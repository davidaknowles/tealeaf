"""Diagnostic local evidence, not a path likelihood or a new eligibility rule."""

from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass

from .junction_benchmark import normalize_starsolo_barcode
from .sashimi import _alignment_junctions


def local_read_count_tensor(feature, subjects, levels, primers, lookup):
    """Sparse marker counts to M by P by two types by two classes.

    The lookup keys are (feature, subject, type, primer) and values have
    included/excluded fields. Missing strata remain zero, not absent subjects.
    Subjects and primers retain the caller's declared order.
    """
    import numpy as np

    subjects, levels, primers = tuple(subjects), tuple(levels), tuple(primers)
    if len(set(subjects)) != len(subjects) or len(levels) != 2 or levels[0] == levels[1] or not primers or len(set(primers)) != len(primers):
        raise ValueError("unique subjects/primers and two different cell-type levels required")
    values = np.zeros((len(subjects), len(primers), 2, 2), dtype=np.int64)
    for u, subject in enumerate(subjects):
        for p, primer in enumerate(primers):
            for c, level in enumerate(levels):
                record = lookup.get((feature, subject, level, primer), {})
                pair = np.asarray([record.get(key, 0) for key in ('included', 'excluded')], dtype=float)
                if not np.isfinite(pair).all() or (pair < 0).any() or (pair != np.floor(pair)).any() or (pair > 2**50).any():
                    raise ValueError('nonnegative exact finite marker counts required')
                values[u, p, c] = pair.astype(np.int64)
    return values


@dataclass(frozen=True)
class LocalReadContrast:
    """One event's zero-based intervals and class-unanimous local features.

    Features are present in every transcript in one event class and absent
    from every transcript in the other. Other transcripts are not ruled out.
    Gene assignment and informative read counts therefore remain descriptive.
    """

    event_id: str
    gene_id: str
    chromosome: str
    strand: str
    gene_start: int
    gene_end: int
    local_start: int
    local_end: int
    gene_exons: tuple
    included_exons: tuple
    excluded_exons: tuple
    included_junctions: tuple
    excluded_junctions: tuple


def _overlaps(blocks, intervals):
    return any(a < d and b > c for a, b in blocks for c, d in intervals)


def event_local_read_contrast(event_id, gene_id, annotation, included, excluded, local_interval):
    """Use annotated exon chains inside an explicit event coordinate envelope.

    The envelope must come from the event definition, not its entire gene or
    fitted effect. Arbitrary within-class isoform-specific features are not
    labeled as evidence for the common event class. Coordinate conversion is
    the caller's responsibility. Returned features do not estimate RNA PSI.
    """
    chains = annotation["transcripts"]
    included, excluded = tuple(included), tuple(excluded)
    start, end = local_interval
    if not included or not excluded or set(included) & set(excluded) or not (set(included) | set(excluded)) <= chains.keys() or not 0 <= start < end:
        raise ValueError("disjoint nonempty annotated classes and a half-open envelope required")
    exon_chains = [tuple(chains[key]) for key in included + excluded]
    gene_exons = tuple(sorted({exon for chain in chains.values() for exon in chain}))
    if not gene_exons or any(not 0 <= a < b for chain in chains.values() for a, b in chain):
        raise ValueError("positive zero-based annotated exon intervals required")
    boundaries = sorted({start, end, *(max(start, min(end, x)) for chain in exon_chains for exon in chain for x in exon)})
    features = [], []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        presence = [any(a <= left and right <= b for a, b in chain) for chain in exon_chains]
        inc, exc = presence[:len(included)], presence[len(included):]
        target = 0 if all(inc) and not any(exc) else 1 if all(exc) and not any(inc) else None
        if target is not None and left < right:
            features[target].append((left, right))
    junction_sets = []
    for chain in exon_chains:
        ordered = sorted(chain)
        junction_sets.append({(a[1], b[0]) for a, b in zip(ordered, ordered[1:]) if a[1] < b[0] and start <= a[1] and b[0] <= end})
    inc_sets, exc_sets = junction_sets[:len(included)], junction_sets[len(included):]
    inc_junctions = set.intersection(*inc_sets) - set.union(*exc_sets)
    exc_junctions = set.intersection(*exc_sets) - set.union(*inc_sets)
    return LocalReadContrast(str(event_id), str(gene_id), annotation["chromosome"], annotation["strand"], min(a for a, _ in gene_exons), max(b for _, b in gene_exons), start, end, gene_exons, tuple(features[0]), tuple(features[1]), tuple(sorted(inc_junctions)), tuple(sorted(exc_junctions)))


def alignment_local_signature(blocks, junctions, event):
    """Five bits, included/excluded exon, included/excluded junction, local.

    No bits imply a gene-exonic but distal alignment when used by the BAM
    collector. Both class bits remain conflicting evidence, never reassigned.
    An exon support bit needs at least one aligned base, not intron overlap.
    """
    blocks, junctions = tuple(blocks), set(junctions)
    return int(_overlaps(blocks, event.included_exons)) | (int(_overlaps(blocks, event.excluded_exons)) << 1) | (int(bool(junctions.intersection(event.included_junctions))) << 2) | (int(bool(junctions.intersection(event.excluded_junctions))) << 3) | (int(_overlaps(blocks, ((event.local_start, event.local_end),))) << 4)


def collect_indexed_local_read_support(bam_path, index_path, barcode_groups, events):
    """Fetch selected genes, OR feature support over a barcode/UB/gene key.

    Only primary, nonsupplementary, NH=1 alignments with CB, UB and a unique
    matching GX gene assignment are used. Keys deduplicate positions/CIGARs
    within this BAM, not across libraries or the upstream quantifier's UMI
    algorithm. Conflicting inclusion/exclusion reads remain a separate class.
    Returned counts are by event, subject, cell type, primer and signature.
    """
    return collect_indexed_library_read_support((bam_path,), (index_path,), barcode_groups, events)


def collect_indexed_library_read_support(bam_paths, index_paths, barcode_groups, events):
    """Union exact CB/UB/gene keys across runs of ONE physical library.

    The caller must supply only that library's barcodes and alignment runs.
    Never union keys across physical libraries. Read support is OR-combined
    across runs as well as positions, so a cross-run conflict stays a conflict.
    Counters expose the per-run key sum and the number removed by union.
    Exact STAR UB union does not reproduce joint upstream UMI error correction.
    """
    import pysam

    bam_paths, index_paths = tuple(bam_paths), tuple(index_paths)
    if not bam_paths or len(bam_paths) != len(index_paths) or len(set(map(str, bam_paths))) != len(bam_paths):
        raise ValueError("distinct runs and one index per run required")

    by_gene = {}
    for event in events:
        by_gene.setdefault(event.gene_id, []).append(event)
    if len({event.event_id for event in events}) != len(events):
        raise ValueError("unique event identities required")
    counts, filters = Counter(), Counter()
    with ExitStack() as stack:
        sources = [stack.enter_context(pysam.AlignmentFile(str(bam), "rb", index_filename=str(index))) for bam, index in zip(bam_paths, index_paths)]
        for gene, local_events in sorted(by_gene.items()):
            chromosomes = {event.chromosome for event in local_events}
            if len(chromosomes) != 1:
                raise ValueError("one chromosome per gene required")
            chromosome = chromosomes.pop()
            left = min(event.gene_start for event in local_events)
            right = max(event.gene_end for event in local_events)
            signatures = {}
            for source in sources:
                if chromosome not in source.references:
                    filters["unavailable_chromosome_events"] += len(local_events)
                    continue
                run_keys = set()
                for alignment in source.fetch(chromosome, left, right):
                    filters["fetched_alignments"] += 1
                    if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                        filters["nonprimary"] += 1
                        continue
                    if not all(alignment.has_tag(tag) for tag in ("NH", "CB", "UB", "GX")):
                        filters["missing_tags"] += 1
                        continue
                    if alignment.get_tag("NH") != 1:
                        filters["multimapped"] += 1
                        continue
                    genes = str(alignment.get_tag("GX")).split(";")
                    if len(genes) != 1 or genes[0].split(".", 1)[0] != gene.split(".", 1)[0]:
                        filters["nonunique_or_other_gene"] += 1
                        continue
                    barcode = normalize_starsolo_barcode(alignment.get_tag("CB"))
                    if barcode not in barcode_groups:
                        filters["unrequested_barcode"] += 1
                        continue
                    blocks = tuple(alignment.get_blocks())
                    junctions = tuple(_alignment_junctions(alignment))
                    for event in local_events:
                        if not _overlaps(blocks, event.gene_exons):
                            continue
                        key = (event.event_id, barcode, str(alignment.get_tag("UB")))
                        run_keys.add(key)
                        signatures[key] = signatures.get(key, 0) | alignment_local_signature(blocks, junctions, event)
                filters["sum_per_run_barcode_umi_event_keys"] += len(run_keys)
            for (event_id, barcode, _), signature in signatures.items():
                counts[(event_id, *barcode_groups[barcode], signature)] += 1
            filters["barcode_umi_event_keys"] += len(signatures)
    filters["cross_run_duplicate_event_keys"] = filters["sum_per_run_barcode_umi_event_keys"] - filters["barcode_umi_event_keys"]
    return counts, filters
