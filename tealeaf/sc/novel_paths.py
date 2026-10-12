"""Unannotated block paths from label-blind splice junction support.

Junction support comes from STAR SJ.out.tab files (read_star_junctions),
summed over all cells of a library, so no cell-type label enters. A
junction inside a block window that is not a junction of any annotated
mature path is novel. novel_block_paths turns supported novel junctions into
exon chains, each a modification of one annotated path:

- one novel junction (d, a) replaces the part of a path between the exon
  holding its donor and the exon holding its acceptor (alternative 5' or 3'
  splice site, new exon skipping or joining, or an intron inside an exon);
  a splice site may extend an exon into the neighbouring intron by at most
  max_extension bases;
- two novel junctions (d, x) and (y, a) with x < y and y - x <= max_exon add
  the new exon [x, y) between the donor and acceptor exons.

Novel junctions outside every annotated block window (for example skipping
an exon that is constitutive in the annotation) get new blocks
(novel_junction_blocks): overlapping junctions are clustered, the gene's
annotated transcripts give the local exon chains between the exon holding the
cluster's first donor and the exon holding its last acceptor, clipped to
shared anchors, and novel_block_paths adds the novel paths.

Junctions are genomic (donor exon end, acceptor exon start), 0-based
half-open as in tealeaf.sc.local_path_reads.
"""

from collections import defaultdict

from .local_path_reads import path_structures


def read_star_junctions(libraries, min_motif=1):
    """{(chromosome, donor end, acceptor start): (unique reads, libraries with reads, strand)}.

    libraries is a list of lists of SJ.out.tab paths, one list per library.
    Junctions with STAR motif below min_motif (0 is non-canonical) are
    skipped. strand is "+", "-" or "." (STAR intron strand 1, 2, 0).
    """
    reads, present, strands = defaultdict(int), defaultdict(set), {}
    for index, files in enumerate(libraries):
        for path in files:
            with open(path) as handle:
                for line in handle:
                    chromosome, start, end, strand, motif, _, unique = line.split("\t")[:7]
                    if int(motif) < min_motif or int(unique) == 0:
                        continue
                    key = (chromosome, int(start) - 1, int(end))
                    reads[key] += int(unique)
                    present[key].add(index)
                    strands[key] = ".+-"[int(strand)]
    return {key: (value, len(present[key]), strands[key]) for key, value in reads.items()}


def _donor_exon(chain, donor, max_extension):
    """(index, extension) of the exon whose end can move to donor: inside the exon, or at most max_extension into the next intron."""
    for index, (start, end) in enumerate(chain):
        following = chain[index + 1][0] if index + 1 < len(chain) else None
        if start < donor <= end:
            return index, 0
        if following is not None and end < donor < following and donor - end <= max_extension:
            return index, donor - end
    return None, 0


def _acceptor_exon(chain, acceptor, max_extension):
    """(index, extension) of the exon whose start can move to acceptor."""
    for index, (start, end) in enumerate(chain):
        previous = chain[index - 1][1] if index > 0 else None
        if start <= acceptor < end:
            return index, 0
        if previous is not None and previous < acceptor < start and start - acceptor <= max_extension:
            return index, start - acceptor
    return None, 0


def _apply(chain, donor, acceptor, middle, max_extension):
    """(chain with donor ... acceptor and optional middle exons spliced in, bases of exon extension), or None."""
    (i, left), (j, right) = _donor_exon(chain, donor, max_extension), _acceptor_exon(chain, acceptor, max_extension)
    if i is None or j is None or i > j or (i == j and (middle or not donor < acceptor)):
        return None
    return tuple(chain[:i]) + ((chain[i][0], donor),) + tuple(middle) + ((acceptor, chain[j][1]),) + tuple(chain[j + 1:]), left + right


def novel_block_paths(paths, support, *, min_reads=20, min_libraries=4, min_fraction=.02, max_extension=300, max_exon=1000, max_novel=8):
    """Novel exon chains for one block, most supported first.

    paths are the block's annotated mature exon chains; support maps
    (donor end, acceptor start) -> (reads, libraries) for junctions of the
    block's chromosome. A novel junction inside the block window is used if it
    has at least min_reads reads in at least min_libraries libraries and at
    least min_fraction of the reads of the block's most used annotated
    junction. Each candidate is applied to one annotated path, preferring a
    path where it replaces a junction sharing its donor or acceptor, then the
    least exon extension, then the path whose least supported junction has
    most reads. A junction that pairs with another into a new exon is not
    used alone with an exon extension. Returns at most max_novel
    (chain, reads) pairs, chains distinct from the annotated paths, ordered by
    the smallest read count of their novel junctions.
    """
    exons, junctions, (low, high) = path_structures(paths)
    annotated = set().union(*junctions)
    top = max((support.get(junction, (0, 0))[0] for junction in annotated), default=0)
    threshold = max(min_reads, min_fraction * top)
    novel = {junction: value[0] for junction, value in support.items() if junction not in annotated and low < junction[0] < junction[1] < high and value[0] >= threshold and value[1] >= min_libraries}
    abundance = [min((support.get(junction, (0, 0))[0] for junction in path_junctions), default=0) for path_junctions in junctions]
    candidates, paired = [], set()
    for (donor, x), first in novel.items():
        for (y, acceptor), second in novel.items():
            if donor < x < y < acceptor and y - x <= max_exon:
                candidates.append(((donor, acceptor), ((x, y),), min(first, second)))
                paired.update({(donor, x), (y, acceptor)})
    candidates += [(junction, (), reads) for junction, reads in novel.items()]
    known, chosen = set(exons), {}
    for (donor, acceptor), middle, reads in candidates:
        best = None
        for index, chain in enumerate(exons):
            result = _apply(chain, donor, acceptor, middle, max_extension)
            if result is None or result[0] in known or (not middle and result[1] and (donor, acceptor) in paired):
                continue
            new, extension = result
            new_junctions = {(a[1], b[0]) for a, b in zip(new, new[1:]) if a[1] < b[0]}
            shares = any(d == donor or a == acceptor for d, a in junctions[index] - new_junctions)
            rank = (shares, -extension, abundance[index], -index)
            if best is None or rank > best[0]:
                best = (rank, new)
        if best is not None and chosen.get(best[1], -1) < reads:
            chosen[best[1]] = reads
    ordered = sorted(chosen.items(), key=lambda item: (-item[1], item[0]))[:max_novel]
    return [([list(exon) for exon in chain], reads) for chain, reads in ordered]


def novel_junction_blocks(transcripts, support, windows, *, min_reads=20, min_libraries=4, min_fraction=.02, max_extension=300, **options):
    """New blocks (annotated paths, novel (chain, reads) list) for one gene.

    transcripts are the gene's annotated exon chains, support maps
    (donor end, acceptor start) -> (reads, libraries, ...) for junctions on
    the gene's chromosome and strand, and windows are the (start, end) spans of
    the gene's annotated blocks. A junction is used if it is in no annotated
    transcript, lies inside the gene, is not inside any window, and has at
    least min_reads reads in min_libraries libraries and min_fraction of the
    reads of the most used annotated junction overlapping it. Junctions with
    overlapping spans form one block. options go to novel_block_paths.
    """
    chains = [sorted(tuple(exon) for exon in chain) for chain in transcripts]
    annotated = {(a[1], b[0]) for chain in chains for a, b in zip(chain, chain[1:]) if a[1] < b[0]}
    low, high = min(chain[0][0] for chain in chains), max(chain[-1][1] for chain in chains)
    used = []
    for (donor, acceptor), value in sorted(support.items()):
        if (donor, acceptor) in annotated or not low < donor < acceptor < high or any(start < donor and acceptor < end for start, end in windows):
            continue
        competing = max((support.get(junction, (0,))[0] for junction in annotated if junction[0] < acceptor and junction[1] > donor), default=0)
        if value[0] >= max(min_reads, min_fraction * competing) and value[1] >= min_libraries:
            used.append((donor, acceptor))
    clusters = []
    for donor, acceptor in used:
        if clusters and donor < clusters[-1][1]:
            clusters[-1][1] = max(clusters[-1][1], acceptor)
        else:
            clusters.append([donor, acceptor])
    blocks = []
    for donor, acceptor in clusters:
        segments = []
        for chain in chains:
            (i, _), (j, _) = _donor_exon(chain, donor, max_extension), _acceptor_exon(chain, acceptor, max_extension)
            if i is not None and j is not None and i <= j:
                segments.append(chain[i:j + 1])
        if not segments:
            continue
        start, end = max(segment[0][0] for segment in segments), min(segment[-1][1] for segment in segments)
        paths = sorted({tuple(((max(a, start), min(b, end)) for a, b in segment)) for segment in segments if segment[0][1] > start and segment[-1][0] < end})
        if not paths or not start < donor or not acceptor < end:
            continue
        novel = novel_block_paths([list(path) for path in paths], support, min_reads=min_reads, min_libraries=min_libraries, min_fraction=min_fraction, max_extension=max_extension, **options)
        if novel:
            blocks.append(([[list(exon) for exon in path] for path in paths], novel))
    return blocks
