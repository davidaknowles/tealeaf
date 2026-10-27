"""Intron retention: annotated intron catalogue, sense A-runs and read placement per intron.

Intron coordinates are 0-based half-open (first intronic base, end), so the
intron between exons [a, b) and [c, d) is [b, c). An intron is clean when no
exon boundary of any transcript (any gene, either strand) lies strictly inside
it: an exon that overlaps an intron without containing it has a boundary
inside, so clean introns overlap no such exon and never partially overlap each
other. retained_by_ri marks introns inside an exon of a retained_intron-type
transcript of the same gene, retained_by_any inside an exon of any transcript
of the gene.

Internal priming. Oligo-dT primers also anneal to A-rich stretches of the
RNA, so poly(dT) reads pile up upstream (5', in transcript orientation) of
genomic A-runs on the sense strand, intronic ones included. a_run_positions
finds the runs and downstream_run_distance gives, for a position, the distance
to the next run downstream on the sense strand.

Differential retention. retention_test fits, for one intron and two cell
types, the random-subject Beta-binomial model of the block path tests
(path_marginal_grid) to retained r and spliced n molecules per subject and
cell type, with the retained fraction as the path proportion and a fixed
row offset, typically the logit of the gene's retained fraction over its
other introns in the same row, so that a gene-wide change in unspliced RNA
is not read as retention of the intron. The likelihood-ratio statistic is
calibrated against fits with cell-type labels swapped at random within
subjects: above its lower tail that null is a chi-square on one degree of
freedom scaled by a constant (null_inflation), and calibrated_p_value divides
the statistic by it.
"""

import gzip

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .differential import _ATTRIBUTE
from .path_marginal import BinaryECPathLikelihood
from .path_marginal_grid import binary_grid_test


def read_gtf_transcripts(path):
    """{transcript_id: {gene_id, chromosome, strand, type, exons}} with sorted 0-based half-open exons."""
    transcripts = {}
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = dict(_ATTRIBUTE.findall(fields[8]))
            record = transcripts.setdefault(attributes["transcript_id"], {"gene_id": attributes["gene_id"], "chromosome": fields[0], "strand": fields[6], "type": attributes.get("transcript_type", ""), "exons": []})
            record["exons"].append((int(fields[3]) - 1, int(fields[4])))
    for record in transcripts.values():
        record["exons"].sort()
    return transcripts


def intron_catalogue(transcripts, genes=None):
    """Unique introns (gene, chromosome, strand, start, end) of the transcripts of genes, with clean and retained labels."""
    introns = {}
    for transcript, record in transcripts.items():
        if genes is not None and record["gene_id"] not in genes:
            continue
        for (_, b), (c, _) in zip(record["exons"], record["exons"][1:]):
            if c > b:
                introns.setdefault((record["gene_id"], record["chromosome"], record["strand"], b, c), []).append(transcript)
    table = pd.DataFrame(list(introns), columns=["gene_id", "chromosome", "strand", "start", "end"])
    table["n_transcripts"] = [len(members) for members in introns.values()]
    exons = pd.DataFrame([(record["chromosome"], a, b, record["gene_id"], record["type"]) for record in transcripts.values() for a, b in record["exons"]], columns=["chromosome", "exon_start", "exon_end", "gene_id", "type"])
    table["clean"] = False
    for chromosome, local in table.groupby("chromosome"):
        on = exons[exons.chromosome.eq(chromosome)]
        bounds = np.sort(np.r_[on.exon_start.to_numpy(), on.exon_end.to_numpy()])
        inside = np.searchsorted(bounds, local.end.to_numpy(), side="left") - np.searchsorted(bounds, local.start.to_numpy(), side="right")
        table.loc[local.index, "clean"] = inside == 0
    pairs = table.reset_index().merge(exons, on=["chromosome", "gene_id"])
    pairs = pairs[(pairs.exon_start <= pairs.start) & (pairs.exon_end >= pairs.end)]
    table["retained_by_any"] = table.index.isin(pairs["index"])
    table["retained_by_ri"] = table.index.isin(pairs.loc[pairs.type.eq("retained_intron"), "index"])
    return table.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)


def a_run_positions(codes, base, window=10, min_count=8):
    """Start offsets of the windows of length window holding at least min_count of base (2-bit code)."""
    hits = np.r_[0, np.cumsum(codes == base, dtype=np.int32)]
    return np.flatnonzero(hits[window:] - hits[:-window] >= min_count)


def downstream_run_distance(positions, runs, plus):
    """Distance from each position to the next sense A-run downstream (inf without one).

    runs are sorted genomic centres of the runs on the gene's sense sequence
    (A-rich windows for plus, T-rich windows for minus genes); downstream is
    increasing for plus and decreasing for minus genes.
    """
    positions = np.asarray(positions)
    distance = np.full(len(positions), np.inf)
    if not len(runs):
        return distance
    if plus:
        index = np.searchsorted(runs, positions, side="left")
        found = index < len(runs)
        distance[found] = runs[index[found]] - positions[found]
    else:
        index = np.searchsorted(runs, positions, side="right") - 1
        found = index >= 0
        distance[found] = positions[found] - runs[index[found]]
    return distance


def contained_introns(starts, ends, intron_starts, intron_ends, overhang=0):
    """(block index, intron index) pairs for introns inside aligned blocks with overhang bases on each side.

    Blocks [starts, ends) and introns are on one chromosome; introns must be
    sorted and pairwise disjoint (clean introns), so the introns inside a block
    form one contiguous index range.
    """
    starts, ends = np.asarray(starts), np.asarray(ends)
    first = np.searchsorted(intron_starts, starts + overhang, side="left")
    last = np.searchsorted(intron_ends, ends - overhang, side="right")
    n = np.clip(last - first, 0, None)
    block = np.repeat(np.arange(len(starts)), n)
    intron = np.repeat(first, n) + np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    return block, intron


def gene_baseline_offset(retained_other, spliced_other, pseudocount=.5):
    """Logit of the retained fraction over a row's other introns."""
    return np.log((np.asarray(retained_other, float) + pseudocount) / (np.asarray(spliced_other, float) + pseudocount))


def retention_test(retained, spliced, subjects, labels, offset=None, nodes=9):
    """Random-subject Beta-binomial test of retention between two labels.

    Rows are (subject, label) aggregates with retained and spliced molecule
    counts; logit of the mean retained fraction is intercept + effect of the
    second label (sorted) + subject effect + offset. Returns binary_grid_test's
    fields; coefficients[1] is the label effect on the logit scale.
    """
    retained, spliced = np.asarray(retained, float), np.asarray(spliced, float)
    components = (np.array([[0., 1., 0.], [0., 0., 1.]]),)
    n = len(retained)
    likelihood = BinaryECPathLikelihood((np.column_stack([retained, spliced]),), components, np.asarray(subjects), np.asarray(labels), np.ones((n, 2)), np.zeros(n))
    return binary_grid_test(likelihood, nodes=nodes, offset=offset)


def null_inflation(null_statistics, levels=(.05, .01, .001)):
    """Scale c of a chi2_1 fitted to null statistics: mean ratio of their upper quantiles to chi2_1's at the given levels."""
    null_statistics = np.asarray(null_statistics, float)
    return float(np.mean([np.quantile(null_statistics, 1 - a) / chi2.isf(a, 1) for a in levels]))


def calibrated_p_value(statistic, inflation):
    """P-value of a likelihood-ratio statistic against the scaled chi2_1 null."""
    return chi2.sf(np.asarray(statistic, float) / inflation, 1)
