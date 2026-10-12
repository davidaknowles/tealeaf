"""Primer-specific positional capture bias for block-local read opportunities.

Read-end profiles f_m(d3) and g_m(d5) give the relative density of reads of
primer m whose 3'-most transcript base lies d3 bases (spliced) from the
transcript 3' end and d5 bases from its 5' end, relative to the middle of the
transcript. They are estimated label-blind from single-transcript genes
(end_profiles) and are 1 beyond their last bin. A read start position on a
block path is weighted by f_m(d3) g_m(d5) averaged over the annotated
transcripts of the path (path_end_distances), and the weighted start
positions give primer-specific read opportunities
(weighted_path_read_opportunities) in place of the uniform counts of
tealeaf.sc.local_path_reads.path_read_opportunities.
"""

from collections import defaultdict

import numpy as np

from .local_path_reads import compatibility_mask, path_structures


def transcript_end_distance(exons, strand, position):
    """Spliced distances (d3, d5) of a genomic base from the transcript ends, or None if not exonic.

    exons are genomic half-open (start, end) intervals of one transcript.
    """
    exons = sorted(exons)
    total = sum(b - a for a, b in exons)
    offset = 0
    for a, b in exons:
        if a <= position < b:
            left = offset + position - a
            return (total - 1 - left, left) if strand == "+" else (left, total - 1 - left)
        offset += b - a
    return None


def end_profiles(gene_counts, bin_width=25, n_bins=80, middle=2000):
    """Read-end profiles from per-gene read-end histograms.

    gene_counts is an iterable of (length, counts3, counts5, middle_count)
    with length the spliced transcript length, counts3 / counts5 read-end
    counts in n_bins bins of bin_width bases from the 3' / 5' end, and
    middle_count the reads ending at least middle bases from both ends.
    Returns (f, g), each n_bins relative densities: summed bin counts over the
    summed expected count at the middle density (ratio of sums over genes).
    """
    observed3, observed5, expected = np.zeros(n_bins), np.zeros(n_bins), 0.
    for length, counts3, counts5, middle_count in gene_counts:
        span = length - 2 * middle
        if span <= 0:
            continue
        observed3 += counts3
        observed5 += counts5
        expected += middle_count / span * bin_width
    if expected <= 0:
        raise ValueError("no gene has a middle region")
    return observed3 / expected, observed5 / expected


def profile_weight(profile, distance, bin_width):
    """Relative density at integer distances (array); 0 for negative distances, 1 beyond the profile."""
    distance = np.asarray(distance)
    index = np.floor_divide(distance, bin_width)
    inside = (distance >= 0) & (index < len(profile))
    return np.where(distance < 0, 0., np.where(inside, np.asarray(profile)[np.clip(index, 0, len(profile) - 1)], 1.))


def path_end_distances(path, transcripts, strand):
    """(D3, D5) per transcript: spliced exonic length beyond the path's 3' and 5' genomic ends.

    path is the block path's genomic exon chain, transcripts a list of
    genomic exon lists of the annotated transcripts of that path.
    """
    low, high = min(a for a, _ in path), max(b for _, b in path)
    distances = []
    for exons in transcripts:
        left = sum(max(0, min(b, low) - a) for a, b in exons)
        right = sum(max(0, b - max(a, high)) for a, b in exons)
        distances.append((right, left) if strand == "+" else (left, right))
    return distances


def weighted_path_read_opportunities(paths, read_length, strand, distances, profiles, bin_width, variable=None, floor=.02):
    """Primer-specific opportunities {primer: {mask: length-S vector}}.

    paths and variable are as in path_read_opportunities. distances[s] is the
    list of (D3, D5) of path s's annotated transcripts (an empty list, e.g.
    for the precursor window, gives uniform weight 1). profiles maps primer
    -> (f, g). Each start position's weight is the transcript-averaged
    f(d3) g(d5) of the read's 3'-most base, read_length - 1 bases downstream
    of the start in transcript orientation, floored at floor.
    """
    structures = path_structures(paths)
    opportunities = {primer: defaultdict(lambda: np.zeros(len(paths))) for primer in profiles}
    for index, chain in enumerate(structures[0]):
        offsets = np.cumsum([0] + [b - a for a, b in chain])
        length = offsets[-1]
        positions = np.arange(1 - read_length, length)
        ends = positions + read_length - 1 if strand == "+" else positions
        weights = {}
        for primer, (f, g) in profiles.items():
            if not distances[index]:
                weights[primer] = np.ones(len(positions))
                continue
            total = np.zeros(len(positions))
            for d3, d5 in distances[index]:
                local3 = d3 + (length - 1 - ends) if strand == "+" else d3 + ends
                local5 = d5 + ends if strand == "+" else d5 + (length - 1 - ends)
                total += profile_weight(f, local3, bin_width) * profile_weight(g, local5, bin_width)
            weights[primer] = np.maximum(total / len(distances[index]), floor)
        for row, position in enumerate(positions):
            first, last = max(position, 0), min(position + read_length, length)
            blocks = []
            for exon, (a, b) in enumerate(chain):
                low, high = max(first, offsets[exon]), min(last, offsets[exon + 1])
                if low < high:
                    blocks.append((a + low - offsets[exon], a + high - offsets[exon]))
            junctions = [(blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1) if blocks[i][1] < blocks[i + 1][0]]
            if variable is not None and not variable.intersection(junctions):
                continue
            mask = compatibility_mask(blocks, junctions, structures)
            if mask:
                for primer in profiles:
                    opportunities[primer][mask][index] += weights[primer][row]
    return {primer: dict(values) for primer, values in opportunities.items()}
