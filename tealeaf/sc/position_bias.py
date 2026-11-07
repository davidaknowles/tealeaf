"""Primer-specific positional capture bias for block-local read opportunities.

Read-end profiles f_m(d3) and g_m(d5) give the relative density of reads of
primer m whose 3'-most transcript base lies d3 bases (spliced) from the
transcript 3' end and d5 bases from its 5' end, relative to the middle of the
transcript. They are estimated label-blind from single-transcript genes
(end_profiles) and keep their last bin's value beyond it; a constant weight
far from both ends cancels within a block. A read start position on a
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


def end_profiles(gene_counts, bin_width=25, n_bins=80, middle=2000, smooth_from=8, window=5, clip=50.):
    """Read-end profiles from per-gene read-end histograms.

    gene_counts is an iterable of (length, counts3, counts5, middle_count)
    with length the spliced transcript length, counts3 / counts5 read-end
    counts in n_bins bins of bin_width bases from the 3' / 5' end, and
    middle_count the reads ending at least middle bases from both ends. Each
    gene's bin counts are divided by its expected count at its middle
    density (ratios clipped at clip) and averaged over genes with weights
    min(expected count, median expected count), so that no single highly
    expressed gene dominates. Bins from smooth_from on are smoothed by a
    centred moving average of window bins. Returns (f, g).
    """
    ratios3, ratios5, expected = [], [], []
    for length, counts3, counts5, middle_count in gene_counts:
        span = length - 2 * middle
        if span <= 0 or middle_count <= 0:
            continue
        value = middle_count / span * bin_width
        expected.append(value)
        ratios3.append(np.clip(np.asarray(counts3, dtype=float) / value, 0., clip))
        ratios5.append(np.clip(np.asarray(counts5, dtype=float) / value, 0., clip))
    if not expected:
        raise ValueError("no gene has reads in a middle region")
    weights = np.minimum(expected, np.median(expected))
    profiles = []
    for ratios in (ratios3, ratios5):
        profile = weights @ np.asarray(ratios) / weights.sum()
        tail = profile[smooth_from:]
        padded = np.r_[np.full(window // 2, tail[0]), tail, np.full(window // 2, tail[-1])]
        profile[smooth_from:] = np.convolve(padded, np.ones(window) / window, mode="valid")
        profiles.append(profile)
    return tuple(profiles)


def profile_weight(profile, distance, bin_width):
    """Relative density at integer distances (array); 0 for negative distances, the last bin's value beyond the profile."""
    distance = np.asarray(distance)
    index = np.clip(np.floor_divide(distance, bin_width), 0, len(profile) - 1)
    return np.where(distance < 0, 0., np.asarray(profile)[index])


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
        _accumulate(opportunities, structures, index, chain, offsets, positions, read_length, weights, variable)
    return {primer: dict(values) for primer, values in opportunities.items()}


def _accumulate(opportunities, structures, index, chain, offsets, positions, read_length, weights, variable):
    """Add each read start's primer weights to the class (mask) of the read it produces on path index."""
    length = offsets[-1]
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
            for primer, weight in weights.items():
                opportunities[primer][mask][index] += weight[row]


def a_run_starts(sequence, minimum=8, window=10):
    """Start positions of A-runs (windows of window bases with at least minimum A, merged) in a sense sequence."""
    is_a = np.frombuffer(sequence.upper().encode(), dtype=np.uint8) == ord("A")
    if len(is_a) < window:
        return np.zeros(0, dtype=int)
    hit = np.convolve(is_a, np.ones(window, dtype=int), "valid") >= minimum
    return np.flatnonzero(hit & ~np.r_[False, hit[:-1]])


def priming_kernel(enrichment, delta, *, smooth=21, low=-300, high=400):
    """Excess read-end density kappa(delta) = max(smoothed enrichment - 1, 0) for low <= delta <= high, else 0.

    enrichment is the observed / uniform read-end density at offset delta = A-run start - read end (transcript
    orientation), from read ends of single-transcript genes (internal priming diagnostic).
    """
    enrichment, delta = np.asarray(enrichment, dtype=float), np.asarray(delta, dtype=int)
    smoothed = np.convolve(enrichment, np.ones(smooth) / smooth, "same")
    keep = (delta >= low) & (delta <= high)
    return delta[keep], np.maximum(smoothed[keep] - 1., 0.)


def priming_weights(sequence, ends, kernel):
    """Weight 1 + sum over A-runs a of kappa(a - e) for read ends e (transcript coordinates of a sense sequence)."""
    delta, kappa = kernel
    ends = np.asarray(ends, dtype=int)
    weights = np.ones(len(ends))
    if not len(ends):
        return weights
    offset = ends.min()
    density = np.zeros(ends.max() - offset + 1)
    for start in a_run_starts(sequence):
        positions = start - delta - offset  # read ends e = a - delta
        inside = (positions >= 0) & (positions < len(density))
        np.add.at(density, positions[inside], kappa[inside])
    return weights + density[ends - offset]


def primed_path_read_opportunities(paths, read_length, strand, sequences, kernel, variable=None, primed="poly(dT)", flat="random hexamer"):
    """Primer-specific opportunities {primer: {mask: length-S vector}} with internal priming of primed (poly(dT)).

    sequences[s] is the sense (transcript-orientation) sequence of path s, its exons concatenated. Each read start
    position of the primed primer is weighted by priming_weights at the read's 3'-most base, so a path whose
    sequence holds A-runs within the kernel's reach of its read ends yields more primed reads; the flat primer
    (random hexamer) has weight 1. A-runs outside the path are not seen.
    """
    structures = path_structures(paths)
    opportunities = {primer: defaultdict(lambda: np.zeros(len(paths))) for primer in (primed, flat)}
    for index, chain in enumerate(structures[0]):
        offsets = np.cumsum([0] + [b - a for a, b in chain])
        length = offsets[-1]
        positions = np.arange(1 - read_length, length)
        ends = positions + read_length - 1 if strand == "+" else length - 1 - positions
        weights = {primed: priming_weights(sequences[index], ends, kernel), flat: np.ones(len(positions))}
        _accumulate(opportunities, structures, index, chain, offsets, positions, read_length, weights, variable)
    return {primer: dict(values) for primer, values in opportunities.items()}
