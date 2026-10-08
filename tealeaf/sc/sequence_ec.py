"""Exact gene-local sequence opportunity controls, not a production aligner."""

from dataclasses import dataclass

import numpy as np

_BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def exact_sequence_words(sequence, length, *, unstranded=True):
    """Yield (zero-based start, collision-free two-bit word) for valid DNA.

    Unstranded words identify reverse complements. Any non-ACGT character
    invalidates its entire window. This enumerates full exact read sequences,
    not alignment seeds, mismatch-tolerant mappings or compound UMI classes.
    """
    if not isinstance(length, (int, np.integer)) or length <= 0:
        raise ValueError("positive integer sequence word length required")
    if length > len(sequence):
        return
    mask, shift = (1 << (2 * int(length))) - 1, 2 * (int(length) - 1)
    forward = reverse = valid = 0
    for end, character in enumerate(sequence.upper()):
        base = _BASES.get(character)
        if base is None:
            forward = reverse = valid = 0
            continue
        forward = ((forward << 2) | base) & mask
        reverse = (reverse >> 2) | ((3 - base) << shift)
        valid += 1
        if valid >= length:
            yield end + 1 - length, min(forward, reverse) if unstranded else forward


@dataclass(frozen=True)
class SequenceECKernels:
    """K exact classes by T RNA-molecule opportunity kernels for two primers.

    DT conditions on one uniform start in the terminal window (or transcript
    if unspecified), RH uses every valid start as a capture opportunity.
    No observed count, expression estimate, label or external outcome enters.
    Class membership only considers the supplied gene's transcript strings.
    """

    masks: tuple[int, ...]
    dt: np.ndarray
    rh: np.ndarray
    dt_starts: np.ndarray
    rh_starts: np.ndarray
    read_length: int
    end_window: int | None


def exact_sequence_ec_kernels(sequences, read_length, *, end_window=None, unstranded=True):
    """Enumerate exact-match classes and count source positions in each class.

    Repeated sequence positions count separately, unlike counting distinct
    classes. This is a working single-read opportunity control, not a claim
    that actual UMI ECs arise from this exact-match, gene-local mapper.
    """
    sequences = tuple(sequences)
    if not sequences or any(not isinstance(value, str) for value in sequences):
        raise ValueError("one or more transcript sequence strings required")
    if not isinstance(read_length, (int, np.integer)) or read_length <= 0:
        raise ValueError("positive integer read length required")
    if end_window is not None and (not isinstance(end_window, (int, np.integer)) or end_window <= 0):
        raise ValueError("positive terminal start window or None required")
    words = {}
    for transcript, sequence in enumerate(sequences):
        bit = 1 << transcript
        for _, word in exact_sequence_words(sequence, read_length, unstranded=unstranded):
            words[word] = words.get(word, 0) | bit
    if not words:
        raise ValueError("no valid exact read starts in the supplied sequences")
    rows = {}
    size = len(sequences)
    for transcript, sequence in enumerate(sequences):
        minimum_dt_start = 0 if end_window is None else max(0, len(sequence) - read_length + 1 - end_window)
        for start, word in exact_sequence_words(sequence, read_length, unstranded=unstranded):
            membership = words[word]
            if membership not in rows:
                rows[membership] = np.zeros((2, size), dtype=np.int64)
            rows[membership][1, transcript] += 1
            if start >= minimum_dt_start:
                rows[membership][0, transcript] += 1
    masks = tuple(sorted(rows))
    opportunities = np.stack([rows[value] for value in masks])
    dt_starts, rh_starts = opportunities.sum(axis=0)
    dt = np.divide(opportunities[:, 0], dt_starts, out=np.zeros((len(masks), size)), where=dt_starts > 0)
    return SequenceECKernels(masks, dt, opportunities[:, 1].astype(float), dt_starts, rh_starts, int(read_length), None if end_window is None else int(end_window))


def project_sequence_ec_kernels(kernels, membership):
    """Return aligned K_observed-by-T kernels and exact-class match indicators.

    Unknown observed classes get zero rows, not invented read opportunities.
    Duplicate observed membership classes are rejected; splitting one class
    into several observed bins requires an additional generative definition.
    """
    membership = np.asarray(membership)
    if membership.ndim != 2 or membership.shape[1] != kernels.dt.shape[1] or not np.isin(membership, (0, 1)).all():
        raise ValueError("aligned binary observed membership required")
    masks = [sum((1 << int(index)) for index in np.flatnonzero(row)) for row in membership]
    if len(set(masks)) != len(masks) or 0 in masks:
        raise ValueError("nonempty unique observed class memberships required")
    lookup = {value: index for index, value in enumerate(kernels.masks)}
    matched = np.array([value in lookup for value in masks])
    dt, rh = (np.zeros(membership.shape) for _ in range(2))
    for row, value in enumerate(masks):
        if value in lookup:
            dt[row], rh[row] = kernels.dt[lookup[value]], kernels.rh[lookup[value]]
    return (dt, rh), matched


def sequence_ec_background_kernels(kernels, membership, fraction):
    """Mix fixed compatibility-uniform background without discarding counts.

    Fraction is a declared per-transcript mixture component, not estimated
    sequencing-error probability. It allows observed classes absent from the
    single-read control. DT molecule exposure is one, RH exposure is its full
    valid-start count. Correct-read mass on unretained classes is censored,
    so columns are NOT renormalized after projecting to the observed rows.
    This is a working contamination model and requires real-data validation.
    """
    if not np.isfinite(fraction) or not 0 <= fraction < 1:
        raise ValueError("background fraction must lie in [0,1)")
    maps, matched = project_sequence_ec_kernels(kernels, membership)
    if fraction == 0:
        return maps, matched
    if (kernels.dt_starts <= 0).any() or (kernels.rh_starts <= 0).any():
        raise ValueError("each transcript needs positive sequence capture exposure")
    membership = np.asarray(membership, dtype=float)
    degrees = membership.sum(axis=0)
    if (degrees <= 0).any():
        raise ValueError("each transcript needs observed compatibility support")
    background = membership / degrees
    return ((1 - fraction) * maps[0] + fraction * background, (1 - fraction) * maps[1] + fraction * background * kernels.rh_starts), matched


def replace_ec_kernel_rows_inplace(designs, ec_indices, transcript_indices, kernels):
    """Replace K-by-T weights without changing any global compatibility.

    Mutates data buffers of caller-owned canonical floating CSR matrices.
    Validate all primers before mutation. Rows must contain exclusively the
    declared transcripts; arbitrary row/column order is allowed. Explicit
    zeros and changed support are rejected. Copy large maps once per builder,
    not once per gene. Count matrices are not arguments and cannot change.
    """
    from scipy.sparse import isspmatrix_csr

    rows, columns = np.asarray(ec_indices), np.asarray(transcript_indices)
    if any(value.ndim != 1 or value.dtype.kind not in "iu" or len(value) == 0 or len(np.unique(value)) != len(value) for value in (rows, columns)):
        raise ValueError("unique nonempty integer EC and transcript indices required")
    if len(designs) == 0 or len(designs) != len(kernels):
        raise ValueError("one kernel per primer required")
    lookup = {int(value): index for index, value in enumerate(columns)}
    patches = []
    for design, kernel in zip(designs, kernels):
        if not isspmatrix_csr(design) or not design.has_canonical_format or design.dtype.kind != "f" or not design.data.flags.writeable:
            raise ValueError("writable canonical floating CSR design required")
        if (rows < 0).any() or (rows >= design.shape[0]).any() or (columns < 0).any() or (columns >= design.shape[1]).any():
            raise ValueError("EC or transcript index outside design")
        kernel = np.asarray(kernel, dtype=float)
        if kernel.shape != (len(rows), len(columns)) or not np.isfinite(kernel).all() or (kernel < 0).any():
            raise ValueError("finite nonnegative aligned K-by-T kernel required")
        for local_row, global_row in enumerate(rows):
            start, end = design.indptr[global_row:global_row + 2]
            global_columns = design.indices[start:end]
            if any(int(value) not in lookup for value in global_columns):
                raise ValueError("EC row contains transcripts outside declared gene")
            local_columns = np.array([lookup[int(value)] for value in global_columns], dtype=int)
            support = np.zeros(len(columns), dtype=bool)
            support[local_columns] = True
            if not np.array_equal(kernel[local_row] > 0, support) or not np.isfinite(design.data[start:end]).all() or (design.data[start:end] <= 0).any():
                raise ValueError("replacement must preserve positive compatibility support")
            values = kernel[local_row, local_columns].astype(design.dtype)
            if not np.isfinite(values).all() or (values <= 0).any():
                raise ValueError("replacement is not representable in design dtype")
            patches.append((design, start, end, values.copy()))
    for design, start, end, values in patches:
        design.data[start:end] = values
