"""Transcript-class event operations, independent of event catalog or dataset."""

import numpy as np

from .ec_glmm import ECGLMMData


def collapse_event_nuisance(data, path_index, baseline):
    """Fix within-class shares for two event paths and an optional outside class."""
    path_index = np.asarray(path_index, dtype=int)
    baseline = np.asarray(baseline, dtype=float)
    if path_index.shape != baseline.shape or baseline.shape != (data.n_isoforms,) or not np.isfinite(baseline).all() or (baseline < 0).any() or baseline.sum() <= 0 or (path_index < -1).any() or (path_index > 1).any():
        raise ValueError("aligned nonnegative transcript weights and binary event paths required")
    groups = [np.flatnonzero(path_index == value) for value in (0, 1, -1)]
    if not len(groups[0]) or not len(groups[1]):
        raise ValueError("both event classes required")
    groups = [group for group in groups if len(group)]
    collapsed_maps = []
    for mapping in data.compatibility:
        columns = []
        for group in groups:
            weights = baseline[group]
            if weights.sum() <= 0:
                weights = np.ones(len(group), dtype=float)
            columns.append(mapping[:, group] @ (weights / weights.sum()))
        collapsed_maps.append(np.column_stack(columns))
    collapsed = ECGLMMData(data.counts, tuple(collapsed_maps), data.design, data.clusters)
    return collapsed, np.asarray([0, 1] + ([-1] if len(groups) == 3 else [])), np.asarray([baseline[group].sum() for group in groups])


def binary_event_information(compatibility, path_index, *, inclusion=.5, event_mass=.7, total_per_type_primer=10000.):
    """Compare interior likelihood information with free versus fixed shares.

    Maps have K_p rows and T transcript columns; path_index has length T.
    Two types share the same positive transcript mixture. Information is for
    their scalar inclusion/exclusion ILR contrast after profiling nuisance.
    This diagnoses local identifiability, not calibration or biological power.
    """
    paths = np.asarray(path_index, dtype=int)
    if not 0 < inclusion < 1 or not 0 < event_mass < 1 or not np.isfinite(total_per_type_primer) or total_per_type_primer <= 0 or not np.any(paths == 0) or not np.any(paths == 1) or np.any((paths < -1) | (paths > 1)):
        raise ValueError("interior binary event proportions, positive totals and both paths required")
    mixture = np.zeros(len(paths))
    mass = event_mass if np.any(paths < 0) else 1.
    for group, value in ((0, mass * inclusion), (1, mass * (1 - inclusion)), (-1, 1 - mass)):
        mask = paths == group
        if mask.any():
            mixture[mask] = value / mask.sum()
    return mixture_event_information(compatibility, paths, mixture, np.full((len(compatibility), 2), total_per_type_primer))


def mixture_event_information(compatibility, path_index, mixture, primer_totals):
    """Expected binary-event contrast information at an explicit T-mixture.

    primer_totals is P by 2, allowing unequal primer/type depth. Transcript
    fractions must be strictly positive so the interior efficient score exists.
    """
    from .path_score_mixed import efficient_shared_path_score

    mixture = np.asarray(mixture, dtype=float)
    paths = np.asarray(path_index, dtype=int)
    totals = np.asarray(primer_totals, dtype=float)
    if mixture.shape != paths.shape or not np.isfinite(mixture).all() or (mixture <= 0).any() or totals.shape != (len(compatibility), 2) or not np.isfinite(totals).all() or (totals < 0).any() or (totals.sum(axis=0) <= 0).any():
        raise ValueError("positive transcript mixture and nonnegative P by 2 primer totals required")
    mixture = mixture / mixture.sum()
    counts = []
    for mapping, depth in zip(compatibility, totals):
        expected = np.asarray(mapping @ mixture).ravel()
        if expected.sum() <= 0:
            if depth.sum() > 0:
                raise ValueError("positive interior primer mass required")
            counts.append(np.zeros((2, len(expected))))
            continue
        counts.append(depth[:, None] * (expected / expected.sum())[None, :])
    data = ECGLMMData(tuple(counts), compatibility, np.ones((2, 1)), np.array([0, 0]))
    _, full, _ = efficient_shared_path_score(data.counts, data.compatibility, np.tile(mixture, (2, 1)), paths, [0, 1], 2)
    collapsed, collapsed_paths, baseline = collapse_event_nuisance(data, paths, mixture)
    _, fixed, _ = efficient_shared_path_score(collapsed.counts, collapsed.compatibility, np.tile(baseline, (2, 1)), collapsed_paths, [0, 1], 2)
    return {"free_transcript_information": float(full[0, 0]), "fixed_mixture_information": float(fixed[0, 0])}
