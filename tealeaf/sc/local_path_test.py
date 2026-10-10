"""Block-local path tests from read-compatibility class counts.

Each S-path block is tested by S one-versus-rest binary random-subject EC
models (tealeaf.sc.path_marginal_grid), with the other paths collapsed by
label-blind pooled shares, and the S p-values are combined by Simes. Two-path
blocks reduce to a single binary test.
"""

import numpy as np

from .differential import helmert_basis
from .local_path_reads import pooled_path_shares
from .path_marginal import BinaryECPathLikelihood
from .path_marginal_grid import binary_grid_test, estimate_primer_offset

PRIMERS = ("poly(dT)", "random hexamer")


def local_class_masks(opportunities, n_paths, anchors):
    """Observed read classes; without anchors the all-paths class is dropped."""
    full = (1 << n_paths) - 1
    return sorted(mask for mask in opportunities if anchors or mask != full)


def local_read_likelihood(lookup, block, opportunities, subjects, labels, levels, anchors, target=0, shares=None):
    """Primer-specific block-local class counts, target path versus the rest.

    The compatibility of class k with the target is its read-opportunity count
    n[k][target]; with the collapsed rest it is sum_j w_j n[k][j] over the
    other paths, w their label-blind pooled shares renormalized. For two paths
    this is the plain binary model. Normalizers are effective local lengths.
    """
    n_paths = len(next(iter(opportunities.values())))
    masks = local_class_masks(opportunities, n_paths, anchors)
    rest = np.ones(n_paths) if shares is None else np.array(shares, dtype=float)
    rest[target] = 0.
    rest /= rest.sum()
    components = np.column_stack([np.zeros(len(masks)), [opportunities[mask][target] for mask in masks], [opportunities[mask] @ rest for mask in masks]])
    if (components[:, 1].sum() == 0) or (components[:, 2].sum() == 0):
        raise ValueError("a path has no read opportunity in the chosen classes")
    keys, rows = [], []
    for subject in np.unique(subjects):
        for label in np.unique(labels[subjects == subject]):
            values = [np.array([lookup.get((block, subject, levels[int(label)], primer, mask), 0.) for mask in masks]) for primer in PRIMERS]
            if sum(value.sum() for value in values) > 0:
                keys.append((subject, label))
                rows.append(values)
    if not rows:
        raise ValueError("no block-local molecules")
    counts = tuple(np.asarray([row[primer] for row in rows]) for primer in range(len(PRIMERS)))
    n = len(keys)
    return BinaryECPathLikelihood(counts, (components, components), np.asarray([key[0] for key in keys]), np.asarray([key[1] for key in keys]), np.ones((n, 2)), np.zeros(n))


def simes(p_values):
    """Simes combination of a block's per-path p-values."""
    ordered = np.sort(np.asarray(p_values, dtype=float))
    return float(min(1., np.min(len(ordered) * ordered / np.arange(1, len(ordered) + 1))))


def fit_local_block(lookup, block, opportunities, n_paths, subjects, labels, levels, *, anchors=True, primer_offset=True, nodes=9):
    """Test one block contrast; returns summary fields and per-level path usage.

    lookup maps (block, subject, cell type, primer, mask) -> molecules, labels
    index levels (0 = level a). Failed target fits enter Simes as p = 1. The
    effect is the vector of per-path standardized proportion differences
    (level b minus level a), zero where a target fit failed.
    """
    masks = local_class_masks(opportunities, n_paths, anchors)
    pooled = {mask: sum(lookup.get((block, subject, levels[int(label)], primer, mask), 0.) for subject, label in set(zip(subjects, labels)) for primer in PRIMERS) for mask in masks}
    shares = pooled_path_shares(pooled, {mask: opportunities[mask] for mask in masks})
    targets = [0] if n_paths == 2 else list(range(n_paths))
    p_values, statistics, offsets, means = [], [], [], {}
    for target in targets:
        try:
            likelihood = local_read_likelihood(lookup, block, opportunities, subjects, labels, levels, anchors, target, shares)
            if primer_offset:
                likelihood.primer_offsets = np.array([0., estimate_primer_offset(likelihood)])
                offsets.append(float(likelihood.primer_offsets[1]))
            fit = binary_grid_test(likelihood, nodes=nodes)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            p_values.append(1.)
            continue
        p_values.append(fit["p_value"])
        statistics.append(fit["statistic"])
        if fit["converged"]:
            means[target] = fit["standardized_means"][:, 0]
            subject_count = fit["n_subjects"]
    if n_paths == 2 and 0 in means:
        means[1] = 1 - means[0]
    effects = np.zeros(n_paths)
    for target, values in means.items():
        effects[target] = values[1] - values[0]
    fields = {"p_value": simes(p_values), "statistic": max(statistics, default=0.), "converged": bool(means), "n_subjects": subject_count if means else 0, "n_path_tests": len(targets), "n_converged_path_tests": sum(target in means for target in targets), "path_shares": shares.tolist(), "primer_offsets": offsets, "mean_difference": (helmert_basis(n_paths).T @ effects).tolist() if means else [], "mean_difference_norm": float(np.linalg.norm(effects)) if means else np.nan}
    usage = [(target, level, float(values[level])) for target, values in means.items() for level in (0, 1)]
    return fields, usage
