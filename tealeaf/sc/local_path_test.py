"""Block-local path tests from read-compatibility class counts.

Each S-path block is tested by S one-versus-rest binary random-subject EC
models (tealeaf.sc.path_marginal_grid), with the other paths collapsed by
label-blind pooled shares, and the S p-values are combined by Simes. Two-path
blocks reduce to a single binary test.
"""

import numpy as np

from .differential import helmert_basis
from .local_path_reads import group_paths, pooled_path_shares, project_mask, project_opportunities
from .path_marginal import BinaryECPathLikelihood
from .path_marginal_grid import binary_grid_test, estimate_primer_offset
from .path_marginal_kernels import profiled_precursor_log_likelihood
from scipy.special import expit

PRIMERS = ("poly(dT)", "random hexamer")


class ProfiledPrecursorLikelihood(BinaryECPathLikelihood):
    """Binary path likelihood with a profiled per-row, per-primer precursor scale.

    Components are K x 3: precursor, target and rest read opportunities.
    Class k has probability proportional to x a_k + psi b_k + (1 - psi) d_k,
    where x >= 0 (precursor molecules relative to mature, in opportunity units)
    is maximized separately for every row, primer and psi. Intronic and
    exon-body reads therefore inform the precursor fraction, and an excess of
    unspliced RNA in one cell type is not read as path usage.
    """

    @property
    def binomial_terms(self):
        return np.full((len(self.subjects), 3), np.nan)

    def primer_log_likelihood(self, row, t, primer):
        """Compiled evaluation; primer_log_likelihood_reference is the numpy version."""
        observed, components = self.counts[primer][row], self.components[primer]
        if observed.sum() <= 0:
            return np.zeros(len(t))
        if components[:, 0].sum() <= 0:
            return self.primer_log_likelihood_reference(row, t, primer)
        seen = observed > 0
        return profiled_precursor_log_likelihood(np.asarray(t, dtype=float), observed[seen], components[seen, 0], components[seen, 1], components[seen, 2], observed.sum(), components[:, 0].sum(), components[:, 1].sum(), components[:, 2].sum())

    def primer_log_likelihood_reference(self, row, t, primer):
        observed, components = self.counts[primer][row], self.components[primer]
        if observed.sum() <= 0:
            return np.zeros(len(t))
        seen = observed > 0
        counts = observed[seen]
        total = observed.sum()
        psi = expit(np.asarray(t, dtype=float))[:, None]
        mature = psi * components[:, 1] + (1 - psi) * components[:, 2]
        precursor, mature_total, mature = components[seen, 0], mature.sum(axis=1), mature[:, seen]
        precursor_total = components[:, 0].sum()
        if precursor_total <= 0:
            return np.log(np.maximum(mature, 1e-300)) @ counts - total * np.log(mature_total)

        grid = np.linspace(-20., 12., 33)
        values = np.stack([np.log(np.maximum(np.exp(u) * precursor + mature, 1e-300)) @ counts - total * np.log(np.exp(u) * precursor_total + mature_total) for u in grid])
        u = grid[np.argmax(values, axis=0)]
        best = values.max(axis=0)
        for _ in range(6):
            scale = np.exp(u)
            inner = scale[:, None] * precursor + mature
            outer = scale * precursor_total + mature_total
            first = scale * ((counts * precursor / inner).sum(axis=1) - total * precursor_total / outer)
            second = first + scale ** 2 * (-(counts * precursor ** 2 / inner ** 2).sum(axis=1) + total * precursor_total ** 2 / outer ** 2)
            step = np.where(second < 0, -first / np.where(second < 0, second, -1.), 0.)
            u = np.clip(u + np.clip(step, -2., 2.), -25., 15.)
        scale = np.exp(u)
        final = np.log(np.maximum(scale[:, None] * precursor + mature, 1e-300)) @ counts - total * np.log(scale * precursor_total + mature_total)
        return np.maximum(final, best)

    def row_log_likelihood(self, row, proportions):
        t = np.log(proportions) - np.log1p(-np.asarray(proportions))
        return sum(self.primer_log_likelihood(row, t, primer) for primer in range(len(self.counts))) - self.offsets[row]


def prepare_block_test(block, entries, opportunities, n_paths, levels, min_share=.02, precursor=False, groups=None, min_subjects=4, min_molecules=20., primer_opportunities=None):
    """Prune a block's paths for one contrast and project its read classes.

    entries are (subject, cell type, primer, mask, molecules) tuples for the
    two levels, with masks over n_paths mature paths then the precursor, and
    opportunities are the matching full-target vectors. Mature paths whose
    label-blind pooled share is below min_share are dropped; masks and
    opportunities are projected onto the retained paths (plus the precursor
    when precursor is set). Raises ValueError if fewer than two paths remain.
    With groups (lists of mature path indices, for junction-only classes),
    paths are first merged by group_paths and kept indexes the retained
    groups, returned as "groups"; precursor must then be False. The
    eligibility screen (min_subjects subjects observed in both levels,
    min_molecules molecules) is repeated on the projected classes, since
    molecules compatible only with a dropped precursor leave the fit.
    primer_opportunities ({primer: opportunities}, e.g. position-weighted,
    with opportunities their mean) are grouped and projected in the same
    way and returned as "primer_opportunities".
    """
    if groups is not None:
        if precursor:
            raise ValueError("path groups are for junction-only classes without a precursor")
        full_paths = n_paths
        entries, opportunities, groups = group_paths(entries, opportunities, groups, n_paths)
        if primer_opportunities is not None:
            regrouped = {primer: group_paths([], values, groups, full_paths) for primer, values in primer_opportunities.items()}
            if any(len(value[2]) != len(groups) for value in regrouped.values()):
                raise ValueError("a path group has no read opportunity for one primer")
            primer_opportunities = {primer: value[1] for primer, value in regrouped.items()}
        n_paths = len(groups)
        if n_paths < 2:
            raise ValueError("fewer than two expressed paths")
    pooled = {}
    for subject, cell_type, _, mask, value in entries:
        pooled[mask] = pooled.get(mask, 0.) + value
    shares = pooled_path_shares(pooled, opportunities)[:n_paths]
    shares = shares / shares.sum()
    kept = [index for index in range(n_paths) if shares[index] >= min_share]
    if len(kept) < 2:
        raise ValueError("fewer than two expressed paths")
    lookup, totals = {}, {}
    for subject, cell_type, primer, mask, value in entries:
        new = project_mask(mask, kept, n_paths, precursor)
        if new:
            lookup[(block, subject, cell_type, primer, new)] = lookup.get((block, subject, cell_type, primer, new), 0.) + value
            totals[(subject, cell_type)] = totals.get((subject, cell_type), 0.) + value
    present = sorted((subject, levels.index(cell_type)) for subject, cell_type in totals if totals[(subject, cell_type)] > 0)
    paired = sum((subject, levels[0]) in totals and (subject, levels[1]) in totals for subject in {subject for subject, _ in present})
    if paired < min_subjects or sum(totals.values()) < min_molecules:
        raise ValueError("too few paired subjects or molecules in the projected classes")
    if primer_opportunities is not None:
        primer_opportunities = {primer: project_opportunities(values, kept, n_paths, precursor) for primer, values in primer_opportunities.items()}
    return {"primer_opportunities": primer_opportunities, "lookup": lookup, "opportunities": project_opportunities(opportunities, kept, n_paths, precursor), "kept": kept, "shares": shares, "subjects": np.array([subject for subject, _ in present]), "labels": np.array([label for _, label in present]), "median_molecules": float(np.median([totals.get((subject, levels[label]), 0.) for subject, label in present])), "groups": groups}


def local_class_masks(opportunities, n_paths, anchors):
    """Observed read classes; without anchors, classes compatible with every
    mature path (with or without the precursor bit) are dropped."""
    full = (1 << n_paths) - 1
    return sorted(mask for mask in opportunities if anchors or (mask & full) != full)


def local_read_likelihood(lookup, block, opportunities, subjects, labels, levels, anchors, target=0, shares=None, n_paths=None, precursor=False, primer_opportunities=None):
    """Primer-specific block-local class counts, target path versus the rest.

    Opportunity vectors have n_paths mature entries, plus a final precursor
    entry when precursor is set. The compatibility of class k with the target
    is n[k][target]; with the collapsed rest it is sum_j w_j n[k][j] over the
    other mature paths, w_j proportional to share_j / L_j, with share the
    label-blind pooled molecule shares (pooled_path_shares) and L_j = sum_k
    n[k][j] over the chosen classes; w is then a per-position abundance on
    the scale of the target's psi, and the rest's molecules split as the
    pooled shares. With the precursor it is n[k][precursor], scaled by a profiled per-row factor.
    For two mature paths without precursor this is the plain binary model.
    With primer_opportunities ({primer: opportunities}), each primer's
    components use its own opportunities over the classes of opportunities.
    """
    n_paths = len(next(iter(opportunities.values()))) - int(precursor) if n_paths is None else n_paths
    masks = local_class_masks(opportunities, n_paths, anchors)
    if shares is None:
        rest = np.ones(n_paths)
    else:
        lengths = np.sum([np.asarray(opportunities[mask][:n_paths], dtype=float) for mask in masks], axis=0)
        rest = np.array(shares[:n_paths], dtype=float) / lengths
    rest[target] = 0.
    rest /= rest.sum()
    def build(values):
        vectors = [np.asarray(values.get(mask, np.zeros(len(opportunities[mask]))), dtype=float) for mask in masks]
        return np.column_stack([[vector[n_paths] if precursor else 0. for vector in vectors], [vector[target] for vector in vectors], [vector[:n_paths] @ rest for vector in vectors]])

    components = tuple(build(primer_opportunities[primer]) for primer in PRIMERS) if primer_opportunities is not None else (build(opportunities),) * len(PRIMERS)
    if any(value[:, 1].sum() == 0 or value[:, 2].sum() == 0 for value in components):
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
    model = ProfiledPrecursorLikelihood if precursor else BinaryECPathLikelihood
    return model(counts, components, np.asarray([key[0] for key in keys]), np.asarray([key[1] for key in keys]), np.ones((n, 2)), np.zeros(n))


def row_compositions(lookup, block, opportunities, subjects, labels, levels, masks, n_paths):
    """Equal-subject mean mature S-path compositions per level, by row EM.

    Each subject/type row (both primers pooled) gets its own effective-length
    EM composition over all supplied paths (precursor included when present,
    then removed and renormalized). Returns a 2 x n_paths array, level a then
    b, averaged over subjects observed in both levels.
    """
    compositions = {}
    for subject in np.unique(subjects):
        for label in np.unique(labels[subjects == subject]):
            counts = {mask: sum(lookup.get((block, subject, levels[int(label)], primer, mask), 0.) for primer in PRIMERS) for mask in masks}
            if sum(counts.values()) > 0:
                shares = pooled_path_shares(counts, {mask: opportunities[mask] for mask in masks})[:n_paths]
                compositions[(subject, int(label))] = shares / shares.sum()
    paired = [subject for subject in np.unique(subjects) if (subject, 0) in compositions and (subject, 1) in compositions]
    if not paired:
        raise ValueError("no subject observed in both levels")
    return np.asarray([np.mean([compositions[(subject, level)] for subject in paired], axis=0) for level in (0, 1)])


def simes(p_values):
    """Simes combination of a block's per-path p-values."""
    ordered = np.sort(np.asarray(p_values, dtype=float))
    return float(min(1., np.min(len(ordered) * ordered / np.arange(1, len(ordered) + 1))))


def fit_local_block(lookup, block, opportunities, n_paths, subjects, labels, levels, *, anchors=True, primer_offset=True, nodes=9, precursor=False, primer_opportunities=None):
    """Test one block contrast; returns summary fields and per-level path usage.

    lookup maps (block, subject, cell type, primer, mask) -> molecules, labels
    index levels (0 = level a). Failed target fits enter Simes as p = 1. The
    effect for two paths is the model's standardized proportion difference;
    for more paths it is the equal-subject mean difference of row EM
    compositions (row_compositions), a joint S-path estimate.
    """
    masks = local_class_masks(opportunities, n_paths, anchors)
    pooled = {mask: sum(lookup.get((block, subject, levels[int(label)], primer, mask), 0.) for subject, label in set(zip(subjects, labels)) for primer in PRIMERS) for mask in masks}
    shares = pooled_path_shares(pooled, {mask: opportunities[mask] for mask in masks})
    targets = [0] if n_paths == 2 else list(range(n_paths))
    p_values, statistics, offsets, means = [], [], [], {}
    for target in targets:
        try:
            likelihood = local_read_likelihood(lookup, block, opportunities, subjects, labels, levels, anchors, target, shares, n_paths=n_paths, precursor=precursor, primer_opportunities=primer_opportunities)
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
    effects = np.zeros(n_paths)
    if n_paths == 2 and 0 in means:
        means[1] = 1 - means[0]
        effects = np.array([means[0][1] - means[0][0], means[1][1] - means[1][0]])
    elif means:
        joint = row_compositions(lookup, block, opportunities, subjects, labels, levels, masks, n_paths)
        effects = joint[1] - joint[0]
        means = {target: joint[:, target] for target in range(n_paths)}
    fields = {"p_value": simes(p_values), "statistic": max(statistics, default=0.), "converged": bool(means), "n_subjects": subject_count if means else 0, "n_path_tests": len(targets), "n_converged_path_tests": sum(target in means for target in targets), "path_shares": shares.tolist(), "primer_offsets": offsets, "mean_difference": (helmert_basis(n_paths).T @ effects).tolist() if means else [], "mean_difference_norm": float(np.linalg.norm(effects)) if means else np.nan}
    usage = [(target, level, float(values[level])) for target, values in means.items() for level in (0, 1)]
    return fields, usage
