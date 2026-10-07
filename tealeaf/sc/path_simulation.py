"""Count-level local-path nulls with fixed observed EC designs and depth."""

import numpy as np
from scipy.special import expit, softmax

from .ec_glmm import ECGLMMData


def simulate_independent_binary_blocks(rng, *, n_subjects=12, gene_depth=100, a_effect=0., b_effect=0., subject_sd=.3):
    """Challenge local A inference with an independent B splice change.

    Four transcripts are A0B0,A0B1,A1B0,A1B1. Actual counts are N by 7,
    N=2*n_subjects, for each of two primer maps with known opportunities.
    Read-origin classes are distinct, even if transcript-compatible sets
    coincide. A-origin junction/exon classes are rows 0:4; other gene
    classes are rows 4:7. Subject A offsets are shared across both types.
    B has a type effect but no subject offsets or residual variation.
    Details include the N by 4 truth, label-blind pooled oracle baseline,
    and aligned subject/type labels. No fitted quantity defines the truth.
    """
    if n_subjects < 1 or gene_depth < 1 or not np.isfinite([a_effect, b_effect, subject_sd]).all() or subject_sd < 0:
        raise ValueError("positive sizes and finite effects/nonnegative subject SD required")
    labels = np.tile([0, 1], n_subjects)
    subjects = np.repeat(np.arange(n_subjects), 2)
    mappings = []
    for exon, first_b_excluded, second_b_excluded in ((1., 4., 3.), (3., 6., 5.)):
        mappings.append(np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [exon, exon, 0, 0], [1, 0, first_b_excluded, 0], [0, 2, 0, second_b_excluded], [1, 1, 1, 1]], dtype=float))
    offset = rng.normal(0, subject_sd, n_subjects)[subjects]
    a = expit(offset + a_effect * (2 * labels - 1))
    b = expit(b_effect * (2 * labels - 1))
    weights = np.column_stack([a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)])
    counts = []
    for primer, mapping in enumerate(mappings):
        masses = weights @ mapping.T
        counts.append(np.array([rng.multinomial(gene_depth * (primer + 1), values / values.sum()) for values in masses]))
    data = ECGLMMData(tuple(counts), tuple(mappings), np.ones((len(labels), 1)), subjects)
    details = {"labels": labels, "subjects": subjects, "path_index": np.array([0, 0, 1, 1]), "observation_weights": weights, "baseline": weights.mean(axis=0), "a_usage": a, "b_usage": b, "true_delta": float(a[labels == 1].mean() - a[labels == 0].mean())}
    return data, details


def simulate_counts(base, baseline, subjects, rng, subject_scale=.5, *, labels=None, path_index=None, residual_concentration=None, return_details=False, within_path_type_scale=None):
    """Generate primer-conditioned EC counts under a conditional-mean null.

    Baseline is length T; subjects/labels are length N; path_index is length T.
    With no residual concentration, each subject has one shared transcript
    mixture across all labels, preserving the original measurement-only null.
    Otherwise, a latent S-path composition is independently drawn for each
    subject/label aggregate from Dirichlet(kappa * subject_mean). Its expected
    composition is identical across labels within each subject. All technical
    observations and both primers of that aggregate share the latent draw.
    Within-path transcript ratios and outside-block mass stay subject-specific
    and constant across labels. Depths, compatibility and missingness are fixed.
    Count totals are rounded to the nearest integer.
    Optional details retain the generating subject and observation mixtures
    for diagnostic oracle fits and conditional repeated-count resampling.
    Optional within_path_type_scale changes transcript ratios within each
    path by a label-specific log-normal tilt shared across subjects, while
    preserving every subject/label path mass and outside-block abundance.
    This tests nuisance changes under a genuine target-path null. Its
    observation truth, not the subject baseline, includes the type tilt.
    """
    baseline, subjects = np.asarray(baseline, float), np.asarray(subjects)
    if baseline.shape != (base.n_isoforms,) or np.any(~np.isfinite(baseline)) or np.any(baseline < 0) or baseline.sum() <= 0:
        raise ValueError("a finite nonnegative transcript baseline is required")
    if subjects.shape != (len(base.counts[0]),) or not np.isfinite(subject_scale) or subject_scale < 0:
        raise ValueError("aligned subject labels and a nonnegative subject scale are required")
    levels, encoded = np.unique(subjects, return_inverse=True)
    offsets = rng.normal(scale=subject_scale, size=(len(levels), len(baseline)))
    weights = softmax(np.log(np.maximum(baseline, 1e-12))[None, :] + offsets, axis=1)
    observation_weights = weights[encoded].copy()
    if residual_concentration is not None:
        if not np.isfinite(residual_concentration) or residual_concentration <= 0:
            raise ValueError("residual concentration must be positive and finite")
        if labels is None or path_index is None:
            raise ValueError("residual path variation requires aligned labels and path indices")
        labels, path_index = np.asarray(labels), np.asarray(path_index, int)
        if labels.shape != subjects.shape or path_index.shape != baseline.shape or not np.any(path_index >= 0):
            raise ValueError("residual path variation requires aligned labels and path indices")
        size = int(path_index.max()) + 1
        if size < 2 or not np.array_equal(np.unique(path_index[path_index >= 0]), np.arange(size)):
            raise ValueError("at least two consecutive path categories required")
        for subject, weight in zip(levels, weights):
            mass = np.asarray([weight[path_index == path].sum() for path in range(size)])
            proportions = mass / mass.sum()
            for label in np.unique(labels[subjects == subject]):
                selected = (subjects == subject) & (labels == label)
                latent = rng.dirichlet(residual_concentration * proportions)
                for path in range(size):
                    positions = path_index == path
                    observation_weights[np.ix_(selected, positions)] = weight[positions] * (latent[path] / proportions[path])
    if within_path_type_scale is not None:
        if not np.isfinite(within_path_type_scale) or within_path_type_scale < 0 or labels is None or path_index is None:
            raise ValueError("within-path type changes require aligned labels/paths and a nonnegative scale")
        labels, path_index = np.asarray(labels), np.asarray(path_index, int)
        if labels.shape != subjects.shape or path_index.shape != baseline.shape or not np.any(path_index >= 0):
            raise ValueError("within-path type changes require aligned labels and paths")
        type_levels, encoded_types = np.unique(labels, return_inverse=True)
        tilts = rng.normal(scale=within_path_type_scale, size=(len(type_levels), len(baseline)))
        for path in np.unique(path_index[path_index >= 0]):
            positions = path_index == path
            mass = observation_weights[:, positions].sum(axis=1, keepdims=True)
            shares = softmax(np.log(np.maximum(observation_weights[:, positions], 1e-300)) + tilts[encoded_types][:, positions], axis=1)
            observation_weights[:, positions] = mass * shares
    generated = resample_counts(base, observation_weights, rng)
    if return_details:
        return generated, {"subject_levels": levels, "subject_weights": weights, "observation_weights": observation_weights}
    return generated


def resample_counts(base, weights, rng):
    """Resample actual EC counts conditional on fixed N by T transcript weights.

    Primer counts use their own compatibility and observed rounded totals.
    This does not redraw biological compositions or alter sample missingness.
    """
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (len(base.counts[0]), base.n_isoforms) or not np.isfinite(weights).all() or np.any(weights < 0) or np.any(weights.sum(axis=1) <= 0):
        raise ValueError("finite nonnegative observation-by-transcript weights required")
    generated = []
    for counts, mapping in zip(base.counts, base.compatibility):
        mass = weights @ np.asarray(mapping).T
        totals = np.rint(np.asarray(counts).sum(axis=1)).astype(int)
        sums = mass.sum(axis=1)
        if np.any((totals > 0) & (sums <= 0)):
            raise ValueError("positive-count observation has no simulated compatibility mass")
        draws = np.zeros_like(np.asarray(counts), dtype=float)
        for index, total in enumerate(totals):
            if total:
                draws[index] = rng.multinomial(total, mass[index] / sums[index])
        generated.append(draws)
    return ECGLMMData(tuple(generated), base.compatibility, base.design, base.clusters)
