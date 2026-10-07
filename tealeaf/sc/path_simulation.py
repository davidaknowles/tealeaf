"""Count-level local-path nulls with fixed observed EC designs and depth."""

import numpy as np
from scipy.special import softmax

from .ec_glmm import ECGLMMData


def simulate_counts(base, baseline, subjects, rng, subject_scale=.5, *, labels=None, path_index=None, residual_concentration=None):
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
    generated = []
    for counts, mapping in zip(base.counts, base.compatibility):
        mass = observation_weights @ np.asarray(mapping).T
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
