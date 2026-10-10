"""Junction-only local path usage for a Tealeaf splice block.

Every variable junction of a path is an equal read opportunity, so a molecule
of path s yields junction reads at a rate proportional to its number of
variable junctions |J_s|. This is the usual junction-count PSI estimator
(inclusion junctions averaged for a cassette) written as a compatibility
likelihood with effective lengths, so multi-path blocks are handled like
Tealeaf's EC likelihood but without exon-body or unspliced reads.
"""

from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize
from scipy.special import softmax

from .sashimi import variable_path_features


def block_path_exons(left_anchor, signature, right_anchor):
    """Full local exon chain of one path: anchors plus its internal exons."""
    exons = [tuple(exon) for exon in signature]
    if left_anchor is not None:
        exons.append(tuple(left_anchor))
    if right_anchor is not None:
        exons.append(tuple(right_anchor))
    return exons


def path_junction_map(paths, strand):
    """Variable junctions and the J x S junction-by-path membership matrix.

    paths is a list of S exon chains with 0-based half-open intervals. Returns
    junction keys (exon end, next exon start) in genomic order and a 0/1
    matrix M with M[j, s] = 1 when variable junction j lies on path s. A path
    without a variable junction has a zero column.
    """
    features = variable_path_features(SimpleNamespace(strand=strand), paths)
    junctions = features.loc[features.feature_type.eq("junction")] if len(features) else features
    keys = [(int(row.start), int(row.end)) for row in junctions.itertuples(index=False)]
    membership = np.zeros((len(keys), len(paths)))
    for row_index, row in enumerate(junctions.itertuples(index=False)):
        membership[row_index, np.asarray(row.path_numbers, dtype=int) - 1] = 1.
    return keys, membership


def junction_identifiable(weights, tolerance=1e-8):
    """True when junction counts identify every path proportion."""
    if weights.shape[0] == 0 or (weights.sum(axis=0) <= 0).any():
        return False
    return np.linalg.matrix_rank(weights, tol=tolerance) == weights.shape[1]


def junction_path_usage(counts, membership, concentration=1.):
    """MAP path proportions psi from junction counts.

    counts is a length-J nonnegative vector, membership the J x S matrix M from
    path_junction_map, and concentration the total Dirichlet pseudo-count
    zeta, so each path receives zeta/S as in Tealeaf's fixed-total prior. The
    objective is sum_j c_j log{(M psi)_j / (1'M psi)} + (zeta/S) sum_s log psi_s.
    """
    counts = np.asarray(counts, dtype=float)
    lengths = membership.sum(axis=0)
    prior, total = concentration / membership.shape[1], counts.sum()
    observed = counts > 0

    def objective(logits):
        psi = softmax(logits)
        rate, length = membership @ psi, lengths @ psi
        value = counts[observed] @ np.log(rate[observed]) - total * np.log(length) + prior * np.log(psi).sum()
        score = membership[observed].T @ (counts[observed] / rate[observed]) - total * lengths / length + prior / psi
        gradient = psi * (score - psi @ score)
        return -value, -gradient

    fit = minimize(objective, np.zeros(membership.shape[1]), jac=True, method="L-BFGS-B", options={"maxiter": 1000, "gtol": 1e-9})
    return softmax(fit.x)
