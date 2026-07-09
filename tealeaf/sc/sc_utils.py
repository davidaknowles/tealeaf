"""Utility functions for the tealeaf single-cell (sc) module.

Includes helpers for reading alevin-fry equivalence-class output, building
sparse EC-to-transcript matrices, and running the EM transcript quantification.
"""

import gzip
import collections

import numpy as np
import scipy.sparse as sp
from scipy.optimize import lsq_linear
from scipy.sparse.linalg import eigsh, svds
from scipy.linalg import svd
from collections import OrderedDict

def smart_open(filename, *args, **kwargs):
    """Open a file, transparently decompressing if the suffix is ``.gz``."""
    return gzip.open(filename, *args, **kwargs) if filename.suffix == ".gz" else open(filename, *args, **kwargs)


def read_alevin_ec(fn):
    """Read ``gene_eqclass.txt.gz`` produced by ``alevin quant --dump-eqclasses``.

    Returns
    -------
    num_genes : int
    num_ec : int
    ecs : collections.OrderedDict
        Mapping from EC index → list of transcript indices.
    """
    ecs = collections.OrderedDict()

    with smart_open(fn) as f: 
        for i,l in enumerate(f):
            if i==0: # first line gives number of features (normally genes but can be transcripts)
                num_genes = int(l.decode().strip())
                continue
            if i==1: # second line gives number of ECs
                num_ec = int(l.decode().strip())
                continue
            l = l.decode().strip().split()
            l = [int(g) for g in l]
            ec_idx = l[-1] # last element of line is EC index
            gene_idx = l[:-1]
            ecs[ec_idx] = gene_idx
    return num_genes, num_ec, ecs




def to_coo(x, shape=None):
    """Build a COO sparse matrix from a list-of-column-id-lists.

    Parameters
    ----------
    x : list[list[int]]
        Each entry ``x[i]`` gives the column indices for non-zero entries in row *i*.
    shape : tuple, optional
        Explicit ``(rows, cols)`` shape; inferred from the data if omitted.
    """
    nnz = sum([len(g) for g in x])

    indices = np.zeros((2,nnz), dtype=int) # cell then EC idx

    nz_idx = 0
    for row_idx,col_ids in enumerate(x): 
        nhere = len(col_ids)
        indices[0,nz_idx:nz_idx+nhere] = row_idx
        indices[1,nz_idx:nz_idx+nhere] = list(col_ids) # might be a set
        nz_idx += nhere
    
    return sp.coo_matrix((np.ones(nnz), indices), shape = shape)


def sparse_sum(x, dim):
    return np.squeeze(np.asarray(x.sum(dim)))


def get_fasta(fasta_file, first_field = False):
    with smart_open(fasta_file) as f:
        F = f.read().split(">")
    dic = OrderedDict()
    for x in F:
        x = x.split("\n")
        seq =  "".join("".join(x[1:]).split("\r"))
        seq=seq.upper()
        if len(x) <= 1: continue
        if first_field:
            dic[x[0].split()[0].strip()] = seq
        else:
            dic[x[0].strip()] = seq

    return dic

def get_transcript_lengths(fasta_file):
    cdna = get_fasta(fasta_file, first_field = True)
    return OrderedDict([(transcript,len(seq)) for transcript,seq in cdna.items()])



def get_feature_weights(features, transcript_lengths_dic,  fragment_size = 300):
    feature_lengths = np.array([transcript_lengths_dic[g] for g in features])
    eff_lens = np.array([ 
        ((g-fragment_size) if (g>fragment_size) else g) 
        for g in feature_lengths ]) # discontinous :(
    return feature_lengths, 1. / eff_lens


def EM(counts, ec_transcript_mat, w, iterations=30):
    """Expectation-Maximisation transcript quantification.

    Parameters
    ----------
    counts : array-like, shape (n_ec,)
        Observed equivalence-class counts for one pseudobulk sample.
    ec_transcript_mat : scipy.sparse.coo_matrix, shape (n_ec, n_transcripts)
        Binary EC-to-transcript membership matrix.
    w : ndarray, shape (n_transcripts,)
        Transcript weight vector (typically 1 / effective_length).
    iterations : int
        Number of EM iterations.

    Returns
    -------
    alpha : ndarray, shape (n_transcripts,)
        Estimated relative abundance for each transcript.
    """
    n_transcripts = len(w)
    alpha = np.full(n_transcripts, 1.0 / n_transcripts)  # uniform initialisation

    alpha_w = ec_transcript_mat.copy()

    for i in range(iterations):
        # alpha_w[e, t] = alpha_t * w_t  (only for transcripts in EC e)
        alpha_w.data = (alpha * w)[ec_transcript_mat.col]
        ec_sums = sparse_sum(alpha_w, 1)
        z = sp.diags(counts / ec_sums) @ alpha_w
        alpha_new = sparse_sum(z, 0)
        alpha_new /= alpha_new.sum()

        if i == iterations - 1:
            print(i, np.mean(np.abs(alpha - alpha_new)))

        alpha = alpha_new

    return alpha


def weighted_ec_transcript_matrix(ec_transcript_mat, w, ec_effective_lengths=None):
    """Build the GLM compatibility matrix ``A`` from EC membership.

    ``docs/glm.tex`` defines ``A[s, t] = u_s / l_t`` for compatible EC/isoform
    pairs.  Alevin-fry EC dumps do not include EC effective lengths, so callers
    pass ``ec_effective_lengths=None`` to use ``u_s = 1``.
    """
    mat = ec_transcript_mat.tocoo(copy=True)
    if ec_effective_lengths is None:
        mat.data = w[mat.col].astype(float, copy=False)
    else:
        u = np.asarray(ec_effective_lengths, dtype=float)
        if len(u) != mat.shape[0]:
            raise ValueError("ec_effective_lengths must have one value per EC row")
        mat.data = u[mat.row] * w[mat.col]
    return mat.tocsr()


def NNLS(counts, ec_transcript_mat, w, ec_effective_lengths=None,
         max_iter=200, tol=1e-6):
    """Non-negative least-squares approximation to EM.

    Solves ``min_phi ||c - A phi||_2^2`` with ``phi >= 0`` for one pseudobulk,
    where ``c`` is the EC count vector normalized to proportions.  The returned
    vector is normalized to sum to one when possible so downstream count scaling
    can reuse the EM code path.
    """
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 0:
        return np.zeros(ec_transcript_mat.shape[1], dtype=float)

    c = counts / total
    A = weighted_ec_transcript_matrix(ec_transcript_mat, w, ec_effective_lengths)
    result = lsq_linear(
        A,
        c,
        bounds=(0, np.inf),
        max_iter=max_iter,
        tol=tol,
        lsmr_tol="auto",
    )
    phi = np.maximum(result.x, 0)
    phi_sum = phi.sum()
    if phi_sum > 0:
        phi /= phi_sum
    return phi


def _svt_nonnegative(x, threshold, rank=None):
    """Apply singular-value soft-thresholding followed by nonnegative projection."""
    min_dim = min(x.shape)
    if rank is not None and 0 < rank < min_dim - 1:
        u, s, vt = svds(x, k=rank)
        order = np.argsort(s)[::-1]
        u, s, vt = u[:, order], s[order], vt[order]
    else:
        u, s, vt = svd(x, full_matrices=False, check_finite=False)

    keep = s > threshold
    if not np.any(keep):
        return np.zeros_like(x)
    shrunk = s[keep] - threshold
    x_new = (u[:, keep] * shrunk) @ vt[keep, :]
    return np.maximum(x_new, 0)


def NNLS_nucnorm(count_matrix, ec_transcript_mat, w, regularization=0.01,
                 ec_effective_lengths=None, max_iter=50, tol=1e-4,
                 svd_rank=50, max_dense_entries=100_000_000):
    """Many-pseudobulk NNLS with nuclear-norm regularization.

    This is a direct reference implementation of the objective in
    ``docs/glm.tex``:

        0.5 * ||C - A Phi||_F^2 + lambda * ||Phi||_*

    using sparse sufficient statistics ``B=A.T@C`` and ``G=A.T@A``.  It stores
    ``Phi`` densely, so callers should keep the pseudobulk/feature matrix
    filtered or raise ``max_dense_entries`` knowingly.
    """
    counts = count_matrix.tocsr()
    n_samples, _ = counts.shape
    A = weighted_ec_transcript_matrix(ec_transcript_mat, w, ec_effective_lengths)
    n_transcripts = A.shape[1]

    dense_entries = n_samples * n_transcripts
    if dense_entries > max_dense_entries:
        raise MemoryError(
            "NNLS_nucnorm would create a dense transcript-by-sample matrix with "
            f"{dense_entries} entries; increase max_dense_entries if intended."
        )

    totals = np.asarray(counts.sum(axis=1)).ravel()
    nonzero = totals > 0
    inv_totals = np.zeros_like(totals, dtype=float)
    inv_totals[nonzero] = 1.0 / totals[nonzero]
    C = counts.T @ sp.diags(inv_totals)

    B = (A.T @ C).toarray()
    G = A.T @ A

    try:
        lipschitz = float(eigsh(G, k=1, return_eigenvectors=False)[0])
    except Exception:
        lipschitz = float(sp.linalg.norm(G))
    if not np.isfinite(lipschitz) or lipschitz <= 0:
        lipschitz = 1.0

    phi = np.zeros((n_transcripts, n_samples), dtype=float)
    y = phi.copy()
    momentum = 1.0
    threshold = regularization / lipschitz

    for _ in range(max_iter):
        grad = G @ y - B
        proposal = y - grad / lipschitz
        phi_new = _svt_nonnegative(proposal, threshold, svd_rank)

        denom = max(np.linalg.norm(phi), 1.0)
        rel_change = np.linalg.norm(phi_new - phi) / denom
        next_momentum = (1 + np.sqrt(1 + 4 * momentum * momentum)) / 2
        y = phi_new + ((momentum - 1) / next_momentum) * (phi_new - phi)
        phi, momentum = phi_new, next_momentum
        if rel_change < tol:
            break

    col_sums = phi.sum(axis=0)
    positive = col_sums > 0
    phi[:, positive] /= col_sums[positive]
    phi[:, ~nonzero] = 0
    return phi.T




