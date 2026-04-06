"""Utility functions for the tealeaf single-cell (sc) module.

Includes helpers for reading alevin-fry equivalence-class output, building
sparse EC-to-transcript matrices, and running the EM transcript quantification.
"""

import gzip
import collections

import numpy as np
import scipy.sparse as sp
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





