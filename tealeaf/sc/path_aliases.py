"""Exact transcript-column equivalence diagnostics, never approximate pruning."""

import numpy as np
from scipy import sparse


def exact_path_column_groups(designs, path_index):
    """Group identical columns across every primer within the same path.

    Arrays are K_p by T, indices length T. Cross-path and outside/path
    aliases must remain distinct since merging them changes the target.
    This helper does not fit a collapsed prior or choose an estimator.
    """
    paths = np.asarray(path_index, dtype=int)
    matrices = tuple(np.asarray(sparse.csr_matrix(value).toarray(), dtype=np.float64) for value in designs)
    if paths.ndim != 1 or not matrices or any(value.ndim != 2 or value.shape[1] != len(paths) or not np.isfinite(value).all() for value in matrices):
        raise ValueError("finite aligned primer maps and path indices required")
    groups = {}
    for transcript, path in enumerate(paths):
        signature = []
        for mapping in matrices:
            column = mapping[:, transcript].copy()
            column[column == 0] = 0.  # Canonicalize signed zero, not near-zero.
            signature.append(column.tobytes())
        key = int(path), tuple(signature)
        groups.setdefault(key, []).append(transcript)
    return tuple(np.asarray(group, dtype=int) for group in groups.values())
