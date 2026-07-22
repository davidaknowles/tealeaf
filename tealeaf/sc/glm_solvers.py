"""Memory-bounded genome-wide GLM solvers for single-cell EC counts.

The scalable solvers keep a transcript-by-rank factor and a cell-by-rank
factor.  Count blocks stay sparse, so neither ``A.T @ C`` nor a full
transcript-by-cell abundance matrix is materialized during fitting.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import time
from typing import Iterator

import numpy as np
import scipy.sparse as sp


SCALABLE_METHODS = {
    "admm", "admm_factorized", "frank_wolfe",
    "frank_wolfe_penalized", "factorized",
}


def _torch():
    try:
        import torch
    except ImportError as exc:
        raise ImportError(
            "Torch GLM methods require the optional dependency. Install tealeaf[glm] "
            "or load a PyTorch environment."
        ) from exc
    return torch


def resolve_device(device: str):
    """Resolve ``auto`` to CUDA when available and otherwise CPU."""
    torch = _torch()
    if device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    target = torch.device(device)
    if target.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("--glm_device cuda was requested but CUDA is unavailable")
    return target


def _csr_to_torch(matrix: sp.spmatrix, device, dtype=None):
    torch = _torch()
    matrix = matrix.tocsr()
    dtype = dtype or torch.float32
    return torch.sparse_csr_tensor(
        torch.as_tensor(matrix.indptr, device=device),
        torch.as_tensor(matrix.indices, device=device),
        torch.as_tensor(matrix.data, dtype=dtype, device=device),
        size=matrix.shape,
        device=device,
    )


def _normalized_counts(counts: sp.csr_matrix):
    """Normalize CSR rows once while sharing the immutable sparse structure."""
    counts = counts.tocsr()
    totals = np.asarray(counts.sum(axis=1)).ravel()
    nonempty = totals > 0
    inverse = np.zeros_like(totals, dtype=np.float32)
    inverse[nonempty] = 1.0 / totals[nonempty]
    values = counts.data.astype(np.float32, copy=True)
    values *= np.repeat(inverse, np.diff(counts.indptr))
    normalized = sp.csr_matrix(
        (values, counts.indices, counts.indptr), shape=counts.shape, copy=False
    )
    return normalized, nonempty


def _coo_components(matrix: sp.spmatrix, nonempty, *, pin=False):
    """Build reusable CPU tensors for a sparse cell block."""
    torch = _torch()
    coo = matrix.tocoo(copy=False)
    indices = torch.as_tensor(
        np.vstack((coo.row, coo.col)), dtype=torch.int64
    )
    values = torch.as_tensor(coo.data, dtype=torch.float32)
    mask = torch.as_tensor(nonempty, dtype=torch.bool)
    if pin:
        indices = indices.pin_memory()
        values = values.pin_memory()
        mask = mask.pin_memory()
    return indices, values, mask, matrix.shape


def _components_to_sparse(components, device, *, non_blocking=False):
    torch = _torch()
    indices, values, mask, shape = components
    indices = indices.to(device, non_blocking=non_blocking)
    values = values.to(device, non_blocking=non_blocking)
    mask = mask.to(device, non_blocking=non_blocking)
    sparse = torch.sparse_coo_tensor(
        indices, values, size=shape, device=device
    ).coalesce()
    return sparse, mask


@dataclass
class GLMResult:
    method: str
    left: object | None
    right: object | None
    dense_phi: object | None
    diagnostics: dict


class SparseGLM:
    """Sparse count access and low-rank objective helpers."""

    def __init__(self, counts, compatibility, device="auto", batch_cells=4096,
                 data_backend="auto"):
        cache_started = time.perf_counter()
        self.torch = _torch()
        self.device = resolve_device(device)
        self.counts = counts.tocsr()
        self.a = _csr_to_torch(compatibility, self.device)
        self.a_t = self.a.transpose(0, 1).to_sparse_coo().coalesce()
        self.n_cells, self.n_ec = self.counts.shape
        self.n_transcripts = compatibility.shape[1]
        self.batch_cells = max(1, int(batch_cells))
        self.normalized_counts, self.nonempty = _normalized_counts(self.counts)
        self.n_nonempty = max(1, int(np.count_nonzero(self.nonempty)))
        self.data_backend = self._resolve_data_backend(data_backend)
        self._cell_batches = self._prepare_cell_batches()
        self.cache_build_seconds = time.perf_counter() - cache_started
        self.gram_lipschitz = self._estimate_gram_lipschitz()
        self.total_count_sq = float(np.square(self.normalized_counts.data).sum())

    def _resolve_data_backend(self, backend):
        allowed = {"auto", "cuda", "pinned", "cpu"}
        if backend not in allowed:
            raise ValueError(f"data_backend must be one of {sorted(allowed)}")
        if self.device.type != "cuda":
            if backend in {"cuda", "pinned"}:
                raise ValueError(f"data_backend={backend} requires a CUDA device")
            return "cpu"
        if backend == "cpu":
            raise ValueError("data_backend=cpu is only valid with a CPU device")
        if backend != "auto":
            return backend
        free_bytes, _ = self.torch.cuda.mem_get_info(self.device)
        # Cached COO uses two int64 coordinates plus one float32 value.
        estimated = int(self.normalized_counts.nnz) * 20 + self.n_cells
        return "cuda" if estimated <= 0.7 * free_bytes else "pinned"

    def _prepare_cell_batches(self):
        batches = []
        pin = self.data_backend == "pinned"
        target = self.device if self.data_backend in {"cuda", "cpu"} else None
        for start in range(0, self.n_cells, self.batch_cells):
            stop = min(start + self.batch_cells, self.n_cells)
            components = _coo_components(
                self.normalized_counts[start:stop], self.nonempty[start:stop],
                pin=pin,
            )
            payload = (
                _components_to_sparse(components, target)
                if target is not None else components
            )
            batches.append((start, stop, payload))
        return batches

    def _estimate_gram_lipschitz(self, iterations=8):
        """Estimate the largest eigenvalue of A.T A with sparse power steps."""
        vector = self.torch.full(
            (self.n_transcripts, 1),
            1.0 / np.sqrt(max(self.n_transcripts, 1)),
            device=self.device,
        )
        for _ in range(int(iterations)):
            image = self.gram_times(vector)
            norm = self.torch.linalg.vector_norm(image)
            if norm.item() == 0:
                return 0.0
            vector = image / norm
        return float((vector * self.gram_times(vector)).sum().item())

    def blocks(self, order=None) -> Iterator[tuple[int, int, object, object]]:
        """Yield cached normalized blocks, optionally in a supplied batch order."""
        indices = range(len(self._cell_batches)) if order is None else order
        if self.data_backend != "pinned":
            for index in indices:
                start, stop, payload = self._cell_batches[int(index)]
                block, nonempty = payload
                yield start, stop, block, nonempty
            return
        # Copies are issued on a separate stream. The generator resumes only
        # after the caller has queued work consuming the previously yielded block.
        stream = self.torch.cuda.Stream(device=self.device)

        def schedule(index):
            start, stop, components = self._cell_batches[int(index)]
            with self.torch.cuda.stream(stream):
                loaded = _components_to_sparse(
                    components, self.device, non_blocking=True
                )
                ready = self.torch.cuda.Event()
                ready.record(stream)
            return (start, stop, *loaded, ready)

        iterator = iter(indices)
        try:
            pending = schedule(next(iterator))
        except StopIteration:
            return
        for index in iterator:
            following = schedule(index)
            start, stop, block, nonempty, ready = pending
            current = self.torch.cuda.current_stream(self.device)
            current.wait_event(ready)
            block._indices().record_stream(current)
            block._values().record_stream(current)
            nonempty.record_stream(current)
            yield start, stop, block, nonempty
            pending = following
        start, stop, block, nonempty, ready = pending
        current = self.torch.cuda.current_stream(self.device)
        current.wait_event(ready)
        block._indices().record_stream(current)
        block._values().record_stream(current)
        nonempty.record_stream(current)
        yield start, stop, block, nonempty

    def a_times(self, value):
        return self.torch.sparse.mm(self.a, value)

    def at_times(self, value):
        return self.torch.sparse.mm(self.a_t, value)

    def gram_times(self, value):
        return self.at_times(self.a_times(value))

    def b_times(self, right):
        """Compute ``A.T @ C @ right`` without forming ``C`` or ``A.T @ C``."""
        ec_cross = self.torch.zeros(
            (self.n_ec, right.shape[1]), device=self.device
        )
        for start, stop, block, _ in self.blocks():
            ec_cross += self.torch.sparse.mm(
                block.transpose(0, 1), right[start:stop]
            )
        return self.at_times(ec_cross)

    def factor_statistics(self, right):
        """Return V'V and C V, sufficient statistics for a fixed cell factor."""
        rank = right.shape[1]
        right_gram = self.torch.zeros((rank, rank), device=self.device)
        ec_cross = self.torch.zeros((self.n_ec, rank), device=self.device)
        for start, stop, block, nonempty in self.blocks():
            right_block = right[start:stop] * nonempty[:, None]
            right_gram += right_block.T @ right_block
            ec_cross += self.torch.sparse.mm(
                block.transpose(0, 1), right_block
            )
        return right_gram, ec_cross

    def loss_from_statistics(self, left, right_gram, ec_cross,
                             regularization=0.0, right_sq=None):
        au = self.a_times(left)
        cross = (au * ec_cross).sum()
        left_gram = left.T @ self.at_times(au)
        loss = 0.5 * (
            self.total_count_sq - 2.0 * cross
            + (left_gram * right_gram).sum()
        )
        if regularization:
            right_sq = right_gram.trace() if right_sq is None else right_sq
            loss += 0.5 * regularization * ((left * left).sum() + right_sq)
        return float(loss.detach().cpu())

    def loss_for_factors(self, left, right, regularization=0.0):
        right_gram, ec_cross = self.factor_statistics(right)
        return self.loss_from_statistics(
            left, right_gram, ec_cross, regularization
        )


def _initial_factors(data: SparseGLM, rank: int, seed: int):
    torch = data.torch
    generator = torch.Generator(device=data.device)
    generator.manual_seed(seed)
    # A coefficient row should initially have O(1) total mass because each
    # response row is normalized to sum to one. Scaling both factors this way
    # avoids an initial prediction mass proportional to the transcript count.
    scale = 1.0 / np.sqrt(max(rank * data.n_transcripts, 1))
    left = torch.rand((data.n_transcripts, rank), generator=generator, device=data.device) * scale
    right = torch.rand((data.n_cells, rank), generator=generator, device=data.device) * scale
    totals = np.asarray(data.counts.sum(axis=1)).ravel()
    right[torch.as_tensor(totals <= 0, device=data.device)] = 0
    return left, right


def _coerce_initial_factors(data, initial_factors, *, rank=None,
                            nonnegative=False, prune_zeros=False):
    """Copy and validate a compact coefficient factorization for continuation."""
    if initial_factors is None:
        return None
    if len(initial_factors) != 2:
        raise ValueError("initial_factors must contain left and right factors")
    torch = data.torch
    left = torch.as_tensor(
        initial_factors[0], dtype=torch.float32, device=data.device
    ).detach().clone()
    right = torch.as_tensor(
        initial_factors[1], dtype=torch.float32, device=data.device
    ).detach().clone()
    if left.ndim != 2 or right.ndim != 2:
        raise ValueError("initial factors must be matrices")
    if left.shape[0] != data.n_transcripts or right.shape[0] != data.n_cells:
        raise ValueError("initial factor dimensions do not match the GLM")
    if left.shape[1] != right.shape[1]:
        raise ValueError("initial factors must have the same rank")
    if rank is not None and left.shape[1] != int(rank):
        raise ValueError("initial factor rank does not match rank")
    if not bool(torch.isfinite(left).all() and torch.isfinite(right).all()):
        raise ValueError("initial factors must be finite")
    if nonnegative:
        left.relu_()
        right.relu_()
    if prune_zeros and left.shape[1]:
        contribution = (
            torch.linalg.vector_norm(left, dim=0)
            * torch.linalg.vector_norm(right, dim=0)
        )
        threshold = torch.finfo(left.dtype).eps * contribution.max().clamp_min(1.0)
        keep = contribution > threshold
        left = left[:, keep]
        right = right[:, keep]
    return left, right


def _objective_stop(diagnostics, objective, previous, stable, iteration,
                    tol, min_iter, patience):
    """Update a patient relative-objective stopping rule."""
    if previous is None:
        relative_change = None
        stable = 0
    else:
        relative_change = abs(previous - objective) / max(abs(previous), 1.0)
        stable = stable + 1 if relative_change <= tol else 0
    diagnostics["objective_relative_change"].append(relative_change)
    stop = iteration + 1 >= min_iter and stable >= patience
    return stop, stable


def _validate_stopping(max_iter, min_iter, patience, tol):
    if int(max_iter) < 1 or int(min_iter) < 1 or int(patience) < 1:
        raise ValueError("max_iter, min_iter, and patience must be positive")
    if float(tol) < 0:
        raise ValueError("tol must be nonnegative")


def _expand_initial_factors(data, rank, seed, initial_factors=None):
    """Initialize a rank path, retaining existing columns and adding small ones."""
    warm_factors = _coerce_initial_factors(
        data, initial_factors, nonnegative=True
    )
    if warm_factors is None:
        return (*_initial_factors(data, rank, seed), 0)
    left, right = warm_factors
    warm_rank = int(left.shape[1])
    if warm_rank > int(rank):
        raise ValueError("warm-start rank exceeds requested rank")
    if warm_rank == int(rank):
        return left, right, warm_rank
    torch = data.torch
    generator = torch.Generator(device=data.device)
    generator.manual_seed(seed + warm_rank)
    added_rank = int(rank) - warm_rank
    scale = 1.0 / np.sqrt(max(int(rank) * data.n_transcripts, 1))
    added_left = torch.rand(
        (data.n_transcripts, added_rank),
        generator=generator,
        device=data.device,
    ) * scale
    added_right = torch.rand(
        (data.n_cells, added_rank),
        generator=generator,
        device=data.device,
    ) * scale
    totals = np.asarray(data.counts.sum(axis=1)).ravel()
    added_right[torch.as_tensor(totals <= 0, device=data.device)] = 0
    return (
        torch.cat((left, added_left), dim=1),
        torch.cat((right, added_right), dim=1),
        warm_rank,
    )


def _balance_factor_columns(left, right):
    """Equalize paired column norms without changing their matrix product."""
    torch = _torch()
    left_norm = torch.linalg.vector_norm(left, dim=0)
    right_norm = torch.linalg.vector_norm(right, dim=0)
    active = (left_norm > 0) & (right_norm > 0)
    scale = torch.ones_like(left_norm)
    scale[active] = torch.sqrt(right_norm[active] / left_norm[active])
    return left * scale, right / scale


def _as_sparse_glm(counts, compatibility, device, batch_cells, data_backend):
    if isinstance(counts, SparseGLM):
        if compatibility is not None:
            raise ValueError("compatibility must be omitted for a prepared SparseGLM")
        return counts
    if compatibility is None:
        raise ValueError("compatibility is required with a raw count matrix")
    return SparseGLM(
        counts, compatibility, device, batch_cells, data_backend=data_backend
    )


def _fista_quadratic_right(initial, gram, linear, step, iterations, active=None):
    """Solve independent nonnegative quadratic rows with FISTA steps."""
    torch = _torch()
    current = initial
    extrapolated = current
    momentum_scale = 1.0
    for _ in range(int(iterations)):
        candidate = torch.relu(
            extrapolated - step * (extrapolated @ gram - linear)
        )
        if active is not None:
            candidate[~active] = 0
        next_scale = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * momentum_scale ** 2))
        extrapolated = candidate + (
            (momentum_scale - 1.0) / next_scale
        ) * (candidate - current)
        current = candidate
        momentum_scale = next_scale
    return current


def _fista_quadratic_left(data, initial, right_gram, linear, step, iterations):
    """Solve the nonnegative transcript-factor quadratic with FISTA steps."""
    torch = data.torch
    current = initial
    extrapolated = current
    momentum_scale = 1.0
    for _ in range(int(iterations)):
        gradient = data.gram_times(extrapolated) @ right_gram - linear
        candidate = torch.relu(extrapolated - step * gradient)
        next_scale = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * momentum_scale ** 2))
        extrapolated = candidate + (
            (momentum_scale - 1.0) / next_scale
        ) * (candidate - current)
        current = candidate
        momentum_scale = next_scale
    return current


def _factorized_exact_epoch(data, left, right, learning_rate, inner_steps=1):
    """One exact alternating epoch with accelerated nonnegative subproblems."""
    torch = data.torch
    au = data.a_times(left)
    gram_left_full = data.at_times(au)
    gram_left = left.T @ gram_left_full
    right_lipschitz = float(torch.linalg.matrix_norm(gram_left, ord=2).item())
    right_step = 1.0 / max(right_lipschitz, 1e-12)
    right_gram = torch.zeros(
        (right.shape[1], right.shape[1]), device=data.device
    )
    ec_cross = torch.zeros((data.n_ec, right.shape[1]), device=data.device)
    for start, stop, block, nonempty in data.blocks():
        right_block = right[start:stop]
        linear = torch.sparse.mm(block, au)
        right_block = _fista_quadratic_right(
            right_block, gram_left, linear, right_step, inner_steps, nonempty
        )
        right[start:stop] = right_block
        right_gram += right_block.T @ right_block
        ec_cross += torch.sparse.mm(block.transpose(0, 1), right_block)
    b_times_right = data.at_times(ec_cross)
    left_lipschitz = (
        data.gram_lipschitz
        * float(torch.linalg.matrix_norm(right_gram, ord=2).item())
    )
    # The sparse power estimate is a lower bound on the true spectral norm;
    # retain a small safety margin while taking a scale-correct Lipschitz step.
    left_step = 1.0 / max(1.1 * left_lipschitz, 1e-12)
    left = _fista_quadratic_left(
        data, left, right_gram, b_times_right, left_step, inner_steps
    )
    objective = data.loss_from_statistics(left, right_gram, ec_cross)
    left, right = _balance_factor_columns(left, right)
    return left, right, objective


def _factorized_minibatch_epoch(data, left, right, learning_rate, generator):
    """Random-reshuffling projected updates of cell and transcript factors."""
    torch = data.torch
    order = torch.randperm(
        len(data._cell_batches), generator=generator, device="cpu"
    ).tolist()
    epoch_right_gram = torch.zeros(
        (right.shape[1], right.shape[1]), device=data.device
    )
    epoch_ec_cross = torch.zeros(
        (data.n_ec, right.shape[1]), device=data.device
    )
    for start, stop, block, nonempty in data.blocks(order):
        au = data.a_times(left)
        gram_left_full = data.at_times(au)
        gram_left = left.T @ gram_left_full
        right_lipschitz = float(
            torch.linalg.matrix_norm(gram_left, ord=2).item()
        )
        right_step = min(
            learning_rate, 1.0 / max(right_lipschitz, 1e-12)
        )
        right_block = right[start:stop]
        gradient = right_block @ gram_left - torch.sparse.mm(block, au)
        right_block = torch.relu(right_block - right_step * gradient)
        right_block[~nonempty] = 0
        right[start:stop] = right_block

        right_gram = right_block.T @ right_block
        ec_cross = torch.sparse.mm(block.transpose(0, 1), right_block)
        b_times_right = data.at_times(ec_cross)
        batch_nonempty = nonempty.sum().clamp_min(1)
        gradient_left = (
            gram_left_full @ right_gram - b_times_right
        ) / batch_nonempty
        left_lipschitz = (
            data.gram_lipschitz
            * float(torch.linalg.matrix_norm(right_gram, ord=2).item())
            / float(batch_nonempty.item())
        )
        left_step = min(
            learning_rate, 1.0 / max(left_lipschitz, 1e-12)
        )
        left = torch.relu(left - left_step * gradient_left)
        epoch_right_gram += right_gram
        epoch_ec_cross += ec_cross
    objective = data.loss_from_statistics(
        left, epoch_right_gram, epoch_ec_cross
    )
    left, right = _balance_factor_columns(left, right)
    return left, right, objective


def fit_factorized(counts, compatibility=None, *, rank=64, max_iter=100,
                   tol=1e-4, learning_rate=0.05, device="auto",
                   batch_cells=4096, seed=0, min_iter=10, patience=5,
                   initial_factors=None, minibatch=False, polish_max_iter=32,
                   data_backend="auto", exact_inner_steps=32):
    """Projected alternating optimization of a genome-wide nonnegative factorization."""
    _validate_stopping(max_iter, min_iter, patience, tol)
    if int(polish_max_iter) < 1:
        raise ValueError("polish_max_iter must be positive")
    if int(exact_inner_steps) < 1:
        raise ValueError("exact_inner_steps must be positive")
    data = _as_sparse_glm(
        counts, compatibility, device, batch_cells, data_backend
    )
    torch = data.torch
    if data.device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(data.device)
    left, right, warm_start_rank = _expand_initial_factors(
        data, rank, seed, initial_factors
    )
    diagnostics = {
        "objective": [], "objective_relative_change": [],
        "device": str(data.device), "rank": rank,
        "gram_lipschitz": data.gram_lipschitz,
        "min_iter": int(min_iter), "patience": int(patience), "tolerance": tol,
        "warm_started": initial_factors is not None,
        "warm_start_rank": warm_start_rank,
        "factor_penalty": None,
        "data_backend": data.data_backend,
        "batch_cells": data.batch_cells,
        "minibatch": bool(minibatch),
        "polish_max_iter": int(polish_max_iter),
        "exact_inner_steps": int(exact_inner_steps),
        "phase": [], "epoch_seconds": [],
        "cache_build_seconds": data.cache_build_seconds,
    }
    previous = None
    stable = 0
    stochastic_stable = 0
    phase = "minibatch" if minibatch else "polish"
    polish_budget = min(
        int(polish_max_iter), max(1, int(max_iter) // 4)
    )
    latest_polish_start = max(0, int(max_iter) - polish_budget)
    generator = torch.Generator(device="cpu")
    generator.manual_seed(seed)

    for iteration in range(int(max_iter)):
        transition_after_epoch = False
        if phase == "minibatch" and iteration >= latest_polish_start:
            phase = "polish"
            previous = None
            stable = 0
        started = time.perf_counter()
        if phase == "minibatch":
            left, right, objective = _factorized_minibatch_epoch(
                data, left, right, learning_rate, generator
            )
        else:
            left, right, objective = _factorized_exact_epoch(
                data, left, right, learning_rate, exact_inner_steps
            )
        diagnostics["phase"].append(phase)
        diagnostics["epoch_seconds"].append(time.perf_counter() - started)
        diagnostics["objective"].append(objective)
        if phase == "minibatch":
            relative_change = (
                None if previous is None else
                abs(previous - objective) / max(abs(previous), 1.0)
            )
            diagnostics["objective_relative_change"].append(relative_change)
            loose_tol = max(10.0 * float(tol), 1e-4)
            stochastic_stable = (
                stochastic_stable + 1
                if relative_change is not None and relative_change <= loose_tol
                else 0
            )
            if stochastic_stable >= 3:
                phase = "polish"
                stable = 0
                transition_after_epoch = True
        else:
            stop, stable = _objective_stop(
                diagnostics, objective, previous, stable, iteration,
                tol, int(min_iter), int(patience),
            )
            if stop:
                diagnostics.update(
                    iterations=iteration + 1, converged=True,
                    convergence_reason="objective_patience",
                )
                break
        previous = None if transition_after_epoch else objective
    else:
        diagnostics.update(
            iterations=int(max_iter), converged=False,
            convergence_reason="max_iterations",
        )
    diagnostics["minibatch_iterations"] = diagnostics["phase"].count("minibatch")
    diagnostics["polish_iterations"] = diagnostics["phase"].count("polish")
    diagnostics["cells_per_second"] = (
        data.n_cells * diagnostics["iterations"]
        / max(sum(diagnostics["epoch_seconds"]), 1e-12)
    )
    diagnostics["peak_cuda_memory_bytes"] = (
        int(torch.cuda.max_memory_allocated(data.device))
        if data.device.type == "cuda" else None
    )
    return GLMResult("factorized", left, right, None, diagnostics)


def fit_factorized_admm(counts, compatibility=None, *, rank=64, regularization=0.01,
                        rho=1.0, max_iter=100, tol=1e-4, learning_rate=0.05,
                        device="auto", batch_cells=4096, seed=0,
                        min_iter=10, patience=5, adaptive_rho=True,
                        rho_update_interval=10, rho_balance=10.0,
                        rho_scale=2.0, initial_factors=None,
                        data_backend="auto"):
    """Factor-sized ADMM split for nonnegative low-rank genome-wide fitting.

    This is deliberately distinct from :func:`fit_admm`: it uses the variational
    factorization and therefore is non-convex, while retaining only ``O(r(T+M))``
    state.
    """
    _validate_stopping(max_iter, min_iter, patience, tol)
    if rho <= 0 or rho_update_interval < 1 or rho_balance <= 1 or rho_scale <= 1:
        raise ValueError("rho adaptation parameters are outside their valid range")
    rho = float(rho)
    data = _as_sparse_glm(
        counts, compatibility, device, batch_cells, data_backend
    )
    torch = data.torch
    if data.device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(data.device)
    warm_factors = _coerce_initial_factors(
        data, initial_factors, rank=rank, nonnegative=True
    )
    left, right = (
        _initial_factors(data, rank, seed)
        if warm_factors is None else warm_factors
    )
    left_copy, right_copy = left.clone(), right.clone()
    left_dual, right_dual = torch.zeros_like(left), torch.zeros_like(right)
    diagnostics = {
        "objective": [], "objective_relative_change": [],
        "primal_relative_residual": [], "dual_relative_residual": [],
        "device": str(data.device), "rank": rank,
        "regularization": float(regularization),
        "rho": rho, "gram_lipschitz": data.gram_lipschitz,
        "rho_history": [], "rho_updates": [],
        "adaptive_rho": bool(adaptive_rho),
        "rho_update_interval": int(rho_update_interval),
        "rho_balance": float(rho_balance), "rho_scale": float(rho_scale),
        "min_iter": int(min_iter), "patience": int(patience), "tolerance": tol,
        "warm_started": warm_factors is not None,
        "warm_start_rank": 0 if warm_factors is None else int(left.shape[1]),
        "data_backend": data.data_backend,
        "batch_cells": data.batch_cells,
        "epoch_seconds": [],
        "cache_build_seconds": data.cache_build_seconds,
    }
    previous = None
    stable = 0

    for iteration in range(int(max_iter)):
        started = time.perf_counter()
        au = data.a_times(left)
        gram_left_full = data.at_times(au)
        gram_left = left.T @ gram_left_full
        right_lipschitz = float(
            torch.linalg.matrix_norm(gram_left, ord=2).item()
        ) + regularization + rho
        right_step = min(learning_rate, 1.0 / max(right_lipschitz, 1e-12))
        ec_cross = torch.zeros((data.n_ec, rank), device=data.device)
        right_gram = torch.zeros((rank, rank), device=data.device)
        primal_sq = torch.zeros((), device=data.device)
        dual_sq = torch.zeros((), device=data.device)
        split_norm_sq = torch.zeros((), device=data.device)
        for start, stop, block, nonempty in data.blocks():
            previous_right_copy = right_copy[start:stop].clone()
            right_block = right[start:stop]
            gradient = right_block @ gram_left - torch.sparse.mm(block, au)
            gradient += regularization * right_block + rho * (right_block - right_copy[start:stop] + right_dual[start:stop])
            right_block = right_block - right_step * gradient
            right[start:stop] = right_block
            right_copy[start:stop] = torch.relu(right_block + right_dual[start:stop])
            right_dual[start:stop] += right_block - right_copy[start:stop]
            right_copy[start:stop][~nonempty] = 0
            primal_sq += (right_block - right_copy[start:stop]).square().sum()
            dual_sq += (
                right_copy[start:stop] - previous_right_copy
            ).square().sum()
            split_norm_sq += right_copy[start:stop].square().sum()
            right_gram += right_copy[start:stop].T @ right_copy[start:stop]
            ec_cross += torch.sparse.mm(
                block.transpose(0, 1), right_copy[start:stop]
            )
        b_times_right = data.at_times(ec_cross)
        # Keep the loss, regularization, and augmented-Lagrangian terms on the
        # same scale when stabilizing this block update by the cell count.
        gradient_left = (
            gram_left_full @ right_gram
            - b_times_right
            + regularization * left
            + rho * (left - left_copy + left_dual)
        ) / data.n_nonempty
        left_lipschitz = (
            data.gram_lipschitz
            * float(torch.linalg.matrix_norm(right_gram, ord=2).item())
            / data.n_nonempty
            + (regularization + rho) / data.n_nonempty
        )
        left_step = min(learning_rate, 1.0 / max(left_lipschitz, 1e-12))
        previous_left_copy = left_copy.clone()
        left = left - left_step * gradient_left
        left_copy = torch.relu(left + left_dual)
        left_dual += left - left_copy
        primal_sq += (left - left_copy).square().sum()
        dual_sq += (left_copy - previous_left_copy).square().sum()
        split_norm_sq += left_copy.square().sum()
        primal_norm = torch.sqrt(primal_sq)
        dual_norm = rho * torch.sqrt(dual_sq)
        split_norm = torch.sqrt(split_norm_sq).clamp_min(1e-12)
        primal_relative = float((primal_norm / split_norm).item())
        dual_relative = float((dual_norm / (rho * split_norm)).item())
        left = left_copy
        right = right_copy

        objective = data.loss_from_statistics(
            left, right_gram, ec_cross, regularization
        )
        diagnostics["epoch_seconds"].append(time.perf_counter() - started)
        diagnostics["objective"].append(objective)
        diagnostics["primal_relative_residual"].append(primal_relative)
        diagnostics["dual_relative_residual"].append(dual_relative)
        diagnostics["rho_history"].append(rho)
        stop, stable = _objective_stop(
            diagnostics, objective, previous, stable, iteration,
            tol, int(min_iter), int(patience),
        )
        if max(primal_relative, dual_relative) > tol:
            stable = 0
            stop = False
        if stop:
            diagnostics.update(
                iterations=iteration + 1, converged=True,
                convergence_reason="objective_patience",
            )
            break
        if adaptive_rho and (iteration + 1) % int(rho_update_interval) == 0:
            primal_residual = float(primal_norm.item())
            dual_residual = float(dual_norm.item())
            previous_rho = rho
            if primal_residual > float(rho_balance) * dual_residual:
                rho *= float(rho_scale)
            elif dual_residual > float(rho_balance) * primal_residual:
                rho /= float(rho_scale)
            if rho != previous_rho:
                dual_rescale = previous_rho / rho
                left_dual *= dual_rescale
                right_dual *= dual_rescale
                diagnostics["rho_updates"].append(
                    {"iteration": iteration + 1, "old": previous_rho, "new": rho}
                )
        previous = objective
    else:
        diagnostics.update(
            iterations=int(max_iter), converged=False,
            convergence_reason="max_iterations",
        )
    diagnostics["final_rho"] = rho
    diagnostics["cells_per_second"] = (
        data.n_cells * diagnostics["iterations"]
        / max(sum(diagnostics["epoch_seconds"]), 1e-12)
    )
    diagnostics["peak_cuda_memory_bytes"] = (
        int(torch.cuda.max_memory_allocated(data.device))
        if data.device.type == "cuda" else None
    )
    return GLMResult("admm_factorized", left, right, None, diagnostics)


def fit_admm(counts, compatibility, *, regularization=0.01, rho=1.0,
             max_iter=100, inner_iter=25, tol=1e-4, max_dense_entries=100_000_000,
             device="auto", batch_cells=4096, adaptive_rho=True,
             rho_update_interval=10, rho_balance=10.0, rho_scale=2.0):
    """Bounded dense convex ADMM reference solver for the nuclear-norm objective."""
    data = SparseGLM(counts, compatibility, device, batch_cells)
    torch = data.torch
    if rho <= 0 or rho_update_interval < 1 or rho_balance <= 1 or rho_scale <= 1:
        raise ValueError("rho adaptation parameters are outside their valid range")
    rho = float(rho)
    dense_entries = data.n_transcripts * data.n_cells
    if dense_entries > max_dense_entries:
        raise MemoryError(
            f"admm needs {dense_entries} dense transcript-by-cell entries; "
            f"the configured cap is {max_dense_entries}. Use admm_factorized instead."
        )
    if data.n_ec * data.n_cells > max_dense_entries:
        raise MemoryError("admm would densify an EC-by-cell matrix above the configured cap")
    count_dense = torch.as_tensor(data.counts.toarray(), dtype=torch.float32, device=data.device)
    totals = count_dense.sum(dim=1, keepdim=True).clamp_min(1.0)
    c = (count_dense / totals).T
    a = data.a.to_dense()
    b = a.T @ c
    gram = a.T @ a
    lipschitz = torch.linalg.matrix_norm(gram, ord=2).item() + rho
    phi = torch.zeros((data.n_transcripts, data.n_cells), device=data.device)
    z = phi.clone()
    dual = phi.clone()
    diagnostics = {
        "objective": [], "device": str(data.device), "rho": rho,
        "regularization": float(regularization),
        "rho_history": [], "rho_updates": [],
        "adaptive_rho": bool(adaptive_rho),
    }

    for iteration in range(int(max_iter)):
        candidate = phi
        for _ in range(int(inner_iter)):
            gradient = gram @ candidate - b + rho * (candidate - z + dual)
            candidate = torch.relu(candidate - gradient / lipschitz)
        phi = candidate
        u, singular, vt = torch.linalg.svd(phi + dual, full_matrices=False)
        shrunk = torch.relu(singular - regularization / rho)
        z_previous = z
        z = (u * shrunk) @ vt
        dual += phi - z
        primal = torch.linalg.vector_norm(phi - z).item()
        dual_residual = rho * torch.linalg.vector_norm(z - z_previous).item()
        residual = c - a @ phi
        objective = 0.5 * float((residual * residual).sum().detach().cpu()) + regularization * float(torch.linalg.svdvals(phi).sum().detach().cpu())
        diagnostics["objective"].append(objective)
        diagnostics["rho_history"].append(rho)
        if max(primal, dual_residual) < tol:
            diagnostics.update(iterations=iteration + 1, converged=True, primal_residual=primal, dual_residual=dual_residual)
            break
        if adaptive_rho and (iteration + 1) % int(rho_update_interval) == 0:
            previous_rho = rho
            if primal > float(rho_balance) * dual_residual:
                rho *= float(rho_scale)
            elif dual_residual > float(rho_balance) * primal:
                rho /= float(rho_scale)
            if rho != previous_rho:
                dual *= previous_rho / rho
                lipschitz = torch.linalg.matrix_norm(gram, ord=2).item() + rho
                diagnostics["rho_updates"].append(
                    {"iteration": iteration + 1, "old": previous_rho, "new": rho}
                )
    else:
        diagnostics.update(iterations=int(max_iter), converged=False)
    diagnostics["final_rho"] = rho
    return GLMResult("admm", None, None, torch.relu(phi), diagnostics)


def fit_frank_wolfe(counts, compatibility, *, rank=64, max_atoms=None,
                    tau=None, max_iter=64,
                    power_iter=3, tol=1e-4, device="auto", batch_cells=4096,
                    seed=0, min_iter=10, patience=5):
    """Streaming nonnegative Frank-Wolfe solver with rank-one atom storage."""
    _validate_stopping(max_iter, min_iter, patience, tol)
    data = SparseGLM(counts, compatibility, device, batch_cells)
    torch = data.torch
    max_atoms = int(rank if max_atoms is None else max_atoms)
    if max_atoms < 1 or int(power_iter) < 1:
        raise ValueError("max_atoms and power_iter must be positive")
    max_iter = min(int(max_iter), max_atoms)
    tau = float(np.sqrt(data.n_cells) if tau is None else tau)
    left = torch.empty((data.n_transcripts, 0), device=data.device)
    right = torch.empty((data.n_cells, 0), device=data.device)
    generator = torch.Generator(device=data.device)
    generator.manual_seed(seed)
    diagnostics = {
        "objective": [], "duality_gap": [], "relative_duality_gap": [],
        "line_search_step": [], "device": str(data.device), "rank": rank,
        "max_atoms": max_atoms, "tau": tau, "power_iterations": int(power_iter),
        "min_iter": int(min_iter), "patience": int(patience), "tolerance": tol,
    }
    initial_gap = None
    stable = 0

    for iteration in range(max_iter):
        if left.shape[1]:
            gram_left = data.gram_times(left)
        else:
            gram_left = None
        vector = torch.rand(data.n_cells, generator=generator, device=data.device)
        vector /= torch.linalg.vector_norm(vector).clamp_min(1e-12)
        for _ in range(int(power_iter)):
            left_vector = torch.zeros(data.n_transcripts, device=data.device)
            for start, stop, block, _ in data.blocks():
                b_block = data.at_times(block.transpose(0, 1)).to_dense()
                if gram_left is not None:
                    b_block -= gram_left @ right[start:stop].T
                left_vector += torch.relu(b_block) @ vector[start:stop]
            left_vector = torch.relu(left_vector)
            left_vector /= torch.linalg.vector_norm(left_vector).clamp_min(1e-12)
            next_vector = torch.zeros_like(vector)
            for start, stop, block, _ in data.blocks():
                b_block = data.at_times(block.transpose(0, 1)).to_dense()
                if gram_left is not None:
                    b_block -= gram_left @ right[start:stop].T
                next_vector[start:stop] = torch.relu(b_block).T @ left_vector
            vector = torch.relu(next_vector)
            vector /= torch.linalg.vector_norm(vector).clamp_min(1e-12)

        # Exact quadratic line search for S - Phi using compact factors.
        atom_left = np.sqrt(tau) * left_vector[:, None]
        atom_right = np.sqrt(tau) * vector[:, None]
        direction_left = torch.cat((atom_left, left), dim=1)
        direction_right = torch.cat((atom_right, -right), dim=1)
        a_direction = data.a_times(direction_left)
        b_dot_direction = torch.zeros((), device=data.device)
        for start, stop, block, _ in data.blocks():
            b_dot_direction += (torch.sparse.mm(block, a_direction) * direction_right[start:stop]).sum()
        quadratic = (direction_left.T @ data.gram_times(direction_left) * (direction_right.T @ direction_right)).sum()
        if left.shape[1]:
            cross = (direction_left.T @ data.gram_times(left) * (direction_right.T @ right)).sum()
        else:
            cross = torch.zeros((), device=data.device)
        gap_tensor = torch.relu(b_dot_direction - cross)
        gap = float(gap_tensor.item())
        initial_gap = gap if initial_gap is None else initial_gap
        relative_gap = gap / max(initial_gap, 1e-12)
        gamma = torch.clamp(gap_tensor / quadratic.clamp_min(1e-12), 0.0, 1.0)
        keep = torch.sqrt(1.0 - gamma)
        add = torch.sqrt(gamma)
        left = torch.cat((left * keep, atom_left * add), dim=1)
        right = torch.cat((right * keep, atom_right * add), dim=1)
        objective = data.loss_for_factors(left, right)
        diagnostics["objective"].append(objective)
        diagnostics["duality_gap"].append(gap)
        diagnostics["relative_duality_gap"].append(relative_gap)
        diagnostics["line_search_step"].append(float(gamma.item()))
        stable = stable + 1 if relative_gap <= tol else 0
        if iteration + 1 >= int(min_iter) and stable >= int(patience):
            diagnostics.update(
                iterations=iteration + 1, converged=True,
                convergence_reason="duality_gap_patience",
            )
            break
    else:
        reason = "atom_capacity" if max_iter == max_atoms else "max_iterations"
        diagnostics.update(
            iterations=max_iter, converged=False, convergence_reason=reason,
        )
    return GLMResult("frank_wolfe", left, right, None, diagnostics)


def _penalized_fw_q_times(data, left, right, vectors, penalty):
    """Apply the negative penalized objective gradient to cell vectors."""
    torch = data.torch
    result = torch.zeros(
        (data.n_transcripts, vectors.shape[1]), device=data.device
    )
    if left.shape[1]:
        result -= data.gram_times(left) @ (right.T @ vectors)
    for start, stop, block, _ in data.blocks():
        vector_block = vectors[start:stop]
        result += data.at_times(
            torch.sparse.mm(block.transpose(0, 1), vector_block)
        )
        if penalty and left.shape[1]:
            phi_block = left @ right[start:stop].T
            negative = torch.minimum(phi_block, torch.zeros_like(phi_block))
            result -= penalty * negative @ vector_block
    return result


def _penalized_fw_qt_times(data, left, right, vectors, penalty):
    """Apply the transpose negative penalized objective gradient."""
    torch = data.torch
    result = torch.zeros((data.n_cells, vectors.shape[1]), device=data.device)
    au = data.a_times(vectors)
    if left.shape[1]:
        result -= right @ (left.T @ data.gram_times(vectors))
    for start, stop, block, _ in data.blocks():
        result[start:stop] += torch.sparse.mm(block, au)
        if penalty and left.shape[1]:
            phi_block = left @ right[start:stop].T
            result[start:stop] -= penalty * (
                torch.minimum(phi_block, torch.zeros_like(phi_block)).T @ vectors
            )
    return result


def _penalized_fw_direction(data, left, right, atom_left, atom_right, penalty):
    """Return candidate gap and a quadratic smoothness bound for S - Phi."""
    torch = data.torch
    direction_left = torch.cat((atom_left, left), dim=1)
    direction_right = torch.cat((atom_right, -right), dim=1)
    a_direction = data.a_times(direction_left)
    b_dot_direction = torch.zeros((), device=data.device)
    penalty_dot_direction = torch.zeros((), device=data.device)
    for start, stop, block, _ in data.blocks():
        b_dot_direction += (
            torch.sparse.mm(block, a_direction) * direction_right[start:stop]
        ).sum()
        if penalty and left.shape[1]:
            phi_block = left @ right[start:stop].T
            direction_block = direction_left @ direction_right[start:stop].T
            penalty_dot_direction += (
                torch.minimum(phi_block, torch.zeros_like(phi_block))
                * direction_block
            ).sum()
    if left.shape[1]:
        cross = (
            direction_left.T @ data.gram_times(left)
            * (direction_right.T @ right)
        ).sum()
    else:
        cross = torch.zeros((), device=data.device)
    direction_gram_left = direction_left.T @ direction_left
    direction_gram_right = direction_right.T @ direction_right
    direction_norm_sq = (direction_gram_left * direction_gram_right).sum()
    prediction_norm_sq = (
        direction_left.T @ data.gram_times(direction_left)
        * direction_gram_right
    ).sum()
    gap = b_dot_direction - cross - penalty * penalty_dot_direction
    curvature = prediction_norm_sq + penalty * direction_norm_sq
    return torch.clamp(gap, min=0.0), curvature.clamp_min(1e-12)


def _negative_factor_stats(data, left, right):
    """Measure entrywise negativity of a compact fitted matrix."""
    torch = data.torch
    negative_sq = torch.zeros((), device=data.device)
    total_sq = torch.zeros((), device=data.device)
    negative_l1 = torch.zeros((), device=data.device)
    total_l1 = torch.zeros((), device=data.device)
    for start in range(0, data.n_cells, data.batch_cells):
        stop = min(start + data.batch_cells, data.n_cells)
        block = left @ right[start:stop].T
        negative = torch.minimum(block, torch.zeros_like(block))
        negative_sq += negative.square().sum()
        total_sq += block.square().sum()
        negative_l1 -= negative.sum()
        total_l1 += block.abs().sum()
    return {
        "negative_squared_norm": float(negative_sq.item()),
        "fitted_squared_norm": float(total_sq.item()),
        "negative_frobenius_fraction": float(
            torch.sqrt(negative_sq / total_sq.clamp_min(1e-24)).item()
        ),
        "negative_l1_fraction": float(
            (negative_l1 / total_l1.clamp_min(1e-24)).item()
        ),
    }


def fit_frank_wolfe_penalized(
    counts,
    compatibility,
    *,
    rank=64,
    max_atoms=None,
    tau=None,
    nonnegative_penalty=1.0,
    max_iter=64,
    power_iter=10,
    tol=1e-4,
    device="auto",
    batch_cells=4096,
    seed=0,
    min_iter=10,
    patience=5,
    initial_factors=None,
):
    """Nuclear-ball Frank-Wolfe with a smooth negative-mass penalty.

    The feasible set is the ordinary nuclear-norm ball, whose signed leading
    singular-vector oracle is well defined. ``nonnegative_penalty`` is a
    dimensionless multiplier of the estimated squared spectral norm of A.
    """
    _validate_stopping(max_iter, min_iter, patience, tol)
    if float(nonnegative_penalty) < 0:
        raise ValueError("nonnegative_penalty must be nonnegative")
    data = SparseGLM(counts, compatibility, device, batch_cells)
    torch = data.torch
    max_atoms = int(rank if max_atoms is None else max_atoms)
    if max_atoms < 1 or int(power_iter) < 1:
        raise ValueError("max_atoms and power_iter must be positive")
    tau = float(np.sqrt(data.n_cells) if tau is None else tau)
    penalty = float(nonnegative_penalty) * data.gram_lipschitz
    warm_factors = _coerce_initial_factors(
        data, initial_factors, prune_zeros=True
    )
    if warm_factors is None:
        left = torch.empty((data.n_transcripts, 0), device=data.device)
        right = torch.empty((data.n_cells, 0), device=data.device)
    else:
        left, right = warm_factors
    initial_atoms = int(left.shape[1])
    if initial_atoms > max_atoms:
        raise ValueError("warm start exceeds max_atoms")
    initial_nuclear_upper = float(
        (
            torch.linalg.vector_norm(left, dim=0)
            * torch.linalg.vector_norm(right, dim=0)
        ).sum().item()
    )
    if initial_nuclear_upper > tau:
        scale = np.sqrt(tau / max(initial_nuclear_upper, 1e-12))
        left *= scale
        right *= scale
    fit_iterations = min(int(max_iter), max_atoms - initial_atoms)
    generator = torch.Generator(device=data.device)
    generator.manual_seed(seed)
    diagnostics = {
        "objective": [], "objective_relative_change": [],
        "candidate_gap": [], "relative_candidate_gap": [],
        "line_search_step": [], "oracle_singular_value": [],
        "oracle_relative_residual": [], "device": str(data.device),
        "rank": rank, "max_atoms": max_atoms, "tau": tau,
        "power_iterations": int(power_iter),
        "nonnegative_penalty_multiplier": float(nonnegative_penalty),
        "nonnegative_penalty": penalty,
        "gram_lipschitz": data.gram_lipschitz,
        "min_iter": int(min_iter), "patience": int(patience), "tolerance": tol,
        "gap_is_certificate": False,
        "warm_started": warm_factors is not None,
        "warm_start_rank": initial_atoms,
        "warm_start_nuclear_upper": initial_nuclear_upper,
    }
    initial_gap = None
    previous = None
    if warm_factors is not None:
        initial_stats = _negative_factor_stats(data, left, right)
        previous = data.loss_for_factors(left, right)
        previous += 0.5 * penalty * initial_stats["negative_squared_norm"]
    stable = 0

    for iteration in range(fit_iterations):
        vector = torch.randn(
            (data.n_cells, 1), generator=generator, device=data.device
        )
        vector /= torch.linalg.vector_norm(vector).clamp_min(1e-12)
        for _ in range(int(power_iter)):
            left_vector = _penalized_fw_q_times(
                data, left, right, vector, penalty
            )
            left_vector /= torch.linalg.vector_norm(left_vector).clamp_min(1e-12)
            vector = _penalized_fw_qt_times(
                data, left, right, left_vector, penalty
            )
            vector /= torch.linalg.vector_norm(vector).clamp_min(1e-12)
        q_vector = _penalized_fw_q_times(data, left, right, vector, penalty)
        singular = torch.linalg.vector_norm(q_vector).clamp_min(1e-12)
        left_vector = q_vector / singular
        qt_vector = _penalized_fw_qt_times(
            data, left, right, left_vector, penalty
        )
        oracle_residual = torch.linalg.vector_norm(
            qt_vector - singular * vector
        ) / singular

        atom_left = np.sqrt(tau) * left_vector
        atom_right = np.sqrt(tau) * vector
        gap_tensor, curvature = _penalized_fw_direction(
            data, left, right, atom_left, atom_right, penalty
        )
        gamma = torch.clamp(gap_tensor / curvature, 0.0, 1.0)
        keep = torch.sqrt(1.0 - gamma)
        add = torch.sqrt(gamma)
        left = torch.cat((left * keep, atom_left * add), dim=1)
        right = torch.cat((right * keep, atom_right * add), dim=1)

        objective = data.loss_for_factors(left, right)
        negative_stats = _negative_factor_stats(data, left, right)
        objective += 0.5 * penalty * negative_stats["negative_squared_norm"]
        gap = float(gap_tensor.item())
        initial_gap = gap if initial_gap is None else initial_gap
        relative_gap = gap / max(initial_gap, 1e-12)
        diagnostics["objective"].append(objective)
        stop, stable = _objective_stop(
            diagnostics, objective, previous, stable, iteration,
            tol, int(min_iter), int(patience),
        )
        diagnostics["candidate_gap"].append(gap)
        diagnostics["relative_candidate_gap"].append(relative_gap)
        diagnostics["line_search_step"].append(float(gamma.item()))
        diagnostics["oracle_singular_value"].append(float(singular.item()))
        diagnostics["oracle_relative_residual"].append(float(oracle_residual.item()))
        if stop:
            diagnostics.update(
                iterations=iteration + 1, converged=True,
                convergence_reason="objective_patience",
            )
            break
        previous = objective
    else:
        reason = (
            "atom_capacity"
            if initial_atoms + fit_iterations == max_atoms
            else "max_iterations"
        )
        diagnostics.update(
            iterations=fit_iterations, converged=False, convergence_reason=reason,
        )
    diagnostics.update(_negative_factor_stats(data, left, right))
    return GLMResult("frank_wolfe_penalized", left, right, None, diagnostics)


def fit_glm(counts, compatibility, method, **kwargs):
    """Dispatch a GLM method using a cells-by-EC sparse count matrix."""
    if method == "factorized":
        allowed = {
            "rank", "max_iter", "tol", "learning_rate", "device",
            "batch_cells", "seed", "min_iter", "patience",
            "initial_factors", "minibatch", "polish_max_iter",
            "data_backend", "exact_inner_steps",
        }
        return fit_factorized(counts, compatibility, **{k: v for k, v in kwargs.items() if k in allowed})
    if method == "admm_factorized":
        allowed = {
            "rank", "regularization", "rho", "max_iter", "tol",
            "learning_rate", "device", "batch_cells", "seed", "min_iter",
            "patience", "adaptive_rho", "rho_update_interval",
            "rho_balance", "rho_scale", "initial_factors", "data_backend",
        }
        return fit_factorized_admm(counts, compatibility, **{k: v for k, v in kwargs.items() if k in allowed})
    if method == "admm":
        allowed = {"regularization", "rho", "max_iter", "inner_iter", "tol", "max_dense_entries", "device", "batch_cells", "adaptive_rho", "rho_update_interval", "rho_balance", "rho_scale"}
        return fit_admm(counts, compatibility, **{k: v for k, v in kwargs.items() if k in allowed})
    if method == "frank_wolfe":
        allowed = {"rank", "max_atoms", "tau", "max_iter", "power_iter", "tol", "device", "batch_cells", "seed", "min_iter", "patience"}
        return fit_frank_wolfe(counts, compatibility, **{k: v for k, v in kwargs.items() if k in allowed})
    if method == "frank_wolfe_penalized":
        allowed = {
            "rank", "max_atoms", "tau", "nonnegative_penalty", "max_iter",
            "power_iter", "tol", "device", "batch_cells", "seed",
            "min_iter", "patience", "initial_factors",
        }
        return fit_frank_wolfe_penalized(
            counts, compatibility, **{k: v for k, v in kwargs.items() if k in allowed}
        )
    raise ValueError(f"Unsupported Torch GLM method: {method}")


def result_to_csr(result: GLMResult, start: int, stop: int, threshold=0.0):
    """Materialize a normalized cell block as SciPy CSR for output only."""
    torch = _torch()
    if result.dense_phi is not None:
        phi = result.dense_phi[:, start:stop].T
    else:
        phi = result.right[start:stop] @ result.left.T
    phi = torch.relu(phi)
    phi /= phi.sum(dim=1, keepdim=True).clamp_min(1e-12)
    values = phi.detach().cpu().numpy()
    if threshold > 0:
        values[values < threshold] = 0
    return sp.csr_matrix(values)


def _factor_diagnostics(result: GLMResult):
    """Summarize whether a fitted representation is finite and nondegenerate."""
    if result.right is None:
        return {}
    right = result.right.detach()
    row_norm = right.square().sum(dim=1).sqrt()
    return {
        "left_finite": bool(result.left.isfinite().all().item()),
        "right_finite": bool(right.isfinite().all().item()),
        "left_norm": float(result.left.square().sum().sqrt().item()),
        "right_norm": float(right.square().sum().sqrt().item()),
        "active_cell_fraction": float((row_norm > 0).float().mean().item()),
    }


def factor_profile_diagnostics(result: GLMResult, *, max_features=512,
                               batch_cells=4096, target_sum=10_000.0):
    """Measure variation after rectification and per-cell profile normalization."""
    if result.left is None or result.right is None:
        return {}
    torch = _torch()
    left = result.left.detach()
    right = result.right.detach()
    n_features = min(max(1, int(max_features)), left.shape[0])
    loading_norms = left.square().sum(dim=1)
    selected = torch.topk(loading_norms, n_features).indices
    selected_left = left.index_select(0, selected)
    sums = np.zeros(n_features, dtype=np.float64)
    squared_sums = np.zeros(n_features, dtype=np.float64)
    active_count = 0
    for start in range(0, right.shape[0], max(1, int(batch_cells))):
        stop = min(start + max(1, int(batch_cells)), right.shape[0])
        abundance = torch.relu(right[start:stop] @ selected_left.T)
        totals = abundance.sum(dim=1)
        active = torch.isfinite(totals) & (totals > 0)
        if not bool(active.any().item()):
            continue
        profile = abundance[active]
        profile = torch.log1p(
            profile * (float(target_sum) / totals[active])[:, None]
        )
        sums += profile.sum(dim=0).detach().cpu().numpy().astype(np.float64)
        squared_sums += (
            profile.square().sum(dim=0).detach().cpu().numpy().astype(np.float64)
        )
        active_count += int(active.sum().item())
    if active_count < 2:
        relative_variance = 0.0
    else:
        means = sums / active_count
        variances = np.maximum(
            (squared_sums - active_count * np.square(means)) / (active_count - 1),
            0.0,
        )
        relative_variance = float(
            variances.sum() / max(np.square(means).sum(), np.finfo(float).tiny)
        )
    return {
        "normalized_profile_active_fraction": float(
            active_count / max(right.shape[0], 1)
        ),
        "normalized_profile_relative_variance": relative_variance,
        "normalized_profile_features": int(n_features),
    }


def write_chunked_result(result: GLMResult, output_prefix, barcodes, features,
                         batch_cells=4096, threshold=1e-8, write_chunks=True):
    """Write compact factors first and optional transcript estimates in blocks."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    n_cells = len(barcodes)
    np.savetxt(f"{output_prefix}glm_rows.txt", np.asarray(barcodes), fmt="%s")
    np.savetxt(f"{output_prefix}glm_cols.txt", np.asarray(features), fmt="%s")
    if result.left is not None:
        np.savez_compressed(
            f"{output_prefix}glm_factors.npz",
            left=result.left.detach().cpu().numpy(),
            right=result.right.detach().cpu().numpy(),
        )
    result.diagnostics.update(_factor_diagnostics(result))
    result.diagnostics.update(
        factor_profile_diagnostics(
            result, batch_cells=batch_cells
        )
    )
    chunk_paths = []
    if write_chunks:
        for index, start in enumerate(range(0, n_cells, int(batch_cells))):
            stop = min(start + int(batch_cells), n_cells)
            filename = f"{output_prefix.name}glm_cells_{index:05d}.npz"
            path = output_prefix.parent / filename
            sp.save_npz(path, result_to_csr(result, start, stop, threshold))
            chunk_paths.append({"path": filename, "start": start, "stop": stop})
    with open(f"{output_prefix}glm_manifest.json", "w") as handle:
        json.dump(
            {
                "method": result.method,
                "chunks_written": bool(write_chunks),
                "chunks": chunk_paths,
                "diagnostics": result.diagnostics,
            },
            handle,
            indent=2,
        )
