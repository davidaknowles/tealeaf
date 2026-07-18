"""Memory-bounded genome-wide GLM solvers for single-cell EC counts.

The scalable solvers keep a transcript-by-rank factor and a cell-by-rank
factor.  Count blocks stay sparse, so neither ``A.T @ C`` nor a full
transcript-by-cell abundance matrix is materialized during fitting.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
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
        torch.as_tensor(matrix.indptr, dtype=torch.int64, device=device),
        torch.as_tensor(matrix.indices, dtype=torch.int64, device=device),
        torch.as_tensor(matrix.data, dtype=dtype, device=device),
        size=matrix.shape,
        device=device,
    )


def _normalized_block(counts: sp.csr_matrix, start: int, stop: int, device):
    """Return a normalized cells-by-EC sparse COO block and nonempty mask."""
    torch = _torch()
    block = counts[start:stop].tocsr(copy=True).astype(np.float32)
    totals = np.asarray(block.sum(axis=1)).ravel()
    nonempty = totals > 0
    inverse = np.zeros_like(totals, dtype=np.float32)
    inverse[nonempty] = 1.0 / totals[nonempty]
    block.data *= np.repeat(inverse, np.diff(block.indptr))
    coo = block.tocoo()
    sparse = torch.sparse_coo_tensor(
        torch.as_tensor(np.vstack((coo.row, coo.col)), dtype=torch.int64, device=device),
        torch.as_tensor(coo.data, dtype=torch.float32, device=device),
        size=block.shape,
        device=device,
    ).coalesce()
    return sparse, torch.as_tensor(nonempty, device=device)


@dataclass
class GLMResult:
    method: str
    left: object | None
    right: object | None
    dense_phi: object | None
    diagnostics: dict


class SparseGLM:
    """Sparse count access and low-rank objective helpers."""

    def __init__(self, counts, compatibility, device="auto", batch_cells=4096):
        self.torch = _torch()
        self.device = resolve_device(device)
        self.counts = counts.tocsr()
        self.a = _csr_to_torch(compatibility, self.device)
        self.a_t = self.a.transpose(0, 1).to_sparse_coo().coalesce()
        self.n_cells, self.n_ec = self.counts.shape
        self.n_transcripts = compatibility.shape[1]
        self.batch_cells = max(1, int(batch_cells))
        self.n_nonempty = max(1, int(np.count_nonzero(np.asarray(self.counts.sum(axis=1)).ravel())))
        self.gram_lipschitz = self._estimate_gram_lipschitz()
        self.total_count_sq = float(np.square(self.counts.data / np.repeat(
            np.maximum(np.asarray(self.counts.sum(axis=1)).ravel(), 1.0),
            np.diff(self.counts.indptr),
        )).sum())

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

    def blocks(self) -> Iterator[tuple[int, int, object, object]]:
        for start in range(0, self.n_cells, self.batch_cells):
            stop = min(start + self.batch_cells, self.n_cells)
            block, nonempty = _normalized_block(self.counts, start, stop, self.device)
            yield start, stop, block, nonempty

    def a_times(self, value):
        return self.torch.sparse.mm(self.a, value)

    def at_times(self, value):
        return self.torch.sparse.mm(self.a_t, value)

    def gram_times(self, value):
        return self.at_times(self.a_times(value))

    def b_times(self, right):
        """Compute ``A.T @ C @ right`` without forming ``C`` or ``A.T @ C``."""
        result = self.torch.zeros((self.n_transcripts, right.shape[1]), device=self.device)
        for start, stop, block, _ in self.blocks():
            result += self.at_times(self.torch.sparse.mm(block.transpose(0, 1), right[start:stop]))
        return result

    def loss_for_factors(self, left, right, regularization=0.0):
        au = self.a_times(left)
        cross = self.torch.zeros((), device=self.device)
        right_gram = self.torch.zeros((right.shape[1], right.shape[1]), device=self.device)
        for start, stop, block, nonempty in self.blocks():
            right_block = right[start:stop] * nonempty[:, None]
            cross += (self.torch.sparse.mm(block, au) * right_block).sum()
            right_gram += right_block.T @ right_block
        left_gram = left.T @ self.gram_times(left)
        loss = 0.5 * (self.total_count_sq - 2.0 * cross + (left_gram * right_gram).sum())
        if regularization:
            loss += 0.5 * regularization * ((left * left).sum() + (right * right).sum())
        return float(loss.detach().cpu())


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


def fit_factorized(counts, compatibility, *, rank=64, regularization=0.01,
                   max_iter=100, tol=1e-4, learning_rate=0.05, device="auto",
                   batch_cells=4096, seed=0, min_iter=10, patience=5):
    """Projected alternating optimization of a genome-wide nonnegative factorization."""
    _validate_stopping(max_iter, min_iter, patience, tol)
    data = SparseGLM(counts, compatibility, device, batch_cells)
    torch = data.torch
    left, right = _initial_factors(data, rank, seed)
    diagnostics = {
        "objective": [], "objective_relative_change": [],
        "device": str(data.device), "rank": rank,
        "regularization": float(regularization),
        "gram_lipschitz": data.gram_lipschitz,
        "min_iter": int(min_iter), "patience": int(patience), "tolerance": tol,
    }
    previous = None
    stable = 0

    for iteration in range(int(max_iter)):
        au = data.a_times(left)
        gram_left = left.T @ data.gram_times(left)
        right_lipschitz = float(
            torch.linalg.matrix_norm(gram_left, ord=2).item()
        ) + regularization
        right_step = min(learning_rate, 1.0 / max(right_lipschitz, 1e-12))
        b_times_right = torch.zeros_like(left)
        right_gram = torch.zeros((rank, rank), device=data.device)
        for start, stop, block, nonempty in data.blocks():
            right_block = right[start:stop]
            gradient = right_block @ gram_left - torch.sparse.mm(block, au)
            gradient += regularization * right_block
            right_block = torch.relu(right_block - right_step * gradient)
            right_block[~nonempty] = 0
            right[start:stop] = right_block
            right_gram += right_block.T @ right_block
            b_times_right += data.at_times(
                torch.sparse.mm(block.transpose(0, 1), right_block)
            )
        # Average the complete block gradient. Averaging only the data term
        # makes the effective nuclear penalty depend on the number of cells.
        gradient_left = (
            data.gram_times(left) @ right_gram
            - b_times_right
            + regularization * left
        ) / data.n_nonempty
        left_lipschitz = (
            data.gram_lipschitz
            * float(torch.linalg.matrix_norm(right_gram, ord=2).item())
            / data.n_nonempty
            + regularization / data.n_nonempty
        )
        left_step = min(learning_rate, 1.0 / max(left_lipschitz, 1e-12))
        left = torch.relu(left - left_step * gradient_left)

        objective = data.loss_for_factors(left, right, regularization)
        diagnostics["objective"].append(objective)
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
        previous = objective
    else:
        diagnostics.update(
            iterations=int(max_iter), converged=False,
            convergence_reason="max_iterations",
        )
    return GLMResult("factorized", left, right, None, diagnostics)


def fit_factorized_admm(counts, compatibility, *, rank=64, regularization=0.01,
                        rho=1.0, max_iter=100, tol=1e-4, learning_rate=0.05,
                        device="auto", batch_cells=4096, seed=0,
                        min_iter=10, patience=5, adaptive_rho=True,
                        rho_update_interval=10, rho_balance=10.0,
                        rho_scale=2.0):
    """Factor-sized ADMM split for nonnegative low-rank genome-wide fitting.

    This is deliberately distinct from :func:`fit_admm`: it uses the variational
    factorization and therefore is non-convex, while retaining only ``O(r(T+M))``
    state.
    """
    _validate_stopping(max_iter, min_iter, patience, tol)
    if rho <= 0 or rho_update_interval < 1 or rho_balance <= 1 or rho_scale <= 1:
        raise ValueError("rho adaptation parameters are outside their valid range")
    rho = float(rho)
    data = SparseGLM(counts, compatibility, device, batch_cells)
    torch = data.torch
    left, right = _initial_factors(data, rank, seed)
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
    }
    previous = None
    stable = 0

    for iteration in range(int(max_iter)):
        au = data.a_times(left)
        gram_left = left.T @ data.gram_times(left)
        right_lipschitz = float(
            torch.linalg.matrix_norm(gram_left, ord=2).item()
        ) + regularization + rho
        right_step = min(learning_rate, 1.0 / max(right_lipschitz, 1e-12))
        b_times_right = torch.zeros_like(left)
        right_gram = torch.zeros((rank, rank), device=data.device)
        primal_sq = 0.0
        dual_sq = 0.0
        split_norm_sq = 0.0
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
            primal_sq += float(
                (right_block - right_copy[start:stop]).square().sum().item()
            )
            dual_sq += float(
                (right_copy[start:stop] - previous_right_copy).square().sum().item()
            )
            split_norm_sq += float(right_copy[start:stop].square().sum().item())
            right_gram += right_copy[start:stop].T @ right_copy[start:stop]
            b_times_right += data.at_times(
                torch.sparse.mm(block.transpose(0, 1), right_copy[start:stop])
            )
        # Keep the loss, regularization, and augmented-Lagrangian terms on the
        # same scale when stabilizing this block update by the cell count.
        gradient_left = (
            data.gram_times(left) @ right_gram
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
        primal_sq += float((left - left_copy).square().sum().item())
        dual_sq += float((left_copy - previous_left_copy).square().sum().item())
        split_norm_sq += float(left_copy.square().sum().item())
        primal_relative = np.sqrt(primal_sq) / max(np.sqrt(split_norm_sq), 1e-12)
        dual_relative = rho * np.sqrt(dual_sq) / max(
            rho * np.sqrt(split_norm_sq), 1e-12
        )
        left = left_copy
        right = right_copy

        objective = data.loss_for_factors(left, right, regularization)
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
            primal_residual = np.sqrt(primal_sq)
            dual_residual = rho * np.sqrt(dual_sq)
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
    max_iter = min(int(max_iter), max_atoms)
    tau = float(np.sqrt(data.n_cells) if tau is None else tau)
    penalty = float(nonnegative_penalty) * data.gram_lipschitz
    left = torch.empty((data.n_transcripts, 0), device=data.device)
    right = torch.empty((data.n_cells, 0), device=data.device)
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
    }
    initial_gap = None
    previous = None
    stable = 0

    for iteration in range(max_iter):
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
        if previous is None:
            relative_change = None
            stable = 0
        else:
            relative_change = abs(previous - objective) / max(abs(previous), 1.0)
            stable = stable + 1 if (
                relative_change <= tol and relative_gap <= tol
            ) else 0
        diagnostics["objective"].append(objective)
        diagnostics["objective_relative_change"].append(relative_change)
        diagnostics["candidate_gap"].append(gap)
        diagnostics["relative_candidate_gap"].append(relative_gap)
        diagnostics["line_search_step"].append(float(gamma.item()))
        diagnostics["oracle_singular_value"].append(float(singular.item()))
        diagnostics["oracle_relative_residual"].append(float(oracle_residual.item()))
        if iteration + 1 >= int(min_iter) and stable >= int(patience):
            diagnostics.update(
                iterations=iteration + 1, converged=True,
                convergence_reason="objective_and_candidate_gap_patience",
            )
            break
        previous = objective
    else:
        reason = "atom_capacity" if max_iter == max_atoms else "max_iterations"
        diagnostics.update(
            iterations=max_iter, converged=False, convergence_reason=reason,
        )
    diagnostics.update(_negative_factor_stats(data, left, right))
    return GLMResult("frank_wolfe_penalized", left, right, None, diagnostics)


def fit_glm(counts, compatibility, method, **kwargs):
    """Dispatch a GLM method using a cells-by-EC sparse count matrix."""
    if method == "factorized":
        allowed = {"rank", "regularization", "max_iter", "tol", "learning_rate", "device", "batch_cells", "seed", "min_iter", "patience"}
        return fit_factorized(counts, compatibility, **{k: v for k, v in kwargs.items() if k in allowed})
    if method == "admm_factorized":
        allowed = {"rank", "regularization", "rho", "max_iter", "tol", "learning_rate", "device", "batch_cells", "seed", "min_iter", "patience", "adaptive_rho", "rho_update_interval", "rho_balance", "rho_scale"}
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
            "min_iter", "patience",
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
