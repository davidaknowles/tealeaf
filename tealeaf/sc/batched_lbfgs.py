"""Independent bounded L-BFGS fits evaluated in one vectorized batch."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BatchedLBFGSResult:
    """Result of separable batched optimization.

    Arrays retain a leading fit dimension. ``evaluations`` counts vectorized
    objective calls, while ``line_search_evaluations`` is recorded per fit.
    """

    parameters: np.ndarray
    values: np.ndarray
    gradients: np.ndarray
    converged: np.ndarray
    failed: np.ndarray
    iterations: np.ndarray
    line_search_evaluations: np.ndarray
    evaluations: int
    scales: np.ndarray


def _projected_gradient(parameters, gradients, lower, upper):
    projected = gradients.copy()
    projected[(parameters <= lower) & (gradients > 0.0)] = 0.0
    projected[(parameters >= upper) & (gradients < 0.0)] = 0.0
    return projected


def _two_loop_direction(gradient, steps, differences):
    """Apply one fit's L-BFGS inverse-Hessian approximation."""
    if not steps:
        return -gradient
    vector = gradient.copy()
    coefficients = []
    for step, difference in zip(reversed(steps), reversed(differences)):
        inverse_curvature = 1.0 / np.dot(difference, step)
        coefficient = inverse_curvature * np.dot(step, vector)
        coefficients.append(coefficient)
        vector -= coefficient * difference
    last_step = steps[-1]
    last_difference = differences[-1]
    scale = np.dot(last_step, last_difference) / np.dot(
        last_difference, last_difference
    )
    vector *= scale
    for (step, difference), coefficient in zip(
        zip(steps, differences), reversed(coefficients)
    ):
        inverse_curvature = 1.0 / np.dot(difference, step)
        beta = inverse_curvature * np.dot(difference, vector)
        vector += step * (coefficient - beta)
    return -vector


def minimize_batched_lbfgs(
    value_and_gradient,
    initial,
    bounds,
    *,
    max_iter=200,
    history_size=10,
    gradient_tolerance=1e-5,
    maximum_line_search_steps=30,
    armijo=1e-4,
    backtrack=0.5,
    minimum_step=1e-12,
):
    """Minimize independent objectives with batched evaluations.

    ``value_and_gradient(parameters)`` must return one value per fit and a
    gradient with the same shape as ``parameters``. Each fit has independent
    L-BFGS history, active state, and projected Armijo backtracking search.
    """
    parameters = np.asarray(initial, dtype=float).copy()
    if parameters.ndim != 2:
        raise ValueError("initial must have shape (fits, parameters)")
    fit_count, parameter_count = parameters.shape
    if len(bounds) != parameter_count:
        raise ValueError("bounds must contain one pair per parameter")
    lower = np.broadcast_to(
        np.asarray([-np.inf if bound[0] is None else bound[0] for bound in bounds]),
        parameters.shape,
    )
    upper = np.broadcast_to(
        np.asarray([np.inf if bound[1] is None else bound[1] for bound in bounds]),
        parameters.shape,
    )
    parameters = np.clip(parameters, lower, upper)
    values, gradients = value_and_gradient(parameters)
    values = np.asarray(values, dtype=float).copy()
    gradients = np.asarray(gradients, dtype=float).copy()
    if values.shape != (fit_count,) or gradients.shape != parameters.shape:
        raise ValueError("objective returned incompatible value or gradient shapes")
    evaluations = 1
    scales = np.maximum(1.0, np.abs(values))
    scaled_values = values / scales
    scaled_gradients = gradients / scales[:, None]
    histories = [[] for _ in range(fit_count)]
    differences = [[] for _ in range(fit_count)]
    converged = np.zeros(fit_count, dtype=bool)
    failed = ~(np.isfinite(values) & np.isfinite(gradients).all(axis=1))
    iterations = np.zeros(fit_count, dtype=int)
    line_search_evaluations = np.zeros(fit_count, dtype=int)
    for _ in range(int(max_iter)):
        projected = _projected_gradient(parameters, scaled_gradients, lower, upper)
        gradient_norm = np.max(np.abs(projected), axis=1)
        converged |= gradient_norm <= float(gradient_tolerance)
        active = ~(converged | failed)
        if not np.any(active):
            break
        directions = np.zeros_like(parameters)
        for fit in np.flatnonzero(active):
            directions[fit] = _two_loop_direction(
                scaled_gradients[fit], histories[fit], differences[fit]
            )
            at_lower = (parameters[fit] <= lower[fit]) & (directions[fit] < 0.0)
            at_upper = (parameters[fit] >= upper[fit]) & (directions[fit] > 0.0)
            directions[fit, at_lower | at_upper] = 0.0
            if np.dot(scaled_gradients[fit], directions[fit]) >= 0.0:
                directions[fit] = -projected[fit]
        step_lengths = np.ones(fit_count, dtype=float)
        pending = active.copy()
        accepted = np.zeros(fit_count, dtype=bool)
        trial_parameters = parameters.copy()
        trial_values = values.copy()
        trial_gradients = gradients.copy()
        for _ in range(int(maximum_line_search_steps)):
            if not np.any(pending):
                break
            proposed = np.clip(
                parameters + step_lengths[:, None] * directions,
                lower,
                upper,
            )
            evaluation_parameters = np.where(pending[:, None], proposed, parameters)
            proposed_values, proposed_gradients = value_and_gradient(
                evaluation_parameters
            )
            evaluations += 1
            proposed_values = np.asarray(proposed_values, dtype=float)
            proposed_gradients = np.asarray(proposed_gradients, dtype=float)
            line_search_evaluations[pending] += 1
            displacement = proposed - parameters
            directional_gain = np.sum(scaled_gradients * displacement, axis=1)
            finite = np.isfinite(proposed_values) & np.isfinite(proposed_gradients).all(
                axis=1
            )
            armijo_ok = (
                proposed_values / scales
                <= scaled_values + float(armijo) * directional_gain
            )
            newly_accepted = pending & finite & armijo_ok
            trial_parameters[newly_accepted] = proposed[newly_accepted]
            trial_values[newly_accepted] = proposed_values[newly_accepted]
            trial_gradients[newly_accepted] = proposed_gradients[newly_accepted]
            accepted |= newly_accepted
            pending &= ~newly_accepted
            step_lengths[pending] *= float(backtrack)
            pending &= step_lengths >= float(minimum_step)
        failed |= active & ~accepted
        for fit in np.flatnonzero(accepted):
            step = trial_parameters[fit] - parameters[fit]
            difference = (trial_gradients[fit] - gradients[fit]) / scales[fit]
            curvature = np.dot(step, difference)
            threshold = 1e-12 * np.linalg.norm(step) * np.linalg.norm(difference)
            if curvature > threshold:
                histories[fit].append(step.copy())
                differences[fit].append(difference.copy())
                if len(histories[fit]) > int(history_size):
                    histories[fit].pop(0)
                    differences[fit].pop(0)
            iterations[fit] += 1
        parameters[accepted] = trial_parameters[accepted]
        values[accepted] = trial_values[accepted]
        gradients[accepted] = trial_gradients[accepted]
        scaled_values[accepted] = values[accepted] / scales[accepted]
        scaled_gradients[accepted] = gradients[accepted] / scales[accepted, None]
    projected = _projected_gradient(parameters, scaled_gradients, lower, upper)
    converged |= np.max(np.abs(projected), axis=1) <= float(gradient_tolerance)
    return BatchedLBFGSResult(
        parameters=parameters,
        values=values,
        gradients=gradients,
        converged=converged,
        failed=failed,
        iterations=iterations,
        line_search_evaluations=line_search_evaluations,
        evaluations=evaluations,
        scales=scales,
    )
