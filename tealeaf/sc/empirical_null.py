"""Apply frozen empirical sign-null pools to independently generated trials."""

import numpy as np


def leave_parent_out_cdf(probabilities, parents, strata, training_probabilities, training_parents, training_strata):
    """Plus-one empirical CDF with all draws of the target parent excluded.

    Target arrays are aligned N-vectors, training arrays aligned K-vectors.
    Repeated target parents are permitted (e.g. distinct actual count draws).
    Training strata stay frozen, rather than changing with a simulated test's
    geometry. An absent stratum gives p1 and zero remaining training draws,
    explicitly identifying unavailable tail calibration. Returns two
    N-vectors, calibrated probabilities and leave-parent-out pool sizes.
    """
    probabilities, training_probabilities = (np.asarray(value, dtype=float) for value in (probabilities, training_probabilities))
    parents, strata, training_parents, training_strata = (np.asarray(value, dtype=object) for value in (parents, strata, training_parents, training_strata))
    if probabilities.ndim != 1 or training_probabilities.ndim != 1 or parents.shape != probabilities.shape or strata.shape != probabilities.shape or training_parents.shape != training_probabilities.shape or training_strata.shape != training_probabilities.shape:
        raise ValueError("aligned one-dimensional target and training arrays required")
    if any(not np.isfinite(value).all() or (value < 0).any() or (value > 1).any() for value in (probabilities, training_probabilities)):
        raise ValueError("finite probabilities in the unit interval required")
    if any(any(not isinstance(value, str) or not value for value in labels) for labels in (parents, strata, training_parents, training_strata)):
        raise ValueError("nonempty string parents and strata required")
    calibrated, counts = np.ones(len(probabilities)), np.zeros(len(probabilities), dtype=int)
    for stratum in np.unique(strata):
        positions = np.flatnonzero(strata == stratum)
        selected = training_strata == stratum
        values, pool_parents = training_probabilities[selected], training_parents[selected]
        pool = np.sort(values)
        wanted = np.unique(parents[positions])
        own_mask = np.isin(pool_parents, wanted)
        own_parents, own_probabilities = pool_parents[own_mask], values[own_mask]
        own = {parent: np.sort(own_probabilities[own_parents == parent]) for parent in wanted}
        for position in positions:
            excluded = own[parents[position]]
            counts[position] = len(pool) - len(excluded)
            numerator = 1 + np.searchsorted(pool, probabilities[position], side="right") - np.searchsorted(excluded, probabilities[position], side="right")
            calibrated[position] = numerator / (1 + counts[position])
    return calibrated, counts
