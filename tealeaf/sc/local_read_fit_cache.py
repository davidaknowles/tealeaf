"""Exact local-count fit reuse, no approximation or hypothesis reduction."""

from collections import OrderedDict
from copy import deepcopy

import numpy as np


class ExactLocalReadFitCache:
    """Bounded worker-local cache under one immutable supplied fit recipe.

    Keys contain the complete int64 tensor and its shape, not its significance,
    direction, gene ID or LR outcome. All independent requested records remain.
    The caller must create a new cache when changing any model/design setting.
    """

    def __init__(self, fit_function, max_entries=512):
        if not isinstance(max_entries, int) or max_entries < 1:
            raise ValueError('positive finite cache capacity required')
        self.fit_function = fit_function
        self.max_entries = max_entries
        self.entries = OrderedDict()
        self.evaluations = 0
        self.hits = 0

    def evaluate(self, counts):
        values = np.asarray(counts)
        if values.ndim != 4 or values.shape[2:] != (2, 2) or values.dtype != np.dtype(np.int64) or (values < 0).any() or (values > 2 ** 50).any():
            raise ValueError('original exact M/P/two-type/two-class int64 counts required')
        key = values.shape, values.tobytes()
        if key in self.entries:
            self.hits += 1
            value = self.entries.pop(key)
            self.entries[key] = value
            return *deepcopy(value), True
        self.evaluations += 1
        try:
            value = self.fit_function(values), ''
        except (ValueError, np.linalg.LinAlgError) as exc:
            value = None, str(exc)
        self.entries[key] = deepcopy(value)
        if len(self.entries) > self.max_entries:
            self.entries.popitem(last=False)
        return *value, False
