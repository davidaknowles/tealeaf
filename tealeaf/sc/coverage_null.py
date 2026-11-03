"""Coverage-only null filter applied to every method in the evaluation.

Each method is rerun on data in which every sample keeps its coverage of the tested feature and the two cell types
share one usage (drawn several times). A test is coverage-driven when the median of its null p-values (columns
p_null0, p_null1, ...) is below `level`. Only tests with a small observed p-value are simulated (targeted null), as only
they can change discoveries or replication; the others are kept.
"""

from __future__ import annotations

import pandas as pd


def null_median(table: pd.DataFrame) -> pd.Series:
    """Median null p-value over the draws (columns p_null<k>)."""
    columns = [column for column in table if column.startswith("p_null")]
    if not columns:
        raise ValueError("no p_null columns")
    return table[columns].median(axis=1)


def coverage_flagged(table: pd.DataFrame, keys, level: float = .05, p_column: str | None = None, threshold: float | None = None) -> pd.DataFrame:
    """Keys of the tests the coverage-only null calls (median null p < level).

    With p_column and threshold, only tests with observed p below threshold count as simulated (the targeted rule
    when the null was computed for more tests).
    """
    keys = list(keys)
    flagged = null_median(table) < level
    if threshold is not None:
        flagged &= pd.to_numeric(table[p_column], errors="coerce") < threshold
    return table.loc[flagged, keys].drop_duplicates().reset_index(drop=True)


def drop_flagged(tests: pd.DataFrame, flagged: pd.DataFrame, keys) -> pd.DataFrame:
    """Rows of tests whose keys are not in flagged."""
    keys = list(keys)
    marked = tests.merge(flagged[keys].drop_duplicates().assign(_flagged=True), on=keys, how="left")
    return marked.loc[marked._flagged.isna()].drop(columns="_flagged").reset_index(drop=True)
