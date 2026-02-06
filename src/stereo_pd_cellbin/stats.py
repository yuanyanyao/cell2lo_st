from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency, mannwhitneyu


def p_to_stars(p: float) -> str:
    """
    Convert a p-value into a star label used in figures.

    Convention:
    - ns:   p >= 0.05
    - *:    p <  0.05
    - **:   p <  0.01
    - ***:  p <  0.001
    - ****: p <  1e-4
    """
    if p is None or np.isnan(p):
        return "na"
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def fisher_2x2(a: int, b: int, c: int, d: int) -> Tuple[float, float]:
    """
    Fisher's exact test for a 2x2 table:
        [[a, b],
         [c, d]]
    Returns (odds_ratio, p_value).
    """
    odds, p = fisher_exact([[a, b], [c, d]])
    return float(odds), float(p)


def chi2_test(table: np.ndarray) -> Tuple[float, float, int]:
    """
    Chi-square test for an RxC contingency table.
    Returns (chi2, p_value, dof).
    """
    chi2, p, dof, _ = chi2_contingency(table)
    return float(chi2), float(p), int(dof)


def mwu(x: np.ndarray, y: np.ndarray, alternative: str = "two-sided") -> float:
    """
    Mann–Whitney U test p-value for comparing two distributions.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size == 0 or y.size == 0:
        return np.nan
    _, p = mannwhitneyu(x, y, alternative=alternative)
    return float(p)
