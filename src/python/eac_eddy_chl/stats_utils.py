"""Statistical helper functions for robust eddy diagnostics."""

from __future__ import annotations

import numpy as np
from scipy.stats import pearsonr, spearmanr, theilslopes, wilcoxon


def safe_wilcoxon(x, min_n=5, alpha=0.05):
    """One-sample Wilcoxon signed-rank test against zero; return p-value or NaN."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < min_n or np.allclose(x, 0.0):
        return np.nan
    try:
        res = wilcoxon(x, y=None, alternative="two-sided", zero_method="wilcox", method="auto", nan_policy="omit")
        return res.pvalue
    except Exception:
        return np.nan


def safe_theilsen(x, y, min_n=5, alpha=0.95, method="separate"):
    """Robust Theil-Sen slope fit for y versus x; return None if invalid."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < min_n or np.nanstd(x) == 0:
        return None
    try:
        return theilslopes(y, x, alpha=alpha, method=method)
    except Exception:
        return None


def safe_spearman(x, y, min_n=5):
    """Spearman correlation with NaN/constant safeguards."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < min_n or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return None
    try:
        return spearmanr(x, y, nan_policy="omit", alternative="two-sided")
    except Exception:
        return None


def safe_pearson(x, y, min_n=5):
    """Pearson correlation with NaN/constant safeguards."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < min_n or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return None
    try:
        return pearsonr(x, y, alternative="two-sided")
    except Exception:
        return None


def median_iqr(x, min_n=5):
    """Return median, 25th percentile, 75th percentile, and sample size."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < min_n:
        return np.nan, np.nan, np.nan, x.size
    med = np.nanmedian(x)
    q25 = np.nanpercentile(x, 25)
    q75 = np.nanpercentile(x, 75)
    return med, q25, q75, x.size
