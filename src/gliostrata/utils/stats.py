"""Statistical helpers for GlioStrata.

Bootstrap confidence intervals, RMST (restricted mean survival time) helpers,
concordance correlation coefficient, and multiple-comparison correction.

These wrap standard formulas — kept here so all tests use the same implementation.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def bootstrap_ci(
    data: np.ndarray | list[float],
    statistic=np.mean,
    n_resamples: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Bootstrap confidence interval for any statistic.

    Returns (point_estimate, lower_bound, upper_bound).

    The lower_bound is what pre-registered thresholds are evaluated against
    per the defensibility stack — point estimate alone is not enough.
    """
    rng = np.random.default_rng(seed)
    arr = np.asarray(data)
    n = len(arr)
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    point = float(statistic(arr))
    resamples = np.empty(n_resamples)
    for i in range(n_resamples):
        sample = rng.choice(arr, size=n, replace=True)
        resamples[i] = statistic(sample)
    alpha = (1 - confidence) / 2
    lower = float(np.quantile(resamples, alpha))
    upper = float(np.quantile(resamples, 1 - alpha))
    return (point, lower, upper)


def concordance_correlation_coefficient(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's concordance correlation coefficient (CCC).

    Measures agreement, not just linear correlation — catches systematic bias
    that Pearson r misses. Used for T7 (predicted-vs-actual tumor kill).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mean_x, mean_y = x.mean(), y.mean()
    var_x, var_y = x.var(), y.var()
    covariance = ((x - mean_x) * (y - mean_y)).mean()
    denom = var_x + var_y + (mean_x - mean_y) ** 2
    if denom == 0:
        return 0.0
    return float(2 * covariance / denom)


def benjamini_hochberg(p_values: list[float], q: float = 0.05) -> list[bool]:
    """Benjamini-Hochberg FDR correction.

    Returns a list of booleans indicating whether each p-value is significant
    after multiple-testing correction.
    """
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    significant = [False] * n
    running_max = 0.0
    for rank, (orig_idx, p) in enumerate(indexed, start=1):
        threshold = q * rank / n
        if p <= threshold:
            significant[orig_idx] = True
            running_max = max(running_max, p)
    # Enforce monotonicity: once significant at any rank, all lower ranks are too
    cumulative = False
    for rank, (orig_idx, p) in enumerate(indexed, start=1):
        if significant[orig_idx]:
            cumulative = True
        significant[orig_idx] = cumulative
    return significant


def rmst_difference(
    time_a: np.ndarray,
    event_a: np.ndarray,
    time_b: np.ndarray,
    event_b: np.ndarray,
    horizon: float,
) -> float:
    """Restricted Mean Survival Time difference (B - A) up to a horizon.

    RMST is the area under the survival curve truncated at `horizon`.
    More robust than median survival — always defined, tighter CIs.

    Uses the Kaplan-Meier estimator (simple implementation).
    For production use, prefer lifelines.RMST or survival::survRM2 (R).
    """
    def _km_rmst(time, event, h):
        # Sort by time
        order = np.argsort(time)
        t = time[order]
        e = event[order]
        # Kaplan-Meier
        unique_t = np.unique(t[t <= h])
        if len(unique_t) == 0:
            return h  # no events before horizon → all survive
        survival = 1.0
        prev_t = 0.0
        area = 0.0
        for ut in unique_t:
            mask = (t == ut)
            n_at_risk = np.sum(t >= ut)
            n_events = np.sum(mask & (e == 1))
            if n_at_risk > 0:
                survival *= (1 - n_events / n_at_risk)
            area += survival * (ut - prev_t)
            prev_t = ut
        # Add tail from last event time to horizon
        area += survival * (h - prev_t)
        return area

    rmst_a = _km_rmst(np.asarray(time_a, dtype=float), np.asarray(event_a, dtype=float), horizon)
    rmst_b = _km_rmst(np.asarray(time_b, dtype=float), np.asarray(event_b, dtype=float), horizon)
    return float(rmst_b - rmst_a)
