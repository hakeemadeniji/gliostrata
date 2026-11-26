"""Statistical helpers for GlioStrata.

Bootstrap confidence intervals, RMST (restricted mean survival time) helpers,
concordance correlation coefficient, and multiple-comparison correction.

These wrap standard formulas -- kept here so all tests use the same implementation.
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
    per the defensibility stack -- point estimate alone is not enough.

    Method: the PERCENTILE bootstrap. ``n_resamples`` resamples of size ``n``
    are drawn with replacement from ``data``, ``statistic`` is applied to each,
    and the bounds are the ``alpha/2`` and ``1 - alpha/2`` empirical quantiles
    where ``alpha = 1 - confidence``. This is NOT the bias-corrected and
    accelerated (BCa) interval and carries no bias/skew correction; for a
    strongly skewed statistic the interval may be shifted relative to BCa.

    Reproducible: the same ``seed`` always yields the same interval.

    Edge cases: an empty input returns ``(nan, nan, nan)``. NaNs are not
    removed -- they propagate through the default ``np.mean``, so a NaN in the
    input surfaces as a NaN interval rather than being silently dropped.
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

    CCC = 2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y)) ** 2)

    Measures agreement, not just linear correlation -- catches systematic bias
    that Pearson r misses. Used for T7 (predicted-vs-actual tumor kill).

    All moments use the population normalisation (ddof=0), consistently
    between numerator and denominator, as in Lin (1989). Mixing ddof=1 in the
    variances with ddof=0 in the covariance (or vice versa) biases the result.

    Returns NaN for mismatched lengths or fewer than 2 points, and 0.0 when the
    denominator vanishes (both inputs constant and equal, i.e. no variance to
    agree about).
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
    after multiple-testing correction, in the ORIGINAL input order.

    Step-up procedure (Benjamini & Hochberg, 1995): sort the m p-values
    ascending, find the LARGEST rank k with p_(k) <= (k / m) * q, and reject
    every hypothesis of rank 1..k -- including ones whose own p-value exceeds
    its own (rank / m) * q threshold. If no rank satisfies the inequality,
    nothing is rejected.
    """
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda pair: pair[1])
    # Largest rank k (1-based) satisfying p_(k) <= (k / n) * q.
    k_max = 0
    for rank, (_orig_idx, p) in enumerate(indexed, start=1):
        if p <= q * rank / n:
            k_max = rank
    significant = [False] * n
    for rank, (orig_idx, _p) in enumerate(indexed, start=1):
        if rank <= k_max:
            significant[orig_idx] = True
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
    More robust than median survival -- always defined, tighter CIs.

    `event` is 1 for an observed event (death/progression) and 0 for a
    right-censored observation. Censoring is handled properly: the
    Kaplan-Meier estimator only drops at observed event times, while censored
    subjects still count in the number at risk up to and including their
    censoring time (the standard convention for ties between an event and a
    censoring at the same time).

    The area is the exact integral of the KM step function, accumulated as
    S(t_{j-1}) * (t_j - t_{j-1}) -- the survival probability BEFORE the drop at
    t_j applies over the interval leading up to t_j. No trapezoid rule: that
    would systematically bias a step function.

    If `horizon` lies beyond the longest follow-up the last KM value is carried
    forward to the horizon (the usual convention; the estimate is then
    extrapolated and its variance understated). An arm with no observed events
    at or before the horizon has RMST == horizon. An empty arm gives NaN.

    Pure numpy -- no lifelines / scikit-survival dependency.
    """
    def _km_rmst(time, event, h):
        if len(time) == 0:
            return float("nan")
        order = np.argsort(time, kind="mergesort")
        t = time[order]
        e = event[order]
        # KM drops only at observed event times at or before the horizon.
        event_times = np.unique(t[(e == 1) & (t <= h)])
        survival = 1.0
        prev_t = 0.0
        area = 0.0
        for ut in event_times:
            n_at_risk = np.sum(t >= ut)
            n_events = np.sum((t == ut) & (e == 1))
            # Area of the interval (prev_t, ut] uses survival BEFORE this drop.
            area += survival * (ut - prev_t)
            if n_at_risk > 0:
                survival *= 1.0 - n_events / n_at_risk
            prev_t = ut
        # Tail from the last event time (or 0) out to the horizon.
        area += survival * (h - prev_t)
        return float(area)

    rmst_a = _km_rmst(np.asarray(time_a, dtype=float), np.asarray(event_a, dtype=float), horizon)
    rmst_b = _km_rmst(np.asarray(time_b, dtype=float), np.asarray(event_b, dtype=float), horizon)
    return float(rmst_b - rmst_a)
