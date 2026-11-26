"""Tests for gliostrata.utils.stats.

This module carries every number in the pre-registered threshold test plan
(T1-T5): every reported metric needs a 1000-resample bootstrap 95% CI, PASS is
decided on the CI lower bound, multiplicity is controlled with
Benjamini-Hochberg, and the headline endpoint is an RMST difference. So the
tests below check each function against its mathematical definition using
hand-computed closed forms and independent reference implementations, not just
against itself.

Regression markers:
  * ``test_bh_does_not_reject_above_largest_passing_rank`` and
    ``test_bh_rejects_everything_below_largest_passing_rank`` pin the BH
    step-up direction (the old code propagated significance UPWARD in rank,
    inverting the procedure).
  * ``test_rmst_single_event_uses_survival_before_the_drop`` and
    ``test_km_rmst_matches_hand_worked_example`` pin the step-function
    integration (the old code applied the post-drop survival to the interval
    leading up to the drop, biasing RMST downward).
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import false_discovery_control

from gliostrata.utils.stats import (
    benjamini_hochberg,
    bootstrap_ci,
    concordance_correlation_coefficient,
    rmst_difference,
)

# --------------------------------------------------------------------------
# Independent reference implementations (deliberately written from the
# textbook definitions rather than by refactoring the code under test).
# --------------------------------------------------------------------------


def ref_ccc(x, y):
    """Lin's CCC from the definition, population moments throughout."""
    x = [float(v) for v in x]
    y = [float(v) for v in y]
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    vx = sum((a - mx) ** 2 for a in x) / n
    vy = sum((b - my) ** 2 for b in y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y)) / n
    return 2 * cov / (vx + vy + (mx - my) ** 2)


def ref_bh(p_values, q):
    """BH step-up via scipy's adjusted p-values (independent of our code)."""
    if len(p_values) == 0:
        return []
    adjusted = false_discovery_control(np.asarray(p_values, dtype=float), method="bh")
    return [bool(v) for v in (adjusted <= q)]


def ref_km_rmst(time, event, horizon):
    """Kaplan-Meier RMST written out longhand, independent of the code."""
    time = [float(t) for t in time]
    event = [float(e) for e in event]
    event_times = sorted({t for t, e in zip(time, event) if e == 1 and t <= horizon})
    surv = 1.0
    prev = 0.0
    area = 0.0
    for ut in event_times:
        at_risk = sum(1 for t in time if t >= ut)
        deaths = sum(1 for t, e in zip(time, event) if t == ut and e == 1)
        area += surv * (ut - prev)
        surv *= 1.0 - deaths / at_risk
        prev = ut
    return area + surv * (horizon - prev)


# ==========================================================================
# bootstrap_ci
# ==========================================================================


class TestBootstrapCI:
    def test_point_estimate_is_the_statistic_of_the_data(self):
        data = [1.0, 2.0, 3.0, 4.0, 5.0]
        point, _lo, _hi = bootstrap_ci(data)
        assert point == pytest.approx(3.0)

    def test_point_estimate_honours_custom_statistic(self):
        data = [1.0, 2.0, 3.0, 4.0, 100.0]
        point, _lo, _hi = bootstrap_ci(data, statistic=np.median)
        assert point == pytest.approx(3.0)
        point_max, _lo, _hi = bootstrap_ci(data, statistic=np.max)
        assert point_max == pytest.approx(100.0)

    def test_custom_statistic_drives_the_bounds_too(self):
        """A resample's statistic, not just the point estimate, must be the
        supplied callable: bounds for max must sit in the data's range."""
        data = list(range(1, 51))
        _p, lo, hi = bootstrap_ci(data, statistic=np.max, n_resamples=200)
        assert 1 <= lo <= hi <= 50
        assert hi == pytest.approx(50.0)  # max of a resample tops out at 50

    def test_ordering_and_containment_of_point_estimate(self):
        rng = np.random.default_rng(0)
        data = rng.normal(10.0, 2.0, size=60)
        point, lo, hi = bootstrap_ci(data)
        assert lo <= point <= hi

    def test_determinism_under_fixed_seed(self):
        data = list(range(40))
        a = bootstrap_ci(data, seed=7)
        b = bootstrap_ci(data, seed=7)
        assert a == b

    def test_different_seeds_give_different_intervals(self):
        data = list(range(40))
        a = bootstrap_ci(data, seed=1)
        b = bootstrap_ci(data, seed=2)
        assert a[0] == pytest.approx(b[0])  # point estimate is seed-free
        assert (a[1], a[2]) != (b[1], b[2])

    def test_wider_confidence_gives_wider_interval(self):
        rng = np.random.default_rng(3)
        data = rng.normal(0.0, 1.0, size=80)
        _p, lo95, hi95 = bootstrap_ci(data, confidence=0.95, seed=11)
        _p, lo99, hi99 = bootstrap_ci(data, confidence=0.99, seed=11)
        _p, lo80, hi80 = bootstrap_ci(data, confidence=0.80, seed=11)
        assert (hi99 - lo99) > (hi95 - lo95) > (hi80 - lo80)
        assert lo99 <= lo95 <= lo80
        assert hi99 >= hi95 >= hi80

    def test_bounds_are_the_requested_percentiles(self):
        """Reproduce the percentile method exactly: with the same seed and the
        same draw sequence, the bounds must be the alpha/2 and 1-alpha/2
        empirical quantiles of the resample distribution."""
        data = np.arange(1.0, 26.0)
        n_resamples, confidence, seed = 300, 0.90, 123
        rng = np.random.default_rng(seed)
        stats_ = np.array(
            [np.mean(rng.choice(data, size=len(data), replace=True)) for _ in range(n_resamples)]
        )
        alpha = (1 - confidence) / 2
        expected_lo = float(np.quantile(stats_, alpha))
        expected_hi = float(np.quantile(stats_, 1 - alpha))
        _p, lo, hi = bootstrap_ci(
            data, n_resamples=n_resamples, confidence=confidence, seed=seed
        )
        assert lo == pytest.approx(expected_lo)
        assert hi == pytest.approx(expected_hi)

    def test_resampling_is_with_replacement_and_size_n(self):
        """With replacement, a resample of a 2-element set can be degenerate,
        so the sample mean of [0, 1] must reach the extremes 0.0 and 1.0.
        Sampling without replacement would always give exactly 0.5."""
        _p, lo, hi = bootstrap_ci([0.0, 1.0], n_resamples=500, confidence=0.999, seed=5)
        assert lo == pytest.approx(0.0)
        assert hi == pytest.approx(1.0)

    def test_empty_input_returns_all_nan(self):
        point, lo, hi = bootstrap_ci([])
        assert math.isnan(point) and math.isnan(lo) and math.isnan(hi)

    def test_single_element_gives_degenerate_interval(self):
        point, lo, hi = bootstrap_ci([4.2], n_resamples=50)
        assert point == pytest.approx(4.2)
        assert lo == pytest.approx(4.2)
        assert hi == pytest.approx(4.2)

    def test_constant_input_gives_zero_width_interval(self):
        point, lo, hi = bootstrap_ci([3.0] * 25, n_resamples=100)
        assert (point, lo, hi) == pytest.approx((3.0, 3.0, 3.0))

    def test_nan_propagates_rather_than_being_silently_dropped(self):
        point, lo, hi = bootstrap_ci([1.0, 2.0, float("nan"), 4.0], n_resamples=50)
        assert math.isnan(point)
        assert math.isnan(lo) or math.isnan(hi)

    def test_nan_safe_statistic_can_be_supplied(self):
        point, lo, hi = bootstrap_ci(
            [1.0, 2.0, float("nan"), 3.0], statistic=np.nanmean, n_resamples=100, seed=9
        )
        assert point == pytest.approx(2.0)
        assert lo <= point <= hi

    def test_interval_covers_the_true_mean_for_a_large_sample(self):
        rng = np.random.default_rng(42)
        data = rng.normal(5.0, 1.0, size=500)
        _p, lo, hi = bootstrap_ci(data, n_resamples=1000, seed=42)
        assert lo < 5.0 < hi

    def test_accepts_numpy_array_and_list_identically(self):
        data = [2.0, 4.0, 6.0, 8.0]
        assert bootstrap_ci(data, seed=4) == bootstrap_ci(np.asarray(data), seed=4)

    def test_default_resample_count_is_the_preregistered_1000(self):
        """T1-T5 specify 1000 resamples; the default must not drift."""
        import inspect

        sig = inspect.signature(bootstrap_ci)
        assert sig.parameters["n_resamples"].default == 1000
        assert sig.parameters["confidence"].default == 0.95


# ==========================================================================
# concordance_correlation_coefficient
# ==========================================================================


class TestCCC:
    def test_hand_computed_value(self):
        """x=[1,2,3,4], y=[1,2,3,5].
        mean_x=2.5, mean_y=2.75; var_x=1.25, var_y=2.1875, cov=1.625
        denom = 1.25 + 2.1875 + 0.0625 = 3.5
        CCC = 2 * 1.625 / 3.5 = 3.25 / 3.5 = 0.9285714285714286
        """
        ccc = concordance_correlation_coefficient([1, 2, 3, 4], [1, 2, 3, 5])
        assert ccc == pytest.approx(0.9285714285714286, rel=1e-12)

    def test_identical_inputs_give_exactly_one(self):
        x = [1.0, 5.0, 2.0, 9.0, -3.0]
        assert concordance_correlation_coefficient(x, x) == pytest.approx(1.0)

    def test_sign_flipped_inputs_centred_give_minus_one(self):
        x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        assert concordance_correlation_coefficient(x, -x) == pytest.approx(-1.0)

    def test_matches_independent_reimplementation(self):
        rng = np.random.default_rng(17)
        for _ in range(25):
            n = int(rng.integers(2, 40))
            x = rng.normal(size=n)
            y = 0.8 * x + rng.normal(scale=0.5, size=n) + rng.normal()
            assert concordance_correlation_coefficient(x, y) == pytest.approx(ref_ccc(x, y))

    def test_ddof_is_internally_consistent(self):
        """A ddof mismatch between covariance and variances is the classic CCC
        bug. Build the ddof=0 and ddof=1 variants explicitly; for n=5 they
        differ, and the implementation must equal the ddof=0 (Lin) form."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 7.0])
        y = np.array([2.0, 2.0, 4.0, 3.0, 9.0])
        n = len(x)
        mx, my = x.mean(), y.mean()
        bias = (mx - my) ** 2
        cov0 = ((x - mx) * (y - my)).sum() / n
        cov1 = ((x - mx) * (y - my)).sum() / (n - 1)
        consistent0 = 2 * cov0 / (x.var(ddof=0) + y.var(ddof=0) + bias)
        consistent1 = 2 * cov1 / (x.var(ddof=1) + y.var(ddof=1) + bias)
        mismatched = 2 * cov0 / (x.var(ddof=1) + y.var(ddof=1) + bias)
        assert consistent0 != pytest.approx(mismatched)
        got = concordance_correlation_coefficient(x, y)
        assert got == pytest.approx(consistent0, rel=1e-12)
        assert got != pytest.approx(mismatched, rel=1e-9)
        # Lin's CCC is scale-invariant in its normalisation only when both
        # moments use the same ddof; the two consistent forms differ because
        # the bias term is unscaled, so pin the population one.
        assert consistent0 != pytest.approx(consistent1)

    def test_bounded_in_minus_one_to_one(self):
        rng = np.random.default_rng(99)
        for _ in range(200):
            n = int(rng.integers(2, 25))
            x = rng.normal(scale=rng.uniform(0.1, 10.0), size=n) + rng.uniform(-20, 20)
            y = rng.normal(scale=rng.uniform(0.1, 10.0), size=n) + rng.uniform(-20, 20)
            ccc = concordance_correlation_coefficient(x, y)
            assert -1.0 - 1e-12 <= ccc <= 1.0 + 1e-12

    def test_never_exceeds_pearson_r_in_magnitude(self):
        """CCC = r * bias-correction factor, and the factor is in (0, 1], so
        |CCC| <= |r|. This is what makes CCC catch systematic bias."""
        rng = np.random.default_rng(5)
        for _ in range(40):
            n = int(rng.integers(3, 30))
            x = rng.normal(size=n)
            y = 1.5 * x + 3.0 + rng.normal(scale=0.3, size=n)
            r = np.corrcoef(x, y)[0, 1]
            ccc = concordance_correlation_coefficient(x, y)
            assert abs(ccc) <= abs(r) + 1e-12

    def test_penalises_additive_bias_where_pearson_does_not(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = x + 10.0
        assert np.corrcoef(x, y)[0, 1] == pytest.approx(1.0)
        assert concordance_correlation_coefficient(x, y) < 0.2

    def test_penalises_scale_bias_where_pearson_does_not(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = 3.0 * x
        assert np.corrcoef(x, y)[0, 1] == pytest.approx(1.0)
        assert concordance_correlation_coefficient(x, y) < 1.0

    def test_symmetric_in_its_arguments(self):
        x = [1.0, 4.0, 2.0, 8.0]
        y = [2.0, 3.0, 2.5, 7.0]
        assert concordance_correlation_coefficient(x, y) == pytest.approx(
            concordance_correlation_coefficient(y, x)
        )

    def test_mismatched_lengths_return_nan(self):
        assert math.isnan(concordance_correlation_coefficient([1, 2, 3], [1, 2]))

    def test_fewer_than_two_points_returns_nan(self):
        assert math.isnan(concordance_correlation_coefficient([1.0], [1.0]))
        assert math.isnan(concordance_correlation_coefficient([], []))

    def test_both_constant_and_equal_returns_zero(self):
        assert concordance_correlation_coefficient([2.0] * 5, [2.0] * 5) == 0.0

    def test_all_zero_inputs_return_zero(self):
        assert concordance_correlation_coefficient([0.0] * 4, [0.0] * 4) == 0.0

    def test_constant_but_unequal_inputs_are_not_the_zero_branch(self):
        """Constant inputs with different means have a non-zero denominator
        (the bias term), so CCC is 0 by the formula, not by the guard."""
        got = concordance_correlation_coefficient([1.0] * 4, [5.0] * 4)
        assert got == pytest.approx(0.0)

    def test_nan_input_propagates(self):
        assert math.isnan(
            concordance_correlation_coefficient([1.0, 2.0, float("nan")], [1.0, 2.0, 3.0])
        )

    def test_accepts_integer_lists(self):
        assert concordance_correlation_coefficient([1, 2, 3], [1, 2, 3]) == pytest.approx(1.0)


# ==========================================================================
# benjamini_hochberg
# ==========================================================================


class TestBenjaminiHochberg:
    def test_bh_does_not_reject_above_largest_passing_rank(self):
        """REGRESSION: p=[0.01, 0.9], q=0.05, m=2.
        Thresholds are 0.025 and 0.05. Only rank 1 passes (0.01 <= 0.025), so
        k_max = 1 and ONLY the 0.01 hypothesis is rejected.
        The old code propagated significance to all HIGHER ranks and returned
        [True, True], declaring p=0.9 significant.
        """
        assert benjamini_hochberg([0.01, 0.9], q=0.05) == [True, False]

    def test_bh_rejects_everything_below_largest_passing_rank(self):
        """REGRESSION: p=[0.04, 0.045], q=0.05, m=2.
        Rank 1: 0.04 > 0.025 (fails on its own). Rank 2: 0.045 <= 0.05 (passes).
        Largest passing rank k_max = 2, so BOTH are rejected -- the step-up
        property. The old code returned [False, True].
        """
        assert benjamini_hochberg([0.04, 0.045], q=0.05) == [True, True]

    def test_step_up_distinguishes_from_per_hypothesis_thresholding(self):
        """p = [0.001, 0.03, 0.031], q=0.05, m=3; thresholds .0167/.0333/.05.
        Rank 2 (0.03 <= 0.0333) and rank 3 (0.031 <= 0.05) pass; rank 1 passes
        too. k_max = 3 so all three reject. Naive per-hypothesis thresholding
        would also reject all here, so pair it with a case where rank 1 fails.
        """
        assert benjamini_hochberg([0.001, 0.03, 0.031], q=0.05) == [True, True, True]
        # Rank 1 fails its own threshold but is still rejected via k_max=3.
        p = [0.02, 0.03, 0.031]
        assert p[0] > 0.05 * 1 / 3  # 0.02 > 0.01667, fails individually
        assert benjamini_hochberg(p, q=0.05) == [True, True, True]

    def test_booleans_map_back_to_original_input_order(self):
        """Shuffled input: the significant flag must follow the p-value, not
        the sorted position."""
        p = [0.9, 0.001, 0.8, 0.002, 0.7]
        got = benjamini_hochberg(p, q=0.05)
        assert got == [False, True, False, True, False]
        # And under a different permutation of the same multiset:
        order = [3, 0, 4, 1, 2]
        permuted = [p[i] for i in order]
        got_perm = benjamini_hochberg(permuted, q=0.05)
        assert got_perm == [got[i] for i in order]

    @pytest.mark.parametrize(
        "p_values",
        [
            [0.01, 0.9],
            [0.04, 0.045],
            [0.001, 0.2, 0.3, 0.4, 0.5],
            [0.9, 0.01, 0.8],
            [0.01, 0.02, 0.03, 0.9, 0.95],
            [0.005, 0.011, 0.02, 0.04, 0.13],
            [0.5],
            [0.049, 0.049, 0.049],
            [0.0, 0.0, 0.5],
            [1.0, 1.0, 1.0],
        ],
    )
    @pytest.mark.parametrize("q", [0.01, 0.05, 0.1, 0.5])
    def test_matches_scipy_false_discovery_control(self, p_values, q):
        assert benjamini_hochberg(p_values, q=q) == ref_bh(p_values, q)

    def test_matches_scipy_on_random_inputs(self):
        rng = np.random.default_rng(2024)
        for _ in range(150):
            m = int(rng.integers(1, 30))
            p = rng.uniform(0.0, 1.0, size=m)
            # Mix in a few strongly significant values so rejections occur.
            if rng.random() < 0.5:
                p[: max(1, m // 4)] = rng.uniform(0.0, 0.01, size=max(1, m // 4))
            q = float(rng.choice([0.01, 0.05, 0.1, 0.2]))
            assert benjamini_hochberg(list(p), q=q) == ref_bh(list(p), q)

    def test_empty_input_returns_empty_list(self):
        assert benjamini_hochberg([], q=0.05) == []

    def test_single_p_value_reduces_to_plain_comparison(self):
        assert benjamini_hochberg([0.04], q=0.05) == [True]
        assert benjamini_hochberg([0.06], q=0.05) == [False]
        assert benjamini_hochberg([0.05], q=0.05) == [True]  # boundary is inclusive

    def test_all_zero_p_values_all_reject(self):
        assert benjamini_hochberg([0.0] * 6, q=0.05) == [True] * 6

    def test_all_ones_p_values_reject_nothing(self):
        assert benjamini_hochberg([1.0] * 6, q=0.05) == [False] * 6

    def test_q_equals_one_rejects_everything(self):
        """With q=1 the largest threshold is 1.0, which every p-value in [0,1]
        meets at rank m, so k_max = m and all hypotheses reject."""
        rng = np.random.default_rng(8)
        p = list(rng.uniform(0.0, 1.0, size=12))
        assert benjamini_hochberg(p, q=1.0) == [True] * 12

    def test_q_equals_zero_rejects_nothing_unless_p_is_zero(self):
        assert benjamini_hochberg([0.1, 0.5, 0.9], q=0.0) == [False] * 3
        assert benjamini_hochberg([0.0, 0.5], q=0.0) == [True, False]

    def test_default_q_is_five_percent(self):
        import inspect

        assert inspect.signature(benjamini_hochberg).parameters["q"].default == 0.05

    def test_ties_are_handled_consistently(self):
        """Equal p-values must get the same verdict regardless of sort order."""
        got = benjamini_hochberg([0.02, 0.02, 0.02, 0.9], q=0.05)
        assert got[:3] == [True, True, True]
        assert got[3] is False

    def test_rejection_set_is_monotone_in_q(self):
        rng = np.random.default_rng(31)
        p = list(rng.uniform(0.0, 0.3, size=15))
        prev = set()
        for q in [0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0]:
            rejected = {i for i, s in enumerate(benjamini_hochberg(p, q=q)) if s}
            assert prev <= rejected
            prev = rejected

    def test_is_at_least_as_powerful_as_bonferroni(self):
        rng = np.random.default_rng(77)
        for _ in range(30):
            m = int(rng.integers(2, 20))
            p = list(rng.uniform(0.0, 0.1, size=m))
            q = 0.05
            bh = benjamini_hochberg(p, q=q)
            bonf = [pv <= q / m for pv in p]
            for i in range(m):
                if bonf[i]:
                    assert bh[i], "BH must reject whatever Bonferroni rejects"

    def test_returns_plain_python_bools(self):
        got = benjamini_hochberg([0.01, 0.9], q=0.05)
        assert all(isinstance(v, bool) for v in got)


# ==========================================================================
# rmst_difference
# ==========================================================================


class TestRMSTDifference:
    def test_rmst_single_event_uses_survival_before_the_drop(self):
        """REGRESSION: one subject with an event at t=5, horizon 10.
        S(t) = 1 for t < 5 and 0 for t >= 5, so RMST = 5.0 exactly.
        Compare against a reference arm of one subject censored at t=20, whose
        RMST is the full horizon 10.0, giving a difference of 10 - 5 = 5.0.
        The old code multiplied the interval (0, 5] by the POST-drop survival
        (0.0) and returned a difference of 10.0, i.e. an RMST of 0 for the
        event arm -- the entire area under the curve was lost.
        """
        event_arm_t, event_arm_e = np.array([5.0]), np.array([1.0])
        survivor_t, survivor_e = np.array([20.0]), np.array([0.0])
        diff = rmst_difference(event_arm_t, event_arm_e, survivor_t, survivor_e, 10.0)
        assert diff == pytest.approx(5.0)

    def test_km_rmst_matches_hand_worked_example(self):
        """REGRESSION / hand-worked Kaplan-Meier, n=6, tau=6.
        times  1  2  3  4  5  6
        event  1  0  1  0  1  0
          t=1: at risk 6, d=1 -> S = 5/6   = 0.8333333
          t=2: censored, no drop
          t=3: at risk 4, d=1 -> S = 0.8333333 * 3/4 = 0.625
          t=4: censored, no drop
          t=5: at risk 2, d=1 -> S = 0.625 * 1/2      = 0.3125
          t=6: censored, no drop
        RMST(6) = 1.0*(1-0) + 0.8333333*(3-1) + 0.625*(5-3) + 0.3125*(6-5)
                = 1.0 + 1.6666667 + 1.25 + 0.3125 = 4.2291666667
        """
        t = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        e = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
        # Reference arm: fully censored beyond the horizon, so RMST == tau == 6.
        ref_t, ref_e = np.array([99.0]), np.array([0.0])
        expected = 1.0 + (5.0 / 6.0) * 2 + 0.625 * 2 + 0.3125
        assert expected == pytest.approx(4.229166666666667)
        diff = rmst_difference(t, e, ref_t, ref_e, 6.0)
        assert 6.0 - diff == pytest.approx(expected, rel=1e-12)
        assert 6.0 - diff == pytest.approx(ref_km_rmst(t, e, 6.0), rel=1e-12)

    def test_censoring_changes_the_estimate(self):
        """A censored observation must NOT be treated as an event. Same times,
        different event flags -> different RMST. (A naive 'everyone who left
        the study died' estimator would give the same answer for both.)"""
        t = np.array([2.0, 4.0, 6.0, 8.0])
        all_events = np.array([1.0, 1.0, 1.0, 1.0])
        some_censored = np.array([1.0, 0.0, 1.0, 0.0])
        ref_t, ref_e = np.array([999.0]), np.array([0.0])
        tau = 10.0
        rmst_all = tau - rmst_difference(t, all_events, ref_t, ref_e, tau)
        rmst_cens = tau - rmst_difference(t, some_censored, ref_t, ref_e, tau)
        assert rmst_all != pytest.approx(rmst_cens)
        # Censoring can only raise the KM curve relative to calling them events.
        assert rmst_cens > rmst_all

    def test_censored_subject_stays_in_the_risk_set_at_its_own_time(self):
        """Tie between an event and a censoring at t=3: the censored subject is
        at risk at t=3, so n_at_risk=3 and S drops to 2/3, not to 1/2."""
        t = np.array([3.0, 3.0, 5.0])
        e = np.array([1.0, 0.0, 0.0])  # event, censored (tied), censored later
        ref_t, ref_e = np.array([999.0]), np.array([0.0])
        tau = 6.0
        rmst = tau - rmst_difference(t, e, ref_t, ref_e, tau)
        expected = 1.0 * 3.0 + (2.0 / 3.0) * (6.0 - 3.0)  # = 3 + 2 = 5.0
        assert rmst == pytest.approx(expected)
        assert rmst == pytest.approx(ref_km_rmst(t, e, tau))

    def test_matches_reference_km_on_random_censored_data(self):
        rng = np.random.default_rng(1234)
        for _ in range(60):
            n = int(rng.integers(1, 30))
            t = np.round(rng.uniform(0.5, 20.0, size=n), 1)
            e = (rng.random(n) < 0.6).astype(float)
            tau = float(rng.choice([5.0, 12.0, 20.0, 30.0]))
            other_t, other_e = np.array([1e6]), np.array([0.0])
            rmst = tau - rmst_difference(t, e, other_t, other_e, tau)
            assert rmst == pytest.approx(ref_km_rmst(t, e, tau), rel=1e-12, abs=1e-12)

    def test_integration_is_a_step_function_not_a_trapezoid(self):
        """Single event at t=5 with horizon 10: the step-function area is 5.0;
        the trapezoid rule over (0, 5] would give 1.0 * 5 / 2 = 2.5 plus a
        0-height tail, i.e. 2.5. Pin the step answer."""
        t, e = np.array([5.0]), np.array([1.0])
        ref_t, ref_e = np.array([1e6]), np.array([0.0])
        rmst = 10.0 - rmst_difference(t, e, ref_t, ref_e, 10.0)
        assert rmst == pytest.approx(5.0)
        assert rmst != pytest.approx(2.5)

    def test_identical_arms_give_zero_difference(self):
        t = np.array([1.0, 3.0, 5.0, 7.0])
        e = np.array([1.0, 0.0, 1.0, 1.0])
        assert rmst_difference(t, e, t, e, 6.0) == pytest.approx(0.0)

    def test_sign_convention_is_b_minus_a(self):
        """Arm B survives longer, so the returned difference must be positive."""
        a_t, a_e = np.array([1.0, 2.0, 3.0]), np.array([1.0, 1.0, 1.0])
        b_t, b_e = np.array([8.0, 9.0, 10.0]), np.array([1.0, 1.0, 1.0])
        tau = 10.0
        assert rmst_difference(a_t, a_e, b_t, b_e, tau) > 0
        assert rmst_difference(b_t, b_e, a_t, a_e, tau) == pytest.approx(
            -rmst_difference(a_t, a_e, b_t, b_e, tau)
        )

    def test_antisymmetry_on_random_arms(self):
        rng = np.random.default_rng(55)
        for _ in range(20):
            a_t = rng.uniform(1, 15, size=10)
            a_e = (rng.random(10) < 0.7).astype(float)
            b_t = rng.uniform(1, 15, size=12)
            b_e = (rng.random(12) < 0.7).astype(float)
            ab = rmst_difference(a_t, a_e, b_t, b_e, 10.0)
            ba = rmst_difference(b_t, b_e, a_t, a_e, 10.0)
            assert ab == pytest.approx(-ba)

    def test_full_censoring_gives_rmst_equal_to_horizon(self):
        """100% censoring: KM never drops, so RMST == tau for both arms."""
        t = np.array([1.0, 2.0, 3.0])
        e = np.zeros(3)
        other_t, other_e = np.array([10.0, 20.0]), np.array([0.0, 0.0])
        assert rmst_difference(t, e, other_t, other_e, 5.0) == pytest.approx(0.0)
        # And against an arm that does have events, the all-censored arm is
        # the maximum possible RMST.
        ev_t, ev_e = np.array([1.0, 2.0, 3.0]), np.array([1.0, 1.0, 1.0])
        assert rmst_difference(ev_t, ev_e, t, e, 5.0) > 0

    def test_no_events_before_horizon_gives_full_horizon(self):
        t = np.array([50.0, 60.0])
        e = np.array([1.0, 1.0])
        other_t, other_e = np.array([70.0]), np.array([1.0])
        # Both arms have all events after tau=10, so both RMSTs are 10.
        assert rmst_difference(t, e, other_t, other_e, 10.0) == pytest.approx(0.0)

    def test_events_after_horizon_are_truncated_away(self):
        """An event at t=20 must not affect RMST at tau=10."""
        ref_t, ref_e = np.array([1e6]), np.array([0.0])
        a = 10.0 - rmst_difference(np.array([20.0]), np.array([1.0]), ref_t, ref_e, 10.0)
        b = 10.0 - rmst_difference(np.array([20.0]), np.array([0.0]), ref_t, ref_e, 10.0)
        assert a == pytest.approx(10.0)
        assert b == pytest.approx(10.0)

    def test_tau_beyond_longest_followup_carries_last_value_forward(self):
        """One event at t=2 (S drops to 0.5) and one censored at t=4, tau=10.
        KM is undefined past t=4; the convention is to carry S=0.5 forward.
        RMST = 1.0*2 + 0.5*(10-2) = 2 + 4 = 6.0
        """
        t = np.array([2.0, 4.0])
        e = np.array([1.0, 0.0])
        ref_t, ref_e = np.array([1e6]), np.array([0.0])
        rmst = 10.0 - rmst_difference(t, e, ref_t, ref_e, 10.0)
        assert rmst == pytest.approx(6.0)
        assert rmst == pytest.approx(ref_km_rmst(t, e, 10.0))

    def test_rmst_is_bounded_by_zero_and_horizon(self):
        rng = np.random.default_rng(606)
        ref_t, ref_e = np.array([1e9]), np.array([0.0])
        for _ in range(50):
            n = int(rng.integers(1, 25))
            t = rng.uniform(0.1, 30.0, size=n)
            e = (rng.random(n) < 0.5).astype(float)
            tau = 15.0
            rmst = tau - rmst_difference(t, e, ref_t, ref_e, tau)
            assert 0.0 <= rmst <= tau + 1e-12

    def test_rmst_is_monotone_nondecreasing_in_horizon(self):
        t = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
        e = np.array([1.0, 0.0, 1.0, 0.0, 1.0])
        ref_t, ref_e = np.array([1e9]), np.array([0.0])
        prev = -1.0
        for tau in [1.0, 2.0, 5.0, 10.0, 20.0]:
            rmst = tau - rmst_difference(t, e, ref_t, ref_e, tau)
            assert rmst >= prev - 1e-12
            prev = rmst

    def test_empty_arm_gives_nan(self):
        t, e = np.array([1.0, 2.0]), np.array([1.0, 1.0])
        empty = np.array([])
        assert math.isnan(rmst_difference(empty, empty, t, e, 5.0))
        assert math.isnan(rmst_difference(t, e, empty, empty, 5.0))

    def test_accepts_lists_and_integer_arrays(self):
        from_lists = rmst_difference([1, 3, 5], [1, 0, 1], [2, 4, 6], [1, 1, 0], 6)
        from_arrays = rmst_difference(
            np.array([1.0, 3.0, 5.0]),
            np.array([1.0, 0.0, 1.0]),
            np.array([2.0, 4.0, 6.0]),
            np.array([1.0, 1.0, 0.0]),
            6.0,
        )
        assert from_lists == pytest.approx(from_arrays)

    def test_unsorted_input_gives_the_same_answer_as_sorted(self):
        t = np.array([7.0, 1.0, 5.0, 3.0])
        e = np.array([1.0, 1.0, 0.0, 1.0])
        order = np.argsort(t)
        ref_t, ref_e = np.array([1e9]), np.array([0.0])
        a = rmst_difference(t, e, ref_t, ref_e, 8.0)
        b = rmst_difference(t[order], e[order], ref_t, ref_e, 8.0)
        assert a == pytest.approx(b)

    def test_all_events_at_the_same_time(self):
        """5 subjects, all with an event at t=4, tau=10: S drops 1 -> 0 at 4,
        so RMST = 4.0."""
        t = np.full(5, 4.0)
        e = np.ones(5)
        ref_t, ref_e = np.array([1e9]), np.array([0.0])
        rmst = 10.0 - rmst_difference(t, e, ref_t, ref_e, 10.0)
        assert rmst == pytest.approx(4.0)

    def test_event_exactly_at_the_horizon_is_included(self):
        """An event at t == tau contributes no area after it, so RMST == tau
        either way; but the drop must be applied, not skipped, so a later
        comparison stays consistent with the reference."""
        t, e = np.array([5.0, 10.0]), np.array([1.0, 1.0])
        ref_t, ref_e = np.array([1e9]), np.array([0.0])
        rmst = 10.0 - rmst_difference(t, e, ref_t, ref_e, 10.0)
        assert rmst == pytest.approx(ref_km_rmst(t, e, 10.0))
        assert rmst == pytest.approx(5.0 + 0.5 * 5.0)  # 1*5 + 0.5*(10-5) = 7.5


# ==========================================================================
# Cross-function integration: the shape of a T1-T5 PASS decision.
# ==========================================================================


def test_preregistered_pass_decision_uses_ci_lower_bound():
    """A metric whose point estimate clears the threshold but whose bootstrap
    lower bound does not must NOT pass."""
    rng = np.random.default_rng(314)
    data = rng.normal(0.72, 0.20, size=30)
    threshold = 0.70
    point, lo, _hi = bootstrap_ci(data, n_resamples=1000, seed=314)
    assert point > threshold
    assert lo < threshold
    assert not (lo >= threshold), "PASS must be decided on the lower bound"


def test_rmst_difference_with_bootstrap_ci_is_deterministic():
    """The headline endpoint: a bootstrap CI around the RMST difference, which
    must be byte-identical across runs with the same seed."""
    rng = np.random.default_rng(1)
    a_t = rng.uniform(1, 24, size=40)
    a_e = (rng.random(40) < 0.7).astype(float)
    b_t = rng.uniform(4, 30, size=40)
    b_e = (rng.random(40) < 0.7).astype(float)

    def stat(idx):
        idx = np.asarray(idx, dtype=int)
        return rmst_difference(a_t[idx], a_e[idx], b_t[idx], b_e[idx], 18.0)

    first = bootstrap_ci(np.arange(40), statistic=stat, n_resamples=100, seed=2)
    second = bootstrap_ci(np.arange(40), statistic=stat, n_resamples=100, seed=2)
    assert first == second
    assert first[1] <= first[0] <= first[2] or math.isnan(first[0])
