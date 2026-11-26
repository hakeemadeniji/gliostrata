"""Tests for `gliostrata.prediction.cox` -- the real, fitted Cox PH model.

The central test is ground-truth recovery: survival times are simulated from a
Cox model with KNOWN coefficients and the fitter must recover them. That is the
only way to prove the partial-likelihood maximiser is actually correct, as
opposed to merely producing plausible-looking numbers.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gliostrata.prediction.cox import (
    CoxPHModel,
    brier_score_ipcw,
    build_design_matrix,
    concordance_index,
    cross_validate_cox,
    fit_clinical_baseline,
)

REAL_COHORT = Path(__file__).resolve().parents[1] / "data/raw/tcga_gbm/clinical_cohort.parquet"

TRUE_BETA = np.array([0.80, -0.50, 0.30])
BASELINE_HAZARD = 0.1


def simulate_cox(
    n: int = 4000,
    beta: np.ndarray = TRUE_BETA,
    *,
    seed: int = 0,
    censor_scale: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Exponential-baseline Cox simulation with a KNOWN beta.

    With h(t|x) = h0 * exp(x'beta) and constant h0, the survival time is
    T = -log(U) / (h0 * exp(x'beta)) -- an exact inverse-CDF draw.

    Returns (X, observed_time, event, true_event_time). Censoring is applied
    independently (administrative/random right censoring), so `event == 0`
    means "still alive when follow-up ended".
    """
    rng = np.random.default_rng(seed)
    beta = np.asarray(beta, dtype=float)
    X = rng.normal(size=(n, len(beta)))
    true_time = -np.log(rng.uniform(size=n)) / (BASELINE_HAZARD * np.exp(X @ beta))
    if censor_scale is None:
        return X, true_time, np.ones(n, dtype=int), true_time
    censor_time = rng.exponential(scale=censor_scale, size=n)
    observed = np.minimum(true_time, censor_time)
    event = (true_time <= censor_time).astype(int)
    return X, observed, event, true_time


# ---------------------------------------------------------------------------
# 1. Ground-truth coefficient recovery
# ---------------------------------------------------------------------------
def test_recovers_known_coefficients_without_censoring():
    X, time, event, _ = simulate_cox(n=6000, seed=1)
    model = CoxPHModel().fit(X, time, event)
    assert model.converged_
    assert np.allclose(model.coef_, TRUE_BETA, atol=0.06), model.coef_
    # Wald CIs from the analytic Hessian should cover the truth.
    summary = model.summary()
    assert np.all(summary["hr_lower"] < np.exp(TRUE_BETA))
    assert np.all(summary["hr_upper"] > np.exp(TRUE_BETA))


def test_recovers_known_coefficients_with_heavy_censoring():
    X, time, event, _ = simulate_cox(n=6000, seed=1, censor_scale=10.0)
    censored_fraction = 1.0 - event.mean()
    assert 0.3 < censored_fraction < 0.7, censored_fraction
    model = CoxPHModel().fit(X, time, event)
    assert model.converged_
    # Censoring costs precision but must not cost correctness.
    assert np.allclose(model.coef_, TRUE_BETA, atol=0.10), model.coef_


def test_recovery_improves_with_sample_size():
    errors = []
    for n in (500, 8000):
        X, time, event, _ = simulate_cox(n=n, seed=7, censor_scale=10.0)
        model = CoxPHModel().fit(X, time, event)
        errors.append(np.abs(model.coef_ - TRUE_BETA).max())
    assert errors[1] < errors[0]


def test_standardization_does_not_change_the_fit():
    X, time, event, _ = simulate_cox(n=2000, seed=3, censor_scale=10.0)
    a = CoxPHModel(standardize=True).fit(X, time, event)
    b = CoxPHModel(standardize=False).fit(X, time, event)
    assert np.allclose(a.coef_, b.coef_, atol=1e-5)


def test_feature_rescaling_rescales_coefficient_inversely():
    """Cox coefficients are per-unit: scaling a column by 10 divides its coef by 10."""
    X, time, event, _ = simulate_cox(n=2000, seed=4, censor_scale=10.0)
    base = CoxPHModel().fit(X, time, event)
    scaled = X.copy()
    scaled[:, 0] *= 10.0
    other = CoxPHModel().fit(scaled, time, event)
    assert other.coef_[0] == pytest.approx(base.coef_[0] / 10.0, rel=1e-4)


# ---------------------------------------------------------------------------
# 2. Right-censoring is handled correctly (and mishandling it biases results)
# ---------------------------------------------------------------------------
def test_censoring_handled_correctly_vs_treating_censored_as_dead():
    """Same underlying data, 0% vs 50% censoring, plus the classic error.

    Correct handling: coefficients stay near the truth under censoring.
    The common error -- treating `vital_status == 0` as a death at the
    censoring time -- attenuates every coefficient toward zero, because
    long-lived (low-risk) patients get assigned artificially short "deaths".
    """
    X, obs_time, event, true_time = simulate_cox(n=6000, seed=11, censor_scale=10.0)
    censored_fraction = 1.0 - event.mean()
    assert censored_fraction > 0.3

    uncensored = CoxPHModel().fit(X, true_time, np.ones(len(X), dtype=int))
    correct = CoxPHModel().fit(X, obs_time, event)
    naive = CoxPHModel().fit(X, obs_time, np.ones(len(X), dtype=int))

    # Correct handling tracks the fully observed fit.
    assert np.allclose(correct.coef_, uncensored.coef_, atol=0.08), (
        correct.coef_,
        uncensored.coef_,
    )

    # The naive fit is measurably and systematically attenuated.
    shrinkage = np.abs(naive.coef_) / np.abs(uncensored.coef_)
    assert np.all(shrinkage < 0.8), shrinkage
    naive_err = np.abs(naive.coef_ - TRUE_BETA).max()
    correct_err = np.abs(correct.coef_ - TRUE_BETA).max()
    assert naive_err > 3 * correct_err, (naive_err, correct_err)
    # Record the magnitude so a regression in censoring handling is obvious.
    assert naive_err > 0.2, naive_err


def test_censored_only_rows_still_enter_the_risk_sets():
    """Dropping censored patients is NOT equivalent to censoring them."""
    X, time, event, _ = simulate_cox(n=3000, seed=12, censor_scale=6.0)
    with_censored = CoxPHModel().fit(X, time, event)
    keep = event == 1
    deaths_only = CoxPHModel().fit(X[keep], time[keep], event[keep])
    # Deleting censored observations throws away information and biases the fit.
    assert not np.allclose(with_censored.coef_, deaths_only.coef_, atol=0.02)
    assert np.abs(with_censored.coef_ - TRUE_BETA).max() < np.abs(
        deaths_only.coef_ - TRUE_BETA
    ).max()


def test_event_column_validated():
    X, time, event, _ = simulate_cox(n=100, seed=1)
    with pytest.raises(ValueError, match="event must be 0"):
        CoxPHModel().fit(X, time, np.full(100, 2))
    with pytest.raises(ValueError, match="zero events"):
        CoxPHModel().fit(X, time, np.zeros(100, dtype=int))


# ---------------------------------------------------------------------------
# 3. C-index properties
# ---------------------------------------------------------------------------
def test_c_index_of_perfect_predictor_is_one():
    rng = np.random.default_rng(0)
    time = rng.uniform(1, 100, size=300)
    event = np.ones(300, dtype=int)
    # risk perfectly inversely ordered with survival time
    assert concordance_index(time, event, -time) == pytest.approx(1.0)
    assert concordance_index(time, event, time) == pytest.approx(0.0)


def test_c_index_of_random_predictor_is_about_half():
    rng = np.random.default_rng(1)
    time = rng.exponential(20, size=2000)
    event = rng.binomial(1, 0.7, size=2000)
    risk = rng.normal(size=2000)
    assert concordance_index(time, event, risk) == pytest.approx(0.5, abs=0.03)


def test_c_index_ties_score_half():
    time = np.array([1.0, 2.0, 3.0])
    event = np.array([1, 1, 1])
    assert concordance_index(time, event, np.zeros(3)) == pytest.approx(0.5)


def test_c_index_only_uses_comparable_pairs():
    """A patient censored early is not scored as if they had died early."""
    time = np.array([5.0, 10.0])
    # patient 0 censored at 5, patient 1 dies at 10 -> only pair (1,0) is
    # non-comparable from 1's side and (0,1) is non-comparable because 0 is
    # censored, so the single comparable pair must come from the death.
    with pytest.raises(ValueError, match="no comparable pairs"):
        concordance_index(time, np.array([0, 0]), np.array([1.0, 0.0]))
    # With the death first, the pair is comparable.
    assert concordance_index(
        np.array([5.0, 10.0]), np.array([1, 0]), np.array([1.0, 0.0])
    ) == pytest.approx(1.0)


def test_fitted_model_beats_chance_on_simulated_data():
    X, time, event, _ = simulate_cox(n=2000, seed=5, censor_scale=10.0)
    model = CoxPHModel().fit(X, time, event)
    c = concordance_index(time, event, model.predict_log_hazard(X))
    assert c > 0.7


# ---------------------------------------------------------------------------
# 4. Baseline survival function (Breslow)
# ---------------------------------------------------------------------------
def test_baseline_survival_is_a_valid_monotone_curve():
    X, time, event, _ = simulate_cox(n=1500, seed=6, censor_scale=10.0)
    model = CoxPHModel().fit(X, time, event)
    grid = np.linspace(0.1, np.quantile(time, 0.9), 40)
    S = model.predict_survival_function(X[:20], grid)
    assert S.shape == (20, 40)
    assert np.all((S >= 0) & (S <= 1))
    assert np.all(np.diff(S, axis=1) <= 1e-12)  # non-increasing in t
    assert model.baseline_cumulative_hazard(0.0) == pytest.approx(0.0)


def test_higher_risk_patients_have_lower_survival_curves():
    X, time, event, _ = simulate_cox(n=2000, seed=8, censor_scale=10.0)
    model = CoxPHModel().fit(X, time, event)
    risk = model.predict_log_hazard(X)
    lo, hi = np.argmin(risk), np.argmax(risk)
    S = model.predict_survival_function(X[[lo, hi]], [6.0, 12.0, 18.0])
    assert np.all(S[0] > S[1])


def test_breslow_baseline_recovers_known_exponential_baseline():
    """With h0 = 0.1 constant, H0(t) should be ~0.1 * t near the mean patient."""
    X, time, event, _ = simulate_cox(n=8000, seed=9, censor_scale=30.0)
    model = CoxPHModel().fit(X, time, event)
    # The baseline is anchored at the mean covariate vector; E[exp(x'b)] > 1
    # for mean-zero x, so the curve at x=mean is the population baseline
    # rescaled. Check linearity in t rather than the absolute level.
    t = np.array([2.0, 4.0, 6.0, 8.0])
    H = model.baseline_cumulative_hazard(t)
    ratios = H[1:] / H[:-1]
    expected = t[1:] / t[:-1]
    assert np.allclose(ratios, expected, rtol=0.2), (ratios, expected)


# ---------------------------------------------------------------------------
# 5. Brier score with IPCW
# ---------------------------------------------------------------------------
def test_brier_score_of_uninformative_half_prediction_is_quarter():
    rng = np.random.default_rng(2)
    time = rng.exponential(20, size=4000)
    event = np.ones(4000, dtype=int)
    bs = brier_score_ipcw(time, event, np.full(4000, 0.5), horizon=18.0)
    assert bs == pytest.approx(0.25, abs=0.01)


def test_brier_score_rewards_the_better_model():
    rng = np.random.default_rng(3)
    n = 3000
    time = rng.exponential(20, size=n)
    event = np.ones(n, dtype=int)
    truth = (time > 18.0).astype(float)
    good = np.clip(truth * 0.9 + 0.05, 0, 1)
    bad = 1.0 - good
    assert brier_score_ipcw(time, event, good, 18.0) < brier_score_ipcw(
        time, event, bad, 18.0
    )


def test_brier_score_ignores_patients_censored_before_the_horizon():
    """Censored-before-horizon patients get zero weight, not a guessed label."""
    time = np.array([5.0, 30.0, 30.0, 30.0])
    event = np.array([0, 1, 1, 1])
    s_bad_for_censored = np.array([0.0, 0.9, 0.9, 0.9])
    s_good_for_censored = np.array([1.0, 0.9, 0.9, 0.9])
    a = brier_score_ipcw(time, event, s_bad_for_censored, 18.0)
    b = brier_score_ipcw(time, event, s_good_for_censored, 18.0)
    assert a == pytest.approx(b)


def test_brier_score_ipcw_reweights_under_censoring():
    """IPCW keeps the estimate near the uncensored truth despite censoring."""
    rng = np.random.default_rng(4)
    n = 6000
    true_time = rng.exponential(20, size=n)
    uncensored = brier_score_ipcw(
        true_time, np.ones(n, dtype=int), np.full(n, 0.5), 18.0
    )
    c = rng.exponential(25, size=n)
    obs = np.minimum(true_time, c)
    ev = (true_time <= c).astype(int)
    censored = brier_score_ipcw(obs, ev, np.full(n, 0.5), 18.0)
    assert censored == pytest.approx(uncensored, abs=0.02)


def test_brier_score_length_validated():
    with pytest.raises(ValueError, match="one entry per patient"):
        brier_score_ipcw(np.array([1.0, 2.0]), np.array([1, 1]), np.array([0.5]), 1.0)


# ---------------------------------------------------------------------------
# 6. Missingness strategy is applied exactly as documented
# ---------------------------------------------------------------------------
def _toy_cohort() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "age": [50.0, 60.0, 70.0, 80.0, 55.0, 65.0],
            "kps": [90.0, np.nan, 70.0, 60.0, np.nan, 80.0],
            "sex": ["M", "F", "M", "M", "F", "M"],
            "treatment_arm": ["ChemoRadiation"] * 4 + [None, "MonoTherapy"],
            "overall_survival_months": [10.0, 20.0, 5.0, 30.0, 15.0, 25.0],
            "vital_status": [1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
        }
    )


def test_median_indicator_imputes_and_adds_indicator():
    df = _toy_cohort()
    X, time, event, report = build_design_matrix(
        df, ["age", "kps"], impute="median_indicator"
    )
    assert report.feature_names == ["age", "kps_missing", "kps"]
    assert report.n_imputed["kps"] == 2
    expected_median = float(np.nanmedian(df["kps"].to_numpy(dtype=float)))
    assert report.imputed_values["kps"] == pytest.approx(expected_median)
    kps_missing = X[:, report.feature_names.index("kps_missing")]
    kps = X[:, report.feature_names.index("kps")]
    assert kps_missing.tolist() == [0, 1, 0, 0, 1, 0]
    assert kps[1] == pytest.approx(expected_median)
    assert kps[4] == pytest.approx(expected_median)
    assert not np.isnan(X).any()


def test_median_only_omits_the_indicator():
    df = _toy_cohort()
    X, _, _, report = build_design_matrix(df, ["age", "kps"], impute="median_only")
    assert report.feature_names == ["age", "kps"]
    assert "kps_missing" not in report.feature_names
    assert X.shape[1] == 2
    assert report.n_imputed["kps"] == 2


def test_impute_drop_refuses_rather_than_silently_dropping():
    df = _toy_cohort()
    with pytest.raises(ValueError, match="requires pre-filtered data"):
        build_design_matrix(df, ["age", "kps"], impute="drop")


def test_categorical_missing_becomes_an_explicit_level():
    df = _toy_cohort()
    _, _, _, report = build_design_matrix(df, ["age", "treatment_arm"])
    assert "treatment_arm_missing" in report.feature_names
    assert "missing" in report.categorical_levels["treatment_arm"]
    # The most frequent level is the reference and gets no dummy.
    assert "treatment_arm_ChemoRadiation" not in report.feature_names


def test_fully_missing_feature_is_refused_not_imputed():
    df = _toy_cohort()
    df["mgmt_methylated"] = np.nan
    with pytest.raises(ValueError, match="entirely missing"):
        build_design_matrix(df, ["age", "mgmt_methylated"])


def test_missing_outcome_rows_are_dropped_and_counted():
    df = _toy_cohort()
    df.loc[0, "vital_status"] = np.nan
    df.loc[1, "overall_survival_months"] = np.nan
    X, time, event, report = build_design_matrix(df, ["age", "kps"])
    assert report.n_input == 6
    assert report.n_used == 4
    assert report.dropped_missing_outcome == 2
    assert len(time) == len(event) == X.shape[0] == 4


def test_design_report_summary_mentions_imputation():
    df = _toy_cohort()
    _, _, _, report = build_design_matrix(df, ["age", "kps", "treatment_arm"])
    text = report.summary()
    assert "median_indicator" in text
    assert "imputed kps" in text
    assert "missing" in text


def test_censoring_flag_survives_design_matrix_construction():
    df = _toy_cohort()
    _, time, event, _ = build_design_matrix(df, ["age"])
    assert event.tolist() == [1, 0, 1, 1, 0, 1]
    assert time.tolist() == [10.0, 20.0, 5.0, 30.0, 15.0, 25.0]


# ---------------------------------------------------------------------------
# 7. Determinism
# ---------------------------------------------------------------------------
def test_fit_is_deterministic():
    X, time, event, _ = simulate_cox(n=1200, seed=21, censor_scale=10.0)
    a = CoxPHModel().fit(X, time, event)
    b = CoxPHModel().fit(X, time, event)
    assert np.array_equal(a.coef_, b.coef_)
    assert np.array_equal(a.se_, b.se_)
    assert np.array_equal(a.baseline_cumhaz_, b.baseline_cumhaz_)
    assert a.log_likelihood_ == b.log_likelihood_


def test_cross_validation_is_deterministic_under_a_fixed_seed():
    X, time, event, _ = simulate_cox(n=800, seed=22, censor_scale=10.0)
    kw = dict(n_splits=5, seed=1234, horizon=18.0, n_boot=120)
    a = cross_validate_cox(X, time, event, **kw)
    b = cross_validate_cox(X, time, event, **kw)
    assert a.c_index == b.c_index
    assert a.c_index_ci == b.c_index_ci
    assert a.brier_score == b.brier_score
    assert a.brier_ci == b.brier_ci
    assert np.array_equal(a.oof_risk, b.oof_risk)


def test_different_seed_changes_folds_but_not_wildly():
    X, time, event, _ = simulate_cox(n=800, seed=23, censor_scale=10.0)
    a = cross_validate_cox(X, time, event, seed=1, n_boot=50)
    b = cross_validate_cox(X, time, event, seed=2, n_boot=50)
    assert a.c_index != b.c_index
    assert abs(a.c_index - b.c_index) < 0.05


# ---------------------------------------------------------------------------
# 8. Cross-validation sanity
# ---------------------------------------------------------------------------
def test_cv_c_index_is_lower_than_in_sample_c_index():
    X, time, event, _ = simulate_cox(n=400, seed=24, censor_scale=10.0)
    cv = cross_validate_cox(X, time, event, n_boot=100)
    full = CoxPHModel().fit(X, time, event)
    in_sample = concordance_index(time, event, full.predict_log_hazard(X))
    assert cv.c_index <= in_sample + 1e-9
    assert np.isfinite(cv.oof_risk).all()
    assert len(cv.fold_c_indices) == 5


def test_cv_confidence_interval_brackets_the_point_estimate():
    X, time, event, _ = simulate_cox(n=600, seed=25, censor_scale=10.0)
    cv = cross_validate_cox(X, time, event, n_boot=200)
    lo, hi = cv.c_index_ci
    assert lo < cv.c_index < hi
    blo, bhi = cv.brier_ci
    assert blo < cv.brier_score < bhi
    assert "C-index" in cv.summary()


def test_cv_on_noise_features_is_near_chance():
    rng = np.random.default_rng(26)
    n = 1000
    X = rng.normal(size=(n, 3))
    time = rng.exponential(20, size=n)  # independent of X
    event = rng.binomial(1, 0.8, size=n)
    cv = cross_validate_cox(X, time, event, n_boot=50)
    assert cv.c_index == pytest.approx(0.5, abs=0.05)


# ---------------------------------------------------------------------------
# 9. Smoke test on the real TCGA-GBM cohort
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not REAL_COHORT.exists(), reason="real TCGA cohort not downloaded")
def test_smoke_real_tcga_gbm_cohort():
    cohort = pd.read_parquet(REAL_COHORT)
    # Canary on the ingested cohort size. If this trips, the L0 filter changed:
    # confirm the change was intended, refit, and update the number here rather
    # than loosening the assertion -- a silent cohort change invalidates every
    # downstream metric and every frozen partition.
    assert len(cohort) == 594

    model, cv, report = fit_clinical_baseline(cohort, n_boot=200, seed=42)

    # Missingness audit: only the 2 patients with no vital_status are dropped.
    assert report.dropped_missing_outcome == 2
    assert report.n_used == 592
    assert report.n_imputed["kps"] == 150  # of the 592 retained
    assert "kps_missing" in report.feature_names
    assert model.converged_
    assert model.n_obs_ == 592
    assert model.n_events_ == 491

    # Performance: better than chance, and honestly short of the >= 0.65
    # EXTERNAL threshold this internal CV cannot speak to.
    lo, hi = cv.c_index_ci
    assert 0.55 < cv.c_index < 0.80
    assert lo < cv.c_index < hi
    assert 0.0 < cv.brier_score < 0.25

    # Biology: in GBM, older age and lower KPS both mean worse survival.
    s = model.summary().set_index("feature")
    assert s.loc["age", "hazard_ratio"] > 1.0, "older age must increase hazard"
    assert s.loc["kps", "hazard_ratio"] < 1.0, "higher KPS must decrease hazard"
    assert s.loc["age", "p"] < 0.05
    assert s.loc["kps", "p"] < 0.05

    # The model emits a real survival curve, not just a score.
    S = model.predict_survival_function(
        np.asarray(model.center_)[None, :], [6.0, 12.0, 18.0, 24.0]
    )[0]
    assert np.all(np.diff(S) <= 0)
    assert 0.0 < S[-1] < S[0] < 1.0


@pytest.mark.skipif(not REAL_COHORT.exists(), reason="real TCGA cohort not downloaded")
def test_real_cohort_fit_is_reproducible():
    cohort = pd.read_parquet(REAL_COHORT)
    m1, cv1, _ = fit_clinical_baseline(cohort, n_boot=100, seed=7)
    m2, cv2, _ = fit_clinical_baseline(cohort, n_boot=100, seed=7)
    assert np.array_equal(m1.coef_, m2.coef_)
    assert cv1.c_index == cv2.c_index
    assert cv1.c_index_ci == cv2.c_index_ci
    assert cv1.brier_score == cv2.brier_score


@pytest.mark.skipif(not REAL_COHORT.exists(), reason="real TCGA cohort not downloaded")
def test_real_cohort_molecular_columns_are_still_empty():
    """Guard: when expression ingestion lands, this test should start failing.

    The clinical-only baseline exists precisely because these are all null.
    """
    cohort = pd.read_parquet(REAL_COHORT)
    for col in ("mgmt_methylated", "idh_mutated", "subtype", "extent_of_resection"):
        assert cohort[col].isna().all(), f"{col} now has data -- refit with it"
