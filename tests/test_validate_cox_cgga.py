"""Unit tests for scripts/validate_cox_cgga.py (no network access).

What is at stake here
---------------------
This script produces the T2 external-validation number. Three things can
go wrong and each would produce a plausible-looking JSON:

  1. CGGA information reaches the fit. The model must be fitted on TCGA
     and scored on CGGA, never the reverse and never jointly.
  2. A missing feature is quietly imputed. CGGA has no KPS at all; that
     must come back as "not evaluable", not as a median.
  3. The threshold flag is wrong. It is true only when the CI LOWER bound
     clears 0.65, not when the point estimate does.

The external encoding is checked against the training encoder itself: the
same rows pushed through `design_external` and through
`build_design_matrix` must give the identical design matrix.

CGGA is NOT read by these tests. It is the external validation set and is
read once, by the script, for the confirmatory evaluation. Every test here
uses synthetic cohorts; the single real-data test uses TCGA (development
data) only, to pin the internal sanity check.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "validate_cox_cgga.py"

REAL_COHORT = ROOT / "data/raw/tcga_gbm/clinical_cohort.parquet"
REAL_JOINED = ROOT / "data/processed/tcga_gbm/molecular_calls_joined.parquet"

#: The committed clinical-baseline C-index (docs/mgmt_idh_calling.md 6.3).
COMMITTED_BASELINE_CINDEX = 0.6735


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_validate_cox_cgga_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def S():
    return _load_module()


# ---------------------------------------------------------------------------
# Synthetic cohorts
# ---------------------------------------------------------------------------
ARMS = ["ChemoRadiation", "MonoTherapy", "NoAdjuvant"]


def _cohort(n=400, seed=0, prefix="T", kps=True):
    """Synthetic patients with real survival effects for every shared
    feature, so a correctly transferred model discriminates."""
    rng = np.random.default_rng(seed)
    age = rng.uniform(30, 75, n)
    sex = rng.choice(["M", "F"], n)
    arm = rng.choice(ARMS, n, p=[0.70, 0.25, 0.05])
    idh = rng.choice([0, 1], n, p=[0.90, 0.10])
    mgmt = rng.choice([0, 1], n, p=[0.55, 0.45])
    k = rng.choice([60.0, 70.0, 80.0, 90.0, 100.0], n) if kps else np.zeros(n)
    eta = (
        0.045 * (age - 55)
        + 0.25 * (sex == "M")
        + 0.20 * (arm != "ChemoRadiation")
        - 0.80 * idh
        - 0.35 * mgmt
        - (0.03 * (k - 80.0) if kps else 0.0)
    )
    # Hazard proportional to exp(eta): mean survival is exp(1.9 - eta) months,
    # so a higher eta shortens survival. Censoring is uniform at 3-90 months.
    t = rng.exponential(np.exp(1.9 - eta))
    censor = rng.uniform(3, 90, n)
    return pd.DataFrame(
        {
            "patient_id": [f"{prefix}{i:04d}" for i in range(n)],
            "age": age,
            "sex": sex,
            "kps": k,
            "treatment_arm": arm,
            "idh_mutated": idh.astype(float),
            "mgmt_methylated": mgmt.astype(float),
            "overall_survival_months": np.minimum(t, censor),
            "vital_status": (t <= censor).astype(int),
        }
    )


def _as_cgga(frame: pd.DataFrame, seed=1) -> pd.DataFrame:
    """Reshape a synthetic cohort to CGGA's dtypes and missingness pattern:
    no KPS at all, nullable molecular calls, a few missing outcomes."""
    rng = np.random.default_rng(seed)
    out = frame.copy()
    out["kps"] = pd.Series([None] * len(out), dtype="object")
    out["extent_of_resection"] = pd.Series([None] * len(out), dtype="object")
    out["sex"] = out["sex"].astype("string")
    out["treatment_arm"] = out["treatment_arm"].mask(rng.random(len(out)) < 0.05)
    out["idh_mutated"] = (
        out["idh_mutated"].mask(rng.random(len(out)) < 0.03).astype("Int64")
    )
    out["mgmt_methylated"] = (
        out["mgmt_methylated"].mask(rng.random(len(out)) < 0.06).astype("Int64")
    )
    out["vital_status"] = out["vital_status"].astype("Int64")
    lost = rng.random(len(out)) < 0.02
    out.loc[lost, "overall_survival_months"] = np.nan
    out.loc[lost, "vital_status"] = pd.NA
    return out


# ---------------------------------------------------------------------------
# Feature availability: CGGA KPS is absent, and that is reported, not imputed
# ---------------------------------------------------------------------------
class _NoPredictModel:
    """Stands in for a fitted model; any scoring attempt is a failure."""

    def predict_log_hazard(self, X):  # pragma: no cover - must not be reached
        raise AssertionError("an unscorable variant must not be scored")


def test_absent_kps_is_reported_not_imputed(S):
    tcga = _cohort(400, seed=0, prefix="T", kps=True)
    cgga = _as_cgga(_cohort(200, seed=1, prefix="C", kps=True))
    fitted = S.fit_variant(tcga, S.VARIANT_1, n_boot=4)
    fitted.model = _NoPredictModel()

    ev = S.evaluate_external(fitted, cgga, n_boot=5)

    assert ev["external_evaluable"] is False
    assert ev["absent_in_cgga"] == ["kps"]
    assert "KPS" in ev["external_reason"] or "kps" in ev["external_reason"]
    assert ev["external_cindex"] is None
    assert ev["external_ci_lower"] is None
    assert ev["clears_t2_threshold"] is False
    assert ev["n_cgga_eval"] == 0
    # Availability is a clean count, and it is zero for KPS only.
    avail = ev["cgga_feature_availability"]
    assert avail["kps"] == 0
    assert all(avail[f] > 0 for f in ("age", "sex", "treatment_arm"))


def test_absent_column_and_all_null_column_are_both_absent(S):
    frame = _cohort(50)
    frame["all_null"] = pd.Series([None] * 50, dtype="object")
    avail = S.feature_availability(frame, ["age", "all_null", "not_a_column"])
    assert avail == {"age": 50, "all_null": 0, "not_a_column": 0}


def test_partial_missingness_is_complete_case_never_imputed(S):
    frame = _as_cgga(_cohort(300, seed=3, prefix="C"))
    kept, dropped = S.complete_cases(frame, S.VARIANT_2)
    assert len(kept) + dropped["n_dropped_missing_feature"] + \
        dropped["n_dropped_missing_outcome"] == len(frame)
    assert kept[list(S.VARIANT_2)].notna().all().all()
    assert kept[S.TIME_COL].notna().all() and kept[S.EVENT_COL].notna().all()


# ---------------------------------------------------------------------------
# Fit-then-evaluate, and no CGGA row ever reaches a fit
# ---------------------------------------------------------------------------
def test_fit_then_evaluate_order_and_fit_receives_only_tcga(S, monkeypatch):
    tcga = _cohort(300, seed=0, prefix="T")
    cgga_rows = _as_cgga(_cohort(150, seed=1, prefix="C"))
    tcga_ids = set(tcga["patient_id"])
    cgga_ids = set(cgga_rows["patient_id"])
    assert not (tcga_ids & cgga_ids)
    events: list[tuple[str, str, set]] = []

    def spy_fit(frame, features, **kw):
        events.append(("fit", tuple(features), set(frame["patient_id"])))
        return S.fit_variant(frame, features, **{**kw, "n_boot": 3})

    def spy_baseline(frame, **kw):
        events.append(("baseline", tuple(S.VARIANT_1), set(frame["patient_id"])))
        return S.committed_baseline_internal(frame, **{**kw, "n_boot": 3})

    real_eval = S.evaluate_external

    def spy_eval(fitted, cgga, **kw):
        events.append(("evaluate", tuple(fitted.features), set(cgga["patient_id"])))
        return real_eval(fitted, cgga, **{**kw, "n_boot": 5})

    monkeypatch.setattr(S, "evaluate_external", spy_eval)
    result = S.run(tcga, cgga_rows, n_boot=5, fit_fn=spy_fit, baseline_fn=spy_baseline)

    # Every evaluation happens after a fit of the same variant.
    seen_fits = set()
    for kind, feats, _ in events:
        if kind in ("fit", "baseline"):
            seen_fits.add(feats)
        else:
            assert kind == "evaluate" and feats in seen_fits

    # Every fit saw exactly the TCGA rows and not one CGGA row.
    for kind, _, ids in events:
        if kind in ("fit", "baseline"):
            assert ids == tcga_ids, f"{kind} was not fitted on TCGA rows alone"
            assert not (ids & cgga_ids)
        else:
            assert ids <= cgga_ids and not (ids & tcga_ids)

    assert events[0][0] == "baseline"           # variant 1 first
    assert ("fit", S.VARIANT_2, tcga_ids) in events
    assert ("fit", S.VARIANT_3, tcga_ids) in events
    assert result["variant_2_transferable"]["n_cgga_eval"] > 0


def test_fitted_variant_records_tcga_patient_ids(S):
    tcga = _cohort(200, seed=0, prefix="T")
    fitted = S.fit_variant(tcga, S.VARIANT_2, n_boot=3)
    assert fitted.fit_patient_ids == set(tcga["patient_id"])
    assert fitted.n_fit == 200 == fitted.n_modelled


# ---------------------------------------------------------------------------
# The external cohort is encoded exactly as the training cohort was
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("features", [("age", "sex", "treatment_arm"),
                                      ("age", "sex", "treatment_arm",
                                       "idh_mutated", "mgmt_methylated")])
def test_external_encoding_equals_training_encoding(S, features):
    from gliostrata.prediction.cox import build_design_matrix

    tcga = _cohort(600, seed=0, prefix="T")
    external = _cohort(150, seed=5, prefix="X")
    fitted = S.fit_variant(tcga, features, n_boot=3)

    X_ext, unseen = S.design_external(external, fitted)
    assert not unseen.any()

    # Same external rows, run through the TRAINING encoder with the
    # training reference levels: must be identical, column for column.
    refs = {
        feat: next(lvl for lvl in levels
                   if f"{feat}_{lvl}" not in fitted.design_columns)
        for feat, levels in fitted.categorical_levels.items()
    }
    X_both, _, _, report = build_design_matrix(
        external, features, impute="drop", reference_levels=refs
    )
    assert list(report.feature_names) == fitted.design_columns
    np.testing.assert_array_equal(X_ext, X_both)


def test_external_encoding_has_training_column_order(S):
    tcga = _cohort(400, seed=0, prefix="T")
    fitted = S.fit_variant(tcga, S.VARIANT_3, n_boot=3)
    X_ext, _ = S.design_external(_cohort(40, seed=9, prefix="X"), fitted)
    assert X_ext.shape == (40, len(fitted.design_columns))
    assert "idh_mutated" in fitted.design_columns
    assert "mgmt_methylated" in fitted.design_columns


def test_unseen_category_is_excluded_and_counted_never_reference_coded(S):
    tcga = _cohort(400, seed=0, prefix="T")
    external = _as_cgga(_cohort(200, seed=2, prefix="C"))
    external.loc[[3, 7], "treatment_arm"] = "Experimental"
    # Make sure those two rows survive every other exclusion.
    external.loc[[3, 7], "overall_survival_months"] = 12.0
    external.loc[[3, 7], S.EVENT_COL] = 1
    fitted = S.fit_variant(tcga, S.VARIANT_2, n_boot=3)

    X, unseen = S.design_external(
        S.complete_cases(external, S.VARIANT_2)[0], fitted)
    # "Experimental" is not silently all-zero (which would mean the
    # reference level): the row is flagged.
    assert int(unseen.sum()) == 2

    ev = S.evaluate_external(fitted, external, n_boot=5)
    n_expected = len(S.complete_cases(external, S.VARIANT_2)[0]) - int(unseen.sum())
    assert ev["n_cgga_excluded"]["n_dropped_category_unseen_in_tcga"] == int(unseen.sum())
    assert ev["n_cgga_eval"] == n_expected
    assert ev["n_cgga_eval"] < 200


def test_numeric_missing_value_is_refused_not_imputed(S):
    tcga = _cohort(200, seed=0, prefix="T")
    fitted = S.fit_variant(tcga, S.VARIANT_2, n_boot=3)
    external = _cohort(20, seed=4, prefix="X")
    external.loc[0, "age"] = np.nan
    with pytest.raises(ValueError, match="never imputed"):
        S.design_external(external, fitted)


# ---------------------------------------------------------------------------
# A variant with no evaluable patient reports null; it does not crash
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("empty", [False, True])
def test_no_evaluable_patients_reports_null(S, empty):
    tcga = _cohort(200, seed=0, prefix="T")
    fitted = S.fit_variant(tcga, S.VARIANT_2, n_boot=3)
    external = _as_cgga(_cohort(100, seed=6, prefix="C"))
    if empty:
        external = external.iloc[0:0]
    else:
        external["overall_survival_months"] = np.nan
        external[S.EVENT_COL] = pd.NA

    ev = S.evaluate_external(fitted, external, n_boot=5)

    assert ev["external_evaluable"] is False
    assert ev["external_cindex"] is None
    assert ev["external_ci_lower"] is None and ev["external_ci_upper"] is None
    assert ev["clears_t2_threshold"] is False
    assert ev["n_cgga_eval"] == 0


# ---------------------------------------------------------------------------
# Bootstrap determinism and the threshold flag
# ---------------------------------------------------------------------------
def test_bootstrap_same_seed_gives_identical_ci(S):
    frame = _cohort(250, seed=11, prefix="B")
    rng = np.random.default_rng(0)
    risk = 0.04 * (frame["age"] - 55) + rng.normal(0, 0.3, len(frame))
    args = (frame[S.TIME_COL].to_numpy(), frame[S.EVENT_COL].to_numpy(), risk)
    first = S.bootstrap_cindex(*args, n_boot=300, seed=42)
    second = S.bootstrap_cindex(*args, n_boot=300, seed=42)
    assert first == second
    assert first[0] <= first[1]
    assert first[2] == 300
    # A different seed draws different resamples, hence (almost surely) a
    # different interval.
    other = S.bootstrap_cindex(*args, n_boot=300, seed=7)
    assert other[:2] != first[:2]


def test_bootstrap_brackets_the_point_estimate(S):
    from gliostrata.prediction.cox import concordance_index

    frame = _cohort(300, seed=12, prefix="B")
    rng = np.random.default_rng(1)
    risk = 0.04 * (frame["age"] - 55) + rng.normal(0, 0.3, len(frame))
    time = frame[S.TIME_COL].to_numpy()
    event = frame[S.EVENT_COL].to_numpy()
    lo, hi, _ = S.bootstrap_cindex(time, event, risk, n_boot=400, seed=42)
    assert lo < concordance_index(time, event, risk) < hi


def test_clears_threshold_requires_the_ci_lower_bound(S):
    # Point estimate above 0.65, lower bound below it: NOT cleared.
    assert S.clears_t2_threshold(0.66, 0.62) is False
    # Point estimate itself below 0.65: NOT cleared even with a wide CI.
    assert S.clears_t2_threshold(0.649, 0.70) is False
    # Both clear.
    assert S.clears_t2_threshold(0.70, 0.66) is True
    # The pre-registration says >= 0.65: exactly on the line counts.
    assert S.clears_t2_threshold(0.65, 0.65) is True
    # One tick below does not.
    assert S.clears_t2_threshold(0.70, 0.6499) is False
    # Nothing to judge.
    assert S.clears_t2_threshold(None, None) is False
    assert S.clears_t2_threshold(0.70, None) is False


def test_verdict_labels_follow_the_point_estimate_and_ci(S):
    def block(c, lo):
        return {"external_cindex": c, "external_ci_lower": lo,
                "external_ci_upper": 0.75, "n_cgga_eval": 300}

    assert S.t2_verdict(block(0.70, 0.66))[0] == S.VERDICT_MET
    assert S.t2_verdict(block(0.66, 0.62))[0] == S.VERDICT_CI
    assert S.t2_verdict(block(0.60, 0.55))[0] == S.VERDICT_SHORT
    assert S.t2_verdict(block(None, None))[0] == S.VERDICT_SHORT
    label, explanation = S.t2_verdict(block(0.60, 0.55))
    assert "0.6000" in explanation and "0.5500" in explanation
    assert label in (S.VERDICT_MET, S.VERDICT_CI, S.VERDICT_SHORT)


# ---------------------------------------------------------------------------
# The pipeline end to end on synthetic data
# ---------------------------------------------------------------------------
def test_fitted_tcga_model_transfers_to_a_synthetic_external_cohort(S):
    """Both cohorts are drawn from the same truth, so a correctly
    transferred model must discriminate externally, and recover the
    planted effect directions."""
    tcga = _cohort(600, seed=0, prefix="T")
    external = _as_cgga(_cohort(500, seed=21, prefix="C"))

    fitted = S.fit_variant(tcga, S.VARIANT_3, n_boot=3)
    ev = S.evaluate_external(fitted, external, n_boot=100, seed=42)

    assert ev["external_cindex"] > 0.60
    assert ev["external_ci_lower"] < ev["external_cindex"] < ev["external_ci_upper"]
    assert ev["n_bootstrap_valid"] == 100
    hr = fitted.hazard_ratios
    assert hr["age"] > 1.0            # older -> worse
    assert hr["idh_mutated"] < 1.0    # IDH-mutant -> better
    assert hr["mgmt_methylated"] < 1.0


def test_molecular_variant_fit_excludes_patients_without_a_call(S):
    """A missing call is not wildtype: those patients are not fitted on."""
    tcga = _cohort(500, seed=0, prefix="T")
    tcga["idh_mutated"] = tcga["idh_mutated"].astype("Float64")
    tcga.loc[:49, "idh_mutated"] = pd.NA
    fitted = S.fit_variant(tcga, S.VARIANT_3, n_boot=3)
    assert fitted.n_fit == 450
    assert fitted.fit_patient_ids == set(tcga.loc[50:, "patient_id"])


def test_one_row_per_patient_keeps_first_row(S):
    frame = pd.DataFrame({
        "patient_id": ["A", "B", "A", "C", "B"],
        "v": [1, 2, 3, 4, 5],
    })
    kept, n_dropped = S.one_row_per_patient(frame)
    assert list(kept["patient_id"]) == ["A", "B", "C"]
    assert list(kept["v"]) == [1, 2, 4]
    assert n_dropped == 2


def test_run_produces_the_specified_contract(S):
    tcga = _cohort(400, seed=0, prefix="T")
    rows = _as_cgga(_cohort(200, seed=8, prefix="C"))
    # A duplicated id, as in the real CGGA.
    rows = pd.concat([rows, rows.iloc[[0, 1]]], ignore_index=True)

    result = S.run(tcga, rows, n_boot=40, seed=42)

    # The three blocks and the top-level keys exist under the exact names.
    v1 = result["variant_1_baseline_committed"]
    v2 = result["variant_2_transferable"]
    v3 = result["variant_3_molecular"]
    for key in ("features", "n_tcga_fit", "internal_cindex",
                "internal_cindex_matches_committed", "external_evaluable",
                "external_reason"):
        assert key in v1
    for block in (v2, v3):
        for key in ("features", "n_tcga_fit", "n_cgga_eval", "internal_cindex",
                    "external_cindex", "external_ci_lower", "external_ci_upper",
                    "clears_t2_threshold"):
            assert key in block
    assert "harmonisation_caveat" in v3
    assert result["t2_threshold"] == 0.65
    assert result["t2_verdict"] in (S.VERDICT_MET, S.VERDICT_CI, S.VERDICT_SHORT)
    assert result["t2_verdict_explanation"].endswith(".")

    # Variant 1: the committed feature list, honestly not evaluable.
    assert v1["features"] == ["age", "sex", "kps", "treatment_arm"]
    assert v1["external_evaluable"] is False
    assert v1["external_reason"] == (
        "CGGA has zero KPS records; the committed baseline cannot be "
        "externally validated"
    )
    assert type(v1["internal_cindex_matches_committed"]) is bool

    # Variants 2 and 3: the flag follows the CI lower bound, not the point.
    for block in (v2, v3):
        assert type(block["clears_t2_threshold"]) is bool
        assert block["clears_t2_threshold"] == (
            block["external_cindex"] >= 0.65
            and block["external_ci_lower"] >= 0.65
        )
        assert block["external_ci_lower"] <= block["external_cindex"] \
            <= block["external_ci_upper"]
    assert v2["features"] == ["age", "sex", "treatment_arm"]
    assert v3["features"] == ["age", "sex", "treatment_arm", "idh_mutated",
                              "mgmt_methylated"]

    # The complete-case molecular variant cannot have more patients than
    # variant 2: it needs everything variant 2 needs, plus two calls.
    assert v3["n_cgga_eval"] <= v2["n_cgga_eval"]
    assert v3["n_tcga_fit"] <= v2["n_tcga_fit"] == 400

    # One row per patient: the duplicated rows are not scored twice.
    assert result["cgga_cohort"]["n_duplicate_rows_dropped"] == 2
    assert result["cgga_cohort"]["n_unique_patients"] == 200
    assert v2["n_cgga_eval"] <= 200

    # The verdict is decided by the headline variant.
    c, lo = v2["external_cindex"], v2["external_ci_lower"]
    expected = (S.VERDICT_MET if lo >= 0.65
                else S.VERDICT_CI if c >= 0.65 else S.VERDICT_SHORT)
    assert result["t2_verdict"] == expected
    assert result["headline_variant"] == "variant_2_transferable"


def test_run_result_is_ascii_json(S):
    tcga = _cohort(200, seed=0, prefix="T")
    result = S.run(tcga, _as_cgga(_cohort(100, seed=3, prefix="C")), n_boot=5)
    text = json.dumps(result, indent=2)
    text.encode("ascii")  # raises if any non-ASCII character slipped in
    assert json.loads(text) == result


def test_dry_run_fits_nothing_and_computes_no_metric(S, monkeypatch):
    def forbidden(*a, **k):  # pragma: no cover - must not be reached
        raise AssertionError("dry run must not fit or score anything")

    monkeypatch.setattr(S, "fit_variant", forbidden)
    monkeypatch.setattr(S, "committed_baseline_internal", forbidden)
    monkeypatch.setattr(S, "evaluate_external", forbidden)

    report = S.dry_run(_cohort(300, seed=0, prefix="T"),
                       _as_cgga(_cohort(150, seed=2, prefix="C")))
    v1, v2, v3 = (report["variants"][k] for k, _ in S.VARIANTS)
    assert v1["externally_evaluable"] is False and v1["absent_in_cgga"] == ["kps"]
    assert v2["externally_evaluable"] and v2["n_tcga_fit"] == 300
    assert v3["n_cgga_evaluable"] <= v2["n_cgga_evaluable"] <= 150


# ---------------------------------------------------------------------------
# The real-data tests: TCGA only (development data). CGGA is never read.
# ---------------------------------------------------------------------------
needs_tcga = pytest.mark.skipif(
    not (REAL_COHORT.is_file() and REAL_JOINED.is_file()),
    reason="TCGA cohort or joined molecular calls not present",
)


@needs_tcga
def test_real_tcga_internal_sanity_check_reproduces_the_committed_baseline(S):
    """Variant 1's internal C-index must be the committed number. If this
    trips, the Cox configuration or feature set has drifted and every T2
    number is meaningless."""
    tcga = S.load_tcga(pd.read_parquet(REAL_COHORT), pd.read_parquet(REAL_JOINED))
    base = S.committed_baseline_internal(tcga, n_boot=20, seed=42)
    assert base["n_fit"] == 594
    assert base["n_modelled"] == 592          # 2 patients lack a vital status
    assert abs(base["internal_cindex"] - COMMITTED_BASELINE_CINDEX) <= 0.005


@needs_tcga
def test_real_tcga_load_keeps_missing_calls_missing(S):
    """A missing IDH / MGMT call stays missing; it is never a wildtype."""
    tcga = S.load_tcga(pd.read_parquet(REAL_COHORT), pd.read_parquet(REAL_JOINED))
    assert len(tcga) == 594
    assert tcga["patient_id"].is_unique
    # Manifest counts from the molecular-calls join: 368 IDH, 415 MGMT.
    assert int(tcga["idh_mutated"].notna().sum()) == 368
    assert int(tcga["mgmt_methylated"].notna().sum()) == 415
    assert int(tcga["idh_mutated"].sum()) == 23   # mutants, not wildtypes
    cc2, _ = S.complete_cases(tcga, S.VARIANT_2)
    cc3, _ = S.complete_cases(tcga, S.VARIANT_3)
    assert len(cc2) == 570 and len(cc3) == 298

