"""Unit tests for scripts/validate_brier_cgga.py (no network access).

What is at stake here
---------------------
This script produces the second T2 external number, the Brier score at 18
months. Four things can go wrong and each would produce a plausible JSON:

  1. A different patient set from the C-index row. The scoring set is built by
     the C-index script's own functions, and a test checks that both scripts
     evaluate the same number of patients on the same cohort.
  2. A verdict rule written in the script. The verdict must come from the
     pre-registered registry (key T2_brier_18mo), and a test pins that row.
  3. A silently re-run evaluation. `main` must refuse to overwrite a result.
  4. A TCGA fit that is not the one the C-index row was computed with. The
     fit-identity check is tested against a stand-in result file.

CGGA is NOT read by these tests. Every cohort below is synthetic. The
synthetic builders are imported from the C-index script's test module.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from gliostrata.evaluation.thresholds import CIBound, Direction, get_threshold
from gliostrata.prediction.cox import brier_score_ipcw
from tests.test_validate_cox_cgga import _as_cgga, _cohort

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "validate_brier_cgga.py"


def _load_module():
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_validate_brier_cgga_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def S():
    return _load_module()


def _fitted(S, tcga, features=None):
    return S.CX.fit_variant(tcga, tuple(features or S.CX.VARIANT_2), n_boot=5)


def _tcga_external(seed=1, n=150):
    """TCGA training frame and a CGGA-shaped external frame, synthetic."""
    return _cohort(300, seed=0), _as_cgga(_cohort(n, seed=seed, prefix="C"), seed=seed)


# ---------------------------------------------------------------------------
# The registered row, and a verdict that comes from the registry
# ---------------------------------------------------------------------------
def test_the_row_is_the_registered_one(S):
    row = get_threshold(S.METRIC_KEY)
    assert row.test_id == "T2"
    assert row.value == 0.25 and row.direction is Direction.LTE
    assert row.ci_bound is CIBound.UPPER
    assert "18" in row.metric and "Brier" in row.metric
    assert S.HORIZON == 18.0


@pytest.mark.parametrize(
    "point_lo_hi_status",
    [
        ((0.20, 0.15, 0.24), "PASS"),        # point and upper bound both <= 0.25
        ((0.20, 0.15, 0.26), "INCONCLUSIVE"),  # point clears, upper bound does not
        ((0.25, 0.20, 0.25), "PASS"),        # exactly on the threshold clears
        ((0.27, 0.22, 0.30), "FAIL"),        # the point estimate itself misses
    ],
)
def test_verdict_is_the_registry_decision(S, point_lo_hi_status):
    (point, lo, hi), status = point_lo_hi_status
    verdict = S.registered_verdict(point, lo, hi)
    assert verdict["status"] == status
    assert verdict["passed"] is (status == "PASS")


# ---------------------------------------------------------------------------
# The metric
# ---------------------------------------------------------------------------
def test_score_equals_the_library_implementation(S):
    rng = np.random.default_rng(3)
    time = rng.uniform(1, 40, 120)
    event = (time < 18).astype(int)
    surv = rng.uniform(0.05, 0.6, 120)
    out = S.brier_with_ci(time, event, surv, 18.0, n_boot=60)
    assert out["score"] == brier_score_ipcw(time, event, surv, 18.0)
    assert out["ci_lower"] <= out["score"] <= out["ci_upper"]
    assert out["n_valid"] == 60


def test_bootstrap_is_seed_deterministic(S):
    rng = np.random.default_rng(5)
    time, event, surv = rng.uniform(1, 40, 80), (rng.random(80) < 0.6), rng.uniform(0, 0.7, 80)
    a = S.brier_with_ci(time, event.astype(int), surv, 18.0, n_boot=40, seed=1)
    b = S.brier_with_ci(time, event.astype(int), surv, 18.0, n_boot=40, seed=1)
    assert (a["ci_lower"], a["ci_upper"]) == (b["ci_lower"], b["ci_upper"])


def test_perfect_and_worst_predictions_without_censoring(S):
    time, event = np.array([5.0, 6.0]), np.array([1, 1])
    assert S.brier_with_ci(time, event, np.zeros(2), 18.0, n_boot=5)["score"] == 0.0
    assert S.brier_with_ci(time, event, np.ones(2), 18.0, n_boot=5)["score"] == 1.0


def test_kaplan_meier_survival_handles_censoring(S):
    # death at 5 (3 of 4 left), censored at 10, death at 20 (2 at risk): 3/4 * 1/2
    time = np.array([5.0, 10.0, 20.0, 30.0])
    event = np.array([1, 0, 1, 1])
    assert S.km_survival_at(time, event, 2.0) == 1.0
    assert S.km_survival_at(time, event, 18.0) == pytest.approx(3 / 4)
    assert S.km_survival_at(time, event, 25.0) == pytest.approx(3 / 8)
    assert S.km_survival_at(time, event, 35.0) == 0.0


def test_null_model_brier(S):
    # everyone dies before the horizon: the null model predicts 0 and is right
    assert S.null_model_brier(np.array([5.0, 6.0, 7.0]), np.array([1, 1, 1]), 18.0) == 0.0
    # half die before, half after: the null predicts 0.5, every term is 0.25
    got = S.null_model_brier(np.array([5.0, 6.0, 30.0, 40.0]), np.array([1, 1, 1, 1]), 18.0)
    assert got == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# Same patients as the C-index row
# ---------------------------------------------------------------------------
def test_scoring_set_is_the_cindex_scripts_patient_set(S):
    tcga, cgga_rows = _tcga_external()
    fitted = _fitted(S, tcga)
    cgga, _ = S.CX.one_row_per_patient(cgga_rows)
    frame, X, excluded = S.scoring_set(fitted, cgga)
    cindex = S.CX.evaluate_external(fitted, cgga, n_boot=5)  # synthetic data only
    assert len(frame) == cindex["n_cgga_eval"] > 0
    assert excluded == cindex["n_cgga_excluded"]
    assert X.shape == (len(frame), len(fitted.design_columns))


def test_unseen_category_is_excluded_and_counted(S):
    tcga, cgga_rows = _tcga_external()
    fitted = _fitted(S, tcga)
    cgga, _ = S.CX.one_row_per_patient(cgga_rows)
    baseline = S.scoring_set(fitted, cgga)[2]["n_dropped_category_unseen_in_tcga"]
    cgga.loc[cgga.index[:2], "treatment_arm"] = "Experimental"
    frame, _, excluded = S.scoring_set(fitted, cgga)
    assert excluded["n_dropped_category_unseen_in_tcga"] == baseline + 2
    assert not (frame["treatment_arm"] == "Experimental").any()


# ---------------------------------------------------------------------------
# The whole run
# ---------------------------------------------------------------------------
def test_run_produces_the_specified_contract(S):
    tcga, cgga_rows = _tcga_external()
    result = S.run(tcga, cgga_rows, n_boot=20, cindex_result="does/not/exist.json")

    assert result["metric_key"] == S.METRIC_KEY and result["horizon_months"] == 18.0
    assert result["headline_variant"] == "variant_2_transferable"
    assert result["cgga_read"]["read_number"] == 2
    assert result["pre_specification"]["rules_fixed_before_first_cgga_brier_evaluation"]

    v1 = result["variant_1_baseline_committed"]
    assert not v1["external_evaluable"] and "kps" in v1["absent_in_cgga"]

    v2 = result["variant_2_transferable"]
    assert v2["external_evaluable"] and v2["n_cgga_eval"] > 0
    assert v2["ci_lower"] <= v2["brier_18mo"] <= v2["ci_upper"]
    assert 0.0 <= v2["brier_18mo"] <= 1.0
    assert set(v2["descriptive_other_horizons"]) == {"brier_12mo", "brier_24mo"}
    assert v2["n_cgga_deaths"] <= v2["n_cgga_eval"]
    assert result["variant_3_molecular"]["harmonisation_caveat"]
    assert result["verdict"]["status"] in {"PASS", "FAIL", "INCONCLUSIVE"}
    assert result["fit_identity_check"]["checked"] is False
    assert result["method"]["n_resamples_requested"] == 20


def test_verdict_follows_the_headline_numbers(S):
    tcga, cgga_rows = _tcga_external()
    result = S.run(tcga, cgga_rows, n_boot=20, cindex_result="does/not/exist.json")
    v2 = result["variant_2_transferable"]
    again = S.registered_verdict(v2["brier_18mo"], v2["ci_lower"], v2["ci_upper"])
    assert again["status"] == result["verdict"]["status"]
    if v2["brier_18mo"] > 0.25:
        assert result["verdict"]["status"] == "FAIL"


def test_result_is_ascii_json(S, tmp_path):
    tcga, cgga_rows = _tcga_external()
    result = S.run(tcga, cgga_rows, n_boot=10, cindex_result="x.json")
    text = json.dumps(result, indent=2) + "\n"
    text.encode("ascii")
    (tmp_path / "r.json").write_text(text, encoding="ascii")


def test_interval_that_excludes_its_point_is_widened_and_flagged(S, monkeypatch):
    tcga, _ = _tcga_external()
    fitted = _fitted(S, tcga)

    def shifted(*args, **kwargs):
        return {"score": 0.20, "ci_lower": 0.21, "ci_upper": 0.30, "n_valid": 10}

    monkeypatch.setattr(S, "brier_with_ci", shifted)
    cgga, _ = S.CX.one_row_per_patient(_as_cgga(_cohort(80, seed=2, prefix="C")))
    block, point, lo, hi = S.evaluate_variant(fitted, cgga, n_boot=5)
    assert block["ci_widened_to_contain_point"] is True
    assert (point, lo, hi) == (0.20, 0.20, 0.30)
    assert S.registered_verdict(point, lo, hi)["status"] == "INCONCLUSIVE"


# ---------------------------------------------------------------------------
# Fit identity against the first run, and the two guards in main()
# ---------------------------------------------------------------------------
def test_fit_identity_check(S, tmp_path):
    blocks = {"variant_2_transferable": {"internal_cindex": 0.65, "n_cgga_eval": 283}}
    path = tmp_path / "first.json"
    assert S.compare_with_cindex_run(path, blocks)["checked"] is False

    path.write_text(json.dumps({"variant_2_transferable":
                                {"internal_cindex": 0.65, "n_cgga_eval": 283}}))
    assert S.compare_with_cindex_run(path, blocks)["matches_cindex_run"] is True

    path.write_text(json.dumps({"variant_2_transferable":
                                {"internal_cindex": 0.66, "n_cgga_eval": 283}}))
    assert S.compare_with_cindex_run(path, blocks)["matches_cindex_run"] is False

    path.write_text(json.dumps({"variant_2_transferable":
                                {"internal_cindex": 0.65, "n_cgga_eval": 240}}))
    assert S.compare_with_cindex_run(path, blocks)["matches_cindex_run"] is False


def _write_inputs(tmp_path):
    tcga = _cohort(120)
    cgga = _as_cgga(_cohort(60, seed=1, prefix="C"))
    tcga.to_parquet(tmp_path / "tcga.parquet")
    tcga[["patient_id", "idh_mutated", "mgmt_methylated"]].to_parquet(
        tmp_path / "joined.parquet")
    cgga.to_parquet(tmp_path / "cgga.parquet")
    return ["--tcga", str(tmp_path / "tcga.parquet"), "--joined",
            str(tmp_path / "joined.parquet"), "--cgga", str(tmp_path / "cgga.parquet"),
            "--cindex-result", str(tmp_path / "first.json")]


def test_dry_run_computes_no_metric(S, monkeypatch, tmp_path, capsys):
    def forbidden(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("a dry run must not compute the Brier score")

    monkeypatch.setattr(S, "brier_score_ipcw", forbidden)
    S.main(_write_inputs(tmp_path) + ["--dry-run", "--out", str(tmp_path / "r.json")])
    assert "no metric computed" in capsys.readouterr().out
    assert not (tmp_path / "r.json").exists()


def test_main_refuses_to_overwrite_an_existing_result(S, monkeypatch, tmp_path):
    out = tmp_path / "r.json"
    out.write_text("{}")

    def forbidden(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("must stop before fitting or scoring anything")

    monkeypatch.setattr(S.CX, "fit_variant", forbidden)
    with pytest.raises(SystemExit, match="Refusing to overwrite"):
        S.main(_write_inputs(tmp_path) + ["--out", str(out)])
    assert out.read_text() == "{}"


def test_main_writes_a_result_when_the_run_is_legitimate(S, tmp_path):
    args = _write_inputs(tmp_path)
    S.main(args + ["--out", str(tmp_path / "r.json"), "--n-boot", "10"])
    result = json.loads((tmp_path / "r.json").read_text(encoding="ascii"))
    assert result["verdict"]["status"] in {"PASS", "FAIL", "INCONCLUSIVE"}
    assert result["fit_identity_check"]["checked"] is False
