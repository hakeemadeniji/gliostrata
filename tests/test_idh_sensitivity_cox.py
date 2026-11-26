"""Unit tests for scripts/idh_sensitivity_cox.py (no network access).

What is at stake here
---------------------
This script produces the number that decides whether the primary Cox result
is reported as robust to excluding the IDH-mutant patients. Three things
can go wrong and each would produce a plausible-looking JSON:

  1. The wrong patients are excluded. A missing IDH call is not wildtype
     and must never be excluded; only idh_mutated == 1 is.
  2. The two fits are not the same model. Anything other than the
     committed baseline, fitted identically on both cohorts, compares
     apples to oranges.
  3. The flag is wrong. The threshold is strictly greater-than 0.01, and
     `delta_exceeds_threshold` must be a real boolean, not 0/1.

A fake fitter is injected for the logic tests, so they run without a model
fit; the committed estimator itself is exercised separately, and on the
real cohort when it is present.
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
SCRIPT_PATH = ROOT / "scripts" / "idh_sensitivity_cox.py"

REAL_COHORT = ROOT / "data/raw/tcga_gbm/clinical_cohort.parquet"
REAL_JOINED = ROOT / "data/processed/tcga_gbm/molecular_calls_joined.parquet"

#: The clinical-baseline C-index observed on the committed 594-patient
#: cohort. A canary, not a fitted target: if this trips, the cohort, the
#: call tables or the baseline changed -- confirm the change was intended
#: and update it here, never loosen the tolerance to keep it green.
COMMITTED_BASELINE_CINDEX = 0.6743


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_idh_sensitivity_cox_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def S():
    return _load_module()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _cohort(n=300, seed=0):
    """Synthetic cohort with a real age / KPS survival effect."""
    rng = np.random.default_rng(seed)
    age = rng.uniform(30, 75, n)
    kps = rng.choice([60.0, 70.0, 80.0, 90.0, 100.0], n)
    t = rng.exponential(np.exp(0.04 * (age - 55) - 0.03 * (kps - 80) + 2.6))
    c = rng.uniform(2, 80, n)
    return pd.DataFrame(
        {
            "patient_id": [f"P{i:04d}" for i in range(n)],
            "age": age,
            "sex": rng.choice(["M", "F"], n),
            "kps": kps,
            "treatment_arm": rng.choice(
                ["ChemoRadiation", "MonoTherapy", "NoAdjuvant"], n
            ),
            "overall_survival_months": np.minimum(t, c),
            "vital_status": (t <= c).astype(float),
        }
    )


def _joined(cohort, mutant_ids, no_call_ids=()):
    """Joined table: 1 for mutants, <NA> for no-call patients, 0 otherwise."""
    status = pd.Series(0, index=cohort["patient_id"], dtype="Int64")
    status.loc[list(mutant_ids)] = 1
    status.loc[list(no_call_ids)] = pd.NA
    return pd.DataFrame(
        {
            "patient_id": status.index,
            "idh_mutated": status.to_numpy(),
            "has_idh_call": status.notna().to_numpy(),
        }
    )


def _fake_fitter(c_full, c_without):
    """Fitter returning fixed C-indices and recording the frames it sees."""
    seen = []

    def fit(frame, features, n_splits, seed, n_boot):
        seen.append(frame)
        c = c_full if len(seen) == 1 else c_without
        return {
            "cindex": c,
            "ci": [c - 0.03, c + 0.03],
            "n_used": int(len(frame)),
            "n_events": int(frame["vital_status"].sum()),
            "fold_c_indices": [c] * n_splits,
            "feature_names": list(features),
        }

    fit.seen = seen
    return fit


# ---------------------------------------------------------------------------
# 1. The right patients are excluded
# ---------------------------------------------------------------------------
def test_only_a_nonnull_mutant_call_marks_a_patient(S):
    cohort = _cohort(20)
    joined = _joined(cohort, ["P0001", "P0002"], no_call_ids=["P0003", "P0004"])
    assert S.idh_mutant_ids(joined) == {"P0001", "P0002"}


def test_a_missing_call_is_not_treated_as_mutant(S):
    cohort = _cohort(20)
    joined = _joined(cohort, [], no_call_ids=["P0003", "P0004"])
    assert S.idh_mutant_ids(joined) == set()


def test_second_fit_never_sees_a_mutant(S):
    cohort = _cohort(40)
    mutants = [f"P{i:04d}" for i in range(5)]
    no_call = [f"P{i:04d}" for i in range(5, 9)]
    joined = _joined(cohort, mutants, no_call_ids=no_call)
    fit = _fake_fitter(0.70, 0.68)

    r = S.analyze(cohort, joined, fitter=fit)

    first, second = fit.seen
    assert len(first) == 40 and set(first["patient_id"]) == set(cohort["patient_id"])
    assert set(second["patient_id"]) == set(cohort["patient_id"]) - set(mutants)
    assert len(second) == 35
    # The no-call patients are NOT excluded: an absent call proves nothing.
    assert set(no_call) <= set(second["patient_id"])
    assert r["n_full"] == 40
    assert r["n_excluded"] == 5
    assert r["n_without_idh"] == 35


def test_both_fits_get_identical_settings(S):
    cohort = _cohort(40)
    joined = _joined(cohort, ["P0000"])
    calls = []

    def fit(frame, features, n_splits, seed, n_boot):
        calls.append((tuple(features), n_splits, seed, n_boot))
        return _fake_fitter(0.7, 0.7)(frame, features, n_splits, seed, n_boot)

    S.analyze(cohort, joined, n_boot=25, seed=7, n_splits=3, fitter=fit)
    assert calls[0] == calls[1] == (tuple(S.DEFAULT_CLINICAL_FEATURES), 3, 7, 25)


def test_bad_inputs_are_refused(S):
    cohort = _cohort(10)
    joined = _joined(cohort, ["P0000"])
    fit = _fake_fitter(0.7, 0.7)
    with pytest.raises(ValueError, match="required column"):
        S.analyze(cohort, joined.drop(columns=["idh_mutated"]), fitter=fit)
    with pytest.raises(ValueError, match="duplicate patient_id"):
        S.analyze(pd.concat([cohort, cohort]), joined, fitter=fit)
    with pytest.raises(ValueError, match="leaves no cohort"):
        S.analyze(cohort, _joined(cohort, list(cohort["patient_id"])), fitter=fit)


# ---------------------------------------------------------------------------
# 2. Delta and the flag
# ---------------------------------------------------------------------------
def test_delta_is_without_minus_full(S):
    cohort = _cohort(30)
    joined = _joined(cohort, ["P0000"])
    r = S.analyze(cohort, joined, fitter=_fake_fitter(0.700, 0.683))
    assert r["delta"] == pytest.approx(-0.017)
    assert r["delta_exceeds_threshold"] is True

    r = S.analyze(cohort, joined, fitter=_fake_fitter(0.700, 0.703))
    assert r["delta"] == pytest.approx(0.003)
    assert r["delta_exceeds_threshold"] is False


def test_threshold_is_strictly_greater_than(S):
    cohort = _cohort(30)
    joined = _joined(cohort, ["P0000"])
    exactly = S.analyze(cohort, joined, fitter=_fake_fitter(0.700, 0.710))
    assert abs(exactly["delta"]) == pytest.approx(0.01)
    assert exactly["delta_exceeds_threshold"] is False, "0.01 is not > 0.01"

    over = S.analyze(cohort, joined, fitter=_fake_fitter(0.700, 0.711))
    assert over["delta_exceeds_threshold"] is True

    # The sign does not matter: a drop counts the same as a rise.
    drop = S.analyze(cohort, joined, fitter=_fake_fitter(0.700, 0.680))
    assert drop["delta"] < 0 and drop["delta_exceeds_threshold"] is True

    looser = S.analyze(
        cohort, joined, threshold=0.05, fitter=_fake_fitter(0.700, 0.711)
    )
    assert looser["threshold"] == 0.05
    assert looser["delta_exceeds_threshold"] is False


def test_result_contract(S):
    cohort = _cohort(30)
    joined = _joined(cohort, ["P0000"])
    r = S.analyze(cohort, joined, fitter=_fake_fitter(0.700, 0.695))
    assert tuple(list(r)[: len(S.RESULT_KEYS)]) == S.RESULT_KEYS
    assert type(r["delta_exceeds_threshold"]) is bool
    assert type(r["n_excluded"]) is int
    assert len(r["cindex_full_ci"]) == len(r["cindex_without_idh_ci"]) == 2
    assert r["threshold"] == 0.01
    assert json.loads(json.dumps(r)) == r


# ---------------------------------------------------------------------------
# 3. The committed estimator itself
# ---------------------------------------------------------------------------
def test_committed_estimator_recovers_the_planted_effect(S):
    cohort = _cohort(300)
    fitted = S.fit_cohort(cohort, n_boot=20)
    assert fitted["cindex"] > 0.65, "age + KPS effect must be recovered"
    lo, hi = fitted["ci"]
    assert lo <= fitted["cindex"] <= hi
    assert fitted["n_used"] == 300


def test_analyze_with_the_committed_estimator_is_deterministic(S):
    cohort = _cohort(200)
    joined = _joined(cohort, [f"P{i:04d}" for i in range(10)])
    a = S.analyze(cohort, joined, n_boot=20)
    b = S.analyze(cohort, joined, n_boot=20)
    assert a == b
    assert a["delta_exceeds_threshold"] in (True, False)


def test_fit_ignores_cohort_columns_outside_the_design(S):
    """The real cohort carries an all-null idh_mutated column of its own.

    The fits must not depend on it: both frames are reduced to the design
    columns before fitting.
    """
    cohort = _cohort(200)
    plain = S.fit_cohort(cohort, n_boot=5)
    decorated = cohort.assign(idh_mutated=pd.Series([pd.NA] * 200, dtype="Int64"))
    again = S.fit_cohort(decorated, n_boot=5)
    assert plain == again


def test_run_writes_the_result_and_is_idempotent(S, tmp_path):
    cohort = _cohort(150)
    joined = _joined(cohort, ["P0000", "P0001"])
    cohort.to_parquet(tmp_path / "cohort.parquet", index=False)
    joined.to_parquet(tmp_path / "joined.parquet", index=False)
    out = tmp_path / "nested" / "result.json"
    r = S.run(tmp_path / "cohort.parquet", tmp_path / "joined.parquet", out, n_boot=20)
    text = out.read_text(encoding="ascii")
    assert json.loads(text) == r
    assert r["n_excluded"] == 2

    S.run(tmp_path / "cohort.parquet", tmp_path / "joined.parquet", out, n_boot=20)
    assert out.read_text(encoding="ascii") == text


# ---------------------------------------------------------------------------
# 4. The real cohort (skipped cleanly when the data is absent)
# ---------------------------------------------------------------------------
real_data = pytest.mark.skipif(
    not (REAL_COHORT.exists() and REAL_JOINED.exists()),
    reason="TCGA-GBM cohort or joined molecular calls not present",
)


@real_data
def test_real_full_cohort_matches_the_committed_baseline(S):
    """The 'full' fit IS the committed baseline; only its packaging is new."""
    cohort = pd.read_parquet(REAL_COHORT)
    joined = pd.read_parquet(REAL_JOINED)
    r = S.analyze(cohort, joined, n_boot=100)

    assert r["n_full"] == 594
    assert r["n_excluded"] == 23
    assert r["n_without_idh"] == 571
    assert r["n_modelled_full"] == 592, "the 2 patients with no vital status drop"
    assert r["n_modelled_without_idh"] == 569
    assert list(r["features"]) == list(S.DEFAULT_CLINICAL_FEATURES)
    assert "age" in r["design_columns"] and "kps" in r["design_columns"]
    assert r["cindex_full"] == pytest.approx(COMMITTED_BASELINE_CINDEX, abs=0.005)
    lo, hi = r["cindex_full_ci"]
    assert lo < r["cindex_full"] < hi
    lo2, hi2 = r["cindex_without_idh_ci"]
    assert lo2 < r["cindex_without_idh"] < hi2
    assert r["delta"] == pytest.approx(
        r["cindex_without_idh"] - r["cindex_full"], abs=2e-6
    )
    assert type(r["delta_exceeds_threshold"]) is bool
    assert r["delta_exceeds_threshold"] == (abs(r["delta"]) > r["threshold"])


@real_data
def test_real_exclusion_matches_the_joined_table(S):
    cohort = pd.read_parquet(REAL_COHORT)
    joined = pd.read_parquet(REAL_JOINED)
    mutants = set(joined.loc[joined["idh_mutated"] == 1, "patient_id"])
    assert len(mutants) == 23
    assert mutants <= set(cohort["patient_id"]), "every mutant is a cohort patient"
    assert S.idh_mutant_ids(joined) == mutants
