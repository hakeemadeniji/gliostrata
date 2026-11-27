"""Tests for scripts/rerun_tcga_internal_cv.py (no network; synthetic cohorts).

What is at stake
----------------
The script reports the internal C-index of the committed baseline under the
earlier whole-frame imputation and under fold-wise imputation. Four things
would make its table misleading, and each has a test:

  1. The "new" medians are not those of each fold's training rows.
  2. The two procedures disagree when there is nothing to impute.
  3. The wrong patients are left out of the IDH sensitivity cohort.
  4. The script reads the held-out cohort, or overwrites a record by accident.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "rerun_tcga_internal_cv.py"
#: A compact feature set: the procedures do not depend on how many there are.
FEATURES = ("age", "kps", "treatment_arm")


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_rerun_tcga_internal_cv_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def script():
    return _load_module()


def _cohort(n=90, seed=0, missing=0.25):
    """Synthetic cohort with a KPS that is sometimes missing."""
    rng = np.random.default_rng(seed)
    kps = rng.choice([60.0, 70.0, 80.0, 90.0], n).astype(float)
    kps[rng.random(n) < missing] = np.nan
    return pd.DataFrame(
        {
            "patient_id": [f"P{i:03d}" for i in range(n)],
            "age": np.round(rng.uniform(30, 75, n), 1),
            "kps": kps,
            "treatment_arm": rng.choice(["ChemoRadiation", "MonoTherapy"], n),
            "overall_survival_months": np.round(rng.exponential(15, n) + 0.5, 3),
            "vital_status": rng.choice([0, 1], n),
        }
    )


def _joined(cohort, mutants=("P001", "P002")):
    return pd.DataFrame(
        {
            "patient_id": cohort["patient_id"],
            "idh_mutated": [1 if p in mutants else 0 for p in cohort["patient_id"]],
        }
    )


def test_the_fold_wise_medians_are_each_training_fold_s(script):
    cohort = _cohort()
    out = script.compare(cohort, _joined(cohort), features=FEATURES, n_boot=20)

    # What the medians must be, recomputed from the folds themselves.
    rows = cohort[list(FEATURES) + ["overall_survival_months", "vital_status"]]
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    expected, imputed = [], []
    kps = rows["kps"].to_numpy(dtype=float)
    for train_idx, test_idx in skf.split(np.zeros(len(rows)), rows["vital_status"]):
        expected.append(float(np.nanmedian(kps[train_idx])))
        imputed.append(int(np.isnan(kps[test_idx]).sum()))

    full = out["full"]
    assert full["new"]["kps_median_by_training_fold"] == expected
    assert full["new"]["kps_imputed_by_validation_fold"] == imputed
    assert full["old"]["kps_median_by_training_fold"] is None
    assert sum(imputed) == full["new"]["kps_imputed"] == full["old"]["kps_imputed"]
    assert full["new"]["kps_median_whole_frame"] == float(np.nanmedian(kps))
    assert full["new"]["n_modelled"] == full["old"]["n_modelled"] == len(rows)


def test_with_nothing_to_impute_the_two_procedures_are_identical(script):
    cohort = _cohort(missing=0.0)
    out = script.compare(cohort, _joined(cohort), features=FEATURES, n_boot=20)
    for label in ("full", "without_idh"):
        old, new = out[label]["old"], out[label]["new"]
        assert out[label]["same_cindex"] is True
        assert old["cindex"] == new["cindex"]
        assert old["ci"] == new["ci"]
        assert old["fold_c"] == new["fold_c"]
        assert new["kps_median_by_training_fold"] == [None] * 5
        assert new["kps_imputed"] == old["kps_imputed"] == 0


def test_the_idh_mutant_patients_are_excluded_from_the_second_cohort(script):
    cohort = _cohort()
    out = script.compare(cohort, _joined(cohort), features=FEATURES, n_boot=10)
    assert out["full"]["n_patients"] == 90
    assert out["without_idh"]["n_patients"] == 88
    assert out["without_idh"]["n_idh_mutant_in_frame"] == 2
    for label in ("full", "without_idh"):
        assert out[label]["delta_cindex_new_minus_old"] == round(
            out[label]["new"]["cindex"] - out[label]["old"]["cindex"], 6
        )


def test_compare_calls_each_procedure_on_both_cohorts(script):
    seen = []

    def fake(tag):
        def fn(frame, features, **kwargs):
            seen.append((tag, len(frame), kwargs["n_boot"]))
            return {"cindex": 0.50 if tag == "old" else 0.52}

        return fn

    cohort = _cohort(n=40)
    out = script.compare(
        cohort,
        _joined(cohort),
        features=FEATURES,
        n_boot=7,
        old_fn=fake("old"),
        new_fn=fake("new"),
    )
    assert sorted(seen) == [("new", 38, 7), ("new", 40, 7), ("old", 38, 7), ("old", 40, 7)]
    assert out["full"]["delta_cindex_new_minus_old"] == 0.02
    assert out["full"]["same_cindex"] is False


def test_a_path_naming_the_other_cohort_is_refused(script, tmp_path):
    with pytest.raises(ValueError, match="TCGA only"):
        script.run("data/raw/cgga/clinical_cohort.parquet", "joined.parquet", tmp_path / "o.json")
    with pytest.raises(ValueError, match="TCGA only"):
        script.refuse_cgga(tmp_path / "cgga_external.parquet")


def test_run_writes_a_deterministic_record_and_refuses_to_overwrite_it(script, tmp_path):
    cohort = _cohort(n=60)
    cohort.to_parquet(tmp_path / "cohort.parquet")
    _joined(cohort).to_parquet(tmp_path / "joined.parquet")
    out = tmp_path / "rerun.json"
    paths = (tmp_path / "cohort.parquet", tmp_path / "joined.parquet", out)

    first = script.run(*paths, features=FEATURES, n_boot=10)
    assert out.is_file() and first["reads_cgga"] is False
    assert len(first["inputs"]["cohort"]["sha256"]) == 64
    with pytest.raises(FileExistsError, match="--overwrite"):
        script.run(*paths, features=FEATURES, n_boot=10)
    assert script.run(*paths, features=FEATURES, n_boot=10, overwrite=True) == first
