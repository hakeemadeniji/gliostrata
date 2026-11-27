"""
Rerun the internal cross-validation of the committed clinical baseline with the
imputation done inside each fold, next to the earlier whole-cohort version.

WHY THIS EXISTS
    gliostrata.prediction.cox promised the TRAINING-fold median for a missing
    numeric feature, but until 2026-11-27 the median was computed on the whole
    cohort before the folds were drawn. `fit_clinical_baseline` now fits the
    design of each fold on that fold's training rows. This script puts the two
    procedures side by side on TCGA-GBM, for the committed baseline (variant 1
    of T2) and for its IDH sensitivity analysis.

WHAT IT DOES
    TCGA only: the clinical cohort and the joined IDH calls. No CGGA file is
    read, and a path containing "cgga" is refused before anything is opened.
    "old" is the earlier procedure, reproduced by the two calls that
    implemented it: `build_design_matrix` on the whole frame, then
    `cross_validate_cox`. "new" is `cross_validate_cohort_cox`, which
    `fit_clinical_baseline` now uses. Same patients, same folds (5-fold,
    stratified on event status, seed 42), same pooled out-of-fold C-index.
    The values already on record are read from the IDH sensitivity result.

WHAT IT DOES NOT DO
    It changes no result file and writes only --out. Variant 1 is never
    evaluated on CGGA, so no external number depends on either procedure.

Usage:
    python scripts/rerun_tcga_internal_cv.py
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

from gliostrata.prediction.cox import (
    DEFAULT_CLINICAL_FEATURES,
    build_design_matrix,
    cross_validate_cohort_cox,
    cross_validate_cox,
)

DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
DEFAULT_JOINED = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
DEFAULT_OUT = Path("data/results/tcga_gbm/foldwise_imputation_rerun.json")
#: The values on record: cindex_full and cindex_without_idh of this file.
RECORDED = Path("data/results/tcga_gbm/idh_sensitivity_cox.json")

N_SPLITS = 5
SEED = 42
N_BOOT = 1000

OLD_IMPUTATION = (
    "median of kps over all modelled patients of the frame, computed before the "
    "folds were drawn (the procedure before 2026-11-27)"
)
NEW_IMPUTATION = (
    "median of kps over the training rows of each cross-validation fold, applied "
    "to that fold's validation rows"
)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def refuse_cgga(path) -> None:
    """This script works on TCGA only: refuse a path that names the other cohort."""
    if "cgga" in str(path).lower():
        raise ValueError(
            f"{path}: this script reads TCGA only and opens no path containing 'cgga'"
        )


def load_idh_script():
    """Load scripts/idh_sensitivity_cox.py by path, for its exclusion rule."""
    path = Path(__file__).resolve().parent / "idh_sensitivity_cox.py"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_idh_sensitivity_cox_script", path
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _frame_columns(frame: pd.DataFrame, features) -> pd.DataFrame:
    """Only the feature and outcome columns: both fits see one frame layout."""
    cols = list(features) + ["overall_survival_months", "vital_status"]
    return frame[cols]


def _block(cv, report, *, imputation, medians_by_fold, imputed_by_fold) -> dict:
    return {
        "cindex": round(float(cv.c_index), 6),
        "ci": [round(float(v), 6) for v in cv.c_index_ci],
        "fold_c": [round(float(c), 6) for c in cv.fold_c_indices],
        "brier_18mo": round(float(cv.brier_score), 6),
        "n_input": int(report.n_input),
        "n_modelled": int(report.n_used),
        "n_events": int(cv.n_events),
        "kps_imputed": int(report.n_imputed.get("kps", 0)),
        "kps_median_whole_frame": report.imputed_values.get("kps"),
        "kps_median_by_training_fold": medians_by_fold,
        "kps_imputed_by_validation_fold": imputed_by_fold,
        "imputation": imputation,
    }


def whole_frame_cv(frame, features, *, n_splits=N_SPLITS, seed=SEED, n_boot=N_BOOT) -> dict:
    """The earlier procedure: one design for the whole frame, then CV on it."""
    x, time, event, report = build_design_matrix(
        _frame_columns(frame, features), features, impute="median_indicator"
    )
    cv = cross_validate_cox(x, time, event, n_splits=n_splits, seed=seed, n_boot=n_boot)
    return _block(
        cv, report, imputation=OLD_IMPUTATION, medians_by_fold=None, imputed_by_fold=None
    )


def foldwise_cv(frame, features, *, n_splits=N_SPLITS, seed=SEED, n_boot=N_BOOT) -> dict:
    """The corrected procedure: each fold's design is fitted on its training rows."""
    cols = _frame_columns(frame, features)
    cv, reports = cross_validate_cohort_cox(
        cols, features, impute="median_indicator", n_splits=n_splits, seed=seed, n_boot=n_boot
    )
    report = build_design_matrix(cols, features, impute="median_indicator")[3]
    return _block(
        cv,
        report,
        imputation=NEW_IMPUTATION,
        medians_by_fold=[r.imputed_values.get("kps") for r in reports],
        imputed_by_fold=[int(r.n_imputed.get("kps", 0)) for r in reports],
    )


# --------------------------------------------------------------------------
# The comparison
# --------------------------------------------------------------------------
def compare(
    cohort: pd.DataFrame,
    joined: pd.DataFrame,
    *,
    features=DEFAULT_CLINICAL_FEATURES,
    idh_script=None,
    n_splits=N_SPLITS,
    seed=SEED,
    n_boot=N_BOOT,
    old_fn=whole_frame_cv,
    new_fn=foldwise_cv,
) -> dict:
    """Old and new procedure, side by side, on the full cohort and without the
    IDH-mutant patients. `old_fn` / `new_fn` are injectable for tests."""
    idh_script = load_idh_script() if idh_script is None else idh_script
    if "patient_id" not in cohort.columns:
        raise ValueError("cohort: required column 'patient_id' is missing")
    ids = cohort["patient_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError("cohort: duplicate patient_id values")
    is_mutant = ids.isin(idh_script.idh_mutant_ids(joined)).to_numpy()
    if (~is_mutant).sum() < 2:
        raise ValueError("excluding the IDH-mutant patients leaves no cohort")

    out: dict = {}
    for label, frame in (("full", cohort), ("without_idh", cohort.loc[~is_mutant])):
        old = old_fn(frame, features, n_splits=n_splits, seed=seed, n_boot=n_boot)
        new = new_fn(frame, features, n_splits=n_splits, seed=seed, n_boot=n_boot)
        out[label] = {
            "n_patients": int(len(frame)),
            "n_idh_mutant_in_frame": int(is_mutant.sum()),
            "old": old,
            "new": new,
            "delta_cindex_new_minus_old": round(new["cindex"] - old["cindex"], 6),
            "same_cindex": bool(new["cindex"] == old["cindex"]),
        }
    return out


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(cohort_path, joined_path, out_path, *, overwrite=False, **kwargs) -> dict:
    """Read the two TCGA tables, compare, write the JSON, and return it."""
    for path in (cohort_path, joined_path, out_path, RECORDED):
        refuse_cgga(path)
    out_path = Path(out_path)
    if out_path.exists() and not overwrite:
        raise FileExistsError(
            f"{out_path} exists. This is a rerun record, not a routine output: pass "
            "--overwrite only for a deliberate, reported rerun."
        )
    cohort = pd.read_parquet(cohort_path)
    joined = pd.read_parquet(joined_path, columns=["patient_id", "idh_mutated"])
    result = compare(cohort, joined, **kwargs)

    recorded = None
    if Path(RECORDED).is_file():
        rec = json.loads(Path(RECORDED).read_text(encoding="utf-8"))
        recorded = {
            "source": RECORDED.as_posix(),
            "cindex_full": rec.get("cindex_full"),
            "cindex_without_idh": rec.get("cindex_without_idh"),
        }
    result = {
        "what": (
            "Internal C-index of the committed clinical baseline on TCGA-GBM: the "
            "whole-frame imputation procedure against fold-wise imputation"
        ),
        "features": list(kwargs.get("features", DEFAULT_CLINICAL_FEATURES)),
        "n_splits": kwargs.get("n_splits", N_SPLITS),
        "seed": kwargs.get("seed", SEED),
        "n_boot": kwargs.get("n_boot", N_BOOT),
        "inputs": {
            "cohort": {"path": Path(cohort_path).as_posix(), "sha256": _sha256(cohort_path)},
            "joined": {"path": Path(joined_path).as_posix(), "sha256": _sha256(joined_path)},
        },
        "recorded": recorded,
        "reads_cgga": False,
        **result,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT))
    parser.add_argument("--joined", default=str(DEFAULT_JOINED))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-splits", type=int, default=N_SPLITS)
    parser.add_argument("--overwrite", action="store_true",
                        help="replace an existing rerun record (a deliberate, reported rerun)")
    args = parser.parse_args(argv)

    print("Internal C-index of the committed baseline (age, sex, kps, treatment_arm)")
    print("TCGA-GBM only. 5-fold stratified CV, pooled out-of-fold. Old against new.\n")
    r = run(args.cohort, args.joined, args.out, overwrite=args.overwrite,
            n_boot=args.n_boot, seed=args.seed, n_splits=args.n_splits)
    recorded = r["recorded"] or {}
    print(f"{'cohort':<22}{'n':>5}{'modelled':>10}   {'old':<10}{'new':<10}{'new-old':<11}"
          f"{'recorded':<10}")
    for label, rec_key in (("full", "cindex_full"), ("without_idh", "cindex_without_idh")):
        b = r[label]
        o, n = b["old"], b["new"]
        print(f"{label:<22}{b['n_patients']:>5}{o['n_modelled']:>10}   "
              f"{o['cindex']:<10.6f}{n['cindex']:<10.6f}"
              f"{b['delta_cindex_new_minus_old']:<+11.6f}"
              f"{str(recorded.get(rec_key, '--')):<10}")
        print(f"{'':<37}old 95% CI {o['ci'][0]:.4f} to {o['ci'][1]:.4f}; "
              f"new 95% CI {n['ci'][0]:.4f} to {n['ci'][1]:.4f}")
        print(f"{'':<37}kps median by training fold: {n['kps_median_by_training_fold']}; "
              f"kps imputed by validation fold: {n['kps_imputed_by_validation_fold']}")
    print(f"\nWrote: {args.out}")


if __name__ == "__main__":
    main()

