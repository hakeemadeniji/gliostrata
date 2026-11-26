"""
Pre-registered IDH sensitivity analysis, part 1: the Cox C-index.

WHY THIS EXISTS
    docs/mgmt_idh_calling.md section 6.3 pre-specified this analysis before
    the IDH calls existed. All 23 IDH-mutant patients stay in the primary
    cohort and the frozen 594-patient partition stands. Under WHO 2021 those
    cases would be classified as astrocytoma, IDH-mutant, so the primary
    estimand is re-reported here with them excluded. This is a robustness
    check on the primary analysis, not a change to it.

WHAT IS COMPARED
    The committed clinical-only baseline, unchanged:
    gliostrata.prediction.cox.fit_clinical_baseline with its default
    features (age, sex, kps, treatment_arm), 5-fold stratified
    cross-validation, and the pooled out-of-fold Harrell C-index. The same
    call is made twice:
        full cohort     every clinical-cohort patient (594)
        without IDH     the IDH-mutant patients removed
    Nothing about the model, the features, the folds or the seed differs
    between the two fits except which patients are in the cohort. It is the
    same estimator tests/test_cox.py::test_smoke_real_tcga_gbm_cohort
    exercises.

    IDH status comes from the joined table written by
    scripts/join_molecular_calls.py. Only idh_mutated == 1 is excluded. A
    patient with no IDH call is NOT excluded: an absent call is not
    evidence of anything, and this analysis removes known mutants only.

CONFIDENCE INTERVALS
    95% percentile bootstrap over patients (1000 resamples by default) of
    the pooled out-of-fold risk scores, exactly as the committed baseline
    computes its own interval. The model is NOT refitted inside each
    resample. Both intervals use the same seed but are computed separately
    on different numbers of patients, so they are not paired resamples.

DELTA
    delta = cindex_without_idh - cindex_full, signed: positive means the
    C-index is HIGHER without the IDH-mutant patients.
    delta_exceeds_threshold = abs(delta) > threshold, with the
    pre-specified threshold of 0.01. Strictly greater-than: a delta of
    exactly 0.01 does not exceed it. The comparison is made on the
    six-decimal value that is published in the JSON, so the flag is
    reproducible from that number alone.

WHAT THIS SCRIPT DOES NOT DO
    It does not refit with IDH as a covariate, does not alter the primary
    cohort, does not touch data/raw, and writes only --out.

Usage:
    python scripts/idh_sensitivity_cox.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from gliostrata.prediction.cox import DEFAULT_CLINICAL_FEATURES, fit_clinical_baseline

DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
DEFAULT_JOINED = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
DEFAULT_OUT = Path("data/results/tcga_gbm/idh_sensitivity_cox.json")

TIME_COL = "overall_survival_months"
EVENT_COL = "vital_status"

#: Pre-specified in docs/mgmt_idh_calling.md section 6.3.
DEFAULT_THRESHOLD = 0.01

#: Result keys, in order. docs/mgmt_idh_calling.md section 6.3 reads these.
RESULT_KEYS = (
    "n_full",
    "n_excluded",
    "cindex_full",
    "cindex_full_ci",
    "cindex_without_idh",
    "cindex_without_idh_ci",
    "delta",
    "delta_exceeds_threshold",
    "threshold",
)


def idh_mutant_ids(joined: pd.DataFrame) -> set[str]:
    """Patients with a non-null IDH-mutant call. A null call is never one."""
    for col in ("patient_id", "idh_mutated"):
        if col not in joined.columns:
            raise ValueError(f"joined: required column {col!r} is missing")
    ids = joined["patient_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError("joined: duplicate patient_id values")
    # eq(1) is <NA> where there is no call; fillna makes that an explicit
    # "not known to be mutant". A missing call never excludes a patient.
    called = joined["idh_mutated"].eq(1).fillna(False).to_numpy(dtype=bool)
    return set(ids[called])


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------
def fit_cohort(
    frame: pd.DataFrame,
    features=DEFAULT_CLINICAL_FEATURES,
    n_splits: int = 5,
    seed: int = 42,
    n_boot: int = 1000,
) -> dict:
    """Fit + cross-validate the committed clinical baseline on one cohort.

    Only the feature and outcome columns are passed, so both fits see an
    identical frame layout whatever else the cohort frame carries (the
    cohort's own all-null `idh_mutated` column, for instance, is never
    in play).
    """
    cols = list(features) + [TIME_COL, EVENT_COL]
    model, cv, report = fit_clinical_baseline(
        frame[cols], features=features, n_splits=n_splits, seed=seed, n_boot=n_boot
    )
    return {
        "cindex": float(cv.c_index),
        "ci": [float(cv.c_index_ci[0]), float(cv.c_index_ci[1])],
        "n_input": int(len(frame)),
        "n_used": int(report.n_used),
        "n_events": int(model.n_events_),
        "fold_c_indices": [float(c) for c in cv.fold_c_indices],
        "feature_names": list(report.feature_names),
    }


def analyze(
    cohort: pd.DataFrame,
    joined: pd.DataFrame,
    *,
    features=DEFAULT_CLINICAL_FEATURES,
    threshold: float = DEFAULT_THRESHOLD,
    n_boot: int = 1000,
    seed: int = 42,
    n_splits: int = 5,
    fitter=None,
) -> dict:
    """Fit the baseline on the full cohort and on the cohort without the
    IDH-mutant patients, and return the comparison.

    `fitter` is injectable so the threshold and exclusion logic can be
    tested without a model fit; production code leaves it as None.
    """
    fit = fit_cohort if fitter is None else fitter
    if "patient_id" not in cohort.columns:
        raise ValueError("cohort: required column 'patient_id' is missing")
    mutants = idh_mutant_ids(joined)
    ids = cohort["patient_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError("cohort: duplicate patient_id values")

    is_mutant = ids.isin(mutants).to_numpy()
    reduced = cohort.loc[~is_mutant]
    if len(reduced) < 2:
        raise ValueError("excluding the IDH-mutant patients leaves no cohort")

    full = fit(cohort, features, n_splits, seed, n_boot)
    without = fit(reduced, features, n_splits, seed, n_boot)

    # Compared on the six-decimal value that is published, so the flag is
    # exactly reproducible from the JSON and a delta of exactly `threshold`
    # does not exceed it. (Raw floats would not be: 0.710 - 0.700 is
    # 0.010000000000000009 in binary arithmetic.)
    delta = round(without["cindex"] - full["cindex"], 6)
    result = {
        "n_full": int(len(cohort)),
        "n_excluded": int(is_mutant.sum()),
        "cindex_full": round(full["cindex"], 6),
        "cindex_full_ci": [round(v, 6) for v in full["ci"]],
        "cindex_without_idh": round(without["cindex"], 6),
        "cindex_without_idh_ci": [round(v, 6) for v in without["ci"]],
        "delta": delta,
        "delta_exceeds_threshold": bool(abs(delta) > threshold),
        "threshold": float(threshold),
        # --- provenance: what was fitted, and on whom ---
        "n_without_idh": int(len(reduced)),
        "n_modelled_full": int(full["n_used"]),
        "n_modelled_without_idh": int(without["n_used"]),
        "n_events_full": int(full["n_events"]),
        "n_events_without_idh": int(without["n_events"]),
        "fold_c_indices_full": full["fold_c_indices"],
        "fold_c_indices_without_idh": without["fold_c_indices"],
        "features": list(features),
        "design_columns": list(full["feature_names"]),
        "n_splits": int(n_splits),
        "n_boot": int(n_boot),
        "seed": int(seed),
        "estimator": (
            "pooled out-of-fold Harrell C-index, 5-fold stratified CV "
            "(gliostrata.prediction.cox.fit_clinical_baseline, defaults)"
        ),
        "ci_method": (
            "95% percentile bootstrap of the pooled out-of-fold risk "
            "scores; the model is not refitted per resample"
        ),
        "delta_definition": "cindex_without_idh - cindex_full",
        "exclusion_rule": (
            "idh_mutated == 1 only; a patient with no IDH call is not excluded"
        ),
    }
    assert tuple(list(result)[: len(RESULT_KEYS)]) == RESULT_KEYS
    return result


def run(
    cohort_path: Path | str,
    joined_path: Path | str,
    out_path: Path | str,
    **kwargs,
) -> dict:
    """Execute the analysis, write the JSON result, and return it."""
    cohort = pd.read_parquet(cohort_path)
    joined = pd.read_parquet(joined_path, columns=["patient_id", "idh_mutated"])
    result = analyze(cohort, joined, **kwargs)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT),
                        help="clinical_cohort.parquet")
    parser.add_argument("--joined", default=str(DEFAULT_JOINED),
                        help="molecular_calls_joined.parquet")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output JSON")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                        help="absolute C-index change that is flagged (default 0.01)")
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-splits", type=int, default=5)
    args = parser.parse_args(argv)

    print("IDH sensitivity analysis (pre-specified) -- Cox C-index")
    print("  Committed clinical baseline, fitted twice: full cohort, and")
    print("  the cohort without the IDH-mutant patients. Robustness check")
    print("  only; the primary cohort is not changed.")
    print()

    r = run(
        args.cohort, args.joined, args.out,
        threshold=args.threshold, n_boot=args.n_boot,
        seed=args.seed, n_splits=args.n_splits,
    )

    def show(label, c, ci):
        print(f"{label:<34} C = {c:.4f}  95% CI [{ci[0]:.4f}, {ci[1]:.4f}]")

    print(f"full cohort                        n = {r['n_full']} "
          f"({r['n_modelled_full']} modelled, {r['n_events_full']} deaths)")
    show("  clinical baseline", r["cindex_full"], r["cindex_full_ci"])
    print(f"without the IDH-mutant patients    n = {r['n_without_idh']} "
          f"({r['n_modelled_without_idh']} modelled, "
          f"{r['n_events_without_idh']} deaths)")
    show("  clinical baseline", r["cindex_without_idh"], r["cindex_without_idh_ci"])
    print()
    print(f"IDH-mutant patients excluded      {r['n_excluded']}")
    print(f"delta (without - full)             {r['delta']:+.4f}")
    verdict = "EXCEEDS" if r["delta_exceeds_threshold"] else "does not exceed"
    print(f"|delta| {verdict} the pre-specified threshold of {r['threshold']:g}")
    print(f"\nWrote: {args.out}")


if __name__ == "__main__":
    main()
