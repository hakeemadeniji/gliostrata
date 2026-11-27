"""
T2: external validation of the Cox survival model on CGGA.

WHAT THIS IS
    The pre-registered T2 external validation (docs/preregistration.md,
    section 2): fit the Cox model on TCGA-GBM, evaluate it once on CGGA,
    and report the external C-index with a bootstrap 95% CI. The
    pre-registered criterion for the C-index row of T2 is an external
    C-index >= 0.65 with the 95% CI lower bound clearing 0.65. (T2 has two
    further rows -- improvement over a 2024 baseline and Brier score at 18
    months -- that this script does not evaluate.)

THE CONTAMINATION RULE
    CGGA is the external validation set. It is read for this single
    confirmatory evaluation and for nothing else. Nothing here tunes a
    hyperparameter, selects a feature, or tries a second configuration:
    the three variants below were fixed before the first CGGA evaluation,
    all three are reported, and none is picked after the fact.

THE THREE PRE-SPECIFIED VARIANTS
    1. Committed baseline: age, sex, kps, treatment_arm. Its internal
       C-index is the committed 0.6735 and is checked against it (+/-0.005).
       CGGA publishes no performance score, so KPS is absent for every
       patient and this variant CANNOT be evaluated externally. That is a
       finding, reported as such; KPS is never imputed.
    2. Transferable: age, sex, treatment_arm -- present in both cohorts.
       This is the pre-specified headline for the verdict.
    3. Molecular: age, sex, treatment_arm, idh_mutated, mgmt_methylated.
       Secondary. CGGA's IDH and MGMT status are its own clinical
       annotations, not the calling pipeline used for TCGA, so this
       comparison is cross-platform (docs/mgmt_idh_calling.md,
       docs/cgga_ingestion.md).

RULES FIXED BEFORE THE FIRST CGGA EVALUATION
    * Evaluation unit: one row per patient. CGGA reuses ids across
      batches; the first row in stored order is kept. Scoring all rows
      would count some patients twice.
    * Missing data: complete cases on the variant's features and on the
      outcome, in BOTH cohorts. Nothing is imputed anywhere.
    * A categorical value never seen in the TCGA fit set excludes that
      patient (counted), rather than being coded as the reference level.
    * Internal C-index: variant 1 is the committed fit_clinical_baseline
      unchanged; variants 2 and 3 use 5-fold stratified CV (pooled
      out-of-fold C-index) on the same complete-case TCGA rows that the
      external model is fitted on.
    * External CI: percentile bootstrap of the CGGA patients around the
      fixed TCGA model (1000 resamples, seed 42); the model is not
      refitted inside a resample.
    * Verdict is decided by variant 2. clears_t2_threshold is true only
      if the point estimate AND the CI lower bound are >= 0.65.

WHAT THIS IS NOT
    The pre-registration reports final T2 numbers on a frozen external
    partition (cgga_external_v1), which does not exist yet. This is the
    single evaluation on the ingested CGGA cohort, not a read of a sealed
    partition. Whether that changes T2's registered status is the
    project lead's decision; docs/preregistration.md is not edited.

Usage:
    python scripts/validate_cox_cgga.py --dry-run   # availability only, no fit, no metric
    python scripts/validate_cox_cgga.py             # the single confirmatory evaluation
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.prediction.cox import (
    CoxPHModel,
    build_design_matrix,
    concordance_index,
    cross_validate_cox,
    fit_clinical_baseline,
)

DEFAULT_TCGA = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
DEFAULT_JOINED = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
DEFAULT_CGGA = Path("data/raw/cgga/clinical_cohort.parquet")
DEFAULT_OUT = Path("data/results/cgga/cox_external_validation.json")

TIME_COL = "overall_survival_months"
EVENT_COL = "vital_status"
CALL_COLUMNS = ("idh_mutated", "mgmt_methylated")

T2_THRESHOLD = 0.65
COMMITTED_INTERNAL_CINDEX = 0.6735
INTERNAL_TOLERANCE = 0.005
N_BOOT = 1000
SEED = 42
N_SPLITS = 5
INTERNAL_N_BOOT = 20  # the internal CI is not reported; kept small

VARIANT_1 = ("age", "sex", "kps", "treatment_arm")
VARIANT_2 = ("age", "sex", "treatment_arm")
VARIANT_3 = ("age", "sex", "treatment_arm", "idh_mutated", "mgmt_methylated")

HARMONISATION_CAVEAT = (
    "CGGA IDH/MGMT calling method is not necessarily identical to TCGA's; "
    "see docs/mgmt_idh_calling.md and docs/cgga_ingestion.md"
)
VERDICT_MET = "met"
VERDICT_CI = "not met - CI lower bound below 0.65"
VERDICT_SHORT = "not met - best variant falls short"


# --------------------------------------------------------------------------
# Cohort preparation (no model, no metric)
# --------------------------------------------------------------------------
def load_tcga(cohort: pd.DataFrame, joined: pd.DataFrame) -> pd.DataFrame:
    """TCGA cohort with the joined IDH / MGMT calls in place of the
    cohort's own all-null placeholder columns."""
    for col in ("patient_id", *CALL_COLUMNS):
        if col not in joined.columns:
            raise ValueError(f"joined molecular calls: required column {col!r} is missing")
    calls = joined[["patient_id", *CALL_COLUMNS]].copy()
    frame = cohort.drop(columns=list(CALL_COLUMNS), errors="ignore")
    frame = frame.merge(calls, on="patient_id", how="left", validate="one_to_one")
    for col in CALL_COLUMNS:
        # A missing call stays missing; it is never a wildtype call.
        frame[col] = pd.to_numeric(frame[col]).astype("Float64")
    return frame


def one_row_per_patient(cgga: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Keep the first row per patient_id in stored order.

    CGGA reuses ids across batches. The rule is fixed in advance and does
    not look at any value: first row in stored order wins.
    """
    ids = cgga["patient_id"].astype(str)
    if ids.duplicated().any():
        dup = sorted(ids[ids.duplicated()].unique())
        print(f"  CGGA: {len(dup)} patient ids appear in more than one batch; "
              f"keeping the first row of each")
    keep = ~ids.duplicated(keep="first").to_numpy()
    return cgga.loc[keep].reset_index(drop=True), int((~keep).sum())


def feature_availability(frame: pd.DataFrame, features) -> dict[str, int]:
    """Non-null count per feature. A zero means the cohort does not have it."""
    return {f: int(frame[f].notna().sum()) if f in frame.columns else 0
            for f in features}


def complete_cases(frame: pd.DataFrame, features) -> tuple[pd.DataFrame, dict[str, int]]:
    """Rows with every feature and the outcome observed; nothing imputed."""
    needed = list(features) + [TIME_COL, EVENT_COL]
    missing_col = [c for c in needed if c not in frame.columns]
    if missing_col:
        raise ValueError(f"required column(s) missing: {missing_col}")
    have = frame[list(features)].notna().all(axis=1)
    outcome = frame[TIME_COL].notna() & frame[EVENT_COL].notna()
    return (
        frame.loc[(have & outcome).to_numpy()].reset_index(drop=True),
        {
            "n_input": int(len(frame)),
            "n_dropped_missing_feature": int((~have).sum()),
            "n_dropped_missing_outcome": int((have & ~outcome).sum()),
        },
    )


def _numeric(frame: pd.DataFrame, col: str) -> np.ndarray:
    values = pd.to_numeric(frame[col]).to_numpy(dtype="float64", na_value=np.nan)
    if np.isnan(values).any():
        raise ValueError(f"column {col!r} has missing values; they are never imputed")
    return values


def _categorical(frame: pd.DataFrame, col: str) -> pd.Series:
    return frame[col].astype("object").map(lambda v: str(v))


# --------------------------------------------------------------------------
# Fitting -- TCGA only. Nothing in this section accepts an external cohort.
# --------------------------------------------------------------------------
@dataclass
class Fitted:
    """A variant fitted on TCGA. Carries what is needed to encode another
    cohort exactly the way the training cohort was encoded."""

    features: tuple
    design_columns: list
    categorical_levels: dict
    model: CoxPHModel
    n_fit: int
    n_modelled: int
    n_events: int
    internal_cindex: float
    hazard_ratios: dict
    fit_patient_ids: frozenset


def fit_variant(tcga: pd.DataFrame, features, *, n_splits=N_SPLITS, seed=SEED,
                n_boot=INTERNAL_N_BOOT) -> Fitted:
    """Complete-case TCGA fit of one variant, plus its internal C-index."""
    features = tuple(features)
    fit_frame, _ = complete_cases(tcga, features)
    for feat in features:
        if feat not in ("sex", "treatment_arm") and feat in fit_frame.columns:
            if pd.api.types.is_numeric_dtype(fit_frame[feat]) or feat in CALL_COLUMNS:
                fit_frame[feat] = _numeric(fit_frame, feat)
    X, time, event, report = build_design_matrix(
        fit_frame, features, impute="drop"
    )
    if report.n_used != len(fit_frame):  # pragma: no cover - guard
        raise AssertionError("complete-case rows were dropped inside the design")
    cv = cross_validate_cox(X, time, event, n_splits=n_splits, seed=seed, n_boot=n_boot)
    model = CoxPHModel().fit(X, time, event, feature_names=report.feature_names)
    summary = model.summary().set_index("feature")
    return Fitted(
        features=features,
        design_columns=list(report.feature_names),
        categorical_levels={k: list(v) for k, v in report.categorical_levels.items()},
        model=model,
        n_fit=int(len(fit_frame)),
        n_modelled=int(report.n_used),
        n_events=int(model.n_events_),
        internal_cindex=float(cv.c_index),
        hazard_ratios={str(k): float(v) for k, v in summary["hazard_ratio"].items()},
        fit_patient_ids=frozenset(fit_frame["patient_id"].astype(str)),
    )


def committed_baseline_internal(tcga: pd.DataFrame, *, n_boot=N_BOOT,
                                seed=SEED) -> dict:
    """Variant 1: the committed fit_clinical_baseline, unchanged.

    This is the exact call tests/test_cox.py makes on the real cohort, so
    its C-index is the committed 0.6735 and nothing else.
    """
    model, cv, report = fit_clinical_baseline(tcga, n_boot=n_boot, seed=seed)
    summary = model.summary().set_index("feature")
    return {
        "n_fit": int(len(tcga)),
        "n_modelled": int(report.n_used),
        "n_events": int(model.n_events_),
        "internal_cindex": float(cv.c_index),
        "hazard_ratios": {str(k): float(v) for k, v in summary["hazard_ratio"].items()},
    }


# --------------------------------------------------------------------------
# External evaluation -- CGGA is scored by a model that is already fixed
# --------------------------------------------------------------------------
def design_external(frame: pd.DataFrame, fitted: Fitted) -> tuple[np.ndarray, np.ndarray]:
    """Encode another cohort exactly as the TCGA fit encoded its own.

    Returns (X, unseen) where `unseen` marks rows carrying a categorical
    value that never occurred in the TCGA fit set. Those rows cannot be
    scored under this model and are excluded by the caller, never coded
    as the reference level.
    """
    columns: dict[str, np.ndarray] = {}
    unseen = np.zeros(len(frame), dtype=bool)
    if len(frame) == 0:
        return np.zeros((0, len(fitted.design_columns))), unseen
    for feat in fitted.features:
        if feat in fitted.categorical_levels:
            raw = _categorical(frame, feat)
            levels = fitted.categorical_levels[feat]
            unseen |= ~raw.isin(set(levels)).to_numpy()
            for lvl in levels:
                name = f"{feat}_{lvl}"
                if name in fitted.design_columns:
                    columns[name] = (raw == lvl).to_numpy(dtype=float)
        elif feat in fitted.design_columns:
            columns[feat] = _numeric(frame, feat)
    if set(columns) != set(fitted.design_columns):  # pragma: no cover - guard
        raise AssertionError("external design does not match the training design")
    X = np.column_stack([columns[c] for c in fitted.design_columns])
    return X, unseen


def bootstrap_cindex(time, event, risk, *, n_boot=N_BOOT, seed=SEED):
    """Percentile bootstrap of the C-index over the evaluated patients.

    The risk scores are fixed (the TCGA model is not refitted); only the
    patient set is resampled.
    """
    rng = np.random.default_rng(seed)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    risk = np.asarray(risk, dtype=float)
    n = len(time)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            values.append(concordance_index(time[idx], event[idx], risk[idx]))
        except ValueError:  # a resample with no comparable pair
            continue
    if len(values) < 2:  # pragma: no cover - degenerate cohort
        return None, None, len(values)
    v = np.asarray(values, dtype=float)
    return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975)), len(values)


def clears_t2_threshold(cindex, ci_lower, threshold=T2_THRESHOLD) -> bool:
    """True only if the point estimate AND the CI lower bound are >= threshold."""
    if cindex is None or ci_lower is None:
        return False
    return bool(cindex >= threshold and ci_lower >= threshold)


def evaluate_external(fitted: Fitted, cgga: pd.DataFrame, *, n_boot=N_BOOT,
                      seed=SEED) -> dict:
    """Score one fitted variant on the external cohort, once."""
    availability = feature_availability(cgga, fitted.features)
    absent = [f for f, n in availability.items() if n == 0]
    out: dict = {"external_evaluable": not absent,
                 "cgga_feature_availability": availability}
    if absent:
        out.update(
            absent_in_cgga=absent,
            external_reason=(
                f"CGGA has zero {', '.join(absent)} records; this variant "
                f"cannot be externally validated"
            ),
            external_cindex=None,
            external_ci_lower=None,
            external_ci_upper=None,
            n_cgga_eval=0,
            clears_t2_threshold=False,
        )
        return out

    eval_frame, dropped = complete_cases(cgga, fitted.features)
    X, unseen = design_external(eval_frame, fitted)
    keep = ~unseen
    eval_frame = eval_frame.loc[keep].reset_index(drop=True)
    X = X[keep]
    out.update(
        n_cgga_excluded={
            **dropped,
            "n_dropped_category_unseen_in_tcga": int(unseen.sum()),
        },
        n_cgga_eval=int(len(eval_frame)),
    )
    if len(eval_frame) == 0:
        out.update(
            external_reason="no CGGA patient is evaluable under this variant",
            external_cindex=None,
            external_ci_lower=None,
            external_ci_upper=None,
            external_evaluable=False,
            clears_t2_threshold=False,
        )
        return out

    risk = fitted.model.predict_log_hazard(X)
    time = _numeric(eval_frame, TIME_COL)
    event = _numeric(eval_frame, EVENT_COL).astype(int)
    c = concordance_index(time, event, risk)
    lo, hi, n_valid = bootstrap_cindex(time, event, risk, n_boot=n_boot, seed=seed)
    out.update(
        n_cgga_deaths=int(event.sum()),
        external_cindex=c,
        external_ci_lower=lo,
        external_ci_upper=hi,
        n_bootstrap_valid=int(n_valid),
        clears_t2_threshold=clears_t2_threshold(c, lo),
    )
    return out


def t2_verdict(headline: dict) -> tuple[str, str]:
    """Verdict on T2's C-index criterion, decided by the headline variant."""
    c = headline.get("external_cindex")
    lo = headline.get("external_ci_lower")
    if c is None or lo is None:
        return (
            VERDICT_SHORT,
            "The transferable variant has no evaluable CGGA patients, so "
            "T2's C-index criterion cannot be met on this cohort.",
        )
    hi = headline["external_ci_upper"]
    n = headline["n_cgga_eval"]
    numbers = f"external C-index {c:.4f} (95% CI {lo:.4f} to {hi:.4f}, n = {n})"
    lead_in = "The transferable model (age, sex, treatment_arm) reaches an "
    if lo >= T2_THRESHOLD:
        return (
            VERDICT_MET,
            lead_in + numbers + "; the lower bound clears 0.65, so T2's "
            "C-index criterion is met on CGGA.",
        )
    if c >= T2_THRESHOLD:
        return (
            VERDICT_CI,
            lead_in + numbers + "; the point estimate is at or above 0.65 "
            "but the lower bound does not clear it.",
        )
    return (
        VERDICT_SHORT,
        lead_in + numbers + "; the point estimate is below 0.65. The "
        "committed baseline itself cannot be evaluated on CGGA because "
        "KPS is absent.",
    )


def _rounded(value, digits=6):
    return None if value is None else round(float(value), digits)


def _evaluable_variant(fitted: Fitted, ev: dict, *, n_boot, seed=SEED,
                       caveat=None) -> dict:
    """Result block for a variant that can be scored externally."""
    block = {
        "features": list(fitted.features),
        "n_tcga_fit": int(fitted.n_fit),
        "n_cgga_eval": int(ev["n_cgga_eval"]),
        "internal_cindex": _rounded(fitted.internal_cindex),
        "external_cindex": _rounded(ev["external_cindex"]),
        "external_ci_lower": _rounded(ev["external_ci_lower"]),
        "external_ci_upper": _rounded(ev["external_ci_upper"]),
        "clears_t2_threshold": bool(ev["clears_t2_threshold"]),
    }
    if caveat is not None:
        block["harmonisation_caveat"] = caveat
    block.update(
        external_evaluable=bool(ev["external_evaluable"]),
        n_tcga_events=int(fitted.n_events),
        n_cgga_deaths=int(ev.get("n_cgga_deaths", 0)),
        n_cgga_excluded=ev.get("n_cgga_excluded"),
        tcga_hazard_ratios={k: _rounded(v, 4) for k, v in fitted.hazard_ratios.items()},
        internal_cindex_method=(
            "5-fold stratified CV on the same complete-case TCGA rows, "
            "pooled out-of-fold Harrell C-index (internal, not a T2 number)"
        ),
        bootstrap={
            "method": "percentile bootstrap of the evaluated CGGA patients "
                      "around the fixed TCGA model; the model is not refitted",
            "n_resamples_requested": int(n_boot),
            "n_resamples_valid": int(ev.get("n_bootstrap_valid", 0)),
            "seed": int(seed),
            "level": 0.95,
        },
    )
    return block


def run(tcga_frame: pd.DataFrame, cgga_rows: pd.DataFrame, *, n_boot=N_BOOT,
        seed=SEED, n_splits=N_SPLITS, fit_fn=None, baseline_fn=None) -> dict:
    """The whole pre-specified evaluation: fit on TCGA, evaluate on CGGA
    once per variant. `fit_fn` / `baseline_fn` are injectable for tests."""
    fit_fn = fit_variant if fit_fn is None else fit_fn
    baseline_fn = committed_baseline_internal if baseline_fn is None else baseline_fn

    # CGGA enters here for schema handling only: one row per patient.
    cgga, n_duplicate_rows = one_row_per_patient(cgga_rows)
    cohort_info = {
        "n_rows_ingested": int(len(cgga_rows)),
        "n_unique_patients": int(len(cgga)),
        "n_duplicate_rows_dropped": n_duplicate_rows,
        "n_with_survival_and_vital_status": int(
            (cgga[TIME_COL].notna() & cgga[EVENT_COL].notna()).sum()
        ),
    }

    # --- variant 1: the committed baseline, internal only -------------
    base = baseline_fn(tcga_frame, n_boot=n_boot, seed=seed)
    v1_fit = Fitted(
        features=VARIANT_1, design_columns=[], categorical_levels={}, model=None,
        n_fit=base["n_fit"], n_modelled=base["n_modelled"],
        n_events=base["n_events"], internal_cindex=base["internal_cindex"],
        hazard_ratios=base["hazard_ratios"], fit_patient_ids=frozenset(),
    )
    v1_ev = evaluate_external(v1_fit, cgga)
    matches = abs(base["internal_cindex"] - COMMITTED_INTERNAL_CINDEX) <= INTERNAL_TOLERANCE
    v1 = {
        "features": list(VARIANT_1),
        "n_tcga_fit": int(base["n_fit"]),
        "internal_cindex": _rounded(base["internal_cindex"]),
        "internal_cindex_matches_committed": bool(matches),
        "external_evaluable": bool(v1_ev["external_evaluable"]),
    }
    if not v1_ev["external_evaluable"]:
        v1["external_reason"] = (
            "CGGA has zero KPS records; the committed baseline cannot be "
            "externally validated"
            if "kps" in v1_ev.get("absent_in_cgga", [])
            else v1_ev["external_reason"]
        )
    v1.update(
        n_tcga_modelled=int(base["n_modelled"]),
        n_tcga_events=int(base["n_events"]),
        committed_internal_cindex=COMMITTED_INTERNAL_CINDEX,
        internal_tolerance=INTERNAL_TOLERANCE,
        absent_in_cgga=v1_ev.get("absent_in_cgga", []),
        cgga_feature_availability=v1_ev["cgga_feature_availability"],
        tcga_hazard_ratios={k: _rounded(v, 4) for k, v in base["hazard_ratios"].items()},
    )

    # --- variants 2 and 3: fit on TCGA, then score CGGA once ----------
    v2_fit = fit_fn(tcga_frame, VARIANT_2, n_splits=n_splits, seed=seed)
    v2 = _evaluable_variant(
        v2_fit, evaluate_external(v2_fit, cgga, n_boot=n_boot, seed=seed),
        n_boot=n_boot, seed=seed,
    )
    v3_fit = fit_fn(tcga_frame, VARIANT_3, n_splits=n_splits, seed=seed)
    v3 = _evaluable_variant(
        v3_fit, evaluate_external(v3_fit, cgga, n_boot=n_boot, seed=seed),
        n_boot=n_boot, seed=seed, caveat=HARMONISATION_CAVEAT,
    )

    verdict, explanation = t2_verdict(v2)
    return {
        "variant_1_baseline_committed": v1,
        "variant_2_transferable": v2,
        "variant_3_molecular": v3,
        "t2_threshold": T2_THRESHOLD,
        "t2_verdict": verdict,
        "t2_verdict_explanation": explanation,
        # --- provenance ---
        "internal_sanity_ok": bool(matches),
        "headline_variant": "variant_2_transferable",
        "pre_specification": {
            "variants_fixed_before_first_cgga_evaluation": True,
            "evaluation_unit": "one row per patient (first row in stored order)",
            "missing_data": "complete cases on features and outcome in both "
                            "cohorts; nothing imputed",
            "unseen_category": "a categorical value absent from the TCGA fit "
                               "set excludes the patient (counted)",
            "verdict_rule": "decided by the headline variant; true only if the "
                            "point estimate and the CI lower bound are >= 0.65",
        },
        "cgga_cohort": cohort_info,
        "scope": [
            "Evaluates the C-index row of T2 only; the improvement-over-2024-"
            "baseline and Brier rows are not evaluated here.",
            "Not a read of a sealed cgga_external_v1 partition, which does not "
            "exist yet; this is the single evaluation on the ingested cohort.",
        ],
    }


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def load_inputs(tcga_path, joined_path, cgga_path):
    tcga = pd.read_parquet(tcga_path)
    joined = pd.read_parquet(joined_path, columns=["patient_id", *CALL_COLUMNS])
    return load_tcga(tcga, joined), pd.read_parquet(cgga_path)


VARIANTS = (
    ("variant_1_baseline_committed", VARIANT_1),
    ("variant_2_transferable", VARIANT_2),
    ("variant_3_molecular", VARIANT_3),
)
CATEGORICAL = ("sex", "treatment_arm")


def dry_run(tcga_frame: pd.DataFrame, cgga_rows: pd.DataFrame) -> dict:
    """Availability and evaluable counts only: no fit, no metric."""
    cgga, n_dup = one_row_per_patient(cgga_rows)
    report = {
        "n_rows_ingested": int(len(cgga_rows)),
        "n_unique_patients": int(len(cgga)),
        "n_duplicate_rows_dropped": n_dup,
        "variants": {},
    }
    for name, feats in VARIANTS:
        avail = feature_availability(cgga, feats)
        absent = [f for f, n in avail.items() if n == 0]
        entry = {"feature_availability": avail, "absent_in_cgga": absent,
                 "externally_evaluable": not absent}
        if not absent:
            tcga_cc, _ = complete_cases(tcga_frame, feats)
            cc, dropped = complete_cases(cgga, feats)
            unseen = 0
            for feat in feats:
                if feat in CATEGORICAL:
                    levels = set(tcga_cc[feat].astype(str))
                    unseen += int((~cc[feat].astype(str).isin(levels)).sum())
            entry.update(
                n_tcga_fit=int(len(tcga_cc)),
                n_cgga_evaluable=int(len(cc) - unseen),
                n_cgga_excluded={**dropped, "n_category_unseen_in_tcga": unseen},
            )
        report["variants"][name] = entry
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--tcga", default=str(DEFAULT_TCGA),
                        help="TCGA clinical_cohort.parquet")
    parser.add_argument("--joined", default=str(DEFAULT_JOINED),
                        help="molecular_calls_joined.parquet")
    parser.add_argument("--cgga", default=str(DEFAULT_CGGA),
                        help="CGGA clinical_cohort.parquet (read once)")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output JSON")
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-splits", type=int, default=N_SPLITS)
    parser.add_argument("--dry-run", action="store_true",
                        help="availability and evaluable counts only; fits "
                             "nothing and computes no metric")
    args = parser.parse_args(argv)

    print("T2 -- external validation of the Cox model on CGGA")
    print("  Three pre-specified variants, fit on TCGA, scored on CGGA once.")
    print("  CGGA is not used for any tuning, feature selection or choice.")
    print()
    tcga_frame, cgga_rows = load_inputs(args.tcga, args.joined, args.cgga)

    if args.dry_run:
        report = dry_run(tcga_frame, cgga_rows)
        print(f"CGGA rows {report['n_rows_ingested']}, unique patients "
              f"{report['n_unique_patients']} "
              f"({report['n_duplicate_rows_dropped']} duplicate rows dropped)")
        for name, entry in report["variants"].items():
            print(f"  {name}")
            print(f"      feature availability: {entry['feature_availability']}")
            print(f"      absent in CGGA: {entry['absent_in_cgga']}")
            if entry["externally_evaluable"]:
                print(f"      TCGA fit n = {entry['n_tcga_fit']}, "
                      f"CGGA evaluable n = {entry['n_cgga_evaluable']}")
        print("\nDry run: nothing fitted, no metric computed.")
        return

    result = run(tcga_frame, cgga_rows, n_boot=args.n_boot, seed=args.seed,
                 n_splits=args.n_splits)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")

    v1, v2, v3 = (result["variant_1_baseline_committed"],
                  result["variant_2_transferable"], result["variant_3_molecular"])
    cohort = result["cgga_cohort"]
    print(f"CGGA: {cohort['n_rows_ingested']} rows, "
          f"{cohort['n_unique_patients']} unique patients "
          f"({cohort['n_duplicate_rows_dropped']} duplicate rows dropped)\n")
    print("1. Committed baseline (age, sex, kps, treatment_arm)")
    print(f"   internal C-index {v1['internal_cindex']:.4f} "
          f"(committed {COMMITTED_INTERNAL_CINDEX}; matches: "
          f"{v1['internal_cindex_matches_committed']})")
    print(f"   external: NOT EVALUABLE -- {v1['external_reason']}\n")

    def show(label, block):
        print(label)
        print(f"   TCGA fit n = {block['n_tcga_fit']}, "
              f"internal C-index {block['internal_cindex']:.4f}")
        if block["external_cindex"] is None:
            print(f"   external: NOT EVALUABLE -- "
                  f"{block.get('external_reason')}\n")
            return
        print(f"   CGGA n = {block['n_cgga_eval']} "
              f"({block['n_cgga_deaths']} deaths)")
        print(f"   external C-index {block['external_cindex']:.4f}  "
              f"95% CI [{block['external_ci_lower']:.4f}, "
              f"{block['external_ci_upper']:.4f}]")
        print(f"   CI lower bound clears {T2_THRESHOLD}: "
              f"{block['clears_t2_threshold']}\n")

    show("2. Transferable (age, sex, treatment_arm) -- the headline", v2)
    show("3. Molecular (age, sex, treatment_arm, IDH, MGMT) -- secondary, "
         "cross-platform", v3)
    print(f"T2 (C-index criterion): {result['t2_verdict']}")
    print(f"  {result['t2_verdict_explanation']}")
    print(f"\nWrote: {args.out}")

    if not result["internal_sanity_ok"]:
        print("\nSANITY CHECK FAILED: variant 1's internal C-index is not within "
              f"+/-{INTERNAL_TOLERANCE} of {COMMITTED_INTERNAL_CINDEX}. The Cox "
              "configuration has drifted; the T2 numbers above are not valid.")
        raise SystemExit(2)


if __name__ == "__main__":
    main()
