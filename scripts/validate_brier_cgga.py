"""
T2: Brier score at 18 months on CGGA -- the second registered row of T2.

WHAT THIS IS
    docs/preregistration.md section 2 registers a Brier score at 18 months
    <= 0.25, decided by the upper bound of the 95% CI. The C-index row was
    evaluated on 2026-11-26 (scripts/validate_cox_cgga.py). This script
    evaluates the Brier row, once, on the same ingested CGGA cohort. It is the
    second read of CGGA and the result file says so.

WHAT IT IS NOT
    Not a re-run of the C-index row: that number is not recomputed here. The
    first run's result file is read only to confirm that the TCGA models are
    the same ones. Not a read of a sealed cgga_external_v1 partition, which
    does not exist yet.

RULES FIXED BEFORE THE FIRST CGGA BRIER EVALUATION
    * Metric: inverse-probability-of-censoring-weighted Brier score (Graf et
      al. 1999), gliostrata.prediction.cox.brier_score_ipcw, with
      Kaplan-Meier censoring weights estimated on the evaluated CGGA patients.
    * Patients: identical to the C-index row. One row per patient, complete
      cases, nothing imputed, an unseen category excludes the patient. Those
      rules are imported from scripts/validate_cox_cgga.py, not copied.
    * Models: the same three variants, fitted on TCGA by the same
      fit_variant. Variant 2 (age, sex, treatment_arm) decides the verdict.
    * CI: percentile bootstrap over the evaluated CGGA patients, 1000
      resamples, seed 42. Model and predictions fixed; censoring weights
      re-estimated inside each resample.
    * Verdict: computed by the pre-registered registry
      (gliostrata.evaluation.thresholds.evaluate, key T2_brier_18mo). PASS
      only if the point estimate AND the CI upper bound are <= 0.25.
    * Descriptive only: Brier at 12 and 24 months (section 3.4 names them)
      and the Brier score of a null model that predicts the cohort's own
      Kaplan-Meier survival for everyone.

Usage:
    python scripts/validate_brier_cgga.py --dry-run   # availability only
    python scripts/validate_brier_cgga.py             # the single evaluation
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.evaluation.thresholds import evaluate
from gliostrata.prediction.cox import brier_score_ipcw

DEFAULT_OUT = Path("data/results/cgga/brier18_external.json")
DEFAULT_CINDEX_RESULT = Path("data/results/cgga/cox_external_validation.json")

METRIC_KEY = "T2_brier_18mo"
HORIZON = 18.0
DESCRIPTIVE_HORIZONS = (12.0, 24.0)
N_BOOT = 1000
SEED = 42
N_SPLITS = 5
INTERNAL_TOLERANCE = 1e-6  # same code, same seed: the refit must be identical
CGGA_READ_NUMBER = 2

_COX_VALIDATOR = Path(__file__).resolve().parent / "validate_cox_cgga.py"


def load_cox_validator():
    """The C-index script, loaded by path (scripts/ is not a package).

    Its cohort rules and its fitter are reused unchanged so that the Brier row
    is scored on exactly the patients the C-index row was scored on.
    """
    spec = importlib.util.spec_from_file_location(
        "gliostrata_validate_cox_cgga_for_brier", _COX_VALIDATOR
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


CX = load_cox_validator()
TIME_COL, EVENT_COL = CX.TIME_COL, CX.EVENT_COL

# --------------------------------------------------------------------------
# Scoring set: the C-index script's rule, called not copied
# --------------------------------------------------------------------------
def scoring_set(fitted, cgga: pd.DataFrame):
    """Patients a variant is scored on, and their design matrix.

    Same sequence as CX.evaluate_external: complete cases, then exclusion of
    a categorical value unseen in the TCGA fit set.
    """
    frame, dropped = CX.complete_cases(cgga, fitted.features)
    X, unseen = CX.design_external(frame, fitted)
    keep = ~unseen
    return (
        frame.loc[keep].reset_index(drop=True),
        X[keep],
        {**dropped, "n_dropped_category_unseen_in_tcga": int(unseen.sum())},
    )


# --------------------------------------------------------------------------
# The metric
# --------------------------------------------------------------------------
def survival_at(fitted, X: np.ndarray, horizon: float) -> np.ndarray:
    """The fixed model's predicted S(horizon | x) for each patient."""
    return np.asarray(
        fitted.model.predict_survival_function(X, horizon), dtype=float
    ).ravel()


def brier_with_ci(time, event, survival_prob, horizon, *, n_boot=N_BOOT, seed=SEED):
    """IPCW Brier score at `horizon` with a percentile bootstrap CI.

    `survival_prob` is held fixed (the TCGA model is not refitted); only the
    patient set is resampled, and the censoring weights are re-estimated on
    each resample because they are computed inside brier_score_ipcw.
    """
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    surv = np.asarray(survival_prob, dtype=float).ravel()
    point = brier_score_ipcw(time, event, surv, horizon)
    rng = np.random.default_rng(seed)
    n = len(time)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        value = brier_score_ipcw(time[idx], event[idx], surv[idx], horizon)
        if np.isfinite(value):
            values.append(value)
    if len(values) < 2:  # pragma: no cover - degenerate cohort
        return {"score": float(point), "ci_lower": None, "ci_upper": None, "n_valid": len(values)}
    v = np.asarray(values, dtype=float)
    return {"score": float(point), "ci_lower": float(np.quantile(v, 0.025)),
            "ci_upper": float(np.quantile(v, 0.975)), "n_valid": len(values)}


def km_event_survival(time, event) -> tuple[np.ndarray, np.ndarray]:
    """Kaplan-Meier S(t) of the outcome distribution: (times, survival)."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event).astype(int)
    order = np.argsort(time, kind="mergesort")
    t_s, e_s = time[order], event[order]
    n = len(t_s)
    times, surv = [], []
    s = 1.0
    k = 0
    while k < n:
        j = k
        while j + 1 < n and t_s[j + 1] == t_s[k]:
            j += 1
        d = int(e_s[k:j + 1].sum())
        if d > 0:
            s *= 1.0 - d / (n - k)
            times.append(t_s[k])
            surv.append(s)
        k = j + 1
    return np.asarray(times, dtype=float), np.asarray(surv, dtype=float)


def km_survival_at(time, event, horizon: float) -> float:
    """Observed Kaplan-Meier survival at `horizon` (1.0 if nobody has died)."""
    times, surv = km_event_survival(time, event)
    if len(times) == 0 or horizon < times[0]:
        return 1.0
    idx = int(np.searchsorted(times, horizon, side="right")) - 1
    return float(surv[idx])


def null_model_brier(time, event, horizon: float) -> float:
    """Brier score of a model that predicts the cohort's own marginal
    Kaplan-Meier survival for every patient. Descriptive reference only."""
    time = np.asarray(time, dtype=float)
    s_null = km_survival_at(time, event, horizon)
    return float(brier_score_ipcw(time, event, np.full(len(time), s_null), horizon))


def registered_verdict(point, ci_lower, ci_upper) -> dict:
    """The verdict, computed by the pre-registered registry for METRIC_KEY."""
    v = evaluate(METRIC_KEY, point, ci_lower, ci_upper)
    out = v.to_dict()
    out["passed"] = bool(v.passed)
    return out


# --------------------------------------------------------------------------
# One variant, scored once
# --------------------------------------------------------------------------
def _r(value, digits=6):
    return None if value is None else round(float(value), digits)


def evaluate_variant(fitted, cgga: pd.DataFrame, *, n_boot=N_BOOT, seed=SEED):
    """Score one TCGA-fitted variant on CGGA. Returns (block, point, lo, hi)."""
    eval_frame, X, excluded = scoring_set(fitted, cgga)
    block = {
        "features": list(fitted.features),
        "n_tcga_fit": int(fitted.n_fit),
        "n_tcga_events": int(fitted.n_events),
        "internal_cindex": _r(fitted.internal_cindex),
        "n_cgga_eval": int(len(eval_frame)),
        "n_cgga_excluded": excluded,
        "external_evaluable": True,
    }
    if len(eval_frame) == 0:
        block.update(external_evaluable=False,
                     external_reason="no CGGA patient is evaluable under this variant")
        return block, None, None, None

    time = CX._numeric(eval_frame, TIME_COL)
    event = CX._numeric(eval_frame, EVENT_COL).astype(int)
    main = brier_with_ci(time, event, survival_at(fitted, X, HORIZON), HORIZON,
                         n_boot=n_boot, seed=seed)
    point, lo, hi = main["score"], main["ci_lower"], main["ci_upper"]
    if lo is not None and not (lo <= point <= hi):
        # A percentile interval that excludes its own point estimate: the
        # registry refuses it, so widen it to contain the point and say so.
        block["ci_widened_to_contain_point"] = True
        lo, hi = min(lo, point), max(hi, point)
    block.update(
        n_cgga_deaths=int(event.sum()),
        n_cgga_at_risk_at_horizon=int((time > HORIZON).sum()),
        observed_km_survival_at_horizon=_r(km_survival_at(time, event, HORIZON)),
        brier_18mo=_r(point), ci_lower=_r(lo), ci_upper=_r(hi),
        n_bootstrap_valid=int(main["n_valid"]),
        null_model_brier_18mo=_r(null_model_brier(time, event, HORIZON)),
    )
    other = {}
    for h in DESCRIPTIVE_HORIZONS:
        d = brier_with_ci(time, event, survival_at(fitted, X, h), h, n_boot=n_boot, seed=seed)
        other[f"brier_{int(h)}mo"] = {
            "score": _r(d["score"]), "ci_lower": _r(d["ci_lower"]), "ci_upper": _r(d["ci_upper"]),
            "null_model_brier": _r(null_model_brier(time, event, h)),
        }
    block["descriptive_other_horizons"] = other
    return block, point, lo, hi


def compare_with_cindex_run(result_path, blocks: dict) -> dict:
    """Confirm the TCGA fits are the ones the C-index row was computed with."""
    out = {"cindex_result_file": str(result_path), "checked": False}
    path = Path(result_path)
    if not path.is_file():
        out["reason"] = "C-index result file not found; fit identity not checked"
        return out
    try:
        first = json.loads(path.read_text(encoding="ascii"))
    except (OSError, ValueError) as exc:
        out["reason"] = f"C-index result file unreadable ({exc}); fit identity not checked"
        return out
    out["checked"] = True
    ok = True
    for key in ("variant_2_transferable", "variant_3_molecular"):
        a, b = first.get(key), blocks.get(key)
        if not a or not b:
            continue
        diff = abs(float(a["internal_cindex"]) - float(b["internal_cindex"]))
        same_n = a.get("n_cgga_eval") == b.get("n_cgga_eval")
        out[key] = {
            "internal_cindex_first_run": a["internal_cindex"],
            "internal_cindex_this_run": b["internal_cindex"],
            "abs_difference": diff,
            "n_cgga_eval_first_run": a.get("n_cgga_eval"),
            "n_cgga_eval_this_run": b.get("n_cgga_eval"),
        }
        ok = ok and same_n and diff <= INTERNAL_TOLERANCE
    out["matches_cindex_run"] = bool(ok)
    return out


def run(tcga_frame, cgga_rows: pd.DataFrame, *, n_boot=N_BOOT, seed=SEED,
        n_splits=N_SPLITS, fit_fn=None, cindex_result=DEFAULT_CINDEX_RESULT) -> dict:
    """The pre-specified evaluation: fit on TCGA, score CGGA once per variant.
    `fit_fn` is injectable for tests."""
    fit_fn = CX.fit_variant if fit_fn is None else fit_fn
    cgga, n_dup = CX.one_row_per_patient(cgga_rows)

    avail = CX.feature_availability(cgga, CX.VARIANT_1)
    absent = [f for f, n in avail.items() if n == 0]
    v1 = {"features": list(CX.VARIANT_1), "external_evaluable": not absent,
          "absent_in_cgga": absent, "cgga_feature_availability": avail}
    if absent:
        v1["external_reason"] = (f"CGGA has zero {', '.join(absent)} records; the "
                                 "committed baseline cannot be externally validated")

    v2, point, lo, hi = evaluate_variant(
        fit_fn(tcga_frame, CX.VARIANT_2, n_splits=n_splits, seed=seed),
        cgga, n_boot=n_boot, seed=seed)
    v3, _, _, _ = evaluate_variant(
        fit_fn(tcga_frame, CX.VARIANT_3, n_splits=n_splits, seed=seed),
        cgga, n_boot=n_boot, seed=seed)
    v3["harmonisation_caveat"] = CX.HARMONISATION_CAVEAT

    if v2["external_evaluable"]:
        verdict = registered_verdict(point, lo, hi)
    else:
        verdict = {"status": "NOT EVALUABLE", "passed": False, "reason": v2["external_reason"]}

    return {
        "metric_key": METRIC_KEY,
        "horizon_months": HORIZON,
        "headline_variant": "variant_2_transferable",
        "variant_1_baseline_committed": v1,
        "variant_2_transferable": v2,
        "variant_3_molecular": v3,
        "verdict": verdict,
        "fit_identity_check": compare_with_cindex_run(
            cindex_result, {"variant_2_transferable": v2, "variant_3_molecular": v3}),
        "cgga_cohort": {
            "n_rows_ingested": int(len(cgga_rows)),
            "n_unique_patients": int(len(cgga)),
            "n_duplicate_rows_dropped": n_dup,
        },
        "cgga_read": {
            "read_number": CGGA_READ_NUMBER,
            "previous_reads": "1: T2 C-index row, 2026-11-26 (cox_external_validation.json)",
        },
        "method": {
            "score": "IPCW Brier score (Graf et al. 1999); Kaplan-Meier censoring "
                     "weights on the evaluated CGGA patients",
            "prediction": "TCGA-fitted Cox model, survival probability at the horizon",
            "bootstrap": "percentile bootstrap over evaluated CGGA patients; model and "
                         "predictions fixed, censoring weights re-estimated per resample",
            "n_resamples_requested": int(n_boot), "seed": int(seed), "level": 0.95,
            "verdict_rule": "gliostrata.evaluation.thresholds.evaluate on the registered "
                            "key; PASS only if the point estimate and the CI upper bound "
                            "are both <= 0.25",
        },
        "pre_specification": {
            "rules_fixed_before_first_cgga_brier_evaluation": True,
            "patients": "identical to the C-index row (rules imported from "
                        "scripts/validate_cox_cgga.py)",
            "descriptive_only": ["brier_12mo", "brier_24mo", "null_model_brier"],
        },
        "scope": [
            "Evaluates the Brier-at-18-months row of T2 only.",
            "Not a read of a sealed cgga_external_v1 partition, which does not exist "
            "yet; this is the evaluation on the ingested cohort.",
        ],
    }


def verdict_text(headline: dict, verdict: dict) -> str:
    """One sentence a reader can check against the numbers."""
    if not headline.get("external_evaluable"):
        return "NOT EVALUABLE: " + str(headline.get("external_reason"))
    return (
        f"The transferable model's Brier score at 18 months is "
        f"{headline['brier_18mo']:.4f} on {headline['n_cgga_eval']} CGGA patients "
        f"(95% CI {headline['ci_lower']:.4f} to {headline['ci_upper']:.4f}); "
        f"{verdict['reason']}. Reported as {verdict['status']}. A model that predicts "
        f"the cohort's own survival for everyone scores "
        f"{headline['null_model_brier_18mo']:.4f} (descriptive)."
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tcga", default=str(CX.DEFAULT_TCGA))
    parser.add_argument("--joined", default=str(CX.DEFAULT_JOINED))
    parser.add_argument("--cgga", default=str(CX.DEFAULT_CGGA),
                        help="CGGA clinical_cohort.parquet (read once for this row)")
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--cindex-result", default=str(DEFAULT_CINDEX_RESULT),
                        help="first run's JSON; used only to confirm the fit identity")
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-splits", type=int, default=N_SPLITS)
    parser.add_argument("--dry-run", action="store_true",
                        help="availability and evaluable counts only; no fit, no metric")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace an existing result file. A second evaluation is "
                             "a re-run and must be a deliberate, reported decision")
    args = parser.parse_args(argv)

    if not args.dry_run and Path(args.out).exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite {args.out}. This row is evaluated once; "
                         "pass --overwrite only as a deliberate, reported re-run.")

    print("T2 -- Brier score at 18 months on CGGA (the second registered row)")
    print("  Same three variants and the same patients as the C-index row.")
    print("  This is CGGA read number", CGGA_READ_NUMBER)
    print()
    tcga_frame, cgga_rows = CX.load_inputs(args.tcga, args.joined, args.cgga)

    if args.dry_run:
        report = CX.dry_run(tcga_frame, cgga_rows)
        print(f"CGGA rows {report['n_rows_ingested']}, unique patients "
              f"{report['n_unique_patients']}")
        for name, entry in report["variants"].items():
            print(f"  {name}: evaluable n = {entry.get('n_cgga_evaluable')}, "
                  f"absent in CGGA: {entry['absent_in_cgga']}")
        print("\nDry run: nothing fitted, no metric computed.")
        return

    result = run(tcga_frame, cgga_rows, n_boot=args.n_boot, seed=args.seed,
                 n_splits=args.n_splits, cindex_result=args.cindex_result)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")

    v1, v2, v3 = (result["variant_1_baseline_committed"],
                  result["variant_2_transferable"], result["variant_3_molecular"])
    print("1. Committed baseline (age, sex, kps, treatment_arm)")
    print(f"   NOT EVALUABLE -- {v1.get('external_reason')}\n")

    def show(label, block):
        print(label)
        if not block["external_evaluable"]:
            print(f"   NOT EVALUABLE -- {block.get('external_reason')}\n")
            return
        print(f"   CGGA n = {block['n_cgga_eval']} ({block['n_cgga_deaths']} deaths, "
              f"{block['n_cgga_at_risk_at_horizon']} at risk at 18 months); observed "
              f"KM survival at 18 months {block['observed_km_survival_at_horizon']:.4f}")
        print(f"   Brier@18mo = {block['brier_18mo']:.4f}  "
              f"95% CI [{block['ci_lower']:.4f}, {block['ci_upper']:.4f}]  "
              f"(null model {block['null_model_brier_18mo']:.4f}, descriptive)")
        for key, d in block["descriptive_other_horizons"].items():
            print(f"   {key} = {d['score']:.4f}  95% CI [{d['ci_lower']:.4f}, "
                  f"{d['ci_upper']:.4f}]  (descriptive)")
        print()

    show("2. Transferable (age, sex, treatment_arm) -- the headline", v2)
    show("3. Molecular (age, sex, treatment_arm, IDH, MGMT) -- secondary", v3)
    print(f"T2 (Brier-at-18-months criterion): {result['verdict']['status']}")
    print(f"  {verdict_text(v2, result['verdict'])}")
    check = result["fit_identity_check"]
    label = "matches" if check.get("matches_cindex_run") else check.get(
        "reason", "MISMATCH")
    print(f"\nFit identity against the C-index run: {label}")
    print(f"Wrote: {args.out}")
    if check.get("checked") and not check.get("matches_cindex_run"):
        raise SystemExit("SANITY CHECK FAILED: the TCGA fits differ from the C-index run.")


if __name__ == "__main__":
    main()
