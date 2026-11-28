"""The four pre-read controls of pre-registration item 11(o).

Item 11(o) requires these before any T4 or T5 read of CGGA, and says that a control which
does not give its stated answer **stops the read**. This script runs all four and exits
non-zero if any fails, so a failure cannot be stepped over quietly.

    (i)   no-modifier     every patient modifier zero -> the optimiser's answer cannot vary
                          by patient, so the distinct-configuration count must be exactly 1.
    (ii)  null-payload    Emax = 0 for every payload -> the benefit term E must be exactly 0
                          for every patient and configuration, and the gain must be the one
                          the toxicity term alone implies, which is negative wherever that
                          term is positive. (Corrected 2026-11-27: the clause first required a
                          gain of exactly 0, which item 11(k)'s toxicity term makes
                          unreachable.)
    (iii) synthetic       the whole pipeline on synthetic patients, no cohort touched.
    (iv)  tcga-rehearsal  the whole T4 and T5 pipeline on the held-out TCGA split, before any
                          CGGA read, so that bugs are found without spending the one
                          registered read. Its numbers are NOT results and are never reported
                          as such.

The survival model is the one T2 already validated: `VARIANT_2` fitted by `fit_variant` from
`scripts/validate_cox_cgga.py`, not a new fit, so the Stupp baseline here is the model item
11(g) registers.

**No CGGA path is opened.** Nothing here reads the sealed cohort, so no read-guard access is
recorded; the guard would refuse it anyway until the DOI is in the config.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.evaluation.partitions import read_manifest
from gliostrata.optimization.bayes_opt import gp_ei_optimize, random_search_optimize
from gliostrata.optimization.config_space import (
    CELL_TYPES,
    PAYLOADS,
    SCHEDULES,
    Configuration,
)
from gliostrata.prediction import response_surface as rs
from gliostrata.simulation import delivery
from gliostrata.simulation.parameters import ParameterSet, load_parameters

ROOT = Path(__file__).resolve().parents[1]
OUT = Path("docs/dry_runs")
TCGA_COHORT = Path("data/processed/tcga_gbm/feature_store.parquet")
TCGA_CLINICAL = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
TCGA_CALLS = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
PARTITION = Path("data/partitions/tcga_gbm_clinical_v2.json")
MODIFIERS = ("immune_tme_slope_per_sd", "mes_like_slope_per_sd")


def _load_validate_cox():
    """The T2 validation module, loaded by path: it is a script, not a package."""
    spec = importlib.util.spec_from_file_location(
        "validate_cox_cgga", ROOT / "scripts" / "validate_cox_cgga.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["validate_cox_cgga"] = module
    spec.loader.exec_module(module)
    return module


class PatientScorer:
    """Turns (patient, configuration) into a predicted RMST, per Design A.

    The patient's Stupp baseline is the survival model's prediction with the arm set to
    radiotherapy plus temozolomide (item 11(a)). The configuration shifts the curve by
    exp(-E + T) (item 11(b)).
    """

    def __init__(self, params: ParameterSet, baseline_survival: np.ndarray,
                 tropism: dict[str, float], modifiers: dict[str, float]):
        self.params = params
        self.baseline = np.asarray(baseline_survival, dtype=float)
        self.tropism = tropism          # per carrier
        self.modifiers = modifiers      # per modifier name
        self.baseline_rmst = rs.rmst(self.baseline, params)

    def effect(self, config: Configuration, *, null_payload: bool = False) -> rs.EffectRecord:
        terms = delivery.delivery_terms(
            dose_cells=config.dose_cells, carrier=config.cell_type,
            tropism=self.tropism[config.cell_type],
            n_injection_sites=config.n_injection_sites, schedule=config.schedule,
            params=self.params)
        return rs.configuration_effect(
            payload=config.payload, carrier=config.cell_type,
            dose_cells=config.dose_cells, n_injection_sites=config.n_injection_sites,
            n_eff=float(terms["n_eff"]), modifiers=self.modifiers,
            sigma=float(terms["schedule_factor"]), cov=float(terms["coverage"]),
            params=self.params, null_payload=null_payload)

    def rmst(self, config: Configuration, *, null_payload: bool = False) -> float:
        effect = self.effect(config, null_payload=null_payload)
        curve = rs.configuration_curve(self.baseline, effect.log_hazard_offset)
        return rs.rmst(curve, self.params)

    def gain(self, config: Configuration, *, null_payload: bool = False) -> float:
        return self.rmst(config, null_payload=null_payload) - self.baseline_rmst


def _distinctness(config: Configuration) -> tuple[str, str, int]:
    """Item 4's registered rule: carrier, payload and NUMBER of sites; dose and schedule
    ignored. Two configurations are the same if and only if these three agree."""
    return (config.cell_type, config.payload, config.n_injection_sites)


def _top_k_keys(result, k: int = 5) -> set[tuple[str, str, int]]:
    """The distinctness keys of a patient's top k configurations.

    T4's registered metric counts distinct configurations **in the top 5** across the
    patients, not the single best per patient (item 2's threshold table, with item 4's
    equivalence rule). Counting only the best would be a different, stricter quantity, so a
    rehearsal that did that would not be rehearsing T4.
    """
    order = np.argsort(np.asarray(result.values, dtype=float))[::-1][:k]
    return {_distinctness(result.configs[int(i)]) for i in order}


def _synthetic_scorer(params: ParameterSet, rng: np.random.Generator,
                      *, flat_modifiers: bool = False) -> PatientScorer:
    grid = rs.time_grid(params)
    hazard = float(rng.uniform(0.03, 0.08))
    baseline = np.exp(-hazard * grid)
    tropism = {carrier: float(rng.normal(0.0, 1.0)) for carrier in CELL_TYPES}
    modifiers = ({name: 0.0 for name in MODIFIERS} if flat_modifiers
                 else {name: float(rng.normal(0.0, 1.0)) for name in MODIFIERS})
    return PatientScorer(params, baseline, tropism, modifiers)


# --- (i) no-modifier ---------------------------------------------------------------------

def control_no_modifier(params: ParameterSet, n_patients: int = 20,
                        budget: int = 60) -> dict:
    """With every modifier at zero the optimiser's answer cannot vary by patient."""
    rng = np.random.default_rng(11)
    chosen = []
    for _ in range(n_patients):
        scorer = _synthetic_scorer(params, rng, flat_modifiers=True)
        # The patient's tropism still varies, so hold it flat too: the control is about the
        # response surface having no patient-specific term left, not about the carrier.
        scorer.tropism = {carrier: 0.0 for carrier in CELL_TYPES}
        result = gp_ei_optimize(scorer.rmst, budget=budget, seed=7)
        chosen.append(_distinctness(result.best_config))
    distinct = sorted(set(chosen))
    ok = len(distinct) == 1
    return {"control": "no_modifier", "must_yield": "distinct-configuration count of 1",
            "n_patients": n_patients, "n_distinct": len(distinct),
            "distinct": [list(d) for d in distinct], "passed": ok}


# --- (ii) null-payload -------------------------------------------------------------------

def control_null_payload(params: ParameterSet, n_patients: int = 20) -> dict:
    """Emax = 0 must give E exactly 0, and a gain equal to the toxicity term's alone."""
    rng = np.random.default_rng(12)
    configs = [Configuration(cell_type=c, payload=p, dose_cells=d, n_injection_sites=n,
                             schedule=s)
               for c in CELL_TYPES for p in PAYLOADS for d in (1_000_000, 50_000_000)
               for n in (1, 8) for s in SCHEDULES
               if (c, p) not in {("HB1_F3_line", "CAR"), ("MSC", "CRAd")}]
    worst_benefit, gains, mismatches = 0.0, [], 0
    for _ in range(n_patients):
        scorer = _synthetic_scorer(params, rng)
        for config in configs:
            effect = scorer.effect(config, null_payload=True)
            worst_benefit = max(worst_benefit, abs(effect.benefit))
            gain = scorer.gain(config, null_payload=True)
            gains.append(gain)
            # The gain must be exactly the one the toxicity term alone implies.
            toxicity_only = rs.rmst(
                rs.configuration_curve(scorer.baseline, effect.toxicity), params
            ) - scorer.baseline_rmst
            if abs(gain - toxicity_only) > 1e-12:
                mismatches += 1
    ok = worst_benefit == 0.0 and mismatches == 0 and max(gains) < 0
    return {"control": "null_payload",
            "must_yield": ("benefit exactly 0; gain equal to the toxicity term's alone, "
                           "negative wherever that term is positive"),
            "n_patients": n_patients, "n_configs": len(configs),
            "max_abs_benefit": worst_benefit, "n_gain_mismatches": mismatches,
            "max_gain": max(gains), "min_gain": min(gains), "passed": bool(ok)}


# --- (iii) synthetic dry run -------------------------------------------------------------

def control_synthetic(params: ParameterSet, n_patients: int = 10, budget: int = 60) -> dict:
    """The whole pipeline end to end on synthetic patients. No cohort is touched."""
    rng = np.random.default_rng(13)
    rows, regret_gp, regret_rs = [], [], []
    for index in range(n_patients):
        scorer = _synthetic_scorer(params, rng)
        gp = gp_ei_optimize(scorer.rmst, budget=budget, seed=100 + index)
        rand = random_search_optimize(scorer.rmst, budget=budget, seed=100 + index)
        best = max(gp.best_value, rand.best_value)
        regret_gp.append(best - gp.best_value)
        regret_rs.append(best - rand.best_value)
        gain = scorer.gain(gp.best_config)
        rows.append({"patient": index, "best_rmst": gp.best_value, "gain_months": gain,
                     "config": gp.best_config.to_dict(),
                     "top5_keys": sorted(_top_k_keys(gp)),
                     "advantage_months": gp.best_value - rand.best_value,
                     "distinct_key": list(_distinctness(gp.best_config))})
    mean_gp, mean_rs = float(np.mean(regret_gp)), float(np.mean(regret_rs))
    reduction = 100.0 * (1.0 - mean_gp / mean_rs) if mean_rs > 0 else float("nan")
    gains = [r["gain_months"] for r in rows]
    return {"control": "synthetic_dry_run",
            "must_yield": "completes and writes a verdict of the registered shape",
            "n_patients": n_patients, "budget": budget,
            "t4_distinct_configs_top5": len(
                {k for r in rows for k in map(tuple, r["top5_keys"])}),
            "t4_distinct_configs_best_only": len({tuple(r["distinct_key"]) for r in rows}),
            "t4_regret_reduction_pct": reduction,
            "t4_regret_saturated": bool(mean_gp == 0.0),
            "t4_mean_advantage_months": float(np.mean([r["advantage_months"] for r in rows])),
            "t4_gp_won_every_patient": bool(all(
                r["advantage_months"] >= 0 for r in rows)),
            "t5_success_fraction": float(np.mean([g >= 3.0 for g in gains])),
            "t5_harm_fraction": float(np.mean([g < 0.0 for g in gains])),
            "mean_gain_months": float(np.mean(gains)),
            # The per-patient values, so a reader and a figure can see the distribution
            # rather than a mean standing in for it.
            "gains_months": [float(g) for g in gains],
            "passed": bool(np.isfinite(mean_gp) and np.isfinite(reduction)),
            "note": "synthetic patients; these numbers are not results"}


# --- (iv) TCGA rehearsal ------------------------------------------------------------------

def control_tcga(params: ParameterSet, budget: int = 60, max_patients: int = 50) -> dict:
    """The whole T4/T5 pipeline on the held-out TCGA split, before any CGGA read."""
    module = _load_validate_cox()
    for path in (TCGA_CLINICAL, TCGA_CALLS, PARTITION):
        if not path.exists():
            return {"control": "tcga_rehearsal", "passed": False,
                    "skipped": f"{path} not found"}
    cohort = pd.read_parquet(TCGA_CLINICAL)
    joined = pd.read_parquet(TCGA_CALLS, columns=["patient_id", *module.CALL_COLUMNS])
    tcga = module.load_tcga(cohort, joined)
    fitted = module.fit_variant(tcga, module.VARIANT_2)

    manifest = read_manifest(PARTITION)
    held_out = [str(p) for p in manifest.test_ids]
    frame = tcga[tcga["patient_id"].astype(str).isin(held_out)].copy()
    frame, _ = module.complete_cases(frame, module.VARIANT_2)
    frame = frame.head(max_patients)
    if frame.empty:
        return {"control": "tcga_rehearsal", "passed": False,
                "skipped": "no complete-case held-out TCGA patients"}

    # Stupp baseline: the same encoding the fit used, with the arm set to RT+TMZ.
    stupp = frame.copy()
    arm_levels = fitted.categorical_levels.get("treatment_arm", [])
    target = next((lvl for lvl in arm_levels if "tmz" in str(lvl).lower()
                   or "temo" in str(lvl).lower()), None)
    if target is None:
        target = next((lvl for lvl in arm_levels if "rt" in str(lvl).lower()), None)
    if target is not None:
        stupp["treatment_arm"] = target
    design, _, _, report = module.build_design_matrix(stupp, module.VARIANT_2, impute="drop")
    design = pd.DataFrame(design, columns=report.feature_names).reindex(
        columns=fitted.design_columns, fill_value=0.0).to_numpy()
    log_hazard = fitted.model.predict_log_hazard(design)

    grid = rs.time_grid(params)
    cumhaz = fitted.model.baseline_cumulative_hazard(grid)
    reference = delivery.load_tropism_reference()
    rng = np.random.default_rng(14)

    rows, regret_gp, regret_rs = [], [], []
    for index, lp in enumerate(log_hazard):
        baseline = np.exp(-cumhaz * np.exp(float(lp)))
        baseline = np.clip(baseline, 0.0, 1.0)
        # No expression for these patients in this rehearsal, so the tropism score is drawn
        # from the frozen reference's own scale. It exercises the code path; it is not a
        # measurement, and it is never reported as one.
        tropism = {carrier: float(rng.normal(0.0, 1.0)) for carrier in CELL_TYPES}
        modifiers = {name: float(rng.normal(0.0, 1.0)) for name in MODIFIERS}
        scorer = PatientScorer(params, baseline, tropism, modifiers)
        gp = gp_ei_optimize(scorer.rmst, budget=budget, seed=200 + index)
        rand = random_search_optimize(scorer.rmst, budget=budget, seed=200 + index)
        best = max(gp.best_value, rand.best_value)
        regret_gp.append(best - gp.best_value)
        regret_rs.append(best - rand.best_value)
        rows.append({"gain_months": scorer.gain(gp.best_config),
                     "baseline_rmst": scorer.baseline_rmst,
                     "top5_keys": sorted(_top_k_keys(gp)),
                     "advantage_months": gp.best_value - rand.best_value,
                     "distinct_key": list(_distinctness(gp.best_config))})
    mean_gp, mean_rs = float(np.mean(regret_gp)), float(np.mean(regret_rs))
    reduction = 100.0 * (1.0 - mean_gp / mean_rs) if mean_rs > 0 else float("nan")
    gains = [r["gain_months"] for r in rows]
    return {"control": "tcga_rehearsal",
            "must_yield": "completes without error; numbers are not results",
            "reference_scale": reference.source,
            "survival_model": {"features": list(module.VARIANT_2),
                               "n_fit": fitted.n_fit,
                               "internal_cindex": fitted.internal_cindex},
            "n_patients": int(len(rows)),
            "mean_baseline_rmst_months": float(np.mean([r["baseline_rmst"] for r in rows])),
            "t4_distinct_configs_top5": len(
                {k for r in rows for k in map(tuple, r["top5_keys"])}),
            "t4_distinct_configs_best_only": len({tuple(r["distinct_key"]) for r in rows}),
            "t4_regret_reduction_pct": reduction,
            "t4_regret_saturated": bool(mean_gp == 0.0),
            "t4_mean_advantage_months": float(np.mean([r["advantage_months"] for r in rows])),
            "t4_gp_won_every_patient": bool(all(
                r["advantage_months"] >= 0 for r in rows)),
            "t5_success_fraction": float(np.mean([g >= 3.0 for g in gains])),
            "t5_harm_fraction": float(np.mean([g < 0.0 for g in gains])),
            "mean_gain_months": float(np.mean(gains)),
            "gains_months": [float(g) for g in gains],
            "passed": bool(np.isfinite(reduction)),
            "note": ("a rehearsal on development data with simulated tropism inputs. NOT a "
                     "result, never reported as one, and no CGGA path was opened.")}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("all", "no-modifier", "null-payload", "synthetic",
                                         "tcga-rehearsal"))
    parser.add_argument("--budget", type=int, default=60)
    parser.add_argument("--out", type=Path, default=OUT / "t4t5_controls.json")
    args = parser.parse_args(argv)

    params = load_parameters()
    print(f"parameters: {params.path} (frozen: {params.frozen})")
    if params.pending():
        raise SystemExit(f"Emax still pending for {params.pending()}; enter the values first")

    runners = {
        "no-modifier": lambda: control_no_modifier(params, budget=args.budget),
        "null-payload": lambda: control_null_payload(params),
        "synthetic": lambda: control_synthetic(params, budget=args.budget),
        "tcga-rehearsal": lambda: control_tcga(params, budget=args.budget),
    }
    steps = list(runners) if args.step == "all" else [args.step]
    results = []
    for name in steps:
        print(f"\n--- {name} ---")
        outcome = runners[name]()
        results.append(outcome)
        print(json.dumps({k: v for k, v in outcome.items()
                          if k not in ("distinct",)}, indent=1, default=str))

    document = {
        "registered_under": "docs/preregistration.md section 2.2 item 11(o)",
        "run_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "parameter_file": params.path.as_posix(),
        "parameter_sha256": params.sha256,
        "frozen": params.frozen,
        "controls": results,
        "all_passed": all(r.get("passed") for r in results),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=1, default=str) + "\n",
                        encoding="utf-8", newline="\n")
    print(f"\nWrote {args.out}")
    if not document["all_passed"]:
        failed = [r["control"] for r in results if not r.get("passed")]
        raise SystemExit(f"CONTROL FAILED: {failed}. Item 11(o): a control that does not give "
                         "its stated answer stops the read.")
    print("All controls passed.")


if __name__ == "__main__":
    main()
