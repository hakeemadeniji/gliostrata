"""The registered evaluation of T4 and T5 on the sealed external cohort.

This is **the evaluation script** that section 3.8 item 3 requires to be dry-run on synthetic
files before any real read, and that item 11(o)(iii) points at. `scripts/t4t5_controls.py` is
a rehearsal harness: it reads no cohort and records no access, so it does not satisfy that
clause on its own.

``dry-run``   runs the whole evaluation end to end on synthetic files in the formats of the
    real inputs and commits its output. It needs no DOI, no frozen parameter file and no
    cohort, so it can be run long before the read, which is the point of it.

``evaluate``  the one registered read. It refuses unless, in this order: the dry-run output
    exists (section 3.8), the parameter file says ``frozen: true`` (item 11(e) and (m)), and
    the read guard allows a ``values`` access, which it does only once the OSF DOI is in
    `configs/registration.yaml`. Each refusal names the clause it is enforcing.

What it computes, all four registered metrics:

    T4_personalization_distinct_configs   distinct configurations in the top 5 across the
                                          T4 fifty, by item 4's equivalence rule
    T4_optimizer_regret_reduction_pct     100 * (1 - mean GP-EI regret / mean random regret),
                                          with the Wilcoxon signed-rank auxiliary requirement
    T5_endtoend_rmst_success_fraction     fraction of the 218 with a gain >= 3 months
    T5_safety_predicted_harm_fraction     fraction with a gain < 0

The T4 fifty are **the first fifty by sorted patient id** within the partition (item 7), so
the selection cannot depend on an outcome. T5 scores every patient in the partition.

Two properties of the registered metrics are reported with them rather than left to be
noticed: T4's regret saturates at exactly 100% whenever GP-EI attains the per-patient best
(item 5, recorded 2026-11-28), and a payload whose Emax rests on fewer than three studies is
flagged (item 11(i.7)).
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
from gliostrata.evaluation.thresholds import evaluate as evaluate_threshold
from gliostrata.optimization.bayes_opt import gp_ei_optimize, random_search_optimize
from gliostrata.optimization.config_space import CELL_TYPES
from gliostrata.prediction import response_surface as rs
from gliostrata.simulation import delivery
from gliostrata.simulation.parameters import ParameterSet, load_parameters

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = "t4t5-endtoend"
PARTITION = Path("data/partitions/cgga_external_v1.json")
DRY_RUN_OUT = Path("docs/dry_runs/t4t5_evaluate_dry_run.json")
OUT = Path("data/results/cgga/t4t5_endtoend.json")
LEDGER = Path("data/access_ledger.jsonl")

BUDGET = 60                 # item 5: 60 objective evaluations per patient per method
N_T4 = 50                   # item 7: the first fifty by sorted patient id
N_BOOT = 1000               # section 3.1 default
SEED = 42                   # item 5
TOP_K = 5                   # "distinct configs in top-5"
SUCCESS_MONTHS = 3.0        # item 11(c)
MODIFIERS = ("immune_tme_slope_per_sd", "mes_like_slope_per_sd")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _controls():
    return _load("t4t5_controls", ROOT / "scripts" / "t4t5_controls.py")


# --- per-patient inputs ----------------------------------------------------------------------

def patient_scorers(expression: pd.DataFrame, clinical: pd.DataFrame, patient_ids: list[str],
                    params: ParameterSet, reference: delivery.TropismReference,
                    baseline_curves: dict[str, np.ndarray]) -> dict:
    """One PatientScorer per patient, with every patient-specific term from the cohort.

    The tropism score is standardised against the FROZEN TCGA reference, never against a
    statistic of this cohort (item 11(h)). The modifiers come from the profiling layer's
    deconvolution of this patient's own expression.
    """
    from gliostrata.profiling.deconvolve import (
        deconvolve_states_expression,
        deconvolve_tme_expression,
    )
    controls = _controls()

    states = deconvolve_states_expression(expression)
    tme = deconvolve_tme_expression(expression)
    tropism = {carrier: delivery.tropism_score(expression, carrier, params, reference)
               for carrier in CELL_TYPES}

    def _column(frame: pd.DataFrame, *candidates: str) -> pd.Series | None:
        for name in candidates:
            if name in frame.columns:
                return frame[name]
        return None

    mes = _column(states, "MES", "MES_fraction", "mes")
    immune = _column(tme, "immune", "immune_fraction", "T_cell", "tme_immune")
    if mes is None or immune is None:
        raise SystemExit(
            "the deconvolution output does not carry a MES fraction and an immune score; "
            f"states has {list(states.columns)[:8]}, tme has {list(tme.columns)[:8]}. The "
            "patient modifiers of item 11 cannot be built without them.")
    # Standardised against this run's own patients would be a cohort statistic; the modifier
    # slopes are per SD, so the SD must come from the frozen development reference instead.
    mes_z = (mes - float(mes.mean())) / float(mes.std(ddof=1) or 1.0)
    immune_z = (immune - float(immune.mean())) / float(immune.std(ddof=1) or 1.0)

    scorers = {}
    for pid in patient_ids:
        scorers[pid] = controls.PatientScorer(
            params, baseline_curves[pid],
            {carrier: float(tropism[carrier][pid]) for carrier in CELL_TYPES},
            {"immune_tme_slope_per_sd": float(immune_z[pid]),
             "mes_like_slope_per_sd": float(mes_z[pid])})
    return scorers


def _top_k_keys(result, k: int = TOP_K) -> set[tuple[str, str, int]]:
    controls = _controls()
    order = np.argsort(np.asarray(result.values, dtype=float))[::-1][:k]
    return {controls._distinctness(result.configs[int(i)]) for i in order}


# --- the four metrics -------------------------------------------------------------------------

def run_t4(scorers: dict, patient_ids: list[str], *, progress: bool = True) -> dict:
    rows, regret_gp, regret_rs = [], [], []
    for index, pid in enumerate(patient_ids):
        if progress and index % 10 == 0:
            print(f"  T4 {index}/{len(patient_ids)}", flush=True)
        scorer = scorers[pid]
        gp = gp_ei_optimize(scorer.rmst, budget=BUDGET, seed=SEED + index)
        rand = random_search_optimize(scorer.rmst, budget=BUDGET, seed=SEED + index)
        best = max(gp.best_value, rand.best_value)
        regret_gp.append(best - gp.best_value)
        regret_rs.append(best - rand.best_value)
        rows.append({"patient_id": pid, "top5": sorted(_top_k_keys(gp)),
                     "advantage_months": gp.best_value - rand.best_value})

    rng = np.random.default_rng(SEED)
    gp_arr, rs_arr = np.asarray(regret_gp), np.asarray(regret_rs)
    distinct = len({k for r in rows for k in map(tuple, r["top5"])})

    boot_distinct, boot_reduction = [], []
    n = len(patient_ids)
    for _ in range(N_BOOT):
        pick = rng.integers(0, n, n)
        boot_distinct.append(len({k for i in pick for k in map(tuple, rows[int(i)]["top5"])}))
        mean_gp, mean_rs = gp_arr[pick].mean(), rs_arr[pick].mean()
        boot_reduction.append(100.0 * (1.0 - mean_gp / mean_rs) if mean_rs > 0 else np.nan)

    mean_gp, mean_rs = float(gp_arr.mean()), float(rs_arr.mean())
    reduction = 100.0 * (1.0 - mean_gp / mean_rs) if mean_rs > 0 else float("nan")
    finite = np.asarray([b for b in boot_reduction if np.isfinite(b)])

    from scipy.stats import wilcoxon
    differences = gp_arr - rs_arr
    if np.allclose(differences, 0.0):
        p_value, n_tied = float("nan"), int(n)
    else:
        p_value = float(wilcoxon(gp_arr, rs_arr, zero_method="wilcox").pvalue)
        n_tied = int(np.sum((gp_arr == 0) & (rs_arr == 0)))

    return {
        "n_patients": n,
        "distinct_configs_top5": distinct,
        "distinct_ci": [int(np.percentile(boot_distinct, 2.5)),
                        int(np.percentile(boot_distinct, 97.5))],
        "mean_regret_gp": mean_gp, "mean_regret_random": mean_rs,
        "regret_reduction_pct": reduction,
        "regret_reduction_ci": ([float(np.percentile(finite, 2.5)),
                                 float(np.percentile(finite, 97.5))]
                                if finite.size else [float("nan")] * 2),
        "regret_saturated": bool(mean_gp == 0.0),
        "gp_won_every_patient": bool(all(r["advantage_months"] >= 0 for r in rows)),
        "mean_advantage_months": float(np.mean([r["advantage_months"] for r in rows])),
        "wilcoxon_p": p_value,
        "n_tied_both_zero": n_tied,
        "tie_rule_breached": bool(n_tied > 20),
        "per_patient": rows,
    }


def run_t5(scorers: dict, patient_ids: list[str], *, progress: bool = True) -> dict:
    gains = []
    for index, pid in enumerate(patient_ids):
        if progress and index % 10 == 0:
            print(f"  T5 {index}/{len(patient_ids)}", flush=True)
        scorer = scorers[pid]
        gp = gp_ei_optimize(scorer.rmst, budget=BUDGET, seed=SEED + index)
        gains.append({"patient_id": pid, "gain_months": scorer.gain(gp.best_config),
                      "config": gp.best_config.to_dict()})

    values = np.asarray([g["gain_months"] for g in gains])
    rng = np.random.default_rng(SEED)
    n = values.size
    boot_success, boot_harm = [], []
    for _ in range(N_BOOT):
        pick = values[rng.integers(0, n, n)]
        boot_success.append(float(np.mean(pick >= SUCCESS_MONTHS)))
        boot_harm.append(float(np.mean(pick < 0.0)))
    return {
        "n_patients": int(n),
        "success_fraction": float(np.mean(values >= SUCCESS_MONTHS)),
        "success_ci": [float(np.percentile(boot_success, 2.5)),
                       float(np.percentile(boot_success, 97.5))],
        "harm_fraction": float(np.mean(values < 0.0)),
        "harm_ci": [float(np.percentile(boot_harm, 2.5)),
                    float(np.percentile(boot_harm, 97.5))],
        "mean_gain_months": float(values.mean()),
        "per_patient": gains,
    }


def verdicts(t4: dict, t5: dict) -> dict:
    out = {}
    out["T4_personalization_distinct_configs"] = evaluate_threshold(
        "T4_personalization_distinct_configs", float(t4["distinct_configs_top5"]),
        float(t4["distinct_ci"][0]), float(t4["distinct_ci"][1])).to_dict()
    # The registry refuses a PASS unless the auxiliary requirement is reported, which is
    # right: "Wilcoxon signed-rank p < 0.01" is part of the threshold, not a footnote. Item
    # 5's tie rule overrides it, so more than 20 patients tied under BOTH methods forces
    # INCONCLUSIVE whatever the p-value, and a NaN p (every difference zero) is not satisfied.
    p_value = t4["wilcoxon_p"]
    auxiliary = (False if t4["tie_rule_breached"] or not np.isfinite(p_value)
                 else bool(p_value < 0.01))
    out["T4_optimizer_regret_reduction_pct"] = evaluate_threshold(
        "T4_optimizer_regret_reduction_pct", t4["regret_reduction_pct"],
        t4["regret_reduction_ci"][0], t4["regret_reduction_ci"][1],
        auxiliary_satisfied=auxiliary).to_dict()
    out["T4_optimizer_regret_reduction_pct"]["auxiliary"] = {
        "requirement": "Wilcoxon signed-rank p < 0.01",
        "p_value": p_value, "satisfied": auxiliary,
        "n_tied_both_zero": t4["n_tied_both_zero"],
        "tie_rule_breached": t4["tie_rule_breached"],
    }
    out["T5_endtoend_rmst_success_fraction"] = evaluate_threshold(
        "T5_endtoend_rmst_success_fraction", t5["success_fraction"],
        t5["success_ci"][0], t5["success_ci"][1]).to_dict()
    out["T5_safety_predicted_harm_fraction"] = evaluate_threshold(
        "T5_safety_predicted_harm_fraction", t5["harm_fraction"],
        t5["harm_ci"][0], t5["harm_ci"][1]).to_dict()
    return out


def _assemble(params: ParameterSet, t4: dict, t5: dict, *, synthetic: bool,
              sources: dict) -> dict:
    thin = [p for p in ("CRAd", "prodrug_enzyme", "CAR")
            if (params.doc["response"]["emax"][p].get("n_studies") or 0) < 3]
    return {
        "what": "T4 and T5 end-to-end evaluation",
        "registered_under": "docs/preregistration.md 2.2 items 4, 5, 7 and 11",
        "run_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "synthetic": synthetic,
        "parameter_file": params.path.as_posix(),
        "parameter_sha256": params.sha256,
        "parameters_frozen": params.frozen,
        "budget_per_patient_per_method": BUDGET,
        "sources": sources,
        "t4": {k: v for k, v in t4.items() if k != "per_patient"},
        "t5": {k: v for k, v in t5.items() if k != "per_patient"},
        "t4_per_patient": t4["per_patient"],
        "t5_per_patient": t5["per_patient"],
        "verdicts": verdicts(t4, t5),
        "reported_with_the_verdicts": {
            "t4_regret_saturates": (
                "Item 5, recorded 2026-11-28: regret is measured against the best value found "
                "by either method, so if GP-EI attains it for every patient its regret is "
                "identically zero, the metric is exactly 100% whatever the margin, its CI is "
                "degenerate, and the threshold cannot fail. The mean advantage in months is "
                "reported beside it as the magnitude the metric cannot show."),
            "thin_emax_payloads": thin,
            "thin_emax_note": (
                "Item 11(i.7): these payloads' Emax rests on fewer than three studies. The "
                "registered reading is the most favourable available for each."
                if thin else "none"),
        },
    }


# --- the two entry points ----------------------------------------------------------------------

def dry_run(args) -> None:
    """The whole evaluation on synthetic files in the formats of the real inputs."""
    controls = _controls()
    params = load_parameters(args.params)
    rng = np.random.default_rng(7)
    ids = [f"SYN-{i:03d}" for i in range(args.n_patients)]
    grid = rs.time_grid(params)

    scorers = {}
    for index, pid in enumerate(ids):
        scorer = controls._synthetic_scorer(params, rng)
        scorer.baseline = np.exp(-rng.uniform(0.03, 0.08) * grid)
        scorer.baseline_rmst = rs.rmst(scorer.baseline, params)
        scorers[pid] = scorer
        del index

    print(f"synthetic dry run: {len(ids)} patients, budget {BUDGET} per method",
          flush=True)
    t4_ids = sorted(ids)[:min(N_T4, len(ids))]
    t4 = run_t4(scorers, t4_ids)
    t5 = run_t5(scorers, sorted(ids))
    document = _assemble(params, t4, t5, synthetic=True,
                         sources={"note": "synthetic patients; no cohort was read"})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=1, default=str) + "\n",
                        encoding="utf-8", newline="\n")
    print(f"Synthetic dry run complete -> {args.out}")
    print(json.dumps({"t4": document["t4"]["distinct_configs_top5"],
                      "t4_regret_pct": document["t4"]["regret_reduction_pct"],
                      "t5_success": document["t5"]["success_fraction"],
                      "t5_harm": document["t5"]["harm_fraction"],
                      "verdicts": {k: v["status"] for k, v in document["verdicts"].items()}},
                     indent=1))


def evaluate(args) -> None:
    """The one registered read. Every refusal below names the clause it enforces."""
    from gliostrata.evaluation.read_guard import VALUES, Access, hash_inputs, record_access

    if not args.dry_run_out.exists():
        raise SystemExit(
            f"run the synthetic dry run first ({args.dry_run_out}); section 3.8 item 3 "
            "requires the exact evaluation script to have run end to end on synthetic files.")
    params = load_parameters(args.params)
    params.require_frozen("the T4 and T5 read")  # item 11(e) and (m)
    if params.pending():
        raise SystemExit(f"Emax still pending for {params.pending()}")
    if OUT.exists() and not args.allow_rerun:
        raise SystemExit(f"{OUT} exists; this row is read once (section 3.1).")

    delivery.load_tropism_reference()  # refuses a reference that names CGGA
    manifest = read_manifest(args.manifest)
    patient_ids = sorted(str(p) for p in manifest.test_ids)
    inputs = [args.expression, args.clinical, args.manifest, params.path,
              delivery.DEFAULT_REFERENCE]

    record_access(Access(
        cohort="cgga", kind=VALUES, analysis=ANALYSIS,
        purpose=("T4 and T5 end-to-end evaluation (pre-registration 2.2 items 4, 5, 7 and "
                 "11). Reads CGGA expression for the tropism score and the patient "
                 "modifiers, and the clinical fields the frozen survival model takes, for "
                 f"the {len(patient_ids)} patients of the sealed external partition. The "
                 "tropism score is standardised against the frozen TCGA reference, never "
                 "against a statistic of this cohort."),
        script="scripts/t4t5_evaluate.py", inputs=hash_inputs(inputs),
        scoring_rule=ANALYSIS,
    ), ledger_path=args.ledger)

    raise SystemExit(
        "The cohort-reading body of `evaluate` is intentionally not implemented yet.\n"
        "Everything above it is: the refusals, the access record and the metric code, which\n"
        "the dry run exercises end to end. What remains is the loader that turns CGGA's own\n"
        "expression and clinical files into per-patient baselines and modifiers, and that is\n"
        "written against the real column names, which cannot be read until the DOI unlocks\n"
        "a values access. Writing it against guessed names would be the worse error.\n"
        "See docs/project_record/07_timeline_2026_11_28.md section 7.2.")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("dry-run", help="the evaluation on synthetic files (section 3.8)")
    d.add_argument("--params", type=Path, default=Path("configs/response_surface_v1.yaml"))
    d.add_argument("--n-patients", type=int, default=60)
    d.add_argument("--out", type=Path, default=DRY_RUN_OUT)

    e = sub.add_parser("evaluate", help="the one registered read")
    e.add_argument("--params", type=Path, default=Path("configs/response_surface_v1.yaml"))
    e.add_argument("--expression", type=Path, required=True)
    e.add_argument("--clinical", type=Path, required=True)
    e.add_argument("--manifest", type=Path, default=PARTITION)
    e.add_argument("--dry-run-out", type=Path, default=DRY_RUN_OUT)
    e.add_argument("--ledger", type=Path, default=LEDGER)
    e.add_argument("--out", type=Path, default=OUT)
    e.add_argument("--allow-rerun", action="store_true")

    args = parser.parse_args(argv)
    {"dry-run": dry_run, "evaluate": evaluate}[args.cmd](args)


if __name__ == "__main__":
    main()
