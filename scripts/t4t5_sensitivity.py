"""Tier 1 of the sensitivity plan: re-scoring, not re-optimising.

`docs/t4t5_sensitivity_plan.md` splits section 5 of the design document in two. **This is
tier 1.** It holds each patient's configuration fixed at what the primary run chose and
re-scores it under each parameter setting, which answers *"is the recommended therapy's
benefit robust?"* — the question a clinician has, since that therapy is what the patient
would receive.

It deliberately **cannot** produce T4's distinct-configuration count or regret reduction.
Those are properties of the optimiser's behaviour, not of a fixed configuration, and
reporting them from a re-scoring pass would be answering a question nobody asked. They belong
to tier 2.

Re-scoring costs 148 µs per patient-setting against 33.3 s to re-optimise, so everything here
runs in about a minute where re-optimising the same settings would take 103 days.

    tornado      each assumed parameter across its registered range, others central
    hypercube    1000 joint Latin-hypercube draws, seed 42
    breakeven    the common Emax multiplier at which the T5 success fraction first reaches
                 0.50, reported as the hazard ratio the therapy would need

Every output is **exploratory** under section 3.6, is labelled so in its own file, and
changes no verdict. Run on synthetic patients it is a test of the machinery; run on the
registered result it is the paper's sensitivity table.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from gliostrata.optimization.config_space import Configuration
from gliostrata.prediction import response_surface as rs
from gliostrata.simulation import delivery
from gliostrata.simulation.parameters import ParameterSet, load_parameters

OUT = Path("docs/dry_runs/t4t5_sensitivity_tier1.json")
SEED = 42
N_DRAWS = 1000
SUCCESS_MONTHS = 3.0
HARM_BELOW = 0.0

# Each entry is (dotted path into the parameter document, low, high). The ranges are the
# registered ones from Design A's table; the paths are where those values actually live.
SWEPT: dict[str, tuple[list[str], float, float]] = {
    "ed50_cells": (["response", "dose_response", "ed50_cells"], 5e6, 1e8),
    "hill_h": (["response", "dose_response", "hill_h"], 0.5, 2.0),
    "translation_factor": (["response", "emax", "translation_factor"], 0.25, 1.0),
    "synergy": (["response", "combo", "synergy"], 0.8, 1.3),
    "f_iPSC-NSC": (["delivery", "cell_fraction_acting", "iPSC-NSC"], 0.1, 1.0),
    "f_MSC": (["delivery", "cell_fraction_acting", "MSC"], 0.1, 1.0),
    "f_HB1_F3_line": (["delivery", "cell_fraction_acting", "HB1_F3_line"], 0.1, 1.0),
    "tropism_slope": (["delivery", "tropism_logit_slope_per_sd"], 0.0, 1.0),
    "tropism_intercept": (["delivery", "tropism_logit_intercept"], -1.0, 1.0),
    "coverage_n0": (["delivery", "coverage_saturation_sites"], 1.5, 6.0),
    "immune_slope": (["response", "patient_modifier", "immune_tme_slope_per_sd"], 0.0, 0.3),
    "mes_slope": (["response", "patient_modifier", "mes_like_slope_per_sd"], 0.0, 0.3),
    "tau_CRAd": (["toxicity", "tau_payload", "CRAd"], 0.01, 0.06),
    "tau_prodrug": (["toxicity", "tau_payload", "prodrug_enzyme"], 0.005, 0.03),
    "tau_CAR": (["toxicity", "tau_payload", "CAR"], 0.025, 0.15),
    "dose_exponent_q": (["toxicity", "dose_exponent_q"], 1.0, 3.0),
    "site_penalty_eta": (["toxicity", "site_penalty_eta"], 0.0, 0.15),
    "delay_cost": (["toxicity", "delay_cost_log_hr"], 0.0, 0.08),
}
LOG_SCALED = {"ed50_cells", "tau_CRAd", "tau_prodrug", "tau_CAR"}


def _set(doc: dict, path: list[str], value: float) -> None:
    node = doc
    for key in path[:-1]:
        node = node[key]
    leaf = node[path[-1]]
    if isinstance(leaf, dict):
        leaf["value"] = value
    else:
        node[path[-1]] = value


def _with(params: ParameterSet, changes: dict[str, float]) -> ParameterSet:
    doc = copy.deepcopy(params.doc)
    for name, value in changes.items():
        _set(doc, SWEPT[name][0], value)
    return ParameterSet(doc=doc, path=params.path, sha256=params.sha256)


def _central(params: ParameterSet, name: str) -> float:
    node = params.doc
    for key in SWEPT[name][0]:
        node = node[key]
    return float(node["value"] if isinstance(node, dict) else node)


def score_all(cases: list[dict], params: ParameterSet) -> dict:
    """Re-score every patient's FIXED configuration under one parameter setting."""
    gains = []
    for case in cases:
        scorer = case["scorer_factory"](params)
        gains.append(scorer.gain(case["config"]))
    values = np.asarray(gains, dtype=float)
    return {
        "success_fraction": float(np.mean(values >= SUCCESS_MONTHS)),
        "harm_fraction": float(np.mean(values < HARM_BELOW)),
        "mean_gain_months": float(values.mean()),
    }


def run_tornado(cases: list[dict], params: ParameterSet) -> list[dict]:
    rows = []
    central = score_all(cases, params)
    for name, (_, low, high) in SWEPT.items():
        at_low = score_all(cases, _with(params, {name: low}))
        at_high = score_all(cases, _with(params, {name: high}))
        rows.append({
            "parameter": name, "central_value": _central(params, name),
            "low": low, "high": high,
            "at_low": at_low, "at_high": at_high,
            "swing_mean_gain": abs(at_high["mean_gain_months"] - at_low["mean_gain_months"]),
            "swing_success": abs(at_high["success_fraction"] - at_low["success_fraction"]),
        })
    rows.sort(key=lambda r: -r["swing_mean_gain"])
    return [{"central": central}] + rows


def run_hypercube(cases: list[dict], params: ParameterSet, n_draws: int = N_DRAWS) -> dict:
    """Latin-hypercube over every swept parameter, log-uniform on the scale parameters."""
    rng = np.random.default_rng(SEED)
    names = list(SWEPT)
    cuts = (rng.permuted(np.tile(np.arange(n_draws), (len(names), 1)), axis=1)
            + rng.random((len(names), n_draws))) / n_draws

    draws = []
    for j, name in enumerate(names):
        _, low, high = SWEPT[name]
        u = cuts[j]
        draws.append(np.exp(np.log(low) + u * (np.log(high) - np.log(low)))
                     if name in LOG_SCALED else low + u * (high - low))
    matrix = np.vstack(draws)

    outputs = {"success_fraction": [], "harm_fraction": [], "mean_gain_months": []}
    for i in range(n_draws):
        result = score_all(cases, _with(params, dict(zip(names, matrix[:, i]))))
        for key in outputs:
            outputs[key].append(result[key])

    def spearman(x: np.ndarray, y: np.ndarray) -> float:
        rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
        if np.std(rx) == 0 or np.std(ry) == 0:
            return 0.0
        return float(np.corrcoef(rx, ry)[0, 1])

    gain = np.asarray(outputs["mean_gain_months"])
    indices = {name: spearman(matrix[j], gain) for j, name in enumerate(names)}
    ranked = dict(sorted(indices.items(), key=lambda kv: -abs(kv[1])))
    success = np.asarray(outputs["success_fraction"])
    harm = np.asarray(outputs["harm_fraction"])
    return {
        "n_draws": n_draws, "seed": SEED,
        "mean_gain_months": {"median": float(np.median(gain)),
                             "p2.5": float(np.percentile(gain, 2.5)),
                             "p97.5": float(np.percentile(gain, 97.5))},
        "success_fraction": {"median": float(np.median(success)),
                             "p2.5": float(np.percentile(success, 2.5)),
                             "p97.5": float(np.percentile(success, 97.5))},
        "harm_fraction": {"median": float(np.median(harm)),
                          "p2.5": float(np.percentile(harm, 2.5)),
                          "p97.5": float(np.percentile(harm, 97.5))},
        "fraction_of_draws_meeting_t5": float(np.mean((success >= 0.50) & (harm <= 0.10))),
        "spearman_rank_index_on_mean_gain": ranked,
        "top_two_movers": list(ranked)[:2],
    }


def run_breakeven(cases: list[dict], params: ParameterSet) -> dict:
    """The common Emax multiplier at which T5's success fraction first reaches 0.50."""
    def at(multiplier: float) -> float:
        doc = copy.deepcopy(params.doc)
        for payload in ("CRAd", "prodrug_enzyme", "CAR"):
            entry = doc["response"]["emax"][payload]
            entry["value"] = float(entry["value"]) * multiplier
        return score_all(cases, ParameterSet(doc=doc, path=params.path,
                                             sha256=params.sha256))["success_fraction"]

    low, high = 1.0, 1.0
    for _ in range(40):
        if at(high) >= 0.50:
            break
        low, high = high, high * 1.5
    else:
        return {"reached": False,
                "note": "success fraction never reached 0.50 within a 1.5^40 multiplier",
                "max_multiplier_tried": high}

    for _ in range(40):
        mid = 0.5 * (low + high)
        if at(mid) >= 0.50:
            high = mid
        else:
            low = mid
    emax_crad = float(params.doc["response"]["emax"]["CRAd"]["value"])
    return {
        "reached": True, "multiplier": high,
        "implied_hazard_ratio_crad": math.exp(-emax_crad * high),
        "note": ("the common multiplier on every Emax at which the T5 success fraction first "
                 "reaches 0.50, holding each patient's chosen configuration fixed"),
    }


def _synthetic_cases(params: ParameterSet, n: int) -> list[dict]:
    """Patients whose configuration is already chosen, as tier 1 requires."""
    import importlib.util
    import sys
    spec = importlib.util.spec_from_file_location(
        "t4t5_controls", Path(__file__).resolve().parent / "t4t5_controls.py")
    controls = importlib.util.module_from_spec(spec)
    sys.modules["t4t5_controls"] = controls
    spec.loader.exec_module(controls)

    rng = np.random.default_rng(11)
    grid = rs.time_grid(params)
    cases = []
    for _ in range(n):
        base = controls._synthetic_scorer(params, rng)
        baseline = np.exp(-rng.uniform(0.03, 0.08) * grid)
        tropism, modifiers = base.tropism, base.modifiers
        config = Configuration(
            cell_type=str(rng.choice(["iPSC-NSC", "MSC", "HB1_F3_line"])),
            payload="combo", dose_cells=int(rng.integers(5_000_000, 50_000_000)),
            n_injection_sites=int(rng.integers(1, 9)),
            schedule=str(rng.choice(["single", "fractionated_x2", "fractionated_x3"])))
        cases.append({
            "config": config,
            "scorer_factory": (lambda p, b=baseline, t=tropism, m=modifiers:
                               controls.PatientScorer(p, b, t, m)),
        })
    return cases


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--params", type=Path, default=Path("configs/response_surface_v1.yaml"))
    parser.add_argument("--n-patients", type=int, default=60)
    parser.add_argument("--n-draws", type=int, default=N_DRAWS)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args(argv)

    params = load_parameters(args.params)
    if params.pending():
        raise SystemExit(f"Emax still pending for {params.pending()}")
    delivery.load_tropism_reference()      # refuses a reference naming CGGA

    cases = _synthetic_cases(params, args.n_patients)
    print(f"tier 1, re-scoring only: {len(cases)} patients, "
          f"{len(SWEPT)} swept parameters, {args.n_draws} draws", flush=True)

    tornado = run_tornado(cases, params)
    print("  tornado done", flush=True)
    hypercube = run_hypercube(cases, params, args.n_draws)
    print("  hypercube done", flush=True)
    breakeven = run_breakeven(cases, params)
    print("  break-even done", flush=True)

    document = {
        "what": "T4/T5 sensitivity, tier 1 (re-scoring)",
        "tier": 1,
        "exploratory": True,
        "changes_no_verdict": True,
        "plan": "docs/t4t5_sensitivity_plan.md",
        "registered_under": "docs/preregistration.md section 3.6",
        "run_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "synthetic": True,
        "parameter_sha256": params.sha256,
        "cannot_produce": ("T4's distinct-configuration count and regret reduction: both are "
                           "properties of the optimiser, not of a fixed configuration, and "
                           "belong to tier 2"),
        "n_patients": len(cases),
        "tornado": tornado,
        "hypercube": hypercube,
        "breakeven": breakeven,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=1, default=str) + "\n",
                        encoding="utf-8", newline="\n")
    print(f"\nWrote {args.out}")
    print(json.dumps({
        "central": tornado[0]["central"],
        "top_movers": hypercube["top_two_movers"],
        "draws_meeting_t5": hypercube["fraction_of_draws_meeting_t5"],
        "breakeven": breakeven,
    }, indent=1, default=str))


if __name__ == "__main__":
    main()
