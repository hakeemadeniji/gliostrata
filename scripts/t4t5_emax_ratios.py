"""Compute ratio_p and Emax_p from the final extraction, exactly as item 11(i) registers it.

    Emax_p  = translation_factor * ln(ratio_p)            (i.1); the factor is assumed
    ratio_p = median over qualifying studies of           (i.2)
              treated median survival / control median survival

Two sets are computed for every payload:

* **Primary**, under the disclosed deviation of item 11(i.7): a study qualifies only if its
  control arm received no agent intended to have anti-tumour activity, that is
  ``untreated``, ``vehicle`` or ``mock_or_reporter_vector``. ``active_therapeutic`` and
  ``uncertain`` do not qualify.
* **Exploratory**, under (i.3) as literally written: any arm the paper designates as its
  control, so every extracted study with a ratio.

If a payload's primary set is empty, (i.7)'s fallback applies and that payload takes the
exploratory set's median, flagged ``literal_fallback``. If both are empty, (i.5) applies:
Emax_p = 0 and the payload is reported as unevidenced.

The script reads the translation factor from `configs/response_surface_v1.yaml` rather than
hard-coding it, and it **does not write** that file: entering a value there is the freeze step
of item 11(e) and (m) and is done deliberately, not as a side effect of a calculation.

No cohort data is touched, so no read-guard access is recorded.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

import yaml

SEARCH_DIR = Path("docs/t4t5_emax_search")
FINAL = SEARCH_DIR / "screening" / "extraction_final.jsonl"
OUT_JSON = SEARCH_DIR / "emax_result.json"
OUT_MD = SEARCH_DIR / "emax_result.md"
PARAMS = Path("configs/response_surface_v1.yaml")
PAYLOADS = ("CRAd", "prodrug_enzyme", "CAR")
PRIMARY_CONTROL_KINDS = frozenset({"untreated", "vehicle", "mock_or_reporter_vector"})
# The 3-month gain bar of item 11(c) corresponds to roughly this hazard ratio on the derived
# Stupp curve (docs/t4t5_model_design.md section 2.1). Reported for context only: it is not
# a threshold and nothing here is compared against it to decide anything.
T5_BREAK_EVEN_HR = 0.59


def _load_final() -> list[dict]:
    if not FINAL.exists():
        raise SystemExit(
            f"{FINAL} not found. It is the adjudicated union of the abstract-level and "
            "full-text extractions; build it before computing a ratio.")
    return [json.loads(line) for line in FINAL.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _translation_factor() -> float:
    doc = yaml.safe_load(PARAMS.read_text(encoding="utf-8"))
    entry = doc["response"]["emax"]["translation_factor"]
    return float(entry["value"] if isinstance(entry, dict) else entry)


def _ratio(row: dict) -> float | None:
    if row.get("ratio") is not None:
        return float(row["ratio"])
    treated, control = row.get("treated_median"), row.get("control_median")
    if treated is None or not control:
        return None
    return float(treated) / float(control)


def _eligible(row: dict) -> tuple[bool, str]:
    """P1 to P5 as confirmed by the extraction, before (i.7)'s control-arm narrowing."""
    if row.get("p2_intracranial") == "no":
        return False, "P2: not an intracranial rodent glioma model"
    if row.get("p2_intracranial") not in ("yes",):
        return False, "P2: not confirmed"
    if row.get("p4_payload_match") != "yes":
        return False, "P4: payload class not confirmed"
    if _ratio(row) is None:
        return False, "P3: no treated-and-control median survival obtainable"
    if row.get("arm_is_combination"):
        return False, "treated arm is a combination, not the payload alone"
    return True, ""


def _summarise(ratios: list[float], factor: float) -> dict:
    if not ratios:
        return {"n": 0, "ratio": None, "emax": None, "hazard_ratio": None}
    ratio = statistics.median(ratios)
    emax = factor * math.log(ratio)
    return {
        "n": len(ratios),
        "ratios": sorted(round(r, 4) for r in ratios),
        "ratio": round(ratio, 6),
        "emax": round(emax, 6),
        "hazard_ratio": round(math.exp(-emax), 6),
    }


def compute() -> dict:
    rows = _load_final()
    factor = _translation_factor()
    result: dict = {
        "registered_under": "docs/preregistration.md 2.2 item 11 (i.1) to (i.7)",
        "computed_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "translation_factor": factor,
        "translation_factor_status": "assumed, declared as assumed wherever reported",
        "n_rows": len(rows),
        "payloads": {},
        "t5_break_even_hazard_ratio_for_context": T5_BREAK_EVEN_HR,
    }

    for payload in PAYLOADS:
        mine = [r for r in rows if r["payload"] == payload]
        eligible, excluded = [], []
        for row in mine:
            ok, why = _eligible(row)
            (eligible if ok else excluded).append(row if ok else {**row, "why": why})

        literal = [_ratio(r) for r in eligible]
        primary_rows = [r for r in eligible
                        if r.get("control_arm_kind") in PRIMARY_CONTROL_KINDS]
        primary = [_ratio(r) for r in primary_rows]

        primary_summary = _summarise([r for r in primary if r is not None], factor)
        literal_summary = _summarise([r for r in literal if r is not None], factor)

        if primary_summary["n"] > 0:
            chosen, rule = primary_summary, "narrowed"
        elif literal_summary["n"] > 0:
            chosen, rule = literal_summary, "literal_fallback"
        else:
            chosen, rule = {"n": 0, "ratio": None, "emax": 0.0, "hazard_ratio": 1.0}, "unevidenced"

        result["payloads"][payload] = {
            "n_carried": len(mine),
            "n_eligible_p1_to_p5": len(eligible),
            "n_excluded": len(excluded),
            "exclusions": sorted({e["why"] for e in excluded}),
            "control_arm_kinds": {
                kind: sum(1 for r in eligible if r.get("control_arm_kind") == kind)
                for kind in sorted({r.get("control_arm_kind") for r in eligible} - {None})
            },
            "primary_narrowed": primary_summary,
            "exploratory_literal": literal_summary,
            "n_dropped_by_narrowing": literal_summary["n"] - primary_summary["n"],
            "control_arm_rule": rule,
            "emax": chosen["emax"],
            "ratio": chosen["ratio"],
            "hazard_ratio": chosen["hazard_ratio"],
            "studies": sorted(
                {r["pmid"] for r in primary_rows} if rule == "narrowed"
                else {r["pmid"] for r in eligible}),
        }
    return result


def _markdown(result: dict) -> str:
    lines = [
        "# Emax, computed under item 11(i)", "",
        f"Computed {result['computed_at_utc']} from "
        "`docs/t4t5_emax_search/screening/extraction_final.jsonl`.", "",
        "`Emax_p = f * ln(ratio_p)`, with `f = "
        f"{result['translation_factor']}` — **assumed**, and declared as assumed wherever a "
        "T4 or T5 number is reported. `ratio_p` is the **median** over qualifying studies, "
        "never one study chosen by hand.", "",
        "**Primary** applies the disclosed deviation of item 11(i.7): the control arm must "
        "have received no agent intended to have anti-tumour activity. **Exploratory** "
        "applies (i.3) as literally written. Where the primary set is empty the payload falls "
        "back to the exploratory set, flagged as such.", "",
        "| Payload | Rule | n | ratio | Emax | implied HR | dropped by narrowing |",
        "|---|---|---|---|---|---|---|",
    ]
    for payload, info in result["payloads"].items():
        hr = info["hazard_ratio"]
        lines.append(
            f"| {payload} | {info['control_arm_rule']} | {info['primary_narrowed']['n']} | "
            f"{info['ratio'] if info['ratio'] is not None else '—'} | "
            f"{info['emax'] if info['emax'] is not None else '—'} | "
            f"{hr if hr is not None else '—'} | {info['n_dropped_by_narrowing']} |")
    lines += [
        "", "## Both sets side by side", "",
        "| Payload | Primary n | Primary HR | Exploratory n | Exploratory HR |",
        "|---|---|---|---|---|",
    ]
    for payload, info in result["payloads"].items():
        p, e = info["primary_narrowed"], info["exploratory_literal"]
        lines.append(f"| {payload} | {p['n']} | {p['hazard_ratio'] or '—'} | "
                     f"{e['n']} | {e['hazard_ratio'] or '—'} |")
    lines += [
        "", "## How each payload's set was reduced", "",
    ]
    for payload, info in result["payloads"].items():
        lines.append(f"**{payload}** — {info['n_carried']} carried past stage 1, "
                     f"{info['n_eligible_p1_to_p5']} eligible under P1 to P5, "
                     f"{info['primary_narrowed']['n']} in the primary set.")
        if info["exclusions"]:
            lines += ["", *(f"  - excluded: {why}" for why in info["exclusions"]), ""]
        if info["control_arm_kinds"]:
            lines.append(f"  - control arms: {info['control_arm_kinds']}")
        lines.append(f"  - studies used: {', '.join(info['studies']) or 'none'}")
        lines.append("")
    lines += [
        f"For context only, the 3-month gain bar of item 11(c) sits near a hazard ratio of "
        f"{result['t5_break_even_hazard_ratio_for_context']}. Nothing above is compared "
        "against it to decide anything: T5's verdict comes from the registry in section 6, "
        "on the sealed cohort, after the freeze.", "",
    ]
    return "\n".join(lines)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("compute",))
    parser.parse_args(argv)
    result = compute()
    OUT_JSON.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8", newline="\n")
    OUT_MD.write_text(_markdown(result), encoding="utf-8", newline="\n")
    for payload, info in result["payloads"].items():
        print(f"{payload:>15}  rule={info['control_arm_rule']:<17} "
              f"n={info['primary_narrowed']['n']}/{info['exploratory_literal']['n']}  "
              f"ratio={info['ratio']}  Emax={info['emax']}  HR={info['hazard_ratio']}")
    print(f"\nWrote {OUT_JSON} and {OUT_MD}")
    print("Nothing was written to configs/response_surface_v1.yaml: entering a value there "
          "is the freeze step of item 11(e) and (m).")


if __name__ == "__main__":
    main()
