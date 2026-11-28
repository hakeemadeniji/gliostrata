"""Adjudicate the two extractions into the final set the Emax ratio is computed from.

Pre-registration item 11 (i.4) requires dual independent screening with disagreements
adjudicated. This script does the adjudication by rule, so that no row's fate depends on a
view taken after the disagreement was visible:

* **A number survives only where both extractors read the same number.** A disagreement
  yields a null and the record is reported as unresolved, never averaged or picked between.
* **The control arm's kind resolves conservatively.** If either extractor says
  ``uncertain`` or ``active_therapeutic``, that is the result, and the study leaves the
  item 11(i.7) primary set. Two different inert kinds (untreated, vehicle,
  mock_or_reporter_vector) do not change eligibility, so the disagreement is recorded and
  the row stays in the primary set.
* **P2 and P4 must be affirmed twice.** "yes" only where both extractors say "yes";
  anything else is ``uncertain`` and the study is excluded. Not-refuted is not confirmed.
* **A record with no full text keeps its abstract-level reading**, with its control arm
  recorded as ``uncertain`` — because an abstract that does not describe the control arm is
  not evidence that the control was inert. Those records are excluded from the primary set
  until their full texts are obtained (`docs/t4t5_emax_search/manual_download_needed.md`).

Writes `docs/t4t5_emax_search/screening/extraction_final.jsonl` and an adjudication summary.
No cohort data is touched, so no read-guard access is recorded.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SCREEN = Path("docs/t4t5_emax_search/screening")
OUT = SCREEN / "extraction_final.jsonl"
SUMMARY = SCREEN / "extraction_final_summary.json"
INERT = ("untreated", "vehicle", "mock_or_reporter_vector")
BLOCKING = ("active_therapeutic", "uncertain")


def _read(name: str) -> list[dict]:
    path = SCREEN / name
    if not path.exists():
        raise SystemExit(f"{path} not found")
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _agree(a, b):
    """The value both extractors read, or None where they differ."""
    return a if a == b else None


def _control_kind(a: str | None, b: str | None) -> tuple[str, bool]:
    """Conservative resolution. Returns (kind, extractors_agreed)."""
    if a == b:
        return (a or "uncertain"), True
    for blocking in BLOCKING:
        if blocking in (a, b):
            return blocking, False
    if a in INERT and b in INERT:
        return a, False  # both inert: eligibility is the same either way
    return "uncertain", False


def _affirmed(a: str | None, b: str | None) -> str:
    return "yes" if a == "yes" and b == "yes" else (
        "no" if "no" in (a, b) else "uncertain")


def finalise() -> None:
    carried = [r for r in _read("adjudication.jsonl") if r["carried_to_stage_2"]]
    abstracts = {(r["pmid"], r["payload"]): r for r in _read("extraction_merged.jsonl")}
    ft_a = {(r["pmid"], r["payload"]): r for r in _read("fulltext_extract_a2.jsonl")}
    ft_b = {(r["pmid"], r["payload"]): r for r in _read("fulltext_extract_b2.jsonl")}

    rows, disagreements = [], []
    for item in carried:
        key = (item["pmid"], item["payload"])
        abstract = abstracts.get(key, {})
        a, b = ft_a.get(key), ft_b.get(key)

        if a and b:
            treated = _agree(a["treated_median"], b["treated_median"])
            control = _agree(a["control_median"], b["control_median"])
            unit = _agree(a["unit"], b["unit"])
            kind, kind_agreed = _control_kind(a["control_arm_kind"], b["control_arm_kind"])
            combination = bool(a.get("arm_is_combination") or b.get("arm_is_combination"))
            row = {
                "pmid": item["pmid"], "payload": item["payload"], "source": "fulltext",
                "treated_median": treated, "control_median": control, "unit": unit,
                "ratio": None,
                "arm_is_combination": combination,
                "control_arm_kind": kind,
                "p2_intracranial": _affirmed(a["p2_intracranial"], b["p2_intracranial"]),
                "p4_payload_match": _affirmed(a["p4_payload_match"], b["p4_payload_match"]),
                "is_primary_research": bool(a.get("is_primary_research", True)
                                            and b.get("is_primary_research", True)),
                "extractors_agreed_on_medians": (
                    a["treated_median"] == b["treated_median"]
                    and a["control_median"] == b["control_median"]),
                "extractors_agreed_on_control_kind": kind_agreed,
                "survival_quote": a.get("survival_quote", "") if treated is not None else "",
                "control_arm_quote": a.get("control_arm_quote", ""),
                "a_note": a.get("note", ""), "b_note": b.get("note", ""),
            }
            if not row["extractors_agreed_on_medians"] or not kind_agreed:
                disagreements.append({
                    "pmid": item["pmid"], "payload": item["payload"],
                    "a": {k: a.get(k) for k in
                          ("treated_median", "control_median", "control_arm_kind")},
                    "b": {k: b.get(k) for k in
                          ("treated_median", "control_median", "control_arm_kind")},
                    "resolved_to": {"treated_median": treated, "control_median": control,
                                    "control_arm_kind": kind},
                })
        else:
            # No full text. The abstract-level reading stands, but an abstract that does not
            # describe the control arm is not evidence that the control was inert.
            row = {
                "pmid": item["pmid"], "payload": item["payload"], "source": "abstract_only",
                "treated_median": abstract.get("treated_median"),
                "control_median": abstract.get("control_median"),
                "unit": abstract.get("unit"),
                "ratio": abstract.get("ratio") if abstract.get("treated_median") is None
                else None,
                "arm_is_combination": False,
                "control_arm_kind": "uncertain",
                "p2_intracranial": _affirmed(abstract.get("p2_a"), abstract.get("p2_b")),
                "p4_payload_match": _affirmed(abstract.get("p4_a"), abstract.get("p4_b")),
                "is_primary_research": True,
                "extractors_agreed_on_medians": bool(abstract.get("extractors_agree")),
                "extractors_agreed_on_control_kind": False,
                "survival_quote": abstract.get("quote", ""),
                "control_arm_quote": "",
                "a_note": abstract.get("a_note", ""), "b_note": abstract.get("b_note", ""),
                "blocked_on": "full text not obtainable automatically; see "
                              "docs/t4t5_emax_search/manual_download_needed.md",
            }
        if not row["is_primary_research"]:
            row["p2_intracranial"] = "no"  # P1 failure: not a primary in vivo study
        rows.append(row)

    with OUT.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")

    full = [r for r in rows if r["source"] == "fulltext"]
    summary = {
        "registered_under": "docs/preregistration.md 2.2 item 11 (i.4) and (i.7)",
        "adjudicated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "adjudication_rules": [
            "a number survives only where both extractors read the same number",
            "the control arm's kind resolves to uncertain or active_therapeutic if either "
            "extractor says so; two different inert kinds do not change eligibility",
            "P2 and P4 are 'yes' only where both extractors affirm them",
            "a record without full text keeps its abstract reading and an uncertain control",
        ],
        "n_rows": len(rows),
        "n_from_fulltext": len(full),
        "n_abstract_only": len(rows) - len(full),
        "n_fulltext_median_disagreements": sum(
            1 for r in full if not r["extractors_agreed_on_medians"]),
        "n_fulltext_control_kind_disagreements": sum(
            1 for r in full if not r["extractors_agreed_on_control_kind"]),
        "control_arm_kinds": dict(Counter(r["control_arm_kind"] for r in rows)),
        "control_arm_kinds_fulltext_only": dict(
            Counter(r["control_arm_kind"] for r in full)),
        "disagreements": disagreements,
    }
    SUMMARY.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")

    print(f"{len(rows)} rows: {len(full)} from full text, "
          f"{len(rows) - len(full)} abstract-only")
    print(f"median disagreements: {summary['n_fulltext_median_disagreements']}, "
          f"control-kind disagreements: {summary['n_fulltext_control_kind_disagreements']}")
    print(f"control arms (full text): {summary['control_arm_kinds_fulltext_only']}")
    print(f"Wrote {OUT} and {SUMMARY}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("finalise",))
    parser.parse_args(argv)
    finalise()


if __name__ == "__main__":
    main()
