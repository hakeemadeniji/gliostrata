"""Validate and merge the manually extracted Emax rows.

Pre-registration item 11(i.8) closed the Emax evidence base over the 108 search records, the
18 open-access full texts, and the full texts the project lead retrieved manually. The
retrieved PDFs are copyrighted and stay outside the repository, so what is committed is the
extraction from them: `docs/t4t5_emax_search/manual_extraction/manual_extraction.csv`, filled
in by hand against that folder's README.

``validate``  checks the filled CSV without changing it, and reports every problem rather
    than stopping at the first: unknown pmid/payload pairs, a median with no unit, a number
    with no supporting quote, a control-arm kind outside the registered five, a ratio outside
    a plausible range, and rows left blank. It exits non-zero if anything is wrong, so a bad
    sheet cannot be merged by accident.

``merge``     rewrites the `abstract_only` rows of `extraction_final.jsonl` from the sheet,
    marking them `source: "fulltext_manual"`. It refuses to run unless `validate` passes. It
    never touches a row that came from the open-access full texts, and it never adds a pmid
    that item 11(i.8) did not close over.

A blank row is kept as `abstract_only` with its control arm `uncertain`, which keeps it out of
the (i.7) primary set. That is the registered handling of a study whose control arm is not
described, and it is not an error.

No cohort data is touched, so no read-guard access is recorded.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SEARCH = Path("docs/t4t5_emax_search")
CSV_PATH = SEARCH / "manual_extraction" / "manual_extraction.csv"
FINAL = SEARCH / "screening" / "extraction_final.jsonl"
REPORT = SEARCH / "manual_extraction" / "validation_report.json"
CONTROL_KINDS = ("untreated", "vehicle", "mock_or_reporter_vector",
                 "active_therapeutic", "uncertain")
TERNARY = ("yes", "no", "uncertain")
UNITS = ("days", "weeks", "months")
# A treated-to-control median-survival ratio outside this range in a rodent glioma study is
# far more likely a transcription slip (a unit mixed up, two arms swapped) than a real effect.
# Flagged for the lead to confirm, never silently corrected or dropped.
RATIO_SANE = (0.5, 12.0)


def _rows() -> list[dict]:
    if not CSV_PATH.exists():
        raise SystemExit(f"{CSV_PATH} not found")
    with CSV_PATH.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _num(value: str | None) -> float | None:
    text = (value or "").strip()
    if not text:
        return None
    try:
        return float(text.replace(",", ""))
    except ValueError:
        return None


def _flag(value: str | None) -> bool:
    return (value or "").strip().upper() in {"TRUE", "T", "YES", "Y", "1"}


def validate(verbose: bool = True) -> dict:
    rows = _rows()
    known = {(r["pmid"], r["payload"]) for r in
             (json.loads(line) for line in FINAL.read_text(encoding="utf-8").splitlines()
              if line.strip())}

    problems: list[dict] = []
    filled = blank = 0

    def bad(row, what):
        problems.append({"pmid": row.get("pmid"), "payload": row.get("payload"),
                         "problem": what})

    for row in rows:
        key = (row.get("pmid", "").strip(), row.get("payload", "").strip())
        if key not in known:
            bad(row, f"pmid/payload {key} is not one item 11(i.8) closed over")
            continue
        treated, control = _num(row.get("treated_median")), _num(row.get("control_median"))
        unit = (row.get("unit") or "").strip().lower()
        kind = (row.get("control_arm_kind") or "").strip()

        if treated is None and control is None and not kind:
            blank += 1
            continue
        filled += 1

        for name, raw, parsed in (("treated_median", row.get("treated_median"), treated),
                                  ("control_median", row.get("control_median"), control)):
            if (raw or "").strip() and parsed is None:
                bad(row, f"{name} {raw!r} is not a number")
            if parsed is not None and parsed <= 0:
                bad(row, f"{name} must be positive, got {parsed}")

        if (treated is not None or control is not None) and unit not in UNITS:
            bad(row, f"unit must be one of {UNITS}, got {unit!r}")
        if (treated is not None or control is not None) \
                and not (row.get("survival_quote") or "").strip():
            bad(row, "a median was recorded with no survival_quote")
        if treated is not None and control:
            ratio = treated / control
            if not RATIO_SANE[0] <= ratio <= RATIO_SANE[1]:
                bad(row, f"ratio {ratio:.2f} is outside {RATIO_SANE}; check the unit and "
                         "which arm is which")
        if kind and kind not in CONTROL_KINDS:
            bad(row, f"control_arm_kind must be one of {CONTROL_KINDS}, got {kind!r}")
        if kind in CONTROL_KINDS[:4] and not (row.get("control_arm_quote") or "").strip():
            bad(row, f"control_arm_kind {kind!r} needs a control_arm_quote; "
                     "'uncertain' is the value that needs no quote")
        for field in ("p2_intracranial", "p4_payload_match"):
            value = (row.get(field) or "").strip().lower()
            if value and value not in TERNARY:
                bad(row, f"{field} must be one of {TERNARY}, got {value!r}")
        if _flag(row.get("arm_is_combination")) and not (row.get("combined_with") or "").strip():
            bad(row, "arm_is_combination is TRUE but combined_with is empty")

    report = {
        "checked_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "csv": CSV_PATH.as_posix(),
        "n_rows": len(rows), "n_filled": filled, "n_blank": blank,
        "n_problems": len(problems), "problems": problems,
        "ok": not problems,
        "control_arm_kinds": dict(Counter(
            (r.get("control_arm_kind") or "").strip() for r in rows
            if (r.get("control_arm_kind") or "").strip())),
    }
    REPORT.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8", newline="\n")
    if verbose:
        print(f"{len(rows)} rows: {filled} filled, {blank} blank")
        print(f"control arms: {report['control_arm_kinds']}")
        for problem in problems:
            print(f"  PROBLEM {problem['pmid']} [{problem['payload']}]: {problem['problem']}")
        print(f"\n{'OK' if report['ok'] else str(len(problems)) + ' PROBLEM(S)'} "
              f"— wrote {REPORT}")
    return report


def merge() -> None:
    report = validate(verbose=False)
    if not report["ok"]:
        raise SystemExit(
            f"{report['n_problems']} problem(s) in {CSV_PATH}; run `validate` and fix them. "
            "Nothing was merged.")
    sheet = {(r["pmid"].strip(), r["payload"].strip()): r for r in _rows()}
    out, updated = [], 0
    for line in FINAL.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        key = (row["pmid"], row["payload"])
        entry = sheet.get(key)
        if row["source"] != "abstract_only" or entry is None:
            out.append(row)
            continue
        treated, control = _num(entry.get("treated_median")), _num(entry.get("control_median"))
        kind = (entry.get("control_arm_kind") or "").strip()
        if treated is None and control is None and not kind:
            out.append(row)  # left blank: stays abstract_only with an uncertain control
            continue
        row.update({
            "source": "fulltext_manual",
            "treated_median": treated, "control_median": control,
            "unit": (entry.get("unit") or "").strip().lower() or None,
            "ratio": None,
            "arm_is_combination": _flag(entry.get("arm_is_combination")),
            "control_arm_kind": kind or "uncertain",
            "p2_intracranial": (entry.get("p2_intracranial") or "uncertain").strip().lower(),
            "p4_payload_match": (entry.get("p4_payload_match") or "uncertain").strip().lower(),
            "is_primary_research": (
                (entry.get("is_primary_research") or "TRUE").strip().upper()
                not in {"FALSE", "F", "NO", "N", "0"}),
            "survival_quote": (entry.get("survival_quote") or "").strip(),
            "control_arm_quote": (entry.get("control_arm_quote") or "").strip(),
            "extracted_by": "project lead, from the manually retrieved full text",
            "pdf_filename": (entry.get("pdf_filename") or "").strip(),
            "manual_note": (entry.get("note") or "").strip(),
        })
        row.pop("blocked_on", None)
        if not row["is_primary_research"]:
            row["p2_intracranial"] = "no"  # P1 failure
        out.append(row)
        updated += 1

    with FINAL.open("w", encoding="utf-8", newline="\n") as handle:
        for row in out:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    print(f"merged {updated} manually extracted rows into {FINAL}")
    print(dict(Counter(r["source"] for r in out)))
    print("Now run: python -I scripts/t4t5_emax_ratios.py compute")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("validate", "merge"))
    args = parser.parse_args(argv)
    if args.step == "validate":
        raise SystemExit(0 if validate()["ok"] else 1)
    merge()


if __name__ == "__main__":
    main()
