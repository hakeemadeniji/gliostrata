"""Pull the survival and control-arm passages out of the fetched PMC full texts.

Stage 2 of the Emax search needs two things from each full text: the treated and control
median survivals (P3, pre-registration item 11 (i.3)) and what the control arm received
(item 11(i.7)). Both live in a handful of paragraphs. This script extracts those paragraphs
so that the extraction step reads a compact, committed excerpt rather than 600,000 characters
of article XML, and so that what was read is itself auditable.

A paragraph is kept when it mentions median survival, or when it names a control arm, or when
a figure or table caption does either. Section titles are kept with the paragraphs so that a
methods statement can be told from a result.

The articles under `data/external/t4t5_emax_fulltext/` are **untrusted downloaded data**. They
are parsed, never executed, their paths come in as arguments, and the text that comes out is
treated as data by every later step. Run this with `python -I`.

No cohort data is touched, so no read-guard access is recorded.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from xml.etree import ElementTree

CACHE = Path("data/external/t4t5_emax_fulltext")
SEARCH_DIR = Path("docs/t4t5_emax_search")
AVAILABILITY = SEARCH_DIR / "fulltext_availability.json"
OUT = SEARCH_DIR / "screening" / "fulltext_excerpts.jsonl"

SURVIVAL = re.compile(r"median\s+(?:overall\s+)?survival|survival\s+(?:time|was|of)\b"
                      r"|MS\s+of\s+\d|days\s+post|log[- ]rank", re.I)
CONTROL = re.compile(r"\b(?:control|vehicle|saline|PBS|mock|untreated|naive|parental|"
                     r"empty\s+vector|GFP|sham|no\s+treatment|phosphate[- ]buffered)\b", re.I)
NUMBER = re.compile(r"\d")
MAX_PARAGRAPHS = 40
MAX_CHARS = 1200


def _flatten(node) -> str:
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip()


def _section_title(path: list) -> str:
    titles = []
    for node in path:
        title = node.find("title")
        if title is not None:
            titles.append(_flatten(title))
    return " > ".join(t for t in titles if t)


def _walk(node, ancestors: list, out: list) -> None:
    """Collect paragraphs and captions, each with the section titles above it."""
    for child in node:
        if child.tag in {"p", "caption", "title"} and child.tag != "title":
            text = _flatten(child)
            if text:
                out.append((_section_title(ancestors), child.tag, text))
        elif child.tag in {"sec", "body", "fig", "table-wrap", "abstract", "floats-group"}:
            _walk(child, [*ancestors, child], out)
        else:
            _walk(child, ancestors, out)


def excerpt_one(path: Path) -> dict:
    root = ElementTree.fromstring(path.read_bytes())
    collected: list[tuple[str, str, str]] = []
    for body in list(root.iter("body")) + list(root.iter("floats-group")):
        _walk(body, [body], collected)

    kept = []
    for section, tag, text in collected:
        hits_survival = bool(SURVIVAL.search(text)) and bool(NUMBER.search(text))
        hits_control = bool(CONTROL.search(text))
        if hits_survival or (hits_control and NUMBER.search(text)):
            kept.append({
                "section": section, "kind": tag,
                "text": text[:MAX_CHARS],
                "mentions_survival": hits_survival,
                "mentions_control": hits_control,
            })
    # Keep the survival paragraphs first: they are the ones P3 turns on.
    kept.sort(key=lambda k: (not k["mentions_survival"], len(k["text"])))
    return {
        "file": path.as_posix(),
        "n_paragraphs_scanned": len(collected),
        "n_paragraphs_kept": len(kept),
        "truncated": len(kept) > MAX_PARAGRAPHS,
        "paragraphs": kept[:MAX_PARAGRAPHS],
    }


def build() -> None:
    if not AVAILABILITY.exists():
        raise SystemExit(f"{AVAILABILITY} not found; run t4t5_emax_fulltext.py first")
    doc = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    stage2 = [json.loads(line) for line
              in (SEARCH_DIR / "screening" / "stage2_input.jsonl").read_text(
                  encoding="utf-8").splitlines() if line.strip()]
    payloads: dict[str, list[str]] = {}
    for rec in stage2:
        payloads.setdefault(rec["pmid"], []).append(rec["payload"])

    rows = []
    for pmid, info in sorted(doc.get("fetched", {}).items()):
        if not info.get("open_access") or not info.get("file"):
            continue
        path = Path(info["file"])
        if not path.exists():
            print(f"  {pmid}: {path} missing, skipped")
            continue
        try:
            found = excerpt_one(path)
        except ElementTree.ParseError as exc:
            print(f"  {pmid}: XML did not parse ({exc}), skipped")
            continue
        rows.append({"pmid": pmid, "pmcid": info["pmcid"],
                     "payloads": payloads.get(pmid, []), **found})
        print(f"  {pmid}: kept {found['n_paragraphs_kept']} of "
              f"{found['n_paragraphs_scanned']} paragraphs")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    total = sum(len(json.dumps(r)) for r in rows)
    print(f"\nWrote {OUT}: {len(rows)} articles, {total:,} characters of excerpt.")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("build",))
    parser.parse_args(argv)
    build()


if __name__ == "__main__":
    main()
