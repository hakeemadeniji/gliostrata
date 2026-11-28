"""Resolve and fetch open-access full texts for the Emax search's stage 2.

Pre-registration item 11 (i.4) requires the extracted ratios to be committed to
`docs/t4t5_emax_search/`. 37 of the 40 records carried past stage 1 need full text, because
their abstracts do not give both arms' median survival, and item 11(i.7) additionally needs
each study's control arm classified. This script finds which of them are in PubMed Central's
open-access subset and downloads those.

``resolve``  elinks PubMed to PMC for every carried record and writes
    ``docs/t4t5_emax_search/fulltext_availability.json``: which records have a PMCID and
    which have none. No article text is fetched.

``fetch``    downloads the PMC XML for the resolved PMCIDs into
    ``data/external/t4t5_emax_fulltext/`` (git-ignored), one file per record, and records in
    the availability file whether a body section actually came back. A PMCID that exists but
    is not in the open-access subset returns metadata only; that is recorded as such rather
    than mistaken for a successful fetch.

Downloaded articles are untrusted data. They land in their own directory, are never executed,
and every parser is run with paths passed as arguments.

No cohort data is touched, so no read-guard access is recorded.
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

import requests

ELINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "gliostrata-prereg"  # NCBI asks for a tool name; the lead's e-mail stays private
SEARCH_DIR = Path("docs/t4t5_emax_search")
STAGE2 = SEARCH_DIR / "screening" / "stage2_input.jsonl"
AVAILABILITY = SEARCH_DIR / "fulltext_availability.json"
CACHE = Path("data/external/t4t5_emax_fulltext")
PAUSE_S = 0.4


def _carried() -> list[dict]:
    if not STAGE2.exists():
        raise SystemExit(f"{STAGE2} not found; run the screening stages first")
    return [json.loads(line) for line in STAGE2.read_text(encoding="utf-8").splitlines()
            if line.strip()]


def _unique_pmids(records: list[dict]) -> list[str]:
    seen: dict[str, None] = {}
    for rec in records:
        seen.setdefault(rec["pmid"], None)
    return list(seen)


def resolve() -> None:
    records = _carried()
    pmids = _unique_pmids(records)
    print(f"{len(records)} carried records, {len(pmids)} unique PMIDs")

    linked: dict[str, str | None] = {}
    for start in range(0, len(pmids), 100):
        chunk = pmids[start:start + 100]
        # One `id=` parameter per record. A comma-separated list makes E-utilities return a
        # single merged linkset for the whole chunk, which cannot be attributed per record.
        resp = requests.get(ELINK, timeout=120, params=[
            ("dbfrom", "pubmed"), ("db", "pmc"), ("retmode", "json"), ("tool", TOOL),
            *(("id", pmid) for pmid in chunk)])
        resp.raise_for_status()
        for linkset in resp.json().get("linksets", []):
            pmid = str(linkset.get("ids", [None])[0])
            pmcid = None
            for db in linkset.get("linksetdbs", []):
                if db.get("linkname") == "pubmed_pmc" and db.get("links"):
                    pmcid = f"PMC{db['links'][0]}"
                    break
            linked[pmid] = pmcid
        time.sleep(PAUSE_S)
    for pmid in pmids:
        linked.setdefault(pmid, None)

    with_pmc = sorted(p for p, v in linked.items() if v)
    without = sorted(p for p, v in linked.items() if not v)
    doc = {
        "registered_under": "docs/preregistration.md section 2.2 item 11 (i.4) and (i.7)",
        "resolved_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_unique_pmids": len(pmids),
        "n_with_pmcid": len(with_pmc),
        "n_without_pmcid": len(without),
        "pmcid_by_pmid": {p: linked[p] for p in pmids},
        "without_pmcid": without,
        "fetched": {},
        "note": ("A PMCID means the article is indexed in PMC, not that its full text is "
                 "open-access. `fetch` records which ones actually returned a body."),
    }
    AVAILABILITY.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"  with a PMCID: {len(with_pmc)}\n  without:      {len(without)}")
    print(f"Wrote {AVAILABILITY}")


def _body_chars(xml_bytes: bytes) -> int:
    """How much article body came back. Zero means metadata only, not open access."""
    try:
        root = ElementTree.fromstring(xml_bytes)
    except ElementTree.ParseError:
        return 0
    return sum(len("".join(body.itertext())) for body in root.iter("body"))


def fetch() -> None:
    if not AVAILABILITY.exists():
        raise SystemExit(f"{AVAILABILITY} not found; run `resolve` first")
    doc = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    CACHE.mkdir(parents=True, exist_ok=True)

    targets = {p: c for p, c in doc["pmcid_by_pmid"].items() if c}
    print(f"fetching {len(targets)} PMC records")
    results: dict[str, dict] = dict(doc.get("fetched", {}))
    for pmid, pmcid in sorted(targets.items()):
        out = CACHE / f"{pmid}_{pmcid}.xml"
        if out.exists():
            results[pmid] = {**results.get(pmid, {}), "pmcid": pmcid,
                             "file": out.as_posix(), "cached": True}
            continue
        resp = requests.get(EFETCH, timeout=180, params={
            "db": "pmc", "id": pmcid.removeprefix("PMC"), "retmode": "xml", "tool": TOOL})
        time.sleep(PAUSE_S)
        if resp.status_code != 200:
            results[pmid] = {"pmcid": pmcid, "http": resp.status_code, "body_chars": 0,
                             "open_access": False}
            print(f"  {pmid} {pmcid}: HTTP {resp.status_code}")
            continue
        chars = _body_chars(resp.content)
        if chars > 0:
            out.write_bytes(resp.content)
        results[pmid] = {
            "pmcid": pmcid, "http": 200, "body_chars": chars,
            "open_access": chars > 0,
            "file": out.as_posix() if chars > 0 else None,
        }
        print(f"  {pmid} {pmcid}: {'body ' + str(chars) + ' chars' if chars else 'METADATA ONLY'}")

    doc["fetched"] = results
    doc["fetched_at_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    doc["n_open_access"] = sum(1 for v in results.values() if v.get("open_access"))
    doc["n_metadata_only"] = sum(1 for v in results.values() if not v.get("open_access"))
    AVAILABILITY.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"\nopen access: {doc['n_open_access']}, metadata only: {doc['n_metadata_only']}, "
          f"no PMCID: {doc['n_without_pmcid']}")
    print(f"Articles cached under {CACHE} (git-ignored). Updated {AVAILABILITY}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("resolve", "fetch"))
    args = parser.parse_args(argv)
    (resolve if args.step == "resolve" else fetch)()


if __name__ == "__main__":
    main()
