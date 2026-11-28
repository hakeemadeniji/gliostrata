"""Run the registered Emax literature search (pre-registration 2.2 item 11 (i.4)).

Item 11(i) fixes Emax_p = 0.5 * ln(ratio_p), with ratio_p the median of the
treated-to-control median-survival ratios across every study qualifying under P1 to P5.
The three queries are written out in (i.4), run **once** each, with no year limit. This
script runs them and writes an append-only log; it refuses to run a second time, because
re-running with a different query is exactly what the registration forbids.

``search``  runs the three registered queries against PubMed esearch, fetches each record's
    title, abstract, journal, year and DOI, and writes ``search.json`` (the queries, PubMed's
    own translation of each, the counts and the PMID lists) and ``records.jsonl`` (one record
    per line, per payload). Nothing is screened here: screening is a separate, dual,
    independent step whose decisions go in ``screening/``.

No cohort data is touched, so no read-guard access is recorded. This reads PubMed only.
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

import requests

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# NCBI asks automated clients to identify the tool. It also asks for an e-mail address;
# the project lead's address is private and is not sent (see .clinerules).
TOOL = "gliostrata-prereg"
OUT = Path("docs/t4t5_emax_search")
RETMAX = 2000
PAUSE_S = 0.4  # under NCBI's 3 requests/second limit without a key

# The registered queries, verbatim from pre-registration 2.2 item 11 (i.4).
QUERIES: dict[str, str] = {
    "CRAd": (
        '("oncolytic adenovirus"[tiab] OR CRAd[tiab] OR "conditionally replicative '
        'adenovirus"[tiab]) AND (glioma[tiab] OR glioblastoma[tiab]) AND "median survival"[tiab]'
    ),
    "prodrug_enzyme": (
        '("cytosine deaminase"[tiab] OR "thymidine kinase"[tiab] OR HSV-TK[tiab] OR '
        'carboxylesterase[tiab] OR prodrug[tiab]) AND (glioma[tiab] OR glioblastoma[tiab]) '
        'AND "median survival"[tiab]'
    ),
    "CAR": (
        '("chimeric antigen receptor"[tiab] OR CAR-T[tiab] OR CAR-NK[tiab]) AND '
        '(glioma[tiab] OR glioblastoma[tiab]) AND "median survival"[tiab]'
    ),
}


def _esearch(query: str) -> tuple[list[str], str, int]:
    resp = requests.get(ESEARCH, timeout=120, params={
        "db": "pubmed", "term": query, "retmax": RETMAX, "retmode": "json", "tool": TOOL})
    resp.raise_for_status()
    result = resp.json()["esearchresult"]
    if "ERROR" in result:
        raise SystemExit(f"PubMed refused the query: {result['ERROR']}")
    return result["idlist"], result.get("querytranslation", ""), int(result["count"])


def _text(node, path: str) -> str:
    found = node.find(path)
    return "".join(found.itertext()).strip() if found is not None else ""


def _efetch(pmids: list[str]) -> list[dict]:
    records: list[dict] = []
    for start in range(0, len(pmids), 200):
        chunk = pmids[start:start + 200]
        resp = requests.post(EFETCH, timeout=300, data={
            "db": "pubmed", "id": ",".join(chunk), "retmode": "xml", "tool": TOOL})
        resp.raise_for_status()
        root = ElementTree.fromstring(resp.content)
        for art in root.iter("PubmedArticle"):
            abstract = " ".join(
                "".join(piece.itertext()).strip()
                for piece in art.iter("AbstractText")).strip()
            # The article's OWN doi, never one from its reference list. `art.iter` walks the
            # whole PubmedArticle, including <ReferenceList>, where every cited paper carries
            # its own <ArticleId IdType="doi">: PMID 37420122 has 40 of them, 39 of which are
            # references. Taking the last match gave that record a 2009 DOI for a 2023 paper.
            # <ELocationID EIdType="doi"> lives inside <Article> and cannot be a reference;
            # the ArticleIdList fallback is scoped to <PubmedData>, outside ReferenceList.
            doi = ""
            eloc = art.find(".//Article/ELocationID[@EIdType='doi']")
            if eloc is not None and eloc.text:
                doi = eloc.text.strip()
            else:
                ident = art.find("./PubmedData/ArticleIdList/ArticleId[@IdType='doi']")
                if ident is not None and ident.text:
                    doi = ident.text.strip()
            records.append({
                "pmid": _text(art, ".//PMID"),
                "title": _text(art, ".//ArticleTitle"),
                "abstract": abstract,
                "journal": _text(art, ".//Journal/ISOAbbreviation"),
                "year": _text(art, ".//JournalIssue/PubDate/Year") or _text(
                    art, ".//JournalIssue/PubDate/MedlineDate"),
                "doi": doi,
                "publication_types": sorted({
                    "".join(p.itertext()).strip() for p in art.iter("PublicationType")}),
            })
        time.sleep(PAUSE_S)
    return records


def search() -> None:
    log, records_path = OUT / "search.json", OUT / "records.jsonl"
    if log.exists() or records_path.exists():
        raise SystemExit(
            f"{log} or {records_path} already exists. The registered search runs ONCE "
            "(pre-registration 2.2 item 11 (i.4)); it is never re-run or overwritten.")
    OUT.mkdir(parents=True, exist_ok=True)
    run_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    per_payload, all_records = {}, []
    for payload, query in QUERIES.items():
        pmids, translation, count = _esearch(query)
        if count > RETMAX:
            raise SystemExit(
                f"{payload}: {count} records exceeds retmax {RETMAX}. Stop: narrowing the "
                "query now would be a post-hoc change the registration forbids. Raise "
                "RETMAX and page instead, which changes no criterion.")
        print(f"{payload}: {count} records")
        time.sleep(PAUSE_S)
        fetched = _efetch(pmids)
        if len(fetched) != len(pmids):
            print(f"  warning: fetched {len(fetched)} of {len(pmids)} records")
        for rec in fetched:
            rec["payload"] = payload
        all_records.extend(fetched)
        per_payload[payload] = {
            "query": query, "query_translation": translation,
            "count": count, "n_fetched": len(fetched), "pmids": pmids,
        }

    overlap = {}
    for payload, info in per_payload.items():
        for other, other_info in per_payload.items():
            if payload < other:
                shared = sorted(set(info["pmids"]) & set(other_info["pmids"]))
                if shared:
                    overlap[f"{payload}|{other}"] = shared

    log.write_text(json.dumps({
        "registered_under": "docs/preregistration.md section 2.2 item 11 (i.1) to (i.6)",
        "database": "pubmed",
        "run_at_utc": run_at,
        "run_once": True,
        "year_limit": None,
        "payloads": per_payload,
        "cross_payload_overlap": overlap,
        "n_records_total": len(all_records),
        "n_unique_pmids": len({r["pmid"] for r in all_records}),
        "screening": "not started; dual independent, see screening/ and README.md",
    }, indent=1) + "\n", encoding="utf-8", newline="\n")
    with records_path.open("w", encoding="utf-8", newline="\n") as handle:
        for rec in all_records:
            handle.write(json.dumps(rec, ensure_ascii=False) + "\n")
    print(f"\nWrote {log} and {records_path}: {len(all_records)} records, "
          f"{len({r['pmid'] for r in all_records})} unique PMIDs.")
    if overlap:
        print(f"Cross-payload overlap: { {k: len(v) for k, v in overlap.items()} }")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("search", "fix-dois"))
    args = parser.parse_args(argv)
    (search if args.step == "search" else fix_dois)()




def fix_dois(ledger_unused=None) -> None:
    """Re-fetch the `doi` field of the committed records and correct it in place.

    The first run took each record's doi from the LAST `<ArticleId IdType="doi">` in its
    PubmedArticle element, which for any record carrying a reference list is a cited paper's
    doi rather than its own. Nothing else in the record is touched: the search, its queries,
    its counts, its PMIDs and every screening decision stand, because none of them used the
    doi field. Only the identifier printed for a human to look the paper up was wrong.
    """
    path = OUT / "records.jsonl"
    records = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
               if line.strip()]
    pmids = sorted({r["pmid"] for r in records})
    print(f"re-fetching dois for {len(pmids)} unique PMIDs")
    corrected = {r["pmid"]: r["doi"] for r in _efetch(pmids)}

    changed = []
    for rec in records:
        was, now = rec["doi"], corrected.get(rec["pmid"], rec["doi"])
        if was != now:
            changed.append({"pmid": rec["pmid"], "was": was, "now": now})
            rec["doi"] = now
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for rec in records:
            handle.write(json.dumps(rec, ensure_ascii=False) + "\n")

    report = OUT / "doi_correction.json"
    report.write_text(json.dumps({
        "corrected_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "cause": ("the doi was read from the last <ArticleId IdType='doi'> in the "
                  "PubmedArticle element, which for a record with a reference list is a "
                  "cited paper's doi, not the article's own"),
        "scope": ("the doi field only; queries, counts, PMIDs and every screening and "
                  "extraction decision are unaffected, none having used the doi"),
        "n_records": len(records), "n_changed": len(changed), "changes": changed,
    }, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"corrected {len(changed)} of {len(records)} records; wrote {report}")


if __name__ == "__main__":
    main()
