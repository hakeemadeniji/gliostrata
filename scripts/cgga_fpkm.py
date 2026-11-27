"""Fetch CGGA's published RSEM FPKM release and check it covers T2's comparator.

Pre-registration section 2.2 item 2: the comparator chosen by the registered search,
PMID 38174685, expects log2(FPKM + 1). Its authors validated it on CGGA's own FPKM
release, so that release is fetched, not computed from our counts.

Two steps, each declared to the read guard before any CGGA byte is touched:

``fetch``  (kind ``checksum``): downloads the two public archives into
    ``data/raw/cgga/_cache/`` and records their sizes and SHA-256 in
    ``data/provenance/cgga_fpkm_sources.json``. Nothing is unpacked or parsed beyond
    the two-byte zip signature. An existing archive is never overwritten.

``check``  (kind ``availability``): reads, from each archive, the header line (sample
    ids) and the first column (gene symbols) only. It reports which of the comparator's
    genes are present and how many ``cgga_external_v1`` patients have a column. No
    expression value is kept or printed, and no clinical field or outcome is opened.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pandas as pd
import requests

from gliostrata.evaluation.read_guard import (
    AVAILABILITY,
    CHECKSUM,
    Access,
    hash_inputs,
    record_access,
)

PORTAL = "https://www.cgga.org.cn/download"
CACHE = Path("data/raw/cgga/_cache")
PROVENANCE = Path("data/provenance/cgga_fpkm_sources.json")
PARTITION = Path("data/partitions/cgga_external_v1.json")
LEDGER = Path("data/access_ledger.jsonl")
FILES = {
    "mRNAseq_693": "download/20200506/CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip",
    "mRNAseq_325": "download/20200506/CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip",
}
TIME = "20200506"
#: PMID 38174685 (Sun and Liu, Folia Neuropathol 2024;62:59-75), as printed.
COMPARATOR_GENES = ("CPPED1", "F3", "G0S2", "MAPK1", "MME", "MMP9", "MPO")


def _url(batch: str) -> str:
    return f"{PORTAL}?file={FILES[batch]}&type={batch}&time={TIME}"


def _path(batch: str) -> Path:
    return CACHE / Path(FILES[batch]).name


def fetch(ledger: Path) -> None:
    targets = [_path(b) for b in FILES]
    existing = [p for p in targets if p.exists()]
    if existing:
        raise SystemExit(f"already present, never overwritten: {existing}")
    entry = record_access(Access(
        cohort="cgga", kind=CHECKSUM, analysis="t2-comparator-fpkm-download",
        purpose=("Download CGGA's public RSEM FPKM release (mRNAseq_693, mRNAseq_325), "
                 "the expression scale the T2 comparator (PMID 38174685) expects. Bytes are "
                 "written and hashed only; nothing is unpacked or parsed."),
        script="scripts/cgga_fpkm.py",
    ), ledger_path=ledger)
    print(f"Access recorded: ledger entry {entry['seq']}.")
    record = {"what": "CGGA RSEM FPKM release, fetched for T2's comparator (pre-registration 2.2 item 2).",
              "retrieved_utc": entry["recorded_at"], "ledger_seq": entry["seq"], "files": {}}
    for batch in FILES:
        resp = requests.get(_url(batch), timeout=900, allow_redirects=True)
        if resp.status_code != 200 or not resp.content.startswith(b"PK"):
            raise SystemExit(f"{batch}: HTTP {resp.status_code}, not a zip archive")
        path = _path(batch)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(resp.content)
        digest = hashlib.sha256(resp.content).hexdigest()
        record["files"][batch] = {"path": path.as_posix(), "url": _url(batch),
                                  "bytes": len(resp.content), "sha256": digest}
        print(f"  {batch}: {len(resp.content):,} bytes, sha256 {digest}")
    PROVENANCE.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"Provenance written to {PROVENANCE}.")


def _names(path: Path) -> tuple[list[str], list[str]]:
    """(gene symbols, sample ids) of the single matrix in one archive."""
    with zipfile.ZipFile(path) as zf:
        # macOS archive metadata (__MACOSX/, ._*) is not data.
        members = [m for m in zf.namelist() if not m.endswith("/")
                   and not m.startswith("__MACOSX/") and not Path(m).name.startswith("._")]
        if len(members) != 1:
            raise SystemExit(f"{path}: expected one member, found {members}")
        with zf.open(members[0]) as handle:
            text = io.TextIOWrapper(handle, encoding="utf-8")
            header = text.readline().rstrip("\n\r").split("\t")
        with zf.open(members[0]) as handle:
            genes = pd.read_csv(handle, sep="\t", usecols=[0], dtype=str).iloc[:, 0]
    return [str(g).strip() for g in genes], [h.strip() for h in header[1:]]


def check(ledger: Path) -> None:
    paths = [_path(b) for b in FILES]
    entry = record_access(Access(
        cohort="cgga", kind=AVAILABILITY, analysis="t2-comparator-fpkm-availability",
        purpose=("Rule E4 for the T2 comparator (PMID 38174685): which of its 7 genes are rows "
                 "of CGGA's FPKM matrices, and how many cgga_external_v1 patients are columns. "
                 "Reads each matrix's header line and first column (gene symbols) only; no "
                 "expression value is kept or printed; no clinical field or outcome."),
        script="scripts/cgga_fpkm.py",
        inputs=hash_inputs(paths + [PARTITION]),
    ), ledger_path=ledger)
    print(f"Access recorded: ledger entry {entry['seq']}.")
    ids = set(json.loads(PARTITION.read_text(encoding="utf-8"))["test_ids"])
    covered: set[str] = set()
    for batch, path in zip(FILES, paths):
        genes, samples = _names(path)
        present = {g: (g in set(genes)) for g in COMPARATOR_GENES}
        dup = sorted({g for g in COMPARATOR_GENES if genes.count(g) > 1})
        hit = ids & set(samples)
        covered |= hit
        print(f"  {batch}: {len(genes):,} genes, {len(samples):,} samples; "
              f"comparator genes present: {sum(present.values())}/7 "
              f"{[g for g, ok in present.items() if not ok] or ''}; duplicated: {dup or 'none'}; "
              f"partition patients with a column: {len(hit)}")
    print(f"Partition patients with an FPKM column in either batch: {len(covered)} of {len(ids)}.")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("fetch", "check"))
    parser.add_argument("--ledger", default=str(LEDGER), type=Path)
    args = parser.parse_args(argv)
    (fetch if args.step == "fetch" else check)(args.ledger)


if __name__ == "__main__":
    main()
