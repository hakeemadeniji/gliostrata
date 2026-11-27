"""Fetch CGGA's mRNA_array_301 expression release and freeze T1's subtype scoring set.

Pre-registration section 2.2 item 9 fixes T1's subtype scoring set by rule: the CGGA
microarray patients (batch mRNA_array_301) that carry a reference Verhaak label. Its
classifier is applied to the gene-level matrix by exact symbol (9(d)). Each step below is
declared to the read guard before any CGGA byte is touched.

``fetch``  (kind ``checksum``): downloads the public gene-level archive into
    ``data/raw/cgga/_cache/`` and records its size and SHA-256 in
    ``data/provenance/cgga_microarray_sources.json``. Nothing is unpacked or parsed beyond
    the zip signature. An existing archive is never overwritten.

``freeze`` (kind ``availability``): reads the matrix's header line (sample ids) and, from
    the cohort file, patient_id, cgga_dataset and whether the reference label is
    *present*. The label itself is never kept, compared or printed. The patients of
    mRNA_array_301 that have a label and an expression column are written as an
    identifiers-only table and sealed as the scoring-set manifest
    ``data/partitions/cgga_t1_subtype_v1.json`` (all-test). Nothing is overwritten.
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
FILE = "download/20200506/CGGA.mRNA_array_301_gene_level.20200506.txt.zip"
TYPE, TIME = "mRNA_array_301_gene_level", "20200506"
ARCHIVE = Path("data/raw/cgga/_cache") / Path(FILE).name
PROVENANCE = Path("data/provenance/cgga_microarray_sources.json")
COHORT = Path("data/raw/cgga/clinical_cohort.parquet")
TABLE = Path("data/processed/cgga/t1_subtype_scoring_set.parquet")
MANIFEST = Path("data/partitions/cgga_t1_subtype_v1.json")
LEDGER = Path("data/access_ledger.jsonl")
LABEL_COLUMN = "subtype"
BATCH = "mRNA_array_301"


def _url() -> str:
    return f"{PORTAL}?file={FILE}&type={TYPE}&time={TIME}"


def fetch(ledger: Path) -> None:
    if ARCHIVE.exists():
        raise SystemExit(f"{ARCHIVE} exists; never overwritten")
    entry = record_access(Access(
        cohort="cgga", kind=CHECKSUM, analysis="t1-microarray-download",
        purpose=("Download CGGA's public mRNA_array_301 gene-level expression release, the "
                 "input of T1's subtype row (pre-registration 2.2 item 9). Bytes are written "
                 "and hashed only; nothing is unpacked or parsed."),
        script="scripts/cgga_microarray.py",
    ), ledger_path=ledger)
    print(f"Access recorded: ledger entry {entry['seq']}.")
    resp = requests.get(_url(), timeout=900, allow_redirects=True)
    if resp.status_code != 200 or not resp.content.startswith(b"PK"):
        raise SystemExit(f"HTTP {resp.status_code}, not a zip archive")
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    ARCHIVE.write_bytes(resp.content)
    digest = hashlib.sha256(resp.content).hexdigest()
    PROVENANCE.write_text(json.dumps({
        "what": "CGGA mRNA_array_301 gene-level expression, fetched for T1's subtype row.",
        "retrieved_utc": entry["recorded_at"], "ledger_seq": entry["seq"],
        "file": {"path": ARCHIVE.as_posix(), "url": _url(), "bytes": len(resp.content),
                 "sha256": digest},
    }, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"  {len(resp.content):,} bytes, sha256 {digest}")


def _header(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as zf:
        members = [m for m in zf.namelist() if not m.endswith("/")
                   and not m.startswith("__MACOSX/") and not Path(m).name.startswith("._")]
        if len(members) != 1:
            raise SystemExit(f"{path}: expected one member, found {members}")
        with zf.open(members[0]) as handle:
            line = io.TextIOWrapper(handle, encoding="utf-8").readline()
    return [h.strip() for h in line.rstrip("\r\n").split("\t")[1:]]


def freeze(ledger: Path) -> None:
    from gliostrata.evaluation.partitions import freeze_partition, read_manifest, verify_partition

    if MANIFEST.exists() or TABLE.exists():
        raise SystemExit("scoring-set manifest or table exists; never overwritten")
    entry = record_access(Access(
        cohort="cgga", kind=AVAILABILITY, analysis="t1-subtype-scoring-set",
        purpose=("Freeze T1's subtype scoring set by its registered rule (pre-registration 2.2 "
                 "item 9): mRNA_array_301 patients with a reference label. Reads the expression "
                 "matrix's header line (sample ids) and, from the cohort file, patient_id, "
                 "cgga_dataset and the PRESENCE of the label (never its value, which is not kept "
                 "or printed). Writes an identifiers-only table and seals it all-test."),
        script="scripts/cgga_microarray.py",
        inputs=hash_inputs([ARCHIVE, COHORT]),
    ), ledger_path=ledger)
    print(f"Access recorded: ledger entry {entry['seq']}.")
    samples = set(_header(ARCHIVE))
    cohort = pd.read_parquet(COHORT, columns=["patient_id", "cgga_dataset", LABEL_COLUMN])
    batch = cohort[cohort["cgga_dataset"] == BATCH]
    has_label = batch[LABEL_COLUMN].notna().to_numpy()
    ids = batch["patient_id"].astype(str).to_numpy()
    del cohort, batch  # the label values go no further
    labelled = {i for i, h in zip(ids, has_label) if h}
    keep = sorted(labelled & samples)
    print(f"  {BATCH}: {len(ids)} rows, {len(set(ids))} ids; with a label: {len(labelled)}; "
          f"expression columns: {len(samples)}; labelled with expression: {len(keep)}")
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate ids within the batch; stop and review")
    table = pd.DataFrame({"patient_id": keep, "cgga_dataset": BATCH})
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(TABLE, index=False)
    _, manifest, written = freeze_partition(
        table, cohort_source=TABLE, manifest_path=MANIFEST, partition_name="cgga_t1_subtype_v1",
        seed=20251126, test_fraction=1.0, seal=True,
        notes=("T1 subtype-accuracy scoring set (pre-registration 2.2 item 9): every "
               f"{BATCH} patient with a reference Verhaak label and an expression column. "
               f"Written by scripts/cgga_microarray.py freeze; access ledger entry {entry['seq']}."))
    report = verify_partition(read_manifest(written), table, cohort_source=TABLE)
    if not report.ok:
        raise SystemExit(f"verification failed:\n{report.describe()}")
    print(f"Sealed {MANIFEST}: {len(manifest.test_ids)} patients; sha256 "
          f"{hashlib.sha256(Path(written).read_bytes()).hexdigest()}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("fetch", "freeze"))
    parser.add_argument("--ledger", default=str(LEDGER), type=Path)
    args = parser.parse_args(argv)
    (fetch if args.step == "fetch" else freeze)(args.ledger)


if __name__ == "__main__":
    main()
