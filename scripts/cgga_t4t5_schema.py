"""Read the column names T4/T5's loader must be written against. Names only, never values.

The T4/T5 evaluation (`scripts/t4t5_evaluate.py`) reads CGGA expression and clinical fields
for the 218 patients of the sealed external partition. That read happens **once**. A loader
written against guessed column names which then fails on the day would cost the registered
read, or force a disclosed re-run, so the names are read here in advance.

This records a ``schema`` access, which the read guard does not gate: only ``values`` is
gated, and only ``values`` needs the OSF DOI. Ledger entries 11, 16 and 17 are the same
shape — availability and schema accesses that read names, never values.

**What is read:** the header line of each expression archive, and the column names of the
clinical table. **What is not read:** any expression value, any clinical field's contents,
any outcome, any identifier. Nothing is printed or written but names, dtypes and counts.
"""

from __future__ import annotations

import argparse
import io
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from gliostrata.evaluation.read_guard import SCHEMA, Access, hash_inputs, record_access

RAW = Path("data/raw/cgga")
CLINICAL = RAW / "clinical_cohort.parquet"
FPKM_325 = RAW / "_cache" / "CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip"
FPKM_693 = RAW / "_cache" / "CGGA.mRNAseq_693.RSEM-genes.20200506.txt.zip"
OUT = Path("docs/cgga_t4t5_schema.json")
LEDGER = Path("data/access_ledger.jsonl")
MAX_SAMPLE_NAMES = 3   # enough to show the id format; never the full roster


def _header(path: Path) -> dict:
    """The first line of the single member of a zip archive. No data row is read."""
    with zipfile.ZipFile(path) as zf:
        members = [m for m in zf.namelist() if not m.endswith("/")
                   and not m.startswith("__MACOSX/") and not Path(m).name.startswith("._")]
        if len(members) != 1:
            raise SystemExit(f"{path}: expected one member, found {members}")
        with zf.open(members[0]) as handle:
            line = io.TextIOWrapper(handle, encoding="utf-8").readline()
    fields = [h.strip() for h in line.rstrip("\r\n").split("\t")]
    return {
        "member": members[0],
        "first_column_name": fields[0],
        "n_columns": len(fields),
        "n_sample_columns": len(fields) - 1,
        "sample_id_examples": fields[1:1 + MAX_SAMPLE_NAMES],
        "note": "header line only; no data row was read",
    }


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ledger", type=Path, default=LEDGER)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args(argv)

    for path in (CLINICAL, FPKM_325, FPKM_693):
        if not path.exists():
            raise SystemExit(f"{path} not found")

    entry = record_access(Access(
        cohort="cgga", kind=SCHEMA, analysis="t4t5-schema",
        purpose=("Read the column names the T4/T5 loader must be written against "
                 "(pre-registration 2.2 item 11): the header line of each RSEM FPKM archive "
                 "and the column names and dtypes of the clinical table. The T4/T5 read "
                 "happens once, and a loader written against guessed names would cost it. "
                 "NO value, outcome or identifier is read: names, dtypes and counts only."),
        script="scripts/cgga_t4t5_schema.py",
        inputs=hash_inputs([CLINICAL, FPKM_325, FPKM_693]),
    ), ledger_path=args.ledger)
    print(f"Access recorded: ledger entry {entry['seq']} (kind=schema, not gated).")

    # Column names and dtypes only. `read_parquet` with columns=[] reads no column data.
    clinical_schema = pd.read_parquet(CLINICAL, columns=[]).dtypes
    empty = pd.read_parquet(CLINICAL).head(0)   # zero rows: dtypes without any value
    document = {
        "what": "CGGA column names for the T4/T5 loader. Names, dtypes and counts only.",
        "registered_under": "docs/preregistration.md 2.2 item 11",
        "read_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "ledger_seq": entry["seq"],
        "access_kind": "schema (ungated; only `values` needs the OSF DOI)",
        "clinical": {
            "path": CLINICAL.as_posix(),
            "n_columns": int(empty.shape[1]),
            "columns": {str(k): str(v) for k, v in empty.dtypes.items()},
        },
        "expression": {
            "mRNAseq_325": _header(FPKM_325),
            "mRNAseq_693": _header(FPKM_693),
        },
        "not_read": ("no expression value, no clinical field contents, no outcome, no "
                     "patient identifier beyond the sample-id format shown"),
    }
    del clinical_schema
    args.out.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")

    print(f"\nclinical: {document['clinical']['n_columns']} columns")
    for name, dtype in document["clinical"]["columns"].items():
        print(f"   {name:<28} {dtype}")
    for batch, info in document["expression"].items():
        print(f"\n{batch}: first column {info['first_column_name']!r}, "
              f"{info['n_sample_columns']} sample columns, "
              f"e.g. {info['sample_id_examples']}")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
