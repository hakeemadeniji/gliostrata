"""
List the CGGA patients eligible for the sealed external partition (T4 and T5).

Protocol tests T4 and T5 draw their patients from a sealed external partition
(docs/preregistration.md, section 2.2 item 7). The registered eligibility rule
is: a CGGA patient with an RNA-seq profile, a survival outcome and a usable
identifier mapping. This script applies that rule and nothing else.

What it does NOT do
-------------------
* It fits nothing and scores nothing.
* It prints no clinical value: only counts of rows and patients.
* It does not open the expression matrix. Expression availability is read from
  the COLUMN NAMES of counts.parquet (its Parquet schema), never its values.
* The table it writes carries identifiers and a batch label only -- no outcome,
  no covariate -- so the file that gets sealed cannot leak a value into
  anything that later reads it.

The counting rule is the ingest's own (``expression_overlap`` in
scripts/download_cgga.py), reused rather than re-implemented, so the two cannot
drift apart.

Read discipline
---------------
Every run, the dry run included, is an access to the held-out cohort. Before
any CGGA file is opened, the run declares itself to the read guard
(``gliostrata.evaluation.read_guard.record_access``) as an ``availability``
access, with the SHA-256 of both inputs, and is appended to the access ledger.
The guard admits an availability access before the deposit; it is recorded all
the same, because the number of looks at the external cohort is itself a
registration commitment. Only the four columns the rule needs are read from the
cohort file.

Usage
-----
    python scripts/cgga_external_eligibility.py --dry-run   # counts only
    python scripts/cgga_external_eligibility.py             # writes the table
    python scripts/cgga_external_eligibility.py --seal data/partitions/cgga_external_v1.json

``--seal`` writes the table, then freezes it as the all-test partition
``cgga_external_v1`` (every eligible patient held out; the CLI's defaults:
seed 20251126, test fraction 1.0) and verifies the manifest against the table,
all under the one ledger entry. It refuses to overwrite a manifest. The
equivalent separate step, which does not record an access, is

    gliostrata freeze-partition --input data/processed/cgga/cgga_external_eligible.parquet --manifest data/partitions/cgga_external_v1.json --name cgga_external_v1 --all-test
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

DEFAULT_COHORT = Path("data/raw/cgga/clinical_cohort.parquet")
DEFAULT_COUNTS = Path("data/raw/cgga/counts.parquet")
DEFAULT_OUT = Path("data/processed/cgga/cgga_external_eligible.parquet")
DEFAULT_LEDGER = Path("data/access_ledger.jsonl")
DEFAULT_PARTITION_NAME = "cgga_external_v1"
#: The freeze-partition CLI's default seed. With every patient held out the seed
#: changes no assignment; it is recorded in the manifest all the same.
DEFAULT_SEED = 20251126
ANALYSIS = "external-partition-eligibility"
#: The cohort columns the eligibility rule needs. Nothing else is read.
COHORT_COLUMNS = ("patient_id", "cgga_dataset", "overall_survival_months", "vital_status")

#: Columns of the table that gets sealed: identifiers and a batch label only.
OUTPUT_COLUMNS = ("patient_id", "cgga_dataset")

_INGEST = Path(__file__).resolve().parent / "download_cgga.py"


def load_ingest_module():
    """The CGGA ingest script, loaded by path (scripts/ is not a package)."""
    spec = importlib.util.spec_from_file_location(
        "gliostrata_cgga_download_script", _INGEST
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def expression_sample_ids(counts_path) -> frozenset[str]:
    """Sample ids that have an expression column, read from the Parquet SCHEMA.

    Only the column names are read (from the file's footer metadata). No
    expression value is loaded.
    """
    import pyarrow.parquet as pq

    return frozenset(str(c) for c in pq.ParquetFile(counts_path).schema_arrow.names)


def read_cohort_columns(cohort_path) -> pd.DataFrame:
    """Only the columns the rule needs (those present), never the whole file."""
    import pyarrow.parquet as pq

    present = set(pq.ParquetFile(cohort_path).schema_arrow.names)
    return pd.read_parquet(cohort_path, columns=[c for c in COHORT_COLUMNS if c in present])


def declare_access(cohort_path, counts_path, ledger_path, *, dry_run: bool,
                   seal: str | None) -> dict:
    """Record the access in the ledger through the read guard, before reading."""
    from gliostrata.evaluation.read_guard import (
        AVAILABILITY,
        Access,
        hash_inputs,
        record_access,
    )

    if dry_run:
        action = "Dry run: counts only, nothing written."
    elif seal:
        action = ("Write the identifiers-only eligibility table, then freeze it as the "
                  f"all-test partition at {Path(seal).as_posix()} and verify the manifest.")
    else:
        action = "Write the identifiers-only eligibility table; nothing is sealed."
    access = Access(
        cohort="cgga",
        kind=AVAILABILITY,
        analysis=ANALYSIS,
        purpose=(
            "Eligibility for the sealed external partition by the registered rule "
            "(pre-registration section 2.2 item 7): an RNA-seq profile, a survival "
            "outcome and a usable identifier. Read: patient_id, cgga_dataset, and "
            "overall_survival_months and vital_status for presence only, from the "
            "cohort file; the column names of counts.parquet (Parquet schema), never "
            f"an expression value. Counts printed, no value. {action} Nothing fitted "
            "or scored."
        ),
        script="scripts/cgga_external_eligibility.py",
        inputs=hash_inputs([cohort_path, counts_path]),
    )
    return record_access(access, ledger_path=ledger_path)


def seal_partition(table: pd.DataFrame, table_path: Path, manifest_path: Path,
                   name: str, *, ledger_seq: int) -> None:
    """Freeze the all-test partition from the written table and verify it."""
    from gliostrata.evaluation.partitions import (
        freeze_partition,
        read_manifest,
        verify_partition,
    )

    _, manifest, written = freeze_partition(
        table,
        cohort_source=table_path,
        manifest_path=manifest_path,
        partition_name=name,
        seed=DEFAULT_SEED,
        test_fraction=1.0,
        seal=True,
        notes=("CGGA external partition, all-test: every patient eligible under "
               "pre-registration section 2.2 item 7 is held out. Primary external test "
               "for T2 (already evaluated before this seal, section 2.3), T4 and T5; not "
               "the scoring set of either T1 row (section 4.5). Written by "
               "scripts/cgga_external_eligibility.py --seal; access ledger entry "
               f"{ledger_seq}."),
    )
    report = verify_partition(read_manifest(written), table, cohort_source=table_path)
    if not report.ok:
        raise SystemExit(f"verification failed after sealing:\n{report.describe()}")
    digest = hashlib.sha256(Path(written).read_bytes()).hexdigest()
    print(f"\nSealed {name}: {len(manifest.test_ids)} test, {len(manifest.train_ids)} "
          f"train -> {written}")
    print(f"  manifest sha256 : {digest}")
    print("  verified        : the table reproduces the manifest")


def eligible_patients(cohort: pd.DataFrame, expr_ids) -> tuple[pd.DataFrame, dict]:
    """Apply the registered eligibility rule and return (table, report).

    A row is eligible when its patient id is present and has an expression
    column, its survival time and vital status are both recorded, and it
    belongs to an RNA-seq batch. One row per patient is kept (first in stored
    order after sorting by batch then id); any collapse is counted.

    The returned table holds ``OUTPUT_COLUMNS`` only.
    """
    for col in ("patient_id", "overall_survival_months", "vital_status"):
        if col not in cohort.columns:
            raise ValueError(
                f"cohort has no '{col}' column; columns are {list(cohort.columns)}."
            )

    ingest = load_ingest_module()
    expr_ids = {str(e) for e in expr_ids}
    # The ingest's own counting rule, fed the same column names.
    overlap = ingest.expression_overlap(cohort, SimpleNamespace(columns=sorted(expr_ids)))

    ids = cohort["patient_id"]
    survival = pd.to_numeric(cohort["overall_survival_months"], errors="coerce")
    has_id = ids.notna()
    has_survival = survival.notna() & cohort["vital_status"].notna()
    has_expression = has_id & ids.astype(str).isin(expr_ids)
    if "cgga_dataset" in cohort.columns:
        is_rna_seq = cohort["cgga_dataset"].isin(ingest.EXPRESSION_DATASETS)
    else:
        # No batch label: the expression column itself is the RNA-seq evidence,
        # because counts.parquet is assembled from the RNA-seq batches only.
        is_rna_seq = pd.Series(True, index=cohort.index)

    eligible = has_id & has_expression & has_survival & is_rna_seq
    rows = cohort.loc[eligible]
    order = [c for c in ("cgga_dataset", "patient_id") if c in rows.columns]
    if order:
        rows = rows.sort_values(order, kind="mergesort")
    n_rows = int(len(rows))
    rows = rows.drop_duplicates(subset="patient_id", keep="first")
    out = rows.loc[:, [c for c in OUTPUT_COLUMNS if c in rows.columns]].copy()
    out = out.sort_values("patient_id", kind="mergesort").reset_index(drop=True)

    by_batch: dict[str, int] = {}
    if "cgga_dataset" in out.columns:
        by_batch = {
            str(k): int(v)
            for k, v in out["cgga_dataset"].value_counts().sort_index().items()
        }
    report = {
        "cohort_rows": int(len(cohort)),
        "cohort_patients": int(ids[has_id].nunique()),
        "rna_seq_rows": int(is_rna_seq.sum()),
        "microarray_rows": int((~is_rna_seq).sum()),
        "with_expression": int(overlap["with_expression"]),
        "with_survival": int(overlap["with_survival"]),
        "with_expression_and_survival": int(overlap["with_expression_and_survival"]),
        "eligible_rows": n_rows,
        "eligible_patients": int(len(out)),
        "duplicate_rows_collapsed": n_rows - int(len(out)),
        "rna_seq_rows_without_expression_column": int((is_rna_seq & ~has_expression).sum()),
        "rna_seq_rows_without_survival": int(
            (is_rna_seq & has_expression & ~has_survival).sum()
        ),
        "by_batch": by_batch,
    }
    return out, report


def print_report(report: dict) -> None:
    """ASCII-only counts. No clinical value is printed."""
    print("CGGA external-partition eligibility (counts only):")
    print(f"  Cohort rows / unique patients      : {report['cohort_rows']} / "
          f"{report['cohort_patients']}")
    print(f"  RNA-seq batch rows                 : {report['rna_seq_rows']}")
    print(f"  Microarray batch rows              : {report['microarray_rows']}")
    print(f"  Patients with an expression column : {report['with_expression']}")
    print("  Patients with survival time AND vital status : "
          f"{report['with_survival']}")
    print(f"  Patients with BOTH                 : {report['with_expression_and_survival']}")
    print("  RNA-seq rows without an expression column    : "
          f"{report['rna_seq_rows_without_expression_column']}")
    print("  RNA-seq rows with expression but no survival : "
          f"{report['rna_seq_rows_without_survival']}")
    print(f"  Eligible rows / eligible patients  : {report['eligible_rows']} / "
          f"{report['eligible_patients']}")
    print(f"  Duplicate rows collapsed           : {report['duplicate_rows_collapsed']}")
    for batch, n in report["by_batch"].items():
        print(f"    {batch}: {n} patients")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="List the CGGA patients eligible for the sealed external partition."
    )
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT), type=Path,
                        help="Mapped CGGA cohort Parquet file.")
    parser.add_argument("--counts", default=str(DEFAULT_COUNTS), type=Path,
                        help="CGGA counts matrix; only its column names are read.")
    parser.add_argument("--out", default=str(DEFAULT_OUT), type=Path,
                        help="Where the eligible table is written (identifiers only).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report the counts only; write nothing.")
    parser.add_argument("--ledger", default=str(DEFAULT_LEDGER), type=Path,
                        help="Access ledger the run is recorded in.")
    parser.add_argument("--seal", default=None, type=Path, metavar="MANIFEST",
                        help="After writing the table, freeze it as the all-test "
                             "partition at MANIFEST and verify it.")
    parser.add_argument("--name", default=DEFAULT_PARTITION_NAME,
                        help="Partition name for --seal.")
    args = parser.parse_args(argv)
    if args.dry_run and args.seal:
        parser.error("--seal writes a manifest; it cannot be combined with --dry-run")
    if args.seal and Path(args.seal).exists():
        raise SystemExit(f"{args.seal} exists; a manifest is never overwritten here "
                         "(docs/preregistration.md section 4.4). Nothing was read.")

    entry = declare_access(args.cohort, args.counts, args.ledger, dry_run=args.dry_run,
                           seal=str(args.seal) if args.seal else None)
    print(f"Access recorded: ledger entry {entry['seq']} ({entry['recorded_at']}), "
          f"kind {entry['kind']}.\n")

    cohort = read_cohort_columns(args.cohort)
    expr_ids = expression_sample_ids(args.counts)
    out, report = eligible_patients(cohort, expr_ids)
    print_report(report)

    if args.dry_run:
        print("\nDRY RUN: nothing was written. Re-run without --dry-run to write the "
              "eligible table.")
        return
    if out.empty:
        raise SystemExit("No eligible patients; refusing to write an empty table.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    digest = hashlib.sha256(args.out.read_bytes()).hexdigest()
    print(f"\nWrote {len(out)} patients -> {args.out}")
    print(f"  columns : {list(out.columns)}")
    print(f"  sha256  : {digest}")
    if not args.seal:
        print("\nNothing is sealed. Re-run with --seal, or use the separate step, which "
              "refuses to overwrite a sealed manifest:")
        print("  gliostrata freeze-partition --input <file above> --manifest "
              "data/partitions/cgga_external_v1.json --name cgga_external_v1 --all-test")
        return
    seal_partition(out, args.out, args.seal, args.name, ledger_seq=entry["seq"])


if __name__ == "__main__":
    main()
