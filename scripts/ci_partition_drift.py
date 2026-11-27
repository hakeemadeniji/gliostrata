"""Partition drift guardrail: the script behind the ``partition-drift`` CI job.

Until 2026-11-27 this logic was written inline in ``.github/workflows/ci.yml``,
where the only test of it read the YAML as text. It is the same code, moved into
a file that both the workflow and ``tests/test_ci_partition_drift.py`` execute.

A sealed train/test split is a pre-registration commitment: the split is a pure
function of (seed, patient_id), so it must keep reproducing as the cohort is
re-ingested. If a patient moves between train and test, every downstream
evaluation is contaminated.

Two levels of check, for each manifest under ``data/partitions/*.json``:

1. ALWAYS, needing only git contents. The manifest is sealed, its id digests
   match its id lists, and re-deriving each patient's split from (seed,
   patient_id) reproduces the frozen assignment.
2. WHEN AVAILABLE. The cohort parquet the manifest names: if it is present the
   script says so, and the job then runs ``gliostrata freeze-partition
   --verify``. ``data/raw/**`` is gitignored, so in CI it is normally absent.

``ACTIVE_MANIFESTS`` is a space-separated list of the manifests that are current.
A manifest not listed is integrity-checked and labelled superseded: v1 (583
patients) is the historical record of the pre-filter-revision cohort, kept
unmodified and deliberately not re-frozen. A listed manifest that does not exist
stops the job.

Exit status is 0 when every manifest verifies and 1 otherwise; stdout says which
and why. Nothing is written and no socket is opened.

Usage:
    ACTIVE_MANIFESTS=data/partitions/tcga_gbm_clinical_v2.json \\
        python scripts/ci_partition_drift.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from gliostrata.evaluation import assign_split, read_manifest, sha256_over_ids

PARTITION_DIR = Path("data/partitions")


def verify(partitions_dir=PARTITION_DIR, active=(), say=print) -> list[str]:
    """Check every manifest in ``partitions_dir``; return the failures.

    ``say`` receives the report, one line at a time. A missing listed manifest,
    or an empty directory, raises ``SystemExit`` with the reason: there is
    nothing to verify, which is not the same as a pass.
    """
    partitions_dir = Path(partitions_dir)
    active = [Path(p) for p in active]
    manifests = sorted(partitions_dir.glob("*.json"))
    if not manifests:
        raise SystemExit(f"No partition manifests found in {partitions_dir}/.")
    for path in active:
        if not path.is_file():
            raise SystemExit(
            f"ACTIVE_MANIFESTS lists {path.as_posix()}, which does not exist."
        )

    failures: list[str] = []
    for path in manifests:
        m = read_manifest(path)
        role = "ACTIVE" if path in active else "superseded (historical record)"
        shown = path.as_posix()  # the same text on every platform
        say(f"\n=== {shown}  [{role}] ===")
        say(f"  partition    : {m.partition_name}")
        say(f"  sealed       : {m.sealed}")
        say(f"  seed         : {m.seed}")
        say(f"  n_train/test : {m.n_train} / {m.n_test}")

        if not m.sealed:
            failures.append(f"{shown}: manifest is UNSEALED (draft)")

        for label, ids, digest in (
            ("train_ids", m.train_ids, m.train_ids_sha256),
            ("test_ids", m.test_ids, m.test_ids_sha256),
            ("all_ids", list(m.train_ids) + list(m.test_ids), m.all_ids_sha256),
        ):
            actual = sha256_over_ids(ids)
            if actual != digest:
                failures.append(
                    f"{shown}: {label} digest mismatch "
                    f"(recorded {digest[:12]}..., recomputed {actual[:12]}...)"
                )

        if m.stratify_by is None:
            drifted = [
                pid
                for expected, ids in (("train", m.train_ids), ("test", m.test_ids))
                for pid in ids
                if assign_split(pid, seed=m.seed, test_fraction=m.test_fraction) != expected
            ]
            if drifted:
                failures.append(
                    f"{shown}: {len(drifted)} id(s) no longer hash to their "
                    f"frozen split, e.g. {drifted[:5]}"
                )
            else:
                say(
                    f"  re-derived   : all {m.n_train + m.n_test} ids reproduce "
                    f"their frozen split from (seed, patient_id)"
                )
        else:
            say(
                f"  re-derived   : SKIPPED (stratified by {m.stratify_by!r}; "
                f"needs the cohort table)"
            )

        cohort = Path(m.cohort_source)
        if path not in active:
            say(
                "  cohort       : not checked -- superseded manifests are frozen "
                "history and are EXPECTED to drift from the current cohort"
            )
        elif cohort.exists():
            say(f"  cohort       : present at {cohort} -> full --verify below")
        else:
            say(
                f"  cohort       : NOT in this checkout ({cohort} is gitignored); "
                f"full --verify not run in CI"
            )

    if failures:
        say("\nPARTITION DRIFT:")
        for failure in failures:
            say(f"  - {failure}")
    else:
        say("\nAll partition manifests verified.")
    return failures


def main(argv: list[str] | None = None) -> int:
    active = os.environ.get("ACTIVE_MANIFESTS", "").split()
    return 1 if verify(active=active) else 0


if __name__ == "__main__":
    # A bad manifest must fail the job, not print a traceback and still exit 0.
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 - the job reports, the log keeps the detail
        print(f"PARTITION DRIFT: could not verify the manifests ({exc})", file=sys.stderr)
        sys.exit(1)
