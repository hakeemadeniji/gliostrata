"""Full cohort verification of the active partition manifests: the script behind
the second step of the ``partition-drift`` CI job.

Until 2026-11-27 this was a shell loop written inline in
``.github/workflows/ci.yml``, where no test could run it. It is the same logic.

For each manifest in ``ACTIVE_MANIFESTS`` (space-separated, as in the first
step), read its ``cohort_source``. If that file exists, run

    gliostrata freeze-partition --input <cohort> --manifest <manifest> --verify

and otherwise print a "Skipping" line (the cohort parquet is gitignored, so in
CI it is normally absent). The first failing verification stops the run with
that command's exit status, as the shell loop did under ``bash -e``.

Held-out guard. A manifest whose ``cohort_source`` lies under
``data/processed/cgga`` or ``data/raw/cgga`` names the sealed external cohort.
Verifying it reads CGGA patient identifiers, which is an access that has to be
recorded in the access ledger (``data/access_ledger.jsonl``). This script
therefore REFUSES such a manifest, when the file is present, unless
``--allow-held-out`` is given. In CI the file is never present, so the manifest
is skipped exactly as before and nothing changes. The guard tests whether the
path exists and never looks at its contents.

Exit status: 0 when every present cohort verifies (or is skipped), 1 when a
manifest is refused, 2 on an unreadable manifest, otherwise the status of the
failing ``freeze-partition`` call.

Usage:
    ACTIVE_MANIFESTS="data/partitions/tcga_gbm_clinical_v2.json" \\
        python scripts/ci_partition_verify.py
"""

from __future__ import annotations

import argparse
import json
import os
import posixpath
import subprocess
import sys
from pathlib import Path

#: Directories (as consecutive path components) that hold the held-out cohort.
HELD_OUT_DIRS = (("data", "processed", "cgga"), ("data", "raw", "cgga"))


def is_held_out(cohort_source) -> bool:
    """True when ``cohort_source`` lies under a held-out CGGA directory.

    A pure string test on the normalised path: nothing is opened or listed.
    Works for relative paths and for absolute ones that contain the same
    consecutive components (so a temporary copy named like the CGGA path is
    treated like it). Case-insensitive, and ``..`` segments are collapsed.
    """
    text = str(cohort_source).replace("\\", "/")
    parts = [p.lower() for p in posixpath.normpath(text).split("/") if p not in ("", ".")]
    return any(
        tuple(parts[i : i + 3]) == held
        for held in HELD_OUT_DIRS
        for i in range(len(parts) - 2)
    )


def read_cohort_source(manifest) -> str:
    with open(manifest, encoding="utf-8") as fh:
        return json.load(fh)["cohort_source"]


def _run_verify(cmd) -> int:
    return subprocess.run(cmd).returncode


def verify_active(manifests, *, allow_held_out=False, say=print, run=_run_verify) -> int:
    """Verify each manifest's cohort where present; return the exit status."""
    status = 0
    for manifest in manifests:
        try:
            cohort = read_cohort_source(manifest)
        except (OSError, ValueError, KeyError) as exc:
            say(f"Cannot read cohort_source from {manifest}: {exc}")
            return 2
        if not os.path.isfile(cohort):
            say(
                f"Skipping: cohort {cohort} is not present in this checkout "
                f"(data/raw/** is gitignored)."
            )
            continue
        if is_held_out(cohort) and not allow_held_out:
            say(
                f"REFUSED: cohort {cohort} for {manifest} is held-out CGGA material. "
                f"Verifying it reads patient identifiers; record the access in "
                f"data/access_ledger.jsonl and re-run with --allow-held-out."
            )
            status = status or 1
            continue
        say(f"Verifying {manifest} against {cohort}")
        code = run(
            [
                "gliostrata", "freeze-partition",
                "--input", cohort,
                "--manifest", str(manifest),
                "--verify",
            ]
        )
        if code != 0:
            return code
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--allow-held-out",
        action="store_true",
        help="verify a manifest whose cohort lies under data/{processed,raw}/cgga "
        "(reads CGGA identifiers; record it in the access ledger)",
    )
    args = parser.parse_args(argv)
    manifests = [Path(p).as_posix() for p in os.environ.get("ACTIVE_MANIFESTS", "").split()]
    return verify_active(manifests, allow_held_out=args.allow_held_out)


if __name__ == "__main__":
    sys.exit(main())
