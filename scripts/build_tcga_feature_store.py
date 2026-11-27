"""
Build the TCGA feature store: a thin wrapper over the CLI.

    python scripts/build_tcga_feature_store.py [OPTIONS]

is exactly

    gliostrata build-feature-store --mode real [OPTIONS]

which runs the real L2 (single-sample subtype calls against the frozen
centroid set; state and TME enrichment scores) and the L3 store on the TCGA
expression input, and writes the store (every cohort patient, one row each),
its label table, its missing-value reasons, its column metadata and a build
report with the flow counts of the named views. Run
`gliostrata harmonize-clinical` first: it writes the L1 clinical and label
tables this reads.

There is one code path: gliostrata.features.build.run_real. The earlier
version of this script carried its own assembly and dropped the cohort's
null placeholder columns itself; that workaround is gone, because the L1
clinical table no longer carries those names and build_feature_store refuses
a null-only column that shadows a computed feature. Any path containing
"cgga" is refused.
"""
from __future__ import annotations

import sys

from gliostrata.cli import cli


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if "--mode" in args:
        raise SystemExit("this wrapper always runs --mode real; do not pass --mode")
    cli.main(args=["build-feature-store", "--mode", "real", *args],
             prog_name="build_tcga_feature_store.py", standalone_mode=True)


if __name__ == "__main__":
    main()
