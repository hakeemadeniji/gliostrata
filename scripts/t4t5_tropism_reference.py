"""Freeze the TCGA reference the tropism score is standardised against.

`gliostrata.simulation.delivery.tropism_score` standardises each carrier's chemoattractant
genes against a **frozen TCGA mean and standard deviation, never a CGGA cohort statistic**
(pre-registration item 11(h); Design A item 2 of `docs/t4t5_model_design.md`). A per-cohort
statistic would let the sealed external cohort leak into the transform, which section 3.1
forbids.

This builds that reference from TCGA counts and writes
`configs/response_surface_v1_tropism_reference.json`, to be checksummed in section 8 at the
freeze alongside the parameter file.

Scale and universe are the project's existing ones, not new choices: symbols are collapsed by
`collapse_to_symbols` and values are `log2(CPM + 1)` from `log2_cpm`, with the library size
summed over the same 33,212-symbol universe the T2 comparator used
(`configs/t2_comparator/gene_universe.txt`). Reusing them means a value here means the same
thing as a value anywhere else in the project.

**Only TCGA is read.** No CGGA path is opened, so no read-guard access is recorded. The script
refuses to overwrite an existing reference, because a frozen reference is frozen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from gliostrata.prediction.geneset_comparator import collapse_to_symbols, log2_cpm

COUNTS = Path("data/raw/tcga_gbm/expression/counts.parquet")
ANNOTATION = Path("data/raw/tcga_gbm/expression/gene_annotation.parquet")
UNIVERSE = Path("configs/t2_comparator/gene_universe.txt")
OUT = Path("configs/response_surface_v1_tropism_reference.json")
PARAMS = Path("configs/response_surface_v1.yaml")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(force: bool = False) -> None:
    if OUT.exists() and not force:
        raise SystemExit(
            f"{OUT} exists and a frozen reference is never overwritten. Pass --force only to "
            "rebuild one that has not yet been checksummed in section 8.")
    for path in (COUNTS, ANNOTATION, UNIVERSE, PARAMS):
        if not path.exists():
            raise SystemExit(f"{path} not found")
    if "cgga" in str(COUNTS).lower():
        raise SystemExit("the reference must be built from TCGA, never from CGGA")

    import yaml
    params = yaml.safe_load(PARAMS.read_text(encoding="utf-8"))
    genes_by_carrier = params["delivery"]["tropism_genes"]
    wanted = sorted({g for genes in genes_by_carrier.values() for g in genes})
    print(f"tropism genes, from the parameter file: {wanted}")

    counts = pd.read_parquet(COUNTS)
    annotation = pd.read_parquet(ANNOTATION)
    universe = [line.strip() for line in UNIVERSE.read_text(encoding="utf-8").splitlines()
                if line.strip()]
    print(f"TCGA counts: {counts.shape[0]:,} ids x {counts.shape[1]} samples")

    by_symbol = collapse_to_symbols(counts, annotation)
    missing = [g for g in wanted if g not in by_symbol.index]
    if missing:
        raise SystemExit(f"tropism genes absent from TCGA after symbol collapse: {missing}")
    expression = log2_cpm(by_symbol, universe)
    print(f"log2(CPM + 1) over the {len(universe):,}-symbol universe: "
          f"{expression.shape[0]:,} genes x {expression.shape[1]} samples")

    rows = {}
    for gene in wanted:
        values = expression.loc[gene].astype(float)
        sd = float(values.std(ddof=1))
        if not sd > 0:
            raise SystemExit(f"{gene} has a non-positive standard deviation on TCGA")
        rows[gene] = {"mean": float(values.mean()), "sd": sd,
                      "min": float(values.min()), "max": float(values.max())}
        print(f"  {gene:<8} mean {rows[gene]['mean']:7.4f}  sd {sd:6.4f}  "
              f"range {rows[gene]['min']:.3f} to {rows[gene]['max']:.3f}")

    document = {
        "what": ("The frozen TCGA reference the T4/T5 tropism score is standardised against "
                 "(pre-registration item 11(h)). Never rebuilt from a CGGA statistic."),
        "registered_under": "docs/preregistration.md section 2.2 item 11(h) and 11(n)",
        "built_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "cohort": "TCGA-GBM",
        "n_samples": int(expression.shape[1]),
        "scale": "log2(CPM + 1)",
        "library_size_universe": {
            "file": UNIVERSE.as_posix(), "n_symbols": len(universe),
            "sha256": _sha256(UNIVERSE),
        },
        "sources": {
            "counts": {"file": COUNTS.as_posix(), "sha256": _sha256(COUNTS)},
            "annotation": {"file": ANNOTATION.as_posix(), "sha256": _sha256(ANNOTATION)},
        },
        "symbol_collapse": "counts of every Ensembl id sharing a GENCODE gene_name are summed",
        "genes_by_carrier": genes_by_carrier,
        "reference": rows,
        "caveat_scale": (
            "These are log2(CPM + 1) statistics from TCGA counts. A later cohort scored on a "
            "different normalisation (CGGA's release is RSEM FPKM) is not on the same scale, "
            "so the standardised score carries a systematic offset. The offset is common to "
            "every patient and so does not reorder them, but it does move the absolute level "
            "of the engrafted dose. Stated here so that it is reported rather than discovered."
        ),
    }
    OUT.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"\nWrote {OUT}\n  sha256 {_sha256(OUT)}")
    print("Not yet registered: its checksum goes in section 8 at the freeze "
          "(item 11(e) and (m)).")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("step", choices=("build",))
    parser.add_argument("--force", action="store_true",
                        help="rebuild a reference that has not yet been checksummed")
    args = parser.parse_args(argv)
    build(force=args.force)


if __name__ == "__main__":
    main()
