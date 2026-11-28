"""Freeze the TCGA reference the tropism score is standardised against.

`gliostrata.simulation.delivery.tropism_score` standardises each carrier's chemoattractant
genes against a **frozen TCGA mean and standard deviation, never a CGGA cohort statistic**
(pre-registration item 11(h); Design A item 2 of `docs/t4t5_model_design.md`). A per-cohort
statistic would let the sealed external cohort leak into the transform, which section 3.1
forbids.

This builds that reference from TCGA counts and writes
`configs/response_surface_v1_tropism_reference.json`, to be checksummed in section 8 at the
freeze alongside the parameter file.

Scale: **log2(FPKM + 1)**, taken from the `fpkm_unstranded` column of the cached GDC files,
with the counts of every Ensembl id sharing a GENCODE `gene_name` summed as elsewhere in the
project. FPKM rather than CPM because the cohort this reference will be used to standardise
(CGGA) publishes RSEM FPKM, and a reference on a different normalisation leaves a systematic
offset in the score: measured on these same 284 samples it is -1.36 SD for the NSC
chemoattractant panel and -1.47 SD for TGFB1, which alone would cut the engrafted-dose factor
by about a third for every patient. Item 11(h) requires a frozen TCGA mean and SD and names no
normalisation, so matching the scoring cohort is an implementation choice inside the registered
design. What remains, and is disclosed rather than fixed, is the second-order difference
between GDC's STAR-derived FPKM and CGGA's RSEM FPKM.

**Only TCGA is read.** No CGGA path is opened, so no read-guard access is recorded. The script
refuses to overwrite an existing reference, because a frozen reference is frozen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path("data/raw/tcga_gbm/expression/_cache")
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
    for path in (CACHE, UNIVERSE, PARAMS):
        if not path.exists():
            raise SystemExit(f"{path} not found")
    if "cgga" in str(CACHE).lower():
        raise SystemExit("the reference must be built from TCGA, never from CGGA")

    import yaml
    params = yaml.safe_load(PARAMS.read_text(encoding="utf-8"))
    genes_by_carrier = params["delivery"]["tropism_genes"]
    wanted = sorted({g for genes in genes_by_carrier.values() for g in genes})
    print(f"tropism genes, from the parameter file: {wanted}")

    files = sorted(CACHE.glob("*.tsv"))
    if not files:
        raise SystemExit(f"no cached GDC files under {CACHE}")
    print(f"reading fpkm_unstranded from {len(files)} cached GDC files")

    columns: list[pd.Series] = []
    for path in files:
        frame = pd.read_csv(path, sep="\t", comment="#",
                            usecols=["gene_name", "fpkm_unstranded"])
        frame = frame.dropna(subset=["fpkm_unstranded"])
        # Sum the ids sharing a GENCODE gene_name, as collapse_to_symbols does for counts.
        columns.append(frame.groupby("gene_name")["fpkm_unstranded"].sum().rename(path.stem))
    fpkm = pd.concat(columns, axis=1)
    missing = [g for g in wanted if g not in fpkm.index]
    if missing:
        raise SystemExit(f"tropism genes absent from TCGA after symbol collapse: {missing}")
    if fpkm.loc[wanted].isna().to_numpy().any():
        raise SystemExit("a tropism gene is missing from at least one sample")
    expression = np.log2(fpkm + 1.0)
    print(f"log2(FPKM + 1): {expression.shape[0]:,} symbols x {expression.shape[1]} samples")

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
        "scale": "log2(FPKM + 1)",
        "fpkm_column": "fpkm_unstranded (GDC STAR-derived)",
        "sources": {
            "cache_dir": CACHE.as_posix(), "n_files": len(files),
            "file_manifest_sha256": _sha256(
                Path("data/raw/tcga_gbm/expression/file_manifest.parquet")),
        },
        "symbol_collapse": "counts of every Ensembl id sharing a GENCODE gene_name are summed",
        "genes_by_carrier": genes_by_carrier,
        "reference": rows,
        "scale_rationale": (
            "FPKM, not CPM, because the cohort this reference standardises (CGGA) publishes "
            "RSEM FPKM. On these same 284 samples a CPM reference offsets the score by -1.36 "
            "SD for the NSC panel and -1.47 SD for TGFB1, cutting the engrafted-dose factor "
            "by about a third for every patient. Item 11(h) names no normalisation, so "
            "matching the scoring cohort is an implementation choice inside the registered "
            "design."
        ),
        "caveat_remaining": (
            "TCGA's FPKM here is GDC's STAR-derived fpkm_unstranded; CGGA's is RSEM FPKM. "
            "Both are length-normalised, so the first-order gene-length offset is removed, "
            "but a second-order pipeline difference remains and is reported, not corrected."
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
