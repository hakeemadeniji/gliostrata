"""
Export the TCGA-GBM bulk expression matrix as a CIBERSORTx mixture file.

CIBERSORTx (https://cibersortx.stanford.edu) is a web tool: it cannot be driven
from this pipeline, so L2 deconvolution runs as a manual step. This script
produces the upload file; the returned fractions are re-ingested separately.

CIBERSORTx mixture-file requirements, and how each is satisfied here:

  tab-delimited plain text      -> .txt, sep="\t"
  genes in rows, samples in cols-> rows are HGNC symbols, cols are patient ids
  first column = gene symbol    -> header cell "GeneSymbol"
  HGNC symbols, not Ensembl     -> mapped via gene_annotation.parquet
  NO duplicate gene symbols     -> collapsed, keeping the gene_id with the
                                   highest mean TPM (deterministic)
  linear space, NOT log2        -> TPM is already linear; nothing is logged
  no missing values             -> verified zero NaN before writing
  no all-zero rows              -> dropped (they carry no signal and inflate
                                   the upload)

Usage:
    python scripts/export_cibersortx.py --out data/processed/cibersortx
    python scripts/export_cibersortx.py --gene-set protein_coding
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_TPM = Path("data/raw/tcga_gbm/expression/tpm.parquet")
DEFAULT_ANNOTATION = Path("data/raw/tcga_gbm/expression/gene_annotation.parquet")
DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")

#: Header of the first column. CIBERSORTx reads column 0 as the gene identifier.
GENE_COLUMN_HEADER = "GeneSymbol"

#: Gene sets selectable with --gene-set.
GENE_SETS = {
    "all": None,
    "protein_coding": ("protein_coding",),
    "protein_coding_lncrna": ("protein_coding", "lncRNA"),
}


def map_symbols(tpm: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Replace the Ensembl gene_id index with HGNC symbols.

    Genes with no symbol are dropped; the caller reports how many.
    """
    symbols = annotation.set_index("gene_id")["gene_name"]
    mapped = tpm.copy()
    mapped.index = mapped.index.map(symbols)
    mapped = mapped[mapped.index.notna()]
    mapped = mapped[mapped.index.astype(str).str.strip() != ""]
    mapped.index.name = GENE_COLUMN_HEADER
    return mapped


def collapse_duplicate_symbols(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Collapse rows sharing a gene symbol, keeping the highest-expressed row.

    Several Ensembl ids map to one HGNC symbol (PAR_Y copies, readthrough
    entries, retired ids). CIBERSORTx requires unique symbols. Keeping the row
    with the greatest mean TPM retains the dominant transcript rather than
    averaging a real signal together with a near-silent duplicate. Ties break on
    the original row order (stable sort), so the result is deterministic and
    independent of the order the rows arrive in.

    The reordering is positional (``iloc``) on purpose. ``frame.loc[labels]``
    with duplicate labels returns EVERY row carrying each label, in frame
    order, which silently reduced this function to "first row in frame order
    wins" -- a real bug caught by ``tests/test_export_cibersortx.py``.
    """
    if not frame.index.duplicated().any():
        return frame.sort_index(kind="mergesort"), 0
    means = frame.to_numpy(dtype=float, copy=False).mean(axis=1)
    order = np.argsort(-means, kind="stable")
    ranked = frame.iloc[order]
    collapsed = ranked[~ranked.index.duplicated(keep="first")]
    dropped = len(frame) - len(collapsed)
    return collapsed.sort_index(kind="mergesort"), dropped


def build_mixture(
    tpm: pd.DataFrame,
    annotation: pd.DataFrame,
    *,
    gene_types: tuple[str, ...] | None = None,
    patients: list[str] | None = None,
    decimals: int = 4,
    renormalize: bool = True,
) -> tuple[pd.DataFrame, dict]:
    """Build the CIBERSORTx mixture matrix and a report of what was done."""
    report: dict = {"input_genes": int(len(tpm)), "input_samples": int(tpm.shape[1])}

    if gene_types is not None:
        keep_ids = annotation.loc[annotation["gene_type"].isin(gene_types), "gene_id"]
        tpm = tpm.loc[tpm.index.intersection(keep_ids)]
        report["after_gene_type_filter"] = int(len(tpm))

    if patients is not None:
        selected = [p for p in tpm.columns if p in set(patients)]
        tpm = tpm[selected]
        report["after_patient_filter"] = int(tpm.shape[1])

    mapped = map_symbols(tpm, annotation)
    report["unmapped_genes_dropped"] = int(len(tpm) - len(mapped))

    collapsed, dropped = collapse_duplicate_symbols(mapped)
    report["duplicate_symbol_rows_dropped"] = int(dropped)

    nonzero = collapsed[collapsed.sum(axis=1) > 0]
    report["all_zero_genes_dropped"] = int(len(collapsed) - len(nonzero))

    column_mass = nonzero.sum(axis=0)
    report["column_sum_min"] = float(column_mass.min())
    report["column_sum_max"] = float(column_mass.max())
    if renormalize:
        # Filtering to a gene subset leaves each sample summing to a different
        # total, because the discarded mass differs per sample. In this cohort
        # that spread is severe and systematic: 121 of 284 samples above the
        # 0.3 misc_RNA fraction carry a median 50.0% of their TPM in misc_RNA,
        # versus a median 0.11% in the other 163 -- a library-prep difference
        # rather than biology. (The 0.3 threshold is a report parameter, not a
        # claim about the underlying distribution.) Deconvolution reads relative
        # abundance, so an uncorrected mixture would hand CIBERSORTx a
        # batch artefact. Rescaling each sample to a common total restores
        # comparability across the retained gene universe. It cannot repair
        # the underlying difference in what was sequenced -- only L1
        # batch correction can address that.
        nonzero = nonzero.div(column_mass, axis=1) * 1e6
        report["renormalized_to"] = 1e6

    rounded = nonzero.round(decimals)
    report["final_genes"] = int(len(rounded))
    report["final_samples"] = int(rounded.shape[1])
    report["nan_values"] = int(rounded.isna().sum().sum())
    report["negative_values"] = int((rounded < 0).sum().sum())
    return rounded, report


def write_mixture(frame: pd.DataFrame, path: Path) -> Path:
    """Write the tab-delimited mixture file CIBERSORTx expects."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=True, index_label=GENE_COLUMN_HEADER)
    return path


def print_report(report: dict, path: Path) -> None:
    """ASCII-only summary (cp1252 consoles choke on non-ASCII)."""
    print("\nCIBERSORTx mixture export:")
    print(f"  Input genes: {report['input_genes']}")
    if "after_gene_type_filter" in report:
        print(f"  After gene_type filter: {report['after_gene_type_filter']}")
    if "after_patient_filter" in report:
        print(f"  After patient filter: {report['after_patient_filter']} samples")
    print(f"  Dropped, no HGNC symbol: {report['unmapped_genes_dropped']}")
    print(f"  Dropped, duplicate symbol: {report['duplicate_symbol_rows_dropped']}")
    print(f"  Dropped, all-zero across samples: {report['all_zero_genes_dropped']}")
    print(f"  Final matrix: {report['final_genes']} genes x {report['final_samples']} samples")
    print(f"  NaN values: {report['nan_values']}  Negative values: {report['negative_values']}")
    print(
        "  Column sums before renormalization: "
        f"{report['column_sum_min']:,.0f} .. {report['column_sum_max']:,.0f}"
    )
    if "renormalized_to" in report:
        print(f"  Rescaled every sample to {report['renormalized_to']:,.0f}")
    size_mb = path.stat().st_size / (1024 * 1024)
    print(f"  Wrote: {path}  ({size_mb:.1f} MB)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpm", default=str(DEFAULT_TPM))
    parser.add_argument("--annotation", default=str(DEFAULT_ANNOTATION))
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT))
    parser.add_argument("--out", default="data/processed/cibersortx", help="Output directory")
    parser.add_argument(
        "--gene-set",
        choices=sorted(GENE_SETS),
        default="protein_coding",
        help="Which genes to include (default: protein_coding)",
    )
    parser.add_argument(
        "--cohort-only",
        action="store_true",
        help="Restrict to patients present in the clinical cohort",
    )
    parser.add_argument("--decimals", type=int, default=4)
    parser.add_argument(
        "--no-renormalize",
        action="store_true",
        help="Skip rescaling each sample to a common total after gene filtering",
    )
    args = parser.parse_args()

    tpm = pd.read_parquet(args.tpm)
    annotation = pd.read_parquet(args.annotation)

    patients = None
    if args.cohort_only:
        patients = list(pd.read_parquet(args.cohort)["patient_id"])

    mixture, report = build_mixture(
        tpm,
        annotation,
        gene_types=GENE_SETS[args.gene_set],
        patients=patients,
        decimals=args.decimals,
        renormalize=not args.no_renormalize,
    )
    suffix = "_cohort" if args.cohort_only else ""
    norm = "_raw" if args.no_renormalize else ""
    out_path = Path(args.out) / f"tcga_gbm_mixture_{args.gene_set}{suffix}{norm}.txt"
    write_mixture(mixture, out_path)
    print_report(report, out_path)


if __name__ == "__main__":
    main()
