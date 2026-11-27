"""Prepare the published Verhaak reference for the single-sample classifier.

One-off preparation from downloaded public files (stage 3b, decision D6). It
needs ``openpyxl`` and ``pyreadr``, which are not project dependencies; the
scripts downstream of it read only what it writes.

Inputs (each under data/external/, git-ignored; provenance with URL and
SHA-256 in data/provenance/subtype_reference_sources.json):
  --verhaak-core   TCGA_unified_CORE_ClaNC840.txt (Verhaak et al. 2010, GDC
                   publication page "gbm_exp"): the 840-gene list, and the
                   labels of the 173 core samples
  --brennan        Molecular_subtype_classification.xlsx (Brennan et al. 2013,
                   GDC page "gbm_2013"): one expression subtype per TCGA sample
  --ceccarelli     TableS1.PatientData.20151020.v3.xlsx (Ceccarelli et al.
                   2016, GDC page "lgggbm_2016"): "Original Subtype" and
                   "Transcriptome Subtype", used as a cross-check only
  --hgnc-rda       hgnc.table.rda from the CRAN package HGNChelper 0.8.16: the
                   map from previous HGNC symbols to approved symbols
  --annotation     data/raw/tcga_gbm/expression/gene_annotation.parquet

Outputs:
  --gene-map       configs/subtype_centroids/clanc840_genes.json (committed):
                   for each of the 840 published symbols, the current symbol it
                   resolves to in the annotation, and how; no patient data
  --labels         data/processed/tcga_gbm/published_subtype_labels.parquet
                   (git-ignored): one row per TCGA patient with a published
                   label

Symbol rule, fixed before any centroid was built: a published symbol present in
the annotation is used as is ("direct"); otherwise it is replaced by its HGNC
approved symbol when exactly one approved symbol for it is in the annotation
("hgnc_previous_symbol"); otherwise it is excluded ("ambiguous" or
"not_in_annotation"). Nothing is guessed.

Every input is parsed as data: nothing in it is executed. No CGGA path is
accepted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

FOUR = ("Classical", "Mesenchymal", "Proneural", "Neural")


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _refuse_cgga(*paths) -> None:
    for p in paths:
        if p is not None and "cgga" in str(p).lower():
            raise ValueError(f"refusing {p}: this script reads TCGA references only")


def read_verhaak_core(path: Path) -> tuple[list[str], dict[str, str]]:
    """The 840 genes in file order, and the core labels keyed by patient id."""
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    header = lines[0].split("\t")
    subtype_row = lines[1].split("\t")
    if header[0] != "CLID" or subtype_row[0] != "Subtype":
        raise ValueError("unexpected layout of the Verhaak core file")
    genes = [ln.split("\t")[0] for ln in lines[2:] if ln.strip()]
    if len(genes) != 840 or len(set(genes)) != 840:
        raise ValueError(f"expected 840 distinct genes, found {len(genes)}")
    labels = {h[:12]: s for h, s in zip(header[2:], subtype_row[2:])}
    return genes, labels


def resolve_genes(genes: list[str], annotation: pd.DataFrame, hgnc: pd.DataFrame) -> list[dict]:
    names = set(annotation["gene_name"].astype(str))
    approved = hgnc.groupby("Symbol")["Approved.Symbol"].apply(lambda s: sorted(set(s)))
    out = []
    for g in genes:
        if g in names:
            out.append({"published_symbol": g, "symbol": g, "resolution": "direct"})
            continue
        hits = [a for a in approved.get(g, []) if a in names]
        if len(hits) == 1:
            out.append({"published_symbol": g, "symbol": hits[0],
                        "resolution": "hgnc_previous_symbol"})
        else:
            out.append({"published_symbol": g, "symbol": None,
                        "resolution": "ambiguous" if hits else "not_in_annotation",
                        "candidates": hits})
    used = [r["symbol"] for r in out if r["symbol"]]
    if len(used) != len(set(used)):
        raise ValueError("two published symbols resolve to one current symbol")
    return out


def read_brennan(path: Path) -> pd.DataFrame:
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    header = list(rows[2][:4])
    if header != ["sample id", "Cluster", "CanCellPaper 2010", "Core"]:
        raise ValueError(f"unexpected Brennan header {header}")
    recs = [r[:4] for r in rows[3:] if r[0]]
    df = pd.DataFrame(recs, columns=["sample_id", "brennan2013_cluster",
                                     "verhaak2010_included", "verhaak2010_core"])
    df["patient_id"] = df["sample_id"].astype(str).str[:12]
    # One patient (TCGA-07-0227) has two rows, both entirely "NA": keep one.
    dup = df[df["patient_id"].duplicated(keep=False)]
    if not (dup["brennan2013_cluster"] == "NA").all():
        raise ValueError("a duplicated Brennan patient carries a label")
    return df.drop_duplicates("patient_id").drop(columns="sample_id")


def read_ceccarelli(path: Path) -> pd.DataFrame:
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    h = list(rows[1])
    i_case, i_study = h.index("Case"), h.index("Study")
    i_orig, i_trans = h.index("Original Subtype"), h.index("Transcriptome Subtype")
    i_idh = h.index("IDH status")
    recs = [(r[i_case], r[i_orig], r[i_trans], r[i_idh]) for r in rows[2:]
            if r[i_case] and r[i_study] == "Glioblastoma multiforme"]
    return pd.DataFrame(recs, columns=["patient_id", "ceccarelli2016_original",
                                       "ceccarelli2016_transcriptome",
                                       "ceccarelli2016_idh_status"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verhaak-core", type=Path,
                    default=Path("data/external/verhaak2010_gdc/TCGA_unified_CORE_ClaNC840.txt"))
    ap.add_argument("--brennan", type=Path,
                    default=Path("data/external/brennan2013_gdc/Molecular_subtype_classification.xlsx"))
    ap.add_argument("--ceccarelli", type=Path,
                    default=Path("data/external/ceccarelli2016_gdc/TableS1.PatientData.20151020.v3.xlsx"))
    ap.add_argument("--hgnc-rda", type=Path,
                    default=Path("data/external/hgnchelper_cran/x/HGNChelper/data/hgnc.table.rda"))
    ap.add_argument("--annotation", type=Path,
                    default=Path("data/raw/tcga_gbm/expression/gene_annotation.parquet"))
    ap.add_argument("--gene-map", type=Path, default=Path("configs/subtype_centroids/clanc840_genes.json"))
    ap.add_argument("--labels", type=Path,
                    default=Path("data/processed/tcga_gbm/published_subtype_labels.parquet"))
    a = ap.parse_args(argv)
    _refuse_cgga(a.verhaak_core, a.brennan, a.ceccarelli, a.hgnc_rda, a.annotation, a.gene_map, a.labels)

    import pyreadr

    genes, core = read_verhaak_core(a.verhaak_core)
    annotation = pd.read_parquet(a.annotation)
    hgnc = pyreadr.read_r(str(a.hgnc_rda))["hgnc.table"]
    resolved = resolve_genes(genes, annotation, hgnc)
    counts = pd.Series([r["resolution"] for r in resolved]).value_counts().to_dict()
    gene_doc = {
        "what": "The Verhaak et al. (2010) ClaNC 840-gene signature resolved to the "
                "GENCODE v36 symbols of gene_annotation.parquet. No patient data.",
        "rule": "direct if the published symbol is in the annotation; else the unique HGNC "
                "approved symbol (HGNChelper 0.8.16 hgnc.table) that is in the annotation; "
                "else excluded (ambiguous or not_in_annotation)",
        "inputs": {
            "verhaak_core": {"file": a.verhaak_core.name, "sha256": _sha256(a.verhaak_core)},
            "hgnc_table": {"file": "HGNChelper_0.8.16.tar.gz:HGNChelper/data/hgnc.table.rda",
                           "sha256": _sha256(a.hgnc_rda)},
            "annotation": {"file": str(a.annotation).replace("\\", "/"),
                           "sha256": _sha256(a.annotation)},
        },
        "counts": {k: int(v) for k, v in sorted(counts.items())},
        "n_usable": sum(1 for r in resolved if r["symbol"]),
        "genes": resolved,
    }
    a.gene_map.parent.mkdir(parents=True, exist_ok=True)
    a.gene_map.write_text(json.dumps(gene_doc, indent=1) + "\n", encoding="utf-8")

    br = read_brennan(a.brennan)
    ce = read_ceccarelli(a.ceccarelli)
    lab = br.merge(ce, on="patient_id", how="outer")
    lab["verhaak2010_core_label"] = lab["patient_id"].map(core)
    lab = lab.sort_values("patient_id").reset_index(drop=True)
    a.labels.parent.mkdir(parents=True, exist_ok=True)
    lab.to_parquet(a.labels, index=False)
    print(f"genes: {gene_doc['counts']}; usable {gene_doc['n_usable']} of 840")
    print(f"labels: {len(lab)} patients -> {a.labels}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
