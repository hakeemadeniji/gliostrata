"""Build and freeze the two Verhaak 840-gene centroid sets on TCGA (stage 3b, D6(b)).

  verhaak840_tcga_v2train  estimated on the v2 TRAINING patients only. Used for
                           every per-patient subtype feature and for any
                           internal evaluation on the v2 test split.
  verhaak840_tcga_all      estimated on every labelled TCGA patient. Reserved
                           for T1.

Scoring population: cohort patients with an expression sample (the columns of
tpm.parquet that are cohort patient ids, in cohort order, as in
scripts/build_tcga_feature_store.py). Training labels: the published Brennan et
al. (2013) expression subtype, kept when it is one of the four Verhaak classes.
"G-CIMP" (a separately reported group in that table) and "NA" are not used for
training; they are reported descriptively.

Genes: the 831 of the 840 published symbols that resolve in the annotation
(configs/subtype_centroids/clanc840_genes.json). Where a symbol has more than
one Ensembl id in the annotation, the id with the highest mean log2(TPM + 1)
over the v2 training patients of the scoring population is frozen, for both
sets, so no test patient influences any choice in the v2train set.

Centroid of a class = the mean, over its training patients, of the per-sample
transform of gliostrata.profiling.centroid (log2(TPM + 1), z-scored across the
frozen genes within each sample).

Writes the two set files and configs/subtype_centroids/CHECKSUMS.json. No
CGGA path is accepted; no survival outcome is read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.profiling.centroid import (
    CHECKSUM_FILE,
    SCALING,
    SET_ALL,
    SET_V2TRAIN,
    CentroidSet,
    centroid_set_to_json,
    sha256_bytes,
    transform_sample,
)
from gliostrata.profiling.verhaak import strip_ensembl_version

FOUR = ("Classical", "Mesenchymal", "Proneural", "Neural")
DEFAULTS = {
    "cohort": Path("data/raw/tcga_gbm/clinical_cohort.parquet"),
    "expression": Path("data/raw/tcga_gbm/expression/tpm.parquet"),
    "annotation": Path("data/raw/tcga_gbm/expression/gene_annotation.parquet"),
    "labels": Path("data/processed/tcga_gbm/published_subtype_labels.parquet"),
    "partition": Path("data/partitions/tcga_gbm_clinical_v2.json"),
    "gene_map": Path("configs/subtype_centroids/clanc840_genes.json"),
    "out_dir": Path("configs/subtype_centroids"),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def refuse_cgga(*paths) -> None:
    for p in paths:
        if "cgga" in str(p).lower():
            raise ValueError(f"refusing {p}: TCGA only")


def scoring_population(cohort: pd.DataFrame, expression: pd.DataFrame) -> list[str]:
    cols = {str(c) for c in expression.columns}
    ids = cohort["patient_id"].astype(str)
    return [p for p in ids if p in cols]


def freeze_gene_ids(
    symbols: list[str], annotation: pd.DataFrame, expression: pd.DataFrame, choose_on: list[str]
) -> tuple[list[str], list[str], dict]:
    """One version-stripped Ensembl id per symbol; duplicates by max mean over ``choose_on``."""
    ann = annotation.assign(stripped=strip_ensembl_version(annotation["gene_id"].astype(str)))
    expr_ids = pd.Index(strip_ensembl_version([str(i) for i in expression.index]))
    pos = pd.Series(np.arange(len(expr_ids)), index=expr_ids)
    genes, ens, dup = [], [], {}
    for s in symbols:
        cands = [e for e in ann.loc[ann["gene_name"] == s, "stripped"] if e in pos.index]
        if not cands:
            continue
        if len(cands) > 1:
            sub = np.log2(expression.iloc[[int(pos[e]) for e in cands]][choose_on].to_numpy(float) + 1)
            means = sub.mean(axis=1)
            pick = cands[int(np.argmax(means))]
            dup[s] = {"candidates": cands, "chosen": pick,
                      "rule": "highest mean log2(TPM+1) over the v2 training patients"}
        else:
            pick = cands[0]
        genes.append(s)
        ens.append(pick)
    return genes, ens, dup


def transformed_matrix(expression: pd.DataFrame, ens: list[str], genes: list[str],
                       patients: list[str]) -> pd.DataFrame:
    expr_ids = pd.Index(strip_ensembl_version([str(i) for i in expression.index]))
    pos = pd.Series(np.arange(len(expr_ids)), index=expr_ids)
    sub = expression.iloc[[int(pos[e]) for e in ens]][patients].to_numpy(dtype=float)
    z = np.column_stack([transform_sample(sub[:, j]) for j in range(sub.shape[1])])
    return pd.DataFrame(z, index=pd.Index(genes, name="gene"), columns=patients)


def make_set(set_id: str, z: pd.DataFrame, labels: pd.Series, genes, ens, meta: dict) -> CentroidSet:
    cen = pd.DataFrame({s: z[labels.index[labels == s]].mean(axis=1) for s in FOUR})
    return CentroidSet(set_id=set_id, genes=tuple(genes), subtypes=FOUR, centroids=cen,
                       scaling=SCALING, source_sha256="", meta=meta, ensembl_ids=tuple(ens))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", type=Path, default=v)
    a = ap.parse_args(argv)
    refuse_cgga(*vars(a).values())

    cohort = pd.read_parquet(a.cohort)
    expression = pd.read_parquet(a.expression)
    annotation = pd.read_parquet(a.annotation)
    labels = pd.read_parquet(a.labels).set_index("patient_id")
    part = json.loads(a.partition.read_text(encoding="utf-8"))
    train_ids, test_ids = set(part["train_ids"]), set(part["test_ids"])
    gene_map = json.loads(a.gene_map.read_text(encoding="utf-8"))
    symbols = [g["symbol"] for g in gene_map["genes"] if g["symbol"]]

    pop = scoring_population(cohort, expression)
    lab = labels["brennan2013_cluster"].reindex(pop)
    labelled = [p for p in pop if lab.get(p) in FOUR]
    train = [p for p in labelled if p in train_ids]
    test = [p for p in labelled if p in test_ids]
    pop_train = [p for p in pop if p in train_ids]

    genes, ens, dup = freeze_gene_ids(symbols, annotation, expression, pop_train)
    z = transformed_matrix(expression, ens, genes, labelled)

    def counts(ids):
        return {s: int((lab.reindex(ids) == s).sum()) for s in FOUR}

    base_meta = {
        "label_source": "Brennan et al. 2013 (Cell 155:462), Molecular_subtype_classification.xlsx, "
                        "column 'Cluster'; four Verhaak classes only (G-CIMP and NA not used)",
        "label_provenance": "mixed: 202 samples of Verhaak et al. 2010 carry the labels of its "
                            "consensus hierarchical clustering; the other 342 were assigned by a "
                            "single-sample GSEA classifier (Brennan 2013). Flag: verhaak2010_included",
        "scoring_population": "cohort patients with an expression sample",
        "n_scoring_population": len(pop),
        "n_scoring_population_by_split": {"train": len(pop_train),
                                          "test": sum(p in test_ids for p in pop)},
        "n_labelled_four_class": len(labelled),
        "n_labelled_by_split": {"train": len(train), "test": len(test)},
        "n_genes": len(genes),
        "gene_map": str(a.gene_map).replace("\\", "/"),
        "gene_map_sha256": sha256_file(a.gene_map),
        "duplicate_symbol_resolution": dup,
        "inputs_sha256": {
            "cohort": sha256_file(a.cohort),
            "expression": sha256_file(a.expression),
            "annotation": sha256_file(a.annotation),
            "labels": sha256_file(a.labels),
            "partition": sha256_file(a.partition),
        },
        "partition": part["partition_name"],
    }
    sets = {
        SET_V2TRAIN: (train, "v2 training patients only; for per-patient features and any "
                             "internal evaluation on the v2 test split"),
        SET_ALL: (train + test, "every labelled TCGA patient; reserved for T1"),
    }
    a.out_dir.mkdir(parents=True, exist_ok=True)
    reg_path = a.out_dir / CHECKSUM_FILE
    registry = json.loads(reg_path.read_text(encoding="utf-8")) if reg_path.is_file() else {}
    for set_id, (ids, use) in sets.items():
        meta = dict(base_meta, use=use, n_training=len(ids), n_training_by_class=counts(ids))
        cs = make_set(set_id, z[ids], lab.reindex(ids), genes, ens, meta)
        text = centroid_set_to_json(cs)
        path = a.out_dir / f"{set_id}.json"
        path.write_bytes(text.encode("utf-8"))
        registry[path.name] = sha256_bytes(text.encode("utf-8"))
        print(f"{set_id}: {len(ids)} training patients {counts(ids)}; {len(genes)} genes; "
              f"sha256 {registry[path.name]}")
    reg_path.write_text(json.dumps(dict(sorted(registry.items())), indent=1) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
