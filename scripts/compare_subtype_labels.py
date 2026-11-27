"""Compare subtype calls with published TCGA-GBM labels (stage 3b, D6(c)). Descriptive only.

Reference labels: Brennan et al. (2013) expression subtype ("Cluster"), the
labels the centroids were trained on. They are of mixed provenance: the 202
samples of Verhaak et al. (2010) carry that paper's consensus-clustering labels
and the rest were assigned by a single-sample GSEA classifier (Brennan 2013).
This is a check on TCGA against TCGA-derived labels, not external validation.

Modes compared:
  cohort_relative          classify_verhaak, defaults, the 283 scored jointly.
                           Labelled cohort-relative: used only for cohort
                           distributions, never as a per-patient feature.
  single_sample_v2train    the frozen v2-training centroids. On v2 TEST
                           patients this is the honest number: none of them
                           helped estimate the centroids. On v2 training
                           patients it is resubstitution.
  single_sample_all        the frozen all-TCGA centroids (reserved for T1).
                           Resubstitution on every labelled patient.
Variant (the one allowed by D6(e)): three classes, Neural dropped, centroids
re-estimated on the same v2 training patients, evaluated on v2 test patients
whose label is Classical, Mesenchymal or Proneural. No tuning round was used:
every choice was fixed before any agreement was computed.

Agreement is over patients whose label is one of the four Verhaak classes;
"G-CIMP" is reported as its own row of each confusion matrix and is not
counted. Intervals: exact (Clopper-Pearson) 95%. No outcome is read; no CGGA
path is accepted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_subtype_centroids import (  # noqa: E402
    FOUR,
    make_set,
    refuse_cgga,
    scoring_population,
    sha256_file,
    transformed_matrix,
)

from gliostrata.profiling.centroid import (  # noqa: E402
    SET_ALL,
    SET_V2TRAIN,
    classify_single_sample,
    load_centroid_set,
)
from gliostrata.profiling.verhaak import classify_verhaak  # noqa: E402

THREE = ("Classical", "Mesenchymal", "Proneural")
DEFAULTS = {
    "cohort": Path("data/raw/tcga_gbm/clinical_cohort.parquet"),
    "expression": Path("data/raw/tcga_gbm/expression/tpm.parquet"),
    "annotation": Path("data/raw/tcga_gbm/expression/gene_annotation.parquet"),
    "joined": Path("data/processed/tcga_gbm/molecular_calls_joined.parquet"),
    "labels": Path("data/processed/tcga_gbm/published_subtype_labels.parquet"),
    "partition": Path("data/partitions/tcga_gbm_clinical_v2.json"),
    "centroid_dir": Path("configs/subtype_centroids"),
    "out": Path("data/results/tcga_gbm/subtype_label_comparison.json"),
}


def exact_ci(k: int, n: int) -> list[float]:
    if n == 0:
        return [float("nan"), float("nan")]
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [round(lo, 4), round(hi, 4)]


def agreement(calls: pd.Series, labels: pd.Series, ids: list[str], classes=FOUR) -> dict:
    ids = [p for p in ids if labels.get(p) in classes]
    k = int((calls.loc[ids] == labels.loc[ids]).sum())
    n = len(ids)
    return {"n": n, "agree": k, "rate": round(k / n, 4) if n else None, "exact_95ci": exact_ci(k, n)}


def confusion(calls: pd.Series, labels: pd.Series, ids: list[str], cols=FOUR) -> dict:
    rows = [r for r in (*FOUR, "G-CIMP") if any(labels.get(p) == r for p in ids)]
    out = {}
    for r in rows:
        sel = [p for p in ids if labels.get(p) == r]
        out[r] = {c: int((calls.loc[sel] == c).sum()) for c in cols}
    return out


def by_idh(calls: pd.Series, idh: pd.Series, ids: list[str], cols=FOUR) -> dict:
    groups = {"IDH-mutant": [p for p in ids if idh.get(p) == 1],
              "IDH-wildtype": [p for p in ids if idh.get(p) == 0],
              "no IDH call": [p for p in ids if idh.get(p) not in (0, 1)]}
    return {g: {"n": len(s), **{c: int((calls.loc[s] == c).sum()) for c in cols}}
            for g, s in groups.items()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", type=Path, default=v)
    a = ap.parse_args(argv)
    refuse_cgga(*vars(a).values())
    cohort = pd.read_parquet(a.cohort)
    expression = pd.read_parquet(a.expression)
    annotation = pd.read_parquet(a.annotation)
    lab_tab = pd.read_parquet(a.labels).set_index("patient_id")
    part = json.loads(a.partition.read_text(encoding="utf-8"))
    train_ids, test_ids = set(part["train_ids"]), set(part["test_ids"])
    j = pd.read_parquet(a.joined)
    idh = (j.assign(patient_id=j["patient_id"].astype(str)).set_index("patient_id")["idh_mutated"]
           .astype("Float64").fillna(-1.0).astype(float))

    pop = scoring_population(cohort, expression)
    labels = lab_tab["brennan2013_cluster"].reindex(pop)
    included = lab_tab["verhaak2010_included"].reindex(pop)
    train = [p for p in pop if p in train_ids]
    test = [p for p in pop if p in test_ids]

    cs_tr = load_centroid_set(a.centroid_dir / f"{SET_V2TRAIN}.json")
    cs_all = load_centroid_set(a.centroid_dir / f"{SET_ALL}.json")
    cr = classify_verhaak(expression[pop], annotation).set_index("patient_id")["subtype"]
    ss_tr = classify_single_sample(expression[pop], annotation, cs_tr).set_index("patient_id")["subtype"]
    ss_all = classify_single_sample(expression[pop], annotation, cs_all).set_index("patient_id")["subtype"]

    # The one variant: three classes, Neural dropped, same v2 training patients.
    tr3 = [p for p in train if labels.get(p) in THREE]
    z3 = transformed_matrix(expression, list(cs_tr.ensembl_ids), list(cs_tr.genes), tr3)
    v3 = make_set("variant_three_class_v2train", z3, labels.reindex(tr3), cs_tr.genes,
                  cs_tr.ensembl_ids, {})
    cen3 = v3.centroids[list(THREE)]
    v3 = type(v3)(set_id=v3.set_id, genes=v3.genes, subtypes=THREE, centroids=cen3,
                  scaling=v3.scaling, source_sha256="not frozen (variant, descriptive only)",
                  meta={}, ensembl_ids=v3.ensembl_ids)
    ss_3 = classify_single_sample(expression[pop], annotation, v3).set_index("patient_id")["subtype"]

    def provenance_split(calls, ids):
        inc = [p for p in ids if included.get(p) == "INCLUDED"]
        new = [p for p in ids if included.get(p) == "NOT"]
        return {"verhaak2010_clustering_label": agreement(calls, labels, inc),
                "brennan2013_ssgsea_label": agreement(calls, labels, new)}

    res = {
        "what": "Agreement of subtype calls with published TCGA labels; descriptive, not validation",
        "reference": "Brennan et al. 2013 'Cluster' (four Verhaak classes; G-CIMP reported separately)",
        "populations": {
            "scoring_population": len(pop),
            "by_split": {"train": len(train), "test": len(test)},
            "labelled_four_class": {"all": sum(labels.get(p) in FOUR for p in pop),
                                    "train": sum(labels.get(p) in FOUR for p in train),
                                    "test": sum(labels.get(p) in FOUR for p in test)},
            "label_counts": {str(k): int(v) for k, v in labels.fillna("no label").value_counts().items()},
        },
        "cohort_relative": {
            "note": "cohort-relative mode: scored jointly over the 283; descriptive distributions only",
            "agreement_all": agreement(cr, labels, pop),
            "agreement_v2_test": agreement(cr, labels, test),
            "confusion_all": confusion(cr, labels, pop),
            "calls_by_idh": by_idh(cr, idh, pop),
            "agreement_by_label_provenance_all": provenance_split(cr, pop),
        },
        "single_sample_v2train": {
            "centroid_set": {"set_id": cs_tr.set_id, "sha256": cs_tr.source_sha256},
            "agreement_v2_test_honest": agreement(ss_tr, labels, test),
            "confusion_v2_test_honest": confusion(ss_tr, labels, test),
            "agreement_v2_train_resubstitution": agreement(ss_tr, labels, train),
            "agreement_by_label_provenance_v2_test": provenance_split(ss_tr, test),
            "calls_by_idh_all_283": by_idh(ss_tr, idh, pop),
            "distribution_all_283": {k: int(v) for k, v in ss_tr.value_counts().sort_index().items()},
        },
        "single_sample_all": {
            "centroid_set": {"set_id": cs_all.set_id, "sha256": cs_all.source_sha256},
            "agreement_resubstitution_all": agreement(ss_all, labels, pop),
            "confusion_resubstitution_all": confusion(ss_all, labels, pop),
            "calls_by_idh_all_283": by_idh(ss_all, idh, pop),
        },
        "variant_three_class": {
            "definition": "Neural dropped; centroids of Classical, Mesenchymal and Proneural "
                          "re-estimated on the same v2 training patients; not frozen, not used",
            "n_training": len(tr3),
            "agreement_v2_test_honest_three_class_labels": agreement(ss_3, labels, test, THREE),
            "confusion_v2_test": confusion(ss_3, labels, test, THREE),
        },
        "tuning_rounds_used": 0,
        "label_source_checks": {},
        "inputs_sha256": {k: sha256_file(getattr(a, k)) for k in
                          ("cohort", "expression", "annotation", "joined", "labels", "partition")},
        "reads_cgga": False,
        "reads_outcome": False,
    }
    core = lab_tab["verhaak2010_core_label"]
    both = core.dropna().index.intersection(lab_tab.index[lab_tab["brennan2013_cluster"].isin([*FOUR, "G-CIMP"])])
    res["label_source_checks"]["verhaak2010_core_vs_brennan2013"] = {
        "n_core_with_brennan_label": int(len(both)),
        "agree": int((core.loc[both] == lab_tab.loc[both, "brennan2013_cluster"]).sum()),
        "disagreements": {f"{c} -> {b}": int(n) for (c, b), n in
                          pd.Series(list(zip(core.loc[both], lab_tab.loc[both, "brennan2013_cluster"])))
                          .loc[lambda s: s.map(lambda t: t[0] != t[1])].value_counts().items()},
    }
    orig = lab_tab["ceccarelli2016_original"]
    both2 = lab_tab.index[lab_tab["brennan2013_cluster"].isin([*FOUR, "G-CIMP"]) & orig.notna() & (orig != "NA")]
    res["label_source_checks"]["ceccarelli2016_original_vs_brennan2013"] = {
        "n_both_labelled": int(len(both2)),
        "agree": int((orig.loc[both2] == lab_tab.loc[both2, "brennan2013_cluster"]).sum()),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
