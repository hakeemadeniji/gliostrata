"""Measure how much a patient's subtype call depends on who else is scored (stage 3b, D6(a)).

Two classifiers, the same experiments:
  cohort_relative  gliostrata.profiling.verhaak.classify_verhaak, defaults:
                   each signature score is standardised across the samples of
                   the call, so the call depends on the batch;
  single_sample    gliostrata.profiling.centroid.classify_single_sample with
                   the frozen set verhaak840_tcga_v2train.

Experiments, each compared with joint scoring of the 283 cohort patients that
have an expression sample:
  1. random halves: the 283 split at random into two halves (141 and 142),
     each half scored on its own; repeated N_REPEATS times with SEED;
  2. IDH-mutant samples excluded: the patients with an IDH call of 1 removed,
     the rest scored together;
  3. the 155 feature-store patients scored together (the stage 3 store rule:
     an expression sample, age, KPS and sex present, and both a MGMT and an
     IDH call);
  4. the 45 of those 155 in the v2 TEST split scored alone, compared with the
     283 joint calls and with the 155 joint calls. In stage 3 those 45 were
     scored jointly with the training patients, which is the leak the review
     found (M6).

SEED is new to this analysis and is not a registered seed: it seeds only the
random halves. Counts only are written; no outcome is read; no CGGA path is
accepted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_subtype_centroids import refuse_cgga, scoring_population, sha256_file  # noqa: E402

from gliostrata.profiling.centroid import (  # noqa: E402
    SET_V2TRAIN,
    classify_single_sample,
    load_centroid_set,
)
from gliostrata.profiling.verhaak import classify_verhaak  # noqa: E402

SEED = 20261127
N_REPEATS = 20
DEFAULTS = {
    "cohort": Path("data/raw/tcga_gbm/clinical_cohort.parquet"),
    "expression": Path("data/raw/tcga_gbm/expression/tpm.parquet"),
    "annotation": Path("data/raw/tcga_gbm/expression/gene_annotation.parquet"),
    "joined": Path("data/processed/tcga_gbm/molecular_calls_joined.parquet"),
    "partition": Path("data/partitions/tcga_gbm_clinical_v2.json"),
    "centroids": Path("configs/subtype_centroids") / f"{SET_V2TRAIN}.json",
    "out": Path("data/results/tcga_gbm/subtype_scoring_set_dependence.json"),
}


def store_population(cohort: pd.DataFrame, joined: pd.DataFrame, pop: list[str]) -> list[str]:
    """The stage 3 store rule, from availability only."""
    c = cohort.assign(patient_id=cohort["patient_id"].astype(str)).set_index("patient_id")
    j = joined.assign(patient_id=joined["patient_id"].astype(str)).set_index("patient_id")
    keep = []
    for p in pop:
        row = c.loc[p]
        if pd.isna(row["age"]) or pd.isna(row["kps"]) or pd.isna(row["sex"]):
            continue
        if p not in j.index or pd.isna(j.loc[p, "mgmt_methylated"]) or pd.isna(j.loc[p, "idh_mutated"]):
            continue
        keep.append(p)
    return keep


def calls(classify, expression, annotation, ids) -> pd.Series:
    out = classify(expression[ids], annotation)
    return out.set_index("patient_id")["subtype"]


def n_changed(a: pd.Series, b: pd.Series) -> int:
    common = a.index.intersection(b.index)
    return int((a.loc[common] != b.loc[common]).sum())


def run(classify, expression, annotation, pop, idh_mut, store, test45) -> dict:
    joint = calls(classify, expression, annotation, pop)
    rng = np.random.default_rng(SEED)
    halves = []
    for _ in range(N_REPEATS):
        perm = list(rng.permutation(pop))
        h1, h2 = perm[: len(pop) // 2], perm[len(pop) // 2:]
        split = pd.concat([calls(classify, expression, annotation, h1),
                           calls(classify, expression, annotation, h2)])
        halves.append(n_changed(joint, split))
    no_idh = [p for p in pop if p not in idh_mut]
    excl = calls(classify, expression, annotation, no_idh)
    joint155 = calls(classify, expression, annotation, store)
    alone45 = calls(classify, expression, annotation, test45)
    return {
        "joint_283_distribution": {k: int(v) for k, v in joint.value_counts().sort_index().items()},
        "random_halves": {
            "n_patients_compared": len(pop),
            "changed_per_repeat": halves,
            "changed_mean": float(np.mean(halves)),
            "changed_min": int(min(halves)),
            "changed_max": int(max(halves)),
        },
        "idh_mutant_excluded": {
            "n_idh_mutant_removed": len(pop) - len(no_idh),
            "n_patients_compared": len(no_idh),
            "changed": n_changed(joint, excl),
        },
        "store_155_joint_vs_283_joint": {"n_patients_compared": len(store),
                                         "changed": n_changed(joint, joint155)},
        "test_45_alone": {
            "n_patients": len(test45),
            "changed_vs_283_joint": n_changed(joint, alone45),
            "changed_vs_155_joint": n_changed(joint155, alone45),
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k}", type=Path, default=v)
    a = ap.parse_args(argv)
    refuse_cgga(*vars(a).values())
    cohort = pd.read_parquet(a.cohort)
    expression = pd.read_parquet(a.expression)
    annotation = pd.read_parquet(a.annotation)
    joined = pd.read_parquet(a.joined)
    part = json.loads(a.partition.read_text(encoding="utf-8"))
    test_ids = set(part["test_ids"])
    cs = load_centroid_set(a.centroids)

    pop = scoring_population(cohort, expression)
    j = joined.assign(patient_id=joined["patient_id"].astype(str)).set_index("patient_id")
    idh = j["idh_mutated"].astype("Float64").fillna(-1.0)
    idh_mut = {p for p in pop if p in idh.index and float(idh.loc[p]) == 1.0}
    store = store_population(cohort, joined, pop)
    test45 = [p for p in store if p in test_ids]

    def single(expr, ann):
        return classify_single_sample(expr, ann, cs)

    result = {
        "what": "Subtype calls that change when the scoring set changes; counts only",
        "seed": SEED,
        "seed_note": "new to this analysis; not a registered seed; seeds only the random halves",
        "n_repeats": N_REPEATS,
        "populations": {"expression_patients": len(pop), "store_rule": len(store),
                        "store_rule_in_v2_test": len(test45),
                        "store_rule_in_v2_train": len(store) - len(test45),
                        "idh_mutant_among_expression_patients": len(idh_mut)},
        "cohort_relative": run(classify_verhaak, expression, annotation, pop, idh_mut, store, test45),
        "single_sample": run(single, expression, annotation, pop, idh_mut, store, test45),
        "single_sample_centroid_set": {"set_id": cs.set_id, "sha256": cs.source_sha256},
        "inputs_sha256": {k: sha256_file(getattr(a, k)) for k in
                          ("cohort", "expression", "annotation", "joined", "partition")},
        "reads_cgga": False,
        "reads_outcome": False,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("populations", "cohort_relative", "single_sample")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
