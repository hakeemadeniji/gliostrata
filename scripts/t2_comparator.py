"""T2 improvement row: choose the comparator gene set, then (after the deposit) evaluate it.

Pre-registration section 2.2 item 2, as amended 2026-11-27 (option B and the
batch-stratified C-index).

``select``: run before the deposit. It reads no CGGA value.
    1. One availability access to CGGA, declared to the read guard first. It reads the
       gene names of ``counts.parquet`` (its ``gene_name`` column only) and nothing else.
    2. The gene universe: symbols present in both TCGA (GENCODE v36 names, Ensembl ids
       summed per name) and CGGA. Rule R2 drops a gene set with any gene outside it.
    3. TCGA only: every eligible gene set is refitted and scored by 5-fold CV on the v2
       training-split patients with RNA-seq and an observed outcome. The best is chosen
       (ties: earlier publication, then lower PMID) and refitted on every TCGA patient
       with RNA-seq and an outcome.
    4. Writes the chosen model (terms, coefficients, universe checksum) to a config that
       is committed and named in the pre-registration before the deposit. It never
       overwrites one.

``evaluate``: run once, after the deposit. It is a value read of the sealed partition
``cgga_external_v1`` under the registered metric ``T2_cindex_improvement_vs_baseline``,
and the read guard refuses it before the deposit. It scores variant 2 (refitted on TCGA
exactly as ``scripts/validate_cox_cgga.py`` does) on the partition's complete cases
against two models:
  * the PRIMARY comparator, which decides the verdict (item 2 (j)): PMID 38174685 applied
    as published, on log2(FPKM + 1) from CGGA's RSEM FPKM release
    (``configs/t2_comparator/primary.json``);
  * the SECONDARY, descriptive comparator: option B's refitted gene set
    (``configs/t2_comparator/chosen.json``), reported without a verdict.
Each comparison is a paired bootstrap of the batch-stratified C-index difference.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.prediction.geneset_comparator import (
    GeneSet,
    choose_gene_set,
    collapse_to_symbols,
    design_matrix,
    fit_gene_set,
    log2_cpm,
    paired_bootstrap_difference,
    score_published,
)

ROOT = Path(__file__).resolve().parents[1]
GENE_SETS = Path("configs/t2_comparator/gene_sets.json")
CHOSEN = Path("configs/t2_comparator/chosen.json")
PRIMARY = Path("configs/t2_comparator/primary.json")
UNIVERSE = Path("configs/t2_comparator/gene_universe.txt")
SELECTION_OUT = Path("data/results/t2_comparator/selection.json")
EVAL_OUT = Path("data/results/t2_comparator/improvement_row.json")
LEDGER = Path("data/access_ledger.jsonl")

TCGA_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
TCGA_JOINED = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
TCGA_COUNTS = Path("data/raw/tcga_gbm/expression/counts.parquet")
TCGA_ANNOTATION = Path("data/raw/tcga_gbm/expression/gene_annotation.parquet")
DEV_PARTITION = Path("data/partitions/tcga_gbm_clinical_v2.json")

CGGA_COUNTS = Path("data/raw/cgga/counts.parquet")
CGGA_COHORT = Path("data/raw/cgga/clinical_cohort.parquet")
CGGA_ELIGIBLE = Path("data/processed/cgga/cgga_external_eligible.parquet")
EXT_PARTITION = Path("data/partitions/cgga_external_v1.json")

ANALYSIS = "T2_cindex_improvement_vs_baseline"
TIME, EVENT = "overall_survival_months", "vital_status"


def _validate_module():
    spec = importlib.util.spec_from_file_location(
        "validate_cox_cgga", ROOT / "scripts" / "validate_cox_cgga.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["validate_cox_cgga"] = mod
    spec.loader.exec_module(mod)
    return mod


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_gene_sets(path: Path = GENE_SETS) -> list[GeneSet]:
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    return [GeneSet(pmid=str(g["pmid"]), first_published=g["first_published"],
                    genes=tuple(g.get("genes", [])),
                    pairs=tuple(tuple(p) for p in g.get("pairs", [])),
                    source=g.get("source", ""), label=g.get("label", ""))
            for g in doc["gene_sets"]]


def tcga_outcomes(cohort: pd.DataFrame) -> pd.DataFrame:
    """Patients with an observed outcome (time > 0 and a vital status)."""
    frame = cohort[["patient_id", TIME, EVENT]].copy()
    frame = frame[frame[TIME].notna() & frame[EVENT].notna()]
    frame = frame[pd.to_numeric(frame[TIME]) > 0]
    frame[EVENT] = pd.to_numeric(frame[EVENT]).astype(int)
    return frame.set_index("patient_id")


def cgga_gene_names(ledger: Path, *, dry_run: bool = False) -> list[str]:
    """The gene names of CGGA's counts matrix, read after declaring the access."""
    from gliostrata.evaluation.read_guard import AVAILABILITY, Access, hash_inputs, record_access

    if not dry_run:
        record_access(Access(
            cohort="cgga", kind=AVAILABILITY, analysis="t2-comparator-gene-availability",
            purpose=("Rule R2 of pre-registration section 2.2 item 2 (option B): which gene "
                     "symbols exist in the CGGA RNA-seq count matrix. Reads only the "
                     "gene_name column of counts.parquet; no sample column, no expression "
                     "value, no clinical field, no outcome."),
            script="scripts/t2_comparator.py",
            inputs=hash_inputs([CGGA_COUNTS]),
        ), ledger_path=ledger)
    import pyarrow.parquet as pq

    names = pq.read_table(CGGA_COUNTS, columns=["gene_name"]).column("gene_name").to_pylist()
    return [str(n) for n in names]


def select(args) -> None:
    out_cfg = Path(args.chosen)
    if out_cfg.exists():
        raise SystemExit(f"{out_cfg} exists; a chosen comparator is never overwritten.")
    gene_sets = load_gene_sets(Path(args.gene_sets))
    print(f"{len(gene_sets)} gene sets passed E1, E2, E5 and R1 (from {args.gene_sets}).")

    tcga_counts = collapse_to_symbols(pd.read_parquet(TCGA_COUNTS), pd.read_parquet(TCGA_ANNOTATION))
    cgga_names = cgga_gene_names(Path(args.ledger))
    universe = sorted(set(tcga_counts.index) & set(cgga_names))
    print(f"Gene universe: {len(universe)} symbols in both TCGA ({len(tcga_counts)}) "
          f"and CGGA ({len(cgga_names)}).")

    eligible, dropped = [], []
    for gs in gene_sets:
        missing = gs.missing_from(universe)
        (dropped if missing else eligible).append((gs, missing))
    for gs, missing in dropped:
        print(f"  R2 drops PMID {gs.pmid}: not in both matrices: {missing}")
    eligible = [gs for gs, _ in eligible]
    if not eligible:
        raise SystemExit("No gene set passes R2: the row is NOT EVALUABLE.")

    expr = log2_cpm(tcga_counts, universe)
    outcomes = tcga_outcomes(pd.read_parquet(TCGA_COHORT))
    train_ids = set(json.loads(DEV_PARTITION.read_text(encoding="utf-8"))["train_ids"])
    with_rna = [p for p in expr.columns if p in outcomes.index]
    train = sorted(p for p in with_rna if p in train_ids)
    print(f"Choice set: {len(train)} v2-training patients with RNA-seq and an outcome.")
    table, chosen = choose_gene_set(
        eligible, expr, train, outcomes.loc[train, TIME].to_numpy(float),
        outcomes.loc[train, EVENT].to_numpy(int), n_boot=args.n_boot)
    print(table.to_string(index=False))

    final = sorted(with_rna)
    model = fit_gene_set(chosen, expr, final, outcomes.loc[final, TIME].to_numpy(float),
                         outcomes.loc[final, EVENT].to_numpy(int))
    UNIVERSE.parent.mkdir(parents=True, exist_ok=True)
    UNIVERSE.write_text("\n".join(universe) + "\n", encoding="utf-8", newline="\n")
    doc = {
        "what": "T2 improvement-row comparator, chosen under pre-registration 2.2 item 2 (option B).",
        "pmid": chosen.pmid, "first_published": chosen.first_published, "label": chosen.label,
        "genes": list(chosen.genes), "pairs": [list(p) for p in chosen.pairs],
        "terms": chosen.term_names, "source": chosen.source,
        "coefficients": {n: float(c) for n, c in zip(chosen.term_names, model.coef_)},
        "expression": "log2(CPM + 1) per sample from STAR read counts; library size over the gene universe",
        "gene_universe_file": UNIVERSE.as_posix(), "gene_universe_sha256": _sha256(UNIVERSE),
        "n_universe": len(universe),
        "final_fit": {"n": len(final), "n_events": int(model.n_events_),
                      "patients": "every TCGA patient with RNA-seq and an observed outcome"},
        "choice": {"n": len(train), "set": "tcga_gbm_clinical_v2 training split with RNA-seq and an outcome",
                   "cv": "5-fold stratified, seed 42, pooled out-of-fold Harrell's C-index",
                   "table": json.loads(table.to_json(orient="records"))},
        "r2_dropped": {gs.pmid: m for gs, m in dropped},
        "inputs_sha256": {p.as_posix(): _sha256(p) for p in
                          (TCGA_COUNTS, TCGA_ANNOTATION, TCGA_COHORT, DEV_PARTITION, Path(args.gene_sets))},
    }
    out_cfg.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")
    SELECTION_OUT.parent.mkdir(parents=True, exist_ok=True)
    SELECTION_OUT.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(f"\nChosen: PMID {chosen.pmid} ({len(chosen.term_names)} terms) -> {out_cfg}")


def read_fpkm(primary: dict, genes, samples) -> pd.DataFrame:
    """log2(FPKM + 1) for the given genes and samples, from CGGA's two FPKM archives."""
    import io
    import zipfile

    blocks = []
    for batch, src in primary["cgga_expression_source"].items():
        path = Path(src["path"])
        if _sha256(path) != src["sha256"]:
            raise SystemExit(f"{path}: checksum differs from the registered release")
        with zipfile.ZipFile(path) as zf:
            member = next(m for m in zf.namelist() if not m.endswith("/")
                          and not m.startswith("__MACOSX/") and not Path(m).name.startswith("._"))
            with zf.open(member) as handle:
                header = io.TextIOWrapper(handle, encoding="utf-8").readline().rstrip("\r\n").split("\t")
            wanted = [header[0]] + [c for c in header[1:] if c.strip() in set(samples)]
            with zf.open(member) as handle:
                frame = pd.read_csv(handle, sep="\t", usecols=wanted, index_col=0)
        frame.index = frame.index.astype(str).str.strip()
        blocks.append(frame.loc[list(genes)])
    fpkm = pd.concat(blocks, axis=1)
    missing = sorted(set(samples) - set(fpkm.columns))
    if missing:
        raise SystemExit(f"{len(missing)} partition patients have no FPKM column")
    return np.log2(fpkm[list(samples)].astype(float) + 1.0)


def evaluate(args) -> None:
    from gliostrata.evaluation.read_guard import VALUES, Access, hash_inputs, record_access
    from gliostrata.evaluation.thresholds import evaluate as verdict_of

    chosen = json.loads(Path(args.chosen).read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    if _sha256(UNIVERSE) != chosen["gene_universe_sha256"]:
        raise SystemExit("gene universe file does not match the chosen config")
    fpkm_paths = [Path(v["path"]) for v in primary["cgga_expression_source"].values()]
    inputs = [CGGA_COUNTS, CGGA_COHORT, CGGA_ELIGIBLE, Path(args.chosen), UNIVERSE, PRIMARY,
              *fpkm_paths]
    record_access(Access(
        cohort="cgga", kind=VALUES, analysis=ANALYSIS,
        purpose=("T2 improvement over the 2024 baseline (pre-registration 2.2 item 2, option B): "
                 "variant 2 against the refitted published gene set, batch-stratified C-index, "
                 "paired bootstrap, on cgga_external_v1."),
        script="scripts/t2_comparator.py", inputs=hash_inputs(inputs),
        partition=EXT_PARTITION.as_posix(),
    ), ledger_path=Path(args.ledger))

    v = _validate_module()
    tcga, cgga_rows = v.load_inputs(TCGA_COHORT, TCGA_JOINED, CGGA_COHORT)
    fitted = v.fit_variant(tcga, v.VARIANT_2)
    eligible = pd.read_parquet(CGGA_ELIGIBLE)
    ids = set(json.loads(EXT_PARTITION.read_text(encoding="utf-8"))["test_ids"])
    if set(eligible["patient_id"].astype(str)) != ids:
        raise SystemExit("eligibility table does not reproduce the sealed partition")
    cgga, _ = v.one_row_per_patient(cgga_rows)
    cgga = cgga[cgga["patient_id"].astype(str).isin(ids)]
    cgga = cgga.merge(eligible[["patient_id", "cgga_dataset"]].rename(
        columns={"cgga_dataset": "batch"}), on="patient_id", how="left", validate="one_to_one")
    frame, flow = v.complete_cases(cgga, v.VARIANT_2)
    x_ours, unseen = v.design_external(frame, fitted)
    frame, x_ours = frame.loc[~unseen].reset_index(drop=True), x_ours[~unseen]
    risk_ours = fitted.model.predict_log_hazard(x_ours)

    universe = UNIVERSE.read_text(encoding="utf-8").split()
    counts = pd.read_parquet(CGGA_COUNTS, columns=sorted(frame["patient_id"].astype(str)))
    expr = log2_cpm(counts, universe)
    gs = GeneSet(pmid=chosen["pmid"], first_published=chosen["first_published"],
                 genes=tuple(chosen["genes"]), pairs=tuple(tuple(p) for p in chosen["pairs"]))
    coef = np.array([chosen["coefficients"][t] for t in gs.term_names])
    risk_cmp = design_matrix(expr, gs, list(frame["patient_id"].astype(str))) @ coef

    pids = list(frame["patient_id"].astype(str))
    risk_primary = score_published(read_fpkm(primary, primary["coefficients"], pids),
                                   primary["coefficients"], pids)

    time_, event_ = frame[TIME].to_numpy(float), frame[EVENT].to_numpy(int)
    strata = frame["batch"].astype(str).to_numpy()
    main_res = paired_bootstrap_difference(time_, event_, risk_ours, risk_primary,
                                           strata=strata, n_boot=1000, seed=42)
    main_res["verdict"] = str(verdict_of(ANALYSIS, main_res["difference"],
                                         ci_low=main_res["ci_low"], ci_high=main_res["ci_high"]))
    main_res["comparator_pmid"] = primary["pmid"]
    secondary = paired_bootstrap_difference(time_, event_, risk_ours, risk_cmp,
                                            strata=strata, n_boot=1000, seed=42)
    secondary["comparator_pmid"] = chosen["pmid"]
    secondary["note"] = "descriptive only; carries no verdict (pre-registration 2.2 item 2 (j))"
    res = {"primary": main_res, "secondary_option_b": secondary,
           "flow": flow | {"n_unseen_category": int(unseen.sum())}}
    EVAL_OUT.parent.mkdir(parents=True, exist_ok=True)
    EVAL_OUT.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(res, indent=1))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select", help="choose the comparator (TCGA only; CGGA gene names)")
    s.add_argument("--gene-sets", default=str(GENE_SETS))
    s.add_argument("--chosen", default=str(CHOSEN))
    s.add_argument("--ledger", default=str(LEDGER))
    s.add_argument("--n-boot", type=int, default=1000)
    e = sub.add_parser("evaluate", help="the one CGGA read, after the deposit")
    e.add_argument("--chosen", default=str(CHOSEN))
    e.add_argument("--ledger", default=str(LEDGER))
    args = parser.parse_args(argv)
    (select if args.cmd == "select" else evaluate)(args)


if __name__ == "__main__":
    main()
