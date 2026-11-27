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
The read also reports the checks of pre-registration section 3.8. None of them changes a
verdict:
  * a proportional-hazards test on each fitted Cox model (variant 2 and option B, both on
    TCGA);
  * a negative control: each score against outcomes permuted within batch.

``dry-run``: the same evaluation, end to end, on synthetic CGGA-shaped files written to a
temporary folder. It runs before the real read, reads no CGGA file, records nothing in the
ledger, and writes its output to ``data/results/t2_comparator/dry_run.json``.
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
    permuted_outcome_control,
    score_published,
)
from gliostrata.prediction.ph_test import proportional_hazards_test

ROOT = Path(__file__).resolve().parents[1]
GENE_SETS = Path("configs/t2_comparator/gene_sets.json")
CHOSEN = Path("configs/t2_comparator/chosen.json")
PRIMARY = Path("configs/t2_comparator/primary.json")
UNIVERSE = Path("configs/t2_comparator/gene_universe.txt")
SELECTION_OUT = Path("data/results/t2_comparator/selection.json")
EVAL_OUT = Path("data/results/t2_comparator/improvement_row.json")
DRY_RUN_OUT = Path("data/results/t2_comparator/dry_run.json")
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


def read_fpkm(sources: dict, genes, samples, *, verify: bool = True) -> pd.DataFrame:
    """log2(FPKM + 1) for the given genes and samples, from CGGA's FPKM archives."""
    import io
    import zipfile

    blocks = []
    for src in sources.values():
        path = Path(src["path"])
        if verify and _sha256(path) != src["sha256"]:
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


def _ph_variant2(v, tcga, fitted) -> dict:
    """PH test of variant 2 on the TCGA rows it was fitted on (section 3.8)."""
    from gliostrata.prediction.cox import build_design_matrix

    fit_frame, _ = v.complete_cases(tcga, v.VARIANT_2)
    x_mat, time_, event_, report = build_design_matrix(fit_frame, v.VARIANT_2, impute="drop")
    if list(report.feature_names) != list(fitted.design_columns):
        raise SystemExit("variant 2 design differs from its fit")
    return proportional_hazards_test(x_mat, time_, event_, fitted.model.coef_,
                                     fitted.design_columns)


def _ph_option_b(chosen: dict, gs: GeneSet, universe) -> dict:
    """PH test of option B's refitted gene set on its TCGA fit (section 3.8)."""
    tcga_counts = collapse_to_symbols(pd.read_parquet(TCGA_COUNTS), pd.read_parquet(TCGA_ANNOTATION))
    expr = log2_cpm(tcga_counts, universe)
    outcomes = tcga_outcomes(pd.read_parquet(TCGA_COHORT))
    final = sorted(p for p in expr.columns if p in outcomes.index)
    if len(final) != chosen["final_fit"]["n"]:
        raise SystemExit("option B's TCGA fit set does not reproduce")
    coef = np.array([chosen["coefficients"][t] for t in gs.term_names])
    return proportional_hazards_test(design_matrix(expr, gs, final), outcomes.loc[final, TIME],
                                     outcomes.loc[final, EVENT], coef, gs.term_names)


def _verdict(verdict_of, res: dict) -> dict:
    """The registered verdict (section 6) for the improvement row.

    ``thresholds.evaluate`` requires the point estimate to lie inside its interval, which
    a percentile bootstrap does not guarantee. If it does not, the section 6 rule is
    applied directly to the same three numbers, and the output says so.
    """
    point, lo, hi = res["difference"], res["ci_low"], res["ci_high"]
    try:
        return verdict_of(ANALYSIS, point, lo, hi).to_dict()
    except ValueError as exc:
        status = ("PASS" if lo >= 0.03 else "INCONCLUSIVE" if point >= 0.03 else "FAIL")
        return {"metric_key": ANALYSIS, "status": status, "point_estimate": point,
                "ci_lower": lo, "ci_upper": hi, "threshold": 0.03,
                "reason": f"section 6 rule applied directly: {exc}"}


def run_evaluation(paths: dict, chosen: dict, primary: dict, *, verify_fpkm: bool = True) -> dict:
    """The registered evaluation. ``paths`` names the CGGA-side files; TCGA is the project's own."""
    from gliostrata.evaluation.thresholds import evaluate as verdict_of

    v = _validate_module()
    tcga, cgga_rows = v.load_inputs(TCGA_COHORT, TCGA_JOINED, paths["cohort"])
    fitted = v.fit_variant(tcga, v.VARIANT_2)
    eligible = pd.read_parquet(paths["eligible"])
    ids = set(json.loads(Path(paths["partition"]).read_text(encoding="utf-8"))["test_ids"])
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
    pids = list(frame["patient_id"].astype(str))

    universe = UNIVERSE.read_text(encoding="utf-8").split()
    import pyarrow.parquet as pq

    counts = (pq.read_table(paths["counts"], columns=["gene_name", *sorted(pids)])
              .to_pandas(ignore_metadata=True).set_index("gene_name"))
    expr = log2_cpm(counts, universe)
    gs = GeneSet(pmid=chosen["pmid"], first_published=chosen["first_published"],
                 genes=tuple(chosen["genes"]), pairs=tuple(tuple(p) for p in chosen["pairs"]))
    coef = np.array([chosen["coefficients"][t] for t in gs.term_names])
    risk_cmp = design_matrix(expr, gs, pids) @ coef
    risk_primary = score_published(read_fpkm(paths["fpkm"], primary["coefficients"], pids,
                                             verify=verify_fpkm),
                                   primary["coefficients"], pids)

    time_, event_ = frame[TIME].to_numpy(float), frame[EVENT].to_numpy(int)
    strata = frame["batch"].astype(str).to_numpy()
    main_res = paired_bootstrap_difference(time_, event_, risk_ours, risk_primary,
                                           strata=strata, n_boot=1000, seed=42)
    main_res["verdict"] = _verdict(verdict_of, main_res)
    main_res["comparator_pmid"] = primary["pmid"]
    secondary = paired_bootstrap_difference(time_, event_, risk_ours, risk_cmp,
                                            strata=strata, n_boot=1000, seed=42)
    secondary["comparator_pmid"] = chosen["pmid"]
    secondary["note"] = "descriptive only; carries no verdict (pre-registration 2.2 item 2 (j))"
    checks = {
        "negative_control": permuted_outcome_control(
            time_, event_, {"variant_2": risk_ours, "primary": risk_primary,
                            "secondary_option_b": risk_cmp}, strata, n_perm=1000, seed=42),
        "ph_test_variant_2_tcga": _ph_variant2(v, tcga, fitted),
        "ph_test_option_b_tcga": _ph_option_b(chosen, gs, universe),
        "note": "section 3.8: reported, never decisive",
    }
    return {"primary": main_res, "secondary_option_b": secondary, "checks": checks,
            "flow": flow | {"n_unseen_category": int(unseen.sum())}}


def evaluate(args) -> None:
    from gliostrata.evaluation.read_guard import VALUES, Access, hash_inputs, record_access

    chosen = json.loads(Path(args.chosen).read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    if _sha256(UNIVERSE) != chosen["gene_universe_sha256"]:
        raise SystemExit("gene universe file does not match the chosen config")
    if not DRY_RUN_OUT.exists():
        raise SystemExit(f"run the synthetic dry run first ({DRY_RUN_OUT}); section 3.8")
    fpkm_paths = [Path(v["path"]) for v in primary["cgga_expression_source"].values()]
    inputs = [CGGA_COUNTS, CGGA_COHORT, CGGA_ELIGIBLE, Path(args.chosen), UNIVERSE, PRIMARY,
              *fpkm_paths]
    record_access(Access(
        cohort="cgga", kind=VALUES, analysis=ANALYSIS,
        purpose=("T2 improvement over the 2024 baseline (pre-registration 2.2 item 2 (j)): "
                 "variant 2 against PMID 38174685 as published (primary) and option B's "
                 "refitted gene set (secondary, descriptive), batch-stratified C-index, paired "
                 "bootstrap, section 3.8 checks, on cgga_external_v1."),
        script="scripts/t2_comparator.py", inputs=hash_inputs(inputs),
        partition=EXT_PARTITION.as_posix(),
    ), ledger_path=Path(args.ledger))
    res = run_evaluation({"cohort": CGGA_COHORT, "eligible": CGGA_ELIGIBLE,
                          "partition": EXT_PARTITION, "counts": CGGA_COUNTS,
                          "fpkm": primary["cgga_expression_source"]}, chosen, primary)
    EVAL_OUT.parent.mkdir(parents=True, exist_ok=True)
    EVAL_OUT.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(res, indent=1))


def write_synthetic_cgga(folder: Path, chosen: dict, primary: dict, *, n: int = 60,
                         seed: int = 7) -> dict:
    """CGGA-shaped synthetic inputs, in the formats of the real files. No CGGA file is read."""
    import zipfile

    rng = np.random.default_rng(seed)
    folder.mkdir(parents=True, exist_ok=True)
    tcga = pd.read_parquet(TCGA_COHORT)
    ids = [f"SYN_{i:04d}" for i in range(n)]
    batch = np.where(np.arange(n) < n * 0.6, "mRNAseq_693", "mRNAseq_325")
    sex = rng.choice(tcga["sex"].dropna().astype(str).unique(), n)
    arm = rng.choice(tcga["treatment_arm"].dropna().astype(str).unique(), n)
    cohort = pd.DataFrame({"patient_id": ids, "age": rng.normal(52, 12, n).round(1),
                           "sex": sex, "treatment_arm": arm,
                           "overall_survival_months": rng.exponential(16, n).round(2) + 0.5,
                           "vital_status": rng.binomial(1, 0.8, n).astype(float),
                           "kps": np.nan, "idh_mutated": np.nan, "mgmt_methylated": np.nan})
    cohort.to_parquet(folder / "clinical_cohort.parquet", index=False)
    pd.DataFrame({"patient_id": ids, "cgga_dataset": batch}).to_parquet(
        folder / "eligible.parquet", index=False)
    (folder / "partition.json").write_text(json.dumps({"test_ids": ids}), encoding="utf-8")
    universe = UNIVERSE.read_text(encoding="utf-8").split()
    lam = rng.gamma(1.0, 30.0, size=(len(universe), 1))
    counts = pd.DataFrame(rng.poisson(lam * np.ones((1, n))), index=universe, columns=ids)
    counts.index.name = "gene_name"
    counts.reset_index().to_parquet(folder / "counts.parquet", index=False)
    fpkm_sources = {}
    genes = list(primary["coefficients"]) + ["GAPDH", "ACTB"]
    for b in ("mRNAseq_693", "mRNAseq_325"):
        cols = [i for i, bb in zip(ids, batch) if bb == b]
        mat = pd.DataFrame(rng.gamma(2.0, 10.0, size=(len(genes), len(cols))), index=genes, columns=cols)
        text = "gene_name\t" + "\t".join(cols) + "\n" + "\n".join(
            g + "\t" + "\t".join(f"{x:.4f}" for x in mat.loc[g]) for g in genes) + "\n"
        zpath = folder / f"CGGA.{b}.RSEM-genes.synthetic.txt.zip"
        with zipfile.ZipFile(zpath, "w") as zf:
            zf.writestr(f"CGGA.{b}.RSEM-genes.synthetic.txt", text)
        fpkm_sources[b] = {"path": str(zpath), "sha256": None}
    return {"cohort": folder / "clinical_cohort.parquet", "eligible": folder / "eligible.parquet",
            "partition": folder / "partition.json", "counts": folder / "counts.parquet",
            "fpkm": fpkm_sources}


def dry_run(args) -> None:
    import tempfile

    chosen = json.loads(Path(args.chosen).read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory() as tmp:
        paths = write_synthetic_cgga(Path(tmp), chosen, primary)
        res = run_evaluation(paths, chosen, primary, verify_fpkm=False)
    res["synthetic"] = True
    DRY_RUN_OUT.parent.mkdir(parents=True, exist_ok=True)
    DRY_RUN_OUT.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print(f"Synthetic dry run complete -> {DRY_RUN_OUT}")
    print(json.dumps({k: res[k] for k in ("primary", "secondary_option_b")}, indent=1)[:1500])


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
    d = sub.add_parser("dry-run", help="the evaluation on synthetic CGGA-shaped files")
    d.add_argument("--chosen", default=str(CHOSEN))
    args = parser.parse_args(argv)
    {"select": select, "evaluate": evaluate, "dry-run": dry_run}[args.cmd](args)


if __name__ == "__main__":
    main()
