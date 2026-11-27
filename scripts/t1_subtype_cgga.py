"""T1 subtype-accuracy row: single-sample Verhaak subtyping of the CGGA microarray patients.

Pre-registration section 2.2 item 9 (a) to (e), section 3.8 and section 6.

Scoring set: the CGGA microarray patients (batch mRNA_array_301) that carry a reference
Verhaak label, frozen as a scoring-set manifest (``test_ids``). Classifier: single-sample
nearest centroid against the frozen set ``verhaak840_tcga_all``
(``gliostrata.profiling.centroid``). Its ``classify_single_sample`` applies log2(TPM + 1)
internally, which microarray values must not go through, so this script implements the
microarray path with the same centroid set and the same Pearson / z-score logic:

  * genes are matched by exact symbol to the set's symbols, no alias resolution;
  * rows sharing a symbol are averaged;
  * values enter on the log2 scale: as distributed (``log2_as_distributed``) or as
    log2(x + 1) (``log2_plus_one``). The scale is read from
    ``configs/t1_subtype/input_scale.json``, decided from the source documentation before
    any value is read; the script refuses without it;
  * per sample, z = (x - mean) / sd over that sample's genes (ddof = 0); the call is the
    centroid with the highest Pearson correlation, an exact tie going to the earlier class.

Metric: accuracy = correct / labelled with an exact two-sided Clopper-Pearson 95% interval;
the lower bound decides against 0.88 through ``thresholds.evaluate``. Fewer than 30
labelled patients scored is NOT EVALUABLE. Section 3.8 negative control: accuracy of the
fixed calls against reference labels permuted 1000 times (seed 42).

``dry-run``: the full evaluation on synthetic files in the formats of the real inputs, in a
temporary folder. No guard, no CGGA file. Writes ``data/results/t1/subtype_dry_run.json``.
``evaluate``: refuses without the dry-run output and the scale config; declares the value
read to the read guard, then runs once and writes ``data/results/t1/subtype_accuracy.json``.

Patients without a reference label, with a label outside the four classes, or outside the
manifest are excluded and counted. A manifest patient with a label but no matrix column is
also excluded and counted (``n_labelled_without_expression``); the floor and the accuracy
use the patients scored.
"""

from __future__ import annotations

import argparse
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.profiling.centroid import CENTROID_DIR, SET_ALL, _pearson, load_centroid_set

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = "T1_subtype_accuracy_cgga"
SCALE_CONFIG = Path("configs/t1_subtype/input_scale.json")
CENTROID_FILE = CENTROID_DIR / f"{SET_ALL}.json"
DRY_RUN_OUT = Path("data/results/t1/subtype_dry_run.json")
EVAL_OUT = Path("data/results/t1/subtype_accuracy.json")
LEDGER = Path("data/access_ledger.jsonl")
BATCH = "mRNA_array_301"
SCALES = ("log2_as_distributed", "log2_plus_one")
MIN_LABELLED = 30
N_PERM, SEED = 1000, 42
LABEL_COL = "subtype"


# ---------------------------------------------------------------- inputs
def load_scale(path: Path = SCALE_CONFIG) -> str:
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"{path} is absent: the input scale must be decided from the "
                         "source documentation and recorded before any value is read")
    scale = json.loads(path.read_text(encoding="utf-8")).get("scale")
    if scale not in SCALES:
        raise SystemExit(f"{path}: 'scale' must be one of {SCALES}, got {scale!r}")
    return scale


def read_matrix(path: Path, samples=None) -> pd.DataFrame:
    """Gene-symbol x sample matrix from a zipped tab-delimited file; duplicate symbols averaged."""
    with zipfile.ZipFile(path) as zf:
        member = next(m for m in zf.namelist() if not m.endswith("/")
                      and not m.startswith("__MACOSX/") and not Path(m).name.startswith("._"))
        with zf.open(member) as handle:
            header = io.TextIOWrapper(handle, encoding="utf-8").readline().rstrip("\r\n").split("\t")
        keep = None if samples is None else set(samples)
        wanted = [header[0]] + [c for c in header[1:] if keep is None or c.strip() in keep]
        with zf.open(member) as handle:
            frame = pd.read_csv(handle, sep="\t", usecols=wanted, index_col=0)
    frame.index = frame.index.astype(str).str.strip()
    frame.columns = [str(c).strip() for c in frame.columns]
    return average_duplicates(frame)


def average_duplicates(frame: pd.DataFrame) -> pd.DataFrame:
    """Rows sharing a symbol are averaged (exact symbol, no alias resolution)."""
    return frame.astype(float).groupby(level=0, sort=False).mean()


def to_log2(frame: pd.DataFrame, scale: str) -> tuple[pd.DataFrame, bool]:
    """Values on the log2 scale, and whether the registered contingency applied.

    Contingency (pre-registration 2.2 item 9, amended 2026-11-27, before any value was
    read): linear intensities cannot be negative. So if ``log2_plus_one`` was decided
    but any value is below 0, the values are already on a log scale and are used as
    distributed. The output records that the contingency applied.
    """
    vals = frame.astype(float)
    if scale == "log2_as_distributed":
        return vals, False
    if scale == "log2_plus_one":
        if (vals.to_numpy() < 0).any():
            return vals, True
        return np.log2(vals + 1.0), False
    raise ValueError(f"unknown scale {scale!r}")


def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t" if path.suffix in (".tsv", ".txt") else ",")


# ---------------------------------------------------------------- classifier
def zscore(x: np.ndarray) -> np.ndarray:
    """Per-sample z over genes, population sd; the same rule as centroid.transform_sample."""
    x = np.asarray(x, dtype=float)
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - x.mean()) / sd


def classify_log2(expr_log2: pd.DataFrame, cs) -> pd.DataFrame:
    """Nearest frozen centroid per sample on log2-scale, symbol-indexed values.

    Each sample is scored on its own: nothing is computed across samples.
    """
    present = [g for g in cs.genes if g in set(expr_log2.index)]
    if len(present) < 2:
        raise ValueError(f"only {len(present)} of the {len(cs.genes)} frozen genes are in the input")
    vals = expr_log2.loc[present].to_numpy(dtype=float)
    cen_all = cs.centroids.loc[present].to_numpy(dtype=float)
    rows = []
    for j, sample in enumerate(expr_log2.columns):
        # Contingency (item 9, amended 2026-11-27): a frozen gene with a missing or
        # non-finite value is dropped for this sample only; the count is reported.
        ok = np.isfinite(vals[:, j])
        if ok.sum() < 2:
            rows.append({"patient_id": str(sample), "subtype": None, "n_genes_used": int(ok.sum())})
            continue
        z = zscore(vals[ok, j])
        cen = cen_all[ok]
        corr = np.array([_pearson(z, cen[:, k]) for k in range(cen.shape[1])])
        call = cs.subtypes[int(np.nanargmax(corr))] if np.isfinite(corr).any() else None
        rows.append({"patient_id": str(sample), "subtype": call, "n_genes_used": int(ok.sum())})
    return pd.DataFrame.from_records(rows, columns=["patient_id", "subtype", "n_genes_used"])


# ---------------------------------------------------------------- metric
def clopper_pearson(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    from scipy.stats import binomtest

    ci = binomtest(int(k), int(n)).proportion_ci(confidence_level=level, method="exact")
    return float(ci.low), float(ci.high)


def permutation_control(calls, labels, *, n_perm: int = N_PERM, seed: int = SEED) -> dict:
    """Accuracy of the fixed calls against reference labels permuted n_perm times (section 3.8)."""
    calls, labels = np.asarray(calls, dtype=object), np.asarray(labels, dtype=object)
    rng = np.random.default_rng(seed)
    acc = np.array([float(np.mean(calls == rng.permutation(labels))) for _ in range(n_perm)])
    return {"n_perm": n_perm, "seed": seed, "mean": float(acc.mean()),
            "range_95": [float(np.percentile(acc, 2.5)), float(np.percentile(acc, 97.5))],
            "note": "section 3.8: reported, never decisive; about 1/4 expected"}


def score_calls(calls: pd.Series, labels: pd.Series, subtypes) -> dict:
    """Accuracy, interval, confusion matrix, verdict and control for aligned calls and labels."""
    from gliostrata.evaluation.thresholds import evaluate as verdict_of

    n = int(len(labels))
    correct = int((calls.to_numpy(dtype=object) == labels.to_numpy(dtype=object)).sum())
    cols = [*subtypes, "no_call"]
    conf = {r: {c: 0 for c in cols} for r in subtypes}
    for ref, call in zip(labels, calls):
        conf[ref]["no_call" if call is None or pd.isna(call) else call] += 1
    res = {"n_scored": n, "n_correct": correct, "confusion_matrix": conf,
           "confusion_rows": "reference label", "confusion_columns": "called subtype"}
    if n < MIN_LABELLED:
        res["verdict"] = {"metric_key": ANALYSIS, "status": "NOT EVALUABLE",
                          "reason": f"{n} labelled patients scored; the floor is {MIN_LABELLED}"}
        res["accuracy"] = correct / n if n else None
        return res
    lo, hi = clopper_pearson(correct, n)
    point = correct / n
    res.update({"accuracy": point, "ci_lower": lo, "ci_upper": hi,
                "interval": "exact two-sided Clopper-Pearson 95%",
                "verdict": verdict_of(ANALYSIS, point, lo, hi).to_dict(),
                "negative_control": permutation_control(calls, labels)})
    return res


def run_evaluation(paths: dict, scale: str, *, label_col: str = LABEL_COL,
                   centroid_path: Path = CENTROID_FILE) -> dict:
    cs = load_centroid_set(centroid_path)
    ids = [str(i) for i in json.loads(Path(paths["manifest"]).read_text(encoding="utf-8"))["test_ids"]]
    clin = read_table(paths["clinical"])
    clin["patient_id"] = clin["patient_id"].astype(str).str.strip()
    # CGGA reuses ids across batches; the reference label lives on the microarray row.
    if "cgga_dataset" in clin.columns:
        clin = clin[clin["cgga_dataset"].astype(str) == BATCH]
    if clin["patient_id"].duplicated().any():
        raise SystemExit("duplicate patient ids within the microarray batch")
    in_manifest = clin[clin["patient_id"].isin(set(ids))]
    n_outside = int(len(clin) - len(in_manifest))
    lab = in_manifest[label_col].astype("string").str.strip()
    has_label = lab.notna() & (lab != "")
    known = has_label & lab.isin(cs.subtypes)
    labels = pd.Series(lab[known].to_numpy(), index=in_manifest.loc[known, "patient_id"].to_numpy())
    expr = read_matrix(Path(paths["matrix"]), samples=set(labels.index))
    have = [p for p in labels.index if p in expr.columns]
    labels = labels.loc[have]
    expr_log2, fallback = to_log2(expr[have], scale)
    called = classify_log2(expr_log2, cs).set_index("patient_id")
    res = score_calls(called.loc[have, "subtype"], labels, cs.subtypes)
    res.update({
        "n_labelled": int(known.sum()),
        "n_genes_used_min": int(called["n_genes_used"].min()) if len(called) else 0,
        "n_genes_used_max": int(called["n_genes_used"].max()) if len(called) else 0,
        "scale_contingency_applied": bool(fallback),
        "n_genes_frozen": len(cs.genes),
        "input_scale": scale, "centroid_set": cs.set_id, "centroid_sha256": cs.source_sha256,
        "exclusions": {
            "n_manifest_ids": len(ids),
            "n_outside_manifest_in_clinical_table": n_outside,
            "n_manifest_without_clinical_row": int(len(set(ids) - set(clin["patient_id"]))),
            "n_without_reference_label": int((~has_label).sum()),
            "n_label_not_a_class": int((has_label & ~known).sum()),
            "n_labelled_without_expression": int(known.sum() - len(have)),
        },
    })
    return res


# ---------------------------------------------------------------- commands
def _write(out: Path, res: dict) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")


def write_synthetic(folder: Path, scale: str, *, n: int = 60, seed: int = 7) -> dict:
    """Synthetic inputs in the formats of the real files. No CGGA file is read."""
    rng = np.random.default_rng(seed)
    folder.mkdir(parents=True, exist_ok=True)
    cs = load_centroid_set(CENTROID_FILE)
    ids = [f"SYN_{i:04d}" for i in range(n)]
    true = rng.choice(cs.subtypes, n)
    genes = list(cs.genes) + [f"EXTRA{i}" for i in range(40)]
    base = np.clip(rng.normal(8.0, 2.0, size=(len(genes), 1)), 3.0, None)
    x = np.tile(base, (1, n)) + rng.normal(0, 0.5, size=(len(genes), n))
    for j, s in enumerate(true):
        x[: len(cs.genes), j] += 1.5 * cs.centroids[s].to_numpy(dtype=float)
    raw = x if scale == "log2_as_distributed" else np.clip(np.exp2(x) - 1.0, 0.0, None)
    mat = pd.DataFrame(raw, index=genes, columns=ids)
    dup = mat.loc[[genes[0], genes[1]]] * 1.01
    mat = pd.concat([mat, dup])
    text = "Gene_Name\t" + "\t".join(ids) + "\n" + "".join(
        g + "\t" + "\t".join(f"{v:.5f}" for v in row) + "\n" for g, row in zip(mat.index, mat.to_numpy()))
    with zipfile.ZipFile(folder / "matrix.txt.zip", "w") as zf:
        zf.writestr("matrix.txt", text)
    labels = pd.Series(true, dtype=object)
    labels[::10] = None  # some patients without a reference label
    extra = [f"SYN_X{i}" for i in range(5)]  # outside the manifest
    pd.DataFrame({"patient_id": ids + extra, LABEL_COL: list(labels) + ["Classical"] * 5}
                 ).to_csv(folder / "clinical.csv", index=False)
    (folder / "manifest.json").write_text(json.dumps({"test_ids": ids}), encoding="utf-8")
    return {"matrix": folder / "matrix.txt.zip", "clinical": folder / "clinical.csv",
            "manifest": folder / "manifest.json"}


def dry_run(args) -> None:
    import tempfile

    scale = load_scale(Path(args.scale_config))
    with tempfile.TemporaryDirectory() as tmp:
        res = run_evaluation(write_synthetic(Path(tmp), scale), scale)
    res["synthetic"] = True
    _write(Path(args.out), res)
    print(f"Synthetic dry run complete -> {args.out}")
    print(json.dumps({k: res[k] for k in ("n_labelled", "n_scored", "n_genes_used_min", "scale_contingency_applied", "accuracy",
                                          "verdict")}, indent=1))


def evaluate(args) -> None:
    from gliostrata.evaluation.read_guard import VALUES, Access, hash_inputs, record_access

    dry = Path(args.dry_run_out)
    if not dry.exists():
        raise SystemExit(f"run the synthetic dry run first ({dry}); section 3.8")
    scale = load_scale(Path(args.scale_config))
    inputs = [Path(args.matrix), Path(args.clinical), Path(args.manifest), Path(args.scale_config),
              CENTROID_FILE]
    record_access(Access(
        cohort="cgga", kind=VALUES, analysis=ANALYSIS,
        purpose=("T1 subtype accuracy (pre-registration 2.2 item 9): single-sample nearest "
                 "centroid on the CGGA microarray patients with a reference Verhaak label, "
                 f"input scale {scale}, exact Clopper-Pearson interval, section 3.8 control."),
        script="scripts/t1_subtype_cgga.py", inputs=hash_inputs(inputs),
        scoring_rule=ANALYSIS,
    ), ledger_path=Path(args.ledger))
    res = run_evaluation({"matrix": args.matrix, "clinical": args.clinical,
                          "manifest": args.manifest}, scale, label_col=args.label_col)
    _write(Path(args.out), res)
    print(json.dumps(res, indent=1))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dry-run", help="the evaluation on synthetic files")
    d.add_argument("--scale-config", default=str(SCALE_CONFIG))
    d.add_argument("--out", default=str(DRY_RUN_OUT))
    e = sub.add_parser("evaluate", help="the one CGGA read, after the deposit")
    e.add_argument("--matrix", required=True, help="zipped tab-delimited microarray matrix")
    e.add_argument("--clinical", required=True, help="table with patient_id and the label")
    e.add_argument("--manifest", required=True, help="frozen scoring-set manifest (test_ids)")
    e.add_argument("--label-col", default=LABEL_COL)
    e.add_argument("--scale-config", default=str(SCALE_CONFIG))
    e.add_argument("--dry-run-out", default=str(DRY_RUN_OUT))
    e.add_argument("--out", default=str(EVAL_OUT))
    e.add_argument("--ledger", default=str(LEDGER))
    args = parser.parse_args(argv)
    {"dry-run": dry_run, "evaluate": evaluate}[args.cmd](args)


if __name__ == "__main__":
    main()
