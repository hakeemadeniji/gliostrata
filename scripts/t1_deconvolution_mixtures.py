"""T1, state-deconvolution row: pseudo-bulk mixtures from the Neftel 2019 atlas.

Implements the specification in ``docs/t1_deconvolution_spec.md``, which was
committed before this script was first run on the atlas. Every constant below
is fixed there; changing one after the run is a deviation and must be reported
as such.

WHAT IT DOES
    1. Loads the Smart-seq2 matrix of GEO GSM3828672 (log2(TPM/10+1); converted
       back to linear TPM) and the per-cell annotation of the Curated Cancer
       Cell Atlas (3CA) release of Neftel et al. 2019.
    2. Keeps the malignant cells that carry all six Neftel meta-module scores
       and gives each one state: MES = max(MES1, MES2), NPC = max(NPC1, NPC2),
       AC, OPC; the state is the highest of the four.
    3. Draws 100 mixing-weight vectors from a flat Dirichlet over the four
       states, and for each one a pseudo-bulk mixture of 500 malignant cells
       (cell counts per state multinomial in the weights, cells drawn without
       replacement within a state; the mixture is the mean linear TPM).
    4. Runs the project's estimator,
       ``gliostrata.profiling.deconvolve.deconvolve_states_expression`` with
       its defaults, once, on the 100 mixtures together.
    5. Scores it as registered (pre-registration 2.2 items 1 and 10): for each
       state the Pearson r between estimated fraction and imposed weight across
       the 100 mixtures; the primary metric is the mean of the four; a
       percentile bootstrap over mixtures (1000 resamples) gives its 95%
       interval; the registry (``gliostrata.evaluation.thresholds``, key
       ``T1_state_deconvolution_pearson_r``) gives the verdict.

WHAT IT DOES NOT DO
    It reads no cohort: no TCGA and no CGGA file, and a path containing "cgga"
    is refused. It opens no network connection. It refuses to overwrite an
    existing result file: the registered analysis is run once.

Usage:
    python scripts/t1_deconvolution_mixtures.py --dry-run   # inputs and counts only
    python scripts/t1_deconvolution_mixtures.py             # the run, once
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import sys
import tarfile
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path

import numpy as np
import pandas as pd

# --- Fixed by docs/t1_deconvolution_spec.md ------------------------------------
SPEC_PATH = "docs/t1_deconvolution_spec.md"
METRIC_KEY = "T1_state_deconvolution_pearson_r"
SEED_MIXTURES = 131928
SEED_BOOTSTRAP = 42
N_MIXTURES = 100
CELLS_PER_MIXTURE = 500
DIRICHLET_ALPHA = (1.0, 1.0, 1.0, 1.0)
N_BOOTSTRAP = 1000
CI_LEVEL = 0.95
#: More undefined bootstrap resamples than this (1% of 1000) make the row NOT EVALUABLE.
MAX_UNDEFINED_RESAMPLES = 10
#: State order everywhere: the estimator's column order.
STATES = ("NPC", "OPC", "AC", "MES")
STATE_COLUMNS = {"NPC": "state_npc", "OPC": "state_opc", "AC": "state_ac", "MES": "state_mes"}
MODULE_SCORES = ("MESlike2", "MESlike1", "AClike", "OPClike", "NPClike1", "NPClike2")
#: The estimator's options for the one pre-declared secondary (descriptive) run.
SECONDARY_FRACTION_METHOD = "proportional"

DEFAULT_EXPRESSION = Path(
    "data/external/neftel2019_geo/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
DEFAULT_ANNOTATION = Path("data/external/neftel2019_3ca/Meta-data_Neftel2019_Brain.tar.gz")
ANNOTATION_MEMBER = "Meta-data_Neftel2019_Brain/SmartSeq2/Cells.csv"
#: SHA-256 of the two inputs as downloaded on 2026-11-27
#: (data/provenance/neftel2019_sources.json). A mismatch stops the run.
EXPECTED_SHA256 = {
    "expression": "73be3319fb0b026dcae9dc33582ced9fef00a8e6470f63df3a7cbbf7daf960f5",
    "annotation": "df01ad311ddf152e08c7ebb3c2c6c9fcf22b119f5a4a9744aa0e81a37c7a7b80",
}
DEFAULT_OUT = Path("data/results/t1/neftel_mixture_deconvolution.json")
DEFAULT_TABLE = Path("data/results/t1/neftel_mixture_deconvolution_mixtures.parquet")


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def refuse_cohort(path: Path | str) -> None:
    """This analysis reads the atlas only; the held-out cohort is never read here."""
    if "cgga" in str(path).lower():
        raise ValueError(f"refusing {path}: T1's deconvolution row reads no cohort")


# --- 1. Inputs -------------------------------------------------------------------
def read_cell_annotation(tar_path: Path | str) -> pd.DataFrame:
    """The 3CA Smart-seq2 cell table, read from the archive in memory."""
    with tarfile.open(tar_path) as archive:
        member = archive.extractfile(ANNOTATION_MEMBER)
        if member is None:
            raise ValueError(f"{tar_path}: {ANNOTATION_MEMBER} is not a file")
        return pd.read_csv(io.BytesIO(member.read()))


def assign_states(cells: pd.DataFrame) -> pd.DataFrame:
    """Eligible cells and their state (spec section 3).

    Eligible: ``malignant == "yes"`` and all six meta-module scores present.
    State: the highest of MES = max(MESlike1, MESlike2), NPC = max(NPClike1,
    NPClike2), AC = AClike, OPC = OPClike. A tie for the highest is not
    assigned and is counted.
    """
    missing = [c for c in ("cell_name", "sample", "malignant", *MODULE_SCORES)
               if c not in cells.columns]
    if missing:
        raise ValueError(f"cell annotation lacks column(s) {missing}")
    has_scores = cells[list(MODULE_SCORES)].notna().all(axis=1)
    eligible = cells[(cells["malignant"].astype(str) == "yes") & has_scores].copy()
    four = pd.DataFrame({
        "NPC": eligible[["NPClike1", "NPClike2"]].max(axis=1),
        "OPC": eligible["OPClike"],
        "AC": eligible["AClike"],
        "MES": eligible[["MESlike1", "MESlike2"]].max(axis=1),
    }, index=eligible.index)[list(STATES)]
    values = four.to_numpy(dtype=float)
    top = values.max(axis=1, keepdims=True)
    n_top = (values == top).sum(axis=1)
    eligible["state"] = [STATES[i] for i in values.argmax(axis=1)]
    eligible["state_tied"] = n_top > 1
    out = eligible[~eligible["state_tied"]][["cell_name", "sample", "state"]]
    out.attrs["n_annotated"] = int(len(cells))
    out.attrs["n_malignant_with_scores"] = int(len(eligible))
    out.attrs["n_tied"] = int(eligible["state_tied"].sum())
    return out.reset_index(drop=True)


def read_linear_tpm(path: Path | str, cells: list[str]) -> pd.DataFrame:
    """Genes x cells, linear TPM (10 * (2**x - 1) of the GEO log2(TPM/10+1))."""
    log_tpm = pd.read_csv(path, sep="\t", index_col=0, usecols=["GENE", *cells],
                          dtype={c: np.float32 for c in cells})
    if log_tpm.index.has_duplicates:
        raise ValueError("expression matrix has duplicate gene symbols")
    missing = sorted(set(cells) - set(log_tpm.columns))
    if missing:
        raise ValueError(f"{len(missing)} eligible cell(s) are not in the matrix")
    log_tpm = log_tpm[cells]
    # One float32 array, transformed in place: the matrix is ~0.6 GB per copy.
    values = log_tpm.to_numpy(dtype=np.float32, copy=True)
    index = log_tpm.index
    del log_tpm
    np.exp2(values, out=values)
    values -= np.float32(1.0)
    values *= np.float32(10.0)
    return pd.DataFrame(values, index=index, columns=cells, copy=False)


# --- 2. Mixtures -------------------------------------------------------------------
def draw_weights(rng: np.random.Generator, n: int = N_MIXTURES) -> np.ndarray:
    """``n x 4`` imposed mixing weights, flat Dirichlet, state order ``STATES``."""
    return rng.dirichlet(np.asarray(DIRICHLET_ALPHA, dtype=float), size=n)


def make_mixtures(tpm: pd.DataFrame, cell_states: pd.DataFrame, weights: np.ndarray,
                  rng: np.random.Generator,
                  cells_per_mixture: int = CELLS_PER_MIXTURE,
                  ) -> tuple[pd.DataFrame, np.ndarray]:
    """Pseudo-bulk mixtures (genes x mixtures, linear TPM) and the realised
    state fractions (mixtures x 4).

    For mixture i: counts ~ Multinomial(cells_per_mixture, weights[i]); for each
    state, that many cells drawn without replacement from the state's pool; the
    mixture is the mean of their linear TPM vectors. Draws are made in mixture
    order, and within a mixture in state order, from the one generator.
    """
    pools = {s: cell_states.loc[cell_states["state"] == s, "cell_name"].to_numpy()
             for s in STATES}
    col_of = {c: i for i, c in enumerate(tpm.columns)}
    values = tpm.to_numpy(dtype=np.float32, copy=False)
    mixtures = np.empty((values.shape[0], len(weights)), dtype=np.float64)
    realised = np.empty((len(weights), len(STATES)), dtype=float)
    for i, w in enumerate(weights):
        counts = rng.multinomial(cells_per_mixture, w)
        chosen: list[int] = []
        for s, k in zip(STATES, counts, strict=True):
            if k > len(pools[s]):
                raise ValueError(f"state {s}: {k} cells needed, {len(pools[s])} in the pool")
            picked = rng.choice(pools[s], size=int(k), replace=False) if k else []
            chosen.extend(col_of[c] for c in picked)
        mixtures[:, i] = values[:, chosen].mean(axis=1, dtype=np.float64)
        realised[i] = counts / float(cells_per_mixture)
    names = [f"mix{i + 1:03d}" for i in range(len(weights))]
    return pd.DataFrame(mixtures, index=tpm.index, columns=names), realised


# --- 3. Scoring ----------------------------------------------------------------------
def per_state_r(estimated: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Pearson r per state (column), across mixtures (rows). NaN if undefined."""
    out = np.full(estimated.shape[1], np.nan)
    for j in range(estimated.shape[1]):
        e, t = estimated[:, j], truth[:, j]
        if np.std(e) > 0 and np.std(t) > 0:
            out[j] = float(np.corrcoef(e, t)[0, 1])
    return out


def bootstrap(estimated: np.ndarray, truth: np.ndarray, rng: np.random.Generator,
              n_boot: int = N_BOOTSTRAP) -> dict:
    """Percentile bootstrap over mixtures of the mean and per-state r.

    The estimates are held fixed (the estimator is not re-run in a resample).
    A resample in which any per-state r is undefined is dropped and counted;
    if more than ``MAX_UNDEFINED_RESAMPLES`` are dropped, the intervals are not
    computed (``None``) and the row is NOT EVALUABLE (spec section 5).
    """
    n = estimated.shape[0]
    means, states, undefined = [], [], 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r = per_state_r(estimated[idx], truth[idx])
        if np.isnan(r).any():
            undefined += 1
            continue
        means.append(float(r.mean()))
        states.append(r)
    lo_q, hi_q = (1 - CI_LEVEL) / 2, 1 - (1 - CI_LEVEL) / 2

    def _ci(values) -> tuple[float, float, bool]:
        lo, hi = (float(v) for v in np.quantile(values, [lo_q, hi_q]))
        clipped = lo < -1 or hi > 1
        return max(lo, -1.0), min(hi, 1.0), clipped

    if undefined > MAX_UNDEFINED_RESAMPLES:
        return {"n_resamples": n_boot, "n_undefined_dropped": undefined,
                "mean_r": None, "per_state": None}
    mean_lo, mean_hi, mean_clip = _ci(means)
    per = np.asarray(states)
    per_ci = {s: _ci(per[:, j]) for j, s in enumerate(STATES)}
    return {"n_resamples": n_boot, "n_undefined_dropped": undefined,
            "mean_r": {"lower": mean_lo, "upper": mean_hi, "clipped": mean_clip},
            "per_state": {s: {"lower": v[0], "upper": v[1], "clipped": v[2]}
                          for s, v in per_ci.items()}}


def registry_verdict(point: float, lower: float, upper: float) -> dict:
    """The registry's verdict. A point outside its percentile interval is
    handled by widening the interval to include it (spec section 5), which can
    only lower the decisive lower bound."""
    from gliostrata.evaluation.thresholds import evaluate

    widened = not (lower <= point <= upper)
    lower, upper = float(lower), float(upper)
    lo, hi = min(lower, point), max(upper, point)
    verdict = evaluate(METRIC_KEY, point, lo, hi).to_dict()
    verdict["interval_widened_to_include_point"] = widened
    return verdict


def final_verdict(primary: dict) -> dict:
    """The registry's verdict, or NOT EVALUABLE where spec section 5 says so:
    a per-state r undefined on the full set, or more than
    ``MAX_UNDEFINED_RESAMPLES`` undefined bootstrap resamples."""
    if primary["mean_r"] is None:
        return {"metric_key": METRIC_KEY, "status": "NOT EVALUABLE",
                "reason": "a per-state Pearson r is undefined on the 100 mixtures "
                          f"(per-state r: {primary['per_state_r']})"}
    boot = primary["bootstrap"]
    if boot["mean_r"] is None:
        return {"metric_key": METRIC_KEY, "status": "NOT EVALUABLE",
                "reason": f"{boot['n_undefined_dropped']} of {boot['n_resamples']} "
                          f"bootstrap resamples have an undefined r (limit "
                          f"{MAX_UNDEFINED_RESAMPLES})"}
    return registry_verdict(primary["mean_r"], boot["mean_r"]["lower"], boot["mean_r"]["upper"])


def estimate(mixtures: pd.DataFrame, *, fraction_method: str | None = None) -> np.ndarray:
    """The project's estimator on the mixtures together; mixtures x 4 in ``STATES`` order."""
    from gliostrata.profiling.deconvolve import deconvolve_states_expression

    kwargs = {} if fraction_method is None else {"fraction_method": fraction_method}
    out = deconvolve_states_expression(mixtures, None, **kwargs)
    out = out.set_index("patient_id").loc[list(mixtures.columns)]
    return out[[STATE_COLUMNS[s] for s in STATES]].to_numpy(dtype=float)


def score(estimated: np.ndarray, weights: np.ndarray, realised: np.ndarray,
          rng: np.random.Generator) -> dict:
    """Primary metric, its interval and the descriptive numbers.

    ``mean_r`` is ``None`` when a per-state r is undefined on the full set; the
    row is then NOT EVALUABLE (spec section 5) and no interval is computed.
    """
    r = per_state_r(estimated, weights)
    defined = not np.isnan(r).any()
    point = float(r.mean()) if defined else None
    boot = bootstrap(estimated, weights, rng) if defined else None
    pooled = float(np.corrcoef(estimated.ravel(), weights.ravel())[0, 1])
    r_real = per_state_r(estimated, realised)
    return {
        "mean_r": point,
        "per_state_r": {s: (None if np.isnan(v) else float(v))
                        for s, v in zip(STATES, r, strict=True)},
        "bootstrap": boot,
        "pooled_r_descriptive": pooled,
        "per_state_r_against_realised_fractions": {
            s: (None if np.isnan(v) else float(v)) for s, v in zip(STATES, r_real, strict=True)},
        "rmse_per_state": {s: float(np.sqrt(np.mean((estimated[:, j] - weights[:, j]) ** 2)))
                           for j, s in enumerate(STATES)},
    }


def _software() -> dict:
    try:
        gs = package_version("gliostrata")
    except PackageNotFoundError:
        gs = "unknown"
    return {"gliostrata": gs, "python": platform.python_version(),
            "numpy": np.__version__, "pandas": pd.__version__}


def run(expression: Path, annotation: Path, out: Path, table: Path, *,
        dry_run: bool = False, check_checksums: bool = True) -> dict:
    for p in (expression, annotation, out, table):
        refuse_cohort(p)
    if not dry_run and out.exists():
        raise FileExistsError(
            f"{out} exists. The registered analysis is run once; a re-run is a new "
            "analysis and must be reported as one (spec section 7)")
    inputs = {}
    for name, path in (("expression", expression), ("annotation", annotation)):
        digest = sha256_file(path)
        if check_checksums and digest != EXPECTED_SHA256[name]:
            raise ValueError(f"{name}: {path} has SHA-256 {digest}, expected "
                             f"{EXPECTED_SHA256[name]}")
        inputs[name] = {"path": str(path).replace("\\", "/"),
                        "bytes": Path(path).stat().st_size, "sha256": digest}

    cell_states = assign_states(read_cell_annotation(annotation))
    counts = {
        "n_annotated_cells": cell_states.attrs["n_annotated"],
        "n_malignant_with_all_scores": cell_states.attrs["n_malignant_with_scores"],
        "n_tied_not_assigned": cell_states.attrs["n_tied"],
        "n_cells_used": int(len(cell_states)),
        "cells_per_state": {s: int((cell_states["state"] == s).sum()) for s in STATES},
        "n_tumours": int(cell_states["sample"].nunique()),
    }
    report = {"what": "T1 state-deconvolution row: Neftel 2019 pseudo-bulk mixtures",
              "spec": SPEC_PATH, "metric_key": METRIC_KEY,
              "run_before_deposit": True, "inputs": inputs, "counts": counts,
              "constants": {"seed_mixtures": SEED_MIXTURES, "seed_bootstrap": SEED_BOOTSTRAP,
                            "n_mixtures": N_MIXTURES, "cells_per_mixture": CELLS_PER_MIXTURE,
                            "dirichlet_alpha": list(DIRICHLET_ALPHA),
                            "n_bootstrap": N_BOOTSTRAP, "ci_level": CI_LEVEL,
                            "states": list(STATES)},
              "software": _software(), "dry_run": dry_run}
    if dry_run:
        return report

    tpm = read_linear_tpm(expression, cell_states["cell_name"].tolist())
    sums = tpm.sum(axis=0).to_numpy()
    counts["n_genes"] = int(tpm.shape[0])
    counts["cell_tpm_sum_min_max"] = [float(sums.min()), float(sums.max())]
    from gliostrata.profiling.deconvolve import NEFTEL_SIGNATURES
    counts["markers_in_atlas"] = {k: [int(sum(g in tpm.index for g in v)), len(v)]
                                  for k, v in NEFTEL_SIGNATURES.items()}

    rng_mix = np.random.default_rng(SEED_MIXTURES)
    weights = draw_weights(rng_mix)
    mixtures, realised = make_mixtures(tpm, cell_states, weights, rng_mix)
    estimated = estimate(mixtures)
    primary = score(estimated, weights, realised, np.random.default_rng(SEED_BOOTSTRAP))
    verdict = final_verdict(primary)
    secondary_est = estimate(mixtures, fraction_method=SECONDARY_FRACTION_METHOD)
    secondary = score(secondary_est, weights, realised, np.random.default_rng(SEED_BOOTSTRAP))

    report.update({
        "primary": {"estimator": "deconvolve_states_expression, defaults "
                                 "(ssgsea, softmax_z, temperature 1.0)", **primary},
        "verdict": verdict,
        "secondary_descriptive": {
            "estimator": f"deconvolve_states_expression, fraction_method="
                         f"{SECONDARY_FRACTION_METHOD!r}; pre-declared, descriptive, "
                         "cannot change the verdict", **secondary},
        "mixture_checks": {"mixture_tpm_sum_min_max": [float(mixtures.sum(axis=0).min()),
                                                       float(mixtures.sum(axis=0).max())],
                           "weights_sum_max_abs_dev": float(np.abs(weights.sum(1) - 1).max())},
    })
    table_frame = pd.DataFrame({"mixture": list(mixtures.columns)})
    for j, s in enumerate(STATES):
        table_frame[f"weight_{s}"] = weights[:, j]
        table_frame[f"realised_{s}"] = realised[:, j]
        table_frame[f"estimate_{s}"] = estimated[:, j]
        table_frame[f"estimate_proportional_{s}"] = secondary_est[:, j]
    table.parent.mkdir(parents=True, exist_ok=True)
    table_frame.to_parquet(table, index=False)
    report["outputs"] = {"mixture_table": {"path": str(table).replace("\\", "/"),
                                           "sha256": sha256_file(table)}}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8", newline="\n")
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION)
    ap.add_argument("--annotation", type=Path, default=DEFAULT_ANNOTATION)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    ap.add_argument("--dry-run", action="store_true",
                    help="check the inputs and print the cell counts; no mixture is made")
    args = ap.parse_args(argv)
    rep = run(args.expression, args.annotation, args.out, args.table, dry_run=args.dry_run)
    print(json.dumps({k: rep[k] for k in ("inputs", "counts")}, indent=2))
    if not args.dry_run:
        p, v = rep["primary"], rep["verdict"]
        if v["status"] == "NOT EVALUABLE":
            print(f"verdict NOT EVALUABLE: {v['reason']} (run before deposit)")
        else:
            lo, hi = p["bootstrap"]["mean_r"]["lower"], p["bootstrap"]["mean_r"]["upper"]
            print(f"mean per-state r = {p['mean_r']:.4f} (95% CI {lo:.4f} to {hi:.4f}); "
                  f"verdict {v['status']} (run before deposit)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
