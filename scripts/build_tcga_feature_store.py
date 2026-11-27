"""
Real-data vertical slice, L0 to L3, on the development cohort (TCGA) only.

WHY THIS EXISTS
    The L2 (profiling) and L3 (feature store) code had only ever run on
    synthetic data. This script runs the committed modules on the real TCGA-GBM
    inputs already on disk and writes a feature store with measured values, so
    that the first three layers are shown to work on real data before any model
    is fitted.

WHAT IT DOES
    1. Loads the clinical cohort (L0), the gene x sample TPM matrix, the gene
       annotation and the joined MGMT / IDH calls.
    2. Restricts the expression matrix to the cohort patients that have an
       expression sample. That set is the scoring set.
    3. Calls, unchanged and with their defaults:
           gliostrata.profiling.verhaak.classify_verhaak
           gliostrata.profiling.deconvolve.deconvolve_states_expression
           gliostrata.profiling.deconvolve.deconvolve_tme_expression
    4. Joins the molecular calls onto the cohort and assembles the store with
           gliostrata.features.store.build_feature_store
    5. Writes the store (--out) and a report of counts, missingness, feature
       ranges and input checksums (--report).

    The harmonised cohort carries all-null placeholder columns for the subtype,
    state, TME and molecular values: its schema was fixed before they existed.
    They are dropped before assembly, because ``build_feature_store`` keeps the
    cohort's column when a name collides with an L2 output. A cohort column
    that carries values under one of those names stops the script.

WHAT IS MEASURED AND WHAT IS NOT
    The 19 feature values are computed from real data. They are not validated:
    the subtype label is a within-cohort relative ranking, so it depends on
    which samples are in the scoring set, and the state and TME "fractions" are
    normalised marker-enrichment scores, not cell proportions. Nothing here
    addresses protocol test T1.

PATIENTS WHO DO NOT ENTER THE STORE, AND WHY
    - No expression sample: excluded, counted in the report.
    - No MGMT or IDH call. The default policy is to exclude the patient and
      count it (--missing-molecular drop). ``zero`` keeps the patient and
      records an absent call as 0, which says "not known to be positive", not
      "measured negative". Which policy the model should use is an open
      decision; both are reported so the choice can be made on numbers.
    - No age, KPS or sex, where other patients have one. The script stops
      (--missing-demographic stop). ``drop`` excludes those patients and counts
      them. Nothing is imputed.
    - A field that is absent for every patient is not measured in this source.
      Resection extent is such a field in the TCGA-GBM clinical table. Those
      patients are kept, the features the L3 code derives from the field are
      constant, and both facts are reported (``not_measured_in_source``,
      ``constant_features``).

WHAT THIS SCRIPT DOES NOT DO
    It reads no CGGA file and opens no network connection; a path containing
    "cgga" is refused. It fits no model and reads no outcome. It modifies
    nothing under data/raw.

Usage:
    python scripts/build_tcga_feature_store.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path

import pandas as pd

from gliostrata.features.store import FEATURE_COLUMNS, build_feature_store
from gliostrata.profiling.deconvolve import (
    deconvolve_states_expression,
    deconvolve_tme_expression,
)
from gliostrata.profiling.verhaak import SUBTYPES, classify_verhaak

DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
DEFAULT_EXPRESSION = Path("data/raw/tcga_gbm/expression/tpm.parquet")
DEFAULT_ANNOTATION = Path("data/raw/tcga_gbm/expression/gene_annotation.parquet")
DEFAULT_JOINED = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
DEFAULT_OUT = Path("data/processed/tcga_gbm/feature_store.parquet")
DEFAULT_REPORT = Path("data/results/tcga_gbm/feature_store_report.json")

COHORT_COLUMNS = (
    "patient_id", "age", "kps", "sex", "extent_of_resection",
    "overall_survival_months", "vital_status",
)
MOLECULAR_COLUMNS = ("patient_id", "mgmt_methylated", "idh_mutated")
DEMOGRAPHIC_COLUMNS = ("age", "kps", "sex", "extent_of_resection")
STATE_COLUMNS = ("state_npc", "state_opc", "state_ac", "state_mes")
TME_COLUMNS = ("tme_immune", "tme_stromal", "tme_endothelial")
#: Cohort column names that an L2 output or the joined calls also use.
SHADOWED_COLUMNS = ("subtype", *STATE_COLUMNS, *TME_COLUMNS,
                    "mgmt_methylated", "idh_mutated")
OUTCOME_COLUMNS = ("overall_survival_months", "vital_status", "treatment_arm")

MISSING_MOLECULAR_POLICIES = ("drop", "zero")
MISSING_DEMOGRAPHIC_POLICIES = ("stop", "drop")


def sha256_file(path: Path | str) -> str:
    """SHA-256 of a file's bytes, read in chunks (the TPM matrix is 70 MB)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def refuse_external(path: Path | str) -> None:
    """This script is for the development cohort only."""
    if "cgga" in str(path).lower():
        raise ValueError(
            f"refusing to read {path}: the held-out cohort needs the read guard "
            "(gliostrata.evaluation.read_guard); this script reads TCGA only"
        )


def _require(frame: pd.DataFrame, columns: tuple[str, ...], label: str) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(f"{label}: required column(s) missing: {missing}")


def _one_row_per_patient(ids: pd.Series, label: str) -> None:
    if ids.duplicated().any():
        raise ValueError(f"{label}: duplicate patient_id values")


def drop_placeholders(cohort: pd.DataFrame, label: str) -> pd.DataFrame:
    """Remove all-null cohort columns that would shadow a computed value.

    ``build_feature_store`` keeps the cohort's column when a name collides
    with an L2 output or a joined call, so a cohort column with values under
    such a name would silently replace the computed one. An all-null column is
    a placeholder and is dropped. A column with values stops the script.
    """
    overlapping = [c for c in cohort.columns if c in SHADOWED_COLUMNS]
    with_values = [c for c in overlapping if cohort[c].notna().any()]
    if with_values:
        raise ValueError(
            f"{label}: column(s) {with_values} carry values and would shadow the "
            "computed calls; remove them from the cohort"
        )
    return cohort.drop(columns=overlapping)


def _call_state(calls: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Split a molecular-call column into (positive, missing).

    A call of 1 is positive. A call of 0, or any other value, is not positive.
    A null (no array, no variant, or no row after the left join) is missing.
    """
    missing = calls.isna()
    positive = calls.eq(1).fillna(False).astype(bool)
    return positive, missing


def _check_calls(frame: pd.DataFrame, label: str, expected: list[str]) -> None:
    _require(frame, ("patient_id",), f"{label} calls")
    got = frame["patient_id"].astype(str)
    _one_row_per_patient(got, f"{label} calls")
    if set(got) != set(expected):
        raise ValueError(f"{label} calls do not match the scoring set")
    if frame.drop(columns=["patient_id"]).isna().to_numpy().any():
        raise ValueError(
            f"{label} calls contain a null value; check that the expression "
            "index and the annotation agree"
        )


def build(
    cohort: pd.DataFrame,
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None,
    joined: pd.DataFrame,
    *,
    classify=None,
    state_scorer=None,
    tme_scorer=None,
    missing_molecular: str = "drop",
    missing_demographic: str = "stop",
) -> tuple[pd.DataFrame, dict]:
    """Assemble the feature store from real inputs and summarise what happened.

    The three scorers are injectable for testing; production code leaves them
    as ``None``, which means the committed functions.
    """
    if missing_molecular not in MISSING_MOLECULAR_POLICIES:
        raise ValueError(f"missing_molecular must be one of {MISSING_MOLECULAR_POLICIES}")
    if missing_demographic not in MISSING_DEMOGRAPHIC_POLICIES:
        raise ValueError(f"missing_demographic must be one of {MISSING_DEMOGRAPHIC_POLICIES}")
    _require(cohort, COHORT_COLUMNS, "cohort")
    _require(joined, MOLECULAR_COLUMNS, "joined")

    cohort = cohort.copy()
    cohort["patient_id"] = cohort["patient_id"].astype(str)
    _one_row_per_patient(cohort["patient_id"], "cohort")
    cohort = drop_placeholders(cohort, "cohort")
    joined = joined.loc[:, list(MOLECULAR_COLUMNS)].copy()
    joined["patient_id"] = joined["patient_id"].astype(str)
    _one_row_per_patient(joined["patient_id"], "joined")

    samples = [str(c) for c in expression.columns]
    if len(set(samples)) != len(samples):
        raise ValueError("expression: duplicate sample columns")
    n_cohort = len(cohort)

    # The scoring set: cohort patients that have an expression sample, in
    # cohort order. Expression samples outside the cohort are not scored.
    has_expression = cohort["patient_id"].isin(set(samples))
    scored_ids = cohort.loc[has_expression, "patient_id"].tolist()
    if not scored_ids:
        raise ValueError("no cohort patient has an expression sample")
    scored_expression = expression[scored_ids]

    classify = classify_verhaak if classify is None else classify
    state_scorer = deconvolve_states_expression if state_scorer is None else state_scorer
    tme_scorer = deconvolve_tme_expression if tme_scorer is None else tme_scorer
    subtypes = classify(scored_expression, annotation)
    states = state_scorer(scored_expression, annotation)
    tme = tme_scorer(scored_expression, annotation)
    _check_calls(subtypes, "subtype", scored_ids)
    _check_calls(states, "state", scored_ids)
    _check_calls(tme, "TME", scored_ids)

    # Molecular calls onto the cohort. A patient without a joined row has no call.
    patients = cohort.loc[has_expression].merge(joined, on="patient_id", how="left")
    n_without_joined_row = int((~patients["patient_id"].isin(set(joined["patient_id"]))).sum())
    mgmt_positive, mgmt_missing = _call_state(patients["mgmt_methylated"])
    idh_positive, idh_missing = _call_state(patients["idh_mutated"])
    molecular_missing = mgmt_missing | idh_missing

    n_missing_demographic = {c: int(patients[c].isna().sum()) for c in DEMOGRAPHIC_COLUMNS}
    # A field absent for every patient is not measured in this source. Those
    # patients are kept: excluding them would exclude the whole cohort and say
    # nothing. The policy applies only where some patients have a value.
    not_measured = [c for c in DEMOGRAPHIC_COLUMNS
                    if n_missing_demographic[c] == len(patients)]
    partially_missing = [c for c in DEMOGRAPHIC_COLUMNS
                         if c not in not_measured and n_missing_demographic[c]]
    demographic_missing = pd.Series(False, index=patients.index)
    for col in partially_missing:
        demographic_missing |= patients[col].isna()
    if missing_demographic == "stop" and demographic_missing.any():
        counts = {c: n_missing_demographic[c] for c in partially_missing}
        raise ValueError(
            f"{int(demographic_missing.sum())} patient(s) have no "
            f"{' or '.join(partially_missing)} {counts}. This script does not "
            "impute; pass --missing-demographic drop to exclude them"
        )

    keep = ~demographic_missing
    if missing_molecular == "drop":
        keep &= ~molecular_missing
    else:
        patients["mgmt_methylated"] = mgmt_positive.astype(int)
        patients["idh_mutated"] = idh_positive.astype(int)
    n_excluded_demographic = int(demographic_missing.sum())
    n_excluded_molecular = 0
    if missing_molecular == "drop":
        n_excluded_molecular = int((molecular_missing & ~demographic_missing).sum())
    patients = patients.loc[keep].copy()
    if patients.empty:
        raise ValueError("no patient has a complete feature set")
    resection_counts = {
        str(k): int(v)
        for k, v in patients["extent_of_resection"].value_counts(dropna=False).items()
    }

    store = build_feature_store(patients, subtypes, states, tme)
    if len(store) != len(patients):
        raise ValueError(f"feature store has {len(store)} rows; expected {len(patients)}")
    store_ids = store["patient_id"].astype(str)
    _one_row_per_patient(store_ids, "feature store")
    if set(store_ids) != set(patients["patient_id"]):
        raise ValueError("feature store patients differ from the eligible patients")
    nulls = store[list(FEATURE_COLUMNS)].isna().sum()
    offending = {c: int(n) for c, n in nulls.items() if n}
    if offending:
        raise ValueError(f"feature store has null values: {offending}")
    try:
        store[list(FEATURE_COLUMNS)].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("feature store has a non-numeric feature column") from exc
    unknown = set(store["subtype"]) - set(SUBTYPES)
    if unknown:
        raise ValueError(f"feature store has unknown subtype labels {sorted(unknown)}")
    constant_features = [c for c in FEATURE_COLUMNS if store[c].nunique(dropna=False) <= 1]
    n_missing_outcome = {c: int(store[c].isna().sum()) for c in OUTCOME_COLUMNS
                         if c in store.columns}

    prob_cols = [f"subtype_prob_{s}" for s in SUBTYPES]
    state_sums = store[list(STATE_COLUMNS)].sum(axis=1)
    tme_sums = store[list(TME_COLUMNS)].sum(axis=1)
    summary = {
        "n_cohort": n_cohort,
        "n_expression_samples": len(samples),
        "n_expression_samples_outside_cohort": len(samples) - len(scored_ids),
        "n_scoring_set": len(scored_ids),
        "n_excluded_without_expression": n_cohort - len(scored_ids),
        "n_excluded_missing_demographic": n_excluded_demographic,
        "n_excluded_missing_molecular": n_excluded_molecular,
        "not_measured_in_source": not_measured,
        "constant_features": constant_features,
        "n_missing_outcome": n_missing_outcome,
        "n_missing": {
            "without_joined_row": n_without_joined_row,
            "mgmt_call": int(mgmt_missing.sum()),
            "idh_call": int(idh_missing.sum()),
            "molecular_any": int(molecular_missing.sum()),
            **n_missing_demographic,
        },
        "n_patients": len(store),
        "n_positive": {
            "mgmt_methylated": int(store["mgmt_methylated"].sum()),
            "idh_mutated": int(store["idh_mutated"].sum()),
        },
        "subtype_counts": {st: int((store["subtype"] == st).sum()) for st in SUBTYPES},
        "median_max_subtype_prob": float(store[prob_cols].max(axis=1).median()),
        "row_sum_range": {
            "state": [float(state_sums.min()), float(state_sums.max())],
            "tme": [float(tme_sums.min()), float(tme_sums.max())],
        },
        "extent_of_resection_counts": resection_counts,
        "feature_ranges": {
            col: {"min": float(store[col].min()),
                  "median": float(store[col].median()),
                  "max": float(store[col].max())}
            for col in FEATURE_COLUMNS
        },
        "policies": {
            "missing_molecular": missing_molecular,
            "missing_demographic": missing_demographic,
        },
    }
    return store, summary


def _package_version() -> str:
    try:
        return package_version("gliostrata")
    except PackageNotFoundError:
        return "unknown"


def run(
    cohort_path: Path | str,
    expression_path: Path | str,
    annotation_path: Path | str,
    joined_path: Path | str,
    out_path: Path | str,
    report_path: Path | str,
    **kwargs,
) -> dict:
    """Read the four inputs, build the store, write it and the report."""
    inputs = {
        "cohort": cohort_path,
        "expression": expression_path,
        "annotation": annotation_path,
        "joined": joined_path,
    }
    for path in inputs.values():
        refuse_external(path)
    cohort = pd.read_parquet(cohort_path)
    expression = pd.read_parquet(expression_path)
    annotation = pd.read_parquet(annotation_path)
    joined = pd.read_parquet(joined_path)
    store, summary = build(cohort, expression, annotation, joined, **kwargs)

    out_path, report_path = Path(out_path), Path(report_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    store.to_parquet(out_path, index=False)
    report = {
        "inputs": {
            name: {
                "path": str(path).replace("\\", "/"),
                "bytes": Path(path).stat().st_size,
                "sha256": sha256_file(path),
            }
            for name, path in inputs.items()
        },
        "software": {
            "gliostrata": _package_version(),
            "python": platform.python_version(),
            "pandas": pd.__version__,
        },
        "scoring_set": "cohort patients with an expression sample",
        "feature_columns": list(FEATURE_COLUMNS),
        "not_validated": [
            "Subtype labels are a within-cohort relative ranking and depend on "
            "the scoring set.",
            "State and TME fractions are normalised marker-enrichment scores, "
            "not cell proportions.",
        ],
        **summary,
        "output": {
            "path": str(out_path).replace("\\", "/"),
            "bytes": out_path.stat().st_size,
            "sha256": sha256_file(out_path),
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="ascii"
    )
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT),
                        help="clinical_cohort.parquet")
    parser.add_argument("--expression", default=str(DEFAULT_EXPRESSION),
                        help="gene x sample TPM matrix")
    parser.add_argument("--annotation", default=str(DEFAULT_ANNOTATION),
                        help="gene_annotation.parquet")
    parser.add_argument("--joined", default=str(DEFAULT_JOINED),
                        help="molecular_calls_joined.parquet")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="feature store (parquet)")
    parser.add_argument("--report", default=str(DEFAULT_REPORT), help="report (JSON)")
    parser.add_argument("--missing-molecular", choices=MISSING_MOLECULAR_POLICIES,
                        default="drop",
                        help="patients with no MGMT or IDH call: exclude and count "
                             "(drop), or record the call as 0 (zero)")
    parser.add_argument("--missing-demographic", choices=MISSING_DEMOGRAPHIC_POLICIES,
                        default="stop",
                        help="patients with no age, KPS, sex or resection: stop, or "
                             "exclude and count")
    args = parser.parse_args(argv)

    print("Real-data vertical slice, L0 to L3 (development cohort, TCGA)")
    print("  Committed L2 and L3 modules, unchanged. Nothing here is validated:")
    print("  subtype labels are a within-cohort ranking; state and TME values")
    print("  are normalised enrichment scores, not cell proportions.")
    print()
    r = run(args.cohort, args.expression, args.annotation, args.joined, args.out,
            args.report, missing_molecular=args.missing_molecular,
            missing_demographic=args.missing_demographic)

    print(f"cohort                              {r['n_cohort']}")
    print(f"  no expression sample, excluded    {r['n_excluded_without_expression']}")
    print(f"  no demographic value, excluded    {r['n_excluded_missing_demographic']}")
    print(f"  no molecular call, excluded       {r['n_excluded_missing_molecular']}")
    print(f"patients in the feature store       {r['n_patients']}")
    print()
    print(f"expression samples                  {r['n_expression_samples']}"
          f"  ({r['n_expression_samples_outside_cohort']} outside the cohort, not scored)")
    print(f"  no joined row                      {r['n_missing']['without_joined_row']}")
    print(f"  no MGMT call                       {r['n_missing']['mgmt_call']}")
    print(f"  no IDH call                        {r['n_missing']['idh_call']}")
    print(f"  no KPS                             {r['n_missing']['kps']}")
    if r["not_measured_in_source"]:
        print(f"not measured in this source:        "
              f"{', '.join(r['not_measured_in_source'])}")
    if r["constant_features"]:
        print(f"{'constant features (' + str(len(r['constant_features'])) + ')':<36}"
              f"{', '.join(r['constant_features'])}")
    print(f"store patients with no vital status  "
          f"{r['n_missing_outcome'].get('vital_status', 0)}")
    print()
    print("subtype (within-cohort label)       n")
    for st, n in r["subtype_counts"].items():
        print(f"  {st:<32}{n}")
    print()
    print(f"policy: missing molecular calls {args.missing_molecular}; "
          f"missing demographic values {args.missing_demographic}")
    print(f"\nWrote: {args.out}\n      {args.report}")


if __name__ == "__main__":
    main()
