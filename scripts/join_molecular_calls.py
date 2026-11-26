"""
Join the IDH and MGMT call tables onto the TCGA-GBM clinical cohort (L2).

WHY THIS EXISTS
    scripts/download_tcga_mutations.py and scripts/download_tcga_methylation.py
    each produce a standalone per-patient call table and deliberately leave
    clinical_cohort.parquet untouched: the ingestion scripts fit nothing and
    join nothing. This script is that separate, reviewed join. It creates
    nothing biological -- every value it writes is copied from a call table.

THE ONE THING THAT MUST NOT GO WRONG
    Both call tables carry a row for EVERY cohort patient, including the
    patients who have no data:
        idh_calls.parquet    idh_mutated = <NA>, call_basis = no_mutation_data
        mgmt_calls.parquet   mgmt_methylated = <NA>, call_basis = no_methylation_data
    So "the patient is in the table" does NOT mean "the patient has a call".
    A patient has a call only when the call value is non-null. This script
    derives has_idh_call / has_mgmt_call from the value, never from row
    presence, and never imputes:
        has_idh_call = False  ->  idh_mutated = <NA>
        has_mgmt_call = False  ->  mgmt_methylated = <NA>, mgmt_score = <NA>
    A missing call is NOT wildtype. Coercing it to 0 invents tested patients.

OUTPUTS (under --out, default data/processed/tcga_gbm)
    molecular_calls_joined.parquet   one row per clinical-cohort patient, in
                                     cohort order: patient_id, idh_mutated
                                     (nullable Int8), mgmt_methylated
                                     (nullable Int8), mgmt_score (float),
                                     has_idh_call (bool), has_mgmt_call (bool)
    molecular_calls_manifest.json    the coverage counts that every downstream
                                     analysis depends on

    Both are derived artifacts (git-ignored). This script regenerates them
    from the IDH and MGMT call tables under data/raw/, which are themselves
    git-ignored outputs of the two download scripts -- so a clean checkout
    reproduces these files only after those scripts have been run.

    clinical_cohort.parquet and everything under data/raw are NEVER modified.
    The exclusion of IDH-mutant patients is NOT done here; it is applied at
    analysis time by scripts/idh_sensitivity_*.py.

Usage:
    python scripts/join_molecular_calls.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")
DEFAULT_IDH_CALLS = Path("data/raw/tcga_gbm/mutations/idh_calls.parquet")
DEFAULT_MGMT_CALLS = Path("data/raw/tcga_gbm/methylation/mgmt_calls.parquet")
DEFAULT_EXPRESSION = Path("data/raw/tcga_gbm/expression/tpm.parquet")
DEFAULT_OUT = Path("data/processed/tcga_gbm")

JOINED_FILENAME = "molecular_calls_joined.parquet"
MANIFEST_FILENAME = "molecular_calls_manifest.json"

#: Columns of molecular_calls_joined.parquet, in order.
JOINED_COLUMNS = (
    "patient_id",
    "idh_mutated",
    "mgmt_methylated",
    "mgmt_score",
    "has_idh_call",
    "has_mgmt_call",
)

#: Manifest keys, in order. Downstream tests and scripts read these by name.
MANIFEST_KEYS = (
    "n_clinical_cohort",
    "n_with_idh_call",
    "n_with_mgmt_call",
    "n_with_both",
    "n_with_neither",
    "n_expression",
    "n_expression_with_idh",
    "n_expression_with_mgmt",
    "n_expression_with_both",
    "n_expression_with_neither",
    "n_idh_mutant",
    "n_idh_mutant_with_expression",
    "n_joined_rows",
)


# --------------------------------------------------------------------------
# Loading and validation
# --------------------------------------------------------------------------
def output_paths(out_dir: Path | str) -> tuple[Path, Path]:
    """Derive the two output paths from the output directory."""
    base = Path(out_dir)
    return base / JOINED_FILENAME, base / MANIFEST_FILENAME


def _check_call_table(calls: pd.DataFrame, value_col: str, label: str) -> None:
    for col in ("patient_id", value_col):
        if col not in calls.columns:
            raise ValueError(f"{label}: required column {col!r} is missing")


def _call_series(calls: pd.DataFrame, value_col: str, label: str) -> pd.Series:
    """Boolean call values as a Series indexed by patient_id.

    Values are 1.0 / 0.0 where a call exists and <NA> where it does not.
    """
    _check_call_table(calls, value_col, label)
    ids = calls["patient_id"].astype(str)
    if ids.duplicated().any():
        dup = sorted(ids[ids.duplicated().unique()])
        raise ValueError(f"{label}: duplicate patient_id values: {dup}")
    values = pd.to_numeric(calls[value_col], errors="raise").astype("Float64")
    bad = values.dropna()
    bad = bad[~bad.isin([0.0, 1.0])]
    if len(bad):
        found = sorted({float(v) for v in bad})
        raise ValueError(
            f"{label}: {value_col} must be 0, 1 or <NA>; found {found}"
        )
    values.index = pd.Index(ids, name="patient_id")
    return values


def _score_series(calls: pd.DataFrame, label: str) -> pd.Series:
    """mgmt_score as a float Series indexed by patient_id."""
    _check_call_table(calls, "mgmt_score", label)
    ids = calls["patient_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError(f"{label}: duplicate patient_id values")
    values = pd.to_numeric(calls["mgmt_score"], errors="raise").astype("Float64")
    values.index = pd.Index(ids, name="patient_id")
    return values


def expression_patient_ids(path: Path | str) -> list[str]:
    """Patient ids that are the sample columns of the gene x sample TPM matrix.

    Read from the parquet schema, so the 60k-gene matrix is never loaded.
    """
    import pyarrow.parquet as pq

    names = pq.ParquetFile(path).schema.names
    return [str(c) for c in names if str(c).startswith("TCGA-")]


# --------------------------------------------------------------------------
# The join
# --------------------------------------------------------------------------
def build_joined(
    cohort: pd.DataFrame, idh_calls: pd.DataFrame, mgmt_calls: pd.DataFrame
) -> pd.DataFrame:
    """Join both call tables onto the clinical cohort.

    One row per cohort patient, in cohort order. Patients who appear in a
    call table but not in the cohort are not emitted. Nothing is imputed:
    where a patient has no call, the call value is <NA> and the has_*_call
    flag is False.
    """
    if "patient_id" not in cohort.columns:
        raise ValueError("cohort: required column 'patient_id' is missing")
    ids = cohort["patient_id"].astype(str)
    if ids.duplicated().any():
        raise ValueError("cohort: duplicate patient_id values")

    idh = _call_series(idh_calls, "idh_mutated", "idh_calls")
    meth = _call_series(mgmt_calls, "mgmt_methylated", "mgmt_calls")
    score = _score_series(mgmt_calls, "mgmt_calls")

    # A score without a call would be a value nobody can interpret; refuse.
    score_without_call = score.notna() & ~meth.reindex(score.index).notna()
    if score_without_call.any():
        bad = sorted(score.index[score_without_call])[:5]
        raise ValueError(
            f"mgmt_calls: mgmt_score is present without a methylation call "
            f"for e.g. {bad}"
        )

    joined = pd.DataFrame({"patient_id": ids.reset_index(drop=True)})
    joined = joined.merge(
        idh.rename("idh_mutated").reset_index(), on="patient_id", how="left"
    )
    joined = joined.merge(
        meth.rename("mgmt_methylated").reset_index(), on="patient_id", how="left"
    )
    joined = joined.merge(
        score.rename("mgmt_score").reset_index(), on="patient_id", how="left"
    )

    joined["idh_mutated"] = pd.array(joined["idh_mutated"], dtype="Int8")
    joined["mgmt_methylated"] = pd.array(joined["mgmt_methylated"], dtype="Int8")
    joined["has_idh_call"] = joined["idh_mutated"].notna().astype(bool)
    joined["has_mgmt_call"] = joined["mgmt_methylated"].notna().astype(bool)
    joined["mgmt_score"] = pd.to_numeric(joined["mgmt_score"]).astype("float64")
    joined["mgmt_score"] = joined["mgmt_score"].where(joined["has_mgmt_call"])

    # Guards on the property this script exists to protect.
    if not (joined["has_idh_call"] == joined["idh_mutated"].notna()).all():
        raise AssertionError("has_idh_call disagrees with idh_mutated")
    if not (joined["has_mgmt_call"] == joined["mgmt_methylated"].notna()).all():
        raise AssertionError("has_mgmt_call disagrees with mgmt_methylated")
    if joined.loc[~joined["has_mgmt_call"], "mgmt_score"].notna().any():
        raise AssertionError("mgmt_score must be <NA> where there is no call")

    return joined.loc[:, list(JOINED_COLUMNS)].reset_index(drop=True)


def build_manifest(
    joined: pd.DataFrame,
    expression_ids,
    idh_calls: pd.DataFrame,
    mgmt_calls: pd.DataFrame,
) -> dict:
    """Coverage counts, under the exact keys downstream code depends on.

    Cohort counts come from the joined table. Expression counts come from
    the call tables directly, because the expression set is not a subset of
    the clinical cohort (TCGA-28-2510 has expression but no cohort row), so
    counting through the 594-row joined table would undercount them.
    """
    expr = pd.Index([str(x) for x in expression_ids], name="patient_id")
    if expr.duplicated().any():
        raise ValueError("expression: duplicate patient columns")

    idh = _call_series(idh_calls, "idh_mutated", "idh_calls")
    meth = _call_series(mgmt_calls, "mgmt_methylated", "mgmt_calls")
    idh_called = set(idh.dropna().index)
    meth_called = set(meth.dropna().index)
    idh_mutant = set(idh[idh == 1.0].index)

    e_idh = expr.isin(idh_called)
    e_mgmt = expr.isin(meth_called)

    n_idh = int(joined["has_idh_call"].sum())
    n_mgmt = int(joined["has_mgmt_call"].sum())
    n_both = int((joined["has_idh_call"] & joined["has_mgmt_call"]).sum())
    n_e_idh = int(e_idh.sum())
    n_e_mgmt = int(e_mgmt.sum())
    n_e_both = int((e_idh & e_mgmt).sum())

    return {
        "n_clinical_cohort": int(len(joined)),
        "n_with_idh_call": n_idh,
        "n_with_mgmt_call": n_mgmt,
        "n_with_both": n_both,
        "n_with_neither": int(len(joined)) - n_idh - n_mgmt + n_both,
        "n_expression": int(len(expr)),
        "n_expression_with_idh": n_e_idh,
        "n_expression_with_mgmt": n_e_mgmt,
        "n_expression_with_both": n_e_both,
        "n_expression_with_neither": int(len(expr)) - n_e_idh - n_e_mgmt + n_e_both,
        "n_idh_mutant": int((joined["idh_mutated"] == 1).sum()),
        "n_idh_mutant_with_expression": int(expr.isin(idh_mutant).sum()),
        "n_joined_rows": int(len(joined)),
    }


def run(
    cohort_path: Path | str,
    idh_path: Path | str,
    mgmt_path: Path | str,
    expression_path: Path | str,
    out_dir: Path | str,
) -> dict:
    """Execute the join, write both outputs, and return the manifest."""
    cohort = pd.read_parquet(cohort_path, columns=["patient_id"])
    idh_calls = pd.read_parquet(idh_path, columns=["patient_id", "idh_mutated"])
    mgmt_calls = pd.read_parquet(
        mgmt_path, columns=["patient_id", "mgmt_methylated", "mgmt_score"]
    )
    expression_ids = expression_patient_ids(expression_path)

    joined = build_joined(cohort, idh_calls, mgmt_calls)
    manifest = build_manifest(joined, expression_ids, idh_calls, mgmt_calls)
    assert tuple(manifest) == MANIFEST_KEYS

    joined_path, manifest_path = output_paths(out_dir)
    joined_path.parent.mkdir(parents=True, exist_ok=True)
    joined.to_parquet(joined_path, index=False)
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="ascii"
    )
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT),
                        help="clinical_cohort.parquet")
    parser.add_argument("--idh-calls", default=str(DEFAULT_IDH_CALLS),
                        help="idh_calls.parquet")
    parser.add_argument("--mgmt-calls", default=str(DEFAULT_MGMT_CALLS),
                        help="mgmt_calls.parquet")
    parser.add_argument("--expression", default=str(DEFAULT_EXPRESSION),
                        help="gene x sample TPM matrix; its columns are the "
                             "expression-covered patients (coverage report only)")
    parser.add_argument("--out", default=str(DEFAULT_OUT),
                        help="Output directory")
    args = parser.parse_args(argv)

    print("Molecular-calls join -- IDH and MGMT calls onto the clinical cohort")
    print("  Joins only. Fits nothing, imputes nothing, excludes nothing;")
    print("  data/raw is not modified.")
    print()

    manifest = run(
        args.cohort, args.idh_calls, args.mgmt_calls, args.expression, args.out
    )

    print(f"clinical cohort           {manifest['n_clinical_cohort']}")
    print(f"  with an IDH call        {manifest['n_with_idh_call']}")
    print(f"  with an MGMT call       {manifest['n_with_mgmt_call']}")
    print(f"  with both               {manifest['n_with_both']}")
    print(f"  with neither            {manifest['n_with_neither']}")
    print(f"  IDH-mutant              {manifest['n_idh_mutant']}")
    print(f"expression patients       {manifest['n_expression']}")
    print(f"  with an IDH call        {manifest['n_expression_with_idh']}")
    print(f"  with an MGMT call       {manifest['n_expression_with_mgmt']}")
    print(f"  with both               {manifest['n_expression_with_both']}")
    print(f"  with neither            {manifest['n_expression_with_neither']}")
    print(f"  IDH-mutant              {manifest['n_idh_mutant_with_expression']}")

    joined_path, manifest_path = output_paths(args.out)
    print(f"\nWrote: {joined_path}")
    print(f"      {manifest_path}")


if __name__ == "__main__":
    main()
