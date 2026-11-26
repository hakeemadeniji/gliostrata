"""
Pre-registered IDH sensitivity analysis, part 2: the subtype distribution.

WHY THIS EXISTS
    docs/mgmt_idh_calling.md section 6.3 pre-specified this analysis before
    the IDH calls existed. All IDH-mutant patients stay in the primary
    cohort. Under WHO 2021 those cases would be classified as astrocytoma,
    IDH-mutant, so the subtype distribution is re-reported here with them
    excluded. This is a robustness check on the primary analysis, not a
    change to it.

WHAT IS COMPARED
    The committed classifier, unchanged:
    gliostrata.profiling.verhaak.classify_verhaak with its defaults, run
    twice on the gene x sample TPM matrix:
        full          every expression sample
        without IDH   the IDH-mutant expression samples removed
    The reported quantity is the called-subtype distribution (counts and
    percentages) each time, and the percentage-point change per subtype.

    IDH status comes from the joined table written by
    scripts/join_molecular_calls.py. Only idh_mutated == 1 is excluded. A
    sample with no IDH call is NOT excluded.

A PROPERTY OF THE CLASSIFIER THAT THIS ANALYSIS EXPOSES
    classify_verhaak labels a sample by the argmax of scores that are
    standardised per signature ACROSS SAMPLES (see _softmax_over_cohort_z:
    "the output is a within-cohort relative ranking"). The label therefore
    depends on which other samples are in the matrix. Removing 16 samples
    can change the label of a sample that was not removed. The aggregate
    percentage-point change hides that, so this script also reports how
    many retained samples change label, and what the change would be if
    every retained sample kept its full-cohort label. That second number
    is a decomposition aid, not a second classifier: labels are always the
    committed classifier's, never frozen by this script in its own output.

    An expression sample with no row in the joined table (it is outside the
    clinical cohort) has no IDH status there. It is kept and listed, never
    silently dropped.

WHAT THIS SCRIPT DOES NOT DO
    It does not modify the classifier, the expression matrix or anything
    under data/raw, and writes only --out.

Usage:
    python scripts/idh_sensitivity_subtype.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from gliostrata.profiling.verhaak import SUBTYPES, classify_verhaak

DEFAULT_EXPRESSION = Path("data/raw/tcga_gbm/expression/tpm.parquet")
DEFAULT_ANNOTATION = Path("data/raw/tcga_gbm/expression/gene_annotation.parquet")
DEFAULT_JOINED = Path("data/processed/tcga_gbm/molecular_calls_joined.parquet")
DEFAULT_OUT = Path("data/results/tcga_gbm/idh_sensitivity_subtype.json")

#: Result keys, in order. docs/mgmt_idh_calling.md section 6.3 reads these.
RESULT_KEYS = (
    "n_expression_full",
    "n_idh_mutant_excluded",
    "n_expression_without_idh",
    "distribution_full",
    "distribution_without_idh",
    "pp_change_per_subtype",
    "max_pp_change",
)


def distribution(calls: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Called-subtype counts and percentages, in the canonical order.

    ``percent`` is rounded to two decimals; ``pp_change_per_subtype`` is
    computed from the rounded percentages so the published table adds up.
    """
    n = len(calls)
    counts = calls["subtype"].value_counts() if n else pd.Series(dtype=int)
    out: dict[str, dict[str, float]] = {}
    for st in SUBTYPES:
        c = int(counts.get(st, 0))
        out[st] = {"n": c, "percent": round(100.0 * c / n, 2) if n else 0.0}
    return out


def pp_change(full: dict, without: dict) -> dict[str, float]:
    """Signed percentage-point change per subtype (without minus full)."""
    return {st: round(without[st]["percent"] - full[st]["percent"], 2) for st in SUBTYPES}


def _check_calls(calls: pd.DataFrame, expected: list[str], label: str) -> None:
    """The classifier must return exactly one label per sample, no more."""
    for col in ("patient_id", "subtype"):
        if col not in calls.columns:
            raise ValueError(f"{label}: classifier output lacks column {col!r}")
    got = calls["patient_id"].astype(str)
    if got.duplicated().any():
        raise ValueError(f"{label}: classifier returned duplicate patient_id")
    if set(got) != set(expected):
        raise ValueError(f"{label}: classifier output does not match the samples")
    if calls["subtype"].isna().any():
        raise ValueError(f"{label}: classifier returned a null subtype")
    unknown = set(calls["subtype"]) - set(SUBTYPES)
    if unknown:
        raise ValueError(f"{label}: unknown subtype labels {sorted(unknown)}")


def analyze(expression: pd.DataFrame, annotation: pd.DataFrame, joined: pd.DataFrame,
            *, classify=None) -> dict:
    """Classify the full expression set and the set without IDH mutants.

    `classify` is injectable for testing; production code leaves it as
    None, which means the committed `classify_verhaak`.
    """
    classify = classify_verhaak if classify is None else classify
    for col in ("patient_id", "idh_mutated"):
        if col not in joined.columns:
            raise ValueError(f"joined: required column {col!r} is missing")
    samples = [str(c) for c in expression.columns]
    if len(set(samples)) != len(samples):
        raise ValueError("expression: duplicate sample columns")
    joined_ids = joined["patient_id"].astype(str)
    if joined_ids.duplicated().any():
        raise ValueError("joined: duplicate patient_id values")

    # eq(1) is <NA> where there is no call; fillna makes that an explicit
    # "not known to be mutant". A missing call never excludes a sample.
    is_mutant = joined["idh_mutated"].eq(1).fillna(False).to_numpy(dtype=bool)
    mutants = set(joined_ids[is_mutant])
    excluded = [s for s in samples if s in mutants]
    retained = [s for s in samples if s not in mutants]
    if not retained:
        raise ValueError("excluding the IDH-mutant samples leaves nothing")
    # Retained samples with no joined row have no IDH status there. They are
    # kept (not known to be mutant) and listed so the reader can see it.
    no_row = [s for s in retained if s not in set(joined_ids)]

    full_calls = classify(expression, annotation)
    without_calls = classify(expression[retained], annotation)
    _check_calls(full_calls, samples, "full")
    _check_calls(without_calls, retained, "without IDH")

    dist_full = distribution(full_calls)
    dist_without = distribution(without_calls)
    change = pp_change(dist_full, dist_without)

    # Decomposition aid: the retained samples under their full-cohort labels.
    full_labels = full_calls.set_index("patient_id")["subtype"].astype(str)
    retained_full = full_labels.loc[retained]
    frozen = distribution(pd.DataFrame({"subtype": retained_full.to_numpy()}))
    change_frozen = pp_change(dist_full, frozen)
    # What the excluded samples were called in the full run. frozen + this
    # rebuilds the full distribution exactly, subtype by subtype.
    excluded_dist = distribution(
        pd.DataFrame({"subtype": full_labels.reindex(excluded).to_numpy()})
    )
    without_labels = (
        without_calls.set_index("patient_id")["subtype"].astype(str).loc[retained]
    )
    n_relabelled = int((retained_full != without_labels).sum())

    result = {
        "n_expression_full": int(len(samples)),
        "n_idh_mutant_excluded": int(len(excluded)),
        "n_expression_without_idh": int(len(retained)),
        "distribution_full": dist_full,
        "distribution_without_idh": dist_without,
        "pp_change_per_subtype": change,
        "max_pp_change": max(abs(v) for v in change.values()),
        # --- sample-level movement that the aggregate change hides ---
        "n_relabelled_among_retained": n_relabelled,
        "distribution_retained_with_full_cohort_labels": frozen,
        "pp_change_per_subtype_with_full_cohort_labels": change_frozen,
        "distribution_excluded_with_full_cohort_labels": excluded_dist,
        # --- provenance ---
        "n_retained_without_joined_row": int(len(no_row)),
        "retained_without_joined_row": no_row,
        "idh_mutant_samples_excluded": excluded,
        "subtype_order": list(SUBTYPES),
        "classifier": (
            "gliostrata.profiling.verhaak.classify_verhaak (defaults); labels "
            "are a within-cohort relative ranking and depend on cohort "
            "composition"
        ),
        "exclusion_rule": (
            "idh_mutated == 1 only; a sample with no IDH call is not excluded"
        ),
    }
    assert tuple(list(result)[: len(RESULT_KEYS)]) == RESULT_KEYS
    return result


def run(
    expression_path: Path | str,
    annotation_path: Path | str,
    joined_path: Path | str,
    out_path: Path | str,
    **kwargs,
) -> dict:
    """Execute the analysis, write the JSON result, and return it."""
    expression = pd.read_parquet(expression_path)
    annotation = pd.read_parquet(annotation_path)
    joined = pd.read_parquet(joined_path, columns=["patient_id", "idh_mutated"])
    result = analyze(expression, annotation, joined, **kwargs)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--expression", default=str(DEFAULT_EXPRESSION),
                        help="gene x sample TPM matrix")
    parser.add_argument("--annotation", default=str(DEFAULT_ANNOTATION),
                        help="gene_annotation.parquet")
    parser.add_argument("--joined", default=str(DEFAULT_JOINED),
                        help="molecular_calls_joined.parquet")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output JSON")
    args = parser.parse_args(argv)

    print("IDH sensitivity analysis (pre-specified) -- subtype distribution")
    print("  Committed Verhaak classifier, run twice: full expression set,")
    print("  and the set without the IDH-mutant samples. Robustness check")
    print("  only; the primary cohort is not changed.")
    print()

    r = run(args.expression, args.annotation, args.joined, args.out)

    print(f"{'subtype':<12}{'full n':>8}{'full %':>9}{'no-IDH n':>11}"
          f"{'no-IDH %':>10}{'change pp':>12}")
    for st in r["subtype_order"]:
        f = r["distribution_full"][st]
        w = r["distribution_without_idh"][st]
        c = r["pp_change_per_subtype"][st]
        print(f"{st:<12}{f['n']:>8}{f['percent']:>9.2f}{w['n']:>11}"
              f"{w['percent']:>10.2f}{c:>+12.2f}")
    print()
    print(f"expression samples                 {r['n_expression_full']}")
    print(f"  IDH-mutant, excluded             {r['n_idh_mutant_excluded']}")
    print(f"  retained                         {r['n_expression_without_idh']}")
    print(f"largest per-subtype change         {r['max_pp_change']:.2f} pp")
    print(f"retained samples with a new label  {r['n_relabelled_among_retained']}"
          "  (labels are a within-cohort ranking)")
    print(f"\nWrote: {args.out}")


if __name__ == "__main__":
    main()
