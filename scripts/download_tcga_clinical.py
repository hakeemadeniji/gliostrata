"""
Download TCGA-GBM clinical data from the GDC REST API and map to the
GlioStrata cohort schema.

Phase A scope: clinical data only. Omic-derived fields (MGMT, IDH, subtype,
Neftel states, TME fractions) are set to pd.NA -- they require Phases B/C
(methylation, mutations, expression, scRNA-seq).

Usage:
    python scripts/download_tcga_clinical.py --out data/raw/tcga_gbm
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import requests

GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"

FIELDS = [
    "case_id",
    "submitter_id",
    "demographic.vital_status",
    "demographic.age_at_index",
    "demographic.sex_at_birth",
    "demographic.days_to_death",
    "diagnoses.primary_diagnosis",
    "diagnoses.tumor_grade",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.age_at_diagnosis",
    "diagnoses.treatments.treatment_type",
    "diagnoses.treatments.treatment_outcome",
    "follow_ups.days_to_follow_up",
    "follow_ups.karnofsky_performance_status",
    "exposures.tobacco_smoking_status",
]


def fetch_cases(project_id: str, size: int) -> list[dict]:
    """Fetch cases for a project from the GDC API."""
    filters = {
        "op": "in",
        "content": {"field": "project.project_id", "value": [project_id]},
    }
    payload = {
        "filters": filters,
        "fields": ",".join(FIELDS),
        "format": "JSON",
        "size": size,
    }
    print(f"Querying GDC API for {project_id} (size={size})...")
    response = requests.post(GDC_CASES_URL, json=payload, timeout=120)
    if response.status_code != 200:
        raise RuntimeError(
            f"GDC API error {response.status_code}: {response.text[:500]}"
        )
    data = response.json()
    api_warnings = data.get("warnings", {})
    if api_warnings:
        print(f"WARNING - GDC API field warnings: {api_warnings}")
    hits = data["data"]["hits"]
    print(f"Fetched {len(hits)} cases")
    return hits


def flatten_case(case: dict) -> dict:
    """Flatten one GDC case record to a single-row dict."""
    demo = case.get("demographic") or {}
    diagnoses = case.get("diagnoses") or []
    follow_ups = case.get("follow_ups") or []
    exposures = case.get("exposures") or []

    diag = diagnoses[0] if diagnoses else {}
    treatments = diag.get("treatments") or []
    exposure = exposures[0] if exposures else {}

    # KPS: scan all follow-ups for the first non-null value. GDC returns it as
    # a string ("80", sometimes "Unknown"); build_cohort coerces to numeric.
    # The follow_ups list is unordered and the KPS-bearing record often has a
    # null days_to_follow_up, so the first non-null value is used.
    kps = None
    for fu in follow_ups:
        v = fu.get("karnofsky_performance_status")
        if v is not None:
            kps = v
            break

    # days_to_follow_up: take the maximum (most recent follow-up).
    days_to_follow_up = None
    for fu in follow_ups:
        v = fu.get("days_to_follow_up")
        if v is not None:
            days_to_follow_up = max(days_to_follow_up or 0, v)

    treatment_types = [
        t.get("treatment_type") for t in treatments if t.get("treatment_type")
    ]
    treatment_outcomes = [
        t.get("treatment_outcome") for t in treatments if t.get("treatment_outcome")
    ]

    return {
        "case_id": case.get("case_id"),
        "submitter_id": case.get("submitter_id"),
        "vital_status_raw": demo.get("vital_status"),
        "age_at_index_raw": demo.get("age_at_index"),
        "sex_at_birth_raw": demo.get("sex_at_birth"),
        "days_to_death": demo.get("days_to_death"),
        "primary_diagnosis": diag.get("primary_diagnosis"),
        "tumor_grade": diag.get("tumor_grade"),
        "ajcc_stage": diag.get("ajcc_pathologic_stage"),
        "age_at_diagnosis_raw": diag.get("age_at_diagnosis"),
        "treatment_types": json.dumps(treatment_types),
        "treatment_outcomes": json.dumps(treatment_outcomes),
        "days_to_follow_up": days_to_follow_up,
        "kps_raw": kps,
        "smoking_status": exposure.get("tobacco_smoking_status"),
    }


#: Label for cases where treatment was recorded but contained neither
#: chemotherapy nor radiation (e.g. surgery only).
ARM_NO_ADJUVANT = "NoAdjuvant"

#: Deterministic order in which filter reasons are evaluated and reported.
FILTER_REASONS = [
    ("contradictory_diagnosis", "primary_diagnosis names a non-glioblastoma histology"),
    ("no_usable_data", "no age, sex, survival or vital status"),
    ("nonpositive_survival", "overall_survival_months <= 0"),
]


def classify_arm(types_json: str):
    """Classify a JSON-encoded list of GDC treatment_type strings.

    Returns:
        "ChemoRadiation" - both chemotherapy and radiation recorded.
        "MonoTherapy"    - chemotherapy or radiation (not both) recorded.
        "NoAdjuvant"     - treatment recorded, but no chemo and no radiation
                           (e.g. surgery only). A real, informative category.
        pd.NA            - treatment was NOT RECORDED (empty/missing/unparseable).
    """
    try:
        parsed = json.loads(types_json)
    except (json.JSONDecodeError, TypeError, ValueError):
        return pd.NA
    if not isinstance(parsed, list):
        return pd.NA
    types = [t.lower() for t in parsed if isinstance(t, str) and t]
    if not types:
        # No treatment entries at all -> treatment was not recorded.
        return pd.NA
    has_chemo = any("pharmaceutical" in t or "chemotherapy" in t for t in types)
    has_rad = any("radiation" in t for t in types)
    if has_chemo and has_rad:
        return "ChemoRadiation"
    if has_chemo or has_rad:
        return "MonoTherapy"
    return ARM_NO_ADJUVANT


#: primary_diagnosis values that record no histology rather than a different
#: one. Every case in this download is enrolled in the TCGA-GBM project, so an
#: absent diagnosis is a data-entry gap, not evidence of another disease.
UNINFORMATIVE_DIAGNOSES = frozenset({"", "unknown", "not reported", "na", "n/a", "none"})


def is_glioblastoma(value) -> bool:
    """True only for an explicit glioblastoma primary diagnosis.

    Case-insensitive substring match on "glioblastoma", so variants such as
    "Glioblastoma, NOS" are kept, while NaN, "Unknown", "Not Reported" and
    other histologies are rejected. This is the STRICT test; the cohort filter
    uses `is_diagnosis_contradictory` instead -- see that function for why.
    """
    if value is None or value is pd.NA:
        return False
    if isinstance(value, float) and pd.isna(value):
        return False
    if not isinstance(value, str):
        return False
    return "glioblastoma" in value.strip().lower()


def is_diagnosis_contradictory(value) -> bool:
    """True only when primary_diagnosis names a histology that is NOT glioblastoma.

    A missing, blank, "Unknown" or "Not Reported" diagnosis is uninformative,
    not contradictory: these cases were enrolled in TCGA-GBM, so the project
    itself establishes the disease and the empty field is a curation gap. They
    are retained. A field naming a different histology (e.g. "Adenocarcinoma,
    NOS") genuinely conflicts with the cohort definition and is dropped.

    This matters because 8 of the patients with an absent diagnosis have RNA-seq
    data; excluding them would discard usable molecular samples on the strength
    of a blank clinical field.
    """
    if value is None or value is pd.NA:
        return False
    if isinstance(value, float) and pd.isna(value):
        return False
    if not isinstance(value, str):
        return False
    text = value.strip().lower()
    if text in UNINFORMATIVE_DIAGNOSES:
        return False
    return "glioblastoma" not in text


def filter_cohort(cohort: pd.DataFrame, raw: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Drop rows that do not belong in a survival-modelling cohort.

    Reasons are evaluated in the fixed order given by FILTER_REASONS. The
    reason masks overlap, so the report records both the number of rows
    matching each reason independently ("matched") and the number of rows
    attributed to each reason under first-match attribution ("dropped"),
    which sums exactly to "dropped_total".

    Args:
        cohort: mapped cohort (output of build_cohort).
        raw: the matching raw download, used for primary_diagnosis.

    Returns:
        (filtered cohort, report dict). Row order is preserved.
    """
    if len(cohort) != len(raw):
        raise ValueError(
            f"cohort/raw length mismatch: {len(cohort)} vs {len(raw)}"
        )

    diagnosis = raw["primary_diagnosis"].to_numpy()
    contradictory = pd.Series(
        [is_diagnosis_contradictory(v) for v in diagnosis], index=cohort.index
    )

    masks = {
        "contradictory_diagnosis": contradictory,
        "no_usable_data": (
            cohort["age"].isna()
            & cohort["sex"].isna()
            & cohort["overall_survival_months"].isna()
            & cohort["vital_status"].isna()
        ),
        "nonpositive_survival": (
            cohort["overall_survival_months"].notna()
            & (pd.to_numeric(cohort["overall_survival_months"], errors="coerce") <= 0)
        ),
    }

    report = {"input_rows": int(len(cohort))}
    already = pd.Series(False, index=cohort.index)
    for key, _desc in FILTER_REASONS:
        mask = masks[key].fillna(False).astype(bool)
        report[f"matched_{key}"] = int(mask.sum())
        report[f"dropped_{key}"] = int((mask & ~already).sum())
        already = already | mask

    report["dropped_total"] = int(already.sum())
    filtered = cohort.loc[~already].copy()
    report["kept_rows"] = int(len(filtered))
    return filtered, report


def print_filter_report(report: dict) -> None:
    """Print an ASCII-only filter report to stdout."""
    print("\nCohort filter report:")
    print(f"  Input rows: {report['input_rows']}")
    for key, desc in FILTER_REASONS:
        print(
            f"  Dropped - {desc}: {report[f'dropped_{key}']}"
            f"  (rows matching this reason: {report[f'matched_{key}']})"
        )
    print(f"  Total unique rows dropped: {report['dropped_total']}")
    print(f"  Final cohort rows: {report['kept_rows']}")


def detect_age_unit(series: pd.Series) -> str:
    """GDC reports age in days; detect to be safe."""
    vals = pd.to_numeric(series, errors="coerce").dropna()
    if len(vals) == 0:
        return "unknown"
    return "days" if vals.median() > 3650 else "years"


def build_cohort(raw: pd.DataFrame) -> pd.DataFrame:
    """Map raw GDC clinical records to the GlioStrata cohort schema."""
    cohort = pd.DataFrame()
    cohort["patient_id"] = raw["submitter_id"]

    # Age
    age_raw = pd.to_numeric(raw["age_at_index_raw"], errors="coerce")
    unit = detect_age_unit(raw["age_at_index_raw"])
    cohort["age"] = (age_raw / 365.25).round(0) if unit == "days" else age_raw.round(0)

    # Sex: GDC uses "male"/"female"/"unknown"; cohort schema uses "M"/"F".
    cohort["sex"] = raw["sex_at_birth_raw"].map(
        {"male": "M", "female": "F", "unknown": pd.NA, None: pd.NA}
    )

    # KPS
    cohort["kps"] = pd.to_numeric(raw["kps_raw"], errors="coerce")

    # Not reliably available in TCGA clinical data
    cohort["extent_of_resection"] = pd.NA

    # Omic-derived fields: Phase A does not have these
    for col in ["mgmt_methylated", "idh_mutated", "subtype"]:
        cohort[col] = pd.NA
    for col in [
        "state_npc", "state_opc", "state_ac", "state_mes",
        "tme_immune", "tme_stromal", "tme_endothelial",
    ]:
        cohort[col] = pd.NA

    # Survival: deceased -> days_to_death; alive -> days_to_follow_up
    days = pd.to_numeric(raw["days_to_death"], errors="coerce")
    days = days.fillna(pd.to_numeric(raw["days_to_follow_up"], errors="coerce"))
    cohort["overall_survival_months"] = (days / 30.4375).round(2)

    # Vital status: 1 = deceased, 0 = alive, NA = unknown
    cohort["vital_status"] = raw["vital_status_raw"].map({"Dead": 1, "Alive": 0})

    # Treatment arm: nullable string dtype so unrecorded treatment is a real
    # missing value (pd.NA) rather than the literal string "None".
    cohort["treatment_arm"] = (
        raw["treatment_types"].apply(classify_arm).astype("string")
    )

    return cohort


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/raw/tcga_gbm", help="Output directory")
    parser.add_argument("--project", default="TCGA-GBM", help="GDC project ID")
    parser.add_argument("--size", type=int, default=2000, help="Max cases to fetch")
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Skip cohort row filtering (produce the unfiltered mapped cohort)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    hits = fetch_cases(args.project, args.size)
    raw = pd.DataFrame([flatten_case(c) for c in hits])

    raw_path = out_dir / "clinical_raw.parquet"
    raw.to_parquet(raw_path, index=False)
    print(f"\nSaved raw clinical data: {raw_path}")
    print(f"  Shape: {raw.shape[0]} rows x {raw.shape[1]} cols")

    print("\nField coverage (non-null count):")
    for col in raw.columns:
        n = raw[col].notna().sum()
        print(f"  {col:30s} {n:5d} / {len(raw)}")

    cohort = build_cohort(raw)
    if args.no_filter:
        print("\nCohort filter report: SKIPPED (--no-filter)")
        print(f"  Final cohort rows: {len(cohort)}")
    else:
        cohort, report = filter_cohort(cohort, raw)
        print_filter_report(report)

    cohort_path = out_dir / "clinical_cohort.parquet"
    cohort.to_parquet(cohort_path, index=False)
    print(f"\nSaved mapped cohort: {cohort_path}")
    print(f"  Shape: {cohort.shape[0]} rows x {cohort.shape[1]} cols")

    print("\nCohort summary:")
    print(
        cohort[
            ["age", "sex", "vital_status", "overall_survival_months", "treatment_arm"]
        ].describe(include="all")
    )


if __name__ == "__main__":
    main()