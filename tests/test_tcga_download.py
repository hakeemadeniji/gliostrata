"""Unit tests for scripts/download_tcga_clinical.py (no network access)."""

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "download_tcga_clinical.py"
)


def _load_module():
    """Load the download script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_tcga_download_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


tcga = _load_module()


# --------------------------------------------------------------------------
# classify_arm
# --------------------------------------------------------------------------

def test_classify_arm_empty_list_is_missing():
    assert pd.isna(tcga.classify_arm(json.dumps([])))


@pytest.mark.parametrize(
    "bad",
    ["", "not json", None, "{}", json.dumps({"a": 1}), json.dumps("Surgery")],
)
def test_classify_arm_malformed_is_missing(bad):
    assert pd.isna(tcga.classify_arm(bad))


def test_classify_arm_surgery_only_is_no_adjuvant():
    assert tcga.classify_arm(json.dumps(["Surgery, NOS"])) == "NoAdjuvant"


def test_classify_arm_chemo_and_radiation():
    value = json.dumps(["Pharmaceutical Therapy, NOS", "Radiation Therapy, NOS"])
    assert tcga.classify_arm(value) == "ChemoRadiation"


def test_classify_arm_chemo_only_is_monotherapy():
    assert tcga.classify_arm(json.dumps(["Chemotherapy"])) == "MonoTherapy"


def test_classify_arm_radiation_only_is_monotherapy():
    value = json.dumps(["Radiation Therapy, NOS", "Surgery, NOS"])
    assert tcga.classify_arm(value) == "MonoTherapy"


def test_classify_arm_never_returns_none_string():
    for value in [json.dumps([]), "garbage", json.dumps(["Surgery, NOS"])]:
        result = tcga.classify_arm(value)
        assert pd.isna(result) or result != "None"


# --------------------------------------------------------------------------
# build_cohort / treatment_arm column
# --------------------------------------------------------------------------

def _raw_row(**overrides):
    row = {
        "case_id": "c",
        "submitter_id": "TCGA-00-0000",
        "vital_status_raw": "Dead",
        "age_at_index_raw": 60,
        "sex_at_birth_raw": "male",
        "days_to_death": 365,
        "primary_diagnosis": "Glioblastoma",
        "tumor_grade": "not reported",
        "ajcc_stage": None,
        "age_at_diagnosis_raw": 21915,
        "treatment_types": json.dumps(
            ["Pharmaceutical Therapy, NOS", "Radiation Therapy, NOS"]
        ),
        "treatment_outcomes": json.dumps([]),
        "days_to_follow_up": None,
        "kps_raw": None,
        "smoking_status": None,
    }
    row.update(overrides)
    return row


def _raw_frame():
    return pd.DataFrame(
        [
            _raw_row(submitter_id="P-chemorad"),
            _raw_row(submitter_id="P-unrecorded", treatment_types=json.dumps([])),
            _raw_row(
                submitter_id="P-surgery",
                treatment_types=json.dumps(["Surgery, NOS"]),
            ),
            _raw_row(submitter_id="P-malformed", treatment_types="not json"),
            _raw_row(
                submitter_id="P-chemo",
                treatment_types=json.dumps(["Chemotherapy"]),
            ),
        ]
    )


def test_treatment_arm_column_has_no_literal_none_string():
    cohort = tcga.build_cohort(_raw_frame())
    assert (cohort["treatment_arm"] == "None").sum() == 0
    assert "None" not in set(cohort["treatment_arm"].dropna())


def test_treatment_arm_isna_counts_unrecorded_rows():
    cohort = tcga.build_cohort(_raw_frame())
    assert cohort["treatment_arm"].isna().sum() == 2  # empty list + malformed
    counts = cohort["treatment_arm"].value_counts(dropna=False).to_dict()
    assert counts["ChemoRadiation"] == 1
    assert counts["NoAdjuvant"] == 1
    assert counts["MonoTherapy"] == 1


# --------------------------------------------------------------------------
# filter_cohort
# --------------------------------------------------------------------------

def _filter_fixture():
    raw = pd.DataFrame(
        [
            _raw_row(submitter_id="KEEP-1"),
            _raw_row(submitter_id="KEEP-2", primary_diagnosis="Glioblastoma, NOS"),
            _raw_row(
                submitter_id="DROP-diagnosis",
                primary_diagnosis="Adenocarcinoma, NOS",
            ),
            _raw_row(submitter_id="KEEP-unknown", primary_diagnosis="Unknown"),
            _raw_row(
                submitter_id="KEEP-notreported", primary_diagnosis="Not Reported"
            ),
            _raw_row(
                submitter_id="DROP-empty",
                primary_diagnosis=None,
                vital_status_raw="Not Reported",
                age_at_index_raw=None,
                sex_at_birth_raw=None,
                days_to_death=None,
            ),
            _raw_row(submitter_id="DROP-zero-surv", days_to_death=0),
        ]
    )
    cohort = tcga.build_cohort(raw)
    return cohort, raw


def test_filter_cohort_keeps_only_valid_rows():
    cohort, raw = _filter_fixture()
    filtered, report = tcga.filter_cohort(cohort, raw)
    # "Unknown"/"Not Reported" diagnoses are retained: TCGA-GBM enrollment
    # establishes the disease, so a blank field is a curation gap, not a
    # different histology. Only the named non-GBM histology is dropped.
    assert list(filtered["patient_id"]) == [
        "KEEP-1",
        "KEEP-2",
        "KEEP-unknown",
        "KEEP-notreported",
    ]
    assert report["kept_rows"] == 4


def test_filter_cohort_report_numbers_add_up():
    cohort, raw = _filter_fixture()
    filtered, report = tcga.filter_cohort(cohort, raw)
    assert report["input_rows"] == len(cohort) == 7
    attributed = sum(report[f"dropped_{k}"] for k, _ in tcga.FILTER_REASONS)
    assert attributed == report["dropped_total"] == 3
    assert report["kept_rows"] + report["dropped_total"] == report["input_rows"]
    assert report["kept_rows"] == len(filtered)
    # Only an explicitly different histology is contradictory; NaN/"Unknown"
    # rows are retained because TCGA-GBM enrollment establishes the disease.
    assert report["matched_contradictory_diagnosis"] == 1
    assert report["dropped_contradictory_diagnosis"] == 1
    assert report["matched_no_usable_data"] == 1
    assert report["dropped_no_usable_data"] == 1
    assert report["matched_nonpositive_survival"] == 1
    assert report["dropped_nonpositive_survival"] == 1


def test_filter_cohort_drops_zero_survival_only_when_nonpositive():
    cohort, raw = _filter_fixture()
    filtered, _ = tcga.filter_cohort(cohort, raw)
    assert (filtered["overall_survival_months"].dropna() > 0).all()


def test_filter_cohort_isolated_no_usable_data_row():
    # Glioblastoma diagnosis but nothing else recorded -> dropped as unusable.
    raw = pd.DataFrame(
        [
            _raw_row(submitter_id="KEEP"),
            _raw_row(
                submitter_id="DROP-empty",
                vital_status_raw="Not Reported",
                age_at_index_raw=None,
                sex_at_birth_raw=None,
                days_to_death=None,
            ),
        ]
    )
    cohort = tcga.build_cohort(raw)
    filtered, report = tcga.filter_cohort(cohort, raw)
    assert list(filtered["patient_id"]) == ["KEEP"]
    assert report["dropped_no_usable_data"] == 1
    assert report["dropped_contradictory_diagnosis"] == 0


def test_filter_cohort_preserves_columns_and_order():
    cohort, raw = _filter_fixture()
    filtered, _ = tcga.filter_cohort(cohort, raw)
    assert list(filtered.columns) == list(cohort.columns)
    again, _ = tcga.filter_cohort(cohort, raw)
    pd.testing.assert_frame_equal(filtered, again)


def test_filter_cohort_length_mismatch_raises():
    cohort, raw = _filter_fixture()
    with pytest.raises(ValueError):
        tcga.filter_cohort(cohort, raw.iloc[:-1])


def test_no_filter_path_leaves_row_count_unchanged():
    # --no-filter simply does not apply filter_cohort.
    cohort, _raw = _filter_fixture()
    assert len(cohort) == 7


def test_print_filter_report_is_ascii(capsys):
    cohort, raw = _filter_fixture()
    _filtered, report = tcga.filter_cohort(cohort, raw)
    tcga.print_filter_report(report)
    out = capsys.readouterr().out
    assert out.isascii()
    assert "Final cohort rows: 4" in out


def test_script_source_is_ascii_only():
    assert SCRIPT_PATH.read_bytes().decode("utf-8").isascii()


def test_kps_uses_gdc_status_field_and_coerces_unknown():
    """GDC names the field karnofsky_performance_status (string-valued).

    Regression: the script once requested karnofsky_performance_score, which the
    API rejects, leaving KPS null for every case.
    """
    assert "follow_ups.karnofsky_performance_status" in tcga.FIELDS
    assert not any("karnofsky_performance_score" in f for f in tcga.FIELDS)

    case = {
        "submitter_id": "TCGA-XX-0001",
        "follow_ups": [
            {"days_to_follow_up": 754},
            {"karnofsky_performance_status": "60", "days_to_follow_up": None},
        ],
    }
    assert tcga.flatten_case(case)["kps_raw"] == "60"

    unknown = {"submitter_id": "TCGA-XX-0002",
               "follow_ups": [{"karnofsky_performance_status": "Unknown"}]}
    assert tcga.flatten_case(unknown)["kps_raw"] == "Unknown"
    # build_cohort turns the numeric string into a number and "Unknown" into NA.
    raw = pd.DataFrame([tcga.flatten_case(case), tcga.flatten_case(unknown)])
    kps = tcga.build_cohort(raw)["kps"]
    assert kps.iloc[0] == 60
    assert pd.isna(kps.iloc[1])


def test_uninformative_diagnosis_is_not_contradictory():
    """A blank or 'Unknown' diagnosis must not exclude a TCGA-GBM case.

    Regression: the original filter required an explicit 'glioblastoma'
    string, which discarded 8 patients who have RNA-seq data purely because
    their clinical diagnosis field was never filled in.
    """
    for keep in [None, float("nan"), "", "Unknown", "unknown", "Not Reported", "  "]:
        assert tcga.is_diagnosis_contradictory(keep) is False, keep
    for keep in ["Glioblastoma", "Glioblastoma, NOS", "GLIOBLASTOMA MULTIFORME"]:
        assert tcga.is_diagnosis_contradictory(keep) is False, keep
    for drop in ["Adenocarcinoma, NOS", "Oligodendroglioma", "Astrocytoma, NOS"]:
        assert tcga.is_diagnosis_contradictory(drop) is True, drop
    # The strict helper still exists and still means what it says.
    assert tcga.is_glioblastoma("Glioblastoma, NOS") is True
    assert tcga.is_glioblastoma("Unknown") is False
