"""Unit tests for scripts/download_cgga.py (no network access).

Every test runs against small synthetic fixtures. Nothing here touches the
real CGGA cohort, and nothing here fits a model or relates a feature to an
outcome -- CGGA is the project's pre-registered external test set
(docs/preregistration.md Sections 2 and 3.3).
"""

import importlib.util
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "download_cgga.py"


def _load_module():
    """Load the download script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_cgga_download_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


cgga = _load_module()


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

#: One synthetic row per interesting case. Column names use the CGGA
#: header spellings so normalise_clinical_columns is exercised too.
RAW_HEADER = [
    "CGGA_ID",
    "PRS_type",
    "Histology",
    "Grade",
    "Gender",
    "Age",
    "OS",
    "Censor (alive=0; dead=1)",
    "Radio_status (treated=1;un-treated=0)",
    "Chemo_status (TMZ treated=1;un-treated=0)",
    "IDH_mutation_status",
    "1p19q_codeletion_status",
    "MGMTp_methylation_status",
]

RAW_ROWS = [
    # keep: primary grade-IV GBM, dead, chemoradiation
    ["CGGA_001", "Primary", "GBM", "WHO IV", "Male", "60", "365", "1",
     "1", "1", "Wildtype", "Non-codel", "methylated"],
    # keep: alive (censored), monotherapy, unmethylated, trailing-space sex
    ["CGGA_002", "Primary", "GBM", "WHO IV", "Female ", "45", "730", "0",
     "1", "0", "Mutant", "Non-codel", "un-methylated"],
    # keep: no adjuvant therapy recorded as given
    ["CGGA_003", "Primary", "GBM", "WHO IV", "Male", "70", "90", "1",
     "0", "0", "Wildtype", "NA", "NA"],
    # keep: treatment arm unknown (one flag missing)
    ["CGGA_004", "Primary", "GBM", "WHO IV", "Female", "55", "200", "1",
     "NA", "1", "NA", "NA", "methylated"],
    # drop: recurrent GBM
    ["CGGA_005", "Recurrent", "rGBM", "WHO IV", "Male", "50", "300", "1",
     "1", "1", "Wildtype", "Non-codel", "methylated"],
    # drop: secondary GBM
    ["CGGA_006", "Secondary", "sGBM", "WHO IV", "Male", "41", "400", "1",
     "1", "1", "Mutant", "Non-codel", "methylated"],
    # drop: lower-grade glioma
    ["CGGA_007", "Primary", "A", "WHO II", "Female", "33", "1500", "0",
     "1", "1", "Mutant", "Codel", "methylated"],
    # drop: non-positive survival
    ["CGGA_008", "Primary", "GBM", "WHO IV", "Male", "64", "0", "1",
     "1", "1", "Wildtype", "Non-codel", "un-methylated"],
    # drop: nothing usable at all
    ["CGGA_009", "Primary", "GBM", "WHO IV", "NA", "NA", "NA", "NA",
     "NA", "NA", "NA", "NA", "NA"],
]


def _raw_frame():
    """Normalised raw frame, as read_clinical would produce it."""
    frame = pd.DataFrame(RAW_ROWS, columns=RAW_HEADER).astype("object")
    frame = frame.replace("NA", None)
    frame.columns = cgga.normalise_clinical_columns(frame.columns)
    frame.insert(1, "cgga_dataset", "mRNAseq_693")
    return frame


# --------------------------------------------------------------------------
# ASCII-only source
# --------------------------------------------------------------------------

def test_script_source_is_ascii_only():
    """The project has crashed on cp1252 consoles; stdout must stay ASCII."""
    raw_bytes = SCRIPT_PATH.read_bytes()
    try:
        raw_bytes.decode("ascii")
    except UnicodeDecodeError as exc:
        offending = raw_bytes[exc.start : exc.end]
        line = raw_bytes[: exc.start].count(b"\n") + 1
        pytest.fail(f"non-ASCII byte {offending!r} at line {line}")


def test_script_uses_ascii_arrows_not_unicode():
    text = SCRIPT_PATH.read_text(encoding="ascii")
    for bad in ("→", "←", "—", "–"):
        assert bad not in text


# --------------------------------------------------------------------------
# Column normalisation
# --------------------------------------------------------------------------

def test_normalise_clinical_columns_canonicalises_known_headers():
    out = cgga.normalise_clinical_columns(RAW_HEADER)
    assert out == [
        "cgga_id", "prs_type", "histology", "grade", "gender", "age_years",
        "os_days", "censor", "radio_status", "chemo_status",
        "idh_mutation_status", "codel_1p19q_status",
        "mgmtp_methylation_status",
    ]


def test_normalise_clinical_columns_handles_1p19q_case_variants():
    """The microarray batch capitalises "Codeletion"; the seq batches do not."""
    lower = cgga.normalise_clinical_columns(["1p19q_codeletion_status"])
    upper = cgga.normalise_clinical_columns(["1p19q_Codeletion_status"])
    assert lower == upper == ["codel_1p19q_status"]


def test_normalise_clinical_columns_keeps_unknown_headers():
    assert cgga.normalise_clinical_columns(["Brand_New_Field"]) == [
        "brand_new_field"
    ]


# --------------------------------------------------------------------------
# Sex mapping
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "value,expected",
    [
        ("Male", "M"),
        ("Female", "F"),
        ("male", "M"),
        ("FEMALE", "F"),
        ("Male ", "M"),
        ("Male  ", "M"),
        ("Female ", "F"),
    ],
)
def test_map_sex_known_values(value, expected):
    assert cgga.map_sex(value) == expected


@pytest.mark.parametrize("value", [None, "", "NA", "Unknown", 1, float("nan")])
def test_map_sex_unknown_is_missing(value):
    assert cgga.map_sex(value) is pd.NA


# --------------------------------------------------------------------------
# IDH / MGMT code mapping
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "value,expected", [("Mutant", 1), ("mutant", 0 + 1), ("Wildtype", 0),
                       ("wildtype", 0)]
)
def test_map_idh_known_values(value, expected):
    assert cgga.map_idh(value) == expected


@pytest.mark.parametrize("value", [None, "NA", "", "Unknown", "Mutantish"])
def test_map_idh_unknown_is_missing(value):
    assert cgga.map_idh(value) is pd.NA


def test_map_mgmt_methylated_is_one():
    assert cgga.map_mgmt("methylated") == 1


def test_map_mgmt_unmethylated_is_zero_not_one():
    """"un-methylated" contains "methylated": a substring match would invert
    every negative. This is the regression guard for that bug."""
    assert cgga.map_mgmt("un-methylated") == 0
    assert cgga.map_mgmt("unmethylated") == 0
    assert cgga.map_mgmt("UN-METHYLATED") == 0


@pytest.mark.parametrize("value", [None, "NA", "", "Unknown"])
def test_map_mgmt_unknown_is_missing(value):
    assert cgga.map_mgmt(value) is pd.NA


# --------------------------------------------------------------------------
# Subtype
# --------------------------------------------------------------------------

def test_map_subtype_passes_label_through():
    assert cgga.map_subtype(" Mesenchymal ") == "Mesenchymal"


@pytest.mark.parametrize("value", [None, "", "NA", "Unknown", 3])
def test_map_subtype_unknown_is_missing(value):
    assert cgga.map_subtype(value) is pd.NA


# --------------------------------------------------------------------------
# Survival: days -> months
# --------------------------------------------------------------------------

def test_days_to_months_uses_the_tcga_divisor():
    out = cgga.days_to_months(pd.Series(["365", "30.4375", "0"]))
    assert out.tolist() == [round(365 / 30.4375, 2), 1.0, 0.0]


def test_days_to_months_is_not_an_identity_or_a_year_conversion():
    """CGGA OS is in DAYS. A 730-day survival is ~24 months, not 730 and not 2."""
    assert cgga.days_to_months(pd.Series([730])).iloc[0] == pytest.approx(23.98, abs=0.01)


def test_days_to_months_non_numeric_is_missing():
    out = cgga.days_to_months(pd.Series(["NA", None, "abc"]))
    assert out.isna().all()


# --------------------------------------------------------------------------
# Censoring flag
# --------------------------------------------------------------------------

def test_map_censor_is_identity_not_inversion():
    """CGGA codes the EVENT in a column called "Censor": 1 = dead, 0 = alive,
    which already matches the project's vital_status convention."""
    assert cgga.map_censor("1") == 1
    assert cgga.map_censor("0") == 0
    assert cgga.map_censor(1) == 1
    assert cgga.map_censor(0) == 0


@pytest.mark.parametrize("value", [None, "NA", "", "2", -1, "unknown"])
def test_map_censor_unknown_is_missing(value):
    assert cgga.map_censor(value) is pd.NA


# --------------------------------------------------------------------------
# Treatment arm
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "radio,chemo,expected",
    [
        ("1", "1", "ChemoRadiation"),
        (1, 1, "ChemoRadiation"),
        ("1", "0", "MonoTherapy"),
        ("0", "1", "MonoTherapy"),
        ("0", "0", "NoAdjuvant"),
    ],
)
def test_classify_arm_known_combinations(radio, chemo, expected):
    assert cgga.classify_arm(radio, chemo) == expected


@pytest.mark.parametrize(
    "radio,chemo",
    [(None, "1"), ("1", None), (None, None), ("NA", "0"), ("2", "1"), ("1", "7")],
)
def test_classify_arm_missing_or_invalid_flag_is_missing(radio, chemo):
    """A missing radio flag with chemo=1 is ambiguous between MonoTherapy and
    ChemoRadiation; guessing would invent an exposure."""
    assert cgga.classify_arm(radio, chemo) is pd.NA


# --------------------------------------------------------------------------
# Filter predicates
# --------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [
    ("Primary", True), ("primary", True), (" Primary ", True),
    ("Recurrent", False), ("Secondary", False), (None, False), (1, False),
])
def test_is_primary(value, expected):
    assert cgga.is_primary(value) is expected


@pytest.mark.parametrize("value,expected", [
    ("GBM", True), ("gbm", True), (" GBM ", True),
    ("rGBM", False), ("sGBM", False), ("A", False), ("AA", False),
    (None, False), (1, False),
])
def test_is_gbm_histology_rejects_recurrent_and_secondary(value, expected):
    """rGBM and sGBM both contain "GBM"; a substring match would keep them."""
    assert cgga.is_gbm_histology(value) is expected


@pytest.mark.parametrize("value,expected", [
    ("WHO IV", True), ("who iv", True), (" WHO  IV ", True),
    ("WHO II", False), ("WHO III", False), (None, False),
])
def test_is_grade_iv(value, expected):
    assert cgga.is_grade_iv(value) is expected


# --------------------------------------------------------------------------
# build_cohort
# --------------------------------------------------------------------------

def test_build_cohort_emits_the_schema_columns_in_order():
    cohort = cgga.build_cohort(_raw_frame())
    assert list(cohort.columns) == cgga.COHORT_COLUMNS + ["cgga_dataset"]


def test_build_cohort_preserves_row_count_and_order():
    raw = _raw_frame()
    cohort = cgga.build_cohort(raw)
    assert len(cohort) == len(raw)
    assert cohort["patient_id"].tolist() == raw["cgga_id"].tolist()


def test_build_cohort_maps_every_field_of_the_first_row():
    cohort = cgga.build_cohort(_raw_frame())
    row = cohort.iloc[0]
    assert row["patient_id"] == "CGGA_001"
    assert row["age"] == 60
    assert row["sex"] == "M"
    assert row["idh_mutated"] == 0
    assert row["mgmt_methylated"] == 1
    assert row["vital_status"] == 1
    assert row["treatment_arm"] == "ChemoRadiation"
    assert row["overall_survival_months"] == pytest.approx(11.99, abs=0.01)


def test_build_cohort_leaves_unmappable_fields_null():
    cohort = cgga.build_cohort(_raw_frame())
    for field in cgga.UNMAPPABLE_FIELDS:
        assert cohort[field].isna().all(), field


def test_build_cohort_subtype_is_null_without_a_tcga_subtypes_column():
    cohort = cgga.build_cohort(_raw_frame())
    assert cohort["subtype"].isna().all()


def test_build_cohort_reads_tcga_subtypes_when_present():
    raw = _raw_frame()
    raw["tcga_subtypes"] = ["Classical"] * len(raw)
    cohort = cgga.build_cohort(raw)
    assert (cohort["subtype"] == "Classical").all()


def test_build_cohort_tolerates_a_missing_optional_column():
    raw = _raw_frame().drop(columns=["mgmtp_methylation_status"])
    cohort = cgga.build_cohort(raw)
    assert cohort["mgmt_methylated"].isna().all()


# --------------------------------------------------------------------------
# filter_cohort
# --------------------------------------------------------------------------

def test_filter_cohort_keeps_only_primary_grade_iv_gbm():
    raw = _raw_frame()
    cohort, report = cgga.filter_cohort(cgga.build_cohort(raw), raw)
    assert cohort["patient_id"].tolist() == [
        "CGGA_001", "CGGA_002", "CGGA_003", "CGGA_004"
    ]
    assert report["kept_rows"] == 4


def test_filter_cohort_report_counts_are_exact():
    raw = _raw_frame()
    _, report = cgga.filter_cohort(cgga.build_cohort(raw), raw)
    assert report["input_rows"] == 9
    # first-match attribution
    assert report["dropped_not_primary"] == 2          # 005 recurrent, 006 secondary
    assert report["dropped_non_glioblastoma"] == 1     # 007 astrocytoma
    assert report["dropped_not_grade_iv"] == 0         # 007 already attributed
    assert report["dropped_no_usable_data"] == 1       # 009
    assert report["dropped_nonpositive_survival"] == 1  # 008
    assert report["dropped_total"] == 5
    assert report["kept_rows"] == 4
    # independent matches, which overlap
    assert report["matched_not_primary"] == 2
    assert report["matched_non_glioblastoma"] == 3     # 005, 006, 007
    assert report["matched_not_grade_iv"] == 1         # 007


def test_filter_cohort_dropped_reasons_sum_to_dropped_total():
    raw = _raw_frame()
    _, report = cgga.filter_cohort(cgga.build_cohort(raw), raw)
    total = sum(report[f"dropped_{k}"] for k, _ in cgga.FILTER_REASONS)
    assert total == report["dropped_total"]
    assert report["kept_rows"] + report["dropped_total"] == report["input_rows"]


def test_filter_cohort_reason_order_is_fixed():
    assert [k for k, _ in cgga.FILTER_REASONS] == [
        "not_primary",
        "non_glioblastoma",
        "not_grade_iv",
        "no_usable_data",
        "nonpositive_survival",
    ]


def test_filter_cohort_rejects_misaligned_inputs():
    raw = _raw_frame()
    cohort = cgga.build_cohort(raw)
    with pytest.raises(ValueError, match="length mismatch"):
        cgga.filter_cohort(cohort.iloc[:3], raw)


def test_filter_cohort_is_deterministic_and_order_preserving():
    raw = _raw_frame()
    first, _ = cgga.filter_cohort(cgga.build_cohort(raw), raw)
    second, _ = cgga.filter_cohort(cgga.build_cohort(raw), raw)
    pd.testing.assert_frame_equal(first, second)
    assert first.index.is_monotonic_increasing


def test_print_filter_report_is_ascii(capsys):
    raw = _raw_frame()
    _, report = cgga.filter_cohort(cgga.build_cohort(raw), raw)
    cgga.print_filter_report(report)
    out = capsys.readouterr().out
    out.encode("ascii")
    assert "Final cohort rows: 4" in out


# --------------------------------------------------------------------------
# Deterministic ordering
# --------------------------------------------------------------------------

def test_concat_clinical_sorts_by_dataset_then_id():
    a = _raw_frame().iloc[[2, 0, 1]].copy()
    b = _raw_frame().iloc[[1, 0]].copy()
    b["cgga_dataset"] = "mRNAseq_325"
    combined = cgga.concat_clinical([a, b])
    assert combined["cgga_dataset"].tolist() == (
        ["mRNAseq_325"] * 2 + ["mRNAseq_693"] * 3
    )
    assert combined["cgga_id"].tolist() == [
        "CGGA_001", "CGGA_002", "CGGA_001", "CGGA_002", "CGGA_003"
    ]


def test_concat_clinical_input_order_does_not_matter():
    a = _raw_frame()
    b = _raw_frame().iloc[::-1].copy()
    b["cgga_dataset"] = "mRNAseq_325"
    pd.testing.assert_frame_equal(
        cgga.concat_clinical([a, b]), cgga.concat_clinical([b, a])
    )


def test_concat_clinical_rejects_empty_input():
    with pytest.raises(ValueError, match="no clinical frames"):
        cgga.concat_clinical([])


# --------------------------------------------------------------------------
# Counts assembly
# --------------------------------------------------------------------------

def _counts(samples, genes):
    return pd.DataFrame(
        {s: list(range(len(genes))) for s in samples},
        index=pd.Index(genes, name="gene_name"),
    )


def test_assemble_counts_intersects_genes_and_sorts_both_axes():
    merged = cgga.assemble_counts(
        {
            "mRNAseq_693": _counts(["S2", "S1"], ["TP53", "EGFR", "MKI67"]),
            "mRNAseq_325": _counts(["S4", "S3"], ["EGFR", "TP53", "EXTRA"]),
        }
    )
    assert list(merged.index) == ["EGFR", "TP53"]
    assert list(merged.columns) == ["S1", "S2", "S3", "S4"]


def test_assemble_counts_rejects_shared_sample_ids():
    with pytest.raises(ValueError, match="share sample ids"):
        cgga.assemble_counts(
            {
                "mRNAseq_693": _counts(["S1"], ["TP53"]),
                "mRNAseq_325": _counts(["S1"], ["TP53"]),
            }
        )


def test_assemble_counts_rejects_empty_input():
    with pytest.raises(ValueError, match="no counts matrices"):
        cgga.assemble_counts({})


def test_assemble_counts_is_order_independent():
    a = _counts(["S1"], ["TP53", "EGFR"])
    b = _counts(["S2"], ["EGFR", "TP53"])
    pd.testing.assert_frame_equal(
        cgga.assemble_counts({"x": a, "y": b}),
        cgga.assemble_counts({"y": b, "x": a}),
    )


# --------------------------------------------------------------------------
# Harmonisation facts
# --------------------------------------------------------------------------

def test_harmonisation_facts_records_symbol_identifier_and_combat_seq():
    facts = cgga.harmonisation_facts(_counts(["S1"], ["TP53", "EGFR"]))
    assert "SYMBOL" in facts["cgga_gene_identifier"]
    assert "Ensembl" in facts["tcga_gene_identifier"]
    assert facts["batch_correction_required"] is True
    assert facts["batch_correction_method"].startswith("ComBat-seq")
    assert facts["cgga_gene_count"] == 2
    assert facts["cgga_sample_count"] == 1
    assert "symbols_in_common" not in facts


def test_harmonisation_facts_computes_mapping_fraction(tmp_path):
    annot = pd.DataFrame(
        {
            "gene_id": ["ENSG1.1", "ENSG2.1", "ENSG3.1", "ENSG4.1"],
            "gene_name": ["TP53", "EGFR", "NOTINCGGA", "TP53"],
            "gene_type": ["protein_coding"] * 3 + ["lncRNA"],
        }
    )
    path = tmp_path / "gene_annotation.parquet"
    annot.to_parquet(path, index=False)
    facts = cgga.harmonisation_facts(_counts(["S1"], ["TP53", "EGFR"]), path)
    assert facts["tcga_gene_count"] == 4
    assert facts["symbols_in_common"] == 2
    assert facts["fraction_cgga_genes_mapped"] == 1.0
    assert facts["fraction_tcga_genes_mapped"] == 0.75
    assert facts["tcga_duplicate_symbol_rows"] == 1
    assert facts["tcga_protein_coding_genes"] == 3
    assert facts["tcga_protein_coding_mapped"] == 2


def test_harmonisation_facts_is_json_and_ascii_serialisable():
    import json

    facts = cgga.harmonisation_facts(_counts(["S1"], ["TP53"]))
    json.dumps(facts, indent=2).encode("ascii")


def test_print_harmonisation_report_is_ascii(capsys):
    cgga.print_harmonisation_report(
        cgga.harmonisation_facts(_counts(["S1"], ["TP53"]))
    )
    capsys.readouterr().out.encode("ascii")


# --------------------------------------------------------------------------
# Duplicate ids / overlap reporting
# --------------------------------------------------------------------------

def test_duplicate_id_report_counts_cross_batch_reuse():
    cohort = pd.DataFrame(
        {
            "patient_id": pd.array(
                ["CGGA_1", "CGGA_2", "CGGA_1", "CGGA_3"], dtype="string"
            ),
            "cgga_dataset": pd.array(
                ["mRNAseq_693", "mRNAseq_693", "mRNA_array_301",
                 "mRNA_array_301"],
                dtype="string",
            ),
        }
    )
    report = cgga.duplicate_id_report(cohort)
    assert report["rows"] == 4
    assert report["unique_patients"] == 3
    assert report["ids_appearing_more_than_once"] == 1
    assert report["rows_beyond_first_per_id"] == 1
    assert report["shared_between_batches"] == {
        "mRNA_array_301 & mRNAseq_693": 1
    }


def test_print_duplicate_id_report_is_ascii(capsys):
    cohort = pd.DataFrame(
        {
            "patient_id": pd.array(["CGGA_1"], dtype="string"),
            "cgga_dataset": pd.array(["mRNAseq_693"], dtype="string"),
        }
    )
    cgga.print_duplicate_id_report(cgga.duplicate_id_report(cohort))
    capsys.readouterr().out.encode("ascii")


def test_expression_overlap_counts_identifiers_only():
    cohort = pd.DataFrame(
        {
            "patient_id": pd.array(["A", "B", "C"], dtype="string"),
            "overall_survival_months": [12.0, None, 5.0],
            "vital_status": pd.array([1, 1, None], dtype="Int64"),
        }
    )
    counts = _counts(["A", "B", "Z"], ["TP53"])
    report = cgga.expression_overlap(cohort, counts)
    assert report["cohort_patients"] == 3
    assert report["with_expression"] == 2
    assert report["with_survival"] == 1
    assert report["with_expression_and_survival"] == 1
    assert report["expression_only_samples"] == 1


def test_expression_overlap_without_counts():
    cohort = pd.DataFrame(
        {
            "patient_id": pd.array(["A"], dtype="string"),
            "overall_survival_months": [12.0],
            "vital_status": pd.array([1], dtype="Int64"),
        }
    )
    report = cgga.expression_overlap(cohort, None)
    assert report["with_expression"] == 0
    assert report["with_expression_and_survival"] == 0


# --------------------------------------------------------------------------
# Download plumbing (no network)
# --------------------------------------------------------------------------

def test_file_url_includes_file_type_and_time():
    url = cgga.file_url(cgga.CGGA_FILES["mRNAseq_693_clinical"])
    assert url.startswith("https://www.cgga.org.cn/download?file=")
    assert "type=mRNAseq_693_clinical" in url
    assert "time=20200506" in url


def test_selected_keys_respects_dataset_and_skip_expression():
    assert cgga.selected_keys(["mRNAseq_693"], skip_expression=False) == [
        "mRNAseq_693_clinical",
        "mRNAseq_693_counts",
    ]
    assert cgga.selected_keys(["mRNAseq_693"], skip_expression=True) == [
        "mRNAseq_693_clinical"
    ]
    assert cgga.selected_keys(["mRNA_array_301"], skip_expression=False) == [
        "mRNA_array_301_clinical"
    ]


def test_selected_keys_is_registry_ordered():
    keys = cgga.selected_keys(
        ["mRNA_array_301", "mRNAseq_325", "mRNAseq_693"], skip_expression=True
    )
    assert keys == [
        "mRNAseq_693_clinical",
        "mRNAseq_325_clinical",
        "mRNA_array_301_clinical",
    ]


def test_expression_datasets_excludes_the_microarray_batch():
    assert "mRNA_array_301" not in cgga.EXPRESSION_DATASETS


def test_pending_keys_skips_cached_archives(tmp_path):
    entry = cgga.CGGA_FILES["mRNAseq_693_clinical"]
    keys = ["mRNAseq_693_clinical", "mRNAseq_325_clinical"]
    assert cgga.pending_keys(keys, tmp_path) == keys
    cgga.cache_path_for(tmp_path, entry).write_bytes(b"PK\x03\x04data")
    assert cgga.pending_keys(keys, tmp_path) == ["mRNAseq_325_clinical"]


def test_pending_keys_treats_a_zero_byte_file_as_missing(tmp_path):
    entry = cgga.CGGA_FILES["mRNAseq_693_clinical"]
    cgga.cache_path_for(tmp_path, entry).write_bytes(b"")
    assert cgga.pending_keys(["mRNAseq_693_clinical"], tmp_path) == [
        "mRNAseq_693_clinical"
    ]


def test_download_files_with_everything_cached_makes_no_request(tmp_path, capsys):
    for key in ["mRNAseq_693_clinical", "mRNAseq_325_clinical"]:
        cgga.cache_path_for(tmp_path, cgga.CGGA_FILES[key]).write_bytes(b"PK\x03\x04")
    report = cgga.download_files(
        ["mRNAseq_693_clinical", "mRNAseq_325_clinical"], tmp_path
    )
    assert report["downloaded"] == 0
    assert report["skipped_cached"] == 2
    capsys.readouterr().out.encode("ascii")


# --------------------------------------------------------------------------
# Archive extraction and parsing
# --------------------------------------------------------------------------

def _make_zip(path, members):
    with zipfile.ZipFile(path, "w") as zf:
        for name, content in members.items():
            zf.writestr(name, content)
    return path


def test_extract_text_member_ignores_macosx_resource_forks(tmp_path):
    archive = _make_zip(
        tmp_path / "a.zip",
        {
            "clinical.txt": "CGGA_ID\tAge\nCGGA_1\t60\n",
            "__MACOSX/._clinical.txt": "junk",
        },
    )
    out = cgga.extract_text_member(archive, tmp_path / "ex")
    assert out.name == "clinical.txt"
    assert "CGGA_1" in out.read_text()


def test_extract_text_member_rejects_ambiguous_archives(tmp_path):
    archive = _make_zip(
        tmp_path / "b.zip", {"one.txt": "a", "two.txt": "b"}
    )
    with pytest.raises(ValueError, match="exactly one .txt member"):
        cgga.extract_text_member(archive, tmp_path / "ex")


def test_read_clinical_handles_cr_only_line_endings(tmp_path):
    """CGGA's 2020-05-06 clinical files use classic-Mac CR-only newlines."""
    path = tmp_path / "clin.txt"
    path.write_bytes(
        b"CGGA_ID\tGender\tOS\tCensor (alive=0; dead=1)\r"
        b"CGGA_1\tMale\t365\t1\r"
        b"CGGA_2\tFemale \t730\t0\r"
    )
    frame = cgga.read_clinical(path, "mRNAseq_693")
    assert len(frame) == 2
    assert list(frame.columns) == [
        "cgga_id", "cgga_dataset", "gender", "os_days", "censor"
    ]
    assert frame["gender"].tolist() == ["Male", "Female"]
    assert (frame["cgga_dataset"] == "mRNAseq_693").all()


def test_read_counts_uses_symbols_and_rejects_duplicates(tmp_path):
    good = tmp_path / "good.txt"
    good.write_text("gene_name\tS1\tS2\nTP53\t1\t2\nEGFR\t3\t4\n")
    frame = cgga.read_counts(good)
    assert frame.index.name == "gene_name"
    assert list(frame.index) == ["TP53", "EGFR"]

    bad = tmp_path / "bad.txt"
    bad.write_text("gene_name\tS1\nTP53\t1\nTP53\t2\n")
    with pytest.raises(ValueError, match="duplicate gene symbols"):
        cgga.read_counts(bad)


# --------------------------------------------------------------------------
# Contamination guards
# --------------------------------------------------------------------------

def test_script_declares_cgga_as_an_external_test_cohort():
    text = SCRIPT_PATH.read_text(encoding="ascii")
    assert "EXTERNAL VALIDATION COHORT" in text


def test_script_contains_no_model_fitting_or_association_code():
    """A grep-level guard: no estimator, metric or correlation call may appear
    in an ingest script for the pre-registered external test set."""
    text = SCRIPT_PATH.read_text(encoding="ascii")
    forbidden = [
        "sklearn",
        "lifelines",
        "sksurv",
        "statsmodels",
        "concordance",
        "c_index",
        "brier",
        "roc_auc",
        "logrank",
        "CoxPH",
        ".fit(",
        ".corr(",
        "pearsonr",
        "spearmanr",
        "groupby(\"vital_status\")",
    ]
    present = [token for token in forbidden if token in text]
    assert present == [], f"contamination-risk tokens in ingest script: {present}"
