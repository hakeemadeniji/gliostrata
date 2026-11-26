"""Unit tests for scripts/download_tcga_expression.py (no network access)."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "download_tcga_expression.py"
)


def _load_module():
    """Load the download script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_tcga_expression_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


expr = _load_module()


# --------------------------------------------------------------------------
# Synthetic STAR-format fixtures
# --------------------------------------------------------------------------

GENES = [
    ("ENSG00000000003.15", "TSPAN6", "protein_coding"),
    ("ENSG00000000005.6", "TNMD", "protein_coding"),
    ("ENSG00000000419.13", "DPM1", "protein_coding"),
    ("ENSG00000284332.1", "MIR1302-2", "miRNA"),
]


def _write_star_tsv(path: Path, counts, tpms, genes=GENES, comment=True) -> Path:
    """Write a minimal GDC STAR counts TSV: comment, header, 4 summary rows, genes."""
    lines = []
    if comment:
        lines.append("# gene-model: GENCODE v36")
    lines.append("\t".join(expr.STAR_COLUMNS))
    for name, value in zip(expr.STAR_SUMMARY_ROWS, [11, 22, 33, 44]):
        lines.append(
            "\t".join([name, "", "", str(value), str(value), str(value), "", "", ""])
        )
    for (gene_id, gene_name, gene_type), count, tpm in zip(genes, counts, tpms):
        lines.append(
            "\t".join(
                [
                    gene_id,
                    gene_name,
                    gene_type,
                    str(count),
                    str(count),
                    str(count),
                    f"{tpm:.4f}",
                    f"{tpm:.4f}",
                    f"{tpm:.4f}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


# --------------------------------------------------------------------------
# read_star_counts
# --------------------------------------------------------------------------

def test_read_star_counts_drops_the_four_summary_rows(tmp_path):
    path = _write_star_tsv(tmp_path / "s1.tsv", [10, 20, 30, 40], [1.0, 2.0, 3.0, 4.0])
    frame = expr.read_star_counts(path)
    assert len(frame) == len(GENES)
    assert not set(frame["gene_id"]) & set(expr.STAR_SUMMARY_ROWS)
    assert list(frame["gene_id"]) == [g[0] for g in GENES]


def test_read_star_counts_keeps_all_star_columns(tmp_path):
    path = _write_star_tsv(tmp_path / "s1.tsv", [1, 2, 3, 4], [0.1, 0.2, 0.3, 0.4])
    frame = expr.read_star_counts(path)
    for column in expr.STAR_COLUMNS:
        assert column in frame.columns


def test_read_star_counts_rejects_missing_columns(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("gene_id\tgene_name\nENSG1\tA\n", encoding="ascii")
    with pytest.raises(ValueError, match="missing STAR columns"):
        expr.read_star_counts(path)


def test_read_star_counts_rejects_duplicate_gene_ids(tmp_path):
    genes = list(GENES) + [GENES[0]]
    path = _write_star_tsv(
        tmp_path / "dup.tsv", [1, 2, 3, 4, 5], [1.0] * 5, genes=genes
    )
    with pytest.raises(ValueError, match="duplicate gene_id"):
        expr.read_star_counts(path)


def test_read_star_counts_tolerates_missing_comment_line(tmp_path):
    path = _write_star_tsv(
        tmp_path / "s1.tsv", [1, 2, 3, 4], [1.0] * 4, comment=False
    )
    assert len(expr.read_star_counts(path)) == len(GENES)


# --------------------------------------------------------------------------
# parse_barcode
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "barcode,patient,code,vial",
    [
        ("TCGA-12-0819-01A-01R-1849-01", "TCGA-12-0819", "01", "A"),
        ("TCGA-06-0208-11A-01R-1849-01", "TCGA-06-0208", "11", "A"),
        ("TCGA-02-0099-02B-02R-A96S-41", "TCGA-02-0099", "02", "B"),
        ("TCGA-32-5222-01A", "TCGA-32-5222", "01", "A"),
        ("TCGA-32-5222", "TCGA-32-5222", None, None),
        ("tcga-32-5222-01a", "TCGA-32-5222", "01", "A"),
    ],
)
def test_parse_barcode(barcode, patient, code, vial):
    parsed = expr.parse_barcode(barcode)
    assert parsed["patient_id"] == patient
    assert parsed["sample_type_code"] == code
    assert parsed["vial"] == vial


@pytest.mark.parametrize("bad", [None, 12, "", "not-a-barcode", float("nan")])
def test_parse_barcode_rejects_garbage(bad):
    parsed = expr.parse_barcode(bad)
    assert parsed == {"patient_id": None, "sample_type_code": None, "vial": None}


def test_parse_barcode_patient_id_matches_cohort_form():
    # clinical_cohort.parquet uses submitter_id, e.g. "TCGA-12-0819".
    assert expr.parse_barcode("TCGA-12-0819-01A-01R-1849-01")["patient_id"] == (
        "TCGA-12-0819"
    )


def test_sample_type_labels():
    assert expr.sample_type_label("01") == "Primary Tumor"
    assert expr.sample_type_label("11") == "Solid Tissue Normal"
    assert expr.sample_type_label("02") == "Recurrent Tumor"
    assert "99" in expr.sample_type_label("99")
    assert expr.sample_type_label(None) == "Unknown"


# --------------------------------------------------------------------------
# flatten_file_record
# --------------------------------------------------------------------------

def _hit(file_id, case, aliquot, sample_type="Primary Tumor"):
    return {
        "file_id": file_id,
        "file_name": f"{file_id}.rna_seq.augmented_star_gene_counts.tsv",
        "file_size": 4200000,
        "md5sum": "abc",
        "analysis": {"workflow_type": "STAR - Counts"},
        "cases": [
            {
                "case_id": f"uuid-{case}",
                "submitter_id": case,
                "samples": [
                    {
                        "submitter_id": aliquot[:16],
                        "sample_type": sample_type,
                        "portions": [
                            {"analytes": [{"aliquots": [{"submitter_id": aliquot}]}]}
                        ],
                    }
                ],
            }
        ],
    }


def test_flatten_file_record_extracts_aliquot_and_patient():
    row = expr.flatten_file_record(
        _hit("f1", "TCGA-12-0819", "TCGA-12-0819-01A-01R-1849-01")
    )
    assert row["patient_id"] == "TCGA-12-0819"
    assert row["aliquot_barcode"] == "TCGA-12-0819-01A-01R-1849-01"
    assert row["sample_type_code"] == "01"
    assert row["workflow_type"] == "STAR - Counts"


def test_flatten_file_record_handles_missing_nesting():
    row = expr.flatten_file_record({"file_id": "f1", "cases": []})
    assert row["patient_id"] is None
    assert row["aliquot_barcode"] is None


# --------------------------------------------------------------------------
# filter_sample_types
# --------------------------------------------------------------------------

def _manifest():
    hits = [
        _hit("f1", "TCGA-12-0819", "TCGA-12-0819-01A-01R-1849-01"),
        _hit("f2", "TCGA-12-0819", "TCGA-12-0819-01B-01R-1849-01"),
        _hit("f3", "TCGA-06-0208", "TCGA-06-0208-01A-01R-1849-01"),
        _hit(
            "f4", "TCGA-06-0208", "TCGA-06-0208-11A-01R-1849-01", "Solid Tissue Normal"
        ),
        _hit("f5", "TCGA-02-0099", "TCGA-02-0099-02A-01R-1849-01", "Recurrent Tumor"),
    ]
    return pd.DataFrame([expr.flatten_file_record(h) for h in hits])


def test_filter_sample_types_keeps_primary_tumor_by_default():
    kept, report = expr.filter_sample_types(
        _manifest(), expr.DEFAULT_SAMPLE_TYPE_CODES
    )
    assert list(kept["file_id"]) == ["f1", "f2", "f3"]
    assert report["breakdown"] == {"01": 3, "02": 1, "11": 1}
    assert report["kept_files"] == 3
    assert report["dropped_files"] == 2


def test_filter_sample_types_can_retain_others():
    kept, report = expr.filter_sample_types(_manifest(), ["01", "02", "11"])
    assert report["kept_files"] == 5
    assert report["dropped_files"] == 0
    assert len(kept) == 5


# --------------------------------------------------------------------------
# resolve_duplicate_aliquots
# --------------------------------------------------------------------------

def test_resolve_duplicate_aliquots_is_deterministic_and_picks_lowest():
    manifest, _ = expr.filter_sample_types(
        _manifest(), expr.DEFAULT_SAMPLE_TYPE_CODES
    )
    deduped, report = expr.resolve_duplicate_aliquots(manifest)
    assert report["input_files"] == 3
    assert report["patients"] == 2
    assert report["patients_with_duplicates"] == 1
    assert report["files_dropped_as_duplicate"] == 1
    # -01A- beats -01B- for TCGA-12-0819.
    picked = dict(zip(deduped["patient_id"], deduped["aliquot_barcode"]))
    assert picked["TCGA-12-0819"] == "TCGA-12-0819-01A-01R-1849-01"
    # Deterministic under row shuffling.
    shuffled = manifest.iloc[::-1].reset_index(drop=True)
    again, again_report = expr.resolve_duplicate_aliquots(shuffled)
    pd.testing.assert_frame_equal(deduped, again)
    assert report == again_report


def test_resolve_duplicate_aliquots_sorts_by_patient_id():
    manifest, _ = expr.filter_sample_types(_manifest(), ["01", "02", "11"])
    deduped, _ = expr.resolve_duplicate_aliquots(manifest)
    assert list(deduped["patient_id"]) == sorted(deduped["patient_id"])


def test_resolve_duplicate_aliquots_empty_manifest():
    deduped, report = expr.resolve_duplicate_aliquots(_manifest().iloc[:0])
    assert deduped.empty
    assert report["patients"] == 0
    assert report["patients_with_duplicates"] == 0


# --------------------------------------------------------------------------
# assemble_matrices
# --------------------------------------------------------------------------

def _three_samples(tmp_path):
    specs = {
        "TCGA-06-0208": ([5, 6, 7, 8], [0.5, 0.6, 0.7, 0.8]),
        "TCGA-02-0099": ([1, 2, 3, 4], [0.1, 0.2, 0.3, 0.4]),
        "TCGA-12-0819": ([9, 10, 11, 12], [0.9, 1.0, 1.1, 1.2]),
    }
    per_sample = {}
    for name, (counts, tpms) in specs.items():
        path = _write_star_tsv(tmp_path / f"{name}.tsv", counts, tpms)
        per_sample[name] = expr.read_star_counts(path)
    return per_sample


def test_assemble_matrices_shape_and_ordering(tmp_path):
    counts, tpm, annot = expr.assemble_matrices(_three_samples(tmp_path))
    assert counts.shape == (len(GENES), 3)
    assert tpm.shape == (len(GENES), 3)
    # Sample columns sorted; genes sorted by gene_id.
    assert list(counts.columns) == ["TCGA-02-0099", "TCGA-06-0208", "TCGA-12-0819"]
    assert list(counts.index) == sorted(g[0] for g in GENES)
    assert counts.index.name == "gene_id"


def test_assemble_matrices_counts_and_tpm_stay_aligned(tmp_path):
    counts, tpm, _ = expr.assemble_matrices(_three_samples(tmp_path))
    assert list(counts.index) == list(tpm.index)
    assert list(counts.columns) == list(tpm.columns)
    # TCGA-02-0099 had counts 1..4 / tpm 0.1..0.4 in GENES order.
    assert counts.loc["ENSG00000000003.15", "TCGA-02-0099"] == 1
    assert tpm.loc["ENSG00000000003.15", "TCGA-02-0099"] == pytest.approx(0.1)
    assert counts.loc["ENSG00000000419.13", "TCGA-12-0819"] == 11
    assert tpm.loc["ENSG00000000419.13", "TCGA-12-0819"] == pytest.approx(1.1)


def test_assemble_matrices_value_dtypes_are_numeric(tmp_path):
    counts, tpm, _ = expr.assemble_matrices(_three_samples(tmp_path))
    assert all(str(d).startswith(("int", "float")) for d in counts.dtypes)
    assert all(str(d).startswith("float") for d in tpm.dtypes)


def test_assemble_matrices_annotation_table(tmp_path):
    _counts, _tpm, annot = expr.assemble_matrices(_three_samples(tmp_path))
    assert list(annot.columns) == expr.GENE_ANNOTATION_COLUMNS
    assert len(annot) == len(GENES)
    assert list(annot["gene_id"]) == sorted(g[0] for g in GENES)
    assert set(annot["gene_type"]) == {"protein_coding", "miRNA"}


def test_assemble_matrices_is_byte_stable_across_reruns(tmp_path):
    per_sample = _three_samples(tmp_path)
    first = expr.assemble_matrices(per_sample)
    shuffled = {k: per_sample[k] for k in reversed(list(per_sample))}
    second = expr.assemble_matrices(shuffled)
    for a, b in zip(first, second):
        pd.testing.assert_frame_equal(a, b)


def test_assemble_matrices_rejects_mismatched_gene_sets(tmp_path):
    per_sample = _three_samples(tmp_path)
    odd_path = _write_star_tsv(
        tmp_path / "odd.tsv",
        [1, 2, 3],
        [0.1, 0.2, 0.3],
        genes=GENES[:2] + [("ENSG99999999.1", "FAKE", "protein_coding")],
    )
    per_sample["TCGA-99-9999"] = expr.read_star_counts(odd_path)
    with pytest.raises(ValueError, match="differs from reference"):
        expr.assemble_matrices(per_sample)


def test_assemble_matrices_rejects_empty_input():
    with pytest.raises(ValueError, match="no samples"):
        expr.assemble_matrices({})


# --------------------------------------------------------------------------
# cache / resume
# --------------------------------------------------------------------------

def test_pending_file_ids_skips_cached_files(tmp_path):
    expr.cache_path_for(tmp_path, "f1").write_text("x", encoding="ascii")
    (tmp_path / "f2.tsv").write_text("", encoding="ascii")  # empty -> refetch
    assert expr.pending_file_ids(["f1", "f2", "f3"], tmp_path) == ["f2", "f3"]


def test_cache_path_for_is_file_id_based(tmp_path):
    assert expr.cache_path_for(tmp_path, "abc").name == "abc.tsv"


# --------------------------------------------------------------------------
# reporting / clinical overlap
# --------------------------------------------------------------------------

def test_reports_are_ascii_only(capsys, tmp_path):
    manifest, st_report = expr.filter_sample_types(
        _manifest(), expr.DEFAULT_SAMPLE_TYPE_CODES
    )
    deduped, dup_report = expr.resolve_duplicate_aliquots(manifest)
    expr.print_metadata_report(_manifest(), total_files=5)
    expr.print_sample_type_report(st_report)
    expr.print_duplicate_report(dup_report)
    clinical = tmp_path / "clinical.parquet"
    pd.DataFrame({"patient_id": ["TCGA-12-0819", "TCGA-99-9999"]}).to_parquet(
        clinical, index=False
    )
    overlap = expr.clinical_overlap(deduped["patient_id"], clinical)
    expr.print_overlap_report(overlap)
    out = capsys.readouterr().out
    assert out.isascii()
    assert "->" in out


def test_clinical_overlap_counts(tmp_path):
    clinical = tmp_path / "clinical.parquet"
    pd.DataFrame(
        {"patient_id": ["TCGA-12-0819", "TCGA-06-0208", "TCGA-99-9999"]}
    ).to_parquet(clinical, index=False)
    report = expr.clinical_overlap(["TCGA-12-0819", "TCGA-11-1111"], clinical)
    assert report["clinical_rows"] == 3
    assert report["overlap"] == 1
    assert report["expression_only"] == 1
    assert report["clinical_only"] == 2


def test_clinical_overlap_missing_file(tmp_path):
    report = expr.clinical_overlap(["TCGA-12-0819"], tmp_path / "nope.parquet")
    assert report["overlap"] is None


# --------------------------------------------------------------------------
# GDC query shape
# --------------------------------------------------------------------------

def test_build_file_filters_targets_star_counts():
    filters = expr.build_file_filters("TCGA-GBM", "STAR - Counts")
    flat = {
        item["content"]["field"]: item["content"]["value"]
        for item in filters["content"]
    }
    assert flat["cases.project.project_id"] == ["TCGA-GBM"]
    assert flat["data_category"] == ["Transcriptome Profiling"]
    assert flat["data_type"] == ["Gene Expression Quantification"]
    assert flat["analysis.workflow_type"] == ["STAR - Counts"]


def test_endpoints_are_the_documented_gdc_urls():
    assert expr.GDC_FILES_URL == "https://api.gdc.cancer.gov/files"
    assert expr.GDC_DATA_URL == "https://api.gdc.cancer.gov/data"


def test_script_source_is_ascii_only():
    assert SCRIPT_PATH.read_bytes().decode("utf-8").isascii()
