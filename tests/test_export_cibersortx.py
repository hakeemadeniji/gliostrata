"""Unit tests for scripts/export_cibersortx.py (no network).

CIBERSORTx is a web tool, so a malformed upload is discovered only after a
manual round trip through the portal. These tests pin the format rules that
would otherwise fail silently or, worse, be accepted and return skewed
fractions.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "export_cibersortx.py"


def _load_module():
    """Load the export script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_export_cibersortx_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


export = _load_module()


def _fixture():
    """Four genes, three samples, with one duplicate symbol and one dead gene."""
    tpm = pd.DataFrame(
        {
            "S1": [10.0, 5.0, 0.0, 1.0, 100.0],
            "S2": [20.0, 1.0, 0.0, 2.0, 200.0],
            "S3": [30.0, 3.0, 0.0, 3.0, 300.0],
        },
        index=pd.Index(
            ["ENSG1.1", "ENSG2.2", "ENSG3.3", "ENSG4.4", "ENSG5.5"], name="gene_id"
        ),
    )
    annotation = pd.DataFrame(
        {
            "gene_id": ["ENSG1.1", "ENSG2.2", "ENSG3.3", "ENSG4.4", "ENSG5.5"],
            # GENEB appears twice: ENSG2 is the louder copy, ENSG4 the quiet one.
            "gene_name": ["GENEA", "GENEB", "DEADGENE", "GENEB", "GENEC"],
            "gene_type": [
                "protein_coding",
                "protein_coding",
                "protein_coding",
                "protein_coding",
                "lncRNA",
            ],
        }
    )
    return tpm, annotation


# ---------------------------------------------------------------------------
# symbol mapping and duplicate collapse
# ---------------------------------------------------------------------------

def test_map_symbols_replaces_index_and_drops_unmapped():
    tpm, annotation = _fixture()
    annotation.loc[annotation["gene_id"] == "ENSG3.3", "gene_name"] = None
    mapped = export.map_symbols(tpm, annotation)
    assert "DEADGENE" not in mapped.index
    assert set(mapped.index) == {"GENEA", "GENEB", "GENEC"}
    assert mapped.index.name == export.GENE_COLUMN_HEADER


def test_collapse_keeps_the_highest_expressed_duplicate():
    """CIBERSORTx rejects duplicate symbols; the louder copy must win.

    Averaging a real signal with a near-silent duplicate would understate the
    gene, which for a marker gene directly biases the returned fractions.
    """
    tpm, annotation = _fixture()
    mapped = export.map_symbols(tpm, annotation)
    collapsed, dropped = export.collapse_duplicate_symbols(mapped)
    assert dropped == 1
    assert not collapsed.index.duplicated().any()
    # ENSG2 (mean 3.0) beats ENSG4 (mean 2.0).
    assert collapsed.loc["GENEB", "S1"] == 5.0


def test_collapse_is_deterministic_and_sorted():
    tpm, annotation = _fixture()
    mapped = export.map_symbols(tpm, annotation)
    first, _ = export.collapse_duplicate_symbols(mapped)
    second, _ = export.collapse_duplicate_symbols(mapped.iloc[::-1])
    assert list(first.index) == sorted(first.index)
    pd.testing.assert_frame_equal(first, second)


def test_collapse_without_duplicates_is_a_passthrough():
    frame = pd.DataFrame({"S1": [1.0, 2.0]}, index=["B", "A"])
    collapsed, dropped = export.collapse_duplicate_symbols(frame)
    assert dropped == 0
    assert list(collapsed.index) == ["A", "B"]


# ---------------------------------------------------------------------------
# build_mixture
# ---------------------------------------------------------------------------

def test_all_zero_genes_are_dropped():
    tpm, annotation = _fixture()
    mixture, report = export.build_mixture(tpm, annotation, gene_types=None)
    assert "DEADGENE" not in mixture.index
    assert report["all_zero_genes_dropped"] == 1


def test_gene_type_filter_selects_only_requested_biotypes():
    tpm, annotation = _fixture()
    mixture, report = export.build_mixture(
        tpm, annotation, gene_types=("protein_coding",)
    )
    assert "GENEC" not in mixture.index  # lncRNA
    assert report["after_gene_type_filter"] == 4


def test_patient_filter_subsets_columns():
    tpm, annotation = _fixture()
    mixture, report = export.build_mixture(
        tpm, annotation, gene_types=None, patients=["S1", "S3", "NOT-A-SAMPLE"]
    )
    assert list(mixture.columns) == ["S1", "S3"]
    assert report["after_patient_filter"] == 2


def test_renormalization_puts_every_sample_on_a_common_total():
    """Filtering genes leaves each sample summing to a different total.

    In the real cohort that spread is severe and systematic: 121 of 284
    samples above the 0.3 misc_RNA fraction carry a median 50.0% of TPM in
    misc_RNA versus a median 0.11% in the other 163, a library-prep
    difference. (The 0.3 threshold is a report parameter, not a claim about
    the underlying distribution.) Deconvolution reads relative abundance, so
    an uncorrected mixture hands CIBERSORTx a batch artefact.
    """
    tpm, annotation = _fixture()
    mixture, report = export.build_mixture(
        tpm, annotation, gene_types=("protein_coding",), renormalize=True
    )
    sums = mixture.sum(axis=0)
    assert sums.min() == pytest.approx(1e6, rel=1e-6)
    assert sums.max() == pytest.approx(1e6, rel=1e-6)
    assert report["renormalized_to"] == 1e6
    # The pre-rescaling spread is recorded, not hidden.
    assert report["column_sum_min"] < report["column_sum_max"]


def test_renormalization_preserves_within_sample_ratios():
    """Rescaling must not change a sample's internal composition."""
    tpm, annotation = _fixture()
    raw, _ = export.build_mixture(
        tpm, annotation, gene_types=None, renormalize=False, decimals=10
    )
    scaled, _ = export.build_mixture(
        tpm, annotation, gene_types=None, renormalize=True, decimals=10
    )
    ratio_raw = raw.loc["GENEA", "S1"] / raw.loc["GENEB", "S1"]
    ratio_scaled = scaled.loc["GENEA", "S1"] / scaled.loc["GENEB", "S1"]
    assert ratio_raw == pytest.approx(ratio_scaled)


def test_no_renormalize_leaves_values_untouched():
    tpm, annotation = _fixture()
    mixture, report = export.build_mixture(
        tpm, annotation, gene_types=None, renormalize=False
    )
    assert "renormalized_to" not in report
    assert mixture.loc["GENEA", "S1"] == 10.0


def test_mixture_has_no_nan_or_negative_values():
    tpm, annotation = _fixture()
    mixture, report = export.build_mixture(tpm, annotation, gene_types=None)
    assert report["nan_values"] == 0
    assert report["negative_values"] == 0


def test_build_mixture_is_deterministic():
    tpm, annotation = _fixture()
    first, r1 = export.build_mixture(tpm, annotation, gene_types=None)
    second, r2 = export.build_mixture(tpm, annotation, gene_types=None)
    pd.testing.assert_frame_equal(first, second)
    assert r1 == r2


# ---------------------------------------------------------------------------
# file format -- the rules CIBERSORTx enforces on upload
# ---------------------------------------------------------------------------

def test_written_file_matches_the_cibersortx_contract(tmp_path):
    tpm, annotation = _fixture()
    mixture, _ = export.build_mixture(tpm, annotation, gene_types=None)
    path = export.write_mixture(mixture, tmp_path / "mix.txt")

    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "GeneSymbol"
    assert header[1:] == ["S1", "S2", "S3"]
    # Tab-delimited, not comma: a CSV is silently mis-parsed by the portal.
    assert "," not in lines[0]
    assert all(len(line.split("\t")) == len(header) for line in lines[1:])

    reloaded = pd.read_csv(path, sep="\t", index_col=0)
    assert not reloaded.index.duplicated().any()
    assert reloaded.isna().sum().sum() == 0


def test_values_stay_in_linear_space(tmp_path):
    """CIBERSORTx requires non-log input; log2 data is accepted but wrong."""
    tpm, annotation = _fixture()
    mixture, _ = export.build_mixture(
        tpm, annotation, gene_types=None, renormalize=False
    )
    # The fixture's largest value is 300; a log2 transform would cap near 8.2.
    assert mixture.values.max() == pytest.approx(300.0)


def test_script_source_is_ascii_only():
    """cp1252 consoles raise UnicodeEncodeError on non-ASCII stdout."""
    assert SCRIPT_PATH.read_bytes().decode("utf-8").isascii()
