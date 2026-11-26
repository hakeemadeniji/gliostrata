"""Unit tests for scripts/download_tcga_mutations.py (no network access).

What is at stake here
---------------------
A false IDH-wildtype call is not a small error. A patient who was never
sequenced, written down as "wildtype", silently inflates the wildtype group and
deflates the mutant fraction -- and under WHO 2021 the IDH-mutant cases are the
ones that may not be glioblastoma at all. So the tests below are organised
around three questions:

  1. Is a genuine hotspot called, and is a look-alike NOT called?
  2. Is "no mutation data" kept distinct from "sequenced, no hotspot"?
  3. Is the per-patient table deterministic and correctly typed?

The fixtures are synthetic MAFs written in the real GDC layout: seven ``#``
comment lines, a 140-column-style header, tab-separated rows, optionally
gzip-compressed.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import sys
import tarfile
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "download_tcga_mutations.py"
)


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_tcga_mutations_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


mut = _load_module()

MAF_COMMENT_LINES = [
    "#version gdc-1.0.0",
    "#annotation.spec gdc-2.0.0-ali",
    "#contigs chr1,chr2,chr3",
    "#sort.order BarcodesAndCoordinates",
    "#filedate 20220517",
    "#normal.aliquot 00000000-0000-0000-0000-000000000000",
    "#tumor.aliquot 11111111-1111-1111-1111-111111111111",
]

MAF_COLUMNS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "HGVSp_Short", "Protein_position", "Amino_acids",
]

SAMPLE_BARCODE = "TCGA-12-0819-01A-01D-0182-08"


def _variant(gene, classification, hgvsp, position, amino_acids):
    """One synthetic MAF row; non-tested columns get placeholder values."""
    row = {c: "" for c in MAF_COLUMNS}
    row.update({
        "Hugo_Symbol": gene,
        "Variant_Classification": classification,
        "Tumor_Sample_Barcode": SAMPLE_BARCODE,
        "HGVSp_Short": hgvsp,
        "Protein_position": position,
        "Amino_acids": amino_acids,
        "Chromosome": "chr2" if gene == "IDH1" else "chr15",
        "Start_Position": 209113112 if gene == "IDH1" else 90631934,
    })
    return row


def _write_maf(path: Path, rows, comment=True, compress=False, barcode=None) -> Path:
    """Write a MAF in the GDC layout: comments, header, rows; maybe gzipped."""
    lines = (MAF_COMMENT_LINES if comment else []) + ["\t".join(MAF_COLUMNS)]
    for row in rows:
        if barcode is not None:
            row = {**row, "Tumor_Sample_Barcode": barcode}
        lines.append("\t".join(str(row[c]) for c in MAF_COLUMNS))
    text = ("\n".join(lines) + "\n").encode("ascii")
    if compress:
        path.write_bytes(gzip.compress(text))
    else:
        path.write_bytes(text)
    return path


# ---------------------------------------------------------------------------
# Barcode parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "barcode,patient,code,vial",
    [
        ("TCGA-12-0819-01A-01D-0182-08", "TCGA-12-0819", "01", "A"),
        ("TCGA-12-0819-01A", "TCGA-12-0819", "01", "A"),
        ("TCGA-02-0003-02A-01D-0182-08", "TCGA-02-0003", "02", "A"),
        ("TCGA-02-0003-10A-01D-0182-08", "TCGA-02-0003", "10", "A"),
        ("TCGA-12-0819", "TCGA-12-0819", None, None),
        ("tcga-12-0819-01b-01r-1849-01", "TCGA-12-0819", "01", "B"),
    ],
)
def test_parse_barcode(barcode, patient, code, vial):
    parsed = mut.parse_barcode(barcode)
    assert parsed == {"patient_id": patient, "sample_type_code": code, "vial": vial}


@pytest.mark.parametrize("bad", [None, 7, "", "TCGA-XX", "TCGA", "12-0819-01A"])
def test_parse_barcode_never_guesses(bad):
    assert set(mut.parse_barcode(bad).values()) == {None}


def test_parse_barcode_keeps_the_patient_but_flags_an_unknown_sample_code():
    """An unknown sample-type code is not invented into a known one.

    The patient id is still recovered (it is well-formed), but the code stays
    "99", is labelled "Other", and is therefore dropped by the default
    primary-tumour filter rather than being called on.
    """
    parsed = mut.parse_barcode("TCGA-12-0819-99A-01D-0182-08")
    assert parsed == {"patient_id": "TCGA-12-0819", "sample_type_code": "99", "vial": "A"}
    assert mut.sample_type_label("99") == "Other (99)"
    assert "99" not in mut.DEFAULT_SAMPLE_TYPE_CODES


def test_flatten_file_record_extracts_patient_and_aliquot():
    hit = {
        "file_id": "f1",
        "file_name": "x.maf.gz",
        "file_size": 10,
        "md5sum": "abc",
        "analysis": {"workflow_type": "Aliquot Ensemble Somatic Variant Merging and Masking"},
        "cases": [{
            "case_id": "c1",
            "submitter_id": "TCGA-12-0819",
            "samples": [{
                "submitter_id": "TCGA-12-0819-01A",
                "sample_type": "Primary Tumor",
                "portions": [{"analytes": [{"aliquots": [{
                    "submitter_id": "TCGA-12-0819-01A-01D-0182-08"}]}]}],
            }],
        }],
    }
    row = mut.flatten_file_record(hit)
    assert row["patient_id"] == "TCGA-12-0819"
    assert row["sample_type_code"] == "01"
    assert row["aliquot_barcode"] == "TCGA-12-0819-01A-01D-0182-08"
    assert row["case_submitter_id"] == "TCGA-12-0819"


def _gdc_hit(file_id, patient, samples):
    """One GDC file record whose case links to `samples` in the given order."""
    return {
        "file_id": file_id, "file_name": f"{file_id}.maf.gz", "file_size": 10,
        "md5sum": "abc",
        "analysis": {"workflow_type": "w"},
        "cases": [{
            "case_id": f"case-{patient}", "submitter_id": patient,
            "samples": [{
                "submitter_id": f"{patient}-{sample_code}A",
                "sample_type": sample_type,
                "portions": [{"analytes": [{"aliquots": [{
                    "submitter_id": f"{patient}-{sample_code}A-01D-0182-08"}]}]}],
            } for sample_code, sample_type in samples],
        }],
    }


NORMAL_FIRST = [("10", "Blood Derived Normal"), ("01", "Primary Tumor")]
TUMOUR_FIRST = [("01", "Primary Tumor"), ("10", "Blood Derived Normal")]


@pytest.mark.parametrize("samples", [NORMAL_FIRST, TUMOUR_FIRST], ids=["normal-first", "tumour-first"])
def test_flatten_file_record_always_attributes_the_file_to_the_tumour(samples):
    """REGRESSION (found on the real GDC metadata, 2026-11-26).

    A masked MAF is linked to the tumour aliquot AND to the matched normal.
    For 239 of the 464 TCGA-GBM MAF files the NORMAL sample is listed first,
    so reading ``samples[0]`` labelled the file "Blood Derived Normal" and the
    primary-tumour filter then threw away mutation data for most of the
    cohort. The sample order in the metadata must not matter.
    """
    row = mut.flatten_file_record(_gdc_hit("f1", "TCGA-12-0819", samples))
    assert row["sample_type_code"] == "01"
    assert row["gdc_sample_type"] == "Primary Tumor"
    assert row["aliquot_barcode"] == "TCGA-12-0819-01A-01D-0182-08"
    assert row["sample_barcode"] == "TCGA-12-0819-01A"
    assert row["patient_id"] == "TCGA-12-0819"


def test_flatten_file_record_keeps_a_recurrent_tumour_labelled_as_recurrent():
    """No primary tumour linked: the file is still a tumour file, but code 02."""
    row = mut.flatten_file_record(
        _gdc_hit("f1", "TCGA-12-0819", [("10", "Blood Derived Normal"), ("02", "Recurrent Tumor")])
    )
    assert row["sample_type_code"] == "02"
    # ...and is therefore excluded by the default primary-only filter.
    frame = pd.DataFrame([row])
    kept, _report = mut.select_sample_types(frame)
    assert kept.empty


def test_flatten_file_record_prefers_primary_over_recurrent_when_both_are_linked():
    row = mut.flatten_file_record(_gdc_hit(
        "f1", "TCGA-12-0819",
        [("02", "Recurrent Tumor"), ("01", "Primary Tumor"), ("10", "Blood Derived Normal")],
    ))
    assert row["sample_type_code"] == "01"


def test_flatten_file_record_ordering_is_independent_of_metadata_sample_order():
    a = mut.flatten_file_record(_gdc_hit("f1", "TCGA-12-0819", NORMAL_FIRST))
    b = mut.flatten_file_record(_gdc_hit("f1", "TCGA-12-0819", TUMOUR_FIRST))
    assert a == b


# ---------------------------------------------------------------------------
# Sample-type selection and duplicate resolution
# ---------------------------------------------------------------------------

def _manifest():
    return pd.DataFrame([
        {"file_id": "f-dup-b", "patient_id": "TCGA-B", "sample_type_code": "01",
         "aliquot_barcode": "TCGA-B-01B-01D-01", "sample_barcode": "TCGA-B-01B"},
        {"file_id": "f-a-rec", "patient_id": "TCGA-A", "sample_type_code": "02",
         "aliquot_barcode": "TCGA-A-02A-01D-01", "sample_barcode": "TCGA-A-02A"},
        {"file_id": "f-a-prim", "patient_id": "TCGA-A", "sample_type_code": "01",
         "aliquot_barcode": "TCGA-A-01A-01D-01", "sample_barcode": "TCGA-A-01A"},
        {"file_id": "f-norm", "patient_id": "TCGA-C", "sample_type_code": "10",
         "aliquot_barcode": "TCGA-C-10A-01D-01", "sample_barcode": "TCGA-C-10A"},
        {"file_id": "f-c-prim", "patient_id": "TCGA-C", "sample_type_code": "01",
         "aliquot_barcode": "TCGA-C-01A-01D-01", "sample_barcode": "TCGA-C-01A"},
        {"file_id": "f-dup-a", "patient_id": "TCGA-B", "sample_type_code": "01",
         "aliquot_barcode": "TCGA-B-01A-01D-01", "sample_barcode": "TCGA-B-01A"},
    ])


def test_select_sample_types_keeps_only_primary_tumour():
    kept, report = mut.select_sample_types(_manifest())
    assert set(kept["sample_type_code"]) == {"01"}
    assert "f-norm" not in set(kept["file_id"]) and "f-a-rec" not in set(kept["file_id"])
    assert report["dropped_files"] == 2 and report["kept_files"] == 4
    assert report["breakdown"] == {"01": 4, "02": 1, "10": 1}


def test_select_sample_types_counts_unparsed_barcodes_separately():
    frame = _manifest()
    frame.loc[0, "sample_type_code"] = None
    _kept, report = mut.select_sample_types(frame)
    assert report["breakdown"]["unparsed"] == 1


def test_resolve_duplicate_aliquots_picks_lowest_vial_and_prefers_primary():
    deduped, report = mut.resolve_duplicate_aliquots(_manifest())
    picked = deduped.set_index("patient_id")["file_id"].to_dict()
    assert picked == {"TCGA-A": "f-a-prim", "TCGA-B": "f-dup-a", "TCGA-C": "f-c-prim"}
    assert report["patients"] == 3
    assert report["patients_with_duplicates"] == 3
    assert report["files_dropped_as_duplicate"] == 3
    assert list(deduped["patient_id"]) == ["TCGA-A", "TCGA-B", "TCGA-C"]


def test_resolve_duplicate_aliquots_is_input_order_independent():
    first, _ = mut.resolve_duplicate_aliquots(_manifest())
    second, _ = mut.resolve_duplicate_aliquots(_manifest().iloc[::-1])
    pd.testing.assert_frame_equal(first, second)


def test_resolve_duplicate_aliquots_on_empty_manifest():
    deduped, report = mut.resolve_duplicate_aliquots(_manifest().iloc[0:0])
    assert deduped.empty and report == {
        "input_files": 0, "patients": 0,
        "patients_with_duplicates": 0, "files_dropped_as_duplicate": 0,
    }

# ---------------------------------------------------------------------------
# Download layer (tar extraction; no network)
# ---------------------------------------------------------------------------

def test_cache_path_and_pending_ids(tmp_path):
    cache = tmp_path / "_cache"
    cache.mkdir()
    (cache / "have.maf.gz").write_bytes(b"x")
    assert mut.cache_path_for(cache, "abc") == cache / "abc.maf.gz"
    assert mut.pending_file_ids(["have", "missing"], cache) == ["missing"]
    # A zero-byte file counts as not downloaded, so a truncated fetch resumes.
    (cache / "empty.maf.gz").write_bytes(b"")
    assert mut.pending_file_ids(["have", "empty"], cache) == ["empty"]


def test_download_batch_unpacks_a_multi_file_tar(tmp_path, monkeypatch):
    payloads = {"id-1": b"fake-maf-1", "id-2": b"fake-maf-2"}

    class _Response:
        status_code = 200
        text = ""

        def __init__(self, content):
            self.content = content

    def fake_post(url, json=None, headers=None, timeout=None):
        assert json == {"ids": ["id-1", "id-2"]}
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            info = tarfile.TarInfo("MANIFEST.txt")
            info.size = 4
            tar.addfile(info, io.BytesIO(b"meta"))
            for file_id, blob in payloads.items():
                info = tarfile.TarInfo(f"{file_id}/{file_id}.masked.maf.gz")
                info.size = len(blob)
                tar.addfile(info, io.BytesIO(blob))
        return _Response(buf.getvalue())

    monkeypatch.setattr(mut.requests, "post", fake_post)
    cache = tmp_path / "_cache"
    assert mut.download_batch(list(payloads), cache) == 2
    # MANIFEST.txt at the archive root is ignored; members land under their id.
    assert (cache / "id-1.maf.gz").read_bytes() == payloads["id-1"]
    assert (cache / "id-2.maf.gz").read_bytes() == payloads["id-2"]
    assert not (cache / "MANIFEST.txt").exists()


def test_download_batch_single_id_writes_the_bare_file(tmp_path, monkeypatch):
    class _Response:
        status_code = 200
        text = ""
        content = b"bare-maf-bytes"

    monkeypatch.setattr(mut.requests, "post", lambda *a, **k: _Response())
    cache = tmp_path / "_cache"
    assert mut.download_batch(["only"], cache) == 1
    assert (cache / "only.maf.gz").read_bytes() == b"bare-maf-bytes"


def test_download_batch_raises_on_http_error_after_retries(tmp_path, monkeypatch):
    class _Response:
        status_code = 500
        text = "upstream exploded"

    calls = []

    def fake_post(*args, **kwargs):
        calls.append(1)
        return _Response()

    monkeypatch.setattr(mut.requests, "post", fake_post)
    with pytest.raises(RuntimeError, match="GDC data error 500"):
        mut.download_batch(["x"], tmp_path, max_attempts=3, retry_wait=0)
    assert len(calls) == 3  # it retried before giving up


def test_download_batch_retries_then_succeeds_on_a_transient_500(tmp_path, monkeypatch):
    """GDC answered 500 to a real 50-file batch; a retry must recover it."""
    class _Response:
        def __init__(self, status, content=b""):
            self.status_code = status
            self.text = ""
            self.content = content

    answers = iter([_Response(500), _Response(200, b"recovered")])
    monkeypatch.setattr(mut.requests, "post", lambda *a, **k: next(answers))
    monkeypatch.setattr(mut.time, "sleep", lambda seconds: None)
    cache = tmp_path / "_cache"
    assert mut.download_batch(["only"], cache, max_attempts=3, retry_wait=5) == 1
    assert (cache / "only.maf.gz").read_bytes() == b"recovered"


def test_download_files_isolates_a_bad_file_and_keeps_going(tmp_path, monkeypatch):
    """One unreachable id must cost that file only, and be reported."""
    def fake_batch(file_ids, cache_dir, timeout=600, max_attempts=3, retry_wait=2.0):
        for file_id in file_ids:
            if file_id == "bad":
                raise RuntimeError("GDC data error 500: internal server error")
            mut.cache_path_for(cache_dir, file_id).write_bytes(b"z")
        return len(file_ids)

    monkeypatch.setattr(mut, "download_batch", fake_batch)
    report = mut.download_files(["good1", "bad", "good2"], tmp_path, batch_size=3)
    assert report["downloaded"] == 2
    assert report["failed"] == ["bad"]
    assert (tmp_path / "good1.maf.gz").exists() and (tmp_path / "good2.maf.gz").exists()


def test_download_files_reports_cached_and_fetched(tmp_path, monkeypatch):
    calls = []

    def fake_batch(file_ids, cache_dir, timeout=600):
        calls.append(list(file_ids))
        for file_id in file_ids:
            mut.cache_path_for(cache_dir, file_id).write_bytes(b"z")
        return len(file_ids)

    monkeypatch.setattr(mut, "download_batch", fake_batch)
    mut.cache_path_for(tmp_path, "cached").write_bytes(b"z")
    report = mut.download_files(["cached", "a", "b"], tmp_path, batch_size=2)
    assert calls == [["a", "b"]]
    assert report["requested"] == 3 and report["skipped_cached"] == 1
    assert report["downloaded"] == 2


# ---------------------------------------------------------------------------
# MAF reading
# ---------------------------------------------------------------------------

def test_read_maf_plain_and_gzipped_are_identical(tmp_path):
    rows = [_variant("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H"),
            _variant("TP53", "Missense_Mutation", "p.R175H", "175", "R/H")]
    plain = mut.read_maf(_write_maf(tmp_path / "a.maf", rows))
    zipped = mut.read_maf(_write_maf(tmp_path / "b.maf.gz", rows, compress=True))
    pd.testing.assert_frame_equal(plain, zipped)
    assert len(plain) == 2
    assert "Hugo_Symbol" in plain.columns


def test_read_maf_header_only_is_an_empty_frame_with_columns(tmp_path):
    """Sequenced with zero somatic variants is a valid wildtype sample."""
    frame = mut.read_maf(_write_maf(tmp_path / "empty.maf", []))
    assert len(frame) == 0
    assert set(mut.MAF_REQUIRED_COLUMNS) <= set(frame.columns)


def test_read_maf_tolerates_a_missing_comment_block(tmp_path):
    rows = [_variant("IDH2", "Missense_Mutation", "p.R172K", "172", "R/K")]
    frame = mut.read_maf(_write_maf(tmp_path / "nocomment.maf", rows, comment=False))
    assert len(frame) == 1


def test_read_maf_rejects_a_file_without_a_header(tmp_path):
    path = tmp_path / "junk.txt"
    path.write_text("not,a,maf\n", encoding="ascii")
    with pytest.raises(ValueError, match="no MAF header"):
        mut.read_maf(path)


def test_read_maf_rejects_a_header_missing_required_columns(tmp_path):
    path = tmp_path / "short.maf"
    path.write_text("Hugo_Symbol\tTumor_Sample_Barcode\nIDH1\tTCGA\n", encoding="ascii")
    with pytest.raises(ValueError, match="missing required column"):
        mut.read_maf(path)

# ---------------------------------------------------------------------------
# Hotspot calling: the look-alikes matter as much as the hits
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "gene,classification,hgvsp,position,amino_acids,expected_label",
    [
        ("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H", "IDH1 p.R132H"),
        ("IDH1", "Missense_Mutation", "p.R132C", "132", "R/C", "IDH1 p.R132C"),
        ("IDH1", "Missense_Mutation", "p.R132G", "132", "R/G", "IDH1 p.R132G"),
        ("IDH1", "Missense_Mutation", "p.R132S", "132", "R/S", "IDH1 p.R132S"),
        ("IDH2", "Missense_Mutation", "p.R172K", "172", "R/K", "IDH2 p.R172K"),
        ("IDH2", "Missense_Mutation", "p.R172M", "172", "R/M", "IDH2 p.R172M"),
        ("IDH2", "Missense_Mutation", "p.R172W", "172", "R/W", "IDH2 p.R172W"),
        # In-frame indels that include the arginine are still neomorphic-competent.
        ("IDH1", "In_Frame_Del", "p.R132_P135delinsHQ", "132-135", "RPPP/HQ",
         "IDH1 p.R132_P135delinsHQ"),
        ("IDH1", "In_Frame_Ins", "p.R132_R133insH", "132", "R/RH", "IDH1 p.R132_R133insH"),
    ],
)
def test_genuine_hotspots_are_called(gene, classification, hgvsp, position, amino_acids,
                                     expected_label):
    row = _variant(gene, classification, hgvsp, position, amino_acids)
    assert mut.is_idh_hotspot_variant(row)
    assert mut.hotspot_label(row) == expected_label


@pytest.mark.parametrize(
    "gene,classification,hgvsp,position,amino_acids,why",
    [
        # A non-hotspot IDH1 variant: same gene, wrong codon.
        ("IDH1", "Missense_Mutation", "p.R100Q", "100", "R/Q", "wrong codon"),
        ("IDH1", "Missense_Mutation", "p.G97D", "97", "G/D", "wrong codon"),
        ("IDH1", "Missense_Mutation", "p.V178I", "178", "V/I", "wrong codon"),
        # Right codon number, but the reference residue is not arginine: the
        # annotation is on a different numbering and this is not the hotspot.
        ("IDH1", "Missense_Mutation", "p.G132D", "132", "G/D", "reference is not R"),
        ("IDH2", "Missense_Mutation", "p.K172N", "172", "K/N", "reference is not R"),
        # Right codon, but a class that cannot create the neomorphic enzyme.
        ("IDH1", "Nonsense_Mutation", "p.R132*", "132", "R/*", "nonsense"),
        ("IDH1", "Frame_Shift_Del", "p.R132fs", "132", "R/X", "frameshift"),
        ("IDH1", "Silent", "", "132", "R/R", "silent"),
        ("IDH2", "Splice_Site", "", "172", "", "splice site"),
        # Hotspot codon number in a gene that is not IDH1/IDH2.
        ("TP53", "Missense_Mutation", "p.R132H", "132", "R/H", "wrong gene"),
        # Hotspot residues swapped between the two genes.
        ("IDH1", "Missense_Mutation", "p.R172K", "172", "R/K", "IDH1 codon 172"),
        ("IDH2", "Missense_Mutation", "p.R132H", "132", "R/H", "IDH2 codon 132"),
    ],
)
def test_look_alikes_are_not_called(gene, classification, hgvsp, position, amino_acids, why):
    row = _variant(gene, classification, hgvsp, position, amino_acids)
    assert not mut.is_idh_hotspot_variant(row), why


def test_call_idh_for_maf_labels_the_hotspot_and_counts_all_variants(tmp_path):
    rows = [
        _variant("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H"),
        _variant("TP53", "Missense_Mutation", "p.R175H", "175", "R/H"),
        _variant("EGFR", "In_Frame_Ins", "p.EE702_EE703insKIP", "702", "EE/EEKIP"),
        # Two reads of the same hotspot must collapse to ONE variant label.
        _variant("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H"),
    ]
    call = mut.call_idh_for_maf(mut.read_maf(_write_maf(tmp_path / "s.maf", rows)))
    assert call["n_somatic_variants"] == 4
    assert call["hotspot_variants"] == ["IDH1 p.R132H"]


def test_call_idh_for_maf_on_a_variant_free_sample_is_wildtype(tmp_path):
    call = mut.call_idh_for_maf(mut.read_maf(_write_maf(tmp_path / "e.maf", [])))
    assert call == {"n_somatic_variants": 0, "hotspot_variants": []}


def test_call_idh_falls_back_to_protein_position_and_amino_acids():
    """Some MAFs leave HGVSp_Short empty; the maf2vcf-style pair must work."""
    row = _variant("IDH1", "Missense_Mutation", "", "132", "R/H")
    assert mut.parse_protein_change("", "132", "R/H") == ("R", 132, "H")
    assert mut.is_idh_hotspot_variant(row)
    assert mut.hotspot_label(row) == "IDH1 p.R132H"


def test_parse_protein_change_handles_missing_and_junk():
    for junk in (None, float("nan"), "", "p.", "p.R13", "p.?132?", "R132", 7):
        assert mut.parse_protein_change(junk) == (None, None, None)
    # Without a prefix.
    assert mut.parse_protein_change("R132H") == ("R", 132, "H")
    assert mut.parse_protein_change("p.R132del") == ("R", 132, "del")

# ---------------------------------------------------------------------------
# The three-state call table
# ---------------------------------------------------------------------------

def _record(pid, hits, n=5):
    return {"patient_id": pid, "hotspot_variants": hits, "n_somatic_variants": n,
            "sample_barcode": f"{pid}-01A", "aliquot_barcode": f"{pid}-01A-01D-01",
            "file_id": f"f-{pid}"}


def test_missing_is_never_coerced_to_wildtype():
    """The central contract: 1 = mutant, 0 = sequenced-and-no-hotspot, NA = no data."""
    records = [_record("TCGA-01-0001", ["IDH1 p.R132H"]),
               _record("TCGA-01-0002", [], n=0)]
    calls = mut.build_calls(records, ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003"])
    row = calls.set_index("patient_id")

    assert row.loc["TCGA-01-0001", "idh_mutated"] == 1
    assert row.loc["TCGA-01-0001", "call_basis"] == mut.BASIS_MUTANT
    assert row.loc["TCGA-01-0001", "idh_variants"] == "IDH1 p.R132H"

    # Sequenced, MAF had zero variants: a genuine wildtype, NOT missing.
    assert row.loc["TCGA-01-0002", "idh_mutated"] == 0
    assert row.loc["TCGA-01-0002", "call_basis"] == mut.BASIS_WILDTYPE
    assert row.loc["TCGA-01-0002", "n_somatic_variants"] == 0

    # Never sequenced: nullable, and every dependent field is nullable too.
    assert pd.isna(row.loc["TCGA-01-0003", "idh_mutated"])
    assert row.loc["TCGA-01-0003", "call_basis"] == mut.BASIS_MISSING
    assert pd.isna(row.loc["TCGA-01-0003", "n_somatic_variants"])
    assert pd.isna(row.loc["TCGA-01-0003", "n_idh_hotspot_variants"])

    assert list(calls.columns) == mut.CALL_COLUMNS
    assert list(calls["patient_id"]) == sorted(calls["patient_id"])


def test_calls_survive_a_parquet_round_trip_without_becoming_zero(tmp_path):
    """The join step reads this table back; NA must still be NA afterwards."""
    calls = mut.build_calls([_record("TCGA-01-0001", [])], ["TCGA-01-0001", "TCGA-01-0009"])
    path = tmp_path / "idh_calls.parquet"
    calls.to_parquet(path, index=False)
    back = pd.read_parquet(path).set_index("patient_id")
    assert str(back["idh_mutated"].dtype) == "Int64"
    assert back.loc["TCGA-01-0001", "idh_mutated"] == 0
    assert pd.isna(back.loc["TCGA-01-0009", "idh_mutated"])
    assert back.loc["TCGA-01-0009", "call_basis"] == mut.BASIS_MISSING


def test_two_samples_for_one_patient_merge_to_mutant_if_any_is_mutant():
    records = [_record("TCGA-01-0001", [], n=12), _record("TCGA-01-0001", ["IDH2 p.R172K"], n=40)]
    calls = mut.build_calls(records)
    assert len(calls) == 1
    row = calls.iloc[0]
    assert row["idh_mutated"] == 1
    assert row["idh_variants"] == "IDH2 p.R172K"
    assert row["n_somatic_variants"] == 40


def test_build_calls_is_order_independent_and_dedups_labels():
    records = [_record("TCGA-02-0002", ["IDH1 p.R132C", "IDH1 p.R132H"], n=9),
               _record("TCGA-01-0001", ["IDH1 p.R132H", "IDH1 p.R132H"], n=3)]
    first = mut.build_calls(records, ["TCGA-03-0003"])
    second = mut.build_calls(records[::-1], ["TCGA-03-0003"])
    pd.testing.assert_frame_equal(first, second)
    assert first.iloc[0]["idh_variants"] == "IDH1 p.R132H"
    assert first.iloc[0]["n_idh_hotspot_variants"] == 1


def test_build_calls_with_nothing_at_all_returns_the_contract_columns():
    calls = mut.build_calls([])
    assert calls.empty and list(calls.columns) == mut.CALL_COLUMNS


def test_build_calls_cohort_only_lists_everyone_as_missing():
    calls = mut.build_calls([], ["TCGA-01-0001", "TCGA-01-0002"])
    assert calls["idh_mutated"].isna().all()
    assert (calls["call_basis"] == mut.BASIS_MISSING).all()

# ---------------------------------------------------------------------------
# Summaries, coverage and the printed report
# ---------------------------------------------------------------------------

def _summary_input():
    records = [_record(f"TCGA-01-{i:04d}", ["IDH1 p.R132H"] if i <= 2 else [], n=10)
               for i in range(1, 41)]
    return mut.build_calls(records, ["TCGA-01-0000"] + [f"TCGA-01-{i:04d}" for i in range(41, 46)])


def test_summarize_calls_counts_the_three_states_separately():
    summary = mut.summarize_calls(_summary_input())
    assert summary["idh_mutant"] == 2
    assert summary["idh_wildtype"] == 38
    assert summary["no_mutation_data"] == 6
    assert summary["patients_with_mutation_data"] == 40
    assert summary["patients_in_table"] == 46
    # The denominator is TESTED patients; missing patients never dilute the rate.
    assert summary["mutant_fraction_of_tested"] == pytest.approx(2 / 40)
    assert summary["mutant_fraction_of_all_rows"] == pytest.approx(2 / 46)
    assert summary["distinct_hotspot_variants"] == ["IDH1 p.R132H"]
    assert summary["mutant_patient_ids"] == ["TCGA-01-0001", "TCGA-01-0002"]
    assert summary["red_flag_high_mutant_rate"] is False


def test_summarize_calls_raises_a_red_flag_on_an_implausibly_high_rate():
    """Primary GBM sits at single-digit percent; 12/40 is a red flag, not a result."""
    records = [_record(f"TCGA-01-{i:04d}", ["IDH1 p.R132H"] if i <= 12 else [], n=10)
               for i in range(1, 41)]
    summary = mut.summarize_calls(mut.build_calls(records))
    assert summary["idh_mutant"] == 12
    assert summary["red_flag_high_mutant_rate"] is True


def test_summarize_calls_on_an_empty_table_has_no_division_error():
    summary = mut.summarize_calls(mut.build_calls([]))
    assert summary["idh_mutant"] == 0 and summary["mutant_fraction_of_tested"] is None
    assert summary["red_flag_high_mutant_rate"] is False


def test_coverage_counts_how_many_of_a_patient_set_gain_a_call():
    calls = mut.build_calls(
        [_record("TCGA-01-0001", ["IDH1 p.R132H"]), _record("TCGA-01-0002", [])],
        ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003"],
    )
    row = mut.coverage(calls, ["TCGA-01-0001", "TCGA-01-0003", "TCGA-01-0002", "TCGA-09-0009"],
                       "expression patients")
    assert row == {
        "set": "expression patients", "patients": 4, "with_idh_call": 2,
        "idh_mutant": 1, "idh_wildtype": 1, "still_missing": 2,
    }


def test_parquet_column_names_reads_the_schema_not_an_empty_selection(tmp_path):
    """Regression: read_parquet(columns=[]) yields NO columns at all."""
    path = tmp_path / "tpm.parquet"
    pd.DataFrame({"TCGA-01-0001": [1.0, 2.0], "TCGA-01-0003": [3.0, 4.0]}).to_parquet(path)
    assert mut.parquet_column_names(path) == ["TCGA-01-0001", "TCGA-01-0003"]
    assert list(pd.read_parquet(path, columns=[]).columns) == []  # the trap


def test_expression_patient_ids_never_count_the_gene_column(tmp_path):
    """Regression: the real tpm.parquet has a `gene_id` column, and it was
    counted as a patient, overstating the expression set by one."""
    path = tmp_path / "tpm.parquet"
    pd.DataFrame({"gene_id": ["EGFR", "IDH1"], "TCGA-01-0001": [1.0, 2.0],
                  "TCGA-01-0003": [3.0, 4.0]}).to_parquet(path, index=False)
    ids, ignored = mut.expression_patient_ids(path)
    assert ids == ["TCGA-01-0001", "TCGA-01-0003"]
    assert ignored == ["gene_id"]

    # A file with no barcode column at all is the wrong file: refuse it.
    other = tmp_path / "clinical_cohort.parquet"
    pd.DataFrame({"patient_id": ["TCGA-01-0001"], "age": [61]}).to_parquet(other, index=False)
    with pytest.raises(RuntimeError, match="not a gene x sample matrix"):
        mut.expression_patient_ids(other)


def test_jsonable_makes_a_report_strict_json():
    payload = {"a": float("nan"), "b": [float("inf"), 1.5], "c": {"d": None},
               "e": __import__("numpy").int64(3), "f": __import__("numpy").float64(0.5)}
    cleaned = mut._jsonable(payload)
    text = json.dumps(cleaned, allow_nan=False)
    assert json.loads(text) == {"a": None, "b": [None, 1.5], "c": {"d": None}, "e": 3, "f": 0.5}


def test_printed_report_is_ascii_and_names_the_who_flag(capsys):
    summary = mut.summarize_calls(_summary_input())
    summary["red_flag_high_mutant_rate"] = True
    mut.print_report(summary, [mut.coverage(mut.build_calls([], ["TCGA-01-0001"]),
                                            ["TCGA-01-0001"], "clinical cohort")])
    out = capsys.readouterr().out
    assert out.isascii()
    assert "no mutation data" in out and "NOT wildtype" in out
    assert "WHO 2021 FLAG -- 2 IDH-mutant case(s): TCGA-01-0001, TCGA-01-0002" in out
    assert "RED FLAG" in out
    assert "IDH1 p.R132H" in out


def test_hotspot_definition_matches_the_who_2021_hotspots():
    assert mut.IDH_HOTSPOT_CODONS == {"IDH1": 132, "IDH2": 172}
    assert mut.HOTSPOT_REFERENCE_RESIDUE == "R"
    assert set(mut.HOTSPOT_VARIANT_CLASSES) == {
        "Missense_Mutation", "In_Frame_Ins", "In_Frame_Del"
    }


def test_script_source_is_ascii_only():
    """cp1252 consoles raise UnicodeEncodeError on non-ASCII stdout."""
    assert SCRIPT_PATH.read_bytes().decode("utf-8").isascii()

# ---------------------------------------------------------------------------
# End to end, with the network replaced by fixtures
# ---------------------------------------------------------------------------

def test_maf_tumour_patient_reads_the_barcode_inside_the_file(tmp_path):
    """Cross-check: the MAF's own Tumor_Sample_Barcode must agree with metadata."""
    rows = [_variant("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H"),
            _variant("TP53", "Missense_Mutation", "p.R175H", "175", "R/H")]
    barcode = "TCGA-44-2655-01A-01D-A273-08"
    maf = mut.read_maf(_write_maf(tmp_path / "a.maf", rows, barcode=barcode))
    assert mut.maf_tumour_patient(maf) == "TCGA-44-2655"
    # No rows -> no barcode to check; that is not an error.
    empty = mut.read_maf(_write_maf(tmp_path / "e.maf", []))
    assert mut.maf_tumour_patient(empty) is None


def test_main_runs_end_to_end_offline(tmp_path, monkeypatch, capsys):
    out_dir = tmp_path / "mutations"

    def fake_meta(project_id=mut.PROJECT_ID, page_size=500):
        # The normal sample is listed FIRST on purpose: the real metadata does that.
        return [_gdc_hit("f1", "TCGA-01-0001", NORMAL_FIRST),
                _gdc_hit("f2", "TCGA-01-0002", TUMOUR_FIRST)]

    def fake_download(file_ids, cache_dir, batch_size=50):
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        _write_maf(cache_dir / "f1.maf.gz",
                   [_variant("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H"),
                    _variant("TP53", "Missense_Mutation", "p.R175H", "175", "R/H")],
                   compress=True, barcode="TCGA-01-0001-01A-01D-0182-08")
        _write_maf(cache_dir / "f2.maf.gz", [], compress=True)
        return {"requested": 2, "skipped_cached": 0, "downloaded": 2, "failed": [],
                "elapsed_seconds": 0.0}

    monkeypatch.setattr(mut, "fetch_file_metadata", fake_meta)
    monkeypatch.setattr(mut, "download_files", fake_download)

    cohort = tmp_path / "clinical_cohort.parquet"
    pd.DataFrame({"patient_id": ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003"]}).to_parquet(cohort)
    expression = tmp_path / "tpm.parquet"
    # Laid out like the real matrix: gene x sample written with index=False, so
    # its first column is gene_id and NOT a patient (it was once counted as one).
    pd.DataFrame({"gene_id": ["ENSG00000141510"], "TCGA-01-0001": [1.0],
                  "TCGA-01-0003": [2.0]}).to_parquet(expression, index=False)

    monkeypatch.setattr(sys, "argv", [
        "download_tcga_mutations.py", "--out", str(out_dir),
        "--cohort", str(cohort), "--expression", str(expression),
    ])
    mut.main()

    calls = pd.read_parquet(out_dir / "idh_calls.parquet").set_index("patient_id")
    assert calls.loc["TCGA-01-0001", "idh_variants"] == "IDH1 p.R132H"
    assert calls.loc["TCGA-01-0002", "idh_mutated"] == 0
    assert calls.loc["TCGA-01-0002", "n_somatic_variants"] == 0
    assert pd.isna(calls.loc["TCGA-01-0003", "idh_mutated"])
    assert calls.loc["TCGA-01-0001", "aliquot_barcode"] == "TCGA-01-0001-01A-01D-0182-08"

    report = json.loads((out_dir / "call_report.json").read_text(encoding="ascii"))
    assert report["summary"]["idh_mutant"] == 1
    assert report["summary"]["no_mutation_data"] == 1
    assert report["coverage"] == [
        {"set": "clinical cohort", "patients": 3, "with_idh_call": 2, "idh_mutant": 1,
         "idh_wildtype": 1, "still_missing": 1},
        {"set": "with expression", "patients": 2, "with_idh_call": 1, "idh_mutant": 1,
         "idh_wildtype": 0, "still_missing": 1},
        {"set": "cohort AND expression", "patients": 2, "with_idh_call": 1, "idh_mutant": 1,
         "idh_wildtype": 0, "still_missing": 1},
    ]
    assert report["expression_non_sample_columns"] == ["gene_id"]
    assert (out_dir / "file_manifest.parquet").exists()

    printed = capsys.readouterr().out
    assert printed.isascii()
    assert "WHO 2021 FLAG -- 1 IDH-mutant case(s): TCGA-01-0001" in printed
    # The clinical cohort file itself is never rewritten by this script.
    assert len(pd.read_parquet(cohort)) == 3


def test_main_refuses_a_maf_whose_own_barcode_belongs_to_another_patient(tmp_path, monkeypatch):
    """Metadata and file content disagree -> stop, do not guess which is right."""
    def fake_meta(project_id=mut.PROJECT_ID, page_size=500):
        return [_gdc_hit("f1", "TCGA-01-0001", NORMAL_FIRST)]

    def fake_download(file_ids, cache_dir, batch_size=50):
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        _write_maf(cache_dir / "f1.maf.gz",
                   [_variant("IDH1", "Missense_Mutation", "p.R132H", "132", "R/H")],
                   compress=True, barcode="TCGA-99-9999-01A-01D-0182-08")
        return {"requested": 1, "skipped_cached": 0, "downloaded": 1, "failed": [],
                "elapsed_seconds": 0.0}

    monkeypatch.setattr(mut, "fetch_file_metadata", fake_meta)
    monkeypatch.setattr(mut, "download_files", fake_download)
    monkeypatch.setattr(sys, "argv", [
        "download_tcga_mutations.py", "--out", str(tmp_path / "x"), "--no-cohort",
    ])
    with pytest.raises(RuntimeError, match="belongs to patient TCGA-99-9999"):
        mut.main()


def test_main_dry_run_downloads_and_writes_nothing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mut, "fetch_file_metadata",
                        lambda project_id=mut.PROJECT_ID, page_size=500: [
                            _gdc_hit("f1", "TCGA-01-0001", NORMAL_FIRST)])

    def must_not_download(*args, **kwargs):
        raise AssertionError("dry-run must not download")

    monkeypatch.setattr(mut, "download_files", must_not_download)
    out_dir = tmp_path / "dry"
    monkeypatch.setattr(sys, "argv", [
        "download_tcga_mutations.py", "--out", str(out_dir), "--dry-run",
    ])
    mut.main()
    printed = capsys.readouterr().out
    assert "--dry-run -> nothing downloaded" in printed
    assert not out_dir.exists()


def test_main_refuses_metadata_with_an_unparseable_barcode(tmp_path, monkeypatch):
    hit = _gdc_hit("f1", "TCGA-01-0001", TUMOUR_FIRST)
    hit["cases"][0]["samples"] = []
    hit["cases"][0]["submitter_id"] = "also-not-a-barcode"
    hit["cases"][0]["submitter_id"] = "also-not-a-barcode"
    monkeypatch.setattr(mut, "fetch_file_metadata",
                        lambda project_id=mut.PROJECT_ID, page_size=500: [hit])
    monkeypatch.setattr(sys, "argv", [
        "download_tcga_mutations.py", "--out", str(tmp_path / "x"), "--dry-run",
    ])
    with pytest.raises(RuntimeError, match="unparseable TCGA barcode"):
        mut.main()

