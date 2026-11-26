"""Unit tests for scripts/download_tcga_methylation.py (no network access).

What is at stake here
---------------------
A false "unmethylated" call is worse than no call: it assigns a patient the
worse-prognosis, less-temozolomide-responsive group on the strength of data
that does not exist. So, as with the IDH script, the tests are organised
around three questions:

  1. Are the two MGMT probes read correctly out of a real-format beta file,
     and is a probe that is present-but-NA treated as missing?
  2. Is "not callable" kept distinct from "unmethylated"?
  3. Does the call depend on the threshold in the way the docstring claims,
     and is it labelled as what it is -- NOT MGMT-STP27?

The fixtures are synthetic beta-value files in the real GDC layout: headerless
``probe<TAB>beta`` text with ``NA`` for unmeasurable probes.
"""

from __future__ import annotations

import importlib.util
import io
import json
import sys
import tarfile
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "download_tcga_methylation.py"
)


def _load_module():
    """Load the script by path (scripts/ is not a package)."""
    assert SCRIPT_PATH.is_file(), f"script not found: {SCRIPT_PATH}"
    spec = importlib.util.spec_from_file_location(
        "gliostrata_tcga_methylation_script", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


meth = _load_module()

P1, P2 = meth.MGMT_PROBES  # cg12434587, cg12981137

#: Filler lines that mimic a real file: ordinary cg probes, rs probes, an NA,
#: and a probe whose id merely STARTS with a real MGMT probe id.
DECOYS = [
    ("cg00000292", "0.776195891063564"),
    ("cg00002426", "0.0680390257317113"),
    ("rs966367", "NA"),
    ("cg27665659", "0.0227891643106534"),
    (P1 + "1", "0.99"),  # a longer id; must NOT be mistaken for P1
]


def _write_beta_file(path: Path, p1="0.80", p2="0.75", *, lines=None) -> Path:
    """Write a headerless probe<TAB>beta file in the GDC layout."""
    if lines is None:
        rows = list(DECOYS) + [(P1, p1), (P2, p2)]
    else:
        rows = lines
    path.write_text("".join(f"{probe}\t{value}\n" for probe, value in rows), encoding="ascii")
    return path


def _record(pid, p1, p2, platform=meth.PLATFORM_27K, file_id=None):
    values = {}
    if p1 is not None:
        values[P1] = p1
    if p2 is not None:
        values[P2] = p2
    return {"patient_id": pid, "platform": platform, "values": values,
            "n_lines": 27578, "sample_barcode": f"{pid}-01A",
            "aliquot_barcode": f"{pid}-01A-01D-01", "file_id": file_id or f"f-{pid}"}


# ---------------------------------------------------------------------------
# Beta parsing and the probe reader
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "text,expected",
    [
        ("0.5", 0.5), ("0", 0.0), ("1.0", 1.0), ("0.806195891063564", 0.806195891063564),
        (" 0.25 ", 0.25), (0.25, 0.25),
    ],
)
def test_parse_beta_accepts_numbers_in_range(text, expected):
    assert meth.parse_beta(text) == pytest.approx(expected)


@pytest.mark.parametrize(
    "text,why",
    [
        ("NA", "GDC's missing marker"), ("", "empty"), (None, "absent"),
        ("nan", "NaN is not a beta"), ("1.2", "above 1"), ("-0.1", "below 0"),
        ("0.5;0.6", "not one number"), ("inf", "infinite"),
    ],
)
def test_parse_beta_rejects_everything_else(text, why):
    """A rejected value is None -- never 0, which would read as unmethylated."""
    assert meth.parse_beta(text) is None, why


def test_read_probe_values_finds_both_probes_among_decoys(tmp_path):
    values, n_lines = meth.read_probe_values(_write_beta_file(tmp_path / "b.txt"))
    assert values == {P1: pytest.approx(0.80), P2: pytest.approx(0.75)}
    assert n_lines == len(DECOYS) + 2


def test_read_probe_values_does_not_match_a_longer_probe_id(tmp_path):
    values, _ = meth.read_probe_values(
        _write_beta_file(tmp_path / "b.txt", lines=[(P1 + "1", "0.99")])
    )
    assert values == {}


def test_read_probe_values_keeps_present_but_na_as_none_and_absent_as_absent(tmp_path):
    """Three different facts, kept different: measured, NA, not on the file."""
    path = _write_beta_file(tmp_path / "b.txt", p1="NA")
    values, _ = meth.read_probe_values(path)
    assert P1 in values and values[P1] is None
    assert values[P2] == pytest.approx(0.75)

    only_one, _ = meth.read_probe_values(
        _write_beta_file(tmp_path / "c.txt", lines=[(P1, "0.10")])
    )
    assert list(only_one) == [P1]  # P2 is absent, not None


def test_read_probe_values_keeps_the_first_occurrence_of_a_probe(tmp_path):
    lines = [(P1, "0.10"), (P2, "0.20"), (P1, "0.90")]
    values, n_lines = meth.read_probe_values(
        _write_beta_file(tmp_path / "dup.txt", lines=lines)
    )
    assert values[P1] == pytest.approx(0.10)
    assert n_lines == 3


def test_read_probe_values_on_an_empty_file(tmp_path):
    path = tmp_path / "empty.txt"
    path.write_text("", encoding="ascii")
    assert meth.read_probe_values(path) == ({}, 0)

# ---------------------------------------------------------------------------
# The score and the threshold
# ---------------------------------------------------------------------------

def test_mgmt_score_is_the_plain_mean_of_both_probes():
    assert meth.mgmt_score({P1: 0.10, P2: 0.30}) == pytest.approx(0.20)
    assert meth.mgmt_score({P1: 0.80, P2: 0.75}) == pytest.approx(0.775)


@pytest.mark.parametrize(
    "values,why",
    [
        ({P1: 0.10}, "one probe absent from the file"),
        ({P2: 0.10}, "the other probe absent"),
        ({P1: 0.10, P2: None}, "one probe present but NA"),
        ({P1: None, P2: None}, "both NA"),
        ({}, "nothing measured"),
    ],
)
def test_mgmt_score_refuses_to_average_over_fewer_probes(values, why):
    """A partial score would be silently comparable to a full one. It is not."""
    assert meth.mgmt_score(values) is None, why


def test_the_call_is_score_greater_or_equal_to_threshold():
    at = meth.build_calls([_record("TCGA-01-0001", 0.35, 0.35)], threshold=0.35)
    assert at.loc[0, "mgmt_methylated"] == 1  # exactly on the threshold
    under = meth.build_calls([_record("TCGA-01-0001", 0.30, 0.39)], threshold=0.35)
    assert under.loc[0, "mgmt_methylated"] == 0  # 0.345 mean


def test_changing_the_threshold_flips_only_the_calls_it_should():
    records = [_record("TCGA-01-0001", 0.80, 0.80),   # 0.80: methylated at any threshold
               _record("TCGA-01-0002", 0.36, 0.36),   # 0.36: methylated only below 0.40
               _record("TCGA-01-0003", 0.05, 0.05)]   # 0.05: unmethylated at any threshold
    low = meth.build_calls(records, threshold=0.30).set_index("patient_id")
    high = meth.build_calls(records, threshold=0.50).set_index("patient_id")
    assert low["mgmt_methylated"].tolist() == [1, 1, 0]
    assert high["mgmt_methylated"].tolist() == [1, 0, 0]
    # The continuous score is threshold-independent -- that is why it is the
    # primary output and the binary call a derived convenience.
    assert low["mgmt_score"].tolist() == high["mgmt_score"].tolist()


@pytest.mark.parametrize("bad", [-0.1, 1.5, 2.0])
def test_build_calls_rejects_a_threshold_outside_zero_to_one(bad):
    with pytest.raises(ValueError, match="threshold must lie in"):
        meth.build_calls([], threshold=bad)


def test_defaults_are_what_the_docstring_claims():
    assert meth.MGMT_PROBES == ("cg12434587", "cg12981137")
    assert 0.0 <= meth.DEFAULT_THRESHOLD <= 1.0
    assert meth.IS_MGMT_STP27 is False
    assert meth.CALL_METHOD == "two_probe_mean_beta"
    # Used on every parse; an undefined name here crashes the run AFTER download.
    assert meth.BETA_MIN == 0.0 and meth.BETA_MAX == 1.0


# ---------------------------------------------------------------------------
# The three-state call table
# ---------------------------------------------------------------------------

def test_not_callable_is_never_coerced_to_unmethylated():
    """The central contract: 1 / 0 / NA, with a reason for every NA."""
    records = [_record("TCGA-01-0001", 0.80, 0.75),          # methylated
               _record("TCGA-01-0002", 0.05, 0.10),          # unmethylated
               _record("TCGA-01-0003", 0.50, None)]          # one probe NA
    calls = meth.build_calls(records, ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003",
                                       "TCGA-01-0004"])
    row = calls.set_index("patient_id")

    assert row.loc["TCGA-01-0001", "mgmt_methylated"] == 1
    assert row.loc["TCGA-01-0001", "call_basis"] == meth.BASIS_CALLED
    assert row.loc["TCGA-01-0001", "mgmt_score"] == pytest.approx(0.775)

    assert row.loc["TCGA-01-0002", "mgmt_methylated"] == 0
    assert row.loc["TCGA-01-0002", "call_basis"] == meth.BASIS_CALLED

    # A file exists, but a required probe is missing: not callable, with the
    # measured probe kept so the loss is visible.
    assert pd.isna(row.loc["TCGA-01-0003", "mgmt_methylated"])
    assert row.loc["TCGA-01-0003", "call_basis"] == meth.BASIS_PROBE_UNAVAILABLE
    assert row.loc["TCGA-01-0003", "cg12434587"] == pytest.approx(0.50)
    assert pd.isna(row.loc["TCGA-01-0003", "mgmt_score"])

    # No file at all.
    assert pd.isna(row.loc["TCGA-01-0004", "mgmt_methylated"])
    assert row.loc["TCGA-01-0004", "call_basis"] == meth.BASIS_MISSING
    assert row.loc["TCGA-01-0004", "call_method"] is None or pd.isna(
        row.loc["TCGA-01-0004", "call_method"])

    assert list(calls.columns) == meth.CALL_COLUMNS
    assert str(calls["mgmt_methylated"].dtype) == "Int64"
    assert list(calls["patient_id"]) == sorted(calls["patient_id"])


def test_every_call_carries_its_method_and_threshold():
    calls = meth.build_calls([_record("TCGA-01-0001", 0.9, 0.9)], threshold=0.42)
    row = calls.iloc[0]
    assert row["call_method"] == "two_probe_mean_beta"
    assert float(row["threshold"]) == pytest.approx(0.42)
    # Nothing in the output may read as the published MGMT-STP27 model.
    assert "stp27" not in str(row["call_method"]).lower()

# ---------------------------------------------------------------------------
# Platform preference and duplicate resolution
# ---------------------------------------------------------------------------

def test_a_patient_profiled_on_both_platforms_gets_the_450k_call():
    records = [_record("TCGA-01-0001", 0.10, 0.10, platform=meth.PLATFORM_27K, file_id="f27"),
               _record("TCGA-01-0001", 0.90, 0.90, platform=meth.PLATFORM_450K, file_id="f450")]
    calls = meth.build_calls(records)
    assert len(calls) == 1
    assert calls.loc[0, "mgmt_score"] == pytest.approx(0.90)
    assert calls.loc[0, "platform"] == "450k"
    assert calls.loc[0, "file_id"] == "f450"


def test_build_calls_is_order_independent():
    records = [_record("TCGA-02-0002", 0.20, 0.20, platform=meth.PLATFORM_450K),
               _record("TCGA-01-0001", 0.70, 0.70, platform=meth.PLATFORM_27K),
               _record("TCGA-02-0002", 0.10, 0.10, platform=meth.PLATFORM_27K)]
    first = meth.build_calls(records, ["TCGA-03-0003"])
    second = meth.build_calls(records[::-1], ["TCGA-03-0003"])
    pd.testing.assert_frame_equal(first, second)
    assert first.loc[0, "platform"] == "27k"   # TCGA-01-0001, 27k only
    assert first.loc[1, "platform"] == "450k"  # TCGA-02-0002, both -> 450k


def test_calls_survive_a_parquet_round_trip_without_becoming_zero(tmp_path):
    calls = meth.build_calls([_record("TCGA-01-0001", 0.05, 0.05)],
                             ["TCGA-01-0001", "TCGA-01-0009"])
    path = tmp_path / "mgmt_calls.parquet"
    calls.to_parquet(path, index=False)
    back = pd.read_parquet(path).set_index("patient_id")
    assert str(back["mgmt_methylated"].dtype) == "Int64"
    assert back.loc["TCGA-01-0001", "mgmt_methylated"] == 0
    assert pd.isna(back.loc["TCGA-01-0009", "mgmt_methylated"])
    assert back.loc["TCGA-01-0009", "call_basis"] == meth.BASIS_MISSING


def test_build_calls_with_nothing_at_all_returns_the_contract_columns():
    calls = meth.build_calls([])
    assert calls.empty and list(calls.columns) == meth.CALL_COLUMNS


def _meth_manifest():
    return pd.DataFrame([
        {"file_id": "f-a27", "patient_id": "TCGA-A", "sample_type_code": "01",
         "platform": meth.PLATFORM_27K, "aliquot_barcode": "TCGA-A-01A-01D-01"},
        {"file_id": "f-a450", "patient_id": "TCGA-A", "sample_type_code": "01",
         "platform": meth.PLATFORM_450K, "aliquot_barcode": "TCGA-A-01A-01D-02"},
        {"file_id": "f-b27", "patient_id": "TCGA-B", "sample_type_code": "01",
         "platform": meth.PLATFORM_27K, "aliquot_barcode": "TCGA-B-01A-01D-01"},
        {"file_id": "f-b27b", "patient_id": "TCGA-B", "sample_type_code": "01",
         "platform": meth.PLATFORM_27K, "aliquot_barcode": "TCGA-B-01B-01D-01"},
        {"file_id": "f-c27", "patient_id": "TCGA-C", "sample_type_code": "02",
         "platform": meth.PLATFORM_27K, "aliquot_barcode": "TCGA-C-02A-01D-01"},
        {"file_id": "f-c450", "patient_id": "TCGA-C", "sample_type_code": "01",
         "platform": meth.PLATFORM_450K, "aliquot_barcode": "TCGA-C-01A-01D-01"},
    ])


def test_resolve_duplicate_files_prefers_primary_then_450k():
    deduped, report = meth.resolve_duplicate_files(_meth_manifest())
    picked = deduped.set_index("patient_id")["file_id"].to_dict()
    assert picked == {"TCGA-A": "f-a450", "TCGA-B": "f-b27", "TCGA-C": "f-c450"}
    assert report["patients"] == 3
    assert report["patients_with_duplicates"] == 3
    # Counted over the manifest AS PASSED IN: TCGA-C has a 27k recurrent and a
    # 450k primary file, so it is on both platforms even though only the 450k
    # file is kept. (main() filters to primary tumours before calling this.)
    assert report["patients_on_both_platforms"] == 2
    assert report["files_dropped_as_duplicate"] == 3
    assert "_priority" not in deduped.columns


def test_resolve_duplicate_files_is_input_order_independent():
    first, _ = meth.resolve_duplicate_files(_meth_manifest())
    second, _ = meth.resolve_duplicate_files(_meth_manifest().iloc[::-1])
    pd.testing.assert_frame_equal(first, second)


def test_resolve_duplicate_files_on_empty_manifest():
    deduped, report = meth.resolve_duplicate_files(_meth_manifest().iloc[0:0])
    assert deduped.empty and report["patients"] == 0


def test_build_file_filters_can_restrict_to_one_platform():
    every = meth.build_file_filters("TCGA-GBM")
    one = meth.build_file_filters("TCGA-GBM", meth.PLATFORM_450K)
    assert len(every["content"]) == 4
    assert len(one["content"]) == 5
    assert one["content"][-1]["content"]["value"] == [meth.PLATFORM_450K]


def test_flatten_file_records_platform_and_tumour_aliquot():
    hit = {
        "file_id": "f1", "file_name": "x.txt", "file_size": 13_000_000, "md5sum": "abc",
        "platform": meth.PLATFORM_450K,
        "cases": [{
            "case_id": "c1", "submitter_id": "TCGA-76-4928",
            "samples": [
                {"submitter_id": "TCGA-76-4928-10A", "sample_type": "Blood Derived Normal",
                 "portions": [{"analytes": [{"aliquots": [
                     {"submitter_id": "TCGA-76-4928-10A-01D-01"}]}]}]},
                {"submitter_id": "TCGA-76-4928-01B", "sample_type": "Primary Tumor",
                 "portions": [{"analytes": [{"aliquots": [
                     {"submitter_id": "TCGA-76-4928-01B-01D-01"}]}]}]},
            ],
        }],
    }
    row = meth.flatten_file_record(hit)
    assert row["platform"] == meth.PLATFORM_450K
    assert row["patient_id"] == "TCGA-76-4928"
    assert row["sample_type_code"] == "01"
    assert row["sample_barcode"] == "TCGA-76-4928-01B"

# ---------------------------------------------------------------------------
# Summaries, probe coverage and the printed report
# ---------------------------------------------------------------------------

def _mixed_calls():
    records = [_record("TCGA-01-0001", 0.90, 0.90, platform=meth.PLATFORM_450K),
               _record("TCGA-01-0002", 0.36, 0.36, platform=meth.PLATFORM_450K),
               _record("TCGA-01-0003", 0.05, 0.05, platform=meth.PLATFORM_27K),
               _record("TCGA-01-0004", 0.50, None, platform=meth.PLATFORM_27K)]
    return meth.build_calls(records, ["TCGA-01-0005"]), records


def test_summarize_calls_counts_the_three_states_separately():
    summary = meth.summarize_calls(_mixed_calls()[0])
    assert summary["patients_in_table"] == 5
    assert summary["methylated"] == 2
    assert summary["unmethylated"] == 1
    assert summary["not_callable"] == 2  # one NA probe + one with no file
    assert summary["callable"] == 3
    # The denominator is CALLABLE patients; the not-callable never dilute it.
    assert summary["methylated_fraction_of_callable"] == pytest.approx(2 / 3)
    assert summary["is_mgmt_stp27"] is False
    assert summary["call_method"] == "two_probe_mean_beta"


def test_summarize_breaks_calls_down_by_platform():
    summary = meth.summarize_calls(_mixed_calls()[0])
    assert summary["by_platform"] == {
        "450k": {"callable": 2, "methylated": 2},
        "27k": {"callable": 1, "methylated": 0},
    }


def test_summarize_reports_threshold_sensitivity_and_a_histogram():
    summary = meth.summarize_calls(_mixed_calls()[0])
    sensitivity = {row["threshold"]: row for row in summary["threshold_sensitivity"]}
    assert set(sensitivity) == set(meth.SENSITIVITY_THRESHOLDS)
    # Callable scores are 0.90, 0.36 and 0.05: the borderline sample flips
    # between 0.35 and 0.40, which is precisely what this table is for.
    assert sensitivity[0.20]["methylated"] == 2
    assert sensitivity[0.35]["methylated"] == 2
    assert sensitivity[0.40]["methylated"] == 1
    assert sensitivity[0.60]["methylated"] == 1
    assert sensitivity[0.35]["fraction_of_callable"] == pytest.approx(2 / 3)
    # The histogram covers every callable patient exactly once.
    assert sum(row["n"] for row in summary["score_histogram"]) == summary["callable"]


def test_summarize_on_an_empty_table_has_no_division_error():
    summary = meth.summarize_calls(meth.build_calls([]))
    assert summary["methylated"] == 0
    assert summary["methylated_fraction_of_callable"] is None
    assert summary["score_histogram"] == []


def test_probe_coverage_distinguishes_measured_na_and_absent():
    records = [
        {"platform": meth.PLATFORM_27K, "values": {P1: 0.4, P2: 0.5}},
        {"platform": meth.PLATFORM_27K, "values": {P1: None, P2: 0.1}},   # P1 is NA
        {"platform": meth.PLATFORM_450K, "values": {P1: 0.2}},            # P2 not on file
    ]
    cov = meth.probe_coverage(records)
    assert cov["27k"]["files"] == 2 and cov["450k"]["files"] == 1
    assert cov["27k"]["probes"][P1] == {"measured": 1, "na": 1, "absent": 0}
    assert cov["450k"]["probes"][P2] == {"measured": 0, "na": 0, "absent": 1}


def test_coverage_counts_how_many_of_a_patient_set_gain_a_call():
    calls = meth.build_calls(
        [_record("TCGA-01-0001", 0.9, 0.9), _record("TCGA-01-0002", 0.0, 0.0)],
        ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003"],
    )
    row = meth.coverage(calls, ["TCGA-01-0001", "TCGA-01-0003", "TCGA-09-0009"], "expression")
    assert row == {"set": "expression", "patients": 3, "with_mgmt_call": 1,
                   "methylated": 1, "still_missing": 2}


def test_printed_report_is_ascii_and_says_it_is_not_stp27(capsys):
    calls, records = _mixed_calls()
    meth.print_report(meth.summarize_calls(calls),
                      [meth.coverage(calls, ["TCGA-01-0001"], "clinical cohort")],
                      meth.probe_coverage(records))
    out = capsys.readouterr().out
    assert out.isascii()
    assert "NOT MGMT-STP27" in out
    assert "not callable" in out and "NOT unmethylated" in out
    assert "<- default" in out
    assert "cg12434587" in out and "cg12981137" in out


def test_script_source_is_ascii_only():
    """cp1252 consoles raise UnicodeEncodeError on non-ASCII stdout."""
    assert SCRIPT_PATH.read_bytes().decode("utf-8").isascii()

# ---------------------------------------------------------------------------
# Download layer (tar extraction; no network) and end to end
# ---------------------------------------------------------------------------

def _tar_bytes(payloads: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        info = tarfile.TarInfo("MANIFEST.txt")
        info.size = 4
        tar.addfile(info, io.BytesIO(b"meta"))
        for file_id, blob in payloads.items():
            info = tarfile.TarInfo(f"{file_id}/{file_id}.methylation_array.level3betas.txt")
            info.size = len(blob)
            tar.addfile(info, io.BytesIO(blob))
    return buf.getvalue()


def test_download_batch_unpacks_a_multi_file_tar(tmp_path, monkeypatch):
    payloads = {"id-1": b"beta-1", "id-2": b"beta-2"}

    class _Response:
        status_code = 200
        text = ""
        content = _tar_bytes(payloads)

    monkeypatch.setattr(meth.requests, "post", lambda *a, **k: _Response())
    cache = tmp_path / "_cache"
    assert meth.download_batch(list(payloads), cache) == 2
    assert (cache / "id-1.txt").read_bytes() == b"beta-1"
    assert (cache / "id-2.txt").read_bytes() == b"beta-2"
    assert not (cache / "MANIFEST.txt").exists()


def test_download_batch_retries_a_transient_500_then_succeeds(tmp_path, monkeypatch):
    class _Response:
        def __init__(self, status, content=b""):
            self.status_code = status
            self.text = ""
            self.content = content

    answers = iter([_Response(500), _Response(500), _Response(200, b"ok")])
    monkeypatch.setattr(meth.requests, "post", lambda *a, **k: next(answers))
    monkeypatch.setattr(meth.time, "sleep", lambda seconds: None)
    cache = tmp_path / "_cache"
    assert meth.download_batch(["x"], cache, max_attempts=3, retry_wait=5) == 1
    assert (cache / "x.txt").read_bytes() == b"ok"


def test_download_files_isolates_a_bad_file_and_keeps_going(tmp_path, monkeypatch):
    # Like the real download_batch, this fake is all-or-nothing per batch: the
    # write only happens after a good response, so a failed batch writes nothing.
    def fake_batch(file_ids, cache_dir, timeout=1800, max_attempts=3, retry_wait=2.0):
        if "bad" in file_ids:
            raise RuntimeError("GDC data error 500: internal server error")
        for file_id in file_ids:
            meth.cache_path_for(cache_dir, file_id).write_bytes(b"z")
        return len(file_ids)

    monkeypatch.setattr(meth, "download_batch", fake_batch)
    report = meth.download_files(["good1", "bad", "good2"], tmp_path, batch_size=3)
    # The good files are fetched on the one-at-a-time retry and counted ONCE.
    assert report["downloaded"] == 2 and report["failed"] == ["bad"]
    assert report["requested"] == 3 and report["skipped_cached"] == 0


def _gdc_hit(file_id, patient, platform):
    return {
        "file_id": file_id, "file_name": f"{file_id}.txt", "file_size": 13_000_000,
        "md5sum": "abc", "platform": platform,
        "cases": [{
            "case_id": f"case-{patient}", "submitter_id": patient,
            "samples": [{
                "submitter_id": f"{patient}-01A", "sample_type": "Primary Tumor",
                "portions": [{"analytes": [{"aliquots": [{
                    "submitter_id": f"{patient}-01A-01D-01"}]}]}],
            }],
        }],
    }

def test_main_runs_end_to_end_offline(tmp_path, monkeypatch, capsys):
    out_dir = tmp_path / "methylation"

    def fake_meta(project_id=meth.PROJECT_ID, page_size=500):
        return [_gdc_hit("f1", "TCGA-01-0001", meth.PLATFORM_450K),
                _gdc_hit("f2", "TCGA-01-0002", meth.PLATFORM_27K),
                _gdc_hit("f3", "TCGA-01-0002", meth.PLATFORM_450K)]  # both platforms

    def fake_download(file_ids, cache_dir, batch_size=4):
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        _write_beta_file(cache_dir / "f1.txt", "0.80", "0.80")
        _write_beta_file(cache_dir / "f2.txt", "0.05", "0.05")
        _write_beta_file(cache_dir / "f3.txt", "0.90", "0.90")
        return {"requested": len(file_ids), "skipped_cached": 0, "downloaded": len(file_ids),
                "failed": [], "elapsed_seconds": 0.0}

    monkeypatch.setattr(meth, "fetch_file_metadata", fake_meta)
    monkeypatch.setattr(meth, "download_files", fake_download)

    cohort = tmp_path / "clinical_cohort.parquet"
    pd.DataFrame({"patient_id": ["TCGA-01-0001", "TCGA-01-0002", "TCGA-01-0003"]}).to_parquet(cohort)
    # Laid out like the real matrix: gene x sample written with index=False, so
    # its first column is gene_id and NOT a patient (it was once counted as one).
    expression = tmp_path / "tpm.parquet"
    pd.DataFrame({"gene_id": ["ENSG00000141510"], "TCGA-01-0001": [1.0],
                  "TCGA-01-0003": [3.0]}).to_parquet(expression, index=False)
    monkeypatch.setattr(sys, "argv", [
        "download_tcga_methylation.py", "--out", str(out_dir), "--cohort", str(cohort),
        "--expression", str(expression),
    ])
    meth.main()

    calls = pd.read_parquet(out_dir / "mgmt_calls.parquet").set_index("patient_id")
    assert calls.loc["TCGA-01-0001", "mgmt_methylated"] == 1
    assert calls.loc["TCGA-01-0001", "platform"] == "450k"
    # Profiled on both platforms: the 450k file (0.90) won over the 27k one (0.05).
    assert calls.loc["TCGA-01-0002", "mgmt_score"] == pytest.approx(0.90)
    assert calls.loc["TCGA-01-0002", "platform"] == "450k"
    assert pd.isna(calls.loc["TCGA-01-0003", "mgmt_methylated"])
    assert calls.loc["TCGA-01-0003", "call_basis"] == meth.BASIS_MISSING

    report = json.loads((out_dir / "mgmt_call_report.json").read_text(encoding="ascii"))
    assert report["is_mgmt_stp27"] is False
    assert report["duplicate_resolution"]["patients_on_both_platforms"] == 1
    assert report["probe_coverage"]["450k"]["probes"][P1]["measured"] == 2
    # gene_id is not a patient: the expression set is 2, not 3.
    assert report["expression_non_sample_columns"] == ["gene_id"]
    assert [(row["set"], row["patients"], row["with_mgmt_call"]) for row in report["coverage"]] == [
        ("clinical cohort", 3, 2), ("with expression", 2, 1), ("cohort AND expression", 2, 1)]

    printed = capsys.readouterr().out
    assert printed.isascii()
    assert "NOT MGMT-STP27" in printed
    assert len(pd.read_parquet(cohort)) == 3  # the cohort file is never rewritten


def test_expression_patient_ids_never_count_the_gene_column(tmp_path):
    """Regression: the real tpm.parquet has a `gene_id` column, and it was
    counted as a patient, overstating the expression set by one."""
    path = tmp_path / "tpm.parquet"
    pd.DataFrame({"gene_id": ["EGFR", "IDH1"], "TCGA-01-0001": [1.0, 2.0],
                  "TCGA-01-0003": [3.0, 4.0]}).to_parquet(path, index=False)
    ids, ignored = meth.expression_patient_ids(path)
    assert ids == ["TCGA-01-0001", "TCGA-01-0003"]
    assert ignored == ["gene_id"]

    # A file with no barcode column at all is the wrong file: refuse it.
    other = tmp_path / "clinical_cohort.parquet"
    pd.DataFrame({"patient_id": ["TCGA-01-0001"], "age": [61]}).to_parquet(other, index=False)
    with pytest.raises(RuntimeError, match="not a gene x sample matrix"):
        meth.expression_patient_ids(other)


def test_main_dry_run_downloads_and_writes_nothing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(meth, "fetch_file_metadata", lambda project_id=meth.PROJECT_ID,
                        page_size=500: [_gdc_hit("f1", "TCGA-01-0001", meth.PLATFORM_27K)])

    def must_not_download(*args, **kwargs):
        raise AssertionError("dry-run must not download")

    monkeypatch.setattr(meth, "download_files", must_not_download)
    out_dir = tmp_path / "dry"
    monkeypatch.setattr(sys, "argv", [
        "download_tcga_methylation.py", "--out", str(out_dir), "--dry-run",
    ])
    meth.main()
    assert "--dry-run -> nothing downloaded" in capsys.readouterr().out
    assert not out_dir.exists()
