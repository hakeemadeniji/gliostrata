"""Checks data/provenance/file_manifest_checksums.json against the files it describes.

The per-modality GDC file manifests are not committed: each one maps every
downloaded file to a patient barcode. What is committed is their size, row
count and SHA-256. A manifest that is present on disk must match all three. A
missing one is skipped, as it is in a clean checkout.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "data" / "provenance" / "file_manifest_checksums.json"
RAW = ROOT / "data" / "raw" / "tcga_gbm"
EXPECTED = {
    "expression/file_manifest.parquet",
    "methylation/file_manifest.parquet",
    "mutations/file_manifest.parquet",
}


def _record() -> dict:
    return json.loads(RECORD.read_text(encoding="utf-8"))


def test_the_record_names_each_manifest_and_why_it_is_not_committed():
    record = _record()
    assert set(record["manifests"]) == EXPECTED
    assert "patient" in record["what"].lower() and "not committed" in record["what"]
    for column in ("patient_id", "case_submitter_id", "aliquot_barcode"):
        assert column in record["columns_in_each"]
    for entry in record["manifests"].values():
        assert len(entry["sha256"]) == 64
        assert entry["bytes"] > 0 and entry["rows"] > 0
    assert set(record["sources"]) == {name.split("/")[0] for name in EXPECTED}


def test_no_file_manifest_is_tracked_by_git():
    """They carry patient barcodes; only their checksums travel with git."""
    out = subprocess.run(
        ["git", "ls-files", "data/raw"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout
    assert "file_manifest" not in out
    assert (ROOT / "data" / "provenance" / "file_manifest_checksums.json").is_file()


def test_each_manifest_on_disk_matches_size_hash_and_row_count():
    checked = 0
    for relative, entry in _record()["manifests"].items():
        path = RAW / relative
        if not path.is_file():
            continue
        data = path.read_bytes()
        assert len(data) == entry["bytes"], f"{relative}: size changed"
        assert hashlib.sha256(data).hexdigest() == entry["sha256"], (
            f"{relative}: content changed since data/provenance was written. Re-create "
            "it with the ingest script and update the record in the same commit."
        )
        assert len(pd.read_parquet(path)) == entry["rows"], relative
        checked += 1
    if not checked:
        pytest.skip("no GDC file manifest on disk (clean checkout)")
