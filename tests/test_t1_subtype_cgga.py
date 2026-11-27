"""Synthetic-only tests for scripts/t1_subtype_cgga.py. No CGGA file is read."""

from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import binomtest

from gliostrata.profiling.centroid import load_centroid_set

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("t1_subtype_cgga", ROOT / "scripts" / "t1_subtype_cgga.py")
t1 = importlib.util.module_from_spec(_spec)
sys.modules["t1_subtype_cgga"] = t1
_spec.loader.exec_module(t1)


@pytest.fixture(scope="module")
def cs():
    return load_centroid_set(t1.CENTROID_FILE)


def _zip(path: Path, frame: pd.DataFrame) -> Path:
    text = "Gene_Name\t" + "\t".join(frame.columns) + "\n" + "".join(
        g + "\t" + "\t".join(f"{v:.6f}" for v in row) + "\n" for g, row in zip(frame.index, frame.to_numpy()))
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("m.txt", text)
    return path


def test_duplicate_symbols_are_averaged(tmp_path):
    frame = pd.DataFrame([[1.0, 4.0], [3.0, 8.0], [5.0, 5.0]], index=["A", "A", "B"], columns=["s1", "s2"])
    out = t1.read_matrix(_zip(tmp_path / "m.zip", frame))
    assert list(out.index) == ["A", "B"]
    assert out.loc["A"].tolist() == [2.0, 6.0]
    assert out.loc["B"].tolist() == [5.0, 5.0]


def test_exact_symbol_match_no_alias_or_case_folding(cs):
    genes = list(cs.genes)
    rng = np.random.default_rng(0)
    frame = pd.DataFrame(rng.normal(8, 1, (len(genes), 2)), index=genes, columns=["a", "b"])
    full = t1.classify_log2(frame, cs)
    assert full["n_genes_used"].iloc[0] == len(genes)
    renamed = frame.copy()
    renamed.index = [genes[0].lower()] + genes[1:]  # case differs: must not match
    renamed.index = [i + "-AS1" if i == genes[1] else i for i in renamed.index]  # alias-like: no match
    part = t1.classify_log2(renamed, cs)
    assert part["n_genes_used"].iloc[0] == len(genes) - 2


def test_both_input_scales(cs):
    out, flag = t1.to_log2(pd.DataFrame({"s": [0.0, 3.0]}), "log2_plus_one")
    assert out["s"].tolist() == [0.0, 2.0] and flag is False
    raw = pd.DataFrame({"s": [0.5, 3.0]})
    out, flag = t1.to_log2(raw, "log2_as_distributed")
    assert out["s"].tolist() == [0.5, 3.0] and flag is False
    with pytest.raises(ValueError):
        t1.to_log2(raw, "natural_log")


def test_negative_values_trigger_the_registered_scale_contingency():
    # Linear intensities cannot be negative: such values are already log-scale.
    neg = pd.DataFrame({"s": [-1.5, 0.2, 4.0]})
    out, flag = t1.to_log2(neg, "log2_plus_one")
    assert flag is True and out["s"].tolist() == [-1.5, 0.2, 4.0]


def test_missing_value_drops_that_gene_for_that_sample_only(cs):
    genes = list(cs.genes)
    rng = np.random.default_rng(3)
    expr = pd.DataFrame(rng.normal(8, 2, size=(len(genes), 2)), index=genes, columns=["a", "b"])
    expr.iloc[0, 0] = np.nan
    out = t1.classify_log2(expr, cs).set_index("patient_id")
    assert out.loc["a", "n_genes_used"] == len(genes) - 1
    assert out.loc["b", "n_genes_used"] == len(genes)
    assert out.loc["a", "subtype"] is not None


def test_scales_give_same_call_on_matching_data(cs, tmp_path):
    for scale in t1.SCALES:
        paths = t1.write_synthetic(tmp_path / scale, scale)
        res = t1.run_evaluation(paths, scale)
        assert res["input_scale"] == scale
        assert res["accuracy"] > 0.9  # synthetic profiles are built from the centroids
        assert res["exclusions"]["n_outside_manifest_in_clinical_table"] == 5
        assert res["exclusions"]["n_without_reference_label"] == 6
        assert res["n_genes_used_min"] == res["n_genes_used_max"] == len(cs.genes)
        assert res["scale_contingency_applied"] is False


@pytest.mark.parametrize(("k", "n"), [(102, 108), (0, 40), (40, 40), (28, 31), (95, 100)])
def test_clopper_pearson_matches_scipy(k, n):
    ci = binomtest(k, n).proportion_ci(confidence_level=0.95, method="exact")
    lo, hi = t1.clopper_pearson(k, n)
    assert lo == pytest.approx(ci.low)
    assert hi == pytest.approx(ci.high)


def test_floor_of_30_is_not_evaluable():
    labels = pd.Series(["Classical"] * 29)
    res = t1.score_calls(labels.copy(), labels, ("Classical", "Mesenchymal", "Proneural", "Neural"))
    assert res["verdict"]["status"] == "NOT EVALUABLE"
    ok = pd.Series(["Classical"] * 30)
    res30 = t1.score_calls(ok.copy(), ok, ("Classical", "Mesenchymal", "Proneural", "Neural"))
    assert res30["verdict"]["status"] == "PASS"
    assert "negative_control" in res30


def test_verdict_uses_lower_bound():
    sub = ("Classical", "Mesenchymal", "Proneural", "Neural")
    labels = pd.Series(["Classical"] * 108)
    calls = labels.copy()
    calls.iloc[:6] = "Neural"  # 102 / 108
    assert t1.score_calls(calls, labels, sub)["verdict"]["status"] == "PASS"
    calls.iloc[:7] = "Neural"  # 101 / 108
    assert t1.score_calls(calls, labels, sub)["verdict"]["status"] == "INCONCLUSIVE"


def test_refuses_without_scale_config(tmp_path):
    with pytest.raises(SystemExit):
        t1.load_scale(tmp_path / "absent.json")
    with pytest.raises(SystemExit):
        t1.main(["dry-run", "--scale-config", str(tmp_path / "absent.json"),
                 "--out", str(tmp_path / "o.json")])
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"scale": "ln"}), encoding="utf-8")
    with pytest.raises(SystemExit):
        t1.load_scale(bad)


def test_evaluate_refuses_without_dry_run(tmp_path):
    cfg = tmp_path / "scale.json"
    cfg.write_text(json.dumps({"scale": "log2_plus_one"}), encoding="utf-8")
    with pytest.raises(SystemExit, match="dry run"):
        t1.main(["evaluate", "--matrix", "m", "--clinical", "c", "--manifest", "j",
                 "--scale-config", str(cfg), "--dry-run-out", str(tmp_path / "none.json"),
                 "--ledger", str(tmp_path / "ledger.jsonl")])
    assert not (tmp_path / "ledger.jsonl").exists()


def test_dry_run_writes_output(tmp_path):
    cfg = tmp_path / "scale.json"
    cfg.write_text(json.dumps({"scale": "log2_plus_one"}), encoding="utf-8")
    out = tmp_path / "dry.json"
    t1.main(["dry-run", "--scale-config", str(cfg), "--out", str(out)])
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert doc["synthetic"] is True
    assert doc["negative_control"]["n_perm"] == 1000
    assert 0.1 < doc["negative_control"]["mean"] < 0.4


def test_per_sample_independence(cs, tmp_path):
    paths = t1.write_synthetic(tmp_path, "log2_as_distributed")
    expr, _ = t1.to_log2(t1.read_matrix(paths["matrix"]), "log2_as_distributed")
    together = t1.classify_log2(expr, cs).set_index("patient_id")["subtype"]
    for sid in expr.columns[:8]:
        alone = t1.classify_log2(expr[[sid]], cs)["subtype"].iloc[0]
        assert alone == together[sid]
