"""T1 state-deconvolution row: scripts/t1_deconvolution_mixtures.py.

What is at stake: the registered T1 deconvolution analysis is run once, so every
rule of its specification (docs/t1_deconvolution_spec.md) must be right before
the run. These tests drive the script's functions on synthetic atlases; the one
real-data test reads only the atlas's cell annotation and skips when it is absent.
No cohort is read.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "t1_deconvolution_mixtures.py"


def _load():
    spec = importlib.util.spec_from_file_location("t1_deconvolution_mixtures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


t1 = _load()


def _cells(rows):
    cols = ["cell_name", "sample", "malignant", *t1.MODULE_SCORES]
    return pd.DataFrame(rows, columns=cols)


# --- 1. state assignment ------------------------------------------------------------
def test_state_is_the_highest_of_four_with_mes_and_npc_merged():
    # module order: MESlike2, MESlike1, AClike, OPClike, NPClike1, NPClike2
    cells = _cells([
        ("c1", "T1", "yes", 0.1, 2.0, 0.5, 0.0, 0.0, 0.0),    # MES via MESlike1
        ("c2", "T1", "yes", 0.0, 0.0, 0.5, 0.0, 0.2, 0.9),    # NPC via NPClike2
        ("c3", "T1", "yes", 0.0, 0.0, 0.1, 0.7, 0.0, 0.0),    # OPC
        ("c4", "T2", "yes", 0.0, 0.0, 1.5, 0.7, 0.0, 0.0),    # AC
        ("c5", "T2", "no", 0.0, 0.0, 3.0, 0.0, 0.0, 0.0),     # not malignant
        ("c6", "T2", "yes", None, 0.0, 3.0, 0.0, 0.0, 0.0),   # a score missing
        ("c7", "T2", None, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0),     # malignant unknown
        ("c8", "T2", "yes", 0.0, 0.0, 1.0, 1.0, 0.0, 0.0),    # AC/OPC tie
    ])
    out = t1.assign_states(cells)
    assert dict(zip(out["cell_name"], out["state"], strict=True)) == {
        "c1": "MES", "c2": "NPC", "c3": "OPC", "c4": "AC"}
    assert out.attrs["n_annotated"] == 8
    assert out.attrs["n_malignant_with_scores"] == 5
    assert out.attrs["n_tied"] == 1


def test_assign_states_refuses_a_table_without_the_module_scores():
    with pytest.raises(ValueError, match="lacks column"):
        t1.assign_states(pd.DataFrame({"cell_name": ["a"], "sample": ["T"], "malignant": ["yes"]}))


# --- 2. weights and mixtures ------------------------------------------------------------
def test_weights_are_a_flat_dirichlet_draw_and_reproducible():
    a = t1.draw_weights(np.random.default_rng(t1.SEED_MIXTURES))
    b = t1.draw_weights(np.random.default_rng(t1.SEED_MIXTURES))
    assert a.shape == (100, 4)
    assert np.allclose(a.sum(axis=1), 1.0) and (a > 0).all()
    assert np.array_equal(a, b)
    assert t1.DIRICHLET_ALPHA == (1.0, 1.0, 1.0, 1.0)


def _toy_atlas(n_per_state=30, n_genes=6, seed=0):
    rng = np.random.default_rng(seed)
    names, states, cols = [], [], []
    for k, s in enumerate(t1.STATES):
        for i in range(n_per_state):
            v = rng.gamma(2.0, 1.0, n_genes)
            v[k] += 50.0
            cols.append(v)
            names.append(f"{s}_{i}")
            states.append(s)
    tpm = pd.DataFrame(np.column_stack(cols), index=[f"G{i}" for i in range(n_genes)],
                       columns=names)
    cell_states = pd.DataFrame({"cell_name": names, "sample": "T", "state": states})
    return tpm, cell_states


def test_a_mixture_is_the_mean_of_its_drawn_cells_and_the_draw_is_reproducible():
    tpm, cs = _toy_atlas()
    w = np.array([[0.7, 0.1, 0.1, 0.1], [0.25, 0.25, 0.25, 0.25]])
    m1, r1 = t1.make_mixtures(tpm, cs, w, np.random.default_rng(5), cells_per_mixture=20)
    m2, r2 = t1.make_mixtures(tpm, cs, w, np.random.default_rng(5), cells_per_mixture=20)
    assert m1.equals(m2) and np.array_equal(r1, r2)
    assert list(m1.columns) == ["mix001", "mix002"]
    assert np.allclose(r1.sum(axis=1), 1.0)
    # Re-derive mixture 1 by replaying the generator.
    rng = np.random.default_rng(5)
    counts = rng.multinomial(20, w[0])
    chosen = []
    for s, k in zip(t1.STATES, counts, strict=True):
        pool = cs.loc[cs["state"] == s, "cell_name"].to_numpy()
        chosen += list(rng.choice(pool, size=int(k), replace=False)) if k else []
    assert np.allclose(m1["mix001"].to_numpy(), tpm[chosen].mean(axis=1).to_numpy())
    assert np.allclose(r1[0], counts / 20)


def test_a_state_pool_too_small_for_the_draw_raises():
    tpm, cs = _toy_atlas(n_per_state=3)
    with pytest.raises(ValueError, match="in the pool"):
        t1.make_mixtures(tpm, cs, np.array([[1.0, 0.0, 0.0, 0.0]]),
                         np.random.default_rng(0), cells_per_mixture=10)


# --- 3. scoring --------------------------------------------------------------------------
def test_per_state_r_is_pearson_per_column_and_nan_when_undefined():
    rng = np.random.default_rng(1)
    truth = rng.dirichlet(np.ones(4), size=50)
    assert np.allclose(t1.per_state_r(truth, truth), 1.0)
    est = truth.copy()
    est[:, 2] = 0.25
    r = t1.per_state_r(est, truth)
    assert np.isnan(r[2]) and np.allclose(r[[0, 1, 3]], 1.0)
    noisy = truth + rng.normal(0, 0.05, truth.shape)
    expected = [np.corrcoef(noisy[:, j], truth[:, j])[0, 1] for j in range(4)]
    assert np.allclose(t1.per_state_r(noisy, truth), expected)


def test_bootstrap_resamples_mixtures_and_is_reproducible():
    rng = np.random.default_rng(2)
    truth = rng.dirichlet(np.ones(4), size=100)
    est = truth + rng.normal(0, 0.05, truth.shape)
    a = t1.bootstrap(est, truth, np.random.default_rng(42), n_boot=200)
    b = t1.bootstrap(est, truth, np.random.default_rng(42), n_boot=200)
    assert a == b
    point = float(t1.per_state_r(est, truth).mean())
    assert a["mean_r"]["lower"] <= point <= a["mean_r"]["upper"]
    assert a["n_undefined_dropped"] == 0
    perfect = t1.bootstrap(truth, truth, np.random.default_rng(42), n_boot=50)
    assert perfect["mean_r"]["lower"] == pytest.approx(1.0)


def test_verdict_comes_from_the_registry():
    assert t1.registry_verdict(0.9, 0.8, 0.95)["status"] == "PASS"
    assert t1.registry_verdict(0.8, 0.7, 0.85)["status"] == "INCONCLUSIVE"
    assert t1.registry_verdict(0.6, 0.5, 0.7)["status"] == "FAIL"
    v = t1.registry_verdict(0.76, 0.77, 0.9)  # point below its percentile interval
    assert v["interval_widened_to_include_point"] is True
    assert v["ci_lower"] == pytest.approx(0.76) and v["status"] == "PASS"
    assert t1.registry_verdict(0.9, 0.8, 0.95)["metric_key"] == t1.METRIC_KEY


# --- 4. the run's guard rails -------------------------------------------------------------
def _write_synthetic_inputs(tmp_path, n_per_state=520, seed=3):
    """A synthetic atlas in the real files' layout, large enough for the real constants."""
    from gliostrata.profiling.deconvolve import NEFTEL_SIGNATURES

    rng = np.random.default_rng(seed)
    markers = {"NPC": NEFTEL_SIGNATURES["NPC-like"], "OPC": NEFTEL_SIGNATURES["OPC-like"],
               "AC": NEFTEL_SIGNATURES["AC-like"], "MES": NEFTEL_SIGNATURES["MES-like"]}
    genes = [g for s in t1.STATES for g in markers[s]] + [f"BG{i}" for i in range(40)]
    rows, cols = [], {}
    for s in t1.STATES:
        for i in range(n_per_state):
            name = f"{s}cell{i}"
            tpm = rng.gamma(1.5, 2.0, len(genes))
            for g in markers[s]:
                tpm[genes.index(g)] *= 8.0
            tpm = tpm / tpm.sum() * 1e6
            cols[name] = np.log2(tpm / 10.0 + 1.0)
            score = {st: (2.0 if st == s else 0.0) for st in t1.STATES}
            rows.append((name, "TUM1", "yes", 0.0, score["MES"], score["AC"], score["OPC"],
                         score["NPC"], 0.0))
    expr = pd.DataFrame(cols, index=pd.Index(genes, name="GENE"))
    gz = tmp_path / "atlas.tsv.gz"
    with gzip.open(gz, "wt") as fh:
        expr.to_csv(fh, sep="\t")
    cells = _cells(rows).to_csv(index=False).encode()
    tar_path = tmp_path / "meta.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        info = tarfile.TarInfo(t1.ANNOTATION_MEMBER)
        info.size = len(cells)
        tar.addfile(info, io.BytesIO(cells))
    return gz, tar_path


def test_the_run_refuses_a_held_out_path_and_a_wrong_checksum(tmp_path):
    gz, tar_path = _write_synthetic_inputs(tmp_path, n_per_state=5)
    with pytest.raises(ValueError, match="reads no cohort"):
        t1.run(tmp_path / "absent" / "cgga" / "x.tsv.gz", tar_path,
               tmp_path / "o.json", tmp_path / "t.parquet", dry_run=True)
    with pytest.raises(ValueError, match="expected"):
        t1.run(gz, tar_path, tmp_path / "o.json", tmp_path / "t.parquet", dry_run=True)


def test_end_to_end_on_a_synthetic_atlas_then_refuses_to_run_again(tmp_path):
    gz, tar_path = _write_synthetic_inputs(tmp_path)
    out, table = tmp_path / "r" / "o.json", tmp_path / "r" / "t.parquet"
    dry = t1.run(gz, tar_path, out, table, dry_run=True, check_checksums=False)
    assert dry["counts"]["cells_per_state"] == {s: 520 for s in t1.STATES}
    assert not out.exists()
    rep = t1.run(gz, tar_path, out, table, check_checksums=False)
    assert rep["run_before_deposit"] is True
    assert rep["counts"]["markers_in_atlas"]["NPC-like"] == [15, 15]
    assert abs(rep["counts"]["cell_tpm_sum_min_max"][0] - 1e6) < 5
    p = rep["primary"]
    assert set(p["per_state_r"]) == set(t1.STATES)
    assert p["mean_r"] == pytest.approx(np.mean(list(p["per_state_r"].values())))
    assert rep["verdict"]["metric_key"] == t1.METRIC_KEY
    assert rep["verdict"]["status"] in {"PASS", "INCONCLUSIVE", "FAIL"}
    assert "secondary_descriptive" in rep
    saved = json.loads(out.read_text(encoding="utf-8"))
    assert saved["primary"]["mean_r"] == p["mean_r"]
    frame = pd.read_parquet(table)
    assert len(frame) == 100 and np.allclose(
        frame[[f"weight_{s}" for s in t1.STATES]].sum(axis=1), 1.0)
    with pytest.raises(FileExistsError, match="run once"):
        t1.run(gz, tar_path, out, table, check_checksums=False)


# --- 5. the real atlas, annotation only --------------------------------------------------
@pytest.mark.skipif(not t1.DEFAULT_ANNOTATION.is_file() and
                    not (ROOT / t1.DEFAULT_ANNOTATION).is_file(),
                    reason="the Neftel 2019 3CA annotation is not on disk")
def test_real_annotation_gives_the_counts_written_into_the_spec():
    path = ROOT / t1.DEFAULT_ANNOTATION
    assert t1.sha256_file(path) == t1.EXPECTED_SHA256["annotation"]
    cs = t1.assign_states(t1.read_cell_annotation(path))
    assert cs.attrs["n_annotated"] == 7930
    assert cs.attrs["n_malignant_with_scores"] == 6568
    assert cs.attrs["n_tied"] == 0
    assert cs["state"].value_counts().to_dict() == {
        "NPC": 1985, "AC": 1928, "MES": 1608, "OPC": 1047}
    assert cs["sample"].nunique() == 26


# --- 6. NOT EVALUABLE rules (spec section 5) -------------------------------------------
def test_an_undefined_full_set_r_makes_the_row_not_evaluable():
    rng = np.random.default_rng(4)
    truth = rng.dirichlet(np.ones(4), size=100)
    est = truth.copy()
    est[:, 1] = 0.25  # a constant estimate: r undefined
    primary = t1.score(est, truth, truth, np.random.default_rng(42))
    assert primary["mean_r"] is None and primary["bootstrap"] is None
    v = t1.final_verdict(primary)
    assert v["status"] == "NOT EVALUABLE" and "undefined" in v["reason"]


def test_too_many_undefined_resamples_make_the_row_not_evaluable():
    truth = np.zeros((100, 4))
    truth[:, :] = 0.25
    truth[0] = [0.4, 0.2, 0.2, 0.2]  # varies in one mixture only: most resamples undefined
    est = truth.copy()
    boot = t1.bootstrap(est, truth, np.random.default_rng(42), n_boot=200)
    assert boot["n_undefined_dropped"] > t1.MAX_UNDEFINED_RESAMPLES
    assert boot["mean_r"] is None
    primary = {"mean_r": 1.0, "per_state_r": {}, "bootstrap": boot}
    v = t1.final_verdict(primary)
    assert v["status"] == "NOT EVALUABLE" and "resamples" in v["reason"]
