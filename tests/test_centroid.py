"""Single-sample nearest-centroid subtype assignment (gliostrata.profiling.centroid).

What is at stake: a per-patient subtype feature must not depend on which other
patients were scored with it (stage 3b, D6), and the centroids it uses must be
the frozen, checksummed ones.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gliostrata.profiling.centroid import (
    CHECKSUM_FILE,
    MODE,
    SCALING,
    SET_ALL,
    SET_V2TRAIN,
    CentroidSet,
    centroid_set_to_json,
    classify_single_sample,
    load_centroid_set,
    sha256_bytes,
    transform_sample,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "subtype_centroids"
SUBTYPES = ("Classical", "Mesenchymal", "Proneural", "Neural")
TPM = ROOT / "data" / "raw" / "tcga_gbm" / "expression" / "tpm.parquet"
ANNOTATION = ROOT / "data" / "raw" / "tcga_gbm" / "expression" / "gene_annotation.parquet"


def _synthetic(n_genes=60, n_per_class=6, seed=0):
    """Expression whose classes differ on disjoint gene blocks; symbols as the index."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    block = n_genes // len(SUBTYPES)
    cols, labels, data = [], [], []
    for k, s in enumerate(SUBTYPES):
        for r in range(n_per_class):
            x = rng.gamma(2.0, 5.0, n_genes)
            x[k * block:(k + 1) * block] *= 20.0
            data.append(x)
            cols.append(f"P{s[0]}{r}")
            labels.append(s)
    expr = pd.DataFrame(np.column_stack(data), index=genes, columns=cols)
    return expr, pd.Series(labels, index=cols)


def _centroid_set(expr, labels, set_id="test_set"):
    z = pd.DataFrame({c: transform_sample(expr[c].to_numpy()) for c in expr.columns}, index=expr.index)
    cen = pd.DataFrame({s: z[labels.index[labels == s]].mean(axis=1) for s in SUBTYPES})
    return CentroidSet(set_id=set_id, genes=tuple(expr.index), subtypes=SUBTYPES, centroids=cen,
                       scaling=SCALING, source_sha256="", meta={"n_training": len(labels)},
                       ensembl_ids=("",) * len(expr.index))


def _write(tmp_path, cs, register=True, name=None):
    path = tmp_path / (name or f"{cs.set_id}.json")
    text = centroid_set_to_json(cs).encode("utf-8")
    path.write_bytes(text)
    if register:
        (tmp_path / CHECKSUM_FILE).write_text(json.dumps({path.name: sha256_bytes(text)}), encoding="utf-8")
    return path


# --- 1. the per-sample transform ----------------------------------------------
def test_transform_is_log2_then_z_across_genes_of_one_sample():
    v = np.array([0.0, 1.0, 3.0, 7.0])
    x = np.log2(v + 1)
    assert np.allclose(transform_sample(v), (x - x.mean()) / x.std(ddof=0))
    assert np.allclose(transform_sample(np.full(5, 4.0)), 0.0)  # flat sample: no division by zero


# --- 2. nearest centroid is correct on data with a known answer --------------
def test_synthetic_samples_go_to_their_own_class(tmp_path):
    train, tl = _synthetic(seed=1)
    test, sl = _synthetic(seed=2)
    cs = load_centroid_set(_write(tmp_path, _centroid_set(train, tl)))
    out = classify_single_sample(test, None, cs).set_index("patient_id")
    assert (out["subtype"] == sl).all()
    assert set(out["subtype_mode"]) == {MODE}
    assert set(out["subtype_centroid_sha256"]) == {cs.source_sha256}
    assert (out["subtype_n_genes_used"] == len(train.index)).all()
    corr = out[[f"subtype_corr_{s}" for s in SUBTYPES]]
    assert ((corr >= -1) & (corr <= 1)).all().all()


# --- 3. a sample's result does not depend on its batch ------------------------
def test_a_sample_with_no_finite_correlation_gets_no_call():
    """A zero-variance sample has no defined correlation with any centroid. Its
    subtype is missing, never the first subtype by default."""
    expr, labels = _synthetic()
    cs = _centroid_set(expr, labels)
    flat = expr.iloc[:, :2].copy()
    flat.iloc[:, 1] = 7.0  # constant over every frozen gene
    out = classify_single_sample(flat, None, cs)
    assert out.loc[0, "subtype"] in SUBTYPES
    assert pd.isna(out.loc[1, "subtype"])
    assert out.loc[1, [f"subtype_corr_{s}" for s in SUBTYPES]].isna().all()


def test_result_is_identical_alone_and_in_any_batch(tmp_path):
    train, tl = _synthetic(seed=3)
    test, _ = _synthetic(seed=4)
    cs = load_centroid_set(_write(tmp_path, _centroid_set(train, tl)))
    full = classify_single_sample(test, None, cs).set_index("patient_id")
    rng = np.random.default_rng(5)
    for sample in test.columns:
        alone = classify_single_sample(test[[sample]], None, cs).set_index("patient_id")
        pd.testing.assert_frame_equal(alone.loc[[sample]], full.loc[[sample]])
    for _ in range(5):
        cols = list(rng.choice(test.columns, size=7, replace=False))
        part = classify_single_sample(test[cols], None, cs).set_index("patient_id")
        pd.testing.assert_frame_equal(part, full.loc[cols])


def test_output_is_deterministic(tmp_path):
    train, tl = _synthetic(seed=6)
    cs = load_centroid_set(_write(tmp_path, _centroid_set(train, tl)))
    a = classify_single_sample(train, None, cs)
    b = classify_single_sample(train, None, cs)
    pd.testing.assert_frame_equal(a, b)


# --- 4. the frozen file is verified -------------------------------------------
def test_altered_centroid_file_is_refused(tmp_path):
    train, tl = _synthetic(seed=7)
    path = _write(tmp_path, _centroid_set(train, tl))
    load_centroid_set(path)
    text = path.read_text(encoding="utf-8").replace("1", "2", 1)
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="does not match the recorded"):
        load_centroid_set(path)


def test_unregistered_or_unverifiable_set_is_refused(tmp_path):
    train, tl = _synthetic(seed=8)
    cs = _centroid_set(train, tl)
    no_registry = _write(tmp_path, cs, register=False)
    with pytest.raises(ValueError, match="no checksum registry"):
        load_centroid_set(no_registry)
    _write(tmp_path, cs, name="other.json")  # registry now lists only other.json
    with pytest.raises(ValueError, match="not listed"):
        load_centroid_set(no_registry)


def test_too_few_frozen_genes_raises(tmp_path):
    train, tl = _synthetic(seed=9)
    cs = load_centroid_set(_write(tmp_path, _centroid_set(train, tl)))
    with pytest.raises(ValueError, match="frozen genes"):
        classify_single_sample(train.rename(index=lambda g: "X" + g), None, cs)


def test_ensembl_input_is_matched_by_the_frozen_id_not_a_batch_rule(tmp_path):
    """With an annotation, genes are found by the Ensembl id frozen in the set."""
    train, tl = _synthetic(seed=10)
    base = _centroid_set(train, tl)
    ens = tuple(f"ENSG{i:011d}" for i in range(len(base.genes)))
    cs = CentroidSet(set_id="ens", genes=base.genes, subtypes=base.subtypes, centroids=base.centroids,
                     scaling=base.scaling, source_sha256="", meta={}, ensembl_ids=ens)
    cs = load_centroid_set(_write(tmp_path, cs))
    expr = train.copy()
    expr.index = [f"{e}.7" for e in ens]  # versioned ids, as in the GDC matrices
    ann = pd.DataFrame({"gene_id": list(expr.index), "gene_name": list(base.genes)})
    a = classify_single_sample(expr, ann, cs)
    b = classify_single_sample(train, None, cs)
    pd.testing.assert_frame_equal(a, b)


# --- 5. the committed sets -----------------------------------------------------
def test_committed_sets_match_their_registered_checksums():
    registry = json.loads((CONFIG_DIR / CHECKSUM_FILE).read_text(encoding="utf-8"))
    assert set(registry) == {f"{SET_V2TRAIN}.json", f"{SET_ALL}.json"}
    for set_id in (SET_V2TRAIN, SET_ALL):
        cs = load_centroid_set(CONFIG_DIR / f"{set_id}.json")
        assert cs.set_id == set_id
        assert cs.subtypes == SUBTYPES
        assert len(cs.genes) == 831 == len(set(cs.ensembl_ids))
        assert cs.scaling == SCALING
    tr = load_centroid_set(CONFIG_DIR / f"{SET_V2TRAIN}.json").meta
    al = load_centroid_set(CONFIG_DIR / f"{SET_ALL}.json").meta
    assert tr["n_labelled_by_split"] == {"train": 170, "test": 76}
    assert tr["n_training"] == 170 and al["n_training"] == 246
    assert tr["partition"] == "tcga_gbm_clinical_v2"


def test_committed_gene_map_accounts_for_all_840_published_genes():
    doc = json.loads((CONFIG_DIR / "clanc840_genes.json").read_text(encoding="utf-8"))
    assert len(doc["genes"]) == 840
    assert doc["counts"] == {"ambiguous": 2, "direct": 737, "hgnc_previous_symbol": 94,
                             "not_in_annotation": 7}
    assert doc["n_usable"] == 831
    assert all((g["symbol"] is None) == (g["resolution"] in ("ambiguous", "not_in_annotation"))
               for g in doc["genes"])


def test_gene_resolution_rule_never_guesses():
    spec = importlib.util.spec_from_file_location("prep", ROOT / "scripts" / "prepare_subtype_reference.py")
    prep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prep)
    ann = pd.DataFrame({"gene_name": ["EGFR", "S1PR1", "HEXIM1", "NEW1"]})
    hgnc = pd.DataFrame({"Symbol": ["EDG1", "EDG1", "OLD1", "OLD1"],
                         "Approved.Symbol": ["S1PR1", "HEXIM1", "NEW1", "NEW1"]})
    out = {r["published_symbol"]: r for r in prep.resolve_genes(["EGFR", "EDG1", "OLD1", "LOC9"], ann, hgnc)}
    assert out["EGFR"]["resolution"] == "direct"
    assert out["EDG1"]["resolution"] == "ambiguous" and out["EDG1"]["symbol"] is None
    assert out["OLD1"] == {"published_symbol": "OLD1", "symbol": "NEW1", "resolution": "hgnc_previous_symbol"}
    assert out["LOC9"]["resolution"] == "not_in_annotation"


# --- 6. real data --------------------------------------------------------------
@pytest.mark.skipif(not (TPM.exists() and ANNOTATION.exists()),
                    reason="TCGA-GBM expression parquet not present")
def test_real_tcga_single_sample_calls_do_not_depend_on_the_batch():
    expr = pd.read_parquet(TPM)
    ann = pd.read_parquet(ANNOTATION)
    cs = load_centroid_set(CONFIG_DIR / f"{SET_V2TRAIN}.json")
    cols = list(expr.columns[:40])
    joint = classify_single_sample(expr[cols], ann, cs).set_index("patient_id")
    assert (joint["subtype_n_genes_used"] == 831).all()
    halves = pd.concat([classify_single_sample(expr[cols[:13]], ann, cs),
                        classify_single_sample(expr[cols[13:]], ann, cs)]).set_index("patient_id")
    pd.testing.assert_frame_equal(halves, joint)
    alone = classify_single_sample(expr[[cols[5]]], ann, cs).set_index("patient_id")
    pd.testing.assert_frame_equal(alone, joint.loc[[cols[5]]])
