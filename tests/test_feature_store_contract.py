"""Tests for the L3 feature-store contract (decisions D1 to D5 of stage 3b).

One test (or a small group) per rule:

  * a column supplied by more than one input raises;
  * a null-only column that shadows a computed feature raises;
  * a missing resection is not "biopsy";
  * a missing molecular call is not 0;
  * call status has three states;
  * KPS is never imputed in the store;
  * the counts of each named view, per split;
  * a stub consumed in real mode raises;
  * a label column in a feature matrix raises;
  * an integration test on a small real TCGA fixture, built at test time from
    the real files and skipped when they are absent (no TCGA value is
    committed).

All other tables are synthetic. No test reads CGGA.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from gliostrata.data.clinical import harmonise_clinical
from gliostrata.features.build import molecular_from_call_tables
from gliostrata.features.contract import (
    CALL_STATUSES,
    MISSING_REASONS,
    ColumnOwnershipError,
    DuplicateSourceError,
    FeatureStoreError,
    LabelInFeatureMatrixError,
    ShadowColumnError,
    StubInRealModeError,
    UnexplainedMissingError,
    sha256_file,
    sibling_paths,
    write_store,
)
from gliostrata.features.store import (
    build_feature_store,
    feature_matrix,
    stub_layer_frame,
)
from gliostrata.features.views import VIEWS, build_view, partition_lookup

ROOT = Path(__file__).resolve().parents[1]
SUBTYPES = ("Classical", "Mesenchymal", "Proneural", "Neural")


# ---------------------------------------------------------------------------
# A small synthetic set of inputs, shaped like the real ones
# ---------------------------------------------------------------------------
def _ids(n):
    return [f"P{i:02d}" for i in range(n)]


def _inputs(n=8):
    """Named inputs for ``n`` patients.

    P00..P{n-3} have an expression sample, the last two do not; one sample
    ("OUT") is outside the cohort. Resection is absent for everyone. KPS is
    missing for P01; treatment arm for P02; vital status for P03. MGMT: P04
    not tested, P05 tested but not callable. IDH: P06 not tested.
    """
    ids = _ids(n)
    expr = ids[: n - 2]
    clinical = pd.DataFrame({
        "patient_id": ids,
        "age": [50.0 + i for i in range(n)],
        "sex": ["M" if i % 2 else "F" for i in range(n)],
        "kps": [80.0 if i != 1 else None for i in range(n)],
        "treatment_arm": ["ChemoRadiation" if i != 2 else None for i in range(n)],
        "extent_of_resection": [None] * n,
    })
    clinical["reason__age"] = None
    clinical["reason__sex"] = None
    clinical["reason__kps"] = ["not_recorded_at_source" if i == 1 else None for i in range(n)]
    clinical["reason__treatment_arm"] = [
        "not_recorded_at_source" if i == 2 else None for i in range(n)]
    clinical["reason__extent_of_resection"] = "not_available_in_cohort"
    labels = pd.DataFrame({
        "patient_id": ids,
        "overall_survival_months": [10.0 + i for i in range(n)],
        "vital_status": [1 if i != 3 else None for i in range(n)],
        "reason__overall_survival_months": None,
        "reason__vital_status": ["recorded_unknown_at_source" if i == 3 else None
                                 for i in range(n)],
    })
    subtypes = pd.DataFrame({"patient_id": expr})
    subtypes["subtype"] = [SUBTYPES[i % 4] for i in range(len(expr))]
    for j, s in enumerate(SUBTYPES):
        subtypes[f"subtype_corr_{s}"] = [0.1 * j + 0.01 * i for i in range(len(expr))]
    subtypes["subtype_mode"] = "single_sample"
    subtypes["subtype_centroid_set"] = "test_set"
    subtypes["subtype_centroid_sha256"] = "0" * 64
    subtypes["subtype_n_genes_used"] = 800
    states = pd.DataFrame({"patient_id": expr})
    for c in ("state_npc", "state_opc", "state_ac", "state_mes"):
        states[c] = 0.25
    states["state_deconvolution_method"] = "signature_ssgsea_softmax_z"
    states["state_deconvolution_version"] = "0.1.0"
    states["is_true_deconvolution"] = False
    tme = pd.DataFrame({"patient_id": expr})
    for c in ("tme_immune", "tme_stromal", "tme_endothelial"):
        tme[c] = 1 / 3
    tme["tme_deconvolution_method"] = "signature_ssgsea_softmax_z"
    tme["tme_deconvolution_version"] = "0.1.0"
    tme["is_true_deconvolution"] = False
    mgmt_status = ["not_tested" if i == 4 else "tested_not_callable" if i == 5 else "called"
                   for i in range(n)]
    idh_status = ["not_tested" if i == 6 else "called" for i in range(n)]
    molecular = pd.DataFrame({
        "patient_id": ids,
        "mgmt_score": [0.5 if s == "called" else None for s in mgmt_status],
        "mgmt_threshold": [0.35 if s != "not_tested" else None for s in mgmt_status],
        "mgmt_methylated": [1 if s == "called" else None for s in mgmt_status],
        "mgmt_call_status": mgmt_status,
        "idh_mutated": [0 if s == "called" else None for s in idh_status],
        "idh_call_status": idh_status,
    })
    return {"clinical": clinical, "labels": labels, "subtypes": subtypes,
            "states": states, "tme": tme, "molecular": molecular,
            "expression_samples": [*expr, "OUT"]}


def _store(mode="real", **overrides):
    kw = _inputs()
    kw.update(overrides)
    return build_feature_store(mode=mode, **kw)


# ---------------------------------------------------------------------------
# D1: every patient, one row; one reason per null cell; labels separate
# ---------------------------------------------------------------------------
def test_the_store_keeps_every_patient_one_row_each():
    store = _store()
    assert store.patient_ids == _ids(8)
    assert store.labels["patient_id"].tolist() == _ids(8)
    assert store.meta["n_expression_samples_outside_cohort"] == 1


def test_every_null_cell_has_exactly_one_reason_from_the_vocabulary():
    store = _store()
    n_null = int(store.features.drop(columns="patient_id").isna().sum().sum()
                 + store.labels.drop(columns="patient_id").isna().sum().sum())
    assert len(store.missing) == n_null
    assert not store.missing.duplicated(["patient_id", "column"]).any()
    assert set(store.missing["reason"]) <= set(MISSING_REASONS)
    # every listed cell really is null, and nothing else is
    for row in store.missing.itertuples():
        frame = store.labels if row.column in ("vital_status",
                                               "overall_survival_months") else store.features
        value = frame.loc[frame["patient_id"] == row.patient_id, row.column].iloc[0]
        assert pd.isna(value)


def test_an_untraced_null_raises():
    kw = _inputs()
    kw["clinical"].loc[0, "age"] = None  # missing, and no reason given
    with pytest.raises(UnexplainedMissingError):
        build_feature_store(mode="real", **kw)


def test_l2_values_of_patients_without_expression_are_coded_no_expression_sample():
    store = _store()
    m = store.missing
    for col in ("subtype", "state_npc", "tme_immune", "subtype_corr_Classical"):
        rows = m[m["column"] == col]
        assert set(rows["patient_id"]) == {"P06", "P07"}
        assert set(rows["reason"]) == {"no_expression_sample"}


def test_every_column_carries_role_layer_source_files_and_class():
    srcs = {"clinical": [{"path": "c.parquet", "sha256": "a" * 64}]}
    store = _store(sources=srcs)
    records = {r["name"]: r for r in store.column_records()}
    assert set(records) == {*store.features.columns, *store.labels.columns} - set()
    for rec in records.values():
        assert rec["role"] in ("id", "feature", "status", "label")
        assert rec["class"] in ("measured", "derived_from_measurement",
                                "clinical_record", "not_available")
        assert rec["source_layer"]
    assert records["age"]["source_files"] == srcs["clinical"]
    assert records["vital_status"]["role"] == "label"
    assert records["treatment_arm"]["role"] == "feature"  # a model feature, not an outcome


def test_outcome_labels_live_only_in_the_label_table():
    store = _store()
    assert "vital_status" not in store.features.columns
    assert "overall_survival_months" not in store.features.columns
    assert list(store.labels.columns) == ["patient_id", "overall_survival_months",
                                          "vital_status"]


# ---------------------------------------------------------------------------
# D5: column ownership
# ---------------------------------------------------------------------------
def test_a_column_supplied_by_two_inputs_raises():
    kw = _inputs()
    kw["clinical"]["state_npc"] = 0.3  # a value, also supplied by the states input
    with pytest.raises(DuplicateSourceError, match="state_npc by clinical, states"):
        build_feature_store(mode="real", **kw)


def test_a_null_only_column_that_shadows_a_computed_feature_raises():
    kw = _inputs()
    kw["clinical"]["subtype"] = None  # the L0 placeholder pattern
    with pytest.raises(ShadowColumnError, match="subtype"):
        build_feature_store(mode="real", **kw)


def test_a_null_only_shadow_raises_even_without_a_duplicate():
    """A null-only computed-feature name in a foreign input is refused on its own."""
    kw = _inputs()
    kw["molecular"]["tme_immune"] = None
    with pytest.raises(ShadowColumnError):
        build_feature_store(mode="real", **kw)


def test_the_l0_cohort_placeholders_are_refused_as_a_clinical_input():
    """Passing the L0 cohort straight in, placeholders and all, cannot work."""
    kw = _inputs()
    clinical = kw["clinical"].copy()
    for col in ("mgmt_methylated", "idh_mutated", "state_npc", "tme_stromal"):
        clinical[col] = None
    kw["clinical"] = clinical
    with pytest.raises(ShadowColumnError):
        build_feature_store(mode="real", **kw)


def test_a_column_an_input_does_not_own_raises():
    kw = _inputs()
    kw["states"]["some_other_value"] = 1.0
    with pytest.raises(ColumnOwnershipError):
        build_feature_store(mode="real", **kw)


def test_harmonise_clinical_emits_no_l2_names_and_refuses_one_with_values():
    cohort = pd.DataFrame({
        "patient_id": ["A", "B"], "age": [50.0, 60.0], "sex": ["M", "F"],
        "kps": [80.0, None], "extent_of_resection": [None, None],
        "treatment_arm": ["ChemoRadiation", None],
        "overall_survival_months": [10.0, 20.0], "vital_status": [1.0, 0.0],
        "subtype": [None, None], "state_npc": [None, None],
        "mgmt_methylated": [None, None], "idh_mutated": [None, None],
    })
    clinical, labels = harmonise_clinical(cohort)
    assert set(clinical.columns) == {"patient_id", "age", "sex", "kps", "treatment_arm",
                                     "extent_of_resection", "reason__age", "reason__sex",
                                     "reason__kps", "reason__treatment_arm",
                                     "reason__extent_of_resection"}
    assert "vital_status" in labels.columns and "vital_status" not in clinical.columns
    bad = cohort.copy()
    bad.loc[0, "state_npc"] = 0.4
    with pytest.raises(ValueError, match="state_npc"):
        harmonise_clinical(bad)
    unknown = cohort.assign(something_new=[1, 2])
    with pytest.raises(ValueError, match="something_new"):
        harmonise_clinical(unknown)


def test_harmonise_clinical_traces_each_missing_value_to_the_raw_record():
    cohort = pd.DataFrame({
        "patient_id": ["A", "B", "C", "D"], "age": [50.0] * 4, "sex": ["M"] * 4,
        "kps": [None, None, None, 90.0], "extent_of_resection": [None] * 4,
        "treatment_arm": [None, "MonoTherapy", "MonoTherapy", "MonoTherapy"],
        "overall_survival_months": [1.0] * 4, "vital_status": [1.0, None, 1.0, 1.0],
    })
    raw = pd.DataFrame({
        "submitter_id": ["A", "B", "C", "D"], "age_at_index_raw": [50.0] * 4,
        "sex_at_birth_raw": ["male"] * 4, "kps_raw": [None, "Unknown", "70", "90"],
        "treatment_types": ["[]", '["Radiation Therapy, NOS"]'] + ['["x"]'] * 2,
        "days_to_death": [30.0] * 4, "days_to_follow_up": [None] * 4,
        "vital_status_raw": ["Dead", "Not Reported", "Dead", "Dead"],
    })
    clinical, labels = harmonise_clinical(cohort, raw)
    r = dict(zip(clinical["patient_id"], clinical["reason__kps"], strict=True))
    assert r["A"] == "not_recorded_at_source"
    assert r["B"] == "recorded_unknown_at_source"
    assert r["C"] == "lost_in_harmonisation"  # "70" at source, missing after
    assert pd.isna(r["D"])
    assert clinical.loc[0, "reason__treatment_arm"] == "not_recorded_at_source"
    assert set(clinical["reason__extent_of_resection"]) == {"not_available_in_cohort"}
    assert labels.loc[1, "reason__vital_status"] == "recorded_unknown_at_source"


# ---------------------------------------------------------------------------
# D3: absent never becomes a value
# ---------------------------------------------------------------------------
def test_missing_resection_is_not_biopsy():
    store = _store()
    col = store.features["extent_of_resection"]
    assert col.isna().all()
    assert "biopsy" not in set(col.dropna().astype(str))
    assert store.columns["extent_of_resection"].klass == "not_available"
    reasons = store.missing[store.missing["column"] == "extent_of_resection"]
    assert set(reasons["reason"]) == {"not_available_in_cohort"}
    # and nothing downstream may turn it into the one-hot reference level
    with pytest.raises(FeatureStoreError):
        stub_layer_frame(store)


def test_a_missing_call_is_not_zero():
    store = _store()
    f = store.features.set_index("patient_id")
    assert str(f["mgmt_methylated"].dtype) == "boolean"
    assert pd.isna(f.loc["P04", "mgmt_methylated"])   # not tested
    assert pd.isna(f.loc["P05", "mgmt_methylated"])   # tested, not callable
    assert pd.isna(f.loc["P06", "idh_mutated"])       # not tested
    assert f["idh_mutated"].dropna().eq(False).all()  # the called ones are wildtype here
    assert int(f["mgmt_methylated"].isna().sum()) == 2
    assert int(f["idh_mutated"].isna().sum()) == 1


def test_call_status_has_three_states_and_keeps_score_and_threshold():
    store = _store()
    f = store.features.set_index("patient_id")
    assert set(f["mgmt_call_status"]) == set(CALL_STATUSES)
    assert f.loc["P05", "mgmt_call_status"] == "tested_not_callable"
    assert f.loc["P04", "mgmt_call_status"] == "not_tested"
    assert f.loc["P00", "mgmt_score"] == pytest.approx(0.5)
    assert f.loc["P00", "mgmt_threshold"] == pytest.approx(0.35)
    m = store.missing
    assert m[(m["column"] == "mgmt_methylated") & (m["patient_id"] == "P05")
             ]["reason"].tolist() == ["tested_not_callable"]
    assert m[(m["column"] == "mgmt_methylated") & (m["patient_id"] == "P04")
             ]["reason"].tolist() == ["not_tested"]


def test_molecular_from_call_tables_maps_the_three_states():
    mgmt = pd.DataFrame({
        "patient_id": ["A", "B", "C", "X"],
        "mgmt_methylated": [1, None, None, 0], "mgmt_score": [0.6, None, None, 0.1],
        "threshold": [0.35, 0.35, None, 0.35],
        "call_basis": ["two_probes_measured", "required_probe_absent_or_NA",
                       "no_methylation_data", "two_probes_measured"],
    })
    idh = pd.DataFrame({
        "patient_id": ["A", "B"], "idh_mutated": [1, 0],
        "call_basis": ["hotspot_variant", "sequenced_no_hotspot"],
    })
    mol, counts = molecular_from_call_tables(["A", "B", "C", "D"], mgmt, idh)
    assert counts["idh"]["assay_manifest_used"] is False
    assert mol["mgmt_call_status"].tolist() == ["called", "tested_not_callable",
                                                "not_tested", "not_tested"]
    assert mol["idh_call_status"].tolist() == ["called", "called", "not_tested",
                                               "not_tested"]
    assert counts["mgmt"]["n_rows_outside_cohort"] == 1
    assert counts["mgmt"]["n_cohort_patients_absent_from_table"] == 1
    assert counts["idh"]["n_cohort_patients_absent_from_table"] == 2
    bad = mgmt.assign(call_basis=["something_else"] * 4)
    with pytest.raises(FeatureStoreError, match="unknown call_basis"):
        molecular_from_call_tables(["A"], bad, idh)


def test_a_patient_assayed_only_on_a_non_primary_sample_is_tested_not_callable():
    """The call tables call the primary tumour only. A patient whose only
    sequenced sample is a recurrence appears there as ``no_mutation_data``; the
    file manifest shows it was tested, so it is ``tested_not_callable``."""
    mgmt = pd.DataFrame({
        "patient_id": ["A", "B"], "mgmt_methylated": [1, None],
        "mgmt_score": [0.6, None], "threshold": [0.35, None],
        "call_basis": ["two_probes_measured", "no_methylation_data"],
    })
    idh = pd.DataFrame({
        "patient_id": ["A", "B"], "idh_mutated": [0, None],
        "call_basis": ["sequenced_no_hotspot", "no_mutation_data"],
    })
    assayed = {"mgmt": {"A"}, "idh": {"A", "B"}}  # B: a recurrent sample only
    mol, counts = molecular_from_call_tables(["A", "B"], mgmt, idh, assayed=assayed)
    assert mol["idh_call_status"].tolist() == ["called", "tested_not_callable"]
    assert mol["mgmt_call_status"].tolist() == ["called", "not_tested"]
    assert counts["idh"]["n_assayed_but_no_primary_call"] == 1
    assert counts["mgmt"]["n_assayed_but_no_primary_call"] == 0
    assert pd.isna(mol.loc[1, "idh_mutated"])  # still no call, never 0


def test_a_call_value_without_called_status_raises():
    kw = _inputs()
    kw["molecular"].loc[4, "mgmt_methylated"] = 0  # value although not tested
    with pytest.raises(FeatureStoreError, match="mgmt_methylated"):
        build_feature_store(mode="real", **kw)


# ---------------------------------------------------------------------------
# D4: KPS is never imputed
# ---------------------------------------------------------------------------
def test_kps_is_never_imputed_in_the_store():
    store = _store()
    f = store.features.set_index("patient_id")
    assert pd.isna(f.loc["P01", "kps"])
    assert int(f["kps"].isna().sum()) == 1
    assert "kps_missing" not in store.features.columns
    reasons = store.missing[store.missing["column"] == "kps"]
    assert reasons["reason"].tolist() == ["not_recorded_at_source"]
    with pytest.raises(FeatureStoreError, match="kps"):
        feature_matrix(store, ["age", "kps"])
    for view in VIEWS.values():
        assert "kps" not in view.required_columns()


# ---------------------------------------------------------------------------
# D2: named views and their counts
# ---------------------------------------------------------------------------
def test_the_counts_of_each_view_with_per_split_counts():
    store = _store()
    split_of = partition_lookup(["P00", "P01", "P02", "P03", "P04"], ["P05", "P06"])
    prof = build_view(store, "profiling", split_of)
    assert prof.patient_ids == _ids(6)
    assert prof.flow[-1]["n_remaining"] == 6
    assert prof.flow[-1]["n_remaining_train"] == 5
    assert prof.flow[-1]["n_remaining_test"] == 1
    assert prof.flow[-1]["n_remaining_not_in_partition"] == 0
    trans = build_view(store, "transferable_model", split_of)
    # P02 no treatment arm, P03 no vital status
    assert trans.patient_ids == ["P00", "P01", "P04", "P05", "P06", "P07"]
    steps = {s["step"]: s for s in trans.flow}
    assert steps["treatment arm recorded"]["n_excluded"] == 1
    assert steps["vital status recorded"]["n_excluded"] == 1
    assert steps["vital status recorded"]["n_excluded_train"] == 1
    assert trans.flow[-1]["n_remaining_not_in_partition"] == 1  # P07
    mol = build_view(store, "molecular_model", split_of)
    # and P04 / P05 lack an MGMT call, P06 an IDH call
    assert mol.patient_ids == ["P00", "P01", "P07"]
    steps = {s["step"]: s for s in mol.flow}
    assert steps["IDH called"]["n_excluded"] == 1
    assert steps["MGMT called"]["n_excluded"] == 2
    assert steps["MGMT called"]["n_missing_overall"] == 2
    # sequential counts add up
    for result in (prof, trans, mol):
        assert 8 - sum(s["n_excluded"] for s in result.flow) == len(result.patient_ids)
    assert list(mol.features.columns) == ["patient_id", "age", "sex", "treatment_arm",
                                          "idh_mutated", "mgmt_methylated"]
    assert mol.labels is not None and "vital_status" in mol.labels.columns
    assert "vital_status" not in mol.features.columns


def test_a_view_never_compares_outcomes():
    """Permuting the outcome values leaves every view's membership unchanged."""
    store = _store()
    kw = _inputs()
    lab = kw["labels"]
    lab["overall_survival_months"] = lab["overall_survival_months"].to_numpy()[::-1]
    permuted = build_feature_store(mode="real", **kw)
    for name in VIEWS:
        assert build_view(store, name).patient_ids == build_view(permuted, name).patient_ids


# ---------------------------------------------------------------------------
# Real mode refuses stubs and cohort-relative subtype outputs
# ---------------------------------------------------------------------------
def test_a_stub_subtype_output_consumed_in_real_mode_raises():
    from gliostrata.profiling.subtype import classify_subtypes

    kw = _inputs()
    fake = pd.DataFrame({"patient_id": kw["subtypes"]["patient_id"],
                         "subtype": kw["subtypes"]["subtype"]})
    stub = classify_subtypes(fake)
    with pytest.raises(StubInRealModeError):
        build_feature_store(mode="real", **{**kw, "subtypes": stub})
    # the same output is accepted, labelled, in synthetic mode
    kw2 = _inputs()
    store = build_feature_store(mode="synthetic", **{**kw2, "subtypes": stub})
    assert store.mode == "synthetic"
    assert store.columns["subtype"].source_layer.startswith("synthetic:")


def test_a_stub_state_output_consumed_in_real_mode_raises():
    kw = _inputs()
    kw["states"]["state_deconvolution_method"] = "stub_passthrough"
    kw["states"]["state_deconvolution_version"] = "0.1.0-stub"
    with pytest.raises(StubInRealModeError, match="stub"):
        build_feature_store(mode="real", **kw)


def test_a_cohort_relative_subtype_output_is_refused_in_real_mode():
    kw = _inputs()
    sub = kw["subtypes"].drop(columns=["subtype_mode", "subtype_centroid_set",
                                       "subtype_centroid_sha256"])
    with pytest.raises(StubInRealModeError, match="single_sample"):
        build_feature_store(mode="real", **{**kw, "subtypes": sub})


def test_real_mode_cli_without_a_centroid_set_fails_loudly_and_writes_nothing(tmp_path):
    from gliostrata.cli import cli

    out = tmp_path / "store.parquet"
    result = CliRunner().invoke(cli, [
        "build-feature-store", "--mode", "real",
        "--centroids", str(tmp_path / "absent.json"), "--out", str(out),
        "--report", str(tmp_path / "report.json"),
    ])
    assert result.exit_code != 0
    assert "centroid" in result.output
    assert not out.exists() and not (tmp_path / "report.json").exists()


def test_the_cli_refuses_a_cgga_path_before_opening_anything(tmp_path):
    """The paths do not exist, so only the refusal, not a missing file, can
    produce the refusal text; and no real CGGA path is ever named here."""
    from gliostrata.cli import cli

    scratch = tmp_path / "scratch"
    for args in (["--mode", "real", "--expression", str(scratch / "cgga" / "x.parquet")],
                 ["--mode", "synthetic", "--input", str(scratch / "CGGA_x.parquet"),
                  "--out", str(tmp_path / "s.parquet")]):
        result = CliRunner().invoke(cli, ["build-feature-store", *args,
                                          "--report", str(tmp_path / "r.json")])
        assert result.exit_code != 0
        assert "refusing" in result.output and "held-out cohort" in result.output
        assert "not found" not in result.output
    result = CliRunner().invoke(cli, [
        "harmonize-clinical", "--cohort", str(scratch / "Cgga" / "c.parquet"),
        "--raw", str(scratch / "r.parquet"),
        "--clinical-out", str(tmp_path / "c.parquet"),
        "--labels-out", str(tmp_path / "l.parquet")])
    assert result.exit_code != 0 and "held-out cohort" in result.output
    assert not list(tmp_path.iterdir())


# ---------------------------------------------------------------------------
# Labels are refused in a feature matrix
# ---------------------------------------------------------------------------
def test_a_label_column_in_a_feature_matrix_raises():
    store = _store()
    with pytest.raises(LabelInFeatureMatrixError):
        feature_matrix(store, ["age", "vital_status"])
    with pytest.raises(LabelInFeatureMatrixError):
        feature_matrix(store.features, ["overall_survival_months"])
    view = build_view(store, "transferable_model")
    matrix = feature_matrix(view.features, ["age"])
    assert list(matrix.columns) == ["age"]


def test_an_outcome_supplied_by_the_clinical_input_is_refused():
    kw = _inputs()
    kw["clinical"]["vital_status"] = 1
    with pytest.raises((DuplicateSourceError, ColumnOwnershipError)):
        build_feature_store(mode="real", **kw)
    kw = _inputs()
    kw["labels"] = kw["labels"].drop(columns=["vital_status", "reason__vital_status"])
    kw["clinical"]["vital_status"] = 1
    with pytest.raises(ColumnOwnershipError, match="outcome labels"):
        build_feature_store(mode="real", **kw)


# ---------------------------------------------------------------------------
# Files: deterministic
# ---------------------------------------------------------------------------
def test_writing_the_store_twice_gives_byte_identical_files(tmp_path):
    a = write_store(_store(), sibling_paths(tmp_path / "a" / "store.parquet"))
    b = write_store(_store(), sibling_paths(tmp_path / "b" / "store.parquet"))
    assert [r["sha256"] for r in a.values()] == [r["sha256"] for r in b.values()]
    assert set(a) == {"features", "labels", "missing", "columns"}


# ---------------------------------------------------------------------------
# Integration: a small real TCGA fixture, built at test time
# ---------------------------------------------------------------------------
REAL = {
    "cohort": ROOT / "data/raw/tcga_gbm/clinical_cohort.parquet",
    "raw": ROOT / "data/raw/tcga_gbm/clinical_raw.parquet",
    "tpm": ROOT / "data/raw/tcga_gbm/expression/tpm.parquet",
    "annotation": ROOT / "data/raw/tcga_gbm/expression/gene_annotation.parquet",
    "mgmt": ROOT / "data/raw/tcga_gbm/methylation/mgmt_calls.parquet",
    "idh": ROOT / "data/raw/tcga_gbm/mutations/idh_calls.parquet",
    "mgmt_manifest": ROOT / "data/raw/tcga_gbm/methylation/file_manifest.parquet",
    "idh_manifest": ROOT / "data/raw/tcga_gbm/mutations/file_manifest.parquet",
    "centroids": ROOT / "configs/subtype_centroids/verhaak840_tcga_v2train.json",
    "partition": ROOT / "data/partitions/tcga_gbm_clinical_v2.json",
}
_HAVE_REAL = all(p.is_file() for p in REAL.values())
needs_real = pytest.mark.skipif(not _HAVE_REAL, reason="real TCGA inputs or the "
                                "frozen centroid set are not on disk")


def _small_real_fixture(tmp_path, n_with=6, n_without=4):
    """Real TCGA rows for a few patients, written under tmp_path (never committed)."""
    cohort = pd.read_parquet(REAL["cohort"])
    raw = pd.read_parquet(REAL["raw"])
    tpm = pd.read_parquet(REAL["tpm"])
    ids = cohort["patient_id"].astype(str)
    with_expr = [p for p in ids if p in set(tpm.columns)][:n_with]
    without = [p for p in ids if p not in set(tpm.columns)][:n_without]
    keep = [p for p in ids if p in set(with_expr) | set(without)]
    paths = {}
    cohort[cohort["patient_id"].isin(keep)].to_parquet(tmp_path / "cohort.parquet")
    raw[raw["submitter_id"].isin(keep)].to_parquet(tmp_path / "raw.parquet")
    tpm[with_expr].to_parquet(tmp_path / "tpm.parquet")
    for name in ("cohort", "raw", "tpm"):
        paths[name] = tmp_path / f"{name}.parquet"
    return paths, keep, with_expr


@needs_real
def test_integration_on_a_small_real_tcga_fixture(tmp_path):
    from gliostrata.cli import cli

    paths, keep, with_expr = _small_real_fixture(tmp_path)
    runner = CliRunner()
    r1 = runner.invoke(cli, [
        "harmonize-clinical", "--cohort", str(paths["cohort"]), "--raw", str(paths["raw"]),
        "--clinical-out", str(tmp_path / "clinical.parquet"),
        "--labels-out", str(tmp_path / "labels.parquet"),
    ])
    assert r1.exit_code == 0, r1.output
    outs = {}
    for run in ("a", "b"):
        out = tmp_path / run / "store.parquet"
        r2 = runner.invoke(cli, [
            "build-feature-store", "--mode", "real",
            "--clinical", str(tmp_path / "clinical.parquet"),
            "--labels", str(tmp_path / "labels.parquet"),
            "--expression", str(paths["tpm"]), "--annotation", str(REAL["annotation"]),
            "--mgmt-calls", str(REAL["mgmt"]), "--idh-calls", str(REAL["idh"]),
            "--mgmt-manifest", str(REAL["mgmt_manifest"]),
            "--idh-manifest", str(REAL["idh_manifest"]),
            "--centroids", str(REAL["centroids"]), "--partition", str(REAL["partition"]),
            "--out", str(out), "--report", str(tmp_path / run / "report.json"),
        ])
        assert r2.exit_code == 0, r2.output
        outs[run] = {k: p.read_bytes() for k, p in sibling_paths(out).items()}
    assert outs["a"] == outs["b"], "a repeat run is not byte-identical"

    store = pd.read_parquet(tmp_path / "a" / "store.parquet")
    assert store["patient_id"].tolist() == keep  # every patient, cohort order
    assert store["expression_available"].sum() == len(with_expr)
    assert store["extent_of_resection"].isna().all()
    profiled = store[store["expression_available"]]
    assert profiled["subtype"].notna().all()
    assert store.loc[~store["expression_available"], "subtype"].isna().all()
    missing = pd.read_parquet(tmp_path / "a" / "store_missing.parquet")
    n_null = int(store.drop(columns="patient_id").isna().sum().sum())
    labels = pd.read_parquet(tmp_path / "a" / "store_labels.parquet")
    n_null += int(labels.drop(columns="patient_id").isna().sum().sum())
    assert len(missing) == n_null
    assert str(labels["vital_status"].dtype) == "Int8"


@needs_real
def test_the_full_tcga_store_reproduces_the_t2_fit_sets(tmp_path):
    """The two model views reproduce T2's pinned fit sets (570 and 298)."""
    import json

    from gliostrata.features.build import run_harmonize_clinical, run_real

    run_harmonize_clinical(str(REAL["cohort"]), str(REAL["raw"]),
                           str(tmp_path / "c.parquet"), str(tmp_path / "l.parquet"))
    report = run_real(
        clinical=str(tmp_path / "c.parquet"), labels=str(tmp_path / "l.parquet"),
        expression=str(REAL["tpm"]), annotation=str(REAL["annotation"]),
        mgmt_calls=str(REAL["mgmt"]), idh_calls=str(REAL["idh"]),
        mgmt_manifest=str(REAL["mgmt_manifest"]), idh_manifest=str(REAL["idh_manifest"]),
        centroids=str(REAL["centroids"]), partition=str(REAL["partition"]),
        out=str(tmp_path / "store.parquet"),
    )
    fp = json.loads((ROOT / "tests/fixtures/tcga_fit_fingerprint.json").read_text())
    assert report["n_patients"] == 594
    assert report["views"]["profiling"]["n"] == 283
    assert report["views"]["transferable_model"]["n"] == fp["variant_2_transferable"]["n_fit"]
    assert report["views"]["molecular_model"]["n"] == fp["variant_3_molecular"]["n_fit"]
    flow = report["views"]["transferable_model"]["flow"]
    assert flow[-1]["n_remaining_train"] + flow[-1]["n_remaining_test"] == 570
    assert report["missing_reason_counts"]["kps"] == {
        "not_recorded_at_source": 149, "recorded_unknown_at_source": 3}
    assert report["missing_reason_counts"]["extent_of_resection"] == {
        "not_available_in_cohort": 594}
    # One cohort patient was sequenced only on a recurrent sample: tested, but
    # no admissible primary-tumour call.
    idh = report["counts"]["molecular"]["idh"]
    assert idh["call_status"] == {"called": 368, "not_tested": 225,
                                  "tested_not_callable": 1}
    assert idh["n_assayed_but_no_primary_call"] == 1
    assert report["counts"]["molecular"]["mgmt"]["call_status"] == {
        "called": 415, "not_tested": 179, "tested_not_callable": 0}


def _real_kwargs(tmp_path):
    from gliostrata.features.build import run_harmonize_clinical

    run_harmonize_clinical(str(REAL["cohort"]), str(REAL["raw"]),
                           str(tmp_path / "c.parquet"), str(tmp_path / "l.parquet"))
    return dict(
        clinical=str(tmp_path / "c.parquet"), labels=str(tmp_path / "l.parquet"),
        expression=str(REAL["tpm"]), annotation=str(REAL["annotation"]),
        mgmt_calls=str(REAL["mgmt"]), idh_calls=str(REAL["idh"]),
        mgmt_manifest=str(REAL["mgmt_manifest"]), idh_manifest=str(REAL["idh_manifest"]),
        centroids=str(REAL["centroids"]), partition=str(REAL["partition"]),
        out=str(tmp_path / "out" / "store.parquet"),
        report=str(tmp_path / "out" / "report.json"),
    )


@needs_real
def test_real_mode_refuses_the_centroid_set_reserved_for_t1(tmp_path):
    """The all-TCGA set includes the v2 test patients; the store must not use it."""
    from gliostrata.features.build import run_real

    kwargs = _real_kwargs(tmp_path)
    kwargs["centroids"] = str(ROOT / "configs/subtype_centroids/verhaak840_tcga_all.json")
    with pytest.raises(FeatureStoreError, match="reserved for T1"):
        run_real(**kwargs)
    assert not (tmp_path / "out").exists()


@needs_real
def test_real_mode_needs_a_matching_l1_provenance_sidecar(tmp_path):
    """The store's provenance reaches back to L0 through the L1 sidecar; a
    missing sidecar, or one that no longer matches the clinical table, stops it."""
    import json

    from gliostrata.features.build import l1_provenance_path, run_real

    kwargs = _real_kwargs(tmp_path)
    side = l1_provenance_path(kwargs["clinical"])
    doc = json.loads(side.read_text(encoding="utf-8"))
    assert set(doc["inputs"]) == {"cohort", "raw"}
    assert doc["inputs"]["cohort"]["sha256"] == sha256_file(REAL["cohort"])
    doc["outputs"]["clinical"]["sha256"] = "0" * 64
    side.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(FeatureStoreError, match="does not match its L1 provenance"):
        run_real(**kwargs)
    side.unlink()
    with pytest.raises(FeatureStoreError, match="no L1 provenance sidecar"):
        run_real(**kwargs)
    assert not (tmp_path / "out").exists()
