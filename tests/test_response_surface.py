"""Tests for the T4/T5 parameter loader, delivery module and response surface.

Nothing here reads a cohort. The committed parameter file is read as data, and the rest runs
on small synthetic inputs. Every test names the registered clause it holds the code to.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from gliostrata.prediction import response_surface as rs
from gliostrata.simulation import delivery
from gliostrata.simulation.parameters import (
    DEFAULT_PATH,
    PENDING,
    ParameterError,
    ParameterSet,
    load_parameters,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def committed() -> ParameterSet:
    """The real registered parameter file, as committed."""
    return load_parameters(ROOT / DEFAULT_PATH)


@pytest.fixture
def resolved(committed) -> ParameterSet:
    """The same file with the three pending Emax values filled, for arithmetic tests only.

    This is a test fixture, never a registered value: it exists so the surface can be
    exercised before the search of item 11(i.4) has entered its numbers.
    """
    import copy

    doc = copy.deepcopy(committed.doc)
    doc["frozen"] = True
    for payload, value in (("CRAd", 0.14), ("prodrug_enzyme", 0.10), ("CAR", 0.05)):
        doc["response"]["emax"][payload] = {"value": value, "status": "assumed"}
    return ParameterSet(doc=doc, path=committed.path, sha256="test-fixture")


@pytest.fixture
def reference() -> delivery.TropismReference:
    genes = ["CXCL12", "VEGFA", "PLAU", "CCL2", "TGFB1"]
    return delivery.TropismReference(
        mean={g: 5.0 for g in genes}, sd={g: 2.0 for g in genes}, source="TCGA log2 CPM, frozen")


# --- the parameter file itself ------------------------------------------------------------

def test_the_committed_file_loads_and_hashes(committed):
    assert len(committed.sha256) == 64
    assert committed.tau_months == 24  # item 11(b)
    assert committed.doc["survival_model"]["variant"] == 2  # item 11(g)


def test_item_11g_and_11l_keep_idh_and_mgmt_out_of_the_model(committed):
    covariates = committed.doc["survival_model"]["covariates"]
    assert covariates == ["age", "sex", "arm"]
    text = yaml.dump(committed.doc["response"]["patient_modifier"])
    assert "IDH" not in text and "MGMT" not in text.replace("MGMT (a mechanistic", "")


def test_the_draft_file_refuses_a_registered_read(committed):
    """Item 11(e) and (m): a T4 or T5 read waits for the freeze amendment."""
    assert committed.frozen is False
    with pytest.raises(ParameterError, match="frozen: false"):
        committed.require_frozen("T5")


def test_a_pending_emax_raises_rather_than_returning_zero(committed):
    """Whatever the committed file holds, a pending entry must never read as zero."""
    import copy

    doc = copy.deepcopy(committed.doc)
    doc["response"]["emax"]["CRAd"] = {"value": None, "status": PENDING, "open": True}
    waiting = ParameterSet(doc=doc, path=committed.path, sha256="test-fixture")
    assert waiting.pending() == ["CRAd"]
    with pytest.raises(ParameterError, match=PENDING):
        waiting.emax("CRAd")
    # The only zero available is the control's, and it has to be asked for by name.
    assert waiting.null_payload_emax() == 0.0


def test_the_committed_emax_values_are_entered_and_exact(committed):
    """Each value must equal factor * ln(ratio) exactly, not a rounded transcription."""
    assert committed.pending() == [], "a payload is still waiting on the search"
    factor = committed.translation_factor
    for payload in ("CRAd", "prodrug_enzyme", "CAR"):
        entry = committed.doc["response"]["emax"][payload]
        assert entry["status"] == "derived_by_rule"
        assert entry["n_studies"] >= 1
        assert committed.emax(payload) == pytest.approx(
            factor * math.log(entry["ratio_used"]), abs=1e-15)


def test_the_committed_values_match_the_result_file(committed):
    """The parameter file must carry what the computation produced, not a retyping of it."""
    result = json.loads(
        (ROOT / "docs/t4t5_emax_search/emax_result.json").read_text(encoding="utf-8"))
    for payload, info in result["payloads"].items():
        assert committed.emax(payload) == info["emax"], (
            f"{payload}: the parameter file and emax_result.json disagree")


def test_emax_refuses_combo_because_combo_is_assembled_per_carrier(committed):
    with pytest.raises(ParameterError, match="single payload"):
        committed.emax("combo")


# --- item 11(j): the per-carrier combo ---------------------------------------------------

def test_every_carrier_has_a_pair_it_can_actually_carry(committed):
    assert committed.combo_pair("iPSC-NSC") == ("CRAd", "prodrug_enzyme")
    assert committed.combo_pair("MSC") == ("prodrug_enzyme", "CAR")
    assert committed.combo_pair("HB1_F3_line") == ("CRAd", "prodrug_enzyme")


def test_a_combo_pair_that_breaks_item_4s_exclusions_is_refused(committed, tmp_path):
    """The contradiction item 11(j) exists to prevent: MSC cannot carry CRAd."""
    import copy

    doc = copy.deepcopy(committed.doc)
    doc["response"]["combo"]["pairs"]["MSC"] = ["CRAd", "prodrug_enzyme"]
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.dump(doc), encoding="utf-8")
    with pytest.raises(ParameterError, match="incompatible"):
        load_parameters(path)


def test_the_tropism_reference_may_not_be_a_cgga_statistic(committed, tmp_path):
    import copy

    doc = copy.deepcopy(committed.doc)
    doc["delivery"]["tropism_reference"]["source"] = "CGGA mRNAseq_693 cohort mean"
    path = tmp_path / "leak.yaml"
    path.write_text(yaml.dump(doc), encoding="utf-8")
    with pytest.raises(ParameterError, match="CGGA"):
        load_parameters(path)
    with pytest.raises(ParameterError, match="CGGA"):
        delivery.TropismReference(mean={"A": 1.0}, sd={"A": 1.0}, source="CGGA cohort")


def test_a_missing_parameter_file_is_an_error_not_a_default(tmp_path):
    with pytest.raises(ParameterError, match="no defaults"):
        load_parameters(tmp_path / "absent.yaml")


# --- delivery ----------------------------------------------------------------------------

def test_tropism_score_uses_the_frozen_reference_not_the_cohort(committed, reference):
    genes = list(committed.tropism_genes("iPSC-NSC"))
    # Two patients, both exactly at the reference mean: the score must be 0, whatever the
    # spread of the two columns is. A per-cohort z-score would give +/-1 instead.
    expr = pd.DataFrame(5.0, index=genes, columns=["p1", "p2"])
    expr.loc[genes[0], "p2"] = 5.0
    score = delivery.tropism_score(expr, "iPSC-NSC", committed, reference)
    assert np.allclose(score.to_numpy(), 0.0)
    expr.loc[genes[0], "p1"] = 7.0  # +1 SD on one of four genes
    score = delivery.tropism_score(expr, "iPSC-NSC", committed, reference)
    assert score["p1"] == pytest.approx(0.25)
    assert score["p2"] == pytest.approx(0.0)


def test_a_missing_gene_drops_for_that_patient_only(committed, reference):
    genes = list(committed.tropism_genes("iPSC-NSC"))
    expr = pd.DataFrame(5.0, index=genes, columns=["p1", "p2"])
    expr.loc[genes[0], "p1"] = np.nan
    score = delivery.tropism_score(expr, "iPSC-NSC", committed, reference)
    assert score.attrs["n_genes_used"] == len(genes)
    assert np.isfinite(score.to_numpy()).all()


def test_a_gene_absent_from_the_matrix_is_reported(committed, reference):
    genes = list(committed.tropism_genes("iPSC-NSC"))
    expr = pd.DataFrame(5.0, index=genes[1:], columns=["p1"])
    score = delivery.tropism_score(expr, "iPSC-NSC", committed, reference)
    assert score.attrs["n_genes_used"] == len(genes) - 1
    assert genes[0] not in score.attrs["genes_used"]
    with pytest.raises(ParameterError, match="none of"):
        delivery.tropism_score(
            pd.DataFrame(5.0, index=["NOTAGENE"], columns=["p1"]),
            "iPSC-NSC", committed, reference)


def test_effective_dose_is_monotone_in_tropism_and_capped_by_the_acting_fraction(committed):
    low = delivery.effective_dose(1e7, "iPSC-NSC", -2.0, committed)
    mid = delivery.effective_dose(1e7, "iPSC-NSC", 0.0, committed)
    high = delivery.effective_dose(1e7, "iPSC-NSC", 5.0, committed)
    assert low < mid < high
    # a_k = 0, so a patient at the reference mean gets exactly half the acting fraction
    assert mid == pytest.approx(1e7 * committed.cell_fraction_acting("iPSC-NSC") * 0.5)
    assert high < 1e7 * committed.cell_fraction_acting("iPSC-NSC")


def test_coverage_saturates_and_rejects_zero_sites(committed):
    cov = delivery.coverage(np.array([1, 2, 4, 8]), committed)
    assert np.all(np.diff(cov) > 0)
    assert cov[-1] < 1.0
    assert cov[0] == pytest.approx(1 - math.exp(-1 / committed.coverage_saturation_sites))
    with pytest.raises(ValueError, match="starts at 1"):
        delivery.coverage(0, committed)


def test_schedule_factors_are_ordered(committed):
    assert (committed.schedule_factor("single")
            < committed.schedule_factor("fractionated_x2")
            < committed.schedule_factor("fractionated_x3"))
    with pytest.raises(ParameterError):
        committed.schedule_factor("fractionated_x9")


# --- the response surface ----------------------------------------------------------------

def test_the_modifier_does_not_apply_where_the_mechanism_does_not(resolved):
    """The immune-TME term applies to CRAd and CAR, not to prodrug_enzyme."""
    mods = {"immune_tme_slope_per_sd": 2.0}
    assert rs.patient_modifier("CRAd", mods, resolved) > 1.0
    assert rs.patient_modifier("prodrug_enzyme", mods, resolved) == pytest.approx(1.0)
    with pytest.raises(ParameterError, match="no patient modifier"):
        rs.patient_modifier("CRAd", {"invented_slope": 1.0}, resolved)


def test_benefit_is_zero_at_zero_dose_and_rises_with_dose(resolved):
    args = ({}, 1.0, 1.0, resolved)
    zero, _ = rs.payload_benefit("CRAd", "iPSC-NSC", 0.0, *args)
    small, _ = rs.payload_benefit("CRAd", "iPSC-NSC", 1e6, *args)
    large, _ = rs.payload_benefit("CRAd", "iPSC-NSC", 1e9, *args)
    assert zero == 0.0
    assert 0 < small < large < resolved.emax("CRAd")


def test_at_ed50_the_occupancy_is_one_half(resolved):
    value, _ = rs.payload_benefit(
        "CRAd", "iPSC-NSC", resolved.ed50_cells, {}, 1.0, 1.0, resolved)
    assert value == pytest.approx(resolved.emax("CRAd") * 0.5)


def test_combo_is_bliss_independent_over_the_carriers_own_pair(resolved):
    n_eff = resolved.ed50_cells
    total, parts = rs.payload_benefit("combo", "MSC", n_eff, {}, 1.0, 1.0, resolved)
    assert parts["pair"] == ("prodrug_enzyme", "CAR")
    each = [rs.payload_benefit(p, "MSC", n_eff, {}, 1.0, 1.0, resolved)[0]
            for p in ("prodrug_enzyme", "CAR")]
    assert total == pytest.approx(sum(each) * resolved.synergy)
    # The same payload name gives a different combo on a different carrier (item 11(j)).
    other, other_parts = rs.payload_benefit("combo", "iPSC-NSC", n_eff, {}, 1.0, 1.0, resolved)
    assert other_parts["pair"] == ("CRAd", "prodrug_enzyme")
    assert other != total


def test_combo_toxicity_is_the_sum_of_the_carriers_pair(resolved):
    msc = resolved.tau_payload("combo", "MSC")
    assert msc == pytest.approx(resolved.tau_payload("prodrug_enzyme", "MSC")
                                + resolved.tau_payload("CAR", "MSC"))
    nsc = resolved.tau_payload("combo", "iPSC-NSC")
    assert nsc == pytest.approx(resolved.tau_payload("CRAd", "iPSC-NSC")
                                + resolved.tau_payload("prodrug_enzyme", "iPSC-NSC"))
    assert msc != nsc  # the pairs differ, so the toxicity does


def test_toxicity_rises_with_dose_and_sites_and_never_vanishes(resolved):
    base = rs.toxicity("CRAd", "iPSC-NSC", 1e6, 1, resolved)
    more_dose = rs.toxicity("CRAd", "iPSC-NSC", 5e7, 1, resolved)
    more_sites = rs.toxicity("CRAd", "iPSC-NSC", 1e6, 8, resolved)
    assert 0 < base < more_dose
    assert base < more_sites
    # Item 11(k): there is no configuration with zero harm, because of the delay cost.
    assert rs.toxicity("CRAd", "iPSC-NSC", 0.0, 1, resolved) == pytest.approx(
        resolved.delay_cost_log_hr)


def test_the_null_payload_control_as_corrected(resolved):
    """Item 11(o)(ii), as corrected: E is exactly 0 and the gain is negative, not 0."""
    effect = rs.configuration_effect(
        "CRAd", "iPSC-NSC", 2.5e7, 4, 1e7, {}, 1.0, 0.8, resolved, null_payload=True)
    assert effect.benefit == 0.0
    assert effect.toxicity > 0
    assert effect.log_hazard_offset == pytest.approx(effect.toxicity)
    grid = rs.time_grid(resolved)
    baseline = np.exp(-0.03 * grid)
    gain = rs.gain_months(baseline, effect, resolved)
    assert gain < 0  # not "exactly 0": the toxicity term stays in the primary surface
    # ...and with the toxicity term's own parameters at zero, the gain would be exactly 0.
    import copy
    doc = copy.deepcopy(resolved.doc)
    for payload in ("CRAd", "prodrug_enzyme", "CAR"):
        doc["toxicity"]["tau_payload"][payload] = {"value": 0.0}
    doc["toxicity"]["delay_cost_log_hr"] = {"value": 0.0}
    no_tox = ParameterSet(doc=doc, path=resolved.path, sha256="test-fixture")
    flat = rs.configuration_effect(
        "CRAd", "iPSC-NSC", 2.5e7, 4, 1e7, {}, 1.0, 0.8, no_tox, null_payload=True)
    assert rs.gain_months(baseline, flat, no_tox) == pytest.approx(0.0, abs=1e-12)


def test_the_no_modifier_control_makes_every_patient_identical(resolved):
    """Item 11(o)(i): with the modifiers at zero, the effect cannot vary by patient."""
    args = ("CRAd", "iPSC-NSC", 1e7)
    a, _ = rs.payload_benefit(*args, {"immune_tme_slope_per_sd": 0.0}, 1.0, 1.0, resolved)
    b, _ = rs.payload_benefit(*args, {"immune_tme_slope_per_sd": 0.0}, 1.0, 1.0, resolved)
    c, _ = rs.payload_benefit(*args, {}, 1.0, 1.0, resolved)
    assert a == b == c


def test_the_curve_transform_is_proportional_hazards(resolved):
    grid = rs.time_grid(resolved)
    baseline = np.exp(-0.05 * grid)
    better = rs.configuration_curve(baseline, -0.3)   # a benefit
    worse = rs.configuration_curve(baseline, +0.3)    # harm
    assert np.all(better[1:] > baseline[1:])
    assert np.all(worse[1:] < baseline[1:])
    assert rs.configuration_curve(baseline, 0.0) == pytest.approx(baseline)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        rs.configuration_curve(np.array([1.2]), 0.0)


def test_rmst_is_the_registered_grid_and_horizon(resolved):
    grid = rs.time_grid(resolved)
    assert grid[0] == 0.0 and grid[-1] == resolved.tau_months
    assert grid.size == int(round(resolved.tau_months / resolved.grid_step_months)) + 1
    # A patient who survives the whole horizon has RMST = tau.
    assert rs.rmst(np.ones_like(grid), resolved) == pytest.approx(resolved.tau_months)
    with pytest.raises(ValueError, match="registered grid"):
        rs.rmst(np.ones(5), resolved)


def test_gain_is_positive_when_the_benefit_outweighs_the_toxicity(resolved):
    grid = rs.time_grid(resolved)
    baseline = np.exp(-0.05 * grid)
    strong = rs.configuration_effect(
        "CRAd", "iPSC-NSC", 1e6, 1, 1e9, {}, 1.15, 0.95, resolved)
    assert strong.benefit > strong.toxicity
    assert rs.gain_months(baseline, strong, resolved) > 0


# --- the frozen tropism reference (item 11(h)) -------------------------------------------

def test_the_committed_tropism_reference_loads_and_is_tcga():
    reference = delivery.load_tropism_reference(
        ROOT / "configs" / "response_surface_v1_tropism_reference.json")
    assert "TCGA" in reference.source
    assert "cgga" not in reference.source.lower()
    assert set(reference.mean) == {"CCL2", "CXCL12", "PLAU", "TGFB1", "VEGFA"}
    assert all(sd > 0 for sd in reference.sd.values())


def test_the_reference_covers_every_carriers_genes(committed):
    reference = delivery.load_tropism_reference(
        ROOT / "configs" / "response_surface_v1_tropism_reference.json")
    for carrier in ("iPSC-NSC", "MSC", "HB1_F3_line"):
        for gene in committed.tropism_genes(carrier):
            assert gene in reference.mean, f"{carrier}'s gene {gene} has no reference"


def test_the_parameter_file_records_the_references_checksum(committed):
    import hashlib
    entry = committed.doc["delivery"]["tropism_reference"]
    path = ROOT / entry["file"]
    assert entry["present"] is True
    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"], (
        "the tropism reference has changed since its checksum was recorded")


def test_a_missing_reference_is_an_error_not_a_cohort_statistic(tmp_path):
    with pytest.raises(ParameterError, match="never be computed from the cohort"):
        delivery.load_tropism_reference(tmp_path / "absent.json")


def test_the_reference_also_covers_the_patient_modifiers():
    """Item 11: the modifier slopes are per SD, so u needs a FROZEN reference too.

    Only the tropism half existed until 2026-11-28. Without this half, the SD would have
    come from the cohort being scored, which is the leakage the reference exists to prevent.
    """
    doc = json.loads(
        (ROOT / "configs" / "response_surface_v1_tropism_reference.json").read_text(
            encoding="utf-8"))
    modifiers = doc.get("modifier_reference") or {}
    assert set(modifiers) == {"immune_tme_slope_per_sd", "mes_like_slope_per_sd"}
    for slope, entry in modifiers.items():
        assert entry["sd"] > 0, f"{slope} has a non-positive SD"
        assert entry["n"] > 0
    assert "tcga" in doc["modifier_source"].lower()
    assert "cgga" not in doc["modifier_source"].lower()


def test_every_modifier_the_surface_uses_has_a_frozen_reference(committed):
    """A modifier with a slope but no reference would be silently standardised on the cohort."""
    doc = json.loads(
        (ROOT / "configs" / "response_surface_v1_tropism_reference.json").read_text(
            encoding="utf-8"))
    referenced = set(doc.get("modifier_reference") or {})
    for name in committed.modifier_names:
        assert name in referenced, (
            f"{name} is used by the response surface but has no frozen reference")
