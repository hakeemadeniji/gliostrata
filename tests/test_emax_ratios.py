"""Tests for scripts/t4t5_emax_ratios.py, which computes the registered Emax.

Every test runs on synthetic extraction rows written to a temporary file. No article, no
cohort and no committed result file is read, and the real parameter file is read only for its
translation factor. Each test names the clause of pre-registration item 11(i) it holds the
code to.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "t4t5_emax_ratios", ROOT / "scripts" / "t4t5_emax_ratios.py")
ratios = importlib.util.module_from_spec(_spec)
sys.modules["t4t5_emax_ratios"] = ratios
_spec.loader.exec_module(ratios)


def _row(**kw) -> dict:
    """A row that qualifies under P1 to P5 and the (i.7) primary rule, unless overridden."""
    base = {
        "pmid": "1", "payload": "CRAd",
        "treated_median": 20.0, "control_median": 10.0, "unit": "days",
        "ratio": None, "arm_is_combination": False,
        "control_arm_kind": "untreated",
        "p2_intracranial": "yes", "p4_payload_match": "yes",
    }
    return {**base, **kw}


@pytest.fixture
def run(tmp_path, monkeypatch):
    """Write synthetic rows and compute, with the output paths redirected to tmp_path."""
    def _run(rows):
        final = tmp_path / "extraction_final.jsonl"
        final.write_text(
            "".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8", newline="\n")
        monkeypatch.setattr(ratios, "FINAL", final)
        monkeypatch.setattr(ratios, "PARAMS", ROOT / "configs" / "response_surface_v1.yaml")
        return ratios.compute()
    return _run


def test_the_translation_factor_comes_from_the_parameter_file(run):
    result = run([_row()])
    assert result["translation_factor"] == 0.5  # item 11(i.1), read not hard-coded
    assert "assumed" in result["translation_factor_status"]


def test_emax_is_the_registered_formula(run):
    result = run([_row(treated_median=20.0, control_median=10.0)])
    info = result["payloads"]["CRAd"]
    assert info["ratio"] == pytest.approx(2.0)
    assert info["emax"] == pytest.approx(0.5 * math.log(2.0))
    assert info["hazard_ratio"] == pytest.approx(math.exp(-0.5 * math.log(2.0)))


def test_the_ratio_is_a_median_not_a_mean(run):
    """Item 11(i.2): the median over qualifying studies. One outlier must not move it."""
    rows = [_row(pmid="a", treated_median=11.0, control_median=10.0),   # 1.1
            _row(pmid="b", treated_median=12.0, control_median=10.0),   # 1.2
            _row(pmid="c", treated_median=100.0, control_median=10.0)]  # 10.0
    info = run(rows)["payloads"]["CRAd"]
    assert info["ratio"] == pytest.approx(1.2)
    assert info["primary_narrowed"]["n"] == 3


def test_an_even_set_takes_the_midpoint(run):
    rows = [_row(pmid="a", treated_median=10.0, control_median=10.0),
            _row(pmid="b", treated_median=20.0, control_median=10.0)]
    assert run(rows)["payloads"]["CRAd"]["ratio"] == pytest.approx(1.5)


@pytest.mark.parametrize("kind", ["untreated", "vehicle", "mock_or_reporter_vector"])
def test_the_three_inert_control_kinds_enter_the_primary_set(run, kind):
    info = run([_row(control_arm_kind=kind)])["payloads"]["CRAd"]
    assert info["control_arm_rule"] == "narrowed"
    assert info["primary_narrowed"]["n"] == 1


@pytest.mark.parametrize("kind", ["active_therapeutic", "uncertain"])
def test_an_active_or_uncertain_control_is_kept_out_of_the_primary_set(run, kind):
    """Item 11(i.7): an uncertain control is excluded rather than argued about."""
    rows = [_row(pmid="a", control_arm_kind=kind, treated_median=30.0, control_median=10.0),
            _row(pmid="b", control_arm_kind="untreated", treated_median=12.0,
                 control_median=10.0)]
    info = run(rows)["payloads"]["CRAd"]
    assert info["control_arm_rule"] == "narrowed"
    assert info["primary_narrowed"]["n"] == 1
    assert info["ratio"] == pytest.approx(1.2)          # the excluded 3.0 does not count
    assert info["exploratory_literal"]["n"] == 2        # but it is still reported
    assert info["n_dropped_by_narrowing"] == 1


def test_the_narrowing_is_reported_beside_the_primary_not_instead_of_it(run):
    rows = [_row(pmid="a", control_arm_kind="untreated", treated_median=12.0,
                 control_median=10.0),
            _row(pmid="b", control_arm_kind="active_therapeutic", treated_median=30.0,
                 control_median=10.0)]
    info = run(rows)["payloads"]["CRAd"]
    assert info["exploratory_literal"]["ratio"] == pytest.approx((1.2 + 3.0) / 2)
    assert info["primary_narrowed"]["ratio"] == pytest.approx(1.2)


def test_an_empty_primary_set_falls_back_to_the_literal_set(run):
    """Item 11(i.7)'s guard: the deviation can never be worse-evidenced than (i.3)."""
    info = run([_row(control_arm_kind="uncertain", treated_median=20.0,
                     control_median=10.0)])["payloads"]["CRAd"]
    assert info["control_arm_rule"] == "literal_fallback"
    assert info["ratio"] == pytest.approx(2.0)
    assert info["emax"] == pytest.approx(0.5 * math.log(2.0))


def test_both_sets_empty_gives_an_unevidenced_payload(run):
    """Item 11(i.5): Emax = 0, and the consequence is that the optimiser never selects it."""
    info = run([_row(p2_intracranial="no")])["payloads"]["CRAd"]
    assert info["control_arm_rule"] == "unevidenced"
    assert info["emax"] == 0.0
    assert info["hazard_ratio"] == 1.0
    assert info["ratio"] is None


@pytest.mark.parametrize(("field", "value", "fragment"), [
    ("p2_intracranial", "no", "P2"),
    ("p2_intracranial", "uncertain", "P2"),
    ("p4_payload_match", "uncertain", "P4"),
    ("p4_payload_match", "no", "P4"),
])
def test_p2_and_p4_must_be_confirmed_not_merely_not_refuted(run, field, value, fragment):
    info = run([_row(**{field: value})])["payloads"]["CRAd"]
    assert info["n_eligible_p1_to_p5"] == 0
    assert any(fragment in why for why in info["exclusions"])


def test_a_combination_arm_is_not_the_payload_alone(run):
    info = run([_row(arm_is_combination=True)])["payloads"]["CRAd"]
    assert info["n_eligible_p1_to_p5"] == 0
    assert any("combination" in why for why in info["exclusions"])


def test_a_missing_median_is_excluded_under_p3(run):
    for row in (_row(treated_median=None), _row(control_median=None),
                _row(control_median=0.0)):
        info = run([row])["payloads"]["CRAd"]
        assert info["n_eligible_p1_to_p5"] == 0
        assert any("P3" in why for why in info["exclusions"])


def test_a_directly_reported_ratio_is_used_when_no_absolutes_exist(run):
    info = run([_row(treated_median=None, control_median=None, ratio=1.8)])["payloads"]["CRAd"]
    assert info["ratio"] == pytest.approx(1.8)


def test_payloads_are_computed_independently(run):
    rows = [_row(pmid="a", payload="CRAd", treated_median=12.0, control_median=10.0),
            _row(pmid="b", payload="CAR", treated_median=30.0, control_median=10.0)]
    payloads = run(rows)["payloads"]
    assert payloads["CRAd"]["ratio"] == pytest.approx(1.2)
    assert payloads["CAR"]["ratio"] == pytest.approx(3.0)
    assert payloads["prodrug_enzyme"]["control_arm_rule"] == "unevidenced"


def test_every_payload_appears_even_with_no_rows(run):
    payloads = run([])["payloads"]
    assert set(payloads) == {"CRAd", "prodrug_enzyme", "CAR"}
    assert all(p["emax"] == 0.0 for p in payloads.values())


def test_the_studies_used_are_named(run):
    rows = [_row(pmid="111", control_arm_kind="untreated"),
            _row(pmid="222", control_arm_kind="uncertain")]
    info = run(rows)["payloads"]["CRAd"]
    assert info["studies"] == ["111"]  # only what the primary set actually used


def test_computing_does_not_touch_the_parameter_file(run):
    params = ROOT / "configs" / "response_surface_v1.yaml"
    before = params.read_bytes()
    run([_row()])
    assert params.read_bytes() == before, (
        "entering a value in the parameter file is the freeze step of item 11(e) and (m); "
        "it must never be a side effect of computing a ratio")


def test_a_missing_extraction_file_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(ratios, "FINAL", tmp_path / "absent.jsonl")
    with pytest.raises(SystemExit, match="adjudicated union"):
        ratios.compute()
