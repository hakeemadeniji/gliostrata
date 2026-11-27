"""Tests for the held-out cohort read guard (gliostrata.evaluation.read_guard).

No test here reads CGGA. The repository ledger (data/access_ledger.jsonl) is
read as text; every access it describes is re-created on temporary files. The
one CGGA-shaped path that appears, `INPUT`, is written by the tests themselves
under `tmp_path`, with two synthetic rows.

What is at stake
----------------
Six things would each produce a plausible-looking run that the
pre-registration does not allow, and each has a test:

  1. A value read before the registration is deposited.
  2. A value read of a metric that was never registered.
  3. A value read with neither a sealed partition nor a registered scoring rule.
  4. A second value read of an analysis that has already been read.
  5. A value read whose scoring rule is free text, is not registered for the
     analysis, or names no frozen scoring-set manifest, an absent one or an
     altered one.
  6. A value read whose input checksums are not the ones on disk.

Further tests keep the ledger and the pre-registration from disagreeing about
how many times the external cohort has been read, and check that entries which
are not accesses are appended but never counted.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from gliostrata.cli import cli
from gliostrata.evaluation import freeze_partition, registered_metrics
from gliostrata.evaluation.read_guard import (
    KINDS,
    LEDGER_PATH,
    NON_ACCESS_KINDS,
    NOTE,
    OBSERVATION,
    REGISTRATION_PATH,
    SCORING_RULES_PATH,
    VALUES,
    Access,
    ReadRefusedError,
    append_ledger_entry,
    check_access,
    format_ledger,
    hash_inputs,
    is_deposited,
    load_ledger,
    load_scoring_rules,
    record_access,
    registration_state,
    summarize,
)

ROOT = Path(__file__).resolve().parents[1]
METRIC = "T4_optimizer_regret_reduction_pct"
MANIFEST = "data/partitions/ext_v1.json"
INPUT = "data/raw/cgga/clinical_cohort.parquet"
RULE = "fixture_rule"

# Synthetic: it exists only in the temporary files a test writes.
FIXTURE_DOI = "10.0000/osf.io/test-fixture"
NOT_DEPOSITED = dict.fromkeys(("osf_registration_doi", "deposited_at", "registered_commit"))


def _input(tmp_path):
    """Write the synthetic held-out input under `tmp_path`; return its checksum."""
    path = tmp_path / INPUT
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"patient_id": ["a", "b"], "value": [1.5, 2.5]}).to_parquet(path)
    return {INPUT: hashlib.sha256(path.read_bytes()).hexdigest()}


def _values(tmp_path, **overrides):
    """A value read that satisfies every gate once the fixtures are in place."""
    fields = dict(
        cohort="cgga",
        kind=VALUES,
        analysis=METRIC,
        purpose="fixture",
        inputs=_input(tmp_path),
        partition=MANIFEST,
    )
    fields.update(overrides)
    return Access(**fields)


def _deposited():
    return {**NOT_DEPOSITED, "osf_registration_doi": FIXTURE_DOI}


def _seal(tmp_path, *, seal=True):
    """Write a (by default sealed) all-test manifest over 40 synthetic ids."""
    cohort = pd.DataFrame({"patient_id": [f"P{i:03d}" for i in range(40)]})
    source = tmp_path / "cohort.parquet"
    cohort.to_parquet(source)
    path = tmp_path / MANIFEST
    freeze_partition(
        cohort,
        cohort_source=source,
        manifest_path=path,
        partition_name="ext_v1",
        seed=20251126,
        test_fraction=1.0,
        seal=seal,
    )
    return path


def _registry(tmp_path, body):
    """Write a scoring-rule registry under the fixture root."""
    path = tmp_path / SCORING_RULES_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(body), encoding="utf-8")
    return path


def _frozen_rule(tmp_path, *, seal=True):
    """A registry whose one rule names the fixture manifest and its true checksum."""
    path = _seal(tmp_path, seal=seal)
    return {
        "schema_version": 1,
        "rules": {
            RULE: {
                "analyses": [METRIC],
                "preregistration": "fixture section",
                "scoring_set": "fixture",
                "scoring_set_manifest": MANIFEST,
                "scoring_set_manifest_sha256": hashlib.sha256(
                    path.read_bytes()
                ).hexdigest(),
            }
        },
    }


def _by_rule(tmp_path, rule=RULE, analysis=METRIC):
    return _values(tmp_path, analysis=analysis, partition=None, scoring_rule=rule)


# ---------------------------------------------------------------------------
# The gates
# ---------------------------------------------------------------------------
def test_only_a_value_access_is_gated(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    for kind in KINDS:
        if kind != VALUES:
            record_access(
                Access("cgga", kind, "look", "counts only"),
                ledger_path=ledger,
                registration_path=tmp_path / "missing.yaml",
                root=tmp_path,
            )
    summary = summarize(load_ledger(ledger))
    assert summary["entries"] == len(KINDS) - 1
    assert summary["value_reads"] == 0
    assert summary["before_deposit"] == len(KINDS) - 1


def test_a_value_read_is_refused_before_the_deposit(tmp_path):
    _seal(tmp_path)
    reasons = check_access(
        _values(tmp_path), ledger=[], registration=NOT_DEPOSITED, root=tmp_path
    )
    assert len(reasons) == 1
    assert "REGISTRATION GATE" in reasons[0]


def test_placeholder_text_does_not_count_as_a_deposit():
    assert not is_deposited({**NOT_DEPOSITED, "osf_registration_doi": "TODO"})
    assert not is_deposited({**NOT_DEPOSITED, "osf_registration_doi": "pending"})
    assert not is_deposited({**NOT_DEPOSITED, "osf_registration_doi": 20261126})
    assert is_deposited(_deposited())


def test_an_unregistered_metric_is_refused(tmp_path):
    _seal(tmp_path)
    reasons = check_access(
        _values(tmp_path, analysis="T9_something_new"),
        ledger=[],
        registration=_deposited(),
        root=tmp_path,
    )
    assert len(reasons) == 1
    assert "not one of the registered metrics" in reasons[0]
    assert "T9_something_new" in reasons[0]


@pytest.mark.parametrize("metric", registered_metrics())
def test_every_registered_metric_passes_the_registration_gate(tmp_path, metric):
    _seal(tmp_path)
    assert (
        check_access(
            _values(tmp_path, analysis=metric),
            ledger=[],
            registration=_deposited(),
            root=tmp_path,
        )
        == []
    )


def test_an_exploratory_read_needs_its_own_preregistration(tmp_path):
    _seal(tmp_path)
    access = _values(tmp_path, analysis="exploratory:calibration-slope")
    reasons = check_access(access, ledger=[], registration=_deposited(), root=tmp_path)
    assert len(reasons) == 1
    assert "exploratory" in reasons[0]
    cited = _values(
        tmp_path,
        analysis="exploratory:calibration-slope",
        exploratory_preregistration="osf.io/example (a citation that exists only in this test)",
    )
    assert check_access(cited, ledger=[], registration=_deposited(), root=tmp_path) == []


def test_partition_gate(tmp_path):
    missing = check_access(
        _values(tmp_path), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert [r for r in missing if "PARTITION GATE" in r and "does not exist" in r]

    _seal(tmp_path, seal=False)
    draft = check_access(_values(tmp_path), ledger=[], registration=_deposited(), root=tmp_path)
    assert [r for r in draft if "PARTITION GATE" in r and "unsealed draft" in r]

    neither = check_access(
        _values(tmp_path, partition=None), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert [r for r in neither if "PARTITION GATE" in r]


def test_a_scoring_rule_in_free_text_no_longer_passes(tmp_path):
    """Until 2026-11-27 any non-empty text opened this gate."""
    _registry(tmp_path, _frozen_rule(tmp_path))
    reasons = check_access(
        _by_rule(tmp_path, rule="docs/preregistration.md section 2.2 item 9"),
        ledger=[],
        registration=_deposited(),
        root=tmp_path,
    )
    assert len(reasons) == 1
    assert "SCORING RULE GATE" in reasons[0] and "not a rule id" in reasons[0]


def test_a_registered_rule_with_its_frozen_manifest_passes(tmp_path):
    _registry(tmp_path, _frozen_rule(tmp_path))
    access = _by_rule(tmp_path)
    assert check_access(access, ledger=[], registration=_deposited(), root=tmp_path) == []

    config = tmp_path / "configs" / "registration.yaml"
    config.parent.mkdir(exist_ok=True)
    config.write_text(f"osf_registration_doi: '{FIXTURE_DOI}'\n", encoding="utf-8")
    entry = record_access(
        access,
        ledger_path=tmp_path / "ledger.jsonl",
        registration_path=config,
        root=tmp_path,
    )
    assert entry["scoring_rule"] == RULE
    assert entry["partition"] is None


def test_an_unknown_scoring_rule_is_refused(tmp_path):
    _registry(tmp_path, _frozen_rule(tmp_path))
    reasons = check_access(
        _by_rule(tmp_path, rule="T9_not_a_rule"),
        ledger=[],
        registration=_deposited(),
        root=tmp_path,
    )
    assert len(reasons) == 1
    assert "SCORING RULE GATE" in reasons[0] and "not a rule id" in reasons[0]

    # With no registry at all the path is closed too, and says why.
    bare = tmp_path / "bare"
    bare.mkdir()
    reasons = check_access(
        _by_rule(bare), ledger=[], registration=_deposited(), root=bare
    )
    assert [r for r in reasons if "SCORING RULE GATE" in r and "does not exist" in r]


def test_a_rule_without_a_frozen_manifest_is_refused(tmp_path):
    body = _frozen_rule(tmp_path)
    body["rules"][RULE]["scoring_set_manifest"] = None
    body["rules"][RULE]["scoring_set_manifest_sha256"] = None
    _registry(tmp_path, body)
    reasons = check_access(
        _by_rule(tmp_path), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "names no frozen scoring-set manifest" in reasons[0]

    # A manifest named but not on disk is the same refusal, stated.
    body["rules"][RULE]["scoring_set_manifest"] = "data/partitions/never_frozen.json"
    body["rules"][RULE]["scoring_set_manifest_sha256"] = "0" * 64
    _registry(tmp_path, body)
    reasons = check_access(
        _by_rule(tmp_path), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "is missing" in reasons[0]


def test_an_altered_scoring_set_manifest_is_refused(tmp_path):
    _registry(tmp_path, _frozen_rule(tmp_path))
    path = tmp_path / MANIFEST
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["sealed"] = False
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    reasons = check_access(
        _by_rule(tmp_path), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "does not hash to the checksum" in reasons[0] and "altered" in reasons[0]

    # Recording the new checksum does not rescue an unsealed draft.
    _registry(tmp_path, _frozen_rule(tmp_path, seal=False))
    reasons = check_access(
        _by_rule(tmp_path), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "unsealed draft" in reasons[0]


def test_a_rule_registered_for_another_analysis_is_refused(tmp_path):
    _registry(tmp_path, _frozen_rule(tmp_path))
    reasons = check_access(
        _by_rule(tmp_path, analysis="T5_endtoend_rmst_success_fraction"),
        ledger=[],
        registration=_deposited(),
        root=tmp_path,
    )
    assert len(reasons) == 1
    assert f"fixes the scoring set of {METRIC}, not of" in reasons[0]


def test_a_malformed_registry_closes_the_path(tmp_path):
    _seal(tmp_path)
    for text in (
        "",
        "rules: {}\n",
        "rules:\n  fixture_rule:\n    analyses: [a]\n"
        "    scoring_set_manifest: data/partitions/ext_v1.json\n",
    ):
        path = tmp_path / SCORING_RULES_PATH
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        with pytest.raises(ValueError):
            load_scoring_rules(path)
        reasons = check_access(
            _by_rule(tmp_path), ledger=[], registration=_deposited(), root=tmp_path
        )
        assert len(reasons) == 1 and "SCORING RULE GATE" in reasons[0]


def test_second_value_read_needs_a_repeat_reason(tmp_path):
    _seal(tmp_path)
    first = {"cohort": "cgga", "analysis": METRIC, "kind": VALUES, "date": "2026-11-26"}
    reasons = check_access(
        _values(tmp_path), ledger=[first], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "ONCE PER ANALYSIS" in reasons[0] and "2026-11-26" in reasons[0]

    again = _values(tmp_path, repeat_reason="disclosed re-run after a documented amendment")
    assert check_access(again, ledger=[first], registration=_deposited(), root=tmp_path) == []

    # Another analysis is not spent by this read, and neither is a dry run.
    dry_run = {"cohort": "cgga", "analysis": METRIC, "kind": "availability"}
    other = _values(tmp_path, analysis="T5_endtoend_rmst_success_fraction")
    assert check_access(other, ledger=[dry_run], registration=_deposited(), root=tmp_path) == []


def test_value_read_without_input_checksums_is_refused(tmp_path):
    _seal(tmp_path)
    reasons = check_access(
        _values(tmp_path, inputs={}), ledger=[], registration=_deposited(), root=tmp_path
    )
    assert len(reasons) == 1
    assert "SHA-256" in reasons[0]


def test_input_checksums_are_recomputed_at_read_time(tmp_path):
    """The ledger must record what is on disk when the read happens."""
    _seal(tmp_path)
    access = _values(tmp_path)
    gates = dict(ledger=[], registration=_deposited(), root=tmp_path)
    assert check_access(access, **gates) == []

    path = tmp_path / INPUT
    path.write_bytes(path.read_bytes() + b"altered")
    reasons = check_access(access, **gates)
    assert len(reasons) == 1
    assert "INPUT GATE" in reasons[0] and "no longer hashes" in reasons[0]

    path.unlink()
    reasons = check_access(access, **gates)
    assert len(reasons) == 1
    assert "INPUT GATE" in reasons[0] and "cannot be recomputed" in reasons[0]

    _input(tmp_path)  # restored: the recorded checksum is the file's again
    assert check_access(access, **gates) == []


def test_all_reasons_are_reported_at_once(tmp_path):
    reasons = check_access(
        _values(tmp_path, partition=None), ledger=[], registration=NOT_DEPOSITED, root=tmp_path
    )
    joined = " | ".join(reasons)
    assert "REGISTRATION GATE" in joined
    assert "PARTITION GATE" in joined


def test_empty_purpose_unknown_cohort_or_unknown_kind_is_refused():
    assert check_access(Access("cgga", "schema", "x", " "))
    assert check_access(Access("cgga", "schema", " ", "p"))
    assert check_access(Access("tcga", "schema", "x", "development data"))
    assert check_access(Access("cgga", "guess", "x", "p"))
    assert check_access(Access("cgga", "schema", "x", "counts only")) == []


def test_the_committed_scoring_rule_registry_is_derived_from_the_preregistration():
    rules = load_scoring_rules(ROOT / SCORING_RULES_PATH)
    metrics = set(registered_metrics())
    assert rules, "an empty registry would close the scoring-rule path silently"
    for rule in rules.values():
        assert set(rule.analyses) <= metrics
        assert rule.preregistration.startswith("docs/preregistration.md")
        assert "section" in rule.preregistration
        assert rule.scoring_set
    # Since 2026-11-27 T1's scoring set is frozen (cgga_t1_subtype_v1). Every frozen
    # manifest must be committed and still hash to the checksum recorded beside it.
    for rule in rules.values():
        if rule.manifest is not None:
            path = ROOT / rule.manifest
            assert path.is_file(), f"{rule.manifest} is named but not committed"
            assert hashlib.sha256(path.read_bytes()).hexdigest() == rule.manifest_sha256


# ---------------------------------------------------------------------------
# Recording
# ---------------------------------------------------------------------------
def test_record_access_appends_and_a_refused_read_is_not_written(tmp_path):
    _seal(tmp_path)
    config = tmp_path / "configs" / "registration.yaml"
    config.parent.mkdir()
    config.write_text(f"osf_registration_doi: '{FIXTURE_DOI}'\n", encoding="utf-8")
    ledger = tmp_path / "ledger.jsonl"
    kwargs = dict(ledger_path=ledger, registration_path=config, root=tmp_path)

    entry = record_access(_values(tmp_path), **kwargs)
    assert entry["seq"] == 0
    assert entry["backfilled"] is False
    assert entry["before_deposit"] is False
    assert "T" in entry["recorded_at"]  # a time of day, not only a date

    second = record_access(
        _values(tmp_path, analysis="T5_endtoend_rmst_success_fraction"), **kwargs
    )
    assert second["seq"] == 1

    with pytest.raises(ReadRefusedError, match="ONCE PER ANALYSIS"):
        record_access(_values(tmp_path), **kwargs)
    assert len(load_ledger(ledger)) == 2


def test_a_value_read_is_refused_while_not_deposited(tmp_path):
    _seal(tmp_path)
    ledger = tmp_path / "ledger.jsonl"
    with pytest.raises(ReadRefusedError, match="REGISTRATION GATE"):
        record_access(
            _values(tmp_path),
            ledger_path=ledger,
            registration_path=tmp_path / "absent.yaml",
            root=tmp_path,
        )
    assert not ledger.exists()


def test_hash_inputs_hashes_bytes_only(tmp_path):
    path = tmp_path / "table.parquet"
    pd.DataFrame({"patient_id": ["a", "b"], "value": [1.5, 2.5]}).to_parquet(path)
    expected = hashlib.sha256(path.read_bytes()).hexdigest()
    assert hash_inputs([path]) == {path.as_posix(): expected}


def test_a_ledger_line_that_was_edited_is_a_finding(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    ledger.write_text('{"seq": 0}\nnot json\n', encoding="utf-8")
    with pytest.raises(ValueError, match="append-only"):
        load_ledger(ledger)
    assert load_ledger(tmp_path / "absent.jsonl") == []


def test_observations_and_notes_are_appended_and_never_count_as_reads(tmp_path):
    """The ledger is append-only: an observation and its follow-up are entries."""
    ledger = tmp_path / "ledger.jsonl"
    absent = tmp_path / "absent.yaml"
    seen = append_ledger_entry(
        OBSERVATION,
        "unexplained-file-access",
        "Seen in file-system metadata; process not identified.",
        details={"files": [{"path": INPUT, "bytes": 25_020}]},
        ledger_path=ledger,
        registration_path=absent,
    )
    closed = append_ledger_entry(
        NOTE,
        "unexplained-file-access",
        "Identified: a workspace-wide text search.",
        ledger_path=ledger,
        registration_path=absent,
    )
    assert (seen["seq"], closed["seq"]) == (0, 1)
    assert len(seen["recorded_at"]) > len(seen["date"]) and "T" in seen["recorded_at"]
    assert seen["before_deposit"] is True

    summary = summarize(load_ledger(ledger))
    assert summary["entries"] == 2
    assert summary["value_reads"] == 0
    assert summary["non_access_entries"] == 2
    assert summary["by_kind"] == {kind: 0 for kind in KINDS}

    # Neither kind is a declared access, so the gates refuse a script that claims it.
    assert check_access(Access("cgga", OBSERVATION, "x", "p"))
    assert check_access(Access("cgga", NOTE, "x", "p"))
    with pytest.raises(ValueError, match="not one of"):
        append_ledger_entry(VALUES, "x", "p", ledger_path=ledger, registration_path=absent)
    with pytest.raises(ValueError, match="required"):
        append_ledger_entry(NOTE, " ", "p", ledger_path=ledger, registration_path=absent)
    assert len(load_ledger(ledger)) == 2

    text = format_ledger(load_ledger(ledger))
    assert text.count("not an access") == 2
    assert "details:" in text


def test_format_ledger_numbers_value_reads_in_order():
    kinds = ("schema", VALUES, "availability", VALUES)
    entries = [
        {"seq": i, "date": "2026-11-26", "kind": k, "analysis": "a", "purpose": "p"}
        for i, k in enumerate(kinds)
    ]
    text = format_ledger(entries)
    assert "value read 1" in text
    assert "value read 2" in text
    assert "Value reads of the external cohort: 2 (a, a)." in text
    assert "empty" in format_ledger([])


# ---------------------------------------------------------------------------
# The repository's own ledger and registration state
# ---------------------------------------------------------------------------
def test_the_repository_is_not_registered_yet():
    state = registration_state(ROOT / REGISTRATION_PATH)
    assert state["config_found"] is True
    assert not is_deposited(state)


def test_the_repository_ledger_agrees_with_the_registration():
    entries = load_ledger(ROOT / LEDGER_PATH)
    assert [entry["seq"] for entry in entries] == list(range(len(entries)))
    assert all(
        entry["cohort"] == "cgga" and entry["kind"] in KINDS + NON_ACCESS_KINDS
        for entry in entries
    )
    assert all(entry["before_deposit"] for entry in entries)

    reads = summarize(entries)
    assert reads["value_reads"] == 2
    assert reads["value_read_analyses"] == ["T2_survival_cindex_cgga", "T2_brier_18mo"]

    cohort_file = "data/raw/cgga/clinical_cohort.parquet"
    digests = {entry["inputs"][cohort_file] for entry in entries if entry["kind"] == VALUES}
    assert len(digests) == 1
    assert len(next(iter(digests))) == 64

    preregistration = (ROOT / "docs" / "preregistration.md").read_text(encoding="utf-8")
    assert "CGGA has been read twice, both for registered rows of T2" in preregistration


def test_every_ledger_entry_from_seven_on_carries_a_time_of_day():
    """Entries 0 to 6 were written with a date only; a note says so."""
    entries = load_ledger(ROOT / LEDGER_PATH)
    dated_only = [entry["seq"] for entry in entries if "T" not in str(entry["recorded_at"])]
    assert dated_only == [0, 1, 2, 3, 4, 5, 6]
    assert all("T" in str(entry["recorded_at"]) for entry in entries if entry["seq"] >= 7)
    notes = [e for e in entries if e["kind"] == NOTE and "time of day" in e["purpose"].lower()]
    assert notes
    assert notes[0]["purpose"].lower().count("time of day not recorded") == 1


def test_registration_config_and_preregistration_header_agree():
    """Deposit day: this fails until the header and the config are both updated."""
    header = (ROOT / "docs" / "preregistration.md").read_text(encoding="utf-8")[:4000]
    doi = registration_state(ROOT / REGISTRATION_PATH)["osf_registration_doi"]
    if doi:
        assert doi in header, "the DOI is in the config but not in the header"
    else:
        assert "TO BE COMPLETED BY AUTHOR: OSF DOI" in header
        assert "not yet deposited" in header


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def test_cli_prints_the_ledger_and_the_registration_state(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    record_access(
        Access("cgga", "schema", "look", "counts only"),
        ledger_path=ledger,
        registration_path=tmp_path / "absent.yaml",
        root=tmp_path,
    )
    result = CliRunner().invoke(cli, ["access-ledger", "--ledger", str(ledger)])
    assert result.exit_code == 0, result.output
    assert "1 entries" in result.output
    assert "#0" in result.output
    assert "counts only" in result.output
    assert "Registration:" in result.output


def test_cli_prints_an_empty_ledger(tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    ledger.write_text("", encoding="utf-8")
    result = CliRunner().invoke(cli, ["access-ledger", "--ledger", str(ledger)])
    assert result.exit_code == 0, result.output
    assert "never been accessed" in result.output

