"""Read discipline for the held-out cohort.

The pre-registration (``docs/preregistration.md``: section 2.2 item 7, section
3.3 and section 4) allows the external cohort to be read once per registered
analysis, on a sealed partition, once the registration is deposited. This
module turns that procedure into three things a script cannot skip by accident:

* **A ledger.** Every access to a held-out cohort is appended, as one JSON
  line, to ``data/access_ledger.jsonl``. Dry runs and checksum-only accesses
  are recorded too. The file is version-controlled and is never edited.
* **Gates.** An access that inspects values -- any feature or outcome -- is
  refused unless the registration is deposited, the analysis is a registered
  metric or carries its own pre-registration, the partition is sealed or the
  scoring set is fixed by a registered rule, this analysis has not been read
  before, and every input's checksum is still the one given for it.
* **Checksums.** A value read records the SHA-256 of every file it opens, and
  the guard recomputes each one at read time.

This is a discipline, not a sandbox: a script has to call :func:`record_access`
before it opens the data. The two T2 scripts predate it and do not; their
reads are back-filled in the ledger and marked as made before deposit.

The scoring-set rule is an id, not a citation. ``scoring_rule`` must be a rule
id from the committed registry (``configs/scoring_rules.yaml``, derived from
the pre-registration), the rule must be registered for the analysis, it must
name a frozen scoring-set manifest, and that manifest's bytes must still hash
to the checksum recorded beside it. A rule whose manifest is not frozen yet
refuses every read.

The ledger also carries lines that are not access declarations: an
``observation`` (an access seen in file-system metadata, made outside this
module) and a ``note`` (a dated statement about earlier entries). Both are
appended with :func:`append_ledger_entry`; neither passes through the gates
and neither counts as a read.

Typical use, in a script that evaluates a registered test on the external
cohort::

    from gliostrata.evaluation.read_guard import Access, hash_inputs, record_access

    record_access(Access(
        cohort="cgga",
        kind="values",
        analysis="T1_subtype_accuracy_cgga",
        purpose="T1 subtype accuracy on the labelled microarray patients",
        script="scripts/validate_subtype_cgga.py",
        inputs=hash_inputs(["data/raw/cgga/clinical_cohort.parquet"]),
        scoring_rule="T1_subtype_accuracy_cgga",
    ))
"""
from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from gliostrata.evaluation.partitions import file_content_sha256, read_manifest
from gliostrata.evaluation.thresholds import UnregisteredMetricError, get_threshold

LEDGER_PATH = Path("data/access_ledger.jsonl")
REGISTRATION_PATH = Path("configs/registration.yaml")
SCORING_RULES_PATH = Path("configs/scoring_rules.yaml")
LEDGER_SCHEMA_VERSION = 1

#: Cohorts whose data are held out. Everything else is development data.
HELD_OUT_COHORTS: tuple[str, ...] = ("cgga",)

#: What an access inspects. Only ``VALUES`` is gated.
CHECKSUM = "checksum"          # file bytes hashed; nothing deserialised
SCHEMA = "schema"              # names, dtypes, null counts, marginal distributions
AVAILABILITY = "availability"  # which patients are evaluable; counts only
VALUES = "values"              # a model, a metric, or any feature against an outcome
KINDS = (CHECKSUM, SCHEMA, AVAILABILITY, VALUES)

#: Ledger lines that are not access declarations. They are appended with
#: :func:`append_ledger_entry`, never through the gates, and never count as a read.
OBSERVATION = "observation"  # an access seen in file-system metadata, made outside this module
NOTE = "note"                # a dated statement about earlier entries
NON_ACCESS_KINDS = (OBSERVATION, NOTE)

EXPLORATORY_PREFIX = "exploratory:"
_REGISTRATION_FIELDS = ("osf_registration_doi", "deposited_at", "registered_commit")
_REQUIRED_TEXT = ("analysis", "purpose")

# A DOI, and only a DOI, opens the registration gate. "TODO" or "pending" do not.
_DOI_PATTERN = re.compile(r"^10\.\d{4,9}/\S+$")


class ReadRefusedError(RuntimeError):
    """An access to the held-out cohort that the gates do not allow."""


@dataclass(frozen=True)
class Access:
    """One access to a held-out cohort, declared before the data are opened.

    ``kind`` is one of :data:`KINDS`. For a ``values`` access, ``partition``
    names the sealed manifest the analysis is scored on, or ``scoring_rule``
    gives the id of the registered rule that fixes the scoring set instead
    (T1's accuracy row is fixed by rule; T2, T4 and T5 by the sealed external
    partition). The id must be in ``configs/scoring_rules.yaml``.
    ``inputs`` maps each file's path to its SHA-256; the guard recomputes it.
    """

    cohort: str
    kind: str
    analysis: str
    purpose: str
    script: str = ""
    inputs: dict[str, str] = field(default_factory=dict)
    partition: str | None = None
    scoring_rule: str | None = None
    exploratory_preregistration: str | None = None
    repeat_reason: str | None = None


def hash_inputs(paths: Iterable[str | Path]) -> dict[str, str]:
    """SHA-256 of each file's bytes, keyed by its forward-slashed path.

    Reads bytes only: no table is deserialised and no value is interpreted.
    """
    hashed: dict[str, str] = {}
    for item in paths:
        path = Path(item)
        hashed[path.as_posix()] = file_content_sha256(path)
    return hashed


def registration_state(path: str | Path = REGISTRATION_PATH) -> dict[str, Any]:
    """The machine-readable registration state (``configs/registration.yaml``).

    A missing file reads as "not deposited", the safe default: a deleted
    config can never open a gate.
    """
    config = Path(path)
    state: dict[str, Any] = dict.fromkeys(_REGISTRATION_FIELDS)
    state["config_found"] = config.is_file()
    if state["config_found"]:
        data = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
        if not isinstance(data, dict):
            raise ValueError(f"{config}: expected a mapping of registration fields")
        for key in _REGISTRATION_FIELDS:
            state[key] = data.get(key)
    return state


def is_deposited(state: Mapping[str, Any]) -> bool:
    """True only when the registration state carries a DOI-shaped identifier."""
    doi = state.get("osf_registration_doi")
    return isinstance(doi, str) and bool(_DOI_PATTERN.match(doi.strip()))


@dataclass(frozen=True)
class ScoringRule:
    """One rule of the committed registry (``configs/scoring_rules.yaml``).

    ``analyses`` are the registered metrics the rule fixes a scoring set for.
    ``manifest`` and ``manifest_sha256`` are ``None`` until the scoring-set
    manifest is frozen; a rule in that state refuses every read.
    """

    rule_id: str
    analyses: tuple[str, ...]
    preregistration: str
    scoring_set: str
    manifest: str | None
    manifest_sha256: str | None


def load_scoring_rules(path: str | Path = SCORING_RULES_PATH) -> dict[str, ScoringRule]:
    """The committed registry of rules that fix a scoring set, keyed by rule id.

    Raises
    ------
    ValueError
        If the registry is missing or malformed. A registry that cannot be read
        closes the scoring-rule path; it never opens it.
    """
    registry_path = Path(path)
    if not registry_path.is_file():
        raise ValueError(f"{registry_path}: the scoring-rule registry does not exist")
    data = yaml.safe_load(registry_path.read_text(encoding="utf-8")) or {}
    rules = data.get("rules") if isinstance(data, dict) else None
    if not isinstance(rules, dict) or not rules:
        raise ValueError(f"{registry_path}: expected a non-empty 'rules' mapping")

    loaded: dict[str, ScoringRule] = {}
    for rule_id, body in rules.items():
        if not isinstance(body, dict):
            raise ValueError(f"{registry_path}: rule {rule_id!r} is not a mapping")
        manifest = body.get("scoring_set_manifest") or None
        checksum = body.get("scoring_set_manifest_sha256") or None
        if (manifest is None) != (checksum is None):
            raise ValueError(
                f"{registry_path}: rule {rule_id!r} needs both a manifest and its "
                "checksum, or neither"
            )
        analyses = tuple(str(a) for a in (body.get("analyses") or []))
        if not analyses:
            raise ValueError(f"{registry_path}: rule {rule_id!r} lists no analysis")
        loaded[str(rule_id)] = ScoringRule(
            rule_id=str(rule_id),
            analyses=analyses,
            preregistration=str(body.get("preregistration", "")),
            scoring_set=str(body.get("scoring_set", "")),
            manifest=manifest,
            manifest_sha256=checksum,
        )
    return loaded


def load_ledger(path: str | Path = LEDGER_PATH) -> list[dict[str, Any]]:
    """Every ledger entry, oldest first. An absent ledger is an empty list."""
    ledger = Path(path)
    if not ledger.is_file():
        return []
    entries: list[dict[str, Any]] = []
    for number, line in enumerate(ledger.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            entries.append(json.loads(line))
        except ValueError as exc:
            raise ValueError(
                f"{ledger}: line {number} is not JSON ({exc}). The ledger is "
                "append-only; an edited line is a finding, not something to repair."
            ) from exc
    return entries


def _scoring_rule_reasons(analysis: str, rule_id: str, root: str | Path) -> list[str]:
    """Why a value read that cites a scoring-set rule is refused.

    A citation in free text fixes nothing. The rule must be in the committed
    registry (``configs/scoring_rules.yaml``, derived from the
    pre-registration), it must be the rule registered for this analysis, it
    must name a frozen scoring-set manifest, and that manifest's bytes must
    still hash to the checksum recorded beside it.
    """
    rule_id = str(rule_id).strip()
    try:
        rules = load_scoring_rules(Path(root) / SCORING_RULES_PATH)
    except ValueError as exc:
        return [f"SCORING RULE GATE: {exc}. No scoring-set rule can be verified."]
    if rule_id not in rules:
        return [
            f"SCORING RULE GATE: {rule_id!r} is not a rule id in "
            f"{SCORING_RULES_PATH.as_posix()}. A rule id from that registry, not a "
            "citation, fixes the scoring set."
        ]
    rule = rules[rule_id]
    if analysis not in rule.analyses:
        return [
            f"SCORING RULE GATE: rule {rule_id} fixes the scoring set of "
            f"{', '.join(rule.analyses)}, not of {analysis}."
        ]
    if rule.manifest is None or rule.manifest_sha256 is None:
        return [
            f"SCORING RULE GATE: rule {rule_id} names no frozen scoring-set manifest "
            f"yet ({SCORING_RULES_PATH.as_posix()}: scoring_set_manifest is empty)."
        ]
    manifest_path = Path(root) / rule.manifest
    if not manifest_path.is_file():
        return [
            f"SCORING RULE GATE: the scoring-set manifest of rule {rule_id} is "
            f"missing ({rule.manifest})."
        ]
    actual = file_content_sha256(manifest_path)
    if actual != rule.manifest_sha256.strip().lower():
        return [
            f"SCORING RULE GATE: scoring-set manifest {rule.manifest} does not hash "
            f"to the checksum recorded for rule {rule_id} (recorded "
            f"{rule.manifest_sha256[:12]}..., recomputed {actual[:12]}...). It was "
            "altered after it was frozen."
        ]
    if not read_manifest(manifest_path).sealed:
        return [
            f"SCORING RULE GATE: scoring-set manifest {rule.manifest} is an "
            "unsealed draft."
        ]
    return []


def _input_reasons(inputs: Mapping[str, str], root: str | Path) -> list[str]:
    """Recompute every input's SHA-256: the ledger records what is on disk now."""
    reasons: list[str] = []
    for path, recorded in inputs.items():
        target = Path(root) / path
        if not target.is_file():
            reasons.append(
                f"INPUT GATE: {path} does not exist under {root}, so its SHA-256 "
                "cannot be recomputed."
            )
            continue
        actual = file_content_sha256(target)
        if actual != str(recorded).strip().lower():
            reasons.append(
                f"INPUT GATE: {path} no longer hashes to the SHA-256 given for it "
                f"(given {str(recorded)[:12]}..., recomputed {actual[:12]}...). "
                "Hash the file in the run that reads it."
            )
    return reasons


def check_access(
    access: Access,
    *,
    ledger: Iterable[dict[str, Any]] | None = None,
    registration: Mapping[str, Any] | None = None,
    root: str | Path = ".",
) -> list[str]:
    """Reasons the access would be refused. An empty list means it is allowed.

    Only a ``values`` access is gated. Checksum, schema and availability
    accesses are always allowed, but they are still recorded: the number of
    looks at the external cohort is itself a registration commitment.
    """
    reasons: list[str] = []
    if access.cohort not in HELD_OUT_COHORTS:
        reasons.append(
            f"{access.cohort!r} is not a held-out cohort ({', '.join(HELD_OUT_COHORTS)})."
        )
    if access.kind not in KINDS:
        reasons.append(f"kind {access.kind!r} is not one of {', '.join(KINDS)}.")
    for name in _REQUIRED_TEXT:
        if not str(getattr(access, name)).strip():
            reasons.append(f"{name} is empty: every access says what it is and why.")
    if access.kind != VALUES or reasons:
        return reasons

    state = registration_state() if registration is None else registration
    if not is_deposited(state):
        reasons.append(
            "REGISTRATION GATE: the pre-registration is not deposited "
            f"({REGISTRATION_PATH} carries no OSF DOI). The deposit must precede "
            "any read of the external test partition by a model."
        )

    if access.analysis.startswith(EXPLORATORY_PREFIX):
        if not (access.exploratory_preregistration or "").strip():
            reasons.append(
                "REGISTRATION GATE: an exploratory analysis of the external cohort "
                "needs its own pre-registration; cite it in exploratory_preregistration."
            )
    else:
        try:
            get_threshold(access.analysis)
        except UnregisteredMetricError:
            reasons.append(
                f"REGISTRATION GATE: {access.analysis!r} is not one of the registered "
                "metrics, and it is not labelled exploratory."
            )

    if access.partition:
        manifest_path = Path(root) / access.partition
        try:
            sealed = read_manifest(manifest_path).sealed
        except FileNotFoundError:
            reasons.append(
                f"PARTITION GATE: manifest {access.partition} does not exist. "
                "Seal it before reading."
            )
        except (ValueError, KeyError, OSError) as exc:
            reasons.append(
                f"PARTITION GATE: manifest {access.partition} is unreadable ({exc})."
            )
        else:
            if not sealed:
                reasons.append(
                    f"PARTITION GATE: manifest {access.partition} is an unsealed draft."
                )
    elif not (access.scoring_rule or "").strip():
        reasons.append(
            "PARTITION GATE: give the sealed partition the analysis is scored on "
            "(partition), or the id of the registered rule that fixes its scoring set "
            f"(scoring_rule; see {SCORING_RULES_PATH.as_posix()})."
        )
    else:
        reasons.extend(_scoring_rule_reasons(access.analysis, access.scoring_rule, root))

    if not access.inputs:
        reasons.append(
            "A value read must record the SHA-256 of every file it opens "
            "(inputs=hash_inputs([...]))."
        )
    else:
        reasons.extend(_input_reasons(access.inputs, root))

    earlier = [
        entry
        for entry in (ledger or [])
        if entry.get("cohort") == access.cohort
        and entry.get("analysis") == access.analysis
        and entry.get("kind") == VALUES
    ]
    if earlier and not (access.repeat_reason or "").strip():
        when = ", ".join(str(entry.get("date")) for entry in earlier)
        reasons.append(
            f"ONCE PER ANALYSIS: {access.analysis} was already read on {when}. "
            "A re-run is a deliberate, reported decision: pass repeat_reason."
        )
    return reasons


def record_access(
    access: Access,
    *,
    ledger_path: str | Path = LEDGER_PATH,
    registration_path: str | Path = REGISTRATION_PATH,
    root: str | Path = ".",
) -> dict[str, Any]:
    """Gate the access, then append it to the ledger.

    Raises
    ------
    ReadRefusedError
        With every reason, if any gate fails. A refused access is not written.
    """
    ledger = load_ledger(ledger_path)
    state = registration_state(registration_path)
    reasons = check_access(access, ledger=ledger, registration=state, root=root)
    if reasons:
        raise ReadRefusedError(
            "Access to the held-out cohort refused:\n  - " + "\n  - ".join(reasons)
        )
    now = datetime.now(timezone.utc)
    entry: dict[str, Any] = {
        "schema_version": LEDGER_SCHEMA_VERSION,
        "seq": 1 + max((int(e.get("seq", -1)) for e in ledger), default=-1),
        "date": now.date().isoformat(),
        "recorded_at": now.isoformat(timespec="seconds"),
        "backfilled": False,
        "before_deposit": not is_deposited(state),
    }
    entry.update(asdict(access))
    target = Path(ledger_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(entry, sort_keys=True, ensure_ascii=True) + "\n")
    return entry


def append_ledger_entry(
    kind: str,
    analysis: str,
    purpose: str,
    *,
    cohort: str = HELD_OUT_COHORTS[0],
    inputs: Mapping[str, str] | None = None,
    script: str = "",
    source: str = "",
    details: Mapping[str, Any] | None = None,
    ledger_path: str | Path = LEDGER_PATH,
    registration_path: str | Path = REGISTRATION_PATH,
) -> dict[str, Any]:
    """Append an observation or a note to the ledger, and return it.

    Neither is an access declaration, so no gate applies; neither is a read,
    so neither counts as one. An ``observation`` records an access seen in
    file-system metadata and made outside this module. A ``note`` is a dated
    statement about earlier entries: a follow-up, a correction, a closing.
    The ledger is append-only, so a note is how an earlier entry is ever
    annotated. Every entry carries a UTC time of day.

    Raises
    ------
    ValueError
        If ``kind`` is neither an observation nor a note, or the text is empty.
    """
    if kind not in NON_ACCESS_KINDS:
        raise ValueError(f"kind {kind!r} is not one of {', '.join(NON_ACCESS_KINDS)}")
    if not str(analysis).strip() or not str(purpose).strip():
        raise ValueError(
            "analysis and purpose are required: every entry says what it is and why"
        )
    ledger = load_ledger(ledger_path)
    now = datetime.now(timezone.utc)
    entry: dict[str, Any] = {
        "schema_version": LEDGER_SCHEMA_VERSION,
        "seq": 1 + max((int(e.get("seq", -1)) for e in ledger), default=-1),
        "date": now.date().isoformat(),
        "recorded_at": now.isoformat(timespec="seconds"),
        "backfilled": False,
        "before_deposit": not is_deposited(registration_state(registration_path)),
        "cohort": cohort,
        "kind": kind,
        "analysis": str(analysis),
        "purpose": str(purpose),
        "script": script,
        "inputs": dict(inputs or {}),
        "partition": None,
        "scoring_rule": None,
        "exploratory_preregistration": None,
        "repeat_reason": None,
        "source": source,
        "details": dict(details or {}),
    }
    json.dumps(entry)  # it must be serialisable before anything is written
    target = Path(ledger_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(entry, sort_keys=True, ensure_ascii=True) + "\n")
    return entry


def summarize(entries: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Counts a reader needs before opening any entry."""
    rows = list(entries)
    value_reads = [entry for entry in rows if entry.get("kind") == VALUES]
    return {
        "entries": len(rows),
        "by_kind": {kind: sum(1 for e in rows if e.get("kind") == kind) for kind in KINDS},
        "by_cohort": {
            cohort: sum(1 for e in rows if e.get("cohort") == cohort)
            for cohort in sorted({str(e.get("cohort")) for e in rows})
        },
        "value_reads": len(value_reads),
        "value_read_analyses": [str(entry.get("analysis")) for entry in value_reads],
        "non_access_entries": sum(
            1 for entry in rows if entry.get("kind") in NON_ACCESS_KINDS
        ),
        "before_deposit": sum(1 for entry in rows if entry.get("before_deposit")),
    }


def format_ledger(entries: Iterable[dict[str, Any]]) -> str:
    """The ledger as plain text, with each value read numbered in order."""
    rows = list(entries)
    if not rows:
        return "Access ledger: empty. The held-out cohort has never been accessed."
    lines = [f"Access ledger: {len(rows)} entries, oldest first."]
    read_number = 0
    for entry in rows:
        flags: list[str] = []
        if entry.get("kind") == VALUES:
            read_number += 1
            flags.append(f"value read {read_number}")
        if entry.get("before_deposit"):
            flags.append("before deposit")
        if entry.get("backfilled"):
            flags.append("back-filled")
        if entry.get("kind") in NON_ACCESS_KINDS:
            flags.append("not an access")
        lines.append(
            f"#{entry.get('seq'):<4}{entry.get('date')}  {str(entry.get('kind')):<13}"
            f"{str(entry.get('analysis')):<44}[{', '.join(flags)}]"
        )
        lines.append(f"      {entry.get('purpose')}")
        if entry.get("details"):
            lines.append(
                "      details: "
                + json.dumps(entry["details"], sort_keys=True)[:240]
            )
        if entry.get("source"):
            lines.append(f"      source: {entry['source']}")
    summary = summarize(rows)
    lines.append(
        f"Value reads of the external cohort: {summary['value_reads']} "
        f"({', '.join(summary['value_read_analyses']) or 'none'})."
    )
    return "\n".join(lines)
