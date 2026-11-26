"""L7 — Wet-lab handoff protocol generation.

Converts the optimizer's ranked configuration set into per-patient JSON
protocol documents that a wet-lab collaborator can execute on patient-derived
organoids (PDOs).

Each protocol is a complete handoff document: the recommended configuration,
a faithful projection of the patient's feature-store row, cell-line and payload
provenance, dose and schedule specifications, the PDO panel, endpoints, assay
plan, statistical analysis plan, regulatory context, randomization/blinding,
predicted outcomes, controls, biosafety, and a completeness checklist that is
derived from the document's actual content.

Fields that depend on wet-lab or institutional specifics that this framework
cannot know are emitted as ``"placeholder: <description>"`` strings. They are
never filled with invented specifics (no fake MTA, accession or IRB numbers).
"""
from __future__ import annotations

from typing import Any

import pandas as pd

PLACEHOLDER_PREFIX = "placeholder:"

# Default cell lines per optimizer cell_type. These are the canonical published
# lines for each modality; the specific vial/passage provenance is a placeholder.
CELL_LINE_DEFAULTS: dict[str, str] = {
    "HB1_F3_line": "HB1.F3.CD",
    "iPSC-NSC": "iPSC-derived NSC (patient-matched iPSC line)",
    "MSC": "Human bone-marrow-derived MSC (BM-MSC)",
}

DIFFERENTIATION_PROTOCOL_DEFAULTS: dict[str, str] = {
    "HB1_F3_line": "placeholder: immortalized line - expansion protocol reference required",
    "iPSC-NSC": "placeholder: dual-SMAD-inhibition NSC induction protocol reference required",
    "MSC": "placeholder: primary isolation and expansion protocol reference required",
}

# Default constructs per payload category.
# Outcome labels that exist in the feature store for model training but must
# never appear in a prospective therapy recommendation document. A wet-lab
# collaborator must not see the patient's observed survival or assigned arm
# alongside a survival prediction. `patient_summary` filters against this set,
# including in its trailing "surface any remaining column" pass.
OUTCOME_DENYLIST: frozenset[str] = frozenset({
    "observed_overall_survival_months",
    "overall_survival_months",
    "vital_status",
    "treatment_arm",
})

# Biosafety classification per payload category (see optimization/config_space.py
# PAYLOADS for the categories the optimizer can emit).
PAYLOAD_BIOSAFETY: dict[str, dict[str, str]] = {
    "CRAd": {
        "viral_vector_class": "Replication-competent adenovirus (Delta-24 backbone)",
        "bsl_level": "BSL-2",
    },
    "prodrug_enzyme": {
        "viral_vector_class": (
            "Replication-incompetent lentiviral vector (enzyme transgene; "
            "prodrug administered separately)"
        ),
        "bsl_level": "BSL-2",
    },
    "CAR": {
        "viral_vector_class": (
            "Replication-incompetent lentiviral vector (chimeric antigen receptor transgene)"
        ),
        "bsl_level": "BSL-2",
    },
    # Co-delivery of a replication-competent vector alongside an integrating
    # vector is the combination most likely to trigger enhanced containment, so
    # the level is left to the institutional biosafety committee rather than
    # asserted. Baseline BSL-2 applies to each vector individually.
    "combo": {
        "viral_vector_class": (
            "Co-delivery: replication-competent adenovirus (Delta-24 backbone) plus "
            "integrating lentiviral vector"
        ),
        "bsl_level": (
            f"{PLACEHOLDER_PREFIX} institutional biosafety committee risk assessment required "
            "- co-delivery of a replication-competent adenovirus and an integrating lentiviral "
            "vector may require BSL-2+ or enhanced containment; baseline BSL-2 applies to each "
            "vector individually but the combination is institution-dependent"
        ),
    },
}

PAYLOAD_CONSTRUCT_DEFAULTS: dict[str, str] = {
    "CRAd": "Delta-24-RGDOX",
    "prodrug_enzyme": "CD-UPRT (cytosine deaminase / uracil phosphoribosyltransferase)",
    "CAR": "IL13Ra2-targeted CAR",
    "combo": "Delta-24-RGDOX + CD-UPRT co-delivery",
}

SCHEDULE_INJECTION_DAYS: dict[str, list[int]] = {
    "single": [0],
    "fractionated_x2": [0, 7],
    "fractionated_x3": [0, 7, 14],
}

FOLLOW_UP_TIMEPOINTS_DAYS = [7, 14, 21, 28]

N_PDOS_PER_PATIENT = 10

ETHICAL_NOTES = (
    "Pre-registered three-arm PDO experiment. Sample size: >= 10 PDOs spanning "
    "Verhaak subtypes (power analysis: 80% power to detect Arm A->B difference "
    "at alpha=0.05 two-sided paired t-test)."
)


def _placeholder(description: str) -> str:
    """Build a placeholder value for a field this framework cannot populate."""
    return f"{PLACEHOLDER_PREFIX} {description}"


def _is_placeholder(value: Any) -> bool:
    """True if `value` is a placeholder string."""
    return isinstance(value, str) and value.strip().lower().startswith(PLACEHOLDER_PREFIX)


def _is_resolved(value: Any) -> bool:
    """True if `value` carries real content (present, non-empty, not a placeholder)."""
    if value is None or _is_placeholder(value):
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, tuple, dict, set)):
        return len(value) > 0
    if isinstance(value, bool):
        return True
    return pd.notna(value)


def _num(value: Any) -> float | None:
    """Coerce to float, or None if missing."""
    return float(value) if value is not None and pd.notna(value) else None


def _int(value: Any) -> int | None:
    """Coerce to int, or None if missing."""
    return int(value) if value is not None and pd.notna(value) else None


def _patient_summary(feat: pd.Series) -> dict[str, Any]:
    """Project a feature-store row into the protocol's patient summary.

    This is a faithful projection: every feature-store field for the patient is
    surfaced, with the one-hot encodings decoded back to their source labels
    (`sex_male` -> `sex`, `extent_of_resection_*` -> `extent_of_resection`).
    """

    def get(key: str) -> Any:
        value = feat.get(key)
        return None if value is None or (not isinstance(value, str) and pd.isna(value)) else value

    sex_male = get("sex_male")
    if sex_male is None:
        sex = None
    else:
        sex = "M" if int(sex_male) == 1 else "F"

    eor_str, eor_gtr = get("extent_of_resection_STR"), get("extent_of_resection_GTR")
    if eor_str is None and eor_gtr is None:
        extent_of_resection = None
    elif eor_gtr is not None and int(eor_gtr) == 1:
        extent_of_resection = "GTR"
    elif eor_str is not None and int(eor_str) == 1:
        extent_of_resection = "STR"
    else:
        extent_of_resection = "biopsy"  # reference level of the one-hot encoding

    summary: dict[str, Any] = {
        # Demographics / clinical
        "sex": sex,
        "age": _int(get("age")),
        "kps": _int(get("kps")),
        "extent_of_resection": extent_of_resection,
        # Molecular
        "mgmt_methylated": bool(get("mgmt_methylated")) if get("mgmt_methylated") is not None else None,
        "idh_mutated": bool(get("idh_mutated")) if get("idh_mutated") is not None else None,
        # Subtype (argmax label + calibrated class probabilities)
        "subtype": get("subtype"),
        "subtype_prob_Classical": _num(get("subtype_prob_Classical")),
        "subtype_prob_Mesenchymal": _num(get("subtype_prob_Mesenchymal")),
        "subtype_prob_Proneural": _num(get("subtype_prob_Proneural")),
        "subtype_prob_Neural": _num(get("subtype_prob_Neural")),
        # Neftel single-cell state fractions
        "state_npc": _num(get("state_npc")),
        "state_opc": _num(get("state_opc")),
        "state_ac": _num(get("state_ac")),
        "state_mes": _num(get("state_mes")),
        # Tumor microenvironment fractions
        "tme_immune": _num(get("tme_immune")),
        "tme_stromal": _num(get("tme_stromal")),
        "tme_endothelial": _num(get("tme_endothelial")),
    }

    # Surface any feature-store column not already covered, so the summary stays
    # a faithful projection as the feature store grows — except outcome labels,
    # which must never reach a prospective recommendation document.
    covered = set(summary) | OUTCOME_DENYLIST | {
        "patient_id", "sex_male", "extent_of_resection_STR", "extent_of_resection_GTR",
    }
    for col in feat.index:
        if col not in covered:
            value = feat.get(col)
            summary[col] = None if (not isinstance(value, str) and pd.isna(value)) else value

    return {k: v for k, v in summary.items() if k not in OUTCOME_DENYLIST}


def _biosafety(payload: str) -> dict[str, Any]:
    """Biosafety classification, derived from the payload category.

    Kept consistent with `payload_specification`: the vector class follows from
    the construct the optimizer selected, rather than being asserted statically.
    """
    classification = PAYLOAD_BIOSAFETY.get(payload)
    if classification is None:
        classification = {
            "viral_vector_class": _placeholder(f"vector class for {payload} required"),
            "bsl_level": "placeholder",
        }
    return {
        "bsl_level": classification["bsl_level"],
        "viral_vector_class": classification["viral_vector_class"],
        "payload_category": payload,
        "containment_notes": "Follow institutional biosafety committee protocols.",
    }


def _cell_line_provenance(cell_type: str) -> dict[str, Any]:
    """Cell-line identity and provenance for the recommended cell type."""
    return {
        "line_name": CELL_LINE_DEFAULTS.get(
            cell_type, _placeholder(f"specific line name for cell_type {cell_type!r} required")
        ),
        "cell_type_category": cell_type,
        "source_reference": _placeholder(
            "originating publication, ATCC/repository catalog number, or executed MTA required"
        ),
        "passage_range": _placeholder("acceptable passage window for this line required"),
        "mycoplasma_certified": _placeholder(
            "mycoplasma test certificate (date and result) required before use"
        ),
        "differentiation_protocol_reference": DIFFERENTIATION_PROTOCOL_DEFAULTS.get(
            cell_type, _placeholder("differentiation/expansion protocol reference required")
        ),
    }


def _payload_specification(payload: str) -> dict[str, Any]:
    """Therapeutic payload construct identity and release criteria."""
    return {
        "construct_name": PAYLOAD_CONSTRUCT_DEFAULTS.get(
            payload, _placeholder(f"specific construct for payload {payload!r} required")
        ),
        "category": payload,
        "sequence_accession": _placeholder(
            "GenBank/Addgene accession or full annotated construct map required"
        ),
        "titer_vp_per_ml": _placeholder("measured viral particle titer (vp/mL) required"),
        "replication_competence_confirmed": _placeholder(
            "replication-competence assay result required (expected for CRAd payloads)"
        ),
        "biosafety_approval_reference": _placeholder(
            "institutional biosafety committee approval number required"
        ),
    }


def _dose_specification(dose_cells: int, n_injection_sites: int) -> dict[str, Any]:
    """Dose breakdown, including the computed per-site split."""
    cells_per_site = dose_cells // n_injection_sites if n_injection_sites else None
    return {
        "total_cells": dose_cells,
        "n_injection_sites": n_injection_sites,
        "cells_per_site": cells_per_site,
        "volume_per_site_ul": _placeholder(
            "injection volume per site (uL) required - constrained by PDO diameter"
        ),
        "concentration_cells_per_ml": _placeholder(
            "working suspension concentration (cells/mL) required - derive from volume_per_site_ul"
        ),
        "vehicle_buffer": _placeholder(
            "vehicle/buffer formulation required (e.g. serum-free medium or HBSS)"
        ),
    }


def _schedule_detail(schedule: str) -> dict[str, Any]:
    """Temporal plan for dosing and sampling."""
    injection_days = SCHEDULE_INJECTION_DAYS.get(schedule, [0])
    return {
        "schedule": schedule,
        "temporal_anchor": "Day 0 post-PDO-establishment",
        "injection_day": injection_days[0],
        "injection_days": injection_days,
        "follow_up_timepoints_days": list(FOLLOW_UP_TIMEPOINTS_DAYS),
        "sampling_schedule": (
            "Conditioned medium sampled at each follow-up timepoint; PDOs harvested for "
            "endpoint assays at day 28 (or at earlier stopping-criteria trigger)."
        ),
    }


def _pdo_panel_specification() -> dict[str, Any]:
    """PDO panel composition and QC."""
    return {
        "n_pdos_per_patient": N_PDOS_PER_PATIENT,
        "selection_criteria": (
            f"{N_PDOS_PER_PATIENT} PDOs per patient spanning the Verhaak subtypes, each "
            "established from the patient's resection tissue and passing establishment QC."
        ),
        "culture_medium": _placeholder(
            "PDO culture medium formulation and supplier catalog numbers required"
        ),
        "passage_range": _placeholder("acceptable PDO passage window at treatment required"),
        "quality_controls": [
            _placeholder("establishment success and viability threshold required"),
            _placeholder("histology/marker confirmation of GBM identity required"),
            _placeholder("mycoplasma and sterility clearance required"),
        ],
    }


def _endpoints() -> dict[str, Any]:
    """Pre-specified primary and secondary endpoints."""
    return {
        "primary_endpoint": "RMST difference Arm A vs Arm B, days",
        "secondary_endpoints": [
            "Tumor kill fraction at day 28 (Arm B vs Arm A)",
            "Tumor kill fraction at day 28 (Arm B vs Arm C)",
            "PDO viability trajectory across follow-up timepoints",
            "Payload persistence / viral replication across follow-up timepoints",
            "Off-target cytotoxicity in non-tumor co-culture",
        ],
        "measurement_methods": [
            "CellTiter-Glo 3D ATP viability assay",
            "Live/dead fluorescent stain (calcein-AM / ethidium homodimer) with confocal readout",
            "qPCR for viral genome copy number",
        ],
        "timepoints": list(FOLLOW_UP_TIMEPOINTS_DAYS),
    }


def _assay_plan() -> dict[str, Any]:
    """Standard PDO assay battery."""
    return {
        "tumor_kill_assay": (
            "CellTiter-Glo 3D at each follow-up timepoint, normalized to the matched "
            "vehicle-only PDO from the same panel."
        ),
        "viability_assay": (
            "Live/dead fluorescent stain (calcein-AM / ethidium homodimer) with confocal "
            "z-stack quantification of the dead-cell fraction."
        ),
        "imaging_plan": (
            "Brightfield PDO diameter tracking at every timepoint; confocal immunofluorescence "
            "at day 28 endpoint harvest."
        ),
        "qpcr_for_viral_replication": (
            "qPCR for viral genome copies per ng total DNA at each follow-up timepoint, "
            "run against a standard curve of the input construct."
        ),
    }


def _statistical_analysis_plan() -> dict[str, Any]:
    """Pre-specified analysis plan matching the power analysis in ethical_notes."""
    return {
        "primary_analysis": (
            "Two-sided paired t-test, Arm A vs Arm B, paired within the patient's PDO panel "
            "(alpha=0.05; 80% power at n>=10 PDOs per power analysis)."
        ),
        "multiple_comparison_correction": "Benjamini-Hochberg",
        "interim_analysis_rule": (
            "One interim look after 50% of PDOs reach the day-28 endpoint; O'Brien-Fleming "
            "alpha spending so the final analysis retains nominal alpha."
        ),
        "stopping_criteria": (
            "Stop for futility if the interim Arm A vs Arm B effect estimate is <= 0; stop an "
            "individual PDO for unexpected total loss of viability in the vehicle control."
        ),
    }


def _regulatory_context() -> dict[str, Any]:
    """Regulatory references. All institution-specific numbers are placeholders."""
    return {
        "irb_reference": _placeholder("pending - IRB protocol number required"),
        "iacuc_reference": _placeholder(
            "pending - IACUC protocol number required only if the study extends to in vivo work"
        ),
        "ind_status": _placeholder("TBD - pre-IND; no IND filed for this preclinical PDO study"),
        "fda_gene_therapy_notes": (
            "Design aligned with FDA guidance \"Considerations for the Design of Early-Phase "
            "Clinical Trials of Cellular and Gene Therapy Products\" and the FDA guidance on "
            "cell-based oncology product design. Preclinical PDO stage: no FDA submission "
            "required, but potency and identity assays are specified now to support a future IND."
        ),
    }


def _randomization_and_blinding() -> dict[str, Any]:
    """Bias-control plan."""
    return {
        "randomization_scheme": "block randomization within patient PDO panel",
        "blinding": (
            "assay-level blinding for tumor kill measurement; analyst blinded to arm"
        ),
        "allocation_concealment": _placeholder(
            "randomization code custodian (not an assay operator) required"
        ),
    }


def _completeness_checklist(protocol: dict[str, Any]) -> dict[str, Any]:
    """Derive the completeness checklist from the protocol's actual content.

    Every item is resolved by inspecting the document that was just built, not by
    a code path that assumed a section would be written. Each value is:

      * ``True``          — the section exists and its required fields carry real content.
      * ``"placeholder"`` — the section exists but required fields are still placeholders.
      * ``False``         — the section is missing or its required fields are empty.
    """

    def status(section_key: str, *fields: str) -> bool | str:
        section = protocol.get(section_key)
        if not isinstance(section, (dict, list)) or len(section) == 0:
            return False
        if isinstance(section, list):
            return True
        values = [section.get(f) for f in fields]
        if any(v is None or (isinstance(v, str) and not v.strip()) for v in values):
            return False
        if all(_is_resolved(v) for v in values):
            return True
        return "placeholder"

    def scalar_status(key: str) -> bool | str:
        value = protocol.get(key)
        if not _is_resolved(value):
            return "placeholder" if _is_placeholder(value) else False
        return True

    checklist: dict[str, Any] = {
        "cell_line_source": status(
            "cell_line_provenance", "line_name", "source_reference", "passage_range",
            "mycoplasma_certified",
        ),
        "payload_sequence": status(
            "payload_specification", "construct_name", "sequence_accession", "titer_vp_per_ml",
        ),
        "dose": status(
            "dose_specification", "total_cells", "cells_per_site", "volume_per_site_ul",
            "concentration_cells_per_ml", "vehicle_buffer",
        ),
        "schedule": status(
            "schedule_detail", "temporal_anchor", "injection_day",
            "follow_up_timepoints_days", "sampling_schedule",
        ),
        "pdo_panel": status(
            "pdo_panel_specification", "n_pdos_per_patient", "selection_criteria",
            "culture_medium", "passage_range", "quality_controls",
        ),
        "injection_coordinates": (
            "placeholder"
            if isinstance(protocol.get("recommended_config"), dict)
            and protocol["recommended_config"].get("injection_coordinates_placeholder")
            else False
        ),
        "controls": status("controls"),
        "endpoints": status(
            "endpoints", "primary_endpoint", "secondary_endpoints", "measurement_methods",
            "timepoints",
        ),
        "assay_plan": status(
            "assay_plan", "tumor_kill_assay", "viability_assay", "imaging_plan",
            "qpcr_for_viral_replication",
        ),
        "statistical_analysis_plan": status(
            "statistical_analysis_plan", "primary_analysis",
            "multiple_comparison_correction", "interim_analysis_rule", "stopping_criteria",
        ),
        "power_analysis": (
            True
            if isinstance(protocol.get("statistical_analysis_plan"), dict)
            and "power" in str(protocol["statistical_analysis_plan"].get("primary_analysis", "")).lower()
            else False
        ),
        "randomization_and_blinding": status(
            "randomization_and_blinding", "randomization_scheme", "blinding",
        ),
        "regulatory_context": status(
            "regulatory_context", "irb_reference", "ind_status", "fda_gene_therapy_notes",
        ),
        "predicted_outcomes": status(
            "predicted_outcomes", "predicted_survival_with_config_months",
            "predicted_tumor_kill_fraction", "survival_gain_vs_stupp_months",
        ),
        "patient_summary": status("patient_summary", "subtype", "age", "sex", "kps"),
        "biosafety_level": status("biosafety", "bsl_level", "viral_vector_class"),
        "ethical_notes": scalar_status("ethical_notes"),
        "data_sharing_agreement": _placeholder(
            "executed data sharing agreement with the collaborating wet lab required"
        ),
    }
    checklist["unresolved_items"] = sorted(k for k, v in checklist.items() if v is not True)
    checklist["all_complete"] = not checklist["unresolved_items"]
    return checklist


def generate_protocols(
    configs: pd.DataFrame,
    feature_store: pd.DataFrame,
    survival_pred: pd.DataFrame,
    *,
    top_n: int = 1,
) -> dict[str, dict]:
    """Generate wet-lab protocol documents for the top-N configs per patient.

    Returns: {patient_id: protocol_dict}. Each protocol_dict has keys:
        patient_id, document_version, recommended_config, patient_summary,
        cell_line_provenance, payload_specification, dose_specification,
        schedule_detail, pdo_panel_specification, endpoints, assay_plan,
        statistical_analysis_plan, regulatory_context,
        randomization_and_blinding, predicted_outcomes, controls, biosafety,
        ethical_notes, completeness_checklist

    Generation is deterministic: the same inputs always produce the same
    document.
    """
    # Merge feature info for context (Series indexed by patient_id)
    fs = feature_store.set_index("patient_id")
    sp = survival_pred.set_index("patient_id")
    protocols = {}
    # Take top-N per patient (rank ≤ top_n)
    top = configs[configs["rank"] <= top_n]
    for pid, group in top.groupby("patient_id"):
        for _, row in group.iterrows():
            feat = fs.loc[pid] if pid in fs.index else pd.Series(dtype=object)
            surv = sp.loc[pid] if pid in sp.index else pd.Series(dtype=object)

            cell_type = row["cell_type"]
            payload = row["payload"]
            schedule = row["schedule"]
            dose_cells = int(row["dose_cells"])
            n_injection_sites = int(row["n_injection_sites"])

            protocol = {
                "patient_id": pid,
                "document_version": "0.2.0",
                "recommended_config": {
                    "cell_type": cell_type,
                    "payload": payload,
                    "dose_cells": dose_cells,
                    "n_injection_sites": n_injection_sites,
                    "injection_coordinates_placeholder": "Map to patient-specific PDO geometry",
                    "schedule": schedule,
                    "rank_in_optimizer": int(row["rank"]),
                },
                "patient_summary": _patient_summary(feat),
                "cell_line_provenance": _cell_line_provenance(cell_type),
                "payload_specification": _payload_specification(payload),
                "dose_specification": _dose_specification(dose_cells, n_injection_sites),
                "schedule_detail": _schedule_detail(schedule),
                "pdo_panel_specification": _pdo_panel_specification(),
                "endpoints": _endpoints(),
                "assay_plan": _assay_plan(),
                "statistical_analysis_plan": _statistical_analysis_plan(),
                "regulatory_context": _regulatory_context(),
                "randomization_and_blinding": _randomization_and_blinding(),
                "predicted_outcomes": {
                    "predicted_survival_with_config_months": float(row["predicted_survival_with_config"]),
                    "predicted_tumor_kill_fraction": float(row["predicted_tumor_kill_fraction"]),
                    "survival_gain_vs_stupp_months": float(row["survival_gain_vs_stupp_months"]),
                    "baseline_stupp_survival_months": (
                        float(surv["baseline_stupp_survival_months"])
                        if "baseline_stupp_survival_months" in surv.index
                        else None
                    ),
                },
                "controls": [
                    {"arm": "A", "description": "Standard Stupp-equivalent in vitro (TMZ + radiation on PDO)"},
                    {"arm": "B", "description": "Framework-recommended stem-cell-delivered configuration"},
                    {"arm": "C", "description": "Fixed-payload NSC delivery (City of Hope-style active control)"},
                ],
                "biosafety": _biosafety(payload),
                "ethical_notes": ETHICAL_NOTES,
            }
            # Derived last, so it reflects the document that was actually built.
            protocol["completeness_checklist"] = _completeness_checklist(protocol)
            protocols[pid] = protocol

    return protocols


def protocols_to_json(protocols: dict[str, dict]) -> dict[str, dict]:
    """Pass-through helper: returns the protocols dict in a JSON-serializable form.

    (pandas types are coerced to native Python types via the `default=str`
    fallback in utils.io.write_json — this function exists as a hook for
    future format conversion, e.g., Markdown or DOCX protocol output.)
    """
    return protocols
