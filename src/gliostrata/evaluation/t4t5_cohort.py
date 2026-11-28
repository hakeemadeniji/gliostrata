"""Turn the sealed external cohort's own files into per-patient inputs for T4 and T5.

Library code rather than script code, so it can be tested without running an evaluation.
`scripts/t4t5_evaluate.py` is the only caller, and it reaches here only after the read guard
has allowed a ``values`` access.

**The modifier-source rule, declared before the read.** The clinical table carries
``tme_immune`` and ``state_mes``, which are the patient modifiers item 11 uses. Whether they
are populated for this cohort is not known here: an ``object`` dtype on a numeric column in
this table usually means every value is null — ``kps`` has exactly that shape and the cohort
records no KPS at all. Rather than depend on which case holds, the rule covers both:

    per patient, the clinical table's own `tme_immune` and `state_mes` where BOTH are
    non-null for that patient; otherwise BOTH are deconvolved from that patient's
    expression. The source used is recorded per patient.

Mixing the two within one patient is refused, because an immune score from the table and a
MES fraction from deconvolution are not on comparable scales. The table is preferred because
it is the cohort's own curated value; deconvolution is the fallback, not the default. Which
source each patient used is reported, so a reader is told rather than left to infer.

**Standardisation is against the frozen TCGA reference, never this cohort.** Both the tropism
score and the modifiers are z-scored against
`configs/response_surface_v1_tropism_reference.json`; a scale taken from the cohort being
scored would let it into its own transform (item 11(h)).
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.optimization.config_space import CELL_TYPES
from gliostrata.prediction import response_surface as rs
from gliostrata.simulation import delivery
from gliostrata.simulation.parameters import ParameterSet

MODIFIER_SOURCE_RULE = (
    "per patient: the clinical table's own tme_immune and state_mes where BOTH are non-null "
    "for that patient; otherwise BOTH are deconvolved from that patient's expression. The "
    "source used is recorded per patient. Mixing the two within one patient is refused, "
    "because the two scales are not comparable."
)
CLINICAL_MODIFIERS = {"immune_tme_slope_per_sd": "tme_immune",
                      "mes_like_slope_per_sd": "state_mes"}


def read_expression(paths: list[Path], wanted: set[str]) -> pd.DataFrame:
    """Gene x sample log2(FPKM + 1) over the wanted samples, from CGGA's RSEM archives.

    The archives are per batch (mRNAseq_325, mRNAseq_693) and a partition spans both, so the
    columns are taken from each and joined. Duplicate gene symbols are averaged, as
    everywhere else in the project.
    """
    frames = []
    for path in paths:
        with zipfile.ZipFile(path) as handle:
            members = [m for m in handle.namelist()
                       if not m.endswith("/") and not m.startswith("__MACOSX/")
                       and not Path(m).name.startswith("._")]
            if len(members) != 1:
                raise ValueError(f"{path}: expected one member, found {members}")
            with handle.open(members[0]) as member:
                frame = pd.read_csv(member, sep="\t", index_col=0)
        frame.index = frame.index.astype(str)
        keep = [c for c in frame.columns if str(c) in wanted]
        if keep:
            frames.append(frame[keep])
    if not frames:
        raise ValueError("none of the requested patients appear in the expression archives")
    joined = pd.concat(frames, axis=1).groupby(level=0).mean()
    return np.log2(joined.astype(float) + 1.0)


def standardise(value: float, slope: str, reference_doc: dict) -> float:
    """z against the frozen TCGA reference. Refuses a modifier that has none."""
    entry = (reference_doc.get("modifier_reference") or {}).get(slope)
    if not entry:
        raise ValueError(
            f"no frozen reference for modifier {slope!r}. Standardising it on the cohort "
            "being scored would be the leakage the reference exists to prevent (item 11(h)).")
    return (float(value) - float(entry["mean"])) / float(entry["sd"])


def baseline_curves(clinical: pd.DataFrame, patient_ids: list[str], params: ParameterSet,
                    module, tcga: pd.DataFrame) -> tuple[dict[str, np.ndarray], str, dict]:
    """Each patient's Stupp baseline: the registered model with the arm forced to RT+TMZ."""
    # The caller supplies the development frame; this module does not guess its path.
    fitted = module.fit_variant(tcga, module.VARIANT_2)
    stupp = clinical.loc[patient_ids].copy()
    levels = fitted.categorical_levels.get("treatment_arm", [])
    target = next((lvl for lvl in levels
                   if "tmz" in str(lvl).lower() or "temo" in str(lvl).lower()), None)
    if target is None:
        raise ValueError(f"no RT+TMZ level among {levels}; item 11(a) cannot be applied")
    stupp["treatment_arm"] = target

    design, _, _, report = module.build_design_matrix(
        stupp.reset_index(), module.VARIANT_2, impute="drop")
    design = pd.DataFrame(design, columns=report.feature_names).reindex(
        columns=fitted.design_columns, fill_value=0.0).to_numpy()
    log_hazard = fitted.model.predict_log_hazard(design)
    grid = rs.time_grid(params)
    cumhaz = fitted.model.baseline_cumulative_hazard(grid)

    curves = {pid: np.clip(np.exp(-cumhaz * np.exp(float(log_hazard[i]))), 0.0, 1.0)
              for i, pid in enumerate(patient_ids)}
    provenance = {"features": list(module.VARIANT_2), "n_fit": int(fitted.n_fit),
                  "arm_level_used": str(target),
                  "internal_cindex": float(fitted.internal_cindex)}
    return curves, str(target), provenance


def load(expression_paths: list[Path], clinical_path: Path, patient_ids: list[str],
         params: ParameterSet, module, tcga: pd.DataFrame,
         scorer_factory) -> tuple[dict, dict]:
    """Build one scorer per patient. Returns (scorers, provenance)."""
    reference = delivery.load_tropism_reference()
    reference_doc = json.loads(delivery.DEFAULT_REFERENCE.read_text(encoding="utf-8"))

    clinical = pd.read_parquet(clinical_path)
    clinical["patient_id"] = clinical["patient_id"].astype(str)
    clinical = clinical[clinical["patient_id"].isin(patient_ids)].set_index("patient_id")
    absent = [p for p in patient_ids if p not in clinical.index]
    if absent:
        raise ValueError(f"{len(absent)} partition patients absent from the clinical table")

    expression = read_expression(expression_paths, set(patient_ids))
    no_expression = [p for p in patient_ids if p not in expression.columns]
    if no_expression:
        raise ValueError(f"{len(no_expression)} partition patients have no expression column")

    genes = sorted({g for carrier in CELL_TYPES for g in params.tropism_genes(carrier)}
                   & set(expression.index))
    tropism = {carrier: delivery.tropism_score(expression.loc[genes], carrier, params,
                                               reference)
               for carrier in CELL_TYPES}

    table = {slope: (clinical[column] if column in clinical.columns
                     else pd.Series(index=clinical.index, dtype="float64"))
             for slope, column in CLINICAL_MODIFIERS.items()}
    from_table = {p for p in patient_ids
                  if all(pd.notna(table[s].get(p)) for s in CLINICAL_MODIFIERS)}

    deconvolved: dict[str, pd.Series] = {}
    if len(from_table) < len(patient_ids):
        from gliostrata.profiling.deconvolve import (
            deconvolve_states_expression,
            deconvolve_tme_expression,
        )
        states = deconvolve_states_expression(expression)
        tme = deconvolve_tme_expression(expression)
        deconvolved = {"mes_like_slope_per_sd": states.get("MES"),
                       "immune_tme_slope_per_sd": tme.get("immune")}
        if any(v is None for v in deconvolved.values()):
            raise ValueError(
                "deconvolution produced no MES fraction or immune score, and "
                f"{len(patient_ids) - len(from_table)} patients need the fallback")

    curves, arm, survival_provenance = baseline_curves(
        clinical, patient_ids, params, module, tcga)

    scorers, source_by_patient = {}, {}
    for pid in patient_ids:
        if pid in from_table:
            raw = {s: float(table[s][pid]) for s in CLINICAL_MODIFIERS}
            source_by_patient[pid] = "clinical_table"
        else:
            raw = {s: float(deconvolved[s][pid]) for s in CLINICAL_MODIFIERS}
            source_by_patient[pid] = "deconvolved_from_expression"
        scorers[pid] = scorer_factory(
            params, curves[pid],
            {carrier: float(tropism[carrier][pid]) for carrier in CELL_TYPES},
            {s: standardise(raw[s], s, reference_doc) for s in CLINICAL_MODIFIERS})

    provenance = {
        "n_patients": len(patient_ids),
        "modifier_source_rule": MODIFIER_SOURCE_RULE,
        "n_modifiers_from_clinical_table": len(from_table),
        "n_modifiers_deconvolved": len(patient_ids) - len(from_table),
        "modifier_source_by_patient": source_by_patient,
        "standardised_against": delivery.DEFAULT_REFERENCE.as_posix(),
        "tropism_genes_used": genes,
        "survival_model": survival_provenance,
        "stupp_arm_level": arm,
        "expression_archives": [p.as_posix() for p in expression_paths],
    }
    return scorers, provenance
