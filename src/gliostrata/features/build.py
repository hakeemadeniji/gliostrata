"""Run the L2 and L3 layers and write the feature store (CLI: build-feature-store).

Two modes, one code path (``gliostrata.features.store.build_feature_store``):

``synthetic`` (default)
    Offline. Reads a harmonised synthetic cohort, runs the three L2 stubs and
    writes a store labelled ``synthetic`` in every file. This is what the tests,
    ``make pipeline`` and the offline Snakemake fixture run use.
``real``
    The development cohort (TCGA) only. Reads the L1 clinical and label tables,
    the TPM matrix and its gene annotation, the MGMT and IDH call tables, a
    frozen subtype centroid set and the frozen development partition. Runs the
    real L2:

    * subtype: single-sample nearest-centroid assignment
      (``gliostrata.profiling.centroid.classify_single_sample``) against the
      frozen centroid set. One sample's call does not depend on the other
      samples scored with it.
    * states and TME: ``gliostrata.profiling.deconvolve`` with its defaults.
      These are cohort-relative (ssGSEA then a softmax over the cohort
      z-score): a patient's value depends on the other patients scored with
      it, and the report says so.

    A stub output, a cohort-relative subtype output, a missing centroid set, or
    any path containing "cgga" stops the build before anything is written.
"""
from __future__ import annotations

import json
import platform
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path

import pandas as pd

from gliostrata.features.contract import (
    CALL_STATUSES,
    COLUMN_CLASSES,
    MISSING_REASONS,
    FeatureStore,
    FeatureStoreError,
    refuse_external,
    sha256_file,
    sibling_paths,
    write_store,
)
from gliostrata.features.store import build_feature_store, build_synthetic_store
from gliostrata.features.views import VIEWS, build_all_views, partition_lookup

#: Default inputs and outputs of the real build (paths relative to the repo).
REAL_DEFAULTS: dict[str, str] = {
    "clinical": "data/interim/tcga_gbm_clinical_harmonized.parquet",
    "labels": "data/interim/tcga_gbm_labels.parquet",
    "expression": "data/raw/tcga_gbm/expression/tpm.parquet",
    "annotation": "data/raw/tcga_gbm/expression/gene_annotation.parquet",
    "mgmt_calls": "data/raw/tcga_gbm/methylation/mgmt_calls.parquet",
    "idh_calls": "data/raw/tcga_gbm/mutations/idh_calls.parquet",
    "mgmt_manifest": "data/raw/tcga_gbm/methylation/file_manifest.parquet",
    "idh_manifest": "data/raw/tcga_gbm/mutations/file_manifest.parquet",
    "centroids": "configs/subtype_centroids/verhaak840_tcga_v2train.json",
    "partition": "data/partitions/tcga_gbm_clinical_v2.json",
    "out": "data/processed/tcga_gbm_feature_store.parquet",
    "labels_out": "data/processed/tcga_gbm_feature_store_labels.parquet",
    "missing_out": "data/processed/tcga_gbm_feature_store_missing.parquet",
    "columns_out": "data/processed/tcga_gbm_feature_store_columns.json",
    "report": "data/results/tcga_gbm/feature_store_build_report.json",
}

#: ``call_basis`` of the call tables -> call status (decision D3).
MGMT_BASIS_TO_STATUS = {
    "two_probes_measured": "called",
    "required_probe_absent_or_NA": "tested_not_callable",
    "no_methylation_data": "not_tested",
}
IDH_BASIS_TO_STATUS = {
    "hotspot_variant": "called",
    "sequenced_no_hotspot": "called",
    "no_mutation_data": "not_tested",
}


def _file_record(path: Path | str) -> dict:
    p = Path(path)
    return {"path": str(path).replace("\\", "/"), "bytes": p.stat().st_size,
            "sha256": sha256_file(p)}


def _software() -> dict:
    try:
        gs = package_version("gliostrata")
    except PackageNotFoundError:
        gs = "unknown"
    import numpy
    import pyarrow
    return {"gliostrata": gs, "python": platform.python_version(),
            "pandas": pd.__version__, "numpy": numpy.__version__,
            "pyarrow": pyarrow.__version__}


def molecular_from_call_tables(patient_ids, mgmt_calls: pd.DataFrame,
                               idh_calls: pd.DataFrame, *,
                               assayed: dict[str, set[str]] | None = None,
                               ) -> tuple[pd.DataFrame, dict]:
    """The molecular input: three-state call status, continuous MGMT score and
    the threshold used, copied from the two L0 call tables. Nothing is imputed.

    A cohort patient absent from a call table is ``not_tested``. A
    ``call_basis`` outside the known vocabulary stops the build.

    ``assayed`` maps ``"mgmt"`` / ``"idh"`` to the patients that have any
    assayed file of that modality in its GDC file manifest, whatever the sample
    type. The call tables call the primary tumour only, so a patient whose only
    assayed sample is not a primary tumour (one IDH case in TCGA-GBM: a
    recurrent sample) has ``no_*_data`` there. Such a patient was tested, but no
    admissible call exists: ``tested_not_callable``, never ``not_tested``.
    """
    ids = [str(p) for p in patient_ids]
    out = pd.DataFrame({"patient_id": ids})
    counts: dict = {}
    for label, table, value_col, basis_map, prefix in (
        ("mgmt", mgmt_calls, "mgmt_methylated", MGMT_BASIS_TO_STATUS, "mgmt"),
        ("idh", idh_calls, "idh_mutated", IDH_BASIS_TO_STATUS, "idh"),
    ):
        for col in ("patient_id", value_col, "call_basis"):
            if col not in table.columns:
                raise FeatureStoreError(f"{label} calls: column {col!r} missing")
        t = table.copy()
        t["patient_id"] = t["patient_id"].astype(str)
        if t["patient_id"].duplicated().any():
            raise FeatureStoreError(f"{label} calls: duplicate patient_id values")
        unknown = sorted(set(t["call_basis"].dropna().astype(str)) - set(basis_map))
        if unknown or t["call_basis"].isna().any():
            raise FeatureStoreError(f"{label} calls: unknown call_basis {unknown}")
        t = t.set_index("patient_id")
        in_table = [p in t.index for p in ids]
        sub = t.reindex(ids)
        status = [basis_map[str(b)] if present else "not_tested"
                  for b, present in zip(sub["call_basis"], in_table, strict=True)]
        tested = (assayed or {}).get(label)
        n_non_primary_only = 0
        if tested is not None:
            for i, pid in enumerate(ids):
                if status[i] == "not_tested" and pid in tested:
                    status[i] = "tested_not_callable"
                    n_non_primary_only += 1
        out[value_col] = sub[value_col].to_numpy()
        out[f"{prefix}_call_status"] = status
        if prefix == "mgmt":
            for col in ("mgmt_score", "threshold"):
                if col not in t.columns:
                    raise FeatureStoreError(f"mgmt calls: column {col!r} missing")
            out["mgmt_score"] = pd.to_numeric(sub["mgmt_score"]).astype("Float64").to_numpy()
            out["mgmt_threshold"] = pd.to_numeric(sub["threshold"]).astype("Float64").to_numpy()
        counts[label] = {
            "n_rows_in_table": len(t),
            "n_rows_outside_cohort": int((~t.index.isin(ids)).sum()),
            "n_cohort_patients_absent_from_table": int(len(ids) - sum(in_table)),
            "call_status": {s: status.count(s) for s in CALL_STATUSES},
            "n_assayed_but_no_primary_call": n_non_primary_only,
            "assay_manifest_used": tested is not None,
        }
    return out, counts


def _check_no_cgga(paths: dict) -> None:
    for value in paths.values():
        if value is not None:
            refuse_external(value)


def _views_report(store: FeatureStore, split_of: dict | None) -> dict:
    views = build_all_views(store, split_of)
    out = {}
    for name, result in views.items():
        view = VIEWS[name]
        out[name] = {
            "description": view.description,
            "requirements": [
                {"step": r.step, "table": r.table, "column": r.column, "rule": r.rule,
                 "reason": r.reason} for r in view.requirements],
            "feature_columns": list(result.features.columns.drop("patient_id")),
            "label_columns": list(view.label_columns),
            "n": len(result.patient_ids),
            "flow": result.flow,
        }
    return out


def _report(store: FeatureStore, inputs: dict, outputs: dict, extra: dict,
            split_of: dict | None, partition: dict | None) -> dict:
    classes: dict[str, list[str]] = {}
    for rec in store.column_records():
        classes.setdefault(rec["class"], []).append(rec["name"])
    return {
        "what": "GlioStrata L3 feature store build: every cohort patient, one row "
                "each; model-ready populations are named views",
        "mode": store.mode,
        "inputs": inputs,
        "software": _software(),
        "n_patients": len(store.features),
        "n_label_rows": len(store.labels),
        **extra,
        "missing_reason_counts": store.reason_counts(),
        "columns_by_class": classes,
        "vocabularies": {"classes": COLUMN_CLASSES, "missing_reasons": MISSING_REASONS,
                         "call_statuses": list(CALL_STATUSES)},
        "scoring": store.meta.get("scoring", {}),
        "partition": partition,
        "views": _views_report(store, split_of),
        "not_validated": [
            "Subtype calls are nearest-centroid assignments against a frozen "
            "centroid set; agreement with published labels is a descriptive "
            "check, not a validation."
            if store.mode == "real" else "Synthetic mode: every value comes from "
            "the offline simulator and the L2 stubs.",
            "State and TME values are normalised marker-enrichment scores, not "
            "cell proportions, and they are cohort-relative.",
        ],
        "outputs": outputs,
    }


def _write_report(report: dict, path: Path | str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _outputs(out: str, labels_out, missing_out, columns_out) -> dict[str, Path]:
    paths = sibling_paths(out)
    if labels_out:
        paths["labels"] = Path(labels_out)
    if missing_out:
        paths["missing"] = Path(missing_out)
    if columns_out:
        paths["columns"] = Path(columns_out)
    return paths


def run_synthetic(input_path: str, out: str, *, labels_out=None, missing_out=None,
                  columns_out=None, report=None) -> FeatureStore:
    """Synthetic mode: simulator cohort -> L2 stubs -> store files."""
    _check_no_cgga({"input": input_path, "out": out, "labels_out": labels_out,
                    "missing_out": missing_out, "columns_out": columns_out,
                    "report": report})
    from gliostrata.profiling.states import deconvolve_states
    from gliostrata.profiling.subtype import classify_subtypes
    from gliostrata.profiling.tme import deconvolve_tme
    from gliostrata.utils.io import read_parquet

    cohort = read_parquet(input_path)
    src = [_file_record(input_path)]
    store = build_synthetic_store(
        cohort, classify_subtypes(cohort), deconvolve_states(cohort),
        deconvolve_tme(cohort),
        sources={k: src for k in ("clinical", "labels", "subtypes", "states", "tme",
                                  "molecular", "expression")},
    )
    paths = _outputs(out, labels_out, missing_out, columns_out)
    written = write_store(store, paths)
    if report:
        _write_report(_report(store, {"input": src[0]}, written, {}, None, None), report)
    return store


def run_real(*, clinical: str, labels: str, expression: str, annotation: str,
             mgmt_calls: str, idh_calls: str, mgmt_manifest: str, idh_manifest: str,
             centroids: str, partition: str,
             out: str, labels_out: str | None = None, missing_out: str | None = None,
             columns_out: str | None = None, report: str | None = None) -> dict:
    """Real mode on TCGA: returns the report (also written to ``report``).

    The per-patient subtype must come from the centroid set estimated on the
    development partition's training split only (decision D6, v2 split): a set
    estimated on other patients, ``verhaak840_tcga_all`` among them, would put
    the test split's patients into every patient's subtype feature. Any other
    set stops the build.
    """
    named = {"clinical": clinical, "labels": labels, "expression": expression,
             "annotation": annotation, "mgmt_calls": mgmt_calls,
             "idh_calls": idh_calls, "mgmt_manifest": mgmt_manifest,
             "idh_manifest": idh_manifest, "centroids": centroids,
             "partition": partition}
    _check_no_cgga({**named, "out": out, "labels_out": labels_out,
                    "missing_out": missing_out, "columns_out": columns_out,
                    "report": report})
    if not Path(centroids).is_file():
        raise FeatureStoreError(
            f"real mode needs a frozen subtype centroid set; none at {centroids}. "
            "Per-patient subtype features come only from the single-sample "
            "classifier (gliostrata.profiling.centroid); nothing was written"
        )
    for name, path in named.items():
        if not Path(path).is_file():
            raise FeatureStoreError(f"real mode: input {name} not found at {path}")

    from gliostrata.evaluation import read_manifest
    from gliostrata.profiling.centroid import (
        SET_V2TRAIN,
        classify_single_sample,
        load_centroid_set,
    )
    from gliostrata.profiling.deconvolve import (
        deconvolve_states_expression,
        deconvolve_tme_expression,
    )

    inputs = {name: _file_record(path) for name, path in named.items()}
    clin = pd.read_parquet(clinical)
    lab = pd.read_parquet(labels)
    tpm = pd.read_parquet(expression)
    ann = pd.read_parquet(annotation)
    mgmt = pd.read_parquet(mgmt_calls)
    idh = pd.read_parquet(idh_calls)
    centroid_set = load_centroid_set(centroids)
    if centroid_set.set_id != SET_V2TRAIN:
        raise FeatureStoreError(
            f"real mode uses only the centroid set {SET_V2TRAIN!r} (estimated on the "
            f"v2 training split); {centroids} is {centroid_set.set_id!r}. The other "
            "set is reserved for T1 (decision D6, v2 split); nothing was written"
        )
    assayed = {
        "mgmt": set(pd.read_parquet(mgmt_manifest, columns=["patient_id"])
                    ["patient_id"].astype(str)),
        "idh": set(pd.read_parquet(idh_manifest, columns=["patient_id"])
                   ["patient_id"].astype(str)),
    }
    upstream = _read_l1_provenance(clinical)
    manifest = read_manifest(partition)

    ids = clin["patient_id"].astype(str).tolist()
    samples = [str(c) for c in tpm.columns]
    if len(set(samples)) != len(samples):
        raise FeatureStoreError("expression: duplicate sample columns")
    profiled = [p for p in ids if p in set(samples)]
    if not profiled:
        raise FeatureStoreError("no cohort patient has an expression sample")
    scored = tpm[profiled]

    subtypes = classify_single_sample(scored, ann, centroid_set)
    states = deconvolve_states_expression(scored, ann)
    tme = deconvolve_tme_expression(scored, ann)
    molecular, mol_counts = molecular_from_call_tables(ids, mgmt, idh, assayed=assayed)

    n = len(profiled)
    rel = (f"cohort-relative: ssGSEA enrichment, then softmax over the z-score "
           f"across the {n} profiled cohort patients (gliostrata.profiling.deconvolve "
           f"defaults). A patient's value depends on the other patients scored "
           f"with it")
    files = {k: inputs[k] for k in inputs}
    store = build_feature_store(
        clinical=clin, labels=lab, subtypes=subtypes, states=states, tme=tme,
        molecular=molecular, expression_samples=samples, mode="real",
        sources={
            "clinical": [files["clinical"], *upstream], "labels": [files["labels"], *upstream],
            "expression": [files["expression"]],
            "subtypes": [files["expression"], files["annotation"], files["centroids"]],
            "states": [files["expression"], files["annotation"]],
            "tme": [files["expression"], files["annotation"]],
            "molecular": [files["mgmt_calls"], files["idh_calls"],
                          files["mgmt_manifest"], files["idh_manifest"]],
        },
        scoring={
            "subtypes": (f"single_sample: nearest centroid against the frozen set "
                         f"{centroid_set.set_id} (sha256 {centroid_set.source_sha256}); "
                         "a sample's call does not depend on the other samples"),
            "states": rel, "tme": rel,
        },
    )
    split_of = partition_lookup(manifest.train_ids, manifest.test_ids)
    paths = _outputs(out, labels_out, missing_out, columns_out)
    written = write_store(store, paths)
    extra = {
        "counts": {
            "n_clinical": len(ids),
            "n_expression_samples": len(samples),
            "n_expression_samples_outside_cohort": len(set(samples) - set(ids)),
            "n_profiled": n,
            "molecular": mol_counts,
        },
        "centroid_set": {"set_id": centroid_set.set_id,
                         "sha256": centroid_set.source_sha256,
                         "n_genes": len(centroid_set.genes),
                         "subtypes": list(centroid_set.subtypes)},
    }
    partition_info = {
        "name": manifest.partition_name, "file_sha256": inputs["partition"]["sha256"],
        "n_train": len(manifest.train_ids), "n_test": len(manifest.test_ids),
        "store_patients_not_in_partition": sum(1 for p in ids if p not in split_of),
    }
    rep = _report(store, inputs, written, extra, split_of, partition_info)
    if report:
        _write_report(rep, report)
    return rep


def l1_provenance_path(clinical_out: Path | str) -> Path:
    """The provenance sidecar written next to the L1 clinical table."""
    p = Path(clinical_out)
    return p.with_name(f"{p.stem}_provenance.json")


def _read_l1_provenance(clinical: Path | str) -> list[dict]:
    """The L0 inputs of the L1 tables, from their sidecar, with checksums.

    The sidecar records the SHA-256 of the L0 files the L1 tables were made
    from, so the store's column metadata reaches back to L0. A sidecar that
    names a different clinical table, or whose recorded output checksum no
    longer matches the clinical table on disk, stops the build.
    """
    side = l1_provenance_path(clinical)
    if not side.is_file():
        raise FeatureStoreError(
            f"real mode: no L1 provenance sidecar at {side}; re-run "
            "`gliostrata harmonize-clinical`")
    doc = json.loads(side.read_text(encoding="utf-8"))
    recorded = doc.get("outputs", {}).get("clinical", {}).get("sha256")
    if recorded != sha256_file(clinical):
        raise FeatureStoreError(
            f"real mode: {clinical} does not match its L1 provenance sidecar {side}; "
            "re-run `gliostrata harmonize-clinical`")
    return [dict(rec) for rec in doc.get("inputs", {}).values()]


def run_harmonize_clinical(cohort: str, raw: str | None, clinical_out: str,
                           labels_out: str) -> dict:
    """L1: write the owned clinical table, the label table and a provenance
    sidecar (``<clinical stem>_provenance.json``) holding the SHA-256 of the L0
    inputs and of the two outputs."""
    from gliostrata.data.clinical import harmonise_clinical, reason_counts

    _check_no_cgga({"cohort": cohort, "raw": raw, "clinical_out": clinical_out,
                    "labels_out": labels_out})
    c = pd.read_parquet(cohort)
    r = pd.read_parquet(raw) if raw else None
    clinical, labels = harmonise_clinical(c, r)
    for frame, path in ((clinical, clinical_out), (labels, labels_out)):
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(p, index=False, engine="pyarrow")
    inputs = {"cohort": _file_record(cohort)}
    if raw:
        inputs["raw"] = _file_record(raw)
    outputs = {"clinical": _file_record(clinical_out), "labels": _file_record(labels_out)}
    side = l1_provenance_path(clinical_out)
    side.write_text(json.dumps({"layer": "L1", "command": "gliostrata harmonize-clinical",
                                "inputs": inputs, "outputs": outputs},
                               indent=2, sort_keys=True) + "\n",
                    encoding="utf-8", newline="\n")
    return {"n_patients": len(clinical),
            "clinical_reasons": reason_counts(clinical),
            "label_reasons": reason_counts(labels),
            "inputs": inputs, "outputs": {**outputs, "provenance": str(side)}}
