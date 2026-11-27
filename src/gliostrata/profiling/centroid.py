"""Single-sample Verhaak subtype assignment by nearest frozen centroid.

WHY THIS EXISTS
    ``gliostrata.profiling.verhaak.classify_verhaak`` standardises each
    signature score across the samples it is given, so a patient's call depends
    on who else is in the call (the scoring-set dependence recorded in stage 3
    and measured in stage 3b). A per-patient feature must not depend on its
    batch. This module assigns a subtype to each sample on its own, against
    centroids that were estimated once on TCGA, frozen to a file and
    checksummed.

THE RULE, EXACTLY
    1. Genes: the published Verhaak et al. (2010) ClaNC 840-gene signature,
       each symbol resolved to one Ensembl gene id that was fixed when the
       centroids were built (``genes`` / ``ensembl_ids`` of the set). Nothing
       about gene choice is decided at scoring time.
    2. Per-sample transform (``SCALING``): x = log2(TPM + 1) over the frozen
       genes present in the input; then z = (x - mean(x)) / sd(x), the mean and
       the population standard deviation (ddof=0) taken over that sample's
       genes only. This is the harmonisation rule registered for T1
       (pre-registration section 2.2 item 9(d)), applied to log2(TPM + 1).
    3. Similarity: Pearson correlation between z and each subtype's centroid,
       over the same genes. The call is the subtype with the highest
       correlation; an exact tie goes to the earlier subtype in ``subtypes``.

WHERE ONE SAMPLE'S RESULT STILL DEPENDS ON OTHER SAMPLES
    * Through the frozen centroids, and so through the training population and
      its labels. That dependence is fixed at build time and recorded in the
      set's ``meta``; it does not change with the batch being scored.
    * Not through its batch: no quantity is computed across the samples of a
      call. Scoring a sample alone or inside any batch gives identical output
      (tested in ``tests/test_centroid.py``).
    * Through the gene intersection only if the input lacks frozen genes: the
      genes used are those of the sample's own matrix, which is the same for
      every sample of one matrix. ``subtype_n_genes_used`` reports it.

HOW THIS DIFFERS FROM VERHAAK ET AL. (2010)
    Verhaak standardised each gene across patients within each platform
    ("unifiedScaled") and classified with ClaNC. Gene-wise standardisation
    across patients is a batch quantity, so it is not used here. Pearson
    correlation is invariant to shifting and scaling either vector, so step 2
    does not change the call; it is applied so that the centroids live on the
    registered T1 scale.

WHAT THIS IS NOT
    Not a validation. The centroids are means of TCGA RNA-seq profiles of
    patients carrying a published label; agreement with those labels on TCGA is
    a descriptive check (resubstitution for the all-TCGA set).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from gliostrata.profiling.verhaak import strip_ensembl_version

__all__ = [
    "SCALING",
    "MODE",
    "CENTROID_DIR",
    "CHECKSUM_FILE",
    "SET_V2TRAIN",
    "SET_ALL",
    "CentroidSet",
    "transform_sample",
    "load_centroid_set",
    "classify_single_sample",
    "centroid_set_to_json",
    "sha256_bytes",
]

SCALING = (
    "per sample: x = log2(TPM + 1) over the frozen genes present; "
    "z = (x - mean(x)) / sd(x), mean and population sd (ddof=0) over that "
    "sample's genes only (pre-registration 2.2 item 9(d)); similarity = Pearson "
    "correlation with each centroid over the same genes"
)
MODE = "single_sample"
CENTROID_DIR = Path(__file__).resolve().parents[3] / "configs" / "subtype_centroids"
CHECKSUM_FILE = "CHECKSUMS.json"
#: Used for every per-patient subtype feature and any internal evaluation on
#: the v2 test split: estimated on v2 training patients only.
SET_V2TRAIN = "verhaak840_tcga_v2train"
#: Reserved for T1: estimated on every labelled TCGA patient.
SET_ALL = "verhaak840_tcga_all"
SCHEMA_VERSION = 1


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@dataclass(frozen=True)
class CentroidSet:
    set_id: str
    genes: tuple[str, ...]
    subtypes: tuple[str, ...]
    centroids: pd.DataFrame
    scaling: str
    source_sha256: str
    meta: dict = field(default_factory=dict)
    #: Version-stripped Ensembl gene id fixed for each entry of ``genes``.
    ensembl_ids: tuple[str, ...] = ()


def transform_sample(values: np.ndarray) -> np.ndarray:
    """log2(TPM + 1), then z-scored across genes within the one sample."""
    x = np.log2(np.asarray(values, dtype=float) + 1.0)
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - x.mean()) / sd


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    den = float(np.sqrt((a * a).sum() * (b * b).sum()))
    return float((a * b).sum() / den) if den > 0 else float("nan")


def centroid_set_to_json(cs: CentroidSet) -> str:
    """Canonical serialisation: fixed key order, fixed float format."""
    doc = {
        "schema_version": SCHEMA_VERSION,
        "set_id": cs.set_id,
        "mode": MODE,
        "scaling": cs.scaling,
        "subtypes": list(cs.subtypes),
        "genes": [
            {"symbol": g, "ensembl_id": e}
            for g, e in zip(cs.genes, cs.ensembl_ids or ("",) * len(cs.genes))
        ],
        "centroids": {
            s: [float(f"{v:.10g}") for v in cs.centroids[s].to_numpy(dtype=float)]
            for s in cs.subtypes
        },
        "meta": cs.meta,
    }
    return json.dumps(doc, indent=1, sort_keys=False, ensure_ascii=True) + "\n"


def _recorded_checksum(path: Path) -> str:
    registry = path.parent / CHECKSUM_FILE
    if not registry.is_file():
        raise ValueError(f"no checksum registry {registry}: refusing an unverified centroid set")
    recorded = json.loads(registry.read_text(encoding="utf-8"))
    entry = recorded.get(path.name)
    if not entry:
        raise ValueError(f"{path.name} is not listed in {registry}: refusing an unknown centroid set")
    return str(entry)


def load_centroid_set(path) -> CentroidSet:
    """Load a frozen centroid set and verify its bytes against the registry."""
    path = Path(path)
    data = path.read_bytes()
    digest = sha256_bytes(data)
    expected = _recorded_checksum(path)
    if digest != expected:
        raise ValueError(
            f"{path.name}: SHA-256 {digest} does not match the recorded {expected}; "
            "the frozen centroid set has been altered"
        )
    doc = json.loads(data.decode("utf-8"))
    if doc.get("schema_version") != SCHEMA_VERSION or doc.get("mode") != MODE:
        raise ValueError(f"{path.name}: unexpected schema or mode")
    genes = tuple(g["symbol"] for g in doc["genes"])
    ens = tuple(g["ensembl_id"] for g in doc["genes"])
    subtypes = tuple(doc["subtypes"])
    cen = pd.DataFrame({s: doc["centroids"][s] for s in subtypes}, index=pd.Index(genes, name="gene"))
    if len(set(genes)) != len(genes) or cen.isna().to_numpy().any():
        raise ValueError(f"{path.name}: duplicate genes or missing centroid values")
    return CentroidSet(
        set_id=str(doc["set_id"]),
        genes=genes,
        subtypes=subtypes,
        centroids=cen,
        scaling=str(doc["scaling"]),
        source_sha256=digest,
        meta=dict(doc.get("meta", {})),
        ensembl_ids=ens,
    )


def _frozen_gene_matrix(
    expression: pd.DataFrame, annotation: pd.DataFrame | None, cs: CentroidSet
) -> pd.DataFrame:
    """Rows of ``expression`` for the frozen genes, chosen without any cross-sample rule."""
    if annotation is not None and cs.ensembl_ids and all(cs.ensembl_ids):
        stripped = pd.Index(strip_ensembl_version([str(i) for i in expression.index]))
        if stripped.has_duplicates:
            raise ValueError("expression index has duplicate Ensembl ids after version stripping")
        by_id = pd.Series(np.arange(len(stripped)), index=stripped)
        rows, genes = [], []
        for g, e in zip(cs.genes, cs.ensembl_ids):
            if e in by_id.index:
                rows.append(int(by_id[e]))
                genes.append(g)
        out = expression.iloc[rows].copy()
        out.index = pd.Index(genes, name="gene")
        return out
    idx = pd.Index([str(i) for i in expression.index])
    if idx.has_duplicates:
        raise ValueError("symbol-indexed expression has duplicate symbols; pass the annotation")
    out = expression.copy()
    out.index = idx
    return out.loc[[g for g in cs.genes if g in set(idx)]]


def classify_single_sample(
    expression: pd.DataFrame,
    annotation: pd.DataFrame | None,
    centroid_set: CentroidSet,
) -> pd.DataFrame:
    """Assign each sample (column) to its nearest frozen centroid, on its own.

    ``expression`` is gene x sample TPM, Ensembl-indexed when ``annotation`` is
    given (the genes are then matched by the Ensembl id frozen in the set) and
    symbol-indexed otherwise. One output row per sample, in column order.
    """
    cs = centroid_set
    sub = _frozen_gene_matrix(expression, annotation, cs)
    n_used = int(sub.shape[0])
    if n_used < 2:
        raise ValueError(f"only {n_used} of the {len(cs.genes)} frozen genes are in the input")
    cen = cs.centroids.loc[list(sub.index)].to_numpy(dtype=float)
    vals = sub.to_numpy(dtype=float)
    if not np.isfinite(vals).all():
        raise ValueError("expression has non-finite values on the frozen genes")
    records = []
    for j, sample in enumerate(sub.columns):
        z = transform_sample(vals[:, j])
        corr = [_pearson(z, cen[:, k]) for k in range(cen.shape[1])]
        arr = np.array(corr, dtype=float)
        # A sample with no finite correlation (zero variance over the frozen
        # genes) has no call. It is never given the first subtype by default.
        call = cs.subtypes[int(np.nanargmax(arr))] if np.isfinite(arr).any() else None
        rec = {"patient_id": str(sample), "subtype": call}
        for s, c in zip(cs.subtypes, corr):
            rec[f"subtype_corr_{s}"] = c
        records.append(rec)
    cols = ["patient_id", "subtype"] + [f"subtype_corr_{s}" for s in cs.subtypes]
    out = pd.DataFrame.from_records(records, columns=cols)
    out["subtype_mode"] = MODE
    out["subtype_centroid_set"] = cs.set_id
    out["subtype_centroid_sha256"] = cs.source_sha256
    out["subtype_n_genes_used"] = n_used
    return out
