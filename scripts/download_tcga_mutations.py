"""
Download TCGA-GBM open-access somatic mutation data (Masked Somatic Mutation,
MAF) from the GDC REST API and call IDH1 / IDH2 hotspot status per patient.

WHY THIS EXISTS
    TCGA-GBM is the development cohort and CGGA is the external validation
    cohort, and they barely share features. CGGA carries IDH status for
    326/333 patients; the TCGA clinical ingest carries none. This script
    closes that gap on the TCGA side (the MGMT half is
    scripts/download_tcga_methylation.py).

WHAT IS CALLED, EXACTLY
    IDH status is defined by missense / in-frame variants at two codons:
        IDH1 R132  (R132H dominant in glioma; also C, G, S, L and others)
        IDH2 R172  (R172K dominant; also M, W, G and others)
    A variant is a hotspot hit only if ALL hold:
      * the gene is IDH1 or IDH2;
      * the annotated reference residue is arginine (R) at codon 132 / 172;
      * the variant class is Missense_Mutation, In_Frame_Ins or In_Frame_Del.
    Nonsense, frameshift, splice-site and silent changes at those codons are
    NOT called: they do not create the neomorphic enzyme activity that
    defines IDH-mutant glioma, and calling them would overstate the rate.

THREE STATES, NOT TWO -- the point of the script
      1      IDH-mutant      a hotspot variant was found in a tumour sample
      0      IDH-wildtype    the tumour WAS sequenced, no hotspot was found
      <NA>   no mutation data  no MAF exists for the patient
    The third state is NOT wildtype. Coercing it to 0 silently shrinks the
    mutant fraction and invents patients who were never tested. It is written
    as a nullable value and counted separately in every report.

WHO 2021 IMPLICATION (reported, not resolved)
    Under WHO 2021 an IDH-mutant grade IV astrocytic tumour is "astrocytoma,
    IDH-mutant, CNS WHO grade 4", NOT glioblastoma. Any IDH-mutant case here
    is therefore possibly mis-grouped by the clinical ingest. This script
    counts and lists those cases; removing them is a pre-registration decision
    for the author, not something an ingestion script may do on its own.

OUTPUTS (under --out, default data/raw/tcga_gbm/mutations)
    file_manifest.parquet   one row per GDC file, with parsed barcodes
    idh_calls.parquet       per-patient calls, INCLUDING patients with no data
    call_report.json        every count printed to stdout, machine-readable
    _cache/<file_id>.maf.gz download cache, separate from the outputs

    clinical_cohort.parquet is NEVER modified. The calls are a standalone
    table to be joined in a separate, reviewed step.

Usage:
    python scripts/download_tcga_mutations.py --dry-run
    python scripts/download_tcga_mutations.py --limit 5
    python scripts/download_tcga_mutations.py
"""

import argparse
import gzip
import io
import json
import re
import tarfile
import time
from pathlib import Path

import pandas as pd
import requests

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"

#: Project this script ingests. TCGA-GBM only; not a general-purpose tool.
PROJECT_ID = "TCGA-GBM"

#: Columns this script reads from a GDC masked MAF. Any absent -> ValueError.
MAF_REQUIRED_COLUMNS = [
    "Hugo_Symbol",
    "Variant_Classification",
    "Tumor_Sample_Barcode",
    "HGVSp_Short",
    "Protein_position",
    "Amino_acids",
]

#: Gene -> codon that defines IDH-mutant glioma (WHO 2021).
IDH_HOTSPOT_CODONS = {"IDH1": 132, "IDH2": 172}

#: The reference residue at both hotspots is arginine. A record whose
#: annotated reference is not R is on a different numbering, not the hotspot.
HOTSPOT_REFERENCE_RESIDUE = "R"

#: Variant classes that create the neomorphic enzyme. Everything else is not.
HOTSPOT_VARIANT_CLASSES = ("Missense_Mutation", "In_Frame_Ins", "In_Frame_Del")

#: TCGA sample-type codes (barcode positions 14-15).
SAMPLE_TYPE_CODES = {
    "01": "Primary Tumor",
    "02": "Recurrent Tumor",
    "05": "Additional New Primary",
    "06": "Metastatic",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
}

#: Calls are made on the primary tumour wherever one exists.
DEFAULT_SAMPLE_TYPE_CODES = ("01",)

#: Preference order used to pick the tumour aliquot among all aliquots linked
#: to a file (lower = preferred). Normal samples (10, 11) never rank here.
TUMOUR_RANK = {"01": 0, "02": 1, "05": 2, "06": 3}

GDC_FILE_FIELDS = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "analysis.workflow_type",
    "cases.submitter_id",
    "cases.case_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.samples.portions.analytes.aliquots.submitter_id",
]

#: Output columns of idh_calls.parquet, in order.
CALL_COLUMNS = [
    "patient_id",
    "idh_mutated",
    "idh_variants",
    "n_idh_hotspot_variants",
    "call_basis",
    "n_somatic_variants",
    "sample_barcode",
    "aliquot_barcode",
    "file_id",
]

#: Value of `call_basis` for each of the three states.
BASIS_MUTANT = "hotspot_variant"
BASIS_WILDTYPE = "sequenced_no_hotspot"
BASIS_MISSING = "no_mutation_data"

DEFAULT_OUT = Path("data/raw/tcga_gbm/mutations")
DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")

_BARCODE_RE = re.compile(
    r"^(?P<patient_id>TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})"
    r"(?:-(?P<sample>(?P<sample_type_code>\d{2})(?P<vial>[A-Z]?)))?"
    r"(?P<rest>-.*)?$"
)

# p.R132H, R132H (no prefix), p.R132_P135delinsHQ, p.R132del, p.R132dup.
_HGVSP_RE = re.compile(r"^p\.([A-Z])(\d+)(.*)$")

# --------------------------------------------------------------------------
# GDC metadata
# --------------------------------------------------------------------------

def build_file_filters(project_id: str = PROJECT_ID) -> dict:
    """GDC `files` filter for a project's open masked somatic MAFs."""
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": [project_id]}},
            {"op": "in", "content": {"field": "data_category",
                                     "value": ["Simple Nucleotide Variation"]}},
            {"op": "in", "content": {"field": "data_type",
                                     "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }


def fetch_file_metadata(project_id: str = PROJECT_ID, page_size: int = 500) -> list[dict]:
    """Fetch all matching file records from the GDC `files` endpoint (paged)."""
    filters = build_file_filters(project_id)
    hits: list[dict] = []
    offset = 0
    total = None
    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(GDC_FILE_FIELDS),
            "format": "JSON",
            "size": str(page_size),
            "from": str(offset),
            "sort": "file_id:asc",
        }
        response = requests.get(GDC_FILES_URL, params=params, timeout=300)
        if response.status_code != 200:
            raise RuntimeError(
                f"GDC API error {response.status_code}: {response.text[:500]}"
            )
        payload = response.json()["data"]
        total = payload["pagination"]["total"]
        hits.extend(payload["hits"])
        offset += page_size
        if offset >= total:
            break
    if len(hits) != total:
        raise RuntimeError(f"paging mismatch: got {len(hits)} of {total} files")
    return hits


def linked_aliquots(case: dict) -> list[dict]:
    """Every aliquot GDC links to a file, across ALL of the case's samples.

    A masked MAF is linked to the tumour aliquot AND to the matched normal
    (blood or solid tissue) that was used to make it. Which sample comes first
    in the metadata is arbitrary -- in TCGA-GBM the normal sample is first for
    239 of 464 MAF files -- so ``samples[0]`` is NOT a safe way to learn what
    kind of tissue a file describes.
    """
    out: list[dict] = []
    for sample in case.get("samples") or []:
        for portion in sample.get("portions") or []:
            for analyte in portion.get("analytes") or []:
                for aliquot in analyte.get("aliquots") or []:
                    if aliquot.get("submitter_id"):
                        out.append({
                            "aliquot_barcode": aliquot["submitter_id"],
                            "sample_barcode": sample.get("submitter_id"),
                            "gdc_sample_type": sample.get("sample_type"),
                        })
    return out


def choose_tumour_aliquot(aliquots: list[dict]) -> dict:
    """Pick the aliquot that describes the tumour, deterministically.

    Rank by tumour type (primary 01, then recurrent 02, then other tumours),
    then by barcode, so reruns are byte-identical. If no tumour aliquot is
    linked, the first barcode in sorted order is kept and its sample-type code
    is what decides whether the file is used at all.
    """
    ranked = []
    for item in aliquots:
        code = parse_barcode(item["aliquot_barcode"])["sample_type_code"]
        ranked.append((TUMOUR_RANK.get(code, 9), str(item["aliquot_barcode"]), item, code))
    if not ranked:
        return {"aliquot_barcode": None, "sample_barcode": None,
                "gdc_sample_type": None, "sample_type_code": None, "vial": None}
    ranked.sort(key=lambda t: (t[0], t[1]))
    _rank, _barcode, item, code = ranked[0]
    return {
        "aliquot_barcode": item["aliquot_barcode"],
        "sample_barcode": item["sample_barcode"],
        "gdc_sample_type": item["gdc_sample_type"],
        "sample_type_code": code,
        "vial": parse_barcode(item["aliquot_barcode"])["vial"],
    }


def flatten_file_record(hit: dict) -> dict:
    """Flatten one GDC file record to a single-row dict with a parsed barcode.

    The patient id is the 12-character prefix of the tumour aliquot barcode,
    in the same form as clinical_cohort.parquet. See
    :func:`linked_aliquots` for why the tumour aliquot is chosen explicitly
    rather than taken as ``samples[0]``.
    """
    cases = hit.get("cases") or []
    case = cases[0] if cases else {}
    chosen = choose_tumour_aliquot(linked_aliquots(case))
    parsed = parse_barcode(chosen["aliquot_barcode"] or case.get("submitter_id"))

    return {
        "file_id": hit.get("file_id"),
        "file_name": hit.get("file_name"),
        "file_size": hit.get("file_size"),
        "md5sum": hit.get("md5sum"),
        "workflow_type": (hit.get("analysis") or {}).get("workflow_type"),
        "case_id": case.get("case_id"),
        "case_submitter_id": case.get("submitter_id"),
        "sample_barcode": chosen["sample_barcode"],
        "aliquot_barcode": chosen["aliquot_barcode"],
        "gdc_sample_type": chosen["gdc_sample_type"],
        "patient_id": parsed["patient_id"],
        "sample_type_code": chosen["sample_type_code"] or parsed["sample_type_code"],
        "vial": parsed["vial"],
    }

# --------------------------------------------------------------------------
# Barcode parsing (same rules as download_tcga_expression.py)
# --------------------------------------------------------------------------

def parse_barcode(barcode) -> dict:
    """Parse a TCGA barcode into patient_id / sample_type_code / vial.

    "TCGA-12-0819-01A-01R-1849-01" -> patient "TCGA-12-0819", code "01".
    Unparseable or missing input yields all-None, never a partial guess.
    """
    none = {"patient_id": None, "sample_type_code": None, "vial": None}
    if not isinstance(barcode, str):
        return dict(none)
    match = _BARCODE_RE.match(barcode.strip().upper())
    if match is None:
        return dict(none)
    return {
        "patient_id": match.group("patient_id"),
        "sample_type_code": match.group("sample_type_code"),
        "vial": match.group("vial") or None,
    }


def sample_type_label(code) -> str:
    """Human-readable label for a TCGA sample-type code."""
    if not isinstance(code, str):
        return "Unknown"
    return SAMPLE_TYPE_CODES.get(code, f"Other ({code})")


def select_sample_types(
    manifest: pd.DataFrame, keep_codes=DEFAULT_SAMPLE_TYPE_CODES
) -> tuple[pd.DataFrame, dict]:
    """Keep only rows whose sample_type_code is in `keep_codes`."""
    keep = list(keep_codes)
    breakdown: dict[str, int] = {}
    for code, count in manifest["sample_type_code"].value_counts(dropna=False).items():
        key = code if isinstance(code, str) else "unparsed"
        breakdown[key] = int(count)
    filtered = manifest[manifest["sample_type_code"].isin(keep)].copy()
    report = {
        "input_files": int(len(manifest)),
        "breakdown": dict(sorted(breakdown.items())),
        "kept_codes": keep,
        "kept_files": int(len(filtered)),
        "dropped_files": int(len(manifest) - len(filtered)),
    }
    return filtered, report


def resolve_duplicate_aliquots(manifest: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Reduce the manifest to ONE file per patient, deterministically.

    Sort by (sample_type_code, aliquot_barcode, file_id) and keep the first,
    so the lowest-numbered vial/portion wins and reruns are byte-identical.
    """
    if manifest.empty:
        return manifest.copy(), {
            "input_files": 0,
            "patients": 0,
            "patients_with_duplicates": 0,
            "files_dropped_as_duplicate": 0,
        }
    ordered = manifest.sort_values(
        ["patient_id", "sample_type_code", "aliquot_barcode", "file_id"],
        kind="mergesort",
        na_position="last",
    )
    counts = ordered.groupby("patient_id", dropna=False).size()
    deduped = ordered.drop_duplicates(subset=["patient_id"], keep="first").copy()
    deduped = deduped.sort_values("patient_id", kind="mergesort").reset_index(drop=True)
    report = {
        "input_files": int(len(manifest)),
        "patients": int(len(deduped)),
        "patients_with_duplicates": int((counts > 1).sum()),
        "files_dropped_as_duplicate": int(len(manifest) - len(deduped)),
    }
    return deduped, report

# --------------------------------------------------------------------------
# Download
# --------------------------------------------------------------------------

def cache_path_for(cache_dir: Path, file_id: str) -> Path:
    """Cache location for one downloaded MAF."""
    return Path(cache_dir) / f"{file_id}.maf.gz"


def pending_file_ids(file_ids, cache_dir: Path) -> list[str]:
    """File ids not already cached on disk (resume support)."""
    pending = []
    for file_id in file_ids:
        path = cache_path_for(cache_dir, file_id)
        if path.is_file() and path.stat().st_size > 0:
            continue
        pending.append(file_id)
    return pending


def _store_response(file_ids: list[str], content: bytes, cache_dir: Path) -> int:
    """Write a /data response into the cache (bare file, or a .tar.gz)."""
    if len(file_ids) == 1:
        cache_path_for(cache_dir, file_ids[0]).write_bytes(content)
        return 1
    written = 0
    wanted = set(file_ids)
    with tarfile.open(fileobj=io.BytesIO(content), mode="r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            parts = Path(member.name).parts
            if len(parts) < 2 or parts[0] not in wanted:
                continue  # MANIFEST.txt at the archive root, or unknown ids
            handle = tar.extractfile(member)
            if handle is None:
                continue
            cache_path_for(cache_dir, parts[0]).write_bytes(handle.read())
            written += 1
    return written


def download_batch(
    file_ids: list[str],
    cache_dir: Path,
    timeout: int = 600,
    max_attempts: int = 3,
    retry_wait: float = 2.0,
) -> int:
    """Download a batch via POST /data and write each member to the cache.

    The endpoint answers a single-id request with the bare file and a
    multi-id request with a .tar.gz whose members are `<file_id>/<name>`.
    (GET is not allowed on /data; it returns 405.) GDC intermittently answers
    500 to a multi-id POST, so the request is retried ``max_attempts`` times
    before the error is raised to the caller, which then isolates bad files.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, int(max_attempts) + 1):
        try:
            response = requests.post(
                GDC_DATA_URL,
                json={"ids": list(file_ids)},
                headers={"Content-Type": "application/json"},
                timeout=timeout,
            )
            if response.status_code == 200:
                return _store_response(file_ids, response.content, cache_dir)
            last_error = RuntimeError(
                f"GDC data error {response.status_code}: {response.text[:500]}"
            )
        except requests.RequestException as exc:
            last_error = RuntimeError(f"GDC request failed: {exc}")
        if attempt < max_attempts:
            print(f"  attempt {attempt}/{max_attempts} failed ({last_error}); "
                  f"retrying in {retry_wait:.0f}s", flush=True)
            time.sleep(retry_wait)
    assert last_error is not None
    raise last_error


def download_files(file_ids: list[str], cache_dir: Path, batch_size: int = 50) -> dict:
    """Download every missing file, batching requests and showing progress.

    A batch that keeps failing is retried ONE FILE AT A TIME, so a single bad
    file id costs that file and not the whole run; ids that still fail are
    returned under ``failed`` and reported, never silently dropped.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    todo = pending_file_ids(file_ids, cache_dir)
    skipped = len(file_ids) - len(todo)
    print(
        f"Download: {len(file_ids)} files requested -> "
        f"{skipped} already cached, {len(todo)} to fetch"
    )
    started = time.time()
    written = 0
    failed: list[str] = []
    for start in range(0, len(todo), batch_size):
        batch = todo[start:start + batch_size]
        try:
            written += download_batch(batch, cache_dir)
        except RuntimeError as exc:
            print(f"  batch failed ({exc}); retrying its files one at a time", flush=True)
            for file_id in batch:
                try:
                    written += download_batch([file_id], cache_dir)
                except RuntimeError as single_error:
                    failed.append(file_id)
                    print(f"  FAILED {file_id}: {single_error}", flush=True)
        done = min(start + batch_size, len(todo))
        print(
            f"  [{done}/{len(todo)}] fetched -> {time.time() - started:6.1f}s elapsed",
            flush=True,
        )
    if failed:
        print(f"WARNING: {len(failed)} file(s) could not be downloaded and are "
              "reported as having no data: " + ", ".join(failed[:10])
              + (" ..." if len(failed) > 10 else ""))
    return {
        "requested": len(file_ids),
        "skipped_cached": skipped,
        "downloaded": written,
        "failed": failed,
        "elapsed_seconds": round(time.time() - started, 1),
    }


# --------------------------------------------------------------------------
# MAF parsing
# --------------------------------------------------------------------------

def read_maf(path) -> pd.DataFrame:
    """Read one masked MAF (plain or gzipped) into a DataFrame.

    Comment lines starting with ``#`` are dropped. A header-only file (a
    sample that was sequenced and had NO somatic variants) returns an empty
    DataFrame WITH the MAF columns -- that is a valid wildtype sample, not an
    error, and must never be confused with a missing file.
    """
    path = Path(path)
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    lines = raw.decode("utf-8", errors="replace").splitlines()
    header_idx = next(
        (i for i, line in enumerate(lines) if line.startswith("Hugo_Symbol")), None
    )
    if header_idx is None:
        raise ValueError(f"{path}: no MAF header line (Hugo_Symbol) found")
    frame = pd.read_csv(io.StringIO("\n".join(lines[header_idx:])), sep="\t", low_memory=False)
    missing = [c for c in MAF_REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"{path}: MAF is missing required column(s): {missing}")
    return frame

# --------------------------------------------------------------------------
# IDH hotspot calling
# --------------------------------------------------------------------------

def parse_protein_change(hgvsp_short, protein_position=None, amino_acids=None) -> tuple:
    """Return ``(ref_residue, codon, alt_text)`` for one variant, or Nones.

    Primary source is ``HGVSp_Short`` (``p.R132H``, ``p.R132_P135delinsHQ``).
    Fallback is the maf2vcf-style pair ``Protein_position`` / ``Amino_acids``
    (``132`` / ``R/H``). Anything unparseable returns three Nones and is
    therefore never called -- an unparsed variant cannot become a hit.
    """
    if isinstance(hgvsp_short, str) and hgvsp_short.strip():
        value = hgvsp_short.strip()
        if not value.startswith("p."):
            value = "p." + value
        match = _HGVSP_RE.match(value)
        if match and match.group(3):
            # An empty alt ("p.R132") is not a substitution; do not call it.
            return match.group(1), int(match.group(2)), match.group(3)
    if (
        isinstance(amino_acids, str)
        and "/" in amino_acids
        and isinstance(protein_position, str)
    ):
        try:
            codon = int(protein_position.split("-")[0].split("/")[0])
        except ValueError:
            return None, None, None
        ref, alt = amino_acids.split("/", 1)
        if len(ref) == 1 and ref.isalpha() and alt:
            return ref, codon, alt
    return None, None, None


def is_idh_hotspot_variant(row) -> bool:
    """True if one MAF row is an IDH1 R132 / IDH2 R172 hotspot variant."""
    gene = str(row.get("Hugo_Symbol"))
    if gene not in IDH_HOTSPOT_CODONS:
        return False
    if str(row.get("Variant_Classification")) not in HOTSPOT_VARIANT_CLASSES:
        return False
    ref, codon, _alt = parse_protein_change(
        row.get("HGVSp_Short"), row.get("Protein_position"), row.get("Amino_acids")
    )
    if codon != IDH_HOTSPOT_CODONS[gene]:
        return False
    return ref == HOTSPOT_REFERENCE_RESIDUE


def hotspot_label(row) -> str:
    """``IDH1 p.R132H``-style label for a row that IS a hotspot hit."""
    gene = str(row.get("Hugo_Symbol"))
    ref, codon, alt = parse_protein_change(
        row.get("HGVSp_Short"), row.get("Protein_position"), row.get("Amino_acids")
    )
    return f"{gene} p.{ref}{codon}{alt}"


def call_idh_for_maf(maf: pd.DataFrame) -> dict:
    """Summarise one sample's MAF into a hotspot call.

    An empty MAF (header only) is a sequenced sample with no somatic variants:
    zero hits, and NOT an error.
    """
    n_variants = int(len(maf))
    if n_variants == 0:
        return {"n_somatic_variants": 0, "hotspot_variants": []}
    gene = maf["Hugo_Symbol"].astype(str)
    cls = maf["Variant_Classification"].astype(str)
    candidates = maf[gene.isin(IDH_HOTSPOT_CODONS) & cls.isin(HOTSPOT_VARIANT_CLASSES)]
    labels = set()
    for _i, row in candidates.iterrows():
        if is_idh_hotspot_variant(row):
            labels.add(hotspot_label(row))
    return {"n_somatic_variants": n_variants, "hotspot_variants": sorted(labels)}


def maf_tumour_patient(maf: pd.DataFrame) -> str | None:
    """Patient id encoded in the MAF's own ``Tumor_Sample_Barcode``.

    Returns None for a header-only MAF (no barcodes) or when the file mixes
    more than one patient (ambiguous, so not usable as a cross-check).
    """
    if len(maf) == 0:
        return None
    barcodes = maf["Tumor_Sample_Barcode"].dropna().astype(str).unique()
    patients = {parse_barcode(b)["patient_id"] for b in barcodes} - {None}
    return next(iter(patients)) if len(patients) == 1 else None


def _call_row(pid: str, labels: list[str], n_variants, record=None) -> dict:
    record = record or {}
    mutant = bool(labels)
    return {
        "patient_id": pid,
        "idh_mutated": 1 if mutant else 0,
        "idh_variants": ";".join(labels),
        "n_idh_hotspot_variants": len(labels),
        "call_basis": BASIS_MUTANT if mutant else BASIS_WILDTYPE,
        "n_somatic_variants": n_variants,
        "sample_barcode": record.get("sample_barcode"),
        "aliquot_barcode": record.get("aliquot_barcode"),
        "file_id": record.get("file_id"),
    }


def build_calls(sample_records: list[dict], cohort_patient_ids=None) -> pd.DataFrame:
    """Assemble per-patient calls, INCLUDING patients with no mutation data.

    ``sample_records`` holds one dict per tumour sample with ``patient_id``,
    ``hotspot_variants``, ``n_somatic_variants`` and sample provenance. A
    patient absent from every record is MISSING (nullable), never wildtype.
    If a patient has more than one sample the calls are merged: mutant if any
    sample carries a hotspot variant.
    """
    by_patient: dict[str, dict] = {}
    for record in sample_records:
        pid = str(record["patient_id"])
        labels = sorted(set(record.get("hotspot_variants") or []))
        n_variants = int(record.get("n_somatic_variants") or 0)
        prev = by_patient.get(pid)
        if prev is not None:
            labels = sorted(
                set(filter(None, str(prev["idh_variants"]).split(";"))) | set(labels)
            )
            n_variants = max(int(prev["n_somatic_variants"]), n_variants)
        by_patient[pid] = _call_row(pid, labels, n_variants, record)

    if cohort_patient_ids is not None:
        for pid in cohort_patient_ids:
            pid = str(pid)
            if pid in by_patient:
                continue
            by_patient[pid] = {
                "patient_id": pid,
                "idh_mutated": pd.NA,
                "idh_variants": None,
                "n_idh_hotspot_variants": pd.NA,
                "call_basis": BASIS_MISSING,
                "n_somatic_variants": pd.NA,
                "sample_barcode": None,
                "aliquot_barcode": None,
                "file_id": None,
            }

    if not by_patient:
        return pd.DataFrame({c: pd.Series(dtype="object") for c in CALL_COLUMNS})
    out = pd.DataFrame([by_patient[p] for p in sorted(by_patient)])[CALL_COLUMNS]
    for col in ("idh_mutated", "n_idh_hotspot_variants", "n_somatic_variants"):
        out[col] = out[col].astype("Int64")
    return out.reset_index(drop=True)

# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

#: A mutant fraction above this among tested patients is treated as a red flag
#: to investigate, not a result to report: the published rate in primary,
#: de novo glioblastoma is single-digit percent.
RED_FLAG_MUTANT_FRACTION = 0.10


def summarize_calls(calls: pd.DataFrame) -> dict:
    """The counts behind every number printed to stdout."""
    n = int(len(calls))
    mutant = int((calls["idh_mutated"] == 1).sum())
    wildtype = int((calls["idh_mutated"] == 0).sum())
    missing = int(calls["idh_mutated"].isna().sum())
    tested = mutant + wildtype
    variants = sorted(
        {v for value in calls["idh_variants"].dropna() for v in str(value).split(";") if v}
    )
    frac_tested = (mutant / tested) if tested else None
    return {
        "patients_in_table": n,
        "patients_with_mutation_data": tested,
        "idh_mutant": mutant,
        "idh_wildtype": wildtype,
        "no_mutation_data": missing,
        "mutant_fraction_of_tested": frac_tested,
        "mutant_fraction_of_all_rows": (mutant / n) if n else None,
        "distinct_hotspot_variants": variants,
        "mutant_patient_ids": sorted(
            calls.loc[calls["idh_mutated"] == 1, "patient_id"].astype(str).tolist()
        ),
        "red_flag_high_mutant_rate": bool(
            frac_tested is not None and frac_tested > RED_FLAG_MUTANT_FRACTION
        ),
    }


def coverage(calls: pd.DataFrame, patient_ids, label: str) -> dict:
    """How many of a given patient set gain an IDH call."""
    ids = {str(p) for p in patient_ids}
    sub = calls[calls["patient_id"].astype(str).isin(ids)]
    has = sub["idh_mutated"].notna()
    return {
        "set": label,
        "patients": len(ids),
        "with_idh_call": int(has.sum()),
        "idh_mutant": int((sub["idh_mutated"] == 1).sum()),
        "idh_wildtype": int((sub["idh_mutated"] == 0).sum()),
        "still_missing": int(len(ids) - has.sum()),
    }


def parquet_column_names(path) -> list[str]:
    """Column names of a parquet file WITHOUT loading its data.

    ``pd.read_parquet(path, columns=[])`` looks equivalent but returns ZERO
    columns, which silently reports an empty patient set. The schema is read
    instead; the pandas index placeholder column is dropped if present.
    """
    try:
        import pyarrow.parquet as pq

        names = list(pq.ParquetFile(str(path)).schema.names)
    except ImportError:  # pragma: no cover - pyarrow ships with pandas here
        names = [str(c) for c in pd.read_parquet(path).columns]
    return [n for n in names if not n.startswith("__index_level_")]


def expression_patient_ids(path) -> tuple[list[str], list[str]]:
    """Patient ids that are the sample columns of a gene x sample matrix.

    The matrix is written with ``index=False``, so its first column
    (``gene_id``) is a gene and not a patient. Treating every column as a
    patient overstates the expression set by one and makes "cohort AND
    expression" disagree with the expression ingest's own report.

    Returns ``(patient_ids, non_sample_columns)``. A file with no barcode
    column at all is the wrong file and is refused.
    """
    ids: list[str] = []
    other: list[str] = []
    for name in parquet_column_names(path):
        text = str(name).strip()
        if _BARCODE_RE.match(text.upper()):
            ids.append(text)
        else:
            other.append(text)
    if not ids:
        raise RuntimeError(
            f"{path}: none of its {len(other)} column(s) is a TCGA barcode -- "
            "this is not a gene x sample matrix"
        )
    return ids, other


def _jsonable(obj):
    """Recursively make a report structure strict-JSON safe (no NaN/inf)."""
    import math

    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        return _jsonable(obj.item())
    return obj

def _pct(num, den) -> str:
    return f"{100.0 * num / den:5.1f}%" if den else "   n/a"


def print_report(summary: dict, coverage_rows: list[dict]) -> None:
    """ASCII-only summary (cp1252 consoles choke on non-ASCII)."""
    tested = summary["patients_with_mutation_data"]
    print("\nIDH1 R132 / IDH2 R172 status -- three states, kept apart:")
    print(f"  IDH-mutant          {summary['idh_mutant']:5d}   "
          f"{_pct(summary['idh_mutant'], tested)} of tested")
    print(f"  IDH-wildtype        {summary['idh_wildtype']:5d}   "
          f"{_pct(summary['idh_wildtype'], tested)} of tested")
    print(f"  no mutation data    {summary['no_mutation_data']:5d}   "
          f"{_pct(summary['no_mutation_data'], summary['patients_in_table'])} of table "
          f"(NOT wildtype; written as <NA>)")
    variants = ", ".join(summary["distinct_hotspot_variants"]) or "none"
    print(f"  distinct hotspot variants: {variants}")

    print("\nCoverage -- how many patients gain a call:")
    print(f"  {'patient set':28}{'n':>6}{'with call':>11}{'mutant':>8}{'missing':>9}")
    for row in coverage_rows:
        print(f"  {row['set']:28}{row['patients']:6d}{row['with_idh_call']:11d}"
              f"{row['idh_mutant']:8d}{row['still_missing']:9d}")

    ids = summary["mutant_patient_ids"]
    if ids:
        shown = ", ".join(ids[:40]) + (" ..." if len(ids) > 40 else "")
        print(f"\nWHO 2021 FLAG -- {len(ids)} IDH-mutant case(s): {shown}")
        print("  An IDH-mutant grade IV astrocytic tumour is 'astrocytoma, IDH-mutant,")
        print("  CNS WHO grade 4', not glioblastoma. These may be mis-grouped by the")
        print("  clinical ingest. Whether to exclude them is a pre-registration")
        print("  decision for the author; this script does not remove them.")
    if summary["red_flag_high_mutant_rate"]:
        print("\nRED FLAG -- mutant fraction exceeds "
              f"{100.0 * RED_FLAG_MUTANT_FRACTION:.0f}% of tested patients. Primary GBM is")
        print("  expected at single-digit percent. Investigate (sample-type mix,")
        print("  parsing, cohort composition) before using these calls.")

# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--project", default=PROJECT_ID)
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Output directory")
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT),
                        help="clinical_cohort.parquet; its patients are listed even "
                             "when no MAF exists for them")
    parser.add_argument("--no-cohort", action="store_true",
                        help="do not read the clinical cohort")
    parser.add_argument("--expression", default="data/raw/tcga_gbm/expression/tpm.parquet",
                        help="gene x sample matrix; its columns are the "
                             "expression-covered patients (coverage report only)")
    parser.add_argument("--limit", type=int, default=None,
                        help="use only the first N patients (after deduplication)")
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--all-sample-types", action="store_true",
                        help="keep recurrent/normal samples too (default: primary only)")
    parser.add_argument("--dry-run", action="store_true", help="report metadata only")
    args = parser.parse_args()

    print("TCGA-GBM mutation ingest -- IDH1 R132 / IDH2 R172 calling")
    print("  Downloads and calls only. Fits nothing and joins nothing;")
    print("  clinical_cohort.parquet is not modified.")

    out_dir = Path(args.out)
    cache_dir = out_dir / "_cache"

    hits = fetch_file_metadata(args.project)
    manifest = pd.DataFrame([flatten_file_record(h) for h in hits])
    manifest = manifest.sort_values("file_id", kind="mergesort").reset_index(drop=True)
    if manifest["patient_id"].isna().any():
        bad = int(manifest["patient_id"].isna().sum())
        raise RuntimeError(f"{bad} file record(s) have an unparseable TCGA barcode")
    total_mb = float(manifest["file_size"].sum()) / 1e6
    print(f"\nGDC metadata: {len(manifest)} open masked MAF files, "
          f"{manifest['patient_id'].nunique()} patients, {total_mb:.1f} MB")
    for code, count in manifest["sample_type_code"].value_counts(dropna=False).items():
        print(f"  sample type {code}: {count:4d} files  ({sample_type_label(code)})")
    if args.dry_run:
        print("\n--dry-run -> nothing downloaded, nothing called")
        return

    if args.all_sample_types:
        kept, sample_report = manifest.copy(), {"kept_codes": "all"}
    else:
        kept, sample_report = select_sample_types(manifest)
    deduped, dedup_report = resolve_duplicate_aliquots(kept)
    print(f"\nPrimary-tumour files kept: {len(kept)} (dropped {len(manifest) - len(kept)})")
    print(f"One file per patient: {dedup_report['patients']} patients "
          f"({dedup_report['files_dropped_as_duplicate']} duplicate files dropped)")
    if args.limit is not None:
        deduped = deduped.head(int(args.limit)).copy()
        print(f"--limit {args.limit}: calling the first {len(deduped)} patient(s) only")

    download_report = download_files(
        deduped["file_id"].tolist(), cache_dir, batch_size=int(args.batch_size)
    )

    records: list[dict] = []
    unavailable: list[str] = []
    for row in deduped.itertuples(index=False):
        path = cache_path_for(cache_dir, row.file_id)
        if not path.is_file() or path.stat().st_size == 0:
            # The download failed. "No mutation data" is then exactly true for
            # this patient; count it below rather than crashing the run.
            unavailable.append(str(row.file_id))
            continue
        maf = read_maf(path)
        inside = maf_tumour_patient(maf)
        if inside is not None and inside != str(row.patient_id):
            raise RuntimeError(
                f"{path}: the MAF's own Tumor_Sample_Barcode belongs to patient "
                f"{inside} but the GDC metadata says {row.patient_id}"
            )
        call = call_idh_for_maf(maf)
        records.append({
            **call,
            "patient_id": row.patient_id,
            "sample_barcode": row.sample_barcode,
            "aliquot_barcode": row.aliquot_barcode,
            "file_id": row.file_id,
        })

    cohort_ids = None
    if not args.no_cohort and Path(args.cohort).exists():
        cohort_ids = pd.read_parquet(args.cohort, columns=["patient_id"])["patient_id"].tolist()
        print(f"\nClinical cohort read: {len(cohort_ids)} patients from {args.cohort}")
    elif not args.no_cohort:
        print(f"\nClinical cohort not found at {args.cohort}: "
              "only patients with a MAF are listed")

    calls = build_calls(records, cohort_ids)
    summary = summarize_calls(calls)

    coverage_rows: list[dict] = []
    expression_ignored: list[str] = []
    if cohort_ids is not None:
        coverage_rows.append(coverage(calls, cohort_ids, "clinical cohort"))
    if args.expression and Path(args.expression).exists():
        expr_ids, expression_ignored = expression_patient_ids(args.expression)
        if expression_ignored:
            print(f"\nExpression matrix: {len(expression_ignored)} non-sample column(s) "
                  f"not counted as patients: {', '.join(expression_ignored[:5])}")
        coverage_rows.append(coverage(calls, expr_ids, "with expression"))
        if cohort_ids is not None:
            both = set(map(str, cohort_ids)) & set(expr_ids)
            coverage_rows.append(coverage(calls, both, "cohort AND expression"))

    out_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_parquet(out_dir / "file_manifest.parquet", index=False)
    calls.to_parquet(out_dir / "idh_calls.parquet", index=False)
    report = {
        "project": args.project,
        "gdc_files": int(len(manifest)),
        "gdc_patients": int(manifest["patient_id"].nunique()),
        "total_mb": round(total_mb, 1),
        "hotspot_definition": {
            "codons": IDH_HOTSPOT_CODONS,
            "reference_residue": HOTSPOT_REFERENCE_RESIDUE,
            "variant_classes": list(HOTSPOT_VARIANT_CLASSES),
        },
        "sample_type_filter": sample_report,
        "duplicate_resolution": dedup_report,
        "download": download_report,
        "files_unavailable_after_download": unavailable,
        "summary": summary,
        "coverage": coverage_rows,
        "expression_non_sample_columns": expression_ignored,
        "limited_to_first_n_patients": args.limit,
        "who_2021_note": (
            "An IDH-mutant grade IV astrocytic tumour is astrocytoma, IDH-mutant, "
            "CNS WHO grade 4 -- not glioblastoma. Mutant cases are listed, not removed."
        ),
    }
    (out_dir / "call_report.json").write_text(
        json.dumps(_jsonable(report), indent=2), encoding="ascii"
    )

    print_report(summary, coverage_rows)
    print(f"\nWrote: {out_dir / 'file_manifest.parquet'}")
    print(f"      {out_dir / 'idh_calls.parquet'}")
    print(f"      {out_dir / 'call_report.json'}")


if __name__ == "__main__":
    main()
