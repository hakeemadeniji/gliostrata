"""
Download TCGA-GBM bulk RNA-seq gene expression (STAR - Counts) from the GDC
REST API and assemble an analysis-ready gene x sample expression matrix
aligned to the GlioStrata clinical cohort.

This unblocks the L2 layer (Verhaak subtype classification, Neftel state
deconvolution, TME deconvolution), which needs an expression matrix that the
clinical-only ingest (scripts/download_tcga_clinical.py) does not provide.

Outputs (under --out, default data/raw/tcga_gbm/expression):
    file_manifest.parquet     one row per GDC file, with parsed barcodes
    gene_annotation.parquet   gene_id, gene_name, gene_type (GENCODE v36)
    counts.parquet            raw counts  (`unstranded`),      gene x sample
    tpm.parquet               TPM         (`tpm_unstranded`),  gene x sample
    sample_metadata.parquet   one row per retained sample column
    _cache/<file_id>.tsv      download cache, separate from the matrices

Usage:
    # metadata only -- report what the cohort actually contains
    python scripts/download_tcga_expression.py --dry-run

    # small test subset
    python scripts/download_tcga_expression.py --limit 5

    # full cohort
    python scripts/download_tcga_expression.py
"""

import argparse
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

#: STAR counts files carry 4 leading summary rows that are not genes.
STAR_SUMMARY_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

#: Columns present in a GDC "augmented_star_gene_counts" TSV.
STAR_COLUMNS = [
    "gene_id",
    "gene_name",
    "gene_type",
    "unstranded",
    "stranded_first",
    "stranded_second",
    "tpm_unstranded",
    "fpkm_unstranded",
    "fpkm_uq_unstranded",
]

GENE_ANNOTATION_COLUMNS = ["gene_id", "gene_name", "gene_type"]

#: Value column written to counts.parquet / tpm.parquet.
COUNTS_COLUMN = "unstranded"
TPM_COLUMN = "tpm_unstranded"

#: TCGA sample-type codes (barcode positions 14-15).
SAMPLE_TYPE_CODES = {
    "01": "Primary Tumor",
    "02": "Recurrent Tumor",
    "05": "Additional New Primary",
    "06": "Metastatic",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
}

#: Sample-type codes kept unless --all-sample-types is given.
DEFAULT_SAMPLE_TYPE_CODES = ("01",)

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

_BARCODE_RE = re.compile(
    r"^(?P<patient_id>TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})"
    r"(?:-(?P<sample>(?P<sample_type_code>\d{2})(?P<vial>[A-Z]?)))?"
    r"(?P<rest>-.*)?$"
)


# --------------------------------------------------------------------------
# GDC metadata
# --------------------------------------------------------------------------

def build_file_filters(project_id: str, workflow_type: str) -> dict:
    """GDC `files` filter for a project's STAR gene-expression quantifications."""
    return {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [project_id],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": ["Transcriptome Profiling"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": ["Gene Expression Quantification"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "analysis.workflow_type",
                    "value": [workflow_type],
                },
            },
        ],
    }


def fetch_file_metadata(
    project_id: str = "TCGA-GBM",
    workflow_type: str = "STAR - Counts",
    page_size: int = 500,
) -> list[dict]:
    """Fetch all matching file records from the GDC `files` endpoint (paged)."""
    filters = build_file_filters(project_id, workflow_type)
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
        warnings = response.json().get("warnings", {})
        if warnings:
            print(f"WARNING - GDC API field warnings: {warnings}")
        total = payload["pagination"]["total"]
        hits.extend(payload["hits"])
        offset += page_size
        if offset >= total:
            break
    if len(hits) != total:
        raise RuntimeError(f"paging mismatch: got {len(hits)} of {total} files")
    return hits


def flatten_file_record(hit: dict) -> dict:
    """Flatten one GDC file record to a single-row dict.

    The aliquot barcode is the file-level TCGA barcode
    (e.g. TCGA-12-0819-01A-01R-1849-01); the sample barcode is its
    patient+sample prefix (TCGA-12-0819-01A).
    """
    cases = hit.get("cases") or []
    case = cases[0] if cases else {}
    samples = case.get("samples") or []
    sample = samples[0] if samples else {}

    aliquot_barcode = None
    for portion in sample.get("portions") or []:
        for analyte in portion.get("analytes") or []:
            for aliquot in analyte.get("aliquots") or []:
                if aliquot.get("submitter_id"):
                    aliquot_barcode = aliquot["submitter_id"]
                    break
            if aliquot_barcode:
                break
        if aliquot_barcode:
            break

    sample_barcode = sample.get("submitter_id")
    barcode = aliquot_barcode or sample_barcode or case.get("submitter_id")
    parsed = parse_barcode(barcode)

    return {
        "file_id": hit.get("file_id"),
        "file_name": hit.get("file_name"),
        "file_size": hit.get("file_size"),
        "md5sum": hit.get("md5sum"),
        "workflow_type": (hit.get("analysis") or {}).get("workflow_type"),
        "case_id": case.get("case_id"),
        "case_submitter_id": case.get("submitter_id"),
        "sample_barcode": sample_barcode,
        "aliquot_barcode": aliquot_barcode,
        "gdc_sample_type": sample.get("sample_type"),
        "patient_id": parsed["patient_id"],
        "sample_type_code": parsed["sample_type_code"],
        "vial": parsed["vial"],
    }


# --------------------------------------------------------------------------
# Barcode parsing
# --------------------------------------------------------------------------

def parse_barcode(barcode) -> dict:
    """Parse a TCGA barcode into patient_id / sample_type_code / vial.

    The patient_id form matches clinical_cohort.parquet's `patient_id`
    (submitter_id), e.g. "TCGA-12-0819". Sample barcodes carry extra
    suffixes -- TCGA-12-0819-01A-01R-1849-01 -> patient "TCGA-12-0819",
    sample_type_code "01" (primary tumor), vial "A".

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


def filter_sample_types(manifest: pd.DataFrame, keep_codes) -> tuple[pd.DataFrame, dict]:
    """Keep only rows whose sample_type_code is in `keep_codes`.

    Returns (filtered manifest, report dict with the full breakdown).
    """
    keep = list(keep_codes)
    breakdown = {}
    for code, count in manifest["sample_type_code"].value_counts(
        dropna=False
    ).items():
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


# --------------------------------------------------------------------------
# Duplicate aliquot resolution
# --------------------------------------------------------------------------

def resolve_duplicate_aliquots(manifest: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Reduce the manifest to one file per patient, deterministically.

    Rule: within a patient, sort by (sample_type_code, aliquot_barcode,
    file_id) ascending and keep the first row. Because TCGA aliquot barcodes
    encode vial/portion/plate in ascending character order, this picks the
    lowest-numbered vial and portion (e.g. -01A-01R- over -01B-01R-), and
    file_id breaks any remaining tie, so reruns are byte-identical.

    Returns (deduplicated manifest sorted by patient_id, report dict).
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
    """Cache location for one downloaded STAR counts file."""
    return Path(cache_dir) / f"{file_id}.tsv"


def pending_file_ids(file_ids, cache_dir: Path) -> list[str]:
    """File ids not already cached on disk (resume support)."""
    pending = []
    for file_id in file_ids:
        path = cache_path_for(cache_dir, file_id)
        if path.is_file() and path.stat().st_size > 0:
            continue
        pending.append(file_id)
    return pending


def download_batch(file_ids: list[str], cache_dir: Path, timeout: int = 600) -> int:
    """Download a batch of files via POST /data and write them to the cache.

    The GDC /data endpoint returns a single .tar.gz for a multi-id request
    (and a bare TSV for a single id), so batching is far cheaper than one
    request per file. Returns the number of files written.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    response = requests.post(
        GDC_DATA_URL,
        json={"ids": list(file_ids)},
        headers={"Content-Type": "application/json"},
        timeout=timeout,
    )
    if response.status_code != 200:
        raise RuntimeError(
            f"GDC data error {response.status_code}: {response.text[:500]}"
        )
    written = 0
    if len(file_ids) == 1:
        cache_path_for(cache_dir, file_ids[0]).write_bytes(response.content)
        return 1
    with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            parts = Path(member.name).parts
            if len(parts) < 2:
                continue  # MANIFEST.txt at the archive root
            file_id = parts[0]
            if file_id not in set(file_ids):
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            cache_path_for(cache_dir, file_id).write_bytes(handle.read())
            written += 1
    return written


def download_files(
    file_ids: list[str],
    cache_dir: Path,
    batch_size: int = 20,
) -> dict:
    """Download every missing file, batching requests and showing progress."""
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
    for start in range(0, len(todo), batch_size):
        batch = todo[start : start + batch_size]
        written += download_batch(batch, cache_dir)
        done = min(start + batch_size, len(todo))
        elapsed = time.time() - started
        print(
            f"  [{done}/{len(todo)}] fetched -> {elapsed:6.1f}s elapsed",
            flush=True,
        )
    return {
        "requested": len(file_ids),
        "skipped_cached": skipped,
        "downloaded": written,
        "elapsed_seconds": round(time.time() - started, 1),
    }


# --------------------------------------------------------------------------
# STAR counts parsing / matrix assembly
# --------------------------------------------------------------------------

def read_star_counts(path) -> pd.DataFrame:
    """Read one GDC STAR counts TSV, dropping the 4 leading summary rows.

    The file begins with a `# gene-model:` comment line, then the header, then
    N_unmapped / N_multimapping / N_noFeature / N_ambiguous, then one row per
    gene. Those 4 rows are alignment summaries, not genes, and MUST be dropped
    before any expression analysis.
    """
    frame = pd.read_csv(path, sep="\t", comment="#", dtype={"gene_id": "string"})
    missing = [c for c in STAR_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"{path}: missing STAR columns {missing}")
    frame = frame[~frame["gene_id"].isin(STAR_SUMMARY_ROWS)].copy()
    frame = frame.reset_index(drop=True)
    if frame["gene_id"].duplicated().any():
        raise ValueError(f"{path}: duplicate gene_id values")
    return frame


def gene_annotation(frame: pd.DataFrame) -> pd.DataFrame:
    """Gene annotation table (gene_id, gene_name, gene_type), gene_id-sorted."""
    annot = frame[GENE_ANNOTATION_COLUMNS].copy()
    annot["gene_id"] = annot["gene_id"].astype("string")
    annot = annot.sort_values("gene_id", kind="mergesort").reset_index(drop=True)
    return annot


def assemble_matrices(per_sample: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build gene x sample counts and TPM matrices plus a gene annotation.

    Args:
        per_sample: {sample column name -> STAR frame from read_star_counts}.

    Returns:
        (counts, tpm, annotation). Both matrices are indexed by gene_id in
        ascending order with identical, sorted column order, so the two stay
        aligned and reruns are byte-identical.
    """
    if not per_sample:
        raise ValueError("no samples to assemble")

    sample_names = sorted(per_sample)
    reference = None
    counts_cols = {}
    tpm_cols = {}
    for name in sample_names:
        frame = per_sample[name].set_index("gene_id").sort_index(kind="mergesort")
        if reference is None:
            reference = frame
        elif not frame.index.equals(reference.index):
            raise ValueError(f"gene_id index of sample {name} differs from reference")
        counts_cols[name] = pd.to_numeric(frame[COUNTS_COLUMN], errors="coerce")
        tpm_cols[name] = pd.to_numeric(frame[TPM_COLUMN], errors="coerce")

    counts = pd.DataFrame(counts_cols, index=reference.index, columns=sample_names)
    tpm = pd.DataFrame(tpm_cols, index=reference.index, columns=sample_names)
    counts.index.name = "gene_id"
    tpm.index.name = "gene_id"

    annot = gene_annotation(reference.reset_index())
    return counts, tpm, annot


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def print_metadata_report(manifest: pd.DataFrame, total_files: int) -> None:
    """ASCII-only report of what the GDC query matched."""
    size_gb = float(pd.to_numeric(manifest["file_size"], errors="coerce").sum()) / 1e9
    print("\nGDC file query report:")
    print(f"  Files matched: {total_files}")
    print(f"  Distinct cases covered: {manifest['case_submitter_id'].nunique()}")
    print(f"  Distinct patients (parsed barcode): {manifest['patient_id'].nunique()}")
    print(f"  Total uncompressed size: {size_gb:.3f} GB")


def print_sample_type_report(report: dict) -> None:
    """ASCII-only sample-type breakdown."""
    print("\nSample-type breakdown:")
    for code, count in report["breakdown"].items():
        print(f"  {code} {sample_type_label(code):24s} -> {count} files")
    kept = ", ".join(f"{c} ({sample_type_label(c)})" for c in report["kept_codes"])
    print(f"  Keeping: {kept}")
    print(
        f"  Files kept: {report['kept_files']} "
        f"-> dropped: {report['dropped_files']}"
    )


def print_duplicate_report(report: dict) -> None:
    """ASCII-only duplicate-aliquot report."""
    print("\nDuplicate aliquot resolution:")
    print(f"  Input files: {report['input_files']}")
    print(f"  Patients with more than one aliquot: {report['patients_with_duplicates']}")
    print(f"  Files dropped as duplicates: {report['files_dropped_as_duplicate']}")
    print(f"  Patients retained: {report['patients']}")
    print(
        "  Rule: lowest (sample_type_code, aliquot_barcode, file_id) wins "
        "-> deterministic across reruns"
    )


def clinical_overlap(patient_ids, clinical_path) -> dict:
    """Overlap between the expression patients and the clinical cohort."""
    path = Path(clinical_path)
    if not path.is_file():
        return {"clinical_rows": None, "overlap": None, "expression_only": None}
    clinical = pd.read_parquet(path)
    clinical_ids = set(clinical["patient_id"].dropna().astype(str))
    expr_ids = set(str(p) for p in patient_ids)
    return {
        "clinical_rows": int(len(clinical)),
        "overlap": len(expr_ids & clinical_ids),
        "expression_only": len(expr_ids - clinical_ids),
        "clinical_only": len(clinical_ids - expr_ids),
    }


def print_overlap_report(report: dict) -> None:
    """ASCII-only clinical-overlap report."""
    print("\nClinical cohort overlap:")
    if report["clinical_rows"] is None:
        print("  clinical_cohort.parquet not found -> overlap not computed")
        return
    print(f"  Clinical cohort rows: {report['clinical_rows']}")
    print(f"  Patients with BOTH clinical and expression: {report['overlap']}")
    print(f"  Expression-only patients: {report['expression_only']}")
    print(f"  Clinical-only patients: {report['clinical_only']}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default="data/raw/tcga_gbm/expression",
        help="Output directory for the assembled matrices",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Download cache (default: <out>/_cache)",
    )
    parser.add_argument(
        "--clinical",
        default="data/raw/tcga_gbm/clinical_cohort.parquet",
        help="Clinical cohort parquet, for the overlap report",
    )
    parser.add_argument("--project", default="TCGA-GBM", help="GDC project ID")
    parser.add_argument(
        "--workflow", default="STAR - Counts", help="GDC analysis.workflow_type"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only download/assemble the first N patients (testing)",
    )
    parser.add_argument(
        "--all-sample-types",
        action="store_true",
        help="Retain non-primary-tumor samples (recurrent, normal, ...)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=20, help="Files per POST /data request"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report the GDC metadata only; download nothing",
    )
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir) if args.cache_dir else out_dir / "_cache"

    print(f"Querying GDC files endpoint for {args.project} / {args.workflow}...")
    hits = fetch_file_metadata(args.project, args.workflow)
    manifest = pd.DataFrame([flatten_file_record(h) for h in hits])
    manifest = manifest.sort_values("file_id", kind="mergesort").reset_index(drop=True)
    print_metadata_report(manifest, total_files=len(hits))

    keep_codes = (
        sorted(c for c in manifest["sample_type_code"].dropna().unique())
        if args.all_sample_types
        else list(DEFAULT_SAMPLE_TYPE_CODES)
    )
    manifest, st_report = filter_sample_types(manifest, keep_codes)
    print_sample_type_report(st_report)

    manifest, dup_report = resolve_duplicate_aliquots(manifest)
    print_duplicate_report(dup_report)

    if args.limit is not None:
        manifest = manifest.head(args.limit).reset_index(drop=True)
        print(f"\n--limit {args.limit} -> restricting to {len(manifest)} patients")

    manifest_path = out_dir / "file_manifest.parquet"
    manifest.to_parquet(manifest_path, index=False)
    print(f"\nSaved file manifest: {manifest_path}  ({len(manifest)} rows)")

    projected_gb = (
        float(pd.to_numeric(manifest["file_size"], errors="coerce").sum()) / 1e9
    )
    print(f"Projected uncompressed download: {projected_gb:.3f} GB")

    if args.dry_run:
        print("\n--dry-run -> no files downloaded, no matrices assembled")
        print_overlap_report(clinical_overlap(manifest["patient_id"], args.clinical))
        return

    dl_report = download_files(
        manifest["file_id"].tolist(), cache_dir, batch_size=args.batch_size
    )
    print(
        f"Downloaded {dl_report['downloaded']} files "
        f"({dl_report['skipped_cached']} cached) in "
        f"{dl_report['elapsed_seconds']}s"
    )

    print("\nParsing STAR counts files...")
    per_sample = {}
    genes_before = None
    for row in manifest.itertuples(index=False):
        frame = read_star_counts(cache_path_for(cache_dir, row.file_id))
        if genes_before is None:
            genes_before = len(frame) + len(STAR_SUMMARY_ROWS)
        per_sample[row.patient_id] = frame
    print(f"  Parsed {len(per_sample)} samples")
    print(
        f"  Rows per file before filtering: {genes_before} -> "
        f"after dropping {len(STAR_SUMMARY_ROWS)} summary rows: "
        f"{genes_before - len(STAR_SUMMARY_ROWS)}"
    )

    counts, tpm, annot = assemble_matrices(per_sample)

    counts_path = out_dir / "counts.parquet"
    tpm_path = out_dir / "tpm.parquet"
    annot_path = out_dir / "gene_annotation.parquet"
    counts.to_parquet(counts_path)
    tpm.to_parquet(tpm_path)
    annot.to_parquet(annot_path, index=False)

    sample_meta = manifest[
        [
            "patient_id",
            "case_submitter_id",
            "sample_barcode",
            "aliquot_barcode",
            "sample_type_code",
            "gdc_sample_type",
            "file_id",
        ]
    ].copy()
    sample_meta_path = out_dir / "sample_metadata.parquet"
    sample_meta.to_parquet(sample_meta_path, index=False)

    print("\nSaved matrices:")
    print(f"  {counts_path}  -> {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"  {tpm_path}  -> {tpm.shape[0]} genes x {tpm.shape[1]} samples")
    print(f"  {annot_path}  -> {annot.shape[0]} genes x {annot.shape[1]} cols")
    print(f"  {sample_meta_path}  -> {sample_meta.shape[0]} rows")
    print(f"  Cache directory: {cache_dir}")

    print("\nGene counts:")
    print(f"  Genes before filtering (incl. summary rows): {genes_before}")
    print(f"  Genes after dropping summary rows: {counts.shape[0]}")
    print("  No expression-level gene filtering applied (kept all GENCODE genes)")

    print_overlap_report(clinical_overlap(counts.columns, args.clinical))


if __name__ == "__main__":
    main()
