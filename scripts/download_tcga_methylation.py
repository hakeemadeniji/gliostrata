"""
Download TCGA-GBM DNA methylation beta values (Illumina 27k and 450k) from the
GDC REST API and call MGMT promoter methylation status per patient.

WHY THIS EXISTS
    TCGA-GBM is the development cohort and CGGA the external validation
    cohort, and they barely share features. CGGA carries MGMT methylation for
    312/333 patients; the TCGA clinical ingest carries none, and MGMT is the
    single most consequential molecular covariate in glioblastoma. This script
    closes that gap on the TCGA side (the IDH half is
    scripts/download_tcga_mutations.py).

THE MODEL USED HERE IS *NOT* MGMT-STP27 -- read this before citing it
    The protocol asks for the published MGMT-STP27 logistic model. Its exact
    coefficients and intercept could NOT be verified against the primary
    publication during development, and a fabricated coefficient is far worse
    than an honest substitute. So this script implements a deliberately
    simple, fully transparent alternative and labels it as such:

        score       = mean beta of cg12434587 and cg12981137
                      (two MGMT promoter probes; both verified present in
                      real files from BOTH platforms)
        methylated  = score >= threshold        (default 0.35)

    The default threshold is an A-PRIORI working value fixed before the cohort
    was downloaded. It is not a published cut-point and it has not been
    calibrated here. The continuous ``mgmt_score`` is the primary output; the
    binary call is a convenience that anyone can re-derive at another
    threshold. The report prints the methylated fraction across a grid of
    thresholds so the sensitivity of the call is visible.

    TO BE VERIFIED BY AUTHOR: the MGMT-STP27 coefficients, intercept and
    decision threshold from the primary publication's supplementary material.
    Once obtained, implement it in ``mgmt_score()`` behind a new method name,
    bump CALL_METHOD, and re-run -- nothing else needs to change.

THREE STATES, NOT TWO
      1      methylated     both probes measured, score >= threshold
      0      unmethylated   both probes measured, score <  threshold
      <NA>   not callable   no beta file for the patient, or a required probe
                             is absent / NA in the file
    The third state is NOT "unmethylated". It is written as a nullable value,
    with the reason, and counted separately in every report.

    A patient with a file on both platforms gets the 450k call (STP27's
    native platform and the denser array); the other file is dropped.

WHAT IS NOT CONTROLLED HERE
    Beta values are used as distributed by GDC (SeSaMe normalisation for 450k,
    legacy level-3 for 27k). No cross-platform harmonisation is applied, so a
    call made on 27k and a call made on 450k are not guaranteed equivalent.
    The platform of every call is written to the output for exactly that
    reason; stratify on it, or include it as a covariate, when using the calls.

OUTPUTS (under --out, default data/raw/tcga_gbm/methylation)
    file_manifest.parquet     one row per GDC file, with parsed barcodes
    mgmt_calls.parquet        per-patient calls, INCLUDING not-callable patients
    mgmt_call_report.json     every count printed to stdout, machine-readable
    _cache/<file_id>.txt      raw beta-value files (the download cache)

    clinical_cohort.parquet is NEVER modified.

Usage:
    python scripts/download_tcga_methylation.py --dry-run
    python scripts/download_tcga_methylation.py --limit 5
    python scripts/download_tcga_methylation.py
    python scripts/download_tcga_methylation.py --threshold 0.30
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

#: Project this script ingests. TCGA-GBM only; not a general-purpose tool.
PROJECT_ID = "TCGA-GBM"

#: GDC platform names. TCGA-GBM was profiled on both.
PLATFORM_27K = "Illumina Human Methylation 27"
PLATFORM_450K = "Illumina Human Methylation 450"
PLATFORMS = (PLATFORM_27K, PLATFORM_450K)

#: Short platform labels written to the output.
PLATFORM_LABELS = {PLATFORM_27K: "27k", PLATFORM_450K: "450k"}

#: Which platform wins when a patient has a file on both (lower = preferred).
PLATFORM_PRIORITY = {PLATFORM_450K: 0, PLATFORM_27K: 1}

#: The two MGMT promoter probes used by the simple score. Both verified
#: present in real files from BOTH platforms (see the run report).
MGMT_PROBES = ("cg12434587", "cg12981137")

#: Valid range of an Illumina beta value (methylated / (methylated + unmethylated)
#: signal). Anything outside it is a parsing or platform failure, not a beta.
BETA_MIN = 0.0
BETA_MAX = 1.0

#: A-priori working threshold on the two-probe mean beta. NOT a published
#: cut-point; see the module docstring.
DEFAULT_THRESHOLD = 0.35

#: Name recorded in the output so nobody mistakes this for MGMT-STP27.
CALL_METHOD = "two_probe_mean_beta"
IS_MGMT_STP27 = False

#: Thresholds at which the methylated fraction is also reported.
SENSITIVITY_THRESHOLDS = (0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60)

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
    "platform",
    "cases.submitter_id",
    "cases.case_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.samples.portions.analytes.aliquots.submitter_id",
]

#: Output columns of mgmt_calls.parquet, in order.
CALL_COLUMNS = [
    "patient_id",
    "mgmt_methylated",
    "mgmt_score",
    "cg12434587",
    "cg12981137",
    "call_basis",
    "call_method",
    "platform",
    "threshold",
    "sample_barcode",
    "aliquot_barcode",
    "file_id",
]

#: Value of `call_basis` for each of the three states.
BASIS_CALLED = "two_probes_measured"
BASIS_PROBE_UNAVAILABLE = "required_probe_absent_or_NA"
BASIS_MISSING = "no_methylation_data"

DEFAULT_OUT = Path("data/raw/tcga_gbm/methylation")
DEFAULT_COHORT = Path("data/raw/tcga_gbm/clinical_cohort.parquet")

_BARCODE_RE = re.compile(
    r"^(?P<patient_id>TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})"
    r"(?:-(?P<sample>(?P<sample_type_code>\d{2})(?P<vial>[A-Z]?)))?"
    r"(?P<rest>-.*)?$"
)


# --------------------------------------------------------------------------
# GDC metadata
# --------------------------------------------------------------------------

def build_file_filters(project_id: str = PROJECT_ID, platform=None) -> dict:
    """GDC `files` filter for open methylation beta-value files."""
    content = [
        {"op": "in", "content": {"field": "cases.project.project_id",
                                 "value": [project_id]}},
        {"op": "in", "content": {"field": "data_category", "value": ["DNA Methylation"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ]
    if platform is not None:
        content.append({"op": "in", "content": {"field": "platform", "value": [platform]}})
    return {"op": "and", "content": content}


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

    ``samples[0]`` is not a safe way to learn what tissue a file describes:
    the order is arbitrary, and a file linked to a tumour/normal pair may list
    the normal first (true for 239 of 464 TCGA-GBM mutation files).
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
    """Pick the aliquot that describes the tumour, deterministically."""
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
    """Flatten one GDC file record to a single-row dict with a parsed barcode."""
    cases = hit.get("cases") or []
    case = cases[0] if cases else {}
    chosen = choose_tumour_aliquot(linked_aliquots(case))
    parsed = parse_barcode(chosen["aliquot_barcode"] or case.get("submitter_id"))

    return {
        "file_id": hit.get("file_id"),
        "file_name": hit.get("file_name"),
        "file_size": hit.get("file_size"),
        "md5sum": hit.get("md5sum"),
        "platform": hit.get("platform"),
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


def resolve_duplicate_files(manifest: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Reduce the manifest to ONE file per patient, deterministically.

    Order: sample_type_code, then PLATFORM_PRIORITY (450k before 27k), then
    aliquot_barcode, then file_id. A patient profiled on both arrays therefore
    keeps the 450k file; everything else is stable across reruns.
    """
    if manifest.empty:
        return manifest.copy(), {
            "input_files": 0, "patients": 0, "patients_with_duplicates": 0,
            "patients_on_both_platforms": 0, "files_dropped_as_duplicate": 0,
        }
    ordered = manifest.assign(
        _priority=manifest["platform"].map(PLATFORM_PRIORITY).fillna(99)
    )
    ordered = ordered.sort_values(
        ["patient_id", "sample_type_code", "_priority", "aliquot_barcode", "file_id"],
        kind="mergesort",
        na_position="last",
    )
    grouped = ordered.groupby("patient_id", dropna=False)
    deduped = ordered.drop_duplicates(subset=["patient_id"], keep="first").drop(
        columns="_priority"
    )
    deduped = deduped.sort_values("patient_id", kind="mergesort").reset_index(drop=True)
    report = {
        "input_files": int(len(manifest)),
        "patients": int(len(deduped)),
        "patients_with_duplicates": int((grouped.size() > 1).sum()),
        "patients_on_both_platforms": int((grouped["platform"].nunique() > 1).sum()),
        "files_dropped_as_duplicate": int(len(manifest) - len(deduped)),
    }
    return deduped, report

# --------------------------------------------------------------------------
# Download
# --------------------------------------------------------------------------

def cache_path_for(cache_dir: Path, file_id: str) -> Path:
    """Cache location for one downloaded beta-value file."""
    return Path(cache_dir) / f"{file_id}.txt"


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
    timeout: int = 1800,
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


def download_files(file_ids: list[str], cache_dir: Path, batch_size: int = 4) -> dict:
    """Download every missing file, batching requests and showing progress.

    The default batch is small because a 450k beta file is ~13 MB. A batch
    that keeps failing is retried ONE FILE AT A TIME, so a single bad file id
    costs that file and not the whole run; ids that still fail are returned
    under ``failed`` and reported, never silently dropped.
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
              "reported as not callable: " + ", ".join(failed[:10])
              + (" ..." if len(failed) > 10 else ""))
    return {
        "requested": len(file_ids),
        "skipped_cached": skipped,
        "downloaded": written,
        "failed": failed,
        "elapsed_seconds": round(time.time() - started, 1),
    }


# --------------------------------------------------------------------------
# Beta-value reading and the MGMT score
# --------------------------------------------------------------------------

def parse_beta(value) -> float | None:
    """Parse one beta value; anything that is not a number in [0, 1] is None.

    GDC writes ``NA`` for an unmeasurable probe. NaN also fails the range
    test, so both land on None rather than being averaged in as 0.
    """
    try:
        beta = float(value)
    except (TypeError, ValueError):
        return None
    if BETA_MIN <= beta <= BETA_MAX:
        return beta
    return None


def read_probe_values(path, probes=MGMT_PROBES) -> tuple[dict, int]:
    """Stream one beta-value file and return ``(values, n_lines)``.

    The files are headerless ``probe<TAB>beta`` text (27,578 lines on 27k,
    ~485,000 on 450k), so the file is streamed once and only the requested
    probes are kept. ``values`` has a key ONLY for probes present in the file;
    a present-but-``NA`` probe maps to ``None``. That lets the caller tell
    "probe not on this platform" apart from "probe measured as NA".
    """
    wanted = set(probes)
    values: dict = {}
    n_lines = 0
    with open(path, "r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            n_lines += 1
            probe, sep, value = line.partition("\t")
            if sep and probe in wanted and probe not in values:
                values[probe] = parse_beta(value.strip())
    return values, n_lines


def mgmt_score(values: dict, probes=MGMT_PROBES) -> float | None:
    """Mean beta over the required probes, or None unless ALL are measured.

    This is the transparent stand-in for MGMT-STP27, NOT that model; see the
    module docstring. A single missing probe makes the sample not callable
    rather than silently averaging over fewer probes.
    """
    betas = [values.get(p) for p in probes]
    if any(b is None for b in betas):
        return None
    return float(sum(betas) / len(betas))

# --------------------------------------------------------------------------
# The three-state call table
# --------------------------------------------------------------------------

def _call_row(pid: str, record: dict, threshold: float) -> dict:
    values = record.get("values") or {}
    score = mgmt_score(values)
    platform = record.get("platform")
    row = {
        "patient_id": pid,
        "mgmt_methylated": pd.NA,
        "mgmt_score": score,
        "cg12434587": values.get("cg12434587"),
        "cg12981137": values.get("cg12981137"),
        "call_basis": BASIS_PROBE_UNAVAILABLE,
        "call_method": CALL_METHOD,
        "platform": PLATFORM_LABELS.get(platform, platform),
        "threshold": float(threshold),
        "sample_barcode": record.get("sample_barcode"),
        "aliquot_barcode": record.get("aliquot_barcode"),
        "file_id": record.get("file_id"),
    }
    if score is not None:
        row["mgmt_methylated"] = int(score >= threshold)
        row["call_basis"] = BASIS_CALLED
    return row


def build_calls(
    sample_records: list[dict],
    cohort_patient_ids=None,
    threshold: float = DEFAULT_THRESHOLD,
) -> pd.DataFrame:
    """Assemble per-patient calls, INCLUDING patients with no methylation data.

    ``sample_records`` holds one dict per beta-value file with ``patient_id``,
    ``platform``, ``values`` (probe -> beta) and sample provenance. If a
    patient has more than one record the preferred platform (450k) wins. A
    patient absent from every record is NOT CALLABLE (nullable), never
    "unmethylated".
    """
    if not BETA_MIN <= float(threshold) <= BETA_MAX:
        raise ValueError(f"threshold must lie in [0, 1]; got {threshold!r}")
    ordered = sorted(
        sample_records,
        key=lambda r: (
            str(r.get("patient_id")),
            PLATFORM_PRIORITY.get(r.get("platform"), 99),
            str(r.get("aliquot_barcode") or ""),
            str(r.get("file_id") or ""),
        ),
    )
    by_patient: dict[str, dict] = {}
    for record in ordered:
        pid = str(record["patient_id"])
        if pid not in by_patient:
            by_patient[pid] = _call_row(pid, record, float(threshold))

    if cohort_patient_ids is not None:
        for pid in cohort_patient_ids:
            pid = str(pid)
            if pid in by_patient:
                continue
            by_patient[pid] = {
                "patient_id": pid,
                "mgmt_methylated": pd.NA,
                "mgmt_score": None,
                "cg12434587": None,
                "cg12981137": None,
                "call_basis": BASIS_MISSING,
                "call_method": None,
                "platform": None,
                "threshold": None,
                "sample_barcode": None,
                "aliquot_barcode": None,
                "file_id": None,
            }

    if not by_patient:
        return pd.DataFrame({c: pd.Series(dtype="object") for c in CALL_COLUMNS})
    out = pd.DataFrame([by_patient[p] for p in sorted(by_patient)])[CALL_COLUMNS]
    out["mgmt_methylated"] = out["mgmt_methylated"].astype("Int64")
    for col in ("mgmt_score", "cg12434587", "cg12981137", "threshold"):
        out[col] = pd.to_numeric(out[col], errors="coerce").astype("Float64")
    return out.reset_index(drop=True)

# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def summarize_calls(calls: pd.DataFrame) -> dict:
    """The counts behind every number printed to stdout."""
    n = int(len(calls))
    methylated = int((calls["mgmt_methylated"] == 1).sum())
    unmethylated = int((calls["mgmt_methylated"] == 0).sum())
    not_callable = int(calls["mgmt_methylated"].isna().sum())
    callable_n = methylated + unmethylated
    scores = calls["mgmt_score"].dropna().astype(float)

    sensitivity = []
    for thr in SENSITIVITY_THRESHOLDS:
        k = int((scores >= thr).sum())
        sensitivity.append({
            "threshold": thr,
            "methylated": k,
            "fraction_of_callable": (k / len(scores)) if len(scores) else None,
        })

    # Where does the score actually sit? A bimodal histogram is the evidence
    # that a threshold is meaningful at all; a flat one says it is not.
    edges = [round(i / 10, 1) for i in range(11)]
    histogram = []
    if len(scores):
        counts = pd.cut(scores, bins=edges, include_lowest=True, right=False)
        for interval, count in counts.value_counts(sort=False).items():
            histogram.append({"bin": str(interval), "n": int(count)})

    by_platform = {}
    if "platform" in calls.columns:
        called = calls[calls["mgmt_methylated"].notna()]
        for label, group in called.groupby(called["platform"].astype(str)):
            by_platform[str(label)] = {
                "callable": int(len(group)),
                "methylated": int((group["mgmt_methylated"] == 1).sum()),
            }
    return {
        "patients_in_table": n,
        "callable": callable_n,
        "methylated": methylated,
        "unmethylated": unmethylated,
        "not_callable": not_callable,
        "methylated_fraction_of_callable": (methylated / callable_n) if callable_n else None,
        "by_call_basis": {str(k): int(v) for k, v in calls["call_basis"].value_counts().items()},
        "by_platform": by_platform,
        "score_histogram": histogram,
        "threshold_sensitivity": sensitivity,
        "is_mgmt_stp27": IS_MGMT_STP27,
        "call_method": CALL_METHOD,
    }


def probe_coverage(sample_records: list[dict]) -> dict:
    """Per platform: how many files carry each required probe, and how many NA."""
    out: dict = {}
    for record in sample_records:
        label = PLATFORM_LABELS.get(record.get("platform"), record.get("platform"))
        entry = out.setdefault(label, {
            "files": 0,
            "probes": {p: {"measured": 0, "na": 0, "absent": 0} for p in MGMT_PROBES},
        })
        entry["files"] += 1
        values = record.get("values") or {}
        for probe in MGMT_PROBES:
            if probe not in values:
                entry["probes"][probe]["absent"] += 1
            elif values[probe] is None:
                entry["probes"][probe]["na"] += 1
            else:
                entry["probes"][probe]["measured"] += 1
    return out


def coverage(calls: pd.DataFrame, patient_ids, label: str) -> dict:
    """How many of a given patient set gain an MGMT call."""
    ids = {str(p) for p in patient_ids}
    sub = calls[calls["patient_id"].astype(str).isin(ids)]
    has = sub["mgmt_methylated"].notna()
    return {
        "set": label,
        "patients": len(ids),
        "with_mgmt_call": int(has.sum()),
        "methylated": int((sub["mgmt_methylated"] == 1).sum()),
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


#: Sanity band for the methylated fraction, taken from the project brief
#: ("roughly 40-50% methylated in TCGA-GBM"). TO BE VERIFIED against a
#: primary source; outside this band the report says so rather than hiding it.
EXPECTED_METHYLATED_RANGE = (0.30, 0.60)


def print_report(summary: dict, coverage_rows: list[dict], probe_cov: dict) -> None:
    """ASCII-only summary (cp1252 consoles choke on non-ASCII)."""
    method = summary["call_method"]
    if not summary["is_mgmt_stp27"]:
        method += "  -- NOT MGMT-STP27 (see module docstring)"
    print(f"\nMGMT promoter methylation -- method: {method}")
    callable_n = summary["callable"]
    print(f"  methylated      {summary['methylated']:5d}   "
          f"{_pct(summary['methylated'], callable_n)} of callable")
    print(f"  unmethylated    {summary['unmethylated']:5d}   "
          f"{_pct(summary['unmethylated'], callable_n)} of callable")
    print(f"  not callable    {summary['not_callable']:5d}   "
          f"{_pct(summary['not_callable'], summary['patients_in_table'])} of table "
          f"(NOT unmethylated; written as <NA>)")
    for basis, count in sorted(summary["by_call_basis"].items()):
        print(f"    basis {basis:36} {count:5d}")

    print("\nBy platform (callable calls only):")
    for label in sorted(summary["by_platform"]):
        row = summary["by_platform"][label]
        print(f"  {label:8} callable {row['callable']:5d}   methylated "
              f"{row['methylated']:5d}  {_pct(row['methylated'], row['callable'])}")

    print("\nScore histogram (mean beta of cg12434587 and cg12981137):")
    peak = max([h["n"] for h in summary["score_histogram"]] or [1])
    for row in summary["score_histogram"]:
        bar = "#" * int(round(40 * row["n"] / max(peak, 1)))
        print(f"  {row['bin']:14} {row['n']:5d}  {bar}")

    print("\nThreshold sensitivity (methylated fraction of callable):")
    for row in summary["threshold_sensitivity"]:
        mark = "   <- default" if abs(row["threshold"] - DEFAULT_THRESHOLD) < 1e-9 else ""
        frac = row["fraction_of_callable"]
        shown = "n/a" if frac is None else f"{100.0 * frac:5.1f}%"
        print(f"  score >= {row['threshold']:.2f}   {row['methylated']:5d}   {shown}{mark}")

    frac = summary["methylated_fraction_of_callable"]
    if frac is not None:
        low, high = EXPECTED_METHYLATED_RANGE
        verdict = "within" if low <= frac <= high else "OUTSIDE"
        print(f"\nSanity check: methylated fraction {_pct(summary['methylated'], callable_n)} is "
              f"{verdict} the expected ~{100*low:.0f}-{100*high:.0f}% band for TCGA-GBM")
        if verdict == "OUTSIDE":
            print("  (expectation per project brief; verify against a primary source)")

    print("\nRequired-probe coverage by platform:")
    for label in sorted(probe_cov, key=str):
        entry = probe_cov[label]
        print(f"  {label}: {entry['files']} files")
        for probe, counts in entry["probes"].items():
            print(f"    {probe}  measured {counts['measured']:4d}   "
                  f"NA {counts['na']:4d}   absent {counts['absent']:4d}")

    print("\nCoverage -- how many patients gain a call:")
    print(f"  {'patient set':28}{'n':>6}{'with call':>11}{'methylated':>11}{'missing':>9}")
    for row in coverage_rows:
        print(f"  {row['set']:28}{row['patients']:6d}{row['with_mgmt_call']:11d}"
              f"{row['methylated']:11d}{row['still_missing']:9d}")

# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--project", default=PROJECT_ID)
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Output directory")
    parser.add_argument("--platform", choices=("all", "27k", "450k"), default="all",
                        help="restrict to one array platform (default: both)")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                        help="score >= threshold is methylated (default: %(default)s; "
                             "an a-priori working value, NOT a published cut-point)")
    parser.add_argument("--cohort", default=str(DEFAULT_COHORT),
                        help="clinical_cohort.parquet; its patients are listed even "
                             "when no beta file exists for them")
    parser.add_argument("--no-cohort", action="store_true",
                        help="do not read the clinical cohort")
    parser.add_argument("--expression", default="data/raw/tcga_gbm/expression/tpm.parquet",
                        help="gene x sample matrix; its columns are the "
                             "expression-covered patients (coverage report only)")
    parser.add_argument("--limit", type=int, default=None,
                        help="use only the first N patients (after deduplication)")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--all-sample-types", action="store_true",
                        help="keep recurrent/normal samples too (default: primary only)")
    parser.add_argument("--dry-run", action="store_true", help="report metadata only")
    args = parser.parse_args()

    print("TCGA-GBM methylation ingest -- MGMT promoter status")
    print(f"  Method: {CALL_METHOD}; MGMT-STP27 = {IS_MGMT_STP27}; "
          f"threshold = {args.threshold}")
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
    if args.platform != "all":
        wanted = PLATFORM_27K if args.platform == "27k" else PLATFORM_450K
        manifest = manifest[manifest["platform"] == wanted].reset_index(drop=True)

    total_mb = float(manifest["file_size"].sum()) / 1e6
    print(f"\nGDC metadata: {len(manifest)} open beta-value files, "
          f"{manifest['patient_id'].nunique()} patients, {total_mb:.0f} MB")
    for platform, count in manifest["platform"].value_counts().items():
        subset = manifest[manifest["platform"] == platform]
        print(f"  {PLATFORM_LABELS.get(platform, platform):5} {count:4d} files, "
              f"{subset['patient_id'].nunique():4d} patients, "
              f"{subset['file_size'].sum() / 1e6:6.0f} MB")
    for code, count in manifest["sample_type_code"].value_counts(dropna=False).items():
        print(f"  sample type {code}: {count:4d} files  ({sample_type_label(code)})")
    if args.dry_run:
        print("\n--dry-run -> nothing downloaded, nothing called")
        return

    if args.all_sample_types:
        kept, sample_report = manifest.copy(), {"kept_codes": "all"}
    else:
        kept, sample_report = select_sample_types(manifest)
    deduped, dedup_report = resolve_duplicate_files(kept)
    print(f"\nPrimary-tumour files kept: {len(kept)} (dropped {len(manifest) - len(kept)})")
    print(f"One file per patient: {dedup_report['patients']} patients; "
          f"{dedup_report['patients_on_both_platforms']} on both platforms "
          f"(450k kept); {dedup_report['files_dropped_as_duplicate']} files dropped")
    if args.limit is not None:
        deduped = deduped.head(int(args.limit)).copy()
        print(f"--limit {args.limit}: calling the first {len(deduped)} patient(s) only")

    # Batch each platform separately: a 27k file is ~0.8 MB, a 450k file ~13 MB.
    download_report = []
    for platform in PLATFORMS:
        ids = deduped.loc[deduped["platform"] == platform, "file_id"].tolist()
        if ids:
            print(f"\n{PLATFORM_LABELS[platform]}:")
            report = download_files(ids, cache_dir, batch_size=int(args.batch_size))
            report["platform"] = PLATFORM_LABELS[platform]
            download_report.append(report)

    records: list[dict] = []
    unavailable: list[str] = []
    for row in deduped.itertuples(index=False):
        path = cache_path_for(cache_dir, row.file_id)
        if not path.is_file() or path.stat().st_size == 0:
            # The download failed. "Not callable" is then exactly true for this
            # patient; count it below rather than crashing the run.
            unavailable.append(str(row.file_id))
            continue
        values, n_lines = read_probe_values(path)
        if not values:
            raise RuntimeError(
                f"{path}: neither MGMT probe found in {n_lines} lines -- "
                "this is not a probe-indexed beta-value file"
            )
        records.append({
            "patient_id": row.patient_id,
            "platform": row.platform,
            "values": values,
            "n_lines": n_lines,
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
              "only patients with a beta file are listed")

    calls = build_calls(records, cohort_ids, threshold=float(args.threshold))
    summary = summarize_calls(calls)
    probe_cov = probe_coverage(records)

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
    calls.to_parquet(out_dir / "mgmt_calls.parquet", index=False)
    report = {
        "project": args.project,
        "platform_filter": args.platform,
        "gdc_files": int(len(manifest)),
        "gdc_patients": int(manifest["patient_id"].nunique()),
        "total_mb": round(total_mb, 1),
        "call_method": CALL_METHOD,
        "is_mgmt_stp27": IS_MGMT_STP27,
        "probes": list(MGMT_PROBES),
        "threshold": float(args.threshold),
        "threshold_note": (
            "A-priori working value, not a published cut-point. MGMT-STP27 "
            "coefficients TO BE VERIFIED BY AUTHOR against the primary publication."
        ),
        "no_cross_platform_harmonisation": True,
        "sample_type_filter": sample_report,
        "duplicate_resolution": dedup_report,
        "download": download_report,
        "files_unavailable_after_download": unavailable,
        "probe_coverage": probe_cov,
        "summary": summary,
        "coverage": coverage_rows,
        "expression_non_sample_columns": expression_ignored,
        "limited_to_first_n_patients": args.limit,
    }
    (out_dir / "mgmt_call_report.json").write_text(
        json.dumps(_jsonable(report), indent=2), encoding="ascii"
    )

    print_report(summary, coverage_rows, probe_cov)
    print(f"\nWrote: {out_dir / 'file_manifest.parquet'}")
    print(f"      {out_dir / 'mgmt_calls.parquet'}")
    print(f"      {out_dir / 'mgmt_call_report.json'}")


if __name__ == "__main__":
    main()
