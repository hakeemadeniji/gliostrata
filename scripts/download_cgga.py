"""
Download the Chinese Glioma Genome Atlas (CGGA) and map it to the GlioStrata
cohort schema.

CGGA is this project's PRE-REGISTERED EXTERNAL VALIDATION COHORT for protocol
tests T1 and T2 (docs/preregistration.md Sections 2 and 3.3). This script
DOWNLOADS AND CHARACTERISES it only. It fits nothing, tunes nothing, selects
nothing, and computes no model metric and no feature-outcome association on
CGGA data. Reading the external partition for anything other than the single
registered confirmatory evaluation contaminates it irreversibly.

Datasets (all open, no registration or login required -- verified by
unauthenticated HTTP GET against www.cgga.org.cn):

    mRNAseq_693     693 RNA-seq samples, counts + FPKM, clinical
    mRNAseq_325     325 RNA-seq samples, counts + FPKM, clinical (disjoint
                    sample ids from mRNAseq_693)
    mRNA_array_301  301 Agilent microarray samples, clinical; the only batch
                    carrying a precomputed TCGA_subtypes label

Outputs (under --out, default data/raw/cgga):
    _cache/<file>.zip             download cache (resume support)
    _cache/extracted/<file>.txt   unpacked raw text, byte-identical to source
    clinical_raw.parquet          every clinical row of every batch, UNFILTERED
    clinical_cohort.parquet       mapped + filtered cohort (GlioStrata schema)
    counts.parquet                gene x sample raw counts, RNA-seq batches
    gene_annotation.parquet       gene_name (CGGA ships SYMBOLS, not Ensembl)
    harmonisation.json            batch-correction and gene-id facts

Usage:
    python scripts/download_cgga.py --dry-run
    python scripts/download_cgga.py
    python scripts/download_cgga.py --no-filter
"""

import argparse
import io
import json
import time
import zipfile
from pathlib import Path

import pandas as pd
import requests

CGGA_DOWNLOAD_URL = "https://www.cgga.org.cn/download"

#: One entry per downloadable CGGA artefact. `file` is the portal's
#: `download/<date>/<name>` path and `type` its download-counter key; both are
#: required by the portal endpoint. Order is fixed so reports are deterministic.
CGGA_FILES = {
    "mRNAseq_693_clinical": {
        "file": "download/20200506/CGGA.mRNAseq_693_clinical.20200506.txt.zip",
        "type": "mRNAseq_693_clinical",
        "time": "20200506",
        "dataset": "mRNAseq_693",
        "kind": "clinical",
    },
    "mRNAseq_693_counts": {
        "file": (
            "download/20220620/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip"
        ),
        "type": "mRNAseq_693_counts",
        "time": "20220620",
        "dataset": "mRNAseq_693",
        "kind": "counts",
    },
    "mRNAseq_325_clinical": {
        "file": "download/20200506/CGGA.mRNAseq_325_clinical.20200506.txt.zip",
        "type": "mRNAseq_325_clinical",
        "time": "20200506",
        "dataset": "mRNAseq_325",
        "kind": "clinical",
    },
    "mRNAseq_325_counts": {
        "file": (
            "download/20220620/CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip"
        ),
        "type": "mRNAseq_325_counts",
        "time": "20220620",
        "dataset": "mRNAseq_325",
        "kind": "counts",
    },
    "mRNA_array_301_clinical": {
        "file": "download/20200506/CGGA.mRNA_array_301_clinical.20200506.txt.zip",
        "type": "mRNA_array_301_clinical",
        "time": "20200506",
        "dataset": "mRNA_array_301",
        "kind": "clinical",
    },
}

#: RNA-seq batches whose counts matrices are assembled into counts.parquet.
#: The microarray batch is intentionally excluded: Agilent intensities are not
#: counts and ComBat-seq is a negative-binomial method for count data.
EXPRESSION_DATASETS = ("mRNAseq_693", "mRNAseq_325")

#: Canonical clinical column names. CGGA spells the 1p19q column with a
#: lower-case "c" in the RNA-seq batches and an upper-case "C" in the
#: microarray batch, and appends the coding legend to three headers.
CLINICAL_COLUMN_ALIASES = {
    "cgga_id": "cgga_id",
    "tcga_subtypes": "tcga_subtypes",
    "prs_type": "prs_type",
    "histology": "histology",
    "grade": "grade",
    "gender": "gender",
    "age": "age_years",
    "os": "os_days",
    "censor (alive=0; dead=1)": "censor",
    "radio_status (treated=1;un-treated=0)": "radio_status",
    "chemo_status (tmz treated=1;un-treated=0)": "chemo_status",
    "idh_mutation_status": "idh_mutation_status",
    "1p19q_codeletion_status": "codel_1p19q_status",
    "mgmtp_methylation_status": "mgmtp_methylation_status",
}

#: Columns of the GlioStrata cohort schema, in order.
COHORT_COLUMNS = [
    "patient_id",
    "age",
    "sex",
    "kps",
    "extent_of_resection",
    "mgmt_methylated",
    "idh_mutated",
    "subtype",
    "state_npc",
    "state_opc",
    "state_ac",
    "state_mes",
    "tme_immune",
    "tme_stromal",
    "tme_endothelial",
    "overall_survival_months",
    "vital_status",
    "treatment_arm",
]

#: Days per month, identical to scripts/download_tcga_clinical.py so the two
#: cohorts' survival columns are on the same scale.
DAYS_PER_MONTH = 30.4375

#: Label for cases where both radiotherapy and chemotherapy were recorded as
#: NOT given. A real, informative category, as in the TCGA ingest.
ARM_NO_ADJUVANT = "NoAdjuvant"

#: Histology codes that denote glioblastoma. CGGA uses a prefix convention:
#: "rGBM" = recurrent GBM, "sGBM" = secondary GBM (progressed from a
#: lower-grade glioma). Only bare "GBM" is a primary de novo glioblastoma, and
#: PRS_type carries the same information independently.
GBM_HISTOLOGY = ("GBM",)

#: WHO grade string for grade IV.
GRADE_IV = "WHO IV"

#: Deterministic order in which filter reasons are evaluated and reported.
FILTER_REASONS = [
    ("not_primary", "PRS_type is not Primary (recurrent or secondary)"),
    ("non_glioblastoma", "Histology is not GBM"),
    ("not_grade_iv", "Grade is not WHO IV"),
    ("no_usable_data", "no age, sex, survival or vital status"),
    ("nonpositive_survival", "overall_survival_months <= 0"),
]

#: Fields of the GlioStrata schema that CGGA simply does not record. Set to
#: pd.NA rather than imputed or approximated from something else.
UNMAPPABLE_FIELDS = {
    "kps": "CGGA publishes no Karnofsky or ECOG performance score.",
    "extent_of_resection": (
        "CGGA publishes no resection-extent field (no GTR/STR/biopsy column)."
    ),
    "state_npc": "Requires scRNA-seq deconvolution (project Phase C).",
    "state_opc": "Requires scRNA-seq deconvolution (project Phase C).",
    "state_ac": "Requires scRNA-seq deconvolution (project Phase C).",
    "state_mes": "Requires scRNA-seq deconvolution (project Phase C).",
    "tme_immune": "Requires TME deconvolution (project Phase C).",
    "tme_stromal": "Requires TME deconvolution (project Phase C).",
    "tme_endothelial": "Requires TME deconvolution (project Phase C).",
}

#: CGGA clinical fields that have NO column in the GlioStrata schema. They are
#: preserved in clinical_raw.parquet rather than dropped or forced into a
#: schema column that means something else.
UNSCHEMAD_CGGA_FIELDS = {
    "codel_1p19q_status": (
        "1p/19q codeletion. No schema column. Not interchangeable with "
        "idh_mutated; near-absent in grade IV anyway."
    ),
    "prs_type": "Primary/Recurrent/Secondary. Consumed by the filter, not stored.",
    "histology": "CGGA histology code. Consumed by the filter, not stored.",
    "grade": "WHO grade string. Consumed by the filter, not stored.",
    "radio_status": "Folded into treatment_arm; the raw 0/1 flag is not stored.",
    "chemo_status": "Folded into treatment_arm; the raw 0/1 flag is not stored.",
}


# --------------------------------------------------------------------------
# Download / cache
# --------------------------------------------------------------------------

def file_url(entry: dict) -> str:
    """Portal download URL for one CGGA_FILES entry."""
    return (
        f"{CGGA_DOWNLOAD_URL}?file={entry['file']}"
        f"&type={entry['type']}&time={entry['time']}"
    )


def zip_name(entry: dict) -> str:
    """Basename of the zip archive for one CGGA_FILES entry."""
    return Path(entry["file"]).name


def cache_path_for(cache_dir, entry: dict) -> Path:
    """Cache location of one downloaded archive."""
    return Path(cache_dir) / zip_name(entry)


def pending_keys(keys, cache_dir) -> list:
    """Keys whose archive is not already cached (resume support)."""
    pending = []
    for key in keys:
        path = cache_path_for(cache_dir, CGGA_FILES[key])
        if path.is_file() and path.stat().st_size > 0:
            continue
        pending.append(key)
    return pending


def download_one(entry: dict, cache_dir, timeout: int = 600) -> int:
    """Download one archive into the cache. Returns bytes written."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    url = file_url(entry)
    response = requests.get(url, timeout=timeout, allow_redirects=True)
    if response.status_code != 200:
        raise RuntimeError(
            f"CGGA download error {response.status_code} for {url}: "
            f"{response.text[:300]}"
        )
    content = response.content
    if not content.startswith(b"PK"):
        raise RuntimeError(
            f"CGGA returned non-zip content for {url} "
            f"(first bytes {content[:16]!r}) -- the portal may be serving an "
            "HTML error or consent page instead of the file"
        )
    path = cache_path_for(cache_dir, entry)
    path.write_bytes(content)
    return len(content)


def download_files(keys, cache_dir) -> dict:
    """Download every missing archive, reporting progress. ASCII-only."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    todo = pending_keys(keys, cache_dir)
    skipped = len(list(keys)) - len(todo)
    print(
        f"Download: {len(list(keys))} archives requested -> "
        f"{skipped} already cached, {len(todo)} to fetch"
    )
    started = time.time()
    total_bytes = 0
    for i, key in enumerate(todo, start=1):
        size = download_one(CGGA_FILES[key], cache_dir)
        total_bytes += size
        print(
            f"  [{i}/{len(todo)}] {key:26s} -> {size / 1e6:8.3f} MB "
            f"({time.time() - started:6.1f}s elapsed)",
            flush=True,
        )
    return {
        "requested": len(list(keys)),
        "skipped_cached": skipped,
        "downloaded": len(todo),
        "bytes": total_bytes,
        "elapsed_seconds": round(time.time() - started, 1),
    }


def extract_text_member(archive_path, extract_dir) -> Path:
    """Extract the single .txt member of a CGGA archive.

    CGGA zips carry a `__MACOSX/` resource-fork directory alongside the real
    payload; those members are skipped. The extracted bytes are written
    verbatim, so the cache is a complete audit trail of what the portal served.
    """
    archive_path = Path(archive_path)
    extract_dir = Path(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as zf:
        members = [
            n
            for n in zf.namelist()
            if n.lower().endswith(".txt")
            and not n.startswith("__MACOSX/")
            and not Path(n).name.startswith(".")
        ]
        if len(members) != 1:
            raise ValueError(
                f"{archive_path.name}: expected exactly one .txt member, "
                f"found {members}"
            )
        name = members[0]
        target = extract_dir / Path(name).name
        if not (target.is_file() and target.stat().st_size > 0):
            target.write_bytes(zf.read(name))
    return target


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def normalise_clinical_columns(columns) -> list:
    """Map CGGA clinical header spellings onto canonical names.

    CGGA headers vary between batches ("1p19q_codeletion_status" vs
    "1p19q_Codeletion_status") and embed the coding legend
    ("Censor (alive=0; dead=1)"). Unknown headers are lower-cased and kept
    rather than dropped, so a new CGGA column shows up in the raw table.
    """
    out = []
    for col in columns:
        key = str(col).strip().lower()
        out.append(CLINICAL_COLUMN_ALIASES.get(key, key))
    return out


def read_clinical(path, dataset: str) -> pd.DataFrame:
    """Read one CGGA clinical TSV into a normalised raw frame.

    CGGA's 2020-05-06 clinical files use classic-Mac CR-only line endings, so
    they must be read with universal newlines (pandas' C parser does this);
    `wc -l` reports 0 lines for them.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str)
    frame.columns = normalise_clinical_columns(frame.columns)
    # CGGA 693 ships trailing-whitespace variants ("Female ", "Male  ").
    # The dtype check covers both pandas 2 (object) and pandas 3 (str).
    for col in frame.columns:
        series = frame[col]
        if pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(
            series
        ):
            frame[col] = series.str.strip()
    frame.insert(1, "cgga_dataset", dataset)
    return frame


def read_counts(path) -> pd.DataFrame:
    """Read one CGGA read-counts matrix (gene SYMBOL x sample).

    The first column is `gene_name`: CGGA ships HGNC-style gene symbols, NOT
    Ensembl gene ids. See harmonisation_facts() for the mapping consequence.
    """
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str).str.strip()
    frame.index.name = "gene_name"
    if frame.index.duplicated().any():
        dupes = frame.index[frame.index.duplicated()].unique().tolist()
        raise ValueError(f"{path}: duplicate gene symbols {dupes[:5]}")
    return frame


def concat_clinical(frames) -> pd.DataFrame:
    """Concatenate per-batch clinical frames, deterministically ordered.

    Sorted by (cgga_dataset, cgga_id) so reruns are byte-identical. CGGA
    sample ids repeat across batches (the microarray batch shares 84 ids with
    mRNAseq_325 and 14 with mRNAseq_693), so rows are NOT deduplicated here --
    clinical_raw.parquet keeps every batch's own record.
    """
    if not frames:
        raise ValueError("no clinical frames to concatenate")
    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined = combined.sort_values(
        ["cgga_dataset", "cgga_id"], kind="mergesort", na_position="last"
    ).reset_index(drop=True)
    return combined


# --------------------------------------------------------------------------
# Schema mapping
# --------------------------------------------------------------------------

def map_sex(value):
    """CGGA Gender -> cohort sex (M/F).

    CGGA 693 contains trailing-whitespace variants ("Female ", "Male  "), so
    the value is stripped before matching. Anything else becomes pd.NA.
    """
    if not isinstance(value, str):
        return pd.NA
    key = value.strip().lower()
    if key in ("male", "m"):
        return "M"
    if key in ("female", "f"):
        return "F"
    return pd.NA


def map_binary_status(value, positive, negative):
    """Map a CGGA status string to 1 / 0 / pd.NA.

    Used for IDH ("Mutant"/"Wildtype") and MGMT promoter
    ("methylated"/"un-methylated"). Note that "un-methylated" contains
    "methylated" as a substring, so negatives are tested FIRST -- a
    substring-match-positive-first implementation silently codes every
    unmethylated case as methylated.
    """
    if not isinstance(value, str):
        return pd.NA
    key = value.strip().lower()
    if key in ("", "na", "nan", "unknown", "n/a", "."):
        return pd.NA
    if key in negative:
        return 0
    if key in positive:
        return 1
    return pd.NA


def map_idh(value):
    """CGGA IDH_mutation_status -> idh_mutated (1 = mutant)."""
    return map_binary_status(value, positive=("mutant",), negative=("wildtype",))


def map_mgmt(value):
    """CGGA MGMTp_methylation_status -> mgmt_methylated (1 = methylated)."""
    return map_binary_status(
        value,
        positive=("methylated",),
        negative=("un-methylated", "unmethylated"),
    )


def map_subtype(value):
    """CGGA TCGA_subtypes -> cohort subtype.

    Present ONLY in the mRNA_array_301 batch, where it is a precomputed
    Verhaak label (Classical / Mesenchymal / Neural / Proneural). It is carried
    through verbatim and is NOT a substitute for the project's own T1
    prediction -- it is the comparison label, so using it as a feature would
    be target leakage.
    """
    if not isinstance(value, str):
        return pd.NA
    key = value.strip()
    if key.lower() in ("", "na", "nan", "unknown", "n/a"):
        return pd.NA
    return key


def days_to_months(days):
    """Convert CGGA OS (DAYS) to months, using the TCGA ingest's divisor.

    Verified empirically, not assumed: the OS column of all three batches
    ranges 19..5159 with a median near 800, which is days, not months
    (5159 months would be 430 years).
    """
    numeric = pd.to_numeric(days, errors="coerce")
    return (numeric / DAYS_PER_MONTH).round(2)


def map_censor(value):
    """CGGA "Censor (alive=0; dead=1)" -> cohort vital_status.

    Verified against the shipped header legend: CGGA's 0/1 coding is already
    the project's convention (1 = deceased, 0 = alive), so this is an identity
    map, NOT an inversion. CGGA calls the column "Censor" while coding the
    EVENT, which is the trap this function exists to document.
    """
    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return pd.NA
    if numeric == 1:
        return 1
    if numeric == 0:
        return 0
    return pd.NA


def classify_arm(radio, chemo):
    """CGGA Radio_status + Chemo_status -> cohort treatment_arm.

    Returns:
        "ChemoRadiation" - both recorded as given.
        "MonoTherapy"    - exactly one recorded as given.
        "NoAdjuvant"     - both recorded as NOT given.
        pd.NA            - either flag missing, so the arm is unknown. A
                           missing radio flag with chemo=1 could be either
                           MonoTherapy or ChemoRadiation, and guessing would
                           invent an exposure.
    """
    r = pd.to_numeric(radio, errors="coerce")
    c = pd.to_numeric(chemo, errors="coerce")
    if pd.isna(r) or pd.isna(c):
        return pd.NA
    if r not in (0, 1) or c not in (0, 1):
        return pd.NA
    if r == 1 and c == 1:
        return "ChemoRadiation"
    if r == 1 or c == 1:
        return "MonoTherapy"
    return ARM_NO_ADJUVANT


def build_cohort(raw: pd.DataFrame) -> pd.DataFrame:
    """Map a normalised raw CGGA clinical frame to the GlioStrata schema.

    Row order and row count are preserved, so filter_cohort can align the
    result against `raw` positionally.
    """
    def col(name):
        if name in raw.columns:
            return raw[name]
        return pd.Series([None] * len(raw), index=raw.index, dtype="object")

    cohort = pd.DataFrame(index=raw.index)
    cohort["patient_id"] = col("cgga_id").astype("string")
    cohort["age"] = pd.to_numeric(col("age_years"), errors="coerce").round(0)
    cohort["sex"] = pd.Series(
        [map_sex(v) for v in col("gender")], index=raw.index, dtype="string"
    )

    for field in ("kps", "extent_of_resection"):
        cohort[field] = pd.NA

    cohort["mgmt_methylated"] = pd.Series(
        [map_mgmt(v) for v in col("mgmtp_methylation_status")],
        index=raw.index,
        dtype="Int64",
    )
    cohort["idh_mutated"] = pd.Series(
        [map_idh(v) for v in col("idh_mutation_status")],
        index=raw.index,
        dtype="Int64",
    )
    cohort["subtype"] = pd.Series(
        [map_subtype(v) for v in col("tcga_subtypes")],
        index=raw.index,
        dtype="string",
    )

    for field in (
        "state_npc",
        "state_opc",
        "state_ac",
        "state_mes",
        "tme_immune",
        "tme_stromal",
        "tme_endothelial",
    ):
        cohort[field] = pd.NA

    cohort["overall_survival_months"] = days_to_months(col("os_days"))
    cohort["vital_status"] = pd.Series(
        [map_censor(v) for v in col("censor")], index=raw.index, dtype="Int64"
    )
    cohort["treatment_arm"] = pd.Series(
        [
            classify_arm(r, c)
            for r, c in zip(col("radio_status"), col("chemo_status"))
        ],
        index=raw.index,
        dtype="string",
    )

    cohort["cgga_dataset"] = col("cgga_dataset").astype("string")
    return cohort[COHORT_COLUMNS + ["cgga_dataset"]]


# --------------------------------------------------------------------------
# Filtering
# --------------------------------------------------------------------------

def is_primary(value) -> bool:
    """True only for CGGA PRS_type == "Primary"."""
    return isinstance(value, str) and value.strip().lower() == "primary"


def is_gbm_histology(value) -> bool:
    """True only for bare "GBM".

    Deliberately NOT a substring match: "rGBM" (recurrent) and "sGBM"
    (secondary, progressed from lower-grade glioma) both contain "GBM" and
    both must be excluded from a primary-GBM external cohort.
    """
    if not isinstance(value, str):
        return False
    return value.strip().upper() in tuple(h.upper() for h in GBM_HISTOLOGY)


def is_grade_iv(value) -> bool:
    """True only for CGGA Grade == "WHO IV"."""
    if not isinstance(value, str):
        return False
    return " ".join(value.strip().upper().split()) == GRADE_IV.upper()


def filter_cohort(cohort: pd.DataFrame, raw: pd.DataFrame):
    """Restrict the mapped cohort to primary WHO grade IV glioblastoma.

    Reasons are evaluated in the fixed order of FILTER_REASONS. The masks
    overlap, so the report records both the rows matching each reason
    independently ("matched") and the rows attributed to it under first-match
    attribution ("dropped"), which sums exactly to "dropped_total" -- the same
    contract as scripts/download_tcga_clinical.py.

    Returns:
        (filtered cohort, report dict). Row order is preserved.
    """
    if len(cohort) != len(raw):
        raise ValueError(f"cohort/raw length mismatch: {len(cohort)} vs {len(raw)}")

    def flags(name, predicate):
        if name not in raw.columns:
            return pd.Series(False, index=cohort.index)
        return pd.Series(
            [predicate(v) for v in raw[name]], index=cohort.index, dtype=bool
        )

    survival = pd.to_numeric(cohort["overall_survival_months"], errors="coerce")
    masks = {
        "not_primary": ~flags("prs_type", is_primary),
        "non_glioblastoma": ~flags("histology", is_gbm_histology),
        "not_grade_iv": ~flags("grade", is_grade_iv),
        "no_usable_data": (
            cohort["age"].isna()
            & cohort["sex"].isna()
            & cohort["overall_survival_months"].isna()
            & cohort["vital_status"].isna()
        ),
        "nonpositive_survival": survival.notna() & (survival <= 0),
    }

    report = {"input_rows": int(len(cohort))}
    already = pd.Series(False, index=cohort.index)
    for key, _desc in FILTER_REASONS:
        mask = masks[key].fillna(False).astype(bool)
        report[f"matched_{key}"] = int(mask.sum())
        report[f"dropped_{key}"] = int((mask & ~already).sum())
        already = already | mask

    report["dropped_total"] = int(already.sum())
    filtered = cohort.loc[~already].copy()
    report["kept_rows"] = int(len(filtered))
    return filtered, report


def print_filter_report(report: dict) -> None:
    """Print an ASCII-only filter report to stdout."""
    print("\nCohort filter report:")
    print(f"  Input rows: {report['input_rows']}")
    for key, desc in FILTER_REASONS:
        print(
            f"  Dropped - {desc}: {report[f'dropped_{key}']}"
            f"  (rows matching this reason: {report[f'matched_{key}']})"
        )
    print(f"  Total unique rows dropped: {report['dropped_total']}")
    print(f"  Final cohort rows: {report['kept_rows']}")


def print_missingness(cohort: pd.DataFrame) -> None:
    """Column-by-column missingness of the mapped cohort. ASCII-only."""
    total = len(cohort)
    print(f"\nCohort missingness ({total} rows):")
    for col in cohort.columns:
        missing = int(cohort[col].isna().sum())
        pct = (100.0 * missing / total) if total else 0.0
        note = " <- never recorded by CGGA" if col in UNMAPPABLE_FIELDS else ""
        print(f"  {col:26s} missing {missing:5d} / {total:5d} ({pct:5.1f}%){note}")


def duplicate_id_report(cohort: pd.DataFrame) -> dict:
    """Cross-batch duplicate CGGA_IDs in the cohort.

    CGGA reuses sample ids across batches: the microarray batch shares ids
    with both RNA-seq batches. One id can therefore appear up to twice in the
    cohort (once as a microarray sample, once as an RNA-seq sample). This
    matters for the external partition, which must collapse duplicate patient
    ids before assignment so a patient cannot straddle a split
    (docs/preregistration.md Section 4.1).
    """
    ids = cohort["patient_id"].astype("string")
    counts = ids.value_counts()
    repeated = counts[counts > 1]
    pairs = {}
    if "cgga_dataset" in cohort.columns:
        for batch_a in sorted(cohort["cgga_dataset"].dropna().unique()):
            set_a = set(
                cohort.loc[cohort["cgga_dataset"] == batch_a, "patient_id"].dropna()
            )
            for batch_b in sorted(cohort["cgga_dataset"].dropna().unique()):
                if batch_b <= batch_a:
                    continue
                set_b = set(
                    cohort.loc[
                        cohort["cgga_dataset"] == batch_b, "patient_id"
                    ].dropna()
                )
                shared = len(set_a & set_b)
                if shared:
                    pairs[f"{batch_a} & {batch_b}"] = shared
    return {
        "rows": int(len(cohort)),
        "unique_patients": int(ids.nunique(dropna=True)),
        "ids_appearing_more_than_once": int(len(repeated)),
        "rows_beyond_first_per_id": int(len(cohort) - ids.nunique(dropna=True)),
        "shared_between_batches": pairs,
    }


def print_duplicate_id_report(report: dict) -> None:
    """ASCII-only duplicate-id report."""
    print("\nPatient-id uniqueness:")
    print(f"  Cohort rows: {report['rows']}")
    print(f"  Unique CGGA_IDs: {report['unique_patients']}")
    print(
        "  IDs appearing in more than one batch: "
        f"{report['ids_appearing_more_than_once']} "
        f"(extra rows: {report['rows_beyond_first_per_id']})"
    )
    for pair, shared in sorted(report["shared_between_batches"].items()):
        print(f"    {pair}: {shared} shared ids")
    print(
        "  Rows are NOT deduplicated here: the per-batch records differ "
        "(platform, and the microarray batch alone carries tcga_subtypes). "
        "Duplicate ids MUST be collapsed before partition assignment."
    )


def print_unmapped_report() -> None:
    """State explicitly what could not be mapped, rather than forcing it."""
    print("\nSchema fields CGGA cannot fill (left as pd.NA, not imputed):")
    for field, reason in UNMAPPABLE_FIELDS.items():
        print(f"  {field:24s} {reason}")
    print("\nCGGA fields with no schema column (kept in clinical_raw.parquet):")
    for field, reason in UNSCHEMAD_CGGA_FIELDS.items():
        print(f"  {field:24s} {reason}")


# --------------------------------------------------------------------------
# Expression assembly / harmonisation
# --------------------------------------------------------------------------

def assemble_counts(per_dataset: dict) -> pd.DataFrame:
    """Join per-batch counts matrices on the intersection of gene symbols.

    The RNA-seq batches are sample-disjoint, so columns simply concatenate.
    Genes are intersected (not unioned) and both axes are sorted, so reruns
    are byte-identical and no gene column is silently NaN-padded.
    """
    if not per_dataset:
        raise ValueError("no counts matrices to assemble")
    names = sorted(per_dataset)
    genes = None
    for name in names:
        idx = set(per_dataset[name].index)
        genes = idx if genes is None else (genes & idx)
    gene_order = sorted(genes)

    blocks = [per_dataset[name].loc[gene_order] for name in names]
    overlap = set()
    seen = set()
    for block in blocks:
        overlap |= seen & set(block.columns)
        seen |= set(block.columns)
    if overlap:
        raise ValueError(
            f"RNA-seq batches share sample ids {sorted(overlap)[:5]} -- "
            "a batch column would be silently duplicated"
        )

    counts = pd.concat(blocks, axis=1)
    counts = counts[sorted(counts.columns)]
    counts.index.name = "gene_name"
    return counts


def harmonisation_facts(counts: pd.DataFrame, tcga_annotation_path=None) -> dict:
    """Facts a cross-cohort comparison must respect, recorded in the output.

    Computed from gene identifiers and sample ids only -- no expression value
    is compared against any CGGA outcome.
    """
    cgga_symbols = set(counts.index.astype(str))
    facts = {
        "cgga_gene_identifier": "gene_name (HGNC-style SYMBOL)",
        "cgga_gene_count": int(counts.shape[0]),
        "cgga_sample_count": int(counts.shape[1]),
        "tcga_gene_identifier": "gene_id (Ensembl, GENCODE v36, versioned)",
        "batch_correction_required": True,
        "batch_correction_method": "ComBat-seq (protocol-specified)",
        "batch_variable": (
            "cgga_dataset for CGGA-internal batches (mRNAseq_693 vs "
            "mRNAseq_325), plus cohort (CGGA vs TCGA) for cross-cohort work"
        ),
        "batch_correction_notes": [
            "CGGA counts.parquet concatenates two independently processed "
            "RNA-seq batches (mRNAseq_693, mRNAseq_325); they are a batch "
            "variable even within CGGA.",
            "CGGA is a Chinese cohort and TCGA-GBM is predominantly European-"
            "ancestry and US-treated, so cohort differences are population and "
            "care-pathway differences, not only technical batch. ComBat-seq "
            "removes the technical component only.",
            "ComBat-seq operates on raw integer counts, which is why "
            "counts.parquet (not the FPKM/RSEM release) is the ingested "
            "matrix.",
            "Batch correction must be fitted on the development cohort and "
            "APPLIED to CGGA; fitting any correction jointly on CGGA and TCGA "
            "lets external-test information reach the training data.",
            "The microarray batch (mRNA_array_301) is clinical-only in this "
            "ingest: Agilent intensities are not counts and ComBat-seq does "
            "not apply to them.",
        ],
        "gene_id_mapping_strategy": (
            "Map TCGA Ensembl gene_id -> gene_name via the GENCODE v36 "
            "annotation already shipped in "
            "data/raw/tcga_gbm/expression/gene_annotation.parquet, then "
            "intersect on symbol. Symbol is the lossy side of the join: "
            "GENCODE gene_name is not unique, so a symbol-keyed join must "
            "resolve one-to-many Ensembl ids per symbol (restricting to "
            "gene_type == protein_coding resolves nearly all of them)."
        ),
    }
    if tcga_annotation_path is not None and Path(tcga_annotation_path).is_file():
        annot = pd.read_parquet(tcga_annotation_path)
        names = annot["gene_name"].astype(str)
        in_cgga = names.isin(cgga_symbols)
        protein = annot["gene_type"] == "protein_coding"
        facts["tcga_gene_count"] = int(len(annot))
        facts["tcga_duplicate_symbol_rows"] = int(names.duplicated().sum())
        facts["symbols_in_common"] = int(len(cgga_symbols & set(names)))
        facts["fraction_cgga_genes_mapped"] = round(
            len(cgga_symbols & set(names)) / max(len(cgga_symbols), 1), 4
        )
        facts["fraction_tcga_genes_mapped"] = round(float(in_cgga.mean()), 4)
        facts["tcga_protein_coding_genes"] = int(protein.sum())
        facts["tcga_protein_coding_mapped"] = int((protein & in_cgga).sum())
        facts["fraction_tcga_protein_coding_mapped"] = round(
            float(in_cgga[protein].mean()), 4
        )
    return facts


def print_harmonisation_report(facts: dict) -> None:
    """ASCII-only harmonisation report."""
    print("\nHarmonisation report (CGGA vs TCGA):")
    print(f"  CGGA gene identifier: {facts['cgga_gene_identifier']}")
    print(f"  CGGA genes x samples: {facts['cgga_gene_count']} x "
          f"{facts['cgga_sample_count']}")
    print(f"  TCGA gene identifier: {facts['tcga_gene_identifier']}")
    if "symbols_in_common" in facts:
        print(f"  TCGA genes: {facts['tcga_gene_count']} "
              f"(duplicate gene_name rows: {facts['tcga_duplicate_symbol_rows']})")
        print(f"  Symbols in common: {facts['symbols_in_common']}")
        print(f"  Fraction of CGGA genes that map: "
              f"{facts['fraction_cgga_genes_mapped']:.4f}")
        print(f"  Fraction of TCGA genes that map: "
              f"{facts['fraction_tcga_genes_mapped']:.4f}")
        print(f"  TCGA protein-coding mapped: "
              f"{facts['tcga_protein_coding_mapped']} / "
              f"{facts['tcga_protein_coding_genes']} "
              f"({facts['fraction_tcga_protein_coding_mapped']:.4f})")
    else:
        print("  TCGA gene annotation not found -> mapping fraction not computed")
    print(f"  Batch correction required: {facts['batch_correction_required']} "
          f"({facts['batch_correction_method']})")
    print(f"  Batch variable: {facts['batch_variable']}")
    for note in facts["batch_correction_notes"]:
        print(f"    - {note}")
    print(f"  Gene-id strategy: {facts['gene_id_mapping_strategy']}")


def expression_overlap(cohort: pd.DataFrame, counts) -> dict:
    """How many cohort patients have both expression and usable survival.

    Counts identifiers only. No expression value is related to any outcome.
    """
    ids = set(cohort["patient_id"].dropna().astype(str))
    expr_ids = set(str(c) for c in counts.columns) if counts is not None else set()
    survival = pd.to_numeric(cohort["overall_survival_months"], errors="coerce")
    has_surv = cohort["patient_id"][
        survival.notna() & cohort["vital_status"].notna()
    ]
    surv_ids = set(has_surv.dropna().astype(str))
    return {
        "cohort_patients": len(ids),
        "with_expression": len(ids & expr_ids),
        "with_survival": len(surv_ids),
        "with_expression_and_survival": len(ids & expr_ids & surv_ids),
        "expression_only_samples": len(expr_ids - ids),
    }


def print_overlap_report(report: dict) -> None:
    """ASCII-only expression/survival overlap report."""
    print("\nExpression and survival availability:")
    print(f"  Cohort patients: {report['cohort_patients']}")
    print(f"  With an expression column: {report['with_expression']}")
    print(f"  With survival time AND vital status: {report['with_survival']}")
    print(
        "  With BOTH expression and survival: "
        f"{report['with_expression_and_survival']}"
    )
    print(
        "  Expression columns not in the filtered cohort "
        f"(other grades/recurrences): {report['expression_only_samples']}"
    )


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def selected_keys(datasets, skip_expression: bool) -> list:
    """CGGA_FILES keys for the requested datasets, in registry order."""
    wanted = set(datasets)
    keys = []
    for key, entry in CGGA_FILES.items():
        if entry["dataset"] not in wanted:
            continue
        if skip_expression and entry["kind"] != "clinical":
            continue
        keys.append(key)
    return keys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/raw/cgga", help="Output directory")
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=["mRNAseq_693", "mRNAseq_325", "mRNA_array_301"],
        choices=["mRNAseq_693", "mRNAseq_325", "mRNA_array_301"],
        help="CGGA batches to ingest",
    )
    parser.add_argument(
        "--cache-dir", default=None, help="Download cache (default: <out>/_cache)"
    )
    parser.add_argument(
        "--tcga-annotation",
        default="data/raw/tcga_gbm/expression/gene_annotation.parquet",
        help="TCGA gene annotation, for the gene-id mapping fraction",
    )
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Skip cohort row filtering (produce the unfiltered mapped cohort)",
    )
    parser.add_argument(
        "--skip-expression",
        action="store_true",
        help="Download and map clinical data only",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report the planned downloads only; fetch nothing",
    )
    args = parser.parse_args()

    out_dir = Path(args.out)
    cache_dir = Path(args.cache_dir) if args.cache_dir else out_dir / "_cache"
    extract_dir = cache_dir / "extracted"

    keys = selected_keys(args.datasets, args.skip_expression)
    print("CGGA ingest -- EXTERNAL VALIDATION COHORT")
    print("  This script downloads and characterises CGGA only.")
    print("  No model is fitted and no feature-outcome association is computed.")
    print(f"\nDatasets requested: {', '.join(args.datasets)}")
    print(f"Archives planned: {len(keys)}")
    for key in keys:
        print(f"  {key:26s} -> {file_url(CGGA_FILES[key])}")

    if args.dry_run:
        print("\n--dry-run -> nothing downloaded")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    dl_report = download_files(keys, cache_dir)
    print(
        f"Fetched {dl_report['downloaded']} archives "
        f"({dl_report['skipped_cached']} cached), "
        f"{dl_report['bytes'] / 1e6:.3f} MB in {dl_report['elapsed_seconds']}s"
    )
    print(f"Cache directory (raw audit trail): {cache_dir}")

    clinical_frames = []
    counts_frames = {}
    for key in keys:
        entry = CGGA_FILES[key]
        text_path = extract_text_member(cache_path_for(cache_dir, entry), extract_dir)
        if entry["kind"] == "clinical":
            frame = read_clinical(text_path, entry["dataset"])
            print(
                f"  Parsed {key:26s} -> {frame.shape[0]:4d} rows x "
                f"{frame.shape[1]} cols"
            )
            clinical_frames.append(frame)
        else:
            frame = read_counts(text_path)
            print(
                f"  Parsed {key:26s} -> {frame.shape[0]:5d} genes x "
                f"{frame.shape[1]} samples"
            )
            counts_frames[entry["dataset"]] = frame

    raw = concat_clinical(clinical_frames)
    raw_path = out_dir / "clinical_raw.parquet"
    raw.to_parquet(raw_path, index=False)
    print(f"\nSaved raw clinical data (UNFILTERED): {raw_path}")
    print(f"  Shape: {raw.shape[0]} rows x {raw.shape[1]} cols")

    print("\nRaw field coverage (non-null count):")
    for col in raw.columns:
        print(f"  {col:28s} {int(raw[col].notna().sum()):5d} / {len(raw)}")

    print("\nRaw grade / histology breakdown (all batches pooled):")
    if "grade" in raw.columns:
        for value, count in sorted(
            raw["grade"].fillna("NA").value_counts().items()
        ):
            print(f"  grade {value:12s} {int(count):5d}")
    if "histology" in raw.columns:
        for value, count in sorted(
            raw["histology"].fillna("NA").value_counts().items()
        ):
            print(f"  histology {value:8s} {int(count):5d}")

    cohort = build_cohort(raw)
    if args.no_filter:
        print("\nCohort filter report: SKIPPED (--no-filter)")
        print(f"  Final cohort rows: {len(cohort)}")
    else:
        cohort, report = filter_cohort(cohort, raw)
        print_filter_report(report)

    cohort = cohort.sort_values(
        ["cgga_dataset", "patient_id"], kind="mergesort"
    ).reset_index(drop=True)
    cohort_path = out_dir / "clinical_cohort.parquet"
    cohort.to_parquet(cohort_path, index=False)
    print(f"\nSaved mapped cohort: {cohort_path}")
    print(f"  Shape: {cohort.shape[0]} rows x {cohort.shape[1]} cols")
    print("\nRows per batch:")
    for value, count in cohort["cgga_dataset"].value_counts().sort_index().items():
        print(f"  {value:16s} {int(count):5d}")

    print_duplicate_id_report(duplicate_id_report(cohort))
    print_missingness(cohort)
    print_unmapped_report()

    counts = None
    if counts_frames:
        counts = assemble_counts(counts_frames)
        counts_path = out_dir / "counts.parquet"
        counts.to_parquet(counts_path)
        annot = pd.DataFrame({"gene_name": counts.index.astype(str)})
        annot_path = out_dir / "gene_annotation.parquet"
        annot.to_parquet(annot_path, index=False)
        print(f"\nSaved counts matrix: {counts_path}")
        print(f"  Shape: {counts.shape[0]} genes x {counts.shape[1]} samples")
        print(f"Saved gene annotation: {annot_path}")

        facts = harmonisation_facts(counts, args.tcga_annotation)
        facts_path = out_dir / "harmonisation.json"
        facts_path.write_text(json.dumps(facts, indent=2), encoding="ascii")
        print(f"Saved harmonisation facts: {facts_path}")
        print_harmonisation_report(facts)
    else:
        print("\n--skip-expression -> no counts matrix assembled")

    print_overlap_report(expression_overlap(cohort, counts))

    print("\nDistributions (parse verification only -- no outcome association):")
    for col in ["age", "overall_survival_months"]:
        values = pd.to_numeric(cohort[col], errors="coerce").dropna()
        if values.empty:
            print(f"  {col:26s} no values")
            continue
        print(
            f"  {col:26s} n={len(values):4d} min={values.min():8.2f} "
            f"median={values.median():8.2f} max={values.max():8.2f}"
        )
    for col in ["sex", "vital_status", "treatment_arm", "idh_mutated",
                "mgmt_methylated", "subtype"]:
        counts_map = cohort[col].value_counts(dropna=False).to_dict()
        rendered = ", ".join(f"{k}={v}" for k, v in sorted(
            counts_map.items(), key=lambda kv: str(kv[0])
        ))
        print(f"  {col:26s} {rendered}")

    print("\nDone. CGGA is a frozen external test cohort: do not fit, tune or")
    print("evaluate anything on it outside the registered confirmatory run.")


if __name__ == "__main__":
    main()
