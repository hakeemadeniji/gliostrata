# CGGA ingestion

How the Chinese Glioma Genome Atlas (CGGA) is downloaded, mapped into the GlioStrata
cohort schema, and filtered to a primary-glioblastoma external validation cohort.

Implementation: `scripts/download_cgga.py`. Tests: `tests/test_cgga.py` (no network).

> **CGGA is the project's pre-registered external test set** for protocol tests T1 and T2
> (`docs/preregistration.md` §2, §3.3). This ingest downloads and characterises it only.
> Nothing in this pipeline fits, tunes or selects anything on CGGA data, and no model
> metric and no feature-outcome association is computed on it. The characterisation below
> is confined to row counts, column dtypes, missingness, value ranges and univariate
> marginal distributions — the checks needed to confirm the files parsed correctly. Looking
> at how any CGGA feature relates to CGGA survival would contaminate the test set
> irreversibly; retraining does not undo it.

---

## 1. What is actually downloadable

Established empirically, not from the portal's prose: an unauthenticated `GET` against
`https://www.cgga.org.cn/download?file=...&type=...&time=...` returns
`HTTP 200`, `Content-Type: application/zip`, and a real ZIP archive. **No registration,
no login, no cookie and no click-through is required** for any of the files below. The
portal's `download?file=` endpoint is a counter that 302-redirects to the static path
`https://www.cgga.org.cn/download/<date>/<name>.zip`; the script uses the counter endpoint
so CGGA's own download statistics stay accurate, and verifies the response begins with the
ZIP magic `PK` so an HTML error or consent page is never silently parsed as data.

The protocol's description of CGGA as ">1,000 samples, bulk RNA-seq + methylation +
clinical, **open with registration**" is correct on sample count and content, and wrong
about registration — these files are open without it.

### Batches ingested

| Batch | Samples | Clinical | Expression | Gene ids | Ingested |
| --- | --- | --- | --- | --- | --- |
| `mRNAseq_693` | 693 | yes, 13 cols | STAR read counts + RSEM FPKM | symbols | counts + clinical |
| `mRNAseq_325` | 325 | yes, 13 cols | STAR read counts + RSEM FPKM | symbols | counts + clinical |
| `mRNA_array_301` | 301 | yes, **14** cols | Agilent, gene + probe level | symbols | **clinical only** |

`mRNA_array_301` is the only batch with a `TCGA_subtypes` column (a precomputed Verhaak
label). Its expression is Agilent intensity, not counts, so ComBat-seq — a
negative-binomial method for count data — does not apply to it; it is therefore ingested
clinical-only. The two RNA-seq batches are **sample-disjoint** (0 shared `CGGA_ID`s).

### Other CGGA datasets, reachable but not ingested

Available from the same open endpoint, listed so nobody re-derives this:
`WESeq_286` (exome + clinical), `methyl_159` (Illumina 27K methylation + clinical),
`microRNA_198`, `scRNAseq_6148` (UMI matrix + clinical), `CCell_4083` (IDH-mutant
multi-omics: EPIC methylation, RSEM expression, proteomics, phosphoproteomics, Seurat
object), `MRI_268` imaging, and 20 non-glioma RNA-seq controls. The methylation release is
**27K and only 159 samples**, with its own clinical file, so it cannot supply MGMT status
for the RNA-seq cohort; MGMT comes from the clinical annotation instead (§3). Raw FASTQ
lives in NGDC/GSA BioProjects `PRJCA001747` (693) and `PRJCA001746` (325) and is a
separate, controlled request. The spatiotemporal transcriptome and WES datasets are listed
on the portal with **no file links at all** — a human must email the CGGA contact address
for those.

---

## 2. The number that matters: GBM-only n

The protocol's ">1,000 samples" is **all glioma grades**. The grade-IV subset is much
smaller, and the primary-GBM subset smaller still.

Pooled across the three ingested batches (1,319 clinical rows):

| | Rows |
| --- | --- |
| All rows | 1,319 |
| WHO II / III / IV / missing | 408 / 391 / 512 / 8 |
| Histology `GBM` (primary) | 333 |
| Histology `rGBM` (recurrent) | 138 |
| Histology `sGBM` (secondary, progressed from LGG) | 41 |
| **Final cohort (primary, `GBM`, WHO IV, positive survival)** | **333 rows / 305 unique patients** |
| of which RNA-seq (has an expression column) | **225** |
| **with BOTH expression and complete survival** | **218** |

Per batch: `mRNAseq_693` → 140, `mRNAseq_325` → 85, `mRNA_array_301` → 108.

**333 rows but 305 patients.** CGGA reuses `CGGA_ID`s across batches: 27 ids appear in both
`mRNA_array_301` and `mRNAseq_325`, and 1 in both `mRNA_array_301` and `mRNAseq_693`. The
per-batch records differ (different platform, and only the microarray batch carries
`tcga_subtypes`), so the ingest keeps both rows and reports the overlap rather than picking
one arbitrarily. **Duplicate ids must be collapsed before partition assignment**, which
`partition_cohort` already does (`docs/preregistration.md` §4.1) — but a downstream
consumer that groups by row rather than patient would double-count 28 patients.

For planning purposes the honest external-validation n is **225** if the cohort is
RNA-seq-only (the realistic case for T1/T2, which need expression), or **218** with
complete survival as well. Not >1,000.

---

## 3. Schema mapping, decision by decision

Target schema: `patient_id, age, sex, kps, extent_of_resection, mgmt_methylated,
idh_mutated, subtype, state_npc, state_opc, state_ac, state_mes, tme_immune, tme_stromal,
tme_endothelial, overall_survival_months, vital_status, treatment_arm`.

| Schema field | CGGA source | Decision |
| --- | --- | --- |
| `patient_id` | `CGGA_ID` | Verbatim. Not globally unique across batches — see §2. |
| `age` | `Age` | **Already years** (range 8–79 pooled), unlike GDC's days. Rounded. No unit conversion. |
| `sex` | `Gender` | `Male`→`M`, `Female`→`F`, case-insensitive. CGGA 693 ships trailing-whitespace variants (`"Female "`, `"Male  "`, 10 rows), so values are stripped; without that they would become `pd.NA`. |
| `kps` | — | **No equivalent.** `pd.NA`. See §4. |
| `extent_of_resection` | — | **No equivalent.** `pd.NA`. See §4. |
| `mgmt_methylated` | `MGMTp_methylation_status` | `methylated`→1, `un-methylated`→0. **Exact match, negatives tested first:** `"un-methylated"` *contains* `"methylated"`, so a substring match would code every unmethylated case as methylated. Regression-tested. |
| `idh_mutated` | `IDH_mutation_status` | `Mutant`→1, `Wildtype`→0. Anything else `pd.NA` — never defaulted to wild-type, even though wild-type is the GBM norm. |
| `subtype` | `TCGA_subtypes` | Verbatim (`Classical`/`Mesenchymal`/`Neural`/`Proneural`). **Only in `mRNA_array_301`**, so 225/333 rows are null. This is the *comparison label* for T1, so using it as a model feature would be target leakage. |
| `state_*`, `tme_*` | — | **No equivalent.** `pd.NA`; require deconvolution (Phase C). |
| `overall_survival_months` | `OS` | **CGGA `OS` is in DAYS — verified, not assumed:** pooled range 19–5,159 with median ≈ 800. 5,159 months would be 430 years. Converted by `/ 30.4375`, the identical divisor used in `scripts/download_tcga_clinical.py`, so the two cohorts are on one scale. |
| `vital_status` | `Censor (alive=0; dead=1)` | **Identity, NOT an inversion — verified from the shipped header legend.** CGGA names the column "Censor" while coding the *event*: 1 = dead, 0 = alive, which already matches this project's `vital_status` (1 = deceased). Inverting it would silently reverse every survival analysis, so the mapping has its own test. |
| `treatment_arm` | `Radio_status` + `Chemo_status` | Both 1 → `ChemoRadiation`; exactly one → `MonoTherapy`; both 0 → `NoAdjuvant`; **either flag missing → `pd.NA`**. A missing radio flag with chemo=1 is genuinely ambiguous between `MonoTherapy` and `ChemoRadiation`, and guessing would invent an exposure. `NoAdjuvant` is a real category, matching the TCGA ingest. |

### Header normalisation

CGGA headers are not consistent across batches: the 1p/19q column is
`1p19q_codeletion_status` in the RNA-seq batches and `1p19q_**C**odeletion_status` in the
microarray batch, and three headers embed their coding legend
(`Censor (alive=0; dead=1)`). `normalise_clinical_columns` canonicalises the known
spellings case-insensitively and lower-cases unknown ones rather than dropping them, so a
new CGGA column appears in `clinical_raw.parquet` instead of vanishing.

### File-format trap

The 2020-05-06 clinical files use **classic-Mac CR-only (`\r`) line endings**. `wc -l`
reports **0 lines** for them. Pandas' C parser handles this; a hand-rolled
`open().readlines()` reader returns one giant line. The archives also carry a
`__MACOSX/` resource-fork directory, whose members are skipped on extraction.

---

## 4. What cannot be mapped

Reported rather than forced:

| Field | Why |
| --- | --- |
| `kps` | **CGGA publishes no performance score** — no Karnofsky, no ECOG, in any of the three clinical files. 100% null. TCGA-GBM does have KPS, so any model that uses `kps` as a feature cannot be externally validated on CGGA at all; it must be either dropped from the T1/T2 feature set or handled by an imputation strategy fixed *before* the external read. **This is a design decision for the lead, not for this ingest.** |
| `extent_of_resection` | No resection-extent column (no GTR/STR/biopsy field). 100% null — same consequence as `kps`. |
| `state_npc`, `state_opc`, `state_ac`, `state_mes` | Require scRNA-seq reference deconvolution (Phase C). CGGA's own scRNA-seq release (`scRNAseq_6148`) is 16 **IDH-mutant** samples, so it is not a GBM-state reference. |
| `tme_immune`, `tme_stromal`, `tme_endothelial` | Require TME deconvolution (Phase C). |
| `subtype` for 225/333 rows | Only the microarray batch carries a Verhaak label. |

CGGA fields with **no schema column**, preserved in `clinical_raw.parquet` rather than
coerced into a column that means something else:

* `1p19q_codeletion_status` — no schema column. Not interchangeable with `idh_mutated`
  (they are distinct markers) and near-absent in grade IV anyway.
* `PRS_type`, `Histology`, `Grade` — consumed by the filter, not stored in the cohort.
* `Radio_status`, `Chemo_status` — folded into `treatment_arm`; the raw 0/1 flags are not
  stored, so `NoAdjuvant` vs `pd.NA` is the only surviving distinction.

---

## 5. Filtering

Same contract as `scripts/download_tcga_clinical.py`: a fixed reason order, both an
independent match count and a first-match-attribution drop count per reason (the drops sum
exactly to the total), the raw download kept unfiltered, and a `--no-filter` flag that
emits the full 1,319-row mapped cohort.

```
Cohort filter report:
  Input rows: 1319
  Dropped - PRS_type is not Primary (recurrent or secondary): 404  (rows matching this reason: 404)
  Dropped - Histology is not GBM: 582  (rows matching this reason: 986)
  Dropped - Grade is not WHO IV: 0  (rows matching this reason: 807)
  Dropped - no age, sex, survival or vital status: 0  (rows matching this reason: 0)
  Dropped - overall_survival_months <= 0: 0  (rows matching this reason: 0)
  Total unique rows dropped: 986
  Final cohort rows: 333
```

Two details worth reading:

* **`rGBM` and `sGBM` both contain the substring `GBM`.** The histology test is an exact
  match on `GBM`, so recurrent (138 rows) and secondary (41 rows) glioblastoma are
  excluded. A substring match would quietly admit 179 non-primary tumours — including
  secondary GBMs, which are biologically IDH-mutant progressed LGGs with a very different
  natural history.
* **`not_grade_iv` drops 0.** Not a bug: `PRS_type` and `Histology` already remove every
  non-grade-IV row, so under first-match attribution grade IV has nothing left to drop.
  It stays in the report as an independent check — its *match* count is 807 — so a future
  CGGA release that puts a `GBM` histology on a non-grade-IV row would be caught.
* **`nonpositive_survival` drops 0.** CGGA's minimum `OS` is 19 days; no zero or negative
  survival exists in this release. The reason is kept so a future release cannot slip one
  through.

---

## 6. Harmonisation: what must happen before any cross-cohort comparison

Written to `data/raw/cgga/harmonisation.json` by the ingest, and reproduced here.

### 6.1 Gene identifiers: symbols vs Ensembl

| | CGGA | TCGA (this repo) |
| --- | --- | --- |
| Identifier | `gene_name`, HGNC-style **symbol** | `gene_id`, **Ensembl**, GENCODE v36, versioned |
| Genes | 55,523 | 60,660 |

An id-mapping step is mandatory. **Measured mapping fractions** (identifier sets only — no
expression value was compared against any outcome):

* Symbols in common: **33,212**
* Fraction of CGGA's 55,523 symbols that map to a TCGA `gene_name`: **0.5982**
* Fraction of TCGA's 60,660 rows whose `gene_name` is in CGGA: **0.5485**
* **Protein-coding: 17,988 / 19,962 = 0.9011**

The ~55% whole-transcriptome figure looks alarming and mostly is not: the loss is
concentrated in pseudogenes, lncRNAs and small RNAs, where symbol nomenclature has churned
hardest between CGGA's older annotation and GENCODE v36. For a protein-coding feature
space — which is what Verhaak subtyping and a survival model use — **90% maps**.

**Strategy:** map TCGA `gene_id` → `gene_name` using the GENCODE v36 annotation already in
`data/raw/tcga_gbm/expression/gene_annotation.parquet`, then intersect on symbol.
Symbol is the lossy side of the join: GENCODE `gene_name` is **not unique** (1,233
duplicate-name rows in the TCGA annotation), so a symbol-keyed join must resolve
one-to-many Ensembl ids per symbol. Restricting to `gene_type == protein_coding` resolves
nearly all of them. If the final feature set needs the non-coding transcriptome, the
remaining ~40% of CGGA symbols need a proper alias resolution (HGNC previous/alias symbols
or an Ensembl ID history lookup), not a plain string join — and that work must be done and
frozen *before* the confirmatory read.

### 6.2 Batch correction (protocol: ComBat-seq)

Required. Four points, in order of how easily they are missed:

1. **CGGA is itself two batches.** `counts.parquet` concatenates `mRNAseq_693` and
   `mRNAseq_325`, processed independently. `cgga_dataset` is carried in
   `clinical_cohort.parquet` precisely so it can be used as a batch covariate *within*
   CGGA. Treating CGGA as one homogeneous cohort is already wrong before TCGA enters.
2. **Counts, not FPKM.** ComBat-seq is a negative-binomial method requiring raw integer
   counts, which is why the ingest takes the STAR read-counts release and not the RSEM
   FPKM release. The FPKM files remain available if a non-count method is ever chosen.
3. **Fit on development, apply to CGGA.** Estimating correction parameters jointly over
   CGGA and TCGA lets external-test information reach the training data. That is leakage
   even though no outcome is involved, and `docs/preregistration.md` §3.2 already lists
   with/without harmonisation as a pre-specified sensitivity analysis — which only works
   if the correction is a transform that can be applied, not a jointly fitted one.
4. **Not all of the difference is technical.** CGGA is a Chinese cohort treated in Chinese
   centres; TCGA-GBM is predominantly European-ancestry and US-treated. ComBat-seq removes
   the technical component. Population and care-pathway differences remain and are a
   limitation to state, not a batch effect to correct away. Note also the age distribution:
   the CGGA primary-GBM cohort's median age is 52 against the TCGA-GBM cohort's 59 — a
   seven-year shift, which is a different population, not a different protocol.
5. **Microarray excluded.** `mRNA_array_301` contributes clinical rows only; Agilent
   intensities are not counts. Using its 108 patients in an expression model would need a
   separate cross-platform harmonisation that is not in scope here.

---

## 7. Outputs and usage

```
data/raw/cgga/
  _cache/<archive>.zip            download cache; resume skips what is present
  _cache/extracted/<file>.txt     unpacked raw text, byte-identical to source
  clinical_raw.parquet            1,319 rows x 15 cols, UNFILTERED, all batches
  clinical_cohort.parquet         333 rows x 19 cols, mapped + filtered
  counts.parquet                  55,523 genes x 1,018 samples (both RNA-seq batches)
  gene_annotation.parquet         gene_name
  harmonisation.json              the facts in §6
```

`counts.parquet` holds **all 1,018** RNA-seq samples, not only the 225 cohort ones — the
raw download stays a complete audit trail and filtering applies to the derived cohort, as
in the TCGA ingest. A consumer must subset it by `clinical_cohort.parquet`'s `patient_id`.

```bash
python scripts/download_cgga.py --dry-run      # report planned downloads, fetch nothing
python scripts/download_cgga.py                # full ingest (~45 MB, resumable)
python scripts/download_cgga.py --no-filter     # unfiltered mapped cohort (1,319 rows)
python scripts/download_cgga.py --skip-expression
python scripts/download_cgga.py --datasets mRNAseq_693
```

Reruns are byte-identical: clinical rows sort by `(cgga_dataset, cgga_id)`, counts sort on
both axes, and genes are intersected rather than unioned so no column is NaN-padded.

---

## 8. Open questions for the lead

Flagged rather than resolved here, because resolving any of them by looking at CGGA
outcomes would contaminate the test set.

1. **`kps` and `extent_of_resection` are 100% missing in CGGA.** If the T1/T2 models use
   either as a feature, external validation is impossible as specified. The feature set or
   the imputation rule must be fixed *before* the confirmatory read, and recorded as a
   pre-registration amendment if it changes anything registered.
2. **Which cohort definition is the registered one?** 225 (RNA-seq only), 218 (RNA-seq with
   complete survival), 305 (unique patients including microarray) or 333 (rows). This must
   be chosen and sealed in `cgga_external_v1` before anything is read. `docs/preregistration.md`
   §4.5 still lists the partition as not created.
3. **Is `mRNA_array_301` in or out?** It is the only source of a `TCGA_subtypes`
   comparison label for T1, but it is a different platform and shares 28 patients with the
   RNA-seq batches. Including it mixes platforms; excluding it removes the precomputed
   subtype label.
4. **Secondary GBM (`sGBM`, 41 rows) is excluded.** Defensible — IDH-mutant progressed
   LGGs — but it is a cohort-definition choice that belongs in the pre-registration.
5. **Non-coding gene mapping.** Should the feature space be restricted to protein-coding
   (90% mapped) or should a full HGNC alias resolution be built first?
6. **Does the 27K methylation batch (159 samples) matter?** It has its own clinical file
   and does not cover the RNA-seq cohort, so MGMT status here comes from the clinical
   annotation, not from methylation data.
