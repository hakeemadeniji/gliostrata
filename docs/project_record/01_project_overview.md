# 1. Project overview

## 1.1 What GlioStrata is

GlioStrata is an open-source research framework that aims to recommend, for an individual
glioblastoma patient, the **stem-cell-delivered therapeutic configuration** (carrier cell
type, payload, dose, injection sites and schedule) most likely to extend survival beyond
the Stupp protocol. It combines public multi-omic GBM cohorts with mechanistic tumour
simulation and Bayesian configuration optimisation.

It is run as a Master's thesis project (Health Informatics, Boston University), with a
computational deliverable and an optional wet-lab extension. Sources: `README.md`,
`docs/protocol.md`.

## 1.2 Research question and hypotheses

**Central question.** Given a patient's tumour multi-omic profile and single-cell state
composition, can a computational framework predict a stem-cell-delivered configuration
that yields a clinically meaningful survival gain over the Stupp baseline?

| Hypothesis | Statement | Registered measure |
| --- | --- | --- |
| H1 (primary) | Patient-specific selection gives a predicted restricted-mean-survival-time gain of at least 3 months in at least half of held-out patients, with a predicted harm rate below 10%. | T5 success fraction >= 0.50 (CI lower); harm fraction <= 0.10 (point) with CI upper <= 0.15 |
| H0 (null) | No configuration chosen by the optimizer reliably beats Stupp. A pre-registered negative result is publishable. | H1 not met |
| H2 (personalisation) | The top five recommendations across 50 held-out patients include at least three distinct payload/site combinations. | T4 distinct configurations >= 3 |
| H3, H4, H5 | The profiling (T1), survival-prediction (T2) and simulation (T3) components each clear their own external threshold. | T1, T2, T3 |

Source: `docs/preregistration.md` section 1.

## 1.3 The eleven registered thresholds and their status

A threshold is met only if the point estimate **and** its decisive 95% confidence bound
both clear it.

| Test | Metric | Requirement | Status (2026-11-26) |
| --- | --- | --- | --- |
| T1 | Subtyping accuracy, CGGA external | >= 0.88 | Not run |
| T1 | State deconvolution, Pearson r | >= 0.75 | Not run |
| T2 | Survival C-index, CGGA external | >= 0.65 | **FAIL**: 0.5385 (CI 0.4994 to 0.5807) |
| T2 | Improvement over a 2024 baseline | >= 0.03 | **Not evaluated**: no comparator named |
| T2 | Brier score at 18 months | <= 0.25 | **FAIL**: 0.2617 (CI 0.2391 to 0.2829) |
| T3 | Tumour-growth RMSE vs PDX | <= 15% | Not run; input data not obtained |
| T3 | NSC-migration r vs in vivo | >= 0.65 | Not run; input data not obtained |
| T4 | Optimizer beats random search | >= 25% regret reduction, Wilcoxon p < 0.01 | Not run |
| T4 | Distinct configurations in top-5 | >= 3 | Not run |
| T5 | End-to-end success fraction | >= 0.50 | Not run |
| T5 | Safety floor, predicted harm | <= 0.10 (point), <= 0.15 (CI upper) | Not run |

Two wet-lab tests, T6 and T7, are outside the pre-registration and have not been started.
Source: `docs/preregistration.md` section 2; `src/gliostrata/evaluation/thresholds.py`.

### How the thresholds are protected

- The table lives in `src/gliostrata/evaluation/thresholds.py`. A SHA-256 over its
  canonical form is asserted by a test, so any change to a number, a direction or a
  decisive bound fails the build. The frozen checksum is
  `96bbcd49dcda20f8a0bc9c8b71427aaa24cf754dc6d8d7b4f7d6a2df3d2874ab`.
- Evaluating a metric that is not registered raises an error; a more favourable metric
  cannot be substituted after results are seen.
- **Which results went through it.** The T2 Brier verdict was computed by the registry. The
  T2 C-index script applies the same registered rule (point estimate and lower bound both
  >= 0.65) with its own comparison and does not call the registry; on these numbers the
  registry gives the same verdict, FAIL. The molecular-calling and sensitivity results in
  `02_findings.md` are descriptive and are not threshold evaluations.

## 1.4 The pipeline and the status of each layer

The pipeline has eight layers. The table says what each does **today**. "Working" means it
does what it claims and is tested. "Implemented, not wired" means real, tested code that
the `gliostrata` command does not call. "Stub" means it returns synthetic or schema-only
output.

| Layer | Purpose | Module | Status today |
| --- | --- | --- | --- |
| L0 | Ingestion | `scripts/download_tcga_*.py`, `scripts/download_cgga.py` | **Working.** Real downloads from the GDC API and CGGA. |
| L1 | Harmonisation | `data/harmonize.py`, `data/qc.py` | Stub. Schema normalisation and QC rows; no batch correction. |
| L2 | Tumour profiling | `profiling/verhaak.py`, `profiling/deconvolve.py` | **Implemented, not wired.** Reads bulk expression. |
| L2 | Tumour profiling | `profiling/subtype.py`, `states.py`, `tme.py` | Stub. Derives calls from the clinical table; this is what the CLI runs. |
| L3 | Feature store | `features/store.py` | Working, but it consumes the stub L2 outputs. |
| L4 | Simulation | `simulation/tumor.py`, `simulation/nsc.py` | Stub. Synthetic curves. |
| L5 | Outcome prediction | `prediction/cox.py` | **Implemented, not wired.** A fitted Cox model; failed T2 externally. |
| L5 | Outcome prediction | `prediction/survival.py` | Stub. This is what the CLI runs. |
| L6 | Optimiser | `optimization/bayesian.py`, `config_space.py` | Working (ranked enumeration, not a Bayesian loop). |
| L6 | Optimiser | `optimization/bayes_opt.py` | **Implemented, not wired.** Gaussian-process / expected-improvement. |
| L7 | Wet-lab handoff | `handoff/protocol.py` | Working. JSON protocol documents; not bench-executable. |
| Evaluation | Pre-registration | `evaluation/thresholds.py`, `partitions.py` | **Working.** Registry, sealed partitions, drift detection. |

The most important consequence: the L0 expression ingest is real and produces count and
TPM matrices, but **nothing downstream in the pipeline reads them**. Every subtype call,
cell-state fraction and TME fraction the pipeline writes today is computed from clinical
covariates, not measured. Source: `docs/roadmap.md`, `docs/reproducibility.md` section 6.

## 1.5 Data in hand

| Cohort | What was obtained | Size |
| --- | --- | --- |
| TCGA-GBM | Clinical cohort | 617 case records, **594 eligible patients** |
| TCGA-GBM | Bulk RNA-seq | **284 samples** (STAR read counts and TPM) |
| TCGA-GBM | IDH mutation calls (from 464 open MAF files) | 373 patients tested; **368 of 594** cohort patients have a call |
| TCGA-GBM | MGMT promoter methylation calls (from 450 beta-value files) | 422 callable; **415 of 594** cohort patients have a call |
| CGGA | Clinical rows, three batches pooled | 1,319 rows, of which **333 are primary GBM**, which are **305 unique patients** |
| CGGA | RNA-seq read counts (two batches) | 76.9 MB matrix; **225 patients** have an expression profile |

Of the 305 CGGA patients, **218** have both an expression profile and a complete survival
outcome and are therefore eligible for the sealed external partition. The CGGA microarray
batch (301 samples) was ingested for its clinical rows only.

Sources: `docs/preregistration.md` section 5, `docs/mgmt_idh_calling.md` section 6,
`docs/cgga_ingestion.md` section 2, and a file listing of `data/raw/` on 2026-11-26.

## 1.6 Sealed partitions

A partition assigns each patient to train or test by a hash of (seed, patient_id), so it
reproduces without shuffling. The manifest is sealed and version-controlled, and a
verification step re-derives every assignment.

| Manifest | Cohort | Patients | Train / test | Seed | Status |
| --- | --- | --- | --- | --- | --- |
| `tcga_gbm_clinical_v2` | TCGA-GBM | 594 | 425 / 169 | 20251126 | **Current** development split. Verifies against the live cohort. |
| `tcga_gbm_clinical_v1` | TCGA-GBM | 583 | 419 / 164 | 20251126 | Superseded. Kept unmodified as the historical record; does not verify against the current cohort. |
| `cgga_external_v1` | CGGA | up to 218 | all test | to be set | **Not created.** |

Both TCGA partitions are unstratified, because the cohort's subtype column is entirely
null. v1 no longer verifies because the ingest filter was corrected and the cohort grew by
11 patients (583 to 594); v2 was frozen as a new manifest rather than by overwriting v1.
Source: manifests under `data/partitions/`; `docs/reproducibility.md` section 4.

## 1.7 Software facts (measured 2026-11-26)

| Property | Value |
| --- | --- |
| Package | `gliostrata` 0.1.0, MIT licence, Python >= 3.10 (developed on 3.14.3) |
| Test suite | **1016 tests in 26 modules: 1015 pass, 1 skipped, 0 failed**, with all outbound sockets blocked |
| The one skip | A Snakemake-dependent workflow test; Snakemake is an optional dependency and was not installed |
| Real-data tests | Skip when the TCGA or CGGA files are absent; they ran here because the data are present |
| Blocking lint (`E9,F63,F7,F82`) | Clean across `src`, `tests`, `scripts` |
| Advisory lint (full project config) | 171 findings: 31 in `src`, 133 in `tests`, 7 in `scripts`. 167 are naming rules that flag names such as `X` and `N`, which are conventional in statistical code |
| Pipeline definition | 9 DVC stages in `dvc.yaml`; 9 Snakemake rule files under `workflow/rules/` |
| Continuous integration | 4 jobs: lint, tests, partition-drift guardrail, Snakemake dry run with an offline fixture run |
| Dependency floors only | No lockfile, no container image |

The advisory lint count rose from 148 to 171 on 2026-11-26 because two new files were
added (20 and 3 findings). They follow the naming convention of their sibling files.

## 1.8 State of the pre-registration

| Item | State |
| --- | --- |
| Registered text | `docs/preregistration.md`, drafted and **not deposited** |
| Paste-ready OSF text | `docs/osf_final_fields.md` (four Sampling-section fields) |
| Threshold registry | In code, checksum-frozen, tested |
| Author details, DOI, registered commit | **Placeholders.** To be completed at deposit |
| Comparator model for T2's improvement row | **Placeholder.** "TO BE COMPLETED BY AUTHOR" |
| Amendments logged | Five dated entries, all 2026-11-26. Four are marked as written after CGGA was read; one (the T1 harmonisation rule) is marked "No" |
| Machine-readable registration state | `configs/registration.yaml`. The OSF DOI is `null` until the deposit; the read guard reads it |
| Access ledger | `data/access_ledger.jsonl`. 7 entries; 2 are value reads, both T2 and both before deposit |

Until the deposit is made the project has no pre-registration in the formal sense, and
T1 to T5 are exploratory. Source: `docs/preregistration.md` header and section 8.
