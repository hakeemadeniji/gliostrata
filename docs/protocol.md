# GlioStrata — Research Protocol (Markdown Edition)

> This is the markdown edition of the Master's Thesis Research Protocol.
> For the formatted .docx and .pdf editions (with cover page, three-line tables,
> and embedded architecture figure), see the `download/` folder at the project root.

---

## Abstract

Glioblastoma (GBM) is the most common and most lethal primary malignant brain tumor
in adults, with median survival near fifteen months despite maximal standard-of-care
resection, radiotherapy, and temozolomide. Stem-cell-delivered therapy — using neural
stem cells or mesenchymal stem cells as tumor-tropic vehicles for oncolytic viruses,
prodrug converting enzymes, or immunomodulatory payloads — has shown promise in
early-phase clinical trials but has been limited by a uniform, fixed-payload design
applied across molecularly heterogeneous tumors. No existing framework predicts, for
an individual patient's tumor, the optimal stem-cell-delivered therapeutic configuration.

This protocol proposes **GlioStrata**, a patient-specific computational framework that
integrates public multi-omic GBM data (TCGA, CGGA, CPTAC, single-cell atlases, GLASS
longitudinal, TCIA radiomics) with mechanistic tumor simulation and Bayesian
configuration optimization to recommend, per patient, the cell type, payload, dose,
injection sites, and schedule most likely to extend survival beyond the Stupp protocol
baseline. The framework is implemented as an eight-layer pipeline (L0–L7): data
ingestion, harmonization, tumor profiling and stratification, unified feature store,
mechanistic simulation, outcome prediction, therapeutic configuration optimization,
and wet-lab handoff with recalibration.

Validation follows a pre-registered, frozen-test-set protocol: training on TCGA-GBM
and external validation on CGGA and GLASS held-out cohorts, with success defined as a
predicted restricted mean survival time (RMST) gain of at least three months over
Stupp in at least fifty percent of held-out patients, with the lower bound of the
bootstrap ninety-five percent confidence interval clearing that threshold. A
ten-point defensibility stack — pre-registration on the Open Science Framework, frozen
test partitions, pre-specified thresholds, bootstrap intervals, multiple-comparison
correction, external validation, calibration, decision-curve analysis, sensitivity
analyses, and full reproducibility release — wraps every test.

The project delivers an open-source pipeline and curated multi-omic GBM feature
store, targeting the Nucleic Acids Research Database issue as the primary publication
venue. A wet-lab extension, in which a collaborating laboratory tests the framework's
top-ranked configurations on a patient-derived organoid panel, runs as a parallel
track and feeds results back into the model via Bayesian recalibration.

**Keywords**: glioblastoma; stem-cell therapy; neural stem cells; patient-specific
modeling; Bayesian optimization; multi-omics integration; tumor subtyping; survival
prediction; physiCell simulation; open-source pipeline; reproducibility; pre-registration

---

## Chapter 1 — Introduction

### 1.1 Project Overview

GlioStrata is a patient-specific computational framework for optimizing stem-cell-
delivered therapy in glioblastoma. The project reframes what is typically a wet-lab-
first therapeutic question — what cell line, what payload, what dose, what delivery
schedule to use — into a computational optimization problem whose output is a ranked
configuration set for each individual patient, together with predicted outcomes and
a wet-lab-ready protocol document that a collaborating laboratory can execute.

### 1.2 Scientific Background and Significance

Glioblastoma is the most common primary malignant brain tumor in adults and carries
one of the worst prognoses in oncology. Despite maximal safe resection followed by
concurrent radiotherapy and temozolomide — the Stupp protocol that has been the
standard of care since 2005 — median overall survival remains approximately fifteen
months, and five-year survival is below seven percent.

Stem-cell-delivered therapy offers a mechanistically distinct approach. Neural stem
cells and mesenchymal stem cells display intrinsic tumor-tropism — they migrate
preferentially toward GBM in response to chemoattractant gradients including SDF-1,
VEGF, uPA, and MCP-1 — and can be engineered to carry therapeutic payloads.

### 1.3 Research Question and Hypothesis

The central research question is: given a patient's tumor multi-omic profile and
single-cell state composition, can a computational framework predict a stem-cell-
delivered therapeutic configuration that yields a clinically meaningful survival
gain over the Stupp protocol baseline?

### 1.4 Project Aims

1. Build a harmonized multi-omic GBM feature store integrating bulk genomics,
   transcriptomics, proteomics, single-cell state fractions, radiomics, and
   longitudinal clinical outcomes.
2. Train and externally validate a tumor profiling pipeline (Verhaak subtype + Neftel
   state fractions + TME) with pre-registered accuracy thresholds.
3. Calibrate a mechanistic simulation engine (PhysiCell tumor growth + NSC migration
   + payload kinetics) against published PDX growth curves and in vivo NSC tracking.
4. Couple simulation to a Bayesian configuration optimizer and external survival
   prediction model, producing per-patient therapeutic recommendations whose predicted
   survival gain clears the pre-registered success threshold.

---

## Chapter 2 — Scientific Background

### 2.1 Current Treatment Landscape

The Stupp protocol (2005) established concomitant and adjuvant temozolomide with
focal radiotherapy as standard of care, improving median survival from 12.1 to 14.6
months. Tumor-treating fields (Stupp et al. 2017) extended survival to ~20.9 months
in selected patients. No systemic therapy has improved overall survival in phase III
trials for unselected GBM patients beyond this.

### 2.2 Stem-Cell-Delivered Therapy

Aboody et al. (2000) first reported NSC tumor-tropism. Subsequent work extended the
platform to cytosine deaminase / 5-FC prodrug conversion, NSC-CRAd oncolytic viruses
(City of Hope), and MSC-delivered payloads. Phase I clinical translation has been
completed for several configurations, but every trial has used a fixed payload + fixed
cell line + fixed delivery protocol applied uniformly across patients.

### 2.3 Computational Frameworks for Precision Oncology

TCGA (2008, 2013), CGGA, CPTAC, GLASS, and TCIA provide the public data substrate.
Published GBM survival models (DeepSurv and variants) typically achieve C-index
0.62–0.70 on internal CV and 0.60–0.65 on external test sets.

### 2.4 Research Gap and Novel Contribution

No existing framework connects tumor molecular profiling directly to the configuration
of a cell-delivered therapy. GlioStrata's contribution is the integration of tumor
subtyping, mechanistic simulation, Bayesian configuration optimization, and survival
prediction into a single open-source pipeline whose output is a ranked configuration
set per patient, validated against pre-registered thresholds.

---

## Chapter 3 — System Architecture

### 3.1 Design Principles

1. Each layer must be independently runnable and publishable.
2. Public-data-only by default — IRB-exempt secondary analysis.
3. Targets BU Shared Computing Cluster (SCC), Singularity, Snakemake.
4. Wet-lab output is a structured artifact (protocol document + predicted-outcome table).

### 3.2 Layered Architecture Overview

See `docs/architecture.md` for the diagram and Table 3-1 (layer responsibilities).

### 3.3 Data Flow

A representative trace: TCGA sample → L1 harmonize → L2 classify as Mesenchymal with
~62% MES-like cells, MGMT unmethylated → L3 feature vector → L4 simulate tumor growth
on patient-specific geometry → L5 predict 14-month survival under Stupp → L6 Bayesian
optimizer recommends iPSC-NSC carrying CRAd-survivin payload, ~2×10⁷ cells, 5 injection
sites, fractionated ×2 → predicted 22-month survival → L7 emits protocol document.

### 3.4 Infrastructure and Tooling

BU SCC for HPC, Singularity containers, Snakemake orchestration, DVC for data
versioning, MLflow for experiment tracking, GitHub Actions for CI.

---

## Chapter 4 — Data Resources and Ingestion (L0–L1)

### 4.1 Public Data Sources

See README.md "Data sources" table for the full catalog.

### 4.2 Harmonization Pipeline

Source-specific ETL modules → canonical sample-ID map → ComBat-seq (bulk) / harmony
(single-cell) for batch correction → Parquet (tabular) + AnnData h5ad (single-cell),
both DVC-versioned.

### 4.3 Quality Control Gates

Three gates: (1) sequencing depth + mapping rate + contamination; (2) tumor purity
(ABSOLUTE for WES, ESTIMATE for RNA-seq), flag < 40%; (3) molecular consistency (e.g.
IDH from WES == IDH from methylation).

### 4.4 Data Access and Ethical Compliance

All public, de-identified, secondary analysis → IRB exempt under 45 CFR 46.104(d)(4).
dbGaP application filed early (4–8 weeks approval window) for any controlled-access TCGA.

---

## Chapter 5 — Tumor Profiling and Stratification (L2–L3)

### 5.1 Verhaak Subtype Classification

Ensemble classifier (random forest + 1D-CNN) trained on TCGA-GBM labeled samples.
5-fold stratified CV on TCGA; final evaluation on frozen CGGA + GLASS.

### 5.2 Neftel Single-Cell State Decomposition

CIBERSORTx trained on the Neftel atlas as signature matrix. MuSiC2 as independent
deconvolver where matched single-cell data are available.

### 5.3 Tumor Microenvironment Deconvolution

CIBERSORTx with LM22 immune signature + MCP-counter for stromal/endothelial. MGMT
methylation from arrays (or inferred from transcriptome). IDH from WES or methylation.

### 5.4 Unified Patient Feature Store

Per-patient feature vector: genomic + transcriptomic + proteomic + radiomic + clinical
+ longitudinal. Parquet + DuckDB, DVC-versioned, with feature manifest.

---

## Chapter 6 — Disease Simulation and Outcome Modeling (L4–L5)

### 6.1 Tumor Growth Simulation

PhysiCell (agent-based). Patient-specific geometry from MRI segmentations via MONAI.
Calibrated against ≥10 PDX series spanning the four Verhaak subtypes.

### 6.2 NSC Migration and Payload Kinetics

Chemotaxis toward SDF-1 / VEGF / uPA / MCP-1 gradients. Validated against HB1.F3
luciferase in vivo imaging. Payload diffusion half-life constrained to 6–72 hours.

### 6.3 Survival Prediction Model

DeepSurv (Cox-PH neural network) trained on TCGA. Calibration via Brier score +
calibration plots at 12, 18, 24 months. Decision-curve analysis for clinical utility.

### 6.4 Calibration Strategy

L4 calibrated against biological ground truth (PDX + in vivo NSC tracking) independent
of patient outcomes. L5 calibrated against patient outcomes on training cohort +
externally validated on frozen test partitions. The two layers interact only at L6.

---

## Chapter 7 — Therapeutic Configuration Optimizer (L6–L7)

### 7.1 Bayesian Optimization Formulation

BoTorch (PyTorch-based). Multi-objective: simulated tumor kill + predicted survival
gain − toxicity penalty. PyMC for parameter uncertainty propagation.

### 7.2 Configuration Space

| Dimension | Domain | Constraint |
|-----------|--------|------------|
| Cell type | iPSC-NSC / MSC / HB1.F3-line | Immunogenicity profile per type |
| Payload | CRAd / prodrug-enzyme / CAR / combo | Compatible with cell type |
| Dose | 1×10⁶ to 5×10⁷ cells | Within delivery-vehicle capacity |
| Injection sites | 1–8 along resection cavity | Coordinates from patient MRI |
| Schedule | single / fractionated ×2 / ×3 | Within standard-of-care window |

### 7.3 Wet-Lab Handoff Protocol Generator

Per top-ranked configuration: cell line source, payload sequence, dose, injection
coordinates, schedule, controls, predicted outcomes, biosafety, ethical notes.
Auto-generated via pandoc + Jinja2 from JSON output. 15-item completeness checklist.

### 7.4 Recalibration Feedback Loop

PDO experimental results update L4 payload-kinetics priors and L5 survival model
via Bayesian updating. Recalibration procedure versioned in DVC.

---

## Chapter 8 — Test Plan and Expected Results

### 8.1 Defensibility Stack (10-Point)

1. Pre-registration on OSF before test data touched
2. Frozen test sets (DVC-locked)
3. Pre-specified thresholds (committed to OSF)
4. Bootstrap 95% CIs on every metric
5. Multiple-comparison correction (Bonferroni / BH-FDR)
6. External validation (frozen external test set)
7. Calibration, not just discrimination (Brier + Hosmer-Lemeshow)
8. Decision-curve analysis (net benefit)
9. Sensitivity analyses (imputation, batch correction, subtype subgroup)
10. Full reproducibility release (code + containers + weights + manifests)

### 8.2 Pre-Registered Computational Thresholds (T1–T5)

| Test | Metric | Threshold (95% CI lower bound) |
|------|--------|-------------------------------|
| T1 | Subtyping accuracy (CGGA external) | ≥ 0.88 |
| T1 | State deconvolution Pearson r | ≥ 0.75 |
| T2 | Survival C-index (CGGA external) | ≥ 0.65 |
| T2 | Beat best-published 2024 baseline by | ≥ 0.03 C-index |
| T2 | Brier score @ 18 months | ≤ 0.25 |
| T3 | Tumor-growth RMSE vs PDX | ≤ 15% of final volume |
| T3 | NSC-migration r vs in vivo | ≥ 0.65 |
| T4 | Optimizer beats random search | ≥ 25% lower regret (Wilcoxon p<0.01) |
| T4 | Personalization variance | ≥ 3 distinct top-5 configs across 50 patients |
| T5 | End-to-end success | ≥ 50% of held-out with ≥ 3-month RMST gain |
| T5 | Safety floor (predicted harm) | ≤ 10% (95% CI upper bound ≤ 15%) |

### 8.3 Wet-Lab Validation Tests (T6–T7)

T6: protocol generator completeness on top-10 predicted configs. T7: 3-arm PDO
experiment (Stupp / framework-recommended / fixed-payload active control) on ≥10 PDOs.

### 8.4 Resource Paper Deliverables

1. Open-source pipeline (Snakemake + Singularity + MIT license + GitHub CI)
2. Curated multi-omic GBM feature store
3. Trained model weights (subtype, deconvolver, survival)
4. Benchmark suite (T1, T3, T4 runnable by reviewers)
5. Protocol templates (T6 as standalone tool)
6. Documentation + tutorials

---

## Chapter 9 — Timeline, Milestones, and Risk Mitigation

### 9.1 Accelerated Schedule

Standard sequential estimate: 12 months. With agentic code generation on parallelizable
layers (L0–L3 + infrastructure): **8 months to paper submission**.

### 9.2 Parallelization Analysis

| Phase | Parallelizable? | Agents |
|-------|----------------|--------|
| L0/L1 ingestion (8 sources) | Yes | 8 |
| L2 subtyping (4 components) | Yes | 4 |
| L3 feature store | Partial | 1 + 4 |
| Infra (Snakemake/Singularity/DVC/CI/MLflow) | Yes | 5 |
| L4 PhysiCell simulation | No | 1 |
| L5 survival model | Partial | 1 |
| L6 Bayesian optimizer | No | 1 |
| L7 protocol generator | Partial | 1 |
| End-to-end integration | No | Human + 1 |

### 9.3 Eight-Month Milestone Map

| Month | Milestone | Tests Completed |
|-------|-----------|----------------|
| M1 | L0–L3 complete | — |
| M2 | L4 simulation + L5 trained | — |
| M3 | L6 optimizer + L7 scaffolding | — |
| M4 | End-to-end integration | — |
| M5 | Frozen-test-set validation | T1–T5 |
| M6 | Analysis + figures + sensitivity | T1–T5 |
| M7 | Resource release prep + paper draft | T1–T5 |
| M8 | Submit to NAR | T1–T5 |

### 9.4 Risk Register

| Risk | Likelihood | Fallback |
|------|-----------|----------|
| Survival C-index < 0.65 external | High | Report achieved C-index; reframe as benchmark |
| Wet-lab collaborator not identified by M3 | Medium | Computational resource paper stands alone |
| dbGaP access > 8 weeks | Medium | Use open-access TCGA subset only |
| End-to-end success threshold missed | Medium | Pre-registered negative result is publishable |
| Reviewer demands additional validation cohort | Medium | Add CPTAC-GBM as third external test |

### 9.5 Publication and Submission Strategy

Primary: **Nucleic Acids Research Database issue** (April submission window for
October publication). Fallback: **Bioinformatics**. Follow-up applied paper (if T7
wet-lab results return): Neuro-Oncology or Cell Reports Medicine.

---

## Chapter 10 — Conclusion and Immediate Next Steps

### 10.1 Summary

GlioStrata is a computationally complete Master's thesis protocol addressing a real
gap in stem-cell-delivered glioblastoma therapy using only public data and standard
Health Informatics skills. The eight-layer architecture, pre-registered thresholds,
and ten-point defensibility stack produce a project whose results will withstand peer
review whether or not every threshold is cleared.

### 10.2 Immediate Decisions Required

1. Draft and deposit the OSF pre-registration before any model touches frozen test partitions.
2. File the dbGaP application in Month 1 (4–8 week approval window).
3. Begin wet-lab collaborator outreach in Month 1 using M1 results as pitch material.
4. Confirm thesis committee structure (internal vs. external co-advisor permitted?).

### 10.3 First Two-Week Sprint

Three concrete artifacts: (1) OSF pre-registration JSON; (2) dbGaP application draft;
(3) L0 ingestion scaffolding (eight Snakemake rules, one per data source).

---

## Appendix A — Abbreviations

See `docs/architecture.md` and README.md for the abbreviation glossary.

## Appendix B — Software Stack

See `README.md` "Project structure" and `pyproject.toml` for the complete software stack.

## Appendix C — References

[1] Stupp R, et al. N Engl J Med. 2005;352(10):987–996.
[2] Stupp R, et al. JAMA. 2017;318(23):2306–2316.
[3] Verhaak RGW, et al. Cancer Cell. 2010;17(1):98–110.
[4] Neftel C, et al. Cell. 2019;178(4):835–849.e21.
[5] Darmanis S, et al. Cell Reports. 2017;21(5):1399–1410.
[6] Aboody KS, et al. PNAS. 2000;97(23):12846–12851.
[7] TCGA Research Network. Nature. 2008;455(7216):1061–1068.
[8] Brennan CW, et al. Cell. 2013;155(2):462–477.
[9] GLASS Consortium. Neuro-Oncology. 2019;21(7):805–812.
[10] CPTAC. Cancer Cell. 2021;39(10):1322–1335.
[11] Zhao Z, et al. Genomics Proteomics Bioinformatics. 2021;19(3):e210792.
[12] Katzman JL, et al. BMC Med Res Methodol. 2018;18:24.
[13] Ghaffarizadeh A, et al. PLoS Comput Biol. 2018;14(2):e1005991.
[14] Baltrušis T, et al. BoTorch. arXiv:1910.06403. 2019.
[15] Menze BH, et al. IEEE TMI. 2015;34(10):1993–2024.
[16] Clark K, et al. J Digital Imaging. 2013;26(6):1045–1057.
[17] Newman AM, et al. Nature Methods. 2015;12(5):453–457.
[18] Weller M, et al. J Clin Oncol. 2010;28(15):e241–e243.

---

*This markdown edition is current as of v0.1.0. The formatted .docx/.pdf editions
in `download/` are the canonical pre-registration documents.*
