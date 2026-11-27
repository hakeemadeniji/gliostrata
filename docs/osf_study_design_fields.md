# OSF registration: Design, Variables and Analysis Plan fields

Paste-ready text for the OSF fields that `docs/osf_final_fields.md` (the Sampling section)
does not cover. Each block matches `docs/preregistration.md` at the commit tagged
`prereg-deposit`; where they differ, the pre-registration is authoritative. The three PDFs
named below are rendered from `docs/osf_supplements/`.

## Additional blinding during research or analysis

```
No analyst blinding. The project has one analyst, so no analysis is run blind to the design
or the outcomes. Unintended influence is limited by procedure instead (pre-registration
section 3.9):
- design choices are made without reference to external outcomes;
- every rule and every script for an external read is committed before that read;
- external test partitions are sealed, and CI re-derives each manifest from
  (seed, patient id);
- a read guard refuses any value read of the external cohort that is not registered,
  deposited and declared;
- every access, of any kind, is recorded in an append-only ledger;
- verdicts come from a registered threshold function, never from judgement;
- before each read, the exact evaluation script is run end to end on synthetic data, and
  its output is committed.
```

## Study design

```
Observational, retrospective secondary analysis of existing de-identified public cohorts,
with an in silico (computational) component. There is no recruitment, no randomisation of
patients and no counterbalancing.

Cohorts: TCGA-GBM is the development cohort (594 patients; sealed split tcga_gbm_clinical_v2,
425 train / 169 test). CGGA is the external cohort: the sealed partition cgga_external_v1
holds 218 patients, all held out, plus 108 labelled microarray patients for T1's subtype
accuracy. A single-cell atlas (Neftel 2019) and PDX/in vivo data serve T1 and T3.

Five tests:
- T1 profiling: subtype accuracy on labelled CGGA patients; state deconvolution on 100
  simulated pseudo-bulk mixtures.
- T2 survival prediction: a TCGA-fitted Cox model scored on CGGA. C-index and Brier are
  already evaluated. A paired within-patient comparison with a published 2024 model is
  still to come.
- T3 simulation fidelity: growth curves against PDX data; NSC migration against in vivo
  tracking.
- T4 optimiser versus random search: paired within-patient, 50 external patients, 60
  evaluations each.
- T5 end-to-end predicted benefit over the Stupp protocol: 218 external patients.

Between-subject throughout, except two paired (within-patient) comparisons: T2's model
comparison and T4's optimiser comparison. T4 and T5 run only after their simulator and
response surface are frozen by a dated amendment (pre-registration 2.2 item 11).
Pre-registration sections 1 to 6; Statistical models PDF.
```

## Manipulated variables

```
None in patients: there is no randomisation or assignment of treatment; every patient
variable is observed.

In silico only (T4/T5), the optimiser varies a therapy configuration. The space is
registered as coded (pre-registration 2.2 item 4):
- carrier cell type: iPSC-NSC, MSC or HB1_F3_line;
- payload: CRAd, prodrug_enzyme, CAR or combo;
- dose: 1 x 10^6 to 5 x 10^7 cells, continuous;
- number of injection sites: 1 to 8;
- schedule: single, fractionated_x2 or fractionated_x3.
HB1_F3_line x CAR and MSC x CRAd are excluded, which leaves 240 discrete combinations plus
dose.

T4 compares two search methods on each patient: GP expected improvement (12 Sobol points
plus 48 iterations) and random search (60 uniform draws).

T1's deconvolution mixtures use simulated mixing weights, drawn from Dirichlet(1,1,1,1),
seed 131928.
```

## Measured variables

```
Full list, with source, coding and missingness: "Measured variables" PDF.

Outcomes:
- overall survival in months (days / 30.4375);
- vital status (1 = death, 0 = censored).

Predictors and covariates:
- age (years);
- sex;
- Karnofsky performance status (TCGA only; CGGA records none);
- treatment arm (ChemoRadiation, MonoTherapy, NoAdjuvant), from radiotherapy and
  chemotherapy records;
- IDH1/2 hotspot mutation call;
- MGMT promoter methylation call;
- bulk tumour RNA-seq: STAR read counts and TPM in TCGA; STAR read counts and RSEM FPKM
  in CGGA.

Reference labels:
- published Verhaak subtype labels (TCGA, for the centroids; CGGA microarray, for T1's
  accuracy);
- single-cell malignant flags and Neftel meta-module scores (T1 deconvolution).

T3: PDX tumour volumes over time, and in vivo NSC tracking. Not yet obtained.

T4/T5: each patient's predicted survival under a configuration, from the frozen simulator.

Variables for subsetting and exclusion:
- histology, WHO grade, primary/recurrent status;
- the CGGA RNA-seq batch (mRNAseq_693 or mRNAseq_325), used as a stratum;
- membership of the development split;
- availability of an expression profile, an outcome or a reference label.
```

**File upload:** `GlioStrata_OSF_measured_variables.pdf`.

## Indices

```
Every derived variable and its formula: "Indices" PDF. The main ones:

- Treatment arm: ChemoRadiation if radiotherapy and chemotherapy were both given,
  MonoTherapy if one was, NoAdjuvant if neither, missing if either flag is missing.
- MGMT call: mean beta of cg12434587 and cg12981137, at least 0.35 = methylated.
- IDH call: an IDH1 R132 or IDH2 R172 hotspot variant.
- Expression scales: log2(TPM+1) (T1); log2(CPM+1) per sample (option B comparator);
  log2(FPKM+1) (published comparator).
- T1 subtype: Pearson correlation of each sample's per-sample z-scored log2(TPM+1) with
  four frozen TCGA centroids (831 Verhaak genes); the highest correlation gives the call.
- T2 risk scores:
  - our model: the Cox linear predictor of age, sex and treatment arm;
  - primary comparator (PMID 38174685): RS = 0.223340457*CPPED1 + 0.1517013*F3
    + 0.046311722*G0S2 - 0.219913858*MAPK1 + 0.192197058*MME + 0.016505054*MMP9
    + 0.414575274*MPO, on log2(FPKM+1).
- T2 metrics:
  - Harrell's C-index, batch-stratified for the improvement row;
  - IPCW Brier score at 18 months;
  - the paired difference C(ours) - C(comparator).
- T5:
  - predicted RMST to 24 months;
  - gain = RMST(configuration) - RMST(Stupp);
  - success = gain of at least 3 months; harm = gain below 0;
  - success fraction and harm fraction.
- T4:
  - per-patient simple regret;
  - % reduction = 100 x (1 - mean regret under GP-EI / mean regret under random search);
  - number of distinct configurations (cell type, payload, number of sites) among the
    top-five recommendations.
- T3:
  - RMSE as a percentage of each curve's own final volume, averaged over curves;
  - Pearson r of simulated against observed migration.
```

**File upload:** `GlioStrata_OSF_indices.pdf`.

## Statistical models

```
Full specification per test: "Statistical models" PDF. In brief:

T1:
- Subtype accuracy: proportion correct among the labelled CGGA microarray patients,
  exact two-sided Clopper-Pearson 95% interval.
- Deconvolution: mean of the four per-state Pearson correlations between estimated and
  imposed fractions over 100 mixtures, bootstrap over mixtures.

T2:
- Model: Cox proportional hazards (Breslow partial likelihood, ridge 1e-6), fitted on TCGA,
  with age (continuous), sex and treatment arm (indicator coding, most frequent level as
  reference).
- The registered four-feature model (with KPS) is NOT EVALUABLE on CGGA, which records no
  KPS.
- Improvement row: paired percentile bootstrap (1000, seed 42) of the batch-stratified
  C-index difference against PMID 38174685, applied as published. A refitted gene-set
  comparator is reported descriptively.

T3:
- RMSE as a percentage of each curve's final volume, bootstrap over curves.
- Pearson r of migration, bootstrap over animals.

T4:
- Percentage regret reduction of GP-EI against random search: bootstrap CI.
- Two-sided Wilcoxon signed-rank test on per-patient regret differences, p < 0.01, with a
  registered tie rule.
- Distinct-configuration count, bootstrap over patients.

T5:
- Success fraction and harm fraction over 218 patients, bootstrap over patients.

Checks (section 3.8), never decisive:
- Grambsch-Therneau proportional-hazards test (rank time transform) on every Cox fit;
- a within-batch permuted-outcome negative control (1000 permutations);
- a synthetic dry run before each read.
- Positive controls: deconvolution of mixtures with known composition; the reproduced
  internal C-index of the committed model.

No other models, interactions or contrasts; subgroup analyses are exploratory only.
```

**File upload:** `GlioStrata_OSF_statistical_models.pdf`.

## Transformations

```
Categorical variables are indicator-coded, with the most frequent TCGA level as the
reference: sex M; treatment arm ChemoRadiation. A CGGA level never seen in TCGA cannot be
scored and is excluded, never recoded.

Survival time is converted from days to months (/30.4375). Vital status is coded
1 = death, 0 = censored.

Expression is transformed per sample, never with cross-sample statistics:
- log2(TPM+1), then a per-sample z-score across genes (T1 subtype);
- log2(CPM+1), with the library size taken over the genes present in both cohorts
  (option B comparator);
- log2(FPKM+1) (published comparator, as published).

The Cox fit standardises covariates internally, and its coefficients are reported in
original units.

No batch correction is applied across cohorts for any confirmatory row. Instead, the
improvement row's C-index compares patients only within the same CGGA batch.
```

## Inference criteria

```
Every threshold is decided by the relevant bound of a two-sided 95% interval:
- percentile bootstrap, 1000 resamples, with the patient as the unit;
- exact Clopper-Pearson for T1 subtype accuracy;
- lower bound for higher-is-better metrics, upper bound for error and harm metrics.

Verdicts:
- PASS: the decisive bound clears the threshold, and every side condition holds (T4's
  Wilcoxon p < 0.01; T5's harm point estimate at most 0.10).
- INCONCLUSIVE: the point estimate clears but the bound does not.
- FAIL: the point estimate misses.
- NOT EVALUABLE: required data or the comparator cannot be obtained.

Tests are two-sided. Multiplicity: Benjamini-Hochberg at FDR 0.05 across the p-values in
the confirmatory family of eleven thresholds. CI-based thresholds are not re-corrected, but
the family size is always reported. Exploratory analyses are corrected separately and never
change a verdict. The thresholds are frozen by SHA-256 (96bbcd49...74ab).
```

## Data inclusion and exclusion

```
TCGA-GBM: all GDC cases, excluding, in order:
1. a non-glioblastoma diagnosis;
2. age, sex, survival time and vital status all missing;
3. overall survival <= 0.
That leaves 594 patients.

CGGA: excluding recurrent or secondary samples, histology other than GBM, grade other than
WHO IV, and survival <= 0. That leaves 305 unique patients (one row per patient, first in
stored order).

External partition: CGGA patients with an RNA-seq profile, a survival outcome and a usable
identifier, 218 in all, all held out. T1's subtype scoring set is the 108 labelled CGGA
microarray patients; below 30 it is NOT EVALUABLE.

Each analysis uses complete cases on its own inputs, and every exclusion is counted and
reported.

Outliers: none removed, trimmed or winsorised. Values are used as recorded at source;
known implausible values (KPS 0 for 3 TCGA patients) are disclosed (section 3.10).
```

## Missing data

```
Missing values are never imputed in the confirmatory analyses on CGGA. Each analysis uses
complete cases on its own inputs, and the numbers excluded are reported.

A missing molecular call is "no call", never wild type or unmethylated. A missing treatment
arm excludes the patient.

The only imputation anywhere is inside the TCGA fit of the committed four-feature model:
KPS median from the training fold, plus a missing indicator. That model is not scored on
CGGA.

Multiple imputation appears only as an exploratory sensitivity analysis inside registered
reads (section 3.6).
```

## Other planned analysis

```
Exploratory, labelled as such, never changing a verdict, and run only inside the reads
already registered (section 3.6):
- complete-case versus multiple imputation;
- with versus without cross-cohort batch correction;
- subgroups by Verhaak subtype, MGMT and IDH;
- seed sensitivity.

Also reported descriptively:
- the refitted best published gene set as a secondary comparator for T2's improvement
  row;
- Brier scores at 12 and 24 months;
- the section 3.8 checks (proportional hazards, negative control);
- per-state correlations for T1;
- any CIBERSORTx run (post hoc).

T6 and T7 (wet-lab) are outside this registration.
```

## Files to upload

The deposit bundle (`C:\Users\adeni\gliostrata_osf_bundle\`, built by
`scripts/osf_bundle.py` from the commit tagged `prereg-deposit`), uploaded under "Data
collection procedures - File upload", plus the three PDFs in their own fields:

| Field | File |
| --- | --- |
| Measured variables - File upload | `GlioStrata_OSF_measured_variables.pdf` |
| Indices - File upload | `GlioStrata_OSF_indices.pdf` |
| Statistical models - File upload | `GlioStrata_OSF_statistical_models.pdf` |
| Data collection procedures - File upload | everything in the bundle (see its `MANIFEST.txt`) |
