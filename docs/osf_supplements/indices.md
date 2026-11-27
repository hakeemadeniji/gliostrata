# GlioStrata: Indices and derived variables

Supplement to the OSF registration ("Indices" field). Formulas as registered in
`docs/preregistration.md` (sections 2.2 and 3) at the commit tagged `prereg-deposit`.

## 1. Clinical derivations

- **Overall survival (months)** = days / 30.4375. The event is death (1); the alternative
  is censoring (0).
- **Treatment arm**, from the radiotherapy and chemotherapy flags:
  - ChemoRadiation: both given;
  - MonoTherapy: exactly one given;
  - NoAdjuvant: neither given;
  - missing: either flag missing. An unknown arm is never guessed.
- **IDH call** = 1 if the tumour carries an IDH1 R132 or IDH2 R172 hotspot variant, 0 if
  it was sequenced without one, and missing if it was not tested.
- **MGMT call** = 1 if the mean beta of probes cg12434587 and cg12981137 is at least 0.35.
  This is a transparent two-probe score, not MGMT-STP27. A missing call stays missing.

## 2. Expression scales

All are computed per sample, from that sample alone:

| Scale | Formula | Used for |
| --- | --- | --- |
| log2(TPM + 1), then per-sample z-score | z = (x - mean_g x) / sd_g x over the sample's genes | T1 subtype |
| log2(CPM + 1) | CPM_g = 10^6 count_g / sum over the shared gene universe (33,212 symbols) | T2 secondary comparator |
| log2(FPKM + 1) | FPKM as distributed by CGGA (RSEM) | T2 primary comparator |

## 3. T1 profiling

- **Subtype call.** Pearson correlation between the sample's z-scored vector and each of
  four frozen centroids (`verhaak840_tcga_all`; 831 genes; means over 246 labelled TCGA
  patients). The call is the class with the highest correlation; ties go to Classical,
  then Mesenchymal, Proneural, Neural.
- **Subtype accuracy** = correct calls / labelled patients.
- **State deconvolution.** Estimated NPC, OPC, AC and MES fractions from
  `deconvolve_states_expression` (defaults). The metric is the mean over the four states
  of Pearson r(estimated, imposed), across 100 mixtures.

## 4. T2 survival

- **Our risk score** (variant 2) = beta_age x age + beta_F x [sex = F]
  + beta_Mono x [arm = MonoTherapy] + beta_NoAdj x [arm = NoAdjuvant]: the Cox linear
  predictor, fitted on TCGA.
- **Primary comparator** (PMID 38174685, applied as published, on log2(FPKM + 1)):
  RS = 0.223340457 CPPED1 + 0.1517013 F3 + 0.046311722 G0S2 - 0.219913858 MAPK1
  + 0.192197058 MME + 0.016505054 MMP9 + 0.414575274 MPO
- **Secondary comparator**: the same seven genes on log2(CPM + 1), refitted with the
  project's Cox estimator on 282 TCGA patients. The coefficients are frozen in
  `configs/t2_comparator/chosen.json`.
- **Harrell's C-index**: concordant comparable pairs / comparable pairs, ties counting 0.5.
  A pair (i, j) is comparable when t_i < t_j and i died. In the **batch-stratified**
  version, only pairs within the same CGGA batch count, summed over batches.
- **Improvement** = C(our score) - C(comparator), on the same patients.
- **Brier score at 18 months**: inverse-probability-of-censoring weighted squared error of
  the predicted survival probability at 18 months.

## 5. T3 simulation

- **Growth error** = RMSE between simulated and observed volumes over each curve's time
  points, divided by that curve's observed final volume (x 100), then averaged over
  curves.
- **Migration agreement** = Pearson r between simulated and observed migration.

## 6. T4 and T5 (pre-registration 2.2 items 4, 5, 11)

- **Predicted RMST(tau = 24 months)** = the area under the patient's predicted survival
  curve from 0 to 24 months.
- **Gain** = RMST(chosen configuration) - RMST(Stupp). Stupp is the same survival model
  with the treatment arm set to radiotherapy plus temozolomide and no stem-cell
  configuration.
- **T5 success fraction** = the share of the 218 patients with a gain of at least 3
  months. **Harm fraction** = the share with a gain below 0.
- **T4 regret** for patient i under method m = the best objective value found for i by
  either method, minus the best value found by m (months of predicted RMST).
- **T4 percentage reduction** = 100 x (1 - mean regret under GP-EI / mean regret under
  random search).
- **Distinct configurations** = the number of distinct (cell type, payload, number of
  injection sites) among the top five recommendations across the 50 patients.
