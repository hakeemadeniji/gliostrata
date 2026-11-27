# GlioStrata: Statistical models and decision rules

Supplement to the OSF registration ("Statistical models" field). Authoritative text:
`docs/preregistration.md`, sections 2 to 6, at the commit tagged `prereg-deposit`.
Thresholds are frozen in `src/gliostrata/evaluation/thresholds.py` (SHA-256 `96bbcd49...74ab`).

## 1. The eleven confirmatory thresholds

| Test | Metric | Threshold | Decisive bound | Interval and unit |
| --- | --- | --- | --- | --- |
| T1 | Subtype accuracy (CGGA microarray, labelled) | >= 0.88 | lower | Clopper-Pearson, two-sided 95%; patient |
| T1 | Mean per-state Pearson r (100 mixtures) | >= 0.75 | lower | bootstrap over mixtures |
| T2 | C-index on CGGA | >= 0.65 | lower | bootstrap over patients |
| T2 | Improvement over the 2024 comparator | >= 0.03 | lower | paired bootstrap over patients |
| T2 | Brier score at 18 months | <= 0.25 | upper | bootstrap over patients |
| T3 | Growth RMSE, % of final volume | <= 15 | upper | bootstrap over curves |
| T3 | NSC migration r | >= 0.65 | lower | bootstrap over animals |
| T4 | Regret reduction, % | >= 25, and Wilcoxon p < 0.01 | lower | bootstrap over patients |
| T4 | Distinct configurations in the top five | >= 3 | lower | bootstrap over patients |
| T5 | Success fraction (gain >= 3 months) | >= 0.50 | lower | bootstrap over patients |
| T5 | Harm fraction (gain < 0) | <= 0.15, and point estimate <= 0.10 | upper | bootstrap over patients |

The bootstrap is a percentile bootstrap: 1000 resamples, fixed recorded seeds (42 unless
stated otherwise). Bounded metrics are clipped to their range.

## 2. Models per test

**T1.**
- *Subtype:* a nearest-centroid classifier with frozen TCGA centroids, applied to each
  sample alone. No fitting on CGGA, and no tuning.
- *Deconvolution:* the project's state deconvolution with default settings, on 100
  mixtures of 500 malignant cells each, drawn under Dirichlet(1,1,1,1) weights (seed
  131928).

**T2.**
- *Estimator:* Cox proportional hazards, Breslow partial likelihood, ridge 1e-6,
  covariates standardised internally. The optimiser runs to at most 200 iterations,
  tolerance 1e-9.
- *Models:*
  - committed: age, sex, KPS, treatment arm. NOT EVALUABLE on CGGA, which records no KPS.
  - transferable (variant 2, the evaluated stand-in): age, sex, treatment arm.
  - molecular (variant 3, secondary): variant 2 plus IDH and MGMT.
- *Coding:* indicator coding with the most frequent level as reference. No interactions.
- *Fitting:* fitted on TCGA complete cases and scored once on CGGA complete cases.
- *C-index and Brier rows:* evaluated before the deposit; FAIL (section 2.3).
- *Improvement row:*
  - our side: variant 2;
  - primary comparator: PMID 38174685 as published;
  - secondary, descriptive comparator: the same genes refitted;
  - metric: batch-stratified Harrell's C, paired bootstrap of the difference, on
    `cgga_external_v1`, after the deposit.

**T3.** The simulation layers are calibrated on the T3 datasets once those are obtained.
The metrics are those of the Indices supplement, with a bootstrap over curves or animals.

**T4.**
- *Methods:* GP expected improvement (Sobol initial design of 12, then 48 iterations)
  against random search (60 uniform draws), with the same 60-evaluation budget, on each
  of the first 50 patients by sorted id in `cgga_external_v1`.
- *Comparison:* paired per patient. The two-sided Wilcoxon signed-rank test is on the
  regret differences.
- *Ties:* patients with zero regret under both methods are dropped by the test. More than
  20 tied patients makes the row INCONCLUSIVE.

**T5.** Success and harm fractions over all 218 partition patients, from the frozen
simulator, response surface and survival model.

**T4 and T5 contingency.** These rows start only after their models are frozen by a dated
amendment made before any T4 or T5 read. They are NOT EVALUABLE if the models are not
frozen by 2027-06-30 (section 2.2 item 11).

## 3. Checks and controls (section 3.8; none decisive)

- **Proportional hazards:** a Grambsch-Therneau test on scaled Schoenfeld residuals, with
  rank-transformed event times, per term and globally, on every Cox fit. Already known: the
  TCGA fit of variant 2 violates proportional hazards in age (p = 0.0002).
- **Negative control:** each risk score's C-index against outcomes permuted within CGGA
  batch, 1000 permutations, seed 42. About 0.5 is expected.
- **Synthetic dry run:** the exact evaluation script runs on synthetic files in the real
  formats, before each read; its output is committed.
- **Positive controls:**
  - the deconvolution of mixtures with known composition;
  - reproducing the committed internal C-index (0.6735).

## 4. Decision rules (section 6)

- **PASS:** the decisive bound clears the threshold, and the side conditions hold.
- **INCONCLUSIVE:** the point estimate clears but the bound does not.
- **FAIL:** the point estimate misses.
- **NOT EVALUABLE:** the required data or comparator cannot be obtained.

Multiplicity is handled by Benjamini-Hochberg at FDR 0.05 across the p-values in the
family. CI-based rows are not re-corrected, and the family size (eleven) is always
reported. Interpretive precedence:
- an H1 pass with an H2 failure is reported as "a single configuration beats Stupp";
- a T5 pass built on a failed T2 is reported with that dependency.

## 5. Not part of the confirmatory analysis

The following are labelled exploratory:
- the section 3.6 sensitivity analyses, run only inside registered reads;
- the descriptive Brier scores at 12 and 24 months;
- the secondary comparator;
- any post hoc analysis.
