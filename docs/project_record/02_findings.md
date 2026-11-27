# 2. Findings

Every result the project has measured, as of 2026-11-26; findings 9 to 13 added on
2026-11-27 (stage 3b). Each section gives the design,
the numbers, what they do and do not show, and where the numbers come from.

**Reading labels.** **Met**: the point estimate and its decisive 95% bound both clear the
threshold. **FAIL**: the point estimate misses. **INCONCLUSIVE**: the point estimate clears
but the decisive bound does not. **NOT EVALUABLE**: cannot be computed on the data.
**Internal**: measured on TCGA-GBM by cross-validation, context only. **External**: measured
on CGGA.

## 2.1 Summary

| # | Finding | Status |
| --- | --- | --- |
| 1 | A Cox survival model fitted on TCGA-GBM shows no established discrimination on CGGA: external C-index 0.5385 (interval 0.4994 to 0.5807) against a threshold of 0.65. | T2 C-index: **not met** |
| 2 | The same model's predicted 18-month survival is less accurate than predicting the cohort average for everyone: Brier 0.2617 against a threshold of 0.25. | T2 Brier: **not met** |
| 3 | The registered four-feature model cannot be tested on CGGA, because CGGA records no performance status (KPS). | **Not evaluable** |
| 4 | Two descriptive steps sit between the registered model's internal C-index and its transfer: about 0.02 between the internal models, about 0.11 on moving to CGGA. Neither difference carries an interval. | Descriptive; no interval |
| 5 | IDH and MGMT status were called from raw TCGA data for 368 and 415 of 594 patients. | Descriptive; pipeline works |
| 6 | Excluding the 23 IDH-mutant patients moves the internal C-index by 0.0079, under the pre-specified limit of 0.01. | Sensitivity: **robust** |
| 7 | The external cohort is 305 patients, not the ">1,000" the protocol assumed; 218 are eligible for T4 and T5. | Cohort fact |
| 8 | T1's accuracy threshold cannot be met on the 28 labelled RNA-seq patients even with a perfect score; the 108-patient microarray set is needed. | Sample-size bound |
| 9 | Every missing KPS (152 of 594) is absent at the GDC source; none was lost in harmonisation. Resection extent is absent at source for all 594. | Data fact |
| 10 | The cohort-relative subtype calls depend on the scoring set: 2 to 18 of 283 calls change when the set is halved, and 5 to 7 of the 45 v2-test store patients change when they are scored alone. The single-sample calls change in none. | Descriptive; measures review item M6 |
| 11 | The single-sample subtype agrees with the published TCGA labels for 56 of 76 v2-test patients (0.737, exact 95% CI 0.623 to 0.831). | Descriptive; not validation |
| 12 | Every reported internal number had the 169 v2-test patients in its training folds. | Disclosure |
| 13 | The state and TME scores are cohort-relative: 6 of 50 patients change their highest state when scored alone. | Descriptive; not covered by D6 |

Two of eleven registered thresholds have been evaluated. Both were not met.

## 2.2 Finding 1: the survival model does not transfer to CGGA (T2, C-index)

**Design.** A Cox proportional-hazards model was fitted on TCGA-GBM and scored once on
CGGA with Harrell's C-index. Three variants were fixed, in the script and its tests, before
CGGA was read, and all three are reported.

- **Patients.** One row per patient (CGGA reuses identifiers: 333 rows are 305 patients).
  Variants 2 and 3, in both cohorts: complete cases on the variant's features and
  the outcome; nothing imputed. Variant 1 is fitted on TCGA only and is imputed
  there (see below the table).
- **Interval.** Percentile bootstrap over the evaluated CGGA patients, 1000 resamples, seed
  42, with the TCGA model held fixed.
- **Decision rule.** Decided by variant 2, as pre-specified. Met only if the point estimate
  and the lower bound of its 95% interval are both >= 0.65.

| Variant | Features | TCGA fit n | Internal C-index | CGGA n (deaths) | External C-index | 95% CI | Result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1. Committed baseline | age, sex, KPS, treatment arm | 594 passed, 592 modelled; KPS set to the median, 80, for 150 with an indicator (HR 1.57); treatment arm missing for 22 (level, HR 1.20) | 0.6735 | not evaluable | -- | -- | Not evaluable: no KPS in CGGA |
| 2. Transferable (decides) | age, sex, treatment arm | 570 | 0.6520 | 283 (239) | **0.5385** | 0.4994 to 0.5807 | **FAIL** |
| 3. Molecular (secondary) | age, sex, treatment arm, IDH, MGMT | 298 | 0.6606 | 241 (204) | 0.5915 | 0.5477 to 0.6341 | FAIL; not comparable with variant 2 |

Internal C-index: 5-fold stratified cross-validation, pooled out-of-fold.

Variant 1's TCGA fit imputes, although nothing is imputed on CGGA: 594 patients were
passed to the fit and 592 were modelled (2 have no vital status); KPS was set to its
TCGA median, 80, for 150 of them, with a missingness indicator that is prognostic in
its own right (HR 1.57); the treatment arm is missing for 22 and is coded as its own
level (HR 1.20); 491 deaths. Variants 2 and 3 impute nothing, internally or externally.
These counts are pinned by `tests/fixtures/tcga_fit_fingerprint.json`.

**What it shows.**
- The whole 95% interval of the headline model lies below 0.65. This is not a near miss
  that a larger sample would be expected to rescue.
- The interval contains 0.50. These data do not establish that the model discriminates in
  CGGA at all. They also do not establish that it cannot: the upper bound, 0.581, is
  consistent with weak discrimination.

**What it does not show.**
- That clinical features lack prognostic value in CGGA. A TCGA-fitted model failing to
  transfer is consistent with age, sex and treatment still mattering there.
- That no Cox model can reach 0.65 on CGGA. One specification, fitted on one cohort, was
  tested. Nothing else was tried, deliberately.
- Why the model failed to transfer. Candidate explanations are listed in section 2.4; none
  was tested, because testing them would mean further reads of the external cohort.

**Variant 3 is a cross-platform comparison.** TCGA's IDH and MGMT calls come from this
project's pipeline; CGGA's are text fields recorded by its contributing hospitals, with no
established assay or threshold. Variant 3 is also scored on a different, smaller set of
patients (241 against 283) and excludes the 20 no-adjuvant patients who were otherwise
evaluable, because the TCGA fit set contained no such patient with both calls. Its 0.053
advantage over variant 2 is therefore not evidence that the molecular markers help.

Source: `docs/cgga_validation.md` sections 1 to 5;
`data/results/cgga/cox_external_validation.json`; `scripts/validate_cox_cgga.py`.

## 2.3 Finding 2: predicted survival is less accurate than the cohort average (T2, Brier)

**Design.** The Brier score at 18 months is the average squared error between each
patient's predicted probability of surviving 18 months and what happened. It runs from 0
(perfect) to 1. Patients censored before 18 months are handled by
inverse-probability-of-censoring weights (Graf et al. 1999). The models and the patients
are the same as in section 2.2: the TCGA fits were confirmed identical to six decimals.
Met only if the point estimate and the upper bound of its 95% interval are both <= 0.25.

This was a **second read of CGGA**, made because this registered row had not been
evaluated. The pre-registration allows one read per registered analysis.

| Variant | CGGA n (deaths) | Alive at 18 months (observed) | Brier at 18 months | 95% CI | Null model | Result |
| --- | --- | --- | --- | --- | --- | --- |
| 2. Transferable (decides) | 283 (239) | 0.457 | **0.2617** | 0.2391 to 0.2829 | 0.2481 | **FAIL** |
| 3. Molecular (secondary) | 241 (204) | 0.478 | 0.2472 | 0.2226 to 0.2729 | 0.2495 | INCONCLUSIVE |

The **null model** predicts the cohort's own observed 18-month survival for every patient.
It is a descriptive reference, not a criterion.

Brier scores at the other horizons named in the pre-registration (descriptive only):

| Horizon | Variant 2 | 95% CI | Null model | Variant 3 | 95% CI | Null model |
| --- | --- | --- | --- | --- | --- | --- |
| 12 months | 0.2319 | 0.2112 to 0.2532 | 0.2277 | 0.2080 | 0.1842 to 0.2339 | 0.2146 |
| 18 months | 0.2617 | 0.2391 to 0.2829 | 0.2481 | 0.2472 | 0.2226 to 0.2729 | 0.2495 |
| 24 months | 0.2254 | 0.1972 to 0.2533 | 0.2202 | 0.2182 | 0.1874 to 0.2500 | 0.2229 |

**What it shows.**
- The headline model's point estimate is above 0.25 and so is the upper bound of its
  interval: FAIL under the registered rule.
- The model scores **worse than the null model at all three horizons** (0.2617 against
  0.2481 at 18 months). Predicting the same survival probability for everyone would have
  been closer to what happened.
- A Brier score above the null model's says that the predictions are worse than a
  constant. It does not say which of two things makes them so. With no resolution,
  the excess over the null model equals the squared offset (the mean predicted
  probability less the observed survival, squared) plus the variance of the
  predictions. Poor ranking, meaning dispersion that does not match outcomes, can
  worsen the score on its own; so can an offset. The two cannot be told apart
  without a calibration or resolution analysis, and neither was done. Which of them
  dominates here is not known, nor in which direction any offset points.

**A caution about the threshold itself.** With 45.7% of patients alive at 18 months, the
uninformative score is about 0.457 x 0.543 = 0.248, just under the registered 0.25. In this
cohort a model could have passed this row while discriminating not at all. The threshold
was fixed from a published range before the external cohort was read, is reported as
registered, and is not re-stated here.

**Not done, by decision (2026-11-27).** No calibration slope or intercept, calibration
plot or decision-curve analysis. Each would be a further read of CGGA. The decision not
to do them for this thesis is recorded in `docs/cgga_validation.md` section 12.5:
calibration is not a registered row, a further read would not change the T2 verdict, and
the offset-versus-dispersion question above is named as a limitation.

Source: `docs/cgga_validation.md` section 12;
`data/results/cgga/brier18_external.json`; `scripts/validate_brier_cgga.py`.

## 2.4 Finding 3 and 4: the registered model is untestable, and where the gap sits

**The registered model cannot be evaluated on the registered cohort.** T2 was registered
for a four-feature model that includes KPS. CGGA publishes no performance status in any of
its three clinical files. KPS was not imputed on CGGA: a TCGA median pasted onto Chinese
patients and scored as an observation would describe the imputation, not the model. This
is a feature-divergence problem, not a performance result.

**Where the gap sits. Descriptive only.** Three C-indices are on record. The two
differences between them carry no interval: neither was bootstrapped, so nothing here
says how much of either is sampling noise.

| Step | C-index | Difference (no interval) |
| --- | --- | --- |
| Committed baseline, TCGA (4 features) | 0.6735 | |
| Transferable model, TCGA (3 features) | 0.6520 | -0.0215 |
| Same 3-feature model, CGGA | 0.5385 | -0.1135 |

The first contrast is not "the effect of KPS". Variant 1 against variant 2 bundles three
differences at once: the KPS value itself; a missing-KPS indicator that is prognostic in
its own right (HR 1.57, for 150 of the 592 modelled patients); and a different fit set,
592 modelled patients with imputation against 570 complete cases. It is therefore not an
estimate of what KPS is worth. The second contrast compares one fitted model on two
cohorts. Read descriptively, the second gap (about 0.11) is about five times the first
(about 0.02). That ordering is consistent with the shortfall lying mainly in the move
between cohorts rather than in the loss of KPS; it is not a tested attribution.

Candidate explanations for the second step, **none tested**:
- A different population and care pathway. CGGA is a Chinese cohort with median age 52
  against TCGA-GBM's 59, and age is one of only three predictors.
- `treatment_arm` may not mean the same thing in the two cohorts, which derive it from
  different source fields. The mix differs: 80 of 333 CGGA rows are mono-therapy (24%)
  against 57 of 594 in TCGA (10%).

Source: `docs/cgga_validation.md` sections 2.1, 4.2 and 5.

## 2.5 Finding 5: molecular calls from raw TCGA data (descriptive)

Two calling pipelines were run against the live GDC API. These are data-preparation
results: they establish that the pipelines work and how many patients carry a call. They
are not external validation.

**IDH mutation status**, from open masked somatic-mutation files.

| | |
| --- | --- |
| Files / patients with a primary-tumour file | 464 files; 373 patients tested |
| IDH-mutant | 23 (6.2% of tested) |
| IDH-wildtype | 350 (93.8% of tested) |
| Cohort patients with a call | 368 of 594 |

The mutant rate is the single-digit figure expected for primary glioblastoma.

**MGMT promoter methylation**, from Illumina beta-value files. This is a two-probe
average score, **not** the clinically validated MGMT-STP27 assay.

| | |
| --- | --- |
| Files / callable patients | 450 files; 422 callable |
| Methylated at the default threshold of 0.35 | 128 (30.3% of callable) |
| Cohort patients with a call | 415 of 594 |

**The binary call depends heavily on the threshold.** Moving the cut from 0.35 to 0.20
changes the call for 48 patients (11% of the callable set). The two platforms also differ:
33.3% methylated on 27k against 24.3% on 450k, with no harmonisation between them. The
default is an a-priori working value, not a published cut-point; it was fixed before the
cohort was downloaded and left unchanged. The continuous score should be preferred
wherever a model allows it.

Of the 594 cohort patients, **314 have both an IDH and an MGMT call**. Of the 284 with
expression, **206 have both** and 12 have neither.

Source: `docs/mgmt_idh_calling.md` section 6.

## 2.6 Finding 6: the primary analysis is robust to excluding IDH-mutant patients

**Why it was done.** Under WHO 2021 an IDH-mutant grade IV tumour is astrocytoma, not
glioblastoma. All 23 IDH-mutant patients stay in the primary cohort, and the sensitivity
analysis repeats the primary analysis without them. The decision rule was recorded before
the IDH calls existed: flag the result if the C-index moves by more than 0.01.

| Cohort | Patients | Modelled | Deaths | Internal C-index | 95% CI |
| --- | --- | --- | --- | --- | --- |
| Full (primary cohort) | 594 | 592 | 491 | 0.6735 | 0.6478 to 0.6994 |
| Without the 23 IDH-mutant patients | 571 | 569 | 481 | 0.6656 | 0.6370 to 0.6913 |

The C-index moved by **0.0079**, below the 0.01 limit: **not flagged**. The spread across
cross-validation folds is several times the size of the change.

**Rerun with fold-wise imputation (2026-11-27).** The fits above imputed KPS with the
median of all modelled patients, taken before the folds were drawn; the code now takes
it from each fold's training rows. Rerun on TCGA only (`scripts/rerun_tcga_internal_cv.py`,
record `data/results/tcga_gbm/foldwise_imputation_rerun.json`): full cohort 0.673544
before and after, 95% CI 0.6478 to 0.6994; without the IDH-mutant patients 0.665650
before and after, 95% CI 0.6370 to 0.6913. Every training fold's median KPS is 80, so the
imputed values are unchanged and so are the predictions. The delta stays at 0.0079.

**Subtype distribution** of the 284 TCGA expression samples, and of the 268 left after
excluding the 16 IDH-mutant samples:

| Subtype | Full | % | Without IDH-mutant | % | Change |
| --- | --- | --- | --- | --- | --- |
| Classical | 88 | 30.99 | 86 | 32.09 | +1.10 pp |
| Mesenchymal | 82 | 28.87 | 77 | 28.73 | -0.14 pp |
| Proneural | 67 | 23.59 | 61 | 22.76 | -0.83 pp |
| Neural | 47 | 16.55 | 44 | 16.42 | -0.13 pp |

No subtype share moved by more than 1.10 percentage points. Two things sit behind that
stability: 14 of the 16 excluded samples were Proneural, and the classifier assigns labels
by ranking within a cohort, so 8 of the 268 retained samples received a different label in
the reduced run. The aggregate is stable; labels of samples near a subtype boundary are
less so.

**Not shown:** no paired test of the C-index difference was pre-specified, and none is
reported.

Source: `docs/mgmt_idh_calling.md` section 6.3;
`data/results/tcga_gbm/idh_sensitivity_cox.json` and `idh_sensitivity_subtype.json`.

## 2.7 Finding 7: the external cohort is much smaller than planned

The protocol describes CGGA as ">1,000 samples". That figure is all glioma grades.

| Step | n |
| --- | --- |
| Clinical rows pooled across three batches | 1,319 |
| Primary GBM, WHO grade IV, positive survival | 333 rows |
| Unique patients (identifiers are reused across batches) | **305** |
| With an RNA-seq expression profile | 225 |
| With expression and a complete survival outcome (**eligible for T4 and T5**) | **218** |

CGGA requires no registration: its files download directly. This corrects the protocol's
"open with registration".

Source: `docs/cgga_ingestion.md` sections 1 and 2.

## 2.8 Finding 8: a sample-size bound that shapes T1

T1's subtyping threshold requires the lower bound of an exact two-sided 95% interval to
reach 0.88.

| Labelled patients | Best possible result | Lower bound | Clears 0.88 |
| --- | --- | --- | --- |
| 28 (RNA-seq with a reference label) | 28 of 28 | 0.8766 | No |
| 30 | 30 of 30 | 0.8843 | Yes |
| 108 (microarray with a reference label) | tolerates up to 6 errors | 0.8830 at 102 of 108 | Yes |

Even a perfect score on the 28 RNA-seq patients cannot meet the threshold. The scoring set
is therefore the 108-patient microarray set, at the cost of a cross-platform comparison.
These bounds are two-sided; the decision memo carried one-sided figures, which are
superseded.

Source: `docs/preregistration.md` section 2.2 item 9; `docs/osf_final_fields.md`.

## 2.9 Findings 9 to 13: the TCGA inputs and the subtype classifier (stage 3b)

All on TCGA only; no CGGA value was read. Sources: `06_evidence_reports.md`, stage 3b;
`data/results/tcga_gbm/subtype_scoring_set_dependence.json`,
`subtype_label_comparison.json` and `feature_store_build_report.json`.

**Finding 9: what is missing is missing at source.** KPS is missing for 152 of the 594
cohort patients: 149 have no KPS in any GDC follow-up record and 3 carry only the string
"Unknown". Of the 592 modelled by variant 1, that is 150 (147 and 3); of the 283 with
expression, 64 (63 and 1). None had a numeric value that harmonisation lost, and a check
against the live GDC API (data release 46.0) found no numeric KPS in any follow-up of any
of the 152. Resection extent is absent at source: the raw clinical table has no resection
field, and no GDC residual-disease, residual-tumour or procedure field holds a value for
any of the 594. The one expression sample outside the cohort, TCGA-28-2510-01A, belongs to
a case whose GDC record is empty. Consequence: imputing KPS fills a value the record never
had (decision D4); the store records resection as not available, never as "biopsy"
(decision D3).

**Finding 10: the cohort-relative subtype depends on the scoring set.** The committed
signature classifier (`profiling/verhaak.py`) standardises each score across the samples
scored together. Measured against joint scoring of the 283 (seed 20261127, 20 random
halvings):

| Change to the scoring set | Calls compared | Cohort-relative: calls changed | Single-sample: calls changed |
| --- | --- | --- | --- |
| Random halves | 283 | mean 9.15 (2 to 18) | 0 |
| The 16 IDH-mutant samples removed | 267 | 10 | 0 |
| The 155 stage-3 store patients scored alone | 155 | 8 | 0 |
| The 45 v2-test store patients scored alone, against the 283-joint calls | 45 | 7 | 0 |
| The same 45, against the 155-joint calls | 45 | 5 | 0 |

The stage 3 store's subtype features therefore leaked the v2 test patients into each
other's values (review item M6). The per-patient subtype is now the single-sample call
(decision D6).

**Finding 11: the single-sample call against published labels.** Reference: Brennan et
al. 2013, whose labels are partly from the 2010 consensus clustering and partly from a
single-sample GSEA classifier. G-CIMP and unlabelled patients are excluded.

| Mode | Patients | Agree | Rate | Exact 95% CI |
| --- | --- | --- | --- | --- |
| Single-sample, centroids from the v2 training split, on the v2 test split | 76 | 56 | 0.737 | 0.623 to 0.831 |
| Cohort-relative, on the same 76 | 76 | 43 | 0.566 | 0.447 to 0.679 |

Of the 76 labels, 65 come from Brennan's ssGSEA classifier and 11 from the 2010
consensus clustering; the single-sample call agrees with 46 of the 65 and 10 of the 11.
This is a check on the cohort the centroids came from, mostly against another
classifier. It is not external validation and says nothing about T1, which is scored on
CGGA microarray patients and has not been run.

**Finding 12: the v2 test split was never held out.** The 169 v2-test patients are
reserved for internal held-out evaluation, but no script reads the split. Every reported
internal number had them in its training folds: the three T2 fits (including the TCGA
fits scored on CGGA), their internal C-indices, the IDH sensitivity runs and the
fold-wise rerun. The pre-registration discloses this for T2 only.

**Finding 13: the state and TME scores are cohort-relative too.** Their default
normalisation standardises each signature across the samples scored. Scoring the first 50
of the 283 alone moves a state fraction by up to 0.0744 and changes the highest state of
6 of the 50. The store's state and TME values were scored over the 283, v2-test patients
included. No decision covers this yet (`05_future_plan.md`, section 5.5).

## 2.10 What has not been shown

A reader should be able to tell a measured result from an absence of one.

- **No pre-registered test has been met.** T2's two evaluated rows failed. T1, T3, T4 and
  T5 have not been run, so nothing can be said about subtyping accuracy, simulation
  fidelity, the optimiser, or end-to-end benefit.
- **No survival benefit has been demonstrated, predicted or otherwise.** T5, the test that
  would measure a predicted benefit, has not been run, and its credibility is bounded by
  T2, which failed.
- **No calibration analysis exists.** No calibration slope, intercept, plot or
  decision-curve analysis has been produced for any model.
- **Nothing is validated against a published reference in L2.** The expression-based
  subtyping and deconvolution code has unit tests but has not been shown to reproduce
  published calls on a benchmark cohort.
- **The pipeline's own "profiling" output is not measurement.** It is derived from the
  clinical table.
- **No end-to-end run on the real cohort has been recorded.**
- **T2's improvement-over-baseline row is unevaluated**, because no comparator model has
  been named.
- **Both reads of CGGA happened without a sealed external partition**, because none
  existed. This departs from the pre-registration's stated procedure and is disclosed in
  its amendment log.

## 2.11 Reproducing the results

```bash
# T2 C-index row (synthetic-data tests, then the evaluation)
python -m pytest tests/test_validate_cox_cgga.py -q
python scripts/validate_cox_cgga.py

# T2 Brier row
python -m pytest tests/test_validate_brier_cgga.py -q
python scripts/validate_brier_cgga.py        # refuses to overwrite an existing result

# Molecular calls, join and the IDH sensitivity analysis
python scripts/download_tcga_mutations.py --out data/raw/tcga_gbm/mutations
python scripts/download_tcga_methylation.py --out data/raw/tcga_gbm/methylation
python scripts/join_molecular_calls.py
python scripts/idh_sensitivity_cox.py
python scripts/idh_sensitivity_subtype.py

# Stage 3b: the store, the centroid sets and the subtype checks (TCGA only)
gliostrata harmonize-clinical
gliostrata build-feature-store --mode real
python scripts/prepare_subtype_reference.py      # needs the data/external/ downloads, openpyxl and pyreadr
python scripts/build_subtype_centroids.py         # rebuilds configs/subtype_centroids/
python scripts/subtype_scoring_set_dependence.py  # finding 10
python scripts/compare_subtype_labels.py          # finding 11
```

- Both T2 scripts are deterministic (seed 42). Re-running one reproduces its numbers and
  is not a new analysis. Anything else done with CGGA is.
- The result files under `data/results/` and the raw data under `data/raw/` are not in
  git. A clean checkout reproduces them only after the download scripts have been run.
- The TCGA ingest queries the live GDC API, not a pinned release. The snapshot in use is
  defined by the recorded file identifiers and checksums.
