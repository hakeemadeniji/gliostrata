# T1, state-deconvolution row: mixture-generation specification

**Status: fixed before the run.** This document, `scripts/t1_deconvolution_mixtures.py`
and `tests/test_t1_deconvolution.py` were committed together, before the script was run
on the atlas. Every number below that the run depends on is a constant in that script.
The run is made once; a change after it is a deviation and is reported as one.

**Run before deposit.** The pre-registration is not deposited (`configs/registration.yaml`:
no OSF DOI). Until it is, T1 to T5 are exploratory (`docs/preregistration.md`, header),
and this row's result will be known before the deposit. The result is labelled "run
before deposit" wherever it is reported, and the deposit must disclose it.

## 1. What is registered, and what this document adds

Registered (`docs/preregistration.md`; `src/gliostrata/evaluation/thresholds.py`):

| Item | Registered text |
| --- | --- |
| Threshold | `T1_state_deconvolution_pearson_r`: >= 0.75 on the **lower** bound of the 95% CI (section 2) |
| Metric | The mean of the four per-state Pearson r values; per-state r reported alongside; the pooled r is descriptive only (section 2.2 item 1) |
| Truth | 100 pseudo-bulk mixtures simulated from the Neftel 2019 atlas (GEO GSE131928; 28 tumours) with imposed Dirichlet random mixing weights; the imposed weights are what the estimates are correlated with; optimistic, and reported as such (item 10) |
| Interval | Percentile bootstrap, 1000 resamples, fixed recorded seed; the resampling unit is the mixture; clipped to the natural range with the clipping reported (section 3.1; item 10) |
| Verdict | Section 6: PASS, INCONCLUSIVE, FAIL, computed by the registry |
| Other estimators | A CIBERSORTx run, if one is done, is post hoc and exploratory (item 10) |

Not registered, and fixed here: the annotation source, the eligible cells, the state rule,
the Dirichlet parameter, the cells per mixture, how a mixture is built, the
normalisation, the gene set, the estimator and its options, the seeds, the treatment of an
interval that excludes its point estimate, and the one secondary analysis.

## 2. Inputs

Recorded in `data/provenance/neftel2019_sources.json`; the script refuses a file whose
SHA-256 differs.

| Input | Source | Bytes | SHA-256 |
| --- | --- | --- | --- |
| Expression | GEO GSM3828672 (series GSE131928), `GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz` | 310,775,398 | `73be3319fb0b026dcae9dc33582ced9fef00a8e6470f63df3a7cbbf7daf960f5` |
| Cell annotation | Curated Cancer Cell Atlas (3CA), `Meta-data_Neftel2019_Brain.tar.gz`, member `SmartSeq2/Cells.csv` | 888,416 | `df01ad311ddf152e08c7ebb3c2c6c9fcf22b119f5a4a9744aa0e81a37c7a7b80` |

- **Why two sources.** The GEO series has no per-cell malignant flag or state; 3CA, the
  Tirosh lab's curated release of the same data, has both a malignant flag and the six
  Neftel meta-module scores for each Smart-seq2 cell. They are used as distributed by
  3CA; whether 3CA recomputed them or copied the authors' values is not established. The Broad Single Cell Portal copy
  (SCP393) refused an anonymous download. Every 3CA cell name matches a GEO column.
- **Scale.** The GEO values are log2(TPM/10 + 1). They are converted back with
  TPM = 10 x (2^x - 1); a cell then sums to 1,000,000. Checked on the first 60 cells
  before this document was written; the run records the minimum and maximum sum over
  every cell it uses.
- **Platform.** Smart-seq2 only (7,930 cells, 28 tumours). The 10x data carry no module
  scores in 3CA and are not used.

What was looked at before this document was written: file sizes and checksums; the
matrix's header, its gene column and the scale of its first 60 cells; the annotation's
columns and value counts; and the cell counts below. No mixture was made and no
estimator was run on atlas data.

## 3. Cells and their state

- **Eligible:** `malignant == "yes"` and all six meta-module scores present (`MESlike1`,
  `MESlike2`, `AClike`, `OPClike`, `NPClike1`, `NPClike2`). 6,568 of 7,930 cells, from
  26 tumours.
- **State:** the highest of four scores, NPC = max(NPClike1, NPClike2), OPC = OPClike,
  AC = AClike, MES = max(MESlike1, MESlike2). This is Neftel's merging of the two NPC and
  the two MES sub-modules; taking the highest of the four is the quadrant of the paper's
  two-dimensional state plot (it compares max(OPC, NPC) with max(AC, MES), then the two
  on the winning side). A tie for the highest is not assigned and is counted: there are
  none.
- **Counts:** NPC 1,985; OPC 1,047; AC 1,928; MES 1,608.
- Cells are pooled across tumours. A mixture is not a single tumour.

## 4. Mixtures

| Item | Value |
| --- | --- |
| Number | 100 (registered) |
| Weights | w ~ Dirichlet(1, 1, 1, 1), flat over the four states, state order NPC, OPC, AC, MES |
| Cells per mixture | 500 malignant cells |
| Cell counts | n ~ Multinomial(500, w) |
| Cells | n_s cells drawn uniformly without replacement from the pool of state s |
| Mixture | the mean of the 500 cells' linear TPM vectors (so it sums to 1,000,000) |
| Truth | the imposed weights w (registered). The realised fractions n/500 are reported as a descriptive secondary |
| Seed | `numpy.random.default_rng(131928)`, one generator: the 100 weight vectors first, then for each mixture in order its counts and its draws, state by state |
| Genes | all 23,686 genes of the matrix, by symbol; nothing is filtered or renamed |
| Normalisation | none beyond linear TPM; no log, no scaling, no batch step |

The seed 131928 is new to this analysis (the GEO series number), not a registered seed.

## 5. Estimator and scoring

- **Estimator (primary).** `gliostrata.profiling.deconvolve.deconvolve_states_expression`
  with its defaults (`score_method="ssgsea"`, `fraction_method="softmax_z"`,
  `temperature=1.0`), applied once to the 100 mixtures together, symbol-indexed, with no
  annotation. This is the method the feature store uses. All 59 of its marker genes are in
  the atlas (NPC 15, OPC 14, AC 16, MES 14). It is a marker-enrichment score, not
  CIBERSORTx, which the protocol names; its default normalisation standardises each state
  across the samples scored, here the 100 mixtures.
- **Metric.** For each state, the Pearson r between the estimated fraction and the
  imposed weight across the 100 mixtures; the primary metric is the mean of the four.
- **Interval.** 1000 resamples of the 100 mixtures with replacement,
  `numpy.random.default_rng(42)`; the estimates are held fixed (the estimator is not
  re-run in a resample, as the T2 bootstraps held the model fixed). For each resample the
  four r and their mean are recomputed; the interval is the 2.5th and 97.5th percentiles,
  clipped to [-1, 1] with the clipping reported. A resample in which any per-state r is
  undefined is dropped and counted; if more than 10 of the 1000 are dropped, no interval
  is computed and the row is NOT EVALUABLE, with the count as the reason.
- **Undefined on the full set.** If a per-state r is undefined on the 100 mixtures (an
  estimate or a weight constant across them), the row is NOT EVALUABLE, with the reason.
- **Verdict.** `gliostrata.evaluation.thresholds.evaluate("T1_state_deconvolution_pearson_r",
  point, lower, upper)`: PASS if the lower bound is >= 0.75; INCONCLUSIVE if the point
  estimate is >= 0.75 and the lower bound is not; FAIL if the point estimate is below
  0.75. If the point estimate lies outside its percentile interval, the interval is
  widened to include it before the registry is called, and that is reported; this can
  only lower the decisive bound.
- **Reported alongside, descriptive:** each per-state r with its interval; the pooled r
  over all 400 pairs; the per-state r against the realised fractions; the RMSE per state.
- **One secondary analysis, declared here and descriptive:** the same estimator with
  `fraction_method="proportional"` (cohort-independent), on the same mixtures, scored
  the same way. It cannot change the verdict. It is declared now because the store's
  state scores are cohort-relative and a decision on that is open
  (`docs/project_record/05_future_plan.md` section 5.5); nothing else is run.

## 6. Why the result is optimistic

Registered as an optimistic test (item 10). Specifically:

1. The truth is simulated: the mixtures are made from the same cells whose states define
   the truth.
2. The mixtures hold malignant cells only. Real tumours carry immune, stromal, vascular
   and normal-brain cells, which the estimator's state panels do not see here.
3. Circularity. A cell's state comes from Neftel's meta-module scores, and the
   estimator's marker panels are anchor genes of the same meta-modules. Agreement is
   partly built in.
4. Smart-seq2 single cells averaged are not bulk RNA-seq of tissue.
5. Mixtures pool cells across tumours, so tumour-specific programmes are averaged out.

The result says how well the estimator recovers known state proportions under these
conditions. It does not say how well it estimates state composition in a TCGA or CGGA
tumour.

## 7. Outputs and the single run

- `data/results/t1/neftel_mixture_deconvolution.json`: inputs with SHA-256, counts,
  constants, the primary result, the verdict, the secondary, software versions.
- `data/results/t1/neftel_mixture_deconvolution_mixtures.parquet`: per mixture, the
  imposed weights, the realised fractions and both estimators' outputs (no cell or patient
  identifiers).
- Both are git-ignored; their SHA-256 go into the project record.
- The script refuses to overwrite the result file. It is run once.
- **An aborted run.** A run that stops with an error before the result file is written,
  having printed no metric, may be repeated unchanged once the cause is fixed. The abort,
  its error and the fix are recorded with the result. A fix that changes any rule of this
  document is a deviation and is reported as one. A re-run, a different
  seed, a different cell count or another estimator after the result is seen is a new,
  exploratory analysis and is reported as such.
- Nothing here reads TCGA or CGGA; a path containing "cgga" is refused.

## 8. Reproducing

```bash
python -m pytest tests/test_t1_deconvolution.py -q
python scripts/t1_deconvolution_mixtures.py --dry-run   # checksums and counts only
python scripts/t1_deconvolution_mixtures.py             # the run (refuses to overwrite)
```
