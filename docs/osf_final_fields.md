# OSF registration: text for the Sampling section

Paste-ready text for the four fields of the OSF "Sampling" section. It matches
`docs/preregistration.md` as of the commit that adds this file. Where the two
differ, the pre-registration is authoritative and this file is wrong.

No placeholders remain in the four fields. The items still marked "TO BE
COMPLETED BY AUTHOR" are all in the pre-registration's header (deposit DOI,
date, authors, repository URL, registered commit) and must be filled at deposit.

## Data collection procedures

```
No new data are generated from human participants. All data are existing, de-identified,
public cohort data downloaded in full and then filtered by fixed rules. There is no
recruitment and no sampling: every eligible patient in each downloaded release is included.

DEVELOPMENT COHORT: TCGA-GBM
Source and frame. All cases in project TCGA-GBM in the NCI Genomic Data Commons (GDC),
queried through the public REST API. Predominantly US-treated. 617 case records returned.
  - Clinical: demographics, diagnosis, treatments, follow-up, Karnofsky performance status.
  - Bulk RNA-seq: 284 samples (STAR read counts and TPM; 60,660 genes).
  - DNA methylation: 450 files, used to call MGMT promoter methylation.
  - Simple nucleotide variants: 464 files, used to call IDH1/IDH2 mutation status.
Eligibility. Cases are excluded, in this fixed order, if (1) the primary diagnosis names a
histology other than glioblastoma; a missing, "Unknown" or "Not Reported" diagnosis is
retained because enrolment in TCGA-GBM establishes the disease; (2) age, sex, survival time
and vital status are all missing; (3) overall survival is <= 0 months. Result: 594 patients.
Development split. Two sealed manifests exist. The one in force is tcga_gbm_clinical_v2:
425 train / 169 test by a hash of (seed, patient_id), seed 20251126. The earlier
tcga_gbm_clinical_v1 (583 patients) is retained unmodified as the historical record.

EXTERNAL COHORT: CGGA
Source and frame. All rows of the three open clinical files on cgga.org.cn (mRNAseq_693,
mRNAseq_325, mRNA_array_301) plus the STAR read-count matrices. Chinese cohort treated in
Chinese centres. 1,319 clinical rows pooled.
Eligibility. Rows are excluded if (1) the sample is recurrent or secondary rather than
primary; (2) histology is not exactly "GBM"; (3) grade is not WHO IV; (4) survival is <= 0.
Result: 333 rows, which are 305 unique patients because CGGA reuses identifiers across
batches. Each patient is analysed once.

ELIGIBILITY FOR THE T4/T5 EXTERNAL PARTITION
CGGA patients with an RNA-seq profile, a survival outcome and a usable identifier mapping.
Counted on 2026-11-26 from identifiers only (the expression matrix's column names, never
its values): 225 patients have an expression column, 218 of them also have a survival time
and vital status, so 218 patients are eligible (133 in mRNAseq_693, 85 in mRNAseq_325).
Every eligible patient is held out; the external cohort has no training split. The manifest
will be sealed after the single T2 read and before any T4 or T5 computation.

DATA STILL TO BE OBTAINED
- Microarray expression for the 108 CGGA patients who carry a reference Verhaak subtype
  label (batch mRNA_array_301). Only their clinical rows are held so far. The file will be
  stored unopened; the harmonisation rule is fixed in advance (below).
- Deconvolution ground truth: 100 pseudo-bulk mixtures simulated from the Neftel 2019
  single-cell atlas (GEO GSE131928; 28 glioblastoma tumours) with imposed Dirichlet random
  mixing weights.
- T3 inputs, not yet obtained. PDX growth curves: the candidates are the Mayo Clinic Brain
  Tumor PDX National Resource (not re-checked) and the Alcaniz 2023 panel
  (doi:10.3389/fonc.2023.1129627). A read-only check of the latter found 26 subcutaneous
  GBM models (15 also orthotopic), control arms and RNA-seq for 23 models, but it
  classifies three molecular subtypes rather than the four Verhaak subtypes, and whether
  its supplementary files carry per-time-point control volumes was not confirmed. In vivo
  NSC tracking: the candidate is Weber et al. 2022 (JoVE, doi:10.3791/63102). That protocol
  is not open access and lists no supplementary data, and its companion paper (Rust et al.
  2022, J Transl Med, doi:10.1186/s12967-022-03610-5) states that raw data are available on
  request and that the imaging could not resolve local cell migration. Neither source is
  recorded as obtained. A source qualifies as a PDX curve set only if it has longitudinal
  volumes, a control arm, and a subtype call or RNA-seq per model.

PREPARATION COMMON TO BOTH COHORTS
Fields are mapped to one schema with documented rules: survival converted to months with
one divisor; treatment arm derived from radiotherapy and chemotherapy records; vital status
coded 1 = death, 0 = censored. A missing molecular call stays missing and is never read as
wildtype. Every downloaded file is recorded with its GDC or CGGA file identifier and MD5
checksum.

Cross-platform harmonisation for T1 (fixed before the microarray matrix is downloaded).
Each sample's expression vector is z-scored across genes (for that sample alone, mean 0 and
variance 1 over its genes), and the RNA-seq-trained classifier is applied to the z-scored
vector. There is no shared reference and no cross-platform batch correction, and nothing in
the rule depends on the external data beyond the sample being scored. This is a project
choice, not the normalisation of Verhaak et al. (2010), who standardised each gene across
patients within each platform.

DURATION
TCGA-GBM and CGGA were downloaded, filtered and frozen on 2026-11-26. T3 inputs must be in
hand by 2027-01-15 and calibration complete by 2027-03-01. The external cohort has been
read once, for the T2 evaluation, and will not be read again for T4 or T5 until the
external partition is sealed.

KNOWN LIMITATION
The TCGA ingest queries the live GDC API rather than a pinned data release. The snapshot
used is defined by the recorded file identifiers and checksums, not by a release number.
The attached files give the full rules: docs/preregistration.md (sections 4 and 5),
docs/cgga_validation.md, and the two sealed development partition manifests.
```

## Sample size

```
The unit of analysis is the patient throughout, except where stated. Sample sizes are what
the public cohorts contain; none was chosen by power calculation.

T2 (survival prediction; ALREADY EVALUATED)
  Development (TCGA-GBM): 594 patients; 425 train and 169 test in the sealed development
  split. External (CGGA): 305 unique patients; 283 evaluable for the three-feature model
  and 241 for the five-feature model that also needs IDH and MGMT calls.

T1 (subtyping and deconvolution; NOT RUN). Two criteria with different units.
  Subtype accuracy: the 108 CGGA microarray patients with a reference Verhaak label (one
    row per patient). Cross-platform. Guard rails fixed in advance: the interval is exact
    binomial; if fewer than 30 labelled patients are available the criterion is NOT
    EVALUABLE. The 28 of these patients who also have RNA-seq are reported descriptively
    only and are not part of the criterion.
  State deconvolution: 100 simulated pseudo-bulk mixtures. The unit is the mixture, not a
    patient; the bootstrap resamples mixtures.

T3 (simulation fidelity; NOT RUN). Two criteria, each on its own unit.
  Tumour growth: PDX growth curves. One curve = one PDX model, control arm, mean volume per
    time point, at least 3 time points. At least 10 curves (protocol minimum), spanning the
    four Verhaak subtypes.
  NSC migration: in vivo tracking series; the resampling unit is the animal.
  Both counts are unknown until the data are in hand. Each is fixed on the date the dataset
  is obtained and recorded in the pre-registration before any calibration is run.

T4 (optimizer versus random search; NOT RUN)
  50 patients, paired: both algorithms run on the same patients with the same budget of 60
  objective evaluations each. Patients are the first 50 by sorted identifier within the
  sealed external partition, so selection does not depend on any outcome. That is 50 x 2 x
  60 = 6,000 simulator evaluations in total.

T5 (end-to-end internal consistency; NOT RUN)
  All 218 eligible held-out patients in the sealed external partition (counted on
  2026-11-26 from identifiers; fixed by the seal).

The T4 and T5 partition is sealed after one external read (T2) and before any T4 or T5
computation. No cohort is equalised across conditions; T4's pairing uses identical patients
by construction.
```

## Sample size rationale

```
There is no power analysis. Sizes are fixed by what the public cohorts contain and by four
rules fixed in advance.

1. Exhaustive inclusion. Every eligible patient in each release is used. No cohort was
   sized, trimmed or augmented after seeing a result.
2. A bound-based decision rule. PASS requires the decisive bound of a 95% CI to clear the
   threshold, not a point estimate. Small samples widen the interval; the registered
   verdict for an estimate on the right side with a bound that does not clear is
   INCONCLUSIVE.
3. One fixed external read per test. There is no sequential look, no sample-size
   re-estimation and no stopping for a favourable result.
4. Where a size would make a criterion undecidable, that is declared now, not after the
   fact, as follows.

T1 subtype accuracy (threshold: exact binomial lower bound >= 0.88, two-sided 95%
interval). On the 28 patients with both RNA-seq and a label, even a perfect 28 of 28
falls short (lower bound 0.8766). On 108 patients the criterion tolerates up to 6 errors
(102 of 108, lower bound 0.8830). An observed accuracy near 0.95 needs roughly 71 to 100
patients to clear, which is why the 108-patient microarray set is the scoring set, at the
cost of a cross-platform comparison. The floor of 30 is the smallest round number at
which a perfect score clears (30 of 30, lower bound 0.8843).

T1 deconvolution (threshold: r >= 0.75). The truth is imposed mixing weights, which is
optimistic: it tests recovery of known proportions, not agreement with a real tissue. 100
mixtures keep the correlation interval narrow enough for the threshold to be decidable. A
CIBERSORTx run, if performed, is a post-hoc sensitivity analysis and is labelled
exploratory.

T4 budget. 60 evaluations per patient per method, including a 12-point initial design, is
the existing code default. It was not tuned on any data, and it will not be calibrated on
the external cohort. Tie rule: a patient whose regret is zero under both methods is a tie
and is dropped by the signed-rank test. At most 20 tied patients (40% of the 50) are
tolerated; beyond that the p-value is still reported but the comparison is reported as
INCONCLUSIVE with the tie count disclosed.

T5. Its credibility is bounded by T2, which was not met (external C-index 0.5385, 95% CI
0.4994 to 0.5807). T5 is therefore reported as an internal consistency check on predicted
outcomes, not as a confirmatory external validation, and no T5 result may be described as a
survival benefit.
```

## Starting and stopping rules

```
No new data are collected, so these rules govern when each analysis starts and when it is
closed.

PILOT PHASE (synthetic data only; COMPLETE)
Pilot work ran on a seeded synthetic fixture and on TCGA-GBM development data only. It
ended when the full test suite passed on the frozen code. One exception is disclosed: the
T2 evaluation (external C-index 0.5385) was run before this registration.

START OF CONFIRMATORY ANALYSIS
A test starts only when all of the following hold: (1) this registration is deposited;
(2) its partition manifest is sealed and the verification step passes; (3) every input the
test needs has been named in advance. A test already run (T2) is reported as conducted and
is not re-run. Order for the external cohort: seal the partition first, then run T4 and
T5. For T1, the harmonisation rule is fixed before the microarray expression matrix is
examined, and it is fixed already.

STOPPING RULES
1. Single evaluation. Each test is evaluated once on its partition and then closed.
2. No interim looks. There is no interim analysis, no early stopping for efficacy or
   futility, and no extension of a cohort because a result was close.
3. Fixed cohorts. The TCGA-GBM and CGGA snapshots are closed as of 2026-11-26. A changed
   source release is frozen as a new, separately named partition; an existing one is never
   overwritten or topped up.
4. T3 dates. PDX growth-curve data and in vivo NSC tracking data must be in hand by
   2027-01-15; calibration must be complete by 2027-03-01. Each criterion is evaluated
   separately: a criterion whose data are not in hand by its date is closed as NOT
   EVALUABLE, with the reason stated, and the corresponding simulation layer is presented
   as uncalibrated. This is distinct from a FAIL.
5. Minimum n. If fewer than 30 CGGA patients carry a reference subtype label at analysis
   time, T1 subtype accuracy is closed as NOT EVALUABLE.
6. Comparator unavailable. T2's comparator criterion is closed as NOT EVALUABLE if the
   comparator cannot be re-run on identical data.
7. No further external reads without registration. CGGA has been read once. Any analysis
   beyond T4 and T5 requires a new pre-registration.
```

## Notes for the person pasting

- Attachments for "Data collection procedures - File upload": `docs/preregistration.md`
  (re-upload: it has changed since the first upload), `src/gliostrata/evaluation/thresholds.py`,
  both development partition manifests, and `docs/cgga_validation.md`.
- The two-sided bounds above differ from the decision memo's figures, which were
  one-sided. The pre-registration carries the two-sided bounds.
- The harmonisation rule is the memo's rule. The sentence attributing it to Verhaak et al.
  (2010) was removed because that paper standardised each gene across patients
  (PMC2818769, Methods).
