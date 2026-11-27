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
Development split. Two sealed development manifests exist. The one in force is tcga_gbm_clinical_v2:
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
was sealed on 2026-11-27, after the T2 evaluation (two reads) and before any T4 or T5
computation: cgga_external_v1, 218 patients, all-test, re-counted by the same rule at the
seal (133 in mRNAseq_693, 85 in mRNAseq_325). The seal read identifiers and the presence of
outcomes only and is recorded in the access ledger.

DATA STILL TO BE OBTAINED
- Microarray expression for the 108 CGGA patients who carry a reference Verhaak subtype
  label (batch mRNA_array_301). Only their clinical rows are held so far. The file will be
  stored unopened; the harmonisation rule is fixed in advance (below).
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

SINGLE-CELL REFERENCE: NEFTEL 2019 ATLAS (OBTAINED 2026-11-27; T1 DECONVOLUTION ONLY)
Source. The Smart-seq2 expression matrix of GEO GSE131928 (sample GSM3828672; 7,930 cells,
28 IDH-wild-type glioblastomas), and the per-cell malignant flag and six Neftel
meta-module scores from the Curated Cancer Cell Atlas (Tirosh lab) release of the same
data, which GEO does not carry. Both files are identified by size and SHA-256.
Use. 6,568 malignant cells with all six scores (26 tumours), each assigned its highest
state (NPC, OPC, AC, MES). They are the cells from which T1's 100 pseudo-bulk mixtures
were simulated, with imposed flat-Dirichlet weights. The mixture-generation rules were
written down and committed before the mixtures were made. The atlas is not a cohort and
holds no outcome; no patient is analysed.

PREPARATION COMMON TO BOTH COHORTS
Fields are mapped to one schema with documented rules: survival converted to months with
one divisor; treatment arm derived from radiotherapy and chemotherapy records; vital status
coded 1 = death, 0 = censored. A missing molecular call stays missing and is never read as
wildtype. Every downloaded file is recorded with its GDC or CGGA file identifier and MD5
checksum.

Cross-platform harmonisation for T1 (fixed before the microarray matrix is downloaded).
Each sample's expression vector is z-scored across genes (for that sample alone, mean 0 and
variance 1 over its genes), and the RNA-seq-trained classifier is applied to the z-scored
vector. The classifier is single-sample nearest centroid (Pearson correlation) over the
Verhaak 840-gene signature, against centroids frozen from 246 labelled TCGA RNA-seq
patients and identified by SHA-256; nothing is tuned and no other candidate is evaluated. There is no shared reference and no cross-platform batch correction, and nothing in
the rule depends on the external data beyond the sample being scored. This is a project
choice, not the normalisation of Verhaak et al. (2010), who standardised each gene across
patients within each platform.

DURATION
TCGA-GBM and CGGA were downloaded, filtered and frozen on 2026-11-26. T3 inputs must be in
hand by 2027-01-15 and calibration complete by 2027-03-01. The external cohort has been
read twice, both for T2 (the C-index row, then the Brier row), and will not be read again
for T4 or T5 until the external partition is sealed.

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

T1 (subtyping NOT RUN; deconvolution RUN BEFORE THE DEPOSIT). Two criteria with different
units.
  Subtype accuracy: the 108 CGGA microarray patients with a reference Verhaak label (one
    row per patient). Cross-platform. Guard rails fixed in advance: the interval is exact
    binomial; if fewer than 30 labelled patients are available the criterion is NOT
    EVALUABLE. The 28 of these patients who also have RNA-seq are reported descriptively
    only and are not part of the criterion.
  State deconvolution: 100 simulated pseudo-bulk mixtures of 500 malignant cells each,
    from the Neftel 2019 atlas. The unit is the mixture, not a patient; the bootstrap
    resamples mixtures. Evaluated on 2026-11-27, before the deposit: mean per-state
    Pearson r 0.9521 (95% CI 0.9430 to 0.9589) against 0.75.

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
  2026-11-26 from identifiers; fixed by the seal of 2026-11-27).

The T4 and T5 partition was sealed on 2026-11-27, after the T2 evaluation (two external
reads) and before any T4 or T5 computation. No cohort is equalised across conditions by construction.
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
3. One fixed external read per registered analysis. There is no sequential look, no sample-size
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

T5. Its credibility is bounded by T2, which was not met: external C-index 0.5385 (95% CI
0.4994 to 0.5807) and Brier score at 18 months 0.2617 (95% CI 0.2391 to 0.2829, against a
threshold of 0.25). T5 is therefore reported as an internal consistency check on predicted
outcomes, not as a confirmatory external validation, and no T5 result may be described as a
survival benefit.
```

## Starting and stopping rules

```
No new data are collected, so these rules govern when each analysis starts and when
it is closed. They apply to each of the eleven registered thresholds separately.

PILOT PHASE (COMPLETE)
Development and pilot work used a seeded synthetic cohort and the TCGA-GBM development
cohort, including the patients in its held-out test split. It ends at the code commit
recorded in this registration. Two things happened before the deposit and are
disclosed: (a) T2 was evaluated on CGGA in two reads (external C-index 0.5385; Brier
score at 18 months 0.2617), without a sealed external partition, because none existed
yet; (b) before the first read, CGGA had been characterised at schema level only (row
counts, missingness, value ranges), with no feature related to any outcome; (c) T1's
state-deconvolution criterion was evaluated on 100 mixtures simulated from the Neftel
2019 atlas (mean per-state r 0.9521, 95% CI 0.9430 to 0.9589, against 0.75), from rules
committed before the mixtures were made. It reads no cohort. These results changed no
threshold. The rules that were still open were fixed before the
deposit and are logged in the pre-registration's amendment log, which marks what was
written after the T2 read.

START OF CONFIRMATORY ANALYSIS
A registered test starts only when all of the following hold: (1) this registration
is deposited and the code commit is recorded; (2) every input the test needs has
been obtained and named in the pre-registration; (3) where the test is scored on a
cohort partition, that partition is sealed and the verification step reproduces it
from (seed, patient_id). This is required for T4 and T5; T1's scoring set and T3's
datasets are fixed by the registered rules instead.
Order for the external cohort (CGGA): seal the external partition first, then run T4
and T5. T1's subtype criterion is scored on the CGGA microarray batch, on the 108
labelled patients fixed in advance; its harmonisation rule is fixed and was fixed before
the matrix is downloaded. T1's deconvolution criterion has already been evaluated (pilot
phase, item c) and is not re-run. T3 is calibrated only after its datasets and their counts are entered in
the pre-registration. T2 is reported as already conducted and is not re-run.

STOPPING RULES
1. Single evaluation. Each threshold is evaluated once, on its sealed partition, and
   is then closed. Nothing is re-run after a null, failing or close result.
2. No interim looks. There is no interim analysis, no early stopping for efficacy or
   futility, no sample-size re-estimation, and no extension of a cohort because a
   result was close.
3. Fixed cohorts. The TCGA-GBM and CGGA snapshots are closed as of 2026-11-26. If a
   source release changes, the new snapshot is frozen as a new, separately named
   partition; an existing manifest is never overwritten or topped up.
4. Fixed budget. Computation stops at the pre-declared budget: 60 objective
   evaluations per patient per method for T4. It is not extended.
5. T3. PDX growth-curve data must be in hand by 2027-01-15 and calibration must be
   complete by 2027-03-01. If the data are not in hand by that date, the tumour-growth
   criterion is closed as NOT EVALUABLE, with the reason stated, and that layer is
   presented as uncalibrated. The NSC-migration criterion is closed the same way if
   in vivo tracking data cannot be obtained. NOT EVALUABLE is distinct from FAIL.
6. Minimum n. If fewer than 30 CGGA patients carry a reference subtype label when T1
   is run, T1 subtype accuracy is closed as NOT EVALUABLE.
7. Ties. If more than 20 of the 50 T4 patients are tied, the p-value is still
   reported but the regret criterion is INCONCLUSIVE whatever its value, and patients
   are not dropped and re-run to remove the tie.
8. Comparator. T2's improvement-over-baseline comparator is chosen by a registered
   eligibility rule and one literature search, both fixed before the search was run.
   If no model is eligible, or the comparator cannot be re-run on the identical frozen
   partition, the criterion is NOT EVALUABLE. It is never scored against a number
   quoted from a paper.
9. External reads. CGGA has been read twice, both for T2 (its C-index row, then its
   Brier row). The only further registered reads of CGGA are T1, T4 and T5. Any other
   analysis of CGGA requires a new pre-registration.
```

## Notes for the person pasting

- Attachments for "Data collection procedures - File upload": `docs/preregistration.md`
  (re-upload: it has changed since the first upload), `src/gliostrata/evaluation/thresholds.py`,
  both development partition manifests, the external manifest
  `data/partitions/cgga_external_v1.json`, `docs/cgga_validation.md`, and
  `docs/t1_deconvolution_spec.md`.
- The two-sided bounds above differ from the decision memo's figures, which were
  one-sided. The pre-registration carries the two-sided bounds.
- The harmonisation rule is the memo's rule. The sentence attributing it to Verhaak et al.
  (2010) was removed because that paper standardised each gene across patients within
  each platform (PMC2818769, Methods).
- "Starting and stopping rules" was revised on 2026-11-26 against the pre-registration.
  The earlier draft said any CGGA analysis beyond T4 and T5 needs a new pre-registration,
  which contradicted the registered T1 read of the CGGA microarray batch. It also said
  pilot work used development data only, without saying that the T2 models were fitted
  on all eligible TCGA patients, held-out test split included. Both are corrected above.
- T2's Brier-score row was evaluated on 2026-11-26 (a second read of CGGA) and is not
  met: 0.2617, 95% CI 0.2391 to 0.2829, against 0.25. Rule 9 lists it among the reads made.
- Revised on 2026-11-27 after T1's deconvolution criterion was evaluated before the
  deposit: the Neftel 2019 atlas moved from "data still to be obtained" to its own
  paragraph in "Data collection procedures"; "Sample size" gives the result; "Starting
  and stopping rules" discloses it as pilot-phase item (c). Upload
  `docs/t1_deconvolution_spec.md` with the other attachments. The pre-registration
  records the result in section 2.4.
- Revised on 2026-11-27 after cgga_external_v1 was sealed: "Data collection
  procedures" and "Sample size" say it is sealed. Upload `data/partitions/cgga_external_v1.json`
  with the two development manifests.
- Revised on 2026-11-27: the T1 classifier is named under "Cross-platform harmonisation
  for T1", and rule 8 now says the T2 comparator is chosen by a registered rule and one
  search (pre-registration section 2.2 item 2). If the search finds an eligible model,
  rule 9 must gain one read before the deposit; if it finds none, rule 9 is unchanged.
  Attach the search log.
