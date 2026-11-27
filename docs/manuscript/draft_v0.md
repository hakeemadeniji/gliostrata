# GlioStrata manuscript: working draft v0

Status: working draft, not for circulation. Two sections that can be written now.
Every number is traced to a repository file in the list at the end (section 3). Places
that need the lead's decision are marked [LEAD: ...].

Standing caveat. At the time of writing `docs/preregistration.md` is a draft with no OSF
DOI and no deposit date (its header). By the project's own rule, T1 to T5 are exploratory
until the deposit. Everything below that reports a result was obtained before the deposit.
[LEAD: update this paragraph with the DOI and deposit date once recorded, and say which
results, if any, came after it. As of this draft none did.]

---

## Part 1. A reproducibility audit of 2024 glioblastoma survival models

[LEAD: the title is the lead's working title. Section 1.4 explains why "reproducibility"
is narrower here than the title suggests. Consider "Which 2024 glioblastoma survival
models can be re-run on an external cohort without refitting? A registered search".]

### Abstract (draft)

We searched PubMed once for 2024 glioblastoma survival models (549 records) and asked
which could be re-run, as published, on a held-out external cohort. Records were screened
against six pre-registered criteria in stages: metadata, title and abstract, full text,
and a second independent review of the closest calls. One model met every criterion
available for assessment. The most common first failure was criterion E3, a complete
formula that states its expression scale: 85 of 549 records, 80 judged on the full text
and 5 not read. Of the 88 papers that reached the full-text stage, 48 used the Chinese
Glioma Genome Atlas (CGGA) in some role, and 21 used it for training, gene selection or a
cut-off. Screening was AI-assisted, used human-defined criteria, and searched one
database. [LEAD: add the abstract's last sentence once the improvement row is read.]

### 1.1 Introduction

Many gene signatures for glioblastoma survival are published each year. Most report a
C-index or an AUC on a validation cohort chosen by the authors. A reader who wants to know
whether a signature works on a new cohort needs the model itself: every coefficient, and
the scale of expression the coefficients were fitted on.

This paper comes out of a pre-registered project, GlioStrata. One of its registered tests
(T2) asks whether the project's survival model beats a "best-published 2024 baseline" by a
C-index of 0.03 or more, with the lower bound of the 95% interval deciding
(`docs/preregistration.md` section 2). The protocol named no baseline. A number quoted
from a paper, computed on another cohort, has no sampling distribution. The registration
therefore requires the comparator to be re-run on the identical frozen partition
(section 2.2 item 2). That in turn needs a published model that can be re-run.

We registered a rule and one search for choosing the comparator, before the search was
run, and committed the log. This paper reports that search. Its result is a count of how
many 2024 models could be re-run as published, and why the rest could not.

The cohort we need to score on is made of two CGGA RNA-seq batches, mRNAseq_325 and
mRNAseq_693. These batches are widely used to "externally validate" signatures. A model
that was trained on, or had its genes or cut-off chosen on, either batch cannot be
re-run on them as a held-out test. We therefore also report how often these batches were
used, and in which role.

### 1.2 Methods

**Registration.** The rule and the search were registered in `docs/preregistration.md`
section 2.2 item 2 and committed before the search was run (commit `aacedb1`, per
`docs/t2_comparator_search/README.md`). The rule was later amended, the same day and
before the deposit, after the search showed that almost no paper stated its expression
scale; the amendment and its reasons are in section 8 of the registration and are
described in section 1.4 below. This paper reports the search under the rule as first
registered.

**Search.** One PubMed search, run once at 2026-11-27T17:48:54Z, with the query:

`glioblastoma[tiab] AND (prognos*[tiab] OR survival[tiab]) AND (signature[tiab] OR
model[tiab] OR nomogram[tiab] OR "risk score"[tiab]) AND 2024[dp]`

PubMed applied the date limit as 2024/01/01 to 2024/12/31 [Date - Publication]. It
returned 549 records. No other database and no second query were used
(`docs/t2_comparator_search/search.json`).

**Criteria.** A model was eligible only if it met all of E1 to E5 (E6 only ranks eligible
models and did not arise).

- E1: peer-reviewed, first published in 2024. The earlier of the electronic and print
  dates must fall in 2024. Preprints and retracted papers fail.
- E2: predicts overall survival in adult WHO grade 4 glioblastoma. Pan-glioma models that
  include grade 4 glioblastoma pass at the abstract stage.
- E3: can be re-run without refitting. Either a complete published formula (every term,
  its coefficient, and the expression scale it expects) or released code with fitted
  weights.
- E4: every input is available for the target cohort: bulk RNA-seq plus age, sex, IDH,
  MGMT, 1p/19q, radiotherapy and chemotherapy (yes or no). Imaging, histology, DNA
  methylation arrays, performance status and extent of resection fail.
- E5: no CGGA mRNAseq_325 or mRNAseq_693 patient was used for training, gene or feature
  selection, coefficient fitting or cut-off choice. Use for validation only passes.

The first criterion a record failed, in the order E1, E2, E4, then E3 and E5 on the full
text, is its "first failing" criterion. Doubt was never resolved in favour of eligibility.
Sources: `docs/preregistration.md` section 2.2 item 2; `docs/t2_comparator_search/protocols/`.

**Screening, by stage.**

1. *Date and type (script).* E1 was read from PubMed metadata by script, for all 549
   records.
2. *Title and abstract (dual screening).* Two independent screeners judged each of the
   480 remaining records against E1 (non-articles), E2 and E4. Four disagreements were
   adjudicated from the abstract.
3. *Full text.* One reviewer per paper judged E3 and E5, and any E2 or E4 failure that
   only the full text showed. The reviewer recorded, for each judgement, a quotation from
   the paper. Every quotation was checked against the paper's text by script, and a
   quotation that was not found counted against the record. Full texts were the open-access
   Europe PMC texts (71 papers), two publisher pages, and 10 PDFs the author obtained
   through a university library.
4. *Second independent review of close calls.* A second reviewer, who had not seen the
   first reviewer's verdicts, re-judged E2 to E5 on 25 papers: the two that passed, six
   with coefficients in a supplement (the supplements were downloaded and read), 16 with
   complete coefficients but a doubtful scale, and one doubtful E5 call. Three calls
   were adjudicated.

**Who screened.** The screening was AI-assisted. The criteria, the written protocols,
the adjudications and the final rulings were the project author's. The screeners and
reviewers were AI agents working from those protocols. [LEAD: confirm this wording, and
supply the model name and version of the reviewers, the dates, and whether the two
independent reviewers were the same model. The repository record says "screeners" and
"reviewers" and does not name them.]

**Counts.** Every count in this paper was computed by a script that reads only
`docs/t2_comparator_search/records.jsonl` (one line per PMID). The script is in the
Appendix. Two groupings are not fields in the file and were coded by hand from its
free-text fields: the role CGGA played (from `e5_detail`) and the sub-reason for an E3
failure (from `reason`). The hand-coded PMID lists are in the script, so a reader can
check or change them. [LEAD: both codings are one reader's, and should be checked
before submission.]

### 1.3 Results

**Flow.** The search returned 549 records.

| Stage | Records in | Failed | Passed on |
| --- | --- | --- | --- |
| Date and type (E1) | 549 | 69 | 480 |
| Title and abstract (E1, E2, E4) | 480 | 392 (1 E1, 330 E2, 61 E4) | 88 |
| Full text: read | 83 | 82 | 1 eligible |
| Full text: not read | 5 | 5 | 0 |

Of the 88 papers at the full-text stage, 83 were read and 5 were not. Three of the five
could not be obtained, even through the library (PMIDs 38959943, 38318849, 37539935). The
other two (39476302, 38662286) were judged on their abstracts, which stated that CGGA
was used for training. All five are recorded as "full text not obtained; E3 not
established".

**One model was eligible.** Sun G, Liu W. *Folia Neuropathologica* 2024;62(1):59-75,
PMID 38174685: a 7-gene risk score with every coefficient in the Results, a stated
training scale (TCGA-GBM, log2(FPKM + 1)), and CGGA mRNAseq_693 used for validation only.
It was one of the ten library PDFs. All 7 genes are present in both CGGA batches, and all
218 patients of the frozen partition have a column in CGGA's FPKM release; this was a
check of gene names and sample identifiers only, not of any value
(`docs/t2_comparator_search/README.md` section 3). So 1 of 549 records (0.18%) was
re-runnable as published. Its screening record is one full-text reviewer and the lead's
own reading. [LEAD: the README says the second review included "the two that passed".
The record for this paper shows no second independent reviewer. Confirm whether it had
one; if not, say so here.]

**Failure reasons by first failing criterion.** Table 1 gives the 548 records that did
not qualify. Counts are derived from `records.jsonl`; the E3 sub-reasons are hand-coded.

*Table 1. First failing criterion of the 548 ineligible records.*

| First failing criterion | Stage at which it was applied | Records | What failed |
| --- | --- | --- | --- |
| E1 | Date and type (69); title and abstract (1) | 70 | 69 outside 2024 or a preprint or retracted, by PubMed metadata [LEAD: the file does not split these reasons; split them from the `reason` field if wanted]; 1 not a research article |
| E2 | Title and abstract | 330 | Not a per-patient survival prediction for adult grade 4 glioblastoma (for example, wet-lab or mouse work, non-survival outcomes, other tumours) |
| E4 | Title and abstract | 61 | Needs an input the target cohort lacks (for example imaging, slide images, methylation arrays, performance status) |
| E3 | Full text | 80 | See Table 2 |
| E3 (not established) | Full text not obtained | 5 | Not read; 3 unobtainable, 2 failed E5 on the abstract |
| E5 | Full text | 2 | CGGA used for gene selection and fitting (1); CGGA used to choose the cut-off (1) |
| **Eligible** | | **1** | |
| **Total** | | **549** | |

*Table 2. The 85 records whose first failing criterion was E3 (hand-coded from the
`reason` field).*

| Sub-reason | Records |
| --- | --- |
| Coefficients available (text or supplement) but the expression scale not stated, or not reproducible (batch-corrected jointly, normalised per cohort) | 31 |
| No usable coefficients or fitted model (figure or chart only, generic formula, request only, classifier not released, unweighted or cohort-relative score) | 44 |
| Not an expression-based risk model that could be applied to CGGA (no prognostic model, clinical or mutation inputs) | 5 |
| Full text not obtained | 5 |

Of the 31 with coefficients but no stated scale, 7 were also flagged for CGGA training or
selection (E5). The other 24 were not. A looser reading of E3 would have kept those 24 in
contention for further checks. We did not assess them further. [LEAD: decide whether to
report this sensitivity as a finding or leave it to the Discussion, where it now sits.]

**Use of CGGA.** The title-and-abstract stage did not record whether CGGA was used, so
this count covers only the 88 papers that reached the full-text stage.

| Role of CGGA | Papers |
| --- | --- |
| Training, gene or feature selection, or cut-off choice (E5 failed) | 21 |
| Validation only | 22 |
| Descriptive analysis of single genes only | 5 |
| Not used, or not known | 40 |
| **Any use** | **48 of 88** |

The 21 are made of 17 flagged by the first full-text reviewer, 1 by adjudication (38378721,
cut-off), and 3 from an abstract or a page summary (39476302, 38662286, 38376452). Of the
22 validation-only papers, 9 did not name the CGGA batch and 1 used the microarray batch
only; so at most 21 of the 22 can have used the RNA-seq batches in question. The
registered search text says "17 papers" used CGGA for training or selection; that is the
first-reviewer count and is one of the three sources above. The "not used" group
includes the three papers whose full text was not obtained.

**Second review.** Twenty-five papers had two independent full-text reviewers
(22 "independent", 3 "adjudicated"). One ruling went against both reviewers: PMID
39435104, a 22 gene-pair model, was ruled ineligible under E3 by the author, because the
paper does not say whether its values were raw counts or length-normalised, and a pair
comparison depends on that. Another, 39630160, was ruled ineligible because the paper
describes its scale in two contradictory ways; the two reviewers had split on it.

### 1.4 Discussion

**Main finding.** Of 549 records from one 2024 PubMed search, one model could be re-run
as published on the external cohort we needed. The dominant barrier at the full-text
stage was not the absence of a model but the absence of a stated scale or of released
coefficients. Eighty-five of 88 papers that reached the full-text stage failed E3 first,
and 44 of those gave no usable coefficients or fitted model.

**What "reproducibility" means here.** The criteria were written to answer one question:
can this model be re-run, without refitting, on two CGGA RNA-seq batches? E4 and E5 are
about that cohort. A paper that fails them may be reproducible on its own data. The
audit does not measure whether a published result replicates. It measures whether a model
is available for use as a fixed comparator on one named cohort.

**Limitations.**

1. *E3 is strict.* It requires, for the data the coefficients were fitted on, the unit
   and every transformation applied before fitting. Values batch-corrected jointly with
   other cohorts, and coefficients readable only from a figure, count as no stated scale.
   The scale rule was made precise after the first full-text pass, then applied the same
   way to every later record. It interprets registered E3 and adds no criterion. The
   effect is large: 31 of the 85 E3 failures gave coefficients but no stated scale, and
   24 of those were not flagged for CGGA use. A researcher willing to infer the scale
   (for example, log2(FPKM + 1) from the pipeline) would have admitted more candidates.
   We did not, because an inferred scale is a choice, and a wrong choice changes the
   score.
2. *Three papers could not be obtained.* PMIDs 38959943, 38318849 and 37539935 are
   recorded as "E3 not established". They are not known to be ineligible. Two more
   papers were judged on their abstracts.
3. *The screening was AI-assisted, with human-defined criteria.* The criteria and the
   rulings are the author's. The screeners and reviewers were AI agents following written
   protocols. [LEAD: confirm.] The quotations were checked against the text by script,
   which shows that a quoted sentence exists, not that it was read correctly. If the two
   independent reviewers were the same model, their agreement is weaker evidence than
   agreement between two humans. The record shows no independent human screening of the
   480 abstracts or 83 full texts. [LEAD: confirm, and consider a human re-screen of a random sample to
   estimate error rates. The record has none.]
4. *Only 25 of 83 full texts had a second reviewer.* The 25 were the closest calls. A
   paper wrongly failed by a single reviewer, on a call that did not look close, was not
   re-read. We have no estimate of that error rate.
5. *One search engine and one query.* PubMed only, with one query that requires the word
   "glioblastoma" in the title or abstract and one of the listed model terms. A paper
   using only "GBM" or "glioma", or indexed under another date, would be missed. No
   second database, no citation chasing and no check of the search's recall were done.
6. *A single publication year.* 2024 is the registered choice. It tells nothing about
   other years.
7. *Coded groupings.* The CGGA-role and E3 sub-reason groupings are one reader's coding
   of free text. The counts in Table 1 and the flow are fields in the file; Table 2 and
   the CGGA table are not.
8. *Rules set while the search ran.* The author's ruling on 39435104 reversed both
   reviewers, and the scale rule was sharpened after the first pass. Both are disclosed
   in `docs/t2_comparator_search/README.md` section 5. Both make eligibility harder.
9. *The comparator rule was amended after the search.* The registration was amended
   the same day to a re-implementation (option B), after the search found no eligible
   model among the records then read, then amended again to restore the rule as first
   registered for the verdict once 38174685 was found among the 13 records that were
   still unread. The reasons are in `docs/preregistration.md` section 8 (entries of
   2026-11-27 on 2.2 item 2). Option B stays as a secondary, descriptive comparison.

**What this does not show.** It does not show that the one eligible model is good. It has
not been scored on CGGA. [LEAD: update when the improvement row is read.] It does not
show that the other models are poor, only that they cannot be used as a fixed comparator
here. It does not estimate how many models would be re-runnable by a method that refits.

**Implication.** Release of the fitted formula with its expression scale, or of code with
weights, is what would have made more of these models usable. Whether that is a general
failing of the 2024 literature, or a feature of this one search and these criteria, this
audit cannot say.

---

## Part 2. Benchmark establishment: frozen external partitions for GBM survival prediction

### Framing

The registered fallback (`docs/preregistration.md` section 7) says that if the external
survival C-index is below 0.65, the contribution is re-framed as benchmark establishment
for external GBM survival prediction on frozen public partitions, and the threshold is
reported as FAIL, not re-stated. The C-index row is FAIL (section 2.3 of the registration;
`docs/cgga_validation.md`). The registration records the decision that section 7's
fallback applies (section 2.3, "Decisions on T2's open items", item 2).

Scope. This part claims a benchmark: frozen partitions, a read discipline, and
pre-registered thresholds, with the results obtained so far. It does not claim a model
that works. The registration states the C-index shortfall was the most likely outcome in
advance (section 5, "Prior expectations", and section 7).

All results below were obtained before the OSF deposit and are marked "before deposit".
Eleven thresholds are registered (section 2 and section 6 rule 3). Three have been
evaluated: two FAIL (both T2) and one PASS (T1 deconvolution, on simulated data). Eight
have not.

### Outline

1. The problem: external validation of GBM survival models, and why frozen partitions
   and a read discipline matter.
2. The frozen partitions and the read discipline.
3. T2 so far: what was measured, the results, and the KPS non-evaluability.
4. T1 state deconvolution (simulated; optimistic).
5. What remains: T2 improvement row, T1 subtype, T3, T4 and T5.
6. Deviations and what was obtained before deposit.
7. Limitations.

### 2.1 The frozen partitions and the read discipline

*Partitions.* A partition is a pure function of `(seed, patient_id)`: a patient goes to
the test side if the SHA-256 of `"{seed}|{patient_id}"` maps below the test fraction. It
is independent of row order, of the number of rows, and of other patients. In the
registration's check on the first TCGA cohort, shuffling the rows four ways gave the
identical split, adding or removing 50 patients moved 0 patients, and the same change
applied to a shuffle-based split moved 173 patients (section 4.1).

| Partition | Cohort | n | Split | Role |
| --- | --- | --- | --- | --- |
| `tcga_gbm_clinical_v2` | TCGA-GBM clinical | 594 | 425 train, 169 test (seed 20251126) | Development. Not an external test. |
| `cgga_external_v1` | CGGA | 218 | 0 train, 218 test (all held out) | External test for T2, T4, T5. Sealed 2026-11-27. |

`cgga_external_v1` holds the 218 CGGA patients with an RNA-seq profile, a survival
outcome and a usable identifier: 133 from mRNAseq_693 and 85 from mRNAseq_325. Its
manifest has SHA-256 `040ba322d77f1f2cc205993d26c2aab63f28de4c86a8f2cc6bd3c056e8c3829e`
(registration section 4.5). CI re-derives all 218 assignments from `(seed, patient_id)`.
The superseded v1 cohort (583 patients) is kept as a record. A GLASS partition was not
created and is not used (section 4.5).

*Read discipline.* From the registration (sections 3.3, 3.9, 4):

- external partitions are read once per registered analysis;
- a sealed manifest refuses to be overwritten, and verification runs before every
  evaluation;
- a read guard refuses any read of values that is not registered, and every access of
  any kind is recorded in an append-only ledger (`data/access_ledger.jsonl`);
- every rule and script for a read is committed before the read, and each verdict comes
  from the registered threshold function, not from judgement;
- a synthetic dry run runs before each read;
- the threshold table is frozen by SHA-256
  (`96bbcd49dcda20f8a0bc9c8b71427aaa24cf754dc6d8d7b4f7d6a2df3d2874ab`).

There is one analyst, so no analysis is blind to design or outcome. The safeguards are
procedural, and the registration says so (section 3.9).

[LEAD: a sentence on the standing rule `.clinerules`, adopted after a workspace-wide
search opened held-out files on 2026-11-27 (access ledger entries 8 and 9), would show the
discipline working; decide whether it belongs in a methods paper or only in the record.]

### 2.2 T2 results so far (before deposit)

The registered four-feature model (age, sex, KPS, treatment arm) cannot be scored on
CGGA, which records no performance status. KPS was not imputed on CGGA. The registered
T2 model is therefore NOT EVALUABLE on the registered external cohort. The two rows
were decided by the nearest evaluable model, variant 2 (age, sex, treatment arm), which
was pre-specified as the headline in the script before the first read. All three
variants were fixed before the first read and all are reported.

| Variant | Features | CGGA n (deaths) | C-index, external (95% CI) | Brier at 18 months (95% CI) | Reported as |
| --- | --- | --- | --- | --- | --- |
| 1. Committed baseline | age, sex, KPS, treatment arm | none | not evaluable | not evaluable | NOT EVALUABLE |
| 2. Transferable (decides) | age, sex, treatment arm | 283 (239) | **0.5385** (0.4994 to 0.5807) | **0.2617** (0.2391 to 0.2829) | FAIL on both rows |
| 3. Molecular (secondary) | variant 2 plus IDH, MGMT | 241 (204) | 0.5915 (0.5477 to 0.6341) | 0.2472 (0.2226 to 0.2729) | C-index below threshold; Brier INCONCLUSIVE |

Thresholds: C-index lower bound at least 0.65; Brier at 18 months upper bound at most
0.25. Intervals are percentile bootstraps over CGGA patients (1000 resamples, seed 42)
with the TCGA-fitted model held fixed. Variants 2 and 3 are scored on different patients
and are not comparable (`docs/cgga_validation.md` section 4.3).

Readings the record supports:

- The C-index interval for variant 2 lies wholly below 0.65 and contains 0.50. The data
  are consistent with no discrimination and with weak discrimination. The record commits
  to neither.
- On the same 283 patients, 123 were alive and under observation at 18 months and the
  observed proportion alive at 18 months was 0.457. A model predicting that value for
  everyone scores a Brier of 0.2481. Variant 2 scores 0.2617. The score is not better than
  a constant prediction. Whether this is an offset or poor ranking was not examined.
- The registered Brier threshold, 0.25, is close to the uninformative score in this
  cohort (0.457 x 0.543 = 0.248), which the record notes.
- For context only, internal five-fold C-indices on TCGA are 0.6735 (variant 1), 0.6520
  (variant 2) and 0.6606 (variant 3). The step from 0.6520 (TCGA) to 0.5385 (CGGA) is
  0.1135 and the step from variant 1 to variant 2 on TCGA is 0.0215. Neither difference
  carries an interval. The 0.0215 bundles the KPS value, a prognostic missing-KPS
  indicator (HR 1.57) and a different fit set. The record does not attribute the
  shortfall to a cause.
- A proportional-hazards test, added to the registration after these reads, found that
  the TCGA fit of variant 2 violates proportional hazards in age (p = 0.0002; global
  p = 0.0006) (registration section 3.8).

[LEAD: confirm the order in which section 3.8 was added relative to the T2 reads: the
record says the finding came from the synthetic dry run before the deposit. It did not
alter any T2 number.]

Not done: calibration slope and intercept, calibration plots and decision-curve analysis.
They were decided against for this project, because each is a further read of CGGA and
none decides a registered row (`docs/cgga_validation.md` section 12.5). The Brier scores
at 12 and 24 months are descriptive only (variant 2: 0.2319 and 0.2254; null model:
0.2277 and 0.2202).

### 2.3 T1 state deconvolution (before deposit; simulated)

T1's state-deconvolution row was run once on 2026-11-27 on 100 pseudo-bulk mixtures
simulated from the Neftel 2019 atlas (GEO GSE131928), from a specification committed
before the run (`docs/t1_deconvolution_spec.md`). The estimator is the project's own
marker-enrichment scorer, not CIBERSORTx. The metric is the mean of the four per-state
Pearson r values, with the interval a percentile bootstrap over mixtures.

| Metric | Value | 95% CI |
| --- | --- | --- |
| Mean of four per-state r (decides; threshold at least 0.75 on the lower bound) | **0.9521** | 0.9430 to 0.9589 |
| NPC-like | 0.9647 | 0.9472 to 0.9768 |
| OPC-like | 0.9528 | 0.9290 to 0.9682 |
| AC-like | 0.9204 | 0.8885 to 0.9447 |
| MES-like | 0.9705 | 0.9574 to 0.9797 |

Verdict from the registry: PASS. The registration calls this test optimistic by design,
and the specification says why. The mixtures hold malignant cells only. The truth is
computed from the same cells. The atlas states and the estimator's marker panels both
come from the Neftel meta-modules, so some agreement is built in. The result says the
estimator recovers simulated state proportions under these conditions. It does not say
how well it estimates composition in a TCGA or CGGA tumour. A secondary, descriptive
analysis with a cohort-independent normalisation gives 0.8541 (95% CI 0.8358 to 0.8750).
No cohort was read for this row. Subtype accuracy, the other T1 row, has not been run, so
T1 as a whole is not reported as met.

### 2.4 What remains

| Row | State | What is needed |
| --- | --- | --- |
| T2 improvement over a 2024 baseline | Not evaluated; one read after the deposit | The comparator is PMID 38174685 (Part 1), applied as published to CGGA's FPKM release. Batch-stratified C-index, paired percentile bootstrap of the difference (1000 resamples, seed 42), lower bound at least 0.03. A refitted gene-set comparator (option B) is reported descriptively and carries no verdict. |
| T1 subtype accuracy | Not run | The 108 CGGA microarray patients with a reference Verhaak label, frozen as their own manifest; nearest-centroid classifier with frozen TCGA centroids; exact Clopper-Pearson interval; threshold at least 0.88 on the lower bound; NOT EVALUABLE below 30 labelled patients. The microarray matrix has not been downloaded. |
| T3 | Not run; data not obtained | PDX growth curves and in vivo NSC tracking. Data in hand by 2027-01-15 and calibration complete by 2027-03-01, or NOT EVALUABLE. |
| T4 and T5 | Not run; contingent | A simulator, a response surface and a survival model, frozen and versioned by a dated amendment before any T4 or T5 read; closed as NOT EVALUABLE if not frozen by 2027-06-30. |

The shortlist document records the expectation that the improvement row will be FAIL,
because clearing a paired difference of +0.03 at the lower bound would need a comparator
at or below chance on CGGA (`docs/t2_comparator_shortlist.md` section 2 point 5). That
expectation does not change the registered analysis. [LEAD: the shortlist was written
before the comparator was named. Decide whether to keep a forward-looking sentence in a
draft that will be read after the result.] T4 and T5 results, if run, will be predicted
RMST gains only. Section 2.2 item 6 of the registration says no T5 result may be called
a survival benefit, and that its credibility is bounded by T2's result.

### 2.5 Deviations and what was known before the deposit

Disclosed in the registration (sections 2.3, 2.4, 5, 8). None changed a threshold, a
decisive bound or a decision rule.

1. *T2 was read twice, on the ingested cohort as a whole, before the external partition
   existed.* The C-index row and then the Brier row (two reads, 2026-11-26). The
   partition was sealed on 2026-11-27 after both. Section 3.3 of the registration says
   final T2 numbers are reported on the frozen partition. They were not recomputed on it.
   The decision of record is that these two reads are the reads of record, and T2 is not
   re-run.
2. *The evaluation set differs from the partition.* The T2 numbers use 283 patients
   (variant 2) and 241 (variant 3) of CGGA's 305 unique patients (333 rows). Patients
   missing a feature or the outcome were excluded, not imputed. The sealed partition has
   218.
3. *T1's deconvolution row was run before the deposit.* Its result was known at deposit.
4. *The T2 comparator rule was amended, twice, after the search was run and before the
   deposit,* as described in Part 1 (limitation 9), with every CGGA access for this
   row recorded in the ledger (entries 11, 13, 16, 17). No CGGA expression value, clinical
   field or outcome has been read for the improvement row.
5. *The CGGA aside.* An unplanned workspace-wide search on 2026-11-27 scanned five
   extracted CGGA text files, returned nothing and was recorded as access ledger entries
   8 and 9. It used no value. [LEAD: confirm you want this in the paper; it is in
   `docs/cgga_validation.md` section 9.]
6. *Open items at the deposit.* The registration header has no DOI and no deposit date
   [LEAD: fill], and the "confirm at deposit" paragraph of section 5 is not yet signed
   off.

### 2.6 Limitations

- Benchmark scope. The partitions are two public cohorts, one development (TCGA-GBM) and
  one external (CGGA, Chinese patients, median age 52 against 59 in TCGA per
  `docs/cgga_validation.md` section 4.2). A shortfall between them is not explained here.
- The T2 model that could be scored uses three clinical variables. It is not the
  project's registered four-feature model. The benchmark number is for that stand-in.
- The sample is small and the intervals are wide: a C-index interval about 0.08 wide.
- Single analyst, no blinding. No calibration analysis.
- The one improvement-row comparator is a published gene set applied as published, on a
  cohort whose RNA-seq batches were processed separately; the registered batch-stratified
  C-index addresses this.

---

## 3. Numbers used, and their sources

Counts in Part 1 are computed by the script in the Appendix from
`docs/t2_comparator_search/records.jsonl` (cited as "records"). "README" is
`docs/t2_comparator_search/README.md`; "prereg" is `docs/preregistration.md`; "cgga" is
`docs/cgga_validation.md`.

| Number | Where used | Source |
| --- | --- | --- |
| 549 records; query; run time 2026-11-27T17:48:54Z | 1.2, 1.3 | `search.json`; records (count) |
| 69 fail E1 at date and type; 480 continue | 1.3 | records (script) |
| 392 fail at title and abstract (1 E1, 330 E2, 61 E4); 88 continue | 1.3 | records (script) |
| 83 full texts read; 5 not read; 3 unobtainable; 2 judged on abstract | 1.3 | records (script); README section 4 |
| E1 total 70; E2 330; E4 61; E3 85 (80 + 5); E5 2; eligible 1 | Table 1 | records (script); README section 3 |
| E3 sub-reasons 31, 44, 5, 5 | Table 2 | records `reason` field, hand-coded (script) |
| 24 and 7 of the 31 (E5 flag) | 1.3, 1.4 | records (script) |
| 48 of 88 any CGGA use; 21, 22, 5, 40 | 1.3 | records `e5_*` fields, hand-coded (script) |
| 17 first-reviewer E5 fails; 1 adjudicated; 3 from abstract or page summary | 1.3 | records (script) |
| 9 of 22 batch not named; 1 microarray only | 1.3 | records `e5_detail` (script); 39309926 |
| 4 title-and-abstract disagreements adjudicated | 1.2 | records (script); README section 2 |
| 71 open-access texts; 2 publisher pages; 10 library PDFs | 1.2 | README sections 2 and 4 |
| Second review: 25 papers (22 + 3); composition 2, 6, 16, 1 | 1.2, 1.3 | records (script); README section 2 |
| 58 single-reviewer full-text records read | 1.4 | records (script) |
| 7 genes in both CGGA batches; 218 patients have a column | 1.3 | README section 3 |
| E3-strictness rulings (39435104, 39630160) | 1.3 | README section 5; records |
| Eleven thresholds; three evaluated | Part 2 framing | prereg section 2, section 6 |
| TCGA v2: 594; 425 / 169; seed 20251126 | 2.1 | prereg section 4.5 |
| `cgga_external_v1`: 218; 133 and 85; manifest SHA-256 `040ba322...829e` | 2.1 | prereg section 4.5; `docs/project_record/06_evidence_reports.md` (sealing) |
| 173 moved by a shuffle split; 0 moved by hashing | 2.1 | prereg section 4.1 |
| Threshold checksum `96bbcd49...874ab` | 2.1 | prereg header and section 2.1 |
| Variant 2: n 283, 239 deaths; C-index 0.5385 (0.4994 to 0.5807); Brier 0.2617 (0.2391 to 0.2829) | 2.2 | prereg section 2.3; cgga sections 1 and 12.1 |
| Variant 3: n 241, 204 deaths; C-index 0.5915 (0.5477 to 0.6341); Brier 0.2472 (0.2226 to 0.2729) | 2.2 | prereg section 2.3; cgga sections 1 and 12.1 |
| Internal C-index 0.6735, 0.6520, 0.6606; differences 0.0215, 0.1135 | 2.2 | prereg section 2.3; cgga section 4.2 |
| 123 alive at 18 months; proportion 0.457; null Brier 0.2481 | 2.2 | cgga section 12.2 |
| 0.457 x 0.543 = 0.248 | 2.2 | cgga section 12.2 |
| Brier at 12 and 24 months: 0.2319, 0.2254; null 0.2277, 0.2202 | 2.2 | cgga section 12.3 |
| KPS HR 1.57; PH test p = 0.0002, global 0.0006 | 2.2 | prereg section 2.3 and 3.8 |
| 305 unique patients in 333 rows | 2.5 | prereg section 2.3; cgga section 8 |
| T1 mean r 0.9521 (0.9430 to 0.9589); per state; secondary 0.8541 (0.8358 to 0.8750) | 2.3 | prereg section 2.4; `docs/project_record/06_evidence_reports.md` |
| 100 mixtures | 2.3 | `docs/t1_deconvolution_spec.md` section 4; prereg section 2.4 |
| 108 labelled microarray patients; floor 30; thresholds 0.88, 0.03, 0.65, 0.25, 0.75 | 2.4, 2.2 | prereg sections 2, 2.2 item 9 |
| T3 dates 2027-01-15 and 2027-03-01; T4/T5 closing date 2027-06-30 | 2.4 | prereg section 2.2 items 3 and 11 |
| Median age 52 (CGGA) against 59 (TCGA) | 2.6 | cgga section 4.2 |
| Ledger entries 8, 9, 11, 13, 16, 17 | 2.5 | prereg section 5; cgga section 9 |

---

## Appendix. The counting script

Run as `python -I -X utf8 audit_counts.py docs/t2_comparator_search/records.jsonl`. It
reads that one file only. The PMID sets `VAL`, `DESC`, `S`, `N` and `X` are hand-coded
from the `e5_detail` and `reason` fields. [LEAD: decide whether to commit this script
under `scripts/` with a test; this draft does not add it to the repository.]

```python
"""Counts for the audit. Reads docs/t2_comparator_search/records.jsonl only."""
import json, sys, collections
C = collections.Counter
R = [json.loads(l) for l in open(sys.argv[1], encoding="utf-8")]
n = len(R); assert n == len({r["pmid"] for r in R})
print("records", n)
for k, v in sorted(C((r["stage"], r["outcome"], r["first_failing"]) for r in R).items(), key=str):
    print("  ", k, v)
dt = sum(r["stage"] == "date_and_type" for r in R)
ta = sum(r["stage"] == "title_abstract" for r in R)
ft = sum(r["stage"] == "full_text" for r in R)
no = sum(r["stage"] == "full_text_not_obtained" for r in R)
print(dt, n - dt, ta, n - dt - ta, ft, no)
f = [r for r in R if r["stage"] in ("full_text", "full_text_not_obtained")]
e5 = {r["pmid"] for r in f if r.get("e5_first_reviewer") == "fail"
      or r.get("e5_from_available_text") == "fail" or r["first_failing"] == "E5"}
VAL = set("39816535 39697739 39630160 39465792 39364731 39309926 39286023 39253064 39162363 39068393 38890658 38841168 38719174 38687369 38674418 38606090 38419527 38365809 38320988 38255198 38214842 39610810".split())
DESC = set("39736337 39262462 38516683 38763954 38322416".split())
S = set("""39896073 39697704 39697698 39630160 39610810 39539547 39465756 39435104 39364731 39322177 39309926 39286023 39253064 39162363 39145097 39136841 39071633 39033204 38978755 38813229 38763954 38665430 38655339 38638676 38566075 38462406 38365809 38339384 38279111 38255198 38214842""".split())
X = set("39857625 39525666 39265619 38733169 38698847".split())
N = set("""39882240 39816535 39736337 39697739 39640533 39640361 39465792 39458959 39444320 39387496 39333603 39262462 39254749 39216044 39216042 39156868 39144445 39072291 39068393 38927654 38890658 38841168 38832259 38817361 38809929 38765057 38719174 38695342 38687369 38674418 38615384 38606090 38595980 38516683 38510250 38495835 38434266 38419527 38376452 38369937 38335383 38322416 38320988 38249808""".split())
e3 = {r["pmid"] for r in R if r["stage"] == "full_text" and r["first_failing"] == "E3"}
assert e3 == S | N | X and not (S & N) and not (S & X) and not (N & X)
print(len(e5), len(VAL), len(DESC), len(S), len(N), len(X), len(S & e5), len(S - e5))
```

Output (verified when this draft was written): 549 records; 69 and 480; 392 and 88; 83
and 5; first-failing totals E1 70, E2 330, E3 85, E4 61, E5 2, eligible 1; E5 flagged 21
(17 first reviewer, 3 abstract or page summary, 1 adjudicated by cut-off, with no overlap between the three); CGGA validation only 22, descriptive only 5, any
use 48 of 88; E3 sub-reasons S 31, N 44, X 5; S flagged for E5 7, not flagged 24.
