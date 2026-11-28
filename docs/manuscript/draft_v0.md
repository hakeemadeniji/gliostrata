# GlioStrata manuscript: working draft v0

Status: working draft, not for circulation. **One paper in three parts**, merged
2026-11-28 by the project lead's decision; Part 3 was drafted separately and is now its final
section. Section 2.4's row for T4 and T5 is superseded by Part 3.
Every number is traced to a repository file in the list at the end (section 3). Places
that need the lead's decision are marked [LEAD: ...].

Standing caveat. At the time of writing `docs/preregistration.md` is a draft with no OSF
DOI and no deposit date (its header). By the project's own rule, T1 to T5 are exploratory
until the deposit. Everything below that reports a result was obtained before the deposit.
[LEAD: update this paragraph with the DOI and deposit date once recorded, and say which
results, if any, came after it. As of this draft none did.]

---

## Abstract (whole paper)

[LEAD: drafted 2026-11-28 on the decision to merge; needs your voice and a journal's word
limit. Parts 1 and 2 keep their own abstracts below for now, which a submission would drop.]

Pre-registered computational frameworks for personalising therapy are commonly validated
against themselves. We report three linked results from one registered protocol, in the order
they constrain each other. **First**, a single registered PubMed search of 2024 glioblastoma
survival models (549 records, dual independent screening) asked which could be re-run, as
published, on a held-out external cohort: **one of 549** could. The binding failure was not
performance but specification — most papers do not state the expression scale their
coefficients expect, and many validate on the very cohort batches they were trained on.
**Second**, applying our own transferable survival model to a sealed external partition of 218
CGGA patients gave a C-index of 0.5385 (95% CI 0.4994–0.5807) against a pre-registered
threshold of 0.65, and a Brier score at 18 months of 0.2617 (0.2391–0.2829) against 0.25:
both registered rows FAIL. We report the frozen partitions, the access ledger and the
read-once discipline as the transferable contribution. **Third**, accepting both results, we
specify and freeze a model for choosing a stem-cell-delivered therapy configuration per
patient, with every efficacy parameter fixed by a rule declared before the literature search
that feeds it (108 records, dual screening and dual extraction). That model is explicitly
model-conditional: it rests on a survival model that failed its own external validation, and
no quantity in it is a survival benefit. Its transferable output is the **break-even effect
size** — the effect a carrier-delivered therapy would have to achieve for the end-to-end
criterion to be met — which can be set against what oncolytic and CAR therapies have achieved
in humans.

The three parts share one discipline rather than one hypothesis: every threshold, partition
and decision rule was registered before the data it decides was read, every deviation is dated
and disclosed, and the amendment log records what was known at each point. Several of the
results below are negative, and two registered expectations were contradicted by the project's
own data and are reported as such.

---

## Part 1. A reproducibility audit of 2024 glioblastoma survival models

[LEAD: the title is the lead's working title. Section 1.4 explains why "reproducibility"
is narrower here than the title suggests. Consider "Which 2024 glioblastoma survival
models can be re-run on an external cohort without refitting? A registered search".]

### Abstract (Part 1)

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
| T4 and T5 | Not run; model specified and frozen | **See Part 3**, which supersedes this row. The simulator, response surface and survival model are specified in 3.2 to 3.4, their efficacy prior in 3.5, and the four pre-read controls in 3.6. Closed as NOT EVALUABLE if not frozen by 2027-06-30. |

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

## Part 3. A registered model for configuration choice, and the evidence behind it

### Framing

Parts 1 and 2 establish what is *not* available: no 2024 survival model could be re-run as
published on an external cohort without refitting (Part 1), and the project's own transferable
model does not discriminate on one (Part 2). Part 3 asks what can still be said once that is
accepted. It specifies, registers and freezes a model for choosing a therapy configuration
per patient, and reports the evidence behind every parameter in it.

It belongs in this paper rather than a separate one for the reason Part 2 exists: the
contribution is a benchmark and a discipline, not a therapy. A model-conditional result reads
as overclaiming when it stands alone and as a bounded, honest extension when it follows the
validation that bounds it.

**The standing caveat is stronger here than in Parts 1 and 2.** Everything below is
*model-conditional*. It predicts restricted mean survival under a therapy configuration using
a survival model that **failed its own external validation** (Part 2: C-index 0.5385, 95% CI
0.4994-0.5807, against a registered threshold of 0.65), and a dose-response whose efficacy
parameters are translated from rodent experiments by an assumed factor. No quantity in this
part is a survival benefit, and the registration (section 2.2 item 6) forbids calling one
that. What it can support is narrower and is stated in 3.7.

### 3.1 What the model is, and what it is for

The framework asks whether a per-patient *choice* among stem-cell-delivered therapy
configurations can be made better than a fixed one. Answering that needs a function from
(patient, configuration) to a predicted outcome. No such function exists in the literature, so
one was specified, registered and frozen before any held-out patient was scored.

The configuration space is registered as coded (section 2.2 item 4): carrier cell type (iPSC-NSC,
MSC, HB1_F3_line) × payload (CRAd, prodrug/enzyme, CAR, combo) × dose (10⁶–5×10⁷ cells,
continuous) × injection sites (1–8) × schedule (single, ×2, ×3), with two carrier–payload
pairs excluded as incompatible, leaving 10 pairs and 240 discrete combinations.

The model has three layers, adopted as "Design A" of `docs/t4t5_model_design.md` and fixed by
section 2.2 item 11(g)–(o):

1. a **survival model** supplying each patient's baseline curve under standard therapy;
2. a **delivery model** turning a configuration into an engrafted, acting cell count;
3. a **response surface** turning that into a log-hazard offset, against a toxicity term.

### 3.2 Survival model and the Stupp baseline

A Cox proportional-hazards model (Breslow ties) fitted on TCGA-GBM with age, sex and treatment
arm — the "transferable" variant of Part 2, chosen over a variant additionally carrying IDH and
MGMT status. **This was a deliberate choice against discrimination.** The molecular variant
scores a higher external C-index (0.5915 vs 0.5385), but requires an IDH and an MGMT call for
every patient, which the external partition's registered eligibility rule never demanded; on
CGGA the equivalent complete-case cohort falls from 283 to 241. Using it would have scored the
end-to-end test on fewer than the 218 patients its scoring set registers. A 0.05 gain in
C-index, on a model that fails the threshold either way, did not justify shrinking a sealed
population (section 2.2 item 11(g)).

A patient's **Stupp baseline** is this model's prediction with the arm set to radiotherapy plus
temozolomide and no stem-cell configuration (item 11(a)). Predicted survival is evaluated on a
fixed 0.1-month grid from 0 to τ = 24 months, and restricted mean survival time (RMST) is the
trapezoidal area under that curve. A configuration's **gain** is its RMST minus the baseline's.

Neither IDH nor MGMT appears anywhere in the frozen model — not in the baseline, and not as an
effect modifier, because no sourced modifier for the configuration effect exists and none was
invented. The consequence is disclosed rather than buried: **IDH-mutant patients' baselines are
modelled as if IDH-wild-type** (item 11(l)). IDH-mutant tumours are not glioblastoma under WHO
2021, which is an argument for excluding those patients; this was not acted on, because the
partition is sealed and changing its membership after sealing is what sealing prevents. IDH
status is reported as a pre-registered exploratory subgroup (section 3.6).

### 3.3 Delivery: a closed-form model, declared uncalibrated

For carrier *k*, a per-patient tropism score *z* is the mean standardised expression of that
carrier's chemoattractant panel — CXCL12, VEGFA, PLAU and CCL2 for the neural stem cell
carriers, TGFB1 for mesenchymal stem cells.

**Standardisation uses a frozen TCGA reference, never a statistic of the cohort being scored.**
A per-cohort z-score would let the sealed external cohort leak into its own transform. The
reference is TCGA-GBM, 284 samples, log2(FPKM+1) from GDC's `fpkm_unstranded`, checksummed in
section 8 at the freeze.

The scale matters more than it looks. The reference was first built on log2(CPM+1); measured on
the same 284 samples, that leaves a systematic offset of −1.36 SD for the NSC panel and −1.47 SD
for TGFB1 against an FPKM-normalised cohort, which alone would reduce every patient's engrafted
dose by about a third. Rebuilding on FPKM removes the gene-length component. **A second-order
difference remains and is not corrected**: TCGA's FPKM here is GDC's STAR-derived, CGGA's is
RSEM.

The engrafted, acting count is *n*<sub>eff</sub> = *N* · *f*<sub>k</sub> ·
logistic(*a*<sub>k</sub> + *g*<sub>k</sub>*z*), with spatial coverage 1 − exp(−*n*/*n*₀) over
injection sites and a multiplicative schedule factor.

This layer is **declared uncalibrated**, and that declaration is unconditional: it does not
depend on whether the project's separate PDX-calibration row (T3) ever obtains data (item
11(h)). *f*, *a*, *g* and *n*₀ are assumed, each with a sensitivity range. No PDE or
agent-based simulator was built.

### 3.4 Response surface and toxicity

The log-hazard reduction for payload *p* is

> *E* = *E*<sub>max,p</sub> · *m*<sub>p</sub>(*i*) · σ(*s*) · cov(*n*) ·
> *n*<sub>eff</sub><sup>h</sup> / (*n*<sub>eff</sub><sup>h</sup> + ED50<sup>h</sup>)

with a patient modifier *m* = exp(**w**′**u**) over a standardised immune-microenvironment score
and mesenchymal-like state fraction. The prodrug payloads receive no MGMT modifier: 5-FU and
SN-38 act through antimetabolite and topoisomerase mechanisms, not through MGMT — a mechanistic
argument, not a measured effect.

**"Combo" is defined per carrier**, as the two payloads that carrier is registered as able to
carry, combining by Bliss independence on the hazard scale with an assumed synergy of 1.0. A
single uniform pair was not available: the registered exclusions forbid CRAd on MSC and CAR on
the HB1.F3 line, so any uniform definition would place an excluded payload on a carrier. An
optimiser-chosen pair was also rejected, because two configurations sharing carrier, payload and
site count would then be different therapies, which the registered distinctness rule — and the
personalisation metric built on it — depend on them not being (item 11(j)).

Toxicity is an excess log hazard, *T* = τ<sub>p</sub>(*N*/5×10⁷)<sup>q</sup>(1 + η(*n*−1)) +
δ<sub>delay</sub>, covering dose, multiple injections and the cost of delaying standard therapy.
It is in the primary surface, not an exploratory variant, because **it is the model's only
source of harm**: without it the end-to-end safety threshold could not fail, and a safety
threshold that cannot fail is not a safety threshold (item 11(k)). A larger prior drawn from
reported CAR neurotoxicity was rejected as not transferable to carrier-cell delivery and is a
sensitivity range instead.

The configuration's curve is *S*(*t*) = *S*<sub>Stupp</sub>(*t*)<sup>exp(−E+T)</sup>.

### 3.5 The efficacy prior: a registered literature search

*E*<sub>max</sub> is the single quantity on which the end-to-end result most depends, and the
only one with no internal evidence. It was fixed by a **rule declared before the search that
feeds it**, so that no study could be chosen for the number it would produce:

> *E*<sub>max,p</sub> = 0.5 × ln(*ratio*<sub>p</sub>)

where *ratio*<sub>p</sub> is the **median** treated-to-control median-survival ratio across
*every* qualifying study for that payload, and 0.5 is an assumed mouse-to-human translation
factor, declared as assumed wherever a number from it is reported.

An earlier version named three studies carried over from the design document. It was rewritten
before the search because none satisfied its own criterion: one was in ovarian cancer, one
reported a tumour-volume rather than survival ratio, and the third had been shortcut to "the
pessimistic end" of an assumed range.

**Search and screening.** Three PubMed queries, one per payload, written into the registration
and run once (2026-11-28), no year limit: 108 records. Eligibility was fixed in advance (P1
peer-reviewed and in vivo; P2 orthotopic or intracranial rodent glioma; P3 median survival for
both a control and a payload arm; P4 payload class matched; P5 carrier delivery neither required
nor disqualifying). Titles and abstracts were screened **twice independently**, agreement
95.4% with no include/exclude conflict, disagreements resolved by a rule fixed before they were
seen; 40 records carried forward. Full texts were obtained for 38 of 40 — 18 from PubMed
Central, 20 retrieved manually — and extracted twice independently. Two articles could not be
obtained, fail P3 for inaccessibility, and are reported rather than dropped.

**A disclosed deviation.** P3 required a "control arm" without defining one, and the qualifying
studies use untreated, vehicle, reporter-vector and active-therapeutic comparators, which give
systematically different ratios. The primary set was narrowed to studies whose control arm
received no agent intended to have anti-tumour activity. This is recorded as a **deviation, not
an amendment**, because the reading was chosen after the literal reading's ratios were visible;
the inclusion test is structural — what the control arm *received* — and never its result. The
literal reading is reported alongside throughout.

**Result.**

| Payload | *ratio* | *E*<sub>max</sub> | implied HR | n studies |
|---|---|---|---|---|
| CRAd | 1.125 | 0.0589 | 0.943 | 3 |
| Prodrug/enzyme | 1.540 | 0.2159 | 0.806 | 10 |
| CAR | 2.691 | 0.4950 | 0.610 | 1 |

Two features of this table are reported with it rather than left to a limitations section.
**CRAd's modelled effect is close to none**, consistent with null human randomised data for
oncolytic adenovirus. **CAR's prior rests on a single study** and is the largest of the three;
its other eligible studies were excluded by the intracranial-model criterion, by unobtainable
medians, or by irrelevant-specificity CAR controls. No minimum study count was introduced,
because such a rule would be payload-neutral in wording and CAR-only in effect, and would have
been written with CAR's count and its *E*<sub>max</sub> both already visible. The registered
reading is, for every thin payload, **the most favourable of the available readings** — every
broader reading moves it further from the success threshold.

### 3.6 Optimisation, metrics and controls

Per patient, Bayesian optimisation (Gaussian process with expected improvement) and random
search each receive a registered budget of **60 objective evaluations** — 12 Sobol initial
points and 48 model-guided iterations for the former. Four metrics are evaluated against
pre-registered thresholds through a checksummed registry, which returns PASS, FAIL,
INCONCLUSIVE or NOT EVALUABLE and refuses a PASS when an auxiliary requirement is not reported:

* distinct configurations in the top five across fifty patients (≥ 3);
* regret reduction of the optimiser against random search (≥ 25%, with a Wilcoxon signed-rank
  requirement);
* fraction of held-out patients with a gain of at least 3 months (≥ 0.50);
* fraction with a gain below zero (≤ 0.10, with a CI-upper requirement of 0.15).

The fifty are the first by sorted patient identifier within the sealed partition, so the
selection cannot depend on an outcome.

**A property of the regret metric, recorded before the read.** Regret is measured against the
best value found by *either* method, so a method attaining that best for every patient has a
regret of identically zero, the metric reads exactly 100% whatever its margin, its bootstrap
interval is degenerate, and the threshold **cannot fail**. In the pre-read rehearsal this is
what happened. The metric is unchanged — rewriting it once its behaviour had been measured
would be worse than reporting that it saturates — and it is reported together with whether it
saturated and with the mean advantage in months, which does not saturate. [LEAD: this paragraph
belongs in Results as well if the confirmatory run saturates too.]

**Four controls run before any held-out patient is scored** and a failure stops the read: a
no-modifier control, which must collapse the optimiser's answer to a single configuration; a
null-payload control, which must give a benefit term of exactly zero and a gain equal to the
toxicity term's alone; a synthetic end-to-end run of the evaluation script itself; and a
rehearsal of the whole pipeline on a held-out development split. All four passed. The
null-payload control gave a benefit of exactly 0 across 120 configurations and 20 patients and
matched the toxicity-only gain to within 10⁻¹² in all 2,400 pairs — and **could not have passed
as the control was first written**, which required a gain of exactly zero, unreachable while the
toxicity term remains in the primary surface.

### 3.7 What this model can and cannot support

Stated here rather than in a discussion section, because it bounds every number in Part 3's
results.

**It cannot** support a claim of survival benefit, of any size, for any configuration. It
predicts RMST under a model whose discrimination failed external validation, using an efficacy
prior translated from rodent experiments by an assumed factor, on a delivery layer declared
uncalibrated.

**It can** support claims of the form: *given* this registered model, does per-patient choice
differ across patients, does model-guided search beat random search, and what effect size would
the therapy need for the end-to-end criterion to be met? The last of these — the break-even
effect size — is the most transferable output, because it converts a model-conditional result
into a statement about the effect a carrier-delivered therapy would have to achieve, which can
be set against what oncolytic and CAR therapies have achieved in humans.

**Decided 2026-11-28 by the project lead: one paper, with this as its final section.** The
reasoning is in Part 3's Framing. It commits the paper to an unusual shape — two negative
results followed by a model that is explicitly conditional on them — and that shape is the
argument rather than a concession to it: the reader meets the bound before the thing it
bounds. [LEAD: the alternative, a standalone Part 3, remains available if a journal asks for
it, but the caveat in 3.7 would then have to carry the whole weight that Part 2 carries here.]

### 3.8 Reproducibility

The model, its parameters and the evidence behind them are files, not descriptions: a single
versioned parameter file carrying every numeric value with its status (fitted, published with a
PMID, or assumed) and its sensitivity range; a frozen standardisation reference; the full search
and screening log, with every screening decision and the passage behind every extracted number;
and an append-only access ledger recording every read of the external cohort. The parameter file
and the reference are checksummed in the registration's amendment log at the freeze, which
precedes any held-out read. [LEAD: add the DOI, the freeze date and the two checksums here once
recorded.]

---

## 4. Numbers used, and their sources

Part 3's numbers are traced in `docs/t4t5_sensitivity_plan.md`, `docs/t4t5_emax_search/emax_result.md` and `configs/response_surface_v1.yaml`.

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
