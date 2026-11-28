# GlioStrata — OSF Pre-Registration Draft

> **Status: DRAFT, not yet deposited.** This document is the OSF-ready pre-registration
> for the computational tests T1–T5 of the GlioStrata protocol. It must be deposited on
> the Open Science Framework, and the deposit DOI recorded in the header below, **before
> any external test partition is read by a model**. Until that happens, the project has
> no pre-registration and T1–T5 are exploratory analyses, whatever the code says.
>
> **TO BE COMPLETED BY AUTHOR: OSF registration DOI** (insert once the registration is
> created; do not fabricate a placeholder DOI).
> **TO BE COMPLETED BY AUTHOR: date of OSF deposit (ISO-8601).**
> Both are entered after the deposit, in one commit, together with `configs/registration.yaml`.
> Every other header field below is complete.

| Field | Value |
| --- | --- |
| Study title | GlioStrata: A Patient-Specific Computational Framework for Optimizing Stem-Cell-Delivered Therapy in Glioblastoma |
| Registration type | Pre-registration of a confirmatory computational analysis (secondary analysis of existing public data) |
| OSF registration DOI | TO BE COMPLETED BY AUTHOR: OSF DOI |
| Date of registration | TO BE COMPLETED BY AUTHOR: ISO-8601 date |
| Authors | Hakeem Adeniji (sole author). MS Health Informatics student, Boston University; MS Psychology, Arizona State University; independent researcher. ORCID: not provided |
| Thesis committee / supervisor | Not yet appointed. Named by a dated amendment (section 8) when appointed |
| Planned start of the confirmatory phase | Immediately after the deposit. The first registered read (T2's improvement row, section 2.2 item 2 (j)) runs once the DOI is recorded in `configs/registration.yaml` |
| Code repository | `https://github.com/hakeemadeniji/gliostrata` (made public at the deposit) |
| Code commit registered | The commit tagged `prereg-deposit` in the repository. Its full SHA is entered on the OSF form and recorded in `configs/registration.yaml` (`registered_commit`) after the deposit; a document cannot contain its own commit's hash |
| Threshold table checksum | `96bbcd49dcda20f8a0bc9c8b71427aaa24cf754dc6d8d7b4f7d6a2df3d2874ab` (SHA-256 over `gliostrata.evaluation.thresholds`) |
| IRB status | Exempt — secondary analysis of existing, de-identified, public data (TCGA, CGGA, CPTAC, published single-cell atlases) |
| Funding / conflicts | None. No external funding; no conflicts of interest |

---

## 1. Research question and hypotheses

### 1.1 Central research question

Given a patient's tumour multi-omic profile and single-cell state composition, can a
computational framework predict a stem-cell-delivered therapeutic configuration — over
the joint space of carrier cell type, payload, dose, injection sites and schedule — that
yields a clinically meaningful survival gain over the Stupp protocol baseline?

### 1.2 Primary hypothesis (H1)

Patient-specific configuration selection produces a **predicted restricted mean survival
time (RMST) gain of at least three months in at least half of held-out test patients**,
with a **harm rate (predicted outcome worse than Stupp) below ten percent**.

Pre-registered operationalisation: metric `T5_endtoend_rmst_success_fraction` ≥ 0.50
(95% CI lower bound) and `T5_safety_predicted_harm_fraction` ≤ 0.10 point estimate with
95% CI upper bound ≤ 0.15.

### 1.3 Null hypothesis (H0)

No configuration in the search space, chosen by the optimizer, reliably exceeds the Stupp
baseline. Accepting H0 is a publishable, pre-registered negative result: it would indicate
that stem-cell-delivered therapy for GBM remains a uniformly applied experimental
intervention with no detectable personalisation signal in currently available public data.

### 1.4 Secondary hypothesis (H2) — the value of personalisation

If the optimizer converges to a single universal recommendation across all patients, the
framework has not personalised therapy; it has rediscovered the existing fixed-payload
approach under another name. H2 states that across a held-out cohort of fifty patients,
the top-five recommended configurations **include at least three distinct configurations**:
distinct in carrier cell type, payload or number of injection sites (section 2.2 item 4).

Pre-registered operationalisation: `T4_personalization_distinct_configs` ≥ 3.

H2 is logically prior to H1 in interpretation: an H1 pass with an H2 failure must be
reported as "a single configuration outperforms Stupp", **not** as evidence for
personalisation. This interpretive rule is registered here precisely so it cannot be
quietly dropped later.

### 1.5 Supporting hypotheses (component validity)

H3 (profiling, T1), H4 (survival prediction, T2) and H5 (mechanistic simulation, T3) state
that each upstream component clears its own pre-registered external threshold. They are
not merely diagnostics: an H1 result built on a component that failed its threshold is
reported as such, with the failing component named in the abstract.

---

## 2. Pre-registered thresholds (Protocol Table 8-2)

PASS requires **the relevant bound of the bootstrap 95% confidence interval** to clear the
threshold — the *lower* bound for "higher is better" metrics, the *upper* bound for the
harm and error metrics. A point estimate that clears the threshold is never sufficient.
This is the bar reviewers enforce and the bar this project is committed to.

| Test | Metric | Requirement | Decisive bound | Justification |
| --- | --- | --- | --- | --- |
| T1 | Subtyping accuracy (CGGA external) | >= 0.88 | 95% CI lower | Best published external subtyping accuracy is ~0.88-0.90. |
| T1 | State deconvolution Pearson r | >= 0.75 | 95% CI lower | CIBERSORTx-class deconvolution typically achieves r 0.70-0.85. |
| T2 | Survival C-index (CGGA external) | >= 0.65 | 95% CI lower | Best published external GBM survival C-index is ~0.65-0.68. |
| T2 | Improvement over best-published 2024 baseline | >= 0.03 | 95% CI lower | Demonstrate improvement, not parity, against the 2024 state of the art. |
| T2 | Brier score @ 18 months | <= 0.25 | 95% CI upper | GBM 18-month Brier baseline is ~0.23-0.27. |
| T3 | Tumor-growth RMSE vs PDX (% of final volume) | <= 15 | 95% CI upper | Standard PDX growth-modeling benchmark tolerance. |
| T3 | NSC-migration r vs in vivo | >= 0.65 | 95% CI lower | In vivo NSC tracking data is noisy, so the bar is set below T1's r. |
| T4 | Optimizer beats random search (% regret reduction) | >= 25 **and** Wilcoxon signed-rank p < 0.01 | 95% CI lower | A statistical test, not a point estimate, must carry the claim. |
| T4 | Personalization variance (distinct configs in top-5 across 50 patients) | >= 3 | 95% CI lower | Proves the optimizer has not converged to a universal recommendation. |
| T5 | End-to-end success (fraction of held-out patients with >= 3-month RMST gain) | >= 0.50 | 95% CI lower | Pre-registered primary success criterion from the working hypothesis. |
| T5 | Safety floor (predicted harm rate) | <= 0.15 (CI upper) **and** <= 0.10 (point estimate) | 95% CI upper | Hard safety stop: harm rate <= 10% with 95% CI upper bound <= 15%. |

### 2.1 Machine-readable form and tamper protection

The table above is generated from `src/gliostrata/evaluation/thresholds.py`, which is the
authoritative, machine-readable copy. Each row records the test id, metric name, comparison
direction, threshold value, which CI bound is decisive, and the protocol's justification.

The module computes a SHA-256 over the canonical serialisation of the whole table, and
`tests/test_evaluation.py::TestThresholdTableIntegrity::test_checksum_matches_the_frozen_constant`
asserts it against the frozen constant
`96bbcd49dcda20f8a0bc9c8b71427aaa24cf754dc6d8d7b4f7d6a2df3d2874ab`. Any edit to any
pre-registered number, direction or CI bound breaks the build. Changing a threshold is
therefore necessarily a deliberate, reviewed, individually attributable commit, and must
also be recorded in §8 (Amendment log) of this document.

Evaluating a metric that is not in the registry raises `UnregisteredMetricError`. Inventing
a more favourable metric after seeing results is not merely discouraged here; it does not
run.

### 2.2 Ambiguities in the protocol table, recorded rather than silently resolved

These items are under-specified in Protocol Table 8-2. They are listed here so that the
resolutions are pre-registered rather than chosen after the data is seen.

1. **T1 state-deconvolution Pearson r — pooled or per-state?** The protocol gives a single
   r for four Neftel states (NPC-, OPC-, AC-, MES-like). *Registered resolution:* the
   primary metric is the mean of the four per-state Pearson r values across patients, with
   per-state r values reported individually alongside; the bootstrap resamples patients,
   not state-patient pairs. The pooled-across-all-fractions r is reported as a secondary
   descriptive number only, because it is inflated by between-state variance.
2. **T2 "beat best-published 2024 baseline by ≥ 0.03" — against what number?** A delta
   against a C-index quoted from a paper computed on a *different* cohort has no sampling
   distribution and no valid CI. *Registered resolution:* the comparator must be re-run on
   the identical frozen partition, and the metric is the paired bootstrap distribution of
   the difference. If the comparator cannot be re-run (no released code or weights), this
   threshold is reported as **not evaluable**, never as a pass. *Comparator (registered 2026-11-27 as a published model re-run without refitting; amended
   the same day, before the deposit, after the registered search had been run; section 8).*
   *As first registered:* a published model chosen by criteria E1 to E6 and one PubMed
   search, applied as published. E1: peer-reviewed, first published in 2024. E2: predicts
   overall survival in adult WHO grade 4 glioblastoma. E3: re-runnable without refitting,
   with every coefficient and the expression scale it expects. E4: every input available
   for `cgga_external_v1`. E5: no CGGA mRNAseq_325 or mRNAseq_693 patient used for training,
   gene or feature selection, coefficient fitting or cut-off choice. E6: the best external
   non-CGGA C-index. The search was run once, at 2026-11-27T17:48:54Z, and its query, protocols
   and per-record log are committed (`docs/t2_comparator_search/`). Of 549 records none was
   eligible, with 13 still awaiting their full text. Nearly every 2024 model with published
   coefficients fails E3: it does not state the expression scale its coefficients expect.
   *As amended, the comparator is a re-implementation under identical conditions*, which
   needs no published scale. (a) *Candidates:* the records of that search, every one judged
   by its committed protocols first, the 13 included. No other search is made. (b) *Eligible
   gene sets:* E1, E2 and E5 as above, and two replacements for E3, E4 and E6. R1: the
   published risk score is a function of bulk-tumour expression only, of genes or of
   within-sample gene pairs, and every gene and pair is named in the main text or in a
   supplementary file read in the review. R2: every gene is present under its symbol in both
   the TCGA expression matrix and the CGGA RNA-seq count matrix. The CGGA gene names are read
   as an availability access, recorded in the access ledger; no CGGA value is read for the
   choice. (c) *Re-implementation, identical for every gene set:* expression is log2(CPM + 1),
   computed for each sample from its STAR read counts, with the library size summed over the
   genes present in both matrices; a gene pair is 1 when its first gene's value exceeds its
   second's and 0 otherwise. The gene set's terms enter `CoxPHModel` (Breslow, ridge 1e-6, the
   estimator of section 2.3) as covariates. No clinical covariate is added and no term is
   selected or dropped. (d) *Choice:* each eligible gene set is scored by 5-fold stratified
   cross-validation (`cross_validate_cox`, seed 42, pooled out-of-fold Harrell's C-index) on
   the patients of the `tcga_gbm_clinical_v2` training split with an RNA-seq profile and an
   observed outcome. The comparator is the gene set with the highest score; ties go to the
   earlier publication. Neither the v2 test split nor CGGA is used for the choice. (e) *Final
   fit:* the chosen gene set is refitted on every TCGA patient with an RNA-seq profile and an
   observed outcome, as variant 2 was fitted on every complete-case TCGA patient. (f) *Our side
   of the comparison* (unchanged): the transferable model (variant 2 of section 2.3: age, sex,
   treatment arm), fitted on TCGA exactly as committed in `scripts/validate_cox_cgga.py`. It
   decided the C-index and Brier rows, because the committed model needs KPS, which CGGA does
   not record. (g) *Partition and procedure:* `cgga_external_v1`, complete cases on the union
   of both models' inputs. For each model, a batch-stratified Harrell's C-index: a pair of
   patients is comparable only when both are in the same CGGA RNA-seq batch
   (`cgga_dataset`, mRNAseq_693 or mRNAseq_325, as recorded in the eligibility table), and
   the concordant and comparable pairs are summed over both batches. The two batches were
   processed separately, so a technical shift between them would penalise an expression
   score but not the clinical model; within-batch pairs remove that asymmetry. Then a paired
   percentile bootstrap of the difference (variant 2 minus the comparator), resampling
   patients from the whole partition, 1000 resamples, seed 42. The lower bound of the 95%
   interval decides, through the registered metric `T2_cindex_improvement_vs_baseline`.
   The unstratified C-index of each model is reported alongside, descriptively. (h) *Outcomes:* if no gene set is eligible,
   this threshold is NOT EVALUABLE and the log is the stated reason. Otherwise the chosen gene
   set is named here by a further amendment before the deposit, and OSF stopping rule 9 and
   section 5 gain exactly this one read. The read runs after the deposit, on the sealed
   partition, through the read guard. (i) *What it measures:* the comparison is with a
   published 2024 *gene set* refitted under the conditions of our own model, not with the
   published model itself, and it is reported as such. The reasoning, the alternatives
   considered and the likely verdict are in `docs/t2_comparator_shortlist.md`. (j) *The search,
   completed, and the comparator named (2026-11-27, before the deposit).* The 13 records
   still unread when option B was adopted were then screened: 10 from library PDFs, 3 not
   obtainable. One record meets E1 to E5 of the rule as first registered: Sun G, Liu W,
   *Folia Neuropathologica* 2024;62(1):59-75, PMID 38174685, a 7-gene neutrophil
   extracellular trap signature. It was trained on TCGA-GBM log2(FPKM + 1), and CGGA
   mRNAseq_693 was used for validation only. The search found no other eligible model,
   so E6 does not arise. Because the rule fixed before the search yields a comparator,
   **that comparator decides the verdict**, applied as published:
   RS = 0.223340457·CPPED1 + 0.1517013·F3 + 0.046311722·G0S2 - 0.219913858·MAPK1 +
   0.192197058·MME + 0.016505054·MMP9 + 0.414575274·MPO, higher RS meaning higher risk. Each
   gene is matched by the symbol printed in the paper, and its value is log2(FPKM + 1) from
   CGGA's own RSEM FPKM release: CGGA.mRNAseq_693 and CGGA.mRNAseq_325
   .RSEM-genes.20200506, the release the authors validated on. That release was fetched
   before the deposit; its SHA-256 values are in `data/provenance/cgga_fpkm_sources.json`
   (access ledger entry 13). An availability check of gene names and sample ids only
   (entry 16) found all 7 genes in both batches and all 218 partition patients. The
   batch-stratified C-index of (g), variant 2 as in (f), and the paired bootstrap of (g)
   apply unchanged. *Option B is kept as a secondary, descriptive comparison:* the gene
   set chosen by (b) to (e) is scored in the same single read. Its difference from variant
   2 is reported with its interval but carries no verdict, and nothing about it changes the
   primary result. The comparison is one read of CGGA, made after the deposit through the
   read guard on the sealed partition. OSF stopping rule 9 and section 5 list it. *Option B's
   choice, made before the deposit:* 33 gene sets passed E1, E2, E5 and R1
   (`configs/t2_comparator/gene_sets.json`), and R2 removed 10 of them. R2 was read from CGGA
   gene names only (access ledger entry 17); the shared gene universe has 33,212 symbols. Of
   the 23 cross-validated on the 200 v2-training patients with RNA-seq (160 deaths), the
   best was the same 7-gene set as the primary comparator, PMID 38174685 (pooled out-of-fold
   C-index 0.615). The secondary comparison is therefore that gene set refitted under the
   project's own conditions, on 282 TCGA patients (228 deaths), frozen in
   `configs/t2_comparator/chosen.json` (SHA-256 `d728965692a140b95fbc2b9743631d36cb1fe551861dc7aa410e98f00fefe6e3`).
3. **T3 tumour-growth RMSE "≤ 15% of final volume" — whose final volume, and over which
   curves?** *Registered resolution:* RMSE is computed over the full time course per PDX
   growth curve, normalised by that curve's own observed final volume, then averaged across
   curves; the bootstrap resamples curves. *The dataset is not named at deposit,
   because no source has been verified; what is registered instead is the rule that
   fixes it.* (a) *Unit:* one curve is one PDX model, control arm, mean volume per
   time point, with at least three time points; the number of curves is the sample
   size. (b) *Qualifying rule:* a source is usable only if it gives longitudinal
   volumes, a control arm, and a subtype call or RNA-seq for every model, and it must
   carry at least ten curves spanning **at least three of the four Verhaak subtypes,
   with at least two curves per represented subtype**. (The requirement was four
   subtypes until 2026-11-27, when it was loosened to three with a per-subtype floor;
   see section 8. Any absent subtype is named in section 8 on the day the data are
   obtained, and the coverage limit is stated wherever T3's verdict is reported, so
   that no claim is made for a subtype that was never calibrated.) (c) *Candidates, none
   verified:* the Mayo Clinic Brain Tumor PDX National Resource and the Alcaniz 2023
   panel. (d) *How the unknowns are fixed:* the dataset and its curve count are
   entered in section 8 on the day the data are obtained, before any model is
   calibrated. (e) *Dates:* a qualifying source must be **secured by 2027-01-08** — an
   open release verified against (b), or a custodian's written agreement to supply one.
   If none is secured by that date, this criterion is declared NOT EVALUABLE on that
   date (item 8) and the tumour-growth layer is presented as uncalibrated; the
   declaration is not deferred past it, because the T4 and T5 freeze (item 11) has to
   know whether the simulation layer is calibrated. If a source is secured: data in
   hand by 2027-01-15; calibration complete by 2027-03-01. If the data are not in hand
   by 2027-01-15, this criterion is likewise NOT EVALUABLE.
4. **T4 personalisation variance is a count with no sampling distribution.** "≥ 3 distinct
   configs in top-5 across 50 patients" over one fixed cohort has no CI unless the cohort
   is itself resampled. *Registered resolution:* patients are bootstrap-resampled and the
   distinct-config count recomputed, giving a genuine CI lower bound. If no resampling is
   performed the CI is reported as degenerate (lower = point = upper) and labelled as such.
   Note also that "distinct configs" needs an equivalence rule: *registered resolution* —
   two configurations are the same if and only if they share carrier cell type, payload and
   number of injection sites, with dose and schedule ignored for this count. (The wording
   "injection-site set" was corrected on 2026-11-27, because the configuration space has a
   number of sites, not a set of locations.) *The configuration space, registered
   2026-11-27 as coded in `src/gliostrata/optimization/config_space.py`:*
   - carrier cell type: iPSC-NSC, MSC or HB1_F3_line;
   - payload: CRAd, prodrug_enzyme, CAR or combo;
   - dose: 1 x 10^6 to 5 x 10^7 cells (continuous);
   - number of injection sites: 1 to 8 (integer);
   - schedule: single, fractionated_x2 or fractionated_x3.

   Two cell type x payload pairs are excluded as incompatible: HB1_F3_line x CAR and
   MSC x CRAd. That leaves 10 cell type x payload pairs, so 240 discrete combinations with
   dose as the one continuous dimension.
5. **T4 regret reduction "≥ 25% lower regret" needs a regret definition.** *Registered
   resolution:* regret is the simple regret against the best configuration found by any
   method in the shared evaluation budget, at a fixed, pre-declared budget of **60
   objective evaluations per patient per method, initial design included**: GP-EI
   spends 12 points on its Sobol initial design and 48 on GP iterations, and random
   search makes 60 uniform draws. 60 is the optimizer code's existing default; it was
   not tuned on any data and will not be recalibrated on the external cohort.
   *Tie rule:* a patient whose regret is zero under both methods is a tie and is
   dropped by the signed-rank test. At most 20 tied patients are tolerated (40% of
   the 50 patients); if more are tied, the p-value is still reported but this
   threshold is INCONCLUSIVE whatever its value, and it is not re-run after dropping
   patients. That is a registered exception to section 6, under which a signed-rank
   p >= 0.01 is FAIL. The comparison is paired per patient and the Wilcoxon
   signed-rank test is two-sided. *Aggregate (registered 2026-11-27):* for patient i, regret_i
   under a method is the best objective value found by either method for i minus that
   method's best value for i, in months of predicted 24-month RMST (item 11). The
   metric is 100 x (1 - mean regret under GP-EI / mean regret under random search) over
   the 50 patients. Its 95% CI is a percentile bootstrap over patients (1000 resamples,
   seed 42), and the Wilcoxon test is on the per-patient differences.
6. **T5 success and safety are defined on *predicted* outcomes.** Both metrics are
   computed from the project's own survival model, so T5 is a statement about internal
   consistency of the framework, not about realised patient benefit. This is registered
   explicitly: no T5 result may be described as a survival benefit, only as a *predicted*
   RMST gain, and its credibility is bounded above by T2's result. The RMST horizon, the
   Stupp baseline, success and harm are defined in item 11.
7. **"Held-out patients" for T5 and the 50-patient cohort for T4 are not named in the
   protocol.** *Registered resolution:* both are drawn from the frozen external test
   partition (§4), never from the development cohort, and the T4 fifty are the first fifty
   by sorted patient id within that partition, so the selection is not outcome-dependent.
   *Eligibility for that partition, registered in advance:* CGGA patients with an
   RNA-seq profile, a survival outcome and a usable identifier mapping (at most 218
   patients). Every eligible patient is held out; the external cohort has no training
   split. The manifest is sealed after the T2 evaluation (two reads, section 2.3) and before any
   T4 or T5 computation, and that ordering is recorded in section 8.
8. **T3 depends on inputs that may never arrive.** NSC-migration validation requires in
   vivo tracking data that is published only in aggregate form in some sources. If the
   required data cannot be obtained, T3 is reported as **not evaluable**, with the reason
   stated, and the framework's simulation layer is presented as uncalibrated. This is
   distinct from a FAIL. *Decision date, registered 2026-11-27:* both T3 rows are
   decided on **2027-01-08**. By that date each row must have a qualifying source
   secured — for the tumour-growth row as item 3(b) and 3(e) define it, for the
   NSC-migration row in vivo tracking data giving per-animal or per-region distances
   rather than a single aggregate — or that row is declared NOT EVALUABLE on that date.
   The two rows are decided independently: one may be calibrated while the other is
   not. The date is earlier than item 3's data deadline on purpose, so that the T4 and
   T5 freeze (item 11) is made knowing T3's status, and it is not extended, because an
   open-ended search would let the decision be taken after seeing which answer suited
   the framework.
9. **T1 subtype accuracy: which CGGA patients are scored, with which interval, and
   below which n the criterion cannot be decided.** CGGA's reference Verhaak labels
   exist only in its microarray batch, so the scoring set decides whether the 0.88
   threshold can be decided at all. *Registered resolution:* (a) *Scoring set:* the
   108 CGGA microarray patients that carry a reference label. The classifier is
   trained on RNA-seq, so this is a cross-platform evaluation and is reported as
   such. The 28 of these patients that also have RNA-seq are reported descriptively
   and are not part of the criterion: on 28 patients even a perfect 28 of 28 falls
   short (exact 95% lower bound 0.8766). (b) *Interval:* the exact binomial
   (Clopper-Pearson) two-sided 95% interval, whose lower bound is the decisive
   bound, replacing the percentile bootstrap for this metric only (section 3.1).
   (c) *Minimum n:* if fewer than 30 labelled patients are available when T1 is run,
   this criterion is NOT EVALUABLE rather than failed. A perfect score first clears
   at 29 patients (0.8806); 30 is registered as the floor: 30 of 30 gives 0.8843 and
   29 of 30 gives 0.8278. Under these bounds the floor is load-bearing: below 29
   patients even a perfect score would not clear. (d) *Harmonisation:* each
   sample's expression vector is z-scored across genes (for that sample alone,
   mean 0 and variance 1 over its genes), and the RNA-seq-trained classifier is
   applied to the z-scored vector. There is no shared reference and no
   cross-platform batch correction, and nothing in the rule depends on the
   external data beyond the sample being scored. This is a project choice, not
   the normalisation of Verhaak et al. (2010), who standardised each gene across
   patients within each platform. The rule was fixed, and entered in section 8,
   before the microarray expression matrix was downloaded. *The classifier (registered
   2026-11-27, before the deposit and before the matrix was downloaded):* single-sample
   nearest centroid, `gliostrata.profiling.centroid.classify_single_sample`. Genes: the
   831 genes of the Verhaak et al. (2010) ClaNC 840-gene signature that were resolved
   when the centroids were built. Centroids: the frozen set `verhaak840_tcga_all`
   (`configs/subtype_centroids/verhaak840_tcga_all.json`, SHA-256
   `c2b65d8ff0be2708bc8086e909f04a16a299823109186d674837c4d46f159bf7`): for each of the four
   classes, the mean of the per-sample z-scored log2(TPM + 1) profiles of the 246 TCGA
   RNA-seq patients carrying a published four-class label (Classical 64, Mesenchymal 82,
   Proneural 48, Neural 52; labels from the Brennan et al. 2013 table, which mixes
   Verhaak 2010 clustering labels with single-sample GSEA calls). Similarity: the
   Pearson correlation between the sample's z-scored vector and each centroid, over the
   genes present; the call is the highest, and an exact tie goes to the earlier class
   in the order Classical, Mesenchymal, Proneural, Neural. Nothing is tuned; no other
   candidate is evaluated; the three-class variant and the set estimated on the v2
   training split (`verhaak840_tcga_v2train`) are not used for T1. *Microarray input:*
   genes are matched by exact symbol to the set's symbols, with no alias resolution;
   rows sharing a symbol are averaged; values enter on the log2 scale, as distributed if
   the source's documentation says they are log2 values, otherwise as log2(x + 1), and
   that scale is decided from the documentation, and recorded, before any value is read.
   The number of genes used is reported. The code applies log2(TPM + 1) to RNA-seq; the
   microarray input rules are implemented in the T1 script, which is committed before
   the read. *Input scale and two format contingencies (amended 2026-11-27, after the OSF
   submission and before any microarray value was read):* CGGA's documentation says the
   intensities were "normalized using GeneSpring GX 11.0" and does not say log2, so the
   values enter as log2(x + 1) (`configs/t1_subtype/input_scale.json`). If any value is
   below 0, the values are already on a log scale, because linear intensities cannot be
   negative, and they are used as distributed; the output records this. A frozen gene
   with a missing value is dropped for that sample only, and the genes used per sample
   are reported. The script is `scripts/t1_subtype_cgga.py`; it reads only the
   `mRNA_array_301` rows of the cohort file, and it refuses the read until its synthetic
   dry run exists. (e) *Margin, stated so the choice of
   scoring set can be checked:* on 108 patients the criterion tolerates up to 6
   errors (102 of 108 gives 0.8830; 101 of 108 gives 0.8710). An observed accuracy
   near 0.95 needs roughly 71 to 100 patients to clear (68 of 71 gives 0.8814; 95 of
   100 gives 0.8872), so a smaller scoring set would be undecidable unless accuracy
   were near-perfect.
10. **T1 state-deconvolution r has no ground truth in the project.** *Registered
    resolution:* the truth is 100 pseudo-bulk mixtures simulated from the Neftel
    2019 single-cell atlas (GEO GSE131928; 28 glioblastoma tumours) with imposed
    Dirichlet random mixing weights; the imposed weights are what the estimated
    fractions are correlated with. For item 1, "patients" here means mixtures and the
    bootstrap resamples mixtures. Because the truth is simulated, this tests
    recovery of known proportions and is an optimistic test; it is reported as such.
    A CIBERSORTx run, if one is done, is a post-hoc sensitivity analysis and is
    labelled exploratory.
11. **T4 and T5 need a model that does not yet exist** (registered 2026-11-27). The
    simulator (L4) and the model that turns a configuration into a predicted survival
    curve (the response surface) are stubs at the time of deposit. *Registered
    resolution:*
    (a) *Stupp baseline.* For each patient, the predicted survival under standard therapy
    is the survival model's prediction with the treatment arm set to radiotherapy plus
    temozolomide, and no stem-cell configuration.
    (b) *Predicted RMST.* The restricted mean survival time to tau = 24 months, the area
    under the predicted survival curve from 0 to 24 months. The gain of a configuration is
    RMST(configuration) minus RMST(Stupp).
    (c) *T5.* A patient is a success if the gain of the configuration chosen for them is at
    least 3 months, and a harm if the gain is below 0.
    (d) *T4.* The optimisation objective is the predicted RMST at tau = 24 months, and
    regret (item 5) is measured in months of that RMST.
    (e) *Contingent specification.* T4 and T5 start only after the simulator, the response
    surface and the survival model they call are frozen and versioned, and entered in
    section 8 with their checksums, by a dated amendment made before any T4 or T5 read of
    CGGA. Such an amendment changes no threshold, bound or decision rule. It is
    confirmatory because it precedes the read, and it is reported as specified after the
    deposit.
    (f) *Closing rule.* If they are not frozen by 2027-06-30, T4 and T5 are closed as NOT
    EVALUABLE, the reason is stated, and H1 and H2 are reported as untested.

    *The model's design, decided 2026-11-27 (see section 8). Design A of
    `docs/t4t5_model_design.md` is adopted; (g) to (o) fix what that document left to the
    project lead. The numeric values themselves are frozen in the parameter file of (n),
    before any T4 or T5 read, under (e).*

    (g) *Survival model (L5).* `CoxPHModel` (Breslow) fitted on TCGA with covariates age,
    sex and treatment arm — **variant 2**, the model that decided T2's rows (section 2.3) —
    with its step baseline cumulative hazard H0(t) evaluated on [0, 24] months. Variant 3
    (plus IDH and MGMT) is **not** used, although its external C-index is higher (0.5915
    against 0.5385): it requires an IDH and an MGMT call for every patient, which the
    partition's eligibility rule (item 7) never required, so T5 would run on fewer than the
    218 patients its scoring set registers. On CGGA as a whole the equivalent complete-case
    drop is 283 to 241. Every one of the 218 has age, sex and arm. Variant 2 failed external
    validation (C-index 0.5385, 95% CI 0.4994 to 0.5807; Brier 0.2617); that failure is
    stated wherever a T4 or T5 number is reported, together with the proportional-hazards
    test of section 3.8.

    (h) *Simulator (L4).* The deterministic closed-form delivery module of Design A item 2,
    named **"L4-lite: closed-form delivery model, uncalibrated, literature-parameterised"**.
    This is what satisfies (e)'s requirement for a frozen simulator. It is **declared
    uncalibrated**, and that declaration stands whatever T3 returns: T3 data would calibrate
    f_k, a_k, g_k and n0 if they arrive by item 8's decision date, but T4 and T5 do not
    require them and are not deferred for them. No PDE or agent-based simulator is built
    under this registration.

    (i) *Efficacy prior.* Emax_p is fixed by a **rule declared before any read**, never by a
    value chosen with a threshold in view: Emax_p = 0.5 x ln(ratio_p), where ratio_p is the
    treated-to-control median-survival ratio of the nearest published preclinical study for
    payload p, and 0.5 is a translation factor on the log scale. The factor is **assumed**,
    and is declared as assumed in every table that reports a T4 or T5 number. The nearest
    study per payload, and the single value the rule yields, are recorded with their PMIDs in
    the parameter file of (n): CRAd from PMID 42179969 (ratios 1.33 to 1.48; the skeptical
    end, 1.33, is used), prodrug_enzyme from PMID 23658244, and CAR at the pessimistic end of
    its range because human CAR data in glioblastoma are uniformly poor (PMID 38216766, PMID
    40451950). *Consequence, stated in advance:* this rule places T5 below its 3-month bar by
    construction — the bar needs a hazard ratio near 0.59, and the rule yields about 0.82 to
    0.87 — so T5's expected verdict is a FAIL that reflects the state of the evidence, not a
    property of the framework. The break-even effect size is reported alongside it as an
    exploratory analysis (section 3.6).

    (j) *What "combo" combines.* The registered payload "combo" is defined **per carrier**,
    as the two payloads that carrier is registered as able to carry (item 4 excludes
    HB1_F3_line x CAR and MSC x CRAd): iPSC-NSC takes CRAd + prodrug_enzyme, MSC takes
    prodrug_enzyme + CAR, HB1_F3_line takes CRAd + prodrug_enzyme. iPSC-NSC is compatible
    with all three payloads, so its pairing is a declared choice. Payload effects combine by
    **Bliss independence on the hazard scale** (E_combo = E_a + E_b) with a synergy factor of
    1.0, assumed. A uniform definition was rejected because it would have placed an excluded
    payload on a carrier that item 4 registers as unable to carry it; an optimiser-chosen
    pair was rejected because two configurations sharing carrier, payload and site count
    could then be different therapies, which item 4's distinctness rule — and so T4's
    distinct-configuration count — depends on them not being. This changes no entry of item
    4: the 10 carrier x payload pairs, the 240 discrete combinations and the exclusion list
    all stand.

    (k) *Toxicity and delay.* The excess-log-hazard term of Design A item 4 is **in the
    primary surface**, not an exploratory variant: T_i(c) = tau_p (N / 5e7)^q (1 + eta (n-1))
    + delta_delay, with the parameters of that document's table (all assumed, each with a
    sensitivity range). It is the only source of harm in the model, so without it T5's harm
    fraction would be 0 by construction and its harm row could not fail; a safety row that
    cannot fail would be reported as vacuous rather than as met. A larger toxicity prior
    drawn from the CAR neurotoxicity rate (56% grade 3, PMID 40451950) was rejected as not
    transferable to carrier-cell delivery; it is a sensitivity range instead.

    (l) *Population.* All 218 patients of `cgga_external_v1`, with **no exclusion by IDH
    status**. Under WHO 2021 an IDH-mutant tumour is astrocytoma rather than glioblastoma,
    which is an argument for excluding these patients; it is not acted on, because the
    partition is sealed at 218 and changing who is in it after the seal is what sealing
    exists to prevent. Under (g) neither IDH nor MGMT appears anywhere in the frozen model —
    not in the baseline, and not as an effect modifier, because no sourced IDH or MGMT
    modifier exists for the configuration effect and one will not be invented. *The cost is
    disclosed:* IDH-mutant patients' Stupp baselines are therefore modelled as if they were
    IDH-wild-type. Subgroup analysis by IDH status is already registered as exploratory in
    section 3.6 and is where this is reported.

    (m) *Freeze date.* The amendment of (e) is entered by **2027-02-28**, four months ahead
    of (f)'s closing date. The slack is deliberate and is not to be spent on scope.

    (n) *Where the parameters live.* Every numeric parameter of the response surface, the
    delivery module and the toxicity term lives in a single versioned file,
    **`configs/response_surface_v1.yaml`**, whose SHA-256 is entered in section 8 at the
    freeze together with the survival model's. Each entry carries its value, its status
    (fitted on TCGA, published with its PMID, or assumed) and its sensitivity range.
    Constants are not kept in the module source, so that the registered checksum cannot be
    broken by an edit that changes no parameter.

    (o) *Controls, extending section 3.8.* Before the CGGA read, and reported with the
    result: (i) a *no-modifier control* — every patient modifier set to zero, which must
    yield a distinct-configuration count of 1, the optimiser's universal answer; (ii) a
    *null-payload control* — Emax = 0 for every payload, which must yield a gain of exactly 0
    for every patient and a harm equal to the toxicity term alone; (iii) the synthetic dry run
    of section 3.8; (iv) a run of the whole T4 and T5 pipeline on the held-out TCGA split,
    before any CGGA read, so that bugs are found without spending the registered read. A
    control that does not give its stated answer stops the read.

### 2.3 Outcome of T2, recorded after the external read (2026-11-26)

> This subsection is a post-registration record of a result.
> The registered specification it refers to is fixed at commit dc37e6c; the git history
> provides the timestamp. It is not part of the registered specification: it changes no
> threshold, no decisive bound and no decision rule. It was written after the CGGA
> external data had been read and is disclosed as such (section 8).

**T2 (external C-index >= 0.65): evaluated 2026-11-26, NOT MET.** Under the decision
rules in section 6, the C-index row of T2 is reported as **FAIL** for the model that
could be evaluated: its point estimate is below the threshold and its whole 95% interval
is below it. The committed four-feature model is **NOT EVALUABLE** on CGGA. This is one
of the eleven pre-registered thresholds in section 2; the other ten are not reported
here. T2's Brier row (18 months, <= 0.25) was evaluated later the same day, on a second read,
and is also reported as FAIL (below). T2's improvement-over-baseline row has not been
evaluated.

**What was measured.** The project's Cox proportional-hazards model was fitted on
TCGA-GBM and scored once on CGGA with Harrell's C-index. The 95% CI is a percentile
bootstrap over the evaluated CGGA patients (1000 resamples, seed 42) with the fitted
model held fixed. Internal C-index is 5-fold cross-validation on TCGA, quoted for
context only; it is not an external result.

| Variant | Features | TCGA fit n | CGGA n | Internal C-index | External C-index | 95% CI | Reported as |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1. Committed baseline | age, sex, KPS, treatment arm | 594 passed, 592 modelled | 0 | 0.6735 | -- | -- | NOT EVALUABLE. CGGA records no performance status, so KPS is missing for every patient. It was not imputed on CGGA. On TCGA the fit imputes: KPS set to the median, 80, for 150 of the 592, with an indicator (HR 1.57); treatment arm missing for 22, coded as its own level (HR 1.20). |
| 2. Transferable model | age, sex, treatment arm | 570 | 283 | 0.6520 | 0.5385 | 0.4994 to 0.5807 | FAIL |
| 3. Molecular model (secondary) | age, sex, treatment arm, IDH, MGMT | 298 | 241 | 0.6606 | 0.5915 | 0.5477 to 0.6341 | Below threshold; not directly comparable to variant 2 |

- Variant 3 was scored on a different, smaller sample of CGGA patients than variant 2
  (241 against 283), so it is not directly comparable and the 0.053 gap between the two
  is not a test of whether molecular features compensate for the loss of KPS. Its IDH
  and MGMT inputs are CGGA's own clinical annotations rather than this project's calls,
  so it is also a cross-platform comparison.
- An external C-index of 0.5385 whose interval contains 0.50 is consistent with
  discrimination no better than chance, and equally consistent with weak discrimination.
  This record commits to neither reading.

**Where the gap sits (descriptive).** Two differences are on record between three
C-indices, and neither carries an interval. The internal C-indices of the committed
baseline and of the transferable model differ by 0.0215 (0.6735 to 0.6520). That contrast
bundles three things: the KPS value; a missing-KPS indicator that is itself prognostic
(HR 1.57); and a different fit set (592 modelled patients against 570). It is not an
estimate of what KPS is worth. Moving the same three-feature model from TCGA to CGGA
differs by a further 0.1135 (0.6520 to 0.5385). Read descriptively the second gap is
about five times the first, which is consistent with the shortfall lying mainly in the
change of cohort rather than in KPS being unavailable in CGGA. That is not a tested
attribution. The mechanism is not investigated, and no cause is asserted here.

**What was not done.**

1. No tuning, feature selection or hyper-parameter search was performed against CGGA.
2. The three variants in the table were pre-specified before the first CGGA evaluation,
   and all three are reported. None was selected after the fact.
3. There were no re-runs after the first evaluation.
4. The mechanism of the external-validation failure was not investigated, because doing
   so would require additional CGGA reads and would contaminate the external partition.

**Pre-registration integrity.**

1. CGGA has been read twice, both for registered rows of T2: the C-index row, and then,
   on a second read, the Brier row. Before the first read it had only been characterised
   at schema level (row counts, missingness, value ranges), with no feature related to any
   outcome. This is the constraint any pre-registered protocol operates under.
2. Any analysis of CGGA beyond the registered tests (T1, T4 and T5) would require a new,
   separately pre-registered exploratory analysis.
3. T1, T3, T4 and T5 are unchanged by this amendment.

**Relation to the registered expectations.** Section 5 registered, in advance, an external
survival C-index below 0.65 as the most likely first-pass outcome, and section 7
pre-committed the response to it, which applies here. The threshold is reported as FAIL
and is not re-stated. Of what that fallback calls for, only the achieved C-index with its
CI is reported; the calibration and decision-curve analyses it also names (sections 3.4
and 3.5) have not been run.

**Deviations disclosed.**

- No sealed `cgga_external_v1` manifest existed when T2 was run (section 4.5 lists it as
  not yet created), so the evaluation used the ingested CGGA cohort as a whole rather
  than a frozen partition as section 3.3 specifies.
- CGGA's 333 rows are 305 unique patients; each patient was scored once. A patient
  missing a feature of the variant being scored, or the outcome, was excluded rather
  than imputed, which gives the CGGA n in the table.

The full account is `docs/cgga_validation.md`; the result file is
`data/results/cgga/cox_external_validation.json`, regenerated by
`scripts/validate_cox_cgga.py`.

**T2 (Brier score at 18 months <= 0.25): evaluated 2026-11-26, NOT MET.** The registered
rule for this row is PASS only if the point estimate and the upper bound of the 95% interval
are both <= 0.25. The transferable model (age, sex, treatment arm) scores 0.2617 (95% CI
0.2391 to 0.2829) on the same 283 CGGA patients, so the row is reported as **FAIL**. The
committed baseline is again **NOT EVALUABLE** (no KPS in CGGA). The molecular model scores
0.2472 (95% CI 0.2226 to 0.2729) on 241 patients: its point estimate is under the threshold
and its upper bound is over it, which is INCONCLUSIVE by the registered rule; it is
secondary and does not decide the row. The TCGA fits are the ones the C-index row was
computed with (internal C-index identical to six decimals).

**What was measured.** The inverse-probability-of-censoring-weighted Brier score at 18
months, from the TCGA-fitted model's predicted survival probability for each CGGA patient.
The 95% CI is a percentile bootstrap over the evaluated CGGA patients (1000 resamples,
seed 42) with the model and its predictions held fixed. Brier scores at 12 and 24 months and
a null-model reference are reported descriptively in `docs/cgga_validation.md` section 12.

**What this is.** A second read of CGGA, made because this registered row had not been
evaluated. Section 3.3 allows one read per registered analysis. The rules for this row were
fixed in `scripts/validate_brier_cgga.py` and its tests before that read; the rules for the
C-index row and the first script were not changed. No threshold, no decisive bound and no
decision rule changed. The read ran, like the first, without a sealed `cgga_external_v1`
manifest. Calibration slope and intercept and decision-curve analysis (section 3.4) were not
done. With both rows evaluated, T2's C-index and Brier criteria are reported as FAIL for the
model that can be evaluated, and its improvement-over-baseline row remains unevaluated.
Full account: `docs/cgga_validation.md` section 12; result file
`data/results/cgga/brier18_external.json`.

**Decisions on T2's open items (2026-11-27, before the deposit).** Each follows from text
already registered, and none changes a threshold, bound or rule.
(1) *Evaluation unit:* the patient, one row per patient (305, not 333 rows), as the
"Sample size" section states for every test.
(2) *Status:* the C-index and Brier rows are FAIL for the model that can be evaluated.
Section 7's fallback applies: the contribution is re-framed as benchmark establishment.
(3) *The two reads:* they are the reads of record for those two rows. T2 is not re-run on
`cgga_external_v1`, and the deviation (no sealed partition at the time) stays disclosed
above.
(4) *KPS:* the committed four-feature model is NOT EVALUABLE on CGGA, which records no
performance status. No other cohort is registered for it.
(5) *Calibration slope and intercept, calibration plots and decision-curve analysis
(sections 3.4 and 3.5):* not done for T2's evaluation, as decided in
`docs/cgga_validation.md` section 12.5. They decide no registered row, and each would be a
further read of CGGA.
The improvement row is the one T2 read still to come (section 2.2 item 2 (j)).

**Anticipated but not pre-registered.** The T2 result leaves open questions: whether the
cross-cohort shift of about 0.113 reflects population, treatment era, cohort composition
or measurement differences; whether a model whose primary features are molecular would
transfer better; and which cohort could test the committed four-feature model at all,
given that CGGA cannot supply KPS. These are questions, not hypotheses, and none is
registered here. Any such analysis would need its own pre-registration, its own stated
hypothesis, and its own decision about whether the remaining CGGA reads are worth
spending. It is not executed here and its results are not anticipated by this document.

### 2.4 Outcome of T1's state-deconvolution row, recorded before deposit (2026-11-27)

> This subsection is a record of a result, not part of the registered specification. It
> changes no threshold, no decisive bound and no decision rule. The row was run before
> this registration was deposited, so the result was known at the deposit and the row is
> disclosed as run before deposit (section 8). It reads no cohort.

**T1 (state deconvolution Pearson r >= 0.75): evaluated 2026-11-27, PASS, run before
deposit.** Under the decision rules in section 6 the row is reported as PASS: the lower
bound of the 95% interval of the registered metric clears the threshold. This is one of
the eleven pre-registered thresholds in section 2. T1's subtype-accuracy row has not been
evaluated; T1 as a whole is not reported as met.

**What was measured.** The registered truth (section 2.2 item 10): 100 pseudo-bulk
mixtures simulated from the Neftel 2019 atlas (GEO GSE131928) with imposed Dirichlet
mixing weights. Everything the registration left open was fixed in
`docs/t1_deconvolution_spec.md` and committed (`c68f556`) before the run, with the script
`scripts/t1_deconvolution_mixtures.py` and its tests: 6,568 malignant Smart-seq2 cells
(26 tumours) whose per-cell malignant flag and six meta-module scores come from the
Curated Cancer Cell Atlas release of the same data, each assigned its highest state
(NPC, OPC, AC, MES, with the two NPC and the two MES sub-modules merged by maximum);
flat Dirichlet weights; 500 cells per mixture, drawn by state; seed 131928. The
estimator is the project's own marker-enrichment scorer
(`gliostrata.profiling.deconvolve.deconvolve_states_expression`, defaults), applied once
to the 100 mixtures; it is not CIBERSORTx. The metric is the mean of the four per-state
Pearson r between estimate and imposed weight (item 1); the interval is a percentile
bootstrap over mixtures, 1000 resamples, seed 42, with the estimates held fixed. The
verdict was computed by the registry (`T1_state_deconvolution_pearson_r`).

| Metric | Value | 95% CI |
| --- | --- | --- |
| Mean of the four per-state r (decides) | **0.9521** | **0.9430 to 0.9589** |
| NPC-like | 0.9647 | 0.9472 to 0.9768 |
| OPC-like | 0.9528 | 0.9290 to 0.9682 |
| AC-like | 0.9204 | 0.8885 to 0.9447 |
| MES-like | 0.9705 | 0.9574 to 0.9797 |
| Pooled r over all 400 pairs (descriptive only) | 0.9506 | |

No interval was clipped and no bootstrap resample was dropped. One secondary analysis
was declared in the specification before the run and is descriptive only: the same
estimator with its cohort-independent normalisation (`fraction_method="proportional"`)
gives 0.8541 (95% CI 0.8358 to 0.8750).

**How to read it.** As item 10 registers, this tests recovery of known proportions and is
optimistic. Specifically: the mixtures contain malignant cells only, with no immune,
stromal or normal-brain cells; the truth is computed from the same cells; and the atlas
states and the estimator's marker panels are both drawn from the Neftel meta-modules, so
part of the agreement is built in. The result says the estimator recovers simulated state
proportions under these conditions. It does not say how well it estimates state
composition in a TCGA or CGGA tumour.

**What was not done.** One run, from the committed specification, with no re-run, no
change of seed, cell count or estimator after the result, and no CIBERSORTx run. Full
account: `docs/project_record/06_evidence_reports.md`, section "T1, state-deconvolution
row"; result file `data/results/t1/neftel_mixture_deconvolution.json`.

---

## 3. Analysis plan

### 3.1 Uncertainty quantification

Every reported metric carries a **bootstrap 95% confidence interval from 1000 resamples**.
The resampling unit is the **patient** (not the sample, row, or state-patient pair), so
that multi-sample patients cannot inflate precision. Bootstrap seeds are fixed and
recorded. Percentile intervals are used; where a metric is bounded (accuracy, fractions,
correlations) the interval is clipped to the metric's natural range and the clipping is
reported. One exception is registered in advance: T1 subtype accuracy uses the exact
binomial (Clopper-Pearson) interval instead (section 2.2, item 9), because a percentile
bootstrap over a sample in which every patient is classified correctly is degenerate
(lower = point = upper = 1.00) and would clear any threshold.

### 3.2 Multiple-comparison correction

The eleven pre-registered thresholds in §2 constitute the confirmatory family.
**Benjamini–Hochberg** correction at FDR = 0.05 is applied across the p-values of the
family where a p-value exists (T4's Wilcoxon test, and any threshold tested via a formal
hypothesis test). Threshold comparisons that are CI-based rather than p-value-based are
not re-corrected — the 95% CI requirement is itself the pre-registered decision rule — but
the number of thresholds evaluated is always reported alongside, so no reader has to infer
the family size. Any subgroup or sensitivity analysis is labelled exploratory and is
corrected separately; exploratory results never change a confirmatory verdict.

### 3.3 External validation

Final numbers for T2 are reported on the **frozen external test partitions** of CGGA
(primary external cohort), never on
TCGA-derived data. Final numbers for T1 are reported on the scoring sets that section 2.2
items 9 and 10 register, not on those partitions: for subtype accuracy, the CGGA
microarray patients that carry a reference Verhaak label, fixed by rule and frozen as
their own manifest before the read; for state deconvolution, the 100 pseudo-bulk mixtures
simulated from the Neftel 2019 atlas. Neither is TCGA-derived. Internal cross-validation
numbers may be reported for context and are always labelled as internal. Model selection, hyper-parameter tuning, feature selection,
threshold calibration and architecture choices use **only** the development cohort's
training split. The external partitions are read exactly once per registered analysis, at
the point of confirmatory evaluation.

### 3.4 Calibration, not just discrimination

Discrimination (C-index) is reported alongside: Brier score at 12, 18 and 24 months (18
months is the pre-registered threshold); calibration plots of predicted vs observed
survival in deciles of predicted risk; calibration slope and intercept; and the
Hosmer–Lemeshow style goodness-of-fit test adapted for censored survival data. A model
that discriminates but is miscalibrated is reported as such, because a miscalibrated model
cannot support the RMST-gain claims in T5.

*As applied to T2 (2026-11-27):* T2 was evaluated before the deposit, and the analyses of
this section were decided not to be run on CGGA for this thesis (section 2.3, "Decisions on
T2's open items", item 5). They apply to any survival model evaluated under a future
registration.

### 3.5 Decision-curve analysis

Net benefit is reported across the clinically plausible range of threshold probabilities
for the survival model, against both treat-all and treat-none strategies, so that the
result is expressed as clinical utility rather than only as a discrimination statistic.
*As applied to T2 (2026-11-27):* not run on CGGA for this thesis, as for section 3.4.

### 3.6 Sensitivity analyses (pre-specified, labelled exploratory)

*Scope (registered 2026-11-27).* These analyses run only inside the reads already
registered: T1, T2's improvement row, T4 and T5. They use the same data, are labelled
exploratory, and never change a verdict. An analysis that would need a further read of
CGGA is not run. T2's C-index and Brier rows, already evaluated, get none.

1. Missing-data imputation: complete-case vs multiple imputation, reported side by side.
2. Batch-effect correction: with and without cross-cohort harmonisation.
3. Subgroup analysis by Verhaak subtype, MGMT methylation status and IDH status.
4. Seed sensitivity: the full confirmatory analysis re-run under the partition seed
   recorded in each manifest, plus two alternative seeds reported as exploratory only. The
   confirmatory verdict is the one from the registered seed; the alternatives exist to
   show the result is not seed-fragile, and a large spread across seeds is reported as a
   limitation rather than suppressed.

### 3.7 Reproducibility release

Code (GitHub plus a Zenodo DOI), container images, trained model weights, the feature
store, and every partition manifest are released with the paper. The threshold-table
checksum and the partition manifests are sufficient for a third party to confirm that the
reported thresholds and partitions are the registered ones.

### 3.8 Model checks and controls (registered 2026-11-27)

None of these changes a model, an analysis or a verdict. Each runs inside a read already
registered, so none adds an access to CGGA.

1. **Proportional hazards.** Every Cox model fitted for a registered read is reported with
   the Grambsch-Therneau test on scaled Schoenfeld residuals, per term and globally.
   Event times are rank-transformed, as in R's `survival::cox.zph` (survival 2.x). The
   implementation is `gliostrata.prediction.ph_test`; it matched lifelines to six decimals
   on a reference case. A violation is reported as a limitation.
   *Already known at registration:* the TCGA fit of variant 2, the model on our side of
   T2's improvement row, violates proportional hazards in age (p = 0.0002; global
   p = 0.0006). This was found when the synthetic dry run of item 3 was first run, before
   the deposit. It changes nothing in the registered analysis.
2. **Negative control.** In each read, each risk score's C-index is also computed against
   outcomes permuted within CGGA batch: 1000 permutations, seed 42, with time and event
   kept together. The mean and 95% range are reported; about 0.5 is expected.
3. **Synthetic dry run.** Before each real read, the exact evaluation script runs end to
   end on synthetic files in the formats of the real inputs, and its output is committed.
   T2's improvement-row script refuses the real read until the dry-run output exists.
4. **Positive controls already in place.** These are not new rows:
   - T1's deconvolution row recovers proportions that are known by construction;
   - T2's internal sanity check reproduces the committed C-index of variant 1 (0.6735).

### 3.9 Blinding and procedural safeguards (registered 2026-11-27)

The project has one analyst, so no analysis is run blind to the design or the outcomes. The
safeguards are procedural:
- every design choice is made without reference to the external outcomes;
- every rule and every script for a read is committed before that read;
- external test partitions are sealed, and a manifest is re-derived from (seed, patient
  id) in CI;
- the read guard refuses any value read that is not registered, deposited and declared,
  and records every access, of any kind, in an append-only ledger;
- each verdict is produced by the registered threshold function, not by judgement;
- the synthetic dry run of section 3.8 runs before each read.

### 3.10 Outliers (registered 2026-11-27)

No observation is removed, trimmed or winsorised for being extreme. Values are used as
recorded at source. Known implausible values are disclosed rather than corrected: three
TCGA patients have KPS 0 at source. KPS enters only the committed four-feature model,
which cannot be scored on CGGA.

---

## 4. Frozen partition procedure

### 4.1 Why hashing, not shuffling

A shuffle-based split (`train_test_split`, `df.sample`) is a function of the entire table:
the random stream is consumed row by row, so a patient's assignment depends on how many
rows precede it and how many rows exist. Change an upstream QC filter, add a cohort, sort
differently, or upgrade pandas/NumPy, and patients move between train and test. Once a
patient has moved from train to test, every model already fitted has seen them, and the
"external" result for that patient is worthless — permanently, because the leak has
already entered every hyper-parameter and model-selection decision.

Partition assignment is therefore a pure function of `(seed, patient_id)`:

```
u = sha256("{seed}|{patient_id}") mapped into [0, 1)
test  if u < test_fraction  else  train
```

Measured on the first TCGA-GBM clinical cohort (the 583 patients frozen as
`tcga_gbm_clinical_v1`; the cohort now in force is 594 patients, see section 4.5),
seed 20251126, test_fraction 0.30): shuffling the input rows four different ways
reproduced the identical 419/164 split; adding 50 unrelated patients moved **0** of the
original patients; and removing 50 patients moved **0** of the remainder. The same
perturbation applied to a `df.sample(frac=0.30)` split moved **173 of them**. That
contrast is the reason this procedure is registered.

Stratification, when used, salts the hash with the stratum label
(`sha256("{seed}|{stratum}|{patient_id}")`) rather than ranking patients within each
stratum. Ranking would give exact per-stratum quotas but would re-introduce dependence on
the patient set, which is the property being eliminated. The registered consequence:
per-stratum test fractions are **approximately**, not exactly, equal to the target
(empirically within ±0.05 at n ≈ 1000 per stratum), and achieved per-stratum fractions are
recorded in each manifest.

Duplicate patient ids are collapsed before assignment, so a multi-sample patient can never
straddle the split.

### 4.2 Manifest contents

Each partition is frozen as a JSON manifest (`data/partitions/*.json`, under version
control) recording: schema version; partition name; seed; target test fraction; patient-id
column; stratification column; the exact patient-id list for each split; a SHA-256 over
each sorted id list and over the union; row counts per split; the cohort source path and a
SHA-256 over **the source file's bytes**; the source row count; achieved per-stratum
balance; an ISO-8601 UTC creation timestamp; a free-text provenance note; and the seal flag.

The id-list hashes let a reviewer confirm the splits without reading the id lists. The source
content hash catches the case the id list cannot see: a cohort file rewritten with the same
patients but altered column values (an edited survival time, a re-imputed covariate).

### 4.3 Verification

`verify_partition` (CLI: `gliostrata freeze-partition --manifest M.json --verify`)
recomputes the partition from the manifest's own recorded parameters and compares it
patient by patient against the frozen lists. It reports precisely what drifted:

* **ids added** — present in the cohort, absent from the frozen partition;
* **ids removed** — in the frozen partition, missing from the cohort;
* **ids moved** — assigned to a different split than the one frozen, with the direction;
  a train→test move is flagged as `CONTAMINATION RISK` because those patients have already
  been seen by every model fitted on this partition;
* **source content hash mismatch** — same patients, different file bytes;
* **row-count change** — reported for information, and explicitly *not* treated as drift,
  because the hash split is row-count independent. This is the case that would have
  silently broken a shuffle-based split.

Verification is run before every evaluation run. A drifted partition blocks evaluation
(the CLI exits non-zero) until the drift is explained in the amendment log.

### 4.4 Sealing

A manifest written with the seal set refuses to be overwritten. Re-running the freeze
against a sealed manifest raises `PartitionSealedError` with an explanation of the
scientific consequence, and nothing is written. Overwriting requires an explicit `force`
flag, and is legitimate **only** as a documented amendment (§8) in which the old manifest
is archived alongside the new one and every model is retrained from scratch. Re-freezing
is never a routine re-run.

An external cohort has no training data, so its seal holds out every patient: a test
fraction of exactly 1.0 (CLI `--all-test`), an empty training list and a test list of
every patient id. It is produced by the same hash rule as any other partition.

### 4.5 Partitions registered at the time of this pre-registration

| Partition | Cohort | n | Split | Seed | Role | Sealed |
| --- | --- | --- | --- | --- | --- | --- |
| `tcga_gbm_clinical_v1` | TCGA-GBM clinical-only | 583 | 419 train / 164 test | 20251126 | **Development.** Internal held-out set for model-selection discipline. **Not** a pre-registered external test. Superseded by `tcga_gbm_clinical_v2`; retained unmodified as the historical record. | Yes (2026-11-26) |
| `tcga_gbm_clinical_v2` | TCGA-GBM clinical | 594 | 425 train / 169 test | 20251126 | **Development.** Supersedes v1 after the cohort was re-ingested (section 5). Unstratified. | Yes (2026-11-26) |
| `cgga_external_v1` | CGGA | 218 | 0 train / 218 test (all-test) | 20251126 (changes no assignment when every patient is held out) | **Primary external test for T2, T4 and T5.** It is not the scoring set of either T1 row (see the note below the table). | **Yes (2026-11-27).** Sealed after the T2 evaluation and before any T4 or T5 computation, as item 7 registers, over the 218 patients eligible by its rule (133 in mRNAseq_693, 85 in mRNAseq_325). Manifest `data/partitions/cgga_external_v1.json`, SHA-256 `040ba322d77f1f2cc205993d26c2aab63f28de4c86a8f2cc6bd3c056e8c3829e`; it records the SHA-256 of the identifiers-only eligibility table (`43da219d917870d9d169f20a6d0a474bc59f27fe78e6ddfc6a84cb1fd25ebfad`), which is git-ignored. The eligibility step was an availability access, recorded in the access ledger (entry 11). T2 was read twice before the seal, on the ingested cohort as a whole (commits 56cc36b and the Brier row the same day): a disclosed deviation from section 3.3; see section 2.3. |
| `glass_external_v1` | GLASS | — | — | — | Longitudinal external test | **Not created, and not used under this registration** (2026-11-27). T2 is reported as already conducted on CGGA and is not re-run, so no registered verdict needs GLASS. Any analysis of GLASS would need its own pre-registration. |

T1 is scored outside this partition, by registered rule. Its subtype-accuracy row uses the
CGGA microarray patients that carry a reference Verhaak label (108 at ingest; section 2.2
item 9), fixed by the rule `T1_subtype_accuracy_cgga` in `configs/scoring_rules.yaml` and
frozen as its own manifest before the read. Of the 108, 28 also have RNA-seq and may also
be in `cgga_external_v1`; the overlap does not change either scoring set. T1's
deconvolution row is scored on 100 pseudo-bulk mixtures simulated from the Neftel 2019
atlas (section 2.2 item 10) and reads no cohort.

Both TCGA partitions are unstratified: the `subtype` column of the clinical cohort file
is entirely null (see §5). When molecular data lands, the stratified partition must be
created as a **new, separately named** partition, never by overwriting either of these.

---

## 5. What is already known and observed at the time of registration

This section exists so that no reader has to take on trust that the test data was untouched.

**Data actually in hand.** TCGA-GBM: `data/raw/tcga_gbm/clinical_cohort.parquet`, 594
eligible patients from 617 case records, **clinical variables** — patient id, age, sex, KPS, overall
survival in months, vital status and treatment arm; 284 bulk RNA-seq samples (STAR read
counts and TPM); and MGMT promoter methylation and IDH mutation calls from the GDC, held
in a separate joined table for 368 patients (IDH), 415 (MGMT) and 314 (both). In the
clinical file itself the columns for extent of resection, MGMT, IDH, Verhaak subtype and
the Neftel and TME fractions are **entirely null**. No proteomic or imaging data has been
ingested. One single-cell dataset is in hand, obtained 2026-11-27 for T1's deconvolution
row: the Neftel 2019 atlas, as the Smart-seq2 expression matrix of GEO GSE131928 (sample
GSM3828672; 7,930 cells from 28 tumours) and the per-cell annotation of the Curated
Cancer Cell Atlas release of the same data (`data/provenance/neftel2019_sources.json`).
The second cohort is CGGA: 305 unique patients in 333 rows, characterised at schema level
and then read twice for T2 (section 2.3). Since then it has been accessed only without
reading values, and every access is recorded in `data/access_ledger.jsonl`:
- the seal's availability check (entry 11);
- the download of CGGA's RSEM FPKM release for T2's comparator, which is in hand and has
  been opened for gene names and sample ids only (entries 13 and 16);
- the gene names of the counts matrix (entry 17).

**Data not yet obtained.** CGGA has been downloaded and read twice for the T2 evaluation
(C-index row, then Brier row); no reads beyond the registered T1, T4 and T5, and the one
read for T2's improvement row registered in section 2.2 item 2 (j), will occur without a new
pre-registration. GLASS has **not** been
downloaded. CPTAC-GBM, any single-cell reference atlas other than the Neftel 2019 atlas,
TCIA imaging, and the PDX growth curves and in vivo NSC tracking data required for T3
have **not** been obtained. The CGGA microarray expression matrix was fetched on 2026-11-27 and stored unopened (ledger entry 21), and T1's 108-patient scoring set is frozen (`cgga_t1_subtype_v1`). The
external partition `cgga_external_v1` was sealed on 2026-11-27 (218 patients, all-test;
§4.5); GLASS is not used under this registration (§3.3, §4.5).

**Analyses already run.** The repository contains a working but *synthetic-data* pipeline:
the L4 simulation layer and the L6 optimiser's objective are stubs that operate on
generated data, exercised by the test suite and a dry-run CLI command. Those runs establish
that the plumbing executes and produce no result about glioblastoma. One registered row
has been evaluated on the Neftel 2019 atlas before this deposit: T1's state-deconvolution
row, on 100 simulated mixtures (section 2.4; PASS, 0.9521, 95% CI 0.9430 to 0.9589).
Analyses on real TCGA-GBM data have also been run; none is an external result, and no number in this
pre-registration is derived from them apart from those in section 2.3 and the centroid-set
counts in section 2.2 item 9(d): descriptive
statistics of the clinical cohort (cohort size, age and survival distributions, censoring
rate); the IDH sensitivity analyses (`scripts/idh_sensitivity_cox.py`,
`scripts/idh_sensitivity_subtype.py`); a concordance check of the deconvolution module's
marker panels on the expression matrix; the T2 fit itself, whose internal C-index is
cross-validated on TCGA; and, in stage 3b (2026-11-27), four more: the estimation of two
frozen subtype centroid sets from TCGA RNA-seq and published labels
(`verhaak840_tcga_v2train` on the 170 labelled v2 training patients; `verhaak840_tcga_all`
on all 246 labelled patients); single-sample subtype calls for the 283 profiled patients,
with a descriptive agreement check against the published labels of the 76 labelled v2-test
patients (56 of 76, 0.737; most of those labels come from another classifier, so this is
not a validation); a measurement of how far cohort-relative subtype and state scores move
when the scoring set changes; and the build of the feature store over all 594 patients,
whose state and TME scores are cohort-relative over the 283 profiled patients, 82 v2-test
patients among them. One more TCGA analysis followed, for option B of T2's comparator
(section 2.2 item 2 (j)): 23 published gene sets, each refitted and scored by 5-fold
cross-validation on the 200 v2-training patients with RNA-seq, then the chosen set
refitted on 282 patients. Then the synthetic dry run of section 3.8, which reads no CGGA
file. Its proportional-hazards tests ran on the two TCGA fits: variant 2 violates
proportional hazards in age (p = 0.0002), and option B's refit does not (global
p = 0.41). The
partition demonstration quoted in §4.1 used patient identifiers and split assignments only
— no outcome variable was examined in connection with any split.

**Models fitted, and on which patients.** No model has been fitted on CGGA: the T2 Cox
variants were fitted on TCGA and only scored on CGGA, once (section 2.3). The development
split's test set has not been held out from those fits. The T2 script fits variants 2 and 3
on every complete-case TCGA patient, and variant 1 on the 592 with an observed outcome
(KPS imputed, see section 2.3), and reports a 5-fold cross-validated internal C-index for
each, so
the 169 patients in the test split of `tcga_gbm_clinical_v2` were in the fit set. The
subtype centroid set reserved for T1 was estimated on all 246 labelled TCGA patients, the
v2 test split included (section 2.2 item 9(d)). The sealed external partition
`cgga_external_v1` was created on 2026-11-27, after the T2 evaluation and before any T4
or T5 computation (section 2.2 item 7; section 4.5). This pre-registration is being written before any
confirmatory evaluation other than T2 and T1's state-deconvolution row (section 2.4).

**Prior expectations.** The thresholds in §2 are taken from published external benchmarks,
not from any GlioStrata result. The protocol's own risk register judges a first-pass
external survival C-index of ≈0.60 — below the 0.65 threshold — to be the single most
likely outcome (§6). That expectation is registered here in advance, so that reporting
0.60 later is a pre-registered outcome rather than a disappointment to be re-framed.

**TO BE COMPLETED BY AUTHOR: confirm, at the moment of OSF deposit, that the above remains
accurate, and amend it if any further data has been obtained or any model trained in the
interim.**

---

## 6. Decision rules: PASS, FAIL, INCONCLUSIVE, NOT EVALUABLE

| Verdict | Condition | How it is reported |
| --- | --- | --- |
| **PASS** | The decisive 95% CI bound clears the threshold, and every registered side condition (T5's point-estimate limit, T4's Wilcoxon test) also holds. | Threshold met. |
| **INCONCLUSIVE** | The point estimate is on the right side of the threshold but the decisive CI bound does not clear it; or a required auxiliary test was not reported. | Threshold **not** met; reported as an imprecise result, with the CI and the sample size that produced it. Never reported as a pass, and never as a trend toward one. |
| **FAIL** | The point estimate itself misses the threshold; or a registered side condition is definitively violated (harm point estimate > 10%, Wilcoxon p ≥ 0.01). | Threshold not met; achieved value reported in full with its CI. |
| **NOT EVALUABLE** | The data required to compute the metric could not be obtained, or a comparator could not be re-run on the frozen partition (see §2.2 items 2, 3, 8, 9). | Reported as not evaluable, with the specific missing input named. Never pooled with PASS. |

Rules that bind regardless of outcome:

1. Every one of the eleven thresholds is reported, including every miss. A threshold is
   never dropped from the results table because it failed.
2. A FAIL or INCONCLUSIVE verdict does not license switching metrics, cohorts, bootstrap
   procedures or CI methods. Any such change after the test data is read is an amendment
   (§8) and is disclosed as a post-hoc deviation in the paper.
3. The total number of confirmatory thresholds (eleven) is stated wherever any subset is
   reported.
4. Interpretive precedence: an H1 pass with an H2 failure is reported as a single
   configuration beating Stupp, not as personalisation (§1.4). A T5 pass built on a failed
   T2 is reported with that dependency stated in the abstract.
5. Exploratory analyses are labelled exploratory in every table and figure caption.

---

## 7. Pre-committed fallbacks (Protocol risk register, Table 9-3)

These fallbacks are committed **now**, before any result is seen, so that a threshold miss
is handled by a registered plan rather than by improvisation.

| Risk | Likelihood | Impact | Pre-committed fallback |
| --- | --- | --- | --- |
| External survival C-index < 0.65 | High | Medium | Report the achieved C-index with full defensibility (CI, calibration, decision-curve analysis) and **re-frame the contribution as benchmark establishment** for external GBM survival prediction on frozen public partitions, rather than as a threshold-beating result. The threshold is reported as FAIL, not re-stated. |
| End-to-end success threshold (T5) missed | Medium | High | The **pre-registered negative result is the publication**: a framework that does not find a personalised configuration beating Stupp on public data is an informative result about the data and the search space, and is reported as such. No post-hoc redefinition of "success" is permitted. |
| Wet-lab collaborator not identified by M3 | Medium | Low | The computational resource paper stands alone; T6/T7 are dropped from the thesis scope and the omission is stated. |
| dbGaP access takes longer than 8 weeks | Medium | Low | Use the open-access TCGA subset only, and state the resulting reduction in cohort size and feature coverage as a limitation. |
| Reviewer demands an additional validation cohort | Medium | Medium | Add CPTAC-GBM as a third external test, as a **new** sealed partition and a disclosed post-registration addition — never by re-freezing an existing partition. |
| Required T3 calibration data unobtainable | Medium | Medium | Report T3 as NOT EVALUABLE, present the simulation layer as uncalibrated, and withdraw any claim that depends on its quantitative accuracy. |

---

## 8. Amendment log

Any change to §2 (thresholds), §3 (analysis plan), §4 (partitions) or §6 (decision rules)
after deposit is recorded here with its date, its reason, and whether it was made before or
after the external test data was read. An amendment made after reading test data is
disclosed as a post-hoc deviation in the paper. An empty log is a claim, and the
threshold-table checksum plus the partition manifests are what make it a checkable one.

| Date | Section | Change | Reason | Test data already read? |
| --- | --- | --- | --- | --- |
| 2026-11-26 | 2.3 (added) | New subsection 2.3 added to record the outcome of protocol test T2 (external C-index >= 0.65 on CGGA), which was evaluated and not met. It changes no threshold, no decisive bound, no decision rule, and no other test's status. Full account: `docs/cgga_validation.md`; result file: `data/results/cgga/cox_external_validation.json`. | Records a result rather than a plan. The committed model is NOT EVALUABLE on CGGA (no KPS); the transferable model scored an external C-index of 0.5385 (95% CI 0.4994 to 0.5807). | **Yes.** Written after the CGGA external read and disclosed as such. The evaluation ran without a sealed `cgga_external_v1` manifest (section 4.5). |
| 2026-11-26 | 2.2 items 3, 5, 7 (completed), 2.2 items 9 and 10 (added), 3.1 (one exception), 4.5 (`cgga_external_v1` row), 6 (cross-reference) | Completed the two items left open at draft stage: the T3 growth-curve rule and dates (item 3) and the T4 evaluation budget and tie rule (item 5). Registered the eligibility and sealing order of the external partition (item 7). Added the T1 subtype-accuracy guard rails (item 9: scoring set, exact interval, minimum n) and the T1 deconvolution ground truth (item 10). No threshold, no decisive bound and no metric changed. | The protocol left these under-specified, and each is a choice that could otherwise be made after seeing the data. The microarray harmonisation method (item 9d) is still open. | **Yes.** CGGA had been ingested and read once for T2 before this entry. None of T1, T3, T4 or T5 has been evaluated, and the counts used here (108 labelled patients, 28 with RNA-seq, 218 eligible) come from the schema-level ingest report, not from any model result. |
| 2026-11-26 | 2.2 item 9(d) completed; item 9(c) one sentence added | Fixed the T1 cross-platform harmonisation rule: each sample's expression vector is z-scored across genes and the RNA-seq-trained classifier applied to it, with no shared reference and no batch correction. Added to item 9(c) that the floor of 30 is load-bearing under two-sided bounds. | The rule has to be fixed before the matrix it will be applied to is seen. It is a project choice, not the normalisation of Verhaak et al. (2010). | **No.** The microarray expression matrix has not been downloaded. |
| 2026-11-26 | 4.1, 4.2, 4.4, 4.5, 5 (corrected); 9 (row added) | Reconciled passages written for the first TCGA cohort (the superseded v1 partition, and the claim that only clinical data was in hand) with the data in hand: 594 patients, `tcga_gbm_clinical_v2` (425/169), 284 RNA-seq samples, joined IDH and MGMT calls. Corrected section 5's claim that no model had been fitted on any test split: the T2 variants were fitted on all complete-case TCGA patients, which include the v2 test split. Recorded the all-test seal in 4.4 and in the section 9 table. | A registration attachment must not contradict the OSF form or the sealed manifests. | **Yes.** Written after the T2 read; changes no threshold, bound or rule. |
| 2026-11-26 | 2.3 (Brier result added); 2.2 item 7, 2.3, 4.5 and 5 ("read once" corrected to two reads) | Recorded the outcome of T2's Brier row (Brier score at 18 months <= 0.25 on CGGA), evaluated on a second read and not met: the transferable model scores 0.2617 (95% CI 0.2391 to 0.2829). The rules for the row (IPCW score, the same patients and models as the C-index row, percentile bootstrap, verdict from the registry) were fixed in `scripts/validate_brier_cgga.py` and its tests before that read. Changes no threshold, no decisive bound, no decision rule and no other test's status. Full account: `docs/cgga_validation.md` section 12; result file: `data/results/cgga/brier18_external.json`. | Records a result rather than a plan. | **Yes.** Written after the second CGGA read and disclosed as such. That read, like the first, ran without a sealed `cgga_external_v1` manifest (section 4.5). |
| 2026-11-27 | 2.3 (table row 1, the "Where the gap sits" paragraph, section 5 "Models fitted"; wording only) | Stated what variant 1's TCGA fit does with missing data, which the row had shown as a bare 594: 594 passed, 592 modelled; KPS set to the median, 80, for 150, with an indicator (HR 1.57); treatment arm missing for 22 (level, HR 1.20). Scoped "not imputed" to CGGA and to variants 2 and 3. Reframed the two internal contrasts as descriptive: neither carries an interval, and the first bundles three differences (the KPS value, the prognostic missing-KPS indicator, and fit sets of 592 against 570). No threshold, no decisive bound and no decision rule changed. | The independent review of commit `c8168a3` found the imputation disclosed in one place and contradicted in another, and an attribution drawn from two differences that were never bootstrapped. | **No.** No data was read. CGGA had been read twice for T2 before this entry. |
| 2026-11-27 | 2.3 (variant 1 internal C-index; note) | The code that computes variant 1's internal C-index imputed KPS with the median of all modelled patients before cross-validation, not the training-fold median its documentation stated. Corrected, and rerun on TCGA only: 0.673544 before and after (95% CI 0.6478 to 0.6994), because every training fold's median KPS is 80; the IDH sensitivity analysis is likewise unchanged (0.665650). The 0.6735 in the table is unchanged. No external number depends on the path: variant 1 is not evaluated on CGGA, and variants 2 and 3 impute nothing. | The review found that the code and its documentation disagreed. | **No.** TCGA only; no CGGA file was read. |
| 2026-11-27 | 3.3 (first sentence), 4.5 (`cgga_external_v1` role; a note below the table) | Section 4.5 named `cgga_external_v1` the primary external test for T1, T2, T4 and T5, and section 3.3 said final T1 numbers are reported on the frozen external test partitions. Both now say what section 2.2 items 9 and 10 already register: T1's subtype accuracy is scored on the labelled CGGA microarray patients (108 at ingest), fixed by rule and frozen as their own manifest, and T1's deconvolution on 100 mixtures simulated from the Neftel 2019 atlas. `cgga_external_v1` is the primary external test for T2, T4 and T5. No threshold, no decisive bound, no decision rule, no scoring set and no partition changed. | The two sections contradicted items 9 and 10: 80 of the 108 labelled patients have no RNA-seq profile and so fall outside the partition's eligibility (item 7), and the deconvolution row reads no cohort. Found in the stage 3b reconciliation (`docs/project_record/06_evidence_reports.md`, stage 3b, Part A item 2). | **No.** No T1 test data has been read: the CGGA microarray expression matrix has not been downloaded, the Neftel atlas is not on disk, and no T1 computation has been run. CGGA had been read twice for T2 before this entry; the change leaves T2's account as it was. |
| 2026-11-27 | 2.4 (added) | Recorded the outcome of T1's state-deconvolution row: PASS, mean per-state Pearson r 0.9521 (95% CI 0.9430 to 0.9589) against 0.75, on the 100 registered mixtures. The mixture-generation specification (`docs/t1_deconvolution_spec.md`) was committed before the run (`c68f556`) and fixes only what section 2.2 items 1 and 10 left open. No threshold, no decisive bound, no metric and no decision rule changed. | Records a result rather than a plan, in the form of section 2.3. | **Yes.** Written after the row was run, on 2026-11-27, before the deposit, and disclosed as such. The row reads the Neftel 2019 atlas, not a cohort; no CGGA or TCGA data was read for it. |
| 2026-11-27 | 5 (data in hand; data not yet obtained; analyses already run; the last sentence of "Models fitted") | Section 5 said no single-cell data had been ingested, that the single-cell reference atlases had not been obtained, and that the registration was written before any confirmatory evaluation other than T2. It now records that the Neftel 2019 atlas was obtained on 2026-11-27 for T1's deconvolution row, and that this row was evaluated before the deposit (section 2.4). | The section exists so that what is in hand at deposit is stated; it had become out of date. | **Yes.** Written after T1's deconvolution row was run (section 2.4). No CGGA data was read. |
| 2026-11-27 | 4.5 (`cgga_external_v1` row); 5 (data not yet obtained) | Sealed the external partition as section 2.2 item 7 registers: an all-test manifest over the 218 CGGA patients eligible by its rule (an RNA-seq profile, a survival outcome, a usable identifier), after the T2 evaluation and before any T4 or T5 computation. The row now gives n, split, seed and seal date, and the SHA-256 of the manifest and of the eligibility table. No threshold, decisive bound, decision rule or eligibility rule changed; the count equals the one registered in advance from identifiers (218). | Item 7 orders the seal after T2 and before T4 and T5; the registration must name the sealed partition. | **Yes.** CGGA had been read twice for T2 before this entry. The seal itself read identifiers, the batch label, the presence of a survival time and vital status, and the column names of the expression matrix: an availability access, recorded in the access ledger (entry 11). No value was inspected. |
| 2026-11-27 | 2.2 item 2 (completed) | Replaced "TO BE COMPLETED BY AUTHOR: the specific 2024 baseline model" with a rule for choosing the comparator: eligibility criteria E1 to E6, one PubMed search with a fixed query, the log to be committed and attached, our side of the comparison (the transferable model, variant 2), the partition (`cgga_external_v1`) and the paired bootstrap. If no model is eligible the threshold is NOT EVALUABLE and CGGA is not read for it; if one is, it is named by a further amendment before the deposit, with exactly one added read. No threshold, decisive bound or decision rule changed. | The protocol names no baseline. A comparator picked by hand cannot show that a stronger one was not passed over, and many 2024 signatures were fitted on the CGGA batches the partition is made of (`docs/t2_comparator_shortlist.md`). Fixed before the search is run. | **Yes.** CGGA had been read twice for T2 before this entry. No comparator has been run, and the registered search had not been run when this entry was committed. Four candidates had been read in a preliminary, unregistered search (`docs/t2_comparator_shortlist.md` section 4). |
| 2026-11-27 | 2.2 item 9(d) (classifier named) | Named the classifier that item 9(d) called "the RNA-seq-trained classifier": single-sample nearest centroid (Pearson) over the 831 resolved genes of the Verhaak 840-gene signature, against the frozen centroid set `verhaak840_tcga_all` (246 labelled TCGA patients, SHA-256 recorded), with a fixed tie order, no tuning and no other candidate; and the microarray input rules (exact symbol match, duplicate symbols averaged, log2 scale decided from the source's documentation before any value is read). The harmonisation rule itself is unchanged. | Decided in stage 3b (project record decision D6) on principle: a sample's call must not depend on the other samples scored. A classifier chosen after the deposit would be a post-hoc choice. | **No.** The CGGA microarray expression matrix has not been downloaded. The centroids were estimated on TCGA only; their agreement with TCGA labels was examined descriptively (section 5). CGGA had been read twice for T2 before this entry. |
| 2026-11-27 | 5 ("Models fitted"; "Analyses already run") | "Models fitted" still said `cgga_external_v1` did not exist; it now gives the seal date, and that the T1 centroid set was estimated on all labelled TCGA patients, the v2 test split included. "Analyses already run" listed four TCGA analyses; it now also lists the four of stage 3b: the two centroid sets, the single-sample subtype calls with their descriptive agreement check (56 of 76 v2-test patients), the scoring-set dependence measurement, and the feature-store build. | The section exists so that what was done before the deposit is stated; both passages had become out of date. Found while preparing `docs/deposit_checklist.md`. | **No.** No data was read. CGGA had been read twice for T2, and once for the seal's availability check, before this entry. |
| 2026-11-27 | 2.2 item 2 (comparator amended) | Kept the registered search and its candidates, and replaced "re-run as published" with a re-implementation under identical conditions. Each eligible 2024 gene set (E1, E2, E5; expression-only score with every gene and pair named; every gene present in TCGA and CGGA) is refitted with our Cox estimator on log2(CPM + 1). The comparator is the one with the best 5-fold cross-validated C-index on the v2 training split, refitted on all TCGA patients with RNA-seq and an outcome. Our side, the partition, the paired bootstrap and the threshold are unchanged. No threshold, decisive bound or decision rule changed. | The registered search, run once, found no model eligible to be re-run as published: almost no 2024 paper states the expression scale its coefficients expect (E3). The re-implementation removes that ambiguity, keeps CGGA out of gene selection and fitting, and compares both models under the same training data and estimator. Chosen by the project lead, knowing that a FAIL is the likely verdict (`docs/t2_comparator_shortlist.md`). | **Yes.** Written after the registered comparator search had been run and found no eligible model (13 records still unscreened). CGGA had been read twice for T2, and once for the seal's availability check, before this entry. No comparator has been fitted or scored; no CGGA value has been read for this row. |
| 2026-11-27 | 2.2 item 2 (g) (batch-stratified C-index) | For the improvement row only, each model's C-index counts a pair of patients as comparable only when both are in the same CGGA RNA-seq batch (mRNAseq_693 or mRNAseq_325), summed over both batches. The unstratified C-index is reported alongside, descriptively. The paired bootstrap, the partition and the threshold are unchanged, and the C-index row's verdict (section 2.3) is not recomputed. | The two batches were processed separately. A technical shift between them would lower the score of an expression-based comparator but not that of the clinical model, which would tilt the comparison for a technical reason. Fixed before any comparator is fitted. | **Yes.** Written after the comparator search and the option B amendment. CGGA had been read twice for T2, and once for the seal's availability check, before this entry. No CGGA value has been read for this row. |
| 2026-11-27 | 2.2 item 2 (j) added; 5 (data not yet obtained) | Recorded the completed comparator search and named the comparator. 1 of 549 records, PMID 38174685 (Sun and Liu, *Folia Neuropathol* 2024; 7 genes; trained on TCGA log2(FPKM + 1); CGGA used for validation only), meets the rule as first registered. It decides the verdict, applied as published, on CGGA's own RSEM FPKM release (fetched before the deposit, ledger 13; genes and patients checked by name only, ledger 16). Option B's refitted gene set becomes a secondary, descriptive comparison in the same read. One CGGA read is added for this row, after the deposit. No threshold, decisive bound or decision rule changed. | Option B was adopted while 13 records were still unread, on the premise that none was eligible. One of those 13 proved eligible under the rule fixed before the search. The lead chose to restore that rule for the verdict and keep B as secondary. | **Yes.** CGGA had been read twice for T2, once for the seal's availability check, and for the FPKM download (checksum) and its gene-name check (availability) before this entry. No CGGA expression value, clinical field or outcome has been read for this row. |
| 2026-11-27 | 2.2 item 2 (j) (option B's choice recorded) | Ran option B's registered choice before the deposit and recorded the result. 33 candidate gene sets were fixed and committed before any fit. R2 removed 10, read from CGGA gene names only (ledger entry 17). 23 were cross-validated on the v2 training split, and the best was the primary comparator's own 7 genes (C-index 0.615). The refitted model is frozen with its checksum. No threshold, bound or rule changed. | Item 2 (h) requires the chosen gene set to be named before the deposit. | **Yes.** As for the previous entry, plus the gene-name availability access of ledger entry 17. No CGGA expression value, clinical field or outcome has been read for this row. |
| 2026-11-27 | 4.5 (no text changed; storage of the three manifests) | The SHA-256 values in section 4.5 are of the manifests as written, with CRLF line endings. Git had stored them with LF line endings, so a checkout outside Windows produced different bytes and failed to reproduce the registered checksums, although the content was identical. The manifests are now marked byte-exact in `.gitattributes` and stored as written. No partition, id, seed or checksum changed. | A registered checksum must reproduce on any platform. Found while preparing the deposit (`docs/deposit_checklist.md` item 2.5); approved by the project lead. | **Yes**, as for the entries above. No data was read for this entry. |
| 2026-11-27 | 2.3 (decisions on T2's open items); 3.3 (GLASS removed); 3.4, 3.5 (as applied to T2); 4.5 (`glass_external_v1` row) | Wrote down five decisions that follow from registered text: the evaluation unit is the patient; T2 is FAIL on its C-index and Brier rows, and section 7's fallback applies; the two pre-seal reads are the reads of record, and T2 is not re-run; the KPS model is NOT EVALUABLE on CGGA; calibration and decision-curve analyses are not run for T2. GLASS is removed from section 3.3: no registered verdict needs it, because T2 is not re-run. No threshold, decisive bound, decision rule or partition changed. | These items were open in `docs/cgga_validation.md` sections 11 and 12.5 and `docs/deposit_checklist.md` items 1.3 to 1.5; a registration must not leave them implicit. | **Yes.** Written after CGGA had been read for T2. No data was read for this entry. |
| 2026-11-27 | Header; 5 (data in hand; analyses already run) | Completed the header: the author and affiliations, a committee not yet appointed, no funding and no conflicts, the start of the confirmatory phase immediately after the deposit, the repository URL, and the registered commit as the tag `prereg-deposit`. Only the DOI and the deposit date remain, and they exist only after the deposit. GLASS was removed from the IRB data list. Section 5 now lists the CGGA accesses that read no values (ledger entries 11, 13, 16, 17) and the TCGA cross-validation of option B. No threshold, bound or rule changed. | The deposit needs a complete header, and section 5 must be true at the moment of deposit. | **Yes**, as for the entries above. No data was read for this entry. |
| 2026-11-27 | 1.4 (H2 wording); 2.2 items 4, 6, 11 (added); 3.6 (scope); 3.8, 3.9, 3.10 (added); 5 | Answered the gaps found against the OSF study-design questions, with the project lead's choices. (1) T4 and T5 get a contingent model specification (item 11): a Stupp baseline; RMST to 24 months, success at a gain of 3 months or more and harm at a gain below 0; and the simulator and response surface frozen by amendment before any T4 or T5 read, or NOT EVALUABLE after 2027-06-30. (2) The configuration space is registered as coded, and the distinctness rule counts the number of injection sites, not a site set. (3) Sensitivity analyses run only inside registered reads. (4) Added a proportional-hazards test, a permutation negative control and a synthetic dry run; the variant 2 PH result on TCGA is disclosed. (5) No outlier removal. (6) Blinding is procedural only, with one analyst. No threshold, decisive bound or decision rule changed. | The OSF study-design form asks for each of these; the registration left them open or implicit. | **Yes**, as for the entries above. No CGGA value was read; the dry run used synthetic files, and its PH tests used TCGA. |
| 2026-11-27 | 2.2 item 5 (aggregate added) | Registered how T4's percentage regret reduction is aggregated: 100 x (1 - mean regret under GP-EI / mean regret under random search) over the 50 patients, with a patient-bootstrap CI. No threshold, bound or rule changed. | The registration defined regret per patient but not its aggregate; chosen by the project lead. | **Yes**, as for the entries above. T4 has not been run. |
| 2026-11-27 | 2.2 item 9 (scoring set frozen; input file named); 5 | **Made after the OSF submission of the registered commit `c7d7a95` (tag `prereg-deposit`), while approval was pending.** Froze T1's subtype scoring set by its registered rule. All 108 eligible CGGA mRNA_array_301 patients carry a reference label and have an expression column; they are sealed all-test as `data/partitions/cgga_t1_subtype_v1.json` (SHA-256 `7adfc09f068b793fa60fed0d191b34a39dc44c2f5b14cf4cc865d84c25d94049`), and `configs/scoring_rules.yaml` now names it. The classifier's input is CGGA's gene-level release (CGGA.mRNA_array_301_gene_level.20200506; SHA-256 in `data/provenance/cgga_microarray_sources.json`); the registered rule matches genes by symbol. No threshold, bound or rule changed. | Section 3.3 and item 9 require the scoring set to be frozen as a manifest before the read. | **No** T1 data read. The microarray file was fetched and hashed (ledger entries 19 to 21; the first attempt failed with HTTP 500). The freeze read only sample ids and whether a label is present, never a label or an expression value (entry 22). |
| 2026-11-27 | 2.2 item 9(d) (input scale decided; two format contingencies) | **Made after the OSF submission of `c7d7a95`, before any microarray value was read.** Decided the input scale from documentation, as registered: log2(x + 1), because CGGA's sources say only that intensities were normalised with GeneSpring GX 11.0. Added two format contingencies so that the one read cannot fail on format: a negative value means the data are already log-scale and are used as distributed; a missing value drops that gene for that sample. No threshold, bound, scoring set or classifier changed. | Item 9(d) requires the scale to be decided from documentation before the read. A read that crashes on an unexpected format would waste the one registered read. | **No.** No T1 value was read: the matrix has been opened only for its header (sample ids). |
| 2026-11-27 | 2.2 item 3(b), 3(e); 2.2 item 8 (decision date added) | **Made after the OSF submission of `c7d7a95`, before any T3 data was obtained or any T3 source was contacted.** Two changes to T3's data rule. (1) The qualifying rule required a source spanning all four Verhaak subtypes; it now requires at least three of the four, with at least two curves per represented subtype and the ten-curve minimum unchanged. Any absent subtype is named in this section when the data arrive, and the coverage limit is stated wherever T3's verdict is reported. (2) Both T3 rows are decided on 2027-01-08: by that date each row must have a qualifying source secured (a verified open release, or a custodian's written agreement), or that row is declared NOT EVALUABLE on that date. The rows are decided independently. Neither T3 threshold (RMSE <= 15%, migration r >= 0.65), nor the RMSE normalisation, nor the bootstrap unit, nor any decisive bound or decision rule changed. | The four-subtype requirement was written without a verified source in hand, and the registered search (`docs/t3_data_search.md`) found no public panel that pairs longitudinal control-arm volumes with subtype calls across all four. Both candidates classify three classes (classical, mesenchymal, proneural) under the Wang 2017 scheme, which has no Neural class at all, so the missing subtype is a property of how modern panels are labelled rather than of any one panel; Wang 2017 dropped it because its signal tracked normal-brain contamination. The ten-curve minimum is unaffected: the Mayo panel's published counts are 31, 21 and 16. A rule that no obtainable dataset can satisfy makes NOT EVALUABLE the only reachable verdict, which is not a finding about the model. The per-subtype floor of two is added at the same time so that the loosened rule cannot be met by a panel that is nominally three-subtype and effectively one. The decision date exists because item 11 freezes the T4 and T5 model before 2027-06-30, and that freeze must be made knowing whether the simulation layer is calibrated; an open-ended search would let T3's status be settled after its convenience was known. Chosen by the project lead. | **No** T3 data read, and none obtained: no PDX growth curve and no NSC-tracking series is on disk, no custodian has been contacted, and the loosened rule was fixed before any request was sent. CGGA had been read twice for T2 and for the availability and checksum accesses of ledger entries 11, 13, 16, 17 and 19 to 22 before this entry; none of them touches T3. |
| 2026-11-27 | 2.2 item 11 (g) to (o) added | **Made after the OSF submission of `c7d7a95`, before any T4 or T5 computation and before any CGGA value was read for either.** Adopted Design A of `docs/t4t5_model_design.md` and fixed the nine choices that document left open, as decided by the project lead. (g) The frozen survival model is variant 2 (age, sex, arm), not variant 3: variant 3 needs an IDH and an MGMT call for every patient, which item 7's eligibility rule never required, so T5 would score fewer than the 218 its scoring set registers. (h) The simulator is a closed-form delivery module, declared uncalibrated whatever T3 returns. (i) The efficacy prior is a rule — Emax_p = 0.5 x ln(treated-to-control median-survival ratio of the nearest preclinical study) — which places T5 below its bar by construction, stated in advance. (j) "combo" is defined per carrier, as the two payloads that carrier may carry, with Bliss independence and synergy 1.0. (k) The toxicity and delay term is in the primary surface, because it is the only source of harm and T5's harm row could not otherwise fail. (l) All 218 patients, no exclusion by IDH status, and neither IDH nor MGMT anywhere in the frozen model; the resulting mis-modelling of IDH-mutant baselines is disclosed and reported in section 3.6. (m) The freeze amendment is entered by 2027-02-28. (n) Parameters live in `configs/response_surface_v1.yaml`, checksummed here at the freeze. (o) Four pre-read controls are added to section 3.8, any failure of which stops the read. No threshold, decisive bound, decision rule, partition or scoring set changed; item 4's configuration space, its two exclusions, its 10 pairs and its 240 combinations are untouched, as are (a) to (f) of item 11. | Item 11(e) requires the model to be frozen by amendment before any T4 or T5 read, and a model cannot be frozen while nine of its choices are open. Each choice is one that could otherwise be made after seeing the data. Three are recorded against the design document's own recommendation, with the reason given in the text: variant 2 over variant 3 (the 218-patient scoring set), the per-carrier "combo" over a uniform one (all three options in the document contradict item 4's exclusion list or its distinctness rule), and IDH and MGMT in no part of the model rather than as invented modifiers. This entry fixes the design; the numeric values and their checksum follow in the freeze amendment under (e) and (m). | **No.** No T4 or T5 computation has been run, and no CGGA value has been read for either: the response surface and the delivery module are stubs, and `configs/response_surface_v1.yaml` does not yet exist. CGGA had been read twice for T2 and for the availability and checksum accesses of ledger entries 11, 13, 16, 17 and 19 to 22 before this entry. The external C-index and Brier values quoted in (g) are T2's, already recorded in section 2.3. |

---

## 9. Correspondence between this document and the code

| This document | Authoritative implementation |
| --- | --- |
| §2 threshold table | `src/gliostrata/evaluation/thresholds.py` (`THRESHOLDS`, checksum `THRESHOLD_TABLE_SHA256`) |
| §6 decision rules | `gliostrata.evaluation.thresholds.evaluate` → `Verdict.status`, `Verdict.reason` |
| §4 partition procedure | `src/gliostrata/evaluation/partitions.py` (`partition_cohort`, `freeze_partition`) |
| §4.2 manifest | `PartitionManifest`; instances in `data/partitions/*.json` |
| §4.3 verification | `verify_partition` → `DriftReport`; CLI `gliostrata freeze-partition --verify` |
| §4.4 sealing | `write_manifest(..., force=...)`, `PartitionSealedError` |
| §4.4 sealing, all-test | `freeze_partition(..., test_fraction=1.0)`; CLI `gliostrata freeze-partition --all-test` |
| All of the above | `tests/test_evaluation.py` |

Where this document and the code disagree, the code's checksum and manifests are the
evidence of what was registered, and the discrepancy is itself a finding to be disclosed.
