# T2 comparator: shortlist and the decision it needs

Status: **prepared 2026-11-27 for the project lead. Nothing here is registered.** It fills
no "TO BE COMPLETED BY AUTHOR" field. Pre-registration section 2.2 item 2 stays open until
the lead decides, and the decision is then entered there as a dated amendment (section 8).

No CGGA file was opened for this document, and no model was run. The candidates were
judged on what their papers say.

## 1. What the registration already fixes

- **The row.** T2 must beat the "best-published 2024 baseline" by at least 0.03 C-index, and
  the lower bound of the 95% CI is the decisive bound (section 2, the threshold table;
  checksum `96bbcd49…74ab`). The protocol names no baseline (`docs/protocol.md`, Table 8-2).
- **The metric** (section 2.2 item 2). The comparator is re-run on the identical frozen
  partition, and the metric is the paired bootstrap distribution of the difference. A C-index
  quoted from a paper is never used.
- **The fallback** (item 2; section 6; OSF stopping rule 8). If the comparator cannot be
  re-run, the row is NOT EVALUABLE, never a pass.

## 2. Five facts that constrain any choice

1. **The registered model cannot be scored on CGGA.** Variant 1 (age, sex, KPS, treatment
   arm) needs KPS, and CGGA records none. The C-index and Brier rows were therefore decided
   by variant 2 (age, sex, treatment arm), the transferable model, which was pre-specified as
   the headline in the T2 script before the first read. Nothing yet says which of our models
   enters the improvement row. That has to be fixed before any read.
2. **The only frozen external partition is `cgga_external_v1`** (218 patients, all with
   RNA-seq; sealed 2026-11-27). T2's C-index row used 283 patients of the ingested cohort,
   with no sealed partition (a disclosed deviation, section 2.3). "Identical frozen
   partition" can now only mean `cgga_external_v1`. So both models would be scored there,
   variant 2 included, on complete cases for the union of their inputs.
3. **A run would be a fourth value read of CGGA, and the registered text does not allow
   it.** OSF stopping rule 9 and section 5 list T1, T4 and T5 as the only further reads,
   and the OSF text says "T2 is reported as already conducted and is not re-run". Running
   the comparator needs an amendment before the deposit that adds exactly one read. After
   the deposit, the read guard admits a value read only under a rule id in
   `configs/scoring_rules.yaml`, so the comparator script would also need its own rule.
4. **The 2024 literature overlaps heavily with the partition.** Many 2024 GBM expression
   signatures were fitted, or had their genes selected, on CGGA mRNAseq_325 or mRNAseq_693.
   Those are the two batches the partition is made of (85 and 133 patients). A comparator
   fitted on them would be scored on its own training patients.
5. **The likely outcome is known in advance. It must not decide the path.** Variant 2 scored
   0.5385 (95% CI 0.4994 to 0.5807) on 283 CGGA patients. Clearing a paired difference of
   +0.03 at the lower bound would, on those numbers, need a comparator at or below chance
   (variant 2's C-index on the 218 is not known and may differ). If a comparator is
   run, FAIL is the expected verdict. That is a reportable result, and it is not a reason to
   prefer NOT EVALUABLE. The path should be chosen on the eligibility rule below, as if the
   outcome were unknown.

## 3. Proposed eligibility rule (to be registered before anything is run)

A model is eligible only if it meets all six criteria.

| # | Criterion | Why |
| --- | --- | --- |
| E1 | Peer-reviewed, first published (online or in print) between 2024-01-01 and 2024-12-31 | The protocol's "2024 baseline" |
| E2 | Predicts overall survival, and applies to adult WHO grade 4 glioblastoma, IDH-wild-type included | The outcome and population of T2 |
| E3 | Can be re-run without refitting: either a complete published formula (every term, its coefficient and the expression scale it expects) or released code with fitted weights | Item 2: "re-run", never refit on the test partition |
| E4 | Every input is available for `cgga_external_v1` from files already on disk: CGGA RNA-seq and CGGA clinical fields (age, sex, IDH, MGMT, 1p/19q, radiotherapy, chemotherapy). Not KPS, extent of resection, methylation arrays, imaging or histology | Otherwise the comparator is not scorable on the identical partition |
| E5 | No CGGA mRNAseq_325 or mRNAseq_693 patient used for training, gene or feature selection, coefficient fitting or cut-off choice. Use for validation only is allowed, and disclosed | A comparator scored on its own training patients is not a re-run on a held-out partition |
| E6 | Among eligible models: the highest C-index the paper reports on an external cohort other than CGGA. Ties go to the earlier publication. If no eligible model reports one, the largest external validation n | Makes "best-published" a rule, not a pick, and keeps CGGA results out of the choice |

**The search, fixed with the rule.** Run it once, on PubMed, over publication date 2024,
with a query such as `(glioblastoma[tiab]) AND (prognos*[tiab] OR survival[tiab]) AND
(signature[tiab] OR model[tiab] OR nomogram[tiab] OR "risk score"[tiab]) AND
2024[dp]`. Screen by title and abstract for E2 and E4, then by full text for E3 and E5.
Commit the full log: query, date run, hit count, and each candidate with its first failing
criterion. Expect hundreds of hits, and screen them all. The point of the log is that a
reader can check that nothing was skipped.

## 4. Candidates checked so far (a targeted search, not the registered one)

These come from a targeted web search on 2026-11-27, with three papers read in more depth.
This is a preview of what the registered search will find. It is not that search.

| Model | Inputs | E1 | E3 | E4 | E5 | Outcome |
| --- | --- | --- | --- | --- | --- | --- |
| Han X et al., "Identification of an immune signature…for IDH wildtype glioblastoma", *J Cancer* 15(19):6452–6467 (2024), doi:10.7150/jca.100144 | 3 genes (FMOD, MXRA5, RAB36). Published formula: 0.159·FMOD + 0.228·MXRA5 + 0.684·RAB36 | yes | yes; the expression scale must be confirmed | yes | **no**: fitted on CGGA_325, with genes selected by intersecting TCGA and CGGA | **Ineligible (E5)** |
| Zhao Z et al., "Development of a web-based tool for estimating individualized survival curves in glioblastoma…", *Clin Transl Oncol* 27(5):2113–2126, doi:10.1007/s12094-024-03739-3 | SEER clinical variables; an 11-gene mRNA model; TME features | online 2024-10-30, issue 2025 | not verified: abstract only; a Shiny app is linked | no: needs marital status, race, laterality, site and surgery | **no**: the mRNA model was built on merged TCGA, CGGA and GEO data | **Ineligible (E4, E5)** |
| Drexler R et al., "A prognostic neural epigenetic signature in high-grade glioma", *Nat Med* 30(6):1622–1635 (2024), doi:10.1038/s41591-024-02969-w (PMC11186787) | DNA methylation arrays (450k/EPIC) | yes | code released; no coefficients in the text read | **no**: the partition has no methylation arrays | yes: no CGGA in its data statement | **Ineligible (E4)** |
| MRI-based 2024 survival models (e.g. the ReSPOND model, code at `CBICA/GBM_Survival_ReSPOND_2024`; Zhang et al., *npj Precis Oncol* 2024, doi:10.1038/s41698-024-00670-2) | Imaging | yes | code released | **no**: no imaging for CGGA | not checked | **Ineligible (E4)** |

How these were read:
- **Han et al.:** the full text on PMC (PMC11540507). It reports no C-index, only
  time-dependent AUC: 0.63 at 1 year and 0.61 at 3 years on TCGA.
- **Zhao et al.:** the Europe PMC abstract only. There is no open full text.
- **Drexler et al.:** the open full text on PMC, partly through a summary. Whether the
  classifier's weights are released was not established.

None of the four candidates is eligible. A wider search may find one. A model trained on
TCGA alone, with a published coefficient table, validated outside CGGA, and needing only
RNA-seq plus age, sex, IDH or MGMT would pass E3 to E5. That kind of model exists in
earlier years. Whether a 2024 one exists is what the registered search answers.

## 5. The decision for the lead

**Path C (recommended): register the rule and the search now, then let the search decide.**
- Before the deposit, amend section 2.2 item 2 to give:
  - the eligibility rule (E1 to E6) and the search;
  - our side of the comparison: variant 2, the transferable model. It is the same stand-in
    that decided the C-index and Brier rows, and it is fitted on TCGA exactly as committed;
  - the partition (`cgga_external_v1`), complete cases on the union of both models' inputs;
  - the procedure: a paired percentile bootstrap of the difference in Harrell's C-index,
    1000 resamples, seed 42, as for the C-index row.
- Then run the search once and commit its log.
- **If no model is eligible**, the row is NOT EVALUABLE under the existing rule, with the
  log as the evidence. No further CGGA read is needed.
- **If one is eligible**, name it in the same amendment, and add exactly one read to OSF
  stopping rule 9 and section 5, run after the deposit under a new rule id in
  `configs/scoring_rules.yaml`.

**Path A: name a comparator without a registered search.** This is faster, but the choice
is a pick, not a rule. A reviewer cannot tell whether a stronger 2024 model was passed
over.

**Path B: declare the row NOT EVALUABLE now, on the ground that the registered model
(variant 1) cannot be scored on CGGA.** This is defensible under section 6, since the KPS
data cannot be obtained. But it treats the improvement row differently from the C-index and
Brier rows, which were decided by variant 2. A reader can fairly ask why the stand-in counts
for two rows and not the third, so the text would have to answer that.

Paths A and C both treat the improvement row like the other two T2 rows, which is why C is
recommended. In every path, the amendment is entered before the deposit and marked
"Test data already read? **Yes**", because CGGA was read twice for T2.
