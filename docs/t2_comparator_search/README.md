# T2 comparator search: log

Status: **closed 2026-11-27. One model is eligible under the rule as first registered:
PMID 38174685** (section 3). Of the 13 papers that had no open full text, 10 were read from
PDFs the author obtained through a university library and 3 could not be obtained. The
eligible paper was one of the 10. Earlier the same day, while those 13 were still unread,
the rule had been amended to a refitted gene set (option B). After the eligible model was
found, the lead chose to restore the rule as first registered for the verdict. Option B is
kept as a secondary, descriptive comparison (pre-registration section 2.2 item 2).

This is the single search registered in pre-registration section 2.2 item 2 (b),
committed in `aacedb1` before the search was run. The registered text is authoritative.
The reasoning behind it is in `docs/t2_comparator_shortlist.md`.

## 1. The search

`search.json` holds the query, PubMed's translation of it, the time it was run
(2026-11-27T17:48:54Z, once) and all 549 PMIDs. PubMed applied the date limit as
2024/01/01 to 2024/12/31 [Date - Publication].

## 2. How each record was judged

| Stage | Criteria | Who | Records |
| --- | --- | --- | --- |
| Date and type | E1: earliest of the electronic and print dates must fall in 2024; preprints and retracted publications fail. Read from PubMed's metadata by script | script | 549 |
| Title and abstract | E1 (non-articles), E2, E4. Protocol: `protocols/1_title_abstract.md` | two independent screeners per record; 4 disagreements adjudicated from the abstract | 480 |
| Full text | E3, E5, plus any E2/E4 failure that only the full text shows. Protocol: `protocols/2_full_text.md` | one reviewer per paper; every quote checked against the text by script | 71 open-access texts (Europe PMC) and 2 publisher pages |
| Second review | All of E2 to E5 again, on the 25 closest calls: the two that passed, the 6 with coefficients in a supplement (supplements downloaded and read), 16 with complete coefficients but a doubtful scale, and 1 doubtful E5 call. Protocol: `protocols/3_second_review.md` | a second, independent reviewer who had not seen the first verdicts; 3 calls adjudicated | 25 |

`records.jsonl` has one line per PMID. Each line gives the title, journal, DOI and
first-publication date, the stage the record stopped at, its first failing criterion,
the reason, and how it was screened.

## 3. Results

| Stopped at | First failing criterion | Records |
| --- | --- | --- |
| Date and type | E1 | 69 |
| Title and abstract | E1 | 1 |
| Title and abstract | E2 | 330 |
| Title and abstract | E4 | 61 |
| Full text | E3 | 80 |
| Full text | E5 | 2 |
| Full text not obtained | E3 (not established) | 5 |
| **Eligible** | | **1** (PMID 38174685) |

At the full-text stage, reviewers also found that 17 papers had used CGGA mRNAseq_325 or
_693 for training or gene selection (E5). Those papers are counted above under E3, their
first failing criterion.

**The eligible model.** Sun G, Liu W. "The neutrophil extracellular traps-related gene
signature predicts the prognosis of glioblastoma multiforme". *Folia Neuropathologica*
2024;62(1):59-75, doi:10.5114/fn.2023.132980. It is a 7-gene LASSO-Cox risk score, trained
on 167 TCGA-GBM patients and validated on 237 CGGA mRNAseq_693 GBM patients:

RS = 0.223340457·CPPED1 + 0.1517013·F3 + 0.046311722·G0S2 - 0.219913858·MAPK1 +
0.192197058·MME + 0.016505054·MMP9 + 0.414575274·MPO

The training expression is TCGA-GBM "RNA-seq data (log2(fpkm+1))" from UCSC Xena.

- **E1 (2024):** first published in print in 2024. Crossref records no earlier online
  date, and created the record on 2024-01-04. The "2023" in the DOI is the year the DOI
  was assigned.
- **E4 (inputs available):** confirmed by an availability check of CGGA's FPKM release
  (access ledger entry 16). All 7 genes are present in both batches, and all 218
  partition patients have a column.
- **E6 (the best of the eligible):** does not apply; it is the only eligible model.

## 4. Records without a full text

- **3 could not be obtained**, even through the author's library: 38959943, 38318849 and
  37539935. They are recorded as "full text not obtained; E3 not established".
- **2 fail E5 on their abstract**, as trained on CGGA: 39476302 and 38662286.

The 10 PDFs that were obtained are kept, unredistributed, in
`data/external/t2_comparator_fulltext/` (git-ignored).

## 5. Decisions made during the screen, disclosed

1. **The scale rule was made precise after the first full-text pass.** It now requires,
   for the data the coefficients were fitted on, the unit and every transformation
   applied before fitting. Two cases count as no stated scale: values batch-corrected
   jointly with other cohorts, and coefficients readable only from a figure. The rule was
   then applied the same way to every remaining record. It interprets registered E3
   ("the expression scale it expects") and adds no criterion.
2. **39435104 is ruled ineligible under E3, against both reviewers.** It is a 22
   gene-pair model (Wang et al., *Heliyon* 2024). Both reviewers passed it, because a
   protocol sentence said gene-pair models need no scale. That holds only for transforms
   applied to every gene alike. Raw counts and length-normalised values (FPKM or TPM)
   reverse a pair comparison whenever the two genes differ in length, and the paper does
   not say which it used. The author confirmed the ruling on 2026-11-27.
3. **39630160:** the paper describes its scale in two contradictory ways, so E3 is not
   met. The two reviewers had split on this record.
4. **38378721:** its cut-off was chosen within the CGGA mRNAseq_325 and _693 cohorts,
   which fails E5.
5. **Doubt is never resolved in favour of eligibility.** This is the registered rule.

## 6. What was not read

- **Unopened supplementary files:** the 412 MB supplementary spreadsheet of 38320988
  (that paper fails E3 for other reasons), and TIFF images.
- **Unread publisher pages:** those whose sites refused automated access.
- **CGGA:** no CGGA file was opened for this search.
