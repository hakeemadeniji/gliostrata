# T2 comparator search: log

Status: **run 2026-11-27, not yet closed.** 13 of the 549 records still await their full
text (section 4). Under the rule as first registered (re-run as published), no model is
eligible so far. Later the same day, before the deposit, the rule was amended (option B):
the same candidates are judged as *gene sets*, to be refitted under the project's own
conditions. Pre-registration section 2.2 item 2 holds the amended rule.

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

## 3. Results so far

| Stopped at | First failing criterion | Records |
| --- | --- | --- |
| Date and type | E1 | 69 |
| Title and abstract | E1 | 1 |
| Title and abstract | E2 | 330 |
| Title and abstract | E4 | 61 |
| Full text | E3 | 71 |
| Full text | E5 | 2 |
| Full text not obtained | E3 (not established) | 15 |
| **Eligible** | | **0** |

At the full-text stage, the first reviewer also found that 16 papers had used CGGA
mRNAseq_325 or _693 for training or gene selection (E5). Those papers are counted above
under E3, their first failing criterion.

## 4. Records still open

Of the 15 papers with no open full text:
- 2 fail E5 on their abstract: 39476302 and 38662286, both trained on CGGA.
- 13 cannot be judged without the full text: 39736337, 39322177, 39216044, 39216042,
  39162363, 38959943, 38719174, 38615384, 38462406, 38369937, 38318849, 38174685 and
  37539935.

The author is obtaining those 13 through a university library. They will be screened by
the same protocols before the outcome is recorded.

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
