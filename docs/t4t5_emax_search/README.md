# The registered Emax search (pre-registration 2.2 item 11 (i.1) to (i.6))

This folder is the audit trail for the one search that sets T4 and T5's efficacy prior.
Nothing here may be revised: item 11(i.6) states that (i.1) to (i.5) are closed once the
search has been run, and it has been.

## What the rule is

    Emax_p  = 0.5 * ln(ratio_p)                       (i.1); 0.5 is assumed
    ratio_p = median over every qualifying study of   (i.2)
              (treated median survival / control median survival)

`ratio_p` is a **median across all qualifying studies**, not one study chosen by hand. That
is the point of the rule: it leaves no discretion over which paper sets the prior.

## Qualifying criteria, fixed before the search (i.3)

| | Criterion |
|---|---|
| P1 | Peer-reviewed and in vivo. |
| P2 | Orthotopic or intracranial rodent glioma or glioblastoma model. |
| P3 | Median survival reported for both a control arm and an arm receiving the payload, in the text or readable from a survival curve. |
| P4 | Payload class matches: CRAd a conditionally replicative oncolytic adenovirus; prodrug_enzyme an enzyme-prodrug pair (cytosine deaminase with 5-FC, HSV thymidine kinase with ganciclovir, carboxylesterase with CPT-11, or equivalent); CAR a chimeric antigen receptor T-cell or NK-cell product. |
| P5 | Carrier delivery is **neither required nor disqualifying**: the delivery module models carriage separately, so a carrier-delivered study qualifies on the same terms as a direct one. |

**Fallback (i.5).** If no study qualifies for a payload, P2 is dropped for that payload
alone and the same rule is applied to the nearest tumour type available; the disease
mismatch is then reported with every number for that payload. If still nothing qualifies,
`Emax_p = 0` and the payload is reported as unevidenced.

## The run

One PubMed query per payload, verbatim from item 11 (i.4), no year limit, run once by
`scripts/t4t5_emax_search.py search` on **2026-11-28T01:19:52Z** (2026-11-27 in the project
lead's local time; the log records UTC). The script refuses to run again.

| Payload | Records |
|---|---|
| CRAd | 10 |
| prodrug_enzyme | 55 |
| CAR | 43 |
| **Total** | **108 records, 107 unique PMIDs** |

One record (PMID 37420122) is returned by both the CAR and the prodrug query. It is screened
once per payload, because P4 is payload-specific and it may qualify for one and not the other.

PubMed's own translation of each query is recorded in `search.json` and was checked against
the registered text: every term stayed `[Title/Abstract]`, with no MeSH expansion added.

## Files

| File | What it holds |
|---|---|
| `search.json` | The three queries, PubMed's translation of each, the counts, the PMID lists, the cross-payload overlap, the run timestamp. |
| `records.jsonl` | One record per line: PMID, title, abstract, journal, year, DOI, publication types, payload. |
| `screening/` | The two independent screens, the adjudication, and the extracted ratios. |

## Screening protocol

Stage 1, titles and abstracts: two independent screens, each deciding INCLUDE / EXCLUDE /
UNCLEAR against P1, P2, P4 and P5, with the failing criterion named and a quoted phrase for
every exclusion. P3 is not decided at this stage, because an abstract that omits median
survival does not prove the paper does.

Stage 2, full text or abstract: for everything not excluded at stage 1, the treated and
control median survivals are extracted with the sentence they came from. A record whose
median survival cannot be found in what is accessible fails P3 and is excluded with that
reason recorded — not silently dropped.

Disagreements between the two screens are adjudicated, and the adjudication is recorded with
its reason. Counts at every stage are reported, so that the path from 108 records to the
final ratio set is reconstructable.
