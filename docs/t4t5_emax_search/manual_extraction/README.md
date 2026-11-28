# Filling in `manual_extraction.csv`

19 PDFs, 20 rows. There are 20 because PMID 37420122 is screened under two payloads
(`prodrug_enzyme` and `CAR`) and each gets its own row — it is one paper read twice, once per
payload.

Open `manual_extraction.csv` in Excel or Google Sheets. `pmid`, `payload`, `title` and `doi`
are filled in; everything else is yours. **Fill the rows in the order they appear**: CRAd
first (3), then CAR (4), then prodrug_enzyme (13). That order is deliberate — CRAd currently
rests on a single study and CAR's primary set is empty, so those seven rows are where the
answer actually moves. If you run out of time, stopping after the first seven still helps.

## Matching a PDF to a row

The filenames are the publishers', not PMIDs. Open a PDF, find its **DOI** on the first page,
and match it to the `doi` column. Put the filename in `pdf_filename` so the file and the row
are linked in the record.

If a PDF turns out not to match any row, or a row has no PDF, leave the row blank and say so
in `note`. Item 11(i.8) covers that case: it is recorded, not patched.

## What goes in each column

| Column | What to write |
|---|---|
| `pdf_filename` | the PDF's filename, e.g. `djv375.pdf` |
| `page` | the page number the survival numbers are on (helps me spot-check) |
| `treated_median` | median survival of the arm that got **this row's payload**, as a number only (`39.5`, not `39.5 days`) |
| `control_median` | median survival of the **control arm in the same experiment**, number only |
| `unit` | `days`, `weeks` or `months` — whatever the paper uses. Both numbers must be in the same unit. |
| `arm_is_combination` | `TRUE` if the treated arm also got something else (radiation, temozolomide, a second agent); otherwise `FALSE` |
| `combined_with` | if `TRUE` above, what it was combined with |
| `control_arm_kind` | one of five values — see below |
| `p2_intracranial` | `yes` if the tumour was put **in the brain** (orthotopic, intracranial, intracerebral, pons, brainstem) of a mouse or rat; `no` if under the skin (subcutaneous, flank); `uncertain` if the paper doesn't say |
| `p4_payload_match` | see below |
| `is_primary_research` | `FALSE` if it's a review, commentary or methods chapter; otherwise `TRUE` |
| `survival_quote` | copy-paste the sentence the two numbers came from |
| `control_arm_quote` | copy-paste the sentence that says what the control arm received |
| `note` | anything odd: several experiments, a number you weren't sure about, a median "not reached" |

### `control_arm_kind` — the decisive one

This is what your deviation (item 11(i.7)) turns on. Judge it by **what the control animals
received**, never by how long they lived:

- `untreated` — nothing at all
- `vehicle` — saline, PBS, buffer, or an empty carrier only
- `mock_or_reporter_vector` — mock transduction, an empty vector, or a vector carrying only a
  marker like GFP, DsRed or luciferase and no therapeutic gene
- `active_therapeutic` — the control got something expected to fight the tumour: plain
  (untransduced) immune cells, an oncolytic virus on its own, radiation, chemotherapy
- `uncertain` — the paper doesn't say, or you can't tell whether what the control got has
  anti-tumour activity of its own

**`uncertain` is a perfectly good answer and is often the right one.** Only the first three
enter the primary analysis, so an honest `uncertain` keeps a study out — which is what the
rule is for. Please don't stretch to one of the first three.

Two cases that came up repeatedly in the automated pass, for consistency:
- An **irrelevant-target CAR** control (e.g. a CD19 CAR-T used against a brain tumour) →
  `uncertain`, unless the paper states it has no anti-tumour effect in that model.
- A vector that carries the **enzyme but no prodrug was given** → `active_therapeutic` if the
  vector itself does something to the tumour, otherwise `uncertain`.

### `p4_payload_match`

Does the treated arm actually use this row's payload type?

- `prodrug_enzyme` → needs the **activating drug named**: ganciclovir, 5-FC, CPT-11, CB1954.
  Enzyme alone, with no prodrug named → `uncertain`.
- `CRAd` → needs a **replicating** oncolytic adenovirus. A replication-defective vector used
  only to deliver a gene → `no`.
- `CAR` → needs a chimeric antigen receptor T-cell or NK-cell product.

## Rules that matter more than completeness

1. **A blank is always acceptable. A guess is not.** Every blank is handled by a registered
   rule; a wrong number would go straight into the registered prior.
2. **Both numbers from the same experiment.** If a paper has several survival experiments,
   take the one it presents as the main one and say in `note` what the others were.
3. **Prefer the payload on its own** over a combination arm.
4. **"Not reached" or "undefined" is a blank**, with `note` saying so.
5. **Don't look at the ratio.** If you find yourself checking whether a number helps or hurts
   T5, that is the moment the registration exists to protect against. Record what the paper
   says.

## When you're done

Save the CSV in place and tell me. I'll validate it (units, ranges, quotes present where
numbers are), fold it into `extraction_final.jsonl`, and recompute. Any row I can't validate
comes back to you as a question rather than being quietly fixed.

## Correction, 2026-11-27: two wrong DOIs, two wrong PDFs

`pdf_filename` is now pre-filled for **17 of the 19** papers, each confirmed by matching the
paper's recorded title against the text of the PDF — not by its DOI, for the reason below.
Those 17 rows are ready to fill in.

**Two rows have no usable PDF, and that is my error.** The `doi` field in `records.jsonl` was
read from the last `<ArticleId IdType="doi">` in each PubMed record, which for any article
carrying a reference list is a *cited* paper's DOI rather than the article's own. PMID
37420122 contains 40 such identifiers, 39 of them references. 37 of the 108 search records
were affected; two of them were in the download list, so two of the PDFs retrieved are the
wrong papers:

| Row | DOI I gave | Correct DOI | What was downloaded instead |
|---|---|---|---|
| 37420122 | `10.1158/1078-0432.CCR-09-0155` | `10.1007/s13402-023-00842-5` | Candolfi 2009, HMGB1 release in glioma killing (`4401.pdf`) |
| 22996231 | `10.2174/156652309788488569` | `10.1007/s13311-012-0144-7` | Tas 2009, NF-κB gene therapy in inflammatory disease (`0002Q.pdf`) |

Neither downloaded paper is among the 108 search records, so **neither may be used** —
they are outside the evidence base item 11(i.8) closed over.

The DOI field was never used by any screening or extraction decision: screening ran on titles
and abstracts, and matching ran on titles. So the only consequence is the two wasted
retrievals. All DOIs in `records.jsonl` and in this sheet are now corrected, and
`scripts/t4t5_emax_search.py fix-dois` plus `doi_correction.json` record every change.
