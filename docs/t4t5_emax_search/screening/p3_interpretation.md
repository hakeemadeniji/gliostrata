# P3's "control arm": an ambiguity found after the search was run

> ## SUPERSEDED AS PRIMARY, the same day (2026-11-27)
>
> This file first recorded the **literal** reading of P3 as primary, with the restricted
> reading as an exploratory check. Later the same day the project lead reversed the two: the
> **restricted** reading is now primary, registered as pre-registration item 11(i.7) and
> logged in section 8 as a **disclosed deviation** from (i.3) — not an amendment to it, since
> (i.6) closed (i.1) to (i.6) when the search ran.
>
> | | Primary | Exploratory |
> |---|---|---|
> | First decided | literal P3 (any designated control) | restricted (untreated only) |
> | **Now** | **restricted: control arm receives no agent intended to have anti-tumour activity** | literal P3, with the count dropped |
>
> If a payload's restricted set is empty, it falls back to the literal set for that payload,
> so the deviation can never be worse-evidenced than the registration itself.
>
> The reasoning below is kept **unedited**, because it is the record of what was decided
> first and why. Where it says the literal reading is primary, read item 11(i.7) instead.
> The reversal was taken openly, with its effect on the verdict unknown but its direction
> known and against the project's own hypothesis. See section 8.

**Status: an interpretation of a closed criterion, not a revision of it.** Pre-registration
item 11(i.6) says that nothing in (i.1) to (i.5) may be revised once the search has been run.
The search ran on 2026-11-28T01:19:52Z. P3 therefore stands exactly as written, and what
follows is a record of how it is applied, not a change to it.

## The ambiguity

P3 requires "median survival reported for both **a control arm** and an arm receiving the
payload". It does not say what a control arm is, and the extracted studies use at least three
different kinds:

| Kind | Example | Effect on the ratio |
|---|---|---|
| Untreated, vehicle or mock | most studies | largest ratio |
| Vector without the activating prodrug | PMID 32093290 (GALV-NTR with and without CB1954) | smaller: isolates the prodrug |
| Parental, untargeted effector cells | PMID 26640245 (parental NK-92 against CAR NK-92) | smaller: isolates the receptor |

The three give systematically different ratios, so the choice is not cosmetic: it moves
`Emax_p`, and with it T5's verdict.

## How it is applied

**Literally.** P3 says "a control arm", not "the untreated arm", so any arm the paper itself
designates as its control satisfies P3. Where a paper reports more than one, the arm the
paper designates as its control is used. Where the abstract does not identify the comparator
at all, the record goes to the full-text stage rather than being recorded or dropped.

A literal reading is what the discipline requires. Narrowing P3 to untreated controls now
would be a revision, and it would be a revision made while its effect on the verdict was
already visible — which is the thing the registration exists to prevent.

## The sensitivity analysis that goes with it

Because the literal reading is the more generous one, the restricted version is reported
beside it as an exploratory analysis under section 3.6, and in every table that carries an
`Emax` value:

* **Primary (registered):** `ratio_p` = the median over all qualifying studies, whatever kind
  of control arm each uses.
* **Exploratory:** `ratio_p` recomputed over untreated-, vehicle- or mock-controlled studies
  only, with the count of studies dropped reported.

The exploratory version is the more conservative of the two for the framework's own
hypothesis, since an active control yields a smaller ratio and therefore a smaller `Emax`.
Reporting both means the reader can see how much of the prior rests on the reading of a word.

## The one extraction disagreement, and how it was resolved

PMID 32109370, prodrug_enzyme. The abstract says "median survival was extended from 31 to
46 days". Extractor A recorded nothing, because the 31-day arm is never identified; extractor
B recorded 46 against 31 and flagged it for verification.

Resolved **in A's favour**, by a mechanical rule declared before the comparison was made: a
number is recorded from an abstract only where both extractors read the same number *and* the
abstract identifies both arms. Otherwise the record goes to the full-text stage. Nothing is
dropped by this rule — 32109370 remains in the set, awaiting its full text — so the
conservative resolution costs no evidence.

## What the abstracts gave, and what this is not

8 of the 40 carried records yield a ratio from the abstract alone; **37 of 40 need full text**,
and the eight below are therefore **not** `ratio_p`. Item 11(i.2) defines `ratio_p` as the
median over *every* qualifying study, and this is a fifth of them.

| PMID | Payload | Treated | Control | Unit | Ratio | Control arm |
|---|---|---|---|---|---|---|
| 19920199 | CRAd | 75.5 | 42 | days | 1.798 | to verify |
| 28148846 | prodrug_enzyme | 62 | 32 | days | 1.938 | unnamed in abstract |
| 30257000 | prodrug_enzyme | 41 | 34 | days | 1.206 | "control group" |
| 31909530 | prodrug_enzyme | 25 | 12 | days | 2.083 | unnamed in abstract |
| 32093290 | prodrug_enzyme | 30 | 17.5 | days | 1.714 | vector without prodrug |
| 33463272 | prodrug_enzyme | 67 | 53.5 | days | 1.252 | to verify |
| 26640245 | CAR | 200.5 | 73 | days | 2.747 | parental NK-92 |
| 39095594 | CAR | — | — | ratio | 2.000 | unstated ("doubled") |

**Two things in this table are worth saying out loud now, because they bear on a prediction
the registration already made.**

1. The ratios run from 1.21 to 2.75, against the 1.33 to 1.48 the design document assumed from
   a single ovarian study. Larger ratios mean a larger `Emax` and a smaller hazard ratio, so
   item 11(i.6)'s advance statement — that the rule places T5 below its 3-month bar — holds by
   a narrower margin than it was written with. It is not yet contradicted: the bar needs a
   hazard ratio near 0.59, and these ratios imply roughly 0.65 to 0.76.
2. CAR's two ratios are the **highest** of the three payloads. The design document wanted CAR
   fixed at "the pessimistic end", which read literally was a hazard ratio of 1.0 and no
   benefit at all. On the mechanical rule it is, so far, the strongest payload. This is the
   clearest vindication of replacing the shortcut with a rule, and it is also a warning: the
   preclinical CAR literature and the human CAR literature in glioblastoma point in opposite
   directions, and item 11(i) uses the preclinical one. Both are reported.
