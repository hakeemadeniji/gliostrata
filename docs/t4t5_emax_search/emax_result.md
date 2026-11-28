# Emax, computed under item 11(i)

Computed 2026-11-28T06:19:11Z from `docs/t4t5_emax_search/screening/extraction_final.jsonl`.

> ## PROVISIONAL — not the registered value
>
> **5 of 40 extracted rows have no full text yet**, so their control arms are `uncertain` and they cannot enter the primary set. Item 11(i.2) defines `ratio_p` as the median over *every* qualifying study; this is the median over those qualifying so far.
>
> Nothing here may be entered in `configs/response_surface_v1.yaml`. That is the freeze step of item 11(e) and (m), and it happens once, when the evidence base is closed — either because the outstanding full texts have been obtained or because they have been recorded as unobtainable. See `manual_download_needed.md`.

`Emax_p = f * ln(ratio_p)`, with `f = 0.5` — **assumed**, and declared as assumed wherever a T4 or T5 number is reported. `ratio_p` is the **median** over qualifying studies, never one study chosen by hand.

**Primary** applies the disclosed deviation of item 11(i.7): the control arm must have received no agent intended to have anti-tumour activity. **Exploratory** applies (i.3) as literally written. Where the primary set is empty the payload falls back to the exploratory set, flagged as such.

| Payload | Rule | n | ratio | Emax | implied HR | dropped by narrowing |
|---|---|---|---|---|---|---|
| CRAd | narrowed | 3 | 1.1250 | 0.0589 | 0.9428 | 1 |
| prodrug_enzyme | narrowed | 9 | 1.5000 | 0.2027 | 0.8165 | 1 |
| CAR | narrowed | 1 | 2.6913 | 0.4950 | 0.6096 | 2 |

## Both sets side by side

| Payload | Primary n | Primary HR | Exploratory n | Exploratory HR |
|---|---|---|---|---|
| CRAd | 3 | 0.9428 | 4 | 0.8483 |
| prodrug_enzyme | 9 | 0.8165 | 10 | 0.8187 |
| CAR | 1 | 0.6096 | 3 | 0.6096 |

## How each payload's set was reduced

**CRAd** — 7 carried past stage 1, 4 eligible under P1 to P5, 3 in the primary set.

  - excluded: P3: no treated-and-control median survival obtainable

  - control arms: {'uncertain': 1, 'vehicle': 3}
  - studies used: 21629227, 21718006, 23821758

**prodrug_enzyme** — 23 carried past stage 1, 10 eligible under P1 to P5, 9 in the primary set.

  - excluded: P2: not an intracranial rodent glioma model
  - excluded: P2: not confirmed
  - excluded: P3: no treated-and-control median survival obtainable
  - excluded: treated arm is a combination, not the payload alone

  - control arms: {'mock_or_reporter_vector': 3, 'uncertain': 1, 'untreated': 3, 'vehicle': 3}
  - studies used: 17409449, 28148846, 29660022, 30257000, 31909530, 32093290, 32439343, 33463272, 35017150

**CAR** — 10 carried past stage 1, 3 eligible under P1 to P5, 1 in the primary set.

  - excluded: P2: not an intracranial rodent glioma model
  - excluded: P2: not confirmed
  - excluded: P3: no treated-and-control median survival obtainable

  - control arms: {'uncertain': 2, 'vehicle': 1}
  - studies used: 26640245

## Against the advance expectation of item 11(i.6)

For context only, the 3-month gain bar of item 11(c) sits near a hazard ratio of 0.59. Item 11(i.6) states in advance that the rule is expected to place T5 below that bar.

On the evidence so far that expectation holds for every payload.

T5's verdict comes from the registry in section 6, on the sealed cohort, after the freeze. It is not computable yet and nothing above anticipates it.
