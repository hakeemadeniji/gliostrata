# Emax, computed under item 11(i)

Computed 2026-11-28T06:25:32Z from `docs/t4t5_emax_search/screening/extraction_final.jsonl`.

> ## The evidence base is CLOSED
>
> Every obtainable article has been extracted. 2 row(s) are the articles item 11(i.8) records as unobtainable; they fail P3, and are counted and reported here rather than quietly dropped.
>
> These are therefore the values for the freeze amendment of item 11(e) and (m). Entering them in `configs/response_surface_v1.yaml` remains a deliberate step, not a side effect of this computation.

`Emax_p = f * ln(ratio_p)`, with `f = 0.5` — **assumed**, and declared as assumed wherever a T4 or T5 number is reported. `ratio_p` is the **median** over qualifying studies, never one study chosen by hand.

**Primary** applies the disclosed deviation of item 11(i.7): the control arm must have received no agent intended to have anti-tumour activity. **Exploratory** applies (i.3) as literally written. Where the primary set is empty the payload falls back to the exploratory set, flagged as such.

| Payload | Rule | n | ratio | Emax | implied HR | dropped by narrowing |
|---|---|---|---|---|---|---|
| CRAd | narrowed | 3 | 1.1250 | 0.0589 | 0.9428 | 1 |
| prodrug_enzyme | narrowed | 10 | 1.5400 | 0.2159 | 0.8058 | 2 |
| CAR | narrowed | 1 | 2.6913 | 0.4950 | 0.6096 | 3 |

## Both sets side by side

| Payload | Primary n | Primary HR | Exploratory n | Exploratory HR |
|---|---|---|---|---|
| CRAd | 3 | 0.9428 | 4 | 0.8483 |
| prodrug_enzyme | 10 | 0.8058 | 12 | 0.8187 |
| CAR | 1 | 0.6096 | 4 | 0.6529 |

## How each payload's set was reduced

**CRAd** — 7 carried past stage 1, 4 eligible under P1 to P5, 3 in the primary set.

  - excluded: P3: no treated-and-control median survival obtainable

  - control arms: {'uncertain': 1, 'vehicle': 3}
  - studies used: 21629227, 21718006, 23821758

**prodrug_enzyme** — 23 carried past stage 1, 12 eligible under P1 to P5, 10 in the primary set.

  - excluded: P2: not an intracranial rodent glioma model
  - excluded: P2: not confirmed
  - excluded: P3: no treated-and-control median survival obtainable
  - excluded: treated arm is a combination, not the payload alone

  - control arms: {'mock_or_reporter_vector': 3, 'uncertain': 2, 'untreated': 3, 'vehicle': 4}
  - studies used: 17409449, 22996231, 28148846, 29660022, 30257000, 31909530, 32093290, 32439343, 33463272, 35017150

**CAR** — 10 carried past stage 1, 4 eligible under P1 to P5, 1 in the primary set.

  - excluded: P2: not an intracranial rodent glioma model
  - excluded: P3: no treated-and-control median survival obtainable

  - control arms: {'uncertain': 3, 'vehicle': 1}
  - studies used: 26640245

## Against the advance expectation of item 11(i.6)

For context only, the 3-month gain bar of item 11(c) sits near a hazard ratio of 0.59. Item 11(i.6) states in advance that the rule is expected to place T5 below that bar.

On the evidence so far that expectation holds for every payload.

T5's verdict comes from the registry in section 6, on the sealed cohort, after the freeze. It is not computable yet and nothing above anticipates it.
