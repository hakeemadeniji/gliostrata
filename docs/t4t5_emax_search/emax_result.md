# Emax, computed under item 11(i)

Computed 2026-11-28T04:54:28Z from `docs/t4t5_emax_search/screening/extraction_final.jsonl`.

> ## PROVISIONAL — not the registered value
>
> **22 of 40 extracted rows have no full text yet**, so their control arms are `uncertain` and they cannot enter the primary set. Item 11(i.2) defines `ratio_p` as the median over *every* qualifying study; this is the median over those qualifying so far.
>
> Nothing here may be entered in `configs/response_surface_v1.yaml`. That is the freeze step of item 11(e) and (m), and it happens once, when the evidence base is closed — either because the outstanding full texts have been obtained or because they have been recorded as unobtainable. See `manual_download_needed.md`.

`Emax_p = f * ln(ratio_p)`, with `f = 0.5` — **assumed**, and declared as assumed wherever a T4 or T5 number is reported. `ratio_p` is the **median** over qualifying studies, never one study chosen by hand.

**Primary** applies the disclosed deviation of item 11(i.7): the control arm must have received no agent intended to have anti-tumour activity. **Exploratory** applies (i.3) as literally written. Where the primary set is empty the payload falls back to the exploratory set, flagged as such.

| Payload | Rule | n | ratio | Emax | implied HR | dropped by narrowing |
|---|---|---|---|---|---|---|
| CRAd | narrowed | 1 | 1.6546 | 0.2518 | 0.7774 | 1 |
| prodrug_enzyme | narrowed | 4 | 1.4162 | 0.1740 | 0.8403 | 0 |
| CAR | literal_fallback | 0 | 3.1674 | 0.5765 | 0.5619 | 2 |

## Both sets side by side

| Payload | Primary n | Primary HR | Exploratory n | Exploratory HR |
|---|---|---|---|---|
| CRAd | 1 | 0.7774 | 2 | 0.7625 |
| prodrug_enzyme | 4 | 0.8403 | 4 | 0.8403 |
| CAR | 0 | — | 2 | 0.5619 |

## How each payload's set was reduced

**CRAd** — 7 carried past stage 1, 2 eligible under P1 to P5, 1 in the primary set.

  - excluded: P2: not confirmed
  - excluded: P3: no treated-and-control median survival obtainable

  - control arms: {'uncertain': 1, 'vehicle': 1}
  - studies used: 21718006

**prodrug_enzyme** — 23 carried past stage 1, 4 eligible under P1 to P5, 4 in the primary set.

  - excluded: P2: not an intracranial rodent glioma model
  - excluded: P2: not confirmed
  - excluded: P3: no treated-and-control median survival obtainable
  - excluded: treated arm is a combination, not the payload alone

  - control arms: {'mock_or_reporter_vector': 3, 'untreated': 1}
  - studies used: 28148846, 32093290, 33463272, 35017150

**CAR** — 10 carried past stage 1, 2 eligible under P1 to P5, 0 in the primary set.

  - excluded: P2: not an intracranial rodent glioma model
  - excluded: P2: not confirmed
  - excluded: P3: no treated-and-control median survival obtainable

  - control arms: {'uncertain': 2}
  - studies used: 26640245, 42573529

## Against the advance expectation of item 11(i.6)

For context only, the 3-month gain bar of item 11(c) sits near a hazard ratio of 0.59. Item 11(i.6) states in advance that the rule is expected to place T5 below that bar.

**On the evidence so far that expectation no longer holds for: CAR.** This is recorded as it stands, not smoothed over: a pre-registered expectation that the data contradicts is a result, and the point of stating it in advance was so that this could be seen.

Three reasons it is premature to call (i.6) refuted, each of which must be reported alongside it:

1. The evidence base is open — see the provisional notice above.
2. A payload on `literal_fallback` has an **empty** primary set, so its figure rests on studies whose control arms could not be classified.
3. An implied hazard ratio is not a verdict. T5 turns on a 0.59 hazard ratio being reached *at the configuration chosen for a patient*, after the toxicity term of item 11(k) is subtracted — and CAR carries the largest toxicity coefficient of the three payloads. Nothing here is compared against the bar to decide anything.

T5's verdict comes from the registry in section 6, on the sealed cohort, after the freeze. It is not computable yet and nothing above anticipates it.
