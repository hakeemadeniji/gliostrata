# T4/T5 sensitivity analyses: the two-tier split

**Status.** A working plan, fixed 2026-11-28, **before any sensitivity analysis has been
run** and before the T4/T5 read. It implements section 5 of `docs/t4t5_model_design.md`,
which registered *what* to vary but not *how*, and that omission turns out to decide whether
the plan takes a minute or a quarter of a year.

**This is not registered text and needs no amendment.** Item 11 adopted Design A's model;
section 3.6 registers that exploratory analyses exist, run only on the data of the registered
read, and never change a verdict. Neither registers the sensitivity plan's internals. It is
written down anyway, in advance, so that the split is a design choice rather than a retrofit
around a slow script.

## 1. The problem the design document left open

Section 5 asks for a 1000-draw Latin hypercube and a tornado over ~20 assumed parameters,
with outputs including T5's success and harm fractions and T4's distinct-config count and
regret reduction. It never says whether a draw **re-scores** the configuration already chosen
or **re-optimises** to find a new one. Measured on this machine, for the 268 patients of a
full run (T4's fifty plus T5's 218):

| | Re-score | Re-optimise | Ratio |
|---|---|---|---|
| One patient, one setting | 148 µs | 33.3 s | **225,000×** |
| Tornado, ~40 settings | 1.6 s | 99 h | |
| Latin hypercube, 1000 draws | **0.7 min** | **2,480 h = 103 days** | |

So the same sentence describes a minute's work or a quarter of a year. The cost is not the
objective, which really is microseconds; it is the Gaussian-process surrogate refitted at each
of the optimiser's 48 iterations.

## 2. The two questions, which are both legitimate

The split is not a compromise for speed. The two costs answer genuinely different questions,
and a plan that conflated them would have answered only one of them whatever it cost:

* **Tier 1 — "is the recommended therapy's benefit robust?"** Hold each patient's
  configuration fixed at what the central parameters chose, and re-score it under each draw.
  This is the question a clinician has: *that* therapy is what the patient would receive, and
  its predicted benefit should not evaporate under a different assumption about, say, the
  acting fraction.
* **Tier 2 — "is the recommendation itself robust?"** Re-optimise, and see whether a different
  configuration wins. This is the question a methodologist has, and it is the only way to get
  T4's two metrics at all: a distinct-configuration count and a regret reduction are
  properties of the optimiser's behaviour and cannot be re-scored.

Tier 1 cannot produce T4's metrics. Tier 2 cannot afford 1000 draws. Each covers what the
other cannot.

## 3. The split, item by item

Section 5's items map as follows. "Outputs" names only what that tier can honestly produce.

| § 5 item | Tier | Settings | Outputs | Cost |
|---|---|---|---|---|
| 1. One-at-a-time tornado | 1 | ~40 (each assumed parameter at both ends of its range) | T5 success, T5 harm, mean gain | ~2 s |
| 2. Global Latin hypercube | 1 | 1000 joint draws, seed 42 | the same, plus Spearman rank indices and the fraction of draws with success ≥ 0.50 and harm ≤ 0.10 | ~1 min |
| 3. Break-even effect size | 1 then 2 | bisection on a common Emax multiplier, then **one** re-optimised confirmation at the multiplier found | the hazard ratio the therapy would need | ~1 s, then ~30 min |
| 4. Prior regimes | 2 | 3 (skeptical as registered, moderate HR 0.80, optimistic HR 0.60) | all four metrics | ~90 min |
| 5. Structural variants | 2 | ~6 (Hill vs linear; combo Bliss vs additive vs best-of-two; modifiers on/off; coverage function; survival variant 2 vs 3) | all four metrics | ~3 h |
| 6. §3.6 items | 1 and 2 | subgroups and complete-case in tier 1; two alternative optimiser seeds in tier 2 | as the tier allows | ~1 h |
| 7. Reporting | — | — | one table per output, verdict beside each | — |

**Why item 5 is tier 2 and item 1 is tier 1.** A structural variant changes the *shape* of
the response surface, so the configuration that maximises it can move; re-scoring a fixed
configuration under a different functional form would silently answer a question nobody
asked. A one-at-a-time parameter shift mostly rescales the surface, so the argmax usually
holds — "usually" is why tier 2 also takes the parameters tier 1 shows matter most.

## 4. Tier 2's patient sample, and why it is not the full 268

Tier 2 runs on **the T4 fifty** — the first fifty by sorted patient id in the sealed
partition, which item 7 already registers as T4's cohort — rather than on all 218. At 33.3 s
per patient that is about 28 minutes per setting instead of two hours, so twelve settings fit
in roughly five and a half hours of unattended background time rather than a full day.

This is a declared limitation, not a hidden one: tier 2's T5 fractions are computed on fifty
patients and are therefore noisier than the primary's 218, and they are reported with that n
beside them. The primary result is unaffected — tier 2 is exploratory and never changes a
verdict (section 3.6).

The alternative, running tier 2 on all 218, costs about 24 hours for twelve settings. That is
affordable if the result turns out to hinge on it, and the decision to stay at fifty is
revisited if tier 1 shows the outputs moving sharply within the ranges.

## 5. What is fixed in advance, so it cannot be chosen later

1. **Tier 1's reference configuration** is whatever the primary registered run chose for that
   patient. It is never re-chosen, and tier 1 never reports a T4 metric.
2. **Tier 2's twelve settings** are: the three prior regimes, the six structural variants, the
   break-even confirmation, and the two parameters with the largest tier-1 Spearman index.
   The last two are data-dependent by construction, and that is stated in the report; the
   *rule* for picking them is fixed here, before any of it has run.
3. **Seeds.** Latin hypercube seed 42, bootstrap seed 42, optimiser seeds as in the primary
   run. The two alternative optimiser seeds of §5 item 6 are 43 and 44.
4. **Every output is labelled exploratory** in every table and caption, and none changes a
   verdict. Section 3.6 scope holds: these use the data of the registered read and nothing
   further is read for them.
5. **Both tiers report the saturation of T4's regret metric** where they report it at all
   (item 5, recorded 2026-11-28), since a sensitivity analysis of a metric pinned at 100% is
   uninformative by construction and should say so rather than appear to confirm stability.

## 6. Order of operations

Tier 1 runs first, in full, because it is a minute's work and it chooses two of tier 2's
settings. Tier 2 then runs as a single unattended batch. Both run **after** the registered
T4/T5 read, on its data, and neither can begin before the freeze.

The code can be written and tested on synthetic patients now, which is the point of fixing
this plan today rather than on the day the DOI arrives.
