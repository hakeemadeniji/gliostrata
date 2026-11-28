# Part 3. Methods: a registered model for configuration choice, and the evidence behind it

Status: working draft, not for circulation. Companion to `draft_v0.md`, whose Parts 1 and 2
cover the comparator search and the external-validation benchmark. Every number below is
traced to a repository file; places needing the lead's decision are marked [LEAD: ...].

**Standing caveat, stronger here than in Parts 1 and 2.** Everything in this part is
*model-conditional*. It predicts restricted mean survival under a therapy configuration using
a survival model that **failed its own external validation** (Part 2: C-index 0.5385, 95% CI
0.4994–0.5807, against a registered threshold of 0.65), and a dose-response whose efficacy
parameters are translated from rodent experiments by an assumed factor. No quantity in this
part is a survival benefit, and the registration (§2.2 item 6) forbids calling one that. What
it can support is narrower and is stated in §3.7.

---

## 3.1 What the model is, and what it is for

The framework asks whether a per-patient *choice* among stem-cell-delivered therapy
configurations can be made better than a fixed one. Answering that needs a function from
(patient, configuration) to a predicted outcome. No such function exists in the literature, so
one was specified, registered and frozen before any held-out patient was scored.

The configuration space is registered as coded (§2.2 item 4): carrier cell type (iPSC-NSC,
MSC, HB1_F3_line) × payload (CRAd, prodrug/enzyme, CAR, combo) × dose (10⁶–5×10⁷ cells,
continuous) × injection sites (1–8) × schedule (single, ×2, ×3), with two carrier–payload
pairs excluded as incompatible, leaving 10 pairs and 240 discrete combinations.

The model has three layers, adopted as "Design A" of `docs/t4t5_model_design.md` and fixed by
§2.2 item 11(g)–(o):

1. a **survival model** supplying each patient's baseline curve under standard therapy;
2. a **delivery model** turning a configuration into an engrafted, acting cell count;
3. a **response surface** turning that into a log-hazard offset, against a toxicity term.

## 3.2 Survival model and the Stupp baseline

A Cox proportional-hazards model (Breslow ties) fitted on TCGA-GBM with age, sex and treatment
arm — the "transferable" variant of Part 2, chosen over a variant additionally carrying IDH and
MGMT status. **This was a deliberate choice against discrimination.** The molecular variant
scores a higher external C-index (0.5915 vs 0.5385), but requires an IDH and an MGMT call for
every patient, which the external partition's registered eligibility rule never demanded; on
CGGA the equivalent complete-case cohort falls from 283 to 241. Using it would have scored the
end-to-end test on fewer than the 218 patients its scoring set registers. A 0.05 gain in
C-index, on a model that fails the threshold either way, did not justify shrinking a sealed
population (§2.2 item 11(g)).

A patient's **Stupp baseline** is this model's prediction with the arm set to radiotherapy plus
temozolomide and no stem-cell configuration (item 11(a)). Predicted survival is evaluated on a
fixed 0.1-month grid from 0 to τ = 24 months, and restricted mean survival time (RMST) is the
trapezoidal area under that curve. A configuration's **gain** is its RMST minus the baseline's.

Neither IDH nor MGMT appears anywhere in the frozen model — not in the baseline, and not as an
effect modifier, because no sourced modifier for the configuration effect exists and none was
invented. The consequence is disclosed rather than buried: **IDH-mutant patients' baselines are
modelled as if IDH-wild-type** (item 11(l)). IDH-mutant tumours are not glioblastoma under WHO
2021, which is an argument for excluding those patients; this was not acted on, because the
partition is sealed and changing its membership after sealing is what sealing prevents. IDH
status is reported as a pre-registered exploratory subgroup (§3.6).

## 3.3 Delivery: a closed-form model, declared uncalibrated

For carrier *k*, a per-patient tropism score *z* is the mean standardised expression of that
carrier's chemoattractant panel — CXCL12, VEGFA, PLAU and CCL2 for the neural stem cell
carriers, TGFB1 for mesenchymal stem cells.

**Standardisation uses a frozen TCGA reference, never a statistic of the cohort being scored.**
A per-cohort z-score would let the sealed external cohort leak into its own transform. The
reference is TCGA-GBM, 284 samples, log2(FPKM+1) from GDC's `fpkm_unstranded`, checksummed in
§8 at the freeze.

The scale matters more than it looks. The reference was first built on log2(CPM+1); measured on
the same 284 samples, that leaves a systematic offset of −1.36 SD for the NSC panel and −1.47 SD
for TGFB1 against an FPKM-normalised cohort, which alone would reduce every patient's engrafted
dose by about a third. Rebuilding on FPKM removes the gene-length component. **A second-order
difference remains and is not corrected**: TCGA's FPKM here is GDC's STAR-derived, CGGA's is
RSEM.

The engrafted, acting count is *n*<sub>eff</sub> = *N* · *f*<sub>k</sub> ·
logistic(*a*<sub>k</sub> + *g*<sub>k</sub>*z*), with spatial coverage 1 − exp(−*n*/*n*₀) over
injection sites and a multiplicative schedule factor.

This layer is **declared uncalibrated**, and that declaration is unconditional: it does not
depend on whether the project's separate PDX-calibration row (T3) ever obtains data (item
11(h)). *f*, *a*, *g* and *n*₀ are assumed, each with a sensitivity range. No PDE or
agent-based simulator was built.

## 3.4 Response surface and toxicity

The log-hazard reduction for payload *p* is

> *E* = *E*<sub>max,p</sub> · *m*<sub>p</sub>(*i*) · σ(*s*) · cov(*n*) ·
> *n*<sub>eff</sub><sup>h</sup> / (*n*<sub>eff</sub><sup>h</sup> + ED50<sup>h</sup>)

with a patient modifier *m* = exp(**w**′**u**) over a standardised immune-microenvironment score
and mesenchymal-like state fraction. The prodrug payloads receive no MGMT modifier: 5-FU and
SN-38 act through antimetabolite and topoisomerase mechanisms, not through MGMT — a mechanistic
argument, not a measured effect.

**"Combo" is defined per carrier**, as the two payloads that carrier is registered as able to
carry, combining by Bliss independence on the hazard scale with an assumed synergy of 1.0. A
single uniform pair was not available: the registered exclusions forbid CRAd on MSC and CAR on
the HB1.F3 line, so any uniform definition would place an excluded payload on a carrier. An
optimiser-chosen pair was also rejected, because two configurations sharing carrier, payload and
site count would then be different therapies, which the registered distinctness rule — and the
personalisation metric built on it — depend on them not being (item 11(j)).

Toxicity is an excess log hazard, *T* = τ<sub>p</sub>(*N*/5×10⁷)<sup>q</sup>(1 + η(*n*−1)) +
δ<sub>delay</sub>, covering dose, multiple injections and the cost of delaying standard therapy.
It is in the primary surface, not an exploratory variant, because **it is the model's only
source of harm**: without it the end-to-end safety threshold could not fail, and a safety
threshold that cannot fail is not a safety threshold (item 11(k)). A larger prior drawn from
reported CAR neurotoxicity was rejected as not transferable to carrier-cell delivery and is a
sensitivity range instead.

The configuration's curve is *S*(*t*) = *S*<sub>Stupp</sub>(*t*)<sup>exp(−E+T)</sup>.

## 3.5 The efficacy prior: a registered literature search

*E*<sub>max</sub> is the single quantity on which the end-to-end result most depends, and the
only one with no internal evidence. It was fixed by a **rule declared before the search that
feeds it**, so that no study could be chosen for the number it would produce:

> *E*<sub>max,p</sub> = 0.5 × ln(*ratio*<sub>p</sub>)

where *ratio*<sub>p</sub> is the **median** treated-to-control median-survival ratio across
*every* qualifying study for that payload, and 0.5 is an assumed mouse-to-human translation
factor, declared as assumed wherever a number from it is reported.

An earlier version named three studies carried over from the design document. It was rewritten
before the search because none satisfied its own criterion: one was in ovarian cancer, one
reported a tumour-volume rather than survival ratio, and the third had been shortcut to "the
pessimistic end" of an assumed range.

**Search and screening.** Three PubMed queries, one per payload, written into the registration
and run once (2026-11-28), no year limit: 108 records. Eligibility was fixed in advance (P1
peer-reviewed and in vivo; P2 orthotopic or intracranial rodent glioma; P3 median survival for
both a control and a payload arm; P4 payload class matched; P5 carrier delivery neither required
nor disqualifying). Titles and abstracts were screened **twice independently**, agreement
95.4% with no include/exclude conflict, disagreements resolved by a rule fixed before they were
seen; 40 records carried forward. Full texts were obtained for 38 of 40 — 18 from PubMed
Central, 20 retrieved manually — and extracted twice independently. Two articles could not be
obtained, fail P3 for inaccessibility, and are reported rather than dropped.

**A disclosed deviation.** P3 required a "control arm" without defining one, and the qualifying
studies use untreated, vehicle, reporter-vector and active-therapeutic comparators, which give
systematically different ratios. The primary set was narrowed to studies whose control arm
received no agent intended to have anti-tumour activity. This is recorded as a **deviation, not
an amendment**, because the reading was chosen after the literal reading's ratios were visible;
the inclusion test is structural — what the control arm *received* — and never its result. The
literal reading is reported alongside throughout.

**Result.**

| Payload | *ratio* | *E*<sub>max</sub> | implied HR | n studies |
|---|---|---|---|---|
| CRAd | 1.125 | 0.0589 | 0.943 | 3 |
| Prodrug/enzyme | 1.540 | 0.2159 | 0.806 | 10 |
| CAR | 2.691 | 0.4950 | 0.610 | 1 |

Two features of this table are reported with it rather than left to a limitations section.
**CRAd's modelled effect is close to none**, consistent with null human randomised data for
oncolytic adenovirus. **CAR's prior rests on a single study** and is the largest of the three;
its other eligible studies were excluded by the intracranial-model criterion, by unobtainable
medians, or by irrelevant-specificity CAR controls. No minimum study count was introduced,
because such a rule would be payload-neutral in wording and CAR-only in effect, and would have
been written with CAR's count and its *E*<sub>max</sub> both already visible. The registered
reading is, for every thin payload, **the most favourable of the available readings** — every
broader reading moves it further from the success threshold.

## 3.6 Optimisation, metrics and controls

Per patient, Bayesian optimisation (Gaussian process with expected improvement) and random
search each receive a registered budget of **60 objective evaluations** — 12 Sobol initial
points and 48 model-guided iterations for the former. Four metrics are evaluated against
pre-registered thresholds through a checksummed registry, which returns PASS, FAIL,
INCONCLUSIVE or NOT EVALUABLE and refuses a PASS when an auxiliary requirement is not reported:

* distinct configurations in the top five across fifty patients (≥ 3);
* regret reduction of the optimiser against random search (≥ 25%, with a Wilcoxon signed-rank
  requirement);
* fraction of held-out patients with a gain of at least 3 months (≥ 0.50);
* fraction with a gain below zero (≤ 0.10, with a CI-upper requirement of 0.15).

The fifty are the first by sorted patient identifier within the sealed partition, so the
selection cannot depend on an outcome.

**A property of the regret metric, recorded before the read.** Regret is measured against the
best value found by *either* method, so a method attaining that best for every patient has a
regret of identically zero, the metric reads exactly 100% whatever its margin, its bootstrap
interval is degenerate, and the threshold **cannot fail**. In the pre-read rehearsal this is
what happened. The metric is unchanged — rewriting it once its behaviour had been measured
would be worse than reporting that it saturates — and it is reported together with whether it
saturated and with the mean advantage in months, which does not saturate. [LEAD: this paragraph
belongs in Results as well if the confirmatory run saturates too.]

**Four controls run before any held-out patient is scored** and a failure stops the read: a
no-modifier control, which must collapse the optimiser's answer to a single configuration; a
null-payload control, which must give a benefit term of exactly zero and a gain equal to the
toxicity term's alone; a synthetic end-to-end run of the evaluation script itself; and a
rehearsal of the whole pipeline on a held-out development split. All four passed. The
null-payload control gave a benefit of exactly 0 across 120 configurations and 20 patients and
matched the toxicity-only gain to within 10⁻¹² in all 2,400 pairs — and **could not have passed
as the control was first written**, which required a gain of exactly zero, unreachable while the
toxicity term remains in the primary surface.

## 3.7 What this model can and cannot support

Stated here rather than in a discussion section, because it bounds every number in Part 3's
results.

**It cannot** support a claim of survival benefit, of any size, for any configuration. It
predicts RMST under a model whose discrimination failed external validation, using an efficacy
prior translated from rodent experiments by an assumed factor, on a delivery layer declared
uncalibrated.

**It can** support claims of the form: *given* this registered model, does per-patient choice
differ across patients, does model-guided search beat random search, and what effect size would
the therapy need for the end-to-end criterion to be met? The last of these — the break-even
effect size — is the most transferable output, because it converts a model-conditional result
into a statement about the effect a carrier-delivered therapy would have to achieve, which can
be set against what oncolytic and CAR therapies have achieved in humans.

[LEAD: consider whether Part 3 should be a separate paper from Parts 1 and 2 or a final section
of one. The honest framing of Part 3 is a *negative or model-conditional* result, and pairing it
with Part 2's benchmark may serve both better than a standalone.]

## 3.8 Reproducibility

The model, its parameters and the evidence behind them are files, not descriptions: a single
versioned parameter file carrying every numeric value with its status (fitted, published with a
PMID, or assumed) and its sensitivity range; a frozen standardisation reference; the full search
and screening log, with every screening decision and the passage behind every extracted number;
and an append-only access ledger recording every read of the external cohort. The parameter file
and the reference are checksummed in the registration's amendment log at the freeze, which
precedes any held-out read. [LEAD: add the DOI, the freeze date and the two checksums here once
recorded.]
