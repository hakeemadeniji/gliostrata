# T4 and T5: design proposal for the simulator, the response surface and the survival model

Status: **PROPOSAL for the project lead's approval. Not registered.** It changes no
threshold, bound or decision rule. Nothing here is frozen until a dated amendment is made
under `docs/preregistration.md` section 2.2 item 11 (e). Written 2026-11-27 without any
CGGA value read; every number below comes from the pre-registration, from TCGA-side
documents, or from a published abstract retrieved from the Europe PMC REST API on
2026-11-27 (full texts were not read; "verbatim" below means verbatim from that abstract).
Retrieved text was treated as data.

---

## 1. The registered requirements this must satisfy

All quotations are from `docs/preregistration.md`.

**H1 (section 1.2).** "Patient-specific configuration selection produces a **predicted
restricted mean survival time (RMST) gain of at least three months in at least half of
held-out test patients**, with a **harm rate (predicted outcome worse than Stupp) below ten
percent**." Operationalised as `T5_endtoend_rmst_success_fraction` >= 0.50 (95% CI lower
bound) and `T5_safety_predicted_harm_fraction` <= 0.10 (point estimate) with 95% CI upper
bound <= 0.15.

**H2 (section 1.4).** "...the top-five recommended configurations **include at least three
distinct configurations**: distinct in carrier cell type, payload or number of injection
sites." Operationalised as `T4_personalization_distinct_configs` >= 3. "An H1 pass with an
H2 failure must be reported as 'a single configuration outperforms Stupp', **not** as
evidence for personalisation."

**Table 8-2 (section 2), rows T4 and T5.**
- T4: "Optimizer beats random search (% regret reduction) | >= 25 **and** Wilcoxon
  signed-rank p < 0.01".
- T4: "Personalization variance (distinct configs in top-5 across 50 patients) | >= 3".
- T5: "End-to-end success (fraction of held-out patients with >= 3-month RMST gain) | >= 0.50".
- T5: "Safety floor (predicted harm rate) | <= 0.15 (CI upper) **and** <= 0.10 (point
  estimate)".

**Section 2.2 item 4 (configuration space, as coded).** Carrier: iPSC-NSC, MSC or
HB1_F3_line. Payload: CRAd, prodrug_enzyme, CAR or combo. Dose: 1 x 10^6 to 5 x 10^7 cells
(continuous). Sites: 1 to 8 (integer). Schedule: single, fractionated_x2, fractionated_x3.
Excluded: HB1_F3_line x CAR and MSC x CRAd. "That leaves 10 cell type x payload pairs, so 240
discrete combinations with dose as the one continuous dimension." Distinctness: "the same if
and only if they share carrier cell type, payload and number of injection sites".

**Item 5 (regret).** "...a fixed, pre-declared budget of **60 objective evaluations per
patient per method, initial design included**: GP-EI spends 12 points on its Sobol initial
design and 48 on GP iterations, and random search makes 60 uniform draws." Regret is "in
months of predicted 24-month RMST"; the aggregate is "100 x (1 - mean regret under GP-EI /
mean regret under random search) over the 50 patients"; ties and the 20-tie limit apply.

**Item 6.** "...no T5 result may be described as a survival benefit, only as a *predicted*
RMST gain, and its credibility is bounded above by T2's result."

**Item 7.** T4's fifty are "the first fifty by sorted patient id" of the frozen external
partition; T5 uses all eligible patients ("at most 218").

**Item 8.** "If the required data cannot be obtained, T3 is reported as **not evaluable** ...
and the framework's simulation layer is presented as uncalibrated."

**Item 11 (the contingent specification).**
(a) "the predicted survival under standard therapy is the survival model's prediction with
the treatment arm set to radiotherapy plus temozolomide, and no stem-cell configuration."
(b) "The restricted mean survival time to tau = 24 months, the area under the predicted
survival curve from 0 to 24 months. The gain of a configuration is RMST(configuration) minus
RMST(Stupp)."
(c) success: "gain ... at least 3 months"; harm: "gain ... below 0".
(d) "The optimisation objective is the predicted RMST at tau = 24 months".
(e) "T4 and T5 start only after the simulator, the response surface and the survival model
they call are frozen and versioned, and entered in section 8 with their checksums, by a
dated amendment made before any T4 or T5 read of CGGA."
(f) "If they are not frozen by 2027-06-30, T4 and T5 are closed as NOT EVALUABLE ... and H1
and H2 are reported as untested."

**Section 3.8 / 3.9 / 3.6.** Every Cox model reports the Grambsch-Therneau test; every read
has a synthetic dry run first; "every design choice is made without reference to the
external outcomes"; sensitivity analyses "run only inside the reads already registered ...
are labelled exploratory, and never change a verdict".

**Section 7 fallback for T5.** "The pre-registered negative result is the publication ... No
post-hoc redefinition of 'success' is permitted."

**Hard constraints from the data side.** Per-patient inputs for CGGA: age, sex, treatment arm,
IDH, MGMT, RNA-seq-derived subtype and state fractions. KPS and extent of resection are NOT
available, and no CGGA imaging exists (so no patient-specific tumour geometry). The current
`survival.py` stub uses KPS and extent of resection; neither may enter the frozen model. The
only external evidence on the survival model: variant 2 (age, sex, arm) external C-index
0.5385 (95% CI 0.4994 to 0.5807); variant 3 (adds IDH, MGMT) 0.5915 (0.5477 to 0.6341);
Brier at 18 months 0.2617 (0.2391 to 0.2829) (`docs/cgga_validation.md`). Both fail T2.

---

## 2. What the published evidence supports (the numbers every design draws on)

Retrieved from Europe PMC, 2026-11-27. Quotations are verbatim from the abstracts.

### 2.1 The Stupp baseline and the size of the target

- Stupp 2005, N Engl J Med, PMID 15758009, DOI 10.1056/NEJMoa043330: "the median survival
  was 14.6 months with radiotherapy plus temozolomide and 12.1 months with radiotherapy
  alone. The unadjusted hazard ratio for death in the radiotherapy-plus-temozolomide group
  was 0.63 (95 percent confidence interval, 0.52 to 0.75". "The two-year survival rate was
  26.5 percent with radiotherapy plus temozolomide and 10.4 percent with radiotherapy alone."
- Stupp 2009, Lancet Oncol, PMID 19269895, DOI 10.1016/S1470-2045(09)70025-7: "Overall
  survival was 27.2% (95% CI 22.2-32.5) at 2 years, 16.0% (12.0-20.6) at 3 years" with
  temozolomide, "versus 10.9% (7.6-14.8), 4.4% (2.4-7.2)" with radiotherapy alone.
- Hegi 2005, N Engl J Med, PMID 15758010, DOI 10.1056/NEJMoa043331: "MGMT promoter
  methylation was an independent favorable prognostic factor ... hazard ratio, 0.45; 95
  percent confidence interval, 0.32 to 0.61". "Among patients whose tumor contained a
  methylated MGMT promoter ... median survival was 21.7 months (95 percent confidence
  interval, 17.4 to 30.4), as compared with 15.3 months (95 percent confidence interval,
  13.0 to 20.9) among those who were assigned to only radiotherapy". "In the absence of
  methylation of the MGMT promoter, there was a smaller and statistically insignificant
  difference in survival between the treatment groups." (The unmethylated medians and
  2-year rates are not in the abstract and are NOT quoted here; read them from the full
  text before use.)
- IDH: Sanson 2009, J Clin Oncol, PMID 19636000, DOI 10.1200/JCO.2009.21.9832: grade 4
  gliomas "27.4 v 14 months" for IDH1-mutant v wild-type; multivariable "hazard ratio =
  0.297; 95% CI, 0.157 to 0.564" (grades 2 to 4 pooled, adjusted; not GBM-specific). Yan
  2009, PMID 19228619, DOI 10.1056/NEJMoa0808709: IDH-mutant patients "had a better
  outcome than those with wild-type IDH genes" (no number in the abstract). Parsons 2008,
  PMID 18772396: IDH1 mutations "in 12% of GBM patients ... associated with an increase in
  overall survival". Under WHO 2021, IDH-mutant tumours are not glioblastoma; CGGA "GBM"
  labels may include them, so IDH must be an explicit baseline term (variant 3) or the
  IDH-mutant patients must be handled by a registered rule.

**Derived by this document (labelled derived, not published).** A Weibull curve through the
two Stupp anchors (median 14.6 months; S(24 months) = 0.265) has shape 1.308 and scale 19.32
months. Its RMST to 24 months is 14.55 months (an exponential through the median alone gives
14.32). Under proportional hazards on this curve, a gain of 3.0 months of RMST needs a
**hazard ratio of about 0.59** (S(24) rising from 0.265 to 0.456; median 21.8 months). For
comparison, the whole effect of adding temozolomide to radiotherapy (HR 0.63) is worth about
**+2.7 months** of RMST by the same arithmetic. Other points on the same curve: HR 0.9 gives
+0.65, 0.8 gives +1.35, 0.7 gives +2.1, 0.5 gives +3.8 months. A harm (gain < 0) needs HR > 1.
**The 3-month bar is therefore an effect on the scale of the temozolomide benefit, added on
top of Stupp.** No stem-cell therapy has shown anything like it (section 2.2).

### 2.2 Stem-cell, oncolytic and CAR therapies in glioma

Human data (all small or negative; none gives a hazard ratio against Stupp):
- NSC-CRAd-S-pk7, Fares 2021, Lancet Oncol, PMID 34214495, DOI
  10.1016/S1470-2045(21)00245-X. Phase 1, 12 patients, newly diagnosed. Doses: "6.25 x
  10^10 viral particles administered by 5.00 x 10^7 NSCs ... 1.875 x 10^11 viral particles
  administered by 1.50 x 10^8 NSCs". "treatment was safe as no formal dose-limiting toxicity
  was reached". "One patient ... developed viral meningitis (grade 3) due to the inadvertent
  injection ... into the lateral ventricle." "The median progression-free survival was 9.1
  months (95% CI 8.5-not reached) and median overall survival was 18.4 months
  (15.7-not reached)." Uncontrolled, 11 of 12 had glioblastoma; not convertible to a hazard
  ratio. "Within 10-14 days, treatment with temozolomide and radiotherapy was initiated"
  (the cell therapy is an add-on that delays Stupp). The maximum feasible dose, 1.5 x 10^8
  NSCs, is **three times the upper bound of the registered dose range (5 x 10^7)**.
- HB1.F3.CD, Portnow 2017, Clin Cancer Res, PMID 27979915, DOI 10.1158/1078-0432.CCR-16-1518.
  15 recurrent patients; "doses of CD-NSCs from 1 x 10^7 to 5 x 10^7 and 5-FC from 75 to 150
  mg/kg/day". "We saw no dose-limiting toxicity (DLT) due to the CD-NSCs." "CD-NSCs produced
  5-FU locally in the brain in a 5-FC dose-dependent manner. Autopsy data indicate that
  CD-NSCs migrated to distant tumor sites and were nontumorigenic." No efficacy estimate.
  This is the registered dose range.
- Toca 511/Toca FC (a retroviral cytosine-deaminase gene therapy, not a stem-cell carrier),
  Cloughesy 2020, JAMA Oncol, PMID 33119048, DOI 10.1001/jamaoncol.2020.3161. Randomised
  phase 2/3, recurrent: "The median OS was 11.10 months for the Toca 511/FC group and 12.22
  months for the control group (hazard ratio, 1.06; 95% CI 0.83, 1.35; P = .62)." The only
  randomised enzyme/prodrug trial in glioma was null despite promising phase 1 data
  ("a mOS of 11.9 months in 56 patients", PMID 30676111).
- DNX-2401 (Delta-24-RGD, directly injected, not carrier-delivered), Lang 2018, J Clin
  Oncol, PMID 29432077, DOI 10.1200/JCO.2017.75.8219: "20% of patients survived > 3 years
  from treatment, and three patients had a >= 95% reduction in the enhancing tumor (12%)"
  (recurrent, single arm, n = 25).
- CAR T (not stem-cell carriers): EGFRvIII CAR T plus pembrolizumab, newly diagnosed, n = 7,
  PMID 38216766: "median overall survival (OS) ... 11.8 months (90% CI, 9.2-14.2 months)";
  "no efficacy". Bivalent EGFR/IL13Ra2 CAR T, intraventricular, recurrent, n = 18, PMID
  40451950, DOI 10.1038/s41591-025-03745-0: "The maximum tolerated dose was determined to be
  2.5 x 10^7 cells. Of the 18 patients, 10 (56%) experienced grade 3 neurotoxicity"; "objective
  radiographic response, 8%"; "Median progression-free survival was 1.9 months".

Preclinical data (mouse, used only to bound an assumed effect):
- HB1.F3.CD plus 5-FC, Aboody 2013, Sci Transl Med, PMID 23658244, DOI
  10.1126/scitranslmed.3005365: "The average tumor volume in mice that received HB1.F3.CD
  NSCs and 5-fluorocytosine was about one-third that of the average volume in control
  mice." "no NSCs outside the brain".
- NSC.CRAd-S-pk7, ovarian peritoneal model (not glioma), PMID 42179969: "increased median
  survival by 33% (from 57 to 76 days) in immunocompetent mice, and 48% (from 54 to 80
  days) in immunodeficient mice"; repeated dosing "by 60% (from 57 to 91 days)". Glioma
  multi-cycle IND-enabling work, PMID 34944938: "an improvement in treatment efficacy is
  achieved in immunocompetent mice by administering multiple treatment cycles
  intracerebrally" (no number in the abstract).
- MSC carriers: TGF-beta mediates MSC homing to glioma stem cells, Shinojima 2013, Cancer
  Res, PMID 23365134, DOI 10.1158/0008-5472.CAN-12-3086: "TGF-beta as a tumor factor that
  attracts BM-hMSCs via TGF-beta receptors"; "suggest that TGF-beta predicts patients in
  whom BM-hMSC delivery will be effective". MSCs "show little tropism for the resting and
  differentiated cancer stem cell-like glioma cells", PMID 24481670. In an immunocompetent
  glioma model, MSCs were found at the tumour in "seven of nine mice" after intra-arterial
  injection, PMID 21847363.

### 2.3 NSC migration, the T3 quantities

- Aboody 2000, PNAS, PMID 11070094, DOI 10.1073/pnas.97.23.12846: NSCs "migrate through
  normal tissue targeting the tumor cells", "NSC migration can be extensive". No distance or
  velocity in the abstract.
- Human embryonic NSC chemotaxis in vitro, PMID 23109183, DOI 10.1002/smll.201202208:
  "average[s] about 1 um min(-1)" (= 60 um/h); "CXCL12 concentrations above 9.0 nM are above
  the minimum activation threshold, while concentrations below 14.7 nM are below the
  saturation threshold"; untreated NSCs "are not able to chemotax within CXCL12 gradients".
  In vitro only.
- Tumour expansion, Swanson review 2003, PMID 14607296, DOI 10.1016/j.jns.2003.06.001: "the
  velocity of expansion is linear with time and varies about 10-fold, from about 4 mm/year
  for low-grade gliomas to about 3 mm/month for high-grade ones."
- **Gap.** I found no verified in vivo NSC migration distance, engraftment fraction or
  post-injection half-life. The stub values in `simulation/nsc.py` (45, 28 and 35 um/h; half
  lives 48 to 72 h; tropism scores 0.70 to 0.85) have no source ("illustrative") and must
  not be reused as literature parameters.

**Bottom line on evidence.** There is no published effect size that converts "carrier x
payload x dose x sites x schedule" into survival, and the best randomised evidence for the
nearest mechanism (Toca 5) is null. Every efficacy parameter in any design is assumed or
translated from mouse data. This is the central fact for the limitations (section 6).

---

## 3. Candidate designs

All three keep the registered baseline: Stupp = the frozen survival model's prediction with
arm = radiotherapy plus temozolomide (item 11 (a)). All treat the cell therapy as an
**add-on** to Stupp, as in the NSC-CRAd-S-pk7 trial (cells at surgery, Stupp 10 to 14 days
later), so the counterfactual curve is the Stupp curve modified by the configuration.

### Design A: literature-parameterised proportional-hazards response surface (transparent)

**Structure.** The survival model supplies the Stupp curve; the response surface is a
log-hazard offset, the L4 simulator is a deterministic closed-form delivery model.

1. *Survival model (L5, frozen).* `CoxPHModel` (Breslow) on TCGA, with a step baseline
   cumulative hazard H0(t): S_stupp,i(t) = exp(-H0(t) exp(lp_i)), lp_i = x_i' beta with x_i =
   (age, sex, arm = RT+TMZ) (variant 2) or plus IDH and MGMT (variant 3). No KPS, no
   resection. Evaluated on [0, 24] months.
2. *Delivery module (L4, replaces the stubs).* For patient i and configuration c = (k, p,
   N, n, s):
   - tropism score z_i,k: for NSC carriers the mean of standardised log2 CPM of CXCL12, VEGFA,
     PLAU and CCL2 (the chemoattractants named in the stub's docstring); for MSC, TGFB1. The
     standardisation uses a frozen TCGA reference (mean and SD), never a CGGA cohort
     statistic, to avoid leakage of the sealed cohort into the transform.
   - engrafted dose: N_eff = N * f_k * logistic(a_k + g_k z_i,k), with f_k the fraction of
     cells surviving to act, a_k the baseline tropism logit, g_k the patient slope.
   - spatial coverage: Cov(n) = 1 - exp(-n / n0).
   - schedule factor: sigma(s) = 1, 1 + r2, 1 + r3 (single, x2, x3).
3. *Response.* Per payload p, the log-hazard reduction is
   E_i(c) = Emax_p * m_p(i) * sigma(s) * Cov(n) * N_eff^h / (N_eff^h + ED50_p^h),
   with a patient modifier m_p(i) = exp(w_p' u_i), u_i = (standardised immune-TME score,
   MES-like fraction, subtype indicators). A "combo" combines by Bliss independence on the
   hazard scale, E_combo = E_a + E_b, times a synergy factor (1.0 assumed). **Which two
   payloads it combines was decided per carrier on 2026-11-27 (section 7, Q4), because a
   single uniform pair would place an excluded payload on a carrier:** iPSC-NSC takes CRAd +
   prodrug_enzyme, MSC takes prodrug_enzyme + CAR, HB1_F3_line takes CRAd + prodrug_enzyme.
4. *Toxicity and delay (this is the only source of harm).* Excess log hazard
   T_i(c) = tau_p * (N / 5e7)^q * (1 + eta (n-1)) + delta_delay, with delta_delay the cost of
   delaying Stupp 10 to 14 days (PMID 34214495), all assumed.
5. *Curve and RMST.* S_c,i(t) = S_stupp,i(t)^{exp(-E_i(c) + T_i(c))};
   RMST_i(c) = integral_0^24 S_c,i(t) dt, by the trapezoid rule on a fixed 0.1-month grid;
   gain_i(c) = RMST_i(c) - RMST_i(Stupp); T4's objective is RMST_i(c).

**Parameters.** Sources are verified published numbers or "assumed" with the sensitivity
range (section 5). "Fitted" means fitted on TCGA (development data), then frozen.

| Parameter | Meaning | Value | Source / status | Range |
| --- | --- | --- | --- | --- |
| H0(t), beta | Stupp curve | fitted | TCGA fit of `CoxPHModel`; plausibility check only: cohort-mean RMST24 about 14.5 vs derived 14.55 from PMID 15758009 anchors | n/a |
| tau | horizon | 24 months | registered (item 11) | fixed |
| MGMT, IDH in baseline | prognostic terms (variant 3) | fitted | literature check: MGMT HR 0.45 (0.32 to 0.61), PMID 15758010; IDH HR 0.297 (0.157 to 0.564), PMID 19636000 | fitted CI |
| Emax_CRAd | max log-hazard reduction | skeptical-translation prior, section 4 | assumed; bounded by mouse median-survival ratios 1.33 to 1.48 (PMID 42179969, ovarian) and null human randomised data (HR 1.06, PMID 33119048) | HR 1.0 to 0.60 |
| Emax_prodrug | same | as above | assumed; mouse tumour volume ratio about 1/3 (PMID 23658244) | HR 1.0 to 0.60 |
| Emax_CAR | same | as above | assumed; no efficacy in new GBM (PMID 38216766), 8% objective response (PMID 40451950) | HR 1.0 to 0.70 |
| ED50_p, h | dose-response midpoint and slope | 2.5e7 cells; h = 1 | assumed; the clinically tested ranges are 1e7 to 5e7 (PMID 27979915) and 5e7 to 1.5e8 (PMID 34214495); CAR MTD 2.5e7 (PMID 40451950) | ED50 5e6 to 1e8; h 0.5 to 2 |
| f_k | fraction of cells acting | iPSC-NSC 0.5, HB1.F3 0.5, MSC 0.3 | assumed; no verified in vivo value (gap) | 0.1 to 1.0 |
| a_k, g_k | tropism intercept and patient slope | a = 0; g = 0.5 per SD | assumed; direction supported for MSC by PMID 23365134 (TGF-beta) and for NSC by PMID 23109183 (CXCL12, in vitro) | g 0 to 1 |
| n0 | coverage saturation | 3 sites | assumed | 1.5 to 6 |
| r2, r3 | schedule gain | 0.10, 0.15 | assumed; direction only from PMID 34944938 | 0 to 0.3 |
| w_p | immune, MES slopes | 0.1 per SD, immune for CRAd and CAR only | assumed | 0 to 0.3 |
| tau_p, q, eta | toxicity | CAR 0.05, CRAd 0.02, prodrug 0.01, combo sum; q = 2; eta = 0.05 | assumed; direction from PMID 40451950 (56% grade 3 neurotoxicity, CAR) and PMID 34214495 (1 of 12 grade 3 meningitis) | tau x0.5 to x3 |
| delta_delay | Stupp delay cost | log HR +0.02 | assumed; delay 10 to 14 days from PMID 34214495 | 0 to 0.08 |

**Which features modulate the effect, and how.**
- *Payload mechanism:* a payload-specific Emax and its patient modifier. CRAd and CAR depend
  on the immune and stromal microenvironment (`tme_immune`, assumed sign); the prodrug
  payloads (5-FU, SN-38) act through DNA-damaging antimetabolites and topoisomerase
  inhibition, not through MGMT, so they receive no MGMT modifier (a mechanistic argument,
  not a measured effect).
- *NSC tropism:* z_i,k from the RNA-seq chemoattractant score (NSC) or TGFB1 (MSC); these
  are real CGGA inputs and the direction of the TGFB1 effect is published (PMID 23365134).
- *MGMT and IDH:* act on the **baseline** only (variant 3), and so on the *gain* only
  through the nonlinearity of RMST in the baseline (a fixed HR yields a larger gain when the
  baseline curve is higher; on the derived Weibull, HR 0.6 gives +3.26 months at scale 0.5
  and +2.21 at scale 1.5 of the Stupp curve, so gain is not monotone in baseline). **As
  registered (variant 2; section 7, Q1 and Q6) neither acts anywhere in the frozen model:**
  the baseline varies only with age and sex, and no IDH or MGMT modifier is added to the
  configuration effect, because none is sourced. The sentence above describes what variant 3
  would have done.
- *Subtype and state fractions:* enter only through the tropism score and the modifier m_p.
  Their effect sizes are assumed; the stub's "affinity" table has no source and is dropped.

**Calibration data needed.** None for the structure. Optional T3 data would calibrate f_k,
a_k, g_k, n0; the survival model needs only TCGA. **T3 data are not required.** Without them
the module is "uncalibrated" (item 8) and the Emax values stay assumed regardless.

**Feasibility and time.** Pure NumPy, about 3 to 5 weeks of work for one analyst (module,
tests, TCGA-only dry run, checksums, amendment text). Each evaluation is microseconds, so
T4 (6000 evaluations) and T5 are trivial. Can be frozen in January 2027.
**[Measured 2026-11-28: wrong. The objective is microseconds; the GP surrogate refit around it
is not. One patient at budget 60, both methods, takes 33.3 s. See section 7.3 of
`docs/project_record/07_timeline_2026_11_28.md` for what this does to the sensitivity plan.]**

### Design B: reduced-order mechanistic simulator (coupled PDE, deterministic)

**Structure.** A 1D radial (or 2D axisymmetric) reaction-diffusion system on a generic
spherical resection cavity (no patient MRI exists for CGGA):
- tumour cells c: dc/dt = div(D div c) + rho c (1 - c) - k_p(n, v) c;
- chemoattractant s: ds/dt = Ds lap s + q_i c - lambda s;
- carrier cells n: dn/dt = div(mu div n) - chi div(n div s) - delta n + source(N, n_sites,
  schedule);
- payload field v (virus replication, or enzyme-converted drug): dv/dt = production(n) +
  amplification(c) - clearance;
- tumour burden at the end of treatment mapped to a log-hazard reduction by
  E_i(c) = alpha * ln(B_untreated / B_treated), then the same Cox composition as Design A.

**Parameters.** rho, D from glioma literature (Swanson: "about 3 mm/month for high-grade
ones", PMID 14607296, sets the front speed 2 sqrt(rho D) but not rho and D separately);
mu, chi from NSC chemotaxis (about 60 um/h in vitro, PMID 23109183; verified only in vitro);
delta, Ds, lambda, q_i, production and clearance rates, virus burst size: assumed, no source
verified. The mapping alpha from burden to hazard is assumed and is the weakest link; it
makes Design B a mechanistic-looking wrapper around the same unsupported efficacy step. The
tropism modulators (q_i from CXCL12, VEGFA, PLAU, CCL2; TGFB1 for MSC) are as in Design A.

**Calibration data.** T3's two inputs: PDX growth curves (RMSE <= 15% of final volume) and
in vivo NSC tracking (r >= 0.65). Without them every mechanistic rate is a prior and the
mechanistic detail adds assumptions, not evidence. T3 data are strongly desirable (item 8
already allows NOT EVALUABLE).

**Feasibility and time.** About 3 to 5 months (solver, nondimensionalisation, verification
against analytic Fisher-KPP fronts, parameter priors, the burden-to-hazard map, tests, dry
run). Each evaluation about 10 to 100 ms. Freezable by April 2027 if work starts at once;
no slack if any component slips.

### Design C: agent-based PhysiCell simulation (the protocol's description)

**Structure.** PhysiCell (Ghaffarizadeh 2018) with tumour cells, carrier agents moving by
chemotaxis in SDF-1, VEGF, uPA and MCP-1 fields, payload kinetics per agent, coupled to
survival as in Design B.

**Parameters.** Several dozen, nearly all assumed or requiring calibration; no literature
set has been verified. **Calibration data:** at least 10 PDX growth series and in vivo
tracking (T3), neither in hand.

**Feasibility and time.** T4 alone needs 50 patients x 60 evaluations x 2 methods = 6000
simulations; T5 needs 218 x 60 = 13,080 more; with the sensitivity analysis well over
10^5. At 5 to 30 minutes per three-dimensional run that is 10^3 to 10^4 CPU-hours before
sensitivity. Setup, calibration (ABC), build checksumming and reproducibility add 9 to 12
months. **Not freezable by 2027-06-30 and not recommended.**

### Comparison

| | A: PH response surface | B: reduced-order PDE | C: PhysiCell ABM |
| --- | --- | --- | --- |
| Evidence-supported structure | none needed; honest about it | partially (growth speed, chemotaxis in vitro) | partially |
| Assumed parameters | about 25, all listed | about 30, plus the burden-to-hazard map | dozens |
| Needs T3 | no | strongly desirable | required |
| Cost per evaluation | microseconds | 10 to 100 ms | minutes |
| Time to freeze | 1 to 2 months | 3 to 5 months | 9 to 12 months |
| Frozen before 2027-06-30 | yes, with slack | at risk | no |
| Transparency to a reviewer | high | medium | low |

---

## 4. Recommendation

**Freeze Design A**, with the survival model and prior regime set as follows.

> Design A was adopted on 2026-11-27. Item 1 below recommended variant 3; **the lead chose
> variant 2**, for the reason given in section 7, Q1. Read section 7's table for what was
> decided; this section is the proposal as written.

1. **Survival model.** Variant 3 (age, sex, arm, IDH, MGMT), if the lead accepts it (decision
   Q1). Reasons: it is the transferable model that carries the two molecular determinants of
   the Stupp baseline for which there are published effect sizes (section 2.1); its external
   C-index (0.5915) is higher than variant 2's (0.5385), though it too fails T2; and it keeps
   IDH-mutant tumours from being treated like IDH-wild-type ones. Variant 2 is the fallback
   if the lead prefers to freeze the model that decided T2. Either way, report the PH test
   (section 3.8) and state that the model failed external validation.
2. **L4.** Replace the stubs with the deterministic delivery module of Design A, named in the
   amendment as "L4-lite: closed-form delivery model, uncalibrated, literature-parameterised"
   (decision Q2 asks the lead to confirm that this satisfies "the simulator (L4)" of item 11).
3. **Efficacy prior.** Register a *rule* for Emax before the read, not a value tuned for the
   verdict: Emax_p is the log of the inverse of the mouse median-survival ratio of the
   nearest published preclinical study, shrunk by a translation factor of 0.5 on the log
   scale (assumed). With the ratios 1.33 to 1.48 (PMID 42179969) this gives HR of about 0.82
   to 0.87 for the oncolytic payload; the prodrug payload is anchored the same way from PMID
   23658244; the CAR payload is set at the pessimistic end because human CAR data are
   uniformly poor. State plainly that the shrinkage factor is assumed and that this
   places T5 below the 3-month bar by construction (the bar needs HR of about 0.59,
   section 2.1). The break-even analysis (section 5) reports how large the effect would
   have to be.
4. **Why A over B.** The 3-month bar depends on one number, the clinical effect size of the
   carrier-delivered payload, which no design can supply and which T3 cannot supply
   either (T3 measures growth in PDX and NSC tracking, not human survival). B and C add
   assumed mechanisms and an assumed burden-to-hazard map without adding evidence; they
   cost three to twelve months and cannot be validated if T3 is NOT EVALUABLE. A is the
   smallest model that satisfies item 11, can be audited line by line, is cheap enough to
   run the full sensitivity analysis, and can be frozen with slack (target: amendment by
   2027-02-28, four months before the deadline). B remains available as a labelled
   post-registration exploratory analysis if T3 data arrive.
5. **Controls (extending section 3.8, as a pre-read amendment text).** (i) a *no-modifier
   control*: all patient modifiers set to zero, which must give H2 distinct = 1 (the
   optimiser's universal answer); (ii) a *null-payload control*: Emax = 0, which must give
   gain = 0 for every patient and harm = the toxicity term only; (iii) the synthetic dry run
   of section 3.8; (iv) the TCGA-only run of the whole T4/T5 pipeline on the held-out TCGA
   split before the CGGA read, to catch bugs without touching sealed data.

> **Superseded 2026-11-28.** Every item in this timeline was finished on 2026-11-28, about
> three months early, so its dates no longer mean anything; see
> `docs/project_record/07_timeline_2026_11_28.md`. Two of this section's estimates were wrong
> in opposite directions: the build estimate was about 10x too long, and the claim above that
> "each evaluation is microseconds, so T4 and T5 are trivial" costed the objective instead of
> the GP surrogate around it. Measured, one patient at the registered budget of 60 takes
> **33.3 s**, so a T4 plus T5 run is about two hours, and section 5's 1000-draw Latin
> hypercube would be **103 days** if each draw re-ran the optimiser rather than re-scoring.

**Timeline (proposed).** Decision on this document by 2026-12-15. Module, tests and TCGA
dry run by 2027-01-31. Draft amendment text, checksums and the proposed sensitivity ranges
by 2027-02-15. Amendment entered in section 8 by 2027-02-28. Buffer to 2027-06-30 for the
lead's review, a supervisor if appointed, and any T2 improvement-row work that precedes the
seal (the manifest is sealed after T2's two reads and before any T4/T5 computation).

---

## 5. Sensitivity-analysis plan (exploratory, section 3.6)

> **Split into two tiers on 2026-11-28, before any of it ran:**
> `docs/t4t5_sensitivity_plan.md`. This section says *what* to vary but never whether a draw
> re-scores the configuration already chosen or re-optimises to find a new one. Measured, the
> difference is **225,000x** per patient-setting, so item 2's 1000 draws are either 0.7 minutes
> or 103 days. Items 1, 2, 3 and part of 6 re-score; items 4, 5 and the rest of 6 re-optimise.


All analyses are labelled exploratory in every table and caption, never change a verdict,
and run only on the data of the registered T4/T5 read (section 3.6 scope: "They use the same
data ... An analysis that would need a further read of CGGA is not run"). Before the read,
the same code runs on TCGA features and on synthetic files to tune nothing, only to test.

1. **One-at-a-time tornado.** Each assumed parameter in the table of Design A moves across
   its stated range with the others at their central values. Outputs: T5 success fraction,
   T5 harm fraction, T4 distinct-config count, T4 regret reduction, and the mean gain.
2. **Global Latin-hypercube sample.** 1000 joint draws (seed 42) of the assumed parameters
   from uniform (log-uniform for the dose and toxicity scales) distributions over the ranges
   in the table; report the distribution of each output, Spearman rank sensitivity indices,
   and the fraction of draws with T5 success >= 0.50 and harm <= 0.10.
3. **Break-even effect size.** The common multiplier on all Emax_p at which the T5 success
   fraction first reaches 0.50, reported as the hazard ratio the carrier-delivered therapy
   would need to achieve. This is the most informative output: it converts a model-conditional
   result into a statement about the effect size that would be needed, to be set against
   Toca 5 (HR 1.06, 0.83 to 1.35) and the temozolomide benefit (HR 0.63).
4. **Prior regimes.** Skeptical (Emax as registered), moderate (HR 0.80) and optimistic (HR
   0.60, about the break-even) side by side. The optimistic regime is shown for transparency
   and is never promoted to the primary.
5. **Structural variants.** Hill vs linear dose response; combo as Bliss vs additive vs
   best-of-two; tropism and immune modifiers on vs off; coverage function; survival model
   variant 2 vs 3; Weibull-tail vs Breslow-step baseline beyond the last event time.
6. **Section 3.6 items that apply.** Complete-case vs imputed inputs; subgroups by Verhaak
   subtype, MGMT and IDH; two alternative seeds for T4's initial design and bootstrap.
7. **Reporting.** A single table per output, in the T4 and T5 report, showing the verdict
   under the registered model beside the range under sensitivity, with the break-even as the
   headline sentence.

---

## 6. Honest limitations

1. **T5 is a property of the model, not of patients.** Both T5 metrics are computed from the
   project's own response surface. The registered wording (item 6) is binding: "predicted
   RMST gain", never a benefit. The surface's efficacy parameters are assumed, so T5 says
   what follows from them.
2. **T5 is bounded by T2's failure.** The survival model's external C-index is 0.5385
   (variant 2) or 0.5915 (variant 3), and its 18-month Brier score is 0.2617, so its
   patient-level ordering is barely better than chance. The individual Stupp baselines are
   therefore weakly informative, and the section 6 rule 4 text applies: "A T5 pass built on a
   failed T2 is reported with that dependency stated in the abstract."
3. **The effect size comes from nowhere human.** No randomised or single-arm study gives a
   hazard ratio for NSC- or MSC-delivered therapy added to Stupp (section 2.2). The best
   randomised evidence for the nearest mechanism is null (Toca 5, HR 1.06). The three-month
   bar corresponds to a hazard ratio of about 0.59, comparable to the whole temozolomide
   effect. With any evidence-anchored prior T5 is expected to FAIL; that outcome is
   consistent with the registered fallback, "The pre-registered negative result is the
   publication".
4. **H2 and T4 are largely design-determined.** Whether the top five contain three distinct
   configurations depends on the assumed patient modifiers: with no modifiers the optimiser
   must return one universal answer (the no-modifier control); with the registered ones the
   personalisation signal is an assumption written into the surface. A passing T4
   "personalization" row would show the surface has interactions, not that patients differ
   in real response. Likewise GP-EI beats random search on a smooth, deterministic,
   low-dimensional surface almost by construction; T4's regret row tests the optimiser, not
   the biology. If most patients tie at zero regret the row is INCONCLUSIVE under item 5.
5. **Harm is zero unless the surface can harm.** With a benefit-only model the harm rate is
   0 by construction, which would be a vacuous pass of the safety row. The toxicity and delay
   terms are assumed and small; the harm result must be read with that in mind.
6. **Little patient-level heterogeneity in the baseline.** Under variant 2 the baseline
   varies with age and sex only; the gain differs across patients through the modifiers and
   the nonlinearity of RMST. The success fraction will therefore tend to be near 0 or near
   1 as the common Emax crosses the break-even, not a smooth value.
7. **The dose range does not reach the tested range of the leading trial.** The registered
   maximum 5 x 10^7 is a third of the NSC-CRAd-S-pk7 maximum feasible dose (1.5 x 10^8); the
   surface cannot be extended without an amendment, and extrapolation beyond tested doses
   is not attempted.
8. **L4 is not the registered protocol's simulator and is uncalibrated.** T3 is likely NOT
   EVALUABLE (item 8, section 7). The simulation layer is presented as uncalibrated, and
   no claim depends on its quantitative accuracy. No patient-specific geometry or KPS exists
   for CGGA.
9. **IDH-mutant "glioblastoma".** CGGA grade-4 labels may include IDH-mutant tumours, which
   are not glioblastoma under WHO 2021. Handling must be registered (decision Q6).
10. **Mouse-to-human translation.** The shrinkage factor 0.5 is an assumption, not a
    measured translation rate.

---

## 7. Decision points for the lead

**All nine were decided by the project lead on 2026-11-27 and registered the same day as
pre-registration item 11 (g) to (o), with a dated section 8 entry.** The options below are
left as written, so that what was chosen can be read against what was offered.

| | Decided | Note |
|---|---|---|
| Q1 | **(b)** variant 2 (age, sex, arm) | Against this document's recommendation. Variant 3 needs an IDH and an MGMT call for every patient, which `cgga_external_v1`'s eligibility rule never required, so T5 would score fewer than the 218 its scoring set registers (the equivalent drop on CGGA as a whole is 283 to 241). The 0.05 C-index gain does not buy that, on a model that fails T2 either way. |
| Q2 | **(a)** closed-form delivery model, uncalibrated | "L4-lite". Declared uncalibrated whatever T3 returns. |
| Q3 | **(a)** pre-declared shrinkage rule | Emax_p = 0.5 x ln(ratio_p). Places T5 below its bar by construction, stated in advance. |
| Q4 | **per carrier** (not (a), (b) or (c)) | None of the three options offered is consistent with item 4: (a) puts CRAd on MSC and (b) puts CAR on HB1_F3_line, both registered exclusions, while (c) makes two configurations with the same carrier, payload and site count different therapies, which item 4's distinctness rule and T4's distinct-config count depend on. combo is therefore the two payloads each carrier may carry: iPSC-NSC CRAd + prodrug, MSC prodrug + CAR, HB1_F3_line CRAd + prodrug. Bliss independence, synergy 1.0. |
| Q5 | **(a)** toxicity term in the primary surface | It is the only source of harm; without it T5's harm row could not fail and a PASS there would be vacuous. |
| Q6 | **(a)** IDH-mutant patients kept | But note the interaction with Q1: variant 2 has no IDH term, so "with IDH in the baseline" is not available. As registered, neither IDH nor MGMT appears anywhere in the frozen model — inventing a modifier slope for either would add an assumed parameter with no source. IDH-mutant baselines are therefore modelled as IDH-wild-type, which is disclosed, and section 3.6's IDH subgroup is where it is reported. |
| Q7 | **(a)** 2027-02-28 | Four months of slack before the registered closing date, not to be spent on scope. |
| Q8 | **wait, and loosen the subtype clause** | Both T3 rows decided 2027-01-08, independently; the four-subtype clause became three-of-four with a per-subtype floor of two. Registered 2026-11-27; see `docs/t3_data_search.md`. |
| Q9 | **(a)** `configs/response_surface_v1.yaml` | Checksummed in section 8 at the freeze. |

**Q1. Which survival model is frozen for T4 and T5?**
- (a) Variant 3: age, sex, arm, IDH, MGMT (recommended; external C-index 0.5915, also fails T2).
- (b) Variant 2: age, sex, arm (the model that decided T2's row; 0.5385).
- (c) Variant 2 as baseline plus a published-HR IDH and MGMT offset (HR 0.45; 0.297) added
  outside the fit. Not recommended: it mixes fitted and literature coefficients.

**Q2. What counts as "the simulator (L4)" for item 11?**
- (a) A closed-form delivery model, declared uncalibrated (recommended).
- (b) A reduced-order PDE simulator (Design B), with the freeze deadline at risk.
- (c) PhysiCell (Design C); not freezable by 2027-06-30.
- (d) Close T4 and T5 as NOT EVALUABLE now, by deciding not to build the model.

**Q3. How is the efficacy prior (Emax) set?**
- (a) A pre-declared rule: preclinical ratio shrunk by 0.5 on the log scale, which places T5
  below the bar (recommended).
- (b) A central moderate prior of HR 0.80, an assumed round number.
- (c) A prior set at the T5 break-even (HR about 0.59) so that T5 can pass. Not recommended:
  it would be a choice made with the threshold in view.
- (d) Register a family of regimes with equal standing and no primary. Not recommended:
  section 6's single-verdict rule needs one primary.

**Q4. How is "combo" defined?** The registered space names "combo" without saying what it
combines.
- (a) CRAd plus prodrug_enzyme, Bliss independence, synergy 1.0 (recommended).
- (b) CAR plus prodrug_enzyme.
- (c) Any two payloads, per-carrier, with the best pair chosen by the optimiser.

**Q5. Is the toxicity and delay term included in the primary surface?**
- (a) Yes, small and assumed, so harm can occur (recommended).
- (b) No: benefit-only surface; the harm row is then trivially met and must be reported as vacuous.
- (c) Yes, and a larger toxicity prior drawn from the CAR neurotoxicity rate (56% grade 3).

**Q6. How are IDH-mutant patients handled?**
- (a) Kept, with IDH in the baseline (recommended).
- (b) Excluded from T4 and T5 as non-glioblastoma under WHO 2021. This changes the eligible
  n below 218 and needs an amendment note on item 7.
- (c) Kept, and reported separately by IDH status as an exploratory subgroup.

**Q7. What is the freeze date?**
- (a) 2027-02-28 (recommended; leaves four months).
- (b) 2027-04-30.
- (c) 2027-06-30 (no slack; not recommended).

**Q8. T3.**
- (a) Declare T3 NOT EVALUABLE now (if no calibration data are in prospect) and say the
  simulation layer is uncalibrated (recommended).
- (b) Spend the extraction effort on PDX and NSC-tracking series first, and decide by
  2027-01-31.

**Q9. Where do the machine-readable parameter files live?**
- (a) A single versioned `configs/response_surface_v1.yaml`, checksummed in section 8
  (recommended).
- (b) Constants in the module, checksummed through the source file.

---

## 8. Evidence gaps (what to find before the freeze)

1. A verified in vivo NSC engraftment fraction, migration distance and post-injection
   survival time for iPSC-NSC, MSC and HB1.F3 (none found in abstracts; the Lancet Oncol
   trial's supplement and Aboody's studies are the first places to read).
2. Hegi 2005's full-text unmethylated medians and 2-year survival rates by arm.
3. Any hazard ratio for NSC-, MSC- or iPSC-delivered therapy in glioma; none exists in the
   abstracts retrieved.
4. Data for T3: at least ten PDX growth curves across the four Verhaak subtypes and in vivo
   NSC tracking with distances.
5. A decision on whether HB1.F3 and iPSC-NSC need separate parameters at all (no evidence
   separates them).
