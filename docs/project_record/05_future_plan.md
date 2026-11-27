# 5. Future plan

As of 2026-11-26. This document separates three kinds of item, because they carry
different weight:

- **Registered commitment.** Written in `docs/preregistration.md`. Changing it needs a
  dated, disclosed amendment.
- **Required before deposit.** Steps that must happen for the registration to be true and
  complete.
- **Proposed.** Engineering direction from `docs/roadmap.md` and the audit. Not registered;
  the project lead decides.

## 5.1 Before anything else: commit, seal, deposit

The order below is a recommendation; the reasons are given so it can be judged.

| Step | Action | Why it comes here |
| --- | --- | --- |
| 1 | Review and commit the working tree. Suggested split: code and packaging; documentation; this record. | The registration records a commit hash, so the registered code must be fixed first. |
| 2 | Fix the partition-drift CI job before any new manifest is committed. | The job treats every manifest other than the development one as "superseded". It would mislabel the external partition when `cgga_external_v1.json` is added. |
| 3 | Decide the two open T2 items: the comparator model, and whether calibration slope, intercept and decision-curve analysis are run. See section 5.4. | Both change what the registration must say. |
| 4 | Write the eligibility table, then seal the external partition (commands below). | T4 and T5 cannot start without it. The pre-registration places the seal after the T2 evaluation and before any T4 or T5 computation. |
| 5 | Commit the manifest and update the `cgga_external_v1` row in pre-registration section 4.5. | The registration must name the sealed partition. |
| 6 | Re-paste all four OSF fields from `docs/osf_final_fields.md`; upload the attachments; fill in the DOI, date, authors and commit hash in the pre-registration header; deposit. | The text on the OSF form predates the latest corrections. |
| 7 | Record the deposit in this folder and in `CHANGELOG.md`. | So the record stays true. |

**Step 4 in full.** The eligibility script fits nothing and scores nothing. It reads
expression availability from column names only and writes identifiers and a batch label.

```bash
python scripts/cgga_external_eligibility.py --dry-run   # counts only; already run: 218 eligible
python scripts/cgga_external_eligibility.py             # writes data/processed/cgga/cgga_external_eligible.parquet
gliostrata freeze-partition --input data/processed/cgga/cgga_external_eligible.parquet \
    --manifest data/partitions/cgga_external_v1.json --name cgga_external_v1 --all-test
gliostrata freeze-partition --input data/processed/cgga/cgga_external_eligible.parquet \
    --manifest data/partitions/cgga_external_v1.json --verify
```

**A gap to close alongside step 5.** `data/processed/**` is git-ignored, so the eligibility
table will not be committed; only the manifest will. Either un-ignore that one file or
commit its checksum.

**A gap to close alongside step 6.** The result files (`data/results/cgga/*.json`,
`data/results/tcga_gbm/idh_sensitivity_*.json`) are git-ignored, yet the pre-registration
cites their paths. Upload them to the OSF as results, or commit their checksums.

## 5.2 The tests not yet run: what each needs

| Test | Registered rule | Scoring set | What is missing today |
| --- | --- | --- | --- |
| T1 subtype accuracy | Exact binomial lower bound >= 0.88, two-sided | 108 CGGA microarray patients with a reference label; below 30 labelled patients it is NOT EVALUABLE | The microarray expression matrix has not been downloaded. The harmonisation rule is fixed. |
| T1 deconvolution | Mean of four per-state Pearson r >= 0.75, CI lower | 100 simulated pseudo-bulk mixtures from the Neftel 2019 atlas | The mixtures have not been simulated. |
| T2 improvement over baseline | Paired difference in C-index >= 0.03, CI lower | Identical frozen partition for both models | **No comparator model is named.** If none can be re-run, the row is NOT EVALUABLE. |
| T3 tumour growth | RMSE <= 15% of final volume, CI upper | At least 10 PDX control-arm curves spanning four subtypes | **No qualifying dataset.** Data due 2027-01-15; calibration by 2027-03-01; otherwise NOT EVALUABLE. |
| T3 NSC migration | r >= 0.65, CI lower | In vivo tracking series | **No qualifying dataset.** Same closing rule. |
| T4 regret reduction | >= 25% and two-sided Wilcoxon p < 0.01 | First 50 patients by sorted id in the sealed external partition; 60 evaluations each | The partition is not sealed; the simulator and the response surface are stubs. |
| T4 distinct configurations | >= 3 in the top five across 50 patients | Same | Same |
| T5 success and safety | >= 0.50 with CI lower; harm <= 0.10 and CI upper <= 0.15 | All 218 eligible held-out patients | Same, and its credibility is bounded by T2, which failed |

**The dependency that matters most.** T4 and T5 score the optimiser against the simulation
and survival layers. Today L4 (simulation) is a stub, and the survival model the optimiser
calls is a stub. Running T4 or T5 now would measure synthetic outputs. Those tests become
meaningful only when a real simulator and a real response surface exist. The
pre-registration already states that T5 is a statement about **predicted** outcomes, and
that no T5 result may be described as a survival benefit.

**Multiplicity.** Benjamini-Hochberg at FDR 0.05 is applied across the family where a
p-value exists. CI-based thresholds are not re-corrected, but the number of thresholds
evaluated is always reported.

**Each of T1 and T4/T5 is one further registered read of CGGA.** Two reads have been made,
both for T2.

**Wet-lab tests T6 and T7** (protocol completeness on the top ten configurations; a
three-arm organoid experiment on at least ten patient-derived organoids) are outside the
pre-registration and have not been started.

## 5.3 Implementation strategy by layer (proposed)

The order follows the recommended implementation order in `docs/roadmap.md`. Nothing in
this section is registered.

| Layer | What exists | What is needed | Note |
| --- | --- | --- | --- |
| L0 ingestion | TCGA clinical, RNA-seq, IDH, MGMT; CGGA RNA-seq | CGGA microarray batch for T1; CPTAC not started; GLASS not downloaded | The one layer that does what it says end to end. |
| L1 harmonisation | Schema normalisation and QC rows | Batch correction for count data; gene-identifier harmonisation; an imputation policy for molecular covariates | The microarray batch is intensity data, so a count-based method does not apply to it. |
| L2 profiling | Expression-based subtyping and deconvolution, tested, not called by the CLI | Wire both into `gliostrata profile`; the pipeline rule and `dvc.yaml` must gain an expression input; validate against published calls | The CIBERSORTx signature-matrix builder is not implemented. |
| L3 feature store | Working | Consume real L2 output | It currently consumes the stub. |
| L5 outcome model | Fitted Cox model; external validation failed | A decision, not more code: see section 5.4 | The stub is what the optimiser calls. |
| L4 simulation | Synthetic curves | PhysiCell wrapper; NSC chemotaxis and payload kinetics; calibration against PDX and tracking data | Not in the pipeline at all. |
| L6 optimiser | Ranked enumeration (in the CLI); GP / expected-improvement (not wired) | Wire the GP optimiser; BoTorch acquisition | T4's registered budget assumes the GP optimiser: 12 initial design points and 48 iterations. |
| L7 handoff | JSON protocol documents | Convert to a word-processor document; make it bench-executable | |
| Infrastructure | CI with 4 jobs; DVC and Snakemake definitions | Lockfile and container image; initialise DVC and configure a remote; pin the GDC data release; a release DOI | Derived data exist on one machine only. |

**A strategic point for L5.** Three routes were listed for the project lead when T2
failed; none has been chosen:
- keep the registration as written and report T2 as failed, re-framing the contribution
  as benchmark establishment (the pre-registered fallback);
- validate the committed four-feature model on a cohort that records performance status
  (GLASS or CPTAC-GBM; neither has been checked for this);
- build a model whose primary features are molecular or expression-based. That is a new
  model and needs its own pre-registered external test on a partition that has not been
  read; CGGA has been.

Changing the 0.65 threshold was listed only to reject it.

## 5.4 Open decisions for the project lead

None of these has been made. Each changes something a document currently says.

| Decision | Options as recorded | Where it is recorded |
| --- | --- | --- |
| Status of T2 after two failed rows | Report as failed, with the registered fallback; or amend the registered feature set | `docs/cgga_validation.md` section 6.3 |
| Comparator model for T2's improvement row | Name a 2024 model with re-runnable code or weights, or accept NOT EVALUABLE | `docs/preregistration.md` section 2.2 item 2 |
| Calibration slope, intercept and decision-curve analysis | Run them as a further, disclosed read with the rules fixed first; or report them as not done | Pre-registration section 7 |
| Whether the two pre-seal reads count as the registered reads | The reads happened without a sealed partition, which the procedure requires | `docs/cgga_validation.md` section 11 |
| Evaluation unit | One row per patient (305) was used, not every row (333) | Same |
| KPS | The committed model cannot be validated on CGGA under any analysis | Same |
| MGMT-STP27 | The published assay's coefficients are unverified; a two-probe score is used and labelled as such | `docs/mgmt_idh_calling.md` section 3.1 |
| Attribution of the T1 harmonisation rule | The decision memo attributes it to Verhaak et al. (2010); the pre-registration records it as a project choice | `docs/osf_final_fields.md`, notes |
| Download of the CGGA microarray matrix | Manual download; name, size and checksum to be sent back | T1 prerequisite |

## 5.5 Known gaps, found and not fixed

Found in the audit or while writing this record. Listed so that nobody has to rediscover
them.

- **CI manifest labelling.** The partition-drift job hard-codes which manifest is current
  and calls every other one "superseded". Fix before committing the external manifest.
- **Eligibility table not committed.** `data/processed/**` is git-ignored.
- **Result files not committed.** `data/results/**` is git-ignored.
- **Unused declared dependencies.** `python-dotenv`, `tqdm` and `joblib` are declared but
  not imported.
- **A deferred consolidation.** `prediction/cox.py` carries its own bootstrap routine with
  a note to consolidate it; changing it would alter random streams.
- **Advisory lint backlog.** 171 findings, mostly naming rules.
- **No lockfile or container.** Two clean installs a month apart can resolve differently.
- **DVC is not initialised.** `.dvc/` does not exist and no remote is configured.
- **The GDC is not pinned to a data release.** A re-run queries today's GDC.
- **Nothing in L2 is validated against a published reference.**
- **No end-to-end run on the real cohort has been recorded.**

## 5.6 Milestones and publication plan, as stated in the protocol

The protocol's eight-month map, unchanged since it was written:

| Month | Milestone | Tests |
| --- | --- | --- |
| M1 | L0 to L3 complete | none |
| M2 | L4 simulation and L5 trained | none |
| M3 | L6 optimiser and L7 scaffolding | none |
| M4 | End-to-end integration | none |
| M5 | Frozen-test-set validation | T1 to T5 |
| M6 | Analysis, figures, sensitivity analyses | T1 to T5 |
| M7 | Resource release and paper draft | T1 to T5 |
| M8 | Submit to Nucleic Acids Research (Database issue) | T1 to T5 |

Fallback venue: *Bioinformatics*. An applied follow-up paper depends on the wet-lab tests.
Source: `docs/protocol.md` Chapter 9.

**Two things to notice about this map.**

- It places validation at M5, after simulation and integration. In fact **T2 was evaluated
  first**, ahead of that order and before deposit.
- The map was written before T2 failed. The registered response to that failure is to
  report the threshold as FAIL and re-frame the contribution as benchmark establishment.
  The map has not been re-baselined against that outcome.
