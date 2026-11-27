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
| 1 | **Done, 2026-11-26.** Review and commit the working tree, in the suggested split: code and packaging (`0c20999`); documentation (`10cfd93`); this record (`8f79b33`). | The registration records a commit hash, so the registered code must be fixed first. |
| 2 | **Done, 2026-11-26.** The partition-drift CI job now lists its active manifests in `ACTIVE_MANIFESTS` (`ci.yml`). Add the external manifest there, in both steps, in the same commit that adds it. | The job treats every manifest other than the development one as "superseded". It would mislabel the external partition when `cgga_external_v1.json` is added. |
| 3 | Name the comparator model for T2's improvement row (section 5.4). The other open T2 item, whether to run calibration slope, intercept and decision-curve analysis, was decided on 2026-11-27: not done for this thesis. | The comparator changes what the registration must say. |
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

### 5.1.1 Three things found on 2026-11-26 that bear on T1

All three were found by reading the registered text against the code and the files on
disk. None changes a number that has been measured.

1. **No classifier matches the registered T1 rule.** Section 2.2 item 9(d) registers
   "the RNA-seq-trained classifier", applied to each sample alone. The only subtype code
   in the repository, `profiling/verhaak.py`, describes itself as "signature-based
   unsupervised scoring, NOT the supervised RF + 1D-CNN ensemble the protocol
   specifies", because the project holds no subtype labels to train on. And
   `classify_verhaak` standardises each signature's scores across the samples it is
   given before taking the argmax, so a sample's label depends on the rest of its batch.
   Applied to the 108 labelled patients as one batch, it would break the registered
   property that "nothing in the rule depends on the external data beyond the sample
   being scored". T1 needs one of two decisions before it can run. Both change the
   registered text, so both are for the lead:
   - obtain Verhaak reference labels for the TCGA RNA-seq samples, and predeclare up to
     three supervised candidates and one tuning round, chosen by cross-validation on
     TCGA only; or
   - amend item 9(d) before the deposit to name a per-sample deterministic rule that
     needs no labels.
2. **Most of the T1 scoring set lies outside the registered external partition.** The
   scoring set is the 108 labelled patients of the microarray batch. The registered
   eligibility for `cgga_external_v1` (item 7) requires an RNA-seq profile, and only 28
   of the 108 have one, so about 80 T1 patients would be scored outside any sealed
   partition. The OSF text already says T1's scoring set is "fixed by the registered
   rules instead", so the procedure is consistent as written. But the `cgga_external_v1`
   row of section 4.5 calls that partition the "Primary external test for T1, T2, T4 and
   T5", which is not accurate for T1, and should be corrected before the deposit.
3. **The microarray expression matrix is not on disk.** `data/raw/cgga/_cache/` holds
   the clinical archive of `mRNA_array_301` and nothing else from that batch. The other
   CGGA files came from the portal's open download endpoint, without registration
   (`docs/cgga_ingestion.md` section 1). Whether the microarray expression file is on
   the same endpoint, and its exact name and checksum, has not been checked.

## 5.2 The tests not yet run: what each needs

| Test | Registered rule | Scoring set | What is missing today |
| --- | --- | --- | --- |
| T1 subtype accuracy | Exact binomial lower bound >= 0.88, two-sided | 108 CGGA microarray patients with a reference label; below 30 labelled patients it is NOT EVALUABLE | The microarray expression matrix has not been downloaded, and no classifier matches the registered rule (section 5.1.1). The harmonisation rule is fixed. |
| T1 deconvolution | Mean of four per-state Pearson r >= 0.75, CI lower | 100 simulated pseudo-bulk mixtures from the Neftel 2019 atlas | **Run 2026-11-27, before the deposit: PASS** (0.9521, CI 0.9430 to 0.9589). Specification `docs/t1_deconvolution_spec.md`, committed before the run (`f5b9fe6`). Exploratory until the deposit, which must disclose it |
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
| Infrastructure | CI with 5 jobs; lock file (`requirements-lock.txt`, 2026-11-27); DVC and Snakemake definitions | Container image; initialise DVC and configure a remote; pin the GDC data release; a release DOI | Derived data exist on one machine only. |

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

Each changes something a document currently says. The five rows that came from the stage 3
slice were decided on 2026-11-27 (section 5.4.1); the others are still open.

| Decision | Options as recorded | Where it is recorded |
| --- | --- | --- |
| Status of T2 after two failed rows | Report as failed, with the registered fallback; or amend the registered feature set | `docs/cgga_validation.md` section 6.3 |
| Comparator model for T2's improvement row | Name a 2024 model with re-runnable code or weights, or accept NOT EVALUABLE | `docs/preregistration.md` section 2.2 item 2 |
| Calibration slope, intercept and decision-curve analysis | **Decided 2026-11-27: not done for this thesis.** Not a registered row; a further read would not change the T2 verdict; the offset-versus-dispersion question is named as a limitation | `docs/cgga_validation.md` section 12.5 |
| Whether the two pre-seal reads count as the registered reads | The reads happened without a sealed partition, which the procedure requires | `docs/cgga_validation.md` section 11 |
| Evaluation unit | One row per patient (305) was used, not every row (333) | Same |
| KPS | The committed model cannot be validated on CGGA under any analysis | Same |
| MGMT-STP27 | The published assay's coefficients are unverified; a two-probe score is used and labelled as such | `docs/mgmt_idh_calling.md` section 3.1 |
| Attribution of the T1 harmonisation rule | The decision memo attributes it to Verhaak et al. (2010); the pre-registration records it as a project choice | `docs/osf_final_fields.md`, notes |
| Download of the CGGA microarray matrix | Manual download; name, size and checksum to be sent back | T1 prerequisite |
| Resection extent | **Decided 2026-11-27**: see 5.4.1 D3. Absent at source for all 594 patients, so it is `not_available_in_cohort`, never "biopsy"; no view requires it. Finding a source remains open | `06_evidence_reports.md`, stage 3 |
| KPS missing for 64 of the 283 patients with expression | **Decided 2026-11-27**: see 5.4.1 D4. Complete cases wherever an analysis needs KPS; median-plus-indicator only as a labelled sensitivity analysis declared before it is run; never imputed in the store | Same |
| No MGMT or IDH call for 78 of the 283 | **Decided 2026-11-27**: see 5.4.1 D3 and D4. A missing call is no call, never 0; call status has three states; the `molecular_model` view requires both calls | Same |
| Scoring set for the subtype, state and TME calls | **Decided 2026-11-27 for the subtype**: see 5.4.1 D6 and D6-v2 (single-sample, centroids from the v2 training split). The state and TME scores are still cohort-relative over the 283 (section 5.5) | Same |
| Placeholder columns in the harmonised cohort | **Decided 2026-11-27**: see 5.4.1 D5. L1 owns its columns; L3 raises on a duplicate source and on a null-only shadow | Same |

### 5.4.1 Decided on 2026-11-27 (stage 3b)

These are binding. They were fixed from availability counts and from principle, never by
comparing outcomes, and before the store they govern was rebuilt. Each replaces an open row
of the table above, marked there.

| # | Decision | Reason |
| --- | --- | --- |
| D1 | The feature store holds all 594 cohort patients, one row each, and never filters. Every column carries its role, source layer, source input, the SHA-256 of the file it came from, and a class: `measured`, `derived_from_measurement`, `clinical_record` or `not_available` (`<store>_columns.json`). Every missing value carries a reason code from a fixed vocabulary (`<store>_missing.parquet`): `not_recorded_at_source`, `recorded_unknown_at_source`, `lost_in_harmonisation`, `missing_in_source_table`, `not_available_in_cohort`, `no_expression_sample`, `scoring_failed`, `not_tested`, `tested_not_callable`. Outcomes are labels, not features: overall survival and vital status are kept in a separate label table keyed by patient id (`<store>_labels.parquet`), a label column supplied to the store as a feature raises, and a label column in a feature matrix raises. A missing outcome is an exclusion reason for a model-ready view, never a feature. Treatment arm is a feature, not an outcome | A store that filters hides who was lost and why; the stage 3 store held 155 of 594 and only its report said why. Provenance per column makes a value traceable to the file it came from. Keeping labels apart makes outcome leakage into features a raised error instead of a convention |
| D2 | Model-ready populations are separate named views built from the store. Each declares its required columns and emits its own flow counts, with per-split counts for `tcga_gbm_clinical_v2`. Three are defined (`features/views.py`): `profiling` (an expression sample; 283), `transferable_model` (age, sex, treatment arm, and an observed outcome; 570), `molecular_model` (the transferable requirements plus an IDH call and an MGMT call; 298) | One store, several populations: each analysis states its population instead of inheriting the store's |
| D3 | Absent never becomes a value. Missing resection is `not_available_in_cohort`, never "biopsy". A missing molecular call is no call, never 0. Call status has three states: `called`, `not_tested`, `tested_not_callable` (on TCGA: MGMT 415, 179, 0; IDH 368, 225, 1; the one is a patient sequenced on a recurrent sample only). The continuous MGMT score and the threshold used (0.35) are kept | The stage 3 store read every patient as "biopsy", because resection is absent at source for all 594 (Part A item 6), and the slice's `zero` option would have read an untested patient as unmethylated or IDH wild type |
| D4 | KPS and the molecular calls are used as complete cases wherever an analysis needs them. Median-plus-indicator is allowed only as a labelled sensitivity analysis, declared before it is run. KPS is never imputed in the store | Every missing KPS is absent at source (`06_evidence_reports.md`, stage 3b, Part A item 4), so imputation fills a value the record never had; the indicator is itself prognostic (HR 1.57) |
| D5 | Column ownership. The harmonised clinical table (`gliostrata harmonize-clinical`, `data/interim/tcga_gbm_clinical_harmonized.parquet`) carries only the columns the clinical layer owns; a column under an L2 or molecular name that holds a value stops it. `build_feature_store` takes each input by explicit name and raises on a column supplied by more than one input and on a null-only column that shadows a computed feature. The script's workaround is removed, so there is one code path | Stage 3 finding 3: the cohort's null placeholders would have silently replaced the L2 outputs. One code path means the rule cannot be bypassed by calling the store without the script |
| D6 | Subtype: per-patient output always uses single-sample assignment (nearest centroid, Pearson, over the published Verhaak 840-gene signature; each sample's log2(TPM+1) z-scored across its own genes, the registered T1 rule 2.2 item 9(d)). The cohort-relative signature scoring is labelled as such and used only for descriptive cohort distributions. Limits: these two candidates plus one variant (three classes, without Neural), one tuning round (none was used), TCGA labels only | Fixed on principle, not on which mode scored better: a patient's feature must not change when other patients enter or leave the scoring set. The cohort-relative mode changed 2 to 18 calls of 283 under random halving and 5 to 7 of the 45 v2-test store patients when they were scored alone (the leak of review item M6); the single-sample mode changed none |
| D6-v2 | Scoring-set decision for the v2 split. Two frozen centroid sets. `verhaak840_tcga_v2train`, estimated on the v2 training split only (170 labelled patients), is the one used for the store's per-patient subtype and for any internal evaluation on the v2 test split. `verhaak840_tcga_all` (246 labelled patients) is reserved for T1. No v2 test patient contributes to the centroids the store uses | The v2 test split is reserved for internal held-out evaluation (pre-registration 4.5); centroids estimated on it would leak it into every feature derived from them |

## 5.5 Known gaps, found and not fixed

Found in the audit or while writing this record. Listed so that nobody has to rediscover
them.

- **CI manifest labelling. Fixed on 2026-11-26.** The partition-drift job used to
  hard-code which manifest was current and label every other one "superseded". It now
  takes its active manifests from `ACTIVE_MANIFESTS` in `ci.yml`.
- **No machine-readable registration state, no single access ledger. Fixed on
  2026-11-26.** `configs/registration.yaml` (OSF DOI `null` until the deposit) and
  `data/access_ledger.jsonl` (every access to CGGA, back-filled with its source).
  `evaluation/read_guard.py` refuses a value read until the registration is deposited,
  the partition is sealed or a registered rule names a frozen scoring-set manifest
  whose checksum still matches, every input's checksum is the one on disk, and the
  analysis has not been read. Since 2026-11-27 the rule must be an id from
  `configs/scoring_rules.yaml`, not free text, and a rule whose manifest is not yet
  frozen refuses every read. It is a discipline, not a sandbox: a script has to call
  it. The two T2 scripts predate it and do not.
- **Eligibility table not committed.** `data/processed/**` is git-ignored.
- **Result files not committed.** `data/results/**` is git-ignored.
  Their SHA-256 checksums as of 2026-11-26 are recorded in `06_evidence_reports.md`.
- **Unused declared dependencies.** `python-dotenv`, `tqdm` and `joblib` are declared but
  not imported.
- **A deferred consolidation.** `prediction/cox.py` carries its own bootstrap routine with
  a note to consolidate it; changing it would alter random streams.
- **Advisory lint backlog.** 170 findings at `baseline-r1` (171 the day before), mostly
  naming rules. Stage R1 added none. After stage 3b: 166 (four removed in files the stage rewrote; none added).
- **No container image.** The lock file exists since 2026-11-27: `requirements-lock.txt`
  pins the 37 installed distributions of the runtime and `dev` closure (Python 3.14.3),
  and CI's `test` job installs from it. There is still no image, and the `workflow`
  extra (`snakemake`, `dvc`) is not locked, because it is not in the environment that
  produced the results.
- **Deferred from stage R1: the extracted plain-text copies of the CGGA files are still
  in the workspace.** `data/raw/cgga/_cache/extracted/` holds five plain-text files next
  to the five `.zip` archives they came from: two expression matrices (47.7 MB and
  103.6 MB) and three clinical tables. They are what a workspace-wide text search
  opened on 2026-11-27T05:48:55Z (ledger entries 8 and 9; `06_evidence_reports.md`,
  stage R1). The intended state is: keep the compressed source, extract into a
  temporary directory only inside a guarded read, and delete it afterwards. It was not
  done in R1 for two reasons. It deletes files under `data/raw/cgga`, and R1 was allowed
  to touch that tree by name, size and time only. And it changes
  `scripts/download_cgga.py` (`extract_text_member` keeps the extraction as an audit
  trail), a change that cannot be validated without a CGGA read. Interim protection:
  `.clineignore` and `.rgignore` exclude both CGGA directories from workspace-wide
  tools, and `.clinerules` forbids workspace-wide searches in this repository.
- **Deferred from stage R1: not verified that the agent's search tool honours the
  ignore files.** The installed tool's file walk was inspected in its bundle and could
  not be determined; a live test was not run because a failing test would have opened
  the held-out files again. The rule in `.clinerules` is the control that does not
  depend on it.
- **Deferred from stage R1: the new CI jobs have not run on GitHub.** `test` now
  installs `requirements-lock.txt`; `test-floor` is new (Python 3.10,
  `requirements-floor.txt`). Neither was executed, because R1 used no network. The
  floor job tests a claim that had never been tested and may be red.
- **Deferred from stage R1: the second partition-drift step** (the shell loop that
  calls `gliostrata freeze-partition --verify`) is still inline in `ci.yml` and
  untested. The first step's logic is now `scripts/ci_partition_drift.py`, and
  `tests/test_ci_partition_drift.py` runs it.
- **DVC is not initialised.** `.dvc/` does not exist and no remote is configured.
- **The GDC is not pinned to a data release.** A re-run queries today's GDC.
- **Nothing in L2 is validated against a published reference.** The single-sample
  subtype agrees with the published TCGA labels for 56 of 76 v2-test patients (0.737), a
  descriptive check on its own cohort (`06_evidence_reports.md`, stage 3b). T1's
  subtype row has not been run. The state scores recover simulated Neftel mixtures (T1's
  deconvolution row, run before deposit), which is not a check against a published
  reference on tumours.
- **No end-to-end run on the real cohort has been recorded.** L0 to L3 run on TCGA
  through `gliostrata build-feature-store --mode real` since 2026-11-27 (stage 3b; all
  594 patients in the store; `scripts/build_tcga_feature_store.py` is now a thin wrapper
  over it). L4 to L7 do not run on real data: `optimize` and `generate-protocol` refuse
  the real store rather than encode an absent value as a value. `dvc repro` on real data
  still stops at the stub `profile` stage.

Found in stage 3b (2026-11-27), not fixed:

- **Three cohort patients have a KPS of 0 at source.** The raw GDC record says "0" for
  each (all three in the v2 training split). A KPS of 0 means dead, which is implausible
  at diagnosis. The store keeps the value as recorded, and it entered variant 1's TCGA
  fit. Whether it is a coding error at source was not checked.
- **The state and TME scores are cohort-relative.** `deconvolve.py`'s default
  normalisation standardises each signature across the samples scored, so a patient's
  value depends on the others: scoring the first 50 of the 283 alone moves a state
  fraction by up to 0.0744 and changes the highest state of 6 of the 50. The store's
  values were scored over all 283 profiled patients, the 82 v2-test patients among them.
  This is a scoring-set dependence of the kind review item M6 found for the subtype, and
  decision D6 does not cover it. The `proportional` option is cohort-independent but is
  not the default. A decision is needed before these values feed any internal evaluation
  on the v2 test split.
- **The v2 test split has never been held out.** Every reported internal number had the
  169 v2-test patients in its training folds: the three T2 fits (including the TCGA fits
  scored on CGGA), their internal C-indices, the IDH sensitivity runs and the fold-wise
  rerun. No script reads `train_ids` or `test_ids`. The pre-registration discloses this
  for T2 only (`06_evidence_reports.md`, stage 3b, Part A item 8).
- **Snakemake has never executed on the development machine**, so neither rule of stage 3b
  has run there. Snakemake is not installed; the `harmonize_clinical` and `build_feature_store` rules were
  checked only through the tests and by running the same `gliostrata` commands by hand. CI's
  `workflow-dry-run` job will be their first run.
- **`cgga_external_v1` is still blocked** (Part A item 1). The real eligibility step has
  not been run; `scripts/cgga_external_eligibility.py` does not call
  `read_guard.record_access`, so running it as it stands would be an unledgered access;
  and its output under `data/processed/**` is git-ignored, so it must be checksummed in
  the commit that adds the manifest. The guard itself would admit the step: it is an
  `availability` access, and only a `values` access is subject to the deposit and seal
  gates. A script must still call it.
- **Corrected on 2026-11-27, after stage 3b: `configs/scoring_rules.yaml` said
  "per-gene z-scored vector"** for the T1 rule, while pre-registration 2.2 item 9(d) says
  each sample is z-scored across its genes. The registry now uses 9(d)'s wording. Its
  rule id, analyses and (null) manifest fields are unchanged, so the read guard's
  behaviour is unchanged.
- **`treatment_arm` takes treatments from the first diagnosis record only**
  (`download_tcga_clinical.py`, `flatten_case`). The live GDC lists 136 patients with a
  "Surgery, NOS" record against 39 in the raw file. Whether this changes `treatment_arm`
  was not checked.
- **Corrected on 2026-11-27, after stage 3b: pre-registration sections 3.3 and 4.5 named
  `cgga_external_v1` as T1's test partition.** It is not: T1's subtype row is scored on
  the 108 labelled CGGA microarray patients (2.2 item 9; only 28 of them have RNA-seq,
  which the partition requires), and its deconvolution row on simulated mixtures (item
  10). Found in stage 3b (Part A item 2), which could not edit the pre-registration;
  applied afterwards at the project lead's request, with the proposed text, and logged
  as a dated amendment (section 8). No threshold, bound, rule, scoring set or partition
  changed. The OSF form text predates it, like the other corrections (section 5.1,
  step 6).

  And in section 3.3: "Final numbers for T2 are reported on the frozen external test
  partitions; T1 is scored on the sets fixed in section 2.2 items 9 and 10."
- **The stale stage 3 store is still on disk.** `data/processed/tcga_gbm/feature_store.parquet`
  (155 patients, git-ignored) sits beside the new `data/processed/tcga_gbm_feature_store*`
  files; nothing reads it.
- **The 0.737 subtype agreement is mostly with another classifier.** 65 of the 76 v2-test
  labels are Brennan 2013 ssGSEA classifier labels; 11 are Verhaak 2010 clustering labels.
- **Resolved 2026-11-27: the Neftel 2019 atlas is on disk** (GEO GSM3828672 and the 3CA
  per-cell annotation; `data/provenance/neftel2019_sources.json`) and T1's deconvolution
  row has been run once, before the deposit (`06_evidence_reports.md`, section "T1,
  state-deconvolution row").
- **Corrected 2026-11-27: the pre-registration and the OSF text now state the atlas and
  the result.** Section 5 records the Neftel 2019 atlas as obtained and T1's deconvolution
  row as evaluated before the deposit; a new section 2.4 records the result in the form of
  section 2.3; both changes are dated amendments marked as written after the run.
  `docs/osf_final_fields.md` gives the atlas its own paragraph, the result under "Sample
  size", and discloses the run as pilot-phase item (c); `docs/t1_deconvolution_spec.md`
  joins the attachments.
- **The atlas inputs carry no stated licence.** NCBI places no restriction on use or
  distribution of GEO data (submitters may claim rights); 3CA states no terms and
  publishes datasets with the original authors' permission. Recommended, not yet done:
  keep both files out of any release (download and verify by checksum instead), cite
  Neftel et al. 2019 and Gavish et al. 2023, and ask 3CA (3CA@weizmann.ac.il) to confirm
  reuse and whether the Neftel malignant flags and state scores are the original
  authors' annotations.
- **`scripts/prepare_subtype_reference.py` needs `openpyxl` and `pyreadr`**, which are in
  neither `pyproject.toml` nor the lock file, to rebuild the published-label table and the gene map from the
  downloads. The committed centroid sets and gene map do not need it.

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
