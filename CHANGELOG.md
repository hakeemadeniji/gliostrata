# Changelog

All notable changes to GlioStrata will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — T4/T5 model design registered (2026-11-27, after the OSF submission)
- Pre-registration item 11 gains (g) to (o): the nine choices `docs/t4t5_model_design.md`
  left to the project lead are decided and registered. Design A is adopted. Frozen survival
  model: variant 2 (age, sex, arm), so that T5 scores all 218 patients. Simulator: a
  closed-form delivery module, declared uncalibrated. Efficacy prior: the rule
  Emax_p = 0.5 x ln(ratio_p), which places T5 below its bar by construction, stated in
  advance. `combo` is defined per carrier. The toxicity and delay term is in the primary
  surface. All 218 patients, with IDH and MGMT in no part of the model, and the cost
  disclosed. Freeze by 2027-02-28. Parameters in `configs/response_surface_v1.yaml`. Four
  pre-read controls added to section 3.8.
- No threshold, decisive bound, decision rule, partition or scoring set changed; item 4's
  configuration space, its two exclusions, its 10 pairs and its 240 combinations stand.
- Three decisions went against the design document's recommendation, each with its reason in
  the text: variant 2 over variant 3, a per-carrier `combo` (all three options offered
  contradicted item 4), and IDH and MGMT in no part of the model rather than as invented
  modifiers.

### Changed — T3 data rule (2026-11-27, after the OSF submission)
- Pre-registration item 3(b): a qualifying PDX source must now span at least three of the
  four Verhaak subtypes, with at least two curves per represented subtype, instead of all
  four. Modern panels label with the Wang 2017 scheme, which has no Neural class, so no
  obtainable dataset could meet the clause as written (`docs/t3_data_search.md`). Any absent
  subtype is named in section 8 when the data arrive, and the coverage limit is reported with
  T3's verdict.
- Pre-registration items 3(e) and 8: both T3 rows are decided on 2027-01-08, independently.
  A row without a secured qualifying source by that date is declared NOT EVALUABLE then, so
  that the T4/T5 freeze is made knowing whether the simulation layer is calibrated.
- Neither T3 threshold changed, nor the RMSE normalisation, nor the bootstrap unit. The
  amendment was made before any T3 source was contacted and before any T3 data was obtained.

### Deposit readiness (2026-11-27)
- Pre-registration header complete except the OSF DOI and the deposit date; the registered
  commit is the tag `prereg-deposit`. Repository: https://github.com/hakeemadeniji/gliostrata
  (private until the deposit). Author email rewritten before the first push
  (`docs/project_record/commit_map.md`).
- CI green on GitHub, the click 8.1 floor included. The sealed manifests are stored
  byte-exact. T2's open items are decided in pre-registration 2.3; GLASS is not used.
  `scripts/osf_bundle.py` assembles the OSF attachments.

### Added — T2 comparator (2026-11-27)
- The registered comparator search (`docs/t2_comparator_search/`) was run once (549 records,
  dual screening, full-text review). It found one model eligible as first registered:
  PMID 38174685, a 7-gene signature. It decides T2's improvement row, applied as published on
  CGGA's own FPKM release (`configs/t2_comparator/primary.json`,
  `data/provenance/cgga_fpkm_sources.json`, `scripts/cgga_fpkm.py`).
- Option B, a refitted best published gene set (`gliostrata.prediction.geneset_comparator`,
  `scripts/t2_comparator.py`), is kept as a descriptive secondary. The improvement row uses a
  batch-stratified C-index. Pre-registration section 2.2 item 2 carries four dated amendments,
  and OSF stopping rule 9 gains one read, after the deposit.

### Added — deposit preparation (2026-11-27)
- `docs/t2_comparator_shortlist.md`: a proposed eligibility rule for T2's 2024 comparator,
  four candidates checked against it (none eligible), and three paths for the lead. Nothing
  registered; no CGGA file opened.
- `docs/deposit_checklist.md`: the ordered steps to an OSF deposit, with two blocking text
  fixes found in pre-registration section 5 and the T1 classifier that section 2.2 item 9(d)
  does not yet name.

### Added — the sealed external partition (2026-11-27)
- `data/partitions/cgga_external_v1.json`: all-test, 218 CGGA patients eligible by the
  registered rule, sealed after the T2 evaluation and before any T4 or T5 computation. It
  records the SHA-256 of the git-ignored eligibility table. CI's partition-drift job lists
  it in `ACTIVE_MANIFESTS`.
- `scripts/cgga_external_eligibility.py` records every run in the access ledger through
  the read guard before opening a CGGA file, reads only the four columns the rule needs,
  and has `--seal`. The seal was one availability access (ledger entry 11).
- Pre-registration section 4.5 names the sealed partition (one dated amendment); section
  5, `docs/osf_final_fields.md` and `docs/cgga_validation.md` updated to match.

### Added — T1 state-deconvolution row (2026-11-27; run before deposit)
- `docs/t1_deconvolution_spec.md`: the mixture-generation specification, committed
  before the run with `scripts/t1_deconvolution_mixtures.py`, its 13 tests and
  `data/provenance/neftel2019_sources.json` (GEO GSM3828672 and the 3CA per-cell
  annotation of the Neftel 2019 atlas; both git-ignored under `data/external/`).
- Result, run once: mean per-state Pearson r 0.9521 (95% CI 0.9430 to 0.9589) on 100
  pseudo-bulk mixtures of malignant cells, against 0.75: PASS by the registry. Run
  before the deposit, so exploratory until it; optimistic by design (malignant cells
  only; the atlas states and the estimator's markers share the Neftel meta-modules).
  Declared secondary (`fraction_method="proportional"`): 0.8541 (0.8358 to 0.8750).
- Pre-registration: a new section 2.4 records the result (in the form of section 2.3), and
  section 5 records the atlas as obtained and the row as evaluated before the deposit;
  two dated amendments. `docs/osf_final_fields.md` updated to match, with the
  specification added to the attachments. The provenance record gains the verified
  terms of both sources.

### Added — pre-registration and the real data layers
- Pre-registration machinery: checksum-frozen threshold registry
  (`evaluation/thresholds.py`), frozen partitions with drift detection
  (`evaluation/partitions.py`), and the `freeze-partition` CLI command,
  including the `--all-test` seal for external cohorts.
- `docs/preregistration.md` (OSF pre-registration, not yet deposited) and
  `docs/osf_final_fields.md` (paste-ready OSF text).
- L0 ingestion: TCGA-GBM bulk RNA-seq, IDH mutation calls and MGMT methylation
  calls (`scripts/download_tcga_*.py`); CGGA RNA-seq ingestion
  (`scripts/download_cgga.py`) and the external-partition eligibility listing
  (`scripts/cgga_external_eligibility.py`).
- L2 from bulk expression, not wired into the pipeline: signature-based Verhaak
  subtype scoring (`profiling/verhaak.py`) and marker-based state/TME scoring
  (`profiling/deconvolve.py`). `scripts/join_molecular_calls.py` and
  `scripts/export_cibersortx.py`.
- L5: fitted Cox proportional-hazards model (`prediction/cox.py`), its T2
  external validation on CGGA (`scripts/validate_cox_cgga.py`) and the IDH
  sensitivity analyses (`scripts/idh_sensitivity_*.py`).
- T2's Brier row (score at 18 months): `scripts/validate_brier_cgga.py`, evaluated once on
  a second read of CGGA and not met (0.2617, 95% CI 0.2391 to 0.2829, against 0.25).
- `docs/project_record/`: a full project record (overview, findings, work log, commit
  history and future plan). Every number in it was cross-checked against its source.
- Read discipline for the held-out cohort: an append-only access ledger
  (`data/access_ledger.jsonl`, with every access to CGGA to date back-filled and
  sourced), a read guard (`evaluation/read_guard.py`) that refuses a value read
  until the registration is deposited and, where one applies, the partition is
  sealed and the analysis not already read, and `gliostrata access-ledger` to
  print the ledger.
- `configs/registration.yaml`: the machine-readable registration state. The OSF
  DOI is `null` until the deposit; a test keeps this file and the
  pre-registration header in agreement.
- L3 from real data (superseded in stage 3b, below: the script is now a thin wrapper
  over `gliostrata build-feature-store --mode real`, and the store keeps all 594
  patients): `scripts/build_tcga_feature_store.py` runs L0 to L3 on the
  TCGA-GBM inputs on disk (155 of 594 patients; every exclusion is counted and
  reported, with the SHA-256 of each input). TCGA only: it reads no CGGA file and
  fits no model. It found that resection extent is absent for the whole cohort and
  that the harmonised cohort's placeholder columns would shadow the L2 outputs in
  `features/store.py`.

- L6: Gaussian-process / expected-improvement optimizer
  (`optimization/bayes_opt.py`).
- Snakemake workflow, DVC pipeline definition and GitHub Actions CI.

### Changed
- The TCGA-GBM clinical cohort is 594 patients (was 583): only contradictory
  diagnoses are dropped. Partition `tcga_gbm_clinical_v2` is the current sealed
  split; `v1` is kept unmodified as a superseded historical record.
- Pre-registration corrected before deposit: the T1 harmonisation rule
  (per-sample z-scoring), two-sided Clopper-Pearson bounds, the T4 tie rule and
  the cohort facts.
- `requests` is declared as a dependency; the download scripts import it.
- Ruff settings moved to `[tool.ruff.lint]`, so ruff no longer warns.
- `gliostrata ingest` and `data/ingest.py` point to the real download scripts
  instead of the roadmap.
- README, `docs/roadmap.md`, `docs/reproducibility.md` and
  `docs/mgmt_idh_calling.md` brought back in line with the code.
- CI partition guardrail: a manifest is active because `ACTIVE_MANIFESTS` lists
  it, not because it is the one manifest that is not v1. A new sealed manifest
  is added to that list in the commit that adds it, so the external partition
  will not be mislabelled "superseded".

### Changed — conditions of the independent review of `b83f9f7` (2026-11-27)
- Read guard hardened: a value read that cites a scoring rule must name a rule id
  from the committed registry `configs/scoring_rules.yaml`; the rule must be
  registered for that analysis and name a frozen scoring-set manifest whose
  checksum still matches; and every input's checksum is recomputed at read time.
  Free text no longer opens the partition gate. Ledger entries that are
  observations or notes can be appended and are never counted as reads.
- The unexplained access to the held-out cohort is identified and recorded: the
  review agent's own workspace-wide text search opened five extracted CGGA text
  files on 2026-11-27 (access ledger entries 8 and 9). No value was used.
  `.clineignore`, `.rgignore` and `.clinerules` keep workspace-wide tools away
  from the two CGGA directories.
- Fold-wise imputation: `fit_clinical_baseline` now fits each cross-validation
  fold's design on that fold's training rows, as its documentation said. The TCGA
  C-index is unchanged (0.673544; 0.665650 without the IDH-mutant patients) and
  the design matrices are byte-identical
  (`tests/fixtures/tcga_fit_fingerprint.json`).
- T2 account corrected: the KPS imputation of variant 1 is stated with its counts
  wherever its fit was shown as 594; the Brier excess is no longer attributed to
  an offset; the internal C-index contrasts are labelled descriptive, with no
  interval.
- CI: the `test` job installs `requirements-lock.txt`; a `test-floor` job runs the
  suite on Python 3.10 from `requirements-floor.txt`; the partition-drift logic is
  `scripts/ci_partition_drift.py`, which tests execute.
- Calibration slope, intercept and decision-curve analysis: decided not to be
  done for this thesis (`docs/cgga_validation.md`, section 12.5).

### Added — stage R1
- `requirements-lock.txt` (37 distributions, Python 3.14.3) and
  `requirements-floor.txt`.
- `scripts/rerun_tcga_internal_cv.py` and its tests; the rerun record
  `data/results/tcga_gbm/foldwise_imputation_rerun.json` (git-ignored).
- `data/provenance/file_manifest_checksums.json`: sizes, row counts and SHA-256 of
  the per-modality GDC file manifests, which carry patient barcodes and are so
  not committed.

### Added — stage 3b: the feature-store contract and the subtype classifier (2026-11-27)
- `gliostrata harmonize-clinical` (`data/clinical.py`): the L1 TCGA clinical table with
  an explicit owned-column list, written to `data/interim/tcga_gbm_clinical_harmonized.parquet`,
  and a separate label table (`data/interim/tcga_gbm_labels.parquet`). It raises on an L2
  or molecular column name that carries a value and on an unknown column; every missing
  value is traced to the raw GDC record for its reason code; nothing is imputed. It also
  writes a provenance sidecar, `data/interim/tcga_gbm_clinical_harmonized_provenance.json`,
  with the SHA-256 of its inputs and outputs; real mode stops without a matching one.
- `gliostrata build-feature-store --mode real|synthetic` (`features/build.py`,
  `features/contract.py`, `features/views.py`): one command for L2 and L3. The store
  holds every cohort patient, one row each, and never filters. It is four files: the
  features, the labels (outcomes are labels, never features), a missing-value reason
  table, and a column-metadata JSON (role, source layer, source file SHA-256, class).
  Model-ready populations are named views (`profiling`, `transferable_model`,
  `molecular_model`), each with its flow counts and per-split counts for
  `tcga_gbm_clinical_v2`. Synthetic mode stays the default and is labelled.
- `profiling/centroid.py`: single-sample subtype assignment, nearest centroid (Pearson)
  over the published Verhaak 2010 840-gene signature (831 genes resolved), each sample's
  log2(TPM+1) z-scored across its own genes. Two frozen centroid sets with checksums in
  `configs/subtype_centroids/`: `verhaak840_tcga_v2train` (the store's per-patient
  subtype) and `verhaak840_tcga_all` (reserved for T1). Training labels: Brennan et al.
  2013; sources and checksums in `data/provenance/subtype_reference_sources.json`.
- `scripts/prepare_subtype_reference.py`, `scripts/build_subtype_centroids.py`,
  `scripts/subtype_scoring_set_dependence.py`, `scripts/compare_subtype_labels.py`.
- `.gitattributes`: the centroid JSON files are byte-exact (no line-ending conversion),
  so their checksums hold on any checkout.
- Tests: `tests/test_feature_store_contract.py`, `tests/test_centroid.py`,
  `tests/test_subtype_scripts.py`.

### Changed — stage 3b (2026-11-27)
- Pre-registration sections 3.3 and 4.5 no longer name `cgga_external_v1` as T1's test
  partition: it is the primary external test for T2, T4 and T5, and T1 is scored on the
  scoring sets of section 2.2 items 9 and 10. Logged as a dated amendment (section 8);
  no threshold, bound, rule, scoring set or partition changed.
- `configs/scoring_rules.yaml`: the T1 rule's description now says what item 9(d)
  registers (each sample z-scored across its own genes), not "per-gene z-scored". Wording
  only; the read guard behaves as before.
- `build_feature_store` takes each input by explicit name and raises on a column
  supplied by more than one input, a null-only column that shadows a computed feature,
  a column its input does not own (an outcome label among them), a stub or
  cohort-relative subtype input in real mode, and a missing value without a reason.
  The stage 3 script's placeholder-dropping workaround is removed.
- Absent never becomes a value: missing resection is `not_available_in_cohort`, never
  "biopsy"; a missing molecular call is no call, never 0; call status has three states
  (called, not tested, tested but not callable); the continuous MGMT score and the
  threshold used are kept; KPS is never imputed in the store. `treatment_arm` is a
  feature, not an outcome.
- Real mode accepts only the centroid set `verhaak840_tcga_v2train`. A sample with no
  finite centroid correlation gets no subtype (`scoring_failed`). A patient sequenced
  only on a non-primary sample is `tested_not_callable` (from the GDC file manifests; one
  IDH case).
- The per-patient subtype is the single-sample call; the cohort-relative
  `subtype_prob_*` columns leave the real store. `optimize` and `generate-protocol`
  refuse the real store rather than encode an absent value as a value.
- `scripts/build_tcga_feature_store.py` is a thin wrapper over
  `gliostrata build-feature-store --mode real`. The DVC stage `build_feature_store` and
  the Snakemake rule (with `source=tcga`) run real mode with the expression input,
  its annotation and the two molecular file manifests; a new `harmonize_clinical` stage and rule precede them.

### Fixed
- `test_argmax_matches_input_subtype` computed an argmax and asserted nothing
  about it; it now asserts the match rate.
- Lint backlog: unused imports, unsorted import blocks, an `== False`
  comparison and an unused variable.

### Added — repository scaffold (v0.1.0)
- Project structure with `src/gliostrata` package layout (L0–L7 modules).
- Click-based CLI (`gliostrata` command) with smoke-test, synthetic-data,
  pipeline-dry-run, and per-layer stub commands.
- Pytest test suite with synthetic-data fixtures; all tests pass on a clean install.
- Configuration files (`configs/*.yaml`) for every layer.
- Snakemake workflow scaffold (`workflow/Snakefile` + `rules/`).
- Documentation: `README.md`, `docs/protocol.md`, `docs/architecture.md`.
- DVC pipeline definition (`dvc.yaml`).
- MIT license.

### Notes
- L0–L7 modules are **stubs** with realistic function signatures and synthetic-data
  implementations where feasible. The actual scientific implementations
  (TCGA ingestion, Verhaak classifier, PhysiCell simulation, BoTorch optimizer,
  DeepSurv survival model, protocol generator) are TODO — see `docs/roadmap.md`.
- Heavy optional dependencies (scanpy, lifelines, botorch, physiCell) are NOT
  required to run the scaffold or the tests.
