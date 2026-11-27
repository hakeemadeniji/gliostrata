# 3. Work log

Every task the project has run, in order. Times are the commit times recorded by git
(local time); the whole committed history falls on 2026-11-26. Tasks without a commit are
marked **uncommitted** and cite the file or output that evidences them.

Tasks that read the external cohort are marked **CGGA read**, because the number of such
reads is a pre-registration commitment.

## Phase A. Scaffold, statistics and the Cox model (01:22 to 01:50)

| Time | Task | Outcome | Commit |
| --- | --- | --- | --- |
| 01:22 | Repository scaffold: package layout L0 to L7, CLI, configs, tests, docs, DVC and Snakemake definitions | 107 files; pipeline runs end to end on synthetic data | `0b5da09` |
| 01:22 | Remove generated data and an Office lock file from git | 36 files removed (35 under `data/manual/`, one lock file); `.gitignore` updated | `e627203` |
| 01:34 | Fix two statistical bugs on the path to the pre-registered claim | `utils/stats.py` fixed; `tests/test_stats.py` added | `cb98002` |
| 01:37 | Fitted Cox proportional-hazards model (Breslow partial likelihood) alongside the stub | `prediction/cox.py`, `tests/test_cox.py` | `8b9264e` |
| 01:40 | TCGA-GBM bulk RNA-seq ingest | `scripts/download_tcga_expression.py` | `82ee3e5` |
| 01:40 | Threshold registry, frozen partitions, `freeze-partition` command, pre-registration draft | `evaluation/` package; `docs/preregistration.md` | `239aca3` |
| 01:48 | Correct the clinical ingest filter: drop only contradictory diagnoses | Cohort 583 to 594; new partition `v2` frozen | `26be9ce` |
| 01:49, 01:50 | Update the Cox smoke-test counts for the 594-patient cohort | `tests/test_cox.py` | `8f45fb2`, `65b15b4` |

**The cohort change, and why it matters.** The filter first dropped patients with a missing
diagnosis; 11 such patients were reinstated because enrolment in TCGA-GBM establishes the
disease. The first sealed partition (583 patients) then no longer verified against the live
cohort. It was kept unmodified as the historical record, and a new manifest (594) was
frozen under a new name. The partition guardrail caught a real change, which is the only
way to know a guardrail works.

## Phase B. Subtyping and the external cohort (09:34 to 09:38)

| Time | Task | Outcome | Commit |
| --- | --- | --- | --- |
| 09:34 | Signature-based Verhaak subtype scoring from real expression | `profiling/verhaak.py`, 32 tests | `a0e0dd6` |
| 09:38 | CGGA external cohort ingest: three clinical files and two RNA-seq count matrices | `scripts/download_cgga.py`; `docs/cgga_ingestion.md` | `e361bed` |

**Established during ingest** (schema-level inspection, no outcome related to a feature):
CGGA needs no registration; its 1,319 clinical rows are 333 primary GBM, which are 305
patients; it publishes no performance status. **CGGA read**: schema only.

## Phase C. Molecular calling, workflow and optimiser (13:36 to 14:37)

| Time | Task | Outcome | Commit |
| --- | --- | --- | --- |
| 13:36 | IDH mutation calling from open GDC mutation files | `scripts/download_tcga_mutations.py`; real run: 373 tested, 23 mutant | `90a9f1c` |
| 14:12 | MGMT promoter methylation calling, with documentation | `scripts/download_tcga_methylation.py`; `docs/mgmt_idh_calling.md`; 9 Snakemake rule files | `fb247ad` |
| 14:16 | Snakemake workflow, DVC pipeline definition, GitHub Actions CI | `ci.yml`, `.dvcignore`, `docs/reproducibility.md`, `tests/test_workflow.py` | `1777030` |
| 14:34 | CIBERSORTx mixture export for manual portal runs | `scripts/export_cibersortx.py` | `0eda7a8` |
| 14:34 | Marker-based cell-state and TME deconvolution | `profiling/deconvolve.py`, 61 tests | `6bc2766` |
| 14:37 | Gaussian-process / expected-improvement optimiser | `optimization/bayes_opt.py`, 67 tests | `183824f` |

Two commits in this phase carry earlier hashes in the editor's checkpoint references
(`e28ba12`, `0115f98`); they were amended before the final history. The hashes above are
the ones on `master`.

## Phase D. Molecular join and the IDH sensitivity analysis (15:21 to 15:22)

| Time | Task | Outcome | Commit |
| --- | --- | --- | --- |
| 15:21 | Join IDH and MGMT calls onto the cohort | `scripts/join_molecular_calls.py` | `320707c` |
| 15:22 | Run the pre-specified IDH sensitivity analysis | C-index moves 0.0079 (under the 0.01 limit); subtype shares move at most 1.10 pp | `de23202` |

## Phase E. T2 external validation (16:05 to 16:31)

| Time | Task | Outcome | Commit |
| --- | --- | --- | --- |
| 16:05 | T2 external validation of the Cox model on CGGA: script, 24 tests (synthetic and TCGA only), then the single evaluation | **CGGA read 1.** C-index 0.5385 (0.4994 to 0.5807): not met. Committed baseline not evaluable (no KPS) | `7f4bd28` |
| 16:31 | Record the T2 outcome in the pre-registration | New subsection 2.3; first amendment-log entry | `857a4cb` |

**Order of work in that task**, so the number can be trusted: schema inspection of CGGA
(column names, null counts, value counts, survival range), then the script and its tests,
then a dry run (counts only), then the single evaluation. No edit to the script after its
final run.

## Phase F. Closing the pre-registration (18:45 to 19:26)

| Time | Task | Outcome | Commit |
| --- | --- | --- | --- |
| 18:45 | Close open pre-registration items: T1 guard rails, T3 rule and dates, T4 budget and tie rule, partition sealing order | `docs/preregistration.md` | `bfd5c0e` |
| 18:46 | State the T1 accuracy margin and the T4 tie denominator explicitly | `docs/preregistration.md` | `fa159de` |
| 19:26 | Fix the T1 harmonisation rule; reconcile the pre-registration with the data in hand | `docs/preregistration.md` | `ee785f6` |
| 19:26 | All-test seal for external cohorts (`freeze-partition --all-test`) | `cli.py`, `partitions.py`, tests | `5677045` |
| 19:26 | Eligibility listing for the external partition, with a dry run | `scripts/cgga_external_eligibility.py`, 8 tests. A dry run counted 218 eligible patients from identifiers only; the table itself has not been written | `384bdf1` |
| 19:26 | Final paste-ready text for the four OSF Sampling fields | `docs/osf_final_fields.md` | `a092a16` |

## Phase G. Full project audit and correction (after 19:26; committed in `0c20999` and `10cfd93`)

A whole-project audit was run in five passes, and everything clearly broken was fixed.

| Task | Method | Outcome |
| --- | --- | --- |
| Inventory and stale-document scan | Every README and docs claim checked against the code | List of stale claims produced |
| Project gates | Partition verification, pipeline dry run, workflow tests | Partition v2 verifies (594 patients, 425 / 169); v1 drift is expected |
| Lint | `ruff` under the project config | Deprecation warning removed; backlog 182 to 148 |
| Dependency check | Declared dependencies against actual imports; `pip check` | `requests` was imported but undeclared: added to three files |
| Python-floor check | Syntax parse at the declared minimum | No use of features above the minimum |
| Test run with sockets blocked | Full suite with outbound network disabled | 964 passed, 1 skipped, 0 failed, 0 blocked (before the Brier tests) |

**Fixes made.**

- **A vacuous test.** `test_argmax_matches_input_subtype` computed an argmax and asserted
  nothing. The strengthened assertion exposed a latent bug (a column name compared with a
  subtype). It now asserts a match rate of at least 80%.
- **Dead code** removed in `optimization/bayesian.py`.
- **Packaging.** `requests` declared; ruff settings moved to the current section so the
  deprecation warning is gone; `environment.yml` changed from a pinned to a floor Python
  version.
- **Messages that pointed nowhere.** The CLI and ingest module now name the real download
  scripts.
- **Stale documentation rewritten**: README, roadmap, reproducibility, architecture,
  molecular-calling document, CHANGELOG, CI and Makefile comments. Test counts, layer
  statuses and "nine DVC stages" were re-derived from the code.
- **OSF text.** The "starting and stopping rules" block rewritten against the
  pre-registration; the T1 accuracy bounds recomputed as two-sided exact intervals, which
  superseded the one-sided figures in the decision memo.

## Phase H. T2 Brier row and this record (2026-11-26; committed in `0c20999`, `10cfd93` and `8f79b33`)

| Task | Order of work | Outcome |
| --- | --- | --- |
| Evaluate T2's Brier-score row | Specification written into `scripts/validate_brier_cgga.py` and 20 synthetic-data tests; tests passed; dry run (counts only); then the single evaluation | **CGGA read 2.** Brier 0.2617 (0.2391 to 0.2829): not met |
| Record the result | Validation document section 12 and access record; pre-registration subsection 2.3 and amendment log; "read once" corrected to two reads; OSF text, README, roadmap, CHANGELOG | Six documents updated |
| Verify | Full suite with sockets blocked; strict lint | 984 passed, 1 skipped, 0 failed; blocking lint clean |
| Write this record | Every number taken from a named file or command | `docs/project_record/` |

**A bug the tests caught before any data were read.** The first version of the Brier
script unpacked a four-value result into two names. The synthetic-data tests failed on it,
it was fixed, and only then was the script run on CGGA.

## Phase I. Secure baseline and read discipline (2026-11-26; after the record was written)

A new session reviewed and committed the work above, then built the read discipline.
Commands and results are in `06_evidence_reports.md`, stages 1 and 2.

| Task | Method | Outcome | Commit |
| --- | --- | --- | --- |
| Review the uncommitted work | Every diff read; secret scan of the diff and the untracked files | Nothing to change: no credential, no patient-level file, no cache, and no change to any threshold, bound, rule or manifest | |
| Verify before committing | Full suite; blocking lint; partition v2; registry checksum | 984 passed, 1 skipped; lint clean; v2 verified (594, 425 / 169); checksum unchanged | |
| Commit the working tree | Three commits, in the split the future plan suggested | Code, tests and packaging; documentation; this record | `0c20999`, `10cfd93`, `8f79b33` |
| Check commit 28 on its own | Throwaway worktree; `test_workflow.py` and `test_cli.py` | 78 passed, 1 skipped; worktree removed. The package code it ran was the main checkout's, not the commit's (corrected in stage R1) | |
| Checksum the CGGA files | PowerShell `Get-FileHash`; file bytes only | **CGGA access: checksum only.** Four files hashed; nothing deserialised | |
| Read discipline | Access ledger; read guard; registration state; CLI command; 31 tests | 7 ledger entries, 2 of them value reads. A value read is refused until the registration is deposited and the partition sealed | `324a09b` |
| CI partition guardrail | `ACTIVE_MANIFESTS` replaces "every manifest but one is superseded" | Run locally in three scenarios; a listed manifest that does not exist stops the job | `324a09b` |
| Full suite with sockets blocked | Python-level socket entry points replaced before pytest | 1015 passed, 1 skipped; 0 blocked connection attempts | |
| Documentation | Stale numbers corrected (lint count, test counts); a broken README sentence repaired; T1 findings recorded (section 5.1.1 of the future plan) | `f3c2a4b` |
| Real-data slice, L0 to L3 | `scripts/build_tcga_feature_store.py` on the TCGA inputs; no network, no CGGA, no outcome value used | 155 of 594 patients in the store; every exclusion counted. Found: resection extent absent for the whole cohort; KPS absent for 64 of the 283 with expression; cohort placeholder columns would shadow the L2 outputs | The commit after `f3c2a4b` |

## Phase J. Close the conditions of the independent review (2026-11-27; stage R1)

An independent review of commit `b83f9f7` returned "accept with conditions": 0 blockers,
6 major findings and 7 minor ones. This stage closes them. Commands and results are in
`06_evidence_reports.md`, stage R1. TCGA only: no CGGA value was read, hashed or
downloaded, and the only model fitting was the fold-wise imputation rerun.

| Task | Method | Outcome | Commit |
| --- | --- | --- | --- |
| Identify the unexplained CGGA access | File metadata only; then the agent's own session log, the access times of files across the repository, and editor and agent logs | It was the review agent's own workspace-wide text search at the start of the stage 1 review: two path-regex queries, 30 seconds, no result. Five extracted text files were scanned for path strings; nothing was deserialised or returned. Ledger entries 8 (observation), 9 (identification) and 10 (times of day; a correction to entry 6) | First commit of stage R1 |
| Prevent recurrence | `.clineignore`, `.rgignore`, `.clinerules`; a rule in `docs/reproducibility.md` | Both CGGA directories are excluded from workspace-wide tools. Removal of the plain-text copies is deferred, with its reason (`05_future_plan.md`, section 5.5) | Same |
| Harden the read guard | A scoring-rule registry (`configs/scoring_rules.yaml`); input checksums recomputed at read time; notes in the ledger | Free text no longer opens the partition gate. 11 more tests | Same |
| Brier explanation and C-index contrasts | `docs/cgga_validation.md`, `02_findings.md` and the pre-registration rewritten | Two possible causes, offset and dispersion, that the data on record cannot tell apart; contrasts labelled descriptive, with no interval | Documentation commit |
| KPS imputation disclosure | The counts stated wherever variant 1's fit was shown as 594; a committed fixture pins them | 594 passed, 592 modelled, 150 imputed at 80, 22 missing treatment arm; `tests/fixtures/tcga_fit_fingerprint.json` | Same |
| Fold-wise imputation | A design fitted per training fold (`fit_design_plan`, `apply_design_plan`); rerun on TCGA | C-index unchanged: 0.673544 and 0.665650 before and after. The three design matrices are byte-identical to before | Same |
| Lock file, floors, CI | `requirements-lock.txt` (37 distributions); `requirements-floor.txt`; CI jobs | `test` installs from the lock file; `test-floor` is new (Python 3.10). Neither has run on GitHub | Same |
| Provenance manifests | Schemas inspected first | They carry patient barcodes, so their sizes, row counts and checksums are committed (`data/provenance/`), not the files | Same |
| Calibration | A decision, recorded | Not done for this thesis (`docs/cgga_validation.md`, section 12.5) | Same |
| CI partition job | Logic moved to `scripts/ci_partition_drift.py`; 10 tests execute it | Missing, superseded, external all-test and tampered manifests | Same |
| Isolated clean-checkout check | Fresh clone; `PYTHONPATH=<clone>/src`; `gliostrata.__file__` printed, in the shell and inside the pytest process | `gliostrata` resolved inside the clone. The first run, at `441c94a`, found a test that failed in every clean checkout (it assumed the git-ignored TCGA cohort was on disk); fixed. The second run, in a new clone at `9216e6f`: 1052 passed, 23 skipped (the real-data tests and the Snakemake test) | `9216e6f` |
| Finish | Full suite with sockets blocked and a CGGA file-access audit hook; blocking and advisory lint | 1074 passed, 1 skipped; 0 blocked connections; no CGGA access; blocking lint clean; advisory lint 170 findings | Tag `baseline-r1` |

## Phase K. The feature-store contract and the subtype classifier (2026-11-27; stage 3b)

TCGA only: no CGGA value was read, downloaded or hashed, no model was fitted on an outcome,
and the pre-registration was not edited. Population rules were fixed from availability
counts. Work ran as three parallel workstreams (a read-only reconciliation, the subtype
classifier, the feature store), then an independent evaluation. Commands and results are
in `06_evidence_reports.md`, stage 3b.

| Task | Method | Outcome | Commit |
| --- | --- | --- | --- |
| Reconcile what is closed (Part A) | Read-only: code, committed docs, the ledger, the TCGA raw files; live GDC queries for the KPS and resection trace | One page, eight items. `cgga_external_v1` not created; pre-registration 4.5 wrong about T1; every missing KPS absent at source; resection absent at source; the 284th sample's case is empty; the v2 test split was in every internal training fold | `45fdc2d`, `2e19939` |
| Decisions D1 to D6 | Entered as decided, with reasons | `05_future_plan.md` section 5.4.1, including the scoring-set decision for the v2 split | The record commit |
| Scoring-set dependence of the subtype (D6a) | `scripts/subtype_scoring_set_dependence.py`; seed 20261127, 20 random halvings | Cohort-relative: 2 to 18 of 283 calls change; 5 to 7 of the 45 v2-test store patients change when scored alone. Single-sample: 0 in every experiment | `45fdc2d` |
| Single-sample classifier (D6b) | `profiling/centroid.py`; published Verhaak 840-gene signature, 831 resolved; log2(TPM+1) z-scored within each sample; Pearson nearest centroid | Two frozen sets with checksums: `verhaak840_tcga_v2train` (170 training patients) and `verhaak840_tcga_all` (246) | `45fdc2d` |
| Published labels (D6c) | Brennan 2013 labels from the GDC publication page; Verhaak 2010 and Ceccarelli 2016 as cross-checks | Agreement on the v2 test split: single-sample 56 of 76 (0.737), cohort-relative 43 of 76 (0.566). Descriptive | `45fdc2d` |
| Clinical harmonisation (D5) | `gliostrata harmonize-clinical`; owned columns; reasons traced to the raw GDC record | 594 patients; a separate label table; nothing imputed | `2e19939` |
| Feature store and views (D1 to D4, C1 to C5) | `gliostrata build-feature-store --mode real`; DVC stage and Snakemake rule; the stage 3 script made a thin wrapper | 594 patients in the store; views 283, 570 and 298 with per-split counts; repeat run byte-identical | `2e19939` |
| Deconvolution report (C6) | The default state and TME scoring on the 283 | No numerical failure; AC poorly identified; cohort-relative (6 of 50 highest states change when scored alone). Not validation | |
| Neftel atlas (Part D) | File listing of `data/`; public sources looked up | Not on disk. Stopped; nothing downloaded | |
| Verify | Full suite with sockets blocked; blocking and advisory lint | 1113 tests: 1112 passed, 1 skipped; 0 blocked connections; no CGGA access. Blocking lint clean; advisory lint 166 | |
| Correct pre-registration 3.3 and 4.5 (after the stage, at the project lead's request) | The text proposed in Part A item 2, applied; a dated amendment in section 8 | `cgga_external_v1` is the primary external test for T2, T4 and T5; T1 is scored on its own registered scoring sets. No threshold, bound, rule, scoring set or partition changed; threshold-table checksum unchanged | `16b6c09` |
| Align the scoring-rule registry with rule 9(d) (after the stage, at the project lead's request) | `configs/scoring_rules.yaml`, the T1 rule's description | "Per-gene z-scored" replaced by 9(d)'s wording: each sample z-scored across its own genes. Rule id, analyses and manifest fields unchanged | The commit after `16b6c09` |

## Phase L. T1's state-deconvolution row (2026-11-27; run before deposit)

| Task | Method | Outcome | Commit |
| --- | --- | --- | --- |
| Obtain the atlas | GEO GSM3828672 (Smart-seq2 matrix); the 3CA per-cell annotation, because GEO has no malignant flag or state and the Broad portal refused an anonymous download | Both checksummed; provenance recorded | `f5b9fe6` |
| Specify before running | `docs/t1_deconvolution_spec.md`; script and 13 tests; independent pre-run review (2 must-fix, 3 should-fix, all applied) | Committed before the run | `f5b9fe6` |
| Run once | `python scripts/t1_deconvolution_mixtures.py`, 11:40:07Z to 11:40:45Z, exit 0 | Mean per-state r 0.9521 (CI 0.9430 to 0.9589): PASS, run before deposit. Secondary (proportional) 0.8541 | |
| Record | This record, the changelog, reproducibility | Pre-registration section 5 and the OSF text noted as out of date, not edited | The commit after `f5b9fe6` |

## Reads of the external cohort, in total

| # | When | What was read | What was computed |
| --- | --- | --- | --- |
| 0 | Before the first script | Column names, dtypes, null counts, value counts, survival range | Nothing |
| 1 | 2026-11-26, Phase E | CGGA clinical cohort | T2 C-index row, once |
| 2 | 2026-11-26, Phase H | CGGA clinical cohort | T2 Brier row, once |

Dry runs before each read printed availability counts only. The CGGA expression values
have never been read by a model: eligibility was counted from column names.

`data/access_ledger.jsonl` is now the record of every access, including the three in
this table. Its later entries: a checksum-only access from Phase I (entry 6); a second
one, after the L0 to L3 slice (entry 7); and three entries that are not accesses
(8 to 10). Entry 8 is an observation: the review agent's workspace-wide text search
opened five extracted CGGA text files on 2026-11-27T05:48:55Z. Entry 9 identifies it.
It scanned them for two path strings and returned nothing; no value was used. Entry 10
records that entries 0 to 6 have no time of day and corrects entry 6 (five files were
hashed, not four). The full test suite reads no CGGA file.
