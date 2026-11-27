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

## Phase G. Full project audit and correction (after 19:26; **uncommitted**)

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

## Phase H. T2 Brier row and this record (2026-11-26; **uncommitted**)

| Task | Order of work | Outcome |
| --- | --- | --- |
| Evaluate T2's Brier-score row | Specification written into `scripts/validate_brier_cgga.py` and 20 synthetic-data tests; tests passed; dry run (counts only); then the single evaluation | **CGGA read 2.** Brier 0.2617 (0.2391 to 0.2829): not met |
| Record the result | Validation document section 12 and access record; pre-registration subsection 2.3 and amendment log; "read once" corrected to two reads; OSF text, README, roadmap, CHANGELOG | Six documents updated |
| Verify | Full suite with sockets blocked; strict lint | 984 passed, 1 skipped, 0 failed; blocking lint clean |
| Write this record | Every number taken from a named file or command | `docs/project_record/` |

**A bug the tests caught before any data were read.** The first version of the Brier
script unpacked a four-value result into two names. The synthetic-data tests failed on it,
it was fixed, and only then was the script run on CGGA.

## Reads of the external cohort, in total

| # | When | What was read | What was computed |
| --- | --- | --- | --- |
| 0 | Before the first script | Column names, dtypes, null counts, value counts, survival range | Nothing |
| 1 | 2026-11-26, Phase E | CGGA clinical cohort | T2 C-index row, once |
| 2 | 2026-11-26, Phase H | CGGA clinical cohort | T2 Brier row, once |

Dry runs before each read printed availability counts only. The CGGA expression values
have never been read by a model: eligibility was counted from column names.
