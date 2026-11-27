# 4. Commit history

**27 commits on `master`**, the only branch. All are dated 2026-11-26; the project's git
history begins that day. Times are as recorded by git (local time). Newest first.

Source: `git log` on 2026-11-26. "Files" lists what the commit added or changed; the
per-commit file lists were read from `git log --name-status`.

| # | Hash | Time | Commit message | Files |
| --- | --- | --- | --- | --- |
| 27 | `a092a16` | 19:26 | Add the final OSF Sampling-section text | `docs/osf_final_fields.md` (new) |
| 26 | `384bdf1` | 19:26 | Add the CGGA external-partition eligibility script with a dry run | `scripts/cgga_external_eligibility.py`, `tests/test_cgga_eligibility.py` |
| 25 | `5677045` | 19:26 | Add the all-test seal for external cohorts (`--all-test`) | `cli.py`, `evaluation/partitions.py`, `tests/test_cli.py`, `tests/test_evaluation.py` |
| 24 | `ee785f6` | 19:26 | Fix the T1 harmonisation rule and reconcile the pre-registration with the data in hand | `docs/preregistration.md` |
| 23 | `fa159de` | 18:46 | State the T1 accuracy margin and the T4 tie denominator explicitly | `docs/preregistration.md` |
| 22 | `bfd5c0e` | 18:45 | Close the open pre-registration items: T1 guard rails, T3 rule, T4 budget | `docs/preregistration.md` |
| 21 | `857a4cb` | 16:31 | Record T2 external validation result in the pre-registration | `docs/cgga_validation.md`, `docs/preregistration.md` |
| 20 | `7f4bd28` | 16:05 | Add T2 external validation of Cox model on CGGA | `docs/cgga_validation.md`, `scripts/validate_cox_cgga.py`, `tests/test_validate_cox_cgga.py` (all new) |
| 19 | `de23202` | 15:22 | Run the pre-registered IDH sensitivity analysis (section 6.3) | `docs/mgmt_idh_calling.md`; `scripts/idh_sensitivity_cox.py`, `scripts/idh_sensitivity_subtype.py` and their two test modules (new) |
| 18 | `320707c` | 15:21 | Add molecular-calls join for the cohort (L2) | `scripts/join_molecular_calls.py`, `tests/test_join_molecular_calls.py` |
| 17 | `183824f` | 14:37 | Add GP-EI Bayesian configuration optimizer (L6) | `optimization/bayes_opt.py`, `tests/test_bayes_opt.py` |
| 16 | `6bc2766` | 14:34 | Add marker-based state and TME deconvolution (L2) | `profiling/deconvolve.py`, `tests/test_deconvolve.py` |
| 15 | `0eda7a8` | 14:34 | Add CIBERSORTx mixture export for manual portal runs | `scripts/export_cibersortx.py`, `tests/test_export_cibersortx.py` |
| 14 | `1777030` | 14:16 | Add Snakemake workflow, DVC pipeline and GitHub Actions CI | `.dvcignore`, `.github/workflows/ci.yml`, `docs/reproducibility.md`, `tests/test_workflow.py` (new); `.gitignore`, `Makefile`, `workflow/Snakefile` (changed) |
| 13 | `fb247ad` | 14:12 | Add TCGA-GBM MGMT promoter methylation calling and DVC stages (L2) | `docs/mgmt_idh_calling.md`, `scripts/download_tcga_methylation.py`, `tests/test_tcga_methylation.py`, 9 files under `workflow/rules/` (new); `dvc.yaml` (changed) |
| 12 | `90a9f1c` | 13:36 | Add TCGA-GBM IDH mutation calling from GDC MAF (L2) | `scripts/download_tcga_mutations.py`, `tests/test_tcga_mutations.py` |
| 11 | `e361bed` | 09:38 | Add CGGA external validation cohort ingestion (L0) | `docs/cgga_ingestion.md`, `scripts/download_cgga.py`, `tests/test_cgga.py` |
| 10 | `a0e0dd6` | 09:34 | Add signature-based Verhaak subtyping from real expression (L2) | `profiling/verhaak.py`, `tests/test_verhaak.py` |
| 9 | `65b15b4` | 01:50 | Finish updating Cox smoke-test counts for the 594-patient cohort | `tests/test_cox.py` |
| 8 | `8f45fb2` | 01:49 | Update Cox cohort-size canary for the 594-patient cohort | `tests/test_cox.py` |
| 7 | `26be9ce` | 01:48 | Drop only contradictory diagnoses, not missing ones (n 583 to 594) | `data/partitions/tcga_gbm_clinical_v2.json` (new); `scripts/download_tcga_clinical.py`, `tests/test_tcga_download.py` |
| 6 | `239aca3` | 01:40 | Add frozen partitions, threshold registry and pre-registration draft | `data/partitions/tcga_gbm_clinical_v1.json`, `docs/preregistration.md`, `evaluation/` (3 files), `tests/test_evaluation.py` (new); `cli.py` |
| 5 | `82ee3e5` | 01:40 | Add TCGA-GBM RNA-seq expression ingestion (L0) | `scripts/download_tcga_expression.py`, `tests/test_tcga_expression.py` |
| 4 | `8b9264e` | 01:37 | Add fitted Cox PH model (L5) alongside the survival stub | `prediction/cox.py`, `tests/test_cox.py` |
| 3 | `cb98002` | 01:34 | Fix two critical statistical bugs on the pre-registered claim path | `utils/stats.py`; `tests/test_stats.py` (new) |
| 2 | `e627203` | 01:22 | Untrack generated data and Office lock file | `.gitignore`; 36 files removed |
| 1 | `0b5da09` | 01:22 | Baseline: GlioStrata L0-L7 scaffold with TCGA clinical ingestion | 107 files (the whole repository) |

Module paths are abbreviated: `profiling/`, `prediction/`, `optimization/`, `evaluation/`
and `utils/` are under `src/gliostrata/`.

## Notes on the history

- **Two commits were amended** before the final history. The editor's checkpoint refs record
  earlier hashes `e28ba12` (IDH calling) and `0115f98` (workflow, DVC and CI) for those
  subjects. The hashes in the table are the ones on `master`.
- **The repository also holds 14 editor checkpoint refs** under `refs/cline/checkpoints/`.
  These are automatic snapshots, not commits on the branch. Counting every commit
  reachable from all refs gives 70, which is not the number of commits on `master`.
- **No commit records the external evaluation results as data.** The result files are
  git-ignored; the commits record the scripts, the tests and the written account.

## Uncommitted work (at the time of writing)

Everything in Phases G and H of the work log is in the working tree, not in a commit.

**Measured on 2026-11-26 with `git status` and `git diff --shortstat`.**

| Group | Files | Lines |
| --- | --- | --- |
| Tracked files modified | 33 | 590 insertions, 291 deletions |
| of which code, tests and packaging | 22 | 72 insertions, 66 deletions |
| of which documentation, CI and Makefile | 11 | the remainder |
| Untracked, new | 8 | not yet counted by git |

**The 22 code, test and packaging files.** `pyproject.toml`, `requirements.txt`,
`environment.yml`; `src/gliostrata/cli.py`, `data/ingest.py`, `data/harmonize.py`,
`features/store.py`, `optimization/bayesian.py`, `optimization/config_space.py`,
`prediction/survival.py`, `simulation/nsc.py`, `utils/stats.py`; five download and
pipeline scripts; `tests/test_data.py`, `test_features.py`, `test_optimization.py`,
`test_profiling.py`, `test_simulation.py`. Most of these are lint fixes of one or two
lines. The substantive ones are `cli.py` and `data/ingest.py` (messages that now point to
the real download scripts), the three packaging files (`requests` declared, Python floor,
ruff settings), `optimization/bayesian.py` (dead line removed) and
`tests/test_profiling.py` (the vacuous test).

**The 11 documentation and CI files.** `README.md`, `CHANGELOG.md`, `Makefile`,
`.github/workflows/ci.yml`, and `docs/architecture.md`, `cgga_validation.md`,
`mgmt_idh_calling.md`, `osf_final_fields.md`, `preregistration.md`, `reproducibility.md`,
`roadmap.md`.

**The 8 untracked files.** `scripts/validate_brier_cgga.py`,
`tests/test_validate_brier_cgga.py`, and the six files of this record.

**Two things follow.**
- The committed code does not contain the Brier evaluation, and the committed
  pre-registration does not record its result. A reader of the repository at `a092a16`
  sees T2's C-index result only.
- The registered commit hash for the OSF deposit should be taken after these changes are
  committed.
