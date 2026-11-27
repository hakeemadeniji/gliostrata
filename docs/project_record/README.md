# GlioStrata project record

**As of 2026-11-26.** This folder is a complete, self-contained record of what the
GlioStrata project is, what it has done, what it has found, and what remains. It is
written for two readers: the project lead, and a researcher who has not seen the
repository before.

## Files

| File | What it holds |
| --- | --- |
| `01_project_overview.md` | The research question, hypotheses, the eleven registered thresholds, the pipeline layers and their status, the data in hand, and the software facts. |
| `02_findings.md` | Every result measured so far, with the method, the numbers, the interpretation and the limits of each. |
| `03_work_log.md` | Every task run, in order, with its evidence and outcome. |
| `04_commit_history.md` | All 27 commits on `master`, with the files each one changed, and the uncommitted work. |
| `05_future_plan.md` | The tests not yet run and what each needs, the steps before the OSF deposit, the implementation strategy per layer, open decisions and known gaps. |
| `06_evidence_reports.md` | One report per stage: changed files, commands and results, test counts, artifact checksums, external reads, unresolved issues, and whether the next gate is satisfied. |

## How to read it

- **Where a number appears, its source is named.** Numbers were taken on 2026-11-26 from
  the result files under `data/results/`, the sealed manifests under `data/partitions/`,
  `git log`, the test runner, or the document cited. None is quoted from memory.
- **Three outcome labels are used**, as the pre-registration defines them. **FAIL**: the
  point estimate missed its threshold. **INCONCLUSIVE**: the point estimate cleared but
  the decisive confidence bound did not. **NOT EVALUABLE**: the criterion cannot be
  computed on the available data. A criterion is **met** only when the point estimate and
  its decisive bound both clear.
- **"Internal"** means measured on TCGA-GBM by cross-validation. It is context, not a
  result about a new cohort. **"External"** means measured on CGGA.

## State of the project in ten lines

1. GlioStrata is a pre-registered, computational framework that aims to recommend a
   stem-cell-delivered therapeutic configuration per glioblastoma patient.
2. The code runs end to end on synthetic data. Real public-cohort ingestion (L0) works
   for TCGA-GBM and CGGA. Layers L0 to L3 have also been run on the real TCGA data as a
   standalone script, `scripts/build_tcga_feature_store.py` (155 of 594 patients).
   Several scientific modules exist but are not wired into the
   pipeline; others are stubs.
3. Two development cohorts are sealed: 594 TCGA-GBM patients (425 train / 169 test).
4. CGGA, the external cohort, holds 305 unique glioblastoma patients, of whom 218 are
   eligible for the sealed external partition. That partition is **not yet sealed**.
5. Of the eleven registered thresholds, **two have been evaluated and both were not met**:
   T2's external C-index (0.5385 against a threshold of 0.65) and T2's Brier score at
   18 months (0.2617 against 0.25).
6. A third registered criterion, T2's improvement over a 2024 baseline, is **not
   evaluated**: no comparator model has been named.
7. T1, T3, T4 and T5 have **not been run**. T3's input data have not been obtained.
8. The supporting analyses are solid: the IDH sensitivity analysis passed its
   pre-specified criterion, and the molecular calling and subtyping ran on real data.
9. The test suite is 1042 tests: 1041 pass and 1 is skipped. It passes with all network
   sockets blocked.
10. The pre-registration is **drafted but not deposited** on the OSF. The work that was
    uncommitted when this was written is committed (section 4, last table).

## Three things a reader should know before citing anything here

- **The project has no pre-registration yet in the legal sense.** Until the OSF deposit is
  made, T1 to T5 are exploratory analyses, whatever the repository says.
- **CGGA has been read twice, both times for T2, and both times before the external
  partition was sealed.** The pre-registration allows one read per registered analysis,
  and both reads are disclosed in its amendment log. Every access is listed in
   `data/access_ledger.jsonl`.
- **Several layers of the pipeline produce synthetic output.** A "subtype" or "cell-state
  fraction" written by `gliostrata profile` today is derived from the clinical table, not
  measured from expression. `01_project_overview.md` marks every layer.
