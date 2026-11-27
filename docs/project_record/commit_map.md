# Commit map: hashes before and after the author-email rewrite (2026-11-27)

Before the first push to GitHub, the author and committer email of every commit on `master`
was set to the project address, at the project lead's request. 39 commits had carried a
personal address. Git hashes cover the author line, so **every commit hash changed**.
Nothing else did. The check below was run on all 64 commits:

- **Unchanged:** each commit's file tree (so every file's bytes), its author and committer
  dates, and its full message are identical before and after (verified by script; 0
  mismatches).
- **Changed:** the email, and therefore the hash.

The records in this repository were written before the rewrite and cite the **old** hashes.
They are left as written: the project record (`docs/project_record/`), `CHANGELOG.md`, the
T2 search log, `docs/cgga_validation.md`, the commit messages themselves, and the
append-only access ledger (`data/access_ledger.jsonl`, which gains a note pointing here).
Use this table to resolve any of them. The pre-registration's references were moved to the
new hashes. The tag `baseline-r1` was re-pointed to the rewritten commit (`cff12c2` became
`8baf99a`).

Dates are the author dates, unchanged by the rewrite. The order of commits, and therefore
every "committed before" claim in the record, is preserved.

| Old | New | Author date | Subject |
| --- | --- | --- | --- |
| `0b5da09` | `2893052` | 2026-11-26T01:22:01-08:00 | Baseline: GlioStrata L0-L7 scaffold with TCGA clinical ingestion |
| `e627203` | `11216ca` | 2026-11-26T01:22:12-08:00 | Untrack generated data and Office lock file |
| `cb98002` | `c6f27eb` | 2026-11-26T01:34:28-08:00 | Fix two critical statistical bugs on the pre-registered claim path |
| `8b9264e` | `01e3f5c` | 2026-11-26T01:37:18-08:00 | Add fitted Cox PH model (L5) alongside the survival stub |
| `82ee3e5` | `bbec583` | 2026-11-26T01:40:23-08:00 | Add TCGA-GBM RNA-seq expression ingestion (L0) |
| `239aca3` | `dc37e6c` | 2026-11-26T01:40:24-08:00 | Add frozen partitions, threshold registry and pre-registration draft |
| `26be9ce` | `524b480` | 2026-11-26T01:48:51-08:00 | Drop only contradictory diagnoses, not missing ones (n 583 -> 594) |
| `8f45fb2` | `71ebb0b` | 2026-11-26T01:49:40-08:00 | Update Cox cohort-size canary for the 594-patient cohort |
| `65b15b4` | `65a84f6` | 2026-11-26T01:50:47-08:00 | Finish updating Cox smoke-test counts for the 594-patient cohort |
| `a0e0dd6` | `13a7a4f` | 2026-11-26T09:34:00-08:00 | Add signature-based Verhaak subtyping from real expression (L2) |
| `e361bed` | `de8d435` | 2026-11-26T09:38:01-08:00 | Add CGGA external validation cohort ingestion (L0) |
| `90a9f1c` | `f997091` | 2026-11-26T13:36:45-08:00 | Add TCGA-GBM IDH mutation calling from GDC MAF (L2) |
| `fb247ad` | `56e4347` | 2026-11-26T14:12:01-08:00 | Add TCGA-GBM MGMT promoter methylation calling and DVC stages (L2) |
| `1777030` | `85b87e4` | 2026-11-26T14:16:15-08:00 | Add Snakemake workflow, DVC pipeline and GitHub Actions CI |
| `0eda7a8` | `cc5bfac` | 2026-11-26T14:34:20-08:00 | Add CIBERSORTx mixture export for manual portal runs |
| `6bc2766` | `31f408e` | 2026-11-26T14:34:57-08:00 | Add marker-based state and TME deconvolution (L2) |
| `183824f` | `cb25513` | 2026-11-26T14:37:10-08:00 | Add GP-EI Bayesian configuration optimizer (L6) |
| `320707c` | `6de90e7` | 2026-11-26T15:21:36-08:00 | Add molecular-calls join for the cohort (L2) |
| `de23202` | `b35e37d` | 2026-11-26T15:22:06-08:00 | Run the pre-registered IDH sensitivity analysis (section 6.3) |
| `7f4bd28` | `56cc36b` | 2026-11-26T16:05:27-08:00 | Add T2 external validation of Cox model on CGGA |
| `857a4cb` | `e7b2fee` | 2026-11-26T16:31:58-08:00 | Record T2 external validation result in the pre-registration |
| `bfd5c0e` | `f4dd592` | 2026-11-26T18:45:46-08:00 | Close the open pre-registration items: T1 guard rails, T3 rule, T4 budget |
| `fa159de` | `25c6230` | 2026-11-26T18:46:38-08:00 | State the T1 accuracy margin and the T4 tie denominator explicitly |
| `ee785f6` | `a952c0e` | 2026-11-26T19:26:20-08:00 | Fix the T1 harmonisation rule and reconcile the pre-registration with the data in hand |
| `5677045` | `47ef757` | 2026-11-26T19:26:20-08:00 | Add the all-test seal for external cohorts (--all-test) |
| `384bdf1` | `01d7823` | 2026-11-26T19:26:20-08:00 | Add the CGGA external-partition eligibility script with a dry run |
| `a092a16` | `33e3e9b` | 2026-11-26T19:26:20-08:00 | Add the final OSF Sampling-section text |
| `0c20999` | `4149e4a` | 2026-11-26T21:58:16-08:00 | Add the T2 Brier row and the audit corrections to code, tests and packaging |
| `10cfd93` | `c14c56a` | 2026-11-26T21:58:34-08:00 | Record the T2 Brier outcome and bring the documentation back in line with the code |
| `8f79b33` | `cc2f341` | 2026-11-26T21:58:34-08:00 | Add the project record: overview, findings, work log, commit history and future plan |
| `324a09b` | `65244a7` | 2026-11-26T22:25:33-08:00 | Add read discipline for the held-out cohort: access ledger, read guard, registration state |
| `f3c2a4b` | `96303dd` | 2026-11-26T22:30:16-08:00 | Reconcile the documentation with the measured state and record the stage evidence |
| `b83f9f7` | `c8168a3` | 2026-11-26T23:10:59-08:00 | Add the real-data slice L0 to L3 on TCGA and record what it found |
| `4fd8964` | `90439e5` | 2026-11-27T01:45:19-08:00 | Record the unexplained CGGA access and harden the read guard |
| `45166f0` | `d008d59` | 2026-11-27T02:02:51-08:00 | Impute KPS inside each cross-validation fold and pin the TCGA fits |
| `f94c0d5` | `a5afb04` | 2026-11-27T02:02:51-08:00 | Lock the environment, add a floor CI job and test the partition-drift job |
| `441c94a` | `3980950` | 2026-11-27T02:03:06-08:00 | Correct the T2 account: KPS imputation, the Brier excess and the C-index contrasts |
| `9216e6f` | `6797b98` | 2026-11-27T02:08:29-08:00 | Let the partition-drift test pass in a clean checkout |
| `cff12c2` | `8baf99a` | 2026-11-27T02:13:47-08:00 | Record stage R1: the review conditions closed or deferred |
| `45fdc2d` | `321592f` | 2026-11-27T03:03:22-08:00 | Add single-sample subtype assignment over the published Verhaak 840-gene signature |
| `2e19939` | `359a4ef` | 2026-11-27T03:03:23-08:00 | Give the feature store an honest contract: every patient, provenance, reasons, views |
| `b0927de` | `4b199a5` | 2026-11-27T03:12:36-08:00 | Record stage 3b: Part A, the decisions, the classifier comparison and the build |
| `16b6c09` | `c47067c` | 2026-11-27T03:24:07-08:00 | Correct the T1 scoring sets in pre-registration sections 3.3 and 4.5 |
| `5dadf9b` | `7d0b5cc` | 2026-11-27T03:25:29-08:00 | Align the T1 scoring rule's wording with pre-registration item 9(d) |
| `f5b9fe6` | `c68f556` | 2026-11-27T03:39:51-08:00 | Fix the T1 deconvolution mixture specification before any run |
| `3bf0094` | `315b340` | 2026-11-27T03:47:16-08:00 | Record T1's state-deconvolution row, run before deposit |
| `947cae8` | `4d917fc` | 2026-11-27T04:02:03-08:00 | Record T1's deconvolution result in the pre-registration and update section 5 |
| `1398cc0` | `0541d3b` | 2026-11-27T04:15:56-08:00 | Mark the 3CA confirmation email as optional |
| `99dc69d` | `83f8b04` | 2026-11-27T04:23:42-08:00 | Route the CGGA eligibility step through the read guard and add --seal |
| `e6cf0d6` | `bc598db` | 2026-11-27T04:40:13-08:00 | Seal the external partition cgga_external_v1 |
| `1ba5762` | `b594e2e` | 2026-11-27T04:41:34-08:00 | Record a post-seal checksum access in the ledger |
| `d5140b5` | `15c21ae` | 2026-11-27T09:44:10-08:00 | Add the T2 comparator shortlist and the OSF deposit checklist |
| `aacedb1` | `730559c` | 2026-11-27T09:48:46-08:00 | Register the T2 comparator rule (path C) and name the T1 classifier |
| `0f7e82a` | `508a10f` | 2026-11-27T10:48:35-08:00 | Record the T2 comparator search run under the registered rule |
| `5f0e0c2` | `8250eff` | 2026-11-27T10:56:03-08:00 | Amend the T2 comparator rule to a refitted published gene set (option B) |
| `215a67f` | `d8cdfd4` | 2026-11-27T10:59:30-08:00 | Add the gene-set re-implementation for T2's comparator (option B) |
| `ecc9630` | `9adcf3f` | 2026-11-27T11:01:54-08:00 | Use a batch-stratified C-index for T2's improvement row |
| `5199628` | `d73b838` | 2026-11-27T11:44:05-08:00 | Close the T2 comparator search: one model eligible as first registered |
| `573472c` | `f6bb39f` | 2026-11-27T11:46:25-08:00 | Name T2's comparator: PMID 38174685 as published; option B secondary |
| `bd01685` | `d74b24c` | 2026-11-27T11:47:47-08:00 | Fix option B's candidate gene sets before any fit |
| `23bb56b` | `04e3c96` | 2026-11-27T11:48:47-08:00 | Freeze option B's secondary comparator for T2's improvement row |
| `be9675b` | `77a7fd5` | 2026-11-27T11:49:27-08:00 | Keep T2 comparator artefacts byte-exact; flag manifest line-ending checksums |
| `887394b` | `895d8d8` | 2026-11-27T11:57:15-08:00 | Store the sealed partition manifests byte-exact (checklist 2.5) |
| `71d6a18` | `64ca8bb` | 2026-11-27T11:59:15-08:00 | Record T2's open decisions, drop GLASS, and record the CGGA terms check |
