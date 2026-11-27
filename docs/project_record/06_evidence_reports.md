# 6. Evidence reports

One report per stage, in the order the stages were run. Each report names the commits
it describes; the report itself is committed after them. Numbers were taken from the
command output of the run they describe, on 2026-11-26, and not from memory.

---

## Stage 1. Secure baseline (2026-11-26)

**Scope.** Inventory the repository, verify the work the previous session left
uncommitted, and commit it.

**Changed files.** Commits `0c20999` (26 files: code, tests, packaging, CI comments),
`10cfd93` (9 documentation files) and `8f79b33` (the six files of this record). The
content is the previous session's work, unchanged: every diff was read before it was
committed.

**Commands and results.**

| Command | Result |
| --- | --- |
| `python -m pytest -q`, before any change | 984 passed, 1 skipped in 209 s (985 collected) |
| `ruff check --select E9,F63,F7,F82 src tests scripts` | clean |
| `ruff check --statistics`, run on `src`, `tests` and `scripts` | 171 findings: 31, 133 and 7 |
| `gliostrata freeze-partition --input data/raw/tcga_gbm/clinical_cohort.parquet --manifest data/partitions/tcga_gbm_clinical_v2.json --verify` | VERIFIED. 594 ids reproduce 425 / 169; id-list and source content hashes match |
| `threshold_table_checksum()` | `96bbcd49dcda20f8a0bc9c8b71427aaa24cf754dc6d8d7b4f7d6a2df3d2874ab`, the frozen constant |
| `gliostrata pipeline-dry-run --n 20` | L0 to L7 ran on synthetic data; no file written |
| secret scan of `git diff` and of the untracked files | two hits, both documentation text saying a dbGaP key belongs in `.env`. No credential, no user-specific path |
| commit `0c20999` in a throwaway worktree: `pytest tests/test_workflow.py tests/test_cli.py` | 78 passed, 1 skipped. Worktree removed afterwards. **Corrected in stage R1:** the package under test was not `0c20999`'s. `gliostrata` is installed in editable mode from the main checkout, and pytest resolved it through that install, so the two test modules of `0c20999` ran against the main checkout's `src/`. Verified: those two modules at that commit. Not verified: that commit's package code in isolation. Stage R1 repeats the check with `PYTHONPATH` set to the checkout and `gliostrata.__file__` printed |

**Review findings on the uncommitted work.**

- No threshold, direction, decisive bound, decision rule or partition manifest changed.
  `thresholds.py` and `data/partitions/` were not in the diff.
- The Brier script takes its verdict from the registry and refuses to overwrite a
  result. It has one fallback the registration does not describe: a percentile interval
  that excludes its own point estimate is widened and flagged. It was not triggered:
  `brier18_external.json` has no `ci_widened_to_contain_point` key.
- Two stale numbers: the Makefile said about 148 lint findings, and
  `docs/reproducibility.md` said 148 and "965 tests". Corrected in stage 2.
- The README had one broken sentence (about optional dependencies), and a CI comment
  said tests read CGGA. No test does. Both corrected in stage 2.

**The pre-registration had changed since `a092a16`, and this report did not say so.**
The independent review of `b83f9f7` found the omission; recorded here on 2026-11-27.
Since `a092a16` the document changed in commit `10cfd93`: 7 hunks, 42 lines added and
12 removed. Three of them alter registered procedure: hunk 1 (section 2.2 item 7, the
manifest is sealed after "the T2 evaluation (two reads)" instead of "the single T2
read"), hunk 5 (the `cgga_external_v1` row of section 4.5, the same correction) and
hunk 6 (section 5: CGGA read twice, and "no reads beyond the registered T1, T4 and T5"
in place of "no further reads"). All three are disclosed in the document's own
amendment log, in the row that hunk 7 adds ("read once" corrected to two reads). The
other hunks record the Brier result. No threshold, no decisive bound and no decision
rule changed.

**Test counts.** 985 collected: 984 passed, 1 skipped. The skip is the
Snakemake-optional workflow test.

**Environment.** Python 3.14.3, the package installed in editable mode from this
repository, in a global environment of 389 packages (no virtual environment). Versions
that matter: numpy 2.4.3, pandas 3.0.1, scipy 1.17.1, scikit-learn 1.8.0, pyarrow
24.0.0, click 8.4.1, PyYAML 6.0.3, pytest 9.1.0, ruff 0.15.17. Not installed:
snakemake, dvc, lifelines, scikit-survival, botorch. A raw `pip freeze` contains one
machine-specific line, so it was kept outside the repository.

**Artifacts and checksums.** SHA-256 of every data, result and partition file outside
the download caches, taken on 2026-11-26 with PowerShell `Get-FileHash`. File bytes
were hashed only; no table was deserialised.

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `data/partitions/tcga_gbm_clinical_v1.json` | 13,666 | `f0a363142f46cd94ae4e940f9eee0a68425d558215e19df530c82922de44509d` |
| `data/partitions/tcga_gbm_clinical_v2.json` | 13,605 | `17adb31d6375df07f86975c2354fc91487119d21f9f2ae4c4d906597e8c55c16` |
| `data/processed/tcga_gbm/molecular_calls_joined.parquet` | 12,297 | `5961d17c70b1f0fc727582e7e8b31bb9fe75eaa33533e9d7e64482983e91c77f` |
| `data/processed/tcga_gbm/molecular_calls_manifest.json` | 386 | `b703c47c841f0c31c49e8e3cfb8f26c15b279ddb84776a0e8639b6ee6d537d52` |
| `data/results/cgga/brier18_external.json` | 5,558 | `aca808f0be72b2b7eba2bf25b3ea93e81cb46a2f4308a7f5adf82c9d76e1e2a1` |
| `data/results/cgga/cox_external_validation.json` | 4,934 | `2f2c48a5c17197f1ddf29d50fbd4ee2afa81d0b497ef98877c9571766470595d` |
| `data/results/tcga_gbm/idh_sensitivity_cox.json` | 1,496 | `dcf451c2a8df4bd1596320ee6141ad833d89ec0bf7fc83b9996f7bf5eea50867` |
| `data/results/tcga_gbm/idh_sensitivity_subtype.json` | 2,534 | `98647e35f76558e44308c20d1d40d5fc383ee646b11198645735f3607dc01a81` |
| `data/raw/cgga/clinical_cohort.parquet` | 15,202 | `34bf666674a3e32f9b1d90f397efcc2205fde6ea79112230c890376f210a7940` |
| `data/raw/cgga/clinical_raw.parquet` | 27,374 | `3ce8c0d26251dc6c000d19961cc7ff7abb006fc7b65cc925a5826b55664bd169` |
| `data/raw/cgga/counts.parquet` | 76,913,732 | `2189b07c1eab5e071e1e8d2adeced0bdb2417c7dbdb72207fd1578ecb13df0b4` |
| `data/raw/cgga/gene_annotation.parquet` | 406,302 | `e88a8838f0bd5049c1da1db6710e3affc2e6aaa62c4033825c8da8bc2385276c` |
| `data/raw/cgga/harmonisation.json` | 2,074 | `fcb24b9bdbc3cbe3df3f36cc1cb6da114151c455e7e641e60a19f123317e4fc6` |
| `data/raw/tcga_gbm/clinical_cohort.parquet` | 18,399 | `a7f77418dc90ffbb327773de699965901324be6127cabf2ef745a12a07aa4bda` |
| `data/raw/tcga_gbm/clinical_raw.parquet` | 51,445 | `cca82de43bc966b44f950bb02177a27268057b0ddefc78dac8e18dd3e5105d1e` |
| `data/raw/tcga_gbm/expression/counts.parquet` | 29,808,370 | `bc78e730b43d9a187d86dca2d604dec37af7736d0afcb98394b41d4ba1826d48` |
| `data/raw/tcga_gbm/expression/file_manifest.parquet` | 61,445 | `ce160dd1d59d13299ce8da7d5dba8268e33cbc24c800fd9ac7bd0ac80b435cf1` |
| `data/raw/tcga_gbm/expression/gene_annotation.parquet` | 1,001,064 | `408434126c8e57807d79b0f44ab7b13ccde25699155d8fd0236f9adadba34960` |
| `data/raw/tcga_gbm/expression/sample_metadata.parquet` | 24,315 | `35391e98858deb13cbc3931eb517a98a340ead98c5c3e709744010acb71781d7` |
| `data/raw/tcga_gbm/expression/tpm.parquet` | 70,154,079 | `6f4cdfcf8dc831bc56554e222ff298d81bdc90c02bdaaac707ca94ce07eebe31` |
| `data/raw/tcga_gbm/methylation/file_manifest.parquet` | 91,573 | `f3c1809ccbf581eb449278119df0c9c187a50acab07f3072bd698c82a2194f90` |
| `data/raw/tcga_gbm/methylation/mgmt_call_report.json` | 4,810 | `0a02b2b3899dbdfd687787c3b13cbb763dc6cb221960b40639b7766028da071d` |
| `data/raw/tcga_gbm/methylation/mgmt_calls.parquet` | 46,690 | `be464959a294d9140721e4fe63092e5adc8e217d9ad602e0afa73f329eae909d` |
| `data/raw/tcga_gbm/mutations/call_report.json` | 2,723 | `4f98e21a0c84edc2541f15c7c2689b2fb22556cd1a35df0c3294adbd6a3043d5` |
| `data/raw/tcga_gbm/mutations/file_manifest.parquet` | 90,677 | `a41acbf0ea2d7ce6c62ab6cac05d2ecabf993576252f4ec3bc3a22d8a393d91f` |
| `data/raw/tcga_gbm/mutations/idh_calls.parquet` | 32,046 | `d8a5ea4dd456e1331bea252bc806d16a28c6d9662b0d961337d01f3012159fc4` |

None of these files is in git (`data/raw/**`, `data/processed/**` and
`data/results/**` are ignored; the two manifests are). The table is what a later run
can be compared against.

**External reads.** One checksum-only access: the bytes of four CGGA files were hashed
(ledger entry 6). Nothing was deserialised and no value was inspected. The full test
suite read no CGGA file: every real-data guard in `tests/` points at a TCGA path.

**Unresolved issues.**

- The pre-registration is not deposited.
- `cgga_external_v1` is not sealed.
- No lockfile and no container image: dependencies are floors.
- Result files and the eligibility table are git-ignored; their checksums are recorded
  here instead.
- The GDC is queried live, not pinned to a data release.

**Is the next gate satisfied?** Yes. The working tree is committed, the suite passes,
the development partition verifies, and the artifact checksums are on record. The next
gate is registration and read discipline (stage 2).

---

## Stage 2. Registration and read discipline (2026-11-26)

**Scope.** Make the read rules enforceable and the registration state
machine-readable, close the CI labelling gap, and reconcile the documentation with
the measured state.

**Changed files.** Commit `324a09b`: `data/access_ledger.jsonl`,
`configs/registration.yaml`, `src/gliostrata/evaluation/read_guard.py`,
`src/gliostrata/cli.py`, `tests/test_read_guard.py`, `.github/workflows/ci.yml`,
`Makefile`. The documentation is in commit `f3c2a4b` (the 32nd): `README.md`,
`CHANGELOG.md`, `docs/cgga_validation.md`, `docs/reproducibility.md`,
`docs/roadmap.md` and this record.

**Commands and results.**

| Command | Result |
| --- | --- |
| `pytest tests/test_read_guard.py` | 31 passed |
| CI partition job, its script extracted from `ci.yml` and run locally, `ACTIVE_MANIFESTS` set to the v2 manifest | v1 labelled superseded, v2 ACTIVE; "All partition manifests verified" |
| the same, with a listed manifest that does not exist | the job stops: "ACTIVE_MANIFESTS lists ... which does not exist" |
| the same, with v1 and v2 both listed | both labelled ACTIVE and integrity-checked |
| `gliostrata access-ledger` | 7 entries; "Value reads of the external cohort: 2"; "Registration: NOT deposited" |
| full suite, Python-level socket entry points replaced before pytest | 1015 passed, 1 skipped in 243 s; 0 blocked connection attempts |
| after the documentation edits: `pytest tests/test_workflow.py tests/test_cli.py tests/test_read_guard.py` | 109 passed, 1 skipped |
| commit `324a09b` in a throwaway worktree, the same three modules | 109 passed, 1 skipped. Worktree removed. **Corrected in stage R1:** the same limitation as the stage 1 worktree check. The three test modules of `324a09b` ran against the main checkout's `src/`, not that commit's package code |
| `ruff check --select E9,F63,F7,F82 src tests scripts` | clean |
| `ruff check` on the three new or changed Python files | clean |
| advisory lint, whole project | 171 findings, as before: the new files add none |

**Test counts.** 1016 collected in 26 modules: 1015 passed, 1 skipped.

**Artifacts.** `data/access_ledger.jsonl` (7 entries) and `configs/registration.yaml`
(OSF DOI `null`) are in git from `324a09b`. No data or result file changed; the stage 1
checksums still apply.

**External reads.** None. The stage wrote no script that opens CGGA, and the test
suite reads none.

**What the guard does, and does not do.** A script declares an access before it opens
the data (`record_access`). Checksum, schema and availability accesses are recorded and
allowed. A value read is refused unless the registration is deposited, the analysis is
one of the eleven registered metrics or carries its own pre-registration, the manifest
is sealed or the scoring set is fixed by a registered rule, the analysis has not been
read before, and the SHA-256 of every input is given. It is cooperative: a script that
does not call it is not stopped. It checks that a manifest exists and is sealed, not
that it verifies against its cohort, so `freeze-partition --verify` is still required.
The two T2 scripts predate it and were not edited.

**Tightened in stage R1 (2026-11-27).** The review of `b83f9f7` found that any
non-empty text passed the scoring-rule gate. That path now needs a rule id from the
committed registry `configs/scoring_rules.yaml`, registered for the analysis, naming
a frozen scoring-set manifest whose bytes still hash to the checksum recorded beside
it. No rule has a frozen manifest yet, so the path admits nothing. The guard also
recomputes every input's checksum at read time and refuses a mismatch. Ledger entries
that are observations or notes can be appended (`append_ledger_entry`); they are
marked "not an access" and never counted as reads.

**Prepared, not run: sealing the external partition.**

```bash
python scripts/cgga_external_eligibility.py
gliostrata freeze-partition --input data/processed/cgga/cgga_external_eligible.parquet --manifest data/partitions/cgga_external_v1.json --name cgga_external_v1 --all-test
gliostrata freeze-partition --input data/processed/cgga/cgga_external_eligible.parquet --manifest data/partitions/cgga_external_v1.json --verify
```

The manifest path then goes into both `ACTIVE_MANIFESTS` entries of `ci.yml`, in the
same commit. The first command is the only one that touches CGGA, and it is an
availability access: identifiers, the presence of a survival time and vital status,
and the column names of `counts.parquet`. Its dry run counted 218 eligible patients on
2026-11-26. Under an all-test seal the assignment does not depend on the seed; the CLI
default, 20251126, is the seed of the other two manifests.

It was not run, for three reasons. A seal refuses to be overwritten, so a wrong seal
leaves a superseded manifest behind. The record lists the cohort definition and the
evaluation unit as open for the lead, and both decide what the seal contains. And
`05_future_plan.md` section 5.1.1 records that the T1 scoring set lies mostly outside
the eligibility this partition registers, which may call for a second, separately
named seal.

**Unresolved issues.**

1. The registration is not deposited. Until it is, T1 to T5 are exploratory and the
   guard refuses every value read.
2. The external partition is not sealed (above).
3. T1 cannot run as registered: no classifier matches the registered rule, and the
   microarray expression matrix is not on disk (`05_future_plan.md` section 5.1.1).
4. The decisions listed for the lead in `05_future_plan.md` section 5.4 are all still
   open, including the comparator for T2's improvement row and whether the calibration
   analyses are run.
5. The eligibility table and the result files are git-ignored.
6. No lockfile and no container image.
7. The GDC is queried live, not pinned to a data release.

**Is the next gate satisfied?** In part. Read discipline is in place and enforced for
any new script. Registration is not: it needs the OSF deposit, which is a human action.
Until then the permitted work is TCGA-only or explicitly synthetic.

---

## Stage 3, first deliverable. Real-data slice, L0 to L3 (2026-11-26)

**Scope.** Run layers L0 to L3 on the real development cohort (TCGA), so that the
first three layers are shown to work on real data. TCGA only: no CGGA value, no
network connection, no outcome used, no model fitted.

**Changed files.** One commit, the one after `f3c2a4b`:
`scripts/build_tcga_feature_store.py`, `tests/test_build_tcga_feature_store.py`,
`data/access_ledger.jsonl` (entry 7) and the documentation (`README.md`,
`CHANGELOG.md`, `docs/roadmap.md`, `docs/reproducibility.md` and this record). No
module under `src/` changed: the L2 and L3 code ran as committed.

**What the script does.** It loads the clinical cohort, the gene x sample TPM
matrix, its annotation and the joined MGMT / IDH calls. It calls `classify_verhaak`,
`deconvolve_states_expression` and `deconvolve_tme_expression`, each with its
defaults, on the cohort patients that have an expression sample. It assembles the
store with `build_feature_store` and writes the store and a report: the counts,
every exclusion by reason, the range of each feature, and the SHA-256 of every
input. A path containing "cgga" is refused before anything is opened.

**Commands and results.**

| Command | Result |
| --- | --- |
| `python scripts/build_tcga_feature_store.py` (default policies) | Stops: "64 patient(s) have no kps". Nothing is written |
| `python scripts/build_tcga_feature_store.py --missing-demographic drop` | 155 patients; store and report written. Repeated after the last code change: same store, same SHA-256 |
| `pytest tests/test_build_tcga_feature_store.py` | 26 passed |
| full suite, sockets blocked | 1041 passed, 1 skipped in 225 s; 0 blocked connection attempts |
| after the last code change: the new module with `test_cli`, `test_workflow` and `test_features` | 112 passed, 1 skipped |
| `ruff check --select E9,F63,F7,F82` on the two new files | clean |
| `ruff check` on the two new files | clean. The project's advisory count stays at 171 |
| the 26 checksums of the stage 1 table, recomputed from the files | 26 of 26 match |

**Test counts.** 1042 collected in 27 modules: 1041 passed, 1 skipped.

**Result on the real cohort.**

| Step | Patients |
| --- | --- |
| Cohort (L0) | 594 |
| With an expression sample. This is the scoring set | 283 |
| Of those, with a KPS | 219 |
| Of those, with both a MGMT and an IDH call. This is the store | 155 |

Of the 283, 54 have no MGMT call, 36 no IDH call and 78 at least one of the two.
None lacks age or sex. In the store, 46 are MGMT-methylated and 11 are IDH-mutant.
No patient lacks a vital status; 3 lack a treatment arm. Subtype labels, which are
a within-cohort ranking: Classical 46, Mesenchymal 44, Proneural 40, Neural 25. The
median of each patient's highest subtype probability is 0.48. Age runs from 21 to
89 (median 61) and KPS from 40 to 100 (median 80). State and TME values lie in
[0, 1] and sum to 1. The two resection features are 0 for everyone.

**Artifacts of this stage.** Both are git-ignored.

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `data/processed/tcga_gbm/feature_store.parquet` | 34,492 | `6cb446e8f3ce50672c7e39f5924fa1746a5740f49edafdbc49f2e12fbf1e219c` |
| `data/results/tcga_gbm/feature_store_report.json` | 5,342 | `a9d23c650748cbb3b551344f862010f3051b47ba6d07195ffc0a12027ec188fb` |

The 26 files of the stage 1 table are unchanged.

**External reads.** No CGGA value was read. The checksum check hashed the bytes of
five files under `data/raw/cgga` (ledger entry 7, dated 2026-11-27 in UTC):
clinical_cohort, clinical_raw, counts, gene_annotation and harmonisation.json.
Nothing was deserialised. The ledger still counts 2 value reads, both for T2.

**Correction to stage 1.** Entry 6 and the stage 1 report say four CGGA files were
hashed. The hash list of that day holds five: `harmonisation.json` was hashed as
well and is in the stage 1 table. It holds ingest facts and no patient value. The
ledger is append-only and entry 6 was not edited; this note is the correction.

**Findings.** Four facts of the real data that the synthetic data never showed.
They are decisions for the lead (`05_future_plan.md` section 5.4), not things this
stage decided.

1. Resection extent is absent for all 594 patients in the TCGA clinical table. The
   L3 code encodes an absent value as 0 in both one-hot columns, which reads as
   "biopsy" for everyone. The script keeps the patients and reports the two features
   as constant (`not_measured_in_source`, `constant_features`).
2. KPS is missing for 152 of the 594 patients, 64 of them among the 283 with
   expression. The default policy stops; the run above used `drop`. Which policy the
   model should use is open.
3. The harmonised cohort carries all-null placeholder columns under the names of
   the L2 outputs and of the molecular calls. `build_feature_store` keeps the
   cohort's column when a name collides, so on this cohort it would return null
   state, TME and subtype values. The script removes all-null placeholders and
   refuses a cohort column that carries values under such a name. L3 itself is
   unchanged.
4. Subtype probabilities depend on the scoring set. A test in the new module shows
   that the same patient gets different probabilities when other patients leave the
   set. The scoring set is therefore recorded: the 283 cohort patients with
   expression. One expression sample lies outside the cohort and is not scored.

**What this is not.** The values are measured, not validated: the subtype label is
a within-cohort ranking, the state and TME values are normalised marker-enrichment
scores, and nothing here addresses T1. L4 to L7 were not run on real data. The CLI
and the Snakemake workflow still run the stub L2; the slice is a script and is not
wired into them.

**Unresolved issues.** Those of stage 2, and:

1. The five decisions added to `05_future_plan.md` section 5.4: resection extent,
   missing KPS, absent molecular calls, the scoring set, and whether L3 should
   itself refuse a shadowing cohort column.
2. The slice is not wired into the CLI, the workflow or DVC.
3. The store covers 155 of 594 patients. Whether that is a usable development
   cohort for L5 is not assessed here.

**Is the next gate satisfied?** In part. The feature store has real values, its
missingness is documented and its feature ranges are reported. The pipeline does
not run end to end on real data: its own L2 is still the stub, and L4 to L7 have
not been run on real data. Registration is not deposited either; that needs the
OSF deposit.

## Stage R1. Conditions of the independent review of `b83f9f7` (2026-11-27)

**Scope.** The review returned "accept with conditions": 0 blockers, 6 major findings
(M1 to M6) and 7 minor ones. This stage closes them or defers them with a reason. TCGA
only: no CGGA value was read, downloaded or hashed, no `validate_*` script was run, and
no OSF action was taken. No threshold, bound, partition, seed or analysis budget changed.
The only model fitting was the fold-wise imputation rerun on TCGA.

**Commits.** `4fd8964` (the access, the guard), `45166f0` (fold-wise imputation, the
fit fingerprint, outcome wording), `f94c0d5` (lock file, floor job, partition-drift
script, provenance checksums), `441c94a` (the T2 account), `9216e6f` (a test that
failed in a clean checkout, found by R6), and the commit that adds this section, tagged
`baseline-r1`. That last commit changes only `docs/project_record/`.

### R1. The unexplained access to the held-out cohort (M5)

**(a) The files.** From file-system metadata only (name, size, times; NTFS last-access
updates are enabled on this machine, `DisableLastAccess = 2`):

| File under `data/raw/cgga/_cache/extracted/` | Bytes | Last access (UTC) | Kind |
| --- | --- | --- | --- |
| `CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt` | 47,679,555 | 05:48:55.811 | expression matrix |
| `CGGA.mRNAseq_325_clinical.20200506.txt` | 25,090 | 05:48:55.906 | clinical |
| `CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt` | 103,559,604 | 05:48:56.050 | expression matrix |
| `CGGA.mRNAseq_693_clinical.20200506.txt` | 52,285 | 05:48:56.485 | clinical |
| `CGGA.mRNA_array_301_clinical.20200506.txt` | 25,020 | 05:48:56.485 | clinical |

The brief named six files. The directory holds five, and no sixth file under
`data/raw/cgga` carries an access time in that window: the five `.zip` archives beside
them were last accessed on 2026-11-26 at their download, and the five top-level files
(four Parquet files and `harmonisation.json`) at 2026-11-27T07:01:26Z, consistent with
the checksum access of ledger entry 7 (recorded at 07:03:46Z). The only other entry in the tree is the `extracted`
directory itself; its access time has since been overwritten by later listings and
cannot show whether it was the sixth item. Re-checked from metadata during the closure of
this stage: still five files, the same sizes and access times. No file is an annotation
table.

**(b) Ledger.** Entry 8 (`kind: observation`, recorded 2026-11-27T08:35:15Z): the paths,
sizes and times above, "process not identified at the time of the observation", no value
use known. Closed by a later entry, not by an edit.

**(c) Investigation, from logs only.** Entry 9 (`kind: note`) names the process: the
review agent's own workspace-wide text-search tool, called at 05:48:54Z with two
path-regex queries at the start of the stage 1 review. It timed out after 30 s
("Search timed out after 30000ms", 05:49:25Z) and returned nothing. Evidence: the
agent's session log (tool call at 05:48:45Z, result at 05:49:25Z) and the access times of
491 repository files that fall in the same path-ordered sweep from 05:48:55Z to
05:49:25Z. The five files were scanned for two path strings; nothing was deserialised,
interpreted or returned, and the test suite did not read them. Re-checked independently
at the closure of this stage, from the agent's session record of that review (session
opened 05:44:26Z): the assistant message carrying the `search_codebase` call is stamped
05:48:45.604Z, and its tool result, stamped 05:49:25.489Z, holds the two queries
(`cgga[^"']{0,60}\.parquet` and `raw/cgga|processed/cgga|raw\\cgga`), each with an empty
result and "Search failed: Search timed out after 30000ms". Not identifiable from the
logs: which of the five files were fully read and which only partly before the timeout.

**(d) Recurrence.** `.clineignore` and `.rgignore` exclude both CGGA directories;
`.clinerules` forbids workspace-wide searches. Deferred: removing the plain-text copies
from the workspace (`05_future_plan.md`, section 5.5, with the reason); and a live check
that the agent's search tool honours the ignore files, not run because a failing check
would open the files again.

### R2 to R5. The T2 account and the fold-wise imputation (M1, M2, minor)

| Item | Where | Result |
| --- | --- | --- |
| R2 Brier explanation | `docs/cgga_validation.md` 12.2; `02_findings.md` | With no resolution, the excess over the null model is the squared offset plus the variance of the predictions; the two cannot be told apart without calibration or resolution analysis, and neither direction was examined |
| R3 C-index contrasts | `02_findings.md`; `cgga_validation.md`; pre-registration 2.3 | Descriptive, no interval. Variant 1 against 2 bundles the KPS value, the missing-KPS indicator (HR 1.57) and the fit set (592 against 570 modelled) |
| R4 KPS disclosure | every table row and paragraph that showed variant 1's fit as 594 | "594 passed, 592 modelled; KPS set to the median, 80, for 150 with an indicator (HR 1.57); treatment arm missing for 22 (level, HR 1.20)". "Nothing imputed" is scoped to variants 2 and 3. Pinned by `tests/fixtures/tcga_fit_fingerprint.json` and `tests/test_tcga_fit_fingerprint.py`. `docs/osf_final_fields.md` carries no variant-1 fit n |
| R5 fold-wise imputation | `prediction/cox.py` (`fit_design_plan`, `apply_design_plan`) | The KPS median now comes from each fold's training rows |

**R5 rerun, TCGA only** (`python scripts/rerun_tcga_internal_cv.py`; 5 folds, seed 42,
1000 bootstrap resamples; record `data/results/tcga_gbm/foldwise_imputation_rerun.json`,
git-ignored):

| Analysis | Patients (modelled) | KPS imputed | C-index, whole-frame median (old) | C-index, training-fold median (new) | 95% CI |
| --- | --- | --- | --- | --- | --- |
| Variant 1, full cohort | 594 (592) | 150 | 0.673544 | 0.673544 | 0.647818 to 0.699370 |
| Variant 1, without the 23 IDH-mutant patients | 571 (569) | 143 | 0.665650 | 0.665650 | 0.637019 to 0.691274 |

Every training fold's median KPS is 80, the same as the whole-frame median, so the
design matrices are byte-identical and the C-index does not move. Both values are kept.
No external number depends on this path: variant 1 is never evaluated on CGGA, and
variants 2 and 3 use `impute="drop"`.

### R6. Isolated clean-checkout check (M4)

The two earlier worktree checks (stage 1, commit `0c20999`; stage 2, commit `324a09b`)
ran the test modules of those commits against the main checkout's `src/`, because
`gliostrata` is installed in editable mode from the main checkout. Both rows above are
reworded to what was verified.

Repeated properly: a fresh `git clone` of this repository into a scratch directory,
run with `PYTHONPATH=<clone>/src` and the same socket block and CGGA audit hook as R9.
The data directories are git-ignored, so the clone has no TCGA data and every
real-data test skips; that is what a clean checkout and CI see.

**First run, at `441c94a`: 1 failed**, 1051 passed, 23 skipped.
`test_the_repository_manifests_pass_with_the_ci_setting` (in
`tests/test_ci_partition_drift.py`, added in this stage) asserted that the TCGA cohort is
on disk. It is git-ignored, so the test would have failed in every clean checkout and in
CI. `gliostrata` resolved inside the clone. Fixed in `9216e6f`: the test now asserts the
line the script prints in each case.

**Second run, at `9216e6f`, in a new clone** (the last commit that changes code; the
tagged commit after it changes only this record):

| Check | Result |
| --- | --- |
| `python -c "import gliostrata; print(gliostrata.__file__)"` | `<clone>\src\gliostrata\__init__.py` |
| `gliostrata.__file__` read inside the pytest process after the run | `<clone>\src\gliostrata\__init__.py` |
| every loaded `gliostrata.*` module whose file lies outside the clone | none |
| full suite in the clone, sockets blocked | 1052 passed, 23 skipped in 217 s; 0 blocked connection attempts; 0 accesses to a CGGA path inside the clone |

The 23 skips are 22 real-data tests whose inputs are git-ignored and therefore absent from a clone, and the Snakemake-optional workflow test. 1052 + 23 = 1075, the same collection as the main checkout. What this verifies: the package code and tests of `9216e6f` on their own, with no data. It does not verify the real-data tests in isolation; those ran on the main checkout (R9).

### R7. Closure (M3)

| Item | Status |
| --- | --- |
| Lock file | **Closed.** `requirements-lock.txt`: 37 distributions, Python 3.14.3 (Windows ARM64), the runtime and `dev` closure. CI's `test` job installs it, then the package with `--no-deps`. The `workflow` extra (`snakemake`, `dvc`) is not locked: it is not in the environment that produced the results |
| Floor job | **Added, not yet run.** `test-floor`: Python 3.10, `requirements-floor.txt`. It has not run on GitHub (no network in R1) and may be red |
| Provenance manifests | **Closed by checksum.** The three per-modality GDC file manifests carry patient barcodes, so they are patient-level linkage and are not committed. `data/provenance/file_manifest_checksums.json` records their bytes, rows and SHA-256 (expression 284 rows, methylation 450, mutations 464); `tests/test_provenance_manifests.py` checks them when the data are present |
| Calibration | **Decided: not done for this thesis.** Not a registered row; a further read of CGGA would not change the T2 verdict; the offset-versus-dispersion question is named as a limitation (`cgga_validation.md` 12.5) |
| CI partition-drift job | **Closed for its first step.** `scripts/ci_partition_drift.py`; `tests/test_ci_partition_drift.py` executes it against a missing, a superseded, an external all-test and a tampered manifest. **Deferred:** the second step (the shell loop around `freeze-partition --verify`) is still inline and untested |
| Record | **Closed.** Stage 1 now cites `10cfd93` and states that the pre-registration changed since `a092a16` (hunks 1, 5 and 6 alter registered procedure, as disclosed amendments). Entries 7 onward carry a UTC time of day; entry 10 marks entries 0 to 6, value reads 1 and 2 among them, "time of day not recorded" |
| Wording | **Closed.** "Reads no outcome" is now "uses no outcome value"; `test_outcome_values_are_carried_and_counted_but_never_used` tests `n_missing_outcome` |

### R8. Guard hardening (minor)

**Closed** in `4fd8964`. A scoring-rule read now needs a rule id from
`configs/scoring_rules.yaml`, registered for the analysis, naming a frozen scoring-set
manifest whose bytes still hash to the recorded checksum; input checksums are recomputed
at read time. Tests: unknown rule refused, rule without a frozen manifest refused,
altered manifest refused, free text refused, checksum mismatch refused. No rule has a
frozen manifest yet, so that path admits nothing.

### R9. Finish

| Command | Result |
| --- | --- |
| full suite on the main checkout at `9216e6f`, Python-level socket entry points replaced and a file-access audit hook on every CGGA path | 1074 passed, 1 skipped in 261 s (and 252 s on the working tree before the commits); 0 blocked connection attempts; 0 opens, listings or stats of a path under the repository's `data/raw/cgga` or `data/processed/cgga` (126 events, all on `tmp_path` fixtures that imitate the layout) |
| `ruff check --select E9,F63,F7,F82 src tests scripts` (blocking) | All checks passed |
| `ruff check src tests scripts` (advisory) | **170 findings**: N803 102, N806 55, N802 7, UP037 3, N812 2, UP035 1. 171 before stage R1 |
| `git tag -a baseline-r1` | on the commit that adds this section |

**Test counts.** 1075 collected in 31 modules: 1074 passed, 1 skipped (the
Snakemake-optional workflow test).

**External reads.** None. The R1 metadata listings read names, sizes and times only;
the ledger still counts 2 value reads, both for T2.

**Left open.** The extracted plain-text CGGA copies (deferred); the live check of the
ignore files (deferred); the `test` and `test-floor` CI jobs have not run on GitHub; the
second partition-drift step is untested; the comparator model for T2's improvement row
(`05_future_plan.md` 5.4); the OSF deposit and the seal of `cgga_external_v1`.
