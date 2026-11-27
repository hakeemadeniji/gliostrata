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
| commit `0c20999` in a throwaway worktree: `pytest tests/test_workflow.py tests/test_cli.py` | 78 passed, 1 skipped. Worktree removed afterwards |

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
`Makefile`. The documentation is in the commit after this report: `README.md`,
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
| commit `324a09b` in a throwaway worktree, the same three modules | 109 passed, 1 skipped. Worktree removed |
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
