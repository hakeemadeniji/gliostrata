# Reproducibility

Point 10 of the GlioStrata defensibility stack is a full reproducibility
release: code, containers, manifests, and a Zenodo DOI. This document is the
honest status report on that goal. It is deliberately written so that the
"not yet reproducible" section is as prominent as the instructions.

**Summary in one line:** the *plumbing* is reproducible end to end today; the
*science* is not, because the pipeline's L1, L2, L5 and L6 steps are stubs. Real
implementations of L2 (subtype and state/TME scoring), L5 (a fitted Cox model) and
L6 (a GP optimizer) exist as modules (§6). Since stage 3b (2026-11-27) the real L2 and
L3 run on TCGA through one CLI command, `gliostrata build-feature-store --mode real`;
the stub `gliostrata profile` and L5 to L7 do not consume its output.

---

## 1. Reproducing from a clean checkout

### 1.1 The offline path (fully reproducible, ~30 seconds, no network)

This is what CI runs and what a reviewer should run first. It exercises the
real L0→L7 code path on a seeded synthetic fixture.

```bash
git clone <repo> && cd gliostrata
python -m venv .venv && . .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev,workflow]"

make test              # the full pytest suite
make dag               # print the real-cohort Snakemake DAG; executes nothing
make pipeline-fixture  # run the whole pipeline on 12 synthetic patients -> .test/
```

To reproduce the environment that produced the results on record, install from the
lock file instead of from the floors:

```bash
pip install -c requirements-lock.txt -e ".[dev,workflow]"
```

`requirements-lock.txt` (2026-11-27) pins the 37 installed distributions of the
project's runtime dependencies and its `dev` extra, as installed with Python 3.14.3.
It is used as a constraints file here so that `snakemake` and `dvc` (the `workflow`
extra), which are not in it, still resolve. CI's `test` job installs it exactly
(`pip install -r requirements-lock.txt`, then the package with `--no-deps`).

`make pipeline-fixture` writes, under `.test/`:

| Layer | Artifact |
| --- | --- |
| L0 | `.test/interim/fixture_synthetic_cohort.parquet` |
| L1 | `.test/interim/fixture_harmonized.parquet` |
| L2 | `.test/processed/fixture_profiled.parquet` |
| L3 | `.test/processed/fixture_feature_store.parquet`, with `_labels.parquet`, `_missing.parquet` and `_columns.json` beside it (all labelled synthetic), and `.test/results/fixture/feature_store_build_report.json` |
| L6 | `.test/results/fixture/configs/ranked_configurations.parquet` |
| L7 | `.test/results/fixture/protocols/protocol_*.json` (12 files) |
| — | `.test/results/fixture/pipeline_complete.txt` |

It is bit-reproducible given the same seed: `gliostrata synthetic-data` is a
pure function of `(n, seed)`, and every downstream layer takes an explicit
`--seed`. Change `n_synthetic=` or `seed=` on the command line to vary it.

There is also a simpler in-memory smoke test with no Snakemake and no files:

```bash
make pipeline     # gliostrata pipeline-dry-run --config configs/default.yaml
```

### 1.2 The real-cohort path (network, ~1.3 GB, hours)

```bash
make ingest           # L0: GDC clinical + bulk RNA-seq. NETWORK. ~1.2 GB cached.
make verify-partition # MANDATORY gate before anything downstream (see §4)
make features-real    # L0 -> L3 on the real cohort
```

Alternatively, `snakemake --cores 4 --snakefile workflow/Snakefile`, which
resolves the complete DAG and runs only what is missing or stale.

This path is **not** reproducible in the strict sense, for two reasons that are
properties of the upstream source, not of this code:

1. The GDC REST API is a live service. `scripts/download_tcga_clinical.py`
   issues a query against today's GDC, not a pinned snapshot. The TCGA-GBM
   cohort has already changed under this repo: the clinical cohort grew from
   583 to 594 rows during a single day of development (see §4).
2. There is no pinned GENCODE/GDC data-release version in the ingest call, so
   a future re-run can pick up a different gene annotation.

Fixing (1) and (2) means recording the GDC data release in the output manifest
and filtering on it. That is not done yet and is a prerequisite for the Zenodo
release.

### 1.3 Pointing the workflow anywhere

Every path is derived from `config`, so nothing is hardcoded:

```bash
snakemake --cores 2 --snakefile workflow/Snakefile \
  --config source=synthetic cohort=myfixture n_synthetic=40 seed=7 \
           with_expression=false verify_partition=false \
           raw_root=.test/raw interim_root=.test/interim \
           processed_root=.test/processed results_root=.test/results
```

Keys and defaults are defined in one place, `workflow/rules/common.smk`. The
important one is `source`:

- `source=tcga` (default) — real GDC ingestion. **The only rules in the entire
  workflow that touch the network are `ingest_tcga_clinical` and
  `ingest_tcga_expression`.**
- `source=synthetic` — offline, seeded fixture. Guaranteed to open no sockets.

Keep `*_root` values **repo-relative**. Snakemake's metadata store base64-encodes
absolute output paths into filenames, and long absolute out-of-tree paths
exceed Windows' path limit and abort the run after the jobs have already
succeeded.

---

## 2. What is versioned where

| Thing | Where | Why |
| --- | --- | --- |
| `src/`, `scripts/`, `tests/`, `configs/` | **git** | Code. |
| `workflow/`, `dvc.yaml`, `Makefile`, `.github/` | **git** | Pipeline definition. |
| `data/partitions/*.json` | **git** | Sealed train/test splits. Kilobytes, must be diff-reviewable, must never be regenerated. See §4. |
| `data/access_ledger.jsonl`, `configs/registration.yaml` | **git** | Every access to the held-out cohort, and the machine-readable registration state. Must be reviewable in a diff. |
| `requirements-lock.txt`, `requirements-floor.txt` | **git** | The environment that produced the results on record, and the lowest versions pyproject accepts. CI's two test jobs install from them. |
| `configs/scoring_rules.yaml` | **git** | The registered rules that may fix a scoring set for a read of the held-out cohort, with the manifest checksum each is tied to. Read by the read guard. |
| `data/provenance/file_manifest_checksums.json` | **git** | Sizes, row counts and SHA-256 of the per-modality GDC file manifests, which are not committed because they carry patient barcodes. |
| `.clineignore`, `.rgignore`, `.clinerules` | **git** | Keep workspace-wide search and indexing tools away from `data/raw/cgga/` and `data/processed/cgga/`. See section 7. |
| `data/raw/**` (`clinical_*.parquet`, `expression/*.parquet`) | **DVC candidate — not yet tracked** | Reconstructible from the GDC, but only against *today's* GDC. See §3. |
| `data/raw/**/_cache/` (1.2 GB of STAR TSVs) | **neither** | Pure download scratch. Excluded in both `.gitignore` and `.dvcignore`. Never track it. |
| `data/interim/**`, `data/processed/**`, `data/results/**` | **neither (DVC candidate)** | Derived. Cheap to recompute from `data/raw/`. |
| `models/**` | **neither (DVC candidate)** | No trained models exist yet. |
| `.test/`, `.snakemake/` | **neither** | Scratch. |

Anything in the "neither" rows is, right now, **on one laptop only.** That is
the single biggest gap between this repo and a reproducibility release.

---

## 3. DVC

### Current state

`dvc.yaml` is a correct, parsing description of the nine real pipeline stages.
It was verified with `dvc stage list` and `dvc dag` in an isolated DVC 3.67.1
repo.

`.dvc/` has **not** been initialised, and there is **no remote configured.**
This was a deliberate decision:

- DVC lives in the optional `workflow` extra and is not installed in the
  development environment, so `dvc init` could not be run and its output could
  not be validated here.
- `dvc init` writes a `.dvc/` directory *and installs git hooks*
  (`pre-commit`, `post-checkout`, `pre-push`). Those hooks change what every
  `git commit` and `git checkout` in the repo does. That is a change the
  project lead should make knowingly, not one that should appear in a commit
  about CI.
- Configuring a remote that cannot be tested here would be worse than
  configuring none: a `dvc push` against an untested remote fails after the
  data has already been staged into the local cache, doubling disk use.

### To initialise it

```bash
pip install -e ".[workflow]"
dvc init
git status            # inspect .dvc/ and the installed hooks before committing
dvc dag               # confirm the stage graph
dvc repro             # NETWORK + ~1.3 GB
```

`.dvcignore` is already in place and is sane for either choice.

### Remote storage recommendation

Track the **derived matrices**, not the download scratch:

| Artifact | Size | Track? |
| --- | --- | --- |
| `expression/tpm.parquet` | 67 MB | **yes** |
| `expression/counts.parquet` | 29 MB | **yes** |
| `expression/{gene_annotation,sample_metadata,file_manifest}.parquet` | < 5 MB | **yes** |
| `clinical_{raw,cohort}.parquet` | < 1 MB | **yes** |
| `expression/_cache/*.tsv` | **1.2 GB** | **never** |

**Storage implication, stated plainly.** DVC is content-addressed and
copy-on-write, not delta-compressed. Every re-ingest that changes a byte of
`tpm.parquet` stores a *complete new 67 MB copy*, in the local `.dvc/cache`
**and** on the remote. The ~100 MB of tracked matrices therefore costs roughly
`100 MB × (number of distinct ingests you keep)`. Ten re-ingests over the
project's life is ~1 GB locally and ~1 GB remote — comfortable. The same
arithmetic applied to `_cache/` would be 12 GB for the same ten re-ingests,
for data that `file_manifest.parquet` already lets you re-download
deterministically. Hence `.dvcignore`.

Concretely, for the Zenodo release: an S3/GCS bucket or a university
institutional store, configured as

```bash
dvc remote add -d origin s3://<bucket>/gliostrata
dvc push
```

Do not use Google Drive as the DVC remote. Its API rate limits make
`dvc pull` of a 100 MB multi-file dataset unreliable, which is precisely the
operation a reviewer will perform.

### Snakemake vs DVC

They are not redundant and they do not fight:

- **Snakemake** executes. It owns parallelism, the `--config` fixture
  indirection, logging, and the dry run.
- **DVC** versions. It owns content hashing, remote storage, and `dvc.lock` as
  the provenance record.

`tests/test_workflow.py` enforces that the two definitions agree, so they
cannot silently drift apart again (which is exactly how the previous
`dvc.yaml` became fiction — see §6).

---

## 4. Frozen partitions — do not regenerate

`data/partitions/*.json` holds sealed train/test splits. They are in **git**,
on purpose, and nothing in `workflow/` or `dvc.yaml` is capable of writing one.

### What they are

A split assigned by `sha256(seed + patient_id)` mapped into `[0, 1)`. Because
it is a pure function of `(seed, patient_id)`, it is stable under row
reordering, cohort growth, and library upgrades — unlike `train_test_split`
with a shuffle, which silently re-assigns *every* patient whenever upstream
filtering changes the row count.

Each manifest records the seed, the test fraction, the full train/test id
lists, SHA-256 digests over those lists, and a content hash of the cohort file
it was frozen against. `sealed: true` makes
`gliostrata freeze-partition` refuse to overwrite it.

### Why they must not be regenerated

Once a model has trained on a partition, moving one patient from train to test
destroys that patient's value as an external test case **permanently**. You
cannot un-see data. A re-freeze that "fixes" a drift warning is not a bug fix;
it is silent test-set contamination, and it is undetectable in the final
numbers.

So: when the cohort changes, **freeze a new, separately named partition.** Do
not overwrite the old one. That is already what happened here — see below.

### Verifying

```bash
make verify-partition    # or: gliostrata freeze-partition --input ... --verify
```

`--verify` writes nothing and exits non-zero on drift. Run it before every
evaluation. The CI job `partition-drift` runs the parts of this check that are
possible from git alone (§5). Those parts are `scripts/ci_partition_drift.py`, which
`tests/test_ci_partition_drift.py` executes against a missing, a superseded, an
external all-test and a tampered manifest.

### Known live drift (as of 2026-11-26)

`tcga_gbm_clinical_v1.json` (583 patients) **does not verify** against the
current `data/raw/tcga_gbm/clinical_cohort.parquet`:

```
Source file content hash changed
11 patient id(s) present in the cohort but not in the frozen partition
Row count changed: manifest recorded 583, cohort now has 594
```

This is the system working correctly. The clinical ingest filter was revised
and the cohort grew by 11 patients. A `tcga_gbm_clinical_v2.json` (594
patients) was correctly frozen as a *new* manifest rather than by overwriting
v1.

**Decision (recorded):** v2 is the current partition and the one
`make verify-partition` checks. v1 is kept, unmodified, as the historical record
of the 583-row cohort: CI still integrity-checks it, never verifies it against
the live cohort, and it should not be used for evaluation. Do not delete or
re-freeze it.

---

## 5. Continuous integration

`.github/workflows/ci.yml`, on every push and pull request, Python **3.14**
(the version this project is developed against; `pyproject.toml` declares
`requires-python = ">=3.10"`), except the `test-floor` job, which uses 3.10.

| Job | What it does |
| --- | --- |
| `lint` | ruff at the version in the lock file. `ruff check --select E9,F63,F7,F82` as a **blocking** gate, plus the full project config as **advisory**. See the note below. |
| `test` | Installs `requirements-lock.txt`, then the package with `--no-deps`; `python -m pytest -q`. Since 2026-11-27. |
| `test-floor` | Python **3.10**: installs `requirements-floor.txt` (the lowest versions pyproject accepts), then the package; `python -m pytest -q`. New on 2026-11-27; it tests the `requires-python` claim and has not yet run on GitHub. |
| `partition-drift` | `python scripts/ci_partition_drift.py`: partition manifest integrity, namely sealed, digests match, and every id re-derives its frozen split from `(seed, patient_id)`. Plus a full `--verify` when the cohort happens to be in the checkout. |
| `workflow-dry-run` | `snakemake -n` on the real-cohort DAG, then a real offline fixture run, then `dvc stage list` / `dvc dag`. |

### CI never touches the network

Verified, not assumed. The suite was re-run with `socket.socket.connect`,
`connect_ex`, `create_connection` and `getaddrinfo` all patched to raise. The
results were **identical** to the normal run — no test attempted a connection.
Re-verified 2026-11-26 on the full suite (1042 tests: 1041 passed, 1 skipped; the
Python-level socket entry points replaced before pytest started) with zero blocked
connection attempts, and again on 2026-11-27 (1075 tests: 1074 passed, 1 skipped,
355 s) with zero blocked attempts and with a file-access tripwire that recorded no
open and no listing of a CGGA path.
Specifically:

- `tests/test_tcga_download.py` and `tests/test_tcga_expression.py` load
  `scripts/download_tcga_*.py` by path and drive the pure functions with
  synthetic fixtures. They assert on the GDC URL *constants* without ever
  requesting them.
- The real-data tests (`test_cox.py`, `test_verhaak.py`, `test_deconvolve.py`,
  `test_join_molecular_calls.py`, `test_idh_sensitivity_cox.py`,
  `test_idh_sensitivity_subtype.py`, `test_validate_cox_cgga.py`) are guarded by
  `skipif` / `skip` conditions on files under `data/raw/` and `data/processed/`.
  Both are gitignored, so they skip in CI by construction.
- The `workflow-dry-run` job uses `source=synthetic`, which routes L0 through
  `gliostrata synthetic-data` and never reaches the two ingest rules.

**No test requires network access.**

### Why the blocking lint set is narrow

The repo currently has 170 ruff findings under its own config (30 in `src`,
133 in `tests`, 7 in `scripts`, measured 2026-11-27; it had 171 the day before),
dominated by the `N803`/`N806` naming rules
(argument and variable names) — which is arguable in statistical code
where `X`, `X_train` and `N` are the conventional names. Making the full set
blocking today would land a red CI on the first commit and train everyone to
ignore it. The blocking set is the "the code is broken" floor (syntax errors,
invalid comparisons, misplaced statements, undefined names) and passes clean.
The advisory step prints the backlog. Widen the blocking set as it is cleared.

---

## 6. What is NOT reproducible

Everything below is **plumbing that works, wrapped around science that is not
yet in the pipeline.** The DAG runs green.

### Stubbed layers

| Layer | Module | What it actually does | What the architecture requires |
| --- | --- | --- | --- |
| **L1** harmonize | `gliostrata.data.harmonize` | Schema normalisation and basic QC row filtering. | Cross-cohort batch correction (ComBat/limma), platform harmonization, imputation of missing molecular covariates. |
| **L2** profile | `gliostrata.profiling.{subtype,states,tme}` | Derives "subtypes", Neftel state fractions and TME composition **from the clinical table**. Reads no expression matrix. | Verhaak classification, Neftel state deconvolution and TME deconvolution from `counts.parquet`/`tpm.parquet`. |
| **L5** survival | `gliostrata.prediction.survival` | Produces survival numbers without a fitted model. This is what the CLI uses. | A validated survival model with calibration and a reported C-index. |
| **L6** optimize | `gliostrata.optimization.bayesian` | Ranked enumeration of configurations. | Bayesian optimization against a validated response surface. |
| **L7** handoff | `gliostrata.handoff.protocol` | Structurally complete JSON documents templated from the above. | Bench-executable protocols. |

Implemented, and called only through `gliostrata build-feature-store --mode real`
(TCGA): `profiling/centroid.py` (single-sample subtype, the per-patient feature) and
`profiling/deconvolve.py` (state and TME scoring) read the expression matrix.
`profiling/verhaak.py` (cohort-relative subtype scoring) is used for descriptive
distributions only. Not called by the pipeline: `prediction/cox.py`, a fitted Cox
model, and `optimization/bayes_opt.py`, a Gaussian-process optimizer, reached through
standalone scripts. `gliostrata profile` is still the stub, and `optimize` and
`generate-protocol` refuse the real store, because they would encode an absent value
as a value. See `docs/roadmap.md`.

**The single most important consequence:** on the synthetic path, every subtype
call, cell-state fraction and TME fraction is computed from clinical covariates by
the stubs; they are not measurements. On TCGA, `build-feature-store --mode real`
reads `tpm.parquet` and its annotation, and the values it writes are measured, not
validated: the subtype is a nearest-centroid assignment, and the state and TME values
are cohort-relative marker-enrichment scores, not cell proportions.

This is why `workflow/rules/profile.smk` (the stub) declares `HARMONIZED` as its
only input, while `build_feature_store` with `source=tcga` takes `tpm.parquet` and
its annotation as inputs and writes their SHA-256 into its report. The DAG is a map of
the code, not of the intent. Neither Snakemake rule had been run on the development
machine at stage 3b (Snakemake is not installed there); CI's `workflow-dry-run` job is
their first run.

Likewise there is no L4 in the pipeline at all, and L5 has no materialised
artifact — `gliostrata optimize` calls `predict_survival` internally and never
writes the predictions out.

### Infrastructure gaps

- **No container.** No `Dockerfile` and no image digest. There is a lock file
  since 2026-11-27 (`requirements-lock.txt`, which CI's `test` job installs
  from), so two clean installs from it resolve the same 37 distributions.
  `requirements.txt` and `environment.yml` still state floors, and `snakemake`
  and `dvc` are not locked. A reproducibility release needs an image digest.
- **No DVC remote.** Every derived artifact exists on one machine. §3.
- **No GDC data-release pinning.** §1.2.
- **No Zenodo DOI**, and nothing to mint one for yet.
- **No end-to-end run on the real cohort recorded anywhere.** L0 to L3 have run on
  TCGA (`build-feature-store --mode real`); L4 to L7 have not, and `dvc repro` on
  real data still stops at the stub `profile` stage. Every real-cohort derived
  artifact is git-ignored and untracked; their checksums are in
  `docs/project_record/06_evidence_reports.md`.
- **Nothing in L2 is validated against a published reference.** The single-sample
  subtype agrees with the published TCGA labels for 56 of 76 v2-test patients
  (0.737), a descriptive check on the cohort its centroids came from, not a
  validation. The state scores recover the imposed proportions of 100 simulated
  Neftel 2019 mixtures (T1's deconvolution row, mean per-state r 0.9521, run before
  deposit; `docs/t1_deconvolution_spec.md`), an optimistic recovery test, not a
  measurement on tumours. The TME scores have no ground truth here. T1's subtype row
  has not been run.

### What *is* genuinely reproducible today

To be fair to the scaffold, the following are real and verified:

1. The offline fixture pipeline, bit-reproducible from a seed (§1.1).
2. The test suite: 1126 tests, 1125 passing and 1 skipped, with all outbound sockets blocked.
3. **L0 ingestion** — TCGA-GBM clinical, bulk RNA-seq, IDH mutation and MGMT
   methylation calls (real GDC API code with cohort filtering, barcode parsing and
   a GENCODE gene annotation), and CGGA RNA-seq. This is the one layer that does
   what it says end to end.
4. **The partition freeze/verify machinery** — seeded hash splitting, sealed
   manifests, content-hash drift detection. It has already caught a real cohort
   change (§4), which is the only way to know a guardrail works.
5. The Snakemake DAG and the DVC stage graph, both machine-validated, and
   tested against the code so they cannot silently rot (`tests/test_workflow.py`).
6. **The TCGA feature store, L1 to L3** — `gliostrata harmonize-clinical`, then
   `gliostrata build-feature-store --mode real`, rebuild it from the files on disk:
   all 594 patients, one row each; a label table; a reason code for every missing
   value; the column metadata; and a report of the three views (profiling 283,
   transferable model 570, molecular model 298) with per-split flow counts and the
   SHA-256 of every input. A repeat run is byte-identical
   (`docs/project_record/06_evidence_reports.md`, stage 3b).
   `scripts/build_tcga_feature_store.py` is a thin wrapper over the second command.
   The values are measured; they are not validated.

---

## 7. Workspace-wide tools and the held-out cohort

The read guard governs scripts (`docs/project_record/06_evidence_reports.md`, stage
2). It does not govern the tools an editor or an AI agent points at the whole
repository. On 2026-11-27T05:48:55Z a workspace-wide text search, run by the review
agent at the start of a review, opened the five extracted CGGA text files under
`data/raw/cgga/_cache/extracted/` and scanned them for two path strings. It returned
nothing and no value was used, but it was an access that no project script made
(access ledger entries 8 and 9).

What is in place since then:

- `.rgignore`: ripgrep's ignore file. It applies to `rg` and to tools built on it
  that keep ripgrep's default ignore handling. `rg --no-ignore` and `rg -uu` skip it.
- `.clineignore`: the agent's ignore file, for its file listing, reading and search
  tools.
- `.clinerules`: standing instructions for any AI agent working in this repository.
  The first one is: never run a workspace-wide search here; scope every search to
  `src`, `tests`, `scripts`, `docs`, `configs` or `workflow`.
- `data/raw/**` and `data/processed/**` are git-ignored, which keeps them out of git
  and out of any tool that respects `.gitignore`.

Not verified: that the installed agent's search tool honours the two ignore files
during its file walk. A live test was not run, because a failing test would open the
held-out files again. The instruction in `.clinerules` is the control that does not
depend on it.

Still open: the plain-text copies themselves. The intended state is to keep only the
compressed sources and to extract into a temporary directory inside a guarded read,
deleting it afterwards. It is deferred, with its reason, in
`docs/project_record/05_future_plan.md`, section 5.5.

