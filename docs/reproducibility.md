# Reproducibility

Point 10 of the GlioStrata defensibility stack is a full reproducibility
release: code, containers, manifests, and a Zenodo DOI. This document is the
honest status report on that goal. It is deliberately written so that the
"not yet reproducible" section is as prominent as the instructions.

**Summary in one line:** the *plumbing* is reproducible end to end today; the
*science* is not, because L1, L2, L5 and L6 are stubs.

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

`make pipeline-fixture` writes, under `.test/`:

| Layer | Artifact |
| --- | --- |
| L0 | `.test/interim/fixture_synthetic_cohort.parquet` |
| L1 | `.test/interim/fixture_harmonized.parquet` |
| L2 | `.test/processed/fixture_profiled.parquet` |
| L3 | `.test/processed/fixture_feature_store.parquet` |
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
snakemake --cores 4 --snakefile workflow/Snakefile   # the whole thing
```

Or, equivalently, `snakemake --cores 4 --snakefile workflow/Snakefile`, which
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

`dvc.yaml` is a correct, parsing description of the seven real pipeline stages.
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
possible from git alone (§5).

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

**Open decision for the project lead:** v1 is now superseded. Either retire it
(with a note recording that nothing was ever trained on it) or document it as
the partition belonging to the 583-row cohort snapshot. Leaving two live
manifests with no stated precedence is how the wrong one gets used. The CI
integrity check deliberately cannot catch this — both manifests are internally
valid — so it needs a human decision.

---

## 5. Continuous integration

`.github/workflows/ci.yml`, on every push and pull request, Python **3.14**
(the version this project is developed against; `pyproject.toml` declares
`requires-python = ">=3.10"`).

| Job | What it does |
| --- | --- |
| `lint` | `ruff check --select E9,F63,F7,F82` as a **blocking** gate, plus the full project config as **advisory**. See the note below. |
| `test` | `pip install -e ".[dev]"` then `python -m pytest -q`. |
| `partition-drift` | Partition manifest integrity: sealed, digests match, and every id re-derives its frozen split from `(seed, patient_id)`. Plus a full `--verify` when the cohort happens to be in the checkout. |
| `workflow-dry-run` | `snakemake -n` on the real-cohort DAG, then a real offline fixture run, then `dvc stage list` / `dvc dag`. |

### CI never touches the network

Verified, not assumed. The suite was re-run with `socket.socket.connect`,
`connect_ex`, `create_connection` and `getaddrinfo` all patched to raise. The
results were **identical** to the normal run — no test attempted a connection.
Specifically:

- `tests/test_tcga_download.py` and `tests/test_tcga_expression.py` load
  `scripts/download_tcga_*.py` by path and drive the pure functions with
  synthetic fixtures. They assert on the GDC URL *constants* without ever
  requesting them.
- The three tests in `tests/test_cox.py` that read the real cohort are guarded
  by `@pytest.mark.skipif(not REAL_COHORT.exists())`. `data/raw/**` is
  gitignored, so they skip in CI by construction.
- The `workflow-dry-run` job uses `source=synthetic`, which routes L0 through
  `gliostrata synthetic-data` and never reaches the two ingest rules.

**No test requires network access.**

### Why the blocking lint set is narrow

The repo currently has ~96 ruff findings under its own config (90 in
`src`+`tests`, 6 in `scripts`), dominated by `N806`
(non-lowercase-variable-in-function) — which is arguable in statistical code
where `X`, `X_train` and `N` are the conventional names. Making the full set
blocking today would land a red CI on the first commit and train everyone to
ignore it. The blocking set is the "the code is broken" floor (syntax errors,
invalid comparisons, misplaced statements, undefined names) and passes clean.
The advisory step prints the backlog. Widen the blocking set as it is cleared.

Separately: `pyproject.toml` still declares `[tool.ruff] select`/`ignore` at
the top level, which modern ruff warns is deprecated in favour of
`[tool.ruff.lint]`. Harmless today, but the warning appears on every run.

---

## 6. What is NOT reproducible

Everything below is **plumbing that works, wrapped around science that does
not exist yet.** The DAG runs green. That is a statement about the DAG.

### Stubbed layers

| Layer | Module | What it actually does | What the architecture requires |
| --- | --- | --- | --- |
| **L1** harmonize | `gliostrata.data.harmonize` | Schema normalisation and basic QC row filtering. | Cross-cohort batch correction (ComBat/limma), platform harmonization, imputation of missing molecular covariates. |
| **L2** profile | `gliostrata.profiling.{subtype,states,tme}` | Derives "subtypes", Neftel state fractions and TME composition **from the clinical table**. Reads no expression matrix. | Verhaak classification, Neftel state deconvolution and TME deconvolution from `counts.parquet`/`tpm.parquet`. |
| **L5** survival | `gliostrata.prediction.survival` | Produces survival numbers without a model fitted on the real cohort. | A validated survival model with calibration and a reported C-index. |
| **L6** optimize | `gliostrata.optimization.bayesian` | Ranked enumeration of configurations. | Bayesian optimization against a validated response surface. |
| **L7** handoff | `gliostrata.handoff.protocol` | Structurally complete JSON documents templated from the above. | Bench-executable protocols. |

**The single most important consequence:** the L0 expression ingest is real and
fully implemented, and it produces `counts.parquet` and `tpm.parquet` — and
**nothing downstream reads them.** Every subtype call, cell-state fraction and
TME fraction in this pipeline today is computed from clinical covariates. They
are not measurements.

This is why `workflow/rules/profile.smk` declares `HARMONIZED` as its only
input and lists the expression matrices as a top-level target instead of an
input edge: an input edge would make the DAG assert a data dependency the code
does not have. The DAG is a map of the code, not of the intent.

Likewise there is no L4 in the pipeline at all, and L5 has no materialised
artifact — `gliostrata optimize` calls `predict_survival` internally and never
writes the predictions out.

### Infrastructure gaps

- **No container.** No `Dockerfile`, no pinned lockfile. `requirements.txt`
  and `environment.yml` use floors (`>=`), not pins, so two clean installs a
  month apart can resolve to different versions. A reproducibility release
  needs a lockfile and an image digest.
- **No DVC remote.** Every derived artifact exists on one machine. §3.
- **No GDC data-release pinning.** §1.2.
- **No Zenodo DOI**, and nothing to mint one for yet.
- **No end-to-end run on the real cohort recorded anywhere.** `data/manual/`
  holds a synthetic run's output; `data/interim/tcga_gbm_harmonized.parquet` is
  the only real-cohort derived artifact present, and it is gitignored and
  untracked.
- **Nothing in L2 is validated against a published reference.** Until L2 reads
  expression and reproduces Verhaak/Neftel calls on a benchmark cohort, there
  is no external check on the pipeline's central claim.

### What *is* genuinely reproducible today

To be fair to the scaffold, the following are real and verified:

1. The offline fixture pipeline, bit-reproducible from a seed (§1.1).
2. The test suite, 371 tests, passing with no network.
3. **L0 TCGA-GBM ingestion**, both clinical and bulk RNA-seq — real GDC API
   code with cohort filtering, barcode parsing and a GENCODE gene annotation.
   This is the one layer that does what it says.
4. **The partition freeze/verify machinery** — seeded hash splitting, sealed
   manifests, content-hash drift detection. It has already caught a real cohort
   change (§4), which is the only way to know a guardrail works.
5. The Snakemake DAG and the DVC stage graph, both machine-validated, and
   tested against the code so they cannot silently rot (`tests/test_workflow.py`).
