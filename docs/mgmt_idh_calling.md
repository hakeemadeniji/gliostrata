# MGMT and IDH calling in TCGA-GBM

How `scripts/download_tcga_mutations.py` (IDH) and
`scripts/download_tcga_methylation.py` (MGMT) turn raw GDC files into
per-patient calls, what those calls do and do not mean, and what is still
unverified. Read this before using either call table in a model.

**Summary in one line:** both scripts produce a *three-state* call
(1 / 0 / not callable) on the primary tumour; the IDH call is the WHO 2021
hotspot definition and is trustworthy as a genotype; the MGMT call is a
**transparent two-probe score, explicitly NOT MGMT-STP27**, and its binary
form depends on an uncalibrated threshold.

| | IDH | MGMT |
| --- | --- | --- |
| Input | Open masked somatic MAF (WXS) | Methylation beta values (27k / 450k) |
| Definition | IDH1 R132 / IDH2 R172 hotspot variant | Mean beta of cg12434587 + cg12981137 >= 0.35 |
| Equals the protocol's method? | Yes | **No** (MGMT-STP27 coefficients unverified) |
| Primary output | `idh_mutated` | `mgmt_score` (continuous) |
| Join to cohort | Separate, reviewed step | Separate, reviewed step |

---

## 1. Why these scripts exist

TCGA-GBM is the development cohort and CGGA the external validation cohort,
and they barely share molecular features. CGGA publishes IDH status for
326/333 patients and MGMT status for 312/333; the TCGA *clinical* ingest
carries neither. MGMT is the covariate most tied to temozolomide response and
IDH status re-defines the diagnosis under WHO 2021, so a model developed on
TCGA without them is missing the two most consequential molecular variables.

Both scripts only **download and call**. Neither fits anything, joins
anything, or modifies `clinical_cohort.parquet`.

---

## 2. IDH status (`scripts/download_tcga_mutations.py`)

### 2.1 Definition

IDH-mutant glioma is defined by variants at two codons:

| Gene | Codon | Common variants |
| --- | --- | --- |
| IDH1 | R132 | R132H (dominant), R132C, R132G, R132S, R132L |
| IDH2 | R172 | R172K (dominant), R172M, R172W, R172G |

A MAF row counts as a hotspot hit only if **all** of the following hold:

- `Hugo_Symbol` is `IDH1` or `IDH2`;
- `Protein_position` is 132 (IDH1) or 172 (IDH2) **and** the annotated
  reference residue in `Amino_acids` is arginine (`R`);
- `Variant_Classification` is `Missense_Mutation`, `In_Frame_Ins` or
  `In_Frame_Del`.

Nonsense, frameshift, splice-site and silent changes at those codons are
**not** called. They do not create the neomorphic enzyme activity that
defines IDH-mutant glioma, and calling them would inflate the mutant rate.
Requiring the reference residue to be `R` guards against a numbering offset:
a record that says position 132 but reference `G` is not this hotspot.

### 2.2 Three states, not two

| `idh_mutated` | `call_basis` | Meaning |
| --- | --- | --- |
| 1 | `hotspot_variant` | a hotspot variant was found in the tumour |
| 0 | `sequenced_no_hotspot` | the tumour **was** sequenced; no hotspot found |
| `<NA>` | `no_mutation_data` | no primary-tumour MAF exists for the patient |

The third state is **not** wildtype. Coercing it to 0 invents patients who
were never tested and silently shrinks the mutant fraction. It is stored as a
nullable `Int64` and counted separately in every report.


### 2.3 Which sample, and which file

- Only **primary tumour** files are used (sample-type code `01`). The 12
  recurrent-tumour MAFs in TCGA-GBM are dropped.
- GDC lists a file's linked samples in an arbitrary order, and for many
  tumour/normal pairs it lists the **normal blood sample first**. Reading
  `samples[0]` therefore labels the file as normal and the primary-tumour
  filter throws it away. The script instead looks at every linked aliquot,
  picks the best tumour one, and then cross-checks that choice against the
  `Tumor_Sample_Barcode` written inside the MAF itself. If metadata and file
  content disagree about the patient, the run stops rather than guessing.

### 2.4 WHO 2021 implication (reported, not resolved)

Under WHO 2021 an IDH-mutant grade IV astrocytic tumour is *astrocytoma,
IDH-mutant, CNS WHO grade 4*, **not glioblastoma**. Every IDH-mutant case in
this "GBM" cohort is therefore a possible mis-grouping inherited from the
TCGA-era diagnosis. The script lists those patients in its report and on
stdout. It deliberately does **not** remove them: whether to exclude them
from the development cohort is a pre-registration decision.

---

## 3. MGMT promoter methylation (`scripts/download_tcga_methylation.py`)

### 3.1 This is NOT MGMT-STP27

The protocol asks for the published MGMT-STP27 logistic model. Its exact
coefficients and intercept **could not be verified against the primary
publication** during development, and an invented coefficient is far worse
than an honest substitute. So the script computes something simpler and
labels it in every output (`call_method = two_probe_mean_beta`,
`is_mgmt_stp27 = false`, and a `NOT MGMT-STP27` line on stdout).

> **TO BE VERIFIED BY AUTHOR:** the MGMT-STP27 coefficients, intercept and
> decision threshold, from the primary publication's supplementary material.
> Once obtained, implement them in `mgmt_score()` behind a new method name,
> change `CALL_METHOD`, and re-run. Nothing else needs to change.

### 3.2 The score that is actually computed

```
score       = mean beta of cg12434587 and cg12981137
methylated  = score >= threshold          (default 0.35)
```

- Both probes are required. If either is absent from the file or is `NA`,
  the patient is **not callable**; the score is never averaged over one
  probe, because a one-probe score would be silently compared with two-probe
  scores.
- The **continuous `mgmt_score` is the primary output**. The binary call is a
  convenience derived from it, and anyone can re-threshold the score.
- The threshold **0.35 is an a-priori working value** fixed before the cohort
  was downloaded. It is not a published cut-point and is not calibrated. The
  report prints the methylated fraction across thresholds 0.20 to 0.60 so the
  sensitivity of the call is visible; see section 6.

### 3.3 Platforms

TCGA-GBM was profiled on two arrays and they are not interchangeable:

| Platform | Patients with data | Notes |
| --- | --- | --- |
| Illumina 450k | 142 | the denser array; MGMT-STP27's native platform |
| Illumina 27k | 287 | older; most of this cohort's calls come from here |
| Both | 5 | the 450k file is used |

Beta values are used exactly as GDC distributes them, with **no
cross-platform harmonisation**. A call made on 27k and a call made on 450k
are therefore not guaranteed to agree. Every call records its `platform` for
this reason: stratify on it, or include it as a covariate, whenever the calls
are used.

### 3.4 Three states, not two

| `mgmt_methylated` | `call_basis` | Meaning |
| --- | --- | --- |
| 1 | `two_probes_measured` | both probes measured, score >= threshold |
| 0 | `two_probes_measured` | both probes measured, score < threshold |
| `<NA>` | `required_probe_absent_or_NA` | a file exists but a probe is missing / `NA` |
| `<NA>` | `no_methylation_data` | no primary-tumour beta file for the patient |

A false "unmethylated" call is worse than no call: it assigns a patient to
the worse-prognosis, less-temozolomide-responsive group on the strength of
data that does not exist.

---

## 4. Outputs

Both scripts write the same three artefacts under `--out`, plus a cache:

| File | Content |
| --- | --- |
| `file_manifest.parquet` | one row per GDC file, with parsed barcodes |
| `idh_calls.parquet` / `mgmt_calls.parquet` | one row per patient, including not-callable patients |
| `call_report.json` / `mgmt_call_report.json` | every count printed to stdout, machine-readable |
| `_cache/` | raw downloads; scratch, never an input to a model |

Default locations: `data/raw/tcga_gbm/mutations/` and
`data/raw/tcga_gbm/methylation/`.

`idh_calls.parquet`: `patient_id`, `idh_mutated` (nullable `Int64`),
`idh_variants` (e.g. `IDH1 p.R132H`), `n_idh_hotspot_variants`, `call_basis`,
`n_somatic_variants`, `sample_barcode`, `aliquot_barcode`, `file_id`.

`mgmt_calls.parquet`: `patient_id`, `mgmt_methylated` (nullable `Int64`),
`mgmt_score`, `cg12434587`, `cg12981137`, `call_basis`, `call_method`,
`platform`, `threshold`, `sample_barcode`, `aliquot_barcode`, `file_id`.

The row set is every patient with a primary-tumour file **plus** every patient
of `clinical_cohort.parquet` (unless `--no-cohort`), so a missing call is an
explicit row with a reason rather than an absent row.

---

## 5. Known limitations and open decisions

1. **MGMT-STP27 is not implemented** (section 3.1). Until the author verifies
   the coefficients, the MGMT call is a two-probe score.
2. **The MGMT threshold is uncalibrated.** 0.35 was fixed a priori. Use the
   continuous score where possible; if a binary call is unavoidable, run the
   sensitivity analysis across the reported threshold grid.
3. **No cross-platform harmonisation.** Beta values are used as distributed
   by GDC, and most of this cohort's calls come from the older 27k array.
4. **An IDH wildtype call means "no hotspot in the open-access MAF".** It is
   not an orthogonal assay result; a low-tumour-purity sample can read as
   wildtype. It is nevertheless the correct genotype definition, and the
   observed rate (6.2% of tested, section 6) is a single-digit figure, the
   level the script's own sanity check expects for primary, de novo GBM.
5. **Recurrent-tumour samples are not used.** A patient whose only sequenced
   sample is recurrent is not callable for IDH.
6. **IDH-mutant cases stay in the primary cohort.** Decided: see the
   pre-specified sensitivity analysis in section 6.3.
7. **Neither table is joined to the cohort.** They are standalone; the join is
   a separate, reviewed step.
8. **Coverage differs between the two markers.** IDH and MGMT are testable on
   different subsets of patients (section 6); a model using both sees only
   their intersection.
9. **Not yet wired into the pipeline.** Neither script is a `dvc.yaml` stage
   or a Snakemake rule, and the call tables are not under DVC. They exist as
   scripts with offline tests; integrating them is a separate step.

---

## 6. Results of the real run

Measured against the live GDC API. Every number below is printed by the
scripts and stored in the `*_report.json` files.

### 6.1 IDH

| | |
| --- | --- |
| Open masked MAFs in TCGA-GBM | 464 files, 374 patients, 18.6 MB |
| Primary-tumour files kept | 452 (12 recurrent-tumour files dropped) |
| Patients with a primary-tumour MAF | 373 |
| IDH-mutant | 23 (6.2% of tested) |
| IDH-wildtype | 350 (93.8% of tested) |
| Hotspot variants observed | IDH1 p.R132H, p.R132C, p.R132G (no IDH2) |

| Patient set | n | with IDH call | mutant | still missing |
| --- | --- | --- | --- | --- |
| clinical cohort | 594 | 368 | 23 | 226 |
| with expression | 284 | 248 | 16 | 36 |
| cohort AND expression | 283 | 247 | 16 | 36 |

Five of the 373 sequenced patients are not in the 594-patient clinical
cohort; that is the whole difference between 373 and 368.

The mutant rate is the single-digit figure expected for primary, de novo
glioblastoma, and well under the script's 10% red-flag line. All 23
IDH-mutant patients are listed in `call_report.json` under
`summary.mutant_patient_ids`, next to `who_2021_note`, because of section 2.4.

### 6.2 MGMT

| | |
| --- | --- |
| Open beta-value files in TCGA-GBM | 450 files, 423 patients, 2.26 GB |
| Primary-tumour files kept | 435 (13 recurrent, 2 solid-tissue-normal dropped) |
| Patients with a primary-tumour file | 422 (5 profiled on both platforms; 450k used) |
| Callable | 422 of 422 |
| Methylated at the default 0.35 | 128 (30.3% of callable) |
| By platform | 27k: 94 / 282 (33.3%) -- 450k: 34 / 140 (24.3%) |

Both probes were **measured in every file on both platforms** (0 `NA`, 0
absent across 282 27k and 140 450k files), so no patient with a file was
lost to probe coverage. Every file downloaded; none failed.

| Patient set | n | with MGMT call | methylated | still missing |
| --- | --- | --- | --- | --- |
| clinical cohort | 594 | 415 | 128 | 179 |
| with expression | 284 | 230 | 68 | 54 |
| cohort AND expression | 283 | 229 | 68 | 54 |

**The binary call is highly threshold-dependent.** Methylated count among
the 422 callable patients:

| threshold | 0.20 | 0.25 | 0.30 | **0.35** | 0.40 | 0.45 | 0.50 | 0.60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| methylated | 176 | 164 | 151 | **128** | 104 | 93 | 84 | 47 |
| % of callable | 41.7 | 38.9 | 35.8 | **30.3** | 24.6 | 22.0 | 19.9 | 11.1 |

The score distribution is bimodal but lopsided: 230 patients sit below 0.10,
then 16 in 0.10-0.20, 83 in 0.20-0.45 and 93 at 0.45 or above. Two things
follow, and both are stated rather than fixed:

- **The default 0.35 does not sit in the distribution's emptiest region.**
  The sparsest gap above the unmethylated mode is 0.10-0.20; 0.35 falls on
  the slope of the middle band, where 23 and 24 patients lie in the two bins
  either side of it. Moving the cut from 0.35 to 0.20 changes the call for
  48 patients (11% of the callable set). The threshold was fixed a priori
  and is **left unchanged** here; this table exists so the choice is visible.
- **The platform split matters.** The 27k call is methylated for 33.3% and
  the 450k call for 24.3%, with no harmonisation between them.

Practical consequence: use the continuous `mgmt_score` wherever the model
allows it, and report the threshold grid alongside any binary result. The
script's own sanity band (30-60%, taken from the project brief and not
verified against a primary source) contains the default-threshold figure at
its lower edge.

### 6.3 Both markers together, and the IDH sensitivity analysis

| Patient set | n | IDH call | MGMT call | **both** | either | neither |
| --- | --- | --- | --- | --- | --- | --- |
| clinical cohort | 594 | 368 | 415 | **314** | 469 | 125 |
| with expression | 284 | 248 | 230 | **206** | 272 | 12 |
| cohort AND expression | 283 | 247 | 229 | **205** | 271 | 12 |

Of the 284 expression patients, **206 gain both an IDH and an MGMT call**,
and only 12 gain neither. Among those 206:

| | MGMT unmethylated | MGMT methylated |
| --- | --- | --- |
| IDH wildtype | 135 | 55 |
| IDH mutant | 10 | 6 |

**Pre-specified sensitivity analysis (decision recorded; not yet run).**
All 23 IDH-mutant patients stay in the primary cohort and the frozen
594-patient partition stands. Under WHO 2021 those cases would be
reclassified as astrocytoma, IDH-mutant (section 2.4), so the primary
analysis is repeated with them excluded as a **robustness check on the
primary analysis, not a change to it**:

- re-report the primary estimand (Cox C-index) and the subtype distribution
  with the 23 excluded (16 of the 284 expression patients);
- a change in C-index larger than **0.01** (absolute) is flagged prominently
  in the main results, not relegated to a supplement; a smaller change is
  reported in one paragraph;
- the exclusion is applied at analysis time from `idh_calls.parquet`.
  Nothing is removed from any ingested table.

This cannot be run yet: the call tables are not joined to the cohort
(section 5, item 7).

---

## 7. Reproducing

```bash
python scripts/download_tcga_mutations.py --dry-run     # metadata only, no download
python scripts/download_tcga_mutations.py               # ~19 MB, about 3 minutes

python scripts/download_tcga_methylation.py --dry-run
python scripts/download_tcga_methylation.py --limit 12  # a few real files
python scripts/download_tcga_methylation.py             # ~2.3 GB, 15-30 minutes
python scripts/download_tcga_methylation.py --threshold 0.30   # re-threshold
```

- Both scripts resume from `_cache/`. Re-running with a different threshold
  re-reads cached files; nothing is downloaded twice.
- The choice of file per patient is a fixed sort order (patient, sample type,
  platform, aliquot barcode, file id), so it is the same on every rerun and
  independent of the order GDC returns.
- `tests/test_tcga_mutations.py` and `tests/test_tcga_methylation.py` run
  offline against synthetic files in the real GDC layout.

---

## 8. Harmonisation with CGGA

The two cohorts do not measure these markers the same way:

| | TCGA-GBM (this document) | CGGA (`scripts/download_cgga.py`) |
| --- | --- | --- |
| IDH | hotspot genotype called from the MAF | published label `Mutant` / `Wildtype` |
| MGMT | two-probe beta score, uncalibrated threshold | published label `methylated` / `un-methylated` |

- **IDH** is comparable at the level of a three-state indicator, on the
  assumption that CGGA's `Mutant` means a hotspot variant. The downloaded
  CGGA files do not say how that label was determined, so the assumption
  cannot be checked here.
- **MGMT is not the same measurement in the two cohorts.** The TCGA label is
  derived from a score at an a-priori threshold; the CGGA label is whatever
  assay the consortium applied. A model trained on the TCGA label and
  validated on the CGGA label is trained and validated on different
  instruments. Any external-validation result that depends on MGMT inherits
  this mismatch and should say so.
- Both CGGA labels are parsed with negatives tested first, because
  `un-methylated` contains `methylated` as a substring.
