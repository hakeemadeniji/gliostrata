# GlioStrata: Measured variables

Supplement to the OSF registration ("Measured variables" field). It summarises section 5
and the data documentation of `docs/preregistration.md` at the commit tagged
`prereg-deposit`; that file is authoritative.

All data are existing, de-identified public data. No variable is manipulated in any
patient.

## 1. Development cohort: TCGA-GBM (NCI Genomic Data Commons)

594 eligible patients, from 617 case records.

| Variable | Source | Units and coding | Availability | Role |
| --- | --- | --- | --- | --- |
| `patient_id` | GDC case submitter id | TCGA barcode (TCGA-XX-XXXX) | all | key |
| `age` | diagnoses, age at diagnosis | years (days / 365.25, rounded) | all | predictor |
| `sex` | demographic | M / F | all | predictor |
| `kps` | Karnofsky performance score | 0 to 100 | missing for 150 of the 592 modelled patients; 0 at source for 3 (disclosed, not corrected) | predictor (committed model only) |
| `treatment_arm` | treatment records of the first diagnosis | ChemoRadiation / MonoTherapy / NoAdjuvant | missing for 23 | predictor |
| `overall_survival_months` | days to death or last follow-up | months (days / 30.4375) | patients with <= 0 excluded | outcome |
| `vital_status` | demographic | 1 = dead, 0 = alive (censored) | missing for 2 | outcome |
| `idh_mutated` | simple nucleotide variants (MAF) | 1 = IDH1 R132 / IDH2 R172 hotspot, 0 = sequenced without hotspot | called for 368; 225 not tested; 1 sequenced on a recurrent sample only (tested, not callable) | predictor (molecular variant) |
| `mgmt_methylated` | 450k / EPIC methylation | 1 = mean beta of cg12434587 and cg12981137 >= 0.35 | called for 415; 179 not tested | predictor (molecular variant) |
| RNA-seq | STAR - Counts (GENCODE v36) | read counts and TPM, 60,660 genes | 284 samples, 283 patients | expression input |
| Published subtype label | Brennan et al. 2013 table | Classical / Mesenchymal / Proneural / Neural | 246 labelled patients with RNA-seq | reference for the subtype centroids |
| Extent of resection | not present at source for any patient | | none | not used |

Split: `tcga_gbm_clinical_v2` (sealed), 425 train / 169 test, assigned by a hash of
(seed 20251126, patient id).

## 2. External cohort: CGGA (Chinese Glioma Genome Atlas)

305 unique patients after the eligibility rules (333 rows). The sealed external partition
`cgga_external_v1` holds 218 of them.

| Variable | Source field | Units and coding | Availability | Role |
| --- | --- | --- | --- | --- |
| `patient_id` | CGGA_ID | CGGA pseudonymous sample code | all | key |
| `cgga_dataset` | batch | mRNAseq_693 / mRNAseq_325 / mRNA_array_301 | all | stratum (RNA-seq batch), eligibility |
| `age` | Age | years | all | predictor |
| `sex` | Gender | M / F | all | predictor |
| `treatment_arm` | Radio_status, Chemo_status | as TCGA; missing if either flag is missing | 283 of 305 with the T2 features and an outcome | predictor |
| `overall_survival_months` | OS (days) | months (days / 30.4375) | <= 0 excluded | outcome |
| `vital_status` | Censor | 1 = dead, 0 = censored | | outcome |
| `idh_mutated`, `mgmt_methylated`, 1p/19q | CGGA's own annotations | 1 / 0 | 241 with IDH and MGMT | covariates (secondary model) |
| Histology, grade, PRS type | clinical file | GBM; WHO IV; primary / recurrent / secondary | all | eligibility only |
| `kps` | none | | null for every patient | not usable |
| RNA-seq counts | Read_Counts-genes (20220620) | STAR read counts, 55,523 gene symbols | 225 patients, of whom 218 have an outcome | expression input |
| RNA-seq FPKM | RSEM-genes (20200506) | RSEM FPKM, about 24,000 genes | all 218 partition patients | primary T2 comparator input |
| Reference Verhaak subtype | TCGA_subtypes | 4 classes | 108 microarray patients | T1 reference label |
| Microarray expression | mRNA_array_301 | normalised intensities | not yet obtained | T1 subtype input |

## 3. Other sources

| Data | Variables | Status |
| --- | --- | --- |
| Neftel et al. 2019 single-cell atlas (GEO GSE131928, sample GSM3828672; 3CA annotations) | per-cell log2(TPM/10+1) expression (23,686 genes x 7,930 cells); malignant flag; six meta-module scores | obtained; T1 deconvolution only |
| PDX growth curves | tumour volume per time point, control arm, subtype or RNA-seq per model | not yet obtained; due 2027-01-15 (pre-registration 2.2 item 3) |
| In vivo NSC tracking | migration measurements per animal | not yet obtained |
| T4/T5 simulation output | each patient's predicted survival curve under a configuration and under the Stupp protocol | produced by the simulator frozen under 2.2 item 11 |

## 4. Variables used for subsetting or exclusion

The following are used only to subset or exclude patients:
- histology, grade and primary/recurrent status (eligibility);
- availability of an expression profile, an outcome or a reference label (scoring sets);
- the development split (v2 train / test);
- the CGGA RNA-seq batch (the stratum of the improvement-row C-index).

Every exclusion is counted and reported.
