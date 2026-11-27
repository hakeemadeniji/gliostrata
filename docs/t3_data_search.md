# T3 data search (2026-11-27)

Scope: the two T3 rows (tumour-growth RMSE vs PDX curves; NSC-migration r vs in vivo
tracking). Registered deadline: data in hand by 2027-01-15. The qualifying rule is in
docs/preregistration.md, section 2.2 item 3. Nothing was downloaded into the repository.
Small public files were fetched to the session scratchpad only (folders `alcaniz`, `vaubel`,
`epmc`, `cm`). Labels: VERIFIED = quote plus URL read in this session; UNVERIFIED = not
confirmed.

## 1. Summary and recommendation

- No source found so far meets the full item 3 rule. The two named candidates both fail
  the "four Verhaak subtypes" clause on what they publish, and neither is confirmed to
  carry per-time-point control volumes in an open file.
- Likely blocker, flagged for the lead: both panels classify **three** subtypes
  (classical, mesenchymal, proneural), using the Wang 2017 scheme, which has no
  Neural class. Any modern GBM PDX panel is likely to do the same. The clause
  "span the four Verhaak subtypes" may therefore be unmeetable by any PDX panel as
  written. I did not change the pre-registration. A decision on this belongs to the
  lead, and any change must be a dated, declared amendment made before the data are
  seen. Otherwise the row closes NOT EVALUABLE, which the protocol already permits.
- Best growth path: ask the Mayo Clinic Brain Tumor PDX National Resource (Table S3B of
  Vaubel et al. holds per-line growth data, and RNA-seq is public under PRJNA548556) and
  the Alcaniz/Hoffmann group (Charite/EPO Berlin; per-time-point control volumes are
  probably in the figures only). Send both requests this week.
- NSC migration: no quantitative in vivo dataset found. The prior verdict stands.
  Best path is a request to the Rust 2022 authors (raw data "on request") and to the
  Aboody lab. Expect NOT EVALUABLE unless one replies with per-animal migration
  distances.
- Not verified, for access reasons: NCI PDMR, CancerModels.org, PRoXE (details below).

## 2. PDX growth candidates

| Candidate | Data that exist | Where | Access | Longitudinal volumes | Control arm | Subtype or RNA-seq per model | Four Verhaak subtypes, >= 10 curves |
|---|---|---|---|---|---|---|---|
| Mayo GBM PDX (Vaubel et al. 2020, Clin Cancer Res) | 96-model panel; time to exceed 2 cm3 per line; engraftment rates; RNA-seq | PMC7056576; Mayo resource site | Models and genomic data "openly available to investigators" at minimal cost; request through the resource (site returned 403 to automated fetch) | UNKNOWN: only summary kinetics seen ("median time to exceed 2 cm3 ranging from 20 to 189 days"); per-line data are in Table S3B, not read | UNKNOWN (untreated growth lines; not a drug-trial control arm) | MET for RNA-seq: "RNAseq Accession PRJNA548556"; subtype calls in Table S6A (not read) | NOT MET as published: "classical (n=31, 46%), mesenchymal (n=21, 31%), or proneural (n=16, 23%)", Wang 2017 scheme; >= 10 curves plausible |
| Alcaniz 2023 (Front Oncol) | 26 s.c. GBM PDX; T/C vs vehicle controls; RNA-seq for 23 | doi:10.3389/fonc.2023.1129627, PMC10126369; RNA-seq at EGA EGAS00001007119 | RNA-seq via EGA (controlled access, data access committee; not open). Supplement is open | NOT confirmed. Supplement is Table_1.xlsx (a mutation list), four image files and figure images; no volume table | MET: "Compared to their respective vehicle treated control group..." | MET for RNA-seq in 23 models: "RNA sequencing... was performed for 23 established GBM PDX models" | NOT MET: "mesenchymal, proneural and classical" |
| NCI PDMR | Models, molecular data | dctd.cancer.gov PDMR page | "The models are made available upon request" | UNVERIFIED: page does not mention growth curves | UNVERIFIED | UNVERIFIED | UNVERIFIED |
| CancerModels.org / PDX Finder | Model metadata | cancermodels.org | API | UNVERIFIED: API returned HTTP 500 on both attempts | UNVERIFIED | UNVERIFIED | UNVERIFIED |
| PRoXE | Leukaemia and lymphoma PDX | proxe.org | n/a | Not a glioma resource (search result: "mainly a leukemia and lymphoma resource") | n/a | n/a | NOT MET (no GBM found) |
| maeve R package, gdc0084_temo_01 | 1411 per-animal, per-day rows, s.c. GBM PDX, drug arms | rdrr.io/github/wfforrest/maeve | Open | MET for this dataset | UNKNOWN (the description lists drug arms only) | UNKNOWN | NOT MET: appears to be a single study (model count unconfirmed) |

Notes.
- Item 2 of Vaubel: the task text cites "Cancer Research"; the paper is Clin Cancer Res
  (doi:10.1158/1078-0432.CCR-19-0909). VERIFIED by Europe PMC.
- The Mayo page reporting that all 96 models are available is VERIFIED via PMC7056576.
  The resource's own site could not be read (403), and I did not work around that.
- Mayo growth data are "time to 2 cm3" style (per line), not necessarily mean volume
  per time point. The item 3 unit (at least 3 time points per curve) may need the raw
  per-animal data from the resource, which means a request.
- Alcaniz: EGA is controlled access, so RNA-seq needs a DAC application (weeks).
  The subtype call is already computed in the paper (single-sample GSEA on published
  gene sets), which may be enough for the "subtype call" arm of the rule.
- PDMR, CancerModels.org and PRoXE: marked UNVERIFIED. PDMR web pages redirected and
  said nothing on growth data; CancerModels returned 500. Retry later (suggested
  2026-12-04) or query the PDX Finder data on GitHub/EBI FTP.
- No GEO/ArrayExpress series pairing RNA-seq with a growth table was found in the
  searches run. One Dryad GBM xenograft set (spatial transcriptomics, doi:10.5061/dryad.wpzgmsbv6)
  samples tumours by growth stage but is not a volume time course (UNVERIFIED beyond the
  search snippet).

## 3. NSC migration candidates

| Candidate | What exists | Access | Verdict |
|---|---|---|---|
| Weber 2022 JoVE, doi:10.3791/63102 | Protocol | Not open access; no supplementary data (from prior notes in docs/osf_final_fields.md, UNVERIFIED this session) | No data |
| Rust 2022 J Transl Med, doi:10.1186/s12967-022-03610-5 | MSC/NSC imaging in glioma | "raw data are available on request"; imaging could not resolve local migration (prior notes, UNVERIFIED this session) | Request, low yield |
| Ahn/Aboody-lab J Neurosci 2008 (doi:10.1523/JNEUROSCI.0296-08.2008) | BLI vs intravital microscopy of human NSCs in glioma | Not checked for supplement | Search snippet: BLI "does not offer the spatial resolution to track single-cell migration in vivo" (UNVERIFIED). Not a distance-over-time dataset |
| MRI-tracked iron-labelled hNSCs (FE-Pro, 7 T) and UCSF 14.1 T MPIO MSC/NSC studies | Longitudinal MRI of cell location | Not checked | UNVERIFIED; figures only likely |

No deposited quantitative NSC migration dataset (GEO, Zenodo, figshare, Dryad) was found.
The searches were limited in number, so absence is not proven.

## 4. Recommended actions with dates

| By | Action |
|---|---|
| 2026-11-30 | Lead decides the four-subtype issue (accept 3-class panels by dated amendment, or accept NOT EVALUABLE). |
| 2026-12-01 | Send email A (Mayo) and email B (Alcaniz group). |
| 2026-12-01 | Send email C (NSC data request: Rust 2022 corresponding author, copy to Aboody lab). |
| 2026-12-04 | Retry CancerModels.org API and PDMR portal; check PDMR for GBM models with growth curves. |
| 2026-12-15 | Chase non-replies once. |
| 2027-01-08 | Go/no-go: if no per-time-point control curves from at least 10 models, record NOT EVALUABLE in section 8. |
| 2027-01-15 | Registered deadline. Data entered in section 8 before any calibration. |

If the Mayo resource supplies per-animal volumes, 10 or more curves is likely reachable
(96 models). The subtype clause is the open question.

## 5. Email drafts

Email A, Mayo Clinic Brain Tumor PDX National Resource

> Subject: Request for control-arm growth data, glioblastoma PDX panel
>
> Dear Dr. Sarkaria and colleagues,
>
> I am [YOUR NAME / AFFILIATION]. We are running a pre-registered computational study of
> mechanistic tumour-growth models and need to test calibration against measured PDX
> growth. We read Vaubel et al. (Clin Cancer Res 2020), where all models and genomic data
> are stated to be openly available.
>
> We would like, for each glioblastoma line, the untreated (control) tumour volumes at each
> measurement time (mean per time point, or per animal), with at least three time points per
> line, together with the subtype call per line (Table S6A) and the RNA-seq sample
> identifiers (PRJNA548556). We would use the data only to compute normalised RMSE of a
> growth model and would cite the resource. Our deadline is 2027-01-15. Could you tell us
> whether these data can be shared and what agreement is needed?
>
> Thank you,
> [YOUR NAME / AFFILIATION]

Email B, Alcaniz 2023 authors

> Subject: Request for vehicle-control tumour volumes, Front Oncol 2023 GBM PDX panel
>
> Dear Dr. Alcaniz and co-authors,
>
> I am [YOUR NAME / AFFILIATION]. We are testing a tumour-growth model against PDX data
> and found your panel of 26 subcutaneous GBM PDX models (Front Oncol 2023;13:1129627).
> The supplement holds mutation tables and figures, but not volumes. Could you share the
> vehicle-control tumour volumes per time point (per animal or mean) for each model, and
> the per-model subtype assignment? We will use them only for a calibration check and will
> credit your work. Our deadline is 2027-01-15.
>
> Thank you,
> [YOUR NAME / AFFILIATION]

Email C, NSC in vivo tracking

> Subject: Request for quantitative in vivo NSC/MSC migration data in glioma models
>
> Dear Dr. [corresponding author, Rust et al. 2022] / Dr. Aboody,
>
> I am [YOUR NAME / AFFILIATION]. We are validating a simulation of therapeutic stem-cell
> migration toward glioma and need quantitative in vivo data: per-animal migration distance
> (or the fraction of cells reaching the tumour) at two or more time points, from MRI, BLI
> or histology. Rust et al. 2022 notes that raw data are available on request. Do you have a
> dataset of this kind that you could share, even if small? Our deadline is 2027-01-15, and
> we would acknowledge the source.
>
> Thank you,
> [YOUR NAME / AFFILIATION]

## 6. Sources read

VERIFIED (read this session):
- Europe PMC REST, search results for doi:10.3389/fonc.2023.1129627 and Vaubel 2020 (title, PMC ids).
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10126369/fullTextXML (Alcaniz full text; data availability: "EGAS00001007119"; "RNA sequencing (RNASeq) of untreated s.c. tumor tissue samples was performed for 23 established GBM PDX models").
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10126369/supplementaryFiles (zip: Table_1.xlsx, sheet "Mutations", plus images).
- https://pmc.ncbi.nlm.nih.gov/articles/PMC7056576/ (Vaubel; read via WebFetch, first 100000 characters; subtype counts, PRJNA548556, availability statement).
- https://dctd.cancer.gov/drug-discovery-development/reagents-materials/pdmr (PDMR; "made available upon request").
- https://rdrr.io/github/wfforrest/maeve/man/gdc0084_temo_01.html (dataset description).

REFUSED or FAILED (recorded, not worked around):
- Mayo resource page https://www.mayo.edu/research/labs/translational-neuro-oncology/mayo-clinic-brain-tumor-patient-derived-xenograft-national-resource/pdx-characteristics: HTTP 403.
- https://www.cancermodels.org/api/... : HTTP 500 (two attempts).
- https://ccgdata.cancer.gov/node/2932: HTTP 503.
- Europe PMC fullTextXML for PMC7056576: HTTP 500; supplementary files: error.

UNVERIFIED (web-search snippets only, not opened): PRoXE scope; Aboody J Neurosci 2008 content;
MRI-tracking studies; Dryad spatial-transcriptomics set; Rust 2022 and Weber 2022 (carried
over from docs/osf_final_fields.md).
