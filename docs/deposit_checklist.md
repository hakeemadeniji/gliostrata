# OSF deposit checklist

Prepared 2026-11-27. These are the steps between the current commit and an OSF deposit
that `configs/registration.yaml` can record. Do them in the order listed. Each item says
who has to do it, and whether it blocks the deposit:

- **Blocks:** the registration would be false or incomplete without it.
- **Should:** a reviewer would reasonably flag it.
- **After:** belongs after the deposit, not before.

Pre-registration sections are cited as "§". "OSF rule n" is a stopping rule in
`docs/osf_final_fields.md`.

## 1. Decisions only the project lead can make

| # | Decision | Status | Where it lands |
| --- | --- | --- | --- |
| 1.1 | **The T2 comparator.** Choose path A, B or C in `docs/t2_comparator_shortlist.md`; C is recommended. Every path also fixes which of our models enters the row. | **Blocks** | §2.2 item 2 (the remaining "TO BE COMPLETED BY AUTHOR"); for path A, or C with a candidate found, also OSF rule 9 and §5 (one added CGGA read) |
| 1.2 | **Name the T1 subtype classifier.** §2.2 item 9(d) says only "the RNA-seq-trained classifier". Decision D6 of the project record chose single-sample nearest centroid (Pearson) over the Verhaak 840-gene signature, and reserved the centroid set `verhaak840_tcga_all` for T1 (SHA-256 `c2b65d8f…9bf7`, `configs/subtype_centroids/CHECKSUMS.json`). None of this is in the pre-registration, and a classifier chosen after the deposit would be a post-hoc choice. | **Blocks** | §2.2 item 9(d), plus an amendment row. The rule in `configs/scoring_rules.yaml` should then cite the centroid checksum. |
| 1.3 | **Calibration and decision curves.** §3.4 and §3.5 say calibration slope, intercept, calibration plots and net benefit are reported. On 2026-11-27 they were decided "not done for this thesis" (`docs/cgga_validation.md` §12.5). §2.3 says they "have not been run" but does not record that decision. | Should | §3.4 and §3.5, or §2.3, with an amendment row |
| 1.4 | **GLASS.** §3.3 says T2's final numbers are reported on the frozen partitions of CGGA *and GLASS*. GLASS has not been downloaded, and `glass_external_v1` is a placeholder. Either keep it, with a date or a NOT EVALUABLE rule like T3's, or drop it by amendment. | Should | §3.3 and §4.5 |
| 1.5 | **Whether the two pre-seal T2 reads count as the registered reads** (project record 5.4). Both are disclosed in §2.3; the decision itself is not recorded. | Should | §2.3 |

## 2. Fixes to the text. I can make these once you confirm.

| # | Fix | Status |
| --- | --- | --- |
| 2.1 | §5, "Models fitted": still says "The sealed external partition `cgga_external_v1` does not exist yet". Missed when the partition was sealed. | **Blocks** (it is false) |
| 2.2 | §5, "Analyses already run": lists four TCGA analyses. Stage 3b added more on real TCGA data: centroids estimated from TCGA labels (v2 training split, and all labelled patients); single-sample subtype calls, including the descriptive 0.737 agreement on 76 v2-test patients; the feature-store build; and cohort-relative state and TME scores over the 283 profiled patients, 82 v2-test patients among them. The section exists "so that no reader has to take on trust" what was done. | **Blocks** |
| 2.3 | After 1.1 to 1.5 and 2.1 to 2.2: bring `docs/osf_final_fields.md` back in line and remove its "Still open" note on the comparator. | **Blocks** |
| 2.4 | Result files are git-ignored but cited. Upload them as OSF attachments: `data/results/cgga/cox_external_validation.json`, `data/results/cgga/brier18_external.json`, `data/results/t1/neftel_mixture_deconvolution.json`. The 2026-11-26 checksums of the first two are in `docs/project_record/06_evidence_reports.md`; confirm they still match before uploading. | Should |

## 3. Repository (you; these steps are public)

| # | Step | Status |
| --- | --- | --- |
| 3.1 | **Create the GitHub repository and push.** There is no remote, and the header's "Code repository" field needs a URL. Only these files under `data/` are tracked: the three partition manifests, the access ledger and three provenance files. `data/raw`, `processed`, `results` and `external` are git-ignored (checked 2026-11-27). | **Blocks** |
| 3.2 | **Before the repository goes public:** `data/partitions/cgga_external_v1.json` lists 218 CGGA patient identifiers, and the ledger names CGGA files. Whether CGGA's terms allow republishing its identifiers has not been checked. Check them, or keep the repository private until checked. The OSF can hold the manifest privately until then. TCGA barcodes are public. | **Blocks** a public push |
| 3.3 | **Let CI run once on GitHub.** None of the jobs has run there. `test-floor` (Python 3.10) tests a claim that has never been tested, and may fail. Registering a commit whose CI fails would be awkward to explain. | Should |
| 3.4 | **Choose the registered commit.** This is the commit that holds the final pre-registration text, after sections 1 and 2, and the code. A document cannot contain its own commit SHA, so the SHA goes on the OSF form and into the post-deposit commit (step 6.1). Tagging it, for example `prereg-deposit`, makes it easy to cite. | **Blocks** |

## 4. Pre-registration header (you)

The header (lines 9 to 28) lists these fields as TO BE COMPLETED BY AUTHOR:

| Field | When |
| --- | --- |
| Authors, ORCIDs, affiliations | Before the registered commit |
| Thesis committee or supervisor | Before the registered commit |
| Funding and conflict-of-interest statement | Before the registered commit |
| Planned start date of the confirmatory analysis phase | Before the registered commit |
| Code repository URL | Before the registered commit (after 3.1) |
| Confirm §5 is still accurate (the last "TO BE COMPLETED" line of §5) | At the moment of deposit |
| Registered commit SHA | After the commit exists; goes on the OSF form and into step 6.1 |
| OSF DOI; date of deposit | After the deposit (step 6.1) |

## 5. The OSF form (you)

1. Re-paste the four Sampling fields from `docs/osf_final_fields.md`, after step 2.3. The
   copy on the form predates the latest corrections.
2. Check that the other OSF sections (hypotheses, design, analysis plan) match the
   pre-registration at the registered commit. `docs/osf_final_fields.md` covers Sampling
   only.
3. Upload the attachments:
   - `docs/preregistration.md` (the registered commit's copy)
   - `src/gliostrata/evaluation/thresholds.py`
   - `data/partitions/tcga_gbm_clinical_v1.json`, `tcga_gbm_clinical_v2.json` and
     `cgga_external_v1.json` (the last one subject to 3.2)
   - `docs/cgga_validation.md`
   - `docs/t1_deconvolution_spec.md`
   - the result files of 2.4
   - from path C: the comparator search log
4. Deposit. Note the DOI and the date.

## 6. Recording the deposit, in one commit (I can do this once you send the DOI)

| # | Step |
| --- | --- |
| 6.1 | Fill `configs/registration.yaml` (`osf_registration_doi`, `deposited_at`, `registered_commit`) and the pre-registration header **in the same commit**. `tests/test_read_guard.py::test_registration_config_and_preregistration_header_agree` checks that the two agree. Only a DOI-shaped value opens the read guard's registration gate, so never enter a placeholder. |
| 6.2 | Add an amendment row in §8 for the header fill (no rule changes), `CHANGELOG.md`, and project record 05 §5.1 step 7. |
| 6.3 | Run the full suite (sockets blocked, CGGA audit hook) and the threshold checksum test. |

## 7. After the deposit, not before

- **T1 subtype read.** Needs the CGGA microarray matrix (to be downloaded manually and
  stored unopened), a frozen scoring-set manifest for the 108 labelled patients, and that
  manifest's checksum in `configs/scoring_rules.yaml`. Until then the guard refuses the
  read.
- **T2 improvement row**, if path A, or C with a candidate, was chosen: one read, under its
  own rule id.
- **T4 and T5.** Need a real simulator and response surface. Both are stubs today.
- **Not needed for the deposit:** the 3CA confirmation email (optional, see
  `data/provenance/neftel2019_sources.json`), DVC, Snakemake, the container image, the
  stale 155-patient store, and `openpyxl`/`pyreadr` in the lock file. These are recorded as
  known gaps in `docs/project_record/05_future_plan.md` §5.5.
