# KPS = 0 at source: three TCGA patients

Facts only. Checked 2026-11-27 against `data/raw/tcga_gbm/clinical_raw.parquet` (617 rows) and
`data/raw/tcga_gbm/clinical_cohort.parquet` (594 rows), reading
`scripts/download_tcga_clinical.py` for the mapping. TCGA only; no CGGA path was touched, no
data or code was changed, nothing was downloaded.

## The three patients

All three are in the v2 training split (`data/partitions/tcga_gbm_clinical_v2.json`).

| patient_id | raw `kps_raw` | cohort `kps` | vital_status_raw | days_to_death | days_to_follow_up | age_at_index_raw | sex | primary_diagnosis | treatment_types |
|---|---|---|---|---|---|---|---|---|---|
| TCGA-74-6573 | `'0'` | 0.0 | Dead | 105.0 | 105.0 | 67.0 | male | Glioblastoma | Chemotherapy; Radiation, External Beam |
| TCGA-76-6283 | `'0'` | 0.0 | Dead | 165.0 | 165.0 | 71.0 | female | Glioblastoma | none recorded (`[]`) |
| TCGA-76-6660 | `'0'` | 0.0 | Dead | 114.0 | 114.0 | 73.0 | male | Glioblastoma | none recorded (`[]`) |

Cohort `overall_survival_months` for them: 3.45, 5.42, 3.75. Case ids (GDC): TCGA-74-6573
8d0d7a7e-8fae-437d-9b39-f2f17d0f5e4f; TCGA-76-6283 e5c9b14e-3bc7-4892-bd51-4cd1599819ec;
TCGA-76-6660 47fd9c51-298c-41ad-9996-135c8a03d691.

## What the files show

1. **The zero is in the raw download, not made by the mapping.** `kps_raw` is a string column.
   Its values and counts: `'80'` 243, missing 171 (null), `'60'` 86, `'100'` 61, `'90'` 16,
   `'40'` 16, `'70'` 13, `'Unknown'` 4, `'0'` 3, `'20'` 2, `'50'` 2. The three zeros are the
   literal string `'0'`, not a null.
2. **A missing KPS is not coerced to 0 by the code.** `flatten_case` keeps a follow-up value only
   if it is not None, and leaves `kps` as None otherwise. `build_cohort` applies
   `pd.to_numeric(raw["kps_raw"], errors="coerce")`, which turns `'Unknown'` (4 patients) and
   nulls into NaN, never into 0. Cohort `kps` has 152 NaN (171 null + 4 `'Unknown'`, less the
   rows the cohort filter dropped) and 3 zeros, which are exactly the three `'0'` strings.
3. **It is the KPS field itself**, `follow_ups.karnofsky_performance_status`; no other field is
   involved.
4. **No duplicates**: no `submitter_id` repeats in the raw table.
5. **All three died**, at 105, 165 and 114 days. For each, `days_to_follow_up` (the maximum over
   follow-ups, per `flatten_case`) equals `days_to_death`, so each has a follow-up record dated on
   the day of death.

## Not determinable from these files

- `flatten_case` stores only the first non-null KPS across a case's follow-ups, and the GDC
  follow-up list is unordered (the script's own comment). The raw parquet does not keep the other
  follow-ups, so whether these cases also carry a non-zero KPS on another follow-up cannot be
  seen here. It would need a re-query of the GDC, which was not done.
- Item 5 above is consistent with a KPS of 0 recorded on a terminal (death-date) follow-up, 0
  being the Karnofsky value for death, rather than a performance status at baseline. The files do
  not state which follow-up the value came from, so this is an inference, not a finding.

## Consequence to be aware of (not acted on)

A KPS of 0 sits outside the range of the other values (lowest non-zero is 20). Any analysis that
takes KPS as a complete-case covariate (decision D4 in `05_future_plan.md`) includes these three
as 0 in the training split. Nothing was changed.
