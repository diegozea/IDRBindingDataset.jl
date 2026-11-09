# Pipeline Run Summary

## Schema refresh — 9 November 2026

Re-ran `09_clustering.jl` (4 threads, 4G global clustering memory limit),
`10_fixed_labels.jl`, and `15_final_dataset.jl` against the checked-in inputs.
The clustering, fixed-label, and final tables now have 2, 14, and 18 columns,
respectively. All 3,199 rows, retained values, and row order match the preceding
artifacts after excluding `global_representative_id` and
`local_representative_id`. This preserves 2,439 final clusters and the split
counts of 2,448 training, 440 validation, and 311 testing proteins. CSV and JSON
sequences and residue labels agree.

## Earlier complete pipeline run

The timings below describe the earlier full pipeline run, before the schema
refresh above.

- Total scripts: `13`
- Passed: `13`
- Failed: `0`
- Total elapsed (s): `612.29`

| Script | Status | Exit Code | Elapsed (s) |
|---|---|---:|---:|
| 01_download_disprot_data.jl | PASS | 0 | 47.28 |
| 02_download_caid_data.jl | PASS | 0 | 7.97 |
| 03_download_sifts.jl | PASS | 0 | 5.74 |
| 05_download_invalid_regions.jl | PASS | 0 | 6.52 |
| 06_download_go_data.jl | PASS | 0 | 9.82 |
| 08_pdb_release_dates.jl | PASS | 0 | 66.00 |
| 09_clustering.jl | PASS | 0 | 144.27 |
| 10_fixed_labels.jl | PASS | 0 | 31.47 |
| 11_download_elm.jl | PASS | 0 | 70.51 |
| 13_mapping_elm_motifs.jl | PASS | 0 | 21.99 |
| 14_elm_motif_clustering.jl | PASS | 0 | 162.51 |
| 15_final_dataset.jl | PASS | 0 | 29.50 |
| 16_check_idr_binding.jl | PASS | 0 | 8.69 |
