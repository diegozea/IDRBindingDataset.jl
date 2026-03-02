# Final Dataset

This folder stores the final per-protein dataset, adding ELM motif-cluster
assignments and replacing missing fixed split labels with `validation`.

## Artifact

- Compressed CSV: `data/dataset/cluster_dataset_labels.csv.gz` (`56953` bytes; `0.054` MiB)
- Generated at (UTC): `2026-03-02 16:51:32 UTC`
- Number of proteins: `3199`
- Number of final clusters: `2439`
- Proteins with at least one ELM cluster assignment: `442`
- Unique ELM clusters assigned across all proteins: `10`

## Data Dictionary

- `protein_elm_clusters`: one or more ELM motif-cluster IDs for each protein, joined by `-` when multiple clusters are present.
- `cluster_dataset_label`: final split label constrained to `training`,
  `validation`, and `testing`; values are derived from
  `cluster_fixed_dataset_label` by keeping `training`, mapping `caid` to
  `testing`, and assigning missing values to `validation`.

## Counts by Dataset Label

- `training`: proteins `2448`, final clusters `1936`, proteins with ELM clusters `350`, unique ELM clusters `9`
- `validation`: proteins `440`, final clusters `400`, proteins with ELM clusters `25`, unique ELM clusters `2`
- `testing`: proteins `311`, final clusters `103`, proteins with ELM clusters `67`, unique ELM clusters `5`

## Source Inputs

- Fixed labels table: `data/derived/fixed/cluster_fixed_labels.csv.gz`
- Protein ELM motif clustering table: `data/derived/elm/elm_motif_clustering.csv.gz`

## Reproducibility

Run from repository root:

```bash
julia scripts/13_assing_missing_labels.jl
```

Complete example:

```bash
julia scripts/13_assing_missing_labels.jl --fixed_labels_path data/derived/fixed/cluster_fixed_labels.csv.gz --elm_motif_clustering_path data/derived/elm/elm_motif_clustering.csv.gz --output_dir data/dataset --output_filename cluster_dataset_labels.csv.gz --missing_label validation
```

## JSON Artifact

- Compressed JSON: `data/dataset/cluster_dataset_labels.json.gz` (`1293876` bytes; `1.234` MiB)
- Generated at (UTC): `2026-03-02 22:21:20 UTC`
- Number of proteins: `3199`

### Content

- Preserves all columns from `cluster_dataset_labels.csv.gz`.
- Adds `sequence`: full protein sequence from DisProt.
- Adds `residue_labels`: per-residue integer vector with `0 = ordered`,
  `1 = disordered`, and `2 = binding`.
- Label precedence is monotonic (`0 -> 1 -> 2`), so residues labeled as
  `2` are never downgraded by later overlapping regions.

### Source Inputs

- Final dataset CSV: `data/dataset/cluster_dataset_labels.csv.gz`
- DisProt regions JSON: `data/raw/disprot/disprot_regions.json.gz`

### Reproducibility

Run from repository root:

```bash
julia scripts/14_json_format.jl
```

Complete example:

```bash
julia scripts/14_json_format.jl --dataset_csv_path data/dataset/cluster_dataset_labels.csv.gz --raw_data_dir data/raw --output_dir data/dataset --output_filename cluster_dataset_labels.json.gz
```

