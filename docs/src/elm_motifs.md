# ELM Motif Workflow

This workflow builds a protein-level ELM motif clustering artifact in two
steps:

1. `scripts/13_mapping_elm_motifs.jl`: map ELM instances to protein regex hits
   (scratch, non-redistributable intermediate).
2. `scripts/14_elm_motif_clustering.jl`: cluster regexes with
   `CompariMotif.jl` and map resulting motif clusters to proteins
   (derived artifact).

Run commands from repository root.

## Step 1: Build Scratch Mapping

Inputs:

- `data/scratch/elm/instances.tsv` (from `scripts/11_download_elm.jl`)
- `data/scratch/elm/elms_index.tsv` (from `scripts/11_download_elm.jl`)
- `data/derived/clustering/cluster_representatives.csv.gz`

Outputs (fixed by script design):

- `data/scratch/elm/mapping/motifs.csv.gz`
- `data/scratch/elm/mapping/README.md`

Command:

```bash
julia scripts/13_mapping_elm_motifs.jl
```

Complete example:

```bash
julia scripts/13_mapping_elm_motifs.jl --instances_path data/scratch/elm/instances.tsv --elms_index_path data/scratch/elm/elms_index.tsv --cluster_representatives_path data/derived/clustering/cluster_representatives.csv.gz --exclude_instance_logic "false positive,true negative"
```

## Step 2: Cluster Regexes And Map Proteins

Input:

- `data/scratch/elm/mapping/motifs.csv.gz`

Outputs:

- scratch regex clusters:
  - `data/scratch/elm/clustering/regex.csv.gz`
- derived protein mapping:
  - `data/derived/elm/elm_motif_clustering.csv.gz`
  - `data/derived/elm/README.md`

Command:

```bash
julia scripts/14_elm_motif_clustering.jl
```

Complete example:

```bash
julia scripts/14_elm_motif_clustering.jl --motifs_path data/scratch/elm/mapping/motifs.csv.gz --regex_cluster_output_dir data/scratch/elm/clustering --regex_cluster_output_filename regex.csv.gz --output_dir data/derived/elm --output_filename elm_motif_clustering.csv.gz --normic_threshold 0.5 --matchic_threshold 1.5
```

## Clustering Protocol

1. Build unique regex list from `motifs.csv.gz`.
2. Run all-vs-all regex comparisons with `CompariMotif.jl`.
3. Keep regex similarity edges with:
   - `NormIC >= 0.5`
   - `MatchIC >= 1.5`
4. Build an undirected similarity graph and assign connected-component IDs.
5. Map regex cluster IDs back to proteins.

## Final Output Schema

`data/derived/elm/elm_motif_clustering.csv.gz` columns:

- `protein_id`: UniProt accession.
- `elm_motif_cluster`: integer connected-component ID of regex motif cluster.

A protein may appear in multiple rows when associated with multiple motif
clusters.

## Notes

- Step `10` intentionally writes into `data/scratch/` to avoid redistribution of
  ELM-derived intermediate data.
- `data/derived/elm/README.md` is regenerated on each `11` run and records
  thresholds, source paths, and summary statistics.
