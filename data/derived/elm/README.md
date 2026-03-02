# ELM Motif Clustering

This folder stores protein-level assignments to ELM motif clusters derived from
regex similarity graph connected components.

## Artifact

- Compressed CSV: `data/derived/elm/elm_motif_clustering.csv.gz` (`1562` bytes)
- Rows (protein-cluster assignments): `462`
- Unique proteins assigned to at least one motif cluster: `442`
- Unique motif clusters assigned to proteins: `10`

## Data Dictionary

- `protein_id`: UniProt accession.
- `elm_motif_cluster`: connected-component cluster number derived from regex
  similarity graph.

## Clustering Protocol

- Build unique regex list from `motifs.csv.gz`.
- Compare regexes pairwise with `CompariMotif.jl` requiring at least
  `2` matched positions.
- Keep similarity edges with `matched=true`, `NormIC >= 0.5`, and
  `MatchIC >= 1.5`.
- Build an undirected graph and assign connected-component IDs.

## Regex Cluster Summary

- Unique regexes clustered: `233`
- Connected components: `10`
- Singleton components: `9`
- Maximum regexes in one cluster: `224`
- Mean regexes per cluster: `23.3`

## Protein Mapping Summary

- Unique proteins mapped: `442`
- Proteins mapped to multiple motif clusters: `18`
- Mean clusters per protein: `1.045`
- Maximum clusters in one protein: `3`

## Source Inputs

- Protein-regex motifs: `data/scratch/elm/mapping/motifs.csv.gz`
- Regex clusters: `data/scratch/elm/clustering/regex.csv.gz`

## Reproducibility

Run from repository root:

```bash
julia scripts/11_elm_motif_clustering.jl
```

Complete example:

```bash
julia scripts/11_elm_motif_clustering.jl --motifs_path data/scratch/elm/mapping/motifs.csv.gz --regex_cluster_output_dir data/scratch/elm/clustering --regex_cluster_output_filename regex.csv.gz --output_dir data/derived/elm --output_filename elm_motif_clustering.csv.gz --normic_threshold 0.5 --matchic_threshold 1.5
```
