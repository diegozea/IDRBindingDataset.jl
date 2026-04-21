# ELM Motif Mapping

This folder stores regex-level ELM motif mappings for proteins in the
cluster-representative table.

## Artifact

- Compressed CSV: `data/scratch/elm/mapping/motifs.csv.gz` (`11088` bytes)
- Row count: `950`
- Unique proteins with at least one mapped motif: `442`
- Unique regex patterns represented: `233`

## Data Dictionary

- `protein_id`: UniProt accession from cluster representatives.
- `regex`: ELM regex pattern (from `elms_index.tsv`) mapped through `ELMIdentifier`.
- `start`: ELM instance start coordinate (1-based, inclusive).
- `end`: ELM instance end coordinate (1-based, inclusive).

## Mapping Decisions

- `InstanceLogic` filtering excludes: `false positive`, `true negative`.
- Protein mapping first uses `Primary_Acc` and then falls back to parsed tokens
  from `Accessions` when needed.
- Only proteins present in `cluster_representatives.csv.gz` are retained.
- Regex assignments are joined by `ELMIdentifier` from `elms_index.tsv`.
- Duplicate output rows by `(protein_id, regex, start, end)` are removed.

## Summary Statistics

- Cluster proteins evaluated: `3199`
- Proteins with motifs: `442` (`13.82%`)
- Proteins without motifs: `2757`
- Average motifs per motif-positive protein: `2.149`
- Median motifs per motif-positive protein: `1.0`
- Maximum motifs for one protein: `17`
- Average motifs per cluster protein (including proteins without motifs): `0.297`

## Input Filtering

- ELM instance rows before filtering: `4277`
- ELM instance rows after filtering: `4171`
- ELM instance rows removed by `InstanceLogic`: `106`

## Source Inputs

- Cluster representatives: `data/derived/clustering/cluster_representatives.csv.gz`
- ELM instances: `data/scratch/elm/instances.tsv`
- ELM classes index: `data/scratch/elm/elms_index.tsv`

## Reproducibility

Run from repository root:

```bash
julia scripts/13_mapping_elm_motifs.jl
```

Complete example:

```bash
julia scripts/13_mapping_elm_motifs.jl --instances_path data/scratch/elm/instances.tsv --elms_index_path data/scratch/elm/elms_index.tsv --cluster_representatives_path data/derived/clustering/cluster_representatives.csv.gz --exclude_instance_logic "false positive,true negative"
```
