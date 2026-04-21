# Fixed Dataset Labels

This folder contains the per-protein fixed cluster labels used as the
deterministic base split:
- `training`: cluster has at least one PDB release date on or before the AF2
  training cutoff (`2018-04-30`), or has PDB mappings with
  missing release dates (conservative assignment).
- `caid`: cluster contains at least one protein in the CAID reference sets
  used for this build.
- `missing`: cluster does not satisfy either rule above.

## Artifact

- Compressed CSV: `data/derived/fixed/cluster_fixed_labels.csv.gz` (`55567` bytes)
- Number of proteins: `3199`
- Number of final clusters: `2439`

## Counts by Fixed Label

### Proteins
- label `training`: `2448`
- label `caid`: `311`
- label `missing`: `440`

### Final clusters
- label `training`: `1936`
- label `caid`: `103`
- label `missing`: `400`

## Source Inputs

- Cluster representatives: `data/derived/clustering/cluster_representatives.csv.gz`
- DisProt regions JSON: `data/raw/disprot/disprot_regions.json.gz`
- PDB release-date map: `data/derived/pdb/uniprot_pdb_release_date_map.csv.gz`
- `source_caid2_binding_fasta_path`: `data/raw/caid_2/binding.fasta.gz`
- `source_caid3_binding_fasta_path`: `data/raw/caid_3/binding.fasta.gz`
