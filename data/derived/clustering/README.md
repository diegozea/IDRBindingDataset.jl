# Sequence Clustering Outputs

Generated at (UTC): 2026-11-09T17:17:09Z

## Outputs
- DisProt full-protein FASTA: data/derived/clustering/disprot_full_proteins.fasta.gz (1082833 bytes, 1.0828 MB)
- Cluster representatives CSV: data/derived/clustering/cluster_representatives.csv.gz (18845 bytes, 0.0188 MB)

The CSV contains one row per protein with only `protein_id` and `final_cluster_id`.
Global and local representative assignments are used internally to compute final clusters.

## Sequence and Cluster Counts
- Protein sequence count: 3199
- Global cluster count: 2655
- Local cluster count: 2830
- Final cluster count: 2439

## Clustering Methods
- Global method: diamond_cluster_recluster
- Global thresholds: identity >= 30.0% and mutual coverage >= 80.0%
- Local method: all_vs_all_blastp
- Local thresholds: identity >= 80.0% on ungapped alignments with length >= 10 aa

## Source Data
- DisProt regions JSON: data/raw/disprot/disprot_regions.json.gz (5158321 bytes, 5.1583 MB)
