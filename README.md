# IDRBindingDataset

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/IDRBindingDataset.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://diegozea.github.io/IDRBindingDataset.jl/dev/)
[![Build Status](https://github.com/diegozea/IDRBindingDataset.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/diegozea/IDRBindingDataset.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

IDRBindingDataset.jl builds a reproducible, FAIR-oriented per-protein dataset for
intrinsically disordered region (IDR) binding research, with deterministic
training/validation/testing assignment and residue-level labels.

The repository is designed so that the full dataset can be regenerated from
raw sources with explicit provenance and auditable transformation criteria.

## Overview

This repository builds a protein-level dataset for IDR-binding analyses starting
from DisProt proteins and their annotated regions. The final dataset contains
one row per DisProt protein accession, keeps the full amino-acid `sequence`, and
stores a per-residue label track in `residue_labels`. These two fields are
positionally aligned, so residue label `i` describes sequence position `i` in
the same protein.

Protein sequences are grouped before any split is assigned. The pipeline first
applies a global DIAMOND clustering/reclustering pass to capture broader
sequence similarity, then a local all-vs-all DIAMOND blastp-like pass to detect
short high-identity matches. Final clusters are defined by the transitive
closure across both passes, so proteins connected by either rule stay together.
Train/validation/test labels are then assigned per final cluster using CAID
reference-protein membership and PDB release dates relative to the
2018-04-30 AlphaFold 2 training-data cutoff. Residue labels are derived from
retained DisProt regions through GO and IDPO ontology terms, CAID-3 invalid
regions are excluded before labeling, and optional ELM motif-cluster
annotations can be added for proteins with mapped motifs.

## Final Dataset Artifacts

The final dataset is generated under `data/dataset/`:

- `data/dataset/cluster_dataset_labels.csv.gz`
- `data/dataset/cluster_dataset_labels.json.gz`

Both artifacts represent the same proteins and split labels, with one final
record per DisProt protein.
`residue_labels` are encoded as:

- JSON: integer vectors
- CSV: compact strings (`0`, `1`, `2`)

In both formats, `sequence` and `residue_labels` are positionally aligned and
have the same length for a given protein.

## End-to-End Pipeline

The default orchestrated pipeline is implemented in `scripts/run_pipeline.jl`
and runs these scripts in order:

1. `01_download_disprot_data.jl`
2. `02_download_caid_data.jl`
3. `03_download_sifts.jl`
4. `05_download_invalid_regions.jl`
5. `06_download_go_data.jl`
6. `08_pdb_release_dates.jl`
7. `09_clustering.jl`
8. `10_fixed_labels.jl`
9. `11_download_elm.jl`
10. `13_mapping_elm_motifs.jl`
11. `14_elm_motif_clustering.jl`
12. `15_final_dataset.jl`
13. `16_check_idr_binding.jl`

Run from repository root:

```bash
julia scripts/run_pipeline.jl
```

## Deterministic Criteria Used To Build The Dataset

### 1. Sequence-Clustering Criteria

Protein-level clustering is built in `src/sequence_clustering.jl`. The pipeline
first runs a global DIAMOND clustering/reclustering pass, then a local
all-vs-all DIAMOND blastp-like pass, and finally defines `final_cluster_id` by
taking the transitive closure across proteins that share either a global or a
local representative.

- Global DIAMOND reclustering: identity >= 30.0 percent and mutual coverage >=
  80.0 percent
- Local all-vs-all DIAMOND blastp-like clustering: identity >= 80.0 percent on
  ungapped alignments with length >= 10 amino acids

These thresholds are documented in:

- `data/derived/clustering/README.md`

The saved clustering table contains only `protein_id` and `final_cluster_id`.
Global and local representative assignments are used internally to calculate
the final clusters; `global_representative_id` and `local_representative_id`
are no longer exported in intermediate or final dataset tables. The existing
clustering function names and `cluster_representatives.csv.gz` filename are
retained. Downstream table builders also omit these two columns when reading
older intermediate tables. CAID experiments retain their separate protein-level
and cluster-level membership columns.

### 2. Split Criteria (Leakage-Aware)

Fixed labels are assigned at cluster level in `src/fixed_labels.jl`.

Fixed cluster label rules:

- `caid` if a cluster contains at least one CAID reference protein
- `training` if a cluster has at least one structure dated on or before
  2018-04-30 (AlphaFold 2.3 training cutoff)
- `training` (conservative) if a cluster has PDB structure mappings but missing
  release dates
- `missing` otherwise

Final split mapping in `src/dataset_finalization.jl`:

- `training` -> `training`
- `caid` -> `testing`
- `missing` -> `validation` (default)

Why this helps leakage control:

- Split assignment is cluster-based rather than independent per protein, so
  proteins grouped by the sequence-clustering rules stay in the same split.

### 3. Residue-Level Label Provenance

Residue labels are generated in `src/dataset_json_export.jl` using ontology
rules defined in `src/disprot_challenge_terms.jl`. For each protein,
`residue_labels` has the same length as `sequence`, and each label refers to
the residue at the same sequence position.

Label values and provenance:

- `0`: default ordered/non-annotated class in the final dataset; the residue
  was not labeled as disorder or binding after all retained annotations were
  applied
- `1`: residue covered by a retained DisProt region whose term is a descendant
  of `IDPO:0000002`
- `2`: residue covered by a retained DisProt region whose term is a descendant
  of `GO:0005488`

Ontology roots used for label classification:

- Disorder: IDPO branch rooted at `IDPO:0000002`
- Binding: GO branch rooted at `GO:0005488` (with explicit handling for key
  descendants including `GO:0005515` and `GO:0003676`)

Additional labeling rules:

- CAID-3 invalid regions are excluded during labeling
- Binding has precedence over disorder on overlapping regions
  (monotonic update: 0 -> 1 -> 2)
- CAID-3 validation checks the mapping `- -> 0`, `0 -> 1`, `1 -> 2`

### 4. ELM Motif Criteria

ELM motif mapping and clustering are generated by:

- `13_mapping_elm_motifs.jl`
- `14_elm_motif_clustering.jl`

Key criteria:

- Exclude ELM instance logic values: `false positive`, `true negative`
- Regex similarity edges require:
  - `NormIC >= 0.5`
  - `MatchIC >= 1.5`
  - at least 2 shared positions

Summary and thresholds are documented in:

- `data/derived/elm/README.md`

## Final Dataset Labels (What They Mean)

### Split Label (`cluster_dataset_label`)

- `training`: cluster assigned to training according to fixed-label protocol
- `validation`: unresolved fixed label (`missing`) mapped to validation
- `testing`: cluster containing at least one CAID reference protein

### Residue Label (`residue_labels`)

`sequence` and `residue_labels` are aligned position by position.

- `0`: default ordered/non-annotated class; residue not labeled as disorder or
  binding in the final dataset
- `1`: residue in a retained DisProt region with a qualifying IDPO disorder
  term
- `2`: residue in a retained DisProt region with a qualifying GO binding term
- If disorder and binding overlap, binding takes precedence.

## Snapshot Statistics (Current Checked-In Artifacts)

The values below summarize the currently checked-in dataset artifacts in this
workspace snapshot (generated on 2026-04-21 UTC).

### Core Dataset Size

| Metric | Value | Source |
|---|---:|---|
| Proteins in final dataset | 3199 | `data/dataset/cluster_dataset_labels.csv.gz` |
| Final clusters | 2439 | `data/dataset/cluster_dataset_labels.csv.gz` |

### Split Distribution

| Split | Proteins | Final clusters |
|---|---:|---:|
| Training | 2448 | 1936 |
| Validation | 440 | 400 |
| Testing | 311 | 103 |

### IDR-Binding Annotation Coverage

| Split | IDR-binding proteins | IDR-binding clusters |
|---|---:|---:|
| Training | 640 | 549 |
| Validation | 137 | 128 |
| Testing | 182 | 97 |
| Total | 959 | 774 |

### ELM Mapping Summary

| Metric | Value | Source |
|---|---:|---|
| Protein-cluster assignment rows | 462 | `data/derived/elm/README.md` |
| Proteins with >=1 ELM cluster | 442 | `data/derived/elm/README.md` |
| Unique motif clusters | 10 | `data/derived/elm/README.md` |

### CAID-3 Residue-Label Validation

From `data/dataset/CAID3_LABELS.md`:

- CAID-3 raw labels are checked with the mapping `- -> 0`, `0 -> 1`, `1 -> 2`
- CAID-3 proteins checked: 52
- CAID-3 proteins passing all checks: 52
- Residues compared: 28263
- Residue mismatches: 0
- Overall result: PASS

### Full Pipeline Status

From `data/dataset/PIPELINE_RUN_SUMMARY.md`:

- Scripts run: 13
- Passed: 13
- Failed: 0
- Total elapsed: 612.29 seconds

## Data Schema (Key Columns)

Final records preserve cluster/fixed-label metadata and add sequence/label
fields.

| Column | Meaning |
|---|---|
| `protein_id` | UniProt accession |
| `final_cluster_id` | Deterministic sequence-cluster identifier |
| `cluster_dataset_label` | Final split (`training`, `validation`, `testing`) |
| `protein_elm_clusters` | Joined ELM motif-cluster IDs (or missing) |
| `sequence` | Full protein sequence (uppercase) |
| `residue_labels` | Per-residue labels (`0`, `1`, `2`) |

Extended column details for final artifacts are documented in:

- `data/dataset/README.md`

## Reproducibility

Run from repository root.

Generate final dataset artifacts only:

```bash
julia scripts/15_final_dataset.jl
```

Validate CAID-3 residue-label consistency:

```bash
julia scripts/16_check_idr_binding.jl
```

Full pipeline:

```bash
julia scripts/run_pipeline.jl
```

CLI help for modified parameters:

```bash
julia scripts/15_final_dataset.jl --help
julia scripts/16_check_idr_binding.jl --help
julia scripts/run_pipeline.jl --help
```

## External Data Sources And Ontologies

Exact run-specific download URLs and timestamps for the current checked-in
snapshot are recorded in generated provenance files under `data/raw/` and, for
ELM downloads, in `data/scratch/elm/README.md`. The main provenance files are:

- `data/raw/disprot/README.md`
- `data/raw/caid_2/README.md`
- `data/raw/caid_3/README.md`
- `data/raw/go/README.md`
- `data/raw/sifts/README.md`
- `data/scratch/elm/README.md`

General upstream resources represented in those provenance files include:

- DisProt API (`disprot_regions.json`) and DisProt ontology mapping files
- IDPO release assets
- CAID reference sets (including CAID-2 and CAID-3)
- Gene Ontology `go-basic.obo`
- SIFTS chain-level mapping data
- ELM instances and ELM motif index

General upstream entry points include:

- https://disprot.org/
- https://caid.idpcentral.org/
- https://purl.obolibrary.org/obo/go/go-basic.obo
- https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
- https://elm.eu.org/

## License

This repository is licensed under MIT. See `LICENSE`.
