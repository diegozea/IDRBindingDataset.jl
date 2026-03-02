# Raw Data Workflow

This page describes how to generate and update files under `data/raw/` using:

- `scripts/01_download_disprot_data.jl`
- `scripts/02_download_caid_data.jl`
- `scripts/03_download_sifts.jl`

It also documents the optional ELM scratch download script for non-redistributable data:

- `scripts/09_download_elm.jl`

Downstream motif mapping from ELM scratch files to cluster proteins is documented in:

- `scripts/10_mapping_elm_motifs.jl`
- `scripts/11_elm_motif_clustering.jl`
- [ELM Motif Workflow](elm_motifs.md)

Run commands from the repository root.

## Prerequisites

- Julia installed.
- Network access to:
  - `https://disprot.org`
  - `https://caid.idpcentral.org`
  - `https://ftp.ebi.ac.uk`
  - `https://elm.eu.org` (with automatic fallback to `http://elm.eu.org` in script `09`)

Both scripts activate and instantiate `scripts/Project.toml` automatically.

## Generate Or Update Raw Data

Run the scripts in this order:

```bash
julia scripts/01_download_disprot_data.jl
julia scripts/02_download_caid_data.jl
julia scripts/03_download_sifts.jl
```

Optional scratch download (ELM instances + motif regex index):

```bash
julia scripts/09_download_elm.jl
```

Optional derived mapping (requires `09` outputs and clustering output):

```bash
julia scripts/10_mapping_elm_motifs.jl
julia scripts/11_elm_motif_clustering.jl
```

## Optional CLI Arguments

DisProt download script (`01`):

```bash
julia scripts/01_download_disprot_data.jl --release 2025_12 --idpo_version 0.3.0
```

- `--release`: DisProt release identifier.
- `--idpo_version`: IDPO version identifier.

CAID download script (`02`):

```bash
julia scripts/02_download_caid_data.jl --raw_data_dir data/raw
```

- `--raw_data_dir`: target root for CAID output folders.

SIFTS download script (`03`):

```bash
julia scripts/03_download_sifts.jl --raw_data_dir data/raw
```

- `--raw_data_dir`: target root for the SIFTS output folder.

ELM scratch download script (`09`):

```bash
julia scripts/09_download_elm.jl --instances_url "https://elm.eu.org/instances.tsv?q=*" --elms_index_url "http://elm.eu.org/elms/elms_index.tsv" --timeout_seconds 180
```

- `--instances_url`: ELM instances endpoint.
- `--elms_index_url`: ELM classes index endpoint (includes motif regexes).
- `--timeout_seconds`: HTTP timeout per attempt.
- Output location is fixed by design to avoid accidental commits and to keep downstream scripts deterministic:
  - `data/scratch/elm/instances.tsv`
  - `data/scratch/elm/elms_index.tsv`

## Expected Output Layout

After a successful run:

```text
data/raw/
|- disprot/
|  |- README.md
|  |- disprot_regions.json.gz
|  |- idpo.json.gz
|  `- mapping_ontology.tsv.gz
|- caid_2/
|  |- README.md
|  `- *.fasta.gz
|- caid_3/
|  |- README.md
|  `- *.fasta.gz
`- sifts/
   |- README.md
   |- pdb_chain_uniprot.csv.gz
   `- pdb_chain_cath_uniprot.csv.gz
```

Optional scratch outputs:

```text
data/scratch/
`- elm/
   |- README.md
   |- instances.tsv
   `- elms_index.tsv
```

## Notes

- `01_download_disprot_data.jl` rewrites DisProt files inside `data/raw/disprot/`.
- `02_download_caid_data.jl` rewrites CAID files inside `data/raw/caid_2/` and `data/raw/caid_3/`.
- `03_download_sifts.jl` rewrites SIFTS files inside `data/raw/sifts/`.
- `09_download_elm.jl` writes ELM scratch files into `data/scratch/elm/` (intended for local, non-redistributable data).
- `10_mapping_elm_motifs.jl` writes scratch mapping files into `data/scratch/elm/mapping/`.
- `11_elm_motif_clustering.jl` writes `data/derived/elm/elm_motif_clustering.csv.gz` and `data/derived/elm/README.md`, plus scratch regex clusters in `data/scratch/elm/clustering/`.
- If a transient network error occurs, rerun the affected script.
