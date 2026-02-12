# Raw Data Workflow

This page describes how to generate and update files under `data/raw/` using:

- `scripts/01_download_disprot_data.jl`
- `scripts/02_download_caid_data.jl`

Run commands from the repository root.

## Prerequisites

- Julia installed.
- Network access to:
  - `https://disprot.org`
  - `https://caid.idpcentral.org`

Both scripts activate and instantiate `scripts/Project.toml` automatically.

## Generate Or Update Raw Data

Run the scripts in this order:

```bash
julia scripts/01_download_disprot_data.jl
julia scripts/02_download_caid_data.jl
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
`- caid_3/
   |- README.md
   `- *.fasta.gz
```

## Notes

- `01_download_disprot_data.jl` rewrites DisProt files inside `data/raw/disprot/`.
- `02_download_caid_data.jl` rewrites CAID files inside `data/raw/caid_2/` and `data/raw/caid_3/`.
- If a transient network error occurs, rerun the affected script.
