# Raw Data Workflow

This page describes how to generate and update files under `data/raw/` using:

- `scripts/01_download_disprot_data.jl`
- `scripts/02_download_caid_data.jl`
- `scripts/03_download_sifts.jl`
- `scripts/14_download_invalid_regions.jl`
- `scripts/15_download_go_data.jl`

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
  - `https://raw.githubusercontent.com` (CAID-3 invalid-region TSV)
  - `https://purl.obolibrary.org` (GO `go-basic.obo`)
  - `https://elm.eu.org` (with automatic fallback to `http://elm.eu.org` in script `09`)

Both scripts activate and instantiate `scripts/Project.toml` automatically.

## Generate Or Update Raw Data

Run the scripts in this order:

```bash
julia scripts/01_download_disprot_data.jl
julia scripts/02_download_caid_data.jl
julia scripts/03_download_sifts.jl
julia scripts/14_download_invalid_regions.jl
julia scripts/15_download_go_data.jl
```

## Full Pipeline Run

To run the default end-to-end numbered pipeline (including dataset generation
and CAID-3 label validation):

```bash
julia scripts/run_pipeline.jl
```

The final validation stage (`scripts/16_check_idr_binding.jl`) is blocking by
default and fails the pipeline when CAID-3 vs dataset labels mismatch.
Use `--allow_mismatch` only for report-only diagnostic runs.

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
julia scripts/01_download_disprot_data.jl --release 2025_12
```

- `--release`: DisProt release identifier.
- IDPO is always downloaded from the latest `BioComputingUP/idpo` GitHub release asset (`idpo.json`).

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

CAID-3 invalid-region download script (`14`):

```bash
julia scripts/14_download_invalid_regions.jl --output_dir data/raw/caid_3 --filename invalid_regions_fullCAID3.tsv --source_url https://raw.githubusercontent.com/BioComputingUP/caid-reference/main/data/disprot/invalid_regions_fullCAID3.tsv
```

- `--output_dir`: target directory for the invalid-region TSV and CAID-3 README updates.
- `--filename`: output filename for the TSV.
- `--source_url`: source URL of the TSV.

GO download script (`15`):

```bash
julia scripts/15_download_go_data.jl --raw_data_dir data/raw
```

- `--raw_data_dir`: target root for `data/raw/go/`.

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
|  |- invalid_regions_fullCAID3.tsv
|  `- *.fasta.gz
|- go/
|  |- README.md
|  `- go-basic.obo.gz
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
- `14_download_invalid_regions.jl` refreshes `data/raw/caid_3/invalid_regions_fullCAID3.tsv`.
- `15_download_go_data.jl` rewrites GO files inside `data/raw/go/`.
- `09_download_elm.jl` writes ELM scratch files into `data/scratch/elm/` (intended for local, non-redistributable data).
- `10_mapping_elm_motifs.jl` writes scratch mapping files into `data/scratch/elm/mapping/`.
- `11_elm_motif_clustering.jl` writes `data/derived/elm/elm_motif_clustering.csv.gz` and `data/derived/elm/README.md`, plus scratch regex clusters in `data/scratch/elm/clustering/`.
- If a transient network error occurs, rerun the affected script.
