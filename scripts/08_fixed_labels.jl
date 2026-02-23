#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         AF2_TRAINING_DATA_CUTOFF_DATE,
                         DISPROT_REGIONS_RELATIVE_PATH,
                         _discover_caid_binding_datasets,
                         _require_caid_dataset_path,
                         build_disprot_protein_annotations,
                         build_fixed_labels_table,
                         build_oldest_pdb_release_dates,
                         clean_directory,
                         read_caid_disprot_ids,
                         read_disprot_entries_from_raw_data,
                         write_fixed_labels_csv_gz,
                         write_fixed_labels_readme

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_RAW_DATA_DIR = normpath(joinpath(REPO_ROOT, "data", "raw"))
const DEFAULT_CLUSTER_REPRESENTATIVES_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "clustering", "cluster_representatives.csv.gz"))
const DEFAULT_PDB_RELEASE_MAP_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "pdb", "uniprot_pdb_release_date_map.csv.gz"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "derived", "fixed"))
const DEFAULT_OUTPUT_FILENAME = "cluster_fixed_labels.csv.gz"
const FIXED_README_FILENAME = "README.md"

"""
    parse_command_line(args)

Parse command-line options for fixed-label dataset generation.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Build the fixed cluster labels table (`cluster_fixed_dataset_label`) from:
- clustering representatives,
- DisProt annotations,
- PDB release-date mappings,
- CAID binding FASTA IDs from all `caid_*` folders found under `raw_data_dir`.

Outputs are written to `data/derived/fixed/`:
- compressed CSV with per-protein fixed labels,
- README.md with source provenance and label-count summary.

Run from the repository root:
  julia scripts/08_fixed_labels.jl --raw_data_dir data/raw --cluster_representatives_path data/derived/clustering/cluster_representatives.csv.gz --pdb_release_map_path data/derived/pdb/uniprot_pdb_release_date_map.csv.gz --output_dir data/derived/fixed --output_filename cluster_fixed_labels.csv.gz

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--raw_data_dir"
        help = "Raw data directory containing disprot/ and one or more caid_*/ folders with binding.fasta.gz."
        arg_type = String
        default = DEFAULT_RAW_DATA_DIR
        "--cluster_representatives_path"
        help = "Compressed CSV path for clustering representatives."
        arg_type = String
        default = DEFAULT_CLUSTER_REPRESENTATIVES_PATH
        "--pdb_release_map_path"
        help = "Compressed CSV path for the UniProt-PDB release-date map."
        arg_type = String
        default = DEFAULT_PDB_RELEASE_MAP_PATH
        "--output_dir"
        help = "Output directory for fixed-label artifacts."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--output_filename"
        help = "Output compressed CSV filename written under `output_dir`."
        arg_type = String
        default = DEFAULT_OUTPUT_FILENAME
    end

    parsed = parse_args(args, settings)
    (
        raw_data_dir = normpath(parsed["raw_data_dir"]),
        cluster_representatives_path = normpath(parsed["cluster_representatives_path"]),
        pdb_release_map_path = normpath(parsed["pdb_release_map_path"]),
        output_dir = normpath(parsed["output_dir"]),
        output_filename = String(parsed["output_filename"])
    )
end

"""
    _relpath_from_repo(path)

Return the path relative to repository root.
"""
function _relpath_from_repo(path::AbstractString)::String
    return relpath(normpath(path), REPO_ROOT)
end

function main()
    options = parse_command_line(ARGS)
    output_csv_path = normpath(joinpath(options.output_dir, options.output_filename))

    # Rebuild this derived artifact from scratch on each run.
    clean_directory(options.output_dir)

    cluster_representatives = DataFrame(CSV.File(options.cluster_representatives_path))
    pdb_release_map = DataFrame(CSV.File(options.pdb_release_map_path))
    disprot_entries = read_disprot_entries_from_raw_data(options.raw_data_dir)
    disprot_annotations = build_disprot_protein_annotations(disprot_entries)
    oldest_pdb_release_dates = build_oldest_pdb_release_dates(pdb_release_map)
    caid_datasets = _discover_caid_binding_datasets(options.raw_data_dir)
    _require_caid_dataset_path(caid_datasets, 2)
    caid_non2_binding_fasta_paths = [dataset.binding_fasta_path
                                     for dataset in caid_datasets
                                     if dataset.number != 2]
    isempty(caid_non2_binding_fasta_paths) && throw(
        ArgumentError(
        "No CAID datasets beyond caid_2 were found under $(options.raw_data_dir). At least one additional caid_* dataset is required."
    ),
    )

    caid_disprot_ids_by_dataset = Dict(
        dataset.number => read_caid_disprot_ids(dataset.binding_fasta_path)
    for dataset in caid_datasets
    )

    fixed_labels_table = build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE
    )

    write_fixed_labels_csv_gz(fixed_labels_table, output_csv_path)

    disprot_regions_path = joinpath(options.raw_data_dir, DISPROT_REGIONS_RELATIVE_PATH)
    readme_path = joinpath(options.output_dir, FIXED_README_FILENAME)
    write_fixed_labels_readme(
        readme_path;
        output_csv_path = _relpath_from_repo(output_csv_path),
        output_csv_filesystem_path = output_csv_path,
        source_cluster_representatives_path = _relpath_from_repo(
            options.cluster_representatives_path
        ),
        source_disprot_regions_path = _relpath_from_repo(disprot_regions_path),
        source_pdb_release_map_path = _relpath_from_repo(options.pdb_release_map_path),
        source_caid_binding_fasta_paths_by_version = Dict(
            dataset.number => _relpath_from_repo(dataset.binding_fasta_path)
        for dataset in caid_datasets
        ),
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE,
        fixed_labels_table = fixed_labels_table
    )

    println(
        "Fixed labels written to $(output_csv_path) with $(nrow(fixed_labels_table)) proteins.",
    )
    return nothing
end

main()
