#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset:
                         build_pdb_release_map_from_raw,
                         clean_directory,
                         write_pdb_release_map_csv_gz,
                         write_pdb_release_readme

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))
const DERIVED_PDB_DIR = normpath(joinpath(@__DIR__, "..", "data", "derived", "pdb"))
const DEFAULT_OUTPUT_FILENAME = "uniprot_pdb_release_date_map.csv.gz"
const DEFAULT_OUTPUT_PATH = joinpath(DERIVED_PDB_DIR, DEFAULT_OUTPUT_FILENAME)
const DEFAULT_BATCH_SIZE = 300

"""
    parse_command_line(args)

Parse command-line arguments for PDB release-date map generation.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Build a derived UniProt-PDB-release-date map from raw DisProt+SIFTS inputs.

Run from the repository root:
  julia scripts/05_pdb_release_dates.jl --raw_data_dir data/raw --output_path data/derived/pdb/uniprot_pdb_release_date_map.csv.gz --batch_size 300

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--raw_data_dir"
        help = "Raw data directory containing disprot/ and sifts/ folders."
        arg_type = String
        default = RAW_DATA_DIR
        "--output_path"
        help = "Output compressed CSV path for the derived UniProt-PDB-release-date map."
        arg_type = String
        default = DEFAULT_OUTPUT_PATH
        "--batch_size"
        help = "Batch size for RCSB GraphQL release-date queries."
        arg_type = Int
        default = DEFAULT_BATCH_SIZE
    end

    parsed = parse_args(args, settings)
    return (
        raw_data_dir = normpath(parsed["raw_data_dir"]),
        output_path = normpath(parsed["output_path"]),
        batch_size = parsed["batch_size"]
    )
end

function main()
    options = parse_command_line(ARGS)

    # Rebuild this derived artifact from scratch on each run.
    clean_directory(DERIVED_PDB_DIR)

    tables = build_pdb_release_map_from_raw(
        options.raw_data_dir;
        batch_size = options.batch_size
    )

    write_pdb_release_map_csv_gz(tables.disprot_pdb_map_with_release, options.output_path)
    readme_path = joinpath(dirname(options.output_path), "README.md")
    write_pdb_release_readme(
        readme_path;
        output_csv_path = options.output_path,
        disprot_pdb_map_with_release = tables.disprot_pdb_map_with_release,
        failed_pdb_lookups = tables.failed_pdb_lookups,
        invalid_merged_pdb_references = tables.invalid_merged_pdb_references,
        batch_size = options.batch_size
    )

    return nothing
end

main()
