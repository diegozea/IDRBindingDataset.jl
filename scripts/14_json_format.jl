#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset: DISPROT_REGIONS_RELATIVE_PATH,
                         build_dataset_json_records,
                         read_disprot_entries_from_raw_data,
                         write_dataset_json_gz,
                         write_dataset_json_readme_section

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_CSV_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "dataset",
    "cluster_dataset_labels.csv.gz"
))
const DEFAULT_RAW_DATA_DIR = normpath(joinpath(REPO_ROOT, "data", "raw"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "dataset"))
const DEFAULT_OUTPUT_FILENAME = "cluster_dataset_labels.json.gz"
const FINAL_DATASET_README_FILENAME = "README.md"

"""
    _relpath_from_repo(path)

Return one normalized path relative to repository root.
"""
function _relpath_from_repo(path::AbstractString)::String
    relpath(normpath(path), REPO_ROOT)
end

"""
    parse_command_line(args)

Parse command-line options for final dataset JSON export.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Create a JSON representation of `cluster_dataset_labels.csv.gz` where each
protein record keeps all original columns and adds:
- `sequence`: full protein sequence from DisProt.
- `residue_labels`: per-residue integer vector:
  - 0 = ordered residue
  - 1 = disordered residue
  - 2 = binding residue (binding overrides disorder)

Run from repository root:
  julia scripts/14_json_format.jl

Complete example:
  julia scripts/14_json_format.jl --dataset_csv_path data/dataset/cluster_dataset_labels.csv.gz --raw_data_dir data/raw --output_dir data/dataset --output_filename cluster_dataset_labels.json.gz

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_csv_path"
        help = "Input compressed CSV with final per-protein dataset labels."
        arg_type = String
        default = DEFAULT_DATASET_CSV_PATH
        "--raw_data_dir"
        help = "Raw data directory containing `disprot/disprot_regions.json.gz`."
        arg_type = String
        default = DEFAULT_RAW_DATA_DIR
        "--output_dir"
        help = "Output directory for generated JSON artifact."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--output_filename"
        help = "Output gzip-compressed JSON filename written under `output_dir`."
        arg_type = String
        default = DEFAULT_OUTPUT_FILENAME
    end

    parsed = parse_args(args, settings)
    return (
        dataset_csv_path = normpath(parsed["dataset_csv_path"]),
        raw_data_dir = normpath(parsed["raw_data_dir"]),
        output_dir = normpath(parsed["output_dir"]),
        output_filename = String(parsed["output_filename"])
    )
end

function main()
    options = parse_command_line(ARGS)
    output_path = normpath(joinpath(options.output_dir, options.output_filename))
    output_readme_path = normpath(joinpath(options.output_dir, FINAL_DATASET_README_FILENAME))
    disprot_regions_path = normpath(joinpath(options.raw_data_dir, DISPROT_REGIONS_RELATIVE_PATH))

    final_dataset_table = DataFrame(CSV.File(options.dataset_csv_path))
    disprot_entries = read_disprot_entries_from_raw_data(options.raw_data_dir)
    json_records = build_dataset_json_records(final_dataset_table, disprot_entries)
    write_dataset_json_gz(json_records, output_path)
    write_dataset_json_readme_section(
        output_readme_path;
        output_json_path = _relpath_from_repo(output_path),
        output_json_filesystem_path = output_path,
        source_dataset_csv_path = _relpath_from_repo(options.dataset_csv_path),
        source_disprot_regions_path = _relpath_from_repo(disprot_regions_path),
        n_proteins = length(json_records),
        script_relpath = _relpath_from_repo(@__FILE__)
    )

    println("JSON dataset written to $(output_path) with $(length(json_records)) proteins.")
    nothing
end

main()
