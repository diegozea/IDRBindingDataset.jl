#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         FINAL_DATASET_LABEL_VALIDATION,
                         build_dataset_csv_table,
                         build_dataset_json_records,
                         build_final_cluster_dataset_table,
                         build_protein_elm_clusters_table,
                         read_disprot_entries_from_raw_data,
                         remove_files,
                         write_dataset_json_gz,
                         write_final_dataset_csv_gz

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_FIXED_LABELS_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "derived",
    "fixed",
    "cluster_fixed_labels.csv.gz"
))
const DEFAULT_ELM_MOTIF_CLUSTERING_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "derived",
    "elm",
    "elm_motif_clustering.csv.gz"
))
const DEFAULT_RAW_DATA_DIR = normpath(joinpath(REPO_ROOT, "data", "raw"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "dataset"))
const DEFAULT_OUTPUT_CSV_FILENAME = "cluster_dataset_labels.csv.gz"
const DEFAULT_OUTPUT_JSON_FILENAME = "cluster_dataset_labels.json.gz"
const DEFAULT_ELM_CLUSTER_SEPARATOR = "-"

"""
    parse_command_line(args)

Parse command-line options for the final dataset generation script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Create both final dataset artifacts in one run:
1) load fixed labels (`cluster_fixed_labels.csv.gz`),
2) aggregate ELM motif-cluster IDs per protein from
   `elm_motif_clustering.csv.gz`,
3) map `cluster_fixed_dataset_label` to final labels:
   - `training` -> `training`,
   - `caid` -> `testing`,
   - missing -> one configurable final label (default: `validation`),
4) add `sequence` and residue-level labels from DisProt,
5) write:
   - CSV (`cluster_dataset_labels.csv.gz`) with `sequence` and
     `residue_labels` as a compact string (`0`, `1`, `2` per residue),
   - JSON (`cluster_dataset_labels.json.gz`) with `residue_labels` as
     per-residue integer vectors.

Run from repository root:
  julia scripts/13_final_dataset.jl

Complete example:
  julia scripts/13_final_dataset.jl --fixed_labels_path data/derived/fixed/cluster_fixed_labels.csv.gz --elm_motif_clustering_path data/derived/elm/elm_motif_clustering.csv.gz --raw_data_dir data/raw --output_dir data/dataset --output_csv_filename cluster_dataset_labels.csv.gz --output_json_filename cluster_dataset_labels.json.gz --missing_label validation

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--fixed_labels_path"
        help = "Input compressed CSV with fixed dataset labels."
        arg_type = String
        default = DEFAULT_FIXED_LABELS_PATH
        "--elm_motif_clustering_path"
        help = "Input compressed CSV with protein-to-ELM motif-cluster rows."
        arg_type = String
        default = DEFAULT_ELM_MOTIF_CLUSTERING_PATH
        "--raw_data_dir"
        help = "Raw data directory containing DisProt source files."
        arg_type = String
        default = DEFAULT_RAW_DATA_DIR
        "--output_dir"
        help = "Output directory for final CSV and JSON dataset artifacts."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--output_csv_filename"
        help = "Output compressed CSV filename written under `output_dir`."
        arg_type = String
        default = DEFAULT_OUTPUT_CSV_FILENAME
        "--output_json_filename"
        help = "Output gzip-compressed JSON filename written under `output_dir`."
        arg_type = String
        default = DEFAULT_OUTPUT_JSON_FILENAME
        "--missing_label"
        help = "Label assigned to proteins whose `cluster_fixed_dataset_label` is missing (must be one of: training, validation, testing)."
        arg_type = String
        default = FINAL_DATASET_LABEL_VALIDATION
        "--elm_cluster_separator"
        help = "Separator used to join multiple ELM cluster IDs for one protein."
        arg_type = String
        default = DEFAULT_ELM_CLUSTER_SEPARATOR
    end

    parsed = parse_args(args, settings)
    (
        fixed_labels_path = normpath(parsed["fixed_labels_path"]),
        elm_motif_clustering_path = normpath(parsed["elm_motif_clustering_path"]),
        raw_data_dir = normpath(parsed["raw_data_dir"]),
        output_dir = normpath(parsed["output_dir"]),
        output_csv_filename = String(parsed["output_csv_filename"]),
        output_json_filename = String(parsed["output_json_filename"]),
        missing_label = String(parsed["missing_label"]),
        elm_cluster_separator = String(parsed["elm_cluster_separator"])
    )
end

function main()
    options = parse_command_line(ARGS)

    output_csv_path = normpath(joinpath(options.output_dir, options.output_csv_filename))
    output_json_path = normpath(joinpath(options.output_dir, options.output_json_filename))

    remove_files([output_csv_path, output_json_path])

    fixed_labels_table = DataFrame(CSV.File(options.fixed_labels_path))
    protein_elm_motif_clusters = DataFrame(CSV.File(options.elm_motif_clustering_path))
    protein_elm_clusters = build_protein_elm_clusters_table(
        protein_elm_motif_clusters;
        separator = options.elm_cluster_separator
    )
    final_dataset_table = build_final_cluster_dataset_table(
        fixed_labels_table,
        protein_elm_clusters;
        missing_label = options.missing_label
    )

    disprot_entries = read_disprot_entries_from_raw_data(options.raw_data_dir)
    json_records = build_dataset_json_records(final_dataset_table, disprot_entries)
    csv_table = build_dataset_csv_table(json_records)

    write_final_dataset_csv_gz(csv_table, output_csv_path)
    write_dataset_json_gz(json_records, output_json_path)

    println(
        "Final dataset artifacts written: CSV=$(output_csv_path), JSON=$(output_json_path), proteins=$(length(json_records)), final_clusters=$(length(unique(final_dataset_table.final_cluster_id)))."
    )

    nothing
end

main()
