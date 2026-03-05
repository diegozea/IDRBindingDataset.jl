#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         FINAL_DATASET_LABEL_VALIDATION,
                         _relpath_from_repo,
                         build_final_cluster_dataset_table,
                         build_protein_elm_clusters_table,
                         clean_directory,
                         write_final_dataset_csv_gz,
                         write_final_dataset_readme

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
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "dataset"))
const DEFAULT_OUTPUT_FILENAME = "cluster_dataset_labels.csv.gz"
const DEFAULT_ELM_CLUSTER_SEPARATOR = "-"
const FINAL_DATASET_README_FILENAME = "README.md"

"""
    parse_command_line(args)

Parse command-line options for final dataset generation.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Create the final per-protein dataset by:
1) loading fixed labels (`cluster_fixed_labels.csv.gz`),
2) aggregating ELM motif-cluster IDs per protein from
   `elm_motif_clustering.csv.gz`,
3) mapping `cluster_fixed_dataset_label` to final labels:
   - `training` -> `training`,
   - `caid` -> `testing`,
   - missing -> one configurable final label (default: `validation`),
4) writing one final compressed CSV and README summary under `data/dataset/`.

Run from repository root:
  julia scripts/13_assing_missing_labels.jl --fixed_labels_path data/derived/fixed/cluster_fixed_labels.csv.gz --elm_motif_clustering_path data/derived/elm/elm_motif_clustering.csv.gz --output_dir data/dataset --output_filename cluster_dataset_labels.csv.gz --missing_label validation

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
        "--output_dir"
        help = "Output directory for the final dataset CSV and README."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--output_filename"
        help = "Output compressed CSV filename written under `output_dir`."
        arg_type = String
        default = DEFAULT_OUTPUT_FILENAME
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
        output_dir = normpath(parsed["output_dir"]),
        output_filename = String(parsed["output_filename"]),
        missing_label = String(parsed["missing_label"]),
        elm_cluster_separator = String(parsed["elm_cluster_separator"])
    )
end

function main()
    options = parse_command_line(ARGS)
    output_csv_path = normpath(joinpath(options.output_dir, options.output_filename))
    output_readme_path = normpath(joinpath(options.output_dir, FINAL_DATASET_README_FILENAME))

    clean_directory(options.output_dir)

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

    write_final_dataset_csv_gz(final_dataset_table, output_csv_path)
    write_final_dataset_readme(
        output_readme_path;
        output_csv_path = _relpath_from_repo(output_csv_path),
        output_csv_filesystem_path = output_csv_path,
        source_fixed_labels_path = _relpath_from_repo(options.fixed_labels_path),
        source_elm_motif_clustering_path = _relpath_from_repo(
            options.elm_motif_clustering_path
        ),
        final_dataset_table = final_dataset_table,
        missing_label = options.missing_label,
        elm_cluster_separator = options.elm_cluster_separator,
        script_relpath = _relpath_from_repo(@__FILE__)
    )

    println(
        "Final dataset written to $(output_csv_path) with $(nrow(final_dataset_table)) proteins and $(length(unique(final_dataset_table.final_cluster_id))) final clusters."
    )

    nothing
end

main()
