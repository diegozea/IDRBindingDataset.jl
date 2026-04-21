#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC,
                         _filter_elm_instances_by_logic,
                         _parse_excluded_instance_logic,
                         _relpath_from_repo,
                         build_cluster_elm_motifs_table,
                         clean_directory,
                         write_elm_motifs_csv_gz,
                         write_elm_motifs_readme

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_INSTANCES_PATH = normpath(joinpath(
    REPO_ROOT, "data", "scratch", "elm", "instances.tsv"))
const DEFAULT_ELMS_INDEX_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "scratch",
    "elm",
    "elms_index.tsv"
))
const DEFAULT_CLUSTER_REPRESENTATIVES_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "derived",
    "clustering",
    "cluster_representatives.csv.gz"
))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "scratch", "elm", "mapping"))
const DEFAULT_OUTPUT_FILENAME = "motifs.csv.gz"
const ELM_MOTIFS_README_FILENAME = "README.md"
const DEFAULT_EXCLUDE_INSTANCE_LOGIC_TEXT = join(DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC, ",")

"""
    parse_command_line(args)

Parse command-line options for ELM motif mapping.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Map ELM motif instances to proteins in `cluster_representatives.csv.gz` and
write one intermediate motif table for downstream processing.

Protocol implemented by this script:
1) load `instances.tsv` and remove rows whose `InstanceLogic` is excluded,
2) map proteins using `Primary_Acc` and fallback tokens in `Accessions`,
3) keep only proteins found in cluster representatives,
4) map `ELMIdentifier` to `Regex` through `elms_index.tsv`,
5) write deterministic CSV and README statistics.

To avoid accidental redistribution of ELM-derived data, output is fixed to:
- `data/scratch/elm/mapping/motifs.csv.gz`
- `data/scratch/elm/mapping/README.md`

Run from repository root:
  julia scripts/13_mapping_elm_motifs.jl --instances_path data/scratch/elm/instances.tsv --elms_index_path data/scratch/elm/elms_index.tsv --cluster_representatives_path data/derived/clustering/cluster_representatives.csv.gz --exclude_instance_logic "false positive,true negative"

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--instances_path"
        help = "ELM instances TSV path (downloaded with scripts/11_download_elm.jl)."
        arg_type = String
        default = DEFAULT_INSTANCES_PATH
        "--elms_index_path"
        help = "ELM classes index TSV path with regex definitions."
        arg_type = String
        default = DEFAULT_ELMS_INDEX_PATH
        "--cluster_representatives_path"
        help = "Cluster representatives compressed CSV path."
        arg_type = String
        default = DEFAULT_CLUSTER_REPRESENTATIVES_PATH
        "--exclude_instance_logic"
        help = "Comma-separated InstanceLogic values to exclude (quotes optional, empty string keeps all rows)."
        arg_type = String
        default = DEFAULT_EXCLUDE_INSTANCE_LOGIC_TEXT
    end

    parsed = parse_args(args, settings)
    excluded_instance_logic = _parse_excluded_instance_logic(parsed["exclude_instance_logic"])

    (
        instances_path = normpath(parsed["instances_path"]),
        elms_index_path = normpath(parsed["elms_index_path"]),
        cluster_representatives_path = normpath(parsed["cluster_representatives_path"]),
        excluded_instance_logic = excluded_instance_logic
    )
end

function main()
    options = parse_command_line(ARGS)

    output_dir = DEFAULT_OUTPUT_DIR
    output_csv_path = normpath(joinpath(output_dir, DEFAULT_OUTPUT_FILENAME))
    output_readme_path = normpath(joinpath(output_dir, ELM_MOTIFS_README_FILENAME))

    # Rebuild this intermediate artifact from scratch on each run.
    clean_directory(output_dir)

    elm_instances = DataFrame(CSV.File(
        options.instances_path;
        delim = '\t',
        quotechar = '"',
        comment = "#"
    ))
    elms_index = DataFrame(CSV.File(
        options.elms_index_path;
        delim = '\t',
        quotechar = '"',
        comment = "#"
    ))
    cluster_representatives = DataFrame(CSV.File(options.cluster_representatives_path))

    filtered_instances = _filter_elm_instances_by_logic(
        elm_instances;
        excluded_instance_logic = options.excluded_instance_logic
    )

    motif_table = build_cluster_elm_motifs_table(
        cluster_representatives,
        filtered_instances,
        elms_index;
        excluded_instance_logic = String[]
    )

    write_elm_motifs_csv_gz(motif_table, output_csv_path)

    cluster_protein_count = length(unique(String.(cluster_representatives.protein_id)))
    write_elm_motifs_readme(
        output_readme_path;
        output_csv_path = _relpath_from_repo(output_csv_path),
        output_csv_filesystem_path = output_csv_path,
        source_cluster_representatives_path = _relpath_from_repo(
            options.cluster_representatives_path
        ),
        source_elm_instances_path = _relpath_from_repo(options.instances_path),
        source_elms_index_path = _relpath_from_repo(options.elms_index_path),
        script_relpath = _relpath_from_repo(@__FILE__),
        excluded_instance_logic = options.excluded_instance_logic,
        input_instance_rows = nrow(elm_instances),
        filtered_instance_rows = nrow(filtered_instances),
        motif_table = motif_table,
        cluster_protein_count = cluster_protein_count
    )

    proteins_with_motifs = length(unique(motif_table.protein_id))
    println(
        "ELM motif mapping completed: proteins_with_motifs=$(proteins_with_motifs), motif_rows=$(nrow(motif_table)), output=$(output_csv_path)."
    )

    nothing
end

main()
