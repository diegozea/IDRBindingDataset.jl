#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using CompariMotif
using DataFrames
using IDRBindingDataset:
                         DEFAULT_ELM_MOTIF_CLUSTER_MATCHIC_THRESHOLD,
                         DEFAULT_ELM_MOTIF_CLUSTER_MIN_SHARED_POSITIONS,
                         DEFAULT_ELM_MOTIF_CLUSTER_NORMIC_THRESHOLD,
                         _ensure_path_not_within_directory,
                         _parse_optional_score_threshold,
                         _relpath_from_repo,
                         build_protein_elm_motif_cluster_table,
                         build_regex_motif_clusters_table,
                         clean_directory,
                         write_elm_motif_clustering_csv_gz,
                         write_elm_motif_clustering_readme,
                         write_regex_motif_clusters_csv_gz

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_MOTIFS_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "scratch",
    "elm",
    "mapping",
    "motifs.csv.gz"
))
const DEFAULT_REGEX_CLUSTER_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT,
    "data",
    "scratch",
    "elm",
    "clustering"
))
const DEFAULT_REGEX_CLUSTER_OUTPUT_FILENAME = "regex.csv.gz"
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "derived", "elm"))
const DEFAULT_OUTPUT_FILENAME = "elm_motif_clustering.csv.gz"
const ELM_MOTIF_CLUSTERING_README_FILENAME = "README.md"
const DEFAULT_CLUSTERING_REFERENCE_DOI = "10.1093/bioinformatics/btv155"

"""
    parse_command_line(args)

Parse command-line options for ELM motif clustering.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Cluster ELM motif regexes using `CompariMotif.jl` and map resulting regex
clusters to proteins.

Inputs:
- `data/scratch/elm/mapping/motifs.csv.gz` (from `scripts/13_mapping_elm_motifs.jl`)

Outputs:
- scratch regex clusters: `data/scratch/elm/clustering/regex.csv.gz`
- derived protein mapping: `data/derived/elm/elm_motif_clustering.csv.gz`
- derived README: `data/derived/elm/README.md`

Default matching settings are aligned with criteria reported in:
Palopoli NT, Lythgow KT, Edwards RJ (2015). Bioinformatics 31(14):2284-2293.
DOI: `$(DEFAULT_CLUSTERING_REFERENCE_DOI)`.

Run from repository root:
  julia scripts/14_elm_motif_clustering.jl --motifs_path data/scratch/elm/mapping/motifs.csv.gz --regex_cluster_output_dir data/scratch/elm/clustering --regex_cluster_output_filename regex.csv.gz --output_dir data/derived/elm --output_filename elm_motif_clustering.csv.gz

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--motifs_path"
        help = "Input motifs compressed CSV path from scripts/13_mapping_elm_motifs.jl."
        arg_type = String
        default = DEFAULT_MOTIFS_PATH
        "--regex_cluster_output_dir"
        help = "Output directory for regex cluster scratch artifacts."
        arg_type = String
        default = DEFAULT_REGEX_CLUSTER_OUTPUT_DIR
        "--regex_cluster_output_filename"
        help = "Filename for regex cluster compressed CSV written under regex_cluster_output_dir."
        arg_type = String
        default = DEFAULT_REGEX_CLUSTER_OUTPUT_FILENAME
        "--output_dir"
        help = "Output directory for protein-level derived clustering artifacts."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--output_filename"
        help = "Filename for protein-level compressed CSV written under output_dir."
        arg_type = String
        default = DEFAULT_OUTPUT_FILENAME
        "--min_shared_positions"
        help = "Minimum number of matched positions required by CompariMotif for a regex comparison."
        arg_type = Int
        default = DEFAULT_ELM_MOTIF_CLUSTER_MIN_SHARED_POSITIONS
        "--normic_threshold"
        help = "Minimum normalized information content (NormIC) to keep a regex similarity edge."
        arg_type = Float64
        default = DEFAULT_ELM_MOTIF_CLUSTER_NORMIC_THRESHOLD
        "--matchic_threshold"
        help = "Minimum match information content (MatchIC) to keep a regex similarity edge."
        arg_type = Float64
        default = DEFAULT_ELM_MOTIF_CLUSTER_MATCHIC_THRESHOLD
        "--score_threshold"
        help = "Optional minimum CompariMotif score to keep a regex similarity edge (empty or 'nothing' disables score filtering)."
        arg_type = String
        default = ""
    end

    parsed = parse_args(args, settings)
    score_threshold = _parse_optional_score_threshold(parsed["score_threshold"])

    (
        motifs_path = normpath(parsed["motifs_path"]),
        regex_cluster_output_dir = normpath(parsed["regex_cluster_output_dir"]),
        regex_cluster_output_filename = String(parsed["regex_cluster_output_filename"]),
        output_dir = normpath(parsed["output_dir"]),
        output_filename = String(parsed["output_filename"]),
        min_shared_positions = Int(parsed["min_shared_positions"]),
        normic_threshold = Float64(parsed["normic_threshold"]),
        matchic_threshold = Float64(parsed["matchic_threshold"]),
        score_threshold = score_threshold
    )
end

function main()
    options = parse_command_line(ARGS)

    options.normic_threshold >= 0 ||
        throw(ArgumentError("--normic_threshold must be >= 0."))
    options.matchic_threshold >= 0 ||
        throw(ArgumentError("--matchic_threshold must be >= 0."))
    options.min_shared_positions >= 1 ||
        throw(ArgumentError("--min_shared_positions must be >= 1."))
    options.score_threshold === nothing || options.score_threshold >= 0 ||
        throw(ArgumentError("--score_threshold must be >= 0 when provided."))

    regex_cluster_output_path = normpath(joinpath(
        options.regex_cluster_output_dir,
        options.regex_cluster_output_filename
    ))
    output_csv_path = normpath(joinpath(options.output_dir, options.output_filename))
    output_readme_path = normpath(joinpath(options.output_dir, ELM_MOTIF_CLUSTERING_README_FILENAME))

    _ensure_path_not_within_directory(
        options.motifs_path,
        options.regex_cluster_output_dir;
        path_label = "--motifs_path",
        directory_label = "--regex_cluster_output_dir"
    )
    clean_directory(options.regex_cluster_output_dir)
    mkpath(options.output_dir)

    motif_table = DataFrame(CSV.File(options.motifs_path))
    regexes = sort(unique(String.(motif_table.regex)))

    compare_options = ComparisonOptions(;
        min_shared_positions = options.min_shared_positions,
        normalized_ic_cutoff = 0.0
    )
    regex_similarity_matrix = compare(regexes, compare_options)
    regex_similarity_table = DataFrame(to_column_table(regex_similarity_matrix))

    regex_cluster_table = build_regex_motif_clusters_table(
        regexes,
        regex_similarity_table;
        normic_threshold = options.normic_threshold,
        matchic_threshold = options.matchic_threshold,
        score_threshold = options.score_threshold
    )
    write_regex_motif_clusters_csv_gz(regex_cluster_table, regex_cluster_output_path)

    protein_cluster_table = build_protein_elm_motif_cluster_table(motif_table, regex_cluster_table)
    write_elm_motif_clustering_csv_gz(protein_cluster_table, output_csv_path)

    write_elm_motif_clustering_readme(
        output_readme_path;
        output_csv_path = _relpath_from_repo(output_csv_path),
        output_csv_filesystem_path = output_csv_path,
        source_motifs_path = _relpath_from_repo(options.motifs_path),
        source_regex_clusters_path = _relpath_from_repo(regex_cluster_output_path),
        min_shared_positions = options.min_shared_positions,
        normic_threshold = options.normic_threshold,
        matchic_threshold = options.matchic_threshold,
        score_threshold = options.score_threshold,
        regex_cluster_table = regex_cluster_table,
        protein_cluster_table = protein_cluster_table,
        script_relpath = _relpath_from_repo(@__FILE__)
    )

    println(
        "ELM motif clustering completed: regexes=$(nrow(regex_cluster_table)), proteins=$(length(unique(protein_cluster_table.protein_id))), rows=$(nrow(protein_cluster_table))."
    )

    nothing
end

main()
