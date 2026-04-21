#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
import IDRBindingDataset
using IDRBindingDataset:
                         CLUSTER_REPRESENTATIVES_FILENAME,
                         CLUSTERING_README_FILENAME,
                         DISPROT_FULL_PROTEIN_FASTA_GZ_FILENAME,
                         DISPROT_REGIONS_RELATIVE_PATH,
                         GLOBAL_CLUSTER_MIN_COVERAGE,
                         GLOBAL_CLUSTER_MIN_IDENTITY,
                         LOCAL_CLUSTER_MIN_IDENTITY,
                         LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH,
                         clean_directory,
                         build_cluster_representatives_table,
                         run_protein_sequence_clustering,
                         write_cluster_representatives_csv,
                         write_disprot_full_protein_fasta_gz,
                         write_sequence_clustering_readme

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "data", "derived", "clustering"))
const DEFAULT_GLOBAL_CLUSTER_MEMORY_LIMIT = "16G"

"""
    parse_command_line(args)

Parse command-line options for DIAMOND-based protein sequence clustering.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Build DisProt full-protein FASTA and run DIAMOND clustering in two modes:
1) global clustering with DIAMOND cluster + recluster (default: >=80% coverage and >=30% identity),
2) local clustering (default: >=80% identity on ungapped regions of length >=10 aa).

Final artifacts are:
- `disprot_full_proteins.fasta.gz` (full DisProt protein sequences clustered),
- `cluster_representatives.csv.gz` (protein ID + global/local representative),
- `README.md` (run summary).

Run from the repository root:
  julia scripts/09_clustering.jl --raw_data_dir data/raw --output_dir data/derived/clustering --threads 8 --global_cluster_memory_limit 16G

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--raw_data_dir"
        help = "Raw data directory containing disprot/disprot_regions.json.gz."
        arg_type = String
        default = RAW_DATA_DIR
        "--output_dir"
        help = "Output directory for final clustering artifacts."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--threads"
        help = "Number of CPU threads for DIAMOND commands."
        arg_type = Int
        default = Sys.CPU_THREADS
        "--global_cluster_memory_limit"
        help = "DIAMOND --memory-limit value for global diamond cluster mode (e.g., 16G)."
        arg_type = String
        default = DEFAULT_GLOBAL_CLUSTER_MEMORY_LIMIT
        "--global_min_coverage"
        help = "Global clustering minimum mutual coverage percentage."
        arg_type = Float64
        default = GLOBAL_CLUSTER_MIN_COVERAGE
        "--global_min_identity"
        help = "Global clustering minimum identity percentage."
        arg_type = Float64
        default = GLOBAL_CLUSTER_MIN_IDENTITY
        "--local_min_identity"
        help = "Local clustering minimum identity percentage."
        arg_type = Float64
        default = LOCAL_CLUSTER_MIN_IDENTITY
        "--local_min_ungapped_length"
        help = "Local clustering minimum ungapped aligned length (aa)."
        arg_type = Int
        default = LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH
    end

    parsed = parse_args(args, settings)

    return (
        raw_data_dir = normpath(parsed["raw_data_dir"]),
        output_dir = normpath(parsed["output_dir"]),
        threads = parsed["threads"],
        global_cluster_memory_limit = parsed["global_cluster_memory_limit"],
        global_min_coverage = parsed["global_min_coverage"],
        global_min_identity = parsed["global_min_identity"],
        local_min_identity = parsed["local_min_identity"],
        local_min_ungapped_length = parsed["local_min_ungapped_length"]
    )
end

function main()
    options = parse_command_line(ARGS)
    repo_root = normpath(joinpath(@__DIR__, ".."))

    clean_directory(options.output_dir)

    disprot_regions_path = joinpath(options.raw_data_dir, DISPROT_REGIONS_RELATIVE_PATH)
    disprot_fasta_output_path = joinpath(options.output_dir, DISPROT_FULL_PROTEIN_FASTA_GZ_FILENAME)
    representatives_output_path = joinpath(options.output_dir, CLUSTER_REPRESENTATIVES_FILENAME)
    readme_output_path = joinpath(options.output_dir, CLUSTERING_README_FILENAME)

    write_disprot_full_protein_fasta_gz(disprot_regions_path, disprot_fasta_output_path)
    results = run_protein_sequence_clustering(
        disprot_fasta_output_path;
        threads = options.threads,
        global_cluster_memory_limit = options.global_cluster_memory_limit,
        global_min_coverage = options.global_min_coverage,
        global_min_identity = options.global_min_identity,
        local_min_identity = options.local_min_identity,
        local_min_ungapped_length = options.local_min_ungapped_length
    )
    representatives = build_cluster_representatives_table(
        results.global_assignments,
        results.local_assignments
    )

    write_cluster_representatives_csv(representatives, representatives_output_path)

    global_cluster_count = length(unique(results.global_assignments.cluster_id))
    local_cluster_count = length(unique(results.local_assignments.cluster_id))
    final_cluster_count = length(unique(representatives.final_cluster_id))

    source_disprot_regions_relpath = relpath(disprot_regions_path, repo_root)
    disprot_fasta_output_relpath = relpath(disprot_fasta_output_path, repo_root)
    representatives_output_relpath = relpath(representatives_output_path, repo_root)

    write_sequence_clustering_readme(
        readme_output_path;
        source_disprot_regions_path = source_disprot_regions_relpath,
        source_disprot_regions_size_bytes = filesize(disprot_regions_path),
        disprot_fasta_path = disprot_fasta_output_relpath,
        disprot_fasta_size_bytes = filesize(disprot_fasta_output_path),
        representatives_path = representatives_output_relpath,
        representatives_size_bytes = filesize(representatives_output_path),
        sequence_count = results.sequence_count,
        global_cluster_count = global_cluster_count,
        local_cluster_count = local_cluster_count,
        final_cluster_count = final_cluster_count,
        global_method = results.global_method,
        local_method = results.local_method,
        global_min_identity = results.global_min_identity,
        global_min_coverage = results.global_min_coverage,
        local_min_identity = results.local_min_identity,
        local_min_ungapped_length = results.local_min_ungapped_length
    )

    println(
        "Clustering completed: global_clusters=$(global_cluster_count), local_clusters=$(local_cluster_count), final_clusters=$(final_cluster_count)"
    )

    return nothing
end

main()
