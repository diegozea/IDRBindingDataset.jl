#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using Printf
using IDRBindingDataset: remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SCRIPTS_DIR = normpath(@__DIR__)
const DEFAULT_SUMMARY_OUTPUT_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "dataset",
    "PIPELINE_RUN_SUMMARY.md"
))
const DEFAULT_PIPELINE_SCRIPTS = [
    "01_download_disprot_data.jl",
    "02_download_caid_data.jl",
    "03_download_sifts.jl",
    "05_download_invalid_regions.jl",
    "06_download_go_data.jl",
    "08_pdb_release_dates.jl",
    "09_clustering.jl",
    "10_fixed_labels.jl",
    "11_download_elm.jl",
    "13_mapping_elm_motifs.jl",
    "14_elm_motif_clustering.jl",
    "15_final_dataset.jl",
    "16_check_idr_binding.jl"
]

"""
    _is_numbered_pipeline_script(filename)

Return `true` when `filename` matches the numbered pipeline script pattern
`NN_name.jl`.
"""
function _is_numbered_pipeline_script(filename::AbstractString)::Bool
    occursin(r"^\d{2}_.+\.jl$", String(filename))
end

"""
    _discover_numbered_pipeline_scripts(scripts_dir)

Discover and sort all numbered pipeline scripts under `scripts_dir`.
"""
function _discover_numbered_pipeline_scripts(
        scripts_dir::AbstractString
)::Vector{String}
    scripts = filter(_is_numbered_pipeline_script, readdir(scripts_dir))
    sort!(scripts)
    isempty(scripts) && throw(
        ArgumentError("No numbered pipeline scripts were found in $(scripts_dir)."),
    )
    scripts
end

"""
    _parse_script_list(value)

Parse one comma-separated script list and return normalized script filenames.
"""
function _parse_script_list(value::AbstractString)::Vector{String}
    scripts = String[]
    for token in split(value, ',')
        script = strip(token)
        isempty(script) && continue
        push!(scripts, script)
    end
    scripts
end

"""
    parse_command_line(args)

Parse command-line options for pipeline orchestration.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Run the full numbered Julia pipeline in `scripts/` from the repository root.

By default, this executes the dependency-aware ordered pipeline defined in
`DEFAULT_PIPELINE_SCRIPTS`.

Run from repository root:
  julia scripts/run_pipeline.jl

Complete example:
  julia scripts/run_pipeline.jl --julia_cmd julia --summary_output_path data/dataset/PIPELINE_RUN_SUMMARY.md --continue_on_error

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--scripts"
        help = "Optional comma-separated script filenames under `scripts/` to run in order. Default: dependency-aware pipeline order."
        arg_type = String
        default = ""
        "--julia_cmd"
        help = "Julia executable used to run each pipeline script."
        arg_type = String
        default = "julia"
        "--summary_output_path"
        help = "Markdown summary output path."
        arg_type = String
        default = DEFAULT_SUMMARY_OUTPUT_PATH
        "--continue_on_error"
        help = "Continue running subsequent scripts after a failure."
        action = :store_true
    end

    parsed = parse_args(args, settings)
    scripts = if isempty(strip(String(parsed["scripts"])))
        discovered_scripts = _discover_numbered_pipeline_scripts(SCRIPTS_DIR)
        missing_from_default = setdiff(discovered_scripts, DEFAULT_PIPELINE_SCRIPTS)
        extras_in_default = setdiff(DEFAULT_PIPELINE_SCRIPTS, discovered_scripts)
        isempty(missing_from_default) || throw(
            ArgumentError(
            "DEFAULT_PIPELINE_SCRIPTS is missing numbered scripts: $(join(sort!(missing_from_default), ", "))."
        ),
        )
        isempty(extras_in_default) || throw(
            ArgumentError(
            "DEFAULT_PIPELINE_SCRIPTS references missing scripts: $(join(sort!(extras_in_default), ", "))."
        ),
        )
        copy(DEFAULT_PIPELINE_SCRIPTS)
    else
        _parse_script_list(String(parsed["scripts"]))
    end
    isempty(scripts) && throw(
        ArgumentError("No scripts were selected for execution."),
    )

    (
        scripts = scripts,
        julia_cmd = String(parsed["julia_cmd"]),
        summary_output_path = normpath(String(parsed["summary_output_path"])),
        continue_on_error = Bool(parsed["continue_on_error"])
    )
end

"""
    _run_one_script(script_filename; julia_cmd)

Run one script in a subprocess and return execution metadata.
"""
function _run_one_script(
        script_filename::AbstractString;
        julia_cmd::AbstractString
)
    script_path = joinpath(SCRIPTS_DIR, script_filename)
    isfile(script_path) || throw(
        ArgumentError("Script does not exist: $(script_path)"),
    )

    script_relpath = joinpath("scripts", script_filename)
    println(">>> Running $(script_relpath)")
    started_ns = time_ns()
    command = Cmd(Cmd([julia_cmd, script_relpath]); dir = REPO_ROOT)
    process = run(command; wait = false)
    wait(process)
    elapsed_seconds = (time_ns() - started_ns) / 1e9
    status = process.exitcode == 0 ? "PASS" : "FAIL"

    println("<<< $(script_relpath): $(status) (exit=$(process.exitcode), elapsed=$(round(elapsed_seconds; digits = 2))s)")
    (
        script = String(script_filename),
        status = status,
        exit_code = Int(process.exitcode),
        elapsed_seconds = elapsed_seconds
    )
end

"""
    _write_summary(path, results)

Write one Markdown execution summary.
"""
function _write_summary(path::AbstractString, results)::Nothing
    mkpath(dirname(path))
    remove_files([path])

    n_pass = count(result -> result.exit_code == 0, results)
    n_fail = length(results) - n_pass
    total_elapsed_seconds = sum(result.elapsed_seconds for result in results)

    lines = String[
    "# Pipeline Run Summary",
    "",
    "- Total scripts: `$(length(results))`",
    "- Passed: `$(n_pass)`",
    "- Failed: `$(n_fail)`",
    "- Total elapsed (s): `$(round(total_elapsed_seconds; digits = 2))`",
    "",
    "| Script | Status | Exit Code | Elapsed (s) |",
    "|---|---|---:|---:|"
]

    for result in results
        push!(
            lines,
            @sprintf("| %s | %s | %d | %.2f |",
                result.script,
                result.status,
                result.exit_code,
                result.elapsed_seconds)
        )
    end

    write(path, join(lines, "\n") * "\n")
    nothing
end

function main()
    options = parse_command_line(ARGS)
    results = NamedTuple{
        (:script, :status, :exit_code, :elapsed_seconds),
        Tuple{String, String, Int, Float64}
    }[]

    for script in options.scripts
        result = _run_one_script(script; julia_cmd = options.julia_cmd)
        push!(results, result)

        if result.exit_code != 0 && !options.continue_on_error
            _write_summary(options.summary_output_path, results)
            throw(ErrorException(
                "Pipeline failed at scripts/$(script). Summary: $(options.summary_output_path)"
            ))
        end
    end

    _write_summary(options.summary_output_path, results)

    n_fail = count(result -> result.exit_code != 0, results)
    if n_fail > 0
        throw(ErrorException(
            "Pipeline finished with $(n_fail) failing script(s). Summary: $(options.summary_output_path)"
        ))
    end

    println("Pipeline completed successfully. Summary: $(options.summary_output_path)")
    nothing
end

main()
