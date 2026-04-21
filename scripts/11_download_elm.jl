#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset:
                         _count_elm_tsv_rows,
                         _download_elm_tsv,
                         _parse_elm_header_metadata,
                         _parse_required_header_int,
                         _write_elm_download_readme,
                         remove_files,
                         utc_now_iso8601

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_SCRATCH_DIR = normpath(joinpath(REPO_ROOT, "data", "scratch"))
const DEFAULT_ELM_SUBDIR = "elm"
const DEFAULT_INSTANCES_FILENAME = "instances.tsv"
const DEFAULT_ELM_INSTANCES_URL = "https://elm.eu.org/instances.tsv?q=*"
const DEFAULT_ELMS_INDEX_FILENAME = "elms_index.tsv"
const DEFAULT_ELMS_INDEX_URL = "http://elm.eu.org/elms/elms_index.tsv"
const DEFAULT_TIMEOUT_SECONDS = 180.0
const ELM_README_FILENAME = "README.md"
const REQUIRED_ELM_INSTANCE_COLUMNS = (
    "Accession",
    "ELMIdentifier",
    "Primary_Acc",
    "Start",
    "End"
)
const REQUIRED_ELM_INDEX_COLUMNS = ("Accession", "ELMIdentifier", "Regex")
const ELM_INSTANCES_COUNT_HEADER_KEY = "NumInstances"
const ELM_INDEX_COUNT_HEADER_KEY = "Num_Classes"

"""
    parse_command_line(args)

Parse command-line options for ELM instances download into scratch storage.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download ELM instances TSV into data/scratch/elm.

The script tries the requested URL first. For `https://elm.eu.org/...` URLs,
it automatically falls back to `http://elm.eu.org/...` if HTTPS is unreachable.
It downloads both:
- ELM instances TSV (`instances.tsv`)
- ELM classes index TSV with motif regexes (`elms_index.tsv`)

Run from the repository root:
  julia scripts/11_download_elm.jl --instances_url "https://elm.eu.org/instances.tsv?q=*" --elms_index_url "http://elm.eu.org/elms/elms_index.tsv" --timeout_seconds 180

Output location is fixed to:
  data/scratch/elm/instances.tsv
  data/scratch/elm/elms_index.tsv

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--instances_url"
        help = "ELM endpoint used to download all motif instances."
        arg_type = String
        default = DEFAULT_ELM_INSTANCES_URL
        "--elms_index_url"
        help = "ELM classes index endpoint used to download motif regex definitions."
        arg_type = String
        default = DEFAULT_ELMS_INDEX_URL
        "--timeout_seconds"
        help = "HTTP timeout per download attempt in seconds."
        arg_type = Float64
        default = DEFAULT_TIMEOUT_SECONDS
    end

    parsed = parse_args(args, settings)
    return (
        instances_url = String(parsed["instances_url"]),
        elms_index_url = String(parsed["elms_index_url"]),
        timeout_seconds = Float64(parsed["timeout_seconds"])
    )
end

function main()
    options = parse_command_line(ARGS)
    options.timeout_seconds > 0 ||
        throw(ArgumentError("--timeout_seconds must be greater than zero."))
    elm_dir = normpath(joinpath(DEFAULT_SCRATCH_DIR, DEFAULT_ELM_SUBDIR))
    instances_path = normpath(joinpath(elm_dir, DEFAULT_INSTANCES_FILENAME))
    classes_index_path = normpath(joinpath(elm_dir, DEFAULT_ELMS_INDEX_FILENAME))
    readme_path = normpath(joinpath(elm_dir, ELM_README_FILENAME))

    mkpath(elm_dir)
    remove_files([instances_path, classes_index_path, readme_path])
    resolved_instances_url = _download_elm_tsv(
        options.instances_url,
        instances_path;
        timeout_seconds = options.timeout_seconds,
        required_columns = REQUIRED_ELM_INSTANCE_COLUMNS,
        dataset_name = "ELM instances",
        row_count_header_key = ELM_INSTANCES_COUNT_HEADER_KEY
    )
    resolved_elms_index_url = _download_elm_tsv(
        options.elms_index_url,
        classes_index_path;
        timeout_seconds = options.timeout_seconds,
        required_columns = REQUIRED_ELM_INDEX_COLUMNS,
        dataset_name = "ELM classes index",
        row_count_header_key = ELM_INDEX_COUNT_HEADER_KEY
    )

    instances_metadata = _parse_elm_header_metadata(instances_path)
    classes_metadata = _parse_elm_header_metadata(classes_index_path)
    instances_rows = _count_elm_tsv_rows(instances_path)
    classes_index_rows = _count_elm_tsv_rows(classes_index_path)
    instances_header_rows = _parse_required_header_int(
        instances_metadata,
        ELM_INSTANCES_COUNT_HEADER_KEY,
        "ELM instances"
    )
    classes_index_header_rows = _parse_required_header_int(
        classes_metadata,
        ELM_INDEX_COUNT_HEADER_KEY,
        "ELM classes index"
    )

    _write_elm_download_readme(
        readme_path;
        instances_file_report = (
            requested_url = options.instances_url,
            resolved_url = resolved_instances_url
        ),
        classes_index_file_report = (
            requested_url = options.elms_index_url,
            resolved_url = resolved_elms_index_url
        ),
        instances_path = instances_path,
        classes_index_path = classes_index_path,
        instances_rows = instances_rows,
        classes_index_rows = classes_index_rows,
        instances_header_rows = instances_header_rows,
        classes_index_header_rows = classes_index_header_rows,
        downloaded_at_utc = utc_now_iso8601()
    )
    println("ELM instances saved to $(instances_path)")
    println("ELM classes index saved to $(classes_index_path)")
    return nothing
end

main()
