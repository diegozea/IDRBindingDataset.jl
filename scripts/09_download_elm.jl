#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using Downloads
using IDRBindingDataset: download_to_file_with_retries, utc_now_iso8601

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_SCRATCH_DIR = normpath(joinpath(REPO_ROOT, "data", "scratch"))
const DEFAULT_ELM_SUBDIR = "elm"
const DEFAULT_INSTANCES_FILENAME = "instances.tsv"
const DEFAULT_ELM_INSTANCES_URL = "https://elm.eu.org/instances.tsv?q=*"
const DEFAULT_ELMS_INDEX_FILENAME = "elms_index.tsv"
const DEFAULT_ELMS_INDEX_URL = "http://elm.eu.org/elms/elms_index.tsv"
const DEFAULT_TIMEOUT_SECONDS = 180.0
const ELM_README_FILENAME = "README.md"
const ELM_HTTP_URL_PREFIX = "http://elm.eu.org/"
const ELM_HTTPS_URL_PREFIX = "https://elm.eu.org/"
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
  julia scripts/09_download_elm.jl --instances_url "https://elm.eu.org/instances.tsv?q=*" --elms_index_url "http://elm.eu.org/elms/elms_index.tsv" --timeout_seconds 180

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

"""
    _write_elm_readme(
        readme_path;
        instances_file_report,
        classes_index_file_report,
        instances_path,
        classes_index_path,
        downloaded_at_utc
    )

Write provenance metadata for a local ELM scratch download.
"""
function _write_elm_readme(
        readme_path::AbstractString;
        instances_file_report::NamedTuple,
        classes_index_file_report::NamedTuple,
        instances_path::AbstractString,
        classes_index_path::AbstractString,
        instances_rows::Int,
        classes_index_rows::Int,
        instances_header_rows::Int,
        classes_index_header_rows::Int,
        downloaded_at_utc::AbstractString
)::Nothing
    instances_relpath = relpath(instances_path, REPO_ROOT)
    classes_index_relpath = relpath(classes_index_path, REPO_ROOT)
    content = """
# ELM Scratch Download

This folder stores local ELM downloads for analysis workflows.
These files are not redistributed with this software due to ELM licensing restrictions.

- Download timestamp (UTC): `$(downloaded_at_utc)`

## Files
- `$(basename(instances_path))`
  - Requested source URL: `$(instances_file_report.requested_url)`
  - Downloaded from URL: `$(instances_file_report.resolved_url)`
  - Local file: `$(instances_relpath)`
  - Size (bytes): `$(filesize(instances_path))`
  - Data rows (CSV): `$(instances_rows)`
  - Header `#$(ELM_INSTANCES_COUNT_HEADER_KEY)`: `$(instances_header_rows)`
- `$(basename(classes_index_path))`
  - Requested source URL: `$(classes_index_file_report.requested_url)`
  - Downloaded from URL: `$(classes_index_file_report.resolved_url)`
  - Local file: `$(classes_index_relpath)`
  - Size (bytes): `$(filesize(classes_index_path))`
  - Data rows (CSV): `$(classes_index_rows)`
  - Header `#$(ELM_INDEX_COUNT_HEADER_KEY)`: `$(classes_index_header_rows)`
"""
    write(readme_path, content)
    return nothing
end

"""
    _build_elm_url_candidates(instances_url)

Build ordered candidate URLs for ELM download attempts.
"""
function _build_elm_url_candidates(instances_url::AbstractString)::Vector{String}
    primary = String(instances_url)
    candidates = String[primary]
    if startswith(primary, ELM_HTTPS_URL_PREFIX)
        fallback = replace(primary, ELM_HTTPS_URL_PREFIX => ELM_HTTP_URL_PREFIX)
        fallback ∉ candidates && push!(candidates, fallback)
    end
    return candidates
end

"""
    _validate_elm_tsv(path, required_columns, dataset_name)

Validate that an ELM TSV has non-empty content, required columns, and
header-reported row counts matching CSV row counts.
"""
function _parse_elm_header_metadata(path::AbstractString)::Dict{String, String}
    metadata = Dict{String, String}()
    open(path, "r") do io
        for line in eachline(io)
            startswith(line, "#") || break
            parsed_line = match(r"^#([^:]+):\s*(.*)\s*$", line)
            isnothing(parsed_line) && continue
            key = strip(String(parsed_line.captures[1]))
            value = strip(String(parsed_line.captures[2]))
            metadata[key] = value
        end
    end
    return metadata
end

function _parse_required_header_int(
        metadata::Dict{String, String},
        key::AbstractString,
        dataset_name::AbstractString
)::Int
    haskey(metadata, key) || throw(
        ArgumentError(
        "Downloaded $(dataset_name) TSV is missing required header metadata key: #$(key)."
    ),
    )
    return try
        parse(Int, metadata[key])
    catch
        throw(
            ArgumentError(
            "Downloaded $(dataset_name) TSV has non-integer header metadata value for #$(key): $(metadata[key])."
        ),
        )
    end
end

_count_elm_tsv_rows(path::AbstractString)::Int = length(
    CSV.File(path; delim = '\t', quotechar = '"', comment = "#")
)

function _validate_elm_tsv(
        path::AbstractString,
        required_columns,
        dataset_name::AbstractString,
        row_count_header_key::AbstractString
)::Nothing
    isfile(path) || throw(ArgumentError("Downloaded file does not exist: $(path)"))
    filesize(path) > 0 || throw(ArgumentError("Downloaded file is empty: $(path)"))

    parsed = CSV.File(path; delim = '\t', quotechar = '"', comment = "#")
    header_fields = String.(propertynames(parsed))
    header_set = Set(header_fields)
    missing_columns = [column for column in required_columns if column ∉ header_set]
    if !isempty(missing_columns)
        missing_text = join(missing_columns, ", ")
        header_text = join(header_fields, ", ")
        throw(
            ArgumentError(
            "Downloaded $(dataset_name) TSV is missing required columns: $(missing_text). Header found: $(header_text)."
        ),
        )
    end

    actual_rows = length(parsed)
    actual_rows == 0 &&
        throw(ArgumentError("Downloaded $(dataset_name) TSV has no data rows: $(path)"))

    metadata = _parse_elm_header_metadata(path)
    expected_rows = _parse_required_header_int(
        metadata,
        row_count_header_key,
        dataset_name
    )

    expected_rows == actual_rows || throw(
        ArgumentError(
        "Downloaded $(dataset_name) TSV row count mismatch: header #$(row_count_header_key)=$(expected_rows), CSV rows=$(actual_rows)."
    ),
    )
    return nothing
end

"""
    _download_elm_tsv(
        source_url,
        output_path;
        timeout_seconds,
        required_columns,
        dataset_name,
        row_count_header_key
    )

Download an ELM TSV using ordered URL candidates and validate the result.
Returns the URL that produced the validated file.
"""
function _download_elm_tsv(
        source_url::AbstractString,
        output_path::AbstractString;
        timeout_seconds::Float64,
        required_columns,
        dataset_name::AbstractString,
        row_count_header_key::AbstractString
)::String
    url_candidates = _build_elm_url_candidates(source_url)
    failures = String[]

    for url_candidate in url_candidates
        downloader = (url,
            destination) -> Downloads.download(
            url,
            destination;
            timeout = timeout_seconds
        )
        delays = startswith(url_candidate, ELM_HTTPS_URL_PREFIX) &&
                 length(url_candidates) > 1 ? ExponentialBackOff(n = 1) :
                 ExponentialBackOff(n = 5)

        try
            println("Trying $(dataset_name) URL: $(url_candidate)")
            download_to_file_with_retries(
                url_candidate,
                output_path;
                delays = delays,
                downloader = downloader
            )
            _validate_elm_tsv(
                output_path,
                required_columns,
                dataset_name,
                row_count_header_key
            )
            return url_candidate
        catch err
            push!(failures, "- $(url_candidate): $(sprint(showerror, err))")
            isfile(output_path) && rm(output_path; force = true)
        end
    end

    failure_text = join(failures, "\n")
    throw(
        ErrorException(
        "Failed to download/validate $(dataset_name) TSV.\n" *
        "Tried URLs:\n$(failure_text)"
    ),
    )
end

function main()
    options = parse_command_line(ARGS)
    options.timeout_seconds > 0 ||
        throw(ArgumentError("--timeout_seconds must be greater than zero."))
    elm_dir = normpath(joinpath(DEFAULT_SCRATCH_DIR, DEFAULT_ELM_SUBDIR))
    instances_path = normpath(joinpath(elm_dir, DEFAULT_INSTANCES_FILENAME))
    classes_index_path = normpath(joinpath(elm_dir, DEFAULT_ELMS_INDEX_FILENAME))

    mkpath(elm_dir)
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

    _write_elm_readme(
        joinpath(elm_dir, ELM_README_FILENAME);
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
