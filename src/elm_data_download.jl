using CSV
using Downloads

const _ELM_HTTP_URL_PREFIX = "http://elm.eu.org/"
const _ELM_HTTPS_URL_PREFIX = "https://elm.eu.org/"

"""
    _write_elm_download_readme(
        readme_path;
        instances_file_report,
        classes_index_file_report,
        instances_path,
        classes_index_path,
        downloaded_at_utc
    )

Write provenance metadata for one local ELM scratch download.
"""
function _write_elm_download_readme(
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
    instances_relpath = _relpath_from_repo(instances_path)
    classes_index_relpath = _relpath_from_repo(classes_index_path)
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
  - Header `#NumInstances`: `$(instances_header_rows)`
- `$(basename(classes_index_path))`
  - Requested source URL: `$(classes_index_file_report.requested_url)`
  - Downloaded from URL: `$(classes_index_file_report.resolved_url)`
  - Local file: `$(classes_index_relpath)`
  - Size (bytes): `$(filesize(classes_index_path))`
  - Data rows (CSV): `$(classes_index_rows)`
  - Header `#Num_Classes`: `$(classes_index_header_rows)`
"""
    write(readme_path, content)
    nothing
end

"""
    _build_elm_url_candidates(source_url)

Build ordered candidate URLs for ELM download attempts.
"""
function _build_elm_url_candidates(source_url::AbstractString)::Vector{String}
    primary = String(source_url)
    candidates = String[primary]
    if startswith(primary, _ELM_HTTPS_URL_PREFIX)
        fallback = replace(primary, _ELM_HTTPS_URL_PREFIX => _ELM_HTTP_URL_PREFIX)
        fallback ∉ candidates && push!(candidates, fallback)
    end
    candidates
end

"""
    _parse_elm_header_metadata(path)

Parse one ELM TSV metadata header into key/value pairs.
"""
function _parse_elm_header_metadata(path::AbstractString)::Dict{String, String}
    metadata = Dict{String, String}()
    open(path, "r") do io
        for line in eachline(io)
            startswith(line, "#") || break
            parsed_line = match(r"^#([^:]+):\s*(.*)\s*$", line)
            isnothing(parsed_line) && continue
            key_capture = parsed_line.captures[1]
            value_capture = parsed_line.captures[2]
            (isnothing(key_capture) || isnothing(value_capture)) && continue

            key = strip(String(key_capture))
            value = strip(String(value_capture))
            metadata[key] = value
        end
    end
    metadata
end

"""
    _parse_required_header_int(metadata, key, dataset_name)

Read and parse one required integer metadata key from one parsed ELM header.
"""
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

    try
        parse(Int, metadata[key])
    catch
        throw(
            ArgumentError(
            "Downloaded $(dataset_name) TSV has non-integer header metadata value for #$(key): $(metadata[key])."
        ),
        )
    end
end

"""
    _count_elm_tsv_rows(path)

Count data rows in one ELM TSV while ignoring comment-prefixed metadata lines.
"""
_count_elm_tsv_rows(path::AbstractString)::Int = length(
    CSV.File(path; delim = '\t', quotechar = '"', comment = "#")
)

"""
    _validate_elm_tsv(path, required_columns, dataset_name, row_count_header_key)

Validate that one ELM TSV has non-empty content, required columns, and
header-reported row counts matching parsed row counts.
"""
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
    nothing
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

Download one ELM TSV using ordered URL candidates and validate the result.
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
        delays = startswith(url_candidate, _ELM_HTTPS_URL_PREFIX) &&
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
