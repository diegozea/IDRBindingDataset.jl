using CodecZlib
using Dates
using Downloads

const _REPO_ROOT_PATH = normpath(joinpath(@__DIR__, ".."))
const _IDPO_LATEST_RELEASE_API_URL = "https://api.github.com/repos/BioComputingUP/idpo/releases/latest"
const _IDPO_RELEASE_ASSET_FILENAME = "idpo.json"

"""
    _relpath_from_repo(path)

Return one normalized path relative to this repository root.
"""
function _relpath_from_repo(path::AbstractString)::String
    relpath(normpath(path), _REPO_ROOT_PATH)
end

"""
    build_disprot_regions_url(release)

Build the DisProt API URL for retrieving all entries of a release, including
ambiguous regions and excluding obsolete ones.
"""
function build_disprot_regions_url(release::AbstractString)::String
    return "https://disprot.org/api/search?release=$(release)&show_ambiguous=true&show_obsolete=false&format=json"
end

"""
    parse_disprot_release(release)

Parse a DisProt release tag in `YYYY_MM` format into a `Date` at the first day
of that month.
"""
function parse_disprot_release(release::AbstractString)::Date
    try
        return Date(release, dateformat"y_m")
    catch err
        if err isa ArgumentError
            throw(ArgumentError("Invalid DisProt release tag '$(release)'. Expected format YYYY_MM with month 01-12."))
        end
        rethrow(err)
    end
end

"""
    _escape_regex_literal(value)

Escape regular-expression metacharacters in `value`.
"""
function _escape_regex_literal(value::AbstractString)::String
    replace(value, r"([\\.^$|?*+()\[\]{}])" => s"\\\1")
end

"""
    _download_url_to_string(url; downloader = Downloads.download)

Download `url` and return the payload as a UTF-8 `String`.
"""
function _download_url_to_string(
        url::AbstractString;
        downloader::Function = Downloads.download
)::String
    tmp_path, tmp_io = mktemp(; cleanup = false)
    close(tmp_io)
    try
        downloader(url, tmp_path)
        return read(tmp_path, String)
    finally
        if isfile(tmp_path)
            rm(tmp_path; force = true)
        end
    end
end

"""
    _extract_json_key_value(json_text, key)

Extract a JSON string value for top-level key `key` using a strict regex match.
Throws an `ErrorException` when the key cannot be found.
"""
function _extract_json_key_value(
        json_text::AbstractString,
        key::AbstractString
)::String
    escaped_key = _escape_regex_literal(key)
    pattern = Regex("\"$(escaped_key)\"\\s*:\\s*\"([^\"]+)\"")
    matched = match(pattern, json_text)
    if matched === nothing
        error("Could not find key \"$(key)\" in JSON payload.")
    end
    return replace(only(matched.captures), "\\/" => "/")
end

"""
    _extract_release_asset_download_url(json_text, asset_filename)

Extract the `browser_download_url` for release asset `asset_filename` from
GitHub release API JSON payload.
"""
function _extract_release_asset_download_url(
        json_text::AbstractString,
        asset_filename::AbstractString
)::String
    escaped_name = _escape_regex_literal(asset_filename)
    pattern = Regex(
        "\"name\"\\s*:\\s*\"$(escaped_name)\"[\\s\\S]*?\"browser_download_url\"\\s*:\\s*\"([^\"]+)\"",
        "s"
    )
    matched = match(pattern, json_text)
    if matched === nothing
        error(
            "GitHub latest-release payload does not include asset \"$(asset_filename)\"."
        )
    end
    return replace(only(matched.captures), "\\/" => "/")
end

"""
    _resolve_latest_idpo_release(; api_url = _IDPO_LATEST_RELEASE_API_URL,
        asset_filename = _IDPO_RELEASE_ASSET_FILENAME,
        downloader = Downloads.download)

Resolve the latest IDPO release tag and download URL for `asset_filename`
from the GitHub Releases API.
"""
function _resolve_latest_idpo_release(;
        api_url::AbstractString = _IDPO_LATEST_RELEASE_API_URL,
        asset_filename::AbstractString = _IDPO_RELEASE_ASSET_FILENAME,
        downloader::Function = Downloads.download
)::NamedTuple{(:release_tag, :asset_url), Tuple{String, String}}
    release_payload = _download_url_to_string(api_url; downloader = downloader)
    release_tag = _extract_json_key_value(release_payload, "tag_name")
    asset_url = _extract_release_asset_download_url(release_payload, asset_filename)
    return (; release_tag, asset_url)
end

"""
    _validate_idpo_obojson_payload(path)

Validate that `path` contains a plausible IDPO OBO-JSON payload (not HTML
fallback content), including the expected `graphs` section and IDPO term IDs.
"""
function _validate_idpo_obojson_payload(path::AbstractString)::Nothing
    payload = read(path, String)
    stripped_payload = strip(payload)
    if isempty(stripped_payload)
        error("IDPO payload is empty.")
    end

    payload_prefix_length = min(length(stripped_payload), 512)
    payload_prefix = lowercase(first(stripped_payload, payload_prefix_length))
    if occursin("<!doctype html", payload_prefix) || occursin("<html", payload_prefix)
        error("IDPO payload appears to be HTML, not JSON.")
    end

    if !startswith(stripped_payload, "{")
        error("IDPO payload does not start with a JSON object.")
    end

    if !occursin(r"\"graphs\"\s*:", payload)
        error("IDPO payload is missing the expected \"graphs\" section.")
    end

    if !occursin(r"IDPO_[0-9]{7}", payload)
        error("IDPO payload is missing expected IDPO term identifiers.")
    end

    return nothing
end

"""
    utc_now_iso8601(; clock = () -> now(UTC))

Return the current UTC timestamp as `yyyy-mm-ddTHH:MM:SSZ`.

The `clock` keyword is injectable to make tests deterministic.
"""
function utc_now_iso8601(; clock::Function = () -> now(UTC))::String
    return string(Dates.format(clock(), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z")
end

"""
    clean_directory(path)

Remove `path` recursively if it exists, then recreate it as an empty directory.
"""
function clean_directory(path::AbstractString)::Nothing
    if isdir(path)
        rm(path; recursive = true, force = true)
    end
    mkpath(path)
    return nothing
end

"""
    _is_path_within_directory(path, directory)

Return `true` when `path` is the same as, or nested under, `directory` after
normalizing both to absolute paths.
"""
function _is_path_within_directory(
        path::AbstractString,
        directory::AbstractString
)::Bool
    normalized_path = normpath(abspath(path))
    normalized_directory = normpath(abspath(directory))

    normalized_path == normalized_directory && return true
    startswith(normalized_path, joinpath(normalized_directory, ""))
end

"""
    _ensure_path_not_within_directory(path, directory;
        path_label = "path", directory_label = "directory")

Throw an `ArgumentError` when `path` is the same as, or nested under,
`directory`.
"""
function _ensure_path_not_within_directory(
        path::AbstractString,
        directory::AbstractString;
        path_label::AbstractString = "path",
        directory_label::AbstractString = "directory"
)::Nothing
    _is_path_within_directory(path, directory) || return nothing

    normalized_path = normpath(path)
    normalized_directory = normpath(directory)
    throw(
        ArgumentError(
        "Refusing to clean $(directory_label)=$(repr(normalized_directory)) because $(path_label)=$(repr(normalized_path)) is inside it. Choose a non-overlapping directory."
    ),
    )
end

"""
    remove_files(paths)

Remove each file path in `paths` if it exists. Parent directories are kept.
"""
function remove_files(paths)::Nothing
    for path in paths
        path_string = String(path)
        if isfile(path_string)
            rm(path_string; force = true)
        end
    end
    return nothing
end

"""
    _attempt_download(url, destination, downloader, logger, log_success)

Run a single atomic download attempt into `destination`.

The payload is first downloaded into a temporary file in the destination folder
and only moved into place on success.
"""
function _attempt_download(
        url::AbstractString,
        destination::AbstractString,
        downloader::Function,
        logger::Function,
        log_success::Bool
)::Nothing
    tmp_path, tmp_io = mktemp(dirname(destination); cleanup = false)
    close(tmp_io)

    try
        # Write to a temporary file first to avoid partial outputs on failure.
        downloader(url, tmp_path)
        mv(tmp_path, destination; force = true)
        if log_success
            logger("Saved $(filesize(destination)) bytes to $(destination)")
        end
        return nothing
    catch err
        if isfile(tmp_path)
            rm(tmp_path; force = true)
        end
        rethrow(err)
    end
end

"""
    download_to_file_with_retries(url, destination;
        delays = ExponentialBackOff(n = 5), downloader = Downloads.download,
        logger = println, log_success = true)

Download `url` into `destination` using atomic writes and bounded retries.
If all attempts fail, the last exception is rethrown.
"""
function download_to_file_with_retries(
        url::AbstractString,
        destination::AbstractString;
        delays = ExponentialBackOff(n = 5),
        downloader::Function = Downloads.download,
        logger::Function = println,
        log_success::Bool = true
)::Nothing
    mkpath(dirname(destination))

    retry(
        () -> _attempt_download(url, destination, downloader, logger, log_success);
        delays = delays,
        check = (_,
            err) -> begin
            logger(
                "Download failed for $(url): $(sprint(showerror, err)). Retrying with exponential backoff...",
            )
            return true
        end
    )()

    return nothing
end

"""
    _bytes_to_mb(bytes)

Convert bytes to decimal megabytes.
"""
_bytes_to_mb(bytes::Integer)::Float64 = bytes / 1_000_000

"""
    _format_mb(bytes)

Convert bytes to decimal megabytes and round to 4 decimal digits.
"""
_format_mb(bytes::Integer)::Float64 = round(_bytes_to_mb(bytes); digits = 4)

"""
    download_to_gzip_file_with_retries(url, destination;
        delays = ExponentialBackOff(n = 5), downloader = Downloads.download,
        logger = println, precompress_validator = nothing)

Download `url`, compress it as gzip, and save it to `destination`.

If `precompress_validator` is provided, it is called with the path to the
downloaded temporary uncompressed file before compression starts.
"""
function download_to_gzip_file_with_retries(
        url::AbstractString,
        destination::AbstractString;
        delays = ExponentialBackOff(n = 5),
        downloader::Function = Downloads.download,
        logger::Function = println,
        precompress_validator::Union{Nothing, Function} = nothing
)::Nothing
    mkpath(dirname(destination))
    raw_tmp_path, raw_tmp_io = mktemp(dirname(destination); cleanup = false)
    close(raw_tmp_io)
    gzip_tmp_path, gzip_tmp_io = mktemp(dirname(destination); cleanup = false)
    close(gzip_tmp_io)

    try
        download_to_file_with_retries(
            url,
            raw_tmp_path;
            delays = delays,
            downloader = downloader,
            logger = logger,
            log_success = false
        )

        if precompress_validator !== nothing
            precompress_validator(raw_tmp_path)
        end

        open(raw_tmp_path, "r") do input_io
            open(gzip_tmp_path, "w") do output_io
                gzip_io = GzipCompressorStream(output_io)
                try
                    while !eof(input_io)
                        write(gzip_io, read(input_io, 1 << 20))
                    end
                finally
                    close(gzip_io)
                end
            end
        end

        mv(gzip_tmp_path, destination; force = true)
        logger(
            "Saved $(filesize(destination)) bytes ($(_format_mb(filesize(destination))) MB compressed) to $(destination)",
        )
        return nothing
    catch err
        rethrow(err)
    finally
        if isfile(raw_tmp_path)
            rm(raw_tmp_path; force = true)
        end
        if isfile(gzip_tmp_path)
            rm(gzip_tmp_path; force = true)
        end
    end
end

"""
    ensure_files_within_size_limit(paths; max_bytes, limit_name, remediation)

Validate that all files in `paths` exist and are smaller than or equal to
`max_bytes`.

If any file is missing or exceeds the configured limit, this function throws
an `ErrorException` with a clear multi-line diagnostic.
"""
function ensure_files_within_size_limit(
        paths;
        max_bytes::Integer,
        limit_name::AbstractString,
        remediation::AbstractString = ""
)::Nothing
    max_bytes > 0 || throw(ArgumentError("`max_bytes` must be greater than zero."))

    missing_paths = String[]
    oversized = Tuple{String, Int64}[]

    for path in paths
        path_string = String(path)
        if !isfile(path_string)
            push!(missing_paths, path_string)
            continue
        end

        size_in_bytes = filesize(path_string)
        if size_in_bytes > max_bytes
            push!(oversized, (path_string, size_in_bytes))
        end
    end

    if !isempty(missing_paths)
        message = join(
            [
                "Cannot validate file size limits because these files are missing:",
                join(["- $(path)" for path in missing_paths], "\n")
            ],
            "\n"
        )
        error(message)
    end

    if !isempty(oversized)
        limit_in_mib = max_bytes / (1024^2)
        details = join(
            ["- $(path): $(size) bytes ($(round(size / (1024^2); digits = 2)) MiB)"
             for (path, size) in oversized],
            "\n"
        )

        lines = [
            "One or more downloaded files exceed the configured Git/GitHub size limit.",
            "Limit ($(limit_name)): $(max_bytes) bytes ($(round(limit_in_mib; digits = 2)) MiB).",
            "Oversized files:",
            details
        ]

        if !isempty(remediation)
            push!(lines, "Suggested remediation: $(remediation)")
        end

        error(join(lines, "\n"))
    end

    return nothing
end

"""
    write_download_report(report_path; status, timestamp_utc, release, idpo_release_tag,
        sources, error_message = nothing)

Write a Markdown report describing the last raw-data download run.

`sources` must be an iterable of `(url, filename)` pairs where each `filename`
is expected to live in `dirname(report_path)`.
"""
function write_download_report(
        report_path::AbstractString;
        status::AbstractString,
        timestamp_utc::AbstractString,
        release::AbstractString,
        idpo_release_tag::AbstractString,
        sources,
        error_message::Union{Nothing, AbstractString} = nothing
)::Nothing
    report_dir = dirname(report_path)
    total_size_bytes = 0

    for (_, filename) in sources
        destination = joinpath(report_dir, filename)
        if isfile(destination)
            total_size_bytes += filesize(destination)
        end
    end

    content = """
# Raw Data Download

- Status: `$(status)`
- Download timestamp (UTC): `$(timestamp_utc)`
- DisProt release: `$(release)`
- IDPO release tag: `$(idpo_release_tag)`
"""

    if error_message !== nothing
        content *= "- Error: `$(error_message)`\n"
    end

    content *= """
- Final size (compressed, total): `$(_format_mb(total_size_bytes)) MB`

## Files
"""

    for (url, filename) in sources
        destination = joinpath(report_dir, filename)
        if isfile(destination)
            size_bytes = filesize(destination)
            content *= """
- `$(filename)` ($(size_bytes) bytes, $(_format_mb(size_bytes)) MB)
  - Source: `$(url)`
"""
        else
            content *= """
- `$(filename)` (missing)
  - Source: `$(url)`
"""
        end
    end

    write(report_path, content)
    return nothing
end
