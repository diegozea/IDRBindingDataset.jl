using CodecZlib
using Dates
using Downloads

"""
    build_disprot_regions_url(release)

Build the DisProt API URL for retrieving all entries of a release, including
ambiguous regions and excluding obsolete ones.
"""
function build_disprot_regions_url(release::AbstractString)::String
    return "https://disprot.org/api/search?release=$(release)&show_ambiguous=true&show_obsolete=false&format=json"
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
    download_to_file_with_retries(url, destination;
        delays = ExponentialBackOff(n = 5), downloader = Downloads.download,
        logger = println, log_success = true)

Download `url` into `destination` using atomic writes and bounded retries.
If all attempts fail, the last exception is rethrown.
"""
function _attempt_download(
    url::AbstractString,
    destination::AbstractString,
    downloader::Function,
    logger::Function,
    log_success::Bool,
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

function download_to_file_with_retries(
    url::AbstractString,
    destination::AbstractString;
    delays = ExponentialBackOff(n = 5),
    downloader::Function = Downloads.download,
    logger::Function = println,
    log_success::Bool = true,
)::Nothing
    mkpath(dirname(destination))

    retry(
        () -> _attempt_download(url, destination, downloader, logger, log_success);
        delays = delays,
        check = (_, err) -> begin
            logger(
                "Download failed for $(url): $(sprint(showerror, err)). Retrying with exponential backoff...",
            )
            return true
        end,
    )()

    return nothing
end

bytes_to_mb(bytes::Integer)::Float64 = bytes / 1_000_000
format_mb(bytes::Integer)::Float64 = round(bytes_to_mb(bytes); digits = 4)

"""
    download_to_gzip_file_with_retries(url, destination;
        delays = ExponentialBackOff(n = 5), downloader = Downloads.download,
        logger = println)

Download `url`, compress it as gzip, and save it to `destination`.
"""
function download_to_gzip_file_with_retries(
    url::AbstractString,
    destination::AbstractString;
    delays = ExponentialBackOff(n = 5),
    downloader::Function = Downloads.download,
    logger::Function = println,
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
            log_success = false,
        )

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
            "Saved $(filesize(destination)) bytes ($(format_mb(filesize(destination))) MB compressed) to $(destination)",
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
    remediation::AbstractString = "",
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
                join(["- $(path)" for path in missing_paths], "\n"),
            ],
            "\n",
        )
        error(message)
    end

    if !isempty(oversized)
        limit_in_mib = max_bytes / (1024^2)
        details = join(
            [
                "- $(path): $(size) bytes ($(round(size / (1024^2); digits = 2)) MiB)"
                for (path, size) in oversized
            ],
            "\n",
        )

        lines = [
            "One or more downloaded files exceed the configured Git/GitHub size limit.",
            "Limit ($(limit_name)): $(max_bytes) bytes ($(round(limit_in_mib; digits = 2)) MiB).",
            "Oversized files:",
            details,
        ]

        if !isempty(remediation)
            push!(lines, "Suggested remediation: $(remediation)")
        end

        error(join(lines, "\n"))
    end

    return nothing
end

"""
    write_download_report(report_path; status, timestamp_utc, release, idpo_version,
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
    idpo_version::AbstractString,
    sources,
    error_message::Union{Nothing, AbstractString} = nothing,
)::Nothing
    report_dir = dirname(report_path)
    total_size_bytes = 0

    open(report_path, "w") do io
        println(io, "# Raw Data Download")
        println(io)
        println(io, "- Status: `$(status)`")
        println(io, "- Download timestamp (UTC): `$(timestamp_utc)`")
        println(io, "- DisProt release: `$(release)`")
        println(io, "- IDPO version: `$(idpo_version)`")

        if error_message !== nothing
            println(io, "- Error: `$(error_message)`")
        end

        for (_, filename) in sources
            destination = joinpath(report_dir, filename)
            if isfile(destination)
                total_size_bytes += filesize(destination)
            end
        end

        println(
            io,
            "- Final size (compressed, total): `$(format_mb(total_size_bytes)) MB`",
        )

        println(io)
        println(io, "## Files")

        for (url, filename) in sources
            destination = joinpath(report_dir, filename)
            if isfile(destination)
                size_bytes = filesize(destination)
                println(
                    io,
                    "- `$(filename)` ($(size_bytes) bytes, $(format_mb(size_bytes)) MB)",
                )
            else
                println(io, "- `$(filename)` (missing)")
            end
            println(io, "  - Source: `$(url)`")
        end
    end

    return nothing
end
