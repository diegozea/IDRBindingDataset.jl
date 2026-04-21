const CAID_FOLDER_README_FILENAME = "README.md"
const CAID_3_INVALID_REGIONS_FILENAME = "invalid_regions_fullCAID3.tsv"
const CAID_3_INVALID_REGIONS_URL = "https://raw.githubusercontent.com/BioComputingUP/caid-reference/main/data/disprot/invalid_regions_fullCAID3.tsv"

const CAID_2_SOURCES = Dict{String, String}(
    "disorder_nox.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/2/disorder_nox.fasta",
    "disorder_pdb.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/2/disorder_pdb.fasta",
    "binding.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/2/binding.fasta",
    "linker.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/2/linker.fasta"
)

const CAID_3_SOURCES = Dict{String, String}(
    "disorder_nox.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/3/disorder_nox.fasta",
    "disorder_pdb.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/3/disorder_pdb.fasta",
    "binding.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/3/binding.fasta",
    "binding_idr.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/3/binding_idr.fasta",
    "linker.fasta" => "https://caid.idpcentral.org/assets/sections/challenge/static/references/3/linker.fasta"
)

const CAID_DATASETS = (
    (folder = "caid_2", challenge = "CAID-2", sources = CAID_2_SOURCES),
    (folder = "caid_3", challenge = "CAID-3", sources = CAID_3_SOURCES)
)

"""
    _caid_fasta_prefix_validator(path)

Validate that a downloaded CAID FASTA file starts with `>DP`.
"""
function _caid_fasta_prefix_validator(path::AbstractString)::Nothing
    open(path, "r") do io
        eof(io) && error("Downloaded CAID FASTA file is empty: $(path)")
        first_line = readline(io)
        startswith(first_line, ">DP") ||
            error("Downloaded CAID FASTA file does not start with >DP: $(path)")
    end
    return nothing
end

"""
    _validate_caid3_invalid_regions_tsv(path)

Validate that one downloaded CAID-3 invalid-regions TSV exists, is non-empty,
and starts with the expected header prefix.
"""
function _validate_caid3_invalid_regions_tsv(path::AbstractString)::Nothing
    isfile(path) || error("Downloaded CAID-3 invalid-regions TSV does not exist: $(path)")

    open(path, "r") do io
        eof(io) && error("Downloaded CAID-3 invalid-regions TSV is empty: $(path)")
        header = chomp(readline(io))
        startswith(header, "disprot_id\tregion_id\tvalid") || error(
            "Downloaded CAID-3 invalid-regions TSV has unexpected header: $(repr(header))",
        )
    end

    return nothing
end

"""
    _caid_entry_size_text(filename, size_in_bytes)

Render one human-readable size string for a CAID raw-data file.
"""
function _caid_entry_size_text(
        filename::AbstractString,
        size_in_bytes::Integer
)::String
    compression_suffix = endswith(lowercase(String(filename)), ".gz") ? " compressed" : ""
    "$(size_in_bytes) bytes, $(_format_mb(size_in_bytes)) MB$(compression_suffix)"
end

"""
    _caid_total_size_label(entries)

Return the total-size label used in a CAID raw-data folder README.
"""
function _caid_total_size_label(entries)::String
    all(entry -> endswith(lowercase(String(entry.filename)), ".gz"), entries) ||
        return "total"
    "compressed, total"
end

"""
    _write_caid_folder_readme(folder_path, challenge, timestamp_utc, entries)

Write `README.md` metadata for a CAID raw-data folder.
"""
function _write_caid_folder_readme(
        folder_path::AbstractString,
        challenge::AbstractString,
        timestamp_utc::AbstractString,
        entries
)::Nothing
    readme_path = joinpath(folder_path, CAID_FOLDER_README_FILENAME)
    total_bytes = sum(filesize(entry.destination) for entry in entries)
    content = """
# $(challenge) Raw Reference Data

- Download timestamp (UTC): `$(timestamp_utc)`
- Final size ($(_caid_total_size_label(entries))): `$(_format_mb(total_bytes)) MB`

## Files
"""

    for entry in entries
        size_in_bytes = filesize(entry.destination)
        content *= """
- `$(entry.filename)` ($(_caid_entry_size_text(entry.filename, size_in_bytes)))
  - Source: `$(entry.url)`
"""
    end

    write(readme_path, content)
    return nothing
end

"""
    _build_caid_folder_entries(folder_path, source_urls_by_filename)

Collect README entries for files present in one CAID raw-data folder using a
known filename-to-source mapping.
"""
function _build_caid_folder_entries(
        folder_path::AbstractString,
        source_urls_by_filename::AbstractDict
)
    entries = NamedTuple{
        (:filename, :url, :destination),
        Tuple{String, String, String}
    }[]

    for filename in sort!(String.(collect(keys(source_urls_by_filename))))
        destination = joinpath(folder_path, filename)
        isfile(destination) || continue

        push!(
            entries,
            (
                filename = filename,
                url = String(source_urls_by_filename[filename]),
                destination = destination
            )
        )
    end

    entries
end

"""
    _caid3_readme_source_urls(; invalid_regions_filename, invalid_regions_url)

Build one filename-to-source mapping used to render the CAID-3 raw-data README.
"""
function _caid3_readme_source_urls(
        ;
        invalid_regions_filename::AbstractString = CAID_3_INVALID_REGIONS_FILENAME,
        invalid_regions_url::AbstractString = CAID_3_INVALID_REGIONS_URL
)::Dict{String, String}
    source_urls = Dict{String, String}()

    for filename in sort!(collect(keys(CAID_3_SOURCES)))
        source_urls[string(filename, ".gz")] = CAID_3_SOURCES[filename]
    end

    source_urls[String(invalid_regions_filename)] = String(invalid_regions_url)
    source_urls
end

"""
    _download_caid_dataset(dataset, raw_data_dir, timestamp_utc; logger = println)

Download one CAID dataset into its folder, validate raw FASTA headers, gzip files,
and write a per-folder `README.md`.
"""
function _download_caid_dataset(
        dataset,
        raw_data_dir::AbstractString,
        timestamp_utc::AbstractString;
        logger::Function = println
)
    folder_path = joinpath(raw_data_dir, dataset.folder)
    logger("Preparing folder $(folder_path)")
    clean_directory(folder_path)

    entries = NamedTuple{
        (:filename, :url, :destination),
        Tuple{String, String, String}
    }[]

    for filename in sort(collect(keys(dataset.sources)))
        url = dataset.sources[filename]
        destination = joinpath(folder_path, string(filename, ".gz"))
        logger("Downloading $(url)")
        download_to_gzip_file_with_retries(
            url,
            destination;
            delays = ExponentialBackOff(n = 5),
            logger = logger,
            precompress_validator = _caid_fasta_prefix_validator
        )

        push!(
            entries,
            (
                filename = basename(destination),
                url = url,
                destination = destination
            )
        )
    end

    _write_caid_folder_readme(folder_path, dataset.challenge, timestamp_utc, entries)
    return (; folder = dataset.folder, challenge = dataset.challenge,
        path = folder_path, entries = entries)
end

"""
    download_caid_reference_data(
        raw_data_dir;
        timestamp_utc = utc_now_iso8601(),
        logger = println,
        datasets = CAID_DATASETS
    )

Download CAID-2 and CAID-3 reference FASTA files into `raw_data_dir`, regenerate
their folder READMEs.
"""
function download_caid_reference_data(
        raw_data_dir::AbstractString;
        timestamp_utc::AbstractString = utc_now_iso8601(),
        logger::Function = println,
        datasets = CAID_DATASETS
)::Nothing
    mkpath(raw_data_dir)
    logger("Downloading CAID reference data into $(raw_data_dir)")

    for dataset in datasets
        logger("Processing $(dataset.challenge)")
        _download_caid_dataset(dataset, raw_data_dir, timestamp_utc; logger = logger)
    end

    logger("CAID reference download completed.")
    return nothing
end

"""
    download_caid3_invalid_regions_data(
        output_dir;
        filename = CAID_3_INVALID_REGIONS_FILENAME,
        source_url = CAID_3_INVALID_REGIONS_URL,
        timestamp_utc = utc_now_iso8601(),
        logger = println,
        downloader = Downloads.download
    )

Download the CAID reference invalid-region curation TSV into `output_dir`,
validate the downloaded file, and regenerate the CAID-3 raw-data `README.md`
for all recognized files currently present in that folder.
"""
function download_caid3_invalid_regions_data(
        output_dir::AbstractString;
        filename::AbstractString = CAID_3_INVALID_REGIONS_FILENAME,
        source_url::AbstractString = CAID_3_INVALID_REGIONS_URL,
        timestamp_utc::AbstractString = utc_now_iso8601(),
        logger::Function = println,
        downloader::Function = Downloads.download
)::Nothing
    mkpath(output_dir)
    destination = joinpath(output_dir, filename)

    logger("Downloading $(source_url)")
    remove_files([destination])

    try
        download_to_file_with_retries(
            source_url,
            destination;
            delays = ExponentialBackOff(n = 5),
            downloader = downloader,
            logger = logger
        )
        _validate_caid3_invalid_regions_tsv(destination)

        entries = _build_caid_folder_entries(
            output_dir,
            _caid3_readme_source_urls(
                ;
                invalid_regions_filename = filename,
                invalid_regions_url = source_url
            )
        )
        _write_caid_folder_readme(output_dir, "CAID-3", timestamp_utc, entries)
    catch
        remove_files([destination])
        rethrow()
    end

    logger("CAID-3 invalid-region TSV download completed.")
    return nothing
end
