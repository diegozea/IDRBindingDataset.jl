const CAID_RAW_README_FILENAME = "README.md"
const CAID_FOLDER_README_FILENAME = "README.md"
const CAID_REFERENCE_SECTION_HEADING = "## CAID Reference Folders"

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
- Final size (compressed, total): `$(_format_mb(total_bytes)) MB`

## Files
"""

    for entry in entries
        size_in_bytes = filesize(entry.destination)
        content *= """
- `$(entry.filename)` ($(size_in_bytes) bytes, $(_format_mb(size_in_bytes)) MB compressed)
  - Source: `$(entry.url)`
"""
    end

    write(readme_path, content)
    return nothing
end

"""
    _build_caid_folders_section(timestamp_utc)

Build the markdown section appended to `data/raw/README.md` describing CAID folders.
"""
function _build_caid_folders_section(timestamp_utc::AbstractString)::String
    return """
## CAID Reference Folders

- `caid_2/`: CAID-2 reference files (`*.fasta.gz`) and per-folder metadata.
- `caid_3/`: CAID-3 reference files (`*.fasta.gz`) and per-folder metadata.
- Last CAID refresh (UTC): `$(timestamp_utc)`
"""
end

"""
    upsert_caid_folders_section(content, timestamp_utc)

Insert or replace the `## CAID Reference Folders` section in `content`.
"""
function upsert_caid_folders_section(
        content::AbstractString,
        timestamp_utc::AbstractString
)::String
    section = _build_caid_folders_section(timestamp_utc)
    section_lines = split(chomp(section), '\n'; keepempty = true)
    content_lines = split(content, '\n'; keepempty = true)

    section_start = findfirst(
        line -> strip(line) == CAID_REFERENCE_SECTION_HEADING,
        content_lines
    )

    if section_start === nothing
        trimmed_content = rstrip(content)
        if isempty(trimmed_content)
            return string(section, "\n")
        end
        return string(trimmed_content, "\n\n", section, "\n")
    end

    section_end = length(content_lines) + 1
    for line_index in (section_start + 1):length(content_lines)
        if startswith(strip(content_lines[line_index]), "## ")
            section_end = line_index
            break
        end
    end

    prefix = content_lines[1:(section_start - 1)]
    suffix = content_lines[section_end:end]
    updated_lines = vcat(prefix, section_lines, suffix)
    updated_content = join(updated_lines, "\n")
    return endswith(updated_content, "\n") ? updated_content : string(updated_content, "\n")
end

"""
    append_caid_mention_to_raw_readme(raw_data_dir, timestamp_utc)

Append or refresh the CAID section in `raw_data_dir/README.md` if it exists.
"""
function append_caid_mention_to_raw_readme(
        raw_data_dir::AbstractString,
        timestamp_utc::AbstractString
)::Nothing
    raw_readme_path = joinpath(raw_data_dir, CAID_RAW_README_FILENAME)
    if !isfile(raw_readme_path)
        @warn "README.md not found at $(raw_readme_path). Skipping CAID folder mention update."
        return nothing
    end

    existing_content = read(raw_readme_path, String)
    updated_content = upsert_caid_folders_section(existing_content, timestamp_utc)
    write(raw_readme_path, updated_content)
    return nothing
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
    download_caid_reference_data(raw_data_dir; timestamp_utc = utc_now_iso8601(), logger = println)

Download CAID-2 and CAID-3 reference FASTA files into `raw_data_dir`, regenerate
their folder READMEs, and update `raw_data_dir/README.md` with a CAID section.
"""
function download_caid_reference_data(
        raw_data_dir::AbstractString;
        timestamp_utc::AbstractString = utc_now_iso8601(),
        logger::Function = println
)::Nothing
    mkpath(raw_data_dir)
    logger("Downloading CAID reference data into $(raw_data_dir)")

    for dataset in CAID_DATASETS
        logger("Processing $(dataset.challenge)")
        _download_caid_dataset(dataset, raw_data_dir, timestamp_utc; logger = logger)
    end

    append_caid_mention_to_raw_readme(raw_data_dir, timestamp_utc)
    logger("CAID reference download completed.")
    return nothing
end
