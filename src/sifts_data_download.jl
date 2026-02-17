const SIFTS_FOLDER_NAME = "sifts"
const SIFTS_FOLDER_README_FILENAME = "README.md"
const SIFTS_SUMMARY_DOC_URL = "https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html"
const SIFTS_DATABASES = (
    MIToS.SIFTS.dbUniProt,
    MIToS.SIFTS.dbCATH
)

"""
    _sifts_mapping_description(database)

Return a human-readable mapping description for one SIFTS summary `database`.
"""
function _sifts_mapping_description(database::Type{<:MIToS.SIFTS.DataBase})::String
    database == MIToS.SIFTS.dbUniProt && return "PDB chain to UniProt mapping"
    database == MIToS.SIFTS.dbCATH &&
        return "PDB chain to CATH mapping with UniProt cross-reference"
    error("Unsupported SIFTS summary database type: $(database)")
end

"""
    _write_sifts_folder_readme(folder_path, timestamp_utc, entries)

Write `README.md` metadata for a SIFTS raw-data folder.
"""
function _write_sifts_folder_readme(
        folder_path::AbstractString,
        timestamp_utc::AbstractString,
        entries
)::Nothing
    readme_path = joinpath(folder_path, SIFTS_FOLDER_README_FILENAME)
    total_bytes = sum(filesize(entry.destination) for entry in entries)
    content = """
# SIFTS Chain-Level Mapping Data

- Download timestamp (UTC): `$(timestamp_utc)`
- Final size (compressed, total): `$(_format_mb(total_bytes)) MB`
- Downloaded with: `MIToS.SIFTS.downloadsifts`
- Source documentation: `$(SIFTS_SUMMARY_DOC_URL)`

## Files
"""

    for entry in entries
        size_in_bytes = filesize(entry.destination)
        content *= """
- `$(entry.filename)` ($(size_in_bytes) bytes, $(_format_mb(size_in_bytes)) MB compressed)
  - Mapping: `$(entry.mapping)`
"""
    end

    write(readme_path, content)
    return nothing
end

"""
    _download_sifts_summary_file(database, folder_path; logger = println)

Download one SIFTS summary file into `folder_path` and return its metadata.
"""
function _download_sifts_summary_file(
        database::Type{<:MIToS.SIFTS.DataBase},
        folder_path::AbstractString;
        logger::Function = println
)
    logger("Downloading SIFTS summary for $(database)")
    downloaded_path = cd(folder_path) do
        MIToS.SIFTS.downloadsifts(database)
    end
    destination = isabspath(downloaded_path) ?
                  downloaded_path :
                  normpath(joinpath(folder_path, downloaded_path))
    filename = basename(destination)
    endswith(filename, ".gz") ||
        error("SIFTS summary download did not produce a compressed file: $(filename)")
    isfile(destination) || error("SIFTS download failed to produce file: $(destination)")
    return (
        filename = filename,
        mapping = _sifts_mapping_description(database),
        destination = destination
    )
end

"""
    download_sifts_reference_data(
        raw_data_dir;
        timestamp_utc = utc_now_iso8601(),
        logger = println
    )

Download SIFTS chain-level summary CSV files for UniProt and CATH into
`joinpath(raw_data_dir, "sifts")`, preserving gzip compression and regenerating
the folder README.
"""
function download_sifts_reference_data(
        raw_data_dir::AbstractString;
        timestamp_utc::AbstractString = utc_now_iso8601(),
        logger::Function = println
)::Nothing
    mkpath(raw_data_dir)
    sifts_dir = joinpath(raw_data_dir, SIFTS_FOLDER_NAME)
    logger("Downloading SIFTS summary data into $(sifts_dir)")
    clean_directory(sifts_dir)

    entries = NamedTuple{
        (:filename, :mapping, :destination),
        Tuple{String, String, String}
    }[]

    for database in SIFTS_DATABASES
        push!(
            entries,
            _download_sifts_summary_file(
                database,
                sifts_dir;
                logger = logger
            )
        )
    end

    _write_sifts_folder_readme(sifts_dir, timestamp_utc, entries)
    logger("SIFTS summary download completed.")
    return nothing
end
