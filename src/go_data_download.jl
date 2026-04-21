using OBOParse

const GO_FOLDER_NAME = "go"
const GO_FOLDER_README_FILENAME = "README.md"
const GO_BASIC_FILENAME = "go-basic.obo.gz"
const GO_BASIC_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

"""
    _validate_go_obo_file(path; parser_loader = OBOParse.load)

Validate that one downloaded GO OBO file exists, is not empty, and can be
parsed by `OBOParse`.
"""
function _validate_go_obo_file(
        path::AbstractString;
        parser_loader::Function = OBOParse.load
)::Nothing
    isfile(path) || error("Downloaded GO OBO file does not exist: $(path)")
    filesize(path) > 0 || error("Downloaded GO OBO file is empty: $(path)")

    ontology = parser_loader(path, "GO")
    hasproperty(ontology, :terms) || error(
        "Parsed GO OBO object does not expose `terms`: $(typeof(ontology)).",
    )
    terms = getproperty(ontology, :terms)
    length(terms) > 0 || error("Parsed GO OBO ontology has no terms: $(path)")
    nothing
end

"""
    _write_go_folder_readme(folder_path, timestamp_utc, filename, source_url, destination)

Write `README.md` metadata for the GO raw-data folder.
"""
function _write_go_folder_readme(
        folder_path::AbstractString,
        timestamp_utc::AbstractString,
        filename::AbstractString,
        source_url::AbstractString,
        destination::AbstractString
)::Nothing
    readme_path = joinpath(folder_path, GO_FOLDER_README_FILENAME)
    size_in_bytes = filesize(destination)
    content = """
# GO Raw Ontology Data

- Download timestamp (UTC): `$(timestamp_utc)`
- Final size (compressed, total): `$(_format_mb(size_in_bytes)) MB`
- Source URL: `$(source_url)`
- Validation: `OBOParse.load(path, "GO")` on downloaded uncompressed content before gzip compression

## Files
- `$(filename)` ($(size_in_bytes) bytes, $(_format_mb(size_in_bytes)) MB compressed)
"""
    write(readme_path, content)
    nothing
end

"""
    download_go_reference_data(
        raw_data_dir;
        filename = GO_BASIC_FILENAME,
        source_url = GO_BASIC_URL,
        timestamp_utc = utc_now_iso8601(),
        logger = println,
        downloader = Downloads.download,
        parser_loader = OBOParse.load
    )

Download `go-basic.obo` into `joinpath(raw_data_dir, "go", "go-basic.obo.gz")`,
validate uncompressed content with `OBOParse` before compression, and regenerate
the folder README.
"""
function download_go_reference_data(
        raw_data_dir::AbstractString;
        filename::AbstractString = GO_BASIC_FILENAME,
        source_url::AbstractString = GO_BASIC_URL,
        timestamp_utc::AbstractString = utc_now_iso8601(),
        logger::Function = println,
        downloader::Function = Downloads.download,
        parser_loader::Function = OBOParse.load
)::Nothing
    mkpath(raw_data_dir)
    go_dir = joinpath(raw_data_dir, GO_FOLDER_NAME)
    logger("Downloading GO ontology data into $(go_dir)")
    clean_directory(go_dir)

    destination = joinpath(go_dir, filename)
    logger("Downloading $(source_url)")

    try
        download_to_gzip_file_with_retries(
            source_url,
            destination;
            delays = ExponentialBackOff(n = 5),
            downloader = downloader,
            logger = logger,
            precompress_validator = path -> _validate_go_obo_file(
                path;
                parser_loader = parser_loader
            )
        )
        _write_go_folder_readme(
            go_dir,
            timestamp_utc,
            filename,
            source_url,
            destination
        )
    catch
        remove_files([destination])
        rethrow()
    end

    logger("GO ontology download completed.")
    nothing
end
