#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse

using IDRBindingDataset:
                         build_disprot_regions_url,
                         remove_files,
                         download_to_gzip_file_with_retries,
                         ensure_files_within_size_limit,
                         _resolve_latest_idpo_release,
                         _validate_idpo_obojson_payload,
                         utc_now_iso8601,
                         write_download_report

const RELEASE = "2025_12"

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))
const DISPROT_DATA_DIR = joinpath(RAW_DATA_DIR, "disprot")
const RAW_README_FILENAME = "README.md"

const DISPROT_REGIONS_FILENAME = "disprot_regions.json.gz"
const IDPO_FILENAME = "idpo.json.gz"
const IDPO_MAPPINGS_FILENAME = "mapping_ontology.tsv.gz"

const IDPO_MAPPINGS_URL = "https://disprot.org/assets/data/mapping_ontology_Disprot-DO_GO_IDPO.tsv"
const IDPO_LATEST_RELEASE_API_URL = "https://api.github.com/repos/BioComputingUP/idpo/releases/latest"
# GitHub blocks files larger than 100 MiB in regular Git repositories.
const GITHUB_SINGLE_FILE_HARD_LIMIT_BYTES = 100 * 1024 * 1024

"""
    parse_command_line(args)

Parse command-line arguments and return selected versions.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download DisProt regions, IDPO ontology, and IDPO mappings into data/raw/disprot.

Run from the repository root:
  julia scripts/01_download_disprot_data.jl --release 2025_12

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--release"
        help = "DisProt release identifier (default: $(RELEASE))."
        arg_type = String
        default = RELEASE
    end

    parsed = parse_args(args, settings)
    return (; release = parsed["release"])
end

function main()
    options = parse_command_line(ARGS)
    release = options.release
    download_timestamp = utc_now_iso8601()
    report_path = joinpath(DISPROT_DATA_DIR, RAW_README_FILENAME)
    idpo_release_tag = "unresolved"
    idpo_url = IDPO_LATEST_RELEASE_API_URL
    sources = (
        (build_disprot_regions_url(release), DISPROT_REGIONS_FILENAME),
        (idpo_url, IDPO_FILENAME),
        (IDPO_MAPPINGS_URL, IDPO_MAPPINGS_FILENAME)
    )
    downloaded_paths = [joinpath(DISPROT_DATA_DIR, filename) for (_, filename) in sources]

    mkpath(DISPROT_DATA_DIR)
    println("Downloading DisProt sources to $(DISPROT_DATA_DIR)")
    println("Using release = $(release)")

    try
        idpo_release = _resolve_latest_idpo_release()
        idpo_release_tag = idpo_release.release_tag
        idpo_url = idpo_release.asset_url

        sources = (
            (build_disprot_regions_url(release), DISPROT_REGIONS_FILENAME),
            (idpo_url, IDPO_FILENAME),
            (IDPO_MAPPINGS_URL, IDPO_MAPPINGS_FILENAME)
        )
        downloaded_paths = [joinpath(DISPROT_DATA_DIR, filename)
                            for (_, filename) in sources]

        println("Using IDPO release tag = $(idpo_release_tag)")
        println("Using IDPO source = $(idpo_url)")
        println("Removing previously downloaded raw data files before download")
        remove_files(vcat(downloaded_paths, [report_path]))

        for (url, filename) in sources
            destination = joinpath(DISPROT_DATA_DIR, filename)
            println("Downloading $(url)")
            download_to_gzip_file_with_retries(
                url,
                destination;
                delays = ExponentialBackOff(n = 5),
                precompress_validator = filename == IDPO_FILENAME ?
                                        _validate_idpo_obojson_payload : nothing
            )
        end

        ensure_files_within_size_limit(
            downloaded_paths;
            max_bytes = GITHUB_SINGLE_FILE_HARD_LIMIT_BYTES,
            limit_name = "GitHub single-file hard limit",
            remediation = "GitHub blocks files larger than 100 MiB. Use Git LFS or store large artifacts in release assets/object storage."
        )

        write_download_report(
            report_path;
            status = "SUCCESS",
            timestamp_utc = download_timestamp,
            release = release,
            idpo_release_tag = idpo_release_tag,
            sources = sources
        )
    catch err
        write_download_report(
            report_path;
            status = "FAILED",
            timestamp_utc = download_timestamp,
            release = release,
            idpo_release_tag = idpo_release_tag,
            sources = sources,
            error_message = sprint(showerror, err)
        )
        rethrow()
    end

    println("All downloads completed.")
end

main()
