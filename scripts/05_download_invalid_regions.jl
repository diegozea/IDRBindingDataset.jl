#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset: download_caid3_invalid_regions_data

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(RAW_DATA_DIR, "caid_3"))
const DEFAULT_FILENAME = "invalid_regions_fullCAID3.tsv"
const DEFAULT_SOURCE_URL = "https://raw.githubusercontent.com/BioComputingUP/caid-reference/main/data/disprot/invalid_regions_fullCAID3.tsv"

"""
    parse_command_line(args)

Parse command-line options for CAID-3 invalid-region TSV download.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download the CAID-3 invalid-region curation TSV into data/raw/caid_3.

Run from the repository root:
  julia scripts/05_download_invalid_regions.jl --output_dir data/raw/caid_3 --filename invalid_regions_fullCAID3.tsv --source_url https://raw.githubusercontent.com/BioComputingUP/caid-reference/main/data/disprot/invalid_regions_fullCAID3.tsv

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--output_dir"
        help = "Directory where the invalid-region TSV and CAID-3 README will be written."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--filename"
        help = "Output filename for the invalid-region TSV."
        arg_type = String
        default = DEFAULT_FILENAME
        "--source_url"
        help = "Source URL for the CAID-3 invalid-region TSV."
        arg_type = String
        default = DEFAULT_SOURCE_URL
    end

    parsed = parse_args(args, settings)
    return (
        output_dir = normpath(parsed["output_dir"]),
        filename = String(parsed["filename"]),
        source_url = String(parsed["source_url"])
    )
end

function main()
    options = parse_command_line(ARGS)
    download_caid3_invalid_regions_data(
        options.output_dir;
        filename = options.filename,
        source_url = options.source_url
    )
    return nothing
end

main()
