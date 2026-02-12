#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset: download_caid_reference_data

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))

"""
    parse_command_line(args)

Parse command-line arguments for CAID raw-data download.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download CAID-2 and CAID-3 annotated FASTA reference files into data/raw/caid_2 and data/raw/caid_3.

Run from the repository root:
  julia scripts/02_download_caid_data.jl --raw_data_dir data/raw

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--raw_data_dir"
        help = "Raw data directory where caid_2 and caid_3 folders will be created."
        arg_type = String
        default = RAW_DATA_DIR
    end

    parsed = parse_args(args, settings)
    return (; raw_data_dir = normpath(parsed["raw_data_dir"]))
end

function main()
    options = parse_command_line(ARGS)
    download_caid_reference_data(options.raw_data_dir)
    return nothing
end

main()
