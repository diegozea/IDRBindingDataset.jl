#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset: download_go_reference_data

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))

"""
    parse_command_line(args)

Parse command-line arguments for GO raw-data download.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download Gene Ontology `go-basic.obo` into data/raw/go and validate it using OBOParse.
The stored artifact is gzip-compressed as `go-basic.obo.gz`.

Run from the repository root:
  julia scripts/15_download_go_data.jl --raw_data_dir data/raw

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--raw_data_dir"
        help = "Raw data directory where the go folder will be created."
        arg_type = String
        default = RAW_DATA_DIR
    end

    parsed = parse_args(args, settings)
    return (; raw_data_dir = normpath(parsed["raw_data_dir"]))
end

function main()
    options = parse_command_line(ARGS)
    download_go_reference_data(options.raw_data_dir)
    nothing
end

main()
