#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset: download_sifts_reference_data

const RAW_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data", "raw"))

function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download SIFTS summary CSV mappings for PDB-chain to UniProt and CATH into data/raw/sifts.

Run from the repository root:
  julia scripts/03_download_sifts.jl --raw_data_dir data/raw

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--raw_data_dir"
        help = "Raw data directory where the sifts folder will be created."
        arg_type = String
        default = RAW_DATA_DIR
    end

    return normpath(parse_args(args, settings)["raw_data_dir"])
end

download_sifts_reference_data(parse_command_line(ARGS))
