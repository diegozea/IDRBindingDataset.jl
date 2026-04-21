#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using IDRBindingDataset:
                         _relpath_from_repo,
                         check_caid3_binding_idr_labels,
                         remove_files,
                         write_caid3_binding_labels_report

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_BINDING_IDR_FASTA_GZ_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "raw",
    "caid_3",
    "binding_idr.fasta.gz"
))
const DEFAULT_DATASET_JSON_GZ_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "dataset",
    "cluster_dataset_labels.json.gz"
))
const DEFAULT_OUTPUT_REPORT_PATH = normpath(joinpath(
    REPO_ROOT,
    "data",
    "dataset",
    "CAID3_LABELS.md"
))

"""
    parse_command_line(args)

Parse command-line options for the CAID-3 binding-IDR label consistency check.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Check CAID-3 binding-IDR residue labels against the final dataset JSON labels.

For proteins in `data/raw/caid_3/binding_idr.fasta.gz`, validate that final
`residue_labels` in `data/dataset/cluster_dataset_labels.json.gz` follow:
- `-` in CAID-3 -> `0` in final dataset
- `0` in CAID-3 -> `1` in final dataset
- `1` in CAID-3 -> `2` in final dataset

The script writes a markdown report summarizing all checks.

Run from repository root:
  julia scripts/16_check_idr_binding.jl

Complete example:
  julia scripts/16_check_idr_binding.jl --binding_idr_fasta_path data/raw/caid_3/binding_idr.fasta.gz --dataset_json_path data/dataset/cluster_dataset_labels.json.gz --output_report_path data/dataset/CAID3_LABELS.md --max_mismatch_examples 10

Optional non-blocking mode:
  julia scripts/16_check_idr_binding.jl --allow_mismatch

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--binding_idr_fasta_path"
        help = "Input CAID-3 binding-IDR FASTA gzip file."
        arg_type = String
        default = DEFAULT_BINDING_IDR_FASTA_GZ_PATH
        "--dataset_json_path"
        help = "Input final dataset JSON gzip file."
        arg_type = String
        default = DEFAULT_DATASET_JSON_GZ_PATH
        "--output_report_path"
        help = "Output markdown report path."
        arg_type = String
        default = DEFAULT_OUTPUT_REPORT_PATH
        "--max_mismatch_examples"
        help = "Maximum mismatch examples reported per protein."
        arg_type = Int
        default = 10
        "--fail_on_mismatch"
        help = "Exit with an error when mismatches are found. Enabled by default; this flag is kept for compatibility."
        action = :store_true
        "--allow_mismatch"
        help = "Do not exit with an error when mismatches are found (report-only mode)."
        action = :store_true
    end

    parsed = parse_args(args, settings)
    fail_on_mismatch = (
        !Bool(parsed["allow_mismatch"]) ||
        Bool(parsed["fail_on_mismatch"])
    )
    (
        binding_idr_fasta_path = normpath(parsed["binding_idr_fasta_path"]),
        dataset_json_path = normpath(parsed["dataset_json_path"]),
        output_report_path = normpath(parsed["output_report_path"]),
        max_mismatch_examples = Int(parsed["max_mismatch_examples"]),
        fail_on_mismatch = fail_on_mismatch
    )
end

function main()
    options = parse_command_line(ARGS)
    mkpath(dirname(options.output_report_path))
    remove_files([options.output_report_path])

    check_result = check_caid3_binding_idr_labels(
        options.binding_idr_fasta_path,
        options.dataset_json_path;
        max_mismatch_examples = options.max_mismatch_examples
    )

    write_caid3_binding_labels_report(
        options.output_report_path;
        binding_idr_fasta_gz_path = _relpath_from_repo(options.binding_idr_fasta_path),
        dataset_json_gz_path = _relpath_from_repo(options.dataset_json_path),
        check_result = check_result,
        script_relpath = _relpath_from_repo(@__FILE__)
    )

    println(
        "CAID-3 binding-IDR check report written to $(options.output_report_path) (result=$(check_result.passed ? "PASS" : "FAIL"), checked=$(check_result.n_proteins_found)/$(check_result.n_binding_idr_proteins), residue_mismatches=$(check_result.residues_mismatched))."
    )

    if options.fail_on_mismatch && !check_result.passed
        throw(ErrorException("CAID-3 binding-IDR label check failed. See $(options.output_report_path)."))
    end

    nothing
end

main()
