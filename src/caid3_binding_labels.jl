using CodecZlib
using JSON

const _CAID3_BINDING_IDR_RAW_LABEL_TO_DATASET_LABEL = Dict{Char, Int}(
    '-' => 0,
    '0' => 1,
    '1' => 2
)
const _CAID3_BINDING_IDR_LABEL_ORDER = ['-', '0', '1']

"""
    _read_gzip_text(path)

Read a gzip-compressed text file and return its full contents as a string.
"""
function _read_gzip_text(path::AbstractString)::String
    isfile(path) || throw(ArgumentError("File does not exist: $(path)"))
    open(path) do io
        read(GzipDecompressorStream(io), String)
    end
end

"""
    _normalize_binding_idr_protein_id(raw_identifier)

Normalize one CAID-3 `binding_idr` record identifier to the first
whitespace-delimited token.
"""
function _normalize_binding_idr_protein_id(raw_identifier::AbstractString)::String
    identifier_text = strip(String(raw_identifier))
    isempty(identifier_text) && throw(
        ArgumentError("Encountered empty protein identifier in CAID-3 binding IDR record."),
    )

    token = split(identifier_text; limit = 2)[1]
    isempty(token) &&
        throw(ArgumentError("Encountered empty protein identifier in CAID-3 binding IDR record."))

    token
end

"""
    _read_caid3_binding_idr_records(binding_idr_fasta_gz_path)

Read CAID-3 `binding_idr.fasta.gz` records and return a dictionary keyed by
DisProt protein ID with values `(sequence, raw_labels)`.

This parser uses `MIToS.MSA.AnnotatedFASTASequences` and expects CAID-style
residue labels in the `"Feature"` annotation line (`-`, `0`, `1`).
"""
function _read_caid3_binding_idr_records(
        binding_idr_fasta_gz_path::AbstractString
)::Dict{String, NamedTuple{(:sequence, :raw_labels), Tuple{String, String}}}
    isfile(binding_idr_fasta_gz_path) || throw(
        ArgumentError("File does not exist: $(binding_idr_fasta_gz_path)"),
    )

    parsed_records = try
        MIToS.MSA.read_file(
            binding_idr_fasta_gz_path,
            MIToS.MSA.AnnotatedFASTASequences
        )
    catch err
        throw(
            ArgumentError(
            "Failed to parse CAID-3 binding IDR FASTA with MIToS.MSA.AnnotatedFASTASequences: $(binding_idr_fasta_gz_path). Underlying error: $(sprint(showerror, err))."
        ),
        )
    end
    isempty(parsed_records) &&
        throw(ArgumentError("CAID-3 binding IDR FASTA is empty: $(binding_idr_fasta_gz_path)"))

    records = Dict{String, NamedTuple{(:sequence, :raw_labels), Tuple{String, String}}}()
    for parsed_record in parsed_records
        parsed_sequence_id = MIToS.MSA.sequence_id(parsed_record)
        original_sequence_id = MIToS.MSA.getannotsequence(parsed_record, "OriginalSeqName", "")
        !isempty(original_sequence_id) && throw(
            ArgumentError(
            "Duplicate protein identifier in CAID-3 binding IDR FASTA: $(original_sequence_id)",
        ),
        )

        protein_id = _normalize_binding_idr_protein_id(parsed_sequence_id)
        haskey(records, protein_id) && throw(
            ArgumentError("Duplicate protein identifier in CAID-3 binding IDR FASTA: $(protein_id)"),
        )

        sequence = join(parsed_record)
        raw_labels = MIToS.MSA.getannotresidue(parsed_record, "Feature", "")
        isempty(sequence) && throw(
            ArgumentError("Protein $(protein_id): sequence is empty in $(binding_idr_fasta_gz_path)."),
        )
        isempty(raw_labels) && throw(
            ArgumentError(
            "Protein $(protein_id): residue annotation feature \"Feature\" is missing or empty in $(binding_idr_fasta_gz_path)."
        ),
        )
        length(sequence) == length(raw_labels) || throw(
            ArgumentError(
            "Protein $(protein_id): sequence length $(length(sequence)) does not match raw label length $(length(raw_labels)) in $(binding_idr_fasta_gz_path)."
        ),
        )

        records[protein_id] = (sequence = sequence, raw_labels = raw_labels)
    end

    records
end

"""
    _map_caid3_binding_raw_labels(raw_labels; protein_id)

Map one raw CAID-3 binding-IDR label string to final dataset residue labels:
- `-` -> `0`
- `0` -> `1`
- `1` -> `2`
"""
function _map_caid3_binding_raw_labels(
        raw_labels::AbstractString;
        protein_id::AbstractString
)::Vector{Int}
    mapped_labels = Vector{Int}(undef, length(raw_labels))
    for (index, raw_label) in enumerate(raw_labels)
        mapped_label = get(_CAID3_BINDING_IDR_RAW_LABEL_TO_DATASET_LABEL, raw_label, nothing)
        mapped_label === nothing && throw(
            ArgumentError(
            "Protein $(protein_id): unexpected CAID-3 raw label $(repr(raw_label)) at position $(index). Allowed labels: $(join(_CAID3_BINDING_IDR_LABEL_ORDER, ", "))."
        ),
        )
        mapped_labels[index] = mapped_label
    end

    mapped_labels
end

"""
    _read_dataset_records_by_disprot_id(dataset_json_gz_path)

Read final dataset JSON records and index them by `protein_disprot_id`.
Records without a non-empty `protein_disprot_id` are skipped.
"""
function _read_dataset_records_by_disprot_id(
        dataset_json_gz_path::AbstractString
)::Dict{String, Dict{String, Any}}
    payload = JSON.parse(_read_gzip_text(dataset_json_gz_path))
    payload isa AbstractVector || throw(
        ArgumentError(
        "Final dataset JSON must be an array of records: $(dataset_json_gz_path)"
    ),
    )

    records_by_disprot_id = Dict{String, Dict{String, Any}}()
    for (record_index, record_any) in enumerate(payload)
        record_any isa AbstractDict || throw(
            ArgumentError(
            "Record $(record_index) in $(dataset_json_gz_path) must be a dictionary-like object."
        ),
        )

        record = Dict{String, Any}(pairs(record_any))
        disprot_id_value = get(record, "protein_disprot_id", missing)
        disprot_id_value isa AbstractString || continue

        disprot_id = strip(String(disprot_id_value))
        isempty(disprot_id) && continue

        haskey(records_by_disprot_id, disprot_id) && throw(
            ArgumentError(
            "Duplicate protein_disprot_id $(disprot_id) found in $(dataset_json_gz_path)."
        ),
        )

        records_by_disprot_id[disprot_id] = record
    end

    records_by_disprot_id
end

"""
    _extract_dataset_residue_labels(record; protein_id)

Extract and validate `residue_labels` from one dataset JSON record.
"""
function _extract_dataset_residue_labels(
        record::AbstractDict;
        protein_id::AbstractString
)::Vector{Int}
    residue_labels_value = get(record, "residue_labels", missing)
    residue_labels_value isa AbstractVector || throw(
        ArgumentError(
        "Protein $(protein_id): dataset record is missing vector field `residue_labels`."
    ),
    )

    residue_labels = Vector{Int}(undef, length(residue_labels_value))
    for (index, value) in enumerate(residue_labels_value)
        value isa Integer || throw(
            ArgumentError(
            "Protein $(protein_id): residue_labels[$(index)] must be an integer, found $(repr(value))."
        ),
        )
        normalized_value = Int(value)
        normalized_value in (0, 1, 2) || throw(
            ArgumentError(
            "Protein $(protein_id): residue_labels[$(index)] must be one of 0, 1, 2, found $(normalized_value)."
        ),
        )
        residue_labels[index] = normalized_value
    end

    residue_labels
end

"""
    _format_mismatch_sample(expected_labels, observed_labels, mismatch_positions; max_examples)

Render mismatch examples as `position:expected!=observed` for reporting.
"""
function _format_mismatch_sample(
        expected_labels::Vector{Int},
        observed_labels::Vector{Int},
        mismatch_positions::Vector{Int};
        max_examples::Int
)::String
    max_examples >= 1 || throw(ArgumentError("`max_examples` must be >= 1."))

    examples = String[]
    for position in mismatch_positions[1:min(length(mismatch_positions), max_examples)]
        push!(
            examples,
            "$(position):$(expected_labels[position])!=$(observed_labels[position])"
        )
    end

    if length(mismatch_positions) > max_examples
        push!(examples, "...")
    end

    join(examples, ", ")
end

"""
    check_caid3_binding_idr_labels(
        binding_idr_fasta_gz_path,
        dataset_json_gz_path;
        max_mismatch_examples = 10
    )

Check CAID-3 `binding_idr` raw labels against final dataset `residue_labels`
for proteins present in `binding_idr.fasta.gz`.

Expected mapping:
- `-` -> `0`
- `0` -> `1`
- `1` -> `2`

Returns a named tuple with summary counts and mismatch details.
"""
function check_caid3_binding_idr_labels(
        binding_idr_fasta_gz_path::AbstractString,
        dataset_json_gz_path::AbstractString;
        max_mismatch_examples::Int = 10
)
    max_mismatch_examples >= 1 || throw(
        ArgumentError("`max_mismatch_examples` must be >= 1."),
    )

    binding_records = _read_caid3_binding_idr_records(binding_idr_fasta_gz_path)
    dataset_records_by_disprot_id = _read_dataset_records_by_disprot_id(dataset_json_gz_path)

    missing_proteins = String[]
    length_mismatches = NamedTuple{
        (:protein_id, :expected_length, :observed_length), Tuple{String, Int, Int}}[]
    label_mismatches = NamedTuple{
        (:protein_id, :mismatch_count, :sample_positions), Tuple{String, Int, String}}[]

    residues_compared = 0
    residues_mismatched = 0
    proteins_found = 0

    for protein_id in sort(collect(keys(binding_records)))
        expected_labels = _map_caid3_binding_raw_labels(
            binding_records[protein_id].raw_labels;
            protein_id = protein_id
        )

        if !haskey(dataset_records_by_disprot_id, protein_id)
            push!(missing_proteins, protein_id)
            continue
        end

        proteins_found += 1
        observed_labels = _extract_dataset_residue_labels(
            dataset_records_by_disprot_id[protein_id];
            protein_id = protein_id
        )

        if length(observed_labels) != length(expected_labels)
            push!(
                length_mismatches,
                (
                    protein_id = protein_id,
                    expected_length = length(expected_labels),
                    observed_length = length(observed_labels)
                )
            )
            continue
        end

        residues_compared += length(expected_labels)
        mismatch_positions = Int[]
        for index in eachindex(expected_labels, observed_labels)
            expected_labels[index] == observed_labels[index] ||
                push!(mismatch_positions, index)
        end

        residues_mismatched += length(mismatch_positions)
        isempty(mismatch_positions) || push!(
            label_mismatches,
            (
                protein_id = protein_id,
                mismatch_count = length(mismatch_positions),
                sample_positions = _format_mismatch_sample(
                    expected_labels,
                    observed_labels,
                    mismatch_positions;
                    max_examples = max_mismatch_examples
                )
            )
        )
    end

    proteins_passed = proteins_found - length(length_mismatches) - length(label_mismatches)
    passed = isempty(missing_proteins) &&
             isempty(length_mismatches) &&
             isempty(label_mismatches)

    (
        n_binding_idr_proteins = length(binding_records),
        n_dataset_records_with_disprot_id = length(dataset_records_by_disprot_id),
        n_proteins_found = proteins_found,
        n_proteins_passed = proteins_passed,
        missing_proteins = missing_proteins,
        length_mismatches = length_mismatches,
        label_mismatches = label_mismatches,
        residues_compared = residues_compared,
        residues_mismatched = residues_mismatched,
        passed = passed
    )
end

"""
    _render_missing_proteins_section(missing_proteins)

Render a markdown section listing proteins missing in the final dataset.
"""
function _render_missing_proteins_section(missing_proteins::Vector{String})::String
    isempty(missing_proteins) && return "## Missing Proteins\n\nNone."

    lines = ["- `$(protein_id)`" for protein_id in missing_proteins]
    "## Missing Proteins\n\n$(join(lines, "\n"))"
end

"""
    _render_length_mismatches_section(length_mismatches)

Render a markdown table for proteins with residue-label length mismatches.
"""
function _render_length_mismatches_section(length_mismatches)::String
    isempty(length_mismatches) && return "## Length Mismatches\n\nNone."

    table_lines = [
        "| protein_disprot_id | expected_length | observed_length |",
        "|---|---:|---:|"
    ]

    for mismatch in length_mismatches
        push!(
            table_lines,
            "| $(mismatch.protein_id) | $(mismatch.expected_length) | $(mismatch.observed_length) |"
        )
    end

    "## Length Mismatches\n\n$(join(table_lines, "\n"))"
end

"""
    _render_label_mismatches_section(label_mismatches)

Render a markdown table for proteins with residue-label value mismatches.
"""
function _render_label_mismatches_section(label_mismatches)::String
    isempty(label_mismatches) && return "## Residue Label Mismatches\n\nNone."

    table_lines = [
        "| protein_disprot_id | mismatched_positions | sample (position:expected!=observed) |",
        "|---|---:|---|"
    ]

    for mismatch in label_mismatches
        push!(
            table_lines,
            "| $(mismatch.protein_id) | $(mismatch.mismatch_count) | $(mismatch.sample_positions) |"
        )
    end

    "## Residue Label Mismatches\n\n$(join(table_lines, "\n"))"
end

"""
    write_caid3_binding_labels_report(
        output_path;
        binding_idr_fasta_gz_path,
        dataset_json_gz_path,
        check_result,
        script_relpath
    )

Write a markdown report for CAID-3 binding-IDR label consistency checks.
"""
function write_caid3_binding_labels_report(
        output_path::AbstractString;
        binding_idr_fasta_gz_path::AbstractString,
        dataset_json_gz_path::AbstractString,
        check_result,
        script_relpath::AbstractString
)::Nothing
    report = """
# CAID-3 Binding IDR Label Check

- Generated (UTC): `$(utc_now_iso8601())`
- Script: `$(script_relpath)`
- CAID-3 `binding_idr` FASTA: `$(binding_idr_fasta_gz_path)`
- Final dataset JSON: `$(dataset_json_gz_path)`

## Checked Mapping

- `-` -> `0`
- `0` -> `1`
- `1` -> `2`

## Summary

- Proteins in CAID-3 `binding_idr`: `$(check_result.n_binding_idr_proteins)`
- Final dataset records with `protein_disprot_id`: `$(check_result.n_dataset_records_with_disprot_id)`
- CAID-3 proteins found in final dataset: `$(check_result.n_proteins_found)`
- CAID-3 proteins passing all checks: `$(check_result.n_proteins_passed)`
- Missing proteins: `$(length(check_result.missing_proteins))`
- Length mismatches: `$(length(check_result.length_mismatches))`
- Residue-label mismatches: `$(length(check_result.label_mismatches))`
- Residues compared: `$(check_result.residues_compared)`
- Mismatched residues: `$(check_result.residues_mismatched)`
- Overall result: `$(check_result.passed ? "PASS" : "FAIL")`

$(_render_missing_proteins_section(check_result.missing_proteins))

$(_render_length_mismatches_section(check_result.length_mismatches))

$(_render_label_mismatches_section(check_result.label_mismatches))
"""

    _write_text_file(output_path, report)
    nothing
end
