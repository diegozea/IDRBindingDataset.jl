using CodecZlib
using DataFrames
using JSON

const _RESIDUE_LABEL_ORDERED = 0
const _RESIDUE_LABEL_DISORDERED = 1
const _RESIDUE_LABEL_BINDING = 2

"""
    _build_disprot_entry_index(disprot_entries)

Build one accession-indexed dictionary from DisProt entries with normalized
sequence and region payload.
"""
function _build_disprot_entry_index(disprot_entries)::Dict{
        String,
        NamedTuple{(:sequence, :regions), Tuple{String, Any}}
}
    entry_index = Dict{String, NamedTuple{(:sequence, :regions), Tuple{String, Any}}}()

    for (entry_number, entry) in enumerate(disprot_entries)
        entry isa AbstractDict || throw(
            ArgumentError("DisProt entry $(entry_number) must be a dictionary-like record."),
        )

        protein_id = _normalize_required_protein_id(
            get(entry, "acc", missing);
            argument_name = "disprot_entries[$(entry_number)].acc"
        )
        haskey(entry_index, protein_id) && throw(
            ArgumentError("Duplicate DisProt accession detected: $(protein_id)."),
        )

        sequence_value = get(entry, "sequence", missing)
        sequence_value isa AbstractString || throw(
            ArgumentError(
            "DisProt entry $(entry_number) ($(protein_id)) is missing a string `sequence` field."
        ),
        )
        sequence = _normalize_disprot_sequence(sequence_value)

        regions = get(entry, "regions", Any[])
        regions isa AbstractVector || throw(
            ArgumentError(
            "DisProt entry $(entry_number) ($(protein_id)) must provide `regions` as an array."
        ),
        )

        entry_index[protein_id] = (sequence = sequence, regions = regions)
    end

    entry_index
end

"""
    _region_has_disorder_annotation(region, term_rules)

Return `true` when one region term ID belongs to the ontology-derived IDPO
disorder challenge set.
"""
function _region_has_disorder_annotation(
        region,
        term_rules::_DisProtChallengeTermRules
)::Bool
    term_id = _region_term_id(region)
    term_id === nothing && return false
    startswith(term_id, "IDPO:") || return false
    term_id in term_rules.idpo_disorder_terms
end

"""
    _region_has_binding_annotation(region, term_rules)

Return `true` when one region term ID belongs to the ontology-derived GO
binding challenge set.
"""
function _region_has_binding_annotation(
        region,
        term_rules::_DisProtChallengeTermRules
)::Bool
    _region_has_idr_binding_annotation(region, term_rules)
end

"""
    _region_bounds(region, sequence_length; protein_id, region_label)

Normalize one region `(start, end)` to a validated 1-based inclusive range.
"""
function _region_bounds(
        region::AbstractDict,
        sequence_length::Integer;
        protein_id::AbstractString,
        region_label::AbstractString
)::UnitRange{Int}
    start_value = get(region, "start", missing)
    end_value = get(region, "end", missing)
    start_value isa Integer || throw(
        ArgumentError(
        "Protein $(protein_id): $(region_label) region start must be an integer, found $(repr(start_value))."
    ),
    )
    end_value isa Integer || throw(
        ArgumentError(
        "Protein $(protein_id): $(region_label) region end must be an integer, found $(repr(end_value))."
    ),
    )

    start_pos = Int(start_value)
    end_pos = Int(end_value)
    start_pos <= end_pos || throw(
        ArgumentError(
        "Protein $(protein_id): $(region_label) region start $(start_pos) is greater than end $(end_pos)."
    ),
    )
    start_pos >= 1 || throw(
        ArgumentError(
        "Protein $(protein_id): $(region_label) region start $(start_pos) is out of bounds for sequence length $(sequence_length)."
    ),
    )
    end_pos <= sequence_length || throw(
        ArgumentError(
        "Protein $(protein_id): $(region_label) region end $(end_pos) is out of bounds for sequence length $(sequence_length)."
    ),
    )

    start_pos:end_pos
end

"""
    _paint_region_labels!(residue_labels, region; target_label, protein_id, region_label)

Assign one label over a validated region range. Existing labels are preserved
when they are greater than `target_label` (for binding precedence).
"""
function _paint_region_labels!(
        residue_labels::Vector{Int},
        region::AbstractDict;
        target_label::Int,
        protein_id::AbstractString,
        region_label::AbstractString
)::Nothing
    bounds = _region_bounds(
        region,
        length(residue_labels);
        protein_id = protein_id,
        region_label = region_label
    )
    for position in bounds
        # Use monotonic precedence: ordered (0) -> disorder (1) -> binding (2),
        # so once a residue reaches 2 it is never downgraded by later regions.
        residue_labels[position] = max(residue_labels[position], target_label)
    end

    nothing
end

"""
    _build_residue_labels(sequence, regions; protein_id)

Build one residue-label vector:
- `0`: ordered
- `1`: disordered
- `2`: binding (overrides disorder)
"""
function _build_residue_labels(
        sequence::AbstractString,
        regions;
        protein_id::AbstractString,
        term_rules::_DisProtChallengeTermRules,
        invalid_region_ids::Union{Nothing, Set{String}} = nothing
)::Vector{Int}
    residue_labels = fill(_RESIDUE_LABEL_ORDERED, length(sequence))

    for region in regions
        region isa AbstractDict || continue
        invalid_region_ids !== nothing &&
            _region_is_marked_invalid(region, invalid_region_ids) && continue

        has_disorder_annotation = _region_has_disorder_annotation(region, term_rules)
        if has_disorder_annotation
            _paint_region_labels!(
                residue_labels,
                region;
                target_label = _RESIDUE_LABEL_DISORDERED,
                protein_id = protein_id,
                region_label = "disorder"
            )
        end

        has_binding_annotation = _region_has_binding_annotation(region, term_rules)
        if has_binding_annotation
            _paint_region_labels!(
                residue_labels,
                region;
                target_label = _RESIDUE_LABEL_BINDING,
                protein_id = protein_id,
                region_label = "binding"
            )
        end
    end

    residue_labels
end

"""
    build_dataset_json_records(
        final_dataset_table,
        disprot_entries;
        term_rules,
        invalid_region_ids = nothing
    )

Build JSON-ready records from one final dataset table by preserving all original
columns and adding:
- `sequence`
- `residue_labels`

Region classification uses ontology-derived GO and IDPO term sets from
`term_rules`.

When `invalid_region_ids` is provided, regions whose `region_id` is present in
that set are skipped.
"""
function build_dataset_json_records(
        final_dataset_table::AbstractDataFrame,
        disprot_entries;
        term_rules::_DisProtChallengeTermRules,
        invalid_region_ids::Union{Nothing, Set{String}} = nothing
)::Vector{NamedTuple}
    _require_columns(final_dataset_table, [:protein_id], "final_dataset_table")

    dataset_columns = collect(propertynames(final_dataset_table))
    output_fields = Tuple(vcat(dataset_columns, [:sequence, :residue_labels]))
    disprot_index = _build_disprot_entry_index(disprot_entries)
    records = NamedTuple[]

    for row in eachrow(final_dataset_table)
        protein_id = _normalize_required_protein_id(
            row.protein_id;
            argument_name = "final_dataset_table.protein_id"
        )
        haskey(disprot_index, protein_id) || throw(
            ArgumentError(
            "Protein $(protein_id) from final_dataset_table is missing in DisProt entries."
        ),
        )

        disprot_entry = disprot_index[protein_id]
        sequence = disprot_entry.sequence
        residue_labels = _build_residue_labels(
            sequence,
            disprot_entry.regions;
            protein_id = protein_id,
            term_rules = term_rules,
            invalid_region_ids = invalid_region_ids
        )
        record_values = tuple(
            (row[column_name] for column_name in dataset_columns)...,
            sequence,
            residue_labels
        )
        push!(records, NamedTuple{output_fields}(record_values))
    end

    records
end

"""
    _residue_label_to_char(value)::Char

Convert one residue label value (`0`, `1`, `2`) to its one-character string
code.
"""
function _residue_label_to_char(value)::Char
    value isa Integer || throw(
        ArgumentError(
        "Residue label values must be integers `0`, `1`, or `2`; found $(repr(value))."
    ),
    )

    label = Int(value)
    label in (_RESIDUE_LABEL_ORDERED, _RESIDUE_LABEL_DISORDERED, _RESIDUE_LABEL_BINDING) ||
        throw(
            ArgumentError(
            "Residue label values must be integers `0`, `1`, or `2`; found $(repr(value))."
        ),
        )

    Char(Int('0') + label)
end

"""
    _encode_residue_labels(labels)::String

Encode one per-residue integer-label vector as a compact string with one
character per residue (`0`, `1`, `2`).
"""
function _encode_residue_labels(labels)::String
    labels isa AbstractVector || throw(
        ArgumentError("`residue_labels` must be a vector-like value."),
    )

    encoded = Vector{Char}(undef, length(labels))
    for (index, label) in enumerate(labels)
        encoded[index] = _residue_label_to_char(label)
    end

    String(encoded)
end

"""
    build_dataset_csv_table(json_records)::DataFrame

Build one CSV-ready table from JSON-ready records. All record fields are
preserved, and `residue_labels` is encoded from integer vectors to compact
strings with one character per residue (`0`, `1`, `2`).
"""
function build_dataset_csv_table(json_records)::DataFrame
    table = DataFrame(json_records)
    _require_columns(table, [:residue_labels], "json_records")

    transform!(
        table,
        :residue_labels => ByRow(_encode_residue_labels) => :residue_labels
    )

    table
end

"""
    _json_ready_value(value)

Recursively convert Julia values to JSON-friendly values, including conversion
of `missing` to `nothing` (`null` in JSON).
"""
function _json_ready_value(value)
    value
end

"""
    _json_ready_value(value::Missing)

Convert `missing` to `nothing` so JSON output encodes as `null`.
"""
function _json_ready_value(::Missing)
    nothing
end

"""
    _json_ready_value(value::NamedTuple)

Convert one named tuple recursively while preserving field names.
"""
function _json_ready_value(value::NamedTuple)
    (; (name => _json_ready_value(getfield(value, name)) for name in keys(value))...)
end

"""
    _json_ready_value(value::AbstractVector)

Convert one vector recursively for JSON serialization.
"""
function _json_ready_value(value::AbstractVector)
    [_json_ready_value(element) for element in value]
end

"""
    _json_ready_value(value::AbstractDict)

Convert one dictionary recursively for JSON serialization using string keys.
"""
function _json_ready_value(value::AbstractDict)
    Dict(string(key) => _json_ready_value(element) for (key, element) in pairs(value))
end

"""
    write_dataset_json_gz(records, output_path)

Write one gzip-compressed JSON array containing `records`.
"""
function write_dataset_json_gz(records, output_path::AbstractString)::Nothing
    mkpath(dirname(output_path))
    json_ready_records = _json_ready_value(records)
    open(GzipCompressorStream, output_path, "w") do io
        write(io, JSON.json(json_ready_records))
        write(io, '\n')
    end

    nothing
end
