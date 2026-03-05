using CodecZlib
using Dates
using JSON

const _DISPROT_DISORDER_TERM_ID = "IDPO:0000002"
const _DISPROT_DISORDER_TRANSITION_TERM_ID_DISORDER_TO_ORDER = "IDPO:0000011"
const _DISPROT_DISORDER_TRANSITION_TERM_ID_DISORDER_TO_DISORDER = "IDPO:0000017"
const _DISPROT_DISORDER_TRANSITION_TERM_ID_ORDER_TO_DISORDER = "IDPO:0000014"
const _DISPROT_DISORDER_TERM_IDS = Set([
    _DISPROT_DISORDER_TERM_ID,
    _DISPROT_DISORDER_TRANSITION_TERM_ID_DISORDER_TO_ORDER,
    _DISPROT_DISORDER_TRANSITION_TERM_ID_DISORDER_TO_DISORDER,
    _DISPROT_DISORDER_TRANSITION_TERM_ID_ORDER_TO_DISORDER
])
const _DISPROT_DISORDER_LABEL_PATTERN = r"\bdisorder\b"i
const _RESIDUE_LABEL_ORDERED = 0
const _RESIDUE_LABEL_DISORDERED = 1
const _RESIDUE_LABEL_BINDING = 2
const _JSON_README_SECTION_HEADER = "## JSON Artifact"
const _JSON_DEFAULT_SCRIPT_RELPATH = "scripts/14_json_format.jl"

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
    _region_has_disorder_annotation(region)

Return `true` when one region is annotated as disorder (`IDPO:0000002`),
one disorder-transition term (`IDPO:0000011`, `IDPO:0000017`,
`IDPO:0000014`), or, for compatibility, when `term_name` includes the token
`disorder` (case-insensitive).
"""
function _region_has_disorder_annotation(region)::Bool
    region isa AbstractDict || return false

    term_id = get(region, "term_id", missing)
    if term_id isa AbstractString && strip(String(term_id)) in _DISPROT_DISORDER_TERM_IDS
        return true
    end

    term_name = get(region, "term_name", missing)
    term_name isa AbstractString || return false
    occursin(_DISPROT_DISORDER_LABEL_PATTERN, String(term_name))
end

"""
    _region_has_binding_annotation(region)

Return `true` when one region is annotated as binding using the same logic as
`build_disprot_protein_annotations`.
"""
function _region_has_binding_annotation(region)::Bool
    _region_has_idr_binding_annotation(region)
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
        protein_id::AbstractString
)::Vector{Int}
    residue_labels = fill(_RESIDUE_LABEL_ORDERED, length(sequence))

    for region in regions
        region isa AbstractDict || continue

        if _region_has_disorder_annotation(region)
            _paint_region_labels!(
                residue_labels,
                region;
                target_label = _RESIDUE_LABEL_DISORDERED,
                protein_id = protein_id,
                region_label = "disorder"
            )
        end

        if _region_has_binding_annotation(region)
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
    build_dataset_json_records(final_dataset_table, disprot_entries)

Build JSON-ready records from one final dataset table by preserving all original
columns and adding:
- `sequence`
- `residue_labels`
"""
function build_dataset_json_records(
        final_dataset_table::AbstractDataFrame,
        disprot_entries
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
            protein_id = protein_id
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

"""
    _strip_json_readme_section(readme_text)

Remove one existing `## JSON Artifact` section from the last occurrence of the
section header to the end of the README text.
"""
function _strip_json_readme_section(readme_text::AbstractString)::String
    section_range = findlast(_JSON_README_SECTION_HEADER, readme_text)
    isnothing(section_range) && return rstrip(readme_text)

    section_start = first(section_range)
    section_start == firstindex(readme_text) && return ""

    prefix_end = prevind(readme_text, section_start)
    rstrip(readme_text[firstindex(readme_text):prefix_end])
end

"""
    _render_dataset_json_readme_section(;
        output_json_path,
        output_json_filesystem_path,
        source_dataset_csv_path,
        source_disprot_regions_path,
        n_proteins,
        generated_at_utc,
        script_relpath
    )

Render one deterministic README markdown section describing the JSON artifact.
"""
function _render_dataset_json_readme_section(;
        output_json_path::AbstractString,
        output_json_filesystem_path::AbstractString,
        source_dataset_csv_path::AbstractString,
        source_disprot_regions_path::AbstractString,
        n_proteins::Integer,
        generated_at_utc::DateTime,
        script_relpath::AbstractString
)::String
    n_proteins >= 0 || throw(ArgumentError("n_proteins must be non-negative."))

    isfile(output_json_filesystem_path) || throw(
        ArgumentError(
        "Final dataset JSON file does not exist for README metadata: $(output_json_filesystem_path)",
    ),
    )

    generated_at_text = _format_generated_at_utc(generated_at_utc)
    output_json_size_bytes = filesize(output_json_filesystem_path)
    output_json_size_mb = output_json_size_bytes / 1024^2

    """
$_JSON_README_SECTION_HEADER

- Compressed JSON: `$(output_json_path)` (`$(output_json_size_bytes)` bytes; `$(round(output_json_size_mb; digits = 3))` MiB)
- Generated at (UTC): `$(generated_at_text)`
- Number of proteins: `$(n_proteins)`

### Content

- Preserves all columns from `cluster_dataset_labels.csv.gz`.
- Adds `sequence`: full protein sequence from DisProt.
- Adds `residue_labels`: per-residue integer vector with `0 = ordered`,
  `1 = disordered`, and `2 = binding`.
- Label precedence is monotonic (`0 -> 1 -> 2`), so residues labeled as
  `2` are never downgraded by later overlapping regions.

### Source Inputs

- Final dataset CSV: `$(source_dataset_csv_path)`
- DisProt regions JSON: `$(source_disprot_regions_path)`

### Reproducibility

Run from repository root:

```bash
julia $(script_relpath)
```

Complete example:

```bash
julia $(script_relpath) --dataset_csv_path data/dataset/cluster_dataset_labels.csv.gz --raw_data_dir data/raw --output_dir data/dataset --output_filename cluster_dataset_labels.json.gz
```
"""
end

"""
    write_dataset_json_readme_section(
        readme_path;
        output_json_path,
        output_json_filesystem_path = output_json_path,
        source_dataset_csv_path,
        source_disprot_regions_path,
        n_proteins,
        generated_at_utc = now(UTC),
        script_relpath = "scripts/14_json_format.jl"
    )

Append one JSON artifact section at the end of a dataset README, replacing any
previous `## JSON Artifact` section written by this pipeline.
"""
function write_dataset_json_readme_section(
        readme_path::AbstractString;
        output_json_path::AbstractString,
        output_json_filesystem_path::AbstractString = output_json_path,
        source_dataset_csv_path::AbstractString,
        source_disprot_regions_path::AbstractString,
        n_proteins::Integer,
        generated_at_utc::DateTime = now(UTC),
        script_relpath::AbstractString = _JSON_DEFAULT_SCRIPT_RELPATH
)::Nothing
    json_section = _render_dataset_json_readme_section(
        output_json_path = output_json_path,
        output_json_filesystem_path = output_json_filesystem_path,
        source_dataset_csv_path = source_dataset_csv_path,
        source_disprot_regions_path = source_disprot_regions_path,
        n_proteins = n_proteins,
        generated_at_utc = generated_at_utc,
        script_relpath = script_relpath
    )

    existing_readme = isfile(readme_path) ? read(readme_path, String) : ""
    base_readme = _strip_json_readme_section(existing_readme)
    merged_readme = isempty(base_readme) ? json_section :
                    string(base_readme, "\n\n", json_section)

    _write_text_file(readme_path, string(merged_readme, "\n"))
end
