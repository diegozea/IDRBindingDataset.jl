using DataFrames
using Statistics

const DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC = ["false positive", "true negative"]

const _REQUIRED_ELM_INSTANCE_COLUMNS = [
    :ELMIdentifier,
    :Primary_Acc,
    :Accessions,
    :Start,
    :End,
    :InstanceLogic
]
const _REQUIRED_ELM_INDEX_COLUMNS = [:ELMIdentifier, :Regex]
const _REQUIRED_CLUSTER_REPRESENTATIVE_COLUMNS = [:protein_id]
const _REQUIRED_ELM_MOTIF_OUTPUT_COLUMNS = [:protein_id, :regex, :start, :end]

"""
    _strip_wrapping_quotes(value)

Trim one string value and remove wrapping single or double quote characters.
"""
function _strip_wrapping_quotes(value::AbstractString)::String
    stripped = strip(String(value))
    strip(stripped, ['"', '\''])
end

"""
    _normalize_logic_values(values)

Normalize `InstanceLogic` filter values by trimming whitespace, removing wrapping
quotes, dropping empty values, and preserving insertion order.
"""
function _normalize_logic_values(values::AbstractVector)::Vector{String}
    normalized = String[]
    seen = Set{String}()

    for raw_value in values
        cleaned_value = _strip_wrapping_quotes(String(raw_value))
        isempty(cleaned_value) && continue
        cleaned_value in seen && continue

        push!(seen, cleaned_value)
        push!(normalized, cleaned_value)
    end

    normalized
end

"""
    _parse_excluded_instance_logic(option_value)

Parse one comma-separated command-line value for `InstanceLogic` exclusions.
Wrapping quote characters are ignored.
"""
function _parse_excluded_instance_logic(option_value::AbstractString)::Vector{String}
    _normalize_logic_values(split(String(option_value), ","))
end

"""
    _coerce_position(value, column_name)

Coerce one motif coordinate value to `Int`.
"""
function _coerce_position(value, column_name::AbstractString)::Int
    value isa Integer && return Int(value)

    text_value = _strip_wrapping_quotes(String(value))
    parsed_value = tryparse(Int, text_value)
    parsed_value === nothing && throw(
        ArgumentError(
        "Invalid $(column_name) value in ELM instances table: $(repr(value)). Expected integer coordinate."
    ),
    )

    parsed_value
end

"""
    _split_accessions(primary_acc, accessions)

Collect one deterministic, de-duplicated accession list from `Primary_Acc` and
the whitespace/comma/semicolon-separated `Accessions` field.
"""
function _split_accessions(primary_acc, accessions)::Vector{String}
    primary_text = _strip_wrapping_quotes(String(coalesce(primary_acc, "")))
    accessions_text = String(coalesce(accessions, ""))

    accession_values = String[]
    seen = Set{String}()

    for accession in vcat([primary_text], split(accessions_text, r"[\s,;]+"))
        cleaned_accession = _strip_wrapping_quotes(accession)
        isempty(cleaned_accession) && continue
        cleaned_accession in seen && continue

        push!(seen, cleaned_accession)
        push!(accession_values, cleaned_accession)
    end

    accession_values
end

"""
    _filter_elm_instances_by_logic(elm_instances; excluded_instance_logic)

Return a copy of `elm_instances` with rows removed when `InstanceLogic` matches
one of the normalized `excluded_instance_logic` values.
"""
function _filter_elm_instances_by_logic(
        elm_instances::AbstractDataFrame;
        excluded_instance_logic::AbstractVector = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC
)::DataFrame
    _require_columns(elm_instances, _REQUIRED_ELM_INSTANCE_COLUMNS, "elm_instances")

    normalized_excluded_values = Set(_normalize_logic_values(excluded_instance_logic))
    isempty(normalized_excluded_values) && return DataFrame(elm_instances)

    filtered_rows = filter(
        row -> _strip_wrapping_quotes(String(row.InstanceLogic)) ∉
               normalized_excluded_values,
        elm_instances
    )

    DataFrame(filtered_rows)
end

"""
    _expand_elm_instances_by_accessions(elm_instances)

Expand ELM instances into one row per accession token mapped to one motif
instance. The returned table includes:
- `protein_id`
- `ELMIdentifier`
- `start`
- `end`
"""
function _expand_elm_instances_by_accessions(elm_instances::AbstractDataFrame)::DataFrame
    _require_columns(elm_instances, _REQUIRED_ELM_INSTANCE_COLUMNS, "elm_instances")

    protein_ids = String[]
    elm_identifiers = String[]
    starts = Int[]
    ends = Int[]

    for row in eachrow(elm_instances)
        elm_identifier = _strip_wrapping_quotes(String(row.ELMIdentifier))
        isempty(elm_identifier) && throw(
            ArgumentError("Found ELM instance with empty ELMIdentifier."),
        )

        start_position = _coerce_position(row.Start, "Start")
        end_position = _coerce_position(row.End, "End")

        for accession in _split_accessions(row.Primary_Acc, row.Accessions)
            push!(protein_ids, accession)
            push!(elm_identifiers, elm_identifier)
            push!(starts, start_position)
            push!(ends, end_position)
        end
    end

    expanded_instances = DataFrame(
        protein_id = protein_ids,
        ELMIdentifier = elm_identifiers,
        start = starts
    )
    expanded_instances[!, :end] = ends

    unique!(expanded_instances)
    sort!(expanded_instances, [:protein_id, :ELMIdentifier, :start, :end])

    expanded_instances
end

"""
    _elm_identifier_to_regex(elm_index)

Build one deterministic mapping from `ELMIdentifier` to `Regex` from one ELM
classes index table. Throws when one identifier maps to conflicting regexes.
"""
function _elm_identifier_to_regex(elm_index::AbstractDataFrame)::Dict{String, String}
    _require_columns(elm_index, _REQUIRED_ELM_INDEX_COLUMNS, "elm_index")

    mapping = Dict{String, String}()
    for row in eachrow(elm_index)
        elm_identifier = _strip_wrapping_quotes(String(row.ELMIdentifier))
        regex_pattern = _strip_wrapping_quotes(String(row.Regex))

        isempty(elm_identifier) && throw(
            ArgumentError("Found ELM class with empty ELMIdentifier."),
        )
        isempty(regex_pattern) && throw(
            ArgumentError("Found ELM class with empty Regex for ELMIdentifier $(elm_identifier)."),
        )

        if haskey(mapping, elm_identifier) && mapping[elm_identifier] != regex_pattern
            throw(
                ArgumentError(
                "Conflicting Regex values for ELMIdentifier $(elm_identifier): " *
                "$(repr(mapping[elm_identifier])) vs $(repr(regex_pattern))."
            ),
            )
        end

        mapping[elm_identifier] = regex_pattern
    end

    mapping
end

"""
    build_cluster_elm_motifs_table(
        cluster_representatives,
        elm_instances,
        elm_index;
        excluded_instance_logic = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC
    )

Build a deterministic table with motif regex mappings for proteins present in
`cluster_representatives`.

Rules:
- remove `elm_instances` rows with excluded `InstanceLogic` values,
- map proteins using `Primary_Acc` and fallback accession tokens in `Accessions`,
- keep only proteins present in `cluster_representatives`,
- join motif regex patterns by `ELMIdentifier` from `elm_index`,
- return unique rows sorted by `protein_id`, `start`, `end`, and `regex`.

The returned columns are:
- `protein_id`
- `regex`
- `start`
- `end`
"""
function build_cluster_elm_motifs_table(
        cluster_representatives::AbstractDataFrame,
        elm_instances::AbstractDataFrame,
        elm_index::AbstractDataFrame;
        excluded_instance_logic::AbstractVector = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC
)::DataFrame
    _require_columns(
        cluster_representatives,
        _REQUIRED_CLUSTER_REPRESENTATIVE_COLUMNS,
        "cluster_representatives"
    )

    filtered_instances = _filter_elm_instances_by_logic(
        elm_instances;
        excluded_instance_logic = excluded_instance_logic
    )
    expanded_instances = _expand_elm_instances_by_accessions(filtered_instances)
    regex_by_identifier = _elm_identifier_to_regex(elm_index)

    motifs_with_regex = copy(expanded_instances)
    motifs_with_regex[!, :regex] = Union{Missing, String}[get(regex_by_identifier, elm_identifier, missing)
                                                          for elm_identifier in
                                                              motifs_with_regex.ELMIdentifier]

    missing_regex_identifiers = String[String(identifier)
                                       for
                                       identifier in
                                       motifs_with_regex[ismissing.(motifs_with_regex.regex), :ELMIdentifier]]
    unique!(missing_regex_identifiers)
    sort!(missing_regex_identifiers)
    if !isempty(missing_regex_identifiers)
        missing_regex_text = join(missing_regex_identifiers, ", ")
        throw(
            ArgumentError(
            "ELM identifiers in instances were not found in the classes index: $(missing_regex_text)",
        ),
        )
    end

    disallowmissing!(motifs_with_regex, :regex)

    cluster_proteins = unique(select(cluster_representatives, :protein_id))
    transform!(
        cluster_proteins,
        :protein_id => ByRow(value -> _strip_wrapping_quotes(String(value))) => :protein_id
    )

    mapped_motifs = innerjoin(cluster_proteins, motifs_with_regex; on = :protein_id)
    motif_table = select(mapped_motifs, :protein_id, :regex, :start, :end)

    unique!(motif_table)
    sort!(motif_table, [:protein_id, :start, :end, :regex])

    motif_table
end

"""
    write_elm_motifs_csv_gz(motif_table, output_path)

Write one gzip-compressed CSV artifact for mapped ELM motifs.
"""
function write_elm_motifs_csv_gz(
        motif_table::AbstractDataFrame,
        output_path::AbstractString
)::Nothing
    _require_columns(motif_table, _REQUIRED_ELM_MOTIF_OUTPUT_COLUMNS, "motif_table")

    _write_csv_gz(motif_table, output_path)
end

"""
    _summarize_elm_motif_table(motif_table, cluster_protein_count)

Compute summary statistics for one mapped ELM motif table.
"""
function _summarize_elm_motif_table(
        motif_table::AbstractDataFrame,
        cluster_protein_count::Int
)::NamedTuple
    _require_columns(motif_table, _REQUIRED_ELM_MOTIF_OUTPUT_COLUMNS, "motif_table")

    proteins_with_motifs = length(unique(motif_table.protein_id))
    motif_rows = nrow(motif_table)
    unique_regex_count = length(unique(motif_table.regex))

    motifs_per_protein = isempty(motif_table) ? Int[] :
                         combine(groupby(motif_table, :protein_id), nrow => :count).count

    average_motifs_per_positive_protein = isempty(motifs_per_protein) ? 0.0 :
                                          mean(motifs_per_protein)
    median_motifs_per_positive_protein = isempty(motifs_per_protein) ? 0.0 :
                                         median(Float64.(motifs_per_protein))
    max_motifs_per_protein = isempty(motifs_per_protein) ? 0 : maximum(motifs_per_protein)

    average_motifs_per_cluster_protein = cluster_protein_count == 0 ? 0.0 :
                                         motif_rows / cluster_protein_count
    proteins_without_motifs = max(cluster_protein_count - proteins_with_motifs, 0)
    proteins_with_motifs_percentage = cluster_protein_count == 0 ? 0.0 :
                                      proteins_with_motifs / cluster_protein_count * 100.0

    (
        motif_rows = motif_rows,
        unique_regex_count = unique_regex_count,
        proteins_with_motifs = proteins_with_motifs,
        proteins_without_motifs = proteins_without_motifs,
        proteins_with_motifs_percentage = proteins_with_motifs_percentage,
        average_motifs_per_positive_protein = average_motifs_per_positive_protein,
        median_motifs_per_positive_protein = median_motifs_per_positive_protein,
        max_motifs_per_protein = max_motifs_per_protein,
        average_motifs_per_cluster_protein = average_motifs_per_cluster_protein
    )
end

"""
    write_elm_motifs_readme(
        output_path;
        output_csv_path,
        output_csv_filesystem_path = output_csv_path,
        source_cluster_representatives_path,
        source_elm_instances_path,
        source_elms_index_path,
        script_relpath = "scripts/10_mapping_elm_motifs.jl",
        excluded_instance_logic = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC,
        input_instance_rows,
        filtered_instance_rows,
        motif_table,
        cluster_protein_count
    )

Write a README for the mapped ELM motif artifact.
"""
function write_elm_motifs_readme(
        output_path::AbstractString;
        output_csv_path::AbstractString,
        output_csv_filesystem_path::AbstractString = output_csv_path,
        source_cluster_representatives_path::AbstractString,
        source_elm_instances_path::AbstractString,
        source_elms_index_path::AbstractString,
        script_relpath::AbstractString = "scripts/10_mapping_elm_motifs.jl",
        excluded_instance_logic::AbstractVector = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC,
        input_instance_rows::Int,
        filtered_instance_rows::Int,
        motif_table::AbstractDataFrame,
        cluster_protein_count::Int
)::Nothing
    _require_columns(motif_table, _REQUIRED_ELM_MOTIF_OUTPUT_COLUMNS, "motif_table")
    isfile(output_csv_filesystem_path) || throw(
        ArgumentError("Motif CSV file does not exist for README metadata: $(output_csv_filesystem_path)"),
    )

    excluded_rows = max(input_instance_rows - filtered_instance_rows, 0)
    normalized_excluded_logic = _normalize_logic_values(excluded_instance_logic)
    normalized_excluded_logic_text = isempty(normalized_excluded_logic) ? "(none)" :
                                     join(
        ["`$(value)`"
         for value in normalized_excluded_logic], ", ")

    summary = _summarize_elm_motif_table(motif_table, cluster_protein_count)
    output_csv_size_bytes = filesize(output_csv_filesystem_path)

    readme_text = """
# ELM Motif Mapping

This folder stores regex-level ELM motif mappings for proteins in the
cluster-representative table.

## Artifact

- Compressed CSV: `$(output_csv_path)` (`$(output_csv_size_bytes)` bytes)
- Row count: `$(summary.motif_rows)`
- Unique proteins with at least one mapped motif: `$(summary.proteins_with_motifs)`
- Unique regex patterns represented: `$(summary.unique_regex_count)`

## Data Dictionary

- `protein_id`: UniProt accession from cluster representatives.
- `regex`: ELM regex pattern (from `elms_index.tsv`) mapped through `ELMIdentifier`.
- `start`: ELM instance start coordinate (1-based, inclusive).
- `end`: ELM instance end coordinate (1-based, inclusive).

## Mapping Decisions

- `InstanceLogic` filtering excludes: $(normalized_excluded_logic_text).
- Protein mapping first uses `Primary_Acc` and then falls back to parsed tokens
  from `Accessions` when needed.
- Only proteins present in `cluster_representatives.csv.gz` are retained.
- Regex assignments are joined by `ELMIdentifier` from `elms_index.tsv`.
- Duplicate output rows by `(protein_id, regex, start, end)` are removed.

## Summary Statistics

- Cluster proteins evaluated: `$(cluster_protein_count)`
- Proteins with motifs: `$(summary.proteins_with_motifs)` (`$(round(summary.proteins_with_motifs_percentage; digits = 2))%`)
- Proteins without motifs: `$(summary.proteins_without_motifs)`
- Average motifs per motif-positive protein: `$(round(summary.average_motifs_per_positive_protein; digits = 3))`
- Median motifs per motif-positive protein: `$(round(summary.median_motifs_per_positive_protein; digits = 3))`
- Maximum motifs for one protein: `$(summary.max_motifs_per_protein)`
- Average motifs per cluster protein (including proteins without motifs): `$(round(summary.average_motifs_per_cluster_protein; digits = 3))`

## Input Filtering

- ELM instance rows before filtering: `$(input_instance_rows)`
- ELM instance rows after filtering: `$(filtered_instance_rows)`
- ELM instance rows removed by `InstanceLogic`: `$(excluded_rows)`

## Source Inputs

- Cluster representatives: `$(source_cluster_representatives_path)`
- ELM instances: `$(source_elm_instances_path)`
- ELM classes index: `$(source_elms_index_path)`

## Reproducibility

Run from repository root:

```bash
julia $(script_relpath)
```

Complete example:

```bash
julia $(script_relpath) --instances_path data/scratch/elm/instances.tsv --elms_index_path data/scratch/elm/elms_index.tsv --cluster_representatives_path data/derived/clustering/cluster_representatives.csv.gz --exclude_instance_logic "false positive,true negative"
```
"""

    _write_text_file(output_path, readme_text)
end
