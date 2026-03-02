using CodecZlib
using CSV
using DataFrames
using Dates

const FINAL_DATASET_LABEL_VALIDATION = "validation"
const _FINAL_DATASET_DEFAULT_SCRIPT_RELPATH = "scripts/13_assing_missing_labels.jl"
const _FINAL_DATASET_DEFAULT_LABEL_ORDER = [
    FIXED_DATASET_LABEL_TRAINING,
    FINAL_DATASET_LABEL_VALIDATION,
    _FIXED_DATASET_LABEL_TESTING
]
const _FINAL_DATASET_ALLOWED_LABELS = Set(
    _FINAL_DATASET_DEFAULT_LABEL_ORDER
)
const _FIXED_TO_FINAL_DATASET_LABEL_MAP = Dict(
    FIXED_DATASET_LABEL_TRAINING => FIXED_DATASET_LABEL_TRAINING,
    FIXED_DATASET_LABEL_CAID => _FIXED_DATASET_LABEL_TESTING,
    _FIXED_DATASET_LABEL_TESTING => _FIXED_DATASET_LABEL_TESTING,
    FINAL_DATASET_LABEL_VALIDATION => FINAL_DATASET_LABEL_VALIDATION
)

"""
    _normalize_elm_cluster_token(value)

Normalize one ELM motif-cluster token value to `Union{Nothing, String}`.
Returns `nothing` for missing or empty values.
"""
function _normalize_elm_cluster_token(value)::Union{Nothing, String}
    ismissing(value) && return nothing

    token = strip(string(value))
    isempty(token) && return nothing
    token
end

"""
    _sort_elm_cluster_tokens(tokens)

Sort one vector of ELM motif-cluster tokens deterministically. Numeric tokens
are sorted numerically when all tokens are valid integers; otherwise tokens are
sorted lexicographically.
"""
function _sort_elm_cluster_tokens(tokens::Vector{String})::Vector{String}
    parsed_tokens = tryparse.(Int, tokens)
    if all(!isnothing, parsed_tokens)
        return sort(tokens; by = token -> parse(Int, token))
    end

    return sort(tokens)
end

"""
    _join_unique_sorted_elm_clusters(values; separator)

Build one deterministic, separator-joined ELM motif-cluster string from
possibly repeated values, or return `missing` when no cluster token is present.
"""
function _join_unique_sorted_elm_clusters(
        values;
        separator::AbstractString
)::Union{Missing, String}
    tokens = String[]
    for value in values
        token = _normalize_elm_cluster_token(value)
        token === nothing && continue
        push!(tokens, token)
    end

    unique!(tokens)
    isempty(tokens) && return missing
    sorted_tokens = _sort_elm_cluster_tokens(tokens)
    join(sorted_tokens, separator)
end

"""
    _normalize_required_protein_id(value; argument_name = "protein_id")

Normalize one protein identifier to a stripped `String` and throw an
`ArgumentError` when the identifier is missing or empty.
"""
function _normalize_required_protein_id(
        value;
        argument_name::AbstractString = "protein_id"
)::String
    (ismissing(value) || isnothing(value)) && throw(
        ArgumentError("$(argument_name) must not contain missing values."),
    )

    protein_id = strip(string(value))
    isempty(protein_id) && throw(
        ArgumentError("$(argument_name) must not contain empty protein IDs."),
    )

    protein_id
end

"""
    build_protein_elm_clusters_table(
        protein_elm_motif_clusters;
        separator = "-",
        output_column = :protein_elm_clusters
    )

Collapse one protein-level ELM motif-cluster table to one row per protein with
joined unique cluster IDs.

Expected input columns:
- `protein_id`
- `elm_motif_cluster`

Returned columns:
- `protein_id`
- `output_column` (default `protein_elm_clusters`)
"""
function build_protein_elm_clusters_table(
        protein_elm_motif_clusters::AbstractDataFrame;
        separator::AbstractString = "-",
        output_column::Symbol = :protein_elm_clusters
)::DataFrame
    isempty(separator) && throw(ArgumentError("separator must not be empty."))
    _require_columns(
        protein_elm_motif_clusters,
        [:protein_id, :elm_motif_cluster],
        "protein_elm_motif_clusters"
    )

    normalized = select(protein_elm_motif_clusters, :protein_id, :elm_motif_cluster)
    transform!(
        normalized,
        :protein_id => ByRow(
            value -> _normalize_required_protein_id(
            value;
            argument_name = "protein_elm_motif_clusters.protein_id"
        )
        ) => :protein_id
    )

    protein_elm_clusters = combine(
        groupby(normalized, :protein_id),
        :elm_motif_cluster => (
            values -> _join_unique_sorted_elm_clusters(values; separator = separator)
        ) => output_column
    )

    sort!(protein_elm_clusters, :protein_id)
    protein_elm_clusters
end

"""
    _normalize_protein_elm_clusters_table(
        protein_elm_clusters;
        elm_cluster_column = :protein_elm_clusters
    )

Normalize one protein ELM cluster table to unique `protein_id` rows while
verifying duplicate rows carry the same ELM cluster assignment.
"""
function _normalize_protein_elm_clusters_table(
        protein_elm_clusters::AbstractDataFrame;
        elm_cluster_column::Symbol = :protein_elm_clusters
)::DataFrame
    _require_columns(
        protein_elm_clusters,
        [:protein_id, elm_cluster_column],
        "protein_elm_clusters"
    )

    normalized = select(
        protein_elm_clusters,
        :protein_id => ByRow(
            value -> _normalize_required_protein_id(
            value;
            argument_name = "protein_elm_clusters.protein_id"
        )
        ) => :protein_id,
        elm_cluster_column
    )

    normalized = combine(
        groupby(normalized, :protein_id),
        elm_cluster_column => (
            values -> _single_nonmissing_or_missing(
            values;
            value_label = string(elm_cluster_column)
        )
        ) => elm_cluster_column
    )
    sort!(normalized, :protein_id)
    normalized
end

"""
    _normalize_dataset_label_token(value)

Normalize one dataset-label token by trimming whitespace and lowercasing.
"""
function _normalize_dataset_label_token(value)::String
    lowercase(strip(string(value)))
end

"""
    _validate_final_dataset_label(value; argument_name)

Validate one final dataset label belongs to:
`training`, `validation`, `testing`.
"""
function _validate_final_dataset_label(
        value::AbstractString;
        argument_name::AbstractString
)::String
    normalized_label = _normalize_dataset_label_token(value)
    normalized_label in _FINAL_DATASET_ALLOWED_LABELS || throw(
        ArgumentError(
        "Invalid $(argument_name): $(repr(value)). Allowed labels are: $(join(sort!(collect(_FINAL_DATASET_ALLOWED_LABELS)), ", "))."
    ),
    )

    normalized_label
end

"""
    _map_fixed_to_final_dataset_label(value; missing_label)

Map one `cluster_fixed_dataset_label` value to the final
`cluster_dataset_label`:
- missing -> `missing_label` (`validation` by default),
- `training` -> `training`,
- `caid` -> `testing`,
- `testing` -> `testing`.
"""
function _map_fixed_to_final_dataset_label(
        value;
        missing_label::AbstractString
)::String
    ismissing(value) && return missing_label

    normalized_fixed_label = _normalize_dataset_label_token(value)
    isempty(normalized_fixed_label) && return missing_label
    haskey(_FIXED_TO_FINAL_DATASET_LABEL_MAP, normalized_fixed_label) || throw(
        ArgumentError(
        "Unsupported cluster_fixed_dataset_label value: $(repr(value)). Supported values are: $(join(sort!(collect(keys(_FIXED_TO_FINAL_DATASET_LABEL_MAP))), ", "))."
    ),
    )

    _FIXED_TO_FINAL_DATASET_LABEL_MAP[normalized_fixed_label]
end

"""
    _ensure_final_dataset_allowed_labels(final_dataset; label_column)

Ensure one final dataset label column only contains:
`training`, `validation`, `testing`.
"""
function _ensure_final_dataset_allowed_labels(
        final_dataset::AbstractDataFrame;
        label_column::Symbol
)::Nothing
    _require_columns(final_dataset, [label_column], "final_dataset")

    unique_labels = Set(_normalize_dataset_label_token.(final_dataset[!, label_column]))
    invalid_labels = setdiff(unique_labels, _FINAL_DATASET_ALLOWED_LABELS)
    isempty(invalid_labels) || throw(
        ArgumentError(
        "Column $(label_column) contains unsupported labels: $(join(sort!(collect(invalid_labels)), ", ")). Allowed labels are: $(join(sort!(collect(_FINAL_DATASET_ALLOWED_LABELS)), ", "))."
    ),
    )

    nothing
end

"""
    build_final_cluster_dataset_table(
        fixed_labels_table,
        protein_elm_clusters;
        missing_label = FINAL_DATASET_LABEL_VALIDATION,
        elm_cluster_column = :protein_elm_clusters
    )

Build the final per-protein dataset by joining ELM motif-cluster assignments and
mapping fixed cluster labels to a constrained final split:
- `training` stays `training`,
- `caid` becomes `testing`,
- missing becomes `missing_label` (`validation` by default).

Expected columns in `fixed_labels_table`:
- `protein_id`
- `final_cluster_id`
- `cluster_fixed_dataset_label`

Expected columns in `protein_elm_clusters`:
- `protein_id`
- `elm_cluster_column`

Adds one new column:
- `cluster_dataset_label`
"""
function build_final_cluster_dataset_table(
        fixed_labels_table::AbstractDataFrame,
        protein_elm_clusters::AbstractDataFrame;
        missing_label::AbstractString = FINAL_DATASET_LABEL_VALIDATION,
        elm_cluster_column::Symbol = :protein_elm_clusters
)::DataFrame
    cleaned_missing_label = _validate_final_dataset_label(
        missing_label;
        argument_name = "missing_label"
    )
    _require_columns(
        fixed_labels_table,
        [:protein_id, :final_cluster_id, :cluster_fixed_dataset_label],
        "fixed_labels_table"
    )

    normalized_elm_clusters = _normalize_protein_elm_clusters_table(
        protein_elm_clusters;
        elm_cluster_column = elm_cluster_column
    )

    final_dataset = leftjoin(copy(fixed_labels_table), normalized_elm_clusters; on = :protein_id)
    transform!(
        final_dataset,
        :cluster_fixed_dataset_label => (
            values -> _map_fixed_to_final_dataset_label.(
            values; missing_label = cleaned_missing_label)
        ) => :cluster_dataset_label
    )
    _ensure_final_dataset_allowed_labels(
        final_dataset;
        label_column = :cluster_dataset_label
    )

    final_dataset
end

"""
    write_final_dataset_csv_gz(final_dataset_table, output_path)

Write one final dataset table as a gzip-compressed CSV file.
"""
function write_final_dataset_csv_gz(
        final_dataset_table::AbstractDataFrame,
        output_path::AbstractString
)::Nothing
    _require_columns(
        final_dataset_table,
        [:protein_id, :final_cluster_id, :cluster_dataset_label],
        "final_dataset_table"
    )

    mkpath(dirname(output_path))
    open(GzipCompressorStream, output_path, "w") do io
        CSV.write(io, final_dataset_table)
    end

    nothing
end

"""
    _split_elm_clusters(value; separator)

Split one separator-joined ELM motif-cluster value into normalized token strings.
Missing/empty values return an empty vector.
"""
function _split_elm_clusters(
        value;
        separator::AbstractString
)::Vector{String}
    token = _normalize_elm_cluster_token(value)
    token === nothing && return String[]

    split_tokens = split(token, separator; keepempty = false)
    normalized_tokens = String[]
    for split_token in split_tokens
        normalized_token = _normalize_elm_cluster_token(split_token)
        normalized_token === nothing && continue
        push!(normalized_tokens, normalized_token)
    end

    normalized_tokens
end

"""
    _count_unique_elm_clusters(values; separator)

Count unique ELM motif-cluster IDs from one collection of separator-joined
cluster strings.
"""
function _count_unique_elm_clusters(
        values;
        separator::AbstractString
)::Int
    unique_tokens = Set{String}()
    for value in values
        for token in _split_elm_clusters(value; separator = separator)
            push!(unique_tokens, token)
        end
    end

    length(unique_tokens)
end

"""
    _protein_has_elm_cluster(value)

Return `true` when one protein ELM cluster assignment is non-missing and non-empty.
"""
function _protein_has_elm_cluster(value)::Bool
    _normalize_elm_cluster_token(value) !== nothing
end

"""
    _final_dataset_label_sort_key(label)

Return a stable sort key for final dataset labels. Known labels are ordered as:
`training`, `validation`, `testing`; unknown labels are sorted after known
labels in lexicographic order.
"""
function _final_dataset_label_sort_key(label::AbstractString)::Tuple{Int, Int, String}
    index = findfirst(==(label), _FINAL_DATASET_DEFAULT_LABEL_ORDER)
    isnothing(index) && return (1, 0, label)
    (0, index, label)
end

"""
    _summarize_final_dataset_by_label(
        final_dataset_table;
        label_column = :cluster_dataset_label,
        elm_cluster_column = :protein_elm_clusters,
        elm_cluster_separator = "-"
    )

Summarize proteins and final clusters by dataset label, including ELM-cluster
coverage statistics.
"""
function _summarize_final_dataset_by_label(
        final_dataset_table::AbstractDataFrame;
        label_column::Symbol = :cluster_dataset_label,
        elm_cluster_column::Symbol = :protein_elm_clusters,
        elm_cluster_separator::AbstractString = "-"
)::DataFrame
    _require_columns(
        final_dataset_table,
        [:final_cluster_id, label_column, elm_cluster_column],
        "final_dataset_table"
    )

    protein_summary = combine(
        groupby(final_dataset_table, label_column),
        nrow => :n_proteins,
        elm_cluster_column => (
            values -> count(_protein_has_elm_cluster, values)
        ) => :n_proteins_with_elm_clusters,
        elm_cluster_column => (
            values -> _count_unique_elm_clusters(
            values; separator = elm_cluster_separator)
        ) => :n_unique_elm_clusters
    )
    rename!(protein_summary, label_column => :dataset_label)

    cluster_summary = combine(
        groupby(final_dataset_table, :final_cluster_id),
        label_column => (
            values -> _single_nonmissing_or_missing(
            values; value_label = string(label_column))
        ) => :dataset_label
    )
    cluster_summary = combine(
        groupby(cluster_summary, :dataset_label),
        nrow => :n_final_clusters
    )

    summary = outerjoin(protein_summary, cluster_summary; on = :dataset_label, matchmissing = :equal)
    summary.n_proteins = coalesce.(summary.n_proteins, 0)
    summary.n_final_clusters = coalesce.(summary.n_final_clusters, 0)
    summary.n_proteins_with_elm_clusters = coalesce.(summary.n_proteins_with_elm_clusters, 0)
    summary.n_unique_elm_clusters = coalesce.(summary.n_unique_elm_clusters, 0)
    transform!(
        summary,
        :dataset_label => ByRow(value -> ismissing(value) ? "missing" :
                                         string(value)) => :dataset_label
    )

    summary._sort_key = _final_dataset_label_sort_key.(summary.dataset_label)
    sort!(summary, :_sort_key)
    select!(
        summary,
        :dataset_label,
        :n_proteins,
        :n_final_clusters,
        :n_proteins_with_elm_clusters,
        :n_unique_elm_clusters
    )
end

"""
    _render_final_dataset_label_summary_lines(summary)

Render markdown bullet lines for per-label final dataset summary rows.
"""
function _render_final_dataset_label_summary_lines(summary::AbstractDataFrame)::String
    _require_columns(
        summary,
        [
            :dataset_label,
            :n_proteins,
            :n_final_clusters,
            :n_proteins_with_elm_clusters,
            :n_unique_elm_clusters
        ],
        "summary"
    )

    lines = String[]
    for row in eachrow(summary)
        push!(
            lines,
            "- `$(row.dataset_label)`: proteins `$(row.n_proteins)`, final clusters `$(row.n_final_clusters)`, proteins with ELM clusters `$(row.n_proteins_with_elm_clusters)`, unique ELM clusters `$(row.n_unique_elm_clusters)`"
        )
    end

    join(lines, "\n")
end

"""
    _format_generated_at_utc(generated_at_utc)

Format one UTC generation timestamp as `YYYY-mm-dd HH:MM:SS UTC`.
"""
function _format_generated_at_utc(generated_at_utc::DateTime)::String
    "$(Dates.format(generated_at_utc, dateformat"yyyy-mm-dd HH:MM:SS")) UTC"
end

"""
    write_final_dataset_readme(
        output_path;
        output_csv_path,
        output_csv_filesystem_path = output_csv_path,
        source_fixed_labels_path,
        source_elm_motif_clustering_path,
        final_dataset_table,
        generated_at_utc = now(UTC),
        missing_label = FINAL_DATASET_LABEL_VALIDATION,
        elm_cluster_column = :protein_elm_clusters,
        elm_cluster_separator = "-",
        script_relpath = "scripts/13_assing_missing_labels.jl"
    )

Write one README for the final dataset artifact generated after assigning
missing labels and joining ELM motif-cluster mappings.
"""
function write_final_dataset_readme(
        output_path::AbstractString;
        output_csv_path::AbstractString,
        output_csv_filesystem_path::AbstractString = output_csv_path,
        source_fixed_labels_path::AbstractString,
        source_elm_motif_clustering_path::AbstractString,
        final_dataset_table::AbstractDataFrame,
        generated_at_utc::DateTime = now(UTC),
        missing_label::AbstractString = FINAL_DATASET_LABEL_VALIDATION,
        elm_cluster_column::Symbol = :protein_elm_clusters,
        elm_cluster_separator::AbstractString = "-",
        script_relpath::AbstractString = _FINAL_DATASET_DEFAULT_SCRIPT_RELPATH
)::Nothing
    cleaned_missing_label = _validate_final_dataset_label(
        missing_label;
        argument_name = "missing_label"
    )
    _require_columns(
        final_dataset_table,
        [:protein_id, :final_cluster_id, :cluster_dataset_label, elm_cluster_column],
        "final_dataset_table"
    )
    _ensure_final_dataset_allowed_labels(
        final_dataset_table;
        label_column = :cluster_dataset_label
    )

    isfile(output_csv_filesystem_path) || throw(
        ArgumentError(
        "Final dataset CSV file does not exist for README metadata: $(output_csv_filesystem_path)",
    ),
    )

    generated_at_text = _format_generated_at_utc(generated_at_utc)
    output_csv_size_bytes = filesize(output_csv_filesystem_path)
    output_csv_size_mb = output_csv_size_bytes / 1024^2
    label_summary = _summarize_final_dataset_by_label(
        final_dataset_table;
        label_column = :cluster_dataset_label,
        elm_cluster_column = elm_cluster_column,
        elm_cluster_separator = elm_cluster_separator
    )
    proteins_with_elm_clusters = count(_protein_has_elm_cluster, final_dataset_table[!, elm_cluster_column])
    unique_elm_clusters = _count_unique_elm_clusters(
        final_dataset_table[!, elm_cluster_column];
        separator = elm_cluster_separator
    )

    readme_text = """
# Final Dataset

This folder stores the final per-protein dataset, adding ELM motif-cluster
assignments and replacing missing fixed split labels with `$(cleaned_missing_label)`.

## Artifact

- Compressed CSV: `$(output_csv_path)` (`$(output_csv_size_bytes)` bytes; `$(round(output_csv_size_mb; digits = 3))` MiB)
- Generated at (UTC): `$(generated_at_text)`
- Number of proteins: `$(nrow(final_dataset_table))`
- Number of final clusters: `$(length(unique(final_dataset_table.final_cluster_id)))`
- Proteins with at least one ELM cluster assignment: `$(proteins_with_elm_clusters)`
- Unique ELM clusters assigned across all proteins: `$(unique_elm_clusters)`

## Data Dictionary

- `protein_elm_clusters`: one or more ELM motif-cluster IDs for each protein, joined by `$(elm_cluster_separator)` when multiple clusters are present.
- `cluster_dataset_label`: final split label constrained to `training`,
  `validation`, and `testing`; values are derived from
  `cluster_fixed_dataset_label` by keeping `training`, mapping `caid` to
  `testing`, and assigning missing values to `$(cleaned_missing_label)`.

## Counts by Dataset Label

$( _render_final_dataset_label_summary_lines(label_summary))

## Source Inputs

- Fixed labels table: `$(source_fixed_labels_path)`
- Protein ELM motif clustering table: `$(source_elm_motif_clustering_path)`

## Reproducibility

Run from repository root:

```bash
julia $(script_relpath)
```

Complete example:

```bash
julia $(script_relpath) --fixed_labels_path data/derived/fixed/cluster_fixed_labels.csv.gz --elm_motif_clustering_path data/derived/elm/elm_motif_clustering.csv.gz --output_dir data/dataset --output_filename cluster_dataset_labels.csv.gz --missing_label $(cleaned_missing_label)
```
"""

    write(output_path, readme_text)
    nothing
end
