using CodecZlib
using CSV
using DataFrames
using Dates

# Split-protocol constant: AlphaFold 2.3 training-data cutoff date.
const AF2_TRAINING_DATA_CUTOFF_DATE = Date(2018, 4, 30)
const FIXED_DATASET_LABEL_TRAINING = "training"
const FIXED_DATASET_LABEL_CAID = "caid"

# Backward-compatible internal aliases used by dataset-stats helpers.
const _AF2_TRAINING_DATA_CUTOFF_DATE = AF2_TRAINING_DATA_CUTOFF_DATE
const _FIXED_DATASET_LABEL_TRAINING = FIXED_DATASET_LABEL_TRAINING
const _FIXED_DATASET_LABEL_CAID = FIXED_DATASET_LABEL_CAID

"""
    _require_columns(df, required_columns, table_name)

Validate that `df` contains all `required_columns`. Throws an `ArgumentError`
listing missing columns when validation fails.
"""
function _require_columns(
        df::DataFrame,
        required_columns,
        table_name::AbstractString
)::Nothing
    existing = Set(Symbol.(names(df)))
    missing_columns = [Symbol(column)
                       for column in required_columns if !(Symbol(column) in existing)]

    isempty(missing_columns) || throw(
        ArgumentError(
        "$(table_name) is missing required columns: $(join(string.(missing_columns), ", "))",
    ),
    )

    return nothing
end

"""
    _minimum_date_or_missing(values)

Return the oldest non-missing date in `values`, or `missing` when no date is
available.
"""
function _minimum_date_or_missing(values)::Union{Missing, Date}
    nonmissing = collect(skipmissing(values))
    isempty(nonmissing) && return missing
    return minimum(nonmissing)
end

"""
    _any_true(values)

Return `true` when at least one value in `values` is a non-missing boolean
true.
"""
function _any_true(values)::Bool
    return any(value -> !ismissing(value) && Bool(value), values)
end

"""
    _normalize_oldest_pdb_release_dates_table(oldest_pdb_release_dates)

Ensure that `oldest_pdb_release_dates` includes:
- `acc`
- `oldest_pdb_release_date`
- `has_pdb_structure`
- `has_missing_pdb_release_date`

If optional boolean columns are absent, conservative defaults are inferred.
"""
function _normalize_oldest_pdb_release_dates_table(
        oldest_pdb_release_dates::DataFrame
)::DataFrame
    _require_columns(
        oldest_pdb_release_dates,
        [:acc, :oldest_pdb_release_date],
        "oldest_pdb_release_dates"
    )

    normalized = copy(oldest_pdb_release_dates)
    has_pdb_structure_present = :has_pdb_structure in Symbol.(names(normalized))
    has_missing_release_present = :has_missing_pdb_release_date in Symbol.(names(normalized))

    if !has_pdb_structure_present
        if has_missing_release_present
            transform!(
                normalized,
                [
                    :oldest_pdb_release_date, :has_missing_pdb_release_date] => ByRow(
                    (release_date,
                    has_missing_release) -> !ismissing(release_date) ||
                                            (!ismissing(has_missing_release) &&
                                             Bool(has_missing_release))
                ) => :has_pdb_structure
            )
        else
            # Conservative inference for legacy 2-column tables (`acc`,
            # `oldest_pdb_release_date`): each row is treated as a known PDB
            # mapping even when the oldest release date is missing.
            normalized.has_pdb_structure = trues(nrow(normalized))
        end
    end

    if !has_missing_release_present
        transform!(
            normalized,
            [
                :oldest_pdb_release_date, :has_pdb_structure] => ByRow((
                release_date, has_pdb_structure) -> !ismissing(has_pdb_structure) &&
                                                    Bool(has_pdb_structure) &&
                                                    ismissing(release_date)
            ) => :has_missing_pdb_release_date
        )
    end

    transform!(
        normalized,
        :has_pdb_structure => ByRow(x -> !ismissing(x) && Bool(x)) => :has_pdb_structure,
        :has_missing_pdb_release_date => ByRow(
            x -> !ismissing(x) && Bool(x)
        ) => :has_missing_pdb_release_date
    )

    return normalized
end

"""
    _fixed_dataset_label_for_cluster(
        cluster_has_caid,
        cluster_has_af2_training_structure,
        cluster_has_pdb_structure,
        cluster_has_missing_pdb_release_date
    )

Return the fixed dataset label for one final cluster according to the split
rules:
- caid overrides all if cluster has CAID members,
- otherwise training if cluster has AF2-pre-cutoff structures,
- otherwise training (conservative) if cluster has PDB structures and missing
  release dates,
- otherwise `missing` (unassigned).
"""
function _fixed_dataset_label_for_cluster(
        cluster_has_caid::Bool,
        cluster_has_af2_training_structure::Bool,
        cluster_has_pdb_structure::Bool,
        cluster_has_missing_pdb_release_date::Bool
)::Union{Missing, String}
    if cluster_has_caid
        return _FIXED_DATASET_LABEL_CAID
    end

    if cluster_has_af2_training_structure ||
       (cluster_has_pdb_structure && cluster_has_missing_pdb_release_date)
        return _FIXED_DATASET_LABEL_TRAINING
    end

    return missing
end

"""
    _build_cluster_split_labels(stats)

Build one row per `final_cluster_id` containing:
- `cluster_in_af2_training_set`
- `cluster_in_caid2`
- `cluster_in_caid3`
- `cluster_fixed_dataset_label`
"""
function _build_cluster_split_labels(stats::DataFrame)::DataFrame
    _require_columns(
        stats,
        [
            :final_cluster_id,
            :protein_in_af2_training_set,
            :protein_has_pdb_structure,
            :protein_has_missing_pdb_release_date,
            :protein_in_caid2,
            :protein_in_caid3
        ],
        "stats"
    )

    cluster_flags = combine(
        groupby(stats, :final_cluster_id),
        :protein_in_af2_training_set => _any_true => :cluster_in_af2_training_set,
        :protein_has_pdb_structure => _any_true => :cluster_has_pdb_structure,
        :protein_has_missing_pdb_release_date => _any_true => :cluster_has_missing_pdb_release_date,
        :protein_in_caid2 => _any_true => :cluster_in_caid2,
        :protein_in_caid3 => _any_true => :cluster_in_caid3
    )

    transform!(
        cluster_flags,
        [
            :cluster_in_caid2, :cluster_in_caid3] => ByRow((
            x, y) -> x || y) => :cluster_has_caid
    )
    transform!(
        cluster_flags,
        [
            :cluster_has_caid,
            :cluster_in_af2_training_set,
            :cluster_has_pdb_structure,
            :cluster_has_missing_pdb_release_date
        ] => ByRow(_fixed_dataset_label_for_cluster) => :cluster_fixed_dataset_label
    )

    return select(
        cluster_flags,
        :final_cluster_id,
        :cluster_in_af2_training_set,
        :cluster_in_caid2,
        :cluster_in_caid3,
        :cluster_fixed_dataset_label
    )
end

"""
    _to_string_set(values)

Normalize an iterable collection of identifiers to a `Set{String}`.
"""
function _to_string_set(values)::Set{String}
    return Set(String.(collect(values)))
end

"""
    _normalize_optional_bool_column!(df, column_name)

Normalize an optional boolean column in-place, replacing `missing` with
`false`.
"""
function _normalize_optional_bool_column!(
        df::DataFrame,
        column_name::Symbol
)::Nothing
    transform!(
        df,
        column_name => ByRow(x -> ismissing(x) ? false : Bool(x)) => column_name
    )
    return nothing
end

"""
    build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid2_disprot_ids,
        caid3_disprot_ids;
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE
    )

Build one per-protein table containing all fixed-label inputs and outputs used
for deterministic split assignment.

Returned columns include:
- all columns from `cluster_representatives`,
- `protein_disprot_id`,
- `protein_has_idr_binding_annotation`,
- `protein_oldest_pdb_release_date`,
- `protein_has_pdb_structure`,
- `protein_has_missing_pdb_release_date`,
- `protein_in_caid2`,
- `protein_in_caid3`,
- `protein_in_af2_training_set`,
- `cluster_in_af2_training_set`,
- `cluster_in_caid2`,
- `cluster_in_caid3`,
- `cluster_fixed_dataset_label`.
"""
function build_fixed_labels_table(
        cluster_representatives::DataFrame,
        disprot_annotations::DataFrame,
        oldest_pdb_release_dates::DataFrame,
        caid2_disprot_ids,
        caid3_disprot_ids;
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE
)::DataFrame
    _require_columns(
        cluster_representatives,
        [:protein_id, :final_cluster_id],
        "cluster_representatives"
    )
    _require_columns(
        disprot_annotations,
        [:acc, :disprot_id, :has_idr_binding_annotation],
        "disprot_annotations"
    )

    caid2_ids = _to_string_set(caid2_disprot_ids)
    caid3_ids = _to_string_set(caid3_disprot_ids)
    normalized_oldest_pdb_release_dates = _normalize_oldest_pdb_release_dates_table(
        oldest_pdb_release_dates
    )

    stats = leftjoin(
        cluster_representatives,
        select(
            disprot_annotations,
            :acc,
            :disprot_id => :protein_disprot_id,
            :has_idr_binding_annotation => :protein_has_idr_binding_annotation
        ),
        on = :protein_id => :acc
    )
    stats = leftjoin(
        stats,
        select(
            normalized_oldest_pdb_release_dates,
            :acc,
            :oldest_pdb_release_date => :protein_oldest_pdb_release_date,
            :has_pdb_structure => :protein_has_pdb_structure,
            :has_missing_pdb_release_date => :protein_has_missing_pdb_release_date
        ),
        on = :protein_id => :acc
    )

    transform!(
        stats,
        :protein_disprot_id => ByRow(x -> !ismissing(x) &&
                                          in(String(x), caid2_ids)) => :protein_in_caid2,
        :protein_disprot_id => ByRow(x -> !ismissing(x) &&
                                          in(String(x), caid3_ids)) => :protein_in_caid3
    )

    _normalize_optional_bool_column!(stats, :protein_has_pdb_structure)
    _normalize_optional_bool_column!(stats, :protein_has_missing_pdb_release_date)
    _normalize_optional_bool_column!(stats, :protein_has_idr_binding_annotation)

    transform!(
        stats,
        :protein_oldest_pdb_release_date => ByRow(
            x -> !ismissing(x) &&
                 x <= af2_training_cutoff_date
        ) => :protein_in_af2_training_set
    )

    cluster_split_labels = _build_cluster_split_labels(stats)
    stats = leftjoin(stats, cluster_split_labels, on = :final_cluster_id)

    cluster_columns = Symbol.(names(cluster_representatives))
    output_columns = unique([
        cluster_columns...,
        :protein_disprot_id,
        :protein_has_idr_binding_annotation,
        :protein_oldest_pdb_release_date,
        :protein_has_pdb_structure,
        :protein_has_missing_pdb_release_date,
        :protein_in_caid2,
        :protein_in_caid3,
        :protein_in_af2_training_set,
        :cluster_in_af2_training_set,
        :cluster_in_caid2,
        :cluster_in_caid3,
        :cluster_fixed_dataset_label
    ])
    stats = select(stats, output_columns)

    sort!(stats, :protein_id)
    return stats
end

"""
    write_fixed_labels_csv_gz(fixed_labels_table, output_path)

Write `fixed_labels_table` to a gzip-compressed CSV file.
"""
function write_fixed_labels_csv_gz(
        fixed_labels_table::DataFrame,
        output_path::AbstractString
)::Nothing
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        gzip_io = GzipCompressorStream(io)
        try
            CSV.write(gzip_io, fixed_labels_table)
        finally
            close(gzip_io)
        end
    end
    return nothing
end

"""
    _normalize_label_for_report(label)

Normalize a split label for reporting, using `"missing"` for `missing` values.
"""
function _normalize_label_for_report(label)::String
    ismissing(label) && return "missing"
    return String(label)
end

"""
    _count_labels(table, label_column)

Count rows by string-normalized values in `label_column`.
"""
function _count_labels(table::DataFrame, label_column::Symbol)::Dict{String, Int}
    _require_columns(table, [label_column], "table")
    counts = Dict{String, Int}()
    for value in table[!, label_column]
        label = _normalize_label_for_report(value)
        counts[label] = get(counts, label, 0) + 1
    end
    return counts
end

"""
    _render_readme_summary_lines(label_counts; prefix)

Render markdown bullet lines for one summary dictionary.
"""
function _render_readme_summary_lines(
        label_counts::Dict{String, Int};
        prefix::AbstractString
)::String
    ordered_labels = ["training", "caid", "missing"]
    lines = String[]

    for label in ordered_labels
        push!(lines, "- $(prefix) `$(label)`: `$(get(label_counts, label, 0))`")
    end

    for label in sort(setdiff(collect(keys(label_counts)), ordered_labels))
        push!(lines, "- $(prefix) `$(label)`: `$(label_counts[label])`")
    end

    return join(lines, "\n")
end

"""
    write_fixed_labels_readme(
        output_path;
        output_csv_path,
        output_csv_filesystem_path = output_csv_path,
        source_cluster_representatives_path,
        source_disprot_regions_path,
        source_pdb_release_map_path,
        source_caid2_binding_fasta_path,
        source_caid3_binding_fasta_path,
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE,
        fixed_labels_table
    )

Write a README for the fixed-label derived dataset.
"""
function write_fixed_labels_readme(
        output_path::AbstractString;
        output_csv_path::AbstractString,
        output_csv_filesystem_path::AbstractString = output_csv_path,
        source_cluster_representatives_path::AbstractString,
        source_disprot_regions_path::AbstractString,
        source_pdb_release_map_path::AbstractString,
        source_caid2_binding_fasta_path::AbstractString,
        source_caid3_binding_fasta_path::AbstractString,
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE,
        fixed_labels_table::DataFrame
)::Nothing
    _require_columns(
        fixed_labels_table,
        [:final_cluster_id, :cluster_fixed_dataset_label],
        "fixed_labels_table"
    )

    unique_clusters = unique(
        select(fixed_labels_table, :final_cluster_id, :cluster_fixed_dataset_label),
        :final_cluster_id
    )
    protein_label_counts = _count_labels(fixed_labels_table, :cluster_fixed_dataset_label)
    cluster_label_counts = _count_labels(unique_clusters, :cluster_fixed_dataset_label)
    isfile(output_csv_filesystem_path) || throw(
        ArgumentError(
        "Fixed-label CSV file does not exist for README metadata: $(output_csv_filesystem_path)",
    ),
    )
    output_csv_size_bytes = filesize(output_csv_filesystem_path)

    readme = """
# Fixed Dataset Labels

This folder contains the per-protein fixed cluster labels used as the
deterministic base split:
- `training`: cluster has at least one PDB release date on or before the AF2
  training cutoff (`$(af2_training_cutoff_date)`), or has PDB mappings with
  missing release dates (conservative assignment).
- `caid`: cluster contains at least one protein in CAID2 or CAID3.
- `missing`: cluster does not satisfy either rule above.

## Artifact

- Compressed CSV: `$(output_csv_path)` (`$(output_csv_size_bytes)` bytes)
- Number of proteins: `$(nrow(fixed_labels_table))`
- Number of final clusters: `$(nrow(unique_clusters))`

## Counts by Fixed Label

### Proteins
$( _render_readme_summary_lines(protein_label_counts; prefix = "label"))

### Final clusters
$( _render_readme_summary_lines(cluster_label_counts; prefix = "label"))

## Source Inputs

- Cluster representatives: `$(source_cluster_representatives_path)`
- DisProt regions JSON: `$(source_disprot_regions_path)`
- PDB release-date map: `$(source_pdb_release_map_path)`
- CAID2 binding FASTA: `$(source_caid2_binding_fasta_path)`
- CAID3 binding FASTA: `$(source_caid3_binding_fasta_path)`
"""

    write(output_path, readme)
    return nothing
end
