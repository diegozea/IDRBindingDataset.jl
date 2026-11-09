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

# ------------------------------------------------------------------
# Core validation and normalization helpers
# ------------------------------------------------------------------

"""
    _require_columns(df, required_columns, table_name)

Validate that `df` contains all `required_columns`. Throws an `ArgumentError`
listing missing columns when validation fails.
"""
function _require_columns(
        df::AbstractDataFrame,
        required_columns::AbstractVector{Symbol},
        table_name::AbstractString
)
    missing_columns = setdiff(required_columns, propertynames(df))

    isempty(missing_columns) || throw(
        ArgumentError(
        "$(table_name) is missing required columns: $(join(string.(missing_columns), ", "))",
    ),
    )

    nothing
end

"""
    _minimum_date_or_missing(values)

Return the oldest non-missing date in `values`, or `missing` when no date is
available.
"""
function _minimum_date_or_missing(values)
    (v = collect(skipmissing(values)); isempty(v) ? missing : minimum(v))
end

"""
    _parse_date_argument(date_text, argument_name)

Parse one date argument in `yyyy-mm-dd` format and throw an `ArgumentError` for
invalid input.
"""
function _parse_date_argument(
        date_text::AbstractString,
        argument_name::AbstractString
)::Date
    parsed = tryparse(Date, String(date_text), dateformat"yyyy-mm-dd")
    parsed === nothing && throw(
        ArgumentError(
        "Invalid $(argument_name): $(repr(date_text)). Expected format yyyy-mm-dd."
    ),
    )
    parsed
end

"""
    _any_true(values)

Return `true` when at least one value in `values` is a non-missing boolean
true.
"""
_any_true(values) = any(skipmissing(values))

"""
    _normalize_oldest_pdb_release_dates_table(oldest_pdb_release_dates)

Ensure that `oldest_pdb_release_dates` includes:
- `acc`
- `oldest_pdb_release_date`
- `has_pdb_structure`
- `has_missing_pdb_release_date`
"""
function _normalize_oldest_pdb_release_dates_table(
        oldest_pdb_release_dates::AbstractDataFrame
)
    _require_columns(
        oldest_pdb_release_dates,
        [:acc, :oldest_pdb_release_date, :has_pdb_structure, :has_missing_pdb_release_date],
        "oldest_pdb_release_dates"
    )

    normalized = copy(oldest_pdb_release_dates)
    _normalize_optional_bool_column!(normalized, :has_pdb_structure)
    _normalize_optional_bool_column!(normalized, :has_missing_pdb_release_date)

    normalized
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
    cluster_has_caid && return _FIXED_DATASET_LABEL_CAID
    (cluster_has_af2_training_structure ||
     (cluster_has_pdb_structure && cluster_has_missing_pdb_release_date)) &&
        return _FIXED_DATASET_LABEL_TRAINING
    missing
end

"""
    _build_cluster_split_labels(stats)

Build one row per `final_cluster_id` containing:
- `cluster_in_af2_training_set`
- one `cluster_in_caid<N>` column per CAID dataset present
- `cluster_fixed_dataset_label`
"""
function _build_cluster_split_labels(
        stats::AbstractDataFrame;
        protein_caid_columns::AbstractVector{Symbol} = Symbol[]
)
    _require_columns(
        stats,
        [
            :final_cluster_id,
            :protein_in_af2_training_set,
            :protein_has_pdb_structure,
            :protein_has_missing_pdb_release_date
        ],
        "stats"
    )
    _require_columns(stats, collect(protein_caid_columns), "stats")

    # Keep a predictable protein->cluster naming convention for CAID flags.
    cluster_caid_columns = Symbol[Symbol(replace(String(protein_column), "protein_" => "cluster_"))
                                  for protein_column in protein_caid_columns]
    cluster_caid_specs = [protein_column => _any_true => cluster_column
                          for (protein_column,
        cluster_column) in zip(protein_caid_columns, cluster_caid_columns)]

    cluster_flags = combine(
        groupby(stats, :final_cluster_id),
        :protein_in_af2_training_set => _any_true => :cluster_in_af2_training_set,
        :protein_has_pdb_structure => _any_true => :cluster_has_pdb_structure,
        :protein_has_missing_pdb_release_date =>
            _any_true => :cluster_has_missing_pdb_release_date,
        cluster_caid_specs...
    )

    # A cluster is CAID-positive if any CAID dataset contributes at least one protein.
    transform!(
        cluster_flags,
        AsTable(cluster_caid_columns) => ByRow(values -> any(values)) => :cluster_has_caid
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

    select(
        cluster_flags,
        :final_cluster_id,
        :cluster_in_af2_training_set,
        cluster_caid_columns...,
        :cluster_fixed_dataset_label
    )
end

# ------------------------------------------------------------------
# CAID dataset/version handling
# ------------------------------------------------------------------

"""
    _to_string_set(values)

Normalize an iterable collection of identifiers to a `Set{String}`.
"""
_to_string_set(values)::Set{String} = Set(String.(values))

"""
    _caid_number_from_membership_column(column_name, prefix)

Extract the CAID dataset number from one membership column name, for example:
- `protein_in_caid4` -> `4`
- `cluster_in_caid4` -> `4`
"""
function _caid_number_from_membership_column(
        column_name::Symbol,
        prefix::AbstractString
)::Int
    text = String(column_name)
    captures = match(Regex("^$(prefix)(\\d+)\$"), text)
    captures === nothing && throw(
        ArgumentError(
        "Invalid CAID membership column name: $(repr(text)). Expected prefix $(repr(prefix)) followed by a number."
    ),
    )
    capture = only(captures.captures)
    capture === nothing && throw(
        ArgumentError("Missing CAID membership column number capture for column: $(repr(text))."),
    )
    parse(Int, String(capture))
end

"""
    _extract_caid_dataset_number(dataset_key)

Extract a CAID dataset number from one dataset key (`Int` or `caid_<N>` string).
"""
function _extract_caid_dataset_number(dataset_key)::Int
    dataset_key isa Integer && return Int(dataset_key)

    text = String(dataset_key)
    captures = match(r"^caid_(\d+)$", text)
    captures === nothing && throw(
        ArgumentError(
        "Invalid CAID dataset key: $(repr(text)). Expected an integer key or a string formatted as caid_<number>."
    ),
    )
    capture = only(captures.captures)
    capture === nothing && throw(
        ArgumentError("Missing CAID dataset number capture for key: $(repr(text))."),
    )
    parse(Int, String(capture))
end

"""
    _discover_caid_binding_datasets(raw_data_dir)

Discover CAID datasets by listing folders under `raw_data_dir` whose names
start with `caid_`. Each discovered folder must contain `binding.fasta.gz`.

Returns a vector of named tuples:
- `name`
- `number`
- `binding_fasta_path`
"""
function _discover_caid_binding_datasets(raw_data_dir::AbstractString)
    isdir(raw_data_dir) ||
        throw(ArgumentError("Raw data directory does not exist: $(raw_data_dir)"))

    dataset_names = sort(filter(
        name -> startswith(name, "caid_") && isdir(joinpath(raw_data_dir, name)),
        readdir(raw_data_dir)
    ))
    isempty(dataset_names) && throw(
        ArgumentError(
        "No CAID dataset folders found under $(raw_data_dir). Expected folders named like caid_2, caid_3."
    ),
    )

    datasets = NamedTuple{
        (:name, :number, :binding_fasta_path), Tuple{String, Int, String}}[]
    for name in dataset_names
        dataset_number = _extract_caid_dataset_number(name)
        binding_fasta_path = normpath(joinpath(raw_data_dir, name, "binding.fasta.gz"))
        isfile(binding_fasta_path) || throw(
            ArgumentError("Missing binding FASTA for $(name): $(binding_fasta_path)"),
        )
        push!(datasets,
            (
                name = String(name),
                number = dataset_number,
                binding_fasta_path = binding_fasta_path
            ))
    end

    sort!(datasets; by = dataset -> dataset.number)
    names_by_number = Dict{Int, Vector{String}}()
    for dataset in datasets
        push!(get!(names_by_number, dataset.number, String[]), dataset.name)
    end
    duplicates = filter(entry -> length(last(entry)) > 1, collect(pairs(names_by_number)))
    isempty(duplicates) || throw(
        ArgumentError(
        "Duplicate CAID dataset version numbers found under $(raw_data_dir): " *
        join(
        ["caid_$(version) => $(join(sort(names), ", "))"
         for (version, names) in sort(duplicates; by = first)],
        "; "
    )
    ),
    )

    datasets
end

"""
    _require_caid_dataset_path(caid_datasets, required_number)

Return the binding FASTA path for one required CAID dataset number.
"""
function _require_caid_dataset_path(
        caid_datasets::AbstractVector,
        required_number::Integer
)::String
    index = findfirst(dataset -> dataset.number == required_number, caid_datasets)
    isnothing(index) || return caid_datasets[index].binding_fasta_path

    discovered = isempty(caid_datasets) ? "(none)" :
                 join((dataset.name for dataset in caid_datasets),
        ", ")
    throw(
        ArgumentError(
        "Required CAID dataset folder not found: caid_$(required_number). Discovered CAID folders: $(discovered)"
    ),
    )
end

"""
    _caid_ids_by_number(caid_disprot_ids_by_dataset)

Normalize CAID dataset ID collections to `Dict{Int, Set{String}}` keyed by
dataset number.
"""
function _caid_ids_by_number(caid_disprot_ids_by_dataset::AbstractDict)::Dict{
        Int, Set{String}}
    ids_by_number = Dict{Int, Set{String}}()
    source_key_by_number = Dict{Int, String}()
    for (dataset_key, dataset_ids) in pairs(caid_disprot_ids_by_dataset)
        dataset_number = _extract_caid_dataset_number(dataset_key)
        haskey(ids_by_number, dataset_number) && throw(
            ArgumentError(
            "Duplicate CAID dataset version $(dataset_number) detected for keys $(repr(source_key_by_number[dataset_number])) and $(repr(dataset_key))."
        ),
        )
        ids_by_number[dataset_number] = _to_string_set(dataset_ids)
        source_key_by_number[dataset_number] = string(dataset_key)
    end

    ids_by_number
end

"""
    _add_caid_membership_columns!(stats, caid_ids_by_number)

Add per-protein CAID membership columns to `stats` and return the list of
protein CAID columns that should be considered for cluster CAID assignment.
"""
function _add_caid_membership_columns!(
        stats::AbstractDataFrame,
        caid_ids_by_number::Dict{Int, Set{String}}
)::Vector{Symbol}
    protein_caid_columns = Symbol[]
    for dataset_number in sort(collect(keys(caid_ids_by_number)))
        column_name = Symbol("protein_in_caid$(dataset_number)")
        dataset_ids = caid_ids_by_number[dataset_number]
        transform!(
            stats,
            :protein_disprot_id =>
                ByRow(
                    value -> !ismissing(value) &&
                             in(String(value), dataset_ids)
                ) => column_name
        )
        push!(protein_caid_columns, column_name)
    end

    sort!(
        unique(protein_caid_columns);
        by = column_name -> _caid_number_from_membership_column(
            column_name,
            "protein_in_caid"
        )
    )
    protein_caid_columns
end

# ------------------------------------------------------------------
# Fixed-label table assembly
# ------------------------------------------------------------------

"""
    _normalize_optional_bool_column!(df, column_name)

Normalize an optional boolean column in-place, replacing `missing` with
`false`.
"""
function _normalize_optional_bool_column!(
        df::AbstractDataFrame,
        column_name::Symbol
)
    column = df[!, column_name]
    if eltype(column) <: Union{Missing, Bool}
        df[!, column_name] = coalesce.(column, false)
    else
        transform!(
            df,
            column_name => ByRow(x -> ismissing(x) ? false : Bool(x)) => column_name
        )
    end
    nothing
end

"""
    build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE
    )

Build one per-protein table containing all fixed-label inputs and outputs used
for deterministic split assignment.

Returned columns include:
- all columns from `cluster_representatives` except the retired
  `global_representative_id` and `local_representative_id` columns,
- `protein_disprot_id`,
- `protein_has_idr_binding_annotation`,
- `protein_oldest_pdb_release_date`,
- `protein_has_pdb_structure`,
- `protein_has_missing_pdb_release_date`,
- one `protein_in_caid<N>` column per provided CAID dataset,
- `protein_in_af2_training_set`,
- `cluster_in_af2_training_set`,
- one `cluster_in_caid<N>` column per provided CAID dataset,
- `cluster_fixed_dataset_label`.
"""
function _build_fixed_labels_table_from_caid_datasets(
        cluster_representatives::AbstractDataFrame,
        disprot_annotations::AbstractDataFrame,
        oldest_pdb_release_dates::AbstractDataFrame,
        caid_ids_by_number::Dict{Int, Set{String}};
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE
)
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

    # Add one protein-level CAID flag column per discovered CAID dataset version.
    protein_caid_columns = _add_caid_membership_columns!(stats, caid_ids_by_number)

    _normalize_optional_bool_column!(stats, :protein_has_pdb_structure)
    _normalize_optional_bool_column!(stats, :protein_has_missing_pdb_release_date)
    _normalize_optional_bool_column!(stats, :protein_has_idr_binding_annotation)

    transform!(
        stats,
        :protein_oldest_pdb_release_date =>
            ByRow(
                x -> coalesce(x <= af2_training_cutoff_date, false)
            ) => :protein_in_af2_training_set
    )

    cluster_split_labels = _build_cluster_split_labels(
        stats;
        protein_caid_columns = protein_caid_columns
    )
    stats = leftjoin(stats, cluster_split_labels, on = :final_cluster_id)

    # Keep cluster CAID columns aligned with the dynamic protein CAID column set.
    cluster_caid_columns = Symbol[Symbol(replace(String(protein_column), "protein_" => "cluster_"))
                                  for protein_column in protein_caid_columns]

    cluster_columns = _dataset_columns_without_representatives(cluster_representatives)
    output_columns = unique([
        cluster_columns...,
        :protein_disprot_id,
        :protein_has_idr_binding_annotation,
        :protein_oldest_pdb_release_date,
        :protein_has_pdb_structure,
        :protein_has_missing_pdb_release_date,
        protein_caid_columns...,
        :protein_in_af2_training_set,
        :cluster_in_af2_training_set,
        cluster_caid_columns...,
        :cluster_fixed_dataset_label
    ])
    stats = select(stats, output_columns)

    sort!(stats, :protein_id)
    stats
end

"""
    build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE
    )

Build fixed labels from many CAID dataset ID collections, keyed by dataset name
(for example `caid_2`, `caid_3`, `caid_4`) or by integer dataset number.
"""
function build_fixed_labels_table(
        cluster_representatives::AbstractDataFrame,
        disprot_annotations::AbstractDataFrame,
        oldest_pdb_release_dates::AbstractDataFrame,
        caid_disprot_ids_by_dataset::AbstractDict;
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE
)
    _build_fixed_labels_table_from_caid_datasets(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        _caid_ids_by_number(caid_disprot_ids_by_dataset);
        af2_training_cutoff_date = af2_training_cutoff_date
    )
end

# ------------------------------------------------------------------
# Artifact/report generation helpers
# ------------------------------------------------------------------

"""
    write_fixed_labels_csv_gz(fixed_labels_table, output_path)

Write `fixed_labels_table` to a gzip-compressed CSV file.
"""
function write_fixed_labels_csv_gz(
        fixed_labels_table::AbstractDataFrame,
        output_path::AbstractString
)
    _write_csv_gz(fixed_labels_table, output_path)
end

"""
    _normalize_caid_source_paths_by_version(; source_caid_binding_fasta_paths_by_version)

Normalize CAID source FASTA paths to `Dict{Int, String}` keyed by CAID version.
"""
function _normalize_caid_source_paths_by_version(;
        source_caid_binding_fasta_paths_by_version::AbstractDict
)::Dict{Int, String}
    paths_by_version = Dict{Int, String}()
    for (version_key, path) in pairs(source_caid_binding_fasta_paths_by_version)
        version = _extract_caid_dataset_number(version_key)
        haskey(paths_by_version, version) && throw(
            ArgumentError(
            "Duplicate CAID source path provided for version $(version) in source_caid_binding_fasta_paths_by_version."
        ),
        )
        paths_by_version[version] = String(path)
    end

    isempty(paths_by_version) && throw(
        ArgumentError("At least one CAID source FASTA path must be provided."),
    )

    paths_by_version
end

"""
    _render_caid_source_lines(paths_by_version)

Render markdown bullet lines for CAID source paths using keys named
`source_caidN_binding_fasta_path`.
"""
function _render_caid_source_lines(paths_by_version::Dict{Int, String})::String
    lines = String[]
    for version in sort(collect(keys(paths_by_version)))
        field_name = "source_caid$(version)_binding_fasta_path"
        push!(lines, "- `$(field_name)`: `$(paths_by_version[version])`")
    end
    join(lines, "\n")
end

"""
    write_fixed_labels_readme(
        output_path;
        output_csv_path,
        output_csv_filesystem_path = output_csv_path,
        source_cluster_representatives_path,
        source_disprot_regions_path,
        source_pdb_release_map_path,
        source_caid_binding_fasta_paths_by_version,
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
        source_caid_binding_fasta_paths_by_version::AbstractDict,
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE,
        fixed_labels_table::AbstractDataFrame
)
    _require_columns(
        fixed_labels_table,
        [:final_cluster_id, :cluster_fixed_dataset_label],
        "fixed_labels_table"
    )

    unique_clusters = unique(
        select(fixed_labels_table, :final_cluster_id, :cluster_fixed_dataset_label),
        :final_cluster_id
    )
    protein_label_counts = _count_labels_by_value(
        fixed_labels_table,
        :cluster_fixed_dataset_label;
        table_name = "fixed_labels_table"
    )
    cluster_label_counts = _count_labels_by_value(
        unique_clusters,
        :cluster_fixed_dataset_label;
        table_name = "unique_clusters"
    )
    isfile(output_csv_filesystem_path) || throw(
        ArgumentError(
        "Fixed-label CSV file does not exist for README metadata: $(output_csv_filesystem_path)",
    ),
    )
    output_csv_size_bytes = filesize(output_csv_filesystem_path)
    caid_source_paths_by_version = _normalize_caid_source_paths_by_version(
        source_caid_binding_fasta_paths_by_version = source_caid_binding_fasta_paths_by_version
    )

    readme = """
# Fixed Dataset Labels

This folder contains the per-protein fixed cluster labels used as the
deterministic base split:
- `training`: cluster has at least one PDB release date on or before the AF2
  training cutoff (`$(af2_training_cutoff_date)`), or has PDB mappings with
  missing release dates (conservative assignment).
- `caid`: cluster contains at least one protein in the CAID reference sets
  used for this build.
- `missing`: cluster does not satisfy either rule above.

## Artifact

- Compressed CSV: `$(output_csv_path)` (`$(output_csv_size_bytes)` bytes)
- Number of proteins: `$(nrow(fixed_labels_table))`
- Number of final clusters: `$(nrow(unique_clusters))`

## Counts by Fixed Label

### Proteins
$( _render_label_count_lines(
    protein_label_counts;
    ordered_labels = ["training", "caid", "missing"],
    prefix = "label"
))

### Final clusters
$( _render_label_count_lines(
    cluster_label_counts;
    ordered_labels = ["training", "caid", "missing"],
    prefix = "label"
))

## Source Inputs

- Cluster representatives: `$(source_cluster_representatives_path)`
- DisProt regions JSON: `$(source_disprot_regions_path)`
- PDB release-date map: `$(source_pdb_release_map_path)`
$( _render_caid_source_lines(caid_source_paths_by_version))
"""

    _write_text_file(output_path, readme)
end
