using DataFrames
using Dates

const _FIXED_DATASET_LABEL_TESTING = "testing"
const _SPLIT_DATASET_LABELS = [
    FIXED_DATASET_LABEL_TRAINING, FIXED_DATASET_LABEL_CAID, _FIXED_DATASET_LABEL_TESTING]

"""
    _header_primary_token(header)

Extract the first whitespace-delimited token from a FASTA header string.
"""
function _header_primary_token(header::AbstractString)::String
    header = strip(String(header))
    isempty(header) &&
        throw(ArgumentError("Empty FASTA header identifier encountered in CAID FASTA file."))
    return split(header; limit = 2)[1]
end

"""
    read_caid_disprot_ids(caid_binding_fasta_gz_path)

Read one CAID FASTA `.gz` file and return the unique DisProt IDs found in FASTA
headers.
"""
function read_caid_disprot_ids(caid_binding_fasta_gz_path::AbstractString)::Set{String}
    isfile(caid_binding_fasta_gz_path) || throw(
        ArgumentError("CAID FASTA file does not exist: $(caid_binding_fasta_gz_path)"),
    )

    ids = Set{String}()
    sequences = MIToS.MSA.read_file(
        caid_binding_fasta_gz_path,
        MIToS.MSA.AnnotatedFASTASequences
    )
    for sequence in sequences
        raw_id = String(MIToS.MSA.sequence_id(sequence))
        push!(ids, _header_primary_token(raw_id))
    end

    return ids
end

"""
    _normalize_pdb_release_date(value)

Normalize one raw PDB release-date value to `Union{Missing, Date}`.

Accepted values:
- `Date`
- `missing`
- strings in `yyyy-mm-dd` format
- strings equal to `missing` (case-insensitive), which are converted to `missing`
"""
function _normalize_pdb_release_date(value)::Union{Missing, Date}
    ismissing(value) && return missing
    value isa Date && return value

    if value isa AbstractString
        text = strip(String(value))
        isempty(text) && return missing
        lowercase(text) == "missing" && return missing

        parsed = tryparse(Date, text)
        parsed === nothing && throw(
            ArgumentError("Invalid PDB release date string: $(repr(text))"),
        )
        return parsed
    end

    throw(ArgumentError("Unsupported PDB release-date value type: $(typeof(value))"))
end

"""
    _count_true(values)

Count values that are non-missing booleans equal to `true`.
"""
function _count_true(values)::Int
    return count(value -> !ismissing(value) && Bool(value), values)
end

"""
    _single_nonmissing_or_missing(values; value_label = "value")

Return the single non-missing value in `values`, or `missing` when all values
are missing. Throws an `ArgumentError` when multiple distinct non-missing
values are found.
"""
function _single_nonmissing_or_missing(
        values;
        value_label::AbstractString = "value"
)
    nonmissing_values = unique(collect(skipmissing(values)))
    if isempty(nonmissing_values)
        return missing
    end
    length(nonmissing_values) == 1 && return only(nonmissing_values)

    throw(
        ArgumentError(
        "Expected one non-missing $(value_label), found multiple values: $(join(string.(nonmissing_values), ", "))"
    ),
    )
end

"""
    _cluster_label_with_disprot_cutoff(
        fixed_label,
        cluster_min_oldest_disprot_release_date,
        disprot_testing_cutoff_date
    )

Resolve one cluster label:
- keep `fixed_label` when present,
- otherwise assign `testing` when cluster minimum oldest DisProt release date is
  after `disprot_testing_cutoff_date`,
- otherwise assign `training`.
"""
function _cluster_label_with_disprot_cutoff(
        fixed_label,
        cluster_min_oldest_disprot_release_date::Union{Missing, Date},
        disprot_testing_cutoff_date::Date
)::String
    !ismissing(fixed_label) && return String(fixed_label)
    if !ismissing(cluster_min_oldest_disprot_release_date) &&
       cluster_min_oldest_disprot_release_date > disprot_testing_cutoff_date
        return _FIXED_DATASET_LABEL_TESTING
    end
    return _FIXED_DATASET_LABEL_TRAINING
end

"""
    _split_label_sort_key(label)

Return a stable order key for split labels (`training`, `caid`,
`testing`), with unrecognized labels sorted after known labels.
"""
function _split_label_sort_key(label)::Int
    label_string = String(label)
    index = findfirst(==(label_string), _SPLIT_DATASET_LABELS)
    isnothing(index) && return length(_SPLIT_DATASET_LABELS) + 1
    return index
end

"""
    _sorted_caid_membership_columns(column_names, prefix)

Return CAID membership columns that match `prefix * <number>`, sorted by CAID
dataset version number.
"""
function _sorted_caid_membership_columns(
        column_names::AbstractVector{Symbol},
        prefix::AbstractString
)::Vector{Symbol}
    pattern = Regex("^$(prefix)(\\d+)\$")
    caid_columns = Symbol[]
    for column_name in column_names
        match(pattern, String(column_name)) === nothing || push!(caid_columns, column_name)
    end
    sort!(
        caid_columns;
        by = column_name -> _caid_number_from_membership_column(column_name, prefix)
    )
    caid_columns
end

"""
    _complete_split_label_summary_rows(summary)

Ensure `summary` contains one row for each split label
(`training`, `caid`, `testing`) with zero-filled count columns when
missing.
"""
function _complete_split_label_summary_rows(summary::DataFrame)::DataFrame
    _require_columns(
        summary,
        [
            :dataset_label,
            :n_proteins,
            :n_final_clusters,
            :n_idr_binding_proteins,
            :n_idr_binding_clusters
        ],
        "summary"
    )

    row_by_label = Dict(String(row.dataset_label) => row for row in eachrow(summary))
    completed = DataFrame(
        dataset_label = String[],
        n_proteins = Int[],
        n_final_clusters = Int[],
        n_idr_binding_proteins = Int[],
        n_idr_binding_clusters = Int[]
    )
    for label in _SPLIT_DATASET_LABELS
        if haskey(row_by_label, label)
            row = row_by_label[label]
            push!(
                completed,
                (
                    dataset_label = label,
                    n_proteins = Int(row.n_proteins),
                    n_final_clusters = Int(row.n_final_clusters),
                    n_idr_binding_proteins = Int(row.n_idr_binding_proteins),
                    n_idr_binding_clusters = Int(row.n_idr_binding_clusters)
                )
            )
        else
            push!(
                completed,
                (
                    dataset_label = label,
                    n_proteins = 0,
                    n_final_clusters = 0,
                    n_idr_binding_proteins = 0,
                    n_idr_binding_clusters = 0
                )
            )
        end
    end

    return completed
end

"""
    build_oldest_pdb_release_dates(pdb_release_map)

Collapse a PDB release-date table to one row per UniProt accession with the
oldest available PDB release date.

Expected input columns:
- `acc`
- `pdb_release_date`

Returns a `DataFrame` with columns:
- `acc`
- `oldest_pdb_release_date`
- `has_pdb_structure`
- `has_missing_pdb_release_date`
"""
function build_oldest_pdb_release_dates(pdb_release_map::DataFrame)::DataFrame
    _require_columns(pdb_release_map, [:acc, :pdb_release_date], "pdb_release_map")

    normalized = select(
        pdb_release_map,
        :acc => ByRow(x -> String(x)) => :acc,
        :pdb_release_date => ByRow(_normalize_pdb_release_date) => :_normalized_release_date
    )

    oldest_by_accession = combine(
        groupby(normalized, :acc),
        :_normalized_release_date => _minimum_date_or_missing => :oldest_pdb_release_date,
        :_normalized_release_date => (values -> any(ismissing, values)) => :has_missing_pdb_release_date
    )
    oldest_by_accession.has_pdb_structure = trues(nrow(oldest_by_accession))

    select!(
        oldest_by_accession,
        :acc,
        :oldest_pdb_release_date,
        :has_pdb_structure,
        :has_missing_pdb_release_date
    )
    sort!(oldest_by_accession, :acc)
    return oldest_by_accession
end

"""
    _normalize_disprot_release_tag(value)

Normalize one DisProt release tag (`YYYY_MM`) to `Union{Missing, Date}`.
"""
function _normalize_disprot_release_tag(value)::Union{Missing, Date}
    ismissing(value) && return missing
    value isa AbstractString || throw(
        ArgumentError("Unsupported DisProt release tag value type: $(typeof(value))"),
    )

    tag = strip(String(value))
    isempty(tag) && return missing
    return parse_disprot_release(tag)
end

"""
    _region_has_idr_binding_annotation(region, term_rules)

Return `true` when one region term ID belongs to the ontology-derived GO
binding challenge set.
"""
function _region_has_idr_binding_annotation(
        region,
        term_rules::_DisProtChallengeTermRules
)::Bool
    term_id = _region_term_id(region)
    term_id === nothing && return false
    startswith(term_id, "GO:") || return false
    term_id in term_rules.go_binding_terms
end

"""
    _disprot_release_dates_for_entry(entry)

Collect all release tags available for one DisProt protein entry (protein-level
and region-level), normalized as `Union{Missing, Date}`.
"""
function _disprot_release_dates_for_entry(entry)::Vector{Union{Missing, Date}}
    release_dates = Union{Missing, Date}[]
    push!(release_dates, _normalize_disprot_release_tag(get(entry, "released", missing)))

    for region in get(entry, "regions", Any[])
        push!(
            release_dates,
            _normalize_disprot_release_tag(get(region, "released", missing))
        )
    end

    return release_dates
end

"""
    build_disprot_release_dates(disprot_entries)

Extract all DisProt release dates from protein-level and region-level release
tags and return unique sorted dates.
"""
function build_disprot_release_dates(disprot_entries)::Vector{Date}
    release_dates = Date[]

    for entry in disprot_entries
        append!(release_dates, skipmissing(_disprot_release_dates_for_entry(entry)))
    end

    unique_release_dates = unique(release_dates)
    sort!(unique_release_dates)
    return unique_release_dates
end

"""
    build_disprot_protein_annotations(
        disprot_entries;
        term_rules,
        invalid_region_ids = nothing
    )

Build one row per DisProt protein with:
- `acc`
- `disprot_id`
- `oldest_disprot_release_date`
- `has_idr_binding_annotation`

Binding is determined from ontology-derived GO term membership in
`term_rules`.

When `invalid_region_ids` is provided, regions with `region_id` in the set are
ignored.
"""
function build_disprot_protein_annotations(
        disprot_entries;
        term_rules::_DisProtChallengeTermRules,
        invalid_region_ids::Union{Nothing, Set{String}} = nothing
)::DataFrame
    rows = DataFrame(
        acc = String[],
        disprot_id = String[],
        oldest_disprot_release_date = Union{Missing, Date}[],
        has_idr_binding_annotation = Bool[]
    )
    seen_accessions = Set{String}()

    for (entry_index, entry) in enumerate(disprot_entries)
        acc = strip(String(get(entry, "acc", "")))
        isempty(acc) && throw(
            ArgumentError("DisProt entry $(entry_index) is missing a non-empty `acc`."),
        )

        acc in seen_accessions && throw(
            ArgumentError("Duplicate DisProt accession detected: $(acc)"),
        )

        disprot_id = strip(String(get(entry, "disprot_id", "")))
        release_dates = _disprot_release_dates_for_entry(entry)
        oldest_disprot_release_date = _minimum_date_or_missing(release_dates)
        has_idr_binding_annotation = false
        for region in get(entry, "regions", Any[])
            invalid_region_ids !== nothing &&
                _region_is_marked_invalid(region, invalid_region_ids) && continue

            has_idr_binding_annotation = _region_has_idr_binding_annotation(region, term_rules)
            has_idr_binding_annotation && break
        end

        push!(seen_accessions, acc)
        push!(
            rows,
            (
                acc = acc,
                disprot_id = disprot_id,
                oldest_disprot_release_date = oldest_disprot_release_date,
                has_idr_binding_annotation = has_idr_binding_annotation
            )
        )
    end

    sort!(rows, :acc)
    return rows
end

"""
    build_dataset_stats_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE
    )

Build the per-protein stats table used for cluster-distribution analysis.

Output columns include all columns from `cluster_representatives`, plus:
- `protein_disprot_id`
- `protein_oldest_disprot_release_date`
- `protein_has_idr_binding_annotation`
- `protein_oldest_pdb_release_date`
- `protein_has_pdb_structure`
- `protein_has_missing_pdb_release_date`
- one `protein_in_caid<N>` column per provided CAID dataset
- `protein_in_af2_training_set`
- `cluster_in_af2_training_set`
- one `cluster_in_caid<N>` column per provided CAID dataset
- `cluster_fixed_dataset_label` (`training`, `caid`, or `missing`)
"""
function build_dataset_stats_table(
        cluster_representatives::DataFrame,
        disprot_annotations::DataFrame,
        oldest_pdb_release_dates::DataFrame,
        caid_disprot_ids_by_dataset::AbstractDict;
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE
)::DataFrame
    _require_columns(
        cluster_representatives,
        [:protein_id, :final_cluster_id],
        "cluster_representatives"
    )
    _require_columns(
        disprot_annotations,
        [:acc, :disprot_id, :oldest_disprot_release_date, :has_idr_binding_annotation],
        "disprot_annotations"
    )
    stats = build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = af2_training_cutoff_date
    )

    stats = leftjoin(
        stats,
        select(
            disprot_annotations,
            :acc,
            :oldest_disprot_release_date => :protein_oldest_disprot_release_date
        ),
        on = :protein_id => :acc
    )

    stats_column_names = Symbol.(names(stats))
    protein_caid_columns = _sorted_caid_membership_columns(
        stats_column_names,
        "protein_in_caid"
    )
    cluster_caid_columns = _sorted_caid_membership_columns(
        stats_column_names,
        "cluster_in_caid"
    )

    cluster_columns = Symbol.(names(cluster_representatives))
    output_columns = unique([
        cluster_columns...,
        :protein_disprot_id,
        :protein_oldest_disprot_release_date,
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
    return stats
end

"""
    assign_missing_cluster_labels_by_disprot_cutoff(
        dataset_stats;
        disprot_testing_cutoff_date
    )

Add a cutoff-based label for clusters with missing fixed labels.

Required `dataset_stats` columns:
- `final_cluster_id`
- `cluster_fixed_dataset_label`
- `protein_oldest_disprot_release_date`

Returned table includes all input columns plus:
- `cluster_min_protein_oldest_disprot_release_date`
- `cluster_disprot_cutoff_dataset_label`

Rules:
- Keep existing `cluster_fixed_dataset_label` values (`training`, `caid`)
  unchanged.
- For clusters with missing fixed labels:
  - assign `testing` when
    `cluster_min_protein_oldest_disprot_release_date > disprot_testing_cutoff_date`;
  - otherwise assign `training`.
"""
function assign_missing_cluster_labels_by_disprot_cutoff(
        dataset_stats::DataFrame;
        disprot_testing_cutoff_date::Date
)::DataFrame
    _require_columns(
        dataset_stats,
        [
            :final_cluster_id,
            :cluster_fixed_dataset_label,
            :protein_oldest_disprot_release_date
        ],
        "dataset_stats"
    )

    cluster_cutoff_labels = combine(
        groupby(dataset_stats, :final_cluster_id),
        :cluster_fixed_dataset_label => (
            values -> _single_nonmissing_or_missing(
            values;
            value_label = "cluster_fixed_dataset_label"
        )
        ) => :cluster_fixed_dataset_label,
        :protein_oldest_disprot_release_date => _minimum_date_or_missing => :cluster_min_protein_oldest_disprot_release_date
    )
    transform!(
        cluster_cutoff_labels,
        [
            :cluster_fixed_dataset_label,
            :cluster_min_protein_oldest_disprot_release_date
        ] => ByRow(
            (fixed_label,
            min_disprot_release_date) -> _cluster_label_with_disprot_cutoff(
            fixed_label,
            min_disprot_release_date,
            disprot_testing_cutoff_date
        )
        ) => :cluster_disprot_cutoff_dataset_label
    )
    select!(
        cluster_cutoff_labels,
        :final_cluster_id,
        :cluster_min_protein_oldest_disprot_release_date,
        :cluster_disprot_cutoff_dataset_label
    )

    return leftjoin(dataset_stats, cluster_cutoff_labels, on = :final_cluster_id)
end

"""
    summarize_dataset_labels_by_binding(
        dataset_stats;
        dataset_label_column
    )

Summarize protein and cluster counts by dataset label, including IDR-binding
annotation counts.

Required `dataset_stats` columns:
- `final_cluster_id`
- `protein_has_idr_binding_annotation`
- `dataset_label_column`

Returned columns:
- `dataset_label`
- `n_proteins`
- `n_final_clusters`
- `n_idr_binding_proteins`
- `n_idr_binding_clusters`

Missing labels are normalized to the string `"missing"` in `dataset_label`.
"""
function summarize_dataset_labels_by_binding(
        dataset_stats::DataFrame;
        dataset_label_column::Symbol
)::DataFrame
    _require_columns(
        dataset_stats,
        [:final_cluster_id, :protein_has_idr_binding_annotation, dataset_label_column],
        "dataset_stats"
    )

    protein_summary = combine(
        groupby(dataset_stats, dataset_label_column),
        nrow => :n_proteins,
        :protein_has_idr_binding_annotation => _count_true => :n_idr_binding_proteins
    )
    rename!(protein_summary, dataset_label_column => :dataset_label)

    cluster_binding = combine(
        groupby(dataset_stats, :final_cluster_id),
        dataset_label_column => (
            values -> _single_nonmissing_or_missing(
            values;
            value_label = string(dataset_label_column)
        )
        ) => :dataset_label,
        :protein_has_idr_binding_annotation => _any_true => :cluster_has_idr_binding_annotation
    )
    cluster_summary = combine(
        groupby(cluster_binding, :dataset_label),
        nrow => :n_final_clusters,
        :cluster_has_idr_binding_annotation => _count_true => :n_idr_binding_clusters
    )

    summary = outerjoin(
        protein_summary,
        cluster_summary,
        on = :dataset_label,
        matchmissing = :equal
    )
    summary.n_proteins = coalesce.(summary.n_proteins, 0)
    summary.n_final_clusters = coalesce.(summary.n_final_clusters, 0)
    summary.n_idr_binding_proteins = coalesce.(summary.n_idr_binding_proteins, 0)
    summary.n_idr_binding_clusters = coalesce.(summary.n_idr_binding_clusters, 0)
    transform!(
        summary,
        :dataset_label => ByRow(value -> _normalize_label_for_report(value)) => :dataset_label
    )
    summary._dataset_label_sort_key = _split_label_sort_key.(summary.dataset_label)
    sort!(summary, :_dataset_label_sort_key)
    return select!(
        summary,
        :dataset_label,
        :n_proteins,
        :n_final_clusters,
        :n_idr_binding_proteins,
        :n_idr_binding_clusters
    )
end

"""
    sweep_disprot_cutoff_label_summaries(
        dataset_stats;
        disprot_testing_cutoff_dates
    )

For each cutoff date, assign missing cluster fixed labels with
`assign_missing_cluster_labels_by_disprot_cutoff` and compute summary counts via
`summarize_dataset_labels_by_binding`.

The output includes one row per (`disprot_testing_cutoff_date`, `dataset_label`)
with labels completed to always include:
`training`, `caid`, `testing`.
"""
function sweep_disprot_cutoff_label_summaries(
        dataset_stats::DataFrame;
        disprot_testing_cutoff_dates
)::DataFrame
    cutoff_dates = collect(disprot_testing_cutoff_dates)
    isempty(cutoff_dates) && throw(
        ArgumentError("`disprot_testing_cutoff_dates` must contain at least one date."),
    )

    sweep_rows = DataFrame(
        disprot_testing_cutoff_date = Date[],
        dataset_label = String[],
        n_proteins = Int[],
        n_final_clusters = Int[],
        n_idr_binding_proteins = Int[],
        n_idr_binding_clusters = Int[]
    )

    for cutoff_date in cutoff_dates
        labeled = assign_missing_cluster_labels_by_disprot_cutoff(
            dataset_stats;
            disprot_testing_cutoff_date = cutoff_date
        )
        summary = summarize_dataset_labels_by_binding(
            labeled;
            dataset_label_column = :cluster_disprot_cutoff_dataset_label
        )
        completed_summary = _complete_split_label_summary_rows(summary)
        completed_summary.disprot_testing_cutoff_date = fill(cutoff_date, nrow(completed_summary))
        append!(sweep_rows, completed_summary; cols = :setequal)
    end

    sweep_rows._dataset_label_sort_key = _split_label_sort_key.(sweep_rows.dataset_label)
    sort!(sweep_rows, [:disprot_testing_cutoff_date, :_dataset_label_sort_key])
    return select!(
        sweep_rows,
        :disprot_testing_cutoff_date,
        :dataset_label,
        :n_proteins,
        :n_final_clusters,
        :n_idr_binding_proteins,
        :n_idr_binding_clusters
    )
end

"""
    build_dataset_stats_table(
        cluster_representatives,
        disprot_entries,
        pdb_release_map,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE,
        term_rules,
        invalid_region_ids = nothing
    )

Convenience overload that derives the DisProt annotation and oldest-PDB tables
from raw DisProt entries and raw PDB release-map rows.
"""
function build_dataset_stats_table(
        cluster_representatives::DataFrame,
        disprot_entries,
        pdb_release_map::DataFrame,
        caid_disprot_ids_by_dataset::AbstractDict;
        af2_training_cutoff_date::Date = AF2_TRAINING_DATA_CUTOFF_DATE,
        term_rules::_DisProtChallengeTermRules,
        invalid_region_ids::Union{Nothing, Set{String}} = nothing
)::DataFrame
    disprot_annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = term_rules,
        invalid_region_ids = invalid_region_ids
    )
    oldest_pdb_release_dates = build_oldest_pdb_release_dates(pdb_release_map)

    return build_dataset_stats_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset;
        af2_training_cutoff_date = af2_training_cutoff_date
    )
end
