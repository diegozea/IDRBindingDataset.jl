using DataFrames
using Statistics

const DEFAULT_ELM_MOTIF_CLUSTER_NORMIC_THRESHOLD = 0.5
const DEFAULT_ELM_MOTIF_CLUSTER_MATCHIC_THRESHOLD = 1.5
const DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD = nothing
# Defaults follow the CompariMotif-based matching criteria described in:
# Palopoli, Lythgow, Edwards (2015), Bioinformatics, doi:10.1093/bioinformatics/btv155.
const DEFAULT_ELM_MOTIF_CLUSTER_MIN_SHARED_POSITIONS = 2

const _REQUIRED_REGEX_SIMILARITY_COLUMNS = [
    :query,
    :search,
    :matched,
    :normalized_ic,
    :match_ic
]
const _REQUIRED_REGEX_CLUSTER_COLUMNS = [:regex, :elm_motif_cluster]
const _REQUIRED_MOTIF_TABLE_COLUMNS_FOR_CLUSTERING = [:protein_id, :regex]
const _REQUIRED_PROTEIN_CLUSTER_COLUMNS = [:protein_id, :elm_motif_cluster]

"""
    _parse_optional_score_threshold(raw_value)

Parse one optional score-threshold CLI value.
Returns `nothing` when the value is empty or `"nothing"` (case-insensitive).
"""
function _parse_optional_score_threshold(
        raw_value::AbstractString
)::Union{Nothing, Float64}
    cleaned_value = strip(String(raw_value))
    isempty(cleaned_value) && return DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD
    lowercase(cleaned_value) == "nothing" &&
        return DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD

    parsed_value = tryparse(Float64, cleaned_value)
    parsed_value === nothing && throw(
        ArgumentError(
        "Invalid --score_threshold value: $(repr(raw_value)). Expected a floating-point number, empty string, or 'nothing'."
    ),
    )
    parsed_value
end

"""
    _coerce_similarity_bool(value, column_name)

Coerce one boolean-like similarity value to `Bool`.
"""
function _coerce_similarity_bool(value, column_name::AbstractString)::Bool
    value isa Bool && return value

    cleaned_text = lowercase(_strip_wrapping_quotes(String(coalesce(value, ""))))
    cleaned_text in ("true", "t", "1", "yes", "y") && return true
    cleaned_text in ("false", "f", "0", "no", "n", "") && return false

    throw(
        ArgumentError(
        "Invalid $(column_name) value in regex similarity table: $(repr(value)). Expected boolean-like value."
    ),
    )
end

"""
    _coerce_similarity_float(value, column_name)

Coerce one numeric similarity value to `Float64`.
"""
function _coerce_similarity_float(value, column_name::AbstractString)::Float64
    value isa Real && return Float64(value)

    cleaned_text = _strip_wrapping_quotes(String(coalesce(value, "")))
    parsed_value = tryparse(Float64, cleaned_text)
    parsed_value === nothing && throw(
        ArgumentError(
        "Invalid $(column_name) value in regex similarity table: $(repr(value)). Expected floating-point value."
    ),
    )

    parsed_value
end

"""
    _normalize_regexes(regexes)

Normalize motif regex strings by trimming whitespace/quotes, dropping empty
entries, preserving insertion order, and deduplicating values.
"""
function _normalize_regexes(regexes::AbstractVector)::Vector{String}
    normalized = String[]
    seen = Set{String}()

    for regex in regexes
        cleaned_regex = _strip_wrapping_quotes(String(regex))
        isempty(cleaned_regex) && continue
        cleaned_regex in seen && continue

        push!(seen, cleaned_regex)
        push!(normalized, cleaned_regex)
    end

    normalized
end

"""
    _build_regex_similarity_edges(
        regex_similarity_table;
        normic_threshold,
        matchic_threshold,
        score_threshold,
        allowed_regexes
    )

Build unique undirected regex edges that pass matching and threshold filters.
"""
function _build_regex_similarity_edges(
        regex_similarity_table::AbstractDataFrame;
        normic_threshold::Float64,
        matchic_threshold::Float64,
        score_threshold::Union{Nothing, Float64} = DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD,
        allowed_regexes::Set{String}
)::Vector{Tuple{String, String}}
    required_similarity_columns = score_threshold === nothing ?
                                  _REQUIRED_REGEX_SIMILARITY_COLUMNS :
                                  [_REQUIRED_REGEX_SIMILARITY_COLUMNS; :score]
    _require_columns(
        regex_similarity_table,
        required_similarity_columns,
        "regex_similarity_table"
    )

    passing_edges = Tuple{String, String}[]
    for row in eachrow(regex_similarity_table)
        query_regex = _strip_wrapping_quotes(String(row.query))
        search_regex = _strip_wrapping_quotes(String(row.search))

        isempty(query_regex) && continue
        isempty(search_regex) && continue
        query_regex in allowed_regexes || continue
        search_regex in allowed_regexes || continue
        query_regex == search_regex && continue

        _coerce_similarity_bool(row.matched, "matched") || continue
        normalized_ic = _coerce_similarity_float(row.normalized_ic, "normalized_ic")
        match_ic = _coerce_similarity_float(row.match_ic, "match_ic")

        (normalized_ic >= normic_threshold && match_ic >= matchic_threshold) || continue
        if score_threshold !== nothing
            score = _coerce_similarity_float(row.score, "score")
            score >= score_threshold || continue
        end

        edge = query_regex <= search_regex ? (query_regex, search_regex) :
               (search_regex, query_regex)
        push!(passing_edges, edge)
    end

    unique!(passing_edges)
    sort!(passing_edges)
    passing_edges
end

"""
    _build_regex_adjacency(regexes, edges)

Build one undirected adjacency map from regex nodes and edge list.
"""
function _build_regex_adjacency(
        regexes::Vector{String},
        edges::Vector{Tuple{String, String}}
)::Dict{String, Set{String}}
    adjacency = Dict{String, Set{String}}(regex => Set{String}() for regex in regexes)

    for (regex_a, regex_b) in edges
        push!(get!(adjacency, regex_a, Set{String}()), regex_b)
        push!(get!(adjacency, regex_b, Set{String}()), regex_a)
    end

    adjacency
end

"""
    _connected_component_ids(regexes, adjacency)

Assign deterministic connected-component IDs (`1, 2, ...`) to regex nodes.
"""
function _connected_component_ids(
        regexes::Vector{String},
        adjacency::Dict{String, Set{String}}
)::Dict{String, Int}
    component_ids = Dict{String, Int}()
    visited = Set{String}()
    next_component_id = 0

    for seed_regex in regexes
        seed_regex in visited && continue

        next_component_id += 1
        stack = String[seed_regex]
        while !isempty(stack)
            current_regex = pop!(stack)
            current_regex in visited && continue

            push!(visited, current_regex)
            component_ids[current_regex] = next_component_id

            neighbors = sort!(collect(get(adjacency, current_regex, Set{String}())))
            append!(stack, Iterators.reverse(neighbors))
        end
    end

    component_ids
end

"""
    build_regex_motif_clusters_table(
        regexes,
        regex_similarity_table;
        normic_threshold = DEFAULT_ELM_MOTIF_CLUSTER_NORMIC_THRESHOLD,
        matchic_threshold = DEFAULT_ELM_MOTIF_CLUSTER_MATCHIC_THRESHOLD,
        score_threshold = DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD
    )

Cluster regex motifs using thresholded similarity edges and connected
components. When `score_threshold` is provided, regex edges must also satisfy
`Score >= score_threshold`.

Returns one deterministic table with columns:
- `regex`
- `elm_motif_cluster`
"""
function build_regex_motif_clusters_table(
        regexes::AbstractVector,
        regex_similarity_table::AbstractDataFrame;
        normic_threshold::Float64 = DEFAULT_ELM_MOTIF_CLUSTER_NORMIC_THRESHOLD,
        matchic_threshold::Float64 = DEFAULT_ELM_MOTIF_CLUSTER_MATCHIC_THRESHOLD,
        score_threshold::Union{Nothing, Float64} = DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD
)::DataFrame
    normic_threshold >= 0 || throw(ArgumentError("normic_threshold must be >= 0."))
    matchic_threshold >= 0 || throw(ArgumentError("matchic_threshold must be >= 0."))
    score_threshold === nothing || score_threshold >= 0 ||
        throw(ArgumentError("score_threshold must be >= 0 when provided."))

    normalized_regexes = _normalize_regexes(regexes)
    unique!(normalized_regexes)
    sort!(normalized_regexes)

    isempty(normalized_regexes) && return DataFrame(
        regex = String[],
        elm_motif_cluster = Int[]
    )

    allowed_regexes = Set(normalized_regexes)
    similarity_edges = _build_regex_similarity_edges(
        regex_similarity_table;
        normic_threshold = normic_threshold,
        matchic_threshold = matchic_threshold,
        score_threshold = score_threshold,
        allowed_regexes = allowed_regexes
    )

    adjacency = _build_regex_adjacency(normalized_regexes, similarity_edges)
    component_ids = _connected_component_ids(normalized_regexes, adjacency)

    clusters = DataFrame(
        regex = normalized_regexes,
        elm_motif_cluster = [component_ids[regex] for regex in normalized_regexes]
    )

    sort!(clusters, [:elm_motif_cluster, :regex])
    clusters
end

"""
    build_protein_elm_motif_cluster_table(motif_table, regex_cluster_table)

Map regex clusters to proteins by joining motif rows on `regex`.

Returns one deterministic table with columns:
- `protein_id`
- `elm_motif_cluster`

A protein may appear in multiple rows when it matches multiple motif clusters.
"""
function build_protein_elm_motif_cluster_table(
        motif_table::AbstractDataFrame,
        regex_cluster_table::AbstractDataFrame
)::DataFrame
    _require_columns(motif_table, _REQUIRED_MOTIF_TABLE_COLUMNS_FOR_CLUSTERING, "motif_table")
    _require_columns(regex_cluster_table, _REQUIRED_REGEX_CLUSTER_COLUMNS, "regex_cluster_table")

    normalized_motifs = select(motif_table, :protein_id, :regex)
    transform!(normalized_motifs,
        :protein_id => ByRow(value -> _strip_wrapping_quotes(String(value))) => :protein_id)
    transform!(normalized_motifs, :regex =>
        ByRow(value -> _strip_wrapping_quotes(String(value))) => :regex)

    normalized_clusters = select(regex_cluster_table, :regex, :elm_motif_cluster)
    transform!(normalized_clusters,
        :regex => ByRow(value -> _strip_wrapping_quotes(String(value))) => :regex)

    protein_cluster_table = innerjoin(normalized_motifs, normalized_clusters; on = :regex)
    protein_cluster_table = select(protein_cluster_table, :protein_id, :elm_motif_cluster)

    unique!(protein_cluster_table)
    sort!(protein_cluster_table, [:protein_id, :elm_motif_cluster])
    protein_cluster_table
end

"""
    write_regex_motif_clusters_csv_gz(regex_cluster_table, output_path)

Write regex-to-cluster assignments to one gzip-compressed CSV file.
"""
function write_regex_motif_clusters_csv_gz(
        regex_cluster_table::AbstractDataFrame,
        output_path::AbstractString
)::Nothing
    _require_columns(regex_cluster_table, _REQUIRED_REGEX_CLUSTER_COLUMNS, "regex_cluster_table")

    _write_csv_gz(regex_cluster_table, output_path)
end

"""
    write_elm_motif_clustering_csv_gz(protein_cluster_table, output_path)

Write protein-to-motif-cluster assignments to one gzip-compressed CSV file.
"""
function write_elm_motif_clustering_csv_gz(
        protein_cluster_table::AbstractDataFrame,
        output_path::AbstractString
)::Nothing
    _require_columns(
        protein_cluster_table,
        _REQUIRED_PROTEIN_CLUSTER_COLUMNS,
        "protein_cluster_table"
    )

    _write_csv_gz(protein_cluster_table, output_path)
end

"""
    _summarize_regex_cluster_table(regex_cluster_table)

Compute summary statistics for one regex cluster table.
"""
function _summarize_regex_cluster_table(regex_cluster_table::AbstractDataFrame)::NamedTuple
    _require_columns(regex_cluster_table, _REQUIRED_REGEX_CLUSTER_COLUMNS, "regex_cluster_table")

    cluster_sizes = isempty(regex_cluster_table) ? Int[] :
                    combine(groupby(regex_cluster_table, :elm_motif_cluster), nrow =>
        :count).count

    (
        regex_count = nrow(regex_cluster_table),
        cluster_count = length(unique(regex_cluster_table.elm_motif_cluster)),
        singleton_cluster_count = count(==(1), cluster_sizes),
        max_cluster_size = isempty(cluster_sizes) ? 0 : maximum(cluster_sizes),
        mean_cluster_size = isempty(cluster_sizes) ? 0.0 : mean(cluster_sizes)
    )
end

"""
    _summarize_protein_cluster_table(protein_cluster_table)

Compute summary statistics for one protein cluster mapping table.
"""
function _summarize_protein_cluster_table(protein_cluster_table::AbstractDataFrame)::NamedTuple
    _require_columns(
        protein_cluster_table,
        _REQUIRED_PROTEIN_CLUSTER_COLUMNS,
        "protein_cluster_table"
    )

    clusters_per_protein = isempty(protein_cluster_table) ? Int[] :
                           combine(groupby(protein_cluster_table, :protein_id), nrow =>
        :count).count

    (
        row_count = nrow(protein_cluster_table),
        protein_count = length(unique(protein_cluster_table.protein_id)),
        cluster_count = length(unique(protein_cluster_table.elm_motif_cluster)),
        proteins_with_multiple_clusters = count(>(1), clusters_per_protein),
        mean_clusters_per_protein = isempty(clusters_per_protein) ? 0.0 :
                                    mean(clusters_per_protein),
        max_clusters_per_protein = isempty(clusters_per_protein) ? 0 :
                                   maximum(clusters_per_protein)
    )
end

"""
    write_elm_motif_clustering_readme(
        output_path;
        output_csv_path,
        output_csv_filesystem_path = output_csv_path,
        source_motifs_path,
        source_regex_clusters_path,
        min_shared_positions = DEFAULT_ELM_MOTIF_CLUSTER_MIN_SHARED_POSITIONS,
        normic_threshold,
        matchic_threshold,
        score_threshold = DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD,
        regex_cluster_table,
        protein_cluster_table,
        script_relpath = "scripts/14_elm_motif_clustering.jl"
    )

Write a README for the derived protein motif-cluster mapping artifact.
"""
function write_elm_motif_clustering_readme(
        output_path::AbstractString;
        output_csv_path::AbstractString,
        output_csv_filesystem_path::AbstractString = output_csv_path,
        source_motifs_path::AbstractString,
        source_regex_clusters_path::AbstractString,
        min_shared_positions::Int = DEFAULT_ELM_MOTIF_CLUSTER_MIN_SHARED_POSITIONS,
        normic_threshold::Float64,
        matchic_threshold::Float64,
        score_threshold::Union{Nothing, Float64} = DEFAULT_ELM_MOTIF_CLUSTER_SCORE_THRESHOLD,
        regex_cluster_table::AbstractDataFrame,
        protein_cluster_table::AbstractDataFrame,
        script_relpath::AbstractString = "scripts/14_elm_motif_clustering.jl"
)::Nothing
    _require_columns(regex_cluster_table, _REQUIRED_REGEX_CLUSTER_COLUMNS, "regex_cluster_table")
    _require_columns(
        protein_cluster_table,
        _REQUIRED_PROTEIN_CLUSTER_COLUMNS,
        "protein_cluster_table"
    )

    isfile(output_csv_filesystem_path) || throw(
        ArgumentError(
        "Protein motif-cluster CSV file does not exist for README metadata: $(output_csv_filesystem_path)",
    ),
    )

    regex_summary = _summarize_regex_cluster_table(regex_cluster_table)
    protein_summary = _summarize_protein_cluster_table(protein_cluster_table)
    output_csv_size_bytes = filesize(output_csv_filesystem_path)
    min_shared_positions > 0 ||
        throw(ArgumentError("min_shared_positions must be >= 1."))
    score_threshold_line = score_threshold === nothing ? "" :
                           "- Keep similarity edges with `Score >= $(score_threshold)`.\n"
    min_shared_positions_arg = min_shared_positions ==
                               DEFAULT_ELM_MOTIF_CLUSTER_MIN_SHARED_POSITIONS ? "" :
                               " --min_shared_positions $(min_shared_positions)"
    score_threshold_arg = score_threshold === nothing ? "" :
                          " --score_threshold $(score_threshold)"

    readme_text = """
# ELM Motif Clustering

This folder stores protein-level assignments to ELM motif clusters derived from
regex similarity graph connected components.

## Artifact

- Compressed CSV: `$(output_csv_path)` (`$(output_csv_size_bytes)` bytes)
- Rows (protein-cluster assignments): `$(protein_summary.row_count)`
- Unique proteins assigned to at least one motif cluster: `$(protein_summary.protein_count)`
- Unique motif clusters assigned to proteins: `$(protein_summary.cluster_count)`

## Data Dictionary

- `protein_id`: UniProt accession.
- `elm_motif_cluster`: connected-component cluster number derived from regex
  similarity graph.

## Clustering Protocol

- Build unique regex list from `motifs.csv.gz`.
- Compare regexes pairwise with `CompariMotif.jl` requiring at least
  `$(min_shared_positions)` matched positions.
- Keep similarity edges with `matched=true`, `NormIC >= $(normic_threshold)`, and
  `MatchIC >= $(matchic_threshold)`.
$(score_threshold_line)- Build an undirected graph and assign connected-component IDs.

## Regex Cluster Summary

- Unique regexes clustered: `$(regex_summary.regex_count)`
- Connected components: `$(regex_summary.cluster_count)`
- Singleton components: `$(regex_summary.singleton_cluster_count)`
- Maximum regexes in one cluster: `$(regex_summary.max_cluster_size)`
- Mean regexes per cluster: `$(round(regex_summary.mean_cluster_size; digits = 3))`

## Protein Mapping Summary

- Unique proteins mapped: `$(protein_summary.protein_count)`
- Proteins mapped to multiple motif clusters: `$(protein_summary.proteins_with_multiple_clusters)`
- Mean clusters per protein: `$(round(protein_summary.mean_clusters_per_protein; digits = 3))`
- Maximum clusters in one protein: `$(protein_summary.max_clusters_per_protein)`

## Source Inputs

- Protein-regex motifs: `$(source_motifs_path)`
- Regex clusters: `$(source_regex_clusters_path)`

## Reproducibility

Run from repository root:

```bash
julia $(script_relpath)
```

Complete example:

```bash
julia $(script_relpath) --motifs_path data/scratch/elm/mapping/motifs.csv.gz --regex_cluster_output_dir data/scratch/elm/clustering --regex_cluster_output_filename regex.csv.gz --output_dir data/derived/elm --output_filename elm_motif_clustering.csv.gz$(min_shared_positions_arg) --normic_threshold $(normic_threshold) --matchic_threshold $(matchic_threshold)$(score_threshold_arg)
```
"""

    _write_text_file(output_path, readme_text)
end
