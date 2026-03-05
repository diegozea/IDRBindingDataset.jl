using CSV
using DataFrames
using DataStructures
using Diamond_jll
using ArgCheck

# ---------------------------------------------------------------------------
# Clustering thresholds and output conventions.
# ---------------------------------------------------------------------------
const GLOBAL_CLUSTER_MIN_COVERAGE = 80.0
const GLOBAL_CLUSTER_MIN_IDENTITY = 30.0
const LOCAL_CLUSTER_MIN_IDENTITY = 80.0
const LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH = 10

const DISPROT_FULL_PROTEIN_FASTA_GZ_FILENAME = "disprot_full_proteins.fasta.gz"
const CLUSTERING_README_FILENAME = "README.md"
const CLUSTER_REPRESENTATIVES_FILENAME = "cluster_representatives.csv.gz"

const _DIAMOND_BLASTP_OUTFMT_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qlen",
    "slen"
]
const _DIAMOND_BLASTP_OUTFMT_COLUMNS = Symbol.(_DIAMOND_BLASTP_OUTFMT_FIELDS)
const _DIAMOND_BLASTP_OUTFMT_TYPES = Dict(
    :qseqid => String,
    :sseqid => String,
    :pident => Float64,
    :length => Int,
    :gapopen => Int,
    :qstart => Int,
    :qend => Int,
    :sstart => Int,
    :send => Int,
    :qlen => Int,
    :slen => Int
)

const _DIAMOND_CLUSTER_PAIR_COLUMNS = [:representative_id, :sequence_id]
const _DIAMOND_CLUSTER_PAIR_TYPES = Dict(
    :representative_id => String,
    :sequence_id => String
)
const _CLUSTER_ASSIGNMENT_COLUMN_NAMES = [
    "cluster_id", "representative_id", "sequence_id", "cluster_size"]
const _CLUSTER_REPRESENTATIVE_COLUMN_NAMES = [
    "protein_id",
    "global_representative_id",
    "local_representative_id",
    "final_cluster_id"
]

const _DisProtProteinSequence = NamedTuple{(:acc, :sequence), Tuple{String, String}}
const _ClusteringThresholds = NamedTuple{
    (:global_min_coverage,
        :global_min_identity,
        :local_min_identity,
        :local_min_ungapped_length),
    Tuple{Float64, Float64, Float64, Int}
}

"""
    _validate_clustering_thresholds(;
        global_min_coverage,
        global_min_identity,
        local_min_identity,
        local_min_ungapped_length)

Validate clustering thresholds and return a normalized named tuple.
"""
function _validate_clustering_thresholds(;
        global_min_coverage::Real,
        global_min_identity::Real,
        local_min_identity::Real,
        local_min_ungapped_length::Integer
)::_ClusteringThresholds
    percent_values = (
        global_min_coverage = Float64(global_min_coverage),
        global_min_identity = Float64(global_min_identity),
        local_min_identity = Float64(local_min_identity)
    )
    h_len = Int(local_min_ungapped_length)

    for (name, value) in pairs(percent_values)
        @argcheck(isfinite(value) && 0.0 <= value <= 100.0,
            "`$(name)` must be a finite percentage in [0, 100] (use 30 for 30%, not 0.30).")
    end
    @argcheck(h_len >= 1, "`local_min_ungapped_length` must be an integer >= 1.")

    return (
        global_min_coverage = percent_values.global_min_coverage,
        global_min_identity = percent_values.global_min_identity,
        local_min_identity = percent_values.local_min_identity,
        local_min_ungapped_length = h_len
    )
end

# ---------------------------------------------------------------------------
# DisProt sequence extraction.
# ---------------------------------------------------------------------------
"""
    _normalize_disprot_sequence(sequence)

Normalize a DisProt protein sequence by trimming whitespace and uppercasing.
Only alphabetic amino-acid letters (`A`-`Z`) are accepted.
"""
function _normalize_disprot_sequence(sequence::AbstractString)::String
    normalized = uppercase(strip(String(sequence)))
    isempty(normalized) &&
        throw(ArgumentError("DisProt entry contains an empty `sequence` value."))

    all(c -> c in 'A':'Z', normalized) || throw(
        ArgumentError("DisProt sequence contains non-letter characters: $(sequence)"),
    )

    return normalized
end

"""
    extract_disprot_full_protein_sequences(disprot_entries)

Extract one `(acc, sequence)` pair per DisProt protein entry.

Validation rules:
- `acc` must be present and non-empty,
- `sequence` must be present and non-empty,
- accessions must be unique,
- sequences are normalized by `_normalize_disprot_sequence`.
"""
function extract_disprot_full_protein_sequences(disprot_entries)::Vector{_DisProtProteinSequence}
    records = _DisProtProteinSequence[]
    seen_accessions = Set{String}()

    for (entry_index, entry) in enumerate(disprot_entries)
        acc_value = get(entry, "acc", missing)
        acc_value isa AbstractString || throw(
            ArgumentError("DisProt entry $(entry_index) is missing a string `acc` field."),
        )

        acc = strip(String(acc_value))
        isempty(acc) &&
            throw(ArgumentError("DisProt entry $(entry_index) contains an empty `acc` value."))
        acc in seen_accessions &&
            throw(ArgumentError("Duplicate DisProt accession detected: $(acc)"))

        sequence_value = get(entry, "sequence", missing)
        sequence_value isa AbstractString || throw(
            ArgumentError("DisProt entry $(entry_index) ($(acc)) is missing a string `sequence` field."),
        )

        sequence = _normalize_disprot_sequence(sequence_value)
        push!(seen_accessions, acc)
        push!(records, (acc = acc, sequence = sequence))
    end

    isempty(records) && throw(ArgumentError("No DisProt protein sequences were found."))
    return records
end

"""
    write_disprot_full_protein_fasta_gz(disprot_regions_path, output_fasta_gz_path; line_width = 80)

Extract full protein sequences from a DisProt JSON file and save them as
compressed FASTA (`.fasta.gz`).

Returns the number of records written.
"""
function write_disprot_full_protein_fasta_gz(
        disprot_regions_path::AbstractString,
        output_fasta_gz_path::AbstractString;
        line_width::Int = 80
)::Int
    isfile(disprot_regions_path) || throw(
        ArgumentError("DisProt regions file does not exist: $(disprot_regions_path)"),
    )

    # Keep this argument for API compatibility with prior revisions.
    # MIToS FASTA sequence writers currently emit one sequence line per record.
    line_width > 0 || throw(ArgumentError("`line_width` must be a positive integer."))

    entries = _read_disprot_entries(disprot_regions_path)
    records = extract_disprot_full_protein_sequences(entries)
    sequences = [MIToS.MSA.AnnotatedSequence(record.acc, record.sequence)
                 for record in records]

    mkpath(dirname(String(output_fasta_gz_path)))
    MIToS.MSA.write_file(String(output_fasta_gz_path), sequences, MIToS.MSA.FASTASequences)

    return length(records)
end

# ---------------------------------------------------------------------------
# DIAMOND command helpers.
# ---------------------------------------------------------------------------
"""
    _diamond_cmd(args)

Build a DIAMOND command by appending `args` to the JLL-provided executable.
"""
function _diamond_cmd(args::AbstractVector{<:AbstractString})::Cmd
    # Interpolate the entire argument vector so each token remains a separate
    # CLI argument (`$args`), instead of concatenating token text.
    return `$(Diamond_jll.diamond()) $args`
end

"""
    _run_checked_command(cmd, context)

Run `cmd` and throw a descriptive error if it exits with a non-zero status.
"""
function _run_checked_command(cmd::Cmd, context::AbstractString)::Nothing
    process = run(Cmd(cmd; ignorestatus = true))
    process.exitcode == 0 || throw(
        ErrorException(
        "DIAMOND command failed for $(context) with exit code $(process.exitcode): $(cmd)"
    ),
    )
    return nothing
end

"""
    _build_all_vs_all_blastp_args(input_fasta_path, db_prefix, blast_output_path; threads)

Build DIAMOND `blastp` arguments for all-vs-all local homology discovery.

Parameter rationale:
- `--evalue 1000`: keep search permissive; final inclusion is controlled by
  explicit Julia-side thresholds (identity, coverage, ungapped length).
- `--max-target-seqs 0`: keep all candidate neighbors for graph edges.
- `--max-hsps 0`: keep alternative local alignments (needed for local homology).
- `--no-self-hits`: suppress trivial self alignments.
- `--ultra-sensitive`: maximize sensitivity for local hit discovery.
"""
function _build_all_vs_all_blastp_args(
        input_fasta_path::AbstractString,
        db_prefix::AbstractString,
        blast_output_path::AbstractString;
        threads::Int
)::Vector{String}
    args = [
        "blastp",
        "--db",
        String(db_prefix),
        "--query",
        String(input_fasta_path),
        "--out",
        String(blast_output_path),
        "--outfmt",
        "6",
        _DIAMOND_BLASTP_OUTFMT_FIELDS...,
        "--evalue",
        "1000",
        "--max-target-seqs",
        "0",
        "--max-hsps",
        "0",
        "--ultra-sensitive",
        "--threads",
        string(threads),
        "--masking",
        "0",
        "--quiet",
        "--no-self-hits"
    ]
    return args
end

"""
    _build_global_cluster_args(input_fasta_path, cluster_output_path; threads, memory_limit)

Build DIAMOND `cluster` arguments for the global criterion.

Parameter rationale:
- `--approx-id`: configured global identity threshold for clustering.
- `--mutual-cover`: configured bidirectional coverage requirement.
- `--memory-limit`: explicit memory control for deterministic runtime behavior.
"""
function _build_global_cluster_args(
        input_fasta_path::AbstractString,
        cluster_output_path::AbstractString;
        threads::Int,
        memory_limit::AbstractString,
        global_min_identity::Real,
        global_min_coverage::Real
)::Vector{String}
    args = [
        "cluster",
        "-d",
        String(input_fasta_path),
        "--out",
        String(cluster_output_path),
        "--approx-id",
        string(Float64(global_min_identity)),
        "--mutual-cover",
        string(Float64(global_min_coverage)),
        "--memory-limit",
        String(memory_limit),
        "--threads",
        string(threads),
        "--masking",
        "0",
        "--quiet"
    ]
    return args
end

"""
    _build_global_recluster_args(input_fasta_path, input_clusters_path, output_clusters_path; threads, memory_limit)

Build DIAMOND `recluster` arguments for the global criterion.
"""
function _build_global_recluster_args(
        input_fasta_path::AbstractString,
        input_clusters_path::AbstractString,
        output_clusters_path::AbstractString;
        threads::Int,
        memory_limit::AbstractString,
        global_min_identity::Real,
        global_min_coverage::Real
)::Vector{String}
    args = [
        "recluster",
        "-d",
        String(input_fasta_path),
        "--clusters",
        String(input_clusters_path),
        "--out",
        String(output_clusters_path),
        "--approx-id",
        string(Float64(global_min_identity)),
        "--mutual-cover",
        string(Float64(global_min_coverage)),
        "--memory-limit",
        String(memory_limit),
        "--threads",
        string(threads),
        "--masking",
        "0",
        "--quiet"
    ]
    return args
end

# ---------------------------------------------------------------------------
# FASTA readers for downstream clustering.
# ---------------------------------------------------------------------------
"""
    _read_fasta_ids(input_fasta_path)

Read unique FASTA sequence identifiers from `input_fasta_path`.
FASTA parsing is delegated to `MIToS.MSA.read_file(..., FASTASequences)`.
"""
function _read_fasta_ids(input_fasta_path::AbstractString)::Vector{String}
    @argcheck(isfile(input_fasta_path),
        "Input FASTA file does not exist: $(input_fasta_path)")

    sequences = MIToS.MSA.read_file(String(input_fasta_path), MIToS.MSA.FASTASequences)
    ids = [strip(String(MIToS.MSA.sequence_id(sequence))) for sequence in sequences]

    isempty(ids) &&
        throw(ErrorException("No FASTA headers were found in $(input_fasta_path)."))
    any(isempty, ids) &&
        throw(ErrorException("Encountered at least one empty FASTA identifier in $(input_fasta_path)."))
    allunique(ids) || throw(ErrorException("Duplicate FASTA identifiers detected."))
    return ids
end

# ---------------------------------------------------------------------------
# Disjoint-set helpers.
# ---------------------------------------------------------------------------
"""
    _sequence_index_map(sequence_ids)

Build a sequence identifier to index mapping for disjoint-set operations.
"""
function _sequence_index_map(sequence_ids::Vector{String})::Dict{String, Int}
    return Dict(sequence_id => sequence_index
    for
    (sequence_index, sequence_id) in enumerate(sequence_ids))
end

"""
    _components_from_disjoint_set!(disjoint_set, sequence_ids)

Build deterministic connected components from `disjoint_set` membership.
"""
function _components_from_disjoint_set!(
        disjoint_set::IntDisjointSet,
        sequence_ids::Vector{String}
)::Vector{Vector{String}}
    components = Dict{Int, Vector{String}}()
    for (sequence_index, sequence_id) in enumerate(sequence_ids)
        root_index = find_root!(disjoint_set, sequence_index)
        push!(get!(components, root_index, String[]), sequence_id)
    end

    sorted_components = [sort!(members) for members in values(components)]
    sort!(sorted_components; by = first)
    return sorted_components
end

"""
    _build_cluster_assignments(mode_prefix, components)

Convert connected components to a deterministic cluster-assignment table with
columns: `cluster_id`, `representative_id`, `sequence_id`, `cluster_size`.
"""
function _build_cluster_assignments(
        mode_prefix::AbstractString,
        components::Vector{Vector{String}}
)::DataFrame
    row_count = sum(length, components)
    cluster_id_column = Vector{String}(undef, row_count)
    representative_id_column = Vector{String}(undef, row_count)
    sequence_id_column = Vector{String}(undef, row_count)
    cluster_size_column = Vector{Int}(undef, row_count)

    row_index = 1
    for (component_index, members) in enumerate(components)
        cluster_id = "$(mode_prefix)$(lpad(string(component_index), 6, '0'))"
        representative_id = first(members)
        cluster_size = length(members)

        for sequence_id in members
            cluster_id_column[row_index] = cluster_id
            representative_id_column[row_index] = representative_id
            sequence_id_column[row_index] = sequence_id
            cluster_size_column[row_index] = cluster_size
            row_index += 1
        end
    end

    return DataFrame(
        cluster_id = cluster_id_column,
        representative_id = representative_id_column,
        sequence_id = sequence_id_column,
        cluster_size = cluster_size_column
    )
end

# ---------------------------------------------------------------------------
# Pairwise-hit filtering for local clustering.
# ---------------------------------------------------------------------------
"""
    _read_diamond_blastp_hits_tsv(blast_output_path)

Read DIAMOND outfmt-6 blastp hits as a typed `DataFrame`.
"""
function _read_diamond_blastp_hits_tsv(blast_output_path::AbstractString)::DataFrame
    return CSV.read(
        String(blast_output_path),
        DataFrame;
        delim = '\t',
        header = _DIAMOND_BLASTP_OUTFMT_COLUMNS,
        types = _DIAMOND_BLASTP_OUTFMT_TYPES,
        normalizenames = false,
        silencewarnings = true
    )
end

"""
    _passes_local_filter(hit)

Check whether a hit satisfies the local criterion:
- identity >= `LOCAL_CLUSTER_MIN_IDENTITY`,
- ungapped alignment (`gapopen == 0`),
- aligned length >= `LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH`.
"""
function _passes_local_filter(
        hit;
        local_min_identity::Real = LOCAL_CLUSTER_MIN_IDENTITY,
        local_min_ungapped_length::Integer = LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH
)::Bool
    return hit.pident >= local_min_identity &&
           hit.gapopen == 0 &&
           hit.length >= local_min_ungapped_length
end

"""
    _local_cluster_assignments_from_hits(sequence_ids, hits)

Construct local cluster assignments by filtering `hits` and grouping accepted
edges with disjoint-set connected components.
"""
function _local_cluster_assignments_from_hits(
        sequence_ids::Vector{String},
        hits;
        local_min_identity::Real = LOCAL_CLUSTER_MIN_IDENTITY,
        local_min_ungapped_length::Integer = LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH
)::DataFrame
    sequence_index_map = _sequence_index_map(sequence_ids)
    local_disjoint_set = IntDisjointSet(length(sequence_ids))

    for hit in hits
        query_index = get(sequence_index_map, hit.qseqid, 0)
        subject_index = get(sequence_index_map, hit.sseqid, 0)
        (query_index == 0 || subject_index == 0 || query_index == subject_index) && continue

        if _passes_local_filter(
            hit;
            local_min_identity = local_min_identity,
            local_min_ungapped_length = local_min_ungapped_length
        )
            union!(local_disjoint_set, query_index, subject_index)
        end
    end

    return _build_cluster_assignments(
        "L",
        _components_from_disjoint_set!(local_disjoint_set, sequence_ids)
    )
end

# ---------------------------------------------------------------------------
# DIAMOND workflows.
# ---------------------------------------------------------------------------
"""
    _run_all_vs_all_blastp(input_fasta_path, blast_output_path; threads)

Run DIAMOND all-vs-all `blastp` for one FASTA file and write tabular output
with fields required by both clustering modes.
"""
function _run_all_vs_all_blastp(
        input_fasta_path::AbstractString,
        blast_output_path::AbstractString;
        threads::Int
)::Nothing
    mktempdir() do working_dir
        db_prefix = joinpath(working_dir, "sequence_db")
        _run_checked_command(
            _diamond_cmd([
                "makedb",
                "--in",
                String(input_fasta_path),
                "--db",
                db_prefix
            ]),
            "makedb"
        )

        blastp_args = _build_all_vs_all_blastp_args(
            input_fasta_path,
            db_prefix,
            blast_output_path;
            threads = threads
        )
        _run_checked_command(_diamond_cmd(blastp_args), "all-vs-all blastp")
    end
    return nothing
end

"""
    _local_cluster_assignments_from_blastp(input_fasta_path, sequence_ids; threads)

Run all-vs-all `blastp`, parse hits, and build local assignments.
"""
function _local_cluster_assignments_from_blastp(
        input_fasta_path::AbstractString,
        sequence_ids::Vector{String};
        threads::Int,
        local_min_identity::Real = LOCAL_CLUSTER_MIN_IDENTITY,
        local_min_ungapped_length::Integer = LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH
)::DataFrame
    return mktempdir() do working_dir
        blast_output_path = joinpath(working_dir, "diamond_all_vs_all_hits.tsv")
        _run_all_vs_all_blastp(input_fasta_path, blast_output_path; threads = threads)
        hits = _read_diamond_blastp_hits_tsv(blast_output_path)
        return _local_cluster_assignments_from_hits(
            sequence_ids,
            eachrow(hits);
            local_min_identity = local_min_identity,
            local_min_ungapped_length = local_min_ungapped_length
        )
    end
end

"""
    _cluster_assignments_from_diamond_cluster_pairs(sequence_ids, pairs)

Convert `(representative, member)` pairs from `diamond cluster` output into a
deterministic assignment table.
"""
function _cluster_assignments_from_diamond_cluster_pairs(
        sequence_ids::Vector{String},
        pairs
)::DataFrame
    sequence_index_map = _sequence_index_map(sequence_ids)
    disjoint_set = IntDisjointSet(length(sequence_ids))

    for (representative, member) in pairs
        representative_index = get(sequence_index_map, representative, 0)
        member_index = get(sequence_index_map, member, 0)
        (representative_index == 0 || member_index == 0 ||
         representative_index == member_index) &&
            continue
        union!(disjoint_set, representative_index, member_index)
    end

    return _build_cluster_assignments(
        "G",
        _components_from_disjoint_set!(disjoint_set, sequence_ids)
    )
end

"""
    _parse_diamond_cluster_pairs(cluster_output_path)

Read representative/member pairs from DIAMOND cluster tabular output.
"""
function _parse_diamond_cluster_pairs(cluster_output_path::AbstractString)
    cluster_pairs_df = CSV.read(
        String(cluster_output_path),
        DataFrame;
        delim = '\t',
        header = _DIAMOND_CLUSTER_PAIR_COLUMNS,
        types = _DIAMOND_CLUSTER_PAIR_TYPES,
        select = [1, 2],
        normalizenames = false,
        silencewarnings = true,
        comment = "#"
    )

    return collect(zip(cluster_pairs_df.representative_id, cluster_pairs_df.sequence_id))
end

"""
    _run_global_diamond_cluster(input_fasta_path, sequence_ids; threads, memory_limit)

Run `diamond cluster` + `diamond recluster` for the global criterion and
convert the refined output pairs to a deterministic assignment table.
"""
function _run_global_diamond_cluster(
        input_fasta_path::AbstractString,
        sequence_ids::Vector{String};
        threads::Int,
        global_min_identity::Real = GLOBAL_CLUSTER_MIN_IDENTITY,
        global_min_coverage::Real = GLOBAL_CLUSTER_MIN_COVERAGE,
        memory_limit::AbstractString = "16G"
)::DataFrame
    return mktempdir() do working_dir
        cluster_output_path = joinpath(working_dir, "diamond_cluster_pairs.tsv")
        recluster_output_path = joinpath(working_dir, "diamond_recluster_pairs.tsv")
        cluster_args = _build_global_cluster_args(
            input_fasta_path,
            cluster_output_path;
            threads = threads,
            memory_limit = memory_limit,
            global_min_identity = global_min_identity,
            global_min_coverage = global_min_coverage
        )

        _run_checked_command(_diamond_cmd(cluster_args), "diamond cluster")
        recluster_args = _build_global_recluster_args(
            input_fasta_path,
            cluster_output_path,
            recluster_output_path;
            threads = threads,
            memory_limit = memory_limit,
            global_min_identity = global_min_identity,
            global_min_coverage = global_min_coverage
        )
        _run_checked_command(_diamond_cmd(recluster_args), "diamond recluster")
        pairs = _parse_diamond_cluster_pairs(recluster_output_path)
        return _cluster_assignments_from_diamond_cluster_pairs(sequence_ids, pairs)
    end
end

# ---------------------------------------------------------------------------
# Public API.
# ---------------------------------------------------------------------------
"""
    run_protein_sequence_clustering(input_fasta_path;
        threads = Sys.CPU_THREADS,
        global_cluster_memory_limit = "16G",
        global_min_coverage = 80.0,
        global_min_identity = 30.0,
        local_min_identity = 80.0,
        local_min_ungapped_length = 10)

    Run protein-level clustering in two modes:
- global mode: `diamond cluster` followed by `diamond recluster` with
  configurable identity and mutual coverage,
- local mode: all-vs-all `diamond blastp` filtered with configurable
  identity and minimum ungapped length.
"""
function run_protein_sequence_clustering(
        input_fasta_path::AbstractString;
        threads::Integer = Sys.CPU_THREADS,
        global_cluster_memory_limit::AbstractString = "16G",
        global_min_coverage::Real = GLOBAL_CLUSTER_MIN_COVERAGE,
        global_min_identity::Real = GLOBAL_CLUSTER_MIN_IDENTITY,
        local_min_identity::Real = LOCAL_CLUSTER_MIN_IDENTITY,
        local_min_ungapped_length::Integer = LOCAL_CLUSTER_MIN_UNGAPPED_LENGTH
)
    thresholds = _validate_clustering_thresholds(
        global_min_coverage = global_min_coverage,
        global_min_identity = global_min_identity,
        local_min_identity = local_min_identity,
        local_min_ungapped_length = local_min_ungapped_length
    )

    @argcheck(threads > 0, "`threads` must be a positive integer.")

    sequence_ids = _read_fasta_ids(input_fasta_path)

    local_assignments = _local_cluster_assignments_from_blastp(
        input_fasta_path,
        sequence_ids;
        threads = Int(threads),
        local_min_identity = thresholds.local_min_identity,
        local_min_ungapped_length = thresholds.local_min_ungapped_length
    )
    global_assignments = _run_global_diamond_cluster(
        input_fasta_path,
        sequence_ids;
        threads = Int(threads),
        global_min_identity = thresholds.global_min_identity,
        global_min_coverage = thresholds.global_min_coverage,
        memory_limit = global_cluster_memory_limit
    )

    return (
        global_assignments = global_assignments,
        local_assignments = local_assignments,
        global_method = "diamond_cluster_recluster",
        local_method = "all_vs_all_blastp",
        sequence_count = length(sequence_ids),
        global_min_coverage = thresholds.global_min_coverage,
        global_min_identity = thresholds.global_min_identity,
        local_min_identity = thresholds.local_min_identity,
        local_min_ungapped_length = thresholds.local_min_ungapped_length
    )
end

"""
    _union_group_members!(disjoint_set, grouped_indices)

Union all indices that share the same group key.
"""
function _union_group_members!(
        disjoint_set::IntDisjointSet,
        grouped_indices::Dict{String, Vector{Int}}
)::Nothing
    for member_indices in values(grouped_indices)
        first_index = first(member_indices)
        for member_index in Iterators.drop(member_indices, 1)
            union!(disjoint_set, first_index, member_index)
        end
    end
    return nothing
end

"""
    _final_cluster_ids(protein_ids, global_representative_ids, local_representative_ids)

Build deterministic final cluster identifiers by taking the transitive closure
of co-membership across global and local representative assignments.
"""
function _final_cluster_ids(
        protein_ids::Vector{String},
        global_representative_ids::Vector{String},
        local_representative_ids::Vector{String}
)::Vector{String}
    @argcheck(length(protein_ids) == length(global_representative_ids),
        "`protein_ids` and `global_representative_ids` must have the same length.")
    @argcheck(length(protein_ids) == length(local_representative_ids),
        "`protein_ids` and `local_representative_ids` must have the same length.")

    protein_index_map = _sequence_index_map(protein_ids)
    disjoint_set = IntDisjointSet(length(protein_ids))
    global_groups = Dict{String, Vector{Int}}()
    local_groups = Dict{String, Vector{Int}}()

    for (protein_id, global_representative_id, local_representative_id) in
        zip(protein_ids, global_representative_ids, local_representative_ids)
        protein_index = protein_index_map[protein_id]
        push!(get!(global_groups, global_representative_id, Int[]), protein_index)
        push!(get!(local_groups, local_representative_id, Int[]), protein_index)
    end

    _union_group_members!(disjoint_set, global_groups)
    _union_group_members!(disjoint_set, local_groups)

    components = _components_from_disjoint_set!(disjoint_set, protein_ids)
    component_id_map = Dict{String, String}()
    for (component_index, members) in enumerate(components)
        cluster_id = "cluster_$(lpad(string(component_index), 6, '0'))"
        for protein_id in members
            component_id_map[protein_id] = cluster_id
        end
    end

    return [component_id_map[protein_id] for protein_id in protein_ids]
end

"""
    build_cluster_representatives_table(global_assignments, local_assignments)

Build a final one-row-per-protein table with columns:
- `protein_id`
- `global_representative_id`
- `local_representative_id`
- `final_cluster_id`
"""
function build_cluster_representatives_table(
        global_assignments::DataFrame,
        local_assignments::DataFrame
)::DataFrame
    @argcheck(names(global_assignments) == _CLUSTER_ASSIGNMENT_COLUMN_NAMES,
        "Global assignments must have columns $(_CLUSTER_ASSIGNMENT_COLUMN_NAMES), got $(names(global_assignments)).")
    @argcheck(names(local_assignments) == _CLUSTER_ASSIGNMENT_COLUMN_NAMES,
        "Local assignments must have columns $(_CLUSTER_ASSIGNMENT_COLUMN_NAMES), got $(names(local_assignments)).")

    global_representatives = select(
        global_assignments,
        :sequence_id => :protein_id,
        :representative_id => :global_representative_id
    )
    local_representatives = select(
        local_assignments,
        :sequence_id => :protein_id,
        :representative_id => :local_representative_id
    )

    @argcheck(nrow(global_representatives) == nrow(local_representatives),
        "Global and local assignments contain different numbers of proteins.")

    representatives = innerjoin(
        global_representatives,
        local_representatives;
        on = :protein_id,
        validate = (true, true)
    )
    @argcheck(nrow(representatives) == nrow(global_representatives),
        "Global and local assignments do not cover the same protein IDs.")

    sort!(representatives, :protein_id)
    representatives.final_cluster_id = _final_cluster_ids(
        representatives.protein_id,
        representatives.global_representative_id,
        representatives.local_representative_id
    )
    return representatives
end

"""
    write_cluster_representatives_csv(representatives, output_path)

Write the final cluster-representative table to CSV. If `output_path` ends with
`.gz`, the CSV is gzip-compressed.
"""
function write_cluster_representatives_csv(
        representatives::DataFrame,
        output_path::AbstractString
)::Nothing
    @argcheck(names(representatives) == _CLUSTER_REPRESENTATIVE_COLUMN_NAMES,
        "Representative table must have columns $(_CLUSTER_REPRESENTATIVE_COLUMN_NAMES), got $(names(representatives)).")

    @argcheck(allunique(representatives.protein_id),
        "Representative table contains duplicated `protein_id` values.")

    mkpath(dirname(String(output_path)))
    compress_output = endswith(lowercase(String(output_path)), ".gz")
    CSV.write(String(output_path), representatives; compress = compress_output)
    return nothing
end

"""
    _bytes_and_mb_string(size_in_bytes)

Render a file size as bytes plus decimal megabytes.
"""
function _bytes_and_mb_string(size_in_bytes::Integer)::String
    return "$(size_in_bytes) bytes, $(_format_mb(size_in_bytes)) MB"
end

"""
    write_sequence_clustering_readme(output_path;
        source_disprot_regions_path,
        source_disprot_regions_size_bytes,
        disprot_fasta_path,
        disprot_fasta_size_bytes,
        representatives_path,
        representatives_size_bytes,
        sequence_count,
        global_cluster_count,
        local_cluster_count,
        final_cluster_count,
        global_method,
        local_method,
        global_min_identity,
        global_min_coverage,
        local_min_identity,
        local_min_ungapped_length,
        generated_at_utc = utc_now_iso8601())

Write a run README for clustering outputs, including methods and thresholds.
"""
function write_sequence_clustering_readme(
        output_path::AbstractString;
        source_disprot_regions_path::AbstractString,
        source_disprot_regions_size_bytes::Integer,
        disprot_fasta_path::AbstractString,
        disprot_fasta_size_bytes::Integer,
        representatives_path::AbstractString,
        representatives_size_bytes::Integer,
        sequence_count::Integer,
        global_cluster_count::Integer,
        local_cluster_count::Integer,
        final_cluster_count::Integer,
        global_method::AbstractString,
        local_method::AbstractString,
        global_min_identity::Real,
        global_min_coverage::Real,
        local_min_identity::Real,
        local_min_ungapped_length::Integer,
        generated_at_utc::AbstractString = utc_now_iso8601()
)::Nothing
    thresholds = _validate_clustering_thresholds(
        global_min_identity = global_min_identity,
        global_min_coverage = global_min_coverage,
        local_min_identity = local_min_identity,
        local_min_ungapped_length = local_min_ungapped_length
    )

    @argcheck(sequence_count >= 1, "`sequence_count` must be >= 1.")
    @argcheck(global_cluster_count >= 1, "`global_cluster_count` must be >= 1.")
    @argcheck(local_cluster_count >= 1, "`local_cluster_count` must be >= 1.")
    @argcheck(final_cluster_count >= 1, "`final_cluster_count` must be >= 1.")
    @argcheck(source_disprot_regions_size_bytes >= 0,
        "`source_disprot_regions_size_bytes` must be >= 0.")
    @argcheck(disprot_fasta_size_bytes >= 0, "`disprot_fasta_size_bytes` must be >= 0.")
    @argcheck(representatives_size_bytes >= 0,
        "`representatives_size_bytes` must be >= 0.")

    readme_text = """
# Sequence Clustering Outputs

Generated at (UTC): $(generated_at_utc)

## Outputs
- DisProt full-protein FASTA: $(disprot_fasta_path) ($(_bytes_and_mb_string(disprot_fasta_size_bytes)))
- Cluster representatives CSV: $(representatives_path) ($(_bytes_and_mb_string(representatives_size_bytes)))

## Sequence and Cluster Counts
- Protein sequence count: $(sequence_count)
- Global cluster count: $(global_cluster_count)
- Local cluster count: $(local_cluster_count)
- Final cluster count: $(final_cluster_count)

## Clustering Methods
- Global method: $(global_method)
- Global thresholds: identity >= $(thresholds.global_min_identity)% and mutual coverage >= $(thresholds.global_min_coverage)%
- Local method: $(local_method)
- Local thresholds: identity >= $(thresholds.local_min_identity)% on ungapped alignments with length >= $(thresholds.local_min_ungapped_length) aa

## Source Data
- DisProt regions JSON: $(source_disprot_regions_path) ($(_bytes_and_mb_string(source_disprot_regions_size_bytes)))
"""

    _write_text_file(output_path, readme_text)
    return nothing
end
