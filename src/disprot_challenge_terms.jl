using CodecZlib
using JSON
using OBOParse

const _GO_BASIC_RELATIVE_PATH = joinpath("go", "go-basic.obo.gz")
const _DISPROT_IDPO_RELATIVE_PATH = joinpath("disprot", "idpo.json.gz")
const _GO_BINDING_ROOT_TERM_ID = "GO:0005488"
const _GO_BINDING_NUCLEIC_ACID_ROOT_TERM_ID = "GO:0003676"
const _GO_BINDING_PROTEIN_ROOT_TERM_ID = "GO:0005515"
const _IDPO_DISORDER_ROOT_TERM_ID = "IDPO:0000002"

"""
    _DisProtChallengeTermRules

Ontology-derived term sets used to classify DisProt regions into challenge
categories without relying on textual heuristics.
"""
struct _DisProtChallengeTermRules
    go_binding_terms::Set{String}
    go_binding_nucleic_acid_terms::Set{String}
    go_binding_protein_terms::Set{String}
    idpo_disorder_terms::Set{String}
end

"""
    _read_gzip_text_file(path)

Read one gzip-compressed UTF-8 text file and return the full payload string.
"""
function _read_gzip_text_file(path::AbstractString)::String
    isfile(path) || throw(ArgumentError("File does not exist: $(path)"))
    open(path) do io
        read(GzipDecompressorStream(io), String)
    end
end

"""
    _region_term_id(region)

Return one normalized region `term_id` or `nothing` when missing/empty.
"""
function _region_term_id(region)::Union{Nothing, String}
    region isa AbstractDict || return nothing
    term_id = get(region, "term_id", missing)
    term_id isa AbstractString || return nothing

    normalized = strip(String(term_id))
    isempty(normalized) && return nothing
    normalized
end

"""
    _go_is_a_children_map_from_obo_gz(go_obo_gz_path)

Build a parent->children map for GO terms from one `go-basic.obo.gz` file using
`OBOParse`.
"""
function _go_is_a_children_map_from_obo_gz(
        go_obo_gz_path::AbstractString
)::Dict{String, Set{String}}
    ontology = mktemp() do obo_path, io
        write(io, _read_gzip_text_file(go_obo_gz_path))
        close(io)
        OBOParse.load(obo_path, "GO")
    end
    hasproperty(ontology, :terms) || throw(
        ArgumentError(
        "Parsed GO ontology object does not expose `terms`: $(typeof(ontology))."
    ),
    )
    terms = getproperty(ontology, :terms)
    terms isa AbstractDict || throw(
        ArgumentError(
        "Parsed GO ontology `terms` must be an associative mapping; found $(typeof(terms))."
    ),
    )

    children = Dict{String, Set{String}}()
    for (child_id_any, term_any) in pairs(terms)
        child_id_any isa AbstractString || continue
        child_id = String(child_id_any)
        startswith(child_id, "GO:") || continue
        hasproperty(term_any, :obsolete) || continue
        getproperty(term_any, :obsolete) && continue
        hasproperty(term_any, :relationships) || continue

        relationships = getproperty(term_any, :relationships)
        relationships isa AbstractDict || continue
        parents = get(relationships, :is_a, Set{String}())
        for parent_id_any in parents
            parent_id_any isa AbstractString || continue
            parent_id = String(parent_id_any)
            startswith(parent_id, "GO:") || continue

            parent_term = get(terms, parent_id, nothing)
            parent_term === nothing && continue
            hasproperty(parent_term, :obsolete) || continue
            getproperty(parent_term, :obsolete) && continue

            push!(get!(children, parent_id, Set{String}()), child_id)
        end
    end

    children
end

"""
    _idpo_curie_from_any_term_identifier(value)

Normalize one IDPO identifier from either CURIE (`IDPO:...`) or OBO IRI
(`.../IDPO_...`) to CURIE form `IDPO:xxxxxxx`.
"""
function _idpo_curie_from_any_term_identifier(
        value
)::Union{Nothing, String}
    value isa AbstractString || return nothing
    text = strip(String(value))
    isempty(text) && return nothing

    if startswith(text, "IDPO:")
        return text
    end

    marker_index = findlast("IDPO_", text)
    marker_index === nothing && return nothing
    suffix = text[(marker_index.stop + 1):end]
    isempty(suffix) && return nothing
    "IDPO:$(suffix)"
end

"""
    _idpo_is_a_children_map_from_json_gz(idpo_json_gz_path)

Build one parent->children map for IDPO terms by parsing `is_a` edges from
`idpo.json.gz` OBO-JSON.
"""
function _idpo_is_a_children_map_from_json_gz(
        idpo_json_gz_path::AbstractString
)::Dict{String, Set{String}}
    payload = JSON.parse(_read_gzip_text_file(idpo_json_gz_path))
    payload isa AbstractDict || throw(
        ArgumentError("IDPO JSON payload must be an object: $(idpo_json_gz_path)"),
    )

    graphs = get(payload, "graphs", missing)
    graphs isa AbstractVector || throw(
        ArgumentError("IDPO JSON payload must include vector key `graphs`."),
    )
    isempty(graphs) && throw(ArgumentError("IDPO JSON payload has no graphs."))

    graph = first(graphs)
    graph isa AbstractDict || throw(
        ArgumentError("IDPO JSON graph entry must be an object."),
    )
    edges = get(graph, "edges", missing)
    edges isa AbstractVector || throw(
        ArgumentError("IDPO JSON graph must include vector key `edges`."),
    )

    children = Dict{String, Set{String}}()
    for edge_any in edges
        edge_any isa AbstractDict || continue
        pred = get(edge_any, "pred", missing)
        pred == "is_a" || continue

        parent = _idpo_curie_from_any_term_identifier(get(edge_any, "obj", missing))
        child = _idpo_curie_from_any_term_identifier(get(edge_any, "sub", missing))
        if parent === nothing || child === nothing
            continue
        end

        push!(get!(children, parent, Set{String}()), child)
    end

    children
end

"""
    _descendants_from_children_map(children_map, root_term_id)

Collect one root term plus all descendants reachable via repeated child edges.
"""
function _descendants_from_children_map(
        children_map::Dict{String, Set{String}},
        root_term_id::AbstractString
)::Set{String}
    descendants = Set{String}([String(root_term_id)])
    stack = String[String(root_term_id)]

    while !isempty(stack)
        parent = pop!(stack)
        for child in get(children_map, parent, Set{String}())
            if child in descendants
                continue
            end
            push!(descendants, child)
            push!(stack, child)
        end
    end

    descendants
end

"""
    _build_disprot_challenge_term_rules(raw_data_dir)

Build ontology-derived challenge term sets from:
- `raw_data_dir/go/go-basic.obo.gz`
- `raw_data_dir/disprot/idpo.json.gz`

Semantics follow the CAID reference notebook:
- GO: `binding`/subtypes are built from `is_a` closure.
- IDPO disorder uses the structural-state branch rooted at `IDPO:0000002`.
"""
function _build_disprot_challenge_term_rules(
        raw_data_dir::AbstractString
)::_DisProtChallengeTermRules
    go_obo_gz_path = joinpath(raw_data_dir, _GO_BASIC_RELATIVE_PATH)
    idpo_json_gz_path = joinpath(raw_data_dir, _DISPROT_IDPO_RELATIVE_PATH)
    isfile(go_obo_gz_path) || throw(
        ArgumentError("GO ontology file is missing: $(go_obo_gz_path)"),
    )
    isfile(idpo_json_gz_path) || throw(
        ArgumentError("IDPO ontology file is missing: $(idpo_json_gz_path)"),
    )

    go_children_map = _go_is_a_children_map_from_obo_gz(go_obo_gz_path)
    idpo_children_map = _idpo_is_a_children_map_from_json_gz(idpo_json_gz_path)

    go_binding_terms = _descendants_from_children_map(
        go_children_map,
        _GO_BINDING_ROOT_TERM_ID
    )
    go_binding_nucleic_acid_terms = _descendants_from_children_map(
        go_children_map,
        _GO_BINDING_NUCLEIC_ACID_ROOT_TERM_ID
    )
    go_binding_protein_terms = _descendants_from_children_map(
        go_children_map,
        _GO_BINDING_PROTEIN_ROOT_TERM_ID
    )
    idpo_disorder_terms = _descendants_from_children_map(
        idpo_children_map,
        _IDPO_DISORDER_ROOT_TERM_ID
    )

    _DisProtChallengeTermRules(
        go_binding_terms,
        go_binding_nucleic_acid_terms,
        go_binding_protein_terms,
        idpo_disorder_terms
    )
end
