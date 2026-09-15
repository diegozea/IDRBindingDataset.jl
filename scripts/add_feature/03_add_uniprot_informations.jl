#!/usr/bin/env julia

"""

Takes as input a CSV with (among others) the columns:
    protein_id, protein_disprot_id, sequence, residue_labels,
    ADAPTED_SEQUENCE, ADAPTED_LABELS, SP_BEG, SP_END, ...

`residue_labels` is a string of concatenated digits, one character per
residue of `sequence` (e.g. "0000111220000..."), not a Python-style list.
`ADAPTED_LABELS`, on the other hand, is in list format: "[0, 0, 1, 2, ...]",
aligned with `ADAPTED_SEQUENCE` (a substring of `sequence`, one row per PDB
entry for a given protein_id). ADAPTED_SEQUENCE is built as
sequence[SP_BEG:SP_END] (1-indexed, inclusive).
Labels: 0 = ordered, 1 = disordered, 2 = disordered-and-binding.

Steps:
  0. Rows with a missing `protein_disprot_id` are dropped.
  1. For each UNIQUE UniProt `protein_id` (a single HTTP request per
     protein, regardless of how many PDB rows are associated with it), the
     corresponding UniProt entry is fetched via BATCH requests.
  2. "Binding site" features are extracted.
  3. Per protein_id, a full-length "corrected" label vector is computed
     (same length as `sequence`): any residue labelled 1 (disorder) that
     overlaps a UniProt binding site becomes 2 (disorder-and-binding).
  4. Adds to the CSV:
       - binding_in_uniprot          :: Bool   (protein-level, repeated)
       - residues_binding_in_uniprot :: Vector{String}  (protein-level)
       - publication_mention_binding :: Vector{Int}      (protein-level)
       - ADAPTED_LABELS_UNIPROT      :: String  (PDB-ROW level: this is the
         sub-vector of the corrected label vector corresponding to
         sequence[SP_BEG:SP_END] for THIS row, i.e. the same range used to
         build ADAPTED_SEQUENCE/ADAPTED_LABELS. Consistency between
         sequence[SP_BEG:SP_END] and ADAPTED_SEQUENCE is checked before
         use, see `sp_range_matches_adapted`.)

Usage:
    julia 03_add_uniprot_informations.jl input.csv output.csv [batch_size]

Dependencies: CSV, DataFrames, HTTP, JSON3, CodecZlib
    julia --project -e 'import Pkg; Pkg.add(["CSV","DataFrames","HTTP","JSON3","CodecZlib"])'
"""

import Pkg


if abspath(PROGRAM_FILE) == @__FILE__
    env_idx = findfirst(==("--julia-env"), ARGS)
    julia_env = env_idx !== nothing && env_idx < length(ARGS) ? ARGS[env_idx+1] : DEFAULT_ENV
    Pkg.activate(julia_env)
end
Pkg.instantiate()


using CSV
using DataFrames
using HTTP
using JSON3
using CodecZlib

const SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
const REQUEST_DELAY = 0.3     # seconds between batch requests
const MAX_RETRIES = 3
const DEFAULT_BATCH_SIZE = 90  # accessions per request (margin under the URL length limit)
const PAGE_SIZE = 500          # UniProt-side size= (max allowed by the API)

# ---------------------------------------------------------------------------
# 1. Label parsing
# ---------------------------------------------------------------------------
"""
Converts the HTTP response body to text, manually decompressing it if the
bytes start with the gzip signature (0x1f 0x8b) — some server/HTTP.jl
combinations let a still-compressed body through even when requesting
Accept-Encoding: identity.
"""
function response_text(resp)
    body = resp.body
    if length(body) >= 2 && body[1] == 0x1f && body[2] == 0x8b
        return String(transcode(GzipDecompressor, body))
    else
        return String(body)
    end
end

"""
Parses `residue_labels`: a string of concatenated digits, one character per
residue (e.g. "0001120..."), NOT a comma-separated "[0, 0, 1]" list. This
mismatch was the cause of the original OverflowError: the previous parser
tried to read the whole string as a single Int64.
"""
function parse_digit_string_labels(s::AbstractString)
    cleaned = strip(s)
    isempty(cleaned) && return Int[]
    return [Int(c) - Int('0') for c in cleaned]
end

# ---------------------------------------------------------------------------
# 2. UniProt batch request: returns Dict{accession_uppercase => entry_json}
#    Handles internal pagination if a batch exceeds PAGE_SIZE results
#    (shouldn't happen with batch_size <= 500, but included as a safeguard).
# ---------------------------------------------------------------------------
function fetch_uniprot_batch(accessions::Vector{String})
    query = join(("accession:$acc" for acc in accessions), " OR ")
    results = Dict{String, Any}()

    cursor_url = SEARCH_URL
    query_params = Dict(
        "query" => query,
        "format" => "json",
        "size" => string(PAGE_SIZE),
    )
    req_headers = ["Accept-Encoding" => "identity"]  # avoids gzip being mis-decompressed by HTTP.jl

    while true
        local resp
        success = false
        for attempt in 1:MAX_RETRIES
            try
                resp = HTTP.get(cursor_url, req_headers; query=query_params, status_exception=false)
                if resp.status == 200
                    success = true
                    break
                else
                    @warn "Réponse HTTP inattendue sur batch" status=resp.status attempt
                end
            catch e
                @warn "Erreur réseau sur batch, nouvelle tentative" attempt exception=e
            end
            sleep(REQUEST_DELAY * attempt)
        end

        if !success
            @warn "Échec définitif pour ce batch" n_accessions=length(accessions)
            return results
        end

        parsed = JSON3.read(response_text(resp))
        if haskey(parsed, :results)
            for entry in parsed.results
                acc = uppercase(string(entry.primaryAccession))
                results[acc] = entry
            end
        end

        # Pagination via the Link header (rel="next"), just in case
        next_link = nothing
        for h in resp.headers
            if lowercase(h[1]) == "link" && occursin("rel=\"next\"", h[2])
                m = match(r"<([^>]+)>", h[2])
                if m !== nothing
                    next_link = m.captures[1]
                end
            end
        end

        if next_link === nothing
            break
        else
            cursor_url = next_link
            query_params = Dict{String,String}()  # already encoded in next_link
        end
    end

    return results
end

# ---------------------------------------------------------------------------
# 3. Extracting binding sites from a UniProt JSON entry
# ---------------------------------------------------------------------------
function extract_binding_sites(entry)
    sites = NamedTuple{(:start, :stop, :n_pubmed), Tuple{Int,Int,Int}}[]
    entry === nothing && return sites
    !haskey(entry, :features) && return sites

    for feature in entry.features
        get(feature, :type, "") == "Binding site" || continue

        loc = feature.location
        start_val = loc.start.value
        end_val = loc.end.value

        n_pubmed = 0
        if haskey(feature, :evidences)
            for ev in feature.evidences
                if get(ev, :source, "") == "PubMed"
                    n_pubmed += 1
                end
            end
        end

        push!(sites, (start=start_val, stop=end_val, n_pubmed=n_pubmed))
    end
    return sites
end

# ---------------------------------------------------------------------------
# 4. Overlap between binding site and residues labelled 1 or 2
# ---------------------------------------------------------------------------
function overlaps_disorder_label(site_start::Int, site_stop::Int, labels::Vector{Int})
    n = length(labels)
    lo = max(site_start, 1)
    hi = min(site_stop, n)
    lo > hi && return false
    for i in lo:hi
        if labels[i] == 1 || labels[i] == 2
            return true
        end
    end
    return false
end

format_residue_range(a::Int, b::Int) = a == b ? string(a) : "$(a)-$(b)"

"""
Reproduces the label vector over the full length of `sequence`, with
disordered residues (label 1) that overlap a UniProt binding site recoded
to 2 (disorder-and-binding). Residues already at 0 are never changed, even
if they fall within a binding site's range.
"""
function labels_with_uniprot_binding(labels::Vector{Int}, sites)
    n = length(labels)
    out = copy(labels)
    for s in sites
        lo = max(s.start, 1)
        hi = min(s.stop, n)
        lo > hi && continue
        for i in lo:hi
            if out[i] == 1
                out[i] = 2
            end
        end
    end
    return out
end

"""
Checks that `sequence[sp_beg:sp_end]` matches `adapted_sequence` (this is
the convention used to build ADAPTED_SEQUENCE/ADAPTED_LABELS from
SP_BEG/SP_END). Returns `true` if consistent, `false` otherwise
(out-of-bounds range, length mismatch, or different content).
"""
function sp_range_matches_adapted(sequence::AbstractString, adapted_sequence::AbstractString,
                                   sp_beg::Int, sp_end::Int)
    n = length(sequence)
    (sp_beg < 1 || sp_end > n || sp_beg > sp_end) && return false
    (sp_end - sp_beg + 1) != length(adapted_sequence) && return false
    return sequence[sp_beg:sp_end] == adapted_sequence
end

# ---------------------------------------------------------------------------
# 5. Main processing
# ---------------------------------------------------------------------------
function enrich_csv(input_path::AbstractString, output_path::AbstractString;
                     batch_size::Int = DEFAULT_BATCH_SIZE)
    df = CSV.read(input_path, DataFrame; missingstring=["", "NA", "NaN"])

    required_cols = ["protein_id", "protein_disprot_id", "sequence", "residue_labels",
                      "ADAPTED_SEQUENCE", "ADAPTED_LABELS", "SP_BEG", "SP_END"]
    for col in required_cols
        @assert col in names(df) "Colonne manquante dans le CSV : $col"
    end

    # --- Filter: drop rows with no DisProt annotation ---
    n_before = nrow(df)
    df = filter(:protein_disprot_id => x -> !ismissing(x) && !isempty(strip(string(x))), df)
    println("Lignes supprimées (protein_disprot_id manquant) : $(n_before - nrow(df))")

    n = nrow(df)
    binding_in_uniprot          = Vector{Union{Missing,Bool}}(missing, n)
    residues_binding_in_uniprot = Vector{Union{Missing,Vector{String}}}(missing, n)
    publication_mention_binding = Vector{Union{Missing,Vector{Int}}}(missing, n)
    adapted_labels_uniprot      = Vector{Union{Missing,String}}(missing, n)

    grouped = groupby(df, :protein_id)
    n_groups = length(grouped)
    println("Nombre de protein_id (UniProt) uniques : $n_groups")

    # --- Step 1: a single UniProt identifier per protein_id ---
    pid_to_accession = Dict{String,String}()   # protein_id -> UniProt accession
    pid_to_rowidx    = Dict{String,Vector{Int}}()

    for subdf in grouped
        pid = string(subdf.protein_id[1])
        accession = if hasproperty(df, :UNIPROT) &&
                       !ismissing(subdf.UNIPROT[1]) &&
                       !isempty(strip(string(subdf.UNIPROT[1])))
            strip(string(subdf.UNIPROT[1]))
        else
            pid
        end
        pid_to_accession[pid] = accession
        pid_to_rowidx[pid] = collect(parentindices(subdf)[1])
    end

    unique_accessions = unique(collect(values(pid_to_accession)))
    println("Nombre d'accessions UniProt uniques à interroger : $(length(unique_accessions))")

    # --- Step 2: batch fetching (one request per unique accession) ---
    accession_to_entry = Dict{String, Any}()
    n_batches = cld(length(unique_accessions), batch_size)
    for (bi, chunk_start) in enumerate(1:batch_size:length(unique_accessions))
        chunk = unique_accessions[chunk_start:min(chunk_start + batch_size - 1, end)]
        println("[Batch $bi/$n_batches] $(length(chunk)) accessions...")
        batch_results = fetch_uniprot_batch(chunk)
        for (acc, entry) in batch_results
            accession_to_entry[acc] = entry
        end
        missing_in_batch = [a for a in chunk if !haskey(batch_results, uppercase(a))]
        if !isempty(missing_in_batch)
            @warn "Accessions non trouvées dans ce batch" missing_in_batch
        end
        sleep(REQUEST_DELAY)
    end

    # --- Step 3: per protein_id, compute the UniProt info once, then
    #     derive per-PDB-row values from it (ADAPTED_LABELS_UNIPROT) ---
    for (gi, subdf) in enumerate(grouped)
        pid = string(subdf.protein_id[1])
        accession = pid_to_accession[pid]
        entry = get(accession_to_entry, uppercase(accession), nothing)

        seq = string(subdf.sequence[1])
        labels = parse_digit_string_labels(string(subdf.residue_labels[1]))

        if length(seq) != length(labels)
            @warn "sequence et residue_labels de longueurs différentes, protein_id ignoré pour l'enrichissement" pid length(seq) length(labels)
        end

        sites = extract_binding_sites(entry)

        kept_residues = String[]
        kept_pubs = Int[]
        corrected_labels = Int[]
        if length(seq) == length(labels)
            for s in sites
                if overlaps_disorder_label(s.start, s.stop, labels)
                    push!(kept_residues, format_residue_range(s.start, s.stop))
                    push!(kept_pubs, s.n_pubmed)
                end
            end
            corrected_labels = labels_with_uniprot_binding(labels, sites)
        end

        has_binding = !isempty(kept_residues)
        for row_idx in pid_to_rowidx[pid]
            binding_in_uniprot[row_idx] = has_binding
            residues_binding_in_uniprot[row_idx] = kept_residues
            publication_mention_binding[row_idx] = kept_pubs

            # --- ADAPTED_LABELS_UNIPROT: specific to THIS PDB row ---
            # ADAPTED_SEQUENCE = sequence[SP_BEG:SP_END] (same convention
            # used to build ADAPTED_LABELS).
            raw_adapted = df.ADAPTED_SEQUENCE[row_idx]
            adapted_seq = ismissing(raw_adapted) ? "" : string(raw_adapted)
            sp_beg_raw = df.SP_BEG[row_idx]
            sp_end_raw = df.SP_END[row_idx]

            if isempty(adapted_seq) || isempty(corrected_labels) ||
               ismissing(sp_beg_raw) || ismissing(sp_end_raw)
                adapted_labels_uniprot[row_idx] = missing
                continue
            end

            sp_beg = Int(sp_beg_raw)
            sp_end = Int(sp_end_raw)

            if !sp_range_matches_adapted(seq, adapted_seq, sp_beg, sp_end)
                @warn "sequence[SP_BEG:SP_END] ne correspond pas à ADAPTED_SEQUENCE, ligne ignorée" pid row_idx sp_beg sp_end
                adapted_labels_uniprot[row_idx] = missing
                continue
            end

            sub = corrected_labels[sp_beg:sp_end]
            adapted_labels_uniprot[row_idx] = "[" * join(sub, ", ") * "]"
        end

        if gi % 200 == 0
            println("  ... $gi/$n_groups protein_id traités")
        end
    end

    df.binding_in_uniprot = binding_in_uniprot
    df.residues_binding_in_uniprot = residues_binding_in_uniprot
    df.publication_mention_binding = publication_mention_binding
    df.ADAPTED_LABELS_UNIPROT = adapted_labels_uniprot

    CSV.write(output_path, df)
    println("\nTerminé. Fichier écrit : $output_path")
    return df
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: julia uniprot_binding_enrich.jl input.csv output.csv [batch_size]")
        exit(1)
    end
    bsize = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : DEFAULT_BATCH_SIZE
    enrich_csv(ARGS[1], ARGS[2]; batch_size=bsize)
end