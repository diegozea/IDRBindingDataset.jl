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
  2. "Binding site" features are extracted (absolute coordinates on the
     full-length `sequence`).
  3. Adds to the CSV:
       - binding_in_uniprot          :: Bool   (protein-level, repeated)
       - residues_binding_in_uniprot :: Vector{String}  (protein-level)
       - publication_mention_binding :: Vector{Int}      (protein-level)
       - ADAPTED_LABELS_UNIPROT      :: String  (PDB-ROW level)

  ADAPTED_LABELS_UNIPROT — CORRECTED CONSTRUCTION (2024 fix):
     Earlier versions of this script rebuilt a full-length "corrected"
     label vector from `sequence`/`residue_labels` and then re-sliced it
     with SP_BEG:SP_END. Any inconsistency at the *protein* level (missing
     `residue_labels`, length mismatch with `sequence`, a mismatched
     `sequence[SP_BEG:SP_END]` vs `ADAPTED_SEQUENCE`, ...) silently made
     the WHOLE row `missing`, even though `ADAPTED_LABELS` itself was
     perfectly valid. That could make ADAPTED_LABELS_UNIPROT shorter than
     ADAPTED_LABELS, or show fewer class-2 residues overall than
     ADAPTED_LABELS — which should be structurally impossible.

     This version instead:
       - parses THIS row's own `ADAPTED_LABELS` (already the correct
         length, already consistent with ADAPTED_SEQUENCE by construction
         upstream),
       - converts each UniProt binding site's ABSOLUTE coordinates
         (on `sequence`) into LOCAL coordinates within this row's adapted
         region using SP_BEG,
       - recodes ONLY residues currently at 1 (disorder) that fall in a
         binding site to 2. Residues at 0 or already at 2 are never
         touched.
     ADAPTED_LABELS_UNIPROT is therefore GUARANTEED to have exactly the
     same length as ADAPTED_LABELS, and its count of class-2 residues can
     only stay equal or increase, never decrease — as long as
     ADAPTED_LABELS itself is present. `sequence`/`residue_labels` are
     only used for the protein-level fields (binding_in_uniprot,
     residues_binding_in_uniprot, publication_mention_binding), which are
     informative flags, not something the row-level label vector depends
     on anymore.

Usage:
    julia 03_add_uniprot_informations.jl input.csv output.csv [batch_size]

Dependencies: CSV, DataFrames, HTTP, JSON3, CodecZlib
    julia --project -e 'import Pkg; Pkg.add(["CSV","DataFrames","HTTP","JSON3","CodecZlib"])'
"""

import Pkg

const DEFAULT_ENV = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDP/scripts/env/env1_020226"

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
residue (e.g. "0001120..."), NOT a comma-separated "[0, 0, 1]" list.

Defensive fix: if `s` is not a clean digit-only string (e.g. the literal
text "missing" coming from `string(missing)` on an actually-missing value
upstream, or any other unexpected content), this returns an EMPTY vector
instead of silently computing garbage numbers from ASCII codes (which was
the original, undetected source of `length(seq) != length(labels)` for
some proteins).
"""
function parse_digit_string_labels(s::AbstractString)
    cleaned = strip(s)
    isempty(cleaned) && return Int[]
    if !all(isdigit, cleaned)
        return Int[]
    end
    return [Int(c) - Int('0') for c in cleaned]
end

"""
Parses `ADAPTED_LABELS` / `ADAPTED_LABELS_UNIPROT`-style strings written as
a Julia/Python list: "[0, 0, 1, 2]". Returns an empty vector for anything
that isn't a well-formed bracketed, comma-separated list of integers.
"""
function parse_adapted_labels_list(s::AbstractString)
    cleaned = strip(s)
    isempty(cleaned) && return Int[]
    i = findfirst('[', cleaned)
    j = findlast(']', cleaned)
    (i === nothing || j === nothing || j <= i) && return Int[]
    inner = strip(cleaned[i+1:j-1])
    isempty(inner) && return Int[]
    parts = split(inner, ",")
    out = Int[]
    for p in parts
        p = strip(p)
        isempty(p) && continue
        v = tryparse(Int, p)
        v === nothing && return Int[]   # contenu invalide -> on ne devine pas, on renvoie vide
        push!(out, v)
    end
    return out
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
#    Coordinates are ABSOLUTE positions on the full-length `sequence`
#    (1-indexed, UniProt convention).
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
#    (used only for the protein-level informative fields)
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

    # --- Diagnostics counters (printed at the end, for sanity-checking) ---
    n_protein_seqlabel_mismatch    = 0  # protein-level: sequence/residue_labels incohérents (n'affecte QUE has_binding/residues_binding_in_uniprot)
    n_rows_no_sites_unchanged      = 0  # ADAPTED_LABELS_UNIPROT == ADAPTED_LABELS (aucun site UniProt disponible / applicable)
    n_rows_recoded                 = 0  # au moins un résidu recodé 1 -> 2
    n_rows_missing_adapted_labels  = 0  # ADAPTED_LABELS lui-même manquant/invalide -> ADAPTED_LABELS_UNIPROT reste missing (cas légitime, pas un bug)

    # --- Step 3: per protein_id ---
    for (gi, subdf) in enumerate(grouped)
        pid = string(subdf.protein_id[1])
        accession = pid_to_accession[pid]
        entry = get(accession_to_entry, uppercase(accession), nothing)
        sites = extract_binding_sites(entry)   # coordonnées ABSOLUES sur `sequence`

        # -- Champs informatifs au niveau protéine --
        # Purement indicatifs : une incohérence ici n'empêche plus jamais le
        # calcul de ADAPTED_LABELS_UNIPROT (voir plus bas).
        seq = ismissing(subdf.sequence[1]) ? "" : string(subdf.sequence[1])
        raw_labels_str = ismissing(subdf.residue_labels[1]) ? "" : string(subdf.residue_labels[1])
        labels = parse_digit_string_labels(raw_labels_str)

        if !isempty(seq) && !isempty(labels) && length(seq) != length(labels)
            @warn "sequence et residue_labels de longueurs différentes (n'affecte que has_binding_in_uniprot / residues_binding_in_uniprot pour cette protéine)" pid length(seq) length(labels)
            n_protein_seqlabel_mismatch += 1
        end

        kept_residues = String[]
        kept_pubs = Int[]
        if length(seq) == length(labels) && !isempty(labels)
            for s in sites
                if overlaps_disorder_label(s.start, s.stop, labels)
                    push!(kept_residues, format_residue_range(s.start, s.stop))
                    push!(kept_pubs, s.n_pubmed)
                end
            end
        end
        has_binding = !isempty(kept_residues)

        for row_idx in pid_to_rowidx[pid]
            binding_in_uniprot[row_idx] = has_binding
            residues_binding_in_uniprot[row_idx] = kept_residues
            publication_mention_binding[row_idx] = kept_pubs

            # --- ADAPTED_LABELS_UNIPROT : construit à partir du ADAPTED_LABELS
            # DE CETTE LIGNE (déjà correct/cohérent), jamais reconstruit depuis
            # sequence/residue_labels. Garantit même longueur, jamais moins de
            # classe 2 que ADAPTED_LABELS.
            raw_adapted_labels = df.ADAPTED_LABELS[row_idx]
            if ismissing(raw_adapted_labels) || isempty(strip(string(raw_adapted_labels)))
                adapted_labels_uniprot[row_idx] = missing
                n_rows_missing_adapted_labels += 1
                continue
            end

            out = parse_adapted_labels_list(string(raw_adapted_labels))
            if isempty(out)
                adapted_labels_uniprot[row_idx] = missing
                n_rows_missing_adapted_labels += 1
                continue
            end

            sp_beg_raw = df.SP_BEG[row_idx]
            sp_end_raw = df.SP_END[row_idx]
            recoded_here = false
            if !ismissing(sp_beg_raw) && !ismissing(sp_end_raw) && !isempty(sites)
                sp_beg = Int(sp_beg_raw)
                sp_end = Int(sp_end_raw)
                if (sp_end - sp_beg + 1) != length(out)
                    @warn "SP_END-SP_BEG+1 ne correspond pas à la longueur de ADAPTED_LABELS pour cette ligne ; correction UniProt appliquée sur la longueur réelle de ADAPTED_LABELS" pid row_idx sp_beg sp_end length(out)
                end
                for s in sites
                    lo = max(s.start, sp_beg)
                    hi = min(s.stop, sp_end)
                    lo > hi && continue
                    for g in lo:hi
                        local_idx = g - sp_beg + 1
                        (local_idx < 1 || local_idx > length(out)) && continue
                        if out[local_idx] == 1
                            out[local_idx] = 2
                            recoded_here = true
                        end
                    end
                end
            end
            recoded_here ? (n_rows_recoded += 1) : (n_rows_no_sites_unchanged += 1)

            adapted_labels_uniprot[row_idx] = "[" * join(out, ", ") * "]"
        end

        if gi % 200 == 0
            println("  ... $gi/$n_groups protein_id traités")
        end
    end

    println()
    println("--- Diagnostics ---")
    println("Protéines avec sequence/residue_labels incohérents (champs protéine seulement) : $n_protein_seqlabel_mismatch")
    println("Lignes ADAPTED_LABELS_UNIPROT recodées (≥1 résidu 1->2)                          : $n_rows_recoded")
    println("Lignes ADAPTED_LABELS_UNIPROT = ADAPTED_LABELS (aucun site applicable)           : $n_rows_no_sites_unchanged")
    println("Lignes ADAPTED_LABELS_UNIPROT manquant (ADAPTED_LABELS lui-même absent/invalide) : $n_rows_missing_adapted_labels")

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
        println("Usage: julia 03_add_uniprot_informations.jl input.csv output.csv [batch_size]")
        exit(1)
    end
    bsize = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : DEFAULT_BATCH_SIZE
    enrich_csv(ARGS[1], ARGS[2]; batch_size=bsize)
end