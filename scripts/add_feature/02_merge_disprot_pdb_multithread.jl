#!/usr/bin/env julia -t 19

#Same as code 1 but use multithread and save at each step 

#Merge Disport informations created before with PDB information merge on the uniprot ID 
#Merge ligne by Chain if only informations that change 
#Create Label at PDB level + add informations about Missing residues 

import Pkg

if abspath(PROGRAM_FILE) == @__FILE__
    env_idx = findfirst(==("--julia-env"), ARGS)
    julia_env = env_idx !== nothing && env_idx < length(ARGS) ? ARGS[env_idx+1] : DEFAULT_ENV
    Pkg.activate(julia_env)
end
Pkg.instantiate()

using CSV
using DataFrames
using CodecZlib
using MIToS
using MIToS.SIFTS
using Dates
using Base.Threads

############################ UTILITAIRES ############################

now_str() = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

# Verrou global pour éviter que deux threads ne téléchargent le même fichier
# SIFTS en même temps (rare, mais possible si le même PDB apparaît dans deux
# groupes protéine différents traités en parallèle).
const DOWNLOAD_LOCK = ReentrantLock()
const JULIA_NUM_THREADS=19
# ----------------------------------------------------------------------
# Parsing / label utilities
# ----------------------------------------------------------------------

"""
    parse_residue_labels(s::String) -> Vector{Int}

Convertit une chaîne de chiffres (ex. `"0102201"`) en vecteur d'entiers.
C'est le format du CSV brut d'entrée.
"""
function parse_residue_labels(s::String)
    return [parse(Int, c) for c in s]
end

"""
    parse_residue_labels_field(s::AbstractString) -> Vector{Int}

Parseur robuste utilisé quand on RELIT un CSV déjà écrit par ce pipeline :
`CSV.write` sérialise une colonne `Vector{Int}` sous la forme
`"[0, 1, 2, 2]"` (avec crochets), et non comme une chaîne de chiffres brute.
Cette fonction gère les deux formats.
"""
function parse_residue_labels_field(s::AbstractString)
    s = strip(s)
    if startswith(s, "[")
        inner = strip(s, ['[', ']'])
        isempty(inner) && return Int[]
        return [parse(Int, strip(x)) for x in split(inner, ",")]
    else
        return parse_residue_labels(String(s))
    end
end
#------------------------------------------------------------------------

function join_disprot_pdb(df_disprot::DataFrame, df_pdb::DataFrame)
    if !("protein_id" in names(df_disprot)) || !("SP_PRIMARY" in names(df_pdb))
        error("Both dataframes must contain ID of uniprot column for joining.")
    end
    df_merged = outerjoin(df_disprot, df_pdb, on=:protein_id => :SP_PRIMARY, makeunique=true)
    return df_merged
end

function merge_by_chain(df_merged::DataFrame)
    autres_cols = setdiff(names(df_merged), ["PDB", "SP_BEG", "SP_END", "CHAIN"])

    df_merged_on_chain = combine(groupby(df_merged, [:PDB, :SP_BEG, :SP_END])) do sdf
        merge(
            (CHAIN = [passmissing(String).(sdf.CHAIN)],),
            NamedTuple(Symbol(c) => [first(sdf[!, c])] for c in autres_cols)
        )
    end

    return df_merged_on_chain
end

#----------------------------------------------------------------
# Diagnostic des résidus manquants via SIFTS (une protéine / un PDB)
#----------------------------------------------------------------

"""
    diagnose_labels_resolved(pdb, chain, adapt_label, sift_folder) -> (n_missing, n_missing_disorder, n_missing_disorder_binding, c_keep)

Version robuste : toute erreur (téléchargement, lecture du XML, comptage)
est capturée et journalisée avec `@warn`/`@error` ; la chaîne concernée est
simplement ignorée plutôt que de faire planter tout le job.
"""
function diagnose_labels_resolved(pdb, chain, adapt_label, sift_folder::String)
    siftsfile = joinpath(sift_folder, lowercase(pdb) * ".xml.gz")

    if !isfile(siftsfile)
         #=
        lock(DOWNLOAD_LOCK) do
            if !isfile(siftsfile)  # double vérification après acquisition du verrou
                @warn "Fichier SIFTS introuvable pour $pdb, tentative de téléchargement..."
                try
                    downloadsifts(uppercase(pdb), filename=siftsfile)
                catch e
                    @error "Échec du téléchargement SIFTS pour $pdb : $(sprint(showerror, e))"
                end
            end
        end
        =#
        if !isfile(siftsfile)
            @warn "Téléchargement impossible pour $pdb, N_MISSING* laissés manquants."
            return (missing, missing, missing, chain)
        end
    end

    best_n_missing = missing
    best_n_missing_disorder = missing
    best_n_missing_disorder_binding = missing
    c_keep = String[]

    for c in chain
        local residue_data
        try
            residue_data = MIToS.SIFTS.read_file(siftsfile, SIFTSXML; chain=String(c))
        catch e
            @warn "Lecture SIFTS impossible pour $pdb chaîne $c : $(sprint(showerror, e)). Chaîne ignorée."
            continue
        end

        n_missing = 0
        n_missing_disorder = 0
        n_missing_disorder_binding = 0
        i = 0

        try
            for res in residue_data
                i += 1
                i > length(adapt_label) && break
                label = adapt_label[i]
                if label == 1 && res.missing
                    n_missing_disorder += 1
                    n_missing += 1
                elseif label == 2 && res.missing
                    n_missing_disorder_binding += 1
                    n_missing += 1
                elseif res.missing
                    n_missing += 1
                end
            end
        catch e
            @warn "Erreur pendant le comptage des résidus manquants pour $pdb chaîne $c : $(sprint(showerror, e)). Chaîne ignorée."
            continue
        end

        if best_n_missing === missing || n_missing < best_n_missing
            best_n_missing = n_missing
            best_n_missing_disorder = n_missing_disorder
            best_n_missing_disorder_binding = n_missing_disorder_binding
            c_keep = [c]
        elseif n_missing == best_n_missing
            push!(c_keep, c)
        end
    end

    return (best_n_missing, best_n_missing_disorder, best_n_missing_disorder_binding, c_keep)
end

#----------------------------------------------------------------
# Traitement d'UNE protéine (un protein_id / SP_PRIMARY unique)
#----------------------------------------------------------------

"""
    adapt_labels_for_protein_group(subdf, sift_folder) -> DataFrame

Équivalent de l'ancienne boucle de `adapt_labels_to_pdb`, mais restreint aux
lignes d'UNE seule protéine (`subdf` = toutes les lignes de `df_merged_on_chain`
partageant le même `protein_id`). Chaque ligne est traitée dans son propre
try/catch : une ligne fautive est journalisée et ignorée, elle ne fait pas
échouer le reste de la protéine.
"""
function adapt_labels_for_protein_group(subdf::AbstractDataFrame, sift_folder::String)
    df_new_info = DataFrame(
        PDB = String[], CHAIN = Vector{String}[], ADAPTED_SEQUENCE = String[],
        ADAPTED_LABELS = Vector{Int}[], COMPLETE = Bool[],
        N_MISSING = Union{Int,Missing}[], RESIDUES_DISORDER_IN_ADAPTED_REGION = Bool[],
        LENGTH_DISORDER_IN_ADAPTED_REGION = Int[], N_MISSING_DISORDER = Union{Int,Missing}[],
        RESIDUES_BINDING_IN_ADAPTED_REGION = Bool[], LENGTH_BINDING_IN_ADAPTED_REGION = Int[],
        N_MISSING_DISORDER_BINDING = Union{Int,Missing}[],
    )

    for row in eachrow(subdf)
        try
            sequence = row.sequence
            residue_labels = row.residue_labels
            residue_labels_array = parse_residue_labels_field(residue_labels)

            if length(residue_labels_array) != length(sequence)
                @warn "protein_id=$(row.protein_id) : longueur residue_labels ≠ longueur sequence. Ligne ignorée."
                continue
            end
            if row.SP_BEG === missing || row.SP_END === missing
                @warn "protein_id=$(row.protein_id) : SP_BEG/SP_END manquant. Ligne ignorée."
                continue
            end

            res_beg = row.SP_BEG
            res_end = row.SP_END
            if res_beg < 1 || res_end > length(sequence) || res_beg > res_end
                @warn "protein_id=$(row.protein_id) : SP_BEG/SP_END invalide (res_beg=$res_beg, res_end=$res_end, len=$(length(sequence))). Ligne ignorée."
                continue
            end

            adapt_seq = sequence[res_beg:res_end]
            adapt_labels = residue_labels_array[res_beg:res_end]

            complete_resolution = length(sequence) == length(adapt_labels)
            disorder = any(==(1), adapt_labels)
            length_disorder = count(==(1), adapt_labels)
            binding = any(==(2), adapt_labels)
            length_binding = count(==(2), adapt_labels)

            (n_missing, n_missing_disorder, n_missing_disorder_binding, c_keep) =
                diagnose_labels_resolved(row.PDB, row.CHAIN, adapt_labels, sift_folder)

            push!(df_new_info, (row.PDB, c_keep, adapt_seq, adapt_labels, complete_resolution,
                                 n_missing, disorder, length_disorder, n_missing_disorder,
                                 binding, length_binding, n_missing_disorder_binding))
        catch e
            pid = hasproperty(row, :protein_id) ? row.protein_id : missing
            @warn "Ligne ignorée suite à une erreur inattendue (protein_id=$pid) : $(sprint(showerror, e))"
            continue
        end
    end

    return outerjoin(subdf, df_new_info, on=[:PDB, :CHAIN], makeunique=true)
end

"""
    safe_adapt_labels_for_protein(protein_id, subdf, sift_folder) -> (result, error_msg)

Encapsule `adapt_labels_for_protein_group` dans un try/catch "attrape-tout" :
si une erreur imprévue plante le traitement d'une protéine entière (bug,
donnée corrompue, etc.), elle est capturée ici, la protéine est marquée en
échec (et donc PAS ajoutée au journal de reprise), mais le reste du batch
et des batches suivants continuent normalement.
"""
function safe_adapt_labels_for_protein(protein_id, subdf::AbstractDataFrame, sift_folder::String)
    try
        result = adapt_labels_for_protein_group(subdf, sift_folder)
        return (result, nothing)
    catch e
        msg = "Erreur fatale pour protein_id=$(protein_id) :\n" *
              sprint(showerror, e) * "\n" *
              sprint(Base.show_backtrace, catch_backtrace())
        return (nothing, msg)
    end
end

#----------------------------------------------------------------
# Orchestrateur : parallélisation + checkpoint/reprise
#----------------------------------------------------------------

"""
    checkpointed_parallel_pipeline(df_merged_on_chain, sift_folder, output_csv;
                                    batch_size=15)

- Regroupe `df_merged_on_chain` par `protein_id` (= SP_PRIMARY unique).
- Traite chaque protéine indépendamment, en parallèle sur les threads
  disponibles (`Threads.@threads`), par lots de `batch_size` protéines
  (10 à 20 recommandé).
- Après chaque lot : écrit les résultats dans `output_csv` (append), note
  les protéines terminées dans `output_csv.done.txt`, et les erreurs dans
  `output_csv.errors.log`.
- Si le script est relancé après une interruption/erreur, les protéines déjà
  listées dans `output_csv.done.txt` sont automatiquement sautées : la
  reprise se fait donc au niveau du dernier lot sauvegardé.
- Une protéine en erreur n'est PAS marquée "done" : elle sera retentée au
  prochain lancement (elle n'est donc jamais perdue silencieusement).
"""
function checkpointed_parallel_pipeline(df_merged_on_chain::DataFrame,
                                         sift_folder::String,
                                         output_csv::String;
                                         batch_size::Int = 15)

    done_log  = output_csv * ".done.txt"
    error_log = output_csv * ".errors.log"

    if !("protein_id" in names(df_merged_on_chain))
        error("La colonne :protein_id est absente de df_merged_on_chain — impossible de grouper par protéine.")
    end

    gdf = groupby(df_merged_on_chain, :protein_id)
    group_keys = collect(keys(gdf))
    protein_ids = [k.protein_id for k in group_keys]

    valid_idx = findall(!ismissing, protein_ids)
    if length(valid_idx) < length(protein_ids)
        @warn "$(length(protein_ids) - length(valid_idx)) groupe(s) sans protein_id ignoré(s)."
    end

    # Reprise : lecture des protéines déjà traitées avec succès
    done_ids = Set{String}()
    if isfile(done_log)
        for line in eachline(done_log)
            line = strip(line)
            isempty(line) || push!(done_ids, line)
        end
        @info "Reprise détectée : $(length(done_ids)) protéine(s) déjà traitée(s), elles seront sautées."
    end

    master_indices = [i for i in valid_idx if !(string(protein_ids[i]) in done_ids)]
    @info "Protéines au total : $(length(protein_ids)) — restantes à traiter : $(length(master_indices))"

    if isempty(master_indices)
        @info "Rien à faire : toutes les protéines ont déjà été traitées."
        return
    end

    n_threads = nthreads()
    if n_threads == 1
        @warn "Julia tourne avec 1 seul thread : le traitement sera SÉRIE. " *
              "Relancez avec `julia --threads=\$SLURM_CPUS_PER_TASK` (ou définissez " *
              "JULIA_NUM_THREADS avant le lancement) pour paralléliser."
    else
        @info "Exécution avec $n_threads threads."
    end

    header_written = isfile(output_csv) && filesize(output_csv) > 0
    write_lock = ReentrantLock()

    for batch_start in 1:batch_size:length(master_indices)
        batch_end = min(batch_start + batch_size - 1, length(master_indices))
        batch_indices = master_indices[batch_start:batch_end]

        results = Vector{Union{DataFrame,Nothing}}(nothing, length(batch_indices))
        errors  = Vector{Union{Nothing,String}}(nothing, length(batch_indices))

        Threads.@threads for j in eachindex(batch_indices)
            i = batch_indices[j]
            pid = protein_ids[i]
            try
                subdf = DataFrame(gdf[i])  # copie locale : sûr en contexte multi-thread
                (res, err) = safe_adapt_labels_for_protein(pid, subdf, sift_folder)
                results[j] = res
                errors[j]  = err
            catch e
                # Filet de sécurité ultime : même une erreur ici (ex. copie du
                # sous-DataFrame) ne doit pas tuer le thread / le job entier.
                errors[j] = "Erreur imprévue avant traitement pour protein_id=$(pid) :\n" *
                            sprint(showerror, e)
            end
        end

        # Partie séquentielle (écriture fichiers) : protégée par verrou, même si
        # à ce stade on est de toute façon hors de la région @threads.
        lock(write_lock) do
            for j in eachindex(batch_indices)
                pid = protein_ids[batch_indices[j]]

                if errors[j] !== nothing
                    open(error_log, "a") do io
                        println(io, "="^80)
                        println(io, "protein_id = $(pid) — $(now_str())")
                        println(io, errors[j])
                    end
                    @error "Échec sur la protéine $(pid) — voir $(error_log). Elle sera retentée au prochain lancement."
                    continue
                end

                res = results[j]
                try
                    CSV.write(output_csv, res; append=header_written, writeheader=!header_written)
                    header_written = true
                catch e
                    open(error_log, "a") do io
                        println(io, "="^80)
                        println(io, "ÉCRITURE CSV échouée pour protein_id = $(pid) — $(now_str())")
                        println(io, sprint(showerror, e))
                    end
                    @error "Impossible d'écrire les résultats de $(pid) dans $(output_csv). Elle sera retentée."
                    continue
                end

                open(done_log, "a") do io
                    println(io, pid)
                end
            end
        end

        @info "Lot $(batch_start)-$(batch_end) / $(length(master_indices)) traité et sauvegardé ($(now_str()))."
    end

    @info "Pipeline terminé. Résultats : $(output_csv). Erreurs éventuelles : $(error_log)."
end

###################################################################

function main()
    sift_folder = "/alpha/database/sift_new-2025-10-27/sift_files"
    input_csv = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl/scripts/csv/cluster_dataset_labels_update.csv"

    df_disprot = try
        CSV.read(input_csv, DataFrame)
    catch e
        error("Impossible de lire le CSV d'entrée $input_csv : $(sprint(showerror, e))")
    end
    @info "df_disprot chargé" nrow(df_disprot)

    output_dir = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl/scripts/add_feature/output"
    isdir(output_dir) || mkpath(output_dir)

    pdb_chain_uniprot_csv = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl/scripts/add_feature/pdb_chain_uniprot.csv.gz"
    if !isfile(pdb_chain_uniprot_csv)
        try
            pdb_chain_uniprot_csv = abspath(downloadsifts(dbUniProt))
        catch e
            error("Fichier de correspondance PDB->UniProt introuvable et téléchargement impossible : $(sprint(showerror, e))")
        end
    end

    df_pdb = try
        CSV.read(pdb_chain_uniprot_csv, DataFrame; header=2)
    catch e
        error("Impossible de lire $pdb_chain_uniprot_csv : $(sprint(showerror, e))")
    end
    @info "df_pdb chargé" nrow(df_pdb)

    df_merged = join_disprot_pdb(df_disprot, df_pdb)
    @info "Fusion DisProt/PDB effectuée" nrow(df_merged)

    df_merged_on_chain = merge_by_chain(df_merged)
    @info "Fusion par chaîne effectuée" nrow(df_merged_on_chain)

    output_csv = joinpath(output_dir, "merged_disprot_pdb_with_adapted_labels.csv")

    # batch_size entre 10 et 20 comme demandé ; ajuster si besoin.
    checkpointed_parallel_pipeline(df_merged_on_chain, sift_folder, output_csv; batch_size=15)

    println("All downloads completed.")
end

try
    main()
catch e
    @error "Erreur fatale non rattrapée dans main()" exception=(e, catch_backtrace())
    exit(1)
end