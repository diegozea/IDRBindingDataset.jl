#!/usr/bin/env julia

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

using ArgParse
using CSV
using DataFrames
using CodecZlib
using MIToS
using MIToS.SIFTS

############################ FUNCTIONS ############################
# ----------------------------------------------------------------------
# Parsing / label utilities
# ----------------------------------------------------------------------

"""
    parse_residue_labels(s::String) -> Vector{Int}

Convert a string of digit characters (e.g. `"0102201"`) into a vector
of integer residue labels. This is the format found in the RAW input
dataset CSV.
"""
function parse_residue_labels(s::String)
    return [parse(Int, c) for c in s]
end

"""
    parse_residue_labels_field(s::AbstractString) -> Vector{Int}

Robust parser used when RE-READING a CSV that this pipeline already
wrote out: `CSV.write` serializes a `Vector{Int}` column as
`"[0, 1, 2, 2]"` (brackets + comma-separated), not as a bare digit
string like the original raw input. This function handles BOTH forms,
so it's safe to use on either a fresh raw CSV or one of our own
annotated/merged CSVs.
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
    # Join DisProt and PDB dataframes on protein_id
    if !("protein_id" in names(df_disprot)) || !("SP_PRIMARY" in names(df_pdb))
        error("Both dataframes must contain ID of uniprot column for joining.")
    end
    df_merged = outerjoin(df_disprot, df_pdb, on=:protein_id => :SP_PRIMARY, makeunique=true)

    return df_merged
end

function merge_by_chain(df_merged::DataFrame)
    # Group by PDB ID, res_beg, res_end and merge CHAIN values into a single row
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

"""
    diagnose_labels_resolved!(df::DataFrame, sift_folder::String) -> DataFrame

Add `N_MISSING_DISORDER::Union{Missing,Int}` and
`N_MISSING_DISORDER_BINDING::Union{Missing,Int}` columns using the
SIFTS XML file for each PDB entry. Binding status is now read PER ROW
from `protein_has_idr_binding_annotation` (instead of a single global
`binding` flag): `N_MISSING_DISORDER_BINDING` stays `missing` for
non-binding rows. If the SIFTS file is not found for a row, both
columns are left `missing` for that row only (no early return / no
rows dropped, unlike the old version).
"""
function diagnose_labels_resolved(pdb, chain,adapt_label, sift_folder::String)
    
    siftsfile = joinpath(sift_folder, lowercase(pdb)* ".xml.gz")

    if !isfile(siftsfile)
        @warn "SIFTS file not found for $pdb, attempting download"
        try 
            downloadsifts(uppercase(pdb), filename=siftsfile)
        catch e
            @error "Failed to download SIFTS file for $pdb: $e"
            return (missing, missing, missing,chain)
        end
        if !isfile(siftsfile)
            @warn "Download failed silently for $pdb, leaving N_MISSING_DISORDER* as missing"
            return (missing, missing, missing,chain)
        end
    end

    best_n_missing = missing
    best_n_missing_disorder = missing
    best_n_missing_disorder_binding = missing
    c_keep = String[]
    for c in chain
        residue_data = MIToS.SIFTS.read_file(siftsfile, SIFTSXML; chain=String(c))
            
        n_missing = 0
        n_missing_disorder = 0
        n_missing_disorder_binding = 0

        i = 0
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
        if best_n_missing === missing || n_missing < best_n_missing
            # nouvelle meilleure chaîne : on repart de zéro
            best_n_missing = n_missing
            best_n_missing_disorder = n_missing_disorder
            best_n_missing_disorder_binding = n_missing_disorder_binding
            c_keep = [c]
        elseif n_missing == best_n_missing
            # égalité avec le best actuel : on garde aussi celle-ci
            push!(c_keep, c)
        end
    end

    return (best_n_missing, best_n_missing_disorder, best_n_missing_disorder_binding,c_keep)
end


function adapt_labels_to_pdb(df_merged_on_chain::DataFrame, sift_folder::String)
    df_new_info = DataFrame(PDB=String[], CHAIN=Vector{String}[], ADAPTED_SEQUENCE=String[], ADAPTED_LABELS=Vector{Int}[],COMPLETE=Bool[],N_MISSING=Union{Int, Missing}[],RESIDUES_DISORDER_IN_ADAPTED_REGION=Bool[],LENGTH_DISORDER_IN_ADAPTED_REGION=Int[],N_MISSING_DISORDER=Union{Int, Missing}[],RESIDUES_BINDING_IN_ADAPTED_REGION=Bool[],LENGTH_BINDING_IN_ADAPTED_REGION=Int[],N_MISSING_DISORDER_BINDING=Union{Int, Missing}[])
    # Adapt label to PDB


    for row in eachrow(df_merged_on_chain)
        sequence = row.sequence
        residue_labels = row.residue_labels
        residue_labels_array = parse_residue_labels_field(residue_labels)

        # Ensure that the length of residue_labels_array matches the length of the sequence
        if length(residue_labels_array) != length(sequence)
            @warn "Length of residue_labels_array does not match length of sequence for protein_id $(row.protein_id). Skipping adaptation."
            continue
        end
        @show length(sequence), length(residue_labels_array)
        if row.SP_BEG === missing || row.SP_END === missing
            @warn "SP_BEG or SP_END is missing for protein_id $(row.protein_id). Skipping adaptation."
            continue
        end
        # Adapt labels to PDB
        res_beg = row.SP_BEG
        res_end = row.SP_END
        if res_beg < 1 || res_end > length(sequence) || res_beg > res_end
            @warn "Invalid SP_BEG or SP_END for protein_id $(row.protein_id). Skipping adaptation."
            continue
        end
        adapt_seq=sequence[res_beg:res_end]
        adapt_labels=residue_labels_array[res_beg:res_end]
        @show length(adapt_seq), length(adapt_labels)

        #Check if complete resolution is available
        if length(sequence) == length(adapt_labels)
            complete_resolution=true
        else
            complete_resolution=false
        end
        # Check disorder
        if any(x -> x == 1, adapt_labels)
            disorder=true
            length_disorder=sum(x -> x == 1, adapt_labels)
        else
            disorder=false
            length_disorder=0
        end

        # Check binding
        if any(x -> x == 2, adapt_labels)
            binding=true
            length_binding=sum(x -> x == 2, adapt_labels)
        else
            binding=false
            length_binding=0
        end

        # Check missing residues in the adapted region
        (n_missing, n_missing_disorder, n_missing_disorder_binding,c_keep) = diagnose_labels_resolved(row.PDB, row.CHAIN,adapt_labels, sift_folder)
        
        push!(df_new_info, (row.PDB, c_keep, adapt_seq, adapt_labels, complete_resolution, n_missing, disorder, length_disorder, n_missing_disorder, binding, length_binding, n_missing_disorder_binding))
    end

    df_complete_info_pdb = outerjoin(df_merged_on_chain, df_new_info, on=[:PDB,:CHAIN], makeunique=true)

    return df_complete_info_pdb
end
###################################################################

"""
    parse_command_line(args)

Parse command-line arguments and return selected versions.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Download DisProt regions, IDPO ontology, and IDPO mappings into data/raw/disprot.

Run from the repository root:
  julia scripts/01_download_disprot_data.jl --input_csv /path/to/input.csv --output_dir /path/to/output_dir --pdb_chain_uniprot_csv /path/to/pdb_chain_uniprot.csv

The script activates and instantiates scripts/Project.toml internally.
""",
    )

    @add_arg_table! settings begin
        "--release"
        help = "DisProt release identifier (default: $(RELEASE))."
        arg_type = String
        default = RELEASE
    end

    parsed = parse_args(args, settings)
    return (; release = parsed["release"])
end

function main()
    #options = parse_command_line(ARGS)
    sift_folder = "/alpha/database/sift_new-2025-10-27/sift_files"
    input_csv = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl/scripts/csv/cluster_dataset_labels_update.csv"
    df_disprot = CSV.read(input_csv, DataFrame)
    @show nrow(df_disprot)
    @show first(df_disprot, 1)

    output_dir = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl/scripts/add_feature/output"
    isdir(output_dir) || mkdir(output_dir)

    pdb_chain_uniprot_csv = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl/scripts/add_feature/pdb_chain_uniprot.csv.gz"
    df_pdb = DataFrame()
    if !isfile(pdb_chain_uniprot_csv)
        try 
            pdb_chain_uniprot_csv = abspath(downloadsifts(dbUniProt))
        catch e
            error("PDB chain to UniProt mapping CSV file not found: $pdb_chain_uniprot_csv")
        end
    end
    df_pdb = CSV.read(pdb_chain_uniprot_csv, DataFrame; header=2)
    @show nrow(df_pdb)
    @show first(df_pdb, 1)

    # Join DisProt and PDB dataframes
    df_merged = join_disprot_pdb(df_disprot, df_pdb)
    @show nrow(df_merged)
    @show first(df_merged, 1)

    df_merged_on_chain = merge_by_chain(df_merged)
    @show nrow(df_merged_on_chain)
    @show first(df_merged_on_chain, 5)

    # Adapt label to pdb 
    df_complete_info_pdb = adapt_labels_to_pdb(df_merged_on_chain, sift_folder)
    @show nrow(df_complete_info_pdb)
    @show first(df_complete_info_pdb, 5)

    CSV.write(joinpath(output_dir, "merged_disprot_pdb_with_adapted_labels.csv"), df_complete_info_pdb)

    println("All downloads completed.")
end

main()
