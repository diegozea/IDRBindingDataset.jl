using TestItems

@testitem "SIFTS databases include only UniProt and CATH summaries" begin
    using Test

    @test IDRBindingDataset.SIFTS_DATABASES == (
        IDRBindingDataset.MIToS.SIFTS.dbUniProt,
        IDRBindingDataset.MIToS.SIFTS.dbCATH
    )
end

@testitem "SIFTS mapping description helper" begin
    using Test
    using IDRBindingDataset: _sifts_mapping_description

    @test _sifts_mapping_description(IDRBindingDataset.MIToS.SIFTS.dbUniProt) ==
          "PDB chain to UniProt mapping"
    @test _sifts_mapping_description(IDRBindingDataset.MIToS.SIFTS.dbCATH) ==
          "PDB chain to CATH mapping with UniProt cross-reference"
end

@testitem "Write SIFTS folder README includes sources and metadata" begin
    using Test
    using IDRBindingDataset: _write_sifts_folder_readme

    tmp_dir = mktempdir()
    uniprot_path = joinpath(tmp_dir, "pdb_chain_uniprot.csv.gz")
    cath_path = joinpath(tmp_dir, "pdb_chain_cath_uniprot.csv.gz")
    write(uniprot_path, "uni")
    write(cath_path, "cath")

    entries = [
        (
            database = "dbUniProt",
            filename = "pdb_chain_uniprot.csv.gz",
            mapping = "PDB chain to UniProt mapping",
            destination = uniprot_path
        ),
        (
            database = "dbCATH",
            filename = "pdb_chain_cath_uniprot.csv.gz",
            mapping = "PDB chain to CATH mapping with UniProt cross-reference",
            destination = cath_path
        )
    ]

    timestamp = "2026-02-17T13:40:00Z"
    _write_sifts_folder_readme(tmp_dir, timestamp, entries)

    readme = read(joinpath(tmp_dir, "README.md"), String)
    @test occursin("# SIFTS Chain-Level Mapping Data", readme)
    @test occursin("Download timestamp (UTC): `$(timestamp)`", readme)
    @test occursin("Downloaded with: `MIToS.SIFTS.downloadsifts`", readme)
    @test occursin("Source documentation:", readme)
    @test occursin("pdb_chain_uniprot.csv.gz", readme)
    @test occursin("pdb_chain_cath_uniprot.csv.gz", readme)
    @test occursin("PDB chain to UniProt mapping", readme)
    @test occursin("PDB chain to CATH mapping with UniProt cross-reference", readme)
end
