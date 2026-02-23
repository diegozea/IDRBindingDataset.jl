using TestItems

@testitem "Build fixed labels table follows AF2, missing-date, and CAID override rules" begin
    using Dates
    using Test
    using IDRBindingDataset:
                             build_disprot_protein_annotations,
                             build_fixed_labels_table,
                             build_oldest_pdb_release_dates

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_B1", "P_C1", "P_D1"],
        global_representative_id = ["P_A1", "P_B1", "P_C1", "P_D1"],
        local_representative_id = ["P_A1", "P_B1", "P_C1", "P_D1"],
        final_cluster_id = ["cluster_A", "cluster_B", "cluster_C", "cluster_D"]
    )

    disprot_entries = Any[
    Dict(
        "acc" => "P_A1", "disprot_id" => "DP_A1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_B1", "disprot_id" => "DP_B1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_C1", "disprot_id" => "DP_C1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_D1", "disprot_id" => "DP_D1", "released" => "2020_01", "regions" => Any[])
]
    disprot_annotations = build_disprot_protein_annotations(disprot_entries)

    pdb_release_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P_A1", "P_B1", "P_C1"],
        pdb_release_date = ["2017-01-01", "missing", "2016-05-01"]
    )
    oldest_pdb = build_oldest_pdb_release_dates(pdb_release_map)

    fixed_labels = build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb,
        Set(["DP_C1"]), # CAID overrides to caid
        Set{String}()
    )

    @test names(fixed_labels) == [
        "protein_id",
        "global_representative_id",
        "local_representative_id",
        "final_cluster_id",
        "protein_disprot_id",
        "protein_has_idr_binding_annotation",
        "protein_oldest_pdb_release_date",
        "protein_has_pdb_structure",
        "protein_has_missing_pdb_release_date",
        "protein_in_caid2",
        "protein_in_caid3",
        "protein_in_af2_training_set",
        "cluster_in_af2_training_set",
        "cluster_in_caid2",
        "cluster_in_caid3",
        "cluster_fixed_dataset_label"
    ]

    by_acc = Dict(row.protein_id => row for row in eachrow(fixed_labels))

    @test by_acc["P_A1"].protein_in_af2_training_set
    @test by_acc["P_A1"].cluster_in_af2_training_set
    @test by_acc["P_A1"].cluster_fixed_dataset_label == "training"

    @test !by_acc["P_B1"].protein_in_af2_training_set
    @test by_acc["P_B1"].protein_has_pdb_structure
    @test by_acc["P_B1"].protein_has_missing_pdb_release_date
    @test by_acc["P_B1"].cluster_fixed_dataset_label == "training"

    @test by_acc["P_C1"].protein_in_af2_training_set
    @test by_acc["P_C1"].cluster_in_caid2
    @test by_acc["P_C1"].cluster_fixed_dataset_label == "caid"

    @test !by_acc["P_D1"].protein_in_af2_training_set
    @test !by_acc["P_D1"].protein_has_idr_binding_annotation
    @test !by_acc["P_D1"].protein_has_pdb_structure
    @test ismissing(by_acc["P_D1"].cluster_fixed_dataset_label)
end

@testitem "Build fixed labels table requires explicit PDB boolean columns" begin
    using Dates
    using Test
    using IDRBindingDataset:
                             build_disprot_protein_annotations,
                             build_fixed_labels_table

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_B1", "P_C1"],
        global_representative_id = ["P_A1", "P_B1", "P_C1"],
        local_representative_id = ["P_A1", "P_B1", "P_C1"],
        final_cluster_id = ["cluster_A", "cluster_B", "cluster_C"]
    )
    disprot_entries = Any[
    Dict(
        "acc" => "P_A1", "disprot_id" => "DP_A1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_B1", "disprot_id" => "DP_B1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_C1", "disprot_id" => "DP_C1", "released" => "2020_01", "regions" => Any[])
]
    disprot_annotations = build_disprot_protein_annotations(disprot_entries)

    # Legacy two-column format: missing explicit boolean columns.
    oldest_pdb_release_dates_legacy = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P_A1", "P_B1"],
        oldest_pdb_release_date = Any[missing, Date(2020, 1, 1)]
    )

    error = try
        build_fixed_labels_table(
            cluster_representatives,
            disprot_annotations,
            oldest_pdb_release_dates_legacy,
            Set{String}(),
            Set{String}()
        )
        nothing
    catch e
        e
    end
    @test error isa ArgumentError
    error_message = sprint(showerror, error)
    @test occursin("oldest_pdb_release_dates is missing required columns", error_message)
    @test occursin("has_pdb_structure", error_message)
    @test occursin("has_missing_pdb_release_date", error_message)
end

@testitem "Write fixed labels CSV.gz and README artifacts" begin
    using Test
    using IDRBindingDataset:
                             AF2_TRAINING_DATA_CUTOFF_DATE,
                             write_fixed_labels_csv_gz,
                             write_fixed_labels_readme

    tmp_dir = mktempdir()
    output_csv_path = joinpath(tmp_dir, "cluster_fixed_labels.csv.gz")
    output_readme_path = joinpath(tmp_dir, "README.md")

    fixed_labels = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P2", "P3"],
        final_cluster_id = ["cluster_1", "cluster_2", "cluster_3"],
        cluster_fixed_dataset_label = Any["training", "caid", missing]
    )

    write_fixed_labels_csv_gz(fixed_labels, output_csv_path)
    @test isfile(output_csv_path)

    reloaded = IDRBindingDataset.DataFrames.DataFrame(
        IDRBindingDataset.CSV.File(output_csv_path)
    )
    @test IDRBindingDataset.DataFrames.nrow(reloaded) == 3
    @test "cluster_fixed_dataset_label" in names(reloaded)

    write_fixed_labels_readme(
        output_readme_path;
        output_csv_path = "data/derived/fixed/cluster_fixed_labels.csv.gz",
        output_csv_filesystem_path = output_csv_path,
        source_cluster_representatives_path = "data/derived/clustering/cluster_representatives.csv.gz",
        source_disprot_regions_path = "data/raw/disprot/disprot_regions.json.gz",
        source_pdb_release_map_path = "data/derived/pdb/uniprot_pdb_release_date_map.csv.gz",
        source_caid2_binding_fasta_path = "data/raw/caid_2/binding.fasta.gz",
        source_caid3_binding_fasta_path = "data/raw/caid_3/binding.fasta.gz",
        af2_training_cutoff_date = AF2_TRAINING_DATA_CUTOFF_DATE,
        fixed_labels_table = fixed_labels
    )
    @test isfile(output_readme_path)

    readme_text = read(output_readme_path, String)
    @test occursin("Fixed Dataset Labels", readme_text)
    @test occursin("$(filesize(output_csv_path))", readme_text)
    @test occursin("label `training`: `1`", readme_text)
    @test occursin("label `caid`: `1`", readme_text)
    @test occursin("label `missing`: `1`", readme_text)
end
