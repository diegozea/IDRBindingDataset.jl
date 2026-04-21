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
    disprot_annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = IDRBindingDataset._DisProtChallengeTermRules(
            Set(["GO:0005488", "GO:0005515", "GO:0003676"]),
            Set(["GO:0003676"]),
            Set(["GO:0005515"]),
            Set([
                "IDPO:0000002",
                "IDPO:0000003",
                "IDPO:0000004",
                "IDPO:0000011",
                "IDPO:0000017",
                "IDPO:0000014"
            ])
        )
    )

    pdb_release_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P_A1", "P_B1", "P_C1"],
        pdb_release_date = ["2017-01-01", "missing", "2016-05-01"]
    )
    oldest_pdb = build_oldest_pdb_release_dates(pdb_release_map)

    fixed_labels = build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb,
        Dict(
            "caid_2" => Set(["DP_C1"]), # CAID overrides to caid
            "caid_3" => Set{String}()
        )
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

@testitem "Build fixed labels table includes only provided CAID versions" begin
    using Dates
    using Test
    using IDRBindingDataset:
                             build_disprot_protein_annotations,
                             build_fixed_labels_table

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_B1"],
        global_representative_id = ["P_A1", "P_B1"],
        local_representative_id = ["P_A1", "P_B1"],
        final_cluster_id = ["cluster_A", "cluster_B"]
    )
    disprot_entries = Any[
    Dict(
        "acc" => "P_A1", "disprot_id" => "DP_A1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_B1", "disprot_id" => "DP_B1", "released" => "2020_01", "regions" => Any[])
]
    disprot_annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = IDRBindingDataset._DisProtChallengeTermRules(
            Set(["GO:0005488", "GO:0005515", "GO:0003676"]),
            Set(["GO:0003676"]),
            Set(["GO:0005515"]),
            Set([
                "IDPO:0000002",
                "IDPO:0000003",
                "IDPO:0000004",
                "IDPO:0000011",
                "IDPO:0000017",
                "IDPO:0000014"
            ])
        )
    )

    oldest_pdb_release_dates = IDRBindingDataset.DataFrames.DataFrame(
        acc = String[],
        oldest_pdb_release_date = Union{Missing, Date}[],
        has_pdb_structure = Bool[],
        has_missing_pdb_release_date = Bool[]
    )

    caid_disprot_ids_by_dataset = Dict(
        "caid_4" => Set(["DP_B1"])
    )

    fixed_labels = build_fixed_labels_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb_release_dates,
        caid_disprot_ids_by_dataset
    )
    by_acc = Dict(row.protein_id => row for row in eachrow(fixed_labels))

    @test ismissing(by_acc["P_A1"].cluster_fixed_dataset_label)

    @test by_acc["P_B1"].protein_in_caid4
    @test by_acc["P_B1"].cluster_in_caid4
    @test by_acc["P_B1"].cluster_fixed_dataset_label == "caid"
    @test "protein_in_caid4" in names(fixed_labels)
    @test "cluster_in_caid4" in names(fixed_labels)
    @test "protein_in_caid2" ∉ names(fixed_labels)
    @test "protein_in_caid3" ∉ names(fixed_labels)
    @test "cluster_in_caid2" ∉ names(fixed_labels)
    @test "cluster_in_caid3" ∉ names(fixed_labels)
end

@testitem "Build fixed labels table rejects duplicate CAID version keys" begin
    using Dates
    using Test
    using IDRBindingDataset:
                             build_disprot_protein_annotations,
                             build_fixed_labels_table

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_B1"],
        global_representative_id = ["P_A1", "P_B1"],
        local_representative_id = ["P_A1", "P_B1"],
        final_cluster_id = ["cluster_A", "cluster_B"]
    )
    disprot_entries = Any[
    Dict(
        "acc" => "P_A1", "disprot_id" => "DP_A1", "released" => "2020_01", "regions" => Any[]),
    Dict(
        "acc" => "P_B1", "disprot_id" => "DP_B1", "released" => "2020_01", "regions" => Any[])
]
    disprot_annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = IDRBindingDataset._DisProtChallengeTermRules(
            Set(["GO:0005488", "GO:0005515", "GO:0003676"]),
            Set(["GO:0003676"]),
            Set(["GO:0005515"]),
            Set([
                "IDPO:0000002",
                "IDPO:0000003",
                "IDPO:0000004",
                "IDPO:0000011",
                "IDPO:0000017",
                "IDPO:0000014"
            ])
        )
    )
    oldest_pdb_release_dates = IDRBindingDataset.DataFrames.DataFrame(
        acc = String[],
        oldest_pdb_release_date = Union{Missing, Date}[],
        has_pdb_structure = Bool[],
        has_missing_pdb_release_date = Bool[]
    )
    caid_disprot_ids_by_dataset = Dict(
        "caid_3" => Set(["DP_A1"]),
        "caid_03" => Set(["DP_B1"])
    )

    error = try
        build_fixed_labels_table(
            cluster_representatives,
            disprot_annotations,
            oldest_pdb_release_dates,
            caid_disprot_ids_by_dataset
        )
        nothing
    catch e
        e
    end
    @test error isa ArgumentError
    @test occursin("Duplicate CAID dataset version 3", sprint(showerror, error))
end

@testitem "Discover CAID binding datasets and require specific versions" begin
    using Test
    using IDRBindingDataset:
                             _discover_caid_binding_datasets,
                             _require_caid_dataset_path

    raw_data_dir = mktempdir()
    mkpath(joinpath(raw_data_dir, "caid_2"))
    write(joinpath(raw_data_dir, "caid_2", "binding.fasta.gz"), "")
    mkpath(joinpath(raw_data_dir, "caid_4"))
    write(joinpath(raw_data_dir, "caid_4", "binding.fasta.gz"), "")
    mkpath(joinpath(raw_data_dir, "not_a_caid_folder"))

    datasets = _discover_caid_binding_datasets(raw_data_dir)
    @test [dataset.name for dataset in datasets] == ["caid_2", "caid_4"]
    @test [dataset.number for dataset in datasets] == [2, 4]

    caid2_path = _require_caid_dataset_path(datasets, 2)
    @test caid2_path == normpath(joinpath(raw_data_dir, "caid_2", "binding.fasta.gz"))

    error = try
        _require_caid_dataset_path(datasets, 3)
        nothing
    catch e
        e
    end
    @test error isa ArgumentError
    @test occursin(
        "Required CAID dataset folder not found: caid_3",
        sprint(showerror, error)
    )
end

@testitem "Discover CAID binding datasets rejects duplicate normalized versions" begin
    using Test
    using IDRBindingDataset: _discover_caid_binding_datasets

    raw_data_dir = mktempdir()
    for folder_name in ("caid_3", "caid_03")
        mkpath(joinpath(raw_data_dir, folder_name))
        write(joinpath(raw_data_dir, folder_name, "binding.fasta.gz"), "")
    end

    error = try
        _discover_caid_binding_datasets(raw_data_dir)
        nothing
    catch e
        e
    end
    @test error isa ArgumentError
    error_message = sprint(showerror, error)
    @test occursin("Duplicate CAID dataset version numbers found", error_message)
    @test occursin("caid_03", error_message)
    @test occursin("caid_3", error_message)
end

@testitem "Parse date argument enforces yyyy-mm-dd format" begin
    using Dates
    using Test
    using IDRBindingDataset: _parse_date_argument

    @test _parse_date_argument("2018-04-30", "--af2_training_cutoff_date") ==
          Date(2018, 4, 30)

    error = try
        _parse_date_argument("30-04-2018", "--af2_training_cutoff_date")
        nothing
    catch e
        e
    end
    @test error isa ArgumentError
    @test occursin("Invalid --af2_training_cutoff_date", sprint(showerror, error))
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
    disprot_annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = IDRBindingDataset._DisProtChallengeTermRules(
            Set(["GO:0005488", "GO:0005515", "GO:0003676"]),
            Set(["GO:0003676"]),
            Set(["GO:0005515"]),
            Set([
                "IDPO:0000002",
                "IDPO:0000003",
                "IDPO:0000004",
                "IDPO:0000011",
                "IDPO:0000017",
                "IDPO:0000014"
            ])
        )
    )

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
            Dict{String, Set{String}}()
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
        source_caid_binding_fasta_paths_by_version = Dict(
            2 => "data/raw/caid_2/binding.fasta.gz",
            3 => "data/raw/caid_3/binding.fasta.gz",
            4 => "data/raw/caid_4/binding.fasta.gz"
        ),
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
    @test occursin("source_caid2_binding_fasta_path", readme_text)
    @test occursin("source_caid3_binding_fasta_path", readme_text)
    @test occursin("source_caid4_binding_fasta_path", readme_text)
end
