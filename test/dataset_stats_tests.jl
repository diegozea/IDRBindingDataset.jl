using TestItems

@testitem "Read CAID DisProt IDs from FASTA headers" begin
    using Test
    using IDRBindingDataset: read_caid_disprot_ids

    tmp_dir = mktempdir()
    input_path = joinpath(tmp_dir, "binding.fasta.gz")
    open(input_path, "w") do io
        gzip_io = IDRBindingDataset.CodecZlib.GzipCompressorStream(io)
        try
            write(gzip_io, ">DP00001 some comment\n")
            write(gzip_io, "SEQUENCE\n")
            write(gzip_io, "01010101\n")
            write(gzip_io, ">DP00002\n")
            write(gzip_io, "SEQUENCE\n")
            write(gzip_io, "11111111\n")
            write(gzip_io, ">DP00001 duplicate\n")
            write(gzip_io, "SEQUENCE\n")
            write(gzip_io, "00000000\n")
        finally
            close(gzip_io)
        end
    end

    ids = read_caid_disprot_ids(input_path)
    @test ids == Set(["DP00001", "DP00002"])
end

@testitem "Build oldest PDB release date table parses strings and missing values" begin
    using Dates
    using Test
    using IDRBindingDataset: build_oldest_pdb_release_dates

    release_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P1", "P1", "P2", "P2", "P3"],
        pdb_release_date = Any["2020-01-01", Date(2019, 5, 3), "missing", "", missing]
    )

    oldest = build_oldest_pdb_release_dates(release_map)
    by_acc = Dict(row.acc => row for row in eachrow(oldest))

    @test by_acc["P1"].oldest_pdb_release_date == Date(2019, 5, 3)
    @test by_acc["P2"].oldest_pdb_release_date === missing
    @test by_acc["P3"].oldest_pdb_release_date === missing
    @test by_acc["P1"].has_pdb_structure
    @test by_acc["P2"].has_pdb_structure
    @test by_acc["P3"].has_pdb_structure
    @test !by_acc["P1"].has_missing_pdb_release_date
    @test by_acc["P2"].has_missing_pdb_release_date
    @test by_acc["P3"].has_missing_pdb_release_date
end

@testitem "Build oldest PDB release date table rejects invalid date strings" begin
    using Test
    using IDRBindingDataset: build_oldest_pdb_release_dates

    release_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P1"],
        pdb_release_date = ["invalid_date"]
    )

    @test_throws ArgumentError build_oldest_pdb_release_dates(release_map)
end

@testitem "Build DisProt protein annotations captures oldest release and binding labels" begin
    using Dates
    using Test
    using IDRBindingDataset: build_disprot_protein_annotations

    disprot_entries = Any[
        Dict(
            "acc" => "P11111",
            "disprot_id" => "DP00001",
            "released" => "2020_04",
            "regions" => Any[
                Dict(
                    "released" => "2021_01",
                    "term_id" => "IDPO:0000002",
                    "term_is_binding" => false
                ),
                Dict(
                    "released" => "2020_09",
                    "term_id" => "GO:0005515",
                    "term_is_binding" => true
                )
            ]
        ),
        Dict(
            "acc" => "P22222",
            "disprot_id" => "DP00002",
            "released" => "2021_03",
            "regions" => Any[
                Dict(
                "released" => "2021_07",
                "term_id" => "GO:0003676"
            )
            ]
        ),
        Dict(
            "acc" => "P33333",
            "disprot_id" => "DP00003",
            "released" => "2022_02",
            "regions" => Any[
                Dict("released" => "2022_08", "term_id" => "IDPO:0000002")
            ]
        )
    ]

    annotations = build_disprot_protein_annotations(
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
    @test annotations.acc == ["P11111", "P22222", "P33333"]

    by_acc = Dict(row.acc => row for row in eachrow(annotations))
    @test by_acc["P11111"].oldest_disprot_release_date == Date(2020, 4, 1)
    @test by_acc["P22222"].oldest_disprot_release_date == Date(2021, 3, 1)
    @test by_acc["P33333"].oldest_disprot_release_date == Date(2022, 2, 1)
    @test by_acc["P11111"].has_idr_binding_annotation
    @test by_acc["P22222"].has_idr_binding_annotation
    @test !by_acc["P33333"].has_idr_binding_annotation
end

@testitem "Build DisProt protein annotations uses ontology rules when provided" begin
    using Test
    using IDRBindingDataset:
                             _DisProtChallengeTermRules,
                             build_disprot_protein_annotations

    term_rules = _DisProtChallengeTermRules(
        Set(["GO:0005488", "GO:0005515"]),
        Set(["GO:0003676"]),
        Set(["GO:0005515"]),
        Set(["IDPO:0000002"])
    )

    disprot_entries = Any[
        Dict(
            "acc" => "P11111",
            "disprot_id" => "DP00001",
            "released" => "2020_04",
            "regions" => Any[
                Dict(
                "term_id" => "GO:0060090",
                "term_name" => "molecular adaptor activity",
                "term_is_binding" => true
            )
            ]
        ),
        Dict(
            "acc" => "P22222",
            "disprot_id" => "DP00002",
            "released" => "2020_07",
            "regions" => Any[
                Dict("term_id" => "GO:0005515", "term_name" => "protein binding")
            ]
        )
    ]

    annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = term_rules
    )
    by_acc = Dict(row.acc => row for row in eachrow(annotations))

    @test !by_acc["P11111"].has_idr_binding_annotation
    @test by_acc["P22222"].has_idr_binding_annotation
end

@testitem "Build DisProt protein annotations can filter invalid region IDs" begin
    using Test
    using IDRBindingDataset:
                             _DisProtChallengeTermRules,
                             build_disprot_protein_annotations

    term_rules = _DisProtChallengeTermRules(
        Set(["GO:0005488", "GO:0005515"]),
        Set(["GO:0003676"]),
        Set(["GO:0005515"]),
        Set(["IDPO:0000002"])
    )

    disprot_entries = Any[
        Dict(
        "acc" => "P11111",
        "disprot_id" => "DP00001",
        "released" => "2020_04",
        "regions" => Any[
            Dict("region_id" => "DP00001r001", "term_id" => "GO:0005515")
        ]
    )
    ]

    annotations = build_disprot_protein_annotations(
        disprot_entries;
        term_rules = term_rules,
        invalid_region_ids = Set(["DP00001r001"])
    )
    @test !only(annotations.has_idr_binding_annotation)
end

@testitem "Build DisProt release dates extracts unique sorted dates" begin
    using Dates
    using Test
    using IDRBindingDataset: build_disprot_release_dates

    disprot_entries = Any[
        Dict(
            "acc" => "P11111",
            "disprot_id" => "DP00001",
            "released" => "2021_03",
            "regions" => Any[
                Dict("released" => "2020_06"),
                Dict("released" => missing)
            ]
        ),
        Dict(
            "acc" => "P22222",
            "disprot_id" => "DP00002",
            "released" => "2020_06",
            "regions" => Any[
                Dict("released" => "2022_01"),
                Dict()
            ]
        )
    ]

    @test build_disprot_release_dates(disprot_entries) == [
        Date(2020, 6, 1),
        Date(2021, 3, 1),
        Date(2022, 1, 1)
    ]
end

@testitem "Build DisProt release dates rejects invalid release tag" begin
    using Test
    using IDRBindingDataset: build_disprot_release_dates

    disprot_entries = Any[
        Dict(
        "acc" => "P11111",
        "disprot_id" => "DP00001",
        "released" => "invalid",
        "regions" => Any[]
    )
    ]

    @test_throws ArgumentError build_disprot_release_dates(disprot_entries)
end

@testitem "Build dataset stats table joins clustering labels and metadata" begin
    using Dates
    using Test
    using IDRBindingDataset:
                             build_dataset_stats_table,
                             build_disprot_protein_annotations,
                             build_oldest_pdb_release_dates

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P11111", "P22222", "P33333"],
        final_cluster_id = ["cluster_000001", "cluster_000001", "cluster_000002"]
    )

    disprot_entries = Any[
        Dict(
            "acc" => "P11111",
            "disprot_id" => "DP00001",
            "released" => "2020_04",
            "regions" => Any[Dict("released" => "2020_04", "term_id" => "IDPO:0000002")]
        ),
        Dict(
            "acc" => "P22222",
            "disprot_id" => "DP00002",
            "released" => "2021_03",
            "regions" => Any[Dict("released" => "2021_06", "term_id" => "GO:0005515")]
        ),
        Dict(
            "acc" => "P33333",
            "disprot_id" => "DP00003",
            "released" => "2022_02",
            "regions" => Any[Dict("released" => "2022_02", "term_id" => "IDPO:0000002")]
        )
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
        acc = ["P11111", "P11111", "P33333"],
        pdb_release_date = ["2019-01-01", "2020-02-02", "missing"]
    )
    oldest_pdb = build_oldest_pdb_release_dates(pdb_release_map)

    stats = build_dataset_stats_table(
        cluster_representatives,
        disprot_annotations,
        oldest_pdb,
        Dict(
            "caid_2" => Set(["DP00001", "DP00003"]),
            "caid_3" => Set(["DP00002"])
        )
    )

    @test names(stats) == [
        "protein_id",
        "final_cluster_id",
        "protein_disprot_id",
        "protein_oldest_disprot_release_date",
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

    by_acc = Dict(row.protein_id => row for row in eachrow(stats))
    @test by_acc["P11111"].protein_oldest_pdb_release_date == Date(2019, 1, 1)
    @test by_acc["P22222"].protein_oldest_pdb_release_date === missing
    @test by_acc["P11111"].protein_has_pdb_structure
    @test !by_acc["P22222"].protein_has_pdb_structure
    @test by_acc["P33333"].protein_has_missing_pdb_release_date
    @test by_acc["P11111"].protein_in_caid2
    @test !by_acc["P11111"].protein_in_caid3
    @test !by_acc["P22222"].protein_in_caid2
    @test by_acc["P22222"].protein_in_caid3
    @test !by_acc["P11111"].protein_in_af2_training_set
    @test !by_acc["P11111"].cluster_in_af2_training_set
    @test by_acc["P11111"].cluster_in_caid2
    @test by_acc["P11111"].cluster_in_caid3
    @test by_acc["P11111"].cluster_fixed_dataset_label == "caid"
    @test by_acc["P33333"].cluster_fixed_dataset_label == "caid"
    @test !by_acc["P33333"].protein_has_idr_binding_annotation

    legacy_clusters = copy(cluster_representatives)
    legacy_clusters.global_representative_id = ["P11111", "P11111", "P33333"]
    legacy_clusters.local_representative_id = ["P11111", "P22222", "P33333"]
    legacy_clusters.source_note = ["alpha", "beta", "gamma"]
    original_legacy_clusters = copy(legacy_clusters)
    legacy_stats = build_dataset_stats_table(
        legacy_clusters,
        disprot_annotations,
        oldest_pdb,
        Dict(
            "caid_2" => Set(["DP00001", "DP00003"]),
            "caid_3" => Set(["DP00002"])
        )
    )

    @test "global_representative_id" ∉ names(legacy_stats)
    @test "local_representative_id" ∉ names(legacy_stats)
    @test isequal(legacy_stats[:, names(stats)], stats)
    @test legacy_stats.source_note == ["alpha", "beta", "gamma"]
    @test isequal(legacy_clusters, original_legacy_clusters)
end

@testitem "Build dataset stats table convenience overload derives helper tables" begin
    using Dates
    using Test
    using IDRBindingDataset: build_dataset_stats_table

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P11111", "P22222"],
        final_cluster_id = ["cluster_000001", "cluster_000002"]
    )
    disprot_entries = Any[
        Dict(
            "acc" => "P11111",
            "disprot_id" => "DP00001",
            "released" => "2020_04",
            "regions" => Any[Dict("released" => "2020_04", "term_id" => "IDPO:0000002")]
        ),
        Dict(
            "acc" => "P22222",
            "disprot_id" => "DP00002",
            "released" => "2020_07",
            "regions" => Any[Dict("released" => "2021_01", "term_id" => "GO:0005515")]
        )
    ]
    pdb_release_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P11111", "P22222", "P22222"],
        pdb_release_date = ["2018-01-01", "2022-01-01", "2020-01-01"]
    )

    stats = build_dataset_stats_table(
        cluster_representatives,
        disprot_entries,
        pdb_release_map,
        Dict(
            "caid_2" => Set(["DP00001"]),
            "caid_3" => Set(["DP00002"])
        );
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

    by_acc = Dict(row.protein_id => row for row in eachrow(stats))
    @test by_acc["P11111"].protein_oldest_pdb_release_date == Date(2018, 1, 1)
    @test by_acc["P22222"].protein_oldest_pdb_release_date == Date(2020, 1, 1)
    @test by_acc["P11111"].protein_in_caid2
    @test by_acc["P22222"].protein_in_caid3
    @test by_acc["P11111"].cluster_fixed_dataset_label == "caid"
    @test by_acc["P22222"].cluster_fixed_dataset_label == "caid"
end

@testitem "Build dataset stats table keeps CAID columns fully dynamic" begin
    using Dates
    using Test
    using IDRBindingDataset: build_dataset_stats_table

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P11111", "P22222"],
        final_cluster_id = ["cluster_000001", "cluster_000002"]
    )
    disprot_entries = Any[
        Dict(
            "acc" => "P11111",
            "disprot_id" => "DP00001",
            "released" => "2020_04",
            "regions" => Any[Dict("released" => "2020_04", "term_id" => "IDPO:0000002")]
        ),
        Dict(
            "acc" => "P22222",
            "disprot_id" => "DP00002",
            "released" => "2020_07",
            "regions" => Any[Dict("released" => "2021_01", "term_id" => "GO:0005515")]
        )
    ]
    pdb_release_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = String[],
        pdb_release_date = String[]
    )

    stats = build_dataset_stats_table(
        cluster_representatives,
        disprot_entries,
        pdb_release_map,
        Dict("caid_4" => Set(["DP00002"]));
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

    @test "protein_in_caid4" in names(stats)
    @test "cluster_in_caid4" in names(stats)
    @test "protein_in_caid2" ∉ names(stats)
    @test "protein_in_caid3" ∉ names(stats)
    @test "cluster_in_caid2" ∉ names(stats)
    @test "cluster_in_caid3" ∉ names(stats)

    by_acc = Dict(row.protein_id => row for row in eachrow(stats))
    @test ismissing(by_acc["P11111"].cluster_fixed_dataset_label)
    @test by_acc["P22222"].cluster_fixed_dataset_label == "caid"
end

@testitem "Assign missing cluster labels by DisProt cutoff" begin
    using Dates
    using Test
    using IDRBindingDataset: assign_missing_cluster_labels_by_disprot_cutoff

    dataset_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_A2", "P_B1", "P_C1", "P_C2", "P_D1"],
        final_cluster_id = [
            "cluster_A", "cluster_A", "cluster_B", "cluster_C", "cluster_C", "cluster_D"],
        cluster_fixed_dataset_label = Any[
            "training", "training", "caid", missing, missing, missing],
        protein_oldest_disprot_release_date = Any[
            Date(2019, 1, 1),
            Date(2020, 1, 1),
            Date(2021, 6, 1),
            Date(2023, 2, 1),
            Date(2022, 6, 1),
            missing
        ]
    )

    labeled = assign_missing_cluster_labels_by_disprot_cutoff(
        dataset_stats;
        disprot_testing_cutoff_date = Date(2022, 1, 1)
    )
    by_acc = Dict(row.protein_id => row for row in eachrow(labeled))

    @test by_acc["P_A1"].cluster_disprot_cutoff_dataset_label == "training"
    @test by_acc["P_B1"].cluster_disprot_cutoff_dataset_label == "caid"
    @test by_acc["P_C1"].cluster_disprot_cutoff_dataset_label == "testing"
    @test by_acc["P_C2"].cluster_disprot_cutoff_dataset_label == "testing"
    @test by_acc["P_D1"].cluster_disprot_cutoff_dataset_label == "training"
    @test by_acc["P_C1"].cluster_min_protein_oldest_disprot_release_date == Date(2022, 6, 1)
    @test by_acc["P_D1"].cluster_min_protein_oldest_disprot_release_date === missing

    legacy_stats = copy(dataset_stats)
    legacy_stats.global_representative_id = copy(dataset_stats.protein_id)
    legacy_stats.local_representative_id = copy(dataset_stats.protein_id)
    original_legacy_stats = copy(legacy_stats)
    legacy_labeled = assign_missing_cluster_labels_by_disprot_cutoff(
        legacy_stats;
        disprot_testing_cutoff_date = Date(2022, 1, 1)
    )
    @test isequal(legacy_labeled, labeled)
    @test isequal(legacy_stats, original_legacy_stats)
end

@testitem "Assign missing cluster labels rejects inconsistent fixed labels inside cluster" begin
    using Dates
    using Test
    using IDRBindingDataset: assign_missing_cluster_labels_by_disprot_cutoff

    dataset_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_A2"],
        final_cluster_id = ["cluster_A", "cluster_A"],
        cluster_fixed_dataset_label = Any["training", "caid"],
        protein_oldest_disprot_release_date = Any[Date(2020, 1, 1), Date(2020, 2, 1)]
    )

    @test_throws ArgumentError assign_missing_cluster_labels_by_disprot_cutoff(
        dataset_stats;
        disprot_testing_cutoff_date = Date(2022, 1, 1)
    )
end

@testitem "Summarize dataset labels by binding counts proteins and clusters" begin
    using Test
    using IDRBindingDataset: summarize_dataset_labels_by_binding

    labeled_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_A2", "P_B1", "P_B2", "P_C1", "P_D1"],
        final_cluster_id = [
            "cluster_A", "cluster_A", "cluster_B", "cluster_B", "cluster_C", "cluster_D"],
        cluster_disprot_cutoff_dataset_label = [
            "training", "training", "caid", "caid", "testing", "testing"],
        protein_has_idr_binding_annotation = [true, false, false, true, false, true]
    )

    summary = summarize_dataset_labels_by_binding(
        labeled_stats;
        dataset_label_column = :cluster_disprot_cutoff_dataset_label
    )
    by_label = Dict(row.dataset_label => row for row in eachrow(summary))

    @test by_label["training"].n_proteins == 2
    @test by_label["training"].n_final_clusters == 1
    @test by_label["training"].n_idr_binding_proteins == 1
    @test by_label["training"].n_idr_binding_clusters == 1

    @test by_label["caid"].n_proteins == 2
    @test by_label["caid"].n_final_clusters == 1
    @test by_label["caid"].n_idr_binding_proteins == 1
    @test by_label["caid"].n_idr_binding_clusters == 1

    @test by_label["testing"].n_proteins == 2
    @test by_label["testing"].n_final_clusters == 2
    @test by_label["testing"].n_idr_binding_proteins == 1
    @test by_label["testing"].n_idr_binding_clusters == 1
end

@testitem "Summarize dataset labels by binding supports missing labels" begin
    using Test
    using IDRBindingDataset: summarize_dataset_labels_by_binding

    labeled_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_B1", "P_C1", "P_C2"],
        final_cluster_id = ["cluster_A", "cluster_B", "cluster_C", "cluster_C"],
        cluster_fixed_dataset_label = Any["training", "caid", missing, missing],
        protein_has_idr_binding_annotation = [true, false, true, false]
    )

    summary = summarize_dataset_labels_by_binding(
        labeled_stats;
        dataset_label_column = :cluster_fixed_dataset_label
    )
    by_label = Dict(row.dataset_label => row for row in eachrow(summary))

    @test by_label["training"].n_proteins == 1
    @test by_label["caid"].n_proteins == 1
    @test by_label["missing"].n_proteins == 2
    @test by_label["missing"].n_final_clusters == 1
    @test by_label["missing"].n_idr_binding_proteins == 1
    @test by_label["missing"].n_idr_binding_clusters == 1
end

@testitem "Summarize dataset labels by binding rejects inconsistent labels within cluster" begin
    using Test
    using IDRBindingDataset: summarize_dataset_labels_by_binding

    labeled_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_A2"],
        final_cluster_id = ["cluster_A", "cluster_A"],
        cluster_disprot_cutoff_dataset_label = ["training", "caid"],
        protein_has_idr_binding_annotation = [false, true]
    )

    @test_throws ArgumentError summarize_dataset_labels_by_binding(
        labeled_stats;
        dataset_label_column = :cluster_disprot_cutoff_dataset_label
    )
end

@testitem "Sweep DisProt cutoff summaries compares labels across dates" begin
    using Dates
    using Test
    using IDRBindingDataset: sweep_disprot_cutoff_label_summaries

    dataset_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1", "P_B1", "P_C1", "P_C2", "P_D1"],
        final_cluster_id = [
            "cluster_A", "cluster_B", "cluster_C", "cluster_C", "cluster_D"],
        cluster_fixed_dataset_label = Any["training", "caid", missing, missing, missing],
        protein_oldest_disprot_release_date = Any[
            Date(2019, 1, 1),
            Date(2021, 6, 1),
            Date(2022, 6, 1),
            Date(2023, 2, 1),
            missing
        ],
        protein_has_idr_binding_annotation = [true, false, true, false, true]
    )

    sweep = sweep_disprot_cutoff_label_summaries(
        dataset_stats;
        disprot_testing_cutoff_dates = [Date(2021, 12, 31), Date(2022, 12, 31)]
    )

    function _cell(date::Date, label::String, column::Symbol)
        row = only(filter(
            r -> r.disprot_testing_cutoff_date == date && r.dataset_label == label,
            eachrow(sweep)
        ))
        return getproperty(row, column)
    end

    @test _cell(Date(2021, 12, 31), "testing", :n_proteins) == 2
    @test _cell(Date(2021, 12, 31), "testing", :n_final_clusters) == 1
    @test _cell(Date(2021, 12, 31), "testing", :n_idr_binding_clusters) == 1

    @test _cell(Date(2022, 12, 31), "testing", :n_proteins) == 0
    @test _cell(Date(2022, 12, 31), "training", :n_proteins) == 4
    @test _cell(Date(2022, 12, 31), "training", :n_idr_binding_proteins) == 3
    @test _cell(Date(2022, 12, 31), "caid", :n_proteins) == 1
end

@testitem "Sweep DisProt cutoff summaries requires at least one cutoff date" begin
    using Dates
    using Test
    using IDRBindingDataset: sweep_disprot_cutoff_label_summaries

    dataset_stats = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_A1"],
        final_cluster_id = ["cluster_A"],
        cluster_fixed_dataset_label = Any[missing],
        protein_oldest_disprot_release_date = Any[missing],
        protein_has_idr_binding_annotation = [false]
    )

    @test_throws ArgumentError sweep_disprot_cutoff_label_summaries(
        dataset_stats;
        disprot_testing_cutoff_dates = Date[]
    )
end
