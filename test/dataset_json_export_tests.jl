using TestItems

@testitem "Build dataset JSON records keeps CSV fields and adds sequence labels" begin
    using Dates
    using Test
    using IDRBindingDataset: build_dataset_json_records

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00001"],
        final_cluster_id = ["cluster_000001"],
        cluster_dataset_label = ["training"],
        protein_elm_clusters = [missing],
        protein_oldest_pdb_release_date = [Date(2020, 1, 1)]
    )

    disprot_entries = Any[
        Dict(
        "acc" => "P00001",
        "sequence" => "acdefg",
        "regions" => Any[
            Dict("term_id" => "IDPO:0000002", "term_name" => "disorder", "start" => 2, "end" =>
                4),
            Dict("term_id" => "GO:0005515", "start" => 3, "end" => 5)
        ]
    )
    ]

    records = build_dataset_json_records(
        final_dataset_table,
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
    @test length(records) == 1
    record = only(records)

    @test record.protein_id == "P00001"
    @test record.final_cluster_id == "cluster_000001"
    @test record.cluster_dataset_label == "training"
    @test ismissing(record.protein_elm_clusters)
    @test record.protein_oldest_pdb_release_date == Date(2020, 1, 1)
    @test record.sequence == "ACDEFG"
    @test record.residue_labels == [0, 1, 2, 2, 2, 0]
end

@testitem "Build dataset JSON records uses ontology term IDs only" begin
    using Test
    using IDRBindingDataset: build_dataset_json_records

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00002"],
        cluster_dataset_label = ["validation"]
    )

    disprot_entries = Any[
        Dict(
        "acc" => "P00002",
        "sequence" => "AAAAA",
        "regions" => Any[
            Dict("term_id" => "IDPO:0000002", "start" => 1, "end" => 2),
            Dict("term_id" => "GO:0003676", "start" => 2, "end" => 4),
            Dict(
                "term_is_binding" => "false",
                "term_name" => "protein binding",
                "start" => 5,
                "end" => 5
            )
        ]
    )
    ]

    records = build_dataset_json_records(
        final_dataset_table,
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
    @test only(records).residue_labels == [1, 2, 2, 2, 0]
end

@testitem "Build dataset JSON records uses ontology term rules when provided" begin
    using Test
    using IDRBindingDataset:
                             _DisProtChallengeTermRules,
                             build_dataset_json_records

    term_rules = _DisProtChallengeTermRules(
        Set(["GO:0005488", "GO:0005515"]),
        Set(["GO:0003676"]),
        Set(["GO:0005515"]),
        Set(["IDPO:0000002", "IDPO:0000003", "IDPO:0000004"])
    )

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00009"],
        cluster_dataset_label = ["validation"]
    )

    disprot_entries = Any[
        Dict(
        "acc" => "P00009",
        "sequence" => "AAAAAA",
        "regions" => Any[
            Dict("term_id" => "IDPO:0000011", "start" => 1, "end" => 1),
            Dict("term_id" => "IDPO:0000003", "start" => 2, "end" => 2),
            Dict(
                "term_id" => "GO:0060090",
                "term_is_binding" => true,
                "term_name" => "molecular adaptor activity",
                "start" => 3,
                "end" => 3
            ),
            Dict("term_id" => "GO:0005515", "start" => 4, "end" => 5),
            Dict("term_name" => "RNA binding", "start" => 6, "end" => 6)
        ]
    )
    ]

    records = build_dataset_json_records(
        final_dataset_table,
        disprot_entries;
        term_rules = term_rules
    )
    @test only(records).residue_labels == [0, 1, 0, 2, 2, 0]
end

@testitem "Build dataset JSON records can filter invalid region IDs" begin
    using Test
    using IDRBindingDataset:
                             _DisProtChallengeTermRules,
                             build_dataset_json_records

    term_rules = _DisProtChallengeTermRules(
        Set(["GO:0005488", "GO:0005515"]),
        Set(["GO:0003676"]),
        Set(["GO:0005515"]),
        Set(["IDPO:0000002"])
    )

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00010"],
        cluster_dataset_label = ["validation"]
    )

    disprot_entries = Any[
        Dict(
        "acc" => "P00010",
        "sequence" => "AAAAAA",
        "regions" => Any[
            Dict(
                "region_id" => "DP00010r001", "term_id" => "IDPO:0000002", "start" =>
                    1, "end" => 6),
            Dict("region_id" => "DP00010r002", "term_id" => "GO:0005515", "start" => 2, "end" =>
                5)
        ]
    )
    ]

    records = build_dataset_json_records(
        final_dataset_table,
        disprot_entries;
        term_rules = term_rules,
        invalid_region_ids = Set(["DP00010r002"])
    )
    @test only(records).residue_labels == [1, 1, 1, 1, 1, 1]
end

@testitem "Build dataset JSON records marks disorder transitions as disordered" begin
    using Test
    using IDRBindingDataset: build_dataset_json_records

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00006"],
        cluster_dataset_label = ["validation"]
    )

    disprot_entries = Any[
        Dict(
        "acc" => "P00006",
        "sequence" => "AAAAAA",
        "regions" => Any[
            Dict("term_id" => "IDPO:0000017", "start" => 1, "end" => 1),
            Dict("term_id" => "IDPO:0000011", "start" => 2, "end" => 3),
            Dict("term_id" => "IDPO:0000014", "start" => 4, "end" => 5),
            Dict("term_id" => "IDPO:9999999", "start" => 6, "end" => 6)
        ]
    )
    ]

    records = build_dataset_json_records(
        final_dataset_table,
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
    @test only(records).residue_labels == [1, 1, 1, 1, 1, 0]
end

@testitem "Build dataset JSON records fails when dataset protein is missing in DisProt" begin
    using Test
    using IDRBindingDataset: build_dataset_json_records

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P_MISSING"],
        cluster_dataset_label = ["testing"]
    )
    disprot_entries = Any[
        Dict("acc" => "P00003", "sequence" => "AAAA", "regions" => Any[])
    ]

    @test_throws ArgumentError build_dataset_json_records(
        final_dataset_table,
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
end

@testitem "Build dataset JSON records validates region bounds" begin
    using Test
    using IDRBindingDataset: build_dataset_json_records

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00004"],
        cluster_dataset_label = ["testing"]
    )
    disprot_entries = Any[
        Dict(
        "acc" => "P00004",
        "sequence" => "AAAA",
        "regions" => Any[
            Dict("term_id" => "IDPO:0000002", "start" => 3, "end" => 6)
        ]
    )
    ]

    @test_throws ArgumentError build_dataset_json_records(
        final_dataset_table,
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
end

@testitem "Write dataset JSON gzip serializes missing values as JSON null" begin
    using JSON
    using Test
    using IDRBindingDataset: build_dataset_json_records, write_dataset_json_gz

    final_dataset_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P00005"],
        cluster_dataset_label = ["training"],
        protein_elm_clusters = [missing]
    )
    disprot_entries = Any[
        Dict(
        "acc" => "P00005",
        "sequence" => "ABCDE",
        "regions" => Any[
            Dict("term_id" => "GO:0005515", "start" => 2, "end" => 3)
        ]
    )
    ]
    records = build_dataset_json_records(
        final_dataset_table,
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

    tmp_dir = mktempdir()
    output_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")
    write_dataset_json_gz(records, output_path)

    payload = open(output_path) do io
        JSON.parse(read(IDRBindingDataset.CodecZlib.GzipDecompressorStream(io), String))
    end

    @test payload isa Vector
    @test length(payload) == 1
    row = only(payload)
    @test row["protein_id"] == "P00005"
    @test row["cluster_dataset_label"] == "training"
    @test row["protein_elm_clusters"] === nothing
    @test row["sequence"] == "ABCDE"
    @test row["residue_labels"] == [0, 2, 2, 0, 0]
end

@testitem "Build dataset CSV table encodes residue labels as compact strings" begin
    using Test
    using IDRBindingDataset: build_dataset_csv_table

    records = [
        (
        protein_id = "P00007",
        cluster_dataset_label = "training",
        sequence = "ABCDE",
        residue_labels = [0, 1, 2, 2, 0],
        protein_elm_clusters = missing
    )
    ]

    table = build_dataset_csv_table(records)
    @test length(table.protein_id) == 1
    @test only(table.protein_id) == "P00007"
    @test only(table.sequence) == "ABCDE"
    @test only(table.residue_labels) == "01220"
    @test ismissing(only(table.protein_elm_clusters))
end

@testitem "Build dataset CSV table validates residue label values" begin
    using Test
    using IDRBindingDataset: build_dataset_csv_table

    invalid_records = [
        (
        protein_id = "P00008",
        sequence = "AB",
        residue_labels = [0, 3]
    )
    ]

    @test_throws ArgumentError build_dataset_csv_table(invalid_records)
end
