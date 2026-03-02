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
            Dict("term_id" => "IDPO:0000002", "term_name" => "disorder", "start" => 2, "end" => 4),
            Dict(
                "term_is_binding" => true,
                "term_name" => "protein binding",
                "start" => 3,
                "end" => 5
            )
        ]
    )
]

    records = build_dataset_json_records(final_dataset_table, disprot_entries)
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

@testitem "Build dataset JSON records supports term-name fallbacks for disorder and binding" begin
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
            Dict("term_name" => "Disorder", "start" => 1, "end" => 2),
            Dict("term_name" => "RNA binding", "start" => 2, "end" => 4),
            Dict(
                "term_is_binding" => "false",
                "term_name" => "protein binding",
                "start" => 5,
                "end" => 5
            )
        ]
    )
]

    records = build_dataset_json_records(final_dataset_table, disprot_entries)
    @test only(records).residue_labels == [1, 2, 2, 2, 0]
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
            Dict("term_name" => "order to disorder", "start" => 4, "end" => 5),
            Dict("term_name" => "order", "start" => 6, "end" => 6)
        ]
    )
]

    records = build_dataset_json_records(final_dataset_table, disprot_entries)
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

    @test_throws ArgumentError build_dataset_json_records(final_dataset_table, disprot_entries)
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

    @test_throws ArgumentError build_dataset_json_records(final_dataset_table, disprot_entries)
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
            Dict("term_name" => "protein binding", "start" => 2, "end" => 3)
        ]
    )
]
    records = build_dataset_json_records(final_dataset_table, disprot_entries)

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

@testitem "Write dataset JSON README section appends and replaces one JSON section" begin
    using Dates
    using Test
    using IDRBindingDataset: write_dataset_json_readme_section

    tmp_dir = mktempdir()
    readme_path = joinpath(tmp_dir, "README.md")
    json_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")

    write(
        readme_path,
        """
# Final Dataset

## Artifact

- Compressed CSV: `data/dataset/cluster_dataset_labels.csv.gz`
"""
    )
    write(json_path, "[]\n")

    write_dataset_json_readme_section(
        readme_path;
        output_json_path = "data/dataset/cluster_dataset_labels.json.gz",
        output_json_filesystem_path = json_path,
        source_dataset_csv_path = "data/dataset/cluster_dataset_labels.csv.gz",
        source_disprot_regions_path = "data/raw/disprot/disprot_regions.json.gz",
        n_proteins = 7,
        generated_at_utc = DateTime(2026, 3, 2, 10, 11, 12),
        script_relpath = "scripts/14_json_format.jl"
    )

    first_render = read(readme_path, String)
    @test occursin("# Final Dataset", first_render)
    @test occursin("## JSON Artifact", first_render)
    @test occursin("- Number of proteins: `7`", first_render)
    @test occursin("- Generated at (UTC): `2026-03-02 10:11:12 UTC`", first_render)
    @test occursin(
        "- DisProt regions JSON: `data/raw/disprot/disprot_regions.json.gz`",
        first_render
    )
    @test length(collect(eachmatch(r"## JSON Artifact", first_render))) == 1

    write_dataset_json_readme_section(
        readme_path;
        output_json_path = "data/dataset/cluster_dataset_labels.json.gz",
        output_json_filesystem_path = json_path,
        source_dataset_csv_path = "data/dataset/cluster_dataset_labels.csv.gz",
        source_disprot_regions_path = "data/raw/disprot/disprot_regions.json.gz",
        n_proteins = 9,
        generated_at_utc = DateTime(2026, 3, 2, 20, 21, 22),
        script_relpath = "scripts/14_json_format.jl"
    )

    second_render = read(readme_path, String)
    @test length(collect(eachmatch(r"## JSON Artifact", second_render))) == 1
    @test occursin("- Number of proteins: `9`", second_render)
    @test !occursin("- Number of proteins: `7`", second_render)
    @test occursin("- Generated at (UTC): `2026-03-02 20:21:22 UTC`", second_render)
    json_section_start = findfirst("## JSON Artifact", second_render)
    @test !isnothing(json_section_start)
    json_section_tail = second_render[first(json_section_start):end]
    @test isnothing(findnext("\n## ", json_section_tail, 2))
end

@testitem "Write dataset JSON README section requires JSON artifact file" begin
    using Test
    using IDRBindingDataset: write_dataset_json_readme_section

    tmp_dir = mktempdir()
    readme_path = joinpath(tmp_dir, "README.md")
    missing_json_path = joinpath(tmp_dir, "missing.json.gz")

    @test_throws ArgumentError write_dataset_json_readme_section(
        readme_path;
        output_json_path = "data/dataset/cluster_dataset_labels.json.gz",
        output_json_filesystem_path = missing_json_path,
        source_dataset_csv_path = "data/dataset/cluster_dataset_labels.csv.gz",
        source_disprot_regions_path = "data/raw/disprot/disprot_regions.json.gz",
        n_proteins = 1
    )
end
