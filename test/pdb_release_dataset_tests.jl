using TestItems

@testitem "Collect DisProt region-level PDB refs uses regex cleaning" begin
    using Test
    using IDRBindingDataset:
                             collect_disprot_region_level_pdb_refs,
                             build_unique_uniprot_pdb_map

    disprot_entries = Any[
        Dict(
        "acc" => "P11111",
        "regions" => Any[
            Dict(
            "cross_refs" => Any[
            Dict("db" => "PDB", "id" => "1abc:a"),
            Dict("db" => "PDB", "id" => "4TTB;PDB:4TTC"),
            Dict("db" => "PDB", "id" => "invalid"),
            Dict("db" => "UNIPROT", "id" => "Q00000")
        ]
        )
        ]
    )
    ]

    refs = collect_disprot_region_level_pdb_refs(disprot_entries)
    @test size(refs, 1) == 3
    @test sort(unique(refs.pdb)) == ["1ABC", "4TTB", "4TTC"]

    uniprot_pdb_map = build_unique_uniprot_pdb_map(refs)
    @test size(uniprot_pdb_map, 1) == 3
    @test names(uniprot_pdb_map) == ["acc", "pdb"]
end

@testitem "Build SIFTS map for DisProt IDs filters and normalizes" begin
    using Test
    using IDRBindingDataset:
                             build_sifts_uniprot_pdb_map_for_disprot,
                             build_disprot_sifts_pdb_comparison

    sifts_uniprot_df = IDRBindingDataset.DataFrames.DataFrame(
        SP_PRIMARY = ["P11111", "P11111", "P22222", "Q99999"],
        PDB = ["1abc", "2def", "3ghi", "4jkl"]
    )

    sifts_map = build_sifts_uniprot_pdb_map_for_disprot(
        sifts_uniprot_df,
        ["P11111", "P22222"]
    )

    @test size(sifts_map, 1) == 3
    @test sort(unique(sifts_map.pdb)) == ["1ABC", "2DEF", "3GHI"]

    disprot_map = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P11111", "P33333"],
        pdb = ["1ABC", "9XYZ"]
    )

    comparison = build_disprot_sifts_pdb_comparison(disprot_map, sifts_map)
    @test Set(comparison.mapping_status) ==
          Set(["both", "only_disprot", "only_sifts"])
end

@testitem "Build release map from comparison classifies unresolved IDs" begin
    using Dates
    using Test
    using IDRBindingDataset: build_pdb_release_map_from_comparison

    comparison = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P11111", "P22222", "P33333", "P44444"],
        pdb = ["1ABC", "1HHB:A", "BADREF", "12TW"],
        mapping_status = ["both", "only_disprot", "only_sifts", "only_sifts"]
    )

    release_fetch_calls = Vector{Set{String}}()
    status_fetch_calls = Vector{Set{String}}()

    function fake_release_fetcher(pdb_ids; batch_size = 300)
        @test batch_size == 300
        ids = Set(String.(pdb_ids))
        push!(release_fetch_calls, ids)
        return Dict{String, Union{Missing, Date}}(
            pdb_id => if pdb_id == "1ABC"
                Date(2020, 1, 1)
            elseif pdb_id == "4HHB"
                Date(1976, 1, 1)
            else
                missing
            end for pdb_id in String.(pdb_ids)
        )
    end

    function fake_status_fetcher(pdb_ids)
        ids = Set(String.(pdb_ids))
        push!(status_fetch_calls, ids)
        return Dict(
            "1HHB" => (
                http_status = 200,
                status = "REMOVED",
                status_code = "OBS",
                replacement_pdb = "4HHB"
            ),
            "12TW" => (
                http_status = 404,
                status = missing,
                status_code = missing,
                replacement_pdb = missing
            )
        )
    end

    output = build_pdb_release_map_from_comparison(
        comparison;
        batch_size = 300,
        release_date_fetcher = fake_release_fetcher,
        holding_status_fetcher = fake_status_fetcher
    )

    @test size(output.invalid_merged_pdb_references, 1) == 1
    @test output.invalid_merged_pdb_references.pdb_reference == ["BADREF"]

    @test release_fetch_calls == [Set(["1ABC", "1HHB", "12TW"]), Set(["4HHB"])]
    @test status_fetch_calls == [Set(["12TW", "1HHB"])]

    @test size(output.disprot_pdb_map_with_release, 1) == 4
    @test size(output.failed_pdb_lookups, 1) == 2

    classifications = Dict(
        row.pdb => row.pdb_lookup_classification
    for row in eachrow(output.disprot_pdb_map_with_release)
    )
    @test classifications["1ABC"] == "current"
    @test classifications["1HHB"] == "obsolete"
    @test classifications["12TW"] == "invalid_or_not_found"
    @test classifications["4HHB"] == "current"

    replacement_rows = [row
                        for row in eachrow(output.disprot_pdb_map_with_release)
                        if row.acc == "P22222" && row.pdb == "4HHB"]
    @test length(replacement_rows) == 1
    @test replacement_rows[1].pdb_release_date == Date(1976, 1, 1)
end

@testitem "Build release map from comparison does not duplicate replacement PDB rows" begin
    using Dates
    using Test
    using IDRBindingDataset: build_pdb_release_map_from_comparison

    comparison = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P11111", "P11111"],
        pdb = ["1HHB", "4HHB"],
        mapping_status = ["only_disprot", "both"]
    )

    function fake_release_fetcher(pdb_ids; batch_size = 300)
        @test batch_size == 300
        return Dict{String, Union{Missing, Date}}(
            pdb_id => (pdb_id == "4HHB" ? Date(1976, 1, 1) : missing)
        for pdb_id in String.(pdb_ids)
        )
    end

    function fake_status_fetcher(pdb_ids)
        @test Set(String.(pdb_ids)) == Set(["1HHB"])
        return Dict(
            "1HHB" => (
            http_status = 200,
            status = "REMOVED",
            status_code = "OBS",
            replacement_pdb = "4HHB"
        )
        )
    end

    output = build_pdb_release_map_from_comparison(
        comparison;
        batch_size = 300,
        release_date_fetcher = fake_release_fetcher,
        holding_status_fetcher = fake_status_fetcher
    )

    @test size(output.disprot_pdb_map_with_release, 1) == 2
    @test count(row -> row.acc == "P11111" && row.pdb == "4HHB",
        eachrow(output.disprot_pdb_map_with_release)) == 1
end

@testitem "Write PDB release map to gzip CSV" begin
    using Dates
    using Test
    using IDRBindingDataset: write_pdb_release_map_csv_gz

    df = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P11111"],
        pdb = ["1ABC"],
        pdb_release_date = [Date(2020, 1, 1)],
        pdb_lookup_classification = ["current"],
        rcsb_http_status = [0],
        rcsb_status = [missing],
        rcsb_status_code = [missing],
        replacement_pdb = [missing]
    )

    tmp_dir = mktempdir()
    output_path = joinpath(tmp_dir, "pdb_release.csv.gz")
    write_pdb_release_map_csv_gz(df, output_path)

    @test isfile(output_path)
    @test read(output_path)[1:2] == UInt8[0x1f, 0x8b]
end

@testitem "Write PDB release README report" begin
    using Dates
    using Test
    using IDRBindingDataset:
                             write_pdb_release_map_csv_gz,
                             write_pdb_release_readme

    df = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P11111", "P22222"],
        pdb = ["1ABC", "1HHB"],
        pdb_release_date = Union{Missing, Date}[Date(2020, 1, 1), missing],
        pdb_lookup_classification = ["current", "obsolete"],
        rcsb_http_status = [0, 200],
        rcsb_status = Union{Missing, String}[missing, "REMOVED"],
        rcsb_status_code = Union{Missing, String}[missing, "OBS"],
        replacement_pdb = Union{Missing, String}[missing, "4HHB"]
    )

    failed = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P22222"],
        pdb = ["1HHB"],
        pdb_release_date = Union{Missing, Date}[missing],
        pdb_lookup_classification = ["obsolete"],
        rcsb_http_status = [200],
        rcsb_status = ["REMOVED"],
        rcsb_status_code = ["OBS"],
        replacement_pdb = ["4HHB"]
    )

    invalid_refs = IDRBindingDataset.DataFrames.DataFrame(
        acc = ["P33333"],
        pdb_reference = ["BADREF"],
        failure_reason = ["invalid_format_no_pdb_code"]
    )

    tmp_dir = mktempdir()
    output_csv = joinpath(tmp_dir, "uniprot_pdb_release_date_map.csv.gz")
    write_pdb_release_map_csv_gz(df, output_csv)

    readme_path = joinpath(tmp_dir, "README.md")
    write_pdb_release_readme(
        readme_path;
        output_csv_path = output_csv,
        disprot_pdb_map_with_release = df,
        failed_pdb_lookups = failed,
        invalid_merged_pdb_references = invalid_refs,
        batch_size = 300,
        generated_at_utc = "2026-02-18T23:59:59Z"
    )

    @test isfile(readme_path)
    content = read(readme_path, String)
    @test occursin("# Derived UniProt-PDB Release-Date Map", content)
    @test occursin("Generated at (UTC): `2026-02-18T23:59:59Z`", content)
    @test occursin("GraphQL batch size: `300`", content)
    @test occursin("Rows in cleaned UniProt-PDB map: `2`", content)
    @test occursin("Rows with missing release date: `1`", content)
    @test occursin("`P22222` | `1HHB`", content)
    @test occursin("raw_reference=`BADREF`", content)
end
