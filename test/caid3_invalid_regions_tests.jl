using TestItems

@testitem "Load CAID-3 invalid region IDs from TSV" begin
    using Test
    using IDRBindingDataset: _load_caid3_invalid_region_ids

    tmp_dir = mktempdir()
    raw_data_dir = joinpath(tmp_dir, "raw")
    caid3_dir = joinpath(raw_data_dir, "caid_3")
    mkpath(caid3_dir)

    tsv_path = joinpath(caid3_dir, "invalid_regions_fullCAID3.tsv")
    write(
        tsv_path,
        """
disprot_id\tregion_id\tvalid\tstart\tend
DP00001\tDP00001r001\t0\t1\t10
DP00001\tDP00001r002\t1\t11\t20
DP00002\tDP00002r003\t0\t21\t30
"""
    )

    invalid_region_ids = _load_caid3_invalid_region_ids(raw_data_dir)
    @test invalid_region_ids == Set(["DP00001r001", "DP00002r003"])
end

@testitem "Load CAID-3 invalid region IDs validates TSV schema and values" begin
    using Test
    using IDRBindingDataset: _load_caid3_invalid_region_ids

    tmp_dir = mktempdir()
    raw_data_dir = joinpath(tmp_dir, "raw")
    caid3_dir = joinpath(raw_data_dir, "caid_3")
    mkpath(caid3_dir)

    bad_header_path = joinpath(caid3_dir, "invalid_regions_fullCAID3.tsv")
    write(bad_header_path, "disprot_id\tregion\tvalue\nDP00001\tDP00001r001\t0\n")
    @test_throws ArgumentError _load_caid3_invalid_region_ids(raw_data_dir)

    write(
        bad_header_path,
        "disprot_id\tregion_id\tvalid\nDP00001\tDP00001r001\t2\n"
    )
    @test_throws ArgumentError _load_caid3_invalid_region_ids(raw_data_dir)
end

@testitem "Region invalid marker checks region_id membership" begin
    using Test
    using IDRBindingDataset: _region_is_marked_invalid

    invalid_region_ids = Set(["DP00001r001"])
    @test _region_is_marked_invalid(
        Dict("region_id" => "DP00001r001"),
        invalid_region_ids
    )
    @test !_region_is_marked_invalid(
        Dict("region_id" => "DP00001r002"),
        invalid_region_ids
    )
    @test !_region_is_marked_invalid(
        Dict("start" => 1, "end" => 2),
        invalid_region_ids
    )
end
