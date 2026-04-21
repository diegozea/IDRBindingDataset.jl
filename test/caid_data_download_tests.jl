using TestItems

@testitem "CAID prefix validator accepts >DP header" begin
    using Test

    tmp_dir = mktempdir()
    raw_path = joinpath(tmp_dir, "ok.fasta")
    write(raw_path, ">DP00001\nACDE\n")

    @test IDRBindingDataset._caid_fasta_prefix_validator(raw_path) === nothing
end

@testitem "CAID prefix validator rejects non-DP header" begin
    using Test

    tmp_dir = mktempdir()
    raw_path = joinpath(tmp_dir, "bad.fasta")
    write(raw_path, ">P1;SEQ1\nACDE\n")

    err = try
        IDRBindingDataset._caid_fasta_prefix_validator(raw_path)
        nothing
    catch caught
        caught
    end

    @test err isa ErrorException
    @test occursin("does not start with >DP", sprint(showerror, err))
end

@testitem "Write CAID folder README includes sources and file metadata" begin
    using Test
    using IDRBindingDataset: _write_caid_folder_readme

    tmp_dir = mktempdir()
    file_path = joinpath(tmp_dir, "binding.fasta.gz")
    write(file_path, "gzip-placeholder")
    entries = [(
        filename = "binding.fasta.gz",
        url = "https://example.org/binding.fasta",
        destination = file_path
    )]

    timestamp = "2026-02-12T12:10:00Z"
    _write_caid_folder_readme(tmp_dir, "CAID-TEST", timestamp, entries)

    readme_path = joinpath(tmp_dir, "README.md")
    content = read(readme_path, String)
    @test occursin("# CAID-TEST Raw Reference Data", content)
    @test occursin("Download timestamp (UTC): `$(timestamp)`", content)
    @test occursin("binding.fasta.gz", content)
    @test occursin("https://example.org/binding.fasta", content)
    @test !occursin("Validation:", content)
    @test !occursin("Validated records:", content)
end

@testitem "CAID-3 invalid-regions TSV validator accepts expected header" begin
    using Test
    using IDRBindingDataset: _validate_caid3_invalid_regions_tsv

    tmp_dir = mktempdir()
    path = joinpath(tmp_dir, "invalid_regions_fullCAID3.tsv")
    write(path, "disprot_id\tregion_id\tvalid\treason\nDP00001\tDP00001_r1\t0\tinvalid\n")

    @test _validate_caid3_invalid_regions_tsv(path) === nothing
end

@testitem "CAID-3 invalid-regions TSV validator rejects unexpected header" begin
    using Test
    using IDRBindingDataset: _validate_caid3_invalid_regions_tsv

    tmp_dir = mktempdir()
    path = joinpath(tmp_dir, "invalid_regions_fullCAID3.tsv")
    write(path, "disprot\tregion\tvalid\nDP00001\tDP00001_r1\t0\n")

    err = try
        _validate_caid3_invalid_regions_tsv(path)
        nothing
    catch caught
        caught
    end

    @test err isa ErrorException
    @test occursin("unexpected header", sprint(showerror, err))
end

@testitem "Write CAID folder README marks compression per file type" begin
    using Test
    using IDRBindingDataset: _write_caid_folder_readme

    tmp_dir = mktempdir()
    gz_path = joinpath(tmp_dir, "binding.fasta.gz")
    tsv_path = joinpath(tmp_dir, "invalid_regions_fullCAID3.tsv")
    write(gz_path, "gzip-placeholder")
    write(tsv_path, "disprot_id\tregion_id\tvalid\nDP00001\tDP00001_r1\t0\n")

    entries = [
        (
            filename = "binding.fasta.gz",
            url = "https://example.org/binding.fasta",
            destination = gz_path
        ),
        (
            filename = "invalid_regions_fullCAID3.tsv",
            url = "https://example.org/invalid_regions_fullCAID3.tsv",
            destination = tsv_path
        )
    ]

    _write_caid_folder_readme(tmp_dir, "CAID-3", "2026-04-19T10:00:00Z", entries)
    content = read(joinpath(tmp_dir, "README.md"), String)
    lines = split(content, '\n')
    gz_line = only(filter(line -> occursin("- `binding.fasta.gz`", line), lines))
    tsv_line = only(filter(line -> occursin("- `invalid_regions_fullCAID3.tsv`", line), lines))

    @test occursin("Final size (total):", content)
    @test !occursin("Final size (compressed, total):", content)
    @test occursin("compressed", gz_line)
    @test !occursin("compressed", tsv_line)
end

@testitem "CAID-3 invalid-regions download writes file and README" begin
    using Test
    using IDRBindingDataset: download_caid3_invalid_regions_data

    tmp_dir = mktempdir()
    existing_gz = joinpath(tmp_dir, "binding.fasta.gz")
    write(existing_gz, "gzip-placeholder")

    function fake_downloader(::AbstractString, path::AbstractString)
        write(path, "disprot_id\tregion_id\tvalid\treason\nDP00001\tDP00001_r1\t0\tinvalid\n")
        return path
    end

    download_caid3_invalid_regions_data(
        tmp_dir;
        timestamp_utc = "2026-04-19T10:00:00Z",
        logger = _ -> nothing,
        downloader = fake_downloader
    )

    tsv_path = joinpath(tmp_dir, "invalid_regions_fullCAID3.tsv")
    readme_path = joinpath(tmp_dir, "README.md")
    content = read(readme_path, String)

    @test isfile(tsv_path)
    @test isfile(readme_path)
    @test occursin("Download timestamp (UTC): `2026-04-19T10:00:00Z`", content)
    @test occursin("Final size (total):", content)
    @test occursin("binding.fasta.gz", content)
    @test occursin("invalid_regions_fullCAID3.tsv", content)
    @test occursin(
        "https://caid.idpcentral.org/assets/sections/challenge/static/references/3/binding.fasta",
        content
    )
    @test occursin(
        "https://raw.githubusercontent.com/BioComputingUP/caid-reference/main/data/disprot/invalid_regions_fullCAID3.tsv",
        content
    )
end

@testitem "CAID-3 invalid-regions download removes invalid output on failure" begin
    using Test
    using IDRBindingDataset: download_caid3_invalid_regions_data

    tmp_dir = mktempdir()
    destination = joinpath(tmp_dir, "invalid_regions_fullCAID3.tsv")

    function bad_downloader(::AbstractString, path::AbstractString)
        write(path, "disprot\tregion\tvalid\nDP00001\tDP00001_r1\t0\n")
        return path
    end

    @test_throws ErrorException download_caid3_invalid_regions_data(
        tmp_dir;
        logger = _ -> nothing,
        downloader = bad_downloader
    )
    @test !isfile(destination)
end
