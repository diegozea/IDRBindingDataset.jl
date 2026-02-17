using TestItems

@testitem "DisProt URL builder" begin
    using Test
    using IDRBindingDataset: build_disprot_regions_url

    expected = "https://disprot.org/api/search?release=2025_12&show_ambiguous=true&show_obsolete=false&format=json"
    @test build_disprot_regions_url("2025_12") == expected
end

@testitem "DisProt release parser accepts YYYY_MM tags" begin
    using Dates
    using Test
    using IDRBindingDataset: parse_disprot_release

    @test parse_disprot_release("2022_03") == Date(2022, 3, 1)
    @test parse_disprot_release("2025_12") == Date(2025, 12, 1)
end

@testitem "DisProt release parser rejects invalid tags" begin
    using Test
    using IDRBindingDataset: parse_disprot_release

    @test_throws ArgumentError parse_disprot_release("2022-03")
    @test_throws ArgumentError parse_disprot_release("2022_13")
    @test_throws ArgumentError parse_disprot_release("not_a_release")
end

@testitem "UTC timestamp formatter" begin
    using Dates
    using Test
    using IDRBindingDataset: utc_now_iso8601

    fake_clock() = DateTime(2026, 2, 11, 15, 30, 45)
    @test utc_now_iso8601(; clock = fake_clock) == "2026-02-11T15:30:45Z"
end

@testitem "File removal helper only removes listed files" begin
    using Test
    using IDRBindingDataset: remove_files

    tmp_dir = mktempdir()
    keep_file = joinpath(tmp_dir, "keep.txt")
    remove_file = joinpath(tmp_dir, "remove.txt")
    write(keep_file, "keep")
    write(remove_file, "remove")

    remove_files([remove_file, joinpath(tmp_dir, "missing.txt")])

    @test isfile(keep_file)
    @test !isfile(remove_file)
end
@testitem "Directory cleanup helper" begin
    using Test
    using IDRBindingDataset: clean_directory

    tmp_dir = mktempdir()
    target_dir = joinpath(tmp_dir, "raw")
    mkpath(joinpath(target_dir, "nested"))
    write(joinpath(target_dir, "old_file.txt"), "old")
    write(joinpath(target_dir, "nested", "old_nested.txt"), "old")

    clean_directory(target_dir)

    @test isdir(target_dir)
    @test isempty(readdir(target_dir))
end

@testitem "Download helper retries and writes destination" begin
    using Test
    using IDRBindingDataset: download_to_file_with_retries

    tmp_dir = mktempdir()
    destination = joinpath(tmp_dir, "downloaded.txt")
    attempts = Ref(0)

    function flaky_downloader(::AbstractString, path::AbstractString)
        attempts[] += 1
        if attempts[] < 3
            error("transient error")
        end
        write(path, "payload")
        return path
    end

    noop_logger(::AbstractString) = nothing

    download_to_file_with_retries(
        "https://example.org/file.txt",
        destination;
        delays = [0.0, 0.0],
        downloader = flaky_downloader,
        logger = noop_logger
    )

    @test attempts[] == 3
    @test read(destination, String) == "payload"
end

@testitem "Download helper throws after max retries" begin
    using Test
    using IDRBindingDataset: download_to_file_with_retries

    tmp_dir = mktempdir()
    destination = joinpath(tmp_dir, "downloaded.txt")
    attempts = Ref(0)

    function failing_downloader(::AbstractString, ::AbstractString)
        attempts[] += 1
        error("always fails")
    end

    @test_throws ErrorException download_to_file_with_retries(
        "https://example.org/file.txt",
        destination;
        delays = [0.0],
        downloader = failing_downloader,
        logger = _ -> nothing
    )
    @test attempts[] == 2
end

@testitem "Download report writer" begin
    using Test
    using IDRBindingDataset: write_download_report

    tmp_dir = mktempdir()
    report_path = joinpath(tmp_dir, "README.md")

    write(joinpath(tmp_dir, "present.json"), "{}")
    sources = (
        ("https://example.org/present.json", "present.json"),
        ("https://example.org/missing.tsv", "missing.tsv")
    )

    write_download_report(
        report_path;
        status = "FAILED",
        timestamp_utc = "2026-02-11T15:00:00Z",
        release = "2025_12",
        idpo_version = "0.3.0",
        sources = sources,
        error_message = "network timeout"
    )

    report = read(report_path, String)
    @test occursin("Status: `FAILED`", report)
    @test occursin("Download timestamp (UTC): `2026-02-11T15:00:00Z`", report)
    @test occursin("DisProt release: `2025_12`", report)
    @test occursin("IDPO version: `0.3.0`", report)
    @test occursin("Final size (compressed, total):", report)
    @test occursin("`present.json`", report)
    @test occursin("`missing.tsv` (missing)", report)
    @test occursin("Error: `network timeout`", report)
end

@testitem "Download to gzip helper writes compressed output" begin
    using Test
    using IDRBindingDataset: download_to_gzip_file_with_retries

    tmp_dir = mktempdir()
    destination = joinpath(tmp_dir, "downloaded.json.gz")

    function fake_downloader(::AbstractString, path::AbstractString)
        write(path, "{\"k\":\"v\"}")
        return path
    end

    download_to_gzip_file_with_retries(
        "https://example.org/payload.json",
        destination;
        delays = [0.0],
        downloader = fake_downloader,
        logger = _ -> nothing
    )

    header = read(destination)[1:2]
    @test header == UInt8[0x1f, 0x8b]
end

@testitem "Download to gzip helper runs precompress validator" begin
    using Test
    using IDRBindingDataset: download_to_gzip_file_with_retries

    tmp_dir = mktempdir()
    destination = joinpath(tmp_dir, "downloaded.txt.gz")
    validator_calls = Ref(0)

    function fake_downloader(::AbstractString, path::AbstractString)
        write(path, "raw-payload")
        return path
    end

    function validator(path::AbstractString)
        validator_calls[] += 1
        @test isfile(path)
        @test read(path, String) == "raw-payload"
        return nothing
    end

    download_to_gzip_file_with_retries(
        "https://example.org/payload.txt",
        destination;
        delays = [0.0],
        downloader = fake_downloader,
        logger = _ -> nothing,
        precompress_validator = validator
    )

    @test validator_calls[] == 1
    @test isfile(destination)
end

@testitem "Download to gzip helper fails when precompress validator fails" begin
    using Test
    using IDRBindingDataset: download_to_gzip_file_with_retries

    tmp_dir = mktempdir()
    destination = joinpath(tmp_dir, "downloaded.txt.gz")

    function fake_downloader(::AbstractString, path::AbstractString)
        write(path, "raw-payload")
        return path
    end

    function failing_validator(::AbstractString)
        error("invalid header")
    end

    @test_throws ErrorException download_to_gzip_file_with_retries(
        "https://example.org/payload.txt",
        destination;
        delays = [0.0],
        downloader = fake_downloader,
        logger = _ -> nothing,
        precompress_validator = failing_validator
    )
    @test !isfile(destination)
end

@testitem "Git/GitHub file-size validator accepts files within limit" begin
    using Test
    using IDRBindingDataset: ensure_files_within_size_limit

    tmp_dir = mktempdir()
    small_file = joinpath(tmp_dir, "small.bin")
    write(small_file, "12345")

    @test ensure_files_within_size_limit(
        [small_file];
        max_bytes = 10,
        limit_name = "test limit"
    ) === nothing
end

@testitem "Git/GitHub file-size validator fails with clear message" begin
    using Test
    using IDRBindingDataset: ensure_files_within_size_limit

    tmp_dir = mktempdir()
    big_file = joinpath(tmp_dir, "too_big.bin")
    write(big_file, "0123456789ABCDEF")

    err = try
        ensure_files_within_size_limit(
            [big_file];
            max_bytes = 8,
            limit_name = "GitHub single-file hard limit",
            remediation = "Use Git LFS."
        )
        nothing
    catch caught
        caught
    end

    @test err isa ErrorException
    message = sprint(showerror, err)
    @test occursin("exceed the configured Git/GitHub size limit", message)
    @test occursin("GitHub single-file hard limit", message)
    @test occursin("too_big.bin", message)
    @test occursin("Use Git LFS.", message)
end
