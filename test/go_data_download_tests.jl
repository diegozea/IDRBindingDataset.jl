using TestItems

@testitem "GO OBO validator accepts minimal parseable ontology" begin
    using Test
    using IDRBindingDataset: _validate_go_obo_file

    minimal_go_obo = """
    format-version: 1.2
    ontology: GO

    [Term]
    id: GO:0003674
    name: molecular_function
    namespace: molecular_function
    """

    tmp_dir = mktempdir()
    path = joinpath(tmp_dir, "go-basic.obo")
    write(path, minimal_go_obo)

    @test _validate_go_obo_file(path) === nothing
end

@testitem "GO OBO validator rejects invalid file content" begin
    using Test
    using IDRBindingDataset: _validate_go_obo_file

    tmp_dir = mktempdir()
    path = joinpath(tmp_dir, "go-basic.obo")
    write(path, "not-an-obo")

    @test_throws Exception _validate_go_obo_file(path)
end

@testitem "GO downloader writes go-basic.obo.gz and README" begin
    using Test
    using IDRBindingDataset: CodecZlib, download_go_reference_data

    minimal_go_obo = """
    format-version: 1.2
    ontology: GO

    [Term]
    id: GO:0003674
    name: molecular_function
    namespace: molecular_function
    """

    tmp_raw_dir = mktempdir()
    timestamp = "2026-04-20T08:00:00Z"

    function fake_downloader(::AbstractString, path::AbstractString)
        write(path, minimal_go_obo)
        return path
    end

    download_go_reference_data(
        tmp_raw_dir;
        timestamp_utc = timestamp,
        logger = _ -> nothing,
        downloader = fake_downloader
    )

    go_dir = joinpath(tmp_raw_dir, "go")
    obo_path = joinpath(go_dir, "go-basic.obo.gz")
    readme_path = joinpath(go_dir, "README.md")
    content = read(readme_path, String)
    decompressed = open(obo_path, "r") do io
        read(CodecZlib.GzipDecompressorStream(io), String)
    end

    @test isfile(obo_path)
    @test isfile(readme_path)
    @test occursin("format-version: 1.2", decompressed)
    @test occursin("Download timestamp (UTC): `$(timestamp)`", content)
    @test occursin("go-basic.obo.gz", content)
    @test occursin("https://purl.obolibrary.org/obo/go/go-basic.obo", content)
    @test occursin("OBOParse.load(path, \"GO\")", content)
end

@testitem "GO downloader removes invalid output when validation fails" begin
    using Test
    using IDRBindingDataset: download_go_reference_data

    tmp_raw_dir = mktempdir()
    go_dir = joinpath(tmp_raw_dir, "go")
    destination = joinpath(go_dir, "go-basic.obo.gz")

    function bad_downloader(::AbstractString, path::AbstractString)
        write(path, "bad-obo-content")
        return path
    end

    @test_throws Exception download_go_reference_data(
        tmp_raw_dir;
        logger = _ -> nothing,
        downloader = bad_downloader
    )
    @test !isfile(destination)
    @test !isfile(joinpath(go_dir, "README.md"))
end
