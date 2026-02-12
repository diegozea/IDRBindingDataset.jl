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
