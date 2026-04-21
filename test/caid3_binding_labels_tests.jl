using TestItems

@testitem "CAID3 binding labels check summarizes missing and mismatched proteins" begin
    using JSON
    using Test
    using IDRBindingDataset:
                             check_caid3_binding_idr_labels

    tmp_dir = mktempdir()
    binding_path = joinpath(tmp_dir, "binding_idr.fasta.gz")
    dataset_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")

    binding_content = """
>DP00001
ACDE
-01-
>DP00002
MK
10
>DP00003
AA
--
"""
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, binding_path, "w") do io
        write(io, binding_content)
    end

    dataset_records = [
        Dict(
            "protein_disprot_id" => "DP00001",
            "residue_labels" => [0, 1, 2, 0]
        ),
        Dict(
            "protein_disprot_id" => "DP00002",
            "residue_labels" => [2, 0]
        )
    ]
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, dataset_path, "w") do io
        write(io, JSON.json(dataset_records))
    end

    result = check_caid3_binding_idr_labels(
        binding_path,
        dataset_path;
        max_mismatch_examples = 5
    )

    @test result.n_binding_idr_proteins == 3
    @test result.n_dataset_records_with_disprot_id == 2
    @test result.n_proteins_found == 2
    @test result.n_proteins_passed == 1
    @test result.missing_proteins == ["DP00003"]
    @test isempty(result.length_mismatches)
    @test length(result.label_mismatches) == 1
    @test first(result.label_mismatches).protein_id == "DP00002"
    @test first(result.label_mismatches).mismatch_count == 1
    @test occursin("2:1!=0", first(result.label_mismatches).sample_positions)
    @test result.residues_compared == 6
    @test result.residues_mismatched == 1
    @test result.passed == false
end

@testitem "CAID3 binding labels check reports length mismatches" begin
    using JSON
    using Test
    using IDRBindingDataset:
                             check_caid3_binding_idr_labels

    tmp_dir = mktempdir()
    binding_path = joinpath(tmp_dir, "binding_idr.fasta.gz")
    dataset_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")

    binding_content = """
>DP10001
ABC
-01
"""
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, binding_path, "w") do io
        write(io, binding_content)
    end

    dataset_records = [
        Dict(
        "protein_disprot_id" => "DP10001",
        "residue_labels" => [0, 1]
    )
    ]
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, dataset_path, "w") do io
        write(io, JSON.json(dataset_records))
    end

    result = check_caid3_binding_idr_labels(binding_path, dataset_path)

    @test result.n_binding_idr_proteins == 1
    @test result.n_proteins_found == 1
    @test result.n_proteins_passed == 0
    @test isempty(result.missing_proteins)
    @test length(result.length_mismatches) == 1
    @test first(result.length_mismatches).protein_id == "DP10001"
    @test first(result.length_mismatches).expected_length == 3
    @test first(result.length_mismatches).observed_length == 2
    @test isempty(result.label_mismatches)
    @test result.residues_compared == 0
    @test result.residues_mismatched == 0
    @test result.passed == false
end

@testitem "CAID3 binding labels report writer renders summary and sections" begin
    using Test
    using IDRBindingDataset:
                             write_caid3_binding_labels_report

    tmp_dir = mktempdir()
    output_report_path = joinpath(tmp_dir, "CAID3_LABELS.md")

    check_result = (
        n_binding_idr_proteins = 2,
        n_dataset_records_with_disprot_id = 10,
        n_proteins_found = 2,
        n_proteins_passed = 2,
        missing_proteins = String[],
        length_mismatches = NamedTuple{
            (:protein_id, :expected_length, :observed_length), Tuple{String, Int, Int}}[],
        label_mismatches = NamedTuple{
            (:protein_id, :mismatch_count, :sample_positions), Tuple{String, Int, String}}[],
        residues_compared = 20,
        residues_mismatched = 0,
        passed = true
    )

    write_caid3_binding_labels_report(
        output_report_path;
        binding_idr_fasta_gz_path = "data/raw/caid_3/binding_idr.fasta.gz",
        dataset_json_gz_path = "data/dataset/cluster_dataset_labels.json.gz",
        check_result = check_result,
        script_relpath = "scripts/16_check_idr_binding.jl"
    )

    report = read(output_report_path, String)
    @test occursin("# CAID-3 Binding IDR Label Check", report)
    @test occursin("`-` -> `0`", report)
    @test occursin("Overall result: `PASS`", report)
    @test occursin("## Missing Proteins", report)
    @test occursin("## Length Mismatches", report)
    @test occursin("## Residue Label Mismatches", report)
end

@testitem "CAID3 binding labels check normalizes identifier to first header token" begin
    using JSON
    using Test
    using IDRBindingDataset:
                             check_caid3_binding_idr_labels

    tmp_dir = mktempdir()
    binding_path = joinpath(tmp_dir, "binding_idr.fasta.gz")
    dataset_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")

    binding_content = """
>DP20001 extra metadata
AA
01
"""
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, binding_path, "w") do io
        write(io, binding_content)
    end

    dataset_records = [
        Dict(
        "protein_disprot_id" => "DP20001",
        "residue_labels" => [1, 2]
    )
    ]
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, dataset_path, "w") do io
        write(io, JSON.json(dataset_records))
    end

    result = check_caid3_binding_idr_labels(binding_path, dataset_path)

    @test result.n_binding_idr_proteins == 1
    @test result.n_proteins_found == 1
    @test result.n_proteins_passed == 1
    @test isempty(result.missing_proteins)
    @test isempty(result.length_mismatches)
    @test isempty(result.label_mismatches)
    @test result.passed == true
end

@testitem "CAID3 binding labels check rejects duplicate protein identifiers" begin
    using JSON
    using Test
    using IDRBindingDataset:
                             check_caid3_binding_idr_labels

    tmp_dir = mktempdir()
    binding_path = joinpath(tmp_dir, "binding_idr.fasta.gz")
    dataset_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")

    binding_content = """
>DP30001
AA
11
>DP30001
AA
00
"""
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, binding_path, "w") do io
        write(io, binding_content)
    end

    dataset_records = [
        Dict(
        "protein_disprot_id" => "DP30001",
        "residue_labels" => [2, 2]
    )
    ]
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, dataset_path, "w") do io
        write(io, JSON.json(dataset_records))
    end

    err = try
        check_caid3_binding_idr_labels(binding_path, dataset_path)
        nothing
    catch caught_err
        caught_err
    end

    @test err isa ArgumentError
    @test occursin("Duplicate protein identifier", sprint(showerror, err))
end

@testitem "CAID3 binding labels check requires Feature annotation" begin
    using JSON
    using Test
    using IDRBindingDataset:
                             check_caid3_binding_idr_labels

    tmp_dir = mktempdir()
    binding_path = joinpath(tmp_dir, "binding_idr.fasta.gz")
    dataset_path = joinpath(tmp_dir, "cluster_dataset_labels.json.gz")

    binding_content = """
# Format:
# >ProteinID
# Sequence
# Binding
>DP40001
AA
11
"""
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, binding_path, "w") do io
        write(io, binding_content)
    end

    dataset_records = [
        Dict(
        "protein_disprot_id" => "DP40001",
        "residue_labels" => [2, 2]
    )
    ]
    open(IDRBindingDataset.CodecZlib.GzipCompressorStream, dataset_path, "w") do io
        write(io, JSON.json(dataset_records))
    end

    err = try
        check_caid3_binding_idr_labels(binding_path, dataset_path)
        nothing
    catch caught_err
        caught_err
    end

    @test err isa ArgumentError
    @test occursin("feature \"Feature\" is missing or empty", sprint(showerror, err))
end
