using TestItems

@testitem "DisProt full protein FASTA is extracted from JSON" begin
    using JSON
    using Test
    using IDRBindingDataset

    tmp_dir = mktempdir()
    input_json_gz = joinpath(tmp_dir, "disprot_regions.json.gz")
    output_fasta_gz = joinpath(tmp_dir, "disprot_full_proteins.fasta.gz")

    payload = Dict(
        "size" => 2,
        "data" => Any[
            Dict("acc" => "P11111", "sequence" => "acdefghik"),
            Dict("acc" => "P22222", "sequence" => "MNAPQRST")
        ]
    )

    open(input_json_gz, "w") do io
        gzip_io = IDRBindingDataset.CodecZlib.GzipCompressorStream(io)
        try
            write(gzip_io, JSON.json(payload))
        finally
            close(gzip_io)
        end
    end

    records_written = IDRBindingDataset.write_disprot_full_protein_fasta_gz(
        input_json_gz,
        output_fasta_gz;
        line_width = 5
    )

    @test records_written == 2
    @test isfile(output_fasta_gz)

    sequences = IDRBindingDataset.MIToS.MSA.read_file(
        output_fasta_gz,
        IDRBindingDataset.MIToS.MSA.FASTASequences
    )
    sequence_by_id = Dict(
        IDRBindingDataset.MIToS.MSA.sequence_id(sequence) => join(sequence)
    for sequence in sequences
    )

    @test length(sequences) == 2
    @test haskey(sequence_by_id, "P11111")
    @test haskey(sequence_by_id, "P22222")
    @test sequence_by_id["P11111"] == "ACDEFGHIK"
    @test sequence_by_id["P22222"] == "MNAPQRST"
end

@testitem "DisProt sequence extraction validates duplicate accessions" begin
    using Test
    using IDRBindingDataset

    entries = Any[
    Dict("acc" => "P11111", "sequence" => "ACDE"),
    Dict(
        "acc" => "P11111", "sequence" => "FGHI")
]

    @test_throws ArgumentError IDRBindingDataset.extract_disprot_full_protein_sequences(entries)
end

@testitem "Sequence clustering local filters enforce thresholds" begin
    using Test
    using IDRBindingDataset

    local_hit = (
        qseqid = "Q1",
        sseqid = "S1",
        pident = 80.0,
        length = 10,
        gapopen = 0,
        qstart = 1,
        qend = 10,
        sstart = 1,
        send = 10,
        qlen = 100,
        slen = 100
    )

    @test IDRBindingDataset._passes_local_filter(local_hit)
    @test !IDRBindingDataset._passes_local_filter(local_hit; local_min_identity = 80.1)
    @test !IDRBindingDataset._passes_local_filter(local_hit; local_min_ungapped_length = 11)
    @test !IDRBindingDataset._passes_local_filter(merge(local_hit, (gapopen = 1,)))
end

@testitem "DIAMOND command builder accepts AbstractString vectors" begin
    using Test
    using IDRBindingDataset

    raw_args = "makedb --in input.fasta.gz --db example_db"
    args = split(raw_args)
    command = IDRBindingDataset._diamond_cmd(args)

    @test command.exec[(end - 4):end] ==
          ["makedb", "--in", "input.fasta.gz", "--db", "example_db"]
end

@testitem "DIAMOND argument builders include requested flags" begin
    using Test
    using IDRBindingDataset

    blastp_args = IDRBindingDataset._build_all_vs_all_blastp_args(
        "input.fasta.gz",
        "db_prefix",
        "output.tsv";
        threads = 4
    )
    @test "gapopen" in blastp_args
    @test "--ultra-sensitive" in blastp_args
    @test "--masking" in blastp_args
    @test "0" in blastp_args

    global_args = IDRBindingDataset._build_global_cluster_args(
        "input.fasta.gz",
        "pairs.tsv";
        threads = 4,
        memory_limit = "16G",
        global_min_identity = 30.0,
        global_min_coverage = 80.0
    )
    @test "-d" in global_args
    @test !("--db" in global_args)
    @test "--masking" in global_args
    @test "0" in global_args

    recluster_args = IDRBindingDataset._build_global_recluster_args(
        "input.fasta.gz",
        "cluster_pairs.tsv",
        "recluster_pairs.tsv";
        threads = 4,
        memory_limit = "16G",
        global_min_identity = 30.0,
        global_min_coverage = 80.0
    )
    @test first(recluster_args) == "recluster"
    @test "-d" in recluster_args
    @test "--clusters" in recluster_args
    @test "cluster_pairs.tsv" in recluster_args
    @test "--masking" in recluster_args
    @test "0" in recluster_args
end

@testitem "Checked command helper raises for non-zero exit" begin
    using Test
    using IDRBindingDataset

    failing_cmd = Cmd(["bash", "-lc", "exit 7"])
    caught = try
        IDRBindingDataset._run_checked_command(failing_cmd, "unit-test")
        nothing
    catch ex
        ex
    end

    @test caught isa ErrorException
    @test occursin("exit code 7", sprint(showerror, caught))
end

@testitem "Local clustering from hits is deterministic and keeps singletons" begin
    using Test
    using IDRBindingDataset

    sequence_ids = ["A", "B", "C", "D"]
    hits = [
        (
            qseqid = "A",
            sseqid = "B",
            pident = 85.0,
            length = 90,
            gapopen = 0,
            qstart = 1,
            qend = 90,
            sstart = 1,
            send = 90,
            qlen = 100,
            slen = 100
        ),
        (
            qseqid = "B",
            sseqid = "C",
            pident = 84.0,
            length = 90,
            gapopen = 0,
            qstart = 1,
            qend = 90,
            sstart = 1,
            send = 90,
            qlen = 100,
            slen = 100
        )
    ]

    assignments = IDRBindingDataset._local_cluster_assignments_from_hits(sequence_ids, hits)
    @test Set(assignments.sequence_id) == Set(sequence_ids)
    @test maximum(assignments.cluster_size) == 3

    cluster_for_a = assignments[assignments.sequence_id .== "A", :cluster_id][1]
    cluster_for_c = assignments[assignments.sequence_id .== "C", :cluster_id][1]
    cluster_for_d = assignments[assignments.sequence_id .== "D", :cluster_id][1]

    @test cluster_for_a == cluster_for_c
    @test cluster_for_d != cluster_for_a
end

@testitem "Cluster pair conversion includes sequences missing from pairs" begin
    using Test
    using IDRBindingDataset

    sequence_ids = ["M2", "X", "R1", "M1"]
    pairs = [("R1", "M1"), ("M1", "M2")]

    assignments = IDRBindingDataset._cluster_assignments_from_diamond_cluster_pairs(sequence_ids, pairs)
    @test Set(assignments.sequence_id) == Set(sequence_ids)

    cluster_for_x = assignments[assignments.sequence_id .== "X", :cluster_size][1]
    cluster_for_m2 = assignments[assignments.sequence_id .== "M2", :cluster_size][1]
    @test cluster_for_x == 1
    @test cluster_for_m2 == 3
end

@testitem "Cluster representative table adds deterministic final cluster IDs" begin
    using Test
    using IDRBindingDataset

    global_assignments = IDRBindingDataset.DataFrames.DataFrame(
        cluster_id = ["G000002", "G000001", "G000002", "G000001"],
        representative_id = ["P3", "P1", "P3", "P1"],
        sequence_id = ["P3", "P1", "P4", "P2"],
        cluster_size = [2, 2, 2, 2]
    )
    local_assignments = IDRBindingDataset.DataFrames.DataFrame(
        cluster_id = ["L000002", "L000001", "L000002", "L000001"],
        representative_id = ["P2", "P1", "P2", "P1"],
        sequence_id = ["P4", "P1", "P2", "P3"],
        cluster_size = [2, 2, 2, 2]
    )

    representatives = IDRBindingDataset.build_cluster_representatives_table(
        global_assignments,
        local_assignments
    )

    @test names(representatives) == [
        "protein_id",
        "global_representative_id",
        "local_representative_id",
        "final_cluster_id"
    ]
    @test representatives.protein_id == ["P1", "P2", "P3", "P4"]
    @test all(representatives.final_cluster_id .== "cluster_000001")
end

@testitem "Cluster representative CSV writer and README writer" begin
    using Test
    using IDRBindingDataset

    tmp_dir = mktempdir()
    representatives_path = joinpath(tmp_dir, "representatives.csv.gz")
    readme_path = joinpath(tmp_dir, "README.md")
    source_path = "data/raw/disprot/disprot_regions.json.gz"

    global_assignments = IDRBindingDataset.DataFrames.DataFrame(
        cluster_id = ["G000001", "G000001", "G000002"],
        representative_id = ["P1", "P1", "P3"],
        sequence_id = ["P1", "P2", "P3"],
        cluster_size = [2, 2, 1]
    )
    local_assignments = IDRBindingDataset.DataFrames.DataFrame(
        cluster_id = ["L000001", "L000002", "L000003"],
        representative_id = ["P1", "P2", "P3"],
        sequence_id = ["P1", "P2", "P3"],
        cluster_size = [1, 1, 1]
    )

    representatives = IDRBindingDataset.build_cluster_representatives_table(
        global_assignments,
        local_assignments
    )
    IDRBindingDataset.write_cluster_representatives_csv(representatives, representatives_path)

    IDRBindingDataset.write_sequence_clustering_readme(
        readme_path;
        source_disprot_regions_path = source_path,
        source_disprot_regions_size_bytes = 123,
        disprot_fasta_path = "data/derived/clustering/disprot_full_proteins.fasta.gz",
        disprot_fasta_size_bytes = 456,
        representatives_path = "data/derived/clustering/cluster_representatives.csv.gz",
        representatives_size_bytes = 789,
        sequence_count = 3,
        global_cluster_count = 2,
        local_cluster_count = 3,
        final_cluster_count = 2,
        global_method = "diamond_cluster_recluster",
        local_method = "all_vs_all_blastp",
        global_min_identity = 30.0,
        global_min_coverage = 80.0,
        local_min_identity = 80.0,
        local_min_ungapped_length = 10,
        generated_at_utc = "2026-02-19T00:00:00Z"
    )

    @test isfile(representatives_path)
    @test isfile(readme_path)

    representatives_df = IDRBindingDataset.CSV.read(
        representatives_path,
        IDRBindingDataset.DataFrames.DataFrame
    )
    @test names(representatives_df) == [
        "protein_id",
        "global_representative_id",
        "local_representative_id",
        "final_cluster_id"
    ]

    readme_contents = read(readme_path, String)
    @test occursin("## Outputs", readme_contents)
    @test !occursin("## Input", readme_contents)
    @test occursin("## Sequence and Cluster Counts", readme_contents)
    @test occursin(
        "DisProt full-protein FASTA: data/derived/clustering/disprot_full_proteins.fasta.gz (456 bytes, 0.0005 MB)",
        readme_contents
    )
    @test occursin(
        "Cluster representatives CSV: data/derived/clustering/cluster_representatives.csv.gz (789 bytes, 0.0008 MB)",
        readme_contents
    )
    @test occursin("DisProt regions JSON: $(source_path) (123 bytes, 0.0001 MB)", readme_contents)
    @test occursin("Global cluster count: 2", readme_contents)
    @test occursin("Local cluster count: 3", readme_contents)
    @test occursin("Final cluster count: 2", readme_contents)
end

@testitem "FASTA ID reader behavior for duplicate and empty IDs" begin
    using Test
    using IDRBindingDataset

    tmp_dir = mktempdir()

    duplicate_fasta_path = joinpath(tmp_dir, "duplicate.fasta")
    write(duplicate_fasta_path, ">P1\nAAAA\n>P2\nBBBB\n>P1\nCCCC\n")
    duplicate_ids = IDRBindingDataset._read_fasta_ids(duplicate_fasta_path)
    @test duplicate_ids == ["P1", "P2", "P1(1)"]

    empty_id_fasta_path = joinpath(tmp_dir, "empty_id.fasta")
    write(empty_id_fasta_path, ">\nAAAA\n")
    @test_throws ErrorException IDRBindingDataset._read_fasta_ids(empty_id_fasta_path)
end
