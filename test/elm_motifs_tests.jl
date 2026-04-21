using TestItems

@testitem "Parse excluded InstanceLogic CLI option" begin
    using Test
    using IDRBindingDataset: _parse_excluded_instance_logic

    parsed = _parse_excluded_instance_logic(
        "false positive, \"true negative\", , false positive"
    )
    @test parsed == ["false positive", "true negative"]
end

@testitem "ELM InstanceLogic filtering normalizes quoted values" begin
    using Test
    using IDRBindingDataset

    instances = IDRBindingDataset.DataFrames.DataFrame(
        ELMIdentifier = ["ELM_A", "ELM_B", "ELM_C", "ELM_D"],
        Primary_Acc = ["P1", "P2", "P3", "P4"],
        Accessions = ["P1", "P2", "P3", "P4"],
        Start = [10, 20, 30, 40],
        End = [11, 21, 31, 41],
        InstanceLogic = ["true positive", "false positive", "\"true negative\"", "unknown"]
    )

    filtered = IDRBindingDataset._filter_elm_instances_by_logic(
        instances;
        excluded_instance_logic = ["false positive", "\"\"true negative\"\""]
    )

    @test IDRBindingDataset.DataFrames.nrow(filtered) == 2
    @test Set(filtered.ELMIdentifier) == Set(["ELM_A", "ELM_D"])
end

@testitem "ELM mapping uses Accessions fallback and default filtering" begin
    using Test
    using IDRBindingDataset

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P2", "P3", "P4"]
    )
    elm_instances = IDRBindingDataset.DataFrames.DataFrame(
        ELMIdentifier = ["ELM_A", "ELM_B", "ELM_C"],
        Primary_Acc = ["P1", "P3", "P4"],
        Accessions = ["P1 P2", "P3 P3", "P4"],
        Start = [10, 20, 30],
        End = [12, 22, 32],
        InstanceLogic = ["true positive", "unknown", "false positive"]
    )
    elm_index = IDRBindingDataset.DataFrames.DataFrame(
        ELMIdentifier = ["ELM_A", "ELM_B", "ELM_C"],
        Regex = ["A..", "B..", "C.."]
    )

    motifs = IDRBindingDataset.build_cluster_elm_motifs_table(
        cluster_representatives,
        elm_instances,
        elm_index
    )

    expected = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P2", "P3"],
        regex = ["A..", "B.."],
        start = [10, 20]
    )
    expected[!, :end] = [12, 22]

    @test names(motifs) == ["protein_id", "regex", "start", "end"]
    @test motifs == expected
end

@testitem "ELM mapping fails when identifiers are missing from ELM index" begin
    using Test
    using IDRBindingDataset

    cluster_representatives = IDRBindingDataset.DataFrames.DataFrame(protein_id = ["P1"])
    elm_instances = IDRBindingDataset.DataFrames.DataFrame(
        ELMIdentifier = ["ELM_MISSING"],
        Primary_Acc = ["P1"],
        Accessions = ["P1"],
        Start = [1],
        End = [3],
        InstanceLogic = ["true positive"]
    )
    elm_index = IDRBindingDataset.DataFrames.DataFrame(
        ELMIdentifier = ["ELM_OTHER"],
        Regex = ["X"]
    )

    err = try
        IDRBindingDataset.build_cluster_elm_motifs_table(
            cluster_representatives,
            elm_instances,
            elm_index
        )
        nothing
    catch ex
        ex
    end

    @test err isa ArgumentError
    @test occursin("ELM_MISSING", sprint(showerror, err))
end

@testitem "ELM motif CSV and README writers produce reproducible artifacts" begin
    using Test
    using IDRBindingDataset

    tmp_dir = mktempdir()
    output_csv = joinpath(tmp_dir, "motifs.csv.gz")
    output_readme = joinpath(tmp_dir, "README.md")

    motif_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P1", "P2"],
        regex = ["A..", "B..", "A.."],
        start = [10, 40, 12]
    )
    motif_table[!, :end] = [12, 42, 14]

    IDRBindingDataset.write_elm_motifs_csv_gz(motif_table, output_csv)
    IDRBindingDataset.write_elm_motifs_readme(
        output_readme;
        output_csv_path = "data/scratch/elm/mapping/motifs.csv.gz",
        output_csv_filesystem_path = output_csv,
        source_cluster_representatives_path = "data/derived/clustering/cluster_representatives.csv.gz",
        source_elm_instances_path = "data/scratch/elm/instances.tsv",
        source_elms_index_path = "data/scratch/elm/elms_index.tsv",
        excluded_instance_logic = ["false positive", "true negative"],
        input_instance_rows = 100,
        filtered_instance_rows = 90,
        motif_table = motif_table,
        cluster_protein_count = 4
    )

    reloaded = IDRBindingDataset.DataFrames.DataFrame(IDRBindingDataset.CSV.File(output_csv))
    readme_text = read(output_readme, String)

    @test IDRBindingDataset.DataFrames.nrow(reloaded) == 3
    @test occursin("# ELM Motif Mapping", readme_text)
    @test occursin("false positive", readme_text)
    @test occursin("Average motifs per motif-positive protein", readme_text)
    @test occursin("julia scripts/13_mapping_elm_motifs.jl", readme_text)
end
