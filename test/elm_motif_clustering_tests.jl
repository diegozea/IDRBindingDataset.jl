using TestItems

@testitem "Parse optional score threshold CLI option" begin
    using Test
    using IDRBindingDataset: _parse_optional_score_threshold

    @test _parse_optional_score_threshold("") === nothing
    @test _parse_optional_score_threshold("nothing") === nothing
    @test _parse_optional_score_threshold("  Nothing  ") === nothing
    @test _parse_optional_score_threshold("0.75") ≈ 0.75
    @test_throws ArgumentError _parse_optional_score_threshold("invalid")
end

@testitem "Regex motif clustering builds thresholded connected components" begin
    using Test
    using IDRBindingDataset

    regexes = ["A", "B", "C", "D"]
    regex_similarity_table = IDRBindingDataset.DataFrames.DataFrame(
        query = ["A", "B", "B", "C", "A", "D", "A"],
        search = ["B", "A", "C", "B", "D", "A", "A"],
        matched = [true, true, true, true, true, true, true],
        normalized_ic = [0.6, 0.6, 0.7, 0.7, 0.55, 0.55, 1.0],
        match_ic = [1.6, 1.6, 1.8, 1.8, 1.4, 1.4, 2.0],
        score = [0.8, 0.8, 0.9, 0.9, 0.4, 0.4, 1.0]
    )

    clusters = IDRBindingDataset.build_regex_motif_clusters_table(
        regexes,
        regex_similarity_table;
        normic_threshold = 0.5,
        matchic_threshold = 1.5
    )

    cluster_by_regex = Dict(clusters.regex .=> clusters.elm_motif_cluster)
    @test cluster_by_regex["A"] == cluster_by_regex["B"]
    @test cluster_by_regex["B"] == cluster_by_regex["C"]
    @test cluster_by_regex["D"] != cluster_by_regex["A"]
    @test length(unique(clusters.elm_motif_cluster)) == 2
end

@testitem "Regex motif clustering applies optional score threshold" begin
    using Test
    using IDRBindingDataset

    regexes = ["A", "B", "C"]
    regex_similarity_table = IDRBindingDataset.DataFrames.DataFrame(
        query = ["A", "B"],
        search = ["B", "C"],
        matched = [true, true],
        normalized_ic = [1.0, 1.0],
        match_ic = [3.0, 3.0],
        score = [0.6, 0.9]
    )

    clusters = IDRBindingDataset.build_regex_motif_clusters_table(
        regexes,
        regex_similarity_table;
        normic_threshold = 0.5,
        matchic_threshold = 1.5,
        score_threshold = 0.7
    )

    cluster_by_regex = Dict(clusters.regex .=> clusters.elm_motif_cluster)
    @test cluster_by_regex["B"] == cluster_by_regex["C"]
    @test cluster_by_regex["A"] != cluster_by_regex["B"]
end

@testitem "Protein mapping keeps one row per protein and motif cluster" begin
    using Test
    using IDRBindingDataset

    motif_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P1", "P1", "P2", "P2"],
        regex = ["A", "B", "D", "C", "C"],
        start = [1, 2, 3, 4, 4]
    )
    motif_table[!, :end] = [1, 2, 3, 4, 4]
    regex_cluster_table = IDRBindingDataset.DataFrames.DataFrame(
        regex = ["A", "B", "C", "D"],
        elm_motif_cluster = [1, 1, 1, 2]
    )

    protein_cluster_table = IDRBindingDataset.build_protein_elm_motif_cluster_table(
        motif_table,
        regex_cluster_table
    )

    expected = Set([("P1", 1), ("P1", 2), ("P2", 1)])
    observed = Set(zip(protein_cluster_table.protein_id, protein_cluster_table.elm_motif_cluster))

    @test observed == expected
    @test IDRBindingDataset.DataFrames.nrow(protein_cluster_table) == 3
end

@testitem "Regex/protein clustering writers generate reproducible files" begin
    using Test
    using IDRBindingDataset

    tmp_dir = mktempdir()
    regex_csv_path = joinpath(tmp_dir, "regex.csv.gz")
    protein_csv_path = joinpath(tmp_dir, "elm_motif_clustering.csv.gz")
    readme_path = joinpath(tmp_dir, "README.md")

    regex_cluster_table = IDRBindingDataset.DataFrames.DataFrame(
        regex = ["A", "B", "C"],
        elm_motif_cluster = [1, 1, 2]
    )
    protein_cluster_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P1", "P2"],
        elm_motif_cluster = [1, 2, 1]
    )

    IDRBindingDataset.write_regex_motif_clusters_csv_gz(regex_cluster_table, regex_csv_path)
    IDRBindingDataset.write_elm_motif_clustering_csv_gz(protein_cluster_table, protein_csv_path)
    IDRBindingDataset.write_elm_motif_clustering_readme(
        readme_path;
        output_csv_path = "data/derived/elm/elm_motif_clustering.csv.gz",
        output_csv_filesystem_path = protein_csv_path,
        source_motifs_path = "data/scratch/elm/mapping/motifs.csv.gz",
        source_regex_clusters_path = "data/scratch/elm/clustering/regex.csv.gz",
        min_shared_positions = 2,
        normic_threshold = 0.5,
        matchic_threshold = 1.5,
        score_threshold = 0.7,
        regex_cluster_table = regex_cluster_table,
        protein_cluster_table = protein_cluster_table
    )

    reloaded_regex = IDRBindingDataset.DataFrames.DataFrame(IDRBindingDataset.CSV.File(regex_csv_path))
    reloaded_protein = IDRBindingDataset.DataFrames.DataFrame(IDRBindingDataset.CSV.File(protein_csv_path))
    readme_text = read(readme_path, String)

    @test IDRBindingDataset.DataFrames.nrow(reloaded_regex) == 3
    @test IDRBindingDataset.DataFrames.nrow(reloaded_protein) == 3
    @test occursin("# ELM Motif Clustering", readme_text)
    @test occursin("at least\n  `2` matched positions", readme_text)
    @test occursin("NormIC >= 0.5", readme_text)
    @test occursin("MatchIC >= 1.5", readme_text)
    @test occursin("Score >= 0.7", readme_text)
    @test occursin("Proteins mapped to multiple motif clusters: `1`", readme_text)
    @test occursin("data/scratch/elm/clustering/regex.csv.gz", readme_text)
end
