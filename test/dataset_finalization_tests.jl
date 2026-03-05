using TestItems

@testitem "Build protein ELM clusters table joins sorted unique cluster IDs" begin
    using Test
    using IDRBindingDataset: build_protein_elm_clusters_table

    protein_elm_motif_clusters = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P1", "P1", "P2", "P2", "P3", "P3"],
        elm_motif_cluster = Any[3, 1, 3, "10", "2", missing, ""]
    )

    protein_elm_clusters = build_protein_elm_clusters_table(protein_elm_motif_clusters)
    by_protein = Dict(row.protein_id => row for row in eachrow(protein_elm_clusters))

    @test names(protein_elm_clusters) == ["protein_id", "protein_elm_clusters"]
    @test by_protein["P1"].protein_elm_clusters == "1-3"
    @test by_protein["P2"].protein_elm_clusters == "2-10"
    @test by_protein["P3"].protein_elm_clusters === missing
end

@testitem "Build protein ELM clusters table rejects missing or empty protein IDs" begin
    using Test
    using IDRBindingDataset: build_protein_elm_clusters_table

    missing_protein_id = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = Any["P1", missing],
        elm_motif_cluster = ["1", "2"]
    )
    @test_throws ArgumentError build_protein_elm_clusters_table(missing_protein_id)

    empty_protein_id = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "   "],
        elm_motif_cluster = ["1", "2"]
    )
    @test_throws ArgumentError build_protein_elm_clusters_table(empty_protein_id)
end

@testitem "Build final cluster dataset table maps caid to testing and missing to validation" begin
    using Test
    using IDRBindingDataset: build_final_cluster_dataset_table

    fixed_labels_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P2", "P3", "P4"],
        final_cluster_id = ["cluster_A", "cluster_B", "cluster_B", "cluster_C"],
        cluster_fixed_dataset_label = Any["training", "caid", missing, "testing"]
    )
    protein_elm_clusters = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P2", "P4"],
        protein_elm_clusters = ["1-2", "3", "4"]
    )

    final_dataset_table = build_final_cluster_dataset_table(
        fixed_labels_table,
        protein_elm_clusters;
        missing_label = "validation"
    )
    by_protein = Dict(row.protein_id => row for row in eachrow(final_dataset_table))

    @test by_protein["P1"].cluster_dataset_label == "training"
    @test by_protein["P2"].cluster_dataset_label == "testing"
    @test by_protein["P3"].cluster_dataset_label == "validation"
    @test by_protein["P4"].cluster_dataset_label == "testing"
    @test by_protein["P1"].cluster_fixed_dataset_label == "training"
    @test by_protein["P2"].cluster_fixed_dataset_label == "caid"
    @test by_protein["P3"].cluster_fixed_dataset_label === missing
    @test by_protein["P1"].protein_elm_clusters == "1-2"
    @test by_protein["P2"].protein_elm_clusters == "3"
    @test by_protein["P3"].protein_elm_clusters === missing
    @test by_protein["P4"].protein_elm_clusters == "4"
    @test Set(final_dataset_table.cluster_dataset_label) ==
          Set(["training", "validation", "testing"])
end

@testitem "Build final cluster dataset table rejects unsupported labels" begin
    using Test
    using IDRBindingDataset: build_final_cluster_dataset_table

    protein_elm_clusters = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1"],
        protein_elm_clusters = ["1"]
    )

    invalid_fixed_labels_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1"],
        final_cluster_id = ["cluster_A"],
        cluster_fixed_dataset_label = Any["holdout"]
    )
    @test_throws ArgumentError build_final_cluster_dataset_table(
        invalid_fixed_labels_table,
        protein_elm_clusters;
        missing_label = "validation"
    )

    valid_fixed_labels_table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1"],
        final_cluster_id = ["cluster_A"],
        cluster_fixed_dataset_label = Any[missing]
    )
    @test_throws ArgumentError build_final_cluster_dataset_table(
        valid_fixed_labels_table,
        protein_elm_clusters;
        missing_label = "holdout"
    )
end
