using TestItems

@testitem "Shared helpers write gzip CSV and plain text outputs" begin
    using Test
    using IDRBindingDataset: _write_csv_gz, _write_text_file

    tmp_dir = mktempdir()

    table_path = joinpath(tmp_dir, "nested", "table.csv.gz")
    table = IDRBindingDataset.DataFrames.DataFrame(
        protein_id = ["P1", "P2"],
        cluster_id = ["cluster_1", "cluster_2"]
    )
    _write_csv_gz(table, table_path)

    @test isfile(table_path)
    parsed = open(table_path, "r") do io
        IDRBindingDataset.CSV.read(
            IDRBindingDataset.CodecZlib.GzipDecompressorStream(io),
            IDRBindingDataset.DataFrames.DataFrame
        )
    end
    @test parsed == table

    text_path = joinpath(tmp_dir, "notes", "README.md")
    _write_text_file(text_path, "hello\nworld\n")

    @test isfile(text_path)
    @test read(text_path, String) == "hello\nworld\n"
end

@testitem "Shared helpers normalize optional text values" begin
    using Test
    using IDRBindingDataset: _normalize_optional_text

    @test _normalize_optional_text(nothing) === missing
    @test _normalize_optional_text(missing) === missing
    @test _normalize_optional_text("   ") === missing
    @test _normalize_optional_text("  Value  ") == "Value"
    @test _normalize_optional_text("alpha"; uppercase_output = true) == "ALPHA"
    @test _normalize_optional_text("BETA"; lowercase_output = true) == "beta"
    @test _normalize_optional_text("  'Quoted'  "; strip_wrapping_quotes = true) == "Quoted"
    @test _normalize_optional_text(""; missing_value = nothing) === nothing
    @test_throws ArgumentError _normalize_optional_text("x";
        uppercase_output = true,
        lowercase_output = true)
end

@testitem "Shared helpers count labels and render summary lines" begin
    using Test
    using IDRBindingDataset:
                             _count_labels_by_value,
                             _render_label_count_lines

    table = IDRBindingDataset.DataFrames.DataFrame(
        cluster_label = Any["training", missing, "caid", "training", "holdout"]
    )
    counts = _count_labels_by_value(table, :cluster_label; table_name = "table")

    @test counts["training"] == 2
    @test counts["caid"] == 1
    @test counts["missing"] == 1
    @test counts["holdout"] == 1

    lines = _render_label_count_lines(
        counts;
        ordered_labels = ["training", "caid", "missing"],
        prefix = "label"
    )

    @test occursin("- label `training`: `2`", lines)
    @test occursin("- label `caid`: `1`", lines)
    @test occursin("- label `missing`: `1`", lines)
    @test occursin("- label `holdout`: `1`", lines)

    @test_throws ArgumentError _count_labels_by_value(table, :unknown_label;
        table_name = "table")
end
