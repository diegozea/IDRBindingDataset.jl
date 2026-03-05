using TestItems

@testitem "ELM URL candidate builder adds HTTP fallback for ELM HTTPS URLs" begin
    using Test
    using IDRBindingDataset: _build_elm_url_candidates

    @test _build_elm_url_candidates("https://elm.eu.org/instances.tsv?q=*") == [
        "https://elm.eu.org/instances.tsv?q=*",
        "http://elm.eu.org/instances.tsv?q=*"
    ]
    @test _build_elm_url_candidates("http://elm.eu.org/elms/elms_index.tsv") == [
        "http://elm.eu.org/elms/elms_index.tsv"
    ]
    @test _build_elm_url_candidates("https://example.org/instances.tsv") == [
        "https://example.org/instances.tsv"
    ]
end

@testitem "ELM header metadata parser and integer key parser" begin
    using Test
    using IDRBindingDataset: _parse_elm_header_metadata, _parse_required_header_int

    tmp_dir = mktempdir()
    tsv_path = joinpath(tmp_dir, "instances.tsv")
    write(
        tsv_path,
        """
#NumInstances: 2
#DownloadedFrom: https://elm.eu.org/instances.tsv?q=*
Accession\tELMIdentifier\tPrimary_Acc\tStart\tEnd
ELME000001\tLIG_TEST\tP1\t1\t2
ELME000002\tLIG_TEST\tP2\t3\t4
"""
    )

    metadata = _parse_elm_header_metadata(tsv_path)
    @test metadata["NumInstances"] == "2"
    @test metadata["DownloadedFrom"] == "https://elm.eu.org/instances.tsv?q=*"
    @test _parse_required_header_int(metadata, "NumInstances", "ELM instances") == 2

    missing_key_error = try
        _parse_required_header_int(metadata, "Num_Classes", "ELM index")
        nothing
    catch err
        err
    end
    @test missing_key_error isa ArgumentError

    bad_value_error = try
        _parse_required_header_int(Dict("NumInstances" => "two"), "NumInstances", "ELM instances")
        nothing
    catch err
        err
    end
    @test bad_value_error isa ArgumentError
end

@testitem "ELM TSV validator accepts valid table and rejects row-count mismatch" begin
    using Test
    using IDRBindingDataset: _validate_elm_tsv

    tmp_dir = mktempdir()
    valid_path = joinpath(tmp_dir, "valid_instances.tsv")
    write(
        valid_path,
        """
#NumInstances: 2
Accession\tELMIdentifier\tPrimary_Acc\tStart\tEnd
ELME000001\tLIG_TEST\tP1\t1\t2
ELME000002\tLIG_TEST\tP2\t3\t4
"""
    )

    @test _validate_elm_tsv(
        valid_path,
        ("Accession", "ELMIdentifier", "Primary_Acc", "Start", "End"),
        "ELM instances",
        "NumInstances"
    ) === nothing

    mismatch_path = joinpath(tmp_dir, "mismatch_instances.tsv")
    write(
        mismatch_path,
        """
#NumInstances: 3
Accession\tELMIdentifier\tPrimary_Acc\tStart\tEnd
ELME000001\tLIG_TEST\tP1\t1\t2
ELME000002\tLIG_TEST\tP2\t3\t4
"""
    )

    err = try
        _validate_elm_tsv(
            mismatch_path,
            ("Accession", "ELMIdentifier", "Primary_Acc", "Start", "End"),
            "ELM instances",
            "NumInstances"
        )
        nothing
    catch caught
        caught
    end

    @test err isa ArgumentError
    @test occursin("row count mismatch", sprint(showerror, err))
end

@testitem "ELM download README writer includes source URLs and row counts" begin
    using Test
    using IDRBindingDataset: _write_elm_download_readme

    tmp_dir = mktempdir()
    instances_path = joinpath(tmp_dir, "instances.tsv")
    index_path = joinpath(tmp_dir, "elms_index.tsv")
    readme_path = joinpath(tmp_dir, "README.md")

    write(instances_path, "instances")
    write(index_path, "index")

    _write_elm_download_readme(
        readme_path;
        instances_file_report = (
            requested_url = "https://elm.eu.org/instances.tsv?q=*",
            resolved_url = "http://elm.eu.org/instances.tsv?q=*"
        ),
        classes_index_file_report = (
            requested_url = "http://elm.eu.org/elms/elms_index.tsv",
            resolved_url = "http://elm.eu.org/elms/elms_index.tsv"
        ),
        instances_path = instances_path,
        classes_index_path = index_path,
        instances_rows = 10,
        classes_index_rows = 4,
        instances_header_rows = 10,
        classes_index_header_rows = 4,
        downloaded_at_utc = "2026-03-05T11:00:00Z"
    )

    content = read(readme_path, String)
    @test occursin("# ELM Scratch Download", content)
    @test occursin("Download timestamp (UTC): `2026-03-05T11:00:00Z`", content)
    @test occursin("Requested source URL", content)
    @test occursin("Downloaded from URL", content)
    @test occursin("Header `#NumInstances`: `10`", content)
    @test occursin("Header `#Num_Classes`: `4`", content)
end
