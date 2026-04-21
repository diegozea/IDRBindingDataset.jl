using TestItems

@testitem "Build DisProt challenge term rules from GO OBO and IDPO OBO-JSON" begin
    using JSON
    using Test
    using IDRBindingDataset: _build_disprot_challenge_term_rules

    tmp_dir = mktempdir()
    raw_data_dir = joinpath(tmp_dir, "raw")
    go_dir = joinpath(raw_data_dir, "go")
    disprot_dir = joinpath(raw_data_dir, "disprot")
    mkpath(go_dir)
    mkpath(disprot_dir)

    go_content = """
format-version: 1.2
ontology: GO

[Term]
id: GO:0005488
name: binding

[Term]
id: GO:0003676
name: nucleic acid binding

[Term]
id: GO:0005515
name: protein binding

[Term]
id: GO:1234567
name: synthetic binding child
is_a: GO:0005488 ! binding

[Term]
id: GO:2345678
name: synthetic nucleic-acid child
is_a: GO:0003676 ! nucleic acid binding

[Term]
id: GO:3456789
name: synthetic protein-binding child
is_a: GO:0005515 ! protein binding

[Term]
id: GO:4567890
name: obsolete binding child
is_obsolete: true
"""

    go_path = joinpath(go_dir, "go-basic.obo.gz")
    open(go_path, "w") do io
        gzip_io = IDRBindingDataset.CodecZlib.GzipCompressorStream(io)
        try
            write(gzip_io, go_content)
        finally
            close(gzip_io)
        end
    end

    idpo_payload = Dict(
        "graphs" => Any[
        Dict(
        "nodes" => Any[
            Dict("id" => "http://purl.obolibrary.org/obo/IDPO_0000002"),
            Dict("id" => "http://purl.obolibrary.org/obo/IDPO_0000003"),
            Dict("id" => "http://purl.obolibrary.org/obo/IDPO_0000004"),
            Dict("id" => "http://purl.obolibrary.org/obo/IDPO_0000010"),
            Dict("id" => "http://purl.obolibrary.org/obo/IDPO_0000011")
        ],
        "edges" => Any[
            Dict(
                "sub" => "http://purl.obolibrary.org/obo/IDPO_0000003",
                "pred" => "is_a",
                "obj" => "http://purl.obolibrary.org/obo/IDPO_0000002"
            ),
            Dict(
                "sub" => "http://purl.obolibrary.org/obo/IDPO_0000004",
                "pred" => "is_a",
                "obj" => "http://purl.obolibrary.org/obo/IDPO_0000002"
            ),
            Dict(
                "sub" => "http://purl.obolibrary.org/obo/IDPO_0000011",
                "pred" => "is_a",
                "obj" => "http://purl.obolibrary.org/obo/IDPO_0000010"
            )
        ]
    )
    ]
    )

    idpo_path = joinpath(disprot_dir, "idpo.json.gz")
    open(idpo_path, "w") do io
        gzip_io = IDRBindingDataset.CodecZlib.GzipCompressorStream(io)
        try
            write(gzip_io, JSON.json(idpo_payload))
        finally
            close(gzip_io)
        end
    end

    rules = _build_disprot_challenge_term_rules(raw_data_dir)

    @test "GO:0005488" in rules.go_binding_terms
    @test "GO:1234567" in rules.go_binding_terms
    @test "GO:4567890" ∉ rules.go_binding_terms

    @test "GO:0003676" in rules.go_binding_nucleic_acid_terms
    @test "GO:2345678" in rules.go_binding_nucleic_acid_terms
    @test "GO:1234567" ∉ rules.go_binding_nucleic_acid_terms

    @test "GO:0005515" in rules.go_binding_protein_terms
    @test "GO:3456789" in rules.go_binding_protein_terms
    @test "GO:1234567" ∉ rules.go_binding_protein_terms

    @test "IDPO:0000002" in rules.idpo_disorder_terms
    @test "IDPO:0000003" in rules.idpo_disorder_terms
    @test "IDPO:0000004" in rules.idpo_disorder_terms
    @test "IDPO:0000011" ∉ rules.idpo_disorder_terms
end

@testitem "Build DisProt challenge term rules validates required input files" begin
    using Test
    using IDRBindingDataset: _build_disprot_challenge_term_rules

    tmp_dir = mktempdir()
    raw_data_dir = joinpath(tmp_dir, "raw")
    mkpath(raw_data_dir)

    @test_throws ArgumentError _build_disprot_challenge_term_rules(raw_data_dir)
end
