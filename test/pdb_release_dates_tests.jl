using TestItems

@testitem "PDB release dates parse API response and fill missing IDs" begin
    using Dates
    using HTTP
    using JSON
    using Test
    using IDRBindingDataset: fetch_pdb_initial_release_dates

    request_count = Ref(0)
    requested_batches = Vector{Vector{String}}()

    function mock_post(url::AbstractString, _, body::AbstractString)
        request_count[] += 1
        @test url == "https://data.rcsb.org/graphql"

        payload = JSON.parse(body)
        batch_ids = String.(payload["variables"]["ids"])
        push!(requested_batches, batch_ids)

        response_payload = Dict(
            "data" => Dict(
            "entries" => Any[
            Dict(
                "rcsb_id" => "1ABC",
                "rcsb_accession_info" => Dict(
                    "initial_release_date" => "2001-02-03T00:00:00+00:00"
                )
            ),
            Dict(
                "rcsb_id" => "2XYZ",
                "rcsb_accession_info" => Dict("initial_release_date" => nothing)
            )
        ]
        )
        )

        return HTTP.Response(
            200,
            ["Content-Type" => "application/json"],
            JSON.json(response_payload)
        )
    end

    release_dates = fetch_pdb_initial_release_dates(
        ["1abc", "2xyz", "9zzz"];
        batch_size = 300,
        max_retries = 2,
        base_sleep = 0.0,
        post_request = mock_post,
        sleeper = _ -> nothing
    )

    @test request_count[] == 1
    @test requested_batches == [["1ABC", "2XYZ", "9ZZZ"]]
    @test release_dates["1ABC"] == Date(2001, 2, 3)
    @test release_dates["2XYZ"] === missing
    @test release_dates["9ZZZ"] === missing
end

@testitem "PDB release dates respect batch_size=300 partitioning" begin
    using Dates
    using HTTP
    using JSON
    using Test
    using IDRBindingDataset: fetch_pdb_initial_release_dates

    batch_sizes = Int[]

    function mock_post(_, _, body::AbstractString)
        payload = JSON.parse(body)
        batch_ids = String.(payload["variables"]["ids"])
        push!(batch_sizes, length(batch_ids))

        entries = [Dict(
                       "rcsb_id" => pdb_id,
                       "rcsb_accession_info" => Dict("initial_release_date" => "1990-01-01")
                   )
                   for pdb_id in batch_ids]

        response_payload = Dict("data" => Dict("entries" => entries))
        return HTTP.Response(
            200,
            ["Content-Type" => "application/json"],
            JSON.json(response_payload)
        )
    end

    pdb_ids = [string(lpad(i, 4, '0')) for i in 1:650]
    release_dates = fetch_pdb_initial_release_dates(
        pdb_ids;
        batch_size = 300,
        max_retries = 1,
        base_sleep = 0.0,
        post_request = mock_post,
        sleeper = _ -> nothing
    )

    @test batch_sizes == [300, 300, 50]
    @test length(release_dates) == 650
    @test release_dates["0001"] == Date(1990, 1, 1)
    @test release_dates["0650"] == Date(1990, 1, 1)
end

@testitem "PDB release dates retry HTTP 429 responses" begin
    using Dates
    using HTTP
    using JSON
    using Test
    using IDRBindingDataset: fetch_pdb_initial_release_dates

    attempts = Ref(0)
    sleep_calls = Float64[]

    function mock_post(_, _, _)
        attempts[] += 1
        if attempts[] == 1
            return HTTP.Response(429, ["Retry-After" => "0.25"], "")
        end

        response_payload = Dict(
            "data" => Dict(
            "entries" => Any[
            Dict(
            "rcsb_id" => "1ABC",
            "rcsb_accession_info" => Dict("initial_release_date" => "1999-12-31")
        )
        ]
        )
        )
        return HTTP.Response(
            200,
            ["Content-Type" => "application/json"],
            JSON.json(response_payload)
        )
    end

    release_dates = fetch_pdb_initial_release_dates(
        ["1ABC"];
        batch_size = 300,
        max_retries = 3,
        base_sleep = 0.0,
        post_request = mock_post,
        sleeper = seconds -> push!(sleep_calls, Float64(seconds))
    )

    @test attempts[] == 2
    @test length(sleep_calls) == 1
    @test sleep_calls[1] ≈ 0.25
    @test release_dates["1ABC"] == Date(1999, 12, 31)
end

@testitem "PDB release dates throw on GraphQL errors" begin
    using HTTP
    using JSON
    using Test
    using IDRBindingDataset: fetch_pdb_initial_release_dates

    function mock_post(_, _, _)
        response_payload = Dict(
            "errors" => Any[
            Dict("message" => "Cannot query field initial_release_date")
        ]
        )
        return HTTP.Response(
            200,
            ["Content-Type" => "application/json"],
            JSON.json(response_payload)
        )
    end

    err = try
        fetch_pdb_initial_release_dates(
            ["1ABC"];
            batch_size = 300,
            max_retries = 1,
            base_sleep = 0.0,
            post_request = mock_post,
            sleeper = _ -> nothing
        )
        nothing
    catch caught
        caught
    end

    @test err isa ErrorException
    @test occursin("RCSB GraphQL error", sprint(showerror, err))
end

@testitem "PDB release dates throw after repeated HTTP failures" begin
    using HTTP
    using Test
    using IDRBindingDataset: fetch_pdb_initial_release_dates

    attempts = Ref(0)
    sleep_calls = Float64[]

    function always_500(_, _, _)
        attempts[] += 1
        return HTTP.Response(500, ["Content-Type" => "application/json"], "")
    end

    err = try
        fetch_pdb_initial_release_dates(
            ["1ABC"];
            batch_size = 300,
            max_retries = 3,
            base_sleep = 0.1,
            post_request = always_500,
            sleeper = seconds -> push!(sleep_calls, Float64(seconds))
        )
        nothing
    catch caught
        caught
    end

    @test err isa ErrorException
    @test occursin("last status = 500", sprint(showerror, err))
    @test attempts[] == 3
    @test sleep_calls == [0.1, 0.2]
end

@testitem "PDB code extractor normalizes and captures canonical IDs" begin
    using Test
    using IDRBindingDataset: extract_pdb_codes

    @test extract_pdb_codes("1DMG:A") == ["1DMG"]
    @test extract_pdb_codes("1ABC_A") == ["1ABC"]
    @test extract_pdb_codes("pdb:2phe") == ["2PHE"]
    @test extract_pdb_codes("4TTB;PDB:4TTC") == ["4TTB", "4TTC"]
    @test extract_pdb_codes("14782") == String[]
    @test extract_pdb_codes("7Q1YY") == String[]
end

@testitem "PDB holding statuses parse removed and not-found IDs" begin
    using HTTP
    using JSON
    using Test
    using IDRBindingDataset: fetch_pdb_holding_statuses

    seen_urls = String[]

    function mock_get(url::AbstractString)
        push!(seen_urls, String(url))

        if endswith(url, "/1HHB")
            body = JSON.json(
                Dict(
                "rcsb_repository_holdings_combined" => Dict(
                "status" => "REMOVED",
                "status_code" => "OBS",
                "id_code_replaced_by_latest" => "4HHB"
            )
            )
            )
            return HTTP.Response(200, ["Content-Type" => "application/json"], body)
        end

        if endswith(url, "/12TW")
            body = JSON.json(Dict("status" => 404, "message" => "No data found"))
            return HTTP.Response(404, ["Content-Type" => "application/json"], body)
        end

        error("Unexpected URL in mock: $(url)")
    end

    statuses = fetch_pdb_holding_statuses(
        ["1hhb", "12tw"];
        max_retries = 1,
        base_sleep = 0.0,
        get_request = mock_get,
        sleeper = _ -> nothing
    )

    @test seen_urls == [
        "https://data.rcsb.org/rest/v1/holdings/status/1HHB",
        "https://data.rcsb.org/rest/v1/holdings/status/12TW"
    ]

    @test statuses["1HHB"].http_status == 200
    @test statuses["1HHB"].status == "REMOVED"
    @test statuses["1HHB"].status_code == "OBS"
    @test statuses["1HHB"].replacement_pdb == "4HHB"

    @test statuses["12TW"].http_status == 404
    @test statuses["12TW"].status === missing
    @test statuses["12TW"].status_code === missing
    @test statuses["12TW"].replacement_pdb === missing
end

@testitem "PDB holding statuses retry HTTP 429 responses" begin
    using HTTP
    using JSON
    using Test
    using IDRBindingDataset: fetch_pdb_holding_statuses

    attempts = Ref(0)
    sleep_calls = Float64[]

    function mock_get(_)
        attempts[] += 1
        if attempts[] == 1
            return HTTP.Response(429, ["Retry-After" => "0.5"], "")
        end

        body = JSON.json(
            Dict(
            "rcsb_repository_holdings_combined" => Dict(
            "status" => "CURRENT",
            "status_code" => "REL"
        )
        )
        )
        return HTTP.Response(200, ["Content-Type" => "application/json"], body)
    end

    statuses = fetch_pdb_holding_statuses(
        ["4HHB"];
        max_retries = 3,
        base_sleep = 0.0,
        get_request = mock_get,
        sleeper = seconds -> push!(sleep_calls, Float64(seconds))
    )

    @test attempts[] == 2
    @test length(sleep_calls) == 1
    @test sleep_calls[1] ≈ 0.5
    @test statuses["4HHB"].status == "CURRENT"
    @test statuses["4HHB"].status_code == "REL"
end

@testitem "PDB holding statuses throw after repeated HTTP failures" begin
    using HTTP
    using Test
    using IDRBindingDataset: fetch_pdb_holding_statuses

    attempts = Ref(0)
    sleep_calls = Float64[]

    function always_500(_)
        attempts[] += 1
        return HTTP.Response(500, ["Content-Type" => "application/json"], "")
    end

    err = try
        fetch_pdb_holding_statuses(
            ["1ABC"];
            max_retries = 3,
            base_sleep = 0.1,
            get_request = always_500,
            sleeper = seconds -> push!(sleep_calls, Float64(seconds))
        )
        nothing
    catch caught
        caught
    end

    @test err isa ErrorException
    @test occursin("RCSB holdings status request failed for 1ABC", sprint(showerror, err))
    @test occursin("last status = 500", sprint(showerror, err))
    @test attempts[] == 3
    @test sleep_calls == [0.1, 0.2]
end

@testitem "PDB lookup classifier distinguishes current obsolete and invalid" begin
    using Dates
    using Test
    using IDRBindingDataset: classify_pdb_lookup_status

    current_status = (
        http_status = 200,
        status = "CURRENT",
        status_code = "REL",
        replacement_pdb = missing
    )
    removed_status = (
        http_status = 200,
        status = "REMOVED",
        status_code = "OBS",
        replacement_pdb = "4HHB"
    )
    missing_status = (
        http_status = 404,
        status = missing,
        status_code = missing,
        replacement_pdb = missing
    )
    unresolved_status = (
        http_status = 200,
        status = missing,
        status_code = missing,
        replacement_pdb = missing
    )

    @test classify_pdb_lookup_status(Date(2024, 1, 1), current_status) == "current"
    @test classify_pdb_lookup_status(missing, removed_status) == "obsolete"
    @test classify_pdb_lookup_status(missing, missing_status) == "invalid_or_not_found"
    @test classify_pdb_lookup_status(missing, current_status) ==
          "current_without_release_date"
    @test classify_pdb_lookup_status(missing, unresolved_status) == "unresolved"
end
