using Dates
using HTTP
using JSON

# ---------------------------------------------------------------------------
# RCSB endpoints, query templates, and request defaults.
# ---------------------------------------------------------------------------
# External API references (official RCSB docs) used by this file:
#
# 1) Data API tutorial: https://data.rcsb.org/
#    - GraphQL endpoint is `https://data.rcsb.org/graphql`.
#    - GraphQL endpoint handles HTTP GET and POST methods.
#    - GraphQL POST body uses `application/json`, must include `query`, and may
#      include `variables`.
#    - The tutorial's "Error handling in GraphQL" section states that the server
#      responds with `200 OK` and returns error details in the response body.
#    - The REST section documents `404` when data/endpoint is not found.
#
# 2) Migration guide: https://data.rcsb.org/migration-guide.html
#    - "Release Status" endpoints:
#      `https://data.rcsb.org/rest/v1/holdings/status/{pdb_id}`
#      `https://data.rcsb.org/rest/v1/holdings/status?ids=...`
#    - Notes section: "GraphQL endpoint allows accessing only current entries."
#
# 3) REST OpenAPI specification:
#    https://data.rcsb.org/redoc/rcsb-restful-api-docs.json
#    - Path `/holdings/status/{entry_id}`:
#      summary "Get the status and status code of a structure.",
#      responses include `200` and `404`.
#    - Path `/holdings/status`:
#      summary "Get the status and status code of a list of entry IDs."
#    - `CombinedEntry` includes:
#      `rcsb_repository_holdings_combined` and
#      `rcsb_repository_holdings_combined_entry_container_identifiers`.
#    - `RcsbRepositoryHoldingsCombined` defines:
#      `status`, `status_code`, and optional `id_code_replaced_by_latest`.
#      The current-entry example (`/holdings/status/4HHB`) includes only
#      `status` and `status_code`; removed-entry payloads can include
#      `id_code_replaced_by_latest` (for example, `/holdings/status/4ACQ`).
#
# 4) Entry JSON schema:
#    https://data.rcsb.org/rest/v1/schema/entry
#    - `rcsb_accession_info.initial_release_date` is documented as
#      "The entry initial release date.".
#
# This module reads:
# - GraphQL field `entries(...).rcsb_accession_info.initial_release_date`.
# - Holdings fields under `rcsb_repository_holdings_combined`.
const RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
const RCSB_HOLDINGS_STATUS_URL = "https://data.rcsb.org/rest/v1/holdings/status"
const RCSB_RELEASE_DATE_QUERY = """
query (\$ids: [String!]!) {
  entries(entry_ids: \$ids) {
    rcsb_id
    rcsb_accession_info { initial_release_date }
  }
}
"""
const RCSB_DEFAULT_BATCH_SIZE = 300
const RCSB_DEFAULT_MAX_RETRIES = 5
const RCSB_DEFAULT_BASE_SLEEP = 0.5
const RCSB_HOLDINGS_NOT_FOUND_HTTP_STATUS = 404

# We are considering the canonical PDB ID format to be a 4-character string starting with 
# a digit, as per RCSB conventions. We do not support the new pdb_* identifiers in this 
# implementation -- for the moment.
const PDB_CODE_CAPTURE_REGEX = r"(?i)(?<![A-Z0-9])(?:PDB:)?([0-9][A-Z0-9]{3})(?![A-Z0-9])"
const PDB_CANONICAL_CODE_REGEX = r"^[0-9][A-Z0-9]{3}$"

# ---------------------------------------------------------------------------
# Shared type aliases for release-date and holdings-status records.
# ---------------------------------------------------------------------------
const PDBHoldingStatus = NamedTuple{
    (:http_status, :status, :status_code, :replacement_pdb),
    Tuple{Int, Union{Missing, String}, Union{Missing, String}, Union{Missing, String}}
}
const PDBReleaseDate = Union{Missing, Date}
const PDBReleaseDateDict = Dict{String, PDBReleaseDate}

# ---------------------------------------------------------------------------
# PDB identifier normalization and extraction helpers.
# ---------------------------------------------------------------------------
"""
    _normalize_pdb_ids(pdb_ids)

Normalize PDB identifiers by trimming whitespace and converting to uppercase.
"""
function _normalize_pdb_ids(
        pdb_ids::AbstractVector{<:AbstractString}
)::Vector{String}
    return uppercase.(strip.(String.(pdb_ids)))
end

"""
    extract_pdb_codes(reference)

Extract canonical 4-character PDB identifiers from one raw `reference` string.

Accepted patterns include raw IDs (`1ABC`), IDs with chain suffixes
(`1ABC:A`, `1ABC_A`), prefixed IDs (`PDB:1ABC`), and concatenated strings
where IDs are token-delimited (for example `1ABC;PDB:2DEF`).
"""
function extract_pdb_codes(reference::AbstractString)::Vector{String}
    extracted = String[]
    for match in eachmatch(PDB_CODE_CAPTURE_REGEX, reference)
        captured_id = match.captures[1]
        captured_id === nothing && continue
        push!(extracted, uppercase(String(captured_id)))
    end
    # `unique` preserves first-seen order, keeping output deterministic.
    return unique(extracted)
end

# ---------------------------------------------------------------------------
# Shared retry/backoff logic for HTTP calls.
# ---------------------------------------------------------------------------
"""
    _retry_after_seconds(response, base_sleep, attempt)

Return the wait duration for an HTTP 429 response. If `Retry-After` is not
present or cannot be parsed, this falls back to exponential backoff.
"""
function _retry_after_seconds(
        response::HTTP.Response,
        base_sleep::Real,
        attempt::Integer
)::Float64
    retry_after = HTTP.header(response, "Retry-After")
    parsed = retry_after === nothing ? nothing : tryparse(Float64, retry_after)
    if parsed === nothing || parsed < 0
        return base_sleep * 2.0^(attempt - 1)
    end
    return parsed
end

"""
    _retry_sleep_seconds_for_status(response, base_sleep, attempt)

Return retry sleep duration for one non-success HTTP response.
"""
function _retry_sleep_seconds_for_status(
        response::HTTP.Response,
        base_sleep::Real,
        attempt::Integer
)::Float64
    # Respect `Retry-After` for 429 responses; otherwise use linear wait.
    if response.status == 429
        return _retry_after_seconds(response, base_sleep, attempt)
    end
    return base_sleep * attempt
end

"""
    _retry_sleep_seconds_for_exception(base_sleep, attempt)

Return retry sleep duration for one transport/connection exception.
"""
function _retry_sleep_seconds_for_exception(
        base_sleep::Real,
        attempt::Integer
)::Float64
    # Connection-level exceptions use exponential backoff.
    return base_sleep * 2.0^(attempt - 1)
end

"""
    _retry_last_status_message(last_status)

Render the terminal HTTP-status message for retry-failure errors.
"""
function _retry_last_status_message(last_status)::String
    return last_status === nothing ? "last status = none" : "last status = $(last_status)"
end

"""
    _request_with_retries(
        request_fn;
        max_retries,
        base_sleep,
        sleeper,
        is_success_response,
        failure_prefix
    )

Run one HTTP request with bounded retries and shared backoff logic.

- Retries on non-success responses and exceptions.
- Uses `Retry-After` when HTTP status is 429.
- Emits consistent terminal errors containing last status and last exception.
"""
function _request_with_retries(
        request_fn::Function;
        max_retries::Integer,
        base_sleep::Real,
        sleeper::Function,
        is_success_response::Function,
        failure_prefix::AbstractString
)::HTTP.Response
    last_error = nothing
    last_status = nothing

    for attempt in 1:max_retries
        try
            response = request_fn()
            if is_success_response(response)
                return response
            end

            last_status = response.status
            attempt == max_retries && break
            sleeper(_retry_sleep_seconds_for_status(response, base_sleep, attempt))
        catch err
            last_error = err
            attempt == max_retries && break
            sleeper(_retry_sleep_seconds_for_exception(base_sleep, attempt))
        end
    end

    status_message = _retry_last_status_message(last_status)
    if last_error !== nothing
        error(
            "$(failure_prefix) after $(max_retries) attempts ($(status_message)): $(sprint(showerror, last_error))",
        )
    end
    error("$(failure_prefix) after $(max_retries) attempts ($(status_message)).")
end

"""
    _post_rcsb_graphql_with_retries(body; max_retries, base_sleep, post_request, sleeper)

POST one GraphQL payload to RCSB with bounded retries.
"""
function _post_rcsb_graphql_with_retries(
        body::AbstractString;
        max_retries::Integer,
        base_sleep::Real,
        post_request::Function,
        sleeper::Function
)::HTTP.Response
    headers = ["Content-Type" => "application/json"]
    return _request_with_retries(
        () -> post_request(RCSB_GRAPHQL_URL, headers, body);
        max_retries = max_retries,
        base_sleep = base_sleep,
        sleeper = sleeper,
        is_success_response = response -> response.status == 200,
        failure_prefix = "RCSB GraphQL request failed"
    )
end

"""
    _default_rcsb_status_get_request(url)

Run one holdings-status GET request without throwing on non-2xx statuses.
"""
function _default_rcsb_status_get_request(url::AbstractString)::HTTP.Response
    return HTTP.get(url; status_exception = false)
end

"""
    _get_rcsb_holding_status_with_retries(pdb_id; max_retries, base_sleep, get_request, sleeper)

Fetch one holdings-status response for `pdb_id`, retrying on transient failures
and HTTP 429 responses.
"""
function _get_rcsb_holding_status_with_retries(
        pdb_id::AbstractString;
        max_retries::Integer,
        base_sleep::Real,
        get_request::Function,
        sleeper::Function
)::HTTP.Response
    request_url = string(RCSB_HOLDINGS_STATUS_URL, "/", pdb_id)
    return _request_with_retries(
        () -> get_request(request_url);
        max_retries = max_retries,
        base_sleep = base_sleep,
        sleeper = sleeper,
        # `/holdings/status/{entry_id}` documents `200` and `404` responses.
        is_success_response = response -> response.status == 200 ||
                                          response.status ==
                                          RCSB_HOLDINGS_NOT_FOUND_HTTP_STATUS,
        failure_prefix = "RCSB holdings status request failed for $(pdb_id)"
    )
end

# ---------------------------------------------------------------------------
# Release-date parsing and GraphQL response conversion.
# ---------------------------------------------------------------------------
"""
    _parse_release_date(raw_date, pdb_id)

Parse an RCSB `initial_release_date` value into `Date`. Returns `missing` when
the API explicitly reports a missing date.
"""
function _parse_release_date(
        raw_date,
        pdb_id::AbstractString
)::Union{Missing, Date}
    raw_date === nothing && return missing

    date_text = String(raw_date)
    length(date_text) >= 10 ||
        error("Invalid release date returned for PDB $(pdb_id): $(date_text)")

    try
        return Date(date_text[1:10])
    catch err
        if err isa ArgumentError
            error("Unparsable release date returned for PDB $(pdb_id): $(date_text)")
        end
        rethrow(err)
    end
end

"""
    _extract_release_dates_from_response(response_body, batch_ids)

Parse one GraphQL response and return a dictionary from PDB ID to release date.
Missing or unknown IDs are represented with `missing`.
"""
function _extract_release_dates_from_response(
        response_body::AbstractString,
        batch_ids::AbstractVector{<:AbstractString}
)::PDBReleaseDateDict
    parsed = JSON.parse(response_body)

    if haskey(parsed, "errors")
        first_error = first(parsed["errors"])
        message = get(first_error, "message", string(first_error))
        error("RCSB GraphQL error: $(message)")
    end

    entries = get(get(parsed, "data", Dict{String, Any}()), "entries", Any[])
    output = PDBReleaseDateDict()

    for entry in entries
        # GraphQL can include explicit `nothing` entries for unresolved IDs.
        entry === nothing && continue

        pdb_id = uppercase(String(get(entry, "rcsb_id", "")))
        isempty(pdb_id) && continue

        accession_info = get(entry, "rcsb_accession_info", nothing)
        raw_date = accession_info === nothing ?
                   nothing :
                   get(accession_info, "initial_release_date", nothing)
        output[pdb_id] = _parse_release_date(raw_date, pdb_id)
    end

    # Ensure every requested ID exists in the output map.
    for pdb_id in batch_ids
        get!(output, String(pdb_id), missing)
    end

    return output
end

# ---------------------------------------------------------------------------
# Public API: initial release dates via GraphQL.
# ---------------------------------------------------------------------------
"""
    fetch_pdb_initial_release_dates(
        pdb_ids;
        batch_size = 300,
        max_retries = 5,
        base_sleep = 0.5
    )

Fetch initial PDB release dates from the RCSB GraphQL API for `pdb_ids`.

Returns `Dict{String,Union{Missing,Date}}` keyed by normalized uppercase PDB ID.
Any unknown ID is included with value `missing`.
"""
function fetch_pdb_initial_release_dates(
        pdb_ids::AbstractVector{<:AbstractString};
        batch_size::Integer = RCSB_DEFAULT_BATCH_SIZE,
        max_retries::Integer = RCSB_DEFAULT_MAX_RETRIES,
        base_sleep::Real = RCSB_DEFAULT_BASE_SLEEP,
        post_request::Function = HTTP.post,
        sleeper::Function = sleep
)::PDBReleaseDateDict
    batch_size > 0 || throw(ArgumentError("`batch_size` must be greater than zero."))
    max_retries > 0 || throw(ArgumentError("`max_retries` must be greater than zero."))
    base_sleep >= 0 || throw(ArgumentError("`base_sleep` must be non-negative."))

    normalized_ids = _normalize_pdb_ids(pdb_ids)
    output = PDBReleaseDateDict()
    sizehint!(output, length(normalized_ids))

    # Preserve explicit empty inputs as `missing` to keep mapping stable.
    for pdb_id in normalized_ids
        isempty(pdb_id) && (output[pdb_id] = missing)
    end

    query_ids = unique(filter(!isempty, normalized_ids))
    isempty(query_ids) && return output

    # Request IDs in fixed-size chunks (project default is 300 IDs per request).
    for batch in Iterators.partition(query_ids, batch_size)
        batch_ids = collect(batch)
        body = JSON.json(Dict(
            "query" => RCSB_RELEASE_DATE_QUERY,
            "variables" => Dict("ids" => batch_ids)
        ))

        response = _post_rcsb_graphql_with_retries(
            body;
            max_retries = max_retries,
            base_sleep = base_sleep,
            post_request = post_request,
            sleeper = sleeper
        )

        batch_output = _extract_release_dates_from_response(
            String(response.body),
            batch_ids
        )
        merge!(output, batch_output)
    end

    return output
end

# ---------------------------------------------------------------------------
# Holdings-status parsing and fallback lookup API.
# ---------------------------------------------------------------------------
"""
    _normalize_optional_string(value; uppercase_output = false)

Normalize an optional JSON field into `Union{Missing,String}`.
"""
function _normalize_optional_string(
        value;
        uppercase_output::Bool = false
)::Union{Missing, String}
    ismissing(value) && return uppercase_output ? "MISSING" : "missing"

    normalized = _normalize_optional_text(
        value;
        uppercase_output = uppercase_output,
        missing_value = nothing
    )
    normalized === nothing && return missing
    normalized
end

"""
    _parse_holding_status_response(response)

Parse one RCSB holdings status response into a typed status tuple.
"""
function _parse_holding_status_response(
        response::HTTP.Response
)::PDBHoldingStatus
    # Data API docs describe `404` for REST missing/not-found responses.
    if response.status == RCSB_HOLDINGS_NOT_FOUND_HTTP_STATUS
        return (
            http_status = RCSB_HOLDINGS_NOT_FOUND_HTTP_STATUS,
            status = missing,
            status_code = missing,
            replacement_pdb = missing
        )
    end

    response.status == 200 ||
        error("Unexpected holdings status HTTP code: $(response.status)")

    parsed = JSON.parse(String(response.body))
    holdings = get(parsed, "rcsb_repository_holdings_combined", Dict{String, Any}())

    # Fields documented in OpenAPI schema `RcsbRepositoryHoldingsCombined`:
    # - `status`
    # - `status_code`
    # - optional `id_code_replaced_by_latest`
    return (
        http_status = 200,
        status = _normalize_optional_string(get(holdings, "status", nothing);
            uppercase_output = true),
        status_code = _normalize_optional_string(get(holdings, "status_code", nothing);
            uppercase_output = true),
        replacement_pdb = _normalize_optional_string(
            get(holdings, "id_code_replaced_by_latest", nothing);
            uppercase_output = true
        )
    )
end

"""
    fetch_pdb_holding_statuses(
        pdb_ids;
        max_retries = 5,
        base_sleep = 0.5
    )

Fetch RCSB holdings status records for canonical 4-character PDB IDs.

Returns `Dict{String,PDBHoldingStatus}` keyed by uppercase PDB ID.
"""
function fetch_pdb_holding_statuses(
        pdb_ids::AbstractVector{<:AbstractString};
        max_retries::Integer = RCSB_DEFAULT_MAX_RETRIES,
        base_sleep::Real = RCSB_DEFAULT_BASE_SLEEP,
        get_request::Function = _default_rcsb_status_get_request,
        sleeper::Function = sleep
)::Dict{String, PDBHoldingStatus}
    max_retries > 0 || throw(ArgumentError("`max_retries` must be greater than zero."))
    base_sleep >= 0 || throw(ArgumentError("`base_sleep` must be non-negative."))

    normalized_ids = unique(filter(!isempty, _normalize_pdb_ids(pdb_ids)))
    output = Dict{String, PDBHoldingStatus}()
    sizehint!(output, length(normalized_ids))

    for pdb_id in normalized_ids
        # Holdings endpoint expects canonical 4-character PDB IDs.
        occursin(PDB_CANONICAL_CODE_REGEX, pdb_id) ||
            throw(ArgumentError("Holding status lookups require canonical 4-character PDB IDs; got '$(pdb_id)'."))

        # Holdings is the fallback only for unresolved release-date lookups.
        response = _get_rcsb_holding_status_with_retries(
            pdb_id;
            max_retries = max_retries,
            base_sleep = base_sleep,
            get_request = get_request,
            sleeper = sleeper
        )
        output[pdb_id] = _parse_holding_status_response(response)
    end

    return output
end

# ---------------------------------------------------------------------------
# Final classification helper used by downstream dataset builders.
# ---------------------------------------------------------------------------
"""
    classify_pdb_lookup_status(release_date, holding_status)

Classify one PDB lookup result into one of:
- `current`
- `obsolete`
- `invalid_or_not_found`
- `current_without_release_date`
- `unresolved`
"""
function classify_pdb_lookup_status(
        release_date::Union{Missing, Date},
        holding_status::NamedTuple{(:http_status, :status, :status_code, :replacement_pdb)}
)::String
    !ismissing(release_date) && return "current"

    holding_status.http_status == RCSB_HOLDINGS_NOT_FOUND_HTTP_STATUS &&
        return "invalid_or_not_found"

    status = holding_status.status
    isequal(status, "REMOVED") && return "obsolete"
    isequal(status, "CURRENT") && return "current_without_release_date"

    return "unresolved"
end
