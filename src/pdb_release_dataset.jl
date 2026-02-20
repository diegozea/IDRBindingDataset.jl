using CodecZlib
using CSV
using DataFrames

# ---------------------------------------------------------------------------
# Input/output file locations and table-shape constants.
# ---------------------------------------------------------------------------
const SIFTS_UNIPROT_RELATIVE_PATH = joinpath("sifts", "pdb_chain_uniprot.csv.gz")
const PDB_LOOKUP_FAILURE_REASON_INVALID_FORMAT = "invalid_format_no_pdb_code"
const CleanedPDBReferenceRow = NamedTuple{
    (:acc, :pdb, :pdb_reference), Tuple{String, String, String}
}

"""
    _read_sifts_uniprot_df(sifts_uniprot_path)

Read the SIFTS UniProt summary file and return it as a `DataFrame`.
"""
function _read_sifts_uniprot_df(sifts_uniprot_path::AbstractString)::DataFrame
    sifts = MIToS.SIFTS.read_file(sifts_uniprot_path, MIToS.SIFTS.SIFTSCSV)
    return DataFrame(sifts.table, sifts.colnames)
end

# ---------------------------------------------------------------------------
# Base DisProt/SIFTS mapping tables.
# ---------------------------------------------------------------------------
"""
    collect_disprot_region_level_pdb_refs(disprot_entries)

Collect region-level DisProt PDB cross references into a tidy table with
columns:
- `acc`
- `pdb`
- `pdb_reference` (raw cross-reference value)

PDB IDs are cleaned using `extract_pdb_codes`.
"""
function collect_disprot_region_level_pdb_refs(disprot_entries)::DataFrame
    rows = DataFrame(
        acc = String[],
        pdb = String[],
        pdb_reference = String[]
    )

    for protein in disprot_entries
        acc = String(get(protein, "acc", ""))
        isempty(acc) && continue

        for region in get(protein, "regions", Any[])
            for ref in get(region, "cross_refs", Any[])
                db = uppercase(strip(String(get(ref, "db", ""))))
                db == "PDB" || continue

                raw_pdb_reference = String(get(ref, "id", ""))
                for pdb in extract_pdb_codes(raw_pdb_reference)
                    push!(rows, (
                        acc = acc,
                        pdb = pdb,
                        pdb_reference = raw_pdb_reference
                    ))
                end
            end
        end
    end

    sort!(rows, [:acc, :pdb, :pdb_reference])
    return rows
end

"""
    build_unique_uniprot_pdb_map(region_level_pdb_refs)

Collapse region-level references into one row per unique `(acc, pdb)` pair.
"""
function build_unique_uniprot_pdb_map(
        region_level_pdb_refs::DataFrame
)::DataFrame
    rows = unique(select(region_level_pdb_refs, :acc, :pdb), [:acc, :pdb])
    sort!(rows, [:acc, :pdb])
    return rows
end

"""
    build_sifts_uniprot_pdb_map_for_disprot(sifts_uniprot_df, disprot_uniprot_ids)

Build a unique SIFTS UniProt-PDB map restricted to DisProt UniProt IDs.
"""
function build_sifts_uniprot_pdb_map_for_disprot(
        sifts_uniprot_df::DataFrame,
        disprot_uniprot_ids
)::DataFrame
    allowed_ids = Set(String.(collect(disprot_uniprot_ids)))

    filtered = subset(
        sifts_uniprot_df,
        :SP_PRIMARY => ByRow(x -> in(String(x), allowed_ids))
    )

    selected = select(
        filtered,
        :SP_PRIMARY => ByRow(x -> String(x)) => :acc,
        :PDB => ByRow(x -> uppercase(strip(String(x)))) => :pdb
    )

    rows = unique(selected, [:acc, :pdb])
    sort!(rows, [:acc, :pdb])
    return rows
end

"""
    build_disprot_sifts_pdb_comparison(disprot_uniprot_pdb_map, sifts_uniprot_pdb_map)

Outer-join DisProt and SIFTS UniProt-PDB maps and add `mapping_status` labels:
- `both`
- `only_disprot`
- `only_sifts`
"""
function build_disprot_sifts_pdb_comparison(
        disprot_uniprot_pdb_map::DataFrame,
        sifts_uniprot_pdb_map::DataFrame
)::DataFrame
    comparison = outerjoin(
        disprot_uniprot_pdb_map,
        sifts_uniprot_pdb_map,
        on = [:acc, :pdb],
        source = :mapping_status
    )

    replace!(
        comparison.mapping_status,
        "left_only" => "only_disprot",
        "right_only" => "only_sifts"
    )
    sort!(comparison, [:mapping_status, :acc, :pdb])
    return comparison
end

# ---------------------------------------------------------------------------
# Helpers to expand and validate raw merged PDB references.
# ---------------------------------------------------------------------------
"""
    _expanded_reference_rows(acc, pdb_reference)

Return cleaned PDB rows extracted from one raw `(acc, pdb_reference)` pair.
"""
function _expanded_reference_rows(
        acc::AbstractString,
        pdb_reference::AbstractString
)::Vector{CleanedPDBReferenceRow}
    reference = String(pdb_reference)
    return [(acc = String(acc), pdb = pdb_code, pdb_reference = reference)
            for pdb_code in extract_pdb_codes(reference)]
end

"""
    _expand_merged_pdb_references(merged_pdb_references)

Expand raw merged references to one row per extracted canonical PDB ID using
DataFrames row-wise transformations.
"""
function _expand_merged_pdb_references(
        merged_pdb_references::DataFrame
)::DataFrame
    # Create one nested vector of cleaned rows per raw `(acc, pdb_reference)` row.
    expanded = transform(
        merged_pdb_references,
        [:acc,
            :pdb] => ByRow((acc,
            pdb_reference) -> _expanded_reference_rows(String(acc), String(pdb_reference))) => :expanded_rows
    )
    # `flatten` + `AsTable` turns nested vectors of named tuples into flat columns.
    cleaned = select(flatten(expanded, :expanded_rows), :expanded_rows => AsTable)
    sort!(cleaned, [:acc, :pdb, :pdb_reference])
    return cleaned
end

"""
    _build_invalid_merged_pdb_references(merged_pdb_references, cleaned_merged_pdb_rows)

Build a table of raw references that contain no extractable canonical PDB ID.
This reuses already-expanded cleaned rows to avoid re-running regex extraction.
"""
function _build_invalid_merged_pdb_references(
        merged_pdb_references::DataFrame,
        cleaned_merged_pdb_rows::DataFrame
)::DataFrame
    merged_reference_keys = select(
        merged_pdb_references,
        :acc,
        :pdb => ByRow(String) => :pdb_reference
    )

    cleaned_reference_keys = unique(
        select(cleaned_merged_pdb_rows, :acc, :pdb_reference),
        [:acc, :pdb_reference]
    )

    # Any raw reference without at least one cleaned PDB code is invalid.
    invalid = antijoin(
        merged_reference_keys,
        cleaned_reference_keys,
        on = [:acc, :pdb_reference]
    )

    transform!(
        invalid,
        :pdb_reference => ByRow(_ -> PDB_LOOKUP_FAILURE_REASON_INVALID_FORMAT) => :failure_reason
    )
    sort!(invalid, [:pdb_reference, :acc])
    return invalid
end

# ---------------------------------------------------------------------------
# Release-date/status annotation helpers.
# ---------------------------------------------------------------------------
"""
    _build_pdb_release_map_row(
        acc,
        pdb,
        pdb_release_dates,
        missing_release_pdb_status,
        default_pdb_status
    )

Build one annotated `(acc, pdb)` output row with release-date and holdings-status
metadata.
"""
function _build_pdb_release_map_row(
        acc::AbstractString,
        pdb::AbstractString,
        pdb_release_dates::AbstractDict,
        missing_release_pdb_status::AbstractDict,
        default_pdb_status::NamedTuple
)
    release_date = get(pdb_release_dates, String(pdb), missing)
    status_record = get(missing_release_pdb_status, String(pdb), default_pdb_status)

    return (
        acc = String(acc),
        pdb = String(pdb),
        pdb_release_date = release_date,
        pdb_lookup_classification = classify_pdb_lookup_status(
            release_date,
            status_record
        ),
        rcsb_http_status = status_record.http_status,
        rcsb_status = status_record.status,
        rcsb_status_code = status_record.status_code,
        replacement_pdb = status_record.replacement_pdb
    )
end

"""
    _replacement_pdb_codes(missing_release_pdb_status)

Return unique replacement PDB IDs extracted from holdings-status records.
"""
function _replacement_pdb_codes(missing_release_pdb_status::AbstractDict)::Vector{String}
    replacement_codes = String[]

    for status_record in values(missing_release_pdb_status)
        ismissing(status_record.replacement_pdb) && continue
        push!(replacement_codes, String(status_record.replacement_pdb))
    end

    return sort(unique(replacement_codes))
end

"""
    _missing_release_pdb_codes(pdb_codes, pdb_release_dates)

Return canonical PDB IDs with missing release-date records.
"""
function _missing_release_pdb_codes(
        pdb_codes::AbstractVector{<:AbstractString},
        pdb_release_dates::AbstractDict
)::Vector{String}
    return sort(unique([String(pdb_code)
                        for pdb_code in pdb_codes
                        if
                        ismissing(get(pdb_release_dates, String(pdb_code), missing))]))
end

"""
    _append_missing_release_statuses_for_codes!(
        missing_release_pdb_status,
        pdb_release_dates,
        pdb_codes;
        holding_status_fetcher
    )

For `pdb_codes` with missing release date, fetch holdings statuses (if not already
present) and merge them into `missing_release_pdb_status`.
"""
function _append_missing_release_statuses_for_codes!(
        missing_release_pdb_status::AbstractDict,
        pdb_release_dates::AbstractDict,
        pdb_codes::AbstractVector{<:AbstractString};
        holding_status_fetcher::Function
)::Nothing
    status_missing_codes = _missing_release_pdb_codes(pdb_codes, pdb_release_dates)
    already_present_codes = Set(String.(collect(keys(missing_release_pdb_status))))
    additional_status_codes = sort(setdiff(status_missing_codes, collect(already_present_codes)))

    # Query holdings only for unresolved IDs that were not fetched before.
    isempty(additional_status_codes) && return nothing
    merge!(
        missing_release_pdb_status,
        holding_status_fetcher(additional_status_codes)
    )
    return nothing
end

"""
    _append_replacement_rows!(
        rows,
        disprot_pdb_map,
        missing_release_pdb_status,
        pdb_release_dates,
        default_pdb_status
    )

Append one row per replacement PDB for obsolete entries, preserving the source
UniProt accession. A replacement row is only added when `(acc, replacement_pdb)`
is not already present.
"""
function _append_replacement_rows!(
        rows::Vector{NamedTuple},
        disprot_pdb_map::DataFrame,
        missing_release_pdb_status::AbstractDict,
        pdb_release_dates::AbstractDict,
        default_pdb_status::NamedTuple
)::Nothing
    # Tracks existing `(acc, pdb)` so replacement rows are additive but never duplicated.
    known_pairs = Set(zip(disprot_pdb_map.acc, disprot_pdb_map.pdb))

    for row in eachrow(disprot_pdb_map)
        status_record = get(missing_release_pdb_status, row.pdb, default_pdb_status)
        replacement_pdb = status_record.replacement_pdb
        ismissing(replacement_pdb) && continue

        replacement_pair = (String(row.acc), String(replacement_pdb))
        replacement_pair in known_pairs && continue
        push!(known_pairs, replacement_pair)

        push!(
            rows,
            _build_pdb_release_map_row(
                row.acc,
                String(replacement_pdb),
                pdb_release_dates,
                missing_release_pdb_status,
                default_pdb_status
            )
        )
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Main builders: comparison-level and raw-input-level dataset generation.
# ---------------------------------------------------------------------------
"""
    build_pdb_release_map_from_comparison(
        disprot_sifts_pdb_comparison;
        batch_size = 300,
        release_date_fetcher = fetch_pdb_initial_release_dates,
        holding_status_fetcher = fetch_pdb_holding_statuses
    )

Build cleaned UniProt-PDB tables and release-date/status annotations starting
from a merged DisProt/SIFTS comparison table.

Returns a named tuple containing intermediate and final tables, including:
- `disprot_pdb_map_with_release`
- `failed_pdb_lookups`
- `invalid_merged_pdb_references`
"""
function build_pdb_release_map_from_comparison(
        disprot_sifts_pdb_comparison::DataFrame;
        batch_size::Integer = 300,
        release_date_fetcher::Function = fetch_pdb_initial_release_dates,
        holding_status_fetcher::Function = fetch_pdb_holding_statuses
)
    # 1) Preserve merged raw references for provenance and diagnostics.
    merged_pdb_references = unique(
        select(disprot_sifts_pdb_comparison, :acc, :pdb),
        [:acc, :pdb]
    )

    # 2) Expand/clean raw reference strings into canonical `(acc, pdb)` rows.
    cleaned_merged_pdb_rows = _expand_merged_pdb_references(merged_pdb_references)

    disprot_pdb_map = unique(
        select(cleaned_merged_pdb_rows, :acc, :pdb),
        [:acc, :pdb]
    )
    sort!(disprot_pdb_map, [:acc, :pdb])

    # 3) Collect raw merged references that could not be parsed into PDB IDs.
    invalid_merged_pdb_references = _build_invalid_merged_pdb_references(
        merged_pdb_references,
        cleaned_merged_pdb_rows
    )

    # 4) Fetch release dates (primary source) for cleaned canonical IDs.
    pdb_release_dates = release_date_fetcher(
        disprot_pdb_map.pdb;
        batch_size = batch_size
    )

    missing_release_pdb_codes = _missing_release_pdb_codes(
        disprot_pdb_map.pdb,
        pdb_release_dates
    )

    # 5) Fetch holdings status (fallback source) for unresolved release dates.
    missing_release_pdb_status = holding_status_fetcher(missing_release_pdb_codes)

    # 6) Add replacement IDs for obsolete PDBs, then annotate replacements too.
    replacement_pdb_codes = _replacement_pdb_codes(missing_release_pdb_status)
    if !isempty(replacement_pdb_codes)
        merge!(
            pdb_release_dates,
            release_date_fetcher(replacement_pdb_codes; batch_size = batch_size)
        )
    end
    _append_missing_release_statuses_for_codes!(
        missing_release_pdb_status,
        pdb_release_dates,
        replacement_pdb_codes;
        holding_status_fetcher = holding_status_fetcher
    )

    default_pdb_status = (
        http_status = 0,
        status = missing,
        status_code = missing,
        replacement_pdb = missing
    )

    # 7) Build the final annotated table and append non-duplicate replacement rows.
    disprot_pdb_map_with_release = let rows = NamedTuple[]
        for row in eachrow(disprot_pdb_map)
            push!(
                rows,
                _build_pdb_release_map_row(
                    row.acc,
                    row.pdb,
                    pdb_release_dates,
                    missing_release_pdb_status,
                    default_pdb_status
                )
            )
        end
        _append_replacement_rows!(
            rows,
            disprot_pdb_map,
            missing_release_pdb_status,
            pdb_release_dates,
            default_pdb_status
        )
        DataFrame(rows)
    end
    sort!(disprot_pdb_map_with_release, [:acc, :pdb])

    missing_release_pdb_codes = _missing_release_pdb_codes(
        disprot_pdb_map_with_release.pdb,
        pdb_release_dates
    )

    failed_pdb_lookups = subset(
        disprot_pdb_map_with_release,
        :pdb_release_date => ByRow(ismissing)
    )
    sort!(failed_pdb_lookups, [:pdb, :acc])

    # 8) Return both final outputs and intermediate tables for notebook analysis.
    return (
        merged_pdb_references = merged_pdb_references,
        cleaned_merged_pdb_rows = cleaned_merged_pdb_rows,
        disprot_pdb_map = disprot_pdb_map,
        invalid_merged_pdb_references = invalid_merged_pdb_references,
        pdb_release_dates = pdb_release_dates,
        missing_release_pdb_codes = missing_release_pdb_codes,
        missing_release_pdb_status = missing_release_pdb_status,
        disprot_pdb_map_with_release = disprot_pdb_map_with_release,
        failed_pdb_lookups = failed_pdb_lookups
    )
end

"""
    build_pdb_release_map_from_raw(
        raw_data_dir;
        batch_size = 300,
        release_date_fetcher = fetch_pdb_initial_release_dates,
        holding_status_fetcher = fetch_pdb_holding_statuses
    )

Build the full UniProt-PDB release-date dataset from repository raw inputs.
"""
function build_pdb_release_map_from_raw(
        raw_data_dir::AbstractString;
        batch_size::Integer = 300,
        release_date_fetcher::Function = fetch_pdb_initial_release_dates,
        holding_status_fetcher::Function = fetch_pdb_holding_statuses
)
    # Resolve expected repository-relative raw input paths.
    disprot_regions_path = joinpath(raw_data_dir, DISPROT_REGIONS_RELATIVE_PATH)
    sifts_uniprot_path = joinpath(raw_data_dir, SIFTS_UNIPROT_RELATIVE_PATH)

    # Build raw DisProt and SIFTS inputs for the comparison table.
    disprot_entries = _read_disprot_entries(disprot_regions_path)
    disprot_uniprot_ids = sort(unique([String(get(protein, "acc", ""))
                                       for protein in disprot_entries
                                       if !isempty(String(get(protein, "acc", "")))]))

    region_level_pdb_refs = collect_disprot_region_level_pdb_refs(disprot_entries)
    disprot_uniprot_pdb_map = build_unique_uniprot_pdb_map(region_level_pdb_refs)

    sifts_uniprot_df = _read_sifts_uniprot_df(sifts_uniprot_path)
    sifts_uniprot_pdb_map = build_sifts_uniprot_pdb_map_for_disprot(
        sifts_uniprot_df,
        disprot_uniprot_ids
    )

    disprot_sifts_pdb_comparison = build_disprot_sifts_pdb_comparison(
        disprot_uniprot_pdb_map,
        sifts_uniprot_pdb_map
    )

    # Delegate release-date/status annotation to the comparison-level builder.
    release_tables = build_pdb_release_map_from_comparison(
        disprot_sifts_pdb_comparison;
        batch_size = batch_size,
        release_date_fetcher = release_date_fetcher,
        holding_status_fetcher = holding_status_fetcher
    )

    # Return both source tables and final derived release table for analysis.
    return (
        disprot_entries = disprot_entries,
        disprot_uniprot_ids = disprot_uniprot_ids,
        region_level_pdb_refs = region_level_pdb_refs,
        disprot_uniprot_pdb_map = disprot_uniprot_pdb_map,
        sifts_uniprot_df = sifts_uniprot_df,
        sifts_uniprot_pdb_map = sifts_uniprot_pdb_map,
        disprot_sifts_pdb_comparison = disprot_sifts_pdb_comparison,
        merged_pdb_references = release_tables.merged_pdb_references,
        cleaned_merged_pdb_rows = release_tables.cleaned_merged_pdb_rows,
        disprot_pdb_map = release_tables.disprot_pdb_map,
        invalid_merged_pdb_references = release_tables.invalid_merged_pdb_references,
        pdb_release_dates = release_tables.pdb_release_dates,
        missing_release_pdb_codes = release_tables.missing_release_pdb_codes,
        missing_release_pdb_status = release_tables.missing_release_pdb_status,
        disprot_pdb_map_with_release = release_tables.disprot_pdb_map_with_release,
        failed_pdb_lookups = release_tables.failed_pdb_lookups
    )
end

# ---------------------------------------------------------------------------
# Output writers and README rendering helpers.
# ---------------------------------------------------------------------------
"""
    write_pdb_release_map_csv_gz(disprot_pdb_map_with_release, output_path)

Write `disprot_pdb_map_with_release` as compressed CSV to `output_path`.
"""
function write_pdb_release_map_csv_gz(
        disprot_pdb_map_with_release::DataFrame,
        output_path::AbstractString
)::Nothing
    mkpath(dirname(output_path))

    open(output_path, "w") do io
        gzip_io = GzipCompressorStream(io)
        try
            CSV.write(gzip_io, disprot_pdb_map_with_release)
        finally
            close(gzip_io)
        end
    end

    return nothing
end

"""
    _markdown_value(value)

Render one value for README markdown output, using `missing` for absent values.
"""
function _markdown_value(value)::String
    return ismissing(value) ? "missing" : string(value)
end

"""
    write_pdb_release_readme(
        readme_path;
        output_csv_path,
        disprot_pdb_map_with_release,
        failed_pdb_lookups,
        invalid_merged_pdb_references,
        batch_size,
        generated_at_utc = utc_now_iso8601()
    )

Write a README report for the derived PDB release-date dataset.
"""
function write_pdb_release_readme(
        readme_path::AbstractString;
        output_csv_path::AbstractString,
        disprot_pdb_map_with_release::DataFrame,
        failed_pdb_lookups::DataFrame,
        invalid_merged_pdb_references::DataFrame,
        batch_size::Integer,
        generated_at_utc::AbstractString = utc_now_iso8601()
)::Nothing
    batch_size > 0 || throw(ArgumentError("`batch_size` must be greater than zero."))
    mkpath(dirname(readme_path))

    output_exists = isfile(output_csv_path)
    output_size_bytes = output_exists ? filesize(output_csv_path) : missing
    missing_unique_pdb = length(unique(failed_pdb_lookups.pdb))

    content = """
# Derived UniProt-PDB Release-Date Map

- Generated at (UTC): `$(generated_at_utc)`
- GraphQL batch size: `$(batch_size)`
- Output file: `$(output_csv_path)`
- Output exists: `$(output_exists)`
- Output size (bytes): `$(_markdown_value(output_size_bytes))`
- Rows in cleaned UniProt-PDB map: `$(size(disprot_pdb_map_with_release, 1))`
- Rows with missing release date: `$(size(failed_pdb_lookups, 1))`
- Unique PDB IDs with missing release date: `$(missing_unique_pdb)`
- Rows with invalid raw PDB references: `$(size(invalid_merged_pdb_references, 1))`

## Failing Cleaned PDB Lookups
"""

    if size(failed_pdb_lookups, 1) == 0
        content *= "\n- None\n"
    else
        for row in eachrow(failed_pdb_lookups)
            content *= """
- `$(row.acc)` | `$(row.pdb)` | classification=`$(_markdown_value(row.pdb_lookup_classification))` | status=`$(_markdown_value(row.rcsb_status))` | status_code=`$(_markdown_value(row.rcsb_status_code))` | replacement=`$(_markdown_value(row.replacement_pdb))`
"""
        end
    end

    content *= "\n## Invalid Raw PDB References\n"
    if size(invalid_merged_pdb_references, 1) == 0
        content *= "\n- None\n"
    else
        for row in eachrow(invalid_merged_pdb_references)
            content *= """
- `$(row.acc)` | raw_reference=`$(row.pdb_reference)` | reason=`$(row.failure_reason)`
"""
        end
    end

    write(readme_path, content)
    return nothing
end
