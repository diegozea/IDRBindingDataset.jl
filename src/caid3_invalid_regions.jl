using CSV

const CAID_3_INVALID_REGIONS_RELATIVE_PATH = joinpath(
    "caid_3",
    CAID_3_INVALID_REGIONS_FILENAME
)

"""
    _region_id(region)

Return one normalized DisProt `region_id` for `region`, or `nothing` when
missing/empty.
"""
function _region_id(region)::Union{Nothing, String}
    region isa AbstractDict || return nothing
    region_id_value = get(region, "region_id", missing)
    region_id_value isa AbstractString || return nothing

    region_id = strip(String(region_id_value))
    isempty(region_id) && return nothing
    region_id
end

"""
    _region_is_marked_invalid(region, invalid_region_ids)

Return `true` when `region` has one `region_id` present in
`invalid_region_ids`.
"""
function _region_is_marked_invalid(
        region,
        invalid_region_ids::Set{String}
)::Bool
    region_id = _region_id(region)
    region_id === nothing && return false
    region_id in invalid_region_ids
end

"""
    _load_caid3_invalid_region_ids(raw_data_dir;
        relative_path = CAID_3_INVALID_REGIONS_RELATIVE_PATH)

Load CAID-3 invalid region IDs from
`joinpath(raw_data_dir, "caid_3", "invalid_regions_fullCAID3.tsv")`.

Rows with `valid == 0` are returned as invalid; rows with `valid == 1` are
treated as valid and excluded from the output set.
"""
function _load_caid3_invalid_region_ids(
        raw_data_dir::AbstractString;
        relative_path::AbstractString = CAID_3_INVALID_REGIONS_RELATIVE_PATH
)::Set{String}
    tsv_path = joinpath(raw_data_dir, relative_path)
    isfile(tsv_path) || throw(
        ArgumentError("CAID-3 invalid-regions TSV does not exist: $(tsv_path)"),
    )
    filesize(tsv_path) > 0 || throw(
        ArgumentError("CAID-3 invalid-regions TSV is empty: $(tsv_path)"),
    )

    rows = try
        CSV.File(
            tsv_path;
            delim = '\t',
            types = String,
            ignoreemptyrows = true
        )
    catch err
        throw(
            ArgumentError(
            "Failed to parse CAID-3 invalid-regions TSV with CSV.jl: $(tsv_path). Underlying error: $(sprint(showerror, err))."
        ),
        )
    end

    column_names = Set(String.(propertynames(rows)))
    ("region_id" in column_names && "valid" in column_names) || throw(
        ArgumentError(
        "CAID-3 invalid-regions TSV must include `region_id` and `valid` columns: $(tsv_path)",
    ),
    )

    invalid_region_ids = Set{String}()
    for (row_index, row) in enumerate(rows)
        line_number = row_index + 1
        region_id_value = row.region_id
        valid_cell = row.valid

        region_id = region_id_value === missing ? "" : strip(String(region_id_value))
        valid_value = valid_cell === missing ? "" : strip(String(valid_cell))
        isempty(region_id) && continue

        if valid_value == "0"
            push!(invalid_region_ids, region_id)
        elseif valid_value == "1"
            continue
        else
            throw(
                ArgumentError(
                "Invalid `valid` value $(repr(valid_value)) at line $(line_number) in $(tsv_path). Expected `0` or `1`.",
            ),
            )
        end
    end

    invalid_region_ids
end
