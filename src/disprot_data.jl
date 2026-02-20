using CodecZlib
using JSON

# ---------------------------------------------------------------------------
# Shared DisProt input locations/readers.
# ---------------------------------------------------------------------------
const DISPROT_REGIONS_RELATIVE_PATH = joinpath("disprot", "disprot_regions.json.gz")

"""
    _read_disprot_entries(disprot_regions_path)

Read a DisProt compressed JSON regions file and return the `data` entries array.
"""
function _read_disprot_entries(disprot_regions_path::AbstractString)
    payload = open(disprot_regions_path) do io
        JSON.parse(read(GzipDecompressorStream(io), String))
    end
    return get(payload, "data", Any[])
end

"""
    read_disprot_entries_from_raw_data(raw_data_dir)

Read DisProt entries from `raw_data_dir/disprot/disprot_regions.json.gz`.
"""
function read_disprot_entries_from_raw_data(raw_data_dir::AbstractString)
    disprot_regions_path = string(joinpath(raw_data_dir, DISPROT_REGIONS_RELATIVE_PATH))
    return _read_disprot_entries(disprot_regions_path)
end
