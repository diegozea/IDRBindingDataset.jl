using CodecZlib
using CSV
using DataFrames

"""
    _write_csv_gz(table, output_path)

Write one table-like object to a gzip-compressed CSV file at `output_path`.
Creates parent directories when needed.
"""
function _write_csv_gz(
        table,
        output_path::AbstractString
)::Nothing
    mkpath(dirname(output_path))
    open(GzipCompressorStream, output_path, "w") do io
        CSV.write(io, table)
    end

    nothing
end

"""
    _write_text_file(output_path, content)

Write plain-text `content` to `output_path`, creating parent directories when
needed.
"""
function _write_text_file(
        output_path::AbstractString,
        content::AbstractString
)::Nothing
    path = String(output_path)
    mkpath(dirname(path))
    write(path, content)

    nothing
end

"""
    _normalize_optional_text(value;
        uppercase_output = false,
        lowercase_output = false,
        strip_wrapping_quotes = false,
        missing_value = missing
    )

Normalize one optional textual value to `String`, returning `missing_value` for
`nothing`, `missing`, and empty strings after trimming.
"""
function _normalize_optional_text(
        value;
        uppercase_output::Bool = false,
        lowercase_output::Bool = false,
        strip_wrapping_quotes::Bool = false,
        missing_value = missing
)::Union{Missing, Nothing, String}
    (uppercase_output && lowercase_output) && throw(
        ArgumentError("Only one of `uppercase_output` or `lowercase_output` may be true."),
    )

    (value === nothing || ismissing(value)) && return missing_value

    text = strip(String(value))
    if strip_wrapping_quotes
        text = strip(strip(text, ['"', '\'']))
    end
    isempty(text) && return missing_value

    uppercase_output && return uppercase(text)
    lowercase_output && return lowercase(text)
    text
end

"""
    _normalize_label_for_report(label; missing_label = "missing")

Normalize one label for reporting, replacing `missing` with `missing_label`.
"""
function _normalize_label_for_report(
        label;
        missing_label::AbstractString = "missing"
)::String
    ismissing(label) && return String(missing_label)
    String(label)
end

"""
    _count_labels_by_value(table, label_column; missing_label = "missing", table_name = "table")

Count rows grouped by normalized values from `label_column`.
"""
function _count_labels_by_value(
        table::AbstractDataFrame,
        label_column::Symbol;
        missing_label::AbstractString = "missing",
        table_name::AbstractString = "table"
)::Dict{String, Int}
    label_column in propertynames(table) || throw(
        ArgumentError(
        "$(table_name) is missing required column: $(label_column)",
    ),
    )

    counts = Dict{String, Int}()
    for value in table[!, label_column]
        label = _normalize_label_for_report(value; missing_label = missing_label)
        counts[label] = get(counts, label, 0) + 1
    end

    counts
end

"""
    _render_label_count_lines(label_counts; ordered_labels, prefix)

Render markdown bullet lines for one label-count dictionary using a stable label
order followed by sorted extra labels.
"""
function _render_label_count_lines(
        label_counts::Dict{String, Int};
        ordered_labels::AbstractVector{<:AbstractString},
        prefix::AbstractString
)::String
    lines = String[]

    for label in ordered_labels
        push!(lines, "- $(prefix) `$(label)`: `$(get(label_counts, String(label), 0))`")
    end

    for label in sort(setdiff(collect(keys(label_counts)), String.(ordered_labels)))
        push!(lines, "- $(prefix) `$(label)`: `$(label_counts[label])`")
    end

    join(lines, "\n")
end
