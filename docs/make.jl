using IDRBindingDataset
using Documenter

DocMeta.setdocmeta!(IDRBindingDataset, :DocTestSetup, :(using IDRBindingDataset); recursive = true)

makedocs(;
    modules = [IDRBindingDataset],
    authors = "Diego Javier Zea <diegozea@gmail.com> and contributors",
    sitename = "IDRBindingDataset.jl",
    format = Documenter.HTML(;
        canonical = "https://diegozea.github.io/IDRBindingDataset.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Raw Data Workflow" => "raw_data.md",
        "ELM Motif Workflow" => "elm_motifs.md"
    ]
)

deploydocs(;
    repo = "github.com/diegozea/IDRBindingDataset.jl",
    devbranch = "main"
)
