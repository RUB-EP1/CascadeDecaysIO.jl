using CascadeDecaysIO
using Documenter

DocMeta.setdocmeta!(
    CascadeDecaysIO,
    :DocTestSetup,
    :(using CascadeDecaysIO);
    recursive = true,
)


makedocs(;
    modules = [CascadeDecaysIO],
    authors = "RUB-EP1 and contributors",
    repo = "https://github.com/RUB-EP1/CascadeDecaysIO.jl/blob/{commit}{path}#{line}",
    sitename = "CascadeDecaysIO.jl",
    doctest = false,
    checkdocs = :none,
    format = Documenter.HTML(;
        canonical = "https://rub-ep1.github.io/CascadeDecaysIO.jl",
        repolink = "https://github.com/RUB-EP1/CascadeDecaysIO.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Writer workflow" => "writer-workflow.md",
        "Schema notes" => "schema.md",
        "API reference" => "api-reference.md",
    ],
)

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/RUB-EP1/CascadeDecaysIO.jl")
end
