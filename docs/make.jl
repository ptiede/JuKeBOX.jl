using JuKeBOX
using Documenter

DocMeta.setdocmeta!(JuKeBOX, :DocTestSetup, :(using JuKeBOX); recursive=true)

makedocs(;
    modules=[JuKeBOX],
    authors="Paul Tiede <ptiede91@gmail.com> and contributors",
    repo="https://github.com/ptiede/JuKeBOX.jl/blob/{commit}{path}#{line}",
    sitename="JuKeBOX.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ptiede.github.io/JuKeBOX.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ptiede/JuKeBOX.jl",
    devbranch="main",
)
