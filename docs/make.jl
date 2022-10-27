using TubePressure
using Documenter

DocMeta.setdocmeta!(TubePressure, :DocTestSetup, :(using TubePressure); recursive=true)

makedocs(;
    modules=[TubePressure],
    authors="Paulo José Saiz Jabardo",
    repo="https://github.com/pjsjipt/TubePressure.jl/blob/{commit}{path}#{line}",
    sitename="TubePressure.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
