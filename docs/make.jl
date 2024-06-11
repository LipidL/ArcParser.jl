using ArcParser
using Documenter

DocMeta.setdocmeta!(ArcParser, :DocTestSetup, :(using ArcParser); recursive=true)

makedocs(;
    modules=[ArcParser],
    authors="lipid <lipid.l@qq.com> and contributors",
    sitename="ArcParser.jl",
    format=Documenter.HTML(;
        canonical="https://LipidL.github.io/ArcParser.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LipidL/ArcParser.jl",
    devbranch="main",
)
