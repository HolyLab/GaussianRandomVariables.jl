using GaussianRandomVariables
using Documenter

DocMeta.setdocmeta!(GaussianRandomVariables, :DocTestSetup, :(using GaussianRandomVariables); recursive=true)

makedocs(;
    modules=[GaussianRandomVariables],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    repo="https://github.com/HolyLab/GaussianRandomVariables.jl/blob/{commit}{path}#{line}",
    sitename="GaussianRandomVariables.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HolyLab.github.io/GaussianRandomVariables.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HolyLab/GaussianRandomVariables.jl",
    devbranch="main",
)
