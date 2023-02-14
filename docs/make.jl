using BoostedDarkMatter
using Documenter

DocMeta.setdocmeta!(BoostedDarkMatter, :DocTestSetup, :(using BoostedDarkMatter); recursive=true)

makedocs(;
    modules=[BoostedDarkMatter],
    authors="Chen Xia <xiachen@sjtu.edu.cn> and contributors",
    repo="https://github.com/physcxia/BoostedDarkMatter.jl/blob/{commit}{path}#{line}",
    sitename="BoostedDarkMatter.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://physcxia.github.io/BoostedDarkMatter.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/physcxia/BoostedDarkMatter.jl",
    devbranch="main",
)
