using KhepriFrame3DD
using Documenter

makedocs(;
    modules=[KhepriFrame3DD],
    authors="António Menezes Leitão <antonio.menezes.leitao@gmail.com>",
    repo="https://github.com/aptmcl/KhepriFrame3DD.jl/blob/{commit}{path}#L{line}",
    sitename="KhepriFrame3DD.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aptmcl.github.io/KhepriFrame3DD.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aptmcl/KhepriFrame3DD.jl",
)
