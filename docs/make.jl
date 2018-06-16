using Documenter, Arpack

makedocs(
    format = :html,
    sitename = "Arpack.jl",
    modules = [Arpack],
    pages = [
        "index.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaLinearAlgebra/Arpack.jl.git",
    target = "build",
    julia  = "0.7-",
    deps = nothing,
    make = nothing
)
