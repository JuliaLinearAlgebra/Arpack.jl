using Documenter, Arpack

makedocs(
    format = Documenter.HTML(
        canonical = "https://julialinearalgebra.github.io/Arpack.jl/stable/",
    ),
    sitename = "Arpack.jl",
    modules = [Arpack],
    pages = [
        "Home" => "index.md",
        "Standard Eigen Decomposition" => "eigs.md",
        "Generalized Eigen Decomposition" => "eigs_gen.md",
        "Singular Value Decomposition" => "svds.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaLinearAlgebra/Arpack.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
