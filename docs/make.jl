using Documenter, Arpack

makedocs(
    format = Documenter.HTML(
        canonical = "https://julialinearalgebra.github.io/Arpack.jl/stable/",
    ),
    sitename = "Arpack.jl",
    modules = [Arpack],
    pages = [
        "Home" => "index.md",
        "eigs" => "eigs.md",
        "eigs (generalized)" => "eigs_gen.md",
        "svds" => "svds.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaLinearAlgebra/Arpack.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
