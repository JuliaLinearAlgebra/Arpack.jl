# Arpack

[![CI](https://github.com/JuliaLinearAlgebra/Arpack.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaLinearAlgebra/Arpack.jl/actions/workflows/ci.yml)
[![][docs-stable-img]][docs-stable-url]

Julia wrapper for the [arpack](https://github.com/opencollab/arpack-ng/) library
designed to solve large-scale eigenvalue problems.

## Installation

Install Arpack.jl through the Julia package manager:
```julia
julia> Pkg.add("Arpack")
```

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://arpack.JuliaLinearAlgebra.org/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://arpack.JuliaLinearAlgebra.org/stable/

# Alternate packages

Users running into issues with this package may want to try [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl) or [ArnoldiMethod.jl](https://github.com/haampie/ArnoldiMethod.jl).
