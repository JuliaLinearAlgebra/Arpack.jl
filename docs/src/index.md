# [Arpack.jl](@id man-arpack)

```@meta
DocTestSetup = :(using Arpack, LinearAlgebra, SparseArrays)
```

This package provides bindings to
[ARPACK](http://www.caam.rice.edu/software/ARPACK/), which can be used
to perform iterative solutions for eigensystems (using [`eigs`](@ref))
or singular value decompositions (using [`svds`](@ref)).

**Notes**

1. The ARPACK Fortran library is not re-entrant. `Arpack.jl` should only be used from one thread in a Julia program.

2. ARPACK uses a random starting vector by default. This causes the phase of the singular vectors to be random (or just the sign, for real values).

```@meta
DocTestSetup = nothing
```
