# [Arpack.jl](@id man-arpack)

```@meta
DocTestSetup = :(using Arpack, LinearAlgebra, SparseArrays)
```

This package provides bindings to [ARPACK](http://www.caam.rice.edu/software/ARPACK/), which
can be used to perform iterative solutions for eigensystems (using [`eigs`](@ref))
or singular value decompositions (using [`svds`](@ref)).

!!! note
    The ARPACK Fortran library is not re-entrant. `Arpack.jl` should only be used from one thread in a Julia program.

```@meta
DocTestSetup = nothing
```
