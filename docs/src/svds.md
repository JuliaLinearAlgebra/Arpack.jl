# [svds](@id man-svds)

```@meta
DocTestSetup = :(using Arpack, LinearAlgebra, SparseArrays, Random)
```
Computes the largest singular values `s` of `A` using implicitly restarted Lanczos
iterations derived from [`eigs`](@ref).

**Inputs**

* `A`: Linear operator whose singular values are desired. `A` may be represented as a
  subtype of `AbstractArray`, e.g., a sparse matrix, or any other type supporting the four
  methods `size(A)`, `eltype(A)`, `A * vector`, and `A' * vector`.
* `nsv`: Number of singular values. Default: 6.
* `ritzvec`: If `true`, return the left and right singular vectors `left_sv` and `right_sv`.
   If `false`, omit the singular vectors. Default: `true`.
* `tol`: tolerance, see [`eigs`](@ref).
* `maxiter`: Maximum number of iterations, see [`eigs`](@ref). Default: 1000.
* `ncv`: Maximum size of the Krylov subspace, see [`eigs`](@ref) (there called `nev`). Default: `2*nsv`.
* `v0`: Initial guess for the first Krylov vector. It may have length `min(size(A)...)`, or 0.

**Outputs**

* `svd`: An `SVD` object containing the left singular vectors, the requested values, and the
  right singular vectors. If `ritzvec = false`, the left and right singular vectors will be
  empty. `U`, `S`, `V` and `Vt` can be obtained from the SVD object with `Z.U`, `Z.S`, `Z.V`
  and `Z.Vt`, where `Z = svds(A)[1]` and `U * Diagonal(S) * Vt` is a low-rank approximation
  of `A` with rank `nsv`. Internally `Vt` is stored and hence `Vt` is more efficient to extract than `V`.
* `nconv`: Number of converged singular values.
* `niter`: Number of iterations.
* `nmult`: Number of matrix--vector products used.
* `resid`: Final residual vector.

# Examples

```jldoctest
julia> Random.seed!(123);

julia> A = Diagonal(1:5);

julia> Z = svds(A, nsv = 2)[1];

julia> Z.U
5×2 Array{Float64,2}:
  0.0           7.80626e-18
  0.0          -0.0
 -1.33227e-16   5.35947e-33
 -6.38552e-17   1.0
 -1.0          -6.38552e-17

julia> Z.S
2-element Array{Float64,1}:
 5.0
 3.999999999999999

julia> Z.Vt
2×5 Array{Float64,2}:
 -2.77556e-17  0.0          -2.22045e-16  -7.9819e-17  -1.0
  3.1225e-17   1.89735e-19   0.0           1.0         -8.32667e-17

julia> Z.V
5×2 Adjoint{Float64,Array{Float64,2}}:
 -2.77556e-17   3.1225e-17
  0.0           1.89735e-19
 -2.22045e-16   0.0
 -7.9819e-17    1.0
 -1.0          -8.32667e-17
```


!!! note "Implementation"
    `svds(A)` is formally equivalent to calling [`eigs`](@ref) to perform implicitly restarted
    Lanczos tridiagonalization on the Hermitian matrix ``A^\\prime A`` or ``AA^\\prime`` such
    that the size is smallest.

```@meta
DocTestSetup = nothing
```
