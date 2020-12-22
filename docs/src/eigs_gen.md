# [Generalized Eigen Decomposition](@id man-eigsgen)

```@meta
DocTestSetup = :(using Arpack, LinearAlgebra, SparseArrays)
```

For the two-input generalized eigensolution version,

`eigs(A, B; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,))) -> (d,[v,],nconv,niter,nmult,resid)`

the following keyword arguments are supported:

* `nev`: Number of eigenvalues
* `ncv`: Number of Krylov vectors used in the computation; should satisfy `nev+1 <= ncv <= n`
  for real symmetric problems and `nev+2 <= ncv <= n` for other problems, where `n` is the
  size of the input matrices `A` and `B`. The default is `ncv = max(20,2*nev+1)`. Note that
  these restrictions limit the input matrix `A` to be of dimension at least 2.
* `which`: type of eigenvalues to compute. See the note below.

| `which` | type of eigenvalues                                                                                                       |
|:--------|:--------------------------------------------------------------------------------------------------------------------------|
| `:LM`   | eigenvalues of largest magnitude (default)                                                                                |
| `:SM`   | eigenvalues of smallest magnitude                                                                                         |
| `:LR`   | eigenvalues of largest real part                                                                                          |
| `:SR`   | eigenvalues of smallest real part                                                                                         |
| `:LI`   | eigenvalues of largest imaginary part (nonsymmetric or complex `A` only)                                                  |
| `:SI`   | eigenvalues of smallest imaginary part (nonsymmetric or complex `A` only)                                                 |
| `:BE`   | compute half of the eigenvalues from each end of the spectrum, biased in favor of the high end. (real symmetric `A` only) |

* `tol`: relative tolerance used in the convergence criterion for eigenvalues, similar to
     `tol` in the [`eigs(A)`](@ref) method for the ordinary eigenvalue
     problem, but effectively for the eigenvalues of ``B^{-1} A`` instead of ``A``.
     See the documentation for the ordinary eigenvalue problem in
     [`eigs(A)`](@ref) and the accompanying note about `tol`.
* `maxiter`: Maximum number of iterations (default = 300)
* `sigma`: Specifies the level shift used in inverse iteration. If `nothing` (default),
  defaults to ordinary (forward) iterations. Otherwise, find eigenvalues close to `sigma`
  using shift and invert iterations.
* `ritzvec`: Returns the Ritz vectors `v` (eigenvectors) if `true`
* `v0`: starting vector from which to start the iterations

`eigs` returns the `nev` requested eigenvalues in `d`, the corresponding Ritz vectors `v`
(only if `ritzvec=true`), the number of converged eigenvalues `nconv`, the number of
iterations `niter` and the number of matrix vector multiplications `nmult`, as well as the
final residual vector `resid`.

We can see the various keywords in action in the following examples:
```jldoctest; filter = r"(1|2)-element Array{(Float64|Complex{Float64}),1}:\n (.|\s)*$"
julia> A = sparse(1.0I, 4, 4); B = Diagonal(1:4);

julia> λ, ϕ = eigs(A, B, nev = 2);

julia> λ
2-element Array{Float64,1}:
 1.0000000000000002
 0.5

julia> A = Diagonal([1, -2im, 3, 4im]); B = sparse(1.0I, 4, 4);

julia> λ, ϕ = eigs(A, B, nev=1, which=:SI);

julia> λ
1-element Array{Complex{Float64},1}:
 -1.5720931501039814e-16 - 1.9999999999999984im

julia> λ, ϕ = eigs(A, B, nev=1, which=:LI);

julia> λ
1-element Array{Complex{Float64},1}:
 0.0 + 4.000000000000002im
```

!!! note
    The `sigma` and `which` keywords interact: the description of eigenvalues searched for by
    `which` do *not* necessarily refer to the eigenvalue problem ``Av = Bv\lambda``, but rather
    the linear operator constructed by the specification of the iteration mode implied by `sigma`.

    | `sigma`         | iteration mode                   | `which` refers to the problem      |
    |:----------------|:---------------------------------|:-----------------------------------|
    | `nothing`       | ordinary (forward)               | ``Av = Bv\lambda``                |
    | real or complex | inverse with level shift `sigma` | ``(A - \sigma B )^{-1}B = v\nu`` |


```@meta
DocTestSetup = nothing
```
