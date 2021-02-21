# This file is a part of Julia. License is MIT: https://julialang.org/license
"""
Arnoldi and Lanczos iteration for computing eigenvalues
"""
module Arpack

# Load in our binary dependencies
using Arpack_jll

using LinearAlgebra: BlasFloat, BlasInt, Diagonal, I, SVD, UniformScaling,
                     checksquare, factorize, ishermitian, issymmetric, mul!,
                     rmul!, qr!
import LinearAlgebra

export eigs, svds

include("libarpack.jl")

## eigs
"""
    eigs(A; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, explicittransform=:auto, v0=zeros((0,)), check=true) -> (d,[v,],nconv,niter,nmult,resid)

Computes eigenvalues `d` of `A` using implicitly restarted Lanczos or Arnoldi iterations for real symmetric or
general nonsymmetric matrices respectively. See [the manual](@ref man-eigs) for more information.

`eigs` returns the `nev` requested eigenvalues in `d`, the corresponding Ritz vectors `v`
(only if `ritzvec=true`), the number of converged eigenvalues `nconv`, the number of
iterations `niter` and the number of matrix vector multiplications `nmult`, as well as the
final residual vector `resid`. The parameter `explicittransform` takes the values `:auto`, `:none`
or `:shiftinvert`, specifying if shift and invert should be explicitly invoked in julia code.

When `check = true`, an error is thrown if maximum number of iterations taken (`info = 1`). This usually means all possible eigenvalues has been found according to ARPACK manual.
When `check = false`, return currently converged eigenvalues when `info = 1`. Only a `@warn` will given.

# Examples
```jldoctest
julia> using LinearAlgebra, Arpack

julia> A = Diagonal(1:4);

julia> λ, ϕ = eigs(A, nev = 2);

julia> λ
2-element Array{Float64,1}:
 3.9999999999999996
 3.000000000000001
```
"""
eigs(A; kwargs...) = eigs(A, I; kwargs...)
eigs(A::AbstractMatrix{<:BlasFloat}, ::UniformScaling; kwargs...) = _eigs(A, I; kwargs...)

eigs(A::AbstractMatrix{T}, B::AbstractMatrix{T}; kwargs...) where {T<:BlasFloat} = _eigs(A, B; kwargs...)
eigs(A::AbstractMatrix{BigFloat}, B::AbstractMatrix...; kwargs...) = throw(MethodError(eigs, Any[A,B,kwargs...]))
eigs(A::AbstractMatrix{BigFloat}, B::UniformScaling; kwargs...) = throw(MethodError(eigs, Any[A,B,kwargs...]))
function eigs(A::AbstractMatrix{T}, ::UniformScaling; kwargs...) where T
    Tnew = typeof(zero(T)/sqrt(one(T)))
    eigs(convert(AbstractMatrix{Tnew}, A), I; kwargs...)
end
function eigs(A::AbstractMatrix, B::AbstractMatrix; kwargs...)
    T = promote_type(eltype(A), eltype(B))
    Tnew = typeof(zero(T)/sqrt(one(T)))
    eigs(convert(AbstractMatrix{Tnew}, A), convert(AbstractMatrix{Tnew}, B); kwargs...)
end
"""
    eigs(A, B; nev=6, ncv=max(20,2*nev+1), which=:LM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,)), check=true) -> (d,[v,],nconv,niter,nmult,resid)

Computes generalized eigenvalues `d` of `A` and `B` using implicitly restarted Lanczos or Arnoldi iterations for real symmetric or general nonsymmetric matrices respectively. See [the manual](@ref man-eigsgen) for more information.

When `check = true`, an error is thrown if maximum number of iterations taken (`info = 1`). This usually means all possible eigenvalues has been found according to ARPACK manual.
When `check = false`, return currently converged eigenvalues when `info = 1`. Only a `@warn` will given.
"""
eigs(A, B; kwargs...) = _eigs(A, B; kwargs...)
function _eigs(A, B;
               nev::Integer=6, ncv::Integer=max(20,2*nev+1), which=:LM,
               tol=0.0, maxiter::Integer=300, sigma=nothing, v0::Vector=zeros(eltype(A),(0,)),
               ritzvec::Bool=true, explicittransform::Symbol=:auto, check::Bool=true)
    n = checksquare(A)

    eigval_postprocess = false; # If we need to shift-and-invert eigvals as postprocessing

    T = eltype(A)
    iscmplx = T <: Complex
    isgeneral = B !== I
    sym = !iscmplx && issymmetric(A) && issymmetric(B)
    nevmax = sym ? n-1 : n-2
    if nevmax <= 0
        throw(ArgumentError("input matrix A is too small. Use eigen instead."))
    end
    if nev > nevmax
        @warn "Adjusting nev from $nev to $nevmax"
        nev = nevmax
    end
    if nev <= 0
        throw(ArgumentError("requested number of eigenvalues (nev) must be ≥ 1, got $nev"))
    end
    ncvmin = nev + (sym ? 1 : 2)
    if ncv < ncvmin
        @warn "Adjusting ncv from $ncv to $ncvmin"
        ncv = ncvmin
    end
    ncv = BlasInt(min(ncv, n))
    bmat = isgeneral ? "G" : "I"
    isshift = sigma !== nothing

    if isa(which,AbstractString)
        @warn "Use symbols instead of strings for specifying which eigenvalues to compute"
        which=Symbol(which)
    end
    if (which != :LM && which != :SM && which != :LR && which != :SR &&
        which != :LI && which != :SI && which != :BE)
        throw(ArgumentError("which must be :LM, :SM, :LR, :SR, :LI, :SI, or :BE, got $(repr(which))"))
    end
    if which == :BE && !sym
        throw(ArgumentError("which=:BE only possible for real symmetric problem"))
    end
    isshift && which == :SM && @warn "Use of :SM in shift-and-invert mode is not recommended, use :LM to find eigenvalues closest to sigma"

    if which==:SM && !isshift # transform into shift-and-invert method with sigma = 0
        isshift=true
        sigma=zero(T)
        which=:LM
    end

    if (explicittransform==:auto)
        # Try to automatically detect if it is good to carry out an explicittransform
        if (isgeneral && (isshift  || which==:LM))
            explicittransform = :shiftinvert
        else
            explicittransform = :none
        end
    end

    if sigma !== nothing && !iscmplx && isa(sigma,Complex)
        throw(ArgumentError("complex shifts for real problems are not yet supported"))
    end
    sigma = isshift ? convert(T,sigma) : zero(T)

    if (explicittransform==:shiftinvert && (which==:LM || which==:LR || which == :LI) && !isgeneral)
        @warn "Explicit transformation with :L* for standard eigenvalue problems has no meaning. Changing to explicittransform=false."
       explicittransform=:none
    end

    sigma0=sigma; # Store for inverted shift-and-invert
    if explicittransform==:shiftinvert
        isgeneral=false
        bmat="I"
        sym=false # Explicit transform destroys symmetry in general
        sigma=zero(T);
        if (isshift) # Try to keep the original meaning of which & sigma
            if (which == :LM)
                which = :SM
            elseif (which == :SM)
                which = :LM
            end
            if (which == :LR)
                which = :SR
            elseif (which == :SR)
                which = :LR
            end
        end

    end

    if !isempty(v0)
        if length(v0) != n
            throw(DimensionMismatch())
        end
        if eltype(v0) != T
            throw(ArgumentError("starting vector must have element type $T, got $(eltype(v0))"))
        end
    end

    whichstr = "LM"
    if which == :BE
        whichstr = "BE"
    end
    if which == :LR
        whichstr = (!sym ? "LR" : "LA")
    end
    if which == :SR
        whichstr = (!sym ? "SR" : "SA")
    end
    if which == :LI
        if !sym
            whichstr = "LI"
        else
            throw(ArgumentError("largest imaginary is meaningless for symmetric eigenvalue problems"))
        end
    end
    if which == :SI
        if !sym
            whichstr = "SI"
        else
            throw(ArgumentError("smallest imaginary is meaningless for symmetric eigenvalue problems"))
        end
    end

    # Refer to ex-*.doc files in ARPACK/DOCUMENTS for calling sequence
    matvecA! = (y, x) -> mul!(y, A, x)
    if !isgeneral || (explicittransform==:shiftinvert)      # Standard problem
        matvecB = x -> x

        if (explicittransform == :none)
            if !isshift         #    Regular mode
                mode       = 1
                solveSI = x->x
            else                #    Shift-invert mode
                mode       = 3
                F = factorize(A - UniformScaling(sigma))
                solveSI = x -> F \ x
            end
        else
            # doing explicit transformation to standard eigprob
            if (which == :LM || which == :LI || which == :LR)
                eigval_postprocess = false # No eigval postprocess necessary the operator is B^{-1}A
                F = factorize(B);
                matvecA! = (y,x) -> (y[:]= F \ (A*x))
            else
                eigval_postprocess = true
                sigma = zero(T);
                F = factorize(sigma0*B - A);
                matvecA! = (y,x) -> (y[:]= F \ (B*x))
            end
            mode = 1;
            solveSI = x -> x;
        end
    else                    # Generalized eigenproblem
        matvecB = x -> B * x
        if !isshift         #    Regular inverse mode
            mode       = 2
            F = factorize(B)
            solveSI = x -> F \ x
        else                #    Shift-invert mode
            mode       = 3
            F = factorize(A - sigma*B)
            solveSI = x -> F \ x
        end
    end

    # Compute the Ritz values and Ritz vectors
    (resid, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, TOL) =
        aupd_wrapper(T, matvecA!, matvecB, solveSI, n, sym, iscmplx, bmat, nev, ncv, whichstr, tol, maxiter, mode, v0, check)
    # Postprocessing to get eigenvalues and eigenvectors
    !check && (iparam[5] < nev) && @warn "nev = $nev, but only $(iparam[5]) found!"
    output = eupd_wrapper(T, n, sym, iscmplx, bmat, check ? nev : iparam[5], whichstr, ritzvec, TOL,
                          resid, ncv, v, ldv, sigma, iparam, ipntr, workd, workl, lworkl, rwork)

    # Issue 10495, 10701: Check that all eigenvalues are converged
    nev = length(output[1])
    nconv = output[ritzvec ? 3 : 2]
    nev ≤ nconv || @warn "Not all wanted Ritz pairs converged. Requested: $nev, converged: $nconv"

    if (eigval_postprocess) # invert the shift-and-inverse
       λ = sigma0 .- 1 ./output[1];
       return (λ, output[2:end]...)
    end

    return output
end


## svds
struct SVDAugmented{T,S} <: AbstractArray{T, 2}
    X::S
    SVDAugmented{T,S}(X::AbstractMatrix) where {T,S} = new(X)
end

function SVDAugmented(A::AbstractMatrix{T}) where T
    Tnew = typeof(zero(T)/sqrt(one(T)))
    Anew = convert(AbstractMatrix{Tnew}, A)
    SVDAugmented{Tnew,typeof(Anew)}(Anew)
end

function LinearAlgebra.mul!(y::StridedVector{T}, A::SVDAugmented{T}, x::StridedVector{T}) where T
    m, mn = size(A.X, 1), length(x)
    mul!( view(y, 1:m), A.X, view(x, m + 1:mn)) # left singular vector
    mul!(view(y, m + 1:mn), adjoint(A.X), view(x, 1:m)) # right singular vector
    return y
end
Base.size(A::SVDAugmented)  = ((+)(size(A.X)...), (+)(size(A.X)...))
LinearAlgebra.ishermitian(A::SVDAugmented) = true

struct AtA_or_AAt{T,S} <: AbstractArray{T, 2}
    A::S
    buffer::Vector{T}
end

function AtA_or_AAt(A)
    T    = eltype(A)
    Tnew = typeof(zero(T)/sqrt(one(T)))
    return AtA_or_AAt{Tnew,typeof(A)}(A, Vector{Tnew}(undef, max(size(A)...)))
end

function LinearAlgebra.mul!(y::StridedVector{T}, A::AtA_or_AAt{T}, x::StridedVector{T}) where T
    if size(A.A, 1) >= size(A.A, 2)
        mul!(A.buffer, A.A, x)
        return mul!(y, adjoint(A.A), A.buffer)
    else
        mul!(A.buffer, adjoint(A.A), x)
        return mul!(y, A.A, A.buffer)
    end
end
Base.size(A::AtA_or_AAt) = ntuple(i -> min(size(A.A)...), Val(2))
LinearAlgebra.ishermitian(s::AtA_or_AAt) = true


svds(A::AbstractMatrix{<:BlasFloat}; kwargs...) = _svds(A; kwargs...)
svds(A::AbstractMatrix{BigFloat}; kwargs...) = throw(MethodError(svds, Any[A, kwargs...]))
function svds(A::AbstractMatrix{T}; kwargs...) where T
    Tnew = typeof(zero(T)/sqrt(one(T)))
    svds(convert(AbstractMatrix{Tnew}, A); kwargs...)
end

"""
    svds(A; nsv=6, ritzvec=true, tol=0.0, maxiter=1000, ncv=2*nsv, v0=zeros((0,))) -> (SVD([left_sv,] s, [right_sv,]), nconv, niter, nmult, resid, check=true)

Computes the largest singular values `s` of `A` using implicitly restarted Lanczos
iterations derived from [`eigs`](@ref). See [the manual](@ref man-svds) for more information.

When `check = true`, an error is thrown if maximum number of iterations taken (`info = 1`). This usually means all possible eigenvalues has been found according to ARPACK manual.
When `check = false`, return currently converged eigenvalues when `info = 1`. Only a `@warn` will given.
"""
svds(A; kwargs...) = _svds(A; kwargs...)
function _orth!(P)
    Q,R = qr!(P)
    _sign(x) = iszero(x) ? one(x) : sign(x)
    rsign = [_sign(R[i,i]) for i in 1:size(R,2)]
    return rmul!(Matrix(Q), Diagonal(rsign))
end
function _svds(X; nsv::Int = 6, ritzvec::Bool = true, tol::Float64 = 0.0, maxiter::Int = 1000, ncv::Int = 2*nsv, v0::Vector=zeros(eltype(X),(0,)), check::Bool=true)
    if nsv < 1
        throw(ArgumentError("number of singular values (nsv) must be ≥ 1, got $nsv"))
    end
    if nsv >= minimum(size(X))
        throw(ArgumentError("number of singular values (nsv) must be < $(minimum(size(X))), got $nsv"))
    end
    m, n = size(X)
    otype = eltype(X)
    if length(v0) ∉ [0,n]
        throw(DimensionMismatch("length of v0, the guess for the starting right Krylov vector, must be 0, or $n, got $(length(v0))"))
    end
    ex    = eigs(AtA_or_AAt(X), I; which = :LM, ritzvec = ritzvec, nev = nsv, tol = tol, maxiter = maxiter, v0=v0, check=check)
    # ind   = [1:2:ncv;]
    # sval  = abs.(ex[1][ind])

    realex1 = real.(ex[1])
    threshold = max(eps(real(otype))*realex1[1], eps(real(otype)))
    firstzero = findfirst(v -> v <= threshold, realex1)
    r = firstzero === nothing ? nsv : firstzero-1 # rank of the decomposition
    realex1[r+1:end] .= zero(real(otype))
    svals = sqrt.(realex1)

    if ritzvec
        # calculating singular vectors
        # left_sv  = sqrt(2) * ex[2][ 1:size(X,1),     ind ] .* sign.(ex[1][ind]')
        if size(X, 1) >= size(X, 2)
            V = ex[2]
            # We cannot assume that X*V is a Matrix even though V is. This is not
            # the case for e.g. LinearMaps.jl so we convert to Matrix explicitly
            U = _orth!(rmul!(convert(Matrix, X*V), Diagonal([inv.(svals[1:r]); ones(nsv-r)])))
        else
            U = ex[2]
            # We cannot assume that X'U is a Matrix even though U is. This is not
            # the case for e.g. LinearMaps.jl so we convert to Matrix explicitly
            V = _orth!(rmul!(convert(Matrix, X'U), Diagonal([inv.(svals[1:r]); ones(nsv-r)])))
        end

        # right_sv = sqrt(2) * ex[2][ size(X,1)+1:end, ind ]
        return (SVD(U, svals, copy(V')), ex[3], ex[4], ex[5], ex[6])
    else
        #The sort is necessary to work around #10329
        return (SVD(zeros(eltype(svals), n, 0),
                    svals,
                    zeros(eltype(svals), 0, m)),
                    ex[2], ex[3], ex[4], ex[5])
    end
end

end # module
