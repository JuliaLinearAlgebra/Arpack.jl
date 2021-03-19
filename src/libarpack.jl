# This file is a part of Julia. License is MIT: https://julialang.org/licenOAse

import LinearAlgebra: BlasInt
using Logging

# A convenient shortcut to show unexpected behavior from libarpack
const ERR_UNEXPECTED_BEHAVIOR = -999

struct XYAUPD_Exception <: Exception
    info::BlasInt
end

const AUPD_ERRORS = [
    (3, "No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV. "),
    (2, "No longer an informational error. Deprecated starting with release 2 of ARPACK."),
    (1, """Maximum number of iterations taken. All possible eigenvalues of OP has been found.
          IPARAM(5) returns the number of wanted converged Ritz values."""),
    (0, "Normal exit."),
    (-1, "N must be positive."),
    (-2, "NEV must be positive."),
    (-3, "NCV-NEV >= 2 and less than or equal to N."),
    (-4, "The maximum number of Arnoldi update iterations allowed must be greater than zero."),
    (-5, " WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'"),
    (-6, "BMAT must be one of 'I' or 'G'."),
    (-7, "Length of private work array WORKL is not sufficient."),
    (-8, "Error return from LAPACK eigenvalue calculation."),
    (-9, "Starting vector is zero."),
    (-10, "IPARAM(7) must be 1,2,3,4."),
    (-11, "IPARAM(7) = 1 and BMAT = 'G' are incompatible."),
    (-12, "IPARAM(1) must be equal to 0 or 1."),
    (-13, "NEV and WHICH = 'BE' are incompatible."),
    (-9999, """Could not build an Arnoldi factorization.
           IPARAM(5) returns the size of the current Arnoldi factorization.
           The user is advised to check that enough workspace and array storage has been allocated.""")
]

function Base.showerror(io::IO, ex::XYAUPD_Exception)
    info = ex.info

    if info == ERR_UNEXPECTED_BEHAVIOR
        @error "XYAUPD_Exception: Undefined error"
    else
      idx = searchsorted(AUPD_ERRORS, info, by=first, rev=true)
      if isempty(idx)
          @error "XYAUPD_Exception: Please check XYAUPD error codes in the ARPACK manual." info
      else
          @error "XYAUPD_Exception: $(last(AUPD_ERRORS[first(idx)]))" info
      end
    end
end

struct XYEUPD_Exception <: Exception
    info::BlasInt
end

const EUPD_ERRORS = [
    (1, """The Schur form computed by LAPACK routine lahqr could not be reordered by LAPACK routine trsen.
       Re-enter subroutine neupd with IPARAM(5)NCV and increase the size of the arrays DR and DI to have dimension at least dimension NCV and allocate at least NCV columns for Z.
       NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs."""),
    (0, "Normal exit."),
    (-1, "N must be positive."),
    (-2, "NEV must be positive."),
    (-3, "NCV-NEV >= 2 and less than or equal to N."),
    (-5, "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'"),
    (-6, "BMAT must be one of 'I' or 'G'."),
    (-7, "Length of private work WORKL array is not sufficient."),
    (-8, """Error return from calculation of a real Schur form.
        "Informational error from LAPACK routine lahqr ."""),
    (-9, """Error return from calculation of eigenvectors.
        "Informational error from LAPACK routine dtrevc."""),
    (-10, "IPARAM(7) must be 1,2,3,4."),
    (-11, "IPARAM(7) = 1 and BMAT = 'G' are incompatible."),
    (-12, "HOWMNY = 'S' not yet implemented"),
    (-13, "HOWMNY must be one of 'A' or 'P' if RVEC = .true."),
    (-14, "DNAUPD  did not find any eigenvalues to sufficient accuracy."),
    (-15, """DNEUPD got a different count of the number of converged Ritz values than NAUPD got.
          This indicates the user probably made an error in passing data from NAUPD to NEUPD or that the data was modified before entering NEUPD""")
]

function Base.showerror(io::IO, ex::XYEUPD_Exception)
    info = ex.info
    if info == ERR_UNEXPECTED_BEHAVIOR
        @error "XYEUPD_Exception: Undefined error"
    else
      idx = searchsorted(EUPD_ERRORS, info, by=first, rev=true)
      if isempty(idx)
          @error "XYEUPD_Exception: Please check XYEUPD error codes in the ARPACK manual." info
      else
          @error "XYEUPD_Exception: $(last(EUPD_ERRORS[first(idx)]))" info
      end
    end
end

## aupd and eupd wrappers

function aupd_wrapper(T, matvecA!::Function, matvecB::Function, solveSI::Function, n::Integer,
                      sym::Bool, cmplx::Bool, bmat,
                      nev::Integer, ncv::Integer, which,
                      tol::Real, maxiter::Integer, mode::Integer, v0::Vector, check::Integer)
    lworkl = cmplx ? ncv * (3*ncv + 5) : (sym ? ncv * (ncv + 8) :  ncv * (3*ncv + 6) )
    TR = cmplx ? T.types[1] : T
    TOL = Ref{TR}(tol)

    v     = Matrix{T}(undef, n, ncv)
    workd = Vector{T}(undef, 3*n)
    workl = Vector{T}(undef, lworkl)
    rwork = cmplx ? Vector{TR}(undef, ncv) : Vector{TR}()

    if isempty(v0)
        resid = Vector{T}(undef, n)
        info  = Ref{BlasInt}(0)
    else
        resid = deepcopy(v0)
        info  = Ref{BlasInt}(1)
    end
    iparam = zeros(BlasInt, 11)
    ipntr  = zeros(BlasInt, (sym && !cmplx) ? 11 : 14)
    ido    = Ref{BlasInt}(0)

    iparam[1] = BlasInt(1)       # ishifts
    iparam[3] = BlasInt(maxiter) # maxiter
    iparam[7] = BlasInt(mode)    # mode

    zernm1 = 0:(n-1)

    while true
        if cmplx
            naupd(ido, bmat, n, which, nev, TOL, resid, ncv, v, n,
                  iparam, ipntr, workd, workl, lworkl, rwork, info)
        elseif sym
            saupd(ido, bmat, n, which, nev, TOL, resid, ncv, v, n,
                  iparam, ipntr, workd, workl, lworkl, info)
        else
            naupd(ido, bmat, n, which, nev, TOL, resid, ncv, v, n,
                  iparam, ipntr, workd, workl, lworkl, info)
        end
        if info[] != 0
            if info[] == 1 && check != 0
                return (resid, v, n, iparam, ipntr, workd, workl, lworkl, rwork, TOL)
            end
            throw(XYAUPD_Exception(info[]))
        end

        x = view(workd, ipntr[1] .+ zernm1)
        y = view(workd, ipntr[2] .+ zernm1)
        if mode == 1  # corresponds to dsdrv1, dndrv1 or zndrv1
            if ido[] == -1 || ido[] == 1
                matvecA!(y, x)
            elseif ido[] == 99
                break
            else
                throw(XYAUPD_Exception(ERR_UNEXPECTED_BEHAVIOR))
            end
        elseif mode == 3 && bmat == "I" # corresponds to dsdrv2, dndrv2 or zndrv2
            if ido[] == -1 || ido[] == 1
                y[:] = solveSI(x)
            elseif ido[] == 99
                break
            else
                throw(XYAUPD_Exception(ERR_UNEXPECTED_BEHAVIOR))
            end
        elseif mode == 2 # corresponds to dsdrv3, dndrv3 or zndrv3
            if ido[] == -1 || ido[] == 1
                matvecA!(y, x)
                if sym
                    x[:] = y    # overwrite as per Remark 5 in dsaupd.f
                end
                y[:] = solveSI(y)
            elseif ido[] == 2
                y[:] = matvecB(x)
            elseif ido[] == 99
                break
            else
                throw(XYAUPD_Exception(ERR_UNEXPECTED_BEHAVIOR))
            end
        elseif mode == 3 && bmat == "G" # corresponds to dsdrv4, dndrv4 or zndrv4
            if ido[] == -1
                y[:] = solveSI(matvecB(x))
            elseif  ido[] == 1
                y[:] = solveSI(view(workd,ipntr[3] .+ zernm1))
            elseif ido[] == 2
                y[:] = matvecB(x)
            elseif ido[] == 99
                break
            else
                throw(XYAUPD_Exception(ERR_UNEXPECTED_BEHAVIOR))
            end
        else
            throw(ArgumentError("ARPACK mode ($mode) not yet supported"))
        end
    end

    return (resid, v, n, iparam, ipntr, workd, workl, lworkl, rwork, TOL)
end

function eupd_wrapper(T, n::Integer, sym::Bool, cmplx::Bool, bmat,
                      nev::Integer, which, ritzvec::Bool,
                      TOL::Ref, resid, ncv::Integer, v, ldv, sigma, iparam, ipntr,
                      workd, workl, lworkl, rwork)
    howmny = "A"
    select = Vector{BlasInt}(undef, ncv)
    info   = Ref{BlasInt}(0)

    dmap = if which == "LM" || which == "SM"
      abs
    elseif which == "LR" || which == "LA" || which == "BE" || which == "SR" || which == "SA"
      real
    elseif which == "LI" || which == "SI"
      abs ∘ imag # ARPACK returns largest,smallest abs(imaginary) (complex pairs come together)
    else
      error("unknown which string $which")
    end

    rev = which[1] == 'L'

    if iparam[7] == 3 # shift-and-invert
      dmap = dmap ∘ (x -> 1 / (x - sigma))
    end

    if cmplx
        d = Vector{T}(undef, nev+1)
        sigmar = Ref{T}(sigma)
        workev = Vector{T}(undef, 2ncv)
        neupd(ritzvec, howmny, select, d, v, ldv, sigmar, workev,
              bmat, n, which, nev, TOL, resid, ncv, v, ldv,
              iparam, ipntr, workd, workl, lworkl, rwork, info)
        if info[] != 0
            throw(XYEUPD_Exception(info[]))
        end

        p = sortperm(d[1:nev], by=dmap, rev=rev)
        return ritzvec ? (d[p], v[1:n, p],iparam[5],iparam[3],iparam[9],resid) : (d[p],iparam[5],iparam[3],iparam[9],resid)
    elseif sym
        d = Vector{T}(undef, nev)
        sigmar = Ref{T}(sigma)
        seupd(ritzvec, howmny, select, d, v, ldv, sigmar,
              bmat, n, which, nev, TOL, resid, ncv, v, ldv,
              iparam, ipntr, workd, workl, lworkl, info)
        if info[] != 0
            throw(XYEUPD_Exception(info[]))
        end

        p = sortperm(d, by=dmap, rev=rev)
        return ritzvec ? (d[p], v[1:n, p],iparam[5],iparam[3],iparam[9],resid) : (d[p],iparam[5],iparam[3],iparam[9],resid)
    else
        dr = Vector{T}(undef, nev+1)
        di = Vector{T}(undef, nev+1)
        fill!(dr,NaN)
        fill!(di,NaN)
        sigmar = Ref{T}(real(sigma))
        sigmai = Ref{T}(imag(sigma))
        workev = Vector{T}(undef, 3*ncv)
        neupd(ritzvec, howmny, select, dr, di, v, ldv, sigmar, sigmai,
              workev, bmat, n, which, nev, TOL, resid, ncv, v, ldv,
              iparam, ipntr, workd, workl, lworkl, info)
        if info[] != 0
            throw(XYEUPD_Exception(info[]))
        end
        evec = complex.(Matrix{T}(undef, n, nev+1), Matrix{T}(undef, n, nev+1))

        j = 1
        while j <= nev
            if di[j] == 0
                evec[:,j] = v[:,j]
            else # For complex conjugate pairs
                evec[:,j]   = v[:,j] + im*v[:,j+1]
                evec[:,j+1] = v[:,j] - im*v[:,j+1]
                j += 1
            end
            j += 1
        end
        if j == nev+1 && !isnan(di[j])
            if di[j] == 0
                evec[:,j] = v[:,j]
                j += 1
            else
                throw(XYEUPD_Exception(ERR_UNEXPECTED_BEHAVIOR))
            end
        end

        d = complex.(dr, di)

        if j == nev+1
            p = sortperm(d[1:nev], by=dmap, rev=rev)
        else
            p = sortperm(d, by=dmap, rev=rev)
            p = p[1:nev]
        end

        return ritzvec ? (d[p], evec[1:n, p],iparam[5],iparam[3],iparam[9],resid) : (d[p],iparam[5],iparam[3],iparam[9],resid)
    end
end

for (T, saupd_name, seupd_name, naupd_name, neupd_name) in
    ((:Float64, :dsaupd_, :dseupd_, :dnaupd_, :dneupd_),
     (:Float32, :ssaupd_, :sseupd_, :snaupd_, :sneupd_))
    @eval begin
        function naupd(ido, bmat, n, evtype, nev, TOL::Ref{$T}, resid::Vector{$T}, ncv, v::Matrix{$T}, ldv,
                       iparam, ipntr, workd::Vector{$T}, workl::Vector{$T}, lworkl, info)
            ccall(($(string(naupd_name)), libarpack), Cvoid,
                  (Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt},
                   Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                  ido, bmat, n, evtype, nev,
                  TOL, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, info, 1, 2)
        end

        function neupd(rvec, howmny, select, dr, di, z, ldz, sigmar, sigmai,
                  workev::Vector{$T}, bmat, n, evtype, nev, TOL::Ref{$T}, resid::Vector{$T}, ncv, v, ldv,
                  iparam, ipntr, workd::Vector{$T}, workl::Vector{$T}, lworkl, info)
            ccall(($(string(neupd_name)), libarpack), Cvoid,
                  (Ref{BlasInt}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ptr{$T}, Ref{BlasInt},
                   Ref{$T}, Ref{$T}, Ptr{$T}, Ptr{UInt8}, Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt},
                   Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                  rvec, howmny, select, dr, di, z, ldz,
                  sigmar, sigmai, workev, bmat, n, evtype, nev,
                  TOL, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, info, 1, 1, 2)
        end

        function saupd(ido, bmat, n, which, nev, TOL::Ref{$T}, resid::Vector{$T}, ncv, v::Matrix{$T}, ldv,
                       iparam, ipntr, workd::Vector{$T}, workl::Vector{$T}, lworkl, info)
            ccall(($(string(saupd_name)), libarpack), Cvoid,
                  (Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt},
                   Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                  ido, bmat, n, which, nev,
                  TOL, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, info, 1, 2)
        end

        function seupd(rvec, howmny, select, d, z, ldz, sigma,
                       bmat, n, evtype, nev, TOL::Ref{$T}, resid::Vector{$T}, ncv, v::Matrix{$T}, ldv,
                       iparam, ipntr, workd::Vector{$T}, workl::Vector{$T}, lworkl, info)
            ccall(($(string(seupd_name)), libarpack), Cvoid,
                  (Ref{BlasInt}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt},
                   Ptr{$T}, Ptr{UInt8}, Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt},
                   Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                  rvec, howmny, select, d, z, ldz,
                  sigma, bmat, n, evtype, nev,
                  TOL, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, info, 1, 1, 2)
        end
    end
end

for (T, TR, naupd_name, neupd_name) in
    ((:ComplexF64, :Float64, :znaupd_, :zneupd_),
     (:ComplexF32, :Float32, :cnaupd_, :cneupd_))
    @eval begin
        function naupd(ido, bmat, n, evtype, nev, TOL::Ref{$TR}, resid::Vector{$T}, ncv, v::Matrix{$T}, ldv,
                       iparam, ipntr, workd::Vector{$T}, workl::Vector{$T}, lworkl,
                       rwork::Vector{$TR}, info)
            ccall(($(string(naupd_name)), libarpack), Cvoid,
                  (Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt},
                   Ptr{$TR}, Ptr{$T}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ptr{$TR}, Ref{BlasInt}, Clong, Clong),
                  ido, bmat, n, evtype, nev,
                  TOL, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, rwork, info, 1, 2)
        end

        function neupd(rvec, howmny, select, d, z, ldz, sigma, workev::Vector{$T},
                       bmat, n, evtype, nev, TOL::Ref{$TR}, resid::Vector{$T}, ncv, v::Matrix{$T}, ldv,
                       iparam, ipntr, workd::Vector{$T}, workl::Vector{$T}, lworkl,
                       rwork::Vector{$TR}, info)
            ccall(($(string(neupd_name)), libarpack), Cvoid,
                  (Ref{BlasInt}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt},
                   Ptr{$T}, Ptr{$T}, Ptr{UInt8}, Ref{BlasInt}, Ptr{UInt8}, Ref{BlasInt},
                   Ptr{$TR}, Ptr{$T}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$T}, Ptr{$T}, Ref{BlasInt}, Ptr{$TR}, Ref{BlasInt}, Clong, Clong, Clong),
                  rvec, howmny, select, d, z, ldz,
                  sigma, workev, bmat, n, evtype, nev,
                  TOL, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, rwork, info, 1, 1, 2)
        end
    end
end
