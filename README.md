# Arpack

[![Build Status](https://travis-ci.org/JuliaLinearAlgebra/Arpack.jl.svg?branch=master)](https://travis-ci.org/JuliaLinearAlgebra/Arpack.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/v6icbqh1xq5y7261?svg=true)](https://ci.appveyor.com/project/andreasnoack/arpack-jl)
[![Coverage Status](https://coveralls.io/repos/JuliaLinearAlgebra/Arpack.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaLinearAlgebra/Arpack.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLinearAlgebra/Arpack.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLinearAlgebra/Arpack.jl?branch=master)
[![][docs-stable-img]][docs-stable-url]
[![][docs-latest-img]][docs-latest-url]

Julia wrapper for the [arpack](https://github.com/opencollab/arpack-ng/) library
designed to solve large scale eigenvalue problems.

## Installation

You can install Arpack.jl through the Julia package manager:
```julia
julia> Pkg.add("Arpack")
```

Arpack.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl)
to automatically install the Arpack binaries.

## Custom Installation

If you get
```
ERROR: LoadError: LibraryProduct(nothing, ["libarpack"], :libarpack, "Prefix(~/.julia/packages/Arpack/cu5By/deps/usr)") is not satisfied, cannot generate deps.jl!
```
when building Arpack, it may be because your Julia installation uses your system
blas for which the symbols are not suffixed by `_64_` while this is required by
the compiled binaries provided by [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl).
This is notably the case on [ArchLinux](https://www.archlinux.org/).
In these case, compile binaries that do not require this suffix as follows.
Download the source of the [v3.5.0 of arpack-ng](https://github.com/opencollab/arpack-ng/releases/tag/3.5.0),
extract it in some `<directory>`, build it and do (note that you may need to
update `cu6By` to match the one printed in the error message printed above).
```
$ cp <directory>/arpack-ng-3.5.0/SRC/.libs/libarpack.so.2.0.0 ~/.julia/packages/Arpack/cu5By/deps/usr/lib/
$ julia -e 'import Pkg; Pkg.build("Arpack")'
  Building Arpack → `~/.julia/packages/Arpack/UiiMc/deps/build.log`
```

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://JuliaLinearAlgebra.github.io/Arpack.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://JuliaLinearAlgebra.github.io/Arpack.jl/stable/

As posted [here](https://discourse.julialang.org/t/resolved-arpack-deps-jl-linux-build-issue/30923) for Deb linux releases there is a less intrusive resolution that only requires symlinking to the distribution library of openblas. For Ubuntu flavours:

```
$ sudo ln -s /usr/lib/x86_64-linux-gnu/libopenblas.so /usr/lib/libopenblas64_.so.0
```
