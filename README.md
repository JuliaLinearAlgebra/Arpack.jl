# Arpack

[![Build Status](https://travis-ci.org/JuliaLinearAlgebra/Arpack.jl.svg?branch=master)](https://travis-ci.org/JuliaLinearAlgebra/Arpack.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/v6icbqh1xq5y7261?svg=true)](https://ci.appveyor.com/project/andreasnoack/arpack-jl)
[![Coverage Status](https://coveralls.io/repos/JuliaLinearAlgebra/Arpack.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaLinearAlgebra/Arpack.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLinearAlgebra/Arpack.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLinearAlgebra/Arpack.jl?branch=master)
[![][docs-stable-img]][docs-stable-url]
[![][docs-latest-img]][docs-latest-url]

If you get
```
ERROR: LoadError: LibraryProduct(nothing, ["libarpack"], :libarpack, "Prefix(~/.julia/packages/Arpack/cu5By/deps/usr)") is not satisfied, cannot generate deps.jl!
```
when building Arpack, download the source of the [v3.5.0 of arpack-ng](https://github.com/opencollab/arpack-ng/releases/tag/3.5.0),
extract it in some `<directory>`, build it and do
```
$ cp <directory>/arpack-ng-3.5.0/SRC/.libs/libarpack.so.2.0.0 ~/.julia/packages/Arpack/cu5By/deps/usr/lib/
$ julia -e 'import Pkg; Pkg.build("Arpack")'
  Building Arpack → `~/.julia/packages/Arpack/UiiMc/deps/build.log`
```

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://JuliaLinearAlgebra.github.io/Arpack.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://JuliaLinearAlgebra.github.io/Arpack.jl/stable/
