using BinaryProvider

# Snarf download_info from `build.jl`
m = Module(:__anon__)
Core.eval(m, quote
	using BinaryProvider
	function write_deps_file(args...; kwargs...); end
    function install(args...; kwargs...); end
	ARGS = ["."]
end)
include_string(m, String(read("build.jl")))
download_info = Core.eval(m, :(download_info))

try
    mkdir("downloads")
catch
end

# Download them all!
for key in keys(download_info)
    download_verify_unpack(
        download_info[key][1],
        download_info[key][2],
        joinpath(@__DIR__, "downloads", triplet(key));
        tarball_path=joinpath(@__DIR__, "downloads", triplet(key)*".tar.gz"),
    )
end

# Get linkage name of openblas for everybody
using ObjectFile
println("default_blaslib_names = Dict(")
for p in keys(download_info)
    libarpack = LibraryProduct(Prefix(joinpath(@__DIR__, "downloads", triplet(p))), "libarpack", :libarpack)
    libarpack_path = locate(libarpack; platform=p)
    blas_name = readmeta(libarpack_path) do oh
        blaslibs = filter(x -> occursin("blas", x), [path(l) for l in DynamicLinks(oh)])
        println("    $(repr(p)) => $(repr(basename(first(blaslibs)))),")
    end
end
println(")")
