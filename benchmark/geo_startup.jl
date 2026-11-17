#!/usr/bin/env julia
# Run each case in a fresh process: julia --project=. --startup-file=no -O2 benchmark/geo_startup.jl transfinite
const USAGE = "Usage: geo_startup.jl transfinite|unstructured [output.toml]"
if ARGS == ["--help"]
    println(USAGE)
    exit()
end
length(ARGS) in (1,2) && ARGS[1] in ("transfinite","unstructured") || error(USAGE)
const OPTIONS = Base.JLOptions()
OPTIONS.compile_enabled == 1 && OPTIONS.opt_level >= 2 &&
    OPTIONS.code_coverage == 0 && OPTIONS.malloc_log == 0 ||
    error("Use --compile=yes -O2 without coverage or allocation instrumentation for timing")
Threads.nthreads() == 1 || error("Run with --threads=1")
using SHA, TOML
function source_hash(root)
    files = sort!([joinpath(dir,file) for (dir,_,files) in walkdir(joinpath(root,"src"))
                   for file in files if endswith(file,".jl")])
    return bytes2hex(sha256(join(relpath(file,root)*":"*bytes2hex(sha256(read(file)))
                                for file in files)))
end
const ROOT = dirname(@__DIR__)
const SOURCE_HASH = source_hash(ROOT)
const PACKAGE_ID = Base.identify_package("Tessella")
PACKAGE_ID === nothing && error("Activate this repository's project with --project=.")
haskey(Base.loaded_modules,PACKAGE_ID) && error("Run the benchmark in a fresh Julia process")
const CACHE_AVAILABLE = Base.isprecompiled(PACKAGE_ID;ignore_loaded=true)
const LOAD_STARTED = time_ns()
using Tessella
const PACKAGE_LOAD_SECONDS = (time_ns()-LOAD_STARTED)/1e9
realpath(pkgdir(Tessella)) == realpath(ROOT) || error("Activate this repository's project")
source_hash(ROOT) == SOURCE_HASH || error("Source changed during package loading")

function mesh_fingerprint(execution)
    isempty(execution.errors) || error(join(execution.errors,"\n"))
    return (coordinates=copy(execution.mesh.coords),
            crc=Tessella.MeshTypes.mesh_crc(execution.mesh).sha)
end
function check_case(execution,kind)
    mesh = execution.mesh
    Tessella.MeshTypes.validate(mesh).ok || error("Invalid mesh")
    all(iszero,mesh.coords[3,:]) && all(x->-64eps(Float64)<=x<=1+64eps(Float64),mesh.coords[1:2,:]) ||
        error("Mesh vertices lie outside the unit square")
    areas = [((mesh.coords[1,b]-mesh.coords[1,a])*(mesh.coords[2,c]-mesh.coords[2,a])-
              (mesh.coords[2,b]-mesh.coords[2,a])*(mesh.coords[1,c]-mesh.coords[1,a]))/2
             for (a,b,c) in eachcol(mesh.tris)]
    all(>(0),areas) && isapprox(sum(areas),1.0;rtol=64eps(Float64)) ||
        error("Mesh does not cover the unit square with positive triangles")
    kind != "transfinite" || (size(mesh.coords,2),size(mesh.tris,2)) == (25,32) ||
        error("Unexpected transfinite mesh size")
    return mesh_fingerprint(execution)
end
function run_case(kind)
    geometry = """
    Point(1)={0,0,0,0.25}; Point(2)={1,0,0,0.25};
    Point(3)={1,1,0,0.25}; Point(4)={0,1,0,0.25};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    """
    kind == "transfinite" && (geometry *= "Transfinite Curve {1,2,3,4}=5; Transfinite Surface {1};\n")
    mktemp() do path,io
        write(io,geometry); close(io)
        first_mesh = @timed Base.inferencebarrier(Tessella.execute_geo)(path;mesh_dim=2)
        println(stderr,"First mesh returned in ",first_mesh.time," seconds")
        execution = first_mesh.value; mesh = execution.mesh
        fingerprint = Base.inferencebarrier(check_case)(execution,kind)
        warm_times = Float64[]
        for _ in 1:7
            warm = @timed Base.inferencebarrier(Tessella.execute_geo)(path;mesh_dim=2)
            isequal(Base.inferencebarrier(mesh_fingerprint)(warm.value),fingerprint) ||
                error("Repeated execution changed the mesh")
            push!(warm_times,warm.time)
        end
        source_hash(ROOT) == SOURCE_HASH || error("Source changed during the benchmark")
        return Dict("case"=>kind,"julia_version"=>string(VERSION),
            "tessella_version"=>string(Base.pkgversion(Tessella)),"source_tree_sha256"=>SOURCE_HASH,
            "compatible_cache_before_load"=>CACHE_AVAILABLE,"package_load_seconds"=>PACKAGE_LOAD_SECONDS,
            "first_mesh_seconds"=>first_mesh.time,"first_mesh_compile_seconds"=>first_mesh.compile_time,
            "first_mesh_bytes"=>first_mesh.bytes,"first_mesh_gc_seconds"=>first_mesh.gctime,
            "warm_median_seconds"=>sort!(warm_times)[4],"nodes"=>size(mesh.coords,2),
            "triangles"=>size(mesh.tris,2),"mesh_crc"=>fingerprint.crc,
            "optimization_level"=>Int(OPTIONS.opt_level),"check_bounds"=>Int(OPTIONS.check_bounds),
            "compiled_modules_mode"=>Int(OPTIONS.use_compiled_modules))
    end
end
result = run_case(ARGS[1])
length(ARGS) == 1 ? TOML.print(stdout,result;sorted=true) :
    open(io->TOML.print(io,result;sorted=true),ARGS[2],"w")
