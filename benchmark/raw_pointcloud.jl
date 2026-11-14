#!/usr/bin/env julia
# Independent benchmark driver; Gmsh is an optional external executable library.
const DRIVER_STARTED = time_ns()
using Libdl, LinearAlgebra, Random, SHA, Statistics, TOML
const COLD = length(ARGS) >= 2 && ARGS[1] == "--cold" ? ARGS[2] : ""
const LATENCY = !isempty(ARGS) && ARGS[1] == "--latency"
const CHECK_DRIVER = ARGS == ["--check-driver"]
const USAGE="Usage: raw_pointcloud.jl [output.toml] | --check-driver | --cold|--latency Gmsh|Tessella 2|3 points output.toml"
if !isempty(ARGS) && ARGS[1] in ("--cold","--latency")
    length(ARGS) == 5 && ARGS[2] in ("Gmsh","Tessella") || error(USAGE)
    case_dimension=tryparse(Int,ARGS[3]);case_points=tryparse(Int,ARGS[4])
    case_dimension in (2,3) && case_points !== nothing &&
        case_dimension+1 <= case_points <= typemax(Int32) || error(USAGE)
elseif !CHECK_DRIVER
    length(ARGS) <= 1 && (isempty(ARGS) || !startswith(ARGS[1],"--")) || error(USAGE)
end
if !CHECK_DRIVER && !LATENCY
    options=Base.JLOptions()
    options.compile_enabled == 1 && options.opt_level >= 2 &&
        options.code_coverage == 0 && options.malloc_log == 0 ||
        error("Use --compile=yes -O2 without coverage or allocation instrumentation for timing")
end
const NEED_GMSH = !LATENCY && COLD != "Tessella"

const LIBGMSH = !NEED_GMSH ? "" : haskey(ENV, "GMSH_LIBRARY") ? ENV["GMSH_LIBRARY"] : Libdl.find_library(
    ["gmsh", "libgmsh", "libgmsh.4.15.2"],
    ["/opt/homebrew/opt/gmsh/lib", "/usr/local/lib"])
NEED_GMSH && isempty(LIBGMSH) && error("Set GMSH_LIBRARY to the Gmsh 4.15.2 shared library")
const HANDLE = NEED_GMSH ? Libdl.dlopen(LIBGMSH) : C_NULL
const SYMBOLS = NEED_GMSH ? Dict(name => Libdl.dlsym(HANDLE, name) for name in
    (:gmshInitialize, :gmshFinalize, :gmshFree, :gmshOptionSetNumber,
     :gmshOptionGetNumber, :gmshOptionGetString, :gmshLoggerGetLastError,
     :gmshAlgorithmTriangulate, :gmshAlgorithmTetrahedralize)) : Dict{Symbol,Ptr{Cvoid}}()

free_gmsh(p) = p == C_NULL ? nothing : ccall(SYMBOLS[:gmshFree], Cvoid, (Ptr{Cvoid},), p)
function gmsh_string(symbol, name="")
    p = Ref{Ptr{Cchar}}(C_NULL); ierr = Ref{Cint}(0)
    try
        if symbol == :gmshOptionGetString
            ccall(SYMBOLS[symbol], Cvoid, (Cstring, Ref{Ptr{Cchar}}, Ref{Cint}), name, p, ierr)
        else
            ccall(SYMBOLS[symbol], Cvoid, (Ref{Ptr{Cchar}}, Ref{Cint}), p, ierr)
        end
        ierr[] == 0 || error("Gmsh string query failed: $name")
        return p[] == C_NULL ? "" : unsafe_string(p[])
    finally
        free_gmsh(p[])
    end
end
function check_gmsh(ierr)
    ierr[] == 0 || error("Gmsh error $(ierr[]): $(gmsh_string(:gmshLoggerGetLastError))")
end
function set_option(name, value)
    ierr = Ref{Cint}(0)
    ccall(SYMBOLS[:gmshOptionSetNumber], Cvoid, (Cstring, Cdouble, Ref{Cint}), name, value, ierr)
    check_gmsh(ierr)
end
function get_option(name)
    value = Ref{Cdouble}(0); ierr = Ref{Cint}(0)
    ccall(SYMBOLS[:gmshOptionGetNumber], Cvoid, (Cstring, Ref{Cdouble}, Ref{Cint}), name, value, ierr)
    check_gmsh(ierr)
    return value[]
end
function initialize_gmsh()
    ierr = Ref{Cint}(0)
    ccall(SYMBOLS[:gmshInitialize], Cvoid,
          (Cint, Ptr{Ptr{Cchar}}, Cint, Cint, Ref{Cint}), 0, C_NULL, 0, 0, ierr)
    check_gmsh(ierr)
    for name in ("General.NumThreads", "Mesh.MaxNumThreads2D", "Mesh.MaxNumThreads3D")
        set_option(name, 1)
    end
    set_option("General.Terminal", 0)
    set_option("Mesh.RandomSeed", 1)
end
function finalize_gmsh()
    ierr = Ref{Cint}(0)
    ccall(SYMBOLS[:gmshFinalize], Cvoid, (Ref{Cint},), ierr)
    check_gmsh(ierr)
end

# Both adapters accept D×N Float64 coordinates and return owned 3×N Float64
# coordinates and (D+1)×M Int32 cells. Adapter allocations are timed.
const EMPTY_CONSTRAINTS = Csize_t[]
function gmsh_mesh(points::Matrix{Float64}, boundary=nothing)
    dim, n = size(points)
    constraint_ids = boundary === nothing ? EMPTY_CONSTRAINTS : Csize_t.(vec(boundary))
    cells = Ref{Ptr{Csize_t}}(C_NULL); ncells = Ref{Csize_t}(0)
    steiner = Ref{Ptr{Cdouble}}(C_NULL); nsteiner = Ref{Csize_t}(0)
    ierr = Ref{Cint}(0)
    try
        if dim == 2
            ccall(SYMBOLS[:gmshAlgorithmTriangulate], Cvoid,
                  (Ptr{Cdouble}, Csize_t, Ref{Ptr{Csize_t}}, Ref{Csize_t}, Ptr{Csize_t}, Csize_t, Ref{Cint}),
                  points, length(points), cells, ncells, constraint_ids, length(constraint_ids), ierr)
        elseif dim == 3
            ccall(SYMBOLS[:gmshAlgorithmTetrahedralize], Cvoid,
                  (Ptr{Cdouble}, Csize_t, Ref{Ptr{Csize_t}}, Ref{Csize_t}, Ref{Ptr{Cdouble}}, Ref{Csize_t}, Ptr{Csize_t}, Csize_t, Ref{Cint}),
                  points, length(points), cells, ncells, steiner, nsteiner, constraint_ids, length(constraint_ids), ierr)
        else
            error("Only 2D and 3D point clouds are supported")
        end
        check_gmsh(ierr)
        ncells[] % (dim + 1) == 0 || error("Malformed Gmsh connectivity")
        nsteiner[] % 3 == 0 || error("Malformed Gmsh Steiner coordinates")
        coords = zeros(3, n + Int(nsteiner[] ÷ 3))
        coords[1:dim, 1:n] = points
        if nsteiner[] > 0
            coords[:, n+1:end] = reshape(unsafe_wrap(Array, steiner[], Int(nsteiner[])), 3, :)
        end
        connectivity = Matrix{Int32}(undef, dim + 1, Int(ncells[] ÷ (dim + 1)))
        if ncells[] > 0
            connectivity[:] = unsafe_wrap(Array, cells[], Int(ncells[]))
        end
        return coords, connectivity
    finally
        free_gmsh(cells[]); free_gmsh(steiner[])
    end
end

# Exact geometric oracle, independent of both meshers' predicates/topology.
const Exact = Rational{BigInt}
function exactdet(input)
    a = copy(input); n = size(a, 1); value = one(Exact)
    for k in 1:n
        pivot = findfirst(i -> a[i, k] != 0, k:n)
        pivot === nothing && return zero(Exact)
        row = k + pivot - 1
        if row != k
            a[k, :], a[row, :] = copy(a[row, :]), copy(a[k, :])
            value = -value
        end
        p = a[k, k]; value *= p
        for i in k+1:n
            factor = a[i, k] / p
            for j in k+1:n
                a[i, j] -= factor * a[k, j]
            end
        end
    end
    return value
end
function simplex_det(p, ids)
    d = size(p, 1)
    exactdet([p[i, ids[j+1]] - p[i, ids[1]] for i in 1:d, j in 1:d])
end
function sphere(p, ids)
    d = size(p, 1)
    a = [2 * (p[j, ids[i+1]] - p[j, ids[1]]) for i in 1:d, j in 1:d]
    b = [sum(p[j, ids[i+1]]^2 - p[j, ids[1]]^2 for j in 1:d) for i in 1:d]
    denominator = exactdet(a)
    center = Exact[]
    for k in 1:d
        replaced = copy(a); replaced[:, k] = b
        push!(center, exactdet(replaced) / denominator)
    end
    return center, sum((center[j] - p[j, ids[1]])^2 for j in 1:d)
end
# Exact determinant in common integer coordinate units; no mesher predicates.
function integer_simplex_det(p, ids)
    a = p[:,ids[2]] - p[:,ids[1]]
    b = p[:,ids[3]] - p[:,ids[1]]
    size(p,1) == 2 && return a[1]*b[2] - a[2]*b[1]
    c = p[:,ids[4]] - p[:,ids[1]]
    return a[1]*(b[2]*c[3]-b[3]*c[2]) - a[2]*(b[1]*c[3]-b[3]*c[1]) +
           a[3]*(b[1]*c[2]-b[2]*c[1])
end
struct OracleUnsupported <: Exception
    message::String
end
Base.showerror(io::IO,err::OracleUnsupported)=print(io,err.message)
function validate_complete(points, coords, cells)
    structural_check(points,coords,cells)
    length(Set(Tuple(points[:,i]) for i in axes(points,2))) == size(coords,2) ||
        throw(OracleUnsupported("Steiner output requires an independent input-hull containment oracle"))
    d = size(points,1)
    size(cells,2) > 0 || error("Empty full-dimensional mesh")
    rational = Rational{BigInt}.(coords[1:d,:])
    denominator_all = maximum(denominator,rational)
    p = map(x->numerator(x)*(denominator_all÷denominator(x)),rational)
    lookup = Dict(Tuple(coords[1:d,i])=>i for i in axes(coords,2))
    originals = [lookup[Tuple(points[:,i])] for i in axes(points,2)]
    faces = Dict{Tuple,Vector{Int}}(); seen = Set{Tuple}(); measure = BigInt(0)
    for (index,cell) in enumerate(eachcol(cells))
        ids = Tuple(Int.(cell)); key = Tuple(sort(collect(ids)))
        length(Set(ids)) == d+1 || error("Repeated simplex vertex")
        key in seen && error("Duplicate simplex")
        push!(seen,key)
        determinant = integer_simplex_det(p,ids)
        determinant != 0 || error("Zero-measure simplex")
        measure += abs(determinant)
        for omitted in 1:d+1
            face = Tuple(sort([ids[i] for i in 1:d+1 if i != omitted]))
            push!(get!(Vector{Int},faces,face),index)
        end
    end
    all(v->length(v) in (1,2),values(faces)) || error("Non-manifold simplex incidence")
    parent = collect(1:size(cells,2))
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]; i = parent[i]
        end
        i
    end
    for (face,adjacent) in faces
        length(adjacent) == 2 || continue
        opposite=[only(i for i in cells[:,cell] if i ∉ face) for cell in adjacent]
        integer_simplex_det(p,(face...,opposite[1])) *
            integer_simplex_det(p,(face...,opposite[2])) < 0 ||
            error("Incident simplices lie on the same side of a shared face")
        parent[root(adjacent[1])] = root(adjacent[2])
    end
    length(Set(root(i) for i in eachindex(parent))) == 1 || error("Disconnected cell complex")
    boundary = [face for (face,adjacent) in faces if length(adjacent) == 1]
    ridges = Dict{Tuple,Int}(); hull_measure = BigInt(0)
    for face in boundary
        signs = [sign(integer_simplex_det(p,(face...,i))) for i in axes(p,2)]
        all(>=(0),signs) || all(<=(0),signs) || error("Boundary face is not a convex-hull supporting plane")
        hull_measure += abs(integer_simplex_det(p,(originals[1],face...)))
        for omitted in 1:d
            ridge = Tuple(sort([face[i] for i in 1:d if i != omitted]))
            ridges[ridge] = get(ridges,ridge,0)+1
        end
    end
    !isempty(boundary) && all(==(2),values(ridges)) || error("Boundary is not closed")
    measure == hull_measure || error("Simplex measure differs from convex-hull measure")
    return Dict("status"=>"pass","method"=>"exact integer convex-hull support, measure and manifold incidence",
        "boundary_faces"=>length(boundary),"cells"=>size(cells,2))
end

function validate_delaunay(points,coords,cells)
    p=Exact.(coords[1:size(points,1),:]);checks=0
    for cell in eachcol(cells)
        ids=Tuple(Int.(cell));center,radius_squared=sphere(p,ids)
        for i in axes(p,2)
            i in ids && continue
            sum((p[j,i]-center[j])^2 for j in axes(p,1)) >= radius_squared ||
                error("Point $i lies strictly inside a simplex circumsphere")
            checks+=1
        end
    end
    return Dict("status"=>"pass","method"=>"rational empty circumspheres","sphere_checks"=>checks)
end
function validate_exact(points,coords,cells)
    completion=validate_complete(points,coords,cells)
    delaunay=validate_delaunay(points,coords,cells)
    merge(completion,delaunay,Dict("method"=>"exact integer completion and rational empty circumspheres"))
end

function tessella_mesh(points::Matrix{Float64}, boundary=nothing)
    if size(points, 1) == 2
        mesh = if boundary === nothing
            Tessella.Mesh2D.triangulate(points[1, :], points[2, :]; rng_seed=1)
        else
            triangulation = Tessella.Mesh2D.constrained_delaunay(points[1, :], points[2, :],
                [Tuple(edge) for edge in eachcol(boundary)]; rng_seed=1)
            Tessella.Mesh2D.to_mesh(triangulation;
                interior=Tessella.Mesh2D.classify_interior(triangulation))
        end
        return mesh.coords, mesh.tris
    end
    # This mode preserves input geometry; default perturb=true is a different task.
    mesh = Tessella.Mesh3D.to_mesh3(Tessella.Mesh3D.delaunay3d(
        points[1, :], points[2, :], points[3, :]; rng_seed=1, perturb=false))
    if boundary !== nothing
        for face in eachcol(boundary)
            vertices=map(i->Tuple(points[:,i]),face)
            if !Tessella.Mesh3D.mesh_covers_triangle3(mesh,vertices...)
                mesh=Tessella.Mesh3D.recover_triangle3(mesh,vertices...)
            end
        end
    end
    return mesh.coords, mesh.tets
end

function cloud(d, n, seed, family)
    rng = MersenneTwister(seed)
    points = rand(rng, d, n)
    if family == "simplex_enclosed"
        points[:, 1:d+1] = hcat(zeros(d), Matrix{Float64}(I, d, d))
        for i in d+2:n
            cuts = sort(rand(rng, d))
            weights = diff(vcat(0.0, cuts, 1.0))
            points[:, i] = 0.8 .* weights[1:d] .+ 0.2 / (d+1)
        end
    elseif family == "thin"
        points[d, :] .*= 1e-10
    elseif family == "small_scale"
        points .*= 1e-9
    elseif family == "large_scale"
        points .*= 1e9
    elseif family == "translated"
        points .+= 1e9
    elseif family != "uniform"
        error("Unknown point family $family")
    end
    return points
end

function structural_check(points, coords, cells)
    d = size(points, 1)
    size(coords, 1) == 3 && size(cells, 1) == d+1 || error("Wrong output dimensions")
    all(isfinite, coords) || error("Non-finite output coordinates")
    d == 2 && any(!=(0), coords[3, :]) && error("Nonzero output z coordinate")
    all(i -> 1 <= i <= size(coords, 2), cells) || error("Invalid node index")
    Set(Int.(vec(cells))) == Set(axes(coords, 2)) || error("Unused output node")
    input = Set(Tuple(points[:, i]) for i in axes(points, 2))
    output = Set(Tuple(coords[1:d, i]) for i in axes(coords, 2))
    issubset(input, output) || error("Input coordinates changed or disappeared")
    length(output) == size(coords, 2) || error("Duplicate output coordinates")
end
function output_hash(coords, cells)
    io = IOBuffer()
    write(io, vec(coords))
    canonical = sort([Tuple(sort(collect(cell))) for cell in eachcol(cells)])
    for cell in canonical
        write(io, collect(cell))
    end
    bytes2hex(sha256(take!(io)))
end

simplex_boundary(d) = d == 2 ? Int32[1 2 3;2 3 1] : Int32[1 1 1 2;3 2 4 3;2 4 3 4]

function paired_measure(points, repetitions, order_rng; boundary=nothing, adapters=nothing,record_geometry=false)
    # Clear preceding exact-oracle garbage outside the meshing measurements.
    GC.gc()
    adapters = adapters === nothing ? (("Gmsh", p->gmsh_mesh(p,boundary)),
                ("Tessella", p->tessella_mesh(p,boundary))) : adapters
    samples = [Float64[], Float64[]]
    outputs = [Tuple{Matrix{Float64},Matrix{Int32}}[], Tuple{Matrix{Float64},Matrix{Int32}}[]]
    errors = [String[], String[]]
    first_calls = zeros(2)
    orders = Vector{String}[]
    # Prime both compiled adapters, retaining these outputs for validation too.
    for i in randperm(order_rng, 2)
        try
            sample = @timed adapters[i][2](points)
            first_calls[i] = sample.time
            push!(outputs[i], sample.value)
        catch err
            err isa InterruptException && rethrow()
            push!(errors[i], sprint(showerror, err))
        end
    end
    for _ in 1:repetitions
        order = randperm(order_rng, 2)
        push!(orders, [adapters[i][1] for i in order])
        for i in order
            try
                sample = @timed adapters[i][2](points)
                push!(samples[i], sample.time)
                push!(outputs[i], sample.value)
            catch err
                err isa InterruptException && rethrow()
                push!(errors[i], sprint(showerror, err))
            end
        end
    end
    rows = Dict{String,Any}[]
    for i in 1:2
        row = Dict{String,Any}("implementation" => adapters[i][1],
            "first_in_case_seconds" => first_calls[i], "warm_samples_seconds" => samples[i],
            "measurement_orders" => orders, "errors" => errors[i])
        if !isempty(samples[i])
            row["warm_median_seconds"] = median(samples[i])
            row["warm_p25_seconds"] = quantile(samples[i], 0.25)
            row["warm_p75_seconds"] = quantile(samples[i], 0.75)
        end
        checked=Set{String}();failures=String[];unsupported=String[];delaunay_failures=String[]
        for (coords, cells) in outputs[i]
            fingerprint = output_hash(coords, cells)
            fingerprint in checked && continue
            push!(checked, fingerprint)
            row["output_points"] = size(coords, 2)
            row["cells"] = size(cells, 2)
            try
                row["completion_check"]=validate_complete(points,coords,cells)
            catch err
                err isa InterruptException && rethrow()
                push!(err isa OracleUnsupported ? unsupported : failures,sprint(showerror,err))
                continue
            end
            if size(points,2) <= 16
                try
                    row["delaunay_check"]=validate_delaunay(points,coords,cells)
                catch err
                    err isa InterruptException && rethrow()
                    push!(delaunay_failures,sprint(showerror,err))
                end
            end
        end
        row["distinct_outputs_checked"] = length(checked)
        row["validation_failures"] = unique(failures)
        row["oracle_limitations"] = unique(unsupported)
        row["delaunay_failures"]=unique(delaunay_failures)
        row["completion_status"]=!isempty(failures) ? "fail" : !isempty(unsupported) ? "unsupported" :
            haskey(row,"completion_check") ? "pass" : "not_checked"
        row["delaunay_status"]=!isempty(delaunay_failures) ? "fail" :
            haskey(row,"delaunay_check") ? "pass" : "not_checked"
        row["validation_scope"] = size(points,2) <= 16 ? "exact_delaunay_and_completion" : "exact_completion_only"
        row["status"] = !isempty(errors[i]) ? "error" : !isempty(failures) ? "returned_invalid" :
            !isempty(unsupported) ? "oracle_unsupported" : !isempty(delaunay_failures) ? "complete_mesh_non_delaunay" :
            size(points,2) <= 16 ? "pass" : "complete_mesh_pass"
        if record_geometry && !isempty(outputs[i])
            coords,cells=outputs[i][1]
            row["geometry"]=Dict("sample"=>"priming call",
                "input_coordinates"=>[collect(column) for column in eachcol(points)],
                "output_coordinates"=>[collect(column) for column in eachcol(coords)],
                "cells"=>[Int.(collect(column)) for column in eachcol(cells)])
        end
        push!(rows, row)
    end
    rows
end

function benchmark_cases()
    cases = vec([(d,n,seed,family) for d in (2,3), n in (16,64,256,1024),
                 seed in (9,1729,2026), family in ("uniform","simplex_enclosed")])
    append!(cases, vec([(d,16,9,family) for d in (2,3),
                        family in ("thin","small_scale","large_scale","translated")]))
    return cases
end
function source_fingerprint()
    root = dirname(@__DIR__)
    paths = sort([joinpath(dir,file) for (dir,_,files) in walkdir(joinpath(root,"src"))
                   for file in files if endswith(file,".jl")])
    hashes = Dict(replace(relpath(path,root),'\\'=>'/') => bytes2hex(sha256(read(path))) for path in paths)
    io = IOBuffer()
    for path in sort(collect(keys(hashes)))
        print(io,path,'\0',hashes[path],'\n')
    end
    return hashes,bytes2hex(sha256(take!(io)))
end
const SOURCE_FILES_SHA256, SOURCE_TREE_SHA256 = source_fingerprint()
const CASES = benchmark_cases()
if isempty(COLD) && !LATENCY
    CASES isa Vector && length(CASES) == 56 || error("Invalid benchmark case vector")
    all(c -> size(cloud(c...)) == (c[1],c[2]), CASES) || error("Invalid cloud dimensions")
end

function metadata()
    root = dirname(@__DIR__)
    return Dict{String,Any}(
        "julia_version" => string(VERSION), "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(), "cpu" => Sys.CPU_NAME,
        "machine" => Sys.MACHINE, "logical_cpu_threads" => Sys.CPU_THREADS,
        "tessella_version" => TOML.parsefile(joinpath(root, "Project.toml"))["version"],
        "git_head" => readchomp(`git -C $root rev-parse HEAD`),
        "git_dirty" => !isempty(readchomp(`git -C $root status --porcelain`)),
        "production_source_dirty" => !isempty(readchomp(`git -C $root status --porcelain -- src`)),
        "source_files_sha256" => SOURCE_FILES_SHA256,
        "kernel_source_sha256" => Dict(path=>SOURCE_FILES_SHA256[path] for path in
            ("src/core/Predicates.jl","src/core/MeshTypes.jl","src/core/GmshLibm.jl",
             "src/meshing/Mesh2D.jl","src/meshing/Mesh3D.jl")),
        "source_tree_sha256" => SOURCE_TREE_SHA256,
        "source_changed_during_load" => source_fingerprint()[2] != SOURCE_TREE_SHA256,
        "driver_sha256" => DRIVER_SHA256,
        "julia_compiled_modules_mode" => Int(Base.JLOptions().use_compiled_modules),
        "julia_compile_enabled" => Int(Base.JLOptions().compile_enabled),
        "julia_check_bounds" => Int(Base.JLOptions().check_bounds),
        "julia_optimization_level" => Int(Base.JLOptions().opt_level),
        "gmsh_default_algorithm_3d" => get_option("Mesh.Algorithm3D"),
        "package_load_seconds" => PACKAGE_LOAD_SECONDS,
        "host_note" => get(ENV, "BENCHMARK_HOST_NOTE", "Shared host; other Julia processes may contend for CPU and memory"),
        "gmsh_version" => gmsh_string(:gmshOptionGetString, "General.Version"),
        "gmsh_build_info" => gmsh_string(:gmshOptionGetString, "General.BuildInfo"),
        "gmsh_library" => Libdl.dlpath(HANDLE),
        "gmsh_library_sha256" => bytes2hex(sha256(read(Libdl.dlpath(HANDLE)))),
        "gmsh_random_factor" => get_option("Mesh.RandomFactor"),
        "gmsh_random_factor_3d" => get_option("Mesh.RandomFactor3D"),
        "gmsh_random_seed" => get_option("Mesh.RandomSeed"),
        "gmsh_threads" => get_option("General.NumThreads"),
        "gmsh_max_threads_2d" => get_option("Mesh.MaxNumThreads2D"),
        "gmsh_max_threads_3d" => get_option("Mesh.MaxNumThreads3D"),
        "omp_num_threads" => get(ENV,"OMP_NUM_THREADS","unset"),
        "tessella_3d_constraints_policy" => "verify coverage; recover uncovered faces",
        "tessella_3d_perturb" => false, "algorithm_seed" => 1,
        "input_seeds" => [9, 1729, 2026], "order_seed" => 8421,
        "repetitions" => 9)
end
function save(output, data)
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        TOML.print(io, data; sorted=true)
    end
end
function cold_run()
    length(ARGS) == 5 || error("Usage: --cold Gmsh|Tessella dimension points output.toml")
    implementation = ARGS[2]; d = parse(Int, ARGS[3]); n = parse(Int, ARGS[4])
    implementation in ("Gmsh", "Tessella") || error("Unknown implementation")
    points = cloud(d, n, 9, "simplex_enclosed")
    initialization_seconds = implementation == "Gmsh" ? @elapsed(initialize_gmsh()) : 0.0
    try
        implementation == "Gmsh" && d == 3 && set_option("Mesh.Algorithm3D",10)
        boundary=simplex_boundary(d)
        f = implementation == "Gmsh" ? gmsh_mesh : tessella_mesh
        sample = @timed f(points,boundary)
        elapsed = (time_ns() - DRIVER_STARTED) / 1e9
        println("MESH_RETURNED");flush(stdout)
        completion=validate_complete(points,sample.value...)
        result=Dict{String,Any}("implementation" => implementation, "dimension" => d,
            "points" => n, "family" => "simplex_enclosed", "julia_version" => string(VERSION),
            "explicit_hull_constraints" => true,"completion_check"=>completion,
            "first_mesh_seconds" => sample.time, "package_load_seconds" => PACKAGE_LOAD_SECONDS,
            "gmsh_initialization_seconds" => initialization_seconds,
            "julia_compiled_modules_mode" => Int(Base.JLOptions().use_compiled_modules),
            "julia_compile_enabled" => Int(Base.JLOptions().compile_enabled),
            "julia_check_bounds" => Int(Base.JLOptions().check_bounds),
            "julia_optimization_level" => Int(Base.JLOptions().opt_level),
            "source_files_sha256" => SOURCE_FILES_SHA256,
            "source_tree_sha256" => SOURCE_TREE_SHA256,
            "driver_sha256" => DRIVER_SHA256,
            "source_changed_during_load" => source_fingerprint()[2] != SOURCE_TREE_SHA256,
            "host_note" => get(ENV,"BENCHMARK_HOST_NOTE","Shared host; other Julia processes may contend for CPU and memory"),
            "driver_through_mesh_seconds" => elapsed,
            "input_sha256" => bytes2hex(sha256(reinterpret(UInt8, vec(points)))))
        if implementation == "Gmsh"
            result["gmsh_version"]=gmsh_string(:gmshOptionGetString,"General.Version")
            result["gmsh_library_sha256"]=bytes2hex(sha256(read(Libdl.dlpath(HANDLE))))
            result["gmsh_algorithm_3d"]=get_option("Mesh.Algorithm3D")
        else
            result["tessella_version"]=string(Base.pkgversion(Tessella))
        end
        save(ARGS[5],result)
    finally
        implementation == "Gmsh" && finalize_gmsh()
    end
end
function latency_run()
    root=dirname(@__DIR__);output=abspath(ARGS[5])
    command=`$(Base.julia_cmd()[1]) --project=$root --startup-file=no --compiled-modules=yes --compile=yes --threads=1 -O2 $(@__FILE__) --cold $(ARGS[2]) $(ARGS[3]) $(ARGS[4]) $output`
    pipe=Pipe();started=time_ns();wall=nothing
    process=run(pipeline(command;stdout=pipe,stderr=stderr);wait=false)
    close(pipe.in)
    for line in eachline(pipe)
        if line == "MESH_RETURNED" && wall === nothing
            wall=(time_ns()-started)/1e9
        else
            println(line)
        end
    end
    wait(process)
    success(process) || error("Cold meshing process failed")
    wall === nothing && error("Cold process returned no mesh timing marker")
    result=TOML.parsefile(output)
    result["fresh_process_through_mesh_seconds"]=wall
    result["fresh_process_total_seconds"]=(time_ns()-started)/1e9
    save(output,result)
    println("Fresh process through returned mesh: ",wall," seconds")
end
function main()
    Threads.nthreads() == 1 || error("Run with JULIA_NUM_THREADS=1")
    BLAS.set_num_threads(1)
    LATENCY && return latency_run()
    !isempty(COLD) && return cold_run()
    initialize_gmsh()
    try
        meta = metadata(); rows = Dict{String,Any}[]; order_rng = MersenneTwister(meta["order_seed"])
        cases = benchmark_cases()
        # Configurations that disable Gmsh random factors are additional small
        # rows, never a replacement for its default-configuration outcomes.
        for zero_random in (false, true)
            zero_random && (set_option("Mesh.RandomFactor",0); set_option("Mesh.RandomFactor3D",0))
            selected = zero_random ? filter(c -> c[2] == 16, cases) : cases
            for (d,n,seed,family) in selected, constrained in (false,true)
                constrained && family != "simplex_enclosed" && continue
                algorithms = d == 3 ? (constrained ? (10,) : (1,10)) : (1,)
                for algorithm in algorithms
                points = cloud(d,n,seed,family)
                boundary = constrained ? simplex_boundary(d) : nothing
                # Gmsh's raw constrained 3D endpoint requires HXT. Retain the
                # default algorithm for every unconstrained case.
                default_algorithm = get_option("Mesh.Algorithm3D")
                set_option("Mesh.Algorithm3D",algorithm)
                record_geometry=n==16 && seed==9 && family=="uniform" && !constrained &&
                    !zero_random && (d==2 || algorithm==10)
                for row in paired_measure(points,meta["repetitions"],order_rng;
                    boundary=boundary,record_geometry=record_geometry)
                    merge!(row, Dict("dimension" => d, "points" => n, "input_seed" => seed,
                        "case" => family, "explicit_hull_constraints" => constrained,
                        "gmsh_algorithm_3d" => get_option("Mesh.Algorithm3D"),
                        "gmsh_zero_random_factors" => zero_random,
                        "input_sha256" => bytes2hex(sha256(reinterpret(UInt8,vec(points))))))
                    push!(rows,row)
                    println(row["implementation"], " d=",d," n=",n," ",family,
                        " seed=",seed," zero_random=",zero_random," ",row["status"],
                        " constrained=",constrained,
                        " median=",get(row,"warm_median_seconds",NaN),
                        " cells=",get(row,"cells",0)," ",row["validation_failures"]," ",row["delaunay_failures"])
                    flush(stdout)
                end
                set_option("Mesh.Algorithm3D",default_algorithm)
                # Keep partial results on interruption without declaring success.
                save(isempty(ARGS) ? "raw_pointcloud.toml" : ARGS[1], Dict("metadata"=>meta,"results"=>rows,"complete"=>false))
                end
            end
        end
        meta["source_changed_during_run"] = source_fingerprint()[2] != SOURCE_TREE_SHA256
        meta["driver_changed_during_run"] = bytes2hex(sha256(read(@__FILE__))) != DRIVER_SHA256
        save(isempty(ARGS) ? "raw_pointcloud.toml" : ARGS[1], Dict("metadata"=>meta,"results"=>rows,"complete"=>true))
    finally
        finalize_gmsh()
    end
end
function check_driver()
    @assert CASES isa Vector && length(CASES) == 56
    initialize_gmsh()
    try
        for d in (2,3)
            points = hcat(zeros(d),Matrix{Float64}(I,d,d))
            @assert validate_exact(points,gmsh_mesh(points)...)["status"] == "pass"
            previous = get_option("Mesh.Algorithm3D")
            d == 3 && set_option("Mesh.Algorithm3D",10)
            @assert validate_exact(points,gmsh_mesh(points,simplex_boundary(d))...)["status"] == "pass"
            set_option("Mesh.Algorithm3D",previous)
        end
        points = [0.0 1 1 0;0 0 1 2]
        coords = vcat(points,zeros(1,4))
        @assert validate_exact(points,coords,Int32[1 1;2 3;3 4])["status"] == "pass"
        rejected = false
        try
            validate_exact(points,coords,Int32[1 2;2 3;4 4])
        catch err
            rejected = occursin("strictly inside",sprint(showerror,err))
        end
        @assert rejected
        function must_reject(points,coords,cells,message)
            try
                validate_complete(points,coords,cells)
            catch err
                occursin(message,sprint(showerror,err)) || rethrow()
                return
            end
            error("Completion oracle accepted $message fixture")
        end
        points=[0.0 1 0 .25;0 0 1 .25];coords=vcat(points,zeros(1,4))
        good=Int32[1 2 3;2 3 1;4 4 4]
        @assert validate_complete(points,coords,good)["status"] == "pass"
        must_reject(points,coords,good[:,1:2],"supporting")
        must_reject(points,coords,hcat(good,good[:,1]),"Duplicate simplex")
        must_reject(points,coords,hcat(good,Int32[1,2,3]),"same side")
        flat_coords=hcat(coords,[.5,0.,0.])
        must_reject(flat_coords[1:2,:],flat_coords,hcat(good,Int32[1,2,5]),"Zero-measure")
        points=hcat(zeros(3),Matrix{Float64}(I,3,3))
        enlarged=hcat(points,[.5,.5,.5])
        must_reject(points,enlarged,Int32[1 2;2 3;3 4;4 5],"Steiner output")
        clipped=[0.0 1 1 .5 0;0 0 .5 1 1]
        square=vcat(hcat(clipped,[1.,1.]),zeros(1,6))
        must_reject(clipped,square,Int32[1 1 1 1;2 3 6 4;3 6 4 5],"Steiner output")
        points=hcat(zeros(2),Matrix{Float64}(I,2,2))
        rows=paired_measure(points,2,MersenneTwister(5);record_geometry=true,adapters=
            (("fixture_a",gmsh_mesh),("fixture_b",gmsh_mesh)))
        @assert length(rows)==2 && all(r->r["status"]=="pass",rows)
        io=IOBuffer();TOML.print(io,Dict("results"=>rows))
        println("Driver setup, public C ABI, exact completion/Delaunay and paired serialization checks passed")
    finally
        finalize_gmsh()
    end
end
const DRIVER_SHA256=bytes2hex(sha256(read(@__FILE__)))
load_started = time_ns()
if !CHECK_DRIVER && !LATENCY && COLD != "Gmsh"
    using Tessella
    realpath(pkgdir(Tessella)) == realpath(dirname(@__DIR__)) ||
        error("Activate this repository's project so source fingerprints identify the loaded Tessella")
end
const PACKAGE_LOAD_SECONDS = (time_ns() - load_started) / 1e9
if abspath(PROGRAM_FILE) == (@__FILE__)
    CHECK_DRIVER ? check_driver() : main()
end
