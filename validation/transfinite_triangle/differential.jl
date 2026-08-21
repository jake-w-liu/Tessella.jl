#!/usr/bin/env julia
# In-memory differential for the specific Gmsh 4.15.2 three-sided structured
# transfinite-triangle algorithm (`Mesh.TransfiniteTri = 1`). No geometry or
# mesh file is created.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: nnodes, node, nsegs, ntris, validate

if !isdefined(Tessella, :TransfiniteTriangle)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src",
                          "TransfiniteTriangle.jl"))
end
using Tessella.TransfiniteTriangle: mesh_transfinite_triangle

const TARGET_GMSH_VERSION = "4.15.2"

function find_gmsh_executable()
    explicit = get(ENV, "GMSH_EXECUTABLE", "")
    if !isempty(explicit)
        isfile(explicit) || error(
            "GMSH_EXECUTABLE does not name a file: $explicit")
        return realpath(explicit)
    end
    executable = Sys.which("gmsh")
    executable !== nothing && return realpath(executable)
    fallback = "/opt/homebrew/bin/gmsh"
    isfile(fallback) && return realpath(fallback)
    error("Gmsh $TARGET_GMSH_VERSION is required; install it or set GMSH_EXECUTABLE")
end

function find_gmsh_api(executable)
    explicit = get(ENV, "GMSH_JULIA_API", "")
    if !isempty(explicit)
        isfile(explicit) || error(
            "GMSH_JULIA_API does not name a file: $explicit")
        return realpath(explicit)
    end
    prefix = dirname(dirname(executable))
    candidates = (joinpath(prefix, "lib", "gmsh.jl"),
                  joinpath(prefix, "lib64", "gmsh.jl"),
                  "/opt/homebrew/lib/gmsh.jl",
                  "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                  "/usr/local/opt/gmsh/lib/gmsh.jl")
    for candidate in candidates
        isfile(candidate) && return realpath(candidate)
    end
    error("could not locate gmsh.jl for $executable; set GMSH_JULIA_API")
end

const GMSH_EXECUTABLE = find_gmsh_executable()
const GMSH_CLI_VERSION = strip(read(`$GMSH_EXECUTABLE --version`, String))
(GMSH_CLI_VERSION == TARGET_GMSH_VERSION ||
 startswith(GMSH_CLI_VERSION, TARGET_GMSH_VERSION * "-")) ||
    error("expected Gmsh $TARGET_GMSH_VERSION, got $GMSH_CLI_VERSION")
const GMSH_API_FILE = find_gmsh_api(GMSH_EXECUTABLE)
include(GMSH_API_FILE)

@inline _distance(a, b) =
    hypot(a[1] - b[1], a[2] - b[2], a[3] - b[3])

function curve_points(curve)
    tags, coordinates, parameters =
        gmsh.model.mesh.getNodes(1, curve, true, true)
    length(tags) == length(parameters) || error(
        "Gmsh curve $curve did not return one parameter per node")
    order = sortperm(parameters)
    return [(coordinates[3index - 2], coordinates[3index - 1],
             coordinates[3index]) for index in order]
end

function canonical_triangles(connectivity)
    length(connectivity) % 3 == 0 || error(
        "triangle connectivity is not divisible by three")
    result = NTuple{3,Int32}[]
    for index in 1:3:length(connectivity)
        values = sort(Int32[connectivity[index], connectivity[index + 1],
                            connectivity[index + 2]])
        push!(result, (values[1], values[2], values[3]))
    end
    sort!(result)
    return result
end

function canonical_segments(connectivity)
    length(connectivity) % 2 == 0 || error(
        "segment connectivity is not divisible by two")
    result = NTuple{2,Int32}[]
    for index in 1:2:length(connectivity)
        first = Int32(connectivity[index])
        second = Int32(connectivity[index + 1])
        push!(result, first < second ? (first, second) : (second, first))
    end
    sort!(result)
    return result
end

function add_curved_triangle(arrangement, count; tilted=false)
    gmsh.clear()
    gmsh.model.add("transfinite_triangle_" * arrangement *
                   (tilted ? "_tilted" : ""))
    ex = (inv(sqrt(2.0)), inv(sqrt(2.0)), 0.0)
    ey = (-inv(sqrt(6.0)), inv(sqrt(6.0)), 2inv(sqrt(6.0)))
    origin = (1.0, -2.0, 3.0)
    coordinates(x, y) = tilted ?
        (origin[1] + x * ex[1] + y * ey[1],
         origin[2] + x * ex[2] + y * ey[2],
         origin[3] + x * ex[3] + y * ey[3]) : (x, y, 0.0)
    add_point(x, y) = gmsh.model.geo.addPoint(coordinates(x, y)...)
    p1 = add_point(0.0, 0.0)
    p12a = add_point(0.8, -0.35)
    p12b = add_point(2.7, -0.55)
    p2 = add_point(4.0, 0.0)
    p23a = add_point(4.2, 0.85)
    p23b = add_point(2.3, 2.5)
    p3 = add_point(0.4, 3.2)
    p31a = add_point(-0.2, 2.35)
    p31b = add_point(-0.35, 0.75)
    curves = (gmsh.model.geo.addSpline([p1, p12a, p12b, p2]),
              gmsh.model.geo.addSpline([p2, p23a, p23b, p3]),
              gmsh.model.geo.addSpline([p3, p31a, p31b, p1]))
    loop = gmsh.model.geo.addCurveLoop(collect(curves))
    surface = gmsh.model.geo.addPlaneSurface([loop])
    for curve in curves
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, count)
    end
    gmsh.model.geo.mesh.setTransfiniteSurface(
        surface, arrangement, [p1, p2, p3])
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    return curves, surface
end

function gmsh_to_tessella_node_map(mesh, surface)
    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    triangle_position = findfirst(==(Int32(2)), types)
    triangle_position === nothing && error(
        "Gmsh emitted no first-order triangles")
    tags = sort!(unique(element_nodes[triangle_position]))
    length(tags) == nnodes(mesh) || error(
        "node-count mismatch: Gmsh $(length(tags)), Tessella $(nnodes(mesh))")
    used = falses(nnodes(mesh))
    mapping = Dict{UInt64,Int32}()
    maximum_error = 0.0
    coordinates = Dict(tag => gmsh.model.mesh.getNode(tag)[1] for tag in tags)
    scale = maximum((maximum(abs, point; init=0.0)
                     for point in values(coordinates)); init=1.0)
    tolerance = 512eps(Float64) * max(scale, 1.0)
    for source in eachindex(tags)
        raw = coordinates[tags[source]]
        point = (raw[1], raw[2], raw[3])
        best = 0
        best_error = Inf
        for destination in 1:nnodes(mesh)
            used[destination] && continue
            candidate_error = _distance(point, node(mesh, destination))
            if candidate_error < best_error
                best = destination
                best_error = candidate_error
            end
        end
        best != 0 && best_error <= tolerance || error(
            "no Tessella node matches Gmsh node $(tags[source]); " *
            "nearest error=$best_error, tolerance=$tolerance")
        mapping[tags[source]] = Int32(best)
        used[best] = true
        maximum_error = max(maximum_error, best_error)
    end
    all(used) || error("some Tessella nodes were not matched to Gmsh nodes")
    return mapping, maximum_error
end

function check_arrangement(arrangement, symbol, count=6; tilted=false)
    curves, surface = add_curved_triangle(arrangement, count; tilted=tilted)
    sides = map(curve_points, curves)
    all(length(side) == count for side in sides) || error(
        "Gmsh did not emit $count nodes on each transfinite curve")
    mesh = mesh_transfinite_triangle(
        sides...; arrangement=symbol, face_tag=21, side_tags=(11, 12, 13))
    validate(mesh).ok || error("Tessella $arrangement patch did not validate")
    divisions = count - 1
    expected = (count * (count + 1) ÷ 2, 3divisions, divisions^2)
    (nnodes(mesh), nsegs(mesh), ntris(mesh)) == expected || error(
        "unexpected Tessella $arrangement counts")
    mapping, maximum_error = gmsh_to_tessella_node_map(mesh, surface)

    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    types == Int32[2] || error(
        "Gmsh $arrangement emitted element types $types instead of only type 2")
    mapped_triangles = Int32[
        mapping[tag] for tag in element_nodes[1]]
    canonical_triangles(mapped_triangles) ==
        canonical_triangles(vec(mesh.tris)) || error(
            "$arrangement triangle connectivity differs from Gmsh")

    mapped_segments = Int32[]
    for curve in curves
        line_types, _, line_nodes = gmsh.model.mesh.getElements(1, curve)
        line_position = findfirst(==(Int32(1)), line_types)
        line_position === nothing && error(
            "Gmsh curve $curve emitted no first-order lines")
        append!(mapped_segments,
                (mapping[tag] for tag in line_nodes[line_position]))
    end
    canonical_segments(mapped_segments) ==
        canonical_segments(vec(mesh.segs)) || error(
            "$arrangement boundary connectivity differs from Gmsh")
    return maximum_error, canonical_triangles(mapped_triangles)
end

gmsh.initialize([GMSH_EXECUTABLE, "-nopopup"], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Smoothing", 1)
    gmsh.option.setNumber("Mesh.QuasiTransfinite", 0)
    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    gmsh.option.setNumber("Mesh.TransfiniteTri", 1)
    api_version = gmsh.option.getString("General.Version")
    (api_version == TARGET_GMSH_VERSION ||
     startswith(api_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh API $TARGET_GMSH_VERSION, got $api_version")

    cases = (("Left", :left), ("Right", :right),
             ("AlternateLeft", :alternate_left),
             ("AlternateRight", :alternate_right))
    errors = Float64[]
    topologies = Vector{NTuple{3,Int32}}[]
    coordinate_samples = 0
    for (arrangement, symbol) in cases
        maximum_error, topology = check_arrangement(arrangement, symbol)
        push!(errors, maximum_error)
        push!(topologies, topology)
        coordinate_samples += 21
    end
    all(topology == topologies[1] for topology in topologies) || error(
        "unrecombined specific triangular topology changed with arrangement")

    # Exercise the minimal patch and two additional triangular-lattice sizes;
    # arrangement parity above uses the nonuniform six-node case.
    for count in (2, 3, 8)
        maximum_error, _ = check_arrangement("Left", :left, count)
        push!(errors, maximum_error)
        coordinate_samples += count * (count + 1) ÷ 2
    end
    maximum_error, _ = check_arrangement("Left", :left, 6; tilted=true)
    push!(errors, maximum_error)
    coordinate_samples += 21
    println("TRANSFINITE_TRIANGLE_DIFFERENTIAL_OK gmsh=$api_version " *
            "arrangements=$(length(cases)) resolutions=4 geometries=2 " *
            "coordinate_samples=$coordinate_samples " *
            "max_node_error=$(maximum(errors)) reference_nodes=21 " *
            "reference_segments=15 reference_triangles=25 " *
            "arrangement_topologies=1")
finally
    gmsh.finalize()
end
