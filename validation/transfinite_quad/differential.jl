#!/usr/bin/env julia
# In-memory differential for recombined four-sided transfinite quadrangle
# patches against pinned Gmsh 4.15.2. No geometry or mesh file is created.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Elements: mixed_crc

if !isdefined(Tessella, :TransfiniteQuad)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "TransfiniteQuad.jl"))
end
using Tessella.TransfiniteQuad: mesh_transfinite_quad_patch

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
gmsh.GMSH_API_VERSION == TARGET_GMSH_VERSION || error(
    "expected Gmsh API $TARGET_GMSH_VERSION, got $(gmsh.GMSH_API_VERSION)")

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

function canonical_segments(connectivity)
    length(connectivity) % 2 == 0 || error(
        "line connectivity is not divisible by two")
    result = NTuple{2,Int32}[]
    for index in 1:2:length(connectivity)
        first = Int32(connectivity[index])
        second = Int32(connectivity[index + 1])
        push!(result, first < second ? (first, second) : (second, first))
    end
    sort!(result)
    return result
end

function canonical_quadrangles(connectivity)
    length(connectivity) % 4 == 0 || error(
        "quadrangle connectivity is not divisible by four")
    result = NTuple{4,Int32}[]
    for index in 1:4:length(connectivity)
        values = sort(Int32[connectivity[index], connectivity[index + 1],
                            connectivity[index + 2], connectivity[index + 3]])
        push!(result, (values[1], values[2], values[3], values[4]))
    end
    sort!(result)
    return result
end

function add_curved_patch(arrangement, counts; tilted=false)
    gmsh.clear()
    gmsh.model.add("transfinite_quad_" * arrangement *
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
    pb1 = add_point(1.0, -0.4)
    pb2 = add_point(2.5, -0.7)
    pb3 = add_point(3.2, -0.2)
    p2 = add_point(4.0, 0.0)
    pr1 = add_point(4.4, 1.0)
    pr2 = add_point(4.2, 2.0)
    p3 = add_point(4.0, 3.0)
    pt1 = add_point(3.0, 3.6)
    pt2 = add_point(1.3, 3.2)
    p4 = add_point(0.0, 3.0)
    pl1 = add_point(-0.3, 2.1)
    pl2 = add_point(-0.6, 0.7)
    curves = (gmsh.model.geo.addSpline([p1, pb1, pb2, pb3, p2]),
              gmsh.model.geo.addSpline([p2, pr1, pr2, p3]),
              gmsh.model.geo.addSpline([p3, pt1, pt2, p4]),
              gmsh.model.geo.addSpline([p4, pl1, pl2, p1]))
    loop = gmsh.model.geo.addCurveLoop(collect(curves))
    surface = gmsh.model.geo.addPlaneSurface([loop])
    for (curve, count) in zip(curves, counts)
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, count)
    end
    gmsh.model.geo.mesh.setTransfiniteSurface(
        surface, arrangement, [p1, p2, p3, p4])
    gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    return curves, surface
end

function gmsh_to_tessella_node_map(mesh, surface)
    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    quadrangle_position = findfirst(==(Int32(3)), types)
    quadrangle_position === nothing && error(
        "Gmsh emitted no first-order quadrangles")
    tags = sort!(unique(element_nodes[quadrangle_position]))
    length(tags) == size(mesh.coords, 2) || error(
        "node-count mismatch: Gmsh $(length(tags)), Tessella $(size(mesh.coords, 2))")
    used = falses(size(mesh.coords, 2))
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
        for destination in axes(mesh.coords, 2)
            used[destination] && continue
            candidate = (mesh.coords[1, destination], mesh.coords[2, destination],
                         mesh.coords[3, destination])
            candidate_error = _distance(point, candidate)
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

function check_case(arrangement, symbol, counts=(5, 4, 5, 4); tilted=false)
    curves, surface = add_curved_patch(arrangement, counts; tilted=tilted)
    sides = map(curve_points, curves)
    all(length(sides[index]) == counts[index] for index in 1:4) || error(
        "Gmsh boundary node counts differ from requested counts $counts")
    mesh = mesh_transfinite_quad_patch(
        sides...; arrangement=symbol, face_tag=21,
        side_tags=(11, 12, 13, 14))
    Tessella.Elements.validate(mesh).ok || error(
        "Tessella $arrangement quadrangle patch did not validate")
    horizontal = counts[1] - 1
    vertical = counts[2] - 1
    expected = (counts[1] * counts[2], 2 * (horizontal + vertical),
                horizontal * vertical)
    (size(mesh.coords, 2), size(mesh.blocks[1].nodes, 2),
     size(mesh.blocks[2].nodes, 2)) == expected || error(
        "unexpected Tessella $arrangement counts")
    mapping, maximum_error = gmsh_to_tessella_node_map(mesh, surface)

    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    types == Int32[3] || error(
        "Gmsh $arrangement emitted element types $types instead of only type 3")
    mapped_quadrangles = Int32[mapping[tag] for tag in element_nodes[1]]
    mapped_quadrangles == vec(mesh.blocks[2].nodes) || error(
        "$arrangement ordered quadrangle connectivity differs from Gmsh")

    mapped_segments = Int32[]
    for curve in curves
        line_types, _, line_nodes = gmsh.model.mesh.getElements(1, curve)
        line_position = findfirst(==(Int32(1)), line_types)
        line_position === nothing && error(
            "Gmsh curve $curve emitted no first-order lines")
        append!(mapped_segments,
                (mapping[tag] for tag in line_nodes[line_position]))
    end
    mapped_segments == vec(mesh.blocks[1].nodes) || error(
        "$arrangement ordered boundary connectivity differs from Gmsh")
    return maximum_error, canonical_quadrangles(mapped_quadrangles),
           mixed_crc(mesh).sha, expected[1]
end

gmsh.initialize([GMSH_EXECUTABLE, "-nopopup"], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Smoothing", 1)
    gmsh.option.setNumber("Mesh.QuasiTransfinite", 0)
    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    api_version = gmsh.option.getString("General.Version")
    (api_version == TARGET_GMSH_VERSION ||
     startswith(api_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh API $TARGET_GMSH_VERSION, got $api_version")

    cases = (("Left", :left), ("Right", :right),
             ("AlternateLeft", :alternate_left),
             ("AlternateRight", :alternate_right))
    errors = Float64[]
    topologies = Vector{NTuple{4,Int32}}[]
    crcs = String[]
    coordinate_samples = 0
    for (arrangement, symbol) in cases
        maximum_error, topology, crc, nodes = check_case(arrangement, symbol)
        push!(errors, maximum_error)
        push!(topologies, topology)
        push!(crcs, crc)
        coordinate_samples += nodes
    end
    all(topology == topologies[1] for topology in topologies) || error(
        "recombined topology changed with triangle arrangement")
    length(unique(crcs)) == 1 || error(
        "Tessella recombined CRC changed with triangle arrangement")

    for counts in ((2, 2, 2, 2), (3, 5, 3, 5), (8, 4, 8, 4))
        maximum_error, _, _, nodes = check_case("Left", :left, counts)
        push!(errors, maximum_error)
        coordinate_samples += nodes
    end
    maximum_error, _, _, nodes = check_case(
        "Left", :left, (5, 4, 5, 4); tilted=true)
    push!(errors, maximum_error)
    coordinate_samples += nodes

    println("TRANSFINITE_QUAD_DIFFERENTIAL_OK gmsh=$api_version " *
            "arrangements=$(length(cases)) resolutions=4 geometries=2 " *
            "coordinate_samples=$coordinate_samples " *
            "max_node_error=$(maximum(errors)) reference_nodes=20 " *
            "reference_segments=14 reference_quadrangles=12 " *
            "arrangement_topologies=1 ordered_connectivity=1")
finally
    gmsh.finalize()
end
