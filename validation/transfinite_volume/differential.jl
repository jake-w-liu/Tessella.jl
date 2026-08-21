#!/usr/bin/env julia
# In-memory differential for the bounded affine six-face transfinite-volume
# implementation against the installed Gmsh 4.15.2 public API. No geometry or
# mesh file is created.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.MeshTypes: nnodes, node, ntris, ntets, validate

if !isdefined(Tessella, :TransfiniteVolume)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "TransfiniteVolume.jl"))
end
using Tessella.TransfiniteVolume: mesh_transfinite_volume

const TARGET_GMSH_VERSION = "4.15.2"

function find_gmsh_executable()
    explicit = get(ENV, "GMSH_EXECUTABLE", "")
    if !isempty(explicit)
        isfile(explicit) || error("GMSH_EXECUTABLE does not name a file: $explicit")
        return realpath(explicit)
    end
    executable = Sys.which("gmsh")
    executable !== nothing && return realpath(executable)
    fallback = "/opt/homebrew/bin/gmsh"
    isfile(fallback) && return realpath(fallback)
    error("Gmsh 4.15.2 is required; install it or set GMSH_EXECUTABLE")
end

function find_gmsh_api(executable)
    explicit = get(ENV, "GMSH_JULIA_API", "")
    if !isempty(explicit)
        isfile(explicit) || error("GMSH_JULIA_API does not name a file: $explicit")
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

function affine_corners(origin, u, v, w)
    add(vectors...) = ntuple(d -> sum(vector[d] for vector in vectors), 3)
    return [origin, add(origin, u), add(origin, u, v), add(origin, v),
            add(origin, w), add(origin, u, w), add(origin, u, v, w),
            add(origin, v, w)]
end

function canonical_cells(connectivity, width::Int, mapping)
    length(connectivity) % width == 0 || error(
        "Gmsh connectivity length is not divisible by $width")
    if width == 3
        result = NTuple{3,Int32}[]
        for index in 1:3:length(connectivity)
            values = sort(Int32[mapping[connectivity[index]],
                                  mapping[connectivity[index + 1]],
                                  mapping[connectivity[index + 2]]])
            push!(result, (values[1], values[2], values[3]))
        end
    elseif width == 4
        result = NTuple{4,Int32}[]
        for index in 1:4:length(connectivity)
            values = sort(Int32[mapping[connectivity[index]],
                                  mapping[connectivity[index + 1]],
                                  mapping[connectivity[index + 2]],
                                  mapping[connectivity[index + 3]]])
            push!(result, (values[1], values[2], values[3], values[4]))
        end
    else
        error("unsupported canonical cell width $width")
    end
    sort!(result)
    return result
end

function canonical_tessella_cells(matrix)
    width = size(matrix, 1)
    if width == 3
        result = NTuple{3,Int32}[]
        for cell in axes(matrix, 2)
            values = sort(matrix[:, cell])
            push!(result, (values[1], values[2], values[3]))
        end
    elseif width == 4
        result = NTuple{4,Int32}[]
        for cell in axes(matrix, 2)
            values = sort(matrix[:, cell])
            push!(result, (values[1], values[2], values[3], values[4]))
        end
    else
        error("unsupported Tessella cell width $width")
    end
    sort!(result)
    return result
end

function add_affine_volume(corners, cells)
    gmsh.clear()
    gmsh.model.add("transfinite_volume")
    points = Int32[gmsh.model.geo.addPoint(point...) for point in corners]
    pairs = ((1, 2), (2, 3), (3, 4), (4, 1),
             (5, 6), (6, 7), (7, 8), (8, 5),
             (1, 5), (2, 6), (3, 7), (4, 8))
    edges = Int32[gmsh.model.geo.addLine(points[a], points[b]) for (a, b) in pairs]

    # Canonical Gmsh face order: vmin, umax, vmax, umin, wmin, wmax.
    signed_loops = ((edges[1], edges[10], -edges[5], -edges[9]),
                    (edges[2], edges[11], -edges[6], -edges[10]),
                    (-edges[3], edges[11], edges[7], -edges[12]),
                    (-edges[4], edges[12], edges[8], -edges[9]),
                    (edges[1], edges[2], edges[3], edges[4]),
                    (edges[5], edges[6], edges[7], edges[8]))
    face_corners = ((1, 2, 6, 5), (2, 3, 7, 6), (4, 3, 7, 8),
                    (1, 4, 8, 5), (1, 2, 3, 4), (5, 6, 7, 8))
    faces = Int32[]
    for loop_edges in signed_loops
        loop = gmsh.model.geo.addCurveLoop(Int32[loop_edges...])
        push!(faces, gmsh.model.geo.addPlaneSurface(Int32[loop]))
    end
    shell = gmsh.model.geo.addSurfaceLoop(
        Int32[faces[1], faces[2], -faces[3], -faces[4], -faces[5], faces[6]])
    volume = gmsh.model.geo.addVolume(Int32[shell])

    nu, nv, nw = cells
    for index in (1, 3, 5, 7)
        gmsh.model.geo.mesh.setTransfiniteCurve(edges[index], nu + 1)
    end
    for index in (2, 4, 6, 8)
        gmsh.model.geo.mesh.setTransfiniteCurve(edges[index], nv + 1)
    end
    for index in (9, 10, 11, 12)
        gmsh.model.geo.mesh.setTransfiniteCurve(edges[index], nw + 1)
    end
    for index in 1:6
        gmsh.model.geo.mesh.setTransfiniteSurface(
            faces[index], "Left", Int32[points[i] for i in face_corners[index]])
    end
    gmsh.model.geo.mesh.setTransfiniteVolume(volume, points)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    return faces, volume
end

@inline distance(a, b) = hypot(a[1] - b[1], a[2] - b[2], a[3] - b[3])

function node_mapping(mesh)
    tags, coordinates, _ = gmsh.model.mesh.getNodes()
    length(coordinates) == 3length(tags) || error("Gmsh returned malformed node coordinates")
    length(tags) == nnodes(mesh) || error(
        "node-count mismatch: Gmsh $(length(tags)), Tessella $(nnodes(mesh))")
    scale = maximum(abs, coordinates; init=1.0)
    # Gmsh obtains straight-curve nodes through its parametric inversion path;
    # on the pinned build this leaves O(1e-12) coordinate residuals even for
    # affine lines. Keep the differential tolerance below 5e-11 at unit scale.
    tolerance = 65_536eps(Float64) * max(scale, 1.0)
    used = falses(nnodes(mesh))
    mapping = Dict{UInt64,Int32}()
    maximum_error = 0.0
    for source in eachindex(tags)
        point = (coordinates[3source - 2], coordinates[3source - 1], coordinates[3source])
        best = 0
        best_error = Inf
        for destination in 1:nnodes(mesh)
            used[destination] && continue
            candidate_error = distance(point, node(mesh, destination))
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

function check_case(corners, cells)
    faces, volume = add_affine_volume(corners, cells)
    mesh = mesh_transfinite_volume(corners, cells)
    validate(mesh).ok || error("Tessella transfinite volume did not validate")
    nu, nv, nw = cells
    expected_nodes = (nu + 1) * (nv + 1) * (nw + 1)
    expected_tets = 6nu * nv * nw
    expected_triangles = 4(nu * nv + nu * nw + nv * nw)
    (nnodes(mesh), ntets(mesh), ntris(mesh)) ==
        (expected_nodes, expected_tets, expected_triangles) ||
        error("unexpected Tessella output counts")
    mapping, maximum_error = node_mapping(mesh)

    types, _, element_nodes = gmsh.model.mesh.getElements(3, volume)
    tetrahedron_position = findfirst(==(Int32(4)), types)
    tetrahedron_position === nothing && error("Gmsh emitted no type-4 tetrahedra")
    length(types) == 1 || error("Gmsh emitted unexpected 3-D element types: $types")
    canonical_cells(element_nodes[tetrahedron_position], 4, mapping) ==
        canonical_tessella_cells(mesh.tets) ||
        error("tetrahedron connectivity differs from Gmsh")

    gmsh_triangles = UInt64[]
    for face in faces
        face_types, _, face_nodes = gmsh.model.mesh.getElements(2, face)
        triangle_position = findfirst(==(Int32(2)), face_types)
        triangle_position === nothing && error("Gmsh face $face emitted no type-2 triangles")
        length(face_types) == 1 || error(
            "Gmsh face $face emitted unexpected 2-D element types: $face_types")
        append!(gmsh_triangles, face_nodes[triangle_position])
    end
    canonical_cells(gmsh_triangles, 3, mapping) ==
        canonical_tessella_cells(mesh.tris) ||
        error("boundary triangle arrangements differ from Gmsh")
    return maximum_error, expected_nodes, expected_tets, expected_triangles
end

gmsh.initialize([GMSH_EXECUTABLE, "-nopopup"], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Smoothing", 0)
    runtime_version = gmsh.option.getString("General.Version")
    (runtime_version == TARGET_GMSH_VERSION ||
     startswith(runtime_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh runtime $TARGET_GMSH_VERSION, got $runtime_version")
    cases = ((affine_corners((0., 0., 0.), (3., 0., 0.),
                             (0., 2., 0.), (0., 0., 1.5)), (3, 2, 2)),
             (affine_corners((1.25, -2., 0.5), (2., 0.4, -0.2),
                             (-0.35, 1.6, 0.25), (0.15, -0.3, 1.4)), (2, 3, 2)))
    maximum_error = 0.0
    total_nodes = 0
    total_tets = 0
    total_triangles = 0
    for (corners, cells) in cases
        error_value, nodes, tets, triangles = check_case(corners, cells)
        maximum_error = max(maximum_error, error_value)
        total_nodes += nodes
        total_tets += tets
        total_triangles += triangles
    end
    println("TRANSFINITE_VOLUME_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$(length(cases)) nodes=$total_nodes tets=$total_tets " *
            "boundary_triangles=$total_triangles max_node_error=$maximum_error")
finally
    gmsh.finalize()
end
