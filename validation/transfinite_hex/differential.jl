#!/usr/bin/env julia
# In-memory differential for fully recombined affine six-face transfinite
# volumes against pinned Gmsh 4.15.2. No geometry or mesh file is created.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella

if !isdefined(Tessella, :TransfiniteHex)
    Base.include(Tessella,
                 joinpath(@__DIR__, "..", "..", "src", "structured",
                          "TransfiniteHex.jl"))
end
using Tessella.TransfiniteHex: mesh_transfinite_hex

const TARGET_GMSH_VERSION = "4.15.2"
const FACE_TAGS = (11, 12, 13, 14, 15, 16)
const VOLUME_TAG = 21
const ARRANGEMENTS = Dict(
    :left => "Left",
    :right => "Right",
    :alternate_left => "AlternateLeft",
    :alternate_right => "AlternateRight")

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
    error("Gmsh $TARGET_GMSH_VERSION is required; install it or set " *
          "GMSH_EXECUTABLE")
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

function affine_corners(origin, u, v, w)
    add(vectors...) = ntuple(d -> sum(vector[d] for vector in vectors), 3)
    return [origin, add(origin, u), add(origin, u, v), add(origin, v),
            add(origin, w), add(origin, u, w), add(origin, u, v, w),
            add(origin, v, w)]
end

function add_affine_volume(corners, cells, arrangement::String)
    gmsh.clear()
    gmsh.model.add("transfinite_hex_" * lowercase(arrangement))
    points = Int32[gmsh.model.geo.addPoint(point...) for point in corners]
    pairs = ((1, 2), (2, 3), (3, 4), (4, 1),
             (5, 6), (6, 7), (7, 8), (8, 5),
             (1, 5), (2, 6), (3, 7), (4, 8))
    edges = Int32[
        gmsh.model.geo.addLine(points[first], points[second])
        for (first, second) in pairs]

    # Canonical order: vmin, umax, vmax, umin, wmin, wmax.
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
        Int32[faces[1], faces[2], -faces[3], -faces[4],
              -faces[5], faces[6]])
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
            faces[index], arrangement,
            Int32[points[corner] for corner in face_corners[index]])
        gmsh.model.geo.mesh.setRecombine(2, faces[index])
    end
    gmsh.model.geo.mesh.setTransfiniteVolume(volume, points)
    gmsh.model.geo.synchronize()
    for index in 1:6
        physical = gmsh.model.addPhysicalGroup(
            2, Int32[faces[index]], FACE_TAGS[index])
        physical == FACE_TAGS[index] || error(
            "Gmsh changed requested face physical tag $(FACE_TAGS[index])")
    end
    physical_volume = gmsh.model.addPhysicalGroup(3, Int32[volume], VOLUME_TAG)
    physical_volume == VOLUME_TAG || error(
        "Gmsh changed requested volume physical tag $VOLUME_TAG")
    gmsh.model.mesh.generate(3)
    return faces, volume
end

@inline distance(a, b) =
    hypot(a[1] - b[1], a[2] - b[2], a[3] - b[3])

function node_mapping(mesh)
    tags, coordinates, _ = gmsh.model.mesh.getNodes()
    length(coordinates) == 3length(tags) || error(
        "Gmsh returned malformed node coordinates")
    length(tags) == size(mesh.coords, 2) || error(
        "node-count mismatch: Gmsh $(length(tags)), " *
        "Tessella $(size(mesh.coords, 2))")
    scale = maximum(abs, coordinates; init=1.0)
    tolerance = 65_536eps(Float64) * max(scale, 1.0)
    used = falses(size(mesh.coords, 2))
    mapping = Dict{UInt64,Int32}()
    maximum_error = 0.0
    for source in eachindex(tags)
        point = (coordinates[3source - 2], coordinates[3source - 1],
                 coordinates[3source])
        best = 0
        best_error = Inf
        for destination in axes(mesh.coords, 2)
            used[destination] && continue
            candidate = (mesh.coords[1, destination],
                         mesh.coords[2, destination],
                         mesh.coords[3, destination])
            candidate_error = distance(point, candidate)
            if candidate_error < best_error
                best = destination
                best_error = candidate_error
            end
        end
        best != 0 && best_error <= tolerance || error(
            "no Tessella node matches Gmsh node $(tags[source]); nearest " *
            "error=$best_error, tolerance=$tolerance")
        mapping[tags[source]] = Int32(best)
        used[best] = true
        maximum_error = max(maximum_error, best_error)
    end
    all(used) || error("some Tessella nodes were not matched to Gmsh nodes")
    return mapping, maximum_error
end

function check_case(corners, cells, arrangement_symbol, arrangement_string)
    faces, volume = add_affine_volume(corners, cells, arrangement_string)
    mesh = mesh_transfinite_hex(
        corners, cells; arrangement=arrangement_symbol,
        face_tags=FACE_TAGS, volume_tag=VOLUME_TAG)
    Tessella.Elements.validate(mesh).ok || error(
        "Tessella $arrangement_string transfinite hex did not validate")
    length(mesh.blocks) == 2 || error("Tessella emitted an unexpected block count")
    quadrangles, hexahedra = mesh.blocks
    (quadrangles.msh, hexahedra.msh) == (3, 5) || error(
        "Tessella emitted unexpected MSH types")

    nu, nv, nw = cells
    expected_nodes = (nu + 1) * (nv + 1) * (nw + 1)
    expected_hexahedra = nu * nv * nw
    expected_quadrangles = 2 * (nu * nv + nu * nw + nv * nw)
    (size(mesh.coords, 2), size(hexahedra.nodes, 2),
     size(quadrangles.nodes, 2)) ==
        (expected_nodes, expected_hexahedra, expected_quadrangles) ||
        error("unexpected Tessella output counts")
    hexahedra.tags == fill(Int32(VOLUME_TAG), expected_hexahedra) || error(
        "Tessella volume tags differ from the requested physical tag")
    expected_face_counts = (nu * nw, nv * nw, nu * nw,
                            nv * nw, nu * nv, nu * nv)
    offset = 0
    for face in 1:6
        span = offset + 1:offset + expected_face_counts[face]
        all(==(Int32(FACE_TAGS[face])), quadrangles.tags[span]) || error(
            "Tessella face $face tags differ from the requested physical tag")
        gmsh.model.getPhysicalGroupsForEntity(2, faces[face]) ==
            Int32[FACE_TAGS[face]] || error(
                "Gmsh face $face physical-tag correspondence changed")
        offset += expected_face_counts[face]
    end
    gmsh.model.getPhysicalGroupsForEntity(3, volume) == Int32[VOLUME_TAG] ||
        error("Gmsh volume physical-tag correspondence changed")

    mapping, maximum_error = node_mapping(mesh)
    types, _, element_nodes = gmsh.model.mesh.getElements(3, volume)
    types == Int32[5] || error(
        "Gmsh emitted volume types $types instead of only type 5")
    mapped_hexahedra = Int32[
        mapping[tag] for tag in element_nodes[1]]
    mapped_hexahedra == vec(hexahedra.nodes) || error(
        "$arrangement_string ordered hexahedron connectivity differs from Gmsh")

    mapped_quadrangles = Int32[]
    for face in faces
        face_types, _, face_nodes = gmsh.model.mesh.getElements(2, face)
        face_types == Int32[3] || error(
            "Gmsh face $face emitted types $face_types instead of only type 3")
        append!(mapped_quadrangles,
                (mapping[tag] for tag in face_nodes[1]))
    end
    mapped_quadrangles == vec(quadrangles.nodes) || error(
        "$arrangement_string ordered boundary connectivity differs from Gmsh")
    return maximum_error, expected_nodes, expected_hexahedra,
           expected_quadrangles
end

gmsh.initialize([GMSH_EXECUTABLE, "-nopopup"], false, false)
try
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Smoothing", 0)
    runtime_version = gmsh.option.getString("General.Version")
    (runtime_version == TARGET_GMSH_VERSION ||
     startswith(runtime_version, TARGET_GMSH_VERSION * "-")) || error(
        "expected Gmsh runtime $TARGET_GMSH_VERSION, got $runtime_version")
    cases = ((affine_corners(
                  (0.0, 0.0, 0.0), (3.0, 0.0, 0.0),
                  (0.0, 2.0, 0.0), (0.0, 0.0, 1.5)), (3, 2, 2)),
             (affine_corners(
                  (1.25, -2.0, 0.5), (2.0, 0.4, -0.2),
                  (-0.35, 1.6, 0.25), (0.15, -0.3, 1.4)), (2, 3, 2)))
    maximum_error = 0.0
    total_nodes = 0
    total_hexahedra = 0
    total_quadrangles = 0
    checked = 0
    for (symbol, name) in sort!(collect(ARRANGEMENTS); by=first)
        for (corners, cells) in cases
            error_value, nodes, hexahedra, quadrangles =
                check_case(corners, cells, symbol, name)
            maximum_error = max(maximum_error, error_value)
            total_nodes += nodes
            total_hexahedra += hexahedra
            total_quadrangles += quadrangles
            checked += 1
        end
    end
    println("TRANSFINITE_HEX_DIFFERENTIAL_OK gmsh=$runtime_version " *
            "cases=$checked nodes=$total_nodes hexahedra=$total_hexahedra " *
            "boundary_quadrangles=$total_quadrangles " *
            "max_node_error=$maximum_error")
finally
    gmsh.finalize()
end
