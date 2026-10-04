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
                 joinpath(@__DIR__, "..", "..", "src", "structured",
                          "TransfiniteTriangle.jl"))
end
using Tessella.TransfiniteTriangle: mesh_transfinite_triangle,
                                    mesh_transfinite_triangle_patch,
                                    mesh_transfinite_triangle_collapsed,
                                    mesh_transfinite_triangle_collapsed_patch

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

function canonical_quadrangles(connectivity)
    length(connectivity) % 4 == 0 || error(
        "quadrangle connectivity is not divisible by four")
    result = NTuple{4,Int32}[]
    for index in 1:4:length(connectivity)
        nodes = ntuple(offset -> Int32(connectivity[index + offset - 1]), 4)
        candidates = NTuple{4,Int32}[]
        for shift in 0:3
            push!(candidates,
                  ntuple(offset -> nodes[mod1(shift + offset, 4)], 4))
        end
        reversed = (nodes[1], nodes[4], nodes[3], nodes[2])
        for shift in 0:3
            push!(candidates,
                  ntuple(offset -> reversed[mod1(shift + offset, 4)], 4))
        end
        push!(result, minimum(candidates))
    end
    sort!(result)
    return result
end

function add_curved_triangle(arrangement, count;
                             tilted=false, recombine=false,
                             counts=nothing, corners=nothing)
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
    per_curve = counts === nothing ? (count, count, count) : counts
    for (curve, nodes) in zip(curves, per_curve)
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, nodes)
    end
    gmsh.model.geo.mesh.setTransfiniteSurface(
        surface, arrangement,
        corners === nothing ? [p1, p2, p3] : corners)
    recombine && gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    return curves, surface
end

# Public CAD-relative winding is a separate contract from the canonical side
# frame used by the low-level interpolation and transfinite-volume kernels.
# Gmsh Generator calls orientMeshGFace after filling, then applies reverseMesh.
function public_orientation_cases()
    result=NamedTuple[]
    for algorithm in (0,1),direction in (-1,1),count in (2,5)
        push!(result,(;algorithm,direction,count))
    end
    for direction in (-1,1)
        for corners in ((1,2,3),(1,3,2))
            push!(result,(;direction,corners))
        end
        push!(result,(;direction,reverse=true))
        push!(result,(;direction,recombine=true,count=5))
        push!(result,(;direction,axis=:xz))
    end
    return result
end

function exact_public_cell_orientation(coordinates,nodes,normal)
    a,b,c=(Rational{BigInt}.(coordinates[:,node]) for node in nodes[1:3])
    u=b.-a;v=c.-a
    cross=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    return sign(sum(cross[d]*Rational{BigInt}(normal[d]) for d in 1:3))
end

function check_public_orientation(case)
    native=Tessella.GeoModel()
    algorithm=get(case,:algorithm,1);count=get(case,:count,2)
    corners=collect(Int,get(case,:corners,()))
    points=get(case,:axis,:xy)===:xy ?
        ((0.,0.,0.),(1.,0.,0.),(0.,1.,0.)) :
        ((0.,0.,0.),(1.,0.,0.),(0.,0.,1.))
    loop=case.direction==1 ? [1,2,3] : [-3,-2,-1]
    gmsh.clear();gmsh.model.add("public_triangle_orientation")
    gmsh.option.setNumber("Mesh.TransfiniteTri",algorithm)
    Tessella.Model.set_transfinite_tri!(native,algorithm)
    for (tag,point) in enumerate(points)
        gmsh.model.geo.addPoint(point...,1.,tag)
        Tessella.Model.add_point!(native,point...;tag,mesh_size=1.)
    end
    for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,1)))
        gmsh.model.geo.addLine(a,b,tag)
        Tessella.Model.add_line!(native,a,b;tag)
        gmsh.model.geo.mesh.setTransfiniteCurve(tag,count)
        Tessella.Model.set_transfinite_curve!(native,tag,count)
    end
    gmsh.model.geo.addCurveLoop(loop,1);gmsh.model.geo.addPlaneSurface([1],1)
    Tessella.Model.add_curve_loop!(native,loop;tag=1)
    Tessella.Model.add_plane_surface!(native,[1];tag=1)
    gmsh.model.geo.mesh.setTransfiniteSurface(1,"Left",corners)
    Tessella.Model.set_transfinite_surface!(native,1,"Left",corners)
    if get(case,:recombine,false)
        gmsh.model.geo.mesh.setRecombine(2,1)
        Tessella.Model.set_recombine!(native,2,1)
    end
    gmsh.model.geo.synchronize()
    if get(case,:reverse,false)
        gmsh.model.mesh.setReverse(2,1,true)
        Tessella.Model.set_reverse!(native,2,1,true)
    end
    gmsh.model.mesh.generate(2)
    mesh=Tessella.Model.mesh_model_surface(native,1)
    normal=get(case,:axis,:xy)===:xy ?
        (0.,0.,Float64(case.direction)) : (0.,Float64(-case.direction),0.)
    gmsh.model.getNormal(1,[0.,0.])==collect(normal) || error(
        "Gmsh public triangle CAD normal differs from the signed analytic plane")
    Tessella.Model.model_normal(native,1,[0.,0.])==collect(normal) || error(
        "native public triangle CAD normal differs from the signed analytic plane")
    expected=get(case,:reverse,false) ? -1 : 1
    tags,xyz,_=gmsh.model.mesh.getNodes()
    coordinates=reshape(xyz,3,:)
    positions=Dict(tag=>index for (index,tag) in enumerate(tags))
    types,_,families=gmsh.model.mesh.getElements(2,1)
    oracle_counts=Dict{Int,Int}()
    for (msh,nodes) in zip(types,families)
        width=gmsh.model.mesh.getElementProperties(msh)[4]
        cells=reshape(nodes,width,:)
        oracle_counts[msh]=size(cells,2)
        for cell in eachcol(cells)
            exact_public_cell_orientation(coordinates,[positions[tag] for tag in cell],normal)==expected ||
                error("Gmsh public triangle winding disagrees with CAD/Reverse on $case")
        end
    end
    blocks=mesh isa Tessella.MeshTypes.Mesh ?
        (Tessella.Elements.ElementBlock(2,mesh.tris),) : mesh.blocks
    native_counts=Dict{Int,Int}()
    for block in blocks
        native_counts[block.msh]=size(block.nodes,2)
        for cell in eachcol(block.nodes)
            exact_public_cell_orientation(mesh.coords,cell,normal)==expected ||
                error("native public triangle winding disagrees with CAD/Reverse on $case")
        end
    end
    native_counts==oracle_counts || error(
        "public triangle family/count mismatch on $case: $native_counts != $oracle_counts")
    return sum(values(native_counts))
end

function gmsh_to_tessella_node_map(mesh, surface)
    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    tags = sort!(unique(reduce(vcat, element_nodes; init=UInt64[])))
    tessella_nodes = size(mesh.coords, 2)
    length(tags) == tessella_nodes || error(
        "node-count mismatch: Gmsh $(length(tags)), Tessella $tessella_nodes")
    used = falses(tessella_nodes)
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
        for destination in 1:tessella_nodes
            used[destination] && continue
            candidate = (mesh.coords[1, destination],
                         mesh.coords[2, destination],
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

function check_recombined_arrangement(
    arrangement, symbol, count=6; tilted=false)
    curves, surface = add_curved_triangle(
        arrangement, count; tilted, recombine=true)
    sides = map(curve_points, curves)
    all(length(side) == count for side in sides) || error(
        "Gmsh did not emit $count nodes on each transfinite curve")
    mesh = mesh_transfinite_triangle_patch(
        sides...; arrangement=symbol, face_tag=21,
        side_tags=(11, 12, 13))
    Tessella.Elements.validate(mesh).ok || error(
        "Tessella recombined $arrangement patch did not validate")
    divisions = count - 1
    expected = (count * (count + 1) ÷ 2, 3divisions, divisions,
                divisions * (divisions - 1) ÷ 2)
    line_block = only(filter(block -> block.msh == 1, mesh.blocks))
    triangle_block = only(filter(block -> block.msh == 2, mesh.blocks))
    quadrangle_blocks = filter(block -> block.msh == 3, mesh.blocks)
    quadrangle_count = isempty(quadrangle_blocks) ? 0 :
                        size(only(quadrangle_blocks).nodes, 2)
    (size(mesh.coords, 2), size(line_block.nodes, 2),
     size(triangle_block.nodes, 2), quadrangle_count) == expected || error(
        "unexpected Tessella recombined $arrangement counts")
    mapping, maximum_error = gmsh_to_tessella_node_map(mesh, surface)

    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    triangle_position = findfirst(==(Int32(2)), types)
    triangle_position === nothing && error(
        "Gmsh recombined $arrangement emitted no first-order triangles")
    mapped_triangles = Int32[
        mapping[tag] for tag in element_nodes[triangle_position]]
    tessella_triangles = canonical_triangles(vec(triangle_block.nodes))
    canonical_triangles(mapped_triangles) == tessella_triangles || error(
        "$arrangement recombined triangle connectivity differs from Gmsh")

    quadrangle_position = findfirst(==(Int32(3)), types)
    mapped_quadrangles = if quadrangle_position === nothing
        Int32[]
    else
        Int32[mapping[tag] for tag in element_nodes[quadrangle_position]]
    end
    tessella_quadrangles = isempty(quadrangle_blocks) ? NTuple{4,Int32}[] :
        canonical_quadrangles(vec(only(quadrangle_blocks).nodes))
    canonical_quadrangles(mapped_quadrangles) == tessella_quadrangles || error(
        "$arrangement quadrangle connectivity differs from Gmsh")

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
        canonical_segments(vec(line_block.nodes)) || error(
            "$arrangement recombined boundary connectivity differs from Gmsh")
    return maximum_error, tessella_triangles, tessella_quadrangles,
           Tessella.Elements.mixed_crc(mesh).sha
end

# Differential for Gmsh's legacy `Mesh.TransfiniteTri = 0` collapsed-grid
# algorithm: `mesh_transfinite_triangle_collapsed` against the same Gmsh
# geometry, with equal counts, an auto-rotated unequal-count boundary, and an
# explicitly pinned collapsed corner.
function check_collapsed(arrangement, symbol, count=6; tilted=false,
                         counts=nothing, corners=nothing, rotate=0)
    curves, surface = add_curved_triangle(
        arrangement, count; tilted=tilted, counts=counts, corners=corners)
    sides = map(curve_points, curves)
    rotate > 0 && (sides = Tuple(sides[mod1(index + rotate, 3)]
                               for index in 1:3))
    mesh = mesh_transfinite_triangle_collapsed(
        sides...; arrangement=symbol, face_tag=21, side_tags=(11, 12, 13),
        allow_corner_rotation=corners isa AbstractVector && isempty(corners))
    validate(mesh).ok || error(
        "Tessella collapsed $arrangement patch did not validate")
    mapping, maximum_error = gmsh_to_tessella_node_map(mesh, surface)

    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    types == Int32[2] || error(
        "Gmsh collapsed $arrangement emitted element types $types")
    mapped_triangles = Int32[
        mapping[tag] for tag in element_nodes[1]]
    canonical_triangles(mapped_triangles) ==
        canonical_triangles(vec(mesh.tris)) || error(
            "collapsed $arrangement triangle connectivity differs from Gmsh")

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
            "collapsed $arrangement boundary connectivity differs from Gmsh")
    return maximum_error
end

# Recombined arm of the same collapsed algorithm (`Mesh.TransfiniteTri = 0`
# plus `Recombine Surface`): upstream keeps the apex fan triangular and emits
# one quadrangle per remaining grid cell, independent of the arrangement.
function check_collapsed_recombined(arrangement, symbol, count=6;
                                    tilted=false, counts=nothing,
                                    corners=nothing, rotate=0)
    curves, surface = add_curved_triangle(
        arrangement, count; tilted=tilted, counts=counts, corners=corners,
        recombine=true)
    sides = map(curve_points, curves)
    rotate > 0 && (sides = Tuple(sides[mod1(index + rotate, 3)]
                               for index in 1:3))
    mesh = mesh_transfinite_triangle_collapsed_patch(
        sides...; arrangement=symbol, face_tag=21, side_tags=(11, 12, 13),
        allow_corner_rotation=corners isa AbstractVector && isempty(corners))
    Tessella.Elements.validate(mesh).ok || error(
        "Tessella recombined collapsed $arrangement patch did not validate")
    mapping, maximum_error = gmsh_to_tessella_node_map(mesh, surface)

    line_block = only(filter(block -> block.msh == 1, mesh.blocks))
    triangle_block = only(filter(block -> block.msh == 2, mesh.blocks))
    quadrangle_blocks = filter(block -> block.msh == 3, mesh.blocks)
    quadrangle_count = isempty(quadrangle_blocks) ? 0 :
                        size(only(quadrangle_blocks).nodes, 2)
    types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    triangle_position = findfirst(==(Int32(2)), types)
    triangle_position === nothing && error(
        "Gmsh recombined collapsed $arrangement emitted no triangles")
    size(triangle_block.nodes, 2) == 0 && error(
        "Tessella recombined collapsed $arrangement emitted no fan triangles")
    mapped_triangles = Int32[
        mapping[tag] for tag in element_nodes[triangle_position]]
    canonical_triangles(mapped_triangles) ==
        canonical_triangles(vec(triangle_block.nodes)) || error(
            "recombined collapsed $arrangement triangle connectivity " *
            "differs from Gmsh")
    length(element_nodes[triangle_position]) ==
        3size(triangle_block.nodes, 2) || error(
            "recombined collapsed $arrangement fan count differs from Gmsh")

    quadrangle_position = findfirst(==(Int32(3)), types)
    mapped_quadrangles = quadrangle_position === nothing ? Int32[] :
        Int32[mapping[tag] for tag in element_nodes[quadrangle_position]]
    tessella_quadrangles = isempty(quadrangle_blocks) ? NTuple{4,Int32}[] :
        canonical_quadrangles(vec(only(quadrangle_blocks).nodes))
    canonical_quadrangles(mapped_quadrangles) == tessella_quadrangles ||
        error("recombined collapsed $arrangement quadrangle connectivity " *
              "differs from Gmsh")

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
        canonical_segments(vec(line_block.nodes)) || error(
            "recombined collapsed $arrangement boundary connectivity " *
            "differs from Gmsh")
    return maximum_error, canonical_triangles(mapped_triangles),
           tessella_quadrangles, Tessella.Elements.mixed_crc(mesh).sha
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

    # Exercise the minimal patch and three additional triangular-lattice sizes;
    # arrangement parity above uses the nonuniform six-node case.
    for count in (2, 3, 8, 10)
        maximum_error, _ = check_arrangement("Left", :left, count)
        push!(errors, maximum_error)
        coordinate_samples += count * (count + 1) ÷ 2
    end
    maximum_error, _ = check_arrangement("Left", :left, 6; tilted=true)
    push!(errors, maximum_error)
    coordinate_samples += 21

    recombined_errors = Float64[]
    recombined_topologies = Tuple{
        Vector{NTuple{3,Int32}},Vector{NTuple{4,Int32}}}[]
    recombined_crcs = String[]
    recombined_coordinate_samples = 0
    for (arrangement, symbol) in cases
        maximum_error, triangles, quadrangles, crc =
            check_recombined_arrangement(arrangement, symbol)
        push!(recombined_errors, maximum_error)
        push!(recombined_topologies, (triangles, quadrangles))
        push!(recombined_crcs, crc)
        recombined_coordinate_samples += 21
    end
    recombined_topologies[3] == recombined_topologies[4] || error(
        "Gmsh AlternateLeft/AlternateRight recombined layouts diverged")
    recombined_topologies[1] != recombined_topologies[2] || error(
        "Gmsh Left/Right recombined layouts unexpectedly coincide")
    for count in (2, 3, 8, 10)
        maximum_error, _, _, _ =
            check_recombined_arrangement("Left", :left, count)
        push!(recombined_errors, maximum_error)
        recombined_coordinate_samples += count * (count + 1) ÷ 2
    end
    maximum_error, _, _, _ = check_recombined_arrangement(
        "Left", :left, 6; tilted=true)
    push!(recombined_errors, maximum_error)
    recombined_coordinate_samples += 21

    # Legacy collapsed-quadrilateral algorithm (Gmsh default
    # Mesh.TransfiniteTri=0): same curved geometries, all four arrangements,
    # several resolutions, one tilted plane, the auto-rotated unequal-count
    # boundary (5,5,8), and an explicitly pinned second collapsed corner.
    gmsh.option.setNumber("Mesh.TransfiniteTri", 0)
    collapsed_errors = Float64[]
    collapsed_topologies = Vector{NTuple{3,Int32}}[]
    collapsed_samples = 0
    for (arrangement, symbol) in cases
        curves, surface = add_curved_triangle(
            arrangement, 6; counts=(6, 6, 6), corners=Int[])
        sides = map(curve_points, curves)
        mesh = mesh_transfinite_triangle_collapsed(
            sides...; arrangement=symbol, face_tag=21,
            side_tags=(11, 12, 13))
        validate(mesh).ok || error(
            "Tessella collapsed $arrangement patch did not validate")
        mapping, maximum_error = gmsh_to_tessella_node_map(mesh, surface)
        types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
        types == Int32[2] || error(
            "Gmsh collapsed $arrangement emitted element types $types")
        topology = canonical_triangles(
            Int32[mapping[tag] for tag in element_nodes[1]])
        topology == canonical_triangles(vec(mesh.tris)) || error(
            "collapsed $arrangement triangle connectivity differs from Gmsh")
        push!(collapsed_errors, maximum_error)
        push!(collapsed_topologies, topology)
        collapsed_samples += nnodes(mesh)
    end
    for count in (2, 4, 7)
        maximum_error = check_collapsed("Left", :left, count)
        push!(collapsed_errors, maximum_error)
        collapsed_samples += 1 + (count - 1) * count
    end
    push!(collapsed_errors,
          check_collapsed("Left", :left, 6; tilted=true))
    collapsed_samples += 31
    # (5,5,8) node counts: Gmsh auto-rotates the collapsed corner to the
    # second junction — the kernel does the same.
    push!(collapsed_errors, check_collapsed(
        "Left", :left, 0; counts=(5, 5, 8), corners=Int[]))
    collapsed_samples += 33
    # Pin the collapsed corner at the third vertex: Gmsh gets corner list
    # (p3,p1,p2) on counts (6,8,8); the kernel gets the same boundary rotated
    # so the pinned junction leads, with rotation disabled.
    push!(collapsed_errors, check_collapsed(
        "Left", :left, 0; counts=(6, 8, 8), rotate=2,
        corners=[7, 1, 4]))
    collapsed_samples += 43

    # Recombined collapsed arm: same geometries under `Recombine Surface`.
    # Upstream's recombined collapsed emission is arrangement-independent —
    # the branch precedes the diagonal dispatch — so Tessella's canonical
    # cells and CRC must be identical across all four arrangements.
    collapsed_recombined_errors = Float64[]
    collapsed_recombined_topologies = Tuple{
        Vector{NTuple{3,Int32}},Vector{NTuple{4,Int32}}}[]
    collapsed_recombined_crcs = String[]
    collapsed_recombined_samples = 0
    for (arrangement, symbol) in cases
        maximum_error, triangles, quadrangles, crc =
            check_collapsed_recombined(arrangement, symbol)
        push!(collapsed_recombined_errors, maximum_error)
        push!(collapsed_recombined_topologies, (triangles, quadrangles))
        push!(collapsed_recombined_crcs, crc)
        collapsed_recombined_samples += 31
    end
    all(topology == collapsed_recombined_topologies[1]
        for topology in collapsed_recombined_topologies) || error(
        "recombined collapsed topology changed with arrangement")
    all(crc == collapsed_recombined_crcs[1]
        for crc in collapsed_recombined_crcs) || error(
        "recombined collapsed mesh changed with arrangement")
    for count in (2, 4, 7)
        maximum_error, _, _, _ = check_collapsed_recombined(
            "Left", :left, count)
        push!(collapsed_recombined_errors, maximum_error)
        collapsed_recombined_samples += 1 + (count - 1) * count
    end
    maximum_error, _, _, _ = check_collapsed_recombined(
        "Left", :left, 6; tilted=true)
    push!(collapsed_recombined_errors, maximum_error)
    collapsed_recombined_samples += 31
    push!(collapsed_recombined_errors, check_collapsed_recombined(
        "Left", :left, 0; counts=(5, 5, 8), corners=Int[])[1])
    collapsed_recombined_samples += 33
    push!(collapsed_recombined_errors, check_collapsed_recombined(
        "Left", :left, 0; counts=(6, 8, 8), rotate=2,
        corners=[7, 1, 4])[1])
    collapsed_recombined_samples += 43
    orientation_cases=public_orientation_cases()
    orientation_cells=sum(check_public_orientation,orientation_cases)
    println("TRANSFINITE_TRIANGLE_DIFFERENTIAL_OK gmsh=$api_version " *
            "arrangements=$(length(cases)) resolutions=5 geometries=2 " *
            "coordinate_samples=$coordinate_samples " *
            "max_node_error=$(maximum(errors)) reference_nodes=21 " *
            "reference_segments=15 reference_triangles=25 " *
            "arrangement_topologies=1 recombined_arrangements=$(length(cases)) " *
            "recombined_resolutions=5 recombined_geometries=2 " *
            "recombined_coordinate_samples=$recombined_coordinate_samples " *
            "recombined_max_node_error=$(maximum(recombined_errors)) " *
            "recombined_reference_triangles=5 " *
            "recombined_reference_quadrangles=10 " *
            "recombined_alternate_topologies=1 " *
            "recombined_crcs=$(join(recombined_crcs, ',')) " *
            "collapsed_arrangements=$(length(cases)) " *
            "collapsed_coordinate_samples=$collapsed_samples " *
            "collapsed_max_node_error=$(maximum(collapsed_errors)) " *
            "collapsed_recombined_arrangements=$(length(cases)) " *
            "collapsed_recombined_samples=$collapsed_recombined_samples " *
            "collapsed_recombined_max_node_error=" *
            "$(maximum(collapsed_recombined_errors)) " *
            "collapsed_recombined_reference_triangles=5 " *
            "collapsed_recombined_reference_quadrangles=20 " *
            "public_orientation_cases=$(length(orientation_cases)) " *
            "public_orientation_cells=$orientation_cells")
finally
    gmsh.finalize()
end
