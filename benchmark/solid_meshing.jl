# Solid meshing benchmark, Tessella side: unstructured tetrahedral meshes of
# CAD solids through the session API. `solid_meshing_gmsh.py` runs the same
# cases through Gmsh's Python API; see README.md for how to run both.
#
#     julia --project=. --startup-file=no -O2 --threads=1 benchmark/solid_meshing.jl [results.toml]
#     julia --project=. --startup-file=no -O2 --threads=1 benchmark/solid_meshing.jl --faces [results.toml]
#
# Every case is meshed twice in the same process and the second run is
# reported, so compilation is excluded. Times are wall-clock seconds for
# `generate(2)` (curves and surfaces) and for the following `generate(3)`.
# `--faces` runs the face-count cases instead: an n-by-n grid of unit plane
# faces at size 0.2, timed for `generate(1)` and the following `generate(2)`.

using Tessella
using Printf

const A = Tessella.API

const CASES = [
    ("box", () -> A.model.occ.add_box(0, 0, 0, 1, 1, 1), (0.1, 0.05, 0.03, 0.02)),
    ("sphere", () -> A.model.occ.add_sphere(0, 0, 0, 1), (0.2, 0.1, 0.06, 0.04)),
    ("cylinder", () -> A.model.occ.add_cylinder(0, 0, 0, 0, 0, 2, 0.5), (0.1, 0.05, 0.03)),
    ("torus", () -> A.model.occ.add_torus(0, 0, 0, 2, 0.5), (0.15, 0.1, 0.07)),
    ("drilled_block", () -> begin
        block = A.model.occ.add_box(0, 0, 0, 2, 2, 2)
        drill = A.model.occ.add_cylinder(1, 1, -1, 0, 0, 4, 0.5)
        A.model.occ.cut([(3, block)], [(3, drill)])
    end, (0.15, 0.08, 0.05)),
]

function tetrahedra()
    types, tags, nodes = A.mesh.get_elements(3)
    for i in eachindex(types)
        types[i] == 4 && return tags[i], nodes[i]
    end
    return UInt64[], UInt64[]
end

function measure(build, size)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        build()
        A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin", size)
        A.option("Mesh.MeshSizeMax", size)
        surface = @elapsed A.mesh.generate(2)
        volume = @elapsed A.mesh.generate(3)
        node_tags, coordinates, _ = A.mesh.get_nodes()
        position = Dict(tag => index for (index, tag) in enumerate(node_tags))
        tags, connectivity = tetrahedra()
        triangle_count = 0
        types, element_tags, _ = A.mesh.get_elements(2)
        for i in eachindex(types)
            types[i] == 2 && (triangle_count += length(element_tags[i]))
        end
        total = 0.0
        for k in 1:4:length(connectivity)
            a, b, c, d = (position[connectivity[k+j]] for j in 0:3)
            p(i) = (coordinates[3i-2], coordinates[3i-1], coordinates[3i])
            pa, pb, pc, pd = p(a), p(b), p(c), p(d)
            u = pb .- pa; v = pc .- pa; w = pd .- pa
            total += (u[1] * (v[2] * w[3] - v[3] * w[2]) - u[2] * (v[1] * w[3] - v[3] * w[1]) +
                      u[3] * (v[1] * w[2] - v[2] * w[1])) / 6
        end
        quality = isempty(tags) ? Float64[] : A.mesh.get_element_qualities(tags, "gamma")
        return (nodes=length(node_tags), triangles=triangle_count, tetrahedra=length(tags),
                volume=total, gamma_min=isempty(quality) ? NaN : minimum(quality),
                gamma_mean=isempty(quality) ? NaN : sum(quality) / length(quality),
                surface_seconds=surface, volume_seconds=volume)
    finally
        A.finalize()
    end
end

function plane_grid(n)
    geo = A.model.geo
    point(i, j) = i * (n + 1) + j + 1
    for i in 0:n, j in 0:n
        geo.add_point(Float64(i), Float64(j), 0.0; tag=point(i, j))
    end
    along(i, j) = i * n + j + 1
    across(i, j) = (n + 1) * n + i * (n + 1) + j + 1
    for i in 0:n, j in 0:n-1
        geo.add_line(point(i, j), point(i, j + 1); tag=along(i, j))
    end
    for i in 0:n-1, j in 0:n
        geo.add_line(point(i, j), point(i + 1, j); tag=across(i, j))
    end
    face = 0
    for i in 0:n-1, j in 0:n-1
        face += 1
        geo.add_curve_loop([along(i, j), across(i, j + 1), -along(i + 1, j), -across(i, j)]; tag=face)
        geo.add_plane_surface([face]; tag=face)
    end
    geo.synchronize()
    return nothing
end

function measure_faces(n)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        plane_grid(n)
        A.option("Mesh.MeshSizeMin", 0.2)
        A.option("Mesh.MeshSizeMax", 0.2)
        curves = @elapsed A.mesh.generate(1)
        surfaces = @elapsed A.mesh.generate(2)
        node_tags, _, _ = A.mesh.get_nodes()
        _, element_tags, _ = A.mesh.get_elements(2)
        return (nodes=length(node_tags), triangles=sum(length, element_tags; init=0),
                curve_seconds=curves, surface_seconds=surfaces)
    finally
        A.finalize()
    end
end

function main(arguments)
    faces = !isempty(arguments) && arguments[1] == "--faces"
    faces && (arguments = arguments[2:end])
    io = isempty(arguments) ? stdout : open(arguments[1], "w")
    try
        println(io, "implementation = \"Tessella $(pkgversion(Tessella))\"")
        println(io, "julia = \"$(VERSION)\"")
        println(io, "threads = $(Threads.nthreads())")
        if faces
            for n in (5, 10, 20)
                measure_faces(n)
                result = measure_faces(n)
                println(io, "\n[[grid]]")
                println(io, "faces = $(n^2)")
                println(io, "nodes = $(result.nodes)")
                println(io, "triangles = $(result.triangles)")
                @printf(io, "curve_seconds = %.4f\n", result.curve_seconds)
                @printf(io, "surface_seconds = %.4f\n", result.surface_seconds)
                flush(io)
            end
            return
        end
        for (name, build, sizes) in CASES, size in sizes
            measure(build, size)                    # compile and warm
            result = measure(build, size)
            println(io, "\n[[case]]")
            println(io, "solid = \"$name\"")
            println(io, "size = $size")
            println(io, "nodes = $(result.nodes)")
            println(io, "triangles = $(result.triangles)")
            println(io, "tetrahedra = $(result.tetrahedra)")
            @printf(io, "volume = %.6f\n", result.volume)
            @printf(io, "gamma_min = %.4f\n", result.gamma_min)
            @printf(io, "gamma_mean = %.4f\n", result.gamma_mean)
            @printf(io, "surface_seconds = %.4f\n", result.surface_seconds)
            @printf(io, "volume_seconds = %.4f\n", result.volume_seconds)
            flush(io)
        end
    finally
        io === stdout || close(io)
    end
end

main(ARGS)
