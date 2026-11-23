module GeoDiscreteParamT13Test

using Test, Tessella

const API = Tessella.API

# `CreateGeometry` on a classified STL-like surface ports Gmsh 4.15.2's
# discreteFace::createGeometry: the harmonic disk map of computeParametrization
# (mean-value weights, longest boundary loop on the unit circle), kept apart
# from the mesh with its own triangles, and the mean-plane coordinates of
# param::checkPlanar for a flat patch. generate(1) then deletes the patch's
# triangles (it is no longer fully discrete) and generate(2) meshes it on the
# stored parametrization, as tutorials/t13.geo does. Every expected value was
# read from the Gmsh 4.15.2 Python wheel.

function classified(kind)
    API.initialize()
    API.option("General.Terminal", 0)
    API.model.add_discrete_entity(2, 1)
    n = 4
    tags = Int[]; xyz = Float64[]
    for j in 0:n, i in 0:n
        if kind === :cylinder
            a = 0.5 * pi * i / n
            append!(xyz, (cos(a), sin(a), j / n))
        else
            append!(xyz, (i / n, j / n, 0.0))
        end
        push!(tags, 1 + i + (n + 1) * j)
    end
    conn = Int[]
    for j in 0:n-1, i in 0:n-1
        a = 1 + i + (n + 1) * j; b = a + 1; c = b + n + 1; d = a + n + 1
        append!(conn, (a, b, c, a, c, d))
    end
    API.mesh.add_nodes(2, 1, tags, xyz)
    API.mesh.add_elements(2, 1, [2], [collect(1:length(conn)÷3)], [conn])
    API.mesh.classify_surfaces(pi / 2, true, false, pi)
    API.mesh.create_geometry()
    return API.model.get_entities(2)[1][2]
end

function node_uv(s, node)
    coord, _ = API.mesh.get_node(node)
    return API.model.get_parametrization(2, s, coord)
end

counts() = [sum(length, API.mesh.get_elements(k)[2]; init=0) for k in 0:2]

@testset "a curved patch gets Gmsh's harmonic parametrization" begin
    try
        s = classified(:cylinder)
        @test API.model.get_entities() == [(0, 1), (1, 3), (2, 2)]
        lo, hi = API.model.get_parametrization_bounds(2, s)
        @test lo == [-1.0, -1.0] && hi == [1.0, 1.0]
        # computeParametrization's (u,v) at mesh nodes (PETSc LU in Gmsh,
        # dense LU here: equal to rounding).
        for (node, uv) in ((1, (1.0000000000000002, 0.0)),
                           (6, (0.9533321816523397, -0.3019234198037482)),
                           (13, (-2.3849088703857873e-17, -1.8188804982411413e-16)),
                           (25, (-1.0000000000000002, -3.865129321289331e-16)))
            @test isapprox(collect(node_uv(s, node)), collect(uv); atol=1e-13)
        end
        @test isapprox(API.model.get_value(2, s, [0.1, 0.2]),
            [0.6857706625742989, 0.7213631198593646, 0.3831629963822848]; atol=1e-12)
        @test isapprox(API.model.get_derivative(2, s, [0.1, 0.2]),
            [0.5561629376175452, -0.37161619411848734, -0.3310099195527826,
             -0.3847620618700165, 0.25708979042332863, -0.41868005831218424]; atol=1e-10)
        API.option("Mesh.MeshSizeMax", 0.25)
        API.mesh.generate(1)
        # The parametrized patch loses its 32 input triangles.
        @test counts() == [1, 52, 0]
        box = API.model.get_bounding_box(2, s)
        @test isapprox(collect(box), [6.123233995736766e-17, 0.0, 0.0, 1.0, 1.0, 1.0]; atol=1e-14)
        API.mesh.generate(2)
        # Gmsh: 242 nodes, 430 triangles (the patch's coordinates went
        # through cos/sin and an LU solve: exact on Windows x86-64 only).
        c = counts()
        @test c[1:2] == [1, 52]
        if Sys.iswindows() && Sys.ARCH === :x86_64
            @test c[3] == 430
        else
            @test abs(c[3] - 430) <= 22
        end
    finally
        API.finalize()
    end
end

@testset "a flat patch is meshed in its mean-plane coordinates" begin
    try
        s = classified(:plate)
        lo, hi = API.model.get_parametrization_bounds(2, s)
        @test isapprox(lo, [0.0, 0.0]; atol=1e-15) && isapprox(hi, [1.0, 1.0]; atol=1e-15)
        for (node, uv) in ((1, (0.0, 0.0)), (6, (0.25, 0.0)), (13, (0.5, 0.5)),
                           (25, (1.0, 1.0)))
            @test isapprox(collect(node_uv(s, node)), collect(uv); atol=1e-14)
        end
        @test isapprox(API.model.get_value(2, s, [0.1, 0.2]), [0.2, 0.1, 0.0]; atol=1e-14)
        API.option("Mesh.MeshSizeMax", 0.25)
        API.mesh.generate(1)
        @test counts() == [1, 40, 0]
        API.mesh.generate(2)
        # Gmsh: 143 nodes, 244 triangles; the curve nodes at exactly 0.1
        # (meshGEdge's arithmetic on discreteEdge::firstDer).
        c = counts()
        @test c[1:2] == [1, 40]
        if Sys.ARCH === :x86_64
            @test c[3] == 244
        else
            @test abs(c[3] - 244) <= 6
        end
        xyz = API.mesh.get_nodes(1, 3)[2]
        @test xyz[1:3] == [0.1, 0.0, 0.0]
    finally
        API.finalize()
    end
end

end # module
