module GeoExtendContextFieldTest

using Test, Tessella

# `Field[] = Extend` reads the listed curves' mesh lines when the first
# surface is meshed (upstream `ExtendField::recomputeCurves`). Both the
# `.geo` executor and the API session used to reject the kind. Gmsh 4.15.2
# meshes this square with 200 nodes / 324 triangles and a mean edge length
# of 0.036 for the triangles within 0.2 of the 0.02-sized bottom curve.
const SQUARE = """
    Point(1)={0,0,0,0.02}; Point(2)={1,0,0,0.02}; Point(3)={1,1,0,0.3}; Point(4)={0,1,0,0.3};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    Mesh.MeshSizeExtendFromBoundary=0;
    """

function generate(fields::AbstractString)
    path = joinpath(mktempdir(), "extend.geo")
    write(path, SQUARE * fields)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(2)
        tags, xyz, _ = A.mesh.get_nodes()
        pos = Dict(tags[i] => (xyz[3i-2], xyz[3i-1]) for i in eachindex(tags))
        _, _, nodes = A.mesh.get_elements(2)
        tri = nodes[1]
        near = Float64[]
        for k in 1:3:length(tri)
            p = [pos[tri[k+j]] for j in 0:2]
            sum(q[2] for q in p) / 3 < 0.2 || continue
            push!(near, sum(hypot((p[j] .- p[mod1(j + 1, 3)])...) for j in 1:3) / 3)
        end
        return length(tags), length(tri) ÷ 3, sum(near) / length(near)
    finally
        A.finalize()
    end
end

@testset "Extend field sizes from the curve mesh" begin
    nodes, triangles, near = generate("""
        Field[1]=Extend; Field[1].CurvesList={1}; Field[1].DistMax=0.5;
        Field[1].SizeMax=0.3; Background Field=1;
        """)
    @test 195 <= nodes <= 205          # Gmsh 4.15.2: 200
    @test 315 <= triangles <= 335      # Gmsh 4.15.2: 324
    @test isapprox(near, 0.036; atol=0.003)
    # Without the field the interior is only boundary-driven (75 triangles
    # near the curve, mean edge 0.066, identical in both programs).
    _, _, plain = generate("")
    @test plain > 0.06
    @test_throws Exception generate("""
        Field[1]=Extend; Field[1].SurfacesList={1}; Background Field=1;
        """)
end

end # module
