module QuadTriNoNewGmshInputsTest
# QuadTriNoNewVerts inputs the bounded catalogs used to refuse and Gmsh 4.15.2
# meshes. Expected values: nodes and surface + volume cells by type, from 20
# fresh Gmsh processes per input, all alike (Gmsh's QuadTri compares vertex
# addresses, so a spread is possible; none was observed on these).
using Test, Tessella
using Tessella.Elements: MixedMesh, msh_dimension
using Tessella.MeshTypes: Mesh, nnodes

const _TRIANGLE = """
Mesh.TransfiniteTri=1;
Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0,1,0,1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
Curve Loop(1)={1,2,3};Plane Surface(1)={1};
Transfinite Curve{:}=2;Transfinite Surface{1};
sweep[]=Extrude{0.0,0.0,1.0}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;};
"""

function _square(n, extra="")
    return """
    Geometry.AutoCoherence=1;
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={1,1,0,1};Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{1,2,3,4}=$n;Transfinite Surface{1} Left;$extra
    sweep[]=Extrude{0.0,0.0,1.0}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;};
    """
end

function _run(source)
    mktempdir() do dir
        path = joinpath(dir, "nonew.geo")
        write(path, source)
        return execute_geo(path; mesh_dim=3)
    end
end

function _counts(mesh)
    counts = Dict{Int,Int}()
    if mesh isa Mesh
        size(mesh.tris, 2) > 0 && (counts[2] = size(mesh.tris, 2))
        size(mesh.tets, 2) > 0 && (counts[4] = size(mesh.tets, 2))
    else
        for block in mesh.blocks
            msh_dimension(block.msh) in (2, 3) || continue
            counts[Int(block.msh)] = get(counts, Int(block.msh), 0) + size(block.nodes, 2)
        end
    end
    return (nnodes(mesh), counts)
end

@testset "QuadTriNoNewVerts inputs Gmsh meshes" begin
    # A rotation about the Source's own normal: flat tetrahedra, as Gmsh
    # emits them (the Point on the axis stays one node).
    rotated = replace(_TRIANGLE, "Extrude{0.0,0.0,1.0}" => "Extrude{{0,0,1},{0,0,0},Pi/6}")
    @test _counts(_run(rotated).mesh) == (9, Dict(2 => 14, 4 => 6))
    # Transfinite cap: meshGFace runs MeshTransfiniteSurface before
    # MeshExtrudedSurface; a one-triangle cap equals the copy.
    @test _counts(_run(_TRIANGLE * "Transfinite Surface{sweep[0]};\n").mesh) ==
          (12, Dict(2 => 20, 4 => 9))
    # One Source swept by two Regions, one on each side.
    both = _TRIANGLE * "Extrude{0,0,-1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;}\n"
    @test _counts(_run(both).mesh) == (21, Dict(2 => 39, 4 => 18))
    # Two triangles rotated in their plane: flat prisms (the arc centre is a
    # Point node of its own, as in Gmsh).
    flat = replace(_square(2), "QuadTriNoNewVerts" => "QuadTriNoNewVerts RecombLaterals",
        "Extrude{0.0,0.0,1.0}" => "Extrude{{0,0,1},{-1,0,0},Pi/6}")
    @test _counts(_run(flat).mesh) == (17, Dict(2 => 4, 3 => 12, 6 => 6))
    # The 2-by-2 quadrangle Source rotated about an axis off the Source goes
    # through the general QuadToTri planner.
    patch = _square(3, "Recombine Surface{1};")
    turned = replace(patch, "Extrude{0.0,0.0,1.0}" => "Extrude{{0,1,0},{-1,0,0},Pi/6}")
    @test _counts(_run(turned).mesh) == (38, Dict(2 => 56, 3 => 4, 4 => 64, 7 => 4))
    # The same Source swept by two Regions.
    shared = patch * "Extrude{0,0,-1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;}\n"
    @test _counts(_run(shared).mesh) == (63, Dict(2 => 112, 3 => 4, 4 => 128, 7 => 8))
    # Two square Regions sharing a transfinite wall: Gmsh meshes the wall
    # transfinite, then QuadToTriLateralRemesh re-triangulates it with the
    # Regions' diagonals; 20 runs gave the two results below.
    squares = """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Point(5)={2,0,0,1};Point(6)={2,1,0,1};
    Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
    Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
    Transfinite Curve{:}=2;
    Transfinite Surface{1};Transfinite Surface{2};
    Recombine Surface{1,2};
    left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
    right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
    Transfinite Surface{20};
    """
    @test _counts(_run(squares).mesh) in ((18, Dict(2 => 32, 3 => 2, 4 => 20, 7 => 2)),
                                          (18, Dict(2 => 32, 3 => 2, 4 => 16, 7 => 4)))
    # A transfinite cap whose diagonals the Region cannot follow: Gmsh stops
    # ("Extruded hexahedron needs subdivision ..."), so does Tessella.
    @test_throws r"do not exactly cover the typed volume boundary" _run(
        patch * "Transfinite Surface{sweep[0]};\n")
    # A reversed top copy places its nodes at the swept generatrix vertices:
    # Gmsh's node at x = 0.4999999999986921 (the inverted parameter gave
    # 0.49999999999869216, one ulp off, and broke the column bijection).
    down = replace(replace(patch, "Curve Loop(1)={1,2,3,4}" => "Curve Loop(1)={-4,-3,-2,-1}"),
        "Extrude{0.0,0.0,1.0}" => "Extrude{0.0,0.0,-1.5}", "Layers{3}" => "Layers{1}")
    mesh = _run(down).mesh
    @test _counts(mesh) == (18, Dict(2 => 24, 3 => 4, 4 => 16, 7 => 4))
    bottom = Set((mesh.coords[1, i], mesh.coords[2, i]) for i in axes(mesh.coords, 2) if mesh.coords[3, i] == 0.0)
    top = Set((mesh.coords[1, i], mesh.coords[2, i]) for i in axes(mesh.coords, 2) if mesh.coords[3, i] == -1.5)
    @test top == bottom && (0.4999999999986921, 0.0) in top
end

end
