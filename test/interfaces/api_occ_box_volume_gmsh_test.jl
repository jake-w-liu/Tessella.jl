# The suite includes every test file into Main; this file keeps its aliases to itself.
module OCCBoxVolumeGmshTests
using Test
using Tessella
const A = Tessella.API

# A Box is filled from its face meshes by Delaunay insertion, as Gmsh fills it.
# It was once filled by bisecting edges down to the size field, which on a unit
# Box at size 0.1 gave 29,120 tetrahedra where Gmsh 4.15.2 has 4,979, and when
# `occ.add_box` began to store the same records as the `.geo` Box that fill
# reached it too. The counts below are Gmsh 4.15.2's for
# `occ.addBox` + `Mesh.MeshSizeMin = Mesh.MeshSizeMax = h` + `generate(3)`.

count_type(dim, type) = (result = A.mesh.get_elements(dim);
    sum((length(result[2][i]) for i in eachindex(result[1]) if result[1][i] == type); init=0))

function counts(build, h)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        build()
        A.option("Mesh.MeshSizeMin", h); A.option("Mesh.MeshSizeMax", h)
        A.mesh.generate(3)
        return (length(A.mesh.get_nodes()[1]), count_type(2, 2), count_type(3, 4))
    finally
        A.finalize()
    end
end

occ_box(dx, dy, dz) = () -> (A.model.occ.add_box(0, 0, 0, dx, dy, dz); A.model.occ.synchronize())
function geo_box()
    path = joinpath(mktempdir(), "box.geo")
    write(path, "SetFactory(\"OpenCASCADE\");\nBox(1) = {0,0,0,1,1,1};\n")
    A.open(path)
end
function two_boxes()
    A.model.occ.add_box(0, 0, 0, 1, 1, 1); A.model.occ.add_box(1, 0, 0, 1, 1, 1)
    A.model.occ.synchronize()
end

@testset "Box volumes have Gmsh's nodes, triangles and tetrahedra" begin
    @test counts(occ_box(1, 1, 1), 0.3) == (143, 264, 387)
    @test counts(occ_box(1, 1, 1), 0.25) == (144, 264, 391)
    # The .geo Box's curve nodes differ from Gmsh's by about 1e-12 (0.2499999999994109
    # where Gmsh has 0.2499999999999997), so its volume is not Gmsh's input. With
    # Gmsh's TetGen pass and insertion order it gives 381 tetrahedra; 387 was a
    # coincidence of the earlier fill. Gmsh: (143, 264, 387).
    @test_broken counts(geo_box, 0.3) == (143, 264, 387)
    @test counts(occ_box(1, 2, 3), 0.4) == (223, 400, 635)
    @test counts(two_boxes, 0.3) == (286, 528, 773)
end

@testset "A Box is not a factory Volume" begin
    A.initialize()
    try
        A.model.occ.add_box(0, 0, 0, 1, 1, 1); A.model.occ.synchronize()
        model = A.CURRENT[]
        @test !Tessella.Model._model_factory_volume_encoding(model, 1)
    finally
        A.finalize()
    end
end
end
