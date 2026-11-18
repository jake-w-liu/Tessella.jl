module APIOCCCurvatureAdjacentFaceSizeChecks
using Test, Tessella
const A=Tessella.API

# `Mesh.MeshSizeFromCurvature` on a curve: upstream `LC_MVertex_CURV` takes
# the larger of the curve's own curvature and `max_surf_curvature`, the
# largest curvature of the faces the curve bounds. The second term was
# missing, so a straight seam of a cylinder ignored the option (7 cells where
# Gmsh 4.15.2 places 11) and a torus seam circle followed its own large
# radius instead of the tube radius (56 cells instead of 60).
function line_cells(build,elements_per_turn)
    A.initialize()
    try
        build();A.model.occ.synchronize()
        A.option("Mesh.MeshSizeFromCurvature",elements_per_turn)
        A.option("Mesh.MeshSizeMax",10)
        A.mesh.generate(1)
        counts=Int[]
        for (_,curve) in A.model.get_entities(1)
            _,tags,_=A.mesh.get_elements(1,curve)
            push!(counts,sum(length,tags;init=0))
        end
        return sort!(counts)
    finally
        A.finalize()
    end
end

@testset "curvature sizing of a curve includes its adjacent faces" begin
    # Cylinder r=0.5, h=2, 16 elements per turn: size 2π·0.5/16 on the two
    # rim circles (16 cells each) and on the straight seam (2/0.196 → 11).
    cylinder=line_cells(()->A.model.occ.add_cylinder(0,0,0,0,0,2,0.5),16)
    @test length(cylinder)==3
    @test cylinder[2:3]==[16,16]
    @test 10<=cylinder[1]<=12                 # Gmsh: 11; was 7
    # Torus R=2, r=0.5, 12 elements per turn: the tube circle gets 12 cells
    # and the outer seam circle (radius 2.5) follows the tube curvature 1/r:
    # 2π·2.5/(2π·0.5/12) = 60 cells.
    torus=line_cells(()->A.model.occ.add_torus(0,0,0,2,0.5),12)
    @test torus[1]==12
    @test 59<=torus[2]<=61                    # Gmsh: 60; was 56
    # Plane faces contribute no curvature: a Box keeps its straight edges.
    box=line_cells(()->A.model.occ.add_box(0,0,0,1,1,1),16)
    @test length(box)==12 && all(==(box[1]),box)
end
end
