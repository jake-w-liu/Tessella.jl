module OCCBooleanSpherePoleChecks
using Test,Tessella
const A=Tessella.API

# A box corner at the sphere center: the three box faces through the corner
# cut the sphere along its equator, its seam meridian and a second meridian
# that runs into the pole. The kept spherical piece is the octant bounded by
# those three arcs and a quarter of the pole's degenerate edge. Counts are
# Gmsh 4.15.2's; the volume is the closed form.
@testset "OCC Boolean sections through a sphere pole" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(0,0,0,2,2,2)
        ball=occ.add_sphere(0,0,0,1.5)
        out,images=occ.intersect([(3,box)],[(3,ball)])
        @test out==[(3,1)]
        @test images==[[(3,1)],[(3,1)]]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(4,7,4,1)
        @test occ.get_mass(3,1)≈π*1.5^3/6 rtol=1e-8
        @test collect(occ.get_center_of_mass(3,1))≈fill(0.5625,3) atol=1e-7
        faces=A.model.get_boundary([(3,1)],false,false,false)
        @test length(faces)==4
        # the spherical face keeps a quarter of the pole edge: 3 arcs + it
        spherical=[t for (_,t) in faces if A.model.get_type(2,t)=="Sphere"]
        @test length(spherical)==1
        edges=A.model.get_boundary([(2,only(spherical))],false,false,false)
        @test length(edges)==4
        spans=[begin
                   lo,hi=A.model.get_parametrization_bounds(1,t)
                   hi[1]-lo[1]
               end for (_,t) in edges if occ.get_mass(1,t)==0]
        @test spans≈[π/2] atol=1e-12
        # every edge of the solid is shared by exactly two faces
        uses=Dict{Int,Int}()
        for face in faces, (_,t) in A.model.get_boundary([face],false,false,false)
            occ.get_mass(1,t)==0 && continue
            uses[t]=get(uses,t,0)+1
        end
        @test length(uses)==6
        @test all(==(2),values(uses))
    finally
        A.finalize()
    end
end
end
