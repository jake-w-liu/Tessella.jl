module OCCBooleanAbuttingChecks
using Test,Tessella
const A=Tessella.API

# Two boxes sharing a whole face: the fragment keeps one copy of the shared
# face and of its four edges. The second box's side faces then reference
# edges owned by the first box, whose curves run the other way — their
# pcurves have to be re-expressed in the surviving edge's parameter, or the
# face's wire encloses a wrong area (5 instead of 4) and the solid a wrong
# volume (8⅓ instead of 8).
@testset "OCC fragment of abutting boxes" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        left=occ.add_box(-1,-1,-1,2,2,2)
        right=occ.add_box(1,-1,-1,2,2,2)
        out,images=occ.fragment([(3,left)],[(3,right)])
        @test out==[(3,1),(3,2)]
        @test images==[[(3,1)],[(3,2)]]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(12,20,11,2)
        @test occ.get_mass(3,1)≈8 rtol=1e-12
        @test occ.get_mass(3,2)≈8 rtol=1e-12
        @test collect(occ.get_center_of_mass(3,1))≈[0.0,0.0,0.0] atol=1e-12
        @test collect(occ.get_center_of_mass(3,2))≈[2.0,0.0,0.0] atol=1e-12
        @test all(t->occ.get_mass(2,t)≈4,last.(A.model.get_entities(2)))
        # each stored pcurve traces its edge: face(pcurve(t)) == edge(t)
        model=A.CURRENT[]
        M=Tessella.Model
        for (surface,loops) in model.surfaces
            geometry=model.surface_geometry[surface]
            for loop in loops,signed_curve in model.loops[loop]
                curve=abs(signed_curve)
                entry=geometry.pcurves[curve]
                pc=signed_curve<0 && entry.rev!==nothing ? entry.rev : entry.fwd
                lo,hi=A.model.get_parametrization_bounds(1,curve)
                for fraction in (0.,.5,1.)
                    t=lo[1]+fraction*(hi[1]-lo[1])
                    uv=M._brep_pc_eval(pc,t)
                    p=A.model.get_value(1,curve,[t])
                    q=A.model.get_value(2,surface,collect(uv))
                    @test maximum(abs,p-q)<=1e-12
                end
            end
        end
        A.option("Mesh.MeshSizeMax",1.0)
        A.mesh.generate(2)
        @test length(A.mesh.get_nodes()[1])>=12
    finally
        A.finalize()
    end
end
end
