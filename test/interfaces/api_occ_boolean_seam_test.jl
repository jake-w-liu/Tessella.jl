module OCCBooleanSeamChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const M=Tessella.Model

@testset "Crossed-cylinder union preserves chart seam topology" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_cylinder(0.,0.,0.,0.,0.,1.,.3;tag=1)
        occ.add_cylinder(.2,0.,0.,0.,1.,0.,.3;tag=2)
        out,images=occ.fuse([(3,1)],[(3,2)];tag=3)
        @test out==[(3,3)]
        @test images==[Tuple{Int,Int}[],Tuple{Int,Int}[]]
        @test Tuple(length(occ.get_entities(d)) for d in 0:3)==(7,11,6,1)
        # Opaque public Gmsh 4.15.2 control for these two finite cylinders.
        @test occ.get_mass(3,3)≈0.5411659937831585 rtol=2e-7
        # Boundary first moments across plane caps and trimmed cylinder
        # side faces (Gmsh 4.15.2 public control for the same union).
        @test collect(occ.get_center_of_mass(3,3))≈
            collect((0.0999999868,0.2553812701,0.2553813159)) atol=2e-7
        # Centroidal inertia: conic plane caps and cylinder side moments.
        @test occ.get_matrix_of_inertia(3,3)≈
            [0.12948926895796575,-0.014187942418562827,0.014187940619115194,
             -0.014187942418562827,0.08288799931771688,0.035699963204375136,
             0.014187940619115194,0.035699963204375136,
             0.08288800851412019] atol=2e-7
        occ.synchronize()
        faces=A.model.get_boundary(out,false,true,false)
        @test length(faces)==6
        @test isempty(A.model.get_boundary(faces,true,true,false))
        model=A.CURRENT[]
        # Fitted sections retain the native 1e-7 tracing accuracy; seam
        # anchors must not create a missing/duplicate pcurve or open wire.
        for (surface,loops) in model.surfaces
            geometry=model.surface_geometry[surface]
            for loop in loops,signed_curve in model.loops[loop]
                curve=abs(signed_curve)
                entry=geometry.pcurves[curve]
                pc=signed_curve<0 && entry.rev!==nothing ? entry.rev : entry.fwd
                @test pc!==nothing
                lo,hi=A.model.get_parametrization_bounds(1,curve)
                for fraction in (0.,.25,.5,.75,1.)
                    t=lo[1]+fraction*(hi[1]-lo[1])
                    uv=M._brep_pc_eval(pc,t)
                    p=A.model.get_value(1,curve,[t])
                    q=A.model.get_value(2,surface,collect(uv))
                    @test norm(p-q)<=1e-7
                end
            end
        end
    finally
        A.finalize()
    end
end
end
