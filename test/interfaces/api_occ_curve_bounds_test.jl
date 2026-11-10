module OCCCurveBoundsChecks
using Test,Tessella
const A=Tessella.API

@testset "Copied CAD curve bounds include interior geometry" begin
    for family in (:add_bezier,:add_bspline)
        A.initialize()
        try
            occ=A.model.occ
            points=[occ.add_point(x,y,z) for (x,y,z) in
                ((0.,0.,0.),(.5,2.,1.),(1.,0.,0.))]
            curve=getproperty(occ,family)(points)
            occ.synchronize()
            for transform in (false,true)
                if transform
                    occ.translate([(1,curve)],3.,-2.,4.)
                    occ.synchronize()
                end
                bounds=A.model.get_bounding_box(1,curve)
                @test bounds==occ.get_bounding_box(1,curve)
                lo,hi=A.model.get_parametrization_bounds(1,curve)
                parameters=collect(range(only(lo),only(hi);length=33))
                xyz=reshape(A.model.get_value(1,curve,parameters),3,:)
                @test all(bounds[axis]<=xyz[axis,row]<=bounds[axis+3]
                    for axis in 1:3 for row in axes(xyz,2))
                @test bounds[5]>bounds[2]
                @test bounds[6]>bounds[3]
            end
        finally
            A.finalize()
        end
    end
end
end
