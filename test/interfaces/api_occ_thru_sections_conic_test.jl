module OCCThruSectionsConicChecks
using Test,Tessella
const A=Tessella.API

# `occ.add_thru_sections` over circular sections threw
# `FieldError: type _ModelOCCNURBSCurve has no field reparam`: the section
# edges are converted to plain rational quadratics, whose angle chart has to
# be built from the edge's own range (thrusections.geo, neuron.geo, pipe.geo,
# twist.geo and t19.geo all stopped there). A ruled loft between two circles
# is an oblique cone frustum with Gmsh 4.15.2's (2,3,3,1) entities; its
# volume πh(r1²+r1r2+r2²)/3 = 0.408407 holds too, but integrating the lofted
# patch takes 19 s interpreted, so the two caps stand in for it.
@testset "OCC thru sections over circles" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        bottom=occ.add_curve_loop([occ.add_circle(0,0,0,0.5)])
        top=occ.add_curve_loop([occ.add_circle(0.1,0.05,1,0.2)])
        out=occ.add_thru_sections([bottom,top];makeSolid=true,makeRuled=true)
        @test out==[(3,1)]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,3,1)
        caps=sort!([occ.get_mass(2,t) for (_,t) in A.model.get_entities(2)
                    if A.model.get_type(2,t)=="Plane"])
        @test caps≈[π*0.04,π*0.25] rtol=1e-12
    finally
        A.finalize()
    end
end
end
