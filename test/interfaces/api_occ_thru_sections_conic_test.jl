module OCCThruSectionsConicChecks
using Test,Tessella
const A=Tessella.API

# `occ.add_thru_sections` over circular sections threw
# `FieldError: type _ModelOCCNURBSCurve has no field reparam`: the section
# edges are converted to plain rational quadratics, whose angle chart has to
# be built from the edge's own range (thrusections.geo, neuron.geo, pipe.geo,
# twist.geo and t19.geo all stopped there). A ruled loft between two circles
# is an oblique cone frustum, volume πh(r1²+r1r2+r2²)/3 — Gmsh 4.15.2 gives
# the same 0.408407 and the same (2,3,3,1) entities.
@testset "OCC thru sections over circles" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        bottom=occ.add_curve_loop([occ.add_circle(0,0,0,0.5)])
        top=occ.add_curve_loop([occ.add_circle(0.1,0.05,1,0.2)])
        out=occ.add_thru_sections([bottom,top];makeSolid=true,makeRuled=true)
        @test out==[(3,1)]
        @test occ.get_mass(3,1)≈π*(0.25+0.1+0.04)/3 rtol=1e-9
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,3,1)
    finally
        A.finalize()
    end
end
end
