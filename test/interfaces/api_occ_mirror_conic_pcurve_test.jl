module OCCMirrorConicPcurveChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# A reflected circle traces its cap plane left-handed and is stored as a
# `neg2d` wrapper. Four constructors wrote `(neg2d=(pc,))` — a parenthesized
# assignment that evaluates to the bare 1-tuple `(pc,)` — instead of the
# NamedTuple `(neg2d=(pc,),)`, so every later `haskey(pc,:neg2d)` threw a
# MethodError: `occ.get_mass` on a mirrored cylinder or cone crashed.
@testset "OCC mirrored conic pcurves" begin
    inner=(circ2d=((0.0,0.0),(1.0,0.0),0.3),)
    for wrapped in (M._model_occ_sew_shift_pc(inner,-1.0,0.0,1e-9,"test"),
                    M._occ_pcurve_reparam((neg2d=((lin2d=((0.,0.),(1.,0.)),),),),
                                          0.5,2.0))
        @test wrapped isa NamedTuple
        @test haskey(wrapped,:neg2d)
        @test wrapped.neg2d isa Tuple
    end
    for (make,volume,area) in (
            (occ->occ.add_cylinder(.5,.5,-.5,0,0,2,.3),π*.09*2,
                2π*.3*2+2π*.09),
            (occ->occ.add_cone(0,0,0,0,0,1,.5,.2),π*(.25+.1+.04)/3,
                π*.7*hypot(1,.3)+π*(.25+.04)))
        for plane in ((1.,0.,0.,.5),(1.,1.,.5,.25))
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                tag=make(occ)
                occ.mirror([(3,tag)],plane...)
                occ.synchronize()
                @test occ.get_mass(3,tag)≈volume rtol=1e-12
                @test sum(occ.get_mass(2,t)
                          for (_,t) in A.model.get_entities(2))≈area rtol=1e-12
            finally
                A.finalize()
            end
        end
    end
end
end
