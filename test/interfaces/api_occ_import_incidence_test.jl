module OCCImportIncidenceChecks
using Test, Tessella
const A=Tessella.API
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

@testset "OCC import returns lower topology occurrences when requested" begin
    for (name,counts) in (("brep_occ_box.step",(48,24,6,1)),
                          ("brep_iges_msbo_tet.igs",(24,12,4,1)),
                          ("iges_occ_curves.igs",(11,5,0,0)))
        A.initialize()
        try
            A.option("General.Terminal",0)
            point=A.model.occ.add_point(-7.,-8.,-9.)
            previous=Set(A.model.occ.get_entities())
            out=A.model.occ.import_shapes(joinpath(FIXTURES,name);highestDimOnly=false)
            @test Tuple(count(dt->first(dt)==d,out) for d in 0:3)==counts
            @test issorted(first.(out);rev=true)
            @test !((0,point) in out)
            @test Set(out)==setdiff(Set(A.model.occ.get_entities()),previous)
            A.model.occ.synchronize()
            @test all(dt->dt in A.model.get_entities(),out)
        finally
            A.finalize()
        end
    end
end
end
