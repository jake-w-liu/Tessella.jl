module APIHomologyEmptyDomainNativeTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
@testset "homology API Cint safety and strict Model queue validation remain atomic" begin
    A.initialize()
    try
        model=A._model_locked()
        for dims in ([true],[Int64(typemax(Cint))+1],[Int64(typemin(Cint))-1],[1.0])
            @test_throws ArgumentError A.mesh.add_homology_request("Homology",Int[],Int[],dims)
            @test isempty(model.meshing.homology_requests)
        end
        @test_throws ArgumentError A.mesh.add_homology_request(:Homology,Int[],Int[],[0])
        @test isempty(model.meshing.homology_requests)
        for (kind,dims) in (("Unknown",[0]),("Homology",[-1]),("Homology",[4]))
            @test_throws ArgumentError Tessella.Model.add_homology_request!(model;kind,dims)
            @test isempty(model.meshing.homology_requests)
        end
    finally A.finalize() end
end

end
