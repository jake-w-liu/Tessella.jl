module DiscreteCurveP2GeometryTests
using Test,Tessella
const RAW_CURVE_API=Tessella.API

@testset "Raw heterogeneous original parameters survive public node rename" begin
    RAW_CURVE_API.initialize();RAW_CURVE_API.option("Mesh.Renumber",0)
    try
        for point in (1,2)
            RAW_CURVE_API.model.add_discrete_entity(0,point)
        end
        RAW_CURVE_API.model.add_discrete_entity(1,11,[1,2])
        RAW_CURVE_API.mesh.add_nodes(0,1,[101],[0.,0.,0.])
        RAW_CURVE_API.mesh.add_nodes(0,2,[102],[1.,0.,0.])
        RAW_CURVE_API.mesh.add_nodes(1,11,[103],[.25,.5,0.])
        RAW_CURVE_API.mesh.add_nodes(1,11,[104],[.75,.5,0.],[1.4])
        RAW_CURVE_API.mesh.add_elements_by_type(11,1,[601,602,603],[101,103,103,104,104,102])
        RAW_CURVE_API.mesh.create_geometry([(1,11)])
        @test isempty(RAW_CURVE_API.mesh.get_node(103)[2])
        @test RAW_CURVE_API.mesh.get_node(104)[2]==[1.4]
        RAW_CURVE_API.mesh.renumber_nodes([101,102,103,104],[900,7,400,300])
        @test isempty(RAW_CURVE_API.mesh.get_node(400)[2])
        @test RAW_CURVE_API.mesh.get_node(300)[2]==[1.4]
    finally
        RAW_CURVE_API.finalize()
    end
end
end
