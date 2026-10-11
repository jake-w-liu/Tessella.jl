module DiscreteCurveGeometryTests
using Test, Tessella
include(joinpath(@__DIR__,"..","..","validation","native_geometry","discrete_curve_provider_primary_v1.jl"))
const _DC_API=Tessella.API
_dc_close(actual,expected)=isapprox(actual,Float64.(expected);atol=2e-11,rtol=0)
function _dc_seed(middle,connectivity,kind=1;supplied=false,labels=[101,102,103],owner_order=[1,2,11])
    _DC_API.initialize()
    _DC_API.option("Mesh.Renumber",0)
    _DC_API.model.add_discrete_entity(0,1)
    _DC_API.model.add_discrete_entity(0,2)
    _DC_API.model.add_discrete_entity(1,11,[1,2])
    for owner in owner_order
        node,coordinate=owner==1 ? (labels[1],[0.,0.,0.]) :
                        owner==2 ? (labels[2],[1.,0.,0.]) : (labels[3],middle)
        _DC_API.mesh.add_nodes(owner==11 ? 1 : 0,owner,[node],coordinate,
                              owner==11 && supplied ? [.7] : Float64[])
    end
    kind==8 && _DC_API.mesh.add_nodes(1,11,[104,105],[9.,8.,7.,-9.,-8.,-7.])
    _DC_API.mesh.add_elements_by_type(1,15,[501],[labels[1]])
    _DC_API.mesh.add_elements_by_type(2,15,[502],[labels[2]])
    _DC_API.mesh.add_elements_by_type(11,kind,[601,602],connectivity)
end
function _dc_check(stage)
    bounds=_DC_API.model.get_parametrization_bounds(1,11)
    @test _dc_close(bounds[1],stage["bounds"]["value"][1])
    @test _dc_close(bounds[2],stage["bounds"]["value"][2])
    for sample in stage["samples"]
        parameter=sample["parameter"]
        @test _dc_close(_DC_API.model.get_value(1,11,[parameter]),sample["point"]["value"])
        @test _dc_close(_DC_API.model.get_derivative(1,11,[parameter]),sample["derivative"]["value"])
        inverse=sample["inverse"]
        inverse===nothing || @test _dc_close(_DC_API.model.get_parametrization(
            1,11,sample["point"]["value"]),inverse["value"])
    end
    for node in stage["nodes"]
        coordinates,parameters,_,_=_DC_API.mesh.get_node(node["tag"])
        @test _dc_close(coordinates,node["coordinates"])
        @test _dc_close(parameters,node["parameters"])
    end
end
@testset "Discrete Curve persistent chain-index geometry" begin
    for row in DISCRETE_CURVE_PROVIDER_ROWS
        @testset "$(row["name"])" begin
            _dc_seed(row["middle"],row["connectivity"],row["kind"];
                     supplied=row["supplied_original_U"])
            try
                _DC_API.mesh.create_geometry([(1,11)])
                _dc_check(row["stages"][1])
                _DC_API.mesh.set_node(103,[.4,.9,.2],row["supplied_original_U"] ? [.7] : Float64[])
                _dc_check(row["stages"][2])
                _DC_API.mesh.create_geometry([(1,11)])
                _dc_check(row["stages"][3])
                _DC_API.mesh.clear([(1,11)])
                _dc_check(row["stages"][4])
            finally
                _DC_API.finalize()
            end
        end
    end
    for (index,row) in enumerate(DISCRETE_CURVE_CHAIN_ROWS)
        @testset "Line sequence and sparse fallback $index" begin
            _dc_seed([.25,.5,0.],row["connectivity"];labels=row["labels"],owner_order=row["add_nodes_owner_order"])
            try
                _DC_API.mesh.create_geometry([(1,11)])
                for sample in row["samples"]
                    @test _dc_close(_DC_API.model.get_value(1,11,[sample["parameter"]]),sample["point"])
                    @test _dc_close(_DC_API.model.get_derivative(1,11,[sample["parameter"]]),sample["derivative"])
                end
            finally
                _DC_API.finalize()
            end
        end
    end
end

end
