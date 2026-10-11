module DiscreteCurveP2GeometryTests
include(joinpath(@__DIR__,"..","..","validation","native_geometry","discrete_curve_geometry_helpers_v1.jl"))
include("../artifacts/discrete_curve_closed_primary_v1.jl")
include("../artifacts/discrete_curve_queries_lifetime_primary_v1.jl")

function _dc_nonfinite_close(actual,expected)
    return length(actual)==length(expected) && all(
        isnan(expected[index]) ? isnan(actual[index]) :
        isapprox(actual[index],expected[index];atol=2e-11,rtol=0) for index in eachindex(expected))
end

@testset "Actual discrete Curve inverse, closest and finite-difference queries" begin
    for row in DISCRETE_CURVE_QUERY_ROWS
        _DC_API.initialize();_DC_API.option("Mesh.Renumber",0)
        try
            _DC_API.model.add_discrete_entity(0,1)
            _DC_API.model.add_discrete_entity(0,2)
            _DC_API.model.add_discrete_entity(1,11,[1,2])
            points=row["points"];labels=row["labels"]
            _DC_API.mesh.add_nodes(0,1,[labels[1]],points[1])
            _DC_API.mesh.add_nodes(0,2,[labels[end]],points[end])
            _DC_API.mesh.add_nodes(1,11,labels[2:end-1],reduce(vcat,points[2:end-1]))
            connectivity=reduce(vcat,[[labels[i],labels[i+1]] for i in 1:length(labels)-1])
            _DC_API.mesh.add_elements_by_type(11,1,collect(601:599+length(labels)),connectivity)
            _DC_API.mesh.create_geometry([(1,11)])
            for sample in row["derivatives"]
                parameter=[sample["u"]]
                @test _dc_close(_DC_API.model.get_derivative(1,11,parameter),sample["first"])
                @test _dc_close(_DC_API.model.get_second_derivative(1,11,parameter),sample["second"])
                @test _dc_nonfinite_close(_DC_API.model.get_curvature(1,11,parameter),sample["curvature"])
            end
            for query in row["queries"]
                @test _dc_close(_DC_API.model.get_parametrization(1,11,query["point"]),query["inverse"])
                point,parameter=_DC_API.model.get_closest_point(1,11,query["point"])
                @test _dc_close(point,query["closest_point"])
                @test _dc_close(parameter,query["closest_parameter"])
            end
        finally
            _DC_API.finalize()
        end
    end
end

function _dc_check_lifetime_stage(stage)
    curve=stage["curve"]
    bounds=_DC_API.model.get_parametrization_bounds(1,curve)
    @test _dc_close(bounds[1],stage["bounds"][1])
    @test _dc_close(bounds[2],stage["bounds"][2])
    for sample in stage["samples"]
        @test _dc_close(_DC_API.model.get_value(1,curve,[sample["u"]]),sample["point"])
        @test _dc_close(_DC_API.model.get_derivative(1,curve,[sample["u"]]),sample["derivative"])
    end
    @test sort(_DC_API.mesh.get_nodes()[1])==sort(UInt64[node["tag"] for node in stage["nodes"]])
    for node in stage["nodes"]
        point,parameter,dimension,owner=_DC_API.mesh.get_node(node["tag"])
        @test _dc_close(point,node["coordinates"])
        @test _dc_close(parameter,node["parameters"])
        @test [dimension,owner]==node["owner"]
    end
    kinds,tags,nodes=_DC_API.mesh.get_elements(1,curve)
    @test kinds==stage["types"]
    @test tags==stage["labels"]
    @test nodes==stage["connectivity"]
end

@testset "Discrete Curve snapshot across public label, topology and mesh lifecycles" begin
    for row in DISCRETE_CURVE_LIFETIME_ROWS
        _dc_seed([.25,.5,0.],[101,103,103,102];supplied=true)
        try
            _DC_API.mesh.create_geometry([(1,11)])
            _dc_check_lifetime_stage(row["stages"][1])
            curve=11;name=row["name"]
            if name=="rekey"
                _DC_API.model.set_tag(1,11,71);curve=71
            elseif name=="node_rename"
                _DC_API.mesh.renumber_nodes([101,102,103],[900,7,400])
            elseif name=="element_rename"
                _DC_API.mesh.renumber_elements([501,502,601,602],[9001,7001,4001,2001])
            elseif name=="clear"
                _DC_API.mesh.clear([(1,11)])
            elseif name=="reverse"
                _DC_API.mesh.reverse([(1,11)])
            elseif name=="remove_recreate"
                _DC_API.model.remove_entities([(1,11)],false)
                _DC_API.model.add_discrete_entity(1,11,[1,2])
                _DC_API.mesh.add_nodes(1,11,[103],[.75,-.25,0.],[.3])
                _DC_API.mesh.add_elements_by_type(11,1,[701,702],[101,103,103,102])
                _DC_API.mesh.create_geometry([(1,11)])
            end
            _dc_check_lifetime_stage(row["stages"][2])
            _DC_API.mesh.create_geometry([(1,curve)])
            _dc_check_lifetime_stage(row["stages"][3])
        finally
            _DC_API.finalize()
        end
    end
end

function _dc_check_raw_stage(stage)
    @test sort(_DC_API.mesh.get_nodes()[1])==sort(UInt64[node["tag"] for node in stage["nodes"]])
    for node in stage["nodes"]
        coordinates,parameters,dimension,owner=_DC_API.mesh.get_node(node["tag"])
        @test _dc_close(coordinates,node["coordinates"])
        @test _dc_close(parameters,node["parameters"])
        @test [dimension,owner]==node["owner"]
    end
    for entity in stage["entities"]
        dimension,owner=entity["entity"]
        tags,coordinates,parameters=_DC_API.mesh.get_nodes(dimension,owner,false,true)
        @test tags==entity["node_tags"]
        @test _dc_close(coordinates,entity["coordinates"])
        @test _dc_close(parameters,entity["parameters"])
        kinds,labels,connectivity=_DC_API.mesh.get_elements(dimension,owner)
        @test kinds==entity["types"]
        @test labels==entity["element_tags"]
        @test connectivity==entity["connectivity"]
    end
end

@testset "Independent full raw Curve state and P2 geometry supports" begin
    for row in DISCRETE_CURVE_STATE_ROWS
        @testset "geometry=$(row["geometry"]) supplied=$(row["supplied"]) order=$(row["order"])" begin
            _dc_seed([.25,0.,0.],[101,103,103,102];supplied=row["supplied"])
            _DC_API.option("Mesh.MeshOnlyEmpty",1)
            try
                for stage in row["stages"]
                    stage["stage"]=="create_geometry" && _DC_API.mesh.create_geometry([(1,11)])
                    stage["stage"]=="set_order2" && _DC_API.mesh.set_order(2)
                    stage["stage"]=="generate2_onlyempty" && _DC_API.mesh.generate(2)
                    _dc_check_raw_stage(stage)
                end
            finally
                _DC_API.finalize()
            end
        end
    end
end

function _dc_seed_admission(row)
    _DC_API.initialize();_DC_API.option("Mesh.Renumber",0)
    for point in (1,2)
        _DC_API.model.add_discrete_entity(0,point)
    end
    _DC_API.model.add_discrete_entity(1,11,row["declared_boundary"])
    if row["name"] in ("no_boundary_curve_owned","declared_boundary_curve_owned")
        _DC_API.mesh.add_nodes(1,11,[101,102,103],[0.,0.,0.,1.,0.,0.,.25,.5,0.])
    else
        _DC_API.mesh.add_nodes(0,1,[101],[0.,0.,0.])
        _DC_API.mesh.add_nodes(0,2,[102],[1.,0.,0.])
        _DC_API.mesh.add_nodes(1,11,[103,104],[.25,.5,0.,.75,.5,0.])
    end
    connectivity=row["connectivity"]
    isempty(connectivity) || _DC_API.mesh.add_elements_by_type(
        11,1,collect(601:600+length(connectivity)÷2),connectivity)
end

@testset "Closed discrete Curve, topology admission and captured curvature" begin
    for row in DISCRETE_CURVE_CLOSED_ROWS
        @testset "$(row["name"])" begin
            _dc_seed_admission(row)
            try
                before_nodes=_DC_API.mesh.get_nodes()
                before_cells=_DC_API.mesh.get_elements()
                if row["created"]["ok"]
                    _DC_API.mesh.create_geometry([(1,11)])
                    @test _DC_API.model.get_boundary([(1,11)],false,true,false)==
                        Tuple{Int,Int}[Tuple(pair) for pair in row["boundary_after"]]
                    for sample in row["samples"]
                        u=sample["parameter"]
                        @test _dc_close(_DC_API.model.get_value(1,11,[u]),sample["value"])
                        @test _dc_close(_DC_API.model.get_derivative(1,11,[u]),sample["derivative"])
                        curvature=_DC_API.model.get_curvature(1,11,[u])
                        expected=sample["curvature"]
                        @test all(isnan(expected[i]) ? isnan(curvature[i]) :
                            isapprox(curvature[i],expected[i];atol=2e-11,rtol=0) for i in eachindex(expected))
                    end
                    if row["order2_nodes"]!==nothing
                        _DC_API.mesh.set_order(2)
                        @test sort(_DC_API.mesh.get_nodes()[1])==
                            sort(UInt64[node["tag"] for node in row["order2_nodes"]])
                        for node in row["order2_nodes"]
                            coordinates,parameters,_,_=_DC_API.mesh.get_node(node["tag"])
                            @test _dc_close(coordinates,node["coordinates"])
                            @test _dc_close(parameters,node["parameters"])
                        end
                    end
                else
                    @test_throws ArgumentError _DC_API.mesh.create_geometry([(1,11)])
                    @test _DC_API.mesh.get_nodes()==before_nodes
                    @test _DC_API.mesh.get_elements()==before_cells
                end
            finally
                _DC_API.finalize()
            end
        end
    end
end

using Test,Tessella
const RAW_CURVE_API=Tessella.API
include("../artifacts/discrete_curve_raw_parameter_rename_primary_v1.jl")
include("../artifacts/discrete_curve_folded_support_primary_v1.jl")
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
        RAW_CURVE_API.mesh.set_order(2)
        for node in RAW_CURVE_PARAMETER_RENAME_ROWS[end]["nodes"]
            actual=RAW_CURVE_API.mesh.get_node(node["tag"])[2]
            @test isapprox(actual,Float64.(node["parameters"]);atol=2e-11,rtol=0)
        end
        stage=FOLDED_DISCRETE_CURVE_STAGE
        kinds,tags,nodes=RAW_CURVE_API.mesh.get_elements(1,11)
        @test kinds==stage["types"]
        @test tags==stage["tags"]
        @test nodes==stage["connectivity"]
        for (tag,point) in stage["coordinates"]
            @test isapprox(RAW_CURVE_API.mesh.get_node(parse(Int,tag))[1],Float64.(point);atol=2e-11,rtol=0)
        end
        for sample in stage["jacobians"]
            jacobian,determinant,point=RAW_CURVE_API.mesh.get_jacobian(stage["folded_element"],[sample["u"],0.,0.])
            @test isapprox(jacobian,Float64.(sample["jacobian"]);atol=2e-11,rtol=0)
            @test isapprox(determinant,Float64.(sample["determinant"]);atol=2e-11,rtol=0)
            @test isapprox(point,Float64.(sample["coordinates"]);atol=2e-11,rtol=0)
        end
    finally
        RAW_CURVE_API.finalize()
    end
end

using Test,Tessella
include("../artifacts/discrete_curve_line8_jacobian_primary_v1.jl")
@testset "Public Line3 weighted Jacobians: folded, zero, spatial and translated" begin
    api=Tessella.API
    api.initialize()
    try
        for row in LINE8_JACOBIAN_ROWS
            coordinates=hcat((Float64.(point) for point in row["points"])...)
            mesh=Tessella.Elements.MixedMesh(coordinates,
                [Tessella.Elements.ElementBlock(8,reshape(Int32[1,2,3],3,1))])
            # Existing classification-free Mixed query entry, as exercised by
            # api_mixed_queries_test. The captured folded case also runs through
            # addNodes/createGeometry/setOrder in the preceding public test.
            lock(api.STATE_LOCK) do
                api._replace_mesh_cache_locked!(api._copy_mesh(mesh))
            end
            local_coordinates=Float64[x for u in row["parameters"] for x in (u,0,0)]
            result=api.mesh.get_jacobian(1,local_coordinates)
            for (actual,key) in zip(result,("jacobian","determinant","coordinates"))
                @test isapprox(actual,Float64.(row[key]);atol=2e-11,rtol=0)
            end
            @test api.mesh.get_jacobians(8,local_coordinates)==result
            @test mesh.coords==coordinates
        end
    finally
        api.finalize()
    end
end

end
