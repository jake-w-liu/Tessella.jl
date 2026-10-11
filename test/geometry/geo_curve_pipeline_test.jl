module SecondOrderGeoCommandsProbeV1
using Tessella,Test
const A=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","geo_curve_pipeline","fixtures.jl"))
function check_result(row)
    expected=row["result"]
    @test Set(A.model.get_entities())==Set(Tuple(pair) for pair in expected["entities"])
    for entity in expected["nodes"]
        dimension=entity["dimension"];owner=entity["owner"]
        tags,coordinates,parameters=A.mesh.get_nodes(dimension,owner,false,true)
        @test tags==UInt64.(entity["tags"])
        @test isapprox(coordinates,Float64.(entity["coordinates"]);atol=2e-11,rtol=0)
        @test isapprox(parameters,Float64.(entity["query_parameters"]);atol=2e-11,rtol=0)
        for node in entity["nodes"]
            point,parameter,dim,tag=A.mesh.get_node(node["tag"])
            @test (dim,tag)==(node["dimension"],node["owner"])
            @test isapprox(point,Float64.(node["coordinates"]);atol=2e-11,rtol=0)
            @test isapprox(parameter,Float64.(node["parameters"]);atol=2e-11,rtol=0)
            if haskey(row,"raw") && row["linear"]==1 && node["tag"]>maximum(v["tag"] for e in row["raw"]["nodes"] for v in e["nodes"])
                @test reinterpret(UInt64,point)==reinterpret(UInt64,Float64.(node["coordinates"]))
                @test isempty(parameter)
            end
        end
    end
    for entity in expected["cells"]
        types,tags,connectivity=A.mesh.get_elements(entity["dimension"],entity["owner"])
        @test types==entity["types"]
        @test tags==[UInt64.(values) for values in entity["tags"]]
        @test connectivity==[UInt64.(values) for values in entity["connectivity"]]
    end
    @test (A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==Tuple(expected["maxima"])
end
@testset "GEO SecondOrderLinear generation, SetOrder and RefineMesh public state" begin
    for row in SECONDORDER_GEO_ROWS
        @testset "$(row["recipe"]) linear=$(row["linear"]) $(row["operation"])" begin
            A.initialize()
            try
                mktemp() do path,io
                    write(io,row["geo"]);close(io);A.open_geo!(path)
                end
                check_result(row)
            finally
                A.finalize()
            end
        end
    end
end
end

module GeoCurveBirthProjectionV1
using Tessella,Test
const A=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","geo_curve_pipeline","fixtures.jl"))
const ARC_ORDER=only(row for row in SECONDORDER_GEO_ROWS if
    row["recipe"]=="arc" && row["linear"]==0 && row["operation"]=="order")
const INITIAL_GEO=replace(ARC_ORDER["geo"],"SetOrder 2;"=>"")

@testset "Native Arc birth UIDs bind separate captured legacy geometry" begin
    A.initialize()
    try
        mktemp() do path,io
            write(io,INITIAL_GEO);close(io);A.open_geo!(path)
        end
        authority=A._model_locked().meshing.mesh_authority
        carrier=A.mesh.get()
        @test carrier isa Tessella.Mesh
        @test authority.mesh===A.LAST_MESH[]
        @test (A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==(5,6)
        @test authority.mesh.entity_data.external_node_tags[
            authority.producer_rows[(1,Int32(11))]]==UInt64[1,4,5,3]
        expected=only(entity for entity in ARC_ORDER["result"]["nodes"] if
            entity["dimension"]==1 && entity["owner"]==11)
        original=[node for node in expected["nodes"] if node["tag"] in (4,5)]
        @test length(original)==2
        geometry_slots=Dict(tag=>slot for (slot,tag) in enumerate(
            authority.mesh.entity_data.external_node_tags[authority.carrier_rows]))
        # Exact pre-existing carrier values are preserved in the failed V8/V9
        # baseline, before any new native birth is bound to these same UIDs.
        legacy=((0.8660254037844387,0.49999999999999994,0.0),
                (0.5000000000000001,0.8660254037844386,0.0))
        for (index,node) in enumerate(original)
            point,parameter,dimension,owner=A.mesh.get_node(node["tag"])
            @test (dimension,owner)==(1,11)
            @test isapprox(point,Float64.(node["coordinates"]);atol=2e-11,rtol=0)
            @test isapprox(parameter,Float64.(node["parameters"]);atol=2e-11,rtol=0)
            @test reinterpret(UInt64,carrier.coords[:,geometry_slots[node["tag"]]])==
                reinterpret(UInt64,collect(legacy[index]))
            @test maximum(abs.(point.-collect(legacy[index])))>2e-11
        end
        @test authority.carrier_cells[1]==UInt64[4,5,6]
    finally
        A.finalize()
    end
    mktemp() do path,io
        write(io,INITIAL_GEO);close(io)
        result=Tessella.execute_geo(path)
        authority=result.model.meshing.mesh_authority
        tags=authority.mesh.entity_data.external_node_tags[authority.carrier_rows]
        slots=Dict(tag=>row for (row,tag) in enumerate(tags))
        @test result.mesh isa Tessella.Mesh
        @test result.mesh.coords[:,slots[4]]==[0.8660254037844387,0.49999999999999994,0.0]
        @test result.mesh.coords[:,slots[5]]==[0.5000000000000001,0.8660254037844386,0.0]
        @test isapprox(result.model.curve_params[11][2:3],
            [0.3333333342552414,0.6666666676645138];atol=2e-11,rtol=0)
    end
end
end
