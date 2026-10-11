module GeoDim01CompletionProbeV2
using Tessella,Test
const A=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","geo_dim01_completion","fixtures.jl"))
function check_entity(entity)
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
    end
    types,tags,connectivity=A.mesh.get_elements(dimension,owner)
    @test types==Int.(entity["types"])
    @test tags==[UInt64.(v) for v in entity["element_tags"]]
    @test connectivity==[UInt64.(v) for v in entity["connectivity"]]
end
@testset "GEO native dim01 factory identities and detached chain publication" begin
    for row in GEO_DIM01_PRODUCER_ROWS
        @testset "$(row["recipe"]) R$(row["renumber"]) $(row["kind"])" begin
            A.initialize()
            try
                mktemp() do path,io
                    write(io,row["geo"]);close(io);A.open_geo!(path)
                end
                expected=row["result"]
                @test Set(A.model.get_entities())==Set((e["dimension"],e["owner"]) for e in expected["entities"])
                for e in expected["entities"];check_entity(e);end
                @test (A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==Tuple(expected["maxima"])
                @test Set(A.model.get_physical_groups())==Set((p["dimension"],p["tag"]) for p in expected["physical"])
                for p in expected["physical"]
                    @test A.model.get_entities_for_physical_group(p["dimension"],p["tag"])==Int.(p["members"])
                    @test A.model.get_physical_name(p["dimension"],p["tag"])==p["name"]
                end
                carrier=A.mesh.get()
                @test carrier isa Tessella.Mesh
                authority=A._model_locked().meshing.mesh_authority
                @test authority.mesh===A.LAST_MESH[]
                @test authority.carrier.coords==carrier.coords
                @test length(authority.carrier_rows)==size(carrier.coords,2)
                # Native public coordinates and the retained geometry carrier are
                # separately bound values of the same emitted vertex identities.
                @test all(A.mesh.get_node(authority.mesh.entity_data.external_node_tags[actual_row])[1]==
                    authority.mesh.coords[:,actual_row] for actual_row in authority.carrier_rows)
                @test allunique(authority.carrier_rows)
                @test length(authority.carrier_cells)==3
                @test length(authority.carrier_cells[1])==size(carrier.segs,2)
                @test isempty(authority.carrier_cells[2]) && isempty(authority.carrier_cells[3])
                for (column,tag) in enumerate(authority.carrier_cells[1])
                    _,nodes,dim,owner=A.mesh.get_element(tag)
                    @test dim==1
                    @test nodes==authority.mesh.entity_data.external_node_tags[authority.carrier_rows[carrier.segs[:,column]]]
                end
            finally
                A.finalize()
            end
        end
    end
end
end
