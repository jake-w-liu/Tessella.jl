module NativeHomologyGeoProbe
using Tessella,Test
include("../../validation/homology_native/geo_public_primary_v1.jl")
const M=Tessella.Model
function execute_row(row)
    source=row["source"]
    return mktemp() do path,stream
        write(stream,source);close(stream)
        Tessella.execute_geo(path)
    end
end
@testset "GEO chains retain sparse actual public identities and completed products" begin
    for row in HOMOLOGY_NATIVE_GEO_ROWS
        @testset "$(row["kind"])" begin
            execution=execute_row(row)
            model=execution.model
            expected=row["immediate"]
            group=only(filter(p->p["dimension"]==0 && p["tag"]==11,expected["physical"]))
            chain_tag=only(group["members"])
            chain=model.discrete[(0,chain_tag)]
            @test reduce(vcat,chain.element_nodes;init=Int32[])==row["expected_chain_node_tags"]
            @test chain.element_tags==only(filter(e->e["dimension"]==0 && e["tag"]==chain_tag,expected["entities"]))["element_tags"][1]
            authority=model.meshing.mesh_authority
            @test authority!==nothing && authority.mesh!==nothing
            actual=authority.mesh;data=actual.entity_data
            @test sort(data.external_node_tags)==sort(UInt64[n["tag"] for n in expected["nodes"]])
            @test UInt64(9) ∉ data.external_node_tags && UInt64(10) in data.external_node_tags
            @test all(tag->tag in data.external_node_tags,reduce(vcat,chain.element_nodes;init=Int32[]))
            parameters=M._model_mesh_runtime_parameters(model,actual)
            positions=Dict(tag=>index for (index,tag) in enumerate(data.external_node_tags))
            for node in expected["nodes"]
                position=positions[UInt64(node["tag"])]
                @test actual.coords[:,position]≈node["coordinates"] atol=2e-11 rtol=0
                @test data.node_entities[position]==(node["dimension"],Int32(node["owner"]))
                @test something(parameters[position],Float64[])≈Float64.(node["parameters"]) atol=2e-11 rtol=0
            end
            for entity in expected["entities"]
                key=(entity["dimension"],entity["tag"])
                @test any(part->(part[1],part[2])==key,execution.actual_mesh_parts)
                part=Tessella.geo_entity_mesh(execution,key...)
                @test part isa Tessella.Elements.MixedMesh
                @test part.entity_data.external_node_tags==entity["node_tags"] ||
                    Set(part.entity_data.external_node_tags)==Set(reduce(vcat,entity["connectivity"];init=Int[]))
                for (msh,tags,nodes) in zip(entity["types"],entity["element_tags"],entity["connectivity"])
                    selected=Int[i for i in eachindex(actual.blocks) if actual.blocks[i].msh==msh]
                    found_tags=UInt64[];found_nodes=UInt64[]
                    for bi in selected,column in axes(actual.blocks[bi].nodes,2)
                        data.block_entities[bi][column]==key[2] || continue
                        append!(found_nodes,data.external_node_tags[actual.blocks[bi].nodes[:,column]])
                        push!(found_tags,data.external_element_tags[bi][column])
                    end
                    @test found_tags==tags
                    @test found_nodes==nodes
                end
            end
            @test sort(M.model_entities(model))==[(e["dimension"],e["tag"]) for e in expected["entities"]]
            @test model.physical[(0,11)]==group["members"]
            @test model.physical_names[(0,11)]==group["name"]
            @test M._model_mesh_product(authority,0,chain_tag).entity_data.external_element_tags[1]==chain.element_tags
        end
    end
end
end
