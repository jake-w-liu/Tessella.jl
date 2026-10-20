module SurfaceConstructorReceiptTests
using Tessella,Test
const M=Tessella.Model
const E=Tessella.Elements
using Tessella.MeshTypes: Mesh
coords=Float64[0 1 1 0;0 0 1 1;0 0 0 0]
uv=Union{Nothing,Vector{Float64}}[[0.,0.],[1.,0.],[1.,1.],[0.,1.]]
tags=UInt64[10,20,30,40];owners=fill((2,Int32(7)),4)
cells=Int32[1 1;2 3;3 4]
metadata=E._native_mixed_entity_data(Dict{Tuple{Int,Int},E.MixedEntity}();
    node_entities=owners,node_parametric=uv,external_node_tags=tags,
    block_entities=[Int32[7,7]],external_element_tags=[UInt64[101,102]])
actual=E.MixedMesh(coords,[E.ElementBlock(2,cells)];entity_data=metadata,
    elementary_entities=[Int32[7,7]])
constructor=Mesh(copy(coords);tris=cells[[1,3,2],:])
production=M._ModelSurfaceProduction(Mesh(copy(coords);tris=copy(cells)),
    tags,owners,uv,falses(4),UInt64(40),UInt64[101,102],UInt64(102),
    Int32[1,2,3,4],constructor)
@testset "operation-owned Surface constructor receipt" begin
    receipt=M._model_surface_constructor_input(production,actual,"test")
    @test receipt.blocks[1].nodes==constructor.tris
    @test reinterpret(UInt64,vec(receipt.coords))==reinterpret(UInt64,vec(actual.coords))
    @test receipt.entity_data.node_entities==actual.entity_data.node_entities
    @test receipt.entity_data.node_parametric==actual.entity_data.node_parametric
    @test receipt.entity_data.external_node_tags==actual.entity_data.external_node_tags
    @test receipt.entity_data.external_element_tags==actual.entity_data.external_element_tags
    authority=M._model_mesh_authority(actual,[(2,7)])
    @test isempty(authority.surface_constructor_inputs)
    bound=M._model_mesh_with_surface_constructor_inputs(authority,Dict((2,Int32(7))=>receipt))
    @test bound.mesh===authority.mesh
    @test only(values(bound.surface_constructor_inputs))===receipt
    @test only(values(M._model_mesh_with_producer_rows(bound,bound.producer_rows).surface_constructor_inputs))===receipt
    @test_throws ArgumentError M._model_mesh_with_surface_constructor_inputs(authority,
        Dict((2,Int32(8))=>receipt))
    for mutation in (:xyz,:uv,:owner,:node_uid,:cell_uid,:cell_owner,:physical,:support)
        broken=deepcopy(receipt)
        if mutation===:xyz;broken.coords[1,1]=nextfloat(broken.coords[1,1])
        elseif mutation===:uv;broken.entity_data.node_parametric[1][1]=nextfloat(0.)
        elseif mutation===:owner;broken.entity_data.node_entities[1]=(1,Int32(7))
        elseif mutation===:node_uid;broken.entity_data.external_node_tags[1]=UInt64(11)
        elseif mutation===:cell_uid;broken.entity_data.external_element_tags[1][2]=UInt64(101)
        elseif mutation===:cell_owner;broken.entity_data.block_entities[1][1]=Int32(8)
        elseif mutation===:physical;broken.blocks[1].tags[1]=Int32(9)
        else;broken.blocks[1].nodes[2,1]=Int32(4)
        end
        @test_throws ArgumentError M._model_mesh_with_surface_constructor_inputs(authority,
            Dict((2,Int32(7))=>broken))
    end
    permutation=[4,3,2,1];inverse=Int32[4,3,2,1]
    remapped=M._ModelSurfaceProduction(production.mesh,
        tags[permutation],owners[permutation],uv[permutation],falses(4),UInt64(40),
        UInt64[101,102],UInt64(102),Int32[1,2,3,4],
        Mesh(coords[:,permutation];tris=inverse[constructor.tris]))
    reordered=M._model_surface_constructor_input(remapped,actual,"test")
    @test reordered.blocks[1].nodes==receipt.blocks[1].nodes
    @test reordered.coords==receipt.coords
    @test actual.blocks[1].nodes==cells
    @test actual.entity_data.external_element_tags==[UInt64[101,102]]
end

end
