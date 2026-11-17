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

function owned_actual(;node_uid=UInt64(40),support=Int32[1 1;2 3;3 4],
        cell_uids=UInt64[201,202],point_shift=0.0)
    coordinates=Float64[0 1 1 0 .25;0 0 1 1 .25;0 0 0 0 0]
    coordinates[1,1]+=point_shift
    metadata=E._native_mixed_entity_data(Dict{Tuple{Int,Int},E.MixedEntity}();
        node_entities=fill((2,Int32(7)),5),
        node_parametric=Union{Nothing,Vector{Float64}}[[0.,0.],[1.,0.],[1.,1.],[0.,1.],[.25,.25]],
        external_node_tags=UInt64[10,20,30,node_uid,999],
        block_entities=[Int32[7,7]],external_element_tags=[cell_uids])
    E.MixedMesh(coordinates,[E.ElementBlock(2,support)];entity_data=metadata,
        elementary_entities=[Int32[7,7]])
end
@testset "Surface first-order transition retains only the verified producer receipt" begin
    original=owned_actual(cell_uids=UInt64[101,102])
    old=M._model_mesh_authority(original,[(2,7)])
    old=M._model_mesh_with_producer_rows(old,Dict((2,Int32(7))=>Int32[1,2,3,4]))
    actual=M._model_mesh_product(old,2,7;caller="transition test",producer_order=true)
    receipt=deepcopy(actual)
    old=M._model_mesh_with_surface_constructor_inputs(old,Dict((2,Int32(7))=>receipt))
    current=M._model_mesh_authority(owned_actual(),[(2,7)])
    mapping=Dict(UInt64(101)=>UInt64(201),UInt64(102)=>UInt64(202))
    before=deepcopy(current.mesh)
    rebound=M._model_mesh_first_order_surface_constructor_inputs(old,current,mapping)
    result=rebound[(2,Int32(7))]
    @test result.entity_data.external_node_tags==UInt64[10,20,30,40]
    @test result.entity_data.external_element_tags==[UInt64[201,202]]
    @test current.mesh.entity_data.external_node_tags==UInt64[10,20,30,40,999]
    @test current.mesh.coords==before.coords
    @test current.mesh.blocks[1].nodes==before.blocks[1].nodes
    @test isempty(current.producer_rows)
    @test old.surface_constructor_inputs[(2,Int32(7))]===receipt
    @test receipt.entity_data.external_element_tags==[UInt64[101,102]]
    @test_throws ArgumentError M._model_mesh_with_surface_constructor_inputs(current,
        Dict((2,Int32(7))=>result))
    for bad in (owned_actual(node_uid=UInt64(41)),
                owned_actual(support=Int32[1 1;2 3;3 5]),
                owned_actual(cell_uids=UInt64[201,203]),
                owned_actual(point_shift=nextfloat(0.0)))
        authority=M._model_mesh_authority(bad,[(2,7)])
        @test_throws ArgumentError M._model_mesh_first_order_surface_constructor_inputs(old,authority,mapping)
    end
    @test_throws ArgumentError M._model_mesh_first_order_surface_constructor_inputs(old,current,
        Dict(UInt64(101)=>UInt64(201)))
    @test_throws ArgumentError M._model_mesh_first_order_surface_constructor_inputs(old,current,
        Dict(UInt64(101)=>UInt64(201),UInt64(102)=>UInt64(201)))
    incomplete=deepcopy(old)
    incomplete.producer_rows[(2,Int32(7))]=Int32[1,2,3]
    @test_throws ArgumentError M._model_mesh_first_order_surface_constructor_inputs(incomplete,current,mapping)
end

@testset "Curve source keeps unused owned rows separate from its producer chain" begin
    model=M.GeoModel();model.curves[1]=(1,2)
    metadata=E._native_mixed_entity_data(Dict{Tuple{Int,Int},E.MixedEntity}();
        node_entities=[(0,Int32(1)),(0,Int32(2)),(1,Int32(1)),(1,Int32(1))],
        node_parametric=Union{Nothing,Vector{Float64}}[nothing,nothing,[.5],[.25]],
        external_node_tags=UInt64[1,2,5,999],block_entities=[Int32[1,1]],
        external_element_tags=[UInt64[1,2]])
    mesh=E.MixedMesh(Float64[0 1 .5 .25;0 0 0 0;0 0 0 0],
        [E.ElementBlock(1,Int32[1 3;3 2])];entity_data=metadata,
        elementary_entities=[Int32[1,1]])
    authority=M._model_mesh_authority(mesh,[(1,1)])
    authority=M._model_mesh_with_producer_rows(authority,
        Dict((1,Int32(1))=>Int32[1,3,2]))
    bound=M._model_mesh_curve_births(authority,model,1,"curve source")
    full=M._model_mesh_curve_births(authority,model,1,"curve source";
        used_only=true,producer_order=false)
    @test !(UInt64(999) in bound.mesh.entity_data.external_node_tags)
    @test UInt64(999) in full.mesh.entity_data.external_node_tags
    @test full.mesh.entity_data.external_node_tags[full.rows]==
          bound.mesh.entity_data.external_node_tags[bound.rows]
    @test_throws ArgumentError M._model_mesh_curve_births(authority,model,1,"curve source";
        producer_order=false)
    points,tags,orphans,orphan_tags=M._extrude_curve_source(model,1,authority,"curve source")
    @test tags==UInt64[1,5,2]
    @test orphans==[(.25,0.,0.)]
    @test orphan_tags==UInt64[999]
end

end
