module NativeHybridQuadConstructorChecks
using Test,TOML,SHA,Tessella
using Tessella.Elements: MixedMesh,ElementBlock,MixedEntity
using Tessella.MeshTypes: nnodes
const Elements=Tessella.Elements
const Model=Tessella.Model
const _Dim01Cell=Tessella.API._Dim01Cell
const _higher_surface_constructor_publication=Tessella.API._higher_surface_constructor_publication
data=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","native_hybrid_quad_constructor_4152_v2.toml"))
@assert data["primary_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"..","artifacts","native_hybrid_quad_birth_4152_v2.json"))))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
ids=UInt64.(data["node_uids"]);coords=reduce(hcat,hexpoint.(data["xyz_hex"]))
node_owners=Tuple{Int,Int32}[(owner[1],Int32(owner[2])) for owner in data["node_owners"]]
parameters=Union{Nothing,Vector{Float64}}[hexpoint(xs) for xs in data["uv_hex"]]
rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(ids))
entities=Dict{Tuple{Int,Int},MixedEntity}()
for owner in union(Set(node_owners),Set((2,Int32(owner)) for owner in data["owners"]))
    points=coords[:,findall(==(owner),node_owners)]
    isempty(points) && (points=coords)
    lo=ntuple(axis->minimum(points[axis,:]),3);hi=ntuple(axis->maximum(points[axis,:]),3)
    entities[(owner[1],Int(owner[2]))]=MixedEntity(owner[1],Int(owner[2]),(lo...,hi...))
end
blocks=ElementBlock[];labels=Vector{UInt64}[];owners=Vector{Int32}[]
for face in sort!(unique(data["owners"]))
    columns=findall(==(face),data["owners"])
    nodes=reduce(hcat,Int32[rows[UInt64(uid)] for uid in data["quad_nodes"][col]] for col in columns)
    push!(blocks,ElementBlock(3,nodes,fill(Int32(0),length(columns))))
    push!(labels,UInt64.(data["quad_uids"][columns]));push!(owners,fill(Int32(face),length(columns)))
end
metadata=Elements._native_mixed_entity_data(entities;node_entities=node_owners,node_parametric=parameters,
    external_node_tags=ids,block_entities=owners,external_element_tags=labels)
actual=MixedMesh(coords,blocks;entity_data=metadata)
authority=Model._model_mesh_authority(actual)
constructors=Dict((2,Int32(face))=>Model._model_mesh_product(authority,2,face;caller="native fixture")
    for face in sort!(unique(data["owners"])))
authority=Model._model_mesh_with_surface_constructor_inputs(authority,constructors)
@testset "current native hybrid constructor binding and fresh final cell events" begin
    snapshot=deepcopy((actual,constructors))
    input=Model._model_native_hybrid_constructor_input(authority,[6,4,2,5,3,1],ids,coords,"hybrid")
    @test [quad.uid for quad in input.quads]==UInt64.(data["quad_uids"])
    @test [quad.nodes for quad in input.quads]==Tuple.(Tuple(UInt64.(x)) for x in data["quad_nodes"])
    @test isempty(input.triangles)
    @test input.input===authority
    @test [input.quads[index].owner[2] for index in input.pyramid_order]==
        Int32[owner for owner in [6,4,2,5,3,1] for _ in findall(==(owner),data["owners"])]
    @test input.coordinates==coords
    @test actual.coords==snapshot[1].coords
    @test actual.entity_data.node_parametric==snapshot[1].entity_data.node_parametric
    @test isempty(authority.done) # A kept actual product has independent provenance.
    @test_throws ArgumentError Model._model_native_hybrid_constructor_input(authority,[1,1],ids,coords,"duplicate face")
    @test_throws ArgumentError Model._model_native_hybrid_constructor_input(authority,[7],ids,coords,"foreign face")
    @test_throws ArgumentError Model._model_native_hybrid_constructor_input(authority,[1],ids[2:end],coords[:,2:end],"missing UID")
    changed=copy(coords);changed[1,1]=nextfloat(changed[1,1])
    @test_throws ArgumentError Model._model_native_hybrid_constructor_input(authority,[1],ids,changed,"Source ULP")
    naked=Model._model_mesh_authority(actual)
    @test_throws ArgumentError Model._model_native_hybrid_constructor_input(naked,[1],ids,coords,"missing receipt")
    edited=deepcopy(constructors);face=(2,Int32(1));part=edited[face]
    for (block,cell_labels) in zip(part.blocks,part.entity_data.external_element_tags)
        columns=reverse(axes(block.nodes,2));block.nodes[:,:]=block.nodes[:,columns];cell_labels[:]=cell_labels[columns]
    end
    bound=Model._model_mesh_with_surface_constructor_inputs(authority,edited)
    reordered=Model._model_native_hybrid_constructor_input(bound,[1],ids,coords,"owned vector")
    @test [quad.uid for quad in reordered.quads]==reverse(UInt64[data["quad_uids"][i] for i in findall(==(1),data["owners"])])
    @test [quad.uid for quad in reordered.quads]!=sort([quad.uid for quad in reordered.quads])
    before_uids=copy(part.entity_data.external_node_tags)
    cell_events=Dict{UInt64,_Dim01Cell}()
    for (owner,constructor) in constructors,(block,cell_labels) in zip(constructor.blocks,constructor.entity_data.external_element_tags),
            (column,uid) in enumerate(cell_labels)
        cell_events[uid]=_Dim01Cell(3,owner,copy(block.nodes[:,column]),uid+UInt64(100000),nothing,Int32(0))
    end
    node_map=Dict(uid=>uid+UInt64(1000) for uid in ids)
    final_actual=deepcopy(actual);final_data=final_actual.entity_data
    final_data.external_node_tags[:]=[node_map[uid] for uid in ids]
    for cell_labels in final_data.external_element_tags
        cell_labels[:]=[cell_events[uid].tag for uid in cell_labels]
    end
    final_authority=Model._model_mesh_authority(final_actual)
    published=_higher_surface_constructor_publication(final_authority,constructors,cell_events,node_map,"fresh final events")
    @test length(published.surface_constructor_inputs)==6
    for (owner,constructor) in published.surface_constructor_inputs
        @test constructor.entity_data.external_node_tags==[node_map[uid] for uid in constructors[owner].entity_data.external_node_tags]
        @test constructor.entity_data.external_element_tags==[
            UInt64[cell_events[uid].tag for uid in block] for block in constructors[owner].entity_data.external_element_tags]
    end
    @test part.entity_data.external_node_tags==before_uids
    @test all(event.source_tag===nothing for event in values(cell_events))
    @test_throws ArgumentError _higher_surface_constructor_publication(final_authority,constructors,cell_events,
        Dict{UInt64,UInt64}(),"missing final node")
    @test_throws ArgumentError _higher_surface_constructor_publication(final_authority,constructors,
        Dict{UInt64,_Dim01Cell}(),node_map,"missing final cell")
    foreign=copy(cell_events);uid=first(keys(foreign));event=foreign[uid]
    foreign[uid]=_Dim01Cell(event.msh,(2,Int32(99)),event.nodes,event.tag,nothing,event.physical)
    @test_throws ArgumentError _higher_surface_constructor_publication(final_authority,constructors,foreign,node_map,"foreign cell")
    @test actual.entity_data.external_node_tags==snapshot[1].entity_data.external_node_tags
    @test actual.entity_data.external_element_tags==snapshot[1].entity_data.external_element_tags
    # A present native live-zero target is distinct from a missing mapping.
    zero_actual=deepcopy(actual);zero_map=Dict(uid=>uid for uid in ids)
    zero_map[first(ids)]=UInt64(0)
    zero_actual.entity_data.external_node_tags[:]=[zero_map[uid] for uid in ids]
    zero_events=Dict{UInt64,_Dim01Cell}()
    first_cell=first(first(actual.entity_data.external_element_tags))
    for (owner,constructor) in constructors,(block,cell_labels) in zip(constructor.blocks,constructor.entity_data.external_element_tags),
            (column,uid) in enumerate(cell_labels)
        zero_events[uid]=_Dim01Cell(3,owner,copy(block.nodes[:,column]),
            uid==first_cell ? UInt64(0) : uid,nothing,Int32(0))
    end
    for cell_labels in zero_actual.entity_data.external_element_tags
        cell_labels[:]=[zero_events[uid].tag for uid in cell_labels]
    end
    zero_authority=Model._model_mesh_authority(zero_actual)
    zero_published=_higher_surface_constructor_publication(zero_authority,constructors,zero_events,zero_map,"native live zero")
    @test UInt64(0) in zero_published.mesh.entity_data.external_node_tags
    @test any(labels->UInt64(0) in labels,zero_published.mesh.entity_data.external_element_tags)
    @test length(zero_published.surface_constructor_inputs)==6
    @test_throws ArgumentError _higher_surface_constructor_publication(zero_authority,constructors,zero_events,
        Dict(uid=>target for (uid,target) in zero_map if uid!=first(ids)),"missing is not live zero")
end
end
