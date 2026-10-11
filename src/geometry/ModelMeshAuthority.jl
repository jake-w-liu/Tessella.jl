# Actual entity products retain storage rank independently of public labels.
# The authority is replaced with its mesh; geometry and connectivity edits must
# rebuild it before a subsequent generation consumes it.
using ..Elements: _msh_primary_nodes,msh_spec
struct _ModelMeshProduct
    owner::Tuple{Int,Int32}
    rows::Vector{Int32}
    blocks::Vector{ElementBlock}
    element_tags::Vector{Vector{UInt64}}
end

include("ModelMeshEmittedProjection.jl")

struct _ModelMeshAuthority
    mesh::Union{Nothing,MixedMesh}
    products::Dict{Tuple{Int,Int32},_ModelMeshProduct}
    done::Set{Tuple{Int,Int32}}
    private_rank::Vector{UInt64}
    native_inputs::Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}
    carrier::Union{Nothing,Mesh,MixedMesh}
    carrier_rows::Vector{Int32}
    producer_rows::Dict{Tuple{Int,Int32},Vector{Int32}}
    carrier_cells::Vector{Vector{Union{Nothing,UInt64}}}
    surface_constructor_inputs::Dict{Tuple{Int,Int32},MixedMesh}
    emitted_geometry::Union{Nothing,_ModelMeshEmittedReceipt}
end

# Every authority reconstruction invalidates the old emitted constructor receipt.
_ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,carrier_cells,surface_constructor_inputs)=
    _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,carrier_cells,surface_constructor_inputs,nothing)

_ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,carrier_cells::Vector{Vector{Union{Nothing,UInt64}}})=
    _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,carrier_cells,Dict{Tuple{Int,Int32},MixedMesh}())

_ModelMeshAuthority(mesh,products,done,private_rank)=_ModelMeshAuthority(
    mesh,products,done,private_rank,Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}())
_ModelMeshAuthority(mesh,products,done,private_rank,native_inputs)=_ModelMeshAuthority(
    mesh,products,done,private_rank,native_inputs,mesh,
    mesh===nothing ? Int32[] : collect(Int32,1:nnodes(mesh)))
_ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows)=
    _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        Dict{Tuple{Int,Int32},Vector{Int32}}())
_ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,producer_rows)=
    _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,Vector{Union{Nothing,UInt64}}[])

# Old UInt64 carriers used zero for missing provenance. A native live-zero
# binding is explicit in the Union vector, or comes from the actual mesh itself.
_model_mesh_cell_binding_type(bindings)=bindings===nothing ? UInt64 :
    any(part->any(block->eltype(block)===Union{Nothing,UInt64},part),bindings) ?
        Union{Nothing,UInt64} : UInt64
_model_mesh_missing_cell_binding(tag,::Type{UInt64})=iszero(tag)
_model_mesh_missing_cell_binding(tag,::Type{Union{Nothing,UInt64}})=isnothing(tag)
_model_mesh_live_cell_tags(tags)=[Union{Nothing,UInt64}[tag for tag in block] for block in tags]
function _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,
        carrier_rows,producer_rows,carrier_cells::Vector{Vector{UInt64}})
    bindings=carrier===mesh ? _model_mesh_live_cell_tags(carrier_cells) :
        [Union{Nothing,UInt64}[tag==0 ? nothing : tag for tag in block] for block in carrier_cells]
    return _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,
        carrier_rows,producer_rows,bindings)
end

function _model_mesh_carrier_blocks(carrier)
    carrier isa Mesh && return ((1,carrier.segs),(2,carrier.tris),(4,carrier.tets))
    carrier isa MixedMesh && return [(block.msh,block.nodes) for block in carrier.blocks]
    return Tuple{Int,Matrix{Int32}}[]
end

function _model_mesh_carrier_cell_tags(mesh,carrier,bindings,caller)
    if bindings===nothing
        return carrier===mesh ? _model_mesh_live_cell_tags(mesh.entity_data.external_element_tags) :
               Vector{Union{Nothing,UInt64}}[]
    end
    isempty(bindings) && return Vector{Union{Nothing,UInt64}}[] # Older carriers have no cell provenance.
    bindings=[Union{Nothing,UInt64}[tag===nothing ? nothing :
        tag==0 && carrier!==mesh && eltype(block)!==Union{Nothing,UInt64} ? nothing : UInt64(tag) for tag in block] for block in bindings]
    blocks=_model_mesh_carrier_blocks(carrier)
    length(bindings)==length(blocks) || throw(ArgumentError(
        "$caller: geometry carrier cell bindings do not align with its blocks"))
    actual=Dict(tag=>msh_dimension(block.msh) for (block,tags) in
        zip(mesh.blocks,mesh.entity_data.external_element_tags) for tag in tags)
    for ((msh,cells),tags) in zip(blocks,bindings)
        length(tags)==size(cells,2) || throw(ArgumentError(
            "$caller: geometry carrier cell bindings do not align with its columns"))
        all(tag->tag===nothing || get(actual,tag,-1)==msh_dimension(msh),tags) ||
            throw(ArgumentError("$caller: geometry carrier cell binding loses its actual identity"))
    end
    return bindings
end

function _model_mesh_rebound_cell_tags(previous,element_tag_map=nothing)
    return [Union{Nothing,UInt64}[tag===nothing ? nothing :
        element_tag_map===nothing ? tag : get(element_tag_map,tag,nothing)
        for tag in tags] for tags in previous.carrier_cells]
end

function _model_mesh_refined_cell_tags(plan,previous;caller="mesh refinement")
    isempty(previous.carrier_cells) && return Vector{Union{Nothing,UInt64}}[]
    previous.carrier isa Mesh && hasproperty(plan,:refine_cell_tags) || return Vector{Union{Nothing,UInt64}}[]
    catalog,_=_model_mesh_primary_catalog(previous.mesh)
    output=Vector{Union{Nothing,UInt64}}[]
    for (bi,(msh,cells)) in enumerate(_model_mesh_carrier_blocks(previous.carrier))
        children=msh==1 ? 2 : msh==2 ? 4 : 8
        tags=Union{Nothing,UInt64}[];sizehint!(tags,Base.checked_mul(children,size(cells,2)))
        for parent in axes(cells,2)
            source=previous.carrier_cells[bi][parent]
            source===nothing && return Vector{Union{Nothing,UInt64}}[]
            actual=get(catalog,source,nothing)
            actual!==nothing && actual[2]==previous.carrier_rows[cells[:,parent]] ||
                return Vector{Union{Nothing,UInt64}}[] # A proxy/permuted parent needs its own child trace.
            produced=get(plan.refine_cell_tags,source,nothing)
            produced!==nothing && length(produced)==children || throw(ArgumentError(
                "$caller: refined carrier parent lacks its actual child events"))
            append!(tags,produced)
        end
        push!(output,tags)
    end
    return output
end

# API supplies the detached identity allocator after both frontends are loaded.
# The operation receives its model, products, options and historical maxima;
# no live API session is consulted by this adapter.
const _MODEL_HIGHER_GENERATOR=Ref{Union{Nothing,Function}}(nothing)
const _MODEL_NATIVE_RECOMBINE_DIAGNOSTICS=Ref{Union{Nothing,NamedTuple}}(nothing)

function _model_mesh_projected_parameter(owner,parameter)
    parameter===nothing && return nothing
    dimension=owner[1]
    length(parameter)>=dimension || return nothing
    return parameter[1:dimension]
end

# Strict mesh metadata projects parameters onto the current owner dimension.
# Original vertex parameters live in the entity records and follow the exact
# public identity when element-based reclassification changes that owner.
function _model_mesh_runtime_parameters(m,mesh::MixedMesh)
    data=mesh.entity_data
    data===nothing && throw(ArgumentError(
        "mesh parameters: actual mesh requires owned entity metadata"))
    parameters=deepcopy(data.node_parametric)
    positions=Dict(tag=>row for (row,tag) in enumerate(data.external_node_tags))
    for ((dimension,tag),record) in _discrete_mesh_records_model(m)
        for (index,node) in enumerate(record.node_tags)
            haskey(record.node_native_params,node) || continue
            row=get(positions,UInt64(node),0)
            row>0 && data.node_entities[row]==(dimension,Int32(tag)) || continue
            stored=_model_mesh_node_parameters(record,index)
            parameters[row]=stored===nothing ? nothing : copy(stored)
        end
    end
    return parameters
end

function _model_mesh_generate_higher(args...;kwargs...)
    generator=_MODEL_HIGHER_GENERATOR[]
    generator===nothing && throw(ArgumentError("higher mesh identity allocator is unavailable"))
    return generator(args...;kwargs...)
end

function _model_mesh_has_native_sources(m)
    for ((dimension,tag),record) in _discrete_mesh_records_model(m)
        haskey(m.discrete,(dimension,tag)) && continue
        dictionary=_model_entity_dictionary(m,dimension)
        haskey(dictionary,tag) || continue
        (!isempty(record.node_tags) || !isempty(record.element_tags)) && return true
    end
    return false
end

function _model_mesh_authority(mesh::MixedMesh,done=();caller="mesh authority",
                               private_rank=collect(UInt64,1:nnodes(mesh)),
                               carrier=mesh,carrier_rows=collect(Int32,1:nnodes(mesh)),
                               carrier_cells=nothing)
    data=mesh.entity_data
    data===nothing && throw(ArgumentError("$caller: actual mesh needs entity metadata"))
    length(data.node_entities)==nnodes(mesh) &&
        length(data.external_node_tags)==nnodes(mesh) &&
        length(data.block_entities)==length(mesh.blocks) &&
        length(data.external_element_tags)==length(mesh.blocks) ||
        throw(ArgumentError("$caller: actual mesh metadata does not align with storage"))
    allunique(data.external_node_tags) || throw(ArgumentError(
        "$caller: actual node labels must be unique"))
    length(private_rank)==nnodes(mesh) && all(>(0),private_rank) &&
        allunique(private_rank) || throw(ArgumentError(
        "$caller: private vertex ranks must be positive, unique and align with actual rows"))
    carrier===nothing ? isempty(carrier_rows) || throw(ArgumentError(
        "$caller: absent geometry carrier has nonempty row bindings")) :
        length(carrier_rows)==nnodes(carrier) && all(row->0<row<=nnodes(mesh),carrier_rows) &&
        allunique(carrier_rows) || throw(ArgumentError(
            "$caller: geometry carrier rows do not align with actual identities"))
    selections=Dict{Tuple{Int,Int32},Vector{Tuple{Int,Vector{Int}}}}()
    owned_rows=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (row,owner) in enumerate(data.node_entities)
        push!(get!(()->Int32[],owned_rows,owner),Int32(row))
        get!(()->Tuple{Int,Vector{Int}}[],selections,owner)
    end
    all_tags=Set{UInt64}()
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock || throw(ArgumentError("$caller: unsupported mesh block"))
        owners=data.block_entities[bi];tags=data.external_element_tags[bi]
        length(owners)==size(block.nodes,2)==length(tags) || throw(ArgumentError(
            "$caller: actual cell metadata does not align with connectivity"))
        for tag in tags
            !(tag in all_tags) || throw(ArgumentError(
                "$caller: actual element labels must be positive, unique and Int32 representable"))
            push!(all_tags,tag)
        end
        by_owner=Dict{Int32,Vector{Int}}()
        for (column,owner) in enumerate(owners)
            push!(get!(()->Int[],by_owner,owner),column)
        end
        for (owner,columns) in by_owner
            key=(msh_dimension(block.msh),owner)
            push!(get!(()->Tuple{Int,Vector{Int}}[],selections,key),(bi,columns))
        end
    end
    products=Dict{Tuple{Int,Int32},_ModelMeshProduct}()
    for (owner,selected) in selections
        referenced=Set(get(owned_rows,owner,Int32[]))
        for (bi,columns) in selected,node in mesh.blocks[bi].nodes[:,columns]
            0<node<=nnodes(mesh) || throw(ArgumentError(
                "$caller: entity $owner references an invalid actual node row"))
            push!(referenced,node)
        end
        # Each product touches its own rows/connectivity. Allocating a global
        # mask and map per entity makes sparse many-entity imports quadratic.
        rows=sort!(collect(referenced);by=row->private_rank[row])
        mapping=Dict(row=>Int32(index) for (index,row) in enumerate(rows))
        blocks=ElementBlock[];tags=Vector{UInt64}[]
        for (bi,columns) in selected
            block=mesh.blocks[bi]
            local_nodes=map(row->mapping[row],block.nodes[:,columns])
            push!(blocks,ElementBlock(block.msh,local_nodes,block.tags[columns]))
            push!(tags,copy(data.external_element_tags[bi][columns]))
        end
        products[owner]=_ModelMeshProduct(owner,rows,blocks,tags)
    end
    completed=Set{Tuple{Int,Int32}}((Int(dim),Int32(tag)) for (dim,tag) in done)
    bindings=_model_mesh_carrier_cell_tags(mesh,carrier,carrier_cells,caller)
    return _ModelMeshAuthority(mesh,products,completed,UInt64.(private_rank),
        Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}(),carrier,Int32.(carrier_rows),
        Dict{Tuple{Int,Int32},Vector{Int32}}(),bindings)
end

# Product metadata retains the declarations it uses and their boundary,
# embedding and parent closure. Copying every declaration for every sparse
# product otherwise dominates import allocation quadratically.
function _model_mesh_product_entities(source,owner,rows)
    data=source.entity_data
    selected=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    queue=Tuple{Int,Int}[(owner[1],Int(owner[2]))]
    append!(queue,[(key[1],Int(key[2])) for key in data.node_entities[rows]])
    visited=Set{Tuple{Int,Int}}();cursor=1
    while cursor<=length(queue)
        key=queue[cursor];cursor+=1
        key in visited && continue
        push!(visited,key)
        entity=get(data.entities,key,nothing)
        if entity===nothing && source.partition_data!==nothing
            entity=get(source.partition_data.entities,key,nothing)
        end
        entity===nothing && continue
        selected[key]=entity
        for tag in entity.boundaries
            push!(queue,(entity.dim-1,abs(Int(tag))))
        end
        for tag in entity.embedded_curves
            push!(queue,(1,Int(tag)))
        end
        entity.parent[1]<0 || push!(queue,(entity.parent[1],Int(entity.parent[2])))
    end
    return selected
end

function _model_mesh_product(authority::_ModelMeshAuthority,dim::Int,tag::Int;
                             caller="mesh authority",producer_order=false)
    product=get(authority.products,(dim,Int32(tag)),nothing)
    product===nothing && return nothing
    source=authority.mesh;data=source.entity_data;rows=product.rows
    blocks=product.blocks
    produced=producer_order ? get(authority.producer_rows,product.owner,nothing) : nothing
    if produced!==nothing
        # Normal PLC input follows the actual face producer vector. NoNew
        # consumers keep the separate private-rank order above.
        inverse=zeros(Int32,nnodes(source));inverse[produced]=Int32.(1:length(produced))
        # Product rows also list owned vertices with no cell incidence (an
        # AddNodes orphan or a retained predecessor row). They have no
        # connectivity to reorder and drop out through rows=produced below;
        # only a used row missing its producer binding invalidates the source.
        for block in blocks,node in block.nodes
            inverse[rows[node]]>0 || throw(ArgumentError(
                "$caller: actual entity $(product.owner) changed without a complete producer row binding"))
        end
        mapping=inverse[rows]
        blocks=ElementBlock[ElementBlock(block.msh,mapping[block.nodes],block.tags)
                            for block in product.blocks]
        rows=produced
    end
    owners=[fill(Int32(tag),size(block.nodes,2)) for block in blocks]
    entities=_model_mesh_product_entities(source,product.owner,rows)
    metadata=Elements._native_mixed_entity_data(entities;
        node_entities=data.node_entities[rows],node_parametric=data.node_parametric[rows],
        external_node_tags=data.external_node_tags[rows],block_entities=owners,
        external_element_tags=product.element_tags)
    return MixedMesh(source.coords[:,rows],blocks;
        physical_names=source.physical_names,entity_data=metadata,elementary_entities=owners)
end

function _model_surface_constructor_binding(actual::MixedMesh,constructor::MixedMesh,
                                             caller::AbstractString)
    actual.entity_data!==nothing && constructor.entity_data!==nothing ||
        throw(ArgumentError("$caller: Surface constructor input lacks actual identity metadata"))
    data=actual.entity_data;receipt=constructor.entity_data
    length(data.node_entities)==length(data.node_parametric)==nnodes(actual) &&
        length(receipt.node_entities)==length(receipt.node_parametric)==nnodes(constructor) ||
        throw(ArgumentError("$caller: Surface constructor node metadata does not align"))
    actual_rows=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    if !(length(actual_rows)==nnodes(actual)==nnodes(constructor) &&
        allunique(receipt.external_node_tags) &&
        all(uid->haskey(actual_rows,uid),receipt.external_node_tags))
        throw(ArgumentError(
            "$caller: Surface constructor input changes the actual node UID set"))
    end
    for (row,uid) in enumerate(receipt.external_node_tags)
        actual_row=actual_rows[uid]
        all(axis->isequal(constructor.coords[axis,row],actual.coords[axis,actual_row]),1:3) &&
            receipt.node_entities[row]==data.node_entities[actual_row] &&
            isequal(receipt.node_parametric[row],data.node_parametric[actual_row]) ||
            throw(ArgumentError("$caller: Surface constructor input changes bound node $uid"))
    end
    cells=Dict{UInt64,Tuple{Int,Int32,Int32,Vector{UInt64}}}()
    length(actual.blocks)==length(data.external_element_tags)==length(data.block_entities) ||
        throw(ArgumentError("$caller: actual Surface cell metadata does not align"))
    for (block,labels,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)
        block.msh in (2,3) && length(labels)==length(owners)==size(block.nodes,2) || throw(ArgumentError(
            "$caller: Surface constructor input requires actual linear Surface cells"))
        for (column,label) in enumerate(labels)
            haskey(cells,label) && throw(ArgumentError(
                "$caller: actual Surface cells repeat constructor UID $label"))
            support=sort!(data.external_node_tags[block.nodes[:,column]])
            cells[label]=(Int(block.msh),owners[column],block.tags[column],support)
        end
    end
    seen=Set{UInt64}()
    length(constructor.blocks)==length(receipt.external_element_tags)==length(receipt.block_entities) ||
        throw(ArgumentError("$caller: Surface constructor cell labels do not align"))
    for (block,labels,owners) in zip(constructor.blocks,receipt.external_element_tags,receipt.block_entities)
        length(labels)==length(owners)==size(block.nodes,2) || throw(ArgumentError(
            "$caller: Surface constructor cell labels do not align"))
        for (column,label) in enumerate(labels)
            label in seen && throw(ArgumentError(
                "$caller: Surface constructor input repeats cell UID $label"))
            push!(seen,label)
            support=sort!(receipt.external_node_tags[block.nodes[:,column]])
            get(cells,label,nothing)==(Int(block.msh),owners[column],block.tags[column],support) || throw(ArgumentError(
                "$caller: Surface constructor input changes cell $label owner or support"))
        end
    end
    length(seen)==length(cells) || throw(ArgumentError(
        "$caller: Surface constructor input omits actual cells"))
    return nothing
end

function _model_surface_constructor_input(production,actual::MixedMesh,
                                           caller::AbstractString;
                                           constructor=production.constructor_mesh)
    constructor===nothing && return nothing
    data=actual.entity_data
    data!==nothing || throw(ArgumentError(
        "$caller: Surface constructor input lacks the authoritative product"))
    length(production.node_tags)==nnodes(constructor) &&
        allunique(production.node_tags) || throw(ArgumentError(
        "$caller: Surface constructor rows do not bind unique actual UIDs"))
    actual_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    all(uid->haskey(actual_rows,uid),production.node_tags) || throw(ArgumentError(
        "$caller: Surface constructor row UID is absent from its actual product"))
    mapping=Int32[actual_rows[uid] for uid in production.node_tags]
    for row in eachindex(mapping)
        bound=mapping[row]
        all(axis->isequal(constructor.coords[axis,row],actual.coords[axis,bound]),1:3) &&
            production.node_owners[row]==data.node_entities[bound] &&
            isequal(production.node_parameters[row],data.node_parametric[bound]) ||
            throw(ArgumentError("$caller: Surface constructor row loses its actual node binding"))
    end
    physical_tags=Dict(uid=>block.tags[column] for (block,uids) in
        zip(actual.blocks,data.external_element_tags) for (column,uid) in enumerate(uids))
    blocks=ElementBlock[];labels=Vector{UInt64}[];cursor=1
    for (msh,nodes) in _surface_dim2_blocks(constructor)
        count=size(nodes,2)
        cursor+count-1<=length(production.element_tags) || throw(ArgumentError(
            "$caller: Surface constructor cells lack actual birth labels"))
        cell_labels=copy(production.element_tags[cursor:cursor+count-1])
        all(uid->haskey(physical_tags,uid),cell_labels) || throw(ArgumentError(
            "$caller: Surface constructor cell UID is absent from its actual product"))
        push!(blocks,ElementBlock(msh,mapping[nodes],Int32[physical_tags[uid] for uid in cell_labels]))
        push!(labels,cell_labels)
        cursor+=count
    end
    cursor==length(production.element_tags)+1 || throw(ArgumentError(
        "$caller: Surface constructor cells do not cover their birth labels"))
    metadata=Elements._native_mixed_entity_data(data.entities;
        node_entities=copy(data.node_entities),node_parametric=deepcopy(data.node_parametric),
        external_node_tags=copy(data.external_node_tags),
        block_entities=deepcopy(data.block_entities),external_element_tags=labels)
    result=MixedMesh(copy(actual.coords),blocks;physical_names=copy(actual.physical_names),
        entity_data=metadata,elementary_entities=deepcopy(actual.elementary_entities))
    _model_surface_constructor_binding(actual,result,caller)
    return result
end

function _model_mesh_with_surface_constructor_inputs(authority::_ModelMeshAuthority,inputs;
        caller="Surface constructor input")
    bound=Dict{Tuple{Int,Int32},MixedMesh}()
    for (owner,constructor) in inputs
        owner[1]==2 && haskey(authority.products,owner) || throw(ArgumentError(
            "$caller: Surface constructor input refers to absent actual entity $owner"))
        # Constructor inputs are bound against the producer-order product at
        # creation time; re-validate against the same producer view so unused
        # owned rows cannot change the actual node set comparison.
        actual=_model_mesh_product(authority,2,Int(owner[2]);caller,producer_order=true)
        _model_surface_constructor_binding(actual,constructor,caller)
        bound[owner]=constructor
    end
    return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
        authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,
        authority.producer_rows,authority.carrier_cells,bound)
end

function _model_mesh_with_empty_surface_product(authority::_ModelMeshAuthority,owner,rows;
                                                caller="mesh authority")
    owner[1]==2 && owner in authority.done &&
        haskey(authority.mesh.entity_data.entities,(2,Int(owner[2]))) || throw(ArgumentError(
            "$caller: empty Surface producer lacks its completed actual declaration"))
    haskey(authority.products,owner) && throw(ArgumentError(
        "$caller: empty Surface producer cannot replace an existing actual product"))
    all(row->0<row<=nnodes(authority.mesh),rows) && allunique(rows) &&
        all(row->authority.mesh.entity_data.node_entities[row][1]<2,rows) ||
        throw(ArgumentError("$caller: empty Surface producer must bind only retained lower vertices"))
    products=copy(authority.products)
    products[owner]=_ModelMeshProduct(owner,Int32.(rows),ElementBlock[],Vector{UInt64}[])
    return _ModelMeshAuthority(authority.mesh,products,authority.done,authority.private_rank,
        authority.native_inputs,authority.carrier,authority.carrier_rows,authority.producer_rows,
        authority.carrier_cells,authority.surface_constructor_inputs)
end

function _model_mesh_with_producer_rows(authority::_ModelMeshAuthority,producer_rows;
                                        caller="mesh authority")
    rows=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (owner,positions) in producer_rows
        haskey(authority.products,owner) || throw(ArgumentError(
            "$caller: producer row binding refers to absent actual entity $owner"))
        all(row->0<row<=nnodes(authority.mesh),positions) && allunique(positions) ||
            throw(ArgumentError("$caller: producer rows for $owner do not bind unique actual vertices"))
        rows[owner]=Int32.(positions)
    end
    return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
        authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,rows,
        authority.carrier_cells)
end

function _model_mesh_with_curve_births(authority::_ModelMeshAuthority,births;
                                      node_tag_map=nothing,caller="mesh authority")
    produced=copy(authority.producer_rows)
    positions=Dict(tag=>Int32(row) for (row,tag) in
                   enumerate(authority.mesh.entity_data.external_node_tags))
    for (owner,tags) in births
        owner[1]==1 && length(tags)>=2 || throw(ArgumentError(
            "$caller: native Curve birth vector must contain its two endpoints"))
        haskey(authority.products,owner) || throw(ArgumentError(
            "$caller: Curve birth vector refers to an absent actual product"))
        rows=Int32[get(positions,node_tag_map===nothing ? tag :
            get(node_tag_map,tag,nothing),Int32(0)) for tag in tags]
        # A closed emitter repeats its initial endpoint as the last incidence.
        # The authority stores each actual vertex once; the Line cells retain
        # the closing incidence, independently of this producer row vector.
        if length(rows)>1 && rows[begin]==rows[end]
            rows=rows[begin:end-1]
        end
        all(>(0),rows) && allunique(rows) || throw(ArgumentError(
            "$caller: native Curve $owner birth vector loses unique actual identities ($tags -> $rows)"))
        produced[owner]=rows
    end
    return _model_mesh_with_producer_rows(authority,produced;caller)
end

function _model_mesh_primary_key(family,tags)
    return (family,Tuple(sort!(collect(tags))))
end

function _model_mesh_primary_catalog(mesh::MixedMesh;include_supports=false)
    data=mesh.entity_data
    catalog=Dict{UInt64,Tuple{Symbol,Vector{Int32}}}()
    identities=Dict{Any,Union{Nothing,UInt64}}()
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock || continue
        family=msh_family(block.msh)
        width=_msh_primary_nodes(family)
        width==0 && continue
        for column in axes(block.nodes,2)
            nodes=include_supports ? block.nodes[:,column] : block.nodes[1:width,column]
            tag=data.external_element_tags[bi][column]
            catalog[tag]=(family,nodes)
            key=_model_mesh_primary_key(family,data.external_node_tags[block.nodes[1:width,column]])
            identities[key]=haskey(identities,key) ? nothing : tag
        end
    end
    return catalog,identities
end

function _model_mesh_rebound_carrier(previous::_ModelMeshAuthority,mesh::MixedMesh;
        node_tag_map=nothing,element_tag_map=nothing)
    carrier=previous.carrier
    carrier!==nothing && previous.mesh!==nothing || return nothing,Int32[]
    carrier isa MixedMesh && isempty(previous.carrier_cells) && return nothing,Int32[]
    old_tags=previous.mesh.entity_data.external_node_tags
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(mesh.entity_data.external_node_tags))
    rows=Int32[get(positions,node_tag_map===nothing ? old_tags[row] :
        get(node_tag_map,old_tags[row],nothing),Int32(0)) for row in previous.carrier_rows]
    all(>(0),rows) && allunique(rows) || return nothing,Int32[]
    inverse=zeros(Int32,nnodes(mesh));inverse[rows]=Int32.(1:length(rows))
    _,old_cells=_model_mesh_primary_catalog(previous.mesh)
    include_supports=carrier isa MixedMesh && any(block->msh_spec(block.msh).order>1,carrier.blocks)
    new_cells,_=_model_mesh_primary_catalog(mesh;include_supports)
    output=Matrix{Int32}[]
    for (bi,(msh,cells)) in enumerate(_model_mesh_carrier_blocks(carrier))
        family=msh_family(msh)
        rebound=similar(cells)
        for column in axes(cells,2)
            old_tag=if isempty(previous.carrier_cells)
                key=_model_mesh_primary_key(family,old_tags[previous.carrier_rows[cells[:,column]]])
                get(old_cells,key,nothing)
            else
                previous.carrier_cells[bi][column]
            end
            old_tag!==nothing || return nothing,Int32[]
            tag=element_tag_map===nothing ? old_tag : get(element_tag_map,old_tag,nothing)
            current=get(new_cells,tag,nothing)
            current!==nothing && current[1]===family || return nothing,Int32[]
            nodes=current[2]
            length(nodes)==size(cells,1) || return nothing,Int32[]
            all(node->inverse[node]>0,nodes) || return nothing,Int32[]
            rebound[:,column]=inverse[nodes]
        end
        push!(output,rebound)
    end
    rebound=if carrier isa Mesh
        coords=copy(mesh.coords[:,rows])
        for (slot,(old_row,new_row)) in enumerate(zip(previous.carrier_rows,rows))
            unchanged=previous.mesh.entity_data.node_entities[old_row]==mesh.entity_data.node_entities[new_row] &&
                isequal(previous.mesh.coords[:,old_row],mesh.coords[:,new_row]) &&
                isequal(previous.mesh.entity_data.node_parametric[old_row],
                        mesh.entity_data.node_parametric[new_row])
            unchanged && (coords[:,slot]=carrier.coords[:,slot])
        end
        Mesh(coords;segs=output[1],tris=output[2],tets=output[3],
             seg_tag=carrier.seg_tag,tri_tag=carrier.tri_tag,tet_tag=carrier.tet_tag)
    else
        blocks=ElementBlock[ElementBlock(block.msh,nodes,block.tags) for
            (block,nodes) in zip(carrier.blocks,output)]
        coords=copy(mesh.coords[:,rows])
        old_data=carrier.entity_data
        data=nothing
        if old_data!==nothing
            parameters=deepcopy(mesh.entity_data.node_parametric[rows])
            for (slot,(old_row,new_row)) in enumerate(zip(previous.carrier_rows,rows))
                unchanged=previous.mesh.entity_data.node_entities[old_row]==mesh.entity_data.node_entities[new_row] &&
                    isequal(previous.mesh.coords[:,old_row],mesh.coords[:,new_row]) &&
                    isequal(previous.mesh.entity_data.node_parametric[old_row],
                            mesh.entity_data.node_parametric[new_row])
                if unchanged
                    coords[:,slot]=carrier.coords[:,slot]
                    parameters[slot]=old_data.node_parametric[slot]
                end
            end
            data=Elements._native_mixed_entity_data(mesh.entity_data.entities;
                node_entities=mesh.entity_data.node_entities[rows],node_parametric=parameters,
                external_node_tags=mesh.entity_data.external_node_tags[rows],
                block_entities=old_data.block_entities,
                external_element_tags=[UInt64[tag for tag in block] for block in
                    _model_mesh_rebound_cell_tags(previous,element_tag_map)])
        end
        MixedMesh(coords,blocks;physical_names=carrier.physical_names,entity_data=data,
            elementary_entities=carrier.elementary_entities,periodic_links=carrier.periodic_links,
            ancillary_sections=carrier.ancillary_sections,data_sections=carrier.data_sections,
            partition_data=carrier.partition_data)
    end
    return rebound,rows
end

function _model_mesh_refined_carrier(plan,carrier::Mesh,tags;caller="mesh refinement")
    trace=Refine._refine_uniform_traced(carrier;require_positive_tets=false)
    length(tags)==nnodes(carrier) || throw(ArgumentError(
        "$caller: original carrier identities do not align with refinement input"))
    rows=Vector{Int32}(undef,nnodes(trace.mesh))
    old_rows=zeros(Int32,count(!iszero,trace.primary_map))
    for (old,compact) in enumerate(trace.primary_map)
        compact==0 && continue
        old_rows[compact]=Int32(old)
        rows[compact]=plan.refine_primary_rows[tags[old]]
    end
    for ((a,b),midpoint) in trace.midpoint_ids
        key=minmax(tags[old_rows[a]],tags[old_rows[b]])
        rows[midpoint]=plan.refine_edge_rows[key]
    end
    allunique(rows) && all(row->0<row<=nnodes(plan.mesh),rows) || throw(ArgumentError(
        "$caller: refined carrier identities do not form a valid injection"))
    # Node provenance alone does not certify a carrier's new cells.
    catalog,actual=_model_mesh_primary_catalog(plan.mesh)
    labels=plan.mesh.entity_data.external_node_tags
    for (family,cells) in ((:lin,trace.mesh.segs),(:tri,trace.mesh.tris),(:tet,trace.mesh.tets))
        for column in axes(cells,2)
            key=_model_mesh_primary_key(family,labels[rows[cells[:,column]]])
            tag=get(actual,key,nothing)
            tag!==nothing || return nothing,Int32[]
            source=rows[cells[:,column]];target=catalog[tag][2]
            permutation=Int[findfirst(==(node),target) for node in source]
            inversions=sum(permutation[i]>permutation[j] for i in eachindex(permutation)
                for j in i+1:length(permutation);init=0)
            iseven(inversions) || return nothing,Int32[]
        end
    end
    refined=trace.mesh
    output=Mesh(plan.mesh.coords[:,rows];segs=refined.segs,tris=refined.tris,tets=refined.tets,
        seg_tag=refined.seg_tag,tri_tag=refined.tri_tag,tet_tag=refined.tet_tag)
    return output,rows
end

function _model_mesh_completed(m,dimension::Int;native_inputs=nothing)
    done=Set{Tuple{Int,Int32}}()
    for curve in keys(m.curve_params)
        push!(done,(1,Int32(curve)))
    end
    if dimension>=2
        for surface in keys(m.surfaces)
            push!(done,(2,Int32(surface)))
        end
    end
    inputs=Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}()
    if native_inputs!==nothing
        for (tag,mesh) in native_inputs
            inputs[(2,Int32(tag))]=mesh
        end
    end
    return _ModelMeshAuthority(nothing,Dict{Tuple{Int,Int32},_ModelMeshProduct}(),done,UInt64[],inputs)
end

function _model_mesh_rebind_authority(previous::_ModelMeshAuthority,mesh::MixedMesh)
    previous.mesh===nothing && return previous
    old=previous.mesh.entity_data
    by_tag=Dict(tag=>previous.private_rank[row] for (row,tag) in enumerate(old.external_node_tags))
    rank=Vector{UInt64}(undef,nnodes(mesh));next=maximum(previous.private_rank;init=UInt64(0))
    for (row,tag) in enumerate(mesh.entity_data.external_node_tags)
        value=get(by_tag,tag,UInt64(0))
        if value==0
            next<typemax(UInt64) || throw(ArgumentError("mesh authority: private vertex rank overflows UInt64"))
            next+=UInt64(1);value=next
        end
        rank[row]=value
    end
    carrier,rows=_model_mesh_rebound_carrier(previous,mesh)
    authority=carrier===nothing ? _model_mesh_authority(mesh,previous.done;private_rank=rank) :
        _model_mesh_authority(mesh,previous.done;private_rank=rank,carrier,carrier_rows=rows,
            carrier_cells=_model_mesh_rebound_cell_tags(previous))
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(mesh.entity_data.external_node_tags))
    produced=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (owner,old_rows) in previous.producer_rows
        haskey(authority.products,owner) || continue
        rebound=Int32[get(positions,old.external_node_tags[row],Int32(0)) for row in old_rows]
        all(>(0),rebound) || continue
        produced[owner]=rebound
    end
    return _model_mesh_with_producer_rows(authority,produced)
end

function _model_mesh_status(m,authority::_ModelMeshAuthority;only_visible=false)
    authority.mesh===nothing && return -1
    have=zeros(Int,4)
    for (owner,product) in authority.products
        entity=(owner[1],Int(owner[2]))
        haskey(m.discrete,entity) && continue
        have[owner[1]+1]+=sum(size(block.nodes,2) for block in product.blocks;init=0)
    end
    have[4]>0 && return 3
    for dim in (2,1)
        have[dim+1]>0 || continue
        entities=(tag for (edim,tag) in model_entities(m) if edim==dim)
        complete=true
        for tag in entities
            only_visible && model_entity_visibility(m,dim,tag)==0 && continue
            if !((dim,Int32(tag)) in authority.done)
                complete=false
                break
            end
        end
        complete && return dim
    end
    any(owner->owner[1]==0 && haskey(m.points,Int(owner[2])),
        authority.mesh.entity_data.node_entities) && return 0
    return -1
end

@inline function _model_mesh_generation_phases(dimension::Int,status::Int)
    dimension in (2,3) || throw(ArgumentError("higher generation requires dimension two or three"))
    status in -1:3 || throw(ArgumentError("invalid actual mesh status"))
    return (one=status<1,two=dimension==2 || status<2,three=dimension==3)
end

function _model_mesh_curve_births(authority,model,curve::Int,caller;used_only=false)
    authority===nothing && return nothing
    mesh=_model_mesh_product(authority,1,curve;caller)
    mesh===nothing && return nothing
    all(block->isempty(block.nodes),mesh.blocks) && return nothing
    all(block->block.msh==1,mesh.blocks) || throw(ArgumentError(
        "$caller: actual Curve[$curve] input must contain linear Line cells"))
    endpoints=get(model.curves,curve,nothing)
    endpoints===nothing && throw(ArgumentError("$caller: actual Curve[$curve] is not native"))
    endpoints[1]!=endpoints[2] || throw(ArgumentError(
        "$caller: actual Curve[$curve] chain must be open"))
    data=mesh.entity_data
    active=used_only ? falses(nnodes(mesh)) : trues(nnodes(mesh))
    if used_only
        for block in mesh.blocks,row in block.nodes
            active[row]=true
        end
    end
    function endpoint_rows(point_tag)
        # Mesh0D chooses its first owned vertex for the Curve endpoint, while
        # an absent Point cell is constructed from the last owned vertex.
        # Bind the Point-owned row already incident to this actual Line chain
        # before consulting Point-cell incidence retained after order pruning.
        owned=findall(row->active[row] && data.node_entities[row]==(0,Int32(point_tag)),eachindex(active))
        isempty(owned) || return owned
        product=_model_mesh_product(authority,0,point_tag;caller)
        point_tags=Set{UInt64}()
        if product!==nothing
            for block in product.blocks,column in axes(block.nodes,2)
                block.msh==15 || continue
                push!(point_tags,product.entity_data.external_node_tags[block.nodes[1,column]])
            end
        end
        return findall(row->active[row] && data.external_node_tags[row] in point_tags,eachindex(active))
    end
    born=get(authority.producer_rows,(1,Int32(curve)),nothing)
    starts,ends=if born===nothing
        endpoint_rows(endpoints[1]),endpoint_rows(endpoints[2])
    else
        # These are actual Curve emitter endpoints, retained through explicit
        # node-tag maps. A Point cell can use a different coincident vertex.
        tags=authority.mesh.entity_data.external_node_tags
        findall(row->active[row] && data.external_node_tags[row]==tags[born[begin]],eachindex(active)),
            findall(row->active[row] && data.external_node_tags[row]==tags[born[end]],eachindex(active))
    end
    length(starts)==length(ends)==1 || throw(ArgumentError(
        "$caller: actual Curve[$curve] chain has ambiguous endpoint identities"))
    neighbors=[Int32[] for _ in 1:nnodes(mesh)]
    for block in mesh.blocks,column in axes(block.nodes,2)
        a,b=block.nodes[:,column]
        a!=b || throw(ArgumentError("$caller: actual Curve[$curve] repeats a Line node"))
        push!(neighbors[a],b);push!(neighbors[b],a)
    end
    first=Int32(only(starts));last=Int32(only(ends))
    all(row->!active[row] || length(neighbors[row])==(row in (first,last) ? 1 : 2),eachindex(neighbors)) ||
        throw(ArgumentError("$caller: actual Curve[$curve] cells do not form one open chain"))
    rows=Int32[first];visited=falses(nnodes(mesh));visited[first]=true
    previous=Int32(0);current=first
    while current!=last
        adjacent=neighbors[current]
        next=adjacent[1]==previous ? adjacent[2] : adjacent[1]
        visited[next] && throw(ArgumentError("$caller: actual Curve[$curve] contains a cycle"))
        push!(rows,next);visited[next]=true;previous=current;current=next
    end
    length(rows)==count(active) || throw(ArgumentError(
        "$caller: actual Curve[$curve] contains disconnected Line cells"))
    return (mesh=mesh,rows=rows)
end

function _model_mesh_curve_input(authority,model,curve::Int,caller)
    source=_model_mesh_curve_births(authority,model,curve,caller)
    source===nothing && return nothing
    return NTuple{3,Float64}[ntuple(axis->source.mesh.coords[axis,row],3) for row in source.rows]
end

function _model_mesh_top_states(top_states,cells::Int,caller)
    top_states===nothing && return nothing
    top_states isa AbstractVector && length(top_states)==cells &&
        all(state->state isa Integer && !(state isa Bool) && state in (1,2),top_states) ||
        throw(ArgumentError("$caller: retained actual cap states do not align with source Quads"))
    Base.require_one_based_indexing(top_states)
    return top_states
end

struct _ModelMeshCapInput
    top_states::Vector{UInt8}
    mesh::MixedMesh
    upper_rows::Vector{Int32}
end

function _model_mesh_cap_input(authority,source_mesh,top::Int,columns,caller)
    authority===nothing && return nothing
    cap=_model_mesh_product(authority,2,top;caller)
    cap===nothing && return nothing
    all(block->isempty(block.nodes),cap.blocks) && return nothing
    source_mesh isa MixedMesh && length(source_mesh.blocks)==1 &&
        only(source_mesh.blocks).msh==3 || throw(ArgumentError(
        "$caller: retained cap requires the actual linear Quad source"))
    length(cap.blocks)==1 && only(cap.blocks).msh==2 || throw(ArgumentError(
        "$caller: retained cap requires one ordered linear Triangle product"))
    source=only(source_mesh.blocks).nodes;triangles=only(cap.blocks).nodes
    size(triangles,2)==Base.checked_mul(2,size(source,2)) || throw(ArgumentError(
        "$caller: retained cap must contain two actual triangles per source Quad"))
    size(columns,2)==nnodes(source_mesh) || throw(ArgumentError(
        "$caller: retained cap source columns do not align with actual source rows"))
    targets=Dict{NTuple{3,Float64},Int32}()
    for row in axes(columns,2)
        point=columns[end,row]
        haskey(targets,point) && throw(ArgumentError(
            "$caller: actual source top columns are coincident"))
        targets[point]=Int32(row)
    end
    spans=ntuple(axis->begin
        bounds=extrema(columns[end,row][axis] for row in axes(columns,2))
        bounds[2]-bounds[1]
    end,3)
    search_axis=argmax(spans)
    column_order=sortperm(collect(axes(columns,2));
        by=row->columns[end,row][search_axis],alg=Base.Sort.MergeSort)
    ordered_values=Float64[columns[end,row][search_axis] for row in column_order]
    projection=zeros(Int32,nnodes(cap));upper=zeros(Int32,nnodes(source_mesh))
    for row in axes(cap.coords,2)
        point=ntuple(axis->cap.coords[axis,row],3)
        target=get(targets,point,Int32(0))
        if target==0
            # This verifies a known analytical column, without deriving the
            # cap vector or a source cell order from a volume's faces.
            tolerance=64eps(Float64)*max(1.0,maximum(abs,point))
            lower=searchsortedfirst(ordered_values,point[search_axis]-tolerance)
            upper_bound=searchsortedlast(ordered_values,point[search_axis]+tolerance)
            for position in lower:upper_bound
                candidate=column_order[position]
                maximum(abs(columns[end,candidate][axis]-point[axis]) for axis in 1:3)<=tolerance || continue
                target==0 || throw(ArgumentError(
                    "$caller: retained cap node has ambiguous source columns"))
                target=Int32(candidate)
            end
        end
        target>0 || throw(ArgumentError("$caller: retained cap node has no actual source column"))
        upper[target]==0 || throw(ArgumentError(
            "$caller: distinct retained cap nodes share a source column"))
        projection[row]=target;upper[target]=Int32(row)
    end
    all(>(0),upper) || throw(ArgumentError("$caller: retained cap omits an actual source column"))
    states=Vector{UInt8}(undef,size(source,2))
    for cell in axes(source,2)
        quad=@view source[:,cell]
        first=projection[triangles[:,2cell-1]];second=projection[triangles[:,2cell]]
        allunique(first) && allunique(second) && Set(vcat(first,second))==Set(quad) ||
            throw(ArgumentError("$caller: actual source/cap vector pairs do not correspond at Quad $cell"))
        diagonal=intersect(Set(first),Set(second))
        state=diagonal==Set((quad[1],quad[3])) ? UInt8(1) :
              diagonal==Set((quad[2],quad[4])) ? UInt8(2) : UInt8(0)
        state>0 || throw(ArgumentError("$caller: retained cap pair has no source diagonal at Quad $cell"))
        states[cell]=state
    end
    return _ModelMeshCapInput(states,cap,upper)
end

# A retained Surface constructor is rebound only by the actual SetOrder1 cell
# events. UID order belongs to the old receipt; current corner winding belongs
# to the current actual product. Neither coordinate sorting nor label offsets
# establish this operation's identity map.
function _model_mesh_first_order_surface_constructor_inputs(previous,current,element_mapping;
        caller="SetOrder1 Surface constructor")
    previous isa _ModelMeshAuthority && current isa _ModelMeshAuthority &&
        previous.mesh!==nothing && current.mesh!==nothing || throw(ArgumentError(
            "$caller: constructor transition lacks its actual before/after authority"))
    expected=UInt64[]
    for constructor in values(previous.surface_constructor_inputs),labels in constructor.entity_data.external_element_tags
        append!(expected,labels)
    end
    allunique(expected) && Set(keys(element_mapping))==Set(expected) &&
        all(key->key isa UInt64,keys(element_mapping)) &&
        all(value->value isa UInt64,values(element_mapping)) &&
        allunique(values(element_mapping)) || throw(ArgumentError(
            "$caller: constructor transition has missing, foreign or duplicate cell events"))
    rebound=Dict{Tuple{Int,Int32},MixedMesh}()
    for (owner,constructor) in previous.surface_constructor_inputs
        before=_model_mesh_product(previous,2,Int(owner[2]);caller)
        after=_model_mesh_product(current,2,Int(owner[2]);caller)
        before!==nothing && after!==nothing || throw(ArgumentError(
            "$caller: constructor transition loses actual Surface $owner"))
        _model_surface_constructor_binding(before,constructor,caller)
        copy_constructor=deepcopy(constructor)
        for labels in copy_constructor.entity_data.external_element_tags
            labels.=[element_mapping[uid] for uid in labels]
        end
        _model_surface_constructor_binding(after,copy_constructor,caller)
        rebound[owner]=copy_constructor
    end
    return rebound
end
