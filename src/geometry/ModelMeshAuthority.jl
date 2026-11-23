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

# Tables derived from an authority's mesh, built when first asked for and kept
# with the authority, so that the producers of a pass share one instead of each
# building its own from every node. They are no part of what an authority
# means: two authorities of one mesh differ at most in what is cached here.
mutable struct _ModelMeshDerived
    label_rows::Union{Nothing,Dict{UInt64,Int}}
    label_max::Union{Nothing,UInt64}
end
_ModelMeshDerived()=_ModelMeshDerived(nothing,nothing)
Base.show(io::IO,::_ModelMeshDerived)=print(io,"_ModelMeshDerived()")

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
    derived::_ModelMeshDerived
end

# A reconstruction over another mesh starts without derived tables.
_ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,carrier_cells,surface_constructor_inputs,emitted_geometry)=
    _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,carrier_rows,
        producer_rows,carrier_cells,surface_constructor_inputs,emitted_geometry,
        _ModelMeshDerived())

# The row of every node label of the actual mesh. Read-only: it is shared by
# everything that holds this authority.
function _model_mesh_label_rows(authority::_ModelMeshAuthority)
    derived=authority.derived
    rows=derived.label_rows
    if rows===nothing
        rows=Dict(tag=>row for (row,tag) in
            enumerate(authority.mesh.entity_data.external_node_tags))
        derived.label_rows=rows
    end
    return rows
end

# The largest node label of the actual mesh, zero for an empty one.
function _model_mesh_label_max(authority::_ModelMeshAuthority)
    derived=authority.derived
    value=derived.label_max
    if value===nothing
        value=maximum(authority.mesh.entity_data.external_node_tags;init=UInt64(0))
        derived.label_max=value
    end
    return value
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
        carrier_rows,producer_rows,carrier_cells::Vector{Vector{UInt64}},
        surface_constructor_inputs=Dict{Tuple{Int,Int32},MixedMesh}())
    bindings=carrier===mesh ? _model_mesh_live_cell_tags(carrier_cells) :
        [Union{Nothing,UInt64}[tag==0 ? nothing : tag for tag in block] for block in carrier_cells]
    return _ModelMeshAuthority(mesh,products,done,private_rank,native_inputs,carrier,
        carrier_rows,producer_rows,bindings,surface_constructor_inputs)
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

# Where the cell with each label sits in `mesh`, as block<<32|column, over the
# blocks `_model_mesh_primary_catalog` lists (a repeated label keeps its last
# cell, as that catalog does). Labels are dense after a renumbering, so the
# table is a vector unless they are sparse; no key and no node vector is made
# per cell.
function _model_mesh_cell_positions(mesh::MixedMesh)
    data=mesh.entity_data
    highest=UInt64(0);lowest=typemax(UInt64);total=0
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock && _msh_primary_nodes(msh_family(block.msh))!=0 || continue
        labels=data.external_element_tags[bi]
        for column in axes(block.nodes,2)
            highest=max(highest,labels[column]);lowest=min(lowest,labels[column])
        end
        total+=size(block.nodes,2)
    end
    dense=total>0 && lowest>0 && highest<=UInt64(4)*UInt64(total)+UInt64(1024)
    table=zeros(UInt64,dense ? Int(highest) : 0)
    sparse=Dict{UInt64,UInt64}()
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock && _msh_primary_nodes(msh_family(block.msh))!=0 || continue
        labels=data.external_element_tags[bi]
        for column in axes(block.nodes,2)
            packed=(UInt64(bi)<<32)|UInt64(column)
            if dense
                table[labels[column]]=packed
            else
                sparse[labels[column]]=packed
            end
        end
    end
    return table,sparse,dense
end

# One carrier block of `_model_mesh_refined_cell_tags`: `false` when a parent
# is not the actual cell of its label on the same rows in the same order.
function _model_mesh_refined_block_tags!(tags::Vector{Union{Nothing,UInt64}},
        sources::AbstractVector,cells::AbstractMatrix{<:Integer},rows::AbstractVector{<:Integer},
        mesh::MixedMesh,table::Vector{UInt64},sparse::Dict{UInt64,UInt64},dense::Bool,
        produced_by,children::Int,caller)
    for parent in axes(cells,2)
        source=sources[parent]
        source===nothing && return false
        packed=dense ? (0<source<=length(table) ? table[source] : UInt64(0)) :
            get(sparse,source,UInt64(0))
        packed==0 && return false
        block=mesh.blocks[Int(packed>>32)]::ElementBlock
        column=Int(packed&0xffffffff)
        width=_msh_primary_nodes(msh_family(block.msh))
        nodes=block.nodes
        width==size(cells,1) || return false
        for slot in 1:width
            nodes[slot,column]==rows[cells[slot,parent]] || return false
        end
        produced=get(produced_by,source,nothing)
        produced!==nothing && length(produced)==children || throw(ArgumentError(
            "$caller: refined carrier parent lacks its actual child events"))
        append!(tags,produced)
    end
    return true
end

# The labels of a refined carrier's cells: each parent's child labels, parent
# after parent, when every parent is the actual cell of its label. A parent is
# found through a table from label to cell position and compared in place;
# the catalog of the whole previous mesh this used to build (a key, three
# vectors and two table entries per cell) cost more than the labels.
function _model_mesh_refined_cell_tags(plan,previous;caller="mesh refinement")
    isempty(previous.carrier_cells) && return Vector{Union{Nothing,UInt64}}[]
    previous.carrier isa Mesh && hasproperty(plan,:refine_cell_tags) || return Vector{Union{Nothing,UInt64}}[]
    # Subdivision into quadrangles or hexahedra leaves no simplex carrier.
    get(plan,:refine_split,0)==0 || return Vector{Union{Nothing,UInt64}}[]
    table,sparse,dense=_model_mesh_cell_positions(previous.mesh)
    output=Vector{Union{Nothing,UInt64}}[]
    for (bi,(msh,cells)) in enumerate(_model_mesh_carrier_blocks(previous.carrier))
        children=msh==1 ? 2 : msh==2 ? 4 : 8
        tags=Union{Nothing,UInt64}[];sizehint!(tags,Base.checked_mul(children,size(cells,2)))
        # A proxy/permuted parent needs its own child trace.
        _model_mesh_refined_block_tags!(tags,previous.carrier_cells[bi],cells,previous.carrier_rows,
            previous.mesh,table,sparse,dense,plan.refine_cell_tags,children,caller) ||
            return Vector{Union{Nothing,UInt64}}[]
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

# One entity's product: the rows it owns or references, in private-rank order,
# and its cells per block. `selected` lists (block, columns) in block order.
#
# `work` is a zeroed table of one entry per mesh row, lent by a caller that
# builds every product of a mesh: rows are then marked and numbered in it
# instead of in a set and a dictionary, and it is zeroed again on return. A
# large product (a Region of half a million cells) spends its time in those
# hashes otherwise. Without `work`, a product costs only its own size, which
# an append to one entity among many needs.
function _model_mesh_owner_product(mesh::MixedMesh,owner,owned,selected,private_rank,caller;
                                   work::Union{Nothing,Vector{Int32}}=nothing)
    data=mesh.entity_data
    if work!==nothing
        rows=Int32[]
        for row in owned
            work[row]==0 || continue
            work[row]=Int32(-1);push!(rows,row)
        end
        for (bi,columns) in selected
            nodes=mesh.blocks[bi].nodes
            for column in columns,slot in axes(nodes,1)
                node=nodes[slot,column]
                0<node<=nnodes(mesh) || throw(ArgumentError(
                    "$caller: entity $owner references an invalid actual node row"))
                work[node]==0 || continue
                work[node]=Int32(-1);push!(rows,node)
            end
        end
        sort!(rows;by=row->private_rank[row])
        for (index,row) in enumerate(rows)
            work[row]=Int32(index)
        end
        blocks=ElementBlock[];tags=Vector{UInt64}[]
        for (bi,columns) in selected
            block=mesh.blocks[bi];nodes=block.nodes
            local_nodes=Matrix{Int32}(undef,size(nodes,1),length(columns))
            for (position,column) in enumerate(columns),slot in axes(nodes,1)
                local_nodes[slot,position]=work[nodes[slot,column]]
            end
            push!(blocks,ElementBlock(Elements._TRUSTED_MIXED_STORAGE,block.msh,local_nodes,
                block.tags[columns]))
            push!(tags,data.external_element_tags[bi][columns])
        end
        for row in rows
            work[row]=Int32(0)
        end
        return _ModelMeshProduct(owner,rows,blocks,tags)
    end
    referenced=Set(owned)
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
    return _ModelMeshProduct(owner,rows,blocks,tags)
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
    # Nodes and cells arrive entity by entity: the table entry of a run of one
    # owner is looked up once.
    run_owner=(-1,Int32(-1));run_rows=Int32[]
    for (row,owner) in enumerate(data.node_entities)
        if owner!=run_owner
            run_owner=owner
            run_rows=get!(()->Int32[],owned_rows,owner)
            get!(()->Tuple{Int,Vector{Int}}[],selections,owner)
        end
        push!(run_rows,Int32(row))
    end
    seen_tags=Elements._SeenTags(data.external_element_tags)
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock || throw(ArgumentError("$caller: unsupported mesh block"))
        owners=data.block_entities[bi];tags=data.external_element_tags[bi]
        length(owners)==size(block.nodes,2)==length(tags) || throw(ArgumentError(
            "$caller: actual cell metadata does not align with connectivity"))
        for tag in tags
            Elements._seen_tag!(seen_tags,tag) && throw(ArgumentError(
                "$caller: actual element labels must be positive, unique and Int32 representable"))
        end
        by_owner=Dict{Int32,Vector{Int}}()
        cell_owner=Int32(-1);cell_columns=Int[]
        for (column,owner) in enumerate(owners)
            if owner!=cell_owner || column==1
                cell_owner=owner
                cell_columns=get!(()->Int[],by_owner,owner)
            end
            push!(cell_columns,column)
        end
        for (owner,columns) in by_owner
            key=(msh_dimension(block.msh),owner)
            push!(get!(()->Tuple{Int,Vector{Int}}[],selections,key),(bi,columns))
        end
    end
    products=Dict{Tuple{Int,Int32},_ModelMeshProduct}()
    work=zeros(Int32,nnodes(mesh))
    for (owner,selected) in selections
        products[owner]=_model_mesh_owner_product(mesh,owner,
            get(owned_rows,owner,Int32[]),selected,private_rank,caller;work)
    end
    completed=Set{Tuple{Int,Int32}}((Int(dim),Int32(tag)) for (dim,tag) in done)
    bindings=_model_mesh_carrier_cell_tags(mesh,carrier,carrier_cells,caller)
    return _ModelMeshAuthority(mesh,products,completed,UInt64.(private_rank),
        Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}(),carrier,Int32.(carrier_rows),
        Dict{Tuple{Int,Int32},Vector{Int32}}(),bindings)
end

include("ModelMeshAuthorityLog.jl")

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

# A retained Surface constructor receipt is a first-order record. An order
# pass elevates the actual product in place, but the receipt still binds its
# first-order skeleton — every elevated cell reduced to its family's corner
# support. Rows an elevated cell references only through a non-corner slot
# (midnodes) leave the skeleton; rows no cell references (owned orphan
# vertices) stay bound. Returns `nothing` when an elevated block has no
# first-order projection.
function _model_mesh_first_order_layout(nrows::Integer,blocks)
    any(block->msh_spec(block.msh).order>1,blocks) || return nothing,blocks
    primary=falses(nrows);incident=falses(nrows);elevated=ElementBlock[]
    for block in blocks
        spec=msh_spec(block.msh)
        for node in block.nodes
            incident[node]=true
        end
        if spec.order>1
            lin=try
                Elements.msh_type(spec.family,1)
            catch err
                err isa ArgumentError || rethrow()
                return nothing,nothing
            end
            width=msh_num_nodes(lin)
            for node in @view block.nodes[1:width,:]
                primary[node]=true
            end
            push!(elevated,ElementBlock(lin,block.nodes[1:width,:],block.tags))
        else
            for node in block.nodes
                primary[node]=true
            end
            push!(elevated,block)
        end
    end
    keep=findall(primary .| .!incident)
    remap=zeros(Int32,nrows)
    for (position,row) in enumerate(keep)
        remap[row]=Int32(position)
    end
    return keep,ElementBlock[ElementBlock(block.msh,remap[block.nodes],block.tags)
                             for block in elevated]
end

function _model_mesh_first_order_skeleton(actual::MixedMesh,caller)
    data=actual.entity_data
    data===nothing && return actual
    keep,blocks=_model_mesh_first_order_layout(nnodes(actual),actual.blocks)
    blocks===nothing && throw(ArgumentError(
        "$caller: Surface constructor input requires actual linear Surface cells"))
    keep===nothing && return actual
    rows=keep
    metadata=Elements._native_mixed_entity_data(data.entities;
        node_entities=data.node_entities[rows],
        node_parametric=data.node_parametric[rows],
        external_node_tags=data.external_node_tags[rows],
        block_entities=data.block_entities,
        external_element_tags=data.external_element_tags)
    return MixedMesh(actual.coords[:,rows],blocks;
        physical_names=actual.physical_names,entity_data=metadata,
        elementary_entities=actual.elementary_entities)
end

# The actual mesh of one entity, cut from the authority's mesh.
#
# The default cut is detached: every array, the entity records and the entity
# table are its own, and its structure is checked. A caller may keep it, hand
# it on or edit it.
#
# `shared=true` is for a caller that reads the cut and drops it. The rows of
# the entity are still sliced into arrays of the cut's own (they are the size
# of the entity's mesh), but everything that does not depend on the rows'
# positions is the authority's own storage: the entity records, each node's
# parameter vector, the cell labels, the physical names, and the blocks when
# the cut keeps the product's order. Its entity table is the one table built
# to select the records, and nothing is checked again: the rows come from a
# mesh that was validated when the authority was built.
#
# What the caller of a shared cut must not do:
#   - edit it in place: not its coordinates and not any of its tables or
#     arrays. Stored meshes ARE edited in place elsewhere
#     (`set_point_coordinates!` writes the coordinates of every mesh in
#     `native_inputs` and `surface_constructor_inputs`, the coordinate
#     publishers write `context.mesh_parts`), so a shared cut must never
#     become a stored one;
#   - keep it beyond the generation that produced the authority: not in a
#     struct field, a plan, `native_inputs`, `surface_constructor_inputs` or
#     `mesh_parts`, and not as a return value. Stored cuts are compared with
#     `===` and used as identity keys, and they are expected to stay as they
#     were while the authority's mesh is edited in place;
#   - hand it to anything that does either.
# A caller that needs any of that takes the default cut. For a Curve or a
# Point read row by row, `_model_mesh_product_rows` builds nothing at all.
function _model_mesh_product(authority::_ModelMeshAuthority,dim::Int,tag::Int;
                             caller="mesh authority",producer_order=false,
                             first_order=false,shared::Bool=false)
    product=get(authority.products,(dim,Int32(tag)),nothing)
    product===nothing && return nothing
    source=authority.mesh;data=source.entity_data;rows=product.rows
    blocks=product.blocks
    if first_order
        keep,blocks=_model_mesh_first_order_layout(length(rows),blocks)
        blocks===nothing && throw(ArgumentError(
            "$caller: actual entity $(product.owner) lacks a first-order projection"))
        keep===nothing || (rows=rows[keep])
    end
    produced=producer_order ? get(authority.producer_rows,product.owner,nothing) : nothing
    if produced!==nothing
        # Normal PLC input follows the actual face producer vector. NoNew
        # consumers keep the separate private-rank order above.
        # The inverse covers the producer rows only: a mesh-sized table per
        # product makes a pass over many entities quadratic.
        checkbounds(source.coords,1,produced)
        inverse=Dict{Int32,Int32}();sizehint!(inverse,length(produced))
        for (index,row) in enumerate(produced)
            inverse[row]=Int32(index)
        end
        # Product rows also list owned vertices with no cell incidence (an
        # AddNodes orphan or a retained predecessor row). They have no
        # connectivity to reorder and drop out through rows=produced below;
        # only a used row missing its producer binding invalidates the source.
        for block in blocks,node in block.nodes
            get(inverse,rows[node],Int32(0))>0 || throw(ArgumentError(
                "$caller: actual entity $(product.owner) changed without a complete producer row binding"))
        end
        mapping=Int32[get(inverse,row,Int32(0)) for row in rows]
        blocks=ElementBlock[ElementBlock(block.msh,mapping[block.nodes],block.tags)
                            for block in blocks]
        rows=produced
    end
    entities=_model_mesh_product_entities(source,product.owner,rows)
    if shared
        cut=Elements.MixedEntityData(Elements._OWNED_MIXED_ENTITY_DATA,entities,
            data.node_entities[rows],data.node_parametric[rows],data.external_node_tags[rows],
            Vector{Int32}[fill(Int32(tag),size(block.nodes,2)) for block in blocks],
            product.element_tags,true)
        return MixedMesh(Elements._TRUSTED_MIXED_STORAGE,source.coords[:,rows],
            MixedElementBlock[block for block in blocks],source.physical_names,cut,
            Vector{Int32}[fill(Int32(tag),size(block.nodes,2)) for block in blocks],
            Elements.MixedPeriodicLink[],Elements.MshAncillarySection[],
            Elements.MshDataSection[],nothing)
    end
    # Every part is cut from the validated actual mesh: copy each once and
    # check the result once, instead of the public constructors' three copies.
    return Elements._detached_mixed_mesh(source.coords[:,rows],
        MixedElementBlock[ElementBlock(Elements._TRUSTED_MIXED_STORAGE,block.msh,
            copy(block.nodes),copy(block.tags)) for block in blocks],
        source.physical_names,entities,data.node_entities[rows],
        Elements._detached_mixed_parameters(@view data.node_parametric[rows]),
        data.external_node_tags[rows],
        Vector{Int32}[fill(Int32(tag),size(block.nodes,2)) for block in blocks],
        Vector{UInt64}[copy(tags) for tags in product.element_tags],
        Vector{Int32}[fill(Int32(tag),size(block.nodes,2)) for block in blocks],true)
end

# An entity's actual product read in place: no mesh is built and nothing is
# checked again, because the pass that produced the entity validated it.
#
#     view=_model_mesh_product_rows(authority,dim,tag)   # `nothing` when absent
#     view.mesh           the authority's mesh; `view.mesh.entity_data` holds
#                         each row's label, owner and parameters
#     view.rows[i]        the row of `view.mesh` that is the product's node i
#     view.blocks         the product's cells; a node index is a position in
#                         `view.rows`
#     view.element_tags   the cells' labels, one vector per block
#
# Node i here is node i of `_model_mesh_product(authority,dim,tag)` in its
# default order, and the blocks and labels are that cut's too. Everything
# returned is the authority's own storage: read it, do not edit or keep it
# beyond the authority.
function _model_mesh_product_rows(authority::_ModelMeshAuthority,dim::Int,tag::Int)
    product=get(authority.products,(dim,Int32(tag)),nothing)
    product===nothing && return nothing
    return (mesh=authority.mesh::MixedMesh,rows=product.rows,blocks=product.blocks,
        element_tags=product.element_tags)
end

# The node labels of one linear triangle or quadrangle in increasing order,
# padded to four: the cell's support without a vector per cell.
@inline function _model_surface_cell_support(tags::Vector{UInt64},nodes::Matrix{Int32},column::Int)
    a=tags[nodes[1,column]];b=tags[nodes[2,column]];c=tags[nodes[3,column]]
    a,b=minmax(a,b);b,c=minmax(b,c);a,b=minmax(a,b)
    size(nodes,1)==3 && return (a,b,c,UInt64(0))
    d=tags[nodes[4,column]]
    c,d=minmax(c,d);b,c=minmax(b,c);a,b=minmax(a,b)
    return (a,b,c,d)
end

# The two halves of the binding check for the usual case, in which the input
# lists the actual product's nodes and cells in the same order: row against
# row and column against column, with no table from label to position. Each
# returns `true` only when the general check below would pass. `false` decides
# nothing: the general check then runs and reports what is wrong.
function _model_surface_nodes_aligned(actual::MixedMesh,constructor::MixedMesh)
    data=actual.entity_data;receipt=constructor.entity_data
    labels=data.external_node_tags
    count=length(labels)
    count==nnodes(actual)==nnodes(constructor) && labels==receipt.external_node_tags || return false
    (issorted(labels;lt=<=) || allunique(labels)) || return false
    for row in 1:count
        all(axis->isequal(constructor.coords[axis,row],actual.coords[axis,row]),1:3) &&
            receipt.node_entities[row]==data.node_entities[row] &&
            isequal(receipt.node_parametric[row],data.node_parametric[row]) || return false
    end
    return true
end

function _model_surface_cells_aligned(actual::MixedMesh,constructor::MixedMesh)
    data=actual.entity_data;receipt=constructor.entity_data
    count=length(actual.blocks)
    count==length(data.external_element_tags)==length(data.block_entities)==
        length(constructor.blocks)==length(receipt.external_element_tags)==
        length(receipt.block_entities) || return false
    previous=UInt64(0);started=false;increasing=true;total=0
    for index in 1:count
        block=actual.blocks[index];other=constructor.blocks[index]
        labels=data.external_element_tags[index];owners=data.block_entities[index]
        other_labels=receipt.external_element_tags[index];other_owners=receipt.block_entities[index]
        block isa ElementBlock && other isa ElementBlock || return false
        block.msh in (2,3) && block.msh==other.msh || return false
        cells=size(block.nodes,2)
        cells==size(other.nodes,2)==length(labels)==length(other_labels)==
            length(owners)==length(other_owners) || return false
        for column in 1:cells
            label=labels[column]
            label==other_labels[column] || return false
            started && label<=previous && (increasing=false)
            started=true;previous=label
            owners[column]==other_owners[column] && block.tags[column]==other.tags[column] &&
                _model_surface_cell_support(data.external_node_tags,block.nodes,column)==
                _model_surface_cell_support(receipt.external_node_tags,other.nodes,column) ||
                return false
        end
        total+=cells
    end
    # Both meshes carry the same labels at the same positions. Increasing
    # labels cannot repeat; birth labels in another order are checked for it.
    increasing && return true
    seen=Set{UInt64}();sizehint!(seen,total)
    for labels in data.external_element_tags,label in labels
        label in seen && return false
        push!(seen,label)
    end
    return true
end

function _model_surface_constructor_binding(actual::MixedMesh,constructor::MixedMesh,
                                             caller::AbstractString)
    actual=_model_mesh_first_order_skeleton(actual,caller)
    actual.entity_data!==nothing && constructor.entity_data!==nothing ||
        throw(ArgumentError("$caller: Surface constructor input lacks actual identity metadata"))
    data=actual.entity_data;receipt=constructor.entity_data
    length(data.node_entities)==length(data.node_parametric)==nnodes(actual) &&
        length(receipt.node_entities)==length(receipt.node_parametric)==nnodes(constructor) ||
        throw(ArgumentError("$caller: Surface constructor node metadata does not align"))
    _model_surface_nodes_aligned(actual,constructor) &&
        _model_surface_cells_aligned(actual,constructor) && return nothing
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
    cells=Dict{UInt64,Tuple{Int,Int32,Int32,NTuple{4,UInt64}}}()
    sizehint!(cells,sum(length,data.external_element_tags;init=0))
    length(actual.blocks)==length(data.external_element_tags)==length(data.block_entities) ||
        throw(ArgumentError("$caller: actual Surface cell metadata does not align"))
    for (block,labels,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)
        block.msh in (2,3) && length(labels)==length(owners)==size(block.nodes,2) || throw(ArgumentError(
            "$caller: Surface constructor input requires actual linear Surface cells"))
        for (column,label) in enumerate(labels)
            haskey(cells,label) && throw(ArgumentError(
                "$caller: actual Surface cells repeat constructor UID $label"))
            support=_model_surface_cell_support(data.external_node_tags,block.nodes,column)
            cells[label]=(Int(block.msh),owners[column],block.tags[column],support)
        end
    end
    seen=Set{UInt64}();sizehint!(seen,length(cells))
    length(constructor.blocks)==length(receipt.external_element_tags)==length(receipt.block_entities) ||
        throw(ArgumentError("$caller: Surface constructor cell labels do not align"))
    for (block,labels,owners) in zip(constructor.blocks,receipt.external_element_tags,receipt.block_entities)
        length(labels)==length(owners)==size(block.nodes,2) || throw(ArgumentError(
            "$caller: Surface constructor cell labels do not align"))
        for (column,label) in enumerate(labels)
            label in seen && throw(ArgumentError(
                "$caller: Surface constructor input repeats cell UID $label"))
            push!(seen,label)
            # Actual cells are linear triangles and quadrangles; any other
            # constructor cell has no equal among them.
            block.msh in (2,3) && get(cells,label,nothing)==(Int(block.msh),owners[column],
                block.tags[column],_model_surface_cell_support(receipt.external_node_tags,
                    block.nodes,column)) || throw(ArgumentError(
                "$caller: Surface constructor input changes cell $label owner or support"))
        end
    end
    length(seen)==length(cells) || throw(ArgumentError(
        "$caller: Surface constructor input omits actual cells"))
    return nothing
end

# `adopt=true`: `actual` is a shared cut made for this call and dropped after
# it. The input then keeps the cut's own arrays (coordinates, node owners and
# labels, cell owners) instead of copying them, shares the authority's entity
# records through the cut's table, and copies only what the cut itself shares
# with the authority and an input must own (the parameter vectors). Its
# structure is not checked separately: the binding check that ends this
# function compares every node and every cell with the validated product.
function _model_surface_constructor_input(production,actual::MixedMesh,
                                           caller::AbstractString;
                                           constructor=production.constructor_mesh,
                                           adopt::Bool=false)
    constructor===nothing && return nothing
    data=actual.entity_data
    data!==nothing || throw(ArgumentError(
        "$caller: Surface constructor input lacks the authoritative product"))
    length(production.node_tags)==nnodes(constructor) &&
        allunique(production.node_tags) || throw(ArgumentError(
        "$caller: Surface constructor rows do not bind unique actual UIDs"))
    # The producer's rows usually are the product's rows, in order: then the
    # row of each (unique) label is its own position.
    mapping=if production.node_tags==data.external_node_tags
        collect(Int32(1):Int32(length(production.node_tags)))
    else
        actual_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
        all(uid->haskey(actual_rows,uid),production.node_tags) || throw(ArgumentError(
            "$caller: Surface constructor row UID is absent from its actual product"))
        Int32[actual_rows[uid] for uid in production.node_tags]
    end
    for row in eachindex(mapping)
        bound=mapping[row]
        all(axis->isequal(constructor.coords[axis,row],actual.coords[axis,bound]),1:3) &&
            production.node_owners[row]==data.node_entities[bound] &&
            isequal(production.node_parameters[row],data.node_parametric[bound]) ||
            throw(ArgumentError("$caller: Surface constructor row loses its actual node binding"))
    end
    # Physical tag by cell label. A constructor block that carries the labels
    # of the product block at its position takes that block's tags as they
    # are; the table from label to tag is built only for a block that does
    # not. (Labels repeated in the product fail the binding check below under
    # either reading.)
    physical_tags=nothing
    blocks=ElementBlock[];labels=Vector{UInt64}[];cursor=1
    for (msh,nodes) in _surface_dim2_blocks(constructor)
        count=size(nodes,2)
        cursor+count-1<=length(production.element_tags) || throw(ArgumentError(
            "$caller: Surface constructor cells lack actual birth labels"))
        cell_labels=production.element_tags[cursor:cursor+count-1]
        position=length(blocks)+1
        tags=if physical_tags===nothing && position<=length(actual.blocks) &&
                position<=length(data.external_element_tags) &&
                length(actual.blocks[position].tags)==count &&
                data.external_element_tags[position]==cell_labels
            Int32[tag for tag in actual.blocks[position].tags]
        else
            physical_tags===nothing && (physical_tags=Dict(uid=>block.tags[column] for (block,uids) in
                zip(actual.blocks,data.external_element_tags) for (column,uid) in enumerate(uids)))
            all(uid->haskey(physical_tags,uid),cell_labels) || throw(ArgumentError(
                "$caller: Surface constructor cell UID is absent from its actual product"))
            Int32[physical_tags[uid] for uid in cell_labels]
        end
        push!(blocks,ElementBlock(msh,mapping[nodes],tags))
        push!(labels,cell_labels)
        cursor+=count
    end
    cursor==length(production.element_tags)+1 || throw(ArgumentError(
        "$caller: Surface constructor cells do not cover their birth labels"))
    elementary=actual.elementary_entities
    result=if elementary===nothing
        metadata=Elements._native_mixed_entity_data(data.entities;
            node_entities=copy(data.node_entities),node_parametric=deepcopy(data.node_parametric),
            external_node_tags=copy(data.external_node_tags),
            block_entities=deepcopy(data.block_entities),external_element_tags=labels)
        MixedMesh(copy(actual.coords),blocks;physical_names=copy(actual.physical_names),
            entity_data=metadata)
    elseif adopt
        MixedMesh(Elements._TRUSTED_MIXED_STORAGE,actual.coords,
            MixedElementBlock[block for block in blocks],
            Elements._copy_physical_names(actual.physical_names,"MixedMesh"),
            Elements.MixedEntityData(Elements._OWNED_MIXED_ENTITY_DATA,data.entities,
                data.node_entities,Elements._detached_mixed_parameters(data.node_parametric),
                data.external_node_tags,data.block_entities,labels,true),
            elementary,Elements.MixedPeriodicLink[],Elements.MshAncillarySection[],
            Elements.MshDataSection[],nothing)
    else
        # The same mesh from one copy of each part: `actual` is a validated
        # product and the blocks above passed their own constructor.
        Elements._detached_mixed_mesh(copy(actual.coords),
            MixedElementBlock[block for block in blocks],actual.physical_names,data.entities,
            copy(data.node_entities),Elements._detached_mixed_parameters(data.node_parametric),
            copy(data.external_node_tags),
            Vector{Int32}[copy(owners) for owners in data.block_entities],labels,
            Vector{Int32}[copy(owners) for owners in elementary],true)
    end
    _model_surface_constructor_binding(actual,result,caller)
    return result
end

function _model_mesh_with_surface_constructor_inputs(authority::_ModelMeshAuthority,inputs;
        caller="Surface constructor input",only=nothing,shared::Bool=false,known::Bool=false,
        validated=nothing,certified::Bool=false,trusted=nothing)
    if shared
        # An authority taken from a log already shares its entity tables across
        # appends, so the pass's own table is bound as it is instead of being
        # copied per Surface. `known` says the other owners were checked
        # against this log before; its products only grow, so they still hold.
        only!==nothing && inputs isa Dict{Tuple{Int,Int32},MixedMesh} || throw(ArgumentError(
            "$caller: a shared Surface constructor table needs the one new owner"))
        for (owner,constructor) in (known ? (haskey(inputs,only) ? (only=>inputs[only],) : ()) : inputs)
            owner[1]==2 && haskey(authority.products,owner) || throw(ArgumentError(
                "$caller: Surface constructor input refers to absent actual entity $owner"))
            # `certified`: the caller built this input from the same product
            # view a moment ago, and that construction ends with this check.
            if owner==only && !certified
                actual=_model_mesh_product(authority,2,Int(owner[2]);caller,
                    producer_order=true,first_order=true,shared=true)
                _model_surface_constructor_binding(actual,constructor,caller)
            end
        end
        return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
            authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,
            authority.producer_rows,authority.carrier_cells,inputs,nothing,authority.derived)
    end
    bound=Dict{Tuple{Int,Int32},MixedMesh}()
    for (owner,constructor) in inputs
        owner[1]==2 && haskey(authority.products,owner) || throw(ArgumentError(
            "$caller: Surface constructor input refers to absent actual entity $owner"))
        # Constructor inputs are bound against the producer-order product at
        # creation time; re-validate against the same producer view so unused
        # owned rows cannot change the actual node set comparison. The
        # first-order view keeps an order-elevated product boundable through
        # its retained corner skeleton. With `only`, just that owner's input
        # is validated: a caller that binds one Surface at a time validates
        # all of them once it is done, rather than all of them every time.
        # `validated` names, per owner, the product, producer rows and input a
        # log authority was already checked with. While they are the very same
        # objects the log has only been appended to since, so the product view
        # and the check's answer are the same too.
        record=validated===nothing ? nothing : get(validated,owner,nothing)
        # `trusted` names owners whose binding the caller carries over from
        # an earlier check of the same cells under a one-to-one relabelling.
        if (only===nothing || owner==only) && !(trusted!==nothing && owner in trusted) &&
                !(record!==nothing && record[1]===authority.products[owner] &&
                record[2]===get(authority.producer_rows,owner,nothing) && record[3]===constructor)
            actual=_model_mesh_product(authority,2,Int(owner[2]);caller,
                producer_order=true,first_order=true,shared=true)
            _model_surface_constructor_binding(actual,constructor,caller)
        end
        bound[owner]=constructor
    end
    return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
        authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,
        authority.producer_rows,authority.carrier_cells,bound,nothing,authority.derived)
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
        authority.carrier_cells,authority.surface_constructor_inputs,nothing,authority.derived)
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
    # The mesh and its products are unchanged, so retained Surface constructor
    # inputs stay bound to the same actual products.
    return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
        authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,rows,
        authority.carrier_cells,authority.surface_constructor_inputs,nothing,authority.derived)
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
    # The carrier was validated when its authority was built, and every child
    # below must be one actual cell of the refined mesh in an even
    # permutation, which no invalid child can be. The trace is run for its
    # node and cell order alone, so it validates neither its input nor its
    # result (1.7 s of 2.2 s on 497,124 tetrahedra).
    trace=Refine._refine_uniform_traced(carrier;require_positive_tets=false,trusted=true)
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
    # Node provenance alone does not certify a carrier's new cells: each must be
    # the only cell of the refined mesh on its primary nodes, and an even
    # permutation of it. Node labels name the rows one to one, so two cells
    # have the same labels exactly when they have the same rows. A cell is
    # looked up among the cells that share its largest row (a refined cell's
    # largest row is a midpoint, which few cells share) and that have the same
    # row sum: nothing is hashed and nothing is allocated per cell.
    labels=plan.mesh.entity_data.external_node_tags
    if allunique(labels)
        function by_largest_row(tables::Vector{Matrix{Int32}},width::Int,node_count::Int)
            starts=zeros(Int32,node_count+2);total=0
            @inbounds for nodes in tables,column in axes(nodes,2)
                high=nodes[1,column]
                for row in 2:width
                    high=max(high,nodes[row,column])
                end
                starts[high+2]+=Int32(1);total+=1
            end
            @inbounds for node in 3:node_count+2
                starts[node]+=starts[node-1]
            end
            table_of=Vector{Int32}(undef,total);column_of=Vector{Int32}(undef,total)
            sum_of=Vector{Int64}(undef,total)
            @inbounds for (table,nodes) in enumerate(tables),column in axes(nodes,2)
                high=nodes[1,column];sum=Int64(high)
                for row in 2:width
                    high=max(high,nodes[row,column]);sum+=nodes[row,column]
                end
                position=starts[high+1]+Int32(1);starts[high+1]=position
                table_of[position]=table;column_of[position]=column;sum_of[position]=sum
            end
            # Afterwards the cells whose largest row is `r` are starts[r]+1:starts[r+1].
            return starts,table_of,column_of,sum_of
        end
        function certified(tables::Vector{Matrix{Int32}},starts::Vector{Int32},
                table_of::Vector{Int32},column_of::Vector{Int32},sum_of::Vector{Int64},
                rows::Vector{Int32},cells::Matrix{Int32},width::Int)
            @inbounds for column in axes(cells,2)
                high=rows[cells[1,column]];sum=Int64(high)
                for row in 2:width
                    node=rows[cells[row,column]]
                    high=max(high,node);sum+=node
                end
                found=0;inversions=0
                for position in starts[high]+1:starts[high+1]
                    sum_of[position]==sum || continue
                    target=tables[table_of[position]];target_column=column_of[position]
                    # Where each carrier node sits in the candidate cell (0: nowhere).
                    p1=0;p2=0;p3=0;p4=0
                    for row in 1:width
                        node=rows[cells[row,column]];spot=0
                        for slot in 1:width
                            if target[slot,target_column]==node
                                spot=slot
                                break
                            end
                        end
                        row==1 ? (p1=spot) : row==2 ? (p2=spot) : row==3 ? (p3=spot) : (p4=spot)
                    end
                    (p1==0 || p2==0 || (width>=3 && p3==0) || (width>=4 && p4==0)) && continue
                    found+=1
                    inversions=(p1>p2)+(width>=3 ? (p1>p3)+(p2>p3) : 0)+
                        (width>=4 ? (p1>p4)+(p2>p4)+(p3>p4) : 0)
                end
                (found==1 && iseven(inversions)) || return false
            end
            return true
        end
        for (family,cells) in ((:lin,trace.mesh.segs),(:tri,trace.mesh.tris),(:tet,trace.mesh.tets))
            isempty(cells) && continue
            width=size(cells,1)
            tables=Matrix{Int32}[block.nodes for block in plan.mesh.blocks
                if block isa ElementBlock && msh_family(block.msh)===family]
            starts,table_of,column_of,sum_of=by_largest_row(tables,width,nnodes(plan.mesh))
            certified(tables,starts,table_of,column_of,sum_of,rows,cells,width) ||
                return nothing,Int32[]
        end
        refined=trace.mesh
        output=Mesh(plan.mesh.coords[:,rows];segs=refined.segs,tris=refined.tris,tets=refined.tets,
            seg_tag=refined.seg_tag,tri_tag=refined.tri_tag,tet_tag=refined.tet_tag)
        return output,rows
    end
    catalog,actual=_model_mesh_primary_catalog(plan.mesh)
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
    result=_model_mesh_with_producer_rows(authority,produced)
    # Constructor receipts stay bound only while the rebound product keeps the
    # same node/cell identities; retain each receipt that still validates and
    # drop the ones the new mesh genuinely invalidated.
    for (owner,mesh) in previous.surface_constructor_inputs
        try
            result=_model_mesh_with_surface_constructor_inputs(result,
                Dict(owner=>mesh))
        catch
            continue
        end
    end
    return result
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

function _model_mesh_curve_births(authority,model,curve::Int,caller;used_only=false,
                                  producer_order=true)
    authority===nothing && return nothing
    # A re-generated Curve leaves its superseded owned vertices in the actual
    # mesh until `_dim01_finish!` drops them. The bound producer rows name the
    # emitted chain, so read the emitted product view — owned orphans must not
    # join the open-chain check. Sweep consumers request the full product and
    # validate its unused owned rows separately.
    mesh=_model_mesh_product(authority,1,curve;caller,producer_order)
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
        product=_model_mesh_product(authority,0,point_tag;caller,shared=true)
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
    positions=Dict(uid=>Int32(row) for (row,uid) in
        enumerate(current.mesh.entity_data.external_node_tags))
    rebound=Dict{Tuple{Int,Int32},MixedMesh}()
    for (owner,constructor) in previous.surface_constructor_inputs
        before=_model_mesh_product(previous,2,Int(owner[2]);caller,
            producer_order=true,first_order=true,shared=true)
        before!==nothing || throw(ArgumentError(
            "$caller: constructor transition loses actual Surface $owner"))
        _model_surface_constructor_binding(before,constructor,caller)
        rows=Int32[get(positions,uid,Int32(0)) for uid in
            before.entity_data.external_node_tags]
        bound=_model_mesh_with_producer_rows(current,Dict(owner=>rows);caller)
        after=_model_mesh_product(bound,2,Int(owner[2]);caller,
            producer_order=true,first_order=true,shared=true)
        after!==nothing || throw(ArgumentError(
            "$caller: constructor transition loses actual Surface $owner"))
        # Only the cell labels change: the rest of the receipt is shared.
        copy_constructor=Elements._relabeled_mixed_mesh(constructor,
            constructor.entity_data.external_node_tags,
            Vector{UInt64}[UInt64[element_mapping[uid] for uid in labels]
                for labels in constructor.entity_data.external_element_tags])
        _model_surface_constructor_binding(after,copy_constructor,caller)
        rebound[owner]=copy_constructor
    end
    return rebound
end
