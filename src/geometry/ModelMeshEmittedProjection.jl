# The emitted constructor receipt authenticates a detached geometry view.
# Publish a fresh receipt after the final authority and emitted parts exist.
# Each later authority replacement needs a new receipt; an old receipt cannot
# authenticate a rebuilt/reclassified/current generation by shape alone.
#
# A receipt is published by every generation and read only when a detached
# part comes back for projection, so publishing it has to be cheap. Two things
# in it were not. The owner of every carrier cell was kept in one table keyed
# by `(msh, rows)` with rows of any length: an untyped key, boxed for each
# cell and hashed through dynamic dispatch, twice per cell. The tables below
# are one per cell arity, so their keys are plain bits. And the state of the
# actual mesh was kept as a deep copy of all of it (a vector per node for the
# parametric coordinates alone) only to be compared with the live mesh later;
# a fingerprint of every word of that state answers the same question.
const _EmittedOwner=Tuple{Int,Int32}
const _EMITTED_NO_OWNER=(-1,Int32(0))

struct _EmittedCellOwners
    two::Dict{Tuple{Int,NTuple{2,Int32}},_EmittedOwner}
    three::Dict{Tuple{Int,NTuple{3,Int32}},_EmittedOwner}
    four::Dict{Tuple{Int,NTuple{4,Int32}},_EmittedOwner}
    five::Dict{Tuple{Int,NTuple{5,Int32}},_EmittedOwner}
    six::Dict{Tuple{Int,NTuple{6,Int32}},_EmittedOwner}
    eight::Dict{Tuple{Int,NTuple{8,Int32}},_EmittedOwner}
    other::Dict{Tuple{Int,Tuple},_EmittedOwner}
end

_EmittedCellOwners()=_EmittedCellOwners(
    Dict{Tuple{Int,NTuple{2,Int32}},_EmittedOwner}(),Dict{Tuple{Int,NTuple{3,Int32}},_EmittedOwner}(),
    Dict{Tuple{Int,NTuple{4,Int32}},_EmittedOwner}(),Dict{Tuple{Int,NTuple{5,Int32}},_EmittedOwner}(),
    Dict{Tuple{Int,NTuple{6,Int32}},_EmittedOwner}(),Dict{Tuple{Int,NTuple{8,Int32}},_EmittedOwner}(),
    Dict{Tuple{Int,Tuple},_EmittedOwner}())

Base.length(owners::_EmittedCellOwners)=length(owners.two)+length(owners.three)+
    length(owners.four)+length(owners.five)+length(owners.six)+length(owners.eight)+
    length(owners.other)

# Every `(msh, rows)` key, untyped: for the rare whole-set comparison only.
function _emitted_owner_keys(owners::_EmittedCellOwners)
    result=Set{Tuple{Int,Tuple}}()
    for table in (owners.two,owners.three,owners.four,owners.five,owners.six,owners.eight,owners.other)
        for key in keys(table)
            push!(result,(key[1],Tuple(key[2])))
        end
    end
    return result
end

@inline _emitted_rows(cells::AbstractMatrix,column::Int,::Val{N}) where {N}=
    ntuple(slot->Int32(cells[slot,column]),Val(N))
@inline _emitted_rows(cells::AbstractMatrix,column::Int,rows,::Val{N}) where {N}=
    ntuple(slot->Int32(rows[cells[slot,column]]),Val(N))

# The owner of carrier cell `column` of `cells`, or `nothing`.
function _emitted_owner(owners::_EmittedCellOwners,msh::Int,cells::AbstractMatrix,column::Int)
    arity=size(cells,1)
    arity==4 && return get(owners.four,(msh,_emitted_rows(cells,column,Val(4))),nothing)
    arity==3 && return get(owners.three,(msh,_emitted_rows(cells,column,Val(3))),nothing)
    arity==2 && return get(owners.two,(msh,_emitted_rows(cells,column,Val(2))),nothing)
    arity==8 && return get(owners.eight,(msh,_emitted_rows(cells,column,Val(8))),nothing)
    arity==6 && return get(owners.six,(msh,_emitted_rows(cells,column,Val(6))),nothing)
    arity==5 && return get(owners.five,(msh,_emitted_rows(cells,column,Val(5))),nothing)
    return get(owners.other,(msh,Tuple(Int32(cells[slot,column]) for slot in 1:arity)),nothing)
end

# Record `owner` for part cell `column`, its rows translated through `rows`.
# `false` when that carrier cell already has an owner.
@inline function _emitted_owner_insert!(table::Dict{K,_EmittedOwner},key::K,owner::_EmittedOwner) where {K}
    before=length(table)
    get!(table,key,owner)
    return length(table)>before
end
function _emitted_owner_insert!(owners::_EmittedCellOwners,msh::Int,cells::AbstractMatrix,
                                column::Int,rows,owner::_EmittedOwner)
    arity=size(cells,1)
    arity==4 && return _emitted_owner_insert!(owners.four,(msh,_emitted_rows(cells,column,rows,Val(4))),owner)
    arity==3 && return _emitted_owner_insert!(owners.three,(msh,_emitted_rows(cells,column,rows,Val(3))),owner)
    arity==2 && return _emitted_owner_insert!(owners.two,(msh,_emitted_rows(cells,column,rows,Val(2))),owner)
    arity==8 && return _emitted_owner_insert!(owners.eight,(msh,_emitted_rows(cells,column,rows,Val(8))),owner)
    arity==6 && return _emitted_owner_insert!(owners.six,(msh,_emitted_rows(cells,column,rows,Val(6))),owner)
    arity==5 && return _emitted_owner_insert!(owners.five,(msh,_emitted_rows(cells,column,rows,Val(5))),owner)
    key=(msh,Tuple(Int32(rows[cells[slot,column]]) for slot in 1:arity))
    return _emitted_owner_insert!(owners.other,key,owner)
end

# Fingerprint of the native state that authenticates row and cell identity:
# two independent 64-bit mixes of every word, in order. Tables whose order is
# not part of their meaning (the entity table) add up the fingerprints of
# their entries. Values are mixed as `isequal` compares them: by value for
# integers, by bits for floating point.
struct _EmittedStateStamp
    nodes::Int
    cells::Int
    first::UInt64
    second::UInt64
end

const _EmittedMix=NTuple{2,UInt64}
const _EMITTED_MIX_SEED=(0xcbf29ce484222325,0x9e3779b97f4a7c15)
@inline function _emitted_mix(state::_EmittedMix,word::UInt64)
    second=(state[2]+word)*0xff51afd7ed558ccd
    return ((state[1]⊻word)*0x00000100000001b3,second⊻(second>>31))
end
_emitted_mix(state::_EmittedMix,value::Bool)=_emitted_mix(state,UInt64(value))
_emitted_mix(state::_EmittedMix,value::Integer)=_emitted_mix(state,value%UInt64)
_emitted_mix(state::_EmittedMix,value::Float64)=_emitted_mix(state,reinterpret(UInt64,value))
_emitted_mix(state::_EmittedMix,value::AbstractFloat)=_emitted_mix(state,Float64(value))
_emitted_mix(state::_EmittedMix,::Nothing)=_emitted_mix(state,0x6e6f7468696e6721)
_emitted_mix(state::_EmittedMix,value::Symbol)=_emitted_mix(state,String(value))
function _emitted_mix(state::_EmittedMix,value::AbstractString)
    state=_emitted_mix(state,UInt64(ncodeunits(value)))
    for unit in codeunits(value)
        state=_emitted_mix(state,UInt64(unit))
    end
    return state
end
function _emitted_mix(state::_EmittedMix,value::Union{Tuple,NamedTuple})
    state=_emitted_mix(state,UInt64(length(value)))
    for item in value
        state=_emitted_mix(state,item)
    end
    return state
end
function _emitted_mix(state::_EmittedMix,value::AbstractArray)
    state=_emitted_mix(state,UInt64(length(value)))
    for item in value
        state=_emitted_mix(state,item)
    end
    return state
end
function _emitted_mix(state::_EmittedMix,value::Array{T}) where {T<:Union{Int32,Int64,UInt64,Float64}}
    state=_emitted_mix(state,UInt64(length(value)))
    @inbounds for index in eachindex(value)
        state=_emitted_mix(state,value[index])
    end
    return state
end
function _emitted_mix(state::_EmittedMix,value::Union{AbstractDict,AbstractSet})
    first=UInt64(0);second=UInt64(0)
    for item in value
        entry=_emitted_mix(_EMITTED_MIX_SEED,item isa Pair ? (item.first,item.second) : item)
        first+=entry[1];second+=entry[2]
    end
    return _emitted_mix(_emitted_mix(_emitted_mix(state,UInt64(length(value))),first),second)
end
# The two record types of a mesh's native state, field by field.
function _emitted_mix(state::_EmittedMix,entity::MixedEntity)
    state=_emitted_mix(state,entity.dim);state=_emitted_mix(state,entity.tag)
    state=_emitted_mix(state,entity.bbox);state=_emitted_mix(state,entity.physical_tags)
    state=_emitted_mix(state,entity.boundaries);state=_emitted_mix(state,entity.embedded_curves)
    state=_emitted_mix(state,entity.partitions)
    return _emitted_mix(state,entity.parent)
end
function _emitted_mix(state::_EmittedMix,link::MixedPeriodicLink)
    state=_emitted_mix(state,link.dim);state=_emitted_mix(state,link.slave_entity)
    state=_emitted_mix(state,link.master_entity);state=_emitted_mix(state,link.affine)
    state=_emitted_mix(state,link.slave_nodes)
    return _emitted_mix(state,link.master_nodes)
end
# Any other value is its fields, in order.
function _emitted_mix(state::_EmittedMix,value::T) where {T}
    isstructtype(T) && fieldcount(T)>0 || throw(ArgumentError(
        "emitted receipt: no fingerprint for a value of type $T"))
    for field in 1:fieldcount(T)
        state=_emitted_mix(state,getfield(value,field))
    end
    return state
end

struct _ModelMeshEmittedReceipt
    actual::MixedMesh
    carrier::Union{Mesh,MixedMesh}
    parts::Dict{Tuple{Int,Int32},Tuple{Union{Mesh,MixedMesh},Vector{Int32}}}
    # The owner of each carrier cell, block by block in carrier order (the
    # tables keyed by rows are needed only to build the receipt).
    cell_owners::Vector{Vector{_EmittedOwner}}
    carrier_rows::Vector{Int32}
    native_state::_EmittedStateStamp
end

_emitted_projection_blocks(mesh::Mesh)=((1,mesh.segs,mesh.seg_tag),
    (2,mesh.tris,mesh.tri_tag),(4,mesh.tets,mesh.tet_tag))
_emitted_projection_blocks(mesh::MixedMesh)=
    ((Int(block.msh),block.nodes,block.tags) for block in mesh.blocks)

function _model_mesh_emitted_receipt(authority,parts,part_rows,caller)
    authority.mesh!==nothing && authority.carrier!==nothing || throw(ArgumentError(
        "$caller: emitted receipt requires an actual authority and its carrier"))
    length(parts)==length(part_rows) || throw(ArgumentError(
        "$caller: emitted constructor parts do not align with their row bindings"))
    bound=Dict{Tuple{Int,Int32},Tuple{Union{Mesh,MixedMesh},Vector{Int32}}}()
    owners=_EmittedCellOwners()
    for ((dimension,tag,part),rows) in zip(parts,part_rows)
        owner=(dimension,Int32(tag))
        haskey(authority.products,owner) && !haskey(bound,owner) || throw(ArgumentError(
            "$caller: emitted constructor part has absent or duplicate ownership"))
        length(rows)==nnodes(part) && allunique(rows) &&
            all(row->0<row<=nnodes(authority.carrier),rows) || throw(ArgumentError(
            "$caller: emitted constructor rows lack a unique carrier binding"))
        isequal(part.coords,authority.carrier.coords[:,rows]) || throw(ArgumentError(
            "$caller: emitted constructor coordinates disagree with their bound carrier"))
        part isa MixedMesh && part.entity_data!==nothing && throw(ArgumentError(
            "$caller: an emitted constructor part claims actual native metadata"))
        all(entry->isempty(entry[2]) || msh_dimension(entry[1])==dimension,
            _emitted_projection_blocks(part)) || throw(ArgumentError(
            "$caller: emitted constructor part contains a foreign cell dimension"))
        bound_rows=length(rows)
        for (msh,cells,tags) in _emitted_projection_blocks(part),column in axes(cells,2)
            for slot in axes(cells,1)
                0<cells[slot,column]<=bound_rows || throw(ArgumentError(
                    "$caller: emitted constructor cell references an absent bound vertex"))
            end
            _emitted_owner_insert!(owners,msh,cells,column,rows,owner) || throw(ArgumentError(
                "$caller: emitted constructors repeat a typed carrier cell"))
        end
        bound[owner]=(deepcopy(part),Int32.(rows))
    end
    # The emitted cells are distinct, so they partition the carrier when it
    # has as many cells and each of them is an emitted one. A carrier that
    # repeats a cell fails the count and takes the set comparison.
    carrier_count=0;covered=true
    for (msh,cells,tags) in _emitted_projection_blocks(authority.carrier)
        carrier_count+=size(cells,2)
        covered || continue
        for column in axes(cells,2)
            if _emitted_owner(owners,msh,cells,column)===nothing
                covered=false
                break
            end
        end
    end
    if !(covered && carrier_count==length(owners))
        carrier_keys=Set{Tuple{Int,Tuple}}((msh,Tuple(Int32.(cells[:,column]))) for (msh,cells,tags) in
            _emitted_projection_blocks(authority.carrier) for column in axes(cells,2))
        _emitted_owner_keys(owners)==carrier_keys || throw(ArgumentError(
            "$caller: emitted constructors do not partition the current carrier"))
    end
    carrier_owners=Vector{_EmittedOwner}[]
    for (msh,cells,tags) in _emitted_projection_blocks(authority.carrier)
        aligned=Vector{_EmittedOwner}(undef,size(cells,2))
        for column in axes(cells,2)
            found=_emitted_owner(owners,msh,cells,column)
            aligned[column]=found===nothing ? _EMITTED_NO_OWNER : found
        end
        push!(carrier_owners,aligned)
    end
    return _ModelMeshEmittedReceipt(authority.mesh,authority.carrier,bound,carrier_owners,
        copy(authority.carrier_rows),_model_mesh_emitted_state(authority.mesh))
end

# Fingerprint only native state that authenticates row/cell identity. Coordinate
# arithmetic of an emitted view never rewrites this actual state.
function _model_mesh_emitted_state(mesh)
    data=mesh.entity_data
    state=_EMITTED_MIX_SEED
    state=_emitted_mix(state,mesh.coords)
    state=_emitted_mix(state,data.node_entities)
    state=_emitted_mix(state,data.node_parametric)
    state=_emitted_mix(state,data.external_node_tags)
    state=_emitted_mix(state,data.entities)
    state=_emitted_mix(state,mesh.periodic_links)
    cells=0
    state=_emitted_mix(state,UInt64(length(mesh.blocks)))
    for (block,owners,tags) in zip(mesh.blocks,data.block_entities,data.external_element_tags)
        cells+=size(block.nodes,2)
        state=_emitted_mix(state,Int(block.msh))
        state=_emitted_mix(state,UInt64(size(block.nodes,1)))
        state=_emitted_mix(state,block.nodes)
        state=_emitted_mix(state,block.tags)
        state=_emitted_mix(state,owners)
        state=_emitted_mix(state,tags)
    end
    return _EmittedStateStamp(size(mesh.coords,2),cells,state[1],state[2])
end

_model_mesh_emitted_state_matches(mesh,state::_EmittedStateStamp)=
    _model_mesh_emitted_state(mesh)==state

function _model_mesh_with_emitted_receipt(authority,parts,rows;caller="emitted geometry")
    receipt=_model_mesh_emitted_receipt(authority,parts,rows,caller)
    return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
        authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,
        authority.producer_rows,authority.carrier_cells,authority.surface_constructor_inputs,receipt)
end

# These branches still require the current generation's exact constructor
# receipt and literal cell incidence below. Other naked inputs retain their
# existing geometric certification path.
function _model_mesh_emitted_scope(m,tag)
    params=_extrude_gate(m,3,tag)
    if params!==nothing
        return params.recombine && params.quad_to_tri===:none
    end
    haskey(m.volumes,tag) && !isempty(m.volumes[tag]) &&
        !haskey(m.meshing.transfinite_volumes,tag) && !_implicit_volume_surface(m,tag)
end

function _model_mesh_emitted_projection(m,mesh,dimension,tag,receipt,caller)
    receipt isa _ModelMeshEmittedReceipt || return nothing
    dimension==3 || return nothing
    _model_mesh_emitted_scope(m,tag) || return nothing
    mesh isa MixedMesh && mesh.entity_data!==nothing && return nothing
    authority=m.meshing.mesh_authority
    authority!==nothing && authority.mesh===receipt.actual &&
        authority.carrier===receipt.carrier || return nothing
    owner=(dimension,Int32(tag))
    constructor=get(receipt.parts,owner,nothing)
    constructor===nothing && return nothing
    part,constructor_rows=constructor
    # Exact constructor equality authenticates the detached naked part.
    # Coordinates verify those row events; they do not supply native identity.
    nnodes(mesh)==nnodes(part) || return nothing
    supplied=collect(_emitted_projection_blocks(mesh))
    expected=collect(_emitted_projection_blocks(part))
    length(supplied)==length(expected) && all(isequal(a,b) for (a,b) in zip(supplied,expected)) ||
        return nothing
    isequal(mesh.coords,part.coords) || throw(ArgumentError(
        "$caller: emitted constructor coordinates no longer match their selected row receipt"))
    isequal(mesh.coords,authority.carrier.coords[:,constructor_rows]) || throw(ArgumentError(
        "$caller: emitted constructor receipt no longer matches its carrier coordinates"))
    authority.carrier_rows==receipt.carrier_rows &&
        _model_mesh_emitted_state_matches(authority.mesh,receipt.native_state) || throw(ArgumentError(
        "$caller: emitted receipt no longer identifies its unchanged actual constructor state"))
    for (owned_part,rows) in values(receipt.parts)
        isequal(owned_part.coords,authority.carrier.coords[:,rows]) || throw(ArgumentError(
            "$caller: emitted closure coordinates no longer match their constructor row receipt"))
    end
    actual=authority.mesh;data=actual.entity_data
    closure=_model_mesh_projection_closure(m,dimension,tag,caller)
    actual_rows=authority.carrier_rows[constructor_rows]
    allunique(actual_rows) && Set(actual_rows)==Set(authority.products[owner].rows) ||
        throw(ArgumentError("$caller: emitted part does not bind the complete selected actual entity"))
    # Verify actual cell UID, owner, family and ordered primary incidence for
    # every usable carrier cell. Proxy/subdivision cells need their own planner.
    native=Dict{UInt64,Tuple{Int,Tuple{Int,Int32},Vector{Int32}}}()
    for (block,owners,uids) in zip(actual.blocks,data.block_entities,data.external_element_tags)
        for column in axes(block.nodes,2)
            uid=uids[column]
            haskey(native,uid) && throw(ArgumentError("$caller: duplicate actual cell UID"))
            native[uid]=(Int(block.msh),(msh_dimension(block.msh),owners[column]),copy(block.nodes[:,column]))
        end
    end
    carrier_cells=collect(_emitted_projection_blocks(authority.carrier))
    length(carrier_cells)==length(authority.carrier_cells) || throw(ArgumentError(
        "$caller: emitted carrier cell blocks lack their actual constructor events"))
    events=Dict{Tuple{Int,Tuple},UInt64}();seen=Set{UInt64}()
    length(carrier_cells)==length(receipt.cell_owners) || throw(ArgumentError(
        "$caller: emitted receipt does not cover the carrier's cell blocks"))
    for (bi,((msh,cells,tags),uids)) in enumerate(zip(carrier_cells,authority.carrier_cells))
        length(uids)==size(cells,2)==length(receipt.cell_owners[bi]) || throw(ArgumentError(
            "$caller: emitted carrier cell columns lack their actual constructor events"))
        for column in axes(cells,2)
            cell_owner=receipt.cell_owners[bi][column]
            cell_owner!=_EMITTED_NO_OWNER || throw(ArgumentError("$caller: carrier cell lacks its emitted constructor owner"))
            cell_owner in closure || continue
            uid=uids[column]
            uid!==nothing || throw(ArgumentError("$caller: emitted closure cell lacks its actual constructor UID"))
            current=get(native,uid,nothing)
            current!==nothing && current[1]==msh && current[2]==cell_owner &&
                current[3]==authority.carrier_rows[cells[:,column]] || throw(ArgumentError(
                "$caller: emitted carrier cell loses its typed actual UID incidence"))
            uid in seen && throw(ArgumentError("$caller: emitted carrier repeats an actual cell event"))
            push!(seen,uid)
            key=(msh,Tuple(cells[:,column]))
            haskey(events,key) && throw(ArgumentError("$caller: emitted carrier repeats a typed constructor cell"))
            events[key]=uid
        end
    end
    selected=Set{UInt64}()
    for (msh,cells,tags) in supplied,column in axes(cells,2)
        uid=get(events,(msh,Tuple(constructor_rows[cells[:,column]])),nothing)
        uid!==nothing && native[uid][2]==owner || throw(ArgumentError(
            "$caller: emitted selected cell lacks its unique owned actual constructor event"))
        push!(selected,uid)
    end
    selected==Set(uid for (uid,cell) in native if cell[2]==owner) || throw(ArgumentError(
        "$caller: emitted selected cells do not cover the current actual entity"))
    dimension==3 || throw(ArgumentError("$caller: this emitted receipt adapter requires the Volume planner"))
    _volume_projection_contract(mesh,caller)
    any(link->(link.dim,link.slave_entity) in closure || (link.dim,link.master_entity) in closure,
        actual.periodic_links) && throw(ArgumentError(
        "$caller: emitted periodic closure requires its affine constructor planner"))
    inverse=Dict(row=>Int32(position) for (position,row) in enumerate(authority.carrier_rows))
    rows=copy(actual_rows);present=Set(rows)
    needed=Set{Int32}()
    for selected_owner in closure
        product=get(authority.products,selected_owner,nothing)
        product===nothing || union!(needed,product.rows)
    end
    for row in sort!(collect(setdiff(needed,present));by=row->authority.private_rank[row])
        haskey(inverse,row) || throw(ArgumentError(
            "$caller: emitted closure loses an actual owned vertex row"))
        push!(rows,row)
    end
    row_map=Dict(row=>Int32(position) for (position,row) in enumerate(rows))
    coords=authority.carrier.coords[:,Int32[inverse[row] for row in rows]]
    blocks=ElementBlock[];owners=Vector{Int32}[]
    for (block,entities,uids) in zip(actual.blocks,data.block_entities,data.external_element_tags)
        columns=Int[column for column in axes(block.nodes,2) if
            msh_dimension(block.msh)<dimension && (msh_dimension(block.msh),entities[column]) in closure]
        isempty(columns) && continue
        for column in columns
            # Constructor-UID attestation is only possible for cells whose
            # owner entity was itself emitted as a part: emitted parts
            # partition the carrier, so their lower cells appear in `seen`.
            # Un-emitted lower entities (and Point cells, which are never
            # carrier cells) are certified by the bound-row, typed-boundary
            # and edge-complex checks below instead.
            cell_owner=(msh_dimension(block.msh),entities[column])
            (cell_owner[1]==0 || !haskey(receipt.parts,cell_owner) ||
                uids[column] in seen) || throw(ArgumentError(
                "$caller: emitted lower cell lacks its literal actual constructor event"))
        end
        nodes=map(row->get(row_map,row,Int32(0)),block.nodes[:,columns])
        all(>(0),nodes) || throw(ArgumentError("$caller: emitted lower cell references an absent bound row"))
        physical=Int32[_model_projection_legacy_tag(_model_projection_physical_tags(m,
            msh_dimension(block.msh),Int(entities[column]))) for column in columns]
        push!(blocks,ElementBlock(block.msh,nodes,physical));push!(owners,entities[columns])
    end
    physical=_model_projection_legacy_tag(_model_projection_physical_tags(m,dimension,tag))
    for (msh,cells,tags) in supplied
        isempty(cells) && continue
        push!(blocks,ElementBlock(msh,cells,fill(physical,size(cells,2))))
        push!(owners,fill(Int32(tag),size(cells,2)))
    end
    # The projection still requires the complete actual typed boundary and
    # every literal lower Line in the emitted Volume's edge complex.
    boundary_surfaces=Set(Int32.(abs.(_model_volume_boundary_surfaces(m,tag,caller))))
    lower_faces=Set{Tuple}();line_edges=Dict{Tuple,Int32}()
    for (block,block_owners) in zip(blocks,owners),column in axes(block.nodes,2)
        cell_dimension=msh_dimension(block.msh)
        if cell_dimension==2 && block_owners[column] in boundary_surfaces
            face=Tuple(sort!(collect(block.nodes[:,column])))
            face in lower_faces && throw(ArgumentError("$caller: emitted lower Surfaces repeat a typed boundary facet"))
            push!(lower_faces,face)
        elseif cell_dimension==1
            edge=Tuple(sort!(collect(block.nodes[:,column])))
            previous=get(line_edges,edge,nothing)
            previous===nothing || previous==block_owners[column] || throw(ArgumentError(
                "$caller: emitted lower Line incidence aliases distinct Curve owners"))
            line_edges[edge]=block_owners[column]
        end
    end
    lower_faces==Set(Tuple(face) for face in first(_volume_cell_boundary_faces(mesh))) ||
        throw(ArgumentError("$caller: emitted lower Surfaces do not cover its exact typed Volume boundary"))
    volume_edges=_volume_edge_set(mesh)
    all(edge->edge in volume_edges,keys(line_edges)) || throw(ArgumentError(
        "$caller: emitted lower Line is absent from its Volume edge complex"))
    # Synthetic projection labels preserve the geometric frontend contract;
    # actual UIDs served solely to attest the constructor and lower incidence.
    entities=Dict{Tuple{Int,Int},MixedEntity}();names=Dict{Tuple{Int,Int},String}()
    for selected_owner in closure
        current=get(data.entities,selected_owner,nothing)
        current!==nothing || throw(ArgumentError("$caller: emitted closure lacks entity metadata"))
        entity_nodes=Set{Int32}(Int32(position) for (position,row) in enumerate(rows)
            if data.node_entities[row]==selected_owner)
        for (block,block_owners) in zip(blocks,owners),column in axes(block.nodes,2)
            (msh_dimension(block.msh),block_owners[column])==selected_owner || continue
            union!(entity_nodes,block.nodes[:,column])
        end
        tags=_model_projection_physical_tags(m,selected_owner[1],Int(selected_owner[2]))
        entities[selected_owner]=_model_mesh_entity_metadata(m,selected_owner[1],Int(selected_owner[2]),
            _model_projection_bbox((;coords),entity_nodes,caller);physical_tags=tags,
            boundaries=copy(current.boundaries),embedded_curves=copy(current.embedded_curves))
        for physical_tag in tags
            key=(selected_owner[1],Int(physical_tag))
            haskey(m.physical_names,key) && (names[key]=m.physical_names[key])
        end
    end
    metadata=MixedEntityData(entities;
        node_entities=data.node_entities[rows],node_parametric=Union{Nothing,Vector{Float64}}[nothing for _ in rows],
        external_node_tags=UInt64.(1:length(rows)),block_entities=owners,
        external_element_tags=_model_projection_external_elements(blocks))
    output=MixedMesh(coords,blocks;physical_names=names,entity_data=metadata,elementary_entities=owners)
    diagnostic=validate(output)
    diagnostic.ok || throw(ArgumentError("$caller: invalid emitted projection — "*join(diagnostic.messages,"; ")))
    return output
end
