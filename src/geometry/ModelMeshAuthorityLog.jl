# Incremental mesh authority.
#
# A Surface pass completes one entity at a time, and each completed product has
# to be visible to the next producer. Rebuilding the whole authority per entity
# costs O(nodes + cells + entities) every time, which is quadratic over the
# pass. The log keeps the same data in append-only storage with spare capacity,
# together with the per-entity indexes `_model_mesh_authority` derives, and
# hands out ordinary `_ModelMeshAuthority` values whose arrays are prefixes of
# that storage. An append costs the rows and cells it adds plus the products of
# the entities it touches.
#
# Contract:
# - Rows and cells are only appended. A caller that edits an existing row
#   discards the log and seeds a new one from a full rebuild.
# - An append the full constructors would reject returns `false` and leaves the
#   log invalid. The caller then runs the full rebuild, which reports the error.
# - Authorities taken from one log share its entity tables (`products`, `done`,
#   `producer_rows`): an older one also sees the entities completed after it.
#   Nothing may keep one across an append and rely on its entity set.
# - The arrays of such an authority alias the log's storage. They are
#   read-only, like the arrays of every authority.
using ..Elements: MshAncillarySection,MshDataSection

# The first `dims` entries of `memory` as an array, without copying where the
# runtime can alias a prefix of a `Memory`.
@static if isdefined(Base,:wrap)
    @inline _model_mesh_prefix(memory::Memory{T},dims::NTuple{N,Int}) where {T,N}=
        Base.wrap(Array,memory,dims)
else
    function _model_mesh_prefix(memory::Memory{T},dims::NTuple{N,Int}) where {T,N}
        output=Array{T,N}(undef,dims)
        count=prod(dims)
        count>0 && copyto!(output,1,memory,1,count)
        return output
    end
end

function _model_mesh_reserve(source::AbstractArray{T},used::Int,capacity::Int) where T
    memory=Memory{T}(undef,capacity)
    used>0 && copyto!(memory,1,source,1,used)
    return memory
end

function _model_mesh_row_numbers(capacity::Int)
    rows=Memory{Int32}(undef,capacity)
    for row in 1:capacity
        rows[row]=Int32(row)
    end
    return rows
end

mutable struct _ModelMeshLog
    valid::Bool
    caller::String
    nodes::Int
    coords::Memory{Float64}
    node_entities::Memory{Tuple{Int,Int32}}
    node_parametric::Memory{Union{Nothing,Vector{Float64}}}
    node_tags::Memory{UInt64}
    rank::Memory{UInt64}
    carrier_rows::Memory{Int32}
    rank_max::UInt64
    node_rows::Dict{UInt64,Int32}
    cells::Int
    types::Vector{Int}
    counts::Vector{Int}
    cell_nodes::Vector{Memory{Int32}}
    cell_physical::Vector{Memory{Int32}}
    cell_owners::Vector{Memory{Int32}}
    cell_elementary::Vector{Memory{Int32}}
    cell_tags::Vector{Memory{UInt64}}
    cell_bindings::Vector{Memory{Union{Nothing,UInt64}}}
    block_index::Dict{Int,Int}
    element_tags::Set{UInt64}
    entities::Dict{Tuple{Int,Int},MixedEntity}
    # Rows each entity owns, and its cells as (block, columns) in block order:
    # the two indexes `_model_mesh_authority` builds from a complete mesh.
    owned::Dict{Tuple{Int,Int32},Vector{Int32}}
    selected::Dict{Tuple{Int,Int32},Vector{Tuple{Int,Vector{Int}}}}
    touched::Vector{Tuple{Int,Int32}}
    touched_set::Set{Tuple{Int,Int32}}
    products::Dict{Tuple{Int,Int32},_ModelMeshProduct}
    done::Set{Tuple{Int,Int32}}
    native_inputs::Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}
    producer_rows::Dict{Tuple{Int,Int32},Vector{Int32}}
end

# Counts of full seeds and incremental appends, for tests of the scaling
# contract: a pass over n independent Surfaces seeds once and appends n times.
const _MODEL_MESH_LOG_SEEDS=Ref(0)
const _MODEL_MESH_LOG_APPENDS=Ref(0)

# Seed a log from a complete authority, or return `nothing` when the authority
# holds anything the log cannot extend: a separate geometry carrier, native
# inputs, mesh sections, or a product that its mesh does not imply (an empty
# Surface product).
function _model_mesh_log(authority::_ModelMeshAuthority;caller="mesh authority")
    mesh=authority.mesh
    mesh===nothing && return nothing
    data=mesh.entity_data
    data!==nothing && data.native_mesh_tags || return nothing
    elementary=mesh.elementary_entities
    (authority.carrier===mesh && elementary!==nothing &&
     isempty(authority.native_inputs) && authority.emitted_geometry===nothing &&
     isempty(authority.surface_constructor_inputs) &&
     isempty(mesh.physical_names) && isempty(mesh.periodic_links) &&
     isempty(mesh.ancillary_sections) && isempty(mesh.data_sections) &&
     mesh.partition_data===nothing &&
     all(block->block isa ElementBlock,mesh.blocks) &&
     allunique(block.msh for block in mesh.blocks)) || return nothing
    count=nnodes(mesh)
    (length(authority.private_rank)==count && length(authority.carrier_rows)==count &&
     length(authority.carrier_cells)==length(mesh.blocks)==length(elementary)) || return nothing
    for row in 1:count
        authority.carrier_rows[row]==row || return nothing
    end
    owned=Dict{Tuple{Int,Int32},Vector{Int32}}()
    selected=Dict{Tuple{Int,Int32},Vector{Tuple{Int,Vector{Int}}}}()
    capacity=min(Int(typemax(Int32)),count+max(count,4096))
    node_rows=Dict{UInt64,Int32}();sizehint!(node_rows,capacity)
    for (row,owner) in enumerate(data.node_entities)
        push!(get!(()->Int32[],owned,owner),Int32(row))
        get!(()->Tuple{Int,Vector{Int}}[],selected,owner)
        node_rows[data.external_node_tags[row]]=Int32(row)
    end
    length(node_rows)==count || return nothing
    types=Int[];counts=Int[];block_index=Dict{Int,Int}()
    cell_nodes=Memory{Int32}[];cell_physical=Memory{Int32}[]
    cell_owners=Memory{Int32}[];cell_elementary=Memory{Int32}[]
    cell_tags=Memory{UInt64}[];cell_bindings=Memory{Union{Nothing,UInt64}}[]
    element_tags=Set{UInt64}();cells=0
    for (bi,block) in enumerate(mesh.blocks)
        owners=data.block_entities[bi];tags=data.external_element_tags[bi]
        columns=size(block.nodes,2)
        (length(owners)==columns==length(tags)==length(elementary[bi]) &&
         isequal(authority.carrier_cells[bi],tags) && elementary[bi]==owners) || return nothing
        dimension=msh_dimension(block.msh)
        for column in 1:columns
            tags[column] in element_tags && return nothing
            push!(element_tags,tags[column])
            selection=get!(()->Tuple{Int,Vector{Int}}[],selected,(dimension,owners[column]))
            if isempty(selection) || last(selection)[1]!=bi
                push!(selection,(bi,Int[column]))
            else
                push!(last(selection)[2],column)
            end
        end
        reserve=columns+max(columns,1024)
        push!(types,block.msh);push!(counts,columns);block_index[block.msh]=bi
        push!(cell_nodes,_model_mesh_reserve(block.nodes,length(block.nodes),size(block.nodes,1)*reserve))
        push!(cell_physical,_model_mesh_reserve(block.tags,columns,reserve))
        push!(cell_owners,_model_mesh_reserve(owners,columns,reserve))
        push!(cell_elementary,_model_mesh_reserve(elementary[bi],columns,reserve))
        push!(cell_tags,_model_mesh_reserve(tags,columns,reserve))
        push!(cell_bindings,_model_mesh_reserve(authority.carrier_cells[bi],columns,reserve))
        cells+=columns
    end
    # Every product must follow from the mesh, so that an append reproduces
    # what a full rebuild of the extended mesh would hold.
    (length(authority.products)==length(selected) &&
     all(owner->haskey(selected,owner),keys(authority.products)) &&
     all(owner->haskey(authority.products,owner),keys(authority.producer_rows))) || return nothing
    _MODEL_MESH_LOG_SEEDS[]+=1
    return _ModelMeshLog(true,String(caller),count,
        _model_mesh_reserve(mesh.coords,3count,3capacity),
        _model_mesh_reserve(data.node_entities,count,capacity),
        _model_mesh_reserve(data.node_parametric,count,capacity),
        _model_mesh_reserve(data.external_node_tags,count,capacity),
        _model_mesh_reserve(authority.private_rank,count,capacity),
        _model_mesh_row_numbers(capacity),
        maximum(authority.private_rank;init=UInt64(0)),node_rows,cells,
        types,counts,cell_nodes,cell_physical,cell_owners,cell_elementary,cell_tags,
        cell_bindings,block_index,element_tags,data.entities,owned,selected,
        Tuple{Int,Int32}[],Set{Tuple{Int,Int32}}(),copy(authority.products),
        copy(authority.done),Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}(),
        copy(authority.producer_rows))
end

function _model_mesh_log_grow_nodes!(log::_ModelMeshLog,needed::Int)
    used=log.nodes
    capacity=min(Int(typemax(Int32)),max(needed,2*length(log.node_tags)))
    log.coords=_model_mesh_reserve(log.coords,3used,3capacity)
    log.node_entities=_model_mesh_reserve(log.node_entities,used,capacity)
    log.node_parametric=_model_mesh_reserve(log.node_parametric,used,capacity)
    log.node_tags=_model_mesh_reserve(log.node_tags,used,capacity)
    log.rank=_model_mesh_reserve(log.rank,used,capacity)
    log.carrier_rows=_model_mesh_row_numbers(capacity)
    return nothing
end

function _model_mesh_log_grow_cells!(log::_ModelMeshLog,bi::Int,needed::Int)
    used=log.counts[bi]
    width=length(log.cell_nodes[bi])÷length(log.cell_tags[bi])
    capacity=max(needed,2*length(log.cell_tags[bi]))
    log.cell_nodes[bi]=_model_mesh_reserve(log.cell_nodes[bi],width*used,width*capacity)
    log.cell_physical[bi]=_model_mesh_reserve(log.cell_physical[bi],used,capacity)
    log.cell_owners[bi]=_model_mesh_reserve(log.cell_owners[bi],used,capacity)
    log.cell_elementary[bi]=_model_mesh_reserve(log.cell_elementary[bi],used,capacity)
    log.cell_tags[bi]=_model_mesh_reserve(log.cell_tags[bi],used,capacity)
    log.cell_bindings[bi]=_model_mesh_reserve(log.cell_bindings[bi],used,capacity)
    return nothing
end

@inline function _model_mesh_log_touch!(log::_ModelMeshLog,owner::Tuple{Int,Int32})
    (!isempty(log.touched) && log.touched[end]==owner) && return nothing
    owner in log.touched_set && return nothing
    push!(log.touched_set,owner);push!(log.touched,owner)
    return nothing
end

# Append one node row. `parameter` is the producer's full parameter vector; the
# stored value is its projection onto the owner's dimension. `false` marks a
# row the full constructors would reject (or one whose rank does not extend
# the private order), and invalidates the log.
function _model_mesh_log_node!(log::_ModelMeshLog,point::NTuple{3,Float64},
                               owner::Tuple{Int,Int32},parameter,tag::UInt64,rank::UInt64)
    log.valid || return false
    log.valid=false
    dimension=owner[1]
    (isfinite(point[1]) && isfinite(point[2]) && isfinite(point[3])) || return false
    0<=dimension<=3 || return false
    (owner[2]>0 || (owner[2]==0 && haskey(log.entities,(dimension,0)))) || return false
    projected=_model_mesh_projected_parameter(owner,parameter)
    projected===nothing || all(isfinite,projected) || return false
    rank>log.rank_max || return false
    haskey(log.node_rows,tag) && return false
    row=log.nodes+1
    row<=typemax(Int32) || return false
    row>length(log.node_tags) && _model_mesh_log_grow_nodes!(log,row)
    base=3*(row-1)
    log.coords[base+1]=point[1];log.coords[base+2]=point[2];log.coords[base+3]=point[3]
    log.node_entities[row]=owner
    log.node_parametric[row]=projected
    log.node_tags[row]=tag
    log.rank[row]=rank
    log.rank_max=rank
    log.node_rows[tag]=Int32(row)
    log.nodes=row
    push!(get!(()->Int32[],log.owned,owner),Int32(row))
    haskey(log.selected,owner) || (log.selected[owner]=Tuple{Int,Vector{Int}}[])
    _model_mesh_log_touch!(log,owner)
    log.valid=true
    return true
end

# Append one cell to the block of its MSH type; a first cell of a new type
# opens a block after the existing ones, as the full grouping does.
function _model_mesh_log_cell!(log::_ModelMeshLog,msh::Int,owner::Int32,
                               nodes::Vector{Int32},physical::Int32,tag::UInt64)
    log.valid || return false
    log.valid=false
    bi=get(log.block_index,msh,0)
    spec=get(Elements._MSH_CATALOG_LOOKUP,msh,nothing)
    spec===nothing && return false
    width=spec.nnodes;dimension=spec.dim
    (width>0 && length(nodes)==width) || return false
    (owner>0 || (owner==0 && haskey(log.entities,(dimension,0)))) || return false
    key=(dimension,owner)
    entity=get(log.entities,(dimension,Int(owner)),nothing)
    entity===nothing && !haskey(log.owned,key) && return false
    projected=entity===nothing || isempty(entity.physical_tags) ? Int32(0) :
        first(entity.physical_tags)
    physical==projected || return false
    for node in nodes
        0<node<=log.nodes || return false
    end
    tag in log.element_tags && return false
    log.cells<typemax(Int32) || return false
    if bi==0
        reserve=1024
        push!(log.types,msh);push!(log.counts,0)
        push!(log.cell_nodes,Memory{Int32}(undef,width*reserve))
        push!(log.cell_physical,Memory{Int32}(undef,reserve))
        push!(log.cell_owners,Memory{Int32}(undef,reserve))
        push!(log.cell_elementary,Memory{Int32}(undef,reserve))
        push!(log.cell_tags,Memory{UInt64}(undef,reserve))
        push!(log.cell_bindings,Memory{Union{Nothing,UInt64}}(undef,reserve))
        bi=length(log.types);log.block_index[msh]=bi
    end
    column=log.counts[bi]+1
    column>length(log.cell_tags[bi]) && _model_mesh_log_grow_cells!(log,bi,column)
    storage=log.cell_nodes[bi];base=width*(column-1)
    for slot in 1:width
        storage[base+slot]=nodes[slot]
    end
    log.cell_physical[bi][column]=physical
    log.cell_owners[bi][column]=owner
    log.cell_elementary[bi][column]=owner
    log.cell_tags[bi][column]=tag
    log.cell_bindings[bi][column]=tag
    log.counts[bi]=column
    log.cells+=1
    push!(log.element_tags,tag)
    selection=get!(()->Tuple{Int,Vector{Int}}[],log.selected,key)
    position=findfirst(entry->entry[1]>=bi,selection)
    if position===nothing
        push!(selection,(bi,Int[column]))
    elseif selection[position][1]==bi
        push!(selection[position][2],column)
    else
        insert!(selection,position,(bi,Int[column]))
    end
    _model_mesh_log_touch!(log,key)
    log.valid=true
    return true
end

function _model_mesh_log_mesh(log::_ModelMeshLog)
    count=log.nodes
    blocks=MixedElementBlock[];owners=Vector{Int32}[];elementary=Vector{Int32}[]
    tags=Vector{UInt64}[]
    for (bi,msh) in enumerate(log.types)
        columns=log.counts[bi]
        width=length(log.cell_nodes[bi])÷length(log.cell_tags[bi])
        push!(blocks,ElementBlock(Elements._TRUSTED_MIXED_STORAGE,msh,
            _model_mesh_prefix(log.cell_nodes[bi],(width,columns)),
            _model_mesh_prefix(log.cell_physical[bi],(columns,))))
        push!(owners,_model_mesh_prefix(log.cell_owners[bi],(columns,)))
        push!(elementary,_model_mesh_prefix(log.cell_elementary[bi],(columns,)))
        push!(tags,_model_mesh_prefix(log.cell_tags[bi],(columns,)))
    end
    data=MixedEntityData(Elements._OWNED_MIXED_ENTITY_DATA,log.entities,
        _model_mesh_prefix(log.node_entities,(count,)),
        _model_mesh_prefix(log.node_parametric,(count,)),
        _model_mesh_prefix(log.node_tags,(count,)),owners,tags,true)
    return MixedMesh(Elements._TRUSTED_MIXED_STORAGE,
        _model_mesh_prefix(log.coords,(3,count)),blocks,Dict{Tuple{Int,Int},String}(),
        data,elementary,MixedPeriodicLink[],MshAncillarySection[],MshDataSection[],nothing)
end

# The authority of everything appended so far. `done` is the pass's set of
# completed entities; `completed` names the entity this append completed, so
# the log's own set follows without a rescan while the sizes agree.
function _model_mesh_log_authority!(log::_ModelMeshLog,done,completed=nothing)
    log.valid || return nothing
    if completed!==nothing && length(done)==length(log.done)+1 &&
            completed in done && !(completed in log.done)
        push!(log.done,completed)
    elseif length(done)!=length(log.done) || (completed!==nothing && !(completed in log.done))
        empty!(log.done)
        for (dimension,tag) in done
            push!(log.done,(Int(dimension),Int32(tag)))
        end
    end
    mesh=_model_mesh_log_mesh(log)
    rank=_model_mesh_prefix(log.rank,(log.nodes,))
    for owner in log.touched
        log.products[owner]=_model_mesh_owner_product(mesh,owner,
            get(log.owned,owner,Int32[]),log.selected[owner],rank,log.caller)
    end
    empty!(log.touched);empty!(log.touched_set)
    bindings=Vector{Union{Nothing,UInt64}}[
        _model_mesh_prefix(log.cell_bindings[bi],(log.counts[bi],)) for bi in eachindex(log.types)]
    return _ModelMeshAuthority(mesh,log.products,log.done,rank,log.native_inputs,mesh,
        _model_mesh_prefix(log.carrier_rows,(log.nodes,)),log.producer_rows,bindings,
        Dict{Tuple{Int,Int32},MixedMesh}())
end

# Bind one entity's producer rows, with the checks `_model_mesh_with_producer_rows`
# applies. `false` leaves the table unchanged.
function _model_mesh_log_producer_rows!(log::_ModelMeshLog,owner::Tuple{Int,Int32},rows::Vector{Int32})
    log.valid && haskey(log.products,owner) || return false
    all(row->0<row<=log.nodes,rows) && allunique(rows) || return false
    log.producer_rows[owner]=rows
    return true
end

# Structural equality of two authorities or any of their parts, for the
# verification mode: the incremental authority against a full rebuild.
function _model_mesh_same(a,b)
    a===b && return true
    typeof(a)===typeof(b) || return false
    if a isa AbstractArray
        size(a)==size(b) || return false
        for index in eachindex(a)
            _model_mesh_same(a[index],b[index]) || return false
        end
        return true
    elseif a isa AbstractDict
        length(a)==length(b) || return false
        for (key,value) in a
            haskey(b,key) && _model_mesh_same(value,b[key]) || return false
        end
        return true
    elseif a isa AbstractSet
        return a==b
    elseif isstructtype(typeof(a)) && fieldcount(typeof(a))>0 && !(a isa Tuple) &&
            !(a isa AbstractString) && !(a isa Number)
        for name in fieldnames(typeof(a))
            _model_mesh_same(getfield(a,name),getfield(b,name)) || return false
        end
        return true
    end
    return isequal(a,b)
end

# Set by tests: every incremental authority is then also rebuilt in full and
# compared, and a difference is an error.
const _MODEL_MESH_LOG_VERIFY=Ref(false)
