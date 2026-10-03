# Native mixed cache support. Dense element tags follow the stored blocks;
# classification and queries retain each block's actual MSH family.
import ..Model
import ..Elements

include("APIMixedP2Certification.jl")

function _canonical_mixed_copy(mesh::MixedMesh)
    order=Int[]
    groups=Dict{Int,Vector{Int}}()
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock || throw(ArgumentError(
            "API: native mixed caches require fixed-connectivity elements"))
        haskey(groups,block.msh) || push!(order,block.msh)
        push!(get!(()->Int[],groups,block.msh),bi)
    end
    blocks=ElementBlock[]
    elementary=mesh.elementary_entities===nothing ? nothing : Vector{Int32}[]
    block_owners=Vector{Int32}[];external_tags=Vector{UInt64}[]
    for msh in order
        indices=groups[msh]
        nodes=hcat((mesh.blocks[i].nodes for i in indices)...)
        tags=vcat((mesh.blocks[i].tags for i in indices)...)
        push!(blocks,ElementBlock(msh,nodes,tags))
        elementary===nothing || push!(elementary,
            vcat((mesh.elementary_entities[i] for i in indices)...))
        if mesh.entity_data!==nothing
            push!(block_owners,vcat((mesh.entity_data.block_entities[i] for i in indices)...))
            push!(external_tags,vcat((mesh.entity_data.external_element_tags[i] for i in indices)...))
        end
    end
    data=mesh.entity_data
    if data!==nothing
        data=Elements.MixedEntityData(data.entities;
            node_entities=data.node_entities,node_parametric=data.node_parametric,
            external_node_tags=data.external_node_tags,
            block_entities=block_owners,external_element_tags=external_tags)
    end
    return MixedMesh(mesh.coords,blocks;physical_names=mesh.physical_names,
        entity_data=data,elementary_entities=elementary,
        periodic_links=mesh.periodic_links,
        ancillary_sections=mesh.ancillary_sections,data_sections=mesh.data_sections,
        partition_data=mesh.partition_data)
end

function _cache_native_blocks(mesh::Mesh)
    return ((1,mesh.segs,mesh.seg_tag),(2,mesh.tris,mesh.tri_tag),
            (4,mesh.tets,mesh.tet_tag))
end

function _cache_native_blocks(mesh::MixedMesh)
    blocks=Tuple{Int,Matrix{Int32},Vector{Int32}}[]
    for block in mesh.blocks
        block isa ElementBlock || throw(ArgumentError(
            "API: generated mixed caches require fixed-connectivity elements"))
        push!(blocks,(block.msh,block.nodes,block.tags))
    end
    return blocks
end

function _cache_catalog(mesh,class=nothing)
    result=Tuple{Int32,Int,Int,Matrix{Int32},Vector{Int32}}[]
    offset=0
    family_offsets=Dict{Int,Int}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        owners=class===nothing ? Int32[] : get(class.cell_entities,msh,Int32[])
        family_offset=get(family_offsets,msh,0)
        if class!==nothing && (family_offset!=0 || length(owners)!=size(cells,2))
            owners=owners[family_offset+1:family_offset+size(cells,2)]
        end
        push!(result,(Int32(msh),msh_dimension(msh),offset,cells,owners))
        offset=Base.checked_add(offset,size(cells,2))
        family_offsets[msh]=family_offset+size(cells,2)
    end
    return result
end

function _mesh_element_offsets(mesh::MixedMesh)
    total=0
    for (_,cells,_) in _cache_native_blocks(mesh)
        total=Base.checked_add(total,size(cells,2))
    end
    return 0,0,total
end

function _mesh_element_block(mesh::MixedMesh,msh::Int)
    offset=0
    for (kind,cells,_) in _cache_native_blocks(mesh)
        kind==msh && return offset,cells
        offset=Base.checked_add(offset,size(cells,2))
    end
    return nothing
end

_high_order_overlay(::MixedMesh)=nothing
_append_p2_midnodes!(tags,coords,params,m::GeoModel,mesh::MixedMesh,
                     dim::Int,tag::Int,include_boundary::Bool)=nothing

function _mesh_element_data(mesh::MixedMesh,dimension::Int)
    types=Int32[];tags=Vector{UInt64}[];nodes=Vector{UInt64}[]
    for (msh,dim,offset,cells,_) in _cache_catalog(mesh)
        (dimension<0 || dimension==dim) || continue
        isempty(cells) && continue
        push!(types,msh);push!(tags,_mesh_dense_tags(offset,size(cells,2)))
        push!(nodes,UInt64.(vec(cells)))
    end
    return types,tags,nodes
end

function _mesh_element_types(mesh::MixedMesh,dimension::Int)
    return sort!(unique!([msh for (msh,dim,_,cells,_) in _cache_catalog(mesh)
        if (dimension<0 || dimension==dim) && !isempty(cells)]))
end

function _mixed_classification(cache,entity,entities,node_entities,boundaries,owners)
    return _MeshClassification(cache,entity,entities,node_entities,boundaries,
        get(owners,1,Int32[]),get(owners,2,Int32[]),get(owners,4,Int32[]),owners)
end

function _classify_cached_mesh(m::GeoModel,mesh::MixedMesh,dim::Int,tag::Int,
                               cache::MixedMesh;_extrude_scope=nothing)
    projected=dim==2 ? model_to_mixed(m,mesh,tag;_extrude_scope=_extrude_scope) :
                      model_to_mixed(m,mesh,dim,tag;_extrude_scope=_extrude_scope)
    data=projected.entity_data
    data===nothing && throw(ErrorException(
        "API.mesh.generate: mixed entity projection has no classification"))
    length(data.node_entities)==nnodes(mesh) || throw(ErrorException(
        "API.mesh.generate: mixed classification changed the node set"))
    boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}(
        (dim,tag)=>Int32[abs(b) for b in entity.boundaries]
        for ((dim,tag),entity) in data.entities)
    ownership=Dict{Tuple{Int,Tuple},Int32}()
    for (bi,block) in enumerate(projected.blocks)
        block isa ElementBlock || continue
        for column in axes(block.nodes,2)
            key=(block.msh,Tuple(sort!(collect(@view block.nodes[:,column]))))
            owner=data.block_entities[bi][column]
            previous=get(ownership,key,owner)
            previous==owner || throw(ErrorException(
                "API.mesh.generate: conflicting mixed cell classification"))
            ownership[key]=owner
        end
    end
    owners=Dict{Int,Vector{Int32}}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        values=get!(()->Int32[],owners,msh)
        for column in axes(cells,2)
            key=(msh,Tuple(sort!(collect(@view cells[:,column]))))
            owner=get(ownership,key,Int32(0))
            owner!=0 || throw(ErrorException(
                "API.mesh.generate: mixed cell has no entity owner"))
            push!(values,owner)
        end
    end
    return _mixed_classification(cache,(dim,Int32(tag)),[(dim,Int32(tag))],
                                copy(data.node_entities),boundaries,owners)
end

function _classification_skeleton(mesh::MixedMesh,reversed::Bool)
    reversed || return mesh
    out=_copy_mesh(mesh)
    for block in out.blocks
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        block.msh in (4,5,6,7) || throw(ArgumentError(
            "API: orientation certification is unavailable for MSH $(block.msh)"))
        Model._extrude_make_positive!(out.coords,block.nodes,block.msh)
    end
    return out
end

function _merge_mixed_entity_meshes(parts,caller;node_entities=nothing)
    positions=Dict{Tuple{Int,Int32,NTuple{3,Float64}},Int32}()
    coordinates=NTuple{3,Float64}[]
    remaps=Vector{Vector{Int32}}(undef,length(parts))
    counts=Dict{Int,Int}()
    for (i,(_,mesh)) in enumerate(parts)
        remap=Vector{Int32}(undef,nnodes(mesh))
        for column in axes(mesh.coords,2)
            p=(mesh.coords[1,column],mesh.coords[2,column],mesh.coords[3,column])
            owner=node_entities===nothing || node_entities[i]===nothing ?
                (-1,Int32(0)) : node_entities[i][column]
            key=(owner[1],owner[2],(p[1]+0.0,p[2]+0.0,p[3]+0.0))
            remap[column]=get!(positions,key) do
                length(coordinates)<typemax(Int32) || throw(ArgumentError(
                    "$caller: merged mesh exceeds the Int32 node limit"))
                push!(coordinates,p)
                Int32(length(coordinates))
            end
        end
        remaps[i]=remap
        for (msh,cells,_) in _cache_native_blocks(mesh)
            counts[msh]=Base.checked_add(get(counts,msh,0),size(cells,2))
        end
    end
    nodes=Dict{Int,Matrix{Int32}}(msh=>Matrix{Int32}(undef,Elements.msh_num_nodes(msh),n)
                                for (msh,n) in counts)
    physical=Dict(msh=>Vector{Int32}(undef,n) for (msh,n) in counts)
    cursors=Dict(msh=>0 for msh in keys(counts))
    keeps=Vector{Vector{BitVector}}(undef,length(parts))
    seen_segments=Set{NTuple{2,Int32}}()
    for (i,(_,mesh)) in enumerate(parts)
        blocks=_cache_native_blocks(mesh)
        keeps[i]=[trues(size(cells,2)) for (_,cells,_) in blocks]
        for (bi,(msh,cells,tags)) in enumerate(blocks)
            target=nodes[msh]
            for column in axes(cells,2)
                if msh==1
                    a=remaps[i][cells[1,column]];b=remaps[i][cells[2,column]]
                    key=minmax(a,b)
                    if key in seen_segments
                        keeps[i][bi][column]=false
                        continue
                    end
                    push!(seen_segments,key)
                end
                slot=cursors[msh]+1
                for row in axes(cells,1)
                    target[row,slot]=remaps[i][cells[row,column]]
                end
                physical[msh][slot]=tags[column]
                cursors[msh]=slot
            end
        end
    end
    coords=Matrix{Float64}(undef,3,length(coordinates))
    for (i,p) in enumerate(coordinates),axis in 1:3
        coords[axis,i]=p[axis]
    end
    blocks=ElementBlock[]
    for msh in sort!(collect(keys(nodes)))
        n=cursors[msh]
        n==0 && continue
        push!(blocks,ElementBlock(msh,nodes[msh][:,1:n],physical[msh][1:n]))
    end
    merged=MixedMesh(coords,blocks)
    diagnostic=validate(merged)
    diagnostic.ok || throw(ErrorException(
        "$caller: merged mixed mesh is invalid — "*join(diagnostic.messages,"; ")))
    return merged,remaps,keeps
end

function _merge_classified_parts(parts,remaps,keeps,merged::MixedMesh,caller)
    entities=Tuple{Int,Int32}[]
    node_entities=fill((0,Int32(0)),nnodes(merged))
    boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}()
    owners=Dict{Int,Vector{Int32}}()
    for (i,(tag,part,class)) in enumerate(parts)
        class===nothing && throw(ErrorException(
            "$caller: mixed generation requires entity classification"))
        push!(entities,class.entity)
        for node in eachindex(class.node_entities)
            slot=remaps[i][node]
            node_entities[slot]=_merged_node_owner(node_entities[slot],class.node_entities[node])
        end
        for (key,value) in class.boundaries
            append!(get!(()->Int32[],boundaries,key),value)
            unique!(boundaries[key])
        end
        for (bi,(msh,_,_,cells,cell_owners)) in enumerate(_cache_catalog(part,class))
            append!(get!(()->Int32[],owners,Int(msh)),cell_owners[keeps[i][bi]])
        end
    end
    return _mixed_classification(merged,entities[1],entities,node_entities,boundaries,owners)
end

function _apply_mesh_order(m::GeoModel,cached::MixedMesh,class,caller)
    order=m.meshing.order
    order in (1,2) || throw(ArgumentError(
        "$caller: mesh order $order is not supported (supported: 1, 2)"))
    order==2 && any(block->msh_spec(block.msh).order==1,cached.blocks) &&
        _mixed_require_linear_cad(m,class,caller)
    if all(block->msh_spec(block.msh).family===:pnt ||
                  msh_spec(block.msh).order==order,cached.blocks)
        if order==2
            replacement,new_class=_mixed_prune_unused(cached,class)
            replacement===cached || _replace_mesh_cache_locked!(replacement,new_class)
        end
        return nothing
    end
    class===nothing && throw(ArgumentError(
        "$caller: changing mixed element order requires mesh classification"))
    linear,linear_class=_mixed_linear_cache(cached,class)
    replacement,new_class=order==1 ? (linear,linear_class) :
        _mixed_quadratic_cache(linear,linear_class,caller)
    _replace_mesh_cache_locked!(replacement,new_class)
    return nothing
end

function _mixed_require_linear_cad(m,class,caller)
    class===nothing && return nothing
    entities=Set{Tuple{Int,Int32}}()
    if hasproperty(class,:mesh)
        for (_,dim,_,cells,owners) in _cache_catalog(class.mesh,class)
            dim>0 || continue
            for owner in owners
                push!(entities,(dim,owner))
            end
            for node in cells
                owner=class.node_entities[node]
                owner[1]>0 && push!(entities,owner)
            end
        end
    else
        # The detached assembler supplies only the active native entities
        # requiring new interpolation nodes, with their boundary closure.
        union!(entities,keys(class.boundaries))
        union!(entities,class.node_entities)
    end
    queue=collect(entities)
    position=1
    while position<=length(queue)
        entity=queue[position]
        position+=1
        entity[1]>0 || continue
        for tag in get(class.boundaries,entity,Int32[])
            boundary=(entity[1]-1,abs(tag))
            boundary in entities && continue
            push!(entities,boundary)
            push!(queue,boundary)
        end
    end
    for (dim,entity) in sort!(collect(entities))
        tag=Int(entity)
        if dim==1 && haskey(m.curves,tag)
            kind=get(m.curve_types,tag,:line)
            kind in (:line,:degenerate) && continue
            throw(ArgumentError("$caller: native mixed quadratic CAD placement is not " *
                "implemented for Curve[$tag] ($kind); curved edge nodes and " *
                "face/volume interior placement require geometry stencils"))
        elseif dim==2 && haskey(m.surfaces,tag)
            kind=get(m.surface_types,tag,:plane)
            kind===:plane && continue
            if kind in (:ruled,:tric) && all(loop->
                    all(curve->get(m.curve_types,abs(curve),:line) in
                        (:line,:degenerate),m.loops[loop]),m.surfaces[tag])
                planar=try
                    parentmodule(@__MODULE__).Model._model_surface_plane(
                        m,tag,caller;allow_ruled=true)
                    true
                catch err
                    err isa InterruptException && rethrow()
                    err isa ArgumentError || rethrow()
                    false
                end
                planar && continue
            end
            throw(ArgumentError("$caller: native mixed quadratic CAD placement is not " *
                "implemented for Surface[$tag] ($kind); curved surface and " *
                "volume interior placement require geometry stencils"))
        end
    end
    return nothing
end

function _mixed_linear_cache(mesh,class)
    used=falses(nnodes(mesh))
    referenced=falses(nnodes(mesh))
    blocks=ElementBlock[]
    owners=Dict{Int,Vector{Int32}}()
    for (bi,(msh,_,_,cells,cell_owners)) in enumerate(_cache_catalog(mesh,class))
        spec=msh_spec(msh)
        target=spec.family===:pnt ? 15 : msh_type(spec.family,1)
        n=Elements.msh_num_nodes(target)
        nodes=cells[1:n,:]
        used[nodes].=true
        referenced[cells].=true
        push!(blocks,ElementBlock(target,nodes,mesh.blocks[bi].tags))
        append!(get!(()->Int32[],owners,target),cell_owners)
    end
    # Keep existing orphan primary nodes, while removing unused high-order
    # interpolation nodes when returning from order two to order one.
    all(block->msh_spec(block.msh).family===:pnt ||
               msh_spec(block.msh).order==1,mesh.blocks) && (used.=true)
    used .|= .!referenced
    mapping=zeros(Int32,nnodes(mesh));mapping[used]=Int32.(1:count(used))
    blocks=[ElementBlock(b.msh,mapping[b.nodes],b.tags) for b in blocks]
    result=_mixed_rebuild_metadata(mesh,blocks;node_order=findall(used))
    return result,_mixed_rebind_class(class,result;
        node_entities=class.node_entities[used],owners=owners)
end

function _mixed_quadratic_cache(input_mesh,input_class,caller)
    mesh,class=_mixed_prune_unused(input_mesh,input_class)
    points=NTuple{3,Float64}[(mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])
                             for i in axes(mesh.coords,2)]
    node_entities=copy(class.node_entities)
    keys=Dict{Tuple,Int32}()
    entity_nodes=Tuple{Tuple{Int,Int32},Set{Int}}[]
    for (entity,_) in class.boundaries
        push!(entity_nodes,(entity,Set(Int.(
            _mesh_entity_nodes_with_boundary(class,entity[1],entity[2])))))
    end
    sort!(entity_nodes;by=first)
    blocks=ElementBlock[];owners=Dict{Int,Vector{Int32}}()
    for (bi,(msh,dim,_,cells,cell_owners)) in enumerate(_cache_catalog(mesh,class))
        family=msh_spec(msh).family
        target=family===:pnt ? 15 : msh_type(family,2)
        reference=Elements.lagrange_nodes(target)
        nprimary=size(cells,1)
        nodes=Matrix{Int32}(undef,size(reference,2),size(cells,2))
        nodes[1:nprimary,:]=cells
        for slot in nprimary+1:size(reference,2)
            weights=MeshFunctionSpaces._first_order_values(Val(family),
                reference[1,slot],reference[2,slot],reference[3,slot],caller,slot)
            support=findall(!iszero,collect(weights))
            for column in axes(cells,2)
                terms=sort!([(cells[j,column],Float64(weights[j])) for j in support])
                key=Tuple(terms)
                node=get!(keys,key) do
                    length(points)<typemax(Int32) || throw(ArgumentError(
                        "$caller: elevated mesh exceeds the Int32 node limit"))
                    p=ntuple(axis->sum(weight*mesh.coords[axis,n] for (n,weight) in terms),3)
                    all(isfinite,p) || throw(ArgumentError(
                        "$caller: elevated node coordinates are non-finite"))
                    push!(points,p)
                    owner=(dim,cell_owners[column])
                    for (candidate,allnodes) in entity_nodes
                        candidate[1]>dim && break
                        all(term->Int(term[1]) in allnodes,terms) || continue
                        owner=candidate
                        break
                    end
                    push!(node_entities,owner)
                    Int32(length(points))
                end
                nodes[slot,column]=node
            end
        end
        push!(blocks,ElementBlock(target,nodes,mesh.blocks[bi].tags))
        append!(get!(()->Int32[],owners,target),cell_owners)
    end
    coords=Matrix{Float64}(undef,3,length(points))
    for (node,p) in enumerate(points),axis in 1:3
        coords[axis,node]=p[axis]
    end
    if mesh.entity_data===nothing
        result=MixedMesh(coords,blocks;physical_names=mesh.physical_names)
    else
        data=mesh.entity_data
        first_new=maximum(data.external_node_tags;init=UInt64(0))
        class.public_tags===nothing || (first_new=max(NODE_TAG_MAX[],first_new))
        count=length(points)-nnodes(mesh)
        first_new<=typemax(Int32)-count || throw(ArgumentError(
            "$caller: elevated public node tags exceed Int32"))
        external=vcat(data.external_node_tags,first_new .+ UInt64.(1:count))
        parameters=vcat(data.node_parametric,
            Union{Nothing,Vector{Float64}}[nothing for _ in 1:count])
        result=_mixed_rebuild_metadata(mesh,blocks;coords,node_entities,
            node_parametric=parameters,external_node_tags=external)
    end
    _api_p2_constructed_mesh_certify(result,caller)
    return result,_mixed_rebind_class(class,result;
        node_entities=node_entities,owners=owners)
end

# Gmsh SetOrderN rebuilds mesh-vertex associations from element connectivity.
# Orphan vertices without a point element cease to belong to the mesh.
function _mixed_prune_unused(mesh,class)
    used=falses(nnodes(mesh))
    for (_,cells,_) in _cache_native_blocks(mesh),node in cells
        used[node]=true
    end
    all(used) && return mesh,class
    selections=[collect(axes(cells,2)) for (_,cells,_) in _cache_native_blocks(mesh)]
    return _mixed_select_columns(mesh,class,selections;node_order=findall(used))
end

# Reversal maps every interpolation node through the same reference reflection
# as the primary vertices, preserving the polynomial map for curved elements.
function _mixed_reverse_permutation(msh)
    spec=msh_spec(msh);reference=Elements.lagrange_nodes(msh)
    spec.family===:pnt && return [1]
    primary=Elements.msh_num_nodes(msh_type(spec.family,1))
    corner_order=Int32.(1:primary)
    _record_reverse_connectivity!(Int32(msh_type(spec.family,1)),corner_order)
    result=Vector{Int}(undef,spec.nnodes)
    for node in 1:spec.nnodes
        weights=MeshFunctionSpaces._first_order_values(Val(spec.family),
            reference[1,node],reference[2,node],reference[3,node],"API.mesh.reverse",node)
        reflected=ntuple(axis->sum(weights[i]*reference[axis,corner_order[i]]
                                   for i in 1:primary),3)
        target=findfirst(i->all(axis->isapprox(reference[axis,i],reflected[axis];
                                      atol=64eps(Float64),rtol=64eps(Float64)),1:3),
                         1:spec.nnodes)
        target===nothing && throw(ArgumentError(
            "API.mesh.reverse: unsupported reference reflection for MSH $msh"))
        result[node]=target
    end
    return result
end

function _mixed_reverse_cache(mesh,class;selected=nothing,tags=nothing)
    blocks=ElementBlock[]
    for (bi,(msh,dim,offset,cells,owners)) in enumerate(_cache_catalog(mesh,class))
        nodes=copy(cells)
        permutation=_mixed_reverse_permutation(msh)
        for column in axes(cells,2)
            (selected===nothing || (dim,owners[column]) in selected) || continue
            (tags===nothing || offset+column in tags) || continue
            nodes[:,column]=cells[permutation,column]
        end
        push!(blocks,ElementBlock(msh,nodes,mesh.blocks[bi].tags))
    end
    result=_mixed_with_coordinates(mesh,mesh.coords;blocks=blocks)
    return result,class===nothing ? nothing : _mixed_rebind_class(class,result)
end

function _mixed_with_coordinates(mesh,coords;blocks=mesh.blocks)
    return MixedMesh(coords,blocks;physical_names=mesh.physical_names,
        entity_data=mesh.entity_data,elementary_entities=mesh.elementary_entities,
        periodic_links=mesh.periodic_links,ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
end

function _mixed_rebuild_metadata(mesh,blocks;node_order=nothing,selections=nothing,
        node_entities=nothing,node_parametric=nothing,external_node_tags=nothing,
        external_element_tags=nothing,block_entities=nothing,coords=nothing,
        periodic_links=nothing)
    order=node_order===nothing ? collect(1:nnodes(mesh)) : node_order
    coordinates=coords===nothing ? mesh.coords[:,order] : coords
    if periodic_links===nothing
        if node_order===nothing || isempty(mesh.periodic_links)
            periodic_links=mesh.periodic_links
        else
            remap=zeros(Int32,nnodes(mesh))
            for (new,old) in enumerate(order)
                remap[old]=Int32(new)
            end
            periodic_links=_mixed_remap_periodic_nodes(mesh,remap,
                "API mesh node selection";drop_removed=true)
        end
    end
    data=mesh.entity_data
    kept_blocks=findall(block->!isempty(block.nodes),blocks)
    if data!==nothing
        stored_data=data
        selected=selections===nothing ? [collect(axes(b.nodes,2)) for b in mesh.blocks] : selections
        entities=node_entities===nothing ? data.node_entities[order] : node_entities
        parameters=node_parametric===nothing ? data.node_parametric[order] : node_parametric
        nodes=external_node_tags===nothing ? data.external_node_tags[order] : external_node_tags
        owners=block_entities===nothing ?
            [stored_data.block_entities[i][selected[i]] for i in eachindex(selected)] : block_entities
        cells=external_element_tags===nothing ?
            [stored_data.external_element_tags[i][selected[i]] for i in eachindex(selected)] : external_element_tags
        data=Elements.MixedEntityData(data.entities;node_entities=entities,
            node_parametric=parameters,external_node_tags=nodes,
            block_entities=owners[kept_blocks],external_element_tags=cells[kept_blocks])
    end
    elementary=data===nothing ? nothing : data.block_entities
    return MixedMesh(coordinates,blocks[kept_blocks];physical_names=mesh.physical_names,
        entity_data=data,elementary_entities=elementary,
        periodic_links=periodic_links,ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
end

function _mixed_affine_cache(mesh,class,matrix,translation,mask=nothing;orientation)
    coords=copy(mesh.coords)
    for node in axes(coords,2)
        (mask===nothing || mask[node]) || continue
        p=mesh.coords[:,node]
        for axis in 1:3
            coords[axis,node]=sum(matrix[axis,j]*p[j] for j in 1:3)+translation[axis]
        end
    end
    result=_mixed_with_coordinates(mesh,coords)
    if mask===nothing && orientation<0
        result,_=_mixed_reverse_cache(result,class)
    end
    return result
end

function _mixed_select_columns(mesh,class,selections;node_order=nothing)
    coords=node_order===nothing ? mesh.coords : mesh.coords[:,node_order]
    mapping=node_order===nothing ? nothing : zeros(Int32,nnodes(mesh))
    node_order===nothing || (mapping[node_order]=Int32.(1:length(node_order)))
    blocks=ElementBlock[];owners=Dict{Int,Vector{Int32}}()
    for (bi,(msh,_,_,cells,cell_owners)) in enumerate(_cache_catalog(mesh,class))
        chosen=selections[bi]
        nodes=cells[:,chosen]
        mapping===nothing || (nodes=mapping[nodes])
        push!(blocks,ElementBlock(msh,nodes,mesh.blocks[bi].tags[chosen]))
        class===nothing || append!(get!(()->Int32[],owners,Int(msh)),cell_owners[chosen])
    end
    result=_mixed_rebuild_metadata(mesh,blocks;node_order,selections)
    new_class=class===nothing ? nothing : _mixed_rebind_class(class,result;
        owners=owners,node_entities=node_order===nothing ? class.node_entities :
                                                    class.node_entities[node_order])
    return result,new_class
end

function _mixed_remove_elements(mesh,class,dim,entity,tags,caller)
    catalog=_cache_catalog(mesh,class)
    _,_,count=_mesh_element_offsets(mesh)
    listed=Set{Int}(_cache_element_index(class.public_tags,value,count,caller) for value in tags)
    if !isempty(listed)
        for tag in listed
            _,_,owner_dim,owner=_mixed_element_record(mesh,class,tag)
            (owner_dim==dim && owner==entity) || throw(ArgumentError(
                "$caller: element $tag is not classified on dimension $dim entity $entity"))
        end
    end
    selections=[findall(column->!(bdim==dim && owners[column]==entity &&
                           (isempty(listed) || offset+column in listed)),axes(cells,2))
        for (_,bdim,offset,cells,owners) in catalog]
    result,new_class=_mixed_select_columns(mesh,class,selections)
    _replace_mesh_cache_locked!(result,new_class;
                               preserve_visibility=class.public_tags!==nothing)
    return nothing
end

function _mixed_reorder(mesh,class,msh,entity,ordering,caller)
    catalog=_cache_catalog(mesh,class)
    selections=[collect(axes(cells,2)) for (_,_,_,cells,_) in catalog]
    block=findfirst(entry->entry[1]==msh,catalog)
    positions=block===nothing ? Int[] : findall(==(Int32(entity)),catalog[block][5])
    isempty(positions) && throw(ArgumentError("$caller: no elements of type $msh on entity $entity"))
    (ordering isa AbstractVector || ordering isa Tuple) || throw(ArgumentError(
        "$caller: ordering must be a vector or tuple"))
    permutation=Int[_mesh_query_integer(value,caller,"ordering entry") for value in ordering]
    sort(permutation)==collect(0:length(positions)-1) || throw(ArgumentError(
        "$caller: ordering must be a permutation of 0:$(length(positions)-1)"))
    selections[block][positions]=positions[permutation.+1]
    replacement,new_class=_mixed_select_columns(mesh,class,selections)
    _replace_mesh_cache_locked!(replacement,new_class;
                               preserve_visibility=class.public_tags!==nothing)
    return nothing
end

function _mixed_renumber_elements(mesh,class,mapping,caller)
    catalog=_cache_catalog(mesh,class)
    selections=Vector{Vector{Int}}()
    for (_,_,offset,cells,_) in catalog
        count=size(cells,2);chosen=Vector{Int}(undef,count)
        for column in 1:count
            target=Int(mapping[offset+column])-offset
            1<=target<=count || throw(ArgumentError(
                "$caller: renumbering cannot move elements across element types"))
            chosen[target]=column
        end
        push!(selections,chosen)
    end
    return _mixed_select_columns(mesh,class,selections)
end

function _mixed_element_record(mesh::MixedMesh,class,tag::Int)
    for (msh,dim,offset,cells,owners) in _cache_catalog(mesh,class)
        column=tag-offset
        1<=column<=size(cells,2) || continue
        return msh,UInt64.(cells[:,column]),dim,Int(owners[column])
    end
    throw(ArgumentError("API.mesh.get_element: unknown element $tag"))
end

function _mixed_rebind_class(class,mesh;node_entities=class.node_entities,
                             owners=class.cell_entities)
    result=_mixed_classification(mesh,class.entity,class.entities,node_entities,
                                class.boundaries,owners)
    class.public_tags===nothing && return result
    return _classification_with_public_tags(result,
        _cache_public_tags(mesh;authority=class.public_tags.authority))
end

# The simplex cache keeps its P2 geometry in an overlay. Lower-dimensional
# generation must ingest the actual cells and midnodes rather than its skeleton.
function _dim01_actual_cache(cached,class;physical_names=Dict{Tuple{Int,Int},String}())
    (cached===nothing || cached isa MixedMesh) && return cached,class
    overlay=_high_order_overlay(cached)
    overlay===nothing && return cached,class
    blocks=ElementBlock[];owners=Dict{Int,Vector{Int32}}()
    for (msh,dim,_,cells,classified) in _cache_catalog(cached,class)
        isempty(cells) && continue
        if msh==_p2_skeleton_etype(overlay)
            msh=_p2_etype(overlay);cells=_p2_cells(overlay)
        end
        physical=dim==1 ? cached.seg_tag : dim==2 ? cached.tri_tag : cached.tet_tag
        push!(blocks,ElementBlock(msh,cells,physical))
        owners[Int(msh)]=classified
    end
    actual=MixedMesh(overlay.coords,blocks;physical_names)
    class===nothing && return actual,nothing
    actual_class=_mixed_classification(actual,class.entity,class.entities,
        vcat(class.node_entities,LAST_MESH_HIGH_ORDER_MIDS[]),class.boundaries,owners)
    return actual,actual_class
end

function _clear_classified_mesh(mesh::MixedMesh,class::_MeshClassification,
                                cleared::Set{Tuple{Int,Int32}})
    catalog=_cache_catalog(mesh,class)
    kept=[findall(owner->!((dim,owner) in cleared),owners)
          for (_,dim,_,_,owners) in catalog]
    removed=Set{Tuple{Int,Int32}}()
    union!(removed,intersect(cleared,Set(class.node_entities)))
    survivors=fill((4,Int32(0)),nnodes(mesh))
    referenced=falses(nnodes(mesh))
    for (bi,(_,dim,_,cells,owners)) in enumerate(catalog)
        for owner in owners
            (dim,owner) in cleared && push!(removed,(dim,owner))
        end
        for column in kept[bi],row in axes(cells,1)
            node=cells[row,column]
            referenced[node]=true
            survivors[node]=min(survivors[node],(dim,owners[column]))
        end
    end
    isempty(removed) && return mesh,class
    keep=BitVector([!(owner in removed) || referenced[node]
                    for (node,owner) in enumerate(class.node_entities)])
    remap=zeros(Int32,nnodes(mesh))
    next=Int32(0)
    for node in eachindex(keep)
        keep[node] && (next+=1;remap[node]=next)
    end
    node_entities=Tuple{Int,Int32}[
        owner in removed ? survivors[node] : owner
        for (node,owner) in enumerate(class.node_entities) if keep[node]]
    blocks=ElementBlock[]
    owners=Dict{Int,Vector{Int32}}()
    for (bi,(msh,dim,_,cells,cell_owners)) in enumerate(catalog)
        block=mesh.blocks[bi]
        push!(blocks,ElementBlock(msh,remap[cells[:,kept[bi]]],block.tags[kept[bi]]))
        append!(get!(()->Int32[],owners,Int(msh)),cell_owners[kept[bi]])
    end
    replacement=_mixed_rebuild_metadata(mesh,blocks;node_order=findall(keep),
        selections=kept,node_entities=node_entities)
    return replacement,_mixed_rebind_class(class,replacement;
                        node_entities=node_entities,owners=owners)
end

function _mixed_refresh_periodic_links(model,mesh,class,caller)
    isempty(mesh.periodic_links) && return mesh,class
    links=Elements.MixedPeriodicLink[]
    for link in mesh.periodic_links
        if link.dim==1 && haskey(model.periodic,(1,Int(link.slave_entity)))
            mapping=_dim01_periodic_nodes(model,mesh,class,Int(link.slave_entity),caller;
                                          include_high_order=true)
            push!(links,Elements.MixedPeriodicLink(1,link.slave_entity,
                mapping.master_entity,mapping.slave_nodes,mapping.master_nodes;
                affine=mapping.affine))
        else
            push!(links,link)
        end
    end
    replacement=_mixed_rebuild_metadata(mesh,mesh.blocks;periodic_links=links)
    return replacement,_mixed_rebind_class(class,replacement)
end
