# Native mixed cache operations. Connectivity keeps its actual MSH family.
function _mixed_remove_duplicate_nodes!(model,cached,pairs,records,caller)
    class=_cached_classification_locked(cached)
    (isempty(pairs) || class!==nothing ||
     all(pair->haskey(model.discrete,pair) || haskey(model.meshing.attached,pair),pairs)) ||
        throw(ArgumentError("$caller: selective duplicate removal requires classification metadata"))
    selected=(isempty(pairs) || class===nothing) ? nothing :
             _mesh_selected_entities(model,class,pairs,caller)
    scan_cache=class!==nothing || isempty(pairs)
    record_pairs=isempty(pairs) ? nothing : collect(Tuple{Int,Int},pairs)
    count=nnodes(cached)
    replacement=collect(Int32,1:count)
    groups=Dict{NTuple{3,Float64},Int32}()
    for node in 1:count
        (scan_cache && (selected===nothing || class.node_entities[node] in selected)) || continue
        point=(cached.coords[1,node],cached.coords[2,node],cached.coords[3,node])
        replacement[node]=get!(groups,point,Int32(node))
    end
    keep=BitVector([replacement[node]==node for node in 1:count])
    remap=zeros(Int32,count)
    next=Int32(0)
    for node in 1:count
        keep[node] && (next+=1;remap[node]=next)
    end
    for node in 1:count
        remap[node]=remap[replacement[node]]
    end
    if !all(keep)
        blocks=ElementBlock[ElementBlock(block.msh,remap[block.nodes],block.tags)
                            for block in cached.blocks]
        mesh=MixedMesh(cached.coords[:,keep],blocks;physical_names=cached.physical_names)
        new_class=class===nothing ? nothing : _mixed_rebind_class(class,mesh;
            node_entities=class.node_entities[keep])
        for (_,_,record) in records,connectivity in record.element_nodes,i in eachindex(connectivity)
            tag=connectivity[i]
            1<=tag<=count && (connectivity[i]=remap[tag])
        end
        _replace_mesh_cache_locked!(mesh,new_class)
        cached=mesh
    end
    isempty(records) && return nothing
    keepers=Dict{NTuple{3,Float64},Tuple{Int32,Bool}}()
    for node in 1:count
        keep[node] || continue
        (scan_cache && (selected===nothing || class.node_entities[node] in selected)) || continue
        point=(cached.coords[1,remap[node]],cached.coords[2,remap[node]],cached.coords[3,remap[node]])
        keepers[point]=(remap[node],true)
    end
    _record_dedup_nodes!(records,record_pairs,keepers)
    return nothing
end

function _mixed_remove_duplicate_elements!(model,cached,pairs,records,caller)
    class=_cached_classification_locked(cached)
    class===nothing && throw(ArgumentError("$caller: duplicate-element removal requires classification metadata"))
    selected=isempty(pairs) ? nothing : _mesh_selected_entities(model,class,pairs,caller)
    record_selected=isempty(pairs) ? nothing : Set{Tuple{Int,Int}}(pairs)
    seen=Set{Tuple{Int32,Int32,Tuple}}()
    blocks=ElementBlock[];owners=Dict{Int,Vector{Int32}}()
    changed=false
    for (bi,(msh,dim,_,cells,cell_owners)) in enumerate(_cache_catalog(cached,class))
        keep=trues(size(cells,2))
        for column in axes(cells,2)
            (selected===nothing || (dim,cell_owners[column]) in selected) || continue
            key=(msh,cell_owners[column],Tuple(sort!(collect(@view cells[:,column]))))
            if key in seen
                keep[column]=false;changed=true
            else
                push!(seen,key)
            end
        end
        push!(blocks,ElementBlock(msh,cells[:,keep],cached.blocks[bi].tags[keep]))
        append!(get!(()->Int32[],owners,Int(msh)),cell_owners[keep])
    end
    if changed
        mesh=MixedMesh(cached.coords,blocks;physical_names=cached.physical_names)
        class=_mixed_rebind_class(class,mesh;owners)
        _replace_mesh_cache_locked!(mesh,class)
        cached=mesh
    end
    for (dim,tag,record) in records
        (record_selected===nothing || (dim,tag) in record_selected) || continue
        seed=nothing
        if !haskey(model.discrete,(dim,tag))
            seed=Vector{Int32}[]
            for (_,dimension,_,cells,cell_owners) in _cache_catalog(cached,class)
                dimension==dim || continue
                for column in axes(cells,2)
                    cell_owners[column]==tag && push!(seed,sort!(collect(@view cells[:,column])))
                end
            end
        end
        _record_dedup_elements!(record,seed)
    end
    return nothing
end

function _partition_elements(mesh::MixedMesh,num_partitions::Int)
    blocks=_cache_native_blocks(mesh)
    total=last(_mesh_element_offsets(mesh))
    incidence=[Int[] for _ in 1:nnodes(mesh)]
    block_ids=Vector{Int}(undef,total);columns=Vector{Int}(undef,total)
    tag=0
    for (bi,(_,cells,_)) in enumerate(blocks),column in axes(cells,2)
        tag+=1;block_ids[tag]=bi;columns[tag]=column
        for row in axes(cells,1)
            push!(incidence[cells[row,column]],tag)
        end
    end
    assignment=zeros(Int32,total)
    quota=cld(total,num_partitions)
    assigned=0;seed=1;part=0
    queue=Int[];neighbors=Int[]
    while assigned<total
        part=min(part+1,num_partitions)
        target=part==num_partitions ? total : min(assigned+quota,total)
        while assigned<target
            empty!(queue)
            while seed<=total && assignment[seed]!=0;seed+=1;end
            seed>total && break
            push!(queue,seed);assignment[seed]=Int32(part);assigned+=1
            head=1
            while head<=length(queue) && assigned<target
                element=queue[head];head+=1
                cells=blocks[block_ids[element]][2];column=columns[element]
                empty!(neighbors)
                for row in axes(cells,1)
                    append!(neighbors,incidence[cells[row,column]])
                end
                sort!(unique!(neighbors))
                for neighbor in neighbors
                    assignment[neighbor]==0 || continue
                    assignment[neighbor]=Int32(part);assigned+=1
                    push!(queue,neighbor)
                    assigned<target || break
                end
            end
        end
    end
    return assignment
end

function _mesh_entity_element_tags(class,model,dimension::Int,entity::Int,cached::MixedMesh)
    result=Int[]
    for (_,dim,offset,_,owners) in _cache_catalog(cached,class)
        dim==dimension || continue
        for (column,owner) in enumerate(owners)
            owner==entity && push!(result,offset+column)
        end
    end
    return result
end
