# Operation-owned topology only. A `nothing` result leaves a legacy imported
# state to its existing fallback; no Cartesian value or tolerance is inspected.
struct _Dim01PeriodicLineCell
    msh::Int
    native_nodes::Vector{Int32}
    directed_nodes::Vector{Int32}
end

function _dim01_periodic_public_tags(mesh,class)
    if mesh isa MixedMesh && mesh.entity_data!==nothing
        return mesh.entity_data.external_node_tags
    end
    class!==nothing && class.public_tags!==nothing || return nothing
    return class.public_tags.node_tags
end

function _dim01_periodic_line_cells(m,mesh,class,curve,caller)
    result=_Dim01PeriodicLineCell[]
    for (msh,dim,_,cells,owners) in _cache_catalog(mesh,class)
        dim==1 || continue
        msh_family(msh)===:lin || return nothing
        for column in axes(cells,2)
            owners[column]==curve || continue
            native=Int32.(cells[:,column])
            length(native)==msh_spec(msh).nnodes &&
                all(row->1<=row<=nnodes(mesh),native) || throw(ArgumentError(
                "$caller: periodic Curve[$curve] has invalid Line connectivity"))
            native[1]!=native[2] || return nothing
            directed=get(m.meshing.reverse,(1,curve),false) ?
                native[_mixed_reverse_permutation(msh)] : copy(native)
            push!(result,_Dim01PeriodicLineCell(Int(msh),native,directed))
        end
    end
    return result
end

function _dim01_periodic_line_rows(cells,include_high_order)
    result=Set{Int32}()
    for cell in cells
        width=include_high_order ? length(cell.native_nodes) : 2
        union!(result,@view cell.native_nodes[1:width])
    end
    return result
end

function _dim01_periodic_existing_uid_pairs(m,mesh,class,slave,master,
        slave_set,master_set,caller;source_mesh=nothing,node_tag_map=nothing)
    tags=_dim01_periodic_public_tags(mesh,class)
    tags===nothing && return nothing
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(tags))
    sources=Union{Mesh,MixedMesh}[mesh]
    source_mesh===nothing || push!(sources,source_mesh)
    authority=m.meshing.mesh_authority
    authority===nothing || authority.mesh===nothing || push!(sources,authority.mesh)
    for source in sources
        source isa MixedMesh && source.entity_data!==nothing || continue
        source_tags=source.entity_data.external_node_tags
        for link in source.periodic_links
            link.dim==1 && link.slave_entity==slave && link.master_entity==master || continue
            pairs=Dict{Int32,Int32}()
            for (srow,mrow) in zip(link.slave_nodes,link.master_nodes)
                1<=srow<=length(source_tags) && 1<=mrow<=length(source_tags) ||
                    throw(ArgumentError("$caller: stored periodic link has invalid node rows"))
                stag=source_tags[srow];mtag=source_tags[mrow]
                if source!==mesh && node_tag_map!==nothing
                    stag=get(node_tag_map,stag,UInt64(0))
                    mtag=get(node_tag_map,mtag,UInt64(0))
                end
                current_slave=get(positions,stag,Int32(0))
                current_master=get(positions,mtag,Int32(0))
                current_slave in slave_set && current_master in master_set || continue
                old=get(pairs,current_slave,Int32(0))
                old==0 || old==current_master || throw(ArgumentError(
                    "$caller: stored periodic UID correspondence is ambiguous"))
                pairs[current_slave]=current_master
            end
            length(pairs)==length(slave_set) || continue
            Set(values(pairs))==master_set || throw(ArgumentError(
                "$caller: stored periodic UID correspondence is not a bijection"))
            return pairs
        end
    end
    return nothing
end

function _dim01_periodic_endpoint_row(m,mesh,class,curve,point,last,
        primary_set,caller;curve_births=nothing,node_tag_map=nothing)
    tags=_dim01_periodic_public_tags(mesh,class)
    tags===nothing && return nothing
    positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(tags))
    owner=(1,Int32(curve))
    if curve_births!==nothing
        births=get(curve_births,owner,nothing)
        if births!==nothing && !isempty(births)
            tag=last ? births[end] : births[begin]
            node_tag_map===nothing || (tag=get(node_tag_map,tag,UInt64(0)))
            row=get(positions,tag,Int32(0))
            row in primary_set && return row
        end
    end
    receipt=get(m.meshing.native_curve_endpoints,point,nothing)
    if receipt!==nothing
        tag=receipt.node
        node_tag_map===nothing || (tag=get(node_tag_map,tag,UInt64(0)))
        row=get(positions,tag,Int32(0))
        row in primary_set && return row
    end
    authority=m.meshing.mesh_authority
    if authority!==nothing && authority.mesh!==nothing
        births=get(authority.producer_rows,owner,nothing)
        if births!==nothing && !isempty(births)
            index=last ? births[end] : births[begin]
            tag=authority.mesh.entity_data.external_node_tags[index]
            node_tag_map===nothing || (tag=get(node_tag_map,tag,UInt64(0)))
            row=get(positions,tag,Int32(0))
            row in primary_set && return row
        end
    end
    candidates=Int32[row for row in primary_set
        if class.node_entities[row]==(0,Int32(point))]
    length(candidates)==1 && return only(candidates)
    empty!(candidates)
    for (msh,dim,_,cells,owners) in _cache_catalog(mesh,class)
        dim==0 && msh==15 || continue
        for column in axes(cells,2)
            owners[column]==point || continue
            row=cells[1,column]
            row in primary_set && !(row in candidates) && push!(candidates,row)
        end
    end
    return length(candidates)==1 ? only(candidates) : nothing
end

function _dim01_periodic_directed_chain(m,mesh,class,curve,cells,caller;
        curve_births=nothing,node_tag_map=nothing)
    endpoints=get(m.curves,curve,nothing)
    endpoints===nothing && return nothing
    primary_set=_dim01_periodic_line_rows(cells,false)
    start=_dim01_periodic_endpoint_row(m,mesh,class,curve,endpoints[1],false,
        primary_set,caller;curve_births,node_tag_map)
    finish=endpoints[1]==endpoints[2] ? start :
        _dim01_periodic_endpoint_row(m,mesh,class,curve,endpoints[2],true,
            primary_set,caller;curve_births,node_tag_map)
    start===nothing || finish===nothing ? (return nothing) : nothing
    outgoing=Dict{Int32,Int}();incoming=Dict{Int32,Int}()
    for (position,cell) in enumerate(cells)
        first,last=cell.directed_nodes[1:2]
        haskey(outgoing,first) && return nothing
        haskey(incoming,last) && return nothing
        outgoing[first]=position;incoming[last]=position
    end
    closed=endpoints[1]==endpoints[2]
    if closed
        Set(keys(outgoing))==primary_set && Set(keys(incoming))==primary_set || return nothing
    else
        haskey(incoming,start) && return nothing
        haskey(outgoing,finish) && return nothing
        union(Set(keys(outgoing)),Set([finish]))==primary_set || return nothing
        union(Set(keys(incoming)),Set([start]))==primary_set || return nothing
    end
    ordered=_Dim01PeriodicLineCell[];visited=Set{Int}();current=start
    while length(ordered)<length(cells)
        position=get(outgoing,current,0)
        position>0 && !(position in visited) || return nothing
        cell=cells[position]
        push!(ordered,cell);push!(visited,position)
        current=cell.directed_nodes[2]
    end
    current==finish || return nothing
    return ordered
end

function _dim01_periodic_add_uid_pair!(pairs,inverse,slave,master,caller)
    old=get(pairs,slave,Int32(0));other=get(inverse,master,Int32(0))
    (old==0 || old==master) && (other==0 || other==slave) || throw(ArgumentError(
        "$caller: periodic Line topology does not establish a node bijection"))
    pairs[slave]=master;inverse[master]=slave
    return nothing
end

function _dim01_periodic_uid_nodes(m,mesh,class,slave::Int,caller;
        include_high_order=false,source_mesh=nothing,node_tag_map=nothing,
        curve_births=nothing)
    constraint=get(m.periodic,(1,slave),nothing)
    constraint===nothing && return nothing
    master=abs(Int(constraint.master_entity))
    class===nothing && return nothing
    _dim01_periodic_public_tags(mesh,class)===nothing && return nothing
    master_cells=_dim01_periodic_line_cells(m,mesh,class,master,caller)
    slave_cells=_dim01_periodic_line_cells(m,mesh,class,slave,caller)
    master_cells===nothing || slave_cells===nothing ? (return nothing) : nothing
    isempty(master_cells) || isempty(slave_cells) ? (return nothing) : nothing
    master_set=_dim01_periodic_line_rows(master_cells,include_high_order)
    slave_set=_dim01_periodic_line_rows(slave_cells,include_high_order)
    pairs=_dim01_periodic_existing_uid_pairs(m,mesh,class,slave,master,
        slave_set,master_set,caller;source_mesh,node_tag_map)
    if pairs===nothing
        master_chain=_dim01_periodic_directed_chain(m,mesh,class,master,master_cells,caller;
            curve_births,node_tag_map)
        slave_chain=_dim01_periodic_directed_chain(m,mesh,class,slave,slave_cells,caller;
            curve_births,node_tag_map)
        master_chain===nothing || slave_chain===nothing ? (return nothing) : nothing
        length(master_chain)==length(slave_chain) || throw(ArgumentError(
            "$caller: periodic Curve[$slave]/Curve[$master] cell counts differ"))
        constraint.reversed && reverse!(slave_chain)
        pairs=Dict{Int32,Int32}();inverse=Dict{Int32,Int32}()
        for (master_cell,slave_cell) in zip(master_chain,slave_chain)
            master_cell.msh==slave_cell.msh || throw(ArgumentError(
                "$caller: corresponding periodic Lines have different native factories"))
            ma,mb=master_cell.directed_nodes[1:2]
            sa,sb=slave_cell.directed_nodes[1:2]
            constraint.reversed && ((sa,sb)=(sb,sa))
            _dim01_periodic_add_uid_pair!(pairs,inverse,sa,ma,caller)
            _dim01_periodic_add_uid_pair!(pairs,inverse,sb,mb,caller)
            if include_high_order
                # Generator.cpp FixPeriodicMesh matches undirected endpoint
                # edges, then copies local interpolation slots verbatim.
                # Reversing these slots would contradict the captured native
                # reversed/order5 correspondence, even when XYZ suggests it.
                for slot in 3:length(master_cell.native_nodes)
                    _dim01_periodic_add_uid_pair!(pairs,inverse,
                        slave_cell.native_nodes[slot],master_cell.native_nodes[slot],caller)
                end
            end
        end
        Set(keys(pairs))==slave_set && Set(values(pairs))==master_set || throw(ArgumentError(
            "$caller: periodic Line topology leaves unmatched node identities"))
    end
    master_nodes=sort!(collect(master_set))
    inverse=Dict(master=>slave for (slave,master) in pairs)
    slave_nodes=Int32[inverse[master] for master in master_nodes]
    return (master_entity=master,slave_nodes=slave_nodes,master_nodes=master_nodes,
            affine=constraint.affine)
end
