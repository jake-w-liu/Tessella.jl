# Lift a recorded Surface vertex bijection through the actual cell stencils.
# All keys contain current primary rows and integer interpolation weights;
# coordinates never establish a vertex or a periodic correspondence.
function _periodic_surface_stencil(corners,coordinate,order,family)
    u,v=coordinate
    weights=family===:tri ? (order-u-v,u,v) :
        ((order-u)*(order-v),u*(order-v),u*v,(order-u)*v)
    divisor=foldl(gcd,weights;init=0)
    divisor>0 || throw(ArgumentError("periodic Surface stencil has no primary support"))
    return Tuple(sort!([(corners[i],weights[i]÷divisor)
        for i in eachindex(corners) if weights[i]!=0];by=first))
end

function _periodic_surface_catalog(mesh,class,surface,caller)
    stencils=Dict{Int32,Tuple}()
    rows_by_stencil=Dict{Tuple,Int32}()
    primary=Set{Int32}()
    topology=Dict{Tuple,Int}()
    native_cells=NamedTuple[]
    for (msh,dim,_,cells,owners) in _cache_catalog(mesh,class)
        dim==2 || continue
        spec=msh_spec(msh)
        any(==(surface),owners) || continue
        spec.family in (:tri,:qua) || throw(ArgumentError(
            "$caller: periodic Surface[$surface] has unsupported cell topology"))
        coordinates=getfield(Elements,Symbol("_",spec.family,"_monomials"))(
            spec.order,spec.serendipity)
        length(coordinates)==size(cells,1) || throw(ArgumentError(
            "$caller: periodic Surface[$surface] has invalid interpolation slots"))
        nprimary=Elements._msh_primary_nodes(spec.family)
        for column in axes(cells,2)
            owners[column]==surface || continue
            corners=Int32[cells[i,column] for i in 1:nprimary]
            length(unique(corners))==nprimary || throw(ArgumentError(
                "$caller: periodic Surface[$surface] has repeated primary identities"))
            union!(primary,corners)
            push!(native_cells,(msh=Int(msh),family=spec.family,order=spec.order,
                serendipity=spec.serendipity,corners=corners,nodes=Int32[cells[:,column]...]))
            edges=Tuple(sort!([minmax(corners[a+1],corners[b+1])
                for (a,b) in _native_order_edges(spec.family)]))
            cellkey=(spec.family,edges)
            topology[cellkey]=get(topology,cellkey,0)+1
            for (slot,coordinate) in enumerate(coordinates)
                row=cells[slot,column]
                1<=row<=nnodes(mesh) || throw(ArgumentError(
                    "$caller: periodic Surface[$surface] has an invalid node row"))
                stencil=_periodic_surface_stencil(corners,coordinate,spec.order,spec.family)
                previous=get(stencils,row,nothing)
                previous===nothing || previous==stencil || throw(ArgumentError(
                    "$caller: periodic Surface[$surface] has conflicting support identities"))
                other=get(rows_by_stencil,stencil,row)
                other==row || throw(ArgumentError(
                    "$caller: periodic Surface[$surface] has duplicate support identities"))
                stencils[row]=stencil;rows_by_stencil[stencil]=row
            end
        end
    end
    return (;stencils,rows_by_stencil,primary,topology,cells=native_cells)
end

function _periodic_surface_lift(mesh,class,link,seeds,caller;
        parent_receipt=nothing,source=nothing,node_tag_map=nothing,element_tag_map=nothing)
    slave=_periodic_surface_catalog(mesh,class,link.slave_entity,caller)
    master=_periodic_surface_catalog(mesh,class,link.master_entity,caller)
    if parent_receipt!==nothing
        parent_receipt isa _SurfacePeriodicRefineReceipt || throw(ArgumentError(
            "$caller: periodic refinement requires an actual parent receipt"))
        _surface_periodic_refine_validate(parent_receipt,source,mesh,class,
            node_tag_map,element_tag_map,caller)
        pairs=Dict(seeds)
        Set(keys(pairs))==Set(keys(slave.stencils)) &&
            Set(values(pairs))==Set(keys(master.stencils)) &&
            length(unique(values(pairs)))==length(pairs) || throw(ArgumentError(
                "$caller: completed periodic parent loses an actual child vertex"))
        return seeds
    end
    isempty(slave.primary) && isempty(master.primary) && return seeds
    roots=Dict{Int32,Int32}()
    for (srow,mrow) in seeds
        srow in slave.primary && mrow in master.primary || continue
        previous=get(roots,srow,mrow)
        previous==mrow || throw(ArgumentError("$caller: periodic primary map is ambiguous"))
        roots[srow]=mrow
    end
    Set(keys(roots))==slave.primary && Set(values(roots))==master.primary &&
        length(unique(values(roots)))==length(roots) || throw(ArgumentError(
            "$caller: periodic Surface link lacks a complete primary UID bijection"))
    mapped_topology=Dict{Tuple,Int}()
    for ((family,edges),count) in slave.topology
        mapped=Tuple(sort!([minmax(roots[a],roots[b]) for (a,b) in edges]))
        key=(family,mapped)
        mapped_topology[key]=get(mapped_topology,key,0)+count
    end
    mapped_topology==master.topology || throw(ArgumentError(
        "$caller: periodic Surface cell incidence differs under its primary UID map"))
    pairs=Dict{Int32,Int32}()
    master_cells=Dict{Tuple,Any}()
    for cell in master.cells
        key=(cell.family,Tuple(sort(cell.corners)))
        haskey(master_cells,key) && throw(ArgumentError(
            "$caller: periodic master has ambiguous cell identities"))
        master_cells[key]=cell
    end
    for cell in slave.cells
        mapped=Int32[roots[row] for row in cell.corners]
        counterpart=get(master_cells,(cell.family,Tuple(sort(mapped))),nothing)
        counterpart!==nothing && counterpart.order==cell.order &&
            counterpart.serendipity==cell.serendipity || throw(ArgumentError(
                "$caller: periodic Surface interpolation topology differs"))
        positions=Int[only(findall(==(row),counterpart.corners)) for row in mapped]
        n=length(positions)
        forward=all(i->positions[i]==mod1(positions[1]+i-1,n),1:n)
        reflected=all(i->positions[i]==mod1(positions[1]-i+1,n),1:n)
        forward || reflected || throw(ArgumentError(
            "$caller: periodic cell primary map has no cyclic orientation"))
        coordinates=getfield(Elements,Symbol("_",cell.family,"_monomials"))(
            cell.order,cell.serendipity)
        # The native periodic HighOrder map keeps its independent primary
        # pairs. New interpolation slots use local ordinal order, with the
        # reference u/v reflection when the matched cell is reversed.
        # Cyclic starting-corner shifts do not rotate this separate map.
        for slot in eachindex(cell.nodes)
            srow=cell.nodes[slot]
            mrow=if slot<=n
                roots[srow]
            else
                coordinate=coordinates[slot]
                target=forward ? slot : only(findall(==((coordinate[2],coordinate[1])),coordinates))
                counterpart.nodes[target]
            end
            previous=get(pairs,srow,mrow)
            previous==mrow || throw(ArgumentError(
                "$caller: periodic shared interpolation slot has conflicting native maps"))
            pairs[srow]=mrow
        end
    end
    Set(values(pairs))==Set(keys(master.stencils)) &&
        length(unique(values(pairs)))==length(pairs) || throw(ArgumentError(
            "$caller: periodic Surface interpolation map is not a bijection"))
    ordered=Pair{Int32,Int32}[]
    seen=Set{Int32}()
    for (srow,mrow) in seeds
        if haskey(pairs,srow)
            pairs[srow]==mrow || throw(ArgumentError(
                "$caller: recorded periodic support disagrees with its cell stencil"))
        end
        srow in seen && continue
        push!(ordered,srow=>mrow);push!(seen,srow)
    end
    for srow in sort!(collect(keys(pairs)))
        srow in seen && continue
        push!(ordered,srow=>pairs[srow]);push!(seen,srow)
    end
    return ordered
end

function _periodic_surface_order_links(source,mesh,class,caller;
        node_tag_map=nothing,element_tag_map=nothing,parent_receipt=nothing)
    source isa MixedMesh || return mesh.periodic_links
    source_tags=_dim01_periodic_public_tags(source,nothing)
    if source_tags===nothing && parent_receipt isa _SurfacePeriodicRefineReceipt
        source===parent_receipt.source || throw(ArgumentError("$caller: periodic parent identity changed"))
        source_tags=parent_receipt.source_class.public_tags.node_tags
    end
    target_tags=_dim01_periodic_public_tags(mesh,class)
    source_tags===nothing && return mesh.periodic_links
    target_tags===nothing && throw(ArgumentError("$caller: periodic Surface UID table is missing"))
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(target_tags))
    links=Elements.MixedPeriodicLink[link for link in mesh.periodic_links if link.dim!=2]
    for link in source.periodic_links
        link.dim==2 || continue
        isempty(link.slave_nodes) && isempty(link.master_nodes) && continue
        seeds=Pair{Int32,Int32}[]
        for (srow,mrow) in zip(link.slave_nodes,link.master_nodes)
            1<=srow<=length(source_tags) && 1<=mrow<=length(source_tags) ||
                throw(ArgumentError("$caller: recorded periodic link has invalid rows"))
            suid=source_tags[srow];muid=source_tags[mrow]
            if node_tag_map!==nothing
                haskey(node_tag_map,suid) && haskey(node_tag_map,muid) || continue
                suid=node_tag_map[suid];muid=node_tag_map[muid]
            end
            s=get(positions,suid,Int32(0));m=get(positions,muid,Int32(0))
            s>0 && m>0 || continue
            push!(seeds,s=>m)
        end
        pairs=_periodic_surface_lift(mesh,class,link,seeds,caller;parent_receipt,source,node_tag_map,element_tag_map)
        push!(links,Elements.MixedPeriodicLink(2,link.slave_entity,link.master_entity,
            Int32[first(pair) for pair in pairs],Int32[last(pair) for pair in pairs];affine=link.affine))
    end
    return links
end
