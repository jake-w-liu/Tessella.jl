# Actual lower Tri3 BRep transport. Interface seams are derived from literal
# node UID incidence, with one authoritative original cell per UID.
const _PARTITION_TET_FACETS=((1,3,2),(1,2,4),(2,3,4),(3,1,4))
const _PARTITION_TRI_EDGES=((1,2),(2,3),(3,1))
_partition_support(uids)=Tuple(sort!(collect(uids)))
_partition_cycle_equal(a,b)=any(a==ntuple(i->b[mod1(i+shift,3)],3) for shift in 0:2)
_partition_cycle_reverse(a,b)=_partition_cycle_equal(a,(b[1],b[3],b[2]))

function _partition_brep_surface_borders(triangles,caller)
    incidences=Dict{Int32,Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}}()
    for tri in triangles,(a,b) in _PARTITION_TRI_EDGES
        edges=get!(()->Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}(),incidences,tri.parent)
        edge=(tri.uids[a],tri.uids[b]);push!(get!(()->Tuple{UInt64,UInt64}[],edges,minmax(edge...)),edge)
    end
    borders=Dict{Int32,Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}}()
    for (parent,edges) in incidences
        border=Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}()
        for (key,hits) in edges
            length(hits)<=2 || throw(ArgumentError("$caller: original Surface has a nonmanifold edge"))
            if length(hits)==2
                hits[1]==reverse(hits[2]) || throw(ArgumentError("$caller: original Surface edge winding disagrees"))
            else
                border[key]=only(hits)
            end
        end
        borders[parent]=border
    end
    return borders
end

function _partition_brep_snapshot(model,cached,class,assignment,caller)
    class.public_tags!==nothing || throw(ArgumentError("$caller: lower BRep needs actual UID authority"))
    _high_order_overlay(cached)===nothing || throw(ArgumentError("$caller: lower BRep partitioning requires first-order cells"))
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    all(row->row[1] in (1,2,4,15),catalog) && any(row->row[1]==4,catalog) || throw(ArgumentError(
        "$caller: lower BRep partitioning supports Point15/Line2/Tri3/Tet4 only"))
    cached isa MixedMesh && (cached.partition_data!==nothing || !isempty(cached.periodic_links)) &&
        throw(ArgumentError("$caller: existing partition or periodic topology is unsupported"))
    isempty(model.periodic) && isempty(model.meshing.periodic_uid_history) ||
        throw(ArgumentError("$caller: periodic lower BRep partitioning is unsupported"))
    all(value->value===nothing || isempty(value),class.node_parametric) ||
        throw(ArgumentError("$caller: parametric lower BRep partitioning is unsupported"))
    public=class.public_tags
    length(public.element_tags)<=typemax(Int32) || throw(ArgumentError("$caller: lower BRep cell capacity exceeded"))
    for (_,record) in Model._discrete_mesh_records_model(model)
        all(uid->haskey(public.node_indices,uid),record.node_tags) &&
        all(uid->haskey(public.element_indices,uid),record.element_tags) || throw(ArgumentError(
            "$caller: lower BRep partitioning cannot omit caller-owned records outside the cache"))
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    _dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    length(assembly.nodes)==nnodes(cached) && length(assembly.cells)==length(public.element_tags) ||
        throw(ArgumentError("$caller: cached and record BRep rows do not agree"))
    entities,_=_dim01_entities(model,model_entities(model),assembly,eachindex(assembly.nodes))
    blocks=ElementBlock[];owners=Vector{Int32}[];uids=Vector{UInt64}[]
    triangles=NamedTuple[];tets=NamedTuple[];lines=NamedTuple[];points=NamedTuple[]
    for (msh,dim,offset,cells,classified) in catalog
        physical=Int32[_dim01_physical(model,(dim,Int(owner))) for owner in classified]
        push!(blocks,ElementBlock(msh,cells,physical));push!(owners,copy(classified))
        push!(uids,copy(public.element_tags[offset+1:offset+size(cells,2)]))
        for column in axes(cells,2)
            rows=Tuple(cells[:,column]);labels=Tuple(public.node_tags[collect(rows)]);index=offset+column
            record=(;rows,uids=labels,index,parent=classified[column],part=assignment[index])
            push!(msh==2 ? triangles : msh==4 ? tets : msh==1 ? lines : points,record)
        end
    end
    parents=sort!(unique(tet.parent for tet in tets));surfaces=sort!(unique(tri.parent for tri in triangles))
    curve_parents=sort!(unique(line.parent for line in lines))
    point_parents=sort!(unique(point.parent for point in points))
    all(parent->parent>0 && haskey(model.discrete,(3,Int(parent))),parents) &&
    all(parent->parent>0 && haskey(model.discrete,(2,Int(parent))),surfaces) &&
    all(parent->parent>0 && haskey(model.discrete,(1,Int(parent))),curve_parents) &&
    all(parent->parent>0 && haskey(model.discrete,(0,Int(parent))),point_parents) || throw(ArgumentError(
        "$caller: lower BRep requires discrete Region, Surface, Curve and Point parents"))
    line_support=Dict{Tuple{UInt64,UInt64},NamedTuple}();curve_edges=Dict{Int32,Vector{NamedTuple}}()
    curve_nodes=Dict{Int32,Set{Int32}}();point_cells=Dict{Int32,NamedTuple}();point_rows=Set{Int32}()
    for line in lines
        key=minmax(line.uids...)
        haskey(line_support,key) && throw(ArgumentError("$caller: original BRep repeats a Line2 support"))
        line_support[key]=line
        push!(get!(()->NamedTuple[],curve_edges,line.parent),(uids=line.uids,rows=line.rows))
        nodes=get!(()->Set{Int32}(),curve_nodes,line.parent);push!(nodes,line.rows[1],line.rows[2])
    end
    for point in points
        haskey(point_cells,point.parent) && throw(ArgumentError(
            "$caller: each declared lower Point needs exactly one actual Point15 cell"))
        point_cells[point.parent]=point
        point.rows[1] in point_rows && throw(ArgumentError("$caller: original BRep repeats a Point15 vertex"))
        push!(point_rows,point.rows[1])
    end
    curve_ends=Dict{Int32,NTuple{2,UInt64}}();used_points=Set{Int32}()
    for parent in curve_parents
        edges=curve_edges[parent];outgoing=Dict{UInt64,Int}();incoming=Dict{UInt64,Int}()
        for (index,edge) in enumerate(edges)
            (haskey(outgoing,edge.uids[1]) || haskey(incoming,edge.uids[2])) && throw(ArgumentError(
                "$caller: declared Curve $parent is not one nonbranching directed chain"))
            outgoing[edge.uids[1]]=index;incoming[edge.uids[2]]=index
        end
        seed=something(findfirst(index->!haskey(incoming,edges[index].uids[1]),eachindex(edges)),1)
        seen=falses(length(edges));index=seed;last=seed
        while index!=0 && !seen[index]
            seen[index]=true;last=index;index=get(outgoing,edges[index].uids[2],0)
        end
        all(seen) || throw(ArgumentError(
            "$caller: declared Curve $parent is not one nonbranching directed chain"))
        boundary=model.discrete[(1,Int(parent))].boundary
        if index==seed
            isempty(boundary) || throw(ArgumentError(
                "$caller: a closed declared Curve $parent cannot have Point boundaries"))
            curve_ends[parent]=(UInt64(0),UInt64(0))
            continue
        end
        curve_ends[parent]=(edges[seed].uids[1],edges[last].uids[2])
        length(boundary)==2 || throw(ArgumentError(
            "$caller: declared Curve $parent needs signed start and end Point boundaries"))
        covered_signs=Set{Int}()
        for (bdim,bsigned) in boundary
            bdim==0 || throw(ArgumentError(
                "$caller: declared Curve $parent bounds require Point entities"))
            point=Int32(abs(bsigned))
            haskey(point_cells,point) || throw(ArgumentError(
                "$caller: declared boundary Point $point has no actual Point15 cell"))
            side=bsigned>0 ? 1 : 2
            side in covered_signs && throw(ArgumentError(
                "$caller: declared Curve $parent repeats a boundary orientation"))
            push!(covered_signs,side);push!(used_points,point)
            node=side==1 ? curve_ends[parent][1] : curve_ends[parent][2]
            point_cells[point].uids[1]==node || throw(ArgumentError(
                "$caller: declared boundary Point $point does not lie on the actual Curve end"))
        end
    end
    surface_borders=_partition_brep_surface_borders(triangles,caller)
    used_curves=Set{Int32}()
    for parent in surfaces
        boundary=model.discrete[(2,Int(parent))].boundary
        isempty(boundary) && continue
        all(entry->entry[1]==1,boundary) && allunique(Int32(abs(entry[2])) for entry in boundary) ||
            throw(ArgumentError(
                "$caller: each declared Surface boundary must list distinct Curve entities"))
        border=surface_borders[parent];covered=Set{Tuple{UInt64,UInt64}}()
        for (_,signed) in boundary
            curve=Int32(abs(signed));haskey(curve_edges,curve) || throw(ArgumentError(
                "$caller: declared boundary Curve $curve has no actual Line2 cells"))
            push!(used_curves,curve)
            for edge in curve_edges[curve]
                key=minmax(edge.uids...);oriented=get(border,key,nothing)
                oriented===nothing && throw(ArgumentError(
                    "$caller: declared Curve $curve cells must lie on the actual Surface border"))
                expected=signed>0 ? edge.uids : (edge.uids[2],edge.uids[1])
                oriented==expected || throw(ArgumentError(
                    "$caller: declared Curve $curve orientation disagrees with its Surface border"))
                key in covered && throw(ArgumentError(
                    "$caller: declared Surface boundary repeats a border Curve cell"))
                push!(covered,key)
            end
        end
        covered==Set(keys(border)) || throw(ArgumentError(
            "$caller: declared lower Surface Curves do not cover its complete border"))
    end
    Set(curve_parents)==used_curves || throw(ArgumentError(
        "$caller: unbound lower Curves are unsupported"))
    Set(point_parents)==used_points || throw(ArgumentError(
        "$caller: unbound lower Points are unsupported"))
    all(row->begin
            dim,tag=class.node_entities[row]
            dim==0 ? (get(point_cells,tag,nothing)!==nothing && point_cells[tag].rows[1]==row) :
            dim==1 ? (haskey(curve_nodes,tag) && row in curve_nodes[tag]) :
            dim==2 ? tag in surfaces : dim==3 && tag in parents
        end,
        unique(vcat([collect(cell.rows) for cell in vcat(triangles,tets,lines,points)]...))) ||
        throw(ArgumentError(
            "$caller: referenced lower BRep vertices have inconsistent entity classification"))
    tet_rows=Matrix{Int32}(undef,4,length(tets))
    for (column,tet) in enumerate(tets);tet_rows[:,column]=collect(tet.rows);end
    validate(Mesh(cached.coords;tets=tet_rows)).ok || throw(ArgumentError(
        "$caller: lower BRep needs positive manifold Tet4 cells"))
    facets=Dict{NTuple{3,UInt64},Vector{NamedTuple}}();facet_order=NTuple{3,UInt64}[]
    for tet in tets,face in _PARTITION_TET_FACETS
        rows=ntuple(i->tet.rows[face[i]],3);labels=ntuple(i->tet.uids[face[i]],3)
        key=_partition_support(labels);haskey(facets,key) || push!(facet_order,key)
        push!(get!(()->NamedTuple[],facets,key),(;rows,uids=labels,tet))
    end
    for hits in values(facets)
        length(hits)<=2 || throw(ArgumentError("$caller: lower BRep facet has too many Tet carriers"))
        length(hits)==2 && !_partition_cycle_reverse(hits[1].uids,hits[2].uids) &&
            throw(ArgumentError("$caller: shared Tet facet has inconsistent winding"))
    end
    triangle_support=Dict{NTuple{3,UInt64},NamedTuple}();lower_edges=Dict{Tuple{UInt64,UInt64},Vector{NamedTuple}}()
    for tri in triangles
        key=_partition_support(tri.uids)
        haskey(triangle_support,key) && throw(ArgumentError("$caller: original BRep repeats a Tri3 support"))
        triangle_support[key]=tri
        hits=get(facets,key,NamedTuple[])
        !isempty(hits) && any(hit->hit.tet.part==tri.part,hits) || throw(ArgumentError(
            "$caller: each lower Tri3 partition must match an actual Tet facet carrier"))
        for (a,b) in _PARTITION_TRI_EDGES
            push!(get!(()->NamedTuple[],lower_edges,minmax(tri.uids[a],tri.uids[b])),
                (uids=(tri.uids[a],tri.uids[b]),rows=(tri.rows[a],tri.rows[b]),tri=tri))
        end
    end
    declared=Dict{Int32,Vector{Int32}}();used=Set{Int32}()
    for parent in parents
        boundary=model.discrete[(3,Int(parent))].boundary
        !isempty(boundary) && all(entry->entry[1]==2 && Int32(abs(entry[2])) in surfaces,boundary) &&
            allunique(abs(entry[2]) for entry in boundary) || throw(ArgumentError(
                "$caller: each Region needs a complete declared Tri3 Surface boundary"))
        expected=Set(key for (key,hits) in facets if
            count(hit->hit.tet.parent==parent,hits)==1)
        covered=Set{NTuple{3,UInt64}}();declared[parent]=Int32[value for (_,value) in boundary]
        for (_,signed) in boundary
            surface=Int32(abs(signed));push!(used,surface)
            for tri in triangles
                tri.parent==surface || continue
                key=_partition_support(tri.uids);key in expected && !(key in covered) ||
                    throw(ArgumentError("$caller: declared Surface does not uniquely cover its Region facets"))
                hit=only(filter(hit->hit.tet.parent==parent,facets[key]))
                oriented=signed>0 ? tri.uids : (tri.uids[1],tri.uids[3],tri.uids[2])
                _partition_cycle_equal(oriented,hit.uids) || throw(ArgumentError(
                    "$caller: declared Surface sign disagrees with its actual Tet facet winding"))
                push!(covered,key)
            end
        end
        covered==expected || throw(ArgumentError("$caller: declared lower BRep omits Region facets"))
    end
    used==Set(surfaces) || throw(ArgumentError("$caller: unbound lower Surfaces are unsupported"))
    return (;entities,blocks,owners,uids,public,triangles,tets,lines,points,
        parents,surfaces,curve_parents,point_parents,curve_edges,point_cells,
        facets,facet_order,triangle_support,lower_edges,declared,surface_borders)
end

function _partition_brep_boundary_certificate(source,children,part_entities,boundaries,
                                               interfaces,curves,caller;topology)
    surface_faces=Dict{Int32,Vector{NamedTuple}}()
    for tri in source.triangles
        child=children[(2,tri.parent,tri.part)]
        push!(get!(()->NamedTuple[],surface_faces,child),tri)
    end
    merge!(surface_faces,interfaces)
    for ((dim,tag),_) in part_entities
        dim==3 || continue
        expected=Dict{NTuple{3,UInt64},NamedTuple}()
        for (key,hits) in source.facets
            selected=[hit for hit in hits if children[(3,hit.tet.parent,hit.tet.part)]==tag]
            length(selected)==1 && (expected[key]=only(selected))
        end
        covered=Set{NTuple{3,UInt64}}()
        for signed in boundaries[(dim,tag)],face in surface_faces[abs(signed)]
            key=_partition_support(face.uids)
            haskey(expected,key) && !(key in covered) || throw(ArgumentError(
                "$caller: child Region boundary repeats or includes a foreign Tet facet"))
            oriented=signed>0 ? face.uids : (face.uids[1],face.uids[3],face.uids[2])
            _partition_cycle_equal(oriented,expected[key].uids) || throw(ArgumentError(
                "$caller: child Region boundary has inconsistent facet winding"))
            push!(covered,key)
        end
        topology && covered!=Set(keys(expected)) && throw(ArgumentError(
            "$caller: child Region topology does not cover its complete Tet boundary"))
    end
    topology || return nothing
    for (tag,faces) in surface_faces
        edges=Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}()
        for face in faces,(a,b) in _PARTITION_TRI_EDGES
            pair=(face.uids[a],face.uids[b]);push!(get!(()->Tuple{UInt64,UInt64}[],edges,minmax(pair...)),pair)
        end
        expected=Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}()
        for (key,hits) in edges
            length(hits)<=2 || throw(ArgumentError("$caller: child Surface has a nonmanifold edge"))
            if length(hits)==2
                hits[1]==reverse(hits[2]) || throw(ArgumentError("$caller: child Surface internal edge winding disagrees"))
            else
                expected[key]=only(hits)
            end
        end
        covered=Set{Tuple{UInt64,UInt64}}()
        for signed in boundaries[(2,Int(tag))],edge in curves[abs(signed)]
            key=minmax(edge.uids...)
            # A declared Curve can span sibling children; only its edges on
            # this child's actual border bind here and must not repeat.
            haskey(expected,key) || continue
            oriented=signed>0 ? edge.uids : reverse(edge.uids)
            expected[key]==oriented && !(key in covered) || throw(ArgumentError(
                "$caller: Surface Curve incidence repeats or disagrees with its actual border"))
            push!(covered,key)
        end
        parent=part_entities[(2,Int(tag))].parent
        inherited=parent[1]==2 ? source.surface_borders[parent[2]] :
            Dict{Tuple{UInt64,UInt64},Tuple{UInt64,UInt64}}()
        all(key->key in covered || get(inherited,key,nothing)==expected[key],keys(expected)) ||
            throw(ArgumentError("$caller: child Surface omits a new seam or changes its original border"))
    end
    return nothing
end

function _partition_brep_chains(edges,caller)
    outgoing=Dict{UInt64,Int}();incoming=Dict{UInt64,Int}()
    for (index,edge) in enumerate(edges)
        haskey(outgoing,edge.uids[1]) && throw(ArgumentError("$caller: branching partition seam is unsupported"))
        outgoing[edge.uids[1]]=index
    end
    for (index,edge) in enumerate(edges)
        haskey(incoming,edge.uids[2]) && throw(ArgumentError("$caller: branching partition seam is unsupported"))
        incoming[edge.uids[2]]=index
    end
    # Open paths begin at literal source edges with no incoming incidence;
    # remaining edges form closed cycles. No coordinate ordering is involved.
    starts=[index for (index,edge) in enumerate(edges) if !haskey(incoming,edge.uids[1])]
    seen=falses(length(edges));chains=NamedTuple[]
    for seed in vcat(starts,collect(eachindex(edges)))
        seen[seed] && continue
        chain=Int[];index=seed
        while !seen[index]
            seen[index]=true;push!(chain,index);next=get(outgoing,edges[index].uids[2],0)
            next==0 && (index=0;break)
            index=next
        end
        closed=index==seed
        (closed || index==0 && seed in starts) || throw(ArgumentError("$caller: partition seam joins another chain"))
        closed && length(chain)<3 && throw(ArgumentError("$caller: partition seam is not a closed polygon"))
        push!(chains,(indices=chain,closed=closed))
    end
    return chains
end

function _partition_brep_production_plan(model,cached,class,count,assignment,caller;
                                         topology::Bool,ghosts::Bool)
    source=_partition_brep_snapshot(model,cached,class,assignment,caller)
    tet_parts=Int32[tet.part for tet in source.tets]
    length(unique(tet_parts))==count || throw(ArgumentError("$caller: every partition must own a Tet4"))
    topology && any(key->length(source.facets[key])==2 &&
        source.facets[key][1].tet.part!=source.facets[key][2].tet.part &&
        haskey(source.triangle_support,key),source.facet_order) && throw(ArgumentError(
        "$caller: a partition interface coincides with an original Tri3; its distinct geometric roles are not implemented"))
    staged=Model._model_mesh_transaction_copy(model)
    part_entities=Dict{Tuple{Int,Int},Elements.MixedEntity}()
    children=Dict{Tuple{Int,Int32,Int32},Int32}()
    interfaces=Dict{Int32,Vector{NamedTuple}}();curves=Dict{Int32,Vector{NamedTuple}}()
    for dim in (2,3)
        cells=dim==2 ? source.triangles : source.tets
        for (parent,part) in sort!(unique((cell.parent,cell.part) for cell in cells))
            tag=Int32(_partition_entity_birth!(staged,dim,caller));children[(dim,parent,part)]=tag
            groups=Model._model_projection_physical_tags(staged,dim,Int(parent))
            staged.entity_physicals[(dim,Int(tag))]=copy(groups)
            for group in groups
                members=staged.physical[(dim,Int(group))];Int(tag) in members || push!(members,Int(tag))
            end
            selected=[cell for cell in cells if cell.parent==parent && cell.part==part]
            rows=unique(vcat([collect(cell.rows) for cell in selected]...))
            part_entities[(dim,Int(tag))]=Elements.MixedEntity(dim,tag,_partition_bounds(cached.coords,rows);
                physical_tags=groups,partitions=Int32[part],parent=(dim,parent))
        end
    end
    curve_children=Dict{Int32,Int32}();point_children=Dict{Int32,Int32}()
    curve_partitions=Dict{Int32,Vector{Int32}}()
    for parent in source.curve_parents
        tag=Int32(_partition_entity_birth!(staged,1,caller));curve_children[parent]=tag
        groups=Model._model_projection_physical_tags(staged,1,Int(parent))
        staged.entity_physicals[(1,Int(tag))]=copy(groups)
        for group in groups
            members=staged.physical[(1,Int(group))];Int(tag) in members || push!(members,Int(tag))
        end
        edges=source.curve_edges[parent]
        parts=sort!(unique(Int32[entry.tri.part for edge in edges
            for entry in source.lower_edges[minmax(edge.uids...)]]))
        curve_partitions[parent]=parts
        rows=unique(reduce(vcat,[collect(edge.rows) for edge in edges]))
        part_entities[(1,Int(tag))]=Elements.MixedEntity(1,tag,_partition_bounds(cached.coords,rows);
            physical_tags=groups,partitions=parts,parent=(1,Int(parent)))
        curves[tag]=edges
    end
    point_declaring=Dict{Int32,Vector{Int32}}()
    for parent in source.curve_parents,(_,signed) in model.discrete[(1,Int(parent))].boundary
        push!(get!(()->Int32[],point_declaring,Int32(abs(signed))),parent)
    end
    for parent in source.point_parents
        tag=Int32(_partition_entity_birth!(staged,0,caller));point_children[parent]=tag
        groups=Model._model_projection_physical_tags(staged,0,Int(parent))
        staged.entity_physicals[(0,Int(tag))]=copy(groups)
        for group in groups
            members=staged.physical[(0,Int(group))];Int(tag) in members || push!(members,Int(tag))
        end
        cell=source.point_cells[parent]
        parts=sort!(unique(reduce(vcat,[curve_partitions[curve] for curve in point_declaring[parent]])))
        part_entities[(0,Int(tag))]=Elements.MixedEntity(0,tag,_partition_bounds(cached.coords,[cell.rows[1]]);
            physical_tags=groups,partitions=parts,parent=(0,Int(parent)))
    end
    flat_owners=Int32[owner for owners in source.owners for owner in owners]
    for cell in source.triangles;flat_owners[cell.index]=children[(2,cell.parent,cell.part)];end
    for cell in source.tets;flat_owners[cell.index]=children[(3,cell.parent,cell.part)];end
    for cell in source.lines;flat_owners[cell.index]=curve_children[cell.parent];end
    for cell in source.points;flat_owners[cell.index]=point_children[cell.parent];end
    blocks=copy(source.blocks);owners=Vector{Int32}[];offset=0
    for block in blocks
        n=size(block.nodes,2);push!(owners,flat_owners[offset+1:offset+n]);offset+=n
    end
    labels=copy(source.uids);maximum_uid=ELEMENT_TAG_MAX[];cell_count=length(source.public.element_tags)
    boundaries=Dict(key=>Int32[] for key in keys(part_entities))
    for tet in source.tets
        child=children[(3,tet.parent,tet.part)]
        for signed in source.declared[tet.parent],tri in source.triangles
            tri.parent==abs(signed) || continue
            hits=source.facets[_partition_support(tri.uids)]
            any(hit->hit.tet.index==tet.index,hits) || continue
            face=children[(2,tri.parent,tri.part)]*sign(signed)
            face in boundaries[(3,Int(child))] || push!(boundaries[(3,Int(child))],face)
        end
    end
    for parent in source.curve_parents
        boundaries[(1,Int(curve_children[parent]))]=Int32[
            sign(signed)*point_children[Int32(abs(signed))]
            for (_,signed) in model.discrete[(1,Int(parent))].boundary]
    end
    seams=Dict{Tuple{Int32,Tuple{Int32,Int32}},Vector{NamedTuple}}()
    seam_order=Tuple{Int32,Tuple{Int32,Int32}}[]
    if topology
        edge_curve=Dict{Tuple{UInt64,UInt64},NamedTuple}()
        for parent in source.curve_parents,edge in source.curve_edges[parent]
            edge_curve[minmax(edge.uids...)]=(uids=edge.uids,tag=curve_children[parent])
        end
        seam_lookup=Dict{Tuple{UInt64,UInt64},NamedTuple}()
        for (key,entry) in edge_curve
            seam_lookup[key]=entry
        end
        grouped=Dict{Tuple{Int32,Int32},Vector{NamedTuple}}()
        for key in source.facet_order
            hits=source.facets[key];length(hits)==2 || continue
            left,right=hits;left.tet.part==right.tet.part && continue
            child=children[(3,left.tet.parent,left.tet.part)];other=children[(3,right.tet.parent,right.tet.part)]
            push!(get!(()->NamedTuple[],grouped,minmax(child,other)),
                (rows=left.rows,uids=left.uids,child=child,other=other))
        end
        for pair in sort!(collect(keys(grouped)))
            faces=grouped[pair]
            for component in _partition_interface_components(faces)
                primary=faces[first(component)].child;tag=Int32(_partition_entity_birth!(staged,2,caller))
                cell_count<=typemax(Int32)-length(component) || throw(ArgumentError("$caller: interface cell capacity exceeded"))
                cell_count+=length(component)
                members=sort!(Int32[part_entities[(3,Int(child))].partitions[1] for child in pair])
                interface_faces=NamedTuple[];nodes=Matrix{Int32}(undef,3,length(component));ids=UInt64[]
                for (column,index) in enumerate(component)
                    face=faces[index];rows=face.child==primary ? face.rows : (face.rows[1],face.rows[3],face.rows[2])
                    uids=Tuple(source.public.node_tags[collect(rows)]);nodes[:,column]=collect(rows)
                    maximum_uid<typemax(UInt64) || throw(ArgumentError("$caller: interface UID capacity exceeded"))
                    maximum_uid+=UInt64(1);push!(ids,maximum_uid);push!(interface_faces,(;rows,uids))
                end
                part_entities[(2,Int(tag))]=Elements.MixedEntity(2,tag,_partition_bounds(cached.coords,unique(vec(nodes)));
                    partitions=members,parent=part_entities[(3,Int(primary))].parent)
                boundaries[(2,Int(tag))]=Int32[];interfaces[tag]=interface_faces
                push!(blocks,ElementBlock(2,nodes));push!(owners,fill(tag,length(ids)));push!(labels,ids)
                for child in pair;push!(boundaries[(3,Int(child))],child==primary ? tag : -tag);end
                edge_hits=Dict{Tuple{UInt64,UInt64},Vector{NamedTuple}}();edge_order=Tuple{UInt64,UInt64}[]
                for face in interface_faces,(a,b) in _PARTITION_TRI_EDGES
                    key=minmax(face.uids[a],face.uids[b]);haskey(edge_hits,key) || push!(edge_order,key)
                    push!(get!(()->NamedTuple[],edge_hits,key),(uids=(face.uids[a],face.uids[b]),rows=(face.rows[a],face.rows[b])))
                end
                for key in edge_order
                    hits=edge_hits[key]
                    length(hits)<=2 || throw(ArgumentError("$caller: nonmanifold interface edge"))
                    if length(hits)==2
                        hits[1].uids==reverse(hits[2].uids) || throw(ArgumentError("$caller: interface edge winding disagrees"))
                        continue
                    end
                    curve_edge=get(edge_curve,key,nothing)
                    if curve_edge!==nothing
                        lower=get(source.lower_edges,key,NamedTuple[])
                        length(lower)==2 &&
                            sort!(unique!(Int32[edge.tri.part for edge in lower]))==members || throw(ArgumentError(
                                "$caller: an interface border on a declared Curve must meet actual lower Tri3 carriers on both partition sides"))
                        seam_lookup[key]=(uids=curve_edge.uids,tag=curve_edge.tag)
                        continue
                    end
                    lower=get(source.lower_edges,key,NamedTuple[])
                    length(lower)==2 &&
                        sort!(Int32[edge.tri.part for edge in lower])==members || throw(ArgumentError(
                            "$caller: interface seam must meet actual lower Tri3 carriers on both partition sides"))
                    reference=argmin([children[(2,edge.tri.parent,edge.tri.part)] for edge in lower])
                    edge=lower[reference];group=(edge.tri.parent,Tuple(members))
                    haskey(seams,group) || push!(seam_order,group)
                    current=get!(()->NamedTuple[],seams,group)
                    any(existing->_partition_support(existing.uids)==_partition_support(edge.uids),current) ||
                        push!(current,(uids=edge.uids,rows=edge.rows))
                end
            end
        end
        for group in seam_order
            edges=seams[group]
            for chain in _partition_brep_chains(edges,caller)
                cycle=chain.indices
                cell_count<=typemax(Int32)-length(cycle) || throw(ArgumentError("$caller: seam cell capacity exceeded"))
                cell_count+=length(cycle)
                tag=Int32(_partition_entity_birth!(staged,1,caller));nodes=Matrix{Int32}(undef,2,length(cycle));ids=UInt64[]
                for (column,index) in enumerate(cycle)
                    edge=edges[index];nodes[:,column]=collect(edge.rows)
                    maximum_uid<typemax(UInt64) || throw(ArgumentError("$caller: seam UID capacity exceeded"))
                    maximum_uid+=UInt64(1);push!(ids,maximum_uid)
                    seam_lookup[minmax(edge.uids...)]=(uids=edge.uids,tag=tag)
                end
                part_entities[(1,Int(tag))]=Elements.MixedEntity(1,tag,_partition_bounds(cached.coords,unique(vec(nodes)));
                    partitions=collect(group[2]),parent=(2,group[1]))
                boundaries[(1,Int(tag))]=Int32[]
                curves[tag]=edges[cycle]
                push!(blocks,ElementBlock(1,nodes));push!(owners,fill(tag,length(ids)));push!(labels,ids)
            end
        end
        for ((dim,tag),entity) in part_entities
            dim==2 || continue
            faces=haskey(interfaces,Int32(tag)) ? interfaces[Int32(tag)] :
                [tri for tri in source.triangles if children[(2,tri.parent,tri.part)]==tag]
            signs=Dict{Int32,Int}()
            for face in faces,(a,b) in _PARTITION_TRI_EDGES
                key=minmax(face.uids[a],face.uids[b]);seam=get(seam_lookup,key,nothing);seam===nothing && continue
                orientation=(face.uids[a],face.uids[b])==seam.uids ? 1 : -1
                haskey(signs,seam.tag) && signs[seam.tag]!=orientation && throw(ArgumentError("$caller: Surface seam has inconsistent orientation"))
                signs[seam.tag]=orientation
            end
            boundaries[(dim,tag)]=Int32[Int32(signs[curve])*curve for curve in sort!(collect(keys(signs)))]
        end
    end
    _partition_brep_boundary_certificate(source,children,part_entities,boundaries,interfaces,curves,caller;topology)
    node_owners=copy(class.node_entities)
    incident=fill((4,Int32(typemax(Int32))),nnodes(cached))
    for (block,classified) in zip(blocks,owners),column in axes(block.nodes,2),row in block.nodes[:,column]
        candidate=(Elements.msh_dimension(block.msh),classified[column])
        isless(candidate,incident[row]) && (incident[row]=candidate)
    end
    for row in eachindex(incident);incident[row][1]<4 && (node_owners[row]=incident[row]);end
    for ((dim,tag),entity) in collect(part_entities)
        boundary=boundaries[(dim,tag)]
        staged.discrete[(dim,tag)].boundary=Tuple{Int,Int}[(dim-1,Int(value)) for value in boundary]
        part_entities[(dim,tag)]=Elements.MixedEntity(dim,tag,entity.bbox;
            physical_tags=entity.physical_tags,boundaries=boundary,partitions=entity.partitions,parent=entity.parent)
    end
    global_entities=copy(source.entities);ghost_entities=Tuple{Int32,Int32}[];ghost_elements=Elements.MixedGhostElement[]
    if ghosts
        node_parts=[Set{Int32}() for _ in 1:nnodes(cached)]
        for tet in source.tets,row in tet.rows;push!(node_parts[row],tet.part);end
        destinations=[Set{Int32}() for _ in source.tets]
        for (index,tet) in enumerate(source.tets),row in tet.rows;union!(destinations[index],node_parts[row]);end
        for (index,tet) in enumerate(source.tets)
            delete!(destinations[index],tet.part);isempty(destinations[index]) && continue
            push!(ghost_elements,Elements.MixedGhostElement(tet.index,tet.part,sort!(collect(destinations[index]))))
        end
        for part in 1:count
            tag=Int32(_partition_entity_birth!(staged,3,caller));push!(ghost_entities,(tag,Int32(part)))
            selected=[tet for (index,tet) in enumerate(source.tets) if Int32(part) in destinations[index]]
            rows=isempty(selected) ? Int32[] : unique(vcat([collect(tet.rows) for tet in selected]...))
            global_entities[(3,Int(tag))]=Elements.MixedEntity(3,tag,_partition_bounds(cached.coords,rows))
        end
    end
    sum(size(block.nodes,2) for block in blocks)<=typemax(Int32) || throw(ArgumentError("$caller: partition cell capacity exceeded"))
    data=Elements._native_mixed_entity_data(global_entities;node_entities=node_owners,
        node_parametric=deepcopy(class.node_parametric),external_node_tags=copy(source.public.node_tags),
        block_entities=owners,external_element_tags=labels)
    partition=Elements.MixedPartitionData(count;entities=part_entities,ghost_entities,ghost_elements)
    target=MixedMesh(cached.coords,blocks;entity_data=data,elementary_entities=owners,
        physical_names=staged.physical_names,partition_data=partition,
        ancillary_sections=cached isa MixedMesh ? cached.ancillary_sections : [],
        data_sections=cached isa MixedMesh ? cached.data_sections : Elements.MshDataSection[])
    validate(target).ok || throw(ArgumentError("$caller: lower BRep partition produced invalid or duplicate cells"))
    _geo_import_entity_receipts!(staged,target)
    identities=Dict(uid=>uid for uid in source.public.node_tags)
    declarations=union(Set(keys(global_entities)),Set(keys(part_entities)))
    _dim01_reconcile_records!(staged,declarations,target,identities;runtime_parameters=data.node_parametric)
    previous=model.meshing.mesh_authority
    ranks=previous!==nothing && previous.mesh===cached ? previous.private_rank : collect(UInt64,1:nnodes(target))
    staged.meshing.mesh_authority=Model._model_mesh_authority(target;caller,private_rank=ranks)
    cell_map=Dict(uid=>uid for uid in source.public.element_tags)
    return (;model=staged,mesh=target,class=_higher_classification_from_data(target,caller),
        authority=Set{Tuple{Int,Int32}}((dim,Int32(tag)) for (dim,tag) in declarations),
        max_node_tag=NODE_TAG_MAX[],max_element_tag=maximum_uid,
        node_tag_map=identities,element_tag_map=cell_map,visibility_tag_map=cell_map)
end
