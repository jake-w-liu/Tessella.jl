# Actual lower Tri3 BRep transport. Interface seams are derived from literal
# node UID incidence, with one authoritative original cell per UID.
const _PARTITION_TET_FACETS=((1,3,2),(1,2,4),(2,3,4),(3,1,4))
const _PARTITION_TRI_EDGES=((1,2),(2,3),(3,1))
_partition_support(uids)=Tuple(sort!(collect(uids)))
_partition_cycle_equal(a,b)=any(a==ntuple(i->b[mod1(i+shift,3)],3) for shift in 0:2)
_partition_cycle_reverse(a,b)=_partition_cycle_equal(a,(b[1],b[3],b[2]))

function _partition_brep_snapshot(model,cached,class,assignment,caller)
    class.public_tags!==nothing || throw(ArgumentError("$caller: lower BRep needs actual UID authority"))
    _high_order_overlay(cached)===nothing || throw(ArgumentError("$caller: lower BRep partitioning requires first-order cells"))
    catalog=[row for row in _cache_catalog(cached,class) if !isempty(row[4])]
    all(row->row[1] in (2,4),catalog) && any(row->row[1]==4,catalog) || throw(ArgumentError(
        "$caller: lower BRep partitioning supports Tri3/Tet4 only, without existing Curve/Point cells"))
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
    triangles=NamedTuple[];tets=NamedTuple[]
    for (msh,dim,offset,cells,classified) in catalog
        physical=Int32[_dim01_physical(model,(dim,Int(owner))) for owner in classified]
        push!(blocks,ElementBlock(msh,cells,physical));push!(owners,copy(classified))
        push!(uids,copy(public.element_tags[offset+1:offset+size(cells,2)]))
        for column in axes(cells,2)
            rows=Tuple(cells[:,column]);labels=Tuple(public.node_tags[collect(rows)]);index=offset+column
            record=(;rows,uids=labels,index,parent=classified[column],part=assignment[index])
            push!(msh==2 ? triangles : tets,record)
        end
    end
    parents=sort!(unique(tet.parent for tet in tets));surfaces=sort!(unique(tri.parent for tri in triangles))
    all(parent->parent>0 && haskey(model.discrete,(3,Int(parent))),parents) &&
    all(parent->parent>0 && haskey(model.discrete,(2,Int(parent))) &&
        isempty(model.discrete[(2,Int(parent))].boundary),surfaces) || throw(ArgumentError(
        "$caller: lower BRep requires discrete Regions and Surfaces without pre-existing Curve boundaries"))
    all(row->class.node_entities[row][1] in (2,3) &&
        class.node_entities[row][2] in (class.node_entities[row][1]==2 ? surfaces : parents),
        unique(vcat([collect(cell.rows) for cell in vcat(triangles,tets)]...))) || throw(ArgumentError(
        "$caller: referenced lower BRep vertices have unsupported Point/Curve classification"))
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
    return (;entities,blocks,owners,uids,public,triangles,tets,parents,surfaces,
        facets,facet_order,triangle_support,lower_edges,declared)
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
            key=minmax(edge.uids...);oriented=signed>0 ? edge.uids : reverse(edge.uids)
            get(expected,key,nothing)==oriented && !(key in covered) || throw(ArgumentError(
                "$caller: Surface Curve incidence repeats or disagrees with its actual border"))
            push!(covered,key)
        end
        covered==Set(keys(expected)) || throw(ArgumentError(
            "$caller: complete Surface boundaries require one supported closed original-Surface seam"))
    end
    return nothing
end

function _partition_brep_cycle(edges,caller)
    outgoing=Dict{UInt64,Int}();incoming=Dict{UInt64,Int}()
    for (index,edge) in enumerate(edges)
        haskey(outgoing,edge.uids[1]) && throw(ArgumentError("$caller: branching partition seam is unsupported"))
        outgoing[edge.uids[1]]=index
    end
    for (index,edge) in enumerate(edges)
        haskey(incoming,edge.uids[2]) && throw(ArgumentError("$caller: branching partition seam is unsupported"))
        incoming[edge.uids[2]]=index
    end
    Set(keys(outgoing))==Set(keys(incoming)) || throw(ArgumentError("$caller: open partition seam is unsupported"))
    seen=falses(length(edges));cycles=Vector{Int}[]
    for seed in eachindex(edges)
        seen[seed] && continue
        cycle=Int[];index=seed
        while !seen[index]
            seen[index]=true;push!(cycle,index);index=outgoing[edges[index].uids[2]]
        end
        index==seed || throw(ArgumentError("$caller: partition seam joins another cycle"))
        length(cycle)>=3 || throw(ArgumentError("$caller: partition seam is not a closed polygon"))
        push!(cycles,cycle)
    end
    return cycles
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
    flat_owners=Int32[owner for owners in source.owners for owner in owners]
    for cell in source.triangles;flat_owners[cell.index]=children[(2,cell.parent,cell.part)];end
    for cell in source.tets;flat_owners[cell.index]=children[(3,cell.parent,cell.part)];end
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
    interfaces=Dict{Int32,Vector{NamedTuple}}();curves=Dict{Int32,Vector{NamedTuple}}()
    seams=Dict{Tuple{Int32,Tuple{Int32,Int32}},Vector{NamedTuple}}()
    if topology
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
                    lower=get(source.lower_edges,key,NamedTuple[])
                    length(lower)==2 && lower[1].tri.parent==lower[2].tri.parent &&
                        sort!(Int32[edge.tri.part for edge in lower])==members || throw(ArgumentError(
                            "$caller: interface seam must meet one original Surface on both partition sides"))
                    reference=argmin([children[(2,edge.tri.parent,edge.tri.part)] for edge in lower])
                    edge=lower[reference];group=(edge.tri.parent,Tuple(members))
                    current=get!(()->NamedTuple[],seams,group)
                    any(existing->_partition_support(existing.uids)==_partition_support(edge.uids),current) ||
                        push!(current,(uids=edge.uids,rows=edge.rows))
                end
            end
        end
        seam_lookup=Dict{Tuple{UInt64,UInt64},NamedTuple}()
        for group in sort!(collect(keys(seams)))
            edges=seams[group]
            for cycle in _partition_brep_cycle(edges,caller)
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
