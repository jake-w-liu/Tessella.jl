# Higher generation uses the same detached identity assembly as dimensions 0/1.
# SetOrder1 and every allocation complete before any session state is published.
include("APIGenerateNoNewIdentity.jl")
function _higher_native_parameters!(assembly,m)
    for position in eachindex(assembly.nodes)
        node=assembly.nodes[position]
        node.parameter===nothing || continue
        dim,tag=node.owner
        native=dim==0 ? haskey(m.points,Int(tag)) : dim==1 ? haskey(m.curves,Int(tag)) :
               dim==2 ? haskey(m.surfaces,Int(tag)) : false
        native || continue
        parameter=dim==0 ? Float64[] :
            Model.model_parametrization(m,dim,Int(tag),collect(node.point))
        assembly.nodes[position]=_Dim01Node(node.tag,node.point,node.owner,parameter)
    end
    return nothing
end

function _higher_node_ranks(assembly,cached,mapping,previous)
    ranks=Dict{UInt64,UInt64}();next=UInt64(0)
    if cached===nothing && previous!==nothing && previous.mesh!==nothing
        # A CAD-only addition can invalidate the public cache while keeping
        # its exact attached products. Their identities retain vertex ranks.
        next=maximum(previous.private_rank;init=UInt64(0))
        for (tag,rank) in zip(previous.mesh.entity_data.external_node_tags,previous.private_rank)
            haskey(assembly.node_positions,tag) || continue
            ranks[tag]=rank;next=max(next,rank)
        end
    end
    if cached!==nothing
        inherited=previous!==nothing && previous.mesh!==nothing &&
                  length(previous.private_rank)==length(mapping)
        for (row,position) in enumerate(mapping)
            rank=inherited ? previous.private_rank[row] : UInt64(row)
            ranks[assembly.nodes[position].tag]=rank
            next=max(next,rank)
        end
    end
    for node in assembly.nodes
        haskey(ranks,node.tag) && continue
        next<typemax(UInt64) || throw(ArgumentError("$(assembly.caller): private vertex rank overflows UInt64"))
        next+=UInt64(1);ranks[node.tag]=next
    end
    return ranks
end

function _higher_complete_ranks!(ranks,assembly)
    next=maximum(values(ranks);init=UInt64(0))
    for node in assembly.nodes
        haskey(ranks,node.tag) && continue
        next<typemax(UInt64) || throw(ArgumentError("$(assembly.caller): private vertex rank overflows UInt64"))
        next+=UInt64(1);ranks[node.tag]=next
    end
    return nothing
end

function _higher_point_incidence_nodes(assembly,m)
    degrees=Dict{Int32,Dict{Int32,Int}}()
    for cell in assembly.cells
        cell.owner[1]==1 && haskey(m.curves,Int(cell.owner[2])) && cell.msh==1 || continue
        counts=get!(()->Dict{Int32,Int}(),degrees,cell.owner[2])
        for node in cell.nodes
            counts[node]=get(counts,node,0)+1
        end
    end
    bindings=Dict{Int32,Int32}()
    for (curve,counts) in degrees
        endpoints=Int32[node for (node,count) in counts if count==1]
        length(endpoints)==2 || continue
        for (endpoint,point_tag) in enumerate(m.curves[Int(curve)])
            # Point ownership is the producer identity even when distinct CAD
            # Points coincide or a large world offset dwarfs the edge length.
            owner=(0,Int32(point_tag))
            matches=Int32[node for node in endpoints if assembly.nodes[node].owner==owner]
            if isempty(matches)
                # SetOrderN may reclassify an endpoint onto its incident Curve.
                # Its stored native parameter still identifies the endpoint.
                lower,upper=Model.model_parametrization_bounds(m,1,Int(curve))
                parameter=endpoint==1 ? only(lower) : only(upper)
                matches=Int32[node for node in endpoints if
                    assembly.nodes[node].owner==(1,curve) &&
                    assembly.nodes[node].parameter!==nothing &&
                    length(assembly.nodes[node].parameter)==1 &&
                    only(assembly.nodes[node].parameter)==parameter]
            end
            if isempty(matches)
                point=m.points[point_tag]
                matches=Int32[node for node in endpoints if assembly.nodes[node].point==point]
            end
            length(matches)==1 || throw(ArgumentError(
                "$(assembly.caller): native Curve[$curve] has ambiguous actual endpoint for Point[$point_tag]"))
            node=only(matches);previous=get(bindings,Int32(point_tag),node)
            previous==node || throw(ArgumentError(
                "$(assembly.caller): incident native Curves use distinct vertices for Point[$point_tag]"))
            bindings[Int32(point_tag)]=node
        end
    end
    return bindings
end

function _higher_assembly_authority(assembly,m,done;ranks=nothing,carrier=nothing,
                                    carrier_rows=Int32[],old_supports=nothing)
    # SetOrder1 deletes old native high-order vertices before extrusion reads
    # the Source. Keep caller primary/orphan rows and every still-used vertex,
    # including supports belonging to an untouched fully discrete cell.
    used=falses(length(assembly.nodes))
    for cell in assembly.cells,node in cell.nodes
        used[node]=true
    end
    kept=Int32[row for row in eachindex(assembly.nodes) if old_supports===nothing ||
        row>length(old_supports) || !old_supports[row] || used[row]]
    mapping=zeros(Int32,length(assembly.nodes))
    coords=Matrix{Float64}(undef,3,length(kept))
    for (position,row) in enumerate(kept)
        mapping[row]=Int32(position)
        node=assembly.nodes[row]
        for axis in 1:3
            coords[axis,position]=node.point[axis]
        end
    end
    types=Int[];groups=Dict{Int,Vector{Int}}()
    for (position,cell) in enumerate(assembly.cells)
        haskey(groups,cell.msh) || push!(types,cell.msh)
        push!(get!(()->Int[],groups,cell.msh),position)
    end
    blocks=ElementBlock[];owners=Vector{Int32}[];tags=Vector{UInt64}[]
    for msh in types
        selected=groups[msh]
        nodes=Matrix{Int32}(undef,Elements.msh_num_nodes(msh),length(selected))
        for (column,position) in enumerate(selected)
            nodes[:,column]=mapping[assembly.cells[position].nodes]
        end
        push!(blocks,ElementBlock(msh,nodes,Int32[assembly.cells[i].physical for i in selected]))
        push!(owners,Int32[assembly.cells[i].owner[2] for i in selected])
        push!(tags,UInt64[assembly.cells[i].tag for i in selected])
    end
    metadata,_=_dim01_entities(m,model_entities(m),assembly,kept)
    data=Elements._native_mixed_entity_data(metadata;
        node_entities=[assembly.nodes[row].owner for row in kept],
        node_parametric=Union{Nothing,Vector{Float64}}[
            Model._model_mesh_projected_parameter(assembly.nodes[row].owner,
                assembly.nodes[row].parameter) for row in kept],
        external_node_tags=[assembly.nodes[row].tag for row in kept],
        block_entities=owners,external_element_tags=tags)
    mesh=MixedMesh(coords,blocks;entity_data=data,elementary_entities=owners)
    rank=ranks===nothing ? UInt64.(kept) :
        UInt64[ranks[assembly.nodes[row].tag] for row in kept]
    return carrier===nothing ? Model._model_mesh_authority(mesh,done;
        caller=assembly.caller,private_rank=rank) : Model._model_mesh_authority(mesh,done;
        caller=assembly.caller,private_rank=rank,carrier,carrier_rows)
end

function _higher_append_surface_source_nodes!(assembly,mesh,production)
    for (tag,point) in production.source_node_updates
        position=get(assembly.node_positions,tag,Int32(0))
        position>0 || throw(ArgumentError(
            "$(assembly.caller): restored Curve identity $tag is absent from the lower producer"))
        previous=assembly.nodes[position]
        previous.owner[1]<2 && all(isfinite,point) || throw(ArgumentError(
            "$(assembly.caller): restored Curve identity $tag has invalid ownership or coordinates"))
        assembly.nodes[position]=_Dim01Node(tag,point,previous.owner,previous.parameter)
    end
    count=nnodes(mesh)
    length(production.node_tags)==length(production.node_owners)==
        length(production.node_parameters)==length(production.retained_nodes)==count ||
        throw(ArgumentError("$(assembly.caller): actual Surface producer rows do not align"))
    allunique(production.node_tags) && all(>(0),production.node_tags) ||
        throw(ArgumentError("$(assembly.caller): actual Surface producer labels must be unique and positive"))
    previous_maximum=assembly.node_max
    production.node_high_water>=previous_maximum &&
        production.node_high_water>=maximum(production.node_tags;init=UInt64(0)) ||
        throw(ArgumentError("$(assembly.caller): actual Surface producer lost consumed node high water"))
    births=production.birth_rows
    length(births)==count && allunique(births) && all(row->1<=row<=count,births) ||
        throw(ArgumentError("$(assembly.caller): actual Surface birth rows are not a complete permutation"))
    mapping=Vector{Int32}(undef,count)
    for position in births
        tag=production.node_tags[position]
        retained=production.retained_nodes[position]
        previous=get(assembly.node_positions,tag,Int32(0))
        if retained
            previous!=0 || throw(ArgumentError(
                "$(assembly.caller): actual Surface source node $tag is absent from the lower producer"))
        else
            previous==0 && tag>previous_maximum || throw(ArgumentError(
                "$(assembly.caller): fresh Surface node $tag collides with consumed source identity"))
        end
        point=ntuple(axis->mesh.coords[axis,position],3)
        owner=production.node_owners[position]
        if retained
            old=assembly.nodes[previous]
            isequal(old.point,point) && old.owner==owner || throw(ArgumentError(
                "$(assembly.caller): actual Surface source node $tag changed coordinates or ownership"))
        end
        mapping[position]=_dim01_add_node!(assembly,tag,point,owner,
            production.node_parameters[position])
    end
    assembly.node_max=production.node_high_water
    return mapping
end

# A native transfinite Region starts with the first Surface's unsigned edge
# loop, not the smallest Point. Keep the carrier's tensor storage, and bind
# the native tensor frame through its actual Point owners and integer rows.
function _higher_transfinite_surface_frame(m,surface,caller)
    loops=m.surfaces[surface]
    length(loops)==1 || return nothing
    curves=copy(m.loops[only(loops)])
    length(curves)==4 || return nothing
    specification=m.meshing.transfinite_surfaces[surface]
    corners=if isempty(specification.corners)
        geometry=get(m.surface_geometry,surface,nothing)
        # Materialized OCC Box planes retain their trimmed-plane pcurves and
        # bounds through affine transforms. Their native wire begins one edge
        # before the declared CAD cycle; signs still own CAD orientation.
        factory_plane=geometry!==nothing && (get(geometry,:primitive,nothing)===:box ||
            (get(geometry,:occ,nothing)===:plane && haskey(geometry,:pcurves) && haskey(geometry,:uvb)))
        if factory_plane
            curves=vcat(curves[end:end],curves[1:end-1])
        end
        pending=collect(2:length(curves));a,b=m.curves[abs(first(curves))]
        frame=Int[a];endpoint=b
        while !isempty(pending)
            index=Model._chained_tri_next(m,curves,pending,endpoint)
            index===nothing && throw(ArgumentError("$caller: native transfinite face is not a closed edge loop"))
            slot=pending[index];deleteat!(pending,index)
            x,y=m.curves[abs(curves[slot])]
            push!(frame,endpoint);endpoint=x==endpoint ? y : x
        end
        endpoint==a || throw(ArgumentError("$caller: native transfinite face is not a closed edge loop"))
        frame
    else
        collect(specification.corners)
    end
    length(corners)==4 || return nothing
    return (corners=Tuple(corners),curves=curves)
end

function _higher_transfinite_hex_frame(m,tag,caller)
    stored=m.meshing.transfinite_volumes[tag]
    isempty(stored) || return Tuple(stored)
    surfaces=Model._model_volume_boundary_surfaces(m,tag,caller)
    length(surfaces)==6 || return nothing
    face=_higher_transfinite_surface_frame(m,abs(first(surfaces)),caller)
    face===nothing && return nothing
    corners=face.corners;curves=face.curves
    face_edges=Set(abs.(curves));remaining=Int[]
    for signed_surface in surfaces,loop in m.surfaces[abs(signed_surface)],curve in m.loops[loop]
        edge=abs(curve)
        edge in face_edges || edge in remaining || push!(remaining,edge)
    end
    opposite=Int[]
    for point in corners
        choices=Int[]
        for edge in remaining
            a,b=m.curves[edge]
            a==point ? push!(choices,b) : b==point ? push!(choices,a) : nothing
        end
        length(choices)==1 || throw(ArgumentError("$caller: native transfinite Region has ambiguous corner incidence"))
        push!(opposite,only(choices))
    end
    return Tuple(vcat(collect(corners),opposite))
end

function _higher_transfinite_quad_rows(m,tag,mesh,point_rows,caller)
    mesh isa MixedMesh && haskey(m.meshing.transfinite_surfaces,tag) &&
        Model._extrude_gate(m,2,tag)===nothing || return nothing
    geometry=get(m.surface_geometry,tag,nothing)
    geometry!==nothing && get(geometry,:occ,nothing)===:plane &&
        haskey(geometry,:pcurves) && haskey(geometry,:uvb) || return nothing
    catalog=_cache_native_blocks(mesh);surfaces=findall(entry->msh_dimension(entry[1])==2,catalog)
    length(surfaces)==1 || return nothing
    block=only(surfaces);msh,cells,_=catalog[block]
    msh==3 && !isempty(cells) || return nothing
    face=_higher_transfinite_surface_frame(m,tag,caller)
    face===nothing && return nothing
    frame=face.corners
    all(point->haskey(point_rows,point),frame) ||
        throw(ArgumentError("$caller: native transfinite face corners lack actual Point rows"))
    corners=sort!(Int[point_rows[point] for point in frame]);npu=corners[2]
    nnodes(mesh)%npu==0 || throw(ArgumentError("$caller: native transfinite face is not a tensor grid"))
    lattice(row)=((row-1)%npu,(row-1)÷npu)
    start=lattice(point_rows[frame[1]])
    ends=(lattice(point_rows[frame[2]]),lattice(point_rows[frame[4]]))
    differences=ntuple(axis->ntuple(k->ends[axis][k]-start[k],2),2)
    all(delta->count(!iszero,delta)==1,differences) ||
        throw(ArgumentError("$caller: native transfinite face frame is not a tensor symmetry"))
    counts=ntuple(axis->sum(abs,differences[axis]),2)
    directions=ntuple(axis->ntuple(k->sign(differences[axis][k]),2),2)
    row_at(i,j)=Int32(1+sum((start[axis]+i*directions[1][axis]+j*directions[2][axis])*
        (axis==1 ? 1 : npu) for axis in 1:2))
    rows=Int32[row_at(i,j) for i in 0:counts[1] for j in 0:counts[2]]
    allunique(rows) && length(rows)==nnodes(mesh) ||
        throw(ArgumentError("$caller: native transfinite face rows are not bijective"))
    source=Dict(Tuple(sort!(collect(cells[:,column])))=>column for column in axes(cells,2))
    nodes=Matrix{Int32}(undef,4,size(cells,2));columns=Int[]
    for i in 0:counts[1]-1,j in 0:counts[2]-1
        cell=Int32[row_at(i,j),row_at(i+1,j),row_at(i+1,j+1),row_at(i,j+1)]
        column=get(source,Tuple(sort(copy(cell))),0)
        column>0 || throw(ArgumentError("$caller: native transfinite face loses its carrier parent"))
        push!(columns,column)
        first=findfirst(==(cell[1]),cells[:,column])
        forward=all(k->cell[k]==cells[mod1(first+k-1,4),column],1:4)
        backward=all(k->cell[k]==cells[mod1(first-k+1,4),column],1:4)
        forward || backward || throw(ArgumentError("$caller: native transfinite face has inconsistent winding"))
        nodes[:,length(columns)]=forward ? cell : cell[[1,4,3,2]]
    end
    length(columns)==size(cells,2) && allunique(columns) ||
        throw(ArgumentError("$caller: native transfinite face cells are not bijective"))
    return (block=block,nodes=nodes,columns=columns,rows=rows)
end

function _higher_transfinite_hex_rows(m,tag,mesh,point_rows,caller)
    mesh isa MixedMesh && haskey(m.meshing.transfinite_volumes,tag) &&
        !(tag in m.meshing.quad_tri) && Model._extrude_gate(m,3,tag)===nothing || return nothing
    catalog=_cache_native_blocks(mesh)
    volumes=findall(entry->msh_dimension(entry[1])==3,catalog)
    length(volumes)==1 || return nothing
    block=only(volumes);msh,cells,_=catalog[block]
    msh==5 && !isempty(cells) || return nothing
    frame=_higher_transfinite_hex_frame(m,tag,caller)
    frame===nothing && return nothing
    length(frame)==8 && all(point->haskey(point_rows,point),frame) ||
        throw(ArgumentError("$caller: native transfinite corners lack actual Point rows"))
    corner_rows=sort!(Int[point_rows[point] for point in frame])
    npu=corner_rows[2];layer=corner_rows[4]
    layer%npu==0 && nnodes(mesh)%layer==0 ||
        throw(ArgumentError("$caller: native transfinite carrier is not a tensor grid"))
    npv=layer÷npu;npw=nnodes(mesh)÷layer
    lattice(row)=((row-1)%npu,((row-1)÷npu)%npv,(row-1)÷layer)
    start=lattice(point_rows[frame[1]])
    ends=(lattice(point_rows[frame[2]]),lattice(point_rows[frame[4]]),lattice(point_rows[frame[5]]))
    differences=ntuple(axis->ntuple(k->ends[axis][k]-start[k],3),3)
    all(delta->count(!iszero,delta)==1,differences) ||
        throw(ArgumentError("$caller: native transfinite corner frame is not a tensor symmetry"))
    counts=ntuple(axis->sum(abs,differences[axis]),3)
    directions=ntuple(axis->ntuple(k->sign(differences[axis][k]),3),3)
    row_at(i,j,k)=Int32(1+sum((start[axis]+i*directions[1][axis]+j*directions[2][axis]+
        k*directions[3][axis])*(axis==1 ? 1 : axis==2 ? npu : layer) for axis in 1:3))
    rows=Int32[row_at(i,j,k) for i in 0:counts[1] for j in 0:counts[2] for k in 0:counts[3]]
    allunique(rows) && length(rows)==nnodes(mesh) ||
        throw(ArgumentError("$caller: native transfinite tensor rows are not bijective"))
    source=Dict(Tuple(sort!(collect(cells[:,column])))=>column for column in axes(cells,2))
    nodes=Matrix{Int32}(undef,8,size(cells,2));columns=Int[]
    reverse=Model.orient3(m.points[frame[1]],m.points[frame[2]],m.points[frame[4]],m.points[frame[5]])>0
    for i in 0:counts[1]-1,j in 0:counts[2]-1,k in 0:counts[3]-1
        cell=Int32[row_at(i,j,k),row_at(i+1,j,k),row_at(i+1,j+1,k),row_at(i,j+1,k),
            row_at(i,j,k+1),row_at(i+1,j,k+1),row_at(i+1,j+1,k+1),row_at(i,j+1,k+1)]
        column=get(source,Tuple(sort(copy(cell))),0)
        column>0 || throw(ArgumentError("$caller: native transfinite cell loses its carrier parent"))
        push!(columns,column)
        # MHexahedron::reverse, applied by GModel::setAllVolumesPositive.
        nodes[:,length(columns)]=reverse ? cell[[3,2,1,4,7,6,5,8]] : cell
    end
    length(columns)==size(cells,2) && allunique(columns) ||
        throw(ArgumentError("$caller: native transfinite cells are not bijective"))
    get(m.meshing.reverse,(3,tag),false) && (nodes=nodes[[3,2,1,4,7,6,5,8],:])
    return (block=block,nodes=nodes,columns=columns,rows=rows)
end

include("APIHybridBatchPublication.jl")

function _higher_append_part!(assembly,m,mesh,class,dimension,tag,discarded;
                              emitted=Set{Tuple{Int,Int32}}(),source_nodes=length(assembly.nodes),
                              allocate_cells=true,geometry_cells=nothing,source_binding=nothing,
                              source_emission=nothing,
                              source_row_tags=nothing,source_row_retained=nothing,
                              source_rows_existing_only=false,source_node_base=nothing,
                              source_node_high_water=nothing,source_cell_tags=nothing,
                              source_cell_high_water=nothing,native_orphan_births=nothing,
                              source_batch_receipt=nothing,
                              source_row_order=nothing,source_volume_production=nothing,
                              emit_cells=true)
    class===nothing && throw(ArgumentError(
        "$(assembly.caller): generated entity ($dimension, $tag) has no ownership metadata"))
    if source_row_order!==nothing
        production=source_volume_production
        dimension==3 && production isa Model._ModelNativeHybridVolumeProduction &&
            production.mesh===mesh && source_emission===nothing &&
            source_row_tags==production.boundary_node_tags && source_cell_tags==production.element_tags &&
            source_node_high_water==production.node_high_water &&
            source_cell_high_water==production.element_high_water &&
            source_row_order==production.registration_rows && length(source_row_order)==nnodes(mesh) &&
            allunique(source_row_order) && all(row->0<row<=nnodes(mesh),source_row_order) ||
            throw(ArgumentError("$(assembly.caller): foreign or incomplete hybrid producer registration receipt"))
        Model._model_native_hybrid_volume_certificate(mesh,source_row_tags,production.constructor,assembly.caller)
    end
    if source_batch_receipt!==nothing
        source_binding===nothing && source_emission===nothing &&
            source_row_retained===nothing && !source_rows_existing_only &&
            allocate_cells && emit_cells || throw(ArgumentError(
            "$(assembly.caller): hybrid batch conflicts with another producer publication"))
        _higher_native_hybrid_batch_binding(source_batch_receipt,assembly,m,mesh,class,dimension,tag,
            source_row_tags,source_node_base,source_node_high_water,source_cell_tags,
            source_cell_high_water,assembly.caller;source_row_order)
    end
    data=mesh isa MixedMesh ? mesh.entity_data : nothing
    if source_row_tags!==nothing
        length(source_row_tags)==nnodes(mesh) || throw(ArgumentError(
            "$(assembly.caller): volume vertex births do not align with producer rows"))
        source_row_retained===nothing || length(source_row_retained)==nnodes(mesh) ||
            throw(ArgumentError("$(assembly.caller): retained volume mask does not cover producer rows"))
        known=source_row_retained===nothing ? filter(!iszero,source_row_tags) :
            source_row_tags[source_row_retained]
        allunique(known) || throw(ArgumentError(
            "$(assembly.caller): volume vertex births are not unique"))
        if source_rows_existing_only
            source_row_retained!==nothing && all(tag->haskey(assembly.node_positions,tag),known) ||
                throw(ArgumentError("$(assembly.caller): retained volume boundary has no actual source identity"))
            source_node_high_water===nothing || source_node_high_water==assembly.node_max ||
                throw(ArgumentError("$(assembly.caller): retained volume boundary changed its allocator clock"))
        elseif source_node_high_water!==nothing
            source_node_base!==nothing && source_node_base<=source_node_high_water &&
                maximum(known;init=UInt64(0))<=source_node_high_water || throw(ArgumentError(
                "$(assembly.caller): volume vertex high-water mark is inconsistent"))
            assembly.node_max=max(assembly.node_max,source_node_high_water)
        end
    end
    point_rows=Dict(Int(owner[2])=>Int32(row) for (row,owner) in enumerate(class.node_entities) if owner[1]==0)
    native=dimension==3 ? _higher_transfinite_hex_rows(m,tag,mesh,point_rows,assembly.caller) :
        dimension==2 ? _higher_transfinite_quad_rows(m,tag,mesh,point_rows,assembly.caller) : nothing
    if source_emission!==nothing
        native=source_emission.block==0 && native!==nothing ?
            merge(native,(rows=source_emission.rows,)) : source_emission
    end
    mapping=if source_binding===nothing
        point_incidence=_higher_point_incidence_nodes(assembly,m)
        candidates=Dict{Tuple{Int,Int32},Vector{Int32}}()
        for (position,node) in enumerate(assembly.nodes)
            position<=source_nodes && node.owner in discarded && continue
            push!(get!(()->Int32[],candidates,node.owner),Int32(position))
        end
        points=Dict{Tuple{Tuple{Int,Int32},NTuple{3,Float64}},Vector{Int32}}()
        ranges=Dict{Tuple{Int,Int32},Tuple{Int,Vector{Float64},Vector{Int32}}}()
        for (owner,positions) in candidates
            for position in positions
                point=assembly.nodes[position].point
                push!(get!(()->Int32[],points,(owner,point)),position)
            end
            spans=ntuple(axis->maximum(assembly.nodes[i].point[axis] for i in positions)-
                               minimum(assembly.nodes[i].point[axis] for i in positions),3)
            selected_axis=argmax(spans)
            ordered=sort(positions;by=i->assembly.nodes[i].point[selected_axis])
            ranges[owner]=(selected_axis,
                Float64[assembly.nodes[i].point[selected_axis] for i in ordered],ordered)
        end
        mapping=Vector{Int32}(undef,nnodes(mesh))
        for position in (native===nothing ?
            source_row_order===nothing ? axes(mesh.coords,2) : source_row_order : native.rows)
            owner=class.node_entities[position]
            point=ntuple(axis->mesh.coords[axis,position],3)
            parameter=data===nothing ? nothing : data.node_parametric[position]
            if source_row_tags!==nothing && (source_row_retained===nothing ?
                source_row_tags[position]!=0 : source_row_retained[position])
                source_tag=source_row_tags[position]
                actual_row=get(assembly.node_positions,source_tag,Int32(0))
                if actual_row==0
                    source_rows_existing_only && throw(ArgumentError(
                        "$(assembly.caller): retained volume boundary vertex $source_tag is absent"))
                    source_node_base!==nothing && source_node_high_water!==nothing &&
                        source_node_base<source_tag<=source_node_high_water || throw(ArgumentError(
                        "$(assembly.caller): volume vertex $source_tag has no actual source or birth receipt"))
                    mapping[position]=_dim01_add_node!(assembly,source_tag,point,owner,parameter)
                    continue
                end
                actual_node=assembly.nodes[actual_row]
                (source_rows_existing_only ? isequal(actual_node.point,point) &&
                    actual_node.owner==owner : actual_node.point==point) || throw(ArgumentError(
                    "$(assembly.caller): volume boundary vertex $source_tag changed coordinates or ownership"))
                mapping[position]=actual_row
                continue
            end
            if owner[1]==0 && haskey(point_incidence,owner[2])
                node=point_incidence[owner[2]]
                endpoint_node=assembly.nodes[node]
                tolerance=64eps(Float64)*max(1.0,maximum(abs,point))
                maximum(abs(endpoint_node.point[axis]-point[axis]) for axis in 1:3)<=tolerance ||
                    throw(ArgumentError("$(assembly.caller): actual Curve endpoint disagrees with generated Point $owner"))
                # The real lower Curve vertex owns this incidence. Its identity
                # survives even if the last P2 prune classified it onto a Curve.
                mapping[position]=node
                continue
            end
            if parameter===nothing && owner[1] in (1,2)
                parameter=Model.model_parametrization(m,owner[1],Int(owner[2]),collect(point))
            end
            matches=Int32[]
            if owner[1]<dimension || owner in emitted
                tolerance=64eps(Float64)*max(1.0,maximum(abs,point))
                append!(matches,get(points,(owner,point),Int32[]))
                if isempty(matches) && haskey(ranges,owner)
                    search_axis,values,search_positions=ranges[owner]
                    lo=searchsortedfirst(values,prevfloat(point[search_axis]-tolerance))
                    hi=searchsortedlast(values,nextfloat(point[search_axis]+tolerance))
                    for slot in lo:hi
                        candidate=search_positions[slot];old=assembly.nodes[candidate]
                        maximum(abs(old.point[a]-point[a]) for a in 1:3)<=tolerance || continue
                        push!(matches,candidate)
                    end
                end
            end
            length(matches)<=1 || throw(ArgumentError(
                "$(assembly.caller): generated boundary on entity $owner has ambiguous source identity"))
            if isempty(matches)
                node=_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),point,owner,parameter)
                push!(get!(()->Int32[],candidates,owner),node)
                push!(get!(()->Int32[],points,(owner,point)),node)
                mapping[position]=node
            else
                node=only(matches);matched_node=assembly.nodes[node]
                mapping[position]=_dim01_add_node!(assembly,matched_node.tag,matched_node.point,owner,parameter)
            end
        end
        mapping
    else
        _higher_append_surface_source_nodes!(assembly,mesh,source_binding)
    end
    if native_orphan_births!==nothing
        used=falses(nnodes(mesh))
        for (_,cells,_) in _cache_native_blocks(mesh),row in cells
            used[row]=true
        end
        for row in eachindex(used)
            !used[row] && class.node_entities[row]==(dimension,Int32(tag)) || continue
            source_row_retained!==nothing && source_row_retained[row] && continue
            push!(native_orphan_births,assembly.nodes[mapping[row]].tag)
        end
    end
    if !emit_cells
        source_binding===nothing && source_cell_tags===nothing || throw(ArgumentError(
            "$(assembly.caller): a node-only stage cannot consume another producer's cell births"))
        push!(emitted,(dimension,Int32(tag)))
        return mapping
    end
    traces=geometry_cells===nothing ? nothing : [Union{Nothing,_Dim01Cell}[
        nothing for _ in axes(cells,2)] for (_,cells,_) in _cache_native_blocks(mesh)]
    source_binding===nothing || source_cell_tags===nothing || throw(ArgumentError(
        "$(assembly.caller): generated cells have conflicting birth providers"))
    cell_births=source_binding!==nothing ? source_binding.element_tags :
        source_cell_tags===nothing ? UInt64[] : source_cell_tags
    cell_high_water=source_binding!==nothing ? source_binding.element_high_water : source_cell_high_water
    if !isempty(cell_births)
        length(cell_births)==sum(size(cells,2) for (_,dim,_,cells,_) in _cache_catalog(mesh,class)
            if dim==dimension) && allunique(cell_births) || throw(ArgumentError(
            "$(assembly.caller): native Surface cell births do not align with emitted cells"))
        cell_high_water!==nothing && cell_high_water>=maximum(cell_births) &&
            (source_batch_receipt!==nothing || minimum(cell_births)>assembly.cell_max) || throw(ArgumentError(
            "$(assembly.caller): native cell births collide with retained identities"))
        source_batch_receipt===nothing || _higher_native_hybrid_batch_cell_births(
            source_batch_receipt,assembly,tag,cell_births,cell_high_water,assembly.caller)
        assembly.cell_max=cell_high_water
    end
    # An exporter can discard every newly constructed Face cell. These
    # deleted objects still consumed their actual native constructor labels.
    if source_binding!==nothing && cell_high_water!==nothing
        assembly.cell_max=max(assembly.cell_max,cell_high_water)
    end
    birth_offset=0
    for (bi,(msh,dim,_,cells,owners)) in enumerate(_cache_catalog(mesh,class))
        dim==dimension || continue
        physical=_cache_native_blocks(mesh)[bi][3]
        for column in axes(cells,2)
            source_column=native!==nothing && bi==native.block ? native.columns[column] : column
            connectivity=native!==nothing && bi==native.block ? native.nodes[:,column] : cells[:,column]
            legacy_physical=data===nothing ? _dim01_physical(m,(dim,Int(owners[source_column]))) :
                            physical[source_column]
            _dim01_add_cell!(assembly,msh,(dim,owners[source_column]),mapping[connectivity],
                !isempty(cell_births) ? cell_births[birth_offset+source_column] :
                allocate_cells ? _dim01_next_tag!(assembly,false) : UInt64(0);
                source_tag=nothing,physical=legacy_physical)
            traces===nothing || (traces[bi][source_column]=last(assembly.cells))
        end
        birth_offset+=size(cells,2)
    end
    push!(emitted,(dimension,Int32(tag)))
    geometry_cells===nothing || (geometry_cells[mesh]=traces)
    return mapping
end

function _higher_analytic_boundary_surfaces(m,dimension,caller)
    skipped=Set{Int}()
    dimension==3 || return skipped
    for volume in keys(m.volumes)
        # A Layers-extruded primitive sweeps through its boundary products —
        # the surfaces classify like an explicit shell's and must mesh.
        (Model._implicit_volume_surface(m,volume) &&
         Model._extrude_gate(m,3,volume)===nothing) || continue
        for signed in Model._model_volume_boundary_surfaces(m,volume,caller)
            surface=abs(signed)
            push!(skipped,surface)
        end
    end
    return skipped
end

# Ordinary Layers constructors weld their emitted positions. Classification
# must bind that emitted position registry to the current lower cell UIDs;
# fitting a curved lateral to a plane cannot supply those ownership facts.
# QuadTri retains its separate planner. Late subdivision supplies certified
# operation-owned lateral incidence, followed by a final actual-cell check.
function _higher_ordinary_sweep_classification(assembly,m,production,input,dimension,tag,caller;
                                               primary_existing=false,planned_laterals=nothing)
    mesh=production.mesh
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: ordinary Layers producer is invalid — "*join(diagnostic.messages,"; ")))
    params=Model._extrude_gate(m,dimension,tag)
    params===nothing && return nothing
    params.quad_to_tri===:none || return nothing
    input!==nothing && input.mesh!==nothing || throw(ArgumentError(
        "$caller: ordinary Layers entity ($dimension,$tag) lacks actual lower products"))
    primary=(dimension,Int32(tag))
    closure=Model._model_mesh_projection_closure(m,dimension,tag,caller)
    actual=input.mesh;data=actual.entity_data
    emitted=production.position_rows
    length(emitted)==nnodes(mesh) && all(row->1<=row<=nnodes(mesh),values(emitted)) &&
        allunique(values(emitted)) || throw(ArgumentError(
            "$caller: ordinary Layers position receipt is not bijective"))
    for (point,row) in emitted
        isequal(point,ntuple(axis->mesh.coords[axis,row],3)) || throw(ArgumentError(
            "$caller: ordinary Layers position receipt changed after emission"))
    end
    # Base source rows are indexed by the renderer, independently of the
    # global position registry. COPIED_ENTITY tops carry transformed columns
    # rather than retaining their source UID on the cap.
    link=m.meshing.extrude_sources[(dimension,tag)]
    source_product=Model._model_mesh_product(input,link[1],abs(link[2]);caller)
    source_product!==nothing &&
        size(production.column_rows,2)==length(production.source_node_tags) &&
        allunique(production.source_node_tags) &&
        Set(production.source_node_tags)==Set(source_product.entity_data.external_node_tags) ||
        throw(ArgumentError("$caller: ordinary Layers entity ($dimension,$tag), Source $link, "*
            "columns $(size(production.column_rows)), emitted UIDs $(production.source_node_tags), "*
            "actual UIDs $(source_product===nothing ? nothing : source_product.entity_data.external_node_tags) "*
            "lack complete actual UID provenance"))
    us=Model._extrude_level_us(params)
    copied=dimension==2 && link[1]==2
    levels=copied ? us[end:end] : us
    size(production.column_rows,1)==length(levels) || throw(ArgumentError(
        "$caller: ordinary Layers source columns have inconsistent layer provenance"))
    spec=m.meshing.extrude_specs[(dimension,tag)]
    transforms=Model._extrude_level_transforms(spec,levels,caller)
    source_positions=Dict(uid=>row for (row,uid) in enumerate(source_product.entity_data.external_node_tags))
    used_source=Set(source_product.entity_data.external_node_tags[row]
        for block in source_product.blocks for row in block.nodes)
    for (column,uid) in enumerate(production.source_node_tags)
        source=get(assembly.node_positions,uid,Int32(0))
        source>0 || throw(ArgumentError("$caller: ordinary Layers column loses source UID $uid"))
        source_row=source_positions[uid];node=assembly.nodes[source]
        node.owner==source_product.entity_data.node_entities[source_row] &&
            isequal(node.point,ntuple(axis->source_product.coords[axis,source_row],3)) ||
            throw(ArgumentError("$caller: ordinary Layers column changed actual Source UID $uid"))
        orphan=(dimension==3 || (dimension==2 && link[1]==1)) && !(uid in used_source)
        orphan && node.owner!=(link[1],Int32(abs(link[2]))) && throw(ArgumentError(
            "$caller: ordinary Layers has foreign unused source UID $uid"))
        for (level,u) in enumerate(levels)
            row=production.column_rows[level,column]
            if orphan && (primary_existing || level in (1,length(levels)))
                row==0 || throw(ArgumentError(
                    "$caller: ordinary Layers orphan source UID $uid fabricated a boundary event"))
                continue
            end
            1<=row<=nnodes(mesh) || throw(ArgumentError(
                "$caller: ordinary Layers source column has an invalid emitted row"))
            point=!copied && level==1 ? node.point :
                Model._extrude_level_at(spec,node.point,u,transforms,level)
            isequal(point,ntuple(axis->mesh.coords[axis,row],3)) || throw(ArgumentError(
                "$caller: ordinary Layers source column does not consume actual UID $uid"))
        end
    end
    node_tags=zeros(UInt64,nnodes(mesh));retained=falses(nnodes(mesh))
    node_owners=fill(primary,nnodes(mesh))
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    source_rows=Dict{UInt64,Int32}()
    function bind_source(row)
        uid=data.external_node_tags[row]
        previous=get(source_rows,uid,Int32(0));previous>0 && return previous
        point=ntuple(axis->actual.coords[axis,row],3)
        target=get(emitted,point,Int32(0))
        target>0 || throw(ArgumentError(
            "$caller: ordinary Layers entity ($dimension,$tag) omits lower cell UID $uid"))
        if retained[target]
            node_tags[target]==uid || throw(ArgumentError(
                "$caller: ordinary Layers emitter welds distinct lower cell UIDs"))
        end
        position=get(assembly.node_positions,uid,Int32(0))
        position>0 || throw(ArgumentError(
            "$caller: ordinary Layers lower cell UID $uid lacks its native record"))
        source=assembly.nodes[position]
        source.owner==data.node_entities[row] && isequal(source.point,point) || throw(ArgumentError(
            "$caller: ordinary Layers lower cell UID $uid changed its actual ownership or coordinates"))
        source.owner in closure || throw(ArgumentError(
            "$caller: ordinary Layers lower cell UID $uid belongs outside its entity closure"))
        node_tags[target]=uid;retained[target]=true;node_owners[target]=source.owner
        parameters[target]=source.parameter===nothing ? nothing : copy(source.parameter)
        source_rows[uid]=target
        return target
    end
    if primary_existing
        dimension==2 || throw(ArgumentError(
            "$caller: existing primary rows are reserved for a lateral cell remesh"))
        primary_product=get(input.products,primary,nothing)
        primary_product!==nothing && primary_product.owner==primary || throw(ArgumentError(
            "$caller: ordinary lateral remesh lacks its actual Surface product"))
        for row in primary_product.rows
            bind_source(row)
        end
        all(retained) || throw(ArgumentError(
            "$caller: ordinary lateral cell remesh invented a new vertex event"))
    end
    edges=Dict{NTuple{2,Int32},Tuple{Int,Int32}}()
    faces=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    quads=Dict{NTuple{4,Int32},Tuple{Int,Int32}}()
    edge_conflicts=Set{NTuple{2,Int32}}();face_conflicts=Set{NTuple{3,Int32}}();quad_conflicts=Set{NTuple{4,Int32}}()
    owners=Dict{Int,Vector{Int32}}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        isempty(cells) && continue
        msh_dimension(msh)==dimension || throw(ArgumentError(
            "$caller: ordinary Layers producer mixes cell dimensions"))
        cell_owners=fill(Int32(tag),size(cells,2))
        append!(get!(()->Int32[],owners,msh),cell_owners)
        _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,msh,cells,cell_owners)
    end
    expected_edges=Set{NTuple{2,Int32}}()
    lower_lines=Set{Tuple{Tuple{Int,Int32},NTuple{2,Int32}}}()
    expected_faces=Set{Tuple}()
    embedded_edges=Set{NTuple{2,Int32}}();embedded_faces=Set{Tuple}()
    # Closed-strip wires list their seam curve with both orientations;
    # its cells are interior support, not emitted boundary.
    boundary_curves=dimension==2 ?
        Set(Model._model_emitted_boundary_curves(m,tag)) : Set{Int}()
    boundary_surfaces=dimension==3 ? Set(abs.(Model._model_volume_boundary_surfaces(m,tag,caller))) : Set{Int}()
    lower=sort!(Tuple{Int,Int32}[owner for owner in closure if 0<owner[1]<dimension])
    for owner in lower
        lower_product=get(input.products,owner,nothing)
        lower_product===nothing && throw(ArgumentError(
            "$caller: ordinary Layers entity ($dimension,$tag) lacks lower product $owner"))
        lower_product.owner==owner && length(lower_product.blocks)==length(lower_product.element_tags) || throw(ArgumentError(
            "$caller: ordinary Layers lower product has inconsistent identity provenance"))
        any(block->!isempty(block.nodes),lower_product.blocks) || throw(ArgumentError(
            "$caller: ordinary Layers entity ($dimension,$tag) lacks actual lower cells for $owner"))
        planned=owner[1]==2 && planned_laterals!==nothing ?
            get(planned_laterals,Int(owner[2]),nothing) : nothing
        if planned!==nothing
            planned.input===input && planned.owner==owner || throw(ArgumentError(
                "$caller: ordinary lateral cell plan has foreign ownership or input"))
            checked=_higher_ordinary_sweep_classification(assembly,m,planned.production,input,
                2,Int(owner[2]),caller;primary_existing=true)
            checked!==nothing && checked.node_tags==planned.node_tags || throw(ArgumentError(
                "$caller: ordinary lateral cell plan changed its actual UID incidence"))
        end
        blocks=planned===nothing ? lower_product.blocks : ElementBlock[
            ElementBlock(msh,cells) for (msh,cells,_) in _cache_native_blocks(planned.production.mesh)
            if msh_dimension(msh)==2 && !isempty(cells)]
        planned_rows=planned===nothing ? nothing :
            Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
        for block in blocks
            msh_dimension(block.msh)==owner[1] || throw(ArgumentError(
                "$caller: ordinary Layers lower product has inconsistent cell ownership"))
            mapped=Matrix{Int32}(undef,size(block.nodes))
            for column in axes(block.nodes,2),slot in axes(block.nodes,1)
                local_row=block.nodes[slot,column]
                actual_row=planned===nothing ? lower_product.rows[local_row] :
                    get(planned_rows,planned.node_tags[local_row],Int32(0))
                actual_row>0 || throw(ArgumentError(
                    "$caller: ordinary lateral cell plan lost an actual node UID"))
                mapped[slot,column]=bind_source(actual_row)
            end
            if owner[1]==1
                block.msh==1 || throw(ArgumentError(
                    "$caller: ordinary Layers requires actual linear Curve cells"))
                # A two-arc closed ring — a vertex generatrix swept
                # through two levels — legitimately lists its only node
                # pair once per arc; any other repeat is an overlap.
                ring2=m.curves[Int(owner[2])][1]==m.curves[Int(owner[2])][2] &&
                      sum(b->size(b.nodes,2),blocks)==2
                for column in axes(mapped,2)
                    edge=minmax(mapped[1,column],mapped[2,column])
                    edge[1]!=edge[2] && (ring2 || !((owner,edge) in lower_lines)) || throw(ArgumentError(
                        "$caller: ordinary Layers lower Curve cells repeat a Line incidence"))
                    push!(lower_lines,(owner,edge))
                    push!(dimension==2 && !(Int(owner[2]) in boundary_curves) ? embedded_edges : expected_edges,edge)
                end
            else
                block.msh in (2,3) || throw(ArgumentError(
                    "$caller: ordinary Layers requires actual linear Surface cells"))
                for column in axes(mapped,2)
                    face=Tuple(sort!(collect(@view mapped[:,column])))
                    target=Int(owner[2]) in boundary_surfaces ? expected_faces : embedded_faces
                    face in target && throw(ArgumentError(
                        "$caller: ordinary Layers lower Surface cells overlap"))
                    push!(target,face)
                end
            end
            _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
                Int(block.msh),mapped,fill(owner[2],size(mapped,2)))
        end
    end
    embedded_points=Int[]
    if dimension==2
        points,_=Model._model_surface_embedding_tags(m,tag,caller)
        append!(embedded_points,points)
    else
        _,_,surfaces,_,_,_=Model._model_volume_embedding_inventory(m,tag,caller)
        for surface in surfaces
            points,_=Model._model_surface_embedding_tags(m,surface,caller)
            append!(embedded_points,points)
        end
        append!(embedded_points,Int[point for (dim,point) in
            get(m.embeds,(3,tag),Tuple{Int,Int}[]) if dim==0])
    end
    used=Set(row for (_,cells,_) in _cache_native_blocks(mesh) for row in cells)
    for point in unique(embedded_points)
        owner=(0,Int32(point));point_product=get(input.products,owner,nothing)
        point_product!==nothing && all(block->block.msh==15,point_product.blocks) || throw(ArgumentError(
            "$caller: ordinary Layers embedded Point[$point] lacks actual Point cells"))
        rows=unique(Int32[point_product.rows[row] for block in point_product.blocks for row in block.nodes])
        length(rows)==1 || throw(ArgumentError(
            "$caller: ordinary Layers embedded Point[$point] has ambiguous actual incidence"))
        bind_source(only(rows)) in used || throw(ArgumentError(
            "$caller: ordinary Layers embedded Point[$point] is absent from emitted cells"))
    end
    _api_check_support_conflicts(edge_conflicts,face_conflicts,quad_conflicts)
    if dimension==2
        boundary,actual_edges=Model._surface_boundary_topology(mesh,caller)
        expected_edges==actual_edges || throw(ArgumentError(
            "$caller: ordinary Layers lower Curve cells do not exactly cover its emitted boundary"))
        actual_support=Set{NTuple{2,Int32}}()
        for (msh,cells,_) in _cache_native_blocks(mesh),column in axes(cells,2),slots in
            Model._surface_cell_edge_slots(size(cells,1))
            push!(actual_support,minmax(cells[slots[1],column],cells[slots[2],column]))
        end
        isempty(intersect(embedded_edges,actual_edges)) && issubset(embedded_edges,actual_support) ||
            throw(ArgumentError("$caller: ordinary Layers embedded Curves do not lie on interior support"))
        all(row->!boundary[row] || retained[row],eachindex(boundary)) || throw(ArgumentError(
            "$caller: ordinary Layers boundary vertex lacks an actual source UID"))
    else
        Model._volume_projection_contract(mesh,caller)
        actual_faces,_=Model._volume_cell_boundary_faces(mesh)
        expected_faces==actual_faces || throw(ArgumentError(
            "$caller: ordinary Layers lower Surface cells do not exactly cover its emitted Volume boundary"))
        isempty(intersect(embedded_faces,actual_faces)) || throw(ArgumentError(
            "$caller: ordinary Layers embedded Surface cannot alias the outer shell"))
        all_faces=Set{Tuple}()
        for (msh,cells,_) in _cache_native_blocks(mesh),column in axes(cells,2),slots in Model._VOLUME_CELL_FACES[msh]
            push!(all_faces,Tuple(sort!(Int32[cells[slot,column] for slot in slots])))
        end
        issubset(embedded_faces,all_faces) || throw(ArgumentError(
            "$caller: ordinary Layers omits an embedded Surface cell from its interior support"))
        actual_edges=Model._volume_edge_set(mesh)
        issubset(expected_edges,actual_edges) || throw(ArgumentError(
            "$caller: ordinary Layers omits an actual Curve cell from its Volume edges"))
    end
    boundaries=Dict(owner=>Int32[abs(child) for child in data.entities[owner].boundaries] for owner in closure)
    class=_mixed_classification(mesh,primary,sort!(collect(closure)),node_owners,boundaries,owners;
        edge_entities=edges,face_entities=faces,quad_entities=quads,node_parametric=parameters)
    return (class=class,node_tags=node_tags,retained=retained)
end

function _higher_surface_production_classification(m,production,tag)
    mesh=production.mesh
    owner=(2,Int32(tag))
    entities=Tuple{Int,Int32}[owner]
    queue=Tuple{Int,Int32}[owner]
    boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}()
    while !isempty(queue)
        current=pop!(queue)
        children=model_boundary(m,[(current[1],Int(current[2]))],false,false,false)
        boundaries[current]=Int32[abs(child[2]) for child in children]
        for child in children
            key=(child[1],Int32(abs(child[2])))
            key in entities && continue
            push!(entities,key);push!(queue,key)
        end
    end
    sort!(entities)
    length(production.node_owners)==nnodes(mesh) || throw(ArgumentError(
        "native Surface producer owner records do not align"))
    cells=Dict{Int,Vector{Int32}}()
    for (msh,nodes,_) in _cache_native_blocks(mesh)
        cells[msh]=fill(Int32(tag),size(nodes,2))
    end
    if mesh isa MixedMesh
        return _mixed_classification(mesh,owner,entities,production.node_owners,
                                     boundaries,cells)
    end
    return _MeshClassification(mesh,owner,entities,production.node_owners,boundaries,
        get(cells,1,Int32[]),get(cells,2,Int32[]),get(cells,4,Int32[]),cells,nothing,
        Dict{NTuple{2,Int32},Tuple{Int,Int32}}(),
        Dict{NTuple{3,Int32},Tuple{Int,Int32}}(),
        Dict{NTuple{4,Int32},Tuple{Int,Int32}}(),production.node_parameters)
end

function _higher_geometry_source_points(assembly,birth)
    birth===nothing && return nothing
    points=Dict{UInt64,NTuple{3,Float64}}()
    for (tag,node) in birth.nodes
        row=get(assembly.node_positions,tag,Int32(0))
        row>0 && assembly.nodes[row].owner==node.owner || continue
        points[tag]=node.point
    end
    return points
end

function _higher_geometry_source_point(assembly,birth,uid)
    birth===nothing && return nothing
    node=get(birth.nodes,uid,nothing)
    node===nothing && return nothing
    row=get(assembly.node_positions,uid,Int32(0))
    row>0 && assembly.nodes[row].owner==node.owner || return nothing
    return node.point
end

function _higher_geometry_emitter_part(mesh,mapping,assembly,birth,geometry_cells)
    birth===nothing && return mesh
    result=mesh
    for (row,position) in enumerate(mapping)
        source=assembly.nodes[position]
        alternate=get(birth.nodes,source.tag,nothing)
        alternate===nothing && continue
        alternate.owner==source.owner || continue
        isequal((result.coords[1,row],result.coords[2,row],result.coords[3,row]),
                alternate.point) && continue
        result===mesh && (result=_copy_mesh(mesh))
        result.coords[:,row].=alternate.point
    end
    result===mesh || (geometry_cells[result]=geometry_cells[mesh])
    return result
end

# Native Mesh2D repeats the tag-ordered face scan: a periodic slave stays
# pending until its master completes. Keep that allocation order in the
# detached operation instead of constructing an untracked recursive master.
function _higher_surface_units(m,units,caller)
    pending=Set(tag for tags in units for tag in tags)
    remaining=copy(units);ordered=Vector{Vector{Int}}()
    while !isempty(remaining)
        deferred=Vector{Vector{Int}}()
        for tags in remaining
            constraint=length(tags)==1 ? get(m.periodic,(2,only(tags)),nothing) : nothing
            if constraint!==nothing && Int(constraint.master_entity) in pending
                push!(deferred,tags)
            else
                push!(ordered,tags);setdiff!(pending,tags)
            end
        end
        length(deferred)<length(remaining) || throw(ArgumentError(
            "$caller: cyclic periodic Surface dependency among " *
            join(sort!(collect(pending)),", ")))
        remaining=deferred
    end
    return ordered
end

# TF Face cells receive their UIDs in the API assembly, after the geometric
# producer returns. Bind the just-created event vectors, including native
# cyclic orientation, rather than expecting preallocated producer labels.
function _higher_transfinite_constructor_input(actual,events,assembly,tag)
    isempty(events) && return nothing
    caller=assembly.caller;owner=(2,Int32(tag));data=actual.entity_data
    data!==nothing && all(cell->cell.owner==owner && cell.msh in (2,3),events) &&
        allunique(cell.tag for cell in events) ||
        throw(ArgumentError("$caller: transfinite constructor lacks its actual Face cell events"))
    for family in (2,3)
        emitted=[(cell.tag,cell.physical,
                  Tuple(assembly.nodes[row].tag for row in cell.nodes))
                 for cell in events if cell.msh==family]
        bound=[(uid,block.tags[column],
                Tuple(data.external_node_tags[block.nodes[:,column]]))
               for (block,uids,owners) in zip(actual.blocks,data.external_element_tags,data.block_entities)
               if block.msh==family
               for (column,uid) in enumerate(uids) if owners[column]==owner[2]]
        emitted==bound || throw(ArgumentError(
            "$caller: transfinite constructor changes its actual Face emission vector"))
    end
    return actual
end

function _higher_surface_phase!(assembly,m,options,entities,discarded,done;
                                size_field=nothing,metric_surface=nothing,
                                source_nodes=length(assembly.nodes),ranks=nothing,
                                old_supports=nothing,geometry_parts=nothing,dimension=2,
                                 geometry_cells=nothing,geometry_birth=nothing,
                                 geometry_alternates=nothing,curve_births=Dict{Tuple{Int,Int32},Vector{UInt64}}(),
                                 periodic_pairs=nothing,native_orphan_births=nothing,native_recombine_events=nothing,
                                 surface_constructor_inputs=nothing)
    # Mesh2D resets all face statuses before OnlyEmpty/OnlyVisible dispatch.
    filter!(owner->owner[1]!=2,done)
    have=Set(cell.owner for cell in assembly.cells)
    retained=Set{Tuple{Int,Int32}}()
    meshed=Int[]
    analytic_boundaries=_higher_analytic_boundary_surfaces(m,dimension,assembly.caller)
    for (dim,tag) in entities
        dim==2 || continue
        owner=(2,Int32(tag))
        if tag in analytic_boundaries
            # These native Volume producers retain their established analytic
            # boundary encoding. Model does not yet admit a curved face as an
            # independent Surface product; do not turn a supported Volume
            # request into that standalone-face diagnostic.
            owner in have && throw(ArgumentError(
                "$(assembly.caller): actual curved Surface[$tag] cells cannot be consumed " *
                "by the analytic primitive Volume producer"))
            continue
        end
        keep=_dim01_fully_discrete(m,(2,tag)) ||
             (options.mesh_only_visible && model_entity_visibility(m,2,tag)==0) ||
             (options.mesh_only_empty && owner in have)
        if keep
            push!(retained,owner)
        else
            push!(meshed,tag);push!(discarded,owner)
        end
    end
    filter!(cell->cell.source_tag===nothing || !(cell.owner in discarded),assembly.cells)
    # The actual products are supplied to extrusion source/cap consumers. A
    # retained surface's cell order is its stored order, including cap pairs.
    ranks===nothing || _higher_complete_ranks!(ranks,assembly)
    input=_higher_assembly_authority(assembly,m,done;ranks,old_supports)
    input=_higher_bind_producer_rows(input,assembly,Tuple[],m.meshing.mesh_authority;curve_births)
    input=Model._model_mesh_with_curve_births(input,curve_births;caller=assembly.caller)
    emitted=Set{Tuple{Int,Int32}}()
    scope=isempty(meshed) ? nothing : Model._extrude_nonew_scope(
        m,assembly.caller;size_field=size_field,_working_model=true,_mesh_input=input)
    periodic_masters=Set(Int(constraint.master_entity) for constraint in
        Model.model_periodic_constraints(m) if constraint.dim==2)
    # A copied top and a later Region consume the just-emitted actual Source,
    # including its edited coordinates and its own constructor row identities.
    ordinary_sources=Set{Int}()
    unmerged_nonew_sources=Set{Int}()
    for volume in keys(m.volumes)
        Model._extrude_nonew_unmerged_b4_incidence(m,volume,assembly.caller) || continue
        push!(unmerged_nonew_sources,abs(m.meshing.extrude_sources[(3,volume)][2]))
    end
    for ((dim,entity),link) in m.meshing.extrude_sources
        dim in (2,3) && link[1]==2 || continue
        params=Model._extrude_gate(m,dim,entity)
        params!==nothing && params.quad_to_tri===:none || continue
        push!(ordinary_sources,abs(link[2]))
    end
    constructor_inputs=surface_constructor_inputs===nothing ?
        Dict{Tuple{Int,Int32},MixedMesh}() : surface_constructor_inputs
    # These Faces have already lost their old cells from the detached assembly.
    # Kept Faces retain their owned constructor vectors; fresh Faces publish
    # new vectors as their actual cell events arrive in this phase.
    for tag in meshed
        delete!(constructor_inputs,(2,Int32(tag)))
    end
    units=_higher_surface_units(m,_generate_units(m,2,sort!(meshed)),assembly.caller)
    for tags in units
        products=metric_surface!==nothing ? [(only(tags),mesh_model_surface(m,only(tags),size_field;
            max_metric_edge=metric_surface.max_metric_edge,max_nodes=metric_surface.max_nodes,
            max_triangles=metric_surface.max_triangles,_mesh_input=input,_mesh1d_options=options,
            _mesh_node_high_water=assembly.node_max,_mesh_element_high_water=assembly.cell_max,
            _return_surface_production=true))] : length(tags)==1 ? [(tags[1],mesh_model_surface(m,tags[1];
            size_field=size_field,_extrude_scope=scope,_mesh_input=input,_mesh1d_options=options,
            _mesh_node_high_water=assembly.node_max,_return_surface_production=true,
            _mesh_element_high_water=assembly.cell_max,
            _geometry_points=_higher_geometry_source_points(assembly,geometry_birth)))] :
            _compound_surface_meshes(m,tags,assembly.caller;size_field=size_field,
                _mesh1d_options=options)
        for (tag,produced) in products
            alternate=produced isa Model._ModelSurfaceGeometryProduction ? produced.geometry : nothing
            source_binding=produced isa Model._ModelSurfaceGeometryProduction ? produced.production :
                produced isa Model._ModelSurfaceProduction ? produced : nothing
            if source_binding!==nothing && native_recombine_events!==nothing
                append!(native_recombine_events,source_binding.native_recombine_events)
            end
            ordinary=produced isa Model._ExtrudeMeshProduction ? produced : nothing
            mesh=source_binding!==nothing ? source_binding.mesh : ordinary!==nothing ? ordinary.mesh : produced
            swept=ordinary===nothing ? nothing : _higher_ordinary_sweep_classification(
                assembly,m,ordinary,input,2,tag,assembly.caller)
            swept===nothing && (swept=_higher_nonew_identity_classification(
                assembly,m,mesh,input,scope,2,tag,assembly.caller))
            class=source_binding!==nothing ? _higher_surface_production_classification(m,source_binding,tag) :
                swept!==nothing ? swept.class : _classify_cached_mesh(
                    m,_classification_skeleton(mesh,false),2,tag,mesh;_extrude_scope=scope)
            first_cell_event=length(assembly.cells)+1
            mapping=_higher_append_part!(assembly,m,mesh,class,2,tag,discarded;
                emitted,source_nodes,geometry_cells,source_binding,
                source_row_tags=swept===nothing ? nothing : swept.node_tags,
                source_row_retained=swept===nothing ? nothing : swept.retained,
                source_rows_existing_only=swept!==nothing,
                native_orphan_births=ordinary===nothing ? nothing : native_orphan_births)
            constraint=get(m.periodic,(2,tag),nothing)
            if periodic_pairs!==nothing && constraint!==nothing && source_binding!==nothing
                master=Int(constraint.master_entity)
                if haskey(m.meshing.transfinite_surfaces,tag)
                    # Explicit slave TF meshing completes independently of the
                    # periodic copy path and replaces the old Face map by empty.
                    periodic_pairs[tag]=(master,Tuple{UInt64,UInt64}[])
                else
                    source=Model._model_periodic_surface_source(m,constraint,input,assembly.caller)
                    source!==nothing && length(source.entity_data.external_node_tags)==length(source_binding.node_tags) ||
                        throw(ArgumentError("$(assembly.caller): periodic Surface copy loses its source UID rows"))
                    periodic_pairs[tag]=(master,collect(zip(source_binding.node_tags,source.entity_data.external_node_tags)))
                end
            end
            if alternate!==nothing && geometry_alternates!==nothing
                geometry_alternates[mesh]=alternate
                geometry_cells[alternate]=geometry_cells[mesh]
                for (row,position) in enumerate(mapping)
                    node=assembly.nodes[position]
                    point=(alternate.coords[1,row],alternate.coords[2,row],alternate.coords[3,row])
                    geometry_birth.nodes[node.tag]=_Dim01Node(node.tag,point,node.owner,node.parameter)
                end
            end
            geometry_parts===nothing || push!(geometry_parts,(tag,mesh,class,mapping))
            push!(done,(2,Int32(tag)))
            # Publish this operation-owned product to downstream periodic
            # copies. The old authority remains immutable and lower source
            # identities retain their private ranks and actual coordinates.
            native_surface=source_binding!==nothing
            if tag in periodic_masters || tag in ordinary_sources || tag in unmerged_nonew_sources || native_surface ||
                    (source_binding!==nothing && !isempty(source_binding.source_node_updates))
                ranks===nothing || _higher_complete_ranks!(ranks,assembly)
                previous_input=input
                input=_higher_assembly_authority(assembly,m,done;ranks,old_supports)
                input=_higher_bind_producer_rows(input,assembly,[(tag,mesh,class,mapping)],
                    previous_input)
                if source_binding!==nothing
                    transfinite_events=source_binding.constructor_mesh===nothing &&
                        isempty(source_binding.element_tags) &&
                        haskey(m.meshing.transfinite_surfaces,tag)
                    actual_product=Model._model_mesh_product(input,2,tag;
                        caller=assembly.caller,
                        producer_order=source_binding.constructor_mesh!==nothing || transfinite_events ||
                            metric_surface!==nothing)
                    if transfinite_events
                        constructor=_higher_transfinite_constructor_input(actual_product,
                            assembly.cells[first_cell_event:end],assembly,tag)
                        constructor===nothing || (constructor_inputs[(2,Int32(tag))]=constructor)
                    elseif source_binding.constructor_mesh===nothing
                        isempty(source_binding.element_tags) ||
                            (constructor_inputs[(2,Int32(tag))]=Model._model_surface_constructor_input(
                                source_binding,actual_product,assembly.caller;constructor=mesh))
                    else
                        constructor_inputs[(2,Int32(tag))]=Model._model_surface_constructor_input(
                            source_binding,actual_product,assembly.caller)
                    end
                end
                input=Model._model_mesh_with_surface_constructor_inputs(input,constructor_inputs;
                    caller=assembly.caller)
                if tag in unmerged_nonew_sources
                    scope=Model._extrude_nonew_scope(m,assembly.caller;size_field=size_field,
                        _working_model=true,_mesh_input=input)
                end
            end
        end
    end
    return input,scope,retained
end

function _higher_mesh_periodic_links(mesh,links)
    return MixedMesh(mesh.coords,mesh.blocks;physical_names=mesh.physical_names,
        entity_data=mesh.entity_data,elementary_entities=mesh.elementary_entities,
        periodic_links=links,ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
end

function _higher_periodic_source_mesh(mesh,pairs)
    mesh isa MixedMesh && !isempty(pairs) || return mesh
    links=Elements.MixedPeriodicLink[link for link in mesh.periodic_links
        if !(link.dim==2 && haskey(pairs,Int(link.slave_entity)))]
    length(links)==length(mesh.periodic_links) && return mesh
    # Regenerated slaves replace their old correspondence. Retained lower and
    # other-dimensional links still pass through the normal UID remapping.
    return _higher_mesh_periodic_links(mesh,links)
end

function _higher_periodic_surface_links(plan,pairs,m,caller)
    isempty(pairs) && return plan
    mesh=plan.mesh
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(mesh.entity_data.external_node_tags))
    links=Elements.MixedPeriodicLink[link for link in mesh.periodic_links
        if !(link.dim==2 && haskey(pairs,Int(link.slave_entity)))]
    for (slave,(master,bindings)) in sort!(collect(pairs);by=first)
        slaves=Int32[];masters=Int32[]
        for (source_slave,source_master) in bindings
            slave_uid=get(plan.node_tag_map,source_slave,nothing)
            master_uid=get(plan.node_tag_map,source_master,nothing)
            slave_row=get(rows,slave_uid,Int32(0));master_row=get(rows,master_uid,Int32(0))
            slave_row>0 && master_row>0 || throw(ArgumentError(
                "$caller: periodic Surface correspondence loses an actual vertex"))
            push!(slaves,slave_row);push!(masters,master_row)
        end
        base=Elements.MixedPeriodicLink(2,slave,master,slaves,masters;
            affine=m.periodic[(2,slave)].affine)
        lifted=isempty(bindings) && haskey(m.meshing.transfinite_surfaces,slave) ?
            Pair{Int32,Int32}[] : _periodic_surface_lift(mesh,plan.class,base,
                collect(zip(slaves,masters)),caller)
        push!(links,Elements.MixedPeriodicLink(2,slave,master,
            Int32[first(pair) for pair in lifted],Int32[last(pair) for pair in lifted];
            affine=base.affine))
    end
    published=_higher_mesh_periodic_links(mesh,links)
    return merge(plan,(mesh=published,class=_mixed_rebind_class(plan.class,published)))
end

function _higher_retained_curves!(assembly,m,entities,first_nodes,high_order,done)
    have=Set(cell.owner for cell in assembly.cells)
    owned=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (position,node) in enumerate(assembly.nodes)
        push!(get!(()->Int32[],owned,node.owner),Int32(position))
    end
    for (dim,tag) in entities
        owner=(dim,Int32(tag))
        dim==1 && owner in done && !(owner in have) || continue
        params=get(m.curve_params,tag,nothing)
        params===nothing && continue
        # A legacy higher-dimensional cache omitted lower Line cells while
        # retaining their actual primary vertices and stored discretization.
        # Complete this skipped phase without renewing those primary vertices.
        _dim01_emit_curve!(assembly,m,tag,params,first_nodes;
                           reuse=true,high_order,owned_nodes=owned)
        push!(have,owner)
    end
    return nothing
end

function _higher_ingest_native_inputs!(assembly,m,authority)
    authority===nothing && return nothing
    have=Set(cell.owner for cell in assembly.cells)
    owners=Set(keys(authority.native_inputs))
    for (owner,mesh) in sort!(collect(authority.native_inputs);by=first)
        owner in have && continue
        # These are the original completed boundary products, retained when a
        # legacy Volume-only cache did not expose its Surface cells. Their
        # vectors are never reconstructed from that volume's boundary faces.
        dim,tag=owner
        class=_classify_cached_mesh(m,_classification_skeleton(mesh,false),dim,Int(tag),mesh)
        _higher_append_part!(assembly,m,mesh,class,dim,Int(tag),Set{Tuple{Int,Int32}}();
                             emitted=owners,allocate_cells=false)
        push!(have,owner)
    end
    return nothing
end

@inline function _higher_triangle_vertex_key(assembly,cell)
    a=assembly.nodes[cell.nodes[1]].tag
    b=assembly.nodes[cell.nodes[2]].tag
    c=assembly.nodes[cell.nodes[3]].tag
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end

function _higher_native_region_surface_order!(assembly,m,tag)
    m.meshing.algorithm_3d in (1,3) || return nothing
    Model._model_factory_volume_encoding(m,tag) && return nothing
    # meshGRegion sorts the real Region's boundary and embedded Surface
    # triangles after boundary recovery using the sorted vertex identifiers.
    # Reorder only detached actual cells; producer rows and the legacy geometry
    # carrier retain their independent insertion order.
    owners=Set((2,Int32(abs(surface))) for surface in
        Model._model_volume_boundary_surfaces(m,tag,assembly.caller))
    for (dim,surface) in get(m.embeds,(3,tag),Tuple{Int,Int}[])
        dim==2 && push!(owners,(2,Int32(surface)))
    end
    positions=Dict{Tuple{Int,Int32},Vector{Int}}()
    for (index,cell) in enumerate(assembly.cells)
        cell.owner in owners && msh_spec(cell.msh).family===:tri || continue
        push!(get!(positions,cell.owner,Int[]),index)
    end
    for rows in values(positions)
        cells=assembly.cells[rows]
        sort!(cells;by=cell->_higher_triangle_vertex_key(assembly,cell))
        assembly.cells[rows]=cells
    end
    return nothing
end

# Native volume rows carry real born UIDs. Bind their lower ownership and
# support incidence from the registered actual products before CAD projection;
# edited source vertices need not lie on the original CAD Point/Curve image.
# A transfinite fill has explicit boundary UIDs and unallocated interior events.
# Its classification consumes those row bindings and the actual lower cells.
# CAD proximity cannot replace either incidence certificate on a warped face.
function _higher_transfinite_volume_classification(assembly,m,mesh,input,tag,
                                                  source_tags,retained,caller;grid_receipt=nothing)
    input!==nothing && input.mesh.entity_data!==nothing || throw(ArgumentError(
        "$caller: transfinite Volume[$tag] lacks actual lower products"))
    length(source_tags)==length(retained)==nnodes(mesh) &&
        all(row->retained[row] ? source_tags[row]>0 : iszero(source_tags[row]),eachindex(source_tags)) &&
        allunique(source_tags[retained]) || throw(ArgumentError(
            "$caller: transfinite Volume[$tag] has inconsistent boundary row identities"))
    certified=_classification_skeleton(mesh,get(m.meshing.reverse,(3,tag),false))
    diagnostic=validate(certified)
    diagnostic.ok || throw(ArgumentError(
        "$caller: transfinite Volume[$tag] has invalid cell geometry — "*join(diagnostic.messages,"; ")))
    Model._volume_projection_contract(certified,caller)
    actual=input.mesh;data=actual.entity_data
    primary=(3,Int32(tag))
    closure=Model._model_mesh_projection_closure(m,3,tag,caller)
    targets=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_tags) if retained[row])
    source_rows=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    node_owners=fill(primary,nnodes(mesh))
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    for (uid,row) in targets
        source=get(source_rows,uid,0)
        position=get(assembly.node_positions,uid,Int32(0))
        source>0 && position>0 || throw(ArgumentError(
            "$caller: transfinite boundary UID $uid lacks its actual lower record"))
        record=assembly.nodes[position]
        owner=data.node_entities[source]
        owner in closure && owner[1]<3 && record.owner==owner &&
            isequal(record.point,ntuple(axis->actual.coords[axis,source],3)) &&
            isequal(record.point,ntuple(axis->mesh.coords[axis,row],3)) || throw(ArgumentError(
                "$caller: transfinite boundary UID $uid changed coordinates or ownership"))
        node_owners[row]=owner
        parameters[row]=record.parameter===nothing ? nothing : copy(record.parameter)
    end
    edges=Dict{NTuple{2,Int32},Tuple{Int,Int32}}()
    faces=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    quads=Dict{NTuple{4,Int32},Tuple{Int,Int32}}()
    edge_conflicts=Set{NTuple{2,Int32}}()
    face_conflicts=Set{NTuple{3,Int32}}()
    quad_conflicts=Set{NTuple{4,Int32}}()
    owners=Dict{Int,Vector{Int32}}()
    for (msh,cells,_) in _cache_native_blocks(mesh)
        isempty(cells) && continue
        msh_dimension(msh)==3 || throw(ArgumentError(
            "$caller: transfinite Volume[$tag] mixes cell dimensions"))
        cell_owners=fill(Int32(tag),size(cells,2))
        append!(get!(()->Int32[],owners,msh),cell_owners)
        _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
                            msh,cells,cell_owners)
    end
    expected_edges=Set{NTuple{2,Int32}}()
    expected_faces=Set{Tuple}()
    embedded_faces=Set{Tuple}()
    lower_lines=Set{Tuple{Tuple{Int,Int32},NTuple{2,Int32}}}()
    boundary_surfaces=Set(abs.(Model._model_volume_boundary_surfaces(m,tag,caller)))
    for owner in sort!(Tuple{Int,Int32}[owner for owner in closure if 0<owner[1]<3])
        product=get(input.products,owner,nothing)
        product!==nothing && product.owner==owner &&
            length(product.blocks)==length(product.element_tags) || throw(ArgumentError(
                "$caller: transfinite Volume[$tag] lacks lower product $owner"))
        any(block->!isempty(block.nodes),product.blocks) || throw(ArgumentError(
            "$caller: transfinite Volume[$tag] lacks actual lower cells for $owner"))
        for block in product.blocks
            msh_dimension(block.msh)==owner[1] || throw(ArgumentError(
                "$caller: transfinite lower product changes cell ownership"))
            mapped=Matrix{Int32}(undef,size(block.nodes))
            for column in axes(block.nodes,2),slot in axes(block.nodes,1)
                source=product.rows[block.nodes[slot,column]]
                uid=data.external_node_tags[source]
                row=get(targets,uid,Int32(0))
                row>0 || throw(ArgumentError(
                    "$caller: transfinite Volume[$tag] omits lower cell UID $uid"))
                mapped[slot,column]=row
            end
            if owner[1]==1
                block.msh==1 || throw(ArgumentError(
                    "$caller: transfinite Volume[$tag] requires actual linear Curve cells"))
                # Two-arc closed rings list their single node pair once
                # per arc; any other repeat is an overlap.
                ring2=m.curves[Int(owner[2])][1]==m.curves[Int(owner[2])][2] &&
                      sum(b->size(b.nodes,2),product.blocks)==2
                for column in axes(mapped,2)
                    edge=minmax(mapped[1,column],mapped[2,column])
                    edge[1]!=edge[2] && (ring2 || !((owner,edge) in lower_lines)) || throw(ArgumentError(
                        "$caller: transfinite lower Curve cells repeat a Line incidence"))
                    push!(lower_lines,(owner,edge));push!(expected_edges,edge)
                end
            else
                block.msh in (2,3) || throw(ArgumentError(
                    "$caller: transfinite Volume[$tag] requires actual linear Surface cells"))
                for column in axes(mapped,2)
                    face=Tuple(sort!(collect(@view mapped[:,column])))
                    destination=Int(owner[2]) in boundary_surfaces ? expected_faces : embedded_faces
                    face in destination && throw(ArgumentError(
                        "$caller: transfinite lower Surface cells overlap"))
                    push!(destination,face)
                end
            end
            _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
                                Int(block.msh),mapped,fill(owner[2],size(mapped,2)))
        end
    end
    _api_check_support_conflicts(edge_conflicts,face_conflicts,quad_conflicts)
    actual_faces,_=Model._volume_cell_boundary_faces(mesh)
    if expected_faces!=actual_faces
        grid_receipt isa Model._ModelTransfiniteGridReceipt &&
            Set(Int(face.owner) for face in grid_receipt.faces)==boundary_surfaces ||
            throw(ArgumentError("$caller: transfinite shell lacks its owned grid receipt"))
        supports=Model._model_transfinite_grid_boundary(
            grid_receipt,input,source_tags,actual_faces,caller;mesh)
        for (face,owner) in supports
            if length(face)==3
                faces[face]=owner
            else
                quads[face]=owner
            end
        end
    end
    isempty(intersect(embedded_faces,actual_faces)) || throw(ArgumentError(
        "$caller: transfinite embedded Surface aliases its outer shell"))
    all_faces=Set{Tuple}()
    for (msh,cells,_) in _cache_native_blocks(mesh),column in axes(cells,2),slots in Model._VOLUME_CELL_FACES[msh]
        push!(all_faces,Tuple(sort!(Int32[cells[slot,column] for slot in slots])))
    end
    issubset(embedded_faces,all_faces) &&
        issubset(expected_edges,Model._volume_edge_set(mesh)) || throw(ArgumentError(
            "$caller: transfinite Volume[$tag] omits actual lower interior support"))
    for (dimension,point) in get(m.embeds,(3,tag),Tuple{Int,Int}[])
        dimension==0 || continue
        product=get(input.products,(0,Int32(point)),nothing)
        product!==nothing && all(block->block.msh==15,product.blocks) || throw(ArgumentError(
            "$caller: transfinite embedded Point[$point] lacks actual Point cells"))
        point_rows=unique(Int32[product.rows[row] for block in product.blocks for row in block.nodes])
        length(point_rows)==1 && haskey(targets,data.external_node_tags[only(point_rows)]) ||
            throw(ArgumentError("$caller: transfinite Volume[$tag] omits its embedded Point[$point]"))
    end
    boundaries=Dict(owner=>Int32[abs(child) for child in data.entities[owner].boundaries] for owner in closure)
    return _mixed_classification(mesh,primary,sort!(collect(closure)),node_owners,boundaries,owners;
        edge_entities=edges,face_entities=faces,quad_entities=quads,node_parametric=parameters)
end


function _higher_hybrid_volume_input_classification(assembly,m,production,input,tag,caller)
    production.input===input || throw(ArgumentError(
        "$caller: hybrid Volume has a foreign actual lower constructor authority"))
    mesh=production.mesh;source_tags=production.boundary_node_tags
    constructor=production.constructor
    Model._model_native_hybrid_constructor_input(input,constructor.boundaries,
        constructor.node_uids,constructor.coordinates,caller;embedded_point=constructor.embedded_point)
    Model._model_native_hybrid_volume_certificate(mesh,source_tags,constructor,caller)
    Model._model_native_hybrid_embedded_point_certificate(mesh,source_tags,constructor,caller)
    length(source_tags)==nnodes(mesh) && allunique(source_tags) && all(!iszero,source_tags) ||
        throw(ArgumentError("$caller: hybrid Volume lacks complete current UID rows"))
    actual=input.mesh;data=actual.entity_data
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_tags))
    old_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    closure=Model._model_mesh_projection_closure(m,3,tag,caller)
    primary=(3,Int32(tag));node_entities=fill(primary,nnodes(mesh))
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    for (row,uid) in enumerate(source_tags)
        source=get(old_rows,uid,Int32(0));source==0 && continue
        owner=data.node_entities[source]
        owner in closure && all(axis->isequal(actual.coords[axis,source],mesh.coords[axis,row]),1:3) ||
            throw(ArgumentError("$caller: hybrid Volume changes or consumes a foreign lower UID $uid"))
        position=get(assembly.node_positions,uid,Int32(0))
        position>0 || throw(ArgumentError("$caller: hybrid lower UID $uid lacks its actual native record"))
        node_entities[row]=owner
        parameter=assembly.nodes[position].parameter
        parameters[row]=parameter===nothing ? nothing : copy(parameter)
    end
    edges=Dict{NTuple{2,Int32},Tuple{Int,Int32}}()
    faces=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    quads=Dict{NTuple{4,Int32},Tuple{Int,Int32}}()
    edge_conflicts=Set{NTuple{2,Int32}}();face_conflicts=Set{NTuple{3,Int32}}();quad_conflicts=Set{NTuple{4,Int32}}()
    owners=Dict{Int,Vector{Int32}}()
    for block in mesh.blocks
        block.msh in (4,7) || throw(ArgumentError("$caller: hybrid Volume has a foreign cell family"))
        ids=fill(Int32(tag),size(block.nodes,2));owners[Int(block.msh)]=ids
        _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
            Int(block.msh),block.nodes,ids)
    end
    for owner in sort!(collect(closure))
        owner[1]<3 || continue
        product=get(input.products,owner,nothing);product===nothing && continue
        for block in product.blocks
            mapped=Matrix{Int32}(undef,size(block.nodes));keep=trues(size(block.nodes,2))
            for column in axes(block.nodes,2),slot in axes(block.nodes,1)
                uid=data.external_node_tags[product.rows[block.nodes[slot,column]]]
                row=get(rows,uid,Int32(0));mapped[slot,column]=row
                if row==0
                    owner[1]==0 || throw(ArgumentError(
                        "$caller: hybrid Volume omits an actual lower cell vertex $uid"))
                    keep[column]=false
                end
            end
            cells=all(keep) ? mapped : mapped[:,keep]
            _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
                Int(block.msh),cells,fill(owner[2],size(cells,2)))
        end
    end
    _api_check_support_conflicts(edge_conflicts,face_conflicts,quad_conflicts)
    boundaries=Dict(owner=>Int32[abs(boundary) for boundary in data.entities[owner].boundaries] for owner in closure)
    return _mixed_classification(mesh,primary,sort!(collect(closure)),node_entities,boundaries,owners;
        edge_entities=edges,face_entities=faces,quad_entities=quads,node_parametric=parameters)
end

function _higher_volume_input_classification(assembly,m,mesh::Mesh,input,tag,
                                            source_tags,caller)
    length(source_tags)==nnodes(mesh) && allunique(source_tags) || throw(ArgumentError(
        "$caller: native Volume[$tag] lacks a unique complete row UID binding"))
    actual=input.mesh;data=actual.entity_data
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_tags))
    old_rows=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    closure=Model._model_mesh_projection_closure(m,3,tag,caller)
    primary=(3,Int32(tag))
    node_entities=fill(primary,nnodes(mesh))
    parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(mesh)]
    for (row,uid) in enumerate(source_tags)
        source=get(old_rows,uid,0)
        source==0 && continue
        owner=data.node_entities[source]
        owner in closure || throw(ArgumentError("$caller: Volume[$tag] UID $uid belongs outside its entity closure"))
        all(axis->isequal(actual.coords[axis,source],mesh.coords[axis,row]),1:3) ||
            throw(ArgumentError("$caller: Volume[$tag] UID $uid changed coordinates during projection"))
        node_entities[row]=owner
        position=get(assembly.node_positions,uid,Int32(0))
        position>0 || throw(ArgumentError("$caller: Volume[$tag] source UID $uid lacks its native record"))
        parameter=assembly.nodes[position].parameter
        parameters[row]=parameter===nothing ? nothing : copy(parameter)
    end
    # Exact lower cell consumption is the admission certificate, independently
    # of coordinate proximity, CAD classification, and constructor row order.
    geometry=(surface=mesh,native_node_tags=source_tags)
    constraints=Model._model_native_volume_input_constraints(m,tag,geometry,input,caller)
    Model._model_certify_native_volume_constraints(mesh,constraints,tag,caller)
    edges=Dict{NTuple{2,Int32},Tuple{Int,Int32}}()
    faces=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    quads=Dict{NTuple{4,Int32},Tuple{Int,Int32}}()
    edge_conflicts=Set{NTuple{2,Int32}}();face_conflicts=Set{NTuple{3,Int32}}();quad_conflicts=Set{NTuple{4,Int32}}()
    tetrahedra=fill(Int32(tag),ntets(mesh))
    _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,4,mesh.tets,tetrahedra)
    for owner in sort!(collect(closure))
        owner[1]<3 || continue
        product=get(input.products,owner,nothing)
        product===nothing && continue
        for block in product.blocks
            mapped=Matrix{Int32}(undef,size(block.nodes))
            keep=trues(size(block.nodes,2))
            for column in axes(block.nodes,2),slot in axes(block.nodes,1)
                uid=data.external_node_tags[product.rows[block.nodes[slot,column]]]
                row=get(rows,uid,Int32(0));mapped[slot,column]=row
                if row==0
                    # Unused Point records remain in the assembly; they are
                    # not added to Volume cells or aliased onto the first Point.
                    owner[1]==0 || throw(ArgumentError("$caller: Volume[$tag] omits an actual lower cell vertex $uid"))
                    keep[column]=false
                end
            end
            cells=all(keep) ? mapped : mapped[:,keep]
            _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,
                Int(block.msh),cells,fill(owner[2],size(cells,2)))
        end
    end
    _api_check_support_conflicts(edge_conflicts,face_conflicts,quad_conflicts)
    boundaries=Dict(owner=>Int32[abs(boundary) for boundary in data.entities[owner].boundaries] for owner in closure)
    owners=Dict{Int,Vector{Int32}}(1=>Int32[],2=>Int32[],4=>tetrahedra)
    return _mixed_classification(mesh,primary,sort!(collect(closure)),node_entities,boundaries,owners;
        edge_entities=edges,face_entities=faces,quad_entities=quads,node_parametric=parameters)
end

function _higher_ordinary_retained_surface_rows(production,input,tag,caller)
    owner=(2,Int32(tag));product=get(input.products,owner,nothing)
    product!==nothing && product.owner==owner || throw(ArgumentError(
        "$caller: lateral subdivision lost its actual Surface vertex storage"))
    lookup=copy(production.position_rows)
    extra=NTuple{3,Float64}[]
    for row in product.rows
        input.mesh.entity_data.node_entities[row]==owner || continue
        point=ntuple(axis->input.mesh.coords[axis,row],3)
        haskey(lookup,point) && continue
        push!(extra,point);lookup[point]=Int32(nnodes(production.mesh)+length(extra))
    end
    isempty(extra) && return production
    appended=Matrix{Float64}(undef,3,length(extra))
    for (column,point) in enumerate(extra),axis in 1:3
        appended[axis,column]=point[axis]
    end
    source=production.mesh;coords=hcat(source.coords,appended)
    mesh=source isa Mesh ? Mesh(coords;segs=source.segs,tris=source.tris,tets=source.tets,
        seg_tag=source.seg_tag,tri_tag=source.tri_tag,tet_tag=source.tet_tag) : MixedMesh(coords,source.blocks)
    return Model._ExtrudeMeshProduction(mesh,lookup,production.source_node_tags,
        production.column_rows,production.cell_events)
end

function _higher_ordinary_lateral_plans(assembly,m,input,extruded,caller)
    plans=Dict{Int,Any}()
    extruded===nothing && return plans
    for tag in sort!(collect(extruded.remesh))
        entry=Model._extrude_entity_params(m,2,tag,caller)
        entry===nothing && throw(ArgumentError(
            "$caller: lateral subdivision plan lost its extrusion recipe"))
        params,spec,link=entry
        params.quad_to_tri===:none || continue
        link[1]==1 || throw(ArgumentError(
            "$caller: lateral subdivision must consume an actual Curve Source"))
        production=Model._extrude_lateral_mesh(m,tag,params,spec,link[2],caller;
            edges=extruded.edges,_mesh_input=input,_return_production=true,_register_orphans=false)
        # The constrained native remesh deletes cells and preserves the whole
        # Surface vertex vector. It does not create new Curve-orphan children.
        production=_higher_ordinary_retained_surface_rows(production,input,tag,caller)
        certified=_higher_ordinary_sweep_classification(assembly,m,production,input,2,tag,caller;
            primary_existing=true)
        certified!==nothing && all(certified.retained) || throw(ArgumentError(
            "$caller: lateral subdivision lacks complete actual UID retention"))
        plans[tag]=(owner=(2,Int32(tag)),input=input,production=production,
            node_tags=certified.node_tags,class=certified.class)
    end
    return plans
end

function _higher_ordinary_volume_events(production,caller)
    expected=Set((msh,column) for (msh,cells,_) in _cache_native_blocks(production.mesh) if msh_dimension(msh)==3
        for column in axes(cells,2))
    events=production.cell_events
    length(events)==length(expected) && allunique(events) && Set(events)==expected ||
        throw(ArgumentError("$caller: ordinary Volume constructor events do not cover its actual cells"))
    return events
end

function _higher_ordinary_volume_cells!(assembly,m,production,class,mapping,geometry_cells)
    mesh=production.mesh
    catalog=_cache_catalog(mesh,class)
    blocks=_cache_native_blocks(mesh)
    events=_higher_ordinary_volume_events(production,assembly.caller)
    tags=Dict{Tuple{Int,Int},UInt64}()
    for event in events
        tags[event]=_dim01_next_tag!(assembly,false)
    end
    traces=geometry_cells===nothing ? nothing : [Union{Nothing,_Dim01Cell}[
        nothing for _ in axes(cells,2)] for (_,cells,_) in blocks]
    for (bi,(msh,dim,_,cells,owners)) in enumerate(catalog)
        dim==3 || continue
        physical=blocks[bi][3]
        for column in axes(cells,2)
            owner=(dim,owners[column])
            value=mesh isa MixedMesh && mesh.entity_data!==nothing ? physical[column] :
                _dim01_physical(m,(dim,Int(owners[column])))
            _dim01_add_cell!(assembly,msh,owner,mapping[cells[:,column]],tags[(msh,column)];
                source_tag=nothing,physical=value)
            traces===nothing || (traces[bi][column]=last(assembly.cells))
        end
    end
    geometry_cells===nothing || (geometry_cells[mesh]=traces)
    return nothing
end

function _higher_ordinary_volume_stage!(assembly,m,input,extruded,meshed,lateral_plans,
        discarded,done,emitted;source_nodes,geometry_parts,geometry_cells,native_orphan_births)
    # All prospective cell and vertex incidence is certified before allocating
    # this stage. The detached operation later certifies the final lower cells.
    prepared=Dict{Int,Any}()
    for tag in sort!(copy(meshed))
        production=extruded.productions[tag]
        certified=_higher_ordinary_sweep_classification(assembly,m,production,input,3,tag,
            assembly.caller;planned_laterals=lateral_plans)
        certified!==nothing || throw(ArgumentError(
            "$(assembly.caller): ordinary subdivision lacks a complete Volume certificate"))
        _higher_ordinary_volume_events(production,assembly.caller)
        if !Model._extrude_gate(m,3,tag).recombine
            count=get(extruded.temporary_elements,tag,nothing)
            count isa Int && count>=0 || throw(ArgumentError(
                "$(assembly.caller): subdivision lost consumed initial Volume cell events"))
        end
        prepared[tag]=(production=production,certificate=certified)
    end
    mappings=Dict{Int,Vector{Int32}}()
    for tag in sort!(copy(meshed))
        item=prepared[tag];certified=item.certificate
        mappings[tag]=_higher_append_part!(assembly,m,item.production.mesh,certified.class,3,tag,discarded;
            emitted,source_nodes,source_row_tags=certified.node_tags,
            source_row_retained=certified.retained,source_rows_existing_only=true,
            native_orphan_births,emit_cells=false)
    end
    # meshGRegionExtruded visits Regions in tag order. Recombined cell births
    # survive; non-recombined initial cells are deleted by subdivision.
    for tag in sort!(copy(meshed))
        item=prepared[tag]
        if Model._extrude_gate(m,3,tag).recombine
            _higher_ordinary_volume_cells!(assembly,m,item.production,item.certificate.class,
                mappings[tag],geometry_cells)
        else
            count=get(extruded.temporary_elements,tag,nothing)
            count!==nothing || throw(ArgumentError(
                "$(assembly.caller): subdivision lost consumed initial Volume cell events"))
            for _ in 1:count
                _dim01_next_tag!(assembly,false)
            end
        end
    end
    # SubdivideExtrudedMesh then allocates final cells for each pending Region.
    for tag in sort!(copy(meshed))
        Model._extrude_gate(m,3,tag).recombine && continue
        item=prepared[tag]
        _higher_ordinary_volume_cells!(assembly,m,item.production,item.certificate.class,
            mappings[tag],geometry_cells)
    end
    for tag in sort!(copy(meshed))
        item=prepared[tag]
        geometry_parts===nothing || push!(geometry_parts,
            (tag,item.production.mesh,item.certificate.class,mappings[tag]))
        push!(done,(3,Int32(tag)))
    end
    return nothing
end

# Fresh Surface cells are actual event objects. Their final labels may change
# during finish! even though they have no retained source_tag mapping entry.
function _higher_surface_constructor_publication(authority,constructors,cells,node_map,caller)
    rebound=Dict{Tuple{Int,Int32},MixedMesh}()
    for (owner,source) in constructors
        constructor=deepcopy(source);data=constructor.entity_data
        data!==nothing || throw(ArgumentError(
            "$caller: Surface constructor publication lacks its actual UID metadata"))
        nodes=UInt64[]
        kept_rows=Int32[]
        for (row,uid) in enumerate(data.external_node_tags)
            final=get(node_map,uid,nothing)
            # finish! removes unreferenced product vertices whose owner was
            # regenerated; such rows cannot support any actual cell.
            final===nothing && continue
            (final>0 || (authority.mesh.entity_data.native_mesh_tags &&
                final in authority.mesh.entity_data.external_node_tags)) || throw(ArgumentError(
                "$caller: Surface constructor publication loses actual node UID $uid"))
            push!(nodes,final)
            push!(kept_rows,Int32(row))
        end
        allunique(nodes) || throw(ArgumentError(
            "$caller: Surface constructor publication aliases actual node UIDs"))
        if length(kept_rows)!=nnodes(constructor)
            row_map=zeros(Int32,nnodes(constructor))
            row_map[kept_rows]=Int32.(1:length(kept_rows))
            blocks=ElementBlock[]
            for block in constructor.blocks
                remapped=map(row->row_map[row],block.nodes)
                all(row->row>0,remapped) || throw(ArgumentError(
                    "$caller: Surface constructor publication loses a supported node UID"))
                push!(blocks,ElementBlock(block.msh,remapped,block.tags))
            end
            constructor=MixedMesh(constructor.coords[:,kept_rows],blocks;
                physical_names=constructor.physical_names,
                entity_data=Elements._native_mixed_entity_data(data.entities;
                    node_entities=data.node_entities[kept_rows],
                    node_parametric=data.node_parametric[kept_rows],
                    external_node_tags=nodes,
                    block_entities=deepcopy(data.block_entities),
                    external_element_tags=deepcopy(data.external_element_tags)),
                elementary_entities=constructor.elementary_entities)
            data=constructor.entity_data
        else
            data.external_node_tags[:]=nodes
        end
        for labels in data.external_element_tags
            final=UInt64[]
            for uid in labels
                event=get(cells,uid,nothing)
                event!==nothing && event.owner==owner && (event.tag>0 ||
                    (authority.mesh.entity_data.native_mesh_tags &&
                     any(labels->event.tag in labels,authority.mesh.entity_data.external_element_tags))) || throw(ArgumentError(
                    "$caller: Surface constructor publication loses its owned actual cell event"))
                push!(final,event.tag)
            end
            allunique(final) || throw(ArgumentError(
                "$caller: Surface constructor publication aliases actual cell UIDs"))
            labels[:]=final
        end
        rebound[owner]=constructor
    end
    return Model._model_mesh_with_surface_constructor_inputs(authority,rebound;caller)
end

function _higher_retained_periodic_surface_pairs!(pairs,m,retained,caller;slave_tags=nothing)
    pairs===nothing && return nothing
    for constraint in Model.model_periodic_constraints(m)
        constraint.dim==2 || continue
        slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
        slave_tags===nothing || slave in slave_tags || continue
        haskey(pairs,slave) && continue
        (2,Int32(slave)) in retained && (2,Int32(master)) in retained || continue
        history=get(m.meshing.periodic_uid_history,(2,slave),nothing)
        history===nothing && continue
        history.master==master && history.affine==constraint.affine || throw(ArgumentError(
            "$caller: retained periodic Surface history changes its declared relation"))
        all(pair->pair.slave_live && pair.master_live,history.primary) || throw(ArgumentError(
            "$caller: retained periodic Surface copy history has retired constructor UIDs"))
        # These are already published live copy events, not a new geometric
        # correspondence. The hybrid admission certificate binds them to the
        # current retained constructor cells and their exact UID incidence.
        pairs[slave]=(master,[(pair.slave,pair.master) for pair in history.primary])
    end
    return nothing
end

function _higher_volume_phase!(assembly,m,options,entities,discarded,done;
                               size_field=nothing,source_nodes=length(assembly.nodes),ranks=nothing,
                               old_supports=nothing,geometry_parts=nothing,geometry_cells=nothing,
                               geometry_birth=nothing,
                               surface_parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[],
                               native_orphan_births=nothing,
                               periodic_surface_pairs=nothing,
                               optimize_pyramids::Float64=0.0,
                               retained_surfaces=Set{Tuple{Int,Int32}}(),
                               surface_constructor_inputs=Dict{Tuple{Int,Int32},MixedMesh}())
    have=Set(cell.owner for cell in assembly.cells)
    meshed=Int[]
    for (dim,tag) in entities
        dim==3 || continue
        owner=(3,Int32(tag))
        visible=!(options.mesh_only_visible && model_entity_visibility(m,3,tag)==0)
        # Native extrusion is scheduled before meshGRegion's OnlyEmpty gate.
        swept=Model._extrude_gate(m,3,tag)!==nothing
        keep=_dim01_fully_discrete(m,(3,tag)) ||
             (!swept && (!visible || (options.mesh_only_empty && owner in have)))
        keep && continue
        push!(meshed,tag);push!(discarded,owner)
    end
    filter!(cell->cell.source_tag===nothing || !(cell.owner in discarded),assembly.cells)
    ranks===nothing || _higher_complete_ranks!(ranks,assembly)
    input_rows=_higher_assembly_authority(assembly,m,done;ranks,old_supports)
    input_bound=_higher_bind_producer_rows(input_rows,assembly,surface_parts,m.meshing.mesh_authority)
    input=Model._model_mesh_with_surface_constructor_inputs(input_bound,surface_constructor_inputs;
        caller=assembly.caller)
    scope=Model._extrude_nonew_scope(m,assembly.caller;tags=meshed,
        size_field,_working_model=true,_mesh_input=input,
        _trace_columns=geometry_parts!==nothing)
    extruded=Model._extrude_volume_pass(m,assembly.caller;
        scope,_mesh_input=input,volume_tags=meshed)
    lateral_plans=_higher_ordinary_lateral_plans(assembly,m,input,extruded,assembly.caller)
    emitted=Set{Tuple{Int,Int32}}()
    staged=extruded!==nothing && any(tag->haskey(extruded.productions,tag) &&
        !Model._extrude_gate(m,3,tag).recombine,meshed)
    if staged
        all(tag->haskey(extruded.productions,tag),meshed) || throw(ArgumentError(
            "$(assembly.caller): ordinary subdivision combined with another Region planner lacks a complete birth plan"))
        _higher_ordinary_volume_stage!(assembly,m,input,extruded,meshed,lateral_plans,
            discarded,done,emitted;source_nodes,geometry_parts,geometry_cells,native_orphan_births)
    end
    hybrid_batch=nothing;hybrid_batch_receipt=nothing
    retained_periodic=periodic_surface_pairs!==nothing && any(constraint->
        constraint.dim==2 && !haskey(periodic_surface_pairs,Int(constraint.slave_entity)) &&
        (2,Int32(constraint.slave_entity)) in retained_surfaces &&
        (2,Int32(constraint.master_entity)) in retained_surfaces,
        Model.model_periodic_constraints(m))
    if retained_periodic && m.meshing.algorithm_3d==1 && !staged && extruded===nothing &&
            length(meshed)==length(m.volumes)==1
        tag=only(meshed)
        if !haskey(m.meshing.transfinite_volumes,tag) && !isempty(m.volumes[tag]) &&
                !Model._model_factory_volume_encoding(m,tag) &&
                Model._model_native_hybrid_volume_input(m,tag,input,assembly.caller)!==nothing
            _higher_retained_periodic_surface_pairs!(periodic_surface_pairs,m,retained_surfaces,assembly.caller)
        end
    end
    if !staged && extruded===nothing && length(meshed)>1 && all(tag->
            !haskey(m.meshing.transfinite_volumes,tag) &&
            !isempty(m.volumes[tag]) && !Model._model_factory_volume_encoding(m,tag) &&
            Model._model_native_hybrid_volume_input(m,tag,input,assembly.caller)!==nothing,meshed)
        if retained_periodic && m.meshing.algorithm_3d==1
            # Only selected hybrid Regions consume retained copy events here.
            # Other retained model relations keep their own planner admission.
            slave_tags=Set(Int(constraint.slave_entity) for tag in meshed
                for constraint in Model._model_volume_periodic_surface_constraints(m,tag,assembly.caller))
            _higher_retained_periodic_surface_pairs!(periodic_surface_pairs,m,retained_surfaces,
                assembly.caller;slave_tags)
        end
        births=Model._ModelNativeVolumeBirths(assembly.node_max,assembly.cell_max)
        hybrid_batch=Model._model_native_hybrid_batch_production(m,meshed,input,options,births,
            assembly.caller;size_field,optimize_pyramids,periodic_surface_pairs)
        hybrid_batch_receipt=_higher_native_hybrid_batch_receipt(assembly,m,hybrid_batch,input,
            meshed,assembly.caller)
    end
    for tag in (staged ? Int[] : sort!(meshed))
        mesh=extruded===nothing ? nothing : get(extruded.parts,tag,nothing)
        source_rows=nothing;source_row_retained=nothing;source_rows_existing_only=false
        source_cells=nothing;source_cell_max=nothing;alternate=nothing;source_emission=nothing;hybrid_production=nothing
        source_node_base=hybrid_batch_receipt===nothing ? assembly.node_max : hybrid_batch_receipt.node_base
        source_node_max=nothing
        if mesh!==nothing && geometry_parts!==nothing && scope!==nothing
            completed=get(scope.volumes,tag,nothing)
            completed===nothing || (alternate=Model._extrude_nonew_b4_geometry(m,tag,
                completed,uid->_higher_geometry_source_point(assembly,geometry_birth,uid),assembly.caller))
        end
        if mesh===nothing
            production=hybrid_batch===nothing ? mesh_model_volume(m,tag;size_field,
                _extrude_scope=scope,_extrude_working_model=true,_mesh_input=input,
                _mesh1d_options=options,_return_volume_production=true,
                _periodic_surface_pairs=periodic_surface_pairs,
                _optimize_pyramids=optimize_pyramids,
                _mesh_node_high_water=assembly.node_max,_mesh_element_high_water=assembly.cell_max,
                _geometry_points=_higher_geometry_source_points(assembly,geometry_birth)) : hybrid_batch.products[tag]
            production isa Union{Model._ModelVolumeGeometryProduction,Model._ModelNativeHybridVolumeProduction} &&
                _higher_native_region_surface_order!(assembly,m,tag)
            mesh=production.mesh
            if production isa Model._ModelNativeHybridVolumeProduction
                hybrid_production=production
            end
            source_node_max=production.node_high_water
            if isempty(production.element_tags)
                if production.temporary_elements>0
                    assembly.node_max=max(assembly.node_max,production.node_high_water)
                    assembly.cell_max=Base.checked_add(assembly.cell_max,production.temporary_elements)
                end
            else
                source_cells=production.element_tags
                source_cell_max=production.element_high_water
            end
            if production isa Model._ModelTransfiniteVolumeProduction
                source_rows=production.boundary_node_tags
                source_row_retained=production.boundary_node_retained
                source_rows_existing_only=true
                alternate=production.geometry
                source_emission=production.emission
            else
                any(!iszero,production.boundary_node_tags) &&
                    (source_rows=production.boundary_node_tags)
            end
        end
        native_rows=mesh isa Mesh && source_rows!==nothing &&
            all(!iszero,source_rows) && allunique(source_rows)
        ordinary=extruded===nothing ? nothing : get(extruded.productions,tag,nothing)
        swept=ordinary===nothing ? nothing : _higher_ordinary_sweep_classification(
            assembly,m,ordinary,input,3,tag,assembly.caller;planned_laterals=lateral_plans)
        swept===nothing && (swept=_higher_nonew_identity_classification(
            assembly,m,mesh,input,scope,3,tag,assembly.caller))
        if swept!==nothing
            source_rows=swept.node_tags;source_row_retained=swept.retained
            source_rows_existing_only=true
        end
        class=if swept!==nothing
            swept.class
        elseif hybrid_production!==nothing
            _higher_hybrid_volume_input_classification(assembly,m,hybrid_production,input,tag,assembly.caller)
        elseif source_rows_existing_only
            _higher_transfinite_volume_classification(assembly,m,mesh,input,tag,
                source_rows,source_row_retained,assembly.caller;grid_receipt=source_emission.grid)
        elseif native_rows
            _higher_volume_input_classification(assembly,m,mesh,input,tag,source_rows,assembly.caller)
        else
            _classify_cached_mesh(m,_classification_skeleton(
                mesh,get(m.meshing.reverse,(3,tag),false)),3,tag,mesh;_extrude_scope=scope,
                _mesh_input=source_rows!==nothing ? input : nothing)
        end
        if source_rows_existing_only && swept===nothing
            class===nothing && throw(ArgumentError("$(assembly.caller): actual TF Volume lacks its shell classification"))
            for row in eachindex(source_rows)
                if source_row_retained[row]
                    actual=get(assembly.node_positions,source_rows[row],Int32(0))
                    actual>0 || throw(ArgumentError("$(assembly.caller): TF boundary source node is absent"))
                    class.node_entities[row]=assembly.nodes[actual].owner
                else
                    # Logical interior events belong to the Region even when
                    # a compact prism stores a point on an aliased diagonal.
                    class.node_entities[row]=(3,Int32(tag))
                end
            end
        end
        mapping=_higher_append_part!(assembly,m,mesh,class,3,tag,discarded;
            emitted,source_nodes,geometry_cells,source_emission,source_row_tags=source_rows,
            source_row_retained,source_rows_existing_only,
            source_node_base,source_node_high_water=source_node_max,
            source_cell_tags=source_cells,source_cell_high_water=source_cell_max,
            source_batch_receipt=hybrid_batch_receipt,
            source_row_order=hybrid_production===nothing ? nothing : hybrid_production.registration_rows,
            source_volume_production=hybrid_production,
            # Recorded TF producers can create unused Region-owned logical rows.
            # Keep their born identities through finish, like ordinary sweep events.
            native_orphan_births=ordinary===nothing && source_emission===nothing ?
                nothing : native_orphan_births)
        if alternate!==nothing
            geometry_cells===nothing || (geometry_cells[alternate]=geometry_cells[mesh])
            if geometry_birth!==nothing
                for (row,position) in enumerate(mapping)
                    node=assembly.nodes[position]
                    point=(alternate.coords[1,row],alternate.coords[2,row],alternate.coords[3,row])
                    geometry_birth.nodes[node.tag]=_Dim01Node(node.tag,point,node.owner,node.parameter)
                end
            end
        end
        if geometry_parts!==nothing
            geometry_mesh,geometry_class,geometry_mapping=
                _higher_volume_geometry_order(alternate===nothing ? mesh : alternate,
                                               class,mapping,m,tag,input,scope;
                                               geometry_cells,caller=assembly.caller)
            push!(geometry_parts,(tag,geometry_mesh,geometry_class,geometry_mapping))
        end
        push!(done,(3,Int32(tag)))
    end
    laterals=extruded===nothing ? Set{Int}() : copy(extruded.remesh)
    if scope!==nothing
        for surface in keys(scope.surfaces)
            link=get(m.meshing.extrude_sources,(2,surface),nothing)
            link!==nothing && link[1]==1 && push!(laterals,surface)
        end
    end
    # Subdivision remeshes lateral cells after allocating the volume. The
    # existing lateral vertex identities remain available to this cell pass.
    filter!(cell->!(cell.owner[1]==2 && Int(cell.owner[2]) in laterals),assembly.cells)
    lateral_owners=Set((2,Int32(tag)) for tag in laterals)
    for tag in sort!(collect(laterals))
        plan=get(lateral_plans,tag,nothing)
        mesh=plan===nothing ? mesh_model_surface(m,tag;size_field,_extrude_scope=scope,
            _extrude_edges=extruded===nothing ? nothing : extruded.edges,_mesh1d_options=options) :
            plan.production.mesh
        swept=plan===nothing ? _higher_nonew_identity_classification(
            assembly,m,mesh,input,scope,2,tag,assembly.caller) : nothing
        class=plan===nothing ? (swept===nothing ? _classify_cached_mesh(
            m,_classification_skeleton(mesh,false),2,tag,mesh;_extrude_scope=scope) : swept.class) : plan.class
        mapping=_higher_append_part!(assembly,m,mesh,class,2,tag,discarded;
            emitted=lateral_owners,source_nodes,geometry_cells,
            source_row_tags=plan===nothing ? (swept===nothing ? nothing : swept.node_tags) : plan.node_tags,
            source_row_retained=plan===nothing ? (swept===nothing ? nothing : swept.retained) : trues(length(plan.node_tags)),
            source_rows_existing_only=plan!==nothing || swept!==nothing)
        push!(surface_parts,(tag,mesh,class,mapping))
        push!(done,(2,Int32(tag)))
    end
    if !isempty(lateral_plans)
        ranks===nothing || _higher_complete_ranks!(ranks,assembly)
        final_input=_higher_assembly_authority(assembly,m,done;ranks,old_supports)
        for (tag,production) in extruded.productions
            # The prospective shell certified allocation order; the actual
            # newly born lower cells must independently certify publication.
            _higher_ordinary_sweep_classification(assembly,m,production,final_input,3,tag,
                assembly.caller)
        end
    end
    return nothing
end

# Actual Source/private ranks drive the extrusion decisions. The geometry view
# separately keeps the original Surface producer's column insertion order. The
# renderer records each real inserted/deduplicated row, so this reorder never
# selects identities by coordinates or assumes one new row per column event.
function _higher_volume_geometry_order(mesh,class,mapping,m,tag,input,scope;
                                       geometry_cells,caller)
    scope===nothing && return mesh,class,mapping
    completed=get(scope.volumes,tag,nothing)
    completed===nothing && return mesh,class,mapping
    trace=completed.column_rows
    isempty(trace) && return mesh,class,mapping
    size(trace)==size(completed.columns) &&
        all(row->0<row<=nnodes(mesh),trace) || throw(ArgumentError(
            "$caller: extrusion column events do not bind actual geometry rows"))
    link=get(m.meshing.extrude_sources,(3,tag),nothing)
    link===nothing && return mesh,class,mapping
    owner=(2,Int32(abs(link[2])))
    produced=get(input.producer_rows,owner,nothing)
    produced===nothing && return mesh,class,mapping
    source=get(completed.surfaces,Int(owner[2]),nothing)
    source isa MixedMesh && source.entity_data!==nothing || throw(ArgumentError(
        "$caller: extrusion Source lacks actual column identities"))
    source_tags=source.entity_data.external_node_tags
    length(source_tags)==size(trace,2) && allunique(source_tags) || throw(ArgumentError(
        "$caller: extrusion Source column identities do not align with emission events"))
    columns=Dict(uid=>column for (column,uid) in enumerate(source_tags))
    produced_tags=input.mesh.entity_data.external_node_tags[produced]
    # A previous producer can be partial after an operation. Keep the current
    # valid emission in that case; it supplies no complete geometry permutation.
    length(produced_tags)==length(source_tags) && all(uid->haskey(columns,uid),produced_tags) ||
        return mesh,class,mapping
    rows=Int32[];sizehint!(rows,nnodes(mesh));seen=falses(nnodes(mesh))
    for uid in produced_tags,level in axes(trace,1)
        row=trace[level,columns[uid]]
        if !seen[row]
            push!(rows,row);seen[row]=true
        end
    end
    for row in eachindex(seen)
        seen[row] || push!(rows,Int32(row))
    end
    length(rows)==nnodes(mesh) && allunique(rows) || throw(ArgumentError(
        "$caller: extrusion geometry row permutation is not bijective"))
    all(i->rows[i]==i,eachindex(rows)) && return mesh,class,mapping
    inverse=zeros(Int32,length(rows));inverse[rows]=Int32.(eachindex(rows))
    geometry=if mesh isa Mesh
        Mesh(mesh.coords[:,rows];segs=inverse[mesh.segs],tris=inverse[mesh.tris],
             tets=inverse[mesh.tets],seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
    else
        mesh.entity_data===nothing || throw(ArgumentError(
            "$caller: extrusion geometry unexpectedly contains public metadata"))
        blocks=ElementBlock[ElementBlock(block.msh,inverse[block.nodes],block.tags)
                            for block in mesh.blocks]
        MixedMesh(mesh.coords[:,rows],blocks;physical_names=mesh.physical_names,
            elementary_entities=mesh.elementary_entities,periodic_links=mesh.periodic_links,
            ancillary_sections=mesh.ancillary_sections,data_sections=mesh.data_sections,
            partition_data=mesh.partition_data)
    end
    geometry_cells===nothing || (geometry_cells[geometry]=geometry_cells[mesh])
    geometry_class=geometry isa MixedMesh ?
        _mixed_rebind_class(class,geometry;node_entities=class.node_entities[rows],
                            primary_map=inverse,node_parametric=class.node_parametric[rows]) :
        _higher_geometry_part_class(geometry,3,tag,class.node_entities[rows])
    return geometry,geometry_class,mapping[rows]
end

function _higher_bind_producer_rows(authority,assembly,parts,previous;node_tag_map=nothing,
                                    curve_births=nothing)
    positions=Dict(tag=>Int32(row) for (row,tag) in
                   enumerate(authority.mesh.entity_data.external_node_tags))
    produced=Dict{Tuple{Int,Int32},Vector{Int32}}()
    replaced=Set((2,Int32(tag)) for (tag,_,_,_) in parts)
    curve_births===nothing || union!(replaced,keys(curve_births))
    if previous!==nothing && previous.mesh!==nothing
        old_tags=previous.mesh.entity_data.external_node_tags
        for (owner,old_rows) in previous.producer_rows
            owner in replaced && continue
            haskey(authority.products,owner) || continue
            rows=Int32[]
            for old_row in old_rows
                source_tag=old_tags[old_row]
                final_tag=node_tag_map===nothing ? source_tag :
                    get(node_tag_map,source_tag,nothing)
                row=get(positions,final_tag,Int32(0))
                row>0 || throw(ArgumentError(
                    "$(assembly.caller): retained producer $owner loses actual node $source_tag"))
                push!(rows,row)
            end
            produced[owner]=rows
        end
    end
    for (tag,mesh,_,mapping) in parts
        rows=Vector{Int32}(undef,length(mapping))
        for (slot,position) in enumerate(mapping)
            source_tag=assembly.nodes[position].tag
            final_tag=node_tag_map===nothing ? source_tag :
                get(node_tag_map,source_tag,nothing)
            row=get(positions,final_tag,Int32(0))
            row>0 || throw(ArgumentError(
                "$(assembly.caller): produced Surface[$tag] loses actual node $source_tag"))
            rows[slot]=row
        end
        authority.mesh.coords[:,rows]==mesh.coords || throw(ArgumentError(
            "$(assembly.caller): produced Surface[$tag] coordinates disagree with actual bound rows"))
        owner=(2,Int32(tag))
        if !haskey(authority.products,owner) && Model._surface_cell_count(mesh)==0
            authority=Model._model_mesh_with_empty_surface_product(authority,owner,rows;
                caller=assembly.caller)
        end
        produced[owner]=rows
    end
    return Model._model_mesh_with_producer_rows(authority,produced;caller=assembly.caller)
end

function _higher_geometry_part_class(mesh,dimension,tag,node_entities)
    segments=fill(Int32(tag),nsegs(mesh));triangles=fill(Int32(tag),ntris(mesh))
    tetrahedra=fill(Int32(tag),ntets(mesh))
    return _MeshClassification(mesh,(dimension,Int32(tag)),[(dimension,Int32(tag))],
        node_entities,Dict{Tuple{Int,Int32},Vector{Int32}}(),segments,triangles,tetrahedra,
        Dict{Int,Vector{Int32}}(1=>segments,2=>triangles,4=>tetrahedra),nothing,
        Dict{NTuple{2,Int32},Tuple{Int,Int32}}(),
        Dict{NTuple{3,Int32},Tuple{Int,Int32}}(),
        Dict{NTuple{4,Int32},Tuple{Int,Int32}}())
end

# Retained Triangle geometry needs local identities, not a metadata projection.
# A local row inverse keeps producer-order extraction proportional to its rows.
function _higher_surface_geometry_input(authority,owner,caller)
    product=get(authority.products,owner,nothing)
    product===nothing && return nothing
    any(block->!isempty(block.nodes),product.blocks) || return nothing
    all(block->block isa ElementBlock && block.msh==2,product.blocks) || return nothing
    source=authority.mesh
    data=source.entity_data
    produced=get(authority.producer_rows,owner,nothing)
    rows=copy(produced===nothing ? product.rows : produced)
    # An AddNodes event can append an owned orphan after the original
    # producer receipt. It has no cell incidence to reorder, but remains an
    # actual source identity in the geometry carrier. Missing used rows must
    # still fail rather than manufacturing a replacement producer order.
    used=Set(product.rows[row] for block in product.blocks for row in block.nodes)
    if produced!==nothing
        for row in product.rows
            row in rows && continue
            !(row in used) && data.node_entities[row]==owner || throw(ArgumentError(
                "$caller: actual entity $owner lost a used producer row"))
            push!(rows,row)
        end
    end
    inverse=Dict{Int32,Int32}();sizehint!(inverse,length(rows))
    for (position,row) in enumerate(rows)
        haskey(inverse,row) && throw(ArgumentError(
            "$caller: actual entity $owner repeats a producer row"))
        inverse[row]=Int32(position)
    end
    all(row->haskey(inverse,row),product.rows) || throw(ArgumentError(
        "$caller: actual entity $owner changed without a complete producer row binding"))
    count=sum(size(block.nodes,2) for block in product.blocks;init=0)
    triangles=Matrix{Int32}(undef,3,count)
    physical=Vector{Int32}(undef,count)
    column=0
    for block in product.blocks,index in axes(block.nodes,2)
        column+=1
        for slot in 1:3
            triangles[slot,column]=inverse[product.rows[block.nodes[slot,index]]]
        end
        physical[column]=block.tags[index]
    end
    geometry=Mesh(source.coords[:,rows];tris=triangles,tri_tag=physical)
    return (mesh=geometry,node_tags=copy(data.external_node_tags[rows]),
        node_entities=copy(data.node_entities[rows]),
        element_tags=[copy(tags) for tags in product.element_tags])
end

function _higher_primal_surface_geometry(mesh)
    any(block->!isempty(block.nodes),mesh.blocks) || return nothing
    all(block->block isa ElementBlock && block.msh==2,mesh.blocks) || return nothing
    triangles=hcat((block.nodes for block in mesh.blocks)...)
    tags=Int32[tag for block in mesh.blocks for tag in block.tags]
    return Mesh(mesh.coords;tris=triangles,tri_tag=tags)
end

function _higher_lower_geometry_parts(assembly,m;sort_points=false,geometry_cells=nothing,
                                       geometry_birth=nothing)
    parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[]
    points=Dict{Int,Vector{Int32}}()
    curves=Dict{Int,Vector{_Dim01Cell}}()
    for (row,node) in enumerate(assembly.nodes)
        node.owner[1]==0 || continue
        push!(get!(()->Int32[],points,Int(node.owner[2])),Int32(row))
    end
    for cell in assembly.cells
        cell.owner[1]==1 && cell.msh==1 || continue
        push!(get!(()->_Dim01Cell[],curves,Int(cell.owner[2])),cell)
    end
    point_order=Int[tag for tag in keys(m.points) if haskey(points,tag)]
    append!(point_order,sort!(Int[tag for tag in keys(points) if !haskey(m.points,tag)]))
    sort_points && sort!(point_order)
    ordered=vcat([(0,tag) for tag in point_order],[(1,tag) for tag in keys(curves)])
    for (dimension,tag) in ordered
        rows=dimension==0 ? points[tag] : unique!(Int32[
            row for cell in curves[tag] for row in cell.nodes])
        if dimension==1
            used=Set(rows)
            append!(rows,Int32[row for (row,node) in enumerate(assembly.nodes)
                if node.owner==(1,Int32(tag)) && !(Int32(row) in used)])
        end
        local_rows=Dict(row=>Int32(position) for (position,row) in enumerate(rows))
        coords=Matrix{Float64}(undef,3,length(rows))
        for (position,row) in enumerate(rows),axis in 1:3
            actual=assembly.nodes[row]
            alternate=geometry_birth===nothing ? actual :
                get(geometry_birth.nodes,actual.tag,actual)
            coords[axis,position]=(alternate.owner==actual.owner ? alternate : actual).point[axis]
        end
        cells=dimension==0 ? _Dim01Cell[] : curves[tag]
        segments=Matrix{Int32}(undef,2,length(cells))
        for (column,cell) in enumerate(cells)
            segments[:,column]=Int32[local_rows[row] for row in cell.nodes]
        end
        mesh=Mesh(coords;segs=segments,seg_tag=Int32[cell.physical for cell in cells])
        class=_higher_geometry_part_class(mesh,dimension,tag,
            [assembly.nodes[row].owner for row in rows])
        push!(parts,(tag,mesh,class,rows))
        geometry_cells===nothing || (geometry_cells[mesh]=[
            Union{Nothing,_Dim01Cell}[cells...],Union{Nothing,_Dim01Cell}[],
            Union{Nothing,_Dim01Cell}[]])
    end
    return parts
end

function _higher_geometry_carrier(plan,assembly,parts,dimension,caller;
                                  geometry_merge=nothing,geometry_cells,geometry_birth=nothing)
    isempty(parts) && return nothing,Int32[],Vector{UInt64}[],nothing,nothing
    geometry_parts=Tuple{Int,Union{Mesh,MixedMesh}}[(tag,mesh) for (tag,mesh,_,_) in parts]
    node_entities=[class.node_entities for (_,_,class,_) in parts]
    node_identities=[mapping for (_,_,_,mapping) in parts]
    live=IdDict(cell=>nothing for cell in assembly.cells)
    cell_identities=Vector{Vector{Union{Nothing,UInt64}}}[]
    for (_,mesh,_,_) in parts
        traces=get(geometry_cells,mesh,nothing)
        traces===nothing && all(isempty(cells) for (_,cells,_) in _cache_native_blocks(mesh)) &&
            (traces=[Union{Nothing,_Dim01Cell}[] for _ in _cache_native_blocks(mesh)])
        traces===nothing && throw(ArgumentError("$caller: produced geometry lacks its cell events"))
        push!(cell_identities,[Union{Nothing,UInt64}[cell===nothing ? nothing : haskey(live,cell) ? cell.tag :
            throw(ArgumentError("$caller: produced geometry cell loses its actual identity"))
            for cell in block] for block in traces])
    end
    if geometry_merge===nothing
        carrier,remaps,_,carrier_cells=_merge_entity_meshes(geometry_parts,caller;node_entities,
            node_identities,cell_identities,
            require_positive_tets=!any(part->part[3].entity[1]==3 &&
                get(plan.model.meshing.reverse,(3,part[1]),false),parts))
    else
        products=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
            (class.entity[1],tag,mesh) for (tag,mesh,class,_) in parts]
        carrier,remaps,carrier_cells=geometry_merge(products,node_entities,node_identities,
            cell_identities,caller)
    end
    # The real native producer merge preserves its geometry storage and CRC;
    # the complete public products remain in the independently bound cache.
    positions=Dict(tag=>Int32(row) for (row,tag) in
                   enumerate(plan.mesh.entity_data.external_node_tags))
    rows=zeros(Int32,nnodes(carrier))
    for ((_,_,_,mapping),merged_rows) in zip(parts,remaps)
        for (position,merged_row) in enumerate(merged_rows)
            source_tag=assembly.nodes[mapping[position]].tag
            final_tag=get(plan.node_tag_map,source_tag,nothing)
            actual_row=get(positions,final_tag,Int32(0))
            actual_row>0 || throw(ArgumentError("$caller: produced geometry vertex loses its actual identity"))
            rows[merged_row]==0 || rows[merged_row]==actual_row || throw(ArgumentError(
                "$caller: merged geometry vertex has conflicting actual identities"))
            rows[merged_row]=actual_row
        end
    end
    all(>(0),rows) && allunique(rows) || throw(ArgumentError(
        "$caller: produced geometry rows do not bind uniquely to actual vertices"))
    for ((_,part,_,mapping),merged_rows) in zip(parts,remaps)
        for (position,merged_row) in enumerate(merged_rows)
            source=assembly.nodes[mapping[position]]
            expected=geometry_birth===nothing ? source :
                get(geometry_birth.nodes,source.tag,source)
            expected.owner==source.owner || (expected=source)
            isequal(carrier.coords[:,merged_row],collect(expected.point)) ||
                throw(ArgumentError("$caller: produced geometry lacks its explicit birth value"))
        end
    end
    emitted=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
        (class.entity[1],tag,mesh) for (tag,mesh,class,_) in parts]
    emitted_rows=Vector{Int32}[collect(Int32,mapping) for mapping in remaps]
    return carrier,rows,carrier_cells,emitted,emitted_rows
end

function _higher_geometry_keep_points(parts,assembly,node_tag_map,caller)
    kept_parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[]
    for (tag,mesh,class,mapping) in parts
        if class.entity[1]!=0
            push!(kept_parts,(tag,mesh,class,mapping))
            continue
        end
        mesh isa Mesh && nsegs(mesh)==0 && ntris(mesh)==0 && ntets(mesh)==0 ||
            throw(ArgumentError("$caller: Point geometry must be a vertex-only producer"))
        kept=Int[]
        for position in eachindex(mapping)
            source=assembly.nodes[mapping[position]].tag
            get(node_tag_map,source,UInt64(0))>0 && push!(kept,position)
        end
        isempty(kept) && continue
        if length(kept)==length(mapping)
            push!(kept_parts,(tag,mesh,class,mapping))
        else
            geometry=Mesh(mesh.coords[:,kept])
            owners=class.node_entities[kept]
            classification=_higher_geometry_part_class(geometry,0,tag,owners)
            push!(kept_parts,(tag,geometry,classification,mapping[kept]))
        end
    end
    return kept_parts
end

function _generate_higher_plan(source_model,cached,cached_class,dimension::Int,
        options::Model._ModelMesh1DOptions,caller;generation_meshing,
        element_order=1,renumber=true,save_all=false,linear=false,incomplete=false,size_field=nothing,
        metric_surface=nothing,
        initial_max_node_tag=UInt64(0),initial_max_element_tag=UInt64(0),
        preserve_geometry=false,geometry_layout=:top,geometry_merge=nothing)
    dimension in (2,3) || throw(ArgumentError("$caller: dim must be two or three"))
    geometry_layout in (:top,:all) || throw(ArgumentError("$caller: invalid geometry layout"))
    m=Model._model_mesh_transaction_copy(source_model)
    _apply_generation_meshing_options!(m,generation_meshing)
    # Each native Model dispatcher admits the ordinary algorithm after its
    # transfinite, extrusion and periodic construction paths. An early API
    # whitelist would reject those specialized paths before their dispatch.
    entities=model_entities(m)
    metric_surface===nothing || (dimension==2 && element_order==1 &&
        [(dim,tag) for (dim,tag) in entities if dim>=2]==[(2,metric_surface.tag)]) ||
        throw(ArgumentError("$caller: metric Surface request must select the sole two-dimensional entity"))
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    source_node_positions=_dim01_ingest!(assembly,m,cached,cached_class;preserve_source_order=true)
    old=m.meshing.mesh_authority
    _higher_native_parameters!(assembly,m)
    # Legacy simplex caches can carry zero proxy physical tags while their
    # current CAD entity belongs to a real group. Entity ownership controls
    # the actual projection, including groups edited since the cache was born.
    if cached!==nothing && (!(cached isa MixedMesh) || cached.entity_data===nothing)
        for cell in assembly.cells
            cell.physical=_dim01_physical(m,(cell.owner[1],Int(cell.owner[2])))
        end
    end
    done=old===nothing ? Set{Tuple{Int,Int32}}() : copy(old.done)
    ranks=_higher_node_ranks(assembly,cached,source_node_positions,old)
    actual=_higher_assembly_authority(assembly,m,done;ranks)
    status=Model._model_mesh_status(m,actual;only_visible=options.mesh_only_visible)
    phases=Model._model_mesh_generation_phases(dimension,status)
    _higher_ingest_native_inputs!(assembly,m,old)
    _higher_native_parameters!(assembly,m)
    _higher_complete_ranks!(ranks,assembly)
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    high_order=if metric_surface===nothing
        _dim01_first_order!(assembly,m)
    else
        # This explicit tensor endpoint retains actual lower cells. Native
        # generation's SetOrder1 factory renews even linear cell identities;
        # admit existing linear input instead of running that reconstruction.
        all(cell->cell.owner[1]==0 || msh_spec(cell.msh).order==1,assembly.cells) ||
            throw(ArgumentError("$caller: metric Surface input must contain only first-order cells"))
        falses(length(assembly.nodes))
    end
    first_order_constructor_inputs=Dict{Tuple{Int,Int32},MixedMesh}()
    if old!==nothing && !phases.one && !isempty(old.surface_constructor_inputs)
        expected=Set(uid for constructor in values(old.surface_constructor_inputs)
            for labels in constructor.entity_data.external_element_tags for uid in labels)
        selected=[cell for cell in assembly.cells if cell.source_tag in expected]
        length(selected)==length(expected) && allunique(cell.source_tag for cell in selected) ||
            throw(ArgumentError("$caller: SetOrder1 loses a retained Surface constructor cell event"))
        mapping=Dict{UInt64,UInt64}(cell.source_tag=>cell.tag for cell in selected)
        current=_higher_assembly_authority(assembly,m,done;ranks,old_supports=high_order)
        first_order_constructor_inputs=Model._model_mesh_first_order_surface_constructor_inputs(
            old,current,mapping;caller)
    end
    surface_parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[]
    geometry_cells=preserve_geometry ?
        IdDict{Any,Vector{Vector{Union{Nothing,_Dim01Cell}}}}() : nothing
    geometry_birth=preserve_geometry ? _dim01_geometry_birth(old) : nothing
    curve_births=Dict{Tuple{Int,Int32},Vector{UInt64}}()
    geometry_alternates=preserve_geometry ? IdDict{Any,Any}() : nothing
    geometry_parts=!preserve_geometry ? nothing : dimension==2 ? surface_parts :
        Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[]
    discarded=Set{Tuple{Int,Int32}}()
    native_orphan_births=Set{UInt64}()
    if phases.one
        for entity in entities
            entity[1]>=2 && !_dim01_fully_discrete(m,entity) &&
                push!(discarded,(entity[1],Int32(entity[2])))
        end
        filter!(cell->!(cell.owner in discarded),assembly.cells)
        filter!(owner->!(owner in discarded),done)
        first_nodes=_dim01_point_phase!(assembly,m,entities)
        union!(discarded,_dim01_curve_phase!(assembly,m,options,entities,first_nodes,high_order;
            done,geometry_birth,curve_births))
    else
        # The legacy Surface/Volume projection omitted its completed Mesh0D
        # products. Recover those once, just as skipped Curve products are
        # recovered below; an actual cache never invents a skipped Point cell.
        legacy_points=old!==nothing && old.mesh===nothing && cached!==nothing
        first_nodes=_dim01_point_phase!(assembly,m,entities;
                                       create_cells=legacy_points,create_nodes=legacy_points)
        _higher_retained_curves!(assembly,m,entities,first_nodes,high_order,done)
    end
    periodic_pairs=Dict{Int,Tuple{Int,Vector{Tuple{UInt64,UInt64}}}}()
    native_recombine_events=Tuple{Symbol,String}[]
    surface_constructor_inputs=first_order_constructor_inputs
    retained_surfaces=Set(keys(first_order_constructor_inputs))
    if phases.two
        for entity in entities
            entity[1]==3 && !_dim01_fully_discrete(m,entity) &&
                push!(discarded,(3,Int32(entity[2])))
        end
        filter!(cell->cell.source_tag===nothing || !(cell.owner in discarded),assembly.cells)
        _,_,retained_surfaces=_higher_surface_phase!(assembly,m,options,entities,discarded,done;
                               size_field,metric_surface,source_nodes=length(high_order),ranks,old_supports=high_order,
                                geometry_parts=surface_parts,dimension,geometry_cells,
                                geometry_birth,geometry_alternates,curve_births,periodic_pairs,native_orphan_births,native_recombine_events,
                                surface_constructor_inputs)
    end
    if phases.three
        _higher_volume_phase!(assembly,m,options,entities,discarded,done;
                             size_field,source_nodes=length(high_order),ranks,old_supports=high_order,
                             geometry_parts,surface_parts,geometry_cells,geometry_birth,native_orphan_births,
                             periodic_surface_pairs=periodic_pairs,
                             optimize_pyramids=Float64(get(generation_meshing,:optimize_pyramids,0.0)),
                             retained_surfaces,
                             surface_constructor_inputs)
    end
    if preserve_geometry && geometry_layout===:all
        # GEO's native carrier has the Surface/Volume products followed by
        # Point/Curve geometry. Reuse each actual emitted Surface vector and
        # the staged lower cells' explicit assembly rows; this does not mesh
        # another product or infer a public identity from coordinates.
        _higher_complete_ranks!(ranks,assembly)
        surfaces=Dict(tag=>(tag,mesh,class,mapping) for (tag,mesh,class,mapping) in surface_parts)
        retained=_higher_assembly_authority(assembly,m,done;ranks,old_supports=high_order)
        retained=_higher_bind_producer_rows(retained,assembly,surface_parts,old)
        retained_cells=Dict(cell.tag=>cell for cell in assembly.cells)
        for owner in keys(retained.products)
            owner[1]==2 && !haskey(surfaces,Int(owner[2])) || continue
            input=_higher_surface_geometry_input(retained,owner,caller)
            input===nothing && continue
            geometry=input.mesh
            mapping=Int32[assembly.node_positions[tag] for tag in input.node_tags]
            class=_higher_geometry_part_class(geometry,2,Int(owner[2]),input.node_entities)
            surfaces[Int(owner[2])]=(Int(owner[2]),geometry,class,mapping)
            geometry_cells[geometry]=[Union{Nothing,_Dim01Cell}[],
                Union{Nothing,_Dim01Cell}[retained_cells[tag] for tags in
                    input.element_tags for tag in tags],
                Union{Nothing,_Dim01Cell}[]]
        end
        volumes=dimension==3 ? geometry_parts :
            Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[]
        geometry_parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[
            surfaces[tag] for tag in sort!(collect(keys(surfaces)))]
        append!(geometry_parts,volumes)
        sort_points=any(volume->Model._extrude_is_nonew(Model._extrude_gate(m,3,volume)),
                        keys(m.volumes))
        append!(geometry_parts,_higher_lower_geometry_parts(assembly,m;
            sort_points,geometry_cells,geometry_birth))
    end
    p2_pass=element_order>1 && any(cell->cell.owner[1]>0,assembly.cells)
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    if element_order>2 || (p2_pass && incomplete)
        _dim01_higher_order!(assembly,m,entities,element_order;
            only_visible=options.mesh_only_visible,linear,incomplete)
    elseif p2_pass
        _dim01_quadratic!(assembly,m,entities;only_visible=options.mesh_only_visible,linear)
    end
    _higher_complete_ranks!(ranks,assembly)
    constructor_cells=Dict{UInt64,_Dim01Cell}()
    for cell in assembly.cells
        cell.owner[1]==2 || continue
        haskey(constructor_cells,cell.tag) && throw(ArgumentError(
            "$caller: Surface constructor publication repeats an actual cell UID"))
        constructor_cells[cell.tag]=cell
    end
    plan=_higher_periodic_surface_links(
        _dim01_finish!(assembly,m,entities,discarded,high_order,renumber,p2_pass;
            save_all=save_all,source_mesh=_higher_periodic_source_mesh(cached,periodic_pairs),
            source_node_positions=source_node_positions,curve_births,native_orphan_births,
            periodic_uid_operation=:defer),
        periodic_pairs,m,caller)
    _periodic_uid_history_generate!(plan,periodic_pairs,caller)
    public_ranks=Dict(plan.node_tag_map[tag]=>rank for (tag,rank) in ranks
                     if haskey(plan.node_tag_map,tag))
    final_ranks=UInt64[public_ranks[tag] for tag in plan.mesh.entity_data.external_node_tags]
    if geometry_parts!==nothing && geometry_layout===:all
        geometry_parts=_higher_geometry_keep_points(geometry_parts,assembly,plan.node_tag_map,caller)
    end
    if geometry_parts!==nothing
        geometry_parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[
            (tag,_higher_geometry_emitter_part(get(geometry_alternates,mesh,mesh),
                mapping,assembly,geometry_birth,geometry_cells),class,mapping)
            for (tag,mesh,class,mapping) in geometry_parts]
    end
    carrier,carrier_rows,carrier_cells,emitted_parts,emitted_rows=geometry_parts===nothing ?
        (nothing,Int32[],Vector{UInt64}[],nothing,nothing) :
        _higher_geometry_carrier(plan,assembly,geometry_parts,dimension,caller;
            geometry_merge,geometry_cells,geometry_birth)
    m.meshing.mesh_authority=carrier===nothing ?
        Model._model_mesh_authority(plan.mesh,done;caller,private_rank=final_ranks) :
        Model._model_mesh_authority(plan.mesh,done;caller,private_rank=final_ranks,
            carrier,carrier_rows,carrier_cells)
    m.meshing.mesh_authority=_higher_bind_producer_rows(m.meshing.mesh_authority,assembly,
         surface_parts,old;node_tag_map=plan.node_tag_map,curve_births)
    m.meshing.mesh_authority=Model._model_mesh_with_curve_births(
        m.meshing.mesh_authority,curve_births;node_tag_map=plan.node_tag_map,caller)
    if element_order==1 && !isempty(surface_constructor_inputs)
        m.meshing.mesh_authority=_higher_surface_constructor_publication(
            m.meshing.mesh_authority,surface_constructor_inputs,constructor_cells,
            plan.node_tag_map,caller)
    end
    if emitted_parts!==nothing && any(emitted_parts) do part
            part[1]==3 || return false
            Model._model_mesh_emitted_scope(m,part[2])
        end
        m.meshing.mesh_authority=Model._model_mesh_with_emitted_receipt(
            m.meshing.mesh_authority,emitted_parts,emitted_rows;caller)
    end
    return merge(plan,(geometry_parts=emitted_parts,geometry_part_rows=emitted_rows,
        native_recombine_events))
end

function _higher_classification_from_data(mesh,caller)
    data=mesh.entity_data
    data===nothing && throw(ArgumentError("$caller: actual products need entity metadata"))
    declarations=mesh.partition_data===nothing ? data.entities :
        merge(data.entities,mesh.partition_data.entities)
    entities=sort!(Tuple{Int,Int32}[(dim,Int32(tag)) for (dim,tag) in keys(declarations)])
    boundaries=Dict((dim,Int32(tag))=>Int32[abs(value) for value in entity.boundaries]
                    for ((dim,tag),entity) in declarations)
    owners=Dict{Int,Vector{Int32}}()
    for (block,values) in zip(mesh.blocks,data.block_entities)
        append!(get!(()->Int32[],owners,block.msh),values)
    end
    authority=union(Set(data.node_entities),Set((msh_dimension(block.msh),owner)
        for (block,values) in zip(mesh.blocks,data.block_entities) for owner in values))
    primary=isempty(entities) ? (0,Int32(0)) : last(entities)
    class=_mixed_classification(mesh,primary,entities,data.node_entities,boundaries,owners)
    return _classification_with_public_tags(class,_cache_public_tags(mesh;authority,caller))
end

function _higher_authority_from_parts(m,merged,parts,remaps,node_owner,caller;
                                     node_max=UInt64(0),element_max=UInt64(0),
                                     producer_order=nothing,producer_parameters=nothing,
                                     producer_coordinates=nothing)
    length(parts)==length(remaps) || throw(ArgumentError("$caller: actual part rows do not align"))
    blocks=ElementBlock[];block_owners=Vector{Int32}[];owners=Dict{Int,Vector{Int32}}()
    entities=Tuple{Int,Int32}[]
    point_rows=Int32[];point_owners=Int32[];seen_points=Set{Tuple{Int,Int32}}()
    for ((dim,tag,part),rows) in zip(parts,remaps)
        length(rows)==nnodes(part) && all(row->0<row<=nnodes(merged),rows) ||
            throw(ArgumentError("$caller: actual entity ($dim,$tag) has invalid merge rows"))
        push!(entities,(dim,Int32(tag)))
        actual_point_cells=any(msh==15 && !isempty(cells) for (msh,cells,_) in _cache_native_blocks(part))
        point_record=dim==0 ? _dim01_record(m,(0,tag)) : nothing
        recorded_point_cells=point_record!==nothing && any(==(15),point_record.element_types)
        if dim==0 && !isempty(rows) && !actual_point_cells && !recorded_point_cells
            key=(tag,rows[end])
            if !(key in seen_points)
                push!(seen_points,key);push!(point_rows,rows[end]);push!(point_owners,Int32(tag))
            end
        end
        for (msh,cells,physical) in _cache_native_blocks(part)
            msh_dimension(msh)==dim || continue
            # Legacy geometry parts can carry zero proxy tags even when the
            # current entity belongs to a Physical group. Bind its actual
            # projection before attaching strict entity metadata.
            projected=fill(_dim01_physical(m,(dim,tag)),size(cells,2))
            push!(blocks,ElementBlock(msh,rows[cells],projected))
            classified=fill(Int32(tag),size(cells,2));push!(block_owners,classified)
            append!(get!(()->Int32[],owners,msh),classified)
        end
    end
    if !isempty(point_rows)
        physical=Int32[_dim01_physical(m,(0,Int(tag))) for tag in point_owners]
        push!(blocks,ElementBlock(15,reshape(point_rows,1,:),physical))
        push!(block_owners,point_owners);append!(get!(()->Int32[],owners,15),point_owners)
    end
    coordinates=producer_coordinates===nothing ? merged.coords : producer_coordinates
    size(coordinates)==size(merged.coords) && all(isfinite,coordinates) ||
        throw(ArgumentError("$caller: native producer coordinates do not align with bound rows"))
    cached=MixedMesh(coordinates,blocks;elementary_entities=block_owners)
    boundaries=Dict(entity=>Int32[abs(tag) for (_,tag) in
        model_boundary(m,[(entity[1],Int(entity[2]))],false,false,false)] for entity in entities)
    node_entities=Tuple{Int,Int32}[(node_owner[row][1],Int32(node_owner[row][2]))
                                for row in 1:nnodes(merged)]
    primary=isempty(entities) ? (0,Int32(0)) : last(entities)
    class=_mixed_classification(cached,primary,entities,node_entities,boundaries,owners)
    if producer_order!==nothing
        # This is the completed native producer schedule supplied by GEO.
        # Merge storage order remains independent of these factory identities.
        for (_,record) in Model._discrete_mesh_records_model(m)
            node_max=max(node_max,maximum(record.node_tags;init=UInt64(0)))
            element_max=max(element_max,maximum(record.element_tags;init=UInt64(0)))
        end
        node_labels=zeros(UInt64,nnodes(merged))
        cell_labels=[zeros(UInt64,size(block.nodes,2)) for block in blocks]
        for (dimension,tag) in producer_order
            owner=(dimension,Int32(tag))
            for ((part_dimension,part_tag,_),rows) in zip(parts,remaps)
                part_dimension==dimension && part_tag==tag || continue
                for row in rows
                    node_entities[row]==owner && node_labels[row]==0 || continue
                    node_max<typemax(UInt64) || throw(ArgumentError("$caller: produced node tags exceed UInt64"))
                    node_max+=UInt64(1);node_labels[row]=node_max
                end
            end
            for (bi,block) in enumerate(blocks)
                msh_dimension(block.msh)==dimension || continue
                for column in axes(block.nodes,2)
                    block_owners[bi][column]==tag && cell_labels[bi][column]==0 || continue
                    element_max<typemax(UInt64) || throw(ArgumentError("$caller: produced element tags exceed UInt64"))
                    element_max+=UInt64(1);cell_labels[bi][column]=element_max
                end
            end
        end
        all(>(0),node_labels) && all(labels->all(>(0),labels),cell_labels) ||
            throw(ArgumentError("$caller: producer schedule omits actual node or cell rows"))
        staged=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
            UInt64(node_max),UInt64(element_max),Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
        # Model-derived carrier metadata uses explicit merge rows. Public
        # identities enter before ingest can encounter the legacy dense slots.
        _dim01_ingest!(staged,m,nothing,nothing)
        metadata,_=_dim01_entities(m,model_entities(m),staged,eachindex(staged.nodes))
        parameters=producer_parameters===nothing ? fill(nothing,nnodes(merged)) :
            producer_parameters
        length(parameters)==nnodes(merged) || throw(ArgumentError(
            "$caller: produced native parameters do not align with geometry rows"))
        data=Elements._native_mixed_entity_data(metadata;node_entities,
            node_parametric=parameters,external_node_tags=node_labels,
            block_entities=block_owners,external_element_tags=cell_labels)
        cached=MixedMesh(coordinates,blocks;entity_data=data,elementary_entities=block_owners)
        class=_higher_classification_from_data(cached,caller)
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(node_max),UInt64(element_max),Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    mapping=_dim01_ingest!(assembly,m,cached,class)
    for cell in assembly.cells
        cell.physical=_dim01_physical(m,(cell.owner[1],Int(cell.owner[2])))
    end
    _higher_native_parameters!(assembly,m)
    ranks=_higher_node_ranks(assembly,cached,mapping,nothing)
    return _higher_assembly_authority(assembly,m,Set(entities);ranks,
        carrier=merged,carrier_rows=mapping)
end

function _higher_frontend_input(m,cached,parts,remaps,node_owner,caller;
                                initial_max_node_tag=UInt64(0),
                                initial_max_element_tag=UInt64(0))
    authority=if cached===nothing
        nothing
    elseif m.meshing.mesh_authority!==nothing &&
            m.meshing.mesh_authority.mesh!==nothing &&
            m.meshing.mesh_authority.carrier===cached &&
            m.meshing.mesh_authority.mesh!==cached
        bound_previous=m.meshing.mesh_authority
        length(bound_previous.carrier_rows)==nnodes(cached) &&
            all(row->0<row<=nnodes(bound_previous.mesh),bound_previous.carrier_rows) &&
            allunique(bound_previous.carrier_rows) || throw(ArgumentError(
                "$caller: geometry carrier loses its bound actual vertex identities"))
        bound_previous
    elseif cached isa MixedMesh && cached.entity_data!==nothing
        cached_previous=m.meshing.mesh_authority
        done=cached_previous===nothing ? Set{Tuple{Int,Int32}}() : cached_previous.done
        ranks=cached_previous!==nothing && cached_previous.mesh!==nothing &&
              length(cached_previous.private_rank)==nnodes(cached) ? cached_previous.private_rank :
              collect(UInt64,1:nnodes(cached))
        cached_previous!==nothing && cached_previous.mesh!==nothing ?
            Model._model_mesh_rebind_authority(cached_previous,cached) :
            Model._model_mesh_authority(cached,done;caller,private_rank=ranks)
    else
        _higher_authority_from_parts(m,cached,parts,remaps,node_owner,caller;
            node_max=initial_max_node_tag,element_max=initial_max_element_tag)
    end
    source=Model._model_mesh_transaction_copy(m)
    source.meshing.mesh_authority=authority
    actual=authority===nothing ? nothing : authority.mesh
    class=actual===nothing ? nothing : _higher_classification_from_data(actual,caller)
    return source,actual,class
end

function _higher_rebind_plan_authority!(plan,previous;
        done=previous===nothing ? Set{Tuple{Int,Int32}}() : previous.done,
        native_inputs=previous===nothing ?
            Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}() : previous.native_inputs,
        carrier=plan.mesh,carrier_rows=collect(Int32,1:nnodes(plan.mesh)),carrier_cells=nothing)
    old_ranks=Dict{UInt64,UInt64}()
    next=UInt64(0)
    if previous!==nothing && previous.mesh!==nothing
        next=maximum(previous.private_rank;init=UInt64(0))
        for (tag,rank) in zip(previous.mesh.entity_data.external_node_tags,previous.private_rank)
            mapped=get(plan.node_tag_map,tag,UInt64(0))
            mapped==0 || (old_ranks[mapped]=rank)
        end
    end
    ranks=Vector{UInt64}(undef,nnodes(plan.mesh))
    for (row,tag) in enumerate(plan.mesh.entity_data.external_node_tags)
        rank=get(old_ranks,tag,UInt64(0))
        if rank==0
            next<typemax(UInt64) || throw(ArgumentError("mesh authority: private vertex rank overflows UInt64"))
            next+=UInt64(1);rank=next
        end
        ranks[row]=rank
    end
    if carrier_cells===nothing && previous!==nothing && carrier===previous.carrier
        carrier_cells=Model._model_mesh_rebound_cell_tags(previous,
            hasproperty(plan,:element_tag_map) ? plan.element_tag_map : nothing)
    end
    authority=Model._model_mesh_authority(plan.mesh,done;
        private_rank=ranks,carrier,carrier_rows,carrier_cells)
    produced=Dict{Tuple{Int,Int32},Vector{Int32}}()
    if previous!==nothing && previous.mesh!==nothing
        positions=Dict(tag=>Int32(row) for (row,tag) in
                       enumerate(plan.mesh.entity_data.external_node_tags))
        for (owner,old_rows) in previous.producer_rows
            haskey(authority.products,owner) || continue
            rows=Int32[get(positions,get(plan.node_tag_map,
                previous.mesh.entity_data.external_node_tags[row],UInt64(0)),Int32(0))
                for row in old_rows]
            all(>(0),rows) || continue # A changed product invalidates its old vector.
            produced[owner]=rows
        end
    end
    plan.model.meshing.mesh_authority=Model._ModelMeshAuthority(authority.mesh,
        authority.products,authority.done,authority.private_rank,copy(native_inputs),
        authority.carrier,authority.carrier_rows,produced,authority.carrier_cells)
    return plan
end

function _higher_capture_dim01_authority!(plan)
    previous=plan.model.meshing.mesh_authority
    # A dimension-one generation deMeshes native higher entities. A zero
    # generation performs the order passes and preserves their actual inputs.
    inputs=plan.mesh_dimension==1 || previous===nothing ?
        Dict{Tuple{Int,Int32},Union{Mesh,MixedMesh}}() : previous.native_inputs
    _higher_rebind_plan_authority!(plan,previous;
        done=plan.mesh_done,native_inputs=inputs,
        carrier=plan.geometry_carrier,
        carrier_cells=plan.mesh.entity_data.external_element_tags)
    plan.model.meshing.mesh_authority=Model._model_mesh_with_curve_births(
        plan.model.meshing.mesh_authority,plan.curve_births;
        node_tag_map=plan.node_tag_map,caller="dimension-one Curve publication")
    return plan
end

function _higher_prepare_api_plan_authority!(plan)
    plan.mesh isa MixedMesh && plan.mesh.entity_data!==nothing || return plan
    previous=plan.model.meshing.mesh_authority
    previous!==nothing && previous.mesh!==nothing && previous.mesh!==plan.mesh || return plan
    hasproperty(plan,:node_tag_map) || return plan
    birth_carrier=hasproperty(plan,:geometry_carrier) &&
        !(previous.carrier isa MixedMesh && previous.carrier.entity_data===nothing)
    carrier,rows=if birth_carrier
        plan.geometry_carrier,collect(Int32,1:nnodes(plan.mesh))
    elseif previous.carrier isa Mesh && hasproperty(plan,:refine_primary_rows)
        tags=previous.mesh.entity_data.external_node_tags[previous.carrier_rows]
        Model._model_mesh_refined_carrier(plan,previous.carrier,tags)
    else
        Model._model_mesh_rebound_carrier(previous,plan.mesh;node_tag_map=plan.node_tag_map,
            element_tag_map=hasproperty(plan,:element_tag_map) ? plan.element_tag_map : nothing)
    end
    cells=birth_carrier ? plan.mesh.entity_data.external_element_tags :
        hasproperty(plan,:refine_cell_tags) ?
        Model._model_mesh_refined_cell_tags(plan,previous) :
        Model._model_mesh_rebound_cell_tags(previous,
            hasproperty(plan,:element_tag_map) ? plan.element_tag_map : nothing)
    return carrier===nothing ? _higher_rebind_plan_authority!(plan,previous) :
        _higher_rebind_plan_authority!(plan,previous;carrier,carrier_rows=rows,carrier_cells=cells)
end

function _higher_recombine_protected(mesh,m,owner)
    curves=Set{Int32}()
    for (dimension,tag) in model_boundary(m,[(2,Int(owner))],false,false,false)
        dimension==1 && push!(curves,Int32(abs(tag)))
    end
    for (dimension,tag) in get(m.embeds,(2,Int(owner)),Tuple{Int,Int}[])
        dimension==1 && push!(curves,Int32(abs(tag)))
    end
    protected=Set{NTuple{2,Int32}}()
    for (bi,block) in enumerate(mesh.blocks)
        msh_spec(block.msh).family===:lin || continue
        for column in axes(block.nodes,2)
            mesh.entity_data.block_entities[bi][column] in curves || continue
            push!(protected,minmax(block.nodes[1,column],block.nodes[2,column]))
        end
    end
    return protected
end

function _higher_recombine_plan(m,actual,caller;algorithm=:greedy,
        initial_max_node_tag=UInt64(0),initial_max_element_tag=UInt64(0),
        max_nodes=typemax(Int32),max_cells=typemax(Int32))
    data=actual.entity_data
    nnodes(actual)<=max_nodes || throw(ArgumentError("$caller: mesh exceeds max_nodes=$max_nodes"))
    maximum_element=max(UInt64(initial_max_element_tag),
        maximum((maximum(tags;init=UInt64(0)) for tags in data.external_element_tags);init=UInt64(0)))
    maximum_node=max(UInt64(initial_max_node_tag),maximum(data.external_node_tags;init=UInt64(0)))
    for (_,record) in Model._discrete_mesh_records_model(m)
        maximum_element=max(maximum_element,maximum(record.element_tags;init=UInt64(0)))
        maximum_node=max(maximum_node,maximum(record.node_tags;init=UInt64(0)))
    end
    groups=Dict{Int32,Vector{Tuple{Int,Int}}}()
    for (bi,block) in enumerate(actual.blocks)
        block.msh==2 || continue
        for column in axes(block.nodes,2)
            owner=data.block_entities[bi][column]
            push!(get!(()->Tuple{Int,Int}[],groups,owner),(bi,column))
        end
    end
    replacements=Dict{Int32,Tuple{Vector{ElementBlock},Vector{Vector{UInt64}}}}()
    changed=false
    for owner in sort!(collect(keys(groups)))
        selected=groups[owner]
        length(selected)>=2 || continue
        triangles=Matrix{Int32}(undef,3,length(selected));physical=Int32[]
        previous=Dict{NTuple{3,Int32},Vector{UInt64}}()
        incidence=Dict{NTuple{2,Int32},Int}()
        for (slot,(bi,column)) in enumerate(selected)
            block=actual.blocks[bi];nodes=block.nodes[:,column]
            triangles[:,slot]=nodes;push!(physical,block.tags[column])
            push!(get!(()->UInt64[],previous,Tuple(sort(nodes))),data.external_element_tags[bi][column])
            for (first,second) in ((1,2),(2,3),(3,1))
                edge=minmax(nodes[first],nodes[second]);incidence[edge]=get(incidence,edge,0)+1
            end
        end
        protected=_higher_recombine_protected(actual,m,owner)
        candidates=0
        for (edge,uses) in incidence
            uses==2 && !(edge in protected) && (candidates+=1)
        end
        angle=get(m.meshing.recombine,(2,Int(owner)),45.0)
        result=recombine_triangles(Mesh(actual.coords;tris=triangles,tri_tag=physical);
            algorithm,protected_edges=protected,recombine_angle=angle)
        result.coords==actual.coords || throw(ArgumentError(
            "$caller: recombination with newly created nodes is not implemented"))
        quadrangles=sum(size(block.nodes,2) for block in result.blocks if block.msh==3;init=0)
        # The pinned first pair pass reserves temporary candidate elements
        # before strict angle screening, then allocates successful quadrangles.
        reservation=UInt64(candidates)+UInt64(quadrangles)
        limit=typemax(UInt64)
        (maximum_element>limit || reservation>limit-maximum_element) &&
            throw(ArgumentError("$caller: recombined element tags exceed UInt64"))
        maximum_element+=UInt64(candidates)
        quadrangles==0 && continue
        changed=true;blocks=ElementBlock[];tags=Vector{UInt64}[]
        for block in result.blocks
            block.msh in (2,3) || continue
            isempty(block.nodes) && continue
            external=UInt64[]
            for nodes in eachcol(block.nodes)
                if block.msh==2
                    source=get(previous,Tuple(sort(nodes)),nothing)
                    source!==nothing && !isempty(source) || throw(ArgumentError(
                        "$caller: a retained triangle lost its source identity"))
                    push!(external,popfirst!(source))
                else
                    # Gmsh's greedy MQuadrangle constructor starts at the
                    # first vertex of the joined MEdge. Its pointer rank is
                    # independent of a public rename and of kernel rotation.
                    first_edge=minmax(nodes[1],nodes[3])
                    second_edge=minmax(nodes[2],nodes[4])
                    joined=get(incidence,first_edge,0)==2 ? first_edge : second_edge
                    get(incidence,joined,0)==2 || throw(ArgumentError(
                        "$caller: quadrangle lost its joined source edge"))
                    ranks=m.meshing.mesh_authority.private_rank
                    first_vertex=ranks[joined[1]]<ranks[joined[2]] ? joined[1] : joined[2]
                    start=0
                    for row in 1:4
                        nodes[row]==first_vertex && (start=row)
                    end
                    original=copy(nodes)
                    for row in 1:4
                        nodes[row]=original[mod1(start+row-1,4)]
                    end
                    maximum_element+=UInt64(1);push!(external,maximum_element)
                end
            end
            push!(blocks,block);push!(tags,external)
        end
        replacements[owner]=(blocks,tags)
    end
    mesh=actual
    if changed
        (any(section->!isempty(section.elements),actual.data_sections) ||
            (actual.partition_data!==nothing && !isempty(actual.partition_data.ghost_elements))) &&
            throw(ArgumentError("$caller: cannot preserve cell-attached MSH metadata while changing triangle blocks"))
        blocks=ElementBlock[];owners=Vector{Int32}[];tags=Vector{UInt64}[]
        emitted=Set{Int32}()
        for (bi,block) in enumerate(actual.blocks)
            if block.msh!=2 || isempty(block.nodes)
                push!(blocks,block);push!(owners,data.block_entities[bi]);push!(tags,data.external_element_tags[bi])
                continue
            end
            kept=Int[column for column in axes(block.nodes,2)
                     if !haskey(replacements,data.block_entities[bi][column])]
            if !isempty(kept)
                push!(blocks,ElementBlock(block.msh,block.nodes[:,kept],block.tags[kept]))
                push!(owners,data.block_entities[bi][kept]);push!(tags,data.external_element_tags[bi][kept])
            end
            for owner in data.block_entities[bi]
                haskey(replacements,owner) && !(owner in emitted) || continue
                push!(emitted,owner);output,external=replacements[owner]
                for (replacement,labels) in zip(output,external)
                    push!(blocks,replacement);push!(owners,fill(owner,size(replacement.nodes,2)));push!(tags,labels)
                end
            end
        end
        sum(size(block.nodes,2) for block in blocks;init=0)<=max_cells ||
            throw(ArgumentError("$caller: mesh exceeds max_cells=$max_cells"))
        metadata=Elements._native_mixed_entity_data(data.entities;node_entities=data.node_entities,
            node_parametric=data.node_parametric,external_node_tags=data.external_node_tags,
            block_entities=owners,external_element_tags=tags)
        mesh=MixedMesh(actual.coords,blocks;physical_names=actual.physical_names,
            entity_data=metadata,elementary_entities=owners,periodic_links=actual.periodic_links,
            ancillary_sections=actual.ancillary_sections,data_sections=actual.data_sections,
            partition_data=actual.partition_data)
    end
    class=_higher_classification_from_data(mesh,caller)
    nodes=Dict(tag=>tag for tag in data.external_node_tags)
    retained=Set(Iterators.flatten(mesh.entity_data.external_element_tags))
    elements=Dict(tag=>tag for tag in Iterators.flatten(data.external_element_tags) if tag in retained)
    covered=Tuple{Int,Int}[(dimension,Int(tag)) for (dimension,tag) in class.public_tags.authority]
    _dim01_reconcile_records!(m,covered,mesh,nodes)
    plan=(model=m,mesh=mesh,class=class,authority=class.public_tags.authority,
        node_tag_map=nodes,element_tag_map=elements,visibility_tag_map=elements,
        max_node_tag=maximum_node,max_element_tag=maximum_element,recombined=changed)
    previous=m.meshing.mesh_authority
    return changed || previous===nothing ? _higher_rebind_plan_authority!(plan,previous) :
        _higher_rebind_plan_authority!(plan,previous;carrier=previous.carrier,
            carrier_rows=previous.carrier_rows)
end

# A completed first-order GEO producer retains its geometry carrier. Apply
# the shared native save/renumber pass to its explicit factory labels, then
# rebind carrier rows and emitted Line cells through those identities.
function _higher_bind_completed_dim01!(m,authority,parts,caller;
        initial_max_node_tag,initial_max_element_tag,renumber,save_all,
        producer_remaps=nothing,producer_geometry_parameters=nothing)
    mesh=authority.mesh
    data=mesh.entity_data
    m.meshing.mesh_authority=authority
    class=_higher_classification_from_data(mesh,caller)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        max(UInt64(initial_max_node_tag),maximum(data.external_node_tags;init=UInt64(0))),
        max(UInt64(initial_max_element_tag),maximum((maximum(tags;init=UInt64(0))
            for tags in data.external_element_tags);init=UInt64(0))),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    mapping=_dim01_ingest!(assembly,m,mesh,class)
    source_cells=Dict{Int32,Vector{UInt64}}()
    for (block,owners,tags) in zip(mesh.blocks,data.block_entities,data.external_element_tags)
        block.msh==1 || continue
        for column in axes(block.nodes,2)
            push!(get!(()->UInt64[],source_cells,owners[column]),tags[column])
        end
    end
    carrier_cells=[UInt64[],UInt64[],UInt64[]]
    for (dimension,tag,part) in parts
        dimension==1 || continue
        tags=get(source_cells,Int32(tag),UInt64[])
        length(tags)==nsegs(part) || throw(ArgumentError(
            "$caller: native Curve[$tag] cell births do not align with its carrier"))
        append!(carrier_cells[1],tags)
    end
    finished=_dim01_finish!(assembly,m,model_entities(m),Set{Tuple{Int,Int32}}(),
        falses(length(assembly.nodes)),renumber,false;save_all,
        source_mesh=mesh,source_node_positions=mapping)
    positions=Dict(tag=>Int32(row) for (row,tag) in
        enumerate(finished.mesh.entity_data.external_node_tags))
    source_tags=data.external_node_tags[authority.carrier_rows]
    rows=Int32[get(positions,get(finished.node_tag_map,tag,nothing),Int32(0))
        for tag in source_tags]
    all(>(0),rows) && allunique(rows) || throw(ArgumentError(
        "$caller: completed native producer loses a carrier vertex identity"))
    finished.mesh.coords[:,rows]==authority.mesh.coords[:,authority.carrier_rows] ||
        throw(ArgumentError("$caller: completed native producer changes its actual birth coordinates"))
    cells=[Union{Nothing,UInt64}[get(finished.element_tag_map,tag,nothing) for tag in tags]
           for tags in carrier_cells]
    all(tags->all(!isnothing,tags),cells) || throw(ArgumentError(
        "$caller: completed native producer loses an emitted Line identity"))
    plan=_higher_rebind_plan_authority!(finished,authority;
        done=authority.done,carrier=authority.carrier,carrier_rows=rows,carrier_cells=cells)
    rebound=plan.model.meshing.mesh_authority
    if producer_remaps!==nothing
        length(parts)==length(producer_remaps) || throw(ArgumentError(
            "$caller: completed producer cell/row bindings do not align"))
        produced=copy(rebound.producer_rows)
        for ((dimension,tag,_),part_rows) in zip(parts,producer_remaps)
            dimension==1 || continue
            produced[(1,Int32(tag))]=Int32[rows[row] for row in part_rows]
        end
        plan.model.meshing.mesh_authority=Model._ModelMeshAuthority(rebound.mesh,
            rebound.products,rebound.done,rebound.private_rank,rebound.native_inputs,
            rebound.carrier,rebound.carrier_rows,produced,rebound.carrier_cells)
    end
    if producer_geometry_parameters!==nothing
        length(producer_geometry_parameters)==length(rows) || throw(ArgumentError(
            "$caller: alternate producer parameters do not align with bound rows"))
        parameters=Union{Nothing,Vector{Float64}}[nothing for _ in 1:nnodes(plan.mesh)]
        for (slot,row) in enumerate(rows)
            value=producer_geometry_parameters[slot]
            parameters[row]=value===nothing ? nothing : copy(value)
        end
        plan=merge(plan,(geometry_node_parametric=parameters,))
    end
    return plan
end

function _higher_geo_renumber_plan(m,mesh,class,operation,caller;
        initial_max_node_tag,initial_max_element_tag,save_all)
    node=operation===:geo_renumber_nodes
    if !node && _native_has_ghost_entities(m)
        plan=_tagged_renumber_plan(m,mesh,class,(),(),:element,caller;
            overlay=nothing,mid_owners=nothing,save_all,
            node_max=initial_max_node_tag,element_max=initial_max_element_tag)
        previous=m.meshing.mesh_authority
        carrier,rows=Model._model_mesh_rebound_carrier(previous,plan.mesh;
            node_tag_map=plan.node_tag_map,element_tag_map=plan.element_tag_map)
        carrier===nothing && return _higher_rebind_plan_authority!(plan,previous)
        cells=Model._model_mesh_rebound_cell_tags(previous,plan.element_tag_map)
        return _higher_rebind_plan_authority!(plan,previous;
            carrier,carrier_rows=rows,carrier_cells=cells)
    end
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        UInt64(initial_max_node_tag),UInt64(initial_max_element_tag),
        Int(typemax(Int32)),Int(typemax(Int32)),String(caller))
    _dim01_ingest!(assembly,m,mesh,class;preserve_source_order=true)
    physical=!save_all && any(!isempty,values(m.physical))
    positions=if node
        sortperm(assembly.nodes;by=value->value.owner,alg=Base.Sort.MergeSort)
    else
        sortperm(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    end
    selected=falses(node ? length(assembly.nodes) : length(assembly.cells))
    if physical
        for (position,cell) in enumerate(assembly.cells)
            _dim01_physical(m,(cell.owner[1],Int(cell.owner[2])))!=0 || continue
            if node
                selected[cell.nodes].=true
            else
                selected[position]=true
            end
        end
    end
    order=UInt64[];sizehint!(order,length(positions))
    for saved in (true,false),position in positions
        physical && selected[position]!=saved && continue
        !physical && !saved && continue
        push!(order,node ? assembly.nodes[position].tag : assembly.cells[position].tag)
    end
    kind=node ? :node : :element
    plan=_tagged_renumber_plan(m,mesh,class,order,collect(UInt64,1:length(order)),
        kind,caller;overlay=nothing,mid_owners=nothing,
        node_max=initial_max_node_tag,element_max=initial_max_element_tag)
    maximum_node=node && physical && !isempty(order) ?
        max(plan.max_node_tag,UInt64(length(order))+UInt64(1)) : plan.max_node_tag
    plan=merge(plan,(max_node_tag=maximum_node,))
    previous=m.meshing.mesh_authority
    carrier,rows=Model._model_mesh_rebound_carrier(previous,plan.mesh;
        node_tag_map=plan.node_tag_map,element_tag_map=plan.element_tag_map)
    carrier===nothing && return _higher_rebind_plan_authority!(plan,previous)
    cells=Model._model_mesh_rebound_cell_tags(previous,plan.element_tag_map)
    return _higher_rebind_plan_authority!(plan,previous;
        carrier,carrier_rows=rows,carrier_cells=cells)
end

function _generate_higher_frontend_plan(m,cached,parts,remaps,node_owner,dimension,options,caller;
        initial_max_node_tag=UInt64(0),initial_max_element_tag=UInt64(0),
        operation=:generate,element_order=1,renumber=true,save_all=false,linear=false,incomplete=false,
        max_nodes=typemax(Int32),max_cells=typemax(Int32),error_callback=nothing,kwargs...)
    if operation===:bind_completed
        source=Model._model_mesh_transaction_copy(m)
        authority=_higher_authority_from_parts(source,cached,parts,remaps,node_owner,caller;
            node_max=initial_max_node_tag,element_max=initial_max_element_tag,
            producer_order=get(kwargs,:producer_order,nothing),
            producer_parameters=get(kwargs,:producer_parameters,nothing),
            producer_coordinates=get(kwargs,:producer_coordinates,nothing))
        if dimension==1 && get(kwargs,:producer_order,nothing)!==nothing
            return _higher_bind_completed_dim01!(source,authority,parts,caller;
                initial_max_node_tag,initial_max_element_tag,renumber,save_all,
                producer_remaps=remaps,
                producer_geometry_parameters=get(kwargs,:producer_geometry_parameters,nothing))
        end
        source.meshing.mesh_authority=authority
        actual=authority.mesh;class=_higher_classification_from_data(actual,caller)
        maximum_node=max(UInt64(initial_max_node_tag),maximum(actual.entity_data.external_node_tags;init=UInt64(0)))
        maximum_element=max(UInt64(initial_max_element_tag),maximum(
            (maximum(tags;init=UInt64(0)) for tags in actual.entity_data.external_element_tags);init=UInt64(0)))
        nodes=Dict(tag=>tag for tag in actual.entity_data.external_node_tags)
        _dim01_reconcile_records!(source,model_entities(source),actual,nodes)
        return (model=source,mesh=actual,class=class,max_node_tag=maximum_node,
            max_element_tag=maximum_element)
    end
    source,actual,class=_higher_frontend_input(m,cached,parts,remaps,node_owner,caller;
        initial_max_node_tag,initial_max_element_tag)
    if operation===:generate
        if dimension in (0,1)
            plan=_generate_dim01_plan(source,actual,class,dimension,options,caller;
                initial_max_node_tag,initial_max_element_tag,element_order,renumber,save_all,linear,incomplete,
                initial_node_counter_wrapped=get(kwargs,:initial_node_counter_wrapped,false),
                initial_element_counter_wrapped=get(kwargs,:initial_element_counter_wrapped,false),
                max_nodes,max_cells,generation_meshing=get(kwargs,:generation_meshing,nothing))
            return _higher_capture_dim01_authority!(plan)
        end
        return _generate_higher_plan(source,actual,class,dimension,options,caller;
            initial_max_node_tag,initial_max_element_tag,element_order,renumber,save_all,linear,incomplete,kwargs...)
    end
    if operation===:native_bipartite_recombine
        native_options=get(kwargs,:native_recombine_options,nothing)
        native_options!==nothing || throw(ArgumentError(
            "$caller: experimental recombination lacks its detached options"))
        return _api_native_recombine_plan(source,actual,class,caller;options=native_options,
            initial_max_node_tag,initial_max_element_tag,max_nodes,max_cells,
            skip_recombination=get(kwargs,:native_recombine_skip,false))
    end
    actual===nothing && throw(ArgumentError("$caller: mesh operation requires actual products"))
    if operation in (:geo_renumber_nodes,:geo_renumber_elements)
        return _higher_geo_renumber_plan(source,actual,class,operation,caller;
            initial_max_node_tag,initial_max_element_tag,save_all)
    elseif operation===:recombine
        return _higher_recombine_plan(source,actual,caller;
            algorithm=get(kwargs,:algorithm,:greedy),initial_max_node_tag,
            initial_max_element_tag,max_nodes,max_cells)
    elseif operation===:homology
        cells,coordinates=Model._native_homology_actual_cells(actual)
        output=Model.compute_homology!(source,cells;coordinates,initial_max_element_tag,error_callback)
        if isempty(output)
            return (model=source,mesh=actual,class=class,
                max_node_tag=max(UInt64(initial_max_node_tag),maximum(
                    actual.entity_data.external_node_tags;init=UInt64(0))),
                max_element_tag=max(UInt64(initial_max_element_tag),maximum(
                    (maximum(tags;init=UInt64(0)) for tags in
                     actual.entity_data.external_element_tags);init=UInt64(0))),
                homology_output=output)
        end
        completed=_reclassify_nodes_plan(source,actual,class,caller;
            node_max=UInt64(initial_max_node_tag),element_max=UInt64(initial_max_element_tag))
        return merge(completed,(homology_output=output,))
    end
    plan=if operation===:refine
        result=_mixed_refine_detached(source,actual,class,max_nodes,max_cells;
            node_max=UInt64(initial_max_node_tag),element_max=UInt64(initial_max_element_tag),
            renumber,save_all,linear_nodes=linear)
        result.plan===nothing && error("$caller: actual refinement lost entity metadata")
        previous=source.meshing.mesh_authority
        carrier_tags=previous.mesh.entity_data.external_node_tags[previous.carrier_rows]
        merge(result.plan,(refine_carrier_tags=carrier_tags,))
    elseif operation===:order
        _set_order_dim01_plan(source,actual,class,element_order,caller;
            initial_max_node_tag,initial_max_element_tag,renumber,save_all,max_nodes,max_cells,linear,incomplete)
    elseif operation in (:renumber_nodes,:renumber_elements)
        kind=operation===:renumber_nodes ? :node : :element
        _tagged_renumber_plan(source,actual,class,(),(),kind,caller;
            overlay=nothing,mid_owners=nothing,
            node_max=initial_max_node_tag,element_max=initial_max_element_tag)
    else
        throw(ArgumentError("$caller: unsupported actual mesh operation $operation"))
    end
    return _higher_rebind_plan_authority!(plan,source.meshing.mesh_authority)
end

Model._MODEL_HIGHER_GENERATOR[]=_generate_higher_frontend_plan

# Read-only native geometry projection. Every support comes from the actual
# parent cell selected by its recorded public identity.
function _higher_geometry_support_key(primary,family,coordinate,order)
    p=order
    weights=if family===:pnt
        (1,)
    elseif family===:lin
        u=coordinate[1]
        (p-u,u)
    elseif family===:tri
        u,v=coordinate
        (p-u-v,u,v)
    elseif family===:tet
        u,v,w=coordinate
        (p-u-v-w,u,v,w)
    elseif family===:qua
        u,v=coordinate
        ((p-u)*(p-v),u*(p-v),u*v,(p-u)*v)
    elseif family===:hex
        u,v,w=coordinate
        ((p-u)*(p-v)*(p-w),u*(p-v)*(p-w),u*v*(p-w),(p-u)*v*(p-w),
         (p-u)*(p-v)*w,u*(p-v)*w,u*v*w,(p-u)*v*w)
    elseif family===:pri
        u,v,w=coordinate
        ((p-u-v)*(p-w),u*(p-w),v*(p-w),(p-u-v)*w,u*w,v*w)
    elseif family===:pyr
        u,v,q=coordinate
        q==0 ? (0,0,0,0,1) : ((q-u)*(q-v),u*(q-v),u*v,(q-u)*v,(p-q)*q)
    else
        throw(ArgumentError("native Mixed geometry: unsupported family $family"))
    end
    # These are exact first-order weights on the integral reference lattice.
    # Normalization also gives the Pyramid apex and base vertices the same
    # representation as every other occurrence of that corner identity.
    divisor=foldl(gcd,weights;init=0)
    terms=Tuple{Int32,Int}[(primary[index],weight÷divisor)
        for (index,weight) in enumerate(weights) if weight!=0]
    sort!(terms;by=first)
    return Tuple(terms)
end

# The legacy carrier owns its interpolation stencil independently of the
# native producer. The caller has already bound this key to the actual public
# node through its parent cell: no coordinate coincidence establishes identity.
function _higher_legacy_carrier_support(model,owner,key,order)
    model===nothing && return false
    if owner[1]==1 && haskey(model.curves,Int(owner[2]))
        return Model._curve_type(model,Int(owner[2])) in (:line,:degenerate)
    end
    order==2 && owner[1]==2 && length(key)==4 || return false
    all(term->term[2]==key[1][2],key) || return false
    surface=Int(owner[2])
    haskey(model.surfaces,surface) && !haskey(model.discrete,(2,surface)) || return false
    geometry=get(model.surface_geometry,surface,nothing)
    return Model._surface_type(model,surface)===:plane &&
        !(geometry!==nothing && hasproperty(geometry,:occ))
end

function _higher_mixed_geometry_view(authority,caller;model=nothing)
    carrier=authority.carrier
    carrier isa MixedMesh && carrier.entity_data===nothing || return carrier
    actual=authority.mesh
    length(authority.carrier_cells)==length(carrier.blocks) || throw(ArgumentError(
        "$caller: native Mixed geometry lacks its actual cell bindings"))
    catalog=Dict{UInt64,Tuple{Int,Vector{Int32}}}()
    for (bi,block) in enumerate(actual.blocks),column in axes(block.nodes,2)
        tag=actual.entity_data.external_element_tags[bi][column]
        catalog[tag]=(block.msh,block.nodes[:,column])
    end
    rows=copy(authority.carrier_rows)
    inverse=Dict(row=>Int32(position) for (position,row) in enumerate(rows))
    points=NTuple{3,Float64}[ntuple(axis->carrier.coords[axis,position],3)
                           for position in axes(carrier.coords,2)]
    blocks=ElementBlock[]
    for (bi,block) in enumerate(carrier.blocks)
        tags=authority.carrier_cells[bi]
        target=isempty(tags) ? block.msh : get(catalog,first(tags),(0,Int32[]))[1]
        target!=0 && all(tag->haskey(catalog,tag) && catalog[tag][1]==target,tags) ||
            throw(ArgumentError("$caller: native Mixed geometry has absent or inconsistent actual parents"))
        family=msh_spec(block.msh).family
        msh_spec(target).family==family || throw(ArgumentError(
            "$caller: native Mixed geometry parent changes element family"))
        primary=Elements.msh_num_nodes(msh_type(family,1))
        size(block.nodes,1)>=primary || throw(ArgumentError(
            "$caller: native Mixed geometry has incomplete primary vertices"))
        spec=msh_spec(target)
        reference=Elements.lagrange_nodes(target)
        lattice=family===:pnt ? [(0,)] : family===:lin ? Elements._line_monomials(spec.order) :
            getfield(Elements,Symbol("_",family,"_monomials"))(spec.order,spec.serendipity)
        nodes=Matrix{Int32}(undef,length(lattice),length(tags))
        supports=Vector{Dict{Tuple,Int32}}(undef,length(tags))
        for (column,tag) in enumerate(tags)
            actual_nodes=catalog[tag][2]
            producer=authority.carrier_rows[block.nodes[1:primary,column]]
            sort(producer)==sort(actual_nodes[1:primary]) || throw(ArgumentError(
                "$caller: native Mixed geometry primary binding disagrees with its actual parent"))
            support=Dict{Tuple,Int32}()
            for (slot,coordinate) in enumerate(lattice)
                key=_higher_geometry_support_key(@view(actual_nodes[1:primary]),family,coordinate,spec.order)
                support[key]=actual_nodes[slot]
            end
            supports[column]=support
        end
        # Match the established native Mixed block/slot/column emission order.
        for (slot,coordinate) in enumerate(lattice),column in eachindex(tags)
            key=_higher_geometry_support_key(authority.carrier_rows[block.nodes[1:primary,column]],
                family,coordinate,spec.order)
            row=get(supports[column],key,Int32(0))
            row>0 || throw(ArgumentError("$caller: native Mixed geometry support loses its actual identity"))
            position=get(inverse,row,Int32(0))
            if position==0
                length(rows)<typemax(Int32) || throw(ArgumentError(
                    "$caller: native Mixed geometry exceeds the Int32 node limit"))
                push!(rows,row);position=Int32(length(rows));inverse[row]=position
                owner=actual.entity_data.node_entities[row]
                point=if _higher_legacy_carrier_support(model,owner,key,spec.order)
                    weights=MeshFunctionSpaces._first_order_values(Val(family),reference[1,slot],
                        reference[2,slot],reference[3,slot],caller,slot)
                    terms=_dim01_support_key(@view(block.nodes[1:primary,column]),weights)
                    ntuple(axis->sum(weight*carrier.coords[axis,node]
                        for (node,weight) in terms),3)
                else
                    ntuple(axis->actual.coords[axis,row],3)
                end
                push!(points,point)
            end
            nodes[slot,column]=position
        end
        push!(blocks,ElementBlock(target,nodes,block.tags))
    end
    coordinates=Matrix{Float64}(undef,3,length(points))
    for position in eachindex(points),axis in 1:3
        coordinates[axis,position]=points[position][axis]
    end
    return MixedMesh(coordinates,blocks;physical_names=carrier.physical_names,
        elementary_entities=carrier.elementary_entities,periodic_links=carrier.periodic_links,
        ancillary_sections=carrier.ancillary_sections,data_sections=carrier.data_sections,
        partition_data=carrier.partition_data)
end
