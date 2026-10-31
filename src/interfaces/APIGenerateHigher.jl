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
        Model._model_native_hybrid_volume_certificate(mesh,source_row_tags,production.constructor,assembly.caller;embedded_pool=production.embedded_pool)
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
            parameter=if data!==nothing
                data.node_parametric[position]
            elseif position<=length(class.node_parametric)
                class.node_parametric[position]
            else
                nothing
            end
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
        implicit=(Model._implicit_volume_surface(m,volume) &&
                  Model._extrude_gate(m,3,volume)===nothing)
        for signed in Model._model_volume_boundary_surfaces(m,volume,caller)
            surface=abs(signed)
            # Analytic primitive Volume producers retain their established
            # boundary encoding. An explicit shell's faces — including
            # OCC-curved tensor/lofted faces, whose standalone production
            # now carries full provenance — mesh as ordinary Surface
            # products upstream of the Region, matching upstream
            # meshGFace-before-meshGRegion so the volume PLC consumes real
            # boundary vertex bindings instead of a coordinate-only re-mesh.
            implicit && push!(skipped,surface)
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
    # Closing rows of a toroidal QuadToTri chain are seam-snapped at the
    # source mesh before emission (`_extrude_volume_columns`,
    # `_extrude_lateral_mesh`), so the certificate must replay the same snap —
    # the raw transform eval differs by evaluator noise from the welded
    # canonical vertex. Copied tops keep the verbatim transform (no snap in
    # `_extrude_top_mesh`).
    seam=if copied
        nothing
    elseif dimension==3
        Model._extrude_quadtri_seam_index(m,tag,caller)
    elseif dimension==2
        begin
            found=nothing
            for region in Model._extrude_quadtri_regions(m,tag,caller)
                found=Model._extrude_quadtri_seam_index(m,region,caller)
                found===nothing || break
            end
            found
        end
    else
        nothing
    end
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
            seam!==nothing && level==length(levels) &&
                (point=Model._extrude_quadtri_snap_seam(seam,point,caller))
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
            "$caller: ordinary Layers entity ($dimension,$tag) omits lower " *
            "cell UID $uid at $point owner $(data.node_entities[row])"))
        if retained[target]
            node_tags[target]==uid || throw(ArgumentError(
                "$caller: ordinary Layers emitter welds distinct lower cell UIDs "*
                "($dimension,$tag) point=$point uid=$uid prior=$(node_tags[target]) "*
                "owner=$(data.node_entities[row])"))
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
    if copied
        # `copyMesh(GFace *from, GFace *to)` carries each copied vertex's
        # stored par() onto the top — the copy's parametric domain IS the
        # source face's (upstream `boundaryLayerSurface` likewise keys its
        # vertex map on source (u,v)). Unbound interior rows therefore
        # inherit the source row's (u,v) bitwise; re-inverting the
        # displaced patch through the generated ruled boundary can fail to
        # converge and is not what upstream stores.
        for column in axes(production.column_rows,2)
            row=production.column_rows[1,column]
            (retained[row] || parameters[row]!==nothing) && continue
            source_row=source_positions[production.source_node_tags[column]]
            parameter=source_product.entity_data.node_parametric[source_row]
            parameter===nothing || (parameters[row]=copy(parameter))
        end
    end
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

# meshGFace writes a recovered boundary chain back onto its GEdge when the
# Face producer subdivides a consumed edge — forced periodic sync vertices
# and vertices sitting exactly on a segment are the producers. The higher
# plan emits a Curve's Line cells before its Surface phase ran, so the
# native producer's ordered UID chains (`production.curve_chains`) reconcile
# here: unchanged chains are left alone, and a grown chain replaces the
# stale dim-1 cells with the actual recovered sub-edges. Fresh element
# identities continue the generation's allocator; the stored parameter
# samples follow so later consumers see the updated discretization.
function _higher_surface_curve_remesh!(assembly,m,chains,curve_births)
    remeshed=Tuple{Int,Int32}[]
    for curve in sort!(collect(keys(chains)))
        owner=(1,curve)
        positions=Int[]
        for (position,cell) in enumerate(assembly.cells)
            cell.owner==owner && cell.msh==1 && push!(positions,position)
        end
        isempty(positions) && continue
        current=UInt64[]
        for (index,position) in enumerate(positions)
            nodes=assembly.cells[position].nodes
            index==1 && push!(current,assembly.nodes[nodes[1]].tag)
            push!(current,assembly.nodes[nodes[2]].tag)
        end
        chain=chains[curve]
        current==chain && continue
        rows=Int32[]
        for tag in chain
            row=get(assembly.node_positions,tag,Int32(0))
            row>0 || throw(ArgumentError(
                "$(assembly.caller): recovered Curve[$curve] vertex $tag " *
                "has no actual node in this generation"))
            push!(rows,row)
        end
        deleteat!(assembly.cells,positions)
        physical=_dim01_physical(m,(1,Int(curve)))
        for index in 1:length(rows)-1
            _dim01_add_cell!(assembly,1,owner,rows[index:index+1],
                _dim01_next_tag!(assembly,false);source_tag=nothing,physical)
        end
        curve_births===nothing || (curve_births[owner]=copy(chain))
        push!(remeshed,owner)
        haskey(m.curves,Int(curve)) || continue
        bounds=Model._model_curve_param_bounds(m,Int(curve),assembly.caller)
        m.curve_params[Int(curve)]=Float64[
            if index==1
                bounds[1]
            elseif index==length(rows)
                bounds[2]
            else
                node=assembly.nodes[row]
                node.owner==owner && node.parameter!==nothing &&
                    !isempty(node.parameter) ? node.parameter[1] :
                    Model.model_parametrization(m,1,Int(curve),
                        collect(node.point))[1]
            end
            for (index,row) in enumerate(rows)]
    end
    return remeshed
end

# Ordered actual vertex positions of one emitted Curve product — the chain
# its Line cells form. `nothing` marks missing, non-linear, branched or
# multi-chain connectivity, which the transport leaves unreconciled.
function _higher_curve_chain_positions(assembly,curve::Int)
    adjacency=Dict{Int32,Vector{Int32}}()
    for cell in assembly.cells
        (cell.owner==(1,Int32(curve)) && cell.msh==1) || continue
        a,b=cell.nodes
        push!(get!(()->Int32[],adjacency,a),b)
        push!(get!(()->Int32[],adjacency,b),a)
    end
    isempty(adjacency) && return nothing
    endpoints=Int32[]
    for (position,neighbors) in adjacency
        length(neighbors)>2 && return nothing
        length(neighbors)==1 && push!(endpoints,position)
    end
    length(endpoints) in (0,2) || return nothing
    start=isempty(endpoints) ? first(keys(adjacency)) : first(endpoints)
    order=Int32[start];previous=Int32(0);current=start
    while true
        next=Int32(0)
        for neighbor in adjacency[current]
            neighbor!=previous && (next=neighbor;break)
        end
        next==0 && break
        push!(order,next);previous,current=current,next
        current==start && break
        length(order)>length(adjacency)+1 && return nothing
    end
    if isempty(endpoints)
        (length(order)==length(adjacency)+1 && last(order)==start) ||
            return nothing
    else
        (length(order)==length(adjacency) && last(order)!=start) ||
            return nothing
    end
    return order
end

# Closed position rings assembled from a patch's boundary Curve chains.
# Closed chains are loops already; open chains join at shared endpoint
# vertices. `nothing` marks a boundary that cannot close into loops.
function _higher_discrete_surface_loops(chains::Vector{Vector{Int32}})
    loops=Vector{Int32}[]
    used=falses(length(chains))
    for (index,chain) in enumerate(chains)
        first(chain)==last(chain) || continue
        used[index]=true
        push!(loops,chain[1:end-1])
    end
    remaining=findall(.!used)
    while !isempty(remaining)
        first_index=popfirst!(remaining);used[first_index]=true
        loop=copy(chains[first_index])
        while first(loop)!=last(loop)
            endpoint=last(loop)
            pick=findfirst(index->first(chains[index])==endpoint ||
                                 last(chains[index])==endpoint,remaining)
            pick===nothing && return nothing
            index=remaining[pick];deleteat!(remaining,pick);used[index]=true
            segment=chains[index]
            first(segment)!=endpoint && (segment=reverse(segment))
            append!(loop,segment[2:end])
        end
        push!(loops,loop[1:end-1])
    end
    return loops
end

# Point-on-segment test within roundoff: used to certify that a produced
# Curve node lies on the retained polyline edge it claims to subdivide.
function _higher_point_on_segment(point,a,b;tol=1e-9)
    abx=b[1]-a[1];aby=b[2]-a[2];abz=b[3]-a[3]
    apx=point[1]-a[1];apy=point[2]-a[2];apz=point[3]-a[3]
    cx=aby*apz-abz*apy;cy=abz*apx-abx*apz;cz=abx*apy-aby*apx
    span=sqrt(abx*abx+aby*aby+abz*abz)
    span>0 || return false
    sqrt(cx*cx+cy*cy+cz*cz)<=tol*span || return false
    t=(apx*abx+apy*aby+apz*abz)/(span*span)
    return -tol<=t<=1+tol
end

# Positions strictly between ring[ia] and ring[ib] walking one direction;
# `forward=false` walks the cyclic ring backwards.
function _higher_chain_arc(ring,ia,ib,forward)
    n=length(ring);segment=Int32[];j=ia
    while j!=ib
        j=forward ? mod1(j+1,n) : mod1(j-1,n)
        j==ib && break
        push!(segment,ring[j])
    end
    return segment
end

# Boundary-edge subdivision transport: when every emitted Curve node lies
# on the retained polyline edge its parameter claims (the discrete Curve
# mesher keeps frame knots as stations), the produced chains are exactly
# subdivisions of the retained boundary edges. Rebuild the patch by
# fanning every retained triangle inside its own footprint — each new
# triangle is a sub-triangle of the retained cell it replaces, so the
# produced patch cannot self-intersect, keeps the retained winding and
# component structure, and every neighbour emits the identical sub-chain
# on the shared Curve. Returns `(triangles, boundary, interior)` in the
# `vcat(boundary,interior)` index space, or `nothing` when a produced node
# sits off the retained polyline (free resampling) or the topology does
# not refine cleanly — the caller then falls back to embedded remeshing.
function _higher_discrete_surface_subdivision(assembly,m,record,curves,chains,owner)
    fail=reason->nothing
    coordinates=Model._model_record_node_coords(m)
    new_position=Dict{UInt64,Int32}()
    subdivided=Dict{NTuple{2,UInt64},Tuple{UInt64,UInt64,Vector{Int32}}}()
    for (curve,chain) in zip(curves,chains)
        curve_record=get(m.discrete,(1,curve),nothing)
        curve_record===nothing && (return fail("no discrete record for curve $curve"))
        edges=NTuple{2,UInt64}[]
        for (index,msh) in enumerate(curve_record.element_types)
            msh_family(msh)===:lin || continue
            nodes=curve_record.element_nodes[index]
            length(nodes)==2 ||
                (return fail("curve $curve non-2 node element"))
            push!(edges,(nodes[1],nodes[2]))
        end
        isempty(edges) && (return fail("curve $curve no edges"))
        old_chain=Model._model_discrete_curve_chain(edges,assembly.caller,curve)
        old_closed=first(old_chain)==last(old_chain)
        old_closed && (old_chain=old_chain[1:end-1])
        length(old_chain)>=2 ||
            (return fail("curve $curve old_chain<2"))
        produced_closed=first(chain)==last(chain)
        ring=produced_closed ? chain[1:end-1] : chain
        length(ring)>=2 ||
            (return fail("curve $curve ring<2"))
        position_of=Dict(position=>index for (index,position) in enumerate(ring))
        length(position_of)==length(ring) ||
            (return fail("curve $curve duplicate ring positions"))
        # The curve mesher places knots by parametric evaluation, so a
        # produced station meant to reproduce a retained polyline vertex
        # can land a few ulps off the stored coordinate. Match each old
        # vertex to its nearest produced ring position within the same
        # 1e-9 relative tolerance `_higher_point_on_segment` certifies
        # stations with; a genuinely dropped vertex stays a hard failure.
        scale=0.0
        for position in ring
            point=assembly.nodes[position].point
            scale=max(scale,abs(point[1]),abs(point[2]),abs(point[3]))
        end
        tol=1e-9*max(scale,1.0)
        for node in old_chain
            haskey(coordinates,node) ||
                (return fail("curve $curve old node $node no coords"))
            target=coordinates[node]
            best=Int32(0);distance=Inf
            for position in ring
                point=assembly.nodes[position].point
                d=(point[1]-target[1])^2+(point[2]-target[2])^2+
                   (point[3]-target[3])^2
                d<distance && (distance=d;best=position)
            end
            distance<=tol*tol ||
                (return fail("curve $curve old node $node at $target not on produced ring"))
            new_position[node]=best
        end
        n=length(ring)
        edge_count=old_closed ? length(old_chain) : length(old_chain)-1
        for index in 1:edge_count
            a=old_chain[index]
            b=old_chain[mod1(index+1,length(old_chain))]
            pa=coordinates[a];pb=coordinates[b]
            ia=position_of[new_position[a]];ib=position_of[new_position[b]]
            segment=nothing
            for forward in (true,false)
                produced_closed || forward==(ia<ib) || continue
                arc=_higher_chain_arc(ring,ia,ib,forward)
                all(position->_higher_point_on_segment(
                        assembly.nodes[position].point,pa,pb),arc) || continue
                segment=arc;break
            end
            segment===nothing &&
                (return fail("curve $curve edge index=$index a=$a b=$b no on-segment arc"))
            # Every produced station is a born Curve vertex the patch must
            # keep referencing — dropping one orphans it against the curve
            # product even when it merely sits close to an endpoint. The
            # emitted polygon therefore keeps all stations; a truly
            # coincident duplicate would surface through the degenerate
            # guard rather than be silently elided here.
            key=a<b ? (a,b) : (b,a)
            haskey(subdivided,key) &&
                (return fail("curve $curve edge $key claimed twice"))
            subdivided[key]=(a,b,segment)
        end
    end
    # Retained cells whose boundary edges are not claimed by any Curve
    # chain cannot be rebuilt by subdivision.
    edge_tris=Dict{NTuple{2,UInt64},Int}()
    retained=NTuple{3,Int32}[]
    for cell in assembly.cells
        (cell.owner==owner && cell.msh==2) || continue
        push!(retained,(cell.nodes[1],cell.nodes[2],cell.nodes[3]))
        for (x,y) in ((cell.nodes[1],cell.nodes[2]),
                      (cell.nodes[2],cell.nodes[3]),
                      (cell.nodes[3],cell.nodes[1]))
            ux=assembly.nodes[x].tag;uy=assembly.nodes[y].tag
            key=ux<uy ? (ux,uy) : (uy,ux)
            edge_tris[key]=get(edge_tris,key,0)+1
        end
    end
    for (key,count) in edge_tris
        count==1 && !haskey(subdivided,key) &&
            (return fail("retained boundary edge $key unclaimed"))
    end
    produced=NTuple{3,Int32}[]
    boundary=Int32[];in_boundary=Set{Int32}()
    for chain in chains,position in chain
        position in in_boundary && continue
        push!(in_boundary,position);push!(boundary,position)
    end
    local_index=Dict{Int32,Int32}(position=>Int32(index)
                                  for (index,position) in enumerate(boundary))
    interior=Int32[]
    for node in record.node_tags
        position=get(assembly.node_positions,node,Int32(0))
        (position==0 || haskey(local_index,position)) && continue
        local_index[position]=Int32(length(boundary)+length(interior)+1)
        push!(interior,position)
    end
    triangles=NTuple{3,Int32}[]
    for cell_nodes in retained
        corners=Int32[]
        chains3=Vector{Int32}[]
        for index in 1:3
            x=cell_nodes[index];y=cell_nodes[mod1(index+1,3)]
            ux=assembly.nodes[x].tag;uy=assembly.nodes[y].tag
            push!(corners,get(new_position,ux,x))
            key=ux<uy ? (ux,uy) : (uy,ux)
            entry=get(subdivided,key,nothing)
            chain=Int32[]
            if entry!==nothing
                a,b,segment=entry
                append!(chain,a==ux ? segment : reverse(segment))
            end
            push!(chains3,chain)
        end
        for tri in _higher_subdivided_cell_triangles(corners,chains3)
            push!(triangles,(local_index[tri[1]],local_index[tri[2]],
                             local_index[tri[3]]))
        end
    end
    # The peel keeps every emitted triangle on a different retained-edge
    # line than its fan apex, so it is non-degenerate by construction;
    # guard anyway so a pathological input falls back to the embedded
    # remesh instead of handing the strict Region certificate a bad cell.
    points=Vector{NTuple{3,Float64}}(undef,length(boundary)+length(interior))
    for (index,position) in enumerate(boundary)
        points[index]=assembly.nodes[position].point
    end
    for (index,position) in enumerate(interior)
        points[length(boundary)+index]=assembly.nodes[position].point
    end
    for tri in triangles
        a,b,c=points[tri[1]],points[tri[2]],points[tri[3]]
        ux,uy,uz=b.-a;vx,vy,vz=c.-a
        cx=uy*vz-uz*vy;cy=uz*vx-ux*vz;cz=ux*vy-uy*vx
        cx*cx+cy*cy+cz*cz==0.0 &&
            (return fail("produced degenerate triangle $tri"))
    end
    return triangles,boundary,interior
end

# Triangulate one retained cell's subdivided boundary polygon. `corners`
# are the produced positions of the cell's three vertices in polygon order
# and `chains[i]` the produced stations on the straight retained edge
# corners[i] → corners[i+1]; every vertex sits on the retained cell's
# convex footprint. A plain fan from the first corner emits a zero-area
# triangle for every consecutive station pair on the apex's own two edges,
# so chained sides are peeled against the next chained side's first
# station instead: each peel covers the sub-polygon (run ++ apex) bounded
# by ring edges plus one interior chord, preserves the retained winding,
# and puts the fan apex on a different retained edge's line than the run —
# no emitted triangle is degenerate. When a single chained side remains it
# fans to the one corner off its line.
function _higher_subdivided_cell_triangles(corners::Vector{Int32},
                                           chains::Vector{Vector{Int32}})
    cs=copy(corners)
    ch=Vector{Int32}[copy(chains[1]),copy(chains[2]),copy(chains[3])]
    triangles=NTuple{3,Int32}[]
    while count(!isempty,ch)>=2
        i=findfirst(k->!isempty(ch[k]) && !isempty(ch[mod1(k+1,3)]),1:3)
        apex=ch[mod1(i+1,3)][1]
        run=Int32[cs[i]];append!(run,ch[i]);push!(run,cs[mod1(i+1,3)])
        for k in 1:length(run)-1
            push!(triangles,(run[k],run[k+1],apex))
        end
        cs=Int32[cs[i],apex,cs[mod1(i+2,3)]]
        ch=Vector{Int32}[Int32[],ch[mod1(i+1,3)][2:end],ch[mod1(i+2,3)]]
    end
    j=findfirst(!isempty,ch)
    if j===nothing
        push!(triangles,(cs[1],cs[2],cs[3]))
    else
        apex=cs[mod1(j+2,3)]
        run=Int32[cs[j]];append!(run,ch[j]);push!(run,cs[mod1(j+1,3)])
        for k in 1:length(run)-1
            push!(triangles,(run[k],run[k+1],apex))
        end
    end
    return triangles
end

# Re-triangulate one `CreateGeometry`-parametrized discrete surface against
# its actual boundary products. The patch is re-embedded on its
# least-squares (u,v) frame and triangulated with the emitted boundary
# chains as constrained loops; every produced cell replaces the retained
# element set under fresh born identities. Returns false when the patch
# cannot be reconciled — the retained cells then carry the unchanged
# product into the strict Region certificate.
function _higher_discrete_surface_transport_one!(assembly,m,tag,record;caller,
                                               geometry_cells=nothing)
    owner=(2,Int32(tag))
    curves=Int[]
    for (dim,bound) in record.boundary
        dim==1 || continue
        push!(curves,abs(bound))
    end
    isempty(curves) && (return nothing,false)
    all(entry->entry==2,record.element_types) || (return nothing,false)
    chains=Vector{Int32}[]
    for curve in curves
        chain=_higher_curve_chain_positions(assembly,curve)
        chain===nothing && (return nothing,false)
        push!(chains,chain)
    end
    # Preferred transport: when the emitted chains refine the retained
    # polyline edge-by-edge, rebuild the patch by boundary-edge
    # subdivision — no retriangulation chords, no possible piercing, and
    # every neighbouring patch emits the identical sub-chain.
    subdivided=_higher_discrete_surface_subdivision(
        assembly,m,record,curves,chains,owner)
    if subdivided!==nothing
        triangles,boundary,interior=subdivided
        width=length(boundary)
        filter!(cell->!(cell.owner==owner && msh_dimension(cell.msh)==2),
                assembly.cells)
        physical=_dim01_physical(m,(2,tag))
        first_emitted=length(assembly.cells)+1
        for tri in triangles
            nodes=Int32[tri[slot]<=width ? boundary[tri[slot]] :
                        interior[tri[slot]-width] for slot in 1:3]
            _dim01_add_cell!(assembly,2,owner,nodes,
                _dim01_next_tag!(assembly,false);source_tag=nothing,physical)
        end
        emitted_cells=assembly.cells[first_emitted:end]
        used=falses(length(boundary)+length(interior))
        for tri in triangles,index in tri
            used[index]=true
        end
        positions=Int32[];remap=zeros(Int32,length(used))
        for index in eachindex(used)
            used[index] || continue
            remap[index]=Int32(length(positions)+1)
            push!(positions,index<=width ? boundary[index] : interior[index-width])
        end
        coords=Matrix{Float64}(undef,3,length(positions))
        for (column,position) in enumerate(positions),axis in 1:3
            coords[axis,column]=assembly.nodes[position].point[axis]
        end
        tri_columns=Matrix{Int32}(undef,3,length(triangles))
        for (column,tri) in enumerate(triangles),slot in 1:3
            tri_columns[slot,column]=remap[tri[slot]]
        end
        mesh=Mesh(coords;tris=tri_columns,
                  tri_tag=fill(physical,length(triangles)))
        geometry_cells===nothing || (geometry_cells[mesh]=
            Vector{Union{Nothing,_Dim01Cell}}[Union{Nothing,_Dim01Cell}[],
            Union{Nothing,_Dim01Cell}[emitted_cells...],
            Union{Nothing,_Dim01Cell}[]])
        class=_higher_geometry_part_class(mesh,2,tag,
            [assembly.nodes[position].owner for position in positions])
        return (tag,mesh,class,positions),true
    end
    position_loops=_higher_discrete_surface_loops(chains)
    (position_loops===nothing || any(loop->length(loop)<3,position_loops)) &&
        (return nothing,false)
    # Retained patch cells: the reference triangulation (winding + shape)
    # and its directed boundary. The produced traversal of a constrained
    # loop follows the loop order on the patch's outer ring and its
    # reverse on hole rings, so each ring is oriented to reproduce the
    # retained boundary traversal — adjacent products then emit their
    # shared Curve in opposite directions and the merged shell stays
    # orientable. Vertices the Curve remesh replaced no longer neighbour
    # each other on the new rings; the vote instead compares the retained
    # boundary's cyclic order restricted to the ring's retained vertices.
    reference=NTuple{3,NTuple{3,Float64}}[]
    ref_dir=Dict{NTuple{2,Int32},NTuple{2,Int32}}();ref_inc=Dict{NTuple{2,Int32},Int}()
    for cell in assembly.cells
        (cell.owner==owner && cell.msh==2) || continue
        push!(reference,ntuple(slot->assembly.nodes[cell.nodes[slot]].point,3))
        for (x,y) in ((cell.nodes[1],cell.nodes[2]),
                      (cell.nodes[2],cell.nodes[3]),
                      (cell.nodes[3],cell.nodes[1]))
            key=x<y ? (x,y) : (y,x)
            ref_inc[key]=get(ref_inc,key,0)+1
            ref_dir[key]=(x,y)
        end
    end
    # The Curve remesh replaces every chain vertex, so the produced and
    # retained rings share no vertices; instead each produced ring edge
    # votes against the retained boundary edge nearest its midpoint —
    # the retained direction along that Curve is the emitted traversal
    # the neighbour's own retained boundary traverses oppositely.
    ret_mid=NTuple{3,Float64}[];ret_dir=NTuple{3,Float64}[]
    for (key,count) in ref_inc
        count==1 || continue
        x,y=ref_dir[key]
        push!(ret_mid,(assembly.nodes[x].point.+assembly.nodes[y].point)./2)
        push!(ret_dir,assembly.nodes[y].point.-assembly.nodes[x].point)
    end
    perimeters=Vector{Float64}(undef,length(position_loops))
    for (index,loop) in enumerate(position_loops)
        total=0.0
        for i in eachindex(loop)
            a=assembly.nodes[loop[i]].point
            b=assembly.nodes[loop[mod1(i+1,length(loop))]].point
            total+=sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2)
        end
        perimeters[index]=total
    end
    outer_loop=argmax(perimeters)
    for (index,loop) in enumerate(position_loops)
        isempty(ret_mid) && break
        vote=0.0
        for i in eachindex(loop)
            a=assembly.nodes[loop[i]].point
            b=assembly.nodes[loop[mod1(i+1,length(loop))]].point
            mid=(a.+b)./2;best=0;dist=Inf
            for r in eachindex(ret_mid)
                d=(ret_mid[r][1]-mid[1])^2+(ret_mid[r][2]-mid[2])^2+
                   (ret_mid[r][3]-mid[3])^2
                d<dist && (dist=d;best=r)
            end
            vote+=ret_dir[best][1]*(b[1]-a[1])+ret_dir[best][2]*(b[2]-a[2])+
                  ret_dir[best][3]*(b[3]-a[3])
        end
        outer=index==outer_loop
        vote==0.0 && continue
        # Outer rings emit the loop order, hole rings its reverse; the
        # loop order should therefore carry sign(retained)·emission.
        outer==(vote>0) && continue
        position_loops[index]=reverse(loop)
    end
    boundary=Int32[];seen=Set{Int32}()
    for loop in position_loops,position in loop
        position in seen || (push!(seen,position);push!(boundary,position))
    end
    local_index=Dict(position=>Int32(index) for (index,position) in enumerate(boundary))
    loops=Vector{Int32}[Int32[local_index[position] for position in loop]
                       for loop in position_loops]
    # Retained vertices that sat on a remeshed Curve were replaced by a
    # produced knot at the identical coordinate; their old position is off
    # the ring yet still listed in node_tags. Passing it as an interior
    # Steiner point makes the patch fold onto its own boundary — the
    # neighbour emits the ring vertex, this patch emits the duplicate —
    # so drop every candidate whose point a loop position already covers.
    interior=Int32[]
    claimed=Set{NTuple{3,Float64}}(ntuple(axis->
            Model._model_projection_coordinate_key(
                assembly.nodes[position].point[axis]),3)
        for position in boundary)
    for uid in record.node_tags
        position=get(assembly.node_positions,uid,Int32(0))
        (position==0 || position in seen) && continue
        push!(seen,position)
        point=ntuple(axis->Model._model_projection_coordinate_key(
                     assembly.nodes[position].point[axis]),3)
        point in claimed && continue
        push!(claimed,point);push!(interior,position)
    end
    # Retained-incidence edges of the patch's own vertex set drive the
    # harmonic embedding fallback for folded projections.
    vertex_index=Dict{UInt64,Int32}()
    for (index,position) in enumerate(boundary)
        vertex_index[assembly.nodes[position].tag]=Int32(index)
    end
    for (index,position) in enumerate(interior)
        vertex_index[assembly.nodes[position].tag]=Int32(length(boundary)+index)
    end
    interior_edges=NTuple{2,Int32}[]
    for nodes in record.element_nodes
        length(nodes)==3 || continue
        for edge in ((1,2),(2,3),(3,1))
            first=get(vertex_index,nodes[edge[1]],Int32(0))
            second=get(vertex_index,nodes[edge[2]],Int32(0))
            (first==0 || second==0) && continue
            push!(interior_edges,(first,second))
        end
    end
    frame=NTuple{3,Float64}[]
    for uid in Model._record_referenced_nodes(record,2)
        position=get(assembly.node_positions,uid,Int32(0))
        position==0 && (return nothing,false)
        push!(frame,assembly.nodes[position].point)
    end
    expected=_higher_record_component_count(record)
    triangles=Model._model_parametrized_discrete_surface_remesh(
        frame,[assembly.nodes[position].point for position in boundary],
        loops,[assembly.nodes[position].point for position in interior],
        interior_edges,reference,expected)
    triangles===nothing && (return nothing,false)
    width=length(boundary)
    filter!(cell->!(cell.owner==owner && msh_dimension(cell.msh)==2),assembly.cells)
    physical=_dim01_physical(m,(2,tag))
    first_emitted=length(assembly.cells)+1
    for tri in triangles
        nodes=Int32[tri[slot]<=width ? boundary[tri[slot]] :
                    interior[tri[slot]-width] for slot in 1:3]
        _dim01_add_cell!(assembly,2,owner,nodes,
            _dim01_next_tag!(assembly,false);source_tag=nothing,physical)
    end
    emitted_cells=assembly.cells[first_emitted:end]
    # The produced product lists only the vertices the emitted cells
    # actually reference, in remesh order (boundary then interior). Rows
    # the triangulation left unused — merged projections, hole interiors —
    # stay owned orphans that `_dim01_finish!` drops through `discarded`;
    # binding this part as producer rows keeps them out of the product
    # view a Region PLC merge reads, exactly like an ordinary remesh.
    used=falses(length(boundary)+length(interior))
    for tri in triangles,index in tri
        used[index]=true
    end
    positions=Int32[];remap=zeros(Int32,length(used))
    for index in eachindex(used)
        used[index] || continue
        remap[index]=Int32(length(positions)+1)
        push!(positions,index<=width ? boundary[index] : interior[index-width])
    end
    coords=Matrix{Float64}(undef,3,length(positions))
    for (column,position) in enumerate(positions),axis in 1:3
        coords[axis,column]=assembly.nodes[position].point[axis]
    end
    tri_columns=Matrix{Int32}(undef,3,length(triangles))
    for (column,tri) in enumerate(triangles),slot in 1:3
        tri_columns[slot,column]=remap[tri[slot]]
    end
    mesh=Mesh(coords;tris=tri_columns,tri_tag=fill(physical,length(triangles)))
    geometry_cells===nothing || (geometry_cells[mesh]=
        Vector{Union{Nothing,_Dim01Cell}}[Union{Nothing,_Dim01Cell}[],
        Union{Nothing,_Dim01Cell}[emitted_cells...],Union{Nothing,_Dim01Cell}[]])
    class=_higher_geometry_part_class(mesh,2,tag,
        [assembly.nodes[position].owner for position in positions])
    return (tag,mesh,class,positions),true
end

# Component count of a record's stored triangle set under shared-edge
# adjacency — the region structure a reconciled product must keep.
function _higher_record_component_count(record)
    edge_tris=Dict{NTuple{2,UInt64},Vector{Int}}()
    for (cell,nodes) in enumerate(record.element_nodes)
        length(nodes)==3 || continue
        for (a,b) in ((nodes[1],nodes[2]),(nodes[2],nodes[3]),(nodes[3],nodes[1]))
            key=a<b ? (a,b) : (b,a)
            push!(get!(()->Int[],edge_tris,key),cell)
        end
    end
    seen=falses(length(record.element_nodes));count=0
    for cell in eachindex(record.element_nodes)
        length(record.element_nodes[cell])==3 || continue
        seen[cell] && continue
        count+=1
        stack=Int[cell];seen[cell]=true
        while !isempty(stack)
            nodes=record.element_nodes[pop!(stack)]
            for (a,b) in ((nodes[1],nodes[2]),(nodes[2],nodes[3]),(nodes[3],nodes[1]))
                key=a<b ? (a,b) : (b,a)
                for neighbor in edge_tris[key]
                    seen[neighbor] || (seen[neighbor]=true;push!(stack,neighbor))
                end
            end
        end
    end
    return count
end

# `discreteFace::mesh` forwards to the ordinary face mesher once the record
# carries a `CreateGeometry` parametrization — the retained element set
# cannot carry the remeshed boundary vertices a Region PLC assembly reads
# through the surface products. Reconcile each such patch against its
# actual boundary chains before the phase builds its input authority;
# visibility/`MeshOnlyEmpty` keeps the upstream skip gates, and an
# unresolvable patch keeps its retained cells so the strict Region
# certificate reports the unchanged missing identity. A reconciled patch
# counts as remeshed: its owner joins `discarded` so `_dim01_finish!`
# drops the unreferenced owned orphans, and its produced part joins
# `parts` so `_higher_bind_producer_rows` binds the used-row product view.
function _higher_discrete_surface_transport!(assembly,m,entities;options,
                                             discarded=nothing,parts=nothing,
                                             geometry_cells=nothing)
    for (dim,tag) in entities
        dim==2 || continue
        record=get(m.discrete,(2,tag),nothing)
        (record===nothing || !record.geometry_parametrized) &&
            (continue)
        options.mesh_only_visible &&
            model_entity_visibility(m,2,tag)==0 && continue
        owner=(2,Int32(tag))
        options.mesh_only_empty &&
            any(cell->cell.owner==owner,assembly.cells) && continue
        part,ok=_higher_discrete_surface_transport_one!(
            assembly,m,tag,record;caller=assembly.caller,geometry_cells)
        ok || continue
        discarded===nothing || push!(discarded,owner)
        parts===nothing || push!(parts,part)
    end
    return nothing
end

# The dim-2 keep/mesh scan `_higher_surface_phase!` runs before its unit
# merge: fully discrete, MeshOnlyVisible-hidden and MeshOnlyEmpty-occupied
# Faces retain their cells; every other Face meshes and sheds its stale
# products.  PAR01 runs the same scan ahead of the Curve phase as well (the
# keep predicates are phase-stable — the Curve phase only appends dim-0/1
# cells, so a Face's `have` membership and visibility do not move) so that
# detached Surfaces can spawn as soon as their boundary Curve products land.
function _higher_surface_mesh_scan(m,entities,options,have,discarded,caller;
                                   dimension=3)
    retained=Set{Tuple{Int,Int32}}()
    meshed=Int[]
    analytic_boundaries=_higher_analytic_boundary_surfaces(m,dimension,caller)
    for (dim,tag) in entities
        dim==2 || continue
        owner=(2,Int32(tag))
        if tag in analytic_boundaries
            # These native Volume producers retain their established analytic
            # boundary encoding. Model does not yet admit a curved face as an
            # independent Surface product; do not turn a supported Volume
            # request into that standalone-face diagnostic.
            owner in have && throw(ArgumentError(
                "$caller: actual curved Surface[$tag] cells cannot be consumed " *
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
    return meshed,retained
end

# The dim-3 keep/mesh scan `_higher_volume_phase!` runs before its merge: a
# Region meshes unless it is fully discrete, MeshOnlyVisible-hidden or a
# MeshOnlyEmpty-occupied non-swept primitive.  PAR01 runs the same scan inside
# the Surface merge loop so eligible Regions can spawn while Surfaces still
# merge; the Region phase re-derives it verbatim before adopting the slots.
function _higher_volume_mesh_tags(m,entities,options,have)
    tags=Int[]
    for (dim,tag) in entities
        dim==3 || continue
        owner=(3,Int32(tag))
        visible=!(options.mesh_only_visible && model_entity_visibility(m,3,tag)==0)
        swept=Model._extrude_gate(m,3,tag)!==nothing
        keep=_dim01_fully_discrete(m,(3,tag)) ||
             (!swept && (!visible || (options.mesh_only_empty && owner in have)))
        keep || push!(tags,tag)
    end
    return tags
end

# A `_Dim01Assembly` view sharing the live node vector but presenting a
# projected cell list.  The Curve phase appends fresh products for a remeshed
# Curve before its end-of-phase compaction drops the stale source cells, so an
# authority built mid-phase must be derived over the projected post-compaction
# cell set — the exact list the Surface phase's own `_higher_assembly_authority`
# call will read.  The node vector is shared (cells only ever select rows) and
# the view is never used for appends.
function _higher_assembly_view(assembly,cells)
    return _Dim01Assembly(assembly.nodes,cells,assembly.node_positions,
        assembly.node_max,assembly.cell_max,assembly.max_nodes,
        assembly.max_cells,assembly.caller,assembly.node_counter_wrapped,
        assembly.cell_counter_wrapped)
end

# The Surface-arm authority a detached produce must see, rebuilt mid-Curve-
# phase over the projected post-compaction cell set.  Mirrors the phase's
# `input` construction exactly (ranks, producer rows with `curve_births`
# replacement, curve-birth row binding).  Returns `nothing` when the rebuild
# throws (a retained producer losing actual rows is only fatal when the
# serial phase itself reaches the same call — skipped spawns simply rerun).
function _higher_surface_spawn_input(assembly,m,done,live_cells;ranks,
        old_supports,surface_parts,curve_births)
    return try
        view=_higher_assembly_view(assembly,live_cells)
        ranks===nothing || _higher_complete_ranks!(ranks,view)
        authority=_higher_assembly_authority(view,m,done;ranks,old_supports)
        authority=_higher_bind_producer_rows(authority,view,
            surface_parts===nothing ? Tuple[] : surface_parts,
            m.meshing.mesh_authority;curve_births)
        Model._model_mesh_with_curve_births(authority,curve_births;
            caller=assembly.caller)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

# The Region-arm authority a detached produce must see, rebuilt mid-Surface-
# merge exactly like the Region phase's `input` (actual lower products,
# producer-row birth bindings for every merged part, then the surface
# constructor-input table).  `nothing` on a rebuilding failure — the serial
# path re-attempts the identical construction at Region-phase entry.
function _higher_volume_spawn_input(assembly,m,done;ranks,old_supports,
        surface_parts,constructor_inputs)
    return try
        ranks===nothing || _higher_complete_ranks!(ranks,assembly)
        authority=_higher_assembly_authority(assembly,m,done;ranks,old_supports)
        authority=_higher_bind_producer_rows(authority,assembly,
            surface_parts===nothing ? Tuple[] : surface_parts,
            m.meshing.mesh_authority)
        Model._model_mesh_with_surface_constructor_inputs(authority,
            constructor_inputs;caller=assembly.caller)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

function _higher_surface_phase!(assembly,m,options,entities,discarded,done;
                                size_field=nothing,metric_surface=nothing,
                                source_nodes=length(assembly.nodes),ranks=nothing,
                                old_supports=nothing,geometry_parts=nothing,dimension=2,
                                 geometry_cells=nothing,geometry_birth=nothing,
                                 geometry_alternates=nothing,curve_births=Dict{Tuple{Int,Int32},Vector{UInt64}}(),
                                 periodic_pairs=nothing,native_orphan_births=nothing,native_recombine_events=nothing,
                                 surface_constructor_inputs=nothing,
                                 surface_forward=nothing,optimize_pyramids::Float64=0.0)
    # Mesh2D resets all face statuses before OnlyEmpty/OnlyVisible dispatch.
    filter!(owner->owner[1]!=2,done)
    have=Set(cell.owner for cell in assembly.cells)
    meshed,retained=_higher_surface_mesh_scan(m,entities,options,have,
        discarded,assembly.caller;dimension)
    filter!(cell->cell.source_tag===nothing || !(cell.owner in discarded),assembly.cells)
    _higher_discrete_surface_transport!(assembly,m,entities;options,
        discarded,parts=geometry_parts,geometry_cells)
    # The actual products are supplied to extrusion source/cap consumers. A
    # retained surface's cell order is its stored order, including cap pairs.
    ranks===nothing || _higher_complete_ranks!(ranks,assembly)
    input=_higher_assembly_authority(assembly,m,done;ranks,old_supports)
    input=_higher_bind_producer_rows(input,assembly,
        geometry_parts===nothing ? Tuple[] : geometry_parts,
        m.meshing.mesh_authority;curve_births)
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
    # PAR01 — `General.NumThreads` detaches runs of independent singleton
    # Surfaces: each eligible unit produces its receipt on a transaction copy
    # under `Threads.@threads`, and the unchanged serial merge below replays or
    # re-runs each slot in unit order.  `nothing` keeps the exact serial path.
    parallel_units=_parallel_surface_prepare(m,units,assembly.caller;
        metric_surface,scope,periodic_masters,ordinary_sources,
        unmerged_nonew_sources,size_field,geometry_birth)
    parallel_results=parallel_units===nothing ? nothing :
        fill!(Vector{Any}(undef,length(units)),nothing)
    parallel_dirty=Set{Tuple{Int,Int32}}()
    parallel_precise=parallel_units===nothing ? false :
        GmshLibm._win_cli_precision()
    # Cross-phase adoption: `_dim01_curve_phase!` may already have spawned
    # detached Surfaces once their boundary Curve products landed.  Adopt
    # exactly the slots this phase's real eligibility scan still admits
    # (footprint/deps/scope are phase-stable, so the adopted work records are
    # interchangeable) and fold its ledger entries into the merge-time dirty
    # set so every mark the workers observed — or missed — is honoured.
    if parallel_units!==nothing && surface_forward!==nothing &&
            surface_forward.units==units
        for (ui,work) in enumerate(parallel_units)
            work===nothing && continue
            spawned=surface_forward.results[ui]
            spawned===nothing && continue
            parallel_results[ui]=spawned
        end
        union!(parallel_dirty,surface_forward.dirty)
        parallel_precise=surface_forward.precise
    end
    parallel_last=parallel_units===nothing ? 0 :
        something(findlast(!isnothing,parallel_units),0)
    # Cross-phase scheduling the other way: precompute the Region eligibility
    # plan so eligible Regions spawn inside the merge loop below while later
    # Surfaces still produce.  The Region phase re-verifies each adopted slot
    # against its real `extruded`/`staged`/hybrid gates and discards whatever
    # it cannot adopt.  `parallel_dirty` is shared with the spawn ledger so a
    # Region spawned now still observes every later Surface mutation.
    volume_forward=nothing
    if dimension==3
        volume_tags=sort!(_higher_volume_mesh_tags(m,entities,options,have))
        volume_plan=isempty(volume_tags) ? nothing :
            _parallel_volume_prepare(m,volume_tags,input,assembly.caller;
                size_field,scope=nothing,extruded=nothing,
                hybrid_batch=nothing,staged=false)
        volume_plan===nothing || (volume_forward=(
            tags=volume_tags,prepared=volume_plan[1],
            pos_draws=volume_plan[2],
            results=fill!(Vector{Any}(undef,length(volume_tags)),nothing),
            dirty=parallel_dirty,precise=GmshLibm._win_cli_precision(),
            size_field=size_field))
    end
    for (ui,tags) in enumerate(units)
        work=parallel_units===nothing ? nothing : parallel_units[ui]
        if work!==nothing && parallel_results[ui]===nothing
            vend=ui
            while vend<length(units) && parallel_units[vend+1]!==nothing
                vend+=1
            end
            _parallel_surface_fill!(parallel_results,ui:vend,m,parallel_units,
                input,options,assembly,parallel_dirty,parallel_precise,
                geometry_birth;size_field)
        end
        refs=ui<parallel_last || volume_forward!==nothing ?
            _parallel_surface_param_refs(m) : nothing
        products=metric_surface!==nothing ? [(only(tags),mesh_model_surface(m,only(tags),size_field;
            max_metric_edge=metric_surface.max_metric_edge,max_nodes=metric_surface.max_nodes,
            max_triangles=metric_surface.max_triangles,_mesh_input=input,_mesh1d_options=options,
            _mesh_node_high_water=assembly.node_max,_mesh_element_high_water=assembly.cell_max,
            _return_surface_production=true,_boundary_regrade=false))] : work!==nothing &&
            parallel_results[ui]!==nothing ?
            _parallel_surface_products(m,work,parallel_results[ui],
                parallel_dirty,assembly;input,scope,options,size_field,
                geometry_birth,caller=assembly.caller) :
            length(tags)==1 ? [(tags[1],mesh_model_surface(m,tags[1];
            size_field=size_field,_extrude_scope=scope,_mesh_input=input,_mesh1d_options=options,
            _mesh_node_high_water=assembly.node_max,_return_surface_production=true,
            _mesh_element_high_water=assembly.cell_max,_boundary_regrade=false,
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
            if source_binding!==nothing
                # A grown recovered chain invalidates every detached produce
                # still waiting on the consumed Curve's stored discretization.
                for owner in _higher_surface_curve_remesh!(
                        assembly,m,source_binding.curve_chains,curve_births)
                    push!(parallel_dirty,owner)
                end
            end
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
                    # A re-meshed Surface keeps its predecessor rows in the
                    # merge-order product (owned vertices with no remaining
                    # cell incidence), while `_model_mesh_with_surface_constructor_inputs`
                    # re-validates every stored constructor against the
                    # producer-order view. Bind the input against that same
                    # view whenever the merge-order product can carry rows the
                    # production did not emit.
                    merged_product=get(input.products,(2,Int32(tag)),nothing)
                    actual_product=Model._model_mesh_product(input,2,tag;
                        caller=assembly.caller,
                        producer_order=source_binding.constructor_mesh!==nothing || transfinite_events ||
                            metric_surface!==nothing ||
                            (merged_product!==nothing &&
                             length(merged_product.rows)!=length(source_binding.node_tags)))
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
        # State this unit left behind that a later detached produce still
        # reads: stored discretizations (writeback, grading, remesh) and
        # relocated source nodes — remesh owners already landed above.
        refs===nothing || _parallel_surface_mark_dirty!(parallel_dirty,m,
            assembly,products,refs)
        # Cross-phase spawn: once every dependency a Region can consult has a
        # product — its boundary Surfaces merged (in `done`), its embedded
        # entities present — detach its produce on a worker now rather than
        # inside the Region phase's own batch arm.  The spawn-time authority
        # rebuilds `input` exactly the way the Region phase will, so a fresh
        # receipt consumes identical boundary/embedded rows.
        if volume_forward!==nothing
            spawn_input=nothing
            for (vi,vwork) in enumerate(volume_forward.prepared)
                vwork===nothing && continue
                volume_forward.results[vi]===nothing || continue
                all(dep->dep[1]!=2 || dep in done ||
                    haskey(input.products,dep),vwork.deps) || continue
                if spawn_input===nothing
                    spawn_input=_higher_volume_spawn_input(assembly,m,done;
                        ranks,old_supports,surface_parts=geometry_parts,
                        constructor_inputs)
                    spawn_input===nothing && break
                end
                _parallel_volume_ready(spawn_input,vwork) || continue
                # `periodic_pairs` is still being written by this merge loop —
                # hand the worker a snapshot taken on this task so it never
                # races a later entry.  Eligible Regions have no periodic
                # boundary anyway, so the copy is never consulted.
                # A predicted RNG consumer additionally needs its serial
                # ledger position: the completed predecessor drawer's
                # reported exit, else the jump-ahead over the predicted
                # draws — `nothing` leaves it to the in-phase batch arm.
                entry_seed=options.random_state[]
                if vwork.may_consume
                    seed=_parallel_volume_spawn_seed(volume_forward.prepared,
                        volume_forward.pos_draws,volume_forward.results,vi,
                        entry_seed)
                    seed===nothing && continue
                    entry_seed=seed
                end
                _parallel_volume_spawn!(volume_forward.results,vi,m,vwork,
                    spawn_input,options,assembly,parallel_dirty,geometry_birth;
                    periodic_surface_pairs=periodic_pairs===nothing ? nothing :
                        copy(periodic_pairs),
                    optimize_pyramids,precise=volume_forward.precise,
                    size_field,entry_seed)
            end
        end
    end
    return input,scope,retained,volume_forward
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
        constructor.node_uids,constructor.coordinates,caller;embedded_point=constructor.embedded_point,
        embedded_curve=constructor.embedded_curve,embedded_surface=constructor.embedded_surface,
        embedded_constraints=constructor.embedded_constraints)
    Model._model_native_hybrid_volume_certificate(mesh,source_tags,constructor,caller;
        embedded_pool=production.embedded_pool)
    Model._model_native_hybrid_embedded_point_certificate(mesh,source_tags,constructor,caller)
    length(source_tags)==nnodes(mesh) && allunique(source_tags) && all(!iszero,source_tags) ||
        throw(ArgumentError("$caller: hybrid Volume lacks complete current UID rows"))
    actual=input.mesh;data=actual.entity_data
    rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_tags))
    old_rows=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    closure=Model._model_mesh_projection_closure(m,3,tag,caller)
    # A joint embedded pool can serve a neighbor-declared entity to this
    # Region (shared interface sub-facets bind its nodes on both sides); its
    # served record node owners join the closure, and a pooled entity's lower
    # cells outside the served set drop instead of demanding local incidence.
    pool=production.embedded_pool
    pool_entities=Set{Tuple{Int,Int32}}()
    if pool!==nothing
        served=Set(get(pool.region_records,tag,UInt64[]))
        for entry in pool.records
            entry.uid in served && push!(closure,entry.node_owner)
        end
        # A pooled entity contributes its complete lower product once: where
        # the pool serves it to another Region its cells drop here, and so do
        # the cells of its lower boundary wire (an embedded sheet's Curves and
        # Points have no served vertices in this Region either).
        pending=NTuple{2,Int}[entity for entity in pool.inventory]
        while !isempty(pending)
            owner=pop!(pending);owner[1]<1 && continue
            key=(owner[1],Int32(owner[2]))
            key in pool_entities && continue
            push!(pool_entities,key)
            haskey(data.entities,key) ||
                throw(ArgumentError("$caller: hybrid embedded pool declares a foreign lower entity"))
            for child in data.entities[key].boundaries
                push!(pending,(owner[1]-1,Int(abs(child))))
            end
        end
    end
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
                    owner[1]==0 || owner in pool_entities || throw(ArgumentError(
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

function _higher_volume_input_classification(assembly,m,mesh::Union{Mesh,MixedMesh},input,tag,
                                            source_tags,caller;
                                            boundary_writeback=nothing)
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
    # A boundary-subdivided mesh carries its write-back receipt: the
    # certificate audits the subdivided product cells it will publish, and
    # the recovered Steiner vertices take their written-back owners.
    geometry=(surface=mesh,native_node_tags=source_tags)
    constraints=boundary_writeback===nothing ?
        Model._model_native_volume_input_constraints(m,tag,geometry,input,caller) :
        boundary_writeback.constraints
    Model._model_certify_native_volume_constraints(mesh,constraints,tag,caller)
    boundary_writeback===nothing || for (row,owner) in boundary_writeback.nodes
        row in eachindex(node_entities) || throw(ArgumentError(
            "$caller: Volume[$tag] write-back vertex $row is outside the mesh"))
        node_entities[row]=owner
        parameters[row]=get(boundary_writeback.node_params,row,nothing)
    end
    edges=Dict{NTuple{2,Int32},Tuple{Int,Int32}}()
    faces=Dict{NTuple{3,Int32},Tuple{Int,Int32}}()
    quads=Dict{NTuple{4,Int32},Tuple{Int,Int32}}()
    edge_conflicts=Set{NTuple{2,Int32}}();face_conflicts=Set{NTuple{3,Int32}}();quad_conflicts=Set{NTuple{4,Int32}}()
    owners=Dict{Int,Vector{Int32}}(1=>Int32[],2=>Int32[])
    if mesh isa Mesh
        tetrahedra=fill(Int32(tag),ntets(mesh))
        _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,quads,quad_conflicts,4,mesh.tets,tetrahedra)
        owners[4]=tetrahedra
    else
        # Recombined hex-dominant part: every dim-3 block contributes its
        # edges to the support map, owned by this Region.
        for block in mesh.blocks
            block isa ElementBlock && msh_dimension(Int(block.msh))==3 || continue
            block_owners=fill(Int32(tag),size(block.nodes,2))
            _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,
                quads,quad_conflicts,Int(block.msh),block.nodes,block_owners)
            owners[Int(block.msh)]=block_owners
        end
    end
    for owner in sort!(collect(closure))
        owner[1]<3 || continue
        wb_cells=boundary_writeback===nothing ? nothing :
            owner[1]==2 ? get(boundary_writeback.surfaces,Int(owner[2]),nothing) :
            owner[1]==1 ? get(boundary_writeback.curves,Int(owner[2]),nothing) : nothing
        if wb_cells!==nothing
            for (msh,conn,_) in wb_cells
                _api_support_cells!(edges,edge_conflicts,faces,face_conflicts,
                    quads,quad_conflicts,msh,reshape(conn,length(conn),1),
                    Int32[owner[2]])
            end
            continue
        end
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
    return _mixed_classification(mesh,primary,sort!(collect(closure)),node_entities,boundaries,owners;
        edge_entities=edges,face_entities=faces,quad_entities=quads,node_parametric=parameters)
end

# A written-back surface part keeps every kept column's coordinate, physical
# and identity trace; the recovered children join as fresh columns so the
# geometry carrier still certifies an explicit cell history. All cell keys
# are sorted Volume-mesh row tuples — the write-back's cell space.
function _higher_volume_wb_rebuild_part(assembly,m,wb,mapping,source_tags,
                                        tag,mesh,class,part_mapping,added_cells,
                                        caller;geometry_cells=nothing)
    cells_wb=get(wb.surfaces,tag,nothing)
    cells_wb===nothing && return nothing
    vrow_of_uid=Dict(uid=>Int32(row) for (row,uid) in enumerate(source_tags))
    pos_to_row=Dict(part_mapping[row]=>Int32(row) for row in eachindex(part_mapping))
    extra=Int32[]
    part_row=Dict{Int32,Int32}()
    for (_,conn,_) in cells_wb,row in conn
        haskey(part_row,row) && continue
        part_row[row]=get(pos_to_row,mapping[Int(row)],Int32(0))
        part_row[row]==0 && push!(extra,row)
    end
    new_rows=Dict{Int32,Int32}()
    for (offset,row) in enumerate(extra)
        part_row[row]=Int32(nnodes(mesh)+offset)
        new_rows[Int32(nnodes(mesh)+offset)]=row
    end
    function part_conn(conn)
        mapped=Int32[part_row[row] for row in conn]
        all(row->row>0,mapped) || throw(ArgumentError(
            "$caller: written-back Surface[$tag] cell loses a produced vertex"))
        return mapped
    end
    node_entities=Vector{Tuple{Int,Int32}}(undef,nnodes(mesh)+length(extra))
    for row in eachindex(part_mapping)
        node_entities[row]=class.node_entities[row]
    end
    new_mapping=Vector{Int32}(undef,nnodes(mesh)+length(extra))
    new_mapping[1:length(part_mapping)].=part_mapping
    coords=Matrix{Float64}(undef,3,nnodes(mesh)+length(extra))
    coords[:,1:nnodes(mesh)].=mesh.coords
    for (part_index,vrow) in new_rows
        position=mapping[Int(vrow)]
        node=assembly.nodes[position]
        new_mapping[part_index]=position
        node_entities[part_index]=get(wb.nodes,vrow,(2,Int32(tag)))
        for axis in 1:3
            coords[axis,part_index]=node.point[axis]
        end
    end
    blocks=_cache_native_blocks(mesh)
    old_index=Dict{Tuple,Tuple{Int,Int32}}()
    for (bi,(msh,cells,_)) in enumerate(blocks)
        msh_dimension(msh)==2 || continue
        for column in axes(cells,2)
            key=Tuple(sort!(Int32[get(vrow_of_uid,assembly.nodes[
                part_mapping[cells[slot,column]]].tag,Int32(0))
                for slot in axes(cells,1)]))
            any(iszero,key) || (old_index[key]=(bi,Int32(column)))
        end
    end
    child_physical=Dict{Tuple,Int32}()
    for (_,conn,_,physical) in get(wb.added,(2,Int32(tag)),
                                   Tuple{Int,Vector{Int32},UInt64,Int32}[])
        child_physical[Tuple(sort(conn))]=physical
    end
    emitted_keys=Set(Tuple(sort!(copy(conn))) for (_,conn,_) in cells_wb)
    old_traces=geometry_cells===nothing ? nothing : get(geometry_cells,mesh,nothing)
    if mesh isa Mesh && all(entry->entry[1]==2,cells_wb)
        tris=Matrix{Int32}(undef,3,0)
        tri_tag=Int32[]
        tri_trace=Union{Nothing,_Dim01Cell}[]
        for (msh,conn,_) in cells_wb
            key=Tuple(sort!(copy(conn)))
            previous=get(old_index,key,nothing)
            tris=hcat(tris,reshape(part_conn(conn),3,1))
            if previous===nothing
                push!(tri_tag,get(child_physical,key,_dim01_physical(m,(2,tag))))
                push!(tri_trace,get(added_cells,((2,Int32(tag)),key),nothing))
            else
                bi,column=previous
                push!(tri_tag,blocks[bi][3][column])
                push!(tri_trace,old_traces===nothing ? nothing :
                      old_traces[bi][column])
            end
        end
        new_mesh=Mesh(coords;segs=mesh.segs,tris=tris,
            seg_tag=mesh.seg_tag,tri_tag=tri_tag)
        new_class=_higher_geometry_part_class(new_mesh,2,tag,node_entities)
        if geometry_cells!==nothing
            delete!(geometry_cells,mesh)
            geometry_cells[new_mesh]=[Union{Nothing,_Dim01Cell}[],tri_trace,
                                      Union{Nothing,_Dim01Cell}[]]
        end
        return (tag,new_mesh,new_class,new_mapping)
    end
    if mesh isa Mesh
        # A written-back quadrangle promotes the Triangle product to a
        # mixed part: cells emit per-msh blocks in receipt order, upstream's
        # `gf->triangles`/`gf->quadrangles` split.
        new_blocks=ElementBlock[];new_traces=Vector{Vector{Union{Nothing,_Dim01Cell}}}()
        group_of=Dict{Int,Int}()
        for (msh,conn,_) in cells_wb
            key=Tuple(sort!(copy(conn)))
            previous=get(old_index,key,nothing)
            gi=get(group_of,msh,0)
            if gi==0
                push!(new_blocks,ElementBlock(msh,Matrix{Int32}(undef,length(conn),0),Int32[]))
                push!(new_traces,Union{Nothing,_Dim01Cell}[])
                gi=length(new_blocks);group_of[msh]=gi
            end
            block=new_blocks[gi]
            block=ElementBlock(msh,hcat(block.nodes,part_conn(conn)),block.tags)
            tagout=previous===nothing ?
                get(child_physical,key,_dim01_physical(m,(2,tag))) :
                blocks[previous[1]][3][previous[2]]
            new_blocks[gi]=ElementBlock(msh,block.nodes,push!(copy(block.tags),tagout))
            push!(new_traces[gi],previous===nothing ?
                get(added_cells,((2,Int32(tag)),key),nothing) :
                (old_traces===nothing ? nothing : old_traces[previous[1]][previous[2]]))
        end
        new_mesh=MixedMesh(coords,new_blocks)
        new_class=_higher_geometry_part_class(new_mesh,2,tag,node_entities)
        if geometry_cells!==nothing
            delete!(geometry_cells,mesh)
            geometry_cells[new_mesh]=new_traces
        end
        return (tag,new_mesh,new_class,new_mapping)
    end
    # MixedMesh part: keep each block's emitted columns; children append as
    # per-msh blocks in receipt order.
    new_blocks=ElementBlock[]
    new_traces=Vector{Vector{Union{Nothing,_Dim01Cell}}}()
    for (bi,(msh,cells,tags)) in enumerate(blocks)
        msh_dimension(msh)==2 || continue
        keep=Int32[]
        for column in axes(cells,2)
            key=Tuple(sort!(Int32[get(vrow_of_uid,assembly.nodes[
                part_mapping[cells[slot,column]]].tag,Int32(0))
                for slot in axes(cells,1)]))
            any(iszero,key) && continue
            key in emitted_keys && push!(keep,Int32(column))
        end
        isempty(keep) && continue
        push!(new_blocks,ElementBlock(msh,cells[:,keep],tags[keep]))
        push!(new_traces,Union{Nothing,_Dim01Cell}[
            old_traces===nothing ? nothing : old_traces[bi][column]
            for column in keep])
    end
    child_groups=Dict{Int,Tuple{Vector{Vector{Int32}},Vector{Int32},
                                Vector{Union{Nothing,_Dim01Cell}}}}()
    child_order=Int[]
    for (msh,conn,_) in cells_wb
        key=Tuple(sort!(copy(conn)))
        haskey(old_index,key) && continue
        msh in (2,3) || throw(ArgumentError(
            "$caller: written-back Surface[$tag] cannot emit non-Triangle/Quadrangle children"))
        group=get(child_groups,msh,nothing)
        if group===nothing
            group=(Vector{Int32}[],Int32[],Union{Nothing,_Dim01Cell}[])
            child_groups[msh]=group;push!(child_order,msh)
        end
        push!(group[1],part_conn(conn))
        push!(group[2],get(child_physical,key,_dim01_physical(m,(2,tag))))
        push!(group[3],get(added_cells,((2,Int32(tag)),key),nothing))
    end
    for msh in child_order
        group=child_groups[msh]
        isempty(group[1]) && continue
        nodes=Matrix{Int32}(undef,length(group[1][1]),length(group[1]))
        for (column,conn) in enumerate(group[1]),slot in eachindex(conn)
            nodes[slot,column]=conn[slot]
        end
        push!(new_blocks,ElementBlock(msh,nodes,group[2]))
        push!(new_traces,group[3])
    end
    new_mesh=MixedMesh(coords,new_blocks;physical_names=mesh.physical_names,
        elementary_entities=mesh.elementary_entities)
    new_class=_higher_geometry_part_class(new_mesh,2,tag,node_entities)
    if geometry_cells!==nothing
        delete!(geometry_cells,mesh)
        geometry_cells[new_mesh]=new_traces
    end
    return (tag,new_mesh,new_class,new_mapping)
end

# The same write-back in constructor label space. A registered Surface
# constructor input is bound against its actual cell UIDs; when the Region's
# recovery replaces boundary cells, the constructor receipt must drop the
# consumed labels and adopt the realized children, or the publication check
# can no longer match its stored labels to live events.
function _higher_volume_wb_rebuild_constructor(assembly,m,wb,mapping,source_tags,
        tag,constructor,added_cells,caller)
    cells_wb=get(wb.surfaces,tag,nothing)
    cells_wb===nothing && return nothing
    data=constructor.entity_data
    data===nothing && return nothing
    owner=(2,Int32(tag))
    function conn_uids(conn)
        uids=UInt64[]
        for row in conn
            position=mapping[Int(row)]
            position==0 && return nothing
            push!(uids,assembly.nodes[position].tag)
        end
        return Tuple(sort(uids))
    end
    emitted_keys=Set{Tuple}()
    child_physical=Dict{Tuple,Int32}()
    for (_,conn,_,physical) in get(wb.added,owner,Tuple{Int,Vector{Int32},UInt64,Int32}[])
        child_physical[Tuple(sort(conn))]=physical
    end
    for (_,conn,_) in cells_wb
        key=conn_uids(conn)
        key===nothing || push!(emitted_keys,key)
    end
    row_of_uid=Dict(uid=>Int32(row) for (row,uid) in enumerate(data.external_node_tags))
    old_index=Set{Tuple}()
    blocks=ElementBlock[];labels=Vector{UInt64}[];owners=Vector{Int32}[]
    for (bi,block) in enumerate(constructor.blocks)
        blabels=data.external_element_tags[bi];bowners=data.block_entities[bi]
        if msh_dimension(block.msh)!=2
            push!(blocks,block);push!(labels,blabels);push!(owners,bowners)
            continue
        end
        keep=Int32[]
        for column in axes(block.nodes,2)
            key=Tuple(sort(data.external_node_tags[block.nodes[:,column]]))
            push!(old_index,key)
            key in emitted_keys && push!(keep,Int32(column))
        end
        isempty(keep) && continue
        push!(blocks,ElementBlock(block.msh,block.nodes[:,keep],block.tags[keep]))
        push!(labels,blabels[keep]);push!(owners,bowners[keep])
    end
    # Recovered children join as fresh columns, like the product rebuild's
    # trailing Triangle block; fresh boundary vertices append to the
    # constructor's node receipt so every child's support binds to an actual
    # UID row.
    extra_uids=UInt64[];extra_entities=Tuple{Int,Int32}[];
        extra_params=Union{Nothing,Vector{Float64}}[];extra_coords=NTuple{3,Float64}[]
    child_nodes=Dict{Int,Vector{Vector{Int32}}}();child_labels=Dict{Int,Vector{UInt64}}()
    child_tags=Dict{Int,Vector{Int32}}();child_order=Int[]
    for (msh,conn,_) in cells_wb
        key=conn_uids(conn)
        key===nothing && continue
        key in old_index && continue
        msh in (2,3) || throw(ArgumentError(
            "$caller: written-back Surface[$tag] constructor cannot emit non-Triangle/Quadrangle children"))
        ckey=Tuple(sort(conn))
        mapped=Int32[]
        for (slot,vrow) in enumerate(conn)
            uid=assembly.nodes[mapping[Int(vrow)]].tag
            row=get(row_of_uid,uid,Int32(0))
            if row==0
                push!(extra_uids,uid)
                push!(extra_entities,get(wb.nodes,Int32(vrow),owner))
                push!(extra_params,nothing)
                push!(extra_coords,assembly.nodes[mapping[Int(vrow)]].point)
                row=Int32(length(data.external_node_tags)+length(extra_uids))
                row_of_uid[uid]=row
            end
            push!(mapped,row)
        end
        haskey(child_nodes,msh) || (child_nodes[msh]=Vector{Int32}[];
            child_labels[msh]=UInt64[];child_tags[msh]=Int32[];
            push!(child_order,msh))
        push!(child_nodes[msh],mapped)
        event=get(added_cells,(owner,ckey),nothing)
        event===nothing && throw(ArgumentError(
            "$caller: written-back Surface[$tag] constructor child lacks its actual event"))
        push!(child_labels[msh],event.tag)
        push!(child_tags[msh],get(child_physical,ckey,_dim01_physical(m,(2,tag))))
    end
    for msh in child_order
        columns=child_nodes[msh]
        isempty(columns) && continue
        nodes=Matrix{Int32}(undef,length(columns[1]),length(columns))
        for (column,conn) in enumerate(columns),slot in eachindex(conn)
            nodes[slot,column]=conn[slot]
        end
        push!(blocks,ElementBlock(msh,nodes,child_tags[msh]))
        push!(labels,child_labels[msh])
        push!(owners,fill(Int32(tag),length(child_labels[msh])))
    end
    isempty(extra_uids) || allunique(extra_uids) || throw(ArgumentError(
        "$caller: written-back Surface[$tag] constructor repeats a fresh node UID"))
    node_tags=vcat(collect(data.external_node_tags),extra_uids)
    node_entities=vcat(collect(data.node_entities),extra_entities)
    node_parametric=vcat(collect(data.node_parametric),extra_params)
    coords=isempty(extra_coords) ? copy(constructor.coords) :
        hcat(constructor.coords,reduce(hcat,(collect(p) for p in extra_coords);init=Matrix{Float64}(undef,3,0)))
    metadata=Elements._native_mixed_entity_data(data.entities;
        node_entities=node_entities,node_parametric=node_parametric,
        external_node_tags=node_tags,block_entities=owners,
        external_element_tags=labels)
    return MixedMesh(coords,blocks;physical_names=constructor.physical_names,
        entity_data=metadata,
        elementary_entities=constructor.elementary_entities===nothing ? nothing : owners)
end

# Commit the volume-boundary subdivision write-back onto the emitted
# Surface/Curve products: drop every consumed source cell, emit the realized
# children under their owning entity, then re-bind the geometry-view parts so
# the published products match the recovered boundary exactly.
function _higher_volume_boundary_writeback!(assembly,m,wb,mapping,source_tags,caller;
                                           geometry_parts=nothing,geometry_cells=nothing,
                                           surface_constructor_inputs=nothing)
    source_tags===nothing && throw(ArgumentError(
        "$caller: native Volume boundary write-back lacks its row UID binding"))
    if !isempty(wb.replaced)
        lookup=Dict{Tuple{Tuple{Int,Int32},Tuple},Int32}()
        for (index,cell) in enumerate(assembly.cells)
            haskey(wb.replaced,cell.owner) || continue
            lookup[(cell.owner,Tuple(sort(cell.nodes)))]=Int32(index)
        end
        drop=Int32[]
        for (owner,conns) in wb.replaced,conn in conns
            mapped=Tuple(sort!(Int32[mapping[Int32(row)] for row in conn]))
            index=get(lookup,(owner,mapped),Int32(0))
            index==0 && throw(ArgumentError(
                "$caller: Volume boundary write-back lost an actual entity " *
                "$owner source cell"))
            push!(drop,index)
        end
        deleteat!(assembly.cells,sort!(drop))
    end
    added_cells=Dict{Tuple{Tuple{Int,Int32},Tuple},_Dim01Cell}()
    for owner in sort!(collect(keys(wb.added)))
        for (msh,conn,label,physical) in wb.added[owner]
            nodes=Int32[mapping[Int32(row)] for row in conn]
            _dim01_add_cell!(assembly,msh,owner,nodes,label;
                source_tag=nothing,physical)
            added_cells[(owner,Tuple(sort(conn)))]=last(assembly.cells)
        end
    end
    if geometry_parts!==nothing
        for (index,(tag,mesh,class,part_mapping)) in enumerate(geometry_parts)
            rebuilt=_higher_volume_wb_rebuild_part(assembly,m,wb,mapping,source_tags,
                tag,mesh,class,part_mapping,added_cells,caller;geometry_cells)
            rebuilt===nothing || (geometry_parts[index]=rebuilt)
        end
    end
    surface_constructor_inputs===nothing && return nothing
    for tag in sort!(collect(keys(wb.surfaces)))
        owner=(2,Int32(tag))
        constructor=get(surface_constructor_inputs,owner,nothing)
        constructor===nothing && continue
        rebuilt=_higher_volume_wb_rebuild_constructor(assembly,m,wb,mapping,source_tags,
            tag,constructor,added_cells,caller)
        rebuilt===nothing || (surface_constructor_inputs[owner]=rebuilt)
    end
    return nothing
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

# Volume-cell publication for a swept part that carries no ordinary
# production — QuadTri and NoNew Regions keep their own planners but emit
# cells during `meshGRegionExtruded` in Region order like recombined
# producers. Synthesize their constructor events from the part's block
# order, the same ordering the unsplit `_higher_append_part!` cell pass
# uses, then replay the certified publication loop.
function _higher_emitted_volume_cells!(assembly,m,mesh,class,mapping,geometry_cells)
    events=Tuple{Int,Int}[]
    for (msh,cells,_) in _cache_native_blocks(mesh)
        msh_dimension(msh)==3 || continue
        for column in axes(cells,2)
            push!(events,(msh,column))
        end
    end
    _higher_ordinary_volume_cells!(assembly,m,
        (mesh=mesh,cell_events=events),class,mapping,geometry_cells)
    return nothing
end

function _higher_ordinary_volume_stage!(assembly,m,input,extruded,meshed,lateral_plans,
        discarded,done,emitted;source_nodes,geometry_parts,geometry_cells,native_orphan_births,
        scope=nothing,geometry_birth=nothing)
    # All prospective cell and vertex incidence is certified before allocating
    # this stage. The detached operation later certifies the final lower cells.
    # Swept Regions without an ordinary production — QuadTri and NoNew
    # volumes — certify through their own planner instead, matching the
    # per-Region classification the unstaged loop applies.
    prepared=Dict{Int,Any}()
    for tag in sort!(copy(meshed))
        production=get(extruded.productions,tag,nothing)
        if production===nothing
            mesh=extruded.parts[tag]
            alternate=nothing
            if geometry_parts!==nothing && scope!==nothing
                completed=get(scope.volumes,tag,nothing)
                completed===nothing || (alternate=Model._extrude_nonew_b4_geometry(
                    m,tag,completed,
                    uid->_higher_geometry_source_point(assembly,geometry_birth,uid),
                    assembly.caller))
            end
            swept=_higher_nonew_identity_classification(
                assembly,m,mesh,input,scope,3,tag,assembly.caller)
            class=swept===nothing ? _classify_cached_mesh(
                m,_classification_skeleton(
                    mesh,get(m.meshing.reverse,(3,tag),false)),3,tag,mesh;
                _extrude_scope=scope) : swept.class
            class===nothing && throw(ArgumentError(
                "$(assembly.caller): staged non-ordinary Volume[$tag] lacks ownership metadata"))
            prepared[tag]=(mesh=mesh,class=class,alternate=alternate,
                node_tags=swept===nothing ? nothing : swept.node_tags,
                retained=swept===nothing ? nothing : swept.retained)
            continue
        end
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
        item=prepared[tag]
        if hasproperty(item,:production)
            certified=item.certificate
            mappings[tag]=_higher_append_part!(assembly,m,item.production.mesh,certified.class,3,tag,discarded;
                emitted,source_nodes,source_row_tags=certified.node_tags,
                source_row_retained=certified.retained,source_rows_existing_only=true,
                native_orphan_births,emit_cells=false)
        else
            mappings[tag]=_higher_append_part!(assembly,m,item.mesh,item.class,3,tag,discarded;
                emitted,source_nodes,source_row_tags=item.node_tags,
                source_row_retained=item.retained,
                source_rows_existing_only=item.node_tags!==nothing,
                native_orphan_births,
                emit_cells=false)
        end
    end
    # meshGRegionExtruded visits Regions in tag order. Recombined cell births
    # survive; non-recombined initial cells are deleted by subdivision.
    for tag in sort!(copy(meshed))
        item=prepared[tag]
        if hasproperty(item,:production)
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
        else
            # QuadTri/NoNew cells are born during the same Region visit and
            # carry no consumed provisional events.
            _higher_emitted_volume_cells!(assembly,m,item.mesh,item.class,
                mappings[tag],geometry_cells)
        end
    end
    # SubdivideExtrudedMesh then allocates final cells for each pending Region.
    for tag in sort!(copy(meshed))
        hasproperty(prepared[tag],:production) || continue
        Model._extrude_gate(m,3,tag).recombine && continue
        item=prepared[tag]
        _higher_ordinary_volume_cells!(assembly,m,item.production,item.certificate.class,
            mappings[tag],geometry_cells)
    end
    for tag in sort!(copy(meshed))
        item=prepared[tag]
        if !hasproperty(item,:production) && item.alternate!==nothing
            geometry_cells===nothing ||
                (geometry_cells[item.alternate]=geometry_cells[item.mesh])
            if geometry_birth!==nothing
                for (row,position) in enumerate(mappings[tag])
                    node=assembly.nodes[position]
                    point=(item.alternate.coords[1,row],
                           item.alternate.coords[2,row],
                           item.alternate.coords[3,row])
                    geometry_birth.nodes[node.tag]=
                        _Dim01Node(node.tag,point,node.owner,node.parameter)
                end
            end
        end
        if geometry_parts!==nothing
            if hasproperty(item,:production)
                push!(geometry_parts,
                    (tag,item.production.mesh,item.certificate.class,mappings[tag]))
            else
                geometry_mesh,geometry_class,geometry_mapping=
                    _higher_volume_geometry_order(
                        item.alternate===nothing ? item.mesh : item.alternate,
                        item.class,mappings[tag],m,tag,input,scope;
                        geometry_cells,caller=assembly.caller)
                push!(geometry_parts,
                    (tag,geometry_mesh,geometry_class,geometry_mapping))
            end
        end
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
                     any(labels->event.tag in labels,authority.mesh.entity_data.external_element_tags))) ||
                    throw(ArgumentError(
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
                               surface_constructor_inputs=Dict{Tuple{Int,Int32},MixedMesh}(),
                               volume_forward=nothing)
    have=Set(cell.owner for cell in assembly.cells)
    # Native extrusion is scheduled before meshGRegion's OnlyEmpty gate.
    meshed=_higher_volume_mesh_tags(m,entities,options,have)
    for tag in meshed
        push!(discarded,(3,Int32(tag)))
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
        # Every scheduled Region must live in the shared sweep: ordinary
        # producers carry a certified production, while the other swept
        # planners (QuadTri, NoNew) publish their part through the same
        # tag-ordered node and cell passes below. A non-swept Region still
        # has no swept birth plan and keeps the existing guard.
        all(tag->haskey(extruded.productions,tag) ||
            haskey(extruded.parts,tag),meshed) || throw(ArgumentError(
            "$(assembly.caller): ordinary subdivision combined with another Region planner lacks a complete birth plan"))
        _higher_ordinary_volume_stage!(assembly,m,input,extruded,meshed,lateral_plans,
            discarded,done,emitted;source_nodes,geometry_parts,geometry_cells,native_orphan_births,
            scope,geometry_birth)
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
            Model._model_native_hybrid_volume_input(m,tag,input,assembly.caller;
                _deferred_embeds=true)!==nothing,meshed)
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
    volume_tags=staged ? Int[] : sort!(meshed)
    # PAR01 — `General.NumThreads` detaches runs of independent Regions: each
    # eligible unit produces its receipt on a transaction copy under
    # `Threads.@threads`, and the unchanged serial merge below replays or
    # re-runs each slot in Region order.  `nothing` keeps the exact serial path.
    parallel_plan=_parallel_volume_prepare(m,volume_tags,input,assembly.caller;
        size_field,scope,extruded,hybrid_batch,staged)
    parallel_units=parallel_plan===nothing ? nothing : parallel_plan[1]
    parallel_draws=parallel_plan===nothing ? nothing : parallel_plan[2]
    parallel_results=parallel_units===nothing ? nothing :
        fill!(Vector{Any}(undef,length(volume_tags)),nothing)
    parallel_dirty=Set{Tuple{Int,Int32}}()
    parallel_precise=parallel_units===nothing ? false :
        GmshLibm._win_cli_precision()
    # Cross-phase adoption: the Surface merge loop may already have spawned
    # detached Regions whose full dependency set had landed.  Adopt exactly
    # the slots this phase's real gates (`extruded`/`staged`/hybrid arms,
    # footprint eligibility) still admit, and thread the shared dirty ledger
    # through — it carries every Curve- and Surface-phase mark the spawned
    # workers must be tested against.
    if parallel_units!==nothing && volume_forward!==nothing &&
            volume_forward.tags==volume_tags
        for (vi,work) in enumerate(parallel_units)
            work===nothing && continue
            spawned=volume_forward.results[vi]
            spawned===nothing && continue
            parallel_results[vi]=spawned
        end
        union!(parallel_dirty,volume_forward.dirty)
        parallel_precise=volume_forward.precise
    end
    parallel_last=parallel_units===nothing ? 0 :
        something(findlast(!isnothing,parallel_units),0)
    for (vi,tag) in enumerate(volume_tags)
        work=parallel_units===nothing ? nothing : parallel_units[vi]
        if work!==nothing && parallel_results[vi]===nothing
            vend=vi
            while vend<length(volume_tags) && parallel_units[vend+1]!==nothing
                vend+=1
            end
            _parallel_volume_fill!(parallel_results,vi:vend,m,parallel_units,
                parallel_draws,input,options,assembly,parallel_dirty,
                parallel_precise,
                geometry_birth;periodic_surface_pairs,optimize_pyramids,
                size_field)
        end
        refs=vi<parallel_last ? _parallel_surface_param_refs(m) : nothing
        mesh=extruded===nothing ? nothing : get(extruded.parts,tag,nothing)
        source_rows=nothing;source_row_retained=nothing;source_rows_existing_only=false
        source_cells=nothing;source_cell_max=nothing;alternate=nothing;source_emission=nothing;hybrid_production=nothing
        wb_receipt=nothing
        source_node_base=hybrid_batch_receipt===nothing ? assembly.node_max : hybrid_batch_receipt.node_base
        source_node_max=nothing
        if mesh!==nothing && geometry_parts!==nothing && scope!==nothing
            completed=get(scope.volumes,tag,nothing)
            completed===nothing || (alternate=Model._extrude_nonew_b4_geometry(m,tag,
                completed,uid->_higher_geometry_source_point(assembly,geometry_birth,uid),assembly.caller))
        end
        if mesh===nothing
            production=hybrid_batch===nothing ? (
                work!==nothing && parallel_results[vi]!==nothing ?
                _parallel_volume_product(m,work,parallel_results[vi],
                    parallel_dirty,assembly,options;input,scope,size_field,
                    periodic_surface_pairs,optimize_pyramids,geometry_birth,
                    caller=assembly.caller) :
                mesh_model_volume(m,tag;size_field,
                _extrude_scope=scope,_extrude_working_model=true,_mesh_input=input,
                _mesh1d_options=options,_return_volume_production=true,
                _periodic_surface_pairs=periodic_surface_pairs,
                _optimize_pyramids=optimize_pyramids,
                _mesh_node_high_water=assembly.node_max,_mesh_element_high_water=assembly.cell_max,
                _boundary_regrade=false,
                _geometry_points=_higher_geometry_source_points(assembly,geometry_birth))) :
                hybrid_batch.products[tag]
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
            wb_receipt=production isa Model._ModelVolumeGeometryProduction ?
                production.boundary_writeback : nothing
        end
        native_rows=(mesh isa Mesh || mesh isa MixedMesh) &&
            source_rows!==nothing &&
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
            _higher_volume_input_classification(assembly,m,mesh,input,tag,source_rows,assembly.caller;
                boundary_writeback=wb_receipt)
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
            # Keep their born identities through finish, like ordinary sweep
            # events — swept non-ordinary Regions share that convention
            # (QuadTri preallocates unused interior centroids upstream).
            native_orphan_births=(ordinary===nothing && swept===nothing &&
                source_emission===nothing &&
                (extruded===nothing || !haskey(extruded.parts,tag))) ?
                nothing : native_orphan_births)
        # The consumed Surface products live in `surface_parts` (the emitted
        # geometry vector only collects Volume parts inside this phase).
        wb_receipt===nothing || _higher_volume_boundary_writeback!(
            assembly,m,wb_receipt,mapping,source_rows,assembly.caller;
            geometry_parts=surface_parts,geometry_cells,surface_constructor_inputs)
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
        # State this Region left behind that a later detached produce still
        # reads: stored Curve discretizations (embedded-curve lazy grading and
        # boundary writeback replayed or written during the merge).
        refs===nothing || _parallel_volume_mark_dirty!(parallel_dirty,m,refs)
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
            _extrude_edges=extruded===nothing ? nothing : extruded.edges,_mesh1d_options=options,
            _boundary_regrade=false) :
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
        if authority.mesh.coords[:,rows]!=mesh.coords
            # Upstream `pos.find` welding: a produced boundary row can bind to
            # a vertex another producer emitted first (the curve parts' T/tops
            # derive from an earlier kernel pass than the strip's own columns)
            # — the appended node then carries the welded position, not the
            # producer's own eval. The same `64eps·max(1,|p|)` rule
            # `_higher_append_part!` welds by decides identity; a bound node
            # beyond it is a real mismatch, not a weld.
            delta=abs.(authority.mesh.coords[:,rows].-mesh.coords)
            bad=count(axes(delta,2)) do j
                tolerance=64eps(Float64)*max(1.0,maximum(@view mesh.coords[:,j]))
                maximum(@view delta[:,j])>tolerance
            end
            bad>0 && begin
                detail=join(("$(j):bound=$(authority.mesh.coords[:,rows[j]]) part=$(mesh.coords[:,j])"
                             for j in axes(delta,2) if any(!=(0.0),delta[:,j])),"\n  ")
                throw(ArgumentError(
                    "$(assembly.caller): produced Surface[$tag] coordinates " *
                    "disagree with actual bound rows " *
                    "($bad/$(size(delta,2)) cols, max |delta|=$(maximum(delta)))\n  $detail"))
            end
        end
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
    # MixedMesh parts (folded recombination products, retained quads) carry
    # their per-family cell counts in the blocks — the simplex vectors stay
    # filled for every dimension the part actually emits.
    nseg_count,ntri_count,ntet_count=if mesh isa MixedMesh
        counts=(0,0,0)
        for block in mesh.blocks
            block isa ElementBlock || continue
            d=msh_dimension(block.msh)
            1<=d<=3 || continue
            counts=Base.setindex(counts,counts[d]+size(block.nodes,2),d)
        end
        counts
    else
        (nsegs(mesh),ntris(mesh),ntets(mesh))
    end
    segments=fill(Int32(tag),nseg_count);triangles=fill(Int32(tag),ntri_count)
    tetrahedra=fill(Int32(tag),ntet_count)
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
    # A recombined transfinite product carries quadrilateral blocks; the
    # geometry view keeps them as a MixedMesh rather than folding the part.
    all(block->block isa ElementBlock && block.msh in (2,3),
        product.blocks) || return nothing
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
    geometry=if all(block->block.msh==2,product.blocks)
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
        Mesh(source.coords[:,rows];tris=triangles,tri_tag=physical)
    else
        blocks=ElementBlock[ElementBlock(block.msh,
            Int32[inverse[product.rows[block.nodes[row,index]]]
                  for row in axes(block.nodes,1), index in axes(block.nodes,2)],
            copy(block.tags)) for block in product.blocks]
        MixedMesh(source.coords[:,rows],blocks)
    end
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
                                       geometry_birth=nothing,discarded=nothing,
                                       high_order=nothing,native_orphan_births=nothing)
    parts=Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}[]
    points=Dict{Int,Vector{Int32}}()
    curves=Dict{Int,Vector{_Dim01Cell}}()
    # A produced part can only list vertices that survive into the actual
    # mesh. `_dim01_finish!` drops unreferenced orphans of discarded
    # (remeshed) owners and old interpolation vertices — the same rule
    # filters the owner-orphan rows below once the caller supplies it.
    referenced=nothing
    if discarded!==nothing || high_order!==nothing || native_orphan_births!==nothing
        referenced=falses(length(assembly.nodes))
        for cell in assembly.cells,row in cell.nodes
            referenced[row]=true
        end
    end
    survives(row,node)=
        (referenced===nothing || referenced[row]) ||
        (native_orphan_births!==nothing && node.tag in native_orphan_births) ||
        ((discarded===nothing || !(node.owner in discarded)) &&
         (high_order===nothing || !(row<=length(high_order) && high_order[row])))
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
            # An orphan kept alive by another entity's incidence still
            # lists: it remains an actual vertex owned by this Curve. A
            # dead orphan of a remeshed Curve does not.
            append!(rows,Int32[row for (row,node) in enumerate(assembly.nodes)
                if node.owner==(1,Int32(tag)) && !(Int32(row) in used) &&
                survives(row,node)])
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
    # Upstream `pos.find` keeps a single vertex store: when a produced
    # boundary row welds onto a vertex an earlier producer emitted, the
    # welded position — not this part's own eval — is the canonical
    # coordinate. Rewrite each part's coords to those bound values before
    # the carrier merge, or two producers' ~1-ulp-different copies of the
    # same logical vertex become two carrier rows bound to one actual
    # vertex (the curve parts' T/tops derive from an earlier kernel pass
    # than the strip's own columns). The birth record names the welded
    # value, matching the `expected.point` audit below.
    canonical=Vector{Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,Vector{Int32}}}(undef,length(parts))
    for (index,(tag,mesh,class,mapping)) in enumerate(parts)
        coords=mesh.coords
        needs=false
        for position in eachindex(mapping)
            source=assembly.nodes[mapping[position]]
            expected=geometry_birth===nothing ? source :
                get(geometry_birth.nodes,source.tag,source)
            expected.owner==source.owner || (expected=source)
            coords[:,position]==collect(expected.point) || (needs=true;break)
        end
        if needs
            new_coords=copy(coords)
            for position in eachindex(mapping)
                source=assembly.nodes[mapping[position]]
                expected=geometry_birth===nothing ? source :
                    get(geometry_birth.nodes,source.tag,source)
                expected.owner==source.owner || (expected=source)
                new_coords[:,position]=collect(expected.point)
            end
            new_mesh=if mesh isa Mesh
                Mesh(new_coords;segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                     seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,
                     tet_tag=mesh.tet_tag)
            else
                MixedMesh(new_coords,mesh.blocks;
                    physical_names=mesh.physical_names,
                    elementary_entities=mesh.elementary_entities)
            end
            traces=get(geometry_cells,mesh,nothing)
            traces===nothing || (geometry_cells[new_mesh]=traces)
            canonical[index]=(tag,new_mesh,class,mapping)
        else
            canonical[index]=(tag,mesh,class,mapping)
        end
    end
    parts=canonical
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
            actual_row>0 || throw(ArgumentError(
                "$caller: produced geometry vertex loses its actual identity"))
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
    if dimension==3
        # `Mesh3D` folds the STORED boundary quads of `recombine3D` Regions at
        # fill time (`quadsToTriangles`), so a face meshed before
        # `Recombine Volume`/`Mesh.Recombine3DAll` was set still reaches the
        # Delaunay PLC all-triangle — and stays split for every later read.
        # `cached` is always the frontend's `authority.mesh` (the retained
        # actual), so the folded replacement is its logical successor.
        folded=Model._recombine3d_fold_authority(m,m.meshing.mesh_authority;caller,
            above=UInt64(initial_max_element_tag))
        if folded!==m.meshing.mesh_authority
            m.meshing.mesh_authority=folded
            cached=folded.mesh
            cached_class=_higher_classification_from_data(cached,caller)
        end
    end
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
    # PAR01 cross-phase scheduling state: `surface_forward` plans the Surface
    # units before Mesh1D so the Curve merge loop can detach each eligible
    # Face the moment its boundary products land; `volume_forward` plans the
    # Region units inside `_higher_surface_phase!` so each eligible Region can
    # detach while later Faces still merge.  Both remain speculative on the
    # consuming phase's real gates, which re-verify every adopted slot.
    surface_forward=nothing
    volume_forward=nothing
    if phases.one
        for entity in entities
            entity[1]>=2 && !_dim01_fully_discrete(m,entity) &&
                push!(discarded,(entity[1],Int32(entity[2])))
        end
        filter!(cell->!(cell.owner in discarded),assembly.cells)
        filter!(owner->!(owner in discarded),done)
        first_nodes=_dim01_point_phase!(assembly,m,entities;
            options=options,size_field=size_field)
        if phases.two
            # The keep/mesh scan is phase-stable (Mesh1D only appends dim-0/1
            # cells), so the unit plan the Surface phase will adopt already
            # exists here; `scope` alone stays speculative and the phase's
            # real `_parallel_surface_prepare` keeps the final word on every
            # adopted slot.
            forward_have=Set(cell.owner for cell in assembly.cells)
            forward_meshed,_=_higher_surface_mesh_scan(m,entities,options,
                forward_have,discarded,assembly.caller;dimension)
            forward_units=_higher_surface_units(m,
                _generate_units(m,2,sort!(forward_meshed)),assembly.caller)
            forward_periodic_masters=Set(Int(constraint.master_entity) for
                constraint in Model.model_periodic_constraints(m) if constraint.dim==2)
            forward_ordinary=Set{Int}();forward_unmerged=Set{Int}()
            for volume in keys(m.volumes)
                Model._extrude_nonew_unmerged_b4_incidence(m,volume,assembly.caller) || continue
                push!(forward_unmerged,abs(m.meshing.extrude_sources[(3,volume)][2]))
            end
            for ((dim,entity),link) in m.meshing.extrude_sources
                dim in (2,3) && link[1]==2 || continue
                params=Model._extrude_gate(m,dim,entity)
                params!==nothing && params.quad_to_tri===:none || continue
                push!(forward_ordinary,abs(link[2]))
            end
            forward_prepared=_parallel_surface_prepare(m,forward_units,
                assembly.caller;metric_surface,scope=nothing,
                periodic_masters=forward_periodic_masters,
                ordinary_sources=forward_ordinary,
                unmerged_nonew_sources=forward_unmerged,
                size_field,geometry_birth)
            forward_prepared===nothing || (surface_forward=(
                units=forward_units,prepared=forward_prepared,
                results=fill!(Vector{Any}(undef,length(forward_units)),nothing),
                dirty=Set{Tuple{Int,Int32}}(),
                precise=GmshLibm._win_cli_precision(),
                ranks,old_supports=high_order,surface_parts,
                size_field=size_field))
        end
        union!(discarded,_dim01_curve_phase!(assembly,m,options,entities,first_nodes,high_order;
            done,geometry_birth,curve_births,surface_forward))
    else
        # The legacy Surface/Volume projection omitted its completed Mesh0D
        # products. Recover those once, just as skipped Curve products are
        # recovered below; an actual cache never invents a skipped Point cell.
        legacy_points=old!==nothing && old.mesh===nothing && cached!==nothing
        first_nodes=_dim01_point_phase!(assembly,m,entities;
                                       create_cells=legacy_points,create_nodes=legacy_points,
                                       options=options,size_field=size_field)
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
        _,_,retained_surfaces,volume_forward=_higher_surface_phase!(assembly,m,options,entities,discarded,done;
                               size_field,metric_surface,source_nodes=length(high_order),ranks,old_supports=high_order,
                                geometry_parts=surface_parts,dimension,geometry_cells,
                                geometry_birth,geometry_alternates,curve_births,periodic_pairs,native_orphan_births,native_recombine_events,
                                surface_constructor_inputs,surface_forward,
                                optimize_pyramids=Float64(get(generation_meshing,:optimize_pyramids,0.0)))
    end
    if phases.three
        _higher_volume_phase!(assembly,m,options,entities,discarded,done;
                             size_field,source_nodes=length(high_order),ranks,old_supports=high_order,
                             geometry_parts,surface_parts,geometry_cells,geometry_birth,native_orphan_births,
                             periodic_surface_pairs=periodic_pairs,
                             optimize_pyramids=Float64(get(generation_meshing,:optimize_pyramids,0.0)),
                             retained_surfaces,
                             surface_constructor_inputs,volume_forward)
    end
    if preserve_geometry && geometry_layout===:all
        # GEO's native carrier has the Surface/Volume products followed by
        # Point/Curve geometry. Reuse each actual emitted Surface vector and
        # the staged lower cells' explicit assembly rows; this does not mesh
        # another product or infer a public identity from coordinates.
        _higher_complete_ranks!(ranks,assembly)
        surfaces=Dict{Int,Tuple{Int,Union{Mesh,MixedMesh},_MeshClassification,
            Vector{Int32}}}(tag=>(tag,mesh,class,mapping)
                            for (tag,mesh,class,mapping) in surface_parts)
        retained=_higher_assembly_authority(assembly,m,done;ranks,old_supports=high_order)
        retained=_higher_bind_producer_rows(retained,assembly,surface_parts,old)
        retained_cells=Dict(cell.tag=>cell for cell in assembly.cells)
        for owner in keys(retained.products)
            owner[1]==2 && !haskey(surfaces,Int(owner[2])) || continue
            input=_higher_surface_geometry_input(retained,owner,caller)
            input===nothing && continue
            geometry=input.mesh
            mapping=Int32[assembly.node_positions[tag] for tag in input.node_tags]
            class=if geometry isa MixedMesh
                _mixed_classification(geometry,(2,Int32(owner[2])),
                    [(2,Int32(owner[2]))],input.node_entities,
                    Dict{Tuple{Int,Int32},Vector{Int32}}(),
                    Dict{Int,Vector{Int32}}(Int(block.msh)=>
                        fill(Int32(owner[2]),size(block.nodes,2))
                        for block in geometry.blocks))
            else
                _higher_geometry_part_class(geometry,2,Int(owner[2]),
                    input.node_entities)
            end
            surfaces[Int(owner[2])]=(Int(owner[2]),geometry,class,mapping)
            geometry_cells[geometry]=geometry isa MixedMesh ?
                Vector{Union{Nothing,_Dim01Cell}}[
                    Union{Nothing,_Dim01Cell}[retained_cells[tag] for tag in tags]
                    for tags in input.element_tags] :
                Vector{Union{Nothing,_Dim01Cell}}[Union{Nothing,_Dim01Cell}[],
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
            sort_points,geometry_cells,geometry_birth,discarded,high_order,
            native_orphan_births))
    end
    p2_pass=element_order>1 && any(cell->cell.owner[1]>0,assembly.cells)
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    # Retained cells above the requested order take the general order path —
    # the same dispatch `set_order` uses — so a discrete entity storing mixed
    # P1/P3 cells converts rather than hitting the P2 pass's first-element guard.
    if element_order>2 || (p2_pass && (incomplete ||
            any(cell->msh_spec(cell.msh).order>2,assembly.cells)))
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
    if previous!==nothing && !isempty(previous.surface_constructor_inputs) &&
            hasproperty(plan,:node_tag_map)
        # Retained Surface constructor receipts carry the previous mesh's node
        # and cell UIDs; rebind them through this plan's tag maps and re-validate
        # against the rebound products so a lost identity cannot slip through.
        element_tag_map=hasproperty(plan,:element_tag_map) ? plan.element_tag_map : nothing
        rebound_constructors=Dict{Tuple{Int,Int32},MixedMesh}()
        for (owner,constructor) in previous.surface_constructor_inputs
            haskey(authority.products,owner) || continue
            copy_constructor=deepcopy(constructor)
            data=copy_constructor.entity_data
            data===nothing && continue
            for (row,tag) in enumerate(data.external_node_tags)
                data.external_node_tags[row]=get(plan.node_tag_map,tag,UInt64(0))
            end
            for labels in data.external_element_tags,(row,tag) in enumerate(labels)
                labels[row]=element_tag_map===nothing ? tag :
                    get(element_tag_map,tag,UInt64(0))
            end
            # The receipt is a first-order record: it rebinds through an
            # order pass onto the product's corner skeleton — an order
            # elevated product keeps it, while a product that was genuinely
            # reshaped (a remesh or import losing a bound row or cell
            # signature) drops it like the pre-rebind behavior did instead of
            # throwing during the strict bind.
            bound=false
            try
                actual=Model._model_mesh_product(
                    plan.model.meshing.mesh_authority,2,Int(owner[2]);
                    caller="higher plan authority Surface constructor rebind",
                    producer_order=true,first_order=true)
                if actual!==nothing
                    Model._model_surface_constructor_binding(actual,copy_constructor,
                        "higher plan authority Surface constructor rebind")
                    bound=true
                end
            catch err
                err isa ArgumentError || rethrow()
            end
            bound && (rebound_constructors[owner]=copy_constructor)
        end
        isempty(rebound_constructors) || (plan.model.meshing.mesh_authority=
            Model._model_mesh_with_surface_constructor_inputs(
                plan.model.meshing.mesh_authority,rebound_constructors;
                caller="higher plan authority Surface constructor rebind"))
    end
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
        # `rebound` keeps its mesh and tag space, so any retained Surface
        # constructor inputs stay bound to the same products.
        plan.model.meshing.mesh_authority=Model._ModelMeshAuthority(rebound.mesh,
            rebound.products,rebound.done,rebound.private_rank,rebound.native_inputs,
            rebound.carrier,rebound.carrier_rows,produced,rebound.carrier_cells,
            rebound.surface_constructor_inputs)
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
                max_nodes,max_cells,generation_meshing=get(kwargs,:generation_meshing,nothing),
                reuse_curve_params=get(kwargs,:reuse_curve_params,false))
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
