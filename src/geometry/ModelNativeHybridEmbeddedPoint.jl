# A selected Point is an actual lower constructor, not an unreferenced row.
struct _ModelNativeHybridEmbeddedPointInput
    model::GeoModel
    region::Int
    input::_ModelMeshAuthority
    seeds::_ModelNativeHybridSizeSeeds
    record::NamedTuple
end

function _model_native_hybrid_embedded_point_input(m,region,input,caller)
    embeds=get(m.embeds,(3,region),NTuple{2,Int}[])
    isempty(embeds) && return nothing
    # Multi-Region and other embedded dimensions retain their existing guards.
    length(m.volumes)==1 && length(embeds)==1 && only(embeds)[1]==0 || return nothing
    seeds=_model_native_hybrid_global_size_seeds(m,input,caller)
    records=[record for record in seeds.records if record.owner==only(embeds)]
    length(records)==1 || throw(ArgumentError(
        "$caller: embedded Point lacks its unique actual Point-cell witness"))
    return _ModelNativeHybridEmbeddedPointInput(m,region,input,seeds,only(records))
end

function _model_native_hybrid_point_source_certificate(witness,input,uids,coordinates,caller)
    witness===nothing && return nothing
    witness isa _ModelNativeHybridEmbeddedPointInput && witness.input===input &&
        get(witness.model.embeds,(3,witness.region),NTuple{2,Int}[])==[witness.record.owner] ||
        throw(ArgumentError("$caller: embedded Point has a foreign or changed constructor"))
    record=witness.record
    record in _model_native_hybrid_seed_records(witness.seeds) &&
        record.node_owner==(0,Int32(record.owner[2])) || throw(ArgumentError(
        "$caller: embedded Point loses its exact current Point-cell/UID/owner binding"))
    row=findfirst(==(record.uid),uids)
    row!==nothing && isequal(Tuple(coordinates[:,row]),record.coordinate) ||
        throw(ArgumentError("$caller: embedded Point UID is missing or changed in the native Source"))
    return nothing
end

function _model_native_hybrid_point_containment_certificate(constructor,caller)
    witness=constructor.embedded_point;witness===nothing && return nothing
    _model_native_hybrid_point_source_certificate(witness,constructor.input,
        constructor.node_uids,constructor.coordinates,caller)
    positions=Dict(uid=>row for (row,uid) in enumerate(constructor.node_uids))
    panels=vcat(constructor.triangles,constructor.quads)
    boundary=Set(uid for panel in panels for uid in panel.nodes)
    witness.record.uid in boundary && throw(ArgumentError(
        "$caller: embedded Point must be strictly inside the hybrid shell"))
    point=witness.record.coordinate
    coordinate(uid)=Tuple(constructor.coordinates[:,positions[uid]])
    for panel in panels
        triangles=length(panel.nodes)==3 ? (panel.nodes,) :
            ((panel.nodes[1],panel.nodes[2],panel.nodes[3]),
             (panel.nodes[1],panel.nodes[3],panel.nodes[4]))
        for triangle in triangles
            a,b,c=map(coordinate,triangle)
            sides=[orient3(a,b,c,coordinate(uid)) for uid in boundary]
            all(isfinite,sides) || throw(ArgumentError(
                "$caller: embedded Point shell halfspaces are not finite"))
            signs=Set(sign(side) for side in sides if side!=0)
            length(signs)==1 || throw(ArgumentError(
                "$caller: embedded Point requires a certified convex hybrid shell"))
            side=orient3(a,b,c,point)
            isfinite(side) && side!=0 && sign(side)==only(signs) || throw(ArgumentError(
                "$caller: embedded Point is outside or on the hybrid shell boundary"))
        end
    end
    return nothing
end

function _model_native_hybrid_embedded_point_certificate(mesh,uids,constructor,caller)
    witness=constructor.embedded_point;witness===nothing && return nothing
    _model_native_hybrid_point_source_certificate(witness,constructor.input,uids,mesh.coords,caller)
    _model_native_hybrid_point_containment_certificate(constructor,caller)
    point=Int32(findfirst(==(witness.record.uid),uids))
    links=Set{NTuple{3,Int32}}()
    for block in mesh.blocks,column in axes(block.nodes,2)
        nodes=block.nodes[:,column];point in nodes || continue
        block.msh==4 || throw(ArgumentError(
            "$caller: embedded Point has a foreign hybrid cell role"))
        face=Tuple(sort!(Int32[node for node in nodes if node!=point]))
        face in links && throw(ArgumentError("$caller: embedded Point repeats a Tet star face"))
        push!(links,face)
    end
    isempty(links) && throw(ArgumentError(
        "$caller: embedded Point UID is unreferenced by the generated volume cells"))
    edges=Dict{Tuple{Int32,Int32},Vector{Int}}();vertices=Set{Int32}();faces=collect(links)
    for (index,face) in enumerate(faces)
        union!(vertices,face)
        for edge in ((face[1],face[2]),(face[1],face[3]),(face[2],face[3]))
            push!(get!(()->Int[],edges,edge),index)
        end
    end
    all(indices->length(indices)==2,values(edges)) &&
        length(vertices)-length(edges)+length(faces)==2 || throw(ArgumentError(
        "$caller: embedded Point lacks a closed interior Tet star"))
    reached=Set([1]);queue=[1];cursor=1
    while cursor<=length(queue)
        current=queue[cursor];cursor+=1
        for indices in values(edges)
            current in indices || continue
            for neighbor in indices
                neighbor in reached || (push!(reached,neighbor);push!(queue,neighbor))
            end
        end
    end
    length(reached)==length(faces) || throw(ArgumentError(
        "$caller: embedded Point has disconnected Tet star components"))
    return nothing
end
