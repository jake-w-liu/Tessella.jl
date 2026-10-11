struct _ModelNativeHybridEmbeddedCurveInput
    model::GeoModel
    region::Int
    curve::Int
    endpoints::Tuple{Int,Int}
    input::_ModelMeshAuthority
    records::Vector{NamedTuple}
    cells::Vector{NamedTuple}
end

function _model_native_hybrid_embedded_curve_input(m,region,input,caller)
    embeds=get(m.embeds,(3,region),NTuple{2,Int}[])
    length(m.volumes)==1 && length(embeds)==1 && only(embeds)[1]==1 || return nothing
    curve=only(embeds)[2]
    haskey(m.curves,curve) || throw(ArgumentError("$caller: embedded Curve has no declared geometry"))
    product=_model_mesh_product(input,1,curve;caller)
    product isa MixedMesh && product.entity_data!==nothing || throw(ArgumentError(
        "$caller: embedded Curve lacks its actual lower discretization"))
    data=product.entity_data;cells=NamedTuple[];used=Set{UInt64}();seen=Set{UInt64}()
    supports=Set{Tuple{UInt64,UInt64}}()
    for (block,labels,owners) in zip(product.blocks,data.external_element_tags,data.block_entities)
        block.msh==1 && size(block.nodes,1)==2 && length(labels)==length(owners)==size(block.nodes,2) ||
            throw(ArgumentError("$caller: embedded Curve needs actual linear Line2 cells"))
        for (column,uid) in enumerate(labels)
            nodes=ntuple(slot->data.external_node_tags[block.nodes[slot,column]],2)
            uid>0 && !(uid in seen) && owners[column]==curve && all(!iszero,nodes) && allunique(nodes) ||
                throw(ArgumentError("$caller: embedded Curve repeats or changes a current cell UID/owner"))
            support=minmax(nodes...)
            support in supports && throw(ArgumentError("$caller: embedded Curve repeats a current segment"))
            push!(seen,uid);push!(supports,support);union!(used,nodes)
            push!(cells,(owner=(1,Int32(curve)),uid=uid,nodes=nodes))
        end
    end
    isempty(cells) && throw(ArgumentError("$caller: embedded Curve has no actual segment cells"))
    endpoint_uids=UInt64[];endpoints=m.curves[curve]
    for point in endpoints
        actual=_model_mesh_product(input,0,point;caller)
        actual isa MixedMesh && actual.entity_data!==nothing && length(actual.blocks)==1 &&
            only(actual.blocks).msh==15 && size(only(actual.blocks).nodes)==(1,1) || throw(ArgumentError(
            "$caller: embedded Curve endpoint lacks its actual Point-cell witness"))
        point_row=only(only(actual.blocks).nodes)
        actual.entity_data.node_entities[point_row]==(0,Int32(point)) &&
            length(only(actual.entity_data.external_element_tags))==1 &&
            only(only(actual.entity_data.external_element_tags))>0 || throw(ArgumentError(
            "$caller: embedded Curve endpoint changes its actual Point-cell UID/owner"))
        push!(endpoint_uids,actual.entity_data.external_node_tags[point_row])
    end
    degrees=Dict(uid=>0 for uid in used)
    for cell in cells,uid in cell.nodes;degrees[uid]+=1;end
    expected_ends=first(endpoint_uids)==last(endpoint_uids) ? Set{UInt64}() : Set(endpoint_uids)
    Set(uid for (uid,degree) in degrees if degree==1)==expected_ends &&
        all(degree->degree in (1,2),values(degrees)) && all(uid->uid in used,endpoint_uids) || throw(ArgumentError(
        "$caller: embedded Curve cells do not form its declared endpoint chain"))
    reached=Set([first(endpoint_uids)]);changed=true
    while changed
        changed=false
        for cell in cells
            any(uid->uid in reached,cell.nodes) || continue
            for uid in cell.nodes
                uid in reached || (push!(reached,uid);changed=true)
            end
        end
    end
    reached==used || throw(ArgumentError("$caller: embedded Curve has disconnected actual segments"))
    positions=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    records=NamedTuple[]
    for uid in sort!(collect(used))
        row=positions[uid];owner=data.node_entities[row]
        (owner==(1,Int32(curve)) || owner[1]==0 && Int(owner[2]) in endpoints) ||
            throw(ArgumentError("$caller: embedded Curve node loses its actual lower owner"))
        push!(records,(uid=uid,coordinate=Tuple(product.coords[:,row]),node_owner=owner))
    end
    return _ModelNativeHybridEmbeddedCurveInput(m,region,curve,endpoints,input,records,cells)
end

function _model_native_hybrid_curve_source_certificate(witness,input,uids,coordinates,caller)
    witness===nothing && return nothing
    witness isa _ModelNativeHybridEmbeddedCurveInput && witness.input===input &&
        get(witness.model.embeds,(3,witness.region),NTuple{2,Int}[])==[(1,witness.curve)] &&
        get(witness.model.curves,witness.curve,nothing)==witness.endpoints || throw(ArgumentError(
        "$caller: embedded Curve has a foreign or changed constructor"))
    fresh=_model_native_hybrid_embedded_curve_input(witness.model,witness.region,input,caller)
    fresh!==nothing && fresh.records==witness.records && fresh.cells==witness.cells || throw(ArgumentError(
        "$caller: embedded Curve changes its actual current UID/owner/cell witness"))
    positions=Dict(uid=>row for (row,uid) in enumerate(uids))
    for record in witness.records
        row=get(positions,record.uid,0)
        row>0 && isequal(Tuple(coordinates[:,row]),record.coordinate) || throw(ArgumentError(
            "$caller: embedded Curve UID/XYZ is missing or changed in the native Source"))
    end
    return nothing
end

function _model_native_hybrid_curve_containment_certificate(constructor,caller)
    witness=constructor.embedded_curve;witness===nothing && return nothing
    _model_native_hybrid_curve_source_certificate(witness,constructor.input,
        constructor.node_uids,constructor.coordinates,caller)
    positions=Dict(uid=>row for (row,uid) in enumerate(constructor.node_uids))
    panels=vcat(constructor.triangles,constructor.quads)
    boundary=Set(uid for panel in panels for uid in panel.nodes)
    coordinate(uid)=Tuple(constructor.coordinates[:,positions[uid]])
    for record in witness.records
        record.uid in boundary && throw(ArgumentError("$caller: embedded Curve must be strictly inside the shell"))
    end
    for panel in panels
        triangles=length(panel.nodes)==3 ? (panel.nodes,) :
            ((panel.nodes[1],panel.nodes[2],panel.nodes[3]),(panel.nodes[1],panel.nodes[3],panel.nodes[4]))
        for triangle in triangles
            a,b,c=map(coordinate,triangle)
            sides=[orient3(a,b,c,coordinate(uid)) for uid in boundary]
            all(isfinite,sides) || throw(ArgumentError("$caller: embedded Curve shell halfspaces are not finite"))
            signs=Set(sign(side) for side in sides if side!=0)
            length(signs)==1 || throw(ArgumentError("$caller: embedded Curve needs a certified convex shell"))
            for record in witness.records
                side=orient3(a,b,c,record.coordinate)
                isfinite(side) && side!=0 && sign(side)==only(signs) || throw(ArgumentError(
                    "$caller: embedded Curve is outside or on the shell boundary"))
            end
        end
    end
    return nothing
end

function _model_native_hybrid_curve_segments(constructor,uids,caller)
    witness=constructor.embedded_curve;witness===nothing && return NTuple{2,Int32}[]
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(uids))
    all(record->haskey(positions,record.uid),witness.records) ||
        throw(ArgumentError("$caller: embedded Curve constraint omits an actual node UID"))
    return NTuple{2,Int32}[ntuple(slot->positions[cell.nodes[slot]],2) for cell in witness.cells]
end

function _model_native_hybrid_embedded_curve_certificate(mesh,uids,constructor,caller)
    witness=constructor.embedded_curve;witness===nothing && return nothing
    _model_native_hybrid_curve_source_certificate(witness,constructor.input,uids,mesh.coords,caller)
    _model_native_hybrid_curve_containment_certificate(constructor,caller)
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(uids))
    for record in witness.records
        _model_native_hybrid_interior_node_star_certificate(mesh,positions[record.uid],caller,"embedded Curve vertex")
    end
    for segment in _model_native_hybrid_curve_segments(constructor,uids,caller)
        links=Set{Tuple{Int32,Int32}}();degrees=Dict{Int32,Int}()
        for block in mesh.blocks,column in axes(block.nodes,2)
            nodes=block.nodes[:,column]
            all(node->node in nodes,segment) || continue
            block.msh==4 || throw(ArgumentError("$caller: embedded Curve segment has a foreign typed cell role"))
            opposite=Int32[node for node in nodes if !(node in segment)]
            length(opposite)==2 || throw(ArgumentError("$caller: embedded Curve segment has degenerate Tet incidence"))
            edge=minmax(opposite...)
            edge in links && throw(ArgumentError("$caller: embedded Curve repeats a Tet edge-link panel"))
            push!(links,edge)
            for node in edge;degrees[node]=get(degrees,node,0)+1;end
        end
        length(links)>=3 && all(==(2),values(degrees)) || throw(ArgumentError(
            "$caller: embedded Curve segment lacks its closed interior Tet edge link"))
        reached=Set([first(first(links))]);changed=true
        while changed
            changed=false
            for edge in links
                any(node->node in reached,edge) || continue
                for node in edge
                    node in reached || (push!(reached,node);changed=true)
                end
            end
        end
        reached==Set(keys(degrees)) || throw(ArgumentError("$caller: embedded Curve has disconnected Tet edge links"))
    end
    return nothing
end
