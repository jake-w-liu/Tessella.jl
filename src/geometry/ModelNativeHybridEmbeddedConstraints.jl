# A complete actual Region inventory, with one lower witness per declared entity.
struct _ModelNativeHybridEmbeddedConstraintsInput
    model::GeoModel
    region::Int
    input::_ModelMeshAuthority
    inventory::Vector{NTuple{2,Int}}
    points::Vector{_ModelNativeHybridEmbeddedPointInput}
    curves::Vector{_ModelNativeHybridEmbeddedCurveInput}
    records::Vector{NamedTuple}
end

function _model_native_hybrid_embedded_constraints_input(m,region,input,caller)
    inventory=copy(get(m.embeds,(3,region),NTuple{2,Int}[]))
    length(m.volumes)==1 && length(inventory)>1 && all(owner->owner[1] in (0,1),inventory) ||
        return nothing
    allunique(inventory) || throw(ArgumentError("$caller: hybrid constraint inventory repeats an entity"))
    points=_ModelNativeHybridEmbeddedPointInput[]
    curves=_ModelNativeHybridEmbeddedCurveInput[]
    for (dim,tag) in inventory
        if dim==0
            push!(points,_model_native_hybrid_embedded_point_input(m,region,input,caller;point=tag))
        else
            push!(curves,_model_native_hybrid_embedded_curve_input(m,region,input,caller;curve=tag))
        end
    end
    records=Dict{UInt64,NamedTuple}()
    function retain(record)
        current=(uid=record.uid,coordinate=record.coordinate,node_owner=record.node_owner)
        prior=get(records,current.uid,nothing)
        prior===nothing || isequal(prior,current) || throw(ArgumentError(
            "$caller: hybrid constraint inventory changes a shared actual node UID/owner/XYZ"))
        records[current.uid]=current
    end
    for point in points;retain(point.record);end
    for curve in curves,record in curve.records;retain(record);end
    return _ModelNativeHybridEmbeddedConstraintsInput(m,region,input,inventory,points,curves,
        NamedTuple[records[uid] for uid in sort!(collect(keys(records)))])
end

function _model_native_hybrid_constraints_source_certificate(receipt,input,uids,coordinates,caller)
    receipt===nothing && return nothing
    receipt isa _ModelNativeHybridEmbeddedConstraintsInput && receipt.input===input &&
        length(receipt.model.volumes)==1 && haskey(receipt.model.volumes,receipt.region) &&
        get(receipt.model.embeds,(3,receipt.region),NTuple{2,Int}[])==receipt.inventory &&
        length(receipt.inventory)>1 && allunique(receipt.inventory) &&
        all(owner->owner[1] in (0,1),receipt.inventory) || throw(ArgumentError(
        "$caller: hybrid constraints lack their complete current Region inventory"))
    owners=vcat([(0,Int(point.record.owner[2])) for point in receipt.points],
                [(1,curve.curve) for curve in receipt.curves])
    length(owners)==length(receipt.inventory) && allunique(owners) &&
        Set(owners)==Set(receipt.inventory) || throw(ArgumentError(
        "$caller: hybrid constraints omit or repeat an actual declared entity"))
    for point in receipt.points
        point.model===receipt.model && point.region==receipt.region && point.inventory==receipt.inventory ||
            throw(ArgumentError("$caller: hybrid Point constraint has a foreign inventory"))
        _model_native_hybrid_point_source_certificate(point,input,uids,coordinates,caller)
    end
    for curve in receipt.curves
        curve.model===receipt.model && curve.region==receipt.region && curve.inventory==receipt.inventory ||
            throw(ArgumentError("$caller: hybrid Curve constraint has a foreign inventory"))
        _model_native_hybrid_curve_source_certificate(curve,input,uids,coordinates,caller)
    end
    fresh=_model_native_hybrid_embedded_constraints_input(receipt.model,receipt.region,input,caller)
    fresh!==nothing && fresh.records==receipt.records || throw(ArgumentError(
        "$caller: hybrid constraint inventory changes its complete actual Source rows"))
    return nothing
end

function _model_native_hybrid_constraints_containment_certificate(constructor,caller)
    receipt=constructor.embedded_constraints;receipt===nothing && return nothing
    _model_native_hybrid_constraints_source_certificate(receipt,constructor.input,
        constructor.node_uids,constructor.coordinates,caller)
    for point in receipt.points
        _model_native_hybrid_point_containment_certificate(merge(constructor,(embedded_point=point,)),caller)
    end
    for curve in receipt.curves
        _model_native_hybrid_curve_containment_certificate(merge(constructor,(embedded_curve=curve,)),caller)
    end
    return nothing
end

function _model_native_hybrid_constraints_segments(constructor,uids,caller)
    receipt=constructor.embedded_constraints;receipt===nothing && return NTuple{2,Int32}[]
    segments=NTuple{2,Int32}[];seen=Set{NTuple{2,Int32}}()
    for curve in receipt.curves
        for edge in _model_native_hybrid_curve_segments(merge(constructor,(embedded_curve=curve,)),uids,caller)
            key=minmax(edge...)
            # The geometric constraint is one edge even when actual lower
            # cells from different declared Curves reference the same support.
            key in seen && continue
            push!(seen,key);push!(segments,edge)
        end
    end
    return segments
end

function _model_native_hybrid_embedded_constraints_certificate(mesh,uids,constructor,caller)
    receipt=constructor.embedded_constraints;receipt===nothing && return nothing
    _model_native_hybrid_constraints_source_certificate(receipt,constructor.input,uids,mesh.coords,caller)
    for point in receipt.points
        _model_native_hybrid_embedded_point_certificate(mesh,uids,
            merge(constructor,(embedded_point=point,)),caller)
    end
    for curve in receipt.curves
        _model_native_hybrid_embedded_curve_certificate(mesh,uids,
            merge(constructor,(embedded_curve=curve,)),caller)
    end
    return nothing
end
