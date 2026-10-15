# `gmsh.model.occ` shape transforms — `OCC_Internals::_transform` (the shared
# translate/rotate/dilate/symmetry/affineTransform path) and `copy` (the
# `_multiBind` recursive rebind). Upstream transforms each bound entity's
# TopoDS subgraph, unbinds the old shapes ownership-aware, and rebinds the
# transformed graph: the top-level shape keeps its tag while every
# re-created sub-entity takes a fresh `getMaxTag+1` tag. Public entities
# churn on the next `synchronize` — the unbound tags flow through
# `kernel.removed` before the new shapes import.

# DFS discovery order over the union subgraph of `roots` — every bound
# sub-shape reachable from the inputs, deduplicated like upstream's single
# transformed compound.
function _model_occ_closure(kernel:: _ModelNativeOCCKernel, roots)
    order=UInt64[]; seen=Set{UInt64}()
    stack=UInt64[reverse(collect(roots))...;]
    while !isempty(stack)
        id=pop!(stack)
        id in seen && continue
        push!(seen,id); push!(order,id)
        children=kernel.shapes[id].children
        for i in length(children):-1:1
            push!(stack,children[i])
        end
    end
    return order
end

# Squared slant of a cone record — its pcurve `v` coordinate is slant length.
@inline _occ_cone_slant(height,dr)=sqrt(height*height+dr*dr)

# Ratio the surface's `v` coordinate stretches by. Planar frames absorb the
# stretch into their stored (un-normalized) basis; spheres and tori
# parametrize in angles; cylinders and cones carry physical length in `v`.
function _model_occ_pcurve_vscale(old,new)
    haskey(old,:occ) || return 1.0
    old.occ===:cylinder && return new.height/old.height
    old.occ===:cone || return 1.0
    sl=_occ_cone_slant(old.height,old.r2-old.r1)
    sl2=_occ_cone_slant(new.height,new.r2-new.r1)
    return sl==0.0 ? 1.0 : sl2/sl
end

function _model_occ_scale_pcurve(pc,vscale)
    pc===nothing && return nothing
    haskey(pc,:lin2d) || return pc
    p0,dir=pc.lin2d
    return (lin2d=((p0[1],vscale*p0[2]),(dir[1],vscale*dir[2])),)
end

function _model_occ_scale_pcurves(pcurves,vscale)
    vscale==1.0 && return pcurves
    result=Dict{Int,NamedTuple}()
    for (tag,e) in pcurves
        result[tag]=(fwd=_model_occ_scale_pcurve(e.fwd,vscale),
                     rev=_model_occ_scale_pcurve(e.rev,vscale))
    end
    return result
end

# Affine image of a standalone NURBS/bezier curve object: control poles map
# exactly under any affine map (weights and knots carry over).
function _model_occ_transform_nurbs_curve(g::_ModelOCCNURBSCurve,t,caller)
    controls=NTuple{3,Float64}[_finite_result(_occ_trsf_apply(t,p),caller)
                               for p in g.curve.controls]
    curve=NURBS.NURBSCurve(g.curve.degree,g.curve.knots,controls,g.curve.weights)
    return _ModelOCCNURBSCurve(g.family,curve)
end

# Compact `primitive` record on a solid — the same family-limited update the
# public `_plan_volume_transform` applies to `m.cylinders`/`m.spheres`/`m.cones`.
function _model_occ_transform_primitive(record,t,caller,what)
    if haskey(record,:radius) && haskey(record,:axis)
        tr=_transform_axis_encoding(t,record.center,record.axis,caller,what)
        return haskey(record,:r1) ?
            (center=tr.center,axis=tr.axis,r1=record.r1*tr.perp,
             r2=record.r2*tr.perp,height=tr.height) :
            (center=tr.center,axis=tr.axis,radius=record.radius*tr.perp,
             height=tr.height)
    elseif haskey(record,:radius)
        s2=_linear_isotropy(t.linear)
        s2===nothing && throw(ArgumentError(
            "$caller: transform is not representable on $what — " *
            "it requires an isotropic linear part"))
        return (center=_finite_result(_occ_trsf_apply(t,record.center),caller),
                radius=record.radius*sqrt(s2))
    end
    return record
end

# Sector record on a partial-sweep solid — the axis-encoded fields
# transform like the matching `primitive` encoding plus the azimuth
# frame (`X`,`Y`) the measure code reads for sector moments; `angle`,
# `lat1`, `lat2` measure against that frame and carry over verbatim.
function _model_occ_transform_sector(rec,t::_AffineTransform,caller,what)
    if rec.occ===:sphere
        s2=_linear_isotropy(t.linear)
        s2===nothing && throw(ArgumentError(
            "$caller: transform is not representable on $what — " *
            "it requires an isotropic linear part"))
        axis=_occ_linear_apply(t,rec.axis)
        X=_occ_linear_apply(t,rec.X); Y=_occ_linear_apply(t,rec.Y)
        return merge(rec,
            (center=_finite_result(_occ_trsf_apply(t,rec.center),caller),
             radius=rec.radius*sqrt(s2),
             axis=_finite_result(axis./sqrt(_dot(axis,axis)),caller),
             X=_finite_result(X./sqrt(_dot(X,X)),caller),
             Y=_finite_result(Y./sqrt(_dot(Y,Y)),caller)))
    end
    tr=_transform_axis_encoding(t,rec.center,rec.axis,caller,what)
    X=_occ_linear_apply(t,rec.X)
    X=_finite_result(X./sqrt(_dot(X,X)),caller)
    Y=_occ_linear_apply(t,rec.Y)
    Y=_finite_result(Y./sqrt(_dot(Y,Y)),caller)
    return rec.occ===:cylinder ?
        merge(rec,(center=tr.center,axis=tr.axis,X=X,Y=Y,
                   radius=rec.radius*tr.perp,height=tr.height)) :
        merge(rec,(center=tr.center,axis=tr.axis,X=X,Y=Y,
                   r1=rec.r1*tr.perp,r2=rec.r2*tr.perp,height=tr.height))
end

# New geometry tuple for one shape under transform `t`. Wire/shell `uses`
# children are remapped separately by the caller; face `pcurves` stay keyed
# by the pre-transform curve tags until the realized tag map exists.
function _model_occ_transform_geometry(m::GeoModel,shape,t,caller,what)
    d=shape.dimension; g=shape.geometry
    d==0 && return merge(g,(point=_finite_result(
        _occ_trsf_apply(t,g.point),caller),))
    if d==1
        rec=g.geometry
        new=if rec.occ===:degenerate
            # The bound face record embedded in the degenerate edge moves
            # with the surface transform; its pcurves are re-keyed below.
            merge(rec,(face=_transform_occ_surface_record(
                rec.face,t,caller,what),))
        elseif rec.occ===:copied_curve
            merge(rec,(curve=_model_occ_transform_nurbs_curve(
                rec.curve,t,caller),))
        else
            _transform_occ_curve_record(m,0,rec,t,caller,what)
        end
        # Embedded face records (degenerate poles, surface_iso seams) carry
        # the same tag-keyed pcurves as bound faces — scale their `v`
        # coordinate like the dim-2 branch.
        if haskey(new,:face) && haskey(new.face,:pcurves)
            vscale=_model_occ_pcurve_vscale(rec.face,new.face)
            new=merge(new,(face=merge(new.face,
                (pcurves=_model_occ_scale_pcurves(
                    new.face.pcurves,vscale),)),))
        end
        return merge(g,(geometry=new,))
    elseif d==2
        rec=g.geometry
        new=_transform_occ_surface_record(rec,t,caller,what)
        if haskey(new,:pcurves)
            new=merge(new,(pcurves=_model_occ_scale_pcurves(
                new.pcurves,_model_occ_pcurve_vscale(rec,new)),))
        end
        return merge(g,(geometry=new,))
    elseif d==3
        sector=get(g,:sector,nothing)
        sector!==nothing && return merge(g,
            (sector=_model_occ_transform_sector(sector,t,caller,what),))
        primitive=get(g,:primitive,nothing)
        primitive===nothing && return g
        return merge(g,(primitive=(occ=primitive.occ,
            record=_model_occ_transform_primitive(
                primitive.record,t,caller,what)),))
    end
    return g
end

# Remap `(child_id, reverse)` use lists and plain child references to the new
# shape identities.
function _model_occ_remap_uses(geometry,remap)
    haskey(geometry,:uses) || return geometry
    uses=Tuple{UInt64,Bool}[(remap[id],rev) for (id,rev) in geometry.uses]
    return merge(geometry,(uses=uses,))
end

# Recursive bind of a new shape tree: the top takes `tag`, each not-yet-bound
# descendant takes the next fresh tag of its dimension — upstream's
# `_bind(shape, tag, true)` exploration order. Shared new shapes already
# bound by an earlier input keep their realized tag.
function _model_occ_bind_tree!(kernel:: _ModelNativeOCCKernel,shapes,id,dim,
                               tag,ctx)
    haskey(kernel.shape_tags,id) && return nothing
    _model_occ_bind_shape!(kernel,dim,tag,id,shapes[id])
    if dim==1
        old=ctx.inverse[id]
        haskey(ctx.oldtags,old) && (ctx.tagmap[ctx.oldtags[old]]=tag)
    end
    for cid in shapes[id].children
        haskey(kernel.shape_tags,cid) && continue
        cdim=shapes[cid].dimension
        ctag=_model_occ_next_tag(kernel,cdim,-1,ctx.caller)
        _model_occ_bind_tree!(kernel,shapes,cid,cdim,ctag,ctx)
    end
    return nothing
end

@inline function _model_occ_rekey_dict(pcurves,tagmap)
    rekeyed=Dict{Int,NamedTuple}()
    for (oldtag,e) in pcurves
        rekeyed[get(tagmap,oldtag,oldtag)]=e
    end
    return rekeyed
end

# Re-key every transformed shape's pcurve dictionaries — bound faces carry
# `rec.pcurves`, embedded faces inside degenerate/surface_iso curve records
# carry `rec.face.pcurves` — from the pre-transform curve tags to the
# realized fresh tags.
function _model_occ_rekey_pcurves!(kernel:: _ModelNativeOCCKernel,newshapes,
                                   tagmap)
    for (id,shape) in newshapes
        if shape.dimension==2
            rec=shape.geometry.geometry
            haskey(rec,:pcurves) || continue
            kernel.shapes[id]=_ModelOCCShape(2,
                (kind=shape.geometry.kind,
                 geometry=merge(rec,(pcurves=_model_occ_rekey_dict(
                     rec.pcurves,tagmap),))),
                shape.children)
        elseif shape.dimension==1
            rec=shape.geometry.geometry
            (haskey(rec,:face) && haskey(rec.face,:pcurves)) || continue
            kernel.shapes[id]=_ModelOCCShape(1,
                (kind=shape.geometry.kind,
                 geometry=merge(rec,(face=merge(rec.face,
                     (pcurves=_model_occ_rekey_dict(
                         rec.face.pcurves,tagmap),)),))),
                shape.children)
        end
    end
    return nothing
end

function _model_occ_transform!(m::GeoModel,dim_tags,t::_AffineTransform,
                               caller::AbstractString)
    kernel=m.occ_kernel
    entries=Tuple{Int,Int,UInt64}[]
    seen=Set{Tuple{Int,Int}}()
    for entry in dim_tags
        d=_model_occ_dimension(first(entry),caller)
        tg=_model_occ_requested_tag(last(entry),caller)
        (d,tg) in seen && throw(ArgumentError(
            "$caller: OCC entity of dimension $d with tag $tg is listed twice"))
        push!(seen,(d,tg))
        id=get(kernel.bindings,(d,tg),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $d with tag $tg"))
        push!(entries,(d,tg,id))
    end
    union=_model_occ_closure(kernel,[id for (_,_,id) in entries])
    remap=Dict{UInt64,UInt64}(zip(
        union,_model_occ_reserve_shapes(kernel,length(union),caller)))
    inverse=Dict{UInt64,UInt64}(v=>k for (k,v) in remap)
    oldtags=Dict{UInt64,Int}(id=>kernel.shape_tags[id][2]
        for id in union if haskey(kernel.shape_tags,id))
    # Worst-case fresh allocations per dimension — every new sub-shape may
    # take its own tag before any of them are bound.
    needed=Dict{Int,Int}()
    for id in union
        dim=kernel.shapes[id].dimension
        needed[dim]=get(needed,dim,0)+1
    end
    for (dim,count) in needed
        kernel.maxima[dim+3]<=typemax(Int32)-count || throw(ArgumentError(
            "$caller: OCC entity tag allocation exceeds Int32"))
    end
    # Every fallible record transform runs before the first bind/unbind —
    # representability failures leave the kernel untouched.
    newshapes=Dict{UInt64,_ModelOCCShape}()
    for id in union
        shape=kernel.shapes[id]
        label=_model_periodic_entity_label(shape.dimension)
        what="OCC $label[$(get(oldtags,id,0))]"
        geometry=_model_occ_remap_uses(
            _model_occ_transform_geometry(m,shape,t,caller,what),remap)
        newshapes[remap[id]]=_ModelOCCShape(shape.dimension,geometry,
            UInt64[remap[c] for c in shape.children])
    end
    # Interleaved unbind/bind, in input order: the ownership-protected
    # recursive unbind mirrors `_unbind(shape, tag, true)`; the recursive
    # rebind mirrors `_bind(out, dim, tag, true)` with fresh sub-tags.
    ctx=(inverse=inverse,oldtags=oldtags,tagmap=Dict{Int,Int}(),caller=caller)
    for (d,tg,id) in entries
        _model_occ_unbind!(kernel,id,true)
        _model_occ_bind_tree!(kernel,newshapes,remap[id],d,tg,ctx)
    end
    _model_occ_rekey_pcurves!(kernel,newshapes,ctx.tagmap)
    kernel.next_shape=max(kernel.next_shape,
                          maximum(collect(values(remap)))+UInt64(1))
    _model_occ_collect_shapes!(kernel)
    return nothing
end

# `OCC_Internals::copy` — a deep copy of the input subgraph rebound entirely
# under fresh tags. `_multiBind(result, -1, out, true, true)` binds only the
# highest present dimension: bare wires contribute their edges and bare
# shells their faces, exactly the TopExp level upstream reports.
function _model_occ_copy!(m::GeoModel,dim_tags,caller::AbstractString)
    kernel=m.occ_kernel
    units=Tuple{Int,UInt64}[]
    for entry in dim_tags
        d=_model_occ_dimension(first(entry),caller)
        tg=_model_occ_requested_tag(last(entry),caller)
        id=get(kernel.bindings,(d,tg),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $d with tag $tg"))
        if d>=0
            push!(units,(d,id))
        elseif d==-1
            append!(units,[(1,c) for c in kernel.shapes[id].children])
        else
            append!(units,[(2,f) for (f,_) in kernel.shapes[id].geometry.uses])
        end
    end
    isempty(units) && return Tuple{Int,Int}[]
    top=maximum(u[1] for u in units)
    units=[u for u in units if u[1]==top]
    union=_model_occ_closure(kernel,[u[2] for u in units])
    remap=Dict{UInt64,UInt64}(zip(
        union,_model_occ_reserve_shapes(kernel,length(union),caller)))
    inverse=Dict{UInt64,UInt64}(v=>k for (k,v) in remap)
    oldtags=Dict{UInt64,Int}(id=>kernel.shape_tags[id][2]
        for id in union if haskey(kernel.shape_tags,id))
    newshapes=Dict{UInt64,_ModelOCCShape}()
    for id in union
        shape=kernel.shapes[id]
        geometry=_model_occ_remap_uses(shape.geometry,remap)
        newshapes[remap[id]]=_ModelOCCShape(shape.dimension,geometry,
            UInt64[remap[c] for c in shape.children])
    end
    ctx=(inverse=inverse,oldtags=oldtags,tagmap=Dict{Int,Int}(),caller=caller)
    out=Tuple{Int,Int}[]
    for (d,id) in units
        nid=remap[id]
        # A copy already bound by an earlier unit's recursion (deduplicated
        # shared sub-shapes, or the same entity listed twice) is not reported
        # again — `_multiBind` only pushes unbound (`!exists`) shapes.
        haskey(kernel.shape_tags,nid) && continue
        tag=_model_occ_next_tag(kernel,d,-1,caller)
        _model_occ_bind_tree!(kernel,newshapes,nid,d,tag,ctx)
        push!(out,(d,tag))
    end
    _model_occ_rekey_pcurves!(kernel,newshapes,ctx.tagmap)
    kernel.next_shape=max(kernel.next_shape,maximum(collect(values(remap)))+
        UInt64(1))
    return out
end

# Public endpoint entry points — each builds its `_AffineTransform` through
# the shared `_affine_*` constructors so `occ.*` and `Transform` agree on the
# gp_Trsf/gp_GTrsf arithmetic.
function _model_occ_translate!(m::GeoModel,dim_tags,dx,dy,dz)
    caller="OCC.translate"
    t=_affine_translation((dx,dy,dz),caller)
    _model_occ_transform!(m,dim_tags,t,caller)
    return nothing
end

function _model_occ_rotate!(m::GeoModel,dim_tags,x,y,z,ax,ay,az,angle)
    caller="OCC.rotate"
    t=_affine_rotation((ax,ay,az),(x,y,z),angle,caller)
    _model_occ_transform!(m,dim_tags,t,caller)
    return nothing
end

function _model_occ_dilate!(m::GeoModel,dim_tags,x,y,z,a,b,c)
    caller="OCC.dilate"
    t=_affine_dilation((x,y,z),(a,b,c),caller)
    _model_occ_transform!(m,dim_tags,t,caller)
    return nothing
end

function _model_occ_mirror!(m::GeoModel,dim_tags,a,b,c,d)
    caller="OCC.mirror"
    t=_affine_symmetry(a,b,c,d,caller)
    _model_occ_transform!(m,dim_tags,t,caller)
    return nothing
end

function _model_occ_affine_transform!(m::GeoModel,dim_tags,matrix)
    caller="OCC.affine_transform"
    matrix isa AbstractVector ||
        throw(ArgumentError("$caller: affine transform must be a vector"))
    entries=Float64[_finite_scalar(v,caller,"affine transform entry")
                    for v in matrix]
    if length(entries)<12
        count=length(entries)
        @warn "$count < 12 entries in affine transform matrix"
        resize!(entries,12)
        entries[count+1:12].=0.0
    end
    linear=(entries[1],entries[2],entries[3],
            entries[5],entries[6],entries[7],
            entries[9],entries[10],entries[11])
    loc=(entries[4],entries[8],entries[12])
    t=_AffineTransform(linear,loc,NTuple{16,Float64}[],
                       (linear=linear,loc=loc))
    _model_occ_transform!(m,dim_tags,t,caller)
    return nothing
end

function _model_occ_copy_entities!(m::GeoModel,dim_tags)
    return _model_occ_copy!(m,dim_tags,"OCC.copy")
end
