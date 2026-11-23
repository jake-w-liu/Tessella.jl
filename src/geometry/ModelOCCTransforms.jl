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

# The shapes under `root` in the order Gmsh's recursive `_bind` meets them
# (`TopExp_Explorer` at each level: a face's edges, an edge's vertices, in
# the kernel's explorer order where it records one); new sub-shapes take
# their tags in this order. => id => walk position.
function _model_occ_walk_rank(kernel::_ModelNativeOCCKernel,root)
    rank=Dict{UInt64,Int}()
    stack=UInt64[root]
    while !isempty(stack)
        id=pop!(stack)
        haskey(rank,id) && continue
        rank[id]=length(rank)+1
        sh=kernel.shapes[id]
        children=sh.children
        if sh.dimension in (1,2) && haskey(kernel.explorer,id)
            # a face's wires take their own rank, its edges come in explorer
            # order rather than through the wires
            if sh.dimension==2
                for w in sh.children
                    haskey(rank,w) || (rank[w]=length(rank)+1)
                end
            end
            children=kernel.explorer[id]
        end
        for i in length(children):-1:1
            push!(stack,children[i])
        end
    end
    return rank
end

# Squared slant of a cone record — its pcurve `v` coordinate is slant length.
@inline _occ_cone_slant(height,dr)=sqrt(height*height+dr*dr)

# Ratio the surface's `v` coordinate stretches by. Planar frames absorb the
# stretch into their stored (un-normalized) basis; spheres and tori
# parametrize in angles; cylinders and cones carry physical length in `v`.
function _model_occ_pcurve_vscale(old,new)
    haskey(old,:occ) || return 1.0
    # Tensor-pole transforms preserve their knot coordinates, including the
    # linear height coordinate of a converted cylinder or cone.
    _model_nurbs_surface(old) && return 1.0
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

# The reparametrization is a parameter-space map — affine pole transforms
# commute with it, so only the homogeneous net needs rebuilding.
function _model_occ_transform_nurbs_curve(g::_ModelOCCReparamCurve,t,caller)
    inner=_model_occ_transform_nurbs_curve(
        _ModelOCCNURBSCurve(g.family,g.curve,g.nets),t,caller)
    return _ModelOCCReparamCurve(inner.family,inner.curve,inner.nets,g.reparam)
end

# Boolean section curve on a tensor surface: the tensor face and the seed
# polyline transform exactly under the affine map, and the pcurve's (u,v)
# coordinates are parameter-space data that carry over verbatim.
function _model_occ_transform_nurbs_curve(g::_ModelOCCSectionCurve,t,caller)
    face=_transform_occ_surface_record(g.face,t,caller,"section curve face")
    controls=NTuple{3,Float64}[_finite_result(_occ_trsf_apply(t,p),caller)
                               for p in g.curve.controls]
    poly=NURBS.NURBSCurve(g.curve.degree,g.curve.knots,controls,g.curve.weights)
    return _ModelOCCSectionCurve(g.family,face,g.pc,poly)
end

# Compact `primitive` record on a solid — the same family-limited update the
# public `_plan_volume_transform` applies to `m.cylinders`/`m.spheres`/`m.cones`.
function _model_occ_transform_primitive(record,t,caller,what)
    if haskey(record,:V)
        # `:revolved` Pappus record — only a mass/centroid carrier
        # (`point + h_c·axis + arm·bisector` is the centroid), so the exact
        # affine image exists for every nonsingular map: the centroid maps
        # with the entity and the volume scales by |det|.
        det=_linear_det(t.linear)
        det==0.0 && throw(ArgumentError(
            "$caller: transform collapses $what"))
        ax=_occ_linear_apply(t,record.axis)
        na=sqrt(_dot(ax,ax))
        na>0 || throw(ArgumentError("$caller: transform collapses $what"))
        new=(point=_finite_result(_occ_trsf_apply(t,record.point),caller),
             axis=_finite_result(ax./na,caller),
             V=record.V*abs(det),h_c=record.h_c*na)
        if haskey(record,:bisector)
            bs=_occ_linear_apply(t,record.bisector)
            nb=sqrt(_dot(bs,bs))
            nb>0 || throw(ArgumentError("$caller: transform collapses $what"))
            new=merge(new,(angle=record.angle,
                bisector=_finite_result(bs./nb,caller),arm=record.arm*nb))
        end
        return new
    end
    if (haskey(record,:radius)||haskey(record,:r1)) && haskey(record,:axis)
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

# GEO-bound volume encodings (`add_box!`/`add_cylinder!`-family records and
# native Boolean snapshots) live on the model, not the kernel — re-derive the
# transformed records so the rebound solid's re-import restores the same
# analytic contract. Boxes travel as the raw encoding: whether the new shell
# still forms an axis-aligned box is a post-materialization check.
function _occ_bound_geo_encodings(m::GeoModel,tag::Int,t,caller)
    enc=NamedTuple()
    haskey(m.cylinders,tag) && (enc=merge(enc,(cylinder=
        _model_occ_transform_primitive(m.cylinders[tag],t,caller,"Cylinder[$tag]"),)))
    haskey(m.spheres,tag) && (enc=merge(enc,(sphere=
        _model_occ_transform_primitive(m.spheres[tag],t,caller,"Sphere[$tag]"),)))
    haskey(m.cones,tag) && (enc=merge(enc,(cone=
        _model_occ_transform_primitive(m.cones[tag],t,caller,"Cone[$tag]"),)))
    haskey(m.box_extents,tag) && (enc=merge(enc,(box=m.box_extents[tag],)))
    haskey(m.booleans,tag) && (enc=merge(enc,(booleans=m.booleans[tag],)))
    haskey(m.boolean_components,tag) &&
        (enc=merge(enc,(boolean_components=m.boolean_components[tag],)))
    if haskey(m.boolean_operands,tag)
        ops=m.boolean_operands[tag]
        enc=merge(enc,(boolean_operands=ops isa Tuple{Mesh,Mesh} ?
            (_transform_mesh_snapshot(ops[1],t,caller),
             _transform_mesh_snapshot(ops[2],t,caller)) :
            merge(ops,(meshes=[_transform_mesh_snapshot(snap,t,caller)
                               for snap in ops.meshes],)),))
    end
    return isempty(fieldnames(typeof(enc))) ? nothing : enc
end

# New geometry tuple for one shape under transform `t`. Wire/shell `uses`
# children are remapped separately by the caller; face `pcurves` stay keyed
# by the pre-transform curve tags until the realized tag map exists.
function _model_occ_transform_geometry(m::GeoModel,shape,t,caller,what;
                                       oldtag::Int=0)
    d=shape.dimension; g=shape.geometry
    d==0 && return merge(g,(point=_finite_result(
        _occ_trsf_apply(t,g.point),caller),))
    if d==1
        rec=g.geometry
        new=_transform_occ_curve_record(m,0,rec,t,caller,what)
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
        if primitive===nothing && oldtag>0
            enc=_occ_bound_geo_encodings(m,oldtag,t,caller)
            enc===nothing && return g
            return merge(g,(geo_encodings=enc,))
        end
        primitive===nothing && return g
        return merge(g,(primitive=(occ=primitive.occ,
            record=_model_occ_transform_primitive(
                primitive.record,t,caller,what)),))
    end
    return g
end

# Does a reflection turn this face's natural normal against the mapped one?
# Sphere, cylinder and cone frames are mapped vector by vector and become
# left-handed; a plane carries its normal along and the torus record negates
# its axis instead, so those two keep their sense.
function _model_occ_reflection_flips(face::_ModelOCCShape)
    rec=face.geometry.geometry
    return rec isa NamedTuple && get(rec,:occ,nothing) in (:sphere,:cylinder,:cone)
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
    elseif dim==0
        old=ctx.inverse[id]
        haskey(ctx.oldtags,old) && (ctx.tagmap0[ctx.oldtags[old]]=tag)
    end
    for cid in shapes[id].children
        haskey(kernel.shape_tags,cid) && continue
        cdim=shapes[cid].dimension
        # Re-bound sub-shapes take their precomputed `getMaxTag+1` tag —
        # upstream emits them in original enumeration order; never-bound
        # shapes take the next fresh tag in walk order.
        ctag=get(ctx.newtags,cid,0)
        ctag==0 && (ctag=_model_occ_next_tag(kernel,cdim,-1,ctx.caller))
        _model_occ_bind_tree!(kernel,shapes,cid,cdim,ctag,ctx)
    end
    return nothing
end

@inline function _model_occ_rekey_dict(pcurves,tagmap,flipped,tol,caller)
    rekeyed=Dict{Int,NamedTuple}()
    for (oldtag,e) in pcurves
        if oldtag in flipped
            # A reflected closed conic stores `p'(t) = T·p(-t)` — its face
            # pcurve columns reparametrize by `s = -t` to stay on the edge's
            # new parameter map (both columns; no use-flag swap).
            e=(fwd=_model_occ_sew_shift_pc(get(e,:fwd,nothing),
                                           -1.0,0.0,tol,caller),
               rev=_model_occ_sew_shift_pc(get(e,:rev,nothing),
                                           -1.0,0.0,tol,caller))
        end
        rekeyed[get(tagmap,oldtag,oldtag)]=e
    end
    return rekeyed
end

# Re-key every transformed shape's pcurve dictionaries — bound faces carry
# `rec.pcurves`, embedded faces inside degenerate/surface_iso curve records
# carry `rec.face.pcurves` — from the pre-transform curve tags to the
# realized fresh tags. `flipped` holds the pre-transform tags of closed
# conic edges whose records a reflection range-negated
# (`_transform_occ_circle`: closed circles keep +T·n with `t0,t1 -> -t1,-t0`).
function _model_occ_rekey_pcurves!(kernel:: _ModelNativeOCCKernel,newshapes,
                                   tagmap;flipped=Set{Int}(),
                                   tagmap0=Dict{Int,Int}(),
                                   tol=_OCC_CONFUSION,caller="OCC.transform")
    for (id,shape) in newshapes
        if shape.dimension==2
            rec=shape.geometry.geometry
            haskey(rec,:pcurves) || haskey(rec,:vertex_uv) ||
                haskey(rec,:plane_rec) || continue
            if haskey(rec,:pcurves)
                rec=merge(rec,(pcurves=_model_occ_rekey_dict(
                    rec.pcurves,tagmap,flipped,tol,caller),))
            end
            # vertex_uv is keyed by vertex tags — the dim-0 map —
            # and a reclassified cap's nested `plane_rec` carries its
            # own edge-keyed pcurve table.
            if haskey(rec,:vertex_uv)
                rec=merge(rec,(vertex_uv=Dict{Int,NTuple{2,Float64}}(
                    get(tagmap0,k,k)=>v for (k,v) in rec.vertex_uv),))
            end
            if haskey(rec,:plane_rec) && haskey(rec.plane_rec,:pcurves)
                rec=merge(rec,(plane_rec=merge(rec.plane_rec,
                    (pcurves=_model_occ_rekey_dict(
                        rec.plane_rec.pcurves,tagmap,flipped,tol,
                        caller),)),))
            end
            kernel.shapes[id]=_ModelOCCShape(2,
                (kind=shape.geometry.kind,geometry=rec),shape.children)
        elseif shape.dimension==1
            rec=shape.geometry.geometry
            (haskey(rec,:face) && haskey(rec.face,:pcurves)) || continue
            kernel.shapes[id]=_ModelOCCShape(1,
                (kind=shape.geometry.kind,
                 geometry=merge(rec,(face=merge(rec.face,
                     (pcurves=_model_occ_rekey_dict(
                         rec.face.pcurves,tagmap,flipped,tol,caller),)),))),
                shape.children)
        end
    end
    return nothing
end

# A transform that is not a similarity (a dilation with unequal factors, a
# general affine matrix) takes a sphere to an ellipsoid and a circle to an
# ellipse: no sphere, cylinder, cone or torus record can hold the image.
# OpenCASCADE converts such shapes to NURBS first (`BRepBuilderAPI_GTransform`)
# and so does this: the faces and edges of the closure become B-splines, on
# which any affine map acts exactly through the poles, and a solid's compact
# primitive record gives way to its shell. Planes are converted too, into
# degree-1 patches: a plane record whose frame took the stretch (X = (2,0,0)
# after a dilation by 2 along x) is no orthonormal frame any more, and what
# meshes a plane assumes one (the cap of a dilated cone had half of Gmsh's
# triangles, with edges 4.6 sizes long). Nothing happens for a similarity.
function _model_occ_affine_prepare!(m::GeoModel,kernel,union,dim_tags,
                                    t::_AffineTransform,caller)
    _linear_isotropy(t.linear)===nothing || return nothing
    curved=false
    for id in union
        shape=kernel.shapes[id]
        g=shape.geometry
        if shape.dimension==2
            rec=g.geometry
            # an imported face on its surface's natural bounds keeps the
            # refusal: its domain bookkeeping does not survive the conversion
            rec isa NamedTuple && haskey(rec,:natural_bounds) && return nothing
            rec isa NamedTuple && !haskey(rec,:patch) &&
                get(rec,:occ,nothing) in (:sphere,:cylinder,:cone,:torus) &&
                (curved=true)
        elseif shape.dimension==1
            rec=g.geometry
            rec isa NamedTuple && get(rec,:occ,nothing) in (:circle,:ellipse) &&
                (curved=true)
        elseif shape.dimension==3
            (haskey(g,:primitive) || haskey(g,:sector)) && (curved=true)
        end
    end
    # A closure of planes is turned into patches when they carry their
    # parameter bounds and their edges are lines. With any other edge (a
    # hyperbola read from a BREP file keeps its own parameter, which the
    # conversion would replace) it keeps its plane records, as before.
    curved || all(id->begin
            shape=kernel.shapes[id]
            rec=shape.dimension in (1,2) ? shape.geometry.geometry : nothing
            if shape.dimension==2
                rec isa NamedTuple && (haskey(rec,:patch) ||
                    (get(rec,:occ,nothing)===:plane && haskey(rec,:uvb)))
            elseif shape.dimension==1
                rec isa NamedTuple && get(rec,:occ,nothing)===:line
            else
                true
            end
        end,union) || return nothing
    _model_occ_convert_to_nurbs!(m,dim_tags,caller)
    # The conversion keeps each quadric's analytic record beside its patch
    # (the shadow a similarity carries along); here the patch is all that
    # survives, so the faces become plain tensor records.
    replaced=Tuple{Any,NamedTuple}[]
    for id in union
        shape=kernel.shapes[id]
        g=shape.geometry
        if shape.dimension==3
            (haskey(g,:primitive) || haskey(g,:sector)) &&
                (kernel.shapes[id]=_ModelOCCShape(3,(kind=:solid,),shape.children))
        elseif shape.dimension==2
            rec=g.geometry
            (rec isa NamedTuple && haskey(rec,:patch) &&
             (haskey(rec,:analytic_chart) ||
              get(rec,:occ,nothing)===:plane)) || continue
            plain=(occ=:tensor_bspline,patch=rec.patch,pcurves=rec.pcurves,
                   reversed=get(rec,:reversed,false),
                   vertex_uv=get(rec,:vertex_uv,Dict{Int,NTuple{2,Float64}}()))
            # A plane keeps its frame for the measures (the patch alone
            # would be measured over its whole rectangle, or its trimmed
            # region by quadrature): the exact area of a stretched cap.
            get(rec,:occ,nothing)===:plane && (plain=merge(plain,(trimmed=true,
                plane_rec=Base.structdiff(rec,(patch=nothing,)))))
            push!(replaced,(rec,plain))
            kernel.shapes[id]=_ModelOCCShape(2,
                merge(g,(kind=:tensor_bspline,geometry=plain)),shape.children)
        end
    end
    # pole and iso edges embed their face's record: follow it
    for id in union
        shape=kernel.shapes[id]
        shape.dimension==1 || continue
        rec=shape.geometry.geometry
        if rec isa NamedTuple &&
                get(rec,:converted_from,nothing) in (:circle,:ellipse,:line)
            # a converted line or conic keeps only its B-spline as well
            kernel.shapes[id]=_ModelOCCShape(1,merge(shape.geometry,
                (geometry=(occ=:copied_curve,curve=rec.curve,
                           t0=rec.t0,t1=rec.t1),)),shape.children)
            continue
        end
        (rec isa NamedTuple && haskey(rec,:face)) || continue
        k=findfirst(pair->pair[1]===rec.face,replaced)
        k===nothing && continue
        kernel.shapes[id]=_ModelOCCShape(1,merge(shape.geometry,
            (geometry=merge(rec,(face=replaced[k][2],)),)),shape.children)
    end
    return nothing
end

# Move the bound vertex `id` from tag `from` to the free tag `to`, in the
# kernel and in what the model has already imported of it: the point record,
# the end points of the curves that use it and the vertex tables of their
# faces. Its tag is then free for another shape.
function _model_occ_rebind_owned_vertex!(m::GeoModel,kernel,id,from,to,caller)
    point=kernel.shapes[id].geometry.point
    _model_occ_retag_shape!(kernel,0,from,to,caller)
    if haskey(kernel.mirrored_points,from)
        delete!(kernel.mirrored_points,from)
        kernel.mirrored_points[to]=point
    end
    if haskey(m.points,from)
        m.points[to]=m.points[from]
        haskey(m.point_size,from) && (m.point_size[to]=m.point_size[from])
        m.next_tag[1]=max(m.next_tag[1],to)
        # the old tag's record goes with the unbind; `synchronize` imports
        # the shape that takes the tag
        delete!(m.points,from);delete!(m.point_size,from)
        delete!(m.point_w,from)
    end
    rename(tags)=map(tag->tag==from ? to : tag,tags)
    for ((dimension,tag),owner) in kernel.bindings
        shape=kernel.shapes[owner]
        if dimension==1 && id in shape.children
            haskey(m.curves,tag) && (m.curves[tag]=rename(m.curves[tag]))
        elseif dimension==2
            record=shape.geometry.geometry
            (record isa NamedTuple && haskey(record,:vertex_uv) &&
             haskey(record.vertex_uv,from)) || continue
            _model_occ_contains_shape(kernel,owner,id) || continue
            table=copy(record.vertex_uv)
            table[to]=pop!(table,from)
            record=merge(record,(vertex_uv=table,))
            kernel.shapes[owner]=_ModelOCCShape(2,
                merge(shape.geometry,(geometry=record,)),shape.children)
            haskey(m.surface_geometry,tag) && (m.surface_geometry[tag]=record)
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
    _model_occ_affine_prepare!(m,kernel,union,
        Tuple{Int,Int}[(d,tg) for (d,tg,_) in entries],t,caller)
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
    # Reflections range-negate closed conic edge records — collect their
    # pre-transform tags so `rekey` can fold `s = -t` into their pcurves.
    flipped=Set{Int}()
    if _linear_det(t.linear)<0
        for id in union
            shape=kernel.shapes[id]
            shape.dimension==1 || continue
            rec=shape.geometry.geometry
            hasproperty(rec,:occ) || continue
            rec.occ in (:circle,:ellipse) || continue
            rec.t1-rec.t0 < _OCC_TWO_PI && continue
            haskey(oldtags,id) && push!(flipped,oldtags[id])
        end
    end
    # Every fallible record transform runs before the first bind/unbind —
    # representability failures leave the kernel untouched.
    newshapes=Dict{UInt64,_ModelOCCShape}()
    for id in union
        shape=kernel.shapes[id]
        label=_model_periodic_entity_label(shape.dimension)
        what="OCC $label[$(get(oldtags,id,0))]"
        geometry=_model_occ_remap_uses(
            _model_occ_transform_geometry(m,shape,t,caller,what;
                oldtag=get(oldtags,id,0)),remap)
        if shape.dimension==-2 && _linear_det(t.linear)<0
            # A reflection maps a quadric's frame vector by vector, which
            # leaves it left-handed: the face's natural normal ∂u×∂v then
            # points against the image of the old one, so the shell has to
            # use the face the other way round (OpenCASCADE reverses the
            # faces of a mirrored shape for the same reason). Without it a
            # mirrored Boolean result measured its curved faces with the
            # wrong sign: box minus cylinder went from 1.7173 to 2.0942.
            geometry=merge(geometry,(uses=Tuple{UInt64,Bool}[
                (new,_model_occ_reflection_flips(kernel.shapes[old]) ? !rev : rev)
                for ((new,rev),(old,_)) in
                    zip(geometry.uses,shape.geometry.uses)],))
        end
        newshapes[remap[id]]=_ModelOCCShape(shape.dimension,geometry,
            UInt64[remap[c] for c in shape.children])
    end
    # Upstream's `_transform` unbinds and rebinds one input at a time:
    # `_unbind(old, tag, true)` releases the input and the sub-shapes nothing
    # else uses (the maxima drop to the surviving entities), then `_bind(new,
    # tag, true)` emits each re-created sub-shape not bound yet at
    # `getMaxTag+1` in original enumeration order. A later input sees the
    # tags the earlier ones took (Gmsh 4.15.2: translating two boxes renumbers
    # the second box's points 25-32, a copied square's points 9-12).
    for id in union
        x=get(kernel.explorer,id,nothing)
        x===nothing || (kernel.explorer[remap[id]]=UInt64[remap[c] for c in x])
    end
    newtags=Dict{UInt64,Int}()
    ctx=(inverse=inverse,oldtags=oldtags,newtags=newtags,
         tagmap=Dict{Int,Int}(),tagmap0=Dict{Int,Int}(),caller=caller)
    closures=Dict{UInt64,Vector{UInt64}}(
        id=>_model_occ_closure(kernel,UInt64[id]) for (_,_,id) in entries)
    walks=Dict{UInt64,Dict{UInt64,Int}}(
        id=>_model_occ_walk_rank(kernel,id) for (_,_,id) in entries)
    for (d,tg,id) in entries
        _model_occ_unbind!(kernel,id,true)
        # `getMaxTag` counts every live GModel entity — including lazily-bound
        # `.geo` records a shared-boundary transform leaves alone — so rebind
        # tags can never overwrite an unrelated public tag.
        for dd in 0:3
            store=dd==0 ? m.points : dd==1 ? m.curves : dd==2 ? m.surfaces :
                m.volumes
            for t in keys(store)
                (dd,t) in kernel.removed && continue
                kernel.maxima[dd+3]=max(kernel.maxima[dd+3],t)
            end
        end
        for t in keys(m.loops)
            (-1,t) in kernel.removed || (kernel.maxima[2]=max(kernel.maxima[2],t))
        end
        for t in keys(m.surface_loops)
            (-2,t) in kernel.removed || (kernel.maxima[1]=max(kernel.maxima[1],t))
        end
        order=Tuple{Int,Int,UInt64}[]
        rank=walks[id]
        for cid in closures[id]
            # the input itself keeps its tag (`_bind(new, tag)`); a point's
            # still-used old vertex needs one (below)
            cid==id && d!=0 && continue
            haskey(oldtags,cid) || continue
            nid=remap[cid]
            (haskey(newtags,nid) || haskey(kernel.shape_tags,nid)) && continue
            push!(order,(kernel.shapes[cid].dimension,rank[cid],nid))
        end
        sort!(order)
        for (od,old,nid) in order
            od>=-2 || continue
            t=_model_occ_next_tag(kernel,od,-1,caller)
            kernel.maxima[od+3]=t
            newtags[nid]=t
        end
        if d==0 && haskey(kernel.shape_tags,id)
            # The ownership-protected unbind left a vertex that a bound edge
            # still uses. Upstream's `_bind` then hands its tag to the moved
            # copy, which is a free vertex; the corner the tag came from
            # stays where it was, part of its solid, and is bound again
            # under the next free tag (Gmsh 4.15.2, `Translate {0,0,1} {
            # Point{1}; }` on a box: nine points, the corner (0,0,1) is
            # Point 9 and Curves 1, 2 and 10 end there). It used to be
            # dragged along with the copy.
            _model_occ_rebind_owned_vertex!(m,kernel,id,tg,
                newtags[remap[id]],caller)
        end
        _model_occ_bind_tree!(kernel,newshapes,remap[id],d,tg,ctx)
    end
    _model_occ_rekey_pcurves!(kernel,newshapes,ctx.tagmap;
        flipped=flipped,tagmap0=ctx.tagmap0,caller=caller)
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
    for id in union
        x=get(kernel.explorer,id,nothing)
        x===nothing || (kernel.explorer[remap[id]]=UInt64[remap[c] for c in x])
    end
    ctx=(inverse=inverse,oldtags=oldtags,newtags=Dict{UInt64,Int}(),
         tagmap=Dict{Int,Int}(),tagmap0=Dict{Int,Int}(),caller=caller)
    out=Tuple{Int,Int}[]
    for (d,id) in units
        nid=remap[id]
        # A copy already bound by an earlier unit's recursion (deduplicated
        # shared sub-shapes, or the same entity listed twice) is not reported
        # again — `_multiBind` only pushes unbound (`!exists`) shapes.
        haskey(kernel.shape_tags,nid) && continue
        tag=_model_occ_next_tag(kernel,d,-1,caller)
        # `BRepBuilderAPI_Copy` keeps the shape's structure, so `_bind`'s
        # walk meets the copied sub-shapes in the order it met the
        # originals: each dimension's copies take the next tags in the
        # order of their originals' tags.
        counters=copy(kernel.maxima)
        order=Tuple{Int,Int,UInt64}[]
        rank=_model_occ_walk_rank(kernel,id)
        for cid in _model_occ_closure(kernel,UInt64[id])
            cid==id && continue
            haskey(oldtags,cid) || continue
            haskey(kernel.shape_tags,remap[cid]) && continue
            haskey(ctx.newtags,remap[cid]) && continue
            push!(order,(kernel.shapes[cid].dimension,rank[cid],remap[cid]))
        end
        sort!(order)
        for (cd,_,cnid) in order
            cd>=0 || continue
            counters[cd+3]+=1
            ctx.newtags[cnid]=counters[cd+3]
        end
        _model_occ_bind_tree!(kernel,newshapes,nid,d,tag,ctx)
        push!(out,(d,tag))
    end
    _model_occ_rekey_pcurves!(kernel,newshapes,ctx.tagmap;
        tagmap0=ctx.tagmap0)
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
