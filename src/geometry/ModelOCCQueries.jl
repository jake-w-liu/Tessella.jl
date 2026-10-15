# `gmsh.model.occ` kernel-side queries — `OCC_Internals::getBoundingBox`,
# `getEntitiesInBoundingBox`, `getCurveLoops`, `getSurfaceLoops` and the
# `occ.mesh.setSize` point attribute. Upstream queries the bound TopoDS
# registry directly (no synchronization), so these operate on
# `kernel.shapes`/`kernel.bindings` — the same graph `synchronize` imports.
# Boxes are exact where the synced `_model_*_bounding_box` paths are exact:
# line endpoints, circle/ellipse/torus extrema, the sphere's ball box,
# NURBS control hulls padded by 1e-7 like `Bnd_Box`, and boundary unions
# elsewhere — which stays exact for plane faces and ruled rim faces whose
# axis-aligned extrema lie on their boundary curves.

# 6-tuple box for one point or an arbitrary point iterator.
@inline _model_occ_box_point(p::NTuple{3,Float64})=
    (p[1],p[2],p[3],p[1],p[2],p[3])

function _model_occ_box_points(points)
    bounds=_MODEL_EMPTY_BOUNDS
    for p in points
        bounds=_model_bounds_include(bounds,(p[1],p[2],p[3]))
    end
    return bounds
end

# NURBS control-polytope box — the `Bnd_Box`-style padded hull the synced
# `_model_entity_bounding_box` reports for tensor surfaces.
function _model_occ_controls_box(controls)
    bounds=_MODEL_EMPTY_BOUNDS
    for p in controls
        bounds=_model_bounds_include(bounds,(p[1],p[2],p[3]))
    end
    return (bounds[1]-1e-7,bounds[2]-1e-7,bounds[3]-1e-7,
            bounds[4]+1e-7,bounds[5]+1e-7,bounds[6]+1e-7)
end

# A partial sphere patch (records carrying `lat1`/`lat2`/`angle`) can
# reach an axis-aligned extremum in the patch interior — the point whose
# radial direction is exactly ±ê — which no boundary edge contains. The
# box is the boundary-edge union plus those interior candidates.
function _occ_sphere_face_bbox(kernel,shape,rec)
    bounds=_MODEL_EMPTY_BOUNDS
    for wire in shape.children
        for edge in kernel.shapes[wire].children
            bounds=_model_bounds_union(bounds,
                _model_occ_shape_bbox(kernel,edge))
        end
    end
    a=get(rec,:angle,_OCC_TWO_PI)
    v1=get(rec,:lat1,-π/2); v2=get(rec,:lat2,π/2)
    for axis in 1:3, s in (-1.0,1.0)
        vc=_myasin(clamp(s*rec.axis[axis],-1.0,1.0))
        (v1-1e-12<=vc<=v2+1e-12) || continue
        uc=mod(_gm_atan2(s*rec.Y[axis],s*rec.X[axis]),_OCC_TWO_PI)
        uc>a+1e-12 && uc<_OCC_TWO_PI-1e-12 && continue
        d=ntuple(k->k==axis ? s : 0.0,3)
        bounds=_model_bounds_union(bounds,_model_occ_box_point(
            ntuple(k->rec.center[k]+rec.radius*d[k],3)))
    end
    return bounds
end

function _model_occ_shape_bbox(kernel,id)
    shape=kernel.shapes[id]
    d=shape.dimension; g=shape.geometry
    d==0 && return _model_occ_box_point(g.point)
    if d==1
        rec=g.geometry
        if rec.occ===:line
            return _model_occ_box_points(
                (_occ_line_point(rec,rec.t0),_occ_line_point(rec,rec.t1)))
        elseif rec.occ===:circle
            return _occ_circle_bounding_box(rec)
        elseif rec.occ===:ellipse
            return _occ_ellipse_bounding_box(rec)
        elseif rec.occ===:degenerate
            tag=kernel.shape_tags[id][2]
            return _model_occ_box_point(_occ_degenerate_point(rec,tag,rec.t0))
        elseif rec.occ===:surface_iso
            face=rec.face
            haskey(face,:patch) || throw(ErrorException(
                "OCC.get_bounding_box: iso curve on a non-NURBS face " *
                "is not supported"))
            surface=face.patch.surface
            controls=rec.axis==1 ?
                @view(surface.controls[:,rec.fixed==first(surface.knots_v) ?
                    1 : size(surface.controls,2)]) :
                @view(surface.controls[rec.fixed==first(surface.knots_u) ?
                    1 : size(surface.controls,1),:])
            return _model_occ_controls_box(controls)
        end
        # :copied_curve — the trimmed sub-range stays inside the full
        # control hull.
        return _model_occ_controls_box(rec.curve.curve.controls)
    elseif d==2
        rec=g.geometry
        rec.occ===:sphere && (haskey(rec,:lat1) || haskey(rec,:lat2) ||
            haskey(rec,:angle)) && return _occ_sphere_face_bbox(kernel,shape,rec)
        rec.occ===:sphere && return (
            rec.center[1]-rec.radius,rec.center[2]-rec.radius,
            rec.center[3]-rec.radius,
            rec.center[1]+rec.radius,rec.center[2]+rec.radius,
            rec.center[3]+rec.radius)
        rec.occ===:torus && return _occ_torus_bounding_box(rec)
        haskey(rec,:patch) &&
            return _model_occ_controls_box(rec.patch.surface.controls)
        # Plane faces and ruled rim faces (cylinder, cone) reach every
        # axis-aligned extremum on a boundary edge — the recursive child
        # union below covers them exactly.
    end
    # Generic recursion: wires → edges, shells → faces, solids → shells,
    # faces → wires. The union is exact whenever the boundary covers the
    # shape's extrema (all current producers).
    bounds=_MODEL_EMPTY_BOUNDS
    for child in shape.children
        bounds=_model_bounds_union(bounds,_model_occ_shape_bbox(kernel,child))
    end
    bounds[1]==Inf && throw(ErrorException(
        "OCC.get_bounding_box: OCC shape has no geometry"))
    return bounds
end

function _model_occ_get_bounding_box(m::GeoModel,dim,tag)
    caller="OCC.get_bounding_box"
    d=_model_occ_dimension(dim,caller)
    t=_model_occ_requested_tag(tag,caller)
    id=get(m.occ_kernel.bindings,(d,t),nothing)
    id===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE entity of dimension $d with tag $t"))
    return _model_occ_shape_bbox(m.occ_kernel,id)
end

# `getEntitiesInBoundingBox` — upstream iterates the bound public maps for
# dims 0..3 (`dim=-1` selects all) and keeps entities whose box is
# CONTAINED in the query box. Result order is dim-ascending like the
# upstream loop, tag-sorted inside each dim for determinism.
function _model_occ_entities_in_bounding_box(
        m::GeoModel,xmin,ymin,zmin,xmax,ymax,zmax,dim)
    caller="OCC.get_entities_in_bounding_box"
    dim isa Integer && !(dim isa Bool) && -1<=dim<=3 || throw(ArgumentError(
        "$caller: OCC dimension must be in -1:3"))
    box=NTuple{6,Float64}(_finite_scalar(v,caller,"bounding box coordinate")
                          for v in (xmin,ymin,zmin,xmax,ymax,zmax))
    kernel=m.occ_kernel
    out=Tuple{Int,Int}[]
    for d in 0:3
        (dim==-1 || dim==d) || continue
        for t in sort!(Int[last(k) for k in keys(kernel.bindings)
                           if first(k)==d])
            b=_model_occ_shape_bbox(kernel,kernel.bindings[(d,t)])
            (b[1]>=box[1] && b[4]<=box[4] && b[2]>=box[2] &&
             b[5]<=box[5] && b[3]>=box[3] && b[6]<=box[6]) &&
                push!(out,(d,t))
        end
    end
    return out
end

# `getCurveLoops` — `TopExp_Explorer(face, TopAbs_WIRE)` then
# `BRepTools_WireExplorer` for traversal-ordered edges. Our face children
# are bound wires and `uses` already stores the explorer order; every wire
# and edge in the kernel is bound, so the on-demand `_bind` upstream
# performs for unbound members never triggers.
function _model_occ_get_curve_loops(m::GeoModel,surface_tag)
    caller="OCC.get_curve_loops"
    tag=_model_occ_requested_tag(surface_tag,caller)
    kernel=m.occ_kernel
    id=get(kernel.bindings,(2,tag),nothing)
    id===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE surface with tag $tag"))
    loop_tags=Int[]; curve_tags=Vector{Int}[]
    for wire_id in kernel.shapes[id].children
        haskey(kernel.shape_tags,wire_id) || throw(ErrorException(
            "$caller: unbound wire inside Surface[$tag]"))
        push!(loop_tags,kernel.shape_tags[wire_id][2])
        push!(curve_tags,[
            begin
                haskey(kernel.shape_tags,edge_id) || throw(ErrorException(
                    "$caller: unbound curve inside Surface[$tag]"))
                kernel.shape_tags[edge_id][2]
            end for (edge_id,_) in kernel.shapes[wire_id].geometry.uses])
    end
    return loop_tags,curve_tags
end

# `getSurfaceLoops` — `TopExp_Explorer(solid, TopAbs_SHELL)` then per-shell
# `TopAbs_FACE` in `uses` order.
function _model_occ_get_surface_loops(m::GeoModel,volume_tag)
    caller="OCC.get_surface_loops"
    tag=_model_occ_requested_tag(volume_tag,caller)
    kernel=m.occ_kernel
    id=get(kernel.bindings,(3,tag),nothing)
    id===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE volume with tag $tag"))
    loop_tags=Int[]; surface_tags=Vector{Int}[]
    for shell_id in kernel.shapes[id].children
        haskey(kernel.shape_tags,shell_id) || throw(ErrorException(
            "$caller: unbound shell inside Volume[$tag]"))
        push!(loop_tags,kernel.shape_tags[shell_id][2])
        push!(surface_tags,[
            begin
                haskey(kernel.shape_tags,face_id) || throw(ErrorException(
                    "$caller: unbound surface inside Volume[$tag]"))
                kernel.shape_tags[face_id][2]
            end for (face_id,_) in kernel.shapes[shell_id].geometry.uses])
    end
    return loop_tags,surface_tags
end

# `occ.mesh.setSize` — `OCC_Internals::setMeshSize`: only dimension-0
# entities are handled and unbound tags are ignored silently. The vertex
# attribute rides the point record's `mesh_size`, which `synchronize`
# already transports to `m.point_size`.
function _model_occ_set_size!(m::GeoModel,dim_tags,size)
    caller="OCC.mesh.set_size"
    s=_finite_scalar(size,caller,"size")
    kernel=m.occ_kernel
    for entry in dim_tags
        d=_model_occ_dimension(first(entry),caller)
        d==0 || continue
        t=_model_occ_requested_tag(last(entry),caller)
        id=get(kernel.bindings,(0,t),nothing)
        id===nothing && continue
        g=kernel.shapes[id].geometry
        kernel.shapes[id]=_ModelOCCShape(0,
            merge(g,(mesh_size=s,)),kernel.shapes[id].children)
        kernel.changed=true
    end
    return nothing
end
