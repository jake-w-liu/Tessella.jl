# `gmsh.model.occ` wire and planar-face producers (`addWire`, `addCurveLoop`,
# `addPlaneSurface`, `addRectangle`). Every producer is failure-atomic: all
# input validation and tag/shape reservations precede the publication loop,
# so a rejected call leaves the registry unchanged. `shape_of`/`tag_of`
# resolvers let the same plane-fit logic run against bound kernel shapes
# (`addPlaneSurface`) and pending pre-publication shapes (`addRectangle`).

# Vertices coincide within OCCT's sewing tolerance (`Precision::Confusion`) —
# `BRepBuilderAPI_MakeWire` chains edges on vertex coincidence.
@inline _model_occ_vertices_coincide(p,q)=
    _occ_modulus(_arc_sub(p,q))<=_OCC_CONFUSION

# `BRepBuilderAPI_MakeWire` head-to-tail chaining: the first edge always goes
# in FORWARD and each later edge orients itself against the running tail; a
# non-chaining edge fails the wire. Returns `(uses,head,tail)` where `uses`
# are `(shape id, reversed)` pairs and `head`/`tail` are the wire's open
# endpoint positions.
function _model_occ_wire_uses(kernel,curve_tags,caller)
    uses=Tuple{UInt64,Bool}[]; head=tail=nothing
    for (i,raw) in enumerate(curve_tags)
        t=abs(_model_occ_requested_tag(raw,caller))
        id=get(kernel.bindings,(1,t),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Curve[$t]"))
        a,b=kernel.shapes[id].children
        pa=kernel.shapes[a].geometry.point; pb=kernel.shapes[b].geometry.point
        rev=false
        if i==1
            head=pa; tail=pb
        elseif _model_occ_vertices_coincide(tail,pa)
            tail=pb
        elseif _model_occ_vertices_coincide(tail,pb)
            rev=true; tail=pa
        else
            throw(ArgumentError("$caller: Could not create wire"))
        end
        push!(uses,(id,rev))
    end
    return uses,head,tail
end

# `OCC_Internals::addWire`: `addCurveLoop` is the same call with
# `checkClosed=true` (`TopoDS_Wire::Closed` on the sewn vertices). A negative
# first curve tag binds `w.Reversed()` — reversed order AND orientation.
function _model_occ_add_wire!(m::GeoModel,curve_tags;tag=-1,check_closed=false)
    caller="OCC.add_wire"
    kernel=m.occ_kernel
    check_closed isa Bool || throw(ArgumentError(
        "$caller: checkClosed must be a boolean"))
    allocated=_model_occ_next_tag(kernel,-1,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(-1,allocated)) && throw(ArgumentError(
        "$caller: OCC wire or curve loop with tag $allocated already exists"))
    isempty(curve_tags) && throw(ArgumentError("$caller: Could not create wire"))
    uses,head,tail=_model_occ_wire_uses(kernel,curve_tags,caller)
    check_closed && !_model_occ_vertices_coincide(head,tail) && throw(
        ArgumentError("$caller: Curve loop is not closed"))
    first(curve_tags)<0 && (uses=[(id,!rev) for (id,rev) in reverse(uses)])
    shape=_ModelOCCShape(-1,(uses=uses,),UInt64[id for (id,_) in uses])
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,-1,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end

_model_occ_add_curve_loop!(m,curve_tags;tag=-1)=
    _model_occ_add_wire!(m,curve_tags;tag=tag,check_closed=true)

# `OCC_Internals::addCurveLoops`: closed wires are searched among the given
# curves — endpoints connect on `Precision::Confusion` coincidence
# (`ShapeAnalysis_FreeBounds` sewing semantics), not only shared vertex
# entities. Each closed walk publishes one curve loop under a fresh tag;
# curves that no closed walk consumes come back unused. Returns
# `(loop_tags, leftover_tags)`.
function _model_occ_curve_loops(m::GeoModel,curve_tags)
    caller="OCC.add_curve_loops"
    kernel=m.occ_kernel
    edges=UInt64[]; labels=Int[]; seen=Set{Int}()
    for raw in curve_tags
        t=abs(_model_occ_requested_tag(raw,caller))
        id=get(kernel.bindings,(1,t),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Curve[$t]"))
        t in seen && continue
        push!(seen,t); push!(edges,id); push!(labels,t)
    end
    isempty(edges) && return Int[],Int[]
    n=length(edges)
    # Cluster coincident endpoints: port 2(e-1)+s of edge e → cluster index.
    coords=NTuple{3,Float64}[]; port=fill(0,2n)
    for e in 1:n, s in 1:2
        vid=kernel.shapes[edges[e]].children[s]
        p=kernel.shapes[vid].geometry.point
        c=findfirst(q->_model_occ_vertices_coincide(p,q),coords)
        c===nothing && (push!(coords,p); c=length(coords))
        port[2(e-1)+s]=c
    end
    incident=[Int[] for _ in coords]
    for e in 1:n, s in 1:2
        push!(incident[port[2(e-1)+s]],e)
    end
    used=falses(n); loops=Vector{Int}[]
    for e0 in 1:n
        used[e0] && continue
        pa,pb=port[2*e0-1],port[2*e0]
        if pa==pb
            used[e0]=true; push!(loops,[labels[e0]]); continue
        end
        # Greedy walk: `pos_of` maps each entered cluster to the chain index
        # just before it, so re-entering any cluster closes the sub-chain
        # (handles lollipop stems and figure-8 shared vertices, not only
        # plain cycles that close at the start).
        chain=[e0]; cur=pb; used[e0]=true
        pos_of=Dict(pa=>0); closed=false
        while true
            if haskey(pos_of,cur)
                k=pos_of[cur]
                used[chain[1:k]].=false
                push!(loops,[labels[e] for e in chain[k+1:end]])
                closed=true; break
            end
            pos_of[cur]=length(chain)
            nxt=0
            for f in incident[cur]
                used[f] || (nxt=f; break)
            end
            nxt==0 && break
            used[nxt]=true; push!(chain,nxt)
            cur= port[2nxt-1]==cur ? port[2nxt] : port[2nxt-1]
        end
        closed || (used[chain].=false)   # open chain — edges stay reusable
    end
    loop_tags=[_model_occ_add_wire!(m,chain;tag=-1,check_closed=true)
               for chain in loops]
    consumed=Set{Int}(Iterators.flatten(loops))
    return loop_tags,[t for t in labels if !(t in consumed)]
end

# Traversal-ordered boundary samples of a wire: vertex positions chained by
# `uses` orientation, with quarter-arc interior points on circle/ellipse
# edges so a curved boundary contributes to the plane fit. Only `:line`,
# `:circle`, and `:ellipse` OCC edges have the analytic pcurves planar
# faces require.
function _model_occ_wire_samples(shape_of,uses,caller)
    pts=NTuple{3,Float64}[]
    for (id,rev) in uses
        edge=shape_of(id); g=edge.geometry.geometry
        g.occ in (:line,:circle,:ellipse) || throw(ArgumentError(
            "$caller: only line, circle, and ellipse OCC curves are " *
            "supported on plane surfaces"))
        a,b=edge.children
        pa=shape_of(a).geometry.point; pb=shape_of(b).geometry.point
        start=rev ? pb : pa; stop=rev ? pa : pb
        # Chained edges repeat the junction vertex, and closed wires repeat
        # the head: dedupe so symmetric boundaries keep a centered mean.
        (isempty(pts) || !_model_occ_vertices_coincide(last(pts),start)) &&
            push!(pts,start)
        if g.occ===:circle
            for k in 1:3
                fraction=rev ? (4-k)/4 : k/4
                push!(pts,_occ_circle_point(g,g.t0+(g.t1-g.t0)*fraction))
            end
        elseif g.occ===:ellipse
            for k in 1:3
                fraction=rev ? (4-k)/4 : k/4
                push!(pts,_occ_ellipse_point(g,g.t0+(g.t1-g.t0)*fraction))
            end
        end
        _model_occ_vertices_coincide(last(pts),stop) || push!(pts,stop)
    end
    length(pts)>1 && _model_occ_vertices_coincide(pts[1],last(pts)) && pop!(pts)
    return pts
end

# Newell's normal is exact for planar polygons and robust on sampled curved
# boundaries; its sign follows the traversal orientation.
function _model_occ_newell_normal(pts)
    nx=ny=nz=0.0; n=length(pts)
    for i in 1:n
        p=pts[i]; q=pts[i==n ? 1 : i+1]
        nx+=(p[2]-q[2])*(p[3]+q[3])
        ny+=(p[3]-q[3])*(p[1]+q[1])
        nz+=(p[1]-q[1])*(p[2]+q[2])
    end
    return (nx,ny,nz)
end

# The 2-D image of a line/circle/ellipse OCC edge on the fitted plane: the
# `lin2d`/`circ2d`/`ellips2d` records `_occ_pcurve_d0` evaluates. The 3-D
# curve's stored `X`/`Y` may be transformed out of the Y=n×X orientation,
# so the 2-D frame is projected (never re-crossed) and must stay
# counterclockwise — the face axis is chosen so every boundary conic
# agrees.
function _model_occ_plane_pcurve(g,center,X,Y)
    if g.occ===:line
        o=_arc_sub(g.origin,center)
        return (lin2d=((_dot3(o,X),_dot3(o,Y)),
                       (_dot3(g.dir,X),_dot3(g.dir,Y))),)
    end
    o=_arc_sub(g.center,center)
    g.occ===:ellipse && return (ellips2d=(
        (_dot3(o,X),_dot3(o,Y)),(_dot3(g.X,X),_dot3(g.X,Y)),g.rx,g.ry),)
    return (circ2d=((_dot3(o,X),_dot3(o,Y)),
                    (_dot3(g.X,X),_dot3(g.X,Y)),g.r),)
end

# `BRepBuilderAPI_MakeFace` on the outer wire, `BRepBuilderAPI_FindPlane`
# semantics: the fitted plane must contain every wire's boundary within
# `Precision::Confusion`. The face axis is the traversal normal, flipped when
# every boundary circle is oppositely wound so their `circ2d` projections
# stay counterclockwise.
function _model_occ_plane_face_geometry(shape_of,tag_of,uses_per_wire,caller)
    # BRepBuilderAPI_MakeFace requires each wire closed on its sewn ends.
    for uses in uses_per_wire
        first_edge=shape_of(first(uses)[1]); last_edge=shape_of(last(uses)[1])
        a0,b0=first_edge.children; aN,bN=last_edge.children
        head=first(uses)[2] ? shape_of(b0).geometry.point :
                              shape_of(a0).geometry.point
        tail=last(uses)[2] ? shape_of(aN).geometry.point :
                             shape_of(bN).geometry.point
        _model_occ_vertices_coincide(head,tail) || throw(ArgumentError(
            "$caller: Could not create surface"))
    end
    outer=_model_occ_wire_samples(shape_of,first(uses_per_wire),caller)
    axis_raw=_model_occ_newell_normal(outer)
    _occ_modulus(axis_raw)>_OCC_CONFUSION || throw(ArgumentError(
        "$caller: Could not create surface"))
    axis=_occ_dir(axis_raw)
    center=(sum(p[1] for p in outer)/length(outer),
            sum(p[2] for p in outer)/length(outer),
            sum(p[3] for p in outer)/length(outer))
    samples=copy(outer)
    for uses in uses_per_wire[2:end]
        append!(samples,_model_occ_wire_samples(shape_of,uses,caller))
    end
    scale=0.0
    for p in samples
        scale=max(scale,_occ_modulus(_arc_sub(p,center)))
    end
    tolerance=max(_OCC_CONFUSION,1e-11*scale)
    for p in samples
        abs(_dot3(_arc_sub(p,center),axis))<=tolerance || throw(ArgumentError(
            "$caller: Could not create surface"))
    end
    # Conic normals must agree with the face axis; a unanimous disagreement
    # flips the fit, a mixed set cannot be projected into a `circ2d`/
    # `ellips2d`.
    signs=Int[]
    for uses in uses_per_wire,(id,_) in uses
        g=shape_of(id).geometry.geometry
        g.occ in (:circle,:ellipse) || continue
        s=_dot3(g.n,axis)
        abs(s)>_OCC_CONFUSION || throw(ArgumentError(
            "$caller: boundary conic is not in the face plane"))
        push!(signs,s>0.0 ? 1 : -1)
    end
    if !isempty(signs)
        if all(==(-1),signs)
            axis=(.-axis)
        elseif !all(==(1),signs)
            throw(ArgumentError(
                "$caller: boundary conics have inconsistent orientations"))
        end
    end
    X,Y=_occ_ax2_frame(axis)
    pcurves=Dict{Int,NamedTuple}()
    umin=vmin=Inf; umax=vmax=-Inf
    for uses in uses_per_wire,(id,_) in uses
        g=shape_of(id).geometry.geometry
        pc=_model_occ_plane_pcurve(g,center,X,Y)
        pcurves[tag_of(id)]=(fwd=pc,rev=nothing)
        lo,hi=_occ_pcurve_fillbndbox(pc,g.t0,g.t1)
        umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
        umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
    end
    return (occ=:plane,center=center,axis=axis,X=X,Y=Y,reversed=false,
            pcurves=pcurves,uvb=([umin,vmin],[umax,vmax]))
end

function _model_occ_add_plane_surface!(m::GeoModel,wire_tags;tag=-1)
    caller="OCC.add_plane_surface"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,2,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(2,allocated)) && throw(ArgumentError(
        "$caller: OCC Surface[$allocated] already exists"))
    isempty(wire_tags) && throw(ArgumentError(
        "$caller: Plane surface requires at least one curve loop"))
    wire_ids=UInt64[]
    for raw in wire_tags
        # All wire tags are > 0 for OCC; upstream takes `abs` for
        # compatibility with the GEO factory's signed loops.
        t=abs(_model_occ_requested_tag(raw,caller))
        id=get(kernel.bindings,(-1,t),nothing)
        id===nothing && throw(ArgumentError("$caller: unknown OCC Wire[$t]"))
        push!(wire_ids,id)
    end
    uses_per_wire=[kernel.shapes[id].geometry.uses for id in wire_ids]
    shape_of=id->kernel.shapes[id]
    tag_of=id->_model_occ_shape_tag(kernel,id,caller)
    geometry=_model_occ_plane_face_geometry(shape_of,tag_of,uses_per_wire,caller)
    shape=_ModelOCCShape(2,(kind=:plane,geometry=geometry),wire_ids)
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,2,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# `OCC_Internals::makeRectangle`: an axis-parallel four-line face, or eight
# edges with quarter-circle corners when `rounded_radius>0`. Sub-entities
# (vertices, edges, the wire) publish alongside the face — the recursive
# upstream `_bind` exposes them too. Corner arcs use the center-form
# geometry `_model_occ_circle_arc_geometry` produces for `addCircleArc`.
function _model_occ_add_rectangle!(m::GeoModel,x,y,z,dx,dy;
                                   tag=-1,rounded_radius=0.0)
    caller="OCC.add_rectangle"
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,2,
        _model_occ_requested_tag(tag,caller),caller)
    haskey(kernel.bindings,(2,allocated)) && throw(ArgumentError(
        "$caller: OCC Surface[$allocated] already exists"))
    ox=_finite_scalar(x,caller,"x"); oy=_finite_scalar(y,caller,"y")
    oz=_finite_scalar(z,caller,"z")
    w=_finite_scalar(dx,caller,"dx"); h=_finite_scalar(dy,caller,"dy")
    (w!=0.0 && h!=0.0) || throw(ArgumentError(
        "$caller: Rectangle with zero width or height"))
    rr=_finite_scalar(rounded_radius,caller,"roundedRadius")
    x1=w>0.0 ? ox : ox+w; x2=w>0.0 ? ox+w : ox
    y1=h>0.0 ? oy : oy+h; y2=h>0.0 ? oy+h : oy
    rounded=rr>0.0
    n_edges=rounded ? 8 : 4
    kernel.maxima[3]<=typemax(Int32)-n_edges || throw(ArgumentError(
        "$caller: OCC vertex tag allocation exceeds Int32"))
    kernel.maxima[4]<=typemax(Int32)-n_edges || throw(ArgumentError(
        "$caller: OCC curve tag allocation exceeds Int32"))
    kernel.maxima[2]<typemax(Int32) || throw(ArgumentError(
        "$caller: OCC wire tag allocation exceeds Int32"))
    vertex_tags=[kernel.maxima[3]+i for i in 1:n_edges]
    edge_tags=[kernel.maxima[4]+i for i in 1:n_edges]
    wire_tag=kernel.maxima[2]+1
    points=NTuple{3,Float64}[]
    if rounded
        r=rr
        append!(points,((x1+r,y1,oz),(x2-r,y1,oz),(x2,y1+r,oz),(x2,y2-r,oz),
                        (x2-r,y2,oz),(x1+r,y2,oz),(x1,y2-r,oz),(x1,y1+r,oz)))
    else
        append!(points,((x1,y1,oz),(x2,y1,oz),(x2,y2,oz),(x1,y2,oz)))
    end
    edge_shapes=_ModelOCCShape[]
    edge_children=Tuple{Int,Int}[]
    if rounded
        r=rr
        centers=NTuple{3,Float64}[(x2-r,y1+r,oz),(x2-r,y2-r,oz),
                                  (x1+r,y2-r,oz),(x1+r,y1+r,oz)]
        for i in 1:4
            a=2i-1; b=2i; nxt=b==8 ? 1 : b+1
            pa=points[a]; pb=points[b]
            difference=_arc_sub(pb,pa); length=_occ_modulus(difference)
            length>_OCC_CONFUSION || throw(ArgumentError(
                "$caller: cannot construct a coincident edge"))
            push!(edge_shapes,_ModelOCCShape(1,
                (kind=:line,geometry=(occ=:line,origin=pa,
                 dir=_occ_dir(difference),t0=0.0,t1=length)),
                UInt64[]))
            push!(edge_children,(a,b))
            geometry,error=_model_occ_circle_arc_geometry(pb,centers[i],
                                                          points[nxt])
            error===nothing || throw(ArgumentError("$caller: $error"))
            push!(edge_shapes,_ModelOCCShape(1,
                (kind=:circle,geometry=geometry),UInt64[]))
            push!(edge_children,(b,nxt))
        end
    else
        for i in 1:4
            a=i; b=i==4 ? 1 : i+1
            pa=points[a]; pb=points[b]
            difference=_arc_sub(pb,pa); length=_occ_modulus(difference)
            length>_OCC_CONFUSION || throw(ArgumentError(
                "$caller: cannot construct a coincident edge"))
            push!(edge_shapes,_ModelOCCShape(1,
                (kind=:line,geometry=(occ=:line,origin=pa,
                 dir=_occ_dir(difference),t0=0.0,t1=length)),
                UInt64[]))
            push!(edge_children,(a,b))
        end
    end
    ids=_model_occ_reserve_shapes(kernel,2n_edges+2,caller)
    pending_shapes=Dict{UInt64,_ModelOCCShape}()
    pending_tags=Dict{UInt64,Int}()
    pending=Tuple{Int,Int,UInt64,_ModelOCCShape}[]
    for (i,point) in enumerate(points)
        id=ids[i]
        shape=_ModelOCCShape(0,(point=point,mesh_size=0.0),UInt64[])
        pending_shapes[id]=shape
        push!(pending,(0,vertex_tags[i],id,shape))
    end
    for i in 1:n_edges
        a,b=edge_children[i]
        id=ids[n_edges+i]
        shape=_ModelOCCShape(1,edge_shapes[i].geometry,
                             UInt64[ids[a],ids[b]])
        pending_shapes[id]=shape; pending_tags[id]=edge_tags[i]
        push!(pending,(1,edge_tags[i],id,shape))
    end
    uses=Tuple{UInt64,Bool}[(ids[n_edges+i],false) for i in 1:n_edges]
    wire_id=ids[2n_edges+1]; face_id=ids[2n_edges+2]
    shape_of=id->pending_shapes[id]
    tag_of=id->pending_tags[id]
    geometry=_model_occ_plane_face_geometry(shape_of,tag_of,(uses,),caller)
    wire=_ModelOCCShape(-1,(uses=uses,),UInt64[first(use) for use in uses])
    face=_ModelOCCShape(2,(kind=:plane,geometry=geometry),UInt64[wire_id])
    pending_shapes[wire_id]=wire; pending_shapes[face_id]=face
    push!(pending,(-1,wire_tag,wire_id,wire),(2,allocated,face_id,face))
    for (dimension,current,id,shape) in pending
        _model_occ_bind_shape!(kernel,dimension,current,id,shape)
    end
    kernel.next_shape=face_id+UInt64(1)
    return allocated
end
