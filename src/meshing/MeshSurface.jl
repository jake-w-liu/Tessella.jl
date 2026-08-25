"""
    MeshSurface

Stage-2 surface meshing (PLAN.md §3, "mesh the flat/patch/coax faces"). Two
routes, both reducing a surface to the exact 2-D machinery of [`Mesh2D`](@ref):

* **Planar faces** — a face bounded by coplanar loops (outer + holes). Boundary
  loops are first discretized under the size field ([`Mesh1D`](@ref)), projected
  to an in-plane orthonormal frame, meshed as a PSLG (constrained Delaunay +
  size-driven Ruppert refinement), then lifted back to 3-D. Isometric ⇒ the 2-D
  quality/size guarantees hold exactly on the surface.

* **Parametric faces** `s(u,v)` — meshed in the parameter rectangle. The dedicated
  cylinder route uses a watertight structured grid and certifies every physical
  chord edge. For a general regular map, the parameter-space target is divided by
  the largest local coordinate stretch (an isotropic approximation), with
  independent physical-edge and physical-area postconditions.

All meshers return a 3-D triangle [`Mesh`](@ref).
"""
module MeshSurface

using ..MeshTypes: Mesh, ntris, nnodes, node, triangle_area, validate
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..Mesh1D: mesh_segment
using ..SizeField: AbstractSizeField, AbstractAnisoField, size_at, ConstantSize,
                   directional_size, metric_edge_length, _metric_curve_increment

export mesh_planar_face, mesh_cylinder_face, mesh_parametric_face
export PlaneFrame, plane_frame, project, lift

@inline function _surface_entity_context(value, caller::AbstractString)
    value === nothing && return nothing
    value isa Tuple && length(value) == 2 &&
        value[1] isa Integer && !(value[1] isa Bool) &&
        value[2] isa Integer && !(value[2] isa Bool) ||
        throw(ArgumentError("$caller: an entity context must be nothing or a " *
                            "(dimension, tag) integer tuple (got $value)"))
    dim=try
        Int(value[1])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: entity dimension exceeds the platform Int range"))
    end
    tag=try
        Int(value[2])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: entity tag exceeds the platform Int range"))
    end
    dim in 0:3 || throw(ArgumentError("$caller: entity dimension must be in 0:3"))
    tag>0 || throw(ArgumentError("$caller: entity tag must be positive"))
    return (dim,tag)
end

@inline _entity_tuple(value)=value isa Tuple && length(value)==2 &&
    value[1] isa Integer && !(value[1] isa Bool) &&
    value[2] isa Integer && !(value[2] isa Bool)

function _surface_boundary_entity(spec, nloops::Int, loop_index::Int,
                                  edge_index::Int, edge_count::Int, p, q,
                                  default, caller::AbstractString)
    spec===nothing && return default
    raw=if _entity_tuple(spec)
        spec
    elseif applicable(spec,loop_index,edge_index,p,q)
        spec(loop_index,edge_index,p,q)
    elseif spec isa AbstractVector
        Base.require_one_based_indexing(spec)
        length(spec)==nloops || throw(ArgumentError(
            "$caller: boundary_entities must have one entry per boundary loop"))
        entry=spec[loop_index]
        if entry===nothing || _entity_tuple(entry)
            entry
        elseif entry isa AbstractVector
            Base.require_one_based_indexing(entry)
            length(entry)==edge_count || throw(ArgumentError(
                "$caller: boundary_entities[$loop_index] must have one entry per loop edge"))
            entry[edge_index]
        else
            throw(ArgumentError(
                "$caller: each boundary_entities entry must be an entity tuple, nothing, or a per-edge vector"))
        end
    else
        throw(ArgumentError(
            "$caller: boundary_entities must be an entity tuple, a per-loop vector, or a callable (loop_index, edge_index, p, q)"))
    end
    return _surface_entity_context(raw,caller)
end

# ── vector helpers (tuples) ─────────────────────────────────────────────────────
@inline _sub(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _norm(a) = hypot(a[1],a[2],a[3])
@inline function _unit(a)
    n = _norm(a)
    (isfinite(n) && n > 0) ||
        throw(ArgumentError("MeshSurface: vector must have finite positive length (got $a)"))
    (a[1]/n,a[2]/n,a[3]/n)
end

@inline function _point3(p, caller::AbstractString)
    count=try
        length(p)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: each point must be an indexable three-coordinate value"))
    end
    count==3 || throw(ArgumentError(
        "$caller: each point must have exactly three coordinates (got $count)"))
    q=ntuple(3) do dimension
        raw=try
            p[dimension]
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: could not read point coordinate $dimension: $(sprint(showerror,err))"))
        end
        raw isa Bool && throw(ArgumentError(
            "$caller: point coordinate $dimension must not be Bool"))
        try
            Float64(raw)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: point coordinate $dimension must be Float64-representable: " *
                sprint(showerror,err)))
        end
    end
    (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) ||
        throw(ArgumentError("$caller: point has non-finite coordinates $q"))
    return q
end

# ── planar frame ────────────────────────────────────────────────────────────────
"""
    PlaneFrame(origin, u, v, normal)

Validated immutable orthonormal frame for a plane. All four arguments are
three-coordinate finite points/vectors; `u × v` must agree with `normal`.
Use [`plane_frame`](@ref) to derive a frame from a boundary loop.
"""
struct PlaneFrame
    o::NTuple{3,Float64}   # origin (a loop vertex)
    u::NTuple{3,Float64}   # in-plane unit axis 1
    v::NTuple{3,Float64}   # in-plane unit axis 2
    n::NTuple{3,Float64}   # unit normal

    function PlaneFrame(o::NTuple{3,Float64},u::NTuple{3,Float64},
                        v::NTuple{3,Float64},n::NTuple{3,Float64})
        all(isfinite,(o...,u...,v...,n...)) || throw(ArgumentError(
            "PlaneFrame: origin and axes must be finite"))
        tolerance=4096eps(Float64)
        for (name,axis) in (("u",u),("v",v),("normal",n))
            abs(_norm(axis)-1)<=tolerance || throw(ArgumentError(
                "PlaneFrame: $name must have unit length"))
        end
        (abs(_dot(u,v))<=tolerance&&abs(_dot(u,n))<=tolerance&&
         abs(_dot(v,n))<=tolerance) || throw(ArgumentError(
            "PlaneFrame: u, v, and normal must be mutually orthogonal"))
        _norm(_sub(_cross(u,v),n))<=tolerance || throw(ArgumentError(
            "PlaneFrame: u × v must agree with normal"))
        new(o,u,v,n)
    end
end

function PlaneFrame(origin,u,v,normal)
    return PlaneFrame(_point3(origin,"PlaneFrame origin"),
                      _point3(u,"PlaneFrame u"),
                      _point3(v,"PlaneFrame v"),
                      _point3(normal,"PlaneFrame normal"))
end

function _loop_points(loop,caller::AbstractString,max_points::Integer=typemax(Int32))
    count=try
        length(loop)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: each loop must be an indexable collection of points"))
    end
    count>=3 || throw(ArgumentError("$caller: a loop needs ≥3 points"))
    count<=max_points || throw(ArgumentError(
        "$caller: boundary point count exceeds Int32 indexing"))
    points=Vector{NTuple{3,Float64}}(undef,count)
    seen=0
    for raw in loop
        seen+=1
        seen<=count || throw(ArgumentError(
            "$caller: loop iteration produced more points than its reported length"))
        points[seen]=_point3(raw,caller)
    end
    seen==count || throw(ArgumentError(
        "$caller: loop iteration produced fewer points than its reported length"))
    @inbounds for index in 1:count
        points[index]!=points[mod1(index+1,count)] || throw(ArgumentError(
            "$caller: loop edge $index has coincident endpoints"))
    end
    return points
end

function _normalized_points(points,caller::AbstractString)
    coordinate_scale=maximum(abs,(point[dimension] for point in points for dimension in 1:3))
    coordinate_scale>0 || throw(ArgumentError("$caller: loop is geometrically degenerate"))
    scaled=map(point->ntuple(d->point[d]/coordinate_scale,3),points)
    origin=scaled[1];span=0.0
    @inbounds for point in scaled,dimension in 1:3
        span=max(span,abs(point[dimension]-origin[dimension]))
    end
    (isfinite(span)&&span>0) || throw(ArgumentError(
        "$caller: loop has no representable coordinate span"))
    normalized=map(point->ntuple(d->(point[d]-origin[d])/span,3),scaled)
    all(point->all(isfinite,point),normalized) || throw(ArgumentError(
        "$caller: normalized loop coordinates are not finite"))
    return normalized
end

function _exact_plane_normal(points)
    R=Rational{BigInt};origin=ntuple(d->R(points[1][d]),3)
    offsets=map(points) do point
        ntuple(d->R(point[d])-origin[d],3)
    end
    nx=zero(R);ny=zero(R);nz=zero(R);count=length(offsets)
    @inbounds for index in 1:count
        point=offsets[index];next=offsets[mod1(index+1,count)]
        nx+=(point[2]-next[2])*(point[3]+next[3])
        ny+=(point[3]-next[3])*(point[1]+next[1])
        nz+=(point[1]-next[1])*(point[2]+next[2])
    end
    exact_normal=(nx,ny,nz);normal_squared=_dot(exact_normal,exact_normal)
    normal_squared>0 || return nothing
    all(point->_dot(exact_normal,point)==0,offsets) || return nothing
    normal=setprecision(BigFloat,256) do
        denominator=sqrt(BigFloat(normal_squared))
        ntuple(d->Float64(BigFloat(exact_normal[d])/denominator),3)
    end
    all(isfinite,normal)&&_norm(normal)>0 || return nothing
    return normal
end

function _frame_from_normal(origin,normal)
    ax = abs(normal[1]) <= abs(normal[2]) ?
         (abs(normal[1]) <= abs(normal[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(normal[2]) <= abs(normal[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    u=_unit(_cross(ax,normal));v=_cross(normal,u)
    return PlaneFrame(origin,u,v,normal)
end

function _plane_frame(points,caller::AbstractString)
    normalized=_normalized_points(points,caller);m=length(points)
    nx=ny=nz=0.0;permanent_x=permanent_y=permanent_z=0.0
    @inbounds for i in 1:m
        p=normalized[i];q=normalized[mod1(i+1,m)]
        term_x=(p[2]-q[2])*(p[3]+q[3])
        term_y=(p[3]-q[3])*(p[1]+q[1])
        term_z=(p[1]-q[1])*(p[2]+q[2])
        nx+=term_x;ny+=term_y;nz+=term_z
        permanent_x+=abs(term_x);permanent_y+=abs(term_y);permanent_z+=abs(term_z)
    end
    length_normal=hypot(nx,ny,nz)
    reliable=abs(nx)>64eps(Float64)*permanent_x||
             abs(ny)>64eps(Float64)*permanent_y||
             abs(nz)>64eps(Float64)*permanent_z
    if reliable&&isfinite(length_normal)&&length_normal>0
        normal=(nx/length_normal,ny/length_normal,nz/length_normal)
        if all(point->abs(_dot(point,normal))<=256eps(Float64),normalized)
            return _frame_from_normal(points[1],normal)
        end
    end
    normal=_exact_plane_normal(points)
    normal===nothing && throw(ArgumentError(
        "$caller: loop is degenerate or is not coplanar"))
    return _frame_from_normal(points[1],normal)
end

function _validate_coplanar_loops(loops,fr::PlaneFrame)
    coordinate_scale=max(abs(fr.o[1]),abs(fr.o[2]),abs(fr.o[3]),
        maximum(abs,(point[dimension] for loop in loops for point in loop for dimension in 1:3)))
    coordinate_scale>0 || throw(ArgumentError(
        "mesh_planar_face: boundary is geometrically degenerate"))
    scaled_origin=ntuple(d->fr.o[d]/coordinate_scale,3);span=0.0
    @inbounds for loop in loops,point in loop,dimension in 1:3
        span=max(span,abs(point[dimension]/coordinate_scale-scaled_origin[dimension]))
    end
    (isfinite(span)&&span>0) || throw(ArgumentError(
        "mesh_planar_face: boundary has no representable coordinate span"))
    tolerance=256eps(Float64)
    @inbounds for (loop_index,loop) in pairs(loops),(point_index,point) in pairs(loop)
        offset=ntuple(d->(point[d]/coordinate_scale-scaled_origin[d])/span,3)
        distance=abs(_dot(offset,fr.n))
        (isfinite(distance)&&distance<=tolerance) || throw(ArgumentError(
            "mesh_planar_face: loop $loop_index point $point_index is not coplanar " *
            "(normalized plane distance $distance exceeds $tolerance)"))
    end
    return nothing
end

"""
    plane_frame(loop) -> PlaneFrame

Construct an orthonormal frame from a finite, nondegenerate coplanar loop.
Coordinates are normalized before Newell evaluation; exact rational fallback
handles representable cancellation cases.
"""
plane_frame(loop)=_plane_frame(_loop_points(loop,"plane_frame"),"plane_frame")

function _project_component(point,origin,axis,caller)
    delta=_sub(point,origin);value=_dot(delta,axis)
    products=ntuple(d->delta[d]*axis[d],3)
    permanent=sum(abs,products)
    if isfinite(value)&&all(isfinite,delta)&&isfinite(permanent)&&
       (permanent==0||abs(value)>32eps(Float64)*permanent)
        return value==0 ? 0.0 : value
    end
    exact=setprecision(BigFloat,256) do
        sum((BigFloat(point[d])-BigFloat(origin[d]))*BigFloat(axis[d]) for d in 1:3)
    end
    result=Float64(exact)
    isfinite(result) || throw(ArgumentError("$caller: projected coordinate is not representable"))
    return result==0 ? 0.0 : result
end

"""Project a finite 3-D point into the two in-plane coordinates of `fr`."""
function project(fr::PlaneFrame,p)
    point=_point3(p,"project")
    return (_project_component(point,fr.o,fr.u,"project"),
            _project_component(point,fr.o,fr.v,"project"))
end

function _affine_coordinate(origin,a,u,b,v,caller)
    value=muladd(a,u,muladd(b,v,origin))
    first=a*u;second=b*v;permanent=abs(origin)+abs(first)+abs(second)
    if isfinite(value)&&isfinite(permanent)&&
       (permanent==0||abs(value)>32eps(Float64)*permanent)
        return value==0 ? 0.0 : value
    end
    exact=setprecision(BigFloat,256) do
        BigFloat(origin)+BigFloat(a)*BigFloat(u)+BigFloat(b)*BigFloat(v)
    end
    result=Float64(exact)
    isfinite(result) || throw(ArgumentError("$caller: output coordinate is not representable"))
    return result==0 ? 0.0 : result
end

function _affine3_coordinate(origin,a,u,b,v,c,w,caller)
    value=muladd(a,u,muladd(b,v,muladd(c,w,origin)))
    first=a*u;second=b*v;third=c*w
    permanent=abs(origin)+abs(first)+abs(second)+abs(third)
    if isfinite(value)&&isfinite(permanent)&&
       (permanent==0||abs(value)>32eps(Float64)*permanent)
        return value==0 ? 0.0 : value
    end
    exact=setprecision(BigFloat,256) do
        BigFloat(origin)+BigFloat(a)*BigFloat(u)+BigFloat(b)*BigFloat(v)+
        BigFloat(c)*BigFloat(w)
    end
    result=Float64(exact)
    isfinite(result) || throw(ArgumentError("$caller: output coordinate is not representable"))
    return result==0 ? 0.0 : result
end

"""Lift finite in-plane coordinates `(a,b)` through `fr` into 3-D."""
function lift(fr::PlaneFrame,a::Real,b::Real)
    aa=_surface_float(a,"lift","a");bb=_surface_float(b,"lift","b")
    (isfinite(aa)&&isfinite(bb)) || throw(ArgumentError(
        "lift: in-plane coordinates must be finite"))
    return ntuple(d->_affine_coordinate(fr.o[d],aa,fr.u[d],bb,fr.v[d],"lift"),3)
end

# ── planar face meshing ─────────────────────────────────────────────────────────
"""
    mesh_planar_face(loops, sf; min_angle_deg=25.0, max_area=Inf,
                     entity=nothing, boundary_entities=nothing) -> Mesh

Mesh a planar face bounded by `loops` (a vector of loops; `loops[1]` is the outer
boundary, the rest are holes). Each loop is a vector of coplanar 3-D points in
order. Boundary edges are size-graded, then the interior is constrained-Delaunay
meshed and Ruppert-refined under `sf`, and lifted back to 3-D.

`entity` is the face context used for interior field samples. `boundary_entities`
optionally supplies the actual curve context for each boundary edge, either as one
entity tuple, one tuple per loop, one tuple per edge within each loop, or a callable
`(loop_index, edge_index, p, q) -> entity`. This distinction is required by Gmsh
`Restrict` and `Constant`: a face does not implicitly classify its boundary curves
when `IncludeBoundary` is false.
"""
function mesh_planar_face(loops::AbstractVector, sf::AbstractSizeField;
                          min_angle_deg::Real=25.0, max_area::Real=Inf,
                          entity=nothing,boundary_entities=nothing)
    isempty(loops) && throw(ArgumentError("mesh_planar_face: no loops"))
    angle=_surface_float(min_angle_deg,"mesh_planar_face","min_angle_deg")
    (isfinite(angle)&&0<=angle<60) || throw(ArgumentError(
        "mesh_planar_face: min_angle_deg must be finite and in [0, 60)"))
    area=_surface_float(max_area,"mesh_planar_face","max_area")
    (!isnan(area)&&area>0) || throw(ArgumentError(
        "mesh_planar_face: max_area must be positive or Inf"))
    face_entity=_surface_entity_context(entity,"mesh_planar_face")
    length(loops)<=typemax(Int32) || throw(ArgumentError(
        "mesh_planar_face: loop count exceeds Int32 indexing"))
    checked_loops=Vector{Vector{NTuple{3,Float64}}}(undef,length(loops))
    point_count=0
    for (index,loop) in enumerate(loops)
        checked_loops[index]=_loop_points(
            loop,"mesh_planar_face",typemax(Int32)-point_count)
        point_count+=length(checked_loops[index])
    end
    fr=_plane_frame(checked_loops[1],"mesh_planar_face")
    _validate_coplanar_loops(checked_loops,fr)
    xs = Float64[]; ys = Float64[]; segs = Tuple{Int,Int}[]
    for (loop_index,loop) in pairs(checked_loops)
        _add_loop!(xs,ys,segs,loop,sf,fr,face_entity,boundary_entities,
                   length(checked_loops),Int(loop_index))
    end
    T = constrained_delaunay(xs, ys, segs)
    sizefn = sf isa AbstractAnisoField ? nothing :
             (a, b) -> size_at(sf, lift(fr, a, b)...,face_entity)
    edgefn=(ax,ay,bx,by)->metric_edge_length(
        sf,lift(fr,ax,ay),lift(fr,bx,by);entity=face_entity)
    interior = refine!(T; min_angle_deg=angle, max_area=area,
                       size=sizefn,edge_metric=edgefn)
    m2 = to_mesh(T; interior=interior)
    return _lift_mesh(m2, fr)
end

# discretize one loop's edges under sf and append its cyclic boundary as segments
function _add_loop!(xs,ys,segs,loop,sf,fr,entity,boundary_entities,nloops,loop_index)
    m = length(loop)
    m >= 3 || throw(ArgumentError("mesh_planar_face: a loop needs ≥3 points"))
    start = length(xs) + 1
    for i in 1:m
        p = _point3(loop[i], "mesh_planar_face")
        q = _point3(loop[mod1(i+1,m)], "mesh_planar_face")
        edge_entity=_surface_boundary_entity(boundary_entities,nloops,loop_index,
            i,m,p,q,entity,"mesh_planar_face")
        pts, _ = mesh_segment(p,q,sf;entity=edge_entity)
        # append all but the last node (shared with the next edge / loop closure)
        length(pts)-1 <= typemax(Int32)-length(xs) ||
            throw(ArgumentError("mesh_planar_face: boundary node count exceeds Int32 indexing"))
        @inbounds for k in 1:length(pts)-1
            a, b = project(fr, pts[k]); push!(xs, a); push!(ys, b)
        end
    end
    stop = length(xs)
    # cyclic segments over this loop's appended nodes
    for i in start:stop
        j = i == stop ? start : i+1
        push!(segs, (i, j))
    end
    return nothing
end

function _lift_mesh(m2::Mesh, fr::PlaneFrame)
    nn = nnodes(m2)
    coords = Matrix{Float64}(undef, 3, nn)
    @inbounds for i in 1:nn
        a = m2.coords[1,i]; b = m2.coords[2,i]
        p = lift(fr, a, b)
        coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]
    end
    out=Mesh(coords; tris=copy(m2.tris), tri_tag=copy(m2.tri_tag))
    d=validate(out)
    d.ok || throw(ErrorException("mesh_planar_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    return out
end

# ── developable cylinder face ───────────────────────────────────────────────────
"""
    mesh_cylinder_face(center, axis, radius, height, sf;
                       entity=nothing, boundary_entities=nothing,
                       max_refine_passes=16, kwargs...) -> Mesh

Mesh the lateral surface of a cylinder: axis start `center`, direction `axis`,
`radius`, `height`. Built as a **graded structured** grid — axial levels are
size-graded (via [`Mesh1D`](@ref) along the axis), the circumference is divided
uniformly into `nθ = max(6, ceil(2πR/h))` sectors, and each quad is split into
two triangles with the seam wrapping (θ_{nθ} ≡ θ_0). Watertight by construction
(no seam to weld), which a Ruppert mesh of the unrolled rectangle cannot
guarantee once the two seam edges refine asymmetrically.

`entity` classifies lateral-face and axial sizing queries. `boundary_entities`
optionally supplies the actual curve contexts for the bottom and top rings,
respectively, as one entity tuple, a two-entry vector, or a callable
`(ring_index, edge_index, p, q) -> entity`. Since the structured surface uses a
single circumferential division count, either ring may conservatively refine
both rings.

The structured grid is rechecked against `sf` at every horizontal, vertical,
and diagonal edge endpoint and midpoint. Violating intervals are subdivided
until the metric-edge bound is certified or `max_refine_passes` is exhausted.
"""
function mesh_cylinder_face(center, axis, radius::Real, height::Real,
                            sf::AbstractSizeField; min_angle_deg::Real=25.0,
                            max_area::Real=Inf, entity=nothing,
                            boundary_entities=nothing,
                            max_refine_passes::Integer=16)
    R = _surface_float(radius,"mesh_cylinder_face","radius")
    H = _surface_float(height,"mesh_cylinder_face","height")
    (isfinite(R) && isfinite(H) && R > 0 && H > 0) ||
        throw(ArgumentError("mesh_cylinder_face: radius and height must be finite and positive (got $radius, $height)"))
    angle=_surface_float(min_angle_deg,"mesh_cylinder_face","min_angle_deg")
    (0<=angle<=atand(inv(sqrt(2.0)))) ||
        throw(ArgumentError("mesh_cylinder_face: min_angle_deg must be in [0, $(atand(inv(sqrt(2.0))))] for the structured right-triangle grid"))
    area=_surface_float(max_area,"mesh_cylinder_face","max_area")
    (!isnan(area)&&area>0) ||
        throw(ArgumentError("mesh_cylinder_face: max_area must be positive or Inf"))
    max_refine_passes isa Bool && throw(ArgumentError(
        "mesh_cylinder_face: max_refine_passes must be a positive integer"))
    passes=try
        Int(max_refine_passes)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_cylinder_face: max_refine_passes exceeds the platform Int limit"))
    end
    passes>0 || throw(ArgumentError(
        "mesh_cylinder_face: max_refine_passes must be positive"))
    c = _point3(center, "mesh_cylinder_face")
    ez = _unit(_point3(axis, "mesh_cylinder_face"))
    face_entity=_surface_entity_context(entity,"mesh_cylinder_face")
    ax = abs(ez[1]) <= abs(ez[2]) ?
         (abs(ez[1]) <= abs(ez[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(ez[2]) <= abs(ez[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    ex = _unit(_cross(ax, ez)); ey = _cross(ez, ex)
    on(θ,z)=begin
        radial_u=R*cos(θ);radial_v=R*sin(θ)
        ntuple(d->_affine3_coordinate(c[d],radial_u,ex[d],radial_v,ey[d],
                                      z,ez[d],"mesh_cylinder_face"),3)
    end
    bottom_entity=_surface_boundary_entity(boundary_entities,2,1,1,1,
        on(0.0,0.0),on(π,0.0),face_entity,"mesh_cylinder_face")
    top_entity=_surface_boundary_entity(boundary_entities,2,2,1,1,
        on(0.0,H),on(π,H),face_entity,"mesh_cylinder_face")
    # axial levels: size-graded along the axis (uses the real 3-D size field)
    axis_end=ntuple(d->_affine_coordinate(c[d],H,ez[d],0.0,0.0,
                                          "mesh_cylinder_face"),3)
    apts, _ = mesh_segment(c,axis_end,sf;
                           entity=face_entity)
    zlev=Float64[_project_component(
        _point3(p,"mesh_cylinder_face axial grading"),c,ez,
        "mesh_cylinder_face axial grading") for p in apts] # 0 .. H, graded
    length(zlev)>=2 || throw(ErrorException(
        "mesh_cylinder_face: axial grading returned fewer than two levels"))
    nz = length(zlev)
    # Circumferential divisions use the smallest sampled target over all axial
    # levels, not only the equator (an axial size field must also refine rings).
    hmin = Inf
    for (j,z) in pairs(zlev),k in 0:11
        θ=2π*k/12
        tangent=(-sin(θ)*ex[1]+cos(θ)*ey[1],
                 -sin(θ)*ex[2]+cos(θ)*ey[2],
                 -sin(θ)*ex[3]+cos(θ)*ey[3])
        ring_entity=j==1 ? bottom_entity : j==length(zlev) ? top_entity : face_entity
        hmin = min(hmin, directional_size(sf,on(θ,z),tangent;
                                          entity=ring_entity))
    end
    mindz=minimum(zlev[i+1]-zlev[i] for i in 1:length(zlev)-1)
    (isfinite(mindz)&&mindz>0) ||
        throw(ArgumentError("mesh_cylinder_face: axial grading produced a non-positive interval"))
    tana=tand(angle)
    chord_target=hmin
    if isfinite(area)
        chord_target=min(chord_target,angle==0 ? sqrt(2area) : sqrt(area/tana))
    end
    angle>0 && (chord_target=min(chord_target,mindz/tana))
    (isfinite(chord_target)&&chord_target>0) ||
        throw(ArgumentError("mesh_cylinder_face: requested controls are below Float64 resolution"))
    sectors = (2π)*(R/chord_target)
    sectors <= prevfloat(Float64(typemax(Int))) ||
        throw(ArgumentError("mesh_cylinder_face: requested circumferential node count exceeds the platform Int limit"))
    ntheta = max(6, ceil(Int, sectors))
    chord=R*(2sinpi(inv(ntheta)))
    (isfinite(chord)&&chord>0) ||
        throw(ArgumentError("mesh_cylinder_face: circumferential chord is below Float64 resolution"))

    # Subdivide every axial interval enough to satisfy the requested right-triangle
    # angle and area controls.  The supported angle ceiling makes adjacent feasible
    # subdivision intervals overlap, so this deterministic ceil rule has no gaps.
    refined_z=Float64[zlev[1]]
    for i in 1:length(zlev)-1
        dz=zlev[i+1]-zlev[i]
        maxdz=angle==0 ? Inf : chord/tana
        isfinite(area) && (maxdz=min(maxdz,2area/chord))
        ratio=dz/maxdz
        (!isfinite(maxdz) || (isfinite(ratio)&&ratio<=prevfloat(Float64(typemax(Int))))) ||
            throw(ArgumentError("mesh_cylinder_face: axial subdivision count exceeds the platform Int limit"))
        nsub=isfinite(maxdz) ? ceil(Int,ratio) : 1
        nsub=max(1,nsub);piece=dz/nsub
        angle>0 && piece < chord*tana*(1-128eps(Float64)) &&
            throw(ErrorException("mesh_cylinder_face: could not satisfy the requested angle bound at Float64 resolution"))
        for k in 1:nsub
            candidate=k==nsub ? zlev[i+1] :
                _param_interp(zlev[i],zlev[i+1],k/nsub)
            candidate>refined_z[end] || throw(ErrorException(
                "mesh_cylinder_face: axial subdivision is below Float64 resolution"))
            push!(refined_z,candidate)
        end
    end
    zlev=refined_z;nz=length(zlev)
    zlev,ntheta=_refine_cylinder_metric_grid(sf,on,zlev,ntheta,R,angle,
        face_entity,bottom_entity,top_entity,passes)
    nz=length(zlev)
    nnout=try
        Base.checked_mul(nz, ntheta)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_cylinder_face: requested mesh dimensions overflow the platform Int limit"))
    end
    nnout<=typemax(Int32) ||
        throw(ArgumentError("mesh_cylinder_face: $nnout nodes exceed the Int32 indexing limit"))
    ntout=try
        Base.checked_mul(Base.checked_mul(2, nz-1), ntheta)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_cylinder_face: requested mesh dimensions overflow the platform Int limit"))
    end
    ntout<=typemax(Int32) ||
        throw(ArgumentError("mesh_cylinder_face: $ntout triangles exceed the Int32 topology limit"))
    # structured nodes: level j (1..nz) × sector i (0..ntheta-1)
    coords = Matrix{Float64}(undef, 3, nnout)
    idx(j, i) = (j-1)*ntheta + (mod(i, ntheta)) + 1
    @inbounds for j in 1:nz, i in 0:ntheta-1
        p = on(2π*i/ntheta, zlev[j]); n = idx(j, i)
        coords[1,n]=p[1]; coords[2,n]=p[2]; coords[3,n]=p[3]
    end
    tris = Matrix{Int32}(undef, 3, ntout); nt = 0
    @inbounds for j in 1:nz-1, i in 0:ntheta-1
        a = idx(j,i); b = idx(j,i+1); cc = idx(j+1,i+1); d = idx(j+1,i)
        nt+=1; tris[1,nt]=a; tris[2,nt]=b; tris[3,nt]=cc
        nt+=1; tris[1,nt]=a; tris[2,nt]=cc; tris[3,nt]=d
    end
    out=Mesh(coords; tris=tris[:, 1:nt])
    d=validate(out)
    d.ok || throw(ErrorException("mesh_cylinder_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    @inbounds for t in 1:ntris(out)
        ia=out.tris[1,t];ib=out.tris[2,t];ic=out.tris[3,t]
        a=node(out,ia);b=node(out,ib);c3=node(out,ic)
        A=triangle_area(a,b,c3);amin=_triangle_min_angle_deg(a,b,c3)
        A<=area*(1+256eps(Float64)) ||
            throw(ErrorException("mesh_cylinder_face: postcondition failed: triangle $t exceeds max_area"))
        amin+8192eps(Float64)*max(1.0,angle)>=angle ||
            throw(ErrorException("mesh_cylinder_face: postcondition failed: triangle $t violates min_angle_deg"))
        for (first_id,last_id,first_point,last_point) in
            ((ia,ib,a,b),(ib,ic,b,c3),(ic,ia,c3,a))
            first_level=div(Int(first_id)-1,ntheta)+1
            last_level=div(Int(last_id)-1,ntheta)+1
            edge_entity=first_level==1&&last_level==1 ? bottom_entity :
                first_level==nz&&last_level==nz ? top_entity : face_entity
            metric_edge_length(sf,first_point,last_point;entity=edge_entity)<=
                1+4096eps(Float64) || throw(ErrorException(
                    "mesh_cylinder_face: postcondition failed: triangle $t has an edge that violates the size field"))
        end
    end
    return out
end

@inline function _cylinder_subdivision_factor(value::Float64)
    value<=1+256eps(Float64) && return 1
    value<=prevfloat(Float64(typemax(Int))) || throw(ArgumentError(
        "mesh_cylinder_face: field sizing requests more subdivisions than the platform Int limit"))
    return ceil(Int,value)
end

function _refine_cylinder_metric_grid(sf,on,zlev::Vector{Float64},ntheta::Int,
                                      radius::Float64,angle::Float64,
                                      face_entity,bottom_entity,top_entity,
                                      maxpasses::Int)
    for pass in 1:maxpasses
        nz=length(zlev)
        axial_factors=ones(Int,nz-1)
        angular_factor=1
        @inbounds for j in 1:nz
            z=zlev[j]
            ring_entity=j==1 ? bottom_entity : j==nz ? top_entity : face_entity
            for i in 0:ntheta-1
                a=on(2π*i/ntheta,z);b=on(2π*(i+1)/ntheta,z)
                angular_factor=max(angular_factor,_cylinder_subdivision_factor(
                    metric_edge_length(sf,a,b;entity=ring_entity)))
            end
        end
        @inbounds for j in 1:nz-1, i in 0:ntheta-1
            a=on(2π*i/ntheta,zlev[j])
            b=on(2π*(i+1)/ntheta,zlev[j])
            c=on(2π*(i+1)/ntheta,zlev[j+1])
            d=on(2π*i/ntheta,zlev[j+1])
            vertical=max(metric_edge_length(sf,a,d;entity=face_entity),
                         metric_edge_length(sf,b,c;entity=face_entity))
            diagonal=metric_edge_length(sf,a,c;entity=face_entity)
            axial_factors[j]=max(axial_factors[j],
                                 _cylinder_subdivision_factor(max(vertical,diagonal)))
            angular_factor=max(angular_factor,
                               _cylinder_subdivision_factor(diagonal))
        end
        next_ntheta=try
            Base.checked_mul(ntheta,angular_factor)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "mesh_cylinder_face: field-refined circumferential count overflows the platform Int limit"))
        end

        # Field-driven subdivisions can change the structured cell aspect ratio.
        # Couple the two directions so the requested physical-space angle bound
        # remains true after every refinement pass.
        if angle>0
            tana=tand(angle)
            chord=2radius*sinpi(inv(ntheta))
            @inbounds for j in eachindex(axial_factors)
                dz=zlev[j+1]-zlev[j]
                axial_factors[j]=max(axial_factors[j],
                    _cylinder_subdivision_factor(
                        (dz*tana/chord)*(1+2048eps(Float64))))
            end
            @inbounds for j in eachindex(axial_factors)
                piece=(zlev[j+1]-zlev[j])/axial_factors[j]
                next_ntheta=_cylinder_sectors_for_chord(
                    radius,(piece/tana)*(1-2048eps(Float64)),next_ntheta)
            end
        end

        next_ntheta==ntheta && all(==(1),axial_factors) && return zlev,ntheta
        pass<maxpasses || throw(ErrorException(
            "mesh_cylinder_face: field metric did not converge within max_refine_passes=$maxpasses"))
        ntheta=next_ntheta
        ntheta<=typemax(Int32) || throw(ArgumentError(
            "mesh_cylinder_face: field-refined circumferential count exceeds Int32 indexing"))
        refined_count=1
        for factor in axial_factors
            refined_count=try
                Base.checked_add(refined_count,factor)
            catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError(
                    "mesh_cylinder_face: field-refined axial count overflows the platform Int limit"))
            end
        end
        refined_count<=typemax(Int32) || throw(ArgumentError(
            "mesh_cylinder_face: field-refined axial count exceeds Int32 indexing"))
        try
            Base.checked_mul(refined_count,ntheta)<=typemax(Int32) ||
                throw(ArgumentError(
                    "mesh_cylinder_face: field-refined node count exceeds Int32 indexing"))
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError && rethrow()
            throw(ArgumentError(
                "mesh_cylinder_face: field-refined node count overflows the platform Int limit"))
        end
        refined=Vector{Float64}(undef,refined_count);out=1;refined[1]=zlev[1]
        @inbounds for j in eachindex(axial_factors)
            n=axial_factors[j];z0=zlev[j];z1=zlev[j+1]
            for k in 1:n
                out+=1
                candidate=k==n ? z1 : _param_interp(z0,z1,k/n)
                candidate>refined[out-1] || throw(ErrorException(
                    "mesh_cylinder_face: field-refined axial subdivision is below Float64 resolution"))
                refined[out]=candidate
            end
        end
        zlev=refined
    end
    throw(ErrorException("mesh_cylinder_face: unreachable metric-refinement state"))
end

function _cylinder_sectors_for_chord(radius::Float64,target::Float64,
                                     current::Int)
    (isfinite(target)&&target>0) || throw(ArgumentError(
        "mesh_cylinder_face: requested angle bound is below Float64 resolution"))
    2radius*sinpi(inv(current))<=target && return current
    ratio=target/(2radius)
    ratio>0 || throw(ArgumentError(
        "mesh_cylinder_face: requested circumferential count exceeds the platform Int limit"))
    estimate=π/asin(min(ratio,1.0))
    estimate<=prevfloat(Float64(typemax(Int))) || throw(ArgumentError(
        "mesh_cylinder_face: requested circumferential count exceeds the platform Int limit"))
    required=max(current,ceil(Int,estimate))
    while 2radius*sinpi(inv(required))>target
        required<typemax(Int) || throw(ArgumentError(
            "mesh_cylinder_face: requested circumferential count exceeds the platform Int limit"))
        required+=1
    end
    return required
end

function _normalized_triangle(a,b,c)
    points=(a,b,c);edge_scale=0.0;finite_edges=true
    @inbounds for i in 1:2,j in i+1:3,dimension in 1:3
        difference=points[j][dimension]-points[i][dimension]
        if isfinite(difference)
            edge_scale=max(edge_scale,abs(difference))
        else
            finite_edges=false
        end
    end
    scaled=points
    if !(finite_edges&&(edge_scale==0||isfinite(inv(edge_scale))))
        coordinate_scale=maximum(abs,(a...,b...,c...))
        coordinate_scale==0 && return nothing
        scaled=ntuple(j->ntuple(d->points[j][d]/coordinate_scale,3),3)
    end
    anchor=scaled[1]
    edge_scale=maximum(abs,(_sub(scaled[j],anchor)[d] for j in 2:3 for d in 1:3))
    edge_scale==0 && return nothing
    return ((0.0,0.0,0.0),
        ntuple(d->(scaled[2][d]-anchor[d])/edge_scale,3),
        ntuple(d->(scaled[3][d]-anchor[d])/edge_scale,3))
end

function _triangle_min_angle_deg(a,b,c)
    normalized=_normalized_triangle(a,b,c)
    normalized===nothing && return 0.0
    sides=sort([_norm(_sub(normalized[1],normalized[2])),
                _norm(_sub(normalized[2],normalized[3])),
                _norm(_sub(normalized[3],normalized[1]))])
    sides[1]>0 || return 0.0
    q=sides./sides[3]
    return acosd(clamp((q[2]^2+1-q[1]^2)/(2q[2]),-1.0,1.0))
end

function _triangle_unit_normal(a,b,c)
    normalized=_normalized_triangle(a,b,c)
    normalized===nothing && return nothing
    normal=_cross(_sub(normalized[2],normalized[1]),
                  _sub(normalized[3],normalized[1]))
    length_normal=_norm(normal)
    (isfinite(length_normal)&&length_normal>0) || return nothing
    return ntuple(d->normal[d]/length_normal,3)
end

function _surface_difference_interval(value,lower,upper,step)
    left=max(lower,value-step);right=min(upper,value+step)
    if left==value&&value>lower
        left=max(lower,prevfloat(value))
    end
    if right==value&&value<upper
        right=min(upper,nextfloat(value))
    end
    right>left || throw(ArgumentError(
        "mesh_parametric_face: finite-difference interval is below parameter resolution"))
    return left,right
end

# ── general parametric face (isotropic-metric approximation) ─────────────────────
"""
    mesh_parametric_face(s, umin, umax, vmin, vmax, sf;
                         min_angle_deg=25.0, max_area=Inf,
                         entity=nothing, boundary_entities=nothing) -> Mesh

Mesh a parametric surface `s(u,v) -> (x,y,z)` over the rectangle
`[umin,umax]×[vmin,vmax]`. Isotropic refinement uses a local-stretch estimate,
and every final chord edge is independently checked against `sf`. `max_area` is
a physical 3-D triangle-area bound, enforced through a conservative chord-length
bound and a postcondition. `min_angle_deg` applies in parameter space; a general
map can distort physical angles. Suitable for regular, gently curved patches;
developable cylinders have the stronger exact route [`mesh_cylinder_face`](@ref).
"""
function mesh_parametric_face(s, umin::Real, umax::Real, vmin::Real, vmax::Real,
                              sf::AbstractSizeField; min_angle_deg::Real=25.0,
                              max_area::Real=Inf, entity=nothing,
                              boundary_entities=nothing)
    face_entity=_surface_entity_context(entity,"mesh_parametric_face")
    angle=_surface_float(min_angle_deg,"mesh_parametric_face","min_angle_deg")
    (isfinite(angle)&&0<=angle<60) || throw(ArgumentError(
        "mesh_parametric_face: min_angle_deg must be finite and in [0, 60)"))
    umin=_surface_float(umin,"mesh_parametric_face","umin")
    umax=_surface_float(umax,"mesh_parametric_face","umax")
    vmin=_surface_float(vmin,"mesh_parametric_face","vmin")
    vmax=_surface_float(vmax,"mesh_parametric_face","vmax")
    (isfinite(umin) && isfinite(umax) && isfinite(vmin) && isfinite(vmax) &&
     umin < umax && vmin < vmax) ||
        throw(ArgumentError("mesh_parametric_face: require finite bounds umin < umax and vmin < vmax"))
    uspan=umax-umin;vspan=vmax-vmin
    (isfinite(uspan)&&isfinite(vspan)&&uspan>0&&vspan>0) || throw(ArgumentError(
        "mesh_parametric_face: parameter spans are not representable"))
    du=1e-6uspan;dv=1e-6vspan
    du>0 || (du=uspan);dv>0 || (dv=vspan)
    applicable(s,umin,vmin) || throw(ArgumentError(
        "mesh_parametric_face: s must be callable as s(u,v)"))
    area=_surface_float(max_area,"mesh_parametric_face","max_area")
    (!isnan(area)&&area>0) || throw(ArgumentError(
        "mesh_parametric_face: max_area must be positive or Inf"))
    area_edge=isfinite(area) ? 2sqrt(area)/sqrt(sqrt(3.0)) : Inf
    (isfinite(area_edge)&&area_edge>0)||!isfinite(area) || throw(ArgumentError(
        "mesh_parametric_face: max_area is below Float64 edge-length resolution"))
    surf(u,v) = _point3(s(u,v), "mesh_parametric_face surface callback")
    # Local coordinate directions and stretch factors (√E, √G) from finite
    # differences. The unit directions also certify mapped triangle orientation.
    surface_jacobian(u,v)=begin
        ul,ur=_surface_difference_interval(u,umin,umax,du)
        vl,vr=_surface_difference_interval(v,vmin,vmax,dv)
        su = _sub(surf(ur,v), surf(ul,v)); sv = _sub(surf(u,vr), surf(u,vl))
        lu=_norm(su);lv=_norm(sv);a=lu/(ur-ul);b=lv/(vr-vl)
        uu=(su[1]/lu,su[2]/lu,su[3]/lu);vv=(sv[1]/lv,sv[2]/lv,sv[3]/lv)
        sinang=_norm(_cross(uu,vv))
        (isfinite(a) && isfinite(b) && isfinite(sinang) && a>0 && b>0 &&
         sinang > 256eps(Float64)) ||
            throw(ArgumentError("mesh_parametric_face: singular or non-finite surface Jacobian at ($u,$v)"))
        (a,b,uu,vv)
    end
    stretch(u,v)=surface_jacobian(u,v)[1:2]
    # Isotropic fields use the local-stretch estimate to select candidates. Every
    # field type is also checked directly on physical chord edges below.
    sizefn = sf isa AbstractAnisoField ? nothing : (u,v) -> begin
        p = surf(u,v); h = size_at(sf,p...,face_entity); (a,b) = stretch(u,v)
        hp=h/max(a,b)
        (isfinite(hp)&&hp>0) || throw(ArgumentError(
            "mesh_parametric_face: parameter-space size is not representable at ($u,$v)"))
        hp
    end
    # PSLG = parameter rectangle boundary, certified in physical chord space.
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    corners=[(umin,vmin),(umax,vmin),(umax,vmax),(umin,vmax)]
    edge_entities=Vector{Union{Nothing,Tuple{Int,Int}}}(undef,4)
    start=1
    for i in 1:4
        p=corners[i]; q=corners[mod1(i+1,4)]
        edge_entity=_surface_boundary_entity(boundary_entities,1,1,i,4,p,q,
            face_entity,"mesh_parametric_face")
        edge_entities[i]=edge_entity
        parameters=_param_boundary_parameters(surf,sf,p,q,edge_entity,area_edge)
        length(parameters)-1<=typemax(Int32)-length(xs) ||
            throw(ArgumentError("mesh_parametric_face: boundary node count exceeds Int32 indexing"))
        for parameter in @view parameters[1:end-1]
            point=_param_lerp(p,q,parameter)
            push!(xs,point[1]);push!(ys,point[2])
        end
    end
    stop=length(xs)
    for i in start:stop; push!(segs, (i, i==stop ? start : i+1)); end
    T=constrained_delaunay(xs,ys,segs)
    bounds=(umin,umax,vmin,vmax)
    edgefn=(u1,v1,u2,v2)->begin
        first=(u1,v1);last=(u2,v2)
        context=_param_surface_entity(first,last,bounds,edge_entities,face_entity)
        _param_edge_score(surf,sf,first,last,context,area_edge)
    end
    interior=refine!(T; min_angle_deg=angle, max_area=Inf,
                     size=sizefn,edge_metric=edgefn)
    m2=to_mesh(T; interior=interior)
    # map parameter mesh to 3-D
    nn=nnodes(m2); coords=Matrix{Float64}(undef,3,nn)
    @inbounds for i in 1:nn
        u=m2.coords[1,i];v=m2.coords[2,i]
        stretch(u,v)
        p=surf(u,v);coords[1,i]=p[1];coords[2,i]=p[2];coords[3,i]=p[3]
    end
    out=Mesh(coords; tris=copy(m2.tris))
    d=validate(out)
    d.ok || throw(ErrorException("mesh_parametric_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    tolerance=1+4096eps(Float64)
    @inbounds for triangle in 1:ntris(out)
        ids=(out.tris[1,triangle],out.tris[2,triangle],out.tris[3,triangle])
        points=(node(out,ids[1]),node(out,ids[2]),node(out,ids[3]))
        triangle_area(points...)<=area*tolerance || throw(ErrorException(
            "mesh_parametric_face: postcondition failed: triangle $triangle exceeds max_area"))
        for edge in ((1,2),(2,3),(3,1))
            first=(m2.coords[1,ids[edge[1]]],m2.coords[2,ids[edge[1]]])
            last=(m2.coords[1,ids[edge[2]]],m2.coords[2,ids[edge[2]]])
            context=_param_surface_entity(first,last,bounds,edge_entities,face_entity)
            metric_edge_length(sf,points[edge[1]],points[edge[2]];entity=context)<=tolerance ||
                throw(ErrorException(
                    "mesh_parametric_face: postcondition failed: triangle $triangle has an edge that violates the size field"))
        end
        u=m2.coords[1,ids[1]]/3+m2.coords[1,ids[2]]/3+m2.coords[1,ids[3]]/3
        v=m2.coords[2,ids[1]]/3+m2.coords[2,ids[2]]/3+m2.coords[2,ids[3]]/3
        _,_,direction_u,direction_v=surface_jacobian(u,v)
        mapped_normal=_triangle_unit_normal(points...)
        if mapped_normal===nothing ||
           _dot(mapped_normal,_cross(direction_u,direction_v))<=256eps(Float64)
            throw(ErrorException(
                "mesh_parametric_face: postcondition failed: triangle $triangle is inverted relative to the surface Jacobian"))
        end
    end
    return out
end

@inline _param_interp(a,b,t)=a==b ? a : muladd(t,b-a,a)
@inline _param_lerp(p,q,t)=(_param_interp(p[1],q[1],t),
                            _param_interp(p[2],q[2],t))

@inline function _param_surface_entity(p,q,bounds,entities,face_entity)
    umin,umax,vmin,vmax=bounds
    p[2]==vmin&&q[2]==vmin && return entities[1]
    p[1]==umax&&q[1]==umax && return entities[2]
    p[2]==vmax&&q[2]==vmax && return entities[3]
    p[1]==umin&&q[1]==umin && return entities[4]
    return face_entity
end

function _param_edge_metric(s,sf,p,q,entity)
    N = 32; L = 0.0; prev = s(p[1], p[2])
    @inbounds for k in 1:N
        point=_param_lerp(p,q,k/N);cur=s(point[1],point[2])
        L += _metric_curve_increment(sf,prev,cur,entity); prev = cur
        isfinite(L) || throw(ArgumentError("mesh_parametric_face: boundary metric length overflowed Float64"))
    end
    return L
end

function _param_edge_score(s,sf,p,q,entity,area_edge)
    a=s(p[1],p[2]);b=s(q[1],q[2])
    metric=metric_edge_length(sf,a,b;entity=entity)
    area_score=isfinite(area_edge) ? _norm(_sub(a,b))/area_edge : 0.0
    (isfinite(area_score)&&area_score>=0) || throw(ArgumentError(
        "mesh_parametric_face: physical edge length is not representable"))
    return max(metric,area_score)
end

function _param_boundary_parameters(s,sf,p,q,entity,area_edge)
    metric=_param_edge_metric(s,sf,p,q,entity)
    initial_score=max(metric,_param_edge_score(s,sf,p,q,entity,area_edge))
    (isfinite(initial_score)&&initial_score>=0&&
     initial_score<=prevfloat(Float64(typemax(Int)))) || throw(ArgumentError(
        "mesh_parametric_face: requested boundary node count exceeds the platform Int limit"))
    intervals=max(1,ceil(Int,initial_score))
    intervals<=typemax(Int32) || throw(ArgumentError(
        "mesh_parametric_face: requested boundary node count exceeds Int32 indexing"))
    parameters=collect(range(0.0,1.0;length=intervals+1))
    for pass in 1:64
        refined=Float64[parameters[1]];changed=false
        @inbounds for index in 1:length(parameters)-1
            first=parameters[index];last=parameters[index+1]
            first_point=_param_lerp(p,q,first);last_point=_param_lerp(p,q,last)
            if _param_edge_score(s,sf,first_point,last_point,entity,area_edge)>1
                midpoint=first/2+last/2
                first<midpoint<last || throw(ErrorException(
                    "mesh_parametric_face: boundary refinement is below parameter resolution"))
                length(refined)<typemax(Int32) || throw(ArgumentError(
                    "mesh_parametric_face: boundary node count exceeds Int32 indexing"))
                push!(refined,midpoint);changed=true
            end
            length(refined)<typemax(Int32) || throw(ArgumentError(
                "mesh_parametric_face: boundary node count exceeds Int32 indexing"))
            push!(refined,last)
        end
        changed || return refined
        pass<64 || throw(ErrorException(
            "mesh_parametric_face: boundary metric did not converge within 64 refinement passes"))
        parameters=refined
    end
    throw(ErrorException("mesh_parametric_face: unreachable boundary-refinement state"))
end

function _surface_float(x::Real,caller::AbstractString,name::AbstractString)
    x isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    try
        return Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be representable as Float64"))
    end
end

end # module MeshSurface
