"""
    MeshSurface

Stage-2 surface meshing (PLAN.md §3, "mesh the flat/patch/coax faces"). Two
routes, both reducing a surface to the exact 2-D machinery of [`Mesh2D`](@ref):

* **Planar faces** — a face bounded by coplanar loops (outer + holes). Boundary
  loops are first discretized under the size field ([`Mesh1D`](@ref)), projected
  to an in-plane orthonormal frame, meshed as a PSLG (constrained Delaunay +
  size-driven Ruppert refinement), then lifted back to 3-D. Isometric ⇒ the 2-D
  quality/size guarantees hold exactly on the surface.

* **Parametric faces** `s(u,v)` — meshed in the parameter rectangle. For a
  **developable** surface (plane, cylinder, cone) the parameter map is (up to a
  constant scale per axis) isometric, so plain 2-D meshing of the scaled
  parameter domain gives exact surface sizing — this is how the coax cylinders
  are meshed. For a general surface the parameter-space size target is divided by
  the local stretch of the first fundamental form (an isotropic metric
  approximation), verified on the sphere by area/where analytic.

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
    dim=Int(value[1]); tag=Int(value[2])
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
        length(spec)==nloops || throw(ArgumentError(
            "$caller: boundary_entities must have one entry per boundary loop"))
        entry=spec[loop_index]
        if entry===nothing || _entity_tuple(entry)
            entry
        elseif entry isa AbstractVector
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
    length(p) >= 3 || throw(ArgumentError("$caller: each point needs three coordinates"))
    q = try
        (Float64(p[1]), Float64(p[2]), Float64(p[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: point coordinates must be Float64-convertible: $(sprint(showerror, err))"))
    end
    (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) ||
        throw(ArgumentError("$caller: point has non-finite coordinates $q"))
    return q
end

# ── planar frame ────────────────────────────────────────────────────────────────
struct PlaneFrame
    o::NTuple{3,Float64}   # origin (a loop vertex)
    u::NTuple{3,Float64}   # in-plane unit axis 1
    v::NTuple{3,Float64}   # in-plane unit axis 2
    n::NTuple{3,Float64}   # unit normal
end

function _validate_coplanar_loops(loops, fr::PlaneFrame)
    scale = 1.0
    for loop in loops, p0 in loop
        p = _point3(p0, "mesh_planar_face")
        dx=p[1]-fr.o[1];dy=p[2]-fr.o[2];dz=p[3]-fr.o[3]
        (isfinite(dx)&&isfinite(dy)&&isfinite(dz)) ||
            throw(ArgumentError("mesh_planar_face: coordinate span overflowed Float64"))
        scale = max(scale, abs(dx),abs(dy),abs(dz))
    end
    tol = 256eps(Float64)*scale
    for (li, loop) in enumerate(loops), (pi, p0) in enumerate(loop)
        p = _point3(p0, "mesh_planar_face")
        dist = abs(_dot(_sub(p, fr.o), fr.n))
        isfinite(dist) || throw(ArgumentError("mesh_planar_face: plane-distance evaluation overflowed Float64"))
        dist <= tol ||
            throw(ArgumentError("mesh_planar_face: loop $li point $pi is not coplanar (plane distance $dist > tolerance $tol)"))
    end
    return nothing
end

"""Orthonormal in-plane frame of a coplanar loop (Newell's method for the normal)."""
function plane_frame(loop)
    m = length(loop)
    m >= 3 || throw(ArgumentError("plane_frame: loop needs ≥3 points"))
    nx=ny=nz=0.0
    @inbounds for i in 1:m
        p = _point3(loop[i], "plane_frame"); q = _point3(loop[mod1(i+1,m)], "plane_frame")
        nx += (p[2]-q[2])*(p[3]+q[3])
        ny += (p[3]-q[3])*(p[1]+q[1])
        nz += (p[1]-q[1])*(p[2]+q[2])
    end
    L = hypot(nx,ny,nz)
    (isfinite(L)&&L>0) || throw(ArgumentError("plane_frame: degenerate or numerically unrepresentable loop"))
    n = (nx/L, ny/L, nz/L)
    # in-plane axis: project the world axis least aligned with n
    ax = abs(n[1]) <= abs(n[2]) ?
         (abs(n[1]) <= abs(n[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(n[2]) <= abs(n[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    u = _unit(_cross(ax, n))
    v = _cross(n, u)          # already unit (n,u orthonormal)
    return PlaneFrame(_point3(loop[1], "plane_frame"), u, v, n)
end

@inline project(fr::PlaneFrame, p) = (_dot(_sub(p, fr.o), fr.u), _dot(_sub(p, fr.o), fr.v))
@inline lift(fr::PlaneFrame, a, b) =
    (fr.o[1]+a*fr.u[1]+b*fr.v[1], fr.o[2]+a*fr.u[2]+b*fr.v[2], fr.o[3]+a*fr.u[3]+b*fr.v[3])

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
    face_entity=_surface_entity_context(entity,"mesh_planar_face")
    fr = plane_frame(loops[1])
    _validate_coplanar_loops(loops, fr)
    xs = Float64[]; ys = Float64[]; segs = Tuple{Int,Int}[]
    for (loop_index,loop) in pairs(loops)
        _add_loop!(xs,ys,segs,loop,sf,fr,face_entity,boundary_entities,
                   length(loops),Int(loop_index))
    end
    T = constrained_delaunay(xs, ys, segs)
    sizefn = sf isa AbstractAnisoField ? nothing :
             (a, b) -> size_at(sf, lift(fr, a, b)...,face_entity)
    edgefn = sf isa AbstractAnisoField ? (ax,ay,bx,by) ->
             metric_edge_length(sf,lift(fr,ax,ay),lift(fr,bx,by);entity=face_entity) : nothing
    interior = refine!(T; min_angle_deg=min_angle_deg, max_area=max_area,
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
uniformly into `nθ = max(6, round(2πR/h))` sectors, and each quad is split into
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
    on(θ, z) = (c[1]+R*cos(θ)*ex[1]+R*sin(θ)*ey[1]+z*ez[1],
                c[2]+R*cos(θ)*ex[2]+R*sin(θ)*ey[2]+z*ez[2],
                c[3]+R*cos(θ)*ex[3]+R*sin(θ)*ey[3]+z*ez[3])
    bottom_entity=_surface_boundary_entity(boundary_entities,2,1,1,1,
        on(0.0,0.0),on(π,0.0),face_entity,"mesh_cylinder_face")
    top_entity=_surface_boundary_entity(boundary_entities,2,2,1,1,
        on(0.0,H),on(π,H),face_entity,"mesh_cylinder_face")
    # axial levels: size-graded along the axis (uses the real 3-D size field)
    apts, _ = mesh_segment(c, (c[1]+H*ez[1], c[2]+H*ez[2], c[3]+H*ez[3]), sf;
                           entity=face_entity)
    zlev = Float64[_dot(_sub(p, c), ez) for p in apts]     # 0 .. H, graded
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
            α=k/nsub
            push!(refined_z,k==nsub ? zlev[i+1] : zlev[i]*(1-α)+zlev[i+1]*α)
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
        a=node(out,out.tris[1,t]);b=node(out,out.tris[2,t]);c3=node(out,out.tris[3,t])
        A=triangle_area(a,b,c3);amin=_triangle_min_angle_deg(a,b,c3)
        A<=area*(1+256eps(Float64)) ||
            throw(ErrorException("mesh_cylinder_face: postcondition failed: triangle $t exceeds max_area"))
        amin+1024eps(Float64)>=angle ||
            throw(ErrorException("mesh_cylinder_face: postcondition failed: triangle $t violates min_angle_deg"))
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
                refined[out]=k==n ? z1 : z0+(z1-z0)*(k/n)
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

function _triangle_min_angle_deg(a,b,c)
    sides=sort([_norm(_sub(a,b)),_norm(_sub(b,c)),_norm(_sub(c,a))])
    sides[1]>0 || return 0.0
    q=sides./sides[3]
    return acosd(clamp((q[2]^2+1-q[1]^2)/(2q[2]),-1.0,1.0))
end

# ── general parametric face (isotropic-metric approximation) ─────────────────────
"""
    mesh_parametric_face(s, umin, umax, vmin, vmax, sf; nu=32, nv=32, kwargs...) -> Mesh

Mesh a parametric surface `s(u,v) -> (x,y,z)` over the rectangle
`[umin,umax]×[vmin,vmax]`. The parameter-space size target is the physical size
divided by the local stretch of the first fundamental form (isotropic metric
approximation via the mean of √E, √G), then the domain is 2-D meshed and mapped
to 3-D. Suitable for gently-curved patches; developable surfaces are exact via
[`mesh_cylinder_face`](@ref).
"""
function mesh_parametric_face(s, umin::Real, umax::Real, vmin::Real, vmax::Real,
                              sf::AbstractSizeField; min_angle_deg::Real=25.0,
                              max_area::Real=Inf, entity=nothing,
                              boundary_entities=nothing)
    face_entity=_surface_entity_context(entity,"mesh_parametric_face")
    umin=_surface_float(umin,"mesh_parametric_face","umin")
    umax=_surface_float(umax,"mesh_parametric_face","umax")
    vmin=_surface_float(vmin,"mesh_parametric_face","vmin")
    vmax=_surface_float(vmax,"mesh_parametric_face","vmax")
    (isfinite(umin) && isfinite(umax) && isfinite(vmin) && isfinite(vmax) &&
     umin < umax && vmin < vmax) ||
        throw(ArgumentError("mesh_parametric_face: require finite bounds umin < umax and vmin < vmax"))
    du = 1e-6*(umax-umin); dv = 1e-6*(vmax-vmin)
    (isfinite(du)&&isfinite(dv)&&du>0&&dv>0) ||
        throw(ArgumentError("mesh_parametric_face: parameter spans are below Float64 finite-difference resolution"))
    surf(u,v) = _point3(s(u,v), "mesh_parametric_face")
    # local stretch factors (√E, √G) from finite differences
    stretch(u,v) = begin
        ul = max(umin, u-du); ur = min(umax, u+du)
        vl = max(vmin, v-dv); vr = min(vmax, v+dv)
        su = _sub(surf(ur,v), surf(ul,v)); sv = _sub(surf(u,vr), surf(u,vl))
        lu=_norm(su);lv=_norm(sv);a=lu/(ur-ul);b=lv/(vr-vl)
        uu=(su[1]/lu,su[2]/lu,su[3]/lu);vv=(sv[1]/lv,sv[2]/lv,sv[3]/lv)
        sinang=_norm(_cross(uu,vv))
        (isfinite(a) && isfinite(b) && isfinite(sinang) && a>0 && b>0 &&
         sinang > 256eps(Float64)) ||
            throw(ArgumentError("mesh_parametric_face: singular or non-finite surface Jacobian at ($u,$v)"))
        (a,b)
    end
    # Isotropic fields use the established local-stretch approximation. Anisotropic
    # fields are checked directly on each physical edge below.
    sizefn = sf isa AbstractAnisoField ? nothing : (u,v) -> begin
        p = surf(u,v); h = size_at(sf,p...,face_entity); (a,b) = stretch(u,v)
        hp=h/max(a,b)
        (isfinite(hp)&&hp>0) || throw(ArgumentError(
            "mesh_parametric_face: parameter-space size is not representable at ($u,$v)"))
        hp
    end
    # PSLG = parameter rectangle boundary, graded in parameter space
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    corners=[(umin,vmin),(umax,vmin),(umax,vmax),(umin,vmax)]
    start=1
    for i in 1:4
        p=corners[i]; q=corners[mod1(i+1,4)]
        edge_entity=_surface_boundary_entity(boundary_entities,1,1,i,4,p,q,
            face_entity,"mesh_parametric_face")
        # metric length along the edge in physical space
        metric = _param_edge_metric(surf,sf,p,q,edge_entity)
        (isfinite(metric)&&metric>0&&metric<=prevfloat(Float64(typemax(Int)))) ||
            throw(ArgumentError("mesh_parametric_face: requested boundary node count exceeds the platform Int limit"))
        n = max(1, round(Int, metric))
        n<=typemax(Int32)-length(xs) ||
            throw(ArgumentError("mesh_parametric_face: boundary node count exceeds Int32 indexing"))
        for k in 0:n-1
            t=k/n; push!(xs, p[1]+t*(q[1]-p[1])); push!(ys, p[2]+t*(q[2]-p[2]))
        end
    end
    stop=length(xs)
    for i in start:stop; push!(segs, (i, i==stop ? start : i+1)); end
    T=constrained_delaunay(xs,ys,segs)
    edgefn=sf isa AbstractAnisoField ? (u1,v1,u2,v2) ->
        metric_edge_length(sf,surf(u1,v1),surf(u2,v2);entity=face_entity) : nothing
    interior=refine!(T; min_angle_deg=min_angle_deg, max_area=max_area,
                     size=sizefn,edge_metric=edgefn)
    m2=to_mesh(T; interior=interior)
    # map parameter mesh to 3-D
    nn=nnodes(m2); coords=Matrix{Float64}(undef,3,nn)
    @inbounds for i in 1:nn
        p=surf(m2.coords[1,i], m2.coords[2,i]); coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]
    end
    out=Mesh(coords; tris=copy(m2.tris))
    d=validate(out)
    d.ok || throw(ErrorException("mesh_parametric_face: produced an invalid surface mesh — "*join(d.messages,"; ")))
    return out
end

@inline function _param_edge_metric(s, sf, p, q, entity)
    N = 32; L = 0.0; prev = s(p[1], p[2])
    @inbounds for k in 1:N
        t = k/N; u = p[1]+t*(q[1]-p[1]); v = p[2]+t*(q[2]-p[2])
        cur = s(u,v)
        L += _metric_curve_increment(sf,prev,cur,entity); prev = cur
        isfinite(L) || throw(ArgumentError("mesh_parametric_face: boundary metric length overflowed Float64"))
    end
    return L
end

function _surface_float(x::Real,caller::AbstractString,name::AbstractString)
    try
        return Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be representable as Float64"))
    end
end

end # module MeshSurface
