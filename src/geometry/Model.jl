"""
    Model

A native geometry/entity kernel: tagged points, curves, loops, surfaces, and
volumes with physical groups. Meshing dispatches to Tessella's certified
simplex and transfinite kernels. This is not OpenCASCADE; unsupported CAD
statements remain explicit blockers.
"""
module Model

using ..MeshTypes: Mesh, validate, nnodes, ntris, ntets
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..SizeField: AbstractSizeField, ConstantSize, size_at
using ..Geometry: box_surface, cylinder_surface, sphere_surface, cone_surface
using ..Mesh3D: tetrahedralize, mesh_boolean, recover_segment3, recover_triangle3
using ..Mesh3D: mesh_covers_segment3, mesh_covers_triangle3

export GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
export add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
export embed!, translate_volume!, dilate_volume!, rotate_volume!
export add_physical_group!, set_physical_name!
export mesh_model_surface, mesh_model_volume, model_entity, model_physical_tags

mutable struct GeoModel
    points::Dict{Int,NTuple{3,Float64}}
    point_size::Dict{Int,Float64}
    curves::Dict{Int,NTuple{2,Int}}
    loops::Dict{Int,Vector{Int}}
    surfaces::Dict{Int,Vector{Int}}
    volumes::Dict{Int,Vector{Int}}
    physical::Dict{Tuple{Int,Int},Vector{Int}}
    physical_names::Dict{Tuple{Int,Int},String}
    box_extents::Dict{Int,NTuple{6,Float64}}
    cylinders::Dict{Int,NamedTuple{(:center,:axis,:radius,:height),
                                   Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64}}}
    spheres::Dict{Int,NamedTuple{(:center,:radius),Tuple{NTuple{3,Float64},Float64}}}
    cones::Dict{Int,NamedTuple{(:center,:axis,:r1,:r2,:height),
                               Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64,Float64}}}
    booleans::Dict{Int,NamedTuple{(:op,:a,:b),Tuple{Symbol,Int,Int}}}
    embeds::Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}
    next_tag::Vector{Int}
end

"""
    GeoModel()

Create an empty native geometry model. Entity tags are positive 32-bit integers;
passing `tag=0` to an `add_*!` function allocates the next tag in that entity
dimension. Geometry is added explicitly and can then be meshed with
[`mesh_model_surface`](@ref) or [`mesh_model_volume`](@ref).
"""
GeoModel() = GeoModel(Dict{Int,NTuple{3,Float64}}(), Dict{Int,Float64}(),
                      Dict{Int,NTuple{2,Int}}(), Dict{Int,Vector{Int}}(),
                      Dict{Int,Vector{Int}}(), Dict{Int,Vector{Int}}(),
                      Dict{Tuple{Int,Int},Vector{Int}}(),
                      Dict{Tuple{Int,Int},String}(),
                      Dict{Int,NTuple{6,Float64}}(),
                      Dict{Int,NamedTuple{(:center,:axis,:radius,:height),
                           Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64}}}(),
                      Dict{Int,NamedTuple{(:center,:radius),Tuple{NTuple{3,Float64},Float64}}}(),
                      Dict{Int,NamedTuple{(:center,:axis,:r1,:r2,:height),
                           Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64,Float64}}}(),
                      Dict{Int,NamedTuple{(:op,:a,:b),Tuple{Symbol,Int,Int}}}(),
                      Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}(),
                      Int[0,0,0,0])

function _tag(value, caller, dim::Int)
    value isa Integer || throw(ArgumentError("$caller: tag must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: tag must not be Bool"))
    value<0 && throw(ArgumentError("$caller: tag must be non-negative"))
    value>typemax(Int32) && throw(ArgumentError("$caller: tag exceeds Int32"))
    return Int(value)
end

function _signed_curve_tag(value, caller)
    value isa Integer || throw(ArgumentError("$caller: curve tag must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: curve tag must not be Bool"))
    value==0 && throw(ArgumentError("$caller: curve tag must be nonzero"))
    (-typemax(Int32)<=value<=typemax(Int32)) || throw(ArgumentError(
        "$caller: curve tag magnitude exceeds Int32"))
    return Int(value)
end

function _dimension(value, caller)
    value isa Integer || throw(ArgumentError("$caller: dimension must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: dimension must not be Bool"))
    (0<=value<=3) || throw(ArgumentError("$caller: dimension must be in 0:3"))
    return Int(value)
end

function _alloc_tag!(m::GeoModel, dim::Int, requested::Int, caller)
    if requested==0
        m.next_tag[dim+1]<typemax(Int32) || throw(ArgumentError(
            "$caller: no automatic tags remain in dimension $dim"))
        m.next_tag[dim+1]+=1
        return m.next_tag[dim+1]
    end
    m.next_tag[dim+1]=max(m.next_tag[dim+1], requested)
    return requested
end

function _alloc_physical_tag(m::GeoModel, dim::Int, requested::Int, caller)
    requested!=0 && return requested
    current=0
    for (d,t) in keys(m.physical)
        d==dim && (current=max(current,t))
    end
    current<typemax(Int32) || throw(ArgumentError(
        "$caller: no automatic physical tags remain in dimension $dim"))
    return current+1
end

function _finite3(x,y,z,caller)
    p=try (Float64(x),Float64(y),Float64(z)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: coordinates must be Float64-representable"))
    end
    all(isfinite,p) || throw(ArgumentError("$caller: coordinates must be finite"))
    return p
end

function _finite_vector3(value, caller, what)
    values=try
        Tuple(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $what must be an iterable with three components"))
    end
    length(values)==3 || throw(ArgumentError(
        "$caller: $what must have exactly three components"))
    return _finite3(values[1],values[2],values[3],caller)
end

function _finite_scalar(value, caller, what)
    result=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $what must be Float64-representable"))
    end
    isfinite(result) || throw(ArgumentError("$caller: $what must be finite"))
    return result
end

function _finite_result(values, caller)
    all(isfinite,values) || throw(ArgumentError(
        "$caller: transformed geometry must remain finite"))
    return values
end

"""
    add_point!(model, x, y, z; tag=0, mesh_size=1.0) -> tag

Add a point with finite coordinates and a positive characteristic mesh size.
`tag=0` requests automatic tag allocation.
"""
function add_point!(m::GeoModel, x, y, z; tag::Integer=0, mesh_size::Real=1.0)
    caller="add_point!"
    p=_finite3(x,y,z,caller)
    h=try Float64(mesh_size) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: mesh_size must be Float64-representable"))
    end
    (isfinite(h) && h>0) || throw(ArgumentError("$caller: mesh_size must be positive"))
    t=_alloc_tag!(m,0,_tag(tag,caller,0),caller)
    haskey(m.points,t) && throw(ArgumentError("$caller: Point[$t] already exists"))
    m.points[t]=p; m.point_size[t]=h
    return t
end

"""
    add_line!(model, start, stop; tag=0) -> tag

Add a straight curve between two distinct existing point tags.
"""
function add_line!(m::GeoModel, a, b; tag::Integer=0)
    caller="add_line!"
    ta=_tag(a,caller,1); tb=_tag(b,caller,1)
    ta==tb && throw(ArgumentError("$caller: line endpoints must be distinct"))
    haskey(m.points,ta) || throw(ArgumentError("$caller: unknown Point[$ta]"))
    haskey(m.points,tb) || throw(ArgumentError("$caller: unknown Point[$tb]"))
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller)
    haskey(m.curves,t) && throw(ArgumentError("$caller: Curve[$t] already exists"))
    m.curves[t]=(ta,tb)
    return t
end

"""
    add_curve_loop!(model, curves; tag=0) -> tag

Add an ordered, closed loop of at least three existing curves. A negative curve
tag traverses that curve in reverse orientation. Consecutive oriented curves
must share endpoints.
"""
function add_curve_loop!(m::GeoModel, curves; tag::Integer=0)
    caller="add_curve_loop!"
    ids=Int[_signed_curve_tag(c,caller) for c in curves]
    length(ids)>=3 || throw(ArgumentError("$caller: a loop needs at least three curves"))
    for id in ids
        haskey(m.curves,abs(id)) || throw(ArgumentError("$caller: unknown Curve[$(abs(id))]"))
    end
    oriented=map(ids) do id
        a,b=m.curves[abs(id)]
        id>0 ? (a,b) : (b,a)
    end
    for i in eachindex(oriented)
        current=oriented[i]
        following=oriented[mod1(i+1,length(oriented))]
        current[2]==following[1] || throw(ArgumentError(
            "$caller: oriented Curve[$(ids[i])] ends at Point[$(current[2])], " *
            "but Curve[$(ids[mod1(i+1,length(ids))])] starts at Point[$(following[1])]"))
    end
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller)
    haskey(m.loops,t) && throw(ArgumentError("$caller: Loop[$t] already exists"))
    m.loops[t]=ids
    return t
end

"""
    add_plane_surface!(model, loops; tag=0) -> tag

Add a planar surface bounded by existing curve loops. The first loop is the
outer boundary and subsequent loops are holes.
"""
function add_plane_surface!(m::GeoModel, loops; tag::Integer=0)
    caller="add_plane_surface!"
    ids=Int[_tag(ℓ,caller,2) for ℓ in loops]
    isempty(ids) && throw(ArgumentError("$caller: need an outer loop"))
    for id in ids
        haskey(m.loops,id) || throw(ArgumentError("$caller: unknown Loop[$id]"))
    end
    t=_alloc_tag!(m,2,_tag(tag,caller,2),caller)
    haskey(m.surfaces,t) && throw(ArgumentError("$caller: Surface[$t] already exists"))
    m.surfaces[t]=ids
    return t
end

"""
    add_box!(model, x, y, z, dx, dy, dz; tag=0) -> tag

Add an axis-aligned box volume with finite origin `(x,y,z)` and positive
extents `(dx,dy,dz)`.
"""
function add_box!(m::GeoModel, xmin, ymin, zmin, dx, dy, dz; tag::Integer=0)
    caller="add_box!"
    origin=_finite3(xmin,ymin,zmin,caller)
    d=_finite3(dx,dy,dz,caller)
    (d[1]>0 && d[2]>0 && d[3]>0) || throw(ArgumentError("$caller: extents must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.box_extents[t]=(origin[1],origin[2],origin[3],d[1],d[2],d[3])
    return t
end

"""
    add_cylinder!(model, x, y, z, dx, dy, dz, radius; tag=0) -> tag

Add a cylinder whose base center is `(x,y,z)` and whose finite, nonzero axis
vector is `(dx,dy,dz)`. `radius` must be finite and positive.
"""
function add_cylinder!(m::GeoModel, x, y, z, dx, dy, dz, radius; tag::Integer=0)
    caller="add_cylinder!"
    c=_finite3(x,y,z,caller); a=_finite3(dx,dy,dz,caller)
    h=hypot(a...)
    (isfinite(h) && h>0) || throw(ArgumentError("$caller: axis must have finite positive length"))
    r=_finite_scalar(radius,caller,"radius")
    r>0 || throw(ArgumentError("$caller: radius must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.cylinders[t]=(center=c, axis=a, radius=r, height=h)
    return t
end

"""
    add_sphere!(model, x, y, z, radius; tag=0) -> tag

Add a sphere with a finite center and finite positive radius.
"""
function add_sphere!(m::GeoModel, x, y, z, radius; tag::Integer=0)
    caller="add_sphere!"
    c=_finite3(x,y,z,caller)
    r=_finite_scalar(radius,caller,"radius")
    r>0 || throw(ArgumentError("$caller: radius must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.spheres[t]=(center=c, radius=r)
    return t
end

"""
    add_cone!(model, x, y, z, dx, dy, dz, r1, r2; tag=0) -> tag

Add a cone or conical frustum along the finite, nonzero axis `(dx,dy,dz)`.
Both radii must be finite and non-negative, and at least one must be positive.
"""
function add_cone!(m::GeoModel, x, y, z, dx, dy, dz, r1, r2; tag::Integer=0)
    caller="add_cone!"
    c=_finite3(x,y,z,caller); a=_finite3(dx,dy,dz,caller)
    h=hypot(a...)
    (isfinite(h) && h>0) || throw(ArgumentError("$caller: axis must have finite positive length"))
    ra=_finite_scalar(r1,caller,"r1")
    rb=_finite_scalar(r2,caller,"r2")
    (ra>=0 && rb>=0 && (ra>0 || rb>0)) || throw(ArgumentError(
        "$caller: radii must be non-negative with at least one positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.cones[t]=(center=c, axis=a, r1=ra, r2=rb, height=h)
    return t
end

"""
    embed!(model, dim, tags, target_dim, target_tag) -> target_tag

Embed existing points or curves in a surface, or existing points, curves, or
surfaces in a volume. The operation is atomic: invalid or duplicate entries do
not add a partial embedding.
"""
function embed!(m::GeoModel, dim, tags, target_dim, target_tag)
    caller="embed!"
    d=_dimension(dim,caller); td=_dimension(target_dim,caller)
    tt=_tag(target_tag,caller,td)
    (d==0 && td==2) || (d==1 && td==2) || (d==0 && td==3) ||
    (d==1 && td==3) || (d==2 && td==3) || throw(ArgumentError(
        "$caller: supported embeddings are Point/Line In Surface and Point/Line/Surface In Volume (got dim=$d in dim=$td)"))
    if td==2
        haskey(m.surfaces,tt) || throw(ArgumentError("$caller: unknown Surface[$tt]"))
    else
        haskey(m.volumes,tt) || throw(ArgumentError("$caller: unknown Volume[$tt]"))
    end
    ents=Int[_tag(raw,caller,d) for raw in tags]
    length(unique(ents))==length(ents) || throw(ArgumentError(
        "$caller: embedding contains duplicate entity tags"))
    existing=get(m.embeds,(td,tt),NTuple{2,Int}[])
    for ent in ents
        if d==0
            haskey(m.points,ent) || throw(ArgumentError("$caller: unknown Point[$ent]"))
        elseif d==1
            haskey(m.curves,ent) || throw(ArgumentError("$caller: unknown Curve[$ent]"))
        else
            haskey(m.surfaces,ent) || throw(ArgumentError("$caller: unknown Surface[$ent]"))
        end
        (d,ent) in existing && throw(ArgumentError(
            "$caller: entity ($d,$ent) already embedded in ($td,$tt)"))
    end
    isempty(ents) || append!(get!(Vector{NTuple{2,Int}},m.embeds,(td,tt)),
                             ((d,ent) for ent in ents))
    return tt
end

function _quarter_turns(angle, caller)
    turns=angle/(π/2)
    (isfinite(turns) && abs(turns)<=typemax(Int)) || throw(ArgumentError(
        "$caller: angle is outside the supported integer-turn range"))
    k=try
        round(Int,turns)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: angle is outside the supported integer-turn range"))
    end
    abs(turns-k)<=1e-9 || throw(ArgumentError(
        "$caller: angle must be an integer multiple of π/2 (got $angle)"))
    return mod(k,4)
end

function _axis_kind(axis, caller)
    L=hypot(axis...)
    L>0 || throw(ArgumentError("$caller: axis must be nonzero"))
    u=(axis[1]/L,axis[2]/L,axis[3]/L)
    if abs(abs(u[1])-1)<=1e-12 && abs(u[2])<=1e-12 && abs(u[3])<=1e-12
        return (:x, u[1]>0 ? 1 : -1)
    elseif abs(abs(u[2])-1)<=1e-12 && abs(u[1])<=1e-12 && abs(u[3])<=1e-12
        return (:y, u[2]>0 ? 1 : -1)
    elseif abs(abs(u[3])-1)<=1e-12 && abs(u[1])<=1e-12 && abs(u[2])<=1e-12
        return (:z, u[3]>0 ? 1 : -1)
    end
    throw(ArgumentError("$caller: only coordinate-axis rotations are implemented"))
end

function _rot90(p, origin, kind, sign, k)
    x,y,z=p[1]-origin[1], p[2]-origin[2], p[3]-origin[3]
    kk=mod(sign*k,4)
    if kind===:z
        q=kk==0 ? (x,y,z) : kk==1 ? (-y,x,z) : kk==2 ? (-x,-y,z) : (y,-x,z)
    elseif kind===:x
        q=kk==0 ? (x,y,z) : kk==1 ? (x,-z,y) : kk==2 ? (x,-y,-z) : (x,z,-y)
    else
        q=kk==0 ? (x,y,z) : kk==1 ? (z,y,-x) : kk==2 ? (-x,y,-z) : (-z,y,x)
    end
    return (q[1]+origin[1], q[2]+origin[2], q[3]+origin[3])
end

function _dilate_point(p, center, s)
    return (center[1]+s*(p[1]-center[1]),
            center[2]+s*(p[2]-center[2]),
            center[3]+s*(p[3]-center[3]))
end

"""
    translate_volume!(model, tag, offset) -> tag

Translate a native primitive volume by the finite three-component `offset`.
The model is unchanged if a translated coordinate is not representable as a
finite `Float64` or the volume has no translatable native encoding.
"""
function translate_volume!(m::GeoModel, tag, offset)
    caller="translate_volume!"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    delta=_finite_vector3(offset,caller,"offset")
    if haskey(m.box_extents,t)
        x0,y0,z0,dx,dy,dz=m.box_extents[t]
        origin=_finite_result((x0+delta[1],y0+delta[2],z0+delta[3]),caller)
        m.box_extents[t]=(origin[1],origin[2],origin[3],dx,dy,dz)
    elseif haskey(m.cylinders,t)
        cyl=m.cylinders[t]
        center=_finite_result((cyl.center[1]+delta[1],cyl.center[2]+delta[2],
                               cyl.center[3]+delta[3]),caller)
        m.cylinders[t]=(center=center,axis=cyl.axis,radius=cyl.radius,height=cyl.height)
    elseif haskey(m.spheres,t)
        sph=m.spheres[t]
        center=_finite_result((sph.center[1]+delta[1],sph.center[2]+delta[2],
                               sph.center[3]+delta[3]),caller)
        m.spheres[t]=(center=center,radius=sph.radius)
    elseif haskey(m.cones,t)
        cone=m.cones[t]
        center=_finite_result((cone.center[1]+delta[1],cone.center[2]+delta[2],
                               cone.center[3]+delta[3]),caller)
        m.cones[t]=(center=center,axis=cone.axis,r1=cone.r1,r2=cone.r2,
                    height=cone.height)
    else
        throw(ArgumentError("$caller: Volume[$t] has no translatable native encoding"))
    end
    return t
end

"""
    dilate_volume!(model, tag, center, scale) -> tag

Dilate a native primitive volume about the three-component finite `center`.
`scale` must be finite and positive.
"""
function dilate_volume!(m::GeoModel, tag, center, scale)
    caller="dilate_volume!"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    c=_finite_vector3(center,caller,"center")
    s=_finite_scalar(scale,caller,"scale")
    s>0 || throw(ArgumentError("$caller: scale must be positive"))
    if haskey(m.box_extents,t)
        x0,y0,z0,dx,dy,dz=m.box_extents[t]
        p0=_dilate_point((x0,y0,z0),c,s)
        transformed=_finite_result((p0[1],p0[2],p0[3],dx*s,dy*s,dz*s),caller)
        m.box_extents[t]=transformed
    elseif haskey(m.cylinders,t)
        cyl=m.cylinders[t]
        center=_dilate_point(cyl.center,c,s)
        radius,height=cyl.radius*s,cyl.height*s
        _finite_result((center...,radius,height),caller)
        m.cylinders[t]=(center=center, axis=cyl.axis, radius=radius, height=height)
    elseif haskey(m.spheres,t)
        sph=m.spheres[t]
        center=_dilate_point(sph.center,c,s); radius=sph.radius*s
        _finite_result((center...,radius),caller)
        m.spheres[t]=(center=center, radius=radius)
    elseif haskey(m.cones,t)
        cone=m.cones[t]
        center=_dilate_point(cone.center,c,s)
        r1,r2,height=cone.r1*s,cone.r2*s,cone.height*s
        _finite_result((center...,r1,r2,height),caller)
        m.cones[t]=(center=center, axis=cone.axis, r1=r1, r2=r2, height=height)
    else
        throw(ArgumentError("$caller: Volume[$t] has no dilatable native encoding"))
    end
    return t
end

"""
    rotate_volume!(model, tag, axis, origin, angle) -> tag

Rotate a native primitive volume about a coordinate-aligned `axis` through
`origin`. The finite angle must be an integer multiple of `π/2`.
"""
function rotate_volume!(m::GeoModel, tag, axis, origin, angle)
    caller="rotate_volume!"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    ax=_finite_vector3(axis,caller,"axis")
    o=_finite_vector3(origin,caller,"origin")
    θ=_finite_scalar(angle,caller,"angle")
    k=_quarter_turns(θ,caller)
    kind,sgn=_axis_kind(ax,caller)
    rot(p)=_rot90(p,o,kind,sgn,k)
    if haskey(m.box_extents,t)
        x0,y0,z0,dx,dy,dz=m.box_extents[t]
        corners=[rot((x0+ix*dx,y0+iy*dy,z0+iz*dz)) for ix in (0,1), iy in (0,1), iz in (0,1)]
        _finite_result(Iterators.flatten(corners),caller)
        xs=sort!(unique([p[1] for p in corners]))
        ys=sort!(unique([p[2] for p in corners]))
        zs=sort!(unique([p[3] for p in corners]))
        (length(xs)==2 && length(ys)==2 && length(zs)==2) || throw(ArgumentError(
            "$caller: rotated box is no longer axis-aligned"))
        m.box_extents[t]=(xs[1],ys[1],zs[1],xs[2]-xs[1],ys[2]-ys[1],zs[2]-zs[1])
    elseif haskey(m.cylinders,t)
        cyl=m.cylinders[t]
        endp=(cyl.center[1]+cyl.axis[1],cyl.center[2]+cyl.axis[2],cyl.center[3]+cyl.axis[3])
        c2=rot(cyl.center); e2=rot(endp)
        axis=(e2[1]-c2[1],e2[2]-c2[2],e2[3]-c2[3])
        _finite_result((c2...,axis...),caller)
        m.cylinders[t]=(center=c2, axis=axis,
                        radius=cyl.radius, height=cyl.height)
    elseif haskey(m.spheres,t)
        sph=m.spheres[t]
        center=rot(sph.center)
        _finite_result(center,caller)
        m.spheres[t]=(center=center, radius=sph.radius)
    elseif haskey(m.cones,t)
        cone=m.cones[t]
        endp=(cone.center[1]+cone.axis[1],cone.center[2]+cone.axis[2],cone.center[3]+cone.axis[3])
        c2=rot(cone.center); e2=rot(endp)
        axis=(e2[1]-c2[1],e2[2]-c2[2],e2[3]-c2[3])
        _finite_result((c2...,axis...),caller)
        m.cones[t]=(center=c2, axis=axis,
                    r1=cone.r1, r2=cone.r2, height=cone.height)
    else
        throw(ArgumentError("$caller: Volume[$t] has no rotatable native encoding"))
    end
    return t
end

"""
    boolean_volumes!(model, op, a, b; tag=0) -> tag

Add a native Boolean volume combining existing volumes `a` and `b`. Supported
operations are `:union`, `:intersection`, and `:difference` (`a \\ b`).
"""
function boolean_volumes!(m::GeoModel, op::Symbol, a, b; tag::Integer=0)
    caller="boolean_volumes!"
    op in (:union,:intersection,:difference) || throw(ArgumentError(
        "$caller: op must be :union, :intersection, or :difference"))
    ta=_tag(a,caller,3); tb=_tag(b,caller,3)
    haskey(m.volumes,ta) || throw(ArgumentError("$caller: unknown Volume[$ta]"))
    haskey(m.volumes,tb) || throw(ArgumentError("$caller: unknown Volume[$tb]"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.booleans[t]=(op=op, a=ta, b=tb)
    return t
end

function _has_entity(m::GeoModel, dim::Int, tag::Int)
    dim==0 && return haskey(m.points,tag)
    dim==1 && return haskey(m.curves,tag)
    dim==2 && return haskey(m.surfaces,tag)
    return haskey(m.volumes,tag)
end

"""
    add_physical_group!(model, dim, tags; tag=0, name="") -> tag

Group one or more existing entities of dimension `dim` under an independent
physical tag. Automatic physical tags use a namespace separate from entity
tags. An optional nonempty `name` is recorded with the group.
"""
function add_physical_group!(m::GeoModel, dim::Integer, tags; tag::Integer=0, name::AbstractString="")
    caller="add_physical_group!"
    d=_dimension(dim,caller)
    ents=Int[_tag(t,caller,d) for t in tags]
    isempty(ents) && throw(ArgumentError("$caller: physical group needs at least one entity"))
    length(unique(ents))==length(ents) || throw(ArgumentError(
        "$caller: physical group contains duplicate entity tags"))
    for ent in ents
        _has_entity(m,d,ent) || throw(ArgumentError(
            "$caller: unknown entity ($d,$ent)"))
    end
    group_name=String(name)
    pt=_alloc_physical_tag(m,d,_tag(tag,caller,d),caller)
    haskey(m.physical,(d,pt)) && throw(ArgumentError("$caller: Physical($d,$pt) already exists"))
    m.physical[(d,pt)]=ents
    isempty(group_name) || (m.physical_names[(d,pt)]=group_name)
    return pt
end

"""
    set_physical_name!(model, dim, tag, name) -> name

Replace the name of an existing physical group.
"""
function set_physical_name!(m::GeoModel, dim::Integer, tag::Integer, name::AbstractString)
    caller="set_physical_name!"
    d=_dimension(dim,caller); t=_tag(tag,caller,d)
    key=(d,t)
    haskey(m.physical,key) || throw(ArgumentError("set_physical_name!: unknown Physical$key"))
    m.physical_names[key]=String(name)
    return name
end

"""
    model_entity(model, dim, tag)

Return the stored native entity representation for `(dim,tag)`, or `nothing`
when that entity does not exist.
"""
function model_entity(m::GeoModel, dim::Integer, tag::Integer)
    caller="model_entity"
    d=_dimension(dim,caller); t=_tag(tag,caller,d)
    d==0 && return get(m.points,t,nothing)
    d==1 && return get(m.curves,t,nothing)
    d==2 && return get(m.surfaces,t,nothing)
    d==3 && return get(m.volumes,t,nothing)
end

"""
    model_physical_tags(model, dim, tag) -> Vector{Int}

Return a copy of the entity tags in a physical group, or an empty vector when
the group does not exist.
"""
function model_physical_tags(m::GeoModel, dim::Integer, tag::Integer)
    caller="model_physical_tags"
    d=_dimension(dim,caller); t=_tag(tag,caller,d)
    return copy(get(m.physical,(d,t),Int[]))
end

function _loop_points(m::GeoModel, loop_id::Int)
    curves=m.loops[loop_id]
    pts=Int[]
    for signed in curves
        a,b=m.curves[abs(signed)]
        signed>0 ? push!(pts,a) : push!(pts,b)
    end
    return pts
end

function _add_surface_point!(xs, ys, index, m::GeoModel, pid::Int, caller)
    haskey(index, pid) && return index[pid]
    haskey(m.points,pid) || throw(ArgumentError("$caller: unknown Point[$pid]"))
    p=m.points[pid]
    abs(p[3])<=1e-12 || throw(ArgumentError(
        "$caller: Point[$pid] is not planar in z=0 (got z=$(p[3]))"))
    push!(xs,p[1]); push!(ys,p[2])
    index[pid]=length(xs)
    return index[pid]
end

function _node_at(mesh::Mesh, p; atol=1e-12)
    @inbounds for i in 1:nnodes(mesh)
        hypot(mesh.coords[1,i]-p[1],mesh.coords[2,i]-p[2],mesh.coords[3,i]-p[3])<=atol && return i
    end
    return 0
end

function _on_segment(x, p, q; atol=1e-12)
    vx,vy=q[1]-p[1], q[2]-p[2]
    wx,wy=x[1]-p[1], x[2]-p[2]
    L2=vx*vx+vy*vy
    L2>0 || return hypot(wx,wy)<=atol
    cross=vx*wy-vy*wx
    abs(cross)<=atol*sqrt(L2) || return false
    t=(wx*vx+wy*vy)/L2
    return -atol<=t<=1+atol
end

function _mesh_covers_segment(mesh::Mesh, p, q; atol=1e-12)
    a=_node_at(mesh,p; atol=atol); b=_node_at(mesh,q; atol=atol)
    a==0 && return false; b==0 && return false
    a==b && return true
    adj=Dict{Int,Vector{Int}}()
    @inbounds for t in 1:ntris(mesh), e in ((1,2),(2,3),(3,1))
        i=Int(mesh.tris[e[1],t]); j=Int(mesh.tris[e[2],t])
        pi=(mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])
        pj=(mesh.coords[1,j],mesh.coords[2,j],mesh.coords[3,j])
        (_on_segment(pi,p,q; atol=atol) && _on_segment(pj,p,q; atol=atol)) || continue
        push!(get!(Vector{Int}, adj, i), j)
        push!(get!(Vector{Int}, adj, j), i)
    end
    seen=falses(nnodes(mesh)); qn=Int[a]; seen[a]=true; head=1
    while head<=length(qn)
        v=qn[head]; head+=1
        v==b && return true
        for u in get(adj, v, Int[])
            seen[u] && continue
            seen[u]=true; push!(qn,u)
        end
    end
    return false
end

"""
    mesh_model_surface(model, tag; min_angle_deg=25.0) -> Mesh

Mesh a native planar surface in `z=0`, including holes and embedded points or
curves. Point characteristic lengths drive refinement. The returned triangle
mesh is validated before it is returned.
"""
function mesh_model_surface(m::GeoModel, tag::Integer; min_angle_deg::Real=25.0)
    caller="mesh_model_surface"
    t=_tag(tag,caller,2)
    haskey(m.surfaces,t) || throw(ArgumentError("$caller: unknown Surface[$t]"))
    loops=m.surfaces[t]
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    index=Dict{Int,Int}()
    for loop_id in loops
        ids=_loop_points(m,loop_id)
        loop_idx=Int[_add_surface_point!(xs,ys,index,m,pid,caller) for pid in ids]
        nloop=length(loop_idx)
        nloop>=3 || throw(ArgumentError("$caller: Loop[$loop_id] needs at least three points"))
        for k in 1:nloop
            push!(segs,(loop_idx[k], loop_idx[mod1(k+1,nloop)]))
        end
    end
    embedded=get(m.embeds,(2,t),NTuple{2,Int}[])
    internal=Tuple{Int,Int}[]
    for (edim,etag) in embedded
        if edim==0
            _add_surface_point!(xs,ys,index,m,etag,caller)
        elseif edim==1
            haskey(m.curves,etag) || throw(ArgumentError("$caller: unknown embedded Curve[$etag]"))
            a,b=m.curves[etag]
            ia=_add_surface_point!(xs,ys,index,m,a,caller)
            ib=_add_surface_point!(xs,ys,index,m,b,caller)
            ia==ib && throw(ArgumentError("$caller: embedded Curve[$etag] has coincident endpoints"))
            push!(internal,(ia,ib))
        else
            throw(ArgumentError("$caller: unsupported embedding dimension $edim"))
        end
    end
    T=constrained_delaunay(xs,ys,segs; internal_segments=internal)
    hmin=minimum(m.point_size[pid] for pid in keys(index))
    sizefn=(x,y)->hmin
    interior=refine!(T; min_angle_deg=min_angle_deg, size=sizefn)
    mesh=to_mesh(T; interior=interior)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntris(mesh)>0 || throw(ErrorException("$caller: Surface[$t] produced no triangles"))
    for (edim,etag) in embedded
        if edim==0
            p=m.points[etag]
            _node_at(mesh,p)==0 && throw(ErrorException(
                "$caller: embedded Point[$etag] at $p is not a mesh node"))
        elseif edim==1
            a,b=m.curves[etag]
            _mesh_covers_segment(mesh, m.points[a], m.points[b]) || throw(ErrorException(
                "$caller: embedded Curve[$etag] is not a chain of mesh edges"))
        end
    end
    return mesh
end

function _volume_surface(m::GeoModel, t::Int)
    if haskey(m.box_extents,t)
        x0,y0,z0,dx,dy,dz=m.box_extents[t]
        return box_surface(x0,x0+dx,y0,y0+dy,z0,z0+dz)
    elseif haskey(m.cylinders,t)
        c=m.cylinders[t]
        return cylinder_surface(c.center,c.axis,c.radius,c.height)
    elseif haskey(m.spheres,t)
        s=m.spheres[t]
        return sphere_surface(s.center,s.radius)
    elseif haskey(m.cones,t)
        c=m.cones[t]
        return cone_surface(c.center,c.axis,c.r1,c.r2,c.height)
    elseif haskey(m.booleans,t)
        spec=m.booleans[t]
        A=_volume_surface(m,spec.a); B=_volume_surface(m,spec.b)
        return mesh_boolean(A,B,spec.op)
    end
    throw(ArgumentError("mesh_model_volume: Volume[$t] has no native solid encoding"))
end

"""
    mesh_model_volume(model, tag) -> Mesh

Mesh a native primitive or Boolean volume, recovering supported embedded
points, curves, and unholed planar sheets. The returned tetrahedral mesh is
validated before it is returned; unsupported solid encodings raise an explicit
error.
"""
function mesh_model_volume(m::GeoModel, tag::Integer)
    caller="mesh_model_volume"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    surface=_volume_surface(m,t)
    extra=NTuple{3,Float64}[]
    lines=NTuple{2,NTuple{3,Float64}}[]
    sheets=NTuple{3,NTuple{3,Float64}}[]
    for (edim,etag) in get(m.embeds,(3,t),NTuple{2,Int}[])
        if edim==0
            haskey(m.points,etag) || throw(ArgumentError("$caller: unknown embedded Point[$etag]"))
            push!(extra, m.points[etag])
        elseif edim==1
            haskey(m.curves,etag) || throw(ArgumentError("$caller: unknown embedded Curve[$etag]"))
            a,b=m.curves[etag]
            push!(lines, (m.points[a], m.points[b]))
            push!(extra, m.points[a]); push!(extra, m.points[b])
        elseif edim==2
            haskey(m.surfaces,etag) || throw(ArgumentError("$caller: unknown embedded Surface[$etag]"))
            loops=m.surfaces[etag]
            isempty(loops) && throw(ArgumentError("$caller: embedded Surface[$etag] has no loop"))
            ids=_loop_points(m, loops[1])
            length(ids)>=3 || throw(ArgumentError(
                "$caller: embedded Surface[$etag] needs at least three points"))
            length(loops)==1 || throw(ArgumentError(
                "$caller: holed sheets In Volume are a blocker"))
            pts=[m.points[pid] for pid in ids]
            for i in 2:length(pts)-1
                push!(sheets, (pts[1], pts[i], pts[i+1]))
            end
        else
            throw(ArgumentError("$caller: unsupported embedding dimension $edim"))
        end
    end
    mesh=isempty(extra) ? tetrahedralize(surface) : tetrahedralize(surface; interior_points=extra)
    for (p,q) in lines
        mesh=recover_segment3(mesh,p,q)
        mesh_covers_segment3(mesh,p,q) || throw(ErrorException(
            "$caller: embedded curve is not a chain of tetrahedron edges"))
    end
    for (a,b,c) in sheets
        mesh=recover_triangle3(mesh,a,b,c)
        mesh_covers_triangle3(mesh,a,b,c) || throw(ErrorException(
            "$caller: embedded sheet is not a union of tetrahedron faces"))
    end
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntets(mesh)>0 || throw(ErrorException("$caller: Volume[$t] produced no tetrahedra"))
    return mesh
end

end # module
