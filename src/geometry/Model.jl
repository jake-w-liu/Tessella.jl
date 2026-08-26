"""
    Model

A native geometry/entity kernel: tagged points, curves, loops, surfaces, and
volumes with physical groups and classified surface/volume mixed-mesh projection.
Meshing dispatches to Tessella's certified simplex and transfinite kernels. This
is not OpenCASCADE; unsupported CAD statements remain explicit blockers.
"""
module Model

using ..MeshTypes: Mesh, validate, nnodes, nsegs, ntris, ntets
using ..Elements: ElementBlock, MixedEntity, MixedEntityData,
                  MixedPeriodicLink, MixedMesh
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..SizeField: AbstractSizeField, ConstantSize, size_at
using ..Geometry: box_surface, cylinder_surface, sphere_surface, cone_surface
using ..Mesh3D: tetrahedralize, mesh_boolean, recover_segment3, recover_triangle3
using ..Mesh3D: mesh_covers_segment3, mesh_covers_triangle3,
                _tet_edge_set, _mesh_covering_faces3, _certify_surface_fill
using ..Periodic: periodic_identify_affine
using ..Transform: _affine_coordinate, _transform_homogeneous

export GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
export add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
export embed!, translate_volume!, dilate_volume!, rotate_volume!
export ModelPeriodicConstraint, set_periodic!, model_periodic_constraints,
       model_periodic_nodes, model_to_mixed
export add_physical_group!, set_physical_name!
export mesh_model_surface, mesh_model_volume, model_entity, model_physical_tags

"""
Owned affine relation between two straight native curves. `affine` maps the
master curve to the slave curve in Gmsh row-major 4×4 order. `reversed` records
whether the master start maps to the slave end.
"""
struct ModelPeriodicConstraint
    dim::Int
    slave_entity::Int32
    master_entity::Int32
    affine::NTuple{16,Float64}
    reversed::Bool
    atol::Float64
end

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
    periodic::Dict{Tuple{Int,Int},ModelPeriodicConstraint}
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
                      Dict{Tuple{Int,Int},ModelPeriodicConstraint}(),
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

function _periodic_entity_tags(values,caller::AbstractString,name::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: $name entities must be a vector or tuple"))
    return Int[_tag(value,caller,1) for value in values]
end

function _model_periodic_tolerance(value,caller::AbstractString)
    value isa Bool && throw(ArgumentError("$caller: atol must not be Bool"))
    value isa Real || throw(ArgumentError("$caller: atol must be real"))
    tolerance=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: atol must be Float64-representable"))
    end
    (isfinite(tolerance) && tolerance>=0) || throw(ArgumentError(
        "$caller: atol must be finite and non-negative"))
    return tolerance
end

function _model_affine_point(coefficients,translation,p,caller,index::Int)
    a11,a21,a31,a12,a22,a32,a13,a23,a33=coefficients
    x,y,z=p
    return (
        _affine_coordinate(0.0,translation[1],a11,a12,a13,x,y,z,
                           0.0,0.0,0.0,caller,index),
        _affine_coordinate(0.0,translation[2],a21,a22,a23,x,y,z,
                           0.0,0.0,0.0,caller,index),
        _affine_coordinate(0.0,translation[3],a31,a32,a33,x,y,z,
                           0.0,0.0,0.0,caller,index),
    )
end

@inline _model_point_distance(a,b)=
    hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3])

function _model_curve_length(m::GeoModel,curve::Int,caller::AbstractString)
    a,b=m.curves[curve];p=m.points[a];q=m.points[b]
    length1=hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3])
    (isfinite(length1) && length1>0) || throw(ArgumentError(
        "$caller: Curve[$curve] must have finite positive geometric length"))
    return length1
end

"""
    set_periodic!(model, dim, slave_entities, master_entities, affine;
                  atol=1e-12)

Persist disjoint affine relations between equally sized lists of straight native
curves. `affine` maps each master curve to its slave in Gmsh row-major 4×4
order. This slice supports dimension 1; surface and volume periodicity are
explicit blockers. Both curves must bound the same planar surface when meshed.
Curve endpoints must be disjoint and agree with the affine map within `atol`.
The update is atomic.
"""
function set_periodic!(m::GeoModel,dim,slave_entities,master_entities,affine;
                       atol=1e-12)
    caller="set_periodic!"
    d=_dimension(dim,caller)
    d==1 || throw(ArgumentError(
        "$caller: only straight Curve periodicity (dimension 1) is implemented"))
    slaves=_periodic_entity_tags(slave_entities,caller,"slave")
    masters=_periodic_entity_tags(master_entities,caller,"master")
    length(slaves)==length(masters) || throw(ArgumentError(
        "$caller: slave and master entity counts differ"))
    isempty(slaves) && throw(ArgumentError(
        "$caller: need at least one slave/master curve pair"))
    length(unique(slaves))==length(slaves) || throw(ArgumentError(
        "$caller: slave curve tags must be unique"))
    participants=vcat(slaves,masters)
    length(unique(participants))==length(participants) || throw(ArgumentError(
        "$caller: every curve may participate in at most one periodic relation"))
    used=Set{Int}()
    for constraint in values(m.periodic)
        push!(used,Int(constraint.slave_entity))
        push!(used,Int(constraint.master_entity))
    end
    overlap=sort!(collect(intersect(Set(participants),used)))
    isempty(overlap) || throw(ArgumentError(
        "$caller: Curve[$(first(overlap))] already participates in a periodic relation"))
    tolerance=_model_periodic_tolerance(atol,caller)
    coefficients,translation,row_major=_transform_homogeneous(
        affine,caller;name="affine transform")
    pending=ModelPeriodicConstraint[]
    for (pair_index,(slave,master)) in enumerate(zip(slaves,masters))
        haskey(m.curves,slave) || throw(ArgumentError(
            "$caller: unknown slave Curve[$slave]"))
        haskey(m.curves,master) || throw(ArgumentError(
            "$caller: unknown master Curve[$master]"))
        _model_curve_length(m,slave,caller)
        _model_curve_length(m,master,caller)
        slave_points=m.curves[slave];master_points=m.curves[master]
        isempty(intersect(Set(slave_points),Set(master_points))) || throw(ArgumentError(
            "$caller: Curve[$slave] and Curve[$master] must have disjoint endpoints"))
        slave_start=m.points[slave_points[1]];slave_stop=m.points[slave_points[2]]
        mapped_start=_model_affine_point(
            coefficients,translation,m.points[master_points[1]],caller,2pair_index-1)
        mapped_stop=_model_affine_point(
            coefficients,translation,m.points[master_points[2]],caller,2pair_index)
        forward=max(_model_point_distance(mapped_start,slave_start),
                    _model_point_distance(mapped_stop,slave_stop))
        reversed=max(_model_point_distance(mapped_start,slave_stop),
                     _model_point_distance(mapped_stop,slave_start))
        mismatch=min(forward,reversed)
        mismatch<=tolerance || throw(ArgumentError(
            "$caller: affine map misses slave Curve[$slave] endpoints by $mismatch"))
        push!(pending,ModelPeriodicConstraint(
            d,Int32(slave),Int32(master),row_major,reversed<forward,tolerance))
    end
    for constraint in pending
        m.periodic[(constraint.dim,Int(constraint.slave_entity))]=constraint
    end
    return nothing
end

"""
    model_periodic_constraints(model) -> Vector{ModelPeriodicConstraint}

Return the model's immutable periodic-curve constraints in deterministic
slave-entity order.
"""
function model_periodic_constraints(m::GeoModel)
    return sort!(collect(values(m.periodic));
                 by=constraint->(constraint.dim,constraint.slave_entity))
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

function _periodic_curve_point(m::GeoModel,curve::Int,parameter::Float64,
                               caller::AbstractString)
    a,b=m.curves[curve];p=m.points[a];q=m.points[b]
    point=ntuple(3) do axis
        _affine_coordinate(
            p[axis],0.0,parameter,0.0,0.0,q[axis],0.0,0.0,
            p[axis],0.0,0.0,caller,curve)
    end
    all(isfinite,point) || throw(ArgumentError(
        "$caller: Curve[$curve] periodic subdivision is not Float64-representable"))
    return point
end

function _surface_curve_parameters(forced,curve::Int,signed::Int)
    parameters=get(forced,curve,nothing)
    parameters===nothing && return signed>0 ? (0.0,) : (1.0,)
    return signed>0 ? @view(parameters[1:end-1]) :
                      Iterators.reverse(@view(parameters[2:end]))
end

function _surface_pslg(m::GeoModel,t::Int,forced,caller::AbstractString)
    xs=Float64[];ys=Float64[];segs=Tuple{Int,Int}[]
    index=Dict{Int,Int}()
    for loop_id in m.surfaces[t]
        loop_idx=Int[]
        for signed in m.loops[loop_id]
            curve=abs(signed);a,b=m.curves[curve]
            for parameter in _surface_curve_parameters(forced,curve,signed)
                vertex=if parameter==0
                    _add_surface_point!(xs,ys,index,m,a,caller)
                elseif parameter==1
                    _add_surface_point!(xs,ys,index,m,b,caller)
                else
                    point=_periodic_curve_point(m,curve,parameter,caller)
                    abs(point[3])<=1e-12 || throw(ArgumentError(
                        "$caller: Curve[$curve] subdivision is not planar in z=0"))
                    push!(xs,point[1]);push!(ys,point[2]);length(xs)
                end
                push!(loop_idx,vertex)
            end
        end
        nloop=length(loop_idx)
        nloop>=3 || throw(ArgumentError(
            "$caller: Loop[$loop_id] needs at least three points"))
        for k in 1:nloop
            push!(segs,(loop_idx[k],loop_idx[mod1(k+1,nloop)]))
        end
    end
    embedded=get(m.embeds,(2,t),NTuple{2,Int}[])
    internal=Tuple{Int,Int}[]
    for (edim,etag) in embedded
        if edim==0
            _add_surface_point!(xs,ys,index,m,etag,caller)
        elseif edim==1
            haskey(m.curves,etag) || throw(ArgumentError(
                "$caller: unknown embedded Curve[$etag]"))
            a,b=m.curves[etag]
            ia=_add_surface_point!(xs,ys,index,m,a,caller)
            ib=_add_surface_point!(xs,ys,index,m,b,caller)
            ia==ib && throw(ArgumentError(
                "$caller: embedded Curve[$etag] has coincident endpoints"))
            push!(internal,(ia,ib))
        else
            throw(ArgumentError(
                "$caller: unsupported embedding dimension $edim"))
        end
    end
    return xs,ys,segs,index,embedded,internal
end

function _surface_boundary_topology(mesh::Mesh,caller::AbstractString)
    edge_counts=Dict{Tuple{Int32,Int32},Int}()
    @inbounds for cell in axes(mesh.tris,2),slots in ((1,2),(2,3),(3,1))
        a=mesh.tris[slots[1],cell];b=mesh.tris[slots[2],cell]
        key=a<b ? (a,b) : (b,a)
        edge_counts[key]=get(edge_counts,key,0)+1
    end
    boundary=falses(nnodes(mesh))
    edges=Set{Tuple{Int32,Int32}}()
    for ((a,b),count) in edge_counts
        count<=2 || throw(ArgumentError(
            "$caller: mesh edge ($a,$b) has $count incident triangles"))
        count==1 || continue
        push!(edges,(a,b))
        boundary[a]=true;boundary[b]=true
    end
    return boundary,edges
end

function _curve_parameter_nodes(m::GeoModel,mesh::Mesh,curve::Int,boundary,
                                boundary_edges,atol::Float64,
                                caller::AbstractString)
    a,b=m.curves[curve];p=m.points[a];q=m.points[b]
    vx=q[1]-p[1];vy=q[2]-p[2]
    length2=muladd(vx,vx,vy*vy)
    (isfinite(length2) && length2>0) || throw(ArgumentError(
        "$caller: Curve[$curve] has an unusable planar length"))
    length1=sqrt(length2)
    scale=max(1.0,hypot(p[1],p[2]),hypot(q[1],q[2]))
    geometric_tolerance=max(atol,128eps(Float64)*scale)
    entries=Tuple{Float64,Int}[]
    @inbounds for node in 1:nnodes(mesh)
        boundary[node] || continue
        wx=mesh.coords[1,node]-p[1];wy=mesh.coords[2,node]-p[2]
        cross=muladd(vx,wy,-vy*wx)
        abs(cross)<=geometric_tolerance*length1 || continue
        parameter=muladd(wx,vx,wy*vy)/length2
        -geometric_tolerance/length1<=parameter<=
            1+geometric_tolerance/length1 || continue
        push!(entries,(clamp(parameter,0.0,1.0),node))
    end
    sort!(entries;by=first)
    length(entries)>=2 || throw(ErrorException(
        "$caller: Curve[$curve] is not represented by two boundary nodes"))
    parameter_tolerance=max(128eps(Float64),geometric_tolerance/length1)
    first(entries)[1]<=parameter_tolerance &&
        1-last(entries)[1]<=parameter_tolerance || throw(ErrorException(
            "$caller: Curve[$curve] boundary chain does not reach both endpoints"))
    for index in 1:(length(entries)-1)
        first_node=Int32(entries[index][2])
        second_node=Int32(entries[index+1][2])
        first_node!=second_node || throw(ErrorException(
            "$caller: Curve[$curve] repeats a boundary node"))
        key=first_node<second_node ? (first_node,second_node) :
                                     (second_node,first_node)
        key in boundary_edges || throw(ErrorException(
            "$caller: Curve[$curve] boundary nodes do not form a mesh-edge chain"))
    end
    return entries,parameter_tolerance
end

function _insert_periodic_parameter!(parameters::Vector{Float64},value::Float64,
                                     tolerance::Float64)
    candidate=clamp(value,0.0,1.0)
    any(existing->abs(existing-candidate)<=tolerance,parameters) && return false
    push!(parameters,candidate);sort!(parameters)
    return true
end

function _surface_periodic_constraints(m::GeoModel,t::Int,
                                       caller::AbstractString)
    boundary_curves=Set{Int}()
    for loop in m.surfaces[t],signed in m.loops[loop]
        push!(boundary_curves,abs(signed))
    end
    constraints=ModelPeriodicConstraint[]
    for constraint in model_periodic_constraints(m)
        slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
        slave_present=slave in boundary_curves
        master_present=master in boundary_curves
        slave_present==master_present || throw(ArgumentError(
            "$caller: periodic Curve[$slave]/Curve[$master] relation " *
            "has only one entity on Surface[$t]"))
        slave_present && push!(constraints,constraint)
    end
    return constraints
end

function _synchronize_periodic_parameters!(forced,m::GeoModel,mesh::Mesh,
                                           constraints)
    boundary,boundary_edges=_surface_boundary_topology(
        mesh,"mesh_model_surface")
    changed=false
    for constraint in constraints
        slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
        master_entries,master_tolerance=_curve_parameter_nodes(
            m,mesh,master,boundary,boundary_edges,constraint.atol,
            "mesh_model_surface")
        slave_entries,slave_tolerance=_curve_parameter_nodes(
            m,mesh,slave,boundary,boundary_edges,constraint.atol,
            "mesh_model_surface")
        tolerance=max(master_tolerance,slave_tolerance)
        common=Float64[first(entry) for entry in master_entries]
        for (parameter,_) in slave_entries
            mapped=constraint.reversed ? 1-parameter : parameter
            _insert_periodic_parameter!(common,mapped,tolerance)
        end
        master_forced=get!(()->Float64[0,1],forced,master)
        slave_forced=get!(()->Float64[0,1],forced,slave)
        for parameter in common
            changed|=_insert_periodic_parameter!(
                master_forced,parameter,tolerance)
            slave_parameter=constraint.reversed ? 1-parameter : parameter
            changed|=_insert_periodic_parameter!(
                slave_forced,slave_parameter,tolerance)
        end
    end
    return changed
end

function _model_periodic_nodes(m::GeoModel,mesh::Mesh,
                               constraint::ModelPeriodicConstraint)
    boundary,boundary_edges=_surface_boundary_topology(
        mesh,"model_periodic_nodes")
    slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
    master_entries,master_tolerance=_curve_parameter_nodes(
        m,mesh,master,boundary,boundary_edges,constraint.atol,
        "model_periodic_nodes")
    slave_entries,slave_tolerance=_curve_parameter_nodes(
        m,mesh,slave,boundary,boundary_edges,constraint.atol,
        "model_periodic_nodes")
    ordered_slave=constraint.reversed ? reverse(slave_entries) : slave_entries
    length(master_entries)==length(ordered_slave) || throw(ErrorException(
        "model_periodic_nodes: Curve[$slave] and Curve[$master] node counts differ"))
    tolerance=max(master_tolerance,slave_tolerance)
    master_nodes=Vector{Int32}(undef,length(master_entries))
    slave_nodes=similar(master_nodes)
    for index in eachindex(master_entries)
        master_parameter,master_node=master_entries[index]
        slave_parameter,slave_node=ordered_slave[index]
        mapped=constraint.reversed ? 1-slave_parameter : slave_parameter
        abs(master_parameter-mapped)<=tolerance || throw(ErrorException(
            "model_periodic_nodes: Curve[$slave]/Curve[$master] parameters differ"))
        master_nodes[index]=Int32(master_node);slave_nodes[index]=Int32(slave_node)
    end
    return (master_entity=master,slave_nodes=slave_nodes,
            master_nodes=master_nodes,affine=constraint.affine)
end

function _model_mapping_matches(mesh::Mesh,constraint::ModelPeriodicConstraint,
                                mapping,caller::AbstractString;exact::Bool)
    coefficients,translation,_=_transform_homogeneous(
        constraint.affine,caller;name="stored affine transform")
    for (master_node,slave_node) in zip(mapping.master_nodes,
                                        mapping.slave_nodes)
        master=(mesh.coords[1,master_node],mesh.coords[2,master_node],
                mesh.coords[3,master_node])
        expected=_model_affine_point(
            coefficients,translation,master,caller,Int(master_node))
        slave=(mesh.coords[1,slave_node],mesh.coords[2,slave_node],
               mesh.coords[3,slave_node])
        if exact
            slave==expected || return false
        elseif _model_point_distance(slave,expected)>constraint.atol
            throw(ArgumentError(
                "$caller: slave node $slave_node is not affine(master node $master_node)"))
        end
    end
    return true
end

function _snap_surface_periodic(m::GeoModel,mesh::Mesh,constraints,
                                caller::AbstractString)
    output=mesh
    for constraint in constraints
        mapping=_model_periodic_nodes(m,output,constraint)
        output=periodic_identify_affine(
            output,constraint.affine,mapping.master_nodes,mapping.slave_nodes;
            atol=constraint.atol)
    end
    for constraint in constraints
        mapping=_model_periodic_nodes(m,output,constraint)
        _model_mapping_matches(
            output,constraint,mapping,caller;exact=true) || throw(ErrorException(
            "$caller: stored periodic constraints do not share an exact boundary-node solution"))
    end
    return output
end

"""
    model_periodic_nodes(model, mesh, dim, slave_entity)

Return the master entity, compact slave/master node arrays, and affine transform
for one meshed straight-curve relation as a named tuple. The mesh must contain a
synchronized boundary discretization produced by [`mesh_model_surface`](@ref).
"""
function model_periodic_nodes(m::GeoModel,mesh::Mesh,dim,slave_entity)
    caller="model_periodic_nodes"
    d=_dimension(dim,caller);slave=_tag(slave_entity,caller,d)
    constraint=get(m.periodic,(d,slave),nothing)
    constraint===nothing && throw(ArgumentError(
        "$caller: no periodic relation for entity ($d,$slave)"))
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(diagnostic.messages,"; ")))
    mapping=_model_periodic_nodes(m,mesh,constraint)
    _model_mapping_matches(mesh,constraint,mapping,caller;exact=false)
    return mapping
end

function _model_projection_physical_tags(m::GeoModel,dim::Int,entity::Int)
    tags=Int32[]
    for ((group_dim,physical),members) in m.physical
        group_dim==dim && entity in members && push!(tags,Int32(physical))
    end
    sort!(tags)
    return tags
end

@inline _model_projection_legacy_tag(tags::Vector{Int32})=
    isempty(tags) ? Int32(0) : first(tags)

function _model_projection_bbox(mesh::Mesh,nodes,caller::AbstractString)
    isempty(nodes) && throw(ErrorException(
        "$caller: internal entity projection has no nodes"))
    first_node=Int(first(nodes))
    xlo=xhi=mesh.coords[1,first_node]
    ylo=yhi=mesh.coords[2,first_node]
    zlo=zhi=mesh.coords[3,first_node]
    for raw_node in Iterators.drop(nodes,1)
        node=Int(raw_node)
        x=mesh.coords[1,node];y=mesh.coords[2,node];z=mesh.coords[3,node]
        xlo=min(xlo,x);xhi=max(xhi,x)
        ylo=min(ylo,y);yhi=max(yhi,y)
        zlo=min(zlo,z);zhi=max(zhi,z)
    end
    return (xlo,ylo,zlo,xhi,yhi,zhi)
end

function _model_projection_surface_curves(m::GeoModel,surface::Int)
    curves=Int[]
    for loop in m.surfaces[surface],signed_curve in m.loops[loop]
        push!(curves,abs(signed_curve))
    end
    return sort!(unique!(curves))
end

function _model_projection_surface_boundaries(m::GeoModel,surface::Int)
    boundaries=Int32[]
    # The selected surface traverses each hole opposite its stored loop direction.
    for (loop_index,loop) in pairs(m.surfaces[surface])
        curves=m.loops[loop]
        if loop_index==1
            append!(boundaries,Int32.(curves))
        else
            for signed_curve in Iterators.reverse(curves)
                push!(boundaries,-Int32(signed_curve))
            end
        end
    end
    return boundaries
end

function _model_projection_point!(point_nodes,node_points,point::Int,node::Int,
                                  caller::AbstractString)
    compact=Int32(node)
    existing=get(point_nodes,point,Int32(0))
    (existing==0 || existing==compact) || throw(ArgumentError(
        "$caller: Point[$point] maps to multiple mesh nodes"))
    other=get(node_points,compact,0)
    (other==0 || other==point) || throw(ArgumentError(
        "$caller: Points[$other] and [$point] map to mesh node $node"))
    point_nodes[point]=compact
    node_points[compact]=point
    return nothing
end

function _model_projection_external_elements(blocks)
    output=Vector{Vector{UInt64}}(undef,length(blocks))
    next_tag=UInt64(1)
    for (block_index,block) in pairs(blocks)
        count=length(block.tags)
        tags=Vector{UInt64}(undef,count)
        for cell in 1:count
            tags[cell]=next_tag
            next_tag+=UInt64(1)
        end
        output[block_index]=tags
    end
    return output
end

function _model_projection_triangle_edges(mesh::Mesh)
    edges=Set{Tuple{Int32,Int32}}()
    @inbounds for cell in axes(mesh.tris,2),slots in ((1,2),(2,3),(3,1))
        first_node=mesh.tris[slots[1],cell]
        second_node=mesh.tris[slots[2],cell]
        push!(edges,first_node<second_node ? (first_node,second_node) :
                                             (second_node,first_node))
    end
    return edges
end

function _model_projection_embedded_curve_nodes(
    m::GeoModel,mesh::Mesh,curve::Int,mesh_edges,
    atol::Float64,caller::AbstractString)
    start_point,stop_point=m.curves[curve]
    first_coordinate=m.points[start_point]
    last_coordinate=m.points[stop_point]
    vx=last_coordinate[1]-first_coordinate[1]
    vy=last_coordinate[2]-first_coordinate[2]
    length2=muladd(vx,vx,vy*vy)
    (isfinite(length2) && length2>0) || throw(ArgumentError(
        "$caller: embedded Curve[$curve] has an unusable planar length"))
    length1=sqrt(length2)
    scale=max(1.0,hypot(first_coordinate[1],first_coordinate[2]),
                    hypot(last_coordinate[1],last_coordinate[2]))
    geometric_tolerance=max(atol,128eps(Float64)*scale)
    parameter_tolerance=max(128eps(Float64),geometric_tolerance/length1)
    entries=Tuple{Float64,Int}[]
    @inbounds for node in 1:nnodes(mesh)
        wx=mesh.coords[1,node]-first_coordinate[1]
        wy=mesh.coords[2,node]-first_coordinate[2]
        cross=muladd(vx,wy,-vy*wx)
        abs(cross)<=geometric_tolerance*length1 || continue
        parameter=muladd(wx,vx,wy*vy)/length2
        -parameter_tolerance<=parameter<=1+parameter_tolerance || continue
        push!(entries,(clamp(parameter,0.0,1.0),node))
    end
    sort!(entries;by=entry->(entry[1],entry[2]))
    length(entries)>=2 || throw(ArgumentError(
        "$caller: embedded Curve[$curve] is not represented by two mesh nodes"))
    first(entries)[1]<=parameter_tolerance &&
        1-last(entries)[1]<=parameter_tolerance || throw(ArgumentError(
        "$caller: embedded Curve[$curve] mesh chain does not reach both endpoints"))
    for index in 1:(length(entries)-1)
        first_parameter,first_raw=entries[index]
        second_parameter,second_raw=entries[index+1]
        second_parameter-first_parameter>parameter_tolerance || throw(ArgumentError(
            "$caller: embedded Curve[$curve] maps multiple mesh nodes to one parameter"))
        first_node=Int32(first_raw);second_node=Int32(second_raw)
        edge=first_node<second_node ? (first_node,second_node) :
                                      (second_node,first_node)
        edge in mesh_edges || throw(ArgumentError(
            "$caller: embedded Curve[$curve] nodes do not form a mesh-edge chain"))
    end
    return entries
end

function _model_projection_embedded_point_node(
    m::GeoModel,mesh::Mesh,point::Int,atol::Float64,caller::AbstractString)
    coordinate=m.points[point]
    scale=max(1.0,hypot(coordinate...))
    tolerance=max(atol,128eps(Float64)*scale)
    matches=Int[]
    @inbounds for node in 1:nnodes(mesh)
        distance=hypot(mesh.coords[1,node]-coordinate[1],
                       mesh.coords[2,node]-coordinate[2],
                       mesh.coords[3,node]-coordinate[3])
        distance<=tolerance && push!(matches,node)
    end
    isempty(matches) && throw(ArgumentError(
        "$caller: embedded Point[$point] is not a mesh node"))
    length(matches)==1 || throw(ArgumentError(
        "$caller: embedded Point[$point] matches $(length(matches)) mesh nodes; " *
        "expected exactly one"))
    return only(matches)
end

function _model_projection_periodic_links(m::GeoModel,mesh::Mesh,point_nodes,
                                          constraints,caller::AbstractString)
    curve_links=MixedPeriodicLink[]
    outgoing=Dict{Int,Vector{Tuple{Int,NTuple{16,Float64}}}}()
    indegree=Dict{Int,Int}()
    periodic_points=Set{Int}()
    for constraint in constraints
        mapping=try
            _model_periodic_nodes(m,mesh,constraint)
        catch err
            err isa InterruptException && rethrow()
            (err isa ArgumentError || err isa ErrorException) || rethrow()
            throw(ArgumentError(
                "$caller: periodic Curve[$(constraint.slave_entity)] mapping " *
                "is incompatible with the input mesh — " * sprint(showerror,err)))
        end
        _model_mapping_matches(mesh,constraint,mapping,caller;exact=true) ||
            throw(ArgumentError(
                "$caller: periodic Curve[$(constraint.slave_entity)] nodes are not exactly snapped"))
        push!(curve_links,MixedPeriodicLink(
            1,constraint.slave_entity,constraint.master_entity,
            mapping.slave_nodes,mapping.master_nodes;affine=constraint.affine))
        slave_start,slave_stop=m.curves[Int(constraint.slave_entity)]
        master_start,master_stop=m.curves[Int(constraint.master_entity)]
        endpoint_pairs=constraint.reversed ?
            ((slave_stop,master_start),(slave_start,master_stop)) :
            ((slave_start,master_start),(slave_stop,master_stop))
        for (slave_point,master_point) in endpoint_pairs
            push!(periodic_points,slave_point);push!(periodic_points,master_point)
            edges=get!(Vector{Tuple{Int,NTuple{16,Float64}}},
                       outgoing,master_point)
            edge=(slave_point,constraint.affine)
            if !(edge in edges)
                push!(edges,edge)
                indegree[slave_point]=get(indegree,slave_point,0)+1
            end
            get!(indegree,master_point,0)
        end
    end
    for edges in values(outgoing)
        sort!(edges;by=edge->(edge[1],edge[2]))
    end

    # MSH periodic metadata permits only one relation per slave entity. Choose
    # a deterministic directed spanning forest through shared periodic corners,
    # matching Gmsh's endpoint-link convention without discarding curve links.
    visited=Set{Int}()
    parents=Dict{Int,Tuple{Int,NTuple{16,Float64}}}()
    roots=sort!(Int[point for point in periodic_points
                    if get(indegree,point,0)==0])
    append!(roots,sort!(collect(setdiff(periodic_points,Set(roots)))))
    for root in roots
        root in visited && continue
        push!(visited,root)
        queue=Int[root];head=1
        while head<=length(queue)
            master_point=queue[head];head+=1
            for (slave_point,affine) in get(
                    outgoing,master_point,
                    Tuple{Int,NTuple{16,Float64}}[])
                slave_point in visited && continue
                push!(visited,slave_point)
                parents[slave_point]=(master_point,affine)
                push!(queue,slave_point)
            end
        end
    end
    visited==periodic_points || throw(ErrorException(
        "$caller: internal periodic endpoint traversal was incomplete"))

    links=MixedPeriodicLink[]
    for slave_point in sort!(collect(keys(parents)))
        master_point,affine=parents[slave_point]
        haskey(point_nodes,slave_point) && haskey(point_nodes,master_point) ||
            throw(ErrorException(
                "$caller: periodic endpoint is absent from the surface projection"))
        push!(links,MixedPeriodicLink(
            0,slave_point,master_point,
            Int32[point_nodes[slave_point]],Int32[point_nodes[master_point]];
            affine=affine))
    end
    append!(links,curve_links)
    return links
end

"""
    model_to_mixed(model, mesh, surface_tag) -> MixedMesh

Project a native planar surface mesh and its geometry ownership into an owned
[`MixedMesh`](@ref). The result contains point, boundary/embedded-line, and
triangle blocks; MSH2 elementary entities; MSH4 point/curve/surface
classification; all projected physical memberships and names; embedded-curve
ownership; and the surface's stored boundary-periodic curve links. It can be
written with [`Tessella.Elements.write_mixed_msh`](@ref) without dropping the
supported metadata. MSH4 retains every physical membership and Gmsh's
curve-in-surface relation. Gmsh 4.15.2 does not serialize Point-In-Surface as an
entity relation, but the classified point node and point element are retained.
MSH2 retains the lowest physical tag as its single legacy element membership.

`mesh` must be a validated, segment-free triangle mesh whose boundary and
embedded chains represent the selected [`mesh_model_surface`](@ref) geometry.
Stored boundary-periodic curves must already be exactly snapped. Periodic
embedded curves are an explicit blocker. Periodic endpoint relations are emitted
as a deterministic spanning forest when curve directions share corners,
satisfying the MSH one-master-per-slave entity constraint while retaining every
boundary-curve link.
"""
function model_to_mixed(m::GeoModel,mesh::Mesh,surface_tag::Integer)
    caller="model_to_mixed"
    surface=_tag(surface_tag,caller,2)
    haskey(m.surfaces,surface) || throw(ArgumentError(
        "$caller: unknown Surface[$surface]"))
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(diagnostic.messages,"; ")))
    nsegs(mesh)==0 || throw(ArgumentError(
        "$caller: input must not contain explicit segment cells"))
    ntets(mesh)==0 || throw(ArgumentError(
        "$caller: input must be a surface mesh without tetrahedra"))
    ntris(mesh)>0 || throw(ArgumentError(
        "$caller: input must contain triangle cells"))
    all(iszero,mesh.tri_tag) || throw(ArgumentError(
        "$caller: input triangle tags must be zero; physical ownership comes from the model"))

    embedded_points=Int[];embedded_curves=Int[]
    for (embedded_dim,embedded_tag) in
            get(m.embeds,(2,surface),NTuple{2,Int}[])
        if embedded_dim==0
            push!(embedded_points,embedded_tag)
        elseif embedded_dim==1
            push!(embedded_curves,embedded_tag)
        else
            throw(ArgumentError(
                "$caller: unsupported Surface[$surface] embedding dimension $embedded_dim"))
        end
    end
    sort!(unique!(embedded_points));sort!(unique!(embedded_curves))

    constraints=_surface_periodic_constraints(m,surface,caller)
    embedded_curve_set=Set(embedded_curves)
    for constraint in model_periodic_constraints(m)
        slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
        (slave in embedded_curve_set || master in embedded_curve_set) &&
            throw(ArgumentError(
                "$caller: periodic embedded curves are not supported " *
                "(Curve[$slave]/Curve[$master])"))
    end
    geometric_tolerance=max(1e-12,
        isempty(constraints) ? 0.0 : maximum(c.atol for c in constraints))
    for node in axes(mesh.coords,2)
        abs(mesh.coords[3,node])<=geometric_tolerance || throw(ArgumentError(
            "$caller: input node $node has z=$(mesh.coords[3,node]), outside " *
            "the z=0 tolerance $geometric_tolerance"))
    end
    curve_tags=_model_projection_surface_curves(m,surface)
    isempty(curve_tags) && throw(ArgumentError(
        "$caller: Surface[$surface] has no boundary curves"))
    any(curve->curve in curve_tags,embedded_curves) && throw(ArgumentError(
        "$caller: a surface curve cannot be both bounding and embedded"))
    all_curve_tags=sort!(unique!(vcat(curve_tags,embedded_curves)))
    for curve in all_curve_tags,point in m.curves[curve]
        abs(m.points[point][3])<=1e-12 || throw(ArgumentError(
            "$caller: Point[$point] is not planar in z=0"))
    end
    for point in embedded_points
        abs(m.points[point][3])<=1e-12 || throw(ArgumentError(
            "$caller: embedded Point[$point] is not planar in z=0"))
    end
    boundary,boundary_edges=_surface_boundary_topology(mesh,caller)
    isempty(boundary_edges) && throw(ArgumentError(
        "$caller: input mesh has no boundary edges"))
    curve_entries=Dict{Int,Vector{Tuple{Float64,Int}}}()
    point_nodes=Dict{Int,Int32}()
    node_points=Dict{Int32,Int}()
    line_cells=NTuple{2,Int32}[]
    line_entities=Int32[]
    claimed_edges=Set{Tuple{Int32,Int32}}()
    for curve in curve_tags
        entries=try
            first(_curve_parameter_nodes(
                m,mesh,curve,boundary,boundary_edges,
                geometric_tolerance,caller))
        catch err
            err isa InterruptException && rethrow()
            err isa ArgumentError && rethrow()
            throw(ArgumentError(
                "$caller: input mesh does not represent Curve[$curve] — " *
                sprint(showerror,err)))
        end
        curve_entries[curve]=entries
        start_point,stop_point=m.curves[curve]
        _model_projection_point!(
            point_nodes,node_points,start_point,entries[1][2],caller)
        _model_projection_point!(
            point_nodes,node_points,stop_point,entries[end][2],caller)
        for index in 1:(length(entries)-1)
            first_node=Int32(entries[index][2])
            second_node=Int32(entries[index+1][2])
            edge=first_node<second_node ? (first_node,second_node) :
                                          (second_node,first_node)
            edge in claimed_edges && throw(ArgumentError(
                "$caller: boundary edge $edge belongs to multiple model curves"))
            push!(claimed_edges,edge)
            push!(line_cells,(first_node,second_node))
            push!(line_entities,Int32(curve))
        end
    end
    claimed_edges==boundary_edges || throw(ArgumentError(
        "$caller: model curves cover $(length(claimed_edges)) of " *
        "$(length(boundary_edges)) mesh boundary edges"))

    mesh_edges=_model_projection_triangle_edges(mesh)
    for curve in embedded_curves
        entries=_model_projection_embedded_curve_nodes(
            m,mesh,curve,mesh_edges,geometric_tolerance,caller)
        curve_entries[curve]=entries
        start_point,stop_point=m.curves[curve]
        _model_projection_point!(
            point_nodes,node_points,start_point,entries[1][2],caller)
        _model_projection_point!(
            point_nodes,node_points,stop_point,entries[end][2],caller)
        for index in 1:(length(entries)-1)
            first_node=Int32(entries[index][2])
            second_node=Int32(entries[index+1][2])
            edge=first_node<second_node ? (first_node,second_node) :
                                          (second_node,first_node)
            edge in claimed_edges && throw(ArgumentError(
                "$caller: mesh edge $edge belongs to multiple model curves"))
            push!(claimed_edges,edge)
            push!(line_cells,(first_node,second_node))
            push!(line_entities,Int32(curve))
        end
    end
    for point in embedded_points
        node=_model_projection_embedded_point_node(
            m,mesh,point,geometric_tolerance,caller)
        _model_projection_point!(
            point_nodes,node_points,point,node,caller)
    end

    node_entities=fill((2,Int32(surface)),nnodes(mesh))
    for (point,node) in point_nodes
        node_entities[node]=(0,Int32(point))
    end
    for curve in all_curve_tags
        entries=curve_entries[curve]
        if length(entries)>2
            for index in 2:(length(entries)-1)
                node=entries[index][2]
                owner=node_entities[node]
                if owner[1]==2
                    node_entities[node]=(1,Int32(curve))
                elseif owner[1]==0
                    # A classified embedded point takes precedence over a curve
                    # that passes through the same mesh node.
                    continue
                elseif owner!=(1,Int32(curve))
                    throw(ArgumentError(
                        "$caller: mesh node $node belongs to Curves[$(owner[2])] and [$curve]"))
                end
            end
        end
    end
    for node in eachindex(boundary)
        boundary[node] && node_entities[node][1]==2 && throw(ArgumentError(
            "$caller: boundary mesh node $node has no model entity"))
    end

    links=_model_projection_periodic_links(
        m,mesh,point_nodes,constraints,caller)

    point_tags=sort!(collect(keys(point_nodes)))
    total_elements=try
        subtotal=Base.checked_add(length(point_tags),length(line_cells))
        Base.checked_add(subtotal,ntris(mesh))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: projected element count overflows Int"))
    end
    total_elements<=typemax(Int32) || throw(ArgumentError(
        "$caller: projected element count exceeds the Int32 MSH limit"))
    point_matrix=reshape(Int32[point_nodes[tag] for tag in point_tags],1,:)
    point_physical=Vector{Int32}(undef,length(point_tags))
    entities=Dict{Tuple{Int,Int},MixedEntity}()
    projected_groups=Set{Tuple{Int,Int}}()
    for (index,point) in pairs(point_tags)
        tags=_model_projection_physical_tags(m,0,point)
        point_physical[index]=_model_projection_legacy_tag(tags)
        union!(projected_groups,((0,Int(tag)) for tag in tags))
        node=point_nodes[point]
        coordinate=(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
        entities[(0,point)]=MixedEntity(
            0,point,coordinate;physical_tags=tags)
    end

    line_matrix=Matrix{Int32}(undef,2,length(line_cells))
    line_physical=Vector{Int32}(undef,length(line_cells))
    curve_physical=Dict{Int,Vector{Int32}}()
    for curve in all_curve_tags
        tags=_model_projection_physical_tags(m,1,curve)
        curve_physical[curve]=tags
        union!(projected_groups,((1,Int(tag)) for tag in tags))
        nodes=Int32[entry[2] for entry in curve_entries[curve]]
        start_point,stop_point=m.curves[curve]
        entities[(1,curve)]=MixedEntity(
            1,curve,_model_projection_bbox(mesh,nodes,caller);
            physical_tags=tags,
            boundaries=Int32[Int32(start_point),-Int32(stop_point)])
    end
    for (cell,(first_node,second_node)) in pairs(line_cells)
        line_matrix[1,cell]=first_node;line_matrix[2,cell]=second_node
        line_physical[cell]=_model_projection_legacy_tag(
            curve_physical[Int(line_entities[cell])])
    end

    surface_physical=_model_projection_physical_tags(m,2,surface)
    union!(projected_groups,((2,Int(tag)) for tag in surface_physical))
    surface_boundaries=_model_projection_surface_boundaries(m,surface)
    entities[(2,surface)]=MixedEntity(
        2,surface,_model_projection_bbox(mesh,axes(mesh.coords,2),caller);
        physical_tags=surface_physical,boundaries=surface_boundaries,
        embedded_curves=embedded_curves)
    triangle_physical=fill(
        _model_projection_legacy_tag(surface_physical),ntris(mesh))

    blocks=ElementBlock[
        ElementBlock(15,point_matrix,point_physical),
        ElementBlock(1,line_matrix,line_physical),
        ElementBlock(2,mesh.tris,triangle_physical),
    ]
    block_entities=Vector{Int32}[
        Int32.(point_tags),
        line_entities,
        fill(Int32(surface),ntris(mesh)),
    ]
    external_node_tags=UInt64.(1:nnodes(mesh))
    external_element_tags=_model_projection_external_elements(blocks)
    node_parametric=Union{Nothing,Vector{Float64}}[
        nothing for _ in 1:nnodes(mesh)]
    data=MixedEntityData(
        entities;node_entities=node_entities,node_parametric=node_parametric,
        external_node_tags=external_node_tags,block_entities=block_entities,
        external_element_tags=external_element_tags)
    names=Dict{Tuple{Int,Int},String}()
    for key in sort!(collect(projected_groups))
        haskey(m.physical_names,key) && (names[key]=m.physical_names[key])
    end
    output=MixedMesh(
        mesh.coords,blocks;physical_names=names,entity_data=data,
        elementary_entities=block_entities,periodic_links=links)
    output_diagnostic=validate(output)
    output_diagnostic.ok || throw(ErrorException(
        "$caller: invalid projected mixed mesh — " *
        join(output_diagnostic.messages,"; ")))
    return output
end

function _model_projection_tet_curve_nodes(
    m::GeoModel,mesh::Mesh,curve::Int,tet_edges,
    caller::AbstractString)
    start_point,stop_point=m.curves[curve]
    first_coordinate=m.points[start_point]
    last_coordinate=m.points[stop_point]
    vx=last_coordinate[1]-first_coordinate[1]
    vy=last_coordinate[2]-first_coordinate[2]
    vz=last_coordinate[3]-first_coordinate[3]
    length2=muladd(vx,vx,muladd(vy,vy,vz*vz))
    (isfinite(length2) && length2>0) || throw(ArgumentError(
        "$caller: Curve[$curve] has an unusable length"))
    length1=sqrt(length2)
    scale=max(1.0,hypot(first_coordinate...),hypot(last_coordinate...))
    geometric_tolerance=max(1e-12,128eps(Float64)*scale)
    parameter_tolerance=max(128eps(Float64),geometric_tolerance/length1)
    entries=Tuple{Float64,Int}[]
    @inbounds for node in 1:nnodes(mesh)
        wx=mesh.coords[1,node]-first_coordinate[1]
        wy=mesh.coords[2,node]-first_coordinate[2]
        wz=mesh.coords[3,node]-first_coordinate[3]
        cross_norm=hypot(vy*wz-vz*wy,vz*wx-vx*wz,vx*wy-vy*wx)
        cross_norm<=geometric_tolerance*length1 || continue
        parameter=muladd(wx,vx,muladd(wy,vy,wz*vz))/length2
        -parameter_tolerance<=parameter<=1+parameter_tolerance || continue
        push!(entries,(clamp(parameter,0.0,1.0),node))
    end
    sort!(entries;by=entry->(entry[1],entry[2]))
    length(entries)>=2 || throw(ArgumentError(
        "$caller: Curve[$curve] is not represented by two mesh nodes"))
    first(entries)[1]<=parameter_tolerance &&
        1-last(entries)[1]<=parameter_tolerance || throw(ArgumentError(
        "$caller: Curve[$curve] mesh chain does not reach both endpoints"))
    for index in 1:(length(entries)-1)
        first_parameter,first_raw=entries[index]
        second_parameter,second_raw=entries[index+1]
        second_parameter-first_parameter>parameter_tolerance || throw(ArgumentError(
            "$caller: Curve[$curve] maps multiple mesh nodes to one parameter"))
        first_node=Int32(first_raw);second_node=Int32(second_raw)
        edge=first_node<second_node ? (first_node,second_node) :
                                      (second_node,first_node)
        edge in tet_edges || throw(ArgumentError(
            "$caller: Curve[$curve] nodes do not form a tetrahedron-edge chain"))
    end
    return entries
end

@inline function _model_projection_face_key(face::NTuple{3,Int32})
    first_node,second_node,third_node=face
    first_node>second_node &&
        ((first_node,second_node)=(second_node,first_node))
    second_node>third_node &&
        ((second_node,third_node)=(third_node,second_node))
    first_node>second_node &&
        ((first_node,second_node)=(second_node,first_node))
    return (first_node,second_node,third_node)
end

function _model_projection_volume_surface_faces!(
    claimed_faces::Set{NTuple{3,Int32}},m::GeoModel,mesh::Mesh,
    surface::Int,caller::AbstractString)
    loops=m.surfaces[surface]
    length(loops)==1 || throw(ArgumentError(
        "$caller: holed embedded Surface[$surface] is not supported"))
    point_tags=_loop_points(m,only(loops))
    length(point_tags)>=3 || throw(ArgumentError(
        "$caller: embedded Surface[$surface] needs at least three points"))
    coordinates=[m.points[point] for point in point_tags]
    local_faces=Set{NTuple{3,Int32}}()
    output=NTuple{3,Int32}[]
    for index in 2:(length(coordinates)-1)
        first_coordinate=coordinates[1]
        second_coordinate=coordinates[index]
        third_coordinate=coordinates[index+1]
        faces,target,covered=_mesh_covering_faces3(
            mesh,first_coordinate,second_coordinate,third_coordinate)
        (isfinite(target) && target>0 && isfinite(covered) &&
         abs(covered-target)<=1e-6*target) || throw(ArgumentError(
            "$caller: embedded Surface[$surface] fan triangle $(index-1) " *
            "is not represented by tetrahedron faces"))
        for face in faces
            key=_model_projection_face_key(face)
            key in local_faces && throw(ArgumentError(
                "$caller: embedded Surface[$surface] fan triangles overlap " *
                "on mesh face $key"))
            key in claimed_faces && throw(ArgumentError(
                "$caller: mesh face $key belongs to multiple embedded surfaces"))
            push!(local_faces,key)
            push!(claimed_faces,key)
            push!(output,face)
        end
    end
    isempty(output) && throw(ArgumentError(
        "$caller: embedded Surface[$surface] has no tetrahedron faces"))
    return output
end

function _model_projection_volume_inventory(
    m::GeoModel,volume::Int,caller::AbstractString)
    point_tags=Int[];curve_tags=Int[];surface_tags=Int[]
    for (embedded_dim,embedded_tag) in
            get(m.embeds,(3,volume),NTuple{2,Int}[])
        if embedded_dim==0
            haskey(m.points,embedded_tag) || throw(ArgumentError(
                "$caller: unknown embedded Point[$embedded_tag]"))
            push!(point_tags,embedded_tag)
        elseif embedded_dim==1
            haskey(m.curves,embedded_tag) || throw(ArgumentError(
                "$caller: unknown embedded Curve[$embedded_tag]"))
            push!(curve_tags,embedded_tag)
        elseif embedded_dim==2
            haskey(m.surfaces,embedded_tag) || throw(ArgumentError(
                "$caller: unknown embedded Surface[$embedded_tag]"))
            isempty(get(m.embeds,(2,embedded_tag),NTuple{2,Int}[])) ||
                throw(ArgumentError(
                    "$caller: nested Point/Line-In-Surface constraints on embedded " *
                    "Surface[$embedded_tag] are not supported"))
            push!(surface_tags,embedded_tag)
            append!(curve_tags,_model_projection_surface_curves(m,embedded_tag))
        else
            throw(ArgumentError(
                "$caller: unsupported Volume[$volume] embedding dimension $embedded_dim"))
        end
    end
    sort!(unique!(curve_tags));sort!(unique!(surface_tags))
    for curve in curve_tags
        start_point,stop_point=m.curves[curve]
        push!(point_tags,start_point);push!(point_tags,stop_point)
    end
    sort!(unique!(point_tags))
    return point_tags,curve_tags,surface_tags
end

function _model_volume_to_mixed(m::GeoModel,mesh::Mesh,volume::Int)
    caller="model_to_mixed"
    haskey(m.volumes,volume) || throw(ArgumentError(
        "$caller: unknown Volume[$volume]"))
    isempty(m.volumes[volume]) || throw(ArgumentError(
        "$caller: explicit boundary-surface projection for Volume[$volume] " *
        "is not supported"))
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(diagnostic.messages,"; ")))
    nsegs(mesh)==0 || throw(ArgumentError(
        "$caller: volume input must not contain explicit segment cells"))
    ntris(mesh)==0 || throw(ArgumentError(
        "$caller: volume input must not contain explicit triangle cells"))
    ntets(mesh)>0 || throw(ArgumentError(
        "$caller: volume input must contain tetrahedron cells"))
    all(iszero,mesh.tet_tag) || throw(ArgumentError(
        "$caller: input tetrahedron tags must be zero; physical ownership comes from the model"))
    isempty(m.periodic) || throw(ArgumentError(
        "$caller: periodic relations are not supported in classified volume projection"))

    domain_surface=_volume_surface(m,volume,caller)
    fills_volume,fill_reason=_certify_surface_fill(domain_surface,mesh)
    fills_volume || throw(ArgumentError(
        "$caller: input tetrahedron mesh does not fill Volume[$volume] — $fill_reason"))

    point_tags,curve_tags,surface_tags=
        _model_projection_volume_inventory(m,volume,caller)
    tet_edges=_tet_edge_set(mesh)
    point_nodes=Dict{Int,Int32}()
    node_points=Dict{Int32,Int}()
    curve_entries=Dict{Int,Vector{Tuple{Float64,Int}}}()
    line_cells=NTuple{2,Int32}[]
    line_entities=Int32[]
    claimed_edges=Set{NTuple{2,Int32}}()
    for curve in curve_tags
        entries=_model_projection_tet_curve_nodes(
            m,mesh,curve,tet_edges,caller)
        curve_entries[curve]=entries
        start_point,stop_point=m.curves[curve]
        _model_projection_point!(
            point_nodes,node_points,start_point,entries[1][2],caller)
        _model_projection_point!(
            point_nodes,node_points,stop_point,entries[end][2],caller)
        for index in 1:(length(entries)-1)
            first_node=Int32(entries[index][2])
            second_node=Int32(entries[index+1][2])
            edge=first_node<second_node ? (first_node,second_node) :
                                          (second_node,first_node)
            edge in claimed_edges && throw(ArgumentError(
                "$caller: mesh edge $edge belongs to multiple model curves"))
            push!(claimed_edges,edge)
            push!(line_cells,(first_node,second_node))
            push!(line_entities,Int32(curve))
        end
    end
    for point in point_tags
        haskey(point_nodes,point) && continue
        node=_model_projection_embedded_point_node(m,mesh,point,1e-12,caller)
        _model_projection_point!(point_nodes,node_points,point,node,caller)
    end

    claimed_faces=Set{NTuple{3,Int32}}()
    surface_cells=NTuple{3,Int32}[]
    surface_entities=Int32[]
    surface_nodes=Dict{Int,Set{Int32}}()
    for surface in surface_tags
        faces=_model_projection_volume_surface_faces!(
            claimed_faces,m,mesh,surface,caller)
        nodes=get!(Set{Int32},surface_nodes,surface)
        for face in faces
            push!(surface_cells,face)
            push!(surface_entities,Int32(surface))
            union!(nodes,face)
        end
    end

    node_entities=fill((3,Int32(volume)),nnodes(mesh))
    for (point,node) in point_nodes
        node_entities[node]=(0,Int32(point))
    end
    for curve in curve_tags
        for (_,node) in curve_entries[curve]
            owner=node_entities[node]
            if owner[1]==3
                node_entities[node]=(1,Int32(curve))
            elseif owner[1]==0 || owner==(1,Int32(curve))
                continue
            else
                throw(ArgumentError(
                    "$caller: mesh node $node belongs to Curves[$(owner[2])] and [$curve]"))
            end
        end
    end
    for (cell,surface) in zip(surface_cells,surface_entities),node in cell
        owner=node_entities[node]
        if owner[1]==3
            node_entities[node]=(2,surface)
        elseif owner[1] in (0,1) || owner==(2,surface)
            continue
        else
            throw(ArgumentError(
                "$caller: mesh node $node belongs to Surfaces[$(owner[2])] and [$surface]"))
        end
    end

    total_elements=0
    for count in (length(point_tags),length(line_cells),
                  length(surface_cells),ntets(mesh))
        total_elements=try Base.checked_add(total_elements,count) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: projected element count overflows Int"))
        end
    end
    total_elements<=typemax(Int32) || throw(ArgumentError(
        "$caller: projected element count exceeds the Int32 MSH limit"))

    entities=Dict{Tuple{Int,Int},MixedEntity}()
    projected_groups=Set{Tuple{Int,Int}}()
    point_physical=Vector{Int32}(undef,length(point_tags))
    for (index,point) in pairs(point_tags)
        tags=_model_projection_physical_tags(m,0,point)
        point_physical[index]=_model_projection_legacy_tag(tags)
        union!(projected_groups,((0,Int(tag)) for tag in tags))
        node=point_nodes[point]
        coordinate=(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
        entities[(0,point)]=MixedEntity(
            0,point,coordinate;physical_tags=tags)
    end

    line_matrix=Matrix{Int32}(undef,2,length(line_cells))
    line_physical=Vector{Int32}(undef,length(line_cells))
    curve_physical=Dict{Int,Vector{Int32}}()
    for curve in curve_tags
        tags=_model_projection_physical_tags(m,1,curve)
        curve_physical[curve]=tags
        union!(projected_groups,((1,Int(tag)) for tag in tags))
        nodes=Int32[entry[2] for entry in curve_entries[curve]]
        start_point,stop_point=m.curves[curve]
        entities[(1,curve)]=MixedEntity(
            1,curve,_model_projection_bbox(mesh,nodes,caller);
            physical_tags=tags,
            boundaries=Int32[Int32(start_point),-Int32(stop_point)])
    end
    for (cell,(first_node,second_node)) in pairs(line_cells)
        line_matrix[1,cell]=first_node;line_matrix[2,cell]=second_node
        line_physical[cell]=_model_projection_legacy_tag(
            curve_physical[Int(line_entities[cell])])
    end

    surface_matrix=Matrix{Int32}(undef,3,length(surface_cells))
    surface_physical=Vector{Int32}(undef,length(surface_cells))
    surface_memberships=Dict{Int,Vector{Int32}}()
    for surface in surface_tags
        tags=_model_projection_physical_tags(m,2,surface)
        surface_memberships[surface]=tags
        union!(projected_groups,((2,Int(tag)) for tag in tags))
        nodes=sort!(collect(surface_nodes[surface]))
        entities[(2,surface)]=MixedEntity(
            2,surface,_model_projection_bbox(mesh,nodes,caller);
            physical_tags=tags,
            boundaries=_model_projection_surface_boundaries(m,surface))
    end
    for (cell,face) in pairs(surface_cells)
        surface_matrix[:,cell].=face
        surface_physical[cell]=_model_projection_legacy_tag(
            surface_memberships[Int(surface_entities[cell])])
    end

    volume_tags=_model_projection_physical_tags(m,3,volume)
    union!(projected_groups,((3,Int(tag)) for tag in volume_tags))
    entities[(3,volume)]=MixedEntity(
        3,volume,_model_projection_bbox(mesh,axes(mesh.coords,2),caller);
        physical_tags=volume_tags)
    tet_physical=fill(_model_projection_legacy_tag(volume_tags),ntets(mesh))

    blocks=ElementBlock[]
    block_entities=Vector{Int32}[]
    if !isempty(point_tags)
        point_matrix=reshape(Int32[point_nodes[tag] for tag in point_tags],1,:)
        push!(blocks,ElementBlock(15,point_matrix,point_physical))
        push!(block_entities,Int32.(point_tags))
    end
    if !isempty(line_cells)
        push!(blocks,ElementBlock(1,line_matrix,line_physical))
        push!(block_entities,line_entities)
    end
    if !isempty(surface_cells)
        push!(blocks,ElementBlock(2,surface_matrix,surface_physical))
        push!(block_entities,surface_entities)
    end
    push!(blocks,ElementBlock(4,mesh.tets,tet_physical))
    push!(block_entities,fill(Int32(volume),ntets(mesh)))

    external_node_tags=UInt64.(1:nnodes(mesh))
    external_element_tags=_model_projection_external_elements(blocks)
    node_parametric=Union{Nothing,Vector{Float64}}[
        nothing for _ in 1:nnodes(mesh)]
    data=MixedEntityData(
        entities;node_entities=node_entities,node_parametric=node_parametric,
        external_node_tags=external_node_tags,block_entities=block_entities,
        external_element_tags=external_element_tags)
    names=Dict{Tuple{Int,Int},String}()
    for key in sort!(collect(projected_groups))
        haskey(m.physical_names,key) && (names[key]=m.physical_names[key])
    end
    output=MixedMesh(
        mesh.coords,blocks;physical_names=names,entity_data=data,
        elementary_entities=block_entities)
    output_diagnostic=validate(output)
    output_diagnostic.ok || throw(ErrorException(
        "$caller: invalid projected mixed mesh — " *
        join(output_diagnostic.messages,"; ")))
    return output
end

"""
    model_to_mixed(model, mesh, entity_dim, entity_tag) -> MixedMesh

Project a validated native surface (`entity_dim=2`) or volume (`entity_dim=3`)
simplex mesh into classified MSH2/MSH4 cells. The three-argument method remains
the surface convenience form. Volume projection first certifies that the tetrahedron
boundary fills the selected native solid, then emits Point/Line/Surface-In-Volume
cells and their entity hierarchy. Gmsh 4.15.2 has no serialized volume-embedding
relation; MSH4 classification and MSH2 cell ownership remain available.
"""
function model_to_mixed(
    m::GeoModel,mesh::Mesh,entity_dim::Integer,entity_tag::Integer)
    caller="model_to_mixed"
    dim=_dimension(entity_dim,caller)
    dim in (2,3) || throw(ArgumentError(
        "$caller: classified projection supports only surface or volume entities"))
    tag=_tag(entity_tag,caller,dim)
    return dim==2 ? model_to_mixed(m,mesh,tag) :
                    _model_volume_to_mixed(m,mesh,tag)
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

function _mesh_model_surface_once(m::GeoModel,t::Int,forced,min_angle_deg,
                                  caller::AbstractString)
    xs,ys,segs,index,embedded,internal=
        _surface_pslg(m,t,forced,caller)
    T=constrained_delaunay(xs,ys,segs; internal_segments=internal)
    hmin=minimum(m.point_size[pid] for pid in keys(index))
    sizefn=(x,y)->hmin
    interior=refine!(T; min_angle_deg=min_angle_deg, size=sizefn)
    mesh=to_mesh(T; interior=interior)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntris(mesh)>0 || throw(ErrorException(
        "$caller: Surface[$t] produced no triangles"))
    return mesh,embedded
end

function _validate_surface_embeddings(m::GeoModel,mesh::Mesh,embedded,
                                      caller::AbstractString)
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
    return nothing
end

"""
    mesh_model_surface(model, tag; min_angle_deg=25.0,
                       max_periodic_passes=8) -> Mesh

Mesh a native planar surface in `z=0`, including holes and embedded points or
curves. Point characteristic lengths drive refinement. Stored straight-curve
periodic relations synchronize their boundary subdivisions through bounded
remeshing passes before slave nodes are certified and snapped. The returned
triangle mesh is validated before it is returned. Relations meeting at a corner
must produce the same exact snapped coordinate.
"""
function mesh_model_surface(m::GeoModel,tag::Integer;min_angle_deg::Real=25.0,
                            max_periodic_passes=8)
    caller="mesh_model_surface"
    t=_tag(tag,caller,2)
    haskey(m.surfaces,t) || throw(ArgumentError(
        "$caller: unknown Surface[$t]"))
    max_periodic_passes isa Bool && throw(ArgumentError(
        "$caller: max_periodic_passes must not be Bool"))
    max_periodic_passes isa Integer || throw(ArgumentError(
        "$caller: max_periodic_passes must be an integer"))
    1<=max_periodic_passes<=64 || throw(ArgumentError(
        "$caller: max_periodic_passes must be in 1:64"))
    npasses=Int(max_periodic_passes)
    constraints=_surface_periodic_constraints(m,t,caller)
    forced=Dict{Int,Vector{Float64}}()
    mesh=nothing;embedded=NTuple{2,Int}[]
    for pass in 1:npasses
        mesh,embedded=_mesh_model_surface_once(
            m,t,forced,min_angle_deg,caller)
        isempty(constraints) && break
        changed=_synchronize_periodic_parameters!(
            forced,m,mesh,constraints)
        if changed
            pass<npasses || throw(ErrorException(
                "$caller: periodic boundary synchronization did not converge " *
                "within $npasses passes"))
            continue
        end
        mesh=_snap_surface_periodic(m,mesh,constraints,caller)
        break
    end
    mesh isa Mesh || throw(ErrorException(
        "$caller: internal surface meshing pass produced no mesh"))
    _validate_surface_embeddings(m,mesh,embedded,caller)
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "$caller: invalid final mesh — "*join(diagnostic.messages,"; ")))
    return mesh
end

function _volume_surface(m::GeoModel,t::Int,caller::AbstractString="mesh_model_volume")
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
        A=_volume_surface(m,spec.a,caller); B=_volume_surface(m,spec.b,caller)
        return mesh_boolean(A,B,spec.op)
    end
    throw(ArgumentError("$caller: Volume[$t] has no native solid encoding"))
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
