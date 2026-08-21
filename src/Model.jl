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
using ..Mesh3D: tetrahedralize, mesh_boolean

export GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
export add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
export embed!, dilate_volume!, rotate_volume!
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
    embeds::Dict{Int,Vector{NTuple{2,Int}}}
    next_tag::Vector{Int}
end

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
                      Dict{Int,Vector{NTuple{2,Int}}}(),
                      Int[0,0,0,0])

function _tag(value, caller, dim::Int)
    value isa Integer || throw(ArgumentError("$caller: tag must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: tag must not be Bool"))
    t=Int(value)
    t<0 && throw(ArgumentError("$caller: tag must be non-negative"))
    t>typemax(Int32) && throw(ArgumentError("$caller: tag exceeds Int32"))
    return t
end

function _alloc_tag!(m::GeoModel, dim::Int, requested::Int, caller)
    if requested==0
        m.next_tag[dim+1]+=1
        return m.next_tag[dim+1]
    end
    m.next_tag[dim+1]=max(m.next_tag[dim+1], requested)
    return requested
end

function _finite3(x,y,z,caller)
    p=try (Float64(x),Float64(y),Float64(z)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: coordinates must be Float64-representable"))
    end
    all(isfinite,p) || throw(ArgumentError("$caller: coordinates must be finite"))
    return p
end

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

function add_curve_loop!(m::GeoModel, curves; tag::Integer=0)
    caller="add_curve_loop!"
    ids=Int[_tag(c,caller,1) for c in curves]
    length(ids)>=3 || throw(ArgumentError("$caller: a loop needs at least three curves"))
    for id in ids
        haskey(m.curves,abs(id)) || throw(ArgumentError("$caller: unknown Curve[$(abs(id))]"))
    end
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller)
    haskey(m.loops,t) && throw(ArgumentError("$caller: Loop[$t] already exists"))
    m.loops[t]=ids
    return t
end

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

function add_cylinder!(m::GeoModel, x, y, z, dx, dy, dz, radius; tag::Integer=0)
    caller="add_cylinder!"
    c=_finite3(x,y,z,caller); a=_finite3(dx,dy,dz,caller)
    h=hypot(a...)
    h>0 || throw(ArgumentError("$caller: axis must have positive length"))
    r=try Float64(radius) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: radius must be Float64-representable"))
    end
    r>0 || throw(ArgumentError("$caller: radius must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.cylinders[t]=(center=c, axis=a, radius=r, height=h)
    return t
end

function add_sphere!(m::GeoModel, x, y, z, radius; tag::Integer=0)
    caller="add_sphere!"
    c=_finite3(x,y,z,caller)
    r=try Float64(radius) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: radius must be Float64-representable"))
    end
    r>0 || throw(ArgumentError("$caller: radius must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.spheres[t]=(center=c, radius=r)
    return t
end

function add_cone!(m::GeoModel, x, y, z, dx, dy, dz, r1, r2; tag::Integer=0)
    caller="add_cone!"
    c=_finite3(x,y,z,caller); a=_finite3(dx,dy,dz,caller)
    h=hypot(a...)
    h>0 || throw(ArgumentError("$caller: axis must have positive length"))
    ra=try Float64(r1) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: r1 must be Float64-representable"))
    end
    rb=try Float64(r2) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: r2 must be Float64-representable"))
    end
    (ra>=0 && rb>=0 && (ra>0 || rb>0)) || throw(ArgumentError(
        "$caller: radii must be non-negative with at least one positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller)
    haskey(m.volumes,t) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.cones[t]=(center=c, axis=a, r1=ra, r2=rb, height=h)
    return t
end

function embed!(m::GeoModel, dim, tags, target_dim, target_tag)
    caller="embed!"
    d=Int(dim); td=Int(target_dim); tt=_tag(target_tag,caller,td)
    d==0 && td==2 || throw(ArgumentError(
        "$caller: only Point In Surface is implemented (got dim=$d in dim=$td)"))
    haskey(m.surfaces,tt) || throw(ArgumentError("$caller: unknown Surface[$tt]"))
    list=get!(Vector{NTuple{2,Int}}, m.embeds, tt)
    for raw in tags
        pt=_tag(raw,caller,0)
        haskey(m.points,pt) || throw(ArgumentError("$caller: unknown Point[$pt]"))
        (0,pt) in list && throw(ArgumentError("$caller: Point[$pt] already embedded in Surface[$tt]"))
        push!(list, (0,pt))
    end
    return tt
end

function _quarter_turns(angle, caller)
    turns=angle/(π/2)
    k=round(Int,turns)
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

function dilate_volume!(m::GeoModel, tag, center, scale)
    caller="dilate_volume!"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    c=_finite3(center[1],center[2],center[3],caller)
    s=try Float64(scale) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: scale must be Float64-representable"))
    end
    (isfinite(s) && s>0) || throw(ArgumentError("$caller: scale must be positive"))
    if haskey(m.box_extents,t)
        x0,y0,z0,dx,dy,dz=m.box_extents[t]
        p0=_dilate_point((x0,y0,z0),c,s)
        m.box_extents[t]=(p0[1],p0[2],p0[3],dx*s,dy*s,dz*s)
    elseif haskey(m.cylinders,t)
        cyl=m.cylinders[t]
        m.cylinders[t]=(center=_dilate_point(cyl.center,c,s), axis=cyl.axis,
                        radius=cyl.radius*s, height=cyl.height*s)
    elseif haskey(m.spheres,t)
        sph=m.spheres[t]
        m.spheres[t]=(center=_dilate_point(sph.center,c,s), radius=sph.radius*s)
    elseif haskey(m.cones,t)
        cone=m.cones[t]
        m.cones[t]=(center=_dilate_point(cone.center,c,s), axis=cone.axis,
                    r1=cone.r1*s, r2=cone.r2*s, height=cone.height*s)
    else
        throw(ArgumentError("$caller: Volume[$t] has no dilatable native encoding"))
    end
    return t
end

function rotate_volume!(m::GeoModel, tag, axis, origin, angle)
    caller="rotate_volume!"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    ax=_finite3(axis[1],axis[2],axis[3],caller)
    o=_finite3(origin[1],origin[2],origin[3],caller)
    θ=try Float64(angle) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: angle must be Float64-representable"))
    end
    isfinite(θ) || throw(ArgumentError("$caller: angle must be finite"))
    k=_quarter_turns(θ,caller)
    kind,sgn=_axis_kind(ax,caller)
    rot(p)=_rot90(p,o,kind,sgn,k)
    if haskey(m.box_extents,t)
        x0,y0,z0,dx,dy,dz=m.box_extents[t]
        corners=[rot((x0+ix*dx,y0+iy*dy,z0+iz*dz)) for ix in (0,1), iy in (0,1), iz in (0,1)]
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
        m.cylinders[t]=(center=c2, axis=(e2[1]-c2[1],e2[2]-c2[2],e2[3]-c2[3]),
                        radius=cyl.radius, height=cyl.height)
    elseif haskey(m.spheres,t)
        sph=m.spheres[t]
        m.spheres[t]=(center=rot(sph.center), radius=sph.radius)
    elseif haskey(m.cones,t)
        cone=m.cones[t]
        endp=(cone.center[1]+cone.axis[1],cone.center[2]+cone.axis[2],cone.center[3]+cone.axis[3])
        c2=rot(cone.center); e2=rot(endp)
        m.cones[t]=(center=c2, axis=(e2[1]-c2[1],e2[2]-c2[2],e2[3]-c2[3]),
                    r1=cone.r1, r2=cone.r2, height=cone.height)
    else
        throw(ArgumentError("$caller: Volume[$t] has no rotatable native encoding"))
    end
    return t
end

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

function add_physical_group!(m::GeoModel, dim::Integer, tags; tag::Integer=0, name::AbstractString="")
    caller="add_physical_group!"
    d=_tag(dim,caller,0)
    0<=d<=3 || throw(ArgumentError("$caller: dimension must be in 0:3"))
    ents=Int[_tag(t,caller,d) for t in tags]
    isempty(ents) && throw(ArgumentError("$caller: physical group needs at least one entity"))
    pt=_alloc_tag!(m,d,_tag(tag,caller,d),caller)
    haskey(m.physical,(d,pt)) && throw(ArgumentError("$caller: Physical($d,$pt) already exists"))
    m.physical[(d,pt)]=ents
    isempty(name) || (m.physical_names[(d,pt)]=String(name))
    return pt
end

function set_physical_name!(m::GeoModel, dim::Integer, tag::Integer, name::AbstractString)
    key=(Int(dim),Int(tag))
    haskey(m.physical,key) || throw(ArgumentError("set_physical_name!: unknown Physical$key"))
    m.physical_names[key]=String(name)
    return name
end

function model_entity(m::GeoModel, dim::Integer, tag::Integer)
    d,t=Int(dim),Int(tag)
    d==0 && return get(m.points,t,nothing)
    d==1 && return get(m.curves,t,nothing)
    d==2 && return get(m.surfaces,t,nothing)
    d==3 && return get(m.volumes,t,nothing)
    throw(ArgumentError("model_entity: dimension must be in 0:3"))
end

model_physical_tags(m::GeoModel, dim::Integer, tag::Integer) =
    get(m.physical,(Int(dim),Int(tag)),Int[])

function _loop_points(m::GeoModel, loop_id::Int)
    curves=m.loops[loop_id]
    pts=Int[]
    for signed in curves
        a,b=m.curves[abs(signed)]
        signed>0 ? push!(pts,a) : push!(pts,b)
    end
    return pts
end

function mesh_model_surface(m::GeoModel, tag::Integer; min_angle_deg::Real=25.0)
    caller="mesh_model_surface"
    t=_tag(tag,caller,2)
    haskey(m.surfaces,t) || throw(ArgumentError("$caller: unknown Surface[$t]"))
    loops=m.surfaces[t]
    xs=Float64[]; ys=Float64[]; segs=Tuple{Int,Int}[]
    for (li,loop_id) in enumerate(loops)
        ids=_loop_points(m,loop_id)
        start=length(xs)+1
        for pid in ids
            p=m.points[pid]
            abs(p[3])<=1e-12 || throw(ArgumentError(
                "$caller: Surface[$t] is not planar in z=0 (got z=$(p[3]))"))
            push!(xs,p[1]); push!(ys,p[2])
        end
        stop=length(xs)
        for i in start:stop
            j=i==stop ? start : i+1
            push!(segs,(i,j))
        end
    end
    embedded=get(m.embeds,t,NTuple{2,Int}[])
    for (edim,etag) in embedded
        edim==0 || throw(ArgumentError(
            "$caller: only Point In Surface embeddings are implemented"))
        haskey(m.points,etag) || throw(ArgumentError("$caller: unknown embedded Point[$etag]"))
        p=m.points[etag]
        abs(p[3])<=1e-12 || throw(ArgumentError(
            "$caller: embedded Point[$etag] is not planar in z=0 (got z=$(p[3]))"))
        push!(xs,p[1]); push!(ys,p[2])
    end
    T=constrained_delaunay(xs,ys,segs)
    hmin=minimum(m.point_size[pid] for pid in _loop_points(m,loops[1]))
    sizefn=(x,y)->hmin
    interior=refine!(T; min_angle_deg=min_angle_deg, size=sizefn)
    mesh=to_mesh(T; interior=interior)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntris(mesh)>0 || throw(ErrorException("$caller: Surface[$t] produced no triangles"))
    for (edim,etag) in embedded
        edim==0 || continue
        p=m.points[etag]
        found=false
        @inbounds for i in 1:nnodes(mesh)
            if hypot(mesh.coords[1,i]-p[1],mesh.coords[2,i]-p[2],mesh.coords[3,i]-p[3])<=1e-12
                found=true; break
            end
        end
        found || throw(ErrorException(
            "$caller: embedded Point[$etag] at $p is not a mesh node"))
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

function mesh_model_volume(m::GeoModel, tag::Integer)
    caller="mesh_model_volume"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    surface=_volume_surface(m,t)
    mesh=tetrahedralize(surface)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntets(mesh)>0 || throw(ErrorException("$caller: Volume[$t] produced no tetrahedra"))
    return mesh
end

end # module
