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
using ..Geometry: box_surface
using ..Mesh3D: tetrahedralize

export GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
export add_box!, add_physical_group!, set_physical_name!
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
    next_tag::Vector{Int}
end

GeoModel() = GeoModel(Dict{Int,NTuple{3,Float64}}(), Dict{Int,Float64}(),
                      Dict{Int,NTuple{2,Int}}(), Dict{Int,Vector{Int}}(),
                      Dict{Int,Vector{Int}}(), Dict{Int,Vector{Int}}(),
                      Dict{Tuple{Int,Int},Vector{Int}}(),
                      Dict{Tuple{Int,Int},String}(),
                      Dict{Int,NTuple{6,Float64}}(),
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
    T=constrained_delaunay(xs,ys,segs)
    hmin=minimum(m.point_size[pid] for pid in _loop_points(m,loops[1]))
    sizefn=(x,y)->hmin
    interior=refine!(T; min_angle_deg=min_angle_deg, size=sizefn)
    mesh=to_mesh(T; interior=interior)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntris(mesh)>0 || throw(ErrorException("$caller: Surface[$t] produced no triangles"))
    return mesh
end

function mesh_model_volume(m::GeoModel, tag::Integer)
    caller="mesh_model_volume"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    extents=get(m.box_extents,t,nothing)
    extents===nothing && throw(ArgumentError("$caller: Volume[$t] has no native box encoding"))
    x0,y0,z0,dx,dy,dz=extents
    surface=box_surface(x0,x0+dx,y0,y0+dy,z0,z0+dz)
    mesh=tetrahedralize(surface)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntets(mesh)>0 || throw(ErrorException("$caller: Volume[$t] produced no tetrahedra"))
    return mesh
end

end # module
