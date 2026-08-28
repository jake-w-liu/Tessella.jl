"""
    Model

A native geometry/entity kernel: tagged points, curves, curve loops, surfaces,
surface loops, and volumes with physical groups and classified surface/volume
mixed-mesh projection.
Meshing dispatches to Tessella's certified simplex and transfinite kernels. This
is not OpenCASCADE; unsupported CAD statements remain explicit blockers.
"""
module Model

using ..MeshTypes: Mesh, validate, nnodes, nsegs, ntris, ntets, boundary_faces,
                   triangle_area, tet_signed_volume, tet_volume
using ..Elements: ElementBlock, MixedEntity, MixedEntityData,
                  MixedPeriodicLink, MixedMesh
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..SizeField: AbstractSizeField, ConstantSize, PostViewField, field_value,
                   size_at
using ..Geometry: box_surface, cylinder_surface, sphere_surface, cone_surface
using ..Mesh3D: tetrahedralize, mesh_boolean, recover_segment3, recover_triangle3
using ..Mesh3D: mesh_covers_segment3, mesh_covers_triangle3,
                _tet_edge_set, _mesh_covering_faces3, _certify_surface_fill
using ..Periodic: periodic_identify_affine
using ..Transform: _affine_coordinate, _transform_homogeneous
using ..Predicates: orient2, orient3

export GeoModel, add_point!, set_point_mesh_size!, add_line!, add_curve_loop!, add_plane_surface!
export add_surface_loop!, add_volume!
export add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
export embed!, translate_volume!, dilate_volume!, rotate_volume!
export ModelPeriodicConstraint, set_periodic!, model_periodic_constraints,
       model_periodic_nodes, model_to_mixed
export add_physical_group!, set_physical_name!, remove_physical_groups!
export remove_physical_name!, model_physical_groups
export model_physical_groups_entities, model_entities_for_physical_group
export model_entities_for_physical_name, model_physical_groups_for_entity
export model_physical_name
export model_entities, model_dimension, model_boundary, model_adjacencies
export mesh_model_surface, mesh_model_volume, model_entity, model_physical_tags

"""
Owned affine relation between two native entities. `affine` maps the master
entity to the slave entity in Gmsh row-major 4×4 order. For dimension 1,
`reversed` records whether the master start maps to the slave end; it is false
for dimension 2.
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
    surface_loops::Dict{Int,Vector{Int}}
    volumes::Dict{Int,Vector{Int}}
    physical::Dict{Tuple{Int,Int},Vector{Int}}
    physical_names::Dict{Tuple{Int,Int},String}
    physical_tag_max::Int
    box_extents::Dict{Int,NTuple{6,Float64}}
    cylinders::Dict{Int,NamedTuple{(:center,:axis,:radius,:height),
                                   Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64}}}
    spheres::Dict{Int,NamedTuple{(:center,:radius),Tuple{NTuple{3,Float64},Float64}}}
    cones::Dict{Int,NamedTuple{(:center,:axis,:r1,:r2,:height),
                               Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64,Float64}}}
    booleans::Dict{Int,NamedTuple{(:op,:a,:b),Tuple{Symbol,Int,Int}}}
    boolean_operands::Dict{Int,Tuple{Mesh,Mesh}}
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
                      Dict{Int,Vector{Int}}(),
                      Dict{Tuple{Int,Int},Vector{Int}}(),
                      Dict{Tuple{Int,Int},String}(),
                      0,
                      Dict{Int,NTuple{6,Float64}}(),
                      Dict{Int,NamedTuple{(:center,:axis,:radius,:height),
                           Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64}}}(),
                      Dict{Int,NamedTuple{(:center,:radius),Tuple{NTuple{3,Float64},Float64}}}(),
                      Dict{Int,NamedTuple{(:center,:axis,:r1,:r2,:height),
                           Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64,Float64}}}(),
                      Dict{Int,NamedTuple{(:op,:a,:b),Tuple{Symbol,Int,Int}}}(),
                      Dict{Int,Tuple{Mesh,Mesh}}(),
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

function _signed_surface_tag(value, caller)
    value isa Integer || throw(ArgumentError("$caller: surface tag must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: surface tag must not be Bool"))
    value==0 && throw(ArgumentError("$caller: surface tag must be nonzero"))
    (-typemax(Int32)<=value<=typemax(Int32)) || throw(ArgumentError(
        "$caller: surface tag magnitude exceeds Int32"))
    return Int(value)
end

function _dimension(value, caller)
    value isa Integer || throw(ArgumentError("$caller: dimension must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: dimension must not be Bool"))
    (0<=value<=3) || throw(ArgumentError("$caller: dimension must be in 0:3"))
    return Int(value)
end

function _query_dimension(value,caller)
    value isa Integer || throw(ArgumentError(
        "$caller: dimension must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: dimension must not be Bool"))
    value==-1 || 0<=value<=3 || throw(ArgumentError(
        "$caller: dimension must be -1 or in 0:3"))
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

function _alloc_physical_tag(m::GeoModel, requested::Int, caller)
    requested!=0 && return requested
    m.physical_tag_max<typemax(Int32) || throw(ArgumentError(
        "$caller: no automatic physical tags remain"))
    return m.physical_tag_max+1
end

function _alloc_surface_loop_tag(m::GeoModel,requested::Int,caller)
    requested!=0 && return requested
    current=isempty(m.surface_loops) ? 0 : maximum(keys(m.surface_loops))
    current<typemax(Int32) || throw(ArgumentError(
        "$caller: no automatic Surface Loop tags remain"))
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
    set_point_mesh_size!(model, points, mesh_size)

Set one finite, positive mesh-size constraint on existing Point tags. The update
is atomic: every tag and the size are validated before the model is changed.
"""
function set_point_mesh_size!(m::GeoModel,points,mesh_size)
    caller="set_point_mesh_size!"
    (points isa AbstractVector || points isa Tuple) || throw(ArgumentError(
        "$caller: points must be a vector or tuple of Point tags"))
    tags=unique(Int[_tag(point,caller,0) for point in points])
    isempty(tags) && throw(ArgumentError(
        "$caller: points must contain at least one Point tag"))
    mesh_size isa Real || throw(ArgumentError(
        "$caller: mesh_size must be a real number"))
    mesh_size isa Bool && throw(ArgumentError(
        "$caller: mesh_size must not be Bool"))
    size=try
        Float64(mesh_size)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: mesh_size must be Float64-representable"))
    end
    (isfinite(size) && size>0) || throw(ArgumentError(
        "$caller: mesh_size must be finite and positive"))
    for tag in tags
        haskey(m.points,tag) || throw(ArgumentError(
            "$caller: unknown Point[$tag]"))
    end
    for tag in tags
        m.point_size[tag]=size
    end
    return nothing
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

function _periodic_entity_tags(values,dim::Int,caller::AbstractString,
                               name::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: $name entities must be a vector or tuple"))
    return Int[_tag(value,caller,dim) for value in values]
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

include("ModelTopologyQueries.jl")

@inline function _model_periodic_entity_label(dim::Int)
    dim==1 && return "Curve"
    dim==2 && return "Surface"
    return "Entity($dim)"
end

function _model_periodic_surface_points(
    m::GeoModel,surface::Int,caller::AbstractString;
    include_embeddings::Bool)
    points=_model_points_of(m,[(2,surface)],caller)
    if include_embeddings
        embedded_points,embedded_curves=
            _model_surface_embedding_tags(m,surface,caller)
        append!(points,embedded_points)
        for curve in embedded_curves
            append!(points,m.curves[curve])
        end
    end
    sort!(unique!(points))
    length(points)>=3 || throw(ArgumentError(
        "$caller: Surface[$surface] needs at least three distinct point tags"))
    coordinates=NTuple{3,Float64}[m.points[point] for point in points]
    _model_surface_projection(coordinates,points,surface,caller)
    return points
end

function _model_periodic_loop_signature(
    m::GeoModel,loop::Int,point_map::Dict{Int,Int})
    points=_loop_points(m,loop)
    edges=NTuple{2,Int}[]
    for index in eachindex(points)
        first_point=point_map[points[index]]
        second_point=point_map[points[mod1(index+1,length(points))]]
        push!(edges,first_point<second_point ?
                    (first_point,second_point) : (second_point,first_point))
    end
    sort!(edges)
    return Tuple(edges)
end

function _model_periodic_surface_topology(
    m::GeoModel,slave::Int,master::Int,point_map::Dict{Int,Int},
    caller::AbstractString;include_embeddings::Bool)
    slave_identity=Dict(point=>point for point in
        _model_periodic_surface_points(
            m,slave,caller;include_embeddings=include_embeddings))
    slave_loops=m.surfaces[slave];master_loops=m.surfaces[master]
    length(slave_loops)==length(master_loops) || throw(ArgumentError(
        "$caller: Surface[$slave] and Surface[$master] loop counts differ"))
    _model_periodic_loop_signature(m,first(master_loops),point_map)==
        _model_periodic_loop_signature(m,first(slave_loops),slave_identity) ||
        throw(ArgumentError(
            "$caller: affine map does not preserve the outer boundary of " *
            "Surface[$master] -> Surface[$slave]"))
    master_holes=sort!([_model_periodic_loop_signature(m,loop,point_map)
                        for loop in Iterators.drop(master_loops,1)])
    slave_holes=sort!([_model_periodic_loop_signature(m,loop,slave_identity)
                       for loop in Iterators.drop(slave_loops,1)])
    master_holes==slave_holes || throw(ArgumentError(
        "$caller: affine map does not preserve the holes of " *
        "Surface[$master] -> Surface[$slave]"))
    include_embeddings || return nothing

    slave_points,slave_curves=
        _model_surface_embedding_tags(m,slave,caller)
    master_points,master_curves=
        _model_surface_embedding_tags(m,master,caller)
    sort!(Int[point_map[point] for point in master_points])==
        sort!(copy(slave_points)) || throw(ArgumentError(
            "$caller: periodic surfaces have different embedded points"))
    function curve_signatures(curves,map)
        signatures=NTuple{2,Int}[]
        for curve in curves
            first_point,second_point=m.curves[curve]
            first_mapped=map[first_point];second_mapped=map[second_point]
            push!(signatures,first_mapped<second_mapped ?
                  (first_mapped,second_mapped) :
                  (second_mapped,first_mapped))
        end
        sort!(signatures)
        return signatures
    end
    curve_signatures(master_curves,point_map)==
        curve_signatures(slave_curves,slave_identity) || throw(ArgumentError(
            "$caller: periodic surfaces have different embedded curves"))
    return nothing
end

function _model_periodic_surface_point_map(
    m::GeoModel,slave::Int,master::Int,affine,atol::Float64,
    caller::AbstractString;include_embeddings::Bool)
    slave_points=_model_periodic_surface_points(
        m,slave,caller;include_embeddings=include_embeddings)
    master_points=_model_periodic_surface_points(
        m,master,caller;include_embeddings=include_embeddings)
    isempty(intersect(Set(slave_points),Set(master_points))) ||
        throw(ArgumentError(
            "$caller: Surface[$slave] and Surface[$master] must have " *
            "disjoint point tags"))
    length(slave_points)==length(master_points) || throw(ArgumentError(
        "$caller: Surface[$slave] and Surface[$master] point counts differ"))
    coefficients,translation,_=_transform_homogeneous(
        affine,caller;name="affine transform")
    available=Set(slave_points)
    point_map=Dict{Int,Int}()
    for (index,master_point) in pairs(master_points)
        expected=_model_affine_point(
            coefficients,translation,m.points[master_point],caller,index)
        matches=Int[slave_point for slave_point in available
                    if _model_point_distance(m.points[slave_point],expected)<=atol]
        length(matches)==1 || throw(ArgumentError(
            "$caller: affine image of Point[$master_point] on Surface[$master] " *
            "matches $(length(matches)) points on Surface[$slave]; expected one"))
        slave_point=only(matches)
        point_map[master_point]=slave_point
        delete!(available,slave_point)
    end
    isempty(available) || throw(ErrorException(
        "$caller: internal periodic surface point matching was incomplete"))
    _model_periodic_surface_topology(
        m,slave,master,point_map,caller;
        include_embeddings=include_embeddings)
    return point_map
end

function _model_periodic_dependency_parents(
    constraints,caller::AbstractString)
    parents=Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for constraint in constraints
        slave=(constraint.dim,Int(constraint.slave_entity))
        master=(constraint.dim,Int(constraint.master_entity))
        haskey(parents,slave) && throw(ArgumentError(
            "$caller: $(_model_periodic_entity_label(slave[1]))[$(slave[2])] " *
            "has more than one periodic master"))
        parents[slave]=master
    end

    state=Dict{Tuple{Int,Int},UInt8}()
    for start in sort!(collect(keys(parents)))
        get(state,start,0x00)==0x00 || continue
        path=Tuple{Int,Int}[]
        positions=Dict{Tuple{Int,Int},Int}()
        current=start
        while haskey(parents,current) && get(state,current,0x00)==0x00
            state[current]=0x01
            push!(path,current)
            positions[current]=length(path)
            current=parents[current]
        end
        if get(state,current,0x00)==0x01
            first_cycle=get(positions,current,0)
            first_cycle>0 || throw(ErrorException(
                "$caller: internal periodic dependency traversal failed"))
            cycle=vcat(path[first_cycle:end],current)
            description=join((
                "$(_model_periodic_entity_label(entity[1]))[$(entity[2])]"
                for entity in cycle)," -> ")
            throw(ArgumentError(
                "$caller: cyclic periodic dependency $description"))
        end
        for entity in path
            state[entity]=0x02
        end
    end
    return parents
end

function _model_periodic_constraint_order(constraints,caller::AbstractString)
    parents=_model_periodic_dependency_parents(constraints,caller)
    depths=Dict{Tuple{Int,Int},Int}()
    for start in sort!(collect(keys(parents)))
        haskey(depths,start) && continue
        path=Tuple{Int,Int}[]
        current=start
        while haskey(parents,current) && !haskey(depths,current)
            push!(path,current)
            current=parents[current]
        end
        depth=get(depths,current,-1)
        for entity in Iterators.reverse(path)
            depth+=1
            depths[entity]=depth
        end
    end
    return sort!(collect(constraints);by=constraint->(
        constraint.dim,
        depths[(constraint.dim,Int(constraint.slave_entity))],
        constraint.slave_entity))
end

"""
    set_periodic!(model, dim, slave_entities, master_entities, affine;
                  atol=1e-12)

Persist affine relations between equally sized lists of straight native curves
(`dim=1`) or planar native surfaces (`dim=2`). `affine` maps each master entity
to its slave in Gmsh row-major 4×4 order. An entity may be the master of multiple
relations or both a slave and a master in an acyclic dependency chain; each
slave has exactly one master. Cycles are explicit blockers. Periodic curves must
belong to the same planar surface when meshed. Periodic surfaces require
disjoint, affine-equivalent boundary point and loop topology; embedded topology
is rechecked when an explicit planar-shell volume is meshed. Volume periodicity
remains an explicit blocker. The update is atomic.
"""
function set_periodic!(m::GeoModel,dim,slave_entities,master_entities,affine;
                       atol=1e-12)
    caller="set_periodic!"
    d=_dimension(dim,caller)
    d in (1,2) || throw(ArgumentError(
        "$caller: only straight Curve and planar Surface periodicity " *
        "(dimensions 1 and 2) are implemented"))
    label=_model_periodic_entity_label(d)
    slaves=_periodic_entity_tags(slave_entities,d,caller,"slave")
    masters=_periodic_entity_tags(master_entities,d,caller,"master")
    length(slaves)==length(masters) || throw(ArgumentError(
        "$caller: slave and master entity counts differ"))
    isempty(slaves) && throw(ArgumentError(
        "$caller: need at least one slave/master $label pair"))
    length(unique(slaves))==length(slaves) || throw(ArgumentError(
        "$caller: slave $label tags must be unique"))
    overlap=sort!(Int[slave for slave in slaves
                      if haskey(m.periodic,(d,slave))])
    isempty(overlap) || throw(ArgumentError(
        "$caller: $label[$(first(overlap))] already has a periodic master"))
    tolerance=_model_periodic_tolerance(atol,caller)
    coefficients,translation,row_major=_transform_homogeneous(
        affine,caller;name="affine transform")
    pending=ModelPeriodicConstraint[]
    for (pair_index,(slave,master)) in enumerate(zip(slaves,masters))
        entities=d==1 ? m.curves : m.surfaces
        haskey(entities,slave) || throw(ArgumentError(
            "$caller: unknown slave $label[$slave]"))
        haskey(entities,master) || throw(ArgumentError(
            "$caller: unknown master $label[$master]"))
        slave!=master || throw(ArgumentError(
            "$caller: slave and master $label tags must differ"))
        reversed=false
        if d==1
            _model_curve_length(m,slave,caller)
            _model_curve_length(m,master,caller)
            slave_points=m.curves[slave];master_points=m.curves[master]
            isempty(intersect(Set(slave_points),Set(master_points))) ||
                throw(ArgumentError(
                    "$caller: Curve[$slave] and Curve[$master] must have " *
                    "disjoint endpoints"))
            slave_start=m.points[slave_points[1]]
            slave_stop=m.points[slave_points[2]]
            mapped_start=_model_affine_point(
                coefficients,translation,m.points[master_points[1]],caller,
                2pair_index-1)
            mapped_stop=_model_affine_point(
                coefficients,translation,m.points[master_points[2]],caller,
                2pair_index)
            forward=max(_model_point_distance(mapped_start,slave_start),
                        _model_point_distance(mapped_stop,slave_stop))
            reverse_error=max(_model_point_distance(mapped_start,slave_stop),
                              _model_point_distance(mapped_stop,slave_start))
            mismatch=min(forward,reverse_error)
            mismatch<=tolerance || throw(ArgumentError(
                "$caller: affine map misses slave Curve[$slave] endpoints " *
                "by $mismatch"))
            reversed=reverse_error<forward
        else
            _model_periodic_surface_point_map(
                m,slave,master,row_major,tolerance,caller;
                include_embeddings=false)
        end
        push!(pending,ModelPeriodicConstraint(
            d,Int32(slave),Int32(master),row_major,reversed,tolerance))
    end
    _model_periodic_dependency_parents(
        vcat(model_periodic_constraints(m),pending),caller)
    for constraint in pending
        m.periodic[(constraint.dim,Int(constraint.slave_entity))]=constraint
    end
    return nothing
end

"""
    model_periodic_constraints(model) -> Vector{ModelPeriodicConstraint}

Return the model's immutable periodic constraints in deterministic dimension and
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

function _validate_surface_loop(m::GeoModel,surfaces,caller::AbstractString)
    unsigned=abs.(surfaces)
    length(unique(unsigned))==length(unsigned) || throw(ArgumentError(
        "$caller: a surface loop cannot repeat a Surface tag"))
    curve_owners=Dict{Int,Set{Int}}()
    curve_incidence=Dict{Int,Int}()
    for surface in unsigned
        haskey(m.surfaces,surface) || throw(ArgumentError(
            "$caller: unknown Surface[$surface]"))
        for loop in m.surfaces[surface]
            haskey(m.loops,loop) || throw(ArgumentError(
                "$caller: Surface[$surface] references unknown Loop[$loop]"))
            for signed_curve in m.loops[loop]
                curve=abs(signed_curve)
                haskey(m.curves,curve) || throw(ArgumentError(
                    "$caller: Loop[$loop] references unknown Curve[$curve]"))
                curve_incidence[curve]=get(curve_incidence,curve,0)+1
                push!(get!(Set{Int},curve_owners,curve),surface)
            end
        end
    end
    isempty(curve_incidence) && throw(ArgumentError(
        "$caller: a surface loop must contain boundary curves"))
    for curve in sort!(collect(keys(curve_incidence)))
        count=curve_incidence[curve]
        owners=curve_owners[curve]
        (count==2 && length(owners)==2) || throw(ArgumentError(
            "$caller: Curve[$curve] must occur once on each of two distinct " *
            "surfaces (found $count occurrences on $(length(owners)) surfaces)"))
    end
    adjacency=Dict(surface=>Set{Int}() for surface in unsigned)
    for owners in values(curve_owners)
        first_surface,second_surface=Tuple(owners)
        push!(adjacency[first_surface],second_surface)
        push!(adjacency[second_surface],first_surface)
    end
    visited=Set{Int}();stack=Int[first(unsigned)]
    while !isempty(stack)
        surface=pop!(stack)
        surface in visited && continue
        push!(visited,surface)
        append!(stack,adjacency[surface])
    end
    length(visited)==length(unsigned) || throw(ArgumentError(
        "$caller: surfaces must form one connected closed shell"))
    return nothing
end

"""
    add_surface_loop!(model, surfaces; tag=0) -> tag

Add one connected closed shell of existing planar surfaces. Signed surface tags
are retained for entity-boundary metadata. Every shell curve must be shared by
exactly two distinct surfaces; open, branched, disconnected, or repeated shells
are rejected before the model is changed.
"""
function add_surface_loop!(m::GeoModel,surfaces;tag::Integer=0)
    caller="add_surface_loop!"
    ids=Int[_signed_surface_tag(surface,caller) for surface in surfaces]
    isempty(ids) && throw(ArgumentError("$caller: need at least one Surface"))
    _validate_surface_loop(m,ids,caller)
    requested=_tag(tag,caller,2)
    t=_alloc_surface_loop_tag(m,requested,caller)
    haskey(m.surface_loops,t) && throw(ArgumentError(
        "$caller: Surface Loop[$t] already exists"))
    m.surface_loops[t]=ids
    return t
end

"""
    add_volume!(model, surface_loops; tag=0) -> tag

Add an explicitly modeled volume. The first surface loop is the exterior shell;
later loops are cavities. Shells must exist, be distinct, and have disjoint
surface entities. Meshing and classified projection also certify that the cavity
shells are disjoint and lie inside the exterior. Native primitive and Boolean
volumes use the same volume tag namespace.
"""
function add_volume!(m::GeoModel,surface_loops;tag::Integer=0)
    caller="add_volume!"
    shells=Int[_tag(shell,caller,3) for shell in surface_loops]
    isempty(shells) && throw(ArgumentError(
        "$caller: need an exterior Surface Loop"))
    all(shell->shell>0,shells) || throw(ArgumentError(
        "$caller: Surface Loop tags must be positive"))
    length(unique(shells))==length(shells) || throw(ArgumentError(
        "$caller: a volume cannot repeat a Surface Loop tag"))
    seen_surfaces=Set{Int}()
    for shell in shells
        haskey(m.surface_loops,shell) || throw(ArgumentError(
            "$caller: unknown Surface Loop[$shell]"))
        _validate_surface_loop(m,m.surface_loops[shell],caller)
        for signed_surface in m.surface_loops[shell]
            surface=abs(signed_surface)
            surface in seen_surfaces && throw(ArgumentError(
                "$caller: Surface[$surface] belongs to multiple shells"))
            push!(seen_surfaces,surface)
        end
    end
    requested=_tag(tag,caller,3)
    requested!=0 && haskey(m.volumes,requested) && throw(ArgumentError(
        "$caller: Volume[$requested] already exists"))
    t=_alloc_tag!(m,3,requested,caller)
    haskey(m.volumes,t) && throw(ArgumentError(
        "$caller: Volume[$t] already exists"))
    m.volumes[t]=shells
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
operations are `:union`, `:intersection`, and `:difference` (`a \\ b`). The
result owns operation-time snapshots of both operands, so later operand changes do
not change the Boolean geometry.
"""
function boolean_volumes!(m::GeoModel, op::Symbol, a, b; tag::Integer=0)
    caller="boolean_volumes!"
    op in (:union,:intersection,:difference) || throw(ArgumentError(
        "$caller: op must be :union, :intersection, or :difference"))
    ta=_tag(a,caller,3); tb=_tag(b,caller,3)
    haskey(m.volumes,ta) || throw(ArgumentError("$caller: unknown Volume[$ta]"))
    haskey(m.volumes,tb) || throw(ArgumentError("$caller: unknown Volume[$tb]"))
    requested=_tag(tag,caller,3)
    requested!=0 && haskey(m.volumes,requested) && throw(ArgumentError(
        "$caller: Volume[$requested] already exists"))
    operand_a=_volume_surface(m,ta,caller)
    operand_b=_volume_surface(m,tb,caller)
    t=_alloc_tag!(m,3,requested,caller)
    haskey(m.volumes,t) && throw(ArgumentError(
        "$caller: automatic Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.booleans[t]=(op=op, a=ta, b=tb)
    m.boolean_operands[t]=(operand_a,operand_b)
    return t
end

function _remove_volume_entity!(m::GeoModel,tag::Int)
    haskey(m.volumes,tag) || return false
    delete!(m.volumes,tag)
    for encoding in (m.box_extents,m.cylinders,m.spheres,m.cones,
                     m.booleans,m.boolean_operands)
        delete!(encoding,tag)
    end
    delete!(m.embeds,(3,tag))
    for key in collect(keys(m.physical))
        key[1]==3 || continue
        entities=m.physical[key]
        tag in entities || continue
        filter!(entity->entity!=tag,entities)
        isempty(entities) || continue
        delete!(m.physical,key)
        delete!(m.physical_names,key)
    end
    return true
end

function _has_entity(m::GeoModel, dim::Int, tag::Int)
    dim==0 && return haskey(m.points,tag)
    dim==1 && return haskey(m.curves,tag)
    dim==2 && return haskey(m.surfaces,tag)
    return haskey(m.volumes,tag)
end

function _physical_group_key(dim,tag,caller)
    d=_dimension(dim,caller)
    t=_tag(tag,caller,d)
    t>0 || throw(ArgumentError("$caller: physical tag must be positive"))
    return (d,t)
end

function _physical_name_is_available(m::GeoModel,dimension::Int,
                                     name::String)
    isempty(name) && return false
    return all(d!=dimension || existing!=name
               for ((d,_),existing) in m.physical_names)
end

"""
    add_physical_group!(model, dim, tags; tag=0, name="") -> tag

Group one or more existing entities of dimension `dim` under an independent
physical tag. Automatic physical tags share one global namespace across entity
dimensions, separate from geometry tags. An optional nonempty `name` is recorded
when no group in the same dimension already uses it; the group is still added when
the requested name is unavailable.
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
    pt=_alloc_physical_tag(m,_tag(tag,caller,d),caller)
    haskey(m.physical,(d,pt)) && throw(ArgumentError("$caller: Physical($d,$pt) already exists"))
    m.physical[(d,pt)]=ents
    _physical_name_is_available(m,d,group_name) &&
        (m.physical_names[(d,pt)]=group_name)
    m.physical_tag_max=max(m.physical_tag_max,pt)
    return pt
end

"""
    set_physical_name!(model, dim, tag, name) -> name

Assign `name` to an unnamed physical group. Names are unique within one entity
dimension. As in Gmsh 4.15.2, an unknown group, an empty name, an already named
group, or a name already used in that dimension is a no-op; use
[`remove_physical_name!`](@ref) before assigning a replacement.
"""
function set_physical_name!(m::GeoModel, dim::Integer, tag::Integer, name::AbstractString)
    caller="set_physical_name!"
    key=_physical_group_key(dim,tag,caller)
    haskey(m.physical,key) || return ""
    existing=get(m.physical_names,key,"")
    isempty(existing) || return existing
    group_name=String(name)
    _physical_name_is_available(m,key[1],group_name) || return ""
    m.physical_names[key]=group_name
    return group_name
end

function _physical_dim_tags(dim_tags,caller)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    result=Tuple{Int,Int}[]
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        push!(result,_physical_group_key(pair[1],pair[2],caller))
    end
    return unique!(result)
end

"""
    remove_physical_groups!(model, dim_tags=()) -> Int

Remove the selected `(dimension, physical_tag)` groups and their names, returning
the number removed. An empty selection removes every group. Unknown valid groups
are ignored, all inputs are checked before mutation, and the global automatic-tag
counter remains monotonic.
"""
function remove_physical_groups!(m::GeoModel,dim_tags=())
    caller="remove_physical_groups!"
    selected=_physical_dim_tags(dim_tags,caller)
    targets=isempty(selected) ? collect(keys(m.physical)) : selected
    removed=0
    for key in targets
        haskey(m.physical,key) || continue
        delete!(m.physical,key)
        delete!(m.physical_names,key)
        removed+=1
    end
    return removed
end

"""
    remove_physical_name!(model, name) -> Int

Remove `name` from every Physical group carrying it without removing any group or
geometry. Return the number of names removed; a missing or empty name is a no-op.
"""
function remove_physical_name!(m::GeoModel,name::AbstractString)
    group_name=String(name)
    isempty(group_name) && return 0
    targets=Tuple{Int,Int}[
        key for (key,existing) in m.physical_names if existing==group_name]
    for key in targets
        delete!(m.physical_names,key)
    end
    return length(targets)
end

"""
    model_physical_groups(model, dim=-1) -> Vector{Tuple{Int,Int}}

Return detached Physical `(dimension, tag)` pairs in deterministic order. `dim=-1`
selects every dimension; `dim=0:3` filters the result.
"""
function model_physical_groups(m::GeoModel,dim=-1)
    dimension=_query_dimension(dim,"model_physical_groups")
    groups=Tuple{Int,Int}[
        key for key in keys(m.physical) if dimension==-1 || key[1]==dimension]
    return sort!(groups)
end

"""
    model_entities_for_physical_group(model, dim, tag) -> Vector{Int}

Return detached, sorted entity tags for an existing Physical group.
"""
function model_entities_for_physical_group(m::GeoModel,dim,tag)
    key=_physical_group_key(dim,tag,"model_entities_for_physical_group")
    haskey(m.physical,key) || throw(ArgumentError(
        "model_entities_for_physical_group: Physical$(key) does not exist"))
    return sort!(copy(m.physical[key]))
end

"""
    model_physical_groups_for_entity(model, dim, tag) -> Vector{Int}

Return the sorted Physical tags containing an existing geometry entity.
"""
function model_physical_groups_for_entity(m::GeoModel,dim,tag)
    caller="model_physical_groups_for_entity"
    dimension=_dimension(dim,caller)
    entity_tag=_tag(tag,caller,dimension)
    entity_tag>0 || throw(ArgumentError("$caller: entity tag must be positive"))
    key=(dimension,entity_tag)
    _has_entity(m,key...) || throw(ArgumentError(
        "$caller: entity $(key) does not exist"))
    groups=Int[physical_tag for ((group_dimension,physical_tag),entities) in m.physical
               if group_dimension==key[1] && key[2] in entities]
    return sort!(groups)
end

"""
    model_physical_name(model, dim, tag) -> String

Return the Physical group name, or an empty string when the group is unnamed or
does not exist.
"""
function model_physical_name(m::GeoModel,dim,tag)
    key=_physical_group_key(dim,tag,"model_physical_name")
    return get(m.physical_names,key,"")
end

"""
    model_entities_for_physical_name(model, name) -> Vector{Tuple{Int,Int}}

Return detached, sorted geometry `(dimension, tag)` pairs belonging to every group
named `name`. The same name can identify one group in each dimension. A missing or
empty name is rejected.
"""
function model_entities_for_physical_name(m::GeoModel,name::AbstractString)
    group_name=String(name)
    isempty(group_name) && throw(ArgumentError(
        "model_entities_for_physical_name: Physical name must not be empty"))
    group_keys=Tuple{Int,Int}[
        key for (key,existing) in m.physical_names if existing==group_name]
    isempty(group_keys) && throw(ArgumentError(
        "model_entities_for_physical_name: Physical name $(repr(group_name)) does not exist"))
    entities=Tuple{Int,Int}[]
    for key in group_keys, entity in m.physical[key]
        push!(entities,(key[1],entity))
    end
    return sort!(unique!(entities))
end

"""
    model_physical_groups_entities(model, dim=-1)

Return the same ordered groups as [`model_physical_groups`](@ref), together with a
detached vector of sorted `(dimension, entity_tag)` members for each group.
"""
function model_physical_groups_entities(m::GeoModel,dim=-1)
    groups=model_physical_groups(m,dim)
    entities=Vector{Vector{Tuple{Int,Int}}}(undef,length(groups))
    for (index,(dimension,physical_tag)) in pairs(groups)
        entities[index]=Tuple{Int,Int}[
            (dimension,tag) for tag in sort!(copy(m.physical[(dimension,physical_tag)]))]
    end
    return groups,entities
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

function _add_surface_point!(xs,ys,mesh_sizes,index,m::GeoModel,pid::Int,caller)
    haskey(index, pid) && return index[pid]
    haskey(m.points,pid) || throw(ArgumentError("$caller: unknown Point[$pid]"))
    p=m.points[pid]
    abs(p[3])<=1e-12 || throw(ArgumentError(
        "$caller: Point[$pid] is not planar in z=0 (got z=$(p[3]))"))
    push!(xs,p[1]); push!(ys,p[2])
    push!(mesh_sizes,m.point_size[pid])
    index[pid]=length(xs)
    return index[pid]
end

function _add_surface_curve_point!(xs,ys,mesh_sizes,index,m::GeoModel,point,
                                   mesh_size::Float64,curve::Int,
                                   caller::AbstractString)
    scale=max(1.0,hypot(point[1],point[2]))
    tolerance=128eps(Float64)*scale
    matched_vertex=0
    matched_point=0
    for point_tag in sort!(collect(keys(index)))
        candidate=m.points[point_tag]
        distance=hypot(candidate[1]-point[1],candidate[2]-point[2],
                       candidate[3]-point[3])
        distance<=tolerance || continue
        vertex=index[point_tag]
        if matched_vertex!=0 && matched_vertex!=vertex
            throw(ArgumentError(
                "$caller: Curve[$curve] subdivision coincides with " *
                "Points[$matched_point] and [$point_tag]"))
        end
        matched_vertex=vertex
        matched_point=point_tag
    end
    if matched_vertex!=0
        mesh_sizes[matched_vertex]=min(mesh_sizes[matched_vertex],mesh_size)
        return matched_vertex
    end
    push!(xs,point[1]);push!(ys,point[2])
    push!(mesh_sizes,mesh_size)
    return length(xs)
end

@inline function _surface_curve_mesh_size(m::GeoModel,curve::Int,
                                          parameter::Float64,
                                          caller::AbstractString)
    a,b=m.curves[curve]
    first_size=m.point_size[a];last_size=m.point_size[b]
    mesh_size=muladd(parameter,last_size-first_size,first_size)
    (isfinite(mesh_size) && mesh_size>0) || throw(ErrorException(
        "$caller: Curve[$curve] has an unrepresentable interpolated Point size"))
    return mesh_size
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
    xs=Float64[];ys=Float64[];mesh_sizes=Float64[]
    segs=Tuple{Int,Int}[]
    index=Dict{Int,Int}()
    for loop_id in m.surfaces[t]
        loop_idx=Int[]
        for signed in m.loops[loop_id]
            curve=abs(signed);a,b=m.curves[curve]
            for parameter in _surface_curve_parameters(forced,curve,signed)
                vertex=if parameter==0
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,a,caller)
                elseif parameter==1
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,b,caller)
                else
                    point=_periodic_curve_point(m,curve,parameter,caller)
                    abs(point[3])<=1e-12 || throw(ArgumentError(
                        "$caller: Curve[$curve] subdivision is not planar in z=0"))
                    _add_surface_curve_point!(
                        xs,ys,mesh_sizes,index,m,point,
                        _surface_curve_mesh_size(m,curve,parameter,caller),
                        curve,caller)
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
            _add_surface_point!(xs,ys,mesh_sizes,index,m,etag,caller)
        elseif edim!=1
            throw(ArgumentError(
                "$caller: unsupported embedding dimension $edim"))
        end
    end
    for (edim,etag) in embedded
        edim==1 || continue
        haskey(m.curves,etag) || throw(ArgumentError(
            "$caller: unknown embedded Curve[$etag]"))
        a,b=m.curves[etag]
        parameters=get(forced,etag,nothing)
        curve_nodes=Int[]
        if parameters===nothing
            push!(curve_nodes,_add_surface_point!(
                xs,ys,mesh_sizes,index,m,a,caller))
            push!(curve_nodes,_add_surface_point!(
                xs,ys,mesh_sizes,index,m,b,caller))
        else
            for parameter in parameters
                vertex=if parameter==0
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,a,caller)
                elseif parameter==1
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,b,caller)
                else
                    point=_periodic_curve_point(m,etag,parameter,caller)
                    abs(point[3])<=1e-12 || throw(ArgumentError(
                        "$caller: embedded Curve[$etag] subdivision is not " *
                        "planar in z=0"))
                    _add_surface_curve_point!(
                        xs,ys,mesh_sizes,index,m,point,
                        _surface_curve_mesh_size(m,etag,parameter,caller),
                        etag,caller)
                end
                push!(curve_nodes,vertex)
            end
        end
        length(curve_nodes)>=2 || throw(ArgumentError(
            "$caller: embedded Curve[$etag] needs two distinct endpoints"))
        for segment_index in 1:(length(curve_nodes)-1)
            first_node=curve_nodes[segment_index]
            second_node=curve_nodes[segment_index+1]
            first_node==second_node && throw(ArgumentError(
                "$caller: embedded Curve[$etag] has coincident subdivision nodes"))
            push!(internal,(first_node,second_node))
        end
    end
    return xs,ys,mesh_sizes,segs,embedded,internal
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

function _curve_parameter_nodes(m::GeoModel,mesh::Mesh,curve::Int,eligible_nodes,
                                eligible_edges,atol::Float64,
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
        eligible_nodes[node] || continue
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
        "$caller: Curve[$curve] is not represented by a two-node mesh-edge chain"))
    parameter_tolerance=max(128eps(Float64),geometric_tolerance/length1)
    first(entries)[1]<=parameter_tolerance &&
        1-last(entries)[1]<=parameter_tolerance || throw(ErrorException(
            "$caller: Curve[$curve] mesh chain does not reach both endpoints"))
    for index in 1:(length(entries)-1)
        first_node=Int32(entries[index][2])
        second_node=Int32(entries[index+1][2])
        first_node!=second_node || throw(ErrorException(
            "$caller: Curve[$curve] repeats a mesh node"))
        key=first_node<second_node ? (first_node,second_node) :
                                     (second_node,first_node)
        key in eligible_edges || throw(ErrorException(
            "$caller: Curve[$curve] nodes do not form a mesh-edge chain"))
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
    _,embedded_curve_tags=_model_surface_embedding_tags(m,t,caller)
    surface_curves=union(boundary_curves,Set(embedded_curve_tags))
    constraints=ModelPeriodicConstraint[]
    for constraint in model_periodic_constraints(m)
        constraint.dim==1 || continue
        slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
        slave_present=slave in surface_curves
        master_present=master in surface_curves
        slave_present==master_present || throw(ArgumentError(
            "$caller: periodic Curve[$slave]/Curve[$master] relation " *
            "has only one entity on Surface[$t]"))
        slave_present && push!(constraints,constraint)
    end
    return constraints
end

function _synchronize_periodic_parameters!(forced,m::GeoModel,mesh::Mesh,
                                           constraints)
    mesh_edges=_model_projection_triangle_edges(mesh)
    all_nodes=trues(nnodes(mesh))
    raw_parameters=Dict{Int,Vector{Float64}}()
    curve_tolerances=Dict{Int,Float64}()
    relations=Tuple{ModelPeriodicConstraint,Float64}[]
    for constraint in constraints
        slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
        master_entries,master_tolerance=_curve_parameter_nodes(
            m,mesh,master,all_nodes,mesh_edges,constraint.atol,
            "mesh_model_surface")
        slave_entries,slave_tolerance=_curve_parameter_nodes(
            m,mesh,slave,all_nodes,mesh_edges,constraint.atol,
            "mesh_model_surface")
        tolerance=max(master_tolerance,slave_tolerance)
        for (curve,entries) in ((master,master_entries),(slave,slave_entries))
            values=get!(()->Float64[],raw_parameters,curve)
            for (parameter,_) in entries
                push!(values,parameter)
            end
            curve_tolerances[curve]=min(
                get(curve_tolerances,curve,Inf),tolerance)
        end
        push!(relations,(constraint,tolerance))
    end

    adjacency=Dict{Int,Vector{Int}}()
    for (constraint,_) in relations
        slave=Int(constraint.slave_entity)
        master=Int(constraint.master_entity)
        push!(get!(()->Int[],adjacency,slave),master)
        push!(get!(()->Int[],adjacency,master),slave)
    end
    component_tolerances=Dict{Int,Float64}()
    visited=Set{Int}()
    for start in sort!(collect(keys(adjacency)))
        start in visited && continue
        component=Int[]
        stack=Int[start]
        tolerance=Inf
        while !isempty(stack)
            curve=pop!(stack)
            curve in visited && continue
            push!(visited,curve)
            push!(component,curve)
            tolerance=min(tolerance,curve_tolerances[curve])
            append!(stack,adjacency[curve])
        end
        for curve in component
            component_tolerances[curve]=tolerance
        end
    end

    parameters=Dict{Int,Vector{Float64}}()
    for curve in sort!(collect(keys(raw_parameters)))
        values=Float64[]
        tolerance=component_tolerances[curve]
        for parameter in raw_parameters[curve]
            _insert_periodic_parameter!(values,parameter,tolerance)
        end
        parameters[curve]=values
    end

    converged=false
    for _ in 1:(length(relations)+1)
        propagated=false
        for (constraint,_) in relations
            slave=Int(constraint.slave_entity)
            master=Int(constraint.master_entity)
            tolerance=component_tolerances[slave]
            master_values=copy(parameters[master])
            slave_values=copy(parameters[slave])
            for parameter in master_values
                mapped=constraint.reversed ? 1-parameter : parameter
                propagated|=_insert_periodic_parameter!(
                    parameters[slave],mapped,tolerance)
            end
            for parameter in slave_values
                mapped=constraint.reversed ? 1-parameter : parameter
                propagated|=_insert_periodic_parameter!(
                    parameters[master],mapped,tolerance)
            end
        end
        if !propagated
            converged=true
            break
        end
    end
    converged || throw(ErrorException(
        "mesh_model_surface: periodic parameter graph did not converge"))

    changed=false
    for curve in sort!(collect(keys(parameters)))
        curve_forced=get!(()->Float64[0,1],forced,curve)
        tolerance=component_tolerances[curve]
        for parameter in parameters[curve]
            changed|=_insert_periodic_parameter!(
                curve_forced,parameter,tolerance)
        end
    end
    return changed
end

function _model_periodic_curve_nodes(m::GeoModel,mesh::Mesh,
                                     constraint::ModelPeriodicConstraint)
    mesh_edges=_model_projection_triangle_edges(mesh)
    all_nodes=trues(nnodes(mesh))
    slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
    master_entries,master_tolerance=_curve_parameter_nodes(
        m,mesh,master,all_nodes,mesh_edges,constraint.atol,
        "model_periodic_nodes")
    slave_entries,slave_tolerance=_curve_parameter_nodes(
        m,mesh,slave,all_nodes,mesh_edges,constraint.atol,
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

function _model_periodic_nodes(m::GeoModel,mesh::Mesh,
                               constraint::ModelPeriodicConstraint)
    constraint.dim==1 && return _model_periodic_curve_nodes(
        m,mesh,constraint)
    constraint.dim==2 && return _model_periodic_surface_nodes(
        m,mesh,constraint)
    throw(ArgumentError(
        "model_periodic_nodes: unsupported periodic dimension $(constraint.dim)"))
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
    ordered=_model_periodic_constraint_order(constraints,caller)
    for constraint in ordered
        mapping=_model_periodic_nodes(m,output,constraint)
        output=periodic_identify_affine(
            output,constraint.affine,mapping.master_nodes,mapping.slave_nodes;
            atol=constraint.atol)
    end
    for constraint in ordered
        mapping=_model_periodic_nodes(m,output,constraint)
        _model_mapping_matches(
            output,constraint,mapping,caller;exact=true) || throw(ErrorException(
            "$caller: stored periodic constraints do not share an exact curve-node solution"))
    end
    return output
end

"""
    model_periodic_nodes(model, mesh, dim, slave_entity)

Return the master entity, compact slave/master node arrays, and affine transform
for one meshed curve or planar boundary-surface relation as a named tuple. Curve
relations require a synchronized boundary or embedded discretization produced by
[`mesh_model_surface`](@ref). Surface relations require a tetrahedron mesh of an
explicit planar-shell volume produced by [`mesh_model_volume`](@ref).
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
    return Int32.(_model_surface_boundary_curves(
        m,surface,"model_to_mixed";orient_holes=true))
end

function _model_surface_embedding_tags(
    m::GeoModel,surface::Int,caller::AbstractString)
    point_tags=Int[];curve_tags=Int[]
    for (embedded_dim,embedded_tag) in
            get(m.embeds,(2,surface),NTuple{2,Int}[])
        if embedded_dim==0
            haskey(m.points,embedded_tag) || throw(ArgumentError(
                "$caller: unknown embedded Point[$embedded_tag]"))
            push!(point_tags,embedded_tag)
        elseif embedded_dim==1
            haskey(m.curves,embedded_tag) || throw(ArgumentError(
                "$caller: unknown embedded Curve[$embedded_tag]"))
            push!(curve_tags,embedded_tag)
        else
            throw(ArgumentError(
                "$caller: unsupported Surface[$surface] embedding dimension " *
                "$embedded_dim"))
        end
    end
    return sort!(unique!(point_tags)),sort!(unique!(curve_tags))
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
ownership; and the surface's stored periodic curve links. It can be
written with [`Tessella.Elements.write_mixed_msh`](@ref) without dropping the
supported metadata. MSH4 retains every physical membership and Gmsh's
curve-in-surface relation. Gmsh 4.15.2 does not serialize Point-In-Surface as an
entity relation, but the classified point node and point element are retained.
MSH2 retains the lowest physical tag as its single legacy element membership.

`mesh` must be a validated, segment-free triangle mesh whose boundary and
embedded chains represent the selected [`mesh_model_surface`](@ref) geometry.
Stored periodic curves, whether boundary or embedded, must already be exactly snapped.
Periodic endpoint relations are emitted as a deterministic spanning forest when
curve directions share corners, satisfying the MSH one-master-per-slave entity
constraint while retaining every curve link.
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

    embedded_points,embedded_curves=
        _model_surface_embedding_tags(m,surface,caller)

    constraints=_surface_periodic_constraints(m,surface,caller)
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

@inline _model_projection_coordinate_key(value::Float64)=
    value==0 ? 0.0 : value

function _model_surface_projection(
    coordinates::Vector{NTuple{3,Float64}},point_tags::Vector{Int},
    surface::Int,caller::AbstractString)
    length(coordinates)>=3 || throw(ArgumentError(
        "$caller: Surface[$surface] needs at least three points"))
    anchor=coordinates[1]
    second_index=findfirst(point->point!=anchor,@view coordinates[2:end])
    second_index===nothing && throw(ArgumentError(
        "$caller: Surface[$surface] has only one distinct coordinate"))
    second_index+=1
    second=coordinates[second_index]
    projection=nothing
    plane_third=nothing
    for candidate in coordinates
        for axes in ((1,2),(1,3),(2,3))
            pa=(anchor[axes[1]],anchor[axes[2]])
            pb=(second[axes[1]],second[axes[2]])
            pc=(candidate[axes[1]],candidate[axes[2]])
            if orient2(pa,pb,pc)!=0
                projection=axes
                plane_third=candidate
                break
            end
        end
        projection===nothing || break
    end
    projection===nothing && throw(ArgumentError(
        "$caller: Surface[$surface] points are collinear"))
    third=plane_third::NTuple{3,Float64}
    for (index,coordinate) in pairs(coordinates)
        orient3(anchor,second,third,coordinate)==0 || throw(ArgumentError(
            "$caller: Point[$(point_tags[index])] is not coplanar with " *
            "Surface[$surface]"))
    end
    return anchor,second,third,projection::Tuple{Int,Int}
end

function _model_planar_surface_mesh(
    m::GeoModel,surface::Int,caller::AbstractString;
    include_embeddings::Bool=true)
    haskey(m.surfaces,surface) || throw(ArgumentError(
        "$caller: unknown Surface[$surface]"))
    point_tags=Int[]
    point_index=Dict{Int,Int}()
    function add_point_tag!(point::Int)
        haskey(m.points,point) || throw(ArgumentError(
            "$caller: Surface[$surface] references unknown Point[$point]"))
        return get!(point_index,point) do
            push!(point_tags,point)
            length(point_tags)
        end
    end
    segments=Tuple{Int,Int}[]
    for loop in m.surfaces[surface]
        haskey(m.loops,loop) || throw(ArgumentError(
            "$caller: Surface[$surface] references unknown Loop[$loop]"))
        loop_points=_loop_points(m,loop)
        length(loop_points)>=3 || throw(ArgumentError(
            "$caller: Surface[$surface] Loop[$loop] needs at least three points"))
        indices=Int[add_point_tag!(point) for point in loop_points]
        for index in eachindex(indices)
            push!(segments,(indices[index],indices[mod1(index+1,length(indices))]))
        end
    end
    internal_segments=Tuple{Int,Int}[]
    if include_embeddings
        embedded_points,embedded_curves=
            _model_surface_embedding_tags(m,surface,caller)
        boundary_curves=Set(_model_projection_surface_curves(m,surface))
        overlap=sort!(collect(intersect(boundary_curves,Set(embedded_curves))))
        isempty(overlap) || throw(ArgumentError(
            "$caller: Surface[$surface] Curve tags $overlap cannot be both " *
            "bounding and embedded"))
        for point in embedded_points
            add_point_tag!(point)
        end
        for curve in embedded_curves
            start_point,stop_point=m.curves[curve]
            first_index=add_point_tag!(start_point)
            second_index=add_point_tag!(stop_point)
            push!(internal_segments,(first_index,second_index))
        end
    end
    coordinates=NTuple{3,Float64}[m.points[point] for point in point_tags]
    _,_,_,projection=_model_surface_projection(
        coordinates,point_tags,surface,caller)
    first_axis,second_axis=projection
    xs=Float64[_model_projection_coordinate_key(point[first_axis])
               for point in coordinates]
    ys=Float64[_model_projection_coordinate_key(point[second_axis])
               for point in coordinates]
    coordinate_map=Dict{Tuple{Float64,Float64},NTuple{3,Float64}}()
    for (index,coordinate) in pairs(coordinates)
        key=(xs[index],ys[index])
        existing=get(coordinate_map,key,nothing)
        (existing===nothing || existing==coordinate) || throw(ArgumentError(
            "$caller: Surface[$surface] projection merges distinct coordinates"))
        coordinate_map[key]=coordinate
    end
    triangulation=constrained_delaunay(
        xs,ys,segments;internal_segments=internal_segments)
    local_mesh=to_mesh(
        triangulation;interior=classify_interior(triangulation))
    ntris(local_mesh)>0 || throw(ArgumentError(
        "$caller: Surface[$surface] has no meshed interior"))
    output_coordinates=Matrix{Float64}(undef,3,nnodes(local_mesh))
    for node in 1:nnodes(local_mesh)
        key=(_model_projection_coordinate_key(local_mesh.coords[1,node]),
             _model_projection_coordinate_key(local_mesh.coords[2,node]))
        coordinate=get(coordinate_map,key,nothing)
        coordinate===nothing && throw(ErrorException(
            "$caller: Surface[$surface] triangulation introduced an " *
            "unmapped coordinate"))
        output_coordinates[:,node].=coordinate
    end
    output=Mesh(output_coordinates;tris=local_mesh.tris)
    diagnostic=validate(output)
    diagnostic.ok || throw(ErrorException(
        "$caller: Surface[$surface] triangulation is invalid — " *
        join(diagnostic.messages,"; ")))
    return output
end

function _model_periodic_surface_mesh(
    m::GeoModel,master_mesh::Mesh,constraint::ModelPeriodicConstraint,
    caller::AbstractString)
    constraint.dim==2 || throw(ErrorException(
        "$caller: internal periodic surface mesh received dimension " *
        "$(constraint.dim)"))
    slave=Int(constraint.slave_entity)
    master=Int(constraint.master_entity)
    point_map=_model_periodic_surface_point_map(
        m,slave,master,constraint.affine,constraint.atol,caller;
        include_embeddings=true)
    master_points=_model_periodic_surface_points(
        m,master,caller;include_embeddings=true)
    coordinate_points=Dict{NTuple{3,Float64},Int}()
    for point in master_points
        coordinate=m.points[point]
        key=ntuple(
            axis->_model_projection_coordinate_key(coordinate[axis]),3)
        haskey(coordinate_points,key) && throw(ArgumentError(
            "$caller: Surface[$master] has multiple point tags at $key"))
        coordinate_points[key]=point
    end
    output_coordinates=Matrix{Float64}(undef,3,nnodes(master_mesh))
    for node in 1:nnodes(master_mesh)
        key=ntuple(axis->_model_projection_coordinate_key(
            master_mesh.coords[axis,node]),3)
        master_point=get(coordinate_points,key,nothing)
        master_point===nothing && throw(ErrorException(
            "$caller: periodic master Surface[$master] mesh introduced " *
            "an unmapped node at $key"))
        slave_point=point_map[master_point]
        output_coordinates[:,node].=m.points[slave_point]
    end
    output=Mesh(output_coordinates;tris=master_mesh.tris)
    diagnostic=validate(output)
    diagnostic.ok || throw(ArgumentError(
        "$caller: synchronized periodic Surface[$slave] mesh is invalid — " *
        join(diagnostic.messages,"; ")))
    return output
end

function _model_loop_position2(point,polygon)
    inside=false
    previous=last(polygon)
    for current in polygon
        orientation=orient2(previous,current,point)
        if orientation==0 &&
                min(previous[1],current[1])<=point[1]<=max(previous[1],current[1]) &&
                min(previous[2],current[2])<=point[2]<=max(previous[2],current[2])
            return 2
        end
        upward=previous[2]<=point[2]<current[2] && orientation>0
        downward=current[2]<=point[2]<previous[2] && orientation<0
        (upward || downward) && (inside=!inside)
        previous=current
    end
    return inside ? 1 : 0
end

function _model_surface_contains2(point,polygons)
    outer=_model_loop_position2(point,first(polygons))
    outer==0 && return false
    outer==2 && return true
    for hole in Iterators.drop(polygons,1)
        position=_model_loop_position2(point,hole)
        position==2 && return true
        position==1 && return false
    end
    return true
end

@inline function _model_mesh_coordinate(mesh::Mesh,node::Integer)
    return (mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
end

@inline function _model_mean3(first::Float64,second::Float64,third::Float64)
    scale=max(abs(first),abs(second),abs(third))
    scale==0 && return 0.0
    value=scale*((first/scale+second/scale+third/scale)/3)
    return _model_projection_coordinate_key(value)
end

function _model_projection_boundary_surface_faces!(
    claimed_faces::Set{NTuple{3,Int32}},mesh_boundary_faces,
    m::GeoModel,mesh::Mesh,surface::Int,caller::AbstractString)
    point_tags=Int[]
    polygons=Vector{NTuple{2,Float64}}[]
    for loop in m.surfaces[surface]
        loop_points=_loop_points(m,loop)
        append!(point_tags,loop_points)
    end
    unique!(point_tags)
    coordinates=NTuple{3,Float64}[m.points[point] for point in point_tags]
    anchor,second,third,projection=_model_surface_projection(
        coordinates,point_tags,surface,caller)
    first_axis,second_axis=projection
    for loop in m.surfaces[surface]
        polygon=NTuple{2,Float64}[]
        for point in _loop_points(m,loop)
            coordinate=m.points[point]
            push!(polygon,
                (_model_projection_coordinate_key(coordinate[first_axis]),
                 _model_projection_coordinate_key(coordinate[second_axis])))
        end
        push!(polygons,polygon)
    end
    anchor2=(anchor[first_axis],anchor[second_axis])
    second2=(second[first_axis],second[second_axis])
    third2=(third[first_axis],third[second_axis])
    target_orientation=orient2(anchor2,second2,third2)
    target_orientation!=0 || throw(ErrorException(
        "$caller: Surface[$surface] projection lost its plane orientation"))
    surface_mesh=_model_planar_surface_mesh(
        m,surface,caller;include_embeddings=true)
    target_area=sum(triangle_area(
        _model_mesh_coordinate(surface_mesh,surface_mesh.tris[1,cell]),
        _model_mesh_coordinate(surface_mesh,surface_mesh.tris[2,cell]),
        _model_mesh_coordinate(surface_mesh,surface_mesh.tris[3,cell]))
        for cell in 1:ntris(surface_mesh))
    (isfinite(target_area) && target_area>0) || throw(ArgumentError(
        "$caller: Surface[$surface] has no finite positive area"))

    output=NTuple{3,Int32}[]
    covered_area=0.0
    for face in sort!(collect(mesh_boundary_faces))
        first_coordinate=_model_mesh_coordinate(mesh,face[1])
        second_coordinate=_model_mesh_coordinate(mesh,face[2])
        third_coordinate=_model_mesh_coordinate(mesh,face[3])
        orient3(anchor,second,third,first_coordinate)==0 || continue
        orient3(anchor,second,third,second_coordinate)==0 || continue
        orient3(anchor,second,third,third_coordinate)==0 || continue
        projected_vertices=ntuple(slot->begin
            coordinate=slot==1 ? first_coordinate :
                       slot==2 ? second_coordinate : third_coordinate
            (_model_projection_coordinate_key(coordinate[first_axis]),
             _model_projection_coordinate_key(coordinate[second_axis]))
        end,3)
        all(vertex->_model_surface_contains2(vertex,polygons),
            projected_vertices) || continue
        centroid=(_model_mean3(projected_vertices[1][1],
                               projected_vertices[2][1],
                               projected_vertices[3][1]),
                  _model_mean3(projected_vertices[1][2],
                               projected_vertices[2][2],
                               projected_vertices[3][2]))
        _model_surface_contains2(centroid,polygons) || continue
        face in claimed_faces && throw(ArgumentError(
            "$caller: mesh face $face belongs to multiple model surfaces"))
        face_orientation=orient2(projected_vertices...)
        face_orientation!=0 || throw(ArgumentError(
            "$caller: boundary Surface[$surface] contains degenerate projected " *
            "mesh face $face"))
        oriented=face_orientation==target_orientation ? face :
                 (face[1],face[3],face[2])
        push!(claimed_faces,face)
        push!(output,oriented)
        covered_area+=triangle_area(
            first_coordinate,second_coordinate,third_coordinate)
    end
    isempty(output) && throw(ArgumentError(
        "$caller: boundary Surface[$surface] has no tetrahedron faces"))
    (isfinite(covered_area) &&
     abs(covered_area-target_area)<=1e-6*target_area) || throw(ArgumentError(
        "$caller: boundary Surface[$surface] mesh area $covered_area does not " *
        "match modeled area $target_area"))
    return output
end

function _model_affine_node_pairs(
    mesh::Mesh,constraint::ModelPeriodicConstraint,master_nodes_raw,
    slave_nodes_raw,caller::AbstractString)
    master_nodes=sort!(unique!(Int[Int(node) for node in master_nodes_raw]))
    slave_nodes=sort!(unique!(Int[Int(node) for node in slave_nodes_raw]);
        by=node->(mesh.coords[1,node],mesh.coords[2,node],
                  mesh.coords[3,node],node))
    length(master_nodes)==length(slave_nodes) || throw(ArgumentError(
        "$caller: periodic $(_model_periodic_entity_label(constraint.dim))" *
        "[$(constraint.slave_entity)] and " *
        "$(_model_periodic_entity_label(constraint.dim))" *
        "[$(constraint.master_entity)] node counts differ"))
    coefficients,translation,_=_transform_homogeneous(
        constraint.affine,caller;name="stored affine transform")
    slave_x=Float64[mesh.coords[1,node] for node in slave_nodes]
    used=Set{Int}()
    mapped_slaves=Vector{Int32}(undef,length(master_nodes))
    mapped_masters=Vector{Int32}(undef,length(master_nodes))
    for (index,master_node) in pairs(master_nodes)
        master=_model_mesh_coordinate(mesh,master_node)
        expected=_model_affine_point(
            coefficients,translation,master,caller,master_node)
        first_candidate=searchsortedfirst(
            slave_x,expected[1]-constraint.atol)
        last_candidate=searchsortedlast(
            slave_x,expected[1]+constraint.atol)
        matches=Int[]
        for position in first_candidate:last_candidate
            position in eachindex(slave_nodes) || continue
            slave_node=slave_nodes[position]
            slave_node in used && continue
            _model_point_distance(
                _model_mesh_coordinate(mesh,slave_node),expected)<=
                constraint.atol && push!(matches,slave_node)
        end
        length(matches)==1 || throw(ArgumentError(
            "$caller: affine image of master node $master_node matches " *
            "$(length(matches)) unused slave nodes; expected one"))
        slave_node=only(matches)
        push!(used,slave_node)
        mapped_masters[index]=Int32(master_node)
        mapped_slaves[index]=Int32(slave_node)
    end
    length(used)==length(slave_nodes) || throw(ErrorException(
        "$caller: internal periodic node matching was incomplete"))
    return (master_entity=Int(constraint.master_entity),
            slave_nodes=mapped_slaves,master_nodes=mapped_masters,
            affine=constraint.affine)
end

function _model_periodic_surface_nodes(
    m::GeoModel,mesh::Mesh,constraint::ModelPeriodicConstraint)
    caller="model_periodic_nodes"
    ntets(mesh)>0 || throw(ArgumentError(
        "$caller: periodic Surface mapping requires a tetrahedron mesh"))
    slave=Int(constraint.slave_entity)
    master=Int(constraint.master_entity)
    boundary=Set(first(boundary_faces(mesh.tets)))
    slave_faces=_model_projection_boundary_surface_faces!(
        Set{NTuple{3,Int32}}(),boundary,m,mesh,slave,caller)
    master_faces=_model_projection_boundary_surface_faces!(
        Set{NTuple{3,Int32}}(),boundary,m,mesh,master,caller)
    slave_nodes,_=_model_projection_face_topology(slave_faces)
    master_nodes,_=_model_projection_face_topology(master_faces)
    mapping=_model_affine_node_pairs(
        mesh,constraint,master_nodes,slave_nodes,caller)
    node_map=Dict(master_node=>slave_node for (master_node,slave_node) in
        zip(mapping.master_nodes,mapping.slave_nodes))
    mapped_master_faces=Set{NTuple{3,Int32}}()
    for face in master_faces
        mapped=ntuple(slot->node_map[face[slot]],3)
        push!(mapped_master_faces,_model_projection_face_key(mapped))
    end
    slave_face_keys=Set(_model_projection_face_key(face)
                        for face in slave_faces)
    mapped_master_faces==slave_face_keys || throw(ArgumentError(
        "$caller: periodic Surface[$slave]/Surface[$master] face topology " *
        "does not match under the affine node map"))
    return mapping
end

function _model_projection_volume_surface_faces!(
    claimed_faces::Set{NTuple{3,Int32}},m::GeoModel,mesh::Mesh,
    surface::Int,caller::AbstractString)
    targets=NTuple{3,NTuple{3,Float64}}[]
    loops=m.surfaces[surface]
    if length(loops)==1
        point_tags=_loop_points(m,only(loops))
        length(point_tags)>=3 || throw(ArgumentError(
            "$caller: Surface[$surface] needs at least three points"))
        coordinates=NTuple{3,Float64}[m.points[point] for point in point_tags]
        for index in 2:(length(coordinates)-1)
            push!(targets,(coordinates[1],coordinates[index],
                           coordinates[index+1]))
        end
    else
        surface_mesh=_model_planar_surface_mesh(
            m,surface,caller;include_embeddings=true)
        for cell in 1:ntris(surface_mesh)
            push!(targets,ntuple(slot->begin
                node=surface_mesh.tris[slot,cell]
                (surface_mesh.coords[1,node],surface_mesh.coords[2,node],
                 surface_mesh.coords[3,node])
            end,3))
        end
    end
    local_faces=Set{NTuple{3,Int32}}()
    output=NTuple{3,Int32}[]
    for (index,(first_coordinate,second_coordinate,third_coordinate)) in
            pairs(targets)
        faces,target,covered=_mesh_covering_faces3(
            mesh,first_coordinate,second_coordinate,third_coordinate)
        (isfinite(target) && target>0 && isfinite(covered) &&
         abs(covered-target)<=1e-6*target) || throw(ArgumentError(
            "$caller: Surface[$surface] triangle $index " *
            "is not represented by tetrahedron faces"))
        for face in faces
            key=_model_projection_face_key(face)
            key in local_faces && throw(ArgumentError(
                "$caller: Surface[$surface] triangles overlap " *
                "on mesh face $key"))
            key in claimed_faces && throw(ArgumentError(
                "$caller: mesh face $key belongs to multiple model surfaces"))
            push!(local_faces,key)
            push!(claimed_faces,key)
            push!(output,face)
        end
    end
    isempty(output) && throw(ArgumentError(
        "$caller: Surface[$surface] has no tetrahedron faces"))
    return output
end

function _model_volume_boundary_surfaces(
    m::GeoModel,volume::Int,caller::AbstractString;
    orient_cavities::Bool=true)
    boundaries=Int[]
    seen=Set{Int}()
    for (shell_index,shell) in pairs(m.volumes[volume])
        haskey(m.surface_loops,shell) || throw(ArgumentError(
            "$caller: Volume[$volume] references unknown Surface Loop[$shell]"))
        surfaces=m.surface_loops[shell]
        _validate_surface_loop(m,surfaces,caller)
        shell_sign=orient_cavities && shell_index>1 ? -1 : 1
        for signed_surface in surfaces
            surface=abs(signed_surface)
            surface in seen && throw(ArgumentError(
                "$caller: Surface[$surface] belongs to multiple Volume[$volume] shells"))
            push!(seen,surface)
            push!(boundaries,shell_sign*signed_surface)
        end
    end
    return boundaries
end

function _model_volume_embedding_inventory(
    m::GeoModel,volume::Int,caller::AbstractString)
    point_tags=Int[];curve_tags=Int[];surface_tags=Int[]
    surface_embedded_points=Dict{Int,Vector{Int}}()
    surface_embedded_curves=Dict{Int,Vector{Int}}()
    boundary_surfaces=_model_volume_boundary_surfaces(m,volume,caller)
    boundary_set=Set(abs.(boundary_surfaces))
    function register_surface!(surface::Int)
        haskey(m.surfaces,surface) || throw(ArgumentError(
            "$caller: unknown Surface[$surface]"))
        nested_points,nested_curves=
            _model_surface_embedding_tags(m,surface,caller)
        boundary_curves=_model_projection_surface_curves(m,surface)
        overlap=sort!(collect(intersect(Set(boundary_curves),Set(nested_curves))))
        isempty(overlap) || throw(ArgumentError(
            "$caller: Surface[$surface] Curve tags $overlap cannot be both " *
            "bounding and embedded"))
        surface_embedded_points[surface]=nested_points
        surface_embedded_curves[surface]=nested_curves
        push!(surface_tags,surface)
        append!(point_tags,nested_points)
        append!(curve_tags,boundary_curves)
        append!(curve_tags,nested_curves)
        return nothing
    end
    for signed_surface in boundary_surfaces
        register_surface!(abs(signed_surface))
    end
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
            embedded_tag in boundary_set && throw(ArgumentError(
                "$caller: boundary Surface[$embedded_tag] cannot also be " *
                "embedded in Volume[$volume]"))
            embedded_tag in surface_tags && throw(ArgumentError(
                "$caller: Volume[$volume] repeats embedded Surface[$embedded_tag]"))
            haskey(m.surfaces,embedded_tag) || throw(ArgumentError(
                "$caller: unknown embedded Surface[$embedded_tag]"))
            register_surface!(embedded_tag)
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
    return point_tags,curve_tags,surface_tags,
           surface_embedded_points,surface_embedded_curves,boundary_surfaces
end

function _model_projection_face_topology(faces)
    nodes=Set{Int32}()
    edges=Set{NTuple{2,Int32}}()
    for (first_node,second_node,third_node) in faces
        push!(nodes,first_node,second_node,third_node)
        for (a,b) in ((first_node,second_node),(second_node,third_node),
                      (third_node,first_node))
            push!(edges,a<b ? (a,b) : (b,a))
        end
    end
    return nodes,edges
end

function _model_projection_validate_nested_surface(
    surface::Int,nodes,edges,point_nodes,curve_entries,
    nested_points,nested_curves,caller::AbstractString)
    for point in nested_points
        node=point_nodes[point]
        node in nodes || throw(ArgumentError(
            "$caller: nested Point[$point] is not a node of embedded " *
            "Surface[$surface]"))
    end
    for curve in nested_curves
        entries=curve_entries[curve]
        for index in 1:(length(entries)-1)
            first_node=Int32(entries[index][2])
            second_node=Int32(entries[index+1][2])
            edge=first_node<second_node ? (first_node,second_node) :
                                          (second_node,first_node)
            edge in edges || throw(ArgumentError(
                "$caller: nested Curve[$curve] edge $edge is not an edge of " *
                "embedded Surface[$surface]"))
        end
    end
    return nothing
end

function _model_periodic_surface_boundary_maps(
    m::GeoModel,constraint::ModelPeriodicConstraint,caller::AbstractString)
    slave=Int(constraint.slave_entity)
    master=Int(constraint.master_entity)
    point_map=_model_periodic_surface_point_map(
        m,slave,master,constraint.affine,constraint.atol,caller;
        include_embeddings=false)
    slave_curves=_model_projection_surface_curves(m,slave)
    master_curves=_model_projection_surface_curves(m,master)
    length(slave_curves)==length(master_curves) || throw(ArgumentError(
        "$caller: periodic surfaces have different boundary-curve counts"))
    slave_signatures=Dict{NTuple{2,Int},Int}()
    for curve in slave_curves
        first_point,second_point=m.curves[curve]
        signature=first_point<second_point ?
            (first_point,second_point) : (second_point,first_point)
        haskey(slave_signatures,signature) && throw(ArgumentError(
            "$caller: Surface[$slave] repeats boundary edge $signature"))
        slave_signatures[signature]=curve
    end
    curve_map=Dict{Int,Int}()
    used=Set{Int}()
    for master_curve in master_curves
        first_point,second_point=m.curves[master_curve]
        first_mapped=point_map[first_point]
        second_mapped=point_map[second_point]
        signature=first_mapped<second_mapped ?
            (first_mapped,second_mapped) :
            (second_mapped,first_mapped)
        slave_curve=get(slave_signatures,signature,nothing)
        slave_curve===nothing && throw(ArgumentError(
            "$caller: affine image of Curve[$master_curve] has no boundary " *
            "curve on Surface[$slave]"))
        slave_curve in used && throw(ArgumentError(
            "$caller: multiple master curves map to Curve[$slave_curve]"))
        curve_map[master_curve]=slave_curve
        push!(used,slave_curve)
    end
    length(used)==length(slave_curves) || throw(ErrorException(
        "$caller: internal periodic boundary-curve matching was incomplete"))
    return point_map,curve_map
end

function _model_periodic_add_relation!(
    relations::Dict{Tuple{Int,Int},Tuple{NTuple{16,Float64},Float64}},
    slave::Int,master::Int,constraint::ModelPeriodicConstraint,
    caller::AbstractString,label::AbstractString)
    key=(slave,master)
    existing=get(relations,key,nothing)
    if existing===nothing
        relations[key]=(constraint.affine,constraint.atol)
    else
        existing[1]==constraint.affine || throw(ArgumentError(
            "$caller: periodic $label[$slave]/$label[$master] is induced " *
            "by inconsistent surface transforms"))
        relations[key]=(existing[1],min(existing[2],constraint.atol))
    end
    return nothing
end

function _model_periodic_spanning_relations(
    relations::Dict{Tuple{Int,Int},Tuple{NTuple{16,Float64},Float64}})
    outgoing=Dict{Int,Vector{Tuple{Int,NTuple{16,Float64},Float64}}}()
    indegree=Dict{Int,Int}()
    entities=Set{Int}()
    for ((slave,master),(affine,atol)) in relations
        push!(get!(Vector{Tuple{Int,NTuple{16,Float64},Float64}},
                   outgoing,master),(slave,affine,atol))
        indegree[slave]=get(indegree,slave,0)+1
        get!(indegree,master,0)
        push!(entities,slave,master)
    end
    for edges in values(outgoing)
        sort!(edges;by=edge->(edge[1],edge[2],edge[3]))
    end
    roots=sort!(Int[entity for entity in entities
                    if get(indegree,entity,0)==0])
    append!(roots,sort!(collect(setdiff(entities,Set(roots)))))
    visited=Set{Int}()
    parents=Dict{Int,Tuple{Int,NTuple{16,Float64},Float64}}()
    for root in roots
        root in visited && continue
        push!(visited,root)
        queue=Int[root];head=1
        while head<=length(queue)
            master=queue[head];head+=1
            for (slave,affine,atol) in get(
                    outgoing,master,
                    Tuple{Int,NTuple{16,Float64},Float64}[])
                slave in visited && continue
                push!(visited,slave)
                parents[slave]=(master,affine,atol)
                push!(queue,slave)
            end
        end
    end
    visited==entities || throw(ErrorException(
        "internal periodic entity traversal was incomplete"))
    return parents
end

function _model_periodic_curve_entry_mapping(
    mesh::Mesh,slave::Int,master::Int,affine,atol::Float64,
    slave_entries,master_entries,caller::AbstractString)
    constraint=ModelPeriodicConstraint(
        1,Int32(slave),Int32(master),affine,false,atol)
    mapping=_model_affine_node_pairs(
        mesh,constraint,last.(master_entries),last.(slave_entries),caller)
    node_map=Dict(master_node=>slave_node for (master_node,slave_node) in
        zip(mapping.master_nodes,mapping.slave_nodes))
    function edge_set(entries,map)
        edges=Set{NTuple{2,Int32}}()
        for index in 1:(length(entries)-1)
            first_node=map[Int32(entries[index][2])]
            second_node=map[Int32(entries[index+1][2])]
            push!(edges,first_node<second_node ?
                        (first_node,second_node) :
                        (second_node,first_node))
        end
        return edges
    end
    identity_map=Dict(Int32(entry[2])=>Int32(entry[2])
                      for entry in slave_entries)
    edge_set(master_entries,node_map)==edge_set(slave_entries,identity_map) ||
        throw(ArgumentError(
            "$caller: periodic Curve[$slave]/Curve[$master] edge topology " *
            "does not match under the affine node map"))
    return mapping
end

function _model_projection_periodic_surface_links(
    m::GeoModel,mesh::Mesh,constraints,point_nodes,curve_entries,
    surface_nodes,caller::AbstractString)
    point_relations=
        Dict{Tuple{Int,Int},Tuple{NTuple{16,Float64},Float64}}()
    curve_relations=
        Dict{Tuple{Int,Int},Tuple{NTuple{16,Float64},Float64}}()
    for constraint in constraints
        point_map,curve_map=
            _model_periodic_surface_boundary_maps(m,constraint,caller)
        for (master,slave) in point_map
            _model_periodic_add_relation!(
                point_relations,slave,master,constraint,caller,"Point")
        end
        for (master,slave) in curve_map
            _model_periodic_add_relation!(
                curve_relations,slave,master,constraint,caller,"Curve")
        end
    end

    links=MixedPeriodicLink[]
    point_parents=_model_periodic_spanning_relations(point_relations)
    for slave in sort!(collect(keys(point_parents)))
        master,affine,atol=point_parents[slave]
        haskey(point_nodes,slave) && haskey(point_nodes,master) ||
            throw(ErrorException(
                "$caller: periodic point entity is absent from projection"))
        constraint=ModelPeriodicConstraint(
            0,Int32(slave),Int32(master),affine,false,atol)
        mapping=(master_entity=master,
                 slave_nodes=Int32[point_nodes[slave]],
                 master_nodes=Int32[point_nodes[master]],affine=affine)
        _model_mapping_matches(mesh,constraint,mapping,caller;exact=false)
        push!(links,MixedPeriodicLink(
            0,slave,master,mapping.slave_nodes,mapping.master_nodes;
            affine=affine))
    end
    curve_parents=_model_periodic_spanning_relations(curve_relations)
    for slave in sort!(collect(keys(curve_parents)))
        master,affine,atol=curve_parents[slave]
        haskey(curve_entries,slave) && haskey(curve_entries,master) ||
            throw(ErrorException(
                "$caller: periodic curve entity is absent from projection"))
        mapping=_model_periodic_curve_entry_mapping(
            mesh,slave,master,affine,atol,
            curve_entries[slave],curve_entries[master],caller)
        push!(links,MixedPeriodicLink(
            1,slave,master,mapping.slave_nodes,mapping.master_nodes;
            affine=affine))
    end
    for constraint in constraints
        slave=Int(constraint.slave_entity)
        master=Int(constraint.master_entity)
        mapping=_model_periodic_surface_nodes(m,mesh,constraint)
        Set(mapping.slave_nodes)==surface_nodes[slave] || throw(ErrorException(
            "$caller: periodic Surface[$slave] mapping omits classified nodes"))
        Set(mapping.master_nodes)==surface_nodes[master] || throw(ErrorException(
            "$caller: periodic Surface[$master] mapping omits classified nodes"))
        _model_mapping_matches(mesh,constraint,mapping,caller;exact=false)
        push!(links,MixedPeriodicLink(
            2,constraint.slave_entity,constraint.master_entity,
            mapping.slave_nodes,mapping.master_nodes;
            affine=constraint.affine))
    end
    return links
end

function _model_volume_to_mixed(m::GeoModel,mesh::Mesh,volume::Int)
    caller="model_to_mixed"
    haskey(m.volumes,volume) || throw(ArgumentError(
        "$caller: unknown Volume[$volume]"))
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

    explicit_geometry=isempty(m.volumes[volume]) ? nothing :
        _model_explicit_volume_geometry(m,volume,caller)
    domain_surface=explicit_geometry===nothing ?
        _volume_surface(m,volume,caller) : explicit_geometry.surface
    fills_volume,fill_reason=_certify_surface_fill(domain_surface,mesh)
    fills_volume || throw(ArgumentError(
        "$caller: input tetrahedron mesh does not fill Volume[$volume] — $fill_reason"))
    explicit_geometry===nothing || _model_certify_explicit_volume_semantics(
        mesh,volume,explicit_geometry.expected_volume,
        explicit_geometry.comparison_scale,caller)

    point_tags,curve_tags,surface_tags,
    surface_embedded_points,surface_embedded_curves,boundary_surfaces=
        _model_volume_embedding_inventory(m,volume,caller)
    periodic_surfaces=
        _model_volume_periodic_surface_constraints(m,volume,caller)
    boundary_surface_tags=abs.(boundary_surfaces)
    boundary_surface_set=Set(boundary_surface_tags)
    projection_surface_tags=vcat(
        boundary_surface_tags,
        Int[surface for surface in surface_tags
            if !(surface in boundary_surface_set)])
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
    surface_edges=Dict{Int,Set{NTuple{2,Int32}}}()
    mesh_boundary_faces=Set(first(boundary_faces(mesh.tets)))
    claimed_boundary_faces=Set{NTuple{3,Int32}}()
    for surface in projection_surface_tags
        faces=if surface in boundary_surface_set
            _model_projection_boundary_surface_faces!(
                claimed_faces,mesh_boundary_faces,m,mesh,surface,caller)
        else
            _model_projection_volume_surface_faces!(
                claimed_faces,m,mesh,surface,caller)
        end
        nodes,edges=_model_projection_face_topology(faces)
        surface_nodes[surface]=nodes
        surface_edges[surface]=edges
        for face in faces
            key=_model_projection_face_key(face)
            if surface in boundary_surface_set
                key in mesh_boundary_faces || throw(ArgumentError(
                    "$caller: boundary Surface[$surface] contains internal " *
                    "tetrahedron face $key"))
                push!(claimed_boundary_faces,key)
            end
            push!(surface_cells,face)
            push!(surface_entities,Int32(surface))
        end
    end
    if !isempty(boundary_surfaces) && claimed_boundary_faces!=mesh_boundary_faces
        missing=sort!(collect(setdiff(mesh_boundary_faces,claimed_boundary_faces)))
        extra=sort!(collect(setdiff(claimed_boundary_faces,mesh_boundary_faces)))
        detail=!isempty(missing) ? "unclassified boundary face $(first(missing))" :
                                  "non-boundary face $(first(extra)) was classified"
        throw(ArgumentError(
            "$caller: explicit Volume[$volume] shell does not exactly cover " *
            "the tetrahedron boundary — $detail"))
    end
    for surface in surface_tags
        _model_projection_validate_nested_surface(
            surface,surface_nodes[surface],surface_edges[surface],
            point_nodes,curve_entries,surface_embedded_points[surface],
            surface_embedded_curves[surface],caller)
    end
    periodic_links=_model_projection_periodic_surface_links(
        m,mesh,periodic_surfaces,point_nodes,curve_entries,
        surface_nodes,caller)

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
            boundaries=_model_projection_surface_boundaries(m,surface),
            embedded_curves=surface_embedded_curves[surface])
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
        physical_tags=volume_tags,boundaries=Int32.(boundary_surfaces))
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
        elementary_entities=block_entities,periodic_links=periodic_links)
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
boundary fills the selected native solid. Explicit volumes additionally require
their modeled surfaces to classify every tetrahedron boundary face exactly once.
The result emits boundary and embedded point/curve/surface cells and their entity
hierarchy, including nested Point/Line-In-Surface constraints. Periodic explicit
volume boundaries retain their surface maps and induced boundary point/curve forest.
Gmsh 4.15.2 has no
serialized volume-embedding relation; MSH4 classification, signed volume boundaries,
its Curve-In-Surface relation, and MSH2 cell ownership remain available.
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

include("SurfacePointSizing.jl")

function _mesh_model_surface_once(m::GeoModel,t::Int,forced,min_angle_deg,
                                  caller::AbstractString)
    xs,ys,mesh_sizes,segs,embedded,internal=
        _surface_pslg(m,t,forced,caller)
    T=constrained_delaunay(xs,ys,segs; internal_segments=internal)
    sizefn=_surface_point_size_field(T,xs,ys,mesh_sizes,t,caller)
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
curves. Point characteristic lengths are linearly interpolated over the
deterministic initial constrained triangulation and drive refinement. Coincident
PSLG inputs use the smaller constraint. Stored
straight-curve periodic relations synchronize boundary or embedded curve subdivisions
across each acyclic dependency graph. Bounded remeshing precedes topology-ordered affine
snapping, so a curve may be both a slave and a downstream master. The returned
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
                "$caller: periodic curve synchronization did not converge " *
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

@inline function _model_compensated_add(
    total::Float64,correction::Float64,value::Float64)
    adjusted=value-correction
    updated=total+adjusted
    return updated,(updated-total)-adjusted
end

function _model_oriented_shell_volume(
    coordinates::Vector{NTuple{3,Float64}},
    faces::Vector{NTuple{3,Int32}},shell::Int,caller::AbstractString)
    isempty(faces) && throw(ArgumentError(
        "$caller: Surface Loop[$shell] has no triangles"))
    incidence=Dict{NTuple{2,Int32},Vector{Tuple{Int32,Bool}}}()
    for (cell,(first_node,second_node,third_node)) in pairs(faces)
        for (start_node,stop_node) in
                ((first_node,second_node),(second_node,third_node),
                 (third_node,first_node))
            edge=start_node<stop_node ? (start_node,stop_node) :
                                        (stop_node,start_node)
            push!(get!(Vector{Tuple{Int32,Bool}},incidence,edge),
                  (Int32(cell),start_node<stop_node))
        end
    end
    adjacency=[Tuple{Int32,Bool}[] for _ in eachindex(faces)]
    for edge in sort!(collect(keys(incidence)))
        entries=incidence[edge]
        length(entries)==2 || throw(ArgumentError(
            "$caller: Surface Loop[$shell] triangle edge $edge has " *
            "incidence $(length(entries)); expected 2"))
        (first_cell,first_direction),(second_cell,second_direction)=entries
        same_direction=first_direction==second_direction
        push!(adjacency[first_cell],(second_cell,same_direction))
        push!(adjacency[second_cell],(first_cell,same_direction))
    end
    seen=falses(length(faces));flipped=falses(length(faces))
    seen[1]=true;stack=Int32[1];visited=0
    while !isempty(stack)
        cell=pop!(stack);visited+=1
        for (neighbor,same_direction) in adjacency[cell]
            required_flip=xor(flipped[cell],same_direction)
            if seen[neighbor]
                flipped[neighbor]==required_flip || throw(ArgumentError(
                    "$caller: Surface Loop[$shell] triangle winding is " *
                    "non-orientable"))
            else
                seen[neighbor]=true
                flipped[neighbor]=required_flip
                push!(stack,neighbor)
            end
        end
    end
    visited==length(faces) || throw(ArgumentError(
        "$caller: Surface Loop[$shell] triangle mesh is disconnected"))

    anchor=coordinates[faces[1][1]]
    signed_total=0.0;signed_correction=0.0
    magnitude_total=0.0;magnitude_correction=0.0
    for (cell,(first_node,second_node,third_node)) in pairs(faces)
        flipped[cell] && ((second_node,third_node)=(third_node,second_node))
        value=tet_signed_volume(
            anchor,coordinates[first_node],coordinates[second_node],
            coordinates[third_node])
        isfinite(value) || throw(ArgumentError(
            "$caller: Surface Loop[$shell] volume is not finite"))
        signed_total,signed_correction=_model_compensated_add(
            signed_total,signed_correction,value)
        magnitude_total,magnitude_correction=_model_compensated_add(
            magnitude_total,magnitude_correction,abs(value))
    end
    volume=abs(signed_total)
    volume>0 || throw(ArgumentError(
        "$caller: Surface Loop[$shell] encloses zero represented volume"))
    return volume,magnitude_total
end

function _model_mesh_volume(mesh::Mesh,caller::AbstractString)
    total=0.0;correction=0.0
    for cell in 1:ntets(mesh)
        value=tet_volume(
            _model_mesh_coordinate(mesh,mesh.tets[1,cell]),
            _model_mesh_coordinate(mesh,mesh.tets[2,cell]),
            _model_mesh_coordinate(mesh,mesh.tets[3,cell]),
            _model_mesh_coordinate(mesh,mesh.tets[4,cell]))
        isfinite(value) || throw(ArgumentError(
            "$caller: tetrahedron $cell has non-finite volume"))
        total,correction=_model_compensated_add(total,correction,value)
    end
    (isfinite(total) && total>0) || throw(ArgumentError(
        "$caller: tetrahedron mesh has no finite positive volume"))
    return total
end

function _model_certify_explicit_volume_semantics(
    mesh::Mesh,volume::Int,expected_volume::Float64,
    comparison_scale::Float64,caller::AbstractString)
    actual_volume=_model_mesh_volume(mesh,caller)
    tolerance=128eps(Float64)*max(comparison_scale,actual_volume)
    abs(actual_volume-expected_volume)<=tolerance || throw(ArgumentError(
        "$caller: Volume[$volume] surface loops do not form disjoint cavity " *
        "shells inside the exterior (mesh volume $actual_volume; " *
        "exterior-minus-cavities volume $expected_volume)"))
    return nothing
end

function _model_explicit_volume_fill(
    surface::Mesh,volume::Int,caller::AbstractString;
    interior_points=nothing)
    try
        return interior_points===nothing ? tetrahedralize(surface) :
               tetrahedralize(surface;interior_points=interior_points)
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        throw(ArgumentError(
            "$caller: explicit Volume[$volume] shell geometry cannot be " *
            "tetrahedralized — $(sprint(showerror,err))"))
    end
end

function _model_volume_periodic_surface_constraints(
    m::GeoModel,volume::Int,caller::AbstractString)
    _,curve_tags,surface_tags,_,_,boundary_surfaces=
        _model_volume_embedding_inventory(m,volume,caller)
    curve_set=Set(curve_tags)
    surface_set=Set(surface_tags)
    boundary_set=Set(abs.(boundary_surfaces))
    constraints=ModelPeriodicConstraint[]
    for constraint in model_periodic_constraints(m)
        slave=Int(constraint.slave_entity)
        master=Int(constraint.master_entity)
        if constraint.dim==1
            slave_present=slave in curve_set
            master_present=master in curve_set
            (slave_present || master_present) || continue
            slave_present==master_present || throw(ArgumentError(
                "$caller: periodic Curve[$slave]/Curve[$master] relation " *
                "has only one entity on Volume[$volume]"))
            throw(ArgumentError(
                "$caller: explicit Volume[$volume] supports planar periodic " *
                "surfaces but not independent periodic curves"))
        elseif constraint.dim==2
            slave_present=slave in surface_set
            master_present=master in surface_set
            (slave_present || master_present) || continue
            slave_present==master_present || throw(ArgumentError(
                "$caller: periodic Surface[$slave]/Surface[$master] relation " *
                "has only one entity on Volume[$volume]"))
            (slave in boundary_set && master in boundary_set) ||
                throw(ArgumentError(
                    "$caller: periodic Surface[$slave]/Surface[$master] must " *
                    "pair boundary surfaces of explicit Volume[$volume]"))
            push!(constraints,constraint)
        else
            throw(ArgumentError(
                "$caller: unsupported periodic dimension $(constraint.dim)"))
        end
    end
    return _model_periodic_constraint_order(constraints,caller)
end

function _model_explicit_volume_geometry(
    m::GeoModel,t::Int,caller::AbstractString)
    boundaries=_model_volume_boundary_surfaces(m,t,caller)
    isempty(boundaries) && throw(ArgumentError(
        "$caller: Volume[$t] has no explicit boundary surfaces"))
    local_meshes=Dict{Int,Mesh}()
    for signed_surface in boundaries
        surface=abs(signed_surface)
        get!(local_meshes,surface) do
            _model_planar_surface_mesh(
                m,surface,caller;include_embeddings=true)
        end
    end
    for constraint in _model_volume_periodic_surface_constraints(m,t,caller)
        slave=Int(constraint.slave_entity)
        master=Int(constraint.master_entity)
        local_meshes[slave]=_model_periodic_surface_mesh(
            m,local_meshes[master],constraint,caller)
    end
    coordinates=NTuple{3,Float64}[]
    node_index=Dict{NTuple{3,Float64},Int32}()
    triangles=NTuple{3,Int32}[]
    seen_faces=Set{NTuple{3,Int32}}()
    shell_faces=Vector{NTuple{3,Int32}}[]
    function global_node(coordinate)
        key=ntuple(axis->_model_projection_coordinate_key(coordinate[axis]),3)
        return get!(node_index,key) do
            length(coordinates)<typemax(Int32) || throw(ArgumentError(
                "$caller: explicit Volume[$t] boundary exceeds the Int32 node limit"))
            push!(coordinates,key)
            Int32(length(coordinates))
        end
    end
    for shell in m.volumes[t]
        faces=NTuple{3,Int32}[]
        for signed_surface in m.surface_loops[shell]
            surface=abs(signed_surface)
            local_mesh=local_meshes[surface]
            local_nodes=Vector{Int32}(undef,nnodes(local_mesh))
            for node in 1:nnodes(local_mesh)
                coordinate=(local_mesh.coords[1,node],local_mesh.coords[2,node],
                            local_mesh.coords[3,node])
                local_nodes[node]=global_node(coordinate)
            end
            for cell in 1:ntris(local_mesh)
                face=(local_nodes[local_mesh.tris[1,cell]],
                      local_nodes[local_mesh.tris[2,cell]],
                      local_nodes[local_mesh.tris[3,cell]])
                key=_model_projection_face_key(face)
                key in seen_faces && throw(ArgumentError(
                    "$caller: explicit Volume[$t] surfaces overlap on " *
                    "triangle $key"))
                push!(seen_faces,key)
                push!(triangles,face)
                push!(faces,face)
            end
        end
        push!(shell_faces,faces)
    end
    isempty(triangles) && throw(ArgumentError(
        "$caller: explicit Volume[$t] boundary has no triangles"))
    length(triangles)<=typemax(Int32) || throw(ArgumentError(
        "$caller: explicit Volume[$t] boundary exceeds the Int32 triangle limit"))
    coordinate_matrix=Matrix{Float64}(undef,3,length(coordinates))
    for (node,coordinate) in pairs(coordinates)
        coordinate_matrix[:,node].=coordinate
    end
    triangle_matrix=Matrix{Int32}(undef,3,length(triangles))
    for (cell,face) in pairs(triangles)
        triangle_matrix[:,cell].=face
    end
    output=Mesh(coordinate_matrix;tris=triangle_matrix)
    diagnostic=validate(output)
    diagnostic.ok || throw(ArgumentError(
        "$caller: explicit Volume[$t] boundary mesh is invalid — " *
        join(diagnostic.messages,"; ")))
    shell_volumes=Float64[]
    comparison_scale=0.0
    for (index,faces) in pairs(shell_faces)
        shell=m.volumes[t][index]
        volume,scale=_model_oriented_shell_volume(
            coordinates,faces,shell,caller)
        push!(shell_volumes,volume)
        comparison_scale+=scale
    end
    isfinite(comparison_scale) || throw(ArgumentError(
        "$caller: Volume[$t] shell-volume scale is not finite"))
    cavity_volume=0.0;cavity_correction=0.0
    for volume in Iterators.drop(shell_volumes,1)
        cavity_volume,cavity_correction=_model_compensated_add(
            cavity_volume,cavity_correction,volume)
    end
    expected_volume=first(shell_volumes)-cavity_volume
    expected_volume>0 || throw(ArgumentError(
        "$caller: Volume[$t] exterior volume $(first(shell_volumes)) is not " *
        "larger than its total cavity volume " *
        "$cavity_volume"))
    return (surface=output,expected_volume=expected_volume,
            comparison_scale=comparison_scale)
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
        haskey(m.boolean_operands,t) || throw(ArgumentError(
            "$caller: Boolean Volume[$t] has no owned operand snapshot; " *
            "recreate the Boolean in a fresh GeoModel"))
        A,B=m.boolean_operands[t]
        return mesh_boolean(A,B,spec.op)
    elseif !isempty(m.volumes[t])
        geometry=_model_explicit_volume_geometry(m,t,caller)
        probe=_model_explicit_volume_fill(geometry.surface,t,caller)
        _model_certify_explicit_volume_semantics(
            probe,t,geometry.expected_volume,geometry.comparison_scale,caller)
        return geometry.surface
    end
    throw(ArgumentError("$caller: Volume[$t] has no native solid encoding"))
end

"""
    mesh_model_volume(model, tag) -> Mesh

Mesh a native primitive, Boolean, or explicitly modeled planar-surface volume,
recovering supported embedded points, curves, and planar sheets with optional holes,
including nested points and curves constrained to those sheets. The returned
tetrahedral mesh is validated before it is returned. For an explicit volume, stored
planar boundary-surface relations synchronize the slave facets from their masters and
certify their affine tetrahedron-boundary node maps. Unsupported solid encodings
raise an explicit error.
"""
function mesh_model_volume(m::GeoModel, tag::Integer)
    caller="mesh_model_volume"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    explicit_geometry=isempty(m.volumes[t]) ? nothing :
        _model_explicit_volume_geometry(m,t,caller)
    periodic_surfaces=explicit_geometry===nothing ? ModelPeriodicConstraint[] :
        _model_volume_periodic_surface_constraints(m,t,caller)
    surface=explicit_geometry===nothing ? _volume_surface(m,t) :
                                         explicit_geometry.surface
    extra=NTuple{3,Float64}[]
    line_tags=Int[]
    sheets=Tuple{Int,NTuple{3,Float64},NTuple{3,Float64},
                       NTuple{3,Float64}}[]
    _,_,surface_tags,surface_embedded_points,surface_embedded_curves,
    boundary_surfaces=
        _model_volume_embedding_inventory(m,t,caller)
    boundary_surface_set=Set(abs.(boundary_surfaces))
    for (edim,etag) in get(m.embeds,(3,t),NTuple{2,Int}[])
        if edim==0
            haskey(m.points,etag) || throw(ArgumentError("$caller: unknown embedded Point[$etag]"))
            push!(extra, m.points[etag])
        elseif edim==1
            haskey(m.curves,etag) || throw(ArgumentError("$caller: unknown embedded Curve[$etag]"))
            a,b=m.curves[etag]
            push!(line_tags,etag)
            push!(extra, m.points[a]); push!(extra, m.points[b])
        elseif edim==2
            for point in surface_embedded_points[etag]
                push!(extra,m.points[point])
            end
            for curve in surface_embedded_curves[etag]
                push!(line_tags,curve)
                a,b=m.curves[curve]
                push!(extra,m.points[a]);push!(extra,m.points[b])
            end
            loops=m.surfaces[etag]
            if length(loops)==1
                ids=_loop_points(m,only(loops))
                length(ids)>=3 || throw(ArgumentError(
                    "$caller: embedded Surface[$etag] needs at least three points"))
                points=NTuple{3,Float64}[m.points[point] for point in ids]
                for index in 2:(length(points)-1)
                    push!(sheets,(etag,points[1],points[index],points[index+1]))
                end
            else
                sheet_mesh=_model_planar_surface_mesh(
                    m,etag,caller;include_embeddings=true)
                for cell in 1:ntris(sheet_mesh)
                    points=ntuple(slot->_model_mesh_coordinate(
                        sheet_mesh,sheet_mesh.tris[slot,cell]),3)
                    push!(sheets,(etag,points...))
                end
            end
        else
            throw(ArgumentError("$caller: unsupported embedding dimension $edim"))
        end
    end
    mesh=if explicit_geometry===nothing
        isempty(extra) ? tetrahedralize(surface) :
                         tetrahedralize(surface;interior_points=extra)
    else
        _model_explicit_volume_fill(
            surface,t,caller;
            interior_points=isempty(extra) ? nothing : extra)
    end
    sort!(unique!(line_tags))
    for curve in line_tags
        a,b=m.curves[curve];p=m.points[a];q=m.points[b]
        mesh=recover_segment3(mesh,p,q)
        mesh_covers_segment3(mesh,p,q) || throw(ErrorException(
            "$caller: embedded Curve[$curve] is not a chain of tetrahedron edges"))
    end
    for (sheet,a,b,c) in sheets
        mesh=recover_triangle3(mesh,a,b,c)
        mesh_covers_triangle3(mesh,a,b,c) || throw(ErrorException(
            "$caller: embedded Surface[$sheet] is not a union of tetrahedron faces"))
    end
    for curve in line_tags
        start_point,stop_point=m.curves[curve]
        mesh_covers_segment3(
            mesh,m.points[start_point],m.points[stop_point]) ||
            throw(ErrorException(
                "$caller: embedded Curve[$curve] is absent from the final " *
                "tetrahedron edge complex"))
    end
    for (sheet,a,b,c) in sheets
        mesh_covers_triangle3(mesh,a,b,c) || throw(ErrorException(
            "$caller: embedded Surface[$sheet] is absent from the final " *
            "tetrahedron face complex"))
    end
    for constraint in periodic_surfaces
        mapping=_model_periodic_surface_nodes(m,mesh,constraint)
        _model_mapping_matches(mesh,constraint,mapping,caller;exact=false)
    end
    tet_edges=_tet_edge_set(mesh)
    point_nodes=Dict{Int,Int32}()
    curve_entries=Dict{Int,Vector{Tuple{Float64,Int}}}()
    final_boundary_faces=Set(first(boundary_faces(mesh.tets)))
    for surface_tag in surface_tags
        nested_points=surface_embedded_points[surface_tag]
        nested_curves=surface_embedded_curves[surface_tag]
        isempty(nested_points) && isempty(nested_curves) && continue
        faces=if surface_tag in boundary_surface_set
            _model_projection_boundary_surface_faces!(
                Set{NTuple{3,Int32}}(),final_boundary_faces,
                m,mesh,surface_tag,caller)
        else
            _model_projection_volume_surface_faces!(
                Set{NTuple{3,Int32}}(),m,mesh,surface_tag,caller)
        end
        nodes,edges=_model_projection_face_topology(faces)
        for point in nested_points
            point_nodes[point]=_model_projection_embedded_point_node(
                m,mesh,point,1e-12,caller)
        end
        for curve in nested_curves
            curve_entries[curve]=_model_projection_tet_curve_nodes(
                m,mesh,curve,tet_edges,caller)
        end
        _model_projection_validate_nested_surface(
            surface_tag,nodes,edges,point_nodes,curve_entries,
            nested_points,nested_curves,caller)
    end
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntets(mesh)>0 || throw(ErrorException("$caller: Volume[$t] produced no tetrahedra"))
    explicit_geometry===nothing || _model_certify_explicit_volume_semantics(
        mesh,t,explicit_geometry.expected_volume,
        explicit_geometry.comparison_scale,caller)
    return mesh
end

end # module
