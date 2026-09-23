"""
    Model

A native geometry/entity kernel: tagged and named points, curves, curve loops,
surfaces, surface loops, and volumes with atomic retagging, dependency-safe
removal, analytical spatial queries, physical groups, and classified
surface/volume mixed-mesh projection. Native entity type and partition-ownership
metadata, explicit Point/Line/Plane evaluation, presentation state, Point-coordinate
updates, and model attributes are available without meshing.
Meshing dispatches to Tessella's certified simplex and transfinite kernels. This
is not OpenCASCADE; unsupported CAD statements remain explicit blockers.
"""
module Model

using ..MeshTypes: Mesh, validate, nnodes, nsegs, ntris, ntets, boundary_faces,
                   triangle_area, tet_signed_volume, tet_volume
using ..Elements: ElementBlock, MixedEntity, MixedEntityData,
                  MixedPeriodicLink, MixedMesh, msh_num_nodes, msh_family,
                  msh_dimension
using ..Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using ..SizeField: AbstractSizeField, ConstantSize, FunctionSize, MinSize,
                   PostViewField, field_value, size_at
using ..Geometry: box_surface, cylinder_surface, sphere_surface, cone_surface
using ..Mesh3D: tetrahedralize, mesh_boolean, recover_segment3, recover_triangle3,
                refine_to_size
using ..Mesh3D: mesh_covers_segment3, mesh_covers_triangle3,
                _tet_edge_set, _mesh_covering_faces3, _certify_surface_fill
using ..Periodic: periodic_identify_affine
using ..TransfiniteVolume: mesh_transfinite_volume
using ..TransfiniteCurve: transfinite_curve_parameters, transfinite_curve_hwall
using ..TransfiniteTriangle: mesh_transfinite_triangle,
                             mesh_transfinite_triangle_collapsed
using ..Transfinite: mesh_transfinite_patch
using ..Transform: _affine_coordinate, _transform_homogeneous,
    _periodic_affine_3x4, _periodic_affine_input
using ..Predicates: orient2, orient3
using ..GmshLibm: _gm_sin, _gm_cos, _gm_tan, _gm_asin, _gm_acos, _gm_atan,
                  _gm_atan2, _gm_pow, _gm_sincos
import LinearAlgebra
using LinearAlgebra: Symmetric, eigen
using Printf: @sprintf

export GeoModel, add_point!, set_point_mesh_size!, add_line!, add_curve_loop!, add_plane_surface!
export add_circle_arc!, add_ellipse_arc!, add_ruled_surface!
export add_spline!, add_bspline!, add_bezier!, add_nurbs!
export add_surface_loop!, add_volume!
export add_box!, add_cylinder!, add_sphere!, add_cone!, add_torus!, boolean_volumes!, boolean_volumes_multi!
export embed!, translate_volume!, dilate_volume!, rotate_volume!
export ModelPeriodicConstraint, set_periodic!, model_periodic_constraints,
       model_periodic_nodes, model_to_mixed
export add_physical_group!, set_physical_name!, remove_physical_groups!
export remove_physical_name!, model_physical_groups
export model_physical_groups_entities, model_entities_for_physical_group
export model_entities_for_physical_name, model_physical_groups_for_entity
export model_physical_name
export model_entities, model_dimension, model_boundary, model_adjacencies
export model_bounding_box, model_entities_in_bounding_box
export model_entity_type, model_entity_properties, model_parent,
       model_number_of_partitions, model_partitions
export model_value, model_derivative, model_second_derivative, model_curvature
export model_principal_curvatures, model_normal, model_parametrization
export model_parametrization_bounds, model_is_inside, model_closest_point
export model_reparametrize_on_surface
export set_entity_name!, remove_entity_name!, model_entity_name, model_set_tag!
export remove_entities!
export set_entity_visibility!, model_entity_visibility
export set_entity_color!, model_entity_color, set_point_coordinates!
export set_model_attribute!, model_attribute, model_attribute_names
export remove_model_attribute!
export mesh_model_surface, mesh_model_volume, model_entity, model_physical_tags

"""
Owned affine relation between two native entities. `affine` maps the master
entity to the slave entity in Gmsh row-major 4×4 order; it is `nothing` for a
dimension-1 relation declared without a transform (Gmsh's orientation-only
`Periodic Curve {slave} = {master}` form, serialized as a zero-row MSH affine
record). For dimension 1, `reversed` records whether the master start maps to
the slave end — inferred from the signed tags for transform-free relations and
from endpoint geometry otherwise; it is false for dimension 2.
"""
struct ModelPeriodicConstraint
    dim::Int
    slave_entity::Int32
    master_entity::Int32
    affine::Union{Nothing,NTuple{16,Float64}}
    reversed::Bool
    atol::Float64
end

"""
Parametrization-free model entity created by `add_discrete_entity!` (Gmsh's
`model.addDiscreteEntity`). Nodes and elements added through the mesh API are
recorded here with their caller-assigned — possibly sparse — tags; `boundary`
stores the declared boundary `(dim, tag)` pairs. The same record also backs the
per-entity `attached` store holding nodes and elements `addNodes`/`addElements`
place on a native entity.
"""
mutable struct DiscreteEntity
    boundary::Vector{NTuple{2,Int}}
    node_tags::Vector{Int32}
    node_coords::Matrix{Float64}
    node_params::Matrix{Float64}
    element_types::Vector{Int32}
    element_tags::Vector{Int32}
    element_nodes::Vector{Vector{Int32}}
    # Parametric coordinates of nodes owned by *other* records with respect to
    # this entity's parametrization (e.g. a classified surface's boundary nodes
    # live on curve records but still carry surface (u,v) parameters after
    # `create_geometry!`).
    aux_params::Dict{Int32,Vector{Float64}}
end

DiscreteEntity() = DiscreteEntity(NTuple{2,Int}[], Int32[],
    Matrix{Float64}(undef,3,0), Matrix{Float64}(undef,0,0),
    Int32[], Int32[], Vector{Int32}[], Dict{Int32,Vector{Float64}}())

Base.:(==)(a::DiscreteEntity,b::DiscreteEntity)=
    all(n->getfield(a,n)==getfield(b,n),fieldnames(DiscreteEntity))

# `Layers`/`Recombine`/`ScaleLast`/`QuadTri*` parameters of a built-in
# `Extrude` shape list, mirrored from Gmsh's `ExtrudeParams::mesh` record and
# attached to every entity the extrusion creates. `layers`/`heights` are empty
# unless a `Layers` parameter ran; `quad_to_tri` is `:none`, `:add_verts`, or
# `:no_new_verts`.
const _GeoExtrudeParams=NamedTuple{
    (:layers,:heights,:scale_last,:recombine,:quad_to_tri,:recomb_laterals),
    Tuple{Vector{Int},Vector{Float64},Bool,Bool,Symbol,Bool}}

"""
Owned per-entity meshing attributes on a [`GeoModel`](@ref), mirroring the
Gmsh `model.mesh` generation-attribute surface. Empty containers mean the
attribute is unset; `remove_constraints!` clears the generation-scoped
categories (transfinite, recombine, smoothing, reverse, algorithm, compound,
outward orientation) while periodic relations, embeddings, and Point sizes
live elsewhere on the model.
"""
mutable struct ModelMeshingAttributes
    transfinite_curves::Dict{Int,NamedTuple{(:num_nodes,:kind,:coef,:reversed),
                                           Tuple{Int,Symbol,Float64,Bool}}}
    transfinite_surfaces::Dict{Int,NamedTuple{(:arrangement,:corners),
                                             Tuple{Symbol,Vector{Int}}}}
    transfinite_volumes::Dict{Int,Vector{Int}}
    recombine::Dict{Tuple{Int,Int},Float64}
    extrude::Dict{Tuple{Int,Int},_GeoExtrudeParams}
    smoothing::Dict{Tuple{Int,Int},Int}
    reverse::Dict{Tuple{Int,Int},Bool}
    algorithm::Dict{Tuple{Int,Int},Int}
    size_at_params::Dict{Tuple{Int,Int},Vector{Tuple{Vector{Float64},Float64}}}
    size_from_boundary::Dict{Tuple{Int,Int},Bool}
    size_callback::Any
    compounds::Vector{Pair{Int,Vector{Int}}}
    outward_orientation::Set{Int}
    degenerated::Set{Int}
    quad_tri::Set{Int}
    order::Int
    transfinite_tri::Int
    attached::Dict{Tuple{Int,Int},DiscreteEntity}
    homology_requests::Vector{NamedTuple{(:kind,:domain,:subdomain,:dims),
        Tuple{String,Vector{Int},Vector{Int},Vector{Int}}}}
end

ModelMeshingAttributes() = ModelMeshingAttributes(
    Dict{Int,NamedTuple{(:num_nodes,:kind,:coef,:reversed),
                        Tuple{Int,Symbol,Float64,Bool}}}(),
    Dict{Int,NamedTuple{(:arrangement,:corners),Tuple{Symbol,Vector{Int}}}}(),
    Dict{Int,Vector{Int}}(),
    Dict{Tuple{Int,Int},Float64}(),
    Dict{Tuple{Int,Int},_GeoExtrudeParams}(),
    Dict{Tuple{Int,Int},Int}(),
    Dict{Tuple{Int,Int},Bool}(),
    Dict{Tuple{Int,Int},Int}(),
    Dict{Tuple{Int,Int},Vector{Tuple{Vector{Float64},Float64}}}(),
    Dict{Tuple{Int,Int},Bool}(),
    nothing,
    Pair{Int,Vector{Int}}[],
    Set{Int}(),
    Set{Int}(),
    Set{Int}(),
    1,
    0,
    Dict{Tuple{Int,Int},DiscreteEntity}(),
    NamedTuple{(:kind,:domain,:subdomain,:dims),
               Tuple{String,Vector{Int},Vector{Int},Vector{Int}}}[])

Base.:(==)(a::ModelMeshingAttributes,b::ModelMeshingAttributes)=
    all(n->getfield(a,n)==getfield(b,n),fieldnames(ModelMeshingAttributes))

mutable struct GeoModel
    points::Dict{Int,NTuple{3,Float64}}
    point_size::Dict{Int,Float64}
    curves::Dict{Int,NTuple{2,Int}}
    # Additional coincident vertices attached to a curve, mirroring the
    # control-point copies Gmsh's `Duplicata` creates for every duplicated
    # curve. Ordinary curves have no entry; transforms and coherence passes
    # treat them as part of the curve.
    curve_control_points::Dict{Int,Vector{Int}}
    # Native geometry kind per curve — :line (default), :circle, :ellipse, or
    # the spline family :spline/:bspline/:bezier/:nurbs. Arcs and splines carry
    # their ordered control points in `curve_control_points` ([start, center,
    # end] or [start, center, major, end] for arcs; the interpolation list for
    # splines); `curve_geometry` stores the `Plane{..}` hint normal on arcs and
    # the float32-rounded knot vector/degree/parameter interval on Nurbs,
    # mirroring Gmsh's `Curve` records.
    curve_types::Dict{Int,Symbol}
    curve_geometry::Dict{Int,NamedTuple}
    loops::Dict{Int,Vector{Int}}
    surfaces::Dict{Int,Vector{Int}}
    # Native geometry kind per surface — :plane (default), :ruled, or :tric.
    # `surface_geometry` carries auxiliary parameters such as the `In Sphere`
    # center point Gmsh stores on surface-filling records.
    surface_types::Dict{Int,Symbol}
    surface_geometry::Dict{Int,NamedTuple}
    surface_loops::Dict{Int,Vector{Int}}
    volumes::Dict{Int,Vector{Int}}
    entity_names::Dict{Tuple{Int,Int},String}
    entity_visibility::Dict{Tuple{Int,Int},Int}
    entity_colors::Dict{Tuple{Int,Int},NTuple{4,Int}}
    attributes::Dict{String,Vector{String}}
    physical::Dict{Tuple{Int,Int},Vector{Int}}
    physical_names::Dict{Tuple{Int,Int},String}
    # `GEntity::physicals` — the signed, insertion-ordered physical tags each
    # entity carries. `.geo` sync rebuilds it (`gmsh_sign(member) * Num`),
    # `add_physical_group!` appends `(t>0 ? tag : -tag)`, and the observable
    # `m.physical` group view derives from its `abs` values.
    entity_physicals::Dict{Tuple{Int,Int},Vector{Int}}
    physical_tag_max::Int
    box_extents::Dict{Int,NTuple{6,Float64}}
    cylinders::Dict{Int,NamedTuple{(:center,:axis,:radius,:height),
                                   Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64}}}
    spheres::Dict{Int,NamedTuple{(:center,:radius),Tuple{NTuple{3,Float64},Float64}}}
    cones::Dict{Int,NamedTuple{(:center,:axis,:r1,:r2,:height),
                               Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64,Float64}}}
    # binary Booleans store (op,a,b); multi-operand results store
    # (op,operands::Vector{Int})
    booleans::Dict{Int,NamedTuple}
    # binary Booleans store (A,B) snapshots; multi-operand results store
    # (meshes::Vector{Mesh}, cell::Union{UInt64,Nothing}) — see
    # `_boolean_result_surface`
    boolean_operands::Dict{Int,Any}
    # extra volume tag → primary tag of a multi-component Boolean result
    boolean_components::Dict{Int,Int}
    periodic::Dict{Tuple{Int,Int},ModelPeriodicConstraint}
    embeds::Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}
    meshing::ModelMeshingAttributes
    discrete::Dict{Tuple{Int,Int},DiscreteEntity}
    next_tag::Vector{Int}
    # `GModel::getName`/`setName` — set by `.geo` `SetName`, reset by
    # `Delete All`/`NewModel`.
    name::String
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
                      Dict{Int,Symbol}(), Dict{Int,NamedTuple}(),
                      Dict{Int,Vector{Int}}(), Dict{Int,Vector{Int}}(),
                      Dict{Int,Symbol}(), Dict{Int,NamedTuple}(),
                      Dict{Int,Vector{Int}}(), Dict{Int,Vector{Int}}(),
                      Dict{Tuple{Int,Int},String}(),
                      Dict{Tuple{Int,Int},Int}(),
                      Dict{Tuple{Int,Int},NTuple{4,Int}}(),
                      Dict{String,Vector{String}}(),
                      Dict{Tuple{Int,Int},Vector{Int}}(),
                      Dict{Tuple{Int,Int},String}(),
                      Dict{Tuple{Int,Int},Vector{Int}}(),
                      0,
                      Dict{Int,NTuple{6,Float64}}(),
                      Dict{Int,NamedTuple{(:center,:axis,:radius,:height),
                           Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64}}}(),
                      Dict{Int,NamedTuple{(:center,:radius),Tuple{NTuple{3,Float64},Float64}}}(),
                      Dict{Int,NamedTuple{(:center,:axis,:r1,:r2,:height),
                           Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64,Float64,Float64}}}(),
                      Dict{Int,NamedTuple{(:op,:a,:b),Tuple{Symbol,Int,Int}}}(),
                      Dict{Int,Tuple{Mesh,Mesh}}(),
                      Dict{Int,Int}(),
                      Dict{Tuple{Int,Int},ModelPeriodicConstraint}(),
                      Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}(),
                      ModelMeshingAttributes(),
                      Dict{Tuple{Int,Int},DiscreteEntity}(),
                      Int[0,0,0,0],
                      "")

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
    (-typemax(Int32)<=value<=typemax(Int32)) || throw(ArgumentError(
        "$caller: curve tag magnitude exceeds Int32"))
    return Int(value)
end

function _signed_surface_tag(value, caller)
    value isa Integer || throw(ArgumentError("$caller: surface tag must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: surface tag must not be Bool"))
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

# `requested == 0` is the public auto-allocation sentinel. `.geo` execution
# additionally needs Gmsh's `addX(int &tag)` convention — `tag < 0` auto-assigns
# while `tag == 0` stays literal — so `literal_zero` inserts at tag 0.
function _alloc_tag!(m::GeoModel, dim::Int, requested::Int, caller;
                     literal_zero::Bool=false)
    if requested==0 && !literal_zero
        # `getMaxTag(dim) + 1` is a C++ `int` increment — a bump past
        # INT32_MAX wraps to INT32_MIN, and `setMaxTag(dim, max(cur, num))`
        # then keeps the larger previous maximum (a wrapped counter stays
        # pinned at INT32_MAX).
        candidate=Int(reinterpret(Int32,(Int64(m.next_tag[dim+1])+1) % UInt32))
        while haskey(m.discrete,(dim,candidate))
            candidate=Int(reinterpret(Int32,(Int64(candidate)+1) % UInt32))
        end
        m.next_tag[dim+1]=max(m.next_tag[dim+1],candidate)
        return candidate
    end
    m.next_tag[dim+1]=max(m.next_tag[dim+1], requested)
    return requested
end

# `GModel::getMaxPhysicalNumber(dim)` — the largest `abs` physical tag carried
# by a dimension-`dim` entity (the synchronized view), or 0 when none exists.
function _max_entity_physical_number(m::GeoModel, dimension::Int)
    found=0
    for ((d,_),pnums) in m.entity_physicals
        d==dimension || continue
        for pnum in pnums
            a=abs(pnum)
            a>found && (found=a)
        end
    end
    return found
end

function _alloc_physical_tag(m::GeoModel, requested::Int, caller;
                             literal_tag::Bool=false, dimension::Int=-1)
    # `.geo` `Physical X(FExpr)` takes a raw `(int)` tag — 0 and negatives are
    # literal (`GEO_Internals::modifyPhysicalGroup` auto-assigns only for the
    # `"name"`-only form, which callers resolve before reaching here).
    literal_tag && return requested
    requested!=0 && return requested
    # `gmsh::model::addPhysicalGroup` auto-assigns
    # `max(getMaxPhysicalNumber(dim), getMaxPhysicalTag()) + 1` — the larger of
    # the synced per-dimension view and the stored global counter.
    current=if dimension<0
        m.physical_tag_max
    else
        max(m.physical_tag_max,_max_entity_physical_number(m,dimension))
    end
    current<typemax(Int32) || throw(ArgumentError(
        "$caller: no automatic physical tags remain"))
    return current+1
end

function _alloc_surface_loop_tag(m::GeoModel,requested::Int,caller;
                                 literal_zero::Bool=false)
    (requested!=0 || literal_zero) && return requested
    current=isempty(m.surface_loops) ? 0 : maximum(keys(m.surface_loops))
    current<typemax(Int32) || throw(ArgumentError(
        "$caller: no automatic Surface Loop tags remain"))
    return current+1
end

# Curve-loop tags live in their own Gmsh namespace (`_maxLineLoopNum`), not the
# curve counter — `Curve Loop() = {1}` after `Circle(7)` still yields loop 1.
function _alloc_curve_loop_tag(m::GeoModel,requested::Int,caller;
                               literal_zero::Bool=false)
    (requested!=0 || literal_zero) && return requested
    current=isempty(m.loops) ? 0 : maximum(keys(m.loops))
    current<typemax(Int32) || throw(ArgumentError(
        "$caller: no automatic Curve Loop tags remain"))
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
function add_point!(m::GeoModel, x, y, z; tag::Integer=0, mesh_size::Real=1.0,
                    _zero_literal::Bool=false)
    caller="add_point!"
    p=_finite3(x,y,z,caller)
    h=try Float64(mesh_size) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: mesh_size must be Float64-representable"))
    end
    (isfinite(h) && h>0) || throw(ArgumentError("$caller: mesh_size must be positive"))
    t=_alloc_tag!(m,0,_tag(tag,caller,0),caller;literal_zero=_zero_literal)
    (haskey(m.points,t) || haskey(m.discrete,(0,t))) && throw(ArgumentError("$caller: Point[$t] already exists"))
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
function add_line!(m::GeoModel, a, b; tag::Integer=0, _zero_literal::Bool=false)
    caller="add_line!"
    ta=_tag(a,caller,1); tb=_tag(b,caller,1)
    ta==tb && throw(ArgumentError("$caller: line endpoints must be distinct"))
    haskey(m.points,ta) || throw(ArgumentError("$caller: unknown Point[$ta]"))
    haskey(m.points,tb) || throw(ArgumentError("$caller: unknown Point[$tb]"))
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller;literal_zero=_zero_literal)
    (haskey(m.curves,t) || haskey(m.discrete,(1,t))) && throw(ArgumentError("$caller: Curve[$t] already exists"))
    m.curves[t]=(ta,tb)
    return t
end

# Periodic slave/master lists may carry signs: for a transform-free curve
# relation Gmsh derives the orientation from `sign(slave*master)`, and every
# lookup resolves `abs(tag)` (`addPeriodicEdge`/`addPeriodicFace`).
function _periodic_signed_tag(value,caller::AbstractString,name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$caller: $name tag must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: $name tag must not be Bool"))
    (-typemax(Int32)<=value<=typemax(Int32)) || throw(ArgumentError(
        "$caller: $name tag magnitude exceeds Int32"))
    return Int(value)
end

function _periodic_entity_tags(values,dim::Int,caller::AbstractString,
                               name::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: $name entities must be a vector or tuple"))
    return Int[_periodic_signed_tag(value,caller,name) for value in values]
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
    _model_require_line_curve(m,curve,caller,"curve length")
    a,b=m.curves[curve];p=m.points[a];q=m.points[b]
    length1=hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3])
    (isfinite(length1) && length1>0) || throw(ArgumentError(
        "$caller: Curve[$curve] must have finite positive geometric length"))
    return length1
end

include("ModelTopologyQueries.jl")
include("ModelIdentity.jl")
include("ModelRemoval.jl")
include("ModelEntityState.jl")
include("ModelSpatialQueries.jl")
include("ModelEntityMetadata.jl")
include("ModelEntityEvaluation.jl")
include("ModelMeshingAttributes.jl")
include("ModelTransforms.jl")
include("ModelCurved.jl")
include("ModelSplines.jl")
include("ModelBoolean.jl")

@inline function _model_periodic_entity_label(dim::Int)
    dim==1 && return "Curve"
    dim==2 && return "Surface"
    dim==3 && return "Volume"
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
    coefficients,translation=_periodic_affine_3x4(affine)
    point_map=Dict{Int,Int}()
    mapped_slaves=Set{Int}()
    for (index,master_point) in pairs(master_points)
        expected=_model_affine_point(
            coefficients,translation,m.points[master_point],caller,index)
        # `GFace::setMeshMaster` scans the slave vertex set in tag order and
        # keeps the first point inside `geom.tolerance*lc` — coincident slave
        # points are not an error until the bijection check below.
        slave_point=0
        dist_min=Inf
        for candidate in slave_points
            distance=_model_point_distance(m.points[candidate],expected)
            dist_min=min(dist_min,distance)
            if distance<atol
                slave_point=candidate
                break
            end
        end
        slave_point==0 && throw(ArgumentError(
            "$caller: no corresponding point $master_point for periodic " *
            "connection of surface $master to $slave (min. distance = " *
            "$dist_min, tolerance = $atol)"))
        point_map[master_point]=slave_point
        push!(mapped_slaves,slave_point)
    end
    length(mapped_slaves)==length(master_points) || throw(ArgumentError(
        "$caller: could not find all point correspondences for the " *
        "periodic connection from surface $master to $slave"))
    _model_periodic_surface_topology(
        m,slave,master,point_map,caller;
        include_embeddings=include_embeddings)
    return point_map
end

# `GFace::setMeshMaster(master, tfo)`'s induced-edge resolution: every slave
# boundary edge (`l_edges`, in loop order) and embedded edge must resolve a
# master counterpart through the vertex correspondence. Directed endpoint
# signatures give the unique forward/backward match; several candidates
# disambiguate by comparing the transformed parametric midpoint inside
# `localbb.diag()*1e-3`, then the transformed sampled bounding box — the same
# order and tolerances as upstream (`GEdge::bounds(true)` samples three
# parameter points, `point(0.5*(low+high))` is the mid). `on_match(slave,
# master)` fires per resolved edge, in upstream's multimap order — slave
# edges sorted by their `(begin, end)` vertex tags, boundary before embedded
# among equal keys — so the executor can emit Gmsh's `Setting curve master`
# progress lines. Returns `(point_map, pairs)` where `pairs` lists
# `slave_curve => master_curve` in that order.
function _model_periodic_surface_edge_map(
        m::GeoModel,slave::Int,master::Int,affine,atol::Float64,
        caller::AbstractString;on_match=nothing)
    function edge_order(surface::Int)
        curves=Int[];seen=Set{Int}()
        for loop in m.surfaces[surface],signed in m.loops[loop]
            curve=abs(signed)
            curve in seen || (push!(seen,curve);push!(curves,curve))
        end
        _,embedded=_model_surface_embedding_tags(m,surface,caller)
        return vcat(curves,embedded)
    end
    slave_edges=edge_order(slave);master_edges=edge_order(master)
    length(slave_edges)==length(master_edges) || throw(ArgumentError(
        "$caller: different number of curves ($(length(slave_edges)) vs " *
        "$(length(master_edges))) for periodic correspondence between " *
        "surfaces $master and $slave"))
    point_map=_model_periodic_surface_point_map(
        m,slave,master,affine,atol,caller;include_embeddings=true)
    slave_to_master=Dict{Int,Int}(slave_point=>master_point
        for (master_point,slave_point) in point_map)
    master_lookup=Dict{NTuple{2,Int},Vector{Int}}()
    for curve in master_edges
        endpoints=m.curves[curve]
        push!(get!(Vector{Int},master_lookup,endpoints),curve)
    end
    coefficients,translation=_periodic_affine_3x4(affine)
    # `localbb`/`localp` are computed once per slave edge upstream — lazily
    # here since the unique-match path never consults them.
    function slave_geometry(edge::Int)
        midpoint=try
            lower,upper=model_parametrization_bounds(m,1,edge)
            _model_curve_point(
                m,edge,0.5*(lower[1]+upper[1]),caller)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        return midpoint
    end
    function edge_bbox(edge::Int)
        lower,upper=model_parametrization_bounds(m,1,edge)
        corners=[_model_curve_point(
            m,edge,lower[1]+(index-1)/2*(upper[1]-lower[1]),caller)
                 for index in 1:3]
        lo=ntuple(axis->minimum(point->point[axis],corners),3)
        hi=ntuple(axis->maximum(point->point[axis],corners),3)
        return lo,hi
    end
    pairs=Pair{Int,Int}[]
    # The multimap orders by (begin, end) vertex tags; equal keys keep
    # insertion order (boundary first, then embedded).
    ordered_slave=sort(slave_edges;by=curve->m.curves[curve])
    for curve in ordered_slave
        s0,s1=m.curves[curve]
        mb0=slave_to_master[s0];mb1=slave_to_master[s1]
        forward_list=get(master_lookup,(mb0,mb1),Int[])
        backward_list=get(master_lookup,(mb1,mb0),Int[])
        master_curve=nothing
        if length(forward_list)==1 && (isempty(backward_list) || mb0==mb1)
            master_curve=forward_list[1]
        elseif length(backward_list)==1 && (isempty(forward_list) || mb0==mb1)
            master_curve=backward_list[1]
        else
            local_lo,local_hi=edge_bbox(curve)
            tolerance=1e-3*_model_point_distance(local_lo,local_hi)
            local_mid=slave_geometry(curve)
            for candidates in (forward_list,backward_list)
                master_curve===nothing || break
                isempty(candidates) && continue
                for candidate in candidates
                    master_mid=try
                        lower,upper=model_parametrization_bounds(m,1,candidate)
                        point=_model_curve_point(
                            m,candidate,0.5*(lower[1]+upper[1]),caller)
                        _model_affine_point(
                            coefficients,translation,point,caller,1)
                    catch err
                        err isa InterruptException && rethrow()
                        nothing
                    end
                    if local_mid!==nothing && master_mid!==nothing &&
                            _model_point_distance(local_mid,master_mid)<tolerance
                        master_curve=candidate
                        break
                    end
                    lo,hi=edge_bbox(candidate)
                    mapped_lo=_model_affine_point(
                        coefficients,translation,lo,caller,1)
                    mapped_hi=_model_affine_point(
                        coefficients,translation,hi,caller,1)
                    if _model_point_distance(mapped_lo,local_lo)<tolerance &&
                            _model_point_distance(mapped_hi,local_hi)<tolerance
                        master_curve=candidate
                        break
                    end
                end
            end
        end
        master_curve===nothing && throw(ArgumentError(
            "$caller: could not find counterpart of curve $curve with end " *
            "points $s0 $s1 (corresponding to curve with end points " *
            "$mb0 $mb1) in surface $master"))
        on_match===nothing || on_match(curve,master_curve)
        push!(pairs,curve=>master_curve)
    end
    return point_map,pairs
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

# `GEdge::setMeshMaster(ge, tfo)`'s endpoint test: transform the master
# endpoints by the affine and compare them with the slave endpoints both
# forward and reversed. Returns `(reversed, mismatch)` — `reversed` is
# `nothing` when neither assignment fits inside `tolerance` (upstream drops the
# relation with an Info diagnostic in that case).
function _model_periodic_curve_orientation(
        m::GeoModel,slave::Int,master::Int,coefficients,translation,
        tolerance::Float64,caller::AbstractString,index::Int=1)
    slave_points=m.curves[slave];master_points=m.curves[master]
    slave_start=m.points[slave_points[1]]
    slave_stop=m.points[slave_points[2]]
    mapped_start=_model_affine_point(
        coefficients,translation,m.points[master_points[1]],caller,index)
    mapped_stop=_model_affine_point(
        coefficients,translation,m.points[master_points[2]],caller,index+1)
    # `GEdge::setMeshMaster`'s distance names: `d01` is slave start against
    # the transformed master end and `d10` slave end against the transformed
    # master start — the Info diagnostic reports them in that order.
    d00=_model_point_distance(mapped_start,slave_start)
    d11=_model_point_distance(mapped_stop,slave_stop)
    d01=_model_point_distance(mapped_stop,slave_start)
    d10=_model_point_distance(mapped_start,slave_stop)
    # `GEdge::setMeshMaster` tries the forward correspondence first, then the
    # reversed one; a symmetric transform resolving both directions always
    # stays forward. On mismatch its Info reports the endpoint distances of
    # the direction preferred by the `d00*d11 < d01*d10` product comparison.
    # Both comparisons are strict like upstream's `d.norm() < tol`.
    reversed=if d00<tolerance && d11<tolerance
        false
    elseif d01<tolerance && d10<tolerance
        true
    else
        nothing
    end
    distances=d00*d11<d01*d10 ? (d00,d11) : (d01,d10)
    return (reversed=reversed,mismatch=max(distances...),distances=distances)
end

"""
    set_periodic!(model, dim, slave_entities, master_entities, affine;
                  atol=1e-12)

Persist affine relations between equally sized lists of straight native curves
(`dim=1`), planar native surfaces (`dim=2`), or native volumes (`dim=3`).
`affine` maps each master entity to its slave in Gmsh row-major 4×4 order —
a 16-entry vector or 4×4 matrix is stored verbatim like
`GEntity::setMeshMaster` (only the first twelve entries are ever applied), a
12-entry vector pads to the canonical homogeneous row.
Entity tags may be signed — lookups resolve `abs(tag)` like Gmsh's
`addPeriodicEdge`/`addPeriodicFace`. For `dim=1`, `affine` may be `nothing`,
recording Gmsh's orientation-only `Periodic Curve {slave} = {master}` relation:
each pair's `reversed` flag is then `sign(slave_tag*master_tag) < 0` and no
geometric endpoint check runs, matching `GEdge::setMeshMaster(source, ori)`.
An entity may be the master of multiple relations or both a slave and a master
in an acyclic dependency chain; each slave has exactly one master. Cycles are
explicit blockers. Periodic curves must belong to the same planar surface when
meshed. Periodic surfaces require disjoint, affine-equivalent boundary point
and loop topology; embedded topology is rechecked when an explicit planar-shell
volume is meshed. Periodic volumes are stored and reported as a Tessella
extension — Gmsh exposes no volume periodicity through either its `.geo`
grammar or `gmsh.model.mesh.setPeriodic` — and constrain no interior mesh:
volume correspondence is carried entirely by the periodic boundary entities.
The update is atomic. `overwrite=true` mirrors `.geo` redeclaration — like
upstream's `setMeshMaster`, a repeated slave tag silently replaces its stored
relation instead of raising.
"""
function set_periodic!(m::GeoModel,dim,slave_entities,master_entities,affine;
                       atol=1e-12,overwrite::Bool=false)
    caller="set_periodic!"
    d=_dimension(dim,caller)
    d in (1,2,3) || throw(ArgumentError(
        "$caller: only Curve, Surface, and Volume periodicity " *
        "(dimensions 1, 2, and 3) are implemented"))
    affine===nothing && d!=1 && throw(ArgumentError(
        "$caller: only dimension-1 relations may omit the affine transform"))
    label=_model_periodic_entity_label(d)
    slaves=_periodic_entity_tags(slave_entities,d,caller,"slave")
    masters=_periodic_entity_tags(master_entities,d,caller,"master")
    length(slaves)==length(masters) || throw(ArgumentError(
        "$caller: slave and master entity counts differ"))
    isempty(slaves) && throw(ArgumentError(
        "$caller: need at least one slave/master $label pair"))
    abs_slaves=abs.(slaves)
    length(unique(abs_slaves))==length(abs_slaves) || throw(ArgumentError(
        "$caller: slave $label tags must be unique"))
    overlap=sort!(Int[slave for slave in abs_slaves
                      if haskey(m.periodic,(d,slave))])
    # `.geo` redeclaration replaces a slave's relation like upstream's
    # `setMeshMaster` — `overwrite` skips the uniqueness check and the write
    # phase below replaces the stored constraint atomically.
    (overwrite || isempty(overlap)) || throw(ArgumentError(
        "$caller: $label[$(first(overlap))] already has a periodic master"))
    tolerance=_model_periodic_tolerance(atol,caller)
    coefficients=translation=nothing
    stored_affine=nothing
    if affine!==nothing
        coefficients,translation,stored_affine=_periodic_affine_input(
            affine,caller;name="affine transform")
    end
    pending=ModelPeriodicConstraint[]
    for (pair_index,(signed_slave,signed_master)) in
            enumerate(zip(slaves,masters))
        slave=abs(signed_slave);master=abs(signed_master)
        entities=d==1 ? m.curves : d==2 ? m.surfaces : m.volumes
        haskey(entities,slave) || throw(ArgumentError(
            "$caller: unknown slave $label[$slave]"))
        haskey(entities,master) || throw(ArgumentError(
            "$caller: unknown master $label[$master]"))
        slave!=master || throw(ArgumentError(
            "$caller: slave and master $label tags must differ"))
        reversed=false
        if d==1
            if coefficients===nothing
                # Gmsh's orientation-only path (`GEdge::setMeshMaster(source,
                # ori)`): no geometric check at all — `masterOrientation` is
                # the sign of `iSource*iTarget` and parameter correspondence
                # does the rest at mesh time.
                reversed=signed_slave*signed_master<0
            else
                _model_curve_length(m,slave,caller)
                _model_curve_length(m,master,caller)
                slave_points=m.curves[slave];master_points=m.curves[master]
                isempty(intersect(Set(slave_points),Set(master_points))) ||
                    throw(ArgumentError(
                        "$caller: Curve[$slave] and Curve[$master] must " *
                        "have disjoint endpoints"))
                orientation=_model_periodic_curve_orientation(
                    m,slave,master,coefficients,translation,tolerance,caller,
                    2pair_index-1)
                orientation.reversed===nothing && throw(ArgumentError(
                    "$caller: affine map misses slave Curve[$slave] endpoints " *
                    "by $(orientation.mismatch)"))
                reversed=orientation.reversed
            end
        elseif d==2
            # Upstream's `GFace::setMeshMaster` resolves the full vertex and
            # (boundary + embedded) edge correspondence at declaration time
            # and stores the induced curve masters — a slave edge without a
            # counterpart aborts the whole relation.
            _model_periodic_surface_edge_map(
                m,slave,master,stored_affine,tolerance,caller)
        end
        push!(pending,ModelPeriodicConstraint(
            d,Int32(slave),Int32(master),stored_affine,reversed,tolerance))
    end
    existing=model_periodic_constraints(m)
    isempty(overlap) ||
        filter!(existing) do constraint
            !((constraint.dim,Int(constraint.slave_entity)) in
              Set((d,slave) for slave in overlap))
        end
    _model_periodic_dependency_parents(vcat(existing,pending),caller)
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

Add a curve loop from signed curve tags. Mirroring Gmsh's `Curve Loop`, the
tags are reordered into connectivity order — each curve's oriented end must
match the next curve's oriented start — so input order is irrelevant. A chain
that dead-ends before consuming every curve fails; closed single-curve loops
(periodic curves) and multi-subloop lists are accepted like Gmsh. Closure of
each boundary chain is certified at mesh time.
"""
function add_curve_loop!(m::GeoModel, curves; tag::Integer=0,
                         _zero_literal::Bool=false)
    caller="add_curve_loop!"
    ids=Int[_signed_curve_tag(c,caller) for c in curves]
    for id in ids
        haskey(m.curves,abs(id)) || throw(ArgumentError("$caller: unknown Curve[$(abs(id))]"))
    end
    ordered=_sort_curve_loop(m,ids,caller)
    t=_alloc_curve_loop_tag(m,_tag(tag,caller,1),caller;
                            literal_zero=_zero_literal)
    haskey(m.loops,t) && throw(ArgumentError("$caller: Loop[$t] already exists"))
    m.loops[t]=ordered
    return t
end

"""
    add_plane_surface!(model, loops; tag=0) -> tag

Add a planar surface bounded by existing curve loops. The first loop is the
outer boundary and subsequent loops are holes.
"""
function add_plane_surface!(m::GeoModel, loops; tag::Integer=0,
                            _zero_literal::Bool=false)
    caller="add_plane_surface!"
    ids=Int[_tag(ℓ,caller,2) for ℓ in loops]
    isempty(ids) && throw(ArgumentError("$caller: need an outer loop"))
    for id in ids
        haskey(m.loops,id) || throw(ArgumentError("$caller: unknown Loop[$id]"))
    end
    t=_alloc_tag!(m,2,_tag(tag,caller,2),caller;literal_zero=_zero_literal)
    (haskey(m.surfaces,t) || haskey(m.discrete,(2,t))) && throw(ArgumentError("$caller: Surface[$t] already exists"))
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
    # Two incidence patterns extend the planar-manifold rule for OCC-style
    # solids: a periodic face's seam curve occurs twice on its own surface in
    # opposite directions, and a degenerate pole/apex edge occurs once. Both
    # curve kinds only arise from materialized primitive construction.
    for curve in sort!(collect(keys(curve_incidence)))
        count=curve_incidence[curve]
        owners=curve_owners[curve]
        count==2 && length(owners)==2 && continue
        a,b=m.curves[curve]
        _curve_type(m,curve)==:degenerate && a==b && continue
        if count==2 && length(owners)==1 &&
                _surface_type(m,only(owners))!=:plane
            surface=only(owners)
            signs=Int[]
            for loop in m.surfaces[surface], signed in m.loops[loop]
                abs(signed)==curve && push!(signs,sign(signed))
            end
            sort!(signs)==[-1,1] && continue
        end
        throw(ArgumentError(
            "$caller: Curve[$curve] must occur once on each of two distinct " *
            "surfaces (found $count occurrences on $(length(owners)) surfaces)"))
    end
    adjacency=Dict(surface=>Set{Int}() for surface in unsigned)
    for owners in values(curve_owners)
        length(owners)==2 || continue
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
are rejected before the model is changed. Two OCC-style exceptions are allowed:
a periodic (non-plane) face may carry a seam curve twice in opposite
directions, and a degenerate pole/apex edge may occur once — both patterns only
arise from materialized primitive construction.
"""
function add_surface_loop!(m::GeoModel,surfaces;tag::Integer=0,
                           _zero_literal::Bool=false,
                           _skip_validation::Bool=false)
    caller="add_surface_loop!"
    ids=Int[_signed_surface_tag(surface,caller) for surface in surfaces]
    isempty(ids) && throw(ArgumentError("$caller: need at least one Surface"))
    # Gmsh's `addSurfaceLoop` stores the member list without a closure check —
    # open shells error later at `Volume` creation. `.geo` execution takes the
    # same deferred path; the public API keeps eager validation.
    _skip_validation || _validate_surface_loop(m,ids,caller)
    requested=_tag(tag,caller,2)
    t=_alloc_surface_loop_tag(m,requested,caller;literal_zero=_zero_literal)
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
function add_volume!(m::GeoModel,surface_loops;tag::Integer=0,
                     _zero_literal::Bool=false,
                     _skip_validation::Bool=false)
    caller="add_volume!"
    shells=Int[_tag(shell,caller,3) for shell in surface_loops]
    isempty(shells) && throw(ArgumentError(
        "$caller: need an exterior Surface Loop"))
    all(shell->shell>=0,shells) || throw(ArgumentError(
        "$caller: Surface Loop tags must be non-negative"))
    # Gmsh's `SetVolumeSurfaces` checks only that each shell and every member
    # surface exists — closure, disjointness, and shell uniqueness are not
    # verified until meshing. `.geo` execution takes the same deferred path.
    if _skip_validation
        for shell in shells
            haskey(m.surface_loops,shell) || throw(ArgumentError(
                "$caller: unknown Surface Loop[$shell]"))
            for signed_surface in m.surface_loops[shell]
                surface=abs(signed_surface)
                (haskey(m.surfaces,surface) ||
                 haskey(m.discrete,(2,surface))) || throw(ArgumentError(
                    "$caller: unknown Surface[$surface]"))
            end
        end
        requested=_tag(tag,caller,3)
        (requested!=0 || _zero_literal) &&
            (haskey(m.volumes,requested) || haskey(m.discrete,(3,requested))) &&
            throw(ArgumentError("$caller: Volume[$requested] already exists"))
        t=_alloc_tag!(m,3,requested,caller;literal_zero=_zero_literal)
        haskey(m.volumes,t) && throw(ArgumentError(
            "$caller: Volume[$t] already exists"))
        m.volumes[t]=shells
        return t
    end
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
    (requested!=0 || _zero_literal) &&
        (haskey(m.volumes,requested) || haskey(m.discrete,(3,requested))) &&
        throw(ArgumentError("$caller: Volume[$requested] already exists"))
    t=_alloc_tag!(m,3,requested,caller;literal_zero=_zero_literal)
    haskey(m.volumes,t) && throw(ArgumentError(
        "$caller: Volume[$t] already exists"))
    m.volumes[t]=shells
    return t
end

"""
    add_box!(model, x, y, z, dx, dy, dz; tag=0) -> tag

Add an axis-aligned box volume with finite origin `(x,y,z)` and positive
extents `(dx,dy,dz)`. Like Gmsh's `addBox`, the box is a real boundary
representation: it owns 8 corner Points, 12 edge Curves, 6 planar Surfaces and
one Surface Loop, all allocated from the shared per-dimension tag namespaces
and queryable through the entity, boundary, and adjacency APIs. Entity order
matches Gmsh 4.15.2's `addBox` numbering, and the shell stores Gmsh's
oriented-boundary signs `[-1, 2, -3, 4, -5, 6]` (each face's stored loop
normal points along its positive coordinate axis). The `box_extents` encoding
is retained alongside the topology for primitive transforms, bounds, and
Boolean operand snapshots; volume meshing itself flows through the explicit
planar-shell path.

Like Gmsh's OCC `addBox` corners, the materialized corner Points carry no
explicit mesh-size constraint; the boundary size field derives their sizes
from incident edges. Assign explicit sizes with `set_point_mesh_size!` or
`setSize` to refine.
"""
function add_box!(m::GeoModel, xmin, ymin, zmin, dx, dy, dz; tag::Integer=0,
                  _zero_literal::Bool=false)
    caller="add_box!"
    origin=_finite3(xmin,ymin,zmin,caller)
    d=_finite3(dx,dy,dz,caller)
    (d[1]>0 && d[2]>0 && d[3]>0) || throw(ArgumentError("$caller: extents must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller;literal_zero=_zero_literal)
    (haskey(m.volumes,t) || haskey(m.discrete,(3,t))) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    x0,y0,z0=origin
    x1=x0+d[1]; y1=y0+d[2]; z1=z0+d[3]
    created_points=Int[];created_curves=Int[]
    created_loops=Int[];created_surfaces=Int[]
    shell=0
    try
        # Gmsh addBox corner order: x outermost, z innermost with zmax first.
        # Each entity is registered for rollback before the next call so a
        # mid-construction failure cannot leak partially built state.
        for corner in ((x0,y0,z1),(x0,y0,z0),(x0,y1,z1),(x0,y1,z0),
                       (x1,y0,z1),(x1,y0,z0),(x1,y1,z1),(x1,y1,z0))
            push!(created_points,add_point!(m,corner...))
        end
        p=created_points
        # Like Gmsh's OCC `addBox` corners, materialized corners carry no
        # explicit mesh-size constraint: `_volume_boundary_size_field` then
        # derives their sizes from incident edges, so a plain box meshes
        # identically to the encoding-only primitive it replaces.
        for corner in p
            delete!(m.point_size,corner)
        end
        # Gmsh addBox edge order and directions.
        for edge in ((p[2],p[1]),(p[1],p[3]),(p[4],p[3]),(p[2],p[4]),
                     (p[6],p[5]),(p[5],p[7]),(p[8],p[7]),(p[6],p[8]),
                     (p[2],p[6]),(p[1],p[5]),(p[4],p[8]),(p[3],p[7]))
            push!(created_curves,add_line!(m,edge...))
        end
        c=created_curves
        # Face loops, each oriented so its stored normal points along the
        # positive coordinate axis (xmin/ymin/zmin faces read inward, matching
        # the natural orientation of Gmsh's OCC faces).
        loop_curve_signs=([c[4],c[3],-c[2],-c[1]],
                          [c[8],c[7],-c[6],-c[5]],
                          [c[1],c[10],-c[5],-c[9]],
                          [c[3],c[12],-c[7],-c[11]],
                          [c[9],c[8],-c[11],-c[4]],
                          [c[10],c[6],-c[12],-c[2]])
        for signed_curves in loop_curve_signs
            push!(created_loops,add_curve_loop!(m,signed_curves))
            push!(created_surfaces,
                  add_plane_surface!(m,[last(created_loops)]))
        end
        shell=add_surface_loop!(m,[-created_surfaces[1],created_surfaces[2],
                                   -created_surfaces[3],created_surfaces[4],
                                   -created_surfaces[5],created_surfaces[6]])
    catch
        shell!=0 && delete!(m.surface_loops,shell)
        for surface in created_surfaces
            delete!(m.surfaces,surface)
        end
        for loop in created_loops
            delete!(m.loops,loop)
        end
        for curve in created_curves
            delete!(m.curves,curve)
        end
        for point in created_points
            delete!(m.points,point);delete!(m.point_size,point)
        end
        rethrow()
    end
    m.volumes[t]=[shell]
    m.box_extents[t]=(origin[1],origin[2],origin[3],d[1],d[2],d[3])
    return t
end

"""
    add_cylinder!(model, x, y, z, dx, dy, dz, radius; tag=0) -> tag

Add a cylinder whose base center is `(x,y,z)` and whose finite, nonzero axis
vector is `(dx,dy,dz)`. `radius` must be finite and positive.

Like Gmsh's OpenCASCADE `addCylinder`, the solid is a real boundary
representation: two rim Points, two closed-circle edges and a seam Line, a
`Cylinder` lateral face plus `Plane` caps, and one Surface Loop, all allocated
from the shared per-dimension tag namespaces and queryable through the entity,
boundary, and adjacency APIs. Entity order, wire signs `[-top,-seam,bottom,
seam]`, and shell signs `[lateral,+top,-bottom]` match Gmsh 4.15.2's OCC
layout. The compact `cylinders` encoding is retained alongside the topology;
volume meshing keeps flowing through the native analytic tessellation rather
than the planar-shell path.
"""
function add_cylinder!(m::GeoModel, x, y, z, dx, dy, dz, radius; tag::Integer=0,
                       _zero_literal::Bool=false)
    caller="add_cylinder!"
    c=_finite3(x,y,z,caller); a=_finite3(dx,dy,dz,caller)
    h=_occ_modulus(a)
    (isfinite(h) && h>0) || throw(ArgumentError("$caller: axis must have finite positive length"))
    r=_finite_scalar(radius,caller,"radius")
    r>0 || throw(ArgumentError("$caller: radius must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller;literal_zero=_zero_literal)
    (haskey(m.volumes,t) || haskey(m.discrete,(3,t))) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    shell=_materialize_cylinder!(m,c,a,r,h)
    m.volumes[t]=[shell]
    m.cylinders[t]=(center=c, axis=a, radius=r, height=h)
    return t
end

"""
    add_sphere!(model, x, y, z, radius; tag=0) -> tag

Add a sphere with a finite center and finite positive radius.

Like Gmsh's OpenCASCADE `addSphere`, the solid materializes two pole Points,
a degenerate edge on each pole, a meridian Circle edge, and one `Sphere` face
bounded by `[-degN,-meridian,degS,meridian]` — queryable through the entity,
boundary, and adjacency APIs. The compact `spheres` encoding is retained for
native meshing and bounds.
"""
function add_sphere!(m::GeoModel, x, y, z, radius; tag::Integer=0,
                     _zero_literal::Bool=false)
    caller="add_sphere!"
    c=_finite3(x,y,z,caller)
    r=_finite_scalar(radius,caller,"radius")
    r>0 || throw(ArgumentError("$caller: radius must be positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller;literal_zero=_zero_literal)
    (haskey(m.volumes,t) || haskey(m.discrete,(3,t))) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    shell=_materialize_sphere!(m,c,r)
    m.volumes[t]=[shell]
    m.spheres[t]=(center=c, radius=r)
    return t
end

"""
    add_cone!(model, x, y, z, dx, dy, dz, r1, r2; tag=0) -> tag

Add a cone or conical frustum along the finite, nonzero axis `(dx,dy,dz)`.
Both radii must be finite and non-negative, and at least one must be positive.

Like Gmsh's OpenCASCADE `addCone`, the boundary representation matches the
cylinder layout except that a zero-radius end collapses to a degenerate apex
edge and loses its Plane cap: a frustum gets `[lateral,+top,-bottom]`, `r1=0`
gets `[lateral,+top]`, and `r2=0` gets `[lateral,-bottom]`. The compact `cones`
encoding is retained for native meshing and bounds.
"""
function add_cone!(m::GeoModel, x, y, z, dx, dy, dz, r1, r2; tag::Integer=0,
                   _zero_literal::Bool=false)
    caller="add_cone!"
    c=_finite3(x,y,z,caller); a=_finite3(dx,dy,dz,caller)
    h=_occ_modulus(a)
    (isfinite(h) && h>0) || throw(ArgumentError("$caller: axis must have finite positive length"))
    ra=_finite_scalar(r1,caller,"r1")
    rb=_finite_scalar(r2,caller,"r2")
    (ra>=0 && rb>=0 && (ra>0 || rb>0)) || throw(ArgumentError(
        "$caller: radii must be non-negative with at least one positive"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller;literal_zero=_zero_literal)
    (haskey(m.volumes,t) || haskey(m.discrete,(3,t))) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    shell=_materialize_cone!(m,c,a,ra,rb,h)
    m.volumes[t]=[shell]
    m.cones[t]=(center=c, axis=a, r1=ra, r2=rb, height=h)
    return t
end

"""
    add_torus!(model, x, y, z, r1, r2; tag=0, angle=2π) -> tag

Add a torus of major radius `r1` and minor radius `r2` centered at `(x,y,z)`,
revolved about the z axis through `angle` radians (`2π` for a full torus).
Both radii must be finite and positive and `angle` finite in `(0, 2π]` —
matching `gmsh.model.occ.addTorus` and the `.geo` `Torus` statement, which
build every torus about the z axis (apply `Rotate` for other axes).

Like Gmsh's OpenCASCADE result, the solid materializes real boundary topology:
a full torus is one rim Point, an outer-equator Circle and a meridian Circle
closed on it, and a single `Torus` face wired `[-equator,+meridian,+equator,
-meridian]`; a partial torus adds a second rim vertex, trims the equator to
`[0,angle]`, closes the two end meridians, and caps the ends with `Plane`
faces under shell `[torus,+start_cap,-end_cap]`. There is no compact encoding:
entity, boundary, evaluation, and bounds queries read the materialized
entities, while volume meshing rejects the non-planar face explicitly.
"""
function add_torus!(m::GeoModel, x, y, z, r1, r2;
                    tag::Integer=0, angle::Real=2π, _zero_literal::Bool=false)
    caller="add_torus!"
    c=_finite3(x,y,z,caller)
    ra=_finite_scalar(r1,caller,"r1")
    rb=_finite_scalar(r2,caller,"r2")
    a=_finite_scalar(angle,caller,"angle")
    (ra>0 && rb>0) || throw(ArgumentError(
        "$caller: radii must be positive"))
    (a>0 && a<=2π) || throw(ArgumentError(
        "$caller: angle must lie in (0, 2π]"))
    t=_alloc_tag!(m,3,_tag(tag,caller,3),caller;literal_zero=_zero_literal)
    (haskey(m.volumes,t) || haskey(m.discrete,(3,t))) && throw(ArgumentError("$caller: Volume[$t] already exists"))
    shell=_materialize_torus!(m,c,ra,rb,a)
    m.volumes[t]=[shell]
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

"""
    remove_embedded!(m::GeoModel, dim_tags, dim=-1)

Remove embedded entities recorded on each listed parent entity. `dim` below 0
removes every embedded dimension; otherwise only embedded entities of that
dimension are dropped. Unknown parents fail explicitly.
"""
function remove_embedded!(m::GeoModel, dim_tags, dim=-1)
    caller="remove_embedded!"
    filter_dim=_query_dimension(dim,caller)
    filter_dim in (-1,0,1,2) || throw(ArgumentError(
        "$caller: embedded dimension filter must be -1 or in 0:2"))
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        pd=_dimension(pair[1],caller)
        pt=_tag(pair[2],caller,pd)
        pd in (2,3) || throw(ArgumentError(
            "$caller: embedded parents must be Surfaces or Volumes " *
            "(got dim=$pd)"))
        parent_dict=pd==2 ? m.surfaces : m.volumes
        haskey(parent_dict,pt) || throw(ArgumentError(
            "$caller: unknown $(pd==2 ? "Surface" : "Volume")[$pt]"))
        entries=get(m.embeds,(pd,pt),NTuple{2,Int}[])
        if filter_dim<0
            delete!(m.embeds,(pd,pt))
        else
            kept=[child for child in entries if child[1]!=filter_dim]
            isempty(kept) ? delete!(m.embeds,(pd,pt)) :
                            (m.embeds[(pd,pt)]=kept)
        end
    end
    return nothing
end

# Every Point owned through the volume's surface loops — used to keep a
# materialized primitive's boundary entities (add_box!) synchronized with its
# compact encoding under the native transforms below.
function _model_volume_owned_points(m::GeoModel,t::Int)
    points=Set{Int}()
    for shell in m.volumes[t], signed_surface in m.surface_loops[shell],
        loop in m.surfaces[abs(signed_surface)], signed_curve in m.loops[loop]
        a,b=m.curves[abs(signed_curve)]
        push!(points,a);push!(points,b)
    end
    return points
end

"""
    translate_volume!(model, tag, offset) -> tag

Translate a native volume by the finite three-component `offset`. Explicit
boundary topology and materialized boundary Points (from `add_box!`) move with
the volume; compact primitive encodings are updated in place. The model is
unchanged if a translated coordinate is not representable as a finite `Float64`.
"""
function translate_volume!(m::GeoModel, tag, offset)
    caller="translate_volume!"
    t=_tag(tag,caller,3)
    transform_entities!(m,_affine_translation(
        _finite_vector3(offset,caller,"offset"),caller),[(3,t)];caller=caller)
    return t
end

"""
    dilate_volume!(model, tag, center, scale) -> tag

Dilate a native volume about the three-component finite `center`. `scale` must
be finite and positive. Explicit boundary topology and materialized boundary
Points move with the volume; compact primitive encodings are updated in place.
"""
function dilate_volume!(m::GeoModel, tag, center, scale)
    caller="dilate_volume!"
    t=_tag(tag,caller,3)
    s=_finite_scalar(scale,caller,"scale")
    s>0 || throw(ArgumentError("$caller: scale must be positive"))
    transform_entities!(m,_affine_dilation(center,s,caller),[(3,t)];caller=caller)
    return t
end

"""
    rotate_volume!(model, tag, axis, origin, angle) -> tag

Rotate a native volume about the finite, nonzero `axis` vector through
`origin` by the finite `angle` (radians). Explicit boundary topology and
materialized boundary Points move with the volume; compact primitive encodings
are updated in place. A materialized `add_box!` volume keeps its
`box_extents` encoding only while the rotated corners still form an
axis-aligned box — otherwise the shell topology survives and the encoding is
dropped.
"""
function rotate_volume!(m::GeoModel, tag, axis, origin, angle)
    caller="rotate_volume!"
    t=_tag(tag,caller,3)
    transform_entities!(m,_affine_rotation(axis,origin,angle,caller),
                        [(3,t)];caller=caller)
    return t
end

"""
    boolean_volumes!(model, op, a, b; tag=0) -> tag

Add a native Boolean volume combining existing volumes `a` and `b`. Supported
operations are `:union`, `:intersection`, and `:difference` (`a \\ b`). The
result owns operation-time snapshots of both operands, so later operand changes do
not change the Boolean geometry. A geometrically empty result (for example a
zero-volume face-touching intersection) binds nothing and returns `0`, matching
Gmsh's OCC kernel which reports an empty output list.
"""
function boolean_volumes!(m::GeoModel, op::Symbol, a, b; tag::Integer=0)
    caller="boolean_volumes!"
    op in (:union,:intersection,:difference) || throw(ArgumentError(
        "$caller: op must be :union, :intersection, or :difference"))
    ta=_tag(a,caller,3); tb=_tag(b,caller,3)
    haskey(m.volumes,ta) || throw(ArgumentError("$caller: unknown Volume[$ta]"))
    haskey(m.volumes,tb) || throw(ArgumentError("$caller: unknown Volume[$tb]"))
    requested=_tag(tag,caller,3)
    requested!=0 && (haskey(m.volumes,requested) || haskey(m.discrete,(3,requested))) && throw(ArgumentError(
        "$caller: Volume[$requested] already exists"))
    operand_a=_volume_surface(m,ta,caller)
    operand_b=_volume_surface(m,tb,caller)
    t=_alloc_tag!(m,3,requested,caller)
    haskey(m.volumes,t) && throw(ArgumentError(
        "$caller: automatic Volume[$t] already exists"))
    m.volumes[t]=Int[]
    m.booleans[t]=(op=op, a=ta, b=tb)
    m.boolean_operands[t]=(operand_a,operand_b)
    try
        _brep_materialize_boolean!(m,t,op,ta,tb,caller)
    catch
        delete!(m.volumes,t)
        delete!(m.booleans,t)
        delete!(m.boolean_operands,t)
        delete!(m.boolean_components,t)
        rethrow()
    end
    if isempty(m.volumes[t])
        # geometrically empty result (zero-volume contact) — OCC/Gmsh bind
        # nothing for it; drop the pre-allocated records. An automatic tag
        # was never claimed by any entity, so the watermark steps back.
        delete!(m.volumes,t)
        delete!(m.booleans,t)
        delete!(m.boolean_operands,t)
        delete!(m.boolean_components,t)
        requested==0 && m.next_tag[4]==t && (m.next_tag[4]-=1)
        return 0
    end
    return t
end

"""
    boolean_volumes_multi!(model, op, objects, tools; tag=0,
                           remove_object=true, remove_tool=true) -> out tags

Multi-operand Boolean over volume lists — the `.geo` `BooleanX{..}{..}` form
mirroring Gmsh's OCC `booleanOperator`. `objects` and `tools` are volume tag
lists (`tools` may be empty); `op` is `:union`, `:intersection`,
`:difference`, or `:fragments`. The arrangement decomposes into cells labeled
by operand membership; each kept cell materializes into its own result
volume — one volume per disconnected piece for a fuse, one per object/tool
membership cell for the other ops — with partition faces shared between
adjacent cells.

Tag/delete semantics follow `occBooleanPreserveNumbering`: an operand whose
only image is itself (no geometric interaction) keeps its tag and stays in
`out`; an operand whose image is a single piece is rebound to its tag when
the operand is removed; all other pieces get fresh tags after the surviving
maximum. `remove_object`/`remove_tool` are the `Delete;` markers of each
operand group. With an explicit `tag`, a multi-piece result is an error (one
tag cannot bind several volumes — OCC reports the same). Returns the output
volume tags in piece order, matching Gmsh's `outDimTags`.
"""
function boolean_volumes_multi!(m::GeoModel,op::Symbol,
        objects::AbstractVector{<:Integer},tools::AbstractVector{<:Integer};
        tag::Integer=0,remove_object::Bool=true,remove_tool::Bool=true,
        caller::AbstractString="boolean_volumes_multi!")
    op in (:union,:intersection,:difference,:fragments) || throw(ArgumentError(
        "$caller: op must be :union, :intersection, :difference, or :fragments"))
    remove_object isa Bool && remove_tool isa Bool || throw(ArgumentError(
        "$caller: remove flags must be Bool"))
    nobj=length(objects)
    nobj>=1 || throw(ArgumentError("$caller: at least one object operand"))
    tags=Int[_tag(t,caller,3) for t in objects]
    append!(tags,Int[_tag(t,caller,3) for t in tools])
    N=length(tags)
    length(unique(tags))==N || throw(ArgumentError(
        "$caller: Boolean operands must be distinct volumes"))
    for t in tags
        haskey(m.volumes,t) || throw(ArgumentError(
            "$caller: unknown Volume[$t]"))
    end
    requested=_tag(tag,caller,3)
    requested!=0 && (haskey(m.volumes,requested) ||
        haskey(m.discrete,(3,requested))) && throw(ArgumentError(
            "$caller: Volume[$requested] already exists"))
    # operand snapshots are operation-time — take them before any mutation
    meshes=[_volume_surface(m,t,caller) for t in tags]
    kept,touched=_brep_collect_faces_n(m,op,tags,nobj,caller)
    if isempty(kept)
        # empty result — OCC binds nothing; operands still get removed per
        # their group's Delete flag (nothing of them remains in the result)
        for i in 1:N
            rem=i<=nobj ? remove_object : remove_tool
            rem && _remove_volume_entity!(m,tags[i])
        end
        return Int[]
    end
    pieces,image,created=_brep_materialize_multi!(
        m,op,tags,kept,touched,N,nobj,caller)
    # operand dispositions (OCC preserve-numbering):
    #   image empty → the operand is deleted from the result — removed when
    #   its group carries Delete
    #   single materialized image + remove → the piece rebinds the tag
    #   unmodified (pseudo) → preserved in place
    try
        unbind=falses(N)
        piece_tag=zeros(Int,length(pieces))
        if requested==0
            for i in 1:N
                rem=i<=nobj ? remove_object : remove_tool
                imgs=image[i]
                if length(imgs)==1 && pieces[imgs[1]].pseudo==i
                    continue                       # Extent 0 — preserved
                elseif isempty(imgs)
                    unbind[i]=rem                  # IsDeleted
                elseif length(imgs)==1 && rem &&
                        count(j->image[j]==imgs &&
                            pieces[imgs[1]].pseudo!=j,1:N)==1
                    # the operand is the unique source whose entire image is
                    # this piece — a 1:1 modification that rebinds the tag
                    # (OCC preserve-numbering); a piece shared by several
                    # sole-image operands (a fuse) is not a bijection
                    unbind[i]=true
                    piece_tag[imgs[1]]=tags[i]
                else
                    unbind[i]=rem
                end
            end
        else
            length(pieces)>1 && throw(ArgumentError(
                "$caller: cannot bind $(length(pieces)) result volumes to " *
                "Volume[$requested]"))
            for i in 1:N
                rem=i<=nobj ? remove_object : remove_tool
                imgs=image[i]
                # an unmodified operand under an explicit tag is re-bound to
                # the result tag, not deleted
                if length(imgs)==1 && pieces[imgs[1]].pseudo==i && !rem
                    continue
                end
                unbind[i]=rem
            end
        end
        if requested!=0
            p=only(pieces)
            for i in 1:N
                (unbind[i] && i!=p.pseudo) || continue
                _remove_volume_entity!(m,tags[i])
            end
            if p.pseudo!=0
                # the sole result is an unmodified operand — it moves to the
                # requested tag when removed; when kept, OCC returns the
                # requested tag in outDimTags but the model keeps the operand
                # in place (no copy is materialized)
                rem_i=p.pseudo<=nobj ? remove_object : remove_tool
                rem_i && model_set_tag!(m,3,tags[p.pseudo],requested)
            else
                _bind_boolean_piece!(m,requested,op,tags,meshes,p,caller)
            end
            return [requested]
        end
        for i in 1:N
            unbind[i] || continue
            _remove_volume_entity!(m,tags[i])
        end
        # rebound pieces bind during the operand pass; fresh pieces take
        # sequential tags after the surviving maximum (OCC _multiBind order)
        for (pi,p) in enumerate(pieces)
            p.pseudo!=0 && continue
            piece_tag[pi]==0 && continue
            _bind_boolean_piece!(m,piece_tag[pi],op,tags,meshes,p,caller)
        end
        fresh=isempty(m.volumes) ? 0 : maximum(keys(m.volumes))
        out=Int[]
        for (pi,p) in enumerate(pieces)
            if p.pseudo!=0
                push!(out,tags[p.pseudo]);continue
            end
            if piece_tag[pi]==0
                fresh+=1
                while haskey(m.volumes,fresh)
                    fresh+=1
                end
                piece_tag[pi]=fresh
                _bind_boolean_piece!(m,fresh,op,tags,meshes,p,caller)
            end
            push!(out,piece_tag[pi])
        end
        # a multi-piece fuse extracts each volume's own component from the
        # shared fused snapshot mesh — link them like the binary path does
        # (preserved operands mesh through their own records)
        bound=[piece_tag[pi] for (pi,p) in enumerate(pieces)
               if p.pseudo==0]
        # pieces sharing one snapshot mesh — all fuse pieces share the fused
        # mesh; cell-op pieces with equal membership masks share a cell mesh —
        # each extracts its own connected components via the link marker
        # (preserved operands mesh through their own records)
        if op===:union
            if length(pieces)>1 && !isempty(bound)
                primary=bound[1]
                m.boolean_components[primary]=primary
                for t in bound[2:end]
                    m.boolean_components[t]=primary
                end
            end
        else
            cell_tags=Dict{UInt64,Vector{Int}}()
            for (pi,p) in enumerate(pieces)
                p.pseudo==0 || continue
                push!(get!(cell_tags,p.srcs,Int[]),piece_tag[pi])
            end
            for shared in values(cell_tags)
                length(shared)>1 || continue
                primary=shared[1]
                m.boolean_components[primary]=primary
                for t in shared[2:end]
                    m.boolean_components[t]=primary
                end
            end
        end
        return out
    catch
        for sh in created.shells
            delete!(m.surface_loops,sh)
        end
        _occ_materialize_rollback!(
            m,created.points,created.curves,created.loops,created.surfaces,0)
        rethrow()
    end
end

# Materialize one result piece's volume records at `t`. `meshes` are the
# operand snapshots; the piece's mesh cell is its membership mask — nothing
# for a fuse piece, whose boundary is the fused mesh's own component.
function _bind_boolean_piece!(m::GeoModel,t::Int,op::Symbol,
        tags::Vector{Int},meshes,p,caller)
    haskey(m.volumes,t) && throw(ArgumentError(
        "$caller: Volume[$t] already exists"))
    _alloc_tag!(m,3,t,caller)
    isempty(p.shells) && throw(ErrorException(
        "$caller: Boolean result piece has no boundary shells"))
    m.volumes[t]=p.shells
    m.booleans[t]=(op=op,operands=copy(tags))
    m.boolean_operands[t]=(meshes=meshes,
                           cell=op===:union ? nothing : p.srcs)
    return t
end

function _has_entity(m::GeoModel, dim::Int, tag::Int)
    haskey(m.discrete,(dim,tag)) && return true
    dim==0 && return haskey(m.points,tag)
    dim==1 && return haskey(m.curves,tag)
    dim==2 && return haskey(m.surfaces,tag)
    return haskey(m.volumes,tag)
end

function _physical_group_key(dim,tag,caller)
    d=_dimension(dim,caller)
    # Lookups take the raw tag — `.geo` can bind names to negative physical
    # tags, and a negative group tag simply resolves to nothing in the
    # observable (abs-keyed) view.
    t=_tag(tag,caller,d)
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
function add_physical_group!(m::GeoModel, dim::Integer, tags; tag::Integer=0,
                             name::AbstractString="", _literal_tag::Bool=false)
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
    requested=if _literal_tag
        tag isa Integer || throw(ArgumentError("$caller: tag must be an integer"))
        tag isa Bool && throw(ArgumentError("$caller: tag must not be Bool"))
        (-typemax(Int32)<=tag<=typemax(Int32)) || throw(ArgumentError(
            "$caller: tag magnitude exceeds Int32"))
        Int(tag)
    else
        _tag(tag,caller,d)
    end
    pt=_alloc_physical_tag(m,requested,caller;literal_tag=_literal_tag,
                           dimension=d)
    haskey(m.physical,(d,pt)) && throw(ArgumentError("$caller: Physical($d,$pt) already exists"))
    m.physical[(d,pt)]=ents
    for ent in ents
        # `GModel::addPhysicalGroup` stores `t>0 ? tag : -tag` on each entity.
        push!(get!(() -> Int[],m.entity_physicals,(d,ent)),
              ent>0 ? pt : -pt)
    end
    _physical_name_is_available(m,d,group_name) &&
        (m.physical_names[(d,pt)]=group_name)
    m.physical_tag_max=max(m.physical_tag_max,pt)
    return pt
end

"""
    set_physical_name!(model, dim, tag, name) -> name

Bind `name` to the physical tag `tag` in dimension `dim`, returning the
effective name bound at that tag. Names are unique within one entity dimension.
As in Gmsh 4.15.2 (`GModel::setPhysicalName`), `tag == 0` resolves to
`maxPhysicalNumber(dim) + 1` and the raw tag is bound whether or not a group
currently carries it; an empty name binds nothing, a name already used in that
dimension binds nothing at the requested tag, and an already named tag keeps
its existing name — use [`remove_physical_name!`](@ref) before assigning a
replacement.
"""
function set_physical_name!(m::GeoModel, dim::Integer, tag::Integer, name::AbstractString)
    caller="set_physical_name!"
    dimension=_dimension(dim,caller)
    tag isa Integer || throw(ArgumentError("$caller: tag must be an integer"))
    tag isa Bool && throw(ArgumentError("$caller: tag must not be Bool"))
    (-typemax(Int32)<=tag<=typemax(Int32)) || throw(ArgumentError(
        "$caller: tag magnitude exceeds Int32"))
    group_name=String(name)
    # `GModel::setPhysicalName` — a name already bound in this dimension wins
    # (`getPhysicalNumber` resolves it and the binding stays put, so nothing is
    # bound at the requested tag); `tag == 0` auto-assigns
    # `getMaxPhysicalNumber(dim) + 1`; the `(dim, tag)` key keeps its existing
    # name when already taken (std::map `insert` never overwrites). The raw
    # tag is bound whether or not a group currently carries it.
    isempty(group_name) ||
        (_physical_name_is_available(m,dimension,group_name) || return "")
    number=tag==0 ? _max_entity_physical_number(m,dimension)+1 : Int(tag)
    key=(dimension,number)
    existing=get(m.physical_names,key,"")
    isempty(existing) || return existing
    isempty(group_name) || (m.physical_names[key]=group_name)
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

Remove the selected `(dimension, physical_tag)` groups and their name bindings,
returning the number of selected tags that carried a group or name binding.
An empty selection clears every entity's physical memberships without erasing
the name table, matching `gmsh::model::removePhysicalGroups`; a targeted
selection erases each tag's name binding even when no group exists there, and
such a name-only erasure counts as a removal. Unknown valid groups are ignored,
all inputs are checked before mutation, and the global automatic-tag counter
remains monotonic.
"""
function remove_physical_groups!(m::GeoModel,dim_tags=())
    caller="remove_physical_groups!"
    selected=_physical_dim_tags(dim_tags,caller)
    removed=0
    if isempty(selected)
        # `gmsh::model::removePhysicalGroups({})` — `resetPhysicalGroups` +
        # `GModel::removePhysicalGroups`: every entity's physical list is
        # cleared but the name table is NOT erased.
        removed=length(m.physical)
        empty!(m.physical)
        for pnums in values(m.entity_physicals)
            empty!(pnums)
        end
        return removed
    end
    for key in selected
        # `GModel::removePhysicalGroup` strips every entity physical whose
        # `abs` equals the tag and erases the raw name binding — both run
        # unconditionally, so a name bound to a groupless tag is erased too.
        # A tag counts as removed when it carried a group or a name binding,
        # so a name-only erasure reports a removal like Gmsh's mutation does.
        had_group=haskey(m.physical,key)
        had_name=haskey(m.physical_names,key)
        had_group && delete!(m.physical,key)
        had_name && delete!(m.physical_names,key)
        (had_group || had_name) && (removed+=1)
        for ((d,e),pnums) in m.entity_physicals
            d==key[1] || continue
            filter!(pnum->abs(pnum)!=key[2],pnums)
        end
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

# Live members of a stored Physical group. `.geo` `Delete` leaves the stale
# member integers in the record (Gmsh's internals keep them and resurrect the
# membership if the tag is re-created), so resolution filters to entities that
# currently exist — the same view `GModel::getPhysicalGroups` produces.
function _physical_live_members(m::GeoModel,dimension::Int,members)
    return Int[tag for tag in members
               if _model_entity_known(m,dimension,tag)]
end

"""
    model_physical_groups(model, dim=-1) -> Vector{Tuple{Int,Int}}

Return detached Physical `(dimension, tag)` pairs in deterministic order. `dim=-1`
selects every dimension; `dim=0:3` filters the result.
"""
function model_physical_groups(m::GeoModel,dim=-1)
    dimension=_query_dimension(dim,"model_physical_groups")
    # A group whose members all resolve to nothing is not listed, matching
    # `GModel::getPhysicalGroups`, which builds the map from live entities.
    groups=Tuple{Int,Int}[
        key for (key,members) in m.physical
        if (dimension==-1 || key[1]==dimension) &&
           !isempty(_physical_live_members(m,key[1],members))]
    return sort!(groups)
end

"""
    model_entities_for_physical_group(model, dim, tag) -> Vector{Int}

Return detached, sorted entity tags for an existing Physical group.
"""
function model_entities_for_physical_group(m::GeoModel,dim,tag)
    key=_physical_group_key(dim,tag,"model_entities_for_physical_group")
    members=get(m.physical,key,nothing)
    # `getEntitiesForPhysicalGroup` reads the *derived* map: a raw group with
    # no live member (empty `={}`, or every member `Delete`d/unresolvable) is
    # absent there and reports "Physical ... does not exist".
    live=members===nothing ? Int[] : _physical_live_members(m,key[1],members)
    isempty(live) && throw(ArgumentError(
        "model_entities_for_physical_group: Physical$(key) does not exist"))
    return sort!(live)
end

"""
    model_physical_groups_for_entity(model, dim, tag) -> Vector{Int}

Return the Physical tags an existing geometry entity carries, in stored order.
Tags are signed like Gmsh's `getPhysicalGroupsForEntity`: an entity added to a
group through a negative `.geo` member or group tag reports the corresponding
negative tag.
"""
function model_physical_groups_for_entity(m::GeoModel,dim,tag)
    caller="model_physical_groups_for_entity"
    dimension=_dimension(dim,caller)
    entity_tag=_tag(tag,caller,dimension)
    entity_tag>=0 || throw(ArgumentError("$caller: entity tag must be non-negative"))
    key=(dimension,entity_tag)
    _has_entity(m,key...) || throw(ArgumentError(
        "$caller: entity $(key) does not exist"))
    return copy(get(m.entity_physicals,key,Int[]))
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
    for key in group_keys
        # `getEntitiesForPhysicalName` matches `abs(entity->physicals[j])`
        # against the *raw* bound tag — a name bound to a negative tag never
        # resolves, and a raw tag of 0 hits the pnum-0 view group.
        key[2]<0 && continue
        for entity in _physical_live_members(
                m,key[1],get(m.physical,key,Int[]))
            push!(entities,(key[1],entity))
        end
    end
    # Gmsh reports "Physical name '...' does not exist" whenever the lookup
    # yields no entities — including a bound name that resolves to nothing.
    isempty(entities) && throw(ArgumentError(
        "model_entities_for_physical_name: Physical name $(repr(group_name)) does not exist"))
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
            (dimension,tag) for tag in sort!(_physical_live_members(
                m,dimension,m.physical[(dimension,physical_tag)]))]
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
    members=get(m.physical,(d,t),nothing)
    members===nothing && return Int[]
    return _physical_live_members(m,d,members)
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

function _add_surface_point!(xs,ys,mesh_sizes,index,m::GeoModel,pid::Int,caller,
                             plane)
    haskey(index, pid) && return index[pid]
    haskey(m.points,pid) || throw(ArgumentError("$caller: unknown Point[$pid]"))
    p=m.points[pid]
    scale=max(1.0,hypot(p...))
    abs(_plane_offset(plane,p))<=1e-12*scale || throw(ArgumentError(
        "$caller: Point[$pid] is not coplanar with the surface " *
        "(off-plane distance $(_plane_offset(plane,p)))"))
    ax=plane.axes
    push!(xs,p[ax[1]]); push!(ys,p[ax[2]])
    # 0.0 marks an unconstrained point; `_surface_pslg` replaces it with the
    # shortest incident boundary edge once the segment table exists.
    push!(mesh_sizes,get(m.point_size,pid,0.0))
    index[pid]=length(xs)
    return index[pid]
end

function _add_surface_curve_point!(xs,ys,mesh_sizes,index,m::GeoModel,point,
                                   mesh_size::Float64,curve::Int,
                                   caller::AbstractString,plane)
    scale=max(1.0,hypot(point...))
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
        existing=mesh_sizes[matched_vertex]
        mesh_sizes[matched_vertex]=existing>0 ? min(existing,mesh_size) : mesh_size
        return matched_vertex
    end
    ax=plane.axes
    push!(xs,point[ax[1]]);push!(ys,point[ax[2]])
    push!(mesh_sizes,mesh_size)
    return length(xs)
end

@inline function _surface_curve_mesh_size(m::GeoModel,curve::Int,
                                          parameter::Float64,
                                          caller::AbstractString)
    a,b=m.curves[curve]
    # A point without an explicit size contributes the curve length, matching
    # the boundary-mesh-derived sizing used for unconstrained vertices.
    p,q=m.points[a],m.points[b]
    edge_length=hypot(q[1]-p[1],q[2]-p[2],q[3]-p[3])
    first_size=get(m.point_size,a,edge_length)
    last_size=get(m.point_size,b,edge_length)
    mesh_size=muladd(parameter,last_size-first_size,first_size)
    (isfinite(mesh_size) && mesh_size>0) || throw(ErrorException(
        "$caller: Curve[$curve] has an unrepresentable interpolated Point size"))
    return mesh_size
end

function _periodic_curve_point(m::GeoModel,curve::Int,parameter::Float64,
                               caller::AbstractString)
    _model_require_line_curve(m,curve,caller,"curve subdivision")
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

function _surface_curve_parameters(m::GeoModel,forced,curve::Int,signed::Int)
    # A `Degenerated` curve meshes to a single edge: it contributes only its
    # first endpoint, regardless of any stored parameter source (transfinite,
    # periodic, or size-at-params).
    curve in m.meshing.degenerated && return signed>0 ? (0.0,) : (1.0,)
    parameters=get(forced,curve,nothing)
    parameters===nothing && return signed>0 ? (0.0,) : (1.0,)
    return signed>0 ? @view(parameters[1:end-1]) :
                      Iterators.reverse(@view(parameters[2:end]))
end

function _surface_pslg(m::GeoModel,t::Int,forced,caller::AbstractString;
                       param_sizes::Dict{Tuple{Int,Float64},Float64}=
                           Dict{Tuple{Int,Float64},Float64}(),
                       plane=_model_surface_plane(m,t,caller;
                                                  allow_ruled=true))
    xs=Float64[];ys=Float64[];mesh_sizes=Float64[]
    segs=Tuple{Int,Int}[]
    index=Dict{Int,Int}()
    for loop_id in m.surfaces[t]
        _verify_loop_closed(m,loop_id,caller,"Surface[$t]")
        loop_idx=Int[]
        for signed in m.loops[loop_id]
            curve=abs(signed)
            _model_require_line_curve(m,curve,caller,"surface meshing")
            a,b=m.curves[curve]
            for parameter in _surface_curve_parameters(m,forced,curve,signed)
                vertex=if parameter==0
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,a,caller,plane)
                elseif parameter==1
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,b,caller,plane)
                else
                    point=_periodic_curve_point(m,curve,parameter,caller)
                    scale=max(1.0,hypot(point...))
                    abs(_plane_offset(plane,point))<=1e-12*scale ||
                        throw(ArgumentError(
                            "$caller: Curve[$curve] subdivision is not " *
                            "coplanar with Surface[$t]"))
                    _add_surface_curve_point!(
                        xs,ys,mesh_sizes,index,m,point,
                        get(param_sizes,(curve,parameter),
                            _surface_curve_mesh_size(
                                m,curve,parameter,caller)),
                        curve,caller,plane)
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
            _add_surface_point!(xs,ys,mesh_sizes,index,m,etag,caller,plane)
        elseif edim!=1
            throw(ArgumentError(
                "$caller: unsupported embedding dimension $edim"))
        end
    end
    for (edim,etag) in embedded
        edim==1 || continue
        haskey(m.curves,etag) || throw(ArgumentError(
            "$caller: unknown embedded Curve[$etag]"))
        _model_require_line_curve(m,etag,caller,"embedded-curve meshing")
        a,b=m.curves[etag]
        parameters=etag in m.meshing.degenerated ? nothing :
                   get(forced,etag,nothing)
        curve_nodes=Int[]
        if parameters===nothing
            push!(curve_nodes,_add_surface_point!(
                xs,ys,mesh_sizes,index,m,a,caller,plane))
            push!(curve_nodes,_add_surface_point!(
                xs,ys,mesh_sizes,index,m,b,caller,plane))
        else
            for parameter in parameters
                vertex=if parameter==0
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,a,caller,plane)
                elseif parameter==1
                    _add_surface_point!(xs,ys,mesh_sizes,index,m,b,caller,plane)
                else
                    point=_periodic_curve_point(m,etag,parameter,caller)
                    scale=max(1.0,hypot(point...))
                    abs(_plane_offset(plane,point))<=1e-12*scale ||
                        throw(ArgumentError(
                            "$caller: embedded Curve[$etag] subdivision is " *
                            "not coplanar with Surface[$t]"))
                    _add_surface_curve_point!(
                        xs,ys,mesh_sizes,index,m,point,
                        get(param_sizes,(etag,parameter),
                            _surface_curve_mesh_size(
                                m,etag,parameter,caller)),
                        etag,caller,plane)
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
    _fill_unsized_surface_vertices!(mesh_sizes,xs,ys,segs,internal)
    return xs,ys,mesh_sizes,segs,embedded,internal
end

# Vertices without an explicit Point size take the shortest incident edge,
# matching `_volume_boundary_size_field`'s boundary-mesh-derived sizing.
# A vertex with no incident edge falls back to the surface bounding scale.
function _fill_unsized_surface_vertices!(
        mesh_sizes::Vector{Float64},xs::Vector{Float64},
        ys::Vector{Float64},segs,internal)
    any(<=(0.0),mesh_sizes) || return nothing
    unsized=falses(length(mesh_sizes))
    for vertex in eachindex(mesh_sizes)
        if mesh_sizes[vertex]<=0.0
            unsized[vertex]=true
            mesh_sizes[vertex]=Inf
        end
    end
    # Vertices with an explicit size keep it verbatim, as in
    # `_volume_boundary_size_field`.
    for (a,b) in Iterators.flatten((segs,internal))
        len=hypot(xs[b]-xs[a],ys[b]-ys[a])
        unsized[a] && len<mesh_sizes[a] && (mesh_sizes[a]=len)
        unsized[b] && len<mesh_sizes[b] && (mesh_sizes[b]=len)
    end
    if any(isinf,mesh_sizes)
        scale=hypot(maximum(xs)-minimum(xs),maximum(ys)-minimum(ys))
        fallback=(isfinite(scale) && scale>0) ? scale : 1.0
        for vertex in eachindex(mesh_sizes)
            isinf(mesh_sizes[vertex]) && (mesh_sizes[vertex]=fallback)
        end
    end
    return nothing
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
    _model_require_line_curve(m,curve,caller,"curve mesh classification")
    a,b=m.curves[curve];p=m.points[a];q=m.points[b]
    vx=q[1]-p[1];vy=q[2]-p[2];vz=q[3]-p[3]
    length2=muladd(vx,vx,muladd(vy,vy,vz*vz))
    (isfinite(length2) && length2>0) || throw(ArgumentError(
        "$caller: Curve[$curve] has an unusable planar length"))
    length1=sqrt(length2)
    scale=max(1.0,hypot(p[1],p[2],p[3]),hypot(q[1],q[2],q[3]))
    geometric_tolerance=max(atol,128eps(Float64)*scale)
    cross_bound=(geometric_tolerance*length1)^2
    entries=Tuple{Float64,Int}[]
    @inbounds for node in 1:nnodes(mesh)
        eligible_nodes[node] || continue
        wx=mesh.coords[1,node]-p[1];wy=mesh.coords[2,node]-p[2]
        wz=mesh.coords[3,node]-p[3]
        cx=vy*wz-vz*wy;cy=vz*wx-vx*wz;cz=vx*wy-vy*wx
        muladd(cx,cx,muladd(cy,cy,cz*cz))<=cross_bound || continue
        parameter=muladd(wx,vx,muladd(wy,vy,wz*vz))/length2
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
    constraint.dim==3 && return (
        master_entity=Int(constraint.master_entity),
        slave_nodes=Int32[],master_nodes=Int32[],affine=constraint.affine)
    throw(ArgumentError(
        "model_periodic_nodes: unsupported periodic dimension $(constraint.dim)"))
end

function _model_mapping_matches(mesh::Mesh,constraint::ModelPeriodicConstraint,
                                mapping,caller::AbstractString;exact::Bool)
    # Orientation-only curve relations carry no affine: the node pairing is
    # parameter-based and was already verified while the mapping was built.
    constraint.affine===nothing && return true
    coefficients,translation=_periodic_affine_3x4(constraint.affine)
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
        # Orientation-only relations have no transform to snap with; the
        # synchronized parameters already pair the nodes on their own curves.
        constraint.affine===nothing && continue
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
explicit planar-shell volume produced by [`mesh_model_volume`](@ref). Volume
relations return empty node arrays with the stored master and affine: as in
Gmsh 4.15.2, volume periodicity is mesh-inert and carries no node
correspondence.
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
    _model_require_line_curve(m,curve,caller,"embedded-curve classification")
    start_point,stop_point=m.curves[curve]
    first_coordinate=m.points[start_point]
    last_coordinate=m.points[stop_point]
    vx=last_coordinate[1]-first_coordinate[1]
    vy=last_coordinate[2]-first_coordinate[2]
    vz=last_coordinate[3]-first_coordinate[3]
    length2=muladd(vx,vx,muladd(vy,vy,vz*vz))
    (isfinite(length2) && length2>0) || throw(ArgumentError(
        "$caller: embedded Curve[$curve] has an unusable planar length"))
    length1=sqrt(length2)
    scale=max(1.0,hypot(first_coordinate...),hypot(last_coordinate...))
    geometric_tolerance=max(atol,128eps(Float64)*scale)
    cross_bound=(geometric_tolerance*length1)^2
    parameter_tolerance=max(128eps(Float64),geometric_tolerance/length1)
    entries=Tuple{Float64,Int}[]
    @inbounds for node in 1:nnodes(mesh)
        wx=mesh.coords[1,node]-first_coordinate[1]
        wy=mesh.coords[2,node]-first_coordinate[2]
        wz=mesh.coords[3,node]-first_coordinate[3]
        cx=vy*wz-vz*wy;cy=vz*wx-vx*wz;cz=vx*wy-vy*wx
        muladd(cx,cx,muladd(cy,cy,cz*cz))<=cross_bound || continue
        parameter=muladd(wx,vx,muladd(wy,vy,wz*vz))/length2
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
    outgoing=Dict{Int,Vector{Tuple{Int,Union{Nothing,NTuple{16,Float64}}}}}()
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
            edges=get!(Vector{Tuple{Int,Union{Nothing,NTuple{16,Float64}}}},
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
        # `nothing` affines compare against tuples through `isless` — map them
        # to a shared sentinel so a corner carrying mixed orientation-only
        # and transform relations still sorts deterministically.
        sort!(edges;by=edge->(edge[1],edge[2]===nothing ? () : edge[2]))
    end

    # MSH periodic metadata permits only one relation per slave entity. Choose
    # a deterministic directed spanning forest through shared periodic corners,
    # matching Gmsh's endpoint-link convention without discarding curve links.
    visited=Set{Int}()
    parents=Dict{Int,Tuple{Int,Union{Nothing,NTuple{16,Float64}}}}()
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
                    Tuple{Int,Union{Nothing,NTuple{16,Float64}}}[])
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

    plane=_model_surface_plane(m,surface,caller;allow_ruled=true)
    constraints=_surface_periodic_constraints(m,surface,caller)
    geometric_tolerance=max(1e-12,
        isempty(constraints) ? 0.0 : maximum(c.atol for c in constraints))
    for node in axes(mesh.coords,2)
        coordinate=(mesh.coords[1,node],mesh.coords[2,node],
                    mesh.coords[3,node])
        offset=_plane_offset(plane,coordinate)
        abs(offset)<=geometric_tolerance || throw(ArgumentError(
            "$caller: input node $node is $offset off the Surface[$surface] " *
            "plane, outside the tolerance $geometric_tolerance"))
    end
    curve_tags=_model_projection_surface_curves(m,surface)
    isempty(curve_tags) && throw(ArgumentError(
        "$caller: Surface[$surface] has no boundary curves"))
    any(curve->curve in curve_tags,embedded_curves) && throw(ArgumentError(
        "$caller: a surface curve cannot be both bounding and embedded"))
    all_curve_tags=sort!(unique!(vcat(curve_tags,embedded_curves)))
    for curve in all_curve_tags,
        point in Iterators.flatten(
            (m.curves[curve],get(m.curve_control_points,curve,Int[])))
        p=m.points[point]
        scale=max(1.0,hypot(p...))
        abs(_plane_offset(plane,p))<=1e-12*scale || throw(ArgumentError(
            "$caller: Point[$point] is not coplanar with Surface[$surface]"))
    end
    for point in embedded_points
        p=m.points[point]
        scale=max(1.0,hypot(p...))
        abs(_plane_offset(plane,p))<=1e-12*scale || throw(ArgumentError(
            "$caller: embedded Point[$point] is not coplanar with " *
            "Surface[$surface]"))
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
    _model_require_line_curve(m,curve,caller,"embedded-curve classification")
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

# The coordinate-axis plane record for one planar surface: `axes` are the two
# kept coordinates of the nondegenerate projection `_model_surface_projection`
# selects, `k` is the dropped axis, `anchor` a boundary point on the plane and
# `normal` the unnormalized plane normal (its `k` component is nonzero because
# the projected anchor triangle is nondegenerate). `norm` is |normal|, used to
# turn dot products into absolute signed distances.
# `allow_ruled` admits `:ruled`/`:tric` surface-filling records for meshing
# paths: Gmsh's translational `Extrude` types every lateral `Surface`, and a
# ruled patch whose boundary certifies planar is planar in fact. Every
# boundary point is offset-checked while the PSLG is built, so a genuinely
# non-planar ruled boundary still fails explicitly there.
function _model_surface_plane(m::GeoModel,surface::Int,caller::AbstractString;
                              allow_ruled::Bool=false)
    haskey(m.surfaces,surface) || throw(ArgumentError(
        "$caller: unknown Surface[$surface]"))
    kind=_surface_type(m,surface)
    if kind!=:plane && !(allow_ruled && kind in (:ruled,:tric))
        _model_require_plane_surface(m,surface,caller,"plane geometry")
    end
    point_tags=Int[]
    for loop in m.surfaces[surface]
        haskey(m.loops,loop) || throw(ArgumentError(
            "$caller: Surface[$surface] references unknown Loop[$loop]"))
        for signed_curve in m.loops[loop]
            curve=abs(signed_curve)
            haskey(m.curves,curve) || throw(ArgumentError(
                "$caller: Loop[$loop] references unknown Curve[$curve]"))
            a,b=m.curves[curve]
            # Arc control points join the fit: the exact coplanarity check then
            # certifies the whole arc lies in the surface's plane.
            for point in Iterators.flatten(
                    ((a,b),get(m.curve_control_points,curve,Int[])))
                haskey(m.points,point) || throw(ArgumentError(
                    "$caller: Curve[$curve] references unknown Point[$point]"))
                push!(point_tags,point)
            end
        end
    end
    unique!(point_tags)
    coordinates=NTuple{3,Float64}[m.points[point] for point in point_tags]
    # An OCC circle edge collapses to a single endpoint; evaluated rim samples
    # stand in for the missing vertices so a cap bounded by one closed circle
    # still yields three non-collinear plane points.
    occ_samples,occ_sample_tags=_occ_surface_samples(m,surface,caller)
    append!(point_tags,occ_sample_tags)
    append!(coordinates,occ_samples)
    anchor,second,third,axes=_model_surface_projection(
        coordinates,point_tags,surface,caller)
    u=(second[1]-anchor[1],second[2]-anchor[2],second[3]-anchor[3])
    v=(third[1]-anchor[1],third[2]-anchor[2],third[3]-anchor[3])
    normal=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    k=axes==(1,2) ? 3 : axes==(1,3) ? 2 : 1
    return (axes=axes,k=k,anchor=anchor,normal=normal,
            norm=sqrt(normal[1]^2+normal[2]^2+normal[3]^2))
end

# Signed distance from `point` to `plane`, in absolute units. On an axis plane
# this reduces to the single dropped-coordinate difference the historical z=0
# checks measured.
@inline function _plane_offset(plane,point::NTuple{3,Float64})
    n=plane.normal;a=plane.anchor
    return (n[1]*(point[1]-a[1])+n[2]*(point[2]-a[2])+
            n[3]*(point[3]-a[3]))/plane.norm
end

# The dropped-axis coordinate of the on-plane point whose kept coordinates are
# `(u, v)` — the exact plane-equation solve. For axis planes the normal's kept
# components vanish and this returns `anchor[k]` bitwise.
@inline function _plane_dropped_coordinate(plane,u::Float64,v::Float64)
    n=plane.normal;a=plane.anchor;ax=plane.axes
    return a[plane.k]-(n[ax[1]]*(u-a[ax[1]])+n[ax[2]]*(v-a[ax[2]]))/n[plane.k]
end

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
        for signed_curve in m.loops[loop]
            _model_require_line_curve(m,abs(signed_curve),caller,
                                      "surface meshing")
        end
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
            _model_require_line_curve(m,curve,caller,
                                      "embedded-curve meshing")
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
    coefficients,translation=_periodic_affine_3x4(constraint.affine)
    output_coordinates=Matrix{Float64}(undef,3,nnodes(master_mesh))
    for node in 1:nnodes(master_mesh)
        key=ntuple(axis->_model_projection_coordinate_key(
            master_mesh.coords[axis,node]),3)
        master_point=get(coordinate_points,key,nothing)
        if master_point===nothing
            # Interior and edge nodes take the transform image, like
            # upstream's per-vertex `SPoint3::transform` copy.
            output_coordinates[:,node].=_model_affine_point(
                coefficients,translation,
                (master_mesh.coords[1,node],master_mesh.coords[2,node],
                 master_mesh.coords[3,node]),caller,node)
        else
            slave_point=point_map[master_point]
            output_coordinates[:,node].=m.points[slave_point]
        end
    end
    tris=master_mesh.tris
    coefficients_det=coefficients[1]*(coefficients[5]*coefficients[9]-
        coefficients[6]*coefficients[8])-
        coefficients[4]*(coefficients[2]*coefficients[9]-
        coefficients[3]*coefficients[8])+
        coefficients[7]*(coefficients[2]*coefficients[6]-
        coefficients[3]*coefficients[5])
    if coefficients_det<0
        # A mirrored transform flips triangle winding; swapping two vertices
        # keeps the copied mesh positively oriented.
        flipped=Matrix{Int32}(undef,3,size(tris,2))
        flipped[1,:].=tris[1,:];flipped[2,:].=tris[3,:];flipped[3,:].=tris[2,:]
        tris=flipped
    end
    output=Mesh(output_coordinates;tris=tris)
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
    coefficients,translation=_periodic_affine_3x4(constraint.affine)
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
    real_boundaries=_model_volume_boundary_surfaces(m,volume,caller)
    boundary_set=Set(abs.(real_boundaries))
    # A materialized curved primitive keeps its boundary entities for queries
    # but meshes through the analytic encoding: the boundary's curves and
    # surfaces are not classified into the tetrahedron mesh. Embedded entities
    # on those surfaces need curved-surface meshing — fail rather than drop.
    boundary_surfaces=if _implicit_volume_surface(m,volume)
        for surface in boundary_set
            nested_points,nested_curves=
                _model_surface_embedding_tags(m,surface,caller)
            (isempty(nested_points) && isempty(nested_curves)) ||
                throw(ArgumentError(
                    "$caller: Surface[$surface] bounds primitive " *
                    "Volume[$volume] and carries embedded entities; curved " *
                    "surface meshing is unsupported"))
        end
        Int[]
    else
        real_boundaries
    end
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
    # The declaration-time resolver already applied `GFace::setMeshMaster`'s
    # correspondence rules (boundary and embedded edges alike); reuse it so
    # projection sees the same induced pairs the declaration validated.
    point_map,pairs=_model_periodic_surface_edge_map(
        m,slave,master,constraint.affine,constraint.atol,caller)
    curve_map=Dict{Int,Int}(master_curve=>slave_curve
                          for (slave_curve,master_curve) in pairs)
    # Upstream's dim-0 records exist only for vertices that are endpoints of
    # a resolved curve — `GEdge::setMeshMaster` sets vertex masters there.
    # Standalone embedded points enter `point_map` for the correspondence
    # check but carry no vertex master, so they emit no point link.
    endpoints=Set{Int}()
    for (slave_curve,_) in pairs
        union!(endpoints,m.curves[slave_curve])
    end
    point_map=Dict{Int,Int}(master_point=>slave_point
        for (master_point,slave_point) in point_map
        if slave_point in endpoints)
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

    explicit_geometry=(isempty(m.volumes[volume]) ||
        _implicit_volume_surface(m,volume)) ? nothing :
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

# Full-3D point-to-segment test: the squared cross-product magnitude measures
# the off-line distance, so no projection-axis choice is needed and the check is
# correct on any coordinate plane or embedded chord.
function _on_segment(x, p, q; atol=1e-12)
    vx=q[1]-p[1];vy=q[2]-p[2];vz=q[3]-p[3]
    wx=x[1]-p[1];wy=x[2]-p[2];wz=x[3]-p[3]
    L2=muladd(vx,vx,muladd(vy,vy,vz*vz))
    L2>0 || return hypot(wx,wy,wz)<=atol
    cx=vy*wz-vz*wy;cy=vz*wx-vx*wz;cz=vx*wy-vy*wx
    muladd(cx,cx,muladd(cy,cy,cz*cz))<=(atol*sqrt(L2))^2 || return false
    t=muladd(wx,vx,muladd(wy,vy,wz*vz))/L2
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
        (_on_segment(pi,p,q; atol=atol) &&
         _on_segment(pj,p,q; atol=atol)) || continue
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
                                  caller::AbstractString;
                                  size_field::Union{Nothing,AbstractSizeField}=
                                      nothing,
                                  param_sizes::Dict{Tuple{Int,Float64},
                                                   Float64}=
                                      Dict{Tuple{Int,Float64},Float64}())
    if haskey(m.meshing.transfinite_surfaces,t)
        return _transfinite_surface_mesh(
            m,t,param_sizes,caller;size_field=size_field),NTuple{2,Int}[]
    end
    plane=_model_surface_plane(m,t,caller;allow_ruled=true)
    xs,ys,mesh_sizes,segs,embedded,internal=
        _surface_pslg(m,t,forced,caller;param_sizes=param_sizes,plane=plane)
    T=constrained_delaunay(xs,ys,segs; internal_segments=internal)
    base=if get(m.meshing.size_from_boundary,(2,t),true)
        _surface_point_size_field(T,xs,ys,mesh_sizes,t,caller)
    else
        lc=_pslg_default_size(xs,ys,caller,t)
        (x,y)->lc
    end
    callback=m.meshing.size_callback
    ax=plane.axes;k=plane.k
    sizefn=if size_field===nothing && callback===nothing
        base
    else
        function sized(x,y)
            h=base(x,y)
            if size_field!==nothing || callback!==nothing
                point=ntuple(3) do axis
                    axis==k ? _plane_dropped_coordinate(plane,x,y) :
                              axis==ax[1] ? x : y
                end
                size_field===nothing ||
                    (h=min(h,size_at(size_field,point[1],point[2],point[3],
                                     (2,t))))
                callback===nothing ||
                    (h=_apply_size_callback(callback,2,t,point[1],point[2],
                                            point[3],h,caller))
            end
            return h
        end
    end
    interior=refine!(T; min_angle_deg=min_angle_deg, size=sizefn)
    mesh=to_mesh(T; interior=interior)
    mesh=_consume_surface_attributes(m,t,mesh,caller)
    # `to_mesh` packs the projected (u,v) coordinates into rows 1,2; scatter
    # them onto the kept axes and solve the dropped coordinate on the plane.
    @inbounds for node in axes(mesh.coords,2)
        u=mesh.coords[1,node];v=mesh.coords[2,node]
        mesh.coords[ax[1],node]=u
        mesh.coords[ax[2],node]=v
        mesh.coords[k,node]=_plane_dropped_coordinate(plane,u,v)
    end
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntris(mesh)>0 || throw(ErrorException(
        "$caller: Surface[$t] produced no triangles"))
    return mesh,embedded
end

# Uniform fallback size for `size_from_boundary=false`: the surface's PSLG
# characteristic length, matching the model-scale default Gmsh falls back to.
function _pslg_default_size(xs,ys,caller::AbstractString,t::Int)
    isempty(xs) && throw(ArgumentError(
        "$caller: Surface[$t] has no boundary points"))
    dx=maximum(xs)-minimum(xs);dy=maximum(ys)-minimum(ys)
    lc=hypot(dx,dy)
    (isfinite(lc) && lc>0) || throw(ArgumentError(
        "$caller: Surface[$t] has a degenerate bounding box"))
    return lc
end

# Gmsh's size callback contract: `(dim, tag, x, y, z, lc) -> size`; a
# non-positive return means "no constraint" and keeps the incoming `lc`.
function _apply_size_callback(callback,dim::Int,tag::Int,x::Float64,
                              y::Float64,z::Float64,lc::Float64,
                              caller::AbstractString)
    h=try
        callback(dim,tag,x,y,z,lc)
    catch err
        err isa InterruptException && rethrow()
        throw(ErrorException(
            "$caller: mesh size callback failed for ($dim,$tag) at " *
            "($x,$y,$z): $(sprint(showerror,err))"))
    end
    h isa Real || throw(ErrorException(
        "$caller: mesh size callback returned a non-numeric value " *
        "$(repr(h)) for ($dim,$tag)"))
    h isa Bool && throw(ErrorException(
        "$caller: mesh size callback returned Bool for ($dim,$tag)"))
    value=Float64(h)
    !isfinite(value) && return lc
    value<=0 && return lc
    return min(lc,value)
end

# Post-refinement attribute passes for a surface mesh: `smoothing` runs that
# many boundary-preserving Laplacian iterations; `reverse` flips triangle
# orientation, matching Gmsh's `setReverse`/`setSmoothing` semantics.
function _consume_surface_attributes(m::GeoModel,t::Int,mesh::Mesh,
                                     caller::AbstractString)
    iterations=get(m.meshing.smoothing,(2,t),0)
    iterations>0 && (mesh=_laplacian_smooth_surface(mesh,iterations,caller,t))
    get(m.meshing.reverse,(2,t),false) &&
        (mesh=_reversed_surface_mesh(mesh,caller,t))
    return mesh
end

function _laplacian_smooth_surface(mesh::Mesh,iterations::Int,
                                   caller::AbstractString,t::Int)
    boundary,_=_surface_boundary_topology(mesh,caller)
    coords=Matrix{Float64}(mesh.coords)
    neighbors=[Int32[] for _ in 1:nnodes(mesh)]
    @inbounds for cell in axes(mesh.tris,2),e in ((1,2),(2,3),(3,1))
        a=mesh.tris[e[1],cell];b=mesh.tris[e[2],cell]
        push!(neighbors[a],b);push!(neighbors[b],a)
    end
    for list in neighbors
        sort!(unique!(list))
    end
    for _ in 1:iterations
        next=copy(coords)
        @inbounds for node in axes(coords,2)
            boundary[node] && continue
            list=neighbors[node]
            isempty(list) && continue
            sx=0.0;sy=0.0
            for other in list
                sx+=coords[1,other];sy+=coords[2,other]
            end
            next[1,node]=sx/length(list)
            next[2,node]=sy/length(list)
        end
        coords=next
    end
    return Mesh(coords;segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,
                tet_tag=mesh.tet_tag)
end

function _reversed_surface_mesh(mesh::Mesh,caller::AbstractString,t::Int)
    tris=Matrix{Int32}(mesh.tris)
    @inbounds for cell in axes(tris,2)
        tris[2,cell],tris[3,cell]=tris[3,cell],tris[2,cell]
    end
    return Mesh(mesh.coords;segs=mesh.segs,tris=tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,
                tet_tag=mesh.tet_tag)
end

# Post-generation attribute passes for a volume mesh: `smoothing` runs that
# many boundary-preserving Laplacian iterations over the tetrahedron edge
# graph; `reverse` flips tetrahedron orientation, matching Gmsh's
# `setReverse`/`setSmoothing` semantics on dim-3 entities.
function _consume_volume_attributes(m::GeoModel,t::Int,mesh::Mesh,
                                    caller::AbstractString)
    iterations=get(m.meshing.smoothing,(3,t),0)
    iterations>0 &&
        (mesh=_laplacian_smooth_volume(mesh,iterations,caller,t))
    get(m.meshing.reverse,(3,t),false) &&
        (mesh=_reversed_volume_mesh(mesh,caller,t))
    return mesh
end

function _laplacian_smooth_volume(mesh::Mesh,iterations::Int,
                                  caller::AbstractString,t::Int)
    boundary=falses(nnodes(mesh))
    @inbounds for face in first(boundary_faces(mesh.tets)),node in face
        boundary[node]=true
    end
    neighbors=[Int32[] for _ in 1:nnodes(mesh)]
    @inbounds for cell in axes(mesh.tets,2),
            edge in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
        a=mesh.tets[edge[1],cell];b=mesh.tets[edge[2],cell]
        push!(neighbors[a],b);push!(neighbors[b],a)
    end
    for list in neighbors
        sort!(unique!(list))
    end
    coords=Matrix{Float64}(mesh.coords)
    for _ in 1:iterations
        next=copy(coords)
        @inbounds for node in axes(coords,2)
            boundary[node] && continue
            list=neighbors[node]
            isempty(list) && continue
            sx=0.0;sy=0.0;sz=0.0
            for other in list
                sx+=coords[1,other];sy+=coords[2,other];sz+=coords[3,other]
            end
            next[1,node]=sx/length(list)
            next[2,node]=sy/length(list)
            next[3,node]=sz/length(list)
        end
        coords=next
    end
    return Mesh(coords;segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,
                tet_tag=mesh.tet_tag)
end

function _reversed_volume_mesh(mesh::Mesh,caller::AbstractString,t::Int)
    tets=Matrix{Int32}(mesh.tets)
    @inbounds for cell in axes(tets,2)
        tets[1,cell],tets[2,cell]=tets[2,cell],tets[1,cell]
    end
    return Mesh(mesh.coords;segs=mesh.segs,tris=mesh.tris,tets=tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,
                tet_tag=mesh.tet_tag)
end

# Transfinite volume fill for an explicit six-surface volume — Gmsh's
# `setTransfiniteVolume`. The 12 boundary edges must all carry transfinite
# curve specs; the 8 corner vertices (points incident to exactly 3 boundary
# edges) are ordered in Gmsh's canonical (s0..s7) order, either from the
# stored corner list or automatically from the boundary edge graph.
function _transfinite_volume_mesh(m::GeoModel,t::Int,caller::AbstractString)
    get(m.embeds,(3,t),NTuple{2,Int}[]) |> isempty || throw(ArgumentError(
        "$caller: transfinite Volume[$t] cannot carry embedded entities"))
    # `TransfQuadTri` selects Gmsh's HAVE_QUADTRI path — hexa/prism elements
    # with boundary-diagonal subdivision at unrecombined faces. The native
    # kernel emits tetrahedra only, so the flag is an explicit blocker here.
    t in m.meshing.quad_tri && throw(ArgumentError(
        "$caller: TransfQuadTri Volume[$t] requires the QuadTri hexahedral " *
        "transfinite algorithm, which Tessella does not implement"))
    boundaries=_model_volume_boundary_surfaces(m,t,caller)
    length(boundaries)==6 || throw(ArgumentError(
        "$caller: transfinite Volume[$t] requires exactly 6 boundary " *
        "surfaces; found $(length(boundaries))"))
    edge_curve=Dict{NTuple{2,Int},Int}()
    neighbors=Dict{Int,Set{Int}}()
    for signed_surface in boundaries
        for loop in m.surfaces[abs(signed_surface)]
            for signed in m.loops[loop]
                curve=abs(signed)
                a,b=m.curves[curve]
                key=a<b ? (a,b) : (b,a)
                if !haskey(edge_curve,key)
                    edge_curve[key]=curve
                    push!(get!(neighbors,a,Set{Int}()),b)
                    push!(get!(neighbors,b,Set{Int}()),a)
                end
            end
        end
    end
    corners=sort!(collect(p for (p,adj) in pairs(neighbors) if
                          length(adj)>=3))
    length(corners)==8 || throw(ArgumentError(
        "$caller: transfinite Volume[$t] requires a hexahedral boundary " *
        "topology (8 corner Points); found $(length(corners))"))
    length(edge_curve)==12 || throw(ArgumentError(
        "$caller: transfinite Volume[$t] requires 12 boundary Curves; " *
        "found $(length(edge_curve))"))
    stored=m.meshing.transfinite_volumes[t]
    ordered=if isempty(stored)
        s0=corners[1]
        near=sort!(collect(neighbors[s0]))
        length(near)==3 || throw(ArgumentError(
            "$caller: transfinite Volume[$t] boundary is not a cube " *
            "edge graph"))
        s1,s3,s4=near
        # The kernel requires orient3(s0,s1,s3,s4) < 0 (positive canonical
        # Gmsh order); swapping the u/v neighbors flips the orientation.
        if orient3(m.points[s0],m.points[s1],m.points[s3],
                   m.points[s4])>0
            s1,s3=s3,s1
        end
        function common(x,y)
            shared=setdiff(neighbors[x]∩neighbors[y],(s0,))
            length(shared)==1 || throw(ArgumentError(
                "$caller: transfinite Volume[$t] boundary is not a cube " *
                "edge graph"))
            return only(shared)
        end
        s2=common(s1,s3);s5=common(s1,s4);s7=common(s3,s4)
        remaining=setdiff(corners,(s0,s1,s2,s3,s4,s5,s7))
        length(remaining)==1 || throw(ArgumentError(
            "$caller: transfinite Volume[$t] boundary is not a cube " *
            "edge graph"))
        (s0,s1,s2,s3,s4,s5,only(remaining),s7)
    else
        length(stored)==8 || throw(ArgumentError(
            "$caller: transfinite Volume[$t] supports only 8-corner blocks"))
        Set(stored)==Set(corners) || throw(ArgumentError(
            "$caller: transfinite Volume[$t] corners must be its 8 " *
            "boundary corner Points"))
        Tuple(stored)
    end
    s=ordered
    function edge_count(a,b)
        curve=get(edge_curve,a<b ? (a,b) : (b,a),0)
        curve==0 && throw(ArgumentError(
            "$caller: transfinite Volume[$t] corner pair ($a,$b) is not " *
            "a boundary edge"))
        spec=get(m.meshing.transfinite_curves,curve,nothing)
        spec===nothing && throw(ArgumentError(
            "$caller: transfinite Volume[$t] requires boundary " *
            "Curve[$curve] to be transfinite"))
        return spec.num_nodes-1
    end
    us=(edge_count(s[1],s[2]),edge_count(s[4],s[3]),
        edge_count(s[5],s[6]),edge_count(s[8],s[7]))
    vs=(edge_count(s[1],s[4]),edge_count(s[2],s[3]),
        edge_count(s[5],s[8]),edge_count(s[6],s[7]))
    ws=(edge_count(s[1],s[5]),edge_count(s[2],s[6]),
        edge_count(s[4],s[8]),edge_count(s[3],s[7]))
    for (direction,family) in (("u",us),("v",vs),("w",ws))
        allequal(family) || throw(ArgumentError(
            "$caller: transfinite Volume[$t] $direction-direction edges " *
            "have mismatched node counts $family"))
    end
    return mesh_transfinite_volume(
        NTuple{3,Float64}[m.points[p] for p in s],
        (us[1],vs[1],ws[1]);volume_tag=t)
end

# Boundary-derived size field for a volume — Gmsh MeshSizeFromBoundary
# semantics: interior sizes extend the sizes prescribed at boundary vertices.
# Sizes come from sized boundary Points coincident with surface-mesh nodes;
# unsized boundary nodes inherit their minimum incident edge length. When no
# sized boundary vertex exists there is nothing to propagate.
function _volume_boundary_size_field(m::GeoModel,t::Int,surface::Mesh,
                                     caller::AbstractString)
    point_size=Dict{NTuple{3,Float64},Float64}()
    for (tag,coordinate) in m.points
        value=get(m.point_size,tag,0.0)
        (isfinite(value) && value>0) || continue
        key=(_model_projection_coordinate_key(coordinate[1]),
             _model_projection_coordinate_key(coordinate[2]),
             _model_projection_coordinate_key(coordinate[3]))
        point_size[key]=min(get(point_size,key,Inf),value)
    end
    isempty(point_size) && return nothing
    values=fill(Inf,nnodes(surface))
    matched=falses(nnodes(surface))
    @inbounds for node in 1:nnodes(surface)
        key=(_model_projection_coordinate_key(surface.coords[1,node]),
             _model_projection_coordinate_key(surface.coords[2,node]),
             _model_projection_coordinate_key(surface.coords[3,node]))
        value=get(point_size,key,0.0)
        value>0 && (values[node]=value;matched[node]=true)
    end
    any(matched) || return nothing
    # Boundary vertices that do not carry a Point size take the shortest
    # incident boundary edge, matching Gmsh's boundary-mesh-derived sizing.
    # Vertices with an explicit Point size keep it verbatim.
    @inbounds for cell in axes(surface.tris,2),edge in ((1,2),(2,3),(3,1))
        a=surface.tris[edge[1],cell];b=surface.tris[edge[2],cell]
        dx=surface.coords[1,b]-surface.coords[1,a]
        dy=surface.coords[2,b]-surface.coords[2,a]
        dz=surface.coords[3,b]-surface.coords[3,a]
        len=hypot(dx,dy,dz)
        !matched[a] && len<values[a] && (values[a]=len)
        !matched[b] && len<values[b] && (values[b]=len)
    end
    # An isolated boundary vertex with no sized neighbor and no incident
    # triangle edge falls back to the surface's bounding scale.
    scale=hypot(maximum(surface.coords[1,:])-minimum(surface.coords[1,:]),
                maximum(surface.coords[2,:])-minimum(surface.coords[2,:]),
                maximum(surface.coords[3,:])-minimum(surface.coords[3,:]))
    fallback=(isfinite(scale) && scale>0) ? scale : 1.0
    @inbounds for node in 1:nnodes(surface)
        isfinite(values[node]) || (values[node]=fallback)
    end
    # PostViewField is an AbstractField; the refiner wants an
    # AbstractSizeField, so adapt through FunctionSize.
    view=PostViewField(surface.coords,values;crop_negative=false)
    return FunctionSize(function(x,y,z)
        field_value(view,x,y,z)
    end)
end

# Compose the effective 3-D size source for `mesh_model_volume`: an explicit
# `size_field`, boundary-derived sizes, and the stored size callback wrap in
# that order (callback sees the incoming `lc` and may tighten it).
function _volume_size_field(m::GeoModel,t::Int,surface::Mesh,
                            size_field::Union{Nothing,AbstractSizeField},
                            caller::AbstractString)
    fields=AbstractSizeField[]
    boundary=_volume_boundary_size_field(m,t,surface,caller)
    boundary===nothing || push!(fields,boundary)
    size_field===nothing || push!(fields,size_field)
    callback=m.meshing.size_callback
    if callback===nothing
        isempty(fields) && return nothing
        return length(fields)==1 ? fields[1] : MinSize(Tuple(fields))
    end
    base=if isempty(fields)
        scale=hypot(maximum(surface.coords[1,:])-minimum(surface.coords[1,:]),
                    maximum(surface.coords[2,:])-minimum(surface.coords[2,:]),
                    maximum(surface.coords[3,:])-minimum(surface.coords[3,:]))
        ConstantSize((isfinite(scale) && scale>0) ? scale : 1.0)
    elseif length(fields)==1
        fields[1]
    else
        MinSize(Tuple(fields))
    end
    return FunctionSize(function(x,y,z)
        _apply_size_callback(callback,3,t,x,y,z,
                             size_at(base,x,y,z),caller)
    end)
end

# Rotate/reverse a transfinite surface's discretized side chains so the first
# pinned corner leads the boundary walk — the model-level half of Gmsh's
# `findTransfiniteCorners` explicit-corner branch. `corners` must contain
# exactly the loop's junction Point tags in forward or reversed cyclic order;
# a reversed order flips the walk as Gmsh's m_vertices reversal does.
function _apply_pinned_surface_corners(m::GeoModel,t::Int,
                                       signed_curves,curve_points,
                                       corners,nside,caller)
    length(corners)==nside || throw(ArgumentError(
        "$caller: transfinite Surface[$t] corner list has $(length(corners)) " *
        "entries for a $nside-curve boundary"))
    junctions=Int[signed>0 ? m.curves[signed][1] : m.curves[-signed][2]
                  for signed in signed_curves]
    sort(corners)==sort(junctions) || throw(ArgumentError(
        "$caller: transfinite Surface[$t] corner list must be the boundary " *
        "junction Points $junctions (got $corners)"))
    start=findfirst(==(corners[1]),junctions)::Int
    if all(i->corners[i]==junctions[mod1(start+i-1,nside)],1:nside)
        return [curve_points[mod1(start+i-1,nside)] for i in 1:nside]
    end
    all(i->corners[i]==junctions[mod1(start-i+1,nside)],1:nside) ||
        throw(ArgumentError(
            "$caller: transfinite Surface[$t] corner list order is " *
            "inconsistent with the curve loop orientation"))
    return [reverse(curve_points[mod1(start-i,nside)]) for i in 1:nside]
end

# Transfinite interpolation for a planar surface whose boundary curves are all
# transfinite — Gmsh's `setTransfiniteSurface`. A 4-sided loop routes through
# `mesh_transfinite_patch`, the structured Coons grid triangulated per quad
# cell with Gmsh's arrangement parity; a 3-sided loop routes through
# `mesh_transfinite_triangle_collapsed` (the legacy `Mesh.TransfiniteTri=0`
# collapsed-quadrilateral algorithm, the default) or `mesh_transfinite_triangle`
# (the compact `TransfiniteTri=1` algorithm, selected by `set_transfinite_tri!`
# and requiring equal node counts on all three sides).
function _transfinite_surface_mesh(m::GeoModel,t::Int,
                                   param_sizes::Dict{Tuple{Int,Float64},
                                                    Float64},
                                   caller::AbstractString;
                                   size_field::Union{Nothing,
                                                     AbstractSizeField}=
                                       nothing)
    spec=m.meshing.transfinite_surfaces[t]
    get(m.embeds,(2,t),NTuple{2,Int}[]) |> isempty || throw(ArgumentError(
        "$caller: transfinite Surface[$t] cannot carry embedded entities"))
    loops=m.surfaces[t]
    length(loops)==1 || throw(ArgumentError(
        "$caller: transfinite Surface[$t] requires exactly one curve loop"))
    # `findVertices` skips `degenerate(0)` curves entirely — a `Degenerated`
    # boundary curve drops out of the side chain, so the surface meshes as
    # `nboundary - ndegenerated`-sided.
    signed_curves=[s for s in m.loops[only(loops)]
                   if !(abs(s) in m.meshing.degenerated)]
    nside=length(signed_curves)
    nside in (3,4) || throw(ArgumentError(
        "$caller: transfinite Surface[$t] requires a 3- or 4-curve boundary " *
        "after skipping degenerated curves (got $nside)"))
    curve_points=Vector{Vector{NTuple{3,Float64}}}(undef,nside)
    for (position,signed) in enumerate(signed_curves)
        curve=abs(signed)
        cspec=get(m.meshing.transfinite_curves,curve,nothing)
        cspec===nothing && throw(ArgumentError(
            "$caller: transfinite Surface[$t] requires boundary Curve[$curve] " *
            "to be transfinite"))
        params=_transfinite_parameters(
            m,cspec.num_nodes,cspec.kind,cspec.coef,caller,curve;
            reversed=cspec.reversed)
        points=[_periodic_curve_point(m,curve,p,caller) for p in params]
        signed<0 && reverse!(points)
        curve_points[position]=points
    end
    isempty(spec.corners) || (curve_points=_apply_pinned_surface_corners(
        m,t,signed_curves,curve_points,spec.corners,nside,caller))
    plane=_model_surface_plane(m,t,caller)
    for p in Iterators.flatten(curve_points)
        scale=max(1.0,hypot(p...))
        abs(_plane_offset(plane,p))<=1e-12*scale || throw(ArgumentError(
            "$caller: transfinite Surface[$t] boundary is not coplanar"))
    end
    # Corner consistency: each side ends where the next begins. The kernels
    # require bitwise-identical shared corners, while `_periodic_curve_point`
    # may differ from the vertex coordinate by an ulp at parameter 1
    # (`p + (q - p)` need not equal `q` exactly) — weld after auditing.
    tolerance=1e-9*max(1.0,maximum(
        p->maximum(abs,p),Iterators.flatten(curve_points)))
    for position in 1:nside
        _points_close(curve_points[position][end],
                      curve_points[mod1(position+1,nside)][1],
                      tolerance) || throw(ArgumentError(
            "$caller: transfinite Surface[$t] boundary corner $position is " *
            "inconsistent"))
    end
    for position in 1:nside
        curve_points[position][end]=curve_points[mod1(position+1,nside)][1]
    end
    if nside==3
        s1,s2,s3=curve_points
        kernel=if m.meshing.transfinite_tri==1
            (length(s1)==length(s2) && length(s1)==length(s3)) ||
                throw(ArgumentError(
                    "$caller: transfinite Surface[$t] has mismatched " *
                    "boundary curve node counts ($(length(s1)), " *
                    "$(length(s2)), $(length(s3)))"))
            mesh_transfinite_triangle(s1,s2,s3;arrangement=spec.arrangement)
        else
            mesh_transfinite_triangle_collapsed(
                s1,s2,s3;arrangement=spec.arrangement,
                allow_corner_rotation=isempty(spec.corners))
        end
        # The entity cache stores the untagged simplex complex; boundary
        # curves are not meshed by generate(2).
        mesh=Mesh(kernel.coords;tris=kernel.tris)
        return _consume_surface_attributes(m,t,mesh,caller)
    end
    bottom,right,top,left=curve_points
    kernel=mesh_transfinite_patch(bottom,right,top,left;
                                  arrangement=spec.arrangement)
    mesh=Mesh(kernel.coords;tris=kernel.tris)
    return _consume_surface_attributes(m,t,mesh,caller)
end

@inline function _points_close(p,q,tolerance)
    return abs(p[1]-q[1])<=tolerance && abs(p[2]-q[2])<=tolerance &&
           abs(p[3]-q[3])<=tolerance
end

function _validate_surface_embeddings(m::GeoModel,mesh::Mesh,embedded,
                                      caller::AbstractString)
    for (edim,etag) in embedded
        if edim==0
            p=m.points[etag]
            _node_at(mesh,p)==0 && throw(ErrorException(
                "$caller: embedded Point[$etag] at $p is not a mesh node"))
        elseif edim==1
            _model_require_line_curve(m,etag,caller,
                                      "embedded-curve verification")
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
                            max_periodic_passes=8,
                            size_field::Union{Nothing,AbstractSizeField}=nothing)
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
    # A slave surface takes its master's mesh verbatim, like upstream's
    # meshGFace copy path — the master is meshed on demand so the relation
    # holds even when the caller never meshed it directly.
    surface_constraint=nothing
    for constraint in model_periodic_constraints(m)
        constraint.dim==2 && Int(constraint.slave_entity)==t || continue
        surface_constraint=constraint
        break
    end
    if surface_constraint!==nothing
        master_mesh=mesh_model_surface(
            m,Int(surface_constraint.master_entity);
            min_angle_deg=min_angle_deg,
            max_periodic_passes=max_periodic_passes,size_field=size_field)
        output=_model_periodic_surface_mesh(
            m,master_mesh,surface_constraint,caller)
        return output
    end
    forced=Dict{Int,Vector{Float64}}()
    param_sizes=_attribute_forced_parameters(m,t,forced,caller)
    mesh=nothing;embedded=NTuple{2,Int}[]
    for pass in 1:npasses
        mesh,embedded=_mesh_model_surface_once(
            m,t,forced,min_angle_deg,caller;size_field=size_field,
            param_sizes=param_sizes)
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
        elseif constraint.dim==3
            continue
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

# A volume carrying a native curved-solid encoding meshes through the analytic
# tessellation even though its OCC boundary topology is materialized — the
# planar explicit-shell path cannot represent Cylinder/Sphere/Cone faces. Box
# volumes deliberately stay explicit: their faces are planar. Boolean results
# likewise keep their materialized boundary for queries but mesh through the
# operand snapshot (`m.booleans`/`m.boolean_operands`).
_implicit_volume_surface(m::GeoModel,t::Int) =
    haskey(m.cylinders,t) || haskey(m.spheres,t) || haskey(m.cones,t) ||
    haskey(m.booleans,t)

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
        haskey(m.boolean_operands,t) || throw(ArgumentError(
            "$caller: Boolean Volume[$t] has no owned operand snapshot; " *
            "recreate the Boolean in a fresh GeoModel"))
        return _boolean_result_surface(m,t,caller)
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
function mesh_model_volume(m::GeoModel, tag::Integer;
                           size_field::Union{Nothing,AbstractSizeField}=nothing)
    caller="mesh_model_volume"
    t=_tag(tag,caller,3)
    haskey(m.volumes,t) || throw(ArgumentError("$caller: unknown Volume[$t]"))
    if haskey(m.meshing.transfinite_volumes,t)
        mesh=_transfinite_volume_mesh(m,t,caller)
        mesh=_consume_volume_attributes(m,t,mesh,caller)
        # model_to_mixed expects an untagged pure tetrahedron complex —
        # boundary faces are re-derived from tet faces downstream and
        # ownership comes from the model, like `tetrahedralize` output.
        mesh=Mesh(mesh.coords;tets=mesh.tets)
        reversed=get(m.meshing.reverse,(3,t),false)
        diagnostic=validate(mesh;require_positive_tets=!reversed)
        diagnostic.ok || throw(ErrorException(
            "$caller: transfinite Volume[$t] produced an invalid mesh — " *
            join(diagnostic.messages,"; ")))
        ntets(mesh)>0 || throw(ErrorException(
            "$caller: Volume[$t] produced no tetrahedra"))
        return mesh
    end
    explicit_geometry=(isempty(m.volumes[t]) ||
        _implicit_volume_surface(m,t)) ? nothing :
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
    effective_field=_volume_size_field(m,t,surface,size_field,caller)
    effective_field===nothing ||
        (mesh=refine_to_size(mesh,effective_field;entity=(3,t),
                             best_effort=true))
    iterations=get(m.meshing.smoothing,(3,t),0)
    iterations>0 &&
        (mesh=_laplacian_smooth_volume(mesh,iterations,caller,t))
    sort!(unique!(line_tags))
    for curve in line_tags
        _model_require_line_curve(m,curve,caller,"embedded-curve recovery")
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
    reversed=get(m.meshing.reverse,(3,t),false)
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    ntets(mesh)>0 || throw(ErrorException("$caller: Volume[$t] produced no tetrahedra"))
    explicit_geometry===nothing || _model_certify_explicit_volume_semantics(
        mesh,t,explicit_geometry.expected_volume,
        explicit_geometry.comparison_scale,caller)
    # `setReverse` flips orientation last so all certifications above run on
    # the positively-oriented complex; structural checks still apply after.
    reversed || return mesh
    mesh=_reversed_volume_mesh(mesh,caller,t)
    post=validate(mesh;require_positive_tets=false)
    post.ok || throw(ErrorException(
        "$caller: reversed mesh is structurally invalid — " *
        join(post.messages,"; ")))
    return mesh
end

# Mesh a dimension-2 compound as a single patch, Gmsh `setCompound` semantics:
# curves shared by two member surfaces become internal constraints rather than
# boundaries, the union is triangulated once, and each resulting triangle is
# attributed back to its member surface by point-in-polygon on that member's
# own boundary. Returns `(member_tag, Mesh)` pairs in member order.
function _compound_surface_meshes(m::GeoModel,members::Vector{Int},
                                  caller::AbstractString;
                                  size_field::Union{Nothing,AbstractSizeField}=
                                      nothing)
    for tag in members
        haskey(m.meshing.transfinite_surfaces,tag) && throw(ArgumentError(
            "$caller: transfinite Surface[$tag] cannot join a compound"))
    end
    counts=Dict{Int,Int}()
    for tag in members,loop_id in m.surfaces[tag],signed in m.loops[loop_id]
        curve=abs(signed)
        counts[curve]=get(counts,curve,0)+1
    end
    any(>(2),values(counts)) && throw(ArgumentError(
        "$caller: compound surfaces $members share a curve more than twice"))
    forced=Dict{Int,Vector{Float64}}()
    param_sizes=Dict{Tuple{Int,Float64},Float64}()
    for (key,specs) in m.meshing.size_at_params
        key[1]==1 || continue
        for (params,value) in specs
            for parameter in params
                param_sizes[(key[2],parameter)]=value
            end
        end
    end
    xs=Float64[];ys=Float64[];mesh_sizes=Float64[]
    boundary_segs=Tuple{Int,Int}[]
    internal_segs=Tuple{Int,Int}[]
    index=Dict{Int,Int}()
    member_polygons=Dict{Int,Vector{Vector{Int}}}()
    # Compound members must share one plane; the first member's plane is the
    # reference every other point is checked against.
    plane=_model_surface_plane(m,first(members),caller;allow_ruled=true)
    for tag in members
        member_polygons[tag]=Vector{Int}[]
        _surface_type(m,tag) in (:plane,:ruled,:tric) ||
            _model_require_plane_surface(m,tag,caller,
                                         "compound surface meshing")
        for loop_id in m.surfaces[tag]
            _verify_loop_closed(m,loop_id,caller,"Surface[$tag]")
            loop_idx=Int[]
            for signed in m.loops[loop_id]
                curve=abs(signed);a,b=m.curves[curve]
                _model_require_line_curve(m,curve,caller,
                                          "compound surface meshing")
                cspec=get(m.meshing.transfinite_curves,curve,nothing)
                (cspec===nothing || curve in m.meshing.degenerated) ||
                    (forced[curve]=collect(_transfinite_parameters(
                        m,cspec.num_nodes,cspec.kind,cspec.coef,
                        caller,curve;reversed=cspec.reversed)))
                loop_segment=Int[]
                for parameter in _surface_curve_parameters(
                        m,forced,curve,signed)
                    vertex=if parameter==0
                        _add_surface_point!(
                            xs,ys,mesh_sizes,index,m,a,caller,plane)
                    elseif parameter==1
                        _add_surface_point!(
                            xs,ys,mesh_sizes,index,m,b,caller,plane)
                    else
                        point=_periodic_curve_point(
                            m,curve,parameter,caller)
                        scale=max(1.0,hypot(point...))
                        abs(_plane_offset(plane,point))<=1e-12*scale ||
                            throw(ArgumentError(
                                "$caller: compound Curve[$curve] subdivision " *
                                "is not coplanar with Surface[" *
                                "$(first(members))]"))
                        _add_surface_curve_point!(
                            xs,ys,mesh_sizes,index,m,point,
                            get(param_sizes,(curve,parameter),
                                _surface_curve_mesh_size(
                                    m,curve,parameter,caller)),
                            curve,caller,plane)
                    end
                    push!(loop_idx,vertex)
                    push!(loop_segment,vertex)
                end
                for k in 1:(length(loop_segment)-1)
                    pair=(loop_segment[k],loop_segment[k+1])
                    pair[1]==pair[2] && throw(ArgumentError(
                        "$caller: compound Loop[$loop_id] has coincident " *
                        "consecutive vertices"))
                    if get(counts,curve,0)>1
                        push!(internal_segs,pair)
                    else
                        push!(boundary_segs,pair)
                    end
                end
            end
            push!(member_polygons[tag],loop_idx)
        end
    end
    _fill_unsized_surface_vertices!(
        mesh_sizes,xs,ys,boundary_segs,internal_segs)
    T=constrained_delaunay(
        xs,ys,boundary_segs;internal_segments=internal_segs)
    base=_surface_point_size_field(
        T,xs,ys,mesh_sizes,first(members),caller)
    callback=m.meshing.size_callback
    sizefn=if size_field===nothing && callback===nothing
        base
    else
        function sized(x,y)
            h=base(x,y)
            size_field===nothing ||
                (h=min(h,size_at(size_field,x,y,0.0,
                                 (2,first(members)))))
            callback===nothing ||
                (h=_apply_size_callback(
                    callback,2,first(members),x,y,0.0,h,caller))
            return h
        end
    end
    interior=refine!(T;min_angle_deg=25.0,size=sizefn)
    mesh=to_mesh(T;interior=interior)
    diag=validate(mesh)
    diag.ok || throw(ErrorException(
        "$caller: invalid compound mesh — "*join(diag.messages,"; ")))
    ntris(mesh)>0 || throw(ErrorException(
        "$caller: compound surfaces $members produced no triangles"))
    # Attribute each triangle to the member surface whose boundary polygon
    # contains its centroid.
    member_tris=Dict{Int,Vector{Int}}(tag=>Int[] for tag in members)
    for cell in axes(mesh.tris,2)
        cx=(mesh.coords[1,mesh.tris[1,cell]]+mesh.coords[1,mesh.tris[2,cell]]+
            mesh.coords[1,mesh.tris[3,cell]])/3
        cy=(mesh.coords[2,mesh.tris[1,cell]]+mesh.coords[2,mesh.tris[2,cell]]+
            mesh.coords[2,mesh.tris[3,cell]])/3
        owner=0
        for tag in members
            inside=false
            for polygon in member_polygons[tag]
                _point_in_polygon(cx,cy,xs,ys,polygon) && (inside=!inside)
            end
            inside && (owner=tag;break)
        end
        owner==0 && throw(ErrorException(
            "$caller: compound triangle $cell is outside every member " *
            "surface polygon"))
        push!(member_tris[owner],cell)
    end
    parts=Tuple{Int,Mesh}[]
    for tag in members
        cells=member_tris[tag]
        isempty(cells) && throw(ErrorException(
            "$caller: compound member Surface[$tag] received no triangles"))
        tris=Matrix{Int32}(mesh.tris[:,cells])
        referenced=falses(nnodes(mesh))
        @inbounds for column in axes(tris,2),slot in 1:3
            referenced[tris[slot,column]]=true
        end
        remap=zeros(Int32,nnodes(mesh))
        next=Int32(0)
        for node in 1:nnodes(mesh)
            referenced[node] || continue
            next+=Int32(1);remap[node]=next
        end
        @inbounds for column in axes(tris,2),slot in 1:3
            tris[slot,column]=remap[tris[slot,column]]
        end
        coordinates=Matrix{Float64}(undef,3,Int(next))
        for node in 1:nnodes(mesh)
            referenced[node] || continue
            coordinates[:,remap[node]].=@view mesh.coords[:,node]
        end
        part=Mesh(coordinates;tris=tris)
        # The member's boundary — including the shared compound seams — is
        # recomputed from its own triangles so downstream classification sees
        # the member's full boundary chain.
        _,member_edges=_surface_boundary_topology(part,caller)
        segs=Matrix{Int32}(undef,2,length(member_edges))
        for (seg,(a,b)) in enumerate(sort!(collect(member_edges)))
            segs[1,seg]=a;segs[2,seg]=b
        end
        part=Mesh(coordinates;segs=segs,tris=part.tris)
        push!(parts,(tag,part))
    end
    return parts
end

# Ray-cast point-in-polygon over the PSLG vertex list `polygon` (indices into
# `xs`/`ys`).
function _point_in_polygon(x,y,xs,ys,polygon)
    inside=false
    n=length(polygon)
    j=n
    for i in 1:n
        xi=xs[polygon[i]];yi=ys[polygon[i]]
        xj=xs[polygon[j]];yj=ys[polygon[j]]
        if (yi>y)!=(yj>y) && x<(xj-xi)*(y-yi)/(yj-yi)+xi
            inside=!inside
        end
        j=i
    end
    return inside
end


end # module
