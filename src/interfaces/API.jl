"""
    API

Gmsh-style model/mesh/option façade over Tessella's native kernels, including
deterministic entity-topology, analytical spatial, native type/partition metadata,
native Point/Line/Plane evaluation and reparametrization, entity presentation and
attribute state, finite Point-coordinate updates, and Physical-group queries,
entity-name, atomic-tag, and dependency-safe removal mutations, Physical-group
mutations, point-local mesh-size constraints, explicit planar surface-loop volumes,
operation-time Boolean operand ownership, and persistent affine relations between
straight periodic boundary or embedded curves and planar periodic volume boundaries.
The session owns atomic uniform refinement, affine coordinate transformation,
complete clearing, and detached Gmsh-shaped bulk node/element retrieval for its
linear-simplex mesh cache, plus deterministic global edge and triangular or
quadrangular face catalogs, fixed-family actual- and explicit-order nodal
reference functions, order-one H1 bases over simplex, Point, Quadrangle,
Hexahedron, and Prism reference families with lowest-order H(curl) bases over
simplexes, orientations,
and node/edge keys. It also
owns a reusable robust AABB locator for
dense element-by-coordinate and reference-coordinate queries, plus scale-robust
named quality queries and Gmsh-shaped forward maps/Jacobians over dense cached
elements. Element type/property lookup and bounded fixed-family reference
quadrature are available without a session. Production meshing is never delegated
to Gmsh.
"""
module API

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_surface_loop!, add_volume!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
using ..Model: _remove_volume_entity!
using ..Model: embed!, remove_embedded!, set_point_mesh_size!
using ..Model: add_physical_group!, set_physical_name!, remove_physical_groups!
using ..Model: remove_physical_name!, model_physical_groups
using ..Model: model_physical_groups_entities, model_entities_for_physical_group
using ..Model: model_entities_for_physical_name, model_physical_groups_for_entity
using ..Model: model_physical_name
using ..Model: model_entities, model_dimension, model_boundary, model_adjacencies
using ..Model: model_bounding_box, model_entities_in_bounding_box
using ..Model: model_entity_type, model_entity_properties, model_parent,
               model_number_of_partitions, model_partitions
using ..Model: model_value, model_derivative, model_second_derivative,
               model_curvature, model_principal_curvatures, model_normal
using ..Model: model_parametrization, model_parametrization_bounds,
               model_is_inside, model_closest_point,
               model_reparametrize_on_surface
using ..Model: set_entity_visibility!, model_entity_visibility
using ..Model: set_entity_color!, model_entity_color, set_point_coordinates!
using ..Model: set_model_attribute!, model_attribute, model_attribute_names,
               remove_model_attribute!
using ..Model: set_entity_name!, remove_entity_name!, model_entity_name, model_set_tag!
using ..Model: remove_entities!
using ..Model: set_periodic!, model_periodic_nodes
using ..Model: mesh_model_surface, mesh_model_volume, model_to_mixed
using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets
using ..MeshEntityTopology: MeshEdgeTopology, MeshFaceTopology,
                            _mesh_edge_topology, _mesh_face_topology,
                            _mesh_edge_topology_for_cells,
                            _mesh_face_topology_for_cells,
                            _mesh_add_edges, _mesh_add_faces,
                            _mesh_edges, _mesh_faces,
                            _mesh_all_edges, _mesh_all_faces,
                            _simplex_edge_patterns, _simplex_face_patterns
using ..MeshPointLocation: SimplexLocator, mesh_element_offsets,
                           mesh_element_block, mesh_element_record,
                           _local_coordinates, _locate_elements,
                           _require_local_coordinates
using ..MeshElementQuality: mesh_element_qualities
using ..MeshQuadrature: mesh_integration_points
using ..MeshReferenceGeometry: mesh_jacobian, mesh_jacobians
using ..MeshFunctionSpaces: MeshFunctionSpaces, mesh_basis_functions,
                            mesh_basis_orientation,
                            mesh_basis_orientations, mesh_key_dimension,
                            mesh_keys, mesh_keys_for_element,
                            mesh_keys_information, mesh_number_of_keys,
                            mesh_number_of_orientations
using ..Elements: msh_spec, msh_type, msh_properties
using ..Refine: refine_uniform
using ..Transform: affine_transform, _transform_gmsh_affine
using ..Optimize: smooth_optimize
using ..GeoExec: execute_geo

export initialize, finalize, option, model, mesh, open_geo!

# Entity classification snapshot for the session mesh cache, derived from the
# canonical `model_to_mixed` projection at generation time. `mesh` is the exact
# cache object the record describes so a cache replaced outside
# `_replace_mesh_cache_locked!` cannot accidentally reuse a stale record.
# `node_entities[i]` is the lowest-dimension entity owning node `i`;
# `seg/tri/tet_entities[c]` is the entity tag owning cell column `c`;
# `boundaries[(dim,tag)]` lists the entity's boundary tags one dimension down.
struct _MeshClassification
    mesh::Mesh
    entity::Tuple{Int,Int32}
    node_entities::Vector{Tuple{Int,Int32}}
    boundaries::Dict{Tuple{Int,Int32},Vector{Int32}}
    seg_entities::Vector{Int32}
    tri_entities::Vector{Int32}
    tet_entities::Vector{Int32}
end

const CURRENT = Ref{Union{Nothing,GeoModel}}(nothing)
const DEFAULT_OPTIONS = Dict{String,Float64}(
    "Mesh.MeshSizeMin"=>0.0,
    "Mesh.MeshSizeMax"=>1.0e22,
    "Mesh.MeshSizeFactor"=>1.0)
const OPTIONS = copy(DEFAULT_OPTIONS)
const LAST_MESH = Ref{Union{Nothing,Mesh}}(nothing)
const LAST_MESH_CLASS = Ref{Union{Nothing,_MeshClassification}}(nothing)
const LAST_MESH_LOCATOR = Ref{Union{Nothing,SimplexLocator}}(nothing)
const LAST_MESH_EDGES = Ref{Union{Nothing,MeshEdgeTopology}}(nothing)
const LAST_MESH_FACES = Ref{Union{Nothing,MeshFaceTopology}}(nothing)
const MODEL_NAME = Ref{String}("unnamed")
const MODEL_FILE_NAME = Ref{String}("")
const ELEMENT_VISIBILITY = Ref{Dict{Int,Int32}}(Dict{Int,Int32}())
const STATE_LOCK = ReentrantLock()

function _replace_mesh_cache_locked!(mesh::Union{Nothing,Mesh},
                                     class::Union{Nothing,_MeshClassification}=
                                         nothing)
    class!==nothing && class.mesh!==mesh && throw(ArgumentError(
        "API: internal classification record does not match the cached mesh"))
    LAST_MESH[]=mesh
    LAST_MESH_CLASS[]=class
    LAST_MESH_LOCATOR[]=nothing
    LAST_MESH_EDGES[]=nothing
    LAST_MESH_FACES[]=nothing
    empty!(ELEMENT_VISIBILITY[])
    return mesh
end

function _cached_classification_locked(mesh::Mesh)
    class=LAST_MESH_CLASS[]
    (class===nothing || class.mesh!==mesh) && return nothing
    return class
end

"""
    initialize()

Start a fresh process-global API session, discarding any prior model or cached
mesh and restoring supported options to their defaults.
"""
function initialize()
    lock(STATE_LOCK) do
        CURRENT[]=GeoModel()
        MODEL_NAME[]="unnamed"
        MODEL_FILE_NAME[]=""
        _replace_mesh_cache_locked!(nothing)
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
    end
    return nothing
end

"""End the API session, discard its model and cached mesh, and restore option defaults."""
function finalize()
    lock(STATE_LOCK) do
        CURRENT[]=nothing
        _replace_mesh_cache_locked!(nothing)
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
    end
    return nothing
end

function _model_locked()
    CURRENT[]===nothing && throw(ArgumentError(
        "API: no current model; call initialize() or model.add(name) first"))
    return CURRENT[]
end

function _with_model(f::Function;invalidate::Bool=false)
    return lock(STATE_LOCK) do
        current=_model_locked()
        result=f(current)
        invalidate && _replace_mesh_cache_locked!(nothing)
        result
    end
end

function _copy_mesh(mesh::Mesh)
    return Mesh(mesh.coords;segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
end

"""
    option(name) -> Float64
    option(name, value) -> Float64

Get or set one supported process-global mesh option in an initialized session.
`MeshSizeMin` is nonnegative; `MeshSizeMax` and `MeshSizeFactor` are positive;
the minimum may not exceed the maximum. Boolean and nonfinite values are
rejected, and a failed update leaves all options unchanged.
"""
function option(name::AbstractString)
    key=String(name)
    return lock(STATE_LOCK) do
        _model_locked()
        haskey(OPTIONS,key) || throw(ArgumentError("API.option: unknown option $name"))
        OPTIONS[key]
    end
end
function option(name::AbstractString, value::Real)
    key=String(name)
    return lock(STATE_LOCK) do
        _model_locked()
        haskey(OPTIONS,key) || throw(ArgumentError("API.option: unknown option $name"))
        value isa Bool && throw(ArgumentError("API.option: value must not be Bool"))
        v=try
            Float64(value)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("API.option: value must be Float64-representable"))
        end
        isfinite(v) || throw(ArgumentError("API.option: value must be finite"))
        if key=="Mesh.MeshSizeMin"
            v>=0 || throw(ArgumentError("API.option: MeshSizeMin must be nonnegative"))
            v<=OPTIONS["Mesh.MeshSizeMax"] || throw(ArgumentError(
                "API.option: MeshSizeMin must not exceed MeshSizeMax"))
        elseif key=="Mesh.MeshSizeMax"
            v>0 || throw(ArgumentError("API.option: MeshSizeMax must be positive"))
            v>=OPTIONS["Mesh.MeshSizeMin"] || throw(ArgumentError(
                "API.option: MeshSizeMax must not be below MeshSizeMin"))
        else
            v>0 || throw(ArgumentError("API.option: MeshSizeFactor must be positive"))
        end
        OPTIONS[key]=v
        v
    end
end

_get_physical_groups(dim=-1)=_with_model() do current
    model_physical_groups(current,dim)
end

_get_entities(dim=-1)=_with_model() do current
    model_entities(current,dim)
end

_get_dimension()=_with_model() do current
    model_dimension(current)
end

# Single-model session state: `get_current` reports the session model name and
# `set_current` only accepts that name — there is no second model to switch
# to. `get_file_name`/`set_file_name` track the model's associated file.
function _model_add(name)
    caller="API.model.add"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    return lock(STATE_LOCK) do
        if CURRENT[]!==nothing
            # `initialize` eagerly creates an unnamed empty model, so the
            # universal `initialize(); model.add("m")` sequence names that
            # fresh model. Any existing content — or a second `add` — is a
            # genuine multi-model request the session cannot serve.
            fresh=isempty(model_entities(CURRENT[])) &&
                  isempty(CURRENT[].physical) &&
                  MODEL_NAME[]=="unnamed"
            fresh || throw(ArgumentError(
                "$caller: the session already owns model " *
                "\"$(MODEL_NAME[])\"; a second model is not supported"))
        else
            CURRENT[]=GeoModel()
        end
        MODEL_NAME[]=String(name)
        MODEL_FILE_NAME[]=""
        _replace_mesh_cache_locked!(nothing)
        return nothing
    end
end

function _model_remove()
    return lock(STATE_LOCK) do
        _model_locked()
        CURRENT[]=nothing
        MODEL_NAME[]="unnamed"
        MODEL_FILE_NAME[]=""
        _replace_mesh_cache_locked!(nothing)
        return nothing
    end
end

function _model_list()
    return lock(STATE_LOCK) do
        return CURRENT[]===nothing ? String[] : String[MODEL_NAME[]]
    end
end

function _get_current()
    return lock(STATE_LOCK) do
        _model_locked()
        return MODEL_NAME[]
    end
end

function _set_current(name)
    caller="API.model.set_current"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    return lock(STATE_LOCK) do
        _model_locked()
        String(name)==MODEL_NAME[] || throw(ArgumentError(
            "$caller: unknown model \"$name\"; the session owns model " *
            "\"$(MODEL_NAME[])\""))
        return nothing
    end
end

function _get_file_name()
    return lock(STATE_LOCK) do
        _model_locked()
        return MODEL_FILE_NAME[]
    end
end

function _set_file_name(name)
    caller="API.model.set_file_name"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    return lock(STATE_LOCK) do
        _model_locked()
        MODEL_FILE_NAME[]=String(name)
        return nothing
    end
end

_get_boundary(dim_tags,combined=true,oriented=false,recursive=false)=
    _with_model() do current
        model_boundary(
            current,dim_tags,combined,oriented,recursive)
    end

_get_adjacencies(dim,tag)=_with_model() do current
    model_adjacencies(current,dim,tag)
end

# `isEntityOrphan` connectivity is the transitive downward boundary closure of
# every entity at the model's highest dimension; embeddings do not connect.
# Classified entities (including implicit primitive or Boolean subentities)
# expand through the classification boundary map; anything else expands through
# explicit model topology, where primitive volumes contribute no children.
function _orphan_boundary_children(model,class,key::Tuple{Int,Int32})
    key[1]==0 && return Int32[]
    if class!==nothing && haskey(class.boundaries,key)
        return class.boundaries[key]
    end
    dictionary=_mesh_entity_dictionary(model,key[1])
    haskey(dictionary,Int(key[2])) || return Int32[]
    (key[1]==3 && isempty(dictionary[Int(key[2])])) && return Int32[]
    return Int32[abs(boundary) for (_,boundary) in model_boundary(
        model,[(key[1],Int(key[2]))],false,false,false)]
end

function _is_entity_orphan(dim,tag)
    caller="API.model.is_entity_orphan"
    return lock(STATE_LOCK) do
        model=_model_locked()
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        entity=_mesh_query_integer(tag,caller,"tag")
        1<=entity<=typemax(Int32) || throw(ArgumentError(
            "$caller: unknown $(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
        cached=LAST_MESH[]
        class=cached===nothing ? nothing :
              _cached_classification_locked(cached)
        known=haskey(_mesh_entity_dictionary(model,dimension),entity) ||
              (class!==nothing &&
               haskey(class.boundaries,(dimension,Int32(entity))))
        known || throw(ArgumentError(
            "$caller: unknown $(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
        highest=model_dimension(model)
        connected=Set{Tuple{Int,Int32}}()
        queue=Tuple{Int,Int32}[(highest,Int32(top)) for top in
            keys(_mesh_entity_dictionary(model,highest))]
        if class!==nothing
            for key in keys(class.boundaries)
                key[1]==highest && push!(queue,key)
            end
        end
        while !isempty(queue)
            key=popfirst!(queue)
            key in connected && continue
            push!(connected,key)
            for child in _orphan_boundary_children(model,class,key)
                push!(queue,(key[1]-1,child))
            end
        end
        return !((dimension,Int32(entity)) in connected)
    end
end

_get_bounding_box(dim,tag)=_with_model() do current
    model_bounding_box(current,dim,tag)
end

_get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,dim=-1)=
    _with_model() do current
        model_entities_in_bounding_box(
            current,xmin,ymin,zmin,xmax,ymax,zmax,dim)
    end

_get_entity_type(dim,tag)=_with_model() do current
    model_entity_type(current,dim,tag)
end

_get_entity_properties(dim,tag)=_with_model() do current
    model_entity_properties(current,dim,tag)
end

_get_value(dim,tag,parametric_coord)=_with_model() do current
    model_value(current,dim,tag,parametric_coord)
end

_get_derivative(dim,tag,parametric_coord)=_with_model() do current
    model_derivative(current,dim,tag,parametric_coord)
end

_get_second_derivative(dim,tag,parametric_coord)=_with_model() do current
    model_second_derivative(current,dim,tag,parametric_coord)
end

_get_curvature(dim,tag,parametric_coord)=_with_model() do current
    model_curvature(current,dim,tag,parametric_coord)
end

_get_principal_curvatures(tag,parametric_coord)=_with_model() do current
    model_principal_curvatures(current,tag,parametric_coord)
end

_get_normal(tag,parametric_coord)=_with_model() do current
    model_normal(current,tag,parametric_coord)
end

_get_parametrization(dim,tag,coord)=_with_model() do current
    model_parametrization(current,dim,tag,coord)
end

_get_parametrization_bounds(dim,tag)=_with_model() do current
    model_parametrization_bounds(current,dim,tag)
end

_is_inside(dim,tag,coord,parametric=false)=_with_model() do current
    model_is_inside(current,dim,tag,coord,parametric)
end

_get_closest_point(dim,tag,coord)=_with_model() do current
    model_closest_point(current,dim,tag,coord)
end

_reparametrize_on_surface(dim,tag,parametric_coord,surface_tag,which=0)=
    _with_model() do current
        model_reparametrize_on_surface(
            current,dim,tag,parametric_coord,surface_tag,which)
    end

_set_visibility(dim_tags,value,recursive=false)=_with_model() do current
    set_entity_visibility!(current,dim_tags,value,recursive)
end

_get_visibility(dim,tag)=_with_model() do current
    model_entity_visibility(current,dim,tag)
end

_set_color(dim_tags,r,g,b,a=255,recursive=false)=_with_model() do current
    set_entity_color!(current,dim_tags,r,g,b,a,recursive)
end

_get_color(dim,tag)=_with_model() do current
    model_entity_color(current,dim,tag)
end

_set_coordinates(tag,x,y,z)=_with_model(invalidate=true) do current
    set_point_coordinates!(current,tag,x,y,z)
end

_set_attribute(name,values)=_with_model() do current
    set_model_attribute!(current,name,values)
end

_get_attribute(name)=_with_model() do current
    model_attribute(current,name)
end

_get_attribute_names()=_with_model() do current
    model_attribute_names(current)
end

_remove_attribute(name)=_with_model() do current
    remove_model_attribute!(current,name)
    nothing
end

_get_parent(dim,tag)=_with_model() do current
    model_parent(current,dim,tag)
end

_get_number_of_partitions()=_with_model() do current
    model_number_of_partitions(current)
end

_get_partitions(dim,tag)=_with_model() do current
    model_partitions(current,dim,tag)
end

_get_entity_name(dim,tag)=_with_model() do current
    model_entity_name(current,dim,tag)
end

function _set_entity_name(dim,tag,name)
    name isa AbstractString || throw(ArgumentError(
        "API.model.set_entity_name: name must be a string"))
    return lock(STATE_LOCK) do
        current=_model_locked()
        before=model_entity_name(current,dim,tag)
        after=set_entity_name!(current,dim,tag,name)
        after==before || _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

function _remove_entity_name(name)
    name isa AbstractString || throw(ArgumentError(
        "API.model.remove_entity_name: name must be a string"))
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_entity_name!(current,name)>0 && _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

_set_tag(dim,tag,new_tag)=_with_model(invalidate=true) do current
    model_set_tag!(current,dim,tag,new_tag)
    nothing
end

function _remove_entities(dim_tags,recursive=false)
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_entities!(current,dim_tags,recursive)>0 &&
            _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

_get_physical_groups_entities(dim=-1)=_with_model() do current
    model_physical_groups_entities(current,dim)
end

_get_entities_for_physical_group(dim,tag)=_with_model() do current
    model_entities_for_physical_group(current,dim,tag)
end

_get_entities_for_physical_name(name)=_with_model() do current
    name isa AbstractString || throw(ArgumentError(
        "API.model.get_entities_for_physical_name: name must be a string"))
    model_entities_for_physical_name(current,name)
end

_get_physical_groups_for_entity(dim,tag)=_with_model() do current
    model_physical_groups_for_entity(current,dim,tag)
end

_get_physical_name(dim,tag)=_with_model() do current
    model_physical_name(current,dim,tag)
end

function _set_physical_name(dim,tag,name)
    name isa AbstractString || throw(ArgumentError(
        "API.model.set_physical_name: name must be a string"))
    return lock(STATE_LOCK) do
        current=_model_locked()
        before=model_physical_name(current,dim,tag)
        after=set_physical_name!(current,dim,tag,name)
        after==before || _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

function _remove_physical_name(name)
    name isa AbstractString || throw(ArgumentError(
        "API.model.remove_physical_name: name must be a string"))
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_physical_name!(current,name)>0 &&
            _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

function _remove_physical_groups(dim_tags=())
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_physical_groups!(current,dim_tags)>0 &&
            _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

"""Gmsh-style geometry-model, topology, and Physical-group operations for the active [`API`](@ref) session."""
module model
using ..API: _with_model, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..API: add_surface_loop!, add_volume!
using ..API: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!, embed!, add_physical_group!
using ..API: _remove_volume_entity!
using ..API: _get_physical_groups, _get_physical_groups_entities
using ..API: _get_entities_for_physical_group, _get_entities_for_physical_name
using ..API: _get_physical_groups_for_entity, _get_physical_name
using ..API: _set_physical_name, _remove_physical_name, _remove_physical_groups
using ..API: _get_entities, _get_dimension, _get_boundary, _get_adjacencies,
             _is_entity_orphan
using ..API: _get_bounding_box, _get_entities_in_bounding_box
using ..API: _get_entity_type, _get_entity_properties, _get_parent,
             _get_number_of_partitions, _get_partitions
using ..API: _get_value, _get_derivative, _get_second_derivative, _get_curvature
using ..API: _get_principal_curvatures, _get_normal, _get_parametrization,
             _get_parametrization_bounds, _is_inside, _get_closest_point
using ..API: _reparametrize_on_surface
using ..API: _set_visibility, _get_visibility, _set_color, _get_color
using ..API: _set_coordinates, _set_attribute, _get_attribute,
             _get_attribute_names, _remove_attribute
using ..API: _get_entity_name, _set_entity_name, _remove_entity_name, _set_tag
using ..API: _remove_entities
using ..API: _get_current, _set_current, _get_file_name, _set_file_name,
             _model_add, _model_remove, _model_list
add_point(x,y,z;tag=0,meshSize=1.0)=_with_model(invalidate=true) do m
    add_point!(m,x,y,z;tag=tag,mesh_size=meshSize)
end
add_line(a,b;tag=0)=_with_model(invalidate=true) do m
    add_line!(m,a,b;tag=tag)
end
add_curve_loop(curves;tag=0)=_with_model(invalidate=true) do m
    add_curve_loop!(m,curves;tag=tag)
end
add_plane_surface(loops;tag=0)=_with_model(invalidate=true) do m
    add_plane_surface!(m,loops;tag=tag)
end
"""Add a connected closed planar-surface shell to the active model."""
add_surface_loop(surfaces;tag=0)=_with_model(invalidate=true) do m
    add_surface_loop!(m,surfaces;tag=tag)
end
"""Add a volume from an exterior surface loop and optional cavity loops."""
add_volume(surface_loops;tag=0)=_with_model(invalidate=true) do m
    add_volume!(m,surface_loops;tag=tag)
end
add_box(x,y,z,dx,dy,dz;tag=0)=_with_model(invalidate=true) do m
    add_box!(m,x,y,z,dx,dy,dz;tag=tag)
end
add_cylinder(x,y,z,dx,dy,dz,r;tag=0)=_with_model(invalidate=true) do m
    add_cylinder!(m,x,y,z,dx,dy,dz,r;tag=tag)
end
add_sphere(x,y,z,r;tag=0)=_with_model(invalidate=true) do m
    add_sphere!(m,x,y,z,r;tag=tag)
end
add_cone(x,y,z,dx,dy,dz,r1,r2;tag=0)=_with_model(invalidate=true) do m
    add_cone!(m,x,y,z,dx,dy,dz,r1,r2;tag=tag)
end
embed(dim,tags,target_dim,target_tag)=_with_model(invalidate=true) do m
    embed!(m,dim,tags,target_dim,target_tag)
end
function _boolean(op,a,b; tag=0)
    return _with_model(invalidate=true) do m
        t=boolean_volumes!(m,op,a,b;tag=tag)
        # OCC cut/fuse/common remove object and tool by default.
        _remove_volume_entity!(m,Int(a))
        _remove_volume_entity!(m,Int(b))
        t
    end
end
boolean_difference(a,b; tag=0)=_boolean(:difference,a,b; tag=tag)
boolean_union(a,b; tag=0)=_boolean(:union,a,b; tag=tag)
boolean_intersection(a,b; tag=0)=_boolean(:intersection,a,b; tag=tag)

"""Return an existing entity's name, or an empty string when it is unnamed or missing."""
get_entity_name(dim,tag)=_get_entity_name(dim,tag)

"""
    set_entity_name(dim, tag, name)

Set or replace the name of an existing positive-tag entity. An empty name removes
the current name; a missing entity is unchanged.
"""
set_entity_name(dim,tag,name)=_set_entity_name(dim,tag,name)

"""Remove `name` from every model entity carrying it."""
remove_entity_name(name)=_remove_entity_name(name)

"""
    set_tag(dim, tag, new_tag)

Atomically move an existing entity to an unused positive tag in the same dimension,
including all live model references, its name, visibility, and color.
"""
set_tag(dim,tag,new_tag)=_set_tag(dim,tag,new_tag)

"""
    remove_entities(dim_tags, recursive=false)

Remove ordered positive `(dimension, entity_tag)` pairs when they are not used by a
surviving boundary or embedding owner. With `recursive=true`, also process explicit
boundary entities down to Points. Missing and still-owned entities are unchanged;
removed presentation state is cleaned.
"""
remove_entities(dim_tags,recursive=false)=
    _remove_entities(dim_tags,recursive)

"""Add a Physical group; `tag=0` allocates globally across entity dimensions."""
add_physical_group(dim,tags;tag=0,name="")=_with_model(invalidate=true) do m
    add_physical_group!(m,dim,tags;tag=tag,name=name)
end

"""Return detached, sorted native model `(dimension, tag)` pairs."""
get_entities(dim=-1)=_get_entities(dim)

"""Return the greatest model-entity dimension, or `-1` for an empty model."""
get_dimension()=_get_dimension()

"""Return the name of the session's current model."""
get_current()=_get_current()

"""
    add(name)

Name the session's model. `initialize` eagerly creates an unnamed empty model,
so a single `add` on that fresh model names it — matching Gmsh's
`model.add(name)` — and `add` after `remove` creates a new model. A second
model on a populated or already-named model fails explicitly.
"""
add(name)=_model_add(name)

"""Remove the session's current model, leaving the session model-less."""
remove()=_model_remove()

"""Return the session's model names — `[get_current()]`, or `[]` after `remove`."""
list()=_model_list()

"""
    set_current(name)

Make `name` the current model. The session owns a single model, so only its own
name — reported by `get_current` — is accepted; anything else throws.
"""
set_current(name)=_set_current(name)

"""Return the file name associated with the current model, or `""` when unset."""
get_file_name()=_get_file_name()

"""Associate `name` with the current model as its file name."""
set_file_name(name)=_set_file_name(name)

"""
Return explicit entity boundaries. `combined` cancels even incidences, `oriented`
retains signed Curve and Surface tags, and `recursive` returns the Point closure.
Primitive and Boolean Volume boundaries are unavailable because their Surface Loop
topology is implicit.
"""
get_boundary(dim_tags,combined=true,oriented=false,recursive=false)=
    _get_boundary(dim_tags,combined,oriented,recursive)

"""
Return detached upward and downward topology-only adjacency tags. Downward queries
for primitive and Boolean Volumes are unavailable because their Surface Loop
topology is implicit.
"""
get_adjacencies(dim,tag)=_get_adjacencies(dim,tag)

"""
    is_entity_orphan(dim, tag) -> Bool

Return `true` when the entity is not connected to any entity of the model's
highest dimension, matching Gmsh 4.15.2's `isEntityOrphan`: connectivity is the
transitive downward boundary closure of every highest-dimension entity, and
mesh embeddings do not connect. Implicit primitive or Boolean boundary entities
resolve through the mesh classification snapshot. Unknown entities fail
explicitly.
"""
is_entity_orphan(dim,tag)=_is_entity_orphan(dim,tag)

"""Return one entity's analytical bounding box; `(-1,-1)` selects the whole model."""
get_bounding_box(dim,tag)=_get_bounding_box(dim,tag)

"""
Return detached, sorted entities whose complete analytical bounding box is inside
the supplied finite box. `dim=-1` selects all dimensions.
"""
get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,dim=-1)=
    _get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,dim)

"""Return the native type of an existing Point, Line, Plane, or Volume entity."""
get_entity_type(dim,tag)=_get_entity_type(dim,tag)

"""Gmsh-compatible synonym for [`get_entity_type`](@ref)."""
get_type(dim,tag)=_get_entity_type(dim,tag)

"""
Return detached native entity-property vectors. Plane reals are `[a,b,c,d]` for
the unit-normal equation `a*x+b*y+c*z=d`; other visible native types return empty
vectors.
"""
get_entity_properties(dim,tag)=_get_entity_properties(dim,tag)

"""Evaluate an explicit Point, straight Line, or Plane parametrization."""
get_value(dim,tag,parametric_coord)=_get_value(dim,tag,parametric_coord)

"""Evaluate first derivatives for an explicit straight Line or Plane."""
get_derivative(dim,tag,parametric_coord)=
    _get_derivative(dim,tag,parametric_coord)

"""Evaluate second derivatives for an explicit straight Line or Plane."""
get_second_derivative(dim,tag,parametric_coord)=
    _get_second_derivative(dim,tag,parametric_coord)

"""Return zero curvature values for an explicit straight Line or Plane."""
get_curvature(dim,tag,parametric_coord)=
    _get_curvature(dim,tag,parametric_coord)

"""Return zero principal curvatures and tangent directions for an explicit Plane."""
get_principal_curvatures(tag,parametric_coord)=
    _get_principal_curvatures(tag,parametric_coord)

"""Return exterior-loop-oriented unit normals for an explicit Plane."""
get_normal(tag,parametric_coord)=_get_normal(tag,parametric_coord)

"""Return orthogonal native parameters for concatenated 3-D coordinates."""
get_parametrization(dim,tag,coord)=_get_parametrization(dim,tag,coord)

"""Return detached native parametric lower and upper bounds."""
get_parametrization_bounds(dim,tag)=
    _get_parametrization_bounds(dim,tag)

"""Count physical or parametric points inside an explicit native entity."""
is_inside(dim,tag,coord,parametric=false)=
    _is_inside(dim,tag,coord,parametric)

"""Project concatenated 3-D coordinates onto an explicit Line or Plane."""
get_closest_point(dim,tag,coord)=_get_closest_point(dim,tag,coord)

"""Map Point or straight-Line parameters into an explicit Plane's parameters."""
reparametrize_on_surface(dim,tag,parametric_coord,surface_tag,which=0)=
    _reparametrize_on_surface(dim,tag,parametric_coord,surface_tag,which)

"""Set an `Int32` visibility value on existing native entities."""
set_visibility(dim_tags,value,recursive=false)=
    _set_visibility(dim_tags,value,recursive)

"""Return an existing native entity's visibility, which defaults to `1`."""
get_visibility(dim,tag)=_get_visibility(dim,tag)

"""Set an RGBA color on existing entities, optionally including their boundaries."""
set_color(dim_tags,r,g,b,a=255,recursive=false)=
    _set_color(dim_tags,r,g,b,a,recursive)

"""Return an existing entity's RGBA color, defaulting to `(0,0,255,0)`."""
get_color(dim,tag)=_get_color(dim,tag)

"""Replace one existing explicit Point's finite coordinates and invalidate the mesh."""
set_coordinates(tag,x,y,z)=_set_coordinates(tag,x,y,z)

"""Set detached string values for a global model attribute."""
set_attribute(name,values)=_set_attribute(name,values)

"""Return a detached model-attribute value list, or an empty list if absent."""
get_attribute(name)=_get_attribute(name)

"""Return model-attribute names in deterministic lexical order."""
get_attribute_names()=_get_attribute_names()

"""Remove a model attribute; a missing name is unchanged."""
remove_attribute(name)=_remove_attribute(name)

"""Return `(-1,-1)`, the partition-parent sentinel, for an existing native entity."""
get_parent(dim,tag)=_get_parent(dim,tag)

"""Return zero because the native geometry model does not own mesh partitions."""
get_number_of_partitions()=_get_number_of_partitions()

"""Return an empty detached partition-membership list for an existing native entity."""
get_partitions(dim,tag)=_get_partitions(dim,tag)

"""Return detached, sorted Physical `(dimension, tag)` pairs; `dim=-1` selects all."""
get_physical_groups(dim=-1)=_get_physical_groups(dim)

"""Return sorted Physical groups and detached `(dimension, entity_tag)` memberships."""
get_physical_groups_entities(dim=-1)=_get_physical_groups_entities(dim)

"""Return detached, sorted entity tags for one existing Physical group."""
get_entities_for_physical_group(dim,tag)=
    _get_entities_for_physical_group(dim,tag)

"""Return detached, sorted entities belonging to groups with the given Physical name."""
get_entities_for_physical_name(name)=_get_entities_for_physical_name(name)

"""Return sorted Physical tags containing one existing geometry entity."""
get_physical_groups_for_entity(dim,tag)=
    _get_physical_groups_for_entity(dim,tag)

"""Return a Physical name, or `""` for an unnamed or missing positive group tag."""
get_physical_name(dim,tag)=_get_physical_name(dim,tag)

"""
Assign a name to an unnamed Physical group. Existing, empty, duplicate, and missing
assignments are no-ops; remove the current name before assigning a replacement.
"""
set_physical_name(dim,tag,name)=_set_physical_name(dim,tag,name)

"""Remove a Physical name from every dimension without removing its groups."""
remove_physical_name(name)=_remove_physical_name(name)

"""Remove selected Physical groups, or every group when `dim_tags` is empty."""
remove_physical_groups(dim_tags=())=_remove_physical_groups(dim_tags)
end

# Project `mesh` onto `entity` through `model_to_mixed` and invert the emitted
# block connectivity back into per-cache-cell entity tags. A projection failure
# never fails generation: the record simply stays absent and entity-filtered
# queries keep their explicit blocker.
function _classify_cached_mesh(m::GeoModel,mesh::Mesh,dim::Int,tag::Int,
                               cache::Mesh)
    classified=try
        dim==2 ? model_to_mixed(m,mesh,tag) : model_to_mixed(m,mesh,dim,tag)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    data=classified.entity_data
    data===nothing && return nothing
    length(data.node_entities)==nnodes(mesh) || return nothing
    boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}()
    for (key,entity) in pairs(data.entities)
        boundaries[(key[1],key[2])]=
            Int32[abs(boundary) for boundary in entity.boundaries]
    end
    owner=Dict{NTuple{4,Int32},Int32}()
    for (block_index,block) in enumerate(classified.blocks)
        block_entities=data.block_entities[block_index]
        for column in axes(block.nodes,2)
            connectivity=sort!(Int32.(vec(block.nodes[:,column])))
            padded=ntuple(
                slot->slot<=length(connectivity) ? connectivity[slot] :
                    Int32(0),4)
            haskey(owner,padded) && return nothing
            owner[padded]=block_entities[column]
        end
    end
    function cell_entities(cells::AbstractMatrix{Int32})
        result=Vector{Int32}(undef,size(cells,2))
        for column in axes(cells,2)
            connectivity=sort!(Int32.(vec(cells[:,column])))
            padded=ntuple(
                slot->slot<=length(connectivity) ? connectivity[slot] :
                    Int32(0),4)
            value=get(owner,padded,Int32(0))
            value==0 && return nothing
            result[column]=value
        end
        return result
    end
    seg_entities=cell_entities(mesh.segs)
    tri_entities=cell_entities(mesh.tris)
    tet_entities=cell_entities(mesh.tets)
    (seg_entities===nothing || tri_entities===nothing ||
     tet_entities===nothing) && return nothing
    return _MeshClassification(
        cache,(dim,Int32(tag)),copy(data.node_entities),boundaries,
        seg_entities,tri_entities,tet_entities)
end

function _generate(dim::Integer)
    dim isa Bool && throw(ArgumentError("API.mesh.generate: dim must not be Bool"))
    dimension=try
        Int(dim)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("API.mesh.generate: dim exceeds the platform Int range"))
    end
    return lock(STATE_LOCK) do
        m=_model_locked()
        generated=if dimension==2
            isempty(m.surfaces) && throw(ArgumentError("API.mesh.generate: no surfaces"))
            length(m.surfaces)==1 || throw(ArgumentError(
                "API.mesh.generate: Mesh 2 with multiple remaining surfaces is a blocker"))
            mesh_model_surface(m,only(keys(m.surfaces)))
        elseif dimension==3
            isempty(m.volumes) && throw(ArgumentError("API.mesh.generate: no volumes"))
            length(m.volumes)==1 || throw(ArgumentError(
                "API.mesh.generate: Mesh 3 with multiple remaining volumes is a blocker — Boolean Delete the operands or mesh a single volume"))
            mesh_model_volume(m,only(keys(m.volumes)))
        else
            throw(ArgumentError("API.mesh.generate: dim must be 2 or 3"))
        end
        entity_tag=Int(only(keys(dimension==2 ? m.surfaces : m.volumes)))
        cache=_copy_mesh(generated)
        class=_classify_cached_mesh(m,generated,dimension,entity_tag,cache)
        _replace_mesh_cache_locked!(cache,class)
        generated
    end
end

function _cached_mesh_locked(caller::AbstractString)
    _model_locked()
    cached=LAST_MESH[]
    cached===nothing && throw(ArgumentError(
        "$caller: no mesh; call API.mesh.generate first"))
    return cached
end

function _cached_mesh_locator_locked(mesh::Mesh)
    locator=LAST_MESH_LOCATOR[]
    if locator===nothing || locator.mesh!==mesh
        replacement=SimplexLocator(mesh)
        LAST_MESH_LOCATOR[]=replacement
        return replacement
    end
    return locator
end

function _get_mesh()
    return lock(STATE_LOCK) do
        _copy_mesh(_cached_mesh_locked("API.mesh.get"))
    end
end

function _mesh_query_integer(value,caller::AbstractString,name::AbstractString)
    value isa Integer || throw(ArgumentError("$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    return try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name exceeds the platform Int range"))
    end
end

function _mesh_query_dimension(value,caller::AbstractString)
    dimension=_mesh_query_integer(value,caller,"dim")
    dimension in -1:3 || throw(ArgumentError(
        "$caller: dim must be -1 or in 0:3"))
    return dimension
end

@inline function _mesh_query_bool(value,caller::AbstractString,
                                  name::AbstractString)
    value isa Bool || throw(ArgumentError("$caller: $name must be Bool"))
    return value
end

function _mesh_element_offsets(mesh::Mesh)
    return mesh_element_offsets(mesh)
end

function _mesh_dense_tags(offset::Int,count::Int)
    count==0 && return UInt64[]
    stop=Base.checked_add(offset,count)
    return UInt64.(offset+1:stop)
end

function _mesh_element_block(mesh::Mesh,element_type::Int)
    return mesh_element_block(mesh,element_type)
end

const _MESH_ENTITY_LABELS=("Point","Curve","Surface","Volume")

@inline function _mesh_entity_dictionary(m::GeoModel,dim::Int)
    return dim==0 ? m.points : dim==1 ? m.curves :
           dim==2 ? m.surfaces : m.volumes
end

# An entity exists when the model knows it or when the canonical
# `model_to_mixed` projection classified nodes/elements on it — the projection
# also enumerates boundary entities the `GeoModel` dictionaries never
# materialize (e.g. `add_box!` stores only the volume).
function _mesh_classified_entity(m::GeoModel,class::_MeshClassification,
                                 dim::Int,tag::Int,caller::AbstractString)
    typemin(Int32)<=tag<=typemax(Int32) || throw(ArgumentError(
        "$caller: unknown $(_MESH_ENTITY_LABELS[dim+1])[$tag]"))
    haskey(class.boundaries,(dim,Int32(tag))) ||
        haskey(_mesh_entity_dictionary(m,dim),tag) || throw(ArgumentError(
            "$caller: unknown $(_MESH_ENTITY_LABELS[dim+1])[$tag]"))
    return nothing
end

# Parses a Gmsh-style `dimTags` selection into concrete `(dim, tag)` pairs.
function _mesh_parse_dim_tags(dim_tags,caller::AbstractString)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    parsed=Tuple{Int,Int}[]
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        dimension,tag=pair
        dimension isa Integer || throw(ArgumentError(
            "$caller: entity dimensions must be integers"))
        dimension isa Bool && throw(ArgumentError(
            "$caller: entity dimensions must not be Bool"))
        0<=dimension<=3 || throw(ArgumentError(
            "$caller: entity dimension $dimension is out of range"))
        tag isa Integer || throw(ArgumentError(
            "$caller: entity tags must be integers"))
        tag isa Bool && throw(ArgumentError(
            "$caller: entity tags must not be Bool"))
        typemin(Int32)<=tag<=typemax(Int32) || throw(ArgumentError(
            "$caller: entity tag $tag is out of range"))
        push!(parsed,(Int(dimension),Int(tag)))
    end
    return parsed
end

# Validates a parsed selection against the classified entity set and returns it
# as a set of `(dim, tag)` pairs. Unknown entities fail explicitly, matching
# Gmsh's "entity does not exist" error.
function _mesh_selected_entities(m::GeoModel,class::_MeshClassification,
                                 pairs,caller::AbstractString)
    selected=Set{Tuple{Int,Int32}}()
    for (dim,tag) in pairs
        _mesh_classified_entity(m,class,dim,tag,caller)
        push!(selected,(dim,Int32(tag)))
    end
    return selected
end

# Entity-filtered cell positions, or `Int[]` for cells of `msh` types that
# cannot exist in the simplex cache.
function _mesh_entity_cell_positions(class::_MeshClassification,msh::Int,
                                     entity::Int)
    cell_entities=msh==1 ? class.seg_entities :
                  msh==2 ? class.tri_entities :
                  msh==4 ? class.tet_entities : Int32[]
    return findall(==(Int32(entity)),cell_entities)
end

# Returns `(msh_type, (offset, cells) or nothing, positions)`. `positions` is
# `nothing` for the complete block and otherwise the 1-based column indices of
# the cells classified on the requested entity.
function _mesh_query_type_block(mesh::Mesh,element_type,tag,
                                caller::AbstractString)
    msh=_mesh_query_integer(element_type,caller,"element_type")
    spec=msh_spec(msh)
    entity=_mesh_query_integer(tag,caller,"tag")
    entity<0 && return msh,_mesh_element_block(mesh,msh),nothing
    class=_cached_classification_locked(mesh)
    class===nothing && throw(ArgumentError(
        "$caller: entity-specific data require mesh classification metadata; " *
        "use a negative tag to query the complete cache"))
    _mesh_classified_entity(_model_locked(),class,spec.dim,entity,caller)
    block=_mesh_element_block(mesh,msh)
    positions=block===nothing ? Int[] :
        _mesh_entity_cell_positions(class,msh,entity)
    return msh,block,positions
end

@inline function _mesh_selected_columns(cells::AbstractMatrix{Int32},
                                        positions)
    return positions===nothing ? axes(cells,2) : positions
end

function _mesh_query_tasks(task,num_tasks,caller::AbstractString)
    task_index=_mesh_query_integer(task,caller,"task")
    task_count=_mesh_query_integer(num_tasks,caller,"num_tasks")
    task_index>=0 || throw(ArgumentError(
        "$caller: task must be nonnegative"))
    task_count>=1 || throw(ArgumentError(
        "$caller: num_tasks must be positive"))
    return task_index,task_count
end

# Gmsh 4.15.2 contiguous-block partition
# (src/common/gmsh.cpp: `begin=(task*count)/numTasks`,
# `end=((task+1)*count)/numTasks` with truncating division). Positions are
# 0-based element slots; the returned 1-based `UnitRange` is already
# intersected with the populated prefix, so `task>=num_tasks` yields `1:0`
# exactly as Gmsh's empty `[begin,end)` range does. Products use `Int128` so
# adversarial task counts cannot wrap before the truncating division.
function _mesh_task_range(count::Int,task_index::Int,task_count::Int)
    (count<=0 || task_index>=task_count) && return 1:0
    first_element=
        Int(div(Int128(task_index)*Int128(count),Int128(task_count)))+1
    last_element=
        Int(div(Int128(task_index+1)*Int128(count),Int128(task_count)))
    last_element=min(last_element,count)
    first_element>last_element && return 1:0
    return first_element:last_element
end

# Select the caller-visible input-tag slice for quality queries, which Gmsh
# partitions over the requested tag vector rather than over cached positions.
# The container contract is checked here with the core message so invalid
# containers fail identically; tag contents are validated by
# `mesh_element_qualities` on the selected slice only, matching Gmsh, which
# leaves entries outside `[begin,end)` untouched.
function _mesh_task_tags(values,task_index::Int,task_count::Int,
                         caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: element_tags must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    return values[_mesh_task_range(length(values),task_index,task_count)]
end

function _mesh_element_family(value,caller::AbstractString)
    value isa AbstractString || throw(ArgumentError(
        "$caller: family_name must be a string"))
    name=String(value)
    occursin('\0',name) && throw(ArgumentError(
        "$caller: family_name must not contain NUL"))
    normalized=lowercase(name)
    normalized=="point" && return :pnt
    normalized=="line" && return :lin
    normalized=="triangle" && return :tri
    normalized=="quadrangle" && return :qua
    normalized=="tetrahedron" && return :tet
    normalized=="hexahedron" && return :hex
    normalized=="prism" && return :pri
    normalized=="pyramid" && return :pyr
    normalized=="trihedron" && return :trih
    throw(ArgumentError(
        "$caller: unknown family_name $(repr(name)); expected Point, Line, " *
        "Triangle, Quadrangle, Tetrahedron, Hexahedron, Prism, Pyramid, or " *
        "Trihedron"))
end

function _get_element_type(family_name,order,serendip=false)
    caller="API.mesh.get_element_type"
    family=_mesh_element_family(family_name,caller)
    polynomial_order=_mesh_query_integer(order,caller,"order")
    incomplete=_mesh_query_bool(serendip,caller,"serendip")
    return Int32(msh_type(
        family,polynomial_order;serendipity=incomplete))
end

function _get_element_properties(element_type)
    caller="API.mesh.get_element_properties"
    msh=_mesh_query_integer(element_type,caller,"element_type")
    properties=msh_properties(msh)
    return properties.name,Int32(properties.dim),Int32(properties.order),
           Int32(properties.num_nodes),properties.local_node_coordinates,
           Int32(properties.num_primary_nodes)
end

function _get_integration_points(element_type,integration_type)
    return mesh_integration_points(
        element_type,integration_type;
        caller="API.mesh.get_integration_points")
end

function _mesh_query_coordinate(value,caller::AbstractString,
                                name::AbstractString)
    value isa Real || throw(ArgumentError(
        "$caller: $name must be a real number"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    converted=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError || err isa MethodError) ||
            rethrow()
        throw(ArgumentError(
            "$caller: $name must be Float64-representable"))
    end
    isfinite(converted) || throw(ArgumentError(
        "$caller: $name must be finite"))
    return converted
end

function _mesh_query_point(x,y,z,caller::AbstractString)
    return (_mesh_query_coordinate(x,caller,"x"),
            _mesh_query_coordinate(y,caller,"y"),
            _mesh_query_coordinate(z,caller,"z"))
end

function _mesh_query_element_tag(mesh::Mesh,value,caller::AbstractString)
    tag=_mesh_query_integer(value,caller,"element_tag")
    _,_,total=_mesh_element_offsets(mesh)
    1<=tag<=total || throw(ArgumentError(
        "$caller: unknown element tag $value; expected a dense tag in 1:$total"))
    return tag
end

function _no_element_at_coordinates(caller::AbstractString,p)
    throw(ArgumentError(
        "$caller: no element found at ($(p[1]), $(p[2]), $(p[3]))"))
end

function _get_element_by_coordinates(x,y,z,dim=-1,strict=false)
    caller="API.mesh.get_element_by_coordinates"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        p=_mesh_query_point(x,y,z,caller)
        dimension=_mesh_query_dimension(dim,caller)
        strict_mode=_mesh_query_bool(strict,caller,"strict")
        locator=_cached_mesh_locator_locked(cached)
        tags=_locate_elements(locator,p,dimension,strict_mode,caller)
        isempty(tags) && _no_element_at_coordinates(caller,p)
        element_tag=first(tags)
        record=mesh_element_record(cached,element_tag)
        coordinates,_,_=
            _local_coordinates(cached,Int(element_tag),p,caller)
        return element_tag,record.element_type,record.node_tags,coordinates...
    end
end

function _get_elements_by_coordinates(x,y,z,dim=-1,strict=false)
    caller="API.mesh.get_elements_by_coordinates"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        p=_mesh_query_point(x,y,z,caller)
        dimension=_mesh_query_dimension(dim,caller)
        strict_mode=_mesh_query_bool(strict,caller,"strict")
        locator=_cached_mesh_locator_locked(cached)
        tags=_locate_elements(locator,p,dimension,strict_mode,caller)
        isempty(tags) && _no_element_at_coordinates(caller,p)
        return tags
    end
end

function _get_local_coordinates_in_element(element_tag,x,y,z)
    caller="API.mesh.get_local_coordinates_in_element"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_element_tag(cached,element_tag,caller)
        p=_mesh_query_point(x,y,z,caller)
        coordinates,_,_=_local_coordinates(cached,tag,p,caller)
        return _require_local_coordinates(coordinates,caller,tag)
    end
end

function _get_element_qualities(element_tags,quality_name="minSICN",
                                task=0,num_tasks=1)
    caller="API.mesh.get_element_qualities"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        selected=_mesh_task_tags(element_tags,task_index,task_count,caller)
        mesh_element_qualities(cached,selected,quality_name)
    end
end

function _get_jacobians(element_type,local_coord,tag=-1,task=0,num_tasks=1)
    caller="API.mesh.get_jacobians"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        if positions!==nothing
            selected=_mesh_task_range(
                length(positions),task_index,task_count)
            return mesh_jacobians(
                cached,msh,local_coord,positions[selected])
        end
        count=block===nothing ? 0 : size(block[2],2)
        mesh_jacobians(
            cached,msh,local_coord,
            _mesh_task_range(count,task_index,task_count))
    end
end

function _get_jacobian(element_tag,local_coord)
    caller="API.mesh.get_jacobian"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_element_tag(cached,element_tag,caller)
        mesh_jacobian(cached,tag,local_coord)
    end
end

function _get_basis_functions(element_type,local_coord,function_space_type,
                              wanted_orientations=Int32[])
    return mesh_basis_functions(
        element_type,local_coord,function_space_type,wanted_orientations;
        caller="API.mesh.get_basis_functions")
end

function _get_number_of_orientations(element_type,function_space_type)
    return mesh_number_of_orientations(
        element_type,function_space_type;
        caller="API.mesh.get_number_of_orientations")
end

function _get_basis_functions_orientation(element_type,function_space_type,
                                          tag=-1,task=0,num_tasks=1)
    caller="API.mesh.get_basis_functions_orientation"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        if positions!==nothing
            selected=_mesh_task_range(
                length(positions),task_index,task_count)
            return mesh_basis_orientations(
                cached,msh,function_space_type,
                positions[selected];caller=caller)
        end
        count=block===nothing ? 0 : size(block[2],2)
        mesh_basis_orientations(
            cached,msh,function_space_type,
            _mesh_task_range(count,task_index,task_count);caller=caller)
    end
end

function _get_basis_functions_orientation_for_element(
    element_tag,function_space_type)
    caller="API.mesh.get_basis_functions_orientation_for_element"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_element_tag(cached,element_tag,caller)
        mesh_basis_orientation(
            cached,tag,function_space_type;caller=caller)
    end
end

function _get_number_of_keys(element_type,function_space_type)
    return mesh_number_of_keys(
        element_type,function_space_type;
        caller="API.mesh.get_number_of_keys")
end

# Hierarchical key queries need the edge catalog whenever the basis owns edge
# functions and the face catalog whenever it owns face functions; both are
# populated incrementally exactly like Gmsh's addMEdge/addMFace calls inside
# getKeys.
function _hierarchical_key_catalog_needs(function_space_type,family::Symbol)
    space=MeshFunctionSpaces._function_space(
        function_space_type,"API.mesh key query")
    space.hierarchical || return false,false
    counts=space.key_dimension==0 ?
        MeshFunctionSpaces._h1_counts(family,space.key_order) :
        MeshFunctionSpaces._hcurl_counts(family,space.key_order)
    return counts.edge>0,counts.face>0
end

function _get_keys(element_type,function_space_type,tag=-1,
                   return_coord=true)
    caller="API.mesh.get_keys"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        needs_edges,needs_faces=_hierarchical_key_catalog_needs(
            function_space_type,msh_spec(msh).family)
        edge_replacement=LAST_MESH_EDGES[]
        face_replacement=LAST_MESH_FACES[]
        selected_cells=nothing
        if block!==nothing
            _,cells=block
            selected_cells=positions===nothing ? cells :
                cells[:,positions]
            if needs_edges
                edge_replacement=_mesh_edge_topology_for_cells(
                    cached,edge_replacement,selected_cells,msh)
            end
            if needs_faces
                face_replacement=_mesh_face_topology_for_cells(
                    cached,face_replacement,selected_cells,msh)
            end
        end
        result=positions===nothing ?
            mesh_keys(
                cached,msh,function_space_type,edge_replacement,
                face_replacement;
                return_coord=return_coord,caller=caller) :
            mesh_keys(
                cached,msh,function_space_type,positions,edge_replacement,
                face_replacement;
                return_coord=return_coord,caller=caller)
        block!==nothing || return result
        needs_edges && (LAST_MESH_EDGES[]=edge_replacement)
        needs_faces && (LAST_MESH_FACES[]=face_replacement)
        return result
    end
end

function _get_keys_for_element(element_tag,function_space_type,
                               return_coord=true)
    caller="API.mesh.get_keys_for_element"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_element_tag(cached,element_tag,caller)
        record=mesh_element_record(cached,tag)
        needs_edges,needs_faces=_hierarchical_key_catalog_needs(
            function_space_type,msh_spec(Int(record.element_type)).family)
        edge_replacement=LAST_MESH_EDGES[]
        face_replacement=LAST_MESH_FACES[]
        cells=reshape(Int32.(record.node_tags),length(record.node_tags),1)
        if needs_edges
            edge_replacement=_mesh_edge_topology_for_cells(
                cached,edge_replacement,cells,Int(record.element_type))
        end
        if needs_faces
            face_replacement=_mesh_face_topology_for_cells(
                cached,face_replacement,cells,Int(record.element_type))
        end
        result=mesh_keys_for_element(
            cached,tag,function_space_type,edge_replacement,face_replacement;
            return_coord=return_coord,caller=caller)
        needs_edges && (LAST_MESH_EDGES[]=edge_replacement)
        needs_faces && (LAST_MESH_FACES[]=face_replacement)
        return result
    end
end

function _get_keys_information(type_keys,entity_keys,element_type,
                               function_space_type)
    return mesh_keys_information(
        type_keys,entity_keys,element_type,function_space_type;
        caller="API.mesh.get_keys_information")
end

# Cell columns classified on any of the selected entities. `owners` carries the
# entity tag per column and `dim` is the owning entities' dimension.
function _mesh_selected_cell_columns(owners::Vector{Int32},dim::Int,
                                     selected::Set{Tuple{Int,Int32}})
    return findall(owner->(dim,owner) in selected,owners)
end

function _create_edges(dim_tags=())
    caller="API.mesh.create_edges"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        if isempty(pairs)
            LAST_MESH_EDGES[]=_mesh_edge_topology(cached,LAST_MESH_EDGES[])
            return nothing
        end
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-selective topology creation requires mesh " *
            "classification metadata; pass an empty collection for the " *
            "complete cache"))
        selected=_mesh_selected_entities(
            _model_locked(),class,pairs,caller)
        replacement=LAST_MESH_EDGES[]
        for (dim,msh,cells,owners) in ((1,1,cached.segs,class.seg_entities),
                                       (2,2,cached.tris,class.tri_entities),
                                       (3,4,cached.tets,class.tet_entities))
            columns=_mesh_selected_cell_columns(owners,dim,selected)
            isempty(columns) && continue
            replacement=_mesh_edge_topology_for_cells(
                cached,replacement,cells[:,columns],msh)
        end
        LAST_MESH_EDGES[]=replacement
        nothing
    end
end

function _create_faces(dim_tags=())
    caller="API.mesh.create_faces"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        if isempty(pairs)
            LAST_MESH_FACES[]=_mesh_face_topology(cached,LAST_MESH_FACES[])
            return nothing
        end
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-selective topology creation requires mesh " *
            "classification metadata; pass an empty collection for the " *
            "complete cache"))
        selected=_mesh_selected_entities(
            _model_locked(),class,pairs,caller)
        replacement=LAST_MESH_FACES[]
        for (dim,msh,cells,owners) in ((2,2,cached.tris,class.tri_entities),
                                       (3,4,cached.tets,class.tet_entities))
            columns=_mesh_selected_cell_columns(owners,dim,selected)
            isempty(columns) && continue
            replacement=_mesh_face_topology_for_cells(
                cached,replacement,cells[:,columns],msh)
        end
        LAST_MESH_FACES[]=replacement
        nothing
    end
end

function _add_edges(edge_tags,edge_nodes)
    caller="API.mesh.add_edges"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        replacement=_mesh_add_edges(
            LAST_MESH_EDGES[],cached,edge_tags,edge_nodes,caller)
        LAST_MESH_EDGES[]=replacement
        nothing
    end
end

function _add_faces(face_type,face_tags,face_nodes)
    caller="API.mesh.add_faces"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        replacement=_mesh_add_faces(
            LAST_MESH_FACES[],cached,face_type,face_tags,face_nodes,caller)
        LAST_MESH_FACES[]=replacement
        nothing
    end
end

function _get_edges(node_tags)
    caller="API.mesh.get_edges"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _mesh_edges(LAST_MESH_EDGES[],cached,node_tags,caller)
    end
end

function _get_faces(face_type,node_tags)
    caller="API.mesh.get_faces"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _mesh_faces(LAST_MESH_FACES[],cached,face_type,node_tags,caller)
    end
end

function _get_all_edges()
    return lock(STATE_LOCK) do
        _cached_mesh_locked("API.mesh.get_all_edges")
        _mesh_all_edges(LAST_MESH_EDGES[])
    end
end

function _get_all_faces(face_type)
    return lock(STATE_LOCK) do
        _cached_mesh_locked("API.mesh.get_all_faces")
        _mesh_all_faces(
            LAST_MESH_FACES[],face_type,"API.mesh.get_all_faces")
    end
end

function _mesh_nodes_for_cells(mesh::Mesh,cells::Matrix{Int32})
    count=length(cells)
    coordinate_count=Base.checked_mul(3,count)
    node_tags=Vector{UInt64}(undef,count)
    coordinates=Vector{Float64}(undef,coordinate_count)
    cursor=0
    @inbounds for cell in axes(cells,2),local_node in axes(cells,1)
        cursor+=1
        node_index=Int(cells[local_node,cell])
        node_tags[cursor]=UInt64(node_index)
        offset=3cursor-3
        coordinates[offset+1]=mesh.coords[1,node_index]
        coordinates[offset+2]=mesh.coords[2,node_index]
        coordinates[offset+3]=mesh.coords[3,node_index]
    end
    return node_tags,coordinates
end

function _mesh_barycenters(mesh::Mesh,cells::AbstractMatrix{Int32},fast::Bool,
                           caller::AbstractString)
    nodes_per_element=size(cells,1)
    count=size(cells,2)
    result=Vector{Float64}(undef,Base.checked_mul(3,count))
    weight=fast ? 1.0 : inv(Float64(nodes_per_element))
    @inbounds for cell in 1:count,axis in 1:3
        value=0.0
        for local_node in 1:nodes_per_element
            value=muladd(weight,
                         mesh.coords[axis,Int(cells[local_node,cell])],value)
        end
        isfinite(value) || throw(ArgumentError(
            "$caller: element $cell coordinate sum is not Float64-representable"))
        result[3cell-3+axis]=value
    end
    return result
end

function _mesh_pattern_nodes(cells::AbstractMatrix{Int32},patterns)
    isempty(patterns) && return UInt64[]
    pattern_width=length(first(patterns))
    per_element=Base.checked_mul(length(patterns),pattern_width)
    count=Base.checked_mul(size(cells,2),per_element)
    result=Vector{UInt64}(undef,count)
    cursor=0
    @inbounds for cell in axes(cells,2),pattern in patterns,
                  local_node in pattern
        cursor+=1
        result[cursor]=UInt64(cells[local_node,cell])
    end
    return result
end

function _mesh_element_data(mesh::Mesh,dimension::Int)
    triangle_offset,tetrahedron_offset,_=_mesh_element_offsets(mesh)
    blocks=((Int32(1),1,0,mesh.segs),
            (Int32(2),2,triangle_offset,mesh.tris),
            (Int32(4),3,tetrahedron_offset,mesh.tets))
    element_types=Int32[]
    element_tags=Vector{Vector{UInt64}}()
    node_tags=Vector{Vector{UInt64}}()
    for (element_type,block_dimension,offset,cells) in blocks
        (dimension<0 || dimension==block_dimension) || continue
        count=size(cells,2)
        count==0 && continue
        push!(element_types,element_type)
        push!(element_tags,_mesh_dense_tags(offset,count))
        push!(node_tags,UInt64.(vec(cells)))
    end
    return element_types,element_tags,node_tags
end

function _mesh_element_types(mesh::Mesh,dimension::Int)
    result=Int32[]
    (dimension<0 || dimension==1) && nsegs(mesh)>0 && push!(result,Int32(1))
    (dimension<0 || dimension==2) && ntris(mesh)>0 && push!(result,Int32(2))
    (dimension<0 || dimension==3) && ntets(mesh)>0 && push!(result,Int32(4))
    return result
end

# Node indices classified on `(dim, tag)`, in node-tag order.
function _mesh_entity_node_positions(class::_MeshClassification,
                                     key::Tuple{Int,Int32})
    return findall(==(key),class.node_entities)
end

# Ordered node indices for one entity: its own classified nodes first, then
# the transitive boundary entities' nodes breadth-first with a visited set,
# matching Gmsh's `includeBoundary` emission order.
function _mesh_entity_nodes_with_boundary(class::_MeshClassification,
                                          dim::Int,tag::Int32)
    order=Int[]
    seen=Set{Tuple{Int,Int32}}()
    queue=Tuple{Int,Int32}[(dim,tag)]
    while !isempty(queue)
        key=popfirst!(queue)
        key in seen && continue
        push!(seen,key)
        append!(order,_mesh_entity_node_positions(class,key))
        key[1]>0 || continue
        for boundary in get(class.boundaries,key,Int32[])
            push!(queue,(key[1]-1,boundary))
        end
    end
    return order
end

# Entities of `dimension` known to the classification: any entity that owns
# nodes or appears in the boundary map.
function _mesh_classified_entities(class::_MeshClassification,dimension::Int)
    tags=Set{Int32}()
    for (entity_dim,entity_tag) in class.node_entities
        entity_dim==dimension && push!(tags,entity_tag)
    end
    for (entity_dim,entity_tag) in keys(class.boundaries)
        entity_dim==dimension && push!(tags,entity_tag)
    end
    return sort!(collect(tags))
end

function _mesh_nodes_payload(mesh::Mesh,indices::Vector{Int})
    count=length(indices)
    node_tags=UInt64.(indices)
    coordinates=Vector{Float64}(undef,Base.checked_mul(3,count))
    @inbounds for position in 1:count,axis in 1:3
        coordinates[3position-3+axis]=mesh.coords[axis,indices[position]]
    end
    return node_tags,coordinates
end

# Flat parametric coordinates of the node positions in `indices` on the entity
# `(dim, tag)`, matching Gmsh's `returnParametricCoord` output: one `u` per node
# on a Line, `(u, v)` per node on a Plane, and nothing for Points, Volumes, or
# all-dimension queries since they own no parametrization. Every classified
# dim-1 entity is a straight Line and every classified dim-2 entity is a
# loop-bounded Plane, so the exact-rational parametrization covers them all.
function _mesh_entity_parameters(model::GeoModel,mesh::Mesh,dim::Int,tag::Int,
                                 indices::Vector{Int})
    (dim==1 || dim==2) || return Float64[]
    coordinates=Vector{Float64}(undef,Base.checked_mul(3,length(indices)))
    @inbounds for position in eachindex(indices),axis in 1:3
        coordinates[3position-3+axis]=mesh.coords[axis,indices[position]]
    end
    return model_parametrization(model,dim,tag,coordinates)
end

function _get_nodes(dim=-1,tag=-1,include_boundary=false,
                    return_parametric_coord=true)
    caller="API.mesh.get_nodes"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        boundary=_mesh_query_bool(include_boundary,caller,"include_boundary")
        parametric=_mesh_query_bool(
            return_parametric_coord,caller,"return_parametric_coord")
        if dimension<0
            return UInt64.(1:nnodes(cached)),vec(copy(cached.coords)),
                Float64[]
        end
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: dimension-specific nodes require mesh classification " *
            "metadata; use dim=-1 for all cached nodes"))
        model=_model_locked()
        indices=Int[]
        parameters=Float64[]
        if entity>=0
            _mesh_classified_entity(model,class,dimension,entity,caller)
            group=boundary ?
                _mesh_entity_nodes_with_boundary(
                    class,dimension,Int32(entity)) :
                _mesh_entity_node_positions(
                    class,(dimension,Int32(entity)))
            append!(indices,group)
            parametric && append!(parameters,_mesh_entity_parameters(
                model,cached,dimension,entity,group))
        else
            for entity_tag in _mesh_classified_entities(class,dimension)
                group=boundary ?
                    _mesh_entity_nodes_with_boundary(
                        class,dimension,entity_tag) :
                    _mesh_entity_node_positions(
                        class,(dimension,entity_tag))
                append!(indices,group)
                parametric && append!(parameters,_mesh_entity_parameters(
                    model,cached,dimension,Int(entity_tag),group))
            end
        end
        node_tags,coordinates=_mesh_nodes_payload(cached,indices)
        return node_tags,coordinates,parameters
    end
end

function _get_node(node_tag)
    caller="API.mesh.get_node"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_integer(node_tag,caller,"node_tag")
        (tag>=1 && tag<=nnodes(cached)) || throw(ArgumentError(
            "$caller: unknown node $tag"))
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: node ownership requires mesh classification " *
            "metadata; generate a mesh so the cache owns entity ownership"))
        dimension,entity=class.node_entities[tag]
        parameters=_mesh_entity_parameters(
            _model_locked(),cached,dimension,Int(entity),[tag])
        return cached.coords[:,tag],parameters,dimension,Int(entity)
    end
end

# Sorted unique node indices on every member entity of Physical(dim, tag):
# each entity's own classified nodes, its transitive boundary entities' nodes,
# and any transitively embedded entities' nodes — matching Gmsh's per-group
# emission, which is a deduplicated ascending tag set.
function _mesh_physical_group_nodes(class::_MeshClassification,model::GeoModel,
                                    dimension::Int,entities)
    indices=Int[]
    seen_nodes=Set{Int}()
    seen_entities=Set{Tuple{Int,Int32}}()
    for entity_tag in entities
        queue=Tuple{Int,Int32}[(dimension,Int32(entity_tag))]
        while !isempty(queue)
            key=popfirst!(queue)
            key in seen_entities && continue
            push!(seen_entities,key)
            for node in _mesh_entity_node_positions(class,key)
                node in seen_nodes && continue
                push!(seen_nodes,node)
                push!(indices,node)
            end
            if key[1]>0
                for boundary in get(class.boundaries,key,Int32[])
                    push!(queue,(key[1]-1,boundary))
                end
            end
            for embedded in get(model.embeds,(key[1],Int(key[2])),
                                NTuple{2,Int}[])
                push!(queue,(embedded[1],Int32(embedded[2])))
            end
        end
    end
    return sort!(indices)
end

function _get_nodes_for_physical_group(dim,tag)
    caller="API.mesh.get_nodes_for_physical_group"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        group=_mesh_query_integer(tag,caller,"tag")
        model=_model_locked()
        typemin(Int32)<=group<=typemax(Int32) || return UInt64[],Float64[]
        entities=get(model.physical,(dimension,group),Int[])
        isempty(entities) && return UInt64[],Float64[]
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: physical-group nodes require mesh classification " *
            "metadata; generate a mesh so the cache owns entity ownership"))
        indices=_mesh_physical_group_nodes(class,model,dimension,entities)
        return _mesh_nodes_payload(cached,indices)
    end
end

function _get_embedded(dim,tag)
    caller="API.mesh.get_embedded"
    return lock(STATE_LOCK) do
        model=_model_locked()
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        entity=_mesh_query_integer(tag,caller,"tag")
        typemin(Int32)<=entity<=typemax(Int32) || throw(ArgumentError(
            "$caller: unknown $(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
        # Like _mesh_classified_entity, entities materialized only through the
        # classification snapshot (implicit primitive boundaries) also resolve.
        known=haskey(_mesh_entity_dictionary(model,dimension),entity)
        if !known
            cached=LAST_MESH[]
            class=cached===nothing ? nothing :
                  _cached_classification_locked(cached)
            known=class!==nothing &&
                  haskey(class.boundaries,(dimension,Int32(entity)))
        end
        known || throw(ArgumentError(
            "$caller: unknown $(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
        return Tuple{Int32,Int32}[Tuple{Int32,Int32}(embedded)
            for embedded in get(model.embeds,(dimension,entity),
                                NTuple{2,Int}[])]
    end
end

# `getSizes`-style lenient pair parsing: well-formedness (pair shape, integer,
# non-Bool members) is still validated, but out-of-range dimensions and tags
# simply report `0.0`, matching Gmsh 4.15.2's silent zeros.
function _mesh_parse_dim_tags_lenient(dim_tags,caller::AbstractString)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    parsed=Tuple{Integer,Integer}[]
    for entry in dim_tags
        pair=entry isa Pair ? (first(entry),last(entry)) :
             (entry isa Tuple && length(entry)==2) ? entry :
             throw(ArgumentError(
                 "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        pair[1] isa Integer || throw(ArgumentError(
            "$caller: entity dimensions must be integers"))
        pair[1] isa Bool && throw(ArgumentError(
            "$caller: entity dimensions must not be Bool"))
        pair[2] isa Integer || throw(ArgumentError(
            "$caller: entity tags must be integers"))
        pair[2] isa Bool && throw(ArgumentError(
            "$caller: entity tags must not be Bool"))
        push!(parsed,pair)
    end
    return parsed
end

function _get_sizes(dim_tags)
    caller="API.mesh.get_sizes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags_lenient(dim_tags,caller)
        return Float64[pair[1]==0 ? get(model.point_size,pair[2],0.0) : 0.0
                       for pair in pairs]
    end
end

function _get_elements(dim=-1,tag=-1)
    caller="API.mesh.get_elements"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        (dimension<0 || entity<0) &&
            return _mesh_element_data(cached,dimension)
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-specific elements require mesh classification " *
            "metadata; use a negative tag to query a complete dimension"))
        _mesh_classified_entity(_model_locked(),class,dimension,entity,caller)
        triangle_offset,tetrahedron_offset,_=_mesh_element_offsets(cached)
        blocks=((Int32(1),1,0,cached.segs,class.seg_entities),
                (Int32(2),2,triangle_offset,cached.tris,class.tri_entities),
                (Int32(4),3,tetrahedron_offset,cached.tets,
                 class.tet_entities))
        element_types=Int32[]
        element_tags=Vector{Vector{UInt64}}()
        node_tags=Vector{Vector{UInt64}}()
        for (element_type,block_dimension,offset,cells,owners) in blocks
            block_dimension==dimension || continue
            positions=findall(==(Int32(entity)),owners)
            isempty(positions) && continue
            push!(element_types,element_type)
            push!(element_tags,UInt64.(offset .+ positions))
            push!(node_tags,UInt64.(vec(@view cells[:,positions])))
        end
        return element_types,element_tags,node_tags
    end
end

function _get_element_types(dim=-1,tag=-1)
    caller="API.mesh.get_element_types"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        (dimension<0 || entity<0) &&
            return _mesh_element_types(cached,dimension)
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-specific element types require mesh " *
            "classification metadata; use a negative tag to query a complete " *
            "dimension"))
        _mesh_classified_entity(_model_locked(),class,dimension,entity,caller)
        result=Int32[]
        dimension==1 && any(==(Int32(entity)),class.seg_entities) &&
            push!(result,Int32(1))
        dimension==2 && any(==(Int32(entity)),class.tri_entities) &&
            push!(result,Int32(2))
        dimension==3 && any(==(Int32(entity)),class.tet_entities) &&
            push!(result,Int32(4))
        return result
    end
end

function _get_element(element_tag)
    caller="API.mesh.get_element"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_integer(element_tag,caller,"element_tag")
        triangle_offset,tetrahedron_offset,total=_mesh_element_offsets(cached)
        (tag>=1 && tag<=total) || throw(ArgumentError(
            "$caller: unknown element $tag"))
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: element classification requires mesh classification " *
            "metadata; generate a mesh so the cache owns entity ownership"))
        if tag<=triangle_offset
            return Int32(1),UInt64.(vec(cached.segs[:,tag])),1,
                Int(class.seg_entities[tag])
        elseif tag<=tetrahedron_offset
            position=tag-triangle_offset
            return Int32(2),UInt64.(vec(cached.tris[:,position])),2,
                Int(class.tri_entities[position])
        end
        position=tag-tetrahedron_offset
        return Int32(4),UInt64.(vec(cached.tets[:,position])),3,
            Int(class.tet_entities[position])
    end
end

function _get_elements_by_type(element_type,tag=-1,task=0,num_tasks=1)
    caller="API.mesh.get_elements_by_type"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return UInt64[],UInt64[]
        offset,cells=block
        columns=_mesh_selected_columns(cells,positions)
        selected=_mesh_task_range(length(columns),task_index,task_count)
        isempty(selected) && return UInt64[],UInt64[]
        chosen=@view columns[selected]
        tags=UInt64.(offset .+ chosen)
        return tags,UInt64.(vec(@view cells[:,chosen]))
    end
end

function _get_nodes_by_element_type(element_type,tag=-1,
                                    return_parametric_coord=true)
    caller="API.mesh.get_nodes_by_element_type"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        parametric=_mesh_query_bool(
            return_parametric_coord,caller,"return_parametric_coord")
        block===nothing && return UInt64[],Float64[],Float64[]
        _,cells=block
        selected=cells[:,_mesh_selected_columns(cells,positions)]
        node_tags,coordinates=_mesh_nodes_for_cells(cached,selected)
        parameters=Float64[]
        if parametric
            class=_cached_classification_locked(cached)
            if class!==nothing
                # Gmsh packs each repeated node's parameters on its owning
                # entity — one `u` for Line owners, `(u, v)` for Plane owners,
                # and nothing for Point or Volume owners — in entry order.
                node_list=Int.(vec(selected))
                groups=Dict{Tuple{Int,Int32},Vector{Int}}()
                for node in unique(node_list)
                    owner=class.node_entities[node]
                    owner[1] in (1,2) || continue
                    push!(get!(groups,owner,Int[]),node)
                end
                node_parameters=Dict{Int,Vector{Float64}}()
                model=_model_locked()
                for ((owner_dim,owner_tag),members) in groups
                    values=_mesh_entity_parameters(
                        model,cached,owner_dim,Int(owner_tag),members)
                    for (index,node) in enumerate(members)
                        node_parameters[node]=values[
                            (index-1)*owner_dim+1:index*owner_dim]
                    end
                end
                for node in node_list
                    append!(parameters,get(
                        node_parameters,node,Float64[]))
                end
            end
        end
        return node_tags,coordinates,parameters
    end
end

function _get_barycenters(element_type,tag,fast,primary,task=0,num_tasks=1)
    caller="API.mesh.get_barycenters"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        fast_mode=_mesh_query_bool(fast,caller,"fast")
        _mesh_query_bool(primary,caller,"primary")
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return Float64[]
        _,cells=block
        columns=_mesh_selected_columns(cells,positions)
        selected=_mesh_task_range(length(columns),task_index,task_count)
        _mesh_barycenters(
            cached,@view(cells[:,columns[selected]]),fast_mode,caller)
    end
end

function _get_element_edge_nodes(element_type,tag=-1,primary=false,
                                 task=0,num_tasks=1)
    caller="API.mesh.get_element_edge_nodes"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        _mesh_query_bool(primary,caller,"primary")
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return UInt64[]
        _,cells=block
        columns=_mesh_selected_columns(cells,positions)
        selected=_mesh_task_range(length(columns),task_index,task_count)
        _mesh_pattern_nodes(
            @view(cells[:,columns[selected]]),_simplex_edge_patterns(msh))
    end
end

function _get_element_face_nodes(element_type,face_type,tag=-1,primary=false,
                                 task=0,num_tasks=1)
    caller="API.mesh.get_element_face_nodes"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block,positions=_mesh_query_type_block(
            cached,element_type,tag,caller)
        face=_mesh_query_integer(face_type,caller,"face_type")
        face in (3,4) || throw(ArgumentError(
            "$caller: face_type must be 3 (triangle) or 4 (quadrangle)"))
        _mesh_query_bool(primary,caller,"primary")
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return UInt64[]
        _,cells=block
        columns=_mesh_selected_columns(cells,positions)
        selected=_mesh_task_range(length(columns),task_index,task_count)
        _mesh_pattern_nodes(
            @view(cells[:,columns[selected]]),_simplex_face_patterns(msh,face))
    end
end

function _get_max_node_tag()
    return lock(STATE_LOCK) do
        UInt64(nnodes(_cached_mesh_locked("API.mesh.get_max_node_tag")))
    end
end

function _get_max_element_tag()
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked("API.mesh.get_max_element_tag")
        _,_,total=_mesh_element_offsets(cached)
        UInt64(total)
    end
end

function _refine(;max_nodes=typemax(Int32),max_cells=typemax(Int32))
    return lock(STATE_LOCK) do
        m=_model_locked()
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "API.mesh.refine: no mesh; call API.mesh.generate first"))
        refined=refine_uniform(
            cached;max_nodes=max_nodes,max_cells=max_cells)
        class=_cached_classification_locked(cached)
        cache=_copy_mesh(refined)
        new_class=class===nothing ? nothing : _classify_cached_mesh(
            m,refined,class.entity[1],Int(class.entity[2]),cache)
        _replace_mesh_cache_locked!(cache,new_class)
        refined
    end
end

# Selective clear on the flat shared-node cache, mirroring Gmsh's per-entity
# mesh ownership: cells classified on a cleared entity are removed; nodes owned
# by such an entity drop when no surviving cell references them and reclassify
# to the lowest-dimension (then lowest-tag) surviving owner otherwise. Entities
# that own no cells in the cache — points, or boundary entities of a higher-
# dimensional generating entity — are no-ops, matching Gmsh 4.15.2's retention
# of boundary meshes and unmeshable vertices.
function _clear_classified_mesh(mesh::Mesh,class::_MeshClassification,
                                cleared::Set{Tuple{Int,Int32}})
    cell_blocks=((1,mesh.segs,class.seg_entities),
                 (2,mesh.tris,class.tri_entities),
                 (3,mesh.tets,class.tet_entities))
    cleared_with_cells=Set{Tuple{Int,Int32}}()
    kept_columns=[Vector{Int}() for _ in cell_blocks]
    for (block,(dim,cells,owners)) in enumerate(cell_blocks)
        kept=kept_columns[block]
        sizehint!(kept,size(cells,2))
        for column in axes(cells,2)
            if (dim,owners[column]) in cleared
                push!(cleared_with_cells,(dim,owners[column]))
            else
                push!(kept,column)
            end
        end
    end
    isempty(cleared_with_cells) && return mesh,class
    count=nnodes(mesh)
    # Lowest-dimension (then lowest-tag) surviving owner referencing each node.
    survivors=Vector{Tuple{Int,Int32}}(undef,count)
    fill!(survivors,(4,Int32(0)))
    for (block,(dim,cells,owners)) in enumerate(cell_blocks)
        for column in kept_columns[block]
            owner=(dim,owners[column])
            for row in axes(cells,1)
                node=Int(cells[row,column])
                owner<survivors[node] && (survivors[node]=owner)
            end
        end
    end
    referenced=falses(count)
    for (block,(dim,cells,owners)) in enumerate(cell_blocks)
        for column in kept_columns[block], row in axes(cells,1)
            referenced[Int(cells[row,column])]=true
        end
    end
    keep=falses(count)
    for node in 1:count
        owner=class.node_entities[node]
        if owner in cleared_with_cells
            keep[node]=referenced[node]
        else
            keep[node]=true
        end
    end
    old_to_new=Vector{Int32}(undef,count)
    index=0
    for node in 1:count
        keep[node] && (index+=1;old_to_new[node]=index)
    end
    index==0 && return nothing,nothing
    coordinates=mesh.coords[:,keep]
    node_entities=Vector{Tuple{Int,Int32}}(undef,index)
    out=0
    for node in 1:count
        keep[node] || continue
        out+=1
        owner=class.node_entities[node]
        node_entities[out]=
            owner in cleared_with_cells ? survivors[node] : owner
    end
    remapped=Vector{Matrix{Int32}}(undef,3)
    cell_entities=Vector{Vector{Int32}}(undef,3)
    for (block,(dim,cells,owners)) in enumerate(cell_blocks)
        kept=kept_columns[block]
        result=Matrix{Int32}(undef,size(cells,1),length(kept))
        for (target,column) in enumerate(kept), row in axes(cells,1)
            result[row,target]=old_to_new[cells[row,column]]
        end
        remapped[block]=result
        cell_entities[block]=owners[kept]
    end
    replacement=Mesh(coordinates;segs=remapped[1],tris=remapped[2],
                     tets=remapped[3],
                     seg_tag=mesh.seg_tag[kept_columns[1]],
                     tri_tag=mesh.tri_tag[kept_columns[2]],
                     tet_tag=mesh.tet_tag[kept_columns[3]])
    record=_MeshClassification(replacement,class.entity,node_entities,
                               class.boundaries,cell_entities[1],
                               cell_entities[2],cell_entities[3])
    return replacement,record
end

function _clear_mesh(dim_tags=())
    caller="API.mesh.clear"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        isempty(pairs) && (_replace_mesh_cache_locked!(nothing);return nothing)
        cached=LAST_MESH[]
        cached===nothing && return nothing
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-selective clearing requires mesh classification " *
            "metadata; pass an empty collection to clear the complete cache"))
        cleared=_mesh_selected_entities(model,class,pairs,caller)
        replacement,record=_clear_classified_mesh(cached,class,cleared)
        replacement===cached ||
            _replace_mesh_cache_locked!(replacement,record)
        nothing
    end
end

function _affine_transform_mesh(affine,dim_tags=())
    caller="API.mesh.affine_transform"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "$caller: no mesh; call API.mesh.generate first"))
        coefficients,translation,_=_transform_gmsh_affine(affine,caller)
        matrix=reshape(collect(coefficients),3,3)
        if isempty(pairs)
            transformed=affine_transform(
                cached,matrix;translation=translation)
        else
            class=_cached_classification_locked(cached)
            class===nothing && throw(ArgumentError(
                "$caller: entity-selective transformation requires mesh " *
                "classification metadata; pass an empty collection to " *
                "transform the complete cached mesh"))
            selected=_mesh_selected_entities(model,class,pairs,caller)
            mask=falses(nnodes(cached))
            for node in eachindex(class.node_entities)
                class.node_entities[node] in selected && (mask[node]=true)
            end
            transformed=affine_transform(
                cached,matrix;translation=translation,node_mask=mask)
        end
        # An affine map preserves connectivity, so the classification stays
        # index-aligned; it is rebound to the transformed cache object.
        class=_cached_classification_locked(cached)
        cache=_copy_mesh(transformed)
        new_class=class===nothing ? nothing : _MeshClassification(
            cache,class.entity,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(cache,new_class)
        transformed
    end
end

# Rebuild the cache after dropping whole element columns from one block. Nodes
# and their entity ownership are retained (Gmsh 4.15.2 keeps orphan nodes);
# only the selected block's connectivity, per-cell owners, and tags shrink.
function _remove_element_columns(mesh::Mesh,class::_MeshClassification,
                                 dimension::Int,dropped::AbstractVector{Bool})
    keep=.!dropped
    replacement=Mesh(mesh.coords;
        segs=dimension==1 ? mesh.segs[:,keep] : mesh.segs,
        tris=dimension==2 ? mesh.tris[:,keep] : mesh.tris,
        tets=dimension==3 ? mesh.tets[:,keep] : mesh.tets,
        seg_tag=dimension==1 ? mesh.seg_tag[keep] : mesh.seg_tag,
        tri_tag=dimension==2 ? mesh.tri_tag[keep] : mesh.tri_tag,
        tet_tag=dimension==3 ? mesh.tet_tag[keep] : mesh.tet_tag)
    record=_MeshClassification(replacement,class.entity,class.node_entities,
        class.boundaries,
        dimension==1 ? class.seg_entities[keep] : class.seg_entities,
        dimension==2 ? class.tri_entities[keep] : class.tri_entities,
        dimension==3 ? class.tet_entities[keep] : class.tet_entities)
    _replace_mesh_cache_locked!(replacement,record)
    return nothing
end

function _remove_elements(dim,tag,element_tags=())
    caller="API.mesh.remove_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        entity=_mesh_query_integer(tag,caller,"tag")
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: removing elements requires mesh classification " *
            "metadata; generate a mesh so the cache owns entity ownership"))
        _mesh_classified_entity(model,class,dimension,entity,caller)
        triangle_offset,tetrahedron_offset,total=_mesh_element_offsets(cached)
        # The simplex cache owns no dimension-0 cells; a Point selection can
        # only no-op on an empty list or reject tags classified elsewhere.
        cells=dimension==0 ? Matrix{Int32}(undef,2,0) :
              dimension==1 ? cached.segs :
              dimension==2 ? cached.tris : cached.tets
        owners=dimension==0 ? Int32[] :
               dimension==1 ? class.seg_entities :
               dimension==2 ? class.tri_entities : class.tet_entities
        (element_tags isa AbstractVector || element_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: element_tags must be a vector or tuple of " *
                "dense element tags"))
        dropped=falses(size(cells,2))
        if isempty(element_tags)
            for column in axes(cells,2)
                owners[column]==Int32(entity) && (dropped[column]=true)
            end
        else
            for value in element_tags
                dense=_mesh_query_integer(value,caller,"element_tags entry")
                (1<=dense<=total) || throw(ArgumentError(
                    "$caller: unknown element $dense; expected a dense tag " *
                    "in 1:$total"))
                block_dimension=dense<=triangle_offset ? 1 :
                                dense<=tetrahedron_offset ? 2 : 3
                position=dense<=triangle_offset ? dense :
                         dense<=tetrahedron_offset ? dense-triangle_offset :
                         dense-tetrahedron_offset
                block_owners=block_dimension==1 ? class.seg_entities :
                             block_dimension==2 ? class.tri_entities :
                             class.tet_entities
                (block_dimension==dimension &&
                 block_owners[position]==Int32(entity)) || throw(ArgumentError(
                    "$caller: element $dense is not classified on " *
                    "$(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
                dropped[position]=true
            end
        end
        any(dropped) &&
            _remove_element_columns(cached,class,dimension,dropped)
        return nothing
    end
end

# Gmsh 4.15.2 `reverse` node-order conventions for first-order simplices:
# segments swap both vertices, triangles swap positions 2 and 3, tetrahedra
# swap positions 1 and 2.
function _reverse_simplex_columns!(cells::AbstractMatrix{Int32},positions)
    for position in positions
        if size(cells,1)==2
            cells[1,position],cells[2,position]=
                cells[2,position],cells[1,position]
        elseif size(cells,1)==3
            cells[2,position],cells[3,position]=
                cells[3,position],cells[2,position]
        else
            cells[1,position],cells[2,position]=
                cells[2,position],cells[1,position]
        end
    end
    return nothing
end

function _reverse_mesh(dim_tags=())
    caller="API.mesh.reverse"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=_cached_mesh_locked(caller)
        class=_cached_classification_locked(cached)
        isempty(pairs) || class===nothing && throw(ArgumentError(
            "$caller: entity-selective reversal requires mesh classification " *
            "metadata; pass an empty collection to reverse the complete cache"))
        selected=isempty(pairs) ? nothing :
                 _mesh_selected_entities(model,class,pairs,caller)
        replacement=_copy_mesh(cached)
        for (dimension,cells,owners) in (
                (1,replacement.segs,class===nothing ? Int32[] :
                    class.seg_entities),
                (2,replacement.tris,class===nothing ? Int32[] :
                    class.tri_entities),
                (3,replacement.tets,class===nothing ? Int32[] :
                    class.tet_entities))
            positions=Int[]
            for column in axes(cells,2)
                (selected===nothing ||
                 (dimension,owners[column]) in selected) &&
                    push!(positions,column)
            end
            _reverse_simplex_columns!(cells,positions)
        end
        # Connectivity is unchanged as a set, so the classification stays
        # index-aligned and is rebound to the reversed cache object.
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        return nothing
    end
end

function _reverse_elements(element_tags)
    caller="API.mesh.reverse_elements"
    return lock(STATE_LOCK) do
        _model_locked()
        cached=_cached_mesh_locked(caller)
        (element_tags isa AbstractVector || element_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: element_tags must be a vector or tuple of dense " *
                "element tags"))
        triangle_offset,tetrahedron_offset,total=_mesh_element_offsets(cached)
        per_block=(Int[],Int[],Int[])
        for value in element_tags
            dense=_mesh_query_integer(value,caller,"element_tags entry")
            (1<=dense<=total) || throw(ArgumentError(
                "$caller: unknown element $dense; expected a dense tag in " *
                "1:$total"))
            push!(per_block[dense<=triangle_offset ? 1 :
                            dense<=tetrahedron_offset ? 2 : 3],
                  dense<=triangle_offset ? dense :
                  dense<=tetrahedron_offset ? dense-triangle_offset :
                  dense-tetrahedron_offset)
        end
        replacement=_copy_mesh(cached)
        _reverse_simplex_columns!(replacement.segs,per_block[1])
        _reverse_simplex_columns!(replacement.tris,per_block[2])
        _reverse_simplex_columns!(replacement.tets,per_block[3])
        class=_cached_classification_locked(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        return nothing
    end
end

# Merge nodes sharing exact coordinates, keeping the lowest tag, then rebuild
# the cache with compacted node numbering and remapped connectivity — matching
# Gmsh 4.15.2's `removeDuplicateNodes`.
function _remove_duplicate_nodes(dim_tags=())
    caller="API.mesh.remove_duplicate_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=_cached_mesh_locked(caller)
        class=_cached_classification_locked(cached)
        isempty(pairs) || class===nothing && throw(ArgumentError(
            "$caller: entity-selective duplicate removal requires mesh " *
            "classification metadata; pass an empty collection to scan the " *
            "complete cache"))
        selected=isempty(pairs) ? nothing :
                 _mesh_selected_entities(model,class,pairs,caller)
        count=nnodes(cached)
        replacement=collect(Int32,1:count)
        groups=Dict{NTuple{3,Float64},Int32}()
        for node in 1:count
            (selected===nothing ||
             class.node_entities[node] in selected) || continue
            key=NTuple{3,Float64}((cached.coords[1,node],
                                   cached.coords[2,node],
                                   cached.coords[3,node]))
            replacement[node]=Int32(get!(groups,key,node))
        end
        keep=Bool[replacement[node]==node for node in 1:count]
        all(keep) && return nothing
        old_to_new=Vector{Int32}(undef,count)
        index=0
        for node in 1:count
            keep[node] && (index+=1;old_to_new[node]=index)
        end
        coordinates=cached.coords[:,keep]
        blocks=Vector{Matrix{Int32}}(undef,3)
        for (block,cells) in enumerate((
                cached.segs,cached.tris,cached.tets))
            result=Matrix{Int32}(undef,size(cells,1),size(cells,2))
            for column in axes(cells,2),row in axes(cells,1)
                result[row,column]=
                    old_to_new[Int(replacement[cells[row,column]])]
            end
            blocks[block]=result
        end
        node_entities=class===nothing ? nothing :
            Tuple{Int,Int32}[class.node_entities[node]
                             for node in 1:count if keep[node]]
        new_mesh=Mesh(coordinates;segs=blocks[1],tris=blocks[2],
                      tets=blocks[3],seg_tag=cached.seg_tag,
                      tri_tag=cached.tri_tag,tet_tag=cached.tet_tag)
        new_class=class===nothing ? nothing : _MeshClassification(
            new_mesh,class.entity,node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(new_mesh,new_class)
        return nothing
    end
end

# Drop cells whose sorted connectivity repeats an earlier cell owned by the
# same entity — Gmsh 4.15.2's `removeDuplicateElements` contract.
function _remove_duplicate_elements(dim_tags=())
    caller="API.mesh.remove_duplicate_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=_cached_mesh_locked(caller)
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: duplicate-element removal requires mesh classification " *
            "metadata; generate a mesh so the cache owns entity ownership"))
        selected=isempty(pairs) ? nothing :
                 _mesh_selected_entities(model,class,pairs,caller)
        changed=false
        new_blocks=Vector{Matrix{Int32}}(undef,3)
        new_owners=Vector{Vector{Int32}}(undef,3)
        new_tags=Vector{Vector{Int32}}(undef,3)
        cell_tags=(cached.seg_tag,cached.tri_tag,cached.tet_tag)
        for (block,(dimension,cells,owners)) in enumerate((
                (1,cached.segs,class.seg_entities),
                (2,cached.tris,class.tri_entities),
                (3,cached.tets,class.tet_entities)))
            keep=trues(size(cells,2))
            seen=Dict{Tuple{Int32,NTuple{4,Int32}},Int}()
            for column in axes(cells,2)
                (selected===nothing ||
                 (dimension,owners[column]) in selected) || continue
                connectivity=sort!(vec(cells[:,column]))
                key=(owners[column],ntuple(
                    slot->slot<=length(connectivity) ? connectivity[slot] :
                        Int32(0),4))
                haskey(seen,key) && (keep[column]=false;changed=true)
                seen[key]=column
            end
            new_blocks[block]=cells[:,keep]
            new_owners[block]=owners[keep]
            new_tags[block]=cell_tags[block][keep]
        end
        changed || return nothing
        replacement=Mesh(cached.coords;segs=new_blocks[1],
                         tris=new_blocks[2],tets=new_blocks[3],
                         seg_tag=new_tags[1],tri_tag=new_tags[2],
                         tet_tag=new_tags[3])
        record=_MeshClassification(replacement,class.entity,
            class.node_entities,class.boundaries,new_owners[1],
            new_owners[2],new_owners[3])
        _replace_mesh_cache_locked!(replacement,record)
        return nothing
    end
end

function _get_duplicate_nodes(dim_tags=())
    caller="API.mesh.get_duplicate_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=_cached_mesh_locked(caller)
        class=_cached_classification_locked(cached)
        selected=nothing
        if !isempty(pairs)
            class===nothing && throw(ArgumentError(
                "$caller: entity-selective duplicate detection requires mesh " *
                "classification metadata; pass an empty collection to scan " *
                "the complete cache"))
            selected=_mesh_selected_entities(model,class,pairs,caller)
        end
        groups=Dict{NTuple{3,Float64},Vector{Int}}()
        for node in axes(cached.coords,2)
            (selected===nothing ||
             class.node_entities[node] in selected) || continue
            push!(get!(groups,Tuple(cached.coords[:,node]),Int[]),node)
        end
        duplicates=UInt64[]
        for members in values(groups)
            length(members)>1 && append!(duplicates,UInt64.(members))
        end
        return sort!(duplicates)
    end
end

# Set one node's Cartesian coordinates — Gmsh 4.15.2 `setNode`. Parametric
# coordinates are not stored: `get_node`/`get_nodes` recompute them from the
# owning entity's geometry, so a nonempty `parametric_coord` is rejected rather
# than silently ignored.
function _set_node(node_tag,coord,parametric_coord=Float64[])
    caller="API.mesh.set_node"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        tag=_mesh_query_integer(node_tag,caller,"node_tag")
        count=nnodes(cached)
        (1<=tag<=count) || throw(ArgumentError(
            "$caller: unknown node $tag; expected a dense tag in 1:$count"))
        (coord isa AbstractVector || coord isa Tuple) ||
            throw(ArgumentError(
                "$caller: coord must be a vector or tuple of 3 coordinates"))
        length(coord)==3 || throw(ArgumentError(
            "$caller: coord must hold exactly 3 coordinates"))
        values=Float64[Float64(c) for c in coord]
        all(isfinite,values) || throw(ArgumentError(
            "$caller: node coordinates must be finite"))
        (parametric_coord isa AbstractVector ||
         parametric_coord isa Tuple) || throw(ArgumentError(
            "$caller: parametric_coord must be a vector or tuple"))
        isempty(parametric_coord) || throw(ArgumentError(
            "$caller: parametric coordinates are recomputed from the owning " *
            "entity; pass an empty parametric_coord"))
        coordinates=copy(cached.coords)
        coordinates[:,tag]=values
        replacement=Mesh(coordinates;segs=cached.segs,tris=cached.tris,
                         tets=cached.tets,seg_tag=cached.seg_tag,
                         tri_tag=cached.tri_tag,tet_tag=cached.tet_tag)
        class=_cached_classification_locked(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        return nothing
    end
end

# Parses an explicit `old_tags`/`new_tags` renumbering pair list into a dense
# permutation of `1:count`, or `nothing` when both lists are empty (Gmsh's
# "renumber continuously" form, already a no-op on the dense cache).
function _mesh_renumber_permutation(old_tags,new_tags,count,caller,label)
    (old_tags isa AbstractVector || old_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: old_tags must be a vector or tuple of $label tags"))
    (new_tags isa AbstractVector || new_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: new_tags must be a vector or tuple of $label tags"))
    length(old_tags)==length(new_tags) || throw(ArgumentError(
        "$caller: old_tags and new_tags must have the same length"))
    isempty(old_tags) && return nothing
    mapping=collect(Int32,1:count)
    seen_new=Set{Int}()
    for (old_value,new_value) in zip(old_tags,new_tags)
        old=_mesh_query_integer(old_value,caller,"old_tags entry")
        new=_mesh_query_integer(new_value,caller,"new_tags entry")
        (1<=old<=count) || throw(ArgumentError(
            "$caller: unknown $label tag $old; expected a dense tag in " *
            "1:$count"))
        new in seen_new && throw(ArgumentError(
            "$caller: new_tags contains tag $new twice"))
        push!(seen_new,new)
        mapping[old]=Int32(new)
    end
    # The flat cache only represents the dense tag space 1:count, so the
    # resulting assignment must stay a permutation of it.
    sort!(collect(Int,mapping))==collect(1:count) || throw(ArgumentError(
        "$caller: renumbering must keep $label tags a permutation of " *
        "1:$count; sparse tags are not representable on the dense cache"))
    return mapping
end

function _renumber_nodes(old_tags=(),new_tags=())
    caller="API.mesh.renumber_nodes"
    return lock(STATE_LOCK) do
        _model_locked()
        cached=_cached_mesh_locked(caller)
        mapping=_mesh_renumber_permutation(
            old_tags,new_tags,nnodes(cached),caller,"node")
        mapping===nothing && return nothing
        order=Vector{Int}(undef,nnodes(cached))
        for node in 1:nnodes(cached)
            order[mapping[node]]=node
        end
        coordinates=cached.coords[:,order]
        inverse=Vector{Int32}(undef,nnodes(cached))
        for node in 1:nnodes(cached)
            inverse[node]=Int32(mapping[node])
        end
        blocks=Vector{Matrix{Int32}}(undef,3)
        for (block,cells) in enumerate((
                cached.segs,cached.tris,cached.tets))
            result=Matrix{Int32}(undef,size(cells,1),size(cells,2))
            for column in axes(cells,2),row in axes(cells,1)
                result[row,column]=inverse[cells[row,column]]
            end
            blocks[block]=result
        end
        class=_cached_classification_locked(cached)
        node_entities=class===nothing ? nothing :
            Tuple{Int,Int32}[class.node_entities[node] for node in order]
        new_mesh=Mesh(coordinates;segs=blocks[1],tris=blocks[2],
                      tets=blocks[3],seg_tag=cached.seg_tag,
                      tri_tag=cached.tri_tag,tet_tag=cached.tet_tag)
        new_class=class===nothing ? nothing : _MeshClassification(
            new_mesh,class.entity,node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(new_mesh,new_class)
        return nothing
    end
end

# Element tags are dense across the segment/triangle/tetrahedron blocks, so a
# renumbering must keep every tag inside its own block range — cross-block
# assignments would change element types and are rejected explicitly.
function _renumber_elements(old_tags=(),new_tags=())
    caller="API.mesh.renumber_elements"
    return lock(STATE_LOCK) do
        _model_locked()
        cached=_cached_mesh_locked(caller)
        triangle_offset,tetrahedron_offset,total=_mesh_element_offsets(cached)
        mapping=_mesh_renumber_permutation(
            old_tags,new_tags,total,caller,"element")
        mapping===nothing && return nothing
        for position in 1:total
            target=Int(mapping[position])
            block_of=position<=triangle_offset ? 1 :
                     position<=tetrahedron_offset ? 2 : 3
            target_block=target<=triangle_offset ? 1 :
                         target<=tetrahedron_offset ? 2 : 3
            block_of==target_block || throw(ArgumentError(
                "$caller: element $position cannot take tag $target; " *
                "renumbering cannot move elements across element types"))
        end
        class=_cached_classification_locked(cached)
        blocks=Vector{Matrix{Int32}}(undef,3)
        owners_out=Vector{Vector{Int32}}(undef,3)
        tags_out=Vector{Vector{Int32}}(undef,3)
        cell_tags=(cached.seg_tag,cached.tri_tag,cached.tet_tag)
        owners_in=(class===nothing ? Int32[] : class.seg_entities,
                   class===nothing ? Int32[] : class.tri_entities,
                   class===nothing ? Int32[] : class.tet_entities)
        for (block,(cells,offset)) in enumerate(zip(
                (cached.segs,cached.tris,cached.tets),
                (0,triangle_offset,tetrahedron_offset)))
            width=size(cells,2)
            order=Vector{Int}(undef,width)
            for column in 1:width
                order[Int(mapping[offset+column])-offset]=column
            end
            blocks[block]=cells[:,order]
            owners_out[block]=isempty(owners_in[block]) ? Int32[] :
                              owners_in[block][order]
            tags_out[block]=cell_tags[block][order]
        end
        replacement=Mesh(cached.coords;segs=blocks[1],tris=blocks[2],
                         tets=blocks[3],seg_tag=tags_out[1],
                         tri_tag=tags_out[2],tet_tag=tags_out[3])
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.node_entities,class.boundaries,
            owners_out[1],owners_out[2],owners_out[3])
        _replace_mesh_cache_locked!(replacement,new_class)
        return nothing
    end
end

# Reorder one entity's cells of one element type — Gmsh 4.15.2's
# `reorderElements`, whose `ordering` holds 0-based source positions:
# `ordering[new_position]` is the column that moves there. Elements keep their
# stored tags and entity ownership.
function _reorder_elements(element_type,tag,ordering)
    caller="API.mesh.reorder_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=_cached_mesh_locked(caller)
        msh=_mesh_query_integer(element_type,caller,"element_type")
        spec=msh_spec(msh)
        dimension=spec.dim
        entity=_mesh_query_integer(tag,caller,"tag")
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: reordering requires mesh classification metadata"))
        _mesh_classified_entity(model,class,dimension,entity,caller)
        cells=dimension==1 ? cached.segs :
              dimension==2 ? cached.tris :
              dimension==3 ? cached.tets :
              Matrix{Int32}(undef,spec.nnodes,0)
        owners=dimension==1 ? class.seg_entities :
               dimension==2 ? class.tri_entities :
               dimension==3 ? class.tet_entities : Int32[]
        positions=findall(==(Int32(entity)),owners)
        (ordering isa AbstractVector || ordering isa Tuple) ||
            throw(ArgumentError(
                "$caller: ordering must be a vector or tuple of 0-based " *
                "source positions"))
        isempty(positions) && throw(ArgumentError(
            "$caller: no elements of type $msh classified on " *
            "$(_MESH_ENTITY_LABELS[dimension+1])[$entity] to reorder"))
        length(ordering)==length(positions) || throw(ArgumentError(
            "$caller: ordering must hold $(length(positions)) entries for the " *
            "elements of type $msh on " *
            "$(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
        permutation=Int[]
        for value in ordering
            source=_mesh_query_integer(value,caller,"ordering entry")
            (0<=source<length(positions)) || throw(ArgumentError(
                "$caller: ordering entry $source is out of range; expected " *
                "0:$(length(positions)-1)"))
            source in permutation && throw(ArgumentError(
                "$caller: ordering repeats position $source"))
            push!(permutation,source)
        end
        reordered=positions[permutation .+ 1]
        blocks=(copy(cached.segs),copy(cached.tris),copy(cached.tets))
        old_tags=(cached.seg_tag,cached.tri_tag,cached.tet_tag)
        cell_tags=(copy(old_tags[1]),copy(old_tags[2]),copy(old_tags[3]))
        replacement_cells=blocks[dimension]
        for (new_slot,column) in enumerate(reordered)
            replacement_cells[:,positions[new_slot]]=cells[:,column]
            cell_tags[dimension][positions[new_slot]]=old_tags[dimension][column]
        end
        replacement=Mesh(cached.coords;segs=blocks[1],tris=blocks[2],
                         tets=blocks[3],seg_tag=cell_tags[1],
                         tri_tag=cell_tags[2],tet_tag=cell_tags[3])
        _replace_mesh_cache_locked!(replacement,
            _MeshClassification(replacement,class.entity,
                class.node_entities,class.boundaries,class.seg_entities,
                class.tri_entities,class.tet_entities))
        return nothing
    end
end

function _remove_embedded(dim_tags,dim=-1)
    caller="API.mesh.remove_embedded"
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    return _with_model(invalidate=true) do current
        remove_embedded!(current,dim_tags,dim)
    end
end

# `optimize` maps Gmsh's default tetrahedral mesh optimizer onto the validated
# boundary-preserving `smooth_optimize` kernel. Connectivity and the
# index-aligned classification snapshot carry over unchanged because only node
# coordinates move and boundary nodes are fixed. Other optimizer names and
# entity-scoped selections are not implemented and fail explicitly.
function _optimize_mesh(method="",force=false,niter=1,dim_tags=())
    caller="API.mesh.optimize"
    method isa AbstractString || throw(ArgumentError(
        "$caller: method must be a string"))
    method=="" || throw(ArgumentError(
        "$caller: unknown or unsupported optimizer \"$method\"; only the " *
        "default tetrahedral optimizer is implemented"))
    force isa Bool || throw(ArgumentError("$caller: force must be Bool"))
    iterations=_mesh_query_integer(niter,caller,"niter")
    iterations>=0 || throw(ArgumentError("$caller: niter must be nonnegative"))
    pairs=_mesh_parse_dim_tags(dim_tags,caller)
    return lock(STATE_LOCK) do
        _model_locked()
        cached=_cached_mesh_locked(caller)
        isempty(pairs) || throw(ArgumentError(
            "$caller: entity-scoped optimization is not implemented; pass an " *
            "empty dim_tags selection"))
        class=_cached_classification_locked(cached)
        smoothed=smooth_optimize(cached;iters=iterations)
        new_class=class===nothing ? nothing : _MeshClassification(
            smoothed,class.entity,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(smoothed,new_class)
        return nothing
    end
end

# Per-element visibility is display state only: Gmsh 4.15.2 stores the raw
# value (default 1), accepts unknown tags silently in both directions, and
# reports 0 for them on read.
function _set_mesh_visibility(element_tags,value)
    caller="API.mesh.set_visibility"
    (element_tags isa AbstractVector || element_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: element_tags must be a vector or tuple of element tags"))
    flag=_mesh_query_integer(value,caller,"value")
    typemin(Int32)<=flag<=typemax(Int32) || throw(ArgumentError(
        "$caller: value $value is out of range"))
    return lock(STATE_LOCK) do
        _model_locked()
        cached=LAST_MESH[]
        total=cached===nothing ? 0 : _mesh_element_offsets(cached)[3]
        for raw in element_tags
            tag=_mesh_query_integer(raw,caller,"element_tags entry")
            # Unknown element tags are dropped like Gmsh, so a later
            # `get_visibility` still reports 0 for them.
            tag in 1:total && (ELEMENT_VISIBILITY[][tag]=Int32(flag))
        end
        return nothing
    end
end

function _get_mesh_visibility(element_tags)
    caller="API.mesh.get_visibility"
    (element_tags isa AbstractVector || element_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: element_tags must be a vector or tuple of element tags"))
    return lock(STATE_LOCK) do
        _model_locked()
        cached=LAST_MESH[]
        total=cached===nothing ? 0 : _mesh_element_offsets(cached)[3]
        values=Int32[]
        for raw in element_tags
            tag=_mesh_query_integer(raw,caller,"element_tags entry")
            push!(values,get(ELEMENT_VISIBILITY[],tag,
                             tag in 1:total ? Int32(1) : Int32(0)))
        end
        return values
    end
end

# `getPeriodic` maps each queried entity to its periodic master tag, or to the
# entity itself when no periodic relation was registered. Unknown entities fail
# like Gmsh's "does not exist" error.
function _get_periodic(dim,tags)
    caller="API.mesh.get_periodic"
    dimension=_mesh_query_integer(dim,caller,"dim")
    dimension in 0:3 || throw(ArgumentError("$caller: dim must be in 0:3"))
    (tags isa AbstractVector || tags isa Tuple) || throw(ArgumentError(
        "$caller: tags must be a vector or tuple of entity tags"))
    return _with_model() do current
        dictionary=_mesh_entity_dictionary(current,dimension)
        masters=Int32[]
        for tag in tags
            (tag isa Integer && !(tag isa Bool) &&
             typemin(Int32)<=tag<=typemax(Int32) &&
             haskey(dictionary,Int(tag))) || throw(ArgumentError(
                "$caller: unknown " *
                "$(_MESH_ENTITY_LABELS[dimension+1])[$tag]"))
            constraint=get(current.periodic,(dimension,Int(tag)),nothing)
            push!(masters,constraint===nothing ? Int32(tag) :
                  constraint.master_entity)
        end
        return masters
    end
end

# Gmsh's `removeConstraints` clears per-entity meshing attributes (transfinite,
# recombine, smoothing, reverse); periodic relations, embeddings, and Point
# mesh sizes are retained — verified against 4.15.2. The native model stores
# none of the cleared attribute kinds, so this is a validated no-op.
function _remove_constraints(dim_tags=())
    caller="API.mesh.remove_constraints"
    pairs=_mesh_parse_dim_tags(dim_tags,caller)
    return _with_model() do current
        for (dimension,tag) in pairs
            haskey(_mesh_entity_dictionary(current,dimension),tag) ||
                throw(ArgumentError(
                    "$caller: unknown " *
                    "$(_MESH_ENTITY_LABELS[dimension+1])[$tag]"))
        end
        return nothing
    end
end

# Reverse Cuthill-McKee over the graph of nodes sharing a selected element.
# Components are processed from the minimum-degree unvisited node and each
# breadth-first frontier enqueues neighbors in increasing degree order; the
# reversed traversal assigns the new dense tags.
function _mesh_rcm_node_order(adjacency)
    count=length(adjacency)
    degree=[length(neighbors) for neighbors in adjacency]
    visited=falses(count)
    order=Int[]
    queue=Int[]
    while length(order)<count
        root=0
        best=typemax(Int)
        for node in 1:count
            if !visited[node] && degree[node]<best
                best=degree[node];root=node
            end
        end
        visited[root]=true
        push!(order,root)
        empty!(queue);push!(queue,root)
        head=1
        while head<=length(queue)
            node=queue[head];head+=1
            pending=sort!(Int[neighbor for neighbor in adjacency[node]
                              if !visited[neighbor]];
                          by=neighbor->degree[neighbor])
            for neighbor in pending
                visited[neighbor]=true
                push!(order,neighbor)
                push!(queue,neighbor)
            end
        end
    end
    return reverse!(order)
end

function _compute_renumbering(method="RCMK",element_tags=())
    caller="API.mesh.compute_renumbering"
    method isa AbstractString || throw(ArgumentError(
        "$caller: method must be a string"))
    method=="RCMK" || throw(ArgumentError(
        "$caller: unknown renumbering method \"$method\""))
    (element_tags isa AbstractVector || element_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: element_tags must be a vector or tuple of element tags"))
    return lock(STATE_LOCK) do
        _model_locked()
        cached=_cached_mesh_locked(caller)
        _,_,total=_mesh_element_offsets(cached)
        selected=Set{Int}()
        for value in element_tags
            tag=_mesh_query_integer(value,caller,"element_tags entry")
            (1<=tag<=total) || throw(ArgumentError(
                "$caller: unknown element $tag"))
            push!(selected,tag)
        end
        triangle_offset,tetrahedron_offset,_=_mesh_element_offsets(cached)
        isempty(element_tags) &&
            (selected=Set{Int}(1:total))
        cells=Tuple{Int,Matrix{Int32}}[]
        size(cached.segs,2)>0 && push!(cells,(0,cached.segs))
        size(cached.tris,2)>0 &&
            push!(cells,(triangle_offset,cached.tris))
        size(cached.tets,2)>0 &&
            push!(cells,(tetrahedron_offset,cached.tets))
        involved=Set{Int32}()
        for (offset,block) in cells
            for column in axes(block,2)
                (offset+column) in selected || continue
                for row in axes(block,1)
                    push!(involved,block[row,column])
                end
            end
        end
        old_tags=sort!(collect(involved))
        local_index=Dict{Int32,Int}(
            node=>index for (index,node) in enumerate(old_tags))
        adjacency=[Set{Int}() for _ in old_tags]
        for (offset,block) in cells
            for column in axes(block,2)
                (offset+column) in selected || continue
                vertices=[local_index[block[row,column]]
                          for row in axes(block,1)]
                for a in vertices,b in vertices
                    a!=b && push!(adjacency[a],b)
                end
            end
        end
        order=_mesh_rcm_node_order(adjacency)
        position=Vector{Int}(undef,length(old_tags))
        for (rank,node) in enumerate(order)
            position[node]=rank
        end
        return UInt64.(old_tags),
               UInt64[position[local_index[node]] for node in old_tags]
    end
end

# The native model owns no mesh partitions: ghosts and unpartitioning are
# empty no-ops, and meshing failures surface through the raised error rather
# than through post-hoc error lists.
function _get_ghost_elements(dim,tag)
    caller="API.mesh.get_ghost_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        entity=_mesh_query_integer(tag,caller,"tag")
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: ghost queries require mesh classification metadata"))
        _mesh_classified_entity(model,class,dimension,entity,caller)
        return UInt64[],Int32[]
    end
end

_unpartition()=lock(STATE_LOCK) do
    _model_locked()
    _cached_mesh_locked("API.mesh.unpartition")
    return nothing
end

_get_last_entity_error()=Tuple{Int32,Int32}[]
_get_last_node_error()=UInt64[]

function _rebuild_node_cache(only_if_necessary=true)
    caller="API.mesh.rebuild_node_cache"
    _mesh_query_bool(only_if_necessary,caller,"only_if_necessary")
    lock(STATE_LOCK) do
        _model_locked()
        _cached_mesh_locked(caller)
    end
    return nothing
end

function _rebuild_element_cache(only_if_necessary=true)
    caller="API.mesh.rebuild_element_cache"
    _mesh_query_bool(only_if_necessary,caller,"only_if_necessary")
    lock(STATE_LOCK) do
        _model_locked()
        _cached_mesh_locked(caller)
    end
    return nothing
end

# Node ownership is derived from the generating entity at classification time
# and maintained through every mutation, so geometric relocation and
# element-based reclassification are identity operations; both validate their
# arguments and the session state like Gmsh.
function _relocate_nodes(dim=-1,tag=-1)
    caller="API.mesh.relocate_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_integer(dim,caller,"dim")
        (-1<=dimension<=3) || throw(ArgumentError(
            "$caller: dim must be -1 or in 0:3"))
        entity=_mesh_query_integer(tag,caller,"tag")
        class=_cached_classification_locked(cached)
        entity>=0 || return nothing
        if entity>0
            if dimension<0
                throw(ArgumentError(
                    "$caller: a nonnegative tag requires a dimension in 0:3"))
            end
            class===nothing && throw(ArgumentError(
                "$caller: entity-selective relocation requires mesh " *
                "classification metadata"))
            _mesh_classified_entity(model,class,dimension,entity,caller)
        end
        return nothing
    end
end

function _reclassify_nodes()
    caller="API.mesh.reclassify_nodes"
    return lock(STATE_LOCK) do
        _model_locked()
        _cached_mesh_locked(caller)
        return nothing
    end
end

function _set_size(dim_tags,size)
    caller="API.mesh.set_size"
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    point_tags=Any[]
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        dimension=pair[1]
        dimension isa Integer || throw(ArgumentError(
            "$caller: entity dimensions must be integers"))
        dimension isa Bool && throw(ArgumentError(
            "$caller: entity dimensions must not be Bool"))
        dimension==0 || throw(ArgumentError(
            "$caller: only dimension-0 Point entities are supported"))
        push!(point_tags,pair[2])
    end
    return _with_model(invalidate=true) do current
        set_point_mesh_size!(current,point_tags,size)
    end
end

function _set_periodic(dim,slave_entities,master_entities,affine;atol=1e-12)
    return _with_model(invalidate=true) do current
        set_periodic!(
            current,dim,slave_entities,master_entities,affine;atol=atol)
    end
end

function _get_periodic_nodes(dim,slave_entity)
    return lock(STATE_LOCK) do
        current=_model_locked()
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "API.mesh.get_periodic_nodes: no mesh"))
        model_periodic_nodes(current,cached,dim,slave_entity)
    end
end

"""Gmsh-style mesh generation, mutation, bulk retrieval, and periodic operations."""
module mesh
using ..API: _generate,_get_mesh,_get_nodes,_get_node,
             _get_nodes_for_physical_group,_get_embedded,_get_sizes,
             _get_elements,_get_element,
             _get_element_types,
             _get_elements_by_type,_get_nodes_by_element_type,_get_barycenters,
             _get_element_edge_nodes,_get_element_face_nodes,
             _get_element_type,_get_element_properties,
             _get_integration_points,
             _get_element_by_coordinates,_get_elements_by_coordinates,
             _get_local_coordinates_in_element,_get_element_qualities,
             _get_jacobians,_get_jacobian,
             _get_basis_functions,_get_basis_functions_orientation,
             _get_basis_functions_orientation_for_element,
             _get_number_of_orientations,_get_keys,_get_keys_for_element,
             _get_number_of_keys,_get_keys_information,
             _create_edges,_create_faces,_get_edges,_get_faces,
             _get_all_edges,_get_all_faces,_add_edges,_add_faces,
             _get_max_node_tag,_get_max_element_tag,
             _refine,_clear_mesh,_affine_transform_mesh,_set_size,_set_periodic,
             _get_periodic_nodes,_remove_elements,_reverse_mesh,
             _reverse_elements,_get_duplicate_nodes,
             _remove_duplicate_nodes,_remove_duplicate_elements,
             _set_node,_renumber_nodes,_renumber_elements,_reorder_elements,
             _remove_embedded,
             _get_ghost_elements,_unpartition,_get_last_entity_error,
             _get_last_node_error,_rebuild_node_cache,_rebuild_element_cache,
             _reclassify_nodes,_relocate_nodes,
             _get_periodic,_remove_constraints,_compute_renumbering,
             _optimize_mesh,_set_mesh_visibility,_get_mesh_visibility
generate(dim::Integer)=_generate(dim)
get()=_get_mesh()

"""
    get_nodes(dim=-1, tag=-1, include_boundary=false,
              return_parametric_coord=true)

Return detached dense `UInt64` node tags, flattened `Float64` coordinates, and
flattened parametric coordinates for the cached mesh. `dim=-1` returns every
cached node and ignores `tag`. A nonnegative `dim` selects nodes classified on
model entities of that dimension — all such entities for a negative `tag`, or the
single `(dim, tag)` entity otherwise. With `include_boundary=true`, nodes
classified on the entity's transitive boundary entities are appended after its
own, matching Gmsh's per-entity emission (boundary nodes can repeat across
entities). With `return_parametric_coord=true` each returned node is followed by
its parameters on the queried entity — one `u` per node for a Line and `(u, v)`
per node for a Plane — matching Gmsh 4.15.2's per-entity parametrization;
Points, Volumes, and `dim=-1` queries own no parametrization and return an
empty parametric vector.
"""
get_nodes(dim=-1,tag=-1,include_boundary=false,return_parametric_coord=true)=
    _get_nodes(dim,tag,include_boundary,return_parametric_coord)

"""
    get_node(node_tag)

Return `(coordinates, parametric_coordinates, entity_dimension, entity_tag)`
for one dense cached node tag, matching Gmsh 4.15.2's `getNode` result order.
Coordinates are a detached 3-vector; parametric coordinates are the node's
parameters on its owning entity — one `u` for a Line owner, `(u, v)` for a
Plane owner, and none for Point or Volume owners. The entity fields come from
the classification snapshot built when the mesh was generated; unknown tags
and caches without classification fail explicitly.
"""
get_node(node_tag)=_get_node(node_tag)

"""
    get_nodes_for_physical_group(dim, tag)

Return detached dense `UInt64` node tags and flattened `Float64` coordinates
for every node on the member entities of Physical group `(dim, tag)`, matching
Gmsh 4.15.2's `getNodesForPhysicalGroup`. Each member contributes its own
classified nodes plus its transitive boundary and embedded entities' nodes;
the result is the deduplicated ascending tag set. An unknown or empty group
returns empty arrays. Caches without classification fail explicitly.
"""
get_nodes_for_physical_group(dim,tag)=_get_nodes_for_physical_group(dim,tag)

"""
    get_embedded(dim, tag)

Return the detached `(dim, tag)` list of entities embedded in entity
`(dim, tag)`, matching Gmsh 4.15.2's `getEmbedded`. Entities without
embeddings return an empty list; unknown entities fail explicitly.
"""
get_embedded(dim,tag)=_get_embedded(dim,tag)

"""
    get_sizes(dim_tags)

Return the detached `Float64` mesh size for each `(dim, tag)` pair, matching
Gmsh 4.15.2's `getSizes`: Points report their assigned mesh size and every
other entity — known, unknown, or out-of-range alike — reports `0.0`. Malformed
`dim_tags` entries fail explicitly.
"""
get_sizes(dim_tags)=_get_sizes(dim_tags)

"""
    get_elements(dim=-1, tag=-1)

Return detached Gmsh-shaped `(element_types, element_tags, node_tags)` arrays for
linear segments (type 1), triangles (type 2), and tetrahedra (type 4) in the
cached mesh. A dimension in `0:3` with negative `tag` filters whole type blocks.
Element tags are dense identifiers derived for the current cache across blocks in
that order. Global queries use `dim=-1`, which ignores `tag`. A nonnegative `tag`
selects only the elements classified on the `(dim, tag)` model entity; unknown
entities fail explicitly. Values outside `-1` or `0:3` are rejected.
"""
get_elements(dim=-1,tag=-1)=_get_elements(dim,tag)

"""
    get_element(element_tag)

Return `(element_type, node_tags, entity_dimension, entity_tag)` for one dense
cached element tag, matching Gmsh 4.15.2's `getElement` result order. The cache
holds only linear-simplex types 1, 2, and 4. The entity fields come from the
classification snapshot built when the mesh was generated; unknown tags and
caches without classification fail explicitly.
"""
get_element(element_tag)=_get_element(element_tag)

"""Return the detached linear-simplex MSH types present in a whole cache/dimension
or, for a nonnegative `tag`, on the `(dim, tag)` model entity."""
get_element_types(dim=-1,tag=-1)=_get_element_types(dim,tag)

"""
    get_element_type(family_name, order, serendip=false)

Return the Gmsh 4.15.2 numeric type for a canonical element family and polynomial
order. Family names are case-insensitive. A requested serendipity type falls back
to the complete type when no distinct incomplete element exists at that order.
Unknown families and unsupported orders fail explicitly instead of returning zero.
"""
get_element_type(family_name,order,serendip=false)=
    _get_element_type(family_name,order,serendip)

"""
    get_element_properties(element_type)

Return `(name, dim, order, num_nodes, local_node_coordinates,
num_primary_nodes)` for one Gmsh element type. Local coordinates are detached and
follow Gmsh node ordering. Verified high-order-prism and trihedron layouts remain
available even where Gmsh 4.15.2's own property call fails.
"""
get_element_properties(element_type)=_get_element_properties(element_type)

"""
    get_integration_points(element_type, integration_type)

Return detached reference coordinates and weights for a fixed-node Point, Line,
Triangle, Quadrangle, Tetrahedron, Hexahedron, Prism, or Pyramid type.
Coordinates are flattened `(u,v,w)` triples. `GaussN` preserves Gmsh 4.15.2's
economical rules, including Triangle through order 20 and Tetrahedron through
order 21, then uses the same tensor transitions as Gmsh. Prism rules combine
the corresponding Triangle and Line rules. `CompositeGaussN` selects bounded
tensor and Duffy rules directly. Trihedra have no integration rule in the pinned
Gmsh release. An omitted `N` means order zero. This query does not require a
model or mesh.
"""
get_integration_points(element_type,integration_type)=
    _get_integration_points(element_type,integration_type)

"""
    get_element_by_coordinates(x, y, z, dim=-1, strict=false)

Return `(element_tag, element_type, node_tags, u, v, w)` for the first cached
linear-simplex element at a finite point. With `dim=-1`, the deterministic first
match has the greatest dimension and then the smallest dense tag. `strict=true`
uses the pinned Gmsh 4.15.2 reference tolerance of `1e-6`; relaxed search widens
that tolerance by decades through `1.0` and stops at the first nonempty level.
No match is an error.
"""
get_element_by_coordinates(x,y,z,dim=-1,strict=false)=
    _get_element_by_coordinates(x,y,z,dim,strict)

"""
    get_elements_by_coordinates(x, y, z, dim=-1, strict=false)

Return all matching dense element tags, ordered by decreasing dimension and then
increasing tag. Search tolerance and no-match behavior are identical to
[`get_element_by_coordinates`](@ref).
"""
get_elements_by_coordinates(x,y,z,dim=-1,strict=false)=
    _get_elements_by_coordinates(x,y,z,dim,strict)

"""
    get_local_coordinates_in_element(element_tag, x, y, z)

Return `(u,v,w)` reference coordinates for one dense cached element tag. Segment
coordinates use `u ∈ [-1,1]`; triangle and tetrahedron coordinates use the standard
unit simplex. Segment and triangle points are orthogonally projected onto the
element span, with unused coordinates returned as exact zeros. The point need not
lie inside the element. A Float64-unrepresentable result fails explicitly.
"""
get_local_coordinates_in_element(element_tag,x,y,z)=
    _get_local_coordinates_in_element(element_tag,x,y,z)

"""
    get_element_qualities(element_tags, quality_name="minSICN",
                          task=0, num_tasks=1)

Return detached `Float64` qualities for dense cached linear-simplex element tags,
in request order. Names follow the documented Gmsh 4.15.2 list. Segment tags are
explicitly unsupported for `minDetJac`, `maxDetJac`, `minSIGE`, and
`minIsotropy`, whose 1-D Gmsh implementations are absent or unreliable.
With `num_tasks > 1`, only the contiguous Gmsh block of requested tags with
0-based positions `begin = (task*count) ÷ num_tasks` through
`end = ((task+1)*count) ÷ num_tasks` is evaluated and returned; unlike Gmsh's
preallocated C++ output, the detached result contains exactly that slice and no
zero padding. Tag validation is slice-scoped, matching Gmsh, which leaves
entries outside the block untouched. `task >= num_tasks` returns an empty
vector without error; a negative `task`, `num_tasks < 1`, or a non-integer
task count fails explicitly. Gmsh 4.15.2 cannot be cross-checked for
`task >= num_tasks` here: its quality loop indexes the request vector
unguarded, so that case reads out of bounds (observed `Unknown element 0`)
where the other six partitioned queries use a guarded range and stay empty.
Tessella deterministically returns the empty slice instead of replicating
that out-of-bounds read.
"""
get_element_qualities(element_tags,quality_name="minSICN",task=0,num_tasks=1)=
    _get_element_qualities(element_tags,quality_name,task,num_tasks)

"""
    get_jacobians(element_type, local_coord, tag=-1, task=0, num_tasks=1)

Return detached `(jacobians, determinants, coordinates)` for every cached element
of one linear-simplex Gmsh type at concatenated `(u,v,w)` points. Results are
ordered by element and then point; each 3×3 Jacobian is flattened by column.
Segment determinants are positive lengths, triangle determinants are positive
area scales, and tetrahedron determinants retain orientation. The cache supports
types 1, 2, and 4. A nonnegative `tag` selects only the elements classified on
the entity of `tag` in the type's own dimension; unknown entities fail
explicitly. With `num_tasks > 1`, only the contiguous Gmsh block of cached elements
with 0-based positions `begin = (task*count) ÷ num_tasks` through
`end = ((task+1)*count) ÷ num_tasks` is evaluated and returned; unlike Gmsh's
preallocated C++ output, the detached result contains exactly that slice and no
zero padding. `task >= num_tasks` returns three empty vectors without error;
a negative `task`, `num_tasks < 1`, or a non-integer task count fails
explicitly.
"""
get_jacobians(element_type,local_coord,tag=-1,task=0,num_tasks=1)=
    _get_jacobians(element_type,local_coord,tag,task,num_tasks)

"""
    get_jacobian(element_tag, local_coord)

Return detached Jacobians, determinants, and mapped physical coordinates for one
dense cached linear-simplex element tag. Layout, reference coordinates, numerical
contracts, and error behavior match [`get_jacobians`](@ref).
"""
get_jacobian(element_tag,local_coord)=
    _get_jacobian(element_tag,local_coord)

"""
    get_basis_functions(element_type, local_coord, function_space_type,
                        wanted_orientations=Int32[])

Return `(num_components, basis_functions, num_orientations)` at concatenated
`(u,v,w)` evaluation points. Unqualified Lagrange and isoparametric aliases use the
input fixed type's actual nodal order and completeness. `LagrangeN` and
`GradLagrangeN` select the complete family basis at order `N`; the catalog covers
orders 0--10, with Hexahedron, Prism, and Pyramid ending at order 9. Order-one
hierarchical H1 functions and gradients cover types 1, 2, and 4, every fixed
Quadrangle, Hexahedron, and Prism type at any Lagrange order, and Point type 15;
their values repeat the reference family's vertex functions once per orientation.
Lowest-order H(curl) functions and curls cover types 1, 2, and 4. Pyramid and
Trihedron hierarchical spaces are rejected: Gmsh 4.15.2 defines no Pyramid
hierarchical family and no Trihedron basis. Hexahedron H1 orientation counts use
8!, which the pinned release's `getBasisFunctions` count confirms; its
`getNumberOfOrientations` metadata query reads uninitialized memory there
(different values across processes) and is not copied. Values use Gmsh's orientation-
then-point-then-function-then-component layout. An empty orientation selection
returns every hierarchical orientation or the sole nodal orientation. This query
does not require a cached mesh.
"""
get_basis_functions(element_type,local_coord,function_space_type,
                    wanted_orientations=Int32[])=
    _get_basis_functions(
        element_type,local_coord,function_space_type,wanted_orientations)

"""
    get_number_of_orientations(element_type, function_space_type)

Return one for Lagrange spaces or the factorial primary-vertex orientation count
for a supported first-order hierarchical space: 1 for Point H1, 2/6/24 for
linear simplexes, 24 for Quadrangle H1, 720 for Prism H1, and 8! for Hexahedron
H1. The Hexahedron count follows the pinned release's verified basis count;
its `getNumberOfOrientations` metadata query reads uninitialized memory
(values differ across processes).
This reference-element query does
not require a cached mesh.
"""
get_number_of_orientations(element_type,function_space_type)=
    _get_number_of_orientations(element_type,function_space_type)

"""
    get_basis_functions_orientation(element_type, function_space_type,
                                    tag=-1, task=0, num_tasks=1)

Return one lexicographic orientation index per cached element of the requested
supported type. Nodal spaces return zeros. Known fixed types absent from the
linear-simplex cache return an empty vector. A nonnegative `tag` selects only
the elements classified on the entity of `tag` in the type's own dimension;
unknown entities fail explicitly. With `num_tasks > 1`, only the
contiguous Gmsh block of cached elements with 0-based positions
`begin = (task*count) ÷ num_tasks` through `end = ((task+1)*count) ÷ num_tasks`
is evaluated and returned; unlike Gmsh's preallocated C++ output, the detached
result contains exactly that slice and no zero padding. `task >= num_tasks`
returns an empty vector without error; a negative `task`, `num_tasks < 1`, or a
non-integer task count fails explicitly. Gmsh 4.15.2 cannot be cross-checked
for `task >= num_tasks` here: its hierarchical orientation loop indexes
per-entity elements unguarded, so that case segfaults the pinned release
where the guarded queries stay empty. Tessella deterministically returns the
empty slice instead.
For hierarchical spaces Gmsh partitions each entity separately while this
entity-free cache partitions the global type block; the two coincide whenever
each queried type lives on one entity.
"""
get_basis_functions_orientation(element_type,function_space_type,
                                tag=-1,task=0,num_tasks=1)=
    _get_basis_functions_orientation(
        element_type,function_space_type,tag,task,num_tasks)

"""
    get_basis_functions_orientation_for_element(element_tag, function_space_type)

Return the lexicographic orientation index for one dense cached element tag.
"""
get_basis_functions_orientation_for_element(element_tag,function_space_type)=
    _get_basis_functions_orientation_for_element(
        element_tag,function_space_type)

"""
    get_keys(element_type, function_space_type, tag=-1, return_coord=true)

Return detached `(type_keys, entity_keys, coordinates)` for every cached element
of one linear-simplex type. Lagrange and order-one H1 keys use dense node tags;
lowest-order H(curl) keys use stable global edge tags and lazily add only the edges
visited by the requested type. Coordinates locate node or edge-midpoint keys and
are omitted when `return_coord=false`. A numeric Lagrange space must match the
stored interpolation-node count; the cache does not synthesize higher-order keys.
A nonnegative `tag` selects only the elements classified on the entity of `tag`
in the type's own dimension; unknown entities fail explicitly.
"""
get_keys(element_type,function_space_type,tag=-1,return_coord=true)=
    _get_keys(element_type,function_space_type,tag,return_coord)

"""
    get_keys_for_element(element_tag, function_space_type, return_coord=true)

Return detached node or edge keys for one dense cached element. Lowest-order
H(curl) calls lazily add only that element's missing edges to the shared catalog.
A numeric Lagrange space must match the stored interpolation-node count.
"""
get_keys_for_element(element_tag,function_space_type,return_coord=true)=
    _get_keys_for_element(element_tag,function_space_type,return_coord)

"""
    get_number_of_keys(element_type, function_space_type)

Return the number of node or edge keys owned by one supported reference element.
Actual- and explicit-order nodal counts cover every fixed family except Trihedron;
hierarchical H1 counts cover linear simplexes, Point type 15, and every fixed
Quadrangle, Hexahedron, and Prism type with the reference family's vertex count.
Lowest-order H(curl) counts cover linear simplexes. This query does not require a cache.
"""
get_number_of_keys(element_type,function_space_type)=
    _get_number_of_keys(element_type,function_space_type)

"""
    get_keys_information(type_keys, entity_keys, element_type,
                         function_space_type)

Return `(entity_dimension, polynomial_order)` for complete element-sized groups
of supported node or edge keys. Actual- and explicit-order nodal metadata covers
every fixed family except Trihedron; hierarchical H1 metadata covers linear
simplexes, Point type 15, and every fixed Quadrangle, Hexahedron, and Prism type.
Point H1 keys report order 0 since the Point vertex function is constant.
Lowest-order H(curl) metadata covers linear simplexes. Key arrays must have equal lengths and the expected type-key value for
the selected space.
"""
get_keys_information(type_keys,entity_keys,element_type,function_space_type)=
    _get_keys_information(
        type_keys,entity_keys,element_type,function_space_type)

"""
    create_edges(dim_tags=())

Create deterministic global identifiers for every unique edge in the cached
linear-simplex mesh, or only for the cells owned by the entities in `dim_tags`.
An empty vector or tuple selects the complete cache. A nonempty selection
requires the classified cache built by [`generate`](@ref); unknown entities
fail explicitly, and listed entities owning no cells contribute nothing, as in
Gmsh 4.15.2. Repeated calls are idempotent.
Edges previously attached with [`add_edges`](@ref) are preserved. Replacing,
refining, transforming, or clearing the cache invalidates the catalog. Automatic
identifiers begin at the current edge count plus one and skip identifiers already
in use.
"""
create_edges(dim_tags=())=_create_edges(dim_tags)

"""
    create_faces(dim_tags=())

Create deterministic global identifiers for every missing triangular face in the
cached linear-simplex mesh while preserving faces attached with
[`add_faces`](@ref), including quadrangles. Lifecycle and entity-selection behavior
match [`create_edges`](@ref). Automatic identifiers use the combined triangle and
quadrangle count and skip identifiers already in use.
"""
create_faces(dim_tags=())=_create_faces(dim_tags)

"""
    add_edges(edge_tags, edge_nodes)

Atomically attach global edges to explicit positive identifiers. `edge_nodes`
contains one node pair per identifier, and each identifier belongs to one edge.
Repeating the same edge/tag association is idempotent; zero or conflicting tags,
repeated or unknown nodes, and malformed batches leave the catalog unchanged. A
later [`create_edges`](@ref) preserves these edges and adds missing simplex edges.
"""
add_edges(edge_tags,edge_nodes)=_add_edges(edge_tags,edge_nodes)

"""
    add_faces(face_type, face_tags, face_nodes)

Atomically attach triangular (`face_type=3`) or quadrangular (`face_type=4`) faces
to explicit positive identifiers shared by both face types. `face_nodes` contains
one complete node group per identifier. Repeating the same face/tag association is
idempotent; zero or conflicting tags, repeated or unknown nodes, and malformed
batches leave the catalog unchanged. A later [`create_faces`](@ref) preserves these
faces and adds missing simplex triangles.
"""
add_faces(face_type,face_tags,face_nodes)=
    _add_faces(face_type,face_tags,face_nodes)

"""
    get_edges(node_tags)

Return detached global edge tags and `Int32` orientations for concatenated node
pairs. Positive orientation follows ascending node tags; reversing a pair returns
the same edge tag and orientation `-1`. The edge catalog must first be populated
with [`add_edges`](@ref) or [`create_edges`](@ref). Incomplete pairs and unknown
nodes or edges fail explicitly.
"""
get_edges(node_tags)=_get_edges(node_tags)

"""
    get_faces(face_type, node_tags)

Return detached global face tags and orientations for concatenated triangular or
quadrangular node groups. Face orientations are zero, matching Gmsh 4.15.2's
public face-topology behavior. The face catalog must first be populated with
[`add_faces`](@ref) or [`create_faces`](@ref); automatic creation adds only
triangles because the cached mesh is linear-simplex.
"""
get_faces(face_type,node_tags)=_get_faces(face_type,node_tags)

"""
    get_all_edges()

Return detached global edge tags and first-encounter node pairs in ascending tag
order. Before either [`add_edges`](@ref) or [`create_edges`](@ref), both arrays are
empty.
"""
get_all_edges()=_get_all_edges()

"""
    get_all_faces(face_type)

Return detached global face tags and first-encounter face nodes in ascending tag
order. Before either [`add_faces`](@ref) or [`create_faces`](@ref), both arrays are
empty. Quadrangular results contain only explicitly attached faces.
"""
get_all_faces(face_type)=_get_all_faces(face_type)

"""
    get_elements_by_type(element_type, tag=-1, task=0, num_tasks=1)

Return detached dense element tags and flattened node tags for one fixed-node Gmsh
element type. The simplex cache can contain only types 1, 2, and 4; another known
fixed-node type returns empty arrays. A nonnegative `tag` selects only the
elements classified on the entity of `tag` in the type's own dimension; unknown
entities fail explicitly.
With `num_tasks > 1`, only the contiguous Gmsh block of cached elements with
0-based positions `begin = (task*count) ÷ num_tasks` through
`end = ((task+1)*count) ÷ num_tasks` is returned; unlike Gmsh's preallocated
C++ output, the detached result contains exactly that slice and no zero
padding. `task >= num_tasks` returns two empty vectors without error; a
negative `task`, `num_tasks < 1`, or a non-integer task count fails explicitly.
"""
get_elements_by_type(element_type,tag=-1,task=0,num_tasks=1)=
    _get_elements_by_type(element_type,tag,task,num_tasks)

"""
    get_nodes_by_element_type(element_type, tag=-1,
                              return_parametric_coord=true)

Return detached node tags and coordinates in per-element connectivity order for
one fixed-node Gmsh element type. Shared nodes consequently appear once per element
use. A nonnegative `tag` selects only the elements classified on the entity of
`tag` in the type's own dimension; unknown entities fail explicitly. With
`return_parametric_coord=true` the third result packs each entry's parameters
on its owning entity — one `u` for a Line owner, `(u, v)` for a Plane owner,
and nothing for Point or Volume owners — matching Gmsh 4.15.2's variable-width
emission order; an unclassified cache returns an empty parametric vector.
"""
get_nodes_by_element_type(element_type,tag=-1,return_parametric_coord=true)=
    _get_nodes_by_element_type(element_type,tag,return_parametric_coord)

"""
    get_barycenters(element_type, tag, fast, primary,
                    task=0, num_tasks=1)

Return detached `x,y,z` barycenters in element order for a cached linear-simplex
type. With `fast=true`, return unnormalized primary-node coordinate sums. All nodes
are primary for types 1, 2, and 4. A nonnegative `tag` selects only the elements
classified on the entity of `tag` in the type's own dimension; unknown entities
fail explicitly.
With `num_tasks > 1`, only the contiguous Gmsh block of cached elements with
0-based positions `begin = (task*count) ÷ num_tasks` through
`end = ((task+1)*count) ÷ num_tasks` is returned; unlike Gmsh's preallocated
C++ output, the detached result contains exactly that slice and no zero
padding. `task >= num_tasks` returns an empty vector without error; a negative
`task`, `num_tasks < 1`, or a non-integer task count fails explicitly.
"""
get_barycenters(element_type,tag,fast,primary,task=0,num_tasks=1)=
    _get_barycenters(element_type,tag,fast,primary,task,num_tasks)

"""
    get_element_edge_nodes(element_type, tag=-1, primary=false,
                           task=0, num_tasks=1)

Return detached edge-node tags in Gmsh local-edge order for every cached element of
one type. The `primary` flag is validated but does not change linear-simplex output.
A nonnegative `tag` selects only the elements classified on the entity of `tag`
in the type's own dimension; unknown entities fail explicitly. With
`num_tasks > 1`, only the contiguous Gmsh block of cached elements with 0-based
positions
`begin = (task*count) ÷ num_tasks` through `end = ((task+1)*count) ÷ num_tasks`
is returned; unlike Gmsh's preallocated C++ output, the detached result contains
exactly that slice and no zero padding. `task >= num_tasks` returns an empty
vector without error; a negative `task`, `num_tasks < 1`, or a non-integer task
count fails explicitly.
"""
get_element_edge_nodes(element_type,tag=-1,primary=false,task=0,num_tasks=1)=
    _get_element_edge_nodes(element_type,tag,primary,task,num_tasks)

"""
    get_element_face_nodes(element_type, face_type, tag=-1, primary=false,
                           task=0, num_tasks=1)

Return detached face-node tags in Gmsh local-face order for every cached element of
one type. `face_type` is 3 for triangles or 4 for quadrangles. The `primary` flag is
validated but does not change linear-simplex output. A nonnegative `tag` selects
only the elements classified on the entity of `tag` in the type's own dimension;
unknown entities fail explicitly. With `num_tasks > 1`, only the contiguous Gmsh block of cached
elements with 0-based positions `begin = (task*count) ÷ num_tasks` through
`end = ((task+1)*count) ÷ num_tasks` is returned; unlike Gmsh's preallocated
C++ output, the detached result contains exactly that slice and no zero padding.
`task >= num_tasks` returns an empty vector without error; a negative `task`,
`num_tasks < 1`, or a non-integer task count fails explicitly.
"""
get_element_face_nodes(element_type,face_type,tag=-1,primary=false,
                       task=0,num_tasks=1)=
    _get_element_face_nodes(
        element_type,face_type,tag,primary,task,num_tasks)

"""Return the greatest dense node tag, or zero when the cached mesh has no nodes."""
get_max_node_tag()=_get_max_node_tag()

"""Return the greatest dense element tag, or zero when the cached mesh has no cells."""
get_max_element_tag()=_get_max_element_tag()

"""
    refine(; max_nodes=typemax(Int32), max_cells=typemax(Int32)) -> Mesh

Uniformly refine the complete cached linear-simplex mesh once. The cache changes
only after successful validation and resource preflight; the returned mesh owns
storage independently of the cache.
"""
refine(;max_nodes=typemax(Int32),max_cells=typemax(Int32))=
    _refine(;max_nodes=max_nodes,max_cells=max_cells)

"""
    clear(dim_tags=())

Clear the complete cached mesh, or only the mesh data owned by the entities in
`dim_tags`, without changing model geometry. An empty vector or tuple selects
the complete cache. A nonempty selection requires the classified cache built by
[`generate`](@ref): cells classified on each listed entity are removed, nodes
they owned drop when no surviving cell references them, and nodes still
referenced reclassify to the lowest-dimension (then lowest-tag) surviving
owner, mirroring Gmsh 4.15.2's per-entity node storage on the flat shared-node
cache. Listed entities owning no cells — points, or boundary entities of a
higher-dimensional generating entity — are no-ops, matching Gmsh's retention
of boundary meshes and unmeshable vertices. Unknown entities fail explicitly,
and clearing is a no-op when no mesh is cached.
"""
clear(dim_tags=())=_clear_mesh(dim_tags)

"""
    affine_transform(affine, dim_tags=()) -> Mesh

Apply a finite nonsingular affine transform to every node in the complete cached
mesh, or only to the nodes owned by the entities in `dim_tags`. `affine` is a
4×4 matrix or exactly 12 or 16 entries in Gmsh row-major order; 12 entries imply
the homogeneous row `(0, 0, 0, 1)`. An empty vector or tuple selects the
complete cache. A nonempty selection requires the classified cache built by
[`generate`](@ref) and moves only nodes classified on the listed entities, as in
Gmsh 4.15.2's per-entity node storage; shared nodes owned by unlisted boundary
entities stay fixed. Under an orientation-reversing matrix only cells whose
every node moved are rewound, since the orientation of a partially transformed
cell is decided by the resulting geometry. The cache changes only after the
transformed mesh validates — stricter than Gmsh, which applies invalid results
— and the returned mesh owns independent storage. Unknown entities fail
explicitly. Model geometry and periodic relations are unchanged.
"""
affine_transform(affine,dim_tags=())=
    _affine_transform_mesh(affine,dim_tags)

"""
    remove_elements(dim, tag, element_tags=[])

Remove the listed dense element tags — or every element on the entity when
`element_tags` is empty — classified on `(dim, tag)`, matching Gmsh 4.15.2's
`removeElements`. Nodes are retained, so dense element tags re-index after
removal. Every listed tag must resolve to an element owned by the entity;
unknown entities, unknown tags, and caches without classification metadata
fail explicitly. A dimension-0 selection is a validated no-op since the
simplex cache owns no Point cells.
"""
remove_elements(dim,tag,element_tags=())=
    _remove_elements(dim,tag,element_tags)

"""
    reverse(dim_tags=())

Reverse the orientation of every element classified on the entities in
`dim_tags` — or the complete cached mesh when empty — using Gmsh 4.15.2's
first-order simplex conventions: segments swap both vertices, triangles swap
the last two, tetrahedra swap the first two. A nonempty selection requires the
classified cache built by [`generate`](@ref); unknown entities fail explicitly.
"""
reverse(dim_tags=())=_reverse_mesh(dim_tags)

"""
    reverse_elements(element_tags)

Reverse the orientation of the listed dense element tags with the same
first-order simplex conventions as [`reverse`](@ref). Unknown tags fail
explicitly.
"""
reverse_elements(element_tags)=_reverse_elements(element_tags)

"""
    get_duplicate_nodes(dim_tags=()) -> Vector{UInt64}

Return the sorted detached tags of coincident nodes owned by the entities in
`dim_tags` — or the complete cache when empty — matching Gmsh 4.15.2's
`getDuplicateNodes`: every node sharing exact coordinates with another scanned
node is reported. A nonempty selection requires classification metadata; an
empty selection scans the cached coordinates directly.
"""
get_duplicate_nodes(dim_tags=())=_get_duplicate_nodes(dim_tags)

"""
    remove_duplicate_nodes(dim_tags=())

Merge nodes sharing exact coordinates, keeping the lowest tag and remapping all
connectivity, matching Gmsh 4.15.2's `removeDuplicateNodes`. An empty selection
scans the complete cache; a nonempty selection merges only nodes owned by the
listed entities and requires classification metadata.
"""
remove_duplicate_nodes(dim_tags=())=_remove_duplicate_nodes(dim_tags)

"""
    remove_duplicate_elements(dim_tags=())

Drop cells whose sorted connectivity repeats an earlier cell owned by the same
entity, matching Gmsh 4.15.2's `removeDuplicateElements`. An empty selection
scans every block; a nonempty selection scans only cells classified on the
listed entities. Requires the classification snapshot built by
[`generate`](@ref).
"""
remove_duplicate_elements(dim_tags=())=_remove_duplicate_elements(dim_tags)

"""
    set_size(dim_tags, size)

Set a finite, positive mesh-size constraint on the dimension-0 Point entities in
`dim_tags`. Every pair is validated before the active model or cached mesh changes.
"""
set_size(dim_tags,size)=_set_size(dim_tags,size)

"""
    set_periodic(dim, slave_entities, master_entities, affine; atol=1e-12)

Store validated straight-curve (`dim=1`), planar-surface (`dim=2`), or volume
(`dim=3`) relations in
the active model and invalidate any cached mesh. `affine` maps each master entity
to its corresponding slave in Gmsh row-major 4×4 order. Each slave has one master;
masters may be reused, and a slave may become a master in an acyclic dependency
chain. Curves must share a planar surface when meshed. Surfaces must be
affine-equivalent boundaries of one explicit planar-shell volume. Volume
relations are stored and reported but mesh-inert, as in Gmsh 4.15.2.
"""
set_periodic(dim,slave_entities,master_entities,affine;atol=1e-12)=
    _set_periodic(dim,slave_entities,master_entities,affine;atol=atol)

"""
    get_periodic_nodes(dim, slave_entity)

Return the master entity, detached slave/master node arrays, and affine transform
for one curve or planar boundary-surface relation in the cached mesh. Stored
volume relations return the master and affine with empty node arrays.
"""
get_periodic_nodes(dim,slave_entity)=
    _get_periodic_nodes(dim,slave_entity)

"""
    get_periodic(dim, tags)

Return the periodic master tag for each entity tag in `tags`, or the entity's
own tag when it has no periodic master — matching Gmsh 4.15.2's `getPeriodic`.
`dim` must be in 0:3 and every tag must name an existing model entity.
"""
get_periodic(dim,tags)=_get_periodic(dim,tags)

"""
    remove_constraints(dim_tags=())

Validate `dim_tags` and remove per-entity meshing attributes. Gmsh 4.15.2's
`removeConstraints` clears transfinite, recombine, smoothing, and reverse
attributes while retaining periodic relations, embeddings, and Point mesh
sizes; the native model stores none of the cleared kinds, so this is a
validated no-op. Unknown entities fail explicitly.
"""
remove_constraints(dim_tags=())=_remove_constraints(dim_tags)

"""
    compute_renumbering(method="RCMK", element_tags=())

Compute a node renumbering for the cached mesh without applying it, matching
Gmsh 4.15.2's `computeRenumbering`. `element_tags` restricts the computation to
the given dense element tags (all elements when empty). Returns `(old_tags,
new_tags)`: the sorted involved node tags and, for each, its new dense tag
under a reverse Cuthill-McKee ordering of the shared-element adjacency graph.
Only `"RCMK"` is supported; other methods and unknown element tags fail
explicitly. Element renumbering is not available, as in Gmsh 4.15.2.
"""
compute_renumbering(method="RCMK",element_tags=())=
    _compute_renumbering(method,element_tags)

"""
    optimize(method="", force=false, niter=1, dim_tags=())

Optimize the cached mesh in place, matching Gmsh 4.15.2's `optimize`. The empty
default method runs the validated boundary-preserving tetrahedral optimizer
(`iters` sweeps); meshes without tetrahedra are unchanged. Other optimizer
names and nonempty `dim_tags` entity scoping are not implemented and fail
explicitly. Connectivity and entity classification are preserved because only
interior node coordinates move.
"""
optimize(method="",force=false,niter=1,dim_tags=())=
    _optimize_mesh(method,force,niter,dim_tags)

"""
    set_visibility(element_tags, value)

Set the display visibility flag of the listed dense element tags to integer
`value`, matching Gmsh 4.15.2's `mesh.setVisibility`. Unknown tags are dropped
silently; the state resets whenever the mesh cache is replaced.
"""
set_visibility(element_tags,value)=_set_mesh_visibility(element_tags,value)

"""
    get_visibility(element_tags)

Return the stored visibility flag for each listed dense element tag — 1 for
known elements never set, 0 for unknown tags — matching Gmsh 4.15.2's
`mesh.getVisibility`. Without a cached mesh every tag is unknown.
"""
get_visibility(element_tags)=_get_mesh_visibility(element_tags)

"""
    set_node(node_tag, coord, parametric_coord=Float64[])

Set the Cartesian coordinates of the cached node `node_tag`, matching Gmsh
4.15.2's `setNode`. `coord` must hold exactly three finite coordinates. The
cache stores no per-node parametric coordinates — they are recomputed from the
owning entity's geometry — so `parametric_coord` must be empty. Unknown tags
fail explicitly.
"""
set_node(node_tag,coord,parametric_coord=Float64[])=
    _set_node(node_tag,coord,parametric_coord)

"""
    renumber_nodes(old_tags=(), new_tags=())

Renumber cached node tags. With no explicit pair lists this is Gmsh 4.15.2's
"renumber continuously" form — already a no-op on the dense cache. Explicit
`old_tags`/`new_tags` pairs must together keep the tag set a permutation of
`1:nnodes` (unlisted tags keep their values); the permutation reorders node
storage and remaps all connectivity. Sparse target tags are not representable
on the dense cache and fail explicitly.
"""
renumber_nodes(old_tags=(),new_tags=())=
    _renumber_nodes(old_tags,new_tags)

"""
    renumber_elements(old_tags=(), new_tags=())

Renumber cached dense element tags. With no explicit pair lists this is a
no-op, matching Gmsh 4.15.2's continuous form. Explicit pairs must keep the
tag set a permutation of `1:nelements` and cannot move an element across
element-type blocks; violations fail explicitly.
"""
renumber_elements(old_tags=(),new_tags=())=
    _renumber_elements(old_tags,new_tags)

"""
    reorder_elements(element_type, tag, ordering)

Reorder the elements of `element_type` classified on the entity `tag` in that
type's own dimension, matching Gmsh 4.15.2's `reorderElements`: `ordering` is
a 0-based source-position permutation — `ordering[new_position]` names the
element that moves there. Dense tags follow positions after the reorder while
each element keeps its stored tag metadata and entity ownership. Entities
owning no elements of the type, malformed orderings, and unclassified caches
fail explicitly.
"""
reorder_elements(element_type,tag,ordering)=
    _reorder_elements(element_type,tag,ordering)

"""
    remove_embedded(dim_tags, dim=-1)

Remove the embedded entities recorded on the parent entities in `dim_tags`,
matching Gmsh 4.15.2's `removeEmbedded`. `dim` below 0 removes every embedded
dimension; `dim` in 0:2 drops only embedded entities of that dimension.
Unknown or non-Surface/Volume parents fail explicitly and leave the model and
cached mesh unchanged; a successful removal invalidates the cached mesh like
every model mutation.
"""
remove_embedded(dim_tags,dim=-1)=_remove_embedded(dim_tags,dim)

"""
    get_ghost_elements(dim, tag) -> (Vector{UInt64}, Vector{Int32})

Return the ghost element tags and partitions of entity `(dim, tag)`. The
native model owns no mesh partitions, so both results are always empty for an
existing entity, matching Gmsh 4.15.2 on an unpartitioned mesh. Requires the
classified cache built by [`generate`](@ref); unknown entities fail
explicitly.
"""
get_ghost_elements(dim,tag)=_get_ghost_elements(dim,tag)

"""
    unpartition()

No-op parity with Gmsh 4.15.2's `unpartition` on an unpartitioned mesh: the
native model owns no mesh partitions, so there is nothing to remove. Requires
a cached mesh like the Gmsh call.
"""
unpartition()=_unpartition()

"""
    get_last_entity_error() -> Vector{Tuple{Int32,Int32}}

Return the entities where the last meshing error occurred. Meshing failures
surface through the raised error rather than post-hoc state, so the result is
always empty — matching Gmsh 4.15.2 after a successful `generate`.
"""
get_last_entity_error()=_get_last_entity_error()

"""
    get_last_node_error() -> Vector{UInt64}

Return the nodes where the last meshing error occurred. Meshing failures
surface through the raised error rather than post-hoc state, so the result is
always empty — matching Gmsh 4.15.2 after a successful `generate`.
"""
get_last_node_error()=_get_last_node_error()

"""
    rebuild_node_cache(only_if_necessary=true)

No-op parity with Gmsh 4.15.2's `rebuildNodeCache`: the cached mesh keeps its
node index coherent through every mutation, so there is nothing to rebuild.
"""
rebuild_node_cache(only_if_necessary=true)=
    _rebuild_node_cache(only_if_necessary)

"""
    rebuild_element_cache(only_if_necessary=true)

No-op parity with Gmsh 4.15.2's `rebuildElementCache`: the cached mesh keeps
its element index coherent through every mutation, so there is nothing to
rebuild.
"""
rebuild_element_cache(only_if_necessary=true)=
    _rebuild_element_cache(only_if_necessary)

"""
    reclassify_nodes()

No-op parity with Gmsh 4.15.2's `reclassifyNodes`: node ownership is derived
from the generating entity when the cache is built and maintained through
every mutation, so element-based reclassification is the identity.
"""
reclassify_nodes()=_reclassify_nodes()

"""
    relocate_nodes(dim=-1, tag=-1)

No-op parity with Gmsh 4.15.2's `relocateNodes`: node ownership derives from
the generating entity and stays geometry-consistent through every supported
mutation, so geometric relocation is the identity. `dim=-1, tag=-1` validates
the whole cache; a nonnegative `tag` additionally requires `dim` in 0:3 and an
existing classified entity.
"""
relocate_nodes(dim=-1,tag=-1)=_relocate_nodes(dim,tag)
end

"""
    open_geo!(path; mesh_dim=0)

Execute a bounded `.geo` file into the initialized API session. Stored model and
mesh state are detached from the returned execution result. Accepted syntax and
blockers follow [`Tessella.GeoExec.execute_geo`](@ref).
"""
function open_geo!(path::AbstractString; mesh_dim::Integer=0)
    return lock(STATE_LOCK) do
        _model_locked()
        result=execute_geo(path;mesh_dim=mesh_dim)
        stored_model=deepcopy(result.model)
        stored_mesh=result.mesh===nothing ? nothing : _copy_mesh(result.mesh)
        CURRENT[]=stored_model
        _replace_mesh_cache_locked!(stored_mesh)
        result
    end
end

end # module
