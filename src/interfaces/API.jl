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
quadrangular face catalogs. It also owns a reusable robust AABB locator for dense
element-by-coordinate and reference-coordinate queries, plus scale-robust named
quality queries and Gmsh-shaped forward maps/Jacobians over dense cached elements.
Element type and property lookup is
available without a session and delegates to the immutable native catalog. Production
meshing is never delegated to Gmsh.
"""
module API

using ..Model: GeoModel, add_point!, add_line!, add_curve_loop!, add_plane_surface!
using ..Model: add_surface_loop!, add_volume!
using ..Model: add_box!, add_cylinder!, add_sphere!, add_cone!, boolean_volumes!
using ..Model: _remove_volume_entity!
using ..Model: embed!, set_point_mesh_size!
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
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets
using ..MeshEntityTopology: MeshEdgeTopology, MeshFaceTopology,
                            _mesh_edge_topology, _mesh_face_topology,
                            _mesh_add_edges, _mesh_add_faces,
                            _mesh_edges, _mesh_faces,
                            _mesh_all_edges, _mesh_all_faces,
                            _simplex_edge_patterns, _simplex_face_patterns
using ..MeshPointLocation: SimplexLocator, mesh_element_offsets,
                           mesh_element_block, mesh_element_record,
                           _local_coordinates, _locate_elements,
                           _require_local_coordinates
using ..MeshElementQuality: mesh_element_qualities
using ..MeshReferenceGeometry: mesh_jacobian, mesh_jacobians
using ..Elements: msh_spec, msh_type, msh_properties
using ..Refine: refine_uniform
using ..Transform: affine_transform, _transform_gmsh_affine
using ..GeoExec: execute_geo

export initialize, finalize, option, model, mesh, open_geo!

const CURRENT = Ref{Union{Nothing,GeoModel}}(nothing)
const DEFAULT_OPTIONS = Dict{String,Float64}(
    "Mesh.MeshSizeMin"=>0.0,
    "Mesh.MeshSizeMax"=>1.0e22,
    "Mesh.MeshSizeFactor"=>1.0)
const OPTIONS = copy(DEFAULT_OPTIONS)
const LAST_MESH = Ref{Union{Nothing,Mesh}}(nothing)
const LAST_MESH_LOCATOR = Ref{Union{Nothing,SimplexLocator}}(nothing)
const LAST_MESH_EDGES = Ref{Union{Nothing,MeshEdgeTopology}}(nothing)
const LAST_MESH_FACES = Ref{Union{Nothing,MeshFaceTopology}}(nothing)
const STATE_LOCK = ReentrantLock()

function _replace_mesh_cache_locked!(mesh::Union{Nothing,Mesh})
    LAST_MESH[]=mesh
    LAST_MESH_LOCATOR[]=nothing
    LAST_MESH_EDGES[]=nothing
    LAST_MESH_FACES[]=nothing
    return mesh
end

"""
    initialize()

Start a fresh process-global API session, discarding any prior model or cached
mesh and restoring supported options to their defaults.
"""
function initialize()
    lock(STATE_LOCK) do
        CURRENT[]=GeoModel()
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
    CURRENT[]===nothing && throw(ArgumentError("API: call initialize() first"))
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

_get_boundary(dim_tags,combined=true,oriented=false,recursive=false)=
    _with_model() do current
        model_boundary(
            current,dim_tags,combined,oriented,recursive)
    end

_get_adjacencies(dim,tag)=_with_model() do current
    model_adjacencies(current,dim,tag)
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
using ..API: _get_entities, _get_dimension, _get_boundary, _get_adjacencies
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
        _replace_mesh_cache_locked!(_copy_mesh(generated))
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

function _mesh_query_type_block(mesh::Mesh,element_type,tag,
                                caller::AbstractString)
    msh=_mesh_query_integer(element_type,caller,"element_type")
    msh_spec(msh)
    entity=_mesh_query_integer(tag,caller,"tag")
    entity<0 || throw(ArgumentError(
        "$caller: entity-specific data require mesh classification metadata; " *
        "use a negative tag to query the complete cache"))
    return msh,_mesh_element_block(mesh,msh)
end

function _mesh_query_tasks(task,num_tasks,caller::AbstractString)
    task_index=_mesh_query_integer(task,caller,"task")
    task_count=_mesh_query_integer(num_tasks,caller,"num_tasks")
    task_index>=0 || throw(ArgumentError(
        "$caller: task must be nonnegative"))
    task_count>=1 || throw(ArgumentError(
        "$caller: num_tasks must be positive"))
    (task_index==0 && task_count==1) || throw(ArgumentError(
        "$caller: only task=0 with num_tasks=1 is supported by the " *
        "detached Julia return arrays"))
    return nothing
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
        _mesh_query_tasks(task,num_tasks,caller)
        mesh_element_qualities(cached,element_tags,quality_name)
    end
end

function _get_jacobians(element_type,local_coord,tag=-1,task=0,num_tasks=1)
    caller="API.mesh.get_jacobians"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,_=_mesh_query_type_block(cached,element_type,tag,caller)
        _mesh_query_tasks(task,num_tasks,caller)
        mesh_jacobians(cached,msh,local_coord)
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

function _mesh_global_topology_selection(dim_tags,caller::AbstractString)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    isempty(dim_tags) || throw(ArgumentError(
        "$caller: entity-selective topology creation requires mesh " *
        "classification metadata; pass an empty collection for the complete cache"))
    return nothing
end

function _create_edges(dim_tags=())
    caller="API.mesh.create_edges"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _mesh_global_topology_selection(dim_tags,caller)
        replacement=_mesh_edge_topology(cached,LAST_MESH_EDGES[])
        LAST_MESH_EDGES[]=replacement
        nothing
    end
end

function _create_faces(dim_tags=())
    caller="API.mesh.create_faces"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _mesh_global_topology_selection(dim_tags,caller)
        replacement=_mesh_face_topology(cached,LAST_MESH_FACES[])
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

function _mesh_barycenters(mesh::Mesh,cells::Matrix{Int32},fast::Bool,
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

function _mesh_pattern_nodes(cells::Matrix{Int32},patterns)
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

function _get_nodes(dim=-1,tag=-1,include_boundary=false,
                    return_parametric_coord=true)
    caller="API.mesh.get_nodes"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        _mesh_query_bool(include_boundary,caller,"include_boundary")
        _mesh_query_bool(
            return_parametric_coord,caller,"return_parametric_coord")
        (dimension<0 && entity<0) || throw(ArgumentError(
            "$caller: entity- or dimension-specific nodes require mesh " *
            "classification metadata; use dim=-1 and a negative tag for all " *
            "cached nodes"))
        return UInt64.(1:nnodes(cached)),vec(copy(cached.coords)),Float64[]
    end
end

function _get_elements(dim=-1,tag=-1)
    caller="API.mesh.get_elements"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        entity<0 || throw(ArgumentError(
            "$caller: entity-specific elements require mesh classification " *
            "metadata; use a negative tag to query a complete dimension"))
        _mesh_element_data(cached,dimension)
    end
end

function _get_element_types(dim=-1,tag=-1)
    caller="API.mesh.get_element_types"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        entity<0 || throw(ArgumentError(
            "$caller: entity-specific element types require mesh " *
            "classification metadata; use a negative tag to query a complete " *
            "dimension"))
        _mesh_element_types(cached,dimension)
    end
end

function _get_elements_by_type(element_type,tag=-1,task=0,num_tasks=1)
    caller="API.mesh.get_elements_by_type"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block=_mesh_query_type_block(cached,element_type,tag,caller)
        _mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return UInt64[],UInt64[]
        offset,cells=block
        return _mesh_dense_tags(offset,size(cells,2)),UInt64.(vec(cells))
    end
end

function _get_nodes_by_element_type(element_type,tag=-1,
                                    return_parametric_coord=true)
    caller="API.mesh.get_nodes_by_element_type"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _,block=_mesh_query_type_block(cached,element_type,tag,caller)
        _mesh_query_bool(
            return_parametric_coord,caller,"return_parametric_coord")
        block===nothing && return UInt64[],Float64[],Float64[]
        _,cells=block
        node_tags,coordinates=_mesh_nodes_for_cells(cached,cells)
        return node_tags,coordinates,Float64[]
    end
end

function _get_barycenters(element_type,tag,fast,primary,task=0,num_tasks=1)
    caller="API.mesh.get_barycenters"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        _,block=_mesh_query_type_block(cached,element_type,tag,caller)
        fast_mode=_mesh_query_bool(fast,caller,"fast")
        _mesh_query_bool(primary,caller,"primary")
        _mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return Float64[]
        _,cells=block
        _mesh_barycenters(cached,cells,fast_mode,caller)
    end
end

function _get_element_edge_nodes(element_type,tag=-1,primary=false,
                                 task=0,num_tasks=1)
    caller="API.mesh.get_element_edge_nodes"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block=_mesh_query_type_block(cached,element_type,tag,caller)
        _mesh_query_bool(primary,caller,"primary")
        _mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return UInt64[]
        _,cells=block
        _mesh_pattern_nodes(cells,_simplex_edge_patterns(msh))
    end
end

function _get_element_face_nodes(element_type,face_type,tag=-1,primary=false,
                                 task=0,num_tasks=1)
    caller="API.mesh.get_element_face_nodes"
    return lock(STATE_LOCK) do
        cached=_cached_mesh_locked(caller)
        msh,block=_mesh_query_type_block(cached,element_type,tag,caller)
        face=_mesh_query_integer(face_type,caller,"face_type")
        face in (3,4) || throw(ArgumentError(
            "$caller: face_type must be 3 (triangle) or 4 (quadrangle)"))
        _mesh_query_bool(primary,caller,"primary")
        _mesh_query_tasks(task,num_tasks,caller)
        block===nothing && return UInt64[]
        _,cells=block
        _mesh_pattern_nodes(cells,_simplex_face_patterns(msh,face))
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
        _model_locked()
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "API.mesh.refine: no mesh; call API.mesh.generate first"))
        refined=refine_uniform(
            cached;max_nodes=max_nodes,max_cells=max_cells)
        _replace_mesh_cache_locked!(_copy_mesh(refined))
        refined
    end
end

function _clear_mesh(dim_tags=())
    return lock(STATE_LOCK) do
        _model_locked()
        (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(
            ArgumentError(
                "API.mesh.clear: dim_tags must be a vector or tuple of " *
                "(dimension, tag) pairs"))
        isempty(dim_tags) || throw(ArgumentError(
            "API.mesh.clear: entity-selective clearing requires mesh " *
            "classification metadata; pass an empty collection to clear the " *
            "complete cached mesh"))
        _replace_mesh_cache_locked!(nothing)
        nothing
    end
end

function _affine_transform_mesh(affine,dim_tags=())
    caller="API.mesh.affine_transform"
    return lock(STATE_LOCK) do
        _model_locked()
        (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(
            ArgumentError(
                "$caller: dim_tags must be a vector or tuple of " *
                "(dimension, tag) pairs"))
        isempty(dim_tags) || throw(ArgumentError(
            "$caller: entity-selective transformation requires mesh " *
            "classification metadata; pass an empty collection to transform " *
            "the complete cached mesh"))
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "$caller: no mesh; call API.mesh.generate first"))
        coefficients,translation,_=_transform_gmsh_affine(affine,caller)
        matrix=reshape(collect(coefficients),3,3)
        transformed=affine_transform(
            cached,matrix;translation=translation)
        _replace_mesh_cache_locked!(_copy_mesh(transformed))
        transformed
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
using ..API: _generate,_get_mesh,_get_nodes,_get_elements,_get_element_types,
             _get_elements_by_type,_get_nodes_by_element_type,_get_barycenters,
             _get_element_edge_nodes,_get_element_face_nodes,
             _get_element_type,_get_element_properties,
             _get_element_by_coordinates,_get_elements_by_coordinates,
             _get_local_coordinates_in_element,_get_element_qualities,
             _get_jacobians,_get_jacobian,
             _create_edges,_create_faces,_get_edges,_get_faces,
             _get_all_edges,_get_all_faces,_add_edges,_add_faces,
             _get_max_node_tag,_get_max_element_tag,
             _refine,_clear_mesh,_affine_transform_mesh,_set_size,_set_periodic,
             _get_periodic_nodes
generate(dim::Integer)=_generate(dim)
get()=_get_mesh()

"""
    get_nodes(dim=-1, tag=-1, include_boundary=false,
              return_parametric_coord=true)

Return detached dense `UInt64` node tags, flattened `Float64` coordinates, and an
empty parametric-coordinate vector for the complete cached mesh. The simplex cache
does not own entity classification or parametric coordinates, so only negative
`tag` with `dim=-1` is supported; the Boolean flags are validated but do not
change a global result.
"""
get_nodes(dim=-1,tag=-1,include_boundary=false,return_parametric_coord=true)=
    _get_nodes(dim,tag,include_boundary,return_parametric_coord)

"""
    get_elements(dim=-1, tag=-1)

Return detached Gmsh-shaped `(element_types, element_tags, node_tags)` arrays for
linear segments (type 1), triangles (type 2), and tetrahedra (type 4) in the
cached mesh. A dimension in `0:3` with negative `tag` filters whole type blocks.
Element tags are dense identifiers derived for the current cache across blocks in
that order. Global queries use `dim=-1`; values outside `-1` or `0:3` are
rejected. Entity-specific queries require unavailable classification metadata.
"""
get_elements(dim=-1,tag=-1)=_get_elements(dim,tag)

"""Return the detached linear-simplex MSH types present in a whole cache/dimension."""
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
Nondefault task partitioning is unavailable for detached Julia return arrays.
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
types 1, 2, and 4. Entity filtering and nondefault task partitioning require
metadata or caller-owned output storage that this Julia cache does not have.
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
    create_edges(dim_tags=())

Create deterministic global identifiers for every unique edge in the cached
linear-simplex mesh. Repeated calls are idempotent. Only whole-cache creation is
available because the cache does not own model-entity classification metadata.
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
fixed-node type returns empty arrays. Entity filtering and nondefault task
partitioning are explicit blockers.
"""
get_elements_by_type(element_type,tag=-1,task=0,num_tasks=1)=
    _get_elements_by_type(element_type,tag,task,num_tasks)

"""
    get_nodes_by_element_type(element_type, tag=-1,
                              return_parametric_coord=true)

Return detached node tags and coordinates in per-element connectivity order for
one fixed-node Gmsh element type. Shared nodes consequently appear once per element
use. The linear-simplex cache has no parametric or entity-classification metadata,
so the parametric result is empty and `tag` must be negative.
"""
get_nodes_by_element_type(element_type,tag=-1,return_parametric_coord=true)=
    _get_nodes_by_element_type(element_type,tag,return_parametric_coord)

"""
    get_barycenters(element_type, tag, fast, primary,
                    task=0, num_tasks=1)

Return detached `x,y,z` barycenters in element order for a cached linear-simplex
type. With `fast=true`, return unnormalized primary-node coordinate sums. All nodes
are primary for types 1, 2, and 4. Entity filtering and nondefault task partitioning
are explicit blockers.
"""
get_barycenters(element_type,tag,fast,primary,task=0,num_tasks=1)=
    _get_barycenters(element_type,tag,fast,primary,task,num_tasks)

"""
    get_element_edge_nodes(element_type, tag=-1, primary=false,
                           task=0, num_tasks=1)

Return detached edge-node tags in Gmsh local-edge order for every cached element of
one type. The `primary` flag is validated but does not change linear-simplex output.
Entity filtering and nondefault task partitioning are explicit blockers.
"""
get_element_edge_nodes(element_type,tag=-1,primary=false,task=0,num_tasks=1)=
    _get_element_edge_nodes(element_type,tag,primary,task,num_tasks)

"""
    get_element_face_nodes(element_type, face_type, tag=-1, primary=false,
                           task=0, num_tasks=1)

Return detached face-node tags in Gmsh local-face order for every cached element of
one type. `face_type` is 3 for triangles or 4 for quadrangles. The `primary` flag is
validated but does not change linear-simplex output. Entity filtering and nondefault
task partitioning are explicit blockers.
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

Clear the complete cached mesh without changing model geometry. An empty vector
or tuple selects the complete cache. Entity-selective clearing is unavailable
until the simplex cache owns entity-classification metadata.
"""
clear(dim_tags=())=_clear_mesh(dim_tags)

"""
    affine_transform(affine, dim_tags=()) -> Mesh

Apply a finite nonsingular affine transform to every node in the complete cached
mesh. `affine` is a 4×4 matrix or exactly 12 or 16 entries in Gmsh row-major
order; 12 entries imply the homogeneous row `(0, 0, 0, 1)`. The cache changes
only after the transformed mesh validates, and the returned mesh owns independent
storage. Entity-selective transforms are unavailable until the simplex cache owns
classification metadata. Model geometry and periodic relations are unchanged.
"""
affine_transform(affine,dim_tags=())=
    _affine_transform_mesh(affine,dim_tags)

"""
    set_size(dim_tags, size)

Set a finite, positive mesh-size constraint on the dimension-0 Point entities in
`dim_tags`. Every pair is validated before the active model or cached mesh changes.
"""
set_size(dim_tags,size)=_set_size(dim_tags,size)

"""
    set_periodic(dim, slave_entities, master_entities, affine; atol=1e-12)

Store validated straight-curve (`dim=1`) or planar-surface (`dim=2`) relations in
the active model and invalidate any cached mesh. `affine` maps each master entity
to its corresponding slave in Gmsh row-major 4×4 order. Each slave has one master;
masters may be reused, and a slave may become a master in an acyclic dependency
chain. Curves must share a planar surface when meshed. Surfaces must be
affine-equivalent boundaries of one explicit planar-shell volume.
"""
set_periodic(dim,slave_entities,master_entities,affine;atol=1e-12)=
    _set_periodic(dim,slave_entities,master_entities,affine;atol=atol)

"""
    get_periodic_nodes(dim, slave_entity)

Return the master entity, detached slave/master node arrays, and affine transform
for one curve or planar boundary-surface relation in the cached mesh.
"""
get_periodic_nodes(dim,slave_entity)=
    _get_periodic_nodes(dim,slave_entity)
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
