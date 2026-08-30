"""
    API

Gmsh-style model/mesh/option façade over Tessella's native kernels, including
deterministic entity-topology, analytical spatial, native type/partition metadata,
and Physical-group queries, entity-name, atomic-tag, and dependency-safe removal
mutations, Physical-group mutations, point-local mesh-size constraints, explicit
planar surface-loop volumes, operation-time Boolean operand ownership, and persistent
affine relations between straight periodic boundary or embedded curves and planar
periodic volume boundaries.
Production meshing is never delegated to Gmsh.
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
               model_is_inside, model_closest_point
using ..Model: set_entity_name!, remove_entity_name!, model_entity_name, model_set_tag!
using ..Model: remove_entities!
using ..Model: set_periodic!, model_periodic_nodes
using ..Model: mesh_model_surface, mesh_model_volume
using ..MeshTypes: Mesh
using ..GeoExec: execute_geo

export initialize, finalize, option, model, mesh, open_geo!

const CURRENT = Ref{Union{Nothing,GeoModel}}(nothing)
const DEFAULT_OPTIONS = Dict{String,Float64}(
    "Mesh.MeshSizeMin"=>0.0,
    "Mesh.MeshSizeMax"=>1.0e22,
    "Mesh.MeshSizeFactor"=>1.0)
const OPTIONS = copy(DEFAULT_OPTIONS)
const LAST_MESH = Ref{Union{Nothing,Mesh}}(nothing)
const STATE_LOCK = ReentrantLock()

"""
    initialize()

Start a fresh process-global API session, discarding any prior model or cached
mesh and restoring supported options to their defaults.
"""
function initialize()
    lock(STATE_LOCK) do
        CURRENT[]=GeoModel()
        LAST_MESH[]=nothing
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
    end
    return nothing
end

"""End the API session, discard its model and cached mesh, and restore option defaults."""
function finalize()
    lock(STATE_LOCK) do
        CURRENT[]=nothing
        LAST_MESH[]=nothing
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
        invalidate && (LAST_MESH[]=nothing)
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
        after==before || (LAST_MESH[]=nothing)
        nothing
    end
end

function _remove_entity_name(name)
    name isa AbstractString || throw(ArgumentError(
        "API.model.remove_entity_name: name must be a string"))
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_entity_name!(current,name)>0 && (LAST_MESH[]=nothing)
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
        remove_entities!(current,dim_tags,recursive)>0 && (LAST_MESH[]=nothing)
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
        after==before || (LAST_MESH[]=nothing)
        nothing
    end
end

function _remove_physical_name(name)
    name isa AbstractString || throw(ArgumentError(
        "API.model.remove_physical_name: name must be a string"))
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_physical_name!(current,name)>0 && (LAST_MESH[]=nothing)
        nothing
    end
end

function _remove_physical_groups(dim_tags=())
    return lock(STATE_LOCK) do
        current=_model_locked()
        remove_physical_groups!(current,dim_tags)>0 && (LAST_MESH[]=nothing)
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
including all live model references and its name.
"""
set_tag(dim,tag,new_tag)=_set_tag(dim,tag,new_tag)

"""
    remove_entities(dim_tags, recursive=false)

Remove ordered positive `(dimension, entity_tag)` pairs when they are not used by a
surviving boundary or embedding owner. With `recursive=true`, also process explicit
boundary entities down to Points. Missing and still-owned entities are unchanged.
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
        LAST_MESH[]=_copy_mesh(generated)
        generated
    end
end

function _get_mesh()
    return lock(STATE_LOCK) do
        _model_locked()
        LAST_MESH[]===nothing && throw(ArgumentError("API.mesh.get: no mesh"))
        _copy_mesh(LAST_MESH[])
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

"""Gmsh-style mesh generation, retrieval, and periodic curve/surface operations."""
module mesh
using ..API: _generate,_get_mesh,_set_size,_set_periodic,_get_periodic_nodes
generate(dim::Integer)=_generate(dim)
get()=_get_mesh()

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
        LAST_MESH[]=stored_mesh
        result
    end
end

end # module
