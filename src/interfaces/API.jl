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

using ..Recombine: recombine_triangles
using ..Elements: ElementBlock, msh_family, msh_dimension, msh_spec
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
using ..Model: set_periodic!, model_periodic_nodes, _model_affine_point,
              _model_cross3
using ..Transform: _transform_homogeneous
using ..Model: mesh_model_surface, mesh_model_volume, model_to_mixed,
              _model_planar_surface_mesh
using ..Model: _model_entity_known, _model_fresh_element_tag,
              _record_append_node!, _record_append_element!,
              set_transfinite_curve!,
              set_transfinite_surface!, set_transfinite_volume!,
              set_transfinite_automatic!, set_recombine!, set_smoothing!,
              set_reverse!, set_algorithm!, set_size_at_parametric_points!,
              set_size_from_boundary!, set_size_callback!, set_compound!,
              set_outward_orientation!, set_order!, remove_constraints!,
              add_discrete_entity!, add_discrete_nodes!, add_discrete_elements!,
              add_homology_request!, clear_homology_requests!,
              model_discrete_entity, model_discrete_entities, DiscreteEntity,
              classify_surfaces!, compute_homology!, create_geometry!,
              create_topology!
using ..MeshTypes: Mesh, nnodes, nsegs, ntris, ntets, validate,
                   boundary_edges
using ..IO: read_stl, GeoParams, GeoFieldSpec, _geo_split_list
using ..SizeField: AbstractSizeField, build_geo_size_field, PostViewField
using ..MeshEntityTopology: MeshEdgeTopology, MeshFaceTopology,
                            _mesh_edge_topology, _mesh_face_topology,
                            _mesh_edge_topology_for_cells,
                            _mesh_face_topology_for_cells,
                            _mesh_add_edges, _mesh_add_faces,
                            _mesh_edges, _mesh_faces, _triangle_key,
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
using ..HighOrder: HighOrder, P2Mesh, P2TriMesh, p2_trimesh, p2_tetmesh
using ..Model: _compound_surface_meshes
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
    entities::Vector{Tuple{Int,Int32}}
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
const MODEL_NAME = Ref{String}("")
const MODEL_FILE_NAME = Ref{String}("")
const ELEMENT_VISIBILITY = Ref{Dict{Int,Int32}}(Dict{Int,Int32}())
const STATE_LOCK = ReentrantLock()

# Session size fields (`model.mesh.field`): each entry owns its kind and its
# Gmsh `Field[tag].Option` values as strings, in first-set order — the same
# representation `GeoFieldSpec` carries, so generation can build them through
# `build_geo_size_field` with no format translation.
mutable struct _ApiField
    kind::String
    options::Dict{String,String}
    order::Vector{String}
end

const FIELD_SPECS = Ref{Dict{Int,_ApiField}}(Dict{Int,_ApiField}())
const FIELD_TAG_MAX = Ref(0)
const BACKGROUND_FIELD = Ref(0)
const BOUNDARY_LAYER_FIELDS = Ref{Vector{Int}}(Int[])

# Partition state attached to a cached mesh: `element_partitions[tag]` is the
# 1-based partition of dense element `tag`. The `mesh` identity anchors the
# record to the cache it describes, mirroring `_MeshClassification`.
struct _MeshPartition
    mesh::Mesh
    num_partitions::Int
    element_partitions::Vector{Int32}
end

const LAST_MESH_PARTITION = Ref{Union{Nothing,_MeshPartition}}(nothing)
const LAST_MESH_HIGH_ORDER = Ref{Union{Nothing,P2Mesh,P2TriMesh}}(nothing)
const LAST_MESH_HIGH_ORDER_MIDS = Ref{Vector{Tuple{Int,Int32}}}(Tuple{Int,Int32}[])
# The exact cache `Mesh` the overlay was elevated from — queries only trust
# the overlay while it describes the live cache object.
const LAST_MESH_HIGH_ORDER_MESH = Ref{Union{Nothing,Mesh}}(nothing)

# Each `_ModelSlot` owns one model in the session model list: its geometry, its
# mesh cache and derived caches, and its element visibility state. The global
# Refs above always mirror the *current* slot; `model.add`/`set_current`/
# `remove` capture and restore them so switching models preserves each model's
# mesh, matching Gmsh's multi-model list semantics (duplicate names allowed).
mutable struct _ModelSlot
    name::String
    file_name::String
    model::GeoModel
    mesh::Union{Nothing,Mesh}
    class::Union{Nothing,_MeshClassification}
    locator::Union{Nothing,SimplexLocator}
    edges::Union{Nothing,MeshEdgeTopology}
    faces::Union{Nothing,MeshFaceTopology}
    visibility::Dict{Int,Int32}
    window_visibility::Dict{Int,Int32}
    window_element_visibility::Dict{Int,Dict{Int,Int32}}
    partition::Union{Nothing,_MeshPartition}
    high_order::Union{Nothing,P2Mesh,P2TriMesh}
    high_order_mids::Vector{Tuple{Int,Int32}}
    fields::Dict{Int,_ApiField}
    background_field::Int
    boundary_layer_fields::Vector{Int}
end

const MODEL_SLOTS = _ModelSlot[]

function _new_slot(name::AbstractString)
    return _ModelSlot(String(name),"",GeoModel(),nothing,nothing,nothing,
                      nothing,nothing,Dict{Int,Int32}(),Dict{Int,Int32}(),
                      Dict{Int,Dict{Int,Int32}}(),nothing,nothing,
                      Tuple{Int,Int32}[],
                      Dict{Int,_ApiField}(),0,Int[])
end

function _current_slot_index_locked()
    model=CURRENT[]
    model===nothing && return 0
    index=findfirst(slot->slot.model===model,MODEL_SLOTS)
    index===nothing && return 0
    return index
end

# Copy the live globals back into the current slot before a switch.
function _slot_capture_locked!()
    index=_current_slot_index_locked()
    index==0 && return nothing
    slot=MODEL_SLOTS[index]
    slot.name=MODEL_NAME[]
    slot.file_name=MODEL_FILE_NAME[]
    slot.model=CURRENT[]
    slot.mesh=LAST_MESH[]
    slot.class=LAST_MESH_CLASS[]
    slot.locator=LAST_MESH_LOCATOR[]
    slot.edges=LAST_MESH_EDGES[]
    slot.faces=LAST_MESH_FACES[]
    slot.partition=LAST_MESH_PARTITION[]
    slot.high_order=LAST_MESH_HIGH_ORDER[]
    slot.high_order_mids=LAST_MESH_HIGH_ORDER_MIDS[]
    slot.fields=FIELD_SPECS[]
    slot.background_field=BACKGROUND_FIELD[]
    slot.boundary_layer_fields=BOUNDARY_LAYER_FIELDS[]
    # `visibility` is a shared Dict object; its contents need no copy.
    return nothing
end

function _slot_load_locked!(slot::_ModelSlot)
    CURRENT[]=slot.model
    MODEL_NAME[]=slot.name
    MODEL_FILE_NAME[]=slot.file_name
    LAST_MESH[]=slot.mesh
    LAST_MESH_CLASS[]=slot.class
    LAST_MESH_LOCATOR[]=slot.locator
    LAST_MESH_EDGES[]=slot.edges
    LAST_MESH_FACES[]=slot.faces
    LAST_MESH_PARTITION[]=slot.partition
    LAST_MESH_HIGH_ORDER[]=slot.high_order
    LAST_MESH_HIGH_ORDER_MESH[]=slot.high_order===nothing ? nothing : slot.mesh
    LAST_MESH_HIGH_ORDER_MIDS[]=slot.high_order_mids
    ELEMENT_VISIBILITY[]=slot.visibility
    FIELD_SPECS[]=slot.fields
    BACKGROUND_FIELD[]=slot.background_field
    BOUNDARY_LAYER_FIELDS[]=slot.boundary_layer_fields
    return nothing
end

function _slot_clear_globals_locked!()
    CURRENT[]=nothing
    MODEL_NAME[]=""
    MODEL_FILE_NAME[]=""
    LAST_MESH[]=nothing
    LAST_MESH_CLASS[]=nothing
    LAST_MESH_LOCATOR[]=nothing
    LAST_MESH_EDGES[]=nothing
    LAST_MESH_FACES[]=nothing
    LAST_MESH_PARTITION[]=nothing
    LAST_MESH_HIGH_ORDER[]=nothing
    LAST_MESH_HIGH_ORDER_MESH[]=nothing
    LAST_MESH_HIGH_ORDER_MIDS[]=Tuple{Int,Int32}[]
    ELEMENT_VISIBILITY[]=Dict{Int,Int32}()
    FIELD_SPECS[]=Dict{Int,_ApiField}()
    BACKGROUND_FIELD[]=0
    BOUNDARY_LAYER_FIELDS[]=Int[]
    return nothing
end

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
    LAST_MESH_PARTITION[]=nothing
    LAST_MESH_HIGH_ORDER[]=nothing
    LAST_MESH_HIGH_ORDER_MESH[]=nothing
    LAST_MESH_HIGH_ORDER_MIDS[]=Tuple{Int,Int32}[]
    empty!(ELEMENT_VISIBILITY[])
    index=_current_slot_index_locked()
    index!=0 && empty!(MODEL_SLOTS[index].window_element_visibility)
    return mesh
end

# Re-elevate a mutated cache to order 2 — used by every mesh mutation that
# preserves Gmsh's "order survives mutation" semantics. Midnode tags are
# re-densified along the new skeleton's edges, consistent with this session's
# positional node-tag model.
function _rebind_high_order!(cache::Mesh,class,caller::AbstractString)
    dimension=ntets(cache)>0 ? 3 : ntris(cache)>0 ? 2 : 0
    dimension==0 && return nothing
    p2=try
        dimension==2 ? p2_trimesh(cache) :
            p2_tetmesh(cache;require_positive_tets=false)
    catch err
        err isa InterruptException && rethrow()
        throw(ErrorException(
            "$caller: order-2 re-elevation failed — "*sprint(showerror,err)))
    end
    entities=class===nothing ? Tuple{Int,Int32}[(dimension,Int32(0))] :
        class.entities
    _set_high_order_overlay(p2,cache,_p2_mid_owners(
        p2,cache,class,dimension,Int.(last.(entities))))
    return nothing
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
        empty!(MODEL_SLOTS)
        push!(MODEL_SLOTS,_new_slot(""))
        _slot_load_locked!(MODEL_SLOTS[1])
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
        FIELD_TAG_MAX[]=0
        empty!(SESSION_VIEWS);VIEW_TAG_MAX[]=0
    end
    return nothing
end

"""End the API session, discard its model and cached mesh, and restore option defaults."""
function finalize()
    lock(STATE_LOCK) do
        empty!(MODEL_SLOTS)
        _slot_clear_globals_locked!()
        empty!(OPTIONS);merge!(OPTIONS,DEFAULT_OPTIONS)
        FIELD_TAG_MAX[]=0
        empty!(SESSION_VIEWS);VIEW_TAG_MAX[]=0
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

# Multi-model session state: the session owns an ordered list of model slots
# (duplicate names are permitted, as in Gmsh). `add` always appends a fresh
# model and selects it; `set_current` selects the first slot with a matching
# name; `remove` deletes the current slot and selects the last remaining one.
# `get_file_name`/`set_file_name` track the current model's associated file.
function _model_add(name)
    caller="API.model.add"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    return lock(STATE_LOCK) do
        _slot_capture_locked!()
        slot=_new_slot(name)
        push!(MODEL_SLOTS,slot)
        _slot_load_locked!(slot)
        return nothing
    end
end

function _model_remove()
    return lock(STATE_LOCK) do
        _model_locked()
        index=_current_slot_index_locked()
        index!=0 || throw(ArgumentError(
            "API.model.remove: the current model is not in the model list"))
        deleteat!(MODEL_SLOTS,index)
        if isempty(MODEL_SLOTS)
            _slot_clear_globals_locked!()
        else
            _slot_load_locked!(MODEL_SLOTS[end])
        end
        return nothing
    end
end

function _model_list()
    return lock(STATE_LOCK) do
        return String[slot.name for slot in MODEL_SLOTS]
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
        target=String(name)
        index=findfirst(slot->slot.name==target,MODEL_SLOTS)
        index===nothing && throw(ArgumentError(
            "$caller: unknown model \"$target\""))
        _slot_capture_locked!()
        _slot_load_locked!(MODEL_SLOTS[index])
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

# `model.set_visibility_per_window` is the global per-window flag — display
# state only, stored on the current slot so it survives model switches.
function _set_visibility_per_window(value,window_index)
    caller="API.model.set_visibility_per_window"
    flag=_mesh_query_integer(value,caller,"value")
    typemin(Int32)<=flag<=typemax(Int32) || throw(ArgumentError(
        "$caller: value $value is out of range"))
    window=_mesh_query_integer(window_index,caller,"window_index")
    window>=0 || throw(ArgumentError(
        "$caller: window_index must be non-negative"))
    return lock(STATE_LOCK) do
        _model_locked()
        index=_current_slot_index_locked()
        index==0 && throw(ArgumentError(
            "$caller: the current model is not in the model list"))
        MODEL_SLOTS[index].window_visibility[window]=Int32(flag)
        return nothing
    end
end

function _add_discrete_entity(dim,tag=0,boundary=())
    caller="API.model.add_discrete_entity"
    dimension=_mesh_query_integer(dim,caller,"dim")
    (boundary isa AbstractVector || boundary isa Tuple) || throw(ArgumentError(
        "$caller: boundary must be a vector or tuple of entity tags or " *
        "(dimension, tag) pairs"))
    pairs=map(boundary) do entry
        if entry isa Integer
            # Gmsh parity: bare tags refer to entities of dimension dim-1.
            return (dimension-1,Int(entry))
        elseif (entry isa Tuple || entry isa AbstractVector ||
                entry isa Pair) && length(entry)==2
            return (Int(first(entry)),Int(last(entry)))
        end
        throw(ArgumentError(
            "$caller: each boundary entry must be an entity tag or a " *
            "(dimension, tag) pair"))
    end
    return lock(STATE_LOCK) do
        add_discrete_entity!(_model_locked(),dim,tag,pairs)
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
    if haskey(model.discrete,key)
        return Int32[abs(pair[2]) for pair in
            model.discrete[key].boundary if pair[1]==key[1]-1]
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
        known=_model_entity_known(model,dimension,entity) ||
              (class!==nothing &&
               haskey(class.boundaries,(dimension,Int32(entity))))
        known || throw(ArgumentError(
            "$caller: unknown $(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
        highest=model_dimension(model)
        connected=Set{Tuple{Int,Int32}}()
        queue=Tuple{Int,Int32}[(highest,Int32(top)) for top in
            keys(_mesh_entity_dictionary(model,highest))]
        for (edim,etag) in keys(model.discrete)
            edim==highest && push!(queue,(edim,Int32(etag)))
        end
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
    partition=LAST_MESH_PARTITION[]
    (partition===nothing || partition.mesh!==LAST_MESH[]) && return 0
    return partition.num_partitions
end

_get_partitions(dim,tag)=_with_model() do current
    model_partitions(current,dim,tag)
    partition=LAST_MESH_PARTITION[]
    (partition===nothing || partition.mesh!==LAST_MESH[]) && return Int[]
    class=LAST_MESH_CLASS[]
    (class===nothing || class.mesh!==partition.mesh) && return Int[]
    cached=partition.mesh
    entity=_mesh_query_integer(tag,"API.model.get_partitions","tag")
    dimension=_mesh_query_integer(dim,"API.model.get_partitions","dim")
    result=Int[]
    for element_tag in _mesh_entity_element_tags(
            class,current,dimension,entity,cached)
        p=Int(partition.element_partitions[element_tag])
        p in result || push!(result,p)
    end
    return sort!(result)
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
             _model_add, _model_remove, _model_list, _set_visibility_per_window,
             _add_discrete_entity
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

Append a new empty model named `name` to the session's model list and make it
current, matching Gmsh's `model.add`. The previously current model keeps its
geometry, mesh, and visibility state; duplicate names are permitted.
"""
add(name)=_model_add(name)

"""
Remove the session's current model from the model list and make the last
remaining model current. Removing the final model leaves the session
model-less until `add` is called.
"""
remove()=_model_remove()

"""Return the session's model names in insertion order."""
list()=_model_list()

"""
    set_current(name)

Make the first model named `name` current. Its geometry, mesh cache, and
element visibility state are restored exactly as it was left; unknown names
throw.
"""
set_current(name)=_set_current(name)

"""Return the file name associated with the current model, or `""` when unset."""
get_file_name()=_get_file_name()

"""Associate `name` with the current model as its file name."""
set_file_name(name)=_set_file_name(name)

"""
    set_visibility_per_window(value, window_index=0)

Store the global display visibility flag for window `window_index`, matching
Gmsh's `model.setVisibilityPerWindow`. This is display state only — there is no
corresponding getter.
"""
set_visibility_per_window(value,window_index=0)=
    _set_visibility_per_window(value,window_index)

"""
    add_discrete_entity(dim, tag=0, boundary=[]) -> tag

Create a parametrization-free entity of dimension `dim`, matching Gmsh's
`model.addDiscreteEntity`. `tag=0` allocates the next tag in `dim`; `boundary`
lists the declared `(dim, tag)` boundary entities used by topology queries.
Nodes and elements are attached with `mesh.add_nodes`, `mesh.add_elements`, or
`mesh.add_elements_by_type`.
"""
add_discrete_entity(dim,tag=0,boundary=[])=_add_discrete_entity(dim,tag,boundary)

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
        cache,(dim,Int32(tag)),[(dim,Int32(tag))],copy(data.node_entities),
        boundaries,seg_entities,tri_entities,tet_entities)
end

# Positively-oriented copy of a reversed volume part for classification:
# reversal swaps the first two nodes of every tet, so the same swap restores
# the orientation `model_to_mixed` certifies while keeping column order
# aligned with the merged cache.
function _classification_skeleton(mesh::Mesh,reversed::Bool)
    (reversed && ntets(mesh)>0) || return mesh
    tets=Matrix{Int32}(mesh.tets)
    @inbounds for cell in axes(tets,2)
        tets[1,cell],tets[2,cell]=tets[2,cell],tets[1,cell]
    end
    return Mesh(mesh.coords;tets=tets)
end

# Merge per-entity classified meshes into one shared-node mesh plus a merged
# `_MeshClassification`, matching Gmsh's `generate(dim)` semantics of meshing
# every remaining top-level entity into a single node-addressable mesh.
# `parts` holds `(entity_tag, part_mesh, part_class)` triples in sorted tag
# order. Nodes shared across parts (boundary nodes on common curves/points)
# unify by rounded coordinate key; on conflicting ownership the
# lower-dimension owner wins, matching `model_to_mixed`'s classification rule.
function _merge_classified_parts(parts,remaps,seg_keeps,merged::Mesh,
                                 caller::AbstractString)
    entities=Tuple{Int,Int32}[]
    node_entities=fill((0,Int32(0)),nnodes(merged))
    boundaries=Dict{Tuple{Int,Int32},Vector{Int32}}()
    seg_entities=Int32[];tri_entities=Int32[];tet_entities=Int32[]
    dim=nothing
    for (tag,_,class) in parts
        class===nothing && return nothing
        dim===nothing && (dim=class.entity[1])
        class.entity[1]==dim || return nothing
        push!(entities,(dim,Int32(tag)))
    end
    for (part,(_,part_mesh,class)) in enumerate(parts)
        remap=remaps[part]
        for (column,entity) in enumerate(class.node_entities)
            owner=_merged_node_owner(node_entities[remap[column]],entity)
            node_entities[remap[column]]=owner
        end
        for (key,tags) in class.boundaries
            merged_tags=get!(()->Int32[],boundaries,key)
            for boundary in tags
                boundary in merged_tags || push!(merged_tags,boundary)
            end
        end
        append!(seg_entities,
                length(class.seg_entities)==length(seg_keeps[part]) ?
                    class.seg_entities[seg_keeps[part]] :
                    class.seg_entities)
        append!(tri_entities,class.tri_entities)
        append!(tet_entities,class.tet_entities)
    end
    return _MeshClassification(
        merged,entities[1],entities,node_entities,boundaries,
        seg_entities,tri_entities,tet_entities)
end

@inline function _merged_node_owner(current,next)
    current[1]==0 && return next
    next[1]==0 && return current
    current[1]!=next[1] && return min(current[1],next[1])==current[1] ?
        current : next
    return min(current,next)
end

# Concatenate disjoint per-entity meshes into one shared-node mesh, unifying
# nodes whose coordinates coincide (boundary nodes on entities common to
# several parts). Columns of `segs`/`tris`/`tets` stay part-ordered so the
# merged per-cell entity vectors remain index-aligned.
function _merge_entity_meshes(parts,caller::AbstractString;
                              require_positive_tets::Bool=true)
    index=Dict{NTuple{3,Float64},Int32}()
    coordinates=NTuple{3,Float64}[]
    remaps=Vector{Vector{Int32}}(undef,length(parts))
    for (part,(tag,mesh)) in enumerate(parts)
        remap=Vector{Int32}(undef,nnodes(mesh))
        for node in 1:nnodes(mesh)
            key=ntuple(axis->round(mesh.coords[axis,node],
                                   sigdigits=12)+0.0,3)
            remap[node]=get!(index,key) do
                push!(coordinates,(mesh.coords[1,node],mesh.coords[2,node],
                                   mesh.coords[3,node]))
                Int32(length(coordinates))
            end
        end
        remaps[part]=remap
    end
    length(coordinates)<=typemax(Int32) || throw(ArgumentError(
        "$caller: merged mesh exceeds the Int32 node limit"))
    counts=zeros(Int,3)
    for (_,mesh) in parts
        counts[1]+=nsegs(mesh);counts[2]+=ntris(mesh);counts[3]+=ntets(mesh)
    end
    segs=Matrix{Int32}(undef,2,counts[1])
    tris=Matrix{Int32}(undef,3,counts[2])
    tets=Matrix{Int32}(undef,4,counts[3])
    seg_tag=Int32[];tri_tag=Int32[];tet_tag=Int32[]
    # Segments on seams shared between parts (e.g. compound members) appear in
    # every part — keep the first occurrence so each curve owns one chain.
    seen_segs=Set{NTuple{2,Int32}}()
    seg_keeps=[falses(nsegs(mesh)) for (_,mesh) in parts]
    seg_count=0
    for (part,(_,mesh)) in enumerate(parts)
        remap=remaps[part]
        for cell in axes(mesh.segs,2)
            a=remap[mesh.segs[1,cell]];b=remap[mesh.segs[2,cell]]
            key=a<b ? (a,b) : (b,a)
            key in seen_segs && continue
            push!(seen_segs,key)
            seg_keeps[part][cell]=true
            seg_count+=1
            segs[1,seg_count]=a;segs[2,seg_count]=b
            push!(seg_tag,mesh.seg_tag[cell])
        end
    end
    tri_count=0;tet_count=0
    for (part,(_,mesh)) in enumerate(parts)
        remap=remaps[part]
        for cell in axes(mesh.tris,2),slot in 1:3
            tris[slot,tri_count+cell]=remap[mesh.tris[slot,cell]]
        end
        for cell in axes(mesh.tets,2),slot in 1:4
            tets[slot,tet_count+cell]=remap[mesh.tets[slot,cell]]
        end
        append!(tri_tag,mesh.tri_tag)
        append!(tet_tag,mesh.tet_tag)
        tri_count+=ntris(mesh);tet_count+=ntets(mesh)
    end
    segs=Matrix{Int32}(segs[:,1:seg_count])
    coordinate_matrix=Matrix{Float64}(undef,3,length(coordinates))
    for (node,coordinate) in enumerate(coordinates)
        coordinate_matrix[:,node].=coordinate
    end
    merged=Mesh(coordinate_matrix;segs=segs,tris=tris,tets=tets,
                seg_tag=seg_tag,tri_tag=tri_tag,tet_tag=tet_tag)
    diag=validate(merged;require_positive_tets=require_positive_tets)
    diag.ok || throw(ErrorException(
        "$caller: merged mesh is invalid — "*join(diag.messages,"; ")))
    return merged,remaps,seg_keeps
end

# Gmsh meshing-algorithm selection consumed by generation: this build's
# kernels are the Delaunay-family surface mesher and the tetrahedral filler —
# the choices `model.mesh.set_algorithm` maps onto. Other algorithm ids are
# rejected loudly rather than silently ignored.
const _GENERATE_SURFACE_ALGORITHMS=(Int32(1),Int32(2),Int32(3))
const _GENERATE_VOLUME_ALGORITHMS=(Int32(1),)

function _validate_generate_algorithms(m::GeoModel,dimension::Int,
                                       caller::AbstractString)
    supported=dimension==2 ? _GENERATE_SURFACE_ALGORITHMS :
               _GENERATE_VOLUME_ALGORITHMS
    entities=dimension==2 ? m.surfaces : m.volumes
    for ((adim,atag),algorithm) in m.meshing.algorithm
        adim==dimension || continue
        haskey(entities,atag) || continue
        Int32(algorithm) in supported || throw(ArgumentError(
            "$caller: mesh algorithm $algorithm on " *
            "$(_MESH_ENTITY_LABELS[dimension+1])[$atag] is not implemented; " *
            "supported: $(join(Int.(supported),", "))"))
    end
    return nothing
end

_p2_cells(p::P2TriMesh)=p.tri6
_p2_cells(p::P2Mesh)=p.tet10
_p2_etype(p::P2TriMesh)=Int32(9)
_p2_etype(p::P2Mesh)=Int32(11)
_p2_skeleton_etype(p::P2TriMesh)=Int32(2)
_p2_skeleton_etype(p::P2Mesh)=Int32(4)
_p2_edge_slots(p::P2TriMesh)=HighOrder._P2TRI_EDGE_SLOTS
_p2_edge_slots(p::P2Mesh)=HighOrder._P2_EDGE_SLOTS

# Elevated elements share the linear skeleton's element-tag space — type 9
# triangles use the triangle offset, type 11 tets the tetrahedron offset.
function _p2_tag_offset(mesh::Mesh,p::P2TriMesh)
    tri_offset,_,_=_mesh_element_offsets(mesh)
    return tri_offset
end

function _p2_tag_offset(mesh::Mesh,p::P2Mesh)
    _,tet_offset,_=_mesh_element_offsets(mesh)
    return tet_offset
end

# Overlay midnode indices (tag - linear node count) owned by `(dim, tag)` —
# a preallocation-free two-pass scan over the stored owner table.
function _p2_entity_mids(dim::Int,tag::Int)
    mids=LAST_MESH_HIGH_ORDER_MIDS[]
    count=0
    @inbounds for owner in mids
        owner==(dim,Int32(tag)) && (count+=1)
    end
    result=Vector{Int}(undef,count)
    position=0
    @inbounds for (index,owner) in enumerate(mids)
        owner==(dim,Int32(tag)) && (result[position+=1]=index)
    end
    return result
end

# Append an overlay's midnodes owned by `(dim, tag)` to a node payload,
# evaluating model parametrization when `parametric` — the same source the
# linear path uses, so midnode params match Gmsh's per-node params.
function _append_p2_midnodes!(node_tags,coordinates,parameters,
                              model::GeoModel,cached::Mesh,dim::Int,tag::Int,
                              parametric::Bool)
    overlay=_high_order_overlay(cached)
    overlay===nothing && return nothing
    mids=_p2_entity_mids(dim,tag)
    isempty(mids) && return nothing
    n0=nnodes(cached)
    for mid in mids
        push!(node_tags,UInt64(n0+mid))
        @inbounds for axis in 1:3
            push!(coordinates,overlay.coords[axis,n0+mid])
        end
    end
    if parametric && dim in (1,2)
        values=Vector{Float64}(undef,Base.checked_mul(3,length(mids)))
        @inbounds for (position,mid) in enumerate(mids),axis in 1:3
            values[3position-3+axis]=overlay.coords[axis,n0+mid]
        end
        append!(parameters,model_parametrization(model,dim,tag,values))
    end
    return nothing
end

# Owner entity for each quadratic midnode, Gmsh-style: a midnode belongs to
# the lowest-dimension entity containing its skeleton edge — approximated by
# endpoint ownership (equal owners win; different-dim owners defer to the
# higher-dim endpoint since a low-dim entity owns only its corner node; equal
# dims on different entities mean an interior edge owned by the cell).
function _p2_mid_owners(p2,mesh::Mesh,class,dimension::Int,entities)
    conn=_p2_cells(p2)
    skeleton=dimension==2 ? mesh.tris : mesh.tets
    cell_entities=class===nothing ? nothing :
        dimension==2 ? class.tri_entities : class.tet_entities
    node_entities=class===nothing ? Tuple{Int,Int32}[] :
        class.node_entities
    fallback=(dimension,Int32(first(entities)))
    owners=fill(fallback,nnodes(p2)-nnodes(mesh))
    @inbounds for cell in axes(conn,2),(slot,i,j) in _p2_edge_slots(p2)
        mid=Int(conn[slot,cell])
        index=mid-nnodes(mesh)
        (index<1 || index>length(owners)) && continue
        owner=if isempty(node_entities)
            fallback
        else
            a=Int(skeleton[i,cell]);b=Int(skeleton[j,cell])
            ea=a<=length(node_entities) ? node_entities[a] : fallback
            eb=b<=length(node_entities) ? node_entities[b] : fallback
            if ea==eb
                ea
            elseif ea[1]==eb[1]
                cell_entities===nothing || cell>length(cell_entities) ?
                    fallback : (dimension,cell_entities[cell])
            else
                ea[1]>eb[1] ? ea : eb
            end
        end
        owners[index]=owner
    end
    return owners
end

# The live order-2 overlay, or `nothing` when the cached mesh is not the
# object the overlay was elevated from — guards every query against a stale
# overlay after an unmanaged cache replacement.
function _high_order_overlay(cached::Union{Nothing,Mesh})
    cached===nothing && return nothing
    overlay=LAST_MESH_HIGH_ORDER[]
    overlay===nothing && return nothing
    LAST_MESH_HIGH_ORDER_MESH[]===cached || return nothing
    return overlay
end

function _set_high_order_overlay(p2,cached::Union{Nothing,Mesh},
                                 mids=Tuple{Int,Int32}[])
    LAST_MESH_HIGH_ORDER[]=p2
    LAST_MESH_HIGH_ORDER_MESH[]=p2===nothing ? nothing : cached
    LAST_MESH_HIGH_ORDER_MIDS[]=mids
    return nothing
end

# Consume `set_order`: order 1 keeps the linear mesh; order 2 elevates the
# top-dimension element family to Gmsh type 9/11 through `p2_trimesh`/
# `p2_tetmesh`, storing the quadratic mesh as a session overlay. The linear
# skeleton remains the canonical cache for classification and topology.
function _apply_mesh_order(m::GeoModel,cached::Mesh,class,caller)
    order=m.meshing.order
    order in (1,2) || throw(ArgumentError(
        "$caller: mesh order $order is not supported (supported: 1, 2)"))
    if order==1
        _set_high_order_overlay(nothing,nothing)
        return nothing
    end
    dimension=ntets(cached)>0 ? 3 : 2
    p2=try
        dimension==2 ? p2_trimesh(cached) :
            p2_tetmesh(cached;require_positive_tets=false)
    catch err
        err isa InterruptException && rethrow()
        throw(ErrorException(
            "$caller: order-2 elevation failed — "*sprint(showerror,err)))
    end
    entities=class===nothing ? Tuple{Int,Int32}[] : class.entities
    isempty(entities) &&
        (entities=Tuple{Int,Int32}[(dimension,Int32(0))])
    _set_high_order_overlay(p2,cached,_p2_mid_owners(
        p2,cached,class,dimension,Int.(last.(entities))))
    return nothing
end

# Per-entity generate units: compounds of dimension `dim` mesh their member
# entities as one patch, everything else meshes independently. Returns a list
# of `(tags, mesh_parts)` work units in first-member-tag order, where
# `mesh_parts` is the list of `(member_tag, Mesh)` pairs the unit produces.
function _generate_units(m::GeoModel,dimension::Int,entities::Vector{Int})
    members_of=Dict{Int,Int}()
    compounds=Tuple{Int,Vector{Int}}[]
    for (compound_dimension,tags) in m.meshing.compounds
        compound_dimension==dimension || continue
        members=sort!(Int.(collect(tags)))
        for tag in members
            tag in entities || throw(ArgumentError(
                "API.mesh.generate: compound member " *
                "$(_MESH_ENTITY_LABELS[dimension+1])[$tag] is not a " *
                "remaining entity"))
            haskey(members_of,tag) && throw(ArgumentError(
                "API.mesh.generate: entity $tag appears in multiple " *
                "dimension-$dimension compounds"))
            members_of[tag]=members
        end
    end
    groups=Vector{Vector{Int}}()
    seen=Set{Int}()
    for tag in entities
        tag in seen && continue
        members=get(members_of,tag,[tag])
        union!(seen,members)
        push!(groups,members)
    end
    sort!(groups;by=first)
    return groups
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
        caller="API.mesh.generate"
        _validate_generate_algorithms(m,dimension,caller)
        size_field=_session_size_field_locked(m)
        parts=Tuple{Int,Mesh}[]
        if dimension==2
            isempty(m.surfaces) && throw(ArgumentError("$caller: no surfaces"))
            units=_generate_units(m,2,sort!(collect(keys(m.surfaces))))
            for tags in units
                if length(tags)==1
                    push!(parts,(tags[1],mesh_model_surface(
                        m,tags[1];size_field=size_field)))
                else
                    for (tag,mesh) in _compound_surface_meshes(
                            m,tags,caller;size_field=size_field)
                        push!(parts,(tag,mesh))
                    end
                end
            end
        elseif dimension==3
            isempty(m.volumes) && throw(ArgumentError("$caller: no volumes"))
            for tag in sort!(collect(keys(m.volumes)))
                push!(parts,(tag,mesh_model_volume(
                    m,tag;size_field=size_field)))
            end
        else
            throw(ArgumentError("$caller: dim must be 2 or 3"))
        end
        generated,remaps,seg_keeps=_merge_entity_meshes(
            parts,caller;require_positive_tets=!any(
                tag->get(m.meshing.reverse,(dimension,tag),false),
                (tag for (tag,_) in parts)))
        classified=[(tag,mesh,_classify_cached_mesh(
                         m,_classification_skeleton(
                             mesh,dimension==3 &&
                                 get(m.meshing.reverse,(3,tag),false)),
                         dimension,tag,mesh)) for (tag,mesh) in parts]
        cache=_copy_mesh(generated)
        class=_merge_classified_parts(classified,remaps,seg_keeps,cache,caller)
        _replace_mesh_cache_locked!(cache,class)
        _apply_mesh_order(m,cache,class,caller)
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
        haskey(_mesh_entity_dictionary(m,dim),tag) ||
        haskey(m.discrete,(dim,tag)) ||
        haskey(m.meshing.attached,(dim,tag)) || throw(ArgumentError(
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

# --- Discrete/attached mesh data -------------------------------------------
# Nodes and elements added through `add_nodes`/`add_elements` live on the model
# (the entity's own record for discrete entities, `meshing.attached` for native
# entities) with caller-assigned — possibly sparse — tags. Queries below merge
# those records with the generated mesh cache; records answer entity-scoped
# queries even when no generated mesh exists, matching Gmsh's behavior on
# discrete models.

"""The [`DiscreteEntity`](@ref) record for `(dim, tag)`, or `nothing`."""
_model_mesh_record(m::GeoModel,dim::Int,tag::Int)=
    get(m.discrete,(dim,tag),get(m.meshing.attached,(dim,tag),nothing))

"""`(dim, tag, record)` triples in sorted order for every record holding data."""
function _discrete_mesh_records(m::GeoModel)
    pairs=Tuple{Int,Int,DiscreteEntity}[]
    for (key,record) in m.discrete
        push!(pairs,(key[1],key[2],record))
    end
    for (key,record) in m.meshing.attached
        push!(pairs,(key[1],key[2],record))
    end
    return sort!(pairs;by=entry->(entry[1],entry[2]))
end

# Global node-tag → (record, position) index across every discrete and
# attached record. Record node tags are globally unique, so connectivity may
# reference nodes owned by a different (boundary) entity's record.
function _record_node_index(m::GeoModel)
    index=Dict{Int32,Tuple{DiscreteEntity,Int}}()
    for (_,_,record) in _discrete_mesh_records(m)
        for (i,tag) in enumerate(record.node_tags)
            index[tag]=(record,i)
        end
    end
    return index
end

# Append a record's nodes to running tag/coordinate/parametric arrays.
function _append_record_nodes!(tags,coordinates,parameters,record,
                               parametric::Bool)
    for i in eachindex(record.node_tags)
        push!(tags,UInt64(record.node_tags[i]))
        append!(coordinates,@view(record.node_coords[:,i]))
    end
    parametric && append!(parameters,vec(record.node_params))
    return nothing
end

# Append nodes of the transitive boundary-entity closure of `record`. Each
# boundary entity contributes once, in depth-first declaration order; nodes
# carry `query_dim` copies of -1.0 as parameters (Gmsh pads lower-dimensional
# parametric coordinates with -1 on `include_boundary` queries).
function _append_record_boundary_nodes!(tags,coordinates,parameters,
                                        m::GeoModel,record,query_dim::Int,
                                        parametric::Bool)
    seen=Set{Tuple{Int,Int}}()
    stack=reverse!(collect(record.boundary))
    while !isempty(stack)
        key=pop!(stack)
        key in seen && continue
        push!(seen,key)
        boundary_record=_model_mesh_record(m,key[1],key[2])
        boundary_record===nothing && continue
        for i in eachindex(boundary_record.node_tags)
            push!(tags,UInt64(boundary_record.node_tags[i]))
            append!(coordinates,@view(boundary_record.node_coords[:,i]))
            parametric && append!(parameters,fill(-1.0,query_dim))
        end
        for entry in Iterators.reverse(boundary_record.boundary)
            push!(stack,entry)
        end
    end
    return nothing
end

# Group a record's elements into `(element_types, element_tags, node_tags)` in
# first-seen type order, matching Gmsh's per-type emission.
function _record_element_blocks(record)
    order=Int32[]
    buckets=Dict{Int32,Tuple{Vector{UInt64},Vector{UInt64}}}()
    for i in eachindex(record.element_tags)
        element_type=record.element_types[i]
        bucket=get!(buckets,element_type,(UInt64[],UInt64[]))
        isempty(bucket[1]) && push!(order,element_type)
        push!(bucket[1],UInt64(record.element_tags[i]))
        append!(bucket[2],UInt64.(record.element_nodes[i]))
    end
    element_types=Int32[]
    element_tags=Vector{UInt64}[]
    node_tags=Vector{UInt64}[]
    for element_type in order
        push!(element_types,element_type)
        push!(element_tags,buckets[element_type][1])
        push!(node_tags,buckets[element_type][2])
    end
    return element_types,element_tags,node_tags
end

# Merge cache and record elements into per-type blocks in Gmsh's emission
# order: entities iterate in (dim, tag) order; each entity's elements join the
# global bucket of their type (first-appearance type ordering).
function _merged_element_data(m::GeoModel,cached,dimension::Int)
    order=Int32[]
    buckets=Dict{Int32,Tuple{Vector{UInt64},Vector{UInt64}}}()
    emit! = function (etype::Int32,tags,nodes)
        bucket=get!(buckets,etype) do
            push!(order,etype)
            (UInt64[],UInt64[])
        end
        append!(bucket[1],tags)
        append!(bucket[2],nodes)
        return nothing
    end
    overlay=_high_order_overlay(cached)
    class=cached===nothing ? nothing : _cached_classification_locked(cached)
    if class===nothing
        # Unclassified cache: flat type blocks first, then record blocks.
        if cached!==nothing
            skeleton_etype=overlay===nothing ? Int32(-1) :
                _p2_skeleton_etype(overlay)
            etypes,etags,enodes=_mesh_element_data(cached,dimension)
            for i in eachindex(etypes)
                etypes[i]==skeleton_etype && continue
                emit!(etypes[i],etags[i],enodes[i])
            end
            if overlay!==nothing &&
                    (dimension<0 || dimension==Int(skeleton_etype))
                conn=_p2_cells(overlay)
                offset=_p2_tag_offset(cached,overlay)
                emit!(_p2_etype(overlay),
                      UInt64.(offset .+ (1:size(conn,2))),
                      UInt64.(vec(conn)))
            end
        end
        for (edim,etag,record) in _discrete_mesh_records(m)
            (dimension<0 || edim==dimension) || continue
            rtypes,rtags,rnodes=_record_element_blocks(record)
            for i in eachindex(rtypes)
                emit!(rtypes[i],rtags[i],rnodes[i])
            end
        end
    else
        tri_offset,tet_offset,_=_mesh_element_offsets(cached)
        for dim in (dimension<0 ? (0:3) : (dimension:dimension))
            tags=Int.(_mesh_classified_entities(class,dim))
            for (edim,etag,_) in _discrete_mesh_records(m)
                edim==dim && push!(tags,etag)
            end
            sort!(unique!(tags))
            for etag in tags
                for (etype,offset,cells) in
                    ((Int32(1),0,cached.segs),
                     (Int32(2),tri_offset,cached.tris),
                     (Int32(4),tet_offset,cached.tets))
                    Int(msh_spec(etype).dim)==dim || continue
                    positions=_mesh_entity_cell_positions(
                        class,Int(etype),etag)
                    isempty(positions) && continue
                    if overlay!==nothing &&
                            _p2_skeleton_etype(overlay)==etype
                        emit!(_p2_etype(overlay),
                              UInt64.(offset .+ positions),
                              UInt64.(vec(@view _p2_cells(overlay)[
                                  :,positions])))
                    else
                        emit!(etype,UInt64.(offset .+ positions),
                              UInt64.(vec(@view cells[:,positions])))
                    end
                end
                record=_model_mesh_record(m,dim,etag)
                if record!==nothing
                    rtypes,rtags,rnodes=_record_element_blocks(record)
                    for i in eachindex(rtypes)
                        emit!(rtypes[i],rtags[i],rnodes[i])
                    end
                end
            end
        end
    end
    element_tags=Vector{UInt64}[]
    node_tags=Vector{UInt64}[]
    for etype in order
        push!(element_tags,buckets[etype][1])
        push!(node_tags,buckets[etype][2])
    end
    return order,element_tags,node_tags
end

function _get_nodes(dim=-1,tag=-1,include_boundary=false,
                    return_parametric_coord=true)
    caller="API.mesh.get_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=LAST_MESH[]
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        boundary=_mesh_query_bool(include_boundary,caller,"include_boundary")
        parametric=_mesh_query_bool(
            return_parametric_coord,caller,"return_parametric_coord")
        if dimension<0
            # Gmsh returns empty arrays when the model holds no mesh rather
            # than erroring — `clear` and unmeshed models behave the same.
            node_tags=UInt64[]
            coordinates=Float64[]
            parameters=Float64[]
            if cached!==nothing
                append!(node_tags,UInt64.(1:nnodes(cached)))
                append!(coordinates,vec(cached.coords))
                overlay=_high_order_overlay(cached)
                if overlay!==nothing
                    n0=nnodes(cached)
                    append!(node_tags,UInt64.(n0+1:nnodes(overlay)))
                    append!(coordinates,
                            vec(@view overlay.coords[:,n0+1:end]))
                end
            end
            # Gmsh never returns parameters on the aggregate (-1,-1) query.
            for (edim,etag,record) in _discrete_mesh_records(model)
                _append_record_nodes!(
                    node_tags,coordinates,parameters,record,false)
                boundary && _append_record_boundary_nodes!(
                    node_tags,coordinates,parameters,model,record,edim,false)
            end
            return node_tags,coordinates,parameters
        end
        if entity>=0
            record=_model_mesh_record(model,dimension,entity)
            if haskey(model.discrete,(dimension,entity))
                node_tags=UInt64[]
                coordinates=Float64[]
                parameters=Float64[]
                _append_record_nodes!(
                    node_tags,coordinates,parameters,record,parametric)
                boundary && _append_record_boundary_nodes!(
                    node_tags,coordinates,parameters,model,record,
                    dimension,parametric)
                return node_tags,coordinates,parameters
            end
            if cached===nothing
                # Unmeshed model: Gmsh answers known entities with empty
                # arrays and rejects unknown tags. Attached records still
                # own their caller-assigned nodes without a cache.
                haskey(_mesh_entity_dictionary(model,dimension),entity) ||
                    haskey(model.meshing.attached,(dimension,entity)) ||
                    throw(ArgumentError(
                        "$caller: unknown " *
                        "$(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
                node_tags=UInt64[];coordinates=Float64[]
                parameters=Float64[]
                record===nothing || _append_record_nodes!(
                    node_tags,coordinates,parameters,record,parametric)
                return node_tags,coordinates,parameters
            end
            class=_cached_classification_locked(cached)
            class===nothing && throw(ArgumentError(
                "$caller: dimension-specific nodes require mesh " *
                "classification metadata; use dim=-1 for all cached nodes"))
            _mesh_classified_entity(model,class,dimension,entity,caller)
            indices=boundary ?
                _mesh_entity_nodes_with_boundary(
                    class,dimension,Int32(entity)) :
                _mesh_entity_node_positions(
                    class,(dimension,Int32(entity)))
            node_tags,coordinates=_mesh_nodes_payload(cached,indices)
            parameters=parametric ?
                _mesh_entity_parameters(model,cached,dimension,entity,indices) :
                Float64[]
            if record!==nothing
                _append_record_nodes!(
                    node_tags,coordinates,parameters,record,parametric)
            end
            _append_p2_midnodes!(node_tags,coordinates,parameters,
                                 model,cached,dimension,entity,parametric)
            return node_tags,coordinates,parameters
        end
        # All entities of `dimension`: classified entities first, then records.
        class=nothing
        node_tags=UInt64[]
        coordinates=Float64[]
        parameters=Float64[]
        if cached!==nothing
            class=_cached_classification_locked(cached)
            class===nothing && throw(ArgumentError(
                "$caller: dimension-specific nodes require mesh " *
                "classification metadata; use dim=-1 for all cached nodes"))
            for entity_tag in _mesh_classified_entities(class,dimension)
                group=boundary ?
                    _mesh_entity_nodes_with_boundary(
                        class,dimension,entity_tag) :
                    _mesh_entity_node_positions(
                        class,(dimension,entity_tag))
                nt,co=_mesh_nodes_payload(cached,group)
                append!(node_tags,nt)
                append!(coordinates,co)
                parametric && append!(parameters,_mesh_entity_parameters(
                    model,cached,dimension,Int(entity_tag),group))
                _append_p2_midnodes!(node_tags,coordinates,parameters,
                                     model,cached,dimension,
                                     Int(entity_tag),parametric)
            end
        end
        for (edim,etag,record) in _discrete_mesh_records(model)
            edim==dimension || continue
            _append_record_nodes!(
                node_tags,coordinates,parameters,record,parametric)
            if boundary && haskey(model.discrete,(edim,etag))
                _append_record_boundary_nodes!(
                    node_tags,coordinates,parameters,model,record,
                    dimension,parametric)
            end
        end
        return node_tags,coordinates,parameters
    end
end

function _get_node(node_tag)
    caller="API.mesh.get_node"
    return lock(STATE_LOCK) do
        model=_model_locked()
        tag=_mesh_query_integer(node_tag,caller,"node_tag")
        # Discrete and attached records carry caller-assigned tags — check them
        # first so a sparse tag resolves even without a generated mesh.
        tag32=0<tag<=typemax(Int32) ? Int32(tag) : Int32(-1)
        for (edim,etag,record) in _discrete_mesh_records(model)
            position=findfirst(==(tag32),record.node_tags)
            position===nothing && continue
            parameters=size(record.node_params,1)>0 ?
                vec(record.node_params[:,position]) : Float64[]
            return vec(record.node_coords[:,position]),parameters,edim,etag
        end
        cached=_cached_mesh_locked(caller)
        overlay=_high_order_overlay(cached)
        total=overlay===nothing ? nnodes(cached) : nnodes(overlay)
        (tag>=1 && tag<=total) || throw(ArgumentError(
            "$caller: unknown node $tag"))
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: node ownership requires mesh classification " *
            "metadata; generate a mesh so the cache owns entity ownership"))
        if tag>nnodes(cached)
            dimension,entity=LAST_MESH_HIGH_ORDER_MIDS[][tag-nnodes(cached)]
            parameters=dimension in (1,2) ? vec(model_parametrization(
                model,dimension,Int(entity),
                vec(overlay.coords[:,tag]))) : Float64[]
            return vec(overlay.coords[:,tag]),parameters,dimension,Int(entity)
        end
        dimension,entity=class.node_entities[tag]
        parameters=_mesh_entity_parameters(
            model,cached,dimension,Int(entity),[tag])
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
        model=_model_locked()
        cached=LAST_MESH[]
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        group=_mesh_query_integer(tag,caller,"tag")
        typemin(Int32)<=group<=typemax(Int32) || return UInt64[],Float64[]
        entities=get(model.physical,(dimension,group),Int[])
        isempty(entities) && return UInt64[],Float64[]
        node_tags=UInt64[]
        coordinates=Float64[]
        if cached!==nothing
            class=_cached_classification_locked(cached)
            class===nothing && throw(ArgumentError(
                "$caller: physical-group nodes require mesh classification " *
                "metadata; generate a mesh so the cache owns entity ownership"))
            indices=_mesh_physical_group_nodes(
                class,model,dimension,entities)
            node_tags,coordinates=_mesh_nodes_payload(cached,indices)
        end
        # Discrete member entities contribute their caller-tagged nodes along
        # with their declared boundary closure, sorted-unique with the
        # classified set.
        merged=Dict{UInt64,NTuple{3,Float64}}()
        for i in eachindex(node_tags)
            merged[node_tags[i]]=(
                coordinates[3i-2],coordinates[3i-1],coordinates[3i])
        end
        seen=Set{Tuple{Int,Int}}()
        queue=Tuple{Int,Int}[(dimension,entity) for entity in entities]
        while !isempty(queue)
            key=popfirst!(queue)
            key in seen && continue
            push!(seen,key)
            record=get(model.discrete,key,nothing)
            record===nothing && continue
            for i in eachindex(record.node_tags)
                merged[UInt64(record.node_tags[i])]=(
                    record.node_coords[1,i],record.node_coords[2,i],
                    record.node_coords[3,i])
            end
            append!(queue,Tuple{Int,Int}.(record.boundary))
        end
        ordered=sort!(collect(keys(merged)))
        return ordered,
               Float64[c for node in ordered for c in merged[node]]
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
        model=_model_locked()
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        if entity>=0 && dimension>=0 &&
                haskey(model.discrete,(dimension,entity))
            return _record_element_blocks(model.discrete[(dimension,entity)])
        end
        cached=LAST_MESH[]
        (dimension<0 || entity<0) &&
            return _merged_element_data(model,cached,dimension)
        if cached===nothing
            # Unmeshed model: Gmsh answers known entities with empty arrays
            # and rejects unknown tags; attached records still own their
            # caller-assigned elements.
            haskey(_mesh_entity_dictionary(model,dimension),entity) ||
                haskey(model.meshing.attached,(dimension,entity)) ||
                throw(ArgumentError(
                    "$caller: unknown " *
                    "$(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
            attached=get(model.meshing.attached,(dimension,entity),nothing)
            return attached===nothing ?
                (Int32[],Vector{UInt64}[],Vector{UInt64}[]) :
                _record_element_blocks(attached)
        end
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-specific elements require mesh classification " *
            "metadata; use a negative tag to query a complete dimension"))
        _mesh_classified_entity(model,class,dimension,entity,caller)
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
        attached=get(model.meshing.attached,(dimension,entity),nothing)
        if attached!==nothing
            rtypes,rtags,rnodes=_record_element_blocks(attached)
            for i in eachindex(rtypes)
                merge_into=findfirst(==(rtypes[i]),element_types)
                if merge_into===nothing
                    push!(element_types,rtypes[i])
                    push!(element_tags,rtags[i])
                    push!(node_tags,rnodes[i])
                else
                    append!(element_tags[merge_into],rtags[i])
                    append!(node_tags[merge_into],rnodes[i])
                end
            end
        end
        return element_types,element_tags,node_tags
    end
end

function _get_element_types(dim=-1,tag=-1)
    caller="API.mesh.get_element_types"
    return lock(STATE_LOCK) do
        model=_model_locked()
        dimension=_mesh_query_dimension(dim,caller)
        entity=_mesh_query_integer(tag,caller,"tag")
        if entity>=0 && dimension>=0 &&
                haskey(model.discrete,(dimension,entity))
            types,_,_=_record_element_blocks(model.discrete[(dimension,entity)])
            return types
        end
        cached=LAST_MESH[]
        if dimension<0 || entity<0
            result=cached===nothing ? Int32[] :
                _mesh_element_types(cached,dimension)
            overlay=_high_order_overlay(cached)
            if overlay!==nothing
                skeleton=_p2_skeleton_etype(overlay)
                for (i,element_type) in enumerate(result)
                    element_type==skeleton &&
                        (result[i]=_p2_etype(overlay))
                end
            end
            for (edim,etag,record) in _discrete_mesh_records(model)
                (dimension<0 || dimension==edim) || continue
                for element_type in record.element_types
                    element_type in result || push!(result,element_type)
                end
            end
            return sort!(result)
        end
        if cached===nothing
            haskey(_mesh_entity_dictionary(model,dimension),entity) ||
                haskey(model.meshing.attached,(dimension,entity)) ||
                throw(ArgumentError(
                    "$caller: unknown " *
                    "$(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
            attached=get(model.meshing.attached,(dimension,entity),nothing)
            return attached===nothing ? Int32[] :
                _record_element_blocks(attached)[1]
        end
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-specific element types require mesh " *
            "classification metadata; use a negative tag to query a complete " *
            "dimension"))
        _mesh_classified_entity(model,class,dimension,entity,caller)
        result=Int32[]
        overlay=_high_order_overlay(cached)
        tri_etype=overlay isa P2TriMesh ? Int32(9) : Int32(2)
        tet_etype=overlay isa P2Mesh ? Int32(11) : Int32(4)
        dimension==1 && any(==(Int32(entity)),class.seg_entities) &&
            push!(result,Int32(1))
        dimension==2 && any(==(Int32(entity)),class.tri_entities) &&
            push!(result,tri_etype)
        dimension==3 && any(==(Int32(entity)),class.tet_entities) &&
            push!(result,tet_etype)
        attached=get(model.meshing.attached,(dimension,entity),nothing)
        if attached!==nothing
            for element_type in attached.element_types
                element_type in result || push!(result,element_type)
            end
        end
        return result
    end
end

function _get_element(element_tag)
    caller="API.mesh.get_element"
    return lock(STATE_LOCK) do
        model=_model_locked()
        tag=_mesh_query_integer(element_tag,caller,"element_tag")
        # Discrete/attached element tags are caller-assigned and sparse —
        # resolve them before the dense generated range.
        tag32=0<tag<=typemax(Int32) ? Int32(tag) : Int32(-1)
        for (edim,etag,record) in _discrete_mesh_records(model)
            position=findfirst(==(tag32),record.element_tags)
            position===nothing && continue
            return record.element_types[position],
                UInt64.(record.element_nodes[position]),edim,etag
        end
        cached=_cached_mesh_locked(caller)
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
            overlay=_high_order_overlay(cached)
            overlay isa P2TriMesh && return Int32(9),
                UInt64.(vec(overlay.tri6[:,position])),2,
                Int(class.tri_entities[position])
            return Int32(2),UInt64.(vec(cached.tris[:,position])),2,
                Int(class.tri_entities[position])
        end
        position=tag-tetrahedron_offset
        overlay=_high_order_overlay(cached)
        overlay isa P2Mesh && return Int32(11),
            UInt64.(vec(overlay.tet10[:,position])),3,
            Int(class.tet_entities[position])
        return Int32(4),UInt64.(vec(cached.tets[:,position])),3,
            Int(class.tet_entities[position])
    end
end

# Element tags and flattened node tags of record elements matching `msh`,
# restricted to entity `tag` when nonnegative (any record dim — Gmsh permits
# elements of arbitrary type on any entity).
function _record_elements_of_type(m::GeoModel,msh::Int32,tag::Int)
    element_tags=UInt64[]
    node_tags=UInt64[]
    edim_filter=tag<0 ? -1 : Int(msh_spec(msh).dim)
    for (edim,etag,record) in _discrete_mesh_records(m)
        (tag<0 || (etag==tag && edim==edim_filter)) || continue
        for i in eachindex(record.element_tags)
            record.element_types[i]==msh || continue
            push!(element_tags,UInt64(record.element_tags[i]))
            append!(node_tags,UInt64.(record.element_nodes[i]))
        end
    end
    return element_tags,node_tags
end

function _get_elements_by_type(element_type,tag=-1,task=0,num_tasks=1)
    caller="API.mesh.get_elements_by_type"
    return lock(STATE_LOCK) do
        model=_model_locked()
        entity=_mesh_query_integer(tag,caller,"tag")
        msh=_mesh_query_integer(element_type,caller,"element_type")
        msh32=try Int32(msh) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: element_type exceeds the Int32 range"))
        end
        spec=msh_spec(msh) # validates the type exists
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        if entity>=0
            found=false
            for dimension in 0:3
                _model_entity_known(model,dimension,entity) &&
                    (found=true;break)
            end
            found || throw(ArgumentError(
                "$caller: unknown entity tag $entity"))
        end
        # Cache columns stay virtual until the task slice is known so a
        # partition only materializes its own share of the data.
        cached=LAST_MESH[]
        cells=nothing
        positions=nothing
        offset=0
        ncache=0
        if cached!==nothing
            class=_cached_classification_locked(cached)
            overlay=_high_order_overlay(cached)
            if overlay!==nothing && Int32(msh)==_p2_etype(overlay)
                offset=_p2_tag_offset(cached,overlay)
                cells=_p2_cells(overlay)
                if entity<0
                    ncache=size(cells,2)
                else
                    class===nothing && throw(ArgumentError(
                        "$caller: entity-specific data require mesh " *
                        "classification metadata; use a negative tag to " *
                        "query the complete cache"))
                    _mesh_classified_entity(
                        model,class,spec.dim,entity,caller)
                    positions=_mesh_entity_cell_positions(
                        class,Int(_p2_skeleton_etype(overlay)),entity)
                    ncache=length(positions)
                end
            else
                block=_mesh_element_block(cached,msh)
                if block!==nothing
                    offset,cells=block
                    if entity<0
                        ncache=size(cells,2)
                    else
                        class===nothing && throw(ArgumentError(
                            "$caller: entity-specific data require mesh " *
                            "classification metadata; use a negative tag to " *
                            "query the complete cache"))
                        haskey(model.discrete,(spec.dim,entity)) ||
                            _mesh_classified_entity(
                                model,class,spec.dim,entity,caller)
                        positions=_mesh_entity_cell_positions(class,msh,entity)
                        ncache=length(positions)
                    end
                end
            end
        end
        rtags,rnodes=_record_elements_of_type(model,msh32,entity)
        total=ncache+length(rtags)
        total==0 && return UInt64[],UInt64[]
        selected=_mesh_task_range(total,task_index,task_count)
        isempty(selected) && return UInt64[],UInt64[]
        first_selected,last_selected=first(selected),last(selected)
        nodes_per=Int(spec.nnodes)
        tags=UInt64[]
        node_tags=UInt64[]
        if first_selected<=ncache
            cache_end=min(last_selected,ncache)
            if positions===nothing
                columns=first_selected:cache_end
                append!(tags,UInt64.(offset .+ columns))
                append!(node_tags,UInt64.(vec(@view cells[:,columns])))
            else
                chosen=@view positions[first_selected:cache_end]
                append!(tags,UInt64.(offset .+ chosen))
                append!(node_tags,UInt64.(vec(@view cells[:,chosen])))
            end
        end
        if last_selected>ncache
            record_first=max(first_selected,ncache+1)-ncache
            record_last=last_selected-ncache
            append!(tags,@view rtags[record_first:record_last])
            append!(node_tags,@view rnodes[
                (record_first-1)*nodes_per+1:record_last*nodes_per])
        end
        return tags,node_tags
    end
end

# Nodes (tags, flat coords, flat params) used by record elements of type `msh`
# on entity `tag` — or every record when `tag` is negative. `msh` is a valid
# Int32 MSH type.
function _record_nodes_of_type(m::GeoModel,msh::Int32,tag::Int,
                               parametric::Bool)
    node_tags=UInt64[]
    coordinates=Float64[]
    parameters=Float64[]
    edim_filter=tag<0 ? -1 : Int(msh_spec(msh).dim)
    index_of=_record_node_index(m)
    for (edim,etag,record) in _discrete_mesh_records(m)
        (tag<0 || (etag==tag && edim==edim_filter)) || continue
        for i in eachindex(record.element_tags)
            record.element_types[i]==msh || continue
            for node in record.element_nodes[i]
                owner,position=index_of[node]
                push!(node_tags,UInt64(node))
                append!(coordinates,@view(owner.node_coords[:,position]))
                if parametric
                    # Params are reported at the element entity's dimension;
                    # a node owned by a different-dim record contributes a
                    # resized (zero-padded) parameter column.
                    column=owner.node_params[:,position]
                    if size(owner.node_params,1)==edim
                        append!(parameters,column)
                    else
                        for row in 1:edim
                            push!(parameters,
                                  row<=size(owner.node_params,1) ?
                                      column[row] : 0.0)
                        end
                    end
                end
            end
        end
    end
    return node_tags,coordinates,parameters
end

function _get_nodes_by_element_type(element_type,tag=-1,
                                    return_parametric_coord=true)
    caller="API.mesh.get_nodes_by_element_type"
    return lock(STATE_LOCK) do
        model=_model_locked()
        parametric=_mesh_query_bool(
            return_parametric_coord,caller,"return_parametric_coord")
        entity=_mesh_query_integer(tag,caller,"tag")
        msh32=try
            Int32(_mesh_query_integer(element_type,caller,"element_type"))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: element_type exceeds the Int32 range"))
        end
        cached=LAST_MESH[]
        node_tags=UInt64[]
        coordinates=Float64[]
        parameters=Float64[]
        if cached!==nothing
            _,block,positions=_mesh_query_type_block(
                cached,element_type,tag,caller)
            if block!==nothing
                _,cells=block
                selected=cells[:,_mesh_selected_columns(cells,positions)]
                nt,co=_mesh_nodes_for_cells(cached,selected)
                append!(node_tags,nt)
                append!(coordinates,co)
                if parametric
                    class=_cached_classification_locked(cached)
                    if class!==nothing
                        # Gmsh packs each repeated node's parameters on its
                        # owning entity — one `u` for Line owners, `(u, v)`
                        # for Plane owners, and nothing for Point or Volume
                        # owners — in entry order.
                        node_list=Int.(vec(selected))
                        groups=Dict{Tuple{Int,Int32},Vector{Int}}()
                        for node in unique(node_list)
                            owner=class.node_entities[node]
                            owner[1] in (1,2) || continue
                            push!(get!(groups,owner,Int[]),node)
                        end
                        node_parameters=Dict{Int,Vector{Float64}}()
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
            end
        end
        rtags,rcoords,rparams=_record_nodes_of_type(
            model,msh32,entity,parametric)
        append!(node_tags,rtags)
        append!(coordinates,rcoords)
        append!(parameters,rparams)
        return node_tags,coordinates,parameters
    end
end

# Node-tag lists of record elements of type `msh` on entity `tag` (all records
# when `tag` is negative).
function _record_cells_of_type(m::GeoModel,msh::Int32,tag::Int)
    cells=Vector{Int32}[]
    edim_filter=tag<0 ? -1 : Int(msh_spec(msh).dim)
    for (edim,etag,record) in _discrete_mesh_records(m)
        (tag<0 || (etag==tag && edim==edim_filter)) || continue
        for i in eachindex(record.element_tags)
            record.element_types[i]==msh && push!(cells,record.element_nodes[i])
        end
    end
    return cells
end

function _get_barycenters(element_type,tag,fast,primary,task=0,num_tasks=1)
    caller="API.mesh.get_barycenters"
    return lock(STATE_LOCK) do
        model=_model_locked()
        entity=_mesh_query_integer(tag,caller,"tag")
        msh=_mesh_query_integer(element_type,caller,"element_type")
        msh_spec(msh)
        msh32=try Int32(msh) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: element_type exceeds the Int32 range"))
        end
        fast_mode=_mesh_query_bool(fast,caller,"fast")
        _mesh_query_bool(primary,caller,"primary")
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        if entity>=0
            known=any(dimension->_model_entity_known(
                          model,dimension,entity),0:3)
            known || throw(ArgumentError(
                "$caller: unknown entity tag $entity"))
        end
        # Count cache and record elements before the task slice so only the
        # selected elements' barycenters are ever computed.
        cached=LAST_MESH[]
        cells=nothing
        columns=nothing
        ncache=0
        if cached!==nothing
            _,block,positions=_mesh_query_type_block(
                cached,element_type,tag,caller)
            if block!==nothing
                _,cells=block
                columns=_mesh_selected_columns(cells,positions)
                ncache=length(columns)
            end
        end
        record_dim=entity<0 ? -1 : Int(msh_spec(msh).dim)
        record_count=0
        for (edim,etag,record) in _discrete_mesh_records(model)
            (entity<0 || (etag==entity && edim==record_dim)) || continue
            record_count+=count(==(msh32),record.element_types)
        end
        selected=_mesh_task_range(
            ncache+record_count,task_index,task_count)
        isempty(selected) && return Float64[]
        first_selected,last_selected=first(selected),last(selected)
        result=Float64[]
        if first_selected<=ncache
            chosen=@view columns[first_selected:min(last_selected,ncache)]
            append!(result,_mesh_barycenters(
                cached,@view(cells[:,chosen]),fast_mode,caller))
        end
        if last_selected>ncache
            record_first=max(first_selected,ncache+1)-ncache
            record_last=last_selected-ncache
            seen=0
            for (edim,etag,record) in _discrete_mesh_records(model)
                (entity<0 || (etag==entity && edim==record_dim)) ||
                    continue
                matching=findall(==(msh32),record.element_types)
                lo=max(1,record_first-seen)
                hi=min(record_last-seen,length(matching))
                if lo>hi
                    seen+=length(matching)
                    continue
                end
                index_of=_record_node_index(model)
                for i in @view matching[lo:hi]
                    nodes=record.element_nodes[i]
                    weight=fast_mode ? 1.0 :
                           inv(Float64(length(nodes)))
                    bx=by=bz=0.0
                    for node in nodes
                        owner,position=index_of[node]
                        bx+=weight*owner.node_coords[1,position]
                        by+=weight*owner.node_coords[2,position]
                        bz+=weight*owner.node_coords[3,position]
                    end
                    all(isfinite,(bx,by,bz)) || throw(ArgumentError(
                        "$caller: element coordinate sum is not " *
                        "Float64-representable"))
                    append!(result,(bx,by,bz))
                end
                seen+=length(matching)
                seen>=record_last && break
            end
        end
        return result
    end
end

# `patterns` applied to record-element node-tag lists — the record equivalent
# of `_mesh_pattern_nodes`.
function _record_pattern_nodes(cells::Vector{Vector{Int32}},patterns)
    isempty(patterns) && return UInt64[]
    result=UInt64[]
    for cell in cells,pattern in patterns,local_node in pattern
        push!(result,UInt64(cell[local_node]))
    end
    return result
end

function _get_element_edge_nodes(element_type,tag=-1,primary=false,
                                 task=0,num_tasks=1)
    caller="API.mesh.get_element_edge_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        entity=_mesh_query_integer(tag,caller,"tag")
        msh=_mesh_query_integer(element_type,caller,"element_type")
        msh32=try Int32(msh) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: element_type exceeds the Int32 range"))
        end
        _mesh_query_bool(primary,caller,"primary")
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        patterns=_simplex_edge_patterns(msh)
        cached=LAST_MESH[]
        cache_columns=0:-1
        cells=nothing
        if cached!==nothing
            _,block,positions=_mesh_query_type_block(
                cached,element_type,tag,caller)
            if block!==nothing
                _,cells=block
                cache_columns=_mesh_selected_columns(cells,positions)
            end
        end
        record_cells=_record_cells_of_type(model,msh32,entity)
        count=length(cache_columns)+length(record_cells)
        isempty(patterns) && return UInt64[]
        selected=_mesh_task_range(count,task_index,task_count)
        isempty(selected) && return UInt64[]
        output=UInt64[]
        cache_count=length(cache_columns)
        cache_selected=selected[1]:min(selected[end],cache_count)
        if !isempty(cache_selected) && cache_selected[1]<=cache_count
            append!(output,_mesh_pattern_nodes(
                @view(cells[:,cache_columns[cache_selected]]),patterns))
        end
        record_first=max(selected[1],cache_count+1)
        record_first<=selected[end] && append!(output,_record_pattern_nodes(
            record_cells[record_first-cache_count:selected[end]-cache_count],
            patterns))
        return output
    end
end

function _get_element_face_nodes(element_type,face_type,tag=-1,primary=false,
                                 task=0,num_tasks=1)
    caller="API.mesh.get_element_face_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        entity=_mesh_query_integer(tag,caller,"tag")
        msh=_mesh_query_integer(element_type,caller,"element_type")
        msh32=try Int32(msh) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: element_type exceeds the Int32 range"))
        end
        face=_mesh_query_integer(face_type,caller,"face_type")
        face in (3,4) || throw(ArgumentError(
            "$caller: face_type must be 3 (triangle) or 4 (quadrangle)"))
        _mesh_query_bool(primary,caller,"primary")
        task_index,task_count=_mesh_query_tasks(task,num_tasks,caller)
        patterns=_simplex_face_patterns(msh,face)
        cached=LAST_MESH[]
        cache_columns=0:-1
        cells=nothing
        if cached!==nothing
            _,block,positions=_mesh_query_type_block(
                cached,element_type,tag,caller)
            if block!==nothing
                _,cells=block
                cache_columns=_mesh_selected_columns(cells,positions)
            end
        end
        record_cells=_record_cells_of_type(model,msh32,entity)
        count=length(cache_columns)+length(record_cells)
        isempty(patterns) && return UInt64[]
        selected=_mesh_task_range(count,task_index,task_count)
        isempty(selected) && return UInt64[]
        output=UInt64[]
        cache_count=length(cache_columns)
        cache_selected=selected[1]:min(selected[end],cache_count)
        if !isempty(cache_selected) && cache_selected[1]<=cache_count
            append!(output,_mesh_pattern_nodes(
                @view(cells[:,cache_columns[cache_selected]]),patterns))
        end
        record_first=max(selected[1],cache_count+1)
        record_first<=selected[end] && append!(output,_record_pattern_nodes(
            record_cells[record_first-cache_count:selected[end]-cache_count],
            patterns))
        return output
    end
end

function _get_max_node_tag()
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=LAST_MESH[]
        overlay=_high_order_overlay(cached)
        maximum_tag=overlay===nothing ?
            (cached===nothing ? 0 : nnodes(cached)) : nnodes(overlay)
        for (_,_,record) in _discrete_mesh_records(model)
            isempty(record.node_tags) ||
                (maximum_tag=max(maximum_tag,Int(maximum(record.node_tags))))
        end
        return UInt64(maximum_tag)
    end
end

function _get_max_element_tag()
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=LAST_MESH[]
        _,_,total=cached===nothing ? (0,0,0) : _mesh_element_offsets(cached)
        maximum_tag=total
        for (_,_,record) in _discrete_mesh_records(model)
            isempty(record.element_tags) ||
                (maximum_tag=max(maximum_tag,Int(maximum(record.element_tags))))
        end
        return UInt64(maximum_tag)
    end
end

# Re-derive the classification after uniform refinement without a geometric
# re-projection: each refined cell inherits its parent cell's owner (children
# are emitted in fixed groups of 2/4/8 per parent), kept nodes keep their
# owner, and each new edge-midpoint node takes the owner of its skeleton edge
# under the same endpoint rule the order-2 overlay uses.
function _inherit_refined_classification(class::_MeshClassification,
                                         refined::Mesh,cache::Mesh)
    mesh=class.mesh
    function inherit(cells::AbstractMatrix{Int32},entities,fanout)
        isempty(entities) && return entities
        result=Vector{Int32}(undef,size(cells,2))
        for cell in axes(cells,2)
            result[cell]=entities[(cell-1)÷fanout+1]
        end
        return result
    end
    seg_entities=inherit(refined.segs,class.seg_entities,2)
    tri_entities=inherit(refined.tris,class.tri_entities,4)
    tet_entities=inherit(refined.tets,class.tet_entities,8)
    node_entities=Vector{Tuple{Int,Int32}}(undef,nnodes(refined))
    assigned=falses(nnodes(refined))
    # `refine_uniform` compacts referenced nodes in order then appends edge
    # midpoints — surviving nodes match by coordinate.
    linear_lookup=Dict{NTuple{3,Int64},Int32}()
    for node in 1:nnodes(mesh)
        key=ntuple(axis->round(Int64,mesh.coords[axis,node]*1e12),3)
        linear_lookup[key]=Int32(node)
    end
    for node in 1:nnodes(refined)
        key=ntuple(axis->round(Int64,refined.coords[axis,node]*1e12),3)
        source=get(linear_lookup,key,Int32(0))
        source==0 && continue
        node_entities[node]=class.node_entities[source]
        assigned[node]=true
    end
    # A midpoint's owner: collect its corner (non-midpoint) neighbors across
    # the child cells it touches — that set is exactly its skeleton edge's two
    # endpoints; apply the same endpoint-ownership rule as `_p2_mid_owners`.
    entity_cells=((1,refined.segs,seg_entities,((1,2),)),
                  (2,refined.tris,tri_entities,((1,2),(2,3),(3,1))),
                  (3,refined.tets,tet_entities,
                   ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))))
    neighbors=[Set{Int}() for _ in 1:nnodes(refined)]
    for (_,cells,entities,edge_patterns) in entity_cells
        isempty(entities) && continue
        for cell in axes(cells,2)
            corners=Int.(cells[:,cell])
            for (i,j) in edge_patterns
                a,b=corners[i],corners[j]
                if assigned[a] != assigned[b]
                    midpoint,endpoint=assigned[a] ? (b,a) : (a,b)
                    push!(neighbors[midpoint],endpoint)
                end
            end
        end
    end
    for (dimension,cells,entities,_) in entity_cells
        for cell in axes(cells,2)
            for corner in Int.(cells[:,cell])
                assigned[corner] && continue
                endpoints=neighbors[corner]
                isempty(endpoints) && continue
                owners=[node_entities[e] for e in endpoints]
                first_owner=owners[1]
                if all(==(first_owner),owners)
                    node_entities[corner]=first_owner
                elseif !all(o->o[1]==first_owner[1],owners)
                    top=argmax(o->o[1],owners)
                    node_entities[corner]=owners[top]
                else
                    node_entities[corner]=(dimension,entities[cell])
                end
                assigned[corner]=true
            end
        end
    end
    fallback=class.entities[1]
    for node in 1:nnodes(refined)
        assigned[node] || (node_entities[node]=fallback)
    end
    return _MeshClassification(cache,class.entity,class.entities,
        node_entities,class.boundaries,seg_entities,tri_entities,
        tet_entities)
end

function _refine(;max_nodes=typemax(Int32),max_cells=typemax(Int32))
    return lock(STATE_LOCK) do
        m=_model_locked()
        cached=LAST_MESH[]
        cached===nothing && throw(ArgumentError(
            "API.mesh.refine: no mesh; call API.mesh.generate first"))
        overlay=_high_order_overlay(cached)
        # `reverse` may legitimately invert tets; refinement preserves
        # orientation, so an inverted input yields an inverted (but otherwise
        # structurally valid) output — bypass only the orientation check.
        refined=refine_uniform(
            cached;max_nodes=max_nodes,max_cells=max_cells,
            require_positive_tets=false)
        class=_cached_classification_locked(cached)
        cache=_copy_mesh(refined)
        new_class=if class===nothing
            nothing
        elseif length(class.entities)<=1
            _classify_cached_mesh(
                m,refined,class.entity[1],Int(class.entity[2]),cache)
        else
            _inherit_refined_classification(class,refined,cache)
        end
        _replace_mesh_cache_locked!(cache,new_class)
        overlay!==nothing &&
            _rebind_high_order!(cache,new_class,"API.mesh.refine")
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
    record=_MeshClassification(replacement,class.entity,class.entities,node_entities,
                               class.boundaries,cell_entities[1],
                               cell_entities[2],cell_entities[3])
    return replacement,record
end

# Drop every node and element stored on a discrete/attached record, keeping
# the entity itself (Gmsh's `clear` removes mesh data, not entities).
function _record_clear!(record)
    empty!(record.node_tags)
    record.node_coords=Matrix{Float64}(undef,3,0)
    record.node_params=Matrix{Float64}(undef,size(record.node_params,1),0)
    empty!(record.element_types)
    empty!(record.element_tags)
    empty!(record.element_nodes)
    return nothing
end

function _clear_mesh(dim_tags=())
    caller="API.mesh.clear"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        if isempty(pairs)
            _replace_mesh_cache_locked!(nothing)
            for (_,_,record) in _discrete_mesh_records(model)
                _record_clear!(record)
            end
            return nothing
        end
        for (dim,tag) in pairs
            record=_model_mesh_record(model,dim,tag)
            record===nothing || _record_clear!(record)
        end
        cached=LAST_MESH[]
        cached===nothing && return nothing
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: entity-selective clearing requires mesh classification " *
            "metadata; pass an empty collection to clear the complete cache"))
        cleared=_mesh_selected_entities(model,class,pairs,caller)
        overlay=_high_order_overlay(cached)
        replacement,record=_clear_classified_mesh(cached,class,cleared)
        if replacement!==cached
            _replace_mesh_cache_locked!(replacement,record)
            overlay!==nothing &&
                _rebind_high_order!(replacement,record,caller)
        end
        nothing
    end
end

function _affine_transform_mesh(affine,dim_tags=())
    caller="API.mesh.affine_transform"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=LAST_MESH[]
        cached===nothing && isempty(model.discrete) &&
            isempty(model.meshing.attached) && throw(ArgumentError(
                "$caller: no mesh; call API.mesh.generate first"))
        coefficients,translation,_=_transform_gmsh_affine(affine,caller)
        matrix=reshape(collect(coefficients),3,3)
        # Discrete and attached records transform with (or instead of) the
        # generated cache: all records on the global form, selected entities
        # on the entity-selective form.
        selected_record_keys=nothing
        if !isempty(pairs)
            selected_record_keys=Set{Tuple{Int,Int}}(
                (Int(p[1]),Int(p[2])) for p in pairs)
        end
        orientation=matrix[1,1]*(matrix[2,2]*matrix[3,3]-
                                matrix[2,3]*matrix[3,2])-
                    matrix[1,2]*(matrix[2,1]*matrix[3,3]-
                                matrix[2,3]*matrix[3,1])+
                    matrix[1,3]*(matrix[2,1]*matrix[3,2]-
                                matrix[2,2]*matrix[3,1])
        for (edim,etag,record) in _discrete_mesh_records(model)
            (selected_record_keys===nothing ||
             (edim,etag) in selected_record_keys) || continue
            record.node_coords=matrix*record.node_coords .+
                               reshape(collect(translation),3,1)
            if orientation<0
                for i in eachindex(record.element_nodes)
                    _record_reverse_connectivity!(
                        record.element_types[i],record.element_nodes[i])
                end
            end
        end
        cached===nothing && return Mesh(Matrix{Float64}(undef,3,0))
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
        overlay=_high_order_overlay(cached)
        class=_cached_classification_locked(cached)
        cache=_copy_mesh(transformed)
        new_class=class===nothing ? nothing : _MeshClassification(
            cache,class.entity,class.entities,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(cache,new_class)
        overlay!==nothing && _rebind_high_order!(cache,new_class,caller)
        transformed
    end
end

# Rebuild the cache after dropping whole element columns from one block. Nodes
# and their entity ownership are retained (Gmsh 4.15.2 keeps orphan nodes);
# only the selected block's connectivity, per-cell owners, and tags shrink.
function _remove_element_columns(mesh::Mesh,class::_MeshClassification,
                                 dimension::Int,dropped::AbstractVector{Bool})
    overlay=_high_order_overlay(mesh)
    keep=.!dropped
    replacement=Mesh(mesh.coords;
        segs=dimension==1 ? mesh.segs[:,keep] : mesh.segs,
        tris=dimension==2 ? mesh.tris[:,keep] : mesh.tris,
        tets=dimension==3 ? mesh.tets[:,keep] : mesh.tets,
        seg_tag=dimension==1 ? mesh.seg_tag[keep] : mesh.seg_tag,
        tri_tag=dimension==2 ? mesh.tri_tag[keep] : mesh.tri_tag,
        tet_tag=dimension==3 ? mesh.tet_tag[keep] : mesh.tet_tag)
    record=_MeshClassification(replacement,class.entity,class.entities,class.node_entities,
        class.boundaries,
        dimension==1 ? class.seg_entities[keep] : class.seg_entities,
        dimension==2 ? class.tri_entities[keep] : class.tri_entities,
        dimension==3 ? class.tet_entities[keep] : class.tet_entities)
    _replace_mesh_cache_locked!(replacement,record)
    overlay!==nothing &&
        _rebind_high_order!(replacement,record,"API.mesh.remove_elements")
    return nothing
end

# Remove caller-tagged elements from a discrete/attached record. An empty
# `element_tags` removes every element on the entity; listed tags must exist.
function _record_remove_elements!(record,element_tags,caller::AbstractString)
    if isempty(element_tags)
        empty!(record.element_types)
        empty!(record.element_tags)
        empty!(record.element_nodes)
        return nothing
    end
    keep=trues(length(record.element_tags))
    for value in element_tags
        entry=_mesh_query_integer(value,caller,"element_tags entry")
        entry32=typemin(Int32)<=entry<=typemax(Int32) ? Int32(entry) :
            Int32(0)
        position=findfirst(==(entry32),record.element_tags)
        position===nothing && throw(ArgumentError(
            "$caller: element $entry is not classified on this entity"))
        keep[position]=false
    end
    record.element_types=record.element_types[keep]
    record.element_tags=record.element_tags[keep]
    record.element_nodes=record.element_nodes[keep]
    return nothing
end

function _remove_elements(dim,tag,element_tags=())
    caller="API.mesh.remove_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        dimension=_mesh_query_integer(dim,caller,"dim")
        dimension in 0:3 || throw(ArgumentError(
            "$caller: dim must be in 0:3"))
        entity=_mesh_query_integer(tag,caller,"tag")
        (element_tags isa AbstractVector || element_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: element_tags must be a vector or tuple of " *
                "element tags"))
        record=get(model.discrete,(dimension,entity),nothing)
        if record!==nothing
            _record_remove_elements!(record,element_tags,caller)
            return nothing
        end
        attached=get(model.meshing.attached,(dimension,entity),nothing)
        if attached!==nothing
            if isempty(element_tags)
                _record_remove_elements!(attached,element_tags,caller)
            else
                listed=Set{Int32}()
                for value in element_tags
                    entry=_mesh_query_integer(value,caller,
                                              "element_tags entry")
                    typemin(Int32)<=entry<=typemax(Int32) &&
                        push!(listed,Int32(entry))
                end
                record_listed=intersect(listed,Set(attached.element_tags))
                if !isempty(record_listed)
                    _record_remove_elements!(
                        attached,[Int(t) for t in record_listed],caller)
                    element_tags=Any[v for v in element_tags
                                     if !(v in record_listed)]
                    # A nonempty selection consumed entirely by record
                    # elements must not fall through to the remove-all
                    # semantics of an empty list.
                    isempty(element_tags) && return nothing
                end
            end
        end
        cached=_cached_mesh_locked(caller)
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

# Reverse one record element's connectivity for `reverse`/`affineTransform`
# with a negative determinant. First-order types follow Gmsh's per-type swaps;
# other types keep the first node and reverse the remainder, which flips the
# orientation of every linear face convention.
function _record_reverse_connectivity!(msh::Int32,nodes::Vector{Int32})
    n=length(nodes)
    n<2 && return nothing
    if msh==Int32(1)                      # 2-node line
        nodes[1],nodes[2]=nodes[2],nodes[1]
    elseif msh==Int32(2)                  # triangle
        nodes[2],nodes[3]=nodes[3],nodes[2]
    elseif msh==Int32(3)                  # quadrangle
        nodes[2],nodes[4]=nodes[4],nodes[2]
    elseif msh==Int32(4)                  # tetrahedron
        nodes[1],nodes[2]=nodes[2],nodes[1]
    elseif msh==Int32(5)                  # hexahedron
        nodes[2],nodes[4]=nodes[4],nodes[2]
        nodes[6],nodes[8]=nodes[8],nodes[6]
    elseif msh==Int32(6)                  # prism
        nodes[2],nodes[3]=nodes[3],nodes[2]
        nodes[5],nodes[6]=nodes[6],nodes[5]
    elseif msh==Int32(7)                  # pyramid
        nodes[2],nodes[4]=nodes[4],nodes[2]
    else
        nodes[2:end]=reverse(nodes[2:end])
    end
    return nothing
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
        cached=LAST_MESH[]
        class=cached===nothing ? nothing :
            _cached_classification_locked(cached)
        isempty(pairs) || (cached===nothing || class===nothing) &&
            any(pair->!haskey(model.discrete,pair) &&
                      !haskey(model.meshing.attached,pair),pairs) &&
            throw(ArgumentError(
                "$caller: entity-selective reversal requires mesh " *
                "classification metadata; pass an empty collection to " *
                "reverse the complete cache"))
        cached===nothing && isempty(model.discrete) &&
            isempty(model.meshing.attached) && throw(ArgumentError(
                "$caller: no mesh; call API.mesh.generate first"))
        selected=(isempty(pairs) || class===nothing) ? nothing :
                 _mesh_selected_entities(model,class,pairs,caller)
        record_keys=isempty(pairs) ? nothing :
                    Set{Tuple{Int,Int}}(pairs)
        for (edim,etag,record) in _discrete_mesh_records(model)
            (record_keys===nothing || (edim,etag) in record_keys) ||
                continue
            for i in eachindex(record.element_nodes)
                _record_reverse_connectivity!(
                    record.element_types[i],record.element_nodes[i])
            end
        end
        cached===nothing && return nothing
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
        overlay=_high_order_overlay(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.entities,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        overlay!==nothing &&
            _rebind_high_order!(replacement,new_class,caller)
        return nothing
    end
end

function _reverse_elements(element_tags)
    caller="API.mesh.reverse_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=LAST_MESH[]
        (element_tags isa AbstractVector || element_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: element_tags must be a vector or tuple of " *
                "element tags"))
        pending=Int[_mesh_query_integer(value,caller,"element_tags entry")
                    for value in element_tags]
        # Caller-assigned record element tags resolve first.
        if !isempty(pending)
            remaining=Int[]
            for dense in pending
                found=false
                dense32=typemin(Int32)<=dense<=typemax(Int32) ?
                        Int32(dense) : Int32(0)
                for (_,_,record) in _discrete_mesh_records(model)
                    position=findfirst(==(dense32),record.element_tags)
                    position===nothing && continue
                    _record_reverse_connectivity!(
                        record.element_types[position],
                        record.element_nodes[position])
                    found=true
                    break
                end
                found || push!(remaining,dense)
            end
            pending=remaining
        end
        isempty(pending) && return nothing
        cached===nothing && throw(ArgumentError(
            "$caller: no mesh; call API.mesh.generate first"))
        triangle_offset,tetrahedron_offset,total=_mesh_element_offsets(cached)
        per_block=(Int[],Int[],Int[])
        for dense in pending
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
        overlay=_high_order_overlay(cached)
        class=_cached_classification_locked(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.entities,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        overlay!==nothing &&
            _rebind_high_order!(replacement,new_class,caller)
        return nothing
    end
end

# Merge record nodes sharing exact coordinates. `groups` maps a coordinate
# key to the surviving tag — record members keep the lowest record tag, while
# a coincident cache member always wins because record connectivity can store
# any tag but dense cells only reference cache positions. `pairs` restricts
# the scan to records on the listed entities.
function _record_dedup_nodes!(records,pairs,cache_keepers)
    selected=pairs===nothing ? nothing :
             Set{Tuple{Int,Int}}(pairs)
    # coord key → (tag, is_cache)
    keepers=Dict{NTuple{3,Float64},Tuple{Int32,Bool}}()
    cache_keepers===nothing || merge!(keepers,cache_keepers)
    remap=Dict{Int32,Int32}()
    deleted=Set{Int32}()
    for (edim,etag,record) in records
        (selected===nothing || (edim,etag) in selected) || continue
        for i in eachindex(record.node_tags)
            key=(record.node_coords[1,i],record.node_coords[2,i],
                 record.node_coords[3,i])
            existing=get(keepers,key,nothing)
            tag=record.node_tags[i]
            if existing===nothing
                keepers[key]=(tag,false)
            elseif existing[2]
                # Cache member wins regardless of tag ordering.
                remap[tag]=existing[1]
                push!(deleted,tag)
            else
                if tag<existing[1]
                    remap[existing[1]]=tag
                    push!(deleted,existing[1])
                    keepers[key]=(tag,false)
                else
                    remap[tag]=existing[1]
                    push!(deleted,tag)
                end
            end
        end
    end
    isempty(deleted) && return false
    for (edim,etag,record) in records
        for connectivity in record.element_nodes
            for i in eachindex(connectivity)
                connectivity[i]=get(remap,connectivity[i],connectivity[i])
            end
        end
        keep=trues(length(record.node_tags))
        for i in eachindex(record.node_tags)
            record.node_tags[i] in deleted && (keep[i]=false)
        end
        all(keep) && continue
        record.node_tags=record.node_tags[keep]
        record.node_coords=record.node_coords[:,keep]
        record.node_params=record.node_params[:,keep]
    end
    return true
end

# Merge nodes sharing exact coordinates, keeping the lowest tag, then rebuild
# the cache with compacted node numbering and remapped connectivity — matching
# Gmsh 4.15.2's `removeDuplicateNodes`.
function _remove_duplicate_nodes(dim_tags=())
    caller="API.mesh.remove_duplicate_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=LAST_MESH[]
        records=_discrete_mesh_records(model)
        if cached===nothing
            isempty(records) && throw(ArgumentError(
                "$caller: no mesh; call API.mesh.generate first"))
            record_pairs=isempty(pairs) ? nothing :
                          collect(Tuple{Int,Int},pairs)
            _record_dedup_nodes!(records,record_pairs,nothing)
            return nothing
        end
        class=_cached_classification_locked(cached)
        (isempty(pairs) || class!==nothing ||
         all(pair->haskey(model.discrete,pair) ||
                  haskey(model.meshing.attached,pair),pairs)) ||
            throw(ArgumentError(
                "$caller: entity-selective duplicate removal requires mesh " *
                "classification metadata; pass an empty collection to scan " *
                "the complete cache"))
        selected=(isempty(pairs) || class===nothing) ? nothing :
                 _mesh_selected_entities(model,class,pairs,caller)
        # An unclassified cache cannot filter by entity: with a nonempty
        # record-only selection no cache node is admissible.
        scan_cache=class!==nothing || isempty(pairs)
        record_pairs=isempty(pairs) ? nothing :
                      collect(Tuple{Int,Int},pairs)
        count=nnodes(cached)
        replacement=collect(Int32,1:count)
        groups=Dict{NTuple{3,Float64},Int32}()
        for node in 1:count
            (scan_cache && (selected===nothing ||
             class.node_entities[node] in selected)) || continue
            key=NTuple{3,Float64}((cached.coords[1,node],
                                   cached.coords[2,node],
                                   cached.coords[3,node]))
            replacement[node]=Int32(get!(groups,key,node))
        end
        keep=Bool[replacement[node]==node for node in 1:count]
        if !all(keep)
            old_to_new=Vector{Int32}(undef,count)
            index=0
            for node in 1:count
                keep[node] && (index+=1;old_to_new[node]=index)
            end
            # Record connectivity referencing dense cache tags follows the
            # same survivor/compaction remap.
            for (_,_,record) in records
                for connectivity in record.element_nodes
                    for i in eachindex(connectivity)
                        tag=connectivity[i]
                        1<=tag<=count || continue
                        connectivity[i]=old_to_new[Int(replacement[tag])]
                    end
                end
            end
            # Surviving cache nodes keep their new dense tags; coincident
            # record nodes merge onto them (dense positions cannot hold
            # caller-assigned tags, so the cache member must win).
            cache_keepers=Dict{NTuple{3,Float64},Tuple{Int32,Bool}}()
            for node in 1:count
                keep[node] || continue
                key=NTuple{3,Float64}((cached.coords[1,node],
                                       cached.coords[2,node],
                                       cached.coords[3,node]))
                cache_keepers[key]=(old_to_new[node],true)
            end
            _record_dedup_nodes!(records,record_pairs,cache_keepers)
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
            overlay=_high_order_overlay(cached)
            new_class=class===nothing ? nothing : _MeshClassification(
                new_mesh,class.entity,class.entities,node_entities,class.boundaries,
                class.seg_entities,class.tri_entities,class.tet_entities)
            _replace_mesh_cache_locked!(new_mesh,new_class)
            overlay!==nothing &&
                _rebind_high_order!(new_mesh,new_class,caller)
            return nothing
        end
        # No cache duplicates: still merge records against unchanged cache
        # coordinates when the selection admits both.
        isempty(records) && return nothing
        cache_keepers=Dict{NTuple{3,Float64},Tuple{Int32,Bool}}()
        for node in 1:count
            (scan_cache && (selected===nothing ||
             class.node_entities[node] in selected)) || continue
            key=NTuple{3,Float64}((cached.coords[1,node],
                                   cached.coords[2,node],
                                   cached.coords[3,node]))
            cache_keepers[key]=(Int32(node),true)
        end
        _record_dedup_nodes!(records,record_pairs,cache_keepers)
        return nothing
    end
end

# Drop cells whose sorted connectivity repeats an earlier cell owned by the
# same entity — Gmsh 4.15.2's `removeDuplicateElements` contract.
# Dedup one record's elements by sorted connectivity — the record analogue of
# the per-entity cache pass. `seed` optionally carries connectivity tuples of
# cache cells already owned by the same entity (attached-record case).
function _record_dedup_elements!(record,seed=nothing)
    seen=Set{Vector{Int32}}()
    seed===nothing || union!(seen,seed)
    keep=trues(length(record.element_tags))
    changed=false
    for i in eachindex(record.element_tags)
        key=sort!(copy(record.element_nodes[i]))
        if key in seen
            keep[i]=false
            changed=true
        else
            push!(seen,key)
        end
    end
    changed || return false
    record.element_types=record.element_types[keep]
    record.element_tags=record.element_tags[keep]
    record.element_nodes=record.element_nodes[keep]
    return true
end

function _remove_duplicate_elements(dim_tags=())
    caller="API.mesh.remove_duplicate_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=LAST_MESH[]
        records=_discrete_mesh_records(model)
        if cached===nothing
            isempty(records) && throw(ArgumentError(
                "$caller: no mesh; call API.mesh.generate first"))
            selected=isempty(pairs) ? nothing :
                      Set{Tuple{Int,Int}}(pairs)
            for (edim,etag,record) in records
                (selected===nothing || (edim,etag) in selected) || continue
                _record_dedup_elements!(record)
            end
            return nothing
        end
        class=_cached_classification_locked(cached)
        class===nothing && throw(ArgumentError(
            "$caller: duplicate-element removal requires mesh " *
            "classification metadata; generate a mesh so the cache owns " *
            "entity ownership"))
        selected=isempty(pairs) ? nothing :
                 _mesh_selected_entities(model,class,pairs,caller)
        record_selected=isempty(pairs) ? nothing :
                        Set{Tuple{Int,Int}}(pairs)
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
                    slot->slot<=length(connectivity) ?
                        connectivity[slot] : Int32(0),4))
                haskey(seen,key) && (keep[column]=false;changed=true)
                seen[key]=column
            end
            new_blocks[block]=cells[:,keep]
            new_owners[block]=owners[keep]
            new_tags[block]=cell_tags[block][keep]
        end
        if changed
            replacement=Mesh(cached.coords;segs=new_blocks[1],
                             tris=new_blocks[2],tets=new_blocks[3],
                             seg_tag=new_tags[1],tri_tag=new_tags[2],
                             tet_tag=new_tags[3])
            overlay=_high_order_overlay(cached)
            new_class=_MeshClassification(replacement,class.entity,class.entities,
                class.node_entities,class.boundaries,new_owners[1],
                new_owners[2],new_owners[3])
            _replace_mesh_cache_locked!(replacement,new_class)
            overlay!==nothing &&
                _rebind_high_order!(replacement,new_class,caller)
        end
        # Record elements dedup per entity; an attached record on a native
        # entity is additionally seeded with that entity's surviving cache
        # connectivity so mixed duplicates collapse.
        class2=LAST_MESH[]===nothing ? nothing :
               _cached_classification_locked(LAST_MESH[])
        for (edim,etag,record) in records
            (record_selected===nothing ||
             (edim,etag) in record_selected) || continue
            seed=nothing
            if !haskey(model.discrete,(edim,etag)) && class2!==nothing
                seed=Vector{Int32}[]
                for (dimension,cells,owners) in (
                        (1,LAST_MESH[].segs,class2.seg_entities),
                        (2,LAST_MESH[].tris,class2.tri_entities),
                        (3,LAST_MESH[].tets,class2.tet_entities))
                    dimension==edim || continue
                    for column in axes(cells,2)
                        owners[column]==Int32(etag) && push!(
                            seed,sort!(vec(cells[:,column])))
                    end
                end
            end
            _record_dedup_elements!(record,seed)
        end
        return nothing
    end
end

function _get_duplicate_nodes(dim_tags=())
    caller="API.mesh.get_duplicate_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        pairs=_mesh_parse_dim_tags(dim_tags,caller)
        cached=LAST_MESH[]
        class=cached===nothing ? nothing :
            _cached_classification_locked(cached)
        selected=nothing
        selected_records=nothing
        if !isempty(pairs)
            (cached===nothing || class===nothing) &&
                any(pair->!haskey(model.discrete,pair) &&
                          !haskey(model.meshing.attached,pair),pairs) &&
                throw(ArgumentError(
                    "$caller: entity-selective duplicate detection requires " *
                    "mesh classification metadata for generated entities"))
            selected=cached===nothing || class===nothing ? Set{Tuple{Int,Int32}}() :
                _mesh_selected_entities(model,class,pairs,caller)
            selected_records=Set{Tuple{Int,Int}}(
                pair for pair in pairs)
        end
        groups=Dict{NTuple{3,Float64},Vector{UInt64}}()
        if cached!==nothing
            for node in axes(cached.coords,2)
                (selected===nothing ||
                 class.node_entities[node] in selected) || continue
                push!(get!(groups,Tuple(cached.coords[:,node]),UInt64[]),
                      UInt64(node))
            end
        end
        for (edim,etag,record) in _discrete_mesh_records(model)
            (selected_records===nothing || (edim,etag) in selected_records) ||
                continue
            for i in eachindex(record.node_tags)
                push!(get!(groups,(
                    record.node_coords[1,i],record.node_coords[2,i],
                    record.node_coords[3,i]),UInt64[]),
                      UInt64(record.node_tags[i]))
            end
        end
        duplicates=UInt64[]
        for members in values(groups)
            length(members)>1 && append!(duplicates,members)
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
        model=_model_locked()
        tag=_mesh_query_integer(node_tag,caller,"node_tag")
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
        # Caller-assigned record nodes resolve before the dense cache so a
        # sparse tag updates its stored coordinate row in place.
        tag32=0<tag<=typemax(Int32) ? Int32(tag) : Int32(-1)
        for (_,_,record) in _discrete_mesh_records(model)
            position=findfirst(==(tag32),record.node_tags)
            position===nothing && continue
            record.node_coords[:,position]=values
            if !isempty(parametric_coord)
                length(parametric_coord)==size(record.node_params,1) ||
                    throw(ArgumentError(
                        "$caller: parametric_coord must hold " *
                        "$(size(record.node_params,1)) coordinates"))
                params=Float64[Float64(c) for c in parametric_coord]
                all(isfinite,params) || throw(ArgumentError(
                    "$caller: parametric coordinates must be finite"))
                record.node_params[:,position]=params
            end
            return nothing
        end
        cached=_cached_mesh_locked(caller)
        overlay=_high_order_overlay(cached)
        if overlay!==nothing && tag>nnodes(cached)
            (tag<=nnodes(overlay)) || throw(ArgumentError(
                "$caller: unknown node $tag; expected a dense tag in " *
                "1:$(nnodes(overlay))"))
            isempty(parametric_coord) || throw(ArgumentError(
                "$caller: parametric coordinates are recomputed from the " *
                "owning entity; pass an empty parametric_coord"))
            overlay.coords[:,tag]=values
            return nothing
        end
        count=nnodes(cached)
        (1<=tag<=count) || throw(ArgumentError(
            "$caller: unknown node $tag; expected a dense tag in 1:$count"))
        isempty(parametric_coord) || throw(ArgumentError(
            "$caller: parametric coordinates are recomputed from the owning " *
            "entity; pass an empty parametric_coord"))
        coordinates=copy(cached.coords)
        coordinates[:,tag]=values
        replacement=Mesh(coordinates;segs=cached.segs,tris=cached.tris,
                         tets=cached.tets,seg_tag=cached.seg_tag,
                         tri_tag=cached.tri_tag,tet_tag=cached.tet_tag)
        overlay=_high_order_overlay(cached)
        class=_cached_classification_locked(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.entities,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        overlay!==nothing &&
            _rebind_high_order!(replacement,new_class,caller)
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

# Remap caller-assigned record tags by an explicit (old, new) pair list.
# `update_nodes` also rewrites record element connectivity that references the
# old tag. New tags must stay positive and must not collide with a tag that is
# not itself being remapped.
function _record_renumber!(records,old_list,new_list,label,caller,
                           update_nodes::Bool,reserved::Int=0)
    isempty(old_list) && return nothing
    record_tags=Dict{Int32,DiscreteEntity}()
    for record in records
        source=label=="node" ? record.node_tags : record.element_tags
        for tag in source
            haskey(record_tags,tag) && throw(ArgumentError(
                "$caller: $label tag $tag is assigned on more than one " *
                "entity"))
            record_tags[tag]=record
        end
    end
    renames=Dict{Int32,Int32}()
    for (old_value,new_value) in zip(old_list,new_list)
        old=_mesh_query_integer(old_value,caller,"old_tags entry")
        new=_mesh_query_integer(new_value,caller,"new_tags entry")
        typemin(Int32)<=old<=typemax(Int32) || continue
        haskey(record_tags,Int32(old)) || continue
        (1<=new<=typemax(Int32)) || throw(ArgumentError(
            "$caller: new_tags entry $new must be a positive Int32 tag"))
        new<=reserved && throw(ArgumentError(
            "$caller: new tag $new collides with an existing $label tag"))
        renames[Int32(old)]=Int32(new)
    end
    for (old,new) in renames
        haskey(record_tags,new) && !haskey(renames,new) &&
            throw(ArgumentError(
                "$caller: new tag $new collides with an existing $label tag"))
    end
    for (old,new) in renames
        record=record_tags[old]
        if label=="node"
            position=findfirst(==(old),record.node_tags)
            record.node_tags[position]=new
            update_nodes || continue
            for connectivity in record.element_nodes
                for i in eachindex(connectivity)
                    connectivity[i]==old && (connectivity[i]=new)
                end
            end
        else
            position=findfirst(==(old),record.element_tags)
            record.element_tags[position]=new
        end
    end
    return nothing
end

function _renumber_nodes(old_tags=(),new_tags=())
    caller="API.mesh.renumber_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        records=[record for (_,_,record) in _discrete_mesh_records(model)]
        cached=LAST_MESH[]
        if cached===nothing
            isempty(records) && throw(ArgumentError(
                "$caller: no mesh; call API.mesh.generate first"))
            # Record-only session: an empty pair list compacts every node to
            # 1:n in entity order; explicit pairs remap sparse tags directly.
            if isempty(old_tags)
                mapping32=Dict{Int32,Int32}()
                next=Int32(1)
                for record in records, tag in record.node_tags
                    mapping32[tag]=next
                    next+=Int32(1)
                end
                for record in records
                    for i in eachindex(record.node_tags)
                        record.node_tags[i]=mapping32[record.node_tags[i]]
                    end
                    for connectivity in record.element_nodes
                        for i in eachindex(connectivity)
                            connectivity[i]=get(
                                mapping32,connectivity[i],connectivity[i])
                        end
                    end
                end
                return nothing
            end
            _record_renumber!(records,old_tags,new_tags,"node",caller,true)
            return nothing
        end
        record_tag_set=Set{Int32}()
        for record in records, tag in record.node_tags
            push!(record_tag_set,tag)
        end
        cache_old=Int[]
        cache_new=Int[]
        record_old=Int[]
        record_new=Int[]
        (old_tags isa AbstractVector || old_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: old_tags must be a vector or tuple of node tags"))
        (new_tags isa AbstractVector || new_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: new_tags must be a vector or tuple of node tags"))
        length(old_tags)==length(new_tags) || throw(ArgumentError(
            "$caller: old_tags and new_tags must have the same length"))
        for (old_value,new_value) in zip(old_tags,new_tags)
            old=_mesh_query_integer(old_value,caller,"old_tags entry")
            new=_mesh_query_integer(new_value,caller,"new_tags entry")
            if Int32(0)<=(old<=typemax(Int32) ? Int32(old) : Int32(0)) &&
               Int32(old) in record_tag_set
                push!(record_old,old)
                push!(record_new,new)
            else
                push!(cache_old,old)
                push!(cache_new,new)
            end
        end
        _record_renumber!(records,record_old,record_new,"node",caller,true,nnodes(cached))
        if isempty(cache_old)
            # Empty-pair form: dense cache tags are already 1:n, so only
            # record nodes compact — to count+1.. in entity order, with
            # record connectivity rewritten through the same mapping.
            isempty(old_tags) || return nothing
            mapping32=Dict{Int32,Int32}()
            next=Int32(nnodes(cached)+1)
            for record in records, tag in record.node_tags
                mapping32[tag]=next
                next+=Int32(1)
            end
            for record in records
                for i in eachindex(record.node_tags)
                    record.node_tags[i]=mapping32[record.node_tags[i]]
                end
                for connectivity in record.element_nodes
                    for i in eachindex(connectivity)
                        connectivity[i]=get(
                            mapping32,connectivity[i],connectivity[i])
                    end
                end
            end
            return nothing
        end
        mapping=_mesh_renumber_permutation(
            cache_old,cache_new,nnodes(cached),caller,"node")
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
        overlay=_high_order_overlay(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            new_mesh,class.entity,class.entities,node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(new_mesh,new_class)
        overlay!==nothing &&
            _rebind_high_order!(new_mesh,new_class,caller)
        return nothing
    end
end

# Element tags are dense across the segment/triangle/tetrahedron blocks, so a
# renumbering must keep every tag inside its own block range — cross-block
# assignments would change element types and are rejected explicitly.
function _renumber_elements(old_tags=(),new_tags=())
    caller="API.mesh.renumber_elements"
    return lock(STATE_LOCK) do
        model=_model_locked()
        records=[record for (_,_,record) in _discrete_mesh_records(model)]
        cached=LAST_MESH[]
        if cached===nothing
            isempty(records) && throw(ArgumentError(
                "$caller: no mesh; call API.mesh.generate first"))
            if isempty(old_tags)
                next=Int32(1)
                for record in records
                    for i in eachindex(record.element_tags)
                        record.element_tags[i]=next
                        next+=Int32(1)
                    end
                end
                return nothing
            end
            _record_renumber!(records,old_tags,new_tags,"element",caller,false)
            return nothing
        end
        triangle_offset,tetrahedron_offset,total=_mesh_element_offsets(cached)
        record_tag_set=Set{Int32}()
        for record in records, etag in record.element_tags
            push!(record_tag_set,etag)
        end
        (old_tags isa AbstractVector || old_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: old_tags must be a vector or tuple of element tags"))
        (new_tags isa AbstractVector || new_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: new_tags must be a vector or tuple of element tags"))
        length(old_tags)==length(new_tags) || throw(ArgumentError(
            "$caller: old_tags and new_tags must have the same length"))
        cache_old=Int[]
        cache_new=Int[]
        record_old=Int[]
        record_new=Int[]
        for (old_value,new_value) in zip(old_tags,new_tags)
            old=_mesh_query_integer(old_value,caller,"old_tags entry")
            new=_mesh_query_integer(new_value,caller,"new_tags entry")
            if 0<old<=typemax(Int32) && Int32(old) in record_tag_set
                push!(record_old,old)
                push!(record_new,new)
            else
                push!(cache_old,old)
                push!(cache_new,new)
            end
        end
        _record_renumber!(records,record_old,record_new,"element",caller,false,total)
        if isempty(cache_old)
            # With records present the empty-pair form compacts record
            # element tags to continue the dense sequence in entity order.
            isempty(old_tags) || return nothing
            next=Int32(total+1)
            for record in records
                for i in eachindex(record.element_tags)
                    record.element_tags[i]=next
                    next+=Int32(1)
                end
            end
            return nothing
        end
        mapping=_mesh_renumber_permutation(
            cache_old,cache_new,total,caller,"element")
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
        overlay=_high_order_overlay(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            replacement,class.entity,class.entities,class.node_entities,class.boundaries,
            owners_out[1],owners_out[2],owners_out[3])
        _replace_mesh_cache_locked!(replacement,new_class)
        overlay!==nothing &&
            _rebind_high_order!(replacement,new_class,caller)
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
        msh=_mesh_query_integer(element_type,caller,"element_type")
        spec=msh_spec(msh)
        dimension=spec.dim
        entity=_mesh_query_integer(tag,caller,"tag")
        record=_model_mesh_record(model,dimension,entity)
        if record!==nothing && haskey(model.discrete,(dimension,entity))
            # Discrete/attached record elements reorder in place: ordering
            # holds 0-based source positions of that type's elements.
            indices=Int[]
            for i in eachindex(record.element_tags)
                record.element_types[i]==Int32(msh) && push!(indices,i)
            end
            (ordering isa AbstractVector || ordering isa Tuple) ||
                throw(ArgumentError(
                    "$caller: ordering must be a vector or tuple of " *
                    "0-based source positions"))
            isempty(indices) && throw(ArgumentError(
                "$caller: no elements of type $msh on " *
                "$(_MESH_ENTITY_LABELS[dimension+1])[$entity] to reorder"))
            length(ordering)==length(indices) || throw(ArgumentError(
                "$caller: ordering must hold $(length(indices)) entries " *
                    "for the elements of type $msh on " *
                    "$(_MESH_ENTITY_LABELS[dimension+1])[$entity]"))
            permutation=Int[]
            for value in ordering
                source=_mesh_query_integer(value,caller,"ordering entry")
                (0<=source<length(indices)) || throw(ArgumentError(
                    "$caller: ordering entry $source is out of range; " *
                    "expected 0:$(length(indices)-1)"))
                source in permutation && throw(ArgumentError(
                    "$caller: ordering repeats position $source"))
                push!(permutation,source)
            end
            chosen=indices[permutation .+ 1]
            tags_snapshot=record.element_tags[chosen]
            types_snapshot=record.element_types[chosen]
            nodes_snapshot=record.element_nodes[chosen]
            for (slot,position) in enumerate(indices)
                record.element_tags[position]=tags_snapshot[slot]
                record.element_types[position]=types_snapshot[slot]
                record.element_nodes[position]=nodes_snapshot[slot]
            end
            return nothing
        end
        cached=_cached_mesh_locked(caller)
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
        overlay=_high_order_overlay(cached)
        new_class=_MeshClassification(replacement,class.entity,class.entities,
            class.node_entities,class.boundaries,class.seg_entities,
            class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(replacement,new_class)
        overlay!==nothing &&
            _rebind_high_order!(replacement,new_class,caller)
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

# Highest caller-assigned tag across every discrete record and attached store.
function _records_max_tag(m::GeoModel,field::Symbol)
    result=0
    for (_,_,record) in _discrete_mesh_records(m)
        values=getfield(record,field)
        isempty(values) || (result=max(result,Int(maximum(values))))
    end
    return result
end

function _import_stl()
    caller="API.mesh.import_stl"
    return lock(STATE_LOCK) do
        m=_model_locked()
        file=MODEL_FILE_NAME[]
        (isempty(file) || !isfile(file)) && return nothing
        lowercase(last(splitext(file)))==".stl" || return nothing
        stl=read_stl(file)
        ntris(stl)==0 && return nothing
        tag=1
        while _model_entity_known(m,2,tag)
            tag+=1
        end
        add_discrete_entity!(m,2,tag)
        cached=LAST_MESH[]
        node_base=cached===nothing ? 0 : nnodes(cached)
        node_base=max(node_base,_records_max_tag(m,:node_tags))
        element_base=0
        if cached!==nothing
            _,_,element_base=_mesh_element_offsets(cached)
        end
        element_base=max(element_base,_records_max_tag(m,:element_tags))
        add_discrete_nodes!(m,2,tag,
                            collect(node_base+1:node_base+nnodes(stl)),
                            vec(stl.coords))
        add_discrete_elements!(m,2,tag,[2],
                               [collect(element_base+1:
                                        element_base+ntris(stl))],
                               [vec(stl.tris) .+ Int32(node_base)])
        return nothing
    end
end

function _create_topology(make_simply_connected=true,export_discrete=true)
    caller="API.mesh.create_topology"
    flag1=_mesh_query_bool(make_simply_connected,caller,"make_simply_connected")
    flag2=_mesh_query_bool(export_discrete,caller,"export_discrete")
    return lock(STATE_LOCK) do
        create_topology!(_model_locked();
                         make_simply_connected=flag1,export_discrete=flag2)
        return nothing
    end
end

function _classify_surfaces(angle=40*pi/180,boundary=true,
                            for_reparametrization=false,curve_angle=pi)
    caller="API.mesh.classify_surfaces"
    for (value,name) in ((angle,"angle"),(curve_angle,"curve_angle"))
        (value isa Real && !(value isa Bool) && isfinite(value)) ||
            throw(ArgumentError("$caller: $name must be a finite real number"))
    end
    flag1=_mesh_query_bool(boundary,caller,"boundary")
    flag2=_mesh_query_bool(for_reparametrization,caller,"for_reparametrization")
    return lock(STATE_LOCK) do
        classify_surfaces!(_model_locked();angle=Float64(angle),
                         boundary=flag1,for_reparametrization=flag2,
                         curve_angle=Float64(curve_angle))
        return nothing
    end
end

function _create_geometry(dim_tags=Tuple{Int,Int}[])
    caller="API.mesh.create_geometry"
    normalized=_mesh_parse_dim_tags(dim_tags,caller)
    return lock(STATE_LOCK) do
        create_geometry!(_model_locked(),
                         isempty(normalized) ? nothing : normalized)
        return nothing
    end
end

function _embed(dim,tags,in_dim,in_tag)
    caller="API.mesh.embed"
    dimension=_mesh_query_integer(dim,caller,"dim")
    target_dim=_mesh_query_integer(in_dim,caller,"inDim")
    target=_mesh_query_integer(in_tag,caller,"inTag")
    (tags isa AbstractVector || tags isa Tuple) || throw(ArgumentError(
        "$caller: tags must be a vector of entity tags"))
    entity_tags=[_mesh_query_integer(entry,caller,"tags") for entry in tags]
    return _with_model(invalidate=true) do current
        embed!(current,dimension,entity_tags,target_dim,target)
    end
end

# `optimize` maps Gmsh's default tetrahedral mesh optimizer onto the validated
# boundary-preserving `smooth_optimize` kernel. Connectivity and the
# index-aligned classification snapshot carry over unchanged because only node
# coordinates move and boundary nodes are fixed. A `dimTags` selection freezes
# every node incident to unselected elements; "Laplace2D"/"Relocate2D" run
# isotropic Laplacian smoothing on a triangle cache. Optimizer names with no
# implemented kernel fail explicitly.
function _optimize_mesh(method="",force=false,niter=1,dim_tags=())
    caller="API.mesh.optimize"
    method isa AbstractString || throw(ArgumentError(
        "$caller: method must be a string"))
    # Method names map onto the implemented kernels: the default quality
    # optimizer is the tetrahedral relocation smoother (Gmsh's "Gmsh" and
    # "Relocate3D" roles), and the 2-D methods run Laplacian smoothing on a
    # triangle cache. Netgen and the high-order optimizers need machinery
    # this build does not have.
    method in ("","Gmsh","Relocate3D","Laplace2D","Relocate2D") ||
        throw(ArgumentError(
            "$caller: unknown or unsupported optimizer \"$method\""))
    force isa Bool || throw(ArgumentError("$caller: force must be Bool"))
    iterations=_mesh_query_integer(niter,caller,"niter")
    iterations>=0 || throw(ArgumentError("$caller: niter must be nonnegative"))
    pairs=_mesh_parse_dim_tags(dim_tags,caller)
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=_cached_mesh_locked(caller)
        for (edim,etag) in pairs
            _model_entity_known(model,edim,etag) || throw(ArgumentError(
                "$caller: unknown " *
                "$(_MESH_ENTITY_LABELS[edim+1])[$etag]"))
        end
        class=_cached_classification_locked(cached)
        movable=_optimize_movable_mask(cached,class,pairs)
        smoothed=if method in ("Laplace2D","Relocate2D")
            _laplacian_smooth_tri_cache(cached,iterations,movable)
        else
            # `reverse` may legitimately leave inverted tets; the smoother's
            # per-star positivity guard already rejects every move there, so
            # bypass only the entry orientation check (structure is still
            # validated) — matching Gmsh, which optimizes without error.
            smooth_optimize(cached;iters=iterations,
                            require_positive_tets=false,movable=movable)
        end
        overlay=_high_order_overlay(cached)
        new_class=class===nothing ? nothing : _MeshClassification(
            smoothed,class.entity,class.entities,class.node_entities,class.boundaries,
            class.seg_entities,class.tri_entities,class.tet_entities)
        _replace_mesh_cache_locked!(smoothed,new_class)
        overlay!==nothing &&
            _rebind_high_order!(smoothed,new_class,caller)
        return nothing
    end
end

# A node is movable under an entity-scoped `optimize` only when every top-dim
# cell incident to it belongs to the selection, so moves never deform
# unselected elements. An empty selection leaves every node movable.
function _optimize_movable_mask(cached::Mesh,class,pairs)
    isempty(pairs) && return nothing
    selected=Set{Tuple{Int,Int32}}(
        (edim,Int32(etag)) for (edim,etag) in pairs)
    inside=falses(nnodes(cached));outside=falses(nnodes(cached))
    if class!==nothing
        for (cells,owners,edim) in
                ((cached.segs,class.seg_entities,1),
                 (cached.tris,class.tri_entities,2),
                 (cached.tets,class.tet_entities,3))
            for cell in axes(cells,2)
                member=(edim,owners[cell]) in selected
                for slot in axes(cells,1)
                    node=cells[slot,cell]
                    member ? (inside[node]=true) : (outside[node]=true)
                end
            end
        end
    end
    return inside .& .!outside
end

# Laplacian smoothing of a triangle cache — the "Laplace2D"/"Relocate2D"
# method roles. Boundary nodes and nodes outside the entity selection stay
# fixed; moves are the isotropic neighbor average in all three coordinates.
function _laplacian_smooth_tri_cache(mesh::Mesh,iterations::Int,
                                     movable)
    boundary=falses(nnodes(mesh))
    @inbounds for edge in first(boundary_edges(mesh.tris)),node in edge
        boundary[node]=true
    end
    neighbors=[Int32[] for _ in 1:nnodes(mesh)]
    @inbounds for cell in axes(mesh.tris,2),e in ((1,2),(2,3),(3,1))
        a=mesh.tris[e[1],cell];b=mesh.tris[e[2],cell]
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
            movable===nothing || movable[node] || continue
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

# The mesh-level variant stores per-element flags for one window. It is display
# state only.
function _set_mesh_visibility_per_window(tag,value,window_index)
    caller="API.mesh.set_visibility_per_window"
    element_tag=_mesh_query_integer(tag,caller,"tag")
    flag=_mesh_query_integer(value,caller,"value")
    typemin(Int32)<=flag<=typemax(Int32) || throw(ArgumentError(
        "$caller: value $value is out of range"))
    window=_mesh_query_integer(window_index,caller,"window_index")
    window>=0 || throw(ArgumentError(
        "$caller: window_index must be non-negative"))
    return lock(STATE_LOCK) do
        _model_locked()
        cached=LAST_MESH[]
        total=cached===nothing ? 0 : _mesh_element_offsets(cached)[3]
        element_tag in 1:total || return nothing
        index=_current_slot_index_locked()
        index==0 && throw(ArgumentError(
            "$caller: the current model is not in the model list"))
        per_window=get!(Dict{Int,Int32},
                        MODEL_SLOTS[index].window_element_visibility,window)
        per_window[element_tag]=Int32(flag)
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
# recombine, smoothing, reverse, algorithm, size channels, compounds, outward
# orientation) for the listed entities — all of them when `dim_tags` is empty;
# periodic relations, embeddings, and Point mesh sizes are retained.
function _remove_constraints(dim_tags=())
    caller="API.mesh.remove_constraints"
    pairs=_mesh_parse_dim_tags(dim_tags,caller)
    return _with_model() do current
        remove_constraints!(current,pairs)
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

# --- Generation-attribute and discrete-entity internals ---------------------
# Thin locked forwards onto the model-layer records in
# geometry/ModelMeshingAttributes.jl; generation consumes them in `_generate`.

function _set_transfinite_curve(tag,num_nodes,mesh_type="Progression",coef=1.0)
    lock(STATE_LOCK) do
        set_transfinite_curve!(_model_locked(),tag,num_nodes,mesh_type,coef)
    end
    return nothing
end

function _set_transfinite_surface(tag,arrangement="Left",corner_tags=Int[])
    lock(STATE_LOCK) do
        set_transfinite_surface!(_model_locked(),tag,arrangement,corner_tags)
    end
    return nothing
end

function _set_transfinite_volume(tag,corner_tags=Int[])
    lock(STATE_LOCK) do
        set_transfinite_volume!(_model_locked(),tag,corner_tags)
    end
    return nothing
end

function _set_transfinite_automatic(dim_tags=(),corner_angle=2.35,
                                    recombine=true)
    lock(STATE_LOCK) do
        set_transfinite_automatic!(_model_locked(),dim_tags,corner_angle,
                                   recombine)
    end
    return nothing
end

function _set_recombine_attribute(dim,tag,angle=45.0)
    lock(STATE_LOCK) do
        set_recombine!(_model_locked(),dim,tag,angle)
    end
    return nothing
end

function _set_smoothing(dim,tag,val)
    lock(STATE_LOCK) do
        set_smoothing!(_model_locked(),dim,tag,val)
    end
    return nothing
end

function _set_reverse_attribute(dim,tag,val=true)
    lock(STATE_LOCK) do
        set_reverse!(_model_locked(),dim,tag,val)
    end
    return nothing
end

function _set_algorithm(dim,tag,val)
    lock(STATE_LOCK) do
        set_algorithm!(_model_locked(),dim,tag,val)
    end
    return nothing
end

function _set_size_at_parametric_points(dim,tag,parametric_coord,sizes)
    lock(STATE_LOCK) do
        set_size_at_parametric_points!(_model_locked(),dim,tag,
                                       parametric_coord,sizes)
    end
    return nothing
end

function _set_size_from_boundary(dim,tag,val)
    lock(STATE_LOCK) do
        set_size_from_boundary!(_model_locked(),dim,tag,val)
    end
    return nothing
end

function _set_size_callback(callback)
    lock(STATE_LOCK) do
        set_size_callback!(_model_locked(),callback)
    end
    return nothing
end

function _remove_size_callback()
    lock(STATE_LOCK) do
        set_size_callback!(_model_locked(),nothing)
    end
    return nothing
end

function _set_compound(dim,tags)
    lock(STATE_LOCK) do
        set_compound!(_model_locked(),dim,tags)
    end
    return nothing
end

function _set_outward_orientation(tag)
    lock(STATE_LOCK) do
        set_outward_orientation!(_model_locked(),tag)
    end
    return nothing
end

# `setOrder` acts immediately on the current mesh in Gmsh 4.15.2 and is also
# stored as a generator attribute consumed by the next `generate`.
function _set_order(order)
    lock(STATE_LOCK) do
        model=_model_locked()
        set_order!(model,order)
        cached=LAST_MESH[]
        cached===nothing && return nothing
        _apply_mesh_order(model,cached,
                          _cached_classification_locked(cached),
                          "API.mesh.set_order")
    end
    return nothing
end

function _add_nodes(dim,tag,node_tags,coord,parametric_coord=Float64[])
    lock(STATE_LOCK) do
        add_discrete_nodes!(_model_locked(),dim,tag,node_tags,coord,
                            parametric_coord)
    end
    return nothing
end

function _add_elements(dim,tag,element_types,element_tags,node_tags)
    lock(STATE_LOCK) do
        add_discrete_elements!(_model_locked(),dim,tag,element_types,
                               element_tags,node_tags)
    end
    return nothing
end

function _add_elements_by_type(tag,element_type,element_tags,node_tags)
    caller="API.mesh.add_elements_by_type"
    mtype=_mesh_query_integer(element_type,caller,"element_type")
    # Gmsh derives the entity dimension from the element type, then looks the
    # entity up by tag ("Surface 1 does not exist" for a type-2 element).
    dim=Int(msh_spec(mtype).dim)
    return lock(STATE_LOCK) do
        current=_model_locked()
        _model_entity_known(current,dim,tag) || throw(ArgumentError(
            "$caller: $(_MESH_ENTITY_LABELS[dim+1]) $tag does not exist"))
        add_discrete_elements!(current,dim,tag,[mtype],[element_tags],
                               [node_tags])
    end
    return nothing
end

function _add_homology_request(kind="Homology",domain_tags=Int[],
                               subdomain_tags=Int[],dims=Int[])
    lock(STATE_LOCK) do
        add_homology_request!(_model_locked();kind=kind,
                              domain_tags=domain_tags,
                              subdomain_tags=subdomain_tags,dims=dims)
    end
    return nothing
end

function _clear_homology_requests()
    lock(STATE_LOCK) do
        clear_homology_requests!(_model_locked())
    end
    return nothing
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
    LAST_MESH_PARTITION[]=nothing
    return nothing
end

# Deterministic element partitioner over the cache dual graph (two elements
# are adjacent when they share a node). Each partition grows breadth-first
# from the lowest-tag unassigned element until it reaches its quota; when the
# frontier disconnects, a fresh seed continues the same partition so every
# element is always assigned.
function _partition_elements(mesh::Mesh,num_partitions::Int)
    total=nsegs(mesh)+ntris(mesh)+ntets(mesh)
    incidence=[Int[] for _ in 1:nnodes(mesh)]
    tag=0
    for cells in (mesh.segs,mesh.tris,mesh.tets)
        for column in axes(cells,2)
            tag+=1
            for row in axes(cells,1)
                push!(incidence[Int(cells[row,column])],tag)
            end
        end
    end
    assignment=zeros(Int32,total)
    quota=cld(total,num_partitions)
    part=0
    assigned=0
    seed=1
    while assigned<total
        part+=1
        part>num_partitions && (part=num_partitions)
        target=part==num_partitions ? total : min(assigned+quota,total)
        queue=Int[]
        while assigned<target
            while seed<=total && assignment[seed]!=0
                seed+=1
            end
            seed>total && break
            push!(queue,seed);assignment[seed]=Int32(part);assigned+=1
            head=1
            while head<=length(queue) && assigned<target
                element=queue[head];head+=1
                cells=element<=nsegs(mesh) ? mesh.segs[:,element] :
                      element<=nsegs(mesh)+ntris(mesh) ?
                          mesh.tris[:,element-nsegs(mesh)] :
                          mesh.tets[:,element-nsegs(mesh)-ntris(mesh)]
                neighbors=Int[]
                for node in cells
                    append!(neighbors,incidence[Int(node)])
                end
                sort!(unique!(neighbors))
                for neighbor in neighbors
                    assignment[neighbor]==0 || continue
                    assignment[neighbor]=Int32(part);assigned+=1
                    push!(queue,neighbor)
                    assigned<target || break
                end
            end
        end
    end
    return assignment
end

function _partition(num_part,element_tags=(),partitions=())
    caller="API.mesh.partition"
    count=_mesh_query_integer(num_part,caller,"num_part")
    count>=0 || throw(ArgumentError("$caller: num_part must be non-negative"))
    (element_tags isa AbstractVector || element_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: element_tags must be a vector of element tags"))
    (partitions isa AbstractVector || partitions isa Tuple) ||
        throw(ArgumentError(
            "$caller: partitions must be a vector of partition tags"))
    length(element_tags)==length(partitions) || throw(ArgumentError(
        "$caller: element_tags and partitions must have matching lengths"))
    return lock(STATE_LOCK) do
        _model_locked()
        cached=_cached_mesh_locked(caller)
        _,_,total=_mesh_element_offsets(cached)
        total>0 || throw(ArgumentError("$caller: the mesh has no elements"))
        if isempty(element_tags)
            count>0 || throw(ArgumentError(
                "$caller: num_part must be positive without an explicit " *
                "element partition list"))
            assignment=_partition_elements(cached,count)
        else
            assignment=zeros(Int32,total)
            supplied=Dict{Int,Int}()
            for (raw_tag,raw_part) in zip(element_tags,partitions)
                tag=_mesh_query_integer(raw_tag,caller,"element_tags")
                1<=tag<=total || throw(ArgumentError(
                    "$caller: element tag $tag is outside 1:$total"))
                part=_mesh_query_integer(raw_part,caller,"partitions")
                1<=part || throw(ArgumentError(
                    "$caller: partition $part must be positive"))
                haskey(supplied,tag) && throw(ArgumentError(
                    "$caller: element tag $tag is listed more than once"))
                supplied[tag]=part
            end
            length(supplied)==total || throw(ArgumentError(
                "$caller: all elements are not partitioned"))
            for (tag,part) in supplied
                assignment[tag]=Int32(part)
            end
            effective=count>0 ? count : maximum(assignment)
            any(>(effective),assignment) && throw(ArgumentError(
                "$caller: partition tags exceed num_part $effective"))
            count=effective
        end
        LAST_MESH_PARTITION[]=_MeshPartition(
            cached,Int(count),assignment)
        return nothing
    end
end

function _mesh_entity_element_tags(class,model,dimension::Int,entity::Int,
                                   cached::Mesh)
    result=Int[]
    tri_offset,tet_offset,_=_mesh_element_offsets(cached)
    for (owner_dim,cells,owners,offset) in
            ((1,cached.segs,class.seg_entities,0),
             (2,cached.tris,class.tri_entities,tri_offset),
             (3,cached.tets,class.tet_entities,tet_offset))
        owner_dim==dimension || continue
        for (column,owner) in enumerate(owners)
            owner==Int32(entity) && push!(result,offset+column)
        end
    end
    return result
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
# Per-entity element list across the simplex cache (dim-matched blocks) and
# the entity's discrete/attached record — the input cells for homology.
function _entity_all_elements(m::GeoModel,cached,class,dim::Int,tag::Int)
    out=Tuple{Int32,Vector{Int32}}[]
    if cached!==nothing && class!==nothing
        _,triangle_offset,tetrahedron_offset=_mesh_element_offsets(cached)
        blocks=dim==1 ? ((cached.segs,class.seg_entities,Int32(1)),) :
               dim==2 ? ((cached.tris,class.tri_entities,Int32(2)),) :
               dim==3 ? ((cached.tets,class.tet_entities,Int32(4)),) : ()
        for (cells,owners,msh_type) in blocks
            for column in axes(cells,2)
                owners[column]==tag || continue
                push!(out,(msh_type,Int32.(cells[:,column])))
            end
        end
    end
    record=_entity_record(m,dim,tag)
    if record!==nothing
        for index in eachindex(record.element_types)
            push!(out,(record.element_types[index],
                       copy(record.element_nodes[index])))
        end
    end
    return out
end

function _compute_homology()
    caller="API.mesh.compute_homology"
    return lock(STATE_LOCK) do
        m=_model_locked()
        isempty(m.meshing.homology_requests) && return Tuple{Int,Int}[]
        cached=LAST_MESH[]
        class=cached===nothing ? nothing :
              _cached_classification_locked(cached)
        cells=Dict{Tuple{Int,Int},Vector{Tuple{Int32,Vector{Int32}}}}()
        if class!==nothing
            for dim in 1:3
                owners=dim==1 ? class.seg_entities :
                       dim==2 ? class.tri_entities : class.tet_entities
                for tag in unique(owners)
                    key=(dim,Int(tag))
                    entries=_entity_all_elements(
                        m,cached,class,dim,Int(tag))
                    isempty(entries) || (cells[key]=entries)
                end
            end
            # Dimension-0 entities own nodes, not elements; their cells come
            # from the closure of incident higher-dimensional cells.
        end
        for (edim,etag,record) in _discrete_mesh_records(m)
            key=(edim,etag)
            entries=get(cells,key,Tuple{Int32,Vector{Int32}}[])
            for index in eachindex(record.element_types)
                push!(entries,(record.element_types[index],
                               copy(record.element_nodes[index])))
            end
            isempty(entries) || (cells[key]=entries)
        end
        isempty(cells) && return Tuple{Int,Int}[]
        return Tuple{Int,Int}[pair for pair in
            compute_homology!(m,cells;
                              requests=m.meshing.homology_requests)]
    end
end

# `getPeriodicKeys` returns one key pair per element-dof incidence: the slave
# node tag, the master node tag it is periodic with, the function-space type
# index of the dof (0 for nodal spaces), and the incidence coordinates.
function _get_periodic_keys(element_type,function_space_type,tag,
                            return_coord=true)
    caller="API.mesh.get_periodic_keys"
    msh_type=_mesh_query_integer(element_type,caller,"elementType")
    spec=try msh_spec(msh_type) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: unknown elementType $msh_type"))
    end
    function_space_type isa AbstractString || throw(ArgumentError(
        "$caller: functionSpaceType must be a string"))
    entity_tag=_mesh_query_integer(tag,caller,"tag")
    flag=_mesh_query_bool(return_coord,caller,"return_coord")
    dimension=spec.dim
    return lock(STATE_LOCK) do
        m=_model_locked()
        entity_tag<=typemax(Int32) || throw(ArgumentError(
            "$caller: tag is out of range"))
        _model_entity_known(m,dimension,Int(entity_tag)) || throw(ArgumentError(
            "$caller: unknown entity ($dimension,$entity_tag)"))
        constraint=get(m.periodic,(dimension,Int(entity_tag)),nothing)
        cached=LAST_MESH[]
        class=cached===nothing ? nothing :
              _cached_classification_locked(cached)
        master_tag=constraint===nothing ? Int32(entity_tag) :
                    Int32(constraint.master_entity)
        pairing=Dict{Int32,Int32}()
        if constraint!==nothing && cached!==nothing && class!==nothing
            mapping=model_periodic_nodes(m,cached,dimension,Int(entity_tag))
            for (slave_node,master_node) in zip(mapping.slave_nodes,
                                                mapping.master_nodes)
                pairing[Int32(slave_node)]=Int32(master_node)
            end
        end
        incidences=Vector{Int32}[]
        if cached!==nothing && class!==nothing
            cells=_mesh_element_block(cached,msh_type)
            if cells!==nothing
                owners=msh_type==1 ? class.seg_entities :
                       msh_type==2 ? class.tri_entities :
                       msh_type==4 ? class.tet_entities : nothing
                if owners!==nothing
                    for column in axes(cells,2)
                        owners[column]==entity_tag || continue
                        push!(incidences,Int32.(cells[:,column]))
                    end
                end
            end
        end
        record=_entity_record(m,dimension,Int(entity_tag))
        if record!==nothing
            for index in eachindex(record.element_types)
                record.element_types[index]==msh_type || continue
                push!(incidences,copy(record.element_nodes[index]))
            end
        end
        if constraint!==nothing && !isempty(incidences)
            master_record=_entity_record(m,dimension,Int(master_tag))
            master_nodes=Dict{NTuple{3,Float64},Int32}()
            if master_record!==nothing
                for column in eachindex(master_record.node_tags)
                    master_nodes[
                        (master_record.node_coords[1,column],
                         master_record.node_coords[2,column],
                         master_record.node_coords[3,column])]=
                        master_record.node_tags[column]
                end
            elseif cached!==nothing && class!==nothing
                for node in eachindex(class.node_entities)
                    class.node_entities[node]!=(dimension,Int32(master_tag)) &&
                        continue
                    master_nodes[(cached.coords[1,node],
                                  cached.coords[2,node],
                                  cached.coords[3,node])]=Int32(node)
                end
            end
            coefficients,translation,_=_transform_homogeneous(
                constraint.affine,caller;name="periodic affine transform")
            for nodes in incidences
                for node in nodes
                    haskey(pairing,node) && continue
                    value=_mesh_node_coords(m,cached,node)
                    value===nothing && continue
                    mapped=_model_affine_point(
                        coefficients,translation,value,caller,Int(node))
                    best=Int32(0);best_distance=Inf
                    for (mc,mtag) in master_nodes
                        distance=sum(k->(mc[k]-mapped[k])^2,1:3)
                        if distance<best_distance
                            best_distance=distance;best=mtag
                        end
                    end
                    paired=best_distance<=max(constraint.atol^2,1e-18) ?
                           best : nothing
                    pairing[node]=paired===nothing ? node : paired
                end
            end
        end
        total=sum(length(nodes) for nodes in incidences)
        entity_keys=Vector{UInt64}(undef,total)
        master_keys=Vector{UInt64}(undef,total)
        coordinates=Float64[]
        flag && sizehint!(coordinates,3*total)
        cursor=0
        for nodes in incidences
            for node in nodes
                cursor+=1
                entity_keys[cursor]=UInt64(node)
                master_keys[cursor]=UInt64(get(pairing,node,node))
                if flag
                    value=_mesh_node_coords(m,cached,node)
                    value===nothing && throw(ArgumentError(
                        "$caller: node $node has no coordinates"))
                    append!(coordinates,value)
                end
            end
        end
        type_keys=zeros(Int32,total)
        return Int32(master_tag),type_keys,copy(type_keys),
               entity_keys,master_keys,coordinates
    end
end

# ---- recombine / split_quadrangles --------------------------------------
#
# Gmsh's `mesh.recombine()` recombines surface triangles into quadrangles on
# every entity flagged through `mesh.set_recombine`; the stored angle is the
# recombination quality angle. `recombine_triangles` produces a `MixedMesh`;
# resulting quadrangles and leftover triangles are stored on the entity's
# discrete record (created on demand for native entities) while the source
# triangles are dropped from the simplex cache or record they came from.

# Collect the (dim=2) triangle elements of `key` across the cache and the
# entity's record. Returns `(tris, sources)` where `tris` are
# `(element_tag, node_tags)` pairs and `sources` identifies where each entry
# came from (`(:cache, column)` or `(:record, index)`).
function _entity_record(m::GeoModel,dim::Int,tag::Int)
    haskey(m.discrete,(dim,tag)) && return m.discrete[(dim,tag)]
    return get(m.meshing.attached,(dim,tag),nothing)
end

function _entity_tri_elements(m::GeoModel,cached,class,dim::Int,tag::Int)
    entries=Tuple{Int32,Vector{Int32}}[]
    cache_cols=Int[]
    if cached!==nothing && class!==nothing
        triangle_offset,_,_=_mesh_element_offsets(cached)
        for column in eachindex(class.tri_entities)
            class.tri_entities[column]==tag || continue
            push!(cache_cols,column)
            push!(entries,(Int32(triangle_offset+column),
                           Int32.(cached.tris[:,column])))
        end
    end
    record=_entity_record(m,dim,tag)
    if record!==nothing
        for index in eachindex(record.element_types)
            record.element_types[index]==2 || continue
            push!(entries,(record.element_tags[index],
                           copy(record.element_nodes[index])))
        end
    end
    return entries,cache_cols,record
end

# Global node coordinate lookup for tag-addressed connectivity: cache nodes
# are dense 1:nnodes, record nodes resolve through the global record index.
function _mesh_node_coords(m::GeoModel,cached,tag::Int32)
    if cached!==nothing && 1<=tag<=size(cached.coords,2)
        return (cached.coords[1,tag],cached.coords[2,tag],cached.coords[3,tag])
    end
    index=_record_node_index(m)
    haskey(index,tag) || return nothing
    record,column=index[tag]
    return (record.node_coords[1,column],record.node_coords[2,column],
            record.node_coords[3,column])
end

function _recombine()
    caller="API.mesh.recombine"
    return lock(STATE_LOCK) do
        m=_model_locked()
        settings=sort!(collect(m.meshing.recombine);by=first)
        isempty(settings) && return nothing
        cached=LAST_MESH[]
        class=cached===nothing ? nothing : _cached_classification_locked(cached)
        for ((dim,tag),angle) in settings
            dim==2 || continue
            _model_entity_known(m,dim,tag) || throw(ArgumentError(
                "$caller: unknown entity ($dim, $tag)"))
            entries,cache_cols,record=_entity_tri_elements(
                m,cached,class,dim,tag)
            length(entries)<2 && continue
            sub_tags=Int32[]
            position=Dict{Int32,Int32}()
            for (_,nodes) in entries
                for node in nodes
                    haskey(position,node) && continue
                    position[node]=Int32(length(sub_tags)+1)
                    push!(sub_tags,node)
                end
            end
            coords=Matrix{Float64}(undef,3,length(sub_tags))
            for (i,node) in enumerate(sub_tags)
                value=_mesh_node_coords(m,cached,node)
                value===nothing && throw(ArgumentError(
                    "$caller: node $node has no coordinates"))
                coords[:,i].=value
            end
            cells=Matrix{Int32}(undef,3,length(entries))
            for (column,(_,nodes)) in enumerate(entries)
                length(nodes)==3 || throw(ArgumentError(
                    "$caller: non-triangle element in recombine scope"))
                cells[:,column].=[position[n] for n in nodes]
            end
            submesh=Mesh(coords;tris=cells)
            recombined=recombine_triangles(submesh)
            target=record
            if target===nothing
                add_discrete_entity!(m,dim,tag)
                target=m.discrete[(dim,tag)]
            end
            extra=size(recombined.coords,2)-length(sub_tags)
            extra<0 && throw(ArgumentError(
                "$caller: recombine dropped nodes unexpectedly"))
            fresh_node_tag=_model_fresh_node_tag(m,cached)
            for i in 1:extra
                column=length(sub_tags)+i
                _record_append_node!(target,fresh_node_tag,
                    [recombined.coords[1,column],recombined.coords[2,column],
                     recombined.coords[3,column]],Float64[])
                sub_tag=fresh_node_tag
                push!(sub_tags,sub_tag)
                fresh_node_tag+=Int32(1)
            end
            element_tag=_model_fresh_element_tag(m)
            for block in recombined.blocks
                block isa ElementBlock || continue
                for column in axes(block.nodes,2)
                    nodes=Int32[sub_tags[block.nodes[row,column]]
                                for row in axes(block.nodes,1)]
                    _record_append_element!(target,block.msh,
                        Int32(element_tag),nodes)
                    element_tag+=1
                end
            end
            isempty(cache_cols) || begin
                dropped=falses(size(cached.tris,2))
                dropped[cache_cols].=true
                _remove_element_columns(cached,class,2,dropped)
            end
            if record!==nothing
                keep=findall(t->t!=2,record.element_types)
                if length(keep)!=length(record.element_types)
                    record.element_types=record.element_types[keep]
                    record.element_tags=record.element_tags[keep]
                    record.element_nodes=record.element_nodes[keep]
                end
            end
        end
        return nothing
    end
end

# Lowest fresh node tag above the dense cache range and every record tag.
function _model_fresh_node_tag(m::GeoModel,cached)
    result=cached===nothing ? 0 : size(cached.coords,2)
    for (_,_,record) in _discrete_mesh_records(m)
        isempty(record.node_tags) ||
            (result=max(result,Int(maximum(record.node_tags))))
    end
    return Int32(result+1)
end

# Corner-angle quality of a quadrangle in [0,1]: 1 for a rectangle, 0 when a
# corner collapses. Each corner scores min(θ, π-θ)·2/π.
function _quad_quality(coords::NTuple{4,NTuple{3,Float64}})
    worst=1.0
    for k in 1:4
        a=coords[mod1(k-1,4)];b=coords[k];c=coords[mod1(k+1,4)]
        u=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
        v=(c[1]-b[1],c[2]-b[2],c[3]-b[3])
        nu=sqrt(u[1]^2+u[2]^2+u[3]^2);nv=sqrt(v[1]^2+v[2]^2+v[3]^2)
        (nu>0 && nv>0) || return 0.0
        cosine=clamp((u[1]*v[1]+u[2]*v[2]+u[3]*v[3])/(nu*nv),-1.0,1.0)
        theta=acos(cosine)
        worst=min(worst,min(theta,pi-theta)*2/pi)
    end
    return worst
end

function _split_quadrangles(quality=1.0,tag=-1)
    caller="API.mesh.split_quadrangles"
    (quality isa Real && !(quality isa Bool) && isfinite(quality)) ||
        throw(ArgumentError("$caller: quality must be a finite real number"))
    limit=Float64(quality)
    tag_value=_mesh_query_integer(tag,caller,"tag")
    return lock(STATE_LOCK) do
        m=_model_locked()
        cached=LAST_MESH[]
        keys=Tuple{Int,Int}[]
        if tag_value<0
            for (edim,etag,record) in _discrete_mesh_records(m)
                edim==2 || continue
                any(t->msh_family(t)===:qua,record.element_types) ||
                    continue
                push!(keys,(edim,etag))
            end
        else
            tag_value<=typemax(Int32) || throw(ArgumentError(
                "$caller: tag is out of range"))
            _model_entity_known(m,2,Int(tag_value)) || throw(ArgumentError(
                "$caller: unknown entity (2, $tag_value)"))
            push!(keys,(2,Int(tag_value)))
        end
        for key in sort!(keys)
            record=_entity_record(m,key[1],key[2])
            record===nothing && continue
            index=_record_node_index(m)
            fresh=_model_fresh_element_tag(m)
            index2=1
            while index2<=length(record.element_types)
                msh_type=Int(record.element_types[index2])
                if msh_family(msh_type)!==:qua
                    index2+=1
                    continue
                end
                nodes=record.element_nodes[index2]
                corners=nodes[1:4]
                points=ntuple(k->begin
                        value=_mesh_node_coords(m,cached,corners[k])
                        value===nothing && throw(ArgumentError(
                            "$caller: node $(corners[k]) has no coordinates"))
                        value
                    end,4)
                if _quad_quality(points)>=limit
                    index2+=1
                    continue
                end
                # Split along the diagonal yielding the better minimum corner
                # angle: (1,3) vs (2,4).
                function tri_min(a,b,c)
                    pa,pb,pc=points[a],points[b],points[c]
                    u=(pa[1]-pb[1],pa[2]-pb[2],pa[3]-pb[3])
                    v=(pc[1]-pb[1],pc[2]-pb[2],pc[3]-pb[3])
                    nu=sqrt(u[1]^2+u[2]^2+u[3]^2)
                    nv=sqrt(v[1]^2+v[2]^2+v[3]^2)
                    (nu>0 && nv>0) || return 0.0
                    return acos(clamp(
                        (u[1]*v[1]+u[2]*v[2]+u[3]*v[3])/(nu*nv),-1.0,1.0))
                end
                a13=min(tri_min(1,2,3),tri_min(3,4,1),
                        tri_min(2,3,1),tri_min(4,1,3))
                a24=min(tri_min(2,3,4),tri_min(4,1,2),
                        tri_min(3,4,2),tri_min(1,2,4))
                keep=nodes[5:end]
                if a13>=a24
                    tri1=vcat(nodes[1:3],keep)
                    tri2=vcat(nodes[[3,4,1]],keep)
                else
                    tri1=vcat(nodes[2:4],keep)
                    tri2=vcat(nodes[[4,1,2]],keep)
                end
                quad_tag=record.element_tags[index2]
                record.element_types[index2]=Int32(2)
                record.element_nodes[index2]=tri1
                _record_append_element!(record,2,Int32(fresh),tri2)
                fresh+=1
                index2+=1
            end
        end
        return nothing
    end
end



# ---- session views --------------------------------------------------------
#
# A session view holds named post-processing data: per-element list data in
# the Gmsh layout (coordinates then components per node), homogeneous
# model data, and numeric/string options. `compute_cross_field` fills three
# views (size H, angle Theta, cross directions); the `view` submodule
# manages them.

mutable struct _SessionView
    name::String
    numbers::Dict{String,Float64}
    strings::Dict{String,String}
    list_types::Vector{String}
    list_counts::Dict{String,Int}
    list_data::Dict{String,Vector{Float64}}
    model_data::Dict{Tuple{Int,String},Tuple{Vector{UInt64},Vector{Float64},Int}}
end

_SessionView(name::AbstractString)=_SessionView(
    String(name),Dict{String,Float64}(),Dict{String,String}(),
    String[],Dict{String,Int}(),Dict{String,Vector{Float64}}(),
    Dict{Tuple{Int,String},Tuple{Vector{UInt64},Vector{Float64},Int}}())

const SESSION_VIEWS = Dict{Int,_SessionView}()
const VIEW_TAG_MAX = Ref(0)

function _fresh_view_tag()
    while true
        VIEW_TAG_MAX[]+=1
        haskey(SESSION_VIEWS,VIEW_TAG_MAX[]) || return VIEW_TAG_MAX[]
    end
end

function _view_add(name,tag=-1)
    caller="API.view.add"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    value=_mesh_query_integer(tag,caller,"tag")
    return lock(STATE_LOCK) do
        allocated=value<0 ? _fresh_view_tag() : Int(value)
        haskey(SESSION_VIEWS,allocated) && throw(ArgumentError(
            "$caller: view tag $allocated already exists"))
        allocated<=typemax(Int32) || throw(ArgumentError(
            "$caller: tag is out of range"))
        SESSION_VIEWS[allocated]=_SessionView(name)
        VIEW_TAG_MAX[]=max(VIEW_TAG_MAX[],allocated)
        return allocated
    end
end

function _session_view(tag,caller::AbstractString)
    value=_mesh_query_integer(tag,caller,"tag")
    view=get(SESSION_VIEWS,Int(value),nothing)
    view===nothing && throw(ArgumentError(
        "$caller: unknown view tag $value"))
    return view
end

function _view_remove(tag)
    lock(STATE_LOCK) do
        view=_session_view(tag,"API.view.remove")
        delete!(SESSION_VIEWS,Int(tag))
        return nothing
    end
end

function _view_get_index(tag)
    value=_mesh_query_integer(tag,"API.view.get_index","tag")
    return lock(STATE_LOCK) do
        view_tags=sort!(collect(keys(SESSION_VIEWS)))
        index=findfirst(==(Int(value)),view_tags)
        index===nothing && return Int32(-1)
        return Int32(index-1)
    end
end

_view_get_tags()=lock(STATE_LOCK) do
    Int32[tag for tag in sort!(collect(keys(SESSION_VIEWS)))]
end

function _view_add_list_data(tag,data_type,num_ele,data)
    caller="API.view.add_list_data"
    data_type isa AbstractString || throw(ArgumentError(
        "$caller: dataType must be a string"))
    count=_mesh_query_integer(num_ele,caller,"numEle")
    count<0 && throw(ArgumentError("$caller: numEle must be non-negative"))
    (data isa AbstractVector || data isa Tuple) || throw(ArgumentError(
        "$caller: data must be a vector"))
    return lock(STATE_LOCK) do
        view=_session_view(tag,caller)
        values=Vector{Float64}(undef,length(data))
        for (i,value) in enumerate(data)
            values[i]=_mesh_query_coordinate(value,caller,"data entry")
        end
        type_key=String(data_type)
        if !(type_key in view.list_types)
            push!(view.list_types,type_key)
        end
        view.list_counts[type_key]=Int(count)
        view.list_data[type_key]=values
        return nothing
    end
end

function _view_get_list_data(tag,return_adaptive=false)
    caller="API.view.get_list_data"
    _mesh_query_bool(return_adaptive,caller,"returnAdaptive")
    return lock(STATE_LOCK) do
        view=_session_view(tag,caller)
        types=copy(view.list_types)
        data=[copy(view.list_data[type_key])
              for type_key in view.list_types]
        return types,data
    end
end

function _view_add_model_data(tag,step,model_name,data_type,tags,data,
                              time=0.0,num_components=-1,partition=0)
    caller="API.view.add_model_data"
    data_type isa AbstractString || throw(ArgumentError(
        "$caller: dataType must be a string"))
    model_name isa AbstractString || throw(ArgumentError(
        "$caller: modelName must be a string"))
    step_value=_mesh_query_integer(step,caller,"step")
    _mesh_query_coordinate(time,caller,"time")
    (tags isa AbstractVector || tags isa Tuple) || throw(ArgumentError(
        "$caller: tags must be a vector"))
    (data isa AbstractVector || data isa Tuple) || throw(ArgumentError(
        "$caller: data must be a vector"))
    length(tags)>0 && length(data)%length(tags)==0 || throw(ArgumentError(
        "$caller: data length must be a multiple of tags"))
    components=_mesh_query_integer(num_components,caller,"numComponents")
    return lock(STATE_LOCK) do
        view=_session_view(tag,caller)
        entity_tags=Vector{UInt64}(undef,length(tags))
        for (i,value) in enumerate(tags)
            entry=_mesh_query_integer(value,caller,"tags entry")
            entry>=1 || throw(ArgumentError(
                "$caller: tags entry must be positive"))
            entity_tags[i]=UInt64(entry)
        end
        values=Vector{Float64}(undef,length(data))
        for (i,value) in enumerate(data)
            values[i]=_mesh_query_coordinate(value,caller,"data entry")
        end
        view.model_data[(Int(step_value),String(data_type))]=
            (entity_tags,values,Int(components))
        return nothing
    end
end

function _view_get_model_data(tag,step)
    caller="API.view.get_model_data"
    step_value=_mesh_query_integer(step,caller,"step")
    return lock(STATE_LOCK) do
        view=_session_view(tag,caller)
        types=String[]
        tags=UInt64[]
        data=Float64[]
        for ((s,data_type),(entity_tags,values,_)) in
                sort!(collect(view.model_data);by=first)
            s==step_value || continue
            push!(types,data_type)
            append!(tags,entity_tags)
            append!(data,values)
        end
        return types,tags,data
    end
end

function _view_set_number(tag,name,value)
    caller="API.view.set_number"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    numeric=_mesh_query_coordinate(value,caller,"value")
    return lock(STATE_LOCK) do
        _session_view(tag,caller).numbers[String(name)]=numeric
        return nothing
    end
end

function _view_get_number(tag,name)
    caller="API.view.get_number"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    return lock(STATE_LOCK) do
        get(_session_view(tag,caller).numbers,String(name),0.0)
    end
end

function _view_set_string(tag,name,value)
    caller="API.view.set_string"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    value isa AbstractString || throw(ArgumentError(
        "$caller: value must be a string"))
    return lock(STATE_LOCK) do
        _session_view(tag,caller).strings[String(name)]=String(value)
        return nothing
    end
end

function _view_get_string(tag,name)
    caller="API.view.get_string"
    name isa AbstractString || throw(ArgumentError(
        "$caller: name must be a string"))
    return lock(STATE_LOCK) do
        get(_session_view(tag,caller).strings,String(name),"")
    end
end

# ---- compute_cross_field --------------------------------------------------
#
# A per-element cross field on the surface mesh: each 2D element's cross is
# the unit direction of its first edge transported across shared edges and
# smoothed by pi/2-symmetric Jacobi relaxation. Three views are created:
# element size H, the angle Theta between the cross and the first edge, and
# the cross direction vectors — matching Gmsh's `mesh.computeCrossField`.

# Tangent-plane unit vector of an element's first edge.
function _crossfield_frame(points)
    a,b,c=points[1],points[2],points[3]
    normal=_model_cross3((b[1]-a[1],b[2]-a[2],b[3]-a[3]),
                         (c[1]-a[1],c[2]-a[2],c[3]-a[3]))
    norm=sqrt(normal[1]^2+normal[2]^2+normal[3]^2)
    norm>0 || return nothing,nothing
    normal=normal./norm
    edge=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
    magnitude=sqrt(edge[1]^2+edge[2]^2+edge[3]^2)
    magnitude>0 || return nothing,nothing
    return [edge[1]/magnitude,edge[2]/magnitude,edge[3]/magnitude],normal
end

# Rotate `direction` into `normal`'s tangent plane, snapped to the pi/2
# representative closest to `reference`.
function _crossfield_transport(direction,normal,reference)
    projected=(direction[1]-normal[1]*(direction[1]*normal[1]+
                                       direction[2]*normal[2]+
                                       direction[3]*normal[3]),
               direction[2]-normal[2]*(direction[1]*normal[1]+
                                       direction[2]*normal[2]+
                                       direction[3]*normal[3]),
               direction[3]-normal[3]*(direction[1]*normal[1]+
                                       direction[2]*normal[2]+
                                       direction[3]*normal[3]))
    magnitude=sqrt(projected[1]^2+projected[2]^2+projected[3]^2)
    magnitude>0 || return copy(reference)
    # explicit arithmetic: closures over the reassigned `projected`/`current`
    # boxed both and allocated on every transported direction
    unit=(projected[1]/magnitude,projected[2]/magnitude,projected[3]/magnitude)
    best=[unit[1],unit[2],unit[3]]
    best_dot=abs(unit[1]*reference[1]+unit[2]*reference[2]+unit[3]*reference[3])
    current=[unit[1],unit[2],unit[3]]
    for _ in 1:3
        rotated=_model_cross3(normal,current)
        rmag=sqrt(rotated[1]^2+rotated[2]^2+rotated[3]^2)
        rmag>0 || break
        current=[rotated[1]/rmag,rotated[2]/rmag,rotated[3]/rmag]
        dot=abs(current[1]*reference[1]+current[2]*reference[2]+
                current[3]*reference[3])
        dot>best_dot && (best=copy(current);best_dot=dot)
    end
    return best
end

function _compute_cross_field()
    caller="API.mesh.compute_cross_field"
    return lock(STATE_LOCK) do
        m=_model_locked()
        cached=LAST_MESH[]
        elements=NTuple{3,NTuple{3,Float64}}[]
        # Surface elements from cache and records, in entity order. A volume
        # cache stores only tets, so fall back to its boundary faces — faces
        # incident to exactly one tet.
        if cached!==nothing && size(cached.tris,2)>0
            for column in axes(cached.tris,2)
                nodes=cached.tris[:,column]
                push!(elements,ntuple(k->(
                    cached.coords[1,nodes[k]],cached.coords[2,nodes[k]],
                    cached.coords[3,nodes[k]]),3))
            end
        elseif cached!==nothing && size(cached.tets,2)>0
            incidence=Dict{NTuple{3,Int32},Int}()
            for column in axes(cached.tets,2)
                for pattern in _simplex_face_patterns(4,3)
                    key=_triangle_key(
                        cached.tets[pattern[1],column],
                        cached.tets[pattern[2],column],
                        cached.tets[pattern[3],column])
                    incidence[key]=get(incidence,key,0)+1
                end
            end
            for (key,count) in incidence
                count==1 || continue
                push!(elements,ntuple(k->(
                    cached.coords[1,key[k]],cached.coords[2,key[k]],
                    cached.coords[3,key[k]]),3))
            end
            sort!(elements)
        end
        index=_record_node_index(m)
        for (edim,_,record) in sort!(collect(_discrete_mesh_records(m));
                                     by=entry->(entry[1],entry[2]))
            edim==2 || continue
            for e in eachindex(record.element_types)
                msh_family(record.element_types[e])===:tri || continue
                nodes=record.element_nodes[e]
                length(nodes)>=3 || continue
                any(k->!haskey(index,nodes[k]),1:3) && continue
                frame=ntuple(3) do k
                    record2,column=index[nodes[k]]
                    (record2.node_coords[1,column],
                     record2.node_coords[2,column],
                     record2.node_coords[3,column])
                end
                push!(elements,frame)
            end
        end
        isempty(elements) && throw(ArgumentError(
            "$caller: no surface elements to compute a cross field on"))
        normals=Vector{NTuple{3,Float64}}(undef,length(elements))
        directions=Vector{Vector{Float64}}(undef,length(elements))
        edge_lookup=Dict{NTuple{2,Int},Vector{Int}}()
        coordinates=NTuple{3,Float64}[]
        node_position=Dict{NTuple{3,Float64},Int}()
        element_nodes=Vector{NTuple{3,Int}}(undef,length(elements))
        for (i,frame) in enumerate(elements)
            direction,normal=_crossfield_frame(frame)
            (direction===nothing) && throw(ArgumentError(
                "$caller: degenerate surface element $i"))
            normals[i]=normal
            directions[i]=direction
            position=ntuple(k->get!(node_position,frame[k]) do
                push!(coordinates,frame[k])
                length(coordinates)
            end,3)
            element_nodes[i]=position
            for k in 1:3
                a=position[k];b=position[mod1(k+1,3)]
                push!(get!(Vector{Int},edge_lookup,
                           a<b ? (a,b) : (b,a)),i)
            end
        end
        # Element adjacency through shared edges.
        adjacency=[Int[] for _ in elements]
        for (_,owners) in edge_lookup
            for x in owners,y in owners
                x!=y && push!(adjacency[x],y)
            end
        end
        # pi/2-symmetric Jacobi smoothing.
        for _ in 1:8
            update=Vector{Vector{Float64}}(undef,length(elements))
            for i in eachindex(elements)
                sumd=copy(directions[i])
                for j in adjacency[i]
                    transported=_crossfield_transport(
                        directions[j],normals[i],directions[i])
                    for k in 1:3
                        sumd[k]+=transported[k]
                    end
                end
                magnitude=sqrt(sumd[1]^2+sumd[2]^2+sumd[3]^2)
                update[i]=magnitude>0 ? sumd./magnitude : directions[i]
            end
            directions=update
        end
        # Emit three views: H (element size), Theta (angle vs first edge),
        # cross (vector per node, replicated element-constant).
        sizes=Float64[];angles=Float64[];vectors=Float64[]
        for (i,frame) in enumerate(elements)
            a,b,c=frame
            lengths=[sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2),
                     sqrt((b[1]-c[1])^2+(b[2]-c[2])^2+(b[3]-c[3])^2),
                     sqrt((c[1]-a[1])^2+(c[2]-a[2])^2+(c[3]-a[3])^2)]
            push!(sizes,sum(lengths)/3)
            direction,normal=_crossfield_frame(frame)
            push!(angles,acos(clamp(
                sum(k->directions[i][k]*direction[k],1:3),-1.0,1.0)))
            for k in 1:3
                append!(vectors,frame[k])
                append!(vectors,directions[i])
            end
        end
        h_tag=_fresh_view_tag()
        SESSION_VIEWS[h_tag]=_SessionView("H")
        _session_view(h_tag,caller).list_types=["ST"]
        _session_view(h_tag,caller).list_counts["ST"]=length(elements)
        h_data=Float64[]
        for (i,frame) in enumerate(elements)
            for node in frame
                append!(h_data,node);push!(h_data,sizes[i])
            end
        end
        _session_view(h_tag,caller).list_data["ST"]=h_data
        theta_tag=_fresh_view_tag()
        SESSION_VIEWS[theta_tag]=_SessionView("Theta")
        theta_data=Float64[]
        for (i,frame) in enumerate(elements)
            for node in frame
                append!(theta_data,node);push!(theta_data,angles[i])
            end
        end
        _session_view(theta_tag,caller).list_types=["ST"]
        _session_view(theta_tag,caller).list_counts["ST"]=length(elements)
        _session_view(theta_tag,caller).list_data["ST"]=theta_data
        cross_tag=_fresh_view_tag()
        SESSION_VIEWS[cross_tag]=_SessionView("cross")
        _session_view(cross_tag,caller).list_types=["VT"]
        _session_view(cross_tag,caller).list_counts["VT"]=length(elements)
        _session_view(cross_tag,caller).list_data["VT"]=vectors
        return Int32[h_tag,theta_tag,cross_tag]
    end
end
# ---- session size fields --------------------------------------------------
#
# `model.mesh.field` state lives in the current model slot. Options are kept
# as strings in first-set order — identical to `GeoFieldSpec` — so generation
# builds the background field through `build_geo_size_field` unchanged.

function _session_field(tag,caller::AbstractString)
    value=_mesh_query_integer(tag,caller,"tag")
    field=get(FIELD_SPECS[],Int(value),nothing)
    field===nothing && throw(ArgumentError(
        "$caller: unknown field tag $value"))
    return field
end

function _field_add(field_type,tag=-1)
    caller="API.field.add"
    field_type isa AbstractString || throw(ArgumentError(
        "$caller: fieldType must be a string"))
    isempty(field_type) && throw(ArgumentError(
        "$caller: fieldType must not be empty"))
    value=_mesh_query_integer(tag,caller,"tag")
    return lock(STATE_LOCK) do
        _model_locked()
        allocated=_field_alloc_tag(value)
        allocated<1 && throw(ArgumentError(
            "$caller: tag must be positive"))
        allocated>typemax(Int32) && throw(ArgumentError(
            "$caller: tag is out of range"))
        haskey(FIELD_SPECS[],allocated) && throw(ArgumentError(
            "$caller: field tag $allocated already exists"))
        FIELD_SPECS[][allocated]=
            _ApiField(String(field_type),Dict{String,String}(),String[])
        FIELD_TAG_MAX[]=max(FIELD_TAG_MAX[],allocated)
        return Int32(allocated)
    end
end

function _field_alloc_tag(value::Int)
    value>=1 && return value
    candidate=FIELD_TAG_MAX[]+1
    while haskey(FIELD_SPECS[],candidate)
        candidate+=1
    end
    return candidate
end

function _field_remove(tag)
    return lock(STATE_LOCK) do
        _model_locked()
        _session_field(tag,"API.field.remove")
        value=Int(_mesh_query_integer(tag,"API.field.remove","tag"))
        delete!(FIELD_SPECS[],value)
        BACKGROUND_FIELD[]==value && (BACKGROUND_FIELD[]=0)
        filter!(!=(value),BOUNDARY_LAYER_FIELDS[])
        return nothing
    end
end

_field_list()=lock(STATE_LOCK) do
    _model_locked()
    Int32[tag for tag in sort!(collect(keys(FIELD_SPECS[])))]
end

function _field_get_type(tag)
    return lock(STATE_LOCK) do
        _model_locked()
        _session_field(tag,"API.field.get_type").kind
    end
end

function _field_set_option(tag,name,value::AbstractString,caller::AbstractString)
    return lock(STATE_LOCK) do
        _model_locked()
        field=_session_field(tag,caller)
        key=String(name)
        haskey(field.options,key) || push!(field.order,key)
        field.options[key]=String(value)
        return nothing
    end
end

function _field_set_number(tag,name,value)
    caller="API.field.set_number"
    name isa AbstractString || throw(ArgumentError(
        "$caller: option must be a string"))
    numeric=_mesh_query_coordinate(value,caller,"value")
    return _field_set_option(tag,name,string(numeric),caller)
end

function _field_get_number(tag,name)
    caller="API.field.get_number"
    name isa AbstractString || throw(ArgumentError(
        "$caller: option must be a string"))
    return lock(STATE_LOCK) do
        _model_locked()
        raw=get(_session_field(tag,caller).options,String(name),nothing)
        raw===nothing && return 0.0
        value=tryparse(Float64,raw)
        return value===nothing ? 0.0 : value
    end
end

function _field_set_string(tag,name,value)
    caller="API.field.set_string"
    name isa AbstractString || throw(ArgumentError(
        "$caller: option must be a string"))
    value isa AbstractString || throw(ArgumentError(
        "$caller: value must be a string"))
    return _field_set_option(tag,name,value,caller)
end

function _field_get_string(tag,name)
    caller="API.field.get_string"
    name isa AbstractString || throw(ArgumentError(
        "$caller: option must be a string"))
    return lock(STATE_LOCK) do
        _model_locked()
        get(_session_field(tag,caller).options,String(name),"")
    end
end

function _field_set_numbers(tag,name,values)
    caller="API.field.set_numbers"
    name isa AbstractString || throw(ArgumentError(
        "$caller: option must be a string"))
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: values must be a vector"))
    items=String[]
    for value in values
        numeric=_mesh_query_coordinate(value,caller,"values entry")
        push!(items,isinteger(numeric) ? string(Int(numeric)) :
              string(numeric))
    end
    return _field_set_option(tag,name,"{"*join(items,",")*"}",caller)
end

function _field_get_numbers(tag,name)
    caller="API.field.get_numbers"
    name isa AbstractString || throw(ArgumentError(
        "$caller: option must be a string"))
    return lock(STATE_LOCK) do
        _model_locked()
        raw=get(_session_field(tag,caller).options,String(name),nothing)
        raw===nothing && return Float64[]
        items=_geo_split_list(raw,caller)
        result=Float64[]
        for item in items
            value=tryparse(Float64,String(strip(item)))
            value===nothing && throw(ArgumentError(
                "$caller: option $name entry $item is not numeric"))
            push!(result,value)
        end
        return result
    end
end

function _field_set_as_background_mesh(tag)
    return lock(STATE_LOCK) do
        _model_locked()
        value=Int(_mesh_query_integer(
            tag,"API.field.set_as_background_mesh","tag"))
        _session_field(value,"API.field.set_as_background_mesh")
        BACKGROUND_FIELD[]=value
        return nothing
    end
end

function _field_set_as_boundary_layer(tag)
    return lock(STATE_LOCK) do
        _model_locked()
        value=Int(_mesh_query_integer(
            tag,"API.field.set_as_boundary_layer","tag"))
        _session_field(value,"API.field.set_as_boundary_layer")
        value in BOUNDARY_LAYER_FIELDS[] ||
            push!(BOUNDARY_LAYER_FIELDS[],value)
        return nothing
    end
end

# Resolve `(dimension, tag)` entity references inside field options to a Mesh,
# as `build_geo_size_field` expects: points become a single-node mesh, curves a
# single-segment mesh, and planar surfaces their boundary PSLG triangulation.
function _field_entity_mesh(m::GeoModel,dim::Int,name::AbstractString)
    tag=tryparse(Int,String(name))
    if tag===nothing
        numeric=tryparse(Float64,String(name))
        (numeric===nothing || !isinteger(numeric) ||
         numeric<typemin(Int) || numeric>typemax(Int)) && return nothing
        tag=Int(numeric)
    end
    if dim==0
        haskey(m.points,tag) || return nothing
        return Mesh(reshape(collect(Float64.(m.points[tag])),3,1))
    elseif dim==1
        haskey(m.curves,tag) || return nothing
        a,b=m.curves[tag]
        coords=hcat(collect(Float64.(m.points[a])),
                    collect(Float64.(m.points[b])))
        return Mesh(coords;segs=reshape(Int32[1;2],2,1))
    elseif dim==2
        haskey(m.surfaces,tag) || return nothing
        return _model_planar_surface_mesh(
            m,tag,"API.field";include_embeddings=false)
    end
    return nothing
end

# Build the model's background size field, or `nothing` when none is selected.
# `context_fields` handles the kinds needing model/view state; `PostView`
# resolves against the session view store.
function _session_size_field_locked(m::GeoModel)
    BACKGROUND_FIELD[]>0 || return nothing
    specs=Dict{Int,GeoFieldSpec}()
    for (tag,field) in FIELD_SPECS[]
        specs[tag]=GeoFieldSpec(tag,field.kind,Dict(field.options),
                                copy(field.order))
    end
    params=GeoParams(OPTIONS["Mesh.MeshSizeMin"],OPTIONS["Mesh.MeshSizeMax"],
                     OPTIONS["Mesh.MeshSizeFactor"],0,
                     Dict{Tuple{Int,Int},String}(),specs,
                     BACKGROUND_FIELD[],copy(BOUNDARY_LAYER_FIELDS[]))
    bounds=try
        model_bounding_box(m,-1,-1)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    function context_fields(spec,config,entities,_params)
        kind=lowercase(spec.kind)
        if kind=="postview"
            view_tag=haskey(spec.options,"ViewTag") ?
                tryparse(Int,strip(spec.options["ViewTag"])) :
                tryparse(Int,strip(get(spec.options,"IField","1")))
            view_tag===nothing && throw(ArgumentError(
                "build_geo_size_field: Field[$(spec.tag)].ViewTag must be an integer"))
            view=get(SESSION_VIEWS,view_tag,nothing)
            view===nothing && throw(ArgumentError(
                "build_geo_size_field: PostView Field[$(spec.tag)] references " *
                "unknown view $view_tag"))
            return _post_view_field(view,spec.tag)
        end
        throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] kind $(spec.kind) " *
            "requires model context that is not available to this session"))
    end
    return build_geo_size_field(
        params,(d,name)->_field_entity_mesh(m,d,name);
        model_bbox=bounds===nothing ? nothing :
                  (bounds[1],bounds[4],bounds[2],bounds[5],bounds[3],bounds[6]),
        context_fields=context_fields)
end

# Node counts for Gmsh list-data element shapes. "S" is ambiguous between
# points (1 node) and tetrahedra (4) — Gmsh disambiguates by data length.
const _VIEW_SHAPE_NODES=Dict{Char,Vector{Int}}(
    'S'=>[1,4],'L'=>[2],'T'=>[3],'Q'=>[4],'H'=>[8],'I'=>[6],'Y'=>[5])

# Wrap a session view's scalar list data as a `PostViewField` size source.
# Each element row holds `nodes×3` coordinates then `nodes` scalar values.
function _post_view_field(view::_SessionView,tag::Int)
    caller="build_geo_size_field"
    coords=Float64[];values=Float64[]
    element_cells=Pair{Symbol,Matrix{Int32}}[]
    for type_key in view.list_types
        length(type_key)==2 || continue
        type_key[1]=='S' || throw(ArgumentError(
            "$caller: Field[$tag] requires scalar view data (got $type_key)"))
        data=view.list_data[type_key]
        count=view.list_counts[type_key]
        candidates=get(_VIEW_SHAPE_NODES,type_key[2],nothing)
        candidates===nothing && throw(ArgumentError(
            "$caller: unsupported view data type $type_key for Field[$tag]"))
        nnodes_elem=findfirst(n->count*(n*3+n)==length(data),candidates)
        nnodes_elem===nothing && throw(ArgumentError(
            "$caller: $type_key data length $(length(data)) is inconsistent " *
            "with numEle=$count"))
        n=candidates[nnodes_elem]
        connectivity=Matrix{Int32}(undef,n,count)
        cursor=0
        for e in 1:count
            first_node=length(values)+1
            for k in 1:n
                append!(coords,data[cursor+1:cursor+3])
                cursor+=3
                connectivity[k,e]=Int32(first_node+k-1)
            end
            append!(values,data[cursor+1:cursor+n])
            cursor+=n
        end
        key=type_key[2]=='S' ? (n==1 ? :points : :tetrahedra) :
            type_key[2]=='L' ? :lines :
            type_key[2]=='T' ? :triangles : type_key[2]=='Q' ? :quadrangles :
            type_key[2]=='H' ? :hexahedra : type_key[2]=='I' ? :prisms :
            :pyramids
        push!(element_cells,key=>connectivity)
    end
    isempty(element_cells) && throw(ArgumentError(
        "$caller: view backing Field[$tag] has no scalar list data"))
    return PostViewField(reshape(coords,3,:),values;
                         NamedTuple(element_cells)...)
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
             _remove_embedded,_embed,
             _get_ghost_elements,_unpartition,_partition,_import_stl,_create_topology,_classify_surfaces,_create_geometry,
              _get_periodic_keys,_recombine,
              _split_quadrangles,_get_last_entity_error,
             _get_last_node_error,_rebuild_node_cache,_rebuild_element_cache,
             _reclassify_nodes,_relocate_nodes,
             _get_periodic,_remove_constraints,_compute_renumbering,
             _optimize_mesh,_set_mesh_visibility,_get_mesh_visibility,
             _set_mesh_visibility_per_window,
             _set_transfinite_curve,_set_transfinite_surface,
             _set_transfinite_volume,_set_transfinite_automatic,
             _set_recombine_attribute,_set_smoothing,_set_reverse_attribute,
             _set_algorithm,_set_size_at_parametric_points,
             _set_size_from_boundary,_set_size_callback,_remove_size_callback,
             _set_compound,_set_outward_orientation,_set_order,
             _add_nodes,_add_elements,_add_elements_by_type,
             _add_homology_request,_clear_homology_requests,
              _compute_homology,_compute_cross_field
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
    set_visibility_per_window(tag, value, window_index=0)

Store the display visibility flag of dense element `tag` for window
`window_index`, matching Gmsh's `mesh.setVisibilityPerWindow`. Unknown element
tags are dropped silently like `set_visibility`; there is no corresponding
getter.
"""
set_visibility_per_window(tag,value,window_index=0)=
    _set_mesh_visibility_per_window(tag,value,window_index)

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
"Partition the current mesh into `num_part` parts (or an explicit element partition list)."
partition(num_part,element_tags=(),partitions=())=_partition(num_part,element_tags,partitions)
"Import the model file's STL triangulation (when available) as discrete mesh data."
import_stl()=_import_stl()
"Create a boundary representation from discrete mesh data (see `mesh.createTopology`)."
create_topology(make_simply_connected=true,export_discrete=true)=
    _create_topology(make_simply_connected,export_discrete)
"Classify the discrete surface mesh into angle-limited patches (see `mesh.classifySurfaces`)."
classify_surfaces(angle=40*pi/180,boundary=true,for_reparametrization=false,
                  curve_angle=pi)=
    _classify_surfaces(angle,boundary,for_reparametrization,curve_angle)
"Recombine flagged surface triangles into quadrangles (see `mesh.recombine`)."
recombine()=_recombine()
"Split low-quality quadrangles into triangles (see `mesh.splitQuadrangles`)."
split_quadrangles(quality=1.0,tag=-1)=_split_quadrangles(quality,tag)
"Compute parametrizations for discrete entities (see `mesh.createGeometry`)."
create_geometry(dim_tags=Tuple{Int,Int}[])=_create_geometry(dim_tags)
"Return periodic function-space keys for an entity (see `mesh.getPeriodicKeys`)."
get_periodic_keys(element_type,function_space_type,tag,return_coord=true)=
    _get_periodic_keys(element_type,function_space_type,tag,return_coord)
"Perform queued homology requests (see `mesh.computeHomology`)."
compute_homology()=_compute_homology()

"""
    compute_cross_field() -> Vector{Int32}

Compute a smoothed per-element cross field on the surface mesh and store it
in three new session views — element size `H`, angle `Theta` between each
cross and its element's first edge, and the `cross` directions — returning
their tags like `mesh.computeCrossField` in Gmsh 4.15.2.
"""
compute_cross_field()=_compute_cross_field()
"Embed entities of dimension `dim` and tags `tags` into entity (`in_dim`, `in_tag`)."
embed(dim,tags,in_dim,in_tag)=_embed(dim,tags,in_dim,in_tag)

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

"""
    set_transfinite_curve(tag, num_nodes, mesh_type="Progression", coef=1.0)

Record a transfinite meshing constraint on `Curve[tag]`, matching Gmsh's
`setTransfiniteCurve`. `mesh_type` is `"Progression"`, `"Bump"`, or `"Beta"`;
the constraint is consumed by `generate`.
"""
set_transfinite_curve(tag,num_nodes,mesh_type="Progression",coef=1.0)=
    _set_transfinite_curve(tag,num_nodes,mesh_type,coef)

"""
    set_transfinite_surface(tag, arrangement="Left", corner_tags=[])

Record a transfinite constraint on `Surface[tag]`, matching Gmsh's
`setTransfiniteSurface`. `arrangement` is `"Left"`, `"Right"`,
`"AlternateLeft"`, or `"AlternateRight"`; `corner_tags` optionally pins 3 or 4
corner Points.
"""
set_transfinite_surface(tag,arrangement="Left",corner_tags=[])=
    _set_transfinite_surface(tag,arrangement,corner_tags)

"""
    set_transfinite_volume(tag, corner_tags=[])

Record a transfinite constraint on `Volume[tag]`, matching Gmsh's
`setTransfiniteVolume`. `corner_tags` optionally pins the 6 or 8 corner Points.
"""
set_transfinite_volume(tag,corner_tags=[])=
    _set_transfinite_volume(tag,corner_tags)

"""
    set_transfinite_automatic(dim_tags=[], corner_angle=2.35, recombine=true)

Apply transfinite constraints automatically, matching Gmsh's
`setTransfiniteAutomatic`: every listed surface (all surfaces when `dim_tags`
is empty) whose boundary resolves to a 3- or 4-cornered patch records a
transfinite surface and matching curve counts; every listed 6-sided volume
records a transfinite volume and processes its faces. Flatter-than-
`corner_angle` patch corners cause the face to be skipped, as upstream.
"""
set_transfinite_automatic(dim_tags=[],corner_angle=2.35,recombine=true)=
    _set_transfinite_automatic(dim_tags,corner_angle,recombine)

"""
    set_recombine(dim, tag, angle=45.0)

Record a recombination constraint on entity `(dim, tag)`, matching Gmsh's
`setRecombine` — generated triangles are recombined into quadrangles up to
`angle` degrees.
"""
set_recombine(dim,tag,angle=45.0)=_set_recombine_attribute(dim,tag,angle)

"""
    set_smoothing(dim, tag, val)

Record `val` Laplace-smoother iterations applied to the mesh of entity
`(dim, tag)` during generation, matching Gmsh's `setSmoothing`.
"""
set_smoothing(dim,tag,val)=_set_smoothing(dim,tag,val)

"""
    set_reverse(dim, tag, val=true)

Record a reverse-orientation constraint on entity `(dim, tag)`, matching
Gmsh's `setReverse` — generated element orientations are flipped relative to
the natural orientation when `val` is true.
"""
set_reverse(dim,tag,val=true)=_set_reverse_attribute(dim,tag,val)

"""
    set_algorithm(dim, tag, val)

Record the meshing algorithm for entity `(dim, tag)`, matching Gmsh's
`setAlgorithm`. `dim` must be 2 or 3; unsupported algorithm numbers are stored
and fail explicitly at `generate` time.
"""
set_algorithm(dim,tag,val)=_set_algorithm(dim,tag,val)

"""
    set_size_at_parametric_points(dim, tag, parametric_coord, sizes)

Record mesh-size constraints at parametric points on `Curve[tag]`, matching
Gmsh's `setSizeAtParametricPoints`. Only `dim == 1` is supported, as upstream.
"""
set_size_at_parametric_points(dim,tag,parametric_coord,sizes)=
    _set_size_at_parametric_points(dim,tag,parametric_coord,sizes)

"""
    set_size_from_boundary(dim, tag, val)

Record whether the interior mesh size of `Surface[tag]` extends its boundary
sizes, matching Gmsh's `setSizeFromBoundary`. Only `dim == 2` is supported, as
upstream.
"""
set_size_from_boundary(dim,tag,val)=_set_size_from_boundary(dim,tag,val)

"""
    set_size_callback(callback)

Install a model-level mesh-size callback `(dim, tag, x, y, z, lc) -> size`,
matching Gmsh's `setSizeCallback`. The generator calls it wherever it would
otherwise prescribe `lc`; returning `lc` is a no-op.
"""
set_size_callback(callback)=_set_size_callback(callback)

"""Remove the model-level mesh-size callback, matching `removeSizeCallback`."""
remove_size_callback()=_remove_size_callback()

"""
    set_compound(dim, tags)

Record a compound meshing constraint — entities `(dim, tags)` are meshed as a
single entity, matching Gmsh's `setCompound`. `dim` must be 1 or 2.
"""
set_compound(dim,tags)=_set_compound(dim,tags)

"""
    set_outward_orientation(tag)

Record that all boundary surfaces of `Volume[tag]` orient outward, matching
Gmsh's `setOutwardOrientation`.
"""
set_outward_orientation(tag)=_set_outward_orientation(tag)

"""
    set_order(order)

Record the global element order, matching Gmsh's `setOrder`. Orders 1 and 2 are
supported by generation; higher values are stored and fail at `generate` time.
"""
set_order(order)=_set_order(order)

"""
    add_nodes(dim, tag, node_tags, coord, parametric_coord=[])

Append nodes to entity `(dim, tag)`, matching Gmsh's `addNodes`. `node_tags`
are unique, strictly positive identifiers; `coord` is flattened x,y,z triples;
`parametric_coord` optionally carries `dim` values per node. Works on discrete
and native entities.
"""
add_nodes(dim,tag,node_tags,coord,parametric_coord=[])=
    _add_nodes(dim,tag,node_tags,coord,parametric_coord)

"""
    add_elements(dim, tag, element_types, element_tags, node_tags)

Append elements to entity `(dim, tag)`, matching Gmsh's `addElements`:
`element_types` lists MSH type numbers; `element_tags[i]` and `node_tags[i]`
hold the tag vector and flattened connectivity of block `i`. Referenced nodes
must already be classified on the entity.
"""
add_elements(dim,tag,element_types,element_tags,node_tags)=
    _add_elements(dim,tag,element_types,element_tags,node_tags)

"""
    add_elements_by_type(tag, element_type, element_tags, node_tags)

Append elements of a single MSH type to the entity owning `tag`, matching
Gmsh's `addElementsByType`. `node_tags` is the flattened connectivity.
"""
add_elements_by_type(tag,element_type,element_tags,node_tags)=
    _add_elements_by_type(tag,element_type,element_tags,node_tags)

"""
    add_homology_request(type="Homology", domain_tags=[], subdomain_tags=[],
                         dims=[])

Queue a (co)homology computation request, matching Gmsh's `addHomologyRequest`.
Requests are consumed by `compute_homology`.
"""
add_homology_request(type="Homology",domain_tags=[],subdomain_tags=[],dims=[])=
    _add_homology_request(type,domain_tags,subdomain_tags,dims)

"""Clear every queued homology request, matching `clearHomologyRequests`."""
clear_homology_requests()=_clear_homology_requests()

"""Mesh size fields (see the Gmsh 4.15.2 `model.mesh.field` module)."""
module field

using ...API: _field_add,_field_remove,_field_list,_field_get_type,
             _field_set_number,_field_get_number,_field_set_string,
             _field_get_string,_field_set_numbers,_field_get_numbers,
             _field_set_as_background_mesh,_field_set_as_boundary_layer

"""Add a field of `fieldType`, returning its tag."""
add(fieldType,tag=-1)=_field_add(fieldType,tag)
"Remove the field with tag `tag`."
remove(tag)=_field_remove(tag)
"Return all field tags in increasing order."
list()=_field_list()
"Return the type name of field `tag`."
get_type(tag)=_field_get_type(tag)
"Set numeric `option` of field `tag` to `value`."
set_number(tag,option,value)=_field_set_number(tag,option,value)
"Return numeric `option` of field `tag` (0 if unset)."
get_number(tag,option)=_field_get_number(tag,option)
"Set string `option` of field `tag` to `value`."
set_string(tag,option,value)=_field_set_string(tag,option,value)
"Return string `option` of field `tag` (empty if unset)."
get_string(tag,option)=_field_get_string(tag,option)
"Set numeric-list `option` of field `tag` to `values`."
set_numbers(tag,option,values)=_field_set_numbers(tag,option,values)
"Return numeric-list `option` of field `tag`."
get_numbers(tag,option)=_field_get_numbers(tag,option)
"Use field `tag` as the background mesh size field."
set_as_background_mesh(tag)=_field_set_as_background_mesh(tag)
"Use field `tag` as a boundary layer field."
set_as_boundary_layer(tag)=_field_set_as_boundary_layer(tag)

end # module field
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
        index=_current_slot_index_locked()
        CURRENT[]=stored_model
        MODEL_FILE_NAME[]=String(path)
        if index!=0
            MODEL_SLOTS[index].model=stored_model
            MODEL_SLOTS[index].file_name=String(path)
        end
        _replace_mesh_cache_locked!(stored_mesh)
        result
    end
end

"""Session post-processing views (see the Gmsh 4.15.2 `view` module)."""
module view

using ..API: _view_add,_view_remove,_view_get_index,_view_get_tags,
             _view_add_list_data,_view_get_list_data,_view_add_model_data,
             _view_get_model_data,_view_set_number,_view_get_number,
             _view_set_string,_view_get_string

"""Add a new empty view named `name`, returning its tag."""
add(name,tag=-1)=_view_add(name,tag)
"Remove the view with tag `tag`."
remove(tag)=_view_remove(tag)
"Return the index of view `tag` in the view list, or -1."
get_index(tag)=_view_get_index(tag)
"Return the tags of all views, in increasing tag order."
get_tags()=_view_get_tags()
"""Add per-element list `data` of `data_type` for `numEle` elements."""
add_list_data(tag,data_type,num_ele,data)=
    _view_add_list_data(tag,data_type,num_ele,data)
"""Return the list data types and data of view `tag`."""
get_list_data(tag,return_adaptive=false)=
    _view_get_list_data(tag,return_adaptive)
"""Add per-entity `data` of `data_type` at `step` to view `tag`."""
add_model_data(tag,step,model_name,data_type,tags,data,time=0.0,
               num_components=-1,partition=0)=
    _view_add_model_data(tag,step,model_name,data_type,tags,data,time,
                         num_components,partition)
"""Return the model data types, entity tags, and values of view `tag` at `step`."""
get_model_data(tag,step)=_view_get_model_data(tag,step)
"Set numeric option `name` of view `tag` to `value`."
set_number(tag,name,value)=_view_set_number(tag,name,value)
"Return numeric option `name` of view `tag` (0 if unset)."
get_number(tag,name)=_view_get_number(tag,name)
"Set string option `name` of view `tag` to `value`."
set_string(tag,name,value)=_view_set_string(tag,name,value)
"Return string option `name` of view `tag` (empty if unset)."
get_string(tag,name)=_view_get_string(tag,name)

end # module view

end # module
