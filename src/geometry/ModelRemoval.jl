function _model_removal_dim_tags(dim_tags,caller::AbstractString)
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dimension, tag) pairs"))
    result=Tuple{Int,Int}[]
    sizehint!(result,length(dim_tags))
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif entry isa Tuple && length(entry)==2
            entry
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dimension, tag) pair"))
        end
        dimension=_dimension(pair[1],caller)
        entity_tag=_tag(pair[2],caller,dimension)
        entity_tag>0 || throw(ArgumentError(
            "$caller: entity tags must be positive"))
        push!(result,(dimension,entity_tag))
    end
    return result
end

@inline function _model_removal_exists(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int)
    return _model_entity_known(m,dimension,tag) &&
           !((dimension,tag) in removed)
end

function _model_removal_surface_curves(
    m::GeoModel,surface::Int,caller::AbstractString)
    loops=get(m.surfaces,surface,nothing)
    loops===nothing && throw(ErrorException(
        "$caller: Surface[$surface] disappeared during removal planning"))
    curves=Int[]
    for loop in loops
        signed_curves=get(m.loops,loop,nothing)
        signed_curves===nothing && throw(ErrorException(
            "$caller: Surface[$surface] references missing Loop[$loop]; " *
            "rebuild the model"))
        for signed_curve in signed_curves
            curve=abs(signed_curve)
            haskey(m.curves,curve) || throw(ErrorException(
                "$caller: Loop[$loop] references missing Curve[$curve]; " *
                "rebuild the model"))
            push!(curves,curve)
        end
    end
    return curves
end

function _model_removal_volume_surfaces(
    m::GeoModel,volume::Int,caller::AbstractString)
    shells=get(m.volumes,volume,nothing)
    shells===nothing && throw(ErrorException(
        "$caller: Volume[$volume] disappeared during removal planning"))
    surfaces=Int[]
    for shell in shells
        signed_surfaces=get(m.surface_loops,shell,nothing)
        signed_surfaces===nothing && throw(ErrorException(
            "$caller: Volume[$volume] references missing Surface Loop[$shell]; " *
            "rebuild the model"))
        for signed_surface in signed_surfaces
            surface=abs(signed_surface)
            haskey(m.surfaces,surface) || throw(ErrorException(
                "$caller: Surface Loop[$shell] references missing " *
                "Surface[$surface]; rebuild the model"))
            push!(surfaces,surface)
        end
    end
    return surfaces
end

function _model_removal_has_boundary_owner(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int,
    caller::AbstractString)
    # A discrete entity's declared boundary keeps its boundary entities owned.
    for ((entity_dimension,entity_tag),entity) in m.discrete
        _model_removal_exists(m,removed,entity_dimension,entity_tag) ||
            continue
        (dimension,tag) in entity.boundary && return true
    end
    if dimension==0
        # A point is owned by every live curve that wires it — as an endpoint
        # or as an arc control point (center, major-axis point) — and by every
        # surviving `In Sphere`/`Using Point` surface constraint.
        for (curve,endpoints) in m.curves
            _model_removal_exists(m,removed,1,curve) || continue
            tag in endpoints && return true
            tag in get(m.curve_control_points,curve,Int[]) && return true
        end
        for (surface,geom) in m.surface_geometry
            _model_removal_exists(m,removed,2,surface) || continue
            hasproperty(geom,:sphere_center) && geom.sphere_center==tag &&
                return true
        end
        return false
    elseif dimension==1
        for surface in keys(m.surfaces)
            _model_removal_exists(m,removed,2,surface) || continue
            tag in _model_removal_surface_curves(m,surface,caller) && return true
        end
    elseif dimension==2
        for volume in keys(m.volumes)
            _model_removal_exists(m,removed,3,volume) || continue
            tag in _model_removal_volume_surfaces(m,volume,caller) && return true
        end
    end
    return false
end

function _model_removal_has_embedding_owner(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int,
    caller::AbstractString)
    source=(dimension,tag)
    for (target,sources) in m.embeds
        target in removed && continue
        _model_removal_exists(m,removed,target...) || throw(ErrorException(
            "$caller: embedding target $target does not exist; rebuild the model"))
        source in sources && return true
    end
    return false
end

function _model_removal_boundary(
    m::GeoModel,dimension::Int,tag::Int,caller::AbstractString)
    discrete=get(m.discrete,(dimension,tag),nothing)
    discrete!==nothing && return Tuple{Int,Int}[
        (bdim,btag) for (bdim,btag) in discrete.boundary]
    if dimension==0
        return Tuple{Int,Int}[]
    elseif dimension==1
        endpoints=get(m.curves,tag,nothing)
        endpoints===nothing && throw(ErrorException(
            "$caller: Curve[$tag] disappeared during removal planning"))
        # The boundary is every vertex the curve wires — endpoints plus arc
        # control points (center, major-axis point).
        owned=unique!(Int[endpoints...,get(m.curve_control_points,tag,Int[])...])
        for point in owned
            haskey(m.points,point) || throw(ErrorException(
                "$caller: Curve[$tag] references missing Point[$point]; " *
                "rebuild the model"))
        end
        return Tuple{Int,Int}[(0,point) for point in owned]
    elseif dimension==2
        return Tuple{Int,Int}[
            (1,curve) for curve in
            _model_removal_surface_curves(m,tag,caller)]
    end
    return Tuple{Int,Int}[
        (2,surface) for surface in
        _model_removal_volume_surfaces(m,tag,caller)]
end

function _model_removal_plan(
    m::GeoModel,requested::Vector{Tuple{Int,Int}},recursive::Bool,
    caller::AbstractString)
    removed=Set{Tuple{Int,Int}}()
    function attempt!(dimension::Int,tag::Int,recurse::Bool)
        _model_removal_exists(m,removed,dimension,tag) || return nothing
        _model_removal_has_boundary_owner(
            m,removed,dimension,tag,caller) && return nothing
        _model_removal_has_embedding_owner(
            m,removed,dimension,tag,caller) && return nothing
        boundary=recurse ?
            _model_removal_boundary(m,dimension,tag,caller) :
            Tuple{Int,Int}[]
        push!(removed,(dimension,tag))
        for (boundary_dimension,boundary_tag) in boundary
            attempt!(boundary_dimension,boundary_tag,true)
        end
        return nothing
    end
    for (dimension,tag) in requested
        attempt!(dimension,tag,recursive)
    end
    return removed
end

function _model_removal_state(
    m::GeoModel,removed::Set{Tuple{Int,Int}};
    keep_dangling_loops::Bool=false,
    keep_stale_physical::Bool=false)
    removed_tags=ntuple(dimension->Set(
        tag for (entity_dimension,tag) in removed
        if entity_dimension==dimension-1),4)

    points=copy(m.points)
    point_size=copy(m.point_size)
    curves=copy(m.curves)
    curve_control_points=copy(m.curve_control_points)
    curve_types=copy(m.curve_types)
    curve_geometry=copy(m.curve_geometry)
    surfaces=copy(m.surfaces)
    surface_types=copy(m.surface_types)
    surface_geometry=copy(m.surface_geometry)
    volumes=copy(m.volumes)
    for tag in removed_tags[1]
        delete!(points,tag);delete!(point_size,tag)
    end
    for tag in removed_tags[2]
        delete!(curves,tag)
        delete!(curve_control_points,tag)
        delete!(curve_types,tag)
        delete!(curve_geometry,tag)
    end
    for tag in removed_tags[3]
        delete!(surfaces,tag)
        delete!(surface_types,tag)
        delete!(surface_geometry,tag)
    end
    for tag in removed_tags[4]
        delete!(volumes,tag)
    end

    loops=copy(m.loops)
    if !keep_dangling_loops
        for (loop,signed_curves) in m.loops
            any(curve->abs(curve) in removed_tags[2],signed_curves) &&
                delete!(loops,loop)
        end
    end
    surface_loops=copy(m.surface_loops)
    if !keep_dangling_loops
        for (loop,signed_surfaces) in m.surface_loops
            any(surface->abs(surface) in removed_tags[3],signed_surfaces) &&
                delete!(surface_loops,loop)
        end
        for (surface,surface_loops_used) in surfaces
            all(loop->haskey(loops,loop),surface_loops_used) || throw(ErrorException(
                "remove_entities!: surviving Surface[$surface] lost a Curve Loop; " *
                "rebuild the model"))
        end
        for (volume,shells) in volumes
            all(shell->haskey(surface_loops,shell),shells) || throw(ErrorException(
                "remove_entities!: surviving Volume[$volume] lost a Surface Loop; " *
                "rebuild the model"))
        end
    end

    entity_names=copy(m.entity_names)
    entity_visibility=copy(m.entity_visibility)
    entity_colors=copy(m.entity_colors)
    for entity in removed
        delete!(entity_names,entity)
        delete!(entity_visibility,entity)
        delete!(entity_colors,entity)
    end

    physical=Dict{Tuple{Int,Int},Vector{Int}}()
    physical_names=copy(m.physical_names)
    if keep_stale_physical
        # `.geo` Delete keeps the physical-group records verbatim: Gmsh's
        # internals retain the stale member integers, so a member tag that is
        # deleted and later re-created resurrects its membership. Query paths
        # filter to live entities; emptied groups stay (with their names).
        for (key,members) in m.physical
            physical[key]=copy(members)
        end
    else
        for (key,members) in m.physical
            retained=Int[member for member in members
                         if !((key[1],member) in removed)]
            if isempty(retained)
                delete!(physical_names,key)
            else
                physical[key]=retained
            end
        end
    end

    embeds=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (target,sources) in m.embeds
        target in removed && continue
        retained=NTuple{2,Int}[source for source in sources
                               if !(source in removed)]
        isempty(retained) || (embeds[target]=retained)
    end

    periodic=Dict{Tuple{Int,Int},ModelPeriodicConstraint}()
    for (key,constraint) in m.periodic
        slave=(constraint.dim,Int(constraint.slave_entity))
        master=(constraint.dim,Int(constraint.master_entity))
        (slave in removed || master in removed) && continue
        periodic[key]=constraint
    end

    box_extents=copy(m.box_extents)
    cylinders=copy(m.cylinders)
    spheres=copy(m.spheres)
    cones=copy(m.cones)
    booleans=copy(m.booleans)
    boolean_operands=copy(m.boolean_operands)
    boolean_components=copy(m.boolean_components)
    for tag in removed_tags[4]
        for encoding in (box_extents,cylinders,spheres,cones,
                         booleans,boolean_operands,boolean_components)
            delete!(encoding,tag)
        end
    end
    # a removed primary's component volumes become self-linked — their own
    # shells and operand snapshots carry everything the link needs
    for (k,v) in collect(boolean_components)
        v in removed_tags[4] && (boolean_components[k]=k)
    end

    discrete=Dict{Tuple{Int,Int},DiscreteEntity}()
    for (key,entity) in m.discrete
        key in removed && continue
        retained=DiscreteEntity(
            NTuple{2,Int}[b for b in entity.boundary if !(b in removed)],
            copy(entity.node_tags),copy(entity.node_coords),
            copy(entity.node_params),copy(entity.element_types),
            copy(entity.element_tags),
            Vector{Int32}[copy(nodes) for nodes in entity.element_nodes],
            Dict{Int32,Vector{Float64}}(k=>copy(v) for (k,v) in
                                        entity.aux_params))
        discrete[key]=retained
    end

    meshing=m.meshing
    for tag in removed_tags[2]
        delete!(meshing.transfinite_curves,tag)
        delete!(meshing.degenerated,tag)
    end
    for tag in removed_tags[3]
        delete!(meshing.transfinite_surfaces,tag)
    end
    for tag in removed_tags[4]
        delete!(meshing.transfinite_volumes,tag)
        delete!(meshing.outward_orientation,tag)
        delete!(meshing.quad_tri,tag)
    end
    for entity in removed
        for store in (meshing.recombine,meshing.smoothing,meshing.reverse,
                      meshing.algorithm,meshing.size_at_params,
                      meshing.size_from_boundary,meshing.attached,
                      meshing.extrude)
            delete!(store,entity)
        end
    end
    filter!(compound->!any(member->(compound.first,member) in removed,
                           compound.second),meshing.compounds)

    return (;points,point_size,curves,curve_control_points,curve_types,
            curve_geometry,loops,surfaces,surface_types,surface_geometry,
            surface_loops,volumes,
            entity_names,entity_visibility,entity_colors,physical,physical_names,
            box_extents,cylinders,spheres,cones,booleans,boolean_operands,
            boolean_components,periodic,embeds,discrete)
end

"""
    remove_entities!(model, dim_tags, recursive=false) -> Int

Atomically process ordered positive `(dimension, entity_tag)` removals and return
the number removed. Missing entities and entities still used as an explicit boundary
or embedding source are skipped. With `recursive=true`, each removed entity's
explicit boundary is processed down to Points; embedded entities are not boundaries.

Entity names, visibility, colors, Physical memberships, embedding targets, periodic
relations, native solid encodings, and Boolean-result snapshots owned by removed
entities are cleaned up. Empty Physical groups and their names are removed,
construction loops that would dangle are discarded, and automatic tag counters
remain monotonic. Materialized Volumes — `add_box!` shells and OCC Boolean
results alike — carry their boundary entities like any explicit shell Volume,
so recursive removal descends through them; unmaterialized primitives have no
boundary entities, so recursive removal stops at the Volume.
"""
function remove_entities!(m::GeoModel,dim_tags,recursive=false)
    caller="remove_entities!"
    recursive isa Bool || throw(ArgumentError(
        "$caller: recursive must be Bool"))
    requested=_model_removal_dim_tags(dim_tags,caller)
    removed=_model_removal_plan(m,requested,recursive,caller)
    isempty(removed) && return 0
    state=_model_removal_state(m,removed)

    m.points=state.points
    m.point_size=state.point_size
    m.curves=state.curves
    m.curve_control_points=state.curve_control_points
    m.curve_types=state.curve_types
    m.curve_geometry=state.curve_geometry
    m.loops=state.loops
    m.surfaces=state.surfaces
    m.surface_types=state.surface_types
    m.surface_geometry=state.surface_geometry
    m.surface_loops=state.surface_loops
    m.volumes=state.volumes
    m.entity_names=state.entity_names
    m.entity_visibility=state.entity_visibility
    m.entity_colors=state.entity_colors
    m.physical=state.physical
    m.physical_names=state.physical_names
    m.box_extents=state.box_extents
    m.cylinders=state.cylinders
    m.spheres=state.spheres
    m.cones=state.cones
    m.booleans=state.booleans
    m.boolean_operands=state.boolean_operands
    m.boolean_components=state.boolean_components
    m.periodic=state.periodic
    m.embeds=state.embeds
    m.discrete=state.discrete
    return length(removed)
end

# Volume removal for consumed operands (Boolean `Delete`, API cut/fuse/common):
# recursive so a materialized primitive's boundary entities go with it, while
# the ownership guards keep children still referenced by surviving parents.
function _remove_volume_entity!(m::GeoModel,tag::Int)
    return remove_entities!(m,[(3,tag)],true)>0
end

# ---------------------------------------------------------------------------
# `.geo` `Delete { ListOfShapes }` / `Recursive Delete { ... }` — the Gmsh
# `GEO_Internals::remove` semantics, which differ from `remove_entities!`:
#
#  * sign handling — vertices compare by `abs(tag)` (`CompareVertex`) and a
#    curve deletion always tries both signs (`DeleteCurve(tag)` +
#    `DeleteCurve(-tag)` for the stored reversed mirror), so `Point{-7}` and
#    `Curve{-9}` still delete their entities; surfaces and volumes compare
#    signed, so `Surface{-1}`/`Volume{-1}` silently match nothing.
#  * refusal is boundary ownership only — a point wired into a live curve
#    (endpoint or arc control point), a curve on a live surface's loops, a
#    surface on a live volume's shells. Embedding sources and sphere-center
#    references do NOT block deletion (the embed record dies with the point;
#    a sphere keeps its resolved center coordinates).
#  * `Recursive` pre-collects the whole transitive boundary into sorted sets
#    per dimension and then deletes top-down (surfaces, curves, points), each
#    child still subject to the ownership refusal against remaining entities.
#  * each successful deletion decrements that dimension's tag counter by one
#    when the deleted tag was the counter (`if(tag == max) max--`), applied in
#    deletion order — it never recomputes the counter from live entities.
#  * curve/surface loop records survive dangling (there is no `Delete Loop`
#    syntax; a later `Plane Surface` on such a loop errors), and physical
#    groups keep their stale member integers, so a re-created tag resurrects
#    its membership. Physical names are never dropped by entity deletion.
#
# Each requested entry is attempted independently in list order; there is no
# atomicity and no error for missing or refused entities.

@inline function _geo_removal_lookup(
    m::GeoModel,removed::Set{Tuple{Int,Int}},dimension::Int,tag::Int)
    if dimension<=1
        tag=abs(tag)
    elseif tag<=0
        return false
    end
    tag==0 && return false
    (dimension,tag) in removed && return false
    return haskey(_model_entity_dictionary(m,dimension),tag)
end

function _geo_removal_owned(m::GeoModel,removed::Set{Tuple{Int,Int}},
                            dimension::Int,tag::Int)
    if dimension==0
        for (curve,endpoints) in m.curves
            _geo_removal_lookup(m,removed,1,curve) || continue
            tag in endpoints && return true
            tag in get(m.curve_control_points,curve,Int[]) && return true
        end
    elseif dimension==1
        for (surface,loops) in m.surfaces
            _geo_removal_lookup(m,removed,2,surface) || continue
            for loop in loops
                signed_curves=get(m.loops,loop,nothing)
                signed_curves===nothing && continue
                any(curve->abs(curve)==tag,signed_curves) && return true
            end
        end
    elseif dimension==2
        for (volume,shells) in m.volumes
            _geo_removal_lookup(m,removed,3,volume) || continue
            for shell in shells
                signed_surfaces=get(m.surface_loops,shell,nothing)
                signed_surfaces===nothing && continue
                any(surface->abs(surface)==tag,signed_surfaces) && return true
            end
        end
    end
    return false
end

function _geo_cascade_collect_surface!(m::GeoModel,surface_tag::Int,
                                       curves::Vector{Int},points::Vector{Int})
    loops=get(m.surfaces,surface_tag,nothing)
    loops===nothing && return nothing
    for loop in loops
        signed_curves=get(m.loops,loop,nothing)
        signed_curves===nothing && continue
        for signed_curve in signed_curves
            curve=abs(signed_curve)
            push!(curves,curve)
            endpoints=get(m.curves,curve,nothing)
            endpoints===nothing && continue
            append!(points,endpoints)
            append!(points,get(m.curve_control_points,curve,Int[]))
        end
    end
    return nothing
end

# Pre-collected recursive boundary, matching Gmsh's `DeleteCurve`/
# `DeleteSurface`/`DeleteVolume`: all boundary entities of the removed entity,
# gathered transitively into per-dimension sorted sets and attempted
# descending-dimension first (surfaces, then curves, then points).
function _geo_removal_cascade(m::GeoModel,dimension::Int,tag::Int)
    dimension==0 && return Tuple{Int,Int}[]
    if dimension==1
        endpoints=get(m.curves,tag,nothing)
        endpoints===nothing && return Tuple{Int,Int}[]
        points=unique!(Int[endpoints...,
                           get(m.curve_control_points,tag,Int[])...])
        sort!(points)
        return Tuple{Int,Int}[(0,point) for point in points]
    end
    surfaces=Int[];curves=Int[];points=Int[]
    if dimension==2
        _geo_cascade_collect_surface!(m,tag,curves,points)
    else
        shells=get(m.volumes,tag,nothing)
        shells===nothing && return Tuple{Int,Int}[]
        for shell in shells
            signed_surfaces=get(m.surface_loops,shell,nothing)
            signed_surfaces===nothing && continue
            for signed_surface in signed_surfaces
                surface=abs(signed_surface)
                push!(surfaces,surface)
                _geo_cascade_collect_surface!(m,surface,curves,points)
            end
        end
        sort!(unique!(surfaces))
    end
    sort!(unique!(curves));sort!(unique!(points))
    result=Tuple{Int,Int}[]
    for surface in surfaces
        push!(result,(2,surface))
    end
    for curve in curves
        push!(result,(1,curve))
    end
    for point in points
        push!(result,(0,point))
    end
    return result
end

"""
    _geo_delete_entities!(model, dim_tags; recursive=false) -> removed

`.geo` `Delete`/`Recursive Delete` — per-entity best-effort removal under
Gmsh's `GEO_Internals::remove` contract (see the module notes above).
`dim_tags` may carry signed tags. Returns the removed `(dim, tag)` pairs in
deletion order so the caller can update tag-allocator counters.
"""
function _geo_delete_entities!(m::GeoModel,dim_tags;recursive::Bool=false)
    removed=Tuple{Int,Int}[]
    inset=Set{Tuple{Int,Int}}()
    function attempt!(dimension::Int,tag::Int)
        _geo_removal_lookup(m,inset,dimension,tag) || return nothing
        key=(dimension,abs(tag))
        _geo_removal_owned(m,inset,dimension,key[2]) && return nothing
        push!(removed,key);push!(inset,key)
        if recursive
            # Gmsh calls the non-recursive deleters on the pre-collected
            # transitive boundary — children are attempted flat.
            for (child_dimension,child_tag) in
                _geo_removal_cascade(m,dimension,key[2])
                _geo_removal_lookup(m,inset,child_dimension,child_tag) ||
                    continue
                _geo_removal_owned(m,inset,child_dimension,child_tag) &&
                    continue
                child_key=(child_dimension,child_tag)
                push!(removed,child_key);push!(inset,child_key)
            end
        end
        return nothing
    end
    for entry in dim_tags
        (dimension,tag)=entry
        attempt!(Int(dimension),Int(tag))
    end
    isempty(removed) && return removed
    state=_model_removal_state(m,inset;keep_dangling_loops=true,
                               keep_stale_physical=true)

    m.points=state.points
    m.point_size=state.point_size
    m.curves=state.curves
    m.curve_control_points=state.curve_control_points
    m.curve_types=state.curve_types
    m.curve_geometry=state.curve_geometry
    m.surfaces=state.surfaces
    m.surface_types=state.surface_types
    m.surface_geometry=state.surface_geometry
    m.volumes=state.volumes
    m.entity_names=state.entity_names
    m.entity_visibility=state.entity_visibility
    m.entity_colors=state.entity_colors
    m.physical=state.physical
    m.physical_names=state.physical_names
    m.box_extents=state.box_extents
    m.cylinders=state.cylinders
    m.spheres=state.spheres
    m.cones=state.cones
    m.booleans=state.booleans
    m.boolean_operands=state.boolean_operands
    m.boolean_components=state.boolean_components
    m.periodic=state.periodic
    m.embeds=state.embeds
    m.discrete=state.discrete
    # loops/surface_loops intentionally untouched — records survive dangling.
    # Counter decrement happens per removal, in order (`if(tag==max) max--`).
    for (dimension,tag) in removed
        tag==m.next_tag[dimension+1] && (m.next_tag[dimension+1]-=1)
    end
    return removed
end

"""
    _geo_reset_model_geometry!(model)

`.geo` `Delete Model` — `GEO_Internals::destroy()` plus `GModel::destroy`:
every entity, loop, physical group, and per-entity meshing attribute is
destroyed and all tag counters reset. Like Gmsh, the compound-mesh multimap
(`m.meshing.compounds`) and the physical-group NAME table
(`m.physical_names`) survive — they live outside the freed entity records.
User variables and the current factory are exec-layer state and are untouched
here.
"""
function _geo_reset_model_geometry!(m::GeoModel)
    empty!(m.points);empty!(m.point_size)
    empty!(m.curves);empty!(m.curve_control_points);empty!(m.curve_types)
    empty!(m.curve_geometry)
    empty!(m.loops)
    empty!(m.surfaces);empty!(m.surface_types);empty!(m.surface_geometry)
    empty!(m.surface_loops)
    empty!(m.volumes)
    empty!(m.entity_names);empty!(m.entity_visibility);empty!(m.entity_colors)
    empty!(m.physical)
    empty!(m.box_extents);empty!(m.cylinders);empty!(m.spheres);empty!(m.cones)
    empty!(m.booleans);empty!(m.boolean_operands);empty!(m.boolean_components)
    empty!(m.periodic);empty!(m.embeds);empty!(m.discrete)
    meshing=m.meshing
    for store in (meshing.transfinite_curves,meshing.degenerated)
        empty!(store)
    end
    empty!(meshing.transfinite_surfaces)
    for store in (meshing.transfinite_volumes,meshing.outward_orientation,
                  meshing.quad_tri)
        empty!(store)
    end
    for store in (meshing.recombine,meshing.smoothing,meshing.reverse,
                  meshing.algorithm,meshing.size_at_params,
                  meshing.size_from_boundary,meshing.attached,meshing.extrude)
        empty!(store)
    end
    empty!(meshing.homology_requests)
    m.next_tag .= 0
    m.physical_tag_max=0
    return nothing
end
