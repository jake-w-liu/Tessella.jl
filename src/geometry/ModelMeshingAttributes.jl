# Per-entity meshing attributes and discrete entities on GeoModel — the model
# half of the Gmsh 4.15.2 `model.mesh` attribute surface. The API wrappers in
# interfaces/API.jl validate and forward here; generation reads these records.

export set_transfinite_curve!, set_transfinite_surface!, set_transfinite_volume!,
       set_transfinite_automatic!, set_recombine!, set_smoothing!, set_reverse!,
       set_algorithm!, set_size_at_parametric_points!, set_size_from_boundary!,
       set_size_callback!, set_compound!, set_outward_orientation!, set_order!,
       set_transfinite_tri!,
       remove_constraints!, add_discrete_entity!, add_discrete_nodes!,
       add_discrete_elements!, add_homology_request!, clear_homology_requests!,
       classify_surfaces!, compute_homology!, create_geometry!,
       create_topology!,
       model_discrete_entity,
       model_discrete_entities

function _mesh_attr_entity!(m::GeoModel,dim::Int,tag::Int,caller::AbstractString;
                            allow_discrete::Bool=false)
    haskey(m.discrete,(dim,tag)) && return allow_discrete ? nothing :
        throw(ArgumentError(
            "$caller: ($dim,$tag) is a discrete entity with no parametrization"))
    entities=_model_entity_dictionary(m,dim)
    haskey(entities,tag) || throw(ArgumentError(
        "$caller: unknown entity ($dim,$tag)"))
    return nothing
end

function _mesh_attr_positive_int(value,caller::AbstractString,name::AbstractString)
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be an integer"))
    result=try Int(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name exceeds the platform Int range"))
    end
    return result
end

function _mesh_attr_finite(value,caller::AbstractString,name::AbstractString)
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    value isa Real || throw(ArgumentError("$caller: $name must be a real number"))
    result=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    isfinite(result) || throw(ArgumentError("$caller: $name must be finite"))
    return result
end

# Gmsh's meshType table (gmsh.cpp `setTransfiniteCurve`): `Power` aliases
# `Progression`, and the `_HWall` kinds interpret `coef` as a signed first-layer
# wall height instead of a distribution coefficient.
const _TRANSFINITE_CURVE_KINDS = Dict{String,Symbol}(
    "progression"=>:progression, "power"=>:progression,
    "bump"=>:bump, "beta"=>:beta,
    "progression_hwall"=>:progression_hwall,
    "bump_hwall"=>:bump_hwall, "beta_hwall"=>:beta_hwall,
    # Grammar-level types 8/9 — Gmsh's `F_Transfinite` has no case for them, so
    # they hit the unknown-type warning and produce a uniform distribution.
    "beta_symmetrical"=>:beta_symmetrical,
    "beta_symmetrical_hwall"=>:beta_symmetrical_hwall)
const _TRANSFINITE_HWALL_KINDS = (:progression_hwall,:bump_hwall,:beta_hwall)
const _TRANSFINITE_ARRANGEMENTS = Dict{String,Symbol}(
    "left"=>:left, "right"=>:right,
    # Gmsh's TransfiniteArrangement maps bare `Alternate` to AlternateRight.
    "alternate"=>:alternate_right,
    "alternateleft"=>:alternate_left, "alternateright"=>:alternate_right)

"""
    set_transfinite_curve!(model, tag, num_nodes, mesh_type="Progression", coef=1.0)

Record a transfinite meshing constraint on `Curve[tag]` — `num_nodes` nodes
distributed by `mesh_type` (`"Progression"`/`"Power"`, `"Bump"`, `"Beta"`, or
the `"Progression_HWall"`/`"Bump_HWall"`/`"Beta_HWall"` wall-height variants,
case-insensitive) with `coef`, matching Gmsh's `setTransfiniteCurve`. For the
ordinary kinds a negative `coef` flips the distribution direction, as upstream;
for the HWall kinds `coef` is the signed first-layer wall height.
"""
function set_transfinite_curve!(m::GeoModel,tag,num_nodes,mesh_type="Progression",
                                coef=1.0)
    caller="set_transfinite_curve!"
    t=_tag(tag,caller,1)
    _mesh_attr_entity!(m,1,t,caller)
    count=_mesh_attr_positive_int(num_nodes,caller,"num_nodes")
    count>=2 || throw(ArgumentError(
        "$caller: num_nodes must be at least 2 (got $count)"))
    mesh_type isa AbstractString || throw(ArgumentError(
        "$caller: mesh_type must be a string"))
    kind=get(_TRANSFINITE_CURVE_KINDS,lowercase(String(mesh_type)),nothing)
    kind===nothing && throw(ArgumentError(
        "$caller: mesh_type must be Progression/Power, Bump, Beta, or a " *
        "*_HWall variant (got \"$mesh_type\")"))
    coefficient=_mesh_attr_finite(coef,caller,"coef")
    # Gmsh's API stores `abs(coef)` for the ordinary types and negates the
    # (signed) type when `coef < 0` — a negative coefficient means a reversed
    # distribution. HWall records keep the signed wall height.
    reversed=false
    if !(kind in _TRANSFINITE_HWALL_KINDS)
        reversed=coefficient<0
        coefficient=abs(coefficient)
    end
    m.meshing.transfinite_curves[t]=(num_nodes=count,kind=kind,coef=coefficient,
                                     reversed=reversed)
    return nothing
end

"""
    set_transfinite_surface!(model, tag, arrangement="Left", corners=Int[])

Record a transfinite constraint on `Surface[tag]`, matching Gmsh's
`setTransfiniteSurface`. `arrangement` is `"Left"`, `"Right"`,
`"AlternateLeft"`, or `"AlternateRight"`; `corners` optionally pins the corner
Point tags (3 or 4 entries, empty = automatic from the boundary loop).
"""
function set_transfinite_surface!(m::GeoModel,tag,arrangement="Left",
                                  corners=Int[])
    caller="set_transfinite_surface!"
    t=_tag(tag,caller,2)
    _mesh_attr_entity!(m,2,t,caller)
    arrangement isa AbstractString || throw(ArgumentError(
        "$caller: arrangement must be a string"))
    mode=get(_TRANSFINITE_ARRANGEMENTS,lowercase(String(arrangement)),nothing)
    mode===nothing && throw(ArgumentError(
        "$caller: arrangement must be Left, Right, AlternateLeft, or " *
        "AlternateRight (got \"$arrangement\")"))
    (corners isa AbstractVector || corners isa Tuple) || throw(ArgumentError(
        "$caller: corners must be a vector or tuple of Point tags"))
    corner_tags=Int[_mesh_attr_positive_int(c,caller,"corners entry")
                    for c in corners]
    if !isempty(corner_tags)
        length(corner_tags) in (3,4) || throw(ArgumentError(
            "$caller: corners must list 3 or 4 Point tags"))
        for point_tag in corner_tags
            haskey(m.points,point_tag) || throw(ArgumentError(
                "$caller: unknown corner Point[$point_tag]"))
        end
    end
    m.meshing.transfinite_surfaces[t]=(arrangement=mode,corners=corner_tags)
    return nothing
end

"""
    set_transfinite_volume!(model, tag, corners=Int[])

Record a transfinite constraint on `Volume[tag]`, matching Gmsh's
`setTransfiniteVolume`. `corners` optionally pins the 6 or 8 corner Point tags.
"""
function set_transfinite_volume!(m::GeoModel,tag,corners=Int[])
    caller="set_transfinite_volume!"
    t=_tag(tag,caller,3)
    _mesh_attr_entity!(m,3,t,caller)
    (corners isa AbstractVector || corners isa Tuple) || throw(ArgumentError(
        "$caller: corners must be a vector or tuple of Point tags"))
    corner_tags=Int[_mesh_attr_positive_int(c,caller,"corners entry")
                    for c in corners]
    if !isempty(corner_tags)
        length(corner_tags) in (5,6,8) || throw(ArgumentError(
            "$caller: corners must list 5, 6, or 8 Point tags"))
        for point_tag in corner_tags
            haskey(m.points,point_tag) || throw(ArgumentError(
                "$caller: unknown corner Point[$point_tag]"))
        end
    end
    m.meshing.transfinite_volumes[t]=corner_tags
    return nothing
end

"""
    set_recombine!(model, dim, tag, angle=45.0)

Record a recombination constraint on entity `(dim, tag)` — triangles in the
generated mesh are recombined into quadrangles up to `angle` degrees, matching
Gmsh's `setRecombine`.
"""
function set_recombine!(m::GeoModel,dim,tag,angle=45.0)
    caller="set_recombine!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    dimension in (1,2,3) || throw(ArgumentError(
        "$caller: recombination applies to curves, surfaces, or volumes"))
    _mesh_attr_entity!(m,dimension,t,caller;allow_discrete=true)
    value=_mesh_attr_finite(angle,caller,"angle")
    0<=value<=90 || throw(ArgumentError(
        "$caller: angle must be in [0, 90] degrees (got $angle)"))
    m.meshing.recombine[(dimension,t)]=value
    return nothing
end

"""
    set_smoothing!(model, dim, tag, val)

Record `val` Laplace-smoother iterations applied to the mesh of entity
`(dim, tag)` during generation, matching Gmsh's `setSmoothing`.
"""
function set_smoothing!(m::GeoModel,dim,tag,val)
    caller="set_smoothing!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    _mesh_attr_entity!(m,dimension,t,caller)
    count=_mesh_attr_positive_int(val,caller,"val")
    m.meshing.smoothing[(dimension,t)]=count
    return nothing
end

"""
    set_reverse!(model, dim, tag, val=true)

Record a reverse-orientation constraint on entity `(dim, tag)`, matching
Gmsh's `setReverse`. Generated element orientations on that entity are flipped
with respect to the natural orientation when `val` is true.
"""
function set_reverse!(m::GeoModel,dim,tag,val=true)
    caller="set_reverse!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    dimension in (1,2,3) || throw(ArgumentError(
        "$caller: reverse applies to curves, surfaces, or volumes"))
    val isa Bool || throw(ArgumentError("$caller: val must be Bool"))
    _mesh_attr_entity!(m,dimension,t,caller)
    if val
        m.meshing.reverse[(dimension,t)]=true
    else
        delete!(m.meshing.reverse,(dimension,t))
    end
    return nothing
end

# The native surface generator is a constrained-Delaunay refinement pipeline —
# Gmsh's algorithm ids 5 (Delaunay) and 6 (Frontal-Delaunay) both select it.
const _SUPPORTED_SURFACE_ALGORITHMS = (5,6)
const _SUPPORTED_VOLUME_ALGORITHMS = (1,)

"""
    set_algorithm!(model, dim, tag, val)

Record the meshing algorithm for entity `(dim, tag)`, matching Gmsh's
`setAlgorithm`. Only `dim == 2` and `dim == 3` are accepted; unsupported
algorithm numbers are stored but fail explicitly at `generate` time, matching
Gmsh's deferred-validation behavior.
"""
function set_algorithm!(m::GeoModel,dim,tag,val)
    caller="set_algorithm!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    dimension in (2,3) || throw(ArgumentError(
        "$caller: set_algorithm is only supported for dim 2 or 3"))
    _mesh_attr_entity!(m,dimension,t,caller)
    number=_mesh_attr_positive_int(val,caller,"val")
    m.meshing.algorithm[(dimension,t)]=number
    return nothing
end

"""
    set_size_at_parametric_points!(model, dim, tag, parametric_coord, sizes)

Record mesh-size constraints at parametric points on entity `(dim, tag)`,
matching Gmsh's `setSizeAtParametricPoints`. Only `dim == 1` is supported, as
upstream: `parametric_coord` lists curve parameters in `[0, 1]` and `sizes` the
matching target sizes — the generator inserts a boundary node at each location
with the recorded size.
"""
function set_size_at_parametric_points!(m::GeoModel,dim,tag,parametric_coord,
                                        sizes)
    caller="set_size_at_parametric_points!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    dimension==1 || throw(ArgumentError(
        "$caller: only dim 1 (curves) is supported"))
    _mesh_attr_entity!(m,dimension,t,caller)
    (parametric_coord isa AbstractVector || parametric_coord isa Tuple) ||
        throw(ArgumentError(
            "$caller: parametric_coord must be a vector or tuple of parameters"))
    (sizes isa AbstractVector || sizes isa Tuple) || throw(ArgumentError(
        "$caller: sizes must be a vector or tuple of real numbers"))
    length(parametric_coord)==length(sizes) || throw(ArgumentError(
        "$caller: parametric_coord and sizes length mismatch " *
        "($(length(parametric_coord)) vs $(length(sizes)))"))
    entries=Tuple{Float64,Float64}[]
    for i in eachindex(parametric_coord)
        parameter=_mesh_attr_finite(parametric_coord[i],caller,
                                    "parametric_coord entry")
        0<=parameter<=1 || throw(ArgumentError(
            "$caller: parametric_coord entry $parameter is outside [0, 1]"))
        size_value=_mesh_attr_finite(sizes[i],caller,"sizes entry")
        size_value>0 || throw(ArgumentError(
            "$caller: sizes entry $size_value must be positive"))
        push!(entries,(parameter,size_value))
    end
    sort!(entries;by=entry->entry[1])
    m.meshing.size_at_params[(1,t)]=[(Float64[p],s) for (p,s) in entries]
    return nothing
end

"""
    set_size_from_boundary!(model, dim, tag, val)

Record whether the interior mesh size of entity `(dim, tag)` extends its
boundary sizes rather than interpolating them, matching Gmsh's
`setSizeFromBoundary`. Only `dim == 2` is supported, as upstream.
"""
function set_size_from_boundary!(m::GeoModel,dim,tag,val)
    caller="set_size_from_boundary!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    dimension==2 || throw(ArgumentError(
        "$caller: only dim 2 (surfaces) is supported"))
    val isa Integer || val isa Bool || throw(ArgumentError(
        "$caller: val must be an integer flag"))
    _mesh_attr_entity!(m,dimension,t,caller)
    if val==0
        delete!(m.meshing.size_from_boundary,(dimension,t))
    else
        m.meshing.size_from_boundary[(dimension,t)]=true
    end
    return nothing
end

"""
    set_size_callback!(model, callback)

Record a model-level mesh-size callback `(dim, tag, x, y, z, lc) -> size`,
matching Gmsh's `setSizeCallback`. The generator calls it wherever it would
otherwise prescribe `lc`; returning `lc` is a no-op. Pass `nothing` to clear.
"""
function set_size_callback!(m::GeoModel,callback)
    caller="set_size_callback!"
    callback===nothing || callback isa Function || throw(ArgumentError(
        "$caller: callback must be a Function or nothing"))
    m.meshing.size_callback=callback
    return nothing
end

"""
    set_compound!(model, dim, tags)

Record a compound meshing constraint — entities `(dim, tags)` are meshed as a
single entity, matching Gmsh's `setCompound`. `dim` must be 1 or 2 and every
listed entity must exist.
"""
function set_compound!(m::GeoModel,dim,tags)
    caller="set_compound!"
    dimension=_dimension(dim,caller)
    dimension in (1,2) || throw(ArgumentError(
        "$caller: only dim 1 (curves) or dim 2 (surfaces) is supported"))
    (tags isa AbstractVector || tags isa Tuple) || throw(ArgumentError(
        "$caller: tags must be a vector or tuple of entity tags"))
    isempty(tags) && throw(ArgumentError("$caller: tags must not be empty"))
    entity_tags=Int[]
    for raw in tags
        t=_mesh_attr_positive_int(raw,caller,"tags entry")
        _mesh_attr_entity!(m,dimension,t,caller)
        push!(entity_tags,t)
    end
    sort!(unique!(entity_tags))
    push!(m.meshing.compounds,dimension=>entity_tags)
    return nothing
end

"""
    set_outward_orientation!(model, tag)

Record the constraint that all boundary surfaces of `Volume[tag]` are oriented
with outward-pointing normals, matching Gmsh's `setOutwardOrientation`. An
existing mesh is reoriented in place.
"""
function set_outward_orientation!(m::GeoModel,tag)
    caller="set_outward_orientation!"
    t=_tag(tag,caller,3)
    _mesh_attr_entity!(m,3,t,caller)
    push!(m.meshing.outward_orientation,t)
    return nothing
end

"""
    remove_constraints!(model, dim_tags=NTuple{2,Int}[])

Clear generation-scoped per-entity meshing attributes — transfinite, recombine,
smoothing, reverse, algorithm, compound, and outward-orientation records, plus
parametric and boundary size channels — matching Gmsh 4.15.2's
`mesh.removeConstraints`. With an empty `dim_tags` every entity's constraints
and the size callback are cleared; otherwise only the listed entities'
records are removed. Periodic relations, embeddings, and Point sizes are
retained.
"""
function remove_constraints!(m::GeoModel,dim_tags=NTuple{2,Int}[])
    caller="remove_constraints!"
    meshing=m.meshing
    if isempty(dim_tags)
        empty!(meshing.transfinite_curves)
        empty!(meshing.transfinite_surfaces)
        empty!(meshing.transfinite_volumes)
        empty!(meshing.recombine)
        empty!(meshing.extrude)
        empty!(meshing.smoothing)
        empty!(meshing.reverse)
        empty!(meshing.algorithm)
        empty!(meshing.size_at_params)
        empty!(meshing.size_from_boundary)
        meshing.size_callback=nothing
        empty!(meshing.compounds)
        empty!(meshing.outward_orientation)
        empty!(meshing.degenerated)
        empty!(meshing.quad_tri)
        return nothing
    end
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dim, tag) pairs"))
    for entry in dim_tags
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif (entry isa Tuple || entry isa AbstractVector) && length(entry)==2
            (entry[1],entry[2])
        else
            throw(ArgumentError(
                "$caller: each dim_tags entry must be a (dim, tag) pair"))
        end
        dimension=_dimension(pair[1],caller)
        tag=_tag(pair[2],caller,dimension)
        _model_entity_known(m,dimension,tag) || throw(ArgumentError(
            "$caller: unknown entity ($dimension,$tag)"))
        dimension==1 && delete!(meshing.transfinite_curves,tag)
        dimension==1 && delete!(meshing.degenerated,tag)
        dimension==2 && delete!(meshing.transfinite_surfaces,tag)
        if dimension==3
            delete!(meshing.transfinite_volumes,tag)
            delete!(meshing.outward_orientation,tag)
            delete!(meshing.quad_tri,tag)
        end
        for store in (meshing.recombine,meshing.smoothing,meshing.reverse,
                      meshing.algorithm,meshing.size_at_params,
                      meshing.size_from_boundary,meshing.extrude)
            delete!(store,(dimension,tag))
        end
        filter!(compound->!(compound.first==dimension && tag in compound.second),
                meshing.compounds)
    end
    return nothing
end

"""
    add_discrete_entity!(model, dim, tag=0, boundary=NTuple{2,Int}[]) -> tag

Create a parametrization-free entity of dimension `dim`, matching Gmsh's
`model.addDiscreteEntity`. `tag=0` allocates the next tag; `boundary` lists
declared `(dim, tag)` boundary entities. The entity coexists with native
entities in queries but has no geometry parametrization.
"""
function add_discrete_entity!(m::GeoModel,dim,tag=0,boundary=NTuple{2,Int}[])
    caller="add_discrete_entity!"
    dimension=_dimension(dim,caller)
    value=_tag(tag,caller,dimension)
    if value==0
        value=m.next_tag[dimension+1]+1
        while haskey(_model_entity_dictionary(m,dimension),value) ||
              haskey(m.discrete,(dimension,value))
            value+=1
        end
        m.next_tag[dimension+1]=value
    end
    (haskey(_model_entity_dictionary(m,dimension),value) ||
     haskey(m.discrete,(dimension,value))) && throw(ArgumentError(
        "$caller: entity ($dimension,$value) already exists"))
    (boundary isa AbstractVector || boundary isa Tuple) || throw(ArgumentError(
        "$caller: boundary must be a vector or tuple of (dim, tag) pairs"))
    pairs=NTuple{2,Int}[]
    for entry in boundary
        pair=if entry isa Pair
            (first(entry),last(entry))
        elseif (entry isa Tuple || entry isa AbstractVector) && length(entry)==2
            (entry[1],entry[2])
        else
            throw(ArgumentError(
                "$caller: each boundary entry must be a (dim, tag) pair"))
        end
        bdim=_dimension(pair[1],caller)
        bdim<dimension || throw(ArgumentError(
            "$caller: boundary dimension $bdim must be below $dimension"))
        btag=_mesh_attr_positive_int(pair[2],caller,"boundary tag")
        _model_entity_known(m,bdim,btag) || throw(ArgumentError(
            "$caller: unknown boundary entity ($bdim,$btag)"))
        push!(pairs,(bdim,btag))
    end
    entity=DiscreteEntity()
    entity.boundary=pairs
    m.discrete[(dimension,value)]=entity
    return value
end

"""Return the [`DiscreteEntity`](@ref) record for `(dim, tag)` or `nothing`."""
model_discrete_entity(m::GeoModel,dim,tag)=get(m.discrete,(Int(dim),Int(tag)),nothing)

"""Sorted `(dim, tag)` pairs of every discrete entity in `model`."""
function model_discrete_entities(m::GeoModel,dim=-1)
    selected=_query_dimension(dim,"model_discrete_entities")
    result=Tuple{Int,Int}[]
    for (dimension,tag) in keys(m.discrete)
        (selected==-1 || selected==dimension) || continue
        push!(result,(dimension,tag))
    end
    return sort!(result)
end

"""
    set_order!(model, order)

Record the global element order for generation, matching Gmsh's `setOrder`.
Orders 1 and 2 are supported; higher values are stored and fail explicitly at
`generate` time.
"""
function set_order!(m::GeoModel,order)
    caller="set_order!"
    value=_mesh_attr_positive_int(order,caller,"order")
    value>=1 || throw(ArgumentError(
        "$caller: order must be at least 1 (got $order)"))
    m.meshing.order=value
    return nothing
end

"""
    set_transfinite_tri!(model, value)

Record the three-sided transfinite surface algorithm for generation, matching
Gmsh's `Mesh.TransfiniteTri` mesh option. `0` (the default) selects the legacy
collapsed-quadrilateral algorithm, which accepts any boundary whose two sides
incident to the collapsed corner carry equal node counts; `1` selects the
compact triangular-lattice algorithm, which requires equal node counts on all
three sides. Other values are rejected.
"""
function set_transfinite_tri!(m::GeoModel,value)
    caller="set_transfinite_tri!"
    value isa Integer || throw(ArgumentError(
        "$caller: transfinite triangle mode must be 0 or 1 (got $value)"))
    value in (0,1) || throw(ArgumentError(
        "$caller: transfinite triangle mode must be 0 or 1 (got $value)"))
    m.meshing.transfinite_tri=Int(value)
    return nothing
end

"""
    set_transfinite_automatic!(model, dim_tags=NTuple{2,Int}[],
                               corner_angle=2.35, recombine=true)

Apply transfinite constraints automatically, matching Gmsh's
`setTransfiniteAutomatic`. For every listed surface (all native surfaces when
`dim_tags` is empty) whose boundary has 3 or 4 corners — vertices where
consecutive boundary segments turn by more than `corner_angle` radians — a
transfinite surface is recorded and matching node counts are written onto the
boundary curves. For every listed 6-faced volume, a transfinite volume is
recorded and its surfaces are processed the same way. When `recombine` is
true, each processed surface also records a 45° recombination constraint.
"""
function set_transfinite_automatic!(m::GeoModel,dim_tags=NTuple{2,Int}[],
                                    corner_angle=2.35,recombine=true)
    caller="set_transfinite_automatic!"
    (dim_tags isa AbstractVector || dim_tags isa Tuple) || throw(ArgumentError(
        "$caller: dim_tags must be a vector or tuple of (dim, tag) pairs"))
    angle=_mesh_attr_finite(corner_angle,caller,"corner_angle")
    0<angle<=pi || throw(ArgumentError(
        "$caller: corner_angle must be in (0, pi] radians (got $corner_angle)"))
    recombine isa Bool || throw(ArgumentError(
        "$caller: recombine must be Bool"))
    selected=Tuple{Int,Int}[]
    if isempty(dim_tags)
        for surface in sort!(collect(keys(m.surfaces)))
            push!(selected,(2,surface))
        end
        for volume in sort!(collect(keys(m.volumes)))
            push!(selected,(3,volume))
        end
    else
        for entry in dim_tags
            pair=if entry isa Pair
                (first(entry),last(entry))
            elseif (entry isa Tuple || entry isa AbstractVector) &&
                   length(entry)==2
                (entry[1],entry[2])
            else
                throw(ArgumentError(
                    "$caller: each dim_tags entry must be a (dim, tag) pair"))
            end
            dimension=_dimension(pair[1],caller)
            t=_tag(pair[2],caller,dimension)
            dimension in (2,3) || throw(ArgumentError(
                "$caller: dim_tags must list surfaces or volumes"))
            _mesh_attr_entity!(m,dimension,t,caller)
            push!(selected,(dimension,t))
        end
    end
    surfaces=Set{Int}()
    for (dimension,t) in selected
        if dimension==2
            push!(surfaces,t)
        else
            shells=get(m.volumes,t,nothing)
            shells===nothing && throw(ArgumentError(
                "$caller: Volume[$t] has no explicit boundary"))
            volume_surfaces=Int[]
            for shell in shells
                signed_surfaces=get(m.surface_loops,shell,nothing)
                signed_surfaces===nothing && throw(ErrorException(
                    "$caller: Volume[$t] references missing Surface " *
                    "Loop[$shell]; rebuild the model"))
                append!(volume_surfaces,abs.(signed_surfaces))
            end
            length(unique(volume_surfaces))==6 || throw(ArgumentError(
                "$caller: Volume[$t] is not 6-sided — transfinite volumes " *
                "require 6 boundary surfaces"))
            set_transfinite_volume!(m,t)
            for surface in volume_surfaces
                push!(surfaces,surface)
            end
        end
    end
    for surface in sort!(collect(surfaces))
        _transfinite_automatic_surface!(m,surface,angle,caller) || continue
        recombine && (m.meshing.recombine[(2,surface)]=45.0)
    end
    return nothing
end

# Ordered boundary point tags of a surface's outer loop — signed curve list
# resolved through `m.curves`, checked for end-to-end connectivity.
function _transfinite_automatic_loop_points(
    m::GeoModel,surface::Int,caller::AbstractString)
    loop_ids=get(m.surfaces,surface,nothing)
    loop_ids===nothing && throw(ArgumentError(
        "$caller: unknown Surface[$surface]"))
    isempty(loop_ids) && return Int[]
    length(loop_ids)==1 || return Int[]
    signed_curves=get(m.loops,loop_ids[1],nothing)
    signed_curves===nothing && throw(ErrorException(
        "$caller: Surface[$surface] references missing Curve " *
        "Loop[$(loop_ids[1])]; rebuild the model"))
    points=Int[]
    for signed_curve in signed_curves
        endpoints=get(m.curves,abs(signed_curve),nothing)
        endpoints===nothing && throw(ErrorException(
            "$caller: Curve Loop[$(loop_ids[1])] references missing " *
            "Curve[$(abs(signed_curve))]; rebuild the model"))
        head,tail=signed_curve>0 ? endpoints : (endpoints[2],endpoints[1])
        if isempty(points)
            push!(points,head,tail)
        else
            points[end]==head || throw(ArgumentError(
                "$caller: Surface[$surface] boundary is not a connected loop"))
            push!(points,tail)
        end
    end
    length(points)>1 && points[end]==points[1] && pop!(points)
    return points
end

# Detect 3 or 4 corner vertices on a surface's boundary loop, then record the
# transfinite surface and matching per-curve node counts.
function _transfinite_automatic_surface!(
    m::GeoModel,surface::Int,corner_angle::Float64,caller::AbstractString)
    points=_transfinite_automatic_loop_points(m,surface,caller)
    n=length(points)
    n>=3 || return false
    coords=NTuple{3,Float64}[]
    for point in points
        coord=get(m.points,point,nothing)
        coord===nothing && throw(ErrorException(
            "$caller: Surface[$surface] boundary references missing " *
            "Point[$point]; rebuild the model"))
        push!(coords,coord)
    end
    # A vertex is a corner when its interior angle does not exceed
    # `corner_angle` — Gmsh ignores faces whose patch corners are flatter than
    # the threshold, so ineligible surfaces are skipped, not rejected.
    corners=Int[]
    for i in 1:n
        previous=coords[mod1(i-1,n)]
        current=coords[i]
        following=coords[mod1(i+1,n)]
        incoming=Float64.(current).-Float64.(previous)
        outgoing=Float64.(following).-Float64.(current)
        ni=sqrt(sum(v->v*v,incoming));no=sqrt(sum(v->v*v,outgoing))
        (ni==0.0 || no==0.0) && continue
        cosine=clamp(sum(incoming.*outgoing)/(ni*no),-1.0,1.0)
        interior=pi-acos(cosine)
        interior<=corner_angle && push!(corners,points[i])
    end
    if n in (3,4)
        length(corners)==n || return false
    else
        length(corners) in (3,4) || return false
    end
    # Distribute one node count per patch edge: equal counts on opposite edges
    # of a quad patch, each edge subdividing its curves proportionally to
    # length. A patch edge runs from one corner to the next along the loop.
    nedges=length(corners)
    edge_lengths=zeros(Float64,nedges)
    edge_curves=[Int[] for _ in 1:nedges]
    # Curve `signed_curves[i]` covers vertices `points[i] -> points[i+1]`
    # because `points` was built in loop order above.
    signed_curves=m.loops[m.surfaces[surface][1]]
    corner_set=Set(corners)
    corner_positions=[i for i in 1:n if points[i] in corner_set]
    for k in 1:nedges
        start=corner_positions[k]
        stop=corner_positions[mod1(k+1,nedges)]
        length_acc=0.0
        i=start
        while i!=stop
            a=coords[i];b=coords[mod1(i+1,n)]
            length_acc+=sqrt(sum((Float64.(b).-Float64.(a)).^2))
            # the curve covering vertex i -> i+1 is the k-th curve index in loop order
            push!(edge_curves[k],signed_curves[i])
            i=mod1(i+1,n)
        end
        edge_lengths[k]=length_acc
    end
    lc=Inf
    for point in points
        lc=min(lc,get(m.point_size,point,Inf))
    end
    isfinite(lc) || (lc=1.0)
    counts=[max(2,round(Int,edge_lengths[k]/lc)+1) for k in 1:nedges]
    if nedges==4
        opposite=max(counts[1],counts[3])
        counts[1]=counts[3]=opposite
        opposite=max(counts[2],counts[4])
        counts[2]=counts[4]=opposite
    end
    for k in 1:nedges
        curves_on_edge=edge_curves[k]
        isempty(curves_on_edge) && continue
        share=max(2,round(Int,counts[k]/length(curves_on_edge)))
        for signed_curve in curves_on_edge
            curve=abs(signed_curve)
            m.meshing.transfinite_curves[curve]=
                (num_nodes=share,kind=:progression,coef=1.0,reversed=false)
        end
    end
    m.meshing.transfinite_surfaces[surface]=
        (arrangement=:left,corners=corners)
    return true
end

# Return the `DiscreteEntity` node/element record for `(dim, tag)` — the
# entity's own record for discrete entities, or a lazily created attachment
# record for native entities (Gmsh stores added mesh data on any entity).
function _discrete_data!(m::GeoModel,dim::Int,tag::Int,caller::AbstractString)
    record=get(m.discrete,(dim,tag),nothing)
    record!==nothing && return record
    _model_entity_known(m,dim,tag) || throw(ArgumentError(
        "$caller: unknown entity ($dim,$tag)"))
    return get!(DiscreteEntity,m.meshing.attached,(dim,tag))
end

# `(key, record)` pairs over both discrete entities and native-entity
# attachments — the model-side counterpart of the API's
# `_discrete_mesh_records`.
function _discrete_mesh_records_model(m::GeoModel)
    pairs=Tuple{Tuple{Int,Int},DiscreteEntity}[]
    for (key,record) in m.discrete
        push!(pairs,(key,record))
    end
    for (key,record) in m.meshing.attached
        push!(pairs,(key,record))
    end
    return pairs
end

"""
    add_discrete_nodes!(model, dim, tag, node_tags, coords, params=Float64[])

Append nodes to entity `(dim, tag)`, matching Gmsh's `mesh.addNodes`.
`node_tags` are unique, strictly positive identifiers; `coords` is flattened
x,y,z triples; `params` optionally carries `dim` parametric coordinates per
node. Duplicate tags update the existing node's coordinates. The entity must
exist — discrete entities use their own record, native entities a meshing
attachment.
"""
function add_discrete_nodes!(m::GeoModel,dim,tag,node_tags,coords,
                             params=Float64[])
    caller="add_discrete_nodes!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    (node_tags isa AbstractVector || node_tags isa Tuple) || throw(ArgumentError(
        "$caller: node_tags must be a vector or tuple"))
    (coords isa AbstractVector || coords isa Tuple) || throw(ArgumentError(
        "$caller: coords must be a flat vector or tuple of x,y,z triples"))
    tags=Int[_mesh_attr_positive_int(t2,caller,"node_tags entry")
             for t2 in node_tags]
    length(coords)==3*length(tags) || throw(ArgumentError(
        "$caller: coords length $(length(coords)) is not 3×" *
        "$(length(tags)) node tags"))
    flat=Float64[_mesh_attr_finite(c,caller,"coords entry") for c in coords]
    pflat=Float64[]
    (params isa AbstractVector || params isa Tuple) || throw(ArgumentError(
        "$caller: params must be a flat vector or tuple"))
    if !isempty(params)
        length(params)==dimension*length(tags) || throw(ArgumentError(
            "$caller: params length $(length(params)) is not $dimension×" *
            "$(length(tags)) node tags"))
        pflat=Float64[_mesh_attr_finite(c,caller,"params entry")
                      for c in params]
    end
    record=_discrete_data!(m,dimension,t,caller)
    for (i,node_tag) in enumerate(tags)
        for (key,other) in _discrete_mesh_records_model(m)
            other===record && continue
            Int32(node_tag) in other.node_tags && throw(ArgumentError(
                "$caller: node tag $node_tag is already owned by entity $key"))
        end
        position=findfirst(==(Int32(node_tag)),record.node_tags)
        if position===nothing
            push!(record.node_tags,Int32(node_tag))
            record.node_coords=hcat(record.node_coords,
                                  reshape(flat[3i-2:3i],3,1))
            if !isempty(pflat)
                column=reshape(pflat[dimension*(i-1)+1:dimension*i],
                               dimension,1)
                if size(record.node_params,1)!=dimension
                    # Earlier nodes carried no parameters — give them zero
                    # columns so the matrix stays aligned with node_tags.
                    record.node_params=zeros(dimension,
                                             length(record.node_tags)-1)
                end
                record.node_params=hcat(record.node_params,column)
            elseif size(record.node_params,2)>0
                record.node_params=hcat(record.node_params,
                                        zeros(dimension,1))
            end
        else
            record.node_coords[:,position].=flat[3i-2:3i]
            if !isempty(pflat)
                if size(record.node_params,1)!=dimension
                    record.node_params=zeros(dimension,
                                             length(record.node_tags))
                end
                record.node_params[:,position].=
                    pflat[dimension*(i-1)+1:dimension*i]
            end
        end
    end
    return nothing
end

"""
    add_discrete_elements!(model, dim, tag, element_types, element_tags,
                           node_tags)

Append elements to entity `(dim, tag)`, matching Gmsh's `mesh.addElements`:
`element_types` lists MSH type numbers and `element_tags[i]`/`node_tags[i]`
the tag vector and flattened connectivity of block `i`. Every referenced node
must already be classified on the entity. Element tags must be unique across
the entity.
"""
function add_discrete_elements!(m::GeoModel,dim,tag,element_types,element_tags,
                                node_tags)
    caller="add_discrete_elements!"
    dimension=_dimension(dim,caller)
    t=_tag(tag,caller,dimension)
    (element_types isa AbstractVector || element_types isa Tuple) ||
        throw(ArgumentError(
            "$caller: element_types must be a vector or tuple of MSH types"))
    (element_tags isa AbstractVector || element_tags isa Tuple) ||
        throw(ArgumentError(
            "$caller: element_tags must be a vector of tag vectors"))
    (node_tags isa AbstractVector || node_tags isa Tuple) || throw(ArgumentError(
        "$caller: node_tags must be a vector of flat connectivity vectors"))
    length(element_tags)==length(element_types) &&
        length(node_tags)==length(element_types) || throw(ArgumentError(
        "$caller: element_types, element_tags, and node_tags must have the " *
        "same block count"))
    record=_discrete_data!(m,dimension,t,caller)
    known=Set(Int.(record.node_tags))
    staged=Tuple{Int32,Int32,Vector{Int32}}[]
    seen=Set{Int32}(record.element_tags)
    for block in eachindex(element_types)
        msh_type=_mesh_attr_positive_int(element_types[block],caller,
                                         "element_types entry")
        spec_nodes=try msh_num_nodes(msh_type) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: unknown MSH element type $msh_type"))
        end
        block_tags=element_tags[block]
        block_nodes=node_tags[block]
        (block_tags isa AbstractVector || block_tags isa Tuple) ||
            throw(ArgumentError(
                "$caller: element_tags[$block] must be a vector of tags"))
        (block_nodes isa AbstractVector || block_nodes isa Tuple) ||
            throw(ArgumentError(
                "$caller: node_tags[$block] must be a flat connectivity " *
                "vector"))
        length(block_nodes)==spec_nodes*length(block_tags) || throw(
            ArgumentError(
                "$caller: node_tags[$block] length $(length(block_nodes)) is " *
                "not $spec_nodes×$(length(block_tags)) element tags"))
        flat_nodes=Int[_mesh_attr_positive_int(v,caller,"node_tags entry")
                       for v in block_nodes]
        for (element_index,raw_tag) in enumerate(block_tags)
            element_tag=_mesh_attr_positive_int(raw_tag,caller,
                                                "element tag")
            Int32(element_tag) in seen && throw(ArgumentError(
                "$caller: duplicate element tag $element_tag"))
            push!(seen,Int32(element_tag))
            connectivity=flat_nodes[
                spec_nodes*(element_index-1)+1:spec_nodes*element_index]
            for node in connectivity
                node in known || throw(ArgumentError(
                    "$caller: element $element_tag references node $node " *
                    "which is not classified on entity ($dimension,$t)"))
            end
            push!(staged,(Int32(msh_type),Int32(element_tag),
                          Int32.(connectivity)))
        end
    end
    for (msh_type,element_tag,connectivity) in staged
        push!(record.element_types,msh_type)
        push!(record.element_tags,element_tag)
        push!(record.element_nodes,connectivity)
    end
    return nothing
end

const _HOMOLOGY_KINDS=("Homology","Cohomology")

"""
    add_homology_request!(model; kind="Homology", domain_tags=Int[],
                          subdomain_tags=Int[], dims=Int[])

Queue a (co)homology computation request, matching Gmsh's
`mesh.addHomologyRequest`. `kind` is `"Homology"` or `"Cohomology"`;
`domain_tags`/`subdomain_tags` list entity tags and `dims` the chain
dimensions to compute.
"""
function add_homology_request!(m::GeoModel;kind="Homology",domain_tags=Int[],
                               subdomain_tags=Int[],dims=Int[])
    caller="add_homology_request!"
    kind isa AbstractString || throw(ArgumentError(
        "$caller: kind must be a string"))
    String(kind) in _HOMOLOGY_KINDS || throw(ArgumentError(
        "$caller: kind must be \"Homology\" or \"Cohomology\" (got \"$kind\")"))
    domain=Int[_mesh_attr_positive_int(t,caller,"domain_tags entry")
               for t in domain_tags]
    subdomain=Int[_mesh_attr_positive_int(t,caller,"subdomain_tags entry")
                  for t in subdomain_tags]
    for t in Iterators.flatten((domain,subdomain))
        0<=t || throw(ArgumentError("$caller: entity tags must be non-negative"))
    end
    dimensions=Int[]
    for raw in dims
        value=_mesh_attr_positive_int(raw,caller,"dims entry")
        0<=value<=3 || throw(ArgumentError(
            "$caller: dims entry $value must be in 0:3"))
        push!(dimensions,value)
    end
    push!(m.meshing.homology_requests,
          (kind=String(kind),domain=domain,subdomain=subdomain,
           dims=dimensions))
    return nothing
end

"""Clear every queued homology request, matching `mesh.clearHomologyRequests`."""
function clear_homology_requests!(m::GeoModel)
    empty!(m.meshing.homology_requests)
    return nothing
end


# ────────────────────────────────────────────────────────────────────────────────
# Topology creation from discrete mesh data
# ────────────────────────────────────────────────────────────────────────────────

# Boundary-facet tables over the *corner* nodes of each fixed-node MSH type
# (1-based local positions). High-order types share their family's corners.
const _RECORD_FACET_FAMILIES = Dict{Symbol,Vector{NTuple{4,Int}}}(
    :pnt => [(1,-1,-1,-1)],
    :lin => [(1,-1,-1,-1),(2,-1,-1,-1)],
    :tri => [(1,2,-1,-1),(2,3,-1,-1),(3,1,-1,-1)],
    :qua => [(1,2,-1,-1),(2,3,-1,-1),(3,4,-1,-1),(4,1,-1,-1)],
    :tet => [(1,3,2,-1),(1,2,4,-1),(1,4,3,-1),(2,3,4,-1)],
    :hex => [(1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)],
    :pri => [(1,3,2,-1),(4,5,6,-1),(1,2,5,4),(2,3,6,5),(3,1,4,6)],
    :pyr => [(1,4,3,2),(1,2,5,-1),(2,3,5,-1),(3,4,5,-1),(4,1,5,-1)],
)

# Facets of one record cell: a vector of node-tag vectors, one per facet.
function _record_cell_facets(msh_type::Int,nodes::Vector{Int32})
    family=msh_family(msh_type)
    table=get(_RECORD_FACET_FAMILIES,family,nothing)
    table===nothing && return Vector{Int32}[]
    facets=Vector{Int32}[]
    for facet in table
        positions=Int[position for position in facet if position>0]
        all(position->position<=length(nodes),positions) || return Vector{Int32}[]
        push!(facets,Int32[nodes[position] for position in positions])
    end
    return facets
end

# Facet key invariant to node ordering within the facet.
_facet_key(facet::Vector{Int32})=Tuple(sort(facet))

# Edge key for facet-adjacency grouping (shared `dim-2` features).
function _facet_subfacets(msh_type::Int,facet::Vector{Int32})
    # Facet adjacency uses shared corner pairs for surface cells and shared
    # endpoints for segments; higher dims use shared edges of the facet.
    if length(facet)<=2
        return [Tuple(sort(facet[i:i])) for i in eachindex(facet)]
    end
    edges=NTuple{2,Int32}[]
    n=length(facet)
    for i in 1:n
        a=facet[i];b=facet[mod1(i+1,n)]
        push!(edges,a<b ? (a,b) : (b,a))
    end
    return edges
end

# `(cells, cell_types)` pairs for every record element living on `dim`.
function _record_dim_cells(record::DiscreteEntity,dim::Int)
    cells=Vector{Int32}[]
    for (i,msh_type) in enumerate(record.element_types)
        msh_dimension(msh_type)==dim || continue
        push!(cells,record.element_nodes[i])
    end
    return cells
end

# Nodes exclusively used by `owned` cells are removed from the record's node
# list and returned as a tag → coordinate/params table for migration.
function _record_extract_nodes!(record::DiscreteEntity,keep::Set{Int32})
    migrated=Dict{Int32,Tuple{Vector{Float64},Vector{Float64}}}()
    mask=trues(length(record.node_tags))
    for (i,tag) in enumerate(record.node_tags)
        tag in keep && continue
        mask[i]=false
        migrated[tag]=(Vector{Float64}(record.node_coords[:,i]),
                       size(record.node_params,1)>0 ?
                           Vector{Float64}(record.node_params[:,i]) :
                           Float64[])
    end
    record.node_tags=record.node_tags[mask]
    record.node_coords=record.node_coords[:,mask]
    size(record.node_params,1)>0 &&
        (record.node_params=record.node_params[:,mask])
    return migrated
end

function _record_append_node!(record::DiscreteEntity,tag::Int32,
                              coords::Vector{Float64},
                              params::Vector{Float64})
    push!(record.node_tags,tag)
    record.node_coords=hcat(record.node_coords,reshape(coords,3,1))
    if size(record.node_params,1)>0 || !isempty(params)
        width=max(size(record.node_params,1),length(params))
        if size(record.node_params,1)!=width
            grown=zeros(width,length(record.node_tags)-1)
            grown[1:size(record.node_params,1),:].=record.node_params
            record.node_params=grown
        end
        column=zeros(width)
        column[1:length(params)].=params
        record.node_params=hcat(record.node_params,column)
    end
    return nothing
end

function _record_append_element!(record::DiscreteEntity,msh_type::Int,
                                 element_tag::Int32,nodes::Vector{Int32})
    push!(record.element_types,Int32(msh_type))
    push!(record.element_tags,element_tag)
    push!(record.element_nodes,copy(nodes))
    return nothing
end

# Fresh element tag above every record and dense-cache tag in the model.
function _model_fresh_element_tag(m::GeoModel)
    result=0
    for (_,record) in _discrete_mesh_records_model(m)
        isempty(record.element_tags) ||
            (result=max(result,Int(maximum(record.element_tags))))
    end
    return result+1
end

# Chain boundary segments (node-tag pairs) into loops or open chains. Each
# entry is `(node sequence, closed)`.
function _chain_boundary_edges(edges::Vector{NTuple{2,Int32}})
    incidence=Dict{Int32,Vector{Int}}()
    for (i,(a,b)) in enumerate(edges)
        push!(get!(Vector{Int},incidence,a),i)
        push!(get!(Vector{Int},incidence,b),i)
    end
    used=falses(length(edges))
    chains=Tuple{Vector{Int32},Bool}[]
    for start in eachindex(edges)
        used[start] && continue
        a,b=edges[start]
        chain=Int32[a,b]
        used[start]=true
        while true
            last_node=chain[end]
            advanced=false
            for edge_index in get(incidence,last_node,Int[])
                used[edge_index] && continue
                u,v=edges[edge_index]
                used[edge_index]=true
                push!(chain,u==last_node ? v : u)
                advanced=true
                break
            end
            advanced || break
        end
        # The walk may have started mid-chain: extend the other direction too.
        while true
            first_node=chain[1]
            advanced=false
            for edge_index in get(incidence,first_node,Int[])
                used[edge_index] && continue
                u,v=edges[edge_index]
                used[edge_index]=true
                pushfirst!(chain,u==first_node ? v : u)
                advanced=true
                break
            end
            advanced || break
        end
        push!(chains,(chain,chain[1]==chain[end]))
    end
    return chains
end

# Group boundary facets into components joined by shared facet edges.
function _group_boundary_facets(facets::Vector{Vector{Int32}})
    n=length(facets)
    n==0 && return Vector{Int}[]
    edge_owner=Dict{NTuple{2,Int32},Vector{Int}}()
    for (i,facet) in enumerate(facets)
        count=length(facet)
        for k in 1:count
            a=facet[k];b=facet[mod1(k+1,count)]
            push!(get!(Vector{Int},edge_owner,a<b ? (a,b) : (b,a)),i)
        end
    end
    components=Vector{Int}[]
    seen=falses(n)
    for start in 1:n
        seen[start] && continue
        component=Int[]
        stack=[start];seen[start]=true
        while !isempty(stack)
            f=pop!(stack);push!(component,f)
            facet=facets[f]
            count=length(facet)
            for k in 1:count
                a=facet[k];b=facet[mod1(k+1,count)]
                for neighbor in get(edge_owner,a<b ? (a,b) : (b,a),Int[])
                    seen[neighbor] && continue
                    seen[neighbor]=true;push!(stack,neighbor)
                end
            end
        end
        push!(components,component)
    end
    return components
end

# Move `tags` out of `record`'s node table into `child`, preserving
# coordinates and parameters. Element connectivity still references the tags —
# record queries resolve node tags through the global record index.
function _record_migrate_nodes!(record::DiscreteEntity,child::DiscreteEntity,
                                tags::Set{Int32})
    keep=setdiff(Set{Int32}(record.node_tags),tags)
    migrated=_record_extract_nodes!(record,keep)
    for (tag,(coords,params)) in sort!(collect(migrated);by=first)
        tag in tags || continue
        _record_append_node!(child,tag,coords,params)
    end
    return nothing
end

# Derive a boundary entity of dimension `dim-1` for `record` at `key`. Returns
# the child entity tag or `nothing` when no boundary exists.
function _record_derive_boundary!(m::GeoModel,dim::Int,tag::Int,
                                  record::DiscreteEntity,caller::AbstractString)
    cells=_record_dim_cells(record,dim)
    isempty(cells) && return nothing
    element_index=0
    incidence=Dict{Any,Vector{Tuple{Int,Int}}}()
    typed_cells=Tuple{Int,Vector{Int32}}[]
    for (i,msh_type) in enumerate(record.element_types)
        msh_dimension(msh_type)==dim || continue
        element_index+=1
        push!(typed_cells,(msh_type,record.element_nodes[i]))
        for (f,facet) in enumerate(
                _record_cell_facets(Int(msh_type),record.element_nodes[i]))
            push!(get!(Vector{Tuple{Int,Int}},incidence,
                       _facet_key(facet)),(element_index,f))
        end
    end
    boundary_facets=Tuple{Int,Int,Vector{Int32}}[]
    for (key2,hits) in incidence
        length(hits)==1 || continue
        cell_index,facet_index=hits[1]
        msh_type,nodes=typed_cells[cell_index]
        facet=_record_cell_facets(msh_type,nodes)[facet_index]
        push!(boundary_facets,(msh_type,facet_index,facet))
    end
    if dim==1 && isempty(boundary_facets) && !isempty(cells)
        # A closed curve loop still gets a seam point (Gmsh semantics).
        node=cells[1][1]
        point_tag=add_discrete_entity!(m,0,0)
        point_record=m.discrete[(0,point_tag)]
        _record_migrate_nodes!(record,point_record,Set{Int32}([node]))
        _record_append_element!(point_record,15,
            Int32(_model_fresh_element_tag(m)),Int32[node])
        record.boundary=NTuple{2,Int}[(0,point_tag)]
        return nothing
    end
    isempty(boundary_facets) && return nothing
    if dim==1
        # Boundary facets of segments are endpoint nodes: each maximal
        # connected boundary-node cluster becomes a point entity.
        seen_nodes=Set{Int32}()
        for (_,_,facet) in boundary_facets
            push!(seen_nodes,facet[1])
        end
        point_tags=Int[]
        for node in sort!(collect(seen_nodes))
            point_tag=add_discrete_entity!(m,0,0)
            point_record=m.discrete[(0,point_tag)]
            _record_migrate_nodes!(record,point_record,Set{Int32}([node]))
            _record_append_element!(point_record,15,
                Int32(_model_fresh_element_tag(m)),Int32[node])
            push!(point_tags,point_tag)
        end
        record.boundary=NTuple{2,Int}[(0,t) for t in point_tags]
        return nothing
    end
    if dim==2
        # Boundary edges chain into loops; each loop becomes one curve plus a
        # seam point (closed) or two endpoint points (open).
        edges=NTuple{2,Int32}[]
        for (_,_,facet) in boundary_facets
            length(facet)==2 && push!(edges,(facet[1],facet[2]))
        end
        isempty(edges) && return nothing
        curve_tags=Int[]
        for (chain,closed) in _chain_boundary_edges(edges)
            curve_tag=add_discrete_entity!(m,1,0)
            curve_record=m.discrete[(1,curve_tag)]
            chain_nodes=Set{Int32}(chain[1:end- (closed ? 1 : 0)])
            _record_migrate_nodes!(record,curve_record,chain_nodes)
            for k in 1:length(chain)-1
                _record_append_element!(curve_record,1,
                    Int32(_model_fresh_element_tag(m)),
                    Int32[chain[k],chain[k+1]])
            end
            point_tags=Int[]
            endpoints=closed ? Int32[chain[1]] : Int32[chain[1],chain[end]]
            for node in endpoints
                point_tag=add_discrete_entity!(m,0,0)
                point_record=m.discrete[(0,point_tag)]
                _record_migrate_nodes!(curve_record,point_record,
                                       Set{Int32}([node]))
                _record_append_element!(point_record,15,
                    Int32(_model_fresh_element_tag(m)),Int32[node])
                push!(point_tags,point_tag)
            end
            curve_record.boundary=NTuple{2,Int}[(0,t) for t in point_tags]
            push!(curve_tags,curve_tag)
        end
        record.boundary=NTuple{2,Int}[(1,t) for t in curve_tags]
        return nothing
    end
    # dim==3: boundary faces group into components, each becoming a surface.
    facets=Vector{Int32}[facet for (_,_,facet) in boundary_facets]
    facet_types=Int[msh_type for (msh_type,_,_) in boundary_facets]
    components=_group_boundary_facets(facets)
    surface_tags=Int[]
    for component in components
        surface_tag=add_discrete_entity!(m,2,0)
        surface_record=m.discrete[(2,surface_tag)]
        nodes=Set{Int32}()
        for index in component
            union!(nodes,facets[index])
        end
        _record_migrate_nodes!(record,surface_record,nodes)
        for index in component
            facet=facets[index]
            child_type=length(facet)==3 ? 2 : 3
            _record_append_element!(surface_record,child_type,
                Int32(_model_fresh_element_tag(m)),facet)
        end
        push!(surface_tags,surface_tag)
    end
    record.boundary=NTuple{2,Int}[(2,t) for t in surface_tags]
    return nothing
end

"""
    create_topology!(model; make_simply_connected=true, export_discrete=true)

Create a boundary representation from discrete mesh data, matching Gmsh's
`mesh.createTopology`: discrete entities that carry elements but no declared
boundary derive boundary entities one dimension down (volumes → surfaces,
surfaces → boundary-loop curves plus seam points, curves → endpoint points).
New entities own their boundary elements and nodes; parent records keep
higher-dimensional connectivity referencing the migrated node tags.
`make_simply_connected` is accepted for Gmsh compatibility — splits of
non-simply-connected entities are not performed. `export_discrete` is
accepted; native CAD entities are retained.
"""
function create_topology!(m::GeoModel;make_simply_connected=true,
                          export_discrete=true)
    caller="create_topology!"
    make_simply_connected isa Bool || throw(ArgumentError(
        "$caller: make_simply_connected must be Bool"))
    export_discrete isa Bool || throw(ArgumentError(
        "$caller: export_discrete must be Bool"))
    for dim in (3,2,1)
        for key in sort!(collect(keys(m.discrete)))
            key[1]==dim || continue
            record=m.discrete[key]
            isempty(record.boundary) || continue
            _record_derive_boundary!(m,dim,key[2],record,caller)
        end
    end
    return nothing
end

# Deduplicated edge set of a patch-boundary graph chained into open chains
# (starting and ending at multi-valent nodes, degree != 2) and closed loops.
# Returns `(chains, multivalent)` where each chain is `(node sequence, closed)`.
function _chain_boundary_graph(edges::Vector{NTuple{2,Int32}})
    incidence=Dict{Int32,Vector{Int}}()
    for (i,(a,b)) in enumerate(edges)
        push!(get!(Vector{Int},incidence,a),i)
        push!(get!(Vector{Int},incidence,b),i)
    end
    used=falses(length(edges))
    multivalent=Set{Int32}()
    for (node,indices) in incidence
        length(indices)==2 || push!(multivalent,node)
    end
    chains=Tuple{Vector{Int32},Bool}[]
    for start_node in sort!(collect(multivalent))
        for edge_index in incidence[start_node]
            used[edge_index] && continue
            used[edge_index]=true
            u,v=edges[edge_index]
            chain=Int32[start_node,u==start_node ? v : u]
            while true
                last_node=chain[end]
                last_node in multivalent && break
                advanced=false
                for next_index in incidence[last_node]
                    used[next_index] && continue
                    used[next_index]=true
                    x,y=edges[next_index]
                    push!(chain,x==last_node ? y : x)
                    advanced=true
                    break
                end
                advanced || push!(multivalent,last_node)
                advanced || break
            end
            push!(chains,(chain,false))
        end
    end
    for start in eachindex(edges)
        used[start] && continue
        used[start]=true
        a,b=edges[start]
        chain=Int32[a,b]
        while chain[end]!=chain[1]
            last_node=chain[end]
            advanced=false
            for next_index in incidence[last_node]
                used[next_index] && continue
                used[next_index]=true
                x,y=edges[next_index]
                push!(chain,x==last_node ? y : x)
                advanced=true
                break
            end
            advanced || break
        end
        push!(chains,(chain,chain[1]==chain[end]))
    end
    return chains,multivalent
end

# Unit normal of a 2D element from its corner coordinates (Newell for quads).
function _classify_normal(coords::Dict{Int32,NTuple{3,Float64}},
                          corners::Vector{Int32})
    normal=zeros(3)
    count=length(corners)
    for k in 1:count
        p=coords[corners[k]];q=coords[corners[mod1(k+1,count)]]
        normal[1]+=(p[2]-q[2])*(p[3]+q[3])
        normal[2]+=(p[3]-q[3])*(p[1]+q[1])
        normal[3]+=(p[1]-q[1])*(p[2]+q[2])
    end
    norm2=sqrt(normal[1]^2+normal[2]^2+normal[3]^2)
    norm2>0 || return normal
    return normal./norm2
end

# Corner node tags (first-order vertices) of a 2D element record entry.
_classify_corners(msh_type::Int,nodes::Vector{Int32})=
    nodes[1:min(length(nodes),msh_family(msh_type)===:qua ? 4 : 3)]

# Edges of a 2D element restricted to its corner polygon.
function _classify_edges(msh_type::Int,nodes::Vector{Int32})
    corners=_classify_corners(msh_type,nodes)
    count=length(corners)
    edges=NTuple{2,Int32}[]
    for k in 1:count
        a=corners[k];b=corners[mod1(k+1,count)]
        push!(edges,a<b ? (a,b) : (b,a))
    end
    return edges
end

# Split each chain at interior nodes where the angle between consecutive
# boundary edges exceeds `curve_angle`; split nodes become open endpoints.
function _split_chains_by_angle(chains::Vector{Tuple{Vector{Int32},Bool}},
                                coords::Dict{Int32,NTuple{3,Float64}},
                                curve_angle::Float64)
    curve_angle>=pi && return chains
    result=Tuple{Vector{Int32},Bool}[]
    for (chain,closed) in chains
        cuts=Int[]
        limit=closed ? length(chain)-1 : length(chain)-1
        for k in 2:limit
            a=coords[chain[k-1]];b=coords[chain[k]];c=coords[chain[k+1]]
            u=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
            v=(c[1]-b[1],c[2]-b[2],c[3]-b[3])
            nu=sqrt(u[1]^2+u[2]^2+u[3]^2);nv=sqrt(v[1]^2+v[2]^2+v[3]^2)
            (nu>0 && nv>0) || continue
            cosine=clamp((u[1]*v[1]+u[2]*v[2]+u[3]*v[3])/(nu*nv),-1.0,1.0)
            acos(cosine)>curve_angle && push!(cuts,k)
        end
        if isempty(cuts)
            push!(result,(chain,closed))
            continue
        end
        start=1
        for cut in cuts
            push!(result,(chain[start:cut],false))
            start=cut
        end
        closed && start==1 || push!(result,(chain[start:end],false))
        if closed && start>1
            # Wrap segment between last cut and first cut through the seam.
            wrapped=vcat(chain[start:end-1],chain[1:cuts[1]])
            pop!(result)
            push!(result,(wrapped,false))
        end
    end
    return result
end

"""
    classify_surfaces!(model; angle, boundary, for_reparametrization,
                       curve_angle)

Classify the discrete surface mesh into patches whose adjacent elements meet
at a dihedral angle below `angle` (radians), matching Gmsh's
`mesh.classifySurfaces`. Every patch becomes a new discrete surface that keeps
its original element tags; each patch-boundary chain (interface or outer
boundary) becomes a discrete curve shared by the patches it borders; chain
endpoints and multi-valent nodes become discrete points owning their node.
Interior chain nodes migrate to their curve record and interior surface nodes
stay on their patch record. Chains split further where consecutive boundary
edges meet above `curve_angle`. `boundary=false` skips curve/point creation;
`for_reparametrization` is accepted for Gmsh compatibility. Source entities
whose elements were fully redistributed are removed.
"""
function classify_surfaces!(m::GeoModel;angle=40*pi/180,boundary=true,
                            for_reparametrization=false,curve_angle=pi)
    caller="classify_surfaces!"
    for (value,name) in ((angle,"angle"),(curve_angle,"curve_angle"))
        (value isa Real && !(value isa Bool) && isfinite(value)) ||
            throw(ArgumentError("$caller: $name must be a finite real number"))
    end
    boundary isa Bool || throw(ArgumentError("$caller: boundary must be Bool"))
    for_reparametrization isa Bool || throw(ArgumentError(
        "$caller: for_reparametrization must be Bool"))
    threshold=cos(Float64(angle))
    cut=Float64(curve_angle)

    coords=Dict{Int32,NTuple{3,Float64}}()
    for (_,record) in _discrete_mesh_records_model(m)
        for column in eachindex(record.node_tags)
            coords[record.node_tags[column]]=
                (record.node_coords[1,column],record.node_coords[2,column],
                 record.node_coords[3,column])
        end
    end

    elements=Tuple{Tuple{Int,Int},Int,Int,Vector{Int32}}[]
    records=Dict{Tuple{Int,Int},DiscreteEntity}()
    for (key,record) in _discrete_mesh_records_model(m)
        records[key]=record
        key[1]==2 || continue
        for index in eachindex(record.element_types)
            msh_dimension(record.element_types[index])==2 || continue
            nodes=record.element_nodes[index]
            corners=_classify_corners(Int(record.element_types[index]),nodes)
            all(haskey(coords,corner) for corner in corners) ||
                throw(ArgumentError(
                    "$caller: element $(record.element_tags[index]) " *
                    "references nodes without coordinates"))
            push!(elements,(key,index,Int(record.element_types[index]),nodes))
        end
    end
    isempty(elements) && return nothing

    normals=Vector{Vector{Float64}}(undef,length(elements))
    edge_owner=Dict{NTuple{2,Int32},Vector{Int}}()
    for (i,(_,_,msh_type,nodes)) in enumerate(elements)
        corners=_classify_corners(msh_type,nodes)
        normals[i]=_classify_normal(coords,corners)
        for edge in _classify_edges(msh_type,nodes)
            push!(get!(Vector{Int},edge_owner,edge),i)
        end
    end

    parent=collect(eachindex(elements))
    function _find(x)
        while parent[x]!=x
            parent[x]=parent[parent[x]]
            x=parent[x]
        end
        return x
    end
    for owners in values(edge_owner)
        for i in 2:length(owners)
            a=_find(owners[1]);b=_find(owners[i])
            a==b && continue
            dot=clamp(normals[a][1]*normals[b][1]+normals[a][2]*normals[b][2]+
                      normals[a][3]*normals[b][3],-1.0,1.0)
            dot>threshold && (parent[b]=a)
        end
    end

    patches=Dict{Int,Vector{Int}}()
    for i in eachindex(elements)
        push!(get!(Vector{Int},patches,_find(i)),i)
    end
    length(patches)<=1 && !boundary && return nothing
    patch_of=[_find(i) for i in eachindex(elements)]
    patch_ids=sort!(collect(keys(patches)))

    # Boundary edges of each patch: edges whose owner set differs from the
    # patch alone (outer boundary or patch interface).
    edge_patches=Dict{NTuple{2,Int32},Set{Int}}()
    for (edge,owners) in edge_owner
        adjacent=Set{Int}(patch_of[i] for i in owners)
        for i in owners
            if length(owners)==1 || length(adjacent)>1
                push!(get!(Set{Int},edge_patches,edge),patch_of[i])
            end
        end
    end

    curves=NTuple{2,Int}[]
    curve_nodes=Dict{Int,Vector{Int32}}()
    curve_patches=Dict{Int,Set{Int}}()
    point_of=Dict{Int32,Int}()
    if boundary
        graph_edges=collect(keys(edge_patches))
        chains,multivalent=_chain_boundary_graph(graph_edges)
        chains=_split_chains_by_angle(chains,coords,cut)
        for (chain,closed) in chains
            curve_tag=add_discrete_entity!(m,1,0)
            curve_record=m.discrete[(1,curve_tag)]
            involved=Set{Int}()
            for k in 1:length(chain)-1
                a,b=chain[k],chain[k+1]
                edge=a<b ? (a,b) : (b,a)
                union!(involved,get(edge_patches,edge,Set{Int}()))
                _record_append_element!(curve_record,1,
                    Int32(_model_fresh_element_tag(m)),Int32[a,b])
            end
            interior=closed ? chain[1:end-1] : chain[2:end-1]
            curve_nodes[curve_tag]=Int32[interior...]
            curve_patches[curve_tag]=involved
            endpoints=closed ? Int32[chain[1]] : Int32[chain[1],chain[end]]
            point_tags=Int[]
            for node in endpoints
                tag=get!(point_of,node) do
                    point_tag=add_discrete_entity!(m,0,0)
                    point_record=m.discrete[(0,point_tag)]
                    _record_append_element!(point_record,15,
                        Int32(_model_fresh_element_tag(m)),Int32[node])
                    point_tag
                end
                push!(point_tags,tag)
            end
            curve_record.boundary=NTuple{2,Int}[(0,t) for t in point_tags]
            push!(curves,(1,curve_tag))
        end
    end

    # Node ownership: points > curves > patches. Interior nodes stay with
    # their patch; boundary chains own their interior nodes.
    point_nodes=Set{Int32}(keys(point_of))
    curve_owner=Dict{Int32,Int}()
    for (tag,nodes) in curve_nodes
        for node in nodes
            node in point_nodes || (curve_owner[node]=tag)
        end
    end

    new_surfaces=Dict{Int,DiscreteEntity}()
    source_keys=Set{Tuple{Int,Int}}()
    for patch in patch_ids
        surface_tag=add_discrete_entity!(m,2,0)
        new_surfaces[patch]=m.discrete[(2,surface_tag)]
    end
    for (i,(key,index,msh_type,nodes)) in enumerate(elements)
        patch=patch_of[i]
        source=records[key]
        target=new_surfaces[patch]
        push!(source_keys,key)
        _record_append_element!(target,msh_type,
            source.element_tags[index],nodes)
        source.element_types[index]=Int32(0)
    end
    for key in source_keys
        record=records[key]
        keep=findall(t->t!=0,record.element_types)
        record.element_types=record.element_types[keep]
        record.element_tags=record.element_tags[keep]
        record.element_nodes=record.element_nodes[keep]
    end

    surface_boundary=Dict{Int,Vector{Int}}(patch=>Int[] for patch in patch_ids)
    for (curve_tag,involved) in curve_patches
        for patch in involved
            push!(surface_boundary[patch],curve_tag)
        end
    end
    for patch in patch_ids
        record=new_surfaces[patch]
        record.boundary=NTuple{2,Int}[(1,t) for t in
                                      sort!(surface_boundary[patch])]
        interior=Set{Int32}()
        for i in patches[patch]
            union!(interior,_classify_corners(elements[i][3],elements[i][4]))
            union!(interior,elements[i][4])
        end
        for node in interior
            node in point_nodes && continue
            haskey(curve_owner,node) && continue
            _record_append_node!(record,node,
                [coords[node][1],coords[node][2],coords[node][3]],
                Float64[])
        end
    end
    for (curve_tag,nodes) in curve_nodes
        record=m.discrete[(1,curve_tag)]
        for node in nodes
            haskey(curve_owner,node) || continue
            _record_append_node!(record,node,
                [coords[node][1],coords[node][2],coords[node][3]],
                Float64[])
        end
    end
    for (node,tag) in point_of
        record=m.discrete[(0,tag)]
        _record_append_node!(record,node,
            [coords[node][1],coords[node][2],coords[node][3]],Float64[])
    end

    # Sources emptied by the redistribution are dropped; attached records on
    # native entities are cleared, discrete entities are removed.
    stale=Tuple{Int,Int}[]
    for key in source_keys
        record=records[key]
        if haskey(m.meshing.attached,key)
            isempty(record.element_types) && delete!(m.meshing.attached,key)
        elseif isempty(record.element_types)
            push!(stale,key)
        else
            record.boundary=NTuple{2,Int}[]
        end
    end
    for key in stale
        boundary_entities=Tuple{Int,Int}[]
        haskey(m.discrete,key) || continue
        for entity in m.discrete[key].boundary
            push!(boundary_entities,entity)
        end
        delete!(m.discrete,key)
        for entity in boundary_entities
            _drop_unreferenced_discrete!(m,entity)
        end
    end
    return nothing
end

# Delete a discrete entity when no remaining entity lists it as a boundary;
# recurse into its own boundary so orphaned chains clean up fully.
function _drop_unreferenced_discrete!(m::GeoModel,key::Tuple{Int,Int})
    haskey(m.discrete,key) || return nothing
    for (_,record) in m.discrete
        key in record.boundary && return nothing
    end
    children=Tuple{Int,Int}[entity for entity in m.discrete[key].boundary]
    delete!(m.discrete,key)
    for child in children
        _drop_unreferenced_discrete!(m,child)
    end
    return nothing
end

# ---- create_geometry ------------------------------------------------------
#
# Gmsh's `mesh.createGeometry` computes a parametrization for discrete
# entities so value/parametrization queries and remeshing work. Tessella
# parametrizes discrete curves by normalized cumulative chord length and
# discrete surfaces by projection onto their least-squares (PCA) plane; the
# resulting per-node parameters are stored on the entity record — on
# `node_params` for nodes the record owns and on `aux_params` for nodes owned
# by boundary entities, mirroring Gmsh's per-entity node parameter storage.

# Node-tag → coordinate lookup across every record of the model.
function _model_record_node_coords(m::GeoModel)
    coords=Dict{Int32,NTuple{3,Float64}}()
    for (_,record) in _discrete_mesh_records_model(m)
        for column in eachindex(record.node_tags)
            coords[record.node_tags[column]]=
                (record.node_coords[1,column],record.node_coords[2,column],
                 record.node_coords[3,column])
        end
    end
    return coords
end

# Every node tag referenced by `record`'s elements of dimension `dim`.
function _record_referenced_nodes(record::DiscreteEntity,dim::Int)
    tags=Set{Int32}()
    for (index,msh_type) in enumerate(record.element_types)
        msh_dimension(msh_type)==dim || continue
        union!(tags,record.element_nodes[index])
    end
    return sort!(collect(tags))
end

# Store `params` for `node` on `target` — into `node_params` when the record
# owns the node, into `aux_params` when a boundary entity owns it.
function _record_store_param!(m::GeoModel,target::DiscreteEntity,dim::Int,
                              node::Int32,params::Vector{Float64},
                              caller::AbstractString)
    column=findfirst(==(node),target.node_tags)
    if column===nothing
        target.aux_params[node]=params
        return nothing
    end
    if size(target.node_params,1)!=dim
        grown=zeros(dim,length(target.node_tags))
        rows=size(target.node_params,1)
        rows>0 && (grown[1:min(rows,dim),:].=target.node_params[1:min(rows,dim),:])
        target.node_params=grown
    end
    target.node_params[:,column].=params
    return nothing
end

# Cumulative-chord parametrization of a discrete curve record's segments.
function _parametrize_discrete_curve!(m::GeoModel,record::DiscreteEntity,
                                      coords::Dict{Int32,NTuple{3,Float64}},
                                      caller::AbstractString)
    edges=NTuple{2,Int32}[]
    for (index,msh_type) in enumerate(record.element_types)
        msh_family(msh_type)===:lin || continue
        nodes=record.element_nodes[index]
        push!(edges,(nodes[1],nodes[2]))
    end
    isempty(edges) && return nothing
    chains,_=_chain_boundary_graph(edges)
    # Order chains deterministically by their smallest endpoint tag.
    sort!(chains,by=entry->minimum(entry[1]))
    lengths=Float64[]
    for (chain,_) in chains
        total=0.0
        for k in 1:length(chain)-1
            a=coords[chain[k]];b=coords[chain[k+1]]
            total+=sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2)
        end
        push!(lengths,total)
    end
    overall=sum(lengths)
    overall>0 || throw(ArgumentError(
        "$caller: discrete curve has zero total length"))
    offset=0.0
    for (chain_index,(chain,_)) in enumerate(chains)
        scale=lengths[chain_index]/overall
        position=0.0
        local_total=lengths[chain_index]
        for k in eachindex(chain)
            if k>1
                a=coords[chain[k-1]];b=coords[chain[k]]
                position+=sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2)
            end
            parameter=offset+(local_total>0 ? scale*position/local_total : 0.0)
            _record_store_param!(m,record,1,chain[k],[parameter],caller)
        end
        offset+=scale
    end
    return nothing
end

# Least-squares-plane parametrization of a discrete surface record.
function _parametrize_discrete_surface!(m::GeoModel,record::DiscreteEntity,
                                        coords::Dict{Int32,NTuple{3,Float64}},
                                        caller::AbstractString)
    tags=_record_referenced_nodes(record,2)
    isempty(tags) && return nothing
    for node in tags
        haskey(coords,node) || throw(ArgumentError(
            "$caller: node $node has no coordinates"))
    end
    centroid=[0.0,0.0,0.0]
    for node in tags
        point=coords[node]
        centroid[1]+=point[1];centroid[2]+=point[2];centroid[3]+=point[3]
    end
    centroid./=length(tags)
    covariance=zeros(3,3)
    for node in tags
        point=coords[node]
        delta=[point[1]-centroid[1],point[2]-centroid[2],point[3]-centroid[3]]
        for i in 1:3,j in 1:3
            covariance[i,j]+=delta[i]*delta[j]
        end
    end
    factorization=eigen(Symmetric(covariance))
    order=sortperm(factorization.values;rev=true)
    first_axis=factorization.vectors[:,order[1]]
    second_axis=factorization.vectors[:,order[2]]
    for node in tags
        point=coords[node]
        delta=(point[1]-centroid[1],point[2]-centroid[2],point[3]-centroid[3])
        u=delta[1]*first_axis[1]+delta[2]*first_axis[2]+delta[3]*first_axis[3]
        v=delta[1]*second_axis[1]+delta[2]*second_axis[2]+
          delta[3]*second_axis[3]
        _record_store_param!(m,record,2,node,[u,v],caller)
    end
    return nothing
end

"""
    create_geometry!(model, dim_tags=nothing)

Compute a parametrization for discrete entities (all discrete entities when
`dim_tags` is `nothing`), matching Gmsh's `mesh.createGeometry`. Discrete
curves receive normalized cumulative-chord parameters; discrete surfaces
receive least-squares-plane parameters covering every node their elements
reference — including nodes owned by boundary entities, stored in the record's
`aux_params` side table.
"""
function create_geometry!(m::GeoModel,dim_tags=nothing)
    caller="create_geometry!"
    targets=Tuple{Int,Int}[]
    if dim_tags===nothing
        for key in sort!(collect(keys(m.discrete)))
            key[1] in (1,2) && push!(targets,key)
        end
    else
        (dim_tags isa AbstractVector || dim_tags isa Tuple) ||
            throw(ArgumentError("$caller: dim_tags must be a vector or tuple"))
        for entry in dim_tags
            (entry isa Tuple || entry isa Pair ||
             (entry isa AbstractVector && length(entry)==2)) &&
                length(entry)==2 || throw(ArgumentError(
                    "$caller: dim_tags entries must be (dim, tag) pairs"))
            dimension=Int(first(entry));tag=Int(last(entry))
            dimension in (1,2) || throw(ArgumentError(
                "$caller: only curve and surface entities can be parametrized"))
            haskey(m.discrete,(dimension,tag)) || throw(ArgumentError(
                "$caller: unknown discrete entity ($dimension,$tag)"))
            push!(targets,(dimension,tag))
        end
    end
    isempty(targets) && return nothing
    coords=_model_record_node_coords(m)
    for (dimension,tag) in targets
        record=m.discrete[(dimension,tag)]
        if dimension==1
            _parametrize_discrete_curve!(m,record,coords,caller)
        else
            _parametrize_discrete_surface!(m,record,coords,caller)
        end
    end
    return nothing
end

# ---- compute_homology ------------------------------------------------------
#
# Gmsh's `mesh.computeHomology` performs the (co)homology computations queued
# through `addHomologyRequest`, representing each generator as a new discrete
# entity holding its representative chain and wrapping it in a physical
# group. Tessella computes simplicial (co)homology over GF(2) on the
# closure of the domain's mesh cells, relative to the subdomain's cell
# closure. Quadrangles are canonically subdivided along the minimal-corner
# diagonal so every complex is simplicial.

# Canonical simplex cells (node-tag tuples, sorted) decomposing one element.
function _homology_simplices(msh_type::Int,nodes::Vector{Int32},
                             caller::AbstractString)
    family=msh_family(msh_type)
    family===:pnt && return NTuple{1,Int32}[(nodes[1],)]
    family===:lin && return NTuple{2,Int32}[Tuple(sort!(Int32[nodes[1],
                                                        nodes[2]]))]
    if family===:tri
        return NTuple{3,Int32}[Tuple(sort!(Int32[nodes[1],nodes[2],
                                                nodes[3]]))]
    end
    if family===:qua
        a,b,c,d=nodes[1],nodes[2],nodes[3],nodes[4]
        corners=Int32[a,b,c,d]
        position=argmin(corners)
        if position==1 || position==3
            return NTuple{3,Int32}[Tuple(sort!(Int32[a,b,c])),
                                   Tuple(sort!(Int32[a,c,d]))]
        end
        return NTuple{3,Int32}[Tuple(sort!(Int32[a,b,d])),
                               Tuple(sort!(Int32[b,c,d]))]
    end
    family===:tet && return NTuple{4,Int32}[Tuple(sort!(Int32[
        nodes[1],nodes[2],nodes[3],nodes[4]]))]
    throw(ArgumentError(
        "$caller: elements of type $msh_type cannot be decomposed " *
        "into simplicial homology cells"))
end

# GF(2) column reduction: eliminate `column` by registered pivot columns
# (`row => pivot column`). Returns the new pivot row, or 0 when the column
# reduced to zero.
function _gf2_eliminate!(column::Set{Int},pivots::Dict{Int,Set{Int}})
    while !isempty(column)
        row=maximum(column)
        pivot=get(pivots,row,nothing)
        pivot===nothing && return row
        symdiff!(column,pivot)
    end
    return 0
end

# Boundary faces (sorted tuples, one vertex dropped each) of a simplex cell.
function _cell_faces(cell::Tuple)
    n=length(cell)
    return Tuple[Tuple(cell[k] for k in 1:n if k!=j) for j in 1:n]
end

# Build the relative simplicial complex: closure of `domain` element
# decompositions minus the closure of `subdomain`'s. `cells` maps each entity
# to `(msh_type, node_tags)` element pairs. Returns per-dimension sorted cell
# lists (0..3).
function _homology_complex(cells::Dict{Tuple{Int,Int},
                                     Vector{Tuple{Int32,Vector{Int32}}}},
                           domain::Set{Tuple{Int,Int}},
                           subdomain::Set{Tuple{Int,Int}},
                           caller::AbstractString)
    dom=[Set{Tuple}() for _ in 0:3]
    sub=[Set{Tuple}() for _ in 0:3]
    for (key,elements) in cells
        target=key in subdomain ? sub : key in domain ? dom : nothing
        target===nothing && continue
        for (msh_type,nodes) in elements
            for simplex in _homology_simplices(Int(msh_type),nodes,caller)
                d=length(simplex)-1
                push!(target[d+1],simplex)
                stack=Tuple[simplex]
                while !isempty(stack)
                    cell=pop!(stack)
                    length(cell)==1 && continue
                    for face in _cell_faces(cell)
                        fd=length(face)-1
                        face in target[fd+1] && continue
                        push!(target[fd+1],face)
                        push!(stack,face)
                    end
                end
            end
        end
    end
    relative=[setdiff(dom[d+1],sub[d+1]) for d in 0:3]
    return [sort!(collect(relative[d+1])) for d in 0:3]
end

# GF(2) nullspace of the boundary map restricted to `cells_k` → `cells_km1`:
# each generator is a set of indices into `cells_k` forming a cycle.
function _homology_cycles(cells_km1::Vector{Tuple},cells_k::Vector{Tuple})
    index=Dict{Tuple,Int}(cell=>i for (i,cell) in enumerate(cells_km1))
    pivots=Dict{Int,Set{Int}}()
    transforms=Dict{Int,Set{Int}}()
    cycles=Vector{Set{Int}}()
    for (j,cell) in enumerate(cells_k)
        column=Set{Int}()
        for face in _cell_faces(cell)
            row=get(index,face,0)
            row!=0 && push!(column,row)
        end
        transform=Set{Int}([j])
        while !isempty(column)
            row=maximum(column)
            pivot=get(pivots,row,nothing)
            if pivot===nothing
                pivots[row]=column
                transforms[row]=transform
                break
            end
            symdiff!(column,pivot)
            symdiff!(transform,transforms[row])
        end
        isempty(column) && push!(cycles,transform)
    end
    return cycles
end

# Pivot structure of the boundary map from (k+1)-cells to k-cells, used to
# quotient cycles by boundaries.
function _homology_boundary_pivots(cells_k::Vector{Tuple},
                                   cells_kp1::Vector{Tuple})
    index=Dict{Tuple,Int}(cell=>i for (i,cell) in enumerate(cells_k))
    pivots=Dict{Int,Set{Int}}()
    for cell in cells_kp1
        column=Set{Int}()
        for face in _cell_faces(cell)
            row=get(index,face,0)
            row!=0 && push!(column,row)
        end
        while !isempty(column)
            row=maximum(column)
            pivot=get(pivots,row,nothing)
            pivot===nothing && (pivots[row]=column;break)
            symdiff!(column,pivot)
        end
    end
    return pivots
end

# Coboundary (transpose) pivots for cohomology representatives: images of
# (k-1)-cells among k-cells.
function _homology_coboundary_pivots(cells_km1::Vector{Tuple},
                                     cells_k::Vector{Tuple})
    pivots=Dict{Int,Set{Int}}()
    for cell in cells_km1
        column=Set{Int}()
        for (j,kcell) in enumerate(cells_k)
            _cell_contains(kcell,cell) && push!(column,j)
        end
        row=_gf2_eliminate!(column,pivots)
        row!=0 && (pivots[row]=column)
    end
    return pivots
end

@inline _cell_contains(cell::Tuple,face::Tuple)=
    length(face)==length(cell)-1 && issubset(face,cell)

# Reduce each cycle against boundary pivots; surviving remainders are the
# homology basis representatives.
function _homology_quotient(cycles::Vector{Set{Int}},
                            pivots::Dict{Int,Set{Int}})
    generators=Set{Int}[]
    for cycle in cycles
        remainder=copy(cycle)
        row=_gf2_eliminate!(remainder,pivots)
        row==0 && continue
        pivots[row]=remainder
        push!(generators,remainder)
    end
    return generators
end

# MSH element type for a k-simplex cell.
_homology_cell_type(cell::Tuple)=(15,1,2,4)[length(cell)]

"""
    compute_homology!(model, cells; requests=meshing.homology_requests)

Perform the queued homology requests (see `add_homology_request!`), matching
Gmsh's `mesh.computeHomology`. `cells` maps `(dim, tag)` entity keys to that
entity's `(msh_type, node_tags)` elements; the simplicial complex is the
closure of the domain entities' cells. Domain and subdomain tags refer to
physical groups of any dimension. Each surviving generator becomes a new
discrete entity of the generator's dimension holding the representative chain,
wrapped in a fresh physical group whose `(dim, tag)` is returned. `"Homology"`
computes relative homology H_k(dom, sub) over GF(2); `"Cohomology"` computes
the matching cocycle representatives on the primal complex. Returns the
`(dim, tag)` pairs of all created physical groups.
"""
function compute_homology!(m::GeoModel,cells::Dict{Tuple{Int,Int},
        Vector{Tuple{Int32,Vector{Int32}}}};
        requests=m.meshing.homology_requests)
    caller="compute_homology!"
    output=Tuple{Int,Int}[]
    for request in requests
        kind=request.kind
        kind in ("Homology","Cohomology") || throw(ArgumentError(
            "$caller: unsupported computation type \"$kind\""))
        domain=Set{Tuple{Int,Int}}()
        for (dim,tag) in keys(m.physical)
            tag in request.domain || continue
            for entity in m.physical[(dim,tag)]
                push!(domain,(dim,entity))
            end
        end
        isempty(domain) && throw(ArgumentError("$caller: domain is empty"))
        subdomain=Set{Tuple{Int,Int}}()
        for (dim,tag) in keys(m.physical)
            tag in request.subdomain || continue
            for entity in m.physical[(dim,tag)]
                push!(subdomain,(dim,entity))
            end
        end
        complex=_homology_complex(cells,domain,subdomain,caller)
        for k in request.dims
            (k isa Integer && !(k isa Bool) && 0<=k<=3) ||
                throw(ArgumentError("$caller: dims entries must be in 0:3"))
            cells_k=complex[Int(k)+1]
            isempty(cells_k) && continue
            lower=Int(k)==0 ? Tuple[] : complex[Int(k)]
            if kind=="Homology"
                cycles=_homology_cycles(lower,cells_k)
                pivots=_homology_boundary_pivots(cells_k,complex[Int(k)+2])
            else
                cycles=_homology_cocycles(cells_k,complex[Int(k)+2])
                pivots=_homology_coboundary_pivots(lower,cells_k)
            end
            for generator in _homology_quotient(cycles,pivots)
                entity_tag=add_discrete_entity!(m,Int(k),0)
                record=m.discrete[(Int(k),entity_tag)]
                for cell_index in sort!(collect(generator))
                    cell=cells_k[cell_index]
                    _record_append_element!(record,
                        _homology_cell_type(cell),
                        Int32(_model_fresh_element_tag(m)),
                        Int32[cell...])
                end
                physical=add_physical_group!(m,Int(k),[entity_tag])
                push!(output,(Int(k),physical))
            end
        end
    end
    return output
end

# Cocycles of degree k: nullspace of the coboundary δ_k : k-cells → (k+1)-cells
# (the transpose incidence). Each generator is a set of k-cell indices.
function _homology_cocycles(cells_k::Vector{Tuple},cells_kp1::Vector{Tuple})
    index=Dict{Tuple,Int}(cell=>i for (i,cell) in enumerate(cells_k))
    incidence=Dict{Int,Set{Int}}()
    for (row,kp1) in enumerate(cells_kp1)
        for face in _cell_faces(kp1)
            column=get(index,face,0)
            column==0 && continue
            push!(get!(Set{Int},incidence,column),row)
        end
    end
    pivots=Dict{Int,Set{Int}}()
    transforms=Dict{Int,Set{Int}}()
    cycles=Vector{Set{Int}}()
    for j in eachindex(cells_k)
        column=copy(get(incidence,j,Set{Int}()))
        transform=Set{Int}([j])
        while !isempty(column)
            row=maximum(column)
            pivot=get(pivots,row,nothing)
            if pivot===nothing
                pivots[row]=column
                transforms[row]=transform
                break
            end
            symdiff!(column,pivot)
            symdiff!(transform,transforms[row])
        end
        isempty(column) && push!(cycles,transform)
    end
    return cycles
end

# ---- generator consumption -------------------------------------------------
#
# The attribute setters above record state; the functions below are the
# canonical paths the generators consume them through.

# Uniform parameter spacing — the distribution Gmsh's `F_Transfinite` emits
# when the density collapses to a constant (`coef <= 0`, `coef == 1`, Beta with
# `coef < 1`, or an unknown signed type falling through to `val = 1.`).
function _transfinite_uniform_parameters(num_nodes::Int)
    parameters=Vector{Float64}(undef,num_nodes)
    for i in 0:num_nodes-1
        parameters[i+1]=i/(num_nodes-1)
    end
    return parameters
end

# Normalized parameters of the `num_nodes` boundary nodes a transfinite curve
# produces, delegated to the differentially validated `TransfiniteCurve` laws.
#
# Gmsh stores a *signed* transfinite type (negative when the `.geo` tag was
# negative — `reversed` here) and a raw coefficient on the GEO record, and
# `F_Transfinite` maps them to a density: Progression/Power spaces
# geometrically with ratio `coef` (sign flips the ratio to `1/coef`), Bump
# clusters at both ends (type sign ignored), Beta applies the beta law
# (negative type mirrors the argument), and the uniform fallbacks above apply
# to the raw `coef`. The HWall kinds (types 5-7) first transform `coef` — a
# signed first-layer wall height — into an ordinary coefficient through Gmsh's
# `newton_get_r`/`bissection_get_*` solves; that transform reads the signed
# type, so a reversed HWall record (negative type) never reaches it and falls
# to the unknown-type warning path — a uniform distribution — like the
# grammar-only `Beta_Symmetrical`/`Beta_Symmetrical_HWall` kinds (types 8/9).
# HWall needs the curve's geometric length, so only straight curves support it.
function _transfinite_parameters(m::GeoModel,num_nodes::Int,kind::Symbol,
                                 coef::Float64,caller::AbstractString,curve::Int;
                                 reversed::Bool=false)
    if !reversed && kind in _TRANSFINITE_HWALL_KINDS
        law=kind===:progression_hwall ? :progression :
            kind===:bump_hwall ? :bump : :beta
        # `F_Transfinite` reads `fabs(coef)` as the wall height and uses the
        # coefficient's sign as the wall side (Bump is symmetric regardless).
        return transfinite_curve_hwall(num_nodes;mesh_type=law,
            wall_height=abs(coef),
            curve_length=_model_curve_length(m,curve,caller),
            orientation=(coef<0 && law!==:bump) ? :end : :start)
    end
    if kind in (:progression,:bump,:beta) &&
       !(coef<=0.0 || coef==1.0 || (kind===:beta && coef<1.0))
        # The module maps a negative coefficient to the reversed distribution,
        # which is exactly what a negative signed type does in `F_Transfinite`.
        return transfinite_curve_parameters(num_nodes;mesh_type=kind,
            coefficient=reversed ? -coef : coef)
    end
    return _transfinite_uniform_parameters(num_nodes)
end

# Seed `forced` curve-parameter lists and the `(curve, parameter) => size`
# override map from stored transfinite and size-at-parametric-point
# attributes. `forced` also carries periodic synchronization entries, which
# take precedence for curves already constrained.
function _attribute_forced_parameters(m::GeoModel,t::Int,
                                      forced::Dict{Int,Vector{Float64}},
                                      caller::AbstractString)
    boundary_curves=Set{Int}()
    for loop_id in m.surfaces[t]
        for signed in m.loops[loop_id]
            push!(boundary_curves,abs(signed))
        end
    end
    for (edim,etag) in get(m.embeds,(2,t),NTuple{2,Int}[])
        edim==1 && push!(boundary_curves,etag)
    end
    param_sizes=Dict{Tuple{Int,Float64},Float64}()
    for (curve,spec) in m.meshing.transfinite_curves
        curve in boundary_curves || continue
        # A `Degenerated` curve collapses to a single edge before the
        # transfinite branch is considered (Gmsh `meshGEdge` checks
        # `degenerate(0)` first), so it overrides every other attribute.
        curve in m.meshing.degenerated && continue
        haskey(forced,curve) && throw(ArgumentError(
            "$caller: Curve[$curve] has both periodic and transfinite " *
            "constraints"))
        forced[curve]=_transfinite_parameters(
            m,spec.num_nodes,spec.kind,spec.coef,caller,curve;
            reversed=spec.reversed)
    end
    for ((edim,etag),entries) in m.meshing.size_at_params
        edim==1 || continue
        etag in boundary_curves || continue
        etag in m.meshing.degenerated && continue
        list=get!(forced,etag,Float64[])
        for (param,size) in entries
            p=only(param)
            p in list || push!(list,p)
            param_sizes[(etag,p)]=size
        end
        sort!(list)
    end
    for (curve,list) in forced
        curve in boundary_curves || continue
        (first(list)==0.0 && last(list)==1.0) || throw(ArgumentError(
            "$caller: Curve[$curve] parameter list must span [0,1]"))
    end
    return param_sizes
end
