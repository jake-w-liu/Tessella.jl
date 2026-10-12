# `gmsh.model.occ.getClosestEntities` — `OCC_Internals::getClosestEntities`.
# Upstream resolves each bound dim-tag, classifies the query point inside
# solid inputs (`BRepClass3d_SolidClassifier`, IN/ON → zero distance at the
# query point), otherwise takes `BRepExtrema_DistShapeShape`'s point on the
# shape, then sorts the (distance, (dim, tag, x, y, z)) keys through an
# `std::set` — duplicate inputs collapse and ties order lexicographically —
# and emits the first `n` (`n <= 0` still emits one) sorted results.
#
# Tessella's pending OCC shapes live in the kernel registry, not the public
# model maps the `_dist_*` projectors read. The materializer below writes the
# queried entities — with their transitive wire/edge/vertex closure and the
# same sign and hole-reversal conventions `_model_occ_synchronize!` uses —
# into a private scratch `GeoModel`, so the verified model-side projection
# and classification machinery applies verbatim without mutating public
# state.

# Materialize one bound shape (transitively) into the scratch model.
# `seen` keeps shared subshapes idempotent. Returns the projection units —
# (dim, tag) pairs `_dist_entity_project` answers — for the input: itself
# for dims 0:3, member edges for a wire, member faces for a shell.
function _model_occ_dist_materialize!(sm::GeoModel,kernel,dim,id,
                                      seen::Set{UInt64},caller)
    dim==0 && begin
        id in seen && return Tuple{Int,Int}[]
        push!(seen,id)
        tag=_model_occ_shape_tag(kernel,id,caller)
        sm.points[tag]=kernel.shapes[id].geometry.point
        return [(0,tag)]
    end
    if dim==1
        id in seen && return Tuple{Int,Int}[]
        push!(seen,id)
        shape=kernel.shapes[id]
        tag=_model_occ_shape_tag(kernel,id,caller)
        a,b=shape.children
        sm.curves[tag]=(_model_occ_shape_tag(kernel,a,caller),
                        _model_occ_shape_tag(kernel,b,caller))
        sm.curve_types[tag]=shape.geometry.kind
        sm.curve_geometry[tag]=shape.geometry.geometry
        _model_occ_dist_materialize!(sm,kernel,0,a,seen,caller)
        _model_occ_dist_materialize!(sm,kernel,0,b,seen,caller)
        return [(1,tag)]
    end
    if dim==-1
        id in seen && return Tuple{Int,Int}[]
        push!(seen,id)
        units=Tuple{Int,Int}[]
        for (edge_id,_) in kernel.shapes[id].geometry.uses
            append!(units,_model_occ_dist_materialize!(
                sm,kernel,1,edge_id,seen,caller))
        end
        return units
    end
    if dim==2
        id in seen && return Tuple{Int,Int}[]
        push!(seen,id)
        shape=kernel.shapes[id]
        tag=_model_occ_shape_tag(kernel,id,caller)
        # Edges must be materialized before `_sort_curve_loop` reads
        # `sm.curves` endpoint pairs.
        for wire_id in shape.children, (edge_id,_) in
            kernel.shapes[wire_id].geometry.uses
            _model_occ_dist_materialize!(sm,kernel,1,edge_id,seen,caller)
        end
        loops=Int[]
        for (index,wire_id) in enumerate(shape.children)
            wire=kernel.shapes[wire_id]
            signed=Int[(reverse ? -1 : 1)*_model_occ_shape_tag(
                           kernel,edge,caller)
                       for (edge,reverse) in wire.geometry.uses]
            index>1 && (signed=[-s for s in Iterators.reverse(signed)])
            ordered=_sort_curve_loop(sm,signed,caller)
            loop=_alloc_curve_loop_tag(sm,0,caller)
            sm.loops[loop]=ordered
            push!(loops,loop)
        end
        sm.surfaces[tag]=loops
        sm.surface_types[tag]=shape.geometry.kind
        sm.surface_geometry[tag]=shape.geometry.geometry
        return [(2,tag)]
    end
    if dim==-2
        id in seen && return Tuple{Int,Int}[]
        push!(seen,id)
        units=Tuple{Int,Int}[]
        for (face_id,_) in kernel.shapes[id].geometry.uses
            append!(units,_model_occ_dist_materialize!(
                sm,kernel,2,face_id,seen,caller))
        end
        return units
    end
    dim==3 || throw(ArgumentError("$caller: unsupported OCC dimension $dim"))
    id in seen && return Tuple{Int,Int}[]
    push!(seen,id)
    shape=kernel.shapes[id]
    tag=_model_occ_shape_tag(kernel,id,caller)
    loops=Int[]
    for shell_id in shape.children
        signed=Int[(reverse ? -1 : 1)*_model_occ_shape_tag(
                       kernel,fid,caller)
                   for (fid,reverse) in
                       kernel.shapes[shell_id].geometry.uses]
        for (fid,_) in kernel.shapes[shell_id].geometry.uses
            _model_occ_dist_materialize!(sm,kernel,2,fid,seen,caller)
        end
        loop=_alloc_surface_loop_tag(sm,0,caller)
        sm.surface_loops[loop]=signed
        push!(loops,loop)
    end
    sm.volumes[tag]=loops
    primitive=get(shape.geometry,:primitive,nothing)
    if primitive!==nothing
        primitive.occ===:cylinder && (sm.cylinders[tag]=primitive.record)
        primitive.occ===:sphere && (sm.spheres[tag]=primitive.record)
        primitive.occ===:cone && (sm.cones[tag]=primitive.record)
    end
    return [(3,tag)]
end

# Solid distance: `_dist_volume_project` minus the degenerate-containment
# throw — inside → (0, q) exactly, otherwise the boundary-face minimum.
# Upstream's `DistShapeShape` never needs containment, so an undecidable
# ray classification still yields the boundary answer instead of
# skipping the entity.
function _model_occ_dist_solid(sm::GeoModel,tag::Int,q,lc,caller)
    inside=_dist_point_in_volume(sm,tag,q,lc,caller)
    inside===true && return 0.0,q
    best=Inf;best_point=(0.0,0.0,0.0)
    for stag in _dist_volume_faces(sm,tag,caller)
        dist,pt=_dist_face_project(sm,stag,q,lc,caller)
        dist<best && (best=dist;best_point=pt)
    end
    isfinite(best) || throw(ErrorException(
        "$caller: Volume[$tag] admits no boundary closest point"))
    return best,best_point
end

# Minimum squared distance and closest point over an entity's projection
# units. A `nothing` return mirrors upstream's `dist.IsDone()` skip.
function _model_occ_dist_units(sm::GeoModel,units,q,lc,caller)
    best=Inf;best_point=(0.0,0.0,0.0)
    for (ud,ut) in units
        sqd,pt=ud==3 ? _model_occ_dist_solid(sm,ut,q,lc,caller) :
                       _dist_entity_project(sm,ud,ut,q,lc,caller)
        sqd<best && (best=sqd;best_point=pt)
    end
    isfinite(best) || return nothing
    return best,best_point
end

function _model_occ_get_closest_entities(m::GeoModel,x,y,z,dim_tags,n)
    caller="OCC.get_closest_entities"
    q=_finite3(x,y,z,caller)
    n isa Integer && !(n isa Bool) || throw(ArgumentError(
        "$caller: n must be an integer"))
    kernel=m.occ_kernel
    entries=Tuple{Int,Int,UInt64}[]
    for entry in dim_tags
        d=_model_occ_dimension(first(entry),caller)
        t=_model_occ_requested_tag(last(entry),caller)
        id=get(kernel.bindings,(d,t),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $d " *
            "with tag $t"))
        push!(entries,(d,t,id))
    end
    sm=GeoModel()
    seen=Set{UInt64}()
    units=Vector{Tuple{Int,Int}}[]
    for (d,_,id) in entries
        push!(units,_model_occ_dist_materialize!(
            sm,kernel,d,id,seen,caller))
    end
    lc=_model_lc(sm)
    results=Set{Tuple{Float64,Tuple{Int,Int,Float64,Float64,Float64}}}()
    for (index,(d,t,_)) in enumerate(entries)
        got=try
            _model_occ_dist_units(sm,units[index],q,lc,caller)
        catch error
            error isa InterruptException && rethrow()
            # Upstream skips entities whose `DistShapeShape` does not
            # converge; the projectors signal that path with ErrorException.
            error isa ErrorException ? nothing : rethrow()
        end
        got===nothing && continue
        sqd,pt=got
        distance=sqrt(sqd)
        all(isfinite,pt) && isfinite(distance) || continue
        push!(results,(distance,(d,t,pt[1],pt[2],pt[3])))
    end
    isempty(results) && return Tuple{Int,Int}[],Float64[],Float64[]
    sorted=sort!(collect(results))
    out_dim_tags=Tuple{Int,Int}[];distances=Float64[];coord=Float64[]
    for (distance,(d,t,px,py,pz)) in Iterators.take(sorted,max(Int(n),1))
        push!(out_dim_tags,(d,t))
        push!(distances,distance)
        append!(coord,(px,py,pz))
    end
    return out_dim_tags,distances,coord
end
