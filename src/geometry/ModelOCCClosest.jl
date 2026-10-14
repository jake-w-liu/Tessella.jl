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
# `seen` keeps shared subshapes idempotent while the projection units —
# (dim, tag) pairs `_dist_entity_project` answers — always enumerate:
# itself for dims 0:3, member edges for a wire, member faces for a shell.
# Re-enumerating a seen shape therefore reports units that pair queries
# (e.g. a wire against one of its own edges) still need.
# `snap` (when given) is a per-materialization list of canonical point
# coordinates: every materialized vertex stores the first coords seen in
# its near-coincidence cluster, so kernel vertices that differ only by
# evaluation dust (~ulp, e.g. engine split points recomputed through a
# different arithmetic path) become bitwise-equal in the scratch model.
# The engine's exact `orient3` coplanarity checks then hold; the kernel
# records themselves keep their own coords. Coordinates closer than the
# dust bound are treated as coincident — distinct points within that
# distance are not representable faithfully anyway.
# Canonicalize the last mantissa bit: section arithmetic can leave a
# corner a few ulps off the grid value its siblings store exactly, and a
# singleton dust vertex has no coincident partner to canonicalize it.
# Fifteen-significant-digit rounding collapses that dust while moving
# any coordinate by well under the geometric coincidence tolerance.
_model_occ_snap_coord(c::Float64)=round(c,sigdigits=15)

function _model_occ_snap_point(canon::Vector{NTuple{3,Float64}},
                               p::NTuple{3,Float64})
    p=_model_occ_snap_coord.(p)
    atol=1e-11*max(1.0,abs(p[1]),abs(p[2]),abs(p[3]))
    for c in canon
        (abs(p[1]-c[1])<=atol && abs(p[2]-c[2])<=atol &&
         abs(p[3]-c[3])<=atol) && return c
    end
    push!(canon,p)
    return p
end

function _model_occ_dist_materialize!(sm::GeoModel,kernel,dim,id,
                                      seen::Set{UInt64},caller;
                                      snap::Union{Nothing,
                                          Vector{NTuple{3,Float64}}}=nothing,
                                      untyped_defaults::Bool=false)
    dim==0 && begin
        tag=_model_occ_shape_tag(kernel,id,caller)
        id in seen || begin
            push!(seen,id)
            p=kernel.shapes[id].geometry.point
            snap===nothing || (p=_model_occ_snap_point(snap,p))
            sm.points[tag]=p
        end
        return [(0,tag)]
    end
    if dim==1
        tag=_model_occ_shape_tag(kernel,id,caller)
        if !(id in seen)
            push!(seen,id)
            shape=kernel.shapes[id]
            a,b=shape.children
            sm.curves[tag]=(_model_occ_shape_tag(kernel,a,caller),
                            _model_occ_shape_tag(kernel,b,caller))
            # `untyped_defaults` leaves :line curves on the implicit
            # default kind like native entities (callers that evaluate
            # sections through the untyped path take it); other kinds
            # always carry their records.
            (untyped_defaults && shape.geometry.kind===:line) || begin
                sm.curve_types[tag]=shape.geometry.kind
                sm.curve_geometry[tag]=shape.geometry.geometry
            end
            _model_occ_dist_materialize!(sm,kernel,0,a,seen,caller;
                snap=snap,untyped_defaults=untyped_defaults)
            _model_occ_dist_materialize!(sm,kernel,0,b,seen,caller;
                snap=snap,untyped_defaults=untyped_defaults)
        end
        return [(1,tag)]
    end
    if dim==-1
        units=Tuple{Int,Int}[]
        for (edge_id,_) in kernel.shapes[id].geometry.uses
            append!(units,_model_occ_dist_materialize!(
                sm,kernel,1,edge_id,seen,caller;
                snap=snap,untyped_defaults=untyped_defaults))
        end
        push!(seen,id)
        return units
    end
    if dim==2
        tag=_model_occ_shape_tag(kernel,id,caller)
        if !(id in seen)
            push!(seen,id)
            shape=kernel.shapes[id]
            # Edges must be materialized before `_sort_curve_loop` reads
            # `sm.curves` endpoint pairs.
            for wire_id in shape.children, (edge_id,_) in
                kernel.shapes[wire_id].geometry.uses
                _model_occ_dist_materialize!(sm,kernel,1,edge_id,seen,
                    caller;snap=snap,untyped_defaults=untyped_defaults)
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
            (untyped_defaults && shape.geometry.kind===:plane) || begin
                sm.surface_types[tag]=shape.geometry.kind
                sm.surface_geometry[tag]=shape.geometry.geometry
            end
        end
        return [(2,tag)]
    end
    if dim==-2
        units=Tuple{Int,Int}[]
        for (face_id,_) in kernel.shapes[id].geometry.uses
            append!(units,_model_occ_dist_materialize!(
                sm,kernel,2,face_id,seen,caller;
                snap=snap,untyped_defaults=untyped_defaults))
        end
        push!(seen,id)
        return units
    end
    dim==3 || throw(ArgumentError("$caller: unsupported OCC dimension $dim"))
    tag=_model_occ_shape_tag(kernel,id,caller)
    if !(id in seen)
        push!(seen,id)
        shape=kernel.shapes[id]
        loops=Int[]
        for shell_id in shape.children
            signed=Int[(reverse ? -1 : 1)*_model_occ_shape_tag(
                           kernel,fid,caller)
                       for (fid,reverse) in
                           kernel.shapes[shell_id].geometry.uses]
            for (fid,_) in kernel.shapes[shell_id].geometry.uses
                _model_occ_dist_materialize!(sm,kernel,2,fid,seen,caller;
                    snap=snap,untyped_defaults=untyped_defaults)
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

# `gmsh.model.occ.get_distance` — `OCC_Internals::getDistance`: bound checks
# on both inputs, then `BRepExtrema_DistShapeShape(shape1, shape2)` yields
# the distance and one attaining point on each shape. Wire and shell inputs
# expand to their member units and the pairwise minimum stands in for the
# shape-level extrema. Volume queries follow the codebase's documented
# convention (inside → 0, the point pair at the interior rep) with a
# boundary-minimum fallback when ray containment is undecidable — matching
# upstream's raw extrema answer in the same corner. An unconverged result
# reports distance −1 like upstream's `!ret → distance = -1`.
const _OCC_DIST_FAILURE=(-1.0,0.0,0.0,0.0,0.0,0.0,0.0)

function _model_occ_dist_pair_units(sm::GeoModel,a,b,lc,caller)
    swapped=a>b
    da,ta,db,tb=swapped ? (b[1],b[2],a[1],a[2]) : (a[1],a[2],b[1],b[2])
    got=nothing
    if db==3
        got=try
            if da==3
                _dist_pair(sm,3,ta,3,tb,lc,caller)
            else
                _dist_pair(sm,da,ta,3,tb,lc,caller)
            end
        catch error
            error isa InterruptException && rethrow()
            error isa ErrorException || rethrow()
            occursin("degenerate",sprint(showerror,error)) || rethrow()
            da==3 ? _dist_shells_min(sm,ta,tb,lc,caller) :
                    _dist_volume_boundary_min(sm,tb,da,ta,lc,caller)
        end
    else
        got=_dist_pair(sm,da,ta,db,tb,lc,caller)
    end
    sqd,pa,pb=got
    return swapped ? (sqd,pb,pa) : (sqd,pa,pb)
end

function _model_occ_get_distance(m::GeoModel,dim1,tag1,dim2,tag2)
    caller="OCC.get_distance"
    kernel=m.occ_kernel
    d1=_model_occ_dimension(dim1,caller)
    t1=_model_occ_requested_tag(tag1,caller)
    id1=get(kernel.bindings,(d1,t1),nothing)
    id1===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE entity of dimension $d1 with tag $t1"))
    d2=_model_occ_dimension(dim2,caller)
    t2=_model_occ_requested_tag(tag2,caller)
    id2=get(kernel.bindings,(d2,t2),nothing)
    id2===nothing && throw(ArgumentError(
        "$caller: Unknown OpenCASCADE entity of dimension $d2 with tag $t2"))
    sm=GeoModel();seen=Set{UInt64}()
    ua=_model_occ_dist_materialize!(sm,kernel,d1,id1,seen,caller)
    ub=_model_occ_dist_materialize!(sm,kernel,d2,id2,seen,caller)
    lc=_model_lc(sm)
    best=Inf;pa=pb=(0.0,0.0,0.0)
    for a in ua, b in ub
        got=try
            _model_occ_dist_pair_units(sm,a,b,lc,caller)
        catch error
            error isa InterruptException && rethrow()
            error isa ErrorException ? nothing : rethrow()
        end
        got===nothing && continue
        sqd,x,y=got
        sqd<best && (best=sqd;pa=x;pb=y)
    end
    (isfinite(best) && all(isfinite,pa) && all(isfinite,pb)) ||
        return _OCC_DIST_FAILURE
    return sqrt(best),pa[1],pa[2],pa[3],pb[1],pb[2],pb[3]
end
