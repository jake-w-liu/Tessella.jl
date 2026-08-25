"""
    SizeField

Scalar and mesh-size fields.  The field graph follows Gmsh's composition model:
generic scalar fields such as [`DistanceField`](@ref) feed size-producing fields
such as [`ThresholdField`](@ref), while [`MinSize`](@ref), [`MaxSize`](@ref) and
[`BoundedSize`](@ref) compose or constrain the result.

`field_value(field, x, y, z)` evaluates any scalar field.  `size_at(field, ...)`
is deliberately restricted to `AbstractSizeField` and guarantees a finite,
strictly-positive mesh size before a meshing kernel can consume it.
"""
module SizeField

using ..MeshTypes: Mesh, validate
using ..IO: GeoParams, GeoFieldSpec, _geo_split_list

export AbstractField, AbstractSizeField, AbstractAnisoField, ConstantSize, FunctionSize
export DistanceField, ThresholdField, BoxField, BallField, CylinderField, FrustumField
export MinSize, MaxSize, BoundedSize
export MinField, MaxField, field_value, size_at, metric_at, build_geo_size_field, GMSH_MAX_SIZE
export Metric3, isotropic_metric, metric_size, metric_eigenvalues
export directional_size, metric_edge_length
export MathEvalField, MathEvalAnisoField, GradientField, LaplacianField, MeanField
export CurvatureField, CurvatureOpField, MaxEigenHessianField, LonLatField, ParametricField
export StructuredField, RestrictField, ConstantField, ExtendField, OctreeField, PostViewField
export MinAnisoField, IntersectAnisoField, AttractorAnisoCurveField
export BoundaryLayerField, AutomaticMeshSizeField, ExternalProcessField
export parse_matheval, eval_matheval, intersection_alauzet, intersection_conserve_mostaniso
export build_geo_boundary_layer_fields

"""Finite sentinel used by Gmsh fields for "do not constrain the mesh size"."""
const GMSH_MAX_SIZE = 1.0e22

"""Base type for scalar fields evaluable with [`field_value`](@ref)."""
abstract type AbstractField end

"""
Base type for scalar fields that meshing kernels may consume through
[`size_at`](@ref). A field can still produce nonpositive intermediate values when
used in a Gmsh-compatible expression graph; `size_at` enforces the finite,
strictly-positive contract at every actual mesh query.
"""
abstract type AbstractSizeField <: AbstractField end

"""Base type for fields that provide a symmetric positive-definite metric."""
abstract type AbstractAnisoField <: AbstractSizeField end

"""
    field_value(field, x, y, z)
    field_value(field, point)

Evaluate a scalar field. Field implementations reject non-finite results; generic
scalar fields may legitimately return zero or a negative value.
"""
function field_value end

@inline function _float_value(x::Real, caller::AbstractString, name::AbstractString)
    x isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    y = try
        Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be representable as Float64"))
    end
    isfinite(y) || throw(ArgumentError("$caller: $name must be finite (got $x)"))
    return y
end

@inline function _positive_value(x::Real, caller::AbstractString, name::AbstractString)
    y = _float_value(x, caller, name)
    y > 0 || throw(ArgumentError("$caller: $name must be positive (got $x)"))
    return y
end

@inline function _point3(p::NTuple{3,Float64}, caller::AbstractString)
    all(isfinite,p) || throw(ArgumentError(
        "$caller: point coordinates must be finite (got $p)"))
    return p
end

@inline function _point3(p, caller::AbstractString)
    n = try
        length(p)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: a point must be an indexable coordinate collection"))
    end
    n >= 2 || throw(ArgumentError("$caller: a point needs at least two coordinates"))
    raw = try
        (p[1], p[2], n >= 3 ? p[3] : 0.0)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: point coordinates must be Float64-representable"))
    end
    any(value -> value isa Bool,raw) && throw(ArgumentError(
        "$caller: point coordinates must not be Bool"))
    q = try
        (Float64(raw[1]),Float64(raw[2]),Float64(raw[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: point coordinates must be Float64-representable"))
    end
    all(isfinite, q) || throw(ArgumentError("$caller: point coordinates must be finite (got $q)"))
    return q
end

@inline function _checked_field_result(value, caller::AbstractString, x, y, z;
                                       positive::Bool=false)
    value isa Real ||
        throw(ArgumentError("$caller: field returned $(typeof(value)), not a real value, at ($x,$y,$z)"))
    value isa Bool && throw(ArgumentError(
        "$caller: field returned Bool, not a numeric value, at ($x,$y,$z)"))
    v = try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: field result $value is not Float64-representable at ($x,$y,$z)"))
    end
    isfinite(v) || throw(ArgumentError("$caller: field returned a non-finite value $value at ($x,$y,$z)"))
    positive && v <= 0 &&
        throw(ArgumentError("$caller: field returned a non-positive mesh size $value at ($x,$y,$z)"))
    return v
end

"""Uniform target size `h` everywhere."""
struct ConstantSize <: AbstractSizeField
    h::Float64
    ConstantSize(h::Real) = new(_positive_value(h, "ConstantSize", "h"))
end

"""Mesh size from an arbitrary callable `f(x,y,z) -> h`."""
struct FunctionSize{F} <: AbstractSizeField
    f::F
end

@inline field_value(sf::ConstantSize, x, y, z) = sf.h
@inline function field_value(sf::FunctionSize, x, y, z)
    return _checked_field_result(sf.f(x, y, z), "FunctionSize", x, y, z; positive=true)
end

# ── Gmsh-compatible geometric fields ──────────────────────────────────────────

function _point_matrix(values, caller::AbstractString)
    if values isa AbstractMatrix
        size(values, 1) in (2, 3) ||
            throw(ArgumentError("$caller: point matrix must have 2 or 3 rows"))
        out = Matrix{Float64}(undef, 3, size(values, 2))
        @inbounds for j in axes(values, 2)
            p = _point3(view(values, :, j), caller)
            out[1,j]=p[1]; out[2,j]=p[2]; out[3,j]=p[3]
        end
        return out
    end
    vals = collect(values)
    out = Matrix{Float64}(undef, 3, length(vals))
    @inbounds for (j, value) in enumerate(vals)
        p = _point3(value, caller)
        out[1,j]=p[1]; out[2,j]=p[2]; out[3,j]=p[3]
    end
    return out
end

function _primitive_matrix(values, arity::Int, caller::AbstractString)
    vals = collect(values)
    out = Matrix{Float64}(undef, 3arity, length(vals))
    @inbounds for (j, value) in enumerate(vals)
        length(value) == arity ||
            throw(ArgumentError("$caller: each primitive needs exactly $arity points"))
        for k in 1:arity
            p = _point3(value[k], caller)
            off = 3(k-1)
            out[off+1,j]=p[1]; out[off+2,j]=p[2]; out[off+3,j]=p[3]
        end
    end
    return out
end

"""
    DistanceField(; points=(), segments=(), triangles=())
    DistanceField(mesh; include_points=true, include_segments=true, include_triangles=true)

Euclidean distance to the union of discrete point, line-segment and triangle
targets.  Curves and surfaces in Gmsh's `Distance` field are internally
discretized before querying; Tessella accepts that discrete representation
directly and evaluates segments and triangles exactly instead of replacing them
with a second point-cloud approximation.

Unlike a mesh-size field, a distance field can legitimately return zero.  Feed it
through [`ThresholdField`](@ref) or [`BoundedSize`](@ref) before passing it to a
mesher.
"""
struct DistanceField <: AbstractField
    points::Matrix{Float64}       # 3 × n
    segments::Matrix{Float64}     # 6 × n: a, b
    triangles::Matrix{Float64}    # 9 × n: a, b, c
    bvh_order::Vector{Int}
    bvh_lo::Vector{NTuple{3,Float64}}
    bvh_hi::Vector{NTuple{3,Float64}}
    bvh_left::Vector{Int}
    bvh_right::Vector{Int}
    bvh_first::Vector{Int}
    bvh_count::Vector{Int}
    function DistanceField(; points=(), segments=(), triangles=())
        p = _point_matrix(points, "DistanceField")
        s = _primitive_matrix(segments, 2, "DistanceField")
        t = _primitive_matrix(triangles, 3, "DistanceField")
        size(p,2)+size(s,2)+size(t,2) > 0 ||
            throw(ArgumentError("DistanceField: need at least one point, segment, or triangle target"))
        tree=_build_distance_bvh(p,s,t)
        new(p,s,t,tree...)
    end
end

const _DISTANCE_BVH_LEAF_SIZE=8

function _build_distance_bvh(points,segments,triangles)
    np=size(points,2);ns=size(segments,2);nt=size(triangles,2)
    n=try
        Base.checked_add(Base.checked_add(np,ns),nt)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("DistanceField: primitive count overflows the platform Int limit"))
    end
    primitive_lo=Matrix{Float64}(undef,3,n)
    primitive_hi=Matrix{Float64}(undef,3,n)
    centroid=Matrix{Float64}(undef,3,n)
    @inbounds for id in 1:n
        if id<=np
            j=id
            for k in 1:3
                primitive_lo[k,id]=points[k,j];primitive_hi[k,id]=points[k,j]
            end
        elseif id<=np+ns
            j=id-np
            for k in 1:3
                a=segments[k,j];b=segments[k+3,j]
                primitive_lo[k,id]=min(a,b);primitive_hi[k,id]=max(a,b)
            end
        else
            j=id-np-ns
            for k in 1:3
                a=triangles[k,j];b=triangles[k+3,j];c=triangles[k+6,j]
                primitive_lo[k,id]=min(a,b,c);primitive_hi[k,id]=max(a,b,c)
            end
        end
        for k in 1:3
            lo=primitive_lo[k,id];hi=primitive_hi[k,id]
            centroid[k,id]=signbit(lo)==signbit(hi) ? lo+(hi-lo)/2 : lo/2+hi/2
        end
    end
    order=collect(1:n)
    lows=NTuple{3,Float64}[];highs=NTuple{3,Float64}[]
    left=Int[];right=Int[];firsts=Int[];counts=Int[]
    function range_bounds(first,last)
        lo=(Inf,Inf,Inf);hi=(-Inf,-Inf,-Inf)
        @inbounds for pos in first:last
            id=order[pos]
            lo=ntuple(k -> min(lo[k],primitive_lo[k,id]),3)
            hi=ntuple(k -> max(hi[k],primitive_hi[k,id]),3)
        end
        return lo,hi
    end
    function build!(first,last)
        node=length(left)+1
        push!(lows,(0.0,0.0,0.0));push!(highs,(0.0,0.0,0.0))
        push!(left,0);push!(right,0);push!(firsts,0);push!(counts,0)
        count=last-first+1
        if count<=_DISTANCE_BVH_LEAF_SIZE
            lows[node],highs[node]=range_bounds(first,last)
            firsts[node]=first;counts[node]=count
            return node
        end
        clo=(Inf,Inf,Inf);chi=(-Inf,-Inf,-Inf)
        @inbounds for pos in first:last
            id=order[pos]
            clo=ntuple(k -> min(clo[k],centroid[k,id]),3)
            chi=ntuple(k -> max(chi[k],centroid[k,id]),3)
        end
        ext=ntuple(k -> chi[k]-clo[k],3)
        axis=ext[1]>=ext[2] ? (ext[1]>=ext[3] ? 1 : 3) :
             (ext[2]>=ext[3] ? 2 : 3)
        sort!(view(order,first:last);by=id -> centroid[axis,id],alg=QuickSort)
        mid=(first+last)>>>1
        l=build!(first,mid);r=build!(mid+1,last)
        left[node]=l;right[node]=r
        lows[node]=ntuple(k -> min(lows[l][k],lows[r][k]),3)
        highs[node]=ntuple(k -> max(highs[l][k],highs[r][k]),3)
        return node
    end
    build!(1,n)
    return order,lows,highs,left,right,firsts,counts
end

function DistanceField(mesh::Mesh; include_points::Bool=true,
                       include_segments::Bool=true, include_triangles::Bool=true)
    diagnostic=validate(mesh;require_positive_tets=false,require_manifold_tris=false)
    diagnostic.ok || throw(ArgumentError(
        "DistanceField: invalid source mesh: $(join(diagnostic.messages, "; "))"))
    points = include_points ? [Tuple(mesh.coords[:,i]) for i in axes(mesh.coords,2)] : ()
    segments = include_segments ?
        [(Tuple(mesh.coords[:,mesh.segs[1,i]]), Tuple(mesh.coords[:,mesh.segs[2,i]]))
         for i in axes(mesh.segs,2)] : ()
    triangles = include_triangles ?
        [(Tuple(mesh.coords[:,mesh.tris[1,i]]), Tuple(mesh.coords[:,mesh.tris[2,i]]),
          Tuple(mesh.coords[:,mesh.tris[3,i]])) for i in axes(mesh.tris,2)] : ()
    return DistanceField(; points=points, segments=segments, triangles=triangles)
end

@inline _sub3(a,b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _addscaled3(a,b,t) = (a[1]+t*b[1], a[2]+t*b[2], a[3]+t*b[3])
@inline _dot3(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _norm3(a) = hypot(a[1],a[2],a[3])

@inline function _distance_segment_kernel(p, a, b)
    ab = _sub3(b,a)
    den = _dot3(ab,ab)
    den == 0 && return _norm3(_sub3(p,a))
    t = clamp(_dot3(_sub3(p,a),ab)/den, zero(den), one(den))
    return _norm3(_sub3(p,_addscaled3(a,ab,t)))
end

function _distance_segment(p, a, b)
    d = _distance_segment_kernel(p,a,b)
    isfinite(d) && return min(Float64(d), GMSH_MAX_SIZE)
    return setprecision(BigFloat, 256) do
        pb=map(BigFloat,p); ab=map(BigFloat,a); bb=map(BigFloat,b)
        min(Float64(_distance_segment_kernel(pb,ab,bb)), GMSH_MAX_SIZE)
    end
end

function _distance_triangle_kernel(p, a, b, c)
    ab=_sub3(b,a); ac=_sub3(c,a)
    cr=(ab[2]*ac[3]-ab[3]*ac[2], ab[3]*ac[1]-ab[1]*ac[3],
        ab[1]*ac[2]-ab[2]*ac[1])
    if _dot3(cr,cr) == 0
        return min(_distance_segment_kernel(p,a,b), _distance_segment_kernel(p,b,c),
                   _distance_segment_kernel(p,c,a))
    end
    ap=_sub3(p,a); d1=_dot3(ab,ap); d2=_dot3(ac,ap)
    (d1 <= 0 && d2 <= 0) && return _norm3(ap)
    bp=_sub3(p,b); d3=_dot3(ab,bp); d4=_dot3(ac,bp)
    (d3 >= 0 && d4 <= d3) && return _norm3(bp)
    vc=d1*d4-d3*d2
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v=d1/(d1-d3)
        return _norm3(_sub3(p,_addscaled3(a,ab,v)))
    end
    cp=_sub3(p,c); d5=_dot3(ab,cp); d6=_dot3(ac,cp)
    (d6 >= 0 && d5 <= d6) && return _norm3(cp)
    vb=d5*d2-d1*d6
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w=d2/(d2-d6)
        return _norm3(_sub3(p,_addscaled3(a,ac,w)))
    end
    va=d3*d6-d5*d4
    if va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0
        bc=_sub3(c,b); w=(d4-d3)/((d4-d3)+(d5-d6))
        return _norm3(_sub3(p,_addscaled3(b,bc,w)))
    end
    denom=inv(va+vb+vc); v=vb*denom; w=vc*denom
    q=_addscaled3(_addscaled3(a,ab,v),ac,w)
    return _norm3(_sub3(p,q))
end

function _distance_triangle(p, a, b, c)
    d = _distance_triangle_kernel(p,a,b,c)
    isfinite(d) && return min(Float64(d), GMSH_MAX_SIZE)
    return setprecision(BigFloat, 256) do
        pb=map(BigFloat,p); ab=map(BigFloat,a); bb=map(BigFloat,b); cb=map(BigFloat,c)
        min(Float64(_distance_triangle_kernel(pb,ab,bb,cb)), GMSH_MAX_SIZE)
    end
end

@inline function _lower_interval_gap(x,a,b)
    d=x<a ? a-x : (x>b ? x-b : 0.0)
    return (isfinite(d)&&d>0) ? prevfloat(d) : d
end

@inline function _distance_box(p,lo,hi)
    # Down-round each positive coordinate gap, then use the L∞ distance. This is
    # a certified lower bound on the Euclidean AABB distance, so pruning cannot
    # discard a primitive that could improve `best`.
    return max(_lower_interval_gap(p[1],lo[1],hi[1]),
               _lower_interval_gap(p[2],lo[2],hi[2]),
               _lower_interval_gap(p[3],lo[3],hi[3]))
end

@inline function _distance_primitive(field::DistanceField,id::Int,p)
    np=size(field.points,2);ns=size(field.segments,2)
    if id<=np
        q=(field.points[1,id],field.points[2,id],field.points[3,id])
        return min(_norm3(_sub3(p,q)),GMSH_MAX_SIZE)
    elseif id<=np+ns
        j=id-np
        a=(field.segments[1,j],field.segments[2,j],field.segments[3,j])
        b=(field.segments[4,j],field.segments[5,j],field.segments[6,j])
        return _distance_segment(p,a,b)
    else
        j=id-np-ns
        a=(field.triangles[1,j],field.triangles[2,j],field.triangles[3,j])
        b=(field.triangles[4,j],field.triangles[5,j],field.triangles[6,j])
        c=(field.triangles[7,j],field.triangles[8,j],field.triangles[9,j])
        return _distance_triangle(p,a,b,c)
    end
end

function _distance_bvh(field::DistanceField,node::Int,p,best::Float64)
    _distance_box(p,field.bvh_lo[node],field.bvh_hi[node])>=best && return best
    count=field.bvh_count[node]
    if count>0
        first=field.bvh_first[node]
        @inbounds for pos in first:first+count-1
            best=min(best,_distance_primitive(field,field.bvh_order[pos],p))
        end
        return best
    end
    left=field.bvh_left[node];right=field.bvh_right[node]
    dl=_distance_box(p,field.bvh_lo[left],field.bvh_hi[left])
    dr=_distance_box(p,field.bvh_lo[right],field.bvh_hi[right])
    if dl<=dr
        dl<best && (best=_distance_bvh(field,left,p,best))
        dr<best && (best=_distance_bvh(field,right,p,best))
    else
        dr<best && (best=_distance_bvh(field,right,p,best))
        dl<best && (best=_distance_bvh(field,left,p,best))
    end
    return best
end

# Return both the nearest primitive distance and its deterministic primitive id.
# This is used by catalog fields that associate data with point samples. Equal
# distances select the lowest original id, independently of BVH traversal order.
function _nearest_distance_bvh(field::DistanceField,node::Int,p,best::Float64,best_id::Int)
    _distance_box(p,field.bvh_lo[node],field.bvh_hi[node])>best && return best,best_id
    count=field.bvh_count[node]
    if count>0
        first=field.bvh_first[node]
        @inbounds for pos in first:first+count-1
            id=field.bvh_order[pos]
            # `_nearest_point` is only defined for point-only fields.  Keep the
            # uncapped distance here: capping every candidate at
            # `GMSH_MAX_SIZE` destroys nearest-sample identity for far queries.
            q=(field.points[1,id],field.points[2,id],field.points[3,id])
            d=_norm3(_sub3(p,q))
            if d<best || (d==best && id<best_id)
                best=d; best_id=id
            end
        end
        return best,best_id
    end
    left=field.bvh_left[node]; right=field.bvh_right[node]
    dl=_distance_box(p,field.bvh_lo[left],field.bvh_hi[left])
    dr=_distance_box(p,field.bvh_lo[right],field.bvh_hi[right])
    if dl<=dr
        dl<=best && ((best,best_id)=_nearest_distance_bvh(field,left,p,best,best_id))
        dr<=best && ((best,best_id)=_nearest_distance_bvh(field,right,p,best,best_id))
    else
        dr<=best && ((best,best_id)=_nearest_distance_bvh(field,right,p,best,best_id))
        dl<=best && ((best,best_id)=_nearest_distance_bvh(field,left,p,best,best_id))
    end
    return best,best_id
end

function _nearest_point_big(field::DistanceField,p)::Tuple{Float64,Int}
    # A Float64 coordinate subtraction can overflow even though both operands
    # are finite.  This slow path is reached only when every candidate distance
    # overflowed; BigFloat then determines the sample identity without changing
    # the public Float64 distance contract.
    best=nothing
    best_id=typemax(Int)
    setprecision(BigFloat,256) do
        px=BigFloat(p[1]); py=BigFloat(p[2]); pz=BigFloat(p[3])
        @inbounds for id in axes(field.points,2)
            dx=px-BigFloat(field.points[1,id])
            dy=py-BigFloat(field.points[2,id])
            dz=pz-BigFloat(field.points[3,id])
            d2=dx*dx+dy*dy+dz*dz
            if best===nothing || d2<best || (d2==best && id<best_id)
                best=d2
                best_id=id
            end
        end
    end
    best_id==typemax(Int) && throw(ErrorException(
        "nearest-point query failed to select a sample"))
    return Inf,(best_id::Int)
end

@inline function _nearest_point(field::DistanceField,p)::Tuple{Float64,Int}
    isempty(field.segments) && isempty(field.triangles) || throw(ArgumentError(
        "nearest-point query requires a point-only DistanceField"))
    # Start at Inf so a finite query farther than GMSH_MAX_SIZE still visits a
    # leaf and returns a valid deterministic sample id.
    best,best_id=_nearest_distance_bvh(field,1,p,Inf,typemax(Int))
    isfinite(best) && return best,best_id
    return _nearest_point_big(field,p)
end

function field_value(field::DistanceField, x, y, z)
    p = (_float_value(x,"DistanceField","x"), _float_value(y,"DistanceField","y"),
         _float_value(z,"DistanceField","z"))
    return _distance_bvh(field,1,p,GMSH_MAX_SIZE)
end

"""
    ThresholdField(input; dist_min=1, dist_max=10, size_min=0.1,
                   size_max=1, sigmoid=false, stop_at_dist_max=false)

Gmsh `Threshold` semantics: `size_min` below `dist_min`, `size_max` above
`dist_max`, and linear (or Gmsh's logistic sigmoid) interpolation in between.
With `stop_at_dist_max=true`, values at and beyond `dist_max` return
[`GMSH_MAX_SIZE`](@ref).
"""
struct ThresholdField{F<:AbstractField} <: AbstractSizeField
    input::F
    dist_min::Float64
    dist_max::Float64
    size_min::Float64
    size_max::Float64
    sigmoid::Bool
    stop_at_dist_max::Bool
end

function ThresholdField(input::F; dist_min::Real=1.0, dist_max::Real=10.0,
                        size_min::Real=0.1, size_max::Real=1.0,
                        sigmoid::Bool=false, stop_at_dist_max::Bool=false) where {F<:AbstractField}
    dmin=_float_value(dist_min,"ThresholdField","dist_min")
    dmax=_float_value(dist_max,"ThresholdField","dist_max")
    dmax>dmin || throw(ArgumentError("ThresholdField: require dist_max > dist_min"))
    hmin=_float_value(size_min,"ThresholdField","size_min")
    hmax=_float_value(size_max,"ThresholdField","size_max")
    return ThresholdField{F}(input,dmin,dmax,hmin,hmax,sigmoid,stop_at_dist_max)
end

@inline function _threshold_value(field::ThresholdField, x, y, z, entity)
    d=_checked_field_result(field_value(field.input,x,y,z,entity),"ThresholdField input",x,y,z)
    field.stop_at_dist_max && d>=field.dist_max && return GMSH_MAX_SIZE
    r=clamp((d-field.dist_min)/(field.dist_max-field.dist_min),0.0,1.0)
    if field.sigmoid
        e=exp(12r-6); r=e/(1+e)
    end
    return muladd(field.size_max-field.size_min,r,field.size_min)
end
field_value(field::ThresholdField,x,y,z)=_threshold_value(field,x,y,z,nothing)
field_value(field::ThresholdField,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_threshold_value(field,x,y,z,entity)

"""
    BoxField(xmin, xmax, ymin, ymax, zmin, zmax; vin, vout, thickness=0)

Gmsh `Box` field.  Returns `vin` inside the closed axis-aligned box and
`vout` outside.  A positive `thickness` linearly interpolates in the exterior
layer using Euclidean distance to the box. Reversed bounds make the interior
empty, while their endpoints still define the transition layer, matching
Gmsh 4.15.2.
"""
struct BoxField <: AbstractSizeField
    xmin::Float64; xmax::Float64
    ymin::Float64; ymax::Float64
    zmin::Float64; zmax::Float64
    vin::Float64; vout::Float64; thickness::Float64
end

function BoxField(xmin::Real,xmax::Real,ymin::Real,ymax::Real,zmin::Real,zmax::Real;
                  vin::Real, vout::Real=GMSH_MAX_SIZE, thickness::Real=0.0)
    xa=_float_value(xmin,"BoxField","xmin"); xb=_float_value(xmax,"BoxField","xmax")
    ya=_float_value(ymin,"BoxField","ymin"); yb=_float_value(ymax,"BoxField","ymax")
    za=_float_value(zmin,"BoxField","zmin"); zb=_float_value(zmax,"BoxField","zmax")
    vi=_float_value(vin,"BoxField","vin"); vo=_float_value(vout,"BoxField","vout")
    th=_float_value(thickness,"BoxField","thickness")
    return BoxField(xa,xb,ya,yb,za,zb,vi,vo,th)
end

function field_value(field::BoxField, x, y, z)
    x=_float_value(x,"BoxField","x"); y=_float_value(y,"BoxField","y")
    z=_float_value(z,"BoxField","z")
    inside=field.xmin<=x<=field.xmax && field.ymin<=y<=field.ymax &&
           field.zmin<=z<=field.zmax
    inside && return field.vin
    if field.thickness>0
        # BoxField::computeDistance projects onto the three endpoint intervals
        # even when the corresponding Min/Max options are reversed.
        dx=max(min(field.xmin,field.xmax)-x,0.0,
               x-max(field.xmin,field.xmax))
        dy=max(min(field.ymin,field.ymax)-y,0.0,
               y-max(field.ymin,field.ymax))
        dz=max(min(field.zmin,field.zmax)-z,0.0,
               z-max(field.zmin,field.zmax))
        d=hypot(dx,dy,dz)
        d<=field.thickness &&
            return muladd(field.vout-field.vin,d/field.thickness,field.vin)
    end
    return field.vout
end

"""
    BallField(center, radius; vin, vout=GMSH_MAX_SIZE, thickness=0)

Gmsh `Ball` field. Returns `vin` strictly inside the sphere and `vout` outside.
For positive `thickness`, the exterior shell from `radius` through
`radius + thickness` is interpolated linearly from `vin` to `vout`. A
non-positive radius has no interior; a non-positive thickness disables the
transition, as in Gmsh 4.15.2.
"""
struct BallField <: AbstractSizeField
    center::NTuple{3,Float64}
    radius::Float64
    vin::Float64
    vout::Float64
    thickness::Float64
end

function BallField(center,radius::Real;vin::Real,vout::Real=GMSH_MAX_SIZE,
                   thickness::Real=0.0)
    c=_point3(center,"BallField")
    r=_float_value(radius,"BallField","radius")
    vi=_float_value(vin,"BallField","vin")
    vo=_float_value(vout,"BallField","vout")
    th=_float_value(thickness,"BallField","thickness")
    return BallField(c,r,vi,vo,th)
end

function field_value(field::BallField,x,y,z)
    p=(_float_value(x,"BallField","x"),_float_value(y,"BallField","y"),
       _float_value(z,"BallField","z"))
    d=hypot(p[1]-field.center[1],p[2]-field.center[2],p[3]-field.center[3])
    d<field.radius && return field.vin
    if field.thickness>0
        outside=d-field.radius
        outside<=field.thickness &&
            return muladd(field.vout-field.vin,outside/field.thickness,field.vin)
    end
    return field.vout
end

"""
    CylinderField(center, axis, radius; vin, vout=GMSH_MAX_SIZE)

Gmsh finite `Cylinder` field. `center` is the cylinder midpoint and `axis` is
the vector from the midpoint to either cap, so the cylinder spans
`center - axis` through `center + axis`. Both the radial and cap tests are
strict. A negative radius is equivalent to its magnitude, and a zero axis
produces an empty cylinder, matching Gmsh 4.15.2.
"""
struct CylinderField <: AbstractSizeField
    center::NTuple{3,Float64}
    unit_axis::NTuple{3,Float64}
    half_length::Float64
    radius::Float64
    vin::Float64
    vout::Float64
end

function CylinderField(center,axis,radius::Real;vin::Real,
                       vout::Real=GMSH_MAX_SIZE)
    c=_point3(center,"CylinderField center")
    a=_point3(axis,"CylinderField axis")
    half=hypot(a...)
    isfinite(half) || throw(ArgumentError(
        "CylinderField: axis must have finite length"))
    r=_float_value(radius,"CylinderField","radius")
    vi=_float_value(vin,"CylinderField","vin")
    vo=_float_value(vout,"CylinderField","vout")
    unit_axis=half>0 ? ntuple(i -> a[i]/half,3) : (0.0,0.0,0.0)
    return CylinderField(c,unit_axis,half,abs(r),vi,vo)
end

function _inside_cylinder(field::CylinderField,p)
    d=ntuple(i -> p[i]-field.center[i],3)
    axial=d[1]*field.unit_axis[1]+d[2]*field.unit_axis[2]+d[3]*field.unit_axis[3]
    radial=ntuple(i -> d[i]-axial*field.unit_axis[i],3)
    r=hypot(radial...)
    if isfinite(axial)&&isfinite(r)
        return abs(axial)<field.half_length && r<field.radius
    end
    return setprecision(BigFloat,256) do
        db=ntuple(i -> BigFloat(p[i])-BigFloat(field.center[i]),3)
        ub=map(BigFloat,field.unit_axis)
        sb=db[1]*ub[1]+db[2]*ub[2]+db[3]*ub[3]
        rb=ntuple(i -> db[i]-sb*ub[i],3)
        abs(sb)<BigFloat(field.half_length) &&
            sqrt(rb[1]^2+rb[2]^2+rb[3]^2)<BigFloat(field.radius)
    end
end

function field_value(field::CylinderField,x,y,z)
    p=(_float_value(x,"CylinderField","x"),
       _float_value(y,"CylinderField","y"),
       _float_value(z,"CylinderField","z"))
    return _inside_cylinder(field,p) ? field.vin : field.vout
end

"""
    FrustumField(p1, p2; inner_r1=0, outer_r1=1, inner_r2=0, outer_r2=1,
                 inner_v1=0.1, outer_v1=1, inner_v2=0.1, outer_v2=1)

Gmsh `Frustum` field. Within the algebraic annular frustum from `p1` to `p2`,
interpolate bilinearly between the mesh sizes on the inner and outer radii at
both endpoints. Radius order is significant and may be reversed. Outside the
axial interval or radial annulus, return [`GMSH_MAX_SIZE`](@ref); coincident
endpoints and an equal inner/outer radius produce no active region.
"""
struct FrustumField <: AbstractSizeField
    p1::NTuple{3,Float64}
    unit_axis::NTuple{3,Float64}
    length::Float64
    inner_r1::Float64
    outer_r1::Float64
    inner_r2::Float64
    outer_r2::Float64
    inner_v1::Float64
    outer_v1::Float64
    inner_v2::Float64
    outer_v2::Float64
end

function FrustumField(p1,p2;inner_r1::Real=0.0,outer_r1::Real=1.0,
                      inner_r2::Real=0.0,outer_r2::Real=1.0,
                      inner_v1::Real=0.1,outer_v1::Real=1.0,
                      inner_v2::Real=0.1,outer_v2::Real=1.0)
    a=_point3(p1,"FrustumField p1");b=_point3(p2,"FrustumField p2")
    axis=ntuple(i -> b[i]-a[i],3); len=hypot(axis...)
    isfinite(len) ||
        throw(ArgumentError("FrustumField: endpoints must have finite separation"))
    r1i=_float_value(inner_r1,"FrustumField","inner_r1")
    r1o=_float_value(outer_r1,"FrustumField","outer_r1")
    r2i=_float_value(inner_r2,"FrustumField","inner_r2")
    r2o=_float_value(outer_r2,"FrustumField","outer_r2")
    v1i=_float_value(inner_v1,"FrustumField","inner_v1")
    v1o=_float_value(outer_v1,"FrustumField","outer_v1")
    v2i=_float_value(inner_v2,"FrustumField","inner_v2")
    v2o=_float_value(outer_v2,"FrustumField","outer_v2")
    unit_axis=len>0 ? ntuple(i -> axis[i]/len,3) : (0.0,0.0,0.0)
    return FrustumField(a,unit_axis,len,r1i,r1o,r2i,r2o,
                        v1i,v1o,v2i,v2o)
end

function field_value(field::FrustumField,x,y,z)
    p=(_float_value(x,"FrustumField","x"),
       _float_value(y,"FrustumField","y"),
       _float_value(z,"FrustumField","z"))
    d=ntuple(i -> p[i]-field.p1[i],3)
    axial=d[1]*field.unit_axis[1]+d[2]*field.unit_axis[2]+d[3]*field.unit_axis[3]
    isfinite(axial) || return GMSH_MAX_SIZE
    u=axial/field.length
    0<=u<=1 || return GMSH_MAX_SIZE
    radialv=ntuple(i -> d[i]-axial*field.unit_axis[i],3)
    radial=hypot(radialv...)
    isfinite(radial) || return GMSH_MAX_SIZE
    ri=muladd(u,field.inner_r2-field.inner_r1,field.inner_r1)
    ro=muladd(u,field.outer_r2-field.outer_r1,field.outer_r1)
    v=(radial-ri)/(ro-ri)
    0<=v<=1 || return GMSH_MAX_SIZE
    vi=muladd(u,field.inner_v2-field.inner_v1,field.inner_v1)
    vo=muladd(u,field.outer_v2-field.outer_v1,field.outer_v1)
    return muladd(v,vo-vi,vi)
end

# ── Composition and final size contract ───────────────────────────────────────

"""Pointwise minimum of one or more scalar fields (Gmsh `Min`)."""
struct MinSize{F<:Tuple} <: AbstractSizeField
    fields::F
end

"""Gmsh `Max`; anisotropic inputs contribute their largest directional size."""
struct MaxSize{F<:Tuple} <: AbstractSizeField
    fields::F
end

function MinSize(fields::Union{Tuple,AbstractVector})
    values=Tuple(fields)
    all(child -> child isa AbstractField,values) ||
        throw(ArgumentError("MinSize: every input must be an AbstractField"))
    return MinSize{typeof(values)}(values)
end


function MaxSize(fields::Union{Tuple,AbstractVector})
    values=Tuple(fields)
    all(child -> child isa AbstractField,values) ||
        throw(ArgumentError("MaxSize: every input must be an AbstractField"))
    return MaxSize{typeof(values)}(values)
end

const MinField = MinSize
const MaxField = MaxSize

@inline _min_value(::Tuple{},x,y,z)=GMSH_MAX_SIZE
@inline function _min_child_value(child::AbstractField,x,y,z,entity)
    if child isa AbstractAnisoField
        eigenvalues=metric_eigenvalues(metric_at(child,x,y,z,entity))
        value=inv(sqrt(maximum(eigenvalues)))
        return _checked_field_result(value,"MinSize anisotropic input",x,y,z)
    end
    return _checked_field_result(field_value(child,x,y,z,entity),"MinSize input",x,y,z)
end
@inline function _min_value(fields::Tuple,x,y,z)
    child=_min_child_value(first(fields),x,y,z,nothing)
    return min(child,_min_value(Base.tail(fields),x,y,z))
end
@inline field_value(field::MinSize,x,y,z)=_min_value(field.fields,x,y,z)

@inline _min_value(::Tuple{},x,y,z,entity)=GMSH_MAX_SIZE
@inline function _min_value(fields::Tuple,x,y,z,entity)
    child=_min_child_value(first(fields),x,y,z,entity)
    return min(child,_min_value(Base.tail(fields),x,y,z,entity))
end
@inline field_value(field::MinSize,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_min_value(field.fields,x,y,z,entity)

@inline _max_value(::Tuple{},x,y,z)=-GMSH_MAX_SIZE
@inline function _max_child_value(child::AbstractField,x,y,z,entity)
    if child isa AbstractAnisoField
        eigenvalues=metric_eigenvalues(metric_at(child,x,y,z,entity))
        value=inv(sqrt(minimum(eigenvalues)))
        return _checked_field_result(value,"MaxSize anisotropic input",x,y,z)
    end
    return _checked_field_result(field_value(child,x,y,z,entity),"MaxSize input",x,y,z)
end
@inline function _max_value(fields::Tuple,x,y,z)
    child=_max_child_value(first(fields),x,y,z,nothing)
    return max(child,_max_value(Base.tail(fields),x,y,z))
end
@inline field_value(field::MaxSize,x,y,z)=_max_value(field.fields,x,y,z)

@inline _max_value(::Tuple{},x,y,z,entity)=-GMSH_MAX_SIZE
@inline function _max_value(fields::Tuple,x,y,z,entity)
    child=_max_child_value(first(fields),x,y,z,entity)
    return max(child,_max_value(Base.tail(fields),x,y,z,entity))
end
@inline field_value(field::MaxSize,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_max_value(field.fields,x,y,z,entity)

"""
    BoundedSize(field; size_min=0, size_max=GMSH_MAX_SIZE, factor=1)

Apply Gmsh's final size policy: clamp the input field into
`[size_min,size_max]`, then multiply by `factor`.  This also turns a generic
scalar field such as `DistanceField` into a valid mesh-size field when
`size_min > 0`. With a zero lower bound, [`size_at`](@ref) still rejects any
non-positive value produced by a generic scalar input at the queried point.
"""
struct BoundedSize{F<:AbstractField} <: AbstractSizeField
    input::F
    size_min::Float64
    size_max::Float64
    factor::Float64
end

function BoundedSize(input::F; size_min::Real=0.0, size_max::Real=GMSH_MAX_SIZE,
                     factor::Real=1.0) where {F<:AbstractField}
    hmin=_float_value(size_min,"BoundedSize","size_min")
    hmin>=0 || throw(ArgumentError("BoundedSize: size_min must be non-negative"))
    hmax=_positive_value(size_max,"BoundedSize","size_max")
    hmax>=hmin || throw(ArgumentError("BoundedSize: require size_max >= size_min"))
    f=_positive_value(factor,"BoundedSize","factor")
    isfinite(hmax*f) || throw(ArgumentError("BoundedSize: size_max * factor must be finite"))
    input isa AbstractAnisoField && return BoundedAniso{F}(input,hmin,hmax,f)
    return BoundedSize{F}(input,hmin,hmax,f)
end

include("SizeFieldCatalog.jl")

"""Directional counterpart returned by `BoundedSize` for anisotropic inputs."""
struct BoundedAniso{F<:AbstractAnisoField} <: AbstractAnisoField
    input::F
    size_min::Float64
    size_max::Float64
    factor::Float64
end

function _bounded_metric(field::BoundedAniso,x,y,z,entity)
    input=metric_at(field.input,x,y,z,entity)
    input isa Metric3 || throw(ArgumentError(
        "BoundedSize: metric_at returned $(typeof(input)), not Metric3"))
    values,vectors=_sym3_eigh(input.m11,input.m22,input.m33,
                              input.m12,input.m13,input.m23)
    all(>(0),values) || throw(ArgumentError(
        "BoundedSize: input metric is not positive definite at ($x,$y,$z)"))
    bounded=ntuple(3) do i
        h=inv(sqrt(values[i]))
        hb=field.factor*clamp(h,field.size_min,field.size_max)
        invh=inv(hb); value=invh*invh
        (isfinite(value) && value>0) || throw(ArgumentError(
            "BoundedSize: bounded anisotropic metric is not representable at ($x,$y,$z)"))
        value
    end
    return metric_from_axes(bounded[1],bounded[2],bounded[3],vectors...)
end
metric_at(field::BoundedAniso,x,y,z)=_bounded_metric(field,x,y,z,nothing)
metric_at(field::BoundedAniso,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_bounded_metric(field,x,y,z,entity)
field_value(field::BoundedAniso,x,y,z)=metric_size(metric_at(field,x,y,z))
field_value(field::BoundedAniso,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=metric_size(metric_at(field,x,y,z,entity))

function field_value(field::BoundedSize,x,y,z)
    value=_checked_field_result(field_value(field.input,x,y,z),"BoundedSize input",x,y,z)
    return field.factor*clamp(value,field.size_min,field.size_max)
end
function field_value(field::BoundedSize,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}
    value=_checked_field_result(field_value(field.input,x,y,z,entity),
                                "BoundedSize input",x,y,z)
    return field.factor*clamp(value,field.size_min,field.size_max)
end

# Resource ownership follows the field graph: closing a composite closes every
# external-process leaf exactly once, even when the same leaf is referenced by
# several MathEval expressions or branches.
function _collect_external_fields!(out::Vector{ExternalProcessField},value)
    if value isa ExternalProcessField
        any(item -> item===value,out) || push!(out,value)
    elseif value isa AbstractField
        @inbounds for i in 1:fieldcount(typeof(value))
            _collect_external_fields!(out,getfield(value,i))
        end
    elseif value isa Tuple
        @inbounds for item in value
            _collect_external_fields!(out,item)
        end
    elseif value isa _MEFieldBinding
        _collect_external_fields!(out,value.field)
    end
    return out
end

function Base.close(field::AbstractField)
    resources=_collect_external_fields!(ExternalProcessField[],field)
    first_error=nothing
    for resource in resources
        try
            close(resource)
        catch err
            err isa InterruptException && rethrow()
            first_error===nothing && (first_error=err)
        end
    end
    first_error===nothing || throw(first_error)
    return nothing
end

function Base.isopen(field::AbstractField)
    resources=_collect_external_fields!(ExternalProcessField[],field)
    return all(isopen,resources)
end

# ── Strict construction from parsed .geo field declarations ──────────────────

function _geo_known_options(spec::GeoFieldSpec, allowed)
    unknown=sort!(collect(setdiff(keys(spec.options),allowed)))
    isempty(unknown) || throw(ArgumentError(
        "build_geo_size_field: unsupported option(s) for Field[$(spec.tag)] $(spec.kind): " *
        join(unknown,", ")))
    return nothing
end

function _geo_float(spec::GeoFieldSpec,name::String,default)
    haskey(spec.options,name) || return default
    raw=strip(spec.options[name]); value=tryparse(Float64,raw)
    (value!==nothing && isfinite(value)) || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be a finite numeric literal (got $raw)"))
    return value::Float64
end

function _geo_float_alias(spec::GeoFieldSpec,name::String,deprecated::String,default)
    option=_geo_last_alias(spec,(name,deprecated))
    return option===nothing ? default : _geo_float(spec,option,default)
end

function _geo_int_alias(spec::GeoFieldSpec,name::String,deprecated::String,default::Int)
    option=_geo_last_alias(spec,(name,deprecated))
    return option===nothing ? default : _geo_int(spec,option;default=default)
end

function _geo_int(spec::GeoFieldSpec,name::String; required::Bool=false,default::Int=0)
    if !haskey(spec.options,name)
        required && throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] is missing $name"))
        return default
    end
    raw=strip(spec.options[name]); numeric=tryparse(Float64,raw)
    (numeric!==nothing && isfinite(numeric)) || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be a finite numeric literal (got $raw)"))
    return try
        trunc(Int,numeric::Float64)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError) || rethrow()
        throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)].$name is outside the platform Int range (got $raw)"))
    end
end

function _geo_bool(spec::GeoFieldSpec,name::String,default::Bool=false)
    haskey(spec.options,name) || return default
    raw=lowercase(strip(spec.options[name]))
    raw=="true" && return true
    raw=="false" && return false
    numeric=tryparse(Float64,raw)
    (numeric!==nothing && isfinite(numeric)) || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be a finite numeric or boolean literal (got $raw)"))
    return !iszero(numeric::Float64)
end

function _geo_list(spec::GeoFieldSpec,name::String; required::Bool=false)
    if !haskey(spec.options,name)
        required && throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] is missing $name"))
        return String[]
    end
    caller="build_geo_size_field: Field[$(spec.tag)].$name"
    return _geo_split_list(spec.options[name],caller)
end

function _geo_infield(spec::GeoFieldSpec;default::Int=1)
    option=_geo_last_alias(spec,("InField","IField"))
    return option===nothing ? default : _geo_int(spec,option;required=true)
end
function _geo_float_string(spec::GeoFieldSpec,name::String,deprecated::String,default::String)
    option=_geo_last_alias(spec,(name,deprecated))
    option===nothing && return default
    return _geo_string_value(spec.options[option],
        "build_geo_size_field: Field[$(spec.tag)].$option")
end
function _geo_string_value(value::AbstractString,caller::AbstractString)
    raw=String(strip(value))
    isempty(raw) && throw(ArgumentError("$caller must be non-empty"))
    leading=Base.first(raw);trailing=Base.last(raw)
    leading_quoted=leading=='"' || leading=='\''
    trailing_quoted=trailing=='"' || trailing=='\''
    if leading_quoted || trailing_quoted
        (leading_quoted && trailing_quoted && leading==trailing) ||
            throw(ArgumentError("$caller has mismatched string quotes"))
    else
        return raw
    end
    first_index=nextind(raw,firstindex(raw));last_index=prevind(raw,lastindex(raw))
    first_index>last_index && return ""
    # Gmsh preserves backslashes and other special characters literally inside
    # both quote styles; these are not Julia/C escape sequences.
    return raw[first_index:last_index]
end
function _geo_int_list_opt(spec::GeoFieldSpec,name::String,deprecated::String="")
    names=isempty(deprecated) ? (name,) : (name,deprecated)
    option=_geo_last_alias(spec,names)
    raws=option===nothing ? String[] : _geo_list(spec,option)
    tags=Int[]
    for raw in raws
        tag=tryparse(Int,raw)
        (tag!==nothing && tag>0) || throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)].$(something(option,name)) needs positive integer tags (got $raw)"))
        push!(tags,tag)
    end
    return tags
end

function _geo_float_list_opt(spec::GeoFieldSpec,name::String,deprecated::String="")
    names=isempty(deprecated) ? (name,) : (name,deprecated)
    option=_geo_last_alias(spec,names)
    option===nothing && return Float64[]
    values=Float64[]
    for raw in _geo_list(spec,option)
        value=tryparse(Float64,raw)
        value!==nothing && isfinite(value) || throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)].$option needs finite numeric literals (got $raw)"))
        push!(values,value)
    end
    return values
end
function _geo_field_tags(spec::GeoFieldSpec,name::String;
                         required::Bool=true,nonempty::Bool=true)
    raws=_geo_list(spec,name;required=required); tags=Int[]
    for raw in raws
        tag=tryparse(Int,raw)
        (tag!==nothing && tag>0) || throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)].$name needs positive field tags (got $raw)"))
        push!(tags,tag)
    end
    nonempty && isempty(tags) && throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name cannot be empty"))
    return tags
end

function _geo_entity_name(raw::AbstractString)
    name=String(strip(raw))
    (startswith(name,"+") || startswith(name,"-")) && (name=String(strip(name[2:end])))
    endswith(name,"[]") && (name=String(strip(name[1:end-2])))
    isempty(name) && throw(ArgumentError("build_geo_size_field: empty geometric entity reference"))
    occursin(r"[{}():]",name) && throw(ArgumentError(
        "build_geo_size_field: unsupported geometric entity expression $raw"))
    return name
end

function _geo_resolve(entities,dim::Int,raw::AbstractString)
    name=_geo_entity_name(raw)
    value = if entities isa Function
        entities(dim,name)
    else
        found=false; result=nothing
        candidates=Any[(dim,name),name]
        numeric=tryparse(Int,name)
        numeric===nothing || append!(candidates,Any[(dim,numeric),numeric])
        for key in candidates
            if haskey(entities,key)
                found=true; result=entities[key]; break
            end
        end
        found ? result : nothing
    end
    value===nothing && throw(ArgumentError(
        "build_geo_size_field: no dimension-$dim geometry supplied for entity $raw"))
    meshes = value isa Mesh ? Mesh[value] :
             value isa AbstractVector && all(v -> v isa Mesh,value) ? Mesh[value...] :
             throw(ArgumentError(
                 "build_geo_size_field: entity $raw must resolve to a Mesh or vector of Mesh values"))
    isempty(meshes) && throw(ArgumentError(
        "build_geo_size_field: entity $raw resolved to an empty mesh collection"))
    return meshes
end

@inline _mesh_point(mesh::Mesh,i) =
    (mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])

function _geo_last_alias(spec::GeoFieldSpec,names::Tuple)
    for name in Iterators.reverse(spec.option_order)
        name in names && haskey(spec.options,name) && return name
    end
    # Defensive fallback for externally constructed four-argument specs whose
    # order vector omitted a key.
    for name in names
        haskey(spec.options,name) && return name
    end
    return nothing
end

function _geo_curve_sample_points(mesh::Mesh,sampling::Int,raw,
                                  max_samples::Int)
    ns=size(mesh.segs,2)
    ns>0 || throw(ArgumentError(
        "build_geo_size_field: curve entity $raw contains no segments"))
    nv=size(mesh.coords,2);adj=[Int[] for _ in 1:nv]
    @inbounds for s in axes(mesh.segs,2)
        a=Int(mesh.segs[1,s]);b=Int(mesh.segs[2,s])
        push!(adj[a],s);push!(adj[b],s)
    end
    used_vertices=Int[i for i in 1:nv if !isempty(adj[i])]
    endpoints=Int[]
    for i in used_vertices
        degree=length(adj[i])
        degree in (1,2) || throw(ArgumentError(
            "build_geo_size_field: curve entity $raw is branching at node $i"))
        degree==1 && push!(endpoints,i)
    end
    length(endpoints) in (0,2) || throw(ArgumentError(
        "build_geo_size_field: curve entity $raw is not one connected open/closed polyline"))
    start=isempty(endpoints) ? minimum(used_vertices) : minimum(endpoints)
    used=falses(ns);order=Int[start];current=start
    for _ in 1:ns
        edge=0
        for candidate in adj[current]
            if !used[candidate] && (edge==0 || candidate<edge)
                edge=candidate
            end
        end
        edge!=0 || throw(ArgumentError(
            "build_geo_size_field: curve entity $raw has disconnected segments"))
        used[edge]=true
        a=Int(mesh.segs[1,edge]);b=Int(mesh.segs[2,edge])
        current=a==current ? b : a
        push!(order,current)
    end
    all(used) || throw(ArgumentError(
        "build_geo_size_field: curve entity $raw has disconnected segments"))
    if isempty(endpoints)
        current==start || throw(ArgumentError(
            "build_geo_size_field: closed curve entity $raw does not form one cycle"))
    else
        current in endpoints || throw(ArgumentError(
            "build_geo_size_field: curve entity $raw has inconsistent segment topology"))
    end
    lengths=Vector{Float64}(undef,ns);total=0.0
    @inbounds for i in 1:ns
        a=_mesh_point(mesh,order[i]);b=_mesh_point(mesh,order[i+1])
        lengths[i]=hypot(b[1]-a[1],b[2]-a[2],b[3]-a[3])
        (isfinite(lengths[i])&&lengths[i]>0) || throw(ArgumentError(
            "build_geo_size_field: curve entity $raw has a zero or non-finite segment"))
        total+=lengths[i]
        isfinite(total) || throw(ArgumentError(
            "build_geo_size_field: curve entity $raw length overflows Float64"))
    end
    count=max(0,sampling-2)
    count<=max_samples || throw(ArgumentError(
        "build_geo_size_field: curve sampling exceeds max_distance_samples=$max_samples"))
    points=Vector{NTuple{3,Float64}}(undef,count);segment=1;cumulative=0.0
    @inbounds for i in 1:count
        target=total*(i/(sampling-1))
        while segment<ns && cumulative+lengths[segment]<target
            cumulative+=lengths[segment];segment+=1
        end
        a=_mesh_point(mesh,order[segment]);b=_mesh_point(mesh,order[segment+1])
        t=(target-cumulative)/lengths[segment]
        points[i]=(a[1]+t*(b[1]-a[1]),a[2]+t*(b[2]-a[2]),
                   a[3]+t*(b[3]-a[3]))
    end
    return points
end

function _geo_surface_sample_points(mesh::Mesh,sampling::Int,raw,
                                    max_samples::Int)
    size(mesh.tris,2)>0 || throw(ArgumentError(
        "build_geo_size_field: surface entity $raw contains no triangles"))
    xmin,xmax=extrema(@view mesh.coords[1,:]);ymin,ymax=extrema(@view mesh.coords[2,:])
    zmin,zmax=extrema(@view mesh.coords[3,:]);diagonal=hypot(xmax-xmin,ymax-ymin,zmax-zmin)
    (isfinite(diagonal)&&diagonal>0) || throw(ArgumentError(
        "build_geo_size_field: surface entity $raw has a zero or non-finite bounding box"))
    max_distance=diagonal/sampling
    (isfinite(max_distance)&&max_distance>0) || throw(ArgumentError(
        "build_geo_size_field: surface sampling scale is not representable"))
    points=NTuple{3,Float64}[]
    @inbounds for t in axes(mesh.tris,2)
        a=_mesh_point(mesh,mesh.tris[1,t]);b=_mesh_point(mesh,mesh.tris[2,t])
        c=_mesh_point(mesh,mesh.tris[3,t])
        longest=max(hypot(b[1]-a[1],b[2]-a[2],b[3]-a[3]),
                    hypot(c[1]-b[1],c[2]-b[2],c[3]-b[3]),
                    hypot(a[1]-c[1],a[2]-c[2],a[3]-c[3]))
        isfinite(longest) || throw(ArgumentError(
            "build_geo_size_field: surface entity $raw has a non-finite triangle edge"))
        ratio=longest/max_distance
        ratio<=prevfloat(Float64(typemax(Int))) || throw(ArgumentError(
            "build_geo_size_field: surface sampling count exceeds the platform Int limit"))
        n=max(floor(Int,ratio),1)
        added=try
            np1=Base.checked_add(n,1)
            iseven(n) ? Base.checked_mul(n÷2,np1) :
                        Base.checked_mul(n,np1÷2)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "build_geo_size_field: surface sampling count overflows the platform Int limit"))
        end
        added<=max_samples && length(points)<=max_samples-added || throw(ArgumentError(
            "build_geo_size_field: surface sampling exceeds max_distance_samples=$max_samples"))
        for iu in 0:n-1,iv in 0:n-iu-1
            u=iu/n;v=iv/n;w=1-u-v
            push!(points,(w*a[1]+u*b[1]+v*c[1],w*a[2]+u*b[2]+v*c[2],
                          w*a[3]+u*b[3]+v*c[3]))
        end
    end
    return points
end

function _geo_distance(spec::GeoFieldSpec,entities,max_samples::Int)
    allowed=Set(["PointsList","CurvesList","SurfacesList","Sampling",
                 "NodesList","EdgesList","FacesList","NNodesByEdge",
                 "NumPointsPerCurve","FieldX","FieldY","FieldZ"])
    _geo_known_options(spec,allowed)
    # Retained by Gmsh 4.15.2 as deprecated integer options, but explicitly
    # unused by DistanceField. Validate their syntax without creating graph
    # dependencies or changing the distance targets.
    for name in ("FieldX","FieldY","FieldZ")
        haskey(spec.options,name) && _geo_int(spec,name;required=true)
    end
    sampling_aliases=("Sampling","NNodesByEdge","NumPointsPerCurve")
    sampling_option=_geo_last_alias(spec,sampling_aliases)
    sampling=sampling_option===nothing ? 20 : _geo_int(spec,sampling_option)
    sampling>0 || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$(something(sampling_option,"Sampling")) must be positive"))
    points=NTuple{3,Float64}[]
    groups=((0,("PointsList","NodesList")),
            (1,("CurvesList","EdgesList")),
            (2,("SurfacesList","FacesList")))
    for (dim,names) in groups
        option=_geo_last_alias(spec,names)
        option===nothing && continue
        for raw in _geo_list(spec,option), mesh in _geo_resolve(entities,dim,raw)
            if dim==0
                length(points)<=max_samples-size(mesh.coords,2) || throw(ArgumentError(
                    "build_geo_size_field: point targets exceed max_distance_samples=$max_samples"))
                append!(points,(_mesh_point(mesh,i) for i in axes(mesh.coords,2)))
            elseif dim==1
                samples=_geo_curve_sample_points(mesh,sampling,raw,max_samples-length(points))
                append!(points,samples)
            else
                samples=_geo_surface_sample_points(mesh,sampling,raw,max_samples-length(points))
                append!(points,samples)
            end
        end
    end
    # Gmsh permits an empty Distance/Attractor target set. Its scalar operator
    # returns MAX_LC, which is also composable through MathEval before the final
    # background-size policy is applied.
    isempty(points) && return ConstantSize(GMSH_MAX_SIZE)
    return DistanceField(;points=points)
end

@inline function _geo_config_positive(value::Real,spec,name)
    value>0 || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be positive"))
    return Float64(value)
end

function _geo_extend_config(spec::GeoFieldSpec)
    curves=_geo_int_list_opt(spec,"CurvesList")
    surfaces=_geo_int_list_opt(spec,"SurfacesList")
    dist_max=_geo_float(spec,"DistMax",1.0)
    size_max=_geo_float(spec,"SizeMax",1.0)
    power=_geo_float(spec,"Power",1.0)
    return (;curves,surfaces,dist_max,size_max,power)
end

function _geo_postview_config(spec::GeoFieldSpec)
    view_index=_geo_int_alias(spec,"ViewIndex","IView",0)
    return (;view_index,view_tag=_geo_int(spec,"ViewTag";default=-1),
            crop_negative_values=_geo_bool(spec,"CropNegativeValues",true),
            use_closest=_geo_bool(spec,"UseClosest",true))
end

function _geo_attractor_aniso_config(spec::GeoFieldSpec)
    sampling_option=_geo_last_alias(spec,
        ("Sampling","NNodesByEdge","NumPointsPerCurve"))
    sampling=sampling_option===nothing ? 20 :
             _geo_int(spec,sampling_option;default=20)
    sampling>0 || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].Sampling must be positive"))
    dist_min=_geo_float_alias(spec,"DistMin","dMin",0.1)
    dist_max=_geo_float_alias(spec,"DistMax","dMax",0.5)
    dist_max>dist_min || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)] requires DistMax > DistMin"))
    size_min_tangent=_geo_config_positive(
        _geo_float_alias(spec,"SizeMinTangent","lMinTangent",0.5),spec,"SizeMinTangent")
    size_max_tangent=_geo_config_positive(
        _geo_float_alias(spec,"SizeMaxTangent","lMaxTangent",0.5),spec,"SizeMaxTangent")
    size_min_normal=_geo_config_positive(
        _geo_float_alias(spec,"SizeMinNormal","lMinNormal",0.05),spec,"SizeMinNormal")
    size_max_normal=_geo_config_positive(
        _geo_float_alias(spec,"SizeMaxNormal","lMaxNormal",0.5),spec,"SizeMaxNormal")
    return (;curves=_geo_int_list_opt(spec,"CurvesList","EdgesList"),sampling,
            dist_min,dist_max,size_min_tangent,size_max_tangent,
            size_min_normal,size_max_normal)
end

function _geo_boundary_layer_config(spec::GeoFieldSpec,
                                    mesh_boundary_layer_fan_elements::Int)
    size=_geo_float_alias(spec,"Size","hwall_n",0.1)
    ratio=_geo_float_alias(spec,"Ratio","ratio",1.1)
    size_far=_geo_float_alias(spec,"SizeFar","hfar",1.0)
    thickness=_geo_float_alias(spec,"Thickness","thickness",0.01)
    beta=_geo_config_positive(_geo_float(spec,"Beta",1.01),spec,"Beta")
    layers=_geo_int(spec,"NbLayers";default=10)
    layers>0 || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].NbLayers must be positive"))
    return (;curves=_geo_int_list_opt(spec,"CurvesList","EdgesList"),
            fan_points=_geo_int_list_opt(spec,"FanPointsList","FanNodesList"),
            fan_sizes=_geo_int_list_opt(spec,"FanPointsSizesList"),
            points=_geo_int_list_opt(spec,"PointsList","NodesList"),
            size,sizes=_geo_float_list_opt(spec,"SizesList","hwall_n_nodes"),
            ratio,size_far,thickness,quads=_geo_int(spec,"Quads";default=0),
            intersect_metrics=_geo_int(spec,"IntersectMetrics";default=0),
            aniso_max=_geo_float(spec,"AnisoMax",1e10),
            beta_law=_geo_int(spec,"BetaLaw";default=0),beta,nb_layers=layers,
            fan_elements=mesh_boundary_layer_fan_elements,
            excluded_surfaces=_geo_int_list_opt(spec,"ExcludedSurfacesList",
                                                 "ExcludedFaceList"))
end

function _geo_automatic_config(spec::GeoFieldSpec)
    filename=haskey(spec.options,"p4estFileToLoad") ?
        _geo_string_value(spec.options["p4estFileToLoad"],
                          "build_geo_size_field: Field[$(spec.tag)].p4estFileToLoad") : ""
    return (;p4est_file=filename,
            points_per_circle=_geo_int(spec,"nPointsPerCircle";
                default=spec.creation_mesh_size_from_curvature==0 ? 20 :
                        spec.creation_mesh_size_from_curvature),
            points_per_gap=_geo_int(spec,"nPointsPerGap";default=0),
            hmin=_geo_float(spec,"hMin",-1.0),hmax=_geo_float(spec,"hMax",-1.0),
            hbulk=_geo_float(spec,"hBulk",-1.0),
            gradation=_geo_float(spec,"gradation",1.1),
            smoothing=_geo_bool(spec,"smoothing",true),
            features=_geo_bool(spec,"features",true))
end

function _geo_context_field(context_fields,spec::GeoFieldSpec,entities,params,config)
    context_fields===nothing && throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)] $(spec.kind) requires a context_fields resolver"))
    value=if applicable(context_fields,spec,config,entities,params)
        context_fields(spec,config,entities,params)
    elseif applicable(context_fields,spec,entities,params)
        context_fields(spec,entities,params)
    elseif context_fields isa AbstractDict
        get(context_fields,spec.tag,nothing)
    else
        throw(ArgumentError(
            "build_geo_size_field: context_fields must be a dictionary or callable " *
            "(spec, config, entities, params)"))
    end
    value isa AbstractField || throw(ArgumentError(
        "build_geo_size_field: context resolver for Field[$(spec.tag)] $(spec.kind) " *
        "must return an AbstractField (got $(typeof(value)))"))
    return value
end

"""
    _build_geo_field(params::GeoParams, entities, root_field; kwargs...) -> AbstractField

Build the parsed `.geo` `Background Field` and apply `Mesh.MeshSizeMin`,
`Mesh.MeshSizeMax`, and `Mesh.MeshSizeFactor` in Gmsh order. `entities` is either
a dictionary or a callable resolving `(dimension, reference)` to a `Mesh` (or a
vector of meshes). Symbolic references such as `surface_group[]` are normalized to
`surface_group`; numeric entity tags are also accepted.

The strict builder covers the coordinate-only catalog directly and delegates the
five model/view-backed kinds through `context_fields`. Unsupported field kinds,
options, expressions, missing entities, cycles, and invalid references raise
`ArgumentError` instead of being ignored.

`base_dir` resolves relative `Structured.FileName` paths, `model_bbox` supplies
the six-coordinate model bounds needed by `Octree`, and the resource keywords
bound structured samples and octree cells. `model_characteristic_length`
overrides the model scale used for Gmsh's default finite-difference step.
Otherwise degenerate bounds are padded with Gmsh 4.15.2's
`FinishUpBoundingBox` rules; `geometry_tolerance` supplies the corresponding
geometry tolerance (with Gmsh's effective minimum of `1e-6`).
`entity_boundaries` and `entity_embedded` map `(dimension, tag)` entities to
their CAD topology for `Restrict` and `Constant` expansion.
`context_fields` is a tag-indexed dictionary override or a
`(spec, config, entities, params)` callable for fields that require model/view
state unavailable in `GeoParams`:
`Extend`, `PostView`, `AttractorAnisoCurve`, `BoundaryLayer`, and
`AutomaticMeshSizeField`. The resolver is responsible for interpreting the
typed, validated `config` and must return an `AbstractField`. The legacy
three-argument callable form remains accepted after the same option validation.
"""
function _gmsh_bbox_characteristic_length(bbox::NTuple{6,Float64},
                                          geometry_tolerance::Real)
    tolerance=_float_value(geometry_tolerance,"build_geo_size_field",
                           "geometry_tolerance")
    tolerance>=0 || throw(ArgumentError(
        "build_geo_size_field: geometry_tolerance must be non-negative"))
    tol=max(1e-6,tolerance)
    dx=bbox[2]-bbox[1];dy=bbox[4]-bbox[3];dz=bbox[6]-bbox[5]
    all(isfinite,(dx,dy,dz)) || throw(ArgumentError(
        "build_geo_size_field: model_bbox span overflows Float64"))
    if dx<tol && dy<tol && dz<tol
        dx=2.0;dy=2.0
    elseif dx<tol && dy<tol
        dx=2dz;dy=2dz
    elseif dx<tol && dz<tol
        dx=2dy
    elseif dy<tol && dz<tol
        dy=2dx
    elseif dx<tol
        dx=2hypot(dy,dz)
    elseif dy<tol
        dy=2hypot(dx,dz)
    end
    all(isfinite,(dx,dy,dz)) || throw(ArgumentError(
        "build_geo_size_field: padded model_bbox span overflows Float64"))
    result=hypot(dx,dy,dz)
    (isfinite(result)&&result>0) || throw(ArgumentError(
        "build_geo_size_field: model characteristic length is not representable"))
    return result
end

function _build_geo_field(params::GeoParams,entities,root_field::Integer;
                          base_dir::AbstractString=pwd(),model_bbox=nothing,
                          max_structured_samples::Integer=10_000_000,
                          max_distance_samples::Integer=10_000_000,
                          max_octree_cells::Integer=1_000_000,
                          model_characteristic_length=nothing,context_fields=nothing,
                          geometry_tolerance=nothing,
                          entity_boundaries=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}(),
                          entity_embedded=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}(),
                          _return_model_lc::Bool=false)
    root=try Int(root_field) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("build_geo_size_field: field tag must fit Int"))
    end
    root>0 || throw(ArgumentError("build_geo_size_field: field tag must be positive"))
    bbox=if model_bbox===nothing
        nothing
    else
        raw=try
            Tuple(model_bbox)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "build_geo_size_field: model_bbox must be a six-coordinate collection"))
        end
        length(raw)==6 || throw(ArgumentError(
            "build_geo_size_field: model_bbox must contain six coordinates"))
        values=ntuple(i->_float_value(raw[i],"build_geo_size_field","model_bbox[$i]"),6)
        (values[2]>=values[1] && values[4]>=values[3] && values[6]>=values[5]) ||
            throw(ArgumentError("build_geo_size_field: model_bbox maxima must not be below minima"))
        values
    end
    model_lc=if model_characteristic_length===nothing
        if bbox===nothing
            1.0
        else
            tolerance=geometry_tolerance===nothing ?
                (isnan(params.geometry_tolerance) ? 1e-8 : params.geometry_tolerance) :
                geometry_tolerance
            _gmsh_bbox_characteristic_length(bbox,tolerance)
        end
    else
        _positive_value(model_characteristic_length,"build_geo_size_field",
                        "model_characteristic_length")
    end
    default_delta=model_lc/1e4
    (isfinite(default_delta) && default_delta>0) || throw(ArgumentError(
        "build_geo_size_field: default finite-difference step is not representable"))
    built=Dict{Int,AbstractField}(); state=Dict{Int,UInt8}()
    function build(tag::Int)
        get(state,tag,0x00)==0x02 && return built[tag]
        get(state,tag,0x00)==0x01 && throw(ArgumentError(
            "build_geo_size_field: cyclic reference involving Field[$tag]"))
        spec=get(params.fields,tag,nothing)
        spec===nothing && throw(ArgumentError(
            "build_geo_size_field: referenced Field[$tag] is not declared"))
        state[tag]=0x01
        kind=lowercase(spec.kind)
        field = if kind in ("distance","attractor")
            max_distance_samples isa Bool && throw(ArgumentError(
                "build_geo_size_field: max_distance_samples must be a positive integer"))
            max_distance_samples>0 || throw(ArgumentError(
                "build_geo_size_field: max_distance_samples must be positive"))
            max_distance_samples<=typemax(Int) || throw(ArgumentError(
                "build_geo_size_field: max_distance_samples exceeds the platform Int limit"))
            _geo_distance(spec,entities,Int(max_distance_samples))
        elseif kind=="threshold"
            _geo_known_options(spec,Set(["InField","DistMin","DistMax","SizeMin",
                                         "SizeMax","Sigmoid","StopAtDistMax",
                                         "IField","LcMin","LcMax"]))
            input_tag=_geo_infield(spec;default=0)
            hmin=_geo_float_alias(spec,"SizeMin","LcMin",0.1)
            hmax=_geo_float_alias(spec,"SizeMax","LcMax",1.0)
            dist_min=_geo_float(spec,"DistMin",1.0)
            dist_max=_geo_float(spec,"DistMax",10.0)
            sigmoid=_geo_bool(spec,"Sigmoid")
            stop_at_dist_max=_geo_bool(spec,"StopAtDistMax")
            input_tag<=0 || input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                ThresholdField(build(input_tag);dist_min=dist_min,dist_max=dist_max,
                    size_min=hmin,size_max=hmax,sigmoid=sigmoid,
                    stop_at_dist_max=stop_at_dist_max)
        elseif kind=="box"
            _geo_known_options(spec,Set(["VIn","VOut","XMin","XMax","YMin","YMax",
                                         "ZMin","ZMax","Thickness"]))
            BoxField(_geo_float(spec,"XMin",0.0),_geo_float(spec,"XMax",0.0),
                     _geo_float(spec,"YMin",0.0),_geo_float(spec,"YMax",0.0),
                     _geo_float(spec,"ZMin",0.0),_geo_float(spec,"ZMax",0.0);
                     vin=_geo_float(spec,"VIn",GMSH_MAX_SIZE),
                     vout=_geo_float(spec,"VOut",GMSH_MAX_SIZE),
                     thickness=_geo_float(spec,"Thickness",0.0))
        elseif kind=="ball"
            _geo_known_options(spec,Set(["VIn","VOut","XCenter","YCenter",
                                         "ZCenter","Radius","Thickness"]))
            BallField((_geo_float(spec,"XCenter",0.0),
                       _geo_float(spec,"YCenter",0.0),
                       _geo_float(spec,"ZCenter",0.0)),
                      _geo_float(spec,"Radius",0.0);
                      vin=_geo_float(spec,"VIn",GMSH_MAX_SIZE),
                      vout=_geo_float(spec,"VOut",GMSH_MAX_SIZE),
                      thickness=_geo_float(spec,"Thickness",0.0))
        elseif kind=="cylinder"
            _geo_known_options(spec,Set(["VIn","VOut","XCenter","YCenter",
                                         "ZCenter","XAxis","YAxis","ZAxis","Radius"]))
            CylinderField((_geo_float(spec,"XCenter",0.0),
                           _geo_float(spec,"YCenter",0.0),
                           _geo_float(spec,"ZCenter",0.0)),
                          (_geo_float(spec,"XAxis",0.0),
                           _geo_float(spec,"YAxis",0.0),
                           _geo_float(spec,"ZAxis",1.0)),
                          _geo_float(spec,"Radius",0.0);
                          vin=_geo_float(spec,"VIn",GMSH_MAX_SIZE),
                          vout=_geo_float(spec,"VOut",GMSH_MAX_SIZE))
        elseif kind=="frustum"
            _geo_known_options(spec,Set(["X1","Y1","Z1","X2","Y2","Z2",
                "InnerR1","OuterR1","InnerR2","OuterR2","InnerV1","OuterV1",
                "InnerV2","OuterV2","R1_inner","R1_outer","R2_inner","R2_outer",
                "V1_inner","V1_outer","V2_inner","V2_outer"]))
            FrustumField((_geo_float(spec,"X1",0.0),_geo_float(spec,"Y1",0.0),
                           _geo_float(spec,"Z1",1.0)),
                          (_geo_float(spec,"X2",0.0),_geo_float(spec,"Y2",0.0),
                           _geo_float(spec,"Z2",0.0));
                inner_r1=_geo_float_alias(spec,"InnerR1","R1_inner",0.0),
                outer_r1=_geo_float_alias(spec,"OuterR1","R1_outer",1.0),
                inner_r2=_geo_float_alias(spec,"InnerR2","R2_inner",0.0),
                outer_r2=_geo_float_alias(spec,"OuterR2","R2_outer",1.0),
                inner_v1=_geo_float_alias(spec,"InnerV1","V1_inner",0.1),
                outer_v1=_geo_float_alias(spec,"OuterV1","V1_outer",1.0),
                inner_v2=_geo_float_alias(spec,"InnerV2","V2_inner",0.1),
                outer_v2=_geo_float_alias(spec,"OuterV2","V2_outer",1.0))
        elseif kind in ("min","max")
            _geo_known_options(spec,Set(["FieldsList"]))
            children=AbstractField[]
            for child in _geo_field_tags(spec,"FieldsList";
                                         required=false,nonempty=false)
                (child==tag || !haskey(params.fields,child)) && continue
                push!(children,build(child))
            end
            kind=="min" ? MinSize(children) : MaxSize(children)
        elseif kind=="matheval"
            _geo_known_options(spec,Set(["F"]))
            expr=haskey(spec.options,"F") ? _geo_string_value(spec.options["F"],
                "build_geo_size_field: Field[$tag].F") : "F2 + Sin(z)"
            ast=parse_matheval(expr)
            for dep in matheval_field_tags(ast); haskey(params.fields,dep) && build(dep); end
            MathEvalField(expr; fields=built)
        elseif kind=="mathevalaniso"
            _geo_known_options(spec,Set(["M11","M22","M33","M12","M13","M23",
                                         "m11","m22","m33","m12","m13","m23"]))
            if isempty(spec.options)
                # Gmsh defaults all six entries to the unresolved expression
                # `F2 + Sin(z)`. Its scalar value is MAX_LC and the singular
                # default metric is treated as an unconstrained background by
                # the 4.15.2 meshers. Preserve that observable no-op without
                # weakening Tessella's strict SPD Metric3 contract.
                ConstantSize(GMSH_MAX_SIZE)
            else
                default_expr="F2 + Sin(z)"
                exprs=(m11=_geo_float_string(spec,"M11","m11",default_expr),
                       m22=_geo_float_string(spec,"M22","m22",default_expr),
                       m33=_geo_float_string(spec,"M33","m33",default_expr),
                       m12=_geo_float_string(spec,"M12","m12",default_expr),
                       m13=_geo_float_string(spec,"M13","m13",default_expr),
                       m23=_geo_float_string(spec,"M23","m23",default_expr))
                for e in exprs
                    for dep in matheval_field_tags(parse_matheval(e))
                        haskey(params.fields,dep) && build(dep)
                    end
                end
                MathEvalAnisoField(; m11=exprs.m11,m22=exprs.m22,m33=exprs.m33,
                                     m12=exprs.m12,m13=exprs.m13,m23=exprs.m23,
                                     fields=built)
            end
        elseif kind=="gradient"
            _geo_known_options(spec,Set(["InField","IField","Kind","Delta"]))
            input_tag=_geo_infield(spec)
            component=_geo_int(spec,"Kind";default=3)
            delta=_geo_float(spec,"Delta",default_delta)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                GradientField(build(input_tag);kind=component,delta=delta)
        elseif kind=="laplacian"
            _geo_known_options(spec,Set(["InField","IField","Delta"]))
            input_tag=_geo_infield(spec);delta=_geo_float(spec,"Delta",default_delta)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                LaplacianField(build(input_tag);delta=delta)
        elseif kind=="mean"
            _geo_known_options(spec,Set(["InField","IField","Delta"]))
            input_tag=_geo_infield(spec);delta=_geo_float(spec,"Delta",default_delta)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                MeanField(build(input_tag);delta=delta)
        elseif kind=="curvature"
            _geo_known_options(spec,Set(["InField","IField","Delta"]))
            input_tag=_geo_infield(spec);delta=_geo_float(spec,"Delta",default_delta)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                CurvatureOpField(build(input_tag);delta=delta)
        elseif kind=="maxeigenhessian"
            _geo_known_options(spec,Set(["InField","IField","Delta"]))
            input_tag=_geo_infield(spec);delta=_geo_float(spec,"Delta",default_delta)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                MaxEigenHessianField(build(input_tag);delta=delta)
        elseif kind=="lonlat"
            _geo_known_options(spec,Set(["InField","IField","FromStereo","RadiusStereo"]))
            input_tag=_geo_infield(spec);from_stereo=_geo_int(spec,"FromStereo")==1
            radius=_geo_float(spec,"RadiusStereo",6371e3)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                LonLatField(build(input_tag);from_stereo=from_stereo,radius=radius)
        elseif kind=="param"
            _geo_known_options(spec,Set(["InField","IField","FX","FY","FZ"]))
            input_tag=_geo_infield(spec)
            fx=haskey(spec.options,"FX") ? _geo_string_value(spec.options["FX"],
                "build_geo_size_field: Field[$tag].FX") : ""
            fy=haskey(spec.options,"FY") ? _geo_string_value(spec.options["FY"],
                "build_geo_size_field: Field[$tag].FY") : ""
            fz=haskey(spec.options,"FZ") ? _geo_string_value(spec.options["FZ"],
                "build_geo_size_field: Field[$tag].FZ") : ""
            if input_tag!=tag && haskey(params.fields,input_tag)
                for expr in (fx,fy,fz)
                    isempty(expr) && continue
                    for dep in matheval_field_tags(parse_matheval(expr))
                        haskey(params.fields,dep) && build(dep)
                    end
                end
                ParametricField(build(input_tag);fx=fx,fy=fy,fz=fz,fields=built)
            else
                ConstantSize(GMSH_MAX_SIZE)
            end
        elseif kind=="structured"
            _geo_known_options(spec,Set(["FileName","TextFormat","SetOutsideValue",
                                         "OutsideValue"]))
            haskey(spec.options,"FileName") || throw(ArgumentError(
                "build_geo_size_field: Structured Field[$tag] is missing FileName"))
            filename=_geo_string_value(spec.options["FileName"],
                "build_geo_size_field: Field[$tag].FileName")
            isempty(filename) && throw(ArgumentError(
                "build_geo_size_field: Field[$tag].FileName must be non-empty"))
            path=isabspath(filename) ? filename : normpath(joinpath(base_dir,filename))
            setoutside=_geo_bool(spec,"SetOutsideValue",false)
            outside=setoutside ? _geo_float(spec,"OutsideValue",GMSH_MAX_SIZE) : nothing
            StructuredField(path;text=_geo_bool(spec,"TextFormat",false),outside=outside,
                            max_samples=max_structured_samples)
        elseif kind=="octree"
            _geo_known_options(spec,Set(["InField"]))
            bbox===nothing && throw(ArgumentError(
                "build_geo_size_field: Octree Field[$tag] requires model_bbox=(xmin,xmax,ymin,ymax,zmin,zmax)"))
            OctreeField(build(_geo_int(spec,"InField";default=1)),bbox...;
                        max_level=4,max_cells=max_octree_cells)
        elseif kind=="externalprocess"
            _geo_known_options(spec,Set(["CommandLine"]))
            if !haskey(spec.options,"CommandLine")
                ConstantSize(GMSH_MAX_SIZE)
            else
                command=_geo_string_value(spec.options["CommandLine"],
                    "build_geo_size_field: Field[$tag].CommandLine")
                isempty(command) ? ConstantSize(GMSH_MAX_SIZE) :
                                   ExternalProcessField(command)
            end
        elseif kind=="restrict"
            _geo_known_options(spec,Set(["InField","IField","PointsList","CurvesList",
                "SurfacesList","VolumesList","IncludeBoundary","VerticesList","EdgesList",
                "FacesList","RegionsList","IncludeEmbedded"]))
            input_tag=_geo_infield(spec)
            input_tag==tag || !haskey(params.fields,input_tag) ?
                ConstantSize(GMSH_MAX_SIZE) :
                RestrictField(build(input_tag);
                    points=_geo_int_list_opt(spec,"PointsList","VerticesList"),
                    curves=_geo_int_list_opt(spec,"CurvesList","EdgesList"),
                    surfaces=_geo_int_list_opt(spec,"SurfacesList","FacesList"),
                    volumes=_geo_int_list_opt(spec,"VolumesList","RegionsList"),
                    include_boundary=_geo_bool(spec,"IncludeBoundary",true),
                    include_embedded=_geo_bool(spec,"IncludeEmbedded",true),
                    entity_boundaries=entity_boundaries,entity_embedded=entity_embedded)
        elseif kind=="constant"
            _geo_known_options(spec,Set(["VIn","VOut","PointsList","CurvesList","SurfacesList",
                "VolumesList","IncludeBoundary","IncludeEmbedded"]))
            ConstantField(; vin=_geo_float(spec,"VIn",GMSH_MAX_SIZE),
                            vout=_geo_float(spec,"VOut",GMSH_MAX_SIZE),
                            points=_geo_int_list_opt(spec,"PointsList"),
                            curves=_geo_int_list_opt(spec,"CurvesList"),
                            surfaces=_geo_int_list_opt(spec,"SurfacesList"),
                            volumes=_geo_int_list_opt(spec,"VolumesList"),
                            include_boundary=_geo_bool(spec,"IncludeBoundary",true),
                            include_embedded=_geo_bool(spec,"IncludeEmbedded",true),
                            entity_boundaries=entity_boundaries,entity_embedded=entity_embedded)
        elseif kind=="minaniso"
            _geo_known_options(spec,Set(["FieldsList"]))
            children=AbstractField[]
            for child in _geo_field_tags(spec,"FieldsList";
                                         required=false,nonempty=false)
                push!(children,(child==tag || !haskey(params.fields,child)) ?
                               _GmshAnisoReferencePlaceholder() : build(child))
            end
            MinAnisoField(children)
        elseif kind=="intersectaniso"
            _geo_known_options(spec,Set(["FieldsList"]))
            children=AbstractField[]
            for child in _geo_field_tags(spec,"FieldsList";
                                         required=false,nonempty=false)
                push!(children,(child==tag || !haskey(params.fields,child)) ?
                               _GmshAnisoReferencePlaceholder() : build(child))
            end
            IntersectAnisoField(children)
        elseif kind=="extend"
            _geo_known_options(spec,Set(["CurvesList","SurfacesList","DistMax",
                                         "SizeMax","Power"]))
            config=_geo_extend_config(spec)
            ((isempty(config.curves) && isempty(config.surfaces)) ||
             config.dist_max<=0 || config.size_max<=0) ?
                ConstantSize(GMSH_MAX_SIZE) :
                _geo_context_field(context_fields,spec,entities,params,config)
        elseif kind=="postview"
            _geo_known_options(spec,Set(["ViewIndex","ViewTag","CropNegativeValues",
                                         "UseClosest","IView"]))
            _geo_context_field(context_fields,spec,entities,params,
                               _geo_postview_config(spec))
        elseif kind=="attractoranisocurve"
            _geo_known_options(spec,Set(["CurvesList","EdgesList","Sampling",
                "NNodesByEdge","NumPointsPerCurve","DistMin","DistMax","dMin","dMax",
                "SizeMinTangent","SizeMaxTangent","SizeMinNormal","SizeMaxNormal",
                "lMinTangent","lMaxTangent","lMinNormal","lMaxNormal"]))
            _geo_context_field(context_fields,spec,entities,params,
                               _geo_attractor_aniso_config(spec))
        elseif kind=="boundarylayer"
            _geo_known_options(spec,Set(["CurvesList","EdgesList","FanPointsList",
                "FanNodesList","FanPointsSizesList","PointsList","NodesList","Size",
                "SizesList","Ratio","SizeFar","Thickness","Quads","IntersectMetrics",
                "AnisoMax","BetaLaw","Beta","NbLayers","ExcludedSurfacesList",
                "ExcludedFaceList","hwall_n","hwall_n_nodes","ratio","hfar","thickness"]))
            config=_geo_boundary_layer_config(
                spec,params.mesh_boundary_layer_fan_elements)
            isempty(config.curves) && isempty(config.points) ?
                ConstantSize(GMSH_MAX_SIZE) :
                _geo_context_field(context_fields,spec,entities,params,config)
        elseif kind=="automaticmeshsizefield"
            _geo_known_options(spec,Set(["p4estFileToLoad","nPointsPerCircle",
                "nPointsPerGap","hMin","hMax","hBulk","gradation","smoothing","features"]))
            _geo_context_field(context_fields,spec,entities,params,
                               _geo_automatic_config(spec))
        else
            throw(ArgumentError(
                "build_geo_size_field: unsupported field kind $(spec.kind) in Field[$tag]"))
        end
        built[tag]=field; state[tag]=0x02
        return field
    end
    field=build(root)
    return _return_model_lc ? (field,model_lc) : field
end

struct _GmshBackgroundSize{F<:AbstractField} <: AbstractSizeField
    input::F
    model_size::Float64
    size_min::Float64
    size_max::Float64
    factor::Float64
end

struct _GmshBackgroundAniso{F<:AbstractAnisoField} <: AbstractAnisoField
    input::F
    model_size::Float64
    size_min::Float64
    size_max::Float64
    factor::Float64
end

function _gmsh_background(params::GeoParams,background::AbstractField,model_lc::Float64)
    hmin=isnan(params.mesh_size_min) ? 0.0 : params.mesh_size_min
    hmax=isnan(params.mesh_size_max) ? GMSH_MAX_SIZE : params.mesh_size_max
    factor=params.mesh_size_factor
    (isfinite(hmin)&&hmin>=0) || throw(ArgumentError(
        "build_geo_size_field: Mesh.MeshSizeMin must be finite and non-negative"))
    (isfinite(hmax)&&hmax>0) || throw(ArgumentError(
        "build_geo_size_field: Mesh.MeshSizeMax must be finite and positive"))
    (isfinite(factor)&&factor>0) || throw(ArgumentError(
        "build_geo_size_field: Mesh.MeshSizeFactor must be finite and positive"))
    hmax>=hmin || throw(ArgumentError(
        "build_geo_size_field: Mesh.MeshSizeMax must be at least Mesh.MeshSizeMin"))
    isfinite(model_lc) && model_lc>0 || throw(ArgumentError(
        "build_geo_size_field: model characteristic length must be finite and positive"))
    isfinite(model_lc*factor) || throw(ArgumentError(
        "build_geo_size_field: model characteristic length times Mesh.MeshSizeFactor overflows Float64"))
    background isa AbstractAnisoField &&
        return _GmshBackgroundAniso(background,model_lc,hmin,hmax,factor)
    return _GmshBackgroundSize(background,model_lc,hmin,hmax,factor)
end

@inline function _gmsh_background_scalar(field::_GmshBackgroundSize,x,y,z,entity)
    background=_checked_field_result(field_value(field.input,x,y,z,entity),
                                     "Gmsh background field",x,y,z)
    lc=clamp(min(field.model_size,background),field.size_min,field.size_max)
    lc>0 || throw(ArgumentError(
        "Gmsh background size: non-positive element size $lc at ($x,$y,$z) " *
        "after Mesh.MeshSizeMin/Max; set a positive lower bound or compose the " *
        "intermediate field into a positive result"))
    return _checked_field_result(lc*field.factor,"Gmsh background size",x,y,z;
                                 positive=true)
end
field_value(field::_GmshBackgroundSize,x,y,z)=
    _gmsh_background_scalar(field,x,y,z,nothing)
field_value(field::_GmshBackgroundSize,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_gmsh_background_scalar(field,x,y,z,entity)

function _gmsh_background_metric(field::_GmshBackgroundAniso,x,y,z,entity)
    base=clamp(field.model_size,field.size_min,field.size_max)
    base>0 || (base=field.model_size)
    background=metric_at(field.input,x,y,z,entity)
    background isa Metric3 || throw(ArgumentError(
        "Gmsh background metric returned $(typeof(background)), not Metric3"))
    metric=_gmsh_metric_intersection(background,isotropic_metric(base))
    if field.factor!=1
        inverse_factor=inv(field.factor)
        scale=inverse_factor*inverse_factor
        isfinite(scale) && scale>0 || throw(ArgumentError(
            "Gmsh background metric factor is not representable"))
        metric=_metric_scale(metric,scale)
    end
    return metric
end
metric_at(field::_GmshBackgroundAniso,x,y,z)=
    _gmsh_background_metric(field,x,y,z,nothing)
metric_at(field::_GmshBackgroundAniso,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=_gmsh_background_metric(field,x,y,z,entity)
field_value(field::_GmshBackgroundAniso,x,y,z)=metric_size(metric_at(field,x,y,z))
field_value(field::_GmshBackgroundAniso,x,y,z,entity::Tuple{T,U}) where
        {T<:Integer,U<:Integer}=metric_size(metric_at(field,x,y,z,entity))

"""
    build_geo_size_field(params::GeoParams, entities; kwargs...) -> AbstractSizeField

Build `params.background_field` with the strict field-graph builder documented by
`_build_geo_field`, then apply `Mesh.MeshSizeMin`, `Mesh.MeshSizeMax`, and
`Mesh.MeshSizeFactor` in Gmsh order.
"""
function build_geo_size_field(params::GeoParams,entities;kwargs...)
    params.background_field>0 || throw(ArgumentError(
        "build_geo_size_field: .geo input has no Background Field"))
    background,model_lc=_build_geo_field(params,entities,params.background_field;
                                         kwargs...,_return_model_lc=true)
    return _gmsh_background(params,background,model_lc)
end

"""
    build_geo_boundary_layer_fields(params, entities; kwargs...) -> Tuple

Build the fields selected by `.geo` `BoundaryLayer Field = ...` declarations.
Boundary-layer selections are independent of `Background Field`, are not wrapped
in global scalar size bounds, and remain in declaration order. Full Gmsh boundary
layer topology requires model-aware objects supplied through the same
`context_fields` resolver accepted by [`build_geo_size_field`](@ref).
"""
function build_geo_boundary_layer_fields(params::GeoParams,entities;kwargs...)
    values=AbstractField[]
    for tag in params.boundary_layer_fields
        spec=get(params.fields,tag,nothing)
        spec===nothing && throw(ArgumentError(
            "build_geo_boundary_layer_fields: Field[$tag] is not declared"))
        lowercase(spec.kind)=="boundarylayer" || throw(ArgumentError(
            "build_geo_boundary_layer_fields: Field[$tag] has kind $(spec.kind), not BoundaryLayer"))
        push!(values,_build_geo_field(params,entities,tag;kwargs...))
    end
    return Tuple(values)
end

"""
    size_at(field, x, y, z)
    size_at(field, point)

Evaluate a mesh-size field and require a finite, strictly positive result.
"""
@inline function size_at(field::AbstractSizeField,x,y,z)
    return _checked_field_result(field_value(field,x,y,z),"size_at",x,y,z;positive=true)
end

@inline function field_value(field::AbstractField,p)
    q=_point3(p,"field_value")
    return field_value(field,q[1],q[2],q[3])
end

@inline function size_at(field::AbstractSizeField,p)
    q=_point3(p,"size_at")
    return size_at(field,q[1],q[2],q[3])
end

end # module SizeField
