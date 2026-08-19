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

using ..MeshTypes: Mesh
using ..IO: GeoParams, GeoFieldSpec

export AbstractField, AbstractSizeField, ConstantSize, FunctionSize
export DistanceField, ThresholdField, BoxField, BallField, CylinderField, FrustumField
export MinSize, MaxSize, BoundedSize
export MinField, MaxField, field_value, size_at, build_geo_size_field, GMSH_MAX_SIZE

"""Finite sentinel used by Gmsh fields for "do not constrain the mesh size"."""
const GMSH_MAX_SIZE = 1.0e22

"""Base type for scalar fields evaluable with [`field_value`](@ref)."""
abstract type AbstractField end

"""
Base type for strictly positive scalar fields that meshing kernels may consume
through [`size_at`](@ref).
"""
abstract type AbstractSizeField <: AbstractField end

"""
    field_value(field, x, y, z)
    field_value(field, point)

Evaluate a scalar field. Field implementations reject non-finite results; generic
scalar fields may legitimately return zero or a negative value.
"""
function field_value end

@inline function _float_value(x::Real, caller::AbstractString, name::AbstractString)
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

@inline function _point3(p, caller::AbstractString)
    n = try
        length(p)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: a point must be an indexable coordinate collection"))
    end
    n >= 2 || throw(ArgumentError("$caller: a point needs at least two coordinates"))
    q = try
        (Float64(p[1]), Float64(p[2]), n >= 3 ? Float64(p[3]) : 0.0)
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
    hmin=_positive_value(size_min,"ThresholdField","size_min")
    hmax=_positive_value(size_max,"ThresholdField","size_max")
    return ThresholdField{F}(input,dmin,dmax,hmin,hmax,sigmoid,stop_at_dist_max)
end

function field_value(field::ThresholdField, x, y, z)
    d=_checked_field_result(field_value(field.input,x,y,z),"ThresholdField input",x,y,z)
    field.stop_at_dist_max && d>=field.dist_max && return GMSH_MAX_SIZE
    r=clamp((d-field.dist_min)/(field.dist_max-field.dist_min),0.0,1.0)
    if field.sigmoid
        e=exp(12r-6); r=e/(1+e)
    end
    return muladd(field.size_max-field.size_min,r,field.size_min)
end

"""
    BoxField(xmin, xmax, ymin, ymax, zmin, zmax; vin, vout, thickness=0)

Gmsh `Box` field.  Returns `vin` inside the closed axis-aligned box and
`vout` outside.  A positive `thickness` linearly interpolates in the exterior
layer using Euclidean distance to the box.
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
    xa<=xb || throw(ArgumentError("BoxField: require xmin <= xmax"))
    ya<=yb || throw(ArgumentError("BoxField: require ymin <= ymax"))
    za<=zb || throw(ArgumentError("BoxField: require zmin <= zmax"))
    vi=_positive_value(vin,"BoxField","vin"); vo=_positive_value(vout,"BoxField","vout")
    th=_float_value(thickness,"BoxField","thickness")
    th>=0 || throw(ArgumentError("BoxField: thickness must be non-negative"))
    return BoxField(xa,xb,ya,yb,za,zb,vi,vo,th)
end

function field_value(field::BoxField, x, y, z)
    x=_float_value(x,"BoxField","x"); y=_float_value(y,"BoxField","y")
    z=_float_value(z,"BoxField","z")
    inside=field.xmin<=x<=field.xmax && field.ymin<=y<=field.ymax &&
           field.zmin<=z<=field.zmax
    inside && return field.vin
    if field.thickness>0
        dx=max(field.xmin-x,0.0,x-field.xmax)
        dy=max(field.ymin-y,0.0,y-field.ymax)
        dz=max(field.zmin-z,0.0,z-field.zmax)
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
`radius + thickness` is interpolated linearly from `vin` to `vout`.
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
    r>=0 || throw(ArgumentError("BallField: radius must be non-negative"))
    vi=_positive_value(vin,"BallField","vin")
    vo=_positive_value(vout,"BallField","vout")
    th=_float_value(thickness,"BallField","thickness")
    th>=0 || throw(ArgumentError("BallField: thickness must be non-negative"))
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
strict, matching Gmsh.
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
    (isfinite(half)&&half>0) ||
        throw(ArgumentError("CylinderField: axis must have finite positive length"))
    r=_float_value(radius,"CylinderField","radius")
    r>=0 || throw(ArgumentError("CylinderField: radius must be non-negative"))
    vi=_positive_value(vin,"CylinderField","vin")
    vo=_positive_value(vout,"CylinderField","vout")
    return CylinderField(c,ntuple(i -> a[i]/half,3),half,r,vi,vo)
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

Gmsh `Frustum` field. Within the annular frustum from `p1` to `p2`, interpolate
bilinearly between the mesh sizes on the inner and outer radii at both endpoints.
Outside the axial interval or radial annulus, return [`GMSH_MAX_SIZE`](@ref).
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
    (isfinite(len)&&len>0) ||
        throw(ArgumentError("FrustumField: endpoints must be distinct with finite separation"))
    r1i=_float_value(inner_r1,"FrustumField","inner_r1")
    r1o=_float_value(outer_r1,"FrustumField","outer_r1")
    r2i=_float_value(inner_r2,"FrustumField","inner_r2")
    r2o=_float_value(outer_r2,"FrustumField","outer_r2")
    (r1i>=0&&r2i>=0) || throw(ArgumentError("FrustumField: inner radii must be non-negative"))
    r1o>r1i || throw(ArgumentError("FrustumField: require outer_r1 > inner_r1"))
    r2o>r2i || throw(ArgumentError("FrustumField: require outer_r2 > inner_r2"))
    v1i=_positive_value(inner_v1,"FrustumField","inner_v1")
    v1o=_positive_value(outer_v1,"FrustumField","outer_v1")
    v2i=_positive_value(inner_v2,"FrustumField","inner_v2")
    v2o=_positive_value(outer_v2,"FrustumField","outer_v2")
    return FrustumField(a,ntuple(i -> axis[i]/len,3),len,r1i,r1o,r2i,r2o,
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

"""Pointwise maximum of one or more scalar fields (Gmsh `Max`)."""
struct MaxSize{F<:Tuple} <: AbstractSizeField
    fields::F
end

function MinSize(fields::Union{Tuple,AbstractVector})
    values=Tuple(fields)
    isempty(values) && throw(ArgumentError("MinSize: need at least one field"))
    all(child -> child isa AbstractField,values) ||
        throw(ArgumentError("MinSize: every input must be an AbstractField"))
    return MinSize{typeof(values)}(values)
end


function MaxSize(fields::Union{Tuple,AbstractVector})
    values=Tuple(fields)
    isempty(values) && throw(ArgumentError("MaxSize: need at least one field"))
    all(child -> child isa AbstractField,values) ||
        throw(ArgumentError("MaxSize: every input must be an AbstractField"))
    return MaxSize{typeof(values)}(values)
end

const MinField = MinSize
const MaxField = MaxSize

@inline _min_value(::Tuple{},x,y,z)=GMSH_MAX_SIZE
@inline function _min_value(fields::Tuple,x,y,z)
    child=_checked_field_result(field_value(first(fields),x,y,z),"MinSize input",x,y,z)
    return min(child,_min_value(Base.tail(fields),x,y,z))
end
@inline field_value(field::MinSize,x,y,z)=_min_value(field.fields,x,y,z)

@inline _max_value(::Tuple{},x,y,z)=-GMSH_MAX_SIZE
@inline function _max_value(fields::Tuple,x,y,z)
    child=_checked_field_result(field_value(first(fields),x,y,z),"MaxSize input",x,y,z)
    return max(child,_max_value(Base.tail(fields),x,y,z))
end
@inline field_value(field::MaxSize,x,y,z)=_max_value(field.fields,x,y,z)

"""
    BoundedSize(field; size_min=0, size_max=GMSH_MAX_SIZE, factor=1)

Apply Gmsh's final size policy: clamp the input field into
`[size_min,size_max]`, then multiply by `factor`.  This also turns a generic
scalar field such as `DistanceField` into a valid mesh-size field when
`size_min > 0`. A zero lower bound is accepted only when the input is already an
`AbstractSizeField`.
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
    (hmin>0 || input isa AbstractSizeField) || throw(ArgumentError(
        "BoundedSize: size_min must be positive for a generic scalar input"))
    hmax=_positive_value(size_max,"BoundedSize","size_max")
    hmax>=hmin || throw(ArgumentError("BoundedSize: require size_max >= size_min"))
    f=_positive_value(factor,"BoundedSize","factor")
    isfinite(hmax*f) || throw(ArgumentError("BoundedSize: size_max * factor must be finite"))
    return BoundedSize{F}(input,hmin,hmax,f)
end

function field_value(field::BoundedSize,x,y,z)
    value=_checked_field_result(field_value(field.input,x,y,z),"BoundedSize input",x,y,z)
    return field.factor*clamp(value,field.size_min,field.size_max)
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
    haskey(spec.options,name) && haskey(spec.options,deprecated) &&
        throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] sets both $name and $deprecated"))
    return haskey(spec.options,name) ? _geo_float(spec,name,default) :
           _geo_float(spec,deprecated,default)
end

function _geo_int(spec::GeoFieldSpec,name::String; required::Bool=false,default::Int=0)
    if !haskey(spec.options,name)
        required && throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] is missing $name"))
        return default
    end
    raw=strip(spec.options[name]); value=tryparse(Int,raw)
    value!==nothing || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be an integer literal (got $raw)"))
    return value::Int
end

function _geo_bool(spec::GeoFieldSpec,name::String,default::Bool=false)
    haskey(spec.options,name) || return default
    raw=lowercase(strip(spec.options[name]))
    raw in ("1","true") && return true
    raw in ("0","false") && return false
    throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be 0/1 or true/false (got $raw)"))
end

function _geo_list(spec::GeoFieldSpec,name::String; required::Bool=false)
    if !haskey(spec.options,name)
        required && throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)] is missing $name"))
        return String[]
    end
    raw=strip(spec.options[name])
    (startswith(raw,"{") && endswith(raw,"}")) || throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name must be a brace-delimited list"))
    body=strip(raw[2:end-1]); isempty(body) && return String[]
    values=strip.(split(body,','))
    any(isempty,values) && throw(ArgumentError(
        "build_geo_size_field: Field[$(spec.tag)].$name contains an empty entry"))
    return values
end

function _geo_field_tags(spec::GeoFieldSpec,name::String)
    raws=_geo_list(spec,name;required=true); tags=Int[]
    for raw in raws
        tag=tryparse(Int,raw)
        (tag!==nothing && tag>0) || throw(ArgumentError(
            "build_geo_size_field: Field[$(spec.tag)].$name needs positive field tags (got $raw)"))
        push!(tags,tag)
    end
    isempty(tags) && throw(ArgumentError(
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

function _geo_distance(spec::GeoFieldSpec,entities)
    allowed=Set(["PointsList","CurvesList","SurfacesList","Sampling",
                 "NodesList","EdgesList","FacesList","NNodesByEdge",
                 "NumPointsPerCurve"])
    _geo_known_options(spec,allowed)
    for name in ("Sampling","NNodesByEdge","NumPointsPerCurve")
        haskey(spec.options,name) && _geo_int(spec,name)>0 || !haskey(spec.options,name) ||
            throw(ArgumentError("build_geo_size_field: Field[$(spec.tag)].$name must be positive"))
    end
    points=NTuple{3,Float64}[]
    segments=Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
    triangles=Tuple{NTuple{3,Float64},NTuple{3,Float64},NTuple{3,Float64}}[]
    groups=((0,("PointsList","NodesList")),
            (1,("CurvesList","EdgesList")),
            (2,("SurfacesList","FacesList")))
    for (dim,names) in groups, option in names, raw in _geo_list(spec,option)
        for mesh in _geo_resolve(entities,dim,raw)
            if dim==0
                append!(points,(_mesh_point(mesh,i) for i in axes(mesh.coords,2)))
            elseif dim==1
                size(mesh.segs,2)>0 || throw(ArgumentError(
                    "build_geo_size_field: curve entity $raw contains no segments"))
                append!(segments,((_mesh_point(mesh,mesh.segs[1,i]),
                                   _mesh_point(mesh,mesh.segs[2,i])) for i in axes(mesh.segs,2)))
            else
                size(mesh.tris,2)>0 || throw(ArgumentError(
                    "build_geo_size_field: surface entity $raw contains no triangles"))
                append!(triangles,((_mesh_point(mesh,mesh.tris[1,i]),
                                    _mesh_point(mesh,mesh.tris[2,i]),
                                    _mesh_point(mesh,mesh.tris[3,i])) for i in axes(mesh.tris,2)))
            end
        end
    end
    isempty(points)&&isempty(segments)&&isempty(triangles) && throw(ArgumentError(
        "build_geo_size_field: Distance Field[$(spec.tag)] has no resolved targets"))
    return DistanceField(;points=points,segments=segments,triangles=triangles)
end

"""
    build_geo_size_field(params::GeoParams, entities) -> AbstractSizeField

Build the parsed `.geo` `Background Field` and apply `Mesh.MeshSizeMin`,
`Mesh.MeshSizeMax`, and `Mesh.MeshSizeFactor` in Gmsh order. `entities` is either
a dictionary or a callable resolving `(dimension, reference)` to a `Mesh` (or a
vector of meshes). Symbolic references such as `surface_group[]` are normalized to
`surface_group`; numeric entity tags are also accepted.

The strict builder currently supports `Distance`/`Attractor`, `Threshold`, `Box`,
`Ball`, `Cylinder`, `Frustum`, `Min`, and `Max`. Unsupported field kinds, options,
expressions, missing entities, cycles, and invalid references raise `ArgumentError`
instead of being ignored.
"""
function build_geo_size_field(params::GeoParams,entities)
    params.background_field>0 || throw(ArgumentError(
        "build_geo_size_field: .geo input has no Background Field"))
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
            _geo_distance(spec,entities)
        elseif kind=="threshold"
            _geo_known_options(spec,Set(["InField","DistMin","DistMax","SizeMin",
                                         "SizeMax","Sigmoid","StopAtDistMax",
                                         "IField","LcMin","LcMax"]))
            haskey(spec.options,"InField") && haskey(spec.options,"IField") &&
                throw(ArgumentError("build_geo_size_field: Field[$tag] sets both InField and IField"))
            input_tag=haskey(spec.options,"InField") ? _geo_int(spec,"InField";required=true) :
                      _geo_int(spec,"IField";required=true)
            haskey(spec.options,"SizeMin") && haskey(spec.options,"LcMin") &&
                throw(ArgumentError("build_geo_size_field: Field[$tag] sets both SizeMin and LcMin"))
            haskey(spec.options,"SizeMax") && haskey(spec.options,"LcMax") &&
                throw(ArgumentError("build_geo_size_field: Field[$tag] sets both SizeMax and LcMax"))
            hmin=haskey(spec.options,"SizeMin") ? _geo_float(spec,"SizeMin",0.1) :
                 _geo_float(spec,"LcMin",0.1)
            hmax=haskey(spec.options,"SizeMax") ? _geo_float(spec,"SizeMax",1.0) :
                 _geo_float(spec,"LcMax",1.0)
            ThresholdField(build(input_tag);
                dist_min=_geo_float(spec,"DistMin",1.0),
                dist_max=_geo_float(spec,"DistMax",10.0),size_min=hmin,size_max=hmax,
                sigmoid=_geo_bool(spec,"Sigmoid"),
                stop_at_dist_max=_geo_bool(spec,"StopAtDistMax"))
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
            children=AbstractField[build(child) for child in _geo_field_tags(spec,"FieldsList")]
            kind=="min" ? MinSize(children) : MaxSize(children)
        else
            throw(ArgumentError(
                "build_geo_size_field: unsupported field kind $(spec.kind) in Field[$tag]"))
        end
        built[tag]=field; state[tag]=0x02
        return field
    end
    background=build(params.background_field)
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
    return BoundedSize(background;size_min=hmin,size_max=hmax,factor=factor)
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
