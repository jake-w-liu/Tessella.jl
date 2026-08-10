"""
    Mesh1D

Stage-2 one-dimensional edge meshing under a size field (PLAN.md §3 "Mesh1D";
gmsh `meshGEdge` in spirit). A curve `γ: [t0,t1] → ℝ³` is discretized so that
consecutive node spacing follows the local target size `h(x)`: nodes are placed
at equal increments of the **metric length** `∫ |γ'(t)| / h(γ(t)) dt`, which is
the standard size-conforming (graded) 1-D mesh.

The curve is supplied as a callable `γ(t) -> (x,y,z)`. Arc length and the metric
integral are evaluated by dense sampling + trapezoidal/piecewise-linear
accumulation, then node parameters are found by inverse interpolation of the
cumulative metric length — an approach whose error is controlled by the sample
count and verified against analytic arc length in the tests.
"""
module Mesh1D

using ..SizeField: AbstractSizeField, size_at, ConstantSize

export mesh_curve, mesh_segment, curve_length, metric_length

@inline _pt3(p) = (Float64(p[1]), Float64(p[2]), length(p) >= 3 ? Float64(p[3]) : 0.0)
@inline _dist3(a, b) = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)

"""
    curve_length(γ; t0=0.0, t1=1.0, nsample=2000) -> Float64

Euclidean arc length of `γ` over `[t0,t1]` by dense polyline sampling.
"""
function curve_length(γ; t0::Real=0.0, t1::Real=1.0, nsample::Integer=2000)
    prev = _pt3(γ(t0)); L = 0.0
    @inbounds for i in 1:nsample
        t = t0 + (t1-t0)*i/nsample
        p = _pt3(γ(t)); L += _dist3(prev, p); prev = p
    end
    return L
end

"""
    metric_length(γ, sf; t0=0.0, t1=1.0, nsample=2000) -> Float64

Metric length `∫ |γ'| / h`, i.e. the ideal number of edges under size field `sf`.
"""
function metric_length(γ, sf::AbstractSizeField; t0::Real=0.0, t1::Real=1.0, nsample::Integer=2000)
    prev = _pt3(γ(t0)); M = 0.0
    @inbounds for i in 1:nsample
        t = t0 + (t1-t0)*i/nsample
        p = _pt3(γ(t))
        hmid = 0.5*(size_at(sf, prev) + size_at(sf, p))
        M += _dist3(prev, p) / hmid
        prev = p
    end
    return M
end

"""
    mesh_curve(γ, sf; t0=0.0, t1=1.0, closed=false, nsample=2000)
        -> (points::Vector{NTuple{3,Float64}}, params::Vector{Float64})

Graded discretization of `γ` under size field `sf`. Nodes are placed at equal
metric-length increments; the number of edges is `max(1, round(metric_length))`
(`max(3, …)` for a closed curve so the loop is non-degenerate). Returns the node
coordinates and their curve parameters (endpoints included; for a closed curve
the last node is *not* duplicated).
"""
function mesh_curve(γ, sf::AbstractSizeField; t0::Real=0.0, t1::Real=1.0,
                    closed::Bool=false, nsample::Integer=2000)
    t0 = Float64(t0); t1 = Float64(t1)
    # cumulative metric length at each sample parameter
    ts = Vector{Float64}(undef, nsample+1)
    cum = Vector{Float64}(undef, nsample+1)
    prev = _pt3(γ(t0)); ts[1] = t0; cum[1] = 0.0
    @inbounds for i in 1:nsample
        t = t0 + (t1-t0)*i/nsample
        p = _pt3(γ(t))
        hmid = 0.5*(size_at(sf, prev) + size_at(sf, p))
        cum[i+1] = cum[i] + _dist3(prev, p)/hmid
        ts[i+1] = t; prev = p
    end
    Mtot = cum[end]
    nedge = max(closed ? 3 : 1, round(Int, Mtot))
    nnode = closed ? nedge : nedge + 1
    pts = Vector{NTuple{3,Float64}}(undef, nnode)
    par = Vector{Float64}(undef, nnode)
    @inbounds for k in 0:nnode-1
        target = Mtot * k / nedge
        tk = _invert_cum(cum, ts, target)
        par[k+1] = tk; pts[k+1] = _pt3(γ(tk))
    end
    return pts, par
end

# parameter at which cumulative metric length == target (piecewise-linear inverse)
@inline function _invert_cum(cum::Vector{Float64}, ts::Vector{Float64}, target::Float64)
    n = length(cum)
    target <= cum[1] && return ts[1]
    target >= cum[n] && return ts[n]
    lo = 1; hi = n                          # binary search for the bracketing interval
    while hi - lo > 1
        mid = (lo+hi) >>> 1
        (cum[mid] <= target) ? (lo = mid) : (hi = mid)
    end
    c0 = cum[lo]; c1 = cum[hi]
    w = c1 == c0 ? 0.0 : (target - c0)/(c1 - c0)
    return ts[lo] + w*(ts[hi] - ts[lo])
end

"""
    mesh_segment(a, b, sf; nsample=2000) -> (points, params)

Convenience: graded mesh of the straight segment `a → b` (3-D points) under `sf`.
"""
function mesh_segment(a, b, sf::AbstractSizeField; nsample::Integer=2000)
    A = _pt3(a); B = _pt3(b)
    γ(t) = (A[1]+t*(B[1]-A[1]), A[2]+t*(B[2]-A[2]), A[3]+t*(B[3]-A[3]))
    return mesh_curve(γ, sf; t0=0.0, t1=1.0, nsample=nsample)
end

end # module Mesh1D
