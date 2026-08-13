"""
    HighOrder

Stage-6 high-order (curved) elements (PLAN.md §4 Stage 6). This module generates
**quadratic (10-node, P2) tetrahedra** from a linear tet [`Mesh`](@ref) by adding
one shared node at the midpoint of every edge, in gmsh's type-11 node order
(4 corners, then the 6 edge nodes on edges `(1,2),(2,3),(3,1),(1,4),(2,4),(3,4)`).

Straight-sided P2 is geometrically identical to the linear mesh (edge nodes are
exact midpoints) — the total volume is unchanged — and is what a P2 FEM solver
consumes. [`curve_to_cylinder!`](@ref) / [`curve_to_surface!`](@ref) then curve the
edge nodes onto the true boundary geometry for curved-boundary accuracy.

Curving is done safely: only mid-nodes on genuine **boundary-surface edges** are
moved (never interior chords — moving those tangles the volume), and every
projection is reverted unless every incident element has a strictly positive
exact cubic-Bernstein Jacobian certificate over the whole reference tetrahedron.
Thus a between-sample fold cannot pass the guard.
"""
module HighOrder

using ..MeshTypes: Mesh, node, tet_volume, boundary_faces, validate
import ..MeshTypes: nnodes, ntets     # extended for P2Mesh
using Printf: @printf

export P2Mesh, p2_tetmesh, p2_volume, write_msh_p2,
       curve_to_cylinder!, curve_to_surface!, p2_min_jacobian

"""
    P2Mesh(coords, tet10)

Quadratic tet mesh: `coords` is `3 × N` (original corners followed by edge-mid
nodes); `tet10` is `10 × ntet` in gmsh type-11 order.
"""
struct P2Mesh
    coords::Matrix{Float64}
    tet10::Matrix{Int32}
    function P2Mesh(coords::AbstractMatrix{<:Real}, tet10::AbstractMatrix{<:Integer})
        size(coords, 1) == 3 || throw(ArgumentError("P2Mesh: coords must be 3 × nnodes"))
        size(tet10, 1) == 10 || throw(ArgumentError("P2Mesh: tet10 must be 10 × ntets"))
        nn = size(coords, 2)
        size(tet10,2)<=typemax(Int32) ||
            throw(ArgumentError("P2Mesh: tetrahedron count exceeds the Int32 topology limit"))
        nn <= typemax(Int32) ||
            throw(ArgumentError("P2Mesh: $nn nodes exceed the Int32 indexing limit"))
        C = try
            Matrix{Float64}(coords)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("P2Mesh: coordinates must be representable as Float64"))
        end
        @inbounds for i in axes(C,2), d in 1:3
            isfinite(C[d,i]) ||
                throw(ArgumentError("P2Mesh: node $i has a non-finite coordinate"))
        end
        T = try
            Matrix{Int32}(tet10)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("P2Mesh: connectivity must fit Int32"))
        end
        @inbounds for t in axes(T,2)
            for k in 1:10
                1 <= T[k,t] <= nn ||
                    throw(ArgumentError("P2Mesh: tet $t references node $(T[k,t]) outside 1:$nn"))
                for j in 1:k-1
                    T[j,t] != T[k,t] ||
                        throw(ArgumentError("P2Mesh: tet $t repeats node $(T[k,t])"))
                end
            end
        end
        if size(T,2)>0
            linear=Mesh(C;tets=Matrix(@view T[1:4,:]))
            diag=validate(linear)
            diag.ok || throw(ArgumentError("P2Mesh: linear corner complex is invalid — " *
                                            join(diag.messages,"; ")))
            slots=((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,2,4),(10,3,4))
            edge_mid=Dict{Tuple{Int32,Int32},Int32}()
            mid_edge=Dict{Int32,Tuple{Int32,Int32}}()
            corners=Set{Int32}(@view T[1:4,:])
            @inbounds for t in axes(T,2),(slot,i,j) in slots
                a=T[i,t];b=T[j,t];key=minmax(a,b);mid=T[slot,t]
                old=get(edge_mid,key,Int32(0))
                (old==0||old==mid) || throw(ArgumentError(
                    "P2Mesh: edge $key uses inconsistent mid-nodes $old and $mid"))
                oldedge=get(mid_edge,mid,nothing)
                (oldedge===nothing||oldedge==key) || throw(ArgumentError(
                    "P2Mesh: mid-node $mid is reused by distinct edges $oldedge and $key"))
                edge_mid[key]=mid;mid_edge[mid]=key
            end
            for mid in keys(mid_edge)
                mid in corners && throw(ArgumentError(
                    "P2Mesh: node $mid is used as both a corner and an edge mid-node"))
            end
        end
        new(C, T)
    end
end

nnodes(p::P2Mesh) = size(p.coords, 2)
ntets(p::P2Mesh) = size(p.tet10, 2)

"""
    p2_tetmesh(m::Mesh) -> P2Mesh

Convert a linear tet mesh to straight-sided quadratic tets, sharing one mid-node
per edge (so the node count grows by exactly the number of unique edges).
"""
function p2_tetmesh(m::Mesh)
    d = validate(m)
    d.ok || throw(ArgumentError("p2_tetmesh: input mesh is invalid — " * join(d.messages, "; ")))
    nn = nnodes(m); nt = ntets(m)
    nt == 0 && return P2Mesh(copy(m.coords), Matrix{Int32}(undef, 10, 0))
    edgeid = Dict{Tuple{Int32,Int32}, Int32}()
    xs=Float64[]; ys=Float64[]; zs=Float64[]
    try Base.checked_mul(6,nt) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("p2_tetmesh: edge-record count overflows Int"))
    end
    function mid(a::Int32, b::Int32)
        key = a < b ? (a, b) : (b, a)
        get!(edgeid, key) do
            nn + length(xs) < typemax(Int32) ||
                throw(ArgumentError("p2_tetmesh: quadratic node count exceeds Int32"))
            pa = node(m, a); pb = node(m, b)
            mx=pa[1]/2+pb[1]/2; my=pa[2]/2+pb[2]/2; mz=pa[3]/2+pb[3]/2
            (isfinite(mx)&&isfinite(my)&&isfinite(mz)) ||
                throw(ArgumentError("p2_tetmesh: edge ($a,$b) has a non-finite midpoint"))
            midp=(mx,my,mz)
            (midp!=pa && midp!=pb) ||
                throw(ArgumentError("p2_tetmesh: edge ($a,$b) midpoint is below Float64 coordinate resolution"))
            push!(xs,mx); push!(ys,my); push!(zs,mz)
            Int32(nn + length(xs))
        end
    end
    tet10 = Matrix{Int32}(undef, 10, nt)
    @inbounds for t in 1:nt
        v1=m.tets[1,t]; v2=m.tets[2,t]; v3=m.tets[3,t]; v4=m.tets[4,t]
        tet10[1,t]=v1; tet10[2,t]=v2; tet10[3,t]=v3; tet10[4,t]=v4
        tet10[5,t]=mid(v1,v2); tet10[6,t]=mid(v2,v3); tet10[7,t]=mid(v3,v1)
        tet10[8,t]=mid(v1,v4); tet10[9,t]=mid(v2,v4); tet10[10,t]=mid(v3,v4)
    end
    N = nn + length(xs)
    coords = Matrix{Float64}(undef, 3, N)
    @inbounds for i in 1:nn; p=node(m,i); coords[1,i]=p[1]; coords[2,i]=p[2]; coords[3,i]=p[3]; end
    @inbounds for i in eachindex(xs); coords[1,nn+i]=xs[i]; coords[2,nn+i]=ys[i]; coords[3,nn+i]=zs[i]; end
    return P2Mesh(coords, tet10)
end

"""
Total isoparametric volume of a P2 mesh. The P2 Jacobian determinant is cubic, so
its exact cubic Bernstein coefficients integrate to their sum divided by 120.
Every element first has to pass the exact global Bernstein positivity certificate;
the positive coefficients are then accumulated with 256-bit arithmetic before the
single `Float64` result conversion. A folded element cannot be hidden by absolute
values or sampled quadrature. For a straight-sided P2 tet this reduces to its
linear corner-tet volume.
"""
function p2_volume(p::P2Mesh)
    return setprecision(BigFloat, 256) do
        total = BigFloat(0)
        @inbounds for t in 1:ntets(p)
            sums, scaleexp = _p2_bernstein_coeffs(p.coords,p.tet10,t)
            all(>(0), sums) ||
                throw(ArgumentError("p2_volume: tet $t lacks a global positive-Jacobian certificate"))
            for i in eachindex(sums)
                total += ldexp(BigFloat(sums[i]) / _B3_NPERM[i], scaleexp) / 120
            end
        end
        v = Float64(total)
        (isfinite(v) && (v > 0 || ntets(p) == 0)) ||
            throw(ArgumentError("p2_volume: total volume is not representable as a positive Float64"))
        v
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# P2 element validity — exact global cubic-Bernstein Jacobian certificate.
# ════════════════════════════════════════════════════════════════════════════════
# Reference tet corners N1=(0,0,0),N2=(1,0,0),N3=(0,1,0),N4=(0,0,1);
# barycentric L1=1-r-s-t, L2=r, L3=s, L4=t.  Shape functions (type-11 order):
#   corner i : Li(2Li-1)  ->  grad = (4Li-1)·grad(Li)
#   edge(a,b): 4·La·Lb    ->  grad = 4·(La·grad(Lb) + Lb·grad(La))
# with grad(L1)=(-1,-1,-1), grad(L2)=(1,0,0), grad(L3)=(0,1,0), grad(L4)=(0,0,1).
@inline function _p2_grads(r::Float64, s::Float64, t::Float64)
    L1=1.0-r-s-t; L2=r; L3=s; L4=t
    g1=(-(4L1-1), -(4L1-1), -(4L1-1))
    g2=(4L2-1, 0.0, 0.0)
    g3=(0.0, 4L3-1, 0.0)
    g4=(0.0, 0.0, 4L4-1)
    # edge(a,b) = 4(La·grad Lb + Lb·grad La)
    g5=(4*(L1*1.0 + L2*(-1.0)), 4*(L2*(-1.0)),           4*(L2*(-1.0)))          # (1,2)
    g6=(4*(L3*1.0),             4*(L2*1.0),               0.0)                    # (2,3)
    g7=(4*(L3*(-1.0)),          4*(L1*1.0 + L3*(-1.0)),   4*(L3*(-1.0)))          # (3,1)
    g8=(4*(L4*(-1.0)),          4*(L4*(-1.0)),            4*(L1*1.0 + L4*(-1.0))) # (1,4)
    g9=(4*(L4*1.0),             0.0,                      4*(L2*1.0))             # (2,4)
    g10=(0.0,                   4*(L4*1.0),               4*(L3*1.0))             # (3,4)
    return (g1,g2,g3,g4,g5,g6,g7,g8,g9,g10)
end

@inline function _detJ(coords, tet10, t::Integer, grads)
    j11=0.0;j12=0.0;j13=0.0;j21=0.0;j22=0.0;j23=0.0;j31=0.0;j32=0.0;j33=0.0
    @inbounds for k in 1:10
        v=tet10[k,t]; x=coords[1,v]; y=coords[2,v]; z=coords[3,v]; g=grads[k]
        j11+=x*g[1]; j12+=x*g[2]; j13+=x*g[3]
        j21+=y*g[1]; j22+=y*g[2]; j23+=y*g[3]
        j31+=z*g[1]; j32+=z*g[2]; j33+=z*g[3]
    end
    return j11*(j22*j33-j23*j32) - j12*(j21*j33-j23*j31) + j13*(j21*j32-j22*j31)
end

# At the four reference vertices J is affine.  Expanding
# det(ΣλᵢJᵢ) by its three columns groups the 4³ mixed determinants into the
# 20 degree-3 tetrahedral Bernstein coefficients.  Positivity of every exact
# coefficient proves det(J)>0 everywhere by the Bernstein convex-hull property.
const _REF_VERTEX_GRADS = let rv=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    ntuple(4) do i
        g=_p2_grads(rv[i]...)
        ntuple(k -> (Int(g[k][1]),Int(g[k][2]),Int(g[k][3])),10)
    end
end
const _B3_MULTI = Tuple((a,b,c,3-a-b-c) for a in 0:3 for b in 0:3-a for c in 0:3-a-b)
const _B3_NPERM = Tuple(6 ÷ (factorial(q[1])*factorial(q[2])*factorial(q[3])*factorial(q[4]))
                         for q in _B3_MULTI)
const _B3_GROUP = let G=Array{Int8}(undef,4,4,4)
    for i in 1:4,j in 1:4,k in 1:4
        q=(count(==(1),(i,j,k)),count(==(2),(i,j,k)),
           count(==(3),(i,j,k)),count(==(4),(i,j,k)))
        G[i,j,k]=Int8(findfirst(==(q),_B3_MULTI))
    end
    G
end

# Exact IEEE-754 decomposition x = mantissa*2^exponent (finite x only).
@inline function _dyadic_parts(x::Float64)
    u=reinterpret(UInt64,x); eb=Int((u>>52)&0x7ff); frac=u&0x000fffffffffffff
    if eb==0
        m=BigInt(frac); e=-1074
    else
        m=BigInt(frac|0x0010000000000000); e=eb-1023-52
    end
    (u>>63)!=0 && (m=-m)
    return m,e
end

function _integer_element_coords(coords,tet10,t::Integer)
    C=Matrix{BigInt}(undef,3,10); exps=Vector{Int}(undef,3)
    @inbounds for d in 1:3
        parts=Vector{Tuple{BigInt,Int}}(undef,10); emin=typemax(Int)
        for k in 1:10
            x=coords[d,tet10[k,t]]; isfinite(x) ||
                throw(ArgumentError("P2 Jacobian: tet $t has a non-finite coordinate"))
            parts[k]=_dyadic_parts(x)
            parts[k][1]!=0 && (emin=min(emin,parts[k][2]))
        end
        emin==typemax(Int) && (emin=0)
        base=parts[1][1]==0 ? BigInt(0) : parts[1][1]<<(parts[1][2]-emin)
        tz=typemax(Int)
        for k in 1:10
            m,e=parts[k]
            q=(m==0 ? BigInt(0) : m<<(e-emin))-base
            C[d,k]=q; q!=0 && (tz=min(tz,trailing_zeros(q)))
        end
        if tz!=typemax(Int) && tz>0
            for k in 1:10; C[d,k] >>= tz; end
            emin += tz
        end
        exps[d]=emin
    end
    return C,exps[1]+exps[2]+exps[3]
end

@inline function _mixed_det(J,i,j,k)
    x1=J[1,i];x2=J[2,i];x3=J[3,i]
    y1=J[4,j];y2=J[5,j];y3=J[6,j]
    z1=J[7,k];z2=J[8,k];z3=J[9,k]
    x1*(y2*z3-y3*z2)-y1*(x2*z3-x3*z2)+z1*(x2*y3-x3*y2)
end

function _p2_bernstein_coeffs(coords,tet10,t::Integer)
    C,scaleexp=_integer_element_coords(coords,tet10,t)
    J=Matrix{BigInt}(undef,9,4)
    @inbounds for q in 1:4,c in 1:3,d in 1:3
        s=BigInt(0)
        for k in 1:10
            g=_REF_VERTEX_GRADS[q][k][c]
            g!=0 && (s += C[d,k]*g)
        end
        J[d+3(c-1),q]=s
    end
    sums=[BigInt(0) for _ in 1:20]
    @inbounds for i in 1:4,j in 1:4,k in 1:4
        g=Int(_B3_GROUP[i,j,k]); sums[g]+=_mixed_det(J,i,j,k)
    end
    return sums,scaleexp
end

function _p2_bernstein_bound(coords,tet10,t::Integer)
    sums,scaleexp=_p2_bernstein_coeffs(coords,tet10,t)
    mini=1
    @inbounds for i in 2:20
        sums[i]*_B3_NPERM[mini] < sums[mini]*_B3_NPERM[i] && (mini=i)
    end
    return all(>(0),sums),sums[mini],_B3_NPERM[mini],scaleexp
end

function _bound_float(n::BigInt,d::Int,e::Int)
    setprecision(BigFloat,128) do
        y=Float64(ldexp(BigFloat(n)/d,e))
        y==0.0 && n>0 && return nextfloat(0.0)
        y==0.0 && n<0 && return -nextfloat(0.0)
        y
    end
end

"""
    p2_min_jacobian(p::P2Mesh) -> Float64

Minimum exact degree-3 Bernstein coefficient of the P2 isoparametric Jacobian
determinant over all elements, returned as `Float64`. It is a conservative global
lower bound: `> 0` formally certifies every element over the whole reference
tetrahedron; `≤ 0` means the mesh is not certified (and may be folded). For a
straight-sided mesh this equals `6·min tet volume > 0`.
"""
function p2_min_jacobian(p::P2Mesh)
    ntets(p) == 0 && return 0.0
    mn = Inf
    for t in 1:ntets(p)
        _,n,dn,e = _p2_bernstein_bound(p.coords,p.tet10,t)
        d = _bound_float(n,dn,e)
        d < mn && (mn = d)
    end
    return mn
end

# ════════════════════════════════════════════════════════════════════════════════
# Curving edge nodes onto boundary geometry (safe: boundary edges only + validity)
# ════════════════════════════════════════════════════════════════════════════════
# mid-node → its two corner endpoints, and mid-node → incident tet list.
function _mid_maps(p::P2Mesh)
    endp = Dict{Int32, Tuple{Int32,Int32}}()
    midtets = Dict{Int32, Vector{Int32}}()
    slots = ((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,2,4),(10,3,4))
    @inbounds for t in 1:ntets(p), (sl,i,j) in slots
        m = p.tet10[sl,t]
        endp[m] = (p.tet10[i,t], p.tet10[j,t])
        push!(get!(() -> Int32[], midtets, m), Int32(t))
    end
    return endp, midtets
end

# undirected corner-edges that lie on a boundary surface face (of the linear tets).
function _boundary_corner_edges(p::P2Mesh)
    bf, _ = boundary_faces(@view p.tet10[1:4, :])
    edges = Set{Tuple{Int32,Int32}}()
    for f in bf
        a,b,c = f[1],f[2],f[3]
        push!(edges, minmax(a,b)); push!(edges, minmax(b,c)); push!(edges, minmax(a,c))
    end
    return edges
end

# bounding-box diagonal → relative displacement scale for the "actually moved" test.
function _bbox_diag(p::P2Mesh)
    lo1=lo2=lo3=Inf; hi1=hi2=hi3=-Inf
    @inbounds for i in 1:nnodes(p)
        x=p.coords[1,i]; y=p.coords[2,i]; z=p.coords[3,i]
        x<lo1 && (lo1=x); x>hi1 && (hi1=x)
        y<lo2 && (lo2=y); y>hi2 && (hi2=y)
        z<lo3 && (lo3=z); z>hi3 && (hi3=z)
    end
    d = hypot(hi1-lo1, hi2-lo2, hi3-lo3)
    isfinite(d) || throw(ArgumentError("P2 curving: mesh bounding-box diagonal is non-finite"))
    return d
end

# Core curving. `qualifies(a,b)::Bool` gates on the two corner endpoints; `project`
# snaps a point onto the surface. Only boundary-surface edges are curved, and any
# projection whose incident elements lack an exact global positive-Jacobian
# certificate is reverted (so the mesh never gains an inverted element). Deterministic (mid-nodes
# processed in ascending id order). Returns the number of nodes actually moved.
function _curve_boundary!(p::P2Mesh, qualifies, project; rtol::Float64=1e-6)
    ntets(p) == 0 && return 0
    @inbounds for t in 1:ntets(p)
        _p2_bernstein_bound(p.coords,p.tet10,t)[1] ||
            throw(ArgumentError("P2 curving: input tet $t lacks a global positive-Jacobian certificate"))
    end
    bnd = _boundary_corner_edges(p)
    endp, midtets = _mid_maps(p)
    movetol = rtol * _bbox_diag(p)
    isfinite(movetol) || throw(ArgumentError("P2 curving: displacement tolerance overflowed Float64"))
    ncurved = 0
    for mid in sort!(collect(keys(endp)))
        a, b = endp[mid]
        (minmax(a,b) in bnd) || continue        # boundary-surface edge only
        qualifies(a, b) || continue             # geometric predicate on endpoints
        @inbounds ox=p.coords[1,mid]; oy=p.coords[2,mid]; oz=p.coords[3,mid]
        raw = project(ox, oy, oz)
        q = try
            (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("P2 curving: project must return three Float64-representable coordinates"))
        end
        (isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3])) || continue
        @inbounds begin p.coords[1,mid]=q[1]; p.coords[2,mid]=q[2]; p.coords[3,mid]=q[3] end
        ok = true
        for t in midtets[mid]
            certified,_,_,_ = _p2_bernstein_bound(p.coords,p.tet10,t)
            certified || (ok=false;break)
        end
        if ok
            hypot(q[1]-ox, q[2]-oy, q[3]-oz) > movetol && (ncurved += 1)
        else
            @inbounds begin p.coords[1,mid]=ox; p.coords[2,mid]=oy; p.coords[3,mid]=oz end
        end
    end
    return ncurved
end

"""
    curve_to_cylinder!(p::P2Mesh, center, axis, radius; rtol=1e-6) -> n_curved

Curve the quadratic mesh onto a cylinder of the given `center`, `axis`, `radius`:
every **boundary-surface** edge whose *both* corner endpoints lie on the lateral
wall (radial distance `radius` from the axis, within `rtol`) has its mid-node
projected radially onto the exact cylinder — turning the straight chord into a
true curved P2 edge. Interior chords are never touched, and any projection that
would invert an incident element is reverted (see module docstring). Flat caps and
axial spokes stay straight. Returns the number of mid-nodes moved.
"""
function curve_to_cylinder!(p::P2Mesh, center, axis, radius::Real; rtol::Real=1e-6)
    R = _input_float(radius, "curve_to_cylinder!: radius")
    c = _input_point3(center, "curve_to_cylinder!: center")
    (isfinite(R) && R > 0) ||
        throw(ArgumentError("curve_to_cylinder!: radius must be finite and positive (got $radius)"))
    (isfinite(c[1]) && isfinite(c[2]) && isfinite(c[3])) ||
        throw(ArgumentError("curve_to_cylinder!: center must be finite"))
    tolrel = _input_float(rtol, "curve_to_cylinder!: rtol")
    (isfinite(tolrel) && tolrel >= 0) ||
        throw(ArgumentError("curve_to_cylinder!: rtol must be finite and non-negative (got $rtol)"))
    ez = _unit3(_input_point3(axis, "curve_to_cylinder!: axis"))
    tol = tolrel * R
    isfinite(tol) || throw(ArgumentError("curve_to_cylinder!: rtol * radius overflowed Float64"))
    @inline function radial(x, y, z)
        rx=x-c[1]; ry=y-c[2]; rz=z-c[3]
        (isfinite(rx) && isfinite(ry) && isfinite(rz)) ||
            throw(ArgumentError("curve_to_cylinder!: coordinate subtraction overflowed Float64"))
        ax = rx*ez[1] + ry*ez[2] + rz*ez[3]
        vx=rx-ax*ez[1]; vy=ry-ax*ez[2]; vz=rz-ax*ez[3]
        d=hypot(vx,vy,vz)
        (isfinite(ax) && isfinite(vx) && isfinite(vy) && isfinite(vz) && isfinite(d)) ||
            throw(ArgumentError("curve_to_cylinder!: radial projection overflowed Float64"))
        (d, ax, vx, vy, vz)
    end
    onwall(v::Integer) = abs(radial(p.coords[1,v], p.coords[2,v], p.coords[3,v])[1] - R) <= tol
    qualifies(a::Integer, b::Integer) = onwall(a) && onwall(b)
    function project(x, y, z)
        d, ax, vx, vy, vz = radial(x, y, z)
        d == 0 && return (Inf, Inf, Inf)          # on the axis — unprojectable, skip
        f = R/d
        (c[1] + ax*ez[1] + vx*f, c[2] + ax*ez[2] + vy*f, c[3] + ax*ez[3] + vz*f)
    end
    return _curve_boundary!(p, qualifies, project; rtol=tolrel)
end

@inline function _unit3(a)
    l = hypot(a[1],a[2],a[3])
    (isfinite(l) && l > 0) ||
        throw(ArgumentError("curve_to_cylinder!: axis must have finite positive length"))
    (a[1]/l, a[2]/l, a[3]/l)
end

function _input_float(x::Real, what::AbstractString)
    y = try
        Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$what must be representable as Float64"))
    end
    return y
end

function _input_point3(p, what::AbstractString)
    try
        return (_input_float(p[1], what), _input_float(p[2], what), _input_float(p[3], what))
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("$what must contain three real coordinates"))
    end
end

"""
    curve_to_surface!(p::P2Mesh, project, on_surface; rtol=1e-6) -> n_curved

Curve the quadratic mesh onto an **arbitrary** target surface — the general form
of [`curve_to_cylinder!`](@ref). `project(x,y,z) -> (x,y,z)` snaps a point onto the
surface; `on_surface(x,y,z) -> Bool` tests whether a *corner* lies on it. A
mid-node is curved only when it sits on a genuine **boundary-surface edge** *and*
both its corner endpoints satisfy `on_surface` — interior chords (whose endpoints
may also lie on the surface) are never moved. Any projection that would invert an
incident element, or that returns a non-finite point, is reverted, so the result
never contains an inverted element. Returns the number of mid-nodes actually moved
(displacement > `rtol` × bounding-box diagonal).
"""
function curve_to_surface!(p::P2Mesh, project, on_surface; rtol::Real=1e-6)
    tolrel = _input_float(rtol, "curve_to_surface!: rtol")
    (isfinite(tolrel) && tolrel >= 0) ||
        throw(ArgumentError("curve_to_surface!: rtol must be finite and non-negative (got $rtol)"))
    function onsurf(v::Integer)
        q = @inbounds on_surface(p.coords[1,v], p.coords[2,v], p.coords[3,v])
        q isa Bool || throw(ArgumentError("curve_to_surface!: on_surface must return Bool"))
        return q
    end
    qualifies(a::Integer, b::Integer) = onsurf(a) && onsurf(b)
    return _curve_boundary!(p, qualifies, (x,y,z) -> project(x,y,z); rtol=tolrel)
end

"""
    write_msh_p2(path, p::P2Mesh; tet_tag=zeros) -> path

Write a quadratic tet mesh as gmsh MSH v2.2 with 10-node (type-11) elements.
"""
function write_msh_p2(path::AbstractString, p::P2Mesh;
                      tet_tag::AbstractVector{<:Integer}=zeros(Int32, ntets(p)))
    length(tet_tag) == ntets(p) || throw(ArgumentError("write_msh_p2: tet_tag length mismatch"))
    tags=Vector{Int32}(undef,ntets(p))
    @inbounds for t in 1:ntets(p)
        tag=tet_tag[t]
        0<=tag<=typemax(Int32) ||
            throw(ArgumentError("write_msh_p2: tet tag $tag must be non-negative and fit Int32"))
        tags[t]=Int32(tag)
    end
    # Refuse folded input before touching the destination.
    @inbounds for t in 1:ntets(p)
        ok,_,_,_=_p2_bernstein_bound(p.coords,p.tet10,t)
        ok || throw(ArgumentError("write_msh_p2: tet $t lacks a global positive-Jacobian certificate"))
    end
    target=abspath(path);parent=dirname(target)
    isdir(parent) || throw(ArgumentError("write_msh_p2: parent directory does not exist: $parent"))
    mktemp(parent) do tmp,io
        println(io, "\$MeshFormat"); println(io, "2.2 0 8"); println(io, "\$EndMeshFormat")
        println(io, "\$Nodes"); println(io, nnodes(p))
        @inbounds for i in 1:nnodes(p)
            @printf(io, "%d %.17g %.17g %.17g\n", i, p.coords[1,i], p.coords[2,i], p.coords[3,i])
        end
        println(io, "\$EndNodes")
        println(io, "\$Elements"); println(io, ntets(p))
        @inbounds for t in 1:ntets(p)
            tag = tags[t]
            print(io, t, " 11 2 ", tag, " ", tag)
            for k in 1:10; print(io, " ", p.tet10[k,t]); end
            println(io)
        end
        println(io, "\$EndElements")
        flush(io);close(io);mv(tmp,target;force=true)
    end
    return path
end

end # module HighOrder
