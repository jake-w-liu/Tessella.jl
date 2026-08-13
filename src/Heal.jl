"""
    Heal

Stage-5 "heal, don't fail" — the *detection* half (PLAN.md §3 "Heal", principle
3). Before a boundary surface is handed to the 3-D mesher it must be a valid
closed manifold; a defect (open edge, non-manifold edge, degenerate/duplicate
triangle, near-coincident vertices, inconsistent orientation) is exactly what
produces "invalid boundary mesh / overlapping facets" downstream. This module
reports those defects precisely so a run either proceeds on a clean surface or
returns an explicit blocker — never a silent bad mesh.

[`surface_diagnostics`](@ref) returns a [`SurfaceReport`](@ref); [`is_meshable`](@ref)
is the go/no-go gate used by the volume mesher.
"""
module Heal

using ..MeshTypes: Mesh, nnodes, ntris, node, triangle_area, boundary_edges, _tri_topology

export SurfaceReport, surface_diagnostics, is_meshable

struct SurfaceReport
    closed::Bool                    # no open (boundary) edges
    manifold::Bool                  # every edge shared by ≤ 2 triangles
    oriented::Bool                  # consistent outward orientation (shared edges reversed)
    n_open_edges::Int
    n_nonmanifold_edges::Int
    n_degenerate_tris::Int          # zero-area
    n_duplicate_tris::Int           # same vertex set
    n_coincident_pairs::Int         # distinct vertices within merge tolerance
    min_edge_length::Float64
    messages::Vector{String}
end

function Base.show(io::IO, r::SurfaceReport)
    print(io, "SurfaceReport(closed=", r.closed, ", manifold=", r.manifold,
          ", oriented=", r.oriented, ", open_edges=", r.n_open_edges,
          ", nonmanifold=", r.n_nonmanifold_edges, ", degenerate=", r.n_degenerate_tris,
          ", duplicate=", r.n_duplicate_tris, ", coincident=", r.n_coincident_pairs,
          ", min_edge=", round(r.min_edge_length; sigdigits=3), ")")
end

@inline _e2(a,b) = a <= b ? (Int32(a),Int32(b)) : (Int32(b),Int32(a))
@inline _f3(a,b,c) = begin
    x=Int32(a);y=Int32(b);z=Int32(c)
    x>y && ((x,y)=(y,x)); y>z && ((y,z)=(z,y)); x>y && ((x,y)=(y,x)); (x,y,z)
end

"""
    surface_diagnostics(m; tol=1e-9) -> SurfaceReport

Analyse a triangle surface mesh for the defects that break volume meshing.
`tol` is the coincident-vertex distance relative to the bounding-box diagonal.
"""
function surface_diagnostics(m::Mesh; tol::Real=1e-9)
    reltol = try
        Float64(tol)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("surface_diagnostics: tol must be a finite non-negative number (got $tol)"))
    end
    (isfinite(reltol) && reltol >= 0.0) ||
        throw(ArgumentError("surface_diagnostics: tol must be finite and non-negative (got $tol)"))

    msgs = String[]
    nt = ntris(m)
    nnonfinite = 0
    @inbounds for i in 1:nnodes(m)
        p = node(m, i)
        (isfinite(p[1]) && isfinite(p[2]) && isfinite(p[3])) || (nnonfinite += 1)
    end
    # edge incidence with directed counts (for orientation) and undirected (manifold)
    dircount = Dict{Tuple{Int32,Int32},Int}()      # directed a→b occurrences
    undcount = Dict{Tuple{Int32,Int32},Int}()      # undirected {a,b} occurrences
    facekeys = Dict{Tuple{Int32,Int32,Int32},Int}()
    ndeg = 0
    minedge = Inf; hasedge = false
    @inbounds for t in 1:nt
        a=m.tris[1,t]; b=m.tris[2,t]; c=m.tris[3,t]
        area = triangle_area(node(m,a),node(m,b),node(m,c))
        if !(isfinite(area) && area > 0.0)
            ndeg += 1
        end
        facekeys[_f3(a,b,c)] = get(facekeys, _f3(a,b,c), 0) + 1
        for (u,v) in ((a,b),(b,c),(c,a))
            dircount[(Int32(u),Int32(v))] = get(dircount, (Int32(u),Int32(v)), 0) + 1
            undcount[_e2(u,v)] = get(undcount, _e2(u,v), 0) + 1
            hasedge = true
            l = _dist(node(m,u), node(m,v)); l < minedge && (minedge = l)
        end
    end
    # open / non-manifold / orientation
    nopen = 0; nnon = 0; noriented_bad = 0
    for (e, c) in undcount
        c == 1 && (nopen += 1)
        c > 2 && (nnon += 1)
        # consistent orientation: for a 2-shared edge, the two directed uses are
        # opposite (a→b once and b→a once) ⇒ dircount[a→b]==1 and [b→a]==1.
        if c == 2
            da = get(dircount, (e[1],e[2]), 0); db = get(dircount, (e[2],e[1]), 0)
            (da == 1 && db == 1) || (noriented_bad += 1)
        end
    end
    ndup = sum(v-1 for v in values(facekeys) if v > 1; init=0)
    ncoin = nnonfinite == 0 ? _count_coincident(m, reltol) : 0

    topok,topreason=_tri_topology(m.tris)
    closed = nt > 0 && nopen == 0
    manifold = nnon == 0 && topok
    oriented = noriented_bad == 0
    nt == 0 && push!(msgs, "surface has no triangles")
    nt > 0 && !closed && push!(msgs, "$nopen open (boundary) edges — surface is not closed")
    nnon==0 || push!(msgs, "$nnon non-manifold edges (shared by >2 triangles)")
    topok || push!(msgs,"non-manifold triangle topology: $topreason")
    (closed && manifold && !oriented) && push!(msgs, "$noriented_bad edges with inconsistent triangle orientation")
    ndeg > 0 && push!(msgs, "$ndeg degenerate or non-finite-area triangles")
    ndup > 0 && push!(msgs, "$ndup duplicate triangles")
    ncoin > 0 && push!(msgs, "$ncoin near-coincident vertex pairs (< tol) — should be merged")
    nnonfinite > 0 && push!(msgs, "$nnonfinite vertices have non-finite coordinates")
    return SurfaceReport(closed, manifold, oriented, nopen, nnon, ndeg, ndup, ncoin,
                         hasedge ? minedge : 0.0, msgs)
end

@inline function _dist(a, b)
    # Direct subtraction preserves local precision for far-origin geometry.  When
    # opposite-sign coordinates make that subtraction overflow, normalize first;
    # a mathematically unrepresentable distance then becomes Inf only on rescaling.
    dx = a[1]-b[1]; dy = a[2]-b[2]; dz = a[3]-b[3]
    # Same-scale, same-sign points (including far-origin local geometry) take this
    # more accurate path; `hypot` scales internally and cannot overflow while the
    # coordinate differences themselves remain finite.
    (isfinite(dx) && isfinite(dy) && isfinite(dz)) && return hypot(dx, dy, dz)
    scale = max(abs(a[1]), abs(a[2]), abs(a[3]),
                abs(b[1]), abs(b[2]), abs(b[3]))
    scale == 0.0 && return 0.0
    return scale * hypot(a[1]/scale - b[1]/scale,
                         a[2]/scale - b[2]/scale,
                         a[3]/scale - b[3]/scale)
end

# count distinct vertex pairs closer than tol·diag (via a hash grid, O(n))
function _count_coincident(m::Mesh, reltol::Float64)
    nn = nnodes(m); nn == 0 && return 0
    # Work in a globally normalized, bbox-relative coordinate system.  Hashing
    # absolute `p/tol` overflows Int for ordinary geometry translated far from the
    # origin; normalization bounds every cell coordinate by roughly 2/eps.
    scale = 0.0
    @inbounds for i in 1:nn
        p = node(m, i)
        scale = max(scale, abs(p[1]), abs(p[2]), abs(p[3]))
    end
    scale == 0.0 && (scale = 1.0)
    lo=(Inf,Inf,Inf); hi=(-Inf,-Inf,-Inf)
    @inbounds for i in 1:nn
        p=node(m,i); q=(p[1]/scale,p[2]/scale,p[3]/scale)
        lo=(min(lo[1],q[1]),min(lo[2],q[2]),min(lo[3],q[3]))
        hi=(max(hi[1],q[1]),max(hi[2],q[2]),max(hi[3],q[3]))
    end
    diag = hypot(hi[1]-lo[1], hi[2]-lo[2], hi[3]-lo[3])
    cellsize = diag*reltol
    if cellsize==0
        counts=Dict{NTuple{3,Float64},Int}()
        pairs=0
        @inbounds for i in 1:nn
            p=node(m,i);key=(p[1]==0 ? 0.0 : p[1],p[2]==0 ? 0.0 : p[2],p[3]==0 ? 0.0 : p[3])
            old=get(counts,key,0);pairs+=old;counts[key]=old+1
        end
        return pairs
    end
    inv = 1.0/cellsize
    if !isfinite(inv) || diag*inv>typemax(Int)-2
        order=sortperm(1:nn;by=i->node(m,i)[1]/scale)
        cnt=0
        @inbounds for ii in eachindex(order)
            i=order[ii];qi=(node(m,i)[1]/scale,node(m,i)[2]/scale,node(m,i)[3]/scale)
            jj=ii-1
            while jj>=1
                j=order[jj];qj=(node(m,j)[1]/scale,node(m,j)[2]/scale,node(m,j)[3]/scale)
                qi[1]-qj[1] < cellsize || break
                hypot(qi[1]-qj[1],qi[2]-qj[2],qi[3]-qj[3]) < cellsize && (cnt+=1)
                jj-=1
            end
        end
        return cnt
    end
    cell = Dict{NTuple{3,Int}, Vector{Int}}()
    sizehint!(cell, nn)
    qcoord(p) = (p[1]/scale-lo[1], p[2]/scale-lo[2], p[3]/scale-lo[3])
    key(q) = (floor(Int,q[1]*inv), floor(Int,q[2]*inv), floor(Int,q[3]*inv))
    cnt = 0
    @inbounds for i in 1:nn
        qi = qcoord(node(m,i)); k = key(qi)
        # check the 27 neighbouring cells for an already-seen coincident vertex
        for dx in -1:1, dy in -1:1, dz in -1:1
            nb = (k[1]+dx, k[2]+dy, k[3]+dz)
            haskey(cell, nb) || continue
            for j in cell[nb]
                qj = qcoord(node(m,j))
                hypot(qi[1]-qj[1], qi[2]-qj[2], qi[3]-qj[3]) < cellsize && (cnt += 1)
            end
        end
        push!(get!(cell, k, Int[]), i)
    end
    return cnt
end

"""
    is_meshable(m; tol=1e-9) -> (ok::Bool, report::SurfaceReport)

Go/no-go gate for volume meshing: the surface must be closed, manifold, oriented,
and free of degenerate/duplicate/coincident defects.
"""
function is_meshable(m::Mesh; tol::Real=1e-9)
    r = surface_diagnostics(m; tol=tol)
    ok = r.closed && r.manifold && r.oriented && r.n_degenerate_tris == 0 &&
         r.n_duplicate_tris == 0 && r.n_coincident_pairs == 0 &&
         all(isfinite, m.coords)
    return ok, r
end

end # module Heal
