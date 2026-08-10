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

using ..MeshTypes: Mesh, nnodes, ntris, node, triangle_area, boundary_edges

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
    msgs = String[]
    nt = ntris(m)
    # edge incidence with directed counts (for orientation) and undirected (manifold)
    dircount = Dict{Tuple{Int32,Int32},Int}()      # directed a→b occurrences
    undcount = Dict{Tuple{Int32,Int32},Int}()      # undirected {a,b} occurrences
    facekeys = Dict{Tuple{Int32,Int32,Int32},Int}()
    ndeg = 0
    minedge = Inf
    @inbounds for t in 1:nt
        a=m.tris[1,t]; b=m.tris[2,t]; c=m.tris[3,t]
        if triangle_area(node(m,a),node(m,b),node(m,c)) == 0
            ndeg += 1
        end
        facekeys[_f3(a,b,c)] = get(facekeys, _f3(a,b,c), 0) + 1
        for (u,v) in ((a,b),(b,c),(c,a))
            dircount[(Int32(u),Int32(v))] = get(dircount, (Int32(u),Int32(v)), 0) + 1
            undcount[_e2(u,v)] = get(undcount, _e2(u,v), 0) + 1
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
    ncoin = _count_coincident(m, tol)

    closed = nopen == 0
    manifold = nnon == 0
    oriented = noriented_bad == 0
    closed || push!(msgs, "$nopen open (boundary) edges — surface is not closed")
    manifold || push!(msgs, "$nnon non-manifold edges (shared by >2 triangles)")
    (closed && manifold && !oriented) && push!(msgs, "$noriented_bad edges with inconsistent triangle orientation")
    ndeg > 0 && push!(msgs, "$ndeg degenerate (zero-area) triangles")
    ndup > 0 && push!(msgs, "$ndup duplicate triangles")
    ncoin > 0 && push!(msgs, "$ncoin near-coincident vertex pairs (< tol) — should be merged")
    return SurfaceReport(closed, manifold, oriented, nopen, nnon, ndeg, ndup, ncoin,
                         minedge == Inf ? 0.0 : minedge, msgs)
end

@inline _dist(a,b) = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)

# count distinct vertex pairs closer than tol·diag (via a hash grid, O(n))
function _count_coincident(m::Mesh, reltol::Real)
    nn = nnodes(m); nn == 0 && return 0
    lo=(Inf,Inf,Inf); hi=(-Inf,-Inf,-Inf)
    @inbounds for i in 1:nn
        p=node(m,i); lo=(min(lo[1],p[1]),min(lo[2],p[2]),min(lo[3],p[3]))
        hi=(max(hi[1],p[1]),max(hi[2],p[2]),max(hi[3],p[3]))
    end
    diag = sqrt((hi[1]-lo[1])^2+(hi[2]-lo[2])^2+(hi[3]-lo[3])^2)
    tol = max(diag*reltol, eps(Float64))
    inv = 1.0/tol
    cell = Dict{NTuple{3,Int}, Vector{Int}}()
    key(p) = (round(Int,p[1]*inv), round(Int,p[2]*inv), round(Int,p[3]*inv))
    cnt = 0
    @inbounds for i in 1:nn
        p = node(m,i); k = key(p)
        # check the 27 neighbouring cells for an already-seen coincident vertex
        for dx in -1:1, dy in -1:1, dz in -1:1
            nb = (k[1]+dx, k[2]+dy, k[3]+dz)
            haskey(cell, nb) || continue
            for j in cell[nb]
                _dist(node(m,i), node(m,j)) < tol && (cnt += 1)
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
         r.n_duplicate_tris == 0 && r.n_coincident_pairs == 0
    return ok, r
end

end # module Heal
