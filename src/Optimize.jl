"""
    Optimize

Stage-4 mesh optimization and quality reporting (PLAN.md §3 "Optimize"). Operates
on finalized tetrahedral [`Mesh`](@ref)es:

* [`mesh_quality`](@ref) — dihedral-angle and radius-edge statistics, minimum tet
  volume, and a sliver count (the honest quality report a solver needs).
* [`smooth_laplacian`](@ref) — boundary-preserving Laplacian smoothing of interior
  nodes, with a positive-volume guard so a move never inverts an incident tet.

Both preserve total volume (boundary fixed) and mesh validity. Smoothing improves
the *mean* dihedral and typically reduces the sliver count; it is **not** min-angle
optimal (Laplacian relaxation can leave or slightly worsen an individual worst
element — genuine sliver removal needs targeted flips/exudation, tracked as
remaining Stage-4 work). All of these properties are verified in the tests.
"""
module Optimize

using ..MeshTypes: Mesh, nnodes, ntets, node, tet_signed_volume, tet_dihedral_extrema,
                   tet_radius_edge, tet_volume, boundary_faces, validate

export mesh_quality, smooth_laplacian, TetQuality

struct TetQuality
    n_tets::Int
    min_dihedral_deg::Float64
    mean_dihedral_deg::Float64
    max_dihedral_deg::Float64
    min_radius_edge::Float64
    mean_radius_edge::Float64
    min_volume::Float64
    n_slivers::Int              # min dihedral < 10° or max dihedral > 170°
end

function Base.show(io::IO, q::TetQuality)
    print(io, "TetQuality(tets=", q.n_tets,
          ", dihedral[min,mean,max]°=(", round(q.min_dihedral_deg;digits=2), ",",
          round(q.mean_dihedral_deg;digits=2), ",", round(q.max_dihedral_deg;digits=2), ")",
          ", radius_edge[min,mean]=(", round(q.min_radius_edge;digits=3), ",",
          round(q.mean_radius_edge;digits=3), ")",
          ", min_vol=", round(q.min_volume;sigdigits=3),
          ", slivers=", q.n_slivers, ")")
end

"""
    mesh_quality(m; sliver_deg=10.0) -> TetQuality

Quality report over a tet mesh: extremal/mean dihedral angles (degrees), radius-
edge ratios, minimum tet volume, and the count of slivers (min dihedral <
`sliver_deg` or max dihedral > `180 − sliver_deg`).
"""
function mesh_quality(m::Mesh; sliver_deg::Real=10.0)
    nt = ntets(m)
    nt == 0 && return TetQuality(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)
    dmin = Inf; dmax = -Inf; dsum = 0.0
    remin = Inf; resum = 0.0
    vmin = Inf; nsliv = 0
    lo = deg2rad(sliver_deg); hi = deg2rad(180 - sliver_deg)
    @inbounds for t in 1:nt
        a=node(m,m.tets[1,t]); b=node(m,m.tets[2,t]); c=node(m,m.tets[3,t]); d=node(m,m.tets[4,t])
        mn, mx = tet_dihedral_extrema(a,b,c,d)
        dmin = min(dmin, mn); dmax = max(dmax, mx); dsum += mn
        (mn < lo || mx > hi) && (nsliv += 1)
        re = tet_radius_edge(a,b,c,d); isfinite(re) && (remin = min(remin, re); resum += re)
        vmin = min(vmin, tet_volume(a,b,c,d))
    end
    return TetQuality(nt, rad2deg(dmin), rad2deg(dsum/nt), rad2deg(dmax),
                      remin, resum/nt, vmin, nsliv)
end

"""
    smooth_laplacian(m; iters=5, relax=1.0) -> Mesh

Boundary-preserving Laplacian smoothing: each interior node is moved toward the
average of its edge-neighbours by factor `relax`, but only if the move keeps every
incident tet positively-oriented (otherwise it is skipped, so validity and total
volume are preserved). Improves the mean dihedral / reduces slivers; not min-angle
optimal. Returns a new smoothed `Mesh`.
"""
function smooth_laplacian(m::Mesh; iters::Integer=5, relax::Real=1.0)
    nn = nnodes(m); nt = ntets(m)
    coords = copy(m.coords)
    isboundary = _boundary_nodes(m)
    neigh = _node_neighbours(m)
    inc = _node_tets(m)
    @inbounds for _ in 1:iters
        for v in 1:nn
            isboundary[v] && continue
            isempty(neigh[v]) && continue
            # target = neighbour centroid
            sx=0.0; sy=0.0; sz=0.0
            for w in neigh[v]; sx+=coords[1,w]; sy+=coords[2,w]; sz+=coords[3,w]; end
            k = length(neigh[v])
            tx = coords[1,v] + relax*(sx/k - coords[1,v])
            ty = coords[2,v] + relax*(sy/k - coords[2,v])
            tz = coords[3,v] + relax*(sz/k - coords[3,v])
            # accept only if no incident tet inverts (positive signed volume kept)
            ox=coords[1,v]; oy=coords[2,v]; oz=coords[3,v]
            coords[1,v]=tx; coords[2,v]=ty; coords[3,v]=tz
            ok = true
            for t in inc[v]
                a=m.tets[1,t];b=m.tets[2,t];c=m.tets[3,t];d=m.tets[4,t]
                pa=(coords[1,a],coords[2,a],coords[3,a]); pb=(coords[1,b],coords[2,b],coords[3,b])
                pc=(coords[1,c],coords[2,c],coords[3,c]); pd=(coords[1,d],coords[2,d],coords[3,d])
                if tet_signed_volume(pa,pb,pc,pd) <= 0; ok=false; break; end
            end
            ok || (coords[1,v]=ox; coords[2,v]=oy; coords[3,v]=oz)
        end
    end
    return Mesh(coords; tets=copy(m.tets), tet_tag=copy(m.tet_tag))
end

# nodes lying on any boundary face (fixed during smoothing)
function _boundary_nodes(m::Mesh)
    isb = falses(nnodes(m))
    bf, _ = boundary_faces(m.tets)
    @inbounds for f in bf; isb[f[1]]=true; isb[f[2]]=true; isb[f[3]]=true; end
    return isb
end

function _node_neighbours(m::Mesh)
    nb = [Set{Int32}() for _ in 1:nnodes(m)]
    @inbounds for t in 1:ntets(m)
        v = (m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
        for i in 1:4, j in 1:4
            i==j && continue
            push!(nb[v[i]], v[j])
        end
    end
    return [collect(s) for s in nb]
end

function _node_tets(m::Mesh)
    inc = [Int32[] for _ in 1:nnodes(m)]
    @inbounds for t in 1:ntets(m), i in 1:4
        push!(inc[m.tets[i,t]], Int32(t))
    end
    return inc
end

end # module Optimize
