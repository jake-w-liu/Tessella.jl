"""
    Optimize

Stage-4 mesh optimization and quality reporting (PLAN.md §3 "Optimize"). Operates
on finalized tetrahedral [`Mesh`](@ref)es:

* [`mesh_quality`](@ref) — dihedral-angle and radius-edge statistics, minimum tet
  volume, and a sliver count (the honest quality report a solver needs).
* [`smooth_laplacian`](@ref) — boundary-preserving Laplacian smoothing of interior
  nodes, with a positive-volume guard so a move never inverts an incident tet.
* [`smooth_odt`](@ref) — boundary-preserving optimal-Delaunay-triangulation (ODT)
  smoothing: interior nodes move to the volume-weighted average of the circumcenters
  of their incident tets, under the same positive-volume guard.
* [`smooth_optimize`](@ref) — boundary-preserving **optimization-based** smoothing
  that maximizes each poor vertex's *worst* local dihedral (a pattern search), so it
  removes needle **and** cap slivers the mean-smoothers leave behind — the geometric
  half of sliver removal, complementing the topological flips in `Mesh3D`.

All three smoothers preserve total volume (boundary fixed) and mesh validity. The
first two lift the *mean* dihedral (relaxation, not min-angle optimal);
`smooth_optimize` targets the *min* angle directly. Slivers whose four vertices are
topology-locked still need [`Mesh3D.optimize_flips!`](@ref) (2-3/3-2 flips) — the two
together are the sliver-removal toolkit; the `mesh_volume(optimize=true)` pipeline
runs flips then optimization smoothing. These properties are verified in the tests.
"""
module Optimize

using ..MeshTypes: Mesh, nnodes, ntets, node, tet_signed_volume, tet_dihedral_extrema,
                   tet_radius_edge, tet_volume, boundary_faces, validate

export mesh_quality, smooth_laplacian, smooth_odt, smooth_optimize, remove_slivers, TetQuality

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
    threshold = float(sliver_deg)
    (isfinite(threshold) && 0 <= threshold < 90) ||
        throw(ArgumentError("mesh_quality: sliver_deg must be finite and in [0, 90) (got $sliver_deg)"))
    nt = ntets(m)
    nt == 0 && return TetQuality(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)
    dmin = Inf; dmax = -Inf; dsum = 0.0
    remin = Inf; resum = 0.0
    vmin = Inf; nsliv = 0
    lo = deg2rad(threshold); hi = deg2rad(180 - threshold)
    @inbounds for t in 1:nt
        a=node(m,m.tets[1,t]); b=node(m,m.tets[2,t]); c=node(m,m.tets[3,t]); d=node(m,m.tets[4,t])
        mn, mx = tet_dihedral_extrema(a,b,c,d)
        dmin = min(dmin, mn); dmax = max(dmax, mx); dsum += mn
        (mn < lo || mx > hi) && (nsliv += 1)
        re = tet_radius_edge(a,b,c,d)
        (isfinite(mn) && isfinite(mx) && isfinite(re)) ||
            throw(ArgumentError("mesh_quality: tet $t has non-finite computed quality; validate the mesh first"))
        remin = min(remin, re); resum += re
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
    iters >= 0 || throw(ArgumentError("smooth_laplacian: iters must be non-negative (got $iters)"))
    r = float(relax)
    (isfinite(r) && 0 <= r <= 1) ||
        throw(ArgumentError("smooth_laplacian: relax must be finite and in [0, 1] (got $relax)"))
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
            tx = coords[1,v] + r*(sx/k - coords[1,v])
            ty = coords[2,v] + r*(sy/k - coords[2,v])
            tz = coords[3,v] + r*(sz/k - coords[3,v])
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
    return Mesh(coords; segs=copy(m.segs), tris=copy(m.tris), tets=copy(m.tets),
                seg_tag=copy(m.seg_tag), tri_tag=copy(m.tri_tag), tet_tag=copy(m.tet_tag))
end

"""
    smooth_odt(m; iters=5) -> Mesh

Boundary-preserving optimal-Delaunay-triangulation (ODT) smoothing. Each interior
node is moved to the **volume-weighted average of the circumcenters** of its
incident tets (the ODT node update), but only if the move keeps every incident tet
positively-oriented (otherwise it is skipped, so validity and total volume are
preserved). Boundary nodes (on any boundary face) stay fixed; near-zero-volume
incident tets are skipped when averaging (their circumcenter is ill-conditioned).
Like [`smooth_laplacian`](@ref) it lifts the *mean* dihedral and is not min-angle
optimal (targeted flips/exudation remain the fix for an individual worst sliver);
it simply uses a different node target. Returns a new smoothed `Mesh`.
"""
function smooth_odt(m::Mesh; iters::Integer=5)
    iters >= 0 || throw(ArgumentError("smooth_odt: iters must be non-negative (got $iters)"))
    nn = nnodes(m)
    coords = copy(m.coords)
    isboundary = _boundary_nodes(m)
    inc = _node_tets(m)
    @inbounds for _ in 1:iters
        for v in 1:nn
            isboundary[v] && continue
            isempty(inc[v]) && continue
            # ODT target: circumcenters of incident tets weighted by their volume
            sx=0.0; sy=0.0; sz=0.0; wsum=0.0
            for t in inc[v]
                a=m.tets[1,t];b=m.tets[2,t];c=m.tets[3,t];d=m.tets[4,t]
                pa=(coords[1,a],coords[2,a],coords[3,a]); pb=(coords[1,b],coords[2,b],coords[3,b])
                pc=(coords[1,c],coords[2,c],coords[3,c]); pd=(coords[1,d],coords[2,d],coords[3,d])
                vol = tet_volume(pa,pb,pc,pd)
                vol <= _odt_voleps(pa,pb,pc,pd) && continue    # skip degenerate tet
                cc = _tet_circumcenter(pa,pb,pc,pd)
                (isfinite(cc[1]) && isfinite(cc[2]) && isfinite(cc[3])) || continue
                sx += vol*cc[1]; sy += vol*cc[2]; sz += vol*cc[3]; wsum += vol
            end
            wsum > 0 || continue                               # no usable incident tet
            tx = sx/wsum; ty = sy/wsum; tz = sz/wsum
            # A symmetric star can reproduce the current point up to a few rounding
            # ulps.  Treat that as an exact no-op so smoothing does not create
            # meaningless coordinate drift (and remains bit-stable on fixed points).
            ox=coords[1,v]; oy=coords[2,v]; oz=coords[3,v]
            move = max(abs(tx-ox), abs(ty-oy), abs(tz-oz))
            movescl = max(abs(ox), abs(oy), abs(oz), abs(tx), abs(ty), abs(tz), 1.0)
            move <= 32eps(Float64)*movescl && continue
            # accept only if no incident tet inverts (positive signed volume kept)
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
    return Mesh(coords; segs=copy(m.segs), tris=copy(m.tris), tets=copy(m.tets),
                seg_tag=copy(m.seg_tag), tri_tag=copy(m.tri_tag), tet_tag=copy(m.tet_tag))
end

"""
    smooth_optimize(m; iters=8, sliver_deg=10.0) -> Mesh

Boundary-preserving **optimization-based** smoothing that targets slivers directly.
Unlike [`smooth_laplacian`](@ref)/[`smooth_odt`](@ref) (which relax toward a
geometric average and lift the *mean* dihedral), this maximizes the **worst** local
dihedral quality of each poor interior vertex's star via a pattern search: the
objective per star is `min over incident tets of min(min_dihedral, π − max_dihedral)`,
so it penalizes both needle (tiny min-angle) and cap (near-flat, huge max-angle)
slivers. Only vertices whose star quality is below `sliver_deg` are moved (targeted),
and a move is accepted only if it *improves* that objective while keeping every
incident tet positively oriented — so validity and total volume are preserved and the
per-star worst angle is monotone non-worsening. Complements the topological
[`Mesh3D.optimize_flips!`](@ref) (flips change connectivity; this repositions nodes).
Returns a new `Mesh`.
"""
function smooth_optimize(m::Mesh; iters::Integer=8, sliver_deg::Real=10.0)
    iters >= 0 || throw(ArgumentError("smooth_optimize: iters must be non-negative (got $iters)"))
    threshold = float(sliver_deg)
    (isfinite(threshold) && 0 <= threshold < 90) ||
        throw(ArgumentError("smooth_optimize: sliver_deg must be finite and in [0, 90) (got $sliver_deg)"))
    nn = nnodes(m)
    coords = copy(m.coords)
    isboundary = _boundary_nodes(m)
    inc = _node_tets(m)
    thr = deg2rad(threshold)
    dirs = ((1.,0.,0.),(-1.,0.,0.),(0.,1.,0.),(0.,-1.,0.),(0.,0.,1.),(0.,0.,-1.))
    @inbounds for _ in 1:iters
        for v in 1:nn
            (isboundary[v] || isempty(inc[v])) && continue
            tv = inc[v]
            q0 = _star_quality(m, coords, tv)
            q0 >= thr && continue                          # targeted: skip already-good stars
            scale = _star_scale(m, coords, tv)
            scale > 0 || continue
            bq = q0; bx=coords[1,v]; by=coords[2,v]; bz=coords[3,v]
            h = 0.25 * scale
            for _pass in 1:24
                improved = false
                for (dx,dy,dz) in dirs
                    coords[1,v]=bx+h*dx; coords[2,v]=by+h*dy; coords[3,v]=bz+h*dz
                    if _star_positive(m, coords, tv)
                        q = _star_quality(m, coords, tv)
                        if q > bq; bq=q; bx=coords[1,v]; by=coords[2,v]; bz=coords[3,v]; improved=true; end
                    end
                end
                coords[1,v]=bx; coords[2,v]=by; coords[3,v]=bz
                improved || (h *= 0.5)
                h < 1e-4 * scale && break
            end
        end
    end
    return Mesh(coords; segs=copy(m.segs), tris=copy(m.tris), tets=copy(m.tets),
                seg_tag=copy(m.seg_tag), tri_tag=copy(m.tri_tag), tet_tag=copy(m.tet_tag))
end

"""
    remove_slivers(m; max_rounds=8, sliver_deg=10.0) -> (Mesh, report)

Converging sliver-exudation driver (geometric route). Repeatedly applies the
targeted optimization-based smoother [`smooth_optimize`](@ref) — re-identifying the
poor vertex stars each round after the previous round's moves reshape their
neighbourhoods — and **accepts a round only if it is valid and strictly reduces the
sliver count**, so the result is monotone-non-worsening and stops at convergence.
Validity, total volume, boundary, and region tags are all preserved. Returns the
improved mesh and a report NamedTuple `(slivers_before, slivers_after,
min_dihedral_before, min_dihedral_after)`. The topological complement (2-3/3-2
flips that repair slivers connectivity-smoothing cannot reach) is
[`Mesh3D.optimize_flips!`](@ref); `mesh_volume(optimize=true)` runs both.
"""
function remove_slivers(m::Mesh; max_rounds::Integer=8, sliver_deg::Real=10.0)
    max_rounds >= 0 ||
        throw(ArgumentError("remove_slivers: max_rounds must be non-negative (got $max_rounds)"))
    q0 = mesh_quality(m; sliver_deg=sliver_deg)
    best = m; bestn = q0.n_slivers
    for _ in 1:max_rounds
        bestn == 0 && break
        cand = smooth_optimize(best; iters=4, sliver_deg=sliver_deg)
        (validate(cand).ok) || break
        q = mesh_quality(cand; sliver_deg=sliver_deg)
        (q.n_slivers < bestn) || break                # converged / no further gain
        best = cand; bestn = q.n_slivers
    end
    qf = mesh_quality(best; sliver_deg=sliver_deg)
    return best, (slivers_before = q0.n_slivers, slivers_after = qf.n_slivers,
                  min_dihedral_before = q0.min_dihedral_deg, min_dihedral_after = qf.min_dihedral_deg)
end

# star quality = min over incident tets of min(min_dihedral, π − max_dihedral): a
# sliver (needle OR cap) drives this toward 0. Maximizing it removes both sliver types.
@inline function _star_quality(m::Mesh, coords, tv)
    q = Inf
    @inbounds for t in tv
        a=m.tets[1,t];b=m.tets[2,t];c=m.tets[3,t];d=m.tets[4,t]
        pa=(coords[1,a],coords[2,a],coords[3,a]); pb=(coords[1,b],coords[2,b],coords[3,b])
        pc=(coords[1,c],coords[2,c],coords[3,c]); pd=(coords[1,d],coords[2,d],coords[3,d])
        dmn, dmx = tet_dihedral_extrema(pa,pb,pc,pd)
        q = min(q, dmn, pi - dmx)
    end
    return q
end

# all incident tets strictly positively oriented at the current coords?
@inline function _star_positive(m::Mesh, coords, tv)
    @inbounds for t in tv
        a=m.tets[1,t];b=m.tets[2,t];c=m.tets[3,t];d=m.tets[4,t]
        pa=(coords[1,a],coords[2,a],coords[3,a]); pb=(coords[1,b],coords[2,b],coords[3,b])
        pc=(coords[1,c],coords[2,c],coords[3,c]); pd=(coords[1,d],coords[2,d],coords[3,d])
        tet_signed_volume(pa,pb,pc,pd) <= 0 && return false
    end
    return true
end

# characteristic length of a vertex star = mean edge length over its incident tets.
@inline function _star_scale(m::Mesh, coords, tv)
    s = 0.0; n = 0
    @inbounds for t in tv
        vs = (m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
        for i in 1:4, j in i+1:4
            a=vs[i]; b=vs[j]
            dx=coords[1,a]-coords[1,b]; dy=coords[2,a]-coords[2,b]; dz=coords[3,a]-coords[3,b]
            s += sqrt(dx*dx+dy*dy+dz*dz); n += 1
        end
    end
    return n == 0 ? 0.0 : s / n
end

# Circumcenter of tet (a,b,c,d): solve |x−a|=|x−b|=|x−c|=|x−d|, denom = 2·A·(B×C) =
# 12·signed volume (same system as MeshTypes.tet_circumradius) but returning the
# absolute point. A degenerate (flat) tet gives denom≈0 → non-finite coords, which
# the caller filters out.
@inline function _tet_circumcenter(a, b, c, d)
    Ax=b[1]-a[1]; Ay=b[2]-a[2]; Az=b[3]-a[3]
    Bx=c[1]-a[1]; By=c[2]-a[2]; Bz=c[3]-a[3]
    Cx=d[1]-a[1]; Cy=d[2]-a[2]; Cz=d[3]-a[3]
    bcx=By*Cz-Bz*Cy; bcy=Bz*Cx-Bx*Cz; bcz=Bx*Cy-By*Cx        # B×C
    cax=Cy*Az-Cz*Ay; cay=Cz*Ax-Cx*Az; caz=Cx*Ay-Cy*Ax        # C×A
    abx=Ay*Bz-Az*By; aby=Az*Bx-Ax*Bz; abz=Ax*By-Ay*Bx        # A×B
    denom=2.0*(Ax*bcx+Ay*bcy+Az*bcz)                         # 12·signed volume
    la=Ax*Ax+Ay*Ay+Az*Az; lb=Bx*Bx+By*By+Bz*Bz; lc=Cx*Cx+Cy*Cy+Cz*Cz
    ox=(la*bcx+lb*cax+lc*abx)/denom
    oy=(la*bcy+lb*cay+lc*aby)/denom
    oz=(la*bcz+lb*caz+lc*abz)/denom
    return (a[1]+ox, a[2]+oy, a[3]+oz)                       # absolute circumcenter
end

# Scale-invariant "degenerate tet" volume threshold: a tet whose volume is a
# vanishing fraction of its longest-edge cube is numerically flat (circumcenter
# blows up). Per-tet, so non-uniform meshes are handled correctly.
@inline function _odt_voleps(a, b, c, d)
    P = (a, b, c, d); e2 = 0.0
    @inbounds for i in 1:4, j in i+1:4
        dx=P[i][1]-P[j][1]; dy=P[i][2]-P[j][2]; dz=P[i][3]-P[j][3]
        l2=dx*dx+dy*dy+dz*dz; l2 > e2 && (e2 = l2)
    end
    return 1e-12 * e2 * sqrt(e2)                             # 1e-12 · (max edge)³
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
