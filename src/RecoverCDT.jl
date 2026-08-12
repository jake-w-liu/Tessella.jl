"""
    RecoverCDT

General robust boundary recovery by **conforming-Delaunay refinement on the exact
kernel** — the piece that closes the one class `Mesh3D.recover_boundary` could not
mesh: a **non-star + reflex** polyhedron (the twisted prism).

`recover_boundary_cdt(surface)` recovers an arbitrary closed PLC surface as a
**conforming** interior tetrahedral mesh (every input facet a union of tet faces, by
exact-arithmetic conformity) or throws an explicit blocker — never a silent
non-conforming mesh (PLAN principle #4).

Algorithm (conforming Delaunay by Gabriel encroachment, Shewchuk — driven directly off
the exact 3-D Delaunay rather than a conservative diametral proxy, which is what makes it
*converge* on thin small-angle features instead of diverging):

 1. exact-Delaunay the current rational vertex set (`ExactMesh3D.delaunay3d_exact`);
 2. any crease **subsegment** that is not a Delaunay edge is encroached → split at its
    exact rational **midpoint**;
 3. once all creases are present, any Delaunay **edge that pierces a coplanar region's
    interior** means that facet-region is not yet recovered → insert the **exact rational
    point where that edge crosses the region plane** (a Steiner point exactly on the facet);
 4. repeat until no crease is missing and no edge pierces a region;
 5. assemble (interior tets kept by centroid ray-cast against the input surface) and pass
    an **exact `Rational{BigInt}` conformity certificate** (crease edges present, no edge
    pierces a region, every boundary face lies in the surface, per-region exact projected
    area matches) — accept only if it passes.

All Steiner coordinates are exact `Rational{BigInt}` (midpoints / edge–plane
intersections of rational points), so they stay exactly on-feature — impossible with a
Float64 kernel. Verified conforming (exact volume + area) on box, genus-1 tunnel, hollow
shell, faceted cylinder, and the **twisted (non-star+reflex) prism**.
"""
module RecoverCDT

using ..ExactMesh3D: delaunay3d_exact
using ..MeshTypes: Mesh, validate, is_closed_manifold, boundary_faces, triangle_area,
                   node, ntris, ntets, nnodes, tet_volume
import ..Mesh3D

export recover_boundary_cdt, mesh_sized_cdt

const RB = Rational{BigInt}
const RBv = NTuple{3,RB}

@inline rsub(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline radd(a,b)=(a[1]+b[1],a[2]+b[2],a[3]+b[3])
@inline rdot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline rcross(a,b)=(a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
@inline rscale(a,s)=(a[1]*s,a[2]*s,a[3]*s)
@inline rsign(x)=x>0 ? 1 : (x<0 ? -1 : 0)
@inline det3(a,b,c)=a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
@inline ekey(a,b)=a<=b ? (a,b) : (b,a)
@inline proj2(p,drop)= drop==1 ? (p[2],p[3]) : (drop==2 ? (p[1],p[3]) : (p[1],p[2]))
@inline function orient2(a,b,c)
    rsign((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))
end
# exact 2D projected (unsigned) area of triangle, dropping axis `drop`
@inline function parea2(a::RBv,b::RBv,c::RBv,drop::Int)
    a2=proj2(a,drop); b2=proj2(b,drop); c2=proj2(c,drop)
    d=(b2[1]-a2[1])*(c2[2]-a2[2])-(b2[2]-a2[2])*(c2[1]-a2[1])
    (d<0 ? -d : d)//2
end

struct Plane; A0::RBv; N::RBv; drop::Int; end
mutable struct Region
    plane::Plane
    tris::Vector{NTuple{3,Int}}       # ORIGINAL facet triangles of the region
    bnd::Set{NTuple{2,Int}}           # current boundary subsegments (atomic)
    intedges::Set{NTuple{2,Int}}      # interior edges of the region's input triangulation
end
@inline side(pl::Plane, x)=rsign(rdot(pl.N, rsub(x, pl.A0)))

function pierces_region(reg::Region, pts, y)
    drop=reg.plane.drop; y2=proj2(y,drop)
    for (a,b,c) in reg.tris
        pa=proj2(pts[a],drop); pb=proj2(pts[b],drop); pc=proj2(pts[c],drop)
        s1=orient2(pa,pb,y2); s2=orient2(pb,pc,y2); s3=orient2(pc,pa,y2)
        (s1!=0 && s2!=0 && s3!=0 && s1==s2 && s2==s3) && return true
        zc=(s1==0)+(s2==0)+(s3==0)
        if zc==1
            nz=filter(!=(0),(s1,s2,s3))
            if length(nz)==2 && nz[1]==nz[2]
                e = s1==0 ? ekey(a,b) : (s2==0 ? ekey(b,c) : ekey(c,a))
                e in reg.intedges && return true
            end
        end
    end
    false
end
# closed (inside or on boundary) point-in-region test (for boundary-in-surface certification)
function in_region_closed(reg::Region, pts, y)
    drop=reg.plane.drop; y2=proj2(y,drop)
    for (a,b,c) in reg.tris
        pa=proj2(pts[a],drop); pb=proj2(pts[b],drop); pc=proj2(pts[c],drop)
        s1=orient2(pa,pb,y2); s2=orient2(pb,pc,y2); s3=orient2(pc,pa,y2)
        hp=(s1>0)||(s2>0)||(s3>0); hn=(s1<0)||(s2<0)||(s3<0)
        (hp && hn) || return true
    end
    false
end

function edge_pierces(reg::Region, pts, p, q)
    pp=pts[p]; qq=pts[q]; sp=side(reg.plane,pp); sq=side(reg.plane,qq)
    (sp!=0 && sq!=0 && sp!=sq) || return (false, pp)
    num=rdot(reg.plane.N, rsub(reg.plane.A0, pp)); den=rdot(reg.plane.N, rsub(qq,pp))
    den==0 && return (false, pp)
    t=num//den
    y=(pp[1]+t*(qq[1]-pp[1]), pp[2]+t*(qq[2]-pp[2]), pp[3]+t*(qq[3]-pp[3]))
    pierces_region(reg,pts,y) ? (true,y) : (false,y)
end

function dedup(surface)
    seen=Dict{NTuple{3,Float64},Int}(); pts=RBv[]; remap=zeros(Int,size(surface.coords,2))
    for i in 1:size(surface.coords,2)
        k=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i]); id=get(seen,k,0)
        if id==0; push!(pts,(RB(k[1]),RB(k[2]),RB(k[3]))); id=length(pts); seen[k]=id; end
        remap[i]=id
    end
    facets=NTuple{3,Int}[]
    for f in 1:size(surface.tris,2)
        a=remap[surface.tris[1,f]];b=remap[surface.tris[2,f]];c=remap[surface.tris[3,f]]
        (a==b||b==c||a==c)&&continue; push!(facets,(a,b,c))
    end
    pts,facets
end
mutable struct UF; p::Vector{Int}; end
uf(n)=UF(collect(1:n)); ufind(u,i)=(while u.p[i]!=i; u.p[i]=u.p[u.p[i]]; i=u.p[i]; end; i); ufunion(u,i,j)=(u.p[ufind(u,i)]=ufind(u,j))

function build_regions(pts, facets)
    nf=length(facets); e2f=Dict{NTuple{2,Int},Vector{Int}}()
    for (fi,(a,b,c)) in enumerate(facets); for e in (ekey(a,b),ekey(b,c),ekey(a,c)); push!(get!(e2f,e,Int[]),fi); end; end
    u=uf(nf)
    for (e,fs) in e2f
        length(fs)==2 || continue
        f1=facets[fs[1]]; f2=facets[fs[2]]; (p,q)=e
        w1=f1[1]; (w1==p||w1==q)&&(w1=f1[2]); (w1==p||w1==q)&&(w1=f1[3])
        w2=f2[1]; (w2==p||w2==q)&&(w2=f2[2]); (w2==p||w2==q)&&(w2=f2[3])
        det3(rsub(pts[q],pts[p]),rsub(pts[w1],pts[p]),rsub(pts[w2],pts[p]))==0 && ufunion(u,fs[1],fs[2])
    end
    comp=Dict{Int,Vector{Int}}(); for fi in 1:nf; push!(get!(comp,ufind(u,fi),Int[]),fi); end
    regions=Region[]
    for (_,fis) in comp
        tris=[facets[fi] for fi in fis]; a,b,c=tris[1]
        A0=pts[a]; N=rcross(rsub(pts[b],pts[a]),rsub(pts[c],pts[a]))
        an=(abs(N[1]),abs(N[2]),abs(N[3])); drop = an[1]>=an[2] ? (an[1]>=an[3] ? 1 : 3) : (an[2]>=an[3] ? 2 : 3)
        inc=Dict{NTuple{2,Int},Int}()
        for (x,y,z) in tris, e in (ekey(x,y),ekey(y,z),ekey(x,z)); inc[e]=get(inc,e,0)+1; end
        bnd=Set{NTuple{2,Int}}(); intl=Set{NTuple{2,Int}}()
        for (e,ci) in inc; ci==1 ? push!(bnd,e) : push!(intl,e); end
        push!(regions, Region(Plane(A0,N,drop), tris, bnd, intl))
    end
    regions
end

# ---- exact certification: does mesh m conform to the input surface, in exact arithmetic? ----
function certify_exact(surface::Mesh, m::Mesh, pts::Vector{RBv}, regions, seg2regs)
    v=validate(m); v.ok || return (false, "invalid mesh: "*join(v.messages,"; "))
    is_closed_manifold(m) || return (false, "not closed manifold")
    key2pid=Dict{NTuple{3,Float64},Int}()
    for i in 1:length(pts); key2pid[(Float64(pts[i][1]),Float64(pts[i][2]),Float64(pts[i][3]))]=i; end
    mid=Vector{Int}(undef, size(m.coords,2))
    for i in 1:size(m.coords,2)
        k=(m.coords[1,i],m.coords[2,i],m.coords[3,i])
        haskey(key2pid,k) || return (false,"mesh node $i has no exact-rational preimage")
        mid[i]=key2pid[k]
    end
    E=Set{NTuple{2,Int}}()
    for t in 1:ntets(m)
        vs=(mid[m.tets[1,t]],mid[m.tets[2,t]],mid[m.tets[3,t]],mid[m.tets[4,t]])
        for i in 1:4,j in i+1:4; push!(E,ekey(vs[i],vs[j])); end
    end
    for e in keys(seg2regs)
        (e in E) || return (false,"crease subsegment $e not a tet edge")
    end
    for (ri,reg) in enumerate(regions), e in E
        pr,_ = edge_pierces(reg, pts, e[1], e[2])
        pr && return (false,"tet edge $e pierces region $ri interior")
    end
    bf = first(boundary_faces(m.tets))
    for f in bf
        a=mid[f[1]]; b=mid[f[2]]; c=mid[f[3]]
        cen=rscale(radd(radd(pts[a],pts[b]),pts[c]), 1//3)
        found=false
        for reg in regions
            (side(reg.plane,pts[a])==0 && side(reg.plane,pts[b])==0 && side(reg.plane,pts[c])==0) || continue
            in_region_closed(reg,pts,cen) && (found=true; break)
        end
        found || return (false,"boundary face ($a,$b,$c) lies outside the input surface")
    end
    # EXACT per-region area conformity: rational projected boundary-face area == input facet area
    for (ri,reg) in enumerate(regions)
        want=zero(RB)
        for (a,b,c) in reg.tris; want += parea2(pts[a],pts[b],pts[c],reg.plane.drop); end
        got=zero(RB)
        for f in bf
            a=mid[f[1]]; b=mid[f[2]]; c=mid[f[3]]
            (side(reg.plane,pts[a])==0 && side(reg.plane,pts[b])==0 && side(reg.plane,pts[c])==0) || continue
            got += parea2(pts[a],pts[b],pts[c],reg.plane.drop)
        end
        want==got || return (false,"region $ri exact area mismatch: want $(Float64(want)) got $(Float64(got))")
    end
    return (true,"conforming (exact)")
end

function _finalize(surface, pts, etets)
    g=Mesh3D._raygrid(surface)
    coordsF=Matrix{Float64}(undef,3,length(pts))
    for i in 1:length(pts); coordsF[1,i]=Float64(pts[i][1]);coordsF[2,i]=Float64(pts[i][2]);coordsF[3,i]=Float64(pts[i][3]); end
    keep=NTuple{4,Int}[]
    for t in etets
        cx=(coordsF[1,t[1]]+coordsF[1,t[2]]+coordsF[1,t[3]]+coordsF[1,t[4]])/4
        cy=(coordsF[2,t[1]]+coordsF[2,t[2]]+coordsF[2,t[3]]+coordsF[2,t[4]])/4
        cz=(coordsF[3,t[1]]+coordsF[3,t[2]]+coordsF[3,t[3]]+coordsF[3,t[4]])/4
        Mesh3D._inside_grid((cx,cy,cz),g) && push!(keep,t)
    end
    tetM=Matrix{Int32}(undef,4,length(keep))
    for (j,t) in enumerate(keep); tetM[1,j]=t[1];tetM[2,j]=t[2];tetM[3,j]=t[3];tetM[4,j]=t[4]; end
    Mesh(coordsF; tets=tetM)
end

"""
    recover_boundary_cdt(surface::Mesh; maxiter=2000, maxpts=6000) -> Mesh

Conforming tetrahedralization of a closed PLC `surface` by exact-kernel conforming-
Delaunay refinement (see the module docstring). Returns a validated, closed-manifold,
**exactly-conforming** tet mesh, or throws an explicit blocker (never a silent
non-conforming mesh). Closes the non-star+reflex class the Float64 `recover_boundary`
cannot; also recovers the supported classes (box / genus-1 tunnel / hollow shell /
faceted cylinder) with 0 Steiner points where they are already Delaunay-conforming.
"""
function recover_boundary_cdt(surface::Mesh; maxiter::Integer=2000, maxpts::Integer=6000)
    m, _, _, _, _ = _recover(surface; maxiter=maxiter, maxpts=maxpts)
    return m
end

# internal: run the refinement and return (mesh, pts, regions, seg2regs, n0) so sizing can reuse it.
function _recover(surface::Mesh; maxiter::Integer=2000, maxpts::Integer=6000)
    pts, facets = dedup(surface)
    length(pts)>=4 || throw(ArgumentError("recover_boundary_cdt: need ≥ 4 distinct vertices"))
    isempty(facets) && throw(ArgumentError("recover_boundary_cdt: surface has no facets"))
    n0=length(pts)
    regions=build_regions(pts,facets)
    seg2regs=Dict{NTuple{2,Int},Vector{Int}}()
    for (ri,reg) in enumerate(regions), e in reg.bnd; push!(get!(seg2regs,e,Int[]),ri); end

    function split_seg!(e)
        a,b=e; mid=rscale(radd(pts[a],pts[b]),1//2); push!(pts,mid); m=length(pts)
        regs=seg2regs[e]; delete!(seg2regs,e)
        for e2 in (ekey(a,m),ekey(m,b)); seg2regs[e2]=copy(regs); end
        for ri in regs; reg=regions[ri]; delete!(reg.bnd,e); push!(reg.bnd,ekey(a,m)); push!(reg.bnd,ekey(m,b)); end
        m
    end

    local etets
    for it in 1:maxiter
        length(pts) > maxpts && error("recover_boundary_cdt: maxpts exceeded ($(length(pts))) — refinement not converging")
        etets = delaunay3d_exact(pts)
        E=Set{NTuple{2,Int}}()
        for t in etets; for i in 1:4,j in i+1:4; push!(E,ekey(t[i],t[j])); end; end
        # 1) missing crease subsegments -> midpoint split
        missing=NTuple{2,Int}[]
        for e in keys(seg2regs); (e in E) || push!(missing,e); end
        if !isempty(missing)
            for e in missing; haskey(seg2regs,e) && split_seg!(e); end
            continue
        end
        # 2) facet piercing -> insert exact piercing point (one per region per iter)
        newpts=RBv[]
        for reg in regions
            for e in E
                pr,y=edge_pierces(reg,pts,e[1],e[2])
                pr && (push!(newpts,y); break)
            end
        end
        if !isempty(newpts)
            for y in unique(newpts); push!(pts,y); end
            continue
        end
        # converged: assemble + EXACT certification (never return a silently non-conforming mesh)
        m = _finalize(surface, pts, etets)
        ok, reason = certify_exact(surface, m, pts, regions, seg2regs)
        ok && return (m, pts, regions, seg2regs, n0)
        error("recover_boundary_cdt: convergence gate passed but exact certification failed: $reason")
    end
    error("recover_boundary_cdt: no convergence in $maxiter iters (verts=$(length(pts)))")
end

"""
    mesh_sized_cdt(surface::Mesh; hmax) -> Mesh

Interior size control on an **arbitrary** (curved / non-star) domain, on top of the
exact conforming-Delaunay recovery. First recovers the boundary conformingly
([`recover_boundary_cdt`](@ref)), then seeds an interior lattice of exact-rational
Steiner points (spacing `hmax/√3`, strictly inside via ray-cast, kept clear of existing
vertices), re-tetrahedralizes with the exact kernel, and **accepts only if the exact
conformity certificate still passes** — otherwise it returns the conforming boundary
mesh unchanged (never a silently non-conforming mesh). Interior edges are driven toward
`≤ hmax`; the surface-facet edges are fixed by the input and bound the achievable size
(refine the surface for finer than that). Always validated + exactly-conforming.
"""
function mesh_sized_cdt(surface::Mesh; hmax::Real)
    hmax > 0 || throw(ArgumentError("mesh_sized_cdt: hmax must be positive (got $hmax)"))
    m0, pts, regions, seg2regs, n0 = _recover(surface)
    g = Mesh3D._raygrid(surface)
    xs=[Float64(p[1]) for p in pts]; ys=[Float64(p[2]) for p in pts]; zs=[Float64(p[3]) for p in pts]
    xlo,xhi=extrema(xs); ylo,yhi=extrema(ys); zlo,zhi=extrema(zs)
    a=float(hmax)/sqrt(3.0)
    nx=max(1,ceil(Int,(xhi-xlo)/a)); ny=max(1,ceil(Int,(yhi-ylo)/a)); nz=max(1,ceil(Int,(zhi-zlo)/a))
    ins2=(0.5*float(hmax))^2
    aug=copy(pts)
    @inbounds for i in 0:nx, j in 0:ny, k in 0:nz
        px=xlo+i*(xhi-xlo)/nx; py=ylo+j*(yhi-ylo)/ny; pz=zlo+k*(zhi-zlo)/nz
        Mesh3D._inside_grid((px,py,pz),g) || continue
        keep=true
        for v in 1:length(aug)
            (px-Float64(aug[v][1]))^2+(py-Float64(aug[v][2]))^2+(pz-Float64(aug[v][3]))^2 < ins2 && (keep=false; break)
        end
        keep && push!(aug,(RB(px),RB(py),RB(pz)))
    end
    length(aug) == length(pts) && return m0                       # nothing to add
    etets = delaunay3d_exact(aug)
    m = _finalize(surface, aug, etets)
    ok, _ = certify_exact(surface, m, aug, regions, seg2regs)
    return ok ? m : m0                                            # sized-or-conforming-baseline (never non-conforming)
end

end # module RecoverCDT
