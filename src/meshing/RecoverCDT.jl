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
import ..Mesh2D

export recover_boundary_cdt, recover_partition_cdt, mesh_sized_cdt

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
        k=(surface.coords[1,i]==0 ? 0.0 : surface.coords[1,i],
           surface.coords[2,i]==0 ? 0.0 : surface.coords[2,i],
           surface.coords[3,i]==0 ? 0.0 : surface.coords[3,i]); id=get(seen,k,0)
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

# Deduplicate a collection of closed surfaces into one exact point/facet registry.
# Coincident triangles (the two copies of a shared material interface) are constraints
# only once.
function dedup_partition(surfaces::AbstractVector{Mesh})
    seen=Dict{NTuple{3,Float64},Int}(); pts=RBv[]
    allfacets=NTuple{3,Int}[]; facekeys=Set{NTuple{3,Int}}()
    for surface in surfaces
        remap=Vector{Int}(undef,size(surface.coords,2))
        for i in axes(surface.coords,2)
            # Canonicalise signed zero: they are the same geometric point and must not
            # become distinct exact Rational vertices.
            k=(surface.coords[1,i]==0 ? 0.0 : surface.coords[1,i],
               surface.coords[2,i]==0 ? 0.0 : surface.coords[2,i],
               surface.coords[3,i]==0 ? 0.0 : surface.coords[3,i])
            id=get(seen,k,0)
            if id==0
                push!(pts,(RB(k[1]),RB(k[2]),RB(k[3])))
                id=length(pts);seen[k]=id
            end
            remap[i]=id
        end
        for f in axes(surface.tris,2)
            tri=(remap[surface.tris[1,f]],remap[surface.tris[2,f]],remap[surface.tris[3,f]])
            q=sort(Int[tri[1],tri[2],tri[3]]);key=(q[1],q[2],q[3])
            if !(key in facekeys)
                push!(facekeys,key);push!(allfacets,tri)
            end
        end
    end
    pts,allfacets
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
    for i in 1:length(pts)
        k=(Float64(pts[i][1]),Float64(pts[i][2]),Float64(pts[i][3]))
        (isfinite(k[1])&&isfinite(k[2])&&isfinite(k[3])) ||
            return (false,"exact point $i is not finite in Float64 output coordinates")
        haskey(key2pid,k) && key2pid[k]!=i &&
            return (false,"distinct exact points $(key2pid[k]) and $i collide in Float64 output coordinates")
        key2pid[k]=i
    end
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
            # count only faces belonging to THIS region — not merely coplanar with its plane.
            # Two DISCONNECTED coplanar regions (a slot / U / plus / multi-component surface) share
            # a geometric plane, so a plane-only filter double-counts the sibling region's faces and
            # spuriously fails the area certificate; the centroid membership test scopes it correctly.
            cen=rscale(radd(radd(pts[a],pts[b]),pts[c]), 1//3)
            in_region_closed(reg,pts,cen) || continue
            got += parea2(pts[a],pts[b],pts[c],reg.plane.drop)
        end
        want==got || return (false,"region $ri exact area mismatch: want $(Float64(want)) got $(Float64(got))")
    end
    return (true,"conforming (exact)")
end

function _finalize(surface, pts, etets)
    length(pts)<=typemax(Int32) ||
        throw(ArgumentError("recover_boundary_cdt: output node count exceeds Int32"))
    g=Mesh3D._raygrid(surface)
    coordsF=Matrix{Float64}(undef,3,length(pts))
    seen=Dict{NTuple{3,Float64},Int}()
    for i in 1:length(pts)
        p=(Float64(pts[i][1]),Float64(pts[i][2]),Float64(pts[i][3]))
        (isfinite(p[1])&&isfinite(p[2])&&isfinite(p[3])) ||
            throw(ErrorException("recover_boundary_cdt: exact point $i cannot be represented by finite Float64 output"))
        haskey(seen,p) && throw(ErrorException(
            "recover_boundary_cdt: distinct exact points $(seen[p]) and $i collapse to one Float64 coordinate"))
        seen[p]=i;coordsF[1,i]=p[1];coordsF[2,i]=p[2];coordsF[3,i]=p[3]
    end
    keep=NTuple{4,Int}[]
    for t in etets
        cx=Float64((pts[t[1]][1]+pts[t[2]][1]+pts[t[3]][1]+pts[t[4]][1])/4)
        cy=Float64((pts[t[1]][2]+pts[t[2]][2]+pts[t[3]][2]+pts[t[4]][2])/4)
        cz=Float64((pts[t[1]][3]+pts[t[2]][3]+pts[t[3]][3]+pts[t[4]][3])/4)
        Mesh3D._inside_grid((cx,cy,cz),g) && push!(keep,t)
    end
    tetM=Matrix{Int32}(undef,4,length(keep))
    for (j,t) in enumerate(keep); tetM[1,j]=t[1];tetM[2,j]=t[2];tetM[3,j]=t[3];tetM[4,j]=t[4]; end
    Mesh(coordsF; tets=tetM)
end

function _float_coords(pts,caller::AbstractString)
    length(pts)<=typemax(Int32) ||
        throw(ArgumentError("$caller: output node count exceeds Int32"))
    coords=Matrix{Float64}(undef,3,length(pts));seen=Dict{NTuple{3,Float64},Int}()
    for i in eachindex(pts)
        p=(Float64(pts[i][1]),Float64(pts[i][2]),Float64(pts[i][3]))
        (isfinite(p[1])&&isfinite(p[2])&&isfinite(p[3])) ||
            throw(ErrorException("$caller: exact point $i cannot be represented by finite Float64 output"))
        haskey(seen,p) && throw(ErrorException(
            "$caller: distinct exact points $(seen[p]) and $i collapse to one Float64 coordinate"))
        seen[p]=i;coords[1,i]=p[1];coords[2,i]=p[2];coords[3,i]=p[3]
    end
    coords
end

function _partition_tets_float(pts;seed::Int=1)
    caller="recover_partition_cdt"
    coords=_float_coords(pts,caller)
    lastreason="no triangulation attempted"
    # First use the exact kernel on the bounded-size dyadic Float64 images.  Unlike
    # repeatedly feeding it arbitrary intersection rationals, operand sizes cannot
    # grow from one recovery iteration to the next; its completeness gate also avoids
    # the physically flat cells possible in the ghost-based Float kernel.
    try
        dyadic=RBv[(RB(coords[1,i]),RB(coords[2,i]),RB(coords[3,i])) for i in axes(coords,2)]
        candidate=delaunay3d_exact(dyadic);etets=NTuple{4,Int}[];used=falses(length(pts));flat=false
        for q0 in candidate
            q=q0;s=det3(rsub(pts[q[2]],pts[q[1]]),rsub(pts[q[3]],pts[q[1]]),rsub(pts[q[4]],pts[q[1]]))
            if s==0;flat=true;break;end
            s<0&&(q=(q[1],q[2],q[4],q[3]))
            push!(etets,q);for v in q;used[v]=true;end
        end
        if !flat&&!isempty(etets)&&all(used)
            tm=Matrix{Int32}(undef,4,length(etets))
            for (j,t) in enumerate(etets),i in 1:4;tm[i,j]=Int32(t[i]);end
            diag=validate(Mesh(coords;tets=tm))
            diag.ok&&return etets
            lastreason="dyadic-exact topology was invalid — "*join(diag.messages,"; ")
        elseif flat
            lastreason="dyadic-exact topology contained a tetrahedron flat in the exact preimage"
        else
            lastreason="dyadic-exact topology was empty or incomplete"
        end
    catch err
        err isa InterruptException&&rethrow()
        (err isa ArgumentError||err isa ErrorException)||rethrow()
        lastreason="dyadic-exact topology failed — "*sprint(showerror,err)
    end
    # The unperturbed topology is preferred.  Cospherical/coplanar Float inputs can
    # make its SoS connectivity physically flat, so bounded deterministic insertion-
    # order trials use the kernel's coordinate perturbation only to choose a topology;
    # all vertices and every acceptance predicate remain the unperturbed exact points.
    for attempt in 0:16
        perturb=attempt!=0;trialseed=attempt==0 ? seed : seed+attempt-1
        try
            T=Mesh3D.delaunay3d(collect(coords[1,:]),collect(coords[2,:]),collect(coords[3,:]);
                                rng_seed=trialseed,perturb=perturb)
            etets=NTuple{4,Int}[];used=falses(length(pts));flat=false
            for ti in eachindex(T.alive)
                (T.alive[ti]&&!Mesh3D._is_ghost_tet(T,ti))||continue
                q=(Int(Mesh3D._vert(T,ti,1)),Int(Mesh3D._vert(T,ti,2)),
                   Int(Mesh3D._vert(T,ti,3)),Int(Mesh3D._vert(T,ti,4)))
                s=det3(rsub(pts[q[2]],pts[q[1]]),rsub(pts[q[3]],pts[q[1]]),rsub(pts[q[4]],pts[q[1]]))
                if s==0;flat=true;break;end
                s<0&&(q=(q[1],q[2],q[4],q[3]))
                push!(etets,q);for v in q;used[v]=true;end
            end
            flat&&(lastreason="attempt $attempt contained an exactly flat tetrahedron";continue)
            isempty(etets)&&(lastreason="attempt $attempt produced no tetrahedra";continue)
            all(used)||(lastreason="attempt $attempt omitted $(count(!,used)) exact points";continue)
            tm=Matrix{Int32}(undef,4,length(etets))
            for (j,t) in enumerate(etets),i in 1:4;tm[i,j]=Int32(t[i]);end
            diag=validate(Mesh(coords;tets=tm))
            if diag.ok;return etets;end
            lastreason="attempt $attempt was invalid — "*join(diag.messages,"; ")
        catch err
            err isa InterruptException&&rethrow()
            (err isa ArgumentError||err isa ErrorException)||rethrow()
            lastreason="attempt $attempt failed — "*sprint(showerror,err)
        end
    end
    throw(ErrorException("$caller: no valid Float-assisted topology in 17 deterministic attempts; $lastreason"))
end

function _certify_partition_exact(surfaces,m,preimage,regions,regionpts=preimage)
    caller="recover_partition_cdt"
    for r in eachindex(surfaces)
        ids=Int[t for t in axes(m.tets,2) if m.tet_tag[t]==Int32(r)]
        isempty(ids)&&return (false,"effective region $r received no tet")
        tm=m.tets[:,ids];part=Mesh(m.coords;tets=tm)
        d=validate(part);d.ok||return (false,"effective region $r is not a tetrahedral manifold — "*join(d.messages,"; "))
        bf=first(boundary_faces(tm))
        for f in bf
            a=Int(f[1]);b=Int(f[2]);c=Int(f[3])
            pa=preimage[a];pb=preimage[b];pc=preimage[c]
            cen=rscale(radd(radd(pa,pb),pc),1//3)
            found=false
            for reg in regions
                (side(reg.plane,pa)==0&&side(reg.plane,pb)==0&&side(reg.plane,pc)==0)||continue
                if in_region_closed(reg,regionpts,cen);found=true;break;end
            end
            found||return (false,"effective region $r has boundary face ($a,$b,$c) outside every input PLC")
        end
    end
    (true,"conforming partition (exact preimage certificate)")
end

function _finalize_partition(surfaces,pts,etets,regions)
    caller="recover_partition_cdt"
    coords=_float_coords(pts,caller)
    grids=[Mesh3D._raygrid(s) for s in surfaces]
    bboxes=[begin
        lo=(Inf,Inf,Inf);hi=(-Inf,-Inf,-Inf)
        for i in axes(s.coords,2)
            lo=(min(lo[1],s.coords[1,i]),min(lo[2],s.coords[2,i]),min(lo[3],s.coords[3,i]))
            hi=(max(hi[1],s.coords[1,i]),max(hi[2],s.coords[2,i]),max(hi[3],s.coords[3,i]))
        end
        (lo,hi)
    end for s in surfaces]
    keep=NTuple{4,Int}[];tags=Int32[]
    sizehint!(keep,length(etets));sizehint!(tags,length(etets))
    for t in etets
        c=((pts[t[1]][1]+pts[t[2]][1]+pts[t[3]][1]+pts[t[4]][1])/4,
           (pts[t[1]][2]+pts[t[2]][2]+pts[t[3]][2]+pts[t[4]][2])/4,
           (pts[t[1]][3]+pts[t[2]][3]+pts[t[3]][3]+pts[t[4]][3])/4)
        cf=(Float64(c[1]),Float64(c[2]),Float64(c[3]))
        (isfinite(cf[1])&&isfinite(cf[2])&&isfinite(cf[3])) ||
            throw(ErrorException("$caller: exact tetrahedron centroid is not finite in Float64"))
        for r in eachindex(surfaces)
            lo,hi=bboxes[r]
            (lo[1]<=cf[1]<=hi[1]&&lo[2]<=cf[2]<=hi[2]&&lo[3]<=cf[3]<=hi[3]) || continue
            if Mesh3D._inside_grid(cf,grids[r])
                push!(keep,t);push!(tags,Int32(r));break
            end
        end
    end
    isempty(keep) && throw(ErrorException("$caller: the classified partition is empty"))
    tets=Matrix{Int32}(undef,4,length(keep))
    for (j,t) in enumerate(keep),i in 1:4;tets[i,j]=Int32(t[i]);end
    out=Mesh(coords;tets=tets,tet_tag=tags)
    diag=validate(out)
    diag.ok || throw(ErrorException("$caller: produced an invalid mesh — "*join(diag.messages,"; ")))
    counts=Dict{Int32,Int}()
    for tag in tags;counts[tag]=get(counts,tag,0)+1;end
    for r in eachindex(surfaces)
        get(counts,Int32(r),0)>0 || throw(ErrorException("$caller: effective region $r received no tet"))
    end
    ok,reason=_certify_partition_exact(surfaces,out,pts,regions)
    ok||throw(ErrorException("$caller: exact partition certification failed — $reason"))
    out
end

@inline _rsub2(a,b)=(a[1]-b[1],a[2]-b[2])
@inline _rcross2(a,b)=a[1]*b[2]-a[2]*b[1]
@inline _project_axis(p,ax)=ax==1 ? (p[2],p[3]) : (ax==2 ? (p[1],p[3]) : (p[1],p[2]))

function _prismatic_axis(surfaces)
    for ax in 1:3
        good=true
        for s in surfaces,f in axes(s.tris,2)
            ids=(s.tris[1,f],s.tris[2,f],s.tris[3,f])
            q=(s.coords[ax,ids[1]],s.coords[ax,ids[2]],s.coords[ax,ids[3]])
            (q[1]==q[2]&&q[2]==q[3])&&continue
            p1=_project_axis((s.coords[1,ids[1]],s.coords[2,ids[1]],s.coords[3,ids[1]]),ax)
            p2=_project_axis((s.coords[1,ids[2]],s.coords[2,ids[2]],s.coords[3,ids[2]]),ax)
            p3=_project_axis((s.coords[1,ids[3]],s.coords[2,ids[3]],s.coords[3,ids[3]]),ax)
            _rcross2(_rsub2(p2,p1),_rsub2(p3,p1))==0||(good=false;break)
        end
        good&&return ax
    end
    0
end

@inline function _on_segment2(p,a,b)
    _rcross2(_rsub2(p,a),_rsub2(b,a))==0&&
        min(a[1],b[1])<=p[1]<=max(a[1],b[1])&&min(a[2],b[2])<=p[2]<=max(a[2],b[2])
end

# Exact planar arrangement + globally conforming prism extrusion for collections of
# co-axial prismatic PLCs (boxes, slots, polygonal cylinders, and shells).  This avoids
# the small-angle nontermination of conforming-Delaunay edge splitting at a feed-through
# while retaining exact preimages for every line/plane intersection.
function _recover_prismatic_partition(surfaces,nmax::Int)
    ax=_prismatic_axis(surfaces);ax==0&&return nothing
    dims=ax==1 ? (2,3) : (ax==2 ? (1,3) : (1,2))
    segset=Set{Tuple{NTuple{2,RB},NTuple{2,RB}}}()
    levels=RB[]
    for s in surfaces
        append!(levels,(RB(s.coords[ax,i]) for i in axes(s.coords,2)))
        for f in axes(s.tris,2)
            v=(s.tris[1,f],s.tris[2,f],s.tris[3,f])
            av=(s.coords[ax,v[1]],s.coords[ax,v[2]],s.coords[ax,v[3]])
            # Cap triangulation diagonals are not cross-section features; retaining
            # them creates artificial near-coincident intersections when a supplied
            # decimal midpoint is not the exact dyadic midpoint of its endpoints.
            (av[1]==av[2]&&av[2]==av[3])&&continue
            for (i,j) in ((1,2),(2,3),(3,1))
                a=(RB(s.coords[dims[1],v[i]]),RB(s.coords[dims[2],v[i]]))
                b=(RB(s.coords[dims[1],v[j]]),RB(s.coords[dims[2],v[j]]))
                a==b&&continue
                if isless(b,a);a,b=b,a;end
                push!(segset,(a,b))
            end
        end
    end
    segs=sort!(collect(segset));nseg=length(segs)
    nseg>0||throw(ErrorException("recover_partition_cdt: prismatic arrangement has no projected segments"))
    pairs=try Base.checked_mul(nseg,nseg-1)÷2 catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("recover_partition_cdt: projected arrangement pair count overflows Int"))
    end
    pairs<=5_000_000||throw(ArgumentError(
        "recover_partition_cdt: projected arrangement has $pairs segment pairs, exceeding the bounded exact-arrangement budget"))
    splits=[Set{NTuple{2,RB}}((a,b)) for (a,b) in segs]
    for i in 1:nseg-1,j in i+1:nseg
        a,b=segs[i];c,d=segs[j];r=_rsub2(b,a);s=_rsub2(d,c);den=_rcross2(r,s)
        if den!=0
            ca=_rsub2(c,a);t=_rcross2(ca,s)/den;u=_rcross2(ca,r)/den
            if 0<=t<=1&&0<=u<=1
                p=(a[1]+t*r[1],a[2]+t*r[2]);push!(splits[i],p);push!(splits[j],p)
            end
        elseif _rcross2(_rsub2(c,a),r)==0
            for p in (a,b,c,d)
                _on_segment2(p,a,b)&&_on_segment2(p,c,d)&&(push!(splits[i],p);push!(splits[j],p))
            end
        end
    end
    p2=NTuple{2,RB}[];pid=Dict{NTuple{2,RB},Int}()
    function getpid(p)
        get!(pid,p) do
            push!(p2,p);length(p2)
        end
    end
    atomic=Set{NTuple{2,Int}}()
    for (si,(a,b)) in enumerate(segs)
        r=_rsub2(b,a);usex=abs(r[1])>=abs(r[2])
        q=sort!(collect(splits[si]);by=p->usex ? (p[1]-a[1])/r[1] : (p[2]-a[2])/r[2])
        for k in 1:length(q)-1
            u=getpid(q[k]);v=getpid(q[k+1]);u==v&&continue
            push!(atomic,u<v ? (u,v) : (v,u))
        end
    end
    length(p2)<=nmax||throw(ArgumentError(
        "recover_partition_cdt: planar arrangement has $(length(p2)) nodes, exceeding maxpts=$nmax"))
    xs=Vector{Float64}(undef,length(p2));ys=similar(xs);fmap=Dict{NTuple{2,Float64},NTuple{2,RB}}()
    for i in eachindex(p2)
        q=(Float64(p2[i][1]),Float64(p2[i][2]))
        (isfinite(q[1])&&isfinite(q[2]))||throw(ErrorException(
            "recover_partition_cdt: projected exact intersection is not finite in Float64"))
        haskey(fmap,q)&&fmap[q]!=p2[i]&&throw(ErrorException(
            "recover_partition_cdt: distinct projected intersections collapse at Float64 point $q"))
        fmap[q]=p2[i];xs[i]=q[1];ys[i]=q[2]
    end
    asegs=sort!(collect(atomic));T2=Mesh2D.constrained_delaunay(xs,ys,asegs)
    m2=Mesh2D.to_mesh(T2)
    size(m2.tris,2)>0||throw(ErrorException("recover_partition_cdt: projected arrangement triangulation is empty"))
    pre2=Vector{NTuple{2,RB}}(undef,size(m2.coords,2))
    for i in eachindex(pre2)
        k=(m2.coords[1,i]==0 ? 0.0 : m2.coords[1,i],m2.coords[2,i]==0 ? 0.0 : m2.coords[2,i])
        haskey(fmap,k)||throw(ErrorException("recover_partition_cdt: 2-D triangulation introduced an unmapped point"))
        pre2[i]=fmap[k]
    end
    sort!(unique!(levels));length(levels)>=2||throw(ErrorException(
        "recover_partition_cdt: prismatic surfaces have fewer than two axial levels"))
    nv=size(m2.coords,2);nl=length(levels)
    nall=try Base.checked_mul(nv,nl) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("recover_partition_cdt: prismatic node count overflows Int"))
    end
    nall<=nmax||throw(ArgumentError(
        "recover_partition_cdt: prismatic partition needs $nall nodes, exceeding maxpts=$nmax"))
    preall=Vector{RBv}(undef,nall)
    for k in 1:nl,i in 1:nv
        q=pre2[i];preall[(k-1)*nv+i]=ax==1 ? (levels[k],q[1],q[2]) :
            (ax==2 ? (q[1],levels[k],q[2]) : (q[1],q[2],levels[k]))
    end
    grids=[Mesh3D._raygrid(s) for s in surfaces]
    tv=NTuple{4,Int}[];tags=Int32[];used=Set{Int}()
    potential=try Base.checked_mul(3,Base.checked_mul(size(m2.tris,2),nl-1)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("recover_partition_cdt: prismatic tetrahedron count overflows Int"))
    end
    sizehint!(tv,potential);sizehint!(tags,potential)
    for k in 1:nl-1,ti in axes(m2.tris,2)
        v=sort!(Int[m2.tris[1,ti],m2.tris[2,ti],m2.tris[3,ti]])
        pce=((pre2[v[1]][1]+pre2[v[2]][1]+pre2[v[3]][1])/3,
             (pre2[v[1]][2]+pre2[v[2]][2]+pre2[v[3]][2])/3)
        pcent=(Float64(pce[1]),Float64(pce[2]));amid=Float64((levels[k]+levels[k+1])/2)
        (isfinite(pcent[1])&&isfinite(pcent[2])&&isfinite(amid)) ||
            throw(ErrorException("recover_partition_cdt: a prismatic classification point is not finite in Float64"))
        p=ax==1 ? (amid,pcent[1],pcent[2]) : (ax==2 ? (pcent[1],amid,pcent[2]) : (pcent[1],pcent[2],amid))
        tag=0
        for r in eachindex(grids);Mesh3D._inside_grid(p,grids[r])&&(tag=r;break);end
        tag==0&&continue
        b1=(k-1)*nv+v[1];b2=(k-1)*nv+v[2];b3=(k-1)*nv+v[3]
        t1=k*nv+v[1];t2=k*nv+v[2];t3=k*nv+v[3]
        for q0 in ((b1,b2,b3,t3),(b1,b2,t3,t2),(b1,t2,t3,t1))
            q=q0;s=det3(rsub(preall[q[2]],preall[q[1]]),rsub(preall[q[3]],preall[q[1]]),rsub(preall[q[4]],preall[q[1]]))
            s==0&&throw(ErrorException("recover_partition_cdt: prismatic subdivision produced a flat tetrahedron"))
            s<0&&(q=(q[1],q[2],q[4],q[3]))
            push!(tv,q);push!(tags,Int32(tag));union!(used,q)
        end
    end
    isempty(tv)&&throw(ErrorException("recover_partition_cdt: prismatic partition is empty"))
    uv=sort!(collect(used));nid=Dict{Int,Int32}(v=>Int32(i) for (i,v) in enumerate(uv))
    pre=preall[uv];coords=_float_coords(pre,"recover_partition_cdt")
    tm=Matrix{Int32}(undef,4,length(tv))
    for (j,q) in enumerate(tv),i in 1:4;tm[i,j]=nid[q[i]];end
    out=Mesh(coords;tets=tm,tet_tag=tags);d=validate(out)
    d.ok||throw(ErrorException("recover_partition_cdt: prismatic partition is invalid — "*join(d.messages,"; ")))
    regionpts,facets=dedup_partition(surfaces);regions=build_regions(regionpts,facets)
    ok,reason=_certify_partition_exact(surfaces,out,pre,regions,regionpts)
    ok||throw(ErrorException("recover_partition_cdt: prismatic exact certificate failed — $reason"))
    out
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

"""
    recover_partition_cdt(surfaces; maxiter=2000, maxpts=6000) -> Mesh

Recover every piecewise-linear interface in an ordered collection of closed surfaces,
then classify each tetrahedron by the first surface containing its centroid. Co-axial
prismatic PLCs use an exact planar arrangement and globally conforming extrusion;
other PLCs use exact-preimage conforming-Delaunay refinement. This implements the same
priority semantics as `tetrahedralize_conforming[_exact]`, including overlapping solids,
and adds exact Steiner points wherever the joint vertex-only Delaunay lacks an interface.
Every effective material region is certified as a tetrahedral manifold whose boundary
is a subcomplex of the input PLC collection; otherwise the function throws.
"""
function recover_partition_cdt(surfaces::AbstractVector{Mesh};
                               maxiter::Integer=2000,maxpts::Integer=6000)
    isempty(surfaces) && throw(ArgumentError("recover_partition_cdt: no surfaces"))
    length(surfaces)<=typemax(Int32) ||
        throw(ArgumentError("recover_partition_cdt: region count exceeds Int32 tags"))
    (1<=maxiter<=typemax(Int)) ||
        throw(ArgumentError("recover_partition_cdt: maxiter must be positive and fit Int (got $maxiter)"))
    (4<=maxpts<=typemax(Int32)) ||
        throw(ArgumentError("recover_partition_cdt: maxpts must be in 4:$(typemax(Int32)) (got $maxpts)"))
    for (r,s) in enumerate(surfaces)
        Mesh3D._require_surface3(s,"recover_partition_cdt region $r")
    end
    nit=Int(maxiter);nmax=Int(maxpts)
    prismatic=_recover_prismatic_partition(surfaces,nmax)
    prismatic===nothing||return prismatic
    pts,facets=dedup_partition(surfaces)
    length(pts)>=4 || throw(ArgumentError("recover_partition_cdt: need at least 4 distinct vertices"))
    length(pts)<=nmax || throw(ArgumentError(
        "recover_partition_cdt: input already has $(length(pts)) points, exceeding maxpts=$nmax"))
    isempty(facets) && throw(ArgumentError("recover_partition_cdt: combined PLC has no facets"))
    pointmap=Dict{RBv,Int}(p=>i for (i,p) in enumerate(pts))
    regions=build_regions(pts,facets)
    isempty(regions) && throw(ArgumentError("recover_partition_cdt: combined PLC has no planar regions"))
    seg2regs=Dict{NTuple{2,Int},Vector{Int}}()
    for (ri,reg) in enumerate(regions),e in reg.bnd
        push!(get!(seg2regs,e,Int[]),ri)
    end

    function split_seg!(e)
        a,b=e;mid=rscale(radd(pts[a],pts[b]),1//2);m=get(pointmap,mid,0)
        if m==0
            length(pts)<nmax || error(
                "recover_partition_cdt: maxpts=$nmax reached while splitting a crease")
            push!(pts,mid);m=length(pts);pointmap[mid]=m
        end
        regs=seg2regs[e];delete!(seg2regs,e)
        for e2 in (ekey(a,m),ekey(m,b))
            dst=get!(seg2regs,e2,Int[])
            for ri in regs;ri in dst||push!(dst,ri);end
        end
        for ri in regs
            reg=regions[ri];delete!(reg.bnd,e)
            push!(reg.bnd,ekey(a,m));push!(reg.bnd,ekey(m,b))
        end
        nothing
    end

    for _ in 1:nit
        # Connectivity is built with the production Float64 kernel for bounded
        # memory, while every recovery and acceptance decision below uses the exact
        # Rational preimages.  The conversion gate rejects point collisions, and the
        # exact certificate—not the approximate Delaunay property—is authoritative.
        etets=_partition_tets_float(pts)
        E=Set{NTuple{2,Int}}()
        for t in etets,i in 1:4,j in i+1:4;push!(E,ekey(t[i],t[j]));end

        missing=NTuple{2,Int}[e for e in keys(seg2regs) if !(e in E)]
        if !isempty(missing)
            for e in missing;haskey(seg2regs,e)&&split_seg!(e);end
            continue
        end

        candidates=RBv[]
        for reg in regions
            for e in E
                pierced,y=edge_pierces(reg,pts,e[1],e[2])
                if pierced;push!(candidates,y);break;end
            end
        end
        if !isempty(candidates)
            fresh=RBv[y for y in unique(candidates) if !haskey(pointmap,y)]
            isempty(fresh) && error(
                "recover_partition_cdt: an interface remains pierced but every exact intersection point already exists")
            length(fresh)<=nmax-length(pts) || error(
                "recover_partition_cdt: adding interface points would exceed maxpts=$nmax")
            for y in fresh;push!(pts,y);pointmap[y]=length(pts);end
            continue
        end

        return _finalize_partition(surfaces,pts,etets,regions)
    end
    error("recover_partition_cdt: no convergence in $nit iterations (verts=$(length(pts)))")
end

# internal: run the refinement and return (mesh, pts, regions, seg2regs, n0) so sizing can reuse it.
function _recover(surface::Mesh; maxiter::Integer=2000, maxpts::Integer=6000)
    Mesh3D._require_surface3(surface,"recover_boundary_cdt")
    (1 <= maxiter <= typemax(Int)) ||
        throw(ArgumentError("recover_boundary_cdt: maxiter must be positive and fit Int (got $maxiter)"))
    (4 <= maxpts <= typemax(Int32)) ||
        throw(ArgumentError("recover_boundary_cdt: maxpts must be in 4:$(typemax(Int32)) (got $maxpts)"))
    nit=Int(maxiter);nmax=Int(maxpts)
    pts, facets = dedup(surface)
    length(pts)>=4 || throw(ArgumentError("recover_boundary_cdt: need ≥ 4 distinct vertices"))
    isempty(facets) && throw(ArgumentError("recover_boundary_cdt: surface has no facets"))
    n0=length(pts)
    n0<=nmax || throw(ArgumentError("recover_boundary_cdt: input already has $n0 points, exceeding maxpts=$nmax"))
    pointmap=Dict{RBv,Int}(p=>i for (i,p) in enumerate(pts))
    regions=build_regions(pts,facets)
    seg2regs=Dict{NTuple{2,Int},Vector{Int}}()
    for (ri,reg) in enumerate(regions), e in reg.bnd; push!(get!(seg2regs,e,Int[]),ri); end

    function split_seg!(e)
        a,b=e;mid=rscale(radd(pts[a],pts[b]),1//2);m=get(pointmap,mid,0)
        if m==0
            length(pts)<nmax || error("recover_boundary_cdt: maxpts=$nmax reached while splitting a crease")
            push!(pts,mid);m=length(pts);pointmap[mid]=m
        end
        regs=seg2regs[e]; delete!(seg2regs,e)
        for e2 in (ekey(a,m),ekey(m,b))
            dst=get!(seg2regs,e2,Int[])
            for ri in regs;ri in dst||push!(dst,ri);end
        end
        for ri in regs; reg=regions[ri]; delete!(reg.bnd,e); push!(reg.bnd,ekey(a,m)); push!(reg.bnd,ekey(m,b)); end
        m
    end

    local etets
    for it in 1:nit
        length(pts) > nmax && error("recover_boundary_cdt: maxpts exceeded ($(length(pts))) — refinement not converging")
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
            fresh=RBv[y for y in unique(newpts) if !haskey(pointmap,y)]
            isempty(fresh) && error(
                "recover_boundary_cdt: a facet remains pierced but every exact intersection point already exists")
            length(fresh)<=nmax-length(pts) ||
                error("recover_boundary_cdt: adding piercing points would exceed maxpts=$nmax")
            for y in fresh;push!(pts,y);pointmap[y]=length(pts);end
            continue
        end
        # converged: assemble + EXACT certification (never return a silently non-conforming mesh)
        m = _finalize(surface, pts, etets)
        ok, reason = certify_exact(surface, m, pts, regions, seg2regs)
        ok && return (m, pts, regions, seg2regs, n0)
        error("recover_boundary_cdt: convergence gate passed but exact certification failed: $reason")
    end
    error("recover_boundary_cdt: no convergence in $nit iters (verts=$(length(pts)))")
end

"""
    mesh_sized_cdt(surface::Mesh; hmax) -> Mesh

Uniform size control on an **arbitrary** (curved / non-star) domain, on top of the
exact conforming-Delaunay recovery. First recover the boundary conformingly
([`recover_boundary_cdt`](@ref)), then apply conforming longest-edge bisection to
every overlong tet edge. Boundary faces are subdivided on their own planes, so the
piecewise-linear PLC is unchanged. The returned mesh is exactly PLC-conforming and
has **every edge `≤ hmax`**; failure to prove either postcondition is an explicit
blocker, never a silent unsized fallback.
"""
function mesh_sized_cdt(surface::Mesh; hmax::Real)
    hm=try Float64(hmax) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_cdt: hmax must be Float64-representable: $(sprint(showerror,err))"))
    end
    (isfinite(hm)&&hm>0) || throw(ArgumentError("mesh_sized_cdt: hmax must be finite and positive (got $hmax)"))
    m0,_,_,_,_=_recover(surface)
    m=Mesh3D.refine_to_size(m0,hm)
    ok,reason=Mesh3D._certify_surface_fill(surface,m)
    ok || throw(ErrorException("mesh_sized_cdt: refined mesh failed the exact PLC certificate — $reason"))
    return m
end

end # module RecoverCDT
