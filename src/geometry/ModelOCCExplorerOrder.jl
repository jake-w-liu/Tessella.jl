# ── `TopExp_Explorer` order of a fillet or chamfer result ────────────────────
#
# Gmsh's `_multiBind` numbers a filleted or chamfered solid's faces, edges and
# vertices in the order `TopExp_Explorer` meets them, and that order is the one
# OpenCASCADE's builder (`ChFi3d_Builder` through `TopOpeBRepBuild`) leaves:
#
# - faces: a breadth-first walk of the shell's faces, starting at the first
#   modified face of the old solid and taking each face's neighbours across
#   its edges in explorer order;
# - a rebuilt face's edges: a breadth-first walk through shared vertices,
#   starting at the first edge of the old face that the blend cut short (in
#   the old face's explorer order; an uncut edge when none was cut); a face
#   the blend did not touch keeps its old order;
# - a cut edge lists its new vertex first (two new vertices: in the old
#   edge's direction); an edge along a blend lists its vertices in the
#   direction of the blended edge, a blend face's end edge from the vertex on
#   the first face of the blended edge;
# - a blend face starts at its end edge at the blended edge's start, then its
#   edges on the two faces of the blended edge.
#
# Checked against Gmsh 4.15.2 on every single and pairwise chamfer and fillet
# of a box's edges (78 cases each) plus the three- and four-edge chains of
# its faces: face, curve and point tags all equal.

_occ_xp_close(a,b)=_occ_modulus(_occ3_sub(a,b))<=1e-7

function _occ_xp_faces(kernel,vid)
    out=UInt64[]
    for sh in kernel.shapes[vid].children
        g=kernel.shapes[sh].geometry
        hasproperty(g,:uses) || continue
        for (f,_) in g.uses
            f in out || push!(out,f)
        end
    end
    return out
end

function _occ_xp_edges(kernel,fid)
    e=get(kernel.explorer,fid,nothing)
    e===nothing || return e
    out=UInt64[]
    for w in kernel.shapes[fid].children
        g=kernel.shapes[w].geometry
        hasproperty(g,:uses) || continue
        for (u,_) in g.uses
            u in out || push!(out,u)
        end
    end
    return out
end

_occ_xp_vertices(kernel,eid)=
    unique(get(kernel.explorer,eid,kernel.shapes[eid].children))

function _occ_xp_plane(kernel,fid)
    g=kernel.shapes[fid].geometry
    hasproperty(g,:geometry) || return nothing
    get(g.geometry,:occ,nothing)===:plane || return nothing
    n=_occ3_scale(g.geometry.axis,1.0/_occ_modulus(g.geometry.axis))
    return (n,_dot3(n,g.geometry.center))
end

_occ_xp_same_plane(p,q)=p!==nothing && q!==nothing &&
    ((_dot3(p[1],q[1])>1-1e-9 && abs(p[2]-q[2])<=1e-7) ||
     (_dot3(p[1],q[1])< -1+1e-9 && abs(p[2]+q[2])<=1e-7))

# What the rebuild needs of the solid before the blend: its explorer orders,
# its planes and its vertex positions (the old shapes may be released).
struct _OCCExplorerSnapshot
    faces::Vector{UInt64}
    fedges::Dict{UInt64,Vector{UInt64}}
    everts::Dict{UInt64,Vector{UInt64}}
    efaces::Dict{UInt64,Vector{UInt64}}
    estart::Dict{UInt64,UInt64}
    vpos::Dict{UInt64,NTuple{3,Float64}}
    plane::Dict{UInt64,Any}
    pgeom::Dict{UInt64,Any}
    frev::Dict{UInt64,Bool}
end

# A blend face: the old edge it blends and its edges along the edge's first
# and second face.
struct _OCCBlendSource
    src::UInt64
    S1::Vector{UInt64}
    S2::Vector{UInt64}
end

function _occ_explorer_snapshot(kernel,vid)
    faces=_occ_xp_faces(kernel,vid)
    fedges=Dict{UInt64,Vector{UInt64}}()
    everts=Dict{UInt64,Vector{UInt64}}()
    efaces=Dict{UInt64,Vector{UInt64}}()
    estart=Dict{UInt64,UInt64}()
    vpos=Dict{UInt64,NTuple{3,Float64}}()
    plane=Dict{UInt64,Any}()
    pgeom=Dict{UInt64,Any}()
    frev=Dict{UInt64,Bool}()
    for sh in kernel.shapes[vid].children
        g=kernel.shapes[sh].geometry
        hasproperty(g,:uses) || continue
        for (f,r) in g.uses
            frev[f]=r
        end
    end
    for f in faces
        plane[f]=_occ_xp_plane(kernel,f)
        plane[f]===nothing || (pgeom[f]=kernel.shapes[f].geometry.geometry)
        fedges[f]=copy(_occ_xp_edges(kernel,f))
        for e in fedges[f]
            push!(get!(efaces,e,UInt64[]),f)
            haskey(everts,e) && continue
            everts[e]=copy(_occ_xp_vertices(kernel,e))
            estart[e]=first(kernel.shapes[e].children)
            for v in kernel.shapes[e].children
                vpos[v]=kernel.shapes[v].geometry.point
            end
        end
    end
    return _OCCExplorerSnapshot(faces,fedges,everts,efaces,estart,vpos,plane,
                                pgeom,frev)
end

# Order the blended solid `vid` as OpenCASCADE's builder does, from the solid
# `old` was taken of and the blended edges `blended` (old edge ids). Records
# the new explorer orders, puts the shell's faces in that order, starts each
# rebuilt wire at its face's first edge, and returns the shapes in the order
# `_multiBind` numbers them (for `_occ_retag_result!`).
function _occ_explorer_rebuild!(kernel,old::_OCCExplorerSnapshot,vid,blended)
    faces=_occ_xp_faces(kernel,vid)
    oldset=Set(old.faces)
    fmatch=Dict{UInt64,Union{Nothing,UInt64}}()
    for f in faces
        if f in oldset
            fmatch[f]=f
        else
            p=_occ_xp_plane(kernel,f)
            fmatch[f]=p===nothing ? nothing :
                _occ_xp_find(g->_occ_xp_same_plane(p,old.plane[g]),
                                     old.faces)
        end
    end
    wedges=Dict{UInt64,Vector{UInt64}}(f=>_occ_xp_edges(kernel,f) for f in faces)
    efaces=Dict{UInt64,Vector{UInt64}}()
    for f in faces, e in wedges[f]
        push!(get!(efaces,e,UInt64[]),f)
    end
    pos(v)=kernel.shapes[v].geometry.point
    verts(e)=unique(kernel.shapes[e].children)
    # new edge => old edge it is (a part of)
    ematch=Dict{UInt64,Union{Nothing,UInt64}}()
    for (e,fs) in efaces
        if haskey(old.everts,e)
            ematch[e]=e; continue
        end
        gs=[fmatch[f] for f in fs]
        ematch[e]=nothing
        (length(gs)==2 && all(!isnothing,gs)) || continue
        for oe in keys(old.efaces)
            Set(old.efaces[oe])==Set(gs) || continue
            ps=[old.vpos[v] for v in old.everts[oe]]
            length(ps)==2 || continue
            a,b=ps; d=_occ3_sub(b,a); L2=_dot3(d,d)
            ok=all(verts(e)) do v
                q=_occ3_sub(pos(v),a); t=_dot3(q,d)/L2
                -1e-9<=t<=1+1e-9 &&
                    _occ_modulus(_occ3_sub(q,_occ3_scale(d,t)))<=1e-7
            end
            ok && (ematch[e]=oe; break)
        end
    end
    isold(e,v)=(oe=ematch[e]; oe!==nothing &&
        any(w->_occ_xp_close(pos(v),old.vpos[w]),old.everts[oe]))
    unsplit(e)=(oe=ematch[e]; oe!==nothing &&
        length(verts(e))==length(old.everts[oe]) && all(v->isold(e,v),verts(e)))
    vorder=Dict{UInt64,Vector{UInt64}}()
    for e in keys(efaces)
        oe=ematch[e]; oe===nothing && continue
        vs=verts(e)
        if unsplit(e)
            vorder[e]=[only(v for v in vs if _occ_xp_close(pos(v),old.vpos[w]))
                       for w in old.everts[oe]]
        else
            olds=[v for v in vs if isold(e,v)]
            news=[v for v in vs if !isold(e,v)]
            s=old.vpos[old.estart[oe]]
            sort!(news,by=v->_occ_modulus(_occ3_sub(pos(v),s)))
            vorder[e]=vcat(news,olds)
        end
    end
    # blend faces: their blended edge, end edges and edges along the blend
    forder=Dict{UInt64,Vector{UInt64}}()
    bsrc=Dict{UInt64,_OCCBlendSource}()
    for f in faces
        fmatch[f]===nothing || continue
        nb=Dict(e=>[fmatch[h] for h in efaces[e] if h!=f] for e in wedges[f])
        near=Set(g for gs in values(nb) for g in gs if g!==nothing)
        src=_occ_xp_find(oe->haskey(old.efaces,oe) &&
            all(g->g in near,old.efaces[oe]),blended)
        if src===nothing
            forder[f]=wedges[f]; continue
        end
        g1,g2=sort(old.efaces[src],by=g->findfirst(==(g),old.faces))
        s=old.vpos[old.estart[src]]
        along(e,g)=g in nb[e]
        S1=[e for e in wedges[f] if along(e,g1)]
        S2=[e for e in wedges[f] if along(e,g2) && !(e in S1)]
        bsrc[f]=_OCCBlendSource(src,S1,S2)
        ends=[e for e in wedges[f] if !(e in S1) && !(e in S2)]
        # which end of the blended edge an end edge sits at, and whether it
        # meets another blend there (a corner) or a face of the solid
        se=old.vpos[only(v for v in old.everts[src] if v!=old.estart[src])]
        dvec=_occ3_sub(se,s)
        atmax(e)=(c=reduce(_occ3_add,(pos(v) for v in verts(e)));
                  _dot3(_occ3_sub(_occ3_scale(c,1/length(verts(e))),s),dvec)>
                      0.5*_dot3(dvec,dvec))
        blend_end(e)=any(isnothing,nb[e])
        cmax=any(e->blend_end(e) && atmax(e),ends)
        cmin=any(e->blend_end(e) && !atmax(e),ends)
        function erank(e)
            b=blend_end(e); mx=atmax(e)
            !b && mx && cmax && return 0     # face cap at a corner, far end
            !b && !mx && !cmin && return 1   # face cap at the free start
            b && !mx && return 2
            !b && !mx && return 3
            b && mx && return 4
            return 5
        end
        sort!(ends,by=erank)
        forder[f]=isempty(ends) ? vcat(S1,S2) :
            vcat(ends[1:1],S1,S2,ends[2:end])
        for e in vcat(S1,S2)
            haskey(vorder,e) && continue
            vorder[e]=sort(verts(e),by=v->_occ_modulus(_occ3_sub(pos(v),s)))
        end
        s1v=Set(v for e in S1 for v in verts(e))
        # `ChFi3d` trims the end curves of a blend: Gmsh types them
        # `TrimmedCurve` (a line, or an arc of a fillet), not Line or Circle
        for e in ends
            g=kernel.shapes[e].geometry
            hasproperty(g.geometry,:native_type) && continue
            kernel.shapes[e]=_ModelOCCShape(1,merge(g,(geometry=merge(
                g.geometry,(native_type="TrimmedCurve",)),)),
                kernel.shapes[e].children)
        end
        for e in ends
            haskey(vorder,e) && continue
            # between two blends: the lower point first (fits every
            # two-blend corner of the box checked)
            blend_end(e) && (vorder[e]=sort(verts(e),
                by=v->map(x->round(x;digits=6)+0.0,pos(v))); continue)
            vorder[e]=sort(verts(e),by=v->v in s1v ? 0 : 1)
        end
    end
    for e in keys(efaces)
        haskey(vorder,e) || (vorder[e]=verts(e))
    end
    modified=Dict{UInt64,Bool}()
    for f in faces
        g=fmatch[f]
        if g===nothing
            modified[f]=true
        else
            es=wedges[f]
            modified[f]=!(length(es)==length(old.fedges[g]) &&
                all(e->unsplit(e) && ematch[e] in old.fedges[g],es))
        end
    end
    for f in faces
        g=fmatch[f]; g===nothing && continue
        es=wedges[f]
        if !modified[f]
            forder[f]=[only(e for e in es if ematch[e]==oe)
                       for oe in old.fedges[g]]
            continue
        end
        kept=UInt64[]
        for oe in old.fedges[g], e in es
            ematch[e]==oe && !unsplit(e) && !(e in kept) && push!(kept,e)
        end
        for oe in old.fedges[g], e in es
            ematch[e]==oe && unsplit(e) && !(e in kept) && push!(kept,e)
        end
        elems=vcat(kept,[e for e in es if !(e in kept)])
        forder[f]=_occ_xp_bfs(elems,e->vorder[e],e->verts(e))
    end
    for f in faces
        fmatch[f]===nothing || continue
        forder[f]=_occ_xp_bfs(forder[f],e->vorder[e],e->verts(e))
    end
    # faces: breadth first from the first modified old face
    rank(f)=(g=fmatch[f]; g===nothing ? typemax(Int) :
             findfirst(==(g),old.faces))
    cands=[f for f in faces if modified[f] && fmatch[f]!==nothing]
    start=isempty(cands) ? first(faces) : cands[argmin(map(rank,cands))]
    fq=UInt64[start]; i=1
    while i<=length(fq)
        for e in forder[fq[i]], h in efaces[e]
            h in fq || push!(fq,h)
        end
        i+=1
    end
    for f in faces
        f in fq || push!(fq,f)
    end
    order=UInt64[]
    seen=Set{UInt64}()
    add!(x)=(x in seen || (push!(seen,x); push!(order,x)))
    for f in fq
        add!(f)
    end
    for f in fq, e in forder[f]
        add!(e)
    end
    for f in fq, e in forder[f], v in vorder[e]
        add!(v)
    end
    for f in faces
        kernel.explorer[f]=forder[f]
    end
    for e in keys(efaces)
        kernel.explorer[e]=vorder[e]
    end
    frev_new=_occ_xp_orient!(kernel,old,faces,oldset,fmatch,ematch,efaces,
                             vorder,bsrc)
    # the shell in explorer order; each rebuilt wire from its first edge
    for sh in kernel.shapes[vid].children
        g=kernel.shapes[sh].geometry
        hasproperty(g,:uses) || continue
        uses=sort([(f,get(frev_new,f,r)) for (f,r) in g.uses],
                  by=u->findfirst(==(u[1]),fq))
        kernel.shapes[sh]=_ModelOCCShape(-2,(uses=uses,),
            UInt64[u[1] for u in uses])
    end
    for f in faces
        f in oldset && continue
        ws=kernel.shapes[f].children
        isempty(ws) && continue
        w=first(ws); wg=kernel.shapes[w].geometry
        hasproperty(wg,:uses) || continue
        k=findfirst(u->u[1]==first(forder[f]),wg.uses)
        (k===nothing || k==1) && continue
        uses=vcat(wg.uses[k:end],wg.uses[1:k-1])
        kernel.shapes[w]=_ModelOCCShape(-1,(uses=uses,),
            UInt64[u[1] for u in uses])
    end
    return order
end

# Edge directions, plane frames and face orientations as OpenCASCADE's
# builder leaves them, where every face of an edge is a plane (a fillet's
# cylinder keeps the pcurves its construction gave it):
# - an edge on an edge of the old solid runs the old edge's way; an edge
#   along a blend or across its end runs from its first explorer vertex;
# - a cut plane keeps the old face's frame and orientation;
# - a chamfer plane has its origin between the chamfer's two lines at the
#   blended edge's start, X along that edge, Y from the line on the first
#   face to the line on the second, and is used reversed when its normal
#   points into the solid;
# - every wire turns counterclockwise about its plane's normal.
# Returns the new shell orientation of each reoriented face.
function _occ_xp_orient!(kernel,old::_OCCExplorerSnapshot,faces,oldset,fmatch,
                         ematch,efaces,vorder,bsrc)
    caller="OCC.fillet"
    pos(v)=kernel.shapes[v].geometry.point
    verts(e)=unique(kernel.shapes[e].children)
    isplane(h)=_occ_xp_plane(kernel,h)!==nothing
    for e in keys(efaces)
        sh=kernel.shapes[e]; g=sh.geometry.geometry
        get(g,:occ,nothing)===:line || continue
        all(isplane,efaces[e]) || continue
        length(unique(sh.children))==2 || continue
        oe=ematch[e]
        want=if oe===nothing
            first(vorder[e])
        else
            s=old.vpos[old.estart[oe]]
            vs=verts(e)
            vs[argmin([_occ_modulus(_occ3_sub(pos(v),s)) for v in vs])]
        end
        first(sh.children)==want && continue
        a,b=sh.children
        g2=merge(g,(origin=_occ3_add(g.origin,_occ3_scale(g.dir,g.t0+g.t1)),
                    dir=_occ3_scale(g.dir,-1.0)))
        kernel.shapes[e]=_ModelOCCShape(1,merge(sh.geometry,(geometry=g2,)),
                                        UInt64[b,a])
        for h in efaces[e], w in kernel.shapes[h].children
            wg=kernel.shapes[w].geometry
            hasproperty(wg,:uses) || continue
            any(u->u[1]==e,wg.uses) || continue
            kernel.shapes[w]=_ModelOCCShape(-1,
                (uses=Tuple{UInt64,Bool}[(u,u==e ? !r : r) for (u,r) in wg.uses],),
                kernel.shapes[w].children)
        end
    end
    frev_new=Dict{UInt64,Bool}()
    for f in faces
        f in oldset && continue
        sh=kernel.shapes[f]; fg=sh.geometry
        get(fg.geometry,:occ,nothing)===:plane || continue
        g=fmatch[f]
        if g!==nothing
            haskey(old.pgeom,g) || continue
            og=old.pgeom[g]
            center,axis,X,Y=og.center,og.axis,og.X,og.Y
            rv=get(og,:reversed,false); r=old.frev[g]
        elseif haskey(bsrc,f) && !isempty(bsrc[f].S1) && !isempty(bsrc[f].S2)
            b=bsrc[f]
            s=old.vpos[old.estart[b.src]]
            se=old.vpos[only(v for v in old.everts[b.src]
                             if v!=old.estart[b.src])]
            U=_occ_dir(_occ3_sub(se,s))
            foot(e)=(p=pos(first(verts(e)));
                     _occ3_sub(p,_occ3_scale(U,_dot3(_occ3_sub(p,s),U))))
            A0=foot(first(b.S1)); B0=foot(first(b.S2))
            center=_occ3_scale(_occ3_add(A0,B0),0.5)
            X=U; Y=_occ_dir(_occ3_sub(B0,A0))
            axis=_occ_cross(X,Y)
            r=_dot3(axis,_occ3_sub(s,center))<0; rv=r
        else
            continue
        end
        pcurves=Dict{Int,NamedTuple}()
        umin=vmin=Inf; umax=vmax=-Inf
        for w in sh.children, (u,_) in kernel.shapes[w].geometry.uses
            eg=kernel.shapes[u].geometry.geometry
            pc=_model_occ_plane_pcurve(eg,center,X,Y)
            eg.occ in (:circle,:ellipse) && _dot3(eg.n,axis)<0 &&
                (pc=(neg2d=(pc,),))
            pcurves[kernel.shape_tags[u][2]]=(fwd=pc,rev=nothing)
            lo,hi=_occ_pcurve_fillbndbox(pc,eg.t0,eg.t1)
            umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
            umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
        end
        newg=merge(fg.geometry,(center=center,axis=axis,X=X,Y=Y,reversed=rv,
                                pcurves=pcurves,uvb=([umin,vmin],[umax,vmax])))
        kernel.shapes[f]=_ModelOCCShape(2,merge(fg,(geometry=newg,)),sh.children)
        frev_new[f]=r
        for w in sh.children
            wg=kernel.shapes[w].geometry
            pts=_model_occ_wire_samples(id->kernel.shapes[id],wg.uses,caller)
            area=0.0
            for i in eachindex(pts)
                p=_occ3_sub(pts[i],center); q=_occ3_sub(pts[mod1(i+1,length(pts))],center)
                area+=_dot3(p,X)*_dot3(q,Y)-_dot3(q,X)*_dot3(p,Y)
            end
            area<0 || continue
            uses=Tuple{UInt64,Bool}[(u,!r2) for (u,r2) in reverse(wg.uses)]
            kernel.shapes[w]=_ModelOCCShape(-1,(uses=uses,),
                UInt64[u[1] for u in uses])
        end
    end
    return frev_new
end

# A `.geo` Box bound into the kernel later (its faces' wires and edges are
# `add_box!`'s, which are OpenCASCADE's): the explorer order of
# `BRepPrim_GWedge`, as `_model_occ_add_box!` records it.
function _occ_explorer_box!(kernel,vid)
    shells=kernel.shapes[vid].children
    length(shells)==1 || return nothing
    faces=_occ_xp_faces(kernel,vid)
    length(faces)==6 || return nothing
    for f in faces
        ws=kernel.shapes[f].children
        (length(ws)==1 && length(kernel.shapes[ws[1]].children)==4) ||
            return nothing
    end
    for f in faces
        us=kernel.shapes[kernel.shapes[f].children[1]].geometry.uses
        kernel.explorer[f]=UInt64[us[j][1] for j in (1,4,3,2)]
        for (e,_) in us
            kernel.explorer[e]=reverse(kernel.shapes[e].children)
        end
    end
    return nothing
end

_occ_xp_find(p,xs)=(i=findfirst(p,xs); i===nothing ? nothing : xs[i])

# Breadth-first walk of a face's edges through their vertices (the order of
# `TopOpeBRepBuild_BlockBuilder`), from the first of `elems`.
function _occ_xp_bfs(elems,vord,vall)
    isempty(elems) && return elems
    q=UInt64[first(elems)]; i=1
    while i<=length(q)
        for v in vord(q[i]), e in elems
            e in q && continue
            v in vall(e) && push!(q,e)
        end
        i+=1
    end
    for e in elems
        e in q || push!(q,e)
    end
    return q
end
