# ── Geometric independence certificate ─────────────────────────────────────
# Strictly separated bounding boxes already prove independence. For pairs whose
# boxes overlap, this certificate decides conforming contact with exact
# predicates: coincident facet pairs with opposed winding weld into shared
# interface Surfaces, coincident vertex/edge-only contacts are admitted
# unwelded (each solid keeps its own coincident boundary node, matching OCC
# import semantics), and every other contact — partial face overlap,
# T-junctions, off-edge segments — is still rejected, as is interpenetration
# or enclosure of one shell inside the other.

const _BREP_R3=NTuple{3,Rational{BigInt}}
@inline _brep_rpoint(p)::_BREP_R3=ntuple(d->Rational{BigInt}(p[d]),3)
@inline _brep_rsub(a::_BREP_R3,b::_BREP_R3)::_BREP_R3=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _brep_rdot(a::_BREP_R3,b::_BREP_R3)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline function _brep_rcross(a::_BREP_R3,b::_BREP_R3)::_BREP_R3
    return (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
end
# Exact orient3 determinant det[a-d,b-d,c-d] over rationalized coordinates.
@inline function _brep_rorient3(a::_BREP_R3,b::_BREP_R3,c::_BREP_R3,d::_BREP_R3)
    return _brep_rdot(_brep_rsub(a,d),_brep_rcross(_brep_rsub(b,d),_brep_rsub(c,d)))
end

@inline function _brep_facet_box(T)
    return ((min(T[1][1],T[2][1],T[3][1]),min(T[1][2],T[2][2],T[3][2]),min(T[1][3],T[2][3],T[3][3])),
            (max(T[1][1],T[2][1],T[3][1]),max(T[1][2],T[2][2],T[3][2]),max(T[1][3],T[2][3],T[3][3])))
end
@inline function _brep_boxes_touch(a,b)
    return !(a[2][1]<b[1][1] || b[2][1]<a[1][1] ||
             a[2][2]<b[1][2] || b[2][2]<a[1][2] ||
             a[2][3]<b[1][3] || b[2][3]<a[1][3])
end
@inline function _brep_box_has(box,p)
    return box[1][1]<=p[1]<=box[2][1] && box[1][2]<=p[2]<=box[2][2] &&
           box[1][3]<=p[3]<=box[2][3]
end

function _brep_shell_facets(surface::Mesh)
    return [ntuple(r->node(surface,surface.tris[r,column]),3)
            for column in 1:ntris(surface)]
end

# 2-D projection helper for coplanar decisions: drop coordinate axis k.
@inline _brep_drop(p,k)=(p[mod1(k+1,3)],p[mod1(k+2,3)])
function _brep_coplanar_drop_axis(T)
    a,b,c=T
    n=((b[2]-a[2])*(c[3]-a[3])-(b[3]-a[3])*(c[2]-a[2]),
       (b[3]-a[3])*(c[1]-a[1])-(b[1]-a[1])*(c[3]-a[3]),
       (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))
    (n[1]==0.0&&n[2]==0.0&&n[3]==0.0) || return argmax(abs.(collect(n)))
    # Degenerate triangle: keep its widest 2-D projection.
    spans=[max(T[1][k],T[2][k],T[3][k])-min(T[1][k],T[2][k],T[3][k]) for k in 1:3]
    return argmin(spans)
end

@inline function _brep_on_seg_2d(a,b,p)
    return min(a[1],b[1])<=p[1]<=max(a[1],b[1]) && min(a[2],b[2])<=p[2]<=max(a[2],b[2])
end
function _brep_segs_intersect_2d(a,b,c,d)
    s1=orient2(a,b,c);s2=orient2(a,b,d)
    s3=orient2(c,d,a);s4=orient2(c,d,b)
    s1*s2<0 && s3*s4<0 && return true
    s1==0 && _brep_on_seg_2d(a,b,c) && return true
    s2==0 && _brep_on_seg_2d(a,b,d) && return true
    s3==0 && _brep_on_seg_2d(c,d,a) && return true
    s4==0 && _brep_on_seg_2d(c,d,b) && return true
    return false
end
# Closed-triangle containment in the plane; false for degenerate triangles
# (their edges carry every real contact through the segment tests).
function _brep_point_in_tri_2d(p,T)
    orient2(T[1],T[2],T[3])==0 && return false
    s1=orient2(T[1],T[2],p);s2=orient2(T[2],T[3],p);s3=orient2(T[3],T[1],p)
    return (s1>=0&&s2>=0&&s3>=0)||(s1<=0&&s2<=0&&s3<=0)
end
# Is point p on the closed triangle T? (T2 is T's dominant-axis projection.)
@inline function _brep_point_on_tri(p,T,T2,drop)
    orient3(T[1],T[2],T[3],p)!=0 && return false
    return _brep_point_in_tri_2d(_brep_drop(p,drop),T2)
end

@inline _brep_rdrop(p,k)=(p[mod1(k+1,3)],p[mod1(k+2,3)])
@inline function _brep_rorient2(a,b,c)
    return sign((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))
end
# Closed-triangle containment on rational coordinates.
function _brep_rpoint_in_tri_2d(p,T)
    _brep_rorient2(T[1],T[2],T[3])==0 && return false
    s1=_brep_rorient2(T[1],T[2],p);s2=_brep_rorient2(T[2],T[3],p)
    s3=_brep_rorient2(T[3],T[1],p)
    return (s1>=0&&s2>=0&&s3>=0)||(s1<=0&&s2<=0&&s3<=0)
end

# Where edge (a,b) pierces the plane of U, kept in exact rationals.
function _brep_edge_plane_pierce(Ur,ra,rb)
    oa=_brep_rorient3(Ur[1],Ur[2],Ur[3],ra)
    ob=_brep_rorient3(Ur[1],Ur[2],Ur[3],rb)
    λ=oa/(oa-ob)
    return ntuple(d->ra[d]+(rb[d]-ra[d])*λ,3)
end

# Where two closed 3-D segments cross at a point interior to both, in exact
# rationals. Endpoint and collinear touches return nothing here.
function _brep_segs_cross_point_3d(a,b,c,d)
    ra,rb,rc,rd=_brep_rpoint(a),_brep_rpoint(b),_brep_rpoint(c),_brep_rpoint(d)
    e1=_brep_rsub(rb,ra);e2=_brep_rsub(rd,rc)
    rhs=_brep_rsub(rc,ra)
    for (i,j) in ((1,2),(1,3),(2,3))
        det=e1[i]*(-e2[j])-e1[j]*(-e2[i])
        det==0 && continue
        μ=(rhs[i]*(-e2[j])-rhs[j]*(-e2[i]))/det
        ν=(e1[i]*rhs[j]-e1[j]*rhs[i])/det
        k=6-i-j
        ra[k]+μ*e1[k]==rc[k]+ν*e2[k] || return nothing
        return (0<μ<1 && 0<ν<1) ? ntuple(x->ra[x]+μ*e1[x],3) : nothing
    end
    return nothing
end
# Do two closed 3-D segments cross at a point interior to both? Endpoint and
# collinear touches are handled by the vertex checks and return false here.
_brep_segs_cross_3d(a,b,c,d)=_brep_segs_cross_point_3d(a,b,c,d)!==nothing

# Is rational point p on the closed rational segment (a,b)?
function _brep_rpoint_on_rseg(p,a,b)
    ab=_brep_rsub(b,a);ap=_brep_rsub(p,a)
    c=_brep_rcross(ab,ap)
    (c[1]==0&&c[2]==0&&c[3]==0) || return false
    t=_brep_rdot(ap,ab);d=_brep_rdot(ab,ab)
    return 0<=t<=d
end

# The contact set of two closed triangles: (:none), (:point,[v]),
# (:segment,[v1,v2]) or :bad for interior crossings, edge/pierce contacts and
# positive-area overlap — anything that cannot be a rim touch.
function _brep_tris_contact(T,U)
    s=ntuple(i->orient3(U[1],U[2],U[3],T[i]),3)
    all(>(0),s)&&return (:none,nothing)
    all(<(0),s)&&return (:none,nothing)
    t=ntuple(j->orient3(T[1],T[2],T[3],U[j]),3)
    all(>(0),t)&&return (:none,nothing)
    all(<(0),t)&&return (:none,nothing)
    dropT=_brep_coplanar_drop_axis(T);dropU=_brep_coplanar_drop_axis(U)
    T2=ntuple(i->_brep_drop(T[i],dropT),3)
    U2=ntuple(j->_brep_drop(U[j],dropU),3)
    pts=Set{NTuple{3,Float64}}()
    for i in 1:3
        s[i]==0 && _brep_point_on_tri(T[i],U,U2,dropU) && push!(pts,T[i])
        t[i]==0 && _brep_point_on_tri(U[i],T,T2,dropT) && push!(pts,U[i])
    end
    if all(==(0),s)
        _brep_tris_overlap_area_2d(T2,U2) && return (:bad,nothing)
        for i in 1:3,j in 1:3
            a,b=T2[i],T2[mod1(i+1,3)];c,d=U2[j],U2[mod1(j+1,3)]
            orient2(a,b,c)*orient2(a,b,d)<0 && orient2(c,d,a)*orient2(c,d,b)<0 &&
                return (:bad,nothing)
        end
    else
        Tr=ntuple(i->_brep_rpoint(T[i]),3)
        Ur=ntuple(j->_brep_rpoint(U[j]),3)
        U2r=ntuple(j->_brep_rdrop(Ur[j],dropU),3)
        T2r=ntuple(i->_brep_rdrop(Tr[i],dropT),3)
        for i in 1:3
            j=i==3 ? 1 : i+1
            if s[i]*s[j]<0
                p=_brep_edge_plane_pierce(Ur,Tr[i],Tr[j])
                _brep_rpoint_in_tri_2d(_brep_rdrop(p,dropU),U2r) &&
                    return (:bad,nothing)
            end
            if t[i]*t[j]<0
                p=_brep_edge_plane_pierce(Tr,Ur[i],Ur[j])
                _brep_rpoint_in_tri_2d(_brep_rdrop(p,dropT),T2r) &&
                    return (:bad,nothing)
            end
        end
        for i in 1:3,j in 1:3
            _brep_segs_cross_3d(T[i],T[mod1(i+1,3)],U[j],U[mod1(j+1,3)]) &&
                return (:bad,nothing)
        end
    end
    length(pts)>2 && return (:bad,nothing)
    length(pts)==2 && return (:segment,collect(pts))
    length(pts)==1 && return (:point,collect(pts))
    return (:none,nothing)
end

# Candidate intersection points of a triangle pair in exact rationals:
# edge-plane pierces inside the other facet, on-plane vertices inside the
# other facet, and proper edge-edge crossings. Coplanar arrangements return
# nothing — they stay on the contact/rejection paths.
function _brep_tri_candidates(T,U)
    s=ntuple(i->orient3(U[1],U[2],U[3],T[i]),3)
    all(==(0),s) && return nothing
    t=ntuple(j->orient3(T[1],T[2],T[3],U[j]),3)
    dropT=_brep_coplanar_drop_axis(T);dropU=_brep_coplanar_drop_axis(U)
    Tr=ntuple(i->_brep_rpoint(T[i]),3)
    Ur=ntuple(j->_brep_rpoint(U[j]),3)
    U2r=ntuple(j->_brep_rdrop(Ur[j],dropU),3)
    T2r=ntuple(i->_brep_rdrop(Tr[i],dropT),3)
    pts=NTuple{3,Rational{BigInt}}[]
    for i in 1:3
        j=i==3 ? 1 : i+1
        if s[i]*s[j]<0
            p=_brep_edge_plane_pierce(Ur,Tr[i],Tr[j])
            _brep_rpoint_in_tri_2d(_brep_rdrop(p,dropU),U2r) && push!(pts,p)
        end
        if t[i]*t[j]<0
            p=_brep_edge_plane_pierce(Tr,Ur[i],Ur[j])
            _brep_rpoint_in_tri_2d(_brep_rdrop(p,dropT),T2r) && push!(pts,p)
        end
        s[i]==0 && _brep_rpoint_in_tri_2d(_brep_rdrop(Tr[i],dropU),U2r) &&
            push!(pts,Tr[i])
        t[i]==0 && _brep_rpoint_in_tri_2d(_brep_rdrop(Ur[i],dropT),T2r) &&
            push!(pts,Ur[i])
    end
    for i in 1:3,j in 1:3
        p=_brep_segs_cross_point_3d(T[i],T[mod1(i+1,3)],U[j],U[mod1(j+1,3)])
        p===nothing || push!(pts,p)
    end
    unique!(pts)
    return pts
end

# Transversal tri-tri crossing: when the two planes differ, each triangle
# straddles the other's plane, and the candidate set spans a segment, return
# its endpoints as rational points. Every other arrangement (point touches,
# edge-in-plane, tangent contact) returns nothing — those stay on the
# contact/rejection paths.
function _brep_tri_crossing(T,U)
    s=ntuple(i->orient3(U[1],U[2],U[3],T[i]),3)
    all(==(0),s)||all(>(0),s)||all(<(0),s) && return nothing
    t=ntuple(j->orient3(T[1],T[2],T[3],U[j]),3)
    all(>(0),t)||all(<(0),t) && return nothing
    # A transversal segment needs each triangle to straddle the other's
    # plane — one-sided sign patterns are tangent/edge-in-plane contact.
    (any(>(0),s)&&any(<(0),s)&&any(>(0),t)&&any(<(0),t)) || return nothing
    pts=_brep_tri_candidates(T,U)
    (pts===nothing || length(pts)<2) && return nothing
    # Every candidate lies on the plane-intersection line; the segment's
    # endpoints are the extremes along its dominant axis.
    Tr=ntuple(i->_brep_rpoint(T[i]),3)
    Ur=ntuple(j->_brep_rpoint(U[j]),3)
    na=_brep_rcross(_brep_rsub(Tr[2],Tr[1]),_brep_rsub(Tr[3],Tr[1]))
    nb=_brep_rcross(_brep_rsub(Ur[2],Ur[1]),_brep_rsub(Ur[3],Ur[1]))
    dir=_brep_rcross(na,nb)
    axis=argmax(abs.((dir[1],dir[2],dir[3])))
    lo=argmin(p->p[axis],pts);hi=argmax(p->p[axis],pts)
    lo==hi && return nothing
    # Every candidate must lie on the segment — a non-collinear contact set
    # is not a single transversal crossing.
    all(p->_brep_rpoint_on_rseg(p,lo,hi),pts) || return nothing
    return lo,hi
end

# Non-axis-aligned probe directions; a ray that never runs parallel to a facet
# plane and pierces only facet interiors gives an unambiguous crossing count.
const _BREP_RAY_DIRS=((1,1,1),(1,2,3),(2,1,3),(1,3,2),(3,1,2),(2,3,1),(3,2,1),
                     (1,1,4),(1,4,1),(4,1,1),(2,3,5),(3,5,2),(5,2,3),(1,4,7),
                     (4,7,1),(7,1,4),(2,5,7),(5,7,2),(7,2,5),(3,4,9),(4,9,3),
                     (9,3,4),(5,6,7),(6,7,5),(7,5,6))
# Möller–Trumbore ray-facet test in exact rationals. :inside/:outside when every
# hit pierces a facet interior; :degenerate when a hit lands on a facet boundary
# (vertex/edge/tangent) so the direction is retried.
function _brep_ray_parity(pr,fr,normals,dr)
    all(n->_brep_rdot(n,dr)!=0,normals) || return :degenerate
    count=0
    for F in fr
        e1=_brep_rsub(F[2],F[1]);e2=_brep_rsub(F[3],F[1])
        pv=_brep_rcross(dr,e2)
        det=_brep_rdot(e1,pv)
        tv=_brep_rsub(pr,F[1])
        u=_brep_rdot(tv,pv)/det
        q=_brep_rcross(tv,e1)
        v=_brep_rdot(dr,q)/det
        t=_brep_rdot(e2,q)/det
        (u>=0&&v>=0&&u+v<=1&&t>0) || continue
        (u>0&&v>0&&u+v<1) || return :degenerate
        count+=1
    end
    return isodd(count) ? :inside : :outside
end
# Classify a point against a closed shell: :onboundary, :inside or :outside.
function _brep_point_vs_shell(p,facets,boxes,caller)
    for (F,box) in zip(facets,boxes)
        _brep_box_has(box,p) || continue
        orient3(F[1],F[2],F[3],p)==0 || continue
        drop=_brep_coplanar_drop_axis(F)
        _brep_point_in_tri_2d(_brep_drop(p,drop),
            ntuple(i->_brep_drop(F[i],drop),3)) && return :onboundary
    end
    pr=_brep_rpoint(p)
    fr=[ntuple(i->_brep_rpoint(F[i]),3) for F in facets]
    normals=[_brep_rcross(_brep_rsub(F[2],F[1]),_brep_rsub(F[3],F[1])) for F in fr]
    for dir in _BREP_RAY_DIRS
        dr=ntuple(d->Rational{BigInt}(dir[d]),3)
        result=_brep_ray_parity(pr,fr,normals,dr)
        result!==:degenerate && return result
    end
    throw(ArgumentError("$caller: cannot certify solid independence; " *
        "no nondegenerate containment probe exists"))
end

# Conforming contact: two facets are coincident when they share the exact same
# vertex positions. Their outward windings must oppose — equal winding means the
# interiors lie on the same side of the support, i.e. interpenetration.
function _brep_coincident_facet(T,U)
    Tr=ntuple(i->_brep_rpoint(T[i]),3)
    Ur=ntuple(j->_brep_rpoint(U[j]),3)
    Set(Tr)==Set(Ur) || return false
    na=_brep_rcross(_brep_rsub(Tr[2],Tr[1]),_brep_rsub(Tr[3],Tr[1]))
    nb=_brep_rcross(_brep_rsub(Ur[2],Ur[1]),_brep_rsub(Ur[3],Ur[1]))
    return _brep_rdot(na,nb)<0
end

# Certify one record pair. Returns (shared, enclosed): the coincident (ai,bi)
# facet pairs that weld into a shared interface, plus :a/:b when one shell is
# strictly enclosed by the other. Throws on any non-conforming contact or
# partial interpenetration.
function _brep_certify_disjoint(sa::Mesh,sb::Mesh,na,nb,caller)
    fa=_brep_shell_facets(sa);fb=_brep_shell_facets(sb)
    ba=[_brep_facet_box(F) for F in fa];bb=[_brep_facet_box(F) for F in fb]
    shared=Tuple{Int,Int}[]
    crossings=Tuple{Int,Int,NTuple{3,Rational{BigInt}},NTuple{3,Rational{BigInt}}}[]
    shared_set=Set{Tuple{Int,Int}}()
    coincident_a=falses(length(fa));coincident_b=falses(length(fb))
    contact_positions=Set{NTuple{3,Float64}}()
    extra=Tuple{Int,Int}[]
    for ai in eachindex(fa),bi in eachindex(fb)
        _brep_boxes_touch(ba[ai],bb[bi]) || continue
        _brep_coincident_facet(fa[ai],fb[bi]) || continue
        push!(shared,(ai,bi));push!(shared_set,(ai,bi))
        coincident_a[ai]=true;coincident_b[bi]=true
    end
    # A non-coincident contact is conforming only at coincident discrete
    # topology: a point common to both facets' vertices welds nothing but is
    # admitted, and a segment must be a shared edge — its endpoints vertices
    # of both facets. Anything else (interior pierces, coplanar area,
    # T-junctions, partial edge overlap) is non-conforming contact.
    for ai in eachindex(fa),bi in eachindex(fb)
        _brep_boxes_touch(ba[ai],bb[bi]) || continue
        (ai,bi) in shared_set && continue
        kind,pts=_brep_tris_contact(fa[ai],fb[bi])
        kind===:none && continue
        legit=false
        if kind===:point
            v=only(pts)
            legit=v in fa[ai] && v in fb[bi]
            legit && push!(contact_positions,v)
        elseif kind===:segment
            v1,v2=pts
            legit=v1 in fa[ai] && v1 in fb[bi] && v2 in fa[ai] && v2 in fb[bi]
            legit && (push!(contact_positions,v1);push!(contact_positions,v2))
        end
        legit && continue
        # A transversal crossing is a fragmentable intersection: record its
        # segment instead of rejecting. Every other non-conforming contact is
        # deferred — a corner or diagonal touch that lies on the collected
        # intersection polyline is conforming by construction; anything else
        # rejects below.
        seg=_brep_tri_crossing(fa[ai],fb[bi])
        seg===nothing || (push!(crossings,(ai,bi,seg[1],seg[2]));continue)
        push!(extra,(ai,bi))
    end
    for (ai,bi) in extra
        cps=_brep_tri_candidates(fa[ai],fb[bi])
        oncurve=cps!==nothing && all(p->any(
            s->_brep_rpoint_on_rseg(p,s[3],s[4]),crossings),cps)
        oncurve || throw(ArgumentError(
            "$caller: cannot certify independent solids $na and $nb; " *
            "their boundaries touch or intersect"))
    end
    (isempty(fa)||isempty(fb)) && return shared,nothing,crossings
    if !isempty(crossings)
        # Crossings mixed with coincident facets on the same pair need a
        # shared-face classification that this slice does not support.
        isempty(shared) || throw(ArgumentError(
            "$caller: cannot certify independent solids $na and $nb; " *
            "their boundaries intersect across coincident facets"))
        return shared,nothing,crossings
    end
    # Every referenced vertex off the shared patch is probed: after the pair
    # sweep no free vertex can lie on the other boundary (non-coincident
    # contact already threw), so each is strictly inside or outside. A shell
    # whose free vertices are all inside encloses nothing — it IS enclosed,
    # admitted as a nested Region; a mix of inside and outside free vertices
    # is interpenetration (a facet hinged on the rim can dip into the other
    # solid). A shell with no free vertex is a duplicate or an embedded
    # subshell, never independent.
    enclosed=nothing
    for (surface,facets,boxes,coincident,direction) in
            ((sa,fb,bb,coincident_a,:a),(sb,fa,ba,coincident_b,:b))
        used=falses(nnodes(surface));referenced=falses(nnodes(surface))
        for column in 1:ntris(surface)
            on_iface=coincident[column]
            for r in 1:3
                row=surface.tris[r,column];referenced[row]=true
                (on_iface || node(surface,row) in contact_positions) &&
                    (used[row]=true)
            end
        end
        probes=0;inside=0
        for row in 1:nnodes(surface)
            (referenced[row] && !used[row]) || continue
            probes+=1
            state=_brep_point_vs_shell(node(surface,row),facets,boxes,caller)
            state===:outside && continue
            state===:inside || throw(ArgumentError(
                "$caller: cannot certify independent solids $na and $nb; " *
                "their boundaries touch or intersect"))
            inside+=1
        end
        # A shell can have no free vertex yet carry non-coincident facets —
        # a region wrapped in shared interfaces, whose facet interiors still
        # decide enclosure. A non-coincident facet's interior cannot touch
        # the other boundary (the pair sweep already classified that), so
        # strict interior probes are decisive; only when every facet is
        # coincident is the shell a duplicate or an embedded subshell.
        if probes==0
            for column in 1:ntris(surface)
                coincident[column] && continue
                probes+=1
                probe=(node(surface,surface.tris[1,column]).+
                    node(surface,surface.tris[2,column]).+
                    node(surface,surface.tris[3,column]))./3
                state=_brep_point_vs_shell(probe,facets,boxes,caller)
                state===:outside && continue
                state===:inside || throw(ArgumentError(
                    "$caller: cannot certify independent solids $na and $nb; " *
                    "their boundaries touch or intersect"))
                inside+=1
            end
            probes==0 && throw(ArgumentError("$caller: cannot certify independent solids $na and $nb; " *
                "their boundaries coincide everywhere"))
        end
        # A shell wholly inside the other encloses a nested Region; a mix of
        # inside and outside free vertices is interpenetration.
        0<inside<probes && throw(ArgumentError(
            "$caller: cannot certify independent solids $na and $nb; " *
            "their boundaries interpenetrate"))
        if inside>0
            enclosed===nothing || throw(ArgumentError(
                "$caller: cannot certify independent solids $na and $nb; " *
                "their boundaries interpenetrate"))
            enclosed=direction
        end
    end
    (enclosed===nothing || (isempty(shared) && isempty(contact_positions))) ||
        throw(ArgumentError("$caller: cannot certify independent solids $na and $nb; " *
            "an enclosed solid cannot share boundary contact"))
    return shared,enclosed,crossings
end

# Does segment (p,q) pass through the strict interior of triangle T2?
# Boundary-only touches return false.
function _brep_seg_in_tri_interior_2d(p,q,T2)
    # A segment lying on an edge line stays on the boundary — it never
    # enters the strict interior, and Float64 materialization makes its
    # midpoint test unreliable (a boundary midpoint lands one ulp inside).
    for i in 1:3
        a,b=T2[i],T2[mod1(i+1,3)]
        orient2(a,b,p)==0 && orient2(a,b,q)==0 && return false
    end
    strict_in(x)=begin
        s1,s2,s3=orient2(T2[1],T2[2],x),orient2(T2[2],T2[3],x),orient2(T2[3],T2[1],x)
        (s1>0&&s2>0&&s3>0)||(s1<0&&s2<0&&s3<0)
    end
    strict_in(p) && return true
    strict_in(q) && return true
    # A chord whose endpoints sit on the cell boundary dips inside exactly
    # when its midpoint does — the triangle is convex.
    strict_in((p .+ q) ./ 2) && return true
    for i in 1:3
        a,b=T2[i],T2[mod1(i+1,3)]
        orient2(a,b,p)*orient2(a,b,q)<0 && orient2(p,q,a)*orient2(p,q,b)<0 &&
            return true
    end
    return false
end

# Is triangle P2 covered by the closed union of coplanar patch triangles G2?
# Every vertex must lie in some patch triangle and no patch rim edge may cut
# the cell's interior — the rim is the support-edge set of single incidence.
function _brep_cell_in_union(P2,G2)
    all(k->any(T2->_brep_point_in_tri_2d(P2[k],T2),G2),1:3) || return false
    counts=Dict{NTuple{2,NTuple{2,Float64}},Int}()
    for T2 in G2,e in 1:3
        edge=minmax(T2[e],T2[mod1(e+1,3)])
        counts[edge]=get(counts,edge,0)+1
    end
    for (edge,count) in counts
        count==1 || continue
        _brep_seg_in_tri_interior_2d(edge[1],edge[2],P2) && return false
    end
    return true
end

# Positive-area 2-D overlap between two projected triangles: strict edge
# crossings or a strict vertex containment. Rim/point touches return false.
function _brep_tris_overlap_area_2d(T2,U2)
    for i in 1:3,j in 1:3
        a,b,c,d=T2[i],T2[mod1(i+1,3)],U2[j],U2[mod1(j+1,3)]
        orient2(a,b,c)*orient2(a,b,d)<0 && orient2(c,d,a)*orient2(c,d,b)<0 &&
            return true
    end
    strict_in(p,T)=begin
        s1,s2,s3=orient2(T[1],T[2],p),orient2(T[2],T[3],p),orient2(T[3],T[1],p)
        (s1>0&&s2>0&&s3>0)||(s1<0&&s2<0&&s3<0)
    end
    orient2(U2[1],U2[2],U2[3])!=0 &&
        any(p->strict_in(p,U2),T2) && return true
    orient2(T2[1],T2[2],T2[3])!=0 &&
        any(p->strict_in(p,T2),U2) && return true
    return false
end

# A rational intersection vertex materialized in Float64. The same geometric
# point computed on either side of a crossing pair yields the same rational,
# so this is also the global dedup key for shared cut vertices.
@inline _brep_ftuple(p)=(Float64(p[1]),Float64(p[2]),Float64(p[3]))

# Cut one closed boundary along its transversal intersection segments. Facets
# the curve misses pass through verbatim; a cut facet is retriangulated in its
# own plane by the 2-D constrained Delaunay kernel — a triangle is its own
# convex hull, so no boundary constraints are needed, only the interior cut
# segments. Emitted pieces keep the source facet's winding, and segment
# endpoints lying on a facet edge are propagated so both incident facets
# split the shared edge identically.
function _brep_cut_surface(facets,segs,caller)
    segends=[x for list in values(segs) for seg in list for x in seg]
    pieces=NTuple{3,NTuple{3,Float64}}[]
    for f in eachindex(facets)
        T=facets[f]
        list=get(segs,f,nothing)
        # The curve can graze an edge without entering this facet: a segment
        # endpoint on a facet edge must split the edge even when the facet
        # carries no segment, so the shared edge is cut identically on both
        # incident facets.
        edgepts=NTuple{3,Rational{BigInt}}[]
        for (a,b) in ((1,2),(2,3),(3,1))
            ra=_brep_rpoint(T[a]);rb=_brep_rpoint(T[b])
            for p in segends
                _brep_rpoint_on_rseg(p,ra,rb) || continue
                p in edgepts || push!(edgepts,p)
            end
        end
        if list===nothing && isempty(edgepts)
            push!(pieces,T);continue
        end
        drop=_brep_coplanar_drop_axis(T)
        p2=function(p)
            d=_brep_drop(p,drop)
            (d[1]==0 ? 0.0 : d[1],d[2]==0 ? 0.0 : d[2])
        end
        pts=NTuple{3,Float64}[T[1],T[2],T[3]]
        # Rational representatives for the same points: Float64 rounding can
        # break exact collinearity, so degeneracy is decided rationally.
        rkey=Dict{NTuple{2,Float64},NTuple{3,Rational{BigInt}}}(
            p2(T[k])=>_brep_rpoint(T[k]) for k in 1:3)
        for p in something(list,Tuple{NTuple{3,Rational{BigInt}},
                NTuple{3,Rational{BigInt}}}[])
            pf,qf=_brep_ftuple(p[1]),_brep_ftuple(p[2])
            pf in pts || push!(pts,pf)
            qf in pts || push!(pts,qf)
            rkey[p2(pf)]=p[1];rkey[p2(qf)]=p[2]
        end
        for p in edgepts
            pf=_brep_ftuple(p);pf in pts || push!(pts,pf)
            rkey[p2(pf)]=p
        end
        xs=Float64[p2(p)[1] for p in pts]
        ys=Float64[p2(p)[2] for p in pts]
        T2=Mesh2D.delaunay2d(xs,ys)
        locate=function(p)
            d=p2(p)
            findfirst(v->T2.x[v]==d[1]&&T2.y[v]==d[2],Int32(1):T2.nreal)
        end
        list===nothing || for (p,q) in list
            vi=locate(_brep_ftuple(p));vj=locate(_brep_ftuple(q))
            (vi===nothing || vj===nothing) && throw(ErrorException(
                "$caller: internal error — cut endpoint missing from the " *
                "facet triangulation"))
            vi!=vj && Mesh2D.insert_segment!(T2,vi,vj;internal=true)
        end
        m=Mesh2D.to_mesh(T2)
        keymap=Dict{NTuple{2,Float64},NTuple{3,Float64}}()
        for p in pts
            keymap[p2(p)]=p
        end
        flip=orient2(_brep_drop(T[1],drop),_brep_drop(T[2],drop),
            _brep_drop(T[3],drop))<0
        for t in 1:ntris(m)
            tri=ntuple(r->keymap[(m.coords[1,m.tris[r,t]],
                m.coords[2,m.tris[r,t]])],3)
            # Cut endpoints collinear on a facet edge can make the kernel
            # emit a sliver that is degenerate in exact arithmetic; it
            # carries no surface and its edges are subsumed by the fine
            # split, so it is dropped.
            r1,r2,r3=ntuple(r->rkey[(m.coords[1,m.tris[r,t]],
                m.coords[2,m.tris[r,t]])],3)
            _brep_rcross(_brep_rsub(r2,r1),_brep_rsub(r3,r1))==(0,0,0) &&
                continue
            push!(pieces,flip ? (tri[1],tri[3],tri[2]) : tri)
        end
    end
    return pieces
end

# Transversal crossings fragment the crossed boundaries: every affected facet
# is cut along its exact rational intersection segments, each resulting piece
# is classified against every other solid, and the pieces regroup into region
# shells — a piece on surface i inside exactly the solids in C borders region
# C∪{i} on its outward side and region C on its inward side. The region shells
# share coincident interface facets, so they re-enter the ordinary assembly
# pipeline where the certificate welds them like any conforming contact.
function _brep_assembly_fragment(records,fill_solids,crossingmap,caller)
    n=length(records)
    facets=[_brep_shell_facets(r.surface) for r in records]
    boxes=[[_brep_facet_box(F) for F in facets[i]] for i in 1:n]
    segs=[Dict{Int,Vector{Tuple{NTuple{3,Rational{BigInt}},
        NTuple{3,Rational{BigInt}}}}}() for _ in 1:n]
    for ((i,j),list) in crossingmap,(ai,bi,p,q) in list
        push!(get!(segs[i],ai,Tuple{NTuple{3,Rational{BigInt}},
            NTuple{3,Rational{BigInt}}}[]),(p,q))
        push!(get!(segs[j],bi,Tuple{NTuple{3,Rational{BigInt}},
            NTuple{3,Rational{BigInt}}}[]),(p,q))
    end
    regions=Dict{Vector{Int},Vector{Tuple{NTuple{3,NTuple{3,Float64}},Bool}}}()
    probes=((1,1,1),(3,1,1),(1,3,1),(1,1,3),(2,2,1),(2,1,2),(1,2,2))
    for i in 1:n
        pieces=_brep_cut_surface(facets[i],segs[i],caller)
        for P in pieces
            # The centroid can land on another surface's facet plane; probe
            # several strict interior points and take the first decisive one.
            # Every probe on the boundary means a coincident piece — the
            # bounded slice rejects it.
            interior=ntuple(k->(probes[k][1].*P[1].+probes[k][2].*P[2].+
                probes[k][3].*P[3])./sum(probes[k]),length(probes))
            C=Int[]
            for j in 1:n
                j==i && continue
                state=:onboundary
                for p in interior
                    st=_brep_point_vs_shell(p,facets[j],boxes[j],caller)
                    st===:onboundary && continue
                    state===:onboundary && (state=st;continue)
                    st===state || throw(ArgumentError(
                        "$caller: inconsistent fragment classification — a " *
                        "cut piece straddles another surface"))
                end
                state===:inside && push!(C,j)
                # :onboundary can only remain if every probe landed on the
                # other surface — a genuinely coincident piece.
                state===:onboundary && throw(ArgumentError(
                    "$caller: cannot assemble solids; a fragment boundary " *
                    "piece lies on another solid's surface"))
            end
            push!(get!(regions,sort!([C;i]),
                Tuple{NTuple{3,NTuple{3,Float64}},Bool}[]),(P,false))
            isempty(C) || push!(get!(regions,C,
                Tuple{NTuple{3,NTuple{3,Float64}},Bool}[]),(P,true))
        end
    end
    pseudo=NamedTuple{(:name,:surface),Tuple{String,Mesh}}[]
    for S in sort!(collect(keys(regions));by=s->(length(s),s))
        tris=regions[S]
        pos2row=Dict{NTuple{3,Float64},Int32}()
        conn=Matrix{Int32}(undef,3,length(tris))
        for (c,(P,rev)) in enumerate(tris)
            for (k,v) in enumerate(rev ? (3,2,1) : (1,2,3))
                pos=P[v]
                row=get!(pos2row,pos) do
                    Int32(length(pos2row)+1)
                end
                conn[k,c]=row
            end
        end
        cm=Matrix{Float64}(undef,3,length(pos2row))
        for (pos,row) in pos2row
            cm[1,row]=pos[1];cm[2,row]=pos[2];cm[3,row]=pos[3]
        end
        name=join((records[k].name for k in S),"∩")
        push!(pseudo,(name=name,surface=Mesh(cm;tris=conn)))
    end
    # Region fragments share coincident interface facets; the joint fill tags
    # every region in one conforming tetrahedralization so the welded
    # interfaces are identical by construction.
    return _brep_assembly(pseudo,fill_solids,caller;joint_fill=fill_solids)
end

# Mesh import retains the source solid partition. These are discrete closed
# boundary entities, not reconstructed analytic CAD Face/Curve parameter maps.
function _brep_assembly(records,fill_solids,caller;joint_fill::Bool=false)
    length(records)>1 || throw(ArgumentError("$caller: assembly requires multiple source solids"))
    length(records)<=typemax(Int32) || throw(ArgumentError("$caller: too many source solids"))
    bounds=[bounding_box(record.surface) for record in records]
    weld_pairs=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    nested=Set{Tuple{Int,Int}}()
    fragment=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int,NTuple{3,Rational{BigInt}},NTuple{3,Rational{BigInt}}}}}()
    for i in 1:length(records)-1,j in i+1:length(records)
        alo,ahi=bounds[i];blo,bhi=bounds[j]
        any(axis->ahi[axis]<blo[axis] || bhi[axis]<alo[axis],1:3) && continue
        shared,enclosed,crossings=_brep_certify_disjoint(records[i].surface,records[j].surface,
            records[i].name,records[j].name,caller)
        isempty(shared) || (weld_pairs[(i,j)]=shared)
        enclosed===:a && push!(nested,(i,j))
        enclosed===:b && push!(nested,(j,i))
        isempty(crossings) || (fragment[(i,j)]=crossings)
    end
    if !isempty(fragment)
        # A crossing pair mixed with a coincident weld or enclosure on ANY
        # pair needs classifications this bounded slice does not support.
        (isempty(weld_pairs) && isempty(nested)) || throw(ArgumentError(
            "$caller: cannot assemble solids; transversal boundary " *
            "intersections combined with coincident or enclosed solids are " *
            "not supported"))
        return _brep_assembly_fragment(records,fill_solids,fragment,caller)
    end
    # Region-tag ordering counts every container of each solid, so the
    # innermost shells tag first. Only direct containment welds an interface:
    # a transitively enclosed solid has no emitted interface with the
    # outermost container (the middle Region's boundary separates them).
    depth=zeros(Int,length(records))
    for (i,j) in nested;depth[i]+=1;end
    for (i,j) in collect(nested)
        any(k->k!=i && k!=j && (i,k) in nested && (k,j) in nested,
            1:length(records)) && delete!(nested,(i,j))
    end
    # Certify every source boundary before any fill. Independent fills leave
    # an enclosed Region overlapped by its container's fill, so nested
    # assemblies take one joint conforming fill — surfaces ordered innermost
    # first so region tags pick the innermost container.
    parts=Mesh[];shells=Matrix{Int32}[]
    if fill_solids && (!isempty(nested) || joint_fill)
        order=sortperm(depth;rev=true)
        joint=tetrahedralize_conforming_exact(
            [records[i].surface for i in order])
        for (r,record) in enumerate(records)
            validate(record.surface).ok || throw(ArgumentError(
                "$caller: invalid source solid $(record.name)"))
            _brep_assembly_shell_certificate(record.surface.coords,
                record.surface.tris,caller)
            pos=findfirst(==(r),order)
            tetcols=[c for c in axes(joint.tets,2) if joint.tet_tag[c]==Int32(pos)]
            rows=sort!(unique!(vec(joint.tets[:,tetcols])))
            localmap=Dict{Int32,Int32}(row=>Int32(i) for (i,row) in enumerate(rows))
            tets=Matrix{Int32}(undef,4,length(tetcols))
            for (c,col) in enumerate(tetcols),k in 1:4
                tets[k,c]=localmap[joint.tets[k,col]]
            end
            part=Mesh(joint.coords[:,rows];tets=tets)
            shell=_brep_tet_shell(part,caller)
            _brep_assembly_shell_certificate(part.coords,shell,caller)
            push!(parts,part);push!(shells,shell)
        end
    else
        for record in records
            validate(record.surface).ok || throw(ArgumentError("$caller: invalid source solid $(record.name)"))
            _brep_assembly_shell_certificate(record.surface.coords,record.surface.tris,caller)
            part=fill_solids ? _filled(record.surface,caller) : record.surface
            shell=fill_solids ? _brep_tet_shell(part,caller) : copy(part.tris)
            _brep_assembly_shell_certificate(part.coords,shell,caller)
            push!(parts,part);push!(shells,shell)
        end
    end
    # Nested interfaces weld only under the joint fill — with fill=false the
    # container's source has no cells on the inner shell, so nested Regions
    # keep their own coincident boundary emission like other coincident
    # contacts.
    effective_nested=fill_solids ? nested : empty(nested)
    (isempty(weld_pairs) && isempty(effective_nested)) ||
        return _brep_assembly_welded(records,parts,shells,
            bounds,fill_solids,weld_pairs,effective_nested,caller)
    node_count=0;cell_count=0
    for (part,shell) in zip(parts,shells)
        node_count=_checked_brep_add(node_count,nnodes(part),0,caller,"assembly node count")
        cell_count=_checked_brep_add(cell_count,size(shell,2),ntets(part),caller,"assembly cell count")
    end
    node_count<=typemax(Int32) && cell_count<=typemax(Int32) ||
        throw(ArgumentError("$caller: assembly exceeds Int32 topology limits"))
    coords=Matrix{Float64}(undef,3,node_count)
    entities=Dict{Tuple{Int,Int},MixedEntity}()
    names=Dict{Tuple{Int,Int},String}()
    owners=Tuple{Int,Int32}[];blocks=ElementBlock[]
    block_owners=Vector{Int32}[];labels=Vector{UInt64}[]
    node_offset=0;cell_offset=UInt64(0)
    for (region,(record,part,shell)) in enumerate(zip(records,parts,shells))
        tag=Int32(region);lo,hi=bounds[region];bbox=(lo...,hi...)
        entities[(2,region)]=MixedEntity(2,tag,bbox;physical_tags=[tag])
        entities[(3,region)]=MixedEntity(3,tag,bbox;physical_tags=[tag],boundaries=[tag])
        names[(2,region)]=record.name*" boundary"
        names[(3,region)]=record.name
        count=nnodes(part);coords[:,node_offset+1:node_offset+count]=part.coords
        boundary=falses(count)
        for row in shell;boundary[row]=true;end
        append!(owners,[(boundary[row] ? 2 : 3,tag) for row in 1:count])
        for (kind,cells) in ((2,shell),(4,part.tets))
            isempty(cells) && continue
            push!(blocks,ElementBlock(kind,cells .+ Int32(node_offset),fill(tag,size(cells,2))))
            push!(block_owners,fill(tag,size(cells,2)))
            push!(labels,UInt64[cell_offset+UInt64(i) for i in 1:size(cells,2)])
            cell_offset+=UInt64(size(cells,2))
        end
        node_offset+=count
    end
    data=MixedEntityData(entities;node_entities=owners,
        node_parametric=fill(nothing,node_count),external_node_tags=UInt64.(1:node_count),
        block_entities=block_owners,external_element_tags=labels)
    result=MixedMesh(coords,blocks;physical_names=names,entity_data=data,
                     elementary_entities=block_owners)
    report=validate(result)
    report.ok || throw(ErrorException("$caller: invalid classified assembly product — " *
        join(report.messages,"; ")))
    return result
end

# Conforming-contact emission: coincident facet pairs weld into one shared
# Surface entity bound by both Regions. Nodes on shared supports deduplicate
# by exact position; every emitted boundary cell inside a shared footprint
# must appear identically on both records' emitted shells.
function _brep_assembly_welded(records,parts,shells,bounds,fill_solids,weld_pairs,nested,caller)
    nrec=length(records)
    # Canonical shared-support key per coincident source facet: the sorted
    # vertex-position triple — identical on both records.
    support_key=Dict{Tuple{Int,Int},NTuple{3,NTuple{3,Float64}}}()
    for (key,list) in weld_pairs,(ai,bi) in list
        i,j=key
        for (r,f) in ((i,ai),(j,bi))
            haskey(support_key,(r,f)) && throw(ArgumentError(
                "$caller: cannot certify conforming contact; " *
                "a boundary facet is shared by three solids"))
            support_key[(r,f)]=Tuple(sort!(collect(ntuple(
                c->node(records[r].surface,records[r].surface.tris[c,f]),3))))
        end
    end
    # One shared Surface entity per welded record pair — coincident facet
    # pairs and nested (enclosed Region) interfaces alike.
    pair_tags=Dict{Tuple{Int,Int},Int32}()
    for key in sort!(collect(union!(collect(keys(weld_pairs)),collect(nested))))
        pair_tags[key]=Int32(nrec+length(pair_tags)+1)
    end
    # Each pair's supports group by coplanarity into plane patches; a fill
    # cell is shared when covered by one patch's union, so fills may
    # retriangulate freely across internal patch edges.
    patches=Dict{Tuple{Int,Int},Vector{Vector{NTuple{3,NTuple{3,Float64}}}}}()
    for (key,list) in weld_pairs,(ai,bi) in list
        support=support_key[(key[1],ai)]
        groups=get!(patches,key,Vector{NTuple{3,NTuple{3,Float64}}}[])
        placed=false
        for group in groups
            S=group[1]
            if all(p->orient3(S[1],S[2],S[3],p)==0,support)
                push!(group,support);placed=true;break
            end
        end
        placed || push!(groups,[support])
    end
    # A nested pair's shared interface is the enclosed Region's whole source
    # boundary; emitted cells on either side are matched by their vertices —
    # the joint fill never inserts new nodes, so every interface vertex is a
    # source vertex of the enclosed solid. Cells on a deeper nested interface
    # (a void inside the enclosed Region) carry different vertices and stay
    # with their own pair.
    nested_verts=Dict{Tuple{Int,Int},Set{NTuple{3,Float64}}}()
    for (i,j) in nested
        set=Set{NTuple{3,Float64}}()
        for column in axes(records[i].surface.tris,2)
            for k in 1:3
                push!(set,node(records[i].surface,records[i].surface.tris[k,column]))
            end
        end
        nested_verts[(i,j)]=set
    end
    # Classify every emitted boundary cell against each pair's interfaces:
    # full closed coverage by a patch union makes it an interface cell;
    # positive-area overlap without coverage means the fill crosses the
    # interface rim. An inner Region's whole boundary is shared.
    emitted_shared=[Dict{Tuple{Int,Int},Vector{NTuple{3,Int}}}()
                    for _ in records]
    for (r,(part,shell)) in enumerate(zip(parts,shells))
        relevant=sort!(Tuple{Int,Int}[key for key in union(keys(patches),
            keys(nested_verts)) if r in key])
        isempty(relevant) && continue
        for column in axes(shell,2)
            cell=(Int(shell[1,column]),Int(shell[2,column]),Int(shell[3,column]))
            P=ntuple(k->node(part,cell[k]),3)
            # Coincident-facet patches classify first — their coverage is
            # exact. A cell on an enclosed Region's welded face carries the
            # inner vertices but borders the weld partner, not the container,
            # so nested pairs claim only the leftovers.
            claimed=false
            for key in relevant
                haskey(nested_verts,key) && continue
                for group in patches[key]
                    S=group[1]
                    all(k->orient3(S[1],S[2],S[3],P[k])==0,1:3) || continue
                    drop=_brep_coplanar_drop_axis(S)
                    P2=ntuple(k->_brep_drop(P[k],drop),3)
                    G2=[ntuple(k->_brep_drop(s[k],drop),3) for s in group]
                    if _brep_cell_in_union(P2,G2)
                        push!(get!(emitted_shared[r],key,NTuple{3,Int}[]),cell)
                        claimed=true;break
                    end
                    any(T2->_brep_tris_overlap_area_2d(P2,T2),G2) && throw(
                        ArgumentError("$caller: conforming contact requires " *
                            "fills to respect the shared interface boundary"))
                end
                claimed && break
            end
            claimed && continue
            for key in relevant
                haskey(nested_verts,key) || continue
                all(k->P[k] in nested_verts[key],1:3) &&
                    push!(get!(emitted_shared[r],key,NTuple{3,Int}[]),cell)
            end
        end
    end
    # Both records must emit the same cells on every shared interface —
    # otherwise their fills triangulate it differently.
    for key in keys(pair_tags)
        i,j=key
        ca=get(emitted_shared[i],key,NTuple{3,Int}[])
        cb=get(emitted_shared[j],key,NTuple{3,Int}[])
        (isempty(ca)||isempty(cb)) && throw(ArgumentError(
            "$caller: fill produced no cell on the shared interface"))
        keys_a=Set(Tuple(sort!(collect(ntuple(k->node(parts[i],cell[k]),3))))
                   for cell in ca)
        keys_b=Set(Tuple(sort!(collect(ntuple(k->node(parts[j],cell[k]),3))))
                   for cell in cb)
        keys_a==keys_b || throw(ArgumentError(
            "$caller: conforming contact requires identical triangulation " *
            "of the shared interface"))
    end
    # Global node rows deduplicate by exact position; a position claimed on a
    # shared cell or support rim belongs to that pair's shared entity.
    pos_owner=Dict{NTuple{3,Float64},Tuple{Int,Int32}}()
    for key in keys(pair_tags)
        i,j=key;tag=pair_tags[key]
        for (ai,bi) in get(weld_pairs,key,Tuple{Int,Int}[])
            support=support_key[(i,ai)]
            for p in support
                pos_owner[p]=min(get(pos_owner,p,(2,tag)),(2,tag))
            end
        end
        for r in (i,j),cell in emitted_shared[r][key],k in 1:3
            p=node(parts[r],cell[k])
            pos_owner[p]=min(get(pos_owner,p,(2,tag)),(2,tag))
        end
    end
    pos2row=Dict{NTuple{3,Float64},Int}()
    global_nodes=NTuple{3,Float64}[]
    owners=Tuple{Int,Int32}[]
    remap=[Vector{Int}(undef,nnodes(part)) for part in parts]
    for (r,(part,shell)) in enumerate(zip(parts,shells))
        tag=Int32(r);boundary=falses(nnodes(part))
        for row in shell;boundary[row]=true;end
        for row in 1:nnodes(part)
            p=node(part,row)
            grow=get!(pos2row,p) do
                owner=get(pos_owner,p,boundary[row] ? (2,tag) : (3,tag))
                push!(global_nodes,p);push!(owners,owner)
                length(global_nodes)
            end
            if boundary[row] && !haskey(pos_owner,p) && owners[grow][1]==3
                owners[grow]=(2,tag)
            end
            remap[r][row]=grow
        end
    end
    node_count=length(global_nodes)
    node_count<=typemax(Int32) ||
        throw(ArgumentError("$caller: assembly exceeds Int32 topology limits"))
    coords=Matrix{Float64}(undef,3,node_count)
    for (row,p) in enumerate(global_nodes);coords[:,row]=collect(p);end
    entities=Dict{Tuple{Int,Int},MixedEntity}()
    names=Dict{Tuple{Int,Int},String}()
    blocks=ElementBlock[];block_owners=Vector{Int32}[];labels=Vector{UInt64}[]
    cell_offset=UInt64(0)
    shared_cells=[Set{NTuple{3,Int}}(Iterators.flatten(values(cell_map)))
                  for cell_map in emitted_shared]
    for (region,(record,part,shell)) in enumerate(zip(records,parts,shells))
        tag=Int32(region);lo,hi=bounds[region];bbox=(lo...,hi...)
        shared_tags=Int32[pair_tags[key] for key in sort!(collect(keys(pair_tags)))
                          if region in key]
        entities[(2,region)]=MixedEntity(2,tag,bbox;physical_tags=[tag])
        entities[(3,region)]=MixedEntity(3,tag,bbox;physical_tags=[tag],
            boundaries=Int32[tag;shared_tags])
        names[(2,region)]=record.name*" boundary"
        names[(3,region)]=record.name
        own=Vector{Int32}[]
        for column in axes(shell,2)
            cell=(Int(shell[1,column]),Int(shell[2,column]),Int(shell[3,column]))
            cell in shared_cells[region] && continue
            push!(own,Int32[remap[region][cell[1]],remap[region][cell[2]],
                            remap[region][cell[3]]])
        end
        nt=ntets(part)
        for (kind,cols,raw) in ((2,length(own),own),(4,nt,nothing))
            cols>0 || continue
            cells=Matrix{Int32}(undef,kind==2 ? 3 : 4,cols)
            if kind==2
                for (c,cell) in enumerate(own);cells[:,c]=cell;end
            else
                for c in 1:nt,k in 1:4
                    cells[k,c]=Int32(remap[region][Int(part.tets[k,c])])
                end
            end
            push!(blocks,ElementBlock(kind,cells,fill(tag,cols)))
            push!(block_owners,fill(tag,cols))
            push!(labels,UInt64[cell_offset+UInt64(c) for c in 1:cols])
            cell_offset+=UInt64(cols)
        end
    end
    # Shared interface cells emit once, from the first record's side — the
    # enclosed Region's own boundary for a nested pair.
    for key in sort!(collect(keys(pair_tags)))
        i,j=key;tag=pair_tags[key]
        emitted=Vector{Int32}[]
        for cell in emitted_shared[i][key]
            push!(emitted,Int32[remap[i][cell[1]],remap[i][cell[2]],
                                remap[i][cell[3]]])
        end
        isempty(emitted) && continue
        cells=Matrix{Int32}(undef,3,length(emitted))
        for (c,cell) in enumerate(emitted);cells[:,c]=cell;end
        push!(blocks,ElementBlock(2,cells,fill(tag,length(emitted))))
        push!(block_owners,fill(tag,length(emitted)))
        push!(labels,UInt64[cell_offset+UInt64(c) for c in 1:length(emitted)])
        cell_offset+=UInt64(length(emitted))
        lo=fill(Inf,3);hi=fill(-Inf,3)
        for cell in emitted
            p=global_nodes[Int(cell[1])]
            lo=min.(lo,p);hi=max.(hi,p)
            for k in 2:3
                p=global_nodes[Int(cell[k])];lo=min.(lo,p);hi=max.(hi,p)
            end
        end
        entities[(2,Int(tag))]=MixedEntity(2,tag,(lo...,hi...);physical_tags=[tag])
        names[(2,Int(tag))]=records[i].name*"/"*records[j].name*" shared boundary"
    end
    cell_offset<=typemax(Int32) ||
        throw(ArgumentError("$caller: assembly exceeds Int32 topology limits"))
    data=MixedEntityData(entities;node_entities=owners,
        node_parametric=fill(nothing,node_count),external_node_tags=UInt64.(1:node_count),
        block_entities=block_owners,external_element_tags=labels)
    result=MixedMesh(coords,blocks;physical_names=names,entity_data=data,
                     elementary_entities=block_owners)
    report=validate(result)
    report.ok || throw(ErrorException("$caller: invalid classified assembly product — " *
        join(report.messages,"; ")))
    return result
end

function _brep_assembly_shell_certificate(coords,shell,caller)
    isempty(shell) && throw(ArgumentError("$caller: source solid has no closed shell"))
    incidences=Dict{Tuple{Int32,Int32},Tuple{Int,Int}}()
    for face in eachcol(shell)
        a,b,c=face
        for (u,v) in ((a,b),(b,c),(c,a))
            key=minmax(u,v);count,orientation=get(incidences,key,(0,0))
            incidences[key]=(count+1,orientation+(u<v ? 1 : -1))
        end
    end
    all(value->value==(2,0),values(incidences)) || throw(ArgumentError(
        "$caller: source solid boundary is not a consistently oriented closed shell"))
    return nothing
end

# Literal positive-tet face traversal supplies outward boundary cells. A Dict
# counts incidence only; the first-encounter vector controls emitted order.
function _brep_tet_shell(mesh,caller)
    counts=Dict{NTuple{3,Int32},Int}();faces=NTuple{3,Int32}[];keys=NTuple{3,Int32}[]
    for cell in eachcol(mesh.tets)
        a,b,c,d=cell
        for face in ((b,c,d),(a,d,c),(a,b,d),(a,c,b))
            key=Tuple(sort(collect(face)))
            if !haskey(counts,key);push!(keys,key);push!(faces,face);end
            counts[key]=get(counts,key,0)+1
            counts[key]<=2 || throw(ArgumentError("$caller: nonmanifold filled solid boundary"))
        end
    end
    selected=[faces[i] for i in eachindex(faces) if counts[keys[i]]==1]
    isempty(selected) && throw(ArgumentError("$caller: filled solid has no exterior shell"))
    shell=Matrix{Int32}(undef,3,length(selected))
    for (column,face) in enumerate(selected);shell[:,column]=collect(face);end
    return shell
end

function _iges_assembly_transform(surface,types,params,bad,transforms,de,caller)
    xf=_iges_de_transform(types,params,bad,transforms,de,caller)
    xf===nothing && return surface
    coordinates=Matrix{Float64}(undef,size(surface.coords))
    for row in axes(coordinates,2)
        coordinates[:,row]=collect(_iges_xf_apply(xf,node(surface,row),true))
    end
    R,_=xf
    determinant=R[1][1]*(R[2][2]*R[3][3]-R[2][3]*R[3][2])-
        R[1][2]*(R[2][1]*R[3][3]-R[2][3]*R[3][1])+R[1][3]*(R[2][1]*R[3][2]-R[2][2]*R[3][1])
    return Mesh(coordinates;tris=determinant<0 ? surface.tris[[1,3,2],:] : surface.tris)
end

function _step_solid_roots(entities)
    breps=[ent for ent in values(entities) if ent.kind in _STEP_BREP_SOLID_KINDS]
    descendants=Set{Int}();queue=Any[ent.args for ent in breps]
    while !isempty(queue)
        value=pop!(queue)
        if value isa Tuple && length(value)==2 && value[1]===:ref
            id=value[2];id in descendants && continue
            push!(descendants,id);ent=get(entities,id,nothing)
            ent===nothing || push!(queue,ent.args)
        elseif value isa AbstractVector
            append!(queue,value)
        elseif value isa StepEntity
            push!(queue,value.args)
        end
    end
    primitives=[ent for ent in values(entities) if !(ent.id in descendants) &&
        ent.kind in ("SPHERE","RIGHT_CIRCULAR_CYLINDER","RIGHT_CIRCULAR_CONE","TOROIDAL_SURFACE")]
    return sort!(vcat(breps,primitives);by=ent->ent.id)
end

function _step_root_surface(ent,entities)
    ent.kind in _STEP_BREP_SOLID_KINDS && return _step_brep_solid(ent,entities,false)
    # Retain placement references, but expose only this root to the existing
    # single-primitive classifier. Selecting a root never depends on Dict order.
    selected=Dict(id=>other for (id,other) in entities if id==ent.id ||
        !(other.kind in ("SPHERE","RIGHT_CIRCULAR_CYLINDER","RIGHT_CIRCULAR_CONE","TOROIDAL_SURFACE")))
    if ent.kind=="SPHERE"
        params=_as_sphere(selected)
        params===nothing || return sphere_surface(params...)
    elseif ent.kind=="RIGHT_CIRCULAR_CYLINDER"
        params=_as_cylinder(selected)
        params===nothing || return cylinder_surface(params...)
    elseif ent.kind=="RIGHT_CIRCULAR_CONE"
        params=_as_cone(selected)
        params===nothing || return cone_surface(params...)
    elseif ent.kind=="TOROIDAL_SURFACE"
        params=_step_torus_params(entities,ent,"import_step")
        return _torus_surface(params...,"import_step")
    end
    throw(ArgumentError("import_step: source solid #$(ent.id) $(ent.kind) is truncated"))
end

function _step_assembly_records(roots,entities)
    for ent in values(entities)
        endswith(ent.kind,"SOLID") && !(ent.kind in _STEP_BREP_SOLID_KINDS) &&
            throw(ArgumentError("import_step: unsupported solid root #$(ent.id) $(ent.kind)"))
    end
    factors=_step_member_unit_factors(entities)
    return [(name="STEP #$(ent.id) $(ent.kind)"*
                 (!isempty(ent.args) && ent.args[1] isa String && !isempty(ent.args[1]) ? " "*ent.args[1] : ""),
             surface=_step_unit_scaled(_step_root_surface(ent,entities),
                 get(factors,ent.id,1.0))) for ent in roots]
end

function _iges_primitive_surface(rec)
    type=_as_int(rec[1],"import_iges","entity type")
    if type==150 && length(rec)>=7
        L,W,H,x,y,z=rec[2:7]
        (L>0 && W>0 && H>0) || throw(ArgumentError("import_iges: Block extents must be positive"))
        return box_surface(x,x+L,y,y+W,z,z+H)
    elseif type==158 && length(rec)>=5
        r,x,y,z=rec[2:5]
        r>0 || throw(ArgumentError("import_iges: Sphere radius must be positive"))
        return sphere_surface((x,y,z),r)
    elseif type==156 && length(rec)>=10
        h,r1,r2,x,y,z,zi,zj,zk=rec[2:10]
        (h>0 && (r1>0 || r2>0) && r1>=0 && r2>=0) || throw(ArgumentError(
            "import_iges: Cone height must be positive and at least one radius must be positive"))
        return cone_surface((x,y,z),(zi,zj,zk),r1,r2,h)
    elseif type==156 && length(rec)>=8
        r,h,x,y,z,zi,zj=rec[2:8];zk=length(rec)>=9 ? rec[9] : 1.0
        (r>0 && h>0) || throw(ArgumentError("import_iges: Cylinder radius/height must be positive"))
        return cylinder_surface((x,y,z),(zi,zj,zk),r,h)
    end
    type in (150,156,158) && throw(ArgumentError("import_iges: primitive type $type is truncated"))
    return nothing
end

function _iges_linked_primitive_entries(types,params,bad)
    ids=Set(de for (de,type) in types if type in (150,156,158))
    for (de,rec) in params
        !isempty(rec) && rec[1] in (150,156,158) && push!(ids,de)
    end
    return [(id=de,record=_iges_rec(params,bad,de,"import_iges")) for de in sort!(collect(ids))]
end

function _iges_assembly_directory_certificate(source,types)
    # Admit the implemented solid and support layouts, plus property records.
    # An unknown geometric entity must not silently disappear from an assembly.
    supported=(100,102,108,110,116,120,123,124,141,142,143,144,150,156,158,
               160,186,190,198,406,502,504,508,510,514)
    all(type->type in supported,values(types)) || throw(ArgumentError(
        "import_iges: assembly contains an unsupported directory-entry entity type"))
    seen=Set{Int}();root_records=Dict{Int,Int}()
    for raw in split(source,r"\r?\n")
        line=length(raw)>=80 ? raw[1:80] : rpad(raw,80)
        if line[73]=='D'
            de=_iges_int_field(line[74:80])
            de===nothing && throw(ArgumentError("import_iges: invalid directory-entry sequence"))
            de>0 || throw(ArgumentError("import_iges: directory-entry sequences must be positive"))
            de in seen && throw(ArgumentError("import_iges: duplicate directory-entry sequence #$de"))
            push!(seen,de)
        elseif line[73]=='P'
            de=_iges_int_field(line[65:72])
            get(types,de,0) in (150,156,158,160,186,198) || continue
            occursin(';',line[1:64]) || continue
            root_records[de]=get(root_records,de,0)+1
            root_records[de]==1 || throw(ArgumentError(
                "import_iges: multiple parameter records claim solid directory entry #$de"))
        end
    end
    return nothing
end

function _iges_assembly_records(primitive_entries,types,params,bad,transforms;frame_context=nothing)
    caller="import_iges";records=NamedTuple[]
    msbo=sort!([de for (de,type) in types if type==186])
    trimmed=sort!([de for (de,type) in types if type==144])
    bounded=sort!([de for (de,type) in types if type==143])
    referenced_tori=Set{Int}()
    if !isempty(msbo)
        isempty(trimmed) && isempty(bounded) || throw(ArgumentError(
            "$caller: mixed MSBO and ungrouped trimmed Face roots cannot be assigned independent solid identities"))
        for de in msbo
            faces,_,_=_iges_brep_faces(types,params,bad,[de],Int[],Int[],caller;frame_context)
            union!(referenced_tori,(face[3] for face in faces))
            surface=_iges_brep_surface(types,params,bad,transforms,[de],Int[],Int[],caller;frame_context)
            surface=_iges_assembly_transform(surface,types,params,bad,transforms,de,caller)
            push!(records,(id=de,name="IGES #$de type 186",surface=surface))
        end
    elseif !isempty(trimmed) || !isempty(bounded)
        # The legacy face collection declares one closed shell. Without an
        # MSBO root there is no source identity for splitting it into solids.
        surface=_iges_brep_surface(types,params,bad,transforms,Int[],trimmed,bounded,caller;frame_context)
        faces,_,_=_iges_brep_faces(types,params,bad,Int[],trimmed,bounded,caller;frame_context)
        union!(referenced_tori,(face[3] for face in faces))
        id=minimum(vcat(trimmed,bounded))
        push!(records,(id=id,name="IGES #$id closed shell",surface=surface))
    end
    for de in sort!([de for (de,type) in types if type in (160,198) && !(de in referenced_tori)])
        push!(records,(id=de,name="IGES #$de type $(types[de])",
            surface=_iges_assembly_transform(_iges_torus_surface(types,params,bad,de,caller),
                types,params,bad,transforms,de,caller)))
    end
    for entry in primitive_entries
        rec=entry.record;type=_as_int(rec[1],caller,"entity type")
        if !(type in (150,156,158))
            isempty(types) && !(type in (124,406)) && throw(ArgumentError(
                "$caller: assembly contains unsupported parameter entity type $type"))
            continue
        end
        entry.id!==nothing && entry.id>0 && isodd(entry.id) || throw(ArgumentError(
            "$caller: multiple recognized solids require positive odd IGES directory-entry pointers"))
        isempty(types) || get(types,entry.id,0)==type || throw(ArgumentError(
            "$caller: primitive parameter type and directory entry #$(entry.id) disagree"))
        surface=_iges_primitive_surface(rec)
        surface=_iges_assembly_transform(surface,types,params,bad,transforms,entry.id,caller)
        push!(records,(id=entry.id,name="IGES #$(entry.id) type $type",surface=surface))
    end
    allunique(record.id for record in records) || throw(ArgumentError(
        "$caller: multiple recognized solids require unique IGES directory-entry pointers"))
    all(record->record.id>0 && isodd(record.id),records) || throw(ArgumentError(
        "$caller: solid roots require positive odd IGES directory-entry identities"))
    for record in records
        haskey(types,record.id) || continue
        parameters=_iges_rec(params,bad,record.id,caller)
        !isempty(parameters) && _as_int(parameters[1],caller,"solid parameter type")==types[record.id] ||
            throw(ArgumentError("$caller: solid parameter type and directory entry #$(record.id) disagree"))
    end
    return sort!(records;by=record->record.id)
end
