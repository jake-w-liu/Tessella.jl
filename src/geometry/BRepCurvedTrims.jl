# ── curved trims: full-tube cylinders + circle-trimmed planar faces ───────────
#
# A STEP closed shell whose faces are PLANEs bounded by LINE wires and closed
# CIRCLE edges, together with full-tube CYLINDRICAL_SURFACE walls, covers the
# drilled-plate / annular-cap / boss-on-plate family: a plate with one or more
# circular through-holes, a hollow cylinder, or a boss on a plate. Every
# CIRCLE edge must be a closed rim (its EDGE_CURVE's two vertices coincide)
# coaxial with a certified tube; every tube rim must be consumed by exactly
# one wire of exactly one planar face; every LINE edge closes by position.
# Partial arcs, B-spline trims, slit walls and non-coaxial circles remain
# explicit blockers routed through the shared curve-kind diagnostic.
#
# Emission welds one canonical rim grid per tube rim — the cap's hole wire and
# the wall's band edge coincide bitwise — so `_dedup_positions` identifies
# shared boundary vertices exactly and the assembled triangle complex is
# verified by a directed shared-edge audit (each undirected edge traversed
# once in each direction) plus the standard positive-volume and `validate`
# certificates. Orientation is recovered geometrically by propagating a flip
# parity across face-adjacent triangles — independent of STEP sense flags —
# then normalizing the single component's signed volume, mirroring
# `_polyhedral_surface`'s component rule.

const _BREP_RIM_SEGMENTS=24

# Traversal-ordered uses of one EDGE_LOOP wire for the tubed-shell
# certificate. Mirrors _step_face_bound_uses for a single bound entity.
function _step_tube_wire_uses(entities,be,curvekinds,caller)
    loop=_deref_entity(entities,be.args[2])
    (loop!==nothing && loop.kind=="EDGE_LOOP") || throw(ArgumentError(
        "$caller: curved-trim certificate requires EDGE_LOOP bounds; saw "*
        "$(loop===nothing ? "a missing loop" : loop.kind)"))
    j=_skip_name(loop.args)
    length(loop.args)>=j || throw(ArgumentError(
        "$caller: EDGE_LOOP #$(loop.id) is truncated"))
    refs=loop.args[j]
    refs isa Vector || throw(ArgumentError(
        "$caller: EDGE_LOOP #$(loop.id) member list is malformed"))
    isempty(refs) && throw(ArgumentError("$caller: EDGE_LOOP #$(loop.id) is empty"))
    uses=_BREP_EDGE_USE[]
    for r in refs
        oe=_deref_entity(entities,r)
        (oe!==nothing && oe.kind=="ORIENTED_EDGE") || throw(ArgumentError(
            "$caller: EDGE_LOOP #$(loop.id) member is not an ORIENTED_EDGE"))
        ecref=nothing
        for a in oe.args
            a isa Tuple && a[1]===:ref && (ecref=a)
        end
        ec=_deref_entity(entities,ecref)
        (ec!==nothing && ec.kind=="EDGE_CURVE") || throw(ArgumentError(
            "$caller: ORIENTED_EDGE #$(oe.id) does not reference an EDGE_CURVE"))
        (v1,v2),_,curve=_step_oriented_edge(entities,oe,curvekinds,caller)
        base=_step_curve_base(entities,curve,0)
        base===nothing && throw(ArgumentError(
            "$caller: EDGE_CURVE geometry of ORIENTED_EDGE #$(oe.id) is unresolvable"))
        push!(uses,(ec=ec,v1=v1,v2=v2,base=base))
    end
    for i in 1:length(uses)-1
        uses[i].v2==uses[i+1].v1 || throw(ArgumentError(
            "$caller: face boundary edges do not form an ordered closed chain"))
    end
    uses[end].v2==uses[1].v1 || throw(ArgumentError(
        "$caller: face boundary loop does not close"))
    return uses
end

# Canonical rim polygon: nθ points on the rim circle at height h of a tube,
# sampled in the tube's own (x,y) frame starting at the x axis. Cap hole wires
# and wall band edges share this list, welding bitwise on dedup.
function _step_tube_rim_grid(c0,z,x,y,R,h,n)
    C=(c0[1]+h*z[1],c0[2]+h*z[2],c0[3]+h*z[3])
    return [ntuple(ax->C[ax]+R*(cos(2π*k/n)*x[ax]+sin(2π*k/n)*y[ax]),3)
            for k in 0:n-1]
end

# Certify one CYLINDRICAL_SURFACE face as a full tube: its bounds are exactly
# two closed coaxial rim circles at distinct heights plus optional axial seam
# LINEs. Returns (c0,z,x,y,R,h_lo,h_hi) or nothing when the face is not a
# full-tube layout (open rim arcs belong to the trimmed-sector certificate);
# throws on certified mismatches.
function _step_tube_certify(entities,face,curvekinds,caller)
    c0,z,x,y,R=_step_cylinder_frame(entities,face.surf,caller)
    tol=1e-9*max(R,1.0); atol=1e-9
    fe=_deref_entity(entities,(:ref,face.id))
    fe===nothing && throw(ArgumentError(
        "$caller: cylindrical face #$(face.id) record is missing"))
    i=_skip_name(fe.args)
    bounds=fe.args[i]
    bounds isa Vector || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) boundary list is malformed"))
    isempty(bounds) && throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has no face boundaries"))
    circles=_BREP_EDGE_USE[]; seams=_BREP_EDGE_USE[]
    for bref in bounds
        be=_deref_entity(entities,bref)
        (be!==nothing && be.kind in _STEP_BREP_BOUND_KINDS) || throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) boundary is not a FACE_BOUND/FACE_OUTER_BOUND"))
        length(be.args)>=3 || throw(ArgumentError(
            "$caller: $(be.kind) #$(be.id) is truncated"))
        for u in _step_tube_wire_uses(entities,be,curvekinds,caller)
            if u.base.kind=="CIRCLE"
                u.v1==u.v2 || return nothing   # open rim arc: sector certificate
                push!(circles,u)
            elseif u.base.kind=="LINE"
                u.v1!=u.v2 || throw(ArgumentError(
                    "$caller: closed LINE edge #$(u.ec.id) of tube face "*
                    "#$(face.id) is not an axial seam"))
                push!(seams,u)
            else
                push!(curvekinds,u.base.kind)
                return nothing
            end
        end
    end
    length(circles)==2 || return nothing
    hs=Float64[]
    for u in circles
        cc,cn,_,_,cr=_step_circle_frame_full(entities,u.base,caller)
        abs(cr-R)<=tol || throw(ArgumentError(
            "$caller: rim CIRCLE #$(u.base.id) radius does not match the "*
            "cylindrical surface"))
        _brep_axis_dist(cc,c0,z)<=tol || throw(ArgumentError(
            "$caller: rim CIRCLE #$(u.base.id) center is off the cylinder axis"))
        abs(_brep_vdot(cn,z))>=1.0-atol || throw(ArgumentError(
            "$caller: rim CIRCLE #$(u.base.id) is not perpendicular to the "*
            "cylinder axis"))
        abs(sqrt(_brep_vdot(_brep_vsub(u.v1,cc),_brep_vsub(u.v1,cc)))-cr)<=tol ||
            throw(ArgumentError(
                "$caller: rim CIRCLE #$(u.base.id) vertex is off the circle"))
        push!(hs,_brep_vdot(_brep_vsub(cc,c0),z))
    end
    abs(hs[1]-hs[2])>tol || return nothing
    h_lo=min(hs[1],hs[2]); h_hi=max(hs[1],hs[2])
    for u in seams
        d=_brep_vsub(u.v2,u.v1)
        abs(_brep_vdot(_brep_vunit(d,caller),z))>=1.0-atol || throw(ArgumentError(
            "$caller: seam edge #$(u.ec.id) of tube face #$(face.id) is not "*
            "parallel to the cylinder axis"))
        for v in (u.v1,u.v2)
            abs(_brep_axis_dist(v,c0,z)-R)<=tol || throw(ArgumentError(
                "$caller: seam edge #$(u.ec.id) endpoint is off the "*
                "cylindrical surface"))
            t=_brep_vdot(_brep_vsub(v,c0),z)
            (abs(t-h_lo)<=tol || abs(t-h_hi)<=tol) ||
                throw(ArgumentError(
                    "$caller: seam edge #$(u.ec.id) of tube face #$(face.id) "*
                    "does not reach either rim circle"))
        end
        t1=_brep_vdot(_brep_vsub(u.v1,c0),z); t2=_brep_vdot(_brep_vsub(u.v2,c0),z)
        (abs(t1-h_lo)<=tol && abs(t2-h_hi)<=tol ||
         abs(t1-h_hi)<=tol && abs(t2-h_lo)<=tol) || throw(ArgumentError(
            "$caller: seam edge #$(u.ec.id) of tube face #$(face.id) does "*
            "not join the two rims"))
    end
    return (c0=c0,z=z,x=x,y=y,R=R,h=(h_lo,h_hi),claimed=falses(2))
end

# A certified full-tube + circle-trimmed planar shell. Returns a Mesh or
# nothing when the face inventory is not this shape.
function _step_tubed_shell(entities,faces,curvekinds,caller)
    length(faces)>=4 || return nothing
    all(f->f.role==1,faces) || return nothing
    all(f->f.surf!==nothing && f.rings!==nothing,faces) || return nothing
    cyl_i=Int[]; plane_i=Int[]
    for (i,f) in enumerate(faces)
        k=f.surf.kind
        if k=="CYLINDRICAL_SURFACE"
            push!(cyl_i,i)
        elseif k=="PLANE"
            push!(plane_i,i)
        else
            return nothing
        end
    end
    isempty(cyl_i) && return nothing

    tubes=[]
    for i in cyl_i
        t=_step_tube_certify(entities,faces[i],curvekinds,caller)
        t===nothing && return nothing
        push!(tubes,(;face=i,t...))
    end

    # ── planar faces: LINE wires and closed-circle rim wires ──────────────
    # per face: (:lines|:circle, ring, tube_index, rim_slot)
    face_wires=Vector{Vector{NamedTuple}}(undef,length(faces))
    for i in plane_i
        f=faces[i]
        p0,pn=_step_plane_frame(entities,f.surf,caller)
        ptol=1e-9*max(_poly_scale(vcat(f.rings...)),1.0)
        fe=_deref_entity(entities,(:ref,f.id))
        fe===nothing && throw(ArgumentError(
            "$caller: face #$(f.id) record is missing"))
        i0=_skip_name(fe.args)
        bounds=fe.args[i0]
        bounds isa Vector || throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) boundary list is malformed"))
        isempty(bounds) && throw(ArgumentError(
            "$caller: $(fe.kind) #$(fe.id) has no face boundaries"))
        wires=NamedTuple[]
        for bref in bounds
            be=_deref_entity(entities,bref)
            (be!==nothing && be.kind in _STEP_BREP_BOUND_KINDS) || throw(ArgumentError(
                "$caller: $(fe.kind) #$(fe.id) boundary is not a FACE_BOUND/FACE_OUTER_BOUND"))
            length(be.args)>=3 || throw(ArgumentError(
                "$caller: $(be.kind) #$(be.id) is truncated"))
            uses=_step_tube_wire_uses(entities,be,curvekinds,caller)
            if length(uses)==1 && only(uses).base.kind=="CIRCLE" &&
               only(uses).v1==only(uses).v2
                u=only(uses)
                cc,cn,_,_,cr=_step_circle_frame_full(entities,u.base,caller)
                match=0; slot=0
                for (ti,t) in enumerate(tubes)
                    ttol=1e-9*max(t.R,1.0)
                    abs(cr-t.R)<=ttol || continue
                    _brep_axis_dist(cc,t.c0,t.z)<=ttol || continue
                    abs(_brep_vdot(cn,t.z))>=1.0-1e-9 || continue
                    h=_brep_vdot(_brep_vsub(cc,t.c0),t.z)
                    for ri in 1:2
                        if abs(h-t.h[ri])<=ttol
                            match==0 || throw(ArgumentError(
                                "$caller: CIRCLE #$(u.base.id) matches "*
                                "multiple certified tube rims"))
                            match=ti; slot=ri
                        end
                    end
                end
                match==0 && (push!(curvekinds,"CIRCLE"); return nothing)
                t=tubes[match]; ttol=1e-9*max(t.R,1.0)
                # A rim claimed by two planar bounds is a T-junction or a
                # duplicated cap: not the tubed-shell layout, so the
                # trimmed-sector certificate reports it instead.
                t.claimed[slot] && return nothing
                # The cap's PLANE must contain the rim circle.
                abs(_brep_vdot(pn,t.z))>=1.0-1e-9 || throw(ArgumentError(
                    "$caller: rim CIRCLE #$(u.base.id) is not parallel to "*
                    "PLANE #$(f.surf.id)"))
                rimC=(t.c0[1]+t.h[slot]*t.z[1],t.c0[2]+t.h[slot]*t.z[2],
                      t.c0[3]+t.h[slot]*t.z[3])
                abs(_brep_vdot(_brep_vsub(rimC,p0),pn))<=ttol || throw(ArgumentError(
                    "$caller: PLANE #$(f.surf.id) does not contain its "*
                    "rim circle #$(u.base.id)"))
                abs(sqrt(_brep_vdot(_brep_vsub(u.v1,cc),
                                    _brep_vsub(u.v1,cc)))-cr)<=ttol ||
                    throw(ArgumentError(
                        "$caller: rim CIRCLE #$(u.base.id) vertex is off the circle"))
                claimed=copy(t.claimed); claimed[slot]=true
                tubes[match]=(;t...,claimed=claimed)
                ring=_step_tube_rim_grid(t.c0,t.z,t.x,t.y,t.R,t.h[slot],
                                         _BREP_RIM_SEGMENTS)
                push!(wires,(;kind=:circle,ring,rim=(match,slot)))
            else
                all(u->u.base.kind=="LINE",uses) ||
                    (push!(curvekinds,first(u.base.kind for u in uses
                                           if u.base.kind!="LINE"));
                     return nothing)
                ring=NTuple{3,Float64}[u.v1 for u in uses]
                _check_coplanar(ring,pn,p0,ptol,f.id,caller)
                push!(wires,(;kind=:lines,ring))
            end
        end
        face_wires[i]=wires
    end
    all(t->all(t.claimed),tubes) || return nothing

    # ── emission ──────────────────────────────────────────────────────────
    pool=NTuple{3,Float64}[]
    items=Tuple{Symbol,UnitRange{Int},Int}[]
    emit!(kind,pts,fi)=(push!(items,(kind,
        length(pool)+1:length(pool)+length(pts),fi));append!(pool,pts);nothing)
    for i in plane_i
        wires=face_wires[i]
        rings=[w.ring for w in wires]
        if length(rings)>1
            # order the containing wire first regardless of bound order, then
            # normalize inner winding to oppose it for the CDT domain.
            nrm=_newell_normal(rings[argmax(length.(rings))])
            drop=argmax(abs.(nrm));axes=filter(!=(drop),collect(1:3))
            projected=[[(p[axes[1]],p[axes[2]]) for p in ring] for ring in rings]
            outer=0
            for k in eachindex(projected)
                inside=false
                for l in eachindex(projected)
                    k==l && continue
                    if _brep_point_in_wire(first(projected[k]),projected[l])
                        inside=true;break
                    end
                end
                if !inside
                    outer==0 || throw(ArgumentError(
                        "$caller: face #$(faces[i].id) has more than one "*
                        "outer boundary wire"))
                    outer=k
                end
            end
            outer==0 && throw(ArgumentError(
                "$caller: face #$(faces[i].id) has no containing outer wire"))
            rings=vcat(rings[outer:outer],
                       rings[setdiff(eachindex(rings),outer)])
            a1=_signed_area2([p[axes[1]] for p in rings[1]],
                             [p[axes[2]] for p in rings[1]])
            for k in 2:length(rings)
                ak=_signed_area2([p[axes[1]] for p in rings[k]],
                                 [p[axes[2]] for p in rings[k]])
                signbit(ak)==signbit(a1) && reverse!(rings[k])
            end
            for tri in _brep_planar_face_holes(rings,caller)
                emit!(:tri,tri,i)
            end
        else
            emit!(:ring,only(rings),i)
        end
    end
    for t in tubes
        lo=_step_tube_rim_grid(t.c0,t.z,t.x,t.y,t.R,t.h[1],_BREP_RIM_SEGMENTS)
        hi=_step_tube_rim_grid(t.c0,t.z,t.x,t.y,t.R,t.h[2],_BREP_RIM_SEGMENTS)
        n=length(lo)
        for k in 1:n
            k2=mod(k,n)+1
            emit!(:quad,(lo[k],lo[k2],hi[k2],hi[k]),t.face)
        end
    end
    nodes,canon=_dedup_positions(pool,caller)
    itris=NTuple{3,Int32}[]; tface=Int[]
    for (kind,rng,fi) in items
        if kind==:ring
            ir=Int32[]; prev=Int32(0)
            for j in rng
                v=canon[j]
                v!=prev && push!(ir,v); prev=v
            end
            length(ir)>1 && first(ir)==last(ir) && pop!(ir)
            n0=length(itris)
            _ear_clip!(itris,ir,nodes,caller)
            append!(tface,fill(fi,length(itris)-n0))
        elseif kind==:quad
            a,b,c,d=(Int32(canon[j]) for j in rng)
            push!(itris,(a,b,c));push!(itris,(a,c,d))
            append!(tface,(fi,fi))
        else
            a,b,c=(Int32(canon[j]) for j in rng)
            push!(itris,(a,b,c));push!(tface,fi)
        end
    end
    _checked_brep_mul(3,length(itris),caller,"tubed shell triangle count")
    length(itris)<=_BREP_MAX_SURFACE_TRIS || throw(ArgumentError(
        "$caller: tubed shell exceeds the surface triangle limit"))

    # Orientation: flip-parity propagation over the face-adjacency graph,
    # then the shared directed-edge/volume/validate closing certificate.
    tris,edge_use=_step_shell_flip_parity(itris,tface,length(faces),
        "tubed shell",caller)
    return _step_shell_audit(tris,nodes,tface,length(faces),edge_use,
        "tubed shell",caller)
end

# Shared closing certificate for emitted curved-trim shells: a directed
# shared-edge audit (every undirected edge traversed once in each direction),
# positive signed-volume normalization and `validate`. `tris` are already
# consistently oriented by the caller's flip-parity propagation over `tface`.
function _step_shell_audit(tris,nodes,tface,nfaces,edge_use,what,caller)
    audit=Dict{Tuple{Int32,Int32},Int}()
    for t in tris,(a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        audit[(a,b)]=get(audit,(a,b),0)+1
    end
    for (key,(t1,t2)) in edge_use
        ab=get(audit,(key[1],key[2]),0);ba=get(audit,(key[2],key[1]),0)
        (ab==1 && ba==1) || throw(ArgumentError(
            "$caller: $what is inconsistently oriented across edge $key"))
    end
    scale=max(_poly_scale(nodes),1.0)
    V=sum(_brep_vdot(nodes[t[1]],
        _brep_vcross(nodes[t[2]],nodes[t[3]]))/6 for t in tris)
    abs(V)>1e-12*scale*scale*scale || throw(ArgumentError(
        "$caller: $what bounds zero volume"))
    V<0 && (tris=map(t->(t[1],t[3],t[2]),tris))
    C=Matrix{Float64}(undef,3,length(nodes))
    for (i,p) in enumerate(nodes); C[:,i]=[p...]; end
    T=Matrix{Int32}(undef,3,length(tris))
    for (i,t) in enumerate(tris); T[:,i]=[t...]; end
    surface=Mesh(C;tris=T)
    diag=validate(surface)
    diag.ok || throw(ErrorException(
        "$caller: $what boundary mesh is invalid — "*
        join(diag.messages,"; ")))
    return surface
end

# Shared flip-parity propagation over the face-adjacency graph implied by the
# undirected edge uses of the emitted triangles: incident faces orient
# consistently iff their flip parities differ exactly when the two edge
# traversals already disagree. Returns (tris, edge_use) where edge_use maps
# each undirected edge to its incident triangle pair for the closing audit.
function _step_shell_flip_parity(itris,tface,nfaces,what,caller)
    edge_use=Dict{Tuple{Int32,Int32},NTuple{2,Int}}()
    edge_dir=Dict{Tuple{Int32,Int32},NTuple{2,Bool}}()
    for (ti,t) in enumerate(itris),(a,b) in
            ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        key=a<b ? (a,b) : (b,a)
        if haskey(edge_use,key)
            prev=edge_use[key]; pd=edge_dir[key]
            edge_use[key]=(prev[1],ti);edge_dir[key]=(pd[1],a==key[1])
        else
            edge_use[key]=(ti,0);edge_dir[key]=(a==key[1],false)
        end
    end
    adj=[Int[] for _ in 1:nfaces]; diff=Dict{Tuple{Int,Int},Bool}()
    for (key,(t1,t2)) in edge_use
        t2==0 && throw(ArgumentError(
            "$caller: $what is not closed — edge $key bounds a "*
            "single triangle"))
        f1=tface[t1];f2=tface[t2]
        f1==f2 && continue
        same=edge_dir[key][1]==edge_dir[key][2]
        k=f1<f2 ? (f1,f2) : (f2,f1)
        want=get(diff,k,nothing)
        if want===nothing
            diff[k]=same
            push!(adj[f1],f2);push!(adj[f2],f1)
        elseif want!=same
            throw(ArgumentError(
                "$caller: $what faces $f1 and $f2 cannot be "*
                "consistently oriented across their shared edges"))
        end
    end
    flip=falses(nfaces);seen=falses(nfaces)
    for seed in 1:nfaces
        seen[seed] && continue
        seen[seed]=true;queue=[seed]
        while !isempty(queue)
            f=popfirst!(queue)
            for g in adj[f]
                k=f<g ? (f,g) : (g,f)
                want=diff[k]
                need=flip[f]⊻want
                if seen[g]
                    flip[g]==need || throw(ArgumentError(
                        "$caller: $what is not orientable — "*
                        "inconsistent winding across face #$g"))
                else
                    flip[g]=need;seen[g]=true;push!(queue,g)
                end
            end
        end
    end
    tris=NTuple{3,Int32}[]
    for (ti,t) in enumerate(itris)
        push!(tris,flip[tface[ti]] ? (t[1],t[3],t[2]) : t)
    end
    return tris,edge_use
end

# ── prismatic shells: two-cap extrusions with line/arc/polyline trims ────────
#
# A STEP closed shell covering a strict prism — the volume swept along one
# axial direction z between two parallel caps — where cap wires admit LINE
# edges, open CIRCLE arcs, closed CIRCLE edges and degree-1 B-spline
# (polyline-equivalent) edges, and wall faces are z-parallel PLANEs (single
# LINE/polyline loops spanning both cap heights) or z-parallel
# CYLINDRICAL_SURFACEs (rim arcs joined by axial seams, both rims covering the
# same angular window). A cylinder face selects the axis; cap planes are
# perpendicular to it at exactly two heights.
#
# The cross-section is certified edge-by-edge: every straight segment — a
# LINE use or one piece of a degree-1 B-spline sub-chain — is shared by
# exactly two uses across the shell, and every circular edge is covered
# exactly twice on its geometric circle (cap wire against wall rim, or two
# coplanar cap faces). Each cylindrical wall's two rims must tile the same
# angular window so the swept band pairs knot-for-knot. Emission welds one
# canonical sampled chain per geometric edge piece — angle-keyed in the
# shared circle frame — so `_dedup_positions` identifies shared boundary
# vertices bitwise, then the shared flip-parity, directed-edge, volume and
# `validate` certificates close the shell. Unsupported edge kinds (ELLIPSE,
# higher-degree or rational splines, TRIMMED_CURVE remnants) and unsupported
# wall surfaces stay explicit blockers through `curvekinds`.

const _BREP_PRISM_ARC_SEGMENTS=48

@inline _posnear(p,q,tol)=abs(p[1]-q[1])<=tol &&
    abs(p[2]-q[2])<=tol && abs(p[3]-q[3])<=tol

# Canonical window of an open CIRCLE-trimmed edge in its own circle frame:
# (center, nr, x, y, r, span, start, stop) where the covered CCW window is
# [theta(start), theta(start)+span] and (start, stop) are the vertex
# positions at the window ends. The EDGE_CURVE same_sense flag selects the
# complement, mirroring _step_cyl_rim_arc's convention.
function _step_arc_span(entities,use,caller)
    ec=use.ec
    length(ec.args)>=3 || throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) is truncated"))
    vs=_step_vertex_point(entities,ec.args[2])
    ve=_step_vertex_point(entities,ec.args[3])
    (vs===nothing || ve===nothing) && throw(ArgumentError(
        "$caller: EDGE_CURVE #$(ec.id) has unresolvable vertex geometry"))
    vs!=ve || throw(ArgumentError(
        "$caller: closed CIRCLE edge #$(ec.id) reached the arc-span "*
        "certificate"))
    sense=length(ec.args)>=5 ? _brep_enum(ec.args[5]) : nothing
    sense===nothing && (sense=true)
    center,nr,xr,yr,cr=_step_circle_frame_full(entities,use.base,caller)
    tol=1e-9*max(cr,1.0)
    for v in (vs,ve)
        abs(sqrt(_brep_vdot(_brep_vsub(v,center),_brep_vsub(v,center)))-cr)<=tol ||
            throw(ArgumentError(
                "$caller: arc EDGE_CURVE #$(ec.id) endpoint is off its circle"))
    end
    ta=atan(_brep_vdot(_brep_vsub(vs,center),yr),
            _brep_vdot(_brep_vsub(vs,center),xr))
    tb=atan(_brep_vdot(_brep_vsub(ve,center),yr),
            _brep_vdot(_brep_vsub(ve,center),xr))
    span=sense ? mod(tb-ta,2π) : mod(ta-tb,2π)
    (span>1e-9 && span<2π-1e-9) || throw(ArgumentError(
        "$caller: arc EDGE_CURVE #$(ec.id) is degenerate"))
    return (center=center,nr=nr,x=xr,y=yr,r=cr,span=span,
            start=(sense ? vs : ve),stop=(sense ? ve : vs))
end

# Degree-1 B-spline trims are exact polylines through the control points; the
# edge's vertex pair must land on a consecutive control sub-chain. Returns
# the sub-chain in control order or nothing for other curve kinds.
function _step_polyline_subchain(entities,use)
    c=_step_nurbs_curve(use.base,entities)
    (c===nothing || c.degree!=1) && return nothing
    ec=use.ec
    length(ec.args)>=3 || throw(ArgumentError(
        "_step_prismatic_shell: EDGE_CURVE #$(ec.id) is truncated"))
    vs=_step_vertex_point(entities,ec.args[2])
    ve=_step_vertex_point(entities,ec.args[3])
    (vs===nothing || ve===nothing) && throw(ArgumentError(
        "_step_prismatic_shell: EDGE_CURVE #$(ec.id) has unresolvable "*
        "vertex geometry"))
    pts=c.controls
    tol=1e-9*max(_poly_scale(pts),1.0)
    a=0;b=0
    for (i,p) in enumerate(pts)
        _posnear(p,vs,tol) && (a=i)
        _posnear(p,ve,tol) && (b=i)
    end
    (a>0 && b>0 && a!=b) || return nothing
    return a<b ? pts[a:b] : pts[b:a]
end

function _step_prismatic_shell(entities,faces,curvekinds,caller)
    length(faces)>=4 || return nothing
    all(f->f.role==1,faces) || return nothing
    all(f->f.surf!==nothing && f.rings!==nothing,faces) || return nothing
    atol=1e-9
    # Extrusion axis: every cylindrical face must be parallel to the first.
    z=nothing; cyl_f=Dict{Int,NamedTuple}()
    for (i,f) in enumerate(faces)
        f.surf.kind=="CYLINDRICAL_SURFACE" || continue
        c0,zc,xc,yc,R=_step_cylinder_frame(entities,f.surf,caller)
        if z===nothing
            z=zc
        else
            abs(_brep_vdot(zc,z))>=1.0-atol || return nothing
        end
        cyl_f[i]=(c0=c0,z=zc,x=xc,y=yc,R=R)
    end
    z===nothing && return nothing
    # Faces: caps are PLANEs perpendicular to z; walls are PLANEs parallel to
    # z or the parallel cylinders. Anything else is not this family.
    cap_i=Int[];wall_i=Int[];cap_hs=Float64[];cap_plane=Dict{Int,Tuple}()
    for (i,f) in enumerate(faces)
        k=f.surf.kind
        if k=="PLANE"
            p0,pn=_step_plane_frame(entities,f.surf,caller)
            pz=_brep_vdot(pn,z)
            if abs(pz)>=1.0-atol
                push!(cap_i,i);push!(cap_hs,_brep_vdot(p0,z))
                cap_plane[i]=(p0,pn)
            elseif abs(pz)<=atol
                push!(wall_i,i)
            else
                return nothing
            end
        elseif k=="CYLINDRICAL_SURFACE"
            push!(wall_i,i)
        else
            push!(curvekinds,"face surface $(_step_surface_name(f.surf))")
            return nothing
        end
    end
    (length(cap_i)>=2 && !isempty(wall_i)) || return nothing
    h0=minimum(cap_hs);h1=maximum(cap_hs)
    h1-h0>atol || return nothing
    hscale=max(abs(h0),abs(h1),1.0);htol=1e-9*hscale
    all(h->abs(h-h0)<=htol || abs(h-h1)<=htol,cap_hs) || return nothing

    # ── per-face classified uses ──────────────────────────────────────────
    # fuse[i][bound] = classified uses in wire order:
    #   (:line,u,nothing) | (:polyline,u,(pts=...,)) |
    #   (:arc,u,arc-span) | (:circle,u,circle-frame)
    fuse=Vector{Vector{Vector{NamedTuple}}}(undef,length(faces))
    for (i,f) in enumerate(faces)
        fe=_deref_entity(entities,(:ref,f.id))
        fe===nothing && throw(ArgumentError(
            "$caller: face #$(f.id) record is missing"))
        bounds=_step_face_bound_uses(entities,fe,curvekinds,caller)
        i in cap_i || length(bounds)==1 || return nothing
        fused=Vector{NamedTuple}[]
        for uses in bounds
            out=NamedTuple[]
            for u in uses
                b=u.base.kind
                if b=="LINE"
                    u.v1!=u.v2 || throw(ArgumentError(
                        "$caller: degenerate LINE edge #$(u.ec.id)"))
                    push!(out,(kind=:line,u=u,aux=nothing))
                elseif b=="CIRCLE"
                    if u.v1==u.v2
                        c,nr,xr,yr,cr=_step_circle_frame_full(entities,
                            u.base,caller)
                        ctol=1e-9*max(cr,1.0)
                        abs(sqrt(_brep_vdot(_brep_vsub(u.v1,c),
                            _brep_vsub(u.v1,c)))-cr)<=ctol || throw(ArgumentError(
                            "$caller: CIRCLE #$(u.base.id) vertex is off "*
                            "the circle"))
                        push!(out,(kind=:circle,u=u,
                            aux=(center=c,nr=nr,x=xr,y=yr,r=cr)))
                    else
                        push!(out,(kind=:arc,u=u,
                            aux=_step_arc_span(entities,u,caller)))
                    end
                elseif b=="B_SPLINE_CURVE_WITH_KNOTS"
                    ch=_step_polyline_subchain(entities,u)
                    ch===nothing && (push!(curvekinds,b);return nothing)
                    push!(out,(kind=:polyline,u=u,aux=(pts=ch,)))
                else
                    push!(curvekinds,b);return nothing
                end
            end
            push!(fused,out)
        end
        fuse[i]=fused
    end

    # Scale reference for positional tolerances.
    allv=NTuple{3,Float64}[]
    for i in eachindex(faces),bound in fuse[i],ur in bound
        push!(allv,ur.u.v1);push!(allv,ur.u.v2)
    end
    pscale=max(_poly_scale(allv),1.0);ptol=1e-9*pscale

    # ── straight-segment incidence: every straight edge exactly twice ─────
    vpos=NTuple{3,Float64}[];vranges=Dict{Tuple{Int,Int,Int},UnitRange{Int}}()
    for (i,bound) in enumerate(fuse),(bi,b) in enumerate(bound),
            (ui,ur) in enumerate(b)
        if ur.kind==:line
            r0=length(vpos)+1;push!(vpos,ur.u.v1,ur.u.v2)
            vranges[(i,bi,ui)]=r0:length(vpos)
        elseif ur.kind==:polyline
            r0=length(vpos)+1;append!(vpos,ur.aux.pts)
            vranges[(i,bi,ui)]=r0:length(vpos)
        end
    end
    vnodes,vcanon=_dedup_positions(vpos,caller)
    segcount=Dict{Tuple{Int32,Int32},Int}()
    for (key,rng) in vranges
        ids=[vcanon[j] for j in rng]
        for k in 1:length(ids)-1
            a,b=ids[k],ids[k+1]
            a==b && throw(ArgumentError(
                "$caller: face #$(faces[key[1]].id) carries a degenerate "*
                "edge segment"))
            e=a<b ? (a,b) : (b,a)
            segcount[e]=get(segcount,e,0)+1
        end
    end
    for (e,c) in segcount
        # Under/over-shared straight edges mean the shell is not this
        # certificate's family; defer to the next certificate's diagnostic.
        c==2 || return nothing
    end

    # ── cap faces: every use in-plane ─────────────────────────────────────
    for i in cap_i
        p0,pn=cap_plane[i]
        for bound in fuse[i],ur in bound
            if ur.kind==:polyline
                for v in ur.aux.pts
                    abs(_brep_vdot(_brep_vsub(v,p0),pn))<=ptol || throw(
                        ArgumentError("$caller: polyline edge #$(ur.u.ec.id) "*
                            "leaves PLANE #$(faces[i].surf.id)"))
                end
            else
                for v in (ur.u.v1,ur.u.v2)
                    abs(_brep_vdot(_brep_vsub(v,p0),pn))<=ptol || throw(
                        ArgumentError("$caller: edge #$(ur.u.ec.id) leaves "*
                            "PLANE #$(faces[i].surf.id)"))
                end
            end
            if ur.kind in (:arc,:circle)
                a=ur.aux
                abs(_brep_vdot(a.nr,pn))>=1.0-atol || throw(ArgumentError(
                    "$caller: CIRCLE #$(ur.u.base.id) is not in PLANE "*
                    "#$(faces[i].surf.id)"))
                abs(_brep_vdot(_brep_vsub(a.center,p0),pn))<=ptol || throw(
                    ArgumentError("$caller: CIRCLE #$(ur.u.base.id) center "*
                        "is off PLANE #$(faces[i].surf.id)"))
            end
        end
    end

    # ── wall faces: single loop, axial or cap-height edges ────────────────
    wallc=Dict{Int,NamedTuple}()   # cyl walls: (c0,z,x,y,R,arcs,pieces,g0i,g1i)
    for i in wall_i
        f=faces[i]
        uses=only(fuse[i])
        if f.surf.kind=="CYLINDRICAL_SURFACE"
            cf=cyl_f[i];tol=1e-9*max(cf.R,1.0)
            arcs=NamedTuple[];axial=false
            for (ui,ur) in enumerate(uses)
                if ur.kind==:arc
                    _step_cyl_rim_arc(entities,ur.u,cf.c0,cf.z,cf.x,cf.y,
                        cf.R,tol,atol,caller)
                    hh=_brep_vdot(ur.aux.center,z)
                    rim=abs(hh-h0)<=htol ? 1 :
                        abs(hh-h1)<=htol ? 2 : return nothing
                    push!(arcs,(rim=rim,ui=ui,ur=ur))
                elseif ur.kind==:circle
                    a=ur.aux
                    abs(a.r-cf.R)<=tol || throw(ArgumentError(
                        "$caller: rim CIRCLE #$(ur.u.base.id) radius does "*
                        "not match the cylindrical surface"))
                    _brep_axis_dist(a.center,cf.c0,cf.z)<=tol || throw(
                        ArgumentError("$caller: rim CIRCLE #$(ur.u.base.id) "*
                            "center is off the cylinder axis"))
                    abs(_brep_vdot(a.nr,cf.z))>=1.0-atol || throw(
                        ArgumentError("$caller: rim CIRCLE #$(ur.u.base.id) "*
                            "is not perpendicular to the cylinder axis"))
                    hh=_brep_vdot(a.center,z)
                    rim=abs(hh-h0)<=htol ? 1 : abs(hh-h1)<=htol ? 2 :
                        return nothing
                    push!(arcs,(rim=rim,ui=ui,ur=ur))
                else
                    pts=ur.kind==:polyline ? ur.aux.pts :
                        (ur.u.v1,ur.u.v2)
                    for v in pts
                        abs(_brep_axis_dist(v,cf.c0,cf.z)-cf.R)<=tol || throw(
                            ArgumentError("$caller: seam edge #$(ur.u.ec.id) "*
                                "vertex is off the cylindrical surface"))
                        t=_brep_vdot(v,z)
                        (abs(t-h0)<=htol || abs(t-h1)<=htol) ||
                            return nothing
                    end
                    for k in 1:length(pts)-1
                        d=_brep_vsub(pts[k+1],pts[k])
                        abs(_brep_vdot(_brep_vunit(d,caller),cf.z))>=1.0-atol ||
                            throw(ArgumentError(
                                "$caller: seam edge #$(ur.u.ec.id) is not "*
                                "parallel to the cylinder axis"))
                    end
                    t1=_brep_vdot(first(pts),z)
                    t2=_brep_vdot(last(pts),z)
                    (abs(t1-h0)<=htol && abs(t2-h1)<=htol ||
                     abs(t1-h1)<=htol && abs(t2-h0)<=htol) || return nothing
                    axial=true
                end
            end
            !isempty(arcs) || return nothing
            wallc[i]=(c0=cf.c0,z=cf.z,x=cf.x,y=cf.y,R=cf.R,arcs=arcs)
        else
            p0,pn=_step_plane_frame(entities,f.surf,caller)
            spans=false
            for ur in uses
                ur.kind in (:line,:polyline) || return nothing
                pts=ur.kind==:polyline ? ur.aux.pts : (ur.u.v1,ur.u.v2)
                for v in pts
                    abs(_brep_vdot(_brep_vsub(v,p0),pn))<=ptol || throw(
                        ArgumentError("$caller: edge #$(ur.u.ec.id) leaves "*
                            "wall PLANE #$(f.surf.id)"))
                end
                axial=all(
                    abs(_brep_vdot(_brep_vunit(_brep_vsub(pts[k+1],pts[k]),
                        caller),z))>=1.0-atol for k in 1:length(pts)-1)
                if axial
                    t1=_brep_vdot(first(pts),z);t2=_brep_vdot(last(pts),z)
                    (abs(t1-h0)<=htol && abs(t2-h1)<=htol ||
                     abs(t1-h1)<=htol && abs(t2-h0)<=htol) || return nothing
                    spans=true
                else
                    ths=[_brep_vdot(v,z) for v in pts]
                    all(t->abs(t-ths[1])<=htol,ths) || return nothing
                    (abs(ths[1]-h0)<=htol || abs(ths[1]-h1)<=htol) ||
                        return nothing
                end
            end
            spans || return nothing
        end
    end

    # ── circular-edge groups + coverage audit ─────────────────────────────
    # One geometric circle per (center, radius, axis-line): every use's
    # canonical coverage window is anchored by vertex position in the group
    # reference frame — the first incident cylinder wall's frame, else the
    # first circle's own — and each elementary arc interval must be covered
    # exactly twice.
    groups=NamedTuple[]
    useg=Dict{Tuple{Int,Int,Int},NamedTuple}()   # (fi,bi,ui) -> group record
    for (i,bound) in enumerate(fuse),(bi,b) in enumerate(bound),
            (ui,ur) in enumerate(b)
        ur.kind in (:arc,:circle) || continue
        a=ur.aux;tol=1e-9*max(a.r,1.0)
        gi=0
        for (k,g) in enumerate(groups)
            _posnear(a.center,g.center,tol) && abs(a.r-g.r)<=tol &&
                abs(_brep_vdot(a.nr,g.nr))>=1.0-atol && (gi=k;break)
        end
        if gi==0
            push!(groups,(center=a.center,r=a.r,nr=a.nr,x=a.x,y=a.y,
                          wall=0,uses=Tuple{Int,Int,Int,NamedTuple}[],
                          knots=Float64[],
                          pieces=Tuple{Float64,Float64}[]))
            gi=length(groups)
        end
        push!(groups[gi].uses,(i,bi,ui,ur))
        useg[(i,bi,ui)]=(gi=gi,)
    end
    # Group reference frames from incident cylinder walls (the first wall
    # wins; split walls on one cylinder share its frame).
    for (i,w) in wallc,a in w.arcs
        gidx=useg[(i,1,a.ui)].gi
        g=groups[gidx]
        g.wall==0 || continue
        groups[gidx]=(;g...,nr=w.z,x=w.x,y=w.y,wall=i)
    end
    # Coverage windows and per-group knot partitions.
    for (gi,g) in enumerate(groups)
        knots=Float64[];entries=NamedTuple[]
        nclosed=0
        for uc in g.uses
            ur=uc[4]
            if ur.kind==:circle
                nclosed+=1
                push!(entries,(fi=uc[1],bi=uc[2],ui=uc[3],ur=ur,
                               lo=0.0,hi=2π,startp=ur.u.v1))
                continue
            end
            a=ur.aux
            sgn=_brep_vdot(a.nr,g.nr)>0
            pos=sgn ? a.start : a.stop
            d=_brep_vsub(pos,g.center)
            lo=mod(atan(_brep_vdot(d,g.y),_brep_vdot(d,g.x)),2π)
            hi=lo+a.span
            push!(entries,(fi=uc[1],bi=uc[2],ui=uc[3],ur=ur,
                           lo=lo,hi=hi,startp=pos))
            push!(knots,lo,mod(hi,2π))
        end
        sort!(knots);bk=Float64[]
        for k in knots
            (isempty(bk) || k-bk[end]>atol) && push!(bk,k)
        end
        (length(bk)>1 && bk[1]+2π-bk[end]<=atol) && pop!(bk)
        pieces=Tuple{Float64,Float64}[]
        for j in eachindex(bk)
            push!(pieces,(bk[j],j<length(bk) ? bk[j+1] : bk[1]+2π))
        end
        groups[gi]=(;g...,knots=bk,pieces=pieces)
        if isempty(pieces)
            # Under/over-shared closed rims are not this certificate's
            # family; defer to the next certificate's diagnostic.
            nclosed==2 || return nothing
            continue
        end
        # Coverage count per elementary interval; remember each use's
        # covered pieces in cyclic order starting at its window.
        cov=zeros(Int,length(pieces))
        covered=Dict{Int,Vector{Int}}()
        for (k,e) in enumerate(entries)
            e.ur.kind==:circle && (cov.+=1;continue)
            lo,hi=e.lo,e.hi;list=Int[]
            for (j,(a,b)) in enumerate(pieces)
                hit=(a>=lo-atol && b<=hi+atol) ||
                    (a+2π>=lo-atol && b+2π<=hi+atol) ||
                    (a-2π>=lo-atol && b-2π<=hi+atol)
                hit && (cov[j]+=1;push!(list,j))
            end
            covered[k]=list
        end
        for (j,c) in enumerate(cov)
            # Intervals must be uncovered (a partial wall's complement) or
            # welded by exactly two uses; an odd count means the shell is
            # not a closed prismatic weld — defer to the next certificate.
            (c==0 || c==2) || return nothing
        end
        for (k,e) in enumerate(entries)
            key=(e.fi,e.bi,e.ui)
            if e.ur.kind==:circle
                useg[key]=(gi=gi,circle=true)
            else
                list=covered[k]
                isempty(list) && return nothing
                j0=findfirst(list) do j
                    (a,b)=pieces[j]
                    abs(a-e.lo)<=atol || abs(a-e.lo-2π)<=atol ||
                        abs(a+2π-e.lo)<=atol
                end
                j0===nothing && return nothing
                rot=vcat(list[j0:end],list[1:j0-1])
                useg[key]=(gi=gi,circle=false,lo=e.lo,hi=e.hi,
                           startp=e.startp,pieces=rot)
            end
        end
    end

    # ── per-wall rim congruence + strip windows ───────────────────────────
    # knots_w: shared partition of the two rim circles in the common owner
    # frame — the union of the two rim groups' knots. Both groups share the
    # first incident wall's frame (split walls on one cylinder use the same
    # frame), so the angular sets are directly comparable.
    for (i,w) in wallc
        rims=Dict{Int,Int}()   # gi -> rim
        for a in w.arcs
            rims[useg[(i,1,a.ui)].gi]=a.rim
        end
        length(rims)==2 || return nothing
        g0i=first(k for (k,v) in rims if v==1)
        g1i=first(k for (k,v) in rims if v==2)
        g0=groups[g0i];g1=groups[g1i]
        kw=sort(union(g0.knots,g1.knots))
        isempty(kw) && (kw=[0.0])
        pieces_w=Tuple{Float64,Float64}[]
        for j in eachindex(kw)
            push!(pieces_w,(kw[j],j<length(kw) ? kw[j+1] : kw[1]+2π))
        end
        covs=[Tuple{Float64,Float64}[],Tuple{Float64,Float64}[]]
        for a in w.arcs
            if a.ur.kind==:circle
                push!(covs[a.rim],(0.0,2π))
            else
                rec=useg[(i,1,a.ui)]
                push!(covs[a.rim],(rec.lo,rec.hi))
            end
        end
        covered_pieces=Tuple{Float64,Float64}[]
        for (a,b) in pieces_w
            inr=[any(a>=lo-atol && b<=hi+atol ||
                     a+2π>=lo-atol && b+2π<=hi+atol for (lo,hi) in covs[r])
                 for r in 1:2]
            inr[1]!=inr[2] && return nothing
            inr[1] && push!(covered_pieces,(a,b))
        end
        wallc[i]=(;w...,pieces=covered_pieces,g0i=g0i,g1i=g1i)
    end

    # ── emission: canonical chains weld shared boundaries bitwise ─────────
    pool=NTuple{3,Float64}[]
    chainmap=Dict{Tuple{Int,Float64,Float64},UnitRange{Int}}()
    function arc_chain!(gi,a,b)
        key=(gi,a,b)
        rng=get(chainmap,key,nothing)
        rng===nothing || return rng
        g=groups[gi];w=b-a
        m=max(1,ceil(Int,_BREP_PRISM_ARC_SEGMENTS*w/(2π)))
        r0=length(pool)+1
        for k in 0:m
            th=a+w*k/m;c=cos(th);s=sin(th)
            push!(pool,(g.center[1]+g.r*(c*g.x[1]+s*g.y[1]),
                        g.center[2]+g.r*(c*g.x[2]+s*g.y[2]),
                        g.center[3]+g.r*(c*g.x[3]+s*g.y[3])))
        end
        rng=r0:length(pool);chainmap[key]=rng
        return rng
    end
    items=Tuple{Symbol,UnitRange{Int},Int}[]
    emit!(kind,pts,fi)=(push!(items,(kind,
        length(pool)+1:length(pool)+length(pts),fi));append!(pool,pts);nothing)

    # Traversal-ordered ring assembly for one wire of face i. A lone closed
    # circle contributes the whole cyclic chain; every other use contributes
    # its chain points in traversal order skipping the first (the previous
    # use's junction node).
    function wire_ring!(bound,i,bi)
        ring=NTuple{3,Float64}[]
        if length(bound)==1 && bound[1].kind==:circle
            rec=useg[(i,bi,1)];g=groups[rec.gi]
            if isempty(g.pieces)
                rng=arc_chain!(rec.gi,0.0,2π)
                append!(ring,pool[rng][1:end-1])
            else
                for (a,b) in g.pieces
                    rng=arc_chain!(rec.gi,a,b)
                    append!(ring,pool[rng][1:end-1])
                end
            end
            return ring
        end
        for (ui,ur) in enumerate(bound)
            if ur.kind==:line
                push!(ring,ur.u.v2)
            elseif ur.kind==:polyline
                ch=ur.aux.pts
                fwd=_posnear(ur.u.v1,first(ch),ptol)
                append!(ring,fwd ? ch[2:end] : reverse(ch[1:end-1]))
            elseif ur.kind==:arc
                rec=useg[(i,bi,ui)]
                fwd=_posnear(ur.u.v1,rec.startp,ptol)
                ps=rec.pieces
                if fwd
                    for j in ps
                        rng=arc_chain!(rec.gi,
                            groups[rec.gi].pieces[j][1],
                            groups[rec.gi].pieces[j][2])
                        append!(ring,pool[rng][2:end])
                    end
                else
                    for j in reverse(ps)
                        rng=arc_chain!(rec.gi,
                            groups[rec.gi].pieces[j][1],
                            groups[rec.gi].pieces[j][2])
                        append!(ring,reverse(pool[rng][1:end-1]))
                    end
                end
            else
                throw(ArgumentError(
                    "$caller: closed CIRCLE edge #$(ur.u.ec.id) shares a "*
                    "wire with other edges — not a prismatic cap"))
            end
        end
        return ring
    end

    for i in cap_i
        rings=[wire_ring!(bound,i,bi) for (bi,bound) in enumerate(fuse[i])]
        any(r->length(r)<3,rings) && throw(ArgumentError(
            "$caller: cap face #$(faces[i].id) has a degenerate wire"))
        if length(rings)>1
            nrm=_newell_normal(rings[argmax(length.(rings))])
            drop=argmax(abs.(nrm));axes=filter(!=(drop),collect(1:3))
            projected=[[(p[axes[1]],p[axes[2]]) for p in ring]
                       for ring in rings]
            outer=0
            for k in eachindex(projected)
                inside=false
                for l in eachindex(projected)
                    k==l && continue
                    if _brep_point_in_wire(first(projected[k]),projected[l])
                        inside=true;break
                    end
                end
                if !inside
                    outer==0 || throw(ArgumentError(
                        "$caller: face #$(faces[i].id) has more than one "*
                        "outer boundary wire"))
                    outer=k
                end
            end
            outer==0 && throw(ArgumentError(
                "$caller: face #$(faces[i].id) has no containing outer wire"))
            rings=vcat(rings[outer:outer],
                       rings[setdiff(eachindex(rings),outer)])
            a1=_signed_area2([p[axes[1]] for p in rings[1]],
                             [p[axes[2]] for p in rings[1]])
            for k in 2:length(rings)
                ak=_signed_area2([p[axes[1]] for p in rings[k]],
                                 [p[axes[2]] for p in rings[k]])
                signbit(ak)==signbit(a1) && reverse!(rings[k])
            end
            for tri in _brep_planar_face_holes(rings,caller)
                emit!(:tri,tri,i)
            end
        else
            emit!(:ring,only(rings),i)
        end
    end
    for i in wall_i
        haskey(wallc,i) && continue
        ring=wire_ring!(only(fuse[i]),i,1)
        length(ring)>=3 || throw(ArgumentError(
            "$caller: wall face #$(faces[i].id) is degenerate"))
        emit!(:ring,ring,i)
    end
    for (i,w) in wallc
        for (a,b) in w.pieces
            r0=arc_chain!(w.g0i,a,b);r1=arc_chain!(w.g1i,a,b)
            p0v=pool[r0];p1v=pool[r1]
            length(p0v)==length(p1v) || throw(ArgumentError(
                "$caller: rim chains on wall face #$(faces[i].id) are "*
                "misaligned"))
            for k in 1:length(p0v)-1
                emit!(:quad,(p0v[k],p0v[k+1],p1v[k+1],p1v[k]),i)
            end
        end
    end
    nodes,canon=_dedup_positions(pool,caller)
    itris=NTuple{3,Int32}[];tface=Int[]
    for (kind,rng,fi) in items
        if kind==:ring
            ir=Int32[];prev=Int32(0)
            for j in rng
                v=canon[j]
                v!=prev && push!(ir,v);prev=v
            end
            length(ir)>1 && first(ir)==last(ir) && pop!(ir)
            n0=length(itris)
            _ear_clip!(itris,ir,nodes,caller)
            append!(tface,fill(fi,length(itris)-n0))
        elseif kind==:quad
            a,b,c,d=(Int32(canon[j]) for j in rng)
            push!(itris,(a,b,c));push!(itris,(a,c,d))
            append!(tface,(fi,fi))
        else
            a,b,c=(Int32(canon[j]) for j in rng)
            push!(itris,(a,b,c));push!(tface,fi)
        end
    end
    _checked_brep_mul(3,length(itris),caller,"prismatic shell triangle count")
    length(itris)<=_BREP_MAX_SURFACE_TRIS || throw(ArgumentError(
        "$caller: prismatic shell exceeds the surface triangle limit"))
    tris,edge_use=_step_shell_flip_parity(itris,tface,length(faces),
        "prismatic shell",caller)
    return _step_shell_audit(tris,nodes,tface,length(faces),edge_use,
        "prismatic shell",caller)
end
