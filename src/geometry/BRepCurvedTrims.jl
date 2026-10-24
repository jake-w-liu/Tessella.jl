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

    # ── orientation: flip-parity propagation over the face-adjacency graph ─
    # edge_tris: undirected edge -> the tris using it, with each tri's
    # traversal direction. Two incident faces orient consistently iff their
    # flip parities differ exactly when the two uses already disagree.
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
    adj=[Int[] for _ in faces]; diff=Dict{Tuple{Int,Int},Bool}()
    for (key,(t1,t2)) in edge_use
        t2==0 && throw(ArgumentError(
            "$caller: tubed shell is not closed — edge $key bounds a "*
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
                "$caller: tubed shell faces $f1 and $f2 cannot be "*
                "consistently oriented across their shared edges"))
        end
    end
    flip=falses(length(faces));seen=falses(length(faces))
    for seed in eachindex(faces)
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
                        "$caller: tubed shell is not orientable — "*
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
    # Directed shared-edge audit: every undirected edge once each way.
    audit=Dict{Tuple{Int32,Int32},Int}()
    for t in tris,(a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        audit[(a,b)]=get(audit,(a,b),0)+1
    end
    for (key,(t1,t2)) in edge_use
        ab=get(audit,(key[1],key[2]),0);ba=get(audit,(key[2],key[1]),0)
        (ab==1 && ba==1) || throw(ArgumentError(
            "$caller: tubed shell is inconsistently oriented across edge $key"))
    end
    scale=max(_poly_scale(nodes),1.0)
    V=sum(_brep_vdot(nodes[t[1]],
        _brep_vcross(nodes[t[2]],nodes[t[3]]))/6 for t in tris)
    abs(V)>1e-12*scale*scale*scale || throw(ArgumentError(
        "$caller: tubed shell bounds zero volume"))
    if V<0
        tris=map(t->(t[1],t[3],t[2]),tris)
    end
    C=Matrix{Float64}(undef,3,length(nodes))
    for (i,p) in enumerate(nodes); C[:,i]=[p...]; end
    T=Matrix{Int32}(undef,3,length(tris))
    for (i,t) in enumerate(tris); T[:,i]=[t...]; end
    surface=Mesh(C;tris=T)
    diag=validate(surface)
    diag.ok || throw(ErrorException(
        "$caller: tubed shell boundary mesh is invalid — "*
        join(diag.messages,"; ")))
    return surface
end
