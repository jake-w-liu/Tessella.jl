# `occ.importShapes` / `OCC_Internals::importShapes` — ISO-10303-21 STEP
# solids imported as native OCC records (points/lines/wires/plane
# faces/shells/volumes) rather than the triangulated `import_step`
# meshes. The BRep module supplies the parser and the face/bound/edge
# graph walk; this file converts each closed shell's rings into shared
# vertex/edge entities — identical vertices and edges publish once, so
# the shell validator's two-surface-per-curve rule holds by
# construction. Bounded to planar faces (STEP `PLANE` or untyped
# `FACE`); curved faces reject loudly until curved STEP→record
# conversion lands.

# Solid-level dedupe state: vertex positions → point tags, unordered
# vertex-tag pairs → (line id, stored stop tag).
struct _OCCImportDedupe
    ppos::Vector{NTuple{3,Float64}}
    ptag::Vector{Int}
    edge::Dict{NTuple{2,Int},Tuple{UInt64,Int}}
end
_OCCImportDedupe()=_OCCImportDedupe(NTuple{3,Float64}[],Int[],
                                   Dict{NTuple{2,Int},Tuple{UInt64,Int}}())

function _occ_import_point!(m,d:: _OCCImportDedupe,p)
    for i in eachindex(d.ppos)
        _model_occ_vertices_coincide(d.ppos[i],p) && return d.ptag[i]
    end
    t=_model_occ_add_point!(m,p[1],p[2],p[3])
    push!(d.ppos,p); push!(d.ptag,t)
    return t
end

# One ring (ordered positions, no duplicated closing point) → wire tag.
# `flip` winds the loop backwards (hole wires are consumed CW relative
# to the face normal by `add_plane_surface` convention checks — STEP
# bounds already carry the correct orientation so `flip` stays false).
function _occ_import_wire!(m,d,kernel,ring,flip,caller)
    verts=[_occ_import_point!(m,d,p) for p in ring]
    uses=Tuple{UInt64,Bool}[]
    K=length(verts)
    for j in 1:K
        a,b=verts[j],verts[mod1(j+1,K)]
        key=a<b ? (a,b) : (b,a)
        eu=get(d.edge,key,nothing)
        if eu===nothing
            et=_model_occ_add_line!(m,a,b)
            eu=(kernel.bindings[(1,et)],b)
            d.edge[key]=eu
        end
        rev=eu[2]!=b
        push!(uses,flip ? (eu[1],!rev) : (eu[1],rev))
    end
    flip && (uses=reverse(uses))
    return _occ_wire_from_uses!(m,uses,caller)
end

# One STEP face: rings from `_step_face_parts` (outer first). Planar
# surfaces only — a non-plane surface record rejects.
function _occ_import_step_face(m,d,kernel,entities,fe,curvekinds,caller)
    rings,surf,sense,straight=BRep._step_face_parts(entities,fe,
                                                  curvekinds,caller)
    rings===nothing && throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) uses an unsupported loop kind"))
    surf===nothing || surf.kind=="PLANE" || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) uses $(surf.kind) — only planar " *
        "faces import currently"))
    straight || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has non-line edges — only " *
        "straight-edge faces import currently"))
    wires=[_occ_import_wire!(m,d,kernel,r,false,caller) for r in rings]
    return _model_occ_add_plane_surface!(m,wires)
end

# One solid root → [(3, tag)]. Each CLOSED_SHELL becomes a surface
# loop; `BREP_WITH_VOIDS`-style later shells carry role −1 (cavities).
function _occ_import_step_solid(m,brep,entities,caller)
    kernel=m.occ_kernel
    length(brep.args)>=2 || throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shell"))
    d=_OCCImportDedupe(); curvekinds=Set{String}()
    shell_tags=Int[]
    for (si,sarg) in enumerate(brep.args[2:end])
        items=sarg isa Vector ? sarg : Any[sarg]
        for item in items
            se=BRep._deref_entity(entities,item)
            se===nothing && throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references a missing shell"))
            se.kind=="CLOSED_SHELL" || throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references $(se.kind) " *
                "#$(se.id) where a CLOSED_SHELL was expected"))
            i=BRep._skip_name(se.args)
            flist=se.args[i]
            flist isa Vector || throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(se.id) face list is malformed"))
            face_tags=Int[]
            for fref in flist
                fe=BRep._deref_entity(entities,fref)
                (fe!==nothing &&
                 fe.kind in BRep._STEP_BREP_FACE_KINDS) || throw(
                    ArgumentError(
                        "$caller: CLOSED_SHELL #$(se.id) member is not a face"))
                push!(face_tags,
                      _occ_import_step_face(m,d,kernel,entities,fe,
                                            curvekinds,caller))
            end
            isempty(face_tags) && throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(se.id) is empty"))
            push!(shell_tags,
                  _model_occ_add_surface_loop!(m,face_tags))
        end
    end
    isempty(shell_tags) && throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shells"))
    return [(3,_model_occ_add_volume!(m,shell_tags))]
end

function _model_occ_import_shapes(m::GeoModel,path;highest_dim_only=true,
                                  caller="OCC.importShapes")
    isfile(path) || throw(ArgumentError("$caller: missing file $path"))
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "$caller: $path is not an ISO-10303-21 STEP file"))
    entities=BRep.parse_step_entities(source)
    breps=sort!([e for e in values(entities)
                 if e.kind in BRep._STEP_BREP_SOLID_KINDS]; by=e->e.id)
    isempty(breps) && throw(ArgumentError(
        "$caller: no MANIFOLD_SOLID_BREP/FACETED_BREP/BREP_WITH_VOIDS " *
        "found in $path"))
    out=Tuple{Int,Int}[]
    for brep in breps
        append!(out,_occ_import_step_solid(m,brep,entities,caller))
    end
    return out
end

# ── `healShapes` ──────────────────────────────────────────────────────────────
#
# `OCC_Internals::healShapes`: weld coincident-but-distinct vertices and
# edges in the target entities' topology within `tolerance`
# (`Precision::Confusion` default), then — when `sewFaces` — regroup
# healed faces into shells and bind solids. Healing works on the stored
# graph: vertex entities cluster by position, edge children rewire to
# the cluster's canonical vertex, and line edges sharing both endpoints
# collapse to one curve entity. Faces/curves are not re-fit — tolerance
# healing moves topology, not geometry.
# Reverse a pcurve onto an oppositely-directed replacement edge:
# `pc'(s)=pc(T−s)` with `T=t0+t1` over the shared `[0,L]` line range —
# lin2d via affine reparam, spl2d by knot-mirror + pole reversal.
function _occ_heal_pcurve_reverse(pc,T::Float64,kernel,caller)
    if hasproperty(pc,:neg2d)
        inner=pc.neg2d[1]
        hasproperty(inner,:lin2d) || throw(ArgumentError(caller*
            ": cannot weld edge — unsupported negated pcurve kind"))
        o,d=inner.lin2d
        return (neg2d=((lin2d=((o[1]-T*d[1],o[2]-T*d[2]),d),),))
    end
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=((o[1]+T*d[1],o[2]+T*d[2]),(-d[1],-d[2])),)
    elseif hasproperty(pc,:spl2d)
        knots,p,cps=pc.spl2d
        return (spl2d=(T .- reverse(knots),p,reverse(cps)),)
    end
    throw(ArgumentError(caller*
        ": cannot weld edge — unsupported pcurve kind for reversal"))
end

function _model_occ_heal_shapes(m::GeoModel,dim_tags;tolerance=-1.0,
                                fixDegenerated=true,fixSmallEdges=true,
                                fixSmallFaces=true,sewFaces=true,
                                makeSolids=true,
                                caller="OCC.healShapes")
    kernel=m.occ_kernel
    tol=tolerance>0.0 ? tolerance : _OCC_CONFUSION
    roots=Tuple{Int,Int}[]
    if isempty(dim_tags)
        for (d,t) in keys(kernel.bindings)
            d>=0 && push!(roots,(d,t))
        end
    else
        for (d,t) in dim_tags
            d=Int(d); t=abs(Int(t))
            haskey(kernel.bindings,(d,t)) || throw(ArgumentError(
                "$caller: unknown OCC entity ($d,$t)"))
            push!(roots,(d,t))
        end
    end
    # Reachable-vertex weld: gather every vertex id under the roots
    # (shapes' children), cluster positions within `tol`, rewire edges.
    vids=UInt64[]
    stack=UInt64[kernel.bindings[r] for r in roots]
    seen=Set{UInt64}()
    while !isempty(stack)
        id=pop!(stack)
        id in seen && continue; push!(seen,id)
        sh=kernel.shapes[id]
        if sh.dimension==0
            push!(vids,id); continue
        end
        append!(stack,sh.children)
        hasproperty(sh.geometry,:uses) || continue
        for (uid,_) in sh.geometry.uses
            uid in seen || push!(stack,uid)
        end
    end
    # Cluster by position — canonical = first member.
    canon=Dict{UInt64,UInt64}()
    reps=UInt64[]
    for id in vids
        p=kernel.shapes[id].geometry.point
        hit=0
        for r in reps
            _occ_modulus(_occ3_sub(kernel.shapes[r].geometry.point,p))<=tol &&
                (hit=r; break)
        end
        if hit==0
            push!(reps,id); canon[id]=id
        else
            canon[id]=hit
        end
    end
    # Rewire edge children to canonical vertices.
    for id in seen
        sh=kernel.shapes[id]
        sh.dimension==1 || continue
        ch=sh.children
        nc=UInt64[get(canon,c,c) for c in ch]
        nc==ch && continue
        kernel.shapes[id]=_ModelOCCShape(1,sh.geometry,nc)
    end
    # Edge dedupe: two lines sharing both endpoints collapse — faces
    # then share the curve entity (needed for shell validation).
    ecanon=Dict{UInt64,UInt64}()
    if fixSmallEdges || sewFaces
        sig=Dict{NTuple{2,UInt64},UInt64}()
        for id in seen
            sh=kernel.shapes[id]
            sh.dimension==1 || continue
            get(sh.geometry.geometry,:occ,nothing)==:line || continue
            a,b=sh.children; key=a<b ? (a,b) : (b,a)
            prev=get(sig,key,0)
            if prev==0
                sig[key]=id
            else
                ecanon[id]=prev
            end
        end
        if !isempty(ecanon)
            # Re-key face pcurves while wires still reference the
            # face's own edges: each pcurve's key moves from the old
            # edge tag to the canonical edge's tag. When the canonical
            # edge runs opposite the dropped one, `fwd` — parametrized
            # along the edge's own direction — must be reversed.
            for id in seen
                sh=kernel.shapes[id]
                sh.dimension==2 || continue
                fr=sh.geometry.geometry
                hasproperty(fr,:pcurves) || continue
                pc=Dict{Int,NamedTuple}()
                for cid in sh.children
                    wsh=kernel.shapes[cid]
                    hasproperty(wsh.geometry,:uses) || continue
                    for (u,_) in wsh.geometry.uses
                        ot=kernel.shape_tags[u][2]
                        nu=get(ecanon,u,u)
                        nt=kernel.shape_tags[nu][2]
                        ent=fr.pcurves[ot]
                        if nu!=u && kernel.shapes[u].children[1]!=
                                kernel.shapes[nu].children[1]
                            eg=kernel.shapes[u].geometry.geometry
                            ent=(fwd=_occ_heal_pcurve_reverse(ent.fwd,
                                     eg.t0+eg.t1,kernel,caller),
                                 rev=nothing)
                        end
                        pc[nt]=ent
                    end
                end
                kernel.shapes[id]=_ModelOCCShape(2,
                    merge(sh.geometry,(geometry=merge(fr,(pcurves=pc,)),)),
                    sh.children)
            end
            # Rewire wire/shell `uses` to canonical edges — flipping
            # `rev` when the dropped edge's stored direction opposes
            # the canonical edge's.
            for id in seen
                sh=kernel.shapes[id]
                hasproperty(sh.geometry,:uses) || continue
                uses=Tuple{UInt64,Bool}[]
                for (u,r) in sh.geometry.uses
                    nu=get(ecanon,u,u)
                    if nu==u
                        push!(uses,(u,r))
                    else
                        och=kernel.shapes[u].children
                        nch=kernel.shapes[nu].children
                        push!(uses,(nu,och[1]==nch[1] ? r : !r))
                    end
                end
                uses==sh.geometry.uses && continue
                kernel.shapes[id]=_ModelOCCShape(sh.dimension,
                    merge(sh.geometry,(uses=uses,)),
                    UInt64[u for (u,_) in uses])
            end
        end
    end
    # Unbind welded-away entities: duplicates no longer exist as
    # distinct geometry once every reference points at a canonical.
    for (id,rep) in canon
        id==rep && continue
        entry=get(kernel.shape_tags,id,nothing)
        entry===nothing && continue
        delete!(kernel.bindings,entry);delete!(kernel.shape_tags,id)
        push!(kernel.removed,entry);kernel.changed=true
    end
    for (id,rep) in ecanon
        id==rep && continue
        entry=get(kernel.shape_tags,id,nothing)
        entry===nothing && continue
        delete!(kernel.bindings,entry);delete!(kernel.shape_tags,id)
        push!(kernel.removed,entry);kernel.changed=true
    end
    return nothing
end
