# ============================================================================
# `occ.addThruSections` — `OCC_Internals::addThruSections`:
# `BRepOffsetAPI_ThruSections` lofting through a sequence of bound wire
# sections (dim −1 shapes). Wires are consumed verbatim — gmsh calls
# `CheckCompatibility(false)`, so edges pair strictly in traversal order
# (no origin/orientation normalization); mismatched member counts reject
# like `BRepFill`'s inconsistent-profile check, and `makeSolid` requires
# closed wires exactly like `ts.AddWire`'s `wire.Closed()` gate.
#
# Every member edge is canonicalized to a directed NURBS curve on the
# normalized internal parameter ũ ∈ [0,1] (knot rescale is affine, hence
# exact; traversal-reversed members are pole/knot-reversed). `makeRuled`
# builds the v-degree-1 ruled patch between each pair of consecutive
# unified nets — the `GeomFill_Generator` construction — producing
# (nsec−1)·ne lateral faces. The default smooth path interpolates all
# section nets in v (`GeomFill_BSplineCurves` skinning in homogeneous
# coordinates) at ChordLength/Centripetal/IsoParametric section
# parameters with degree min(3, nsec−1) capped by `maxDegree`, producing
# one face per edge column. Unification elevates lower-degree members and
# merges knot vectors through exact homogeneous knot insertion, so every
# lateral patch interpolates its boundary edges exactly.
#
# Seam edges join corresponding section vertices: straight line records on
# ruled faces, the surface's v-direction boundary interpolant (a
# `copied_curve` B-spline) on smooth faces. `makeSolid` caps the first/last
# sections through `_model_occ_plane_face_geometry` — the upstream
# `BRepBuilderAPI_MakeFace` contract, so nonplanar caps reject — and binds
# the whole shell as a volume. Section edges keep their own t→ũ pcurves:
# affine `lin2d` for line/B-spline members, the `tan2d` record for conics
# whose t→ũ map is the normalized tangent half-angle reparam. The loft
# never mutates the registry before every patch and record has staged.
# ============================================================================

# ── homogeneous-pole helpers ────────────────────────────────────────────────

function _loft_hcurve(c::NURBS.NURBSCurve)
    n=length(c.controls)
    Pw=Vector{NTuple{4,Float64}}(undef,n)
    @inbounds for i in 1:n
        w=c.weights[i]
        Pw[i]=(w*c.controls[i][1],w*c.controls[i][2],w*c.controls[i][3],w)
    end
    return Pw
end

function _loft_dehom(Pw,caller)
    n=length(Pw)
    controls=Vector{NTuple{3,Float64}}(undef,n)
    weights=Vector{Float64}(undef,n)
    @inbounds for i in 1:n
        w=Pw[i][4]
        (isfinite(w) && w>0.0) || throw(ArgumentError(
            "$caller: loft produced a non-positive control weight"))
        controls[i]=(Pw[i][1]/w,Pw[i][2]/w,Pw[i][3]/w)
        weights[i]=w
    end
    return controls,weights
end

# ── NURBS curve parameter surgery (all locus-exact) ─────────────────────────

# Same curve, reversed traversal direction, unchanged parameter range.
function _loft_nc_reverse(c::NURBS.NURBSCurve)
    s0=c.knots[1]; s1=c.knots[end]; n=length(c.controls)
    knots=Float64[s0+s1-c.knots[end+1-k] for k in 1:length(c.knots)]
    controls=NTuple{3,Float64}[c.controls[n+1-i] for i in 1:n]
    weights=Float64[c.weights[n+1-i] for i in 1:n]
    return NURBS.NURBSCurve(c.degree,knots,controls,weights)
end

# Affine rescale of the knot domain to [0,1] — exact reparametrization.
function _loft_nc_rescale01(c::NURBS.NURBSCurve)
    s0=c.knots[1]; s1=c.knots[end]
    Δ=s1-s0
    Δ>0.0 || throw(ArgumentError("addThruSections: degenerate knot range"))
    knots=Float64[(k-s0)/Δ for k in c.knots]
    return NURBS.NURBSCurve(c.degree,knots,c.controls,c.weights)
end

# Boehm single-knot insertion on the homogeneous pole net — the curve (as a
# parametrized locus) is unchanged; the knot vector gains `u` once.
function _loft_nc_insert(c::NURBS.NURBSCurve,u::Float64)
    p=c.degree; U=c.knots; n=length(c.controls)
    Pw=_loft_hcurve(c)
    # Span index: U[k] <= u < U[k+1] with k in [p+1, n].
    k=clamp(searchsortedlast(U,u),p+1,n)
    U2=Vector{Float64}(undef,length(U)+1)
    copyto!(U2,1,U,1,k); U2[k+1]=u; copyto!(U2,k+2,U,k+1,length(U)-k)
    P2=Vector{NTuple{4,Float64}}(undef,n+1)
    @inbounds for i in 1:k-p
        P2[i]=Pw[i]
    end
    @inbounds for i in k-p+1:k
        den=U[i+p]-U[i]
        α=den>0.0 ? (u-U[i])/den : 0.0
        a,b=Pw[i-1],Pw[i]
        P2[i]=(fma(α,b[1],fma(-α,a[1],a[1])),
               fma(α,b[2],fma(-α,a[2],a[2])),
               fma(α,b[3],fma(-α,a[3],a[3])),
               fma(α,b[4],fma(-α,a[4],a[4])))
    end
    @inbounds for i in k+1:n+1
        P2[i]=Pw[i-1]
    end
    controls,weights=_loft_dehom(P2,"addThruSections")
    return NURBS.NURBSCurve(p,U2,controls,weights)
end

# Raise one degree: drive interior knots to full multiplicity (Bézier split),
# elevate each segment (Prautzsch Q_i = (i/q)·P_{i−1} + (1−i/q)·P_i), rejoin
# at multiplicity q. Exact — the curve locus is unchanged.
function _loft_nc_elevate1(c::NURBS.NURBSCurve)
    p=c.degree
    # Drive every distinct interior knot to multiplicity p.
    interior=c.knots[p+2:length(c.knots)-p-1]
    for u in unique(interior)
        mult=count(==(u),c.knots)
        for _ in mult+1:p
            c=_loft_nc_insert(c,u)
        end
    end
    U=c.knots; Pw=_loft_hcurve(c)
    # Full-multiplicity interior knots split the pole net into Bézier
    # segments of p+1 poles each.
    nseg=(length(Pw)-1)÷p
    nseg*p+1==length(Pw) || throw(ArgumentError(
        "addThruSections: inconsistent Bézier decomposition"))
    q=p+1
    H=Vector{NTuple{4,Float64}}(undef,nseg*q+1)
    @inbounds for s in 0:nseg-1
        base=s*p; dst=s*q
        H[dst+1]=Pw[base+1]
        for i in 1:p
            α=i/q
            a=Pw[base+i]; b=Pw[base+i+1]
            H[dst+i+1]=(fma(α,b[1],fma(-α,a[1],a[1])),
                        fma(α,b[2],fma(-α,a[2],a[2])),
                        fma(α,b[3],fma(-α,a[3],a[3])),
                        fma(α,b[4],fma(-α,a[4],a[4])))
        end
        s==nseg-1 && (H[dst+q+1]=Pw[base+p+1])
    end
    # Knots: ends at multiplicity q+1, each distinct interior knot at q.
    U2=Vector{Float64}(undef,2*(q+1)+(nseg-1)*q)
    U2[1:q+1].=U[1]; U2[end-q:end].=U[end]
    o=q+1
    for u in unique(@view U[q+2:end-q-1])
        U2[o+1:o+q].=u; o+=q
    end
    controls,weights=_loft_dehom(H,"addThruSections")
    return NURBS.NURBSCurve(q,U2,controls,weights)
end

function _loft_nc_elevate_to(c::NURBS.NURBSCurve,degree)
    while c.degree<degree
        c=_loft_nc_elevate1(c)
    end
    return c
end

# Elevate all curves to the common maximum degree, then insert knots until
# every curve carries the merged knot vector (max multiplicity per distinct
# value). All surgery is locus-preserving, so the returned curves coincide
# pointwise with the inputs on the shared [0,1] domain.
function _loft_nc_unify(curves)
    pmax=maximum(c.degree for c in curves)
    cs=[c.degree<pmax ? _loft_nc_elevate_to(c,pmax) : c for c in curves]
    mults=Dict{Float64,Int}()
    for c in cs
        seen=Dict{Float64,Int}()
        for u in @view c.knots[pmax+2:length(c.knots)-pmax-1]
            seen[u]=get(seen,u,0)+1
        end
        for (u,m) in seen
            mults[u]=max(get(mults,u,0),m)
        end
    end
    out=Vector{NURBS.NURBSCurve}(undef,length(cs))
    for (i,c) in pairs(cs)
        for (u,m) in mults
            have=count(==(u),c.knots)
            for _ in have+1:m
                c=_loft_nc_insert(c,u)
            end
        end
        out[i]=c
    end
    return out
end

# ── section edge canonicalization ───────────────────────────────────────────
#
# `_loft_edge_nc(geometry, rev, caller)` → (nc, umap) with `nc` the directed
# NURBS curve on ũ∈[0,1] and `umap` describing the edge's own t→ũ map:
#   (:lin,o,d)          — ũ = o + d·t
#   (:tan,inner,u0,us)  — ũ = u0 + us·inner(t),  inner a _ModelOCCConicReparam

function _loft_edge_nc(geometry,rev,caller)
    if geometry.occ===:line
        cv=_model_occ_nurbs_line(geometry,caller)
        nc=cv.curve; inner=nothing; tan=false
    elseif geometry.occ===:circle || geometry.occ===:ellipse
        cv=_model_occ_nurbs_conic(geometry,caller)
        nc=cv.curve; inner=cv.reparam; tan=true
    elseif hasproperty(geometry,:curve)
        cv=geometry.curve
        if cv isa _ModelOCCReparamCurve
            nc=cv.curve; inner=cv.reparam; tan=true
        elseif cv isa _ModelOCCNURBSCurve
            nc=cv.curve; inner=nothing; tan=false
        else
            throw(ArgumentError(
                "$caller: unsupported edge curve kind $(geometry.occ)"))
        end
    else
        throw(ArgumentError(
            "$caller: unsupported edge curve kind $(geometry.occ)"))
    end
    s0,s1=nc.knots[1],nc.knots[end]
    inv=1.0/(s1-s0)
    if rev
        nc=_loft_nc_reverse(nc)
        # traversal position of edge-t: ũ = (s1 − s(t))/Δs
        umap=tan ? (:tan,inner,s1*inv,-inv) : (:lin,s1*inv,-inv)
    else
        umap=tan ? (:tan,inner,-s0*inv,inv) : (:lin,-s0*inv,inv)
    end
    return _loft_nc_rescale01(nc),umap
end

# ── v-direction interpolation (GeomFill_BSplineCurves skinning) ─────────────

# Section parameters per upstream ParType: IsoParametric → uniform [0,1];
# ChordLength (the default) → normalized cumulative chord lengths between
# section barycenters; Centripetal → sqrt-chord normalization.
function _loft_vparams(centers,parametrization,caller)
    n=length(centers)
    if parametrization=="" || parametrization=="ChordLength"
        power=1.0
    elseif parametrization=="Centripetal"
        power=0.5
    elseif parametrization=="IsoParametric"
        return Float64[(i-1)/(n-1) for i in 1:n]
    else
        throw(ArgumentError(
            "$caller: unknown parametrization '$parametrization'"))
    end
    v=Float64[0.0]
    for j in 2:n
        d2=0.0
        for k in 1:3
            d2+=(centers[j][k]-centers[j-1][k])^2
        end
        push!(v,last(v)+sqrt(d2)^power)
    end
    total=last(v)
    total>0.0 || throw(ArgumentError(
        "$caller: through sections have coincident centers"))
    return v./total
end

# Interpolation knots for q-degree clamped interpolation at parameters vbar
# (averaging method — the GeomFill/BSplCLib convention).
function _loft_vknots(vbar,q)
    n=length(vbar)
    knots=Vector{Float64}(undef,2*(q+1)+max(n-1-q,0))
    knots[1:q+1].=vbar[1]
    knots[end-q:end].=vbar[end]
    o=q+1
    for i in 2:n-q
        s=0.0
        for j in i:i+q-1
            s+=vbar[j]
        end
        knots[o+1]=s/q; o+=1
    end
    return knots
end

# All degree-q basis values at v ∈ [knots[s], knots[s+1]) — The NURBS Book
# A2.2 BasisFuns — scattered into an n-length row (active indices s−q..s).
function _loft_basis_row(knots,q,s,v,n)
    left=Vector{Float64}(undef,q+1)
    right=Vector{Float64}(undef,q+1)
    N=zeros(Float64,q+1); N[1]=1.0
    @inbounds for j in 1:q
        left[j]=v-knots[s-j+1]
        right[j]=knots[s+j]-v
        saved=0.0
        for r in 0:j-1
            den=right[r+1]+left[j-r]
            tmp=den>0.0 ? N[r+1]/den : 0.0
            N[r+1]=saved+right[r+1]*tmp
            saved=left[j-r]*tmp
        end
        N[j+1]=saved
    end
    row=zeros(Float64,n)
    for r in 0:q
        idx=s-q+r
        1<=idx<=n && (row[idx]=N[r+1])
    end
    return row
end

# Solve the n×n clamped interpolation A·X = R in homogeneous space for every
# u-pole column; `columns[j]` is section j's homogeneous net (length nu).
# Returns (rows, knots) — rows[i][j] is the interpolated homogeneous pole at
# u-index i, section row j.
function _loft_vinterp(columns,vbar,q,caller)
    n=length(vbar)
    nu=length(columns[1])
    knots=_loft_vknots(vbar,q)
    nb=length(knots)-q-1
    A=Matrix{Float64}(undef,n,n)
    for j in 1:n
        s=clamp(searchsortedlast(knots,vbar[j]),q+1,nb)
        A[j,:]=_loft_basis_row(knots,q,s,vbar[j],n)
    end
    F=try
        LinearAlgebra.lu(A)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        throw(ArgumentError(
            "$caller: through-sections interpolation is singular"))
    end
    rows=Vector{Vector{NTuple{4,Float64}}}(undef,nu)
    R=Matrix{Float64}(undef,n,4)
    for i in 1:nu
        for j in 1:n
            h=columns[j][i]
            R[j,1]=h[1]; R[j,2]=h[2]; R[j,3]=h[3]; R[j,4]=h[4]
        end
        X=F\R
        all(isfinite,X) || throw(ArgumentError(
            "$caller: through-sections interpolation is singular"))
        rows[i]=NTuple{4,Float64}[
            (X[j,1],X[j,2],X[j,3],X[j,4]) for j in 1:n]
    end
    return rows,knots
end

# ── edge/face record builders ───────────────────────────────────────────────

# `tan2d` pcurve: u(t) = u0 + us·inner(t) where `inner` is the edge's conic
# tangent half-angle reparam; v is constant `vv`.
_loft_tan2d(inner,u0,us,vv)=(tan2d=(inner,u0,us,vv),)
_loft_lin2d(o,d)=(lin2d=(o,d),)

function _loft_umap_pcurve(umap,vv)
    if umap[1]===:lin
        return _loft_lin2d((umap[2],vv),(umap[3],0.0))
    end
    return _loft_tan2d(umap[2],umap[3],umap[4],vv)
end

# Seam pcurve: edge-t → (u_side, v_lo + t·Δv/L) — `t` is the seam edge's own
# parameter over [t0,t1]; the face's v runs [v_lo, v_lo+Δv].
function _loft_seam_pcurve(u_side,t0,t1,v_lo,Δv)
    scale=Δv/(t1-t0)
    return _loft_lin2d((u_side,v_lo-t0*scale),(0.0,scale))
end

# ── driver ──────────────────────────────────────────────────────────────────
#
# `_model_occ_thru_sections!(m, tag, wire_tags; make_solid, make_ruled,
#   max_degree, continuity, parametrization, smoothing)` → outDimTags.
# Publication is staged: every shape rides one pending vector and binds only
# after all geometry has been constructed, so a failure leaves the kernel
# untouched.

function _model_occ_thru_sections!(m::GeoModel,tag,wire_tags;
                                   make_solid=true,make_ruled=false,
                                   max_degree=-1,continuity="",
                                   parametrization="",smoothing=false)
    caller="OCC.addThruSections"
    kernel=m.occ_kernel
    requested=_model_occ_requested_tag(tag,caller)
    dim=make_solid ? 3 : 2
    nsec=length(wire_tags)
    nsec>=2 || throw(ArgumentError(
        "$caller: ThruSections require at least 2 wires"))
    make_solid isa Bool && make_ruled isa Bool ||
        throw(ArgumentError("$caller: makeSolid/makeRuled must be Bool"))
    smoothing isa Bool || throw(ArgumentError(
        "$caller: smoothing must be Bool"))
    smoothing && throw(ArgumentError(
        "$caller: smoothing loft is not implemented"))
    max_degree isa Integer && !(max_degree isa Bool) || throw(ArgumentError(
        "$caller: maxDegree must be an integer"))
    continuity in ("","C0","C1","C2","C3","CN","G1","G2") ||
        throw(ArgumentError("$caller: unknown continuity '$continuity'"))
    parametrization in ("","ChordLength","Centripetal","IsoParametric") ||
        throw(ArgumentError(
            "$caller: unknown parametrization '$parametrization'"))
    allocated=dim==3 ?
        _model_occ_next_tag(kernel,3,requested,caller) :
        -1   # surface mode: every face gets its own tag; `tag` binds the first
    if dim==3
        haskey(kernel.bindings,(3,allocated)) && throw(ArgumentError(
            "$caller: OCC Volume[$allocated] already exists"))
    elseif requested>=0
        haskey(kernel.bindings,(2,requested)) && throw(ArgumentError(
            "$caller: OCC Surface[$requested] already exists"))
    end

    # ── Stage 1: resolve wires, read member uses (no mutation) ──
    uses_per_wire=Vector{Vector{Tuple{UInt64,Bool}}}(undef,nsec)
    for j in 1:nsec
        wt=wire_tags[j]
        wt isa Integer || throw(ArgumentError(
            "$caller: wire tags must be integers"))
        id=get(kernel.bindings,(-1,Int(wt)),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE wire or curve loop with tag $wt"))
        uses=kernel.shapes[id].geometry.uses
        isempty(uses) && throw(ArgumentError(
            "$caller: wire tag $wt has no edges"))
        uses_per_wire[j]=uses
    end
    ne=length(uses_per_wire[1])
    for j in 2:nsec
        length(uses_per_wire[j])==ne || throw(ArgumentError(
            "$caller: inconsistent number of edges in through sections"))
    end

    # Traversal-start vertex of each member edge; the traversal-end vertex of
    # the last edge detects closedness by shared shape identity (an OCC wire
    # is closed exactly when its end vertex IS the start vertex).
    vstart=[Vector{UInt64}(undef,ne) for _ in 1:nsec]
    vend=Vector{UInt64}(undef,nsec)
    closed=Vector{Bool}(undef,nsec)
    for j in 1:nsec
        for k in 1:ne
            eid,rev=uses_per_wire[j][k]
            a,b=kernel.shapes[eid].children
            vstart[j][k]=rev ? b : a
        end
        eid,rev=uses_per_wire[j][ne]
        a,b=kernel.shapes[eid].children
        vend[j]=rev ? a : b
        closed[j]=vend[j]==vstart[j][1]
    end
    if make_solid
        all(closed) || throw(ArgumentError(
            "$caller: Making solid requires closed wires"))
    end
    all(closed) || all(.!closed) || throw(ArgumentError(
        "$caller: cannot loft mixed open and closed sections"))
    ncols=closed[1] ? ne : ne+1
    vcol(j,c)= c<=ne ? vstart[j][c] : vend[j]

    # ── Stage 2: canonical directed nets + pcurve descriptors ──
    nets=[Vector{NURBS.NURBSCurve}(undef,ne) for _ in 1:nsec]
    umaps=[Vector{Any}(undef,ne) for _ in 1:nsec]
    centers=Vector{NTuple{3,Float64}}(undef,nsec)
    for j in 1:nsec
        cx=cy=cz=0.0; cnt=0
        for k in 1:ne
            eid,rev=uses_per_wire[j][k]
            g=kernel.shapes[eid].geometry.geometry
            nc,um=_loft_edge_nc(g,rev,caller)
            nets[j][k]=nc; umaps[j][k]=um
            for vid in kernel.shapes[eid].children
                p=kernel.shapes[vid].geometry.point
                cx+=p[1]; cy+=p[2]; cz+=p[3]; cnt+=1
            end
        end
        cnt>0 || throw(ArgumentError("$caller: wire has no vertices"))
        centers[j]=(cx/cnt,cy/cnt,cz/cnt)
    end

    # ── Stage 3: column unification + patches ──
    vbar=_loft_vparams(centers,parametrization,caller)
    qv=make_ruled ? 1 : min(3,nsec-1)
    if max_degree>0
        qv=min(qv,Int(max_degree))
        qv<1 && (qv=1)
    end
    # col_patches[k]: ruled → Vector of per-gap patches; smooth → single patch
    col_patches=Vector{Any}(undef,ne)
    col_vknots=Vector{Any}(undef,ne)
    for k in 1:ne
        unified=_loft_nc_unify([nets[j][k] for j in 1:nsec])
        pu=unified[1].degree
        knotsU=unified[1].knots
        nu=length(unified[1].controls)
        if make_ruled || qv==1
            patches=Vector{NURBS.NURBSSurface}(undef,nsec-1)
            for j in 1:nsec-1
                A=_loft_hcurve(unified[j]); B=_loft_hcurve(unified[j+1])
                ctrl=Matrix{NTuple{3,Float64}}(undef,nu,2)
                wts=Matrix{Float64}(undef,nu,2)
                for i in 1:nu
                    wa=A[i][4]; wb=B[i][4]
                    ctrl[i,1]=(A[i][1]/wa,A[i][2]/wa,A[i][3]/wa)
                    ctrl[i,2]=(B[i][1]/wb,B[i][2]/wb,B[i][3]/wb)
                    wts[i,1]=wa; wts[i,2]=wb
                end
                patches[j]=NURBS.NURBSSurface(pu,1,knotsU,
                    Float64[vbar[j],vbar[j],vbar[j+1],vbar[j+1]],ctrl,wts)
            end
            col_patches[k]=patches; col_vknots[k]=nothing
        else
            columns=[_loft_hcurve(unified[j]) for j in 1:nsec]
            rows,vknots=_loft_vinterp(columns,vbar,qv,caller)
            ctrl=Matrix{NTuple{3,Float64}}(undef,nu,nsec)
            wts=Matrix{Float64}(undef,nu,nsec)
            for i in 1:nu, j in 1:nsec
                h=rows[i][j]; w=h[4]
                ctrl[i,j]=(h[1]/w,h[2]/w,h[3]/w); wts[i,j]=w
            end
            col_patches[k]=NURBS.NURBSSurface(pu,qv,knotsU,vknots,ctrl,wts)
            col_vknots[k]=vknots
        end
    end

    # ── Stage 4: allocation (ids + tags through the sweep cursor) ──
    sw=_model_occ_sweep_setup((0.0,0.0,0.0),caller,kernel)
    # Seam edges: ruled → one per (gap, column); smooth → one per column.
    ruled=make_ruled || qv==1
    nseams=ruled ? (nsec-1)*ncols : ncols
    nfaces=ruled ? (nsec-1)*ne : ne
    seam_ids=UInt64[_model_occ_sweep_id!(sw,kernel) for _ in 1:nseams]
    seam_tags=Int[_model_occ_sweep_tag!(sw,kernel,1) for _ in 1:nseams]
    wire_ids=UInt64[_model_occ_sweep_id!(sw,kernel) for _ in 1:nfaces+2]
    wire_tags=Int[_model_occ_sweep_tag!(sw,kernel,-1) for _ in 1:nfaces+2]
    face_ids=UInt64[_model_occ_sweep_id!(sw,kernel) for _ in 1:nfaces+2]
    face_tags=Vector{Int}(undef,nfaces+2)
    # `tag` binds the solid (makeSolid) or the first lateral face. Bump the
    # dim-2 cursor past a requested face tag so later draws cannot collide.
    first_face_slot=1
    if !make_solid && requested>=0
        face_tags[1]=requested
        sw.cursor[5]=max(sw.cursor[5],requested)
        first_face_slot=2
    end
    for i in first_face_slot:nfaces+2
        face_tags[i]=_model_occ_sweep_tag!(sw,kernel,2)
    end
    if make_solid
        shell_id=_model_occ_sweep_id!(sw,kernel)
        shell_tag=_model_occ_sweep_tag!(sw,kernel,-2)
        solid_id=_model_occ_sweep_id!(sw,kernel)
        solid_tag=allocated
    end

    # ── Stage 5: records ──
    loft_centroid=(sum(c[1] for c in centers)/nsec,
                   sum(c[2] for c in centers)/nsec,
                   sum(c[3] for c in centers)/nsec)
    # seam_at(j,c): seam edge (id,tag) at gap j (ruled) or column c (smooth).
    seam_at(j,c)=ruled ? ((j-1)*ncols+c) : c
    # v-span helpers
    v_lo_of(j)=ruled ? vbar[j] : 0.0
    v_hi_of(j)=ruled ? vbar[j+1] : 1.0

    # Seam edge records.
    if ruled
        for j in 1:nsec-1, c in 1:ncols
            va=vcol(j,c); vb=vcol(j+1,c)
            pa=kernel.shapes[va].geometry.point
            pb=kernel.shapes[vb].geometry.point
            rec=_loft_line_record(pa,pb,caller)
            idx=seam_at(j,c)
            _model_occ_sweep_bind!(sw,1,seam_tags[idx],seam_ids[idx],
                _ModelOCCShape(1,(kind=:line,geometry=rec),UInt64[va,vb]))
        end
    else
        for c in 1:ncols
            va=vcol(1,c); vb=vcol(nsec,c)
            # The v-boundary interpolant of the adjacent column's patch:
            # vertex column c is ũ=0 of edge-column c (c<=ne) or ũ=1 of
            # edge-column ne (open wire's last column).
            if c<=ne
                patch=col_patches[c]; iu=1
            else
                patch=col_patches[ne]; iu=size(patch.controls,1)
            end
            poles=[patch.controls[iu,j] for j in 1:nsec]
            wts=[patch.weights[iu,j] for j in 1:nsec]
            nc=NURBS.NURBSCurve(patch.degree_v,patch.knots_v,poles,wts)
            cv=_ModelOCCNURBSCurve(:bspline,nc)
            lo,hi=_model_occ_curve_bounds(cv)
            rec=(occ=:copied_curve,curve=cv,t0=lo,t1=hi)
            idx=seam_at(1,c)
            _model_occ_sweep_bind!(sw,1,seam_tags[idx],seam_ids[idx],
                _ModelOCCShape(1,(kind=:tensor_bspline,geometry=rec),
                               UInt64[va,vb]))
        end
    end

    # Lateral faces.
    shape_of=id->_model_occ_sweep_shape_of(sw,kernel,id)
    tag_of=id->_model_occ_sweep_tag_of(sw,kernel,id,caller)
    centroid=loft_centroid
    face_of=Dict{Tuple{Int,Int},Int}()   # (gap_or_0, col) -> face_ids index
    lat_rev=Dict{Int,Bool}()             # face_ids index -> outward flip
    fi=0
    lat_order=Tuple{Int,Int}[]
    for k in 1:ne
        ngaps=ruled ? nsec-1 : 1
        for j in 1:ngaps
            fi+=1
            patch=ruled ? col_patches[k][j] : col_patches[k]
            ulo,ulo2=v_lo_of(j),v_hi_of(j)
            # Wire: bottom edge (ũ inc), right seam up, top edge reversed,
            # left seam down — CCW in (ũ,v).
            eid_lo,rev_lo=uses_per_wire[j][k]
            jhi=ruled ? j+1 : nsec
            eid_hi,rev_hi=uses_per_wire[jhi][k]
            kR= k==ne ? (closed[1] ? 1 : ne+1) : k+1
            seamR=ruled ? seam_at(j,kR) : seam_at(1,kR)
            seamL=ruled ? seam_at(j,k) : seam_at(1,k)
            uses=Tuple{UInt64,Bool}[(eid_lo,rev_lo),
                (seam_ids[seamR],false),(eid_hi,!rev_hi),
                (seam_ids[seamL],true)]
            # Coplanar straight-section laterals classify as planes — OCCT
            # likewise emits Geom_Plane for coplanar ruled sections, and
            # plane faces mesh through the standard (Delaunay-conforming)
            # planar path instead of the curved-patch tessellator. The
            # degree-1×1 gate is exact: a bilinear patch is planar iff its
            # boundary is, and higher-degree patches can bow out of a
            # coplanar boundary's plane.
            plangeo=patch.degree_u==1 && patch.degree_v==1 ? try
                _model_occ_plane_face_geometry(
                    shape_of,tag_of,[uses],caller)
            catch err
                err isa ArgumentError ? nothing : rethrow()
            end : nothing
            if plangeo!==nothing
                planar_rev=false
                if make_solid
                    mid=_occ_surface_point(
                        (patch=_ModelNURBSSurface(:bspline,patch),),
                        0.5,(ulo+ulo2)/2)
                    # Same mechanism as the caps: `reversed` flips the stored
                    # normal; reversing the wire uses would be cancelled by
                    # conic axis normalization.
                    planar_rev=_dot3(plangeo.axis,_arc_sub(mid,centroid))<0
                    planar_rev && (plangeo=(;plangeo...,reversed=true))
                end
                _model_occ_sweep_bind!(sw,-1,wire_tags[fi],wire_ids[fi],
                    _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses]))
                _model_occ_sweep_bind!(sw,2,face_tags[fi],face_ids[fi],
                    _ModelOCCShape(2,(kind=:plane,geometry=plangeo),
                                   UInt64[wire_ids[fi]]))
                face_of[(ruled ? j : 0,k)]=fi
                lat_rev[fi]=planar_rev
                push!(lat_order,(ruled ? j : 0,k))
                continue
            end
            # pcurves: bottom edge t→(ũ(t),v_lo), top t→(ũ(t),v_hi),
            # seams t→(u_side, v(t)).
            pc=Dict{Int,NamedTuple}()
            etag_lo=_model_occ_shape_tag(kernel,eid_lo,caller)
            etag_hi=_model_occ_shape_tag(kernel,eid_hi,caller)
            pc[etag_lo]=(fwd=_loft_umap_pcurve(umaps[j][k],ulo),rev=nothing)
            pc[etag_hi]=(fwd=_loft_umap_pcurve(umaps[jhi][k],ulo2),rev=nothing)
            for (sidx,u_side,reversed_use) in
                    ((seamL,0.0,true),(seamR,1.0,false))
                rec=sw.shapes[seam_ids[sidx]].geometry.geometry
                pcv=_loft_seam_pcurve(u_side,rec.t0,rec.t1,ulo,ulo2-ulo)
                stag=seam_tags[sidx]
                # Closed sections reuse one seam edge on both uv columns —
                # `fwd`/`rev` are selected by the use sign, like the swept
                # cylinder's two-sided seam encoding.
                ent=get(pc,stag,(fwd=nothing,rev=nothing))
                pc[stag]=reversed_use ? (fwd=ent.fwd,rev=pcv) :
                                        (fwd=pcv,rev=ent.rev)
            end
            # vertex_uv for `_occ_vertex_reparam_on_face` on patch faces.
            vuv=Dict{Int,NTuple{2,Float64}}()
            kR_=kR
            for (vid,uv) in ((vcol(j,k),(0.0,ulo)),
                             (vcol(j,kR_),(1.0,ulo)),
                             (vcol(jhi,k),(0.0,ulo2)),
                             (vcol(jhi,kR_),(1.0,ulo2)))
                vuv[_model_occ_shape_tag(kernel,vid,caller)]=uv
            end
            # Bake the outward normal into the record: volume assembly uses
            # `abs(signed_surface)`, so shell use-signs cannot reorient the
            # PLC — the stored `reversed` flag is what the curved mesher's
            # winding alignment (via `_brep_surface_normal`) honors.
            face_reversed=false
            if make_solid
                pmid=_occ_surface_point(
                    (patch=_ModelNURBSSurface(:bspline,patch),),
                    0.5,(ulo+ulo2)/2)
                _,fdu,fdv=_occ_surface_d1(
                    (patch=_ModelNURBSSurface(:bspline,patch),),
                    0.5,(ulo+ulo2)/2)
                fn=_occ_cross(fdu,fdv)
                face_reversed=_dot3(fn,_arc_sub(pmid,centroid))<0
            end
            geometry=(occ=:tensor_bspline,
                      patch=_ModelNURBSSurface(:bspline,patch),
                      pcurves=pc,vertex_uv=vuv,reversed=face_reversed)
            _model_occ_sweep_bind!(sw,-1,wire_tags[fi],wire_ids[fi],
                _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses]))
            _model_occ_sweep_bind!(sw,2,face_tags[fi],face_ids[fi],
                _ModelOCCShape(2,(kind=:tensor_bspline,geometry=geometry),
                               UInt64[wire_ids[fi]]))
            face_of[(ruled ? j : 0,k)]=fi
            lat_rev[fi]=face_reversed
            push!(lat_order,(ruled ? j : 0,k))
        end
    end

    # Caps + shell + solid when makeSolid.
    out=Tuple{Int,Int}[]
    if make_solid
        # Cap faces bind the section wires verbatim; each must be planar.
        cap_uses=[uses_per_wire[1],uses_per_wire[nsec]]
        caps=Tuple{Int,UInt64,Bool}[]
        for (which,uses) in enumerate(cap_uses)
            geo=_model_occ_plane_face_geometry(shape_of,tag_of,[uses],caller)
            outward=which==1 ?
                _arc_sub(centers[1],centroid) :
                _arc_sub(centers[nsec],centroid)
            # `_model_occ_plane_face_geometry` normalizes `axis` back to the
            # boundary conics' intrinsic normal, so reversing the wire uses
            # cannot bake an outward normal — the `reversed` record flag is
            # the mechanism (`_occ_surface_normal` and the winding alignment
            # honor it), matching `_model_box_plane_geometry`'s convention.
            cap_rev=_dot3(geo.axis,outward)<0
            cap_rev && (geo=(;geo...,reversed=true))
            fi+=1
            _model_occ_sweep_bind!(sw,-1,wire_tags[fi],wire_ids[fi],
                _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses]))
            _model_occ_sweep_bind!(sw,2,face_tags[fi],face_ids[fi],
                _ModelOCCShape(2,(kind=:plane,geometry=geo),
                               UInt64[wire_ids[fi]]))
            push!(caps,(fi,face_ids[fi],cap_rev))
        end
        # The measure path orients each face's canonical normal by the shell
        # use-sign alone (σ·axis / σ·(du×dv)), so every use carries that
        # face's outward flip — `rev` equals the record's `reversed` flag
        # here, the same pairing `occ.extrude` emits on its source cap.
        shell_uses=Tuple{UInt64,Bool}[]
        for (fi_,fid,cap_rev) in caps
            push!(shell_uses,(fid,cap_rev))
        end
        for (j,k) in lat_order
            fid=face_ids[face_of[(j,k)]]
            push!(shell_uses,(fid,lat_rev[face_of[(j,k)]]))
        end
        _model_occ_sweep_bind!(sw,-2,shell_tag,shell_id,
            _ModelOCCShape(-2,(uses=shell_uses,),
                           UInt64[u[1] for u in shell_uses]))
        _model_occ_sweep_bind!(sw,3,solid_tag,solid_id,
            _ModelOCCShape(3,(kind=:solid,),UInt64[shell_id]))
        push!(out,(3,solid_tag))
    else
        # Surface mode: `_multiBind` explores FACE — one (2,tag) per lateral,
        # in creation order; `tag` already bound to the first face.
        for (j,k) in lat_order
            push!(out,(2,face_tags[face_of[(j,k)]]))
        end
    end

    for (dimension,t,id,shape) in sw.pending
        _model_occ_bind_shape!(kernel,dimension,t,id,shape)
    end
    isempty(sw.pending) ||
        (kernel.next_shape=max(kernel.next_shape,sw.next_id))
    return out
end

function _loft_line_record(a,b,caller)
    d=_arc_sub(b,a)
    len=_occ_modulus(d)
    len>0.0 || throw(ArgumentError(
        "$caller: loft collapsed on coincident section vertices"))
    return (occ=:line,origin=a,dir=_occ_dir(d),t0=0.0,t1=len)
end

# ── `addSurfaceFilling` / `addBSplineFilling` ─────────────────────────────────
#
# `BRepOffsetAPI_MakeFilling` contract: a bounded surface whose border is the
# given wire. Planar wires publish the exact `:plane` record — upstream's
# filling degenerates to the containing plane likewise. Non-planar 4-edge
# wires become a Coons/Gordon-blended tensor B-spline patch: the two ruled
# families plus the corner bilinear interpolant, sampled on the unified
# boundary Greville grid and re-interpolated in the same knot spaces so the
# boundary is reproduced exactly. Member edges keep `lin2d`/`tan2d`/`tan2dv`
# pcurves on the patch rectangle.

# Homogeneous evaluation of a rescaled NURBSCurve at u ∈ [0,1] via
# `_loft_basis_row` — returns (xw,yw,zw,w) so the Coons blend below stays
# exact for rational boundaries (conic arcs).
function _fill_nc_eval(nc,u)
    knots=nc.knots; p=nc.degree; n=length(nc.controls)
    u=clamp(u,knots[1],knots[end])
    s=clamp(searchsortedlast(knots,u),p+1,length(knots)-p-1)
    row=_loft_basis_row(knots,p,s,u,n)
    X=Y=Z=W=0.0
    @inbounds for i in 1:n
        b=row[i]
        b==0.0 && continue
        w=nc.weights[i]; c=nc.controls[i]
        X+=b*w*c[1]; Y+=b*w*c[2]; Z+=b*w*c[3]; W+=b*w
    end
    W!=0.0 || throw(ArgumentError(
        "OCC.addSurfaceFilling: degenerate filling boundary"))
    return (X,Y,Z,W)
end

# u/v edge umap flip: t ↦ 1−ū(t) — `(:lin,off,sc)` → `(:lin,1−off,−sc)` and
# the conic `(:tan,inner,off,sc)` flips identically in the outer map.
_fill_flip_umap(um)=um[1]===:lin ? (:lin,1.0-um[2],-um[3]) :
                                    (:tan,um[2],1.0-um[3],-um[4])

# v-side pcurve for edge with umap `um` on column `uu`: t ↦ (uu, v0+vs·t)
# for linear edges, (uu, v0+vs·inner(t)) for conic edges.
function _fill_vmap_pcurve(um,uu)
    if um[1]===:lin
        return _loft_lin2d((uu,um[2]),(0.0,um[3]))
    end
    return (tan2dv=(um[2],uu,um[3],um[4]),)
end

# Unify helper: Greville abscissae of a clamped knot vector.
function _fill_greville(knots,p,n)
    out=Vector{Float64}(undef,n)
    for i in 1:n
        s=0.0
        for k in 1:p
            s+=knots[i+k]
        end
        out[i]=s/p
    end
    return out
end

# Coons tensor patch for a 4-edge closed wire. Edges in traversal order:
# e1 bottom (u fwd), e2 right (v fwd), e3 top (u bwd), e4 left (v bwd).
function _fill_coons_geometry(kernel,uses,tag_of,caller)
    ncs=Vector{NURBS.NURBSCurve}(undef,4); ums=Vector{Any}(undef,4)
    for k in 1:4
        eid,rev=uses[k]
        g=kernel.shapes[eid].geometry.geometry
        ncs[k],ums[k]=_loft_edge_nc(g,rev,caller)
    end
    nc_bot,nc_top=ncs[1],_loft_nc_reverse(ncs[3])
    nc_rgt,nc_lft=ncs[2],_loft_nc_reverse(ncs[4])
    uc=_loft_nc_unify([nc_bot,nc_top]); vc=_loft_nc_unify([nc_lft,nc_rgt])
    pu=uc[1].degree; knotsU=uc[1].knots; nu=length(uc[1].controls)
    qv=vc[1].degree; knotsV=vc[1].knots; nv=length(vc[1].controls)
    û=_fill_greville(knotsU,pu,nu); v̂=_fill_greville(knotsV,qv,nv)
    cb=u->_fill_nc_eval(nc_bot,u); ct=u->_fill_nc_eval(nc_top,u)
    cl=v->_fill_nc_eval(nc_lft,v); cr=v->_fill_nc_eval(nc_rgt,v)
    P00=cb(0.0); P10=cb(1.0); P01=ct(0.0); P11=ct(1.0)
    # Coons blend in homogeneous space — the corner bilinear term subtracts
    # the doubly-counted corner patch in (xw,yw,zw,w) too, keeping the
    # boundary curves' rational parametrization exactly representable.
    G=Matrix{NTuple{4,Float64}}(undef,nu,nv)
    for i in 1:nu, j in 1:nv
        u=û[i]; v=v̂[j]
        pb=cb(u); pt=ct(u); pl=cl(v); pr=cr(v)
        w00=(1-u)*(1-v); w10=u*(1-v); w01=(1-u)*v; w11=u*v
        G[i,j]=ntuple(d->(1-v)*pb[d]+v*pt[d]+(1-u)*pl[d]+u*pr[d]-
                (w00*P00[d]+w10*P10[d]+w01*P01[d]+w11*P11[d]),4)
    end
    # Tensor interpolation: C = A_u⁻¹·G·A_v⁻ᵀ with basis rows at Greville pts.
    Au=Matrix{Float64}(undef,nu,nu); Av=Matrix{Float64}(undef,nv,nv)
    for i in 1:nu
        s=clamp(searchsortedlast(knotsU,û[i]),pu+1,length(knotsU)-pu-1)
        Au[i,:].=_loft_basis_row(knotsU,pu,s,û[i],nu)
    end
    for j in 1:nv
        s=clamp(searchsortedlast(knotsV,v̂[j]),qv+1,length(knotsV)-qv-1)
        Av[j,:].=_loft_basis_row(knotsV,qv,s,v̂[j],nv)
    end
    # Interpolate the homogeneous grid; dehomogenize into poles+weights.
    ctrl=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    wts=Matrix{Float64}(undef,nu,nv)
    H=[(Au\Float64[G[i,j][d] for i in 1:nu, j in 1:nv])/transpose(Av)
        for d in 1:4]                     # per-dim A_u⁻¹·G·A_v⁻ᵀ
    for i in 1:nu, j in 1:nv
        w=H[4][i,j]
        w>0.0 || throw(ArgumentError(
            "OCC.addSurfaceFilling: degenerate filling patch"))
        ctrl[i,j]=(H[1][i,j]/w,H[2][i,j]/w,H[3][i,j]/w)
        wts[i,j]=w
    end
    patch=NURBS.NURBSSurface(pu,qv,knotsU,knotsV,ctrl,wts)
    # pcurves — bottom e1 u-trace at v=0, right e2 v-trace at u=1, top e3
    # u-trace at v=1 with flipped direction, left e4 v-trace at u=0 flipped.
    pc=Dict{Int,NamedTuple}()
    eids=[uses[k][1] for k in 1:4]
    pc[tag_of(eids[1])]=(fwd=_loft_umap_pcurve(ums[1],0.0),rev=nothing)
    pc[tag_of(eids[2])]=(fwd=_fill_vmap_pcurve(ums[2],1.0),rev=nothing)
    pc[tag_of(eids[3])]=(fwd=_loft_umap_pcurve(
        _fill_flip_umap(ums[3]),1.0),rev=nothing)
    pc[tag_of(eids[4])]=(fwd=_fill_vmap_pcurve(
        _fill_flip_umap(ums[4]),0.0),rev=nothing)
    # vertex_uv: traversal corners.
    vuv=Dict{Int,NTuple{2,Float64}}()
    for (k,uv) in ((1,(0.0,0.0)),(2,(1.0,0.0)),(3,(1.0,1.0)),(4,(0.0,1.0)))
        eid,rev=uses[k]
        a,b=kernel.shapes[eid].children
        vid=rev ? b : a
        vuv[_model_occ_shape_tag(kernel,vid,caller)]=uv
    end
    return (occ=:tensor_bspline,
            patch=_ModelNURBSSurface(:bspline,patch),
            pcurves=pc,vertex_uv=vuv,reversed=false)
end

# `OCC_Internals::addSurfaceFilling` — the filling surface's tag binds dim 2;
# sub-entities stay owned by the input wire (no recursive rebind).
function _model_occ_surface_filling(m::GeoModel,tag,wire_tag;
                                    caller="OCC.addSurfaceFilling")
    kernel=m.occ_kernel
    requested=_model_occ_requested_tag(tag,caller)
    allocated=_model_occ_next_tag(kernel,2,requested,caller)
    haskey(kernel.bindings,(2,allocated)) && throw(ArgumentError(
        "$caller: OCC Surface[$allocated] already exists"))
    wt=abs(_model_occ_requested_tag(wire_tag,caller))
    wid=get(kernel.bindings,(-1,wt),nothing)
    wid===nothing && throw(ArgumentError("$caller: unknown OCC Wire[$wt]"))
    uses=kernel.shapes[wid].geometry.uses
    length(uses)>=3 || throw(ArgumentError(
        "$caller: filling wire must have at least 3 edges"))
    shape_of=id->kernel.shapes[id]
    tag_of=id->_model_occ_shape_tag(kernel,id,caller)
    geometry=try
        _model_occ_plane_face_geometry(shape_of,tag_of,[uses],caller)
    catch err
        err isa ArgumentError ? nothing : rethrow()
    end
    if geometry===nothing
        length(uses)==4 || throw(ArgumentError(
            "$caller: non-planar filling requires a 4-edge wire"))
        geometry=_fill_coons_geometry(kernel,uses,tag_of,caller)
    end
    kind=haskey(geometry,:patch) ? :tensor_bspline : :plane
    shape=_ModelOCCShape(2,(kind=kind,geometry=geometry),UInt64[wid])
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,2,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end
