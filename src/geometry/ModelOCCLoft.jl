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
            fi+=1
            _model_occ_sweep_bind!(sw,-1,wire_tags[fi],wire_ids[fi],
                _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses]))
            # The planar mesher admits a cap only when its boundary is
            # BITWISE coplanar (`orient3 == 0` on every sample — a rotated
            # transport frame leaves ~ulp residuals the exact predicate
            # rejects even though the cap is genuinely planar).  When the
            # fit plane leaves any residual, publish the same disk as a
            # trimmed degree-1 tensor patch over the identical (X,Y) frame
            # so the curved mesher lifts its pcurve boundary instead — the
            # geometry is truthful, just classified differently.
            cap_samples=_model_occ_wire_samples(shape_of,uses,caller)
            bitwise=all(
                p->_dot3(_arc_sub(p,geo.center),geo.axis)==0.0,cap_samples)
            if bitwise
                cap_rev && (geo=(;geo...,reversed=true))
                _model_occ_sweep_bind!(sw,2,face_tags[fi],face_ids[fi],
                    _ModelOCCShape(2,(kind=:plane,geometry=geo),
                                   UInt64[wire_ids[fi]]))
            else
                ulo,uhi=geo.uvb[1][1],geo.uvb[2][1]
                vlo,vhi=geo.uvb[1][2],geo.uvb[2][2]
                ctrl=Matrix{NTuple{3,Float64}}(undef,2,2)
                wts=ones(Float64,2,2)
                for (i,uu) in enumerate((ulo,uhi)),
                        (j,vv) in enumerate((vlo,vhi))
                    p=geo.center .+ uu.*geo.X .+ vv.*geo.Y
                    ctrl[i,j]=(p[1],p[2],p[3])
                end
                patch=NURBS.NURBSSurface(1,1,Float64[ulo,ulo,uhi,uhi],
                    Float64[vlo,vlo,vhi,vhi],ctrl,wts)
                vuv=Dict{Int,NTuple{2,Float64}}()
                for (eid,_) in uses, vid in shape_of(eid).children
                    o=_arc_sub(shape_of(vid).geometry.point,geo.center)
                    vuv[tag_of(vid)]=(_dot3(o,geo.X),_dot3(o,geo.Y))
                end
                # `plane_rec` keeps the exact fitted-plane record for the
                # measures path — the trimmed-CDT quadrature cannot reach
                # plane-shoelace precision on conic wires, and the
                # reclassification is a meshing admission, not a geometry
                # claim.
                tgeo=(occ=:tensor_bspline,
                      patch=_ModelNURBSSurface(:bspline,patch),
                      pcurves=geo.pcurves,vertex_uv=vuv,reversed=cap_rev,
                      trimmed=true,plane_rec=geo)
                _model_occ_sweep_bind!(sw,2,face_tags[fi],face_ids[fi],
                    _ModelOCCShape(2,(kind=:tensor_bspline,geometry=tgeo),
                                   UInt64[wire_ids[fi]]))
            end
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

# ── `addPipe` ─────────────────────────────────────────────────────────────────
#
# `BRepOffsetAPI_MakePipe` semantics: a wire or face profile swept along a
# spine wire. Upstream sweeps the profile exactly; here the profile is
# transported rigidly along discrete parallel-transport frames — stations
# densify with spine curvature (π/16 on conics, 16 spans on free-form) —
# and the transported section wires loft through `addThruSections`. The
# approximation converges to the true pipe under station refinement and is
# exact for straight spines (two-circle stations yield an exact ruled
# rational cylinder).

function _pipe_edge_point(m,g,etag,u)
    g.occ===:line && return g.origin .+ u .* g.dir
    g.occ===:circle && return _occ_circle_point(g,u)
    g.occ===:ellipse && return _occ_ellipse_point(g,u)
    g.occ===:copied_curve && return _model_occ_curve_d2(g.curve,u)[1]
    g.occ===:surface_iso && return _model_surface_iso_point(g,u)
    throw(ArgumentError(
        "OCC.addPipe: unsupported spine edge kind $(g.occ)"))
end

# Station parameter cuts on one spine edge, ends included.
function _pipe_edge_stations(g)
    if g.occ===:line
        return Float64[g.t0,g.t1]
    elseif g.occ in (:circle,:ellipse)
        k=max(4,ceil(Int,abs(g.t1-g.t0)/(π/16)))
        return Float64[g.t0+(g.t1-g.t0)*i/k for i in 0:k]
    end
    return Float64[g.t0+(g.t1-g.t0)*i/16 for i in 0:16]
end

# Row-major 3×3 product A·B over NTuple{9} layouts.
_pipe_matmul(A,B)=ntuple(k->begin
    r=(k-1)÷3; c=(k-1)%3
    A[3r+1]*B[c+1]+A[3r+2]*B[c+4]+A[3r+3]*B[c+7]
end,9)

# Minimal rotation taking unit a → unit b (axis = a×b, 180° fallback to any
# perpendicular). Row-major NTuple{9}.
function _pipe_min_rotation(a,b)
    ax=_occ_cross(a,b); c=_dot3(a,b)
    s=_occ_modulus(ax)
    if s<1e-12
        c>0.0 && return (1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0)
        # anti-parallel: rotate π about any perpendicular
        p=abs(a[1])<0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
        ax=_occ_dir(_occ_cross(a,p))
        x,y,z=ax
        return (1-2*(y*y+z*z), 2*x*y, 2*x*z,
                2*x*y, 1-2*(x*x+z*z), 2*y*z,
                2*x*z, 2*y*z, 1-2*(x*x+y*y))
    end
    x,y,z=ax./s
    return _occ_rotation_matrix((x,y,z),atan(s,c))
end

function _model_occ_add_pipe!(m::GeoModel,tag,profile,wire_tag)
    caller="OCC.addPipe"
    kernel=m.occ_kernel
    make_solid=false
    profile_uses=Vector{Tuple{UInt64,Bool}}()
    if profile isa Tuple && profile[1]==2
        fid=get(kernel.bindings,(2,profile[2]),nothing)
        fid===nothing && throw(ArgumentError(
            "$caller: unknown OCC Surface[$(profile[2])]"))
        wires=kernel.shapes[fid].children
        isempty(wires) && throw(ArgumentError(
            "$caller: profile face has no boundary"))
        profile_uses=kernel.shapes[wires[1]].geometry.uses
        make_solid=true
    else
        wt=profile isa Tuple ? abs(Int(profile[2])) : abs(Int(profile))
        wid=get(kernel.bindings,(-1,wt),nothing)
        wid===nothing && throw(ArgumentError(
            "$caller: unknown OCC Wire[$wt]"))
        profile_uses=kernel.shapes[wid].geometry.uses
    end
    isempty(profile_uses) && throw(ArgumentError(
        "$caller: profile has no edges"))
    wt=abs(_model_occ_requested_tag(wire_tag,caller))
    swid=get(kernel.bindings,(-1,wt),nothing)
    swid===nothing && throw(ArgumentError("$caller: unknown OCC Wire[$wt]"))
    suses=kernel.shapes[swid].geometry.uses
    isempty(suses) && throw(ArgumentError("$caller: spine wire has no edges"))
    # ── spine stations (position + unit tangent), traversal order ──
    P=NTuple{3,Float64}[]; T=NTuple{3,Float64}[]
    for (eid,rev) in suses
        g=kernel.shapes[eid].geometry.geometry
        etag=_model_occ_shape_tag(kernel,eid,caller)
        params=_pipe_edge_stations(g)
        rev && reverse!(params)
        for u in params
            p=_pipe_edge_point(m,g,etag,u)
            d=_occ_curve_first_der(m,g,etag,u,caller)
            tn=_occ_modulus(d)
            tn>_OCC_CONFUSION || continue
            dir=_occ_dir(rev ? (.-d) : d)
            # A station coincident with the previous (shared edge end)
            # only refreshes the outgoing tangent for frame transport.
            if !isempty(P) && _model_occ_vertices_coincide(p,last(P))
                T[end]=dir
                continue
            end
            push!(P,p); push!(T,dir)
        end
    end
    nsec=length(P)
    nsec>=2 || throw(ArgumentError(
        "$caller: spine yields fewer than 2 sections"))
    # ── cumulative parallel-transport frames ──
    R=Vector{NTuple{9,Float64}}(undef,nsec)
    R[1]=(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0)
    for i in 2:nsec
        Ri=_pipe_min_rotation(T[i-1],T[i])
        R[i]=_pipe_matmul(Ri,R[i-1])
    end
    P0=P[1]
    # ── transported section wires ──
    ne=length(profile_uses)
    # Per section the loft binds one shape per profile edge, one per DISTINCT
    # profile vertex (a closed wire shares `ne`; an open chain carries `ne+1`)
    # plus the section wire itself — the old `2ne+1` budget only covered the
    # closed case.
    nverts=length(Set{UInt64}(vid for (peid,_) in profile_uses
                            for vid in kernel.shapes[peid].children))
    ids=_model_occ_reserve_shapes(kernel,nsec*(nverts+ne+1),caller)
    next=1
    section_tags=Int[]
    ecache=Dict{UInt64,UInt64}()   # profile vertex id → station vertex id
    for i in 1:nsec
        off=(P[i][1]-_occ_matvec(R[i],P0)[1],
             P[i][2]-_occ_matvec(R[i],P0)[2],
             P[i][3]-_occ_matvec(R[i],P0)[3])
        t=_AffineTransform(R[i],off,NTuple{16,Float64}[],
                           (linear=R[i],loc=off))
        empty!(ecache)
        uses=Vector{Tuple{UInt64,Bool}}(undef,ne)
        for (k,(peid,prev)) in enumerate(profile_uses)
            pshape=kernel.shapes[peid]
            newg=_model_occ_transform_geometry(m,pshape,t,caller,
                                               "profile edge")
            kids=UInt64[]
            for vid in pshape.children
                nv_id=get!(ecache,vid) do
                    vshape=kernel.shapes[vid]
                    vg=_model_occ_transform_geometry(m,vshape,t,caller,
                                                     "profile vertex")
                    t0=_model_occ_next_tag(kernel,0,-1,caller)
                    _model_occ_bind_shape!(kernel,0,t0,ids[next],
                        _ModelOCCShape(0,vg,UInt64[]))
                    next+=1
                    ids[next-1]
                end
                push!(kids,nv_id)
            end
            etag=_model_occ_next_tag(kernel,1,-1,caller)
            _model_occ_bind_shape!(kernel,1,etag,ids[next],
                _ModelOCCShape(1,newg,kids))
            next+=1
            uses[k]=(ids[next-1],prev)
        end
        wtag=_model_occ_next_tag(kernel,-1,-1,caller)
        _model_occ_bind_shape!(kernel,-1,wtag,ids[next],
            _ModelOCCShape(-1,(uses=uses,),UInt64[u[1] for u in uses]))
        next+=1
        push!(section_tags,wtag)
    end
    kernel.next_shape=ids[next-1]+UInt64(1)
    return _model_occ_thru_sections!(m,tag,section_tags;
        make_solid=make_solid,make_ruled=false)
end

# ── `addTrimmedSurface` ───────────────────────────────────────────────────────
#
# `BRepBuilderAPI_MakeFace(surface, wire)` semantics: the wire's edges are
# projected onto the base surface's canonical domain, producing a bounded
# face on that surface. Edges projecting linearly in uv collapse to exact
# `lin2d` pcurves; the rest become `spl2d` — cubic B-spline interpolants
# through unwrapped (period-corrected) projection samples. Without a
# surfaceTag the wire resolves through the plane path.

# (u,v) parameter periods of a trim support: the torus additionally wraps
# in v, and tensor patches report their knot-domain width in whichever
# direction the seam is geometrically closed (`_brep_uperiod` covers u).
function _trim_uv_periods(face_g)
    face_g.occ===:plane && return 0.0,0.0
    uperiod=_brep_uperiod(face_g)
    vperiod=if face_g.occ===:torus
        _OCC_TWO_PI
    elseif face_g.occ in (:tensor_bspline,:tensor_bezier)
        _model_nurbs_patch_seam_period(face_g,2)
    else
        0.0
    end
    return uperiod,vperiod
end

# Parameter-domain limits used to detect off-domain trim projections:
# tensor patches are defined only over their knot domain; bounded OCC
# faces carry `uvb`. Returns `(lo,hi)` or `nothing` when unknown.
function _trim_domain_bounds(face_g)
    face_g.occ in (:tensor_bspline,:tensor_bezier) &&
        return _model_nurbs_surface_bounds(face_g)
    return get(face_g,:uvb,nothing)
end

# Project edge `eid` onto base record `face_g`; returns a lin2d or spl2d
# pcurve over the edge's own parameter range.
function _trim_project_pcurve(m,kernel,face_g,eid,caller)
    g=kernel.shapes[eid].geometry.geometry
    etag=_model_occ_shape_tag(kernel,eid,caller)
    t0,t1=g.t0,g.t1
    ns=17
    ts=Float64[t0+(t1-t0)*i/(ns-1) for i in 0:ns-1]
    uvs=Vector{NTuple{2,Float64}}(undef,ns)
    period,period_v=_trim_uv_periods(face_g)
    prev=nothing
    for (i,u) in enumerate(ts)
        P=_pipe_edge_point(m,g,etag,u)
        pr=_occ_surface_closest(face_g,P,1.0)
        pr===nothing && throw(ArgumentError(
            "$caller: Curve[$etag] does not project onto the surface"))
        uv=(pr[1],pr[2])
        # Off-surface edges legitimately project (3-D chords between
        # on-surface points become spl2d traces), but a wire whose
        # projection misses the domain must not silently clamp onto its
        # boundary: reject only when the closest uv is pinned at a
        # NON-periodic domain bound while the edge point stays off the
        # surface. Planes are unbounded — their `uvb` is a sampling box,
        # never a trim limit.
        if face_g.occ!==:plane &&
           _occ_modulus(_arc_sub(pr[3],P))>1e-7*max(1.0,maximum(abs,P))
            db=_trim_domain_bounds(face_g)
            if db!==nothing
                dlo,dhi=db
                for d in 1:2
                    (d==1 ? period : period_v)>0.0 && continue
                    span=dhi[d]-dlo[d]; span>0 || continue
                    (uv[d]-dlo[d]<1e-9*span || dhi[d]-uv[d]<1e-9*span) &&
                        throw(ArgumentError(
                            "$caller: Curve[$etag] does not project within the surface domain"))
                end
            end
        end
        if prev!==nothing
            period>0.0 && (uv=(uv[1]+period*round((prev[1]-uv[1])/period),
                               uv[2]))
            period_v>0.0 && (uv=(uv[1],
                uv[2]+period_v*round((prev[2]-uv[2])/period_v)))
        end
        prev=uv; uvs[i]=uv
    end
    # Exact linear reduction: every sample must sit on the t-interpolated
    # chord — a straight uv trace needs d·t linear in the edge parameter.
    u0,uN=uvs[1],uvs[end]
    du=(uN[1]-u0[1])/(t1-t0),(uN[2]-u0[2])/(t1-t0)
    scale=max(1.0,maximum(abs(u[1]-u0[1]) for u in uvs),
              maximum(abs(u[2]-u0[2]) for u in uvs))
    lin=true
    for i in 2:ns-1
        pu=u0[1]+du[1]*(ts[i]-t0); pv=u0[2]+du[2]*(ts[i]-t0)
        (abs(uvs[i][1]-pu)>1e-9*scale || abs(uvs[i][2]-pv)>1e-9*scale) &&
            (lin=false; break)
    end
    # `lin2d` evaluates `o+d·t` at the ABSOLUTE edge parameter — store
    # the uv at t=0, not at t0 (imported edges can start at t0<0).
    lin && return (lin2d=((u0[1]-du[1]*t0,u0[2]-du[2]*t0),du),)
    return (spl2d=_pchip_bspline_2d(ts,uvs),)
end

# Monotone (Fritsch–Carlson) piecewise-cubic Hermite interpolant through
# (ts[i], uvs[i]) converted to a C1 cubic B-spline: interior knots carry
# multiplicity 2, the three Bezier poles per segment collapse the shared
# endpoint. Monotonicity keeps the trace inside the sampled uv range — a
# global interpolant can overshoot past the surface's v-window.
function _pchip_bspline_2d(ts,uvs)
    n=length(ts)
    us=[p[1] for p in uvs]; vs=[p[2] for p in uvs]
    function slopes(ys)
        m=Vector{Float64}(undef,n)
        d=[(ys[i+1]-ys[i])/(ts[i+1]-ts[i]) for i in 1:n-1]
        # PCHIP endpoints: one-sided extrapolation clamped to the first
        # secant sign, capped at 3|d1| when the secants disagree.
        function end_slope(h1,h2,d1,d2)
            e=((2*h1+h2)*d1-h1*d2)/(h1+h2)
            sign(e)!=sign(d1) && return 0.0
            (sign(d1)!=sign(d2) && abs(e)>3*abs(d1)) && return 3*d1
            return e
        end
        if n==2
            m[1]=d[1]; m[2]=d[1]
        else
            m[1]=end_slope(ts[2]-ts[1],ts[3]-ts[2],d[1],d[2])
            m[n]=end_slope(ts[n]-ts[n-1],ts[n-1]-ts[n-2],
                           d[n-1],d[n-2])
        end
        for i in 2:n-1
            if d[i-1]*d[i]<=0.0
                m[i]=0.0
            else
                w1=2*(ts[i+1]-ts[i])+(ts[i]-ts[i-1])
                w2=(ts[i+1]-ts[i])+2*(ts[i]-ts[i-1])
                m[i]=(w1+w2)/(w1/d[i-1]+w2/d[i])
            end
        end
        return m
    end
    mu=slopes(us); mv=slopes(vs)
    return _hermite_bspline_layout(ts,us,vs,mu,mv)
end

# Shared piecewise-Bezier layout for Hermite-family 2-D interpolants:
# interior knots at full multiplicity p=3 (C0 in the B-spline
# parameterization, tangent-continuous in geometry), junction pole stored
# once per side.
function _hermite_bspline_layout(ts,us,vs,mu,mv)
    n=length(ts)
    knots=Float64[]; cps=NTuple{2,Float64}[]
    for i in 1:n-1
        h=ts[i+1]-ts[i]
        p0=(us[i],vs[i]); p3=(us[i+1],vs[i+1])
        p1=(us[i]+mu[i]*h/3,vs[i]+mv[i]*h/3)
        p2=(us[i+1]-mu[i+1]*h/3,vs[i+1]-mv[i+1]*h/3)
        i>1 && append!(knots,[ts[i],ts[i],ts[i]])
        i==1 && append!(knots,[ts[1],ts[1],ts[1],ts[1]])
        push!(cps,p0); push!(cps,p1); push!(cps,p2)
        i==n-1 && (push!(cps,p3);
                   append!(knots,[ts[end],ts[end],ts[end],ts[end]]))
    end
    return (knots,3,cps)
end

# Cubic Hermite interpolant through (ts[i], uvs[i]) with PRESCRIBED
# per-knot derivatives `duv[i]=(du/dt,dv/dt)` — a NaN entry falls back to
# the centred secant slope so callers can mark points where the exact
# tangent is unavailable (e.g. a vanishing level gradient). Exact
# tangents converge several orders faster than the estimated pchip
# slopes on a strongly curved trace.
function _hermite_bspline_2d(ts,uvs,duv)
    n=length(ts)
    us=[p[1] for p in uvs]; vs=[p[2] for p in uvs]
    function fallback(i)
        lo=i==1 ? 1 : i-1
        hi=i==n ? n : i+1
        h=ts[hi]-ts[lo]
        return h>0.0 ? ((us[hi]-us[lo])/h,(vs[hi]-vs[lo])/h) : (0.0,0.0)
    end
    mu=Vector{Float64}(undef,n); mv=Vector{Float64}(undef,n)
    for i in 1:n
        d=duv[i]
        if isnan(d[1]) || isnan(d[2])
            mu[i],mv[i]=fallback(i)
        else
            mu[i],mv[i]=d
        end
    end
    return _hermite_bspline_layout(ts,us,vs,mu,mv)
end

function _model_occ_trimmed_surface(m::GeoModel,tag,wire_tag,surface_tag;
                                    caller="OCC.addTrimmedSurface")
    kernel=m.occ_kernel
    requested=_model_occ_requested_tag(tag,caller)
    allocated=_model_occ_next_tag(kernel,2,requested,caller)
    haskey(kernel.bindings,(2,allocated)) && throw(ArgumentError(
        "$caller: OCC Surface[$allocated] already exists"))
    wt=abs(_model_occ_requested_tag(wire_tag,caller))
    wid=get(kernel.bindings,(-1,wt),nothing)
    wid===nothing && throw(ArgumentError("$caller: unknown OCC Wire[$wt]"))
    uses=kernel.shapes[wid].geometry.uses
    isempty(uses) && throw(ArgumentError("$caller: wire has no edges"))
    shape_of=id->kernel.shapes[id]
    tag_of=id->_model_occ_shape_tag(kernel,id,caller)
    shape=nothing
    if surface_tag isa Integer && surface_tag>=0
        st=Int(surface_tag)
        fid=get(kernel.bindings,(2,st),nothing)
        fid===nothing && throw(ArgumentError(
            "$caller: unknown OCC Surface[$st]"))
        base=kernel.shapes[fid].geometry
        base_g=base.geometry
        base_g.occ in (:cylinder,:sphere,:cone,:torus,:plane,
                       :tensor_bspline,:tensor_bezier) ||
            throw(ArgumentError(
                "$caller: Surface[$st] has no parametric OCC record"))
        # The base face's own member edges carry pcurves already — reuse
        # them when the trim wire re-uses those edges; otherwise project.
        pc=Dict{Int,NamedTuple}()
        for (eid,_) in uses
            etag=tag_of(eid)
            own=get(base_g,:pcurves,nothing)
            own_pc=own===nothing ? nothing : get(own,etag,nothing)
            pc[etag]=own_pc===nothing ?
                (fwd=_trim_project_pcurve(m,kernel,base_g,eid,caller),
                 rev=nothing) : own_pc
        end
        vuv=Dict{Int,NTuple{2,Float64}}()
        for (eid,_) in uses, vid in kernel.shapes[eid].children
            vtag=tag_of(vid)
            haskey(vuv,vtag) && continue
            pt=kernel.shapes[vid].geometry.point
            pr=_occ_surface_closest(base_g,pt,1.0)
            pr===nothing && throw(ArgumentError(
                "$caller: wire vertex does not project onto the surface"))
            vuv[vtag]=(pr[1],pr[2])
        end
        geometry=merge(base_g,(pcurves=pc,vertex_uv=vuv,trimmed=true))
        if base_g.occ===:plane
            # `uvb` is the canonical sampling domain for plane records —
            # rebox it on the new boundary so stale bounds cannot clip the
            # trim.
            umin=vmin=Inf; umax=vmax=-Inf
            for (eid,_) in uses
                eg=kernel.shapes[eid].geometry.geometry
                lo,hi=_occ_pcurve_fillbndbox(pc[tag_of(eid)].fwd,
                                             eg.t0,eg.t1)
                umin=min(umin,lo[1]); vmin=min(vmin,lo[2])
                umax=max(umax,hi[1]); vmax=max(vmax,hi[2])
            end
            geometry=merge(geometry,(uvb=([umin,vmin],[umax,vmax]),))
        end
        shape=_ModelOCCShape(2,merge(base,(geometry=geometry,)),UInt64[wid])
    else
        geometry=_model_occ_plane_face_geometry(shape_of,tag_of,[uses],caller)
        shape=_ModelOCCShape(2,(kind=:plane,geometry=geometry),UInt64[wid])
    end
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,2,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# Bind a closed wire from a precomputed `uses` list — same storage as
# `_model_occ_add_wire!` but the caller already knows each edge's `rev`
# flag (shared-edge assemblies where forced-forward chaining cannot
# traverse the loop).
function _occ_wire_from_uses!(m::GeoModel,uses,caller)
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,-1,-1,caller)
    shape=_ModelOCCShape(-1,(uses=uses,),UInt64[id for (id,_) in uses])
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,-1,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# ── `addThickSolid` ───────────────────────────────────────────────────────────
#
# `BRepOffsetAPI_MakeThickSolid`: hollows a planar-faced solid by offsetting
# every non-excluded face by `offset` (positive = outward, negative =
# inward) and capping the excluded faces' openings. The result is built
# directly as boundary shells — no boolean is needed once the two parallel
# bodies are known:
#
#   offset < 0  outer wall = original kept faces, inner wall = offset-body
#               faces wound CW about n̂ (wall outward = −n̂), rim = original
#               excluded wire + offset-body polygon as hole
#   offset > 0  outer wall = offset-body faces (+n̂), inner wall = original
#               kept faces consumed with flipped shell rev, rim = offset
#               polygon as outer wire + reversed original wire as hole
#
# Restricted to planar-faced solids — curved faces reject loudly (analytic
# curved offsets remain with the 3-D fillet backlog).
function _model_occ_thick_solid(m::GeoModel,tag,volume_tag,exclude_tags,
                                offset;caller="OCC.addThickSolid")
    kernel=m.occ_kernel
    vt=abs(_model_occ_requested_tag(volume_tag,caller))
    vid=get(kernel.bindings,(3,vt),nothing)
    vid===nothing && throw(ArgumentError("$caller: unknown OCC Volume[$vt]"))
    offset==0.0 && throw(ArgumentError("$caller: offset must be non-zero"))
    excl=Set(abs(_model_occ_requested_tag(t,caller)) for t in exclude_tags)
    # Per-face-use outward planes n̂·x ≤ b. Plane duplicates (two uses of
    # the same face tag can only appear once in a closed shell, so fuses
    # of distinct faces sharing a plane stay separate entries) keep the
    # full face-use bookkeeping for shell reuse.
    planes=NTuple{3,Float64}[]; bs=Float64[]; kept=Bool[]
    fuses=Tuple{UInt64,Bool}[]; ftags=Int[]
    in_volume=Set{Int}()
    for shell_id in kernel.shapes[vid].children
        for (fid,rev) in kernel.shapes[shell_id].geometry.uses
            fg=kernel.shapes[fid].geometry
            get(fg.geometry,:occ,nothing)!==:plane && throw(ArgumentError(
                "$caller: addThickSolid needs planar faces (got occ=" *
                "$(get(fg.geometry,:occ,nothing)))"))
            ftag=_model_occ_shape_tag(kernel,fid,caller)
            push!(in_volume,ftag)
            σ=rev ? -1.0 : 1.0
            n=_occ3_scale(fg.geometry.axis,σ)
            get(fg.geometry,:reversed,false) && (n=(.-n))
            n=_occ3_scale(n,1.0/_occ_modulus(n))
            push!(planes,n); push!(bs,_dot3(n,fg.geometry.center))
            push!(kept,!(ftag in excl))
            push!(fuses,(fid,rev)); push!(ftags,ftag)
        end
    end
    for t in excl
        t in in_volume || throw(ArgumentError(
            "$caller: Surface[$t] is not a face of Volume[$vt]"))
    end
    # Offset half-spaces: kept planes shift outward by `offset`; excluded
    # planes stay and bound the opening.
    bb=Float64[kept[i] ? bs[i]+offset : bs[i] for i in eachindex(bs)]
    P=length(planes)
    P<4 && throw(ArgumentError("$caller: offset body degenerates"))
    # Orientation sanity: the claimed-original halfspace body must be
    # non-empty — a flipped face normal empties or inverts it.
    # Boolean-produced shells store faces inside-out wholesale
    # (rev+reversed bookkeeping), so a fully empty body retries once
    # with every plane flipped.
    orig_verts=_occ_halfspace_verts(planes,bs)
    flipped=false
    if isempty(orig_verts)
        planes=[(.-n) for n in planes]; bs=(.-bs)
        bb=Float64[kept[i] ? bs[i]+offset : bs[i] for i in eachindex(bs)]
        orig_verts=_occ_halfspace_verts(planes,bs)
        flipped=true
    end
    isempty(orig_verts) && throw(ArgumentError(
        "$caller: Volume[$vt] is not oriented outward"))
    # Offset-body vertices: triple intersections inside every shifted
    # half-space, deduplicated within Confusion (4+ planes can share a
    # corner).
    verts=_occ_halfspace_verts(planes,bb)
    isempty(verts) && throw(ArgumentError(
        "$caller: offset produces an empty body"))
    # Halfspace-intersection offsets are exact only for convex bodies —
    # a nonconvex input leaves phantom corners in the halfspace body.
    # Reject any enumerated corner that is not a real vertex of the
    # source solid (collected from its wire endpoints).
    realv=NTuple{3,Float64}[]
    for (fid,_) in fuses
        for wid in kernel.shapes[fid].children
            for (eid,_) in kernel.shapes[wid].geometry.uses
                for vid2 in kernel.shapes[eid].children
                    push!(realv,kernel.shapes[vid2].geometry.point)
                end
            end
        end
    end
    for x in orig_verts
        any(v->_model_occ_vertices_coincide(v,x),realv) || throw(
            ArgumentError("$caller: needs a convex planar solid — " *
                          "nonconvex offsets are not supported"))
    end
    if offset<0
        for x in verts, q in 1:P
            _dot3(planes[q],x)>bs[q]+1e-6 && throw(ArgumentError(
                "$caller: Volume[$vt] is not oriented outward"))
        end
    else
        for x in orig_verts, q in 1:P
            kept[q] || continue
            _dot3(planes[q],x)>bb[q]+1e-6 && throw(ArgumentError(
                "$caller: Volume[$vt] is not oriented outward"))
        end
    end
    # Per-plane polygon (member vertices sorted CCW about n̂); shared
    # points/lines publish once so every curve lands on two surfaces.
    ptags=[_model_occ_add_point!(m,x[1],x[2],x[3]) for x in verts]
    edge_of=Dict{NTuple{2,Int},Tuple{UInt64,Int}}()  # (a,b)→(line id, stored stop vert)
    polygon(i)=begin
        on=Int[]
        for (vi,x) in enumerate(verts)
            abs(_dot3(planes[i],x)-bb[i])<1e-7 && push!(on,vi)
        end
        length(on)<3 && return on
        n=planes[i]
        ax=abs(n[1])<0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
        e1=_occ3_scale(_occ_cross(n,ax),1.0)
        e1=_occ3_scale(e1,1.0/_occ_modulus(e1))
        e2=_occ_cross(n,e1)
        c=_occ3_scale(reduce(_occ3_add,(verts[v] for v in on)),
                      1.0/length(on))
        sort!(on,by=v->atan(
            _dot3(_occ3_sub(verts[v],c),e2),
            _dot3(_occ3_sub(verts[v],c),e1)))
        on
    end
    function poly_wire(on,flip)   # CCW-about-n̂ `on`; flip winds it CW
        uses=Tuple{UInt64,Bool}[]
        K=length(on)
        for j in 1:K
            a,b=on[j],on[mod1(j+1,K)]
            key=a<b ? (a,b) : (b,a)
            eu=get(edge_of,key,nothing)
            if eu===nothing
                et=_model_occ_add_line!(m,ptags[a],ptags[b])
                eu=(kernel.bindings[(1,et)],b)   # stored stop = verts[b]
                edge_of[key]=eu
            end
            # traversal must END at verts[b]; rev when the stored stop is
            # the other end. `flip` reverses both order and orientation.
            rev=eu[2]!=b
            push!(uses,flip ? (eu[1],!rev) : (eu[1],rev))
        end
        flip && (uses=reverse(uses))
        _occ_wire_from_uses!(m,uses,caller)
    end
    polys=[polygon(i) for i in 1:P]
    # Winding rule: outer-wall/rim outer wires run CCW about n̂, cavity
    # walls CW — all of it mirrored when the input read inside-out (the
    # result then carries the same convention as its source).
    sgn=flipped ? -1 : 1
    cw(role)=role*sgn<0
    off_face=Dict{Int,Int}()        # plane index -> new face tag
    for i in 1:P
        length(polys[i])<3 && continue
        kept[i] || continue
        w=poly_wire(polys[i],cw(offset<0 ? -1 : 1))
        off_face[i]=_model_occ_add_plane_surface!(m,[w])
    end
    # Rim faces on excluded planes: outer wire + hole wire on the same
    # plane, wound oppositely. `wires` come from the *outer* side of the
    # wall (offset>0 → offset polygon, offset<0 → the original wire).
    rim_face=Int[]
    for i in 1:P
        kept[i] && continue
        fid,_=fuses[i]
        outer0=kernel.shapes[fid].children[1]          # original outer wire
        outer0_uses=kernel.shapes[outer0].geometry.uses
        if offset<0
            length(polys[i])>=3 || throw(ArgumentError(
                "$caller: excluded Surface[$(ftags[i])] does not meet " *
                "the offset body"))
            # hole wound opposite the (reused) original outer wire —
            # CW about n̂ normally, CCW when the source read inside-out
            w=poly_wire(polys[i],!flipped)
            face=_model_occ_add_plane_surface!(m,
                [_model_occ_shape_tag(kernel,outer0,caller),w])
        else
            length(polys[i])>=3 || throw(ArgumentError(
                "$caller: excluded Surface[$(ftags[i])] does not meet " *
                "the offset body"))
            wo=poly_wire(polys[i],cw(1))               # outer wire
            wh=_occ_wire_from_uses!(m,
                flipped ? collect(outer0_uses) :
                    [(id,!r) for (id,r) in reverse(outer0_uses)],caller)
            face=_model_occ_add_plane_surface!(m,[wo,wh])
        end
        push!(rim_face,face)
    end
    # Shell assembly. `add_surface_loop` bakes rev=false, so records are
    # bound directly wherever flipped uses are needed. With no opening the
    # result is two disconnected shells (outer wall + cavity); an opening
    # joins everything into one shell through the rim faces.
    shell_tags=Int[]
    if isempty(excl)
        # outer shell: offset>0 → offset faces (+n̂); offset<0 → original
        outer_uses=Tuple{UInt64,Bool}[]
        if offset>0
            for (_,ft) in off_face
                push!(outer_uses,(kernel.bindings[(2,ft)],false))
            end
        else
            append!(outer_uses,fuses)
        end
        # cavity shell: offset>0 → original faces flipped (−n̂ into the
        # void); offset<0 → offset faces (axis already −n̂).
        inner_uses=Tuple{UInt64,Bool}[]
        if offset>0
            for (fid,rev) in fuses
                push!(inner_uses,(fid,!rev))
            end
        else
            for (_,ft) in off_face
                push!(inner_uses,(kernel.bindings[(2,ft)],false))
            end
        end
        push!(shell_tags,_occ_shell_from_uses!(m,outer_uses,caller))
        push!(shell_tags,_occ_shell_from_uses!(m,inner_uses,caller))
    else
        shell_uses=Tuple{UInt64,Bool}[]
        if offset<0
            for i in 1:P
                kept[i] && push!(shell_uses,fuses[i])
            end
            for (_,ft) in off_face
                push!(shell_uses,(kernel.bindings[(2,ft)],false))
            end
        else
            for (_,ft) in off_face
                push!(shell_uses,(kernel.bindings[(2,ft)],false))
            end
            for i in 1:P                        # cavity walls: original
                kept[i] && push!(shell_uses,(fuses[i][1],!fuses[i][2]))
            end
        end
        for ft in rim_face
            push!(shell_uses,(kernel.bindings[(2,ft)],false))
        end
        length(shell_uses)<4 && throw(ArgumentError(
            "$caller: thick shell has too few faces"))
        push!(shell_tags,_occ_shell_from_uses!(m,shell_uses,caller))
    end
    _model_occ_add_volume!(m,shell_tags;tag)
end

# Vertices of a half-space intersection {x : n̂ᵢ·x ≤ bᵢ} — triple-plane
# intersections kept inside every half-space, deduplicated within
# Confusion.
function _occ_halfspace_verts(planes,bs)
    P=length(planes)
    verts=NTuple{3,Float64}[]
    for i in 1:P-2, j in i+1:P-1, k in j+1:P
        n1,n2,n3=planes[i],planes[j],planes[k]
        d=_dot3(n1,_occ_cross(n2,n3))
        abs(d)<1e-14 && continue
        x=_occ3_scale(
            _occ3_add(_occ3_add(
                _occ3_scale(_occ_cross(n2,n3),bs[i]),
                _occ3_scale(_occ_cross(n3,n1),bs[j])),
                _occ3_scale(_occ_cross(n1,n2),bs[k])),
            1.0/d)
        all(isfinite,x) || continue
        all(_dot3(planes[q],x)<=bs[q]+1e-7 for q in 1:P) || continue
        any(v->_occ_modulus(_occ3_sub(v,x))<=1e-7,verts) && continue
        push!(verts,(x[1],x[2],x[3]))
    end
    return verts
end

# ── `defeature` ───────────────────────────────────────────────────────────────
#
# `BRepAlgoAPI_Defeature`: remove the listed faces of a planar-faced
# solid and refill the gap — the classic use is dropping a boss or
# pocket's faces so the surrounding face planes close over it. The
# defeatured body is the half-space intersection of the *kept* faces'
# outward planes: exact whenever that result is convex (boss-on-box →
# the plain box; the pocket case likewise). Newly intersected corners
# are published through the same shared-edge path `addThickSolid` uses;
# kept faces are re-bound fresh from their polygons (half-space
# membership can grow their extent past the original wire).
function _model_occ_defeature(m::GeoModel,volume_tag,remove_tags;
                              caller="OCC.defeature")
    kernel=m.occ_kernel
    vt=abs(_model_occ_requested_tag(volume_tag,caller))
    vid=get(kernel.bindings,(3,vt),nothing)
    vid===nothing && throw(ArgumentError("$caller: unknown OCC Volume[$vt]"))
    rem=Set(abs(_model_occ_requested_tag(t,caller)) for t in remove_tags)
    isempty(rem) && throw(ArgumentError(
        "$caller: no faces given to remove"))
    planes=NTuple{3,Float64}[]; bs=Float64[]; in_volume=Set{Int}()
    for shell_id in kernel.shapes[vid].children
        for (fid,rev) in kernel.shapes[shell_id].geometry.uses
            fg=kernel.shapes[fid].geometry
            get(fg.geometry,:occ,nothing)!==:plane && throw(ArgumentError(
                "$caller: defeature needs planar faces (got occ=" *
                "$(get(fg.geometry,:occ,nothing)))"))
            ftag=_model_occ_shape_tag(kernel,fid,caller)
            push!(in_volume,ftag)
            ftag in rem && continue
            σ=rev ? -1.0 : 1.0
            n=_occ3_scale(fg.geometry.axis,σ)
            get(fg.geometry,:reversed,false) && (n=(.-n))
            n=_occ3_scale(n,1.0/_occ_modulus(n))
            push!(planes,n); push!(bs,_dot3(n,fg.geometry.center))
        end
    end
    for t in rem
        t in in_volume || throw(ArgumentError(
            "$caller: Surface[$t] is not a face of Volume[$vt]"))
    end
    P=length(planes)
    P<4 && throw(ArgumentError(
        "$caller: removing $(length(rem)) faces leaves no closed body"))
    verts=_occ_halfspace_verts(planes,bs)
    if isempty(verts)
        # Boolean-produced shells store faces inside-out (rev+reversed
        # bookkeeping): retry with every plane flipped.
        planes=[(.-n) for n in planes]; bs=(.-bs)
        verts=_occ_halfspace_verts(planes,bs)
        isempty(verts) && throw(ArgumentError(
            "$caller: kept faces bound an empty body"))
    end
    # Old sub-shape ids become orphans once the volume rebinds.
    old_ids=UInt64[]
    stack=copy(kernel.shapes[vid].children)
    while !isempty(stack)
        x=pop!(stack); x in old_ids && continue
        push!(old_ids,x); append!(stack,kernel.shapes[x].children)
        hasproperty(kernel.shapes[x].geometry,:uses) &&
            append!(stack,(u for (u,_) in
                           kernel.shapes[x].geometry.uses))
    end
    shell=_occ_solid_from_planes!(m,planes,bs,verts,caller)
    # `OCC_Internals::defeature` mutates the volume in place — rebind the
    # same tag to the healed shell rather than allocating a new volume.
    sid=kernel.bindings[(-2,shell)]
    kernel.shapes[vid]=_ModelOCCShape(3,(kind=:solid,),UInt64[sid])
    _occ_unbind_orphans!(kernel,Set(old_ids))
    _occ_invalidate_volume_tables!(m,vt)
    return vt
end

# A volume whose shell was rebuilt in place keeps its kernel binding
# but its imported model records — the `m.volumes` surface-loop list
# plus any primitive or boolean encodings — still describe the
# pre-rebuild boundary.  Drop them so the next `synchronize`
# re-imports the actual shell instead of serving stale topology to
# `get_boundary` or the analytic mesh paths.
function _occ_invalidate_volume_tables!(m::GeoModel,tag::Int)
    delete!(m.volumes,tag)
    for tab in (m.cylinders,m.spheres,m.cones,m.box_extents,
                m.booleans,m.boolean_operands,m.boolean_components)
        delete!(tab,tag)
    end
    return nothing
end

# Unbind entities in `dead` that nothing outside `dead` references —
# safe post-rebind cleanup for consumed sub-shapes.
function _occ_unbind_orphans!(kernel,dead::Set{UInt64})
    referenced=Set{UInt64}()
    for (sid,sh) in kernel.shapes
        sid in dead && continue
        union!(referenced,sh.children)
        hasproperty(sh.geometry,:uses) &&
            union!(referenced,(u for (u,_) in sh.geometry.uses))
    end
    # A dead-listed shape still referenced by a survivor (e.g. a wire
    # shared with a face that was not rebuilt) keeps protecting its
    # own references — propagate to a fixpoint.
    grown=true
    while grown
        grown=false
        for sid in intersect(dead,referenced)
            sh=kernel.shapes[sid]
            for x in sh.children
                x in referenced || (push!(referenced,x); grown=true)
            end
            hasproperty(sh.geometry,:uses) && for (u,_) in sh.geometry.uses
                u in referenced || (push!(referenced,u); grown=true)
            end
        end
    end
    for id in dead
        id in referenced && continue
        entry=get(kernel.shape_tags,id,nothing)
        entry===nothing && continue
        _model_occ_entity_map_unbind!(
            kernel.entity_maps[first(entry)+3],last(entry))
        delete!(kernel.bindings,entry);delete!(kernel.shape_tags,id)
        push!(kernel.removed,entry)
    end
    kernel.changed=true
    return nothing
end

# Build a closed solid from an outward half-space set: per-plane
# polygons of the enumerated body vertices, shared point/line
# entities, one shell. Returns the surface-loop tag.
function _occ_solid_from_planes!(m::GeoModel,planes,bs,verts,caller)
    kernel=m.occ_kernel
    ptags=[_model_occ_add_point!(m,x[1],x[2],x[3]) for x in verts]
    edge_of=Dict{NTuple{2,Int},Tuple{UInt64,Int}}()
    face_tags=Int[]
    for i in eachindex(planes)
        on=Int[]
        for (vi,x) in enumerate(verts)
            abs(_dot3(planes[i],x)-bs[i])<1e-7 && push!(on,vi)
        end
        length(on)<3 && continue
        n=planes[i]
        ax=abs(n[1])<0.9 ? (1.0,0.0,0.0) : (0.0,1.0,0.0)
        e1=_occ3_scale(_occ_cross(n,ax),1.0)
        e1=_occ3_scale(e1,1.0/_occ_modulus(e1))
        e2=_occ_cross(n,e1)
        c=_occ3_scale(reduce(_occ3_add,(verts[v] for v in on)),
                      1.0/length(on))
        sort!(on,by=v->atan(
            _dot3(_occ3_sub(verts[v],c),e2),
            _dot3(_occ3_sub(verts[v],c),e1)))
        uses=Tuple{UInt64,Bool}[]
        for j in eachindex(on)
            a,b=on[j],on[mod1(j+1,length(on))]
            key=a<b ? (a,b) : (b,a)
            eu=get(edge_of,key,nothing)
            if eu===nothing
                et=_model_occ_add_line!(m,ptags[a],ptags[b])
                eu=(kernel.bindings[(1,et)],b)
                edge_of[key]=eu
            end
            push!(uses,(eu[1],eu[2]!=b))
        end
        w=_occ_wire_from_uses!(m,uses,caller)
        push!(face_tags,_model_occ_add_plane_surface!(m,[w]))
    end
    length(face_tags)<4 && throw(ArgumentError(
        "$caller: kept faces do not close a solid"))
    return _model_occ_add_surface_loop!(m,face_tags)
end

# Bind a surface loop from precomputed `(face id, rev)` uses — parallel to
# `_model_occ_add_surface_loop!` for shells whose member faces need
# explicit orientation flags (cavity walls reusing source faces).
function _occ_shell_from_uses!(m::GeoModel,uses,caller)
    kernel=m.occ_kernel
    allocated=_model_occ_next_tag(kernel,-2,-1,caller)
    shape=_ModelOCCShape(-2,(uses=uses,),UInt64[id for (id,_) in uses])
    id=only(_model_occ_reserve_shapes(kernel,1,caller))
    _model_occ_bind_shape!(kernel,-2,allocated,id,shape)
    kernel.next_shape=id+UInt64(1)
    return allocated
end

# ── `chamfer` ─────────────────────────────────────────────────────────────────
#
# `BRepFilletAPI_MakeChamfer` on planar-faced solids: each chamfered
# edge contributes one new half-space — the plane through the lines
# parallel to the edge at distances d1/d2 inside its two adjacent
# faces — so the chamfered body is exactly the original half-space
# intersection plus those planes (same enumeration as `defeature`).
# The plane must clip only the edge's corner region: its midpoint is
# forced onto the removed side, and any *other* original vertex on the
# removed side means the chamfer overruns a neighboring corner — a
# loud rejection like upstream's "chamfer could not be computed".
function _model_occ_chamfer(m::GeoModel,volume_tags,curve_tags,
                            distance1,distance2;
                            caller="OCC.chamfer")
    kernel=m.occ_kernel
    out=Tuple{Int,Int}[]
    for volume_tag in volume_tags
        vt=abs(_model_occ_requested_tag(volume_tag,caller))
        vid=get(kernel.bindings,(3,vt),nothing)
        vid===nothing && throw(ArgumentError(
            "$caller: unknown OCC Volume[$vt]"))
        # Outward planes + edge→adjacent-face map.
        planes=NTuple{3,Float64}[]; bs=Float64[]
        fnormals=Dict{UInt64,NTuple{3,Float64}}()
        fcenters=Dict{UInt64,NTuple{3,Float64}}()
        edge_faces=Dict{UInt64,Vector{UInt64}}()
        for shell_id in kernel.shapes[vid].children
            for (fid,rev) in kernel.shapes[shell_id].geometry.uses
                fg=kernel.shapes[fid].geometry
                get(fg.geometry,:occ,nothing)!==:plane &&
                    throw(ArgumentError(
                        "$caller: chamfer needs planar faces (got " *
                        "occ=$(get(fg.geometry,:occ,nothing)))"))
                σ=rev ? -1.0 : 1.0
                n=_occ3_scale(fg.geometry.axis,σ)
                get(fg.geometry,:reversed,false) && (n=(.-n))
                n=_occ3_scale(n,1.0/_occ_modulus(n))
                fnormals[fid]=n; fcenters[fid]=fg.geometry.center
                push!(planes,n); push!(bs,_dot3(n,fg.geometry.center))
                for wid in kernel.shapes[fid].children
                    wsh=kernel.shapes[wid]
                    hasproperty(wsh.geometry,:uses) || continue
                    for (u,_) in wsh.geometry.uses
                        push!(get!(edge_faces,u,UInt64[]),fid)
                    end
                end
            end
        end
        overs=_occ_halfspace_verts(planes,bs)
        if isempty(overs)
            planes=[(.-n) for n in planes]; bs=(.-bs)
            fnormals=Dict(k=>(.-v) for (k,v) in fnormals)
            overs=_occ_halfspace_verts(planes,bs)
            isempty(overs) && throw(ArgumentError(
                "$caller: Volume[$vt] faces bound an empty body"))
        end
        # Chamfer half-space per requested edge.
        nplanes=copy(planes); nbs=copy(bs)
        for curve_tag in curve_tags
            et=abs(_model_occ_requested_tag(curve_tag,caller))
            eid=get(kernel.bindings,(1,et),nothing)
            eid===nothing && throw(ArgumentError(
                "$caller: unknown OCC Curve[$et]"))
            fids=get(edge_faces,eid,nothing)
            (fids===nothing || length(fids)!=2) && throw(ArgumentError(
                "$caller: Curve[$et] is not an edge shared by two " *
                "faces of Volume[$vt]"))
            eg=kernel.shapes[eid].geometry.geometry
            get(eg,:occ,nothing)==:line || throw(ArgumentError(
                "$caller: chamfer edge must be a line " *
                "(Curve[$et] is occ=$(get(eg,:occ,nothing)))"))
            t=_occ3_scale(eg.dir,1.0/_occ_modulus(eg.dir))
            em=_occ3_add(eg.origin,
                         _occ3_scale(eg.dir,(eg.t0+eg.t1)/2))
            fA,fB=fids
            nA,nB=fnormals[fA],fnormals[fB]
            # In-face inward directions perpendicular to the edge.
            uA=_occ_cross(t,nA); uB=_occ_cross(t,nB)
            uA=_occ3_scale(uA,1.0/_occ_modulus(uA))
            uB=_occ3_scale(uB,1.0/_occ_modulus(uB))
            _dot3(uA,_occ3_sub(fcenters[fA],em))<0 && (uA=(.-uA))
            _dot3(uB,_occ3_sub(fcenters[fB],em))<0 && (uB=(.-uB))
            qA=_occ3_add(em,_occ3_scale(uA,distance1))
            qB=_occ3_add(em,_occ3_scale(uB,distance2))
            w=_occ3_sub(qB,qA)
            nc=_occ_cross(t,w)
            nl=_occ_modulus(nc)
            nl<1e-12 && throw(ArgumentError(
                "$caller: degenerate chamfer on Curve[$et]"))
            nc=_occ3_scale(nc,1.0/nl)
            bc=_dot3(nc,qA)
            # Corner side must violate the kept half-space.
            _dot3(nc,em)>bc || (nc=(.-nc); bc=-bc)
            # No original vertex other than the edge's own endpoints
            # may cross the cut — else the chamfer overruns a corner.
            ends=[kernel.shapes[c].geometry.point
                  for c in kernel.shapes[eid].children]
            for x in overs
                _dot3(nc,x)>bc+1e-7 || continue
                any(e->_occ_modulus(_occ3_sub(e,x))<=1e-7,ends) &&
                    continue
                throw(ArgumentError(
                    "$caller: chamfer on Curve[$et] overruns a " *
                    "neighboring corner of Volume[$vt]"))
            end
            push!(nplanes,nc); push!(nbs,bc)
        end
        verts=_occ_halfspace_verts(nplanes,nbs)
        isempty(verts) && throw(ArgumentError(
            "$caller: chamfer leaves an empty body"))
        # Old sub-shape ids become orphans once the volume rebinds —
        # collected here, unbound below if nothing else uses them.
        old_ids=UInt64[]
        stack=copy(kernel.shapes[vid].children)
        while !isempty(stack)
            x=pop!(stack); x in old_ids && continue
            push!(old_ids,x); append!(stack,kernel.shapes[x].children)
            hasproperty(kernel.shapes[x].geometry,:uses) &&
                append!(stack,(u for (u,_) in
                               kernel.shapes[x].geometry.uses))
        end
        shell=_occ_solid_from_planes!(m,nplanes,nbs,verts,caller)
        sid=kernel.bindings[(-2,shell)]
        kernel.shapes[vid]=_ModelOCCShape(3,(kind=:solid,),UInt64[sid])
        _occ_unbind_orphans!(kernel,Set(old_ids))
        _occ_invalidate_volume_tables!(m,vt)
        kernel.changed=true
        push!(out,(3,vt))
    end
    return out
end


# ── circle-edge fillets (ChFi3d on coaxial surfaces of revolution) ────────────
#
# A full-circle fillet edge shared by two faces of revolution that share the
# edge's axis gets a rolling-ball fillet which is itself a surface of
# revolution: the torus whose meridian is the 2-D fillet of the two faces'
# meridian profiles in the (r,z) half-plane — r = distance from the edge
# axis, z = height along +axis measured from the edge centre. Plane,
# cylinder and cone faces give line profiles; sphere and torus faces give
# circle profiles. The tube centre is the intersection of the two profiles
# offset by the fillet radius toward the material; the feet of
# perpendiculars are the two tangent circles that replace the rim on each
# kept face.

@inline _fillet_dot2(a,b)=a[1]*b[1]+a[2]*b[2]
@inline _fillet_sub2(a,b)=(a[1]-b[1],a[2]-b[2])
@inline _fillet_sc2(a,s)=(a[1]*s,a[2]*s)
@inline _fillet_add2(a,b)=(a[1]+b[1],a[2]+b[2])
@inline _fillet_perp(v)=(-v[2],v[1])
@inline _fillet_norm2(v)=sqrt(_fillet_dot2(v,v))
@inline _fillet_cross(a,b)=a[1]*b[2]-a[2]*b[1]

# Point of a kernel edge record at parameter t (strip-overrun sampling).
function _occ_fillet_edge_point(g,t)
    k=get(g,:occ,nothing)
    k===:line && return _occ_line_point(g,t)
    k===:circle && return _occ_circle_point(g,t)
    k===:copied_curve && return first(_model_occ_curve_d2(g.curve,t))
    throw(ArgumentError(
        "OCC.fillet: fillet strip sampling needs line or circle boundary " *
        "edges (got $k)"))
end

# Meridian profile of the face record `sg` along its rim with the filleted
# circle edge. `σ` is the shell-outward orientation of the stored record
# (use-rev × `reversed` flag). All coordinates live in the rim's (r,z)
# half-plane where the rim sits at (R,0): `:line` profiles carry the
# profile tangent `d`, `:circ` profiles carry centre `O` and `rad`; `nin`
# is the material-inward unit normal of the profile at the rim.
function _occ_fillet_rim_profile(sg,σ,axis,o,R,et,caller)
    rim=(R,0.0)
    tol=1e-7*max(1.0,R)
    function onaxis(p)
        d=_arc_sub(p,o)
        _occ_modulus(_arc_sub(d,_occ3_scale(axis,_dot3(d,axis))))<tol
    end
    k=get(sg,:occ,nothing)
    if k===:plane
        abs(_dot3(sg.axis,axis))>1-1e-9 || throw(ArgumentError(
            "$caller: Curve[$et] rim lies in a plane not perpendicular " *
            "to the edge axis"))
        abs(_dot3(_arc_sub(o,sg.center),sg.axis))<tol || throw(
            ArgumentError(
                "$caller: Curve[$et] rim does not lie in its incident " *
                "plane"))
        nz=σ*_dot3(sg.axis,axis)
        return (kind=:line,P=rim,d=(1.0,0.0),nin=(0.0,-nz),sg=sg,σ=σ)
    elseif k===:cylinder
        (abs(_dot3(sg.axis,axis))>1-1e-9 && onaxis(sg.center)) || throw(
            ArgumentError("$caller: Curve[$et] touches a cylinder that " *
                          "is not coaxial with the edge"))
        abs(sg.radius-R)<tol || throw(ArgumentError(
            "$caller: Curve[$et] is not a rim of its incident cylinder"))
        return (kind=:line,P=rim,
                d=(0.0,sign(_dot3(sg.axis,axis))),nin=(-σ,0.0),sg=sg,σ=σ)
    elseif k===:cone
        (abs(_dot3(sg.axis,axis))>1-1e-9 && onaxis(sg.center)) || throw(
            ArgumentError("$caller: Curve[$et] touches a cone that is " *
                          "not coaxial with the edge"))
        hr=_dot3(_arc_sub(o,sg.center),sg.axis)
        sa=atan((sg.r2-sg.r1)/sg.height)
        σc=sign(_dot3(sg.axis,axis))
        abs(sg.r1+hr*tan(sa)-R)<tol || throw(ArgumentError(
            "$caller: Curve[$et] is not a rim of its incident cone"))
        cs,ss=cos(sa),sin(sa)
        return (kind=:line,P=rim,d=(ss,σc*cs),nin=(-σ*cs,σ*σc*ss),
                sg=sg,σ=σ)
    elseif k===:sphere
        onaxis(sg.center) || throw(ArgumentError(
            "$caller: Curve[$et] touches a sphere that is not coaxial " *
            "with the edge"))
        hz=_dot3(_arc_sub(sg.center,o),axis)
        abs(R*R+hz*hz-sg.radius*sg.radius)<tol*max(1.0,sg.radius) ||
            throw(ArgumentError(
                "$caller: Curve[$et] is not a rim of its incident sphere"))
        return (kind=:circ,P=rim,O=(0.0,hz),rad=sg.radius,
                nin=(-σ*R/sg.radius,σ*hz/sg.radius),sg=sg,σ=σ)
    elseif k===:torus
        onaxis(sg.center) || throw(ArgumentError(
            "$caller: Curve[$et] touches a torus that is not coaxial " *
            "with the edge"))
        hc=_dot3(_arc_sub(sg.center,o),axis); dr=R-sg.r1
        abs(dr*dr+hc*hc-sg.r2*sg.r2)<tol*max(1.0,sg.r2) || throw(
            ArgumentError(
                "$caller: Curve[$et] is not a rim of its incident torus"))
        return (kind=:circ,P=rim,O=(sg.r1,hc),rad=sg.r2,
                nin=(-σ*dr/sg.r2,σ*hc/sg.r2),sg=sg,σ=σ)
    end
    throw(ArgumentError(
        "$caller: Curve[$et] is a circle edge shared by an occ=$k face " *
        "— only coaxial plane/cylinder/cone/sphere/torus rims are " *
        "supported"))
end

# Offset-intersection solve: the tube centre C of the radius-ρ fillet is a
# point at inward distance ρ from both profiles, inside the notch (the
# wedge where both inward normals apply) and with its tangent feet on the
# domain side of each rim. Nearest valid candidate to the rim wins.
function _occ_fillet_meridian(pA,pB,R,ρ,et,caller)
    rim=(R,0.0); tol=1e-9*max(1.0,R,ρ)
    nA,nB=pA.nin,pB.nin
    cands=NTuple{2,Float64}[]
    if pA.kind===:line && pB.kind===:line
        det=_fillet_cross(nA,nB)
        abs(det)>1e-12 || throw(ArgumentError(
            "$caller: Curve[$et] has parallel rim profiles — no " *
            "fillet centre"))
        b1=_fillet_dot2(nA,rim)+ρ; b2=_fillet_dot2(nB,rim)+ρ
        push!(cands,((b1*nB[2]-b2*nA[2])/det,
                     (nA[1]*b2-nB[1]*b1)/det))
    else
        offs=Tuple{NTuple{2,Float64},Float64}[]
        for p in (pA,pB)
            p.kind===:line && continue
            er=_fillet_sc2(_fillet_sub2(p.P,p.O),1.0/p.rad)
            s=_fillet_dot2(p.nin,er)>=0 ? 1.0 : -1.0
            Rp=p.rad+s*ρ
            Rp>tol || throw(ArgumentError(
                "$caller: fillet radius $ρ swallows the profile circle " *
                "on Curve[$et]"))
            push!(offs,(p.O,Rp))
        end
        if pA.kind!==:circ || pB.kind!==:circ
            L=pA.kind===:line ? pA : pB
            Oc,Rp=only(offs)
            Q=_fillet_add2(rim,_fillet_sc2(L.nin,ρ))
            rel=_fillet_sub2(Q,Oc)
            b=_fillet_dot2(rel,L.d); c=_fillet_dot2(rel,rel)-Rp*Rp
            disc=b*b-c
            disc>=-tol || throw(ArgumentError(
                "$caller: Curve[$et] fillet radius $ρ does not reach " *
                "the profile circle"))
            for sgn in (1.0,-1.0)
                h=sqrt(max(0.0,disc))
                push!(cands,_fillet_add2(Q,_fillet_sc2(L.d,-b+sgn*h)))
            end
        else
            (O1,Rp1),(O2,Rp2)=offs
            dv=_fillet_sub2(O2,O1); d=_fillet_norm2(dv)
            d>1e-12 || throw(ArgumentError(
                "$caller: Curve[$et] has concentric rim profiles — no " *
                "isolated fillet centre"))
            a=(d*d+Rp1*Rp1-Rp2*Rp2)/(2d)
            h2=Rp1*Rp1-a*a
            h2>=-tol || throw(ArgumentError(
                "$caller: Curve[$et] profile circles cannot reach a " *
                "radius-$ρ fillet"))
            h=sqrt(max(0.0,h2))
            e=_fillet_sc2(dv,1.0/d)
            base=_fillet_add2(O1,_fillet_sc2(e,a))
            m=_fillet_perp(e)
            push!(cands,_fillet_add2(base,_fillet_sc2(m,h)),
                         _fillet_add2(base,_fillet_sc2(m,-h)))
        end
    end
    best=nothing; bestd=Inf
    for C in cands
        (_fillet_dot2(_fillet_sub2(C,rim),nA)>tol &&
         _fillet_dot2(_fillet_sub2(C,rim),nB)>tol) || continue
        Ts=NTuple{2,Float64}[]; legs=Float64[]; ok=true
        for (p,nother) in ((pA,nB),(pB,nA))
            local T,leg
            if p.kind===:line
                T=_fillet_add2(p.P,_fillet_sc2(p.d,
                    _fillet_dot2(_fillet_sub2(C,p.P),p.d)))
                leg=_fillet_dot2(_fillet_sub2(T,p.P),p.d)*
                    sign(_fillet_dot2(p.d,nother))
            else
                rel=_fillet_sub2(C,p.O); l=_fillet_norm2(rel)
                l>tol || (ok=false; break)
                T=_fillet_add2(p.O,_fillet_sc2(rel,p.rad/l))
                er=_fillet_sc2(_fillet_sub2(p.P,p.O),1.0/p.rad)
                s=sign(_fillet_dot2(_fillet_perp(er),nother))
                leg=atan(_fillet_cross(er,_fillet_sc2(rel,1.0/l)),
                         _fillet_dot2(er,_fillet_sc2(rel,1.0/l)))*s*p.rad
            end
            leg>tol || (ok=false; break)
            push!(Ts,T); push!(legs,leg)
        end
        ok || continue
        dd=_fillet_dot2(_fillet_sub2(C,rim),_fillet_sub2(C,rim))
        dd<bestd && (best=(C,Ts[1],Ts[2],legs[1],legs[2]); bestd=dd)
    end
    best===nothing && throw(ArgumentError(
        "$caller: Curve[$et] fillet radius $ρ does not fit the rim " *
        "notch — concave or over-running fillets are unsupported"))
    return best
end

# (domain-arclength s, inward distance δ) of meridian point `m` against
# profile `p`'s rim — `nother` is the companion face's inward normal
# (the domain side of the profile keeps inside the other's inward side).
function _occ_fillet_strip_coords(p,rim,nother,m)
    if p.kind===:line
        s=_fillet_dot2(_fillet_sub2(m,rim),p.d)*
            sign(_fillet_dot2(p.d,nother))
        return s,_fillet_dot2(_fillet_sub2(m,rim),p.nin)
    end
    er=_fillet_sc2(_fillet_sub2(rim,p.O),1.0/p.rad)
    rel=_fillet_sub2(m,p.O); l=_fillet_norm2(rel)
    em=l>0 ? _fillet_sc2(rel,1.0/l) : er
    s=sign(_fillet_dot2(_fillet_perp(er),nother))*
        atan(_fillet_cross(er,em),_fillet_dot2(er,em))*p.rad
    δ=sign(_fillet_dot2(p.nin,er))*(l-p.rad)
    return s,δ
end

# Strip-overrun certificate for one kept face: every other boundary of the
# face must stay strictly outside the strip the fillet removes — the
# domain-arclength interval (0, leg) from the rim to its tangent point.
# `rim_eid` is the edge being replaced (its own original position sits at
# s=0 by construction); `extra_ms` carries the tangent meridian points of
# other filleted rims on the same face so colliding fillets reject.
# Parametric seam edges (two uses inside the same face) are interior to
# the parametrization, not material boundaries — skipped. A boundary
# landing exactly on the tangent circle would leave a degenerate remnant.
function _occ_fillet_overrun(kernel,fid,rim_eid,p,rim,nother,leg,extra_ms,
                             axis,o,face_wires,et,caller)
    tol=1e-8*max(1.0,leg)
    counts=Dict{UInt64,Int}()
    for wid in face_wires[fid], (u,_) in kernel.shapes[wid].geometry.uses
        counts[u]=get(counts,u,0)+1
    end
    ftag=_model_occ_shape_tag(kernel,fid,caller)
    function check(m)
        s,δ=_occ_fillet_strip_coords(p,rim,nother,m)
        abs(s-leg)<=tol && throw(ArgumentError(
            "$caller: Curve[$et] fillet tangent lands on another " *
            "boundary of Surface[$ftag]"))
        (tol<s<leg-tol && δ>-tol) && throw(ArgumentError(
            "$caller: Curve[$et] fillet overruns a boundary of " *
            "Surface[$ftag]"))
    end
    for m in extra_ms
        check(m)
    end
    for wid in face_wires[fid], (u,_) in kernel.shapes[wid].geometry.uses
        u==rim_eid && continue
        get(counts,u,0)>=2 && continue
        g2=kernel.shapes[u].geometry.geometry
        coaxial=get(g2,:occ,nothing)===:circle &&
            abs(_dot3(g2.n,axis))>1-1e-9 && begin
            dc=_arc_sub(g2.center,o)
            _occ_modulus(_arc_sub(dc,_occ3_scale(axis,_dot3(dc,axis))))<
                1e-7*max(1.0,_occ_modulus(dc))
        end
        if coaxial
            check((g2.r,_dot3(_arc_sub(g2.center,o),axis)))
        else
            for i in 0:12
                t=g2.t0+(g2.t1-g2.t0)*i/12
                pt=_occ_fillet_edge_point(g2,t)
                d=_arc_sub(pt,o)
                check((_occ_modulus(_arc_sub(d,
                    _occ3_scale(axis,_dot3(d,axis)))),_dot3(d,axis)))
            end
        end
    end
    return nothing
end

# Plan one circle-edge fillet: rim profiles → meridian offset solve.
# Returns the spec the rebuild consumes; every fallible step precedes any
# kernel mutation (overrun certificates run in a separate pass once every
# filleted rim is known).
function _occ_fillet_circle_plan(kernel,eg,eid,fids,frev_of,vert_pos,
                                 et,r,caller)
    abs(eg.t1-eg.t0-_OCC_TWO_PI)<1e-6 || throw(ArgumentError(
        "$caller: fillet Curve[$et] must be a full 2π circle edge — " *
        "partial rims are unsupported"))
    axis=_occ_dir(eg.n); o=eg.center; R=eg.r
    fA,fB=fids
    sgA=kernel.shapes[fA].geometry.geometry
    sgB=kernel.shapes[fB].geometry.geometry
    # Outward normals through the shell-use sign alone — the measure path
    # (`_model_occ_get_mass`) certifies σ·axis/σ·(du×dv); the record's
    # `reversed` flag mirrors the same TopAbs orientation and must not be
    # multiplied in a second time.
    σA=frev_of[fA] ? -1.0 : 1.0
    σB=frev_of[fB] ? -1.0 : 1.0
    pA=_occ_fillet_rim_profile(sgA,σA,axis,o,R,et,caller)
    pB=_occ_fillet_rim_profile(sgB,σB,axis,o,R,et,caller)
    C2,TA,TB,legA,legB=_occ_fillet_meridian(pA,pB,R,r,et,caller)
    r_c,z_c=C2
    sc=max(1.0,R,r)
    r_c>1e-7*sc || throw(ArgumentError(
        "$caller: fillet on Curve[$et] has its tube centre on the " *
        "axis — spindle fillets are unsupported"))
    TA[1]>1e-7*sc && TB[1]>1e-7*sc || throw(ArgumentError(
        "$caller: fillet on Curve[$et] tangent collapses to a pole " *
        "on the axis — unsupported"))
    rimvid=kernel.shapes[eid].children[1]
    return (axis=axis,o=o,R=R,C2=C2,rho=r,rim=(R,0.0),
            rimvid=rimvid,rim_vpos=vert_pos[rimvid],
            fA=fA,fB=fB,sgA=sgA,sgB=sgB,pA=pA,pB=pB,
            TA=TA,TB=TB,legA=legA,legB=legB)
end

# Torus frame/v-range for one planned circle fillet: X_t is the seam
# azimuth (the rim vertex's radial direction); the v interval spans the
# two tangent meridians through the rim direction, canonicalized into the
# [0,2π] period when it fits.
function _occ_fillet_circle_frames(s)
    d=_arc_sub(s.rim_vpos,s.o)
    X_t=_occ_dir(_arc_sub(d,_occ3_scale(s.axis,_dot3(d,s.axis))))
    r_c,z_c=s.C2
    vof(T)=atan(T[2]-z_c,T[1]-r_c)
    vA,vB=vof(s.TA),vof(s.TB); vrim=atan(-z_c,s.R-r_c)
    al,ah=minmax(vA,vB)
    I=al<=vrim<=ah ? (al,ah) : (ah,al+_OCC_TWO_PI)
    mid=(I[1]+I[2])/2
    k=round((π-mid)/_OCC_TWO_PI)
    lo,hi=I[1]+k*_OCC_TWO_PI,I[2]+k*_OCC_TWO_PI
    function repv(T)
        v=vof(T)
        for q in (v,v+_OCC_TWO_PI,v-_OCC_TWO_PI)
            lo-1e-9<=q<=hi+1e-9 && return q
        end
        return v
    end
    vA,vB=repv(s.TA),repv(s.TB)
    vlo,vhi=minmax(vA,vB)
    # Spindle certificate: r_c + ρ·cos v must stay positive over the
    # fillet meridian — otherwise the torus tube self-intersects on the
    # axis and the patch folds through its own normal.
    cmin=min(cos(vlo),cos(vhi))
    for q in (π,π+_OCC_TWO_PI,π-_OCC_TWO_PI)
        vlo<q<vhi && (cmin=-1.0)
    end
    r_c+s.rho*cmin>1e-9*s.rho || throw(ArgumentError(
        "OCC.fillet: fillet tube radius swallows its own axis — " *
        "spindle fillets are unsupported"))
    return merge(s,(X_t=X_t,vlo=vlo,vhi=vhi,vA=vA,vB=vB))
end

# (u0, du, v_t) of a fresh tangent-circle edge `nc_eid` on face record `sg`:
# u0 lands on the [0,2π] period boundary for the seam azimuth, du is the
# azimuth direction sign, v_t the face's own v-coordinate at the tangent.
function _occ_fillet_face_uv(sg,nc_eid,kernel,caller)
    ng=kernel.shapes[nc_eid].geometry.geometry
    u0=atan(_dot3(ng.X,sg.Y),_dot3(ng.X,sg.X))
    du=sign(_dot3(ng.n,sg.axis))
    u0n=mod(u0,_OCC_TWO_PI)
    cen=ng.center; rt=ng.r
    vt=if sg.occ===:cone
        h=_dot3(_arc_sub(cen,sg.center),sg.axis)
        L=sqrt(sg.height*sg.height+(sg.r2-sg.r1)*(sg.r2-sg.r1))
        h*L/sg.height
    elseif sg.occ===:sphere
        atan(_dot3(_arc_sub(cen,sg.center),sg.axis),rt)
    elseif sg.occ===:torus
        atan(_dot3(_arc_sub(cen,sg.center),sg.axis),rt-sg.r1)
    else
        _dot3(_arc_sub(cen,sg.center),sg.axis)
    end
    return u0n,du,vt
end

# Publish every planned circle fillet: tangent circles replace the rims on
# the kept faces (wires/pcurves/vertex_uv re-bound on fresh face shapes,
# plane faces regenerated through `add_plane_surface`), meridian seam arcs
# and torus fillet faces bound directly. `new_face2` collects
# fid ⇒ (new fid, shell-rev) for the shell rewrite.
function _occ_fillet_circle_rebuild!(m,kernel,circ_specs,wires_of,frev_of,
                                    fnormals,new_face2,shell_uses,caller)
    isempty(circ_specs) && return nothing
    rim_repl=Dict{Tuple{UInt64,UInt64},NamedTuple}()
    vmap=Dict{Tuple{UInt64,UInt64},UInt64}()
    for (eid,s0) in circ_specs
        s=_occ_fillet_circle_frames(s0)
        circ_specs[eid]=s
        for (fid,T,sg) in ((s.fA,s.TA,s.sgA),(s.fB,s.TB,s.sgB))
            cen3=_occ3_add(s.o,_occ3_scale(s.axis,T[2]))
            nct=_model_occ_add_circle!(m,cen3[1],cen3[2],cen3[3],T[1];
                zaxis=collect(s.axis),
                xaxis=collect(_arc_sub(s.rim_vpos,cen3)))
            nc_eid=kernel.bindings[(1,nct)]
            nc_vid=kernel.shapes[nc_eid].children[1]
            u0,du,vt=sg.occ===:plane ? (0.0,1.0,0.0) :
                _occ_fillet_face_uv(sg,nc_eid,kernel,caller)
            rim_repl[(fid,eid)]=(eid_new=nc_eid,etag=nct,vid=nc_vid,
                vtag=_model_occ_shape_tag(kernel,nc_vid,caller),
                u0=u0,du=du,vt=vt)
            vmap[(fid,s.rimvid)]=nc_vid
        end
    end
    # Kept faces: replace rim edges in their wires, retarget edges that
    # shared the rim vertex (meridian seams), swap rim pcurve/vertex_uv
    # entries, then rebind each face under a fresh tag.
    dirty=Set{UInt64}(fid for (fid,_) in keys(rim_repl))
    for fid in dirty
        fgsh=kernel.shapes[fid]; fg=fgsh.geometry.geometry
        ispl=fg.occ===:plane
        ftag=_model_occ_shape_tag(kernel,fid,caller)
        if !ispl
            hasproperty(fg,:pcurves) || throw(ArgumentError(
                "$caller: Surface[$ftag] carries no pcurves"))
        end
        pco=ispl ? nothing : copy(fg.pcurves)
        vuvo=hasproperty(fg,:vertex_uv) ? copy(fg.vertex_uv) : nothing
        new_wtags=Int[]; new_wids=UInt64[]
        for wid in wires_of[fid]
            uses=kernel.shapes[wid].geometry.uses
            changed=false; nuses=Tuple{UInt64,Bool}[]
            for (u,rev) in uses
                repl=get(rim_repl,(fid,u),nothing)
                if repl!==nothing
                    push!(nuses,(repl.eid_new,rev)); changed=true
                    continue
                end
                ch=kernel.shapes[u].children
                c1=get(vmap,(fid,ch[1]),ch[1])
                c2=get(vmap,(fid,ch[2]),ch[2])
                if c1==ch[1] && c2==ch[2]
                    push!(nuses,(u,rev)); continue
                end
                changed=true
                g2=kernel.shapes[u].geometry.geometry
                tA=_occ_2df_param(g2,kernel.shapes[c1].geometry.point)
                tB=_occ_2df_param(g2,kernel.shapes[c2].geometry.point)
                tA<tB-1e-12 || throw(ArgumentError(
                    "$caller: fillet truncation inverts an edge of " *
                    "Surface[$ftag]"))
                rec=merge(g2,(t0=tA,t1=tB))
                hasproperty(rec,:native_geometry) && (rec=merge(rec,
                    (native_geometry=merge(rec.native_geometry,
                        (t0=tA,t1=tB)),)))
                neid=only(_model_occ_reserve_shapes(kernel,1,caller))
                netag=_model_occ_next_tag(kernel,1,-1,caller)
                _model_occ_bind_shape!(kernel,1,netag,neid,
                    _ModelOCCShape(1,(kind=kernel.shapes[u].geometry.kind,
                                    geometry=rec),UInt64[c1,c2]))
                kernel.next_shape=neid+UInt64(1)
                push!(nuses,(neid,rev))
                if pco!==nothing
                    oldt=_model_occ_shape_tag(kernel,u,caller)
                    haskey(pco,oldt) && (pco[netag]=pop!(pco,oldt))
                end
            end
            if changed
                wt=_occ_wire_from_uses!(m,nuses,caller)
                push!(new_wtags,wt); push!(new_wids,kernel.bindings[(-1,wt)])
            else
                push!(new_wtags,_model_occ_shape_tag(kernel,wid,caller))
                push!(new_wids,wid)
            end
        end
        if pco!==nothing
            for ((ffid,ueid),repl) in rim_repl
                ffid===fid || continue
                delete!(pco,_model_occ_shape_tag(kernel,ueid,caller))
                pco[repl.etag]=(fwd=(lin2d=((repl.u0,repl.vt),
                                            (repl.du,0.0)),),rev=nothing)
            end
        end
        if vuvo!==nothing
            for ((ffid,ueid),repl) in rim_repl
                ffid===fid || continue
                ovt=_model_occ_shape_tag(kernel,
                    kernel.shapes[ueid].children[1],caller)
                haskey(vuvo,ovt) && delete!(vuvo,ovt)
                vuvo[repl.vtag]=(repl.u0,repl.vt)
            end
        end
        if ispl
            nft=_model_occ_add_plane_surface!(m,new_wtags)
            nfid=kernel.bindings[(2,nft)]
            nga=kernel.shapes[nfid].geometry.geometry.axis
            # Outward use-sign σ·axis (rev only) — same convention the
            # measure path certifies.
            σo=frev_of[fid] ? -1.0 : 1.0
            new_face2[fid]=(nfid,_dot3(nga,fg.axis)*σo<0)
        else
            newgeo=merge(fg,(pcurves=pco,))
            vuvo!==nothing && (newgeo=merge(newgeo,(vertex_uv=vuvo,)))
            nft=_model_occ_next_tag(kernel,2,-1,caller)
            nfid=only(_model_occ_reserve_shapes(kernel,1,caller))
            _model_occ_bind_shape!(kernel,2,nft,nfid,
                _ModelOCCShape(2,(kind=fgsh.geometry.kind,geometry=newgeo),
                               new_wids))
            kernel.next_shape=nfid+UInt64(1)
            new_face2[fid]=(nfid,frev_of[fid])
        end
    end
    # Torus fillet faces — one per circle edge, 4-edge periodic wire like
    # the canonical cylinder lateral.
    for (eid,s) in circ_specs
        a=s.axis; X_t=s.X_t; ρ=s.rho; r_c,z_c=s.C2
        flo,fhi=s.vA<=s.vB ? (s.fA,s.fB) : (s.fB,s.fA)
        rlo=rim_repl[(flo,eid)]; rhi=rim_repl[(fhi,eid)]
        Y_t=_occ_cross(a,X_t)
        C3=_occ3_add(s.o,_occ3_scale(a,z_c))
        Cm=_occ3_add(C3,_occ3_scale(X_t,r_c))
        seid=only(_model_occ_reserve_shapes(kernel,1,caller))
        setag=_model_occ_next_tag(kernel,1,-1,caller)
        _model_occ_bind_shape!(kernel,1,setag,seid,
            _ModelOCCShape(1,(kind=:circle,
                geometry=(occ=:circle,center=Cm,n=_occ_cross(X_t,a),
                          X=X_t,Y=a,r=ρ,t0=s.vlo,t1=s.vhi)),
                UInt64[rlo.vid,rhi.vid]))
        kernel.next_shape=seid+UInt64(1)
        pc=Dict{Int,NamedTuple}(
            rlo.etag=>(fwd=(lin2d=((0.0,s.vlo),(1.0,0.0)),),rev=nothing),
            rhi.etag=>(fwd=(lin2d=((0.0,s.vhi),(1.0,0.0)),),rev=nothing),
            setag=>(fwd=(lin2d=((_OCC_TWO_PI,0.0),(0.0,1.0)),),
                    rev=(lin2d=((0.0,0.0),(0.0,1.0)),)))
        wuses=Tuple{UInt64,Bool}[(rhi.eid_new,true),(seid,true),
                                 (rlo.eid_new,false),(seid,false)]
        wt=_occ_wire_from_uses!(m,wuses,caller)
        wid2=kernel.bindings[(-1,wt)]
        nft=_model_occ_next_tag(kernel,2,-1,caller)
        nfid=only(_model_occ_reserve_shapes(kernel,1,caller))
        _model_occ_bind_shape!(kernel,2,nft,nfid,
            _ModelOCCShape(2,(kind=:torus,
                geometry=(occ=:torus,center=C3,axis=a,X=X_t,Y=Y_t,
                          r1=r_c,r2=ρ,angle=_OCC_TWO_PI,pcurves=pc)),
                UInt64[wid2]))
        kernel.next_shape=nfid+UInt64(1)
        push!(shell_uses,(nfid,false))
    end
    return nothing
end

# ── circle-edge fillets on general faces (ChFi3d pipe blends) ────────────
#
# A full-circle rim that fails the coaxial-revolution profile test — a
# `tensor_bspline` loft wall, a trimmed cap patch, an off-axis analytic
# face — still admits the rolling-ball fillet whenever the ball centre
# stays azimuth-monotone: at rim azimuth θ the meridian section carries
# one tube centre C(θ) = o + r̂(θ)·r_b + axis·z_b solving, per incident
# face k ∈ {A,B},
#   (C − S_k(u_k,v_k))·S_u = 0,   (C − S_k(u_k,v_k))·S_v = 0,
#   |C − S_k| = ρ,
# the ChFi3d walking-section equations. The fillet face is the pipe
# surface between the two tangent curves, stored as a rational
# `tensor_bspline` patch: periodic cubic in the rim azimuth, one
# degree-2 rational span across the meridian arc — the station
# iso-arcs are exact circle arcs, the tangent boundary curves are
# interpolated exactly (they are the shared edges), and the interior
# approximates the true canal surface like OCC's own B-spline fillet.

# Number of azimuth stations for the swept fillet net (retried doubled
# when the mid-station accuracy certificate fails).
const _OCC_FILLET_SWEEP_N = 24

# Periodic cubic B-spline interpolation: solve the circulant system
# (Q_{i−1} + 4Q_i + Q_{i+1})/6 = D_i for homogeneous poles `data`
# sampled at uniformly spaced knots over one period.
function _fillet_periodic_poles(data::Vector{NTuple{4,Float64}},caller)
    n=length(data)
    n>=4 || throw(ArgumentError(
        "$caller: too few stations for the fillet sweep"))
    A=zeros(Float64,n,n)
    for i in 1:n
        A[i,i]=4.0/6
        A[i,i==1 ? n : i-1]=1.0/6
        A[i,i==n ? 1 : i+1]=1.0/6
    end
    F=try
        LinearAlgebra.lu(A)
    catch err
        err isa LinearAlgebra.SingularException || rethrow()
        throw(ArgumentError(
            "$caller: fillet sweep interpolation is singular"))
    end
    R=Matrix{Float64}(undef,n,4)
    for i in 1:n
        d=data[i]
        R[i,1]=d[1]; R[i,2]=d[2]; R[i,3]=d[3]; R[i,4]=d[4]
    end
    X=F\R
    all(isfinite,X) || throw(ArgumentError(
        "$caller: fillet sweep interpolation is singular"))
    return NTuple{4,Float64}[(X[i,1],X[i,2],X[i,3],X[i,4]) for i in 1:n]
end

# Non-clamped uniform knot vector reproducing a degree-3 periodic
# B-spline on [0,2π] with `n` distinct poles wrapped 3 deep: knots
# (j−3)·2π/n for j = 0..n+6, active domain [0,2π].
_fillet_periodic_knots(n::Int)=[(j-3)*(2π/n) for j in 0:n+6]

# Rim pcurve evaluation on the edge's own periodic parameter: the rim
# is a closed 2π circle, so the uv trace is defined mod the edge span
# (analytic pcurves extrapolate harmlessly, spl2d ones must stay in
# range).
_rim_pc_uv(pc,eg,θ)=_occ_pcurve_d0(pc,
    eg.t0+mod(θ-eg.t0,eg.t1-eg.t0))

# Incident-face descriptor for one swept fillet. `perdim`/`period` name
# the periodic parameter direction (1=u, 2=v, 0=none); `seam_uv` is the
# iso value of the face's parametric seam edge in that direction when a
# doubled wire use exists; `seam_raw` keeps the stored (fwd-slot) iso
# value for vertex_uv bookkeeping.
function _occ_fillet_sweep_face(kernel,sg,fid,σ,et,wires_of,caller)
    k=get(sg,:occ,nothing)
    (_model_nurbs_surface(sg) ||
     k in (:plane,:cylinder,:cone,:sphere,:torus)) || throw(ArgumentError(
        "$caller: Curve[$et] is a circle edge shared by an occ=$k " *
        "face — only coaxial plane/cylinder/cone/sphere/torus rims " *
        "and faces carrying surface evaluators are supported"))
    hasproperty(sg,:pcurves) || throw(ArgumentError(
        "$caller: Curve[$et] rim face has no pcurves — surface " *
        "blends are unsupported"))
    ent=get(sg.pcurves,et,nothing)
    ent===nothing && throw(ArgumentError(
        "$caller: Curve[$et] carries no pcurve on its rim face — " *
        "surface blends are unsupported"))
    pc=ent.fwd===nothing ? ent.rev : ent.fwd
    lo,hi=_occ_surface_bounds(sg)
    perdim=0; period=0.0
    if _model_nurbs_surface(sg)
        for (d,p) in ((1,_model_nurbs_patch_seam_period(sg,1)),
                      (2,_model_nurbs_patch_seam_period(sg,2)))
            p>0.0 && (perdim=d; period=p)
        end
    elseif k in (:cylinder,:cone,:sphere,:torus)
        perdim=1; period=get(sg,:angle,_OCC_TWO_PI)
    end
    # Parametric seams: edge records used twice inside this face's
    # wires bound the periodic parametrization, not the material.
    counts=Dict{UInt64,Int}()
    for wid in wires_of[fid], (u,_) in kernel.shapes[wid].geometry.uses
        counts[u]=get(counts,u,0)+1
    end
    seam_uv=NaN; seam_raw=NaN
    if perdim!=0
        for wid in wires_of[fid], (u,_) in kernel.shapes[wid].geometry.uses
            get(counts,u,0)>=2 || continue
            ug=kernel.shapes[u].geometry.geometry
            uetag=_model_occ_shape_tag(kernel,u,caller)
            pent=get(sg.pcurves,uetag,nothing)
            pent===nothing && continue
            ppc=pent.fwd===nothing ? pent.rev : pent.fwd
            ppc===nothing && continue
            uv=_occ_pcurve_d0(ppc,ug.t0)
            seam_raw=perdim==1 ? uv[1] : uv[2]
            seam_uv=seam_raw - period*floor(seam_raw/period)
            seam_uv>period-1e-12 && (seam_uv=0.0)
            break
        end
    end
    return (sg=sg,fid=fid,σ=σ,pc=pc,lo=lo,hi=hi,
            ctr=((lo[1]+hi[1])/2,(lo[2]+hi[2])/2),
            perdim=perdim,period=period,seam_uv=seam_uv,
            seam_raw=seam_raw)
end

# Meridional inward-normal projection at rim azimuth θ on face `info`:
# the 3-D inward normal at the rim point projected ⊥ rim tangent —
# seeds the 2-D profile solve and certifies the tangent leg direction.
function _occ_fillet_sweep_frame(info,uv,th,rh,axis)
    _,Su,Sv=_occ_surface_d1(info.sg,uv[1],uv[2])
    fn=_occ_cross(Su,Sv)
    fl=_occ_modulus(fn)
    fl>1e-12 || return nothing
    fn=_occ3_scale(fn,1.0/fl)
    n_in=_occ3_scale(fn,-info.σ)
    nin3=_arc_sub(n_in,_occ3_scale(th,_dot3(n_in,th)))
    l=_occ_modulus(nin3)
    l>1e-9 || return nothing
    nin3=_occ3_scale(nin3,1.0/l)
    d3=_occ_dir(_occ_cross(th,nin3))
    iud=_fillet_sub2(info.ctr,uv)
    i3=_occ3_add(_occ3_scale(Su,iud[1]),_occ3_scale(Sv,iud[2]))
    _dot3(d3,i3)<0 && (d3=(.-d3))
    return nin3,d3
end

# One rolling-ball section at rim azimuth θ: Newton on
# x = (r_b−R, z_b, u_A, v_A, u_B, v_B). `x0` is the previous station's
# solution (continuation) or the profile seed.
function _occ_fillet_sweep_station(o,axis,rh,R,ρ,infoA,infoB,x,et,caller)
    s,z,uA,vA,uB,vB=x
    ftol=1e-11*max(1.0,ρ)
    for it in 1:80
        P1,Su1,Sv1,Suu1,Svv1,Suv1=_occ_surface_d2(infoA.sg,uA,vA)
        P2,Su2,Sv2,Suu2,Svv2,Suv2=_occ_surface_d2(infoB.sg,uB,vB)
        C=_occ3_add(_occ3_add(o,_occ3_scale(rh,R+s)),
                    _occ3_scale(axis,z))
        rv1=_arc_sub(C,P1); rv2=_arc_sub(C,P2)
        f1=_dot3(rv1,Su1); f2=_dot3(rv1,Sv1)
        f3=_dot3(rv1,rv1)-ρ*ρ
        f4=_dot3(rv2,Su2); f5=_dot3(rv2,Sv2)
        f6=_dot3(rv2,rv2)-ρ*ρ
        max(abs(f1),abs(f2),abs(f3),abs(f4),abs(f5),abs(f6))<ftol &&
            return (s,z,uA,vA,uB,vB)
        J=Float64[
            _dot3(rh,Su1)  _dot3(axis,Su1)  (-_dot3(Su1,Su1)+_dot3(rv1,Suu1))  (-_dot3(Sv1,Su1)+_dot3(rv1,Suv1))  0.0  0.0;
            _dot3(rh,Sv1)  _dot3(axis,Sv1)  (-_dot3(Su1,Sv1)+_dot3(rv1,Suv1))  (-_dot3(Sv1,Sv1)+_dot3(rv1,Svv1))  0.0  0.0;
            2*_dot3(rv1,rh)  2*_dot3(rv1,axis)  (-2*_dot3(rv1,Su1))  (-2*_dot3(rv1,Sv1))  0.0  0.0;
            _dot3(rh,Su2)  _dot3(axis,Su2)  0.0  0.0  (-_dot3(Su2,Su2)+_dot3(rv2,Suu2))  (-_dot3(Sv2,Su2)+_dot3(rv2,Suv2));
            _dot3(rh,Sv2)  _dot3(axis,Sv2)  0.0  0.0  (-_dot3(Su2,Sv2)+_dot3(rv2,Suv2))  (-_dot3(Sv2,Sv2)+_dot3(rv2,Svv2));
            2*_dot3(rv2,rh)  2*_dot3(rv2,axis)  0.0  0.0  (-2*_dot3(rv2,Su2))  (-2*_dot3(rv2,Sv2))]
        δ=try
            J\-[f1,f2,f3,f4,f5,f6]
        catch err
            err isa LinearAlgebra.SingularException || rethrow()
            throw(ArgumentError(
                "$caller: Curve[$et] fillet section solve is singular"))
        end
        # Damped step: retry with halved increments while the combined
        # residual grows.
        scale=1.0
        accepted=false
        for _ in 1:12
            s2=s+scale*δ[1]; z2=z+scale*δ[2]
            uA2=uA+scale*δ[3]; vA2=vA+scale*δ[4]
            uB2=uB+scale*δ[5]; vB2=vB+scale*δ[6]
            uA2=_fillet_sweep_wrapdim(infoA,1,uA2)
            vA2=_fillet_sweep_wrapdim(infoA,2,vA2)
            uB2=_fillet_sweep_wrapdim(infoB,1,uB2)
            vB2=_fillet_sweep_wrapdim(infoB,2,vB2)
            Q1,Qu1,Qv1=_occ_surface_d1(infoA.sg,uA2,vA2)
            Q2,Qw1,Qw2=_occ_surface_d1(infoB.sg,uB2,vB2)
            C2=_occ3_add(_occ3_add(o,_occ3_scale(rh,R+s2)),
                         _occ3_scale(axis,z2))
            rw1=_arc_sub(C2,Q1); rw2=_arc_sub(C2,Q2)
            res=max(abs(_dot3(rw1,Qu1)),abs(_dot3(rw1,Qv1)),
                    abs(_dot3(rw1,rw1)-ρ*ρ),
                    abs(_dot3(rw2,Qw1)),abs(_dot3(rw2,Qw2)),
                    abs(_dot3(rw2,rw2)-ρ*ρ))
            if res<max(abs(f1),abs(f2),abs(f3),abs(f4),abs(f5),abs(f6)) ||
               scale<1e-4
                s,z,uA,vA,uB,vB=s2,z2,uA2,vA2,uB2,vB2
                accepted=true
                break
            end
            scale*=0.5
        end
        accepted || throw(ArgumentError(
            "$caller: Curve[$et] fillet section solve diverged — " *
            "no rolling-ball fit"))
    end
    throw(ArgumentError(
        "$caller: Curve[$et] fillet section solve did not converge — " *
        "no rolling-ball fit"))
end

# Wrap/clamp a surface parameter inside the face's active domain:
# periodic directions wrap by whole periods, others clamp to the
# stored bounds padded by half a span so seeds can recover.
function _fillet_sweep_wrapdim(info,dim,x)
    lo,hi=info.lo[dim],info.hi[dim]
    if info.perdim==dim
        return lo+mod(x-lo,info.period)
    end
    span=max(hi-lo,1e-12)
    return clamp(x,lo-0.5*span,hi+0.5*span)
end

# Profile-seed for the station solve: the local 2-D line-profile
# meridian solve plus a linear surface projection of each tangent.
function _occ_fillet_sweep_seed(o,axis,X,Y,eg,R,ρ,infoA,infoB,θ,et,
                                caller)
    rh=_occ3_add(_occ3_scale(X,cos(θ)),_occ3_scale(Y,sin(θ)))
    th=_occ3_add(_occ3_scale(X,-sin(θ)),_occ3_scale(Y,cos(θ)))
    profs=Any[]; uvs=NTuple{2,Float64}[]; tans=Any[]
    for info in (infoA,infoB)
        uv=_rim_pc_uv(info.pc,eg,θ)
        fr=_occ_fillet_sweep_frame(info,uv,th,rh,axis)
        fr===nothing && throw(ArgumentError(
            "$caller: Curve[$et] rim face is tangent-degenerate — " *
            "no fillet section"))
        nin3,d3=fr
        nin2=(_dot3(nin3,rh),_dot3(nin3,axis))
        d2=(_dot3(d3,rh),_dot3(d3,axis))
        push!(profs,(kind=:line,P=(R,0.0),
            d=_fillet_sc2(d2,1.0/_fillet_norm2(d2)),
            nin=_fillet_sc2(nin2,1.0/_fillet_norm2(nin2))))
        Pr,Su,Sv=_occ_surface_d1(info.sg,uv[1],uv[2])
        push!(uvs,uv); push!(tans,(Pr,Su,Sv))
    end
    C2,TA2,TB2,_,_=_occ_fillet_meridian(profs[1],profs[2],R,ρ,et,caller)
    out=[C2[1]-R,C2[2]]
    for k in 1:2
        uv=uvs[k]; T_2=(TA2,TB2)[k]; (Pr,Su,Sv)=tans[k]
        T3=_occ3_add(_occ3_add(o,_occ3_scale(rh,T_2[1])),
                     _occ3_scale(axis,T_2[2]))
        rel=_arc_sub(T3,Pr)
        a11=_dot3(Su,Su); a12=_dot3(Su,Sv); a22=_dot3(Sv,Sv)
        det=a11*a22-a12*a12
        du=dv=0.0
        if abs(det)>1e-24
            du=(_dot3(rel,Su)*a22-_dot3(rel,Sv)*a12)/det
            dv=(_dot3(rel,Sv)*a11-_dot3(rel,Su)*a12)/det
        end
        push!(out,uv[1]+du,uv[2]+dv)
    end
    return (out[1],out[2],out[3],out[4],out[5],out[6])
end

# Station record assembled per azimuth: ball centre, tangent points,
# tangent uv, meridian arc angle, and the tangent leg in uv.
function _occ_fillet_sweep_stations(kernel,o,axis,X,Y,eg,R,ρ,infoA,
                                  infoB,phis,et,caller)
    n=length(phis)
    out=Vector{NamedTuple}(undef,n)
    x=nothing; uwa=nothing; uwb=nothing
    for (i,θ) in enumerate(phis)
        rh=_occ3_add(_occ3_scale(X,cos(θ)),_occ3_scale(Y,sin(θ)))
        th=_occ3_add(_occ3_scale(X,-sin(θ)),_occ3_scale(Y,cos(θ)))
        x===nothing &&
            (x=_occ_fillet_sweep_seed(o,axis,X,Y,eg,R,ρ,infoA,infoB,θ,
                                      et,caller))
        s,z,uA,vA,uB,vB=_occ_fillet_sweep_station(o,axis,rh,R,ρ,
                                                infoA,infoB,x,et,
                                                caller)
        # Unwrap the periodic tangent coordinates along the sweep so the
        # series stay continuous across the face seam.
        for (info,u,v,sym) in ((infoA,uA,vA,:A),(infoB,uB,vB,:B))
            info.perdim==0 && continue
            raw=info.perdim==1 ? u : v
            prev=sym===:A ? uwa : uwb
            unw=prev===nothing ? raw :
                raw+info.period*round((prev-raw)/info.period)
            sym===:A ? (uwa=unw) : (uwb=unw)
        end
        C=_occ3_add(_occ3_add(o,_occ3_scale(rh,R+s)),
                    _occ3_scale(axis,z))
        TA=_occ_surface_point(infoA.sg,uA,vA)
        TB=_occ_surface_point(infoB.sg,uB,vB)
        uvA=uwa===nothing ? (uA,vA) :
            (infoA.perdim==1 ? (uwa,vA) : (uA,uwa))
        uvB=uwb===nothing ? (uB,vB) :
            (infoB.perdim==1 ? (uwb,vB) : (uB,uwb))
        # Validity: ball on the material-inward side of each face,
        # tangent inside the domain on the rim's interior side, tube
        # centre off the axis, and a nondegenerate meridian arc.
        for (info,P,u,v) in ((infoA,TA,uA,vA),(infoB,TB,uB,vB))
            _,Su,Sv=_occ_surface_d1(info.sg,u,v)
            n_in=_occ3_scale(_occ_cross(Su,Sv),-info.σ)
            _dot3(_arc_sub(C,P),n_in)>0.0 || throw(ArgumentError(
                "$caller: Curve[$et] fillet centre sits outside the " *
                "material notch — concave fillets are unsupported"))
        end
        for (info,uvT,rawu,rawv) in
                ((infoA,uvA,uA,vA),(infoB,uvB,uB,vB))
            uvR=_rim_pc_uv(info.pc,eg,θ)
            info.perdim==1 && (uvR=(uvT[1],uvR[2]))
            info.perdim==2 && (uvR=(uvR[1],uvT[2]))
            iud=_fillet_sub2(info.ctr,uvR)
            leg=_fillet_dot2(_fillet_sub2(uvT,uvR),iud)
            leg>1e-9*max(1.0,_fillet_norm2(iud)) || throw(ArgumentError(
                "$caller: Curve[$et] fillet tangent leaves the face " *
                "domain — no rolling-ball fit"))
            (info.lo[1]-1e-9<=rawu<=info.hi[1]+1e-9 &&
             info.lo[2]-1e-9<=rawv<=info.hi[2]+1e-9) || throw(
                ArgumentError(
                    "$caller: Curve[$et] fillet tangent overruns the " *
                    "face domain"))
        end
        sc=max(1.0,R,ρ)
        R+s>1e-7*sc || throw(ArgumentError(
            "$caller: fillet on Curve[$et] has its tube centre on " *
            "the axis — spindle fillets are unsupported"))
        e1=_arc_sub(TA,C); e2=_arc_sub(TB,C)
        dφ=acos(clamp(_dot3(e1,e2)/(ρ*ρ),-1.0,1.0))
        (dφ>0.05 && dφ<π-0.05) || throw(ArgumentError(
            "$caller: fillet on Curve[$et] has a degenerate meridian " *
            "arc — unsupported"))
        out[i]=(θ=θ,s=s,z=z,C=C,TA=TA,TB=TB,uvA=uvA,uvB=uvB,
                uA=uA,vA=vA,uB=uB,vB=vB,dφ=dφ)
        x=(s,z,uA,vA,uB,vB)
    end
    return out
end


# Strip-overrun certificate for one swept face: the uv strip between
# the rim pcurve and the tangent pcurve must contain no other material
# boundary, and no other filleted rim's tangent may enter it.
function _occ_fillet_sweep_overrun(kernel,fid,eid,spec,others,wires_of,
                                 et,caller)
    info=fid==spec.infoA.fid ? spec.infoA : spec.infoB
    sg=info.sg
    rim_uv=spec.rim_uv[fid]; tan_uv=spec.tan_uv[fid]
    # Strip polygon in unwrapped uv: rim trace forward, tangent back.
    poly=NTuple{2,Float64}[]
    append!(poly,rim_uv)
    for i in length(tan_uv):-1:1
        push!(poly,tan_uv[i])
    end
    ftag=_model_occ_shape_tag(kernel,fid,caller)
    counts=Dict{UInt64,Int}()
    for wid in wires_of[fid], (u,_) in kernel.shapes[wid].geometry.uses
        counts[u]=get(counts,u,0)+1
    end
    function inside(q)
        # Even-odd ray cast at v→+∞ against the closed polyline.
        c=false
        n=length(poly)
        for i in 1:n
            a=poly[i]; b=poly[i==n ? 1 : i+1]
            ((a[2]>q[2]) != (b[2]>q[2])) || continue
            xin=a[1]+(q[2]-a[2])*(b[1]-a[1])/(b[2]-a[2])
            xin>q[1] && (c=!c)
        end
        return c
    end
    function stripd(q)
        d=Inf
        for i in 1:length(tan_uv)-1
            a=tan_uv[i]; b=tan_uv[i+1]
            e=_fillet_sub2(b,a); t=_fillet_dot2(_fillet_sub2(q,a),e)
            L=_fillet_dot2(e,e)
            t=L>0 ? clamp(t/L,0.0,1.0) : 0.0
            p=_fillet_add2(a,_fillet_sc2(e,t))
            d=min(d,_fillet_norm2(_fillet_sub2(q,p)))
        end
        return d
    end
    function check(q,what)
        qq=q
        if info.perdim!=0
            uv=info.perdim==1 ? qq[1] : qq[2]
            ref=info.perdim==1 ? rim_uv[1][1] : rim_uv[1][2]
            uv+=info.period*round((ref-uv)/info.period)
            qq=info.perdim==1 ? (uv,qq[2]) : (qq[1],uv)
        end
        inside(qq) && throw(ArgumentError(
            "$caller: Curve[$et] fillet overruns $what of " *
            "Surface[$ftag]"))
        stripd(qq)<=1e-9*max(1.0,spec.R) && throw(ArgumentError(
            "$caller: Curve[$et] fillet tangent lands on $what of " *
            "Surface[$ftag]"))
    end
    for wid in wires_of[fid], (u,_) in kernel.shapes[wid].geometry.uses
        u==eid && continue
        get(counts,u,0)>=2 && continue
        uetag=_model_occ_shape_tag(kernel,u,caller)
        pent=get(sg.pcurves,uetag,nothing)
        pent===nothing && continue
        ppc=pent.fwd===nothing ? pent.rev : pent.fwd
        ppc===nothing && continue
        ug=kernel.shapes[u].geometry.geometry
        for i in 0:16
            t=ug.t0+(ug.t1-ug.t0)*i/16
            check(_occ_pcurve_d0(ppc,t),"a boundary")
        end
    end
    for (eid2,s2) in others
        eid2===eid && continue
        for (ofid,tan2) in ((s2.infoA.fid,s2.tan_uv[s2.infoA.fid]),
                          (s2.infoB.fid,s2.tan_uv[s2.infoB.fid]))
            ofid==fid || continue
            for q in tan2
                check(q,"another fillet tangent")
            end
        end
    end
    return nothing
end

# Build one closed tangent edge from the unwrapped periodic stations:
# degree-3 periodic interpolation wrapped into the non-clamped uniform
# knot convention, parameter domain [0,2π] starting at the seam
# crossing `φ_s`.
function _occ_fillet_sweep_edge(points)
    n=length(points)
    poles=_fillet_periodic_poles(
        NTuple{4,Float64}[(p[1],p[2],p[3],1.0) for p in points],
        "OCC.fillet")
    controls=NTuple{3,Float64}[]
    # Same abscissa convention as the fillet face: controls[k] is
    # interpolated at t=(k−2)·h, so rotate the net one pole to put
    # station 1 at t=0 (the minted seam vertex).
    for i in 0:n+2
        q=poles[mod(i-1,n)+1]
        w=q[4]
        push!(controls,(q[1]/w,q[2]/w,q[3]/w))
    end
    curve=_ModelOCCNURBSCurve(:bspline,NURBS.NURBSCurve(
        3,_fillet_periodic_knots(n),controls,
        ones(Float64,length(controls))))
    return (occ=:copied_curve,curve=curve,t0=0.0,t1=2π)
end

# Plan one swept circle fillet: stations → tangent curves → the seam
# crossing and the raw rebuild data.  Everything stays a pure
# computation until the caller's rebuild phase mutates the kernel.
function _occ_fillet_sweep_plan(kernel,eg,eid,fids,frev_of,vert_pos,
                                wires_of,et,r,caller)
    abs(eg.t1-eg.t0-_OCC_TWO_PI)<1e-6 || throw(ArgumentError(
        "$caller: fillet Curve[$et] must be a full 2π circle edge — " *
        "partial rims are unsupported"))
    axis=_occ_dir(eg.n); o=eg.center; R=eg.r
    X=eg.X; Y=_occ_circle_y(eg)
    fA,fB=fids
    σA=frev_of[fA] ? -1.0 : 1.0
    σB=frev_of[fB] ? -1.0 : 1.0
    infoA=_occ_fillet_sweep_face(kernel,kernel.shapes[fA].geometry.geometry,
                                 fA,σA,et,wires_of,caller)
    infoB=_occ_fillet_sweep_face(kernel,kernel.shapes[fB].geometry.geometry,
                                 fB,σB,et,wires_of,caller)
    # Detection pass: stations over the rim period, seam-crossing
    # azimuths per face, and the unwrapped tangent traces.
    N=_OCC_FILLET_SWEEP_N
    phis=ntuple(i->(i-1)*_OCC_TWO_PI/N,N)
    sts=_occ_fillet_sweep_stations(kernel,o,axis,X,Y,eg,R,r,
                                   infoA,infoB,phis,et,caller)
    tan_uv=Dict{UInt64,Vector{NTuple{2,Float64}}}()
    rim_uv=Dict{UInt64,Vector{NTuple{2,Float64}}}()
    for info in (infoA,infoB)
        ser=NTuple{2,Float64}[]
        rser=NTuple{2,Float64}[]
        for (i,st) in enumerate(sts)
            uv=info.fid==fA ? st.uvA : st.uvB
            push!(ser,uv)
            ruv=_rim_pc_uv(info.pc,eg,st.θ)
            # Unwrap the rim's periodic coordinate into the tangent's
            # frame so the strip polygon is well-defined.
            if info.perdim==1
                uu=ruv[1]+info.period*round((uv[1]-ruv[1])/info.period)
                ruv=(uu,ruv[2])
            elseif info.perdim==2
                vv=ruv[2]+info.period*round((uv[2]-ruv[2])/info.period)
                ruv=(ruv[1],vv)
            end
            push!(rser,ruv)
        end
        tan_uv[info.fid]=ser; rim_uv[info.fid]=rser
        # Winding certificate: on a periodic face the unwrapped trace
        # must wind monotonically and close onto itself after one
        # period — the stations cover N−1 of N hops, so the check adds
        # the virtual closing hop back to the first sample one period
        # along.
        if info.perdim!=0
            d=info.perdim
            vals=[q[d] for q in ser]
            span=vals[end]-vals[1]
            sgn=sign(span)
            sgn==0.0 && throw(ArgumentError(
                "$caller: Curve[$et] fillet tangent does not wrap " *
                "the periodic face — unsupported"))
            for i in 2:length(vals)
                (vals[i]-vals[i-1])*sgn>-1e-9*info.period || throw(
                    ArgumentError(
                        "$caller: Curve[$et] fillet tangent winds " *
                        "non-monotonically on the periodic face"))
            end
            close=vals[1]+sgn*info.period-vals[end]
            close*sgn>-1e-9*info.period || throw(ArgumentError(
                "$caller: Curve[$et] fillet tangent overshoots one " *
                "period — unsupported"))
            abs(close)<0.75*info.period || throw(ArgumentError(
                "$caller: Curve[$et] fillet tangent does not wind " *
                "once around the periodic face — unsupported"))
        end
    end
    # Seam crossing azimuth φ_s: any face with a parametric seam fixes
    # the edge's start parameter; two seams must agree.
    φs=0.0; have_s=false
    for info in (infoA,infoB)
        (info.perdim==0 || isnan(info.seam_uv)) && continue
        d=info.perdim
        ser=tan_uv[info.fid]
        a=[q[d] for q in ser]
        # target ≡ seam_uv (mod period) interior to [a[1], a[1]+span]
        span=a[end]-a[1]
        sgn=sign(span)
        k=round((a[1]+0.5*sgn*info.period-info.seam_uv)/info.period)
        target=info.seam_uv+k*info.period
        sgn*(target-a[1])>0 && sgn*(a[end]-target)>0 ||
            (target-=sgn*info.period)
        sgn*(target-a[1])>=0 && sgn*(a[end]-target)>=0 ||
            throw(ArgumentError(
                "$caller: Curve[$et] fillet tangent misses the face " *
                "seam"))
        # locate the bracketing interval
        ib=findfirst(i->(a[i]-target)*sgn<=0 &&
                        (a[i+1]-target)*sgn>=0,1:length(a)-1)
        ib===nothing && throw(ArgumentError(
            "$caller: Curve[$et] fillet tangent misses the face " *
            "seam"))
        φlo=sts[ib].θ; φhi=sts[ib+1].θ
        # tighten the crossing on the raw (wrapped) coordinate by
        # bisection against the mod-period value
        info_r=info.fid==fA ? infoB : infoA
        function gapf(θ)
            st2=_occ_fillet_sweep_stations(kernel,o,axis,X,Y,eg,R,r,
                                           info,info_r,(θ,),et,caller)
            # `info` is the first argument, so uvA is its tangent uv.
            uv=st2[1].uvA
            raw=info.perdim==1 ? uv[1] : uv[2]
            uvw=raw+info.period*round((target-raw)/info.period)
            return (uvw-target)*sgn
        end
        for _ in 1:60
            φm=(φlo+φhi)/2; gm=gapf(φm)
            abs(gm)<1e-11 && break
            gm<0 ? (φlo=φm) : (φhi=φm)
        end
        φloc=(φlo+φhi)/2
        if !have_s
            φs=φloc; have_s=true
        else
            δ=mod(φloc-φs+π,_OCC_TWO_PI)-π
            abs(δ)<1e-6 || throw(ArgumentError(
                "$caller: Curve[$et] fillet tangents cross " *
                "misaligned face seams — unsupported"))
        end
    end
    # Build pass: stations on the seam-aligned grid, shared by the edge
    # interpolation and the fillet net.
    phis2=ntuple(i->φs+(i-1)*_OCC_TWO_PI/N,N)
    sts2=_occ_fillet_sweep_stations(kernel,o,axis,X,Y,eg,R,r,
                                    infoA,infoB,phis2,et,caller)
    rimvid=kernel.shapes[eid].children[1]
    return (axis=axis,o=o,R=R,rho=r,rim=(R,0.0),X=X,Y=Y,
            rimvid=rimvid,rim_vpos=vert_pos[rimvid],
            fA=fA,fB=fB,infoA=infoA,infoB=infoB,phs=φs,
            sts=sts2,tan_uv=tan_uv,rim_uv=rim_uv)
end

# Publish the swept fillets: tangent NURBS edges replace the rims on
# the kept faces (wires/pcurves/vertex_uv re-bound on fresh shapes —
# plane records demote to their trimmed-tensor form so the curved
# boundary keeps a surface parametrization), meridian seam arcs and
# tensor pipe faces bound like the torus path.
function _occ_fillet_sweep_rebuild!(m,kernel,specs,wires_of,frev_of,
                                   new_face2,shell_uses,caller)
    isempty(specs) && return nothing
    rim_repl=Dict{Tuple{UInt64,UInt64},NamedTuple}()
    vmap=Dict{Tuple{UInt64,UInt64},UInt64}()
    sweep_faces=Dict{Tuple{UInt64,UInt64},NamedTuple}()
    for (eid,s) in specs
        # Tangent edge records, one per incident face; both are
        # parameterized t ∈ [0,2π] starting at the seam azimuth φ_s.
        st0=s.sts[1]
        for (fid,key) in ((s.fA,:TA),(s.fB,:TB))
            pts=NTuple{3,Float64}[getproperty(st,key) for st in s.sts]
            rec=_occ_fillet_sweep_edge(pts)
            # Mint the vertex first — add_point internally reserves its
            # shape id against kernel.next_shape, which must already be
            # past every edge bound so far.
            vtx=getproperty(st0,key)
            vtag=_model_occ_add_point!(m,vtx[1],vtx[2],vtx[3])
            vnid=kernel.bindings[(0,vtag)]
            eidn=only(_model_occ_reserve_shapes(kernel,1,caller))
            etag=_model_occ_next_tag(kernel,1,-1,caller)
            _model_occ_bind_shape!(kernel,1,etag,eidn,
                _ModelOCCShape(1,(kind=:tensor_bspline,geometry=rec),
                               UInt64[vnid,vnid]))
            kernel.next_shape=eidn+UInt64(1)
            info=fid==s.fA ? s.infoA : s.infoB
            uv0=getproperty(st0,fid==s.fA ? :uvA : :uvB)
            # Vertex uv on this face: the seam vertex sits at both ends
            # of the closed tangent trace — store the trailing domain
            # bound like the loft's u=hi convention, or the raw uv when
            # the face is not periodic.
            vuv=uv0
            if info.perdim!=0
                d=info.perdim
                ser=s.tan_uv[fid]
                sgn=sign(ser[end][d]-ser[1][d])
                uend=sgn>=0 ? info.hi[d] : info.lo[d]
                vuv=d==1 ? (uend,uv0[2]) : (uv0[1],uend)
            end
            sweep_faces[(fid,eid)]=(eid_new=eidn,etag=etag,vid=vnid,
                                    vtag=vtag,vuv=vuv)
            vmap[(fid,s.rimvid)]=vnid
        end
    end
    # Kept faces: swap rim uses for the tangent edges, truncate seam
    # edges at the tangent vertices, rekey pcurves/vertex_uv, rebind.
    dirty=Set{UInt64}(fid for (fid,_) in keys(sweep_faces))
    for fid in dirty
        fgsh=kernel.shapes[fid]; fg=fgsh.geometry.geometry
        ftag=_model_occ_shape_tag(kernel,fid,caller)
        ispl=fg.occ===:plane
        pco=hasproperty(fg,:pcurves) ? copy(fg.pcurves) :
            Dict{Int,NamedTuple}()
        vuvo=hasproperty(fg,:vertex_uv) ? copy(fg.vertex_uv) :
            Dict{Int,NTuple{2,Float64}}()
        new_wids=UInt64[]
        # A parametric seam appears twice in the same wire — both uses
        # must reference the ONE truncated edge entity.
        truncmap=Dict{UInt64,UInt64}()
        for wid in wires_of[fid]
            uses=kernel.shapes[wid].geometry.uses
            changed=false; nuses=Tuple{UInt64,Bool}[]
            for (u,rev) in uses
                repl=get(sweep_faces,(fid,u),nothing)
                if repl!==nothing
                    push!(nuses,(repl.eid_new,rev)); changed=true
                    continue
                end
                ch=kernel.shapes[u].children
                c1=get(vmap,(fid,ch[1]),ch[1])
                c2=get(vmap,(fid,ch[2]),ch[2])
                if c1==ch[1] && c2==ch[2]
                    push!(nuses,(u,rev)); continue
                end
                changed=true
                neid=get(truncmap,u,nothing)
                if neid===nothing
                    g2=kernel.shapes[u].geometry.geometry
                    tA=_occ_2df_param(g2,
                        kernel.shapes[c1].geometry.point)
                    tB=_occ_2df_param(g2,
                        kernel.shapes[c2].geometry.point)
                    tA<tB-1e-12 || throw(ArgumentError(
                        "$caller: fillet truncation inverts an edge " *
                        "of Surface[$ftag]"))
                    rec=merge(g2,(t0=tA,t1=tB))
                    hasproperty(rec,:native_geometry) && (rec=merge(rec,
                        (native_geometry=merge(rec.native_geometry,
                            (t0=tA,t1=tB)),)))
                    neid=only(_model_occ_reserve_shapes(kernel,1,caller))
                    netag=_model_occ_next_tag(kernel,1,-1,caller)
                    _model_occ_bind_shape!(kernel,1,netag,neid,
                        _ModelOCCShape(1,
                            (kind=kernel.shapes[u].geometry.kind,
                             geometry=rec),UInt64[c1,c2]))
                    kernel.next_shape=neid+UInt64(1)
                    truncmap[u]=neid
                    oldt=_model_occ_shape_tag(kernel,u,caller)
                    haskey(pco,oldt) && (pco[netag]=pop!(pco,oldt))
                end
                push!(nuses,(neid,rev))
            end
            if changed
                wt=_occ_wire_from_uses!(m,nuses,caller)
                push!(new_wids,kernel.bindings[(-1,wt)])
            else
                push!(new_wids,wid)
            end
        end
        # Tangent-edge pcurves: spl2d through the station uv trace —
        # pchip keeps the trace inside the sampled window.
        for ((ffid,ueid),repl) in sweep_faces
            ffid===fid || continue
            s=specs[ueid]
            info=fid==s.fA ? s.infoA : s.infoB
            ser=NTuple{2,Float64}[]
            n=length(s.sts)
            for (i,st) in enumerate(s.sts)
                uv=fid==s.fA ? st.uvA : st.uvB
                push!(ser,uv)
            end
            # Periodic faces: remap the unwrapped trace so the seam
            # vertex starts on the domain bound matching the winding
            # direction, then runs one full period — the same layout as
            # the rim's own tan2d/lin2d pcurve on the face.  The trace
            # itself carries the winding sign, so the map is a pure
            # translation `ud = u_start + (uw − uw_0)`.
            if info.perdim!=0
                d=info.perdim
                sgn=sign(ser[end][d]-ser[1][d])
                u0=sgn>=0 ? info.lo[d] : info.hi[d]
                q0=ser[1][d]
                ser=[d==1 ? (u0+q[1]-q0,q[2]) : (q[1],u0+q[2]-q0)
                     for q in ser]
            end
            ts=[s.sts[i].θ-s.phs for i in 1:n]
            # Closing sample at t=2π: the same vertex as t=0 but on the
            # far side of the periodic trace — continue the wind rather
            # than snapping back to the first sample.
            if info.perdim!=0
                d=info.perdim
                sgn=sign(ser[end][d]-ser[1][d])
                qe=ser[1]
                qe=d==1 ? (qe[1]+sgn*info.period,qe[2]) :
                          (qe[1],qe[2]+sgn*info.period)
                push!(ser,qe)
            else
                push!(ser,ser[1])
            end
            push!(ts,2π)
            pc=(spl2d=_pchip_bspline_2d(ts,ser),)
            delete!(pco,_model_occ_shape_tag(kernel,ueid,caller))
            pco[repl.etag]=(fwd=pc,rev=nothing)
            ovt=_model_occ_shape_tag(kernel,
                kernel.shapes[ueid].children[1],caller)
            haskey(vuvo,ovt) && delete!(vuvo,ovt)
            vuvo[repl.vtag]=repl.vuv
        end
        # Rebind: planes demote to the trimmed degree-1 tensor patch
        # over the stored frame (the loft cap fallback), other records
        # keep their geometry with the new pcurve/vertex_uv tables.
        if ispl
            ulo,uhi=fg.uvb[1][1],fg.uvb[2][1]
            vlo,vhi=fg.uvb[1][2],fg.uvb[2][2]
            ctrl=Matrix{NTuple{3,Float64}}(undef,2,2)
            wts=ones(Float64,2,2)
            for (i,uu) in enumerate((ulo,uhi)),
                    (j,vv) in enumerate((vlo,vhi))
                p=fg.center .+ uu.*fg.X .+ vv.*fg.Y
                ctrl[i,j]=(p[1],p[2],p[3])
            end
            patch=NURBS.NURBSSurface(1,1,Float64[ulo,ulo,uhi,uhi],
                Float64[vlo,vlo,vhi,vhi],ctrl,wts)
            # plane_rec keeps the analytic frame for the measure fast
            # path, but its pcurve table must be the rebuilt one — the
            # wire's new tangent edges carry fresh tags.
            newgeo=(occ=:tensor_bspline,
                patch=_ModelNURBSSurface(:bspline,patch),
                pcurves=pco,vertex_uv=vuvo,
                reversed=get(fg,:reversed,false),
                trimmed=true,
                plane_rec=merge(fg,(pcurves=pco,)))
            nft=_model_occ_next_tag(kernel,2,-1,caller)
            nfid=only(_model_occ_reserve_shapes(kernel,1,caller))
            _model_occ_bind_shape!(kernel,2,nft,nfid,
                _ModelOCCShape(2,(kind=:tensor_bspline,geometry=newgeo),
                               new_wids))
            kernel.next_shape=nfid+UInt64(1)
            new_face2[fid]=(nfid,frev_of[fid])
        else
            # The rebuilt wire bounds a sub-region of the patch's knot
            # domain (the strip between the two tangent traces, split by
            # the truncated seam) — the measure and uv-mesher paths must
            # take the trimmed-CDT route, not the full-domain integral.
            newgeo=merge(fg,(pcurves=pco,trimmed=true))
            newgeo=merge(newgeo,(vertex_uv=vuvo,))
            # An already-demoted cap carries `plane_rec` for the measure
            # fast path; its nested pcurve table is keyed by the same
            # edge tags and must be rebuilt alongside the face's own.
            if hasproperty(newgeo,:plane_rec) &&
                    hasproperty(newgeo.plane_rec,:pcurves)
                newgeo=merge(newgeo,(plane_rec=
                    merge(newgeo.plane_rec,(pcurves=pco,)),))
            end
            nft=_model_occ_next_tag(kernel,2,-1,caller)
            nfid=only(_model_occ_reserve_shapes(kernel,1,caller))
            _model_occ_bind_shape!(kernel,2,nft,nfid,
                _ModelOCCShape(2,(kind=fgsh.geometry.kind,
                                  geometry=newgeo),new_wids))
            kernel.next_shape=nfid+UInt64(1)
            new_face2[fid]=(nfid,frev_of[fid])
        end
    end
    # Pipe fillet faces — one per swept rim, 4-use periodic wire like
    # the torus fillet face.
    for (eid,s) in specs
        N=length(s.sts)
        h=_OCC_TWO_PI/N
        # Meridian arc control rows in homogeneous space.
        rows=[NTuple{4,Float64}[],NTuple{4,Float64}[],
              NTuple{4,Float64}[]]
        dphi0=s.sts[1].dφ
        for st in s.sts
            w=cos(st.dφ/2)
            e1=_occ_dir(_arc_sub(st.TA,st.C))
            e2=_occ_dir(_arc_sub(st.TB,st.C))
            mid=_occ_dir(_occ3_add(e1,e2))
            cp=_occ3_add(st.C,_occ3_scale(mid,s.rho/w))
            push!(rows[1],(st.TA[1],st.TA[2],st.TA[3],1.0))
            push!(rows[2],(cp[1]*w,cp[2]*w,cp[3]*w,w))
            push!(rows[3],(st.TB[1],st.TB[2],st.TB[3],1.0))
        end
        nets=[_fillet_periodic_poles(rows[j],caller) for j in 1:3]
        # Abscissa convention: S(u) interpolates controls[k] at
        # u=(k−2)·h, so the net is rotated one pole to put station 1 —
        # the seam azimuth φ_s — at u=0, where the seam arc's pcurve
        # claims it.
        controls=Matrix{NTuple{3,Float64}}(undef,N+3,3)
        wmat=Matrix{Float64}(undef,N+3,3)
        for i in 0:N+2, j in 1:3
            q=nets[j][mod(i-1,N)+1]
            controls[i+1,j]=(q[1]/q[4],q[2]/q[4],q[3]/q[4])
            wmat[i+1,j]=q[4]
        end
        surf=NURBS.NURBSSurface(3,2,_fillet_periodic_knots(N),
            Float64[0,0,0,1,1,1],controls,wmat)
        rA=sweep_faces[(s.fA,eid)]; rB=sweep_faces[(s.fB,eid)]
        # Meridian seam arc at φ_s: circle centred on the ball with the
        # tangent direction frame.
        B0=s.sts[1].C
        xa=_occ_dir(_arc_sub(s.sts[1].TA,B0))
        ya_raw=_arc_sub(s.sts[1].TB,B0)
        ya=_occ_dir(_arc_sub(ya_raw,
            _occ3_scale(xa,_dot3(ya_raw,xa))))
        seam_rec=(occ=:circle,center=B0,n=_occ_cross(xa,ya),
                  X=xa,Y=ya,r=s.rho,t0=0.0,t1=dphi0)
        seid=only(_model_occ_reserve_shapes(kernel,1,caller))
        setag=_model_occ_next_tag(kernel,1,-1,caller)
        _model_occ_bind_shape!(kernel,1,setag,seid,
            _ModelOCCShape(1,(kind=:circle,geometry=seam_rec),
                           UInt64[rA.vid,rB.vid]))
        kernel.next_shape=seid+UInt64(1)
        etagA=rA.etag; etagB=rB.etag
        pc=Dict{Int,NamedTuple}(
            etagA=>(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing),
            etagB=>(fwd=(lin2d=((0.0,1.0),(1.0,0.0)),),rev=nothing),
            setag=>(fwd=(lin2d=((_OCC_TWO_PI,0.0),(0.0,1.0/dphi0)),),
                    rev=(lin2d=((0.0,0.0),(0.0,1.0/dphi0)),)))
        wuses=Tuple{UInt64,Bool}[(rB.eid_new,true),(seid,true),
                                 (rA.eid_new,false),(seid,false)]
        wt=_occ_wire_from_uses!(m,wuses,caller)
        wid2=kernel.bindings[(-1,wt)]
        # Outward sign: the raw Su×Sv should point away from the ball
        # centre — the pipe surface's material side faces the tube.
        bare=(patch=_ModelNURBSSurface(:bspline,surf),)
        pmid=_occ_surface_point(bare,Float64(π),0.5)
        _,fdu,fdv=_occ_surface_d1(bare,Float64(π),0.5)
        fn=_occ_cross(fdu,fdv)
        Bout=s.sts[div(N,2)+1].C
        frev=_dot3(fn,_arc_sub(pmid,Bout))<0
        vuv=Dict{Int,NTuple{2,Float64}}(
            rA.vtag=>(_OCC_TWO_PI,0.0),rB.vtag=>(_OCC_TWO_PI,1.0))
        nft=_model_occ_next_tag(kernel,2,-1,caller)
        nfid=only(_model_occ_reserve_shapes(kernel,1,caller))
        _model_occ_bind_shape!(kernel,2,nft,nfid,
            _ModelOCCShape(2,(kind=:tensor_bspline,
                geometry=(occ=:tensor_bspline,
                    patch=_ModelNURBSSurface(:bspline,surf),
                    pcurves=pc,vertex_uv=vuv,reversed=frev)),
                UInt64[wid2]))
        kernel.next_shape=nfid+UInt64(1)
        push!(shell_uses,(nfid,frev))
    end
    return nothing
end

function _model_occ_fillet(m::GeoModel,volume_tags,curve_tags,radii;
                           caller="OCC.fillet")
    kernel=m.occ_kernel
    rs=radii isa AbstractVector ? Float64[Float64(r) for r in radii] :
        fill(Float64(radii),length(curve_tags))
    length(rs)==length(curve_tags) || throw(ArgumentError(
        "$caller: radii and curveTags must have the same length"))
    all(>(0.0),rs) || throw(ArgumentError(
        "$caller: fillet radii must be positive"))
    out=Tuple{Int,Int}[]
    for volume_tag in volume_tags
        vt=abs(_model_occ_requested_tag(volume_tag,caller))
        vid=get(kernel.bindings,(3,vt),nothing)
        vid===nothing && throw(ArgumentError(
            "$caller: unknown OCC Volume[$vt]"))
        fnormals=Dict{UInt64,NTuple{3,Float64}}()
        fcenters=Dict{UInt64,NTuple{3,Float64}}()
        edge_faces=Dict{UInt64,Vector{UInt64}}()
        vert_edges=Dict{UInt64,Vector{UInt64}}()
        vert_pos=Dict{UInt64,NTuple{3,Float64}}()
        wire_of=Dict{UInt64,UInt64}()
        wires_of=Dict{UInt64,Vector{UInt64}}()
        frev_of=Dict{UInt64,Bool}()
        planes=NTuple{3,Float64}[]; bs=Float64[]
        shell_uses=copy(kernel.shapes[
            kernel.shapes[vid].children[1]].geometry.uses)
        isplane=Dict{UInt64,Bool}()
        for (fid,rev) in shell_uses
            frev_of[fid]=rev
            fg=kernel.shapes[fid].geometry
            isplane[fid]=get(fg.geometry,:occ,nothing)===:plane
            if isplane[fid]
                σ=rev ? -1.0 : 1.0
                n=_occ3_scale(fg.geometry.axis,σ)
                get(fg.geometry,:reversed,false) && (n=(.-n))
                n=_occ3_scale(n,1.0/_occ_modulus(n))
                fnormals[fid]=n; fcenters[fid]=fg.geometry.center
                push!(planes,n); push!(bs,_dot3(n,fg.geometry.center))
            end
            wlist=kernel.shapes[fid].children
            wires_of[fid]=wlist
            wire_of[fid]=wlist[1]
            for wid in wlist, (u,_) in kernel.shapes[wid].geometry.uses
                push!(get!(edge_faces,u,UInt64[]),fid)
                for c in kernel.shapes[u].children
                    push!(get!(vert_edges,c,UInt64[]),u)
                    vert_pos[c]=kernel.shapes[c].geometry.point
                end
            end
        end
        # The half-space boundedness/orientation certificate needs every
        # face planar — a curved shell's planar subset needn't bound a
        # region, so the check (and the inward-flip retry) stays limited
        # to the all-planar case it was designed for.
        if all(values(isplane))
            overs=_occ_halfspace_verts(planes,bs)
            if isempty(overs)
                planes=[(.-n) for n in planes]; bs=(.-bs)
                fnormals=Dict(k=>(.-v) for (k,v) in fnormals)
                overs=_occ_halfspace_verts(planes,bs)
                isempty(overs) && throw(ArgumentError(
                    "$caller: Volume[$vt] faces bound an empty body"))
            end
        end
        fil_verts=Set{UInt64}()
        trunc=Dict{Tuple{UInt64,UInt64},NTuple{3,Float64}}()
        tangent=Dict{Tuple{UInt64,UInt64},Dict{UInt64,NTuple{3,Float64}}}()
        arc_at=Dict{Tuple{UInt64,UInt64},NTuple{3,Float64}}()
        fillet_of=Dict{UInt64,NamedTuple}()
        circ_specs=Dict{UInt64,NamedTuple}()
        sweep_specs=Dict{UInt64,NamedTuple}()
        for (ct,r) in zip(curve_tags,rs)
            et=abs(_model_occ_requested_tag(ct,caller))
            eid=get(kernel.bindings,(1,et),nothing)
            eid===nothing && throw(ArgumentError(
                "$caller: unknown OCC Curve[$et]"))
            fids=get(edge_faces,eid,nothing)
            (fids===nothing || length(fids)!=2) && throw(ArgumentError(
                "$caller: Curve[$et] is not an edge shared by two " *
                "faces of Volume[$vt]"))
            eg=kernel.shapes[eid].geometry.geometry
            eocc=get(eg,:occ,nothing)
            # Parametric seams bound the same face twice — a periodic
            # parametrization line or a degenerate meridian arc, not a
            # material edge. OCC's fillet builder skips them silently
            # (a smooth G1 edge carries no blend); mirror that for any
            # curve kind.
            if fids[1]==fids[2]
                continue
            end
            if eocc===:circle
                for vv in unique(kernel.shapes[eid].children)
                    vv in fil_verts && throw(ArgumentError(
                        "$caller: filleted edges share a vertex — " *
                        "corner blends are unsupported"))
                    push!(fil_verts,vv)
                end
                # Coaxial revolution faces get the exact torus path;
                # anything else with evaluators/pcurves (tensor
                # lofts, trimmed caps, off-axis analytics) takes the
                # swept rolling-ball path.
                anal(k)=k in (:plane,:cylinder,:cone,:sphere,:torus)
                kA=get(kernel.shapes[fids[1]].geometry.geometry,:occ,
                       nothing)
                kB=get(kernel.shapes[fids[2]].geometry.geometry,:occ,
                       nothing)
                if anal(kA) && anal(kB)
                    circ_specs[eid]=_occ_fillet_circle_plan(
                        kernel,eg,eid,fids,frev_of,vert_pos,et,r,
                        caller)
                else
                    sweep_specs[eid]=_occ_fillet_sweep_plan(
                        kernel,eg,eid,fids,frev_of,vert_pos,
                        wires_of,et,r,caller)
                end
                continue
            end
            eocc===:line || throw(ArgumentError(
                "$caller: fillet edge must be a line or a full circle " *
                "on coaxial revolution faces (Curve[$et] is occ=$eocc)"))
            t=_occ3_scale(eg.dir,1.0/_occ_modulus(eg.dir))
            va,vb=kernel.shapes[eid].children
            pa,pb=vert_pos[va],vert_pos[vb]
            em=_occ3_scale(_occ3_add(pa,pb),0.5)
            L=_occ_modulus(_occ3_sub(pb,pa))
            fA,fB=fids
            for f in (fA,fB)
                focc=get(kernel.shapes[f].geometry.geometry,:occ,nothing)
                isplane[f] || throw(ArgumentError(
                    "$caller: Curve[$et] touches a non-planar face " *
                    "(occ=$focc) — surface blends (ChFi3d) are unsupported"))
            end
            nA,nB=fnormals[fA],fnormals[fB]
            φ=acos(clamp(_dot3(nA,nB),-1.0,1.0))
            δ=π-φ
            (δ>0.05 && δ<π-0.05) || throw(ArgumentError(
                "$caller: degenerate dihedral on Curve[$et]"))
            leg=r/tan(δ/2)
            h=r/sin(δ/2)
            nb=_occ3_add(nA,nB)
            bdir=_occ3_scale(nb,-1.0/_occ_modulus(nb))
            pin=_occ3_add(em,_occ3_scale(bdir,0.01))
            (_dot3(nA,pin)<=_dot3(nA,fcenters[fA])+1e-9 &&
             _dot3(nB,pin)<=_dot3(nB,fcenters[fB])+1e-9) ||
                throw(ArgumentError(
                    "$caller: Curve[$et] is a concave edge — " *
                    "concave fillets are unsupported"))
            uA=_occ_cross(t,nA); uB=_occ_cross(t,nB)
            uA=_occ3_scale(uA,1.0/_occ_modulus(uA))
            uB=_occ3_scale(uB,1.0/_occ_modulus(uB))
            _dot3(uA,_occ3_sub(fcenters[fA],em))<0 && (uA=(.-uA))
            _dot3(uB,_occ3_sub(fcenters[fB],em))<0 && (uB=(.-uB))
            for vv in (va,vb)
                vv in fil_verts && throw(ArgumentError(
                    "$caller: filleted edges share a vertex — " *
                    "corner blends are unsupported"))
            end
            push!(fil_verts,va); push!(fil_verts,vb)
            ends=NamedTuple[]
            for vv in (va,vb)
                pv=vert_pos[vv]
                sibs=unique(s for s in vert_edges[vv] if s!=eid)
                length(sibs)==2 || throw(ArgumentError(
                    "$caller: Curve[$et] endpoint is not a simple " *
                    "corner ($(length(sibs)) sibling edges)"))
                iA=findfirst(s->fA in edge_faces[s],sibs)
                iA===nothing && throw(ArgumentError(
                    "$caller: broken corner topology on Curve[$et]"))
                sA=sibs[iA]; sB=sibs[3-iA]
                fC=only(f for f in edge_faces[sA] if f!=fA)
                isplane[fC] || throw(ArgumentError(
                    "$caller: Curve[$et] corner touches a non-planar " *
                    "face — surface blends (ChFi3d) are unsupported"))
                nC=fnormals[fC]
                abs(_dot3(t,nC))>1-1e-6 || throw(ArgumentError(
                    "$caller: fillet end face must be orthogonal " *
                    "to the edge"))
                oA=only(c for c in kernel.shapes[sA].children
                        if c!=vv)
                oB=only(c for c in kernel.shapes[sB].children
                        if c!=vv)
                dA=_occ3_sub(vert_pos[oA],pv)
                dB=_occ3_sub(vert_pos[oB],pv)
                lA=_occ_modulus(dA); lB=_occ_modulus(dB)
                dA=_occ3_scale(dA,1.0/lA); dB=_occ3_scale(dB,1.0/lB)
                (leg<lA-1e-9 && leg<lB-1e-9) || throw(ArgumentError(
                    "$caller: fillet on Curve[$et] overruns an " *
                    "adjacent edge"))
                tA=_occ3_add(pv,_occ3_scale(dA,leg))
                tB=_occ3_add(pv,_occ3_scale(dB,leg))
                s=_occ3_add(pv,_occ3_scale(bdir,h))
                trunc[(sA,vv)]=tA; trunc[(sB,vv)]=tB
                arc_at[(fC,vv)]=s
                push!(ends,(v=vv,tA=tA,tB=tB,s=s,fC=fC,
                            sA=sA,sB=sB))
            end
            tangent[(fA,eid)]=Dict(e.v=>e.tA for e in ends)
            tangent[(fB,eid)]=Dict(e.v=>e.tB for e in ends)
            fillet_of[eid]=(t=t,pa=pa,pb=pb,va=va,vb=vb,r=r,
                            bdir=bdir,fA=fA,fB=fB,ends=ends)
        end
        # A same-face seam edge is skipped per-edge (OCC `Add` treats
        # it as unsuitable); when nothing requested was suitable, OCC
        # reports "no suitable edges for chamfer or fillet" — the
        # guard mirrors that, covering both all-seam lists and the
        # empty input.
        isempty(fillet_of) && isempty(circ_specs) &&
            isempty(sweep_specs) && throw(ArgumentError(
            "$caller: no valid edges to fillet"))
        # Strip-overrun certificates for the circle fillets, now that
        # every filleted rim on each face is known.
        for (eid,s) in circ_specs
            et=_model_occ_shape_tag(kernel,eid,caller)
            for (fid,p,other,T,leg) in ((s.fA,s.pA,s.pB,s.TA,s.legA),
                                        (s.fB,s.pB,s.pA,s.TB,s.legB))
                extra=NTuple{2,Float64}[]
                for (eid2,s2) in circ_specs
                    eid2===eid && continue
                    s2.fA===fid && push!(extra,s2.TA)
                    s2.fB===fid && push!(extra,s2.TB)
                end
                _occ_fillet_overrun(kernel,fid,eid,p,s.rim,other.nin,leg,
                    extra,s.axis,s.o,wires_of,et,caller)
            end
        end
        # uv-strip overruns for the swept fillets: every other boundary
        # of each kept face and every other tangent trace on it must
        # stay clear of the strip the fillet removes.
        for (eid,s) in sweep_specs
            et=_model_occ_shape_tag(kernel,eid,caller)
            for fid in (s.fA,s.fB)
                _occ_fillet_sweep_overrun(kernel,fid,eid,s,sweep_specs,
                    wires_of,et,caller)
            end
        end
        # ── rebuild geometry ─────────────────────────────────────────
        ptag=Dict{NTuple{3,Float64},Int}()
        function _pt(p)
            get!(ptag,p) do
                _model_occ_add_point!(m,p[1],p[2],p[3])
            end
        end
        line_id=Dict{NTuple{2,UInt64},UInt64}()
        function _line(p0,p1)
            i0=kernel.bindings[(0,_pt(p0))]
            i1=kernel.bindings[(0,_pt(p1))]
            key=i0<i1 ? (i0,i1) : (i1,i0)
            id=get!(line_id,key) do
                lt=_model_occ_add_line!(m,_pt(p0),_pt(p1))
                kernel.bindings[(1,lt)]
            end
            # rev: stored children order vs requested traversal
            id, kernel.shapes[id].children[1]!=i0
        end
        arc_id=Dict{Tuple{UInt64,UInt64,UInt64},UInt64}()
        function _arc(p0,p1,cen)
            i0=kernel.bindings[(0,_pt(p0))]
            i1=kernel.bindings[(0,_pt(p1))]
            ic=kernel.bindings[(0,_pt(cen))]
            key=i0<i1 ? (i0,i1,ic) : (i1,i0,ic)
            id=get!(arc_id,key) do
                at=_model_occ_add_circle_arc!(m,_pt(p0),_pt(cen),
                                              _pt(p1))
                kernel.bindings[(1,at)]
            end
            id, kernel.shapes[id].children[1]!=i0
        end
        new_face=Dict{UInt64,Int}()
        for (fid,_) in shell_uses
            wid=wire_of[fid]
            uses=kernel.shapes[wid].geometry.uses
            dirty=false
            for (u,_) in uses
                haskey(tangent,(fid,u)) && (dirty=true; break)
                for c in kernel.shapes[u].children
                    (haskey(trunc,(u,c)) || haskey(arc_at,(fid,c))) &&
                        (dirty=true; break)
                end
                dirty && break
            end
            dirty || continue
            focc=get(kernel.shapes[fid].geometry.geometry,:occ,nothing)
            isplane[fid] || throw(ArgumentError(
                "$caller: fillet truncates an edge of a non-planar face " *
                "(occ=$focc) — surface blends are unsupported"))
            segs=Any[]
            prev_end=nothing
            for (u,rev) in uses
                ch=kernel.shapes[u].children
                oa,ob=rev ? (ch[2],ch[1]) : (ch[1],ch[2])
                if haskey(tangent,(fid,u))
                    p0=tangent[(fid,u)][oa]
                    p1=tangent[(fid,u)][ob]
                else
                    p0=get(trunc,(u,oa),vert_pos[oa])
                    p1=get(trunc,(u,ob),vert_pos[ob])
                end
                if prev_end!==nothing &&
                        _occ_modulus(_occ3_sub(p0,prev_end))>1e-9
                    s=arc_at[(fid,oa)]
                    push!(segs,(:arc,prev_end,p0,s))
                end
                push!(segs,(:line,p0,p1))
                prev_end=p1
            end
            # Wrap-around gap: last segment end → first segment start.
            # The junction vertex is the endpoint shared by the first
            # and last wire edges.
            if length(segs)>1 && segs[1][1]==:line
                p0=segs[1][2]
                if _occ_modulus(_occ3_sub(p0,prev_end))>1e-9
                    c1=kernel.shapes[uses[1][1]].children
                    cN=kernel.shapes[uses[end][1]].children
                    vj=only(v for v in c1 if v in cN)
                    s=arc_at[(fid,vj)]
                    pushfirst!(segs,(:arc,prev_end,p0,s))
                end
            end
            wuses=Tuple{UInt64,Bool}[]
            for seg in segs
                if seg[1]==:line
                    lid,lr=_line(seg[2],seg[3])
                    push!(wuses,(lid,lr))
                else
                    aid,ar=_arc(seg[2],seg[3],seg[4])
                    push!(wuses,(aid,ar))
                end
            end
            wt=_occ_wire_from_uses!(m,wuses,caller)
            new_face[fid]=_model_occ_add_plane_surface!(m,[wt])
        end
        # Circle-edge fillets: conflict with the planar line-edge rebuild
        # is a shared dirty face — reject loudly rather than tangle.
        for (eid,s) in circ_specs, fid in (s.fA,s.fB)
            haskey(new_face,fid) && throw(ArgumentError(
                "$caller: fillet mixes line and circle edges on one " *
                "face — unsupported"))
        end
        for (eid,s) in sweep_specs, fid in (s.fA,s.fB)
            haskey(new_face,fid) && throw(ArgumentError(
                "$caller: fillet mixes line and circle edges on one " *
                "face — unsupported"))
        end
        new_face2=Dict{UInt64,Tuple{UInt64,Bool}}()
        _occ_fillet_circle_rebuild!(m,kernel,circ_specs,wires_of,
            frev_of,fnormals,new_face2,shell_uses,caller)
        _occ_fillet_sweep_rebuild!(m,kernel,sweep_specs,wires_of,
            frev_of,new_face2,shell_uses,caller)
        for (eid,spec) in fillet_of
            e0,e1=spec.ends
            lidA,dirA=_line(e0.tA,e1.tA)
            lidB,dirB=_line(e0.tB,e1.tB)
            a0,dir0=_arc(e0.tA,e0.tB,e0.s)
            a1,dir1=_arc(e1.tA,e1.tB,e1.s)
            s0=e0.s
            X=_occ3_scale(_occ3_sub(e0.tA,s0),1.0/spec.r)
            # The record frame must stay right-handed (X×Y=axis): when
            # t_B lies on the -Y side, flip the axis rather than Y so
            # the sweep reads positive.
            axisc=spec.t
            Y=_occ_crossed(axisc,X)
            db=_occ3_sub(e0.tB,s0)
            θB=atan(_dot3(db,Y),_dot3(db,X))
            if θB<0
                axisc=(.-axisc); Y=(.-Y)
                θB=-θB
            end
            L=_occ_modulus(_occ3_sub(spec.pb,spec.pa))
            pc=Dict{Int,NamedTuple}()
            tagA=_model_occ_shape_tag(kernel,lidA,caller)
            tagB=_model_occ_shape_tag(kernel,lidB,caller)
            # v is the axial coordinate along `axisc`; the stored line
            # direction may be either way (shared cached lines), so the
            # pcurve v(t)=dot(origin−s0,axisc)+dot(dir,axisc)·t.
            for (lid,u_c) in ((lidA,0.0),(lidB,θB))
                lg=kernel.shapes[lid].geometry.geometry
                v0=_dot3(_occ3_sub(lg.origin,s0),axisc)
                dv=_dot3(lg.dir,axisc)
                ltag=_model_occ_shape_tag(kernel,lid,caller)
                pc[ltag]=(fwd=(lin2d=((u_c,v0),(0.0,dv)),),rev=nothing)
            end
            # Arc edge parameters are global circle angles — affine-remap
            # u so the t_A endpoint sits at u=0 and t_B at u=θB.
            for (aid,tA,tB) in ((a0,e0.tA,e0.tB),(a1,e1.tA,e1.tB))
                ag=kernel.shapes[aid].geometry.geometry
                p0c=_occ_circle_point(ag,ag.t0)
                a_off=_occ_modulus(_occ3_sub(p0c,tA))<
                      _occ_modulus(_occ3_sub(p0c,tB)) ? ag.t0 : ag.t1
                a_op=ag.t0+ag.t1-a_off
                du=a_op>a_off ? 1.0 : -1.0
                vv=_dot3(_occ3_sub(tA,s0),axisc)
                atag=_model_occ_shape_tag(kernel,aid,caller)
                pc[atag]=(fwd=(lin2d=((-du*a_off,vv),(du,0.0)),),
                          rev=nothing)
            end
            # uv loop (0,0)→(θB,0)→(θB,L)→(0,L): T_A up, arc1 across
            # the top, T_B back down, arc0 back across the bottom.
            wuses=Tuple{UInt64,Bool}[(lidA,dirA),(a1,dir1),
                                     (lidB,!dirB),(a0,!dir0)]
            wt=_occ_wire_from_uses!(m,wuses,caller)
            wid2=kernel.bindings[(-1,wt)]
            ftag=_model_occ_next_tag(kernel,2,-1,caller)
            fid2=only(_model_occ_reserve_shapes(kernel,1,caller))
            _model_occ_bind_shape!(kernel,2,ftag,fid2,
                _ModelOCCShape(2,(kind=:cylinder,
                    geometry=(occ=:cylinder,center=s0,axis=axisc,
                              X=X,Y=Y,radius=spec.r,height=L,
                              angle=θB,pcurves=pc)),
                    UInt64[wid2]))
            kernel.next_shape=fid2+UInt64(1)
            push!(shell_uses,(fid2,false))
        end
        old_ids=UInt64[]
        stack=copy(kernel.shapes[vid].children)
        while !isempty(stack)
            x=pop!(stack); x in old_ids && continue
            push!(old_ids,x); append!(stack,kernel.shapes[x].children)
            hasproperty(kernel.shapes[x].geometry,:uses) &&
                append!(stack,(u for (u,_) in
                               kernel.shapes[x].geometry.uses))
        end
        sus=Tuple{UInt64,Bool}[]
        for (fid,rev) in shell_uses
            if haskey(new_face,fid)
                nf=kernel.bindings[(2,new_face[fid])]
                na=kernel.shapes[nf].geometry.geometry.axis
                push!(sus,(nf,_dot3(na,fnormals[fid])<0))
            elseif haskey(new_face2,fid)
                nf,rv=new_face2[fid]
                push!(sus,(nf,rv))
            else
                push!(sus,(fid,rev))
            end
        end
        shell=_occ_shell_from_uses!(m,sus,caller)
        sid=kernel.bindings[(-2,shell)]
        kernel.shapes[vid]=_ModelOCCShape(3,(kind=:solid,),UInt64[sid])
        _occ_unbind_orphans!(kernel,Set(old_ids))
        _occ_invalidate_volume_tables!(m,vt)
        kernel.changed=true
        push!(out,(3,vt))
    end
    return out
end
