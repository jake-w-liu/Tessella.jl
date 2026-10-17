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
    ids=_model_occ_reserve_shapes(kernel,nsec*(2ne+1),caller)
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

# Project edge `eid` onto base record `face_g`; returns a lin2d or spl2d
# pcurve over the edge's own parameter range.
function _trim_project_pcurve(m,kernel,face_g,eid,caller)
    g=kernel.shapes[eid].geometry.geometry
    etag=_model_occ_shape_tag(kernel,eid,caller)
    t0,t1=g.t0,g.t1
    ns=17
    ts=Float64[t0+(t1-t0)*i/(ns-1) for i in 0:ns-1]
    uvs=Vector{NTuple{2,Float64}}(undef,ns)
    period=face_g.occ===:plane ? 0.0 : _brep_uperiod(face_g)
    prev=nothing
    for (i,u) in enumerate(ts)
        P=_pipe_edge_point(m,g,etag,u)
        pr=_occ_surface_closest(face_g,P,1.0)
        pr===nothing && throw(ArgumentError(
            "$caller: Curve[$etag] does not project onto the surface"))
        uv=(pr[1],pr[2])
        if prev!==nothing && period>0.0
            uv=(uv[1]+period*round((prev[1]-uv[1])/period),uv[2])
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
    lin && return (lin2d=(u0,du),)
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
    knots=Float64[]; cps=NTuple{2,Float64}[]
    for i in 1:n-1
        h=ts[i+1]-ts[i]
        p0=(us[i],vs[i]); p3=(us[i+1],vs[i+1])
        p1=(us[i]+mu[i]*h/3,vs[i]+mv[i]*h/3)
        p2=(us[i+1]-mu[i+1]*h/3,vs[i+1]-mv[i+1]*h/3)
        # Piecewise-Bezier layout: interior knots at full multiplicity p=3
        # (C0 in the B-spline parameterization, tangent-continuous in
        # geometry), junction pole stored once per side.
        i>1 && append!(knots,[ts[i],ts[i],ts[i]])
        i==1 && append!(knots,[ts[1],ts[1],ts[1],ts[1]])
        push!(cps,p0); push!(cps,p1); push!(cps,p2)
        i==n-1 && (push!(cps,p3);
                   append!(knots,[ts[end],ts[end],ts[end],ts[end]]))
    end
    return (knots,3,cps)
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
    ptags=[_model_occ_add_point!(m,x[1],x[2],x[3]) for x in verts]
    edge_of=Dict{NTuple{2,Int},Tuple{UInt64,Int}}()
    face_tags=Int[]
    for i in 1:P
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
    shell=_model_occ_add_surface_loop!(m,face_tags)
    # `OCC_Internals::defeature` mutates the volume in place — rebind the
    # same tag to the healed shell rather than allocating a new volume.
    sid=kernel.bindings[(-2,shell)]
    kernel.shapes[vid]=_ModelOCCShape(3,(kind=:solid,),UInt64[sid])
    return vt
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
