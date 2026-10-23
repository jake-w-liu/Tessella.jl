# `occ.importShapes` / `Merge` — IGES (.igs/.iges) import into the native
# OCC kernel. IGES carries no shared topology: a bounded face
# (143/144/114) or an MSBO (186 → 514 shell → 510 face → 508 loop →
# 504 edge list → 502 vertex list) stores its own boundary curves, and
# two coincident 3-D curves stay distinct entities. Rather than
# duplicating the STEP import machinery, this file materializes each
# IGES graph as a synthesized STEP-entity table — the same
# `EDGE_CURVE`/`EDGE_LOOP`/`ADVANCED_FACE`/`CLOSED_SHELL`/
# `MANIFOLD_SOLID_BREP` shapes — so the importer's edge records
# (line/circle/degenerate/NURBS), surface records (plane/cylinder/cone/
# sphere/torus/tensor_bspline), projected pcurves, and the shell
# validator's two-surface-per-curve certificate all apply unchanged.
#
# Binding scope matches OCC's `importShapes` result: standalone bounded
# faces and bare surfaces bind independently (unwelded, one face per
# IGES face entity), while inside one MSBO coincident vertices and
# shared edge-list members weld — the same contract STEP's shared
# EDGE_CURVEs give the shell validator. Duplicated seam curves inside
# one boundary (the coincident arc pairs OCC writes on a closed face)
# weld through a canonical geometry signature so the wire's seam pair
# shares one edge like STEP's SEAM_CURVE.
#
# Entity coverage is the IGES BRep core: 100 arcs, 110 lines, 126
# rational B-spline curves, 102 composites; 108/190 planes, 198 tori,
# 120 surfaces of revolution (parallel-line → cylinder, coplanar line
# crossing the axis once → cone, axial-plane circle → torus, every
# other generatrix or partial sweep → the exact NURBS revolve), 122
# tabulated cylinders (line → plane strip, arc extruded along its
# plane normal → cylinder, otherwise the exact NURBS extrusion), and
# 128 rational B-spline surfaces; 124 rigid-transform chains;
# 141/142/143/144/114 boundaries; 502/504/508/510/514/186 manifold
# topology; standalone 116 points; and CSG primitives 150/156/158/160
# bound through the OCC primitive adders. Conics (104), copious data
# (106), parametric splines (112), ruled (118) and offset (140)
# surfaces are explicit blockers.

# ── synthesized entity-graph state ────────────────────────────────────

mutable struct _OCCIgesCtx
    caller::String
    entities::Dict{Int,BRep.StepEntity}
    next_id::Int
    types::Dict{Int,Int}
    params::Dict{Int,Vector{Float64}}
    bad::Set{Int}
    transforms::Dict{Int,Int}
    forms::Dict{Int,Int}
    xf_cache::Dict{Int,Any}                  # de => (R,t)|nothing
    geom::Dict{Int,Int}                      # curve de => entity id
    geom_ends::Dict{Int,
        Tuple{NTuple{3,Float64},NTuple{3,Float64}}}
    surf::Dict{Int,Int}                      # surface de => entity id
    cp::Dict{NTuple{3,Float64},Int}          # pos => CARTESIAN_POINT
    dir::Dict{NTuple{3,Float64},Int}         # dir => DIRECTION
    plc::Dict{Any,Int}                       # frame => AXIS2_PLACEMENT_3D
    referenced::Set{Int}                     # des consumed by roots
end

function _OCCIgesCtx(types,params,bad,transforms,forms,caller)
    maxde=isempty(types) ? 0 : maximum(keys(types))
    return _OCCIgesCtx(caller,Dict{Int,BRep.StepEntity}(),maxde+1,
        types,params,bad,transforms,forms,
        Dict{Int,Any}(),Dict{Int,Int}(),
        Dict{Int,Tuple{NTuple{3,Float64},NTuple{3,Float64}}}(),
        Dict{Int,Int}(),Dict{NTuple{3,Float64},Int}(),
        Dict{NTuple{3,Float64},Int}(),Dict{Any,Int}(),Set{Int}())
end

# Per-binding-unit synthesis dedupe: one IGES face binds independently,
# an MSBO's whole graph shares vertices/edges (STEP shared-entity
# semantics). `snappable` marks boundary edges whose junction vertex may
# still be rewritten bitwise onto its successor's start.
struct _OCCIgesSynth
    vert::Dict{NTuple{3,Float64},Int}        # pos => VERTEX_POINT
    ecmap::Dict{Any,Int}                     # edge key => EDGE_CURVE
    snappable::Set{Int}
end
_OCCIgesSynth()=_OCCIgesSynth(Dict{NTuple{3,Float64},Int}(),
                              Dict{Any,Int}(),Set{Int}())

function _occ_iges_ent!(ctx::_OCCIgesCtx,kind::String,args)
    id=ctx.next_id;ctx.next_id+=1
    ent=BRep.StepEntity(id,kind,Any[args...])
    ctx.entities[id]=ent
    return ent
end

function _occ_iges_cp!(ctx::_OCCIgesCtx,p::NTuple{3,Float64})
    id=get(ctx.cp,p,nothing);id!==nothing && return id
    ent=_occ_iges_ent!(ctx,"CARTESIAN_POINT",("",[p[1],p[2],p[3]]))
    ctx.cp[p]=ent.id
    return ent.id
end

function _occ_iges_dir!(ctx::_OCCIgesCtx,d::NTuple{3,Float64})
    id=get(ctx.dir,d,nothing);id!==nothing && return id
    ent=_occ_iges_ent!(ctx,"DIRECTION",("",[d[1],d[2],d[3]]))
    ctx.dir[d]=ent.id
    return ent.id
end

function _occ_iges_plc!(ctx::_OCCIgesCtx,loc,axis,x)
    key=(loc,axis,x)
    id=get(ctx.plc,key,nothing);id!==nothing && return id
    ent=_occ_iges_ent!(ctx,"AXIS2_PLACEMENT_3D",
        ("",(:ref,_occ_iges_cp!(ctx,loc)),
         (:ref,_occ_iges_dir!(ctx,axis)),
         (:ref,_occ_iges_dir!(ctx,x))))
    ctx.plc[key]=ent.id
    return ent.id
end

function _occ_iges_vert!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                        p::NTuple{3,Float64})
    id=get(synth.vert,p,nothing);id!==nothing && return id
    ent=_occ_iges_ent!(ctx,"VERTEX_POINT",
        ("",(:ref,_occ_iges_cp!(ctx,p))))
    synth.vert[p]=ent.id
    return ent.id
end

# ── type-124 transform chains ─────────────────────────────────────────

# An entity's DE field-7 transform applies first; the 124's own
# transform field composes outward (IGES 5.3 §4.21) — the same order
# `BRep._iges_plane_uv_transform` uses for the UV-chart path.
function _occ_iges_xf(ctx::_OCCIgesCtx,de::Int)
    haskey(ctx.xf_cache,de) && return ctx.xf_cache[de]
    matrix=get(ctx.transforms,de,0)
    matrix==0 && (ctx.xf_cache[de]=nothing;return nothing)
    combined=nothing;owner=de;visited=Set{Int}()
    while matrix!=0
        matrix in visited && throw(ArgumentError(
            "$(ctx.caller): IGES transformation matrices contain a " *
            "cycle at #$matrix"))
        push!(visited,matrix);push!(ctx.referenced,matrix)
        transform=BRep._iges_de_transform(ctx.types,ctx.params,ctx.bad,
            ctx.transforms,owner,ctx.caller)
        if combined===nothing
            combined=transform
        else
            rotation,_=transform
            child_rotation,child_translation=combined
            composed_rotation=ntuple(row->ntuple(col->sum(
                rotation[row][k]*child_rotation[k][col] for k in 1:3),
                3),3)
            composed_translation=BRep._iges_xf_apply(transform,
                child_translation,true)
            all(isfinite,(composed_rotation[1]...,composed_rotation[2]...,
                          composed_rotation[3]...,
                          composed_translation...)) || throw(ArgumentError(
                "$(ctx.caller): IGES transform composition is non-finite"))
            combined=(composed_rotation,composed_translation)
        end
        owner=matrix;matrix=get(ctx.transforms,matrix,0)
    end
    ctx.xf_cache[de]=combined
    return combined
end

_occ_iges_pt(ctx,de,p)=BRep._iges_xf_apply(_occ_iges_xf(ctx,de),p,true)
_occ_iges_dv(ctx,de,d)=BRep._iges_xf_apply(_occ_iges_xf(ctx,de),d,false)

# ── curve geometry entities ───────────────────────────────────────────

# Full knot vector → (multiplicities, distinct knots).
function _occ_iges_knot_mults(U::Vector{Float64})
    mults=Int[];knots=Float64[]
    i=1
    while i<=length(U)
        j=i
        while j<length(U) && U[j+1]==U[i];j+=1;end
        push!(mults,j-i+1);push!(knots,U[i]);i=j+1
    end
    return mults,knots
end

# NURBSCurve → B_SPLINE_CURVE_WITH_KNOTS or a rational COMPLEX entity.
function _occ_iges_bspline_curve_ent!(ctx::_OCCIgesCtx,nc::NURBS.NURBSCurve)
    prefs=Any[(:ref,_occ_iges_cp!(ctx,p)) for p in nc.controls]
    mults,knots=_occ_iges_knot_mults(nc.knots)
    if all(w->w==1.0,nc.weights)
        ent=_occ_iges_ent!(ctx,"B_SPLINE_CURVE_WITH_KNOTS",
            ("",nc.degree,prefs,".UNSPECIFIED.",".F.",".F.",
             mults,knots,".UNSPECIFIED."))
        return ent.id
    end
    bc=BRep.StepEntity(0,"B_SPLINE_CURVE",
        Any[nc.degree,prefs,".UNSPECIFIED.",".F.",".F."])
    bk=BRep.StepEntity(0,"B_SPLINE_CURVE_WITH_KNOTS",
        Any[mults,knots,".UNSPECIFIED."])
    br=BRep.StepEntity(0,"RATIONAL_B_SPLINE_CURVE",
        Any[Float64[w for w in nc.weights]])
    ent=_occ_iges_ent!(ctx,"COMPLEX",Any[bc,bk,br])
    return ent.id
end

# NURBSSurface → B_SPLINE_SURFACE_WITH_KNOTS or a rational COMPLEX
# entity. STEP control rows are u-major (outer index = u index).
function _occ_iges_bspline_surf_ent!(ctx::_OCCIgesCtx,
                                    ns::NURBS.NURBSSurface)
    nu,nv=size(ns.controls)
    rows=Any[Any[(:ref,_occ_iges_cp!(ctx,ns.controls[i,j]))
               for j in 1:nv] for i in 1:nu]
    umults,uknots=_occ_iges_knot_mults(ns.knots_u)
    vmults,vknots=_occ_iges_knot_mults(ns.knots_v)
    if all(w->w==1.0,ns.weights)
        ent=_occ_iges_ent!(ctx,"B_SPLINE_SURFACE_WITH_KNOTS",
            ("",ns.degree_u,ns.degree_v,rows,".UNSPECIFIED.",".F.",".F.",
             ".F.",umults,vmults,uknots,vknots,".UNSPECIFIED."))
        return ent.id
    end
    bs=BRep.StepEntity(0,"B_SPLINE_SURFACE",
        Any[ns.degree_u,ns.degree_v,rows,".UNSPECIFIED.",".F.",".F.",
            ".F."])
    bk=BRep.StepEntity(0,"B_SPLINE_SURFACE_WITH_KNOTS",
        Any[umults,vmults,uknots,vknots,".UNSPECIFIED."])
    br=BRep.StepEntity(0,"RATIONAL_B_SPLINE_SURFACE",
        Any[[[ns.weights[i,j] for j in 1:nv] for i in 1:nu]])
    ent=_occ_iges_ent!(ctx,"COMPLEX",Any[bs,bk,br])
    return ent.id
end

# An IGES 100 arc → NURBSCurve (quadratic rational segments, ≤ π/2 each).
function _occ_iges_arc_nurbs(C,X,n,r::Float64,ta::Float64,tb::Float64)
    caller="IGES arc"
    Y=BRep._brep_vcross(n,X)
    span=tb-ta
    span>0 || throw(ArgumentError("$caller: non-positive arc span"))
    nseg=max(1,ceil(Int,span/(π/2)))
    step=span/nseg
    wc=cos(step/2)
    controls=NTuple{3,Float64}[];weights=Float64[];knots=Float64[]
    for k in 0:nseg-1
        a0,a1=ta+step*k,ta+step*(k+1);am=(a0+a1)/2
        p0=(C[1]+r*(cos(a0)*X[1]+sin(a0)*Y[1]),
            C[2]+r*(cos(a0)*X[2]+sin(a0)*Y[2]),
            C[3]+r*(cos(a0)*X[3]+sin(a0)*Y[3]))
        pm=(C[1]+(r/wc)*(cos(am)*X[1]+sin(am)*Y[1]),
            C[2]+(r/wc)*(cos(am)*X[2]+sin(am)*Y[2]),
            C[3]+(r/wc)*(cos(am)*X[3]+sin(am)*Y[3]))
        p1=(C[1]+r*(cos(a1)*X[1]+sin(a1)*Y[1]),
            C[2]+r*(cos(a1)*X[2]+sin(a1)*Y[2]),
            C[3]+r*(cos(a1)*X[3]+sin(a1)*Y[3]))
        if k==0
            push!(controls,p0,pm,p1);append!(weights,(1.0,wc,1.0))
            append!(knots,(a0,a0,a0))
        else
            push!(controls,pm,p1);append!(weights,(wc,1.0))
            append!(knots,(a0,a0))
        end
    end
    append!(knots,(tb,tb,tb))
    return NURBS.NURBSCurve(2,knots,controls,weights)
end

# Simple curve de → (geometry entity id, transformed stored endpoints).
# 102 expands at the edge layer so it is rejected here.
function _occ_iges_curve_geom!(ctx::_OCCIgesCtx,de::Int)
    id=get(ctx.geom,de,nothing)
    id!==nothing && return id,ctx.geom_ends[de]
    push!(ctx.referenced,de)
    t=get(ctx.types,de,0);caller=ctx.caller
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    if t==110
        length(rec)>=7 || throw(ArgumentError(
            "$caller: IGES line #$de is truncated"))
        p1=_occ_iges_pt(ctx,de,(rec[2],rec[3],rec[4]))
        p2=_occ_iges_pt(ctx,de,(rec[5],rec[6],rec[7]))
        delta=BRep._brep_vsub(p2,p1)
        len=sqrt(BRep._brep_vdot(delta,delta))
        dir=len>0 ? BRep._brep_vunit(delta,caller) : (1.0,0.0,0.0)
        vec=_occ_iges_ent!(ctx,"VECTOR",
            ("",(:ref,_occ_iges_dir!(ctx,dir)),max(len,1.0)))
        ent=_occ_iges_ent!(ctx,"LINE",
            ("",(:ref,_occ_iges_cp!(ctx,p1)),(:ref,vec.id)))
        ends=(p1,p2)
    elseif t==100
        length(rec)>=8 || throw(ArgumentError(
            "$caller: IGES arc #$de is truncated"))
        zl,xc,yc,xs,ys,xe,ye=rec[2],rec[3],rec[4],
            rec[5],rec[6],rec[7],rec[8]
        all(isfinite,(zl,xc,yc,xs,ys,xe,ye)) || throw(ArgumentError(
            "$caller: IGES arc #$de is non-finite"))
        C=_occ_iges_pt(ctx,de,(xc,yc,zl))
        n=BRep._brep_vunit(_occ_iges_dv(ctx,de,(0.0,0.0,1.0)),caller)
        S=_occ_iges_pt(ctx,de,(xs,ys,zl))
        E=_occ_iges_pt(ctx,de,(xe,ye,zl))
        rs=BRep._brep_vsub(S,C);r=sqrt(BRep._brep_vdot(rs,rs))
        r>0 || throw(ArgumentError(
            "$caller: IGES arc #$de has zero radius"))
        X=BRep._brep_vunit(rs,caller)
        plc=_occ_iges_plc!(ctx,C,n,X)
        ent=_occ_iges_ent!(ctx,"CIRCLE",("",(:ref,plc),r))
        ends=(S,E)
    elseif t==126
        nc=BRep._iges_nurbs_curve(rec)
        xf=_occ_iges_xf(ctx,de)
        if xf!==nothing
            controls=NTuple{3,Float64}[
                BRep._iges_xf_apply(xf,p,true) for p in nc.controls]
            nc=NURBS.NURBSCurve(nc.degree,nc.knots,controls,nc.weights)
        end
        # V0/V1 follow the control block, before the unit normal.
        n=length(nc.controls);nknots=length(nc.knots)
        i=8+nknots+n+3n
        lo,hi=nc.knots[1],nc.knots[end]
        v0=length(rec)>=i+1 ? rec[i] : lo
        v1=length(rec)>=i+1 ? rec[i+1] : hi
        (isfinite(v0) && isfinite(v1)) || throw(ArgumentError(
            "$caller: IGES B-spline curve #$de has a non-finite range"))
        v0=clamp(v0,lo,hi);v1=clamp(v1,lo,hi)
        v1>v0 || throw(ArgumentError(
            "$caller: IGES B-spline curve #$de has an empty range"))
        ends=(NURBS.nurbs_eval(nc,v0),NURBS.nurbs_eval(nc,v1))
        id=_occ_iges_bspline_curve_ent!(ctx,nc)
        ctx.geom[de]=id;ctx.geom_ends[de]=ends
        return id,ends
    else
        throw(ArgumentError(
            "$caller: IGES curve entity #$de has unsupported type $t"))
    end
    ctx.geom[de]=ent.id;ctx.geom_ends[de]=ends
    return ent.id,ends
end

# A generatrix de → NURBSCurve for the 120/122 tensor fallback.
function _occ_iges_generatrix_nurbs(ctx::_OCCIgesCtx,de::Int)
    t=get(ctx.types,de,0);caller=ctx.caller
    if t==110 || t==126
        t==126 || begin
            rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
            length(rec)>=7 || throw(ArgumentError(
                "$caller: IGES line #$de is truncated"))
            push!(ctx.referenced,de)
            p1=_occ_iges_pt(ctx,de,(rec[2],rec[3],rec[4]))
            p2=_occ_iges_pt(ctx,de,(rec[5],rec[6],rec[7]))
            d=BRep._brep_vsub(p2,p1)
            len=sqrt(BRep._brep_vdot(d,d))
            len>0 || throw(ArgumentError(
                "$caller: IGES line #$de has zero length"))
            return NURBS.NURBSCurve(1,[0.0,0.0,len,len],
                [p1,p2],[1.0,1.0])
        end
        rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
        push!(ctx.referenced,de)
        nc=BRep._iges_nurbs_curve(rec)
        xf=_occ_iges_xf(ctx,de)
        xf===nothing && return nc
        controls=NTuple{3,Float64}[
            BRep._iges_xf_apply(xf,p,true) for p in nc.controls]
        return NURBS.NURBSCurve(nc.degree,nc.knots,controls,nc.weights)
    elseif t==100
        rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
        length(rec)>=8 || throw(ArgumentError(
            "$caller: IGES arc #$de is truncated"))
        push!(ctx.referenced,de)
        zl,xc,yc,xs,ys,xe,ye=rec[2],rec[3],rec[4],
            rec[5],rec[6],rec[7],rec[8]
        C=_occ_iges_pt(ctx,de,(xc,yc,zl))
        n=BRep._brep_vunit(_occ_iges_dv(ctx,de,(0.0,0.0,1.0)),caller)
        S=_occ_iges_pt(ctx,de,(xs,ys,zl))
        E=_occ_iges_pt(ctx,de,(xe,ye,zl))
        rs=BRep._brep_vsub(S,C);r=sqrt(BRep._brep_vdot(rs,rs))
        r>0 || throw(ArgumentError(
            "$caller: IGES arc #$de has zero radius"))
        X=BRep._brep_vunit(rs,caller)
        Y=BRep._brep_vcross(n,X)
        re=BRep._brep_vsub(E,C)
        closed=sqrt(BRep._brep_vdot(BRep._brep_vsub(E,S),
                    BRep._brep_vsub(E,S)))<=1e-9*max(1.0,r)
        ta=0.0
        tb=atan(BRep._brep_vdot(re,Y),BRep._brep_vdot(re,X))
        tb<=ta && (tb+= closed ? 2π : 2π-mod(ta-tb,2π))
        return _occ_iges_arc_nurbs(C,X,n,r,ta,tb)
    end
    throw(ArgumentError(
        "$caller: IGES generatrix entity #$de has unsupported type $t"))
end

# ── surface entities ──────────────────────────────────────────────────

# Exact revolution of a NURBSCurve generatrix about an axis over
# [sa,ea] (quadratic rational arcs in the sweep direction).
function _occ_iges_revolve_surface(gen::NURBS.NURBSCurve,A0,axis,
                                   sa::Float64,ea::Float64,caller)
    span=ea-sa
    (isfinite(span) && span>0) || throw(ArgumentError(
        "$caller: surface of revolution has a non-positive sweep"))
    nseg=max(1,ceil(Int,span/(π/2)))
    step=span/nseg
    wc=cos(step/2)
    nu=2nseg+1;nv=length(gen.controls)
    C=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    W=Matrix{Float64}(undef,nu,nv)
    for j in 1:nv
        P=gen.controls[j];w=gen.weights[j]
        axial=BRep._brep_vdot(BRep._brep_vsub(P,A0),axis)
        foot=(A0[1]+axis[1]*axial,A0[2]+axis[2]*axial,A0[3]+axis[3]*axial)
        rad=BRep._brep_vsub(P,foot)
        rho=sqrt(BRep._brep_vdot(rad,rad))
        udir=rho>0 ? BRep._brep_vunit(rad,caller) : (0.0,0.0,0.0)
        tdir=BRep._brep_vcross(axis,udir)
        for k in 0:nseg-1
            a0=sa+step*k;a1=sa+step*(k+1);am=(a0+a1)/2
            C[2k+1,j]=(foot[1]+rho*(cos(a0)*udir[1]+sin(a0)*tdir[1]),
                       foot[2]+rho*(cos(a0)*udir[2]+sin(a0)*tdir[2]),
                       foot[3]+rho*(cos(a0)*udir[3]+sin(a0)*tdir[3]))
            C[2k+2,j]=(foot[1]+(rho/wc)*(cos(am)*udir[1]+sin(am)*tdir[1]),
                       foot[2]+(rho/wc)*(cos(am)*udir[2]+sin(am)*tdir[2]),
                       foot[3]+(rho/wc)*(cos(am)*udir[3]+sin(am)*tdir[3]))
            W[2k+1,j]=w;W[2k+2,j]=w*wc
        end
        C[nu,j]=(foot[1]+rho*(cos(ea)*udir[1]+sin(ea)*tdir[1]),
                 foot[2]+rho*(cos(ea)*udir[2]+sin(ea)*tdir[2]),
                 foot[3]+rho*(cos(ea)*udir[3]+sin(ea)*tdir[3]))
        W[nu,j]=w
    end
    ku=Float64[]
    append!(ku,(sa,sa,sa))
    for k in 1:nseg-1
        append!(ku,(sa+step*k,sa+step*k))
    end
    append!(ku,(ea,ea,ea))
    return NURBS.NURBSSurface(2,gen.degree,ku,gen.knots,C,W)
end

# Exact extrusion of a NURBSCurve generatrix along displacement `d`.
function _occ_iges_extrude_surface(gen::NURBS.NURBSCurve,d)
    n=length(gen.controls)
    C=Matrix{NTuple{3,Float64}}(undef,n,2)
    W=Matrix{Float64}(undef,n,2)
    for i in 1:n
        P=gen.controls[i]
        C[i,1]=P
        C[i,2]=(P[1]+d[1],P[2]+d[2],P[3]+d[3])
        W[i,1]=gen.weights[i];W[i,2]=gen.weights[i]
    end
    return NURBS.NURBSSurface(gen.degree,1,gen.knots,
        [0.0,0.0,1.0,1.0],C,W)
end

function _occ_iges_plane_ent!(ctx::_OCCIgesCtx,loc,axis,x)
    plc=_occ_iges_plc!(ctx,loc,axis,x)
    return _occ_iges_ent!(ctx,"PLANE",("",(:ref,plc))).id
end

# Point (116) / direction (123) readers — each record's own DE
# transform applies to its stored data.
function _occ_iges_point116(ctx::_OCCIgesCtx,de::Int)
    get(ctx.types,de,0)==116 || throw(ArgumentError(
        "$(ctx.caller): referenced entity #$de is not a Point " *
        "(type 116)"))
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,ctx.caller)
    length(rec)>=4 || throw(ArgumentError(
        "$(ctx.caller): IGES point #$de is truncated"))
    return _occ_iges_pt(ctx,de,(rec[2],rec[3],rec[4]))
end

function _occ_iges_dir123(ctx::_OCCIgesCtx,de::Int)
    get(ctx.types,de,0)==123 || throw(ArgumentError(
        "$(ctx.caller): referenced entity #$de is not a Direction " *
        "(type 123)"))
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,ctx.caller)
    length(rec)>=4 || throw(ArgumentError(
        "$(ctx.caller): IGES direction #$de is truncated"))
    return BRep._brep_vunit(_occ_iges_dv(ctx,de,
        (rec[2],rec[3],rec[4])),ctx.caller)
end

# Line (110) record → (p1,p2) transformed — the 120 axis reader.
function _occ_iges_line_ends(ctx::_OCCIgesCtx,de::Int)
    get(ctx.types,de,0)==110 || throw(ArgumentError(
        "$(ctx.caller): referenced entity #$de is not a Line " *
        "(type 110)"))
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,ctx.caller)
    length(rec)>=7 || throw(ArgumentError(
        "$(ctx.caller): IGES line #$de is truncated"))
    return _occ_iges_pt(ctx,de,(rec[2],rec[3],rec[4])),
           _occ_iges_pt(ctx,de,(rec[5],rec[6],rec[7]))
end

function _occ_iges_surface_ent!(ctx::_OCCIgesCtx,de::Int;
                                bounded_root::Bool=false)
    id=get(ctx.surf,de,nothing);id!==nothing && return id
    push!(ctx.referenced,de)
    t=get(ctx.types,de,0);caller=ctx.caller
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    if t==108
        length(rec)>=5 || throw(ArgumentError(
            "$caller: IGES plane #$de is truncated"))
        nrm=(rec[2],rec[3],rec[4]);dval=rec[5]
        (all(isfinite,nrm) && isfinite(dval)) || throw(ArgumentError(
            "$caller: IGES plane #$de is non-finite"))
        nlen=sqrt(BRep._brep_vdot(nrm,nrm))
        nlen>0 || throw(ArgumentError(
            "$caller: IGES plane #$de has a zero normal"))
        n=(nrm[1]/nlen,nrm[2]/nlen,nrm[3]/nlen)
        loc=(n[1]*dval/nlen,n[2]*dval/nlen,n[3]*dval/nlen)
        xf=_occ_iges_xf(ctx,de)
        loc=BRep._iges_xf_apply(xf,loc,true)
        n=BRep._brep_vunit(BRep._iges_xf_apply(xf,n,false),caller)
        id=_occ_iges_plane_ent!(ctx,loc,n,BRep._brep_default_x(n,caller))
    elseif t==190
        length(rec)>=3 || throw(ArgumentError(
            "$caller: IGES plane surface #$de is truncated"))
        pde=BRep._iges_param_int(rec,2,caller,"plane surface point")
        nde=BRep._iges_param_int(rec,3,caller,"plane surface normal")
        loc=_occ_iges_point116(ctx,pde)
        n=_occ_iges_dir123(ctx,nde)
        x=length(rec)>=4 && rec[4]!=0 ?
            BRep._brep_ortho(_occ_iges_dir123(ctx,
                BRep._iges_param_int(rec,4,caller,
                    "plane surface reference direction")),n,caller) :
            BRep._brep_default_x(n,caller)
        xf=_occ_iges_xf(ctx,de)
        loc=BRep._iges_xf_apply(xf,loc,true)
        n=BRep._brep_vunit(BRep._iges_xf_apply(xf,n,false),caller)
        x=BRep._brep_ortho(BRep._iges_xf_apply(xf,x,false),n,caller)
        id=_occ_iges_plane_ent!(ctx,loc,n,x)
    elseif t==198
        length(rec)>=5 || throw(ArgumentError(
            "$caller: IGES toroidal surface #$de is truncated"))
        pde=BRep._iges_param_int(rec,2,caller,"torus center")
        ade=BRep._iges_param_int(rec,3,caller,"torus axis")
        C=_occ_iges_point116(ctx,pde)
        n=_occ_iges_dir123(ctx,ade)
        R,r=rec[4],rec[5]
        (isfinite(R) && isfinite(r) && R>0 && r>0) || throw(ArgumentError(
            "$caller: toroidal surface #$de radii must be positive"))
        r<R || throw(ArgumentError(
            "$caller: toroidal surface #$de is not a ring torus"))
        x=length(rec)>=6 && rec[6]!=0 ?
            BRep._brep_ortho(_occ_iges_dir123(ctx,
                BRep._iges_param_int(rec,6,caller,
                    "torus reference direction")),n,caller) :
            BRep._brep_default_x(n,caller)
        xf=_occ_iges_xf(ctx,de)
        C=BRep._iges_xf_apply(xf,C,true)
        n=BRep._brep_vunit(BRep._iges_xf_apply(xf,n,false),caller)
        x=BRep._brep_ortho(BRep._iges_xf_apply(xf,x,false),n,caller)
        plc=_occ_iges_plc!(ctx,C,n,x)
        id=_occ_iges_ent!(ctx,"TOROIDAL_SURFACE",
            ("",(:ref,plc),R,r)).id
    elseif t==120
        id=_occ_iges_surf120!(ctx,de,rec;bounded_root=bounded_root)
    elseif t==122
        id=_occ_iges_surf122!(ctx,de,rec;bounded_root=bounded_root)
    elseif t==128
        ns=BRep._iges_nurbs_surface(rec)
        xf=_occ_iges_xf(ctx,de)
        if xf!==nothing
            controls=Matrix{NTuple{3,Float64}}(undef,
                size(ns.controls,1),size(ns.controls,2))
            for idx in eachindex(ns.controls)
                controls[idx]=BRep._iges_xf_apply(xf,
                    ns.controls[idx],true)
            end
            ns=NURBS.NURBSSurface(ns.degree_u,ns.degree_v,
                ns.knots_u,ns.knots_v,controls,ns.weights)
        end
        id=_occ_iges_bspline_surf_ent!(ctx,ns)
    else
        throw(ArgumentError(
            "$caller: IGES surface entity #$de has unsupported " *
            "type $t"))
    end
    ctx.surf[de]=id
    return id
end

# Type 120 surface of revolution: a full 2π sweep classifies into the
# analytic records when the generatrix allows (line parallel to the
# axis → cylinder, coplanar line crossing the axis once → cone, circle
# in an axial plane → torus); every other generatrix or partial sweep
# produces the exact NURBS revolve. A bare root (`bounded_root`) binds
# the analytic record only when it is closed in both directions (a
# torus); every other generatrix takes the tensor revolve whose
# natural domain is the (sweep × generatrix) box.
function _occ_iges_surf120!(ctx::_OCCIgesCtx,de::Int,rec;
                            bounded_root::Bool=false)
    caller=ctx.caller
    length(rec)>=5 || throw(ArgumentError(
        "$caller: surface of revolution #$de is truncated"))
    lde=BRep._iges_param_int(rec,2,caller,"revolution axis")
    gde=BRep._iges_param_int(rec,3,caller,"revolution generatrix")
    sa,ea=rec[4],rec[5]
    (isfinite(sa) && isfinite(ea) && ea>sa) || throw(ArgumentError(
        "$caller: surface of revolution #$de has an invalid sweep"))
    p1,p2=_occ_iges_line_ends(ctx,lde)
    axis=BRep._brep_vunit(BRep._brep_vsub(p2,p1),caller)
    A0=p1
    full=abs(ea-sa-2π)<=1e-9*max(1.0,abs(ea-sa))
    gt=get(ctx.types,gde,0)
    tol=1e-9
    if full && gt==110 && !bounded_root
        g1,g2=_occ_iges_line_ends(ctx,gde)
        gd=BRep._brep_vsub(g2,g1)
        glen=sqrt(BRep._brep_vdot(gd,gd))
        glen>0 || throw(ArgumentError(
            "$caller: revolution #$de generatrix has zero length"))
        ghat=(gd[1]/glen,gd[2]/glen,gd[3]/glen)
        cr=BRep._brep_vcross(ghat,axis)
        cn=sqrt(BRep._brep_vdot(cr,cr))
        rel=BRep._brep_vsub(g1,A0)
        if cn<=tol*glen
            # Generatrix parallel to the axis → cylinder.
            axial=BRep._brep_vdot(rel,axis)
            rad=BRep._brep_vsub(rel,(axis[1]*axial,axis[2]*axial,
                                    axis[3]*axial))
            rho=sqrt(BRep._brep_vdot(rad,rad))
            rho>tol || throw(ArgumentError(
                "$caller: revolution #$de generatrix lies on its axis"))
            loc=(A0[1]+axis[1]*axial,A0[2]+axis[2]*axial,
                 A0[3]+axis[3]*axial)
            X=BRep._brep_vunit(rad,caller)
            plc=_occ_iges_plc!(ctx,loc,axis,X)
            return _occ_iges_ent!(ctx,"CYLINDRICAL_SURFACE",
                ("",(:ref,plc),rho)).id
        end
        ndir=(cr[1]/cn,cr[2]/cn,cr[3]/cn)
        skew=abs(BRep._brep_vdot(rel,ndir))
        if skew<=tol*max(1.0,sqrt(BRep._brep_vdot(rel,rel)))
            # Coplanar non-parallel — a cone whose apex is the
            # generatrix/axis intersection.
            tstar=BRep._brep_vdot(
                BRep._brep_vcross(BRep._brep_vsub(A0,g1),axis),cr)/(cn*cn)
            apex=(g1[1]+ghat[1]*tstar,g1[2]+ghat[2]*tstar,
                  g1[3]+ghat[3]*tstar)
            h1=BRep._brep_vdot(BRep._brep_vsub(g1,apex),axis)
            h2=BRep._brep_vdot(BRep._brep_vsub(g2,apex),axis)
            sgn=h1>tol && h2>tol ? 1.0 :
                h1<-tol && h2<-tol ? -1.0 : 0.0
            if sgn!=0.0
                ax=sgn>0 ? axis : (-axis[1],-axis[2],-axis[3])
                alpha=atan(cn,abs(BRep._brep_vdot(ghat,axis)))
                mid=BRep._brep_vsub(
                    ((g1[1]+g2[1])/2,(g1[2]+g2[2])/2,(g1[3]+g2[3])/2),
                    apex)
                X=BRep._brep_ortho(mid,ax,caller)
                plc=_occ_iges_plc!(ctx,apex,ax,X)
                return _occ_iges_ent!(ctx,"CONICAL_SURFACE",
                    ("",(:ref,plc),0.0,alpha)).id
            end
        end
    elseif full && gt==100
        grec=BRep._iges_rec(ctx.params,ctx.bad,gde,caller)
        length(grec)>=8 || throw(ArgumentError(
            "$caller: generatrix arc #$gde is truncated"))
        zl,xc,yc=grec[2],grec[3],grec[4]
        C=_occ_iges_pt(ctx,gde,(xc,yc,zl))
        an=BRep._brep_vunit(_occ_iges_dv(ctx,gde,(0.0,0.0,1.0)),caller)
        S=_occ_iges_pt(ctx,gde,(grec[5],grec[6],zl))
        E=_occ_iges_pt(ctx,gde,(grec[7],grec[8],zl))
        rad=BRep._brep_vsub(S,C);rl=sqrt(BRep._brep_vdot(rad,rad))
        rl>0 || throw(ArgumentError(
            "$caller: generatrix arc #$gde has zero radius"))
        ctol=1e-7*max(1.0,rl,sqrt(BRep._brep_vdot(C,C)))
        # A bare-root natural wire covers the whole tube, so the torus
        # classification requires a closed generatrix; a partial arc
        # takes the tensor revolve (its own bounded v-domain).
        closed=sqrt(BRep._brep_vdot(BRep._brep_vsub(E,S),
                    BRep._brep_vsub(E,S)))<=ctol
        if (!bounded_root || closed) &&
           abs(BRep._brep_vdot(axis,an))<=1e-7 &&
           abs(BRep._brep_vdot(BRep._brep_vsub(A0,C),an))<=ctol
            # Circle in an axial plane → torus (a partial arc sweeps a
            # toroidal band over the same surface record).
            rel=BRep._brep_vsub(C,A0)
            axial=BRep._brep_vdot(rel,axis)
            radial=BRep._brep_vsub(rel,(axis[1]*axial,axis[2]*axial,
                                        axis[3]*axial))
            R=sqrt(BRep._brep_vdot(radial,radial))
            if R>rl+ctol
                center=(A0[1]+axis[1]*axial,A0[2]+axis[2]*axial,
                        A0[3]+axis[3]*axial)
                X=BRep._brep_vunit(radial,caller)
                plc=_occ_iges_plc!(ctx,center,axis,X)
                return _occ_iges_ent!(ctx,"TOROIDAL_SURFACE",
                    ("",(:ref,plc),R,rl)).id
            end
        end
    end
    gen=_occ_iges_generatrix_nurbs(ctx,gde)
    ns=_occ_iges_revolve_surface(gen,A0,axis,sa,ea,caller)
    return _occ_iges_bspline_surf_ent!(ctx,ns)
end

# Type 122 tabulated cylinder: a line generatrix gives a plane strip,
# an arc extruded along its plane normal gives a cylinder; everything
# else is the exact NURBS extrusion of the generatrix. A bare root
# always binds the extruded tensor — its natural domain is bounded
# where the analytic plane/cylinder records are not.
function _occ_iges_surf122!(ctx::_OCCIgesCtx,de::Int,rec;
                            bounded_root::Bool=false)
    caller=ctx.caller
    length(rec)>=5 || throw(ArgumentError(
        "$caller: tabulated cylinder #$de is truncated"))
    gde=BRep._iges_param_int(rec,2,caller,
        "tabulated cylinder generatrix")
    disp=(rec[3],rec[4],rec[5])
    all(isfinite,disp) || throw(ArgumentError(
        "$caller: tabulated cylinder #$de has a non-finite direction"))
    dlen=sqrt(BRep._brep_vdot(disp,disp))
    dlen>0 || throw(ArgumentError(
        "$caller: tabulated cylinder #$de has a zero direction"))
    dhat=(disp[1]/dlen,disp[2]/dlen,disp[3]/dlen)
    gt=get(ctx.types,gde,0)
    tol=1e-9
    if !bounded_root && gt==110
        g1,g2=_occ_iges_line_ends(ctx,gde)
        gd=BRep._brep_vsub(g2,g1)
        glen=sqrt(BRep._brep_vdot(gd,gd))
        glen>0 || throw(ArgumentError(
            "$caller: tabulated cylinder #$de generatrix has zero " *
            "length"))
        ghat=(gd[1]/glen,gd[2]/glen,gd[3]/glen)
        n=BRep._brep_vcross(ghat,dhat)
        nn=sqrt(BRep._brep_vdot(n,n))
        nn>tol || throw(ArgumentError(
            "$caller: tabulated cylinder #$de sweeps along its own " *
            "generatrix — degenerate surface"))
        n=(n[1]/nn,n[2]/nn,n[3]/nn)
        return _occ_iges_plane_ent!(ctx,g1,n,ghat)
    elseif !bounded_root && gt==100
        grec=BRep._iges_rec(ctx.params,ctx.bad,gde,caller)
        length(grec)>=8 || throw(ArgumentError(
            "$caller: generatrix arc #$gde is truncated"))
        zl,xc,yc=grec[2],grec[3],grec[4]
        C=_occ_iges_pt(ctx,gde,(xc,yc,zl))
        an=BRep._brep_vunit(_occ_iges_dv(ctx,gde,(0.0,0.0,1.0)),caller)
        S=_occ_iges_pt(ctx,gde,(grec[5],grec[6],zl))
        rad=BRep._brep_vsub(S,C);rl=sqrt(BRep._brep_vdot(rad,rad))
        rl>0 || throw(ArgumentError(
            "$caller: generatrix arc #$gde has zero radius"))
        if abs(BRep._brep_vdot(dhat,an))>=1.0-1e-7
            sgn=BRep._brep_vdot(dhat,an)>0 ? 1.0 : -1.0
            ax=(sgn*dhat[1],sgn*dhat[2],sgn*dhat[3])
            X=BRep._brep_vunit(rad,caller)
            plc=_occ_iges_plc!(ctx,C,ax,X)
            return _occ_iges_ent!(ctx,"CYLINDRICAL_SURFACE",
                ("",(:ref,plc),rl)).id
        end
    end
    gen=_occ_iges_generatrix_nurbs(ctx,gde)
    ns=_occ_iges_extrude_surface(gen,disp)
    return _occ_iges_bspline_surf_ent!(ctx,ns)
end

# ── edges and boundary loops ──────────────────────────────────────────

# Canonical coincidence signature so a boundary's duplicated seam
# curves weld into one EDGE_CURVE (STEP's shared-edge semantics).
# Splines key on their entity id only — coincident NURBS pairs are not
# hashed.
function _occ_iges_edge_signature(ctx::_OCCIgesCtx,de::Int,pa,pb)
    t=get(ctx.types,de,0)
    lo,hi=pa<pb ? (pa,pb) : (pb,pa)
    if t==110
        return (:line,lo,hi)
    elseif t==100
        ent=ctx.entities[ctx.geom[de]]
        plc=BRep._deref_entity(ctx.entities,ent.args[2])
        C=BRep._deref_point(ctx.entities,plc.args[2])
        n=BRep._deref_direction(ctx.entities,plc.args[3])
        r=ent.args[3]
        return (:circle,C,r,abs(n[1]),abs(n[2]),abs(n[3]),lo,hi)
    end
    return (:entity,de)
end

# One (edge_curve_id, rev) use for a simple curve de bounded by pa→pb
# (transformed positions). `key` overrides signature dedupe.
function _occ_iges_edge_use!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                            de::Int,pa,pb;key=nothing,rev::Bool=false,
                            snappable::Bool=true)
    ekey=key===nothing ?
        _occ_iges_edge_signature(ctx,de,pa,pb) : key
    ecid=get(synth.ecmap,ekey,nothing)
    if ecid===nothing
        cent,_=_occ_iges_curve_geom!(ctx,de)
        ec=_occ_iges_ent!(ctx,"EDGE_CURVE",
            ("",(:ref,_occ_iges_vert!(ctx,synth,pa)),
             (:ref,_occ_iges_vert!(ctx,synth,pb)),(:ref,cent),".T."))
        ecid=ec.id;synth.ecmap[ekey]=ecid
        snappable && push!(synth.snappable,ecid)
    end
    return (ecid,rev)
end

# A synthetic LINE edge for curves that have no entity (UV-chart
# boundaries mapped through a Form-1 plane frame).
function _occ_iges_line_edge_use!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                  pa,pb,key)
    ecid=get(synth.ecmap,key,nothing)
    if ecid===nothing
        delta=BRep._brep_vsub(pb,pa)
        len=sqrt(BRep._brep_vdot(delta,delta))
        dir=len>0 ? BRep._brep_vunit(delta,ctx.caller) : (1.0,0.0,0.0)
        vec=_occ_iges_ent!(ctx,"VECTOR",
            ("",(:ref,_occ_iges_dir!(ctx,dir)),max(len,1.0)))
        cent=_occ_iges_ent!(ctx,"LINE",
            ("",(:ref,_occ_iges_cp!(ctx,pa)),(:ref,vec.id)))
        ec=_occ_iges_ent!(ctx,"EDGE_CURVE",
            ("",(:ref,_occ_iges_vert!(ctx,synth,pa)),
             (:ref,_occ_iges_vert!(ctx,synth,pb)),(:ref,cent.id),".T."))
        ecid=ec.id;synth.ecmap[key]=ecid
    end
    return (ecid,false)
end

# A vertex list (502) entry's point, under the list's own transform.
function _occ_iges_502_point(ctx::_OCCIgesCtx,vde::Int,idx::Int)
    get(ctx.types,vde,0)==502 || throw(ArgumentError(
        "$(ctx.caller): edge references #$vde which is not a Vertex " *
        "List (type 502)"))
    push!(ctx.referenced,vde)
    rec=BRep._iges_rec(ctx.params,ctx.bad,vde,ctx.caller)
    base=3+3*(idx-1)
    (idx>=1 && length(rec)>=base+2) || throw(ArgumentError(
        "$(ctx.caller): vertex list #$vde is truncated at index $idx"))
    return _occ_iges_pt(ctx,vde,(rec[base],rec[base+1],rec[base+2]))
end

function _occ_iges_curve102_members(ctx::_OCCIgesCtx,de::Int)
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,ctx.caller)
    n=BRep._iges_param_int(rec,2,ctx.caller,
        "composite member count")
    (n>=1 && length(rec)>=2+n) || throw(ArgumentError(
        "$(ctx.caller): composite curve #$de is empty or truncated"))
    return Int[BRep._iges_param_int(rec,2+k,ctx.caller,
        "composite member") for k in 1:n]
end

# A curve's stored endpoints — nested composites resolve through their
# first/last members.
function _occ_iges_curve_member_ends(ctx::_OCCIgesCtx,de::Int)
    get(ctx.types,de,0)==102 ||
        return _occ_iges_curve_geom!(ctx,de)[2]
    members=_occ_iges_curve102_members(ctx,de)
    isempty(members) && throw(ArgumentError(
        "$(ctx.caller): composite curve #$de is empty"))
    return _occ_iges_curve_member_ends(ctx,first(members))[1],
           _occ_iges_curve_member_ends(ctx,last(members))[2]
end

# A 504 edge-list member → uses (a composite member expands).
function _occ_iges_edge504_uses!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                ede::Int,ndx::Int)
    push!(ctx.referenced,ede)
    erec=BRep._iges_rec(ctx.params,ctx.bad,ede,ctx.caller)
    base=3+5*(ndx-1)
    (ndx>=1 && length(erec)>=base+4) || throw(ArgumentError(
        "$(ctx.caller): edge list #$ede is truncated at edge $ndx"))
    cde=BRep._iges_param_int(erec,base,ctx.caller,"edge curve")
    svp=BRep._iges_param_int(erec,base+1,ctx.caller,
        "edge start vertex list")
    svi=BRep._iges_param_int(erec,base+2,ctx.caller,"edge start vertex")
    tvp=BRep._iges_param_int(erec,base+3,ctx.caller,
        "edge end vertex list")
    tvi=BRep._iges_param_int(erec,base+4,ctx.caller,"edge end vertex")
    sv=_occ_iges_502_point(ctx,svp,svi)
    tv=_occ_iges_502_point(ctx,tvp,tvi)
    if get(ctx.types,cde,0)==102
        uses=Tuple{Int,Bool}[]
        members=_occ_iges_curve102_members(ctx,cde)
        prev=nothing
        for (mi,mde) in enumerate(members)
            pa,pb=_occ_iges_curve_member_ends(ctx,mde)
            pa=prev===nothing ? (mi==1 ? sv : pa) : prev
            pb=mi==length(members) ? tv : pb
            push!(uses,_occ_iges_edge_use!(ctx,synth,mde,pa,pb;
                key=(:e504,ede,ndx,mi),snappable=false))
            prev=pb
        end
        return uses
    end
    return [_occ_iges_edge_use!(ctx,synth,cde,sv,tv;
        key=(:e504,ede,ndx),snappable=false)]
end

# UV-only boundary (a 142 with no model curve) → straight model-space
# edges through a Form-1 190 chart — the `_iges_plane_uv_*` semantics of
# the mesh path, mapped to EDGE_CURVEs instead of a vertex ring.
function _occ_iges_uv_wire_uses!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                bde::Int,pde::Int,surf_de::Int)
    caller=ctx.caller
    st=get(ctx.types,surf_de,0)
    st==190 || throw(ArgumentError(
        "$caller: UV-only boundary #$bde needs a Form-1 plane surface " *
        "base — entity #$surf_de has type $st"))
    context=(;forms=ctx.forms,transforms=ctx.transforms)
    origin,x,y=BRep._iges_plane_uv_frame(ctx.types,ctx.params,ctx.bad,
        surf_de,caller,context)
    segments=BRep._iges_plane_uv_segments(ctx.types,ctx.params,ctx.bad,
        pde,caller,context,Set{Int}())
    uses=Tuple{Int,Bool}[]
    prev_end=nothing
    for (a,b) in segments
        pa=(origin[1]+a[1]*x[1]+a[2]*y[1],
            origin[2]+a[1]*x[2]+a[2]*y[2],
            origin[3]+a[1]*x[3]+a[2]*y[3])
        pb=(origin[1]+b[1]*x[1]+b[2]*y[1],
            origin[2]+b[1]*x[2]+b[2]*y[2],
            origin[3]+b[1]*x[3]+b[2]*y[3])
        prev_end!==nothing && (pa=prev_end)
        prev_end=pb
        push!(uses,_occ_iges_line_edge_use!(ctx,synth,pa,pb,
            (:uv,bde,length(uses)+1)))
    end
    return uses
end

# The traversal start/end vertex positions of one (ecid,rev) use.
function _occ_iges_use_ends(ctx,use)
    ecid,rev=use
    ec=ctx.entities[ecid]
    pa=BRep._step_vertex_point(ctx.entities,ec.args[2])
    pb=BRep._step_vertex_point(ctx.entities,ec.args[3])
    return rev ? (pb,pa) : (pa,pb)
end
_occ_iges_use_end(ctx,use)=_occ_iges_use_ends(ctx,use)[2]

# Simple curve de → uses, snapping pa onto `prev_end` when the junction
# is within tolerance (files that round-trip endpoints separately).
function _occ_iges_member_uses!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                               de::Int;prev_end=nothing)
    t=get(ctx.types,de,0);caller=ctx.caller
    if t==102
        uses=Tuple{Int,Bool}[]
        pe=prev_end
        for mde in _occ_iges_curve102_members(ctx,de)
            more=_occ_iges_member_uses!(ctx,synth,mde;prev_end=pe)
            pe=isempty(more) ? pe : _occ_iges_use_end(ctx,more[end])
            append!(uses,more)
        end
        return uses
    elseif t in (110,100,126)
        _,(pa,pb)=_occ_iges_curve_geom!(ctx,de)
        if prev_end!==nothing
            tol=1e-7*max(1.0,BRep._poly_scale(
                NTuple{3,Float64}[prev_end,pa,pb]))
            BRep._pt_close(prev_end,pa,tol) && (pa=prev_end)
        end
        return [_occ_iges_edge_use!(ctx,synth,de,pa,pb)]
    end
    throw(ArgumentError(
        "$caller: IGES boundary references entity #$de — only Line " *
        "(110), Arc (100), B-spline (126) and Composite (102) curves " *
        "bind as OCC edges (saw type $t)"))
end

# Resolve one boundary/loop entity de to ordered (ecid,rev) uses.
function _occ_iges_boundary_uses!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                 de::Int;expected_surface::Int=0,
                                 prev_end=nothing)
    push!(ctx.referenced,de)
    t=get(ctx.types,de,0);caller=ctx.caller
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    if t==142
        length(rec)>=6 || throw(ArgumentError(
            "$caller: curve-on-surface #$de is truncated"))
        surf=BRep._iges_param_int(rec,3,caller,
            "curve-on-surface surface")
        (expected_surface==0 || surf==expected_surface) ||
            throw(ArgumentError(
            "$caller: curve-on-surface #$de binds surface #$surf, " *
            "expected #$expected_surface"))
        pde=BRep._iges_param_int(rec,4,caller,"parametric curve")
        cde=BRep._iges_param_int(rec,5,caller,"model-space curve")
        cde!=0 && return _occ_iges_boundary_uses!(ctx,synth,cde;
            prev_end=prev_end)
        pde!=0 || throw(ArgumentError(
            "$caller: curve-on-surface #$de carries neither a " *
            "model-space nor a parametric curve"))
        return _occ_iges_uv_wire_uses!(ctx,synth,de,pde,surf)
    elseif t==141
        length(rec)>=5 || throw(ArgumentError(
            "$caller: boundary #$de is truncated"))
        if expected_surface!=0
            BRep._iges_param_int(rec,4,caller,"boundary surface")==
                expected_surface || throw(ArgumentError(
                "$caller: boundary #$de does not bind its declared " *
                "face surface"))
        end
        nb=BRep._iges_param_int(rec,5,caller,"boundary curve count")
        uses=Tuple{Int,Bool}[];i=6
        for _ in 1:nb
            length(rec)>=i+2 || throw(ArgumentError(
                "$caller: boundary #$de is truncated"))
            mde=BRep._iges_param_int(rec,i,caller,"boundary curve")
            mflag=BRep._iges_param_int(rec,i+1,caller,
                "boundary curve orientation")
            npc=BRep._iges_param_int(rec,i+2,caller,
                "boundary parameter-curve count")
            npc>=0 || throw(ArgumentError(
                "$caller: boundary #$de has a negative parameter " *
                "count"))
            i+=3+npc
            mflag in (1,2) || throw(ArgumentError(
                "$caller: boundary #$de uses orientation flag $mflag"))
            part=_occ_iges_boundary_uses!(ctx,synth,mde;
                prev_end=isempty(uses) ? prev_end :
                    _occ_iges_use_end(ctx,uses[end]))
            mflag==2 && (part=reverse!([(e,!r) for (e,r) in part]))
            append!(uses,part)
        end
        isempty(uses) && throw(ArgumentError(
            "$caller: boundary #$de is empty"))
        return uses
    elseif t==508
        length(rec)>=3 || throw(ArgumentError(
            "$caller: loop #$de is truncated"))
        n=BRep._iges_param_int(rec,2,caller,"loop edge count")
        uses=Tuple{Int,Bool}[];i=3
        for _ in 1:n
            length(rec)>=i+4 || throw(ArgumentError(
                "$caller: loop #$de is truncated"))
            typ=BRep._iges_param_int(rec,i,caller,"loop entry type")
            ede=BRep._iges_param_int(rec,i+1,caller,"loop edge list")
            ndx=BRep._iges_param_int(rec,i+2,caller,"loop edge index")
            of=BRep._iges_param_int(rec,i+3,caller,
                "loop edge orientation")
            npc=BRep._iges_param_int(rec,i+4,caller,
                "loop parameter-curve count")
            npc>=0 || throw(ArgumentError(
                "$caller: loop #$de has a negative parameter count"))
            i+=5+2*npc
            length(rec)>=i-1 || throw(ArgumentError(
                "$caller: loop #$de is truncated"))
            typ==0 || throw(ArgumentError(
                "$caller: loop #$de uses non-edge entries " *
                "(type $typ)"))
            get(ctx.types,ede,0)==504 || throw(ArgumentError(
                "$caller: loop #$de edge #$ede is not an Edge List " *
                "(type 504)"))
            of in (0,1) || throw(ArgumentError(
                "$caller: loop #$de uses edge orientation $of"))
            for (ecid,srev) in
                    _occ_iges_edge504_uses!(ctx,synth,ede,ndx)
                push!(uses,(ecid,of==0 ? !srev : srev))
            end
        end
        isempty(uses) && throw(ArgumentError(
            "$caller: loop #$de is empty"))
        return uses
    elseif t in (102,110,100,126)
        return _occ_iges_member_uses!(ctx,synth,de;prev_end=prev_end)
    end
    throw(ArgumentError(
        "$caller: IGES boundary entity #$de has unsupported type $t"))
end

# Snap the closing junction of a uses list onto its first vertex when
# they coincide within tolerance but differ bitwise — the ring path's
# exact chain check compares tuples with `==`.
function _occ_iges_snap_close!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                             uses::Vector{Tuple{Int,Bool}})
    n=length(uses)
    n<2 && return uses
    first_start=_occ_iges_use_ends(ctx,uses[1])[1]
    last_end=_occ_iges_use_ends(ctx,uses[end])[2]
    last_end==first_start && return uses
    tol=1e-7*max(1.0,BRep._poly_scale(
        NTuple{3,Float64}[first_start,last_end]))
    BRep._pt_close(last_end,first_start,tol) || return uses
    ecid,rev=uses[end]
    ecid in synth.snappable || return uses
    ec=ctx.entities[ecid]
    pa=BRep._step_vertex_point(ctx.entities,ec.args[2])
    pb=BRep._step_vertex_point(ctx.entities,ec.args[3])
    # The traversal end is `pb` for a forward use, `pa` reversed.
    na,nb=rev ? (first_start,pb) : (pa,first_start)
    ec2=_occ_iges_ent!(ctx,"EDGE_CURVE",
        ("",(:ref,_occ_iges_vert!(ctx,synth,na)),
         (:ref,_occ_iges_vert!(ctx,synth,nb)),
         (:ref,ec.args[4][2]),".T."))
    uses[end]=(ec2.id,rev)
    return uses
end

# Ring positions for a uses list — the containment certificate
# `BRep._brep_face_outer_wire_order` consumes.
function _occ_iges_uses_ring(ctx,uses)
    ring=NTuple{3,Float64}[]
    for use in uses
        push!(ring,_occ_iges_use_ends(ctx,use)[1])
    end
    return ring
end

# Synthesize FACE_(OUTER_)BOUND/EDGE_LOOP/ORIENTED_EDGE entities around
# a uses list; returns the bound entity id.
function _occ_iges_loop_bound!(ctx::_OCCIgesCtx,uses;outer::Bool)
    oes=Any[]
    for (ecid,rev) in uses
        oe=_occ_iges_ent!(ctx,"ORIENTED_EDGE",
            ("",nothing,nothing,(:ref,ecid),rev ? ".F." : ".T."))
        push!(oes,(:ref,oe.id))
    end
    loop=_occ_iges_ent!(ctx,"EDGE_LOOP",("",oes))
    return _occ_iges_ent!(ctx,
        outer ? "FACE_OUTER_BOUND" : "FACE_BOUND",
        ("",(:ref,loop.id),".T.")).id
end

# Synthesize one ADVANCED_FACE over `surf_ent` from per-loop uses
# lists; `outer` is the designated outer loop index (0 → resolve by
# containment). `sense` is the ADVANCED_FACE same_sense flag.
function _occ_iges_synth_face!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                               surf_ent::Int,loop_uses;outer::Int=1,
                               sense::Bool=true)
    for uses in loop_uses
        _occ_iges_snap_close!(ctx,synth,uses)
    end
    if outer==0 && length(loop_uses)>1
        rings=[_occ_iges_uses_ring(ctx,u) for u in loop_uses]
        ordered=BRep._brep_face_outer_wire_order(rings,ctx.caller)
        idx=findfirst(==(first(ordered)),rings)
        outer=idx===nothing ? 1 : idx
    end
    bounds=Any[]
    for (i,uses) in enumerate(loop_uses)
        isempty(uses) && throw(ArgumentError(
            "$(ctx.caller): IGES face has an empty boundary loop"))
        push!(bounds,(:ref,_occ_iges_loop_bound!(ctx,uses;
            outer=i==outer)))
    end
    fe=_occ_iges_ent!(ctx,"ADVANCED_FACE",
        ("",bounds,(:ref,surf_ent),sense ? ".T." : ".F."))
    return fe
end

# ── natural-domain wires ──────────────────────────────────────────────

# Iso-curve extraction on a knot-plane boundary: collapsing the control
# net in homogeneous space at `value` is exact.
function _occ_iges_iso_curve(ns::NURBS.NURBSSurface,fixed::Symbol,
                             value::Float64)
    nu,nv=size(ns.controls)
    if fixed===:u
        basis=[NURBS.bspline_basis(i,ns.degree_u,ns.knots_u,value)
               for i in 1:nu]
        controls=NTuple{3,Float64}[];weights=Float64[]
        for j in 1:nv
            hp=(0.0,0.0,0.0);hw=0.0
            for i in 1:nu
                b=basis[i]*ns.weights[i,j]
                hp=(hp[1]+b*ns.controls[i,j][1],
                    hp[2]+b*ns.controls[i,j][2],
                    hp[3]+b*ns.controls[i,j][3])
                hw+=b
            end
            hw>0 || throw(ArgumentError(
                "IGES import: degenerate iso-curve weight"))
            push!(controls,(hp[1]/hw,hp[2]/hw,hp[3]/hw))
            push!(weights,hw)
        end
        return NURBS.NURBSCurve(ns.degree_v,ns.knots_v,controls,weights)
    else
        basis=[NURBS.bspline_basis(j,ns.degree_v,ns.knots_v,value)
               for j in 1:nv]
        controls=NTuple{3,Float64}[];weights=Float64[]
        for i in 1:nu
            hp=(0.0,0.0,0.0);hw=0.0
            for j in 1:nv
                b=basis[j]*ns.weights[i,j]
                hp=(hp[1]+b*ns.controls[i,j][1],
                    hp[2]+b*ns.controls[i,j][2],
                    hp[3]+b*ns.controls[i,j][3])
                hw+=b
            end
            hw>0 || throw(ArgumentError(
                "IGES import: degenerate iso-curve weight"))
            push!(controls,(hp[1]/hw,hp[2]/hw,hp[3]/hw))
            push!(weights,hw)
        end
        return NURBS.NURBSCurve(ns.degree_u,ns.knots_u,controls,weights)
    end
end

# Two NURBS curves coincide geometrically when sampled positions agree
# — welds a periodic patch's iso-curve pair into one seam edge.
function _occ_iges_curves_coincide(a::NURBS.NURBSCurve,b::NURBS.NURBSCurve)
    tol=1e-9*max(1.0,
        maximum(sqrt(BRep._brep_vdot(p,p)) for p in a.controls))
    for k in 0:8
        ta=a.knots[1]+(a.knots[end]-a.knots[1])*k/8
        tb=b.knots[1]+(b.knots[end]-b.knots[1])*k/8
        d=BRep._brep_vsub(NURBS.nurbs_eval(a,ta),
                          NURBS.nurbs_eval(b,tb))
        BRep._brep_vdot(d,d)<=tol*tol || return false
    end
    return true
end

# Natural (u,v)-box boundary of a tensor_bspline surface entity —
# four iso-curves wound CCW in parameter space.
function _occ_iges_natural_uses!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                surf_ent_id::Int)
    ent=ctx.entities[surf_ent_id]
    ns=BRep._step_nurbs_surface(ent,ctx.entities)
    ns===nothing && throw(ArgumentError(
        "$(ctx.caller): natural-boundary synthesis needs a B-spline " *
        "surface (saw $(ent.kind) #$surf_ent_id)"))
    u0,u1=ns.knots_u[1],ns.knots_u[end]
    v0,v1=ns.knots_v[1],ns.knots_v[end]
    function viso(value,pa,pb)
        nc=_occ_iges_iso_curve(ns,:v,value)
        cid=_occ_iges_bspline_curve_ent!(ctx,nc)
        _occ_iges_ent!(ctx,"EDGE_CURVE",
            ("",(:ref,_occ_iges_vert!(ctx,synth,pa)),
             (:ref,_occ_iges_vert!(ctx,synth,pb)),
             (:ref,cid),".T.")).id
    end
    function uiso(value,pa,pb)
        nc=_occ_iges_iso_curve(ns,:u,value)
        cid=_occ_iges_bspline_curve_ent!(ctx,nc)
        _occ_iges_ent!(ctx,"EDGE_CURVE",
            ("",(:ref,_occ_iges_vert!(ctx,synth,pa)),
             (:ref,_occ_iges_vert!(ctx,synth,pb)),
             (:ref,cid),".T.")).id
    end
    p(u,v)=NURBS.nurbs_eval(ns,u,v)
    e_v0=viso(v0,p(u0,v0),p(u1,v0))
    e_v1=viso(v1,p(u0,v1),p(u1,v1))
    e_u0=uiso(u0,p(u0,v0),p(u0,v1))
    e_u1=uiso(u1,p(u1,v0),p(u1,v1))
    # Periodic seam welding — coincident iso-curve pairs collapse to
    # one shared edge used in both directions.
    ncv(id)=BRep._step_nurbs_curve(
        ctx.entities[ctx.entities[id].args[4][2]],ctx.entities)
    cv0,cv1,cu0,cu1=ncv(e_v0),ncv(e_v1),ncv(e_u0),ncv(e_u1)
    cu0!==nothing && cu1!==nothing &&
        _occ_iges_curves_coincide(cu0,cu1) && (e_u1=e_u0)
    cv0!==nothing && cv1!==nothing &&
        _occ_iges_curves_coincide(cv0,cv1) && (e_v1=e_v0)
    # CCW in (u,v): v0 forward, u1 forward, v1 reversed, u0 reversed.
    return Tuple{Int,Bool}[(e_v0,false),(e_u1,false),
                            (e_v1,true),(e_u0,true)]
end

# Natural wire of a full torus: the rim (outer equator) and meridian
# (u=0 tube circle), each closed circle used in both directions — the
# same seam-pair structure the importer's periodic-wire handling
# expects.
function _occ_iges_torus_natural!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                  surf_ent_id::Int)
    caller=ctx.caller
    ent=ctx.entities[surf_ent_id]
    refs=[a for a in ent.args if a isa Tuple && a[1]===:ref]
    plc=BRep._deref_entity(ctx.entities,refs[1])
    loc,n,X,Y=_occ_import_placement(ctx.entities,plc,caller)
    R=Float64(ent.args[3]);r=Float64(ent.args[4])
    # Outer equator circle: center on the axis at the torus center,
    # radius R+r, in the X-Y plane.
    vp=(loc[1]+(R+r)*X[1],loc[2]+(R+r)*X[2],loc[3]+(R+r)*X[3])
    rim_ent=_occ_iges_ent!(ctx,"CIRCLE",
        ("",(:ref,_occ_iges_plc!(ctx,loc,n,X)),R+r))
    vid=_occ_iges_vert!(ctx,synth,vp)
    rim=_occ_iges_ent!(ctx,"EDGE_CURVE",
        ("",(:ref,vid),(:ref,vid),(:ref,rim_ent.id),".T."))
    # u=0 meridian: circle of radius r about the tube center
    # loc+R·X, in the plane spanned by X and n (normal n×X).
    merC=(loc[1]+R*X[1],loc[2]+R*X[2],loc[3]+R*X[3])
    tdir=BRep._brep_vcross(n,X)
    mer_ent=_occ_iges_ent!(ctx,"CIRCLE",
        ("",(:ref,_occ_iges_plc!(ctx,merC,tdir,X)),r))
    vm=(merC[1]+r*X[1],merC[2]+r*X[2],merC[3]+r*X[3])
    vid2=_occ_iges_vert!(ctx,synth,vm)
    mer=_occ_iges_ent!(ctx,"EDGE_CURVE",
        ("",(:ref,vid2),(:ref,vid2),(:ref,mer_ent.id),".T."))
    # Wire: rim fwd, meridian fwd, rim rev, meridian rev.
    return Tuple{Int,Bool}[(rim.id,false),(mer.id,false),
                            (rim.id,true),(mer.id,true)]
end

# ── face/solid entry resolution ───────────────────────────────────────

# A 144 trimmed surface → (surf_de, loop_uses | :natural, outer_index).
function _occ_iges_face144_loops!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                 de::Int)
    caller=ctx.caller
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    length(rec)>=5 || throw(ArgumentError(
        "$caller: trimmed surface #$de is truncated"))
    surf_de=BRep._iges_param_int(rec,2,caller,"trimmed surface base")
    n1=BRep._iges_param_int(rec,3,caller,"outer-boundary flag")
    n2=BRep._iges_param_int(rec,4,caller,"inner-boundary count")
    n2>=0 || throw(ArgumentError(
        "$caller: trimmed surface #$de has a negative inner count"))
    length(rec)>=5+n2 || throw(ArgumentError(
        "$caller: trimmed surface #$de is truncated"))
    get(ctx.types,surf_de,0)!=0 || throw(ArgumentError(
        "$caller: trimmed surface #$de references missing surface " *
        "#$surf_de"))
    n1 in (0,1) || throw(ArgumentError(
        "$caller: trimmed surface #$de outer-boundary flag $n1 " *
        "is invalid"))
    if n1==0
        n2==0 || throw(ArgumentError(
            "$caller: trimmed surface #$de has inner wires without " *
            "an outer wire"))
        return surf_de,:natural,1
    end
    pto=BRep._iges_param_int(rec,5,caller,"outer boundary")
    pto!=0 || throw(ArgumentError(
        "$caller: trimmed surface #$de has no outer boundary"))
    loops=Vector{Tuple{Int,Bool}}[
        _occ_iges_boundary_uses!(ctx,synth,pto;
            expected_surface=surf_de)]
    for k in 1:n2
        pti=BRep._iges_param_int(rec,5+k,caller,"inner boundary")
        pti!=0 || throw(ArgumentError(
            "$caller: trimmed surface #$de has a null inner boundary"))
        push!(loops,_occ_iges_boundary_uses!(ctx,synth,pti;
            expected_surface=surf_de))
    end
    return surf_de,loops,1
end

# A 143/114 bounded surface → (surf_de, loop_uses, outer_index).
function _occ_iges_bounded_loops!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                 de::Int)
    caller=ctx.caller
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    length(rec)>=4 || throw(ArgumentError(
        "$caller: bounded surface #$de is truncated"))
    surf_de=BRep._iges_param_int(rec,3,caller,"bounded surface base")
    nb=BRep._iges_param_int(rec,4,caller,
        "bounded surface boundary count")
    (nb>=1 && length(rec)>=4+nb) || throw(ArgumentError(
        "$caller: bounded surface #$de has no boundaries"))
    get(ctx.types,surf_de,0)!=0 || throw(ArgumentError(
        "$caller: bounded surface #$de references missing surface " *
        "#$surf_de"))
    loops=Vector{Tuple{Int,Bool}}[]
    for k in 1:nb
        bde=BRep._iges_param_int(rec,4+k,caller,
            "bounded surface boundary")
        push!(loops,_occ_iges_boundary_uses!(ctx,synth,bde;
            expected_surface=surf_de))
    end
    return surf_de,loops,length(loops)>1 ? 0 : 1
end

# A 510 face → (surf_de, loop_uses, outer_index).
function _occ_iges_face510_loops!(ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                                 de::Int)
    caller=ctx.caller
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    length(rec)>=4 || throw(ArgumentError(
        "$caller: face #$de is truncated"))
    surf_de=BRep._iges_param_int(rec,2,caller,"face surface")
    n=BRep._iges_param_int(rec,3,caller,"face loop count")
    of=BRep._iges_param_int(rec,4,caller,"face outer-loop flag")
    n>=1 || throw(ArgumentError("$caller: face #$de has no loops"))
    n==1 || of in (0,1) || throw(ArgumentError(
        "$caller: face #$de has an invalid outer-loop flag $of"))
    length(rec)>=4+n || throw(ArgumentError(
        "$caller: face #$de is truncated"))
    get(ctx.types,surf_de,0)!=0 || throw(ArgumentError(
        "$caller: face #$de references missing surface #$surf_de"))
    loops=Vector{Tuple{Int,Bool}}[]
    for k in 1:n
        lde=BRep._iges_param_int(rec,4+k,caller,"face loop")
        push!(loops,_occ_iges_boundary_uses!(ctx,synth,lde;
            expected_surface=surf_de))
    end
    # Flag 1 → the first loop is the designated outer; flag 0 →
    # resolve by containment (matches `BRep._iges_face510`).
    outer=n==1 ? 1 : (of==1 ? 1 : 0)
    return surf_de,loops,outer
end

# Bind one IGES face through a synthesized ADVANCED_FACE.
function _occ_iges_bind_face!(m,ctx::_OCCIgesCtx,synth::_OCCIgesSynth,
                            d::_OCCImportDedupe,st::_OCCImportStructural,
                            surf_de::Int,loops,outer::Int,sense::Bool,
                            curvekinds)
    caller=ctx.caller
    surf_ent=_occ_iges_surface_ent!(ctx,surf_de)
    if loops===:natural
        ent=ctx.entities[surf_ent]
        uses=if ent.kind in ("B_SPLINE_SURFACE_WITH_KNOTS","COMPLEX")
            _occ_iges_natural_uses!(ctx,synth,surf_ent)
        elseif ent.kind=="TOROIDAL_SURFACE"
            _occ_iges_torus_natural!(ctx,synth,surf_ent)
        else
            throw(ArgumentError(
                "$caller: IGES face on surface #$surf_de keeps its " *
                "natural boundary — unsupported on a $(ent.kind) base"))
        end
        loops=[uses];outer=1
    end
    fe=_occ_iges_synth_face!(ctx,synth,surf_ent,loops;outer=outer,
        sense=sense)
    ftag,_=_occ_import_step_face(m,d,st,m.occ_kernel,ctx.entities,
        ctx.entities[fe.id],curvekinds,caller)
    return ftag
end

# One MSBO (186) → synthesized BREP → `_occ_import_step_solid` →
# [(3,tag)]. All shells/faces share one synth table so coincident
# vertices and edge-list members weld like shared STEP topology.
function _occ_iges_bind_msbo!(m,ctx::_OCCIgesCtx,de::Int)
    caller=ctx.caller
    push!(ctx.referenced,de)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    nv=BRep._iges_param_int(rec,4,caller,"MSBO void-shell count")
    (nv>=0 && length(rec)>=4+2nv) || throw(ArgumentError(
        "$caller: MSBO #$de is truncated"))
    synth=_OCCIgesSynth()
    shells=Tuple{Int,Int}[]
    push!(shells,(BRep._iges_param_int(rec,2,caller,"MSBO shell"),
                  BRep._iges_param_int(rec,3,caller,
                      "MSBO shell orientation")))
    for k in 1:nv
        push!(shells,(BRep._iges_param_int(rec,4+2k-1,caller,
                          "MSBO void shell"),
                      BRep._iges_param_int(rec,4+2k,caller,
                          "MSBO void orientation")))
    end
    shell_ents=Int[]
    for (sde,sflag) in shells
        push!(ctx.referenced,sde)
        get(ctx.types,sde,0)==514 || throw(ArgumentError(
            "$caller: MSBO #$de references #$sde which is not a " *
            "Shell (type 514)"))
        sflag in (0,1) || throw(ArgumentError(
            "$caller: MSBO #$de uses shell orientation $sflag"))
        srec=BRep._iges_rec(ctx.params,ctx.bad,sde,caller)
        nf=BRep._iges_param_int(srec,2,caller,"shell face count")
        (nf>=1 && length(srec)>=2+2nf) || throw(ArgumentError(
            "$caller: shell #$sde is empty or truncated"))
        face_refs=Any[]
        for k in 1:nf
            fde=BRep._iges_param_int(srec,2+2k-1,caller,"shell face")
            fof=BRep._iges_param_int(srec,2+2k,caller,
                "shell face orientation")
            get(ctx.types,fde,0)==510 || throw(ArgumentError(
                "$caller: shell #$sde member #$fde is not a Face " *
                "(type 510)"))
            fof in (0,1) || throw(ArgumentError(
                "$caller: shell #$sde uses face orientation $fof"))
            flip=xor(fof==1,sflag==1)
            surf_de,loops,outer=_occ_iges_face510_loops!(ctx,synth,fde)
            surf_ent=_occ_iges_surface_ent!(ctx,surf_de)
            fe=_occ_iges_synth_face!(ctx,synth,surf_ent,loops;
                outer=outer,sense=!flip)
            push!(face_refs,(:ref,fe.id))
        end
        se=_occ_iges_ent!(ctx,"CLOSED_SHELL",("",face_refs))
        push!(shell_ents,se.id)
    end
    brep=length(shell_ents)==1 ?
        _occ_iges_ent!(ctx,"MANIFOLD_SOLID_BREP",
            ("",(:ref,shell_ents[1]))) :
        _occ_iges_ent!(ctx,"BREP_WITH_VOIDS",
            ("",(:ref,shell_ents[1]),
             Any[(:ref,i) for i in shell_ents[2:end]]))
    return _occ_import_step_solid(m,brep,ctx.entities,caller)
end

# One standalone IGES face (144/143/114 or unreferenced 510) →
# [(2,tag)] — an independent binding unit matching OCC's
# one-face-per-FACE transfer.
function _occ_iges_bind_face_root!(m,ctx::_OCCIgesCtx,de::Int)
    t=get(ctx.types,de,0)
    synth=_OCCIgesSynth()
    if t==144
        surf_de,loops,outer=_occ_iges_face144_loops!(ctx,synth,de)
    elseif t==510
        surf_de,loops,outer=_occ_iges_face510_loops!(ctx,synth,de)
    else
        surf_de,loops,outer=_occ_iges_bounded_loops!(ctx,synth,de)
    end
    d=_OCCImportDedupe();st=_OCCImportStructural()
    curvekinds=Set{String}()
    ftag=_occ_iges_bind_face!(m,ctx,synth,d,st,surf_de,loops,outer,
        true,curvekinds)
    return [(2,ftag)]
end

# A bare surface root (120/122/128/190/198 not consumed by any face)
# → its natural-domain face where bounded. Unbounded surfaces (bare
# planes, cylinders/cones without a bounded domain) are explicit
# blockers.
function _occ_iges_bind_surface_root!(m,ctx::_OCCIgesCtx,de::Int)
    caller=ctx.caller
    t=get(ctx.types,de,0)
    surf_ent=_occ_iges_surface_ent!(ctx,de;bounded_root=true)
    ent=ctx.entities[surf_ent]
    synth=_OCCIgesSynth()
    uses=if ent.kind in ("B_SPLINE_SURFACE_WITH_KNOTS","COMPLEX")
        _occ_iges_natural_uses!(ctx,synth,surf_ent)
    elseif ent.kind=="TOROIDAL_SURFACE"
        _occ_iges_torus_natural!(ctx,synth,surf_ent)
    else
        throw(ArgumentError(
            "$caller: bare IGES surface #$de of type $t has no " *
            "bounded natural domain — only B-spline and toroidal " *
            "surfaces bind without a bounded-face wrapper"))
    end
    fe=_occ_iges_synth_face!(ctx,synth,surf_ent,[uses];outer=1,
        sense=true)
    d=_OCCImportDedupe();st=_OCCImportStructural()
    curvekinds=Set{String}()
    ftag,_=_occ_import_step_face(m,d,st,m.occ_kernel,ctx.entities,
        ctx.entities[fe.id],curvekinds,caller)
    return [(2,ftag)]
end

# Standalone curve roots (110/100/126/102 not referenced by any face
# or surface) → one edge each; composites bind member-by-member.
function _occ_iges_bind_curve_root!(m,ctx::_OCCIgesCtx,de::Int)
    t=get(ctx.types,de,0)
    if t==102
        out=Tuple{Int,Int}[]
        for mde in _occ_iges_curve102_members(ctx,de)
            append!(out,_occ_iges_bind_curve_root!(m,ctx,mde))
        end
        return out
    end
    synth=_OCCIgesSynth()
    d=_OCCImportDedupe();st=_OCCImportStructural()
    curvekinds=Set{String}()
    _,(pa,pb)=_occ_iges_curve_geom!(ctx,de)
    (ecid,_)=_occ_iges_edge_use!(ctx,synth,de,pa,pb)
    eid,=_occ_import_step_edge!(m,d,st,ctx.entities,(:ref,ecid),
        curvekinds,ctx.caller)
    tag=_model_occ_shape_tag(m.occ_kernel,eid,ctx.caller)
    return [(1,tag)]
end

# Standalone vertex root (unreferenced 116) → one point.
function _occ_iges_bind_point_root!(m,ctx::_OCCIgesCtx,de::Int)
    p=_occ_iges_point116(ctx,de)
    return [(0,_model_occ_add_point!(m,p[1],p[2],p[3]))]
end

# CSG primitive roots → OCC volume primitives.
function _occ_iges_bind_primitive_root!(m,ctx::_OCCIgesCtx,de::Int)
    caller=ctx.caller
    push!(ctx.referenced,de)
    t=get(ctx.types,de,0)
    rec=BRep._iges_rec(ctx.params,ctx.bad,de,caller)
    if t==150
        length(rec)>=7 || throw(ArgumentError(
            "$caller: IGES block #$de is truncated"))
        L,W,H=rec[2],rec[3],rec[4]
        x,y,z=rec[5],rec[6],rec[7]
        (L>0 && W>0 && H>0) || throw(ArgumentError(
            "$caller: block #$de extents must be positive"))
        corner=_occ_iges_pt(ctx,de,(x,y,z))
        ext=(L,W,H)
        xf=_occ_iges_xf(ctx,de)
        if xf!==nothing
            R,_=xf
            # An axis-aligned OCC box survives only signed-permutation
            # rotations: local axis i maps to model axis perm[i] with
            # sign s_i; extents permute, and negative images shift the
            # stored min corner along that model axis.
            perm=zeros(Int,3);sgn=zeros(Int,3)
            for i in 1:3
                col=(R[1][i],R[2][i],R[3][i])
                nz=[j for j in 1:3 if abs(col[j])>0.5]
                (length(nz)==1 &&
                 abs(abs(col[nz[1]])-1.0)<=1e-9) || throw(ArgumentError(
                    "$caller: block #$de rotation is not " *
                    "axis-aligned"))
                perm[i]=nz[1];sgn[i]=col[nz[1]]>0 ? 1 : -1
            end
            sort(perm)==[1,2,3] || throw(ArgumentError(
                "$caller: block #$de rotation is not axis-aligned"))
            new_ext=ntuple(j->ext[findfirst(==(j),perm)],3)
            for i in 1:3
                sgn[i]<0 || continue
                corner=ntuple(k->corner[k]-
                    (k==perm[i] ? ext[i] : 0.0),3)
            end
            ext=new_ext
        end
        tag=_model_occ_add_box!(m,corner[1],corner[2],corner[3],
            ext[1],ext[2],ext[3])
        return [(3,tag)]
    elseif t==158
        length(rec)>=5 || throw(ArgumentError(
            "$caller: IGES sphere #$de is truncated"))
        r=rec[2];r>0 || throw(ArgumentError(
            "$caller: sphere #$de radius must be positive"))
        c=_occ_iges_pt(ctx,de,(rec[3],rec[4],rec[5]))
        return [(3,_model_occ_add_sphere!(m,c[1],c[2],c[3],r))]
    elseif t==156
        if length(rec)>=10
            h,r1,r2=rec[2],rec[3],rec[4]
            b=_occ_iges_pt(ctx,de,(rec[5],rec[6],rec[7]))
            dir=_occ_iges_dv(ctx,de,(rec[8],rec[9],rec[10]))
        elseif length(rec)>=8
            # Short cylinder form: radius, height, base center, axis.
            r1=r2=rec[2];h=rec[3]
            b=_occ_iges_pt(ctx,de,(rec[4],rec[5],rec[6]))
            dir=_occ_iges_dv(ctx,de,(rec[7],rec[8],
                length(rec)>=9 ? rec[9] : 1.0))
        else
            throw(ArgumentError("$caller: IGES cone #$de is truncated"))
        end
        (h>0 && r1>=0 && r2>=0 && (r1>0 || r2>0)) || throw(ArgumentError(
            "$caller: cone #$de has invalid dimensions"))
        dl=sqrt(BRep._brep_vdot(dir,dir))
        dl>0 || throw(ArgumentError(
            "$caller: cone #$de has a zero axis"))
        if r1==r2
            tag=_model_occ_add_cylinder!(m,b[1],b[2],b[3],
                dir[1]/dl*h,dir[2]/dl*h,dir[3]/dl*h,r1)
        else
            tag=_model_occ_add_cone!(m,b[1],b[2],b[3],
                dir[1]/dl*h,dir[2]/dl*h,dir[3]/dl*h,r1,r2)
        end
        return [(3,tag)]
    elseif t==160
        length(rec)>=9 || throw(ArgumentError(
            "$caller: IGES torus #$de is truncated"))
        R,r=rec[2],rec[3]
        (R>0 && r>0 && r<R) || throw(ArgumentError(
            "$caller: torus #$de is not a ring torus"))
        c=_occ_iges_pt(ctx,de,(rec[4],rec[5],rec[6]))
        z=_occ_iges_dv(ctx,de,(rec[7],rec[8],rec[9]))
        tag=_model_occ_add_torus!(m,c[1],c[2],c[3],R,r;
            zAxis=Float64[z[1],z[2],z[3]])
        return [(3,tag)]
    end
    throw(ArgumentError(
        "$caller: IGES primitive entity #$de has unsupported type $t"))
end

# ── entry ─────────────────────────────────────────────────────────────

const _OCC_IGES_FACE_TYPES=(144,143,114)
const _OCC_IGES_SURFACE_ROOTS=(128,120,122,190,198)
const _OCC_IGES_CURVE_ROOTS=(110,100,126,102)
const _OCC_IGES_PRIMITIVE_ROOTS=(150,156,158,160)

"""
Import an IGES file's topology into the OCC kernel — the `.igs`/`.iges`
counterpart of [`_model_occ_import_shapes`](@ref) STEP import. MSBO
(186) solids bind as volumes, bounded faces (143/144/114/510) as
surfaces, bare supported surfaces as natural-domain faces, unreferenced
curves/points as edges/points, and CSG primitives (150/156/158/160) as
primitive volumes. `highest_dim_only` keeps only the top-dimensional
products, matching OCC's `importShapes` contract.
"""
function _model_occ_import_iges(m::GeoModel,path;highest_dim_only=true,
                                caller="OCC.importShapes")
    isfile(path) || throw(ArgumentError("$caller: missing file $path"))
    source=read(path,String)
    types,params,bad,transforms,frame_context=
        BRep._iges_entity_records(source;with_frame_context=true)
    isempty(types) && throw(ArgumentError(
        "$caller: $path has no IGES directory entries"))
    kernel=m.occ_kernel
    for d in 0:3
        kernel.maxima[d+3]=max(kernel.maxima[d+3],m.next_tag[d+1])
    end
    isempty(m.loops) ||
        (kernel.maxima[2]=max(kernel.maxima[2],
                              maximum(keys(m.loops))))
    isempty(m.surface_loops) ||
        (kernel.maxima[1]=max(kernel.maxima[1],
                              maximum(keys(m.surface_loops))))
    ctx=_OCCIgesCtx(types,params,bad,transforms,frame_context.forms,
        caller)
    out=Tuple{Int,Int,Int}[]   # (root de, dim, tag)
    des=sort!(collect(keys(types)))
    for de in des
        get(types,de,0)==186 || continue
        for (dim,tag) in _occ_iges_bind_msbo!(m,ctx,de)
            push!(out,(de,dim,tag))
        end
    end
    # 510 faces consumed by an MSBO-referenced shell are not roots.
    consumed510=Set{Int}()
    for de in ctx.referenced
        get(types,de,0)==514 || continue
        srec=get(params,de,nothing)
        srec===nothing && continue
        nf=Int(srec[2])
        for k in 1:nf
            push!(consumed510,Int(srec[2+2k-1]))
        end
    end
    for de in des
        t=get(types,de,0)
        (t in _OCC_IGES_FACE_TYPES ||
         (t==510 && !(de in consumed510))) || continue
        for (dim,tag) in _occ_iges_bind_face_root!(m,ctx,de)
            push!(out,(de,dim,tag))
        end
    end
    for de in des
        t=get(types,de,0)
        (t in _OCC_IGES_SURFACE_ROOTS && !(de in ctx.referenced)) ||
            continue
        for (dim,tag) in _occ_iges_bind_surface_root!(m,ctx,de)
            push!(out,(de,dim,tag))
        end
    end
    for de in des
        t=get(types,de,0)
        (t in _OCC_IGES_CURVE_ROOTS && !(de in ctx.referenced)) ||
            continue
        for (dim,tag) in _occ_iges_bind_curve_root!(m,ctx,de)
            push!(out,(de,dim,tag))
        end
    end
    for de in des
        t=get(types,de,0)
        (t in _OCC_IGES_PRIMITIVE_ROOTS && !(de in ctx.referenced)) ||
            continue
        for (dim,tag) in _occ_iges_bind_primitive_root!(m,ctx,de)
            push!(out,(de,dim,tag))
        end
    end
    for de in des
        get(types,de,0)==116 && !(de in ctx.referenced) || continue
        for (dim,tag) in _occ_iges_bind_point_root!(m,ctx,de)
            push!(out,(de,dim,tag))
        end
    end
    isempty(out) && throw(ArgumentError(
        "$caller: no supported IGES geometry in $path (saw entity " *
        "types $(sort!(unique(collect(values(types))))))"))
    result=Tuple{Int,Int}[(dim,tag) for (_,dim,tag) in
                          sort!(out;by=r->r[1])]
    if highest_dim_only
        top=maximum(dim for (dim,_) in result)
        result=Tuple{Int,Int}[(dim,tag) for (dim,tag) in result
                              if dim==top]
    end
    return result
end
