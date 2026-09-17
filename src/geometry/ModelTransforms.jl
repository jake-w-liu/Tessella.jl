# General affine entity transforms, `Duplicata` copies, and post-transform
# coherence for the native `GeoModel`.
#
# Semantics mirror the Gmsh built-in kernel's `ApplicationOnShapes` +
# `ReplaceAllDuplicates` pipeline (`Geo.cpp`): a transform mutates the shared
# vertices of each listed entity in place, recursively through boundary
# topology, with every vertex moved at most once per statement; a coincident
# entity merge then restores coherence. `Duplicata` deep-copies entities at
# their current positions (fresh tags allocated through Gmsh's shared
# `Geometry.OldNewReg` counter) and is merge-free on its own — the copies are
# normally collapsed or spread by the enclosing transform statement.

const _COHERENCE_RTOL = 1e-8

# A fused affine map plus the exact sequence of 4×4 Gmsh `ApplicationOnShapes`
# matrices that realize it. `linear`/`offset` describe the map semantically
# (encoding representability checks); `steps` are the per-pass matrices Gmsh
# applies to each owned vertex, in order — `Rotate` is three passes
# (translate by -origin, rotate, translate by +origin), so coordinates match
# Gmsh's rounding bit-for-bit. `occ` carries the gp_Trsf/gp_GTrsf equivalent —
# (linear, loc) applied as M·p + loc — used for OCC-geometry records and
# OCC-owned vertices, where `occ.rotate`'s Rodrigues matrix rounds
# differently from Gmsh's SetRotationMatrix.
struct _AffineTransform
    linear::NTuple{9,Float64}   # row-major 3×3
    offset::NTuple{3,Float64}
    steps::Vector{NTuple{16,Float64}}  # row-major 4×4, homogeneous last column
    occ::@NamedTuple{linear::NTuple{9,Float64},loc::NTuple{3,Float64}}
end

# `gp_XYZ::Multiplied(gp_Mat)` row — x·M(i,1) unrounded into the rounded
# y·M(i,2), then z·M(i,3) unrounded into the partial (verified against the
# compiled libTKMath `fmadd` chain).
@inline _occ_matvec(L::NTuple{9,Float64}, v::NTuple{3,Float64}) =
    (fma(L[3],v[3],fma(L[1],v[1],L[2]*v[2])),
     fma(L[6],v[3],fma(L[4],v[1],L[5]*v[2])),
     fma(L[9],v[3],fma(L[7],v[1],L[8]*v[2])))

# `gp_Trsf::Transforms`: (M·p)·scale + loc. The uniform-scale factor is
# already folded into `occ.linear` for gp_GTrsf paths (dilate/affine), so the
# application is one M·p followed by the componentwise loc add.
@inline _occ_trsf_apply(t::_AffineTransform, p::NTuple{3,Float64}) =
    _occ_matvec(t.occ.linear,p) .+ t.occ.loc

@inline _occ_linear_apply(t::_AffineTransform, v::NTuple{3,Float64}) =
    _occ_matvec(t.occ.linear,v)

# `gp_Mat::SetRotation` — Rodrigues form I + sin·K + (1−cos)·K² with K the
# cross matrix of the normalized axis, built through OCCT's own pass order
# (K·s, += I, += K²·(1−cos)). The diagonal of K² contracts as
# −fma(X,X,Y·Y), and the libm calls go through the platform libm shims.
function _occ_rotation_matrix(axis::NTuple{3,Float64}, angle::Float64)
    mod=sqrt(fma(axis[3],axis[3],fma(axis[1],axis[1],axis[2]*axis[2])))
    A,B,C=axis[1]/mod,axis[2]/mod,axis[3]/mod
    s,c=_gm_sincos(angle)
    K=((0.0,-C,B),(C,0.0,-A),(-B,A,0.0))
    K2=((-fma(C,C,B*B),A*B,A*C),
        (A*B,-fma(A,A,C*C),B*C),
        (A*C,B*C,-fma(A,A,B*B)))
    omc=1.0-c
    return ntuple(9) do k
        i=(k-1)÷3+1; j=(k-1)%3+1
        (K[i][j]*s + (i==j ? 1.0 : 0.0)) + K2[i][j]*omc
    end
end

# `vecmat4x4` from Geo.cpp: res[i] = Σ_j mat[i][j]·vec[j] with vec=(x,y,z,1),
# accumulated in order. The shipped Gmsh binary rounds each product-add
# separately here (the `+=` loop does not fuse), while `SetRotationMatrix`'s
# matrix products do fuse — verified against the 4.15.2 binary.
@inline function _gmsh_matvec4x4(mat::NTuple{16,Float64}, p::NTuple{3,Float64})
    v=(p[1],p[2],p[3],1.0)
    return ntuple(3) do i
        acc=0.0
        for j in 1:4
            acc+=mat[4*(i-1)+j]*v[j]
        end
        acc
    end
end

@inline function _affine_apply_steps(t::_AffineTransform, p::NTuple{3,Float64})
    for step in t.steps
        p=_gmsh_matvec4x4(step,p)
    end
    return p
end

@inline function _affine_apply(t::_AffineTransform, p::NTuple{3,Float64})
    (a,b,c,d,e,f,g,h,i)=t.linear
    return (fma(c,p[3],fma(a,p[1],b*p[2]))+t.offset[1],
            fma(f,p[3],fma(d,p[1],e*p[2]))+t.offset[2],
            fma(i,p[3],fma(g,p[1],h*p[2]))+t.offset[3])
end

@inline function _linear_apply(t::_AffineTransform, v::NTuple{3,Float64})
    (a,b,c,d,e,f,g,h,i)=t.linear
    return (fma(c,v[3],fma(a,v[1],b*v[2])),
            fma(f,v[3],fma(d,v[1],e*v[2])),
            fma(i,v[3],fma(g,v[1],h*v[2])))
end

@inline _dot(a::NTuple{3,Float64},b::NTuple{3,Float64}) =
    fma(a[3],b[3],fma(a[1],b[1],a[2]*b[2]))

@inline function _entity_label(dim::Int)
    dim==0 && return "Point"
    return _model_periodic_entity_label(dim)
end

@inline _gmsh_translation_step(d::NTuple{3,Float64}) =
    (1.0,0.0,0.0,d[1], 0.0,1.0,0.0,d[2], 0.0,0.0,1.0,d[3], 0.0,0.0,0.0,1.0)

function _affine_translation(delta, caller)
    d=_finite_vector3(delta,caller,"translation delta")
    return _AffineTransform((1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),d,
                            [_gmsh_translation_step(d)],
                            (linear=(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),
                             loc=d))
end

function _affine_dilation(center, scale, caller)
    c=_finite_vector3(center,caller,"dilation center")
    scales=scale isa Real ?
        let s=_finite_scalar(scale,caller,"dilation scale"); (s,s,s) end :
        _finite_vector3(scale,caller,"dilation scales")
    (sx,sy,sz)=scales
    # SetDilatationMatrix: mat[i][3] = T[i]*(1-scale_i).
    step=(sx,0.0,0.0,c[1]*(1.0-sx),
          0.0,sy,0.0,c[2]*(1.0-sy),
          0.0,0.0,sz,c[3]*(1.0-sz),
          0.0,0.0,0.0,1.0)
    return _AffineTransform((sx,0.0,0.0,0.0,sy,0.0,0.0,0.0,sz),
                            (c[1]*(1-sx),c[2]*(1-sy),c[3]*(1-sz)),[step],
                            (linear=(sx,0.0,0.0,0.0,sy,0.0,0.0,0.0,sz),
                             loc=(c[1]*(1-sx),c[2]*(1-sy),c[3]*(1-sz))))
end

# `norme`/`prodve` from Gmsh's Numeric.h: normalization multiplies by the
# reciprocal; the cross product uses Gmsh's component order. `fma` replicates
# the fused multiply-add the compiled sum-of-products expressions produce
# (`norm3`/`prosca`: `fma(c,c,fma(a,a,b*b))`; `prodve`: the first product of
# each component fused, the second rounded).
@inline function _gmsh_norme!(v::Vector{Float64})
    mod=sqrt(fma(v[1],v[1],fma(v[2],v[2],v[3]*v[3])))
    if mod!=0.0
        inv=1.0/mod
        v[1]*=inv;v[2]*=inv;v[3]*=inv
    end
    return v
end

@inline _gmsh_prodve(a::Vector{Float64},b::Vector{Float64}) =
    [fma(a[2],b[3],-(a[3]*b[2])), fma(-a[1],b[3],a[3]*b[1]),
     fma(a[1],b[2],-(a[2]*b[1]))]

# SetRotationMatrix from Geo.cpp: orthonormal basis (axis, t1, t2) built by
# Gmsh's GramSchmidt, an in-basis rotation about the axis, then
# plan'·rot·plan with Gmsh's triple-loop order.
function _gmsh_rotation_matrix(axis::NTuple{3,Float64}, angle::Float64)
    axe=[axis[1],axis[2],axis[3]]
    if axe[1]!=0.0
        t1=[0.0,1.0,0.0];t2=[0.0,0.0,1.0]
    elseif axe[2]!=0.0
        t1=[1.0,0.0,0.0];t2=[0.0,0.0,1.0]
    else
        t1=[1.0,0.0,0.0];t2=[0.0,1.0,0.0]
    end
    # GramSchmidt(axe, t1, t2): v1=axis, v2=t1, v3=t2.
    _gmsh_norme!(axe)
    t1=_gmsh_norme!(_gmsh_prodve(t2,axe))
    t2=_gmsh_norme!(_gmsh_prodve(axe,t1))
    plan=(axe,t1,t2)   # rows
    s,c=_gm_sincos(angle)
    rot=((1.0,0.0,0.0),(0.0,c,-s),(0.0,s,c))
    interm=ntuple(3) do i
        ntuple(3) do j
            acc=0.0
            for k in 1:3
                acc=fma(plan[k][i],rot[k][j],acc)   # invplan[i][k] = plan[k][i]
            end
            acc
        end
    end
    mat=ntuple(3) do i
        ntuple(3) do j
            acc=0.0
            for k in 1:3
                acc=fma(interm[i][k],plan[k][j],acc)
            end
            acc
        end
    end
    return (mat[1][1],mat[1][2],mat[1][3],0.0,
            mat[2][1],mat[2][2],mat[2][3],0.0,
            mat[3][1],mat[3][2],mat[3][3],0.0,
            0.0,0.0,0.0,1.0)
end

function _affine_rotation(axis, origin, angle, caller)
    a=_finite_vector3(axis,caller,"rotation axis")
    o=_finite_vector3(origin,caller,"rotation origin")
    θ=_finite_scalar(angle,caller,"rotation angle")
    n=sqrt(a[1]^2+a[2]^2+a[3]^2)
    n>0 || throw(ArgumentError("$caller: rotation axis must be nonzero"))
    rstep=_gmsh_rotation_matrix(a,θ)
    L=(rstep[1],rstep[2],rstep[3],rstep[5],rstep[6],rstep[7],
       rstep[9],rstep[10],rstep[11])
    Ro=_gmsh_matvec4x4(rstep,o)
    # gp_Trsf::SetRotation: loc = A1.Location − M·A1.Location with the
    # Rodrigues-built M.
    Mo=_occ_rotation_matrix(a,θ)
    Mo_o=_occ_matvec(Mo,o)
    return _AffineTransform(L,(o[1]-Ro[1],o[2]-Ro[2],o[3]-Ro[3]),
                            [_gmsh_translation_step((-o[1],-o[2],-o[3])),
                             rstep,
                             _gmsh_translation_step(o)],
                            (linear=Mo,
                             loc=(o[1]-Mo_o[1],o[2]-Mo_o[2],o[3]-Mo_o[3])))
end

# Reflection across the plane `A*x + B*y + C*z + D = 0` (Gmsh `Symmetry`). A
# degenerate zero normal takes Gmsh's `p -> 1e-12` floor and acts as identity.
function _affine_symmetry(a, b, c, d, caller)
    A=_finite_scalar(a,caller,"symmetry plane coefficient A")
    B=_finite_scalar(b,caller,"symmetry plane coefficient B")
    C=_finite_scalar(c,caller,"symmetry plane coefficient C")
    D=_finite_scalar(d,caller,"symmetry plane coefficient D")
    p=fma(C,C,fma(A,A,B*B))
    p==0.0 && (p=1e-12)
    F=-2.0/p
    step=(fma(A*A,F,1.0), A*B*F, A*C*F, A*D*F,
          A*B*F, fma(B*B,F,1.0), B*C*F, B*D*F,
          A*C*F, B*C*F, fma(C*C,F,1.0), C*D*F,
          0.0,0.0,0.0,1.0)
    L=(step[1],step[2],step[3],step[5],step[6],step[7],step[9],step[10],step[11])
    # occ.symmetrize routes through gp_GTrsf with the identical vectorial
    # part and translation column.
    return _AffineTransform(L,(step[4],step[8],step[12]),[step],
                            (linear=L,loc=(step[4],step[8],step[12])))
end

"""
    transform_entities!(model, transform::_AffineTransform, entities; caller) -> entities

Apply an affine transform in place to the shared vertices of each listed
`(dimension, tag)` entity. Boundary topology is traversed recursively —
surfaces contribute their loops' curves' endpoints and volumes their shells' —
and each point is moved at most once per call, matching Gmsh's per-statement
transformed-point set. Compact volume encodings are updated when the transform
is representable (isotropic for spheres; axis-perpendicular isotropic for
cylinders and cones; recomputed corners for materialized `add_box!` volumes)
and are discarded with the shell topology preserved otherwise. Boolean operand
snapshots are transformed with the result; snapshots recorded for other
Booleans stay untouched, matching Gmsh's operation-time snapshot semantics. A
coincident-entity coherence merge runs afterwards. The call is atomic:
non-representable transforms and non-finite results leave the model unchanged.
"""
function transform_entities!(m::GeoModel, t::_AffineTransform,
                             entities::AbstractVector{<:Tuple{Integer,Integer}};
                             caller::AbstractString="transform_entities!")
    normalized=NTuple{2,Int}[]
    for (dim,tag) in entities
        d=_dimension(dim,caller); tg=_tag(tag,caller,d)
        haskey(m.discrete,(d,tg)) && throw(ArgumentError(
            "$caller: discrete $(_entity_label(d))[$tg] is not transformable"))
        push!(normalized,(d,tg))
    end
    move=Int[]; seen=Set{Int}(); volumes=Int[]
    for (d,tg) in normalized
        for p in _entity_point_tags(m,d,tg,caller)
            p in seen || (push!(seen,p); push!(move,p))
        end
        d==3 && push!(volumes,tg)
    end
    plans=[_plan_volume_transform(m,tg,t,caller) for tg in volumes]
    geometry_plans=_plan_occ_geometry_transforms(m,normalized,t,caller)
    # Vertices on OCC-geometried edges transform through the gp_Trsf path
    # (one M·p + loc pass) rather than Gmsh's ApplicationOnShapes steps.
    occ_pts=Set{Int}()
    for (ct,g) in m.curve_geometry
        hasproperty(g,:occ) || continue
        a,b=m.curves[ct]
        push!(occ_pts,a); push!(occ_pts,b)
    end
    coords=[(tag,_finite_result(
                tag in occ_pts ? _occ_trsf_apply(t,m.points[tag]) :
                                 _affine_apply_steps(t,m.points[tag]),caller))
            for tag in move]
    for (tag,p) in coords
        m.points[tag]=p
    end
    for plan in plans
        _apply_volume_plan!(m,plan)
    end
    _apply_occ_geometry_plans!(m,geometry_plans)
    _reconcile_box_encodings!(m,seen)
    _reconcile_curved_encodings!(m,seen)
    coherence!(m)
    return normalized
end

# Curves and surfaces owned by one entity — the closure that carries OCC
# geometry records under a transform.
function _entity_closure_geometry(m::GeoModel, d::Int, tag::Int)
    curves=Int[]; surfaces=Int[]
    if d==1
        push!(curves,tag)
    elseif d==2
        push!(surfaces,tag)
        for l in m.surfaces[tag], c in m.loops[l]
            push!(curves,abs(c))
        end
    elseif d==3
        for shell in m.volumes[tag], signed_surface in m.surface_loops[shell]
            surface=abs(signed_surface)
            push!(surfaces,surface)
            for l in m.surfaces[surface], c in m.loops[l]
                push!(curves,abs(c))
            end
        end
    end
    return curves,surfaces
end

# Stage OCC curve/surface geometry record updates for the transformed entities.
# An OCC circle can only stay circular under an in-plane similarity; anything
# else is rejected atomically like the compact volume encodings.
function _plan_occ_geometry_transforms(m::GeoModel, normalized, t,
                                       caller::AbstractString)
    curve_plans=Tuple{Int,NamedTuple}[]
    surface_plans=Tuple{Int,NamedTuple}[]
    seen_curves=Set{Int}(); seen_surfaces=Set{Int}()
    for (d,tag) in normalized
        curves,surfaces=_entity_closure_geometry(m,d,tag)
        for curve in curves
            curve in seen_curves && continue
            push!(seen_curves,curve)
            g=_occ_geometry(m,curve)
            g===nothing && continue
            what="Curve[$curve]"
            push!(curve_plans,(curve,
                _transform_occ_curve_record(m,curve,g,t,caller,what)))
        end
        for surface in surfaces
            surface in seen_surfaces && continue
            push!(seen_surfaces,surface)
            g=get(m.surface_geometry,surface,nothing)
            (g===nothing || !hasproperty(g,:occ)) && continue
            push!(surface_plans,(surface,
                _transform_occ_surface_record(g,t,caller,"Surface[$surface]")))
        end
    end
    return (curve_plans,surface_plans)
end

# The transformed OCC curve record for `curve`, or `g` unchanged when the
# record carries no coordinates (`:degenerate` stores only its parameter
# range).
function _transform_occ_curve_record(m::GeoModel, curve::Int, g,
                                     t::_AffineTransform, caller, what)
    if g.occ===:circle
        return _transform_occ_circle(g,t,caller,what)
    elseif g.occ===:line
        # Geom_TrimmedCurve::Transformed keeps its [t0,t1] bounds; the basis
        # gp_Lin moves as Location' = M·P + loc, Direction' = M·D — the
        # direction is not renormalized, so arc-length parameters stay
        # consistent under the same stretch as the endpoints.
        dir=_finite_result(_occ_linear_apply(t,g.dir),caller)
        return (occ=:line,origin=_occ_trsf_apply(t,g.origin),dir=dir,
                t0=g.t0,t1=g.t1)
    end
    return g
end

function _transform_occ_surface_record(g, t::_AffineTransform, caller, what)
    if g.occ===:sphere
        s2=_linear_isotropy(t.linear)
        s2===nothing && throw(ArgumentError(
            "$caller: transform is not representable on $what — " *
            "it requires an isotropic linear part"))
        axis=_occ_linear_apply(t,g.axis)
        X=_occ_linear_apply(t,g.X)
        Y=_occ_linear_apply(t,g.Y)
        return (occ=:sphere,
                center=_finite_result(
                    _occ_trsf_apply(t,g.center),caller),
                radius=g.radius*sqrt(s2),
                axis=_finite_result(axis./sqrt(_dot(axis,axis)),caller),
                X=_finite_result(X./sqrt(_dot(X,X)),caller),
                Y=_finite_result(Y./sqrt(_dot(Y,Y)),caller),
                pcurves=g.pcurves)
    elseif g.occ in (:cylinder,:cone)
        axis_vector=g.axis .* g.height
        tr=_transform_axis_encoding(
            t,g.center,axis_vector,caller,what)
        n=(tr.axis[1]/tr.height,tr.axis[2]/tr.height,
           tr.axis[3]/tr.height)
        X=_occ_linear_apply(t,g.X)
        X=_finite_result(X./sqrt(_dot(X,X)),caller)
        Y=_occ_linear_apply(t,g.Y)
        Y=_finite_result(Y./sqrt(_dot(Y,Y)),caller)
        return g.occ===:cylinder ?
            (occ=:cylinder,center=tr.center,axis=n,X=X,Y=Y,
             radius=g.radius*tr.perp,height=tr.height,
             pcurves=g.pcurves) :
            (occ=:cone,center=tr.center,axis=n,X=X,Y=Y,
             r1=g.r1*tr.perp,r2=g.r2*tr.perp,height=tr.height,
             pcurves=g.pcurves)
    elseif g.occ===:torus
        # A torus only survives a similarity: the equator and the tube
        # must stay circular, which needs an isotropic linear part.
        s2=_linear_isotropy(t.linear)
        s2===nothing && throw(ArgumentError(
            "$caller: transform is not representable on $what — " *
            "it requires an isotropic linear part"))
        axis=_occ_linear_apply(t,g.axis)
        # Under a reflection the derived in-plane direction flips
        # handedness; storing −T·axis keeps Y′ = T·Y so
        # p′(u,v) = T·p(u,−v) — the [0,angle] sweep still runs from
        # the start radial to the end radial (v traverses the closed
        # tube, so negating it describes the same surface).
        _linear_det(t.linear)<0 && (axis=.-axis)
        X=_occ_linear_apply(t,g.X)
        Y=_occ_linear_apply(t,g.Y)
        return (occ=:torus,
                center=_finite_result(
                    _occ_trsf_apply(t,g.center),caller),
                axis=_finite_result(axis./sqrt(_dot(axis,axis)),caller),
                X=_finite_result(X./sqrt(_dot(X,X)),caller),
                Y=_finite_result(Y./sqrt(_dot(Y,Y)),caller),
                r1=g.r1*sqrt(s2),r2=g.r2*sqrt(s2),angle=g.angle,
                pcurves=g.pcurves)
    elseif g.occ===:plane
        # `gp_Ax3::Transformed` — every direction takes the raw vectorial
        # part (unrenormalized: a `gp_GTrsf` scale stretches the frame and
        # the parametrization absorbs it, so the pcurves — and therefore
        # `uvb` — carry over verbatim).
        return (occ=:plane,
                center=_finite_result(
                    _occ_trsf_apply(t,g.center),caller),
                axis=_occ_linear_apply(t,g.axis),
                X=_occ_linear_apply(t,g.X),
                Y=_occ_linear_apply(t,g.Y),
                reversed=g.reversed,pcurves=g.pcurves,uvb=g.uvb)
    end
    return g
end

function _apply_occ_geometry_plans!(m::GeoModel, plans)
    curve_plans,surface_plans=plans
    for (tag,record) in curve_plans
        m.curve_geometry[tag]=record
    end
    for (tag,record) in surface_plans
        m.surface_geometry[tag]=record
    end
    return nothing
end

# Transform an OCC circle record: the frame must remain orthonormal, which is
# exactly the condition that the transformed curve is still a circle.
function _transform_occ_circle(g, t::_AffineTransform, caller, what)
    center=_finite_result(_occ_trsf_apply(t,g.center),caller)
    X=_occ_linear_apply(t,g.X)
    Y=_occ_linear_apply(t,_occ_circle_y(g))
    sx=sqrt(_dot(X,X)); sy=sqrt(_dot(Y,Y))
    (sx>0 && sy>0) || throw(ArgumentError(
        "$caller: transform collapses $what"))
    tol=1e-12*max(1.0,sx*sx,sy*sy)
    (abs(sx*sx-sy*sy)<=tol && abs(_dot(X,Y))<=tol) || throw(ArgumentError(
        "$caller: transform is not representable on $what — " *
        "it would warp the circle into an ellipse"))
    n=_occ_linear_apply(t,g.n)
    nl=sqrt(_dot(n,n))
    nl>0 || throw(ArgumentError("$caller: transform collapses $what"))
    # Under a reflection the circle's derived Y = n×X flips handedness. A
    # closed circle keeps +T·n (its endpoints coincide and cap normals must
    # still match the solid axis), the negated range carrying the reversed
    # traversal p'(t) = T·p(−t) — hence −T·Y; a trimmed arc must keep t0 on
    # its start vertex, so it stores −T·n and keeps its range — p'(t) =
    # T·p(t) in the flipped frame, hence +T·Y.
    t0,t1=g.t0,g.t1
    if _linear_det(t.linear)<0
        if g.t1-g.t0<_OCC_TWO_PI
            n=.-n
        else
            t0,t1=-g.t1,-g.t0
            Y=.-Y
        end
    end
    return (occ=:circle,center=center,
            n=(n[1]/nl,n[2]/nl,n[3]/nl),
            X=(X[1]/sx,X[2]/sx,X[3]/sx),
            Y=(Y[1]/sy,Y[2]/sy,Y[3]/sy),r=g.r*sx,
            t0=t0,t1=t1)
end

# After independent sub-entity moves, keep a materialized curved primitive's
# encoding only while its boundary entities still satisfy it — the same
# reconcile-or-drop contract as `_reconcile_box_encodings!`.
function _reconcile_curved_encodings!(m::GeoModel, moved::Set{Int})
    for dict in (:cylinders,:spheres,:cones)
        store=getfield(m,dict)
        for tag in collect(keys(store))
            isempty(m.volumes[tag]) && continue
            owned=_model_volume_owned_points(m,tag)
            any(p->p in moved,owned) || continue
            _materialized_curved_consistent(m,tag) && continue
            delete!(store,tag)
        end
    end
    return nothing
end

# Point tags owned by one entity — the transform's atomic unit. Parent entities
# recurse through boundary topology; duplicate endpoints are deduplicated by the
# caller's per-statement set.
function _entity_point_tags(m::GeoModel, d::Int, tag::Int, caller)
    if d==0
        haskey(m.points,tag) || throw(ArgumentError("$caller: unknown Point[$tag]"))
        return (tag,)
    elseif d==1
        haskey(m.curves,tag) || throw(ArgumentError("$caller: unknown Curve[$tag]"))
        a,b=m.curves[tag]
        cps=get(m.curve_control_points,tag,Int[])
        return (a,b,cps...)
    elseif d==2
        haskey(m.surfaces,tag) || throw(ArgumentError(
            "$caller: unknown Surface[$tag]"))
        pts=Int[]
        for l in m.surfaces[tag]
            for c in m.loops[l]
                a,b=m.curves[abs(c)]
                push!(pts,a,b)
                append!(pts,get(m.curve_control_points,abs(c),Int[]))
            end
        end
        return pts
    else
        haskey(m.volumes,tag) || throw(ArgumentError("$caller: unknown Volume[$tag]"))
        pts=collect(_model_volume_owned_points(m,tag))
        for shell in m.volumes[tag], signed_surface in m.surface_loops[shell],
            loop in m.surfaces[abs(signed_surface)], signed_curve in m.loops[loop]
            append!(pts,get(m.curve_control_points,abs(signed_curve),Int[]))
        end
        return pts
    end
end

# A staged volume-encoding update, validated against pre-transform state and
# applied only after every point move commits.
function _plan_volume_transform(m::GeoModel, tag::Int, t::_AffineTransform,
                                caller)
    if haskey(m.cylinders,tag)
        rec=m.cylinders[tag]
        tr=_transform_axis_encoding(t,rec.center,rec.axis,caller,"Cylinder[$tag]")
        return (dict=:cylinders,tag=tag,
                record=(center=tr.center,axis=tr.axis,radius=rec.radius*tr.perp,
                        height=tr.height))
    elseif haskey(m.spheres,tag)
        rec=m.spheres[tag]
        s2=_linear_isotropy(t.linear)
        s2===nothing && throw(ArgumentError(
            "$caller: transform is not representable on implicit Sphere[$tag] — " *
            "it requires an isotropic linear part"))
        return (dict=:spheres,tag=tag,
                record=(center=_finite_result(
                            _occ_trsf_apply(t,rec.center),caller),
                        radius=rec.radius*sqrt(s2)))
    elseif haskey(m.cones,tag)
        rec=m.cones[tag]
        tr=_transform_axis_encoding(t,rec.center,rec.axis,caller,"Cone[$tag]")
        return (dict=:cones,tag=tag,
                record=(center=tr.center,axis=tr.axis,r1=rec.r1*tr.perp,
                        r2=rec.r2*tr.perp,height=tr.height))
    elseif haskey(m.booleans,tag)
        rec=m.boolean_operands[tag]
        record=if rec isa Tuple{Mesh,Mesh}
            (_transform_mesh_snapshot(rec[1],t,caller),
             _transform_mesh_snapshot(rec[2],t,caller))
        else
            merge(rec,(meshes=[_transform_mesh_snapshot(snap,t,caller)
                               for snap in rec.meshes],))
        end
        return (dict=:boolean_operands,tag=tag,record=record)
    else
        return (dict=nothing,tag=tag,record=nothing)
    end
end

function _apply_volume_plan!(m::GeoModel, plan)
    plan.dict===nothing && return nothing
    getfield(m,plan.dict)[plan.tag]=plan.record
    return nothing
end

function _transform_mesh_snapshot(mesh::Mesh, t::_AffineTransform, caller)
    coords=Matrix{Float64}(undef,3,size(mesh.coords,2))
    for i in axes(mesh.coords,2)
        coords[:,i].=_finite_result(
            _occ_trsf_apply(
                t,(mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])),
            caller)
    end
    return Mesh(coords; segs=mesh.segs,tris=mesh.tris,tets=mesh.tets,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
end

function _linear_det(L::NTuple{9,Float64})
    return L[1]*(L[5]*L[9]-L[6]*L[8])-L[2]*(L[4]*L[9]-L[6]*L[7])+
           L[3]*(L[4]*L[8]-L[5]*L[7])
end

# Squared isotropic scale of the linear part, or `nothing` when its Gram matrix
# is not a scaled identity.
function _linear_isotropy(L::NTuple{9,Float64})
    cols=((L[1],L[4],L[7]),(L[2],L[5],L[8]),(L[3],L[6],L[9]))
    d=(_dot(cols[1],cols[1]),_dot(cols[2],cols[2]),_dot(cols[3],cols[3]))
    s2=sum(d)/3
    s2>0 || return nothing
    tol=1e-12*max(1.0,s2)
    (abs(d[1]-s2)<=tol && abs(d[2]-s2)<=tol && abs(d[3]-s2)<=tol &&
     abs(_dot(cols[1],cols[2]))<=tol && abs(_dot(cols[1],cols[3]))<=tol &&
     abs(_dot(cols[2],cols[3]))<=tol) || return nothing
    return s2
end

# Transform an axis-encoded primitive: the axis is a free vector (linear part
# only, length = height) and the cross-section radii scale by the transform's
# stretch perpendicular to it — representable only when that stretch is
# isotropic.
function _transform_axis_encoding(t::_AffineTransform, center, axis, caller, what)
    a=_occ_linear_apply(t,axis)
    n=sqrt(_dot(a,a))
    n>0 || throw(ArgumentError("$caller: transform collapses the axis of $what"))
    (u,v)=_perp_basis(a ./ n)
    uL=_occ_linear_apply(t,u); vL=_occ_linear_apply(t,v)
    g=(_dot(uL,uL),_dot(uL,vL),_dot(vL,vL))
    s2=(g[1]+g[3])/2
    s2>0 || throw(ArgumentError("$caller: transform collapses $what"))
    tol=1e-12*max(1.0,s2)
    (abs(g[1]-s2)<=tol && abs(g[3]-s2)<=tol && abs(g[2])<=tol) || throw(ArgumentError(
        "$caller: transform is not representable on $what — it would warp the cross-section"))
    return (center=_finite_result(_occ_trsf_apply(t,center),caller),
            axis=_finite_result(a,caller),perp=sqrt(s2),height=n)
end

function _perp_basis(n::NTuple{3,Float64})
    ax,ay,az=abs.(n)
    seed=ax<=ay ? (ax<=az ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
                  (ay<=az ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    u=(n[2]*seed[3]-n[3]*seed[2],n[3]*seed[1]-n[1]*seed[3],
       n[1]*seed[2]-n[2]*seed[1])
    u=u ./ sqrt(_dot(u,u))
    v=(n[2]*u[3]-n[3]*u[2],n[3]*u[1]-n[1]*u[3],n[1]*u[2]-n[2]*u[1])
    return (u,v)
end

# Materialized `add_box!` volumes keep their `box_extents` encoding only while
# the shell still owns exactly the encoded corners (`_model_volume_bounds`
# validates this). After a point move, re-derive extents when the moved corners
# still form an axis-aligned box; otherwise drop the encoding so queries and
# meshing fall back to the explicit shell rather than reading stale bounds.
function _reconcile_box_encodings!(m::GeoModel, moved::Set{Int})
    for (tag,ext) in collect(m.box_extents)
        owned=sort!(collect(_model_volume_owned_points(m,tag)))
        any(p->p in moved,owned) || continue
        x0,y0,z0,dx,dy,dz=ext
        scale=max(1.0,abs(x0),abs(y0),abs(z0),abs(dx),abs(dy),abs(dz))
        tol=1e-9*scale
        keep=false
        if length(owned)==8
            mapped=[m.points[p] for p in owned]
            lo=ntuple(k->minimum(getindex.(mapped,k)),3)
            hi=ntuple(k->maximum(getindex.(mapped,k)),3)
            ext2=hi .- lo
            corners=[(lo[1]+ix*ext2[1],lo[2]+iy*ext2[2],lo[3]+iz*ext2[3])
                     for ix in (0,1),iy in (0,1),iz in (0,1)]
            keep=all(e->e>0,ext2) && length(unique(mapped))==8 && all(
                c->any(p->_points_close(p,c,tol),mapped),corners)
            keep && (m.box_extents[tag]=(lo[1],lo[2],lo[3],
                                         ext2[1],ext2[2],ext2[3]))
        end
        keep || delete!(m.box_extents,tag)
    end
    return nothing
end

"""
    coherence!(model; tol) -> Bool

Merge coincident entities the way Gmsh's `ReplaceAllDuplicates` does after a
transform: points within `tol` (relative to the model's coordinate scale) fuse
into the lowest tag, then curves sharing an endpoint pair, then surfaces whose
flattened boundary-curve sequences match. References are rewired through
physical groups, embeddings, periodic relations, discrete boundaries, compound
members, and transfinite corner lists. Automatic tag counters reset to the
largest surviving tag per merged dimension, matching `Geometry.OldNewReg`.
"""
function coherence!(m::GeoModel; tol::Real=_COHERENCE_RTOL)
    maps=_coherence_merge!(m,tol)
    return !isempty(maps.points) || !isempty(maps.curves) ||
           !isempty(maps.surfaces)
end

# The coordinate-scaled tolerance Gmsh's `ComparePosition`-style coincidence
# checks share with the merge passes.
function _coherence_eps(m::GeoModel, tol::Real=_COHERENCE_RTOL)
    scale=isempty(m.points) ? 1.0 :
          max(1.0,maximum(p->maximum(abs.(p)),values(m.points)))
    return Float64(tol)*scale
end

# Run the three merge passes and return each dimension's dropped=>survivor
# tag map so callers can resolve post-merge entity tags.
function _coherence_merge!(m::GeoModel, tol::Real=_COHERENCE_RTOL)
    eps=_coherence_eps(m,tol)
    return (points=_merge_points!(m,eps),
            curves=_merge_curves!(m),
            surfaces=_merge_surfaces!(m))
end

"""
    merge_vertices!(model, tags; caller)

Gmsh's `Coherence Point{...}` / `mergeVertices`: snap every listed point to
the first tag's coordinates, then run the global coherence merge. Unknown
tags raise an error before any position is changed.
"""
function merge_vertices!(m::GeoModel, tags;
                         caller::AbstractString="merge_vertices!")
    length(tags)<2 && return nothing
    haskey(m.points,tags[1]) || throw(ArgumentError(
        "$caller: unknown Point[$(tags[1])]"))
    for t in tags[2:end]
        haskey(m.points,t) || throw(ArgumentError("$caller: unknown Point[$t]"))
    end
    target=m.points[tags[1]]
    for t in tags[2:end]
        m.points[t]=target
    end
    coherence!(m)
    return nothing
end

# Returns the dropped=>survivor tag map; empty when nothing merged.
function _merge_points!(m::GeoModel,eps)
    mapping=Dict{Int,Int}()
    isempty(m.points) && return mapping
    grid=Dict{NTuple{3,Int},Vector{Int}}()
    for tag in sort!(collect(keys(m.points)))
        p=m.points[tag]
        cell=ntuple(i->floor(Int,p[i]/eps),3)
        keep=0
        for dx in -1:1, dy in -1:1, dz in -1:1
            for other in get(grid,(cell[1]+dx,cell[2]+dy,cell[3]+dz),Int[])
                if _points_close(p,m.points[other],eps)
                    keep=other; break
                end
            end
            keep!=0 && break
        end
        if keep==0
            push!(get!(grid,cell,Int[]),tag)
        else
            mapping[tag]=keep
        end
    end
    isempty(mapping) && return mapping
    for (drop,keep) in sort!(collect(mapping))
        _rewire_entity_refs!(m,0,drop,keep)
        for (c,(a,b)) in collect(m.curves)
            (a==drop || b==drop) &&
                (m.curves[c]=(a==drop ? keep : a, b==drop ? keep : b))
        end
        for (c,cps) in collect(m.curve_control_points)
            rewired=[p==drop ? keep : p for p in cps]
            rewired==cps || (m.curve_control_points[c]=rewired)
        end
        delete!(m.points,drop); delete!(m.point_size,drop)
        _drop_entity_state!(m,0,drop)
    end
    m.next_tag[1]=isempty(m.points) ? 0 : maximum(keys(m.points))
    return mapping
end

function _merge_curves!(m::GeoModel)
    # Gmsh keeps signed records: `+c` and its reversed record `-c` are
    # distinct entries, and `CompareTwoCurves` compares them directed
    # (beg/end, or the control-point list when present). A copy created in
    # the reversed direction therefore merges into the surviving *reversed*
    # record, so the dropped tag maps to `-keep`. Group by the canonical
    # undirected key, then resolve each drop's sign from its direction.
    seen=Dict{Tuple{Symbol,Int,Vector{Int},Any},Int}()
    dir=Dict{Int,Vector{Int}}()
    mapping=Dict{Int,Int}()
    for tag in sort!(collect(keys(m.curves)))
        cps=get(m.curve_control_points,tag,Int[])
        a,b=m.curves[tag]
        dkey=isempty(cps) ? Int[a,b] : copy(cps)
        # `CompareTwoCurves` distinguishes by record type (Line/Circle/Ellipse
        # — Gmsh's CIRC and CIRC_INV records are equivalent, which Tessella's
        # unsigned types already express) and control-point count first, so
        # both are part of the undirected key; the directed vector is
        # canonicalized against its own reversal. OCC curves fold their stored
        # geometry in too: two coincident closed circles sharing a seam vertex
        # only describe the same edge when center, axis, and radius agree.
        occ=_occ_geometry(m,tag)
        signature=occ===nothing ? nothing :
                  occ.occ===:circle ? (occ.center,occ.n,occ.r) :
                  occ.occ===:line ? (occ.t0,occ.t1) : nothing
        ukey=(_curve_type(m,tag),length(cps),min(dkey,reverse(dkey)),
              signature)
        if haskey(seen,ukey)
            keep=seen[ukey]
            mapping[tag]=dkey==dir[keep] ? keep : -keep
        else
            seen[ukey]=tag
            dir[tag]=dkey
        end
    end
    isempty(mapping) && return mapping
    for (drop,keep) in sort!(collect(mapping))
        _rewire_entity_refs!(m,1,drop,abs(keep))
        for l in keys(m.loops)
            m.loops[l]=[v==drop ? keep : (v==-drop ? -keep : v)
                        for v in m.loops[l]]
        end
        delete!(m.curves,drop)
        # The dropped curve's control points stay behind as unreferenced
        # vertices, matching Gmsh's destroyed-copy behavior.
        delete!(m.curve_control_points,drop)
        _drop_entity_state!(m,1,drop)
    end
    m.next_tag[2]=isempty(m.curves) ? 0 : maximum(keys(m.curves))
    return mapping
end

function _merge_surfaces!(m::GeoModel)
    # `CompareTwoSurfaces` runs `Compare2Lists` over the generatrices with
    # `CompareAbsCurve`: both lists are sorted first and compared by
    # `abs(num)`, so equality is a multiset of absolute curve tags —
    # insensitive to loop order, starting edge, and orientation.
    seen=Dict{Vector{Int},Int}()
    mapping=Dict{Int,Int}()
    for tag in sort!(collect(keys(m.surfaces)))
        flat=Int[]
        for l in m.surfaces[tag]; append!(flat,abs.(m.loops[l])); end
        sort!(flat)
        haskey(seen,flat) ? (mapping[tag]=seen[flat]) : (seen[flat]=tag)
    end
    isempty(mapping) && return mapping
    for (drop,keep) in sort!(collect(mapping))
        _rewire_entity_refs!(m,2,drop,keep)
        for sl in keys(m.surface_loops)
            m.surface_loops[sl]=[v==drop ? keep : (v==-drop ? -keep : v)
                                 for v in m.surface_loops[sl]]
        end
        orphan_loops=m.surfaces[drop]
        delete!(m.surfaces,drop)
        _drop_entity_state!(m,2,drop)
        referenced=Set{Int}()
        for loops in values(m.surfaces); union!(referenced,loops); end
        for l in orphan_loops
            l in referenced || delete!(m.loops,l)
        end
    end
    m.next_tag[3]=isempty(m.surfaces) ? 0 : maximum(keys(m.surfaces))
    return mapping
end

# Replace references to the dropped entity in every `(dim, tag)`-indexed or
# member-list store. Signed topology memberships (curve loops, surface loops)
# are rewired by the dedicated merge passes above.
function _rewire_entity_refs!(m::GeoModel, dim::Int, drop::Int, keep::Int)
    for (key,members) in m.physical
        key[1]==dim && drop in members &&
            (m.physical[key]=unique!([v==drop ? keep : v for v in members]))
    end
    embeds=Dict{Tuple{Int,Int},Vector{NTuple{2,Int}}}()
    for (key,entities) in m.embeds
        new_key=key==(dim,drop) ? (dim,keep) : key
        merged=unique!([e==(dim,drop) ? (dim,keep) : e for e in entities])
        if haskey(embeds,new_key)
            append!(embeds[new_key],merged); unique!(embeds[new_key])
        else
            embeds[new_key]=merged
        end
    end
    m.embeds=embeds
    periodic=Dict{Tuple{Int,Int},ModelPeriodicConstraint}()
    for (key,c) in m.periodic
        c.dim==dim || (periodic[key]=c; continue)
        Int(c.slave_entity)==drop && continue   # slave entity is gone
        master=Int(c.master_entity)==drop ? keep : Int(c.master_entity)
        periodic[key]=ModelPeriodicConstraint(c.dim,c.slave_entity,
                                              Int32(master),c.affine,
                                              c.reversed,c.atol)
    end
    m.periodic=periodic
    discrete=Dict{Tuple{Int,Int},DiscreteEntity}()
    for (key,ent) in m.discrete
        new_key=key==(dim,drop) ? (dim,keep) : key
        ent.boundary=unique!([b==(dim,drop) ? (dim,keep) : b
                              for b in ent.boundary])
        discrete[new_key]=ent
    end
    m.discrete=discrete
    attached=Dict{Tuple{Int,Int},DiscreteEntity}()
    for (key,ent) in m.meshing.attached
        new_key=key==(dim,drop) ? (dim,keep) : key
        ent.boundary=unique!([b==(dim,drop) ? (dim,keep) : b
                              for b in ent.boundary])
        attached[new_key]=ent
    end
    m.meshing.attached=attached
    m.meshing.compounds=Pair{Int,Vector{Int}}[
        cdim==dim ? cdim=>unique!([t==drop ? keep : t for t in ctags]) :
                    cdim=>copy(ctags)
        for (cdim,ctags) in m.meshing.compounds]
    if dim==0
        for (s,a) in m.meshing.transfinite_surfaces
            drop in a.corners && (m.meshing.transfinite_surfaces[s]=
                (arrangement=a.arrangement,
                 corners=[c==drop ? keep : c for c in a.corners]))
        end
        for v in keys(m.meshing.transfinite_volumes)
            corners=m.meshing.transfinite_volumes[v]
            drop in corners && (m.meshing.transfinite_volumes[v]=
                [c==drop ? keep : c for c in corners])
        end
    end
    for i in eachindex(m.meshing.homology_requests)
        r=m.meshing.homology_requests[i]
        (drop in r.domain || drop in r.subdomain || drop in r.dims) &&
            (m.meshing.homology_requests[i]=
                (kind=r.kind,
                 domain=[t==drop ? keep : t for t in r.domain],
                 subdomain=[t==drop ? keep : t for t in r.subdomain],
                 dims=[t==drop ? keep : t for t in r.dims]))
    end
    return nothing
end

function _drop_entity_state!(m::GeoModel, dim::Int, tag::Int)
    key=(dim,tag)
    delete!(m.entity_names,key)
    delete!(m.entity_visibility,key)
    delete!(m.entity_colors,key)
    delete!(m.meshing.recombine,key)
    delete!(m.meshing.extrude,key)
    delete!(m.meshing.smoothing,key)
    delete!(m.meshing.reverse,key)
    delete!(m.meshing.algorithm,key)
    delete!(m.meshing.size_at_params,key)
    delete!(m.meshing.size_from_boundary,key)
    delete!(m.meshing.attached,key)
    if dim==1
        delete!(m.meshing.transfinite_curves,tag)
        delete!(m.curve_control_points,tag)
        delete!(m.curve_types,tag)
        delete!(m.curve_geometry,tag)
    elseif dim==2
        delete!(m.meshing.transfinite_surfaces,tag)
        delete!(m.surface_types,tag)
        delete!(m.surface_geometry,tag)
    elseif dim==3
        delete!(m.meshing.transfinite_volumes,tag)
        delete!(m.meshing.outward_orientation,tag)
    end
    return nothing
end

# ── Duplicata ────────────────────────────────────────────────────────────────

# Gmsh's shared `Geometry.OldNewReg` tag source: the maximum of the
# curve/surface/volume counters plus every live curve loop, surface loop,
# physical, and non-point discrete tag.
function _geo_newreg_tag(m::GeoModel)
    mx=max(m.next_tag[2],m.next_tag[3],m.next_tag[4],m.physical_tag_max)
    for t in keys(m.loops); mx=max(mx,t); end
    for (d,t) in keys(m.discrete); d>0 && (mx=max(mx,t)); end
    for t in keys(m.surface_loops); mx=max(mx,t); end
    return mx+1
end

_geo_newreg_alloc!(m::GeoModel,dim::Int,caller) =
    _alloc_tag!(m,dim,_geo_newreg_tag(m),caller)

"""
    duplicate_entities!(model, entities; caller) -> copies

Deep-copy each listed native entity at its current position, mirroring Gmsh's
`Duplicata`: copies carry no physical memberships or names, curve copies
allocate their two coincident endpoint vertices plus Gmsh's transient beg/end
vertex pair, and surface/volume copies recurse through boundary topology with
tags drawn from the `OldNewReg` counter. Coincident results stay unmerged —
an enclosing transform's coherence pass collapses them exactly like Gmsh.
"""
function duplicate_entities!(m::GeoModel,
                             entities::AbstractVector{<:Tuple{Integer,Integer}};
                             caller::AbstractString="duplicate_entities!")
    out=NTuple{2,Int}[]
    # One copy map per call, consulted only by OCC-geometry entities: like
    # `OCC_Internals.copy`, a duplicated edge/vertex is shared by every copied
    # face that references it, keeping periodic wires consistent.
    memo=(Dict{Int,Int}(),Dict{Int,Int}(),Dict{Int,Int}())
    for (dim,tag) in entities
        d=_dimension(dim,caller); tg=_tag(tag,caller,d)
        haskey(m.discrete,(d,tg)) && throw(ArgumentError(
            "$caller: discrete $(_entity_label(d))[$tg] cannot be duplicated"))
        push!(out,(d,_duplicate_entity!(m,d,tg,caller,memo)))
    end
    return out
end

function _duplicate_entity!(m::GeoModel, d::Int, tag::Int, caller, memo)
    if d==0
        haskey(m.points,tag) || throw(ArgumentError("$caller: unknown Point[$tag]"))
        return _fresh_point_copy!(m,tag,caller)
    elseif d==1
        return _duplicate_curve!(m,tag,caller;memo=memo)
    elseif d==2
        return _duplicate_surface!(m,tag,caller,memo)
    else
        return _duplicate_volume!(m,tag,caller,memo)
    end
end

function _fresh_point_copy!(m::GeoModel, src::Int, caller, pmemo=nothing)
    pmemo!==nothing && haskey(pmemo,src) && return pmemo[src]
    t=_alloc_tag!(m,0,0,caller)
    m.points[t]=m.points[src]
    haskey(m.point_size,src) && (m.point_size[t]=m.point_size[src])
    pmemo!==nothing && (pmemo[src]=t)
    return t
end

# True when the surface's own geometry or any boundary curve is an OCC record —
# such copies must share edge/vertex copies across the whole operation.
function _occ_surface_entity(m::GeoModel, src::Int)
    haskey(m.surface_geometry,src) && return true
    for l in m.surfaces[src], c in m.loops[l]
        _occ_geometry(m,abs(c))!==nothing && return true
    end
    return false
end

function _occ_volume_entity(m::GeoModel, src::Int)
    for sl in m.volumes[src], s in m.surface_loops[sl]
        _occ_surface_entity(m,abs(s)) && return true
    end
    return false
end

# Gmsh's `DuplicateCurve` allocates the curve tag first (shared `NEWREG`
# counter), then one coincident copy per control point — stored on the copy as
# its `curve_control_points` — and finally the wired beg/end vertex copies.
# For a plain line the source control list is its endpoint pair, so each
# curve copy burns four point tags: two orphans, then the wired pair.
# `reversed=true` mirrors `DuplicateCurve` applied to Gmsh's reversed-curve
# record: control points copy in reversed order and the copy is wired
# end-to-begin (beg copy = copy of the source's end vertex).
function _duplicate_curve!(m::GeoModel, src::Int, caller; reversed::Bool=false,
                           memo=nothing)
    haskey(m.curves,src) || throw(ArgumentError("$caller: unknown Curve[$src]"))
    a,b=m.curves[src]
    occ=_occ_geometry(m,src)
    # OCC copies share through `memo`; a reversed copy is a distinct edge.
    !reversed && occ!==nothing && memo!==nothing && haskey(memo[2],src) &&
        return memo[2][src]
    # OCC curves carry their parametrization in `curve_geometry`, not control
    # points — duplicating must not invent any.
    source_cps=occ===nothing ? get(m.curve_control_points,src,Int[a,b]) : Int[]
    geometry=haskey(m.curve_geometry,src) ? m.curve_geometry[src] : nothing
    if reversed
        (a,b)=(b,a)
        # Gmsh's `CreateReversedCurve` inverts the control list except for
        # ellipses, where the center and major-axis points keep their
        # positions: [end, center, major, start].
        source_cps=_curve_type(m,src)==:ellipse && !isempty(source_cps) ?
            Int[source_cps[4],source_cps[2],source_cps[3],source_cps[1]] :
            reverse(source_cps)
        # A reversed OCC circle traverses its forward image at -t, so the
        # second in-plane direction flips sign: p_rev(t) = p_fwd(-t) =
        # center + r(cos t·X - sin t·Y).
        occ!==nothing && occ.occ===:circle &&
            (geometry=(occ=:circle,center=occ.center,
                       n=(-occ.n[1],-occ.n[2],-occ.n[3]),X=occ.X,
                       Y=(-occ.Y[1],-occ.Y[2],-occ.Y[3]),r=occ.r,
                       t0=-occ.t1,t1=-occ.t0))
        # `CreateReversedCurve` mirrors the Nurbs knot vector in place
        # (k'[deg+N-i] = k[i]) and stores `ubeg' = 1-uend`, `uend' = 1-ubeg`.
        geometry!==nothing && hasproperty(geometry,:knots) &&
            (geometry=(knots=reverse(geometry.knots),deg=geometry.deg,
                       ubeg=1.0-geometry.uend,uend=1.0-geometry.ubeg))
    end
    t=_geo_newreg_alloc!(m,1,caller)
    isempty(source_cps) || (m.curve_control_points[t]=
        [_fresh_point_copy!(m,c,caller) for c in source_cps])
    pmemo=occ===nothing || memo===nothing ? nothing : memo[1]
    pa=_fresh_point_copy!(m,a,caller,pmemo)
    # A closed OCC edge's endpoints are one shared vertex, not two copies.
    pb=a==b ? pa : _fresh_point_copy!(m,b,caller,pmemo)
    m.curves[t]=(pa,pb)
    haskey(m.curve_types,src) && (m.curve_types[t]=m.curve_types[src])
    geometry===nothing || (m.curve_geometry[t]=geometry)
    !reversed && occ!==nothing && memo!==nothing && (memo[2][src]=t)
    return t
end

function _duplicate_surface!(m::GeoModel, src::Int, caller, memo=nothing)
    haskey(m.surfaces,src) || throw(ArgumentError("$caller: unknown Surface[$src]"))
    occ=_occ_surface_entity(m,src)
    occ && memo!==nothing && haskey(memo[3],src) && return memo[3][src]
    t=_geo_newreg_alloc!(m,2,caller)
    loops=Int[]
    for (i,l) in enumerate(m.surfaces[src])
        # OCC copies share each memoized edge and keep the orientation sign in
        # the loop (`BRepBuilderAPI_Copy` preserves the wire's orientation
        # flags). The built-in `DuplicateSurface` instead resolves each
        # generatrix record first: a `-c` member duplicates the reversed
        # record (reversed data — mirrored Nurbs knots, `[1-uend,1-ubeg]`)
        # and the copy enters the new loop under a positive tag.
        curves=occ ?
            [sign(c)*_duplicate_curve!(m,abs(c),caller;memo=memo)
             for c in m.loops[l]] :
            [_duplicate_curve!(m,abs(c),caller;reversed=c<0)
             for c in m.loops[l]]
        lt=(i==1 && !haskey(m.loops,t)) ? t :
           ((isempty(m.loops) ? 0 : maximum(keys(m.loops)))+1)
        m.loops[lt]=curves
        push!(loops,lt)
    end
    m.surfaces[t]=loops
    haskey(m.surface_types,src) && (m.surface_types[t]=m.surface_types[src])
    haskey(m.surface_geometry,src) &&
        (m.surface_geometry[t]=m.surface_geometry[src])
    occ && memo!==nothing && (memo[3][src]=t)
    return t
end

function _duplicate_volume!(m::GeoModel, src::Int, caller, memo=nothing)
    haskey(m.volumes,src) || throw(ArgumentError("$caller: unknown Volume[$src]"))
    occ=_occ_volume_entity(m,src)
    t=_geo_newreg_alloc!(m,3,caller)
    if !isempty(m.volumes[src])
        shells=Int[]
        for (i,sl) in enumerate(m.volumes[src])
            surfs=[sign(s)*_duplicate_surface!(m,abs(s),caller,
                                               occ ? memo : nothing)
                   for s in m.surface_loops[sl]]
            slt=(i==1 && !haskey(m.surface_loops,t)) ? t :
                ((isempty(m.surface_loops) ? 0 : maximum(keys(m.surface_loops)))+1)
            m.surface_loops[slt]=surfs
            push!(shells,slt)
        end
        m.volumes[t]=shells
    else
        m.volumes[t]=Int[]
    end
    haskey(m.box_extents,src) && (m.box_extents[t]=m.box_extents[src])
    haskey(m.cylinders,src) && (m.cylinders[t]=m.cylinders[src])
    haskey(m.spheres,src) && (m.spheres[t]=m.spheres[src])
    haskey(m.cones,src) && (m.cones[t]=m.cones[src])
    if haskey(m.booleans,src)
        m.booleans[t]=m.booleans[src]
        m.boolean_operands[t]=m.boolean_operands[src]
        # the copy owns its own shells — a component-linked source copies as a
        # self-linked component volume, not a second reference to src's link
        haskey(m.boolean_components,src) && (m.boolean_components[t]=t)
    end
    return t
end

# ── Extrude ────────────────────────────────────────────────────────────────
#
# Gmsh built-in translational `Extrude {dx,dy,dz} { ... }` (`ExtrudeShapes`
# plus `ExtrudePoint`/`ExtrudeCurve`/`ExtrudeSurface` in Geo.cpp). Every
# listed entity is deep-copied and translated; the copy is the "top" and the
# swept topology is the "body": a point yields a connecting curve, a curve a
# ruled lateral surface, and a surface a volume bounded by the source
# (orientation -1), the top (+1), and one lateral surface per generatrix
# (signed by the generatrix's sign). The flat result list holds `[top, body]`
# per input entity.
#
# Tag allocation follows Gmsh's counter order exactly. A curve copy burns one
# shared `NEWREG` tag plus four point tags (two control-point orphans, then
# the wired pair); a generatrix extrusion burns its own curve copy, the two
# endpoint copies, the connecting curves, and the lateral surface; a surface
# extrusion additionally burns one transient surface tag (Gmsh's
# backward-compatible re-tag of the top surface) and allocates the volume
# from the dedicated volume counter (the `oldNewreg=0` hack). Coincident
# copies collapse in a single `ReplaceAllDuplicates` pass at the end of each
# top-level entity's extrusion; a fully coincident point or curve extrusion
# returns the source tag and leaves the unmerged copies behind, exactly like
# Gmsh.

# Gmsh's `ExtrudeParams::geo`: the kernel `type` plus its operands. `:translate`
# carries the displacement `T`; `:rotate` carries the revolved `axis`/`origin`/
# `angle` and `rot`, the precomputed Gmsh rotation transform (three sequential
# 4×4 applications, bit-identical to `ApplyTransformationToPoint`).
_extrude_spec_translate(delta::NTuple{3,Float64}) = (type=:translate,T=delta)
function _extrude_spec_rotate(axis, origin, angle, caller)
    a=_finite_vector3(axis,caller,"rotation axis")
    o=_finite_vector3(origin,caller,"rotation origin")
    θ=_finite_scalar(angle,caller,"rotation angle")
    return (type=:rotate,rot=_affine_rotation(a,o,θ,caller),
            axis=a,origin=o,angle=θ)
end

# The copy transform applied to chapeau points: the displacement for
# `:translate`, the translate/rotate/translate matrix sequence for `:rotate`.
_extrude_move(spec,p) = spec.type===:rotate ?
    _affine_apply_steps(spec.rot,p) : p .+ spec.T

"""
    extrude_entities!(model, entities, delta; params=nothing, caller) -> tags

Translate-extrude every listed `(dim, tag)` entity — Points, Curves, and
Surfaces only — returning the flat `[top, body]` tag list Gmsh's built-in
`Extrude` produces. `params` is a `_GeoExtrudeParams` record (or `nothing`)
attached to every created entity. All inputs are validated before any entity
is created, so unknown or unsupported entities leave the model unchanged.
"""
function extrude_entities!(m::GeoModel,
                           entities::AbstractVector{<:Tuple{Integer,Integer}},
                           delta::NTuple{3,Float64};
                           params::Union{Nothing,_GeoExtrudeParams}=nothing,
                           return_lateral::Bool=true,
                           caller::AbstractString="extrude_entities!")
    all(isfinite,delta) || throw(ArgumentError(
        "$caller: extrusion delta must be finite"))
    return _extrude_entities!(m,entities,_extrude_spec_translate(delta);
                             params=params,return_lateral=return_lateral,
                             caller=caller)
end

"""
    revolve_entities!(model, entities, axis, origin, angle; ...) -> tags

Rotate-extrude every listed `(dim, tag)` entity by `angle` radians about the
axis through `origin` with direction `axis` — Gmsh's `revolve` (`Extrude
{{axis}, {point}, angle}`) with built-in-kernel semantics: swept vertices
produce `Circle` arcs wired `[start, axis-center, end]`, curve extrusions
produce ruled (or triangular) lateral surfaces, and surface extrusions
produce volumes. Tag allocation, signed-generatrix rules, `out` list layout,
and merge semantics match `ExtrudeShapes(ROTATE, ...)`.
"""
function revolve_entities!(m::GeoModel,
                           entities::AbstractVector{<:Tuple{Integer,Integer}},
                           axis, origin, angle;
                           params::Union{Nothing,_GeoExtrudeParams}=nothing,
                           return_lateral::Bool=true,
                           caller::AbstractString="revolve_entities!")
    spec=_extrude_spec_rotate(axis,origin,angle,caller)
    return _extrude_entities!(m,entities,spec;
                             params=params,return_lateral=return_lateral,
                             caller=caller)
end

function _extrude_entities!(m::GeoModel,
                            entities::AbstractVector{<:Tuple{Integer,Integer}},
                            spec;
                            params::Union{Nothing,_GeoExtrudeParams}=nothing,
                            return_lateral::Bool=true,
                            caller::AbstractString="extrude_entities!")
    normalized=NTuple{2,Int}[]
    for (dim,tag) in entities
        d=_dimension(dim,caller)
        t=Int(tag)
        t==0 && throw(ArgumentError(
            "$caller: $(_entity_label(d)) tags must be nonzero"))
        tg=_tag(abs(t),caller,d)
        haskey(m.discrete,(d,tg)) && throw(ArgumentError(
            "$caller: discrete $(_entity_label(d))[$tg] cannot be extruded"))
        d==3 && throw(ArgumentError(
            "$caller: impossible to extrude Volume[$tg]"))
        store=d==0 ? m.points : d==1 ? m.curves : m.surfaces
        haskey(store,tg) || throw(ArgumentError(
            "$caller: unknown $(_entity_label(d))[$tg]"))
        # Gmsh keys the signed curve/surface records on the absolute tag; the
        # sign selects the reversed record for curves and feeds extrusion
        # metadata for surfaces. Point signs resolve to the absolute tag.
        push!(normalized,(d,d==0 ? tg : t))
    end
    out=Int[]
    for (d,tg) in normalized
        if d==0
            _extrude_point!(m,tg,spec,out,params,caller)
        elseif d==1
            _extrude_curve!(m,tg,spec,out,params,return_lateral,caller)
        else
            _extrude_surface!(m,tg,spec,out,params,return_lateral,caller)
        end
    end
    return out
end

# `ExtrudePoint`: copy the vertex, transform the copy, and wire a connecting
# curve — unless the copy still coincides with the source, in which case the
# source tag is returned and the orphan copy stays unmerged (the coherence
# pass is skipped entirely). `final=false` suppresses the merge for the
# nested endpoint extrusions inside `ExtrudeCurve`.
function _extrude_point_copy!(m::GeoModel, src::Int, spec,
                              params, caller)
    chapeau=_fresh_point_copy!(m,src,caller)
    m.points[chapeau]=_finite_result(
        _extrude_move(spec,m.points[chapeau]),caller)
    _points_close(m.points[chapeau],m.points[src],_coherence_eps(m)) &&
        return (nothing,chapeau)
    curve=_geo_newreg_alloc!(m,1,caller)
    m.curves[curve]=(src,chapeau)
    if spec.type===:rotate
        # `MSH_SEGM_CIRC`: the swept arc is a three-point circle
        # [start, axis-center, end]; the center is the source's orthogonal
        # projection onto the axis, allocated *after* the curve tag exactly
        # like `DuplicateVertex` inside `ExtrudePoint`. `CreateCurve` seeds
        # `Circle.n=(0,0,1)` — the `EndCurve` fallback normal for degenerate
        # (half-turn) arcs.
        m.curve_types[curve]=:circle
        ax=_arc_norme(spec.axis)
        o=spec.origin;p=m.points[src]
        d=_arc_dot((p[1]-o[1],p[2]-o[2],p[3]-o[3]),ax)
        center=_alloc_tag!(m,0,0,caller)
        m.points[center]=(fma(d,ax[1],o[1]),fma(d,ax[2],o[2]),fma(d,ax[3],o[3]))
        m.curve_control_points[curve]=Int[src,center,chapeau]
        m.curve_geometry[curve]=(n=(0.0,0.0,1.0),)
    end
    params!==nothing && (m.meshing.extrude[(1,curve)]=params)
    return (curve,chapeau)
end

function _extrude_point!(m::GeoModel, src::Int, spec,
                         out::Vector{Int}, params, caller)
    (curve,chapeau)=_extrude_point_copy!(m,src,spec,params,caller)
    if curve===nothing
        push!(out,src)
        return nothing
    end
    maps=_coherence_merge!(m)
    push!(out,get(maps.points,chapeau,chapeau))
    haskey(m.curves,curve) && push!(out,curve)
    return nothing
end

# `ExtrudeCurve`: duplicate the curve, transform the copy's point set, extrude
# the (possibly signed) endpoints into connecting curves, and wire the lateral
# surface. Returns `(surf, chapeau)` or `nothing` when both connecting curves
# collapsed — Gmsh's `if(!CurveBeg && !CurveEnd) return ic` path.
function _extrude_curve_lateral!(m::GeoModel, c::Int, spec,
                                 params, caller)
    src=abs(c)
    a,b=m.curves[src]
    # For a negative generatrix Gmsh extrudes the reversed-curve record:
    # begin/end swap and the copy is wired end-to-begin.
    sbeg,send=c>0 ? (a,b) : (b,a)
    chapeau=_duplicate_curve!(m,src,caller;reversed=c<0)
    params!==nothing && (m.meshing.extrude[(1,chapeau)]=params)
    _transform_curve_points!(m,chapeau,spec,caller)
    _extrude_occ_transform!(m,1,chapeau,spec,caller)
    (cbeg,_)=_extrude_point_copy!(m,sbeg,spec,params,caller)
    (cend,_)=_extrude_point_copy!(m,send,spec,params,caller)
    (cbeg===nothing && cend===nothing) && return nothing
    return (_extrude_lateral_surface!(m,c,chapeau,cbeg,cend,params,caller),
            chapeau)
end

function _extrude_curve!(m::GeoModel, c::Int, spec,
                         out::Vector{Int}, params, return_lateral::Bool,
                         caller)
    result=_extrude_curve_lateral!(m,c,spec,params,caller)
    if result===nothing
        # `if(!CurveBeg && !CurveEnd) return ic` — Gmsh returns the signed
        # input tag and leaves the unmerged copy behind.
        push!(out,c)
        return nothing
    end
    (surf,chapeau)=result
    maps=_coherence_merge!(m)
    top=get(maps.curves,chapeau,chapeau)
    push!(out,top)
    if haskey(m.surfaces,surf)
        push!(out,surf)
        # `Geometry.ExtrudeReturnLateralEntities` (default on) appends the
        # lateral generatrices besides the source and top curves. Gmsh's
        # filter compares `abs(record tag)` against the *signed* input, so a
        # negative generatrix keeps its own `-tag` in the result list.
        return_lateral || return nothing
        for g in m.loops[only(m.surfaces[surf])]
            (abs(g)==c || abs(g)==abs(top)) || push!(out,g)
        end
    end
    return nothing
end

# The lateral surface of a curve extrusion — generatrices `pc, CurveEnd,
# -chapeau, -CurveBeg`, or Gmsh's triangular three-edge variants when one
# connecting curve collapsed (`MSH_SURF_TRIC`).
function _extrude_lateral_surface!(m::GeoModel, src::Int, chapeau::Int,
                                   cbeg, cend, params, caller)
    surf=_geo_newreg_alloc!(m,2,caller)
    # `src` is the signed generatrix record: `pc` for `ic>0`, the reversed
    # record for `ic<0` — the loop carries the source with its sign.
    generatrices=if cbeg===nothing
        Int[src,cend,-chapeau]
    elseif cend===nothing
        Int[-chapeau,-cbeg,src]
    else
        Int[src,cend,-chapeau,-cbeg]
    end
    lt=_geo_derived_loop_tag(m,surf)
    m.loops[lt]=generatrices
    m.surfaces[surf]=Int[lt]
    # `ExtrudeCurve`: a collapsed side gives a three-edge `MSH_SURF_TRIC`
    # patch, otherwise the four-generatrix `MSH_SURF_REGL` ruled surface.
    m.surface_types[surf]=length(generatrices)==4 ? :ruled : :tric
    params!==nothing && (m.meshing.extrude[(2,surf)]=params)
    return surf
end

# `ExtrudeSurface`: burn the transient duplicate-surface tag, duplicate every
# generatrix into the top copy, allocate the volume on the dedicated volume
# counter, sweep each generatrix into a lateral surface, transform the top,
# re-tag it through a fresh `NEWSURFACE`, and wire the shell.
function _extrude_surface!(m::GeoModel, is::Int, spec,
                           out::Vector{Int}, params, return_lateral::Bool,
                           caller)
    tag=abs(is)   # `ps = FindSurface(std::abs(is))`; the sign is metadata-only
    _geo_newreg_alloc!(m,2,caller)   # transient chapeau tag — burned
    top_loops=[Int[_duplicate_curve!(m,abs(c),caller;reversed=c<0)
                   for c in m.loops[l]] for l in m.surfaces[tag]]
    if params!==nothing
        for copies in top_loops, c in copies
            m.meshing.extrude[(1,c)]=params
        end
    end
    vol=_alloc_tag!(m,3,0,caller)
    params!==nothing && (m.meshing.extrude[(3,vol)]=params)
    laterals=Int[]
    for l in m.surfaces[tag], c in m.loops[l]
        result=_extrude_curve_lateral!(m,c,spec,params,caller)
        result===nothing && continue
        push!(laterals,c<0 ? -result[1] : result[1])
    end
    for copies in top_loops, c in copies
        _transform_curve_points!(m,c,spec,caller)
        _extrude_occ_transform!(m,1,c,spec,caller)
    end
    top=_geo_newreg_alloc!(m,2,caller)
    loop_tags=Int[]
    for (i,copies) in enumerate(top_loops)
        lt=i==1 ? _geo_derived_loop_tag(m,top) : _geo_next_loop_tag(m)
        m.loops[lt]=copies
        push!(loop_tags,lt)
    end
    m.surfaces[top]=loop_tags
    # The top is `DuplicateSurface` output: it inherits the source's record
    # type (`Plane Surface` extrudes to a `Plane` cap) and geometry metadata.
    # OCC geometry records move with the copy — an unmodified record would
    # answer queries at the source's coordinates.
    haskey(m.surface_types,tag) && (m.surface_types[top]=m.surface_types[tag])
    haskey(m.surface_geometry,tag) &&
        (m.surface_geometry[top]=_extrude_occ_surface_transform(
            m.surface_geometry[tag],spec,caller,"Surface[$top]"))
    params!==nothing && (m.meshing.extrude[(2,top)]=params)
    slt=haskey(m.surface_loops,vol) ? _geo_next_surface_loop_tag(m) : vol
    m.surface_loops[slt]=vcat(-tag,top,laterals)
    m.volumes[vol]=Int[slt]
    maps=_coherence_merge!(m)
    top_out=get(maps.surfaces,top,top)
    push!(out,top_out)
    push!(out,vol)
    # `Geometry.ExtrudeReturnLateralEntities` (default on) appends the
    # volume's boundary surfaces besides the source and top — the same
    # signed comparison Gmsh applies, so `is<0` keeps the `+tag` source
    # record in the result list. Boundary entries carry orientation signs;
    # the appended record tags are positive.
    return_lateral || return nothing
    for s in m.surface_loops[slt]
        side=get(maps.surfaces,abs(s),abs(s))
        (side==is || side==top_out) || push!(out,side)
    end
    return nothing
end

# Derived curve loops are keyed by the owning surface's tag when it is free,
# else the next unused loop tag — matching the `Curve Loop(N)` numbering in
# Gmsh's unrolled `.geo` output.
_geo_derived_loop_tag(m::GeoModel, surface_tag::Int) =
    haskey(m.loops,surface_tag) ? _geo_next_loop_tag(m) : surface_tag
_geo_next_loop_tag(m::GeoModel) =
    (isempty(m.loops) ? 0 : maximum(keys(m.loops)))+1
_geo_next_surface_loop_tag(m::GeoModel) =
    (isempty(m.surface_loops) ? 0 : maximum(keys(m.surface_loops)))+1

# Transform every vertex a curve owns — the wired endpoints plus its
# `curve_control_points` copies (Gmsh's `ApplyTransformationToCurve` reaches
# both).
function _transform_curve_points!(m::GeoModel, curve::Int, spec, caller)
    a,b=m.curves[curve]
    pts=Int[a,b]
    append!(pts,get(m.curve_control_points,curve,Int[]))
    for p in unique!(pts)
        m.points[p]=_finite_result(_extrude_move(spec,m.points[p]),caller)
    end
    return nothing
end

# A duplicated OCC curve record moves with its copy: a translate/rotate
# extrusion is always rigid, so the transform stays representable. Built-in
# arc records carry only a fallback normal and need no update.
function _extrude_occ_transform!(m::GeoModel, dim::Int, tag::Int, spec, caller)
    if dim==1
        g=_occ_geometry(m,tag)
        g===nothing && return nothing
        m.curve_geometry[tag]=_transform_occ_curve_record(
            m,tag,g,_extrude_affine(spec),caller,"Curve[$tag]")
    else
        g=get(m.surface_geometry,tag,nothing)
        (g===nothing || !hasproperty(g,:occ)) && return nothing
        m.surface_geometry[tag]=_transform_occ_surface_record(
            g,_extrude_affine(spec),caller,"Surface[$tag]")
    end
    return nothing
end

# The `_AffineTransform` equivalent of the extrusion's copy transform —
# only OCC geometry records read it (chapeau points use `_extrude_move`).
_extrude_affine(spec) = spec.type===:rotate ? spec.rot :
    _affine_translation(spec.T)

_extrude_occ_surface_transform(g, spec, caller, what) =
    (g===nothing || !hasproperty(g,:occ)) ? g :
    _transform_occ_surface_record(g,_extrude_affine(spec),caller,what)
