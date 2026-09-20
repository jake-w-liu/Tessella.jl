# ── Spline-family curves ─────────────────────────────────────────────────────
#
# Gmsh built-in-kernel `Spline`, `BSpline`, `Bezier`, and `Nurbs` records
# (`MSH_SEGM_SPLN`/`BSPLN`/`BEZIER`/`NURBS`; all report geomType "Nurb").
# Ordered control-point tags live in `curve_control_points`; `Nurbs` stores its
# float32-rounded knot vector, inferred degree, and raw parameter interval in
# `curve_geometry`, mirroring `Curve::k`/`degre`/`ubeg`/`uend`. Evaluation,
# derivatives, and projection port `GeoInterpolation.cpp`'s `InterpolateCurve`
# and `GEdge`'s refineProjection/XYZToU/bounds scan.
#
# A curve is closed (Gmsh's `beg == end` periodic branch) exactly when its first
# and last control-point tags coincide. Curved curves remain excluded from the
# straight-line meshing paths via `_model_require_line_curve`.

const _SPLINE_CURVE_TYPES = (:spline, :bspline, :bezier, :nurbs)

# `CreateCurve(tag, Typ, Order, Liste, Knots, -1, -1, 0., 1.)` — beg/end are the
# first/last control-point vertices, and unknown tags are construction errors.
function _add_spline_family!(m::GeoModel, points, kind::Symbol,
                             tag::Integer, caller::AbstractString,
                             what::AbstractString;
                             literal_zero::Bool=false)
    pts=Int[_tag(p,caller,1) for p in points]
    length(pts)>=2 || throw(ArgumentError(
        "$caller: $what curve requires at least 2 control points"))
    for p in pts
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: unknown control point $p in Curve $tag"))
    end
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller;literal_zero=literal_zero)
    (haskey(m.curves,t) || haskey(m.discrete,(1,t))) && throw(ArgumentError(
        "$caller: Curve[$t] already exists"))
    m.curves[t]=(first(pts),last(pts))
    m.curve_control_points[t]=pts
    m.curve_types[t]=kind
    return t
end

"""
    add_spline!(model, points; tag=0) -> tag

Add a Catmull-Rom `Spline` through the ordered control Points, mirroring
`Spline(tag) = {p1,...,pn}` (`MSH_SEGM_SPLN`, uniform knots, extrapolated
ghost ends, cyclic when the first and last tags coincide).
"""
add_spline!(m::GeoModel, points; tag::Integer=0, _zero_literal::Bool=false) =
    _add_spline_family!(m,points,:spline,tag,"add_spline!","Spline";
                       literal_zero=_zero_literal)

"""
    add_bspline!(model, points; tag=0) -> tag

Add a uniform clamped `BSpline` over the ordered control Points, mirroring
`BSpline(tag) = {p1,...,pn}` (`MSH_SEGM_BSPLN`, Gmsh's `InterpolateUBS`
piecewise-cubic construction).
"""
add_bspline!(m::GeoModel, points; tag::Integer=0, _zero_literal::Bool=false) =
    _add_spline_family!(m,points,:bspline,tag,"add_bspline!","BSpline";
                       literal_zero=_zero_literal)

"""
    add_bezier!(model, points; tag=0) -> tag

Add a `Bezier` curve over the ordered control Points, mirroring
`Bezier(tag) = {p1,...,pn}` (`MSH_SEGM_BEZIER`, De Casteljau interpolation).
"""
add_bezier!(m::GeoModel, points; tag::Integer=0, _zero_literal::Bool=false) =
    _add_spline_family!(m,points,:bezier,tag,"add_bezier!","Bezier";
                       literal_zero=_zero_literal)

"""
    add_nurbs!(model, points, knots; tag=0) -> tag

Add a `Nurbs` curve mirroring `Nurbs(tag) = {p1,...,pn} Knots {k0,...} Order o`
(`MSH_SEGM_NURBS`). The built-in kernel ignores `Order` and infers the degree
as `length(knots) - length(points) - 1`; an empty `knots` builds a `BSpline`
record exactly like `GEO_Internals::addBSpline`. Knots must be finite and
nondecreasing, and the inferred degree must lie in
`0 <= degree <= length(points)-1` — outside that range Gmsh's `findSpan`/
`basisFuns` read the control-point list out of bounds, so Tessella rejects the
record instead (a degree-zero Nurbs is the piecewise-constant curve
`findSpan`/`basisFuns` produce). The parameter interval is
`[knots[1],knots[end]]` (`Curve::ubeg`/`uend`) and `getValue` evaluates raw
knot coordinates.
"""
function add_nurbs!(m::GeoModel, points, knots; tag::Integer=0,
                    _zero_literal::Bool=false)
    caller="add_nurbs!"
    knots isa Bool && throw(ArgumentError(
        "$caller: knots must not be Bool"))
    knot_vector=Float64[]
    for (index,k) in enumerate(knots)
        k isa Bool && throw(ArgumentError(
            "$caller: knot $(index-1) must not be Bool"))
        value=try
            Float64(k)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: knot $(index-1) must be real: $(sprint(showerror,err))"))
        end
        isfinite(value) || throw(ArgumentError(
            "$caller: knot $(index-1) is not finite"))
        push!(knot_vector,value)
    end
    # `GEO_Internals::addBSpline` — an empty `Knots {}` builds the plain
    # `MSH_SEGM_BSPLN` record, not a Nurbs.
    isempty(knot_vector) && return _add_spline_family!(
        m,points,:bspline,tag,caller,"BSpline";literal_zero=_zero_literal)
    pts=Int[_tag(p,caller,1) for p in points]
    length(pts)>=2 || throw(ArgumentError(
        "$caller: Nurbs curve requires at least 2 control points"))
    for p in pts
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: unknown control point $p in Curve $tag"))
    end
    degree=length(knot_vector)-length(pts)-1
    degree>=0 || throw(ArgumentError(
        "$caller: Nurbs needs at least length(points)+1 knots " *
        "(degree $degree is not a curve)"))
    # `findSpan` returns spans in `[deg, N-1]` and `basisFuns` plus the
    # control-point accumulation index `Control_Points[span-deg+i]`; an
    # inferred degree above `N-1` (more than 2N knots) sends that index
    # negative — upstream `List_Read` silently reads out of bounds — so the
    # malformed record is rejected at construction like a negative degree.
    degree<=length(pts)-1 || throw(ArgumentError(
        "$caller: Nurbs with $(length(pts)) control points needs at most " *
        "$(2*length(pts)) knots (inferred degree $degree exceeds the polygon)"))
    for index in 2:length(knot_vector)
        knot_vector[index]>=knot_vector[index-1] || throw(ArgumentError(
            "$caller: Nurbs knots must be nondecreasing"))
    end
    t=_alloc_tag!(m,1,_tag(tag,caller,1),caller;literal_zero=_zero_literal)
    (haskey(m.curves,t) || haskey(m.discrete,(1,t))) && throw(ArgumentError(
        "$caller: Curve[$t] already exists"))
    m.curves[t]=(first(pts),last(pts))
    m.curve_control_points[t]=pts
    m.curve_types[t]=:nurbs
    # `CreateCurve` keeps `ubeg`/`uend` as the raw doubles but stores the knot
    # vector as `float` — both conversions are preserved bit-for-bit.
    m.curve_geometry[t]=(knots=Float64.(Float32.(knot_vector)),
                         deg=degree,ubeg=first(knot_vector),
                         uend=last(knot_vector))
    return t
end

# Resolved spline-family geometry — control points re-read from `m.points` on
# every call exactly like `_arc_geometry`, so transforms and point moves need
# no cache invalidation.
function _spline_geometry(m::GeoModel, tag::Int, caller::AbstractString)
    kind=_curve_type(m,tag)
    kind in _SPLINE_CURVE_TYPES || throw(ArgumentError(
        "$caller: Curve[$tag] is a $(_curve_type_name(kind)), not a " *
        "spline-family curve"))
    cptags=get(m.curve_control_points,tag,nothing)
    cptags===nothing && throw(ErrorException(
        "$caller: Curve[$tag] has no spline control points; rebuild the model"))
    cps=NTuple{3,Float64}[]
    sizehint!(cps,length(cptags))
    for p in cptags
        haskey(m.points,p) || throw(ArgumentError(
            "$caller: Curve[$tag] references unknown Point[$p]"))
        push!(cps,m.points[p])
    end
    record=get(m.curve_geometry,tag,nothing)
    if kind===:nurbs
        (record===nothing || !hasproperty(record,:knots)) && throw(
            ErrorException(
                "$caller: Curve[$tag] has no Nurbs knot record; rebuild the " *
                "model"))
        return (kind=kind,cps=cps,closed=first(cptags)==last(cptags),
                knots=record.knots,deg=record.deg,
                ubeg=record.ubeg,uend=record.uend)
    end
    return (kind=kind,cps=cps,closed=first(cptags)==last(cptags),
            knots=nothing,deg=0,ubeg=0.0,uend=1.0)
end

# `InterpolateCatmullRom` — uniform-knot Catmull-Rom with extrapolated ghost
# endpoints and cyclic controls when `beg == end`.
function _catmull_rom_point(cps::Vector{NTuple{3,Float64}},closed::Bool,
                            u::Float64)
    n=length(cps)
    i=trunc(Int,(n-1)*u)
    i<0 && (i=0)
    i>=n-1 && (i=n-2)
    t1=i/(n-1)
    t2=(i+1)/(n-1)
    t=(u-t1)/(t2-t1)
    v1=cps[i+1];v2=cps[i+2]
    v0=if i==0
        closed ? cps[n-1] : (2.0*v1[1]-v2[1],2.0*v1[2]-v2[2],2.0*v1[3]-v2[3])
    else
        cps[i]
    end
    v3=if i==n-2
        closed ? cps[2] : (2.0*v2[1]-v1[1],2.0*v2[2]-v1[2],2.0*v2[3]-v1[3])
    else
        cps[i+3]
    end
    # `fma` placement replicates the contracted assembly of
    # `InterpolateCatmullRom` (`s` coefficients and the `s·v` accumulation);
    # verified bit-for-bit against the pinned binary on dense sweeps.
    t2_=t*t
    t3=t*t2_
    s0=fma(-0.5,t,fma(-0.5,t3,t2_))
    s1=fma(1.5,t3,-(2.5*t2_))+1.0
    s2=fma(0.5,t,fma(-1.5,t3,2.0*t2_))
    s3=fma(0.5,t3,-(0.5*t2_))
    return (fma(s3,v3[1],fma(s2,v2[1],fma(s0,v0[1],s1*v1[1]))),
            fma(s3,v3[2],fma(s2,v2[2],fma(s0,v0[2],s1*v1[2]))),
            fma(s3,v3[3],fma(s2,v2[3],fma(s0,v0[3],s1*v1[3]))))
end

# `InterpolateBezier` — De Casteljau on the full control list.
function _bezier_point(cps::Vector{NTuple{3,Float64}},u::Float64)
    buffer=copy(cps)
    n=length(buffer)
    w=1.0-u
    while n>1
        n-=1
        for i in 1:n
            c1=buffer[i];c2=buffer[i+1]
            # `fma` matches the contracted `(1-u)*c1 + u*c2` in the compiled
            # `InterpolateBezier` inner loop.
            buffer[i]=(fma(w,c1[1],u*c2[1]),fma(w,c1[2],u*c2[2]),
                       fma(w,c1[3],u*c2[3]))
        end
    end
    return buffer[1]
end

# `InterpolateUBS` — piecewise cubics; the stored interior matrix is
# `matbs/6.0` and the boundary spans use the hardcoded endpoint matrices,
# with the right-extremity reflection trick before matrix selection.
const _UBS_MAT2 = ((0.0,0.0,0.0,0.0),(0.0,0.0,0.0,0.0),
                   (-1.0,1.0,0.0,0.0),(1.0,0.0,0.0,0.0))
const _UBS_MAT3 = ((0.0,0.0,0.0,0.0),(1.0,-2.0,1.0,0.0),
                   (-2.0,2.0,0.0,0.0),(1.0,0.0,0.0,0.0))
const _UBS_MAT4 = ((-1.0,3.0,-3.0,1.0),(3.0,-6.0,3.0,0.0),
                   (-3.0,3.0,0.0,0.0),(1.0,0.0,0.0,0.0))
const _UBS_MAT5 = ((-1.0,7.0/4.0,-1.0,1.0/4.0),(3.0,-4.5,1.5,0.0),
                   (-3.0,3.0,0.0,0.0),(1.0,0.0,0.0,0.0))
const _UBS_MAT6_1 = ((-1.0,7.0/4.0,-11.0/12.0,2.0/12.0),(3.0,-4.5,1.5,0.0),
                     (-3.0,3.0,0.0,0.0),(1.0,0.0,0.0,0.0))
const _UBS_MAT6_2 = ((-1.0/4.0,7.0/12.0,-7.0/12.0,1.0/4.0),
                     (3.0/4.0,-5.0/4.0,1.0/2.0,0.0),
                     (-3.0/4.0,1.0/4.0,1.0/2.0,0.0),
                     (1.0/4.0,7.0/12.0,1.0/6.0,0.0))
const _UBS_MATEXT_1 = ((-1.0,7.0/4.0,-11.0/12.0,1.0/6.0),(3.0,-4.5,1.5,0.0),
                       (-3.0,3.0,0.0,0.0),(1.0,0.0,0.0,0.0))
# `mat7_2` and `matext_2` carry identical entries in `InterpolateUBS`, so the
# `iCurve == 1` boundary matrix is shared for `NbControlPoints >= 7`.
const _UBS_MAT7_2 = ((-1.0/4.0,7.0/12.0,-1.0/2.0,1.0/6.0),
                     (3.0/4.0,-5.0/4.0,1.0/2.0,0.0),
                     (-3.0/4.0,1.0/4.0,1.0/2.0,0.0),
                     (1.0/4.0,7.0/12.0,1.0/6.0,0.0))
const _UBS_MAT_INTERIOR = ((-1.0/6.0,3.0/6.0,-3.0/6.0,1.0/6.0),
                           (3.0/6.0,-6.0/6.0,3.0/6.0,0.0),
                           (-3.0/6.0,0.0,3.0/6.0,0.0),
                           (1.0/6.0,4.0/6.0,1.0/6.0,0.0))

# `InterpolateCubicSpline` — per-axis matrix blend `vec[i] = Σ_j mat[i,j]·v[j]`
# then `V = Σ_j T[j]·vec[j]` with T = (t³,t²,t,1). The `fma` chains replicate
# the compiled `+=` contractions (the `T[3] == 1` term folds to a plain add).
function _ubs_cubic(v::NTuple{4,NTuple{3,Float64}},t::Float64,mat)
    t2=t*t
    t3=t*t2
    return ntuple(3) do axis
        vec=ntuple(4) do i
            fma(mat[i][4],v[4][axis],
                fma(mat[i][3],v[3][axis],
                    fma(mat[i][2],v[2][axis],
                        fma(mat[i][1],v[1][axis],0.0))))
        end
        # The `V += T[j]*vec[j]` loop vectorizes in the pinned binary into
        # plain mul+add pairs — no fma.
        t3*vec[1]+t2*vec[2]+t*vec[3]+vec[4]
    end
end

function _ubs_point(cps::Vector{NTuple{3,Float64}},closed::Bool,u::Float64)
    n=length(cps)
    ncurves=max(n+(closed ? -1 : -3),1)
    seg=floor(Int,u*ncurves)
    icurve=seg==ncurves ? seg-1 : seg   # `iCurve == NbCurves` — u = 1
    t1=icurve/ncurves
    t2=(icurve+1)/ncurves
    t=(u-t1)/(t2-t1)
    v=ntuple(4) do i
        k=if closed
            mod(i-2+icurve,n-1)   # (iCurve-1+i) mod (N-1), 0-based
        else
            clamp(icurve+i-1,0,n-1)
        end
        cps[k+1]
    end
    mat=_UBS_MAT_INTERIOR
    if !closed
        mi=icurve
        if (n>6 && mi>=ncurves-2) || (n>4 && mi==ncurves-1)
            v=(v[4],v[3],v[2],v[1])
            t=1.0-t
            mi=ncurves-1-mi
        end
        if mi==0
            mat=n==2 ? _UBS_MAT2 : n==3 ? _UBS_MAT3 : n==4 ? _UBS_MAT4 :
                n==5 ? _UBS_MAT5 : n==6 ? _UBS_MAT6_1 : _UBS_MATEXT_1
        elseif mi==1
            mat=n==6 ? _UBS_MAT6_2 : _UBS_MAT7_2  # mat7_2 == matext_2
        end
    end
    return _ubs_cubic(v,t,mat)
end

# `findSpan`/`basisFuns`/`InterpolateNurbs` — non-rational B-spline evaluation
# on the float32-rounded knot vector (`Curve::k` is `float*`).
# `findSpan` — the upstream check is `u <= U[0]`, which never exits the binary
# search when an unclamped vector leaves `U[0] < u < U[deg]`. The `U[deg]`
# bound below returns the first span for that input instead; on every
# (clamped) vector where `U[0] == U[deg]` — and on every input where upstream
# terminates — the result is identical.
function _nurbs_findspan(u::Float64,deg::Int,n::Int,knots::Vector{Float64})
    u>=knots[n+1] && return n-1
    u<=knots[deg+1] && return deg
    low=deg;high=n+1;mid=(low+high)÷2
    while u<knots[mid+1] || u>=knots[mid+2]
        if u<knots[mid+1]
            high=mid
        else
            low=mid
        end
        mid=(low+high)÷2
    end
    return mid
end

function _nurbs_basisfuns(u::Float64,span0::Int,deg::Int,
                          knots::Vector{Float64})
    left=Vector{Float64}(undef,deg+1)
    right=Vector{Float64}(undef,deg+1)
    N=Vector{Float64}(undef,deg+1)
    N[1]=1.0
    for j in 1:deg
        left[j+1]=u-knots[span0-j+2]
        right[j+1]=knots[span0+j+1]-u
        saved=0.0
        for r in 0:j-1
            temp=N[r+1]/(right[r+2]+left[j-r+1])
            # `fma` matches the compiled `saved + right[r+1]*temp`.
            N[r+1]=fma(right[r+2],temp,saved)
            saved=left[j-r+1]*temp
        end
        N[j+1]=saved
    end
    return N
end

function _nurbs_point(g,u::Float64)
    cps=g.cps;knots=g.knots;deg=g.deg
    span0=_nurbs_findspan(u,deg,length(cps),knots)
    N=_nurbs_basisfuns(u,span0,deg,knots)
    # The `p += Nb[i]*v` accumulation is `fmla`-fused in the pinned binary.
    x=y=z=0.0
    for i in deg:-1:0
        p=cps[span0-deg+i+1]
        x=fma(N[i+1],p[1],x)
        y=fma(N[i+1],p[2],y)
        z=fma(N[i+1],p[3],z)
    end
    return (x,y,z)
end

# `InterpolateCurve` dispatch on `c->Typ` for the spline family.
function _spline_point(g,u::Float64)
    kind=g.kind
    if kind===:spline
        return _catmull_rom_point(g.cps,g.closed,u)
    elseif kind===:bezier
        return _bezier_point(g.cps,u)
    elseif kind===:bspline
        return _ubs_point(g.cps,g.closed,u)
    end
    return _nurbs_point(g,u)
end

# `InterpolateCurve(der=1/2)` — finite differences with the hardcoded [0,1]
# bound checks (`eps1 = u<eps ? 0 : eps`), which apply verbatim even when a
# Nurbs parameter range is not [0,1].
function _spline_first_derivative_fd(g,u::Float64)
    eps=1e-8
    eps1=u<eps ? 0.0 : eps
    eps2=u>1.0-eps ? 0.0 : eps
    p0=_spline_point(g,u-eps1)
    p1=_spline_point(g,u+eps2)
    inv=1.0/(eps1+eps2)
    return ((p1[1]-p0[1])*inv,(p1[2]-p0[2])*inv,(p1[3]-p0[3])*inv)
end

function _spline_second_derivative_fd(g,u::Float64)
    eps=1e-8
    eps1=u<eps ? 0.0 : eps
    eps2=u>1.0-eps ? 0.0 : eps
    d0=_spline_first_derivative_fd(g,u-eps1)
    d1=_spline_first_derivative_fd(g,u+eps2)
    inv=1.0/(eps1+eps2)
    return ((d1[1]-d0[1])*inv,(d1[2]-d0[2])*inv,(d1[3]-d0[3])*inv)
end

# `GEdge::curvature` — same FD derivatives and cross/norm chain as arcs.
function _spline_curvature(g,u::Float64)
    d1=_spline_first_derivative_fd(g,u)
    d2=_spline_second_derivative_fd(g,u)
    cr=_occ_cross(d1,d2)
    return sqrt(_sqlen(cr))*_gm_pow(1.0/sqrt(_sqlen(d1)),3.0)
end

# `GEdge::refineProjection` on a spline — damped Newton on the FD derivative,
# parameter clamped to `parBounds` = [ubeg,uend].
function _spline_refine_projection(
    g,q::NTuple{3,Float64},u::Float64,relax::Float64,tol::Float64,lc::Float64)
    maxDist=tol*lc
    dPQ=_arc_sub(_spline_point(g,u),q)
    err=sqrt(_sqlen(dPQ))
    iter=0
    while (iter+=1)<=25 && err>maxDist
        der=_spline_first_derivative_fd(g,u)
        du=_dot3(dPQ,der)/_dot3(der,der)
        du<tol && sqrt(_sqlen(dPQ))>maxDist && (du=1.0)
        unew=fma(-relax,du,u)
        # `std::min(uMax,std::max(uMin,uNew))` — NaN input clamps to `uMin`.
        unew=g.ubeg<unew ? unew : g.ubeg
        u=unew<g.uend ? unew : g.uend
        dPQ=_arc_sub(_spline_point(g,u),q)
        err=sqrt(_sqlen(dPQ))
    end
    return err<=maxDist,u,err
end

@inline _spline_dist(g,q::NTuple{3,Float64},t::Float64) =
    sqrt(_sqlen(_arc_sub(q,_spline_point(g,t))))

# `GEdge::XYZToU` on a spline: 21 evenly spaced seeds over [ubeg,uend] through
# `refineProjection`, then a relaxed retry from the last seed's refined value.
function _spline_xyz_to_u(
    g,q::NTuple{3,Float64},lc::Float64;relax::Float64=1.0)
    errors=Dict{Float64,Float64}()
    umin=g.ubeg;umax=g.uend
    u_try=umin; err=Inf
    step=(umax-umin)/20.0
    for i in 0:20
        u_try=fma(step,Float64(i),umin)
        ok,u_try,err=_spline_refine_projection(g,q,u_try,relax,1e-8,lc)
        ok && return true,u_try
        errors[err]=u_try
    end
    if relax>0.1
        ok,u_try=_spline_xyz_to_u(g,q,lc;relax=0.75*relax)
        ok && return true,u_try
        errors[_spline_dist(g,q,u_try)]=u_try
    end
    return false,errors[minimum(keys(errors))]
end

# `goldenSectionSearch` on a spline — same recursive minimizer over [ubeg,uend].
function _spline_golden_section(
    g,q::NTuple{3,Float64},x1::Float64,x2::Float64,x3::Float64,tau::Float64)
    golden2=2.0-(1.0+sqrt(5.0))/2.0
    x4=fma(golden2,x3-x2,x2)
    abs(x3-x1)<tau*(abs(x2)+abs(x4)) && return (x3+x1)/2
    return _spline_dist(g,q,x4)<_spline_dist(g,q,x2) ?
        _spline_golden_section(g,q,x2,x4,x3,tau) :
        _spline_golden_section(g,q,x4,x2,x1,tau)
end

# `GEdge::closestPoint` on a spline: 100-sample scan over [ubeg,uend] for the
# bracket, then the recursive golden-section search. Returns (t, point).
function _spline_closest(g,q::NTuple{3,Float64})
    tmin=g.ubeg;tmax=g.uend
    dt=(tmax-tmin)/99.0
    dmin=1e22; topt=tmin
    for i in 0:99
        t=fma(Float64(i),dt,tmin)
        d=_spline_dist(g,q,t)
        d<dmin && (topt=t; dmin=d)
    end
    t=topt==tmin ?
        _spline_golden_section(g,q,topt,topt+dt/2.0,topt+dt,1e-9) :
        topt==tmax ?
        _spline_golden_section(g,q,topt-dt,topt-dt/2.0,topt,1e-9) :
        _spline_golden_section(g,q,topt-dt,topt,topt+dt,1e-9)
    return t,_spline_point(g,t)
end

# `GEdge::bounds(fast=false)` — a 10-sample scan of `point(t)` over
# `parBounds` = [ubeg,uend]. Gmsh's reported curve box is exactly this scan
# (it neither includes control points nor solves extrema), so the box is a
# bit-exact port rather than an exact-extrema computation.
function _spline_bounding_box(g)
    lo=fill(Inf,3);hi=fill(-Inf,3)
    for i in 0:9
        t=fma(Float64(i)/9.0,g.uend-g.ubeg,g.ubeg)
        p=_spline_point(g,t)
        for axis in 1:3
            # `SBoundingBox3d::operator+=` extends through `<`/`>` tests, so a
            # NaN coordinate fails both comparisons and the sample is skipped
            # rather than polluting the box — plain `min`/`max` would
            # propagate the NaN instead.
            p[axis]<lo[axis] && (lo[axis]=p[axis])
            p[axis]>hi[axis] && (hi[axis]=p[axis])
        end
    end
    return (lo[1],lo[2],lo[3],hi[1],hi[2],hi[3])
end
