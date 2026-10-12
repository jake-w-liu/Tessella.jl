# Native center-form OCC CircleArc geometry. The plane fixes the circle axis;
# gce_MakeCirc chooses its own reference direction, then the trim reverses it.
function _model_occ_circle_arc_parameter(center,X,n,point)
    raw=_arc_sub(point,center)
    _sqlen(raw)<floatmin(Float64) && return 0.0
    projected=_occ_crosscross(n,raw,n)
    _sqlen(projected)<floatmin(Float64) && return 0.0
    angle=_occ_angle_with_ref(X,projected,n)
    return angle < -1e-16 ? angle+2π : angle<0.0 ? 0.0 : angle
end

function _model_occ_circle_arc_geometry(first::NTuple{3,Float64},
                                         center::NTuple{3,Float64},
                                         last::NTuple{3,Float64})
    v1=_arc_sub(center,first);v2=_arc_sub(last,first)
    normal=_occ_cross(v1,v2)
    plane_normal=_occ_modulus(normal)<floatmin(Float64) ?
        (0.0,0.0,1.0) : _occ_dir(normal)
    X,original_Y=_occ_ax2_frame(plane_normal)
    Y=ntuple(axis->-original_Y[axis],3)
    n=ntuple(axis->-plane_normal[axis],3)
    radius=_occ_modulus(_arc_sub(first,center))
    lower=_model_occ_circle_arc_parameter(center,X,n,first)
    upper=_model_occ_circle_arc_parameter(center,X,n,last)
    lower==upper && return nothing,"OpenCASCADE exception Geom_TrimmedCurve::U1 == U2"
    upper<lower && (upper+=2π)
    # BRepBuilderAPI_MakeEdge checks the supplied vertices against the circle.
    end_radius=_occ_modulus(_arc_sub(last,center))
    radius>0.0 && isfinite(radius) && abs(end_radius-radius)<=_OCC_CONFUSION ||
        return nothing,"Could not create circle arc"
    native=_model_native_occ_circle_geometry(first,center,last)
    native.t0==native.t1 && return nothing,"OpenCASCADE exception Geom_TrimmedCurve::U1 == U2"
    return (occ=:circle,center=center,n=n,X=X,Y=Y,r=radius,t0=native.t0,t1=native.t1,
            native_geometry=native),nothing
end

# `gce_MakeCirc(P1,P2,P3)` + `ElCLib::Parameter` trim for the
# `addCircleArc(..., center=false)` form: the circle through all three
# points, trimmed start→end through the middle. The plane axis
# (middle-start)×(end-start) orients the parameterization so the middle
# point's parameter always lies in (0, parameter(end)) — upstream's
# `Sense` computation on `Geom_TrimmedCurve` is therefore always the
# forward trim.
function _model_occ_circle_arc_through_geometry(first::NTuple{3,Float64},
                                                middle::NTuple{3,Float64},
                                                last::NTuple{3,Float64})
    v12=_arc_sub(middle,first); v13=_arc_sub(last,first)
    w=_occ_cross(v12,v13)
    wmod2=_dot3(w,w)
    wmod2<=_OCC_CONFUSION*_OCC_CONFUSION &&
        return nothing,"Could not build circle through three points"
    n=_occ_dir(w)
    # Circumcenter of the triangle in its plane:
    # C = first + (|v12|²·(v13×w) − |v13|²·(v12×w)) / (2·|w|²)
    o=_occ_mul(_arc_sub(_occ_mul(_occ_cross(v13,w),_dot3(v12,v12)),
                        _occ_mul(_occ_cross(v12,w),_dot3(v13,v13))),
               0.5/wmod2)
    center=_add3(first,o)
    r=_occ_modulus(_arc_sub(center,first))
    r>_OCC_CONFUSION ||
        return nothing,"Could not build circle through three points"
    X=_occ_dir(_arc_sub(first,center))
    Y=_occ_crossed(n,X)
    upper=_model_occ_circle_arc_parameter(center,X,n,last)
    upper<=0.0 &&
        return nothing,"OpenCASCADE exception Geom_TrimmedCurve::U1 == U2"
    return (occ=:circle,center=center,n=n,X=X,Y=Y,r=r,t0=0.0,t1=upper),nothing
end

# OCCT 7.8.1's center-form circle construction and Geom_Conic::Reverse.
# The latter calls gp_Ax2::SetDirection, projecting X a second time. Keep
# this arithmetic with the actual constructor; older primitive records retain
# the arithmetic used by their original construction and transform protocols.
@inline _model_native_occ_modulus(v)=sqrt((v[1]*v[1]+v[2]*v[2])+v[3]*v[3])
@inline function _model_native_occ_direction(v)
    length=_model_native_occ_modulus(v)
    return (v[1]/length,v[2]/length,v[3]/length)
end
@inline _model_native_occ_cross(a,b)=
    (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
@inline function _model_native_occ_setx(n,d)
    X=_model_native_occ_direction(_model_native_occ_cross(
        n,_model_native_occ_cross(d,n)))
    Y=_model_native_occ_direction(_model_native_occ_cross(n,X))
    return X,Y
end

@inline _model_native_occ_dot(a,b)=(a[1]*b[1]+a[2]*b[2])+a[3]*b[3]

function _model_native_occ_circle_parameter(center,X,n,point)
    raw=_arc_sub(point,center)
    _model_native_occ_dot(raw,raw)<floatmin(Float64) && return 0.0
    projected=_model_native_occ_cross(n,_model_native_occ_cross(raw,n))
    _model_native_occ_dot(projected,projected)<floatmin(Float64) && return 0.0
    # X and n are already gp_Dir values; only the projected XYZ converts.
    direction=_model_native_occ_direction(projected)
    product=_model_native_occ_cross(X,direction)
    cosine=_model_native_occ_dot(X,direction)
    sine=_model_native_occ_modulus(product)
    angle=(-0.70710678118655<cosine<0.70710678118655) ?
        _gm_acos(cosine) : cosine<0.0 ? π-_gm_asin(sine) : _gm_asin(sine)
    _model_native_occ_dot(product,n)<0.0 && (angle=-angle)
    return angle < -1e-16 ? angle+2π : angle<0.0 ? 0.0 : angle
end

function _model_native_occ_adjust_circle_period(lower,upper,precision)
    period=2π
    lower-=floor(lower/period)*period
    period-lower<precision && (lower-=period)
    upper-=floor((upper-lower)/period)*period
    upper-lower<precision && (upper+=period)
    return lower,upper
end

function _model_native_occ_circle_geometry(first,center,last)
    normal=_model_native_occ_cross(_arc_sub(center,first),_arc_sub(last,first))
    n=_model_native_occ_modulus(normal)<floatmin(Float64) ?
        (0.0,0.0,1.0) : _model_native_occ_direction(normal)
    a,b,c=n; aa,bb,cc=abs(a),abs(b),abs(c)
    raw=if bb<=aa && bb<=cc
        aa>cc ? (-c,0.0,a) : (c,0.0,-a)
    elseif aa<=bb && aa<=cc
        bb>cc ? (0.0,-c,b) : (0.0,c,-b)
    else
        aa>bb ? (-b,a,0.0) : (b,-a,0.0)
    end
    X,_=_model_native_occ_setx(n,_model_native_occ_direction(raw))
    alpha1=_model_native_occ_circle_parameter(center,X,n,first)
    alpha2=_model_native_occ_circle_parameter(center,X,n,last)
    if alpha1==alpha2
        lower=upper=alpha1
    else
        # Geom_TrimmedCurve adjusts before Reverse. MakeEdge projects the
        # start vertex to the reversed last bound, then adjusts once more.
        trim1,trim2=_model_native_occ_adjust_circle_period(
            alpha1,alpha2,min(abs(alpha2-alpha1)/2,1e-9))
        lower,upper=_model_native_occ_adjust_circle_period(2π-trim1,2π-trim2,1e-9)
    end
    n=(-n[1],-n[2],-n[3])
    X,Y=_model_native_occ_setx(n,X)
    radius=_model_native_occ_modulus(_arc_sub(first,center))
    return (occ=:circle,center=center,n=n,X=X,Y=Y,r=radius,t0=lower,t1=upper,
            native_arithmetic=true)
end

@inline function _model_native_occ_circle_point(g,u)
    sine,cosine=_model_high_order_curve_sincos(u)
    x,y=g.r*cosine,g.r*sine
    return (x*g.X[1]+y*g.Y[1]+g.center[1],
            x*g.X[2]+y*g.Y[2]+g.center[2],
            x*g.X[3]+y*g.Y[3]+g.center[3])
end

@inline function _model_native_occ_circle_derivative(g,u)
    sine,cosine=_model_high_order_curve_sincos(u)
    x,y=g.r*cosine,g.r*sine
    return ((-y)*g.X[1]+x*g.Y[1],(-y)*g.X[2]+x*g.Y[2],
            (-y)*g.X[3]+x*g.Y[3])
end

@inline function _model_native_occ_circle_second_derivative(g,u)
    sine,cosine=_model_high_order_curve_sincos(u)
    x,y=g.r*cosine,g.r*sine
    return ((-x)*g.X[1]+(-y)*g.Y[1],(-x)*g.X[2]+(-y)*g.Y[2],
            (-x)*g.X[3]+(-y)*g.Y[3])
end

function _model_native_occ_circle_curvature(g,u)
    d1=_model_native_occ_circle_derivative(g,u)
    d2=_model_native_occ_circle_second_derivative(g,u)
    dd1=(d1[1]*d1[1]+d1[2]*d1[2])+d1[3]*d1[3]
    dd2=(d2[1]*d2[1]+d2[2]*d2[2])+d2[3]*d2[3]
    epsilon=1e-15; tolerance=epsilon*epsilon
    dd2<=tolerance && return epsilon
    cross=_model_native_occ_cross(d1,d2)
    numerator=(cross[1]*cross[1]+cross[2]*cross[2])+cross[3]*cross[3]
    curvature=numerator/dd1/dd2<=tolerance ? 0.0 : sqrt(numerator)/dd1/sqrt(dd1)
    return curvature<=epsilon ? epsilon : curvature
end

# `gp_Dir::IsParallel(other, AngularTolerance)` — |angle| ≤ tol or
# |π − angle| ≤ tol; `u` is unit.
function _model_occ_dir_parallel(u,v,tol)
    ang=_gm_atan2(_occ_modulus(_occ_cross(u,v)),_dot3(u,v))
    return ang<=tol || π-ang<=tol
end

# `OCC_Internals::addEllipseArc` — the ellipse through the start and end
# points whose major axis direction is center→major. `u` is that direction;
# `v` is the first endpoint's off-axis component direction (the other's
# when the first is parallel). The conic's semi-axes solve the two-point
# system, then the `gp_Ax2` frame and the `ElCLib::Parameter` trim pick the
# sub-π traversal — forward when `mod2π(α2−α1) < π`, the reversed frame
# (negated Y and n, mirrored parameter range) otherwise.
function _model_occ_ellipse_arc_geometry(first::NTuple{3,Float64},
                                         center::NTuple{3,Float64},
                                         major::NTuple{3,Float64},
                                         last::NTuple{3,Float64})
    x1=_arc_sub(first,center); x2=_arc_sub(last,center)
    udir=_arc_sub(major,center)
    _occ_modulus(udir)<=floatmin(Float64) &&
        return nothing,"OpenCASCADE exception gp_Dir() : input vector has null modulus",nothing
    u=_occ_dir(udir)
    vraw=if !_model_occ_dir_parallel(u,x1,1e-6)
        _arc_sub(x1,_occ_mul(u,_dot3(x1,u)))
    elseif !_model_occ_dir_parallel(u,x2,1e-6)
        _arc_sub(x2,_occ_mul(u,_dot3(x2,u)))
    else
        return nothing,"Cannot create ellipse arc with start and end point " *
                       "on the major axis",nothing
    end
    _occ_modulus(vraw)<=floatmin(Float64) &&
        return nothing,"OpenCASCADE exception gp_Dir() : input vector has null modulus",nothing
    v=_occ_dir(vraw)
    x1u=_dot3(x1,u)^2; x1v=_dot3(x1,v)^2
    x2u=_dot3(x2,u)^2; x2v=_dot3(x2,v)^2
    (abs(x1u-x2u)<_OCC_CONFUSION || abs(x1v-x2v)<_OCC_CONFUSION) &&
        return nothing,"Cannot create ellipse arc with start and end point " *
                       "symmetric with respect to major or minor axis",nothing
    a2=(x1v*x2u-x1u*x2v)/(x1v-x2v)
    b2=(x1u*x2v-x1v*x2u)/(x1u-x2u)
    (a2<=0.0 || b2<=0.0) &&
        return nothing,"Invalid radii during creation of ellipse arc",nothing
    warn=false
    if a2>=b2
        a,b=sqrt(a2),sqrt(b2)
        n=_occ_crossed(u,v); X=u; Y=v
    else
        warn=true
        a,b=sqrt(b2),sqrt(a2)
        n=_occ_crossed(v,u); X=v; Y=u
    end
    alpha1=_gm_atan2(_dot3(x1,Y)/b,_dot3(x1,X)/a)
    alpha2=_gm_atan2(_dot3(x2,Y)/b,_dot3(x2,X)/a)
    span=_angle_02pi(alpha2-alpha1)
    if span<π
        t0,t1=alpha1,alpha1+span
    else
        # The ≤π traversal is clockwise: flip Y and n so the stored frame's
        # parameter mirrors the ellipse's — p′(u) = p(−u) — then re-anchor
        # the range so u(t0) = −α1 maps back to the start point.
        span=_angle_02pi(alpha1-alpha2)
        Y=ntuple(i->-Y[i],3); n=ntuple(i->-n[i],3)
        t0,t1=-alpha1,-alpha1+span
    end
    return (occ=:ellipse,center=center,n=n,X=X,Y=Y,rx=a,ry=b,t0=t0,t1=t1),
           nothing,(warn ? "Major radius smaller than minor radius" : nothing)
end
