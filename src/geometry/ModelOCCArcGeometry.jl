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
