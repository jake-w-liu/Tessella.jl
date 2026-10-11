# Windows native GEO frame arithmetic used by actual high-order placement.
# Existing geometry carriers keep their established arithmetic policy.
@inline function _model_high_order_arc_normalize(vector::NTuple{3,Float64})
    magnitude=sqrt(vector[1]*vector[1]+vector[2]*vector[2]+vector[3]*vector[3])
    magnitude==0.0 && return vector
    reciprocal=1.0/magnitude
    return ntuple(axis->vector[axis]*reciprocal,3)
end
@inline _model_high_order_arc_cross(a::NTuple{3,Float64},b::NTuple{3,Float64})=
    (a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
@inline _model_high_order_arc_dot(a::NTuple{3,Float64},b::NTuple{3,Float64})=
    a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
function _model_high_order_arc_frame(cps::Vector{NTuple{3,Float64}}, kind::Symbol,
                       stored_n::NTuple{3,Float64}, tag::Int,
                       caller::AbstractString)
    length(cps) in (3,4) || throw(ErrorException(
        "$caller: Curve[$tag] has $(length(cps)) control points; " *
        "an arc needs 3 or 4"))
    has_major=length(cps)==4
    P0,C=cps[1],cps[2]
    P2=has_major ? cps[4] : cps[3]
    raw0=_arc_sub(P0,C); raw2=_arc_sub(P2,C)
    u=_model_high_order_arc_normalize(raw0); e=_model_high_order_arc_normalize(raw2)
    n=_model_high_order_arc_cross(u,e)
    if sqrt(n[1]*n[1]+n[2]*n[2]+n[3]*n[3])<1e-15
        n=_model_high_order_arc_normalize(stored_n)
    else
        n=_model_high_order_arc_normalize(n)
        # Gmsh rechecks the normalized vector before falling back to the
        # stored `Plane{..}` normal.
        if abs(n[1])<1e-5 && abs(n[2])<1e-5 && abs(n[3])<1e-5
            n=_model_high_order_arc_normalize(stored_n)
        end
    end
    v=_model_high_order_arc_normalize(_model_high_order_arc_cross(n,u))
    x0=_model_high_order_arc_dot(raw0,u); y0=_model_high_order_arc_dot(raw0,v)
    x2=_model_high_order_arc_dot(raw2,u); y2=_model_high_order_arc_dot(raw2,v)
    R=sqrt(x0*x0+y0*y0); R2=sqrt(x2*x2+y2*y2)
    (R==0.0 || R2==0.0) && throw(ArgumentError(
        "$caller: zero radius in circle or ellipse with tag $tag"))
    if !has_major && abs((R-R2)/(R+R2))>0.1
        throw(ArgumentError(
            "$caller: control points of circle with tag $tag are not " *
            "cocircular (R1=$R, R2=$R2)"))
    end
    local A1, A3, A4, f1, f2
    if has_major
        raw3=_arc_sub(cps[3],C)
        x3=_model_high_order_arc_dot(raw3,u); y3=_model_high_order_arc_dot(raw3,v)
        A4=_angle_02pi(_myatan2(y3,x3))
        s,c=_model_high_order_curve_sincos(A4)
        x1=x0*c+y0*s; y1=-x0*s+y0*c
        xe=x2*c+y2*s; ye=-x2*s+y2*c
        # Literal sys2x2 solve, with pre-rounded squared entries.
        det=(x1*x1)*(ye*ye)-(xe*xe)*(y1*y1)
        matnorm=(x1*x1)^2+(ye*ye)^2+(y1*y1)^2+(xe*xe)^2
        s0=s1=0.0
        if !(matnorm==0.0 || abs(det)/matnorm<1e-16)
            ud=1.0/det
            s0=(1.0*(ye*ye)-(y1*y1)*1.0)*ud
            s1=((x1*x1)*1.0-(xe*xe)*1.0)*ud
        end
        (s0<=0.0 || s1<=0.0) && throw(ArgumentError(
            "$caller: ellipse with tag $tag is wrong"))
        f1=sqrt(1.0/s0); f2=sqrt(1.0/s1)
        A1= x1<0.0 ? -_myasin(y1/f2)+A4+π : _myasin(y1/f2)+A4
        A3= xe<0.0 ? -_myasin(ye/f2)+A4+π : _myasin(ye/f2)+A4
    else
        A1=_myatan2(y0,x0); A3=_myatan2(y2,x2); A4=0.0; f1=f2=R
    end
    A1=_angle_02pi(A1); A3=_angle_02pi(A3)
    A1>=A3 && (A3+=2π)
    span=A3-A1
    span>1.01π && throw(ArgumentError(
        "$caller: circle or ellipse arc $tag greater than Pi (angle=$span)"))
    return (center=C, u=u, v=v, n=n, t1=A1, t2=A3, f1=f1, f2=f2, incl=A4)
end

function _model_high_order_arc_geometry(model::GeoModel,tag::Int,caller)
    @static if !Sys.iswindows()
        return _arc_geometry(model,tag,caller)
    end
    controls=get(model.curve_control_points,tag,nothing)
    controls===nothing && throw(ArgumentError("$caller: Curve[$tag] has no arc control points"))
    coordinates=NTuple{3,Float64}[model.points[point] for point in controls]
    stored=get(model.curve_geometry,tag,nothing)
    normal=stored===nothing ? (0.0,0.0,1.0) : stored.n
    return _model_high_order_arc_frame(coordinates,_curve_type(model,tag),normal,tag,caller)
end
