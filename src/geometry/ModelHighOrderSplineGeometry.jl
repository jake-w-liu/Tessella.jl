# Literal native GEO spline arithmetic for high-order placement.
# Established geometry renderers retain their own compiler-contraction policy.
function _model_high_order_catmull_rom_point(cps::Vector{NTuple{3,Float64}},closed::Bool,
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
    t2_=t*t
    t3=t*t2_
    s0=((-0.5)*(t)+(((-0.5)*(t3)+(t2_))))
    s1=((1.5)*(t3)+(-(2.5*t2_)))+1.0
    s2=((0.5)*(t)+(((-1.5)*(t3)+(2.0*t2_))))
    s3=((0.5)*(t3)+(-(0.5*t2_)))
    return (((s3)*(v3[1])+(((s2)*(v2[1])+(((s0)*(v0[1])+(s1*v1[1])))))),
            ((s3)*(v3[2])+(((s2)*(v2[2])+(((s0)*(v0[2])+(s1*v1[2])))))),
            ((s3)*(v3[3])+(((s2)*(v2[3])+(((s0)*(v0[3])+(s1*v1[3])))))))
end

function _model_high_order_bezier_point(cps::Vector{NTuple{3,Float64}},u::Float64)
    buffer=copy(cps)
    n=length(buffer)
    w=1.0-u
    while n>1
        n-=1
        for i in 1:n
            c1=buffer[i];c2=buffer[i+1]
            buffer[i]=(((w)*(c1[1])+(u*c2[1])),((w)*(c1[2])+(u*c2[2])),
                       ((w)*(c1[3])+(u*c2[3])))
        end
    end
    return buffer[1]
end



function _model_high_order_ubs_cubic(v::NTuple{4,NTuple{3,Float64}},t::Float64,mat)
    t2=t*t
    t3=t*t2
    return ntuple(3) do axis
        vec=ntuple(4) do i
            ((mat[i][4])*(v[4][axis])+(
                ((mat[i][3])*(v[3][axis])+(
                    ((mat[i][2])*(v[2][axis])+(
                        ((mat[i][1])*(v[1][axis])+(0.0))))))))
        end
        t3*vec[1]+t2*vec[2]+t*vec[3]+vec[4]
    end
end

function _model_high_order_ubs_point(cps::Vector{NTuple{3,Float64}},closed::Bool,u::Float64)
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
    return _model_high_order_ubs_cubic(v,t,mat)
end

function _model_high_order_nurbs_findspan(u::Float64,deg::Int,n::Int,knots::Vector{Float64})
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

function _model_high_order_nurbs_basisfuns(u::Float64,span0::Int,deg::Int,
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
            N[r+1]=((right[r+2])*(temp)+(saved))
            saved=left[j-r+1]*temp
        end
        N[j+1]=saved
    end
    return N
end

function _model_high_order_nurbs_point(g,u::Float64)
    cps=g.cps;knots=g.knots;deg=g.deg
    span0=_model_high_order_nurbs_findspan(u,deg,length(cps),knots)
    N=_model_high_order_nurbs_basisfuns(u,span0,deg,knots)
    x=y=z=0.0
    for i in deg:-1:0
        p=cps[span0-deg+i+1]
        x=((N[i+1])*(p[1])+(x))
        y=((N[i+1])*(p[2])+(y))
        z=((N[i+1])*(p[3])+(z))
    end
    return (x,y,z)
end

function _model_high_order_spline_point(g,u::Float64)
    kind=g.kind
    if kind===:spline
        return _model_high_order_catmull_rom_point(g.cps,g.closed,u)
    elseif kind===:bezier
        return _model_high_order_bezier_point(g.cps,u)
    elseif kind===:bspline
        return _model_high_order_ubs_point(g.cps,g.closed,u)
    end
    return _model_high_order_nurbs_point(g,u)
end
