# Literal floating-point arithmetic used by the native Surface insertion
# algorithms. These operations are separate from the robust generic CDT
# predicates and constructions.
@inline function _gmsh_sys2x2(a,b,c,d,u,v)
    magnitude=a*a+d*d+b*b+c*c
    determinant=a*d-c*b
    if magnitude==0.0 || abs(determinant)/magnitude<1e-16
        return (0.0,0.0)
    end
    reciprocal=1.0/determinant
    return ((u*d-b*v)*reciprocal,(a*v-c*u)*reciprocal)
end

@inline function _gmsh_circumcenter_xy(a,b,c)
    x1=a[1];x2=b[1];x3=c[1]
    y1=a[2];y2=b[2];y3=c[2]
    denominator=2.0*(y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2))
    denominator==0.0 && return (-99999.0,-99999.0)
    a1=x1*x1+y1*y1
    a2=x2*x2+y2*y2
    a3=x3*x3+y3*y3
    return ((a1*(y3-y2)+a2*(y1-y3)+a3*(y2-y1))/denominator,
            (a1*(x2-x3)+a2*(x3-x1)+a3*(x1-x2))/denominator)
end

@inline _gmsh_dot3(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _gmsh_sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _gmsh_cross3(a,b)=(a[2]*b[3]-a[3]*b[2],
                          -a[1]*b[3]+a[3]*b[1],
                          a[1]*b[2]-a[2]*b[1])
@inline _gmsh_norm3(a)=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])

@inline function _gmsh_normalize3(a)
    magnitude=_gmsh_norm3(a)
    magnitude==0.0 && return a
    reciprocal=1.0/magnitude
    return (a[1]*reciprocal,a[2]*reciprocal,a[3]*reciprocal)
end

@inline function _gmsh_circumcenter_frame(a,b,c)
    v1=_gmsh_sub3(b,a)
    v2=_gmsh_sub3(c,a)
    vx=v1
    vy=v2
    vz=_gmsh_cross3(vx,vy)
    vy=_gmsh_cross3(vz,vx)
    vx=_gmsh_normalize3(vx)
    vy=_gmsh_normalize3(vy)
    vz=_gmsh_normalize3(vz)
    pb=(_gmsh_dot3(v1,vx),_gmsh_dot3(v1,vy))
    pc=(_gmsh_dot3(v2,vx),_gmsh_dot3(v2,vy))
    center=_gmsh_circumcenter_xy((0.0,0.0),pb,pc)
    return vx,vy,pb,pc,center
end

@inline function _gmsh_circumcenter_xyz(a,b,c)
    vx,vy,_,_,center=_gmsh_circumcenter_frame(a,b,c)
    return (a[1]+center[1]*vx[1]+center[2]*vy[1],
            a[2]+center[1]*vx[2]+center[2]*vy[2],
            a[3]+center[1]*vx[3]+center[2]*vy[3])
end

@inline function _gmsh_circumcenter_xyz_parameters(a,b,c)
    _,_,pb,pc,center=_gmsh_circumcenter_frame(a,b,c)
    return _gmsh_sys2x2(pb[1],pc[1],pb[2],pc[2],center...)
end

@inline function _gmsh_circumcenter_metric(pa,pb,pc,metric)
    a,b,d=metric
    m11=2.0*a*(pa[1]-pb[1])+2.0*b*(pa[2]-pb[2])
    m12=2.0*d*(pa[2]-pb[2])+2.0*b*(pa[1]-pb[1])
    m21=2.0*a*(pa[1]-pc[1])+2.0*b*(pa[2]-pc[2])
    m22=2.0*d*(pa[2]-pc[2])+2.0*b*(pa[1]-pc[1])
    rhs1=a*(pa[1]*pa[1]-pb[1]*pb[1])+
         d*(pa[2]*pa[2]-pb[2]*pb[2])+
         2.0*b*(pa[1]*pa[2]-pb[1]*pb[2])
    rhs2=a*(pa[1]*pa[1]-pc[1]*pc[1])+
         d*(pa[2]*pa[2]-pc[2]*pc[2])+
         2.0*b*(pa[1]*pa[2]-pc[1]*pc[2])
    center=_gmsh_sys2x2(m11,m12,m21,m22,rhs1,rhs2)
    x=center[1]-pa[1];y=center[2]-pa[2]
    radius2=x*x*a+y*y*d+2.0*x*y*b
    return center,radius2
end

@inline function _gmsh_metric_distance2(a,b,metric)
    x=a[1]-b[1];y=a[2]-b[2]
    return x*x*metric[1]+y*y*metric[3]+2.0*x*y*metric[2]
end

@inline function _gmsh_barycentric2(a,b,c,p)
    u,v=_gmsh_sys2x2(b[1]-a[1],c[1]-a[1],b[2]-a[2],c[2]-a[2],
                     p[1]-a[1],p[2]-a[2])
    return (1.0-u-v,u,v)
end

@inline _gmsh_flat_xyz(p)=(p[1],p[2],0.0)

@inline function _gmsh_qr_column_norm(qr,column,first)
    value=0.0
    for row in first:3
        value+=qr[row,column]*qr[row,column]
    end
    return sqrt(value)
end

# fullMatrix<double>::luSolve uses Eigen's column-pivoted Householder QR in
# the pinned native build. Keep its pivot, rank and column-update order for
# the three unknowns of intersectCurveSurface; singular QR is not a rejection.
function _gmsh_qr_solve3!(qr,rhs,solution,h,norms,direct,permutation)
    for column in 1:3
        norms[column]=_gmsh_qr_column_norm(qr,column,1)
        direct[column]=norms[column]
        permutation[column]=column
    end
    initial_norm=max(norms[1],norms[2],norms[3])
    threshold=(initial_norm*eps(Float64))^2/3.0
    downgrade=sqrt(eps(Float64))
    nonzero=3
    for k in 1:3
        best=k
        for column in k+1:3
            norms[column]>norms[best] && (best=column)
        end
        if nonzero==3 && norms[best]*norms[best]<threshold*(4-k)
            nonzero=k-1
        end
        if best!=k
            for row in 1:3
                qr[row,k],qr[row,best]=qr[row,best],qr[row,k]
            end
            norms[k],norms[best]=norms[best],norms[k]
            direct[k],direct[best]=direct[best],direct[k]
            permutation[k],permutation[best]=permutation[best],permutation[k]
        end
        tail=0.0
        for row in k+1:3
            tail+=qr[row,k]*qr[row,k]
        end
        first=qr[k,k]
        if tail<=floatmin(Float64)
            h[k]=0.0
            beta=first
            for row in k+1:3
                qr[row,k]=0.0
            end
        else
            beta=sqrt(first*first+tail)
            first>=0.0 && (beta=-beta)
            for row in k+1:3
                qr[row,k]/=first-beta
            end
            h[k]=(beta-first)/beta
        end
        qr[k,k]=beta
        for column in k+1:3
            if h[k]!=0.0
                temporary=0.0
                for row in k+1:3
                    temporary+=qr[row,k]*qr[row,column]
                end
                temporary+=qr[k,column]
                qr[k,column]-=h[k]*temporary
                for row in k+1:3
                    qr[row,column]-=(h[k]*qr[row,k])*temporary
                end
            end
            if norms[column]!=0.0
                value=abs(qr[k,column])/norms[column]
                value=(1.0+value)*(1.0-value)
                value=value<0.0 ? 0.0 : value
                relative=norms[column]/direct[column]
                reliability=value*(relative*relative)
                if reliability<=downgrade
                    direct[column]=_gmsh_qr_column_norm(qr,column,k+1)
                    norms[column]=direct[column]
                else
                    norms[column]*=sqrt(value)
                end
            end
        end
    end
    for k in 1:nonzero
        if k==3
            rhs[k]*=1.0-h[k]
        elseif h[k]!=0.0
            temporary=0.0
            for row in k+1:3
                temporary+=qr[row,k]*rhs[row]
            end
            temporary+=rhs[k]
            rhs[k]-=h[k]*temporary
            for row in k+1:3
                rhs[row]-=(h[k]*qr[row,k])*temporary
            end
        end
    end
    for row in nonzero:-1:1
        if rhs[row]!=0.0
            rhs[row]/=qr[row,row]
            for previous in 1:row-1
                rhs[previous]-=rhs[row]*qr[previous,row]
            end
        end
    end
    fill!(solution,0.0)
    for row in 1:nonzero
        solution[permutation[row]]=rhs[row]
    end
    return nothing
end

@inline function _gmsh_plane_circle_residual!(residual,x,middle,tangent,distance)
    cosine=cos(x[3])
    residual[1]=x[1]-(middle[1]+tangent[1]*cosine*distance)
    residual[2]=x[2]-(middle[2]+tangent[2]*cosine*distance)
    residual[3]=0.0-sin(x[3])*distance
    return nothing
end

function _gmsh_plane_circle_newton(fast,middle,tangent,distance)
    x=Float64[fast[1],fast[2],0.0]
    residual=zeros(Float64,3)
    perturbed=similar(residual)
    step=similar(residual)
    qr=Matrix{Float64}(undef,3,3)
    h=similar(residual)
    norms=similar(residual)
    direct=similar(residual)
    permutation=Vector{Int}(undef,3)
    for _ in 1:100
        _gmsh_norm3(x)>1e6 && return false,(x[1],x[2],x[3])
        _gmsh_plane_circle_residual!(residual,x,middle,tangent,distance)
        # Pinned newton_fd breaks and reports false on an exactly zero
        # residual. The caller publishes its Newton result only on true.
        all(iszero,residual) && break
        for column in 1:3
            difference=1e-4*abs(x[column])
            difference==0.0 && (difference=1e-4)
            x[column]+=difference
            _gmsh_plane_circle_residual!(perturbed,x,middle,tangent,distance)
            for row in 1:3
                qr[row,column]=(perturbed[row]-residual[row])/difference
            end
            x[column]-=difference
        end
        _gmsh_qr_solve3!(qr,residual,step,h,norms,direct,permutation)
        for row in 1:3
            x[row]-=step[row]
        end
        _gmsh_norm3(step)<1e-6 && return true,(x[1],x[2],x[3])
    end
    return false,(x[1],x[2],x[3])
end

function _gmsh_plane_circle_intersection(fast,middle,tangent,distance)
    residual=(fast[1]-(middle[1]+tangent[1]*distance),
              fast[2]-(middle[2]+tangent[2]*distance),0.0)
    _gmsh_norm3(residual)<distance*1e-8 && return fast
    success,point=_gmsh_plane_circle_newton(fast,middle,tangent,distance)
    return success ? (point[1],point[2]) : fast
end

# optimalPointFrontalB on a genuinely warped chart: the frontal circle lives
# in world space — radius `d` around the edge's world midpoint in the seed
# triangle's tangent plane, `middle + d·(n2·cos t + n1·sin t)` — and
# `intersectCurveSurface` Newton-solves `lift(u,v) == circle(t)` for
# `(u,v,t)` from the fast UV guess. The UV-plane proxy cannot express that
# circle: it applied the metric-normalized `distance` as a raw UV radius and
# overshot by the domain's metric scale. A degenerate tangent normal or a
# failed Newton keeps the fast point, matching upstream's `return true`.
@inline function _gmsh_surface_circle_residual!(residual,x,lift,middle,n2,n1,d)
    s=lift(x[1],x[2])
    cs=cos(x[3]);sn=sin(x[3])
    residual[1]=s[1]-(middle[1]+d*(n2[1]*cs+n1[1]*sn))
    residual[2]=s[2]-(middle[2]+d*(n2[2]*cs+n1[2]*sn))
    residual[3]=s[3]-(middle[3]+d*(n2[3]*cs+n1[3]*sn))
    return nothing
end

function _gmsh_surface_circle_intersection(fast,lift,wa,wb,wc)
    middle=((wa[1]+wb[1])*0.5,(wa[2]+wb[2])*0.5,(wa[3]+wb[3])*0.5)
    e=_gmsh_sub3(wb,wa)
    g=_gmsh_sub3(wc,middle)
    n1=_gmsh_cross3(e,g)
    n1n=_gmsh_norm3(n1)
    n1n<1e-12 && return fast
    n1=(n1[1]/n1n,n1[2]/n1n,n1[3]/n1n)
    n2=_gmsh_normalize3(_gmsh_cross3(n1,e))
    # The circle's world radius is the fast point's own world distance to the
    # edge midpoint — exactly the `L` optimalPointFrontal returned upstream.
    wf=lift(fast[1],fast[2])
    d=_gmsh_norm3(_gmsh_sub3(wf,middle))
    (isfinite(d) && d>0) || return fast
    x=Float64[fast[1],fast[2],0.0]
    residual=zeros(3);perturbed=similar(residual);step=similar(residual)
    qr=Matrix{Float64}(undef,3,3);h=similar(residual)
    norms=similar(residual);direct=similar(residual)
    permutation=Vector{Int}(undef,3)
    _gmsh_surface_circle_residual!(residual,x,lift,middle,n2,n1,d)
    _gmsh_norm3(residual)<d*1e-8 && return (fast[1],fast[2])
    for _ in 1:100
        _gmsh_norm3(x)>1e6 && return fast
        _gmsh_surface_circle_residual!(residual,x,lift,middle,n2,n1,d)
        all(iszero,residual) && break
        for column in 1:3
            difference=1e-4*abs(x[column])
            difference==0.0 && (difference=1e-4)
            x[column]+=difference
            _gmsh_surface_circle_residual!(perturbed,x,lift,middle,n2,n1,d)
            for row in 1:3
                qr[row,column]=(perturbed[row]-residual[row])/difference
            end
            x[column]-=difference
        end
        _gmsh_qr_solve3!(qr,residual,step,h,norms,direct,permutation)
        for row in 1:3
            x[row]-=step[row]
        end
        _gmsh_norm3(step)<1e-6 && return (x[1],x[2])
    end
    return fast
end
