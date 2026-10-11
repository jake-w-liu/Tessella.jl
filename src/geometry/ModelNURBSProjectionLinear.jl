# Small operation-owned linear solves used by the tensor projection root. The
# traversal and arithmetic follow OCCT 7.8.1 math_Recipes, math_Gauss and math_SVD.
_model_nurbs_root_sign(a,b)=b>=0.0 ? abs(a) : -abs(a)
function _model_nurbs_root_pythag(a,b)
    at,bt=abs(a),abs(b)
    if at>bt
        ct=bt/at
        return at*sqrt(1.0+ct*ct)
    elseif !iszero(bt)
        ct=at/bt
        return bt*sqrt(1.0+ct*ct)
    end
    return 0.0
end

function _model_nurbs_root_gauss(matrix,rhs,minimum_pivot)
    a=copy(matrix);b=copy(rhs);n=length(b)
    scales=zeros(n);indices=zeros(Int,n)
    for i in 1:n
        biggest=0.0
        for j in 1:n
            temp=abs(a[i,j])
            temp>biggest && (biggest=temp)
        end
        biggest<=minimum_pivot && return nothing
        scales[i]=1.0/biggest
    end
    for j in 1:n
        for i in 1:j-1
            total=a[i,j]
            for k in 1:i-1
                total-=a[i,k]*a[k,j]
            end
            a[i,j]=total
        end
        biggest=0.0;maximum_row=0
        for i in j:n
            total=a[i,j]
            for k in 1:j-1
                total-=a[i,k]*a[k,j]
            end
            a[i,j]=total
            value=scales[i]*abs(total)
            value<biggest && continue
            biggest=value;maximum_row=i
        end
        if j!=maximum_row
            for k in 1:n
                a[maximum_row,k],a[j,k]=a[j,k],a[maximum_row,k]
            end
            scales[maximum_row]=scales[j]
        end
        indices[j]=maximum_row
        abs(a[j,j])<=minimum_pivot && return nothing
        if j!=n
            reciprocal=1.0/a[j,j]
            for i in j+1:n
                a[i,j]*=reciprocal
            end
        end
    end
    first_nonzero=0
    for i in 1:n
        pivot=indices[i];total=b[pivot];b[pivot]=b[i]
        if first_nonzero!=0
            for j in first_nonzero:i-1
                total-=a[i,j]*b[j]
            end
        elseif !iszero(total)
            first_nonzero=i
        end
        b[i]=total
    end
    for i in n:-1:1
        total=b[i]
        for j in i+1:n
            total-=a[i,j]*b[j]
        end
        b[i]=total/a[i,i]
    end
    return b
end

function _model_nurbs_root_svd(matrix,rhs)
    a=copy(matrix);m,n=size(a)
    m>=n || throw(ArgumentError("tensor projection SVD requires at least as many equations as unknowns"))
    w=zeros(n);v=zeros(n,n);rv1=zeros(n)
    anorm=0.0;g=0.0;scale=0.0;l=0
    for i in 1:n
        l=i+1;rv1[i]=scale*g;g=0.0;s=0.0;scale=0.0
        for k in i:m
            aki=a[k,i]
            scale+=aki>0.0 ? aki : -aki
        end
        if !iszero(scale)
            for k in i:m
                a[k,i]/=scale;s+=a[k,i]*a[k,i]
            end
            f=a[i,i];g=-_model_nurbs_root_sign(sqrt(s),f)
            h=f*g-s;a[i,i]=f-g
            if i!=n
                for j in l:n
                    s=0.0
                    for k in i:m
                        s+=a[k,i]*a[k,j]
                    end
                    f=s/h
                    for k in i:m
                        a[k,j]+=f*a[k,i]
                    end
                end
            end
            for k in i:m
                a[k,i]*=scale
            end
        end
        w[i]=scale*g;g=0.0;s=0.0;scale=0.0
        if i!=n
            for k in l:n
                aik=a[i,k]
                scale+=aik>0.0 ? aik : -aik
            end
            if !iszero(scale)
                for k in l:n
                    a[i,k]/=scale;s+=a[i,k]*a[i,k]
                end
                f=a[i,l];g=-_model_nurbs_root_sign(sqrt(s),f)
                h=f*g-s;a[i,l]=f-g
                for k in l:n
                    rv1[k]=a[i,k]/h
                end
                if i!=m
                    for j in l:m
                        s=0.0
                        for k in l:n
                            s+=a[j,k]*a[i,k]
                        end
                        for k in l:n
                            a[j,k]+=s*rv1[k]
                        end
                    end
                end
                for k in l:n
                    a[i,k]*=scale
                end
            end
        end
        aw=abs(w[i]);ar=rv1[i]
        ar=ar>0.0 ? aw+ar : aw-ar
        anorm<ar && (anorm=ar)
    end
    for i in n:-1:1
        if i<n
            if !iszero(g)
                for j in l:n
                    v[j,i]=(a[i,j]/a[i,l])/g
                end
                for j in l:n
                    s=0.0
                    for k in l:n
                        s+=a[i,k]*v[k,j]
                    end
                    for k in l:n
                        v[k,j]+=s*v[k,i]
                    end
                end
            end
            for j in l:n
                v[i,j]=0.0;v[j,i]=0.0
            end
        end
        v[i,i]=1.0;g=rv1[i];l=i
    end
    for i in n:-1:1
        l=i+1;g=w[i]
        for j in l:n
            a[i,j]=0.0
        end
        if !iszero(g)
            g=1.0/g
            if i!=n
                for j in l:n
                    s=0.0
                    for k in l:m
                        s+=a[k,i]*a[k,j]
                    end
                    f=(s/a[i,i])*g
                    for k in i:m
                        a[k,j]+=f*a[k,i]
                    end
                end
            end
            for j in i:m
                a[j,i]*=g
            end
        else
            for j in i:m
                a[j,i]=0.0
            end
        end
        a[i,i]+=1.0
    end
    for k in n:-1:1
        for iteration in 1:30
            flag=true;nm=0;l=k
            for candidate in k:-1:1
                l=candidate;nm=l-1
                if abs(rv1[l])+anorm==anorm
                    flag=false;break
                end
                nm==0 && return nothing
                abs(w[nm])+anorm==anorm && break
            end
            if flag
                c=0.0;s=1.0
                for i in l:k
                    f=s*rv1[i]
                    if abs(f)+anorm!=anorm
                        g=w[i];h=_model_nurbs_root_pythag(f,g);w[i]=h;h=1.0/h
                        c=g*h;s=-f*h
                        for j in 1:m
                            y=a[j,nm];z=a[j,i]
                            a[j,nm]=y*c+z*s;a[j,i]=z*c-y*s
                        end
                    end
                end
            end
            z=w[k]
            if l==k
                if z<0.0
                    w[k]=-z
                    for j in 1:n
                        v[j,k]=-v[j,k]
                    end
                end
                break
            end
            iteration==30 && return nothing
            x=w[l];nm=k-1;y=w[nm];g=rv1[nm];h=rv1[k]
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=_model_nurbs_root_pythag(f,1.0)
            f=((x-z)*(x+z)+h*((y/(f+_model_nurbs_root_sign(g,f)))-h))/x
            c=1.0;s=1.0
            for j in l:nm
                i=j+1;g=rv1[i];y=w[i];h=s*g;g=c*g
                z=_model_nurbs_root_pythag(f,h);rv1[j]=z;c=f/z;s=h/z
                f=x*c+g*s;g=g*c-x*s;h=y*s;y=y*c
                for jj in 1:n
                    x=v[jj,j];z=v[jj,i]
                    v[jj,j]=x*c+z*s;v[jj,i]=z*c-x*s
                end
                z=_model_nurbs_root_pythag(f,h);w[j]=z
                if !iszero(z)
                    z=1.0/z;c=f*z;s=h*z
                end
                f=c*g+s*y;x=c*y-s*g
                for jj in 1:m
                    y=a[jj,j];z=a[jj,i]
                    a[jj,j]=y*c+z*s;a[jj,i]=z*c-y*s
                end
            end
            rv1[l]=0.0;rv1[k]=f;w[k]=x
        end
    end
    minimum=1.0e-6*maximum(w)
    for i in 1:n
        w[i]<minimum && (w[i]=0.0)
    end
    temp=zeros(n);result=zeros(n)
    for j in 1:n
        total=0.0
        if !iszero(w[j])
            for i in 1:m
                total+=a[i,j]*rhs[i]
            end
            total/=w[j]
        end
        temp[j]=total
    end
    for j in 1:n
        total=0.0
        for jj in 1:n
            total+=v[j,jj]*temp[jj]
        end
        result[j]=total
    end
    return result
end
