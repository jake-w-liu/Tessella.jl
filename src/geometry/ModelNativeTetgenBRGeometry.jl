# Floating-point geometry of tetgenBR.cxx, operation for operation: the
# dihedral cosines decide every quality flip and sliver split, so they must be
# the same bits as TetGen's (lu_decmp/lu_solve with scaled partial pivoting,
# facenormal with its shortest-edge pivot).

# lu_decmp(lu, 3, ps, &d, 0) on rows r1..r3; returns (ok, rows, ps).
function _lu3(r1::NTuple{3,Float64},r2::NTuple{3,Float64},r3::NTuple{3,Float64})
    rows=(r1,r2,r3)
    scales=(0.0,0.0,0.0)
    for i in 1:3
        biggest=0.0
        for j in 1:3
            tempf=abs(rows[i][j])
            biggest<tempf && (biggest=tempf)
        end
        biggest!=0.0 || return false,rows,(1,2,3)
        scales=Base.setindex(scales,1.0/biggest,i)
    end
    ps=(1,2,3)
    pivotindex=1
    for k in 1:2
        biggest=0.0
        for i in k:3
            tempf=abs(rows[ps[i]][k])*scales[ps[i]]
            if biggest<tempf
                biggest=tempf;pivotindex=i
            end
        end
        biggest==0.0 && return false,rows,ps
        if pivotindex!=k
            j=ps[k];ps=Base.setindex(ps,ps[pivotindex],k);ps=Base.setindex(ps,j,pivotindex)
        end
        pivot=rows[ps[k]][k]
        for i in k+1:3
            row=rows[ps[i]]
            mult=row[k]/pivot
            row=Base.setindex(row,mult,k)
            if mult!=0.0
                for j in k+1:3
                    row=Base.setindex(row,row[j]-mult*rows[ps[k]][j],j)
                end
            end
            rows=Base.setindex(rows,row,ps[i])
        end
    end
    return rows[ps[3]][3]!=0.0,rows,ps
end
# lu_solve(lu, 3, ps, b, 0) for b = the j-th unit vector.
function _lusolve3(rows,ps,b::NTuple{3,Float64})
    X=(0.0,0.0,0.0)
    for i in 1:3
        dot=0.0
        for j in 1:i-1
            dot+=rows[ps[i]][j]*X[j]
        end
        X=Base.setindex(X,b[ps[i]]-dot,i)
    end
    for i in 3:-1:1
        dot=0.0
        for j in i+1:3
            dot+=rows[ps[i]][j]*X[j]
        end
        X=Base.setindex(X,(X[i]-dot)/rows[ps[i]][i],i)
    end
    return X
end
@inline function _cross(v1,v2)
    return (v1[2]*v2[3]-v2[2]*v1[3],-(v1[1]*v2[3]-v2[1]*v1[3]),v1[1]*v2[2]-v2[1]*v1[2])
end
# facenormal(pa, pb, pc, n, pivot=1, NULL)
function facenormal(pa,pb,pc)
    v1=(pb[1]-pa[1],pb[2]-pa[2],pb[3]-pa[3])
    v2=(pa[1]-pc[1],pa[2]-pc[2],pa[3]-pc[3])
    v3=(pc[1]-pb[1],pc[2]-pb[2],pc[3]-pb[3])
    L1=dot3(v1,v1);L2=dot3(v2,v2);L3=dot3(v3,v3)
    pv1,pv2=if L1<L2
        L2<L3 ? (v1,v2) : (v3,v1)
    else
        L1<L3 ? (v1,v2) : (v2,v3)
    end
    n=_cross(pv1,pv2)
    return (-n[1],-n[2],-n[3])
end
@inline function _normalized(n)
    len=sqrt(dot3(n,n))
    return len,(len!=0.0 ? (n[1]/len,n[2]/len,n[3]/len) : n)
end
const _ALL180=(-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
# tetalldihedral(pa, pb, pc, pd, cosdd, &cosmaxd, NULL) -> (valid, cosdd, cosmaxd)
function tetalldihedral(pa,pb,pc,pd)
    A1=(pa[1]-pd[1],pa[2]-pd[2],pa[3]-pd[3])
    A2=(pb[1]-pd[1],pb[2]-pd[2],pb[3]-pd[3])
    A3=(pc[1]-pd[1],pc[2]-pd[2],pc[3]-pd[3])
    ok,rows,ps=_lu3(A1,A2,A3)
    vol=0.0
    N1=N2=N3=N4=(0.0,0.0,0.0)
    if ok
        vol=abs(rows[ps[1]][1]*rows[ps[2]][2]*rows[ps[3]][3])/6.0
        N1=_lusolve3(rows,ps,(1.0,0.0,0.0))
        N2=_lusolve3(rows,ps,(0.0,1.0,0.0))
        N3=_lusolve3(rows,ps,(0.0,0.0,1.0))
        N4=(-N1[1]-N2[1]-N3[1],-N1[2]-N2[2]-N3[2],-N1[3]-N2[3]-N3[3])
    end
    if vol>0
        l,N1=_normalized(N1);l==0.0 && (vol=0.0)
        l,N2=_normalized(N2);l==0.0 && (vol=0.0)
        l,N3=_normalized(N3);l==0.0 && (vol=0.0)
        l,N4=_normalized(N4);l==0.0 && (vol=0.0)
    end
    if vol<=0
        N1=facenormal(pc,pb,pd);N2=facenormal(pa,pc,pd)
        N3=facenormal(pb,pa,pd);N4=facenormal(pa,pb,pc)
        l,N1=_normalized(N1);l==0.0 && return false,_ALL180,-1.0
        l,N2=_normalized(N2);l==0.0 && return false,_ALL180,-1.0
        l,N3=_normalized(N3);l==0.0 && return false,_ALL180,-1.0
        l,N4=_normalized(N4);l==0.0 && return false,_ALL180,-1.0
    end
    N=(N1,N2,N3,N4)
    cosdd=_ALL180;cosmaxd=0.0
    for (i,(f1,f2)) in enumerate(((1,2),(2,3),(3,4),(1,4),(3,1),(2,4)))
        cosd=-dot3(N[f1],N[f2])
        cosd< -1.0 && (cosd=-1.0)
        cosd>1.0 && (cosd=1.0)
        cosdd=Base.setindex(cosdd,cosd,i)
        cosmaxd=i==1 ? cosd : (cosd<cosmaxd ? cosd : cosmaxd)
    end
    return true,cosdd,cosmaxd
end
@inline tetalldihedral(m::TG,a,b,c,d)=tetalldihedral(m.xyz[a],m.xyz[b],m.xyz[c],m.xyz[d])
