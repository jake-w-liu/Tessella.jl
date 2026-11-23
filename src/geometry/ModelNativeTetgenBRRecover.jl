# Boundary recovery by flips (recoversegments/recoveredgebyflips,
# recoversubfaces/recoverfacebyflips, removefacebyflips) and the
# triangle-edge intersection classifier they use. Recovery that needs
# Steiner points, and coplanar intersections (tri_edge_2d), throw
# `Unsupported`.

const DISJOINT=0;const SHAREVERT=2;const SHAREEDGE=3;const SHAREFACE=4
const TOUCHEDGE=5;const TOUCHFACE=6;const ACROSSVERT=7;const ACROSSEDGE=8;const ACROSSFACE=9

# tri_edge_tail(); returns (count, types, pos) with count 0, 1 or 2.
function tri_edge_tail(m,A,B,C,P,Q,R,sP,sQ,level)
    types=(0,0);pos=(0,0,0,0)
    if sP<0
        sQ<0 && return 0,types,pos
        if sQ>0
            U=(A,B,C);V=(P,Q,R);pu=(0,1,2);pv=(0,1,2);z1=0
        else
            U=(A,B,C);V=(P,Q,R);pu=(0,1,2);pv=(0,1,2);z1=1
        end
    elseif sP>0
        if sQ<0
            U=(A,B,C);V=(Q,P,R);pu=(0,1,2);pv=(1,0,2);z1=0
        elseif sQ>0
            return 0,types,pos
        else
            U=(B,A,C);V=(P,Q,R);pu=(1,0,2);pv=(0,1,2);z1=1
        end
    else
        if sQ<0
            U=(A,B,C);V=(Q,P,R);pu=(0,1,2);pv=(1,0,2);z1=1
        elseif sQ>0
            U=(B,A,C);V=(Q,P,R);pu=(1,0,2);pv=(1,0,2);z1=1
        else
            return tri_edge_2d(m,A,B,C,P,Q,R,level)
        end
    end
    s1=orient3d(m,U[1],U[2],V[1],V[2]);s1<0 && return 0,types,pos
    s2=orient3d(m,U[2],U[3],V[1],V[2]);s2<0 && return 0,types,pos
    s3=orient3d(m,U[3],U[1],V[1],V[2]);s3<0 && return 0,types,pos
    level==0 && return 1,types,pos
    t=0;p0=0;p1=0
    if z1==0
        if s1>0
            if s2>0
                s3>0 ? (t=ACROSSFACE;p0=3) : (t=ACROSSEDGE;p0=pu[3])
            else
                s3>0 ? (t=ACROSSEDGE;p0=pu[2]) : (t=ACROSSVERT;p0=pu[3])
            end
        else
            if s2>0
                s3>0 ? (t=ACROSSEDGE;p0=pu[1]) : (t=ACROSSVERT;p0=pu[1])
            elseif s3>0
                t=ACROSSVERT;p0=pu[2]
            end
        end
    else
        p1=pv[2]
        if s1>0
            if s2>0
                s3>0 ? (t=TOUCHFACE;p0=0) : (t=TOUCHEDGE;p0=pu[3])
            else
                s3>0 ? (t=TOUCHEDGE;p0=pu[2]) : (t=SHAREVERT;p0=pu[3])
            end
        else
            if s2>0
                s3>0 ? (t=TOUCHEDGE;p0=pu[1]) : (t=SHAREVERT;p0=pu[1])
            elseif s3>0
                t=SHAREVERT;p0=pu[2]
            end
        end
    end
    return 2,(t,DISJOINT),(p0,p1,0,0)
end
# tri_edge_2d(): [P,Q] and the triangle [A,B,C] are coplanar; R lifts the
# plane (`R == 0` is TetGen's NULL: a point above A along the normal). The
# slots a branch leaves unset stay DISJOINT/0 (TetGen leaves them as the
# caller's array had them), and the z1 == 2, s3 < 0, s2 > 0 branches
# overwrite types[0]/pos[0..1] as TetGen does.
function tri_edge_2d(m,A,B,C,P,Q,R,level)
    xyz(p)=p isa NTuple{3,Float64} ? p : m.xyz[p]
    o3(a,b,c,d)=Float64(orient3(xyz(a),xyz(b),xyz(c),xyz(d)))
    if R==0
        a=m.xyz[A];b=m.xyz[B];c=m.xyz[C]
        n=facenormal(a,b,c)
        len=sqrt(dot3(n,n))
        len==0 && return 0,(DISJOINT,DISJOINT),(0,0,0,0)
        inv_len=1.0/len
        n=(n[1]*inv_len,n[2]*inv_len,n[3]*inv_len)
        len=distance(a,b)
        len+=distance(b,c)
        len+=distance(c,a)
        len*=1.0/3.0
        R=(a[1]+len*n[1],a[2]+len*n[2],a[3]+len*n[3])
    elseif R==DUMMY
        throw(Unsupported("tri_edge_2d: the lift point is the dummy point"))
    end
    sA=o3(P,Q,R,A);sB=o3(P,Q,R,B);sC=o3(P,Q,R,C)
    I3=(0,1,2);ST=(2,0,1);STST=(1,2,0);I2=(0,1,2);SL=(1,0,2)
    sel(p)=p==(0,1,2) ? (A,B,C) : p==(2,0,1) ? (C,A,B) : (B,C,A)
    vsel(p)=p==(0,1,2) ? (P,Q,R) : (Q,P,R)
    pu,pv,z1=if sA<0
        if sB<0
            sC<0 && return 0,(DISJOINT,DISJOINT),(0,0,0,0)
            sC>0 ? (I3,I2,0) : (I3,I2,1)
        elseif sB>0
            sC<0 ? (ST,I2,0) : sC>0 ? (STST,SL,0) : (ST,I2,2)
        else
            sC<0 ? (ST,I2,1) : sC>0 ? (STST,SL,2) : (STST,SL,3)
        end
    elseif sA>0
        if sB<0
            sC<0 ? (STST,I2,0) : sC>0 ? (ST,SL,0) : (ST,SL,2)
        elseif sB>0
            sC>0 && return 0,(DISJOINT,DISJOINT),(0,0,0,0)
            sC<0 ? (I3,SL,0) : (I3,SL,1)
        else
            sC<0 ? (STST,I2,2) : sC>0 ? (ST,SL,1) : (STST,I2,3)
        end
    else
        if sB<0
            sC<0 ? (STST,I2,1) : sC>0 ? (I3,I2,2) : (ST,SL,3)
        elseif sB>0
            sC<0 ? (I3,SL,2) : sC>0 ? (STST,SL,1) : (ST,I2,3)
        else
            sC<0 ? (I3,SL,3) : sC>0 ? (I3,I2,3) : (I3,I2,4)
        end
    end
    U=sel(pu);V=vsel(pv)
    s1=o3(U[1],U[3],R,V[2])
    s2=o3(U[2],U[3],R,V[1])
    s1>0 && return 0,(DISJOINT,DISJOINT),(0,0,0,0)
    s2<0 && return 0,(DISJOINT,DISJOINT),(0,0,0,0)
    level==0 && return 1,(DISJOINT,DISJOINT),(0,0,0,0)
    types=[DISJOINT,DISJOINT];pos=[0,0,0,0]
    function set!(t0,p0,p1,t1=nothing,p2=0,p3=0)
        types[1]=t0;pos[1]=p0;pos[2]=p1
        t1===nothing || (types[2]=t1;pos[3]=p2;pos[4]=p3)
    end
    done()=(4,(types[1],types[2]),(pos[1],pos[2],pos[3],pos[4]))
    if z1==1
        if s1==0
            set!(SHAREVERT,pu[3],pv[2],DISJOINT)
        elseif s2==0
            set!(SHAREVERT,pu[3],pv[1],DISJOINT)
        else
            set!(ACROSSVERT,pu[3],pv[1],DISJOINT)
        end
        return done()
    end
    s3=o3(U[1],U[3],R,V[1])
    s4=o3(U[2],U[3],R,V[2])
    if z1==0
        if s1<0
            if s3>0
                if s4>0
                    set!(ACROSSEDGE,pu[3],pv[1],TOUCHFACE,3,pv[2])
                elseif s4==0
                    set!(ACROSSEDGE,pu[3],pv[1],TOUCHEDGE,pu[2],pv[2])
                else
                    set!(ACROSSEDGE,pu[3],pv[1],ACROSSEDGE,pu[2],pv[1])
                end
            elseif s3==0
                if s4>0
                    set!(TOUCHEDGE,pu[3],pv[1],TOUCHFACE,3,pv[2])
                elseif s4==0
                    set!(TOUCHEDGE,pu[3],pv[1],TOUCHEDGE,pu[2],pv[2])
                else
                    set!(TOUCHEDGE,pu[3],pv[1],ACROSSEDGE,pu[2],pv[1])
                end
            else
                if s2>0
                    if s4>0
                        set!(TOUCHFACE,3,pv[1],TOUCHFACE,3,pv[2])
                    elseif s4==0
                        set!(TOUCHFACE,3,pv[1],TOUCHEDGE,pu[2],pv[2])
                    else
                        set!(TOUCHFACE,3,pv[1],ACROSSEDGE,pu[2],pv[1])
                    end
                else
                    set!(TOUCHEDGE,pu[2],pv[1],DISJOINT)
                end
            end
        else
            set!(TOUCHEDGE,pu[3],pv[2],DISJOINT)
        end
    elseif z1==2
        if s1<0
            if s3>0
                if s4>0
                    set!(ACROSSVERT,pu[1],pv[1],TOUCHFACE,3,pv[2])
                elseif s4==0
                    set!(ACROSSVERT,pu[1],pv[1],TOUCHEDGE,pu[2],pv[2])
                else
                    set!(ACROSSVERT,pu[1],pv[1],ACROSSEDGE,pu[2],pv[1])
                end
            elseif s3==0
                if s4>0
                    set!(SHAREVERT,pu[1],pv[1],TOUCHFACE,3,pv[2])
                elseif s4==0
                    set!(SHAREVERT,pu[1],pv[1],TOUCHEDGE,pu[2],pv[2])
                else
                    set!(SHAREVERT,pu[1],pv[1],ACROSSEDGE,pu[2],pv[1])
                end
            else
                if s2>0
                    # TetGen writes the second record over the first here.
                    if s4>0
                        set!(TOUCHFACE,3,pv[2])
                    elseif s4==0
                        set!(TOUCHEDGE,pu[2],pv[2])
                    else
                        set!(ACROSSEDGE,pu[2],pv[1])
                    end
                else
                    set!(TOUCHEDGE,pu[2],pv[1],DISJOINT)
                end
            end
        else
            set!(SHAREVERT,pu[1],pv[2],DISJOINT)
        end
    elseif z1==3
        if s1<0
            if s3>0
                if s4>0
                    set!(ACROSSVERT,pu[1],pv[1],TOUCHEDGE,pu[1],pv[2])
                elseif s4==0
                    set!(ACROSSVERT,pu[1],pv[1],SHAREVERT,pu[2],pv[2])
                else
                    set!(ACROSSVERT,pu[1],pv[1],ACROSSVERT,pu[2],pv[1])
                end
            elseif s3==0
                if s4>0
                    set!(SHAREVERT,pu[1],pv[1],TOUCHEDGE,pu[1],pv[2])
                elseif s4==0
                    set!(SHAREEDGE,pu[1],pv[1],DISJOINT)
                else
                    set!(SHAREVERT,pu[1],pv[1],ACROSSVERT,pu[2],pv[1])
                end
            else
                if s2>0
                    if s4>0
                        set!(TOUCHEDGE,pu[1],pv[1],TOUCHEDGE,pu[1],pv[2])
                    elseif s4==0
                        set!(TOUCHEDGE,pu[1],pv[1],SHAREVERT,pu[2],pv[2])
                    else
                        set!(TOUCHEDGE,pu[1],pv[1],ACROSSVERT,pu[2],pv[1])
                    end
                else
                    set!(SHAREVERT,pu[2],pv[1],DISJOINT)
                end
            end
        else
            set!(SHAREVERT,pu[1],pv[2],DISJOINT)
        end
    end
    return done()
end

function tri_edge_test(m,A,B,C,P,Q,R,level)
    sP=orient3d(m,A,B,C,P);sQ=orient3d(m,A,B,C,Q)
    return tri_edge_tail(m,A,B,C,P,Q,R,sP,sQ,level)
end

# The segment/facet constraint branches of checkflipeligibility().
function _checkconstraints(m,fliptype,pa,pb,pc,pd,pe,fc)
    rejflag=false
    if fc.seg[1]!=0
        if fliptype==1
            for p in (pa,pb,pc)
                rejflag && break
                p==DUMMY && continue
                intflag,types,poss=tri_edge_test(m,pe,pd,p,fc.seg[1],fc.seg[2],0,1)
                if intflag==2
                    types[1]==ACROSSFACE && (rejflag=true)
                    types[1]==ACROSSEDGE && poss[1]==0 && (rejflag=true)
                elseif intflag==4
                    types[1]==ACROSSEDGE && poss[1]==0 && (rejflag=true)
                end
            end
        elseif fliptype==2 && pc!=DUMMY
            intflag,types,poss=tri_edge_test(m,pa,pb,pc,fc.seg[1],fc.seg[2],0,1)
            if intflag==2
                types[1]==ACROSSFACE && (rejflag=true)
            elseif intflag==4
                types[1]==ACROSSEDGE && (rejflag=true)
            end
        end
    end
    if fc.fac[1]!=0 && !rejflag && fliptype==1
        intflag,types,poss=tri_edge_test(m,fc.fac[1],fc.fac[2],fc.fac[3],pe,pd,0,1)
        if intflag==2
            (types[1]==ACROSSFACE || types[1]==ACROSSEDGE) && (rejflag=true)
        elseif intflag==4
            for i in 1:2
                (types[i]==ACROSSFACE || types[i]==ACROSSEDGE) && (rejflag=true;break)
            end
        end
    end
    return rejflag
end

function removefacebyflips!(m,flipface,fc)
    F=m.F3
    F[1]=flipface;F[2]=fsym(m,flipface)
    pa=org(m,F[1]);pb=dest(m,F[1]);pc=apex(m,F[1]);pd=oppo(m,F[1]);pe=oppo(m,F[2])
    flipedge=flipface;reducflag=false
    if orient3d(m,pa,pb,pd,pe)>0
        if orient3d(m,pb,pc,pd,pe)>0
            if orient3d(m,pc,pa,pd,pe)>0
                reducflag=true
            else
                flipedge=eprev(flipface)
            end
        else
            flipedge=enext(flipface)
        end
    end
    if reducflag
        flip23!(m,F,0,0,fc)
        return true
    end
    nn,_=removeedgebyflips!(m,flipedge,fc)
    return nn==2
end

function recoveredgebyflips!(m,startpt,endpt,fullsearch,fc=Flip())
    _flipreset!(fc)
    fc.seg=(startpt,endpt)
    fc.checkflipeligibility=true
    while true
        searchtet,dir=finddirection!(m,point2tetorg(m,startpt),endpt)
        if dir===:acrossvert
            dest(m,searchtet)==endpt && return true,searchtet
            throw(Unsupported("recoveredgebyflips: a vertex lies on a segment"))
        end
        searchtet=enextesym(searchtet)
        if dir===:acrossface
            m.checksubfaceflag && issubface(m,searchtet) &&
                throw(Unsupported("recoveredgebyflips: a segment crosses a facet"))
            removefacebyflips!(m,searchtet,fc) && continue
        elseif dir===:acrossedge
            m.checksubsegflag && issubseg(m,searchtet) &&
                throw(Unsupported("recoveredgebyflips: two segments intersect"))
            nn,_=removeedgebyflips!(m,searchtet,fc)
            nn==2 && continue
        end
        if fullsearch
            searchtet,dir=finddirection!(m,point2tetorg(m,startpt),endpt)
            searchtet=enextesym(searchtet)
            success=false
            neightet=TNULL;pos=0
            while true
                searchtet=fsym(m,searchtet)
                if dir===:acrossface || dir==ACROSSFACE
                    neightet=searchtet
                    j=neightet.ver&3
                    d=DISJOINT
                    for i in j+1:j+3
                        neightet=Tri(neightet.tet,i%4)
                        intflag,types,poss=tri_edge_test(m,org(m,neightet),dest(m,neightet),
                            apex(m,neightet),startpt,endpt,oppo(m,neightet),1)
                        if intflag!=0
                            d=types[1];pos=poss[1];break
                        else
                            d=DISJOINT;pos=0
                        end
                    end
                    d==DISJOINT && throw(Unsupported("recoveredgebyflips: lost the segment"))
                    dir=d
                elseif dir===:acrossedge || dir==ACROSSEDGE
                    d=DISJOINT
                    while true
                        for i in 1:2
                            neightet=i==1 ? enextesym(searchtet) : eprevesym(searchtet)
                            intflag,types,poss=tri_edge_test(m,org(m,neightet),dest(m,neightet),
                                apex(m,neightet),startpt,endpt,oppo(m,neightet),1)
                            if intflag!=0
                                d=types[1];pos=poss[1];break
                            else
                                d=DISJOINT;pos=0
                            end
                        end
                        d!=DISJOINT && break
                        searchtet=fnext(m,searchtet)
                    end
                    dir=d
                else
                    throw(Unsupported("recoveredgebyflips: unexpected direction"))
                end
                for _ in 1:pos;neightet=enext(neightet);end
                if dir==SHAREVERT
                    org(m,neightet)==endpt && break
                    throw(Unsupported("recoveredgebyflips: unexpected shared vertex"))
                end
                searchtet=neightet
                bf=(org(m,searchtet),dest(m,searchtet),apex(m,searchtet),oppo(m,searchtet))
                if dir==ACROSSFACE
                    m.checksubfaceflag && issubface(m,searchtet) &&
                        throw(Unsupported("recoveredgebyflips: a segment crosses a facet"))
                    if removefacebyflips!(m,searchtet,fc)
                        success=true;break
                    end
                elseif dir==ACROSSEDGE
                    m.checksubsegflag && issubseg(m,searchtet) &&
                        throw(Unsupported("recoveredgebyflips: two segments intersect"))
                    nn,searchtet=removeedgebyflips!(m,searchtet,fc)
                    if nn==2
                        success=true;break
                    end
                else
                    throw(Unsupported("recoveredgebyflips: a vertex lies on a segment"))
                end
                if searchtet.tet==0 || isdeadtet(m,searchtet) || org(m,searchtet)!=bf[1] ||
                        dest(m,searchtet)!=bf[2] || apex(m,searchtet)!=bf[3] ||
                        oppo(m,searchtet)!=bf[4]
                    t,dir1=finddirection!(m,point2tetorg(m,bf[1]),bf[2])
                    searchtet=TNULL
                    if dir1===:acrossvert && dest(m,t)==bf[2]
                        spintet=t
                        while true
                            if apex(m,spintet)==bf[3]
                                searchtet=spintet;break
                            end
                            spintet=fnext(m,spintet)
                            spintet.tet==t.tet && break
                        end
                        if searchtet.tet!=0 && oppo(m,searchtet)!=bf[4]
                            searchtet=fsym(m,searchtet)
                            oppo(m,searchtet)!=bf[4] && (searchtet=TNULL)
                        end
                    end
                    searchtet.tet==0 && break
                end
            end
            success && continue
        end
        break
    end
    return false,TNULL
end

function _bondsegment!(m,sseg,t)
    sstbond1!(m,sseg,t)
    spintet=t
    while true
        tssbond1!(m,spintet,sseg)
        spintet=fnext(m,spintet)
        spintet.tet==t.tet && break
    end
end

function recoversegments!(m,misseglist,fullsearch)
    S=m.subsegs
    fc=Flip()
    while !isempty(m.subsegstack)
        sseg=pop!(m.subsegstack)
        sstpivot1(m,sseg).tet!=0 && continue
        startpt=sorg(S,sseg);endpt=sdest(S,sseg)
        ok,t=recoveredgebyflips!(m,startpt,endpt,false,fc)
        ok || ((ok,t)=recoveredgebyflips!(m,endpt,startpt,false,fc))
        !ok && fullsearch && ((ok,t)=recoveredgebyflips!(m,startpt,endpt,true,fc))
        if ok
            _bondsegment!(m,sseg,t)
        else
            push!(misseglist,sseg)
        end
    end
end

function recoverfacebyflips!(m,pa,pb,pc,fc=Flip())
    _flipreset!(fc)
    fc.fac=(pa,pb,pc)
    fc.checkflipeligibility=true
    fac=(pa,pb,pc)
    for i in 0:2
        while true
            searchtet,_=finddirection!(m,point2tetorg(m,fac[i+1]),fac[mod(i+1,3)+1])
            spintet=searchtet
            while true
                if apex(m,spintet)==fac[mod(i+2,3)+1]
                    searchtet=spintet
                    for _ in i:-1:1;searchtet=eprev(searchtet);end
                    return true,searchtet
                end
                spintet=fnext(m,spintet)
                spintet.tet==searchtet.tet && break
            end
            flipedge=TNULL
            spintet=searchtet
            while true
                pd=apex(m,spintet);pe=oppo(m,spintet)
                if pd!=DUMMY && pe!=DUMMY
                    intflag,types,poss=tri_edge_test(m,pa,pb,pc,pd,pe,0,1)
                    if intflag>0
                        intflag==2 || throw(Unsupported("recoverfacebyflips: a coplanar edge"))
                        flipedge=edestoppo(spintet)
                        if types[1]==ACROSSFACE || types[1]==ACROSSEDGE
                            issubseg(m,flipedge) &&
                                throw(Unsupported("recoverfacebyflips: a segment crosses a facet"))
                            chkface=flipedge
                            while true
                                issubface(m,chkface) && break
                                chkface=fsym(m,chkface)
                                chkface.tet==flipedge.tet && break
                            end
                            issubface(m,chkface) &&
                                throw(Unsupported("recoverfacebyflips: two facets intersect"))
                        else
                            throw(Unsupported("recoverfacebyflips: a vertex touches a facet"))
                        end
                        break
                    end
                end
                spintet=fnext(m,spintet)
                spintet.tet==searchtet.tet &&
                    throw(Unsupported("recoverfacebyflips: no crossing edge"))
            end
            nn,_=removeedgebyflips!(m,flipedge,fc)
            nn==2 && continue
            break
        end
    end
    return false,TNULL
end

function recoversubfaces!(m,misshlist)
    F=m.subfaces
    fc=Flip()
    while !isempty(m.subfacstack)
        searchsh=pop!(m.subfacstack)
        shdead(F,searchsh) && continue
        stpivot(m,searchsh).tet!=0 && continue
        for _ in 1:3
            seg=sspivot(F,searchsh)
            seg.sh==0 && throw(Unsupported("a subface edge without a segment"))
            sstpivot1(m,seg).tet==0 && throw(Unsupported("an unrecovered segment of a subface"))
            searchsh=senext(searchsh)
        end
        ok,searchtet=recoverfacebyflips!(m,sorg(F,searchsh),sdest(F,searchsh),sapex(F,searchsh),fc)
        if ok
            tsbond!(m,searchtet,searchsh)
            tsbond!(m,fsym(m,searchtet),sesym(searchsh))
        else
            push!(misshlist,searchsh)
        end
    end
end
