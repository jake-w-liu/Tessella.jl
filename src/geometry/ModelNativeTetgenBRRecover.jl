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
            throw(Unsupported("tri_edge_2d: a coplanar triangle-edge intersection"))
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
    F=Vector{Tri}(undef,3)
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

function recoveredgebyflips!(m,startpt,endpt,fullsearch)
    fc=Flip()
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
    while !isempty(m.subsegstack)
        sseg=pop!(m.subsegstack)
        sstpivot1(m,sseg).tet!=0 && continue
        startpt=sorg(S,sseg);endpt=sdest(S,sseg)
        ok,t=recoveredgebyflips!(m,startpt,endpt,false)
        ok || ((ok,t)=recoveredgebyflips!(m,endpt,startpt,false))
        !ok && fullsearch && ((ok,t)=recoveredgebyflips!(m,startpt,endpt,true))
        if ok
            _bondsegment!(m,sseg,t)
        else
            push!(misseglist,sseg)
        end
    end
end

function recoverfacebyflips!(m,pa,pb,pc)
    fc=Flip()
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
        ok,searchtet=recoverfacebyflips!(m,sorg(F,searchsh),sdest(F,searchsh),sapex(F,searchsh))
        if ok
            tsbond!(m,searchtet,searchsh)
            tsbond!(m,fsym(m,searchtet),sesym(searchsh))
        else
            push!(misshlist,searchsh)
        end
    end
end
