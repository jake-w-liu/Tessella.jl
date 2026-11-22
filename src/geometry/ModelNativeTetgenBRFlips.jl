# Elementary and composite flips of tetgenBR.cxx (flip23, flip32, flipnm,
# flipnm_post, checkflipeligibility, removeedgebyflips, lawsonflip3d) and the
# walks that find edges (finddirection, getedge, gettetrahedron). A flip takes
# its three tetrahedra as F[o+1], F[o+2], F[o+3] and returns its results there,
# like TetGen's `triface *fliptets`.

@inline function enqueuetetrahedron!(m,t)
    if !marktest2ed(m,t)
        marktest2!(m,t)
        push!(m.badtetrahedrons,t)
    end
end
@inline function flippush!(m,t)
    if !facemarked(m,t)
        markface!(m,t)
        push!(m.flipstack,t)
    end
end

function flip23!(m,F,o,hullflag,fc)
    dummyflag=0
    if hullflag>0
        if oppo(m,F[o+2])==DUMMY
            F[o+1],F[o+2]=F[o+2],F[o+1]
            dummyflag=-1
        elseif org(m,F[o+1])==DUMMY
            dummyflag=1
            F[o+1]=enext(F[o+1]);F[o+2]=eprev(F[o+2])
        elseif dest(m,F[o+1])==DUMMY
            dummyflag=2
            F[o+1]=eprev(F[o+1]);F[o+2]=enext(F[o+2])
        end
    end
    f0=F[o+1];f1=F[o+2]
    pa=org(m,f0);pb=dest(m,f0);pc=apex(m,f0);pd=oppo(m,f0);pe=oppo(m,f1)
    t1=fnext(m,f0);f0=enext(f0);t2=fnext(m,f0);f0=enext(f0);t3=fnext(m,f0);f0=enext(f0)
    top=(t1,t2,t3)
    b1=fnext(m,f1);f1=eprev(f1);b2=fnext(m,f1);f1=eprev(f1);b3=fnext(m,f1);f1=eprev(f1)
    bot=(b1,b2,b3)
    f0=Tri(f0.tet,11);f1=Tri(f1.tet,11)
    m.tmark[f0.tet]=0;m.tmark[f1.tet]=0
    m.checksubsegflag && (for i in 1:6;m.tseg[6*(f0.tet-1)+i]=0;m.tseg[6*(f1.tet-1)+i]=0;end)
    m.checksubfaceflag && (for i in 1:4;m.tsub[4*(f0.tet-1)+i]=0;m.tsub[4*(f1.tet-1)+i]=0;end)
    f2=maketetrahedron!(m)
    if hullflag>0
        if pd!=DUMMY
            setvertices!(m,f0,pe,pd,pa,pb)
            setvertices!(m,f1,pe,pd,pb,pc)
            if pc!=DUMMY
                setvertices!(m,f2,pe,pd,pc,pa)
            else
                setvertices!(m,f2,pd,pe,pa,pc);f2=esym(f2)
            end
        else
            setvertices!(m,f0,pa,pb,pe,pd)
            setvertices!(m,f1,pb,pc,pe,pd)
            setvertices!(m,f2,pc,pa,pe,pd)
            f0=enext(eprevesym(f0));f1=enext(eprevesym(f1));f2=enext(eprevesym(f2))
        end
    else
        setvertices!(m,f0,pe,pd,pa,pb)
        setvertices!(m,f1,pe,pd,pb,pc)
        setvertices!(m,f2,pe,pd,pc,pa)
    end
    new=(f0,f1,f2)
    for i in 1:3
        bond!(m,esym(new[i]),new[mod(i,3)+1])
    end
    for i in 1:3
        bond!(m,eorgoppo(new[i]),top[i])
    end
    for i in 1:3
        bond!(m,edestoppo(new[i]),bot[i])
    end
    if m.checksubsegflag
        for i in 1:3
            if issubseg(m,top[i])
                seg=tsspivot1(m,top[i]);nf=eorgoppo(new[i])
                tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
            end
        end
        for i in 1:3
            cas=eprev(top[i])
            if issubseg(m,cas)
                seg=tsspivot1(m,cas)
                nf=enext(new[i]);tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
                nf=eprev(esym(new[mod(i+1,3)+1]));tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
            end
        end
        for i in 1:3
            cas=enext(bot[i])
            if issubseg(m,cas)
                seg=tsspivot1(m,cas)
                nf=eprev(new[i]);tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
                nf=enext(esym(new[mod(i+1,3)+1]));tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
            end
        end
    end
    if m.checksubfaceflag
        for i in 1:3
            if issubface(m,top[i])
                sh=tspivot(m,top[i]);tsbond!(m,eorgoppo(new[i]),sesym(sh))
            end
        end
        for i in 1:3
            if issubface(m,bot[i])
                sh=tspivot(m,bot[i]);tsbond!(m,edestoppo(new[i]),sesym(sh))
            end
        end
    end
    if (fc.chkencflag&4)!=0
        for i in 1:3;enqueuetetrahedron!(m,new[i]);end
    end
    setpoint2tet!(m,pa,Tri(f0.tet,0));setpoint2tet!(m,pb,Tri(f0.tet,0))
    setpoint2tet!(m,pc,Tri(f1.tet,0));setpoint2tet!(m,pd,Tri(f0.tet,0))
    setpoint2tet!(m,pe,Tri(f0.tet,0))
    if hullflag>0 && dummyflag!=0
        if dummyflag==-1
            f0=esym(f0);f1=esym(f1);f2=esym(f2)
            f1,f2=f2,f1
        elseif dummyflag==1
            f0,f1,f2=f2,f0,f1
        else
            f0,f1,f2=f1,f2,f0
        end
    end
    F[o+1]=f0;F[o+2]=f1;F[o+3]=f2
    if fc.enqflag>0
        for i in 1:3;flippush!(m,eprevesym(F[o+i]));end
        if fc.enqflag>1
            for i in 1:3;flippush!(m,enextesym(F[o+i]));end
        end
    end
    return nothing
end

function flip32!(m,F,o,hullflag,fc)
    dummyflag=0
    f0=F[o+1];f1=F[o+2];f2=F[o+3]
    if hullflag>0
        if org(m,f0)==DUMMY
            f0=esym(f0);f1=esym(f1);f2=esym(f2)
            f1,f2=f2,f1
            dummyflag=-1
        elseif apex(m,f0)==DUMMY
            dummyflag=1
            f0,f1,f2=f1,f2,f0
        elseif apex(m,f1)==DUMMY
            dummyflag=2
            f0,f1,f2=f2,f0,f1
        end
    end
    pa=apex(m,f0);pb=apex(m,f1);pc=apex(m,f2);pd=dest(m,f0);pe=org(m,f0)
    old=(f0,f1,f2)
    top=(fsym(m,eorgoppo(f0)),fsym(m,eorgoppo(f1)),fsym(m,eorgoppo(f2)))
    bot=(fsym(m,edestoppo(f0)),fsym(m,edestoppo(f1)),fsym(m,edestoppo(f2)))
    if m.checksubfaceflag
        for t in old
            issubface(m,t) && throw(Unsupported(
                "flip32 across interior subfaces requires the surface 2-to-2 flip"))
        end
    end
    f0=Tri(f0.tet,11);f1=Tri(f1.tet,11)
    m.tmark[f0.tet]=0;m.tmark[f1.tet]=0
    m.checksubsegflag && (for i in 1:6;m.tseg[6*(f0.tet-1)+i]=0;m.tseg[6*(f1.tet-1)+i]=0;end)
    m.checksubfaceflag && (for i in 1:4;m.tsub[4*(f0.tet-1)+i]=0;m.tsub[4*(f1.tet-1)+i]=0;end)
    tetrahedrondealloc!(m,f2.tet)
    if hullflag>0 && pc==DUMMY
        setvertices!(m,f0,pb,pa,pd,pc)
        setvertices!(m,f1,pa,pb,pe,pc)
        f0=esym(f0);f1=esym(f1)
    else
        setvertices!(m,f0,pa,pb,pc,pd)
        setvertices!(m,f1,pb,pa,pc,pe)
    end
    bond!(m,f0,f1)
    for i in 1:3
        bond!(m,esym(f0),top[i]);f0=enext(f0)
    end
    for i in 1:3
        bond!(m,esym(f1),bot[i]);f1=eprev(f1)
    end
    if m.checksubsegflag
        for i in 1:3
            if issubseg(m,top[i])
                seg=tsspivot1(m,top[i])
                tssbond1!(m,f0,seg);sstbond1!(m,seg,f0)
                tssbond1!(m,f1,seg);sstbond1!(m,seg,f1)
            end
            f0=enext(f0);f1=eprev(f1)
        end
        for i in 1:3
            nf=eprev(esym(f0));cas=enext(top[i])
            if issubseg(m,cas)
                seg=tsspivot1(m,cas);tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
            end
            f0=enext(f0)
        end
        for i in 1:3
            nf=enext(esym(f1));cas=eprev(bot[i])
            if issubseg(m,cas)
                seg=tsspivot1(m,cas);tssbond1!(m,nf,seg);sstbond1!(m,seg,nf)
            end
            f1=eprev(f1)
        end
    end
    if m.checksubfaceflag
        for i in 1:3
            if issubface(m,top[i])
                sh=tspivot(m,top[i]);tsbond!(m,esym(f0),sesym(sh))
            end
            f0=enext(f0)
        end
        for i in 1:3
            if issubface(m,bot[i])
                sh=tspivot(m,bot[i]);tsbond!(m,esym(f1),sesym(sh))
            end
            f1=eprev(f1)
        end
    end
    if (fc.chkencflag&4)!=0
        enqueuetetrahedron!(m,f0);enqueuetetrahedron!(m,f1)
    end
    setpoint2tet!(m,pa,Tri(f0.tet,0));setpoint2tet!(m,pb,Tri(f0.tet,0))
    setpoint2tet!(m,pc,Tri(f0.tet,0));setpoint2tet!(m,pd,Tri(f0.tet,0))
    setpoint2tet!(m,pe,Tri(f1.tet,0))
    if hullflag>0 && dummyflag!=0
        if dummyflag==-1
            f0,f1=f1,f0
        elseif dummyflag==1
            f0=eprev(f0);f1=enext(f1)
        else
            f0=enext(f0);f1=eprev(f1)
        end
    end
    F[o+1]=f0;F[o+2]=f1
    if fc.enqflag>0
        flippush!(m,enextesym(f0));flippush!(m,eprevesym(f1))
        if fc.enqflag>1
            flippush!(m,eprevesym(f0));flippush!(m,enextesym(f1))
            flippush!(m,esym(f0));flippush!(m,esym(f1))
        end
    end
    return nothing
end

function checkflipeligibility(m,fliptype,pa,pb,pc,pd,pe,level,edgepivot,fc)
    rejflag=(fc.seg[1]!=0 || fc.fac[1]!=0) && _checkconstraints(m,fliptype,pa,pb,pc,pd,pe,fc)
    (fc.remove_large_angle && !rejflag) || return rejflag
    if fliptype==1
        if pc!=DUMMY && pe!=DUMMY && pd!=DUMMY
            rejflag=_rejectangle!(m,fc,pe,pd,pb,pc) || _rejectangle!(m,fc,pe,pd,pc,pa)
        end
    elseif fliptype==2
        if pa!=DUMMY && pb!=DUMMY && pc!=DUMMY
            if level==0
                rejflag=_rejectangle!(m,fc,pa,pb,pc,pd) || _rejectangle!(m,fc,pb,pa,pc,pe)
            elseif edgepivot==1
                rejflag=_rejectangle!(m,fc,pb,pa,pc,pe)
            else
                rejflag=_rejectangle!(m,fc,pa,pb,pc,pd)
            end
        end
    end
    return rejflag
end
# A new tetrahedron must not have a dihedral angle as large as the one being
# removed; the largest new angle is recorded in fc.cosdihed_out.
function _rejectangle!(m,fc,p,q,r,s)
    _,_,cosmaxd=tetalldihedral(m,p,q,r,s)
    diff=cosmaxd-fc.cosdihed_in
    abs(diff/fc.cosdihed_in)<m.epsilon && (diff=0.0)
    diff<=0 && return true
    cosmaxd<fc.cosdihed_out && (fc.cosdihed_out=cosmaxd)
    return false
end

# finddirection(): walk from org(t) toward `endpt`; returns (t, kind) with
# kind in (:acrossvert, :acrossedge, :acrossface).
function finddirection!(m,t,endpt)
    pa=org(m,t)
    if tvert(m,t.tet,7)==DUMMY
        n=rawnb(m,t.tet,3)
        s=n.tet
        t=tvert(m,s,4)==pa ? Tri(s,11) : tvert(m,s,5)==pa ? Tri(s,3) :
          tvert(m,s,6)==pa ? Tri(s,7) : Tri(s,0)
    end
    pb=dest(m,t)
    pb==endpt && return t,:acrossvert
    pc=apex(m,t)
    if pc==endpt
        return eprevesym(t),:acrossvert
    end
    for _ in 1:1_000_000
        pd=oppo(m,t)
        if pd==endpt
            return enext(esym(t)),:acrossvert
        end
        if pd==DUMMY
            m.nonconvex && return t,:acrossface
            throw(Unsupported("finddirection reached the hull of a convex mesh"))
        end
        hori=orient3d(m,pa,pb,pc,endpt)
        rori=orient3d(m,pb,pa,pd,endpt)
        lori=orient3d(m,pa,pc,pd,endpt)
        move=if hori>0
            if rori>0
                if lori>0
                    s=randomnation!(m,3)
                    s==0 ? 0 : s==1 ? 1 : 2
                else
                    randomnation!(m,2)!=0 ? 0 : 1
                end
            else
                lori>0 ? (randomnation!(m,2)!=0 ? 0 : 2) : 0
            end
        else
            if rori>0
                lori>0 ? (randomnation!(m,2)!=0 ? 1 : 2) : 1
            elseif lori>0
                2
            else
                if hori==0
                    rori==0 && return t,:acrossvert
                    lori==0 && return eprevesym(t),:acrossvert
                    return t,:acrossedge
                end
                if rori==0
                    lori==0 && return enext(esym(t)),:acrossvert
                    return enext(esym(t)),:acrossedge
                end
                lori==0 && return eprevesym(t),:acrossedge
                return t,:acrossface
            end
        end
        if move==1
            t=fnext(m,t)
        elseif move==2
            t=enext(fnext(m,eprev(t)))
        else
            t=enext(fsym(m,t))
        end
        pb=dest(m,t);pc=apex(m,t)
    end
    throw(Unsupported("finddirection did not terminate"))
end

function getedge!(m,e1,e2,tedge)
    if !isdeadtet(m,tedge)
        if org(m,tedge)==e1
            dest(m,tedge)==e2 && return true,tedge
        elseif org(m,tedge)==e2
            dest(m,tedge)==e1 && return true,esym(tedge)
        end
    end
    tedge,_=finddirection!(m,point2tetorg(m,e1),e2)
    dest(m,tedge)==e2 && return true,tedge
    tedge,_=finddirection!(m,point2tetorg(m,e2),e1)
    dest(m,tedge)==e1 && return true,esym(tedge)
    searchtet=enextesym(point2tetorg(m,e1))
    for _ in 1:3
        apex(m,searchtet)==e2 && return true,eorgoppo(searchtet)
        searchtet=enext(searchtet)
    end
    neightet=esym(fnext(m,searchtet))
    apex(m,neightet)==e2 && return true,eorgoppo(neightet)
    list=m.cavebdrylist
    infect!(m,searchtet);push!(list,searchtet)
    infect!(m,neightet);push!(list,neightet)
    done=false
    i=1
    while i<=length(list) && !done
        searchtet=list[i]
        for _ in 1:2
            done && break
            searchtet=enext(searchtet)
            neightet=fnext(m,searchtet)
            if !infected(m,neightet)
                neightet=esym(neightet)
                if apex(m,neightet)==e2
                    tedge=eorgoppo(neightet);done=true
                else
                    infect!(m,neightet);push!(list,neightet)
                end
            end
        end
        i+=1
    end
    for t in list;uninfect!(m,t);end
    empty!(list)
    return done,tedge
end

function gettetrahedron!(m,pa,pb,pc,pd,searchtet)
    ok,searchtet=getedge!(m,pa,pb,searchtet)
    if ok
        spintet=searchtet
        while true
            if apex(m,spintet)==pc
                searchtet=spintet;break
            end
            spintet=fnext(m,spintet)
            spintet.tet==searchtet.tet && break
        end
        if apex(m,searchtet)==pc
            oppo(m,searchtet)==pd && return true,searchtet
            searchtet=fsym(m,searchtet)
            oppo(m,searchtet)==pd && return true,searchtet
        end
    end
    return false,searchtet
end

function flipnm!(m,abtets,n,level,abedgepivot,fc)
    pa=org(m,abtets[1]);pb=dest(m,abtets[1])
    fliptets=m.F3
    if n>3
        reflexlinkedgecount=0
        i=0
        while i<n
            cur=abtets[i+1];prev=abtets[mod(i-1,n)+1]
            if m.checksubfaceflag && issubface(m,cur)
                i+=1;continue
            end
            if elemcounter(m,cur)>1 || elemcounter(m,prev)>1
                i+=1;continue
            end
            pc=apex(m,cur);pd=apex(m,abtets[mod(i+1,n)+1]);pe=apex(m,prev)
            if pd==DUMMY || pe==DUMMY
                i+=1;continue
            end
            reducflag=false;ori=0.0
            hullflag=pc==DUMMY ? 1 : 0
            hulledgeflag=false
            if hullflag==0
                ori=orient3d(m,pb,pc,pd,pe)
                if ori>0
                    ori=orient3d(m,pc,pa,pd,pe)
                    if ori>0
                        ori=orient3d(m,pa,pb,pd,pe)
                        if ori>0
                            reducflag=true
                        elseif ori==0 && n==4
                            reducflag=true
                            pf=apex(m,abtets[mod(i+2,n)+1]);hulledgeflag=pf==DUMMY
                        end
                    end
                end
                reducflag || (reflexlinkedgecount+=1)
            else
                if n==4
                    pf=apex(m,abtets[mod(i+2,n)+1])
                    ori=orient3d(m,pd,pe,pf,pa)
                    if ori<0
                        ori=orient3d(m,pe,pd,pf,pb)
                        if ori<0
                            reducflag=true;ori=0.0;hulledgeflag=true
                        end
                    end
                end
            end
            if reducflag && m.nonconvex && hulledgeflag
                found,_=getedge!(m,pe,pd,TNULL)
                found && (reducflag=false)
            end
            if reducflag
                rejflag=fc.checkflipeligibility &&
                    checkflipeligibility(m,1,pa,pb,pc,pd,pe,level,abedgepivot,fc)
                if !rejflag
                    fliptets[1]=cur;fliptets[2]=fsym(m,cur)
                    flip23!(m,fliptets,0,hullflag,fc)
                    fliptets[1]=edestoppo(fliptets[1])
                    increaseelemcounter!(m,fliptets[1])
                    abtets[mod(i-1,n)+1]=fliptets[1]
                    for j in i:n-2
                        abtets[j+1]=abtets[j+2]
                    end
                    abtets[n]=Tri(Int(pc),(1<<4)|(i<<6))
                    if fc.collectnewtets
                        push!(m.cavetetlist,fliptets[2]);push!(m.cavetetlist,fliptets[3])
                    end
                    nn=flipnm!(m,abtets,n-1,level,abedgepivot,fc)
                    nn==2 && return nn
                    if fc.unflip || ori==0
                        fliptets[1]=edestoppo(abtets[mod(i-1,n-1)+1])
                        fliptets[2]=fnext(m,fliptets[1])
                        fliptets[3]=fnext(m,fliptets[2])
                        flip32!(m,fliptets,0,hullflag,fc)
                        increaseelemcounter!(m,fliptets[1]);increaseelemcounter!(m,fliptets[2])
                        for j in n-2:-1:i
                            abtets[j+2]=abtets[j+1]
                        end
                        abtets[mod(i-1,n)+1]=esym(fliptets[2])
                        abtets[i+1]=fliptets[1]
                        nn+=1
                        fc.collectnewtets && resize!(m.cavetetlist,length(m.cavetetlist)-2)
                    end
                    fc.unflip || return nn
                end
            end
            i+=1
        end
        if reflexlinkedgecount>0 &&
                ((m.fliplinklevel<0 && level<m.autofliplinklevel) ||
                 (m.fliplinklevel>=0 && level<m.fliplinklevel))
            i=0
            while i<n
                cur=abtets[i+1];prev=abtets[mod(i-1,n)+1]
                if elemcounter(m,cur)>1 || elemcounter(m,prev)>1
                    i+=1;continue
                end
                pc=apex(m,cur)
                if pc==DUMMY
                    i+=1;continue
                end
                pd=apex(m,abtets[mod(i+1,n)+1]);pe=apex(m,prev)
                if pd==DUMMY || pe==DUMMY
                    i+=1;continue
                end
                edgepivot=0;flipedge=TNULL
                if orient3d(m,pb,pc,pd,pe)<=0
                    flipedge=enext(cur);edgepivot=1
                end
                if edgepivot==0 && orient3d(m,pc,pa,pd,pe)<=0
                    flipedge=eprev(cur);edgepivot=2
                end
                if edgepivot==0
                    i+=1;continue
                end
                if m.checksubsegflag && issubseg(m,flipedge)
                    i+=1;continue
                end
                flipedge=esym(flipedge)
                n1=0;j=0;spintet=flipedge
                while true
                    n1+=1;j+=elemcounter(m,spintet)
                    spintet=fnext(m,spintet)
                    spintet.tet==flipedge.tet && break
                end
                n1<3 && throw(Unsupported("flipnm: an edge star with fewer than three tetrahedra"))
                if j>2 || (m.flipstarsize>0 && n1>m.flipstarsize)
                    i+=1;continue
                end
                tmpabtets=Vector{Tri}(undef,n1)
                j=0;spintet=flipedge
                while true
                    tmpabtets[j+1]=spintet
                    increaseelemcounter!(m,spintet)
                    j+=1
                    spintet=fnext(m,spintet)
                    spintet.tet==flipedge.tet && break
                end
                nn=flipnm!(m,tmpabtets,n1,level+1,edgepivot,fc)
                if nn==2
                    spintet=edgepivot==1 ? enext(esym(enext(tmpabtets[1]))) :
                                           eprev(esym(eprev(tmpabtets[2])))
                    increaseelemcounter!(m,spintet)
                    abtets[mod(i-1,n)+1]=spintet
                    for j in i:n-2
                        abtets[j+1]=abtets[j+2]
                    end
                    push!(m.stars,tmpabtets)
                    abtets[n]=Tri(length(m.stars),edgepivot|(1<<5)|(i<<6)|(n1<<19))
                    tmpabtets[1]=Tri(Int(pc),1<<5)
                    nn=flipnm!(m,abtets,n-1,level,abedgepivot,fc)
                    nn==2 && return nn
                    if fc.unflip
                        if edgepivot==1
                            tmpabtets[1]=eprev(esym(eprev(abtets[mod(i-1,n-1)+1])))
                            tmpabtets[2]=fsym(m,tmpabtets[1])
                        else
                            tmpabtets[2]=enext(esym(enext(abtets[mod(i-1,n-1)+1])))
                            tmpabtets[1]=fsym(m,tmpabtets[2])
                        end
                        flipnm_post!(m,tmpabtets,n1,2,edgepivot,fc)
                        for j in n-2:-1:i
                            abtets[j+2]=abtets[j+1]
                        end
                        if edgepivot==1
                            f0=esym(enext(tmpabtets[2]));f1=eprev(esym(tmpabtets[1]))
                        else
                            f0=esym(eprev(tmpabtets[2]));f1=enext(esym(tmpabtets[1]))
                        end
                        increaseelemcounter!(m,f0);increaseelemcounter!(m,f1)
                        abtets[mod(i-1,n)+1]=f0
                        abtets[i+1]=f1
                        nn+=1
                    end
                    fc.unflip || return nn
                else
                    fc.unflip || flipnm_post!(m,tmpabtets,n1,nn,edgepivot,fc)
                    for j in 1:nn
                        decreaseelemcounter!(m,tmpabtets[j])
                    end
                end
                i+=1
            end
        end
    else
        hullflag=0
        if apex(m,abtets[1])==DUMMY
            pc=apex(m,abtets[2]);pd=apex(m,abtets[3]);pe=apex(m,abtets[1]);hullflag=1
        elseif apex(m,abtets[2])==DUMMY
            pc=apex(m,abtets[3]);pd=apex(m,abtets[1]);pe=apex(m,abtets[2]);hullflag=2
        else
            pc=apex(m,abtets[1]);pd=apex(m,abtets[2]);pe=apex(m,abtets[3])
            hullflag=pe==DUMMY ? 3 : 0
        end
        reducflag=false;rejflag=false
        if hullflag==0
            if orient3d(m,pd,pc,pe,pa)<0 && orient3d(m,pc,pd,pe,pb)<0
                reducflag=true
            end
        else
            throw(Unsupported("flipnm: a 3-to-2 flip at a hull edge"))
        end
        if reducflag
            if m.checksubfaceflag
                nn=0;edgepivot=-1
                for j in 0:2
                    if issubface(m,abtets[j+1])
                        nn+=1
                    else
                        edgepivot=j
                    end
                end
                if nn==1
                    rejflag=true
                elseif nn==2
                    spintet=eorgoppo(abtets[mod(edgepivot+1,3)+1])
                    if issubface(m,spintet)
                        rejflag=true
                    elseif issubface(m,esym(spintet))
                        rejflag=true
                    end
                elseif nn==3
                    throw(Unsupported("flipnm: three subfaces at one edge"))
                end
            end
            if !rejflag && fc.checkflipeligibility
                rejflag=checkflipeligibility(m,2,pc,pd,pe,pb,pa,level,abedgepivot,fc)
            end
            if !rejflag
                flip32!(m,abtets,0,hullflag,fc)
                if fc.remove_ndelaunay_edge && level==0
                    # tetprismvol() is zero in Gmsh (orient4d is a stub), so
                    # the objective never improves: flip back.
                    if fc.tetprism_vol_sum>=0.0 || abs(fc.tetprism_vol_sum)<fc.bak_tetprism_vol
                        flip23!(m,abtets,0,hullflag,fc)
                        for j in 1:3;increaseelemcounter!(m,abtets[j]);end
                        return 3
                    end
                end
                if fc.collectnewtets
                    if level==0
                        push!(m.cavetetlist,abtets[1]);push!(m.cavetetlist,abtets[2])
                    else
                        push!(m.cavetetlist,abedgepivot==1 ? abtets[2] : abtets[1])
                    end
                end
                return 2
            end
        end
    end
    return n
end

function flipnm_post!(m,abtets,n,nn,abedgepivot,fc)
    fliptets=m.F3
    if nn==2
        if fc.unflip
            flip23!(m,abtets,0,1,fc)
            if fc.collectnewtets
                resize!(m.cavetetlist,length(m.cavetetlist)-(abedgepivot==0 ? 2 : 1))
            end
        end
        nn+=1
    end
    for i in nn:n-1
        fliptype=(abtets[i+1].ver>>4)&3
        if fliptype==1
            t=abtets[i+1].ver>>6
            if fc.unflip
                fliptets[1]=enext(esym(eprev(abtets[mod(t-1,i)+1])))
                fliptets[2]=fnext(m,fliptets[1])
                fliptets[3]=fnext(m,fliptets[2])
                flip32!(m,fliptets,0,1,fc)
                for j in i-1:-1:t
                    abtets[j+2]=abtets[j+1]
                end
                abtets[mod(t-1,i+1)+1]=esym(fliptets[2])
                abtets[t+1]=fliptets[1]
                fc.collectnewtets && resize!(m.cavetetlist,length(m.cavetetlist)-2)
            end
        elseif fliptype==2
            tmpabtets=m.stars[abtets[i+1].tet]
            n1=(abtets[i+1].ver>>19)&8191
            edgepivot=abtets[i+1].ver&3
            t=(abtets[i+1].ver>>6)&8191
            if fc.unflip
                if edgepivot==1
                    tmpabtets[1]=eprev(esym(eprev(abtets[mod(t-1,i)+1])))
                    tmpabtets[2]=fsym(m,tmpabtets[1])
                else
                    tmpabtets[2]=enext(esym(enext(abtets[mod(t-1,i)+1])))
                    tmpabtets[1]=fsym(m,tmpabtets[2])
                end
                flipnm_post!(m,tmpabtets,n1,2,edgepivot,fc)
                for j in i-1:-1:t
                    abtets[j+2]=abtets[j+1]
                end
                if edgepivot==1
                    f0=esym(enext(tmpabtets[2]));f1=eprev(esym(tmpabtets[1]))
                else
                    f0=esym(eprev(tmpabtets[2]));f1=enext(esym(tmpabtets[1]))
                end
                abtets[mod(t-1,i+1)+1]=f0
                abtets[t+1]=f1
            else
                flipnm_post!(m,tmpabtets,n1,2,edgepivot,fc)
            end
        end
    end
    return 1
end

function removeedgebyflips!(m,flipedge,fc)
    m.checksubsegflag && issubseg(m,flipedge) && return 0,flipedge
    n=0;spintet=flipedge
    while true
        n+=1;spintet=fnext(m,spintet)
        spintet.tet==flipedge.tet && break
    end
    n<3 && throw(Unsupported("removeedgebyflips: an inverted edge star"))
    m.flipstarsize>0 && n>m.flipstarsize && return 0,flipedge
    abtets=Vector{Tri}(undef,n)
    spintet=flipedge;i=0
    while true
        abtets[i+1]=spintet;setelemcounter!(m,spintet,1);i+=1
        spintet=fnext(m,spintet)
        spintet.tet==flipedge.tet && break
    end
    nn=flipnm!(m,abtets,n,0,0,fc)
    if nn>2
        for i in 1:nn;setelemcounter!(m,abtets[i],0);end
        flipedge=abtets[1]
    end
    bak=fc.unflip;fc.unflip=false
    flipnm_post!(m,abtets,n,nn,0,fc)
    fc.unflip=bak
    empty!(m.stars)
    return nn,flipedge
end

function lawsonflip3d!(m,fc)
    F=Vector{Tri}(undef,5)
    totalcount=0;sliver_peels=0;loopcount=0
    while true
        flipcount=0
        while !isempty(m.flipstack)
            F[1]=pop!(m.flipstack)
            isdeadtet(m,F[1]) && continue
            facemarked(m,F[1]) || continue
            unmarkface!(m,F[1])
            ishulltet(m,F[1]) && continue
            F[2]=fsym(m,F[1])
            if ishulltet(m,F[2])
                if m.nonconvex
                    checksh=tspivot(m,F[1])
                    checksh.sh==0 && throw(Unsupported("lawsonflip3d: a hull face without a subface"))
                    for _ in 1:3
                        if !isshsubseg(m,checksh)
                            casingout=spivot(m.subfaces,checksh)
                            sorg(m.subfaces,checksh)!=sdest(m.subfaces,casingout) &&
                                (casingout=sesym(casingout))
                            neightet=stpivot(m,casingout)
                            if neightet.tet==F[1].tet
                                hulltet=fsym(m,edestoppo(neightet))
                                if oppo(m,hulltet)==DUMMY
                                    pe=org(m,neightet)
                                    pointtype(m,pe) in (FREEFACETVERTEX,FREESEGVERTEX) &&
                                        throw(Unsupported("lawsonflip3d: removevertexbyflips"))
                                else
                                    hulltet=fsym(m,eorgoppo(neightet))
                                    if oppo(m,hulltet)==DUMMY
                                        pd=dest(m,neightet)
                                        pointtype(m,pd) in (FREEFACETVERTEX,FREESEGVERTEX) &&
                                            throw(Unsupported("lawsonflip3d: removevertexbyflips"))
                                    else
                                        F[1]=neightet;F[2]=fnext(m,F[1]);F[3]=fnext(m,F[2])
                                        flip32!(m,F,0,1,fc)
                                        sliver_peels+=1
                                    end
                                end
                                break
                            end
                        end
                        checksh=senext(checksh)
                    end
                end
                continue
            end
            m.checksubfaceflag && issubface(m,F[1]) && continue
            pts=F[2].tet
            sign=insphere_s(m,tvert(m,pts,4),tvert(m,pts,5),tvert(m,pts,6),tvert(m,pts,7),oppo(m,F[1]))
            if sign<0
                pd=oppo(m,F[1]);pe=oppo(m,F[2])
                len3=distance(m.xyz[pd],m.xyz[pe]);len3=len3*len3*len3
                ori=0.0
                for _ in 1:3
                    ori=orient3d(m,org(m,F[1]),dest(m,F[1]),pd,pe)
                    if ori>0
                        F[3]=esym(F[1]);F[4]=esym(F[2])
                        if issubface(m,F[3]) || issubface(m,F[4])
                            vol=orient3dfast(m.xyz[org(m,F[1])],m.xyz[dest(m,F[1])],m.xyz[pd],m.xyz[pe])
                            abs(vol)/len3<m.epsilon && (ori=0.0)
                        end
                    end
                    ori<=0 && break
                    F[1]=enext(F[1]);F[2]=eprev(F[2])
                end
                if ori>0
                    flip23!(m,F,0,0,fc)
                    flipcount+=1
                    continue
                else
                    m.checksubsegflag && issubseg(m,F[1]) && continue
                    F[1]=esym(F[1])
                    for i in 1:3
                        F[i+1]=fnext(m,F[i])
                    end
                    if F[4].tet==F[1].tet
                        flip32!(m,F,0,0,fc)
                        flipcount+=1
                        continue
                    else
                        F[5]=fnext(m,F[4])
                        if F[5].tet==F[1].tet
                            if m.nonconvex && apex(m,F[4])==DUMMY
                                ori=0.0
                            end
                            if ori==0
                                F[1]=esym(F[1])
                                flip23!(m,F,0,0,fc)
                                F[2]=fnext(m,F[4]);F[3]=fnext(m,F[2])
                                flip32!(m,F,1,apex(m,F[4])==DUMMY ? 1 : 0,fc)
                                flipcount+=1
                                continue
                            end
                        end
                    end
                end
                push!(m.unflipqueue,BadFace(F[1],org(m,F[1]),dest(m,F[1]),apex(m,F[1])))
            end
        end
        totalcount+=flipcount
        isempty(m.unflipqueue) && break
        flipcount==0 && break
        for bface in m.unflipqueue
            t=bface.tt
            if !isdeadtet(m,t) && org(m,t)==bface.forg && dest(m,t)==bface.fdest &&
                    apex(m,t)==bface.fapex
                flippush!(m,t)
            end
        end
        empty!(m.unflipqueue)
        loopcount>10_000 && throw(Unsupported("lawsonflip3d does not terminate"))
        loopcount+=1
    end
    return totalcount+sliver_peels
end
