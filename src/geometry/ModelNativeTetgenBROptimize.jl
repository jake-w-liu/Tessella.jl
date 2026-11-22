# optimizemesh() with Gmsh's options (-Y without -q: optlevel 2, optscheme 7,
# largest allowed dihedral 179 degrees, sliver and smoothing bounds 179.999
# degrees): improvequalitybyflips, improvequalitybysmoothing (Steiner points
# only), removeslivers (splitsliver + insertpoint into the edge star).

@inline _cosdeg(d)=cos(d/180.0*Float64(pi))

function _badface(m,slot)
    a=tvert(m,slot,4);b=tvert(m,slot,5);c=tvert(m,slot,6);d=tvert(m,slot,7)
    _,cent,key=tetalldihedral(m,a,b,c,d)
    return BadFace(Tri(slot,11),key,cent,a,b,c,d)
end

function improvequalitybyflips!(m)
    flipqueue=BadFace[];nextflipqueue=BadFace[]
    bak=(m.autofliplinklevel,m.fliplinklevel,m.flipstarsize)
    m.autofliplinklevel=1;m.fliplinklevel=-1;m.flipstarsize=10
    fc=Flip()
    fc.remove_large_angle=true;fc.unflip=true;fc.collectnewtets=true
    fc.checkflipeligibility=true
    totalremcount=0
    flipqueue,m.unflipqueue=m.unflipqueue,flipqueue
    while !isempty(flipqueue)
        remcount=0
        while !isempty(flipqueue)
            for k in eachindex(flipqueue)
                bface=flipqueue[k]
                ok,tt=gettetrahedron!(m,bface.forg,bface.fdest,bface.fapex,bface.foppo,bface.tt)
                ok || continue
                key=bface.key;cent=bface.cent
                forg=bface.forg;fdest=bface.fdest;fapex=bface.fapex;foppo=bface.foppo
                if tt.ver!=11
                    forg=tvert(m,tt.tet,4);fdest=tvert(m,tt.tet,5)
                    fapex=tvert(m,tt.tet,6);foppo=tvert(m,tt.tet,7)
                    _,cent,key=tetalldihedral(m,forg,fdest,fapex,foppo)
                    tt=Tri(tt.tet,11)
                end
                if key==0
                    _,cent,key=tetalldihedral(m,tvert(m,tt.tet,4),tvert(m,tt.tet,5),
                        tvert(m,tt.tet,6),tvert(m,tt.tet,7))
                end
                remflag=false
                for i in 1:6
                    remflag && break
                    if cent[i]<m.cosmaxdihed
                        tt=Tri(tt.tet,EDGE2VER[i])
                        fc.cosdihed_in=cent[i];fc.cosdihed_out=0.0
                        n,tt=removeedgebyflips!(m,tt,fc)
                        if n==2
                            remflag=true
                            if fc.cosdihed_out<m.cosmaxdihed
                                for parytet in m.cavetetlist
                                    isdeadtet(m,parytet) && continue
                                    tvert(m,parytet.tet,7)==DUMMY && continue
                                    nb=_badface(m,parytet.tet)
                                    nb.key<m.cosmaxdihed && push!(nextflipqueue,nb)
                                end
                            end
                            empty!(m.cavetetlist)
                            remcount+=1
                        end
                    end
                end
                if !remflag
                    push!(m.unflipqueue,BadFace(tt,key,cent,forg,fdest,fapex,foppo))
                end
            end
            empty!(flipqueue)
            flipqueue,nextflipqueue=nextflipqueue,flipqueue
        end
        totalremcount+=remcount
        if !isempty(m.unflipqueue)
            m.autofliplinklevel>=2 && break
            m.autofliplinklevel+=1
        end
        flipqueue,m.unflipqueue=m.unflipqueue,flipqueue
    end
    m.autofliplinklevel,m.fliplinklevel,m.flipstarsize=bak
    return totalremcount
end

# smoothpoint(): returns the number of accepted moves and the final position.
function smoothpoint!(m,start::NTuple{3,Float64},list::Vector{Tri},ccw::Bool,opm::Opt)
    numdirs=min(length(list),opm.numofsearchdirs)
    opm.imprval=opm.initval
    iter=0
    bestpt=startpt=start
    minval=0.0
    while true
        oldval=opm.imprval
        for i in 0:numdirs-1
            k=randomnation!(m,length(list)-i)
            face=list[k+1]
            pa=m.xyz[org(m,face)];pb=m.xyz[dest(m,face)];pc=m.xyz[apex(m,face)]
            fcent=ntuple(j->(pa[j]+pb[j]+pc[j])/3.0,3)
            nextpt=ntuple(j->startpt[j]+opm.searchstep*(fcent[j]-startpt[j]),3)
            j=0
            while j<length(list)
                f=list[j+1]
                a,b=ccw ? (org(m,f),dest(m,f)) : (dest(m,f),org(m,f))
                pa=m.xyz[a];pb=m.xyz[b];pc=m.xyz[apex(m,f)]
                ori=Float64(orient3(pa,pb,pc,nextpt))
                val=0.0
                if ori<0.0
                    if opm.max_min_volume
                        val=-orient3dfast(pa,pb,pc,nextpt)
                    elseif opm.min_max_dihedangle
                        _,_,maxcosd=tetalldihedral(pa,pb,pc,nextpt)
                        maxcosd< -1 && (maxcosd=-1.0)
                        val=maxcosd+1.0
                    end
                else
                    if opm.max_min_volume
                        val=-orient3dfast(pa,pb,pc,nextpt)
                    else
                        break
                    end
                end
                val<=opm.imprval && break
                minval=j==0 ? val : (val<minval ? val : minval)
                j+=1
            end
            if j==length(list)
                opm.imprval=minval
                bestpt=nextpt
            end
            jj=length(list)-i-1
            list[k+1],list[jj+1]=list[jj+1],list[k+1]
        end
        diff=opm.imprval-oldval
        if diff>0.0
            startpt=bestpt
            iter+=1
            opm.maxiter>0 && iter>=opm.maxiter && break
        else
            break
        end
    end
    iter>0 && (opm.smthiter=iter)
    return iter,startpt
end

function getvertexstar!(m,searchpt,tetlist)
    searchtet=enextesym(point2tetorg(m,searchpt))
    infect!(m,searchtet);push!(tetlist,searchtet)
    neightet=esym(fsym(m,esym(searchtet)))
    infect!(m,neightet);push!(tetlist,neightet)
    i=1
    while i<=length(tetlist)
        searchtet=tetlist[i]
        for _ in 1:2
            searchtet=enext(searchtet)
            neightet=fsym(m,esym(searchtet))
            if !infected(m,neightet)
                neightet=esym(neightet)
                infect!(m,neightet);push!(tetlist,neightet)
            end
        end
        i+=1
    end
    for t in tetlist;uninfect!(m,t);end
    return length(tetlist)
end

function improvequalitybysmoothing!(m,opm)
    flipqueue=BadFace[]
    flipqueue,m.unflipqueue=m.unflipqueue,flipqueue
    totalsmtcount=0;iter=0
    while !isempty(flipqueue)
        smtcount=0
        for k in eachindex(flipqueue)
            bface=flipqueue[k]
            ok,tt=gettetrahedron!(m,bface.forg,bface.fdest,bface.fapex,bface.foppo,bface.tt)
            (ok && !marktested(m,tt)) || continue
            ppt=(tvert(m,tt.tet,4),tvert(m,tt.tet,5),tvert(m,tt.tet,6),tvert(m,tt.tet,7))
            _,_,key=tetalldihedral(m,ppt...)
            key<m.cossmtdihed || continue
            smtflag=false
            opm.initval=key+1.0
            for i in 1:4
                smtflag && break
                pointtype(m,ppt[i])==FREEVOLVERTEX || continue
                getvertexstar!(m,ppt[i],m.cavetetlist)
                opm.searchstep=0.001
                it,newpt=smoothpoint!(m,m.xyz[ppt[i]],m.cavetetlist,true,opm)
                smtflag=it>0
                if smtflag
                    m.xyz[ppt[i]]=newpt
                    while opm.smthiter==opm.maxiter
                        opm.searchstep*=10.0
                        opm.initval=opm.imprval
                        opm.smthiter=0
                        it,newpt=smoothpoint!(m,m.xyz[ppt[i]],m.cavetetlist,true,opm)
                        it>0 && (m.xyz[ppt[i]]=newpt)
                    end
                    smtcount+=1
                    if opm.imprval-1.0<m.cossmtdihed
                        for parytet in m.cavetetlist
                            marktested(m,parytet) && continue
                            q=(tvert(m,parytet.tet,4),tvert(m,parytet.tet,5),
                               tvert(m,parytet.tet,6),tvert(m,parytet.tet,7))
                            _,_,qkey=tetalldihedral(m,q...)
                            if qkey<m.cossmtdihed
                                marktest!(m,parytet)
                                push!(m.unflipqueue,BadFace(Tri(parytet.tet,11),0.0,
                                    ntuple(_->0.0,6),q[1],q[2],q[3],q[4]))
                            end
                        end
                    end
                end
                empty!(m.cavetetlist)
            end
            if !smtflag
                marktest!(m,tt)
                push!(m.unflipqueue,BadFace(Tri(tt.tet,11),0.0,ntuple(_->0.0,6),ppt...))
            end
        end
        empty!(flipqueue)
        for b in m.unflipqueue;unmarktest!(m,b.tt);end
        totalsmtcount+=smtcount
        smtcount==0 && break
        iter+=1
        iter==2 && break
        flipqueue,m.unflipqueue=m.unflipqueue,flipqueue
    end
    return totalsmtcount
end

# insertpoint() for a Steiner point whose cavity is the given infected star
# (ivf.iloc = INSTAR, no Bowyer-Watson growth, no boundary splitting).
function insertpoint_instar!(m,insertpt,chkencflag)
    for cavetet in m.caveoldtetlist,j in 0:3
        neightet=rawnb(m,cavetet.tet,j)
        infected(m,neightet) || push!(m.cavebdrylist,Tri(neightet.tet,EPIVOT[neightet.ver+1]))
    end
    S=m.subsegs;F=m.subfaces
    if m.checksubsegflag
        for cavetet in m.caveoldtetlist,j in 1:6
            x=m.tseg[6*(cavetet.tet-1)+j]
            x==0 && continue
            seg=sdec(x)
            if !sinfected(S,seg)
                sinfect!(S,seg);push!(m.cavetetseglist,seg)
            end
        end
        for seg in m.cavetetseglist;suninfect!(S,seg);end
    end
    if m.checksubfaceflag
        for cavetet in m.caveoldtetlist,j in 1:4
            x=m.tsub[4*(cavetet.tet-1)+j]
            x==0 && continue
            sh=sdec(x)
            if !sinfected(F,sh)
                sinfect!(F,sh);push!(m.cavetetshlist,sh)
            end
        end
        for sh in m.cavetetshlist;suninfect!(F,sh);end
    end
    if m.checksubsegflag
        i=1
        while i<=length(m.cavetetseglist)
            seg=m.cavetetseglist[i]
            if !smarktested(S,seg)
                j=0;k=0;neineitet=TNULL
                neightet=sstpivot1(m,seg);spintet=neightet
                while true
                    j+=1
                    if !infected(m,spintet)
                        neineitet=spintet
                    else
                        k+=1
                    end
                    spintet=fnext(m,spintet)
                    spintet.tet==neightet.tet && break
                end
                if k==0
                    m.cavetetseglist[i]=m.cavetetseglist[end];pop!(m.cavetetseglist)
                    continue
                elseif k<j
                    sstbond1!(m,seg,neineitet)
                else
                    throw(Unsupported("insertpoint: a segment inside the Steiner cavity"))
                end
            else
                sinfect!(S,seg)
            end
            i+=1
        end
    end
    if m.checksubfaceflag
        i=1
        while i<=length(m.cavetetshlist)
            sh=m.cavetetshlist[i]
            if !smarktested(F,sh)
                k=0;checksh=SNULL
                for _ in 1:2
                    neightet=stpivot(m,sh)
                    if !infected(m,neightet)
                        checksh=sh
                    else
                        k+=1
                    end
                    sh=sesym(sh)
                end
                if k==0
                    m.cavetetshlist[i]=m.cavetetshlist[end];pop!(m.cavetetshlist)
                    continue
                elseif k==1
                    m.cavetetshlist[i]=checksh
                else
                    throw(Unsupported("insertpoint: a subface inside the Steiner cavity"))
                end
            else
                sinfect!(F,sh)
            end
            i+=1
        end
    end
    newtet=TNULL
    for i in eachindex(m.cavebdrylist)
        neightet=m.cavebdrylist[i]
        unmarktest!(m,neightet)
        oldtet=fsym(m,neightet)
        newtet=maketetrahedron!(m)
        if apex(m,neightet)!=DUMMY
            setorg!(m,newtet,dest(m,neightet));setdest!(m,newtet,org(m,neightet))
            setapex!(m,newtet,apex(m,neightet));setoppo!(m,newtet,insertpt)
        else
            setorg!(m,newtet,org(m,neightet));setdest!(m,newtet,dest(m,neightet))
            setapex!(m,newtet,insertpt);setoppo!(m,newtet,DUMMY)
            newtet=esym(newtet)
        end
        bond!(m,newtet,neightet)
        m.cavebdrylist[i]=oldtet
    end
    setpoint2tet!(m,insertpt,Tri(newtet.tet,0))
    empty!(m.cavetetlist)
    for i in eachindex(m.cavebdrylist)
        oldtet=m.cavebdrylist[i]
        neightet=fsym(m,oldtet)
        newtet=fsym(m,neightet)
        for _ in 1:3
            neightet=esym(newtet)
            if m.nb[_nbi(neightet)]==0
                spintet=oldtet
                while true
                    spintet=fnext(m,spintet)
                    infected(m,spintet) || break
                end
                newneitet=esym(fsym(m,spintet))
                bond!(m,neightet,newneitet)
            end
            setpoint2tet!(m,org(m,newtet),Tri(newtet.tet,0))
            newtet=enext(newtet);oldtet=enext(oldtet)
        end
        m.cavebdrylist[i]=newtet
    end
    if m.checksubfaceflag
        for i in eachindex(m.cavetetshlist)
            sh=m.cavetetshlist[i]
            sinfected(F,sh) && continue
            neightet=stpivot(m,sh)
            spintet=fsym(m,neightet)
            sh=sesym(sh);m.cavetetshlist[i]=sh
            tsbond!(m,spintet,sh)
        end
    end
    if m.checksubsegflag
        for seg in m.cavetetseglist
            sinfected(S,seg) && continue
            neightet=sstpivot1(m,seg);spintet=neightet
            while true
                tssbond1!(m,spintet,seg)
                spintet=fnext(m,spintet)
                spintet.tet==neightet.tet && break
            end
        end
    end
    if (chkencflag&4)!=0
        for t in m.cavebdrylist;enqueuetetrahedron!(m,t);end
    end
    for t in m.caveoldtetlist
        tetrahedrondealloc!(m,t.tet)
    end
    empty!(m.caveoldtetlist);empty!(m.cavebdrylist);empty!(m.cavetetlist)
    empty!(m.cavetetseglist);empty!(m.caveencseglist)
    empty!(m.cavetetshlist);empty!(m.caveencshlist)
    return true
end

function splitsliver!(m,slitet,cosd,chkencflag)
    searchtet=edestoppo(slitet)
    issubseg(m,searchtet) && return false
    spintet=searchtet;n=0
    while true
        ishulltet(m,spintet) && break
        n+=1
        spintet=fnext(m,spintet)
        spintet.tet==searchtet.tet && break
    end
    ishulltet(m,spintet) && return false
    abtets=Vector{Tri}(undef,n)
    spintet=searchtet
    for i in 1:n
        abtets[i]=spintet;spintet=fnext(m,spintet)
    end
    for i in 1:n
        push!(m.cavetetlist,esym(eprev(abtets[i])))
        push!(m.cavetetlist,esym(enext(abtets[i])))
    end
    pa=m.xyz[org(m,abtets[1])];pb=m.xyz[dest(m,abtets[1])]
    smtpt=ntuple(i->0.5*(pa[i]+pb[i]),3)
    opm=Opt()
    opm.min_max_dihedangle=true
    opm.initval=cosd+1.0
    opm.numofsearchdirs=20
    opm.searchstep=0.001
    opm.maxiter=100
    success,newpt=smoothpoint!(m,smtpt,m.cavetetlist,true,opm)
    success>0 && (smtpt=newpt)
    if success>0
        while opm.smthiter==opm.maxiter
            opm.searchstep*=10.0
            opm.initval=opm.imprval
            opm.smthiter=0
            it,newpt=smoothpoint!(m,smtpt,m.cavetetlist,true,opm)
            it>0 && (smtpt=newpt)
        end
    end
    empty!(m.cavetetlist)
    success>0 || return false
    steinerpt=makepoint!(m,FREEVOLVERTEX)
    m.xyz[steinerpt]=smtpt
    for i in 1:n
        infect!(m,abtets[i]);push!(m.caveoldtetlist,abtets[i])
    end
    insertpoint_instar!(m,steinerpt,chkencflag)
    m.steiner+=1
    return true
end

function removeslivers!(m,chkencflag)
    flipqueue=BadFace[]
    flipqueue,m.unflipqueue=m.unflipqueue,flipqueue
    totalsptcount=0;iter=0
    while !isempty(flipqueue)
        sptcount=0
        for k in eachindex(flipqueue)
            bface=flipqueue[k]
            ok,tt=gettetrahedron!(m,bface.forg,bface.fdest,bface.fapex,bface.foppo,bface.tt)
            bface=BadFace(tt,bface.key,bface.cent,bface.forg,bface.fdest,bface.fapex,bface.foppo)
            ok || continue
            if bface.key==0 || tt.ver!=11
                _,cent,key=tetalldihedral(m,tvert(m,tt.tet,4),tvert(m,tt.tet,5),
                    tvert(m,tt.tet,6),tvert(m,tt.tet,7))
                bface=BadFace(tt,key,cent,bface.forg,bface.fdest,bface.fapex,bface.foppo)
            end
            bface.key<m.cosslidihed || continue
            split=false
            for j in 1:6
                if bface.cent[j]<m.cosslidihed
                    if splitsliver!(m,Tri(tt.tet,EDGE2VER[j]),bface.cent[j],chkencflag)
                        sptcount+=1;split=true;break
                    end
                end
            end
            if split
                for parytet in m.badtetrahedrons
                    unmarktest2!(m,parytet)
                    nb=_badface(m,parytet.tet)
                    if nb.key<m.cosslidihed
                        push!(m.unflipqueue,BadFace(Tri(parytet.tet,11),nb.key,nb.cent,
                            nb.forg,nb.fdest,nb.fapex,nb.foppo))
                    end
                end
                empty!(m.badtetrahedrons)
            else
                push!(m.unflipqueue,bface)
            end
        end
        empty!(flipqueue)
        totalsptcount+=sptcount
        sptcount==0 && break
        iter+=1
        iter==2 && break
        flipqueue,m.unflipqueue=m.unflipqueue,flipqueue
    end
    return totalsptcount
end

function optimizemesh!(m)
    optpasses=(1<<2)-1
    m.cosmaxdihed=_cosdeg(179.0)
    m.cossmtdihed=_cosdeg(179.999)
    m.cosslidihed=_cosdeg(179.999)
    for slot in 1:m.tmax
        (tvert(m,slot,4)==0 || tvert(m,slot,7)==DUMMY) && continue
        b=_badface(m,slot)
        b.key<m.cosmaxdihed && push!(m.unflipqueue,b)
    end
    improvequalitybyflips!(m)
    if !isempty(m.unflipqueue)
        opm=Opt()
        opm.min_max_dihedangle=true
        opm.numofsearchdirs=10
        opm.maxiter=30
        chkencflag=4
        iter=0
        while iter<optpasses
            smtcount=0;sptcount=0;remcount=0
            smtcount+=improvequalitybysmoothing!(m,opm)
            smtcount>0 && (remcount=improvequalitybyflips!(m))
            if !isempty(m.unflipqueue)
                sptcount+=removeslivers!(m,chkencflag)
                sptcount>0 && (remcount=improvequalitybyflips!(m))
            end
            if !isempty(m.unflipqueue)
                remcount>0 || break
                iter+=1
            else
                break
            end
        end
        empty!(m.badtetrahedrons)
    end
    empty!(m.unflipqueue)
    return nothing
end

"""
    BR.region(kernel, facets; seed=1, epsilon=1e-12, facet_overlap=0.1)

Run Gmsh's TetGen/BR pass on one Region: reconstruct the kernel's Delaunay
tetrahedralization (its eight enclosing-box points included), recover the
boundary `facets` (native point rows, in Gmsh's GFace and triangle order),
carve the exterior, recover Delaunayness and optimize. Returns the tetrahedra
in pool order and the Steiner points in point-pool order, as
`(cells, steiner, seed)`. A boundary that needs Steiner points or flips to be
recovered throws `Unsupported`.
"""
function region(kernel,facets;seed=1,epsilon=1e-12,facet_overlap=0.1)
    m=_emptytg(kernel.points,length(kernel.points),seed,epsilon)
    _reconstruct!(m,kernel)
    _surface!(m,facets,facet_overlap/180.0*Float64(pi))
    _recoverboundary!(m)
    _carveholes!(m)
    _suppresssteinerpoints!(m)
    _recoverdelaunay!(m)
    optimizemesh!(m)
    keep=_regioncomponents(m)
    cells=NTuple{4,Int32}[]
    used=falses(length(m.xyz))
    for slot in 1:m.tmax
        keep[slot] || continue
        cell=(tvert(m,slot,4),tvert(m,slot,5),tvert(m,slot,6),tvert(m,slot,7))
        push!(cells,cell)
        for p in cell;used[p]=true;end
    end
    steiner=Tuple{Int32,NTuple{3,Float64}}[]
    for p in m.ninput+1:length(m.xyz)
        t=pointtype(m,p)
        used[p] && t in (FREESEGVERTEX,FREEFACETVERTEX,FREEVOLVERTEX) && push!(steiner,(Int32(p),m.xyz[p]))
    end
    return (cells=cells,steiner=steiner,seed=m.seed,mesh=m)
end

# TetGen keeps every tetrahedron enclosed by subfaces, holes included; Gmsh's
# insertVerticesInRegion then classifies the face-connected components and
# deletes those whose bounding faces are not the Region's ("void region").
# For one Region every subface bounds its own component on one side, so a
# component is the Region's when it touches all subfaces.
function _regioncomponents(m)
    comp=zeros(Int32,m.tmax)
    sides=Int[]
    stack=Int[]
    for slot in 1:m.tmax
        (tvert(m,slot,4)==0 || tvert(m,slot,7)==DUMMY || comp[slot]!=0) && continue
        push!(sides,0);c=Int32(length(sides))
        comp[slot]=c;push!(stack,slot)
        while !isempty(stack)
            s=pop!(stack)
            for f in 0:3
                if m.tsub[4*(s-1)+f+1]!=0
                    sides[c]+=1
                    continue
                end
                n=rawnb(m,s,f).tet
                (n==0 || comp[n]!=0 || tvert(m,n,7)==DUMMY) && continue
                comp[n]=c;push!(stack,n)
            end
        end
    end
    keep=falses(m.tmax)
    found=false
    for slot in 1:m.tmax
        c=comp[slot]
        if c!=0 && sides[c]==m.subfaces.items
            keep[slot]=true;found=true
        end
    end
    found || throw(Unsupported("no tetrahedra are bounded by the whole Region boundary"))
    return keep
end
