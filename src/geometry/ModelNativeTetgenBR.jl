# SPDX-FileCopyrightText: 2015-2016 Weierstrass Institute for Applied Analysis and Stochastics
# SPDX-FileCopyrightText: 2026 Jake Liu
# SPDX-License-Identifier: MIT AND GPL-2.0-or-later
# Port of the part of Gmsh 4.15.2's TetGen/BR (src/mesh/tetgenBR.cxx and
# meshGRegionBoundaryRecovery.cpp) that runs on a Region after its Delaunay
# tetrahedralization: mesh reconstruction, boundary recovery when the boundary
# is already present, carveholes, recoverdelaunay and optimizemesh (flips,
# sliver splitting and Steiner-point smoothing).
# Julia adaptation: Tessella contributors, 2026.
# Preserve the upstream agreement notice in LICENSES/NOTICES.md#tetgenbr.
# Licenses apply to their respective portions; see LICENSES/README.md.

# The port keeps TetGen's memory layout where it decides results: tetrahedra,
# subfaces, segments and points live in pools whose slots are reused last-freed
# first, traversals run in slot order, and every random draw comes from the
# same linear congruential generator in the same order. A tetrahedron handle is
# a slot and a version 0..11; an encoded handle is `slot<<4 | version`. A
# shell-face handle is a slot and a version 0..5, encoded `slot<<3 | version`.
# Branches of TetGen that this region pass never reaches with a closed,
# already-conforming boundary throw `Unsupported`, and the caller falls back.
module BR
using ....Predicates: orient3, insphere
import ..Unsupported

# ---------------------------------------------------------------------------
# Tables (tetgenmesh::inittables), 0-based versions, accessed as T[v+1].
const ESYMTBL=(9,6,11,4,3,7,1,5,10,0,8,2)
const ORGPIVOT=(7,7,5,5,6,4,4,6,5,6,7,4)
const DESTPIVOT=(6,4,4,6,5,6,7,4,7,7,5,5)
const APEXPIVOT=(5,6,7,4,7,7,5,5,6,4,4,6)
const OPPOPIVOT=(4,5,6,7,4,5,6,7,4,5,6,7)
const VER2EDGE=(0,1,2,3,3,5,1,5,4,0,4,2)
const EDGE2VER=(0,1,2,3,8,5)
const EPIVOT=(4,5,2,11,4,5,2,11,4,5,2,11)
const SNEXTPIVOT=(2,5,4,1,0,3)
const SORGPIVOT=(3,4,4,5,5,3)
const SDESTPIVOT=(4,3,5,4,3,5)
const SAPEXPIVOT=(5,5,3,3,4,4)
const ENEXTTBL=ntuple(i->(i-1+4)%12,12)
const EPREVTBL=ntuple(i->(i-1+8)%12,12)
const ENEXTESYMTBL=ntuple(i->ESYMTBL[ENEXTTBL[i]+1],12)
const EPREVESYMTBL=ntuple(i->ESYMTBL[EPREVTBL[i]+1],12)
const EORGOPPOTBL=ntuple(i->EPREVTBL[ESYMTBL[ENEXTTBL[i]+1]+1],12)
const EDESTOPPOTBL=ntuple(i->ENEXTTBL[ESYMTBL[EPREVTBL[i]+1]+1],12)
const BONDTBL=[(j&3)+(((i&12)+(j&12))%12) for i in 0:11,j in 0:11]
const FSYMTBL=[(j+12-(i&12))%12 for i in 0:11,j in 0:11]
const FACEPIVOT1=ntuple(i->ESYMTBL[i]&3,12)
const FACEPIVOT2=[FSYMTBL[ESYMTBL[i+1]+1,j+1] for i in 0:11,j in 0:11]
function _shtables(pivot::Bool)
    ts=zeros(Int,12,6);st=zeros(Int,12,6)
    for i in 0:11,j in 0:5
        even=(j&1)==0
        if even==pivot
            soffset=(6-((i&12)>>1))%6;toffset=(12-((j&6)<<1))%12
        else
            soffset=(i&12)>>1;toffset=(j&6)<<1
        end
        ts[i+1,j+1]=(j&1)+(((j&6)+soffset)%6)
        st[i+1,j+1]=(i&3)+(((i&12)+toffset)%12)
    end
    return ts,st
end
const TSBONDTBL,STBONDTBL=_shtables(true)
const TSPIVOTTBL,STPIVOTTBL=_shtables(false)

# Vertex types (enum verttype).
const UNUSEDVERTEX=0;const DUPLICATEDVERTEX=1;const RIDGEVERTEX=2
const FACETVERTEX=4;const VOLVERTEX=5;const FREESEGVERTEX=6
const FREEFACETVERTEX=7;const FREEVOLVERTEX=8;const DEADVERTEX=10

const DUMMY=Int32(-1)   # tetgenmesh::dummypoint

struct Tri
    tet::Int
    ver::Int
end
const TNULL=Tri(0,0)
struct Sh
    sh::Int
    shver::Int
end
const SNULL=Sh(0,0)
@inline tenc(t::Tri)=(t.tet<<4)|t.ver
@inline tdec(x::Int)=Tri(x>>4,x&15)
@inline senc(s::Sh)=(s.sh<<3)|s.shver
@inline sdec(x::Int)=Sh(x>>3,x&7)

# A memory pool of shell faces: eleven words per item (sh[0..10]), the shell
# mark and the flag word at shmarkindex+1.
mutable struct ShellPool
    w::Vector{Int}
    mark::Vector{Int32}
    flag::Vector{Int32}
    dead::Int        # top of the dead-item stack (0: empty)
    items::Int
    maxitems::Int
end
ShellPool()=ShellPool(Int[],Int32[],Int32[],0,0,0)

struct BadFace
    tt::Tri
    key::Float64
    cent::NTuple{6,Float64}
    forg::Int32
    fdest::Int32
    fapex::Int32
    foppo::Int32
end
BadFace(tt::Tri,forg,fdest,fapex)=BadFace(tt,0.0,ntuple(_->0.0,6),forg,fdest,fapex,Int32(0))

mutable struct Flip
    enqflag::Int
    chkencflag::Int
    unflip::Bool
    collectnewtets::Bool
    remove_ndelaunay_edge::Bool
    bak_tetprism_vol::Float64
    tetprism_vol_sum::Float64
    remove_large_angle::Bool
    cosdihed_in::Float64
    cosdihed_out::Float64
    checkflipeligibility::Bool
    seg::NTuple{2,Int32}
    fac::NTuple{3,Int32}
end
Flip()=Flip(0,0,false,false,false,0.0,0.0,false,0.0,0.0,false,(0,0),(0,0,0))
# A fresh flipconstraints, reusing `fc` (TetGen makes one on the stack per call).
function _flipreset!(fc::Flip)
    fc.enqflag=0;fc.chkencflag=0;fc.unflip=false;fc.collectnewtets=false
    fc.remove_ndelaunay_edge=false;fc.bak_tetprism_vol=0.0;fc.tetprism_vol_sum=0.0
    fc.remove_large_angle=false;fc.cosdihed_in=0.0;fc.cosdihed_out=0.0
    fc.checkflipeligibility=false;fc.seg=(0,0);fc.fac=(0,0,0)
    return fc
end

mutable struct Opt
    max_min_volume::Bool
    min_max_dihedangle::Bool
    initval::Float64
    imprval::Float64
    numofsearchdirs::Int
    searchstep::Float64
    maxiter::Int
    smthiter::Int
end
Opt()=Opt(false,false,0.0,0.0,10,0.01,-1,0)

mutable struct TG
    # Points: coordinates, pointmark, flag word (type<<8 | marks), point2tet.
    xyz::Vector{NTuple{3,Float64}}
    pmark::Vector{Int32}
    pflag::Vector{Int32}
    p2tet::Vector{Int}
    pdead::Int
    pitems::Int
    ninput::Int
    # Tetrahedra: neighbours, vertices (tet[4..7]), segment and subface
    # handles per edge/face, element marker.
    nb::Vector{Int}
    tv::Vector{Int32}
    tseg::Vector{Int}
    tsub::Vector{Int}
    tmark::Vector{Int32}
    tdead::Int
    titems::Int
    tmax::Int
    subfaces::ShellPool
    subsegs::ShellPool
    seed::UInt32
    nonconvex::Bool
    checksubsegflag::Bool
    checksubfaceflag::Bool
    autofliplinklevel::Int
    fliplinklevel::Int
    flipstarsize::Int
    epsilon::Float64
    cosmaxdihed::Float64
    cossmtdihed::Float64
    cosslidihed::Float64
    # Working lists (arraypools of tetgenmesh).
    flipstack::Vector{Tri}
    unflipqueue::Vector{BadFace}
    cavetetlist::Vector{Tri}
    cavebdrylist::Vector{Tri}
    caveoldtetlist::Vector{Tri}
    cavetetseglist::Vector{Sh}
    cavetetshlist::Vector{Sh}
    caveencseglist::Vector{Sh}
    caveencshlist::Vector{Sh}
    badtetrahedrons::Vector{Tri}
    subsegstack::Vector{Sh}
    subfacstack::Vector{Sh}
    subvertstack::Vector{Int32}
    cavetetvertlist::Vector{Int32}
    stars::Vector{Vector{Tri}}
    steiner::Int
    # Scratch handles for one flip at a time (flipnm, flipnm_post and
    # removefacebyflips use theirs only around their own flip calls).
    F3::Vector{Tri}
end

# ---------------------------------------------------------------------------
# Tetrahedron primitives.
@inline _nbi(t)=4*(t.tet-1)+(t.ver&3)+1
@inline _tvi(slot,pivot)=4*(slot-1)+pivot-3
@inline tvert(m,slot,pivot)=m.tv[_tvi(slot,pivot)]
@inline org(m,t)=m.tv[_tvi(t.tet,ORGPIVOT[t.ver+1])]
@inline dest(m,t)=m.tv[_tvi(t.tet,DESTPIVOT[t.ver+1])]
@inline apex(m,t)=m.tv[_tvi(t.tet,APEXPIVOT[t.ver+1])]
@inline oppo(m,t)=m.tv[_tvi(t.tet,OPPOPIVOT[t.ver+1])]
@inline setorg!(m,t,p)=(m.tv[_tvi(t.tet,ORGPIVOT[t.ver+1])]=p)
@inline setdest!(m,t,p)=(m.tv[_tvi(t.tet,DESTPIVOT[t.ver+1])]=p)
@inline setapex!(m,t,p)=(m.tv[_tvi(t.tet,APEXPIVOT[t.ver+1])]=p)
@inline setoppo!(m,t,p)=(m.tv[_tvi(t.tet,OPPOPIVOT[t.ver+1])]=p)
@inline function setvertices!(m,t,a,b,c,d)
    setorg!(m,t,a);setdest!(m,t,b);setapex!(m,t,c);setoppo!(m,t,d)
end
@inline enext(t)=Tri(t.tet,ENEXTTBL[t.ver+1])
@inline eprev(t)=Tri(t.tet,EPREVTBL[t.ver+1])
@inline esym(t)=Tri(t.tet,ESYMTBL[t.ver+1])
@inline enextesym(t)=Tri(t.tet,ENEXTESYMTBL[t.ver+1])
@inline eprevesym(t)=Tri(t.tet,EPREVESYMTBL[t.ver+1])
@inline eorgoppo(t)=Tri(t.tet,EORGOPPOTBL[t.ver+1])
@inline edestoppo(t)=Tri(t.tet,EDESTOPPOTBL[t.ver+1])
@inline function bond!(m,t1,t2)
    m.nb[_nbi(t1)]=(t2.tet<<4)|BONDTBL[t1.ver+1,t2.ver+1]
    m.nb[_nbi(t2)]=(t1.tet<<4)|BONDTBL[t2.ver+1,t1.ver+1]
    return nothing
end
@inline dissolve!(m,t)=(m.nb[_nbi(t)]=0)
@inline function fsym(m,t)
    x=m.nb[_nbi(t)]
    n=tdec(x)
    return Tri(n.tet,FSYMTBL[t.ver+1,n.ver+1])
end
@inline function fnext(m,t)
    x=m.nb[4*(t.tet-1)+FACEPIVOT1[t.ver+1]+1]
    n=tdec(x)
    return Tri(n.tet,FACEPIVOT2[t.ver+1,n.ver+1])
end
# decode(t.tet[i], neightet) for a raw face index i in 0..3.
@inline rawnb(m,slot,i)=tdec(m.nb[4*(slot-1)+i+1])
@inline ishulltet(m,t)=tvert(m,t.tet,7)==DUMMY
@inline isdeadtet(m,t)=t.tet==0 || tvert(m,t.tet,4)==0
@inline infect!(m,t)=(m.tmark[t.tet]|=Int32(1))
@inline uninfect!(m,t)=(m.tmark[t.tet]&=~Int32(1))
@inline infected(m,t)=(m.tmark[t.tet]&Int32(1))!=0
@inline marktest!(m,t)=(m.tmark[t.tet]|=Int32(2))
@inline unmarktest!(m,t)=(m.tmark[t.tet]&=~Int32(2))
@inline marktested(m,t)=(m.tmark[t.tet]&Int32(2))!=0
@inline markface!(m,t)=(m.tmark[t.tet]|=Int32(4)<<(t.ver&3))
@inline unmarkface!(m,t)=(m.tmark[t.tet]&=~(Int32(4)<<(t.ver&3)))
@inline facemarked(m,t)=(m.tmark[t.tet]&(Int32(4)<<(t.ver&3)))!=0
@inline marktest2!(m,t)=(m.tmark[t.tet]|=Int32(4096))
@inline unmarktest2!(m,t)=(m.tmark[t.tet]&=~Int32(4096))
@inline marktest2ed(m,t)=(m.tmark[t.tet]&Int32(4096))!=0
@inline elemcounter(m,t)=m.tmark[t.tet]>>16
@inline function setelemcounter!(m,t,value)
    c=m.tmark[t.tet]&Int32(65535)
    m.tmark[t.tet]=c|(Int32(value)<<16)
end
@inline increaseelemcounter!(m,t)=setelemcounter!(m,t,elemcounter(m,t)+1)
@inline decreaseelemcounter!(m,t)=setelemcounter!(m,t,elemcounter(m,t)-1)

# Subface/segment attachments of a tetrahedron.
@inline _subi(t)=4*(t.tet-1)+(t.ver&3)+1
@inline _segi(t)=6*(t.tet-1)+VER2EDGE[t.ver+1]+1
@inline issubface(m,t)=m.tsub[_subi(t)]!=0
@inline issubseg(m,t)=m.tseg[_segi(t)]!=0
@inline function tspivot(m,t)
    x=m.tsub[_subi(t)]
    x==0 && return SNULL
    s=sdec(x)
    return Sh(s.sh,TSPIVOTTBL[t.ver+1,s.shver+1])
end
@inline tsspivot1(m,t)=(x=m.tseg[_segi(t)];x==0 ? SNULL : sdec(x))
@inline function tsbond!(m,t,s)
    m.tsub[_subi(t)]=(s.sh<<3)|TSBONDTBL[t.ver+1,s.shver+1]
    m.subfaces.w[11*(s.sh-1)+10+(s.shver&1)]=(t.tet<<4)|STBONDTBL[t.ver+1,s.shver+1]
    return nothing
end
@inline tssbond1!(m,t,s)=(m.tseg[_segi(t)]=senc(s))
@inline tsdissolve!(m,t)=(m.tsub[_subi(t)]=0)
@inline function _cleartetattachments!(m,slot)
    for i in 1:4;m.tsub[4*(slot-1)+i]=0;end
    for i in 1:6;m.tseg[6*(slot-1)+i]=0;end
end

# Shell-face primitives on a pool `p`.
@inline sw(p,s,i)=p.w[11*(s.sh-1)+i+1]
@inline setsw!(p,s,i,x)=(p.w[11*(s.sh-1)+i+1]=x)
@inline sorg(p,s)=Int32(sw(p,s,SORGPIVOT[s.shver+1]))
@inline sdest(p,s)=Int32(sw(p,s,SDESTPIVOT[s.shver+1]))
@inline sapex(p,s)=Int32(sw(p,s,SAPEXPIVOT[s.shver+1]))
@inline function setshvertices!(p,s,a,b,c)
    setsw!(p,s,SORGPIVOT[s.shver+1],Int(a))
    setsw!(p,s,SDESTPIVOT[s.shver+1],Int(b))
    setsw!(p,s,SAPEXPIVOT[s.shver+1],Int(c))
end
@inline sesym(s)=Sh(s.sh,s.shver⊻1)
@inline senext(s)=Sh(s.sh,SNEXTPIVOT[s.shver+1])
@inline senext2(s)=Sh(s.sh,SNEXTPIVOT[SNEXTPIVOT[s.shver+1]+1])
@inline spivot(p,s)=(x=sw(p,s,s.shver>>1);x==0 ? SNULL : sdec(x))
@inline sbond1!(p,s1,s2)=setsw!(p,s1,s1.shver>>1,senc(s2))
@inline sdissolve!(p,s)=setsw!(p,s,s.shver>>1,0)
@inline sspivot(p,s)=(x=sw(p,s,6+(s.shver>>1));x==0 ? SNULL : sdec(x))
# ssbond: the subface `s` (in subfaces) and the segment `e` (in subsegs).
@inline function ssbond!(m,s,e)
    setsw!(m.subfaces,s,6+(s.shver>>1),senc(e))
    setsw!(m.subsegs,e,0,senc(s))
end
@inline function stpivot(m,s)
    x=sw(m.subfaces,s,9+(s.shver&1))
    x==0 && return TNULL
    t=tdec(x)
    return Tri(t.tet,STPIVOTTBL[t.ver+1,s.shver+1])
end
@inline sstbond1!(m,e,t)=setsw!(m.subsegs,e,9,tenc(t))
@inline sstpivot1(m,e)=(x=sw(m.subsegs,e,9);x==0 ? TNULL : tdec(x))
@inline sinfected(p,s)=(p.flag[s.sh]&Int32(1))!=0
@inline sinfect!(p,s)=(p.flag[s.sh]|=Int32(1))
@inline suninfect!(p,s)=(p.flag[s.sh]&=~Int32(1))
@inline smarktested(p,s)=(p.flag[s.sh]&Int32(2))!=0
@inline isshsubseg(m,s)=sw(m.subfaces,s,6+(s.shver>>1))!=0
@inline shdead(p,s)=sw(p,s,3)==0

function makeshellface!(p::ShellPool)
    if p.dead!=0
        slot=p.dead
        p.dead=tdecsh(p.w[11*(slot-1)+1])
    else
        p.maxitems+=1;slot=p.maxitems
        resize!(p.w,11*slot);push!(p.mark,0);push!(p.flag,0)
    end
    for i in 1:11;p.w[11*(slot-1)+i]=0;end
    p.mark[slot]=0;p.flag[slot]=0
    p.items+=1
    return Sh(slot,0)
end
# A dead item keeps the next dead slot in its first word, as memorypool does.
@inline tdecsh(x::Int)=x>>3
function shellfacedealloc!(p::ShellPool,slot::Int)
    p.w[11*(slot-1)+4]=0                 # dyingsh[3] = NULL
    p.w[11*(slot-1)+1]=p.dead<<3         # memorypool link in word 0
    p.dead=slot
    p.items-=1
    return nothing
end

# ---------------------------------------------------------------------------
# Pools of tetrahedra and points.
function maketetrahedron!(m)
    if m.tdead!=0
        slot=m.tdead
        m.tdead=m.nb[4*(slot-1)+1]>>4
    else
        m.tmax+=1;slot=m.tmax
        append!(m.nb,(0,0,0,0));append!(m.tv,(Int32(0),Int32(0),Int32(0),Int32(0)))
        append!(m.tseg,(0,0,0,0,0,0));append!(m.tsub,(0,0,0,0));push!(m.tmark,Int32(0))
    end
    for i in 1:4
        m.nb[4*(slot-1)+i]=0;m.tv[4*(slot-1)+i]=0
    end
    _cleartetattachments!(m,slot)
    m.tmark[slot]=0
    m.titems+=1
    return Tri(slot,11)
end
function tetrahedrondealloc!(m,slot::Int)
    m.tv[4*(slot-1)+1]=0                 # tet[4] = NULL
    m.nb[4*(slot-1)+1]=m.tdead<<4        # memorypool link in word 0
    m.tdead=slot
    m.titems-=1
    return nothing
end
@inline pointtype(m,p)=m.pflag[p]>>8
@inline setpointtype!(m,p,v)=(m.pflag[p]=(Int32(v)<<8)+(m.pflag[p]&Int32(255)))
function makepoint!(m,vtype)
    if m.pdead!=0
        p=m.pdead
        m.pdead=m.p2tet[p]               # the dead link (see pointdealloc!)
    else
        push!(m.xyz,(0.0,0.0,0.0));push!(m.pmark,0);push!(m.pflag,0);push!(m.p2tet,0)
        p=length(m.xyz)
    end
    m.pitems+=1
    m.p2tet[p]=0
    m.pmark[p]=Int32(m.pitems)           # firstnumber = 1
    m.pflag[p]=0
    setpointtype!(m,p,vtype)
    return Int32(p)
end
function pointdealloc!(m,p)
    setpointtype!(m,p,DEADVERTEX)
    m.p2tet[p]=m.pdead
    m.pdead=Int(p)
    m.pitems-=1
    return nothing
end
@inline setpoint2tet!(m,p,t::Tri)=(p==DUMMY || (m.p2tet[p]=tenc(t)))
function point2tetorg(m,p)
    t=tdec(m.p2tet[p])
    t.tet==0 && throw(Unsupported("point2tet of a point without a tetrahedron"))
    s=t.tet
    tvert(m,s,4)==p && return Tri(s,11)
    tvert(m,s,5)==p && return Tri(s,3)
    tvert(m,s,6)==p && return Tri(s,7)
    return Tri(s,0)
end

# ---------------------------------------------------------------------------
# Predicates on point ids.
@inline _xyz(m,p)=m.xyz[p]
@inline orient3d(m,a,b,c,d)=Float64(orient3(m.xyz[a],m.xyz[b],m.xyz[c],m.xyz[d]))
@inline function orient3d(m,a,b,c,d::NTuple{3,Float64})
    Float64(orient3(m.xyz[a],m.xyz[b],m.xyz[c],d))
end
function orient3dfast(pa,pb,pc,pd)
    adx=pa[1]-pd[1];bdx=pb[1]-pd[1];cdx=pc[1]-pd[1]
    ady=pa[2]-pd[2];bdy=pb[2]-pd[2];cdy=pc[2]-pd[2]
    adz=pa[3]-pd[3];bdz=pb[3]-pd[3];cdz=pc[3]-pd[3]
    return adx*(bdy*cdz-bdz*cdy)+bdx*(cdy*adz-cdz*ady)+cdx*(ady*bdz-adz*bdy)
end
function insphere_s(m,a,b,c,d,e)
    sign=insphere(m.xyz[a],m.xyz[b],m.xyz[c],m.xyz[d],m.xyz[e])
    sign!=0 && return Float64(sign)
    pt=[a,b,c,d,e];swaps=0;n=5
    while true
        count=0;n-=1
        for i in 1:n
            if m.pmark[pt[i]]>m.pmark[pt[i+1]]
                pt[i],pt[i+1]=pt[i+1],pt[i];count+=1
            end
        end
        swaps+=count
        count>0 || break
    end
    oriA=orient3d(m,pt[2],pt[3],pt[4],pt[5])
    if oriA!=0
        isodd(swaps) && (oriA=-oriA)
        return oriA
    end
    oriB=-orient3d(m,pt[1],pt[3],pt[4],pt[5])
    oriB==0 && throw(Unsupported("insphere_s: unresolved symbolic perturbation"))
    isodd(swaps) && (oriB=-oriB)
    return oriB
end
@inline function distance(p,q)
    return sqrt((q[1]-p[1])*(q[1]-p[1])+(q[2]-p[2])*(q[2]-p[2])+(q[3]-p[3])*(q[3]-p[3]))
end
@inline dot3(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

function randomnation!(m,choices::Integer)
    choices=UInt32(choices)
    if choices>=UInt32(714025)
        newrandom=(m.seed*UInt32(1366)+UInt32(150889))%UInt32(714025)
        m.seed=(newrandom*UInt32(1366)+UInt32(150889))%UInt32(714025)
        newrandom=newrandom*(choices÷UInt32(714025))+m.seed
        return Int(newrandom>=choices ? newrandom-choices : newrandom)
    end
    m.seed=(m.seed*UInt32(1366)+UInt32(150889))%UInt32(714025)
    return Int(m.seed%choices)
end

include("ModelNativeTetgenBRGeometry.jl")
include("ModelNativeTetgenBRFlips.jl")
include("ModelNativeTetgenBRRecover.jl")
include("ModelNativeTetgenBRMesh.jl")
include("ModelNativeTetgenBROptimize.jl")
end
