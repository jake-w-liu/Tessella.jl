# Hierarchical H1 and H(curl) reference bases matching Gmsh 4.15.2's
# HierarchicalBasisH1*/HierarchicalBasisHcurl* implementations, including the
# Solin edge/face conventions, orientation sign tables, and key metadata.
#
# The orthogonal polynomials below are verbatim Horner-form ports of
# OrthogonalPoly.cpp so results track the pinned release's rounding closely.

function _gmsh_lobatto(order::Int,x::Float64)
    xsquare=x*x
    order==0 && return 0.5*(1.0-x)
    order==1 && return 0.5*(1.0+x)
    if order==2
        L=(-1.0+xsquare)
        return L*0.5*sqrt(3.0/2.0)
    elseif order==3
        L=x*(-1.0+xsquare)
        return L*0.5*sqrt(5.0/2.0)
    elseif order==4
        L=1.0+xsquare*(-6.0+5.0*xsquare)
        return L*1.0/8.0*sqrt(7.0/2.0)
    elseif order==5
        L=x*(3.0+xsquare*(-10.0+7.0*xsquare))
        return L*3.0/8.0*2.0^-0.5
    elseif order==6
        L=-1.0+xsquare*(15.0+xsquare*(-35.0+21.0*xsquare))
        return L*1.0/16.0*sqrt(11.0/2.0)
    elseif order==7
        L=x*(-5.0+xsquare*(35.0+xsquare*(-63.0+33.0*xsquare)))
        return L*1.0/16.0*sqrt(13.0/2.0)
    elseif order==8
        L=5.0+xsquare*(-140.0+xsquare*(630.0+xsquare*(-924.0+429.0*xsquare)))
        return L*1.0/128.0*sqrt(15.0/2.0)
    elseif order==9
        L=x*(35.0+xsquare*(-420.0+xsquare*(1386.0+xsquare*(-1716.0+715.0*xsquare))))
        return L*1.0/128.0*sqrt(17.0/2.0)
    elseif order==10
        L=-7.0+xsquare*(315.0+xsquare*(-2310.0+xsquare*(6006.0+xsquare*(-6435.0+2431.0*xsquare))))
        return L*1.0/256.0*sqrt(19.0/2.0)
    elseif order==11
        L=x*(-63.0+xsquare*(1155.0+xsquare*(-6006.0+xsquare*(12870.0+xsquare*(-12155.0+4199.0*xsquare)))))
        return L*1.0/256.0*sqrt(21.0/2.0)
    elseif order==12
        L=21.0+xsquare*(-1386.0+xsquare*(15015.0+xsquare*(-60060.0+xsquare*(109395.0+xsquare*(-92378.0+29393.0*xsquare)))))
        return L*1.0/1024.0*sqrt(23.0/2.0)
    elseif order==13
        L=x*(231.0+xsquare*(-6006.0+xsquare*(45045.0+xsquare*(-145860.0+xsquare*(230945.0+xsquare*(-176358.0+52003.0*xsquare))))))
        return L*5.0/1024.0*2.0^-0.5
    elseif order==14
        L=-33.0+xsquare*(3003.0+xsquare*(-45045.0+xsquare*(255255.0+xsquare*(-692835.0+xsquare*(969969.0+xsquare*(-676039.0+185725.0*xsquare))))))
        return L*3.0/2048.0*sqrt(3.0/2.0)
    elseif order==15
        L=x*(-429.0+xsquare*(15015.0+xsquare*(-153153.0+xsquare*(692835.0+xsquare*(-1616615.0+xsquare*(2028117.0+xsquare*(-1300075.0+334305.0*xsquare)))))))
        return L*1.0/2048.0*sqrt(29.0/2.0)
    end
    throw(ArgumentError("Lobatto functions are written for orders =< 15"))
end

function _gmsh_dlobatto(order::Int,x::Float64)
    xsquare=x*x
    order==0 && return -0.5
    order==1 && return 0.5
    if order==2
        dL=2.0*x
        return dL*0.5*sqrt(3.0/2.0)
    elseif order==3
        dL=-1.0+3.0*xsquare
        return dL*0.5*sqrt(5.0/2.0)
    elseif order==4
        dL=x*(-12.0+20.0*xsquare)
        return dL*1.0/8.0*sqrt(7.0/2.0)
    elseif order==5
        dL=3.0+xsquare*(-30.0+35.0*xsquare)
        return dL*3.0/8.0*2.0^-0.5
    elseif order==6
        dL=x*(30.0+xsquare*(-140.0+126.0*xsquare))
        return dL*1.0/16.0*sqrt(11.0/2.0)
    elseif order==7
        dL=-5.0+xsquare*(105.0+xsquare*(-315.0+231.0*xsquare))
        return dL*1.0/16.0*sqrt(13.0/2.0)
    elseif order==8
        dL=x*(-280.0+xsquare*(2520.0+xsquare*(-5544.0+3432.0*xsquare)))
        return dL*1.0/128.0*sqrt(15.0/2.0)
    elseif order==9
        dL=35.0+xsquare*(-1260.0+xsquare*(6930.0+xsquare*(-12012.0+6435.0*xsquare)))
        return dL*1.0/128.0*sqrt(17.0/2.0)
    elseif order==10
        dL=x*(630.0+xsquare*(-9240.0+xsquare*(36036.0+xsquare*(-51480.0+24310.0*xsquare))))
        return dL*1.0/256.0*sqrt(19.0/2.0)
    elseif order==11
        dL=-63.0+xsquare*(3465.0+xsquare*(-30030.0+xsquare*(90090.0+xsquare*(-109395.0+46189.0*xsquare))))
        return dL*1.0/256.0*sqrt(21.0/2.0)
    elseif order==12
        dL=x*(-2772.0+xsquare*(60060.0+xsquare*(-360360.0+xsquare*(875160.0+xsquare*(-923780.0+352716.0*xsquare)))))
        return dL*1.0/1024.0*sqrt(23.0/2.0)
    elseif order==13
        dL=231.0+xsquare*(-18018.0+xsquare*(225225.0+xsquare*(-1021020.0+xsquare*(2078505.0+xsquare*(-1939938.0+676039.0*xsquare)))))
        return dL*5.0/1024.0*2.0^-0.5
    elseif order==14
        dL=x*(6006.0+xsquare*(-180180.0+xsquare*(1531530.0+xsquare*(-5542680.0+xsquare*(9699690.0+xsquare*(-8112468.0+2600150.0*xsquare))))))
        return dL*3.0/2048.0*sqrt(3.0/2.0)
    elseif order==15
        dL=-429.0+xsquare*(45045.0+xsquare*(-765765.0+xsquare*(4849845.0+xsquare*(-14549535.0+xsquare*(22309287.0+xsquare*(-16900975.0+5014575.0*xsquare))))))
        return dL*1.0/2048.0*sqrt(29.0/2.0)
    end
    throw(ArgumentError("Lobatto functions are written for orders =< 15"))
end

function _gmsh_kernel(order::Int,x::Float64)
    xsquare=x*x
    order==0 && return -sqrt(6.0)
    order==1 && return -x*sqrt(10.0)
    if order==2
        phi=1.0-5.0*xsquare
        return phi*0.5*sqrt(7.0/2.0)
    elseif order==3
        phi=x*(3.0-7.0*xsquare)
        return phi*3.0/2.0*2.0^-0.5
    elseif order==4
        phi=-1.0+xsquare*(14.0-21.0*xsquare)
        return phi*1.0/4.0*sqrt(11.0/2.0)
    elseif order==5
        phi=x*(-5.0+xsquare*(30.0-33.0*xsquare))
        return phi*1.0/4.0*sqrt(13.0/2.0)
    elseif order==6
        phi=5.0+xsquare*(-135.0+xsquare*(495.0-429.0*xsquare))
        return phi*1.0/32.0*sqrt(15.0/2.0)
    elseif order==7
        phi=x*(35.0+xsquare*(-385.0+xsquare*(1001.0-715.0*xsquare)))
        return phi*1.0/32.0*sqrt(17.0/2.0)
    elseif order==8
        phi=-7.0+xsquare*(308.0+xsquare*(-2002.0+xsquare*(4004.0-2431.0*xsquare)))
        return phi*1.0/64.0*sqrt(19.0/2.0)
    elseif order==9
        phi=x*(-63.0+xsquare*(1092.0+xsquare*(-4914.0+xsquare*(7956.0-4199.0*xsquare))))
        return phi*1.0/64.0*sqrt(21.0/2.0)
    elseif order==10
        phi=21.0+xsquare*(-1365.0+xsquare*(13650.0+xsquare*(-46410.0+xsquare*(62985.0-29393.0*xsquare))))
        return phi*1.0/256.0*sqrt(23.0/2.0)
    elseif order==11
        phi=x*(231.0+xsquare*(-5775.0+xsquare*(39270.0+xsquare*(-106590.0+xsquare*(124355.0-52003.0*xsquare)))))
        return phi*5.0/256.0*2.0^-0.5
    elseif order==12
        phi=-33.0+xsquare*(2970.0+xsquare*(-42075.0+xsquare*(213180.0+xsquare*(-479655.0+xsquare*(490314.0-185725.0*xsquare)))))
        return phi*3.0/512.0*sqrt(3.0/2.0)
    elseif order==13
        phi=x*(-429.0+xsquare*(14586.0+xsquare*(-138567.0+xsquare*(554268.0+xsquare*(-1062347.0+xsquare*(965770.0-334305.0*xsquare))))))
        return phi*1.0/512.0*sqrt(29.0/2.0)
    end
    throw(ArgumentError("Lobatto functions are written for orders =< 15"))
end

function _gmsh_dkernel(order::Int,x::Float64)
    xsquare=x*x
    order==0 && return 0.0
    order==1 && return -sqrt(10.0)
    if order==2
        dphi=-10.0*x
        return dphi*0.5*sqrt(7.0/2.0)
    elseif order==3
        dphi=3.0-21.0*xsquare
        return dphi*3.0/2.0*2.0^-0.5
    elseif order==4
        dphi=x*(28.0-84.0*xsquare)
        return dphi*1.0/4.0*sqrt(11.0/2.0)
    elseif order==5
        dphi=-5.0+xsquare*(90.0-165.0*xsquare)
        return dphi*1.0/4.0*sqrt(13.0/2.0)
    elseif order==6
        dphi=x*(-270.0+xsquare*(1980.0-2574.0*xsquare))
        return dphi*1.0/32.0*sqrt(15.0/2.0)
    elseif order==7
        dphi=35.0+xsquare*(-1155.0+xsquare*(5005.0-5005.0*xsquare))
        return dphi*1.0/32.0*sqrt(17.0/2.0)
    elseif order==8
        dphi=x*(616.0+xsquare*(-8008.0+xsquare*(24024.0-19448.0*xsquare)))
        return dphi*1.0/64.0*sqrt(19.0/2.0)
    elseif order==9
        dphi=-63.0+xsquare*(3276.0+xsquare*(-24570.0+xsquare*(55692.0-37791.0*xsquare)))
        return dphi*1.0/64.0*sqrt(21.0/2.0)
    elseif order==10
        dphi=x*(-2730.0+xsquare*(54600.0+xsquare*(-278460.0+xsquare*(503880.0-293930.0*xsquare))))
        return dphi*1.0/256.0*sqrt(23.0/2.0)
    elseif order==11
        dphi=231.0+xsquare*(-17325.0+xsquare*(196350.0+xsquare*(-746130.0+xsquare*(1119195.0-572033.0*xsquare))))
        return dphi*5.0/256.0*2.0^-0.5
    elseif order==12
        dphi=x*(5940.0+xsquare*(-168300.0+xsquare*(1279080.0+xsquare*(-3837240.0+xsquare*(4903140.0-2228700.0*xsquare)))))
        return dphi*3.0/512.0*sqrt(3.0/2.0)
    elseif order==13
        dphi=-429.0+xsquare*(43758.0+xsquare*(-692835.0+xsquare*(3879876.0+xsquare*(-9561123.0+xsquare*(10623470.0-4345965.0*xsquare)))))
        return dphi*1.0/512.0*sqrt(29.0/2.0)
    end
    throw(ArgumentError("Lobatto functions are written for orders =< 15"))
end

function _gmsh_legendre(order::Int,x::Float64)
    xsquare=x*x
    order==0 && return 1.0
    order==1 && return x
    order==2 && return 1.5*xsquare-0.5
    order==3 && return 0.5*x*(5.0*xsquare-3.0)
    if order==4
        L=3.0+xsquare*(35.0*xsquare-30.0)
        return 1.0/8.0*L
    elseif order==5
        L=x*(xsquare*(63.0*xsquare-70.0)+15.0)
        return 1.0/8.0*L
    elseif order==6
        L=((231.0*xsquare-315.0)*xsquare+105.0)*xsquare-5.0
        return 1.0/16.0*L
    elseif order==7
        L=x*(((429.0*xsquare-693.0)*xsquare+315.0)*xsquare-35.0)
        return 1.0/16.0*L
    elseif order==8
        L=(((6435.0*xsquare-12012.0)*xsquare+6930.0)*xsquare-1260.0)*xsquare+35.0
        return 1.0/128.0*L
    elseif order==9
        L=((((12155.0*xsquare-25740.0)*xsquare+18018.0)*xsquare-4620.0)*xsquare+315.0)*x
        return 1.0/128.0*L
    elseif order==10
        L=((((46189.0*xsquare-109395.0)*xsquare+90090.0)*xsquare-30030.0)*xsquare+3465.0)*xsquare-63.0
        return 1.0/256.0*L
    end
    throw(ArgumentError("Legendre functions are written for orders =< 10"))
end

function _gmsh_dlegendre(order::Int,x::Float64)
    xsquare=x*x
    order==0 && return 0.0
    order==1 && return 1.0
    order==2 && return 3.0*x
    order==3 && return 0.5*(15.0*xsquare-3.0)
    if order==4
        dL=x*(140.0*xsquare-60.0)
        return 1.0/8.0*dL
    elseif order==5
        dL=15.0+xsquare*(315.0*xsquare-210.0)
        return 1.0/8.0*dL
    elseif order==6
        dL=x*(210.0+xsquare*(1386.0*xsquare-1260.0))
        return 1.0/16.0*dL
    elseif order==7
        dL=((xsquare*3003.0-3465.0)*xsquare+945.0)*xsquare-35.0
        return 1.0/16.0*dL
    elseif order==8
        dL=x*(((51480.0*xsquare-72072.0)*xsquare+27720.0)*xsquare-2520.0)
        return 1.0/128.0*dL
    elseif order==9
        dL=315.0+xsquare*(-13860.0+xsquare*(90090.0+xsquare*(-180180.0+109395.0*xsquare)))
        return 1.0/128.0*dL
    elseif order==10
        dL=x*(6930.0+xsquare*(-120120.0+xsquare*(540540.0+xsquare*(-875160.0+461890.0*xsquare))))
        return 1.0/256.0*dL
    end
    throw(ArgumentError("Legendre functions are written for orders =< 10"))
end

# ---------------------------------------------------------------------------
# Forward-mode value+gradient factor. Every hierarchical basis function is a
# product of linear barycentric factors and kernel/Lobatto factors of linear
# arguments, so a tiny product rule reproduces Gmsh's hand-coded gradients.
# ---------------------------------------------------------------------------
struct _BF
    v::Float64
    g::NTuple{3,Float64}
end
@inline _BF(v::Float64)=_BF(v,(0.0,0.0,0.0))
@inline function Base.:(*)(a::_BF,b::_BF)
    return _BF(a.v*b.v,(a.g[1]*b.v+a.v*b.g[1],
                        a.g[2]*b.v+a.v*b.g[2],
                        a.g[3]*b.v+a.v*b.g[3]))
end
@inline Base.:(-)(a::_BF,b::_BF)=_BF(a.v-b.v,(a.g[1]-b.g[1],
                                            a.g[2]-b.g[2],
                                            a.g[3]-b.g[3]))
@inline Base.:(*)(s::Float64,a::_BF)=_BF(s*a.v,(s*a.g[1],s*a.g[2],s*a.g[3]))
@inline Base.:(*)(s::Integer,a::_BF)=Float64(s)*a
@inline Base.:(*)(a::_BF,s::Integer)=Float64(s)*a
@inline _ker(k::Int,x::_BF)=
    _BF(_gmsh_kernel(k,x.v),
        (_gmsh_dkernel(k,x.v)*x.g[1],
         _gmsh_dkernel(k,x.v)*x.g[2],
         _gmsh_dkernel(k,x.v)*x.g[3]))
@inline _lob(k::Int,x::_BF)=
    _BF(_gmsh_lobatto(k,x.v),
        (_gmsh_dlobatto(k,x.v)*x.g[1],
         _gmsh_dlobatto(k,x.v)*x.g[2],
         _gmsh_dlobatto(k,x.v)*x.g[3]))

# Solin vertex tables (0-based) per reference family.
const _SOLIN_EDGES=Dict{Symbol,NTuple}(
    :lin=>((0,1),),
    :tri=>((0,1),(1,2),(2,0)),
    :qua=>((0,1),(1,2),(3,2),(0,3)),
    :tet=>((0,1),(1,2),(2,0),(0,3),(2,3),(1,3)),
    :hex=>((0,1),(0,3),(0,4),(1,2),(1,5),(3,2),
           (2,6),(3,7),(4,5),(4,7),(5,6),(7,6)),
    :pri=>((0,1),(0,2),(0,3),(1,2),(1,4),
           (2,5),(3,4),(3,5),(4,5)))
const _SOLIN_FACES=Dict{Symbol,Tuple{Vararg{NTuple{N,Int} where N}}}(
    :tri=>((0,1,2),),
    :qua=>((0,1,3,2),),
    :tet=>((0,1,2),(0,1,3),(0,2,3),(1,2,3)),
    :hex=>((0,1,3,2),(0,1,4,5),(0,3,4,7),
           (1,2,5,6),(3,2,7,6),(4,5,7,6)),
    :pri=>((0,1,3,4),(0,2,3,5),(1,2,4,5),
           (0,1,2),(3,4,5)))

# Per-edge function count for a uniform-order hierarchical basis.
@inline _h1_edge_funcs(::Int,p::Int)=p-1
@inline _hcurl_edge_funcs(::Int,p::Int)=p+1

# (offset, count, :tri|:quad) blocks inside the concatenated face table.
function _h1_face_blocks(family::Symbol,p::Int)
    quad_count=max(p-1,0)*max(p-1,0)
    tri_count=max(p-1,0)*max(p-2,0)÷2
    blocks=Tuple{Int,Int,Symbol}[]
    offset=0
    faces=get(_SOLIN_FACES,family,())
    for face in faces
        count=length(face)==3 ? tri_count : quad_count
        push!(blocks,(offset,count,length(face)==3 ? :tri : :quad))
        offset+=count
    end
    return blocks
end

function _h1_counts(family::Symbol,p::Int)
    nv=family===:pnt ? 1 :
       family===:lin ? 2 :
       family===:tri ? 3 :
       (family===:tet || family===:qua) ? 4 :
       family===:hex ? 8 : 6
    ne=length(get(_SOLIN_EDGES,family,()))
    nefuncs=ne*max(p-1,0)
    blocks=_h1_face_blocks(family,p)
    nffuncs=sum(b[2] for b in blocks;init=0)
    nbubble=family===:tet ? max((p-1)*(p-2)*(p-3),0)÷6 :
            family===:hex ? max(p-1,0)^3 :
            family===:pri ? max((p-1)*(p-2)*(p-1),0)÷2 : 0
    return (vertex=nv,edge=nefuncs,face=nffuncs,bubble=nbubble)
end

# ---------------------------------------------------------------------------
# H1 canonical (orientation-0) generators. Gmsh's reference coordinates:
#   lin [-1,1]; tri/tet [0,1] simplices; qua/hex [-1,1]^2/^3; pri [0,1]^2x[-1,1]
# Each returns (vertex, edge, face, bubble) _BF vectors in Gmsh's function
# order. p>=1 for every family except :pnt, which ignores the order.
# ---------------------------------------------------------------------------
function _h1_tables(::Val{:pnt},p::Int,u::Float64,v::Float64,w::Float64)
    return _BF[_BF(1.0)],_BF[],_BF[],_BF[]
end

function _h1_tables(::Val{:lin},p::Int,u::Float64,v::Float64,w::Float64)
    l1=_BF(0.5*(1.0+u),(0.5,0.0,0.0))
    l2=_BF(0.5*(1.0-u),(-0.5,0.0,0.0))
    e=_BF[l1*l2*_ker(k,l1-l2) for k in 0:p-2]
    return _BF[l2,l1],e,_BF[],_BF[]
end

function _h1_tables(::Val{:tri},p::Int,u::Float64,v::Float64,w::Float64)
    l1=_BF(v,(0.0,1.0,0.0))
    l2=_BF(1.0-u-v,(-1.0,-1.0,0.0))
    l3=_BF(u,(1.0,0.0,0.0))
    vb=_BF[l2,l3,l1]
    e=_BF[l3*l2*_ker(k,l3-l2) for k in 0:p-2]
    append!(e,_BF[l1*l3*_ker(k,l1-l3) for k in 0:p-2])
    append!(e,_BF[l2*l1*_ker(k,l2-l1) for k in 0:p-2])
    product=l1*l2*l3
    f=_BF[]
    sizehint!(f,(p-1)*(p-2)÷2)
    for n1 in 0:p-3
        k1=_ker(n1,l3-l2)
        for n2 in 0:p-3-n1
            push!(f,product*k1*_ker(n2,l2-l1))
        end
    end
    return vb,e,f,_BF[]
end

function _h1_tables(::Val{:tet},p::Int,u::Float64,v::Float64,w::Float64)
    l1=_BF(v,(0.0,1.0,0.0))
    l2=_BF(1.0-u-v-w,(-1.0,-1.0,-1.0))
    l3=_BF(u,(1.0,0.0,0.0))
    l4=_BF(w,(0.0,0.0,1.0))
    vb=_BF[l2,l3,l1,l4]
    subs=(l3-l2,l1-l3,l1-l2,l4-l2,l4-l1,l4-l3)
    eprods=(l3*l2,l1*l3,l2*l1,l4*l2,l4*l1,l4*l3)
    e=_BF[]
    sizehint!(e,6*(p-1))
    for edge in 1:6
        for k in 0:p-2
            sign=(edge==3 && isodd(k)) ? -1.0 : 1.0
            push!(e,sign*eprods[edge]*_ker(k,subs[edge]))
        end
    end
    faceprods=(eprods[1]*l1,eprods[1]*l4,eprods[3]*l4,eprods[2]*l4)
    targets=((1,3),(1,4),(3,4),(2,6))
    f=_BF[]
    sizehint!(f,2*(p-1)*(p-2))
    for face in 1:4
        t1,t2=targets[face]
        for n1 in 0:p-3
            k1=_ker(n1,subs[t1])
            for n2 in 0:p-3-n1
                sign=isodd(n2) ? -1.0 : 1.0
                push!(f,sign*faceprods[face]*k1*_ker(n2,subs[t2]))
            end
        end
    end
    product=l1*l2*l3*l4
    b=_BF[]
    sizehint!(b,max(p-1,0)*max(p-2,0)*max(p-3,0)÷6)
    for n1 in 0:p-4
        k1=_ker(n1,subs[3])
        for n2 in 0:p-4-n1
            k2=_ker(n2,subs[1])
            for n3 in 0:p-4-n1-n2
                push!(b,product*k1*k2*_ker(n3,subs[4]))
            end
        end
    end
    return vb,e,f,b
end

function _h1_tables(::Val{:qua},p::Int,u::Float64,v::Float64,w::Float64)
    l1=_BF(0.5*(1.0+u),(0.5,0.0,0.0))
    l2=_BF(0.5*(1.0-u),(-0.5,0.0,0.0))
    l3=_BF(0.5*(1.0+v),(0.0,0.5,0.0))
    l4=_BF(0.5*(1.0-v),(0.0,-0.5,0.0))
    ub=_BF(u,(1.0,0.0,0.0))
    vbb=_BF(v,(0.0,1.0,0.0))
    vb=_BF[l2*l4,l1*l4,l1*l3,l2*l3]
    e=_BF[]
    sizehint!(e,4*(p-1))
    for k in 2:p
        push!(e,l4*_lob(k,ub))
    end
    for k in 2:p
        push!(e,l1*_lob(k,vbb))
    end
    for k in 2:p
        push!(e,l3*_lob(k,ub))
    end
    for k in 2:p
        push!(e,l2*_lob(k,vbb))
    end
    f=_BF[]
    sizehint!(f,(p-1)*(p-1))
    for n1 in 2:p
        k1=_lob(n1,ub)
        for n2 in 2:p
            push!(f,k1*_lob(n2,vbb))
        end
    end
    return vb,e,f,_BF[]
end

function _h1_tables(::Val{:hex},p::Int,u::Float64,v::Float64,w::Float64)
    l=_BF[_BF(0.5*(1.0+u),(0.5,0.0,0.0)),
          _BF(0.5*(1.0-u),(-0.5,0.0,0.0)),
          _BF(0.5*(1.0+v),(0.0,0.5,0.0)),
          _BF(0.5*(1.0-v),(0.0,-0.5,0.0)),
          _BF(0.5*(1.0+w),(0.0,0.0,0.5)),
          _BF(0.5*(1.0-w),(0.0,0.0,-0.5))]
    product=_BF[l[4]*l[6],l[2]*l[6],l[2]*l[4],l[1]*l[6],
                l[4]*l[1],l[3]*l[6],l[3]*l[1],l[3]*l[2],
                l[4]*l[5],l[5]*l[2],l[5]*l[1],l[5]*l[3]]
    vb=_BF[l[2]*product[1],l[1]*product[1],l[1]*product[6],
           l[2]*product[6],l[2]*product[9],l[1]*product[9],
           l[1]*product[12],l[2]*product[12]]
    ub=_BF(u,(1.0,0.0,0.0))
    vbb=_BF(v,(0.0,1.0,0.0))
    wb=_BF(w,(0.0,0.0,1.0))
    lku=_BF[_lob(k,ub) for k in 2:p]
    lkv=_BF[_lob(k,vbb) for k in 2:p]
    lkw=_BF[_lob(k,wb) for k in 2:p]
    dirs=(lku,lkv,lkw)
    edge_dirs=(1,2,3,2,3,1,3,3,1,2,2,1)
    e=_BF[]
    sizehint!(e,12*(p-1))
    for edge in 1:12
        lk=dirs[edge_dirs[edge]]
        for k in 1:p-1
            push!(e,lk[k]*product[edge])
        end
    end
    face_lambda=(6,4,2,1,3,5)
    face_dirs=((1,2),(1,3),(2,3),(2,3),(1,3),(1,2))
    f=_BF[]
    sizehint!(f,6*(p-1)*(p-1))
    for face in 1:6
        lk1,lk2=dirs[face_dirs[face][1]],dirs[face_dirs[face][2]]
        transverse=l[face_lambda[face]]
        for n1 in 1:p-1
            for n2 in 1:p-1
                push!(f,transverse*lk1[n1]*lk2[n2])
            end
        end
    end
    b=_BF[]
    sizehint!(b,(p-1)^3)
    for n1 in 1:p-1
        for n2 in 1:p-1
            for n3 in 1:p-1
                push!(b,lku[n1]*lkv[n2]*lkw[n3])
            end
        end
    end
    return vb,e,f,b
end

function _h1_tables(::Val{:pri},p::Int,u::Float64,v::Float64,w::Float64)
    l1=_BF(v,(0.0,1.0,0.0))
    l2=_BF(1.0-u-v,(-1.0,-1.0,0.0))
    l3=_BF(u,(1.0,0.0,0.0))
    l4=_BF(0.5*(1.0+w),(0.0,0.0,0.5))
    l5=_BF(0.5*(1.0-w),(0.0,0.0,-0.5))
    vb=_BF[l2*l5,l3*l5,l1*l5,l2*l4,l4*l3,l1*l4]
    product=_BF[vb[1]*l3,vb[1]*l1,vb[1]*l4,vb[2]*l1,vb[2]*l4,
                vb[3]*l4,vb[4]*l3,vb[4]*l1,vb[5]*l1]
    subs=(l3-l2,l1-l2,l4-l5,l1-l3,l2-l1)
    wb=_BF(w,(0.0,0.0,1.0))
    edge_phi=(1,2,3,4,3,3,1,2,4)
    e=_BF[]
    sizehint!(e,9*(p-1))
    for edge in 1:9
        arg=subs[edge_phi[edge]]
        for k in 0:p-2
            push!(e,product[edge]*_ker(k,arg))
        end
    end
    # Quadrilateral faces (Solin faces 0-2) come first, then triangular
    # faces (Solin faces 3-4), matching Gmsh's packed face table.
    qlam=(product[1]*l4,product[8]*l5,product[9]*l5)
    qphis=((1,3),(2,3),(4,3))
    f=_BF[]
    sizehint!(f,3*(p-1)*(p-1)+(p-1)*(p-2))
    for face in 1:3
        i1,i2=qphis[face]
        for n1 in 0:p-2
            k1=_ker(n1,subs[i1])
            for n2 in 0:p-2
                push!(f,qlam[face]*k1*_ker(n2,subs[i2]))
            end
        end
    end
    tlams=(product[1]*l1,product[8]*l3)
    for face in 1:2
        lam=tlams[face]
        for n1 in 0:p-3
            k1=_ker(n1,subs[1])
            for n2 in 0:p-3-n1
                push!(f,lam*k1*_ker(n2,subs[5]))
            end
        end
    end
    product_bubble=l1*l2*l3
    b=_BF[]
    sizehint!(b,(p-1)*(p-2)*(p-1)÷2)
    for n1 in 0:p-3
        k1=_ker(n1,subs[1])
        for n2 in 0:p-3-n1
            k2=_ker(n2,subs[5])
            for n3 in 2:p
                push!(b,product_bubble*k1*k2*_lob(n3,wb))
            end
        end
    end
    return vb,e,f,b
end

# ---------------------------------------------------------------------------
# Orientation machinery (ports of orientEdgeFunctionsForNegativeFlag,
# orientOneFace, orientFace, addAllOrientedFaceFunctions, and
# MFace::getOrientationFlagForFace / numberOrientation*Face).
# ---------------------------------------------------------------------------

# Flag triples driving the all-orientation tables, in Gmsh's case order.
const _QUAD_ORIENTATION_FLAGS=
    ((1,1,1),(-1,1,1),(1,-1,1),(-1,-1,1),
     (1,1,-1),(-1,1,-1),(1,-1,-1),(-1,-1,-1))
const _TRI_ORIENTATION_FLAGS=
    ((0,1),(1,1),(2,1),(0,-1),(1,-1),(2,-1))

@inline function _number_orientation_quad(f1::Int,f2::Int,f3::Int)
    f1==1 && f2==1 && f3==1 && return 0
    f1==-1 && f2==1 && f3==1 && return 1
    f1==1 && f2==-1 && f3==1 && return 2
    f1==-1 && f2==-1 && f3==1 && return 3
    f1==1 && f2==1 && f3==-1 && return 4
    f1==-1 && f2==1 && f3==-1 && return 5
    f1==1 && f2==-1 && f3==-1 && return 6
    return 7
end

@inline function _number_orientation_tri(f1::Int,f2::Int)
    f1==0 && f2==1 && return 0
    f1==1 && f2==1 && return 1
    f1==2 && f2==1 && return 2
    f1==0 && f2==-1 && return 3
    f1==1 && f2==-1 && return 4
    return 5
end

# MFace::getOrientationFlagForFace on the Solin face whose element-vertex
# numbers are `nums` (permuted vertex numbers in face-local order).
function _tri_face_flags(nums::NTuple{3,Int})
    si=sortperm([nums[1],nums[2],nums[3]])
    if nums[si[1]]==nums[1]
        return 0,nums[si[2]]==nums[2] ? 1 : -1
    elseif nums[1+1]==nums[si[1]]
        return 1,nums[1]==nums[si[3]] ? 1 : -1
    end
    return 2,nums[2]==nums[si[3]] ? 1 : -1
end

function _quad_face_flags(nums::NTuple{4,Int})
    si=sortperm([nums[1],nums[2],nums[3],nums[4]])
    c=findfirst(i->nums[si[1]]==nums[i],1:4)-1
    opposed=3-c
    num_opposed=nums[opposed+1]
    axis1A=nums[si[1]]
    axis1B=nums[si[2]]==num_opposed ? nums[si[3]] : nums[si[2]]
    if axis1A==nums[1]
        return (1,1,axis1B==nums[2] ? 1 : -1)
    elseif axis1A==nums[2]
        return (-1,1,axis1B==nums[1] ? 1 : -1)
    elseif axis1A==nums[3]
        return (1,-1,axis1B==nums[4] ? 1 : -1)
    end
    return (-1,-1,axis1B==nums[3] ? 1 : -1)
end

# Permute a triangular face's vertex lambdas per (flag1, flag2).
@inline function _permute_tri_lambdas(lam::NTuple{3,_BF},f1::Int,f2::Int)
    a,b,c=lam
    f1==1 && f2==-1 && return (b,a,c)
    f1==0 && f2==-1 && return (a,c,b)
    f1==2 && f2==-1 && return (c,b,a)
    f1==1 && f2==1 && return (b,c,a)
    f1==2 && f2==1 && return (c,a,b)
    return lam
end

# Recompute one triangular face's functions under (flag1,flag2): the product
# factor multiplies kernel functions of the two face coordinate differences.
function _oriented_h1_tri_face!(out::Vector{_BF},offset::Int,pf::Int,
                                lam::NTuple{3,_BF},extra::_BF,
                                f1::Int,f2::Int)
    (f1==0 && f2==1) && return out # identity: canonical block stays
    l=_permute_tri_lambdas(lam,f1,f2)
    subs1=l[2]-l[1]
    subs2=l[1]-l[3]
    product=l[1]*l[2]*l[3]*extra
    phi2=_BF[_ker(n,subs2) for n in 0:max(pf-3,-1)]
    iterator=offset
    for n1 in 0:pf-3
        k1=_ker(n1,subs1)
        for n2 in 0:pf-3-n1
            iterator+=1
            out[iterator]=product*k1*phi2[n2+1]
        end
    end
    return out
end

# Recompute one quadrilateral face's functions under (flag1,flag2,flag3).
# `transverse` is the face-normal factor, `arg1`/`arg2` the two in-face
# polynomial arguments, and `kernel` selects kernel polynomials (prism) over
# Lobatto polynomials (quadrangle, hexahedron).
function _oriented_h1_quad_face!(out::Vector{_BF},offset::Int,
                                 pf1::Int,pf2::Int,
                                 transverse::_BF,arg1::_BF,arg2::_BF,
                                 kernel::Bool,f1::Int,f2::Int,f3::Int)
    (f1==1 && f2==1 && f3==1) && return out # identity
    poly1=kernel ? (k,x)->_ker(k,x) : (k,x)->_lob(k,x)
    poly2=poly1
    if f3==1
        iterator=offset
        for it1 in 2:pf1
            i1=isodd(it1) ? f1 : 1
            for it2 in 2:pf2
                i2=isodd(it2) ? f2 : 1
                iterator+=1
                out[iterator]=i1*i2*out[iterator]
            end
        end
    else
        # Lobatto arguments are degrees 2..pf; kernel arguments are Gmsh's
        # 0-based EvalKernelFunction indices 0..pf-2.
        shift=kernel ? 2 : 0
        lk1=_BF[poly1(k-shift,arg1) for k in 2:pf1]
        lk2=_BF[poly2(k-shift,arg2) for k in 2:pf2]
        iterator=offset
        for it1 in 2:pf2
            i1=isodd(it1) ? f2 : 1
            for it2 in 2:pf1
                i2=isodd(it2) ? f1 : 1
                iterator+=1
                out[iterator]=i1*i2*transverse*lk1[it2-1]*lk2[it1-1]
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Per-family face recompute ingredients for orientOneFace.
# Tri faces return (lambdas in Solin vertex order, extra factor).
# Quad faces return (transverse factor, var1, var2, use_kernel).
# ---------------------------------------------------------------------------
@inline function _lams_tri(u::Float64,v::Float64)
    return (_BF(v,(0.0,1.0,0.0)),
            _BF(1.0-u-v,(-1.0,-1.0,0.0)),
            _BF(u,(1.0,0.0,0.0)))
end
@inline function _lams_tet(u::Float64,v::Float64,w::Float64)
    return (_BF(v,(0.0,1.0,0.0)),
            _BF(1.0-u-v-w,(-1.0,-1.0,-1.0)),
            _BF(u,(1.0,0.0,0.0)),
            _BF(w,(0.0,0.0,1.0)))
end
@inline function _lams_qua(u::Float64,v::Float64)
    return (_BF(0.5*(1.0+u),(0.5,0.0,0.0)),
            _BF(0.5*(1.0-u),(-0.5,0.0,0.0)),
            _BF(0.5*(1.0+v),(0.0,0.5,0.0)),
            _BF(0.5*(1.0-v),(0.0,-0.5,0.0)))
end
@inline function _lams_hex(u::Float64,v::Float64,w::Float64)
    return (_BF(0.5*(1.0+u),(0.5,0.0,0.0)),
            _BF(0.5*(1.0-u),(-0.5,0.0,0.0)),
            _BF(0.5*(1.0+v),(0.0,0.5,0.0)),
            _BF(0.5*(1.0-v),(0.0,-0.5,0.0)),
            _BF(0.5*(1.0+w),(0.0,0.0,0.5)),
            _BF(0.5*(1.0-w),(0.0,0.0,-0.5)))
end
@inline function _lams_pri(u::Float64,v::Float64,w::Float64)
    return (_BF(v,(0.0,1.0,0.0)),
            _BF(1.0-u-v,(-1.0,-1.0,0.0)),
            _BF(u,(1.0,0.0,0.0)),
            _BF(0.5*(1.0+w),(0.0,0.0,0.5)),
            _BF(0.5*(1.0-w),(0.0,0.0,-0.5)))
end

# Tri-face lambdas in Solin face-vertex order plus the extra product factor.
function _h1_tri_face_desc(::Val{:tri},face::Int,u,v,w)
    l=_lams_tri(u,v)
    return (l[2],l[3],l[1]),_BF(1.0)
end
function _h1_tri_face_desc(::Val{:tet},face::Int,u,v,w)
    l=_lams_tet(u,v,w)
    desc=((l[2],l[3],l[1]),(l[2],l[3],l[4]),(l[2],l[1],l[4]),
          (l[3],l[1],l[4]))
    return desc[face],_BF(1.0)
end
function _h1_tri_face_desc(::Val{:pri},face::Int,u,v,w)
    l=_lams_pri(u,v,w)
    lam=(l[2],l[3],l[1])
    return lam,face==4 ? l[5] : l[4]
end

# Quad-face transverse factor, in-face arguments, and polynomial kind.
function _h1_quad_face_desc(::Val{:qua},face::Int,u,v,w)
    ub=_BF(u,(1.0,0.0,0.0))
    vbb=_BF(v,(0.0,1.0,0.0))
    return _BF(1.0),ub,vbb,false
end
function _h1_quad_face_desc(::Val{:hex},face::Int,u,v,w)
    l=_lams_hex(u,v,w)
    ub=_BF(u,(1.0,0.0,0.0))
    vbb=_BF(v,(0.0,1.0,0.0))
    wb=_BF(w,(0.0,0.0,1.0))
    transverse=(l[6],l[4],l[2],l[1],l[3],l[5])
    vars=((ub,vbb),(ub,wb),(vbb,wb),(vbb,wb),(ub,wb),(ub,vbb))
    var1,var2=vars[face]
    return transverse[face],var1,var2,false
end
function _h1_quad_face_desc(::Val{:pri},face::Int,u,v,w)
    l=_lams_pri(u,v,w)
    prods=(l[2]*l[3]*l[4]*l[5],l[2]*l[1]*l[4]*l[5],l[1]*l[3]*l[4]*l[5])
    var1=(l[3]-l[2],l[1]-l[2],l[1]-l[3])
    var2=l[4]-l[5]
    return prods[face],var1[face],var2,true
end

# ---------------------------------------------------------------------------
# Drivers: canonical + negative-edge + all-orientation face tables per point,
# then the Gmsh orientation loop (edge flags from Solin edge vertex order,
# face flags from the MFace algorithm).
# ---------------------------------------------------------------------------
const _H1_FAMILIES=(:pnt,:lin,:tri,:qua,:tet,:hex,:pri)

@inline function _h1_tables(family::Symbol,p::Int,u,v,w)
    family===:pnt && return _h1_tables(Val(:pnt),p,u,v,w)
    family===:lin && return _h1_tables(Val(:lin),p,u,v,w)
    family===:tri && return _h1_tables(Val(:tri),p,u,v,w)
    family===:qua && return _h1_tables(Val(:qua),p,u,v,w)
    family===:tet && return _h1_tables(Val(:tet),p,u,v,w)
    family===:hex && return _h1_tables(Val(:hex),p,u,v,w)
    return _h1_tables(Val(:pri),p,u,v,w)
end

# Negative-edge-flag table: flip functions at odd local index within each
# edge's contiguous block (the shared H1 rule for every family).
function _h1_negative_edge_table(e::Vector{_BF},edge_funcs::Int)
    out=copy(e)
    cursor=0
    while cursor<length(out)
        for k in 1:edge_funcs
            isodd(k-1) && (out[cursor+k]=-1.0*out[cursor+k])
        end
        cursor+=edge_funcs
    end
    return out
end

function _h1_orient_one_face!(work::Vector{_BF},family::Symbol,p::Int,
                              face::Int,block::Tuple{Int,Int,Symbol},
                              f1::Int,f2::Int,f3::Int,
                              u::Float64,v::Float64,w::Float64)
    offset,_,kind=block
    if kind===:tri
        lam,extra=_h1_tri_face_desc(Val(family),face,u,v,w)
        _oriented_h1_tri_face!(work,offset,p,lam,extra,f1,f2)
    else
        transverse,arg1,arg2,kernel=
            _h1_quad_face_desc(Val(family),face,u,v,w)
        _oriented_h1_quad_face!(
            work,offset,p,p,transverse,arg1,arg2,kernel,f1,f2,f3)
    end
    return work
end

# All oriented face tables for one evaluation point:
# quadAll[o][i] = concatenated quad-face functions under quad orientation o,
# triAll[o][i] = concatenated tri-face functions under tri orientation o.
function _h1_all_oriented_faces(family::Symbol,p::Int,
                                u::Float64,v::Float64,w::Float64,
                                f::Vector{_BF})
    blocks=_h1_face_blocks(family,p)
    quad_total=0
    tri_total=0
    for block in blocks
        block[3]===:quad ? (quad_total+=block[2]) : (tri_total+=block[2])
    end
    quad_all=Vector{Vector{_BF}}(undef,8)
    for o in 1:8
        work=copy(f)
        f1,f2,f3=_QUAD_ORIENTATION_FLAGS[o]
        for (face,block) in enumerate(blocks)
            block[3]===:quad || continue
            _h1_orient_one_face!(
                work,family,p,face,block,f1,f2,f3,u,v,w)
        end
        quad_all[o]=work[1:quad_total]
    end
    tri_all=Vector{Vector{_BF}}(undef,6)
    for o in 1:6
        work=copy(f)
        f1,f2=_TRI_ORIENTATION_FLAGS[o]
        for (face,block) in enumerate(blocks)
            block[3]===:tri || continue
            _h1_orient_one_face!(
                work,family,p,face,block,f1,f2,1,u,v,w)
        end
        tri_all[o]=work[quad_total+1:end]
    end
    return quad_all,tri_all,quad_total
end

# Evaluate the whole oriented hierarchy for one point and one vertex
# permutation, writing function-major results into `out` (1:values,
# 2..4:gradient components).
function _h1_write_point!(out::Vector{Float64},cursor::Int,
                          vt::Vector{_BF},et::Vector{_BF},
                          eneg::Vector{_BF},ftab::Vector{_BF},
                          bt::Vector{_BF},
                          quad_all,tri_all,quad_total::Int,
                          family::Symbol,p::Int,perm::Vector{Int},
                          components::Int)
    blocks=_h1_face_blocks(family,p)
    edge_funcs=p-1
    # edge selection by Solin edge flag
    ecopy=et
    edges=get(_SOLIN_EDGES,family,())
    oriented=false
    eout=et
    if !isempty(edges)
        eout=copy(et)
        for (edge,(a,b)) in enumerate(edges)
            flag=perm[b+1]<perm[a+1] ? -1 : 1
            if flag==-1
                oriented=true
                lo=(edge-1)*edge_funcs
                for k in 1:edge_funcs
                    eout[lo+k]=eneg[lo+k]
                end
            end
        end
    end
    # face selection by MFace flags
    fout=ftab
    if !isempty(blocks)
        fout=copy(ftab)
        faces=_SOLIN_FACES[family]
        for (face,verts) in enumerate(faces)
            offset,count,kind=blocks[face]
            count==0 && continue
            if kind===:tri
                nums=ntuple(i->perm[verts[i]+1],3)
                f1,f2=_tri_face_flags(nums)
                io=_number_orientation_tri(f1,f2)
                src=tri_all[io+1]
                base=offset-quad_total
                for k in 1:count
                    fout[offset+k]=src[base+k]
                end
            else
                nums=ntuple(i->perm[verts[i]+1],4)
                f1,f2,f3=_quad_face_flags(nums)
                io=_number_orientation_quad(f1,f2,f3)
                src=quad_all[io+1]
                for k in 1:count
                    fout[offset+k]=src[offset+k]
                end
            end
        end
    end
    if components==1
        for f in vt
            out[cursor+=1]=f.v
        end
        for f in eout
            out[cursor+=1]=f.v
        end
        for f in fout
            out[cursor+=1]=f.v
        end
        for f in bt
            out[cursor+=1]=f.v
        end
    else
        for f in vt
            for c in 1:3
                out[cursor+=1]=f.g[c]
            end
        end
        for f in eout
            for c in 1:3
                out[cursor+=1]=f.g[c]
            end
        end
        for f in fout
            for c in 1:3
                out[cursor+=1]=f.g[c]
            end
        end
        for f in bt
            for c in 1:3
                out[cursor+=1]=f.g[c]
            end
        end
    end
    return cursor
end

# ---------------------------------------------------------------------------
# getKeysInfo tables: (functionType, order) per basis function in output order.
# ---------------------------------------------------------------------------
function _h1_keys_info(family::Symbol,p::Int)
    nv=family===:pnt ? 1 :
       family===:lin ? 2 : family===:tri ? 3 :
       (family===:tet || family===:qua) ? 4 : family===:hex ? 8 : 6
    counts=_h1_counts(family,p)
    info=Vector{Tuple{Int32,Int32}}(undef,
        counts.vertex+counts.edge+counts.face+counts.bubble)
    it=1
    for _ in 1:nv
        info[it]=(Int32(0),Int32(family===:pnt ? 0 : 1))
        it+=1
    end
    nedge=length(get(_SOLIN_EDGES,family,()))
    for _ in 1:nedge, i in 2:p
        info[it]=(Int32(1),Int32(i))
        it+=1
    end
    faces=get(_SOLIN_FACES,family,())
    for face in faces
        if length(face)==4
            for n1 in 2:p, n2 in 2:p
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
        else
            for n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0)
                info[it]=(Int32(2),Int32(n1+n2+1))
                it+=1
            end
        end
    end
    if family===:tet
        for n1 in 1:max(p-3,0), n2 in 1:max(p-2-n1,0),
            n3 in 1:max(p-1-n1-n2,0)
            info[it]=(Int32(3),Int32(n1+n2+n3+1))
            it+=1
        end
    elseif family===:hex
        for n1 in 2:p, n2 in 2:p, n3 in 2:p
            info[it]=(Int32(3),Int32(max(n1,n2,n3)))
            it+=1
        end
    elseif family===:pri
        for n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0), n3 in 2:p
            info[it]=(Int32(3),Int32(max(n1+n2+1,n3)))
            it+=1
        end
    end
    return info
end

# typeKeys per function position: vertex funcs take 0; per edge the funcs take
# 1..per-edge count; each face takes the per-face range starting after the
# per-edge count; bubbles take the remaining upper range.
function _h1_type_keys(family::Symbol,p::Int)
    counts=_h1_counts(family,p)
    per_edge=max(p-1,0)
    blocks=_h1_face_blocks(family,p)
    per_quad=0
    per_tri=0
    for block in blocks
        if block[3]===:quad
            per_quad=block[2]
        else
            per_tri=block[2]
        end
    end
    const1=per_edge+1
    const2=const1+per_quad
    const3=const1+per_tri
    const4=counts.bubble+max(const2,const3)
    keys=Vector{Int32}(undef,
        counts.vertex+counts.edge+counts.face+counts.bubble)
    it=1
    for _ in 1:counts.vertex
        keys[it]=Int32(0)
        it+=1
    end
    for _ in 1:length(get(_SOLIN_EDGES,family,())), k in 1:per_edge
        keys[it]=Int32(k)
        it+=1
    end
    for block in blocks
        limit=block[3]===:quad ? const2 : const3
        for k in const1:limit-1
            keys[it]=Int32(k)
            it+=1
        end
    end
    for k in max(const3,const2):const4-1
        keys[it]=Int32(k)
        it+=1
    end
    return keys
end

# Full oriented H1/GradH1 basis in Gmsh's output layout for the requested
# orientations: orientation-major, then point, then function, then component.
function _h1_basis_eval(family::Symbol,p::Int,gradient::Bool,
                        coords::Vector{Float64},point_count::Int,
                        orientations::Vector{Int},result::Vector{Float64})
    components=gradient ? 3 : 1
    counts=_h1_counts(family,p)
    nfuncs=counts.vertex+counts.edge+counts.face+counts.bubble
    blocks=_h1_face_blocks(family,p)
    nvertices=counts.vertex
    edges=get(_SOLIN_EDGES,family,())
    per_edge=p-1
    # per-point canonical tables and auxiliary orientation tables
    vtabs=Vector{Vector{_BF}}(undef,point_count)
    etabs=Vector{Vector{_BF}}(undef,point_count)
    enegs=Vector{Vector{_BF}}(undef,point_count)
    ftabs=Vector{Vector{_BF}}(undef,point_count)
    btabs=Vector{Vector{_BF}}(undef,point_count)
    quads=Vector{Any}(undef,point_count)
    tris=Vector{Any}(undef,point_count)
    quadtotals=Vector{Int}(undef,point_count)
    @inbounds for q in 1:point_count
        u=coords[3q-2]
        v=coords[3q-1]
        w=coords[3q]
        vt,et,ft,bt=_h1_tables(family,p,u,v,w)
        vtabs[q]=vt
        etabs[q]=et
        enegs[q]=isempty(et) ? et : _h1_negative_edge_table(et,per_edge)
        ftabs[q]=ft
        btabs[q]=bt
        if isempty(ft)
            quads[q]=nothing
            tris[q]=nothing
            quadtotals[q]=0
        else
            qa,ta,qt=_h1_all_oriented_faces(family,p,u,v,w,ft)
            quads[q]=qa
            tris[q]=ta
            quadtotals[q]=qt
        end
    end
    cursor=0
    for orientation in orientations
        perm=_orientation_permutation(nvertices,orientation)
        for q in 1:point_count
            cursor=_h1_write_point!(
                result,cursor,vtabs[q],etabs[q],enegs[q],ftabs[q],
                btabs[q],quads[q],tris[q],quadtotals[q],
                family,p,perm,components)
        end
    end
    return result
end

# ---------------------------------------------------------------------------
# H(curl) generators. Basis functions are reference 3-vectors; counts follow
# the Gmsh constructors (uniform order p):
#   lin: e=p+1
#   tri: e=3(p+1); f=3(p-1)+(p-1)(p-2) (0 for p=0)
#   qua: e=4(p+1); f=2p(p+1)
#   tet: e=6(p+1); f=12(p-1)+4(p-2)(p-1); b=(p-1)(p-2)(p-3)/2+2(p-2)(p-1)
#   hex: e=12(p+1); f=12p(p+1); b=3p^2(p+1)
#   pri: e=9(p+1); f=6p(p+1)+6(p-1)+2(p-1)(p-2);
#        b=3(p-1)p+(p-1)(p-2)p+(p-1)p(p+1)/2
# ---------------------------------------------------------------------------
const _V3=NTuple{3,Float64}

function _hcurl_counts(family::Symbol,p::Int)
    ne=length(get(_SOLIN_EDGES,family,()))
    e=ne*(p+1)
    if family===:tri
        f=p>=1 ? 3*(p-1)+(p-1)*(p-2) : 0
        return (vertex=0,edge=e,face=f,bubble=0)
    elseif family===:qua
        return (vertex=0,edge=e,face=2*p*(p+1),bubble=0)
    elseif family===:tet
        f=p>=1 ? 12*(p-1)+4*(p-2)*(p-1) : 0
        b=p>=1 ? (p-1)*(p-2)*(p-3)÷2+2*(p-2)*(p-1) : 0
        return (vertex=0,edge=e,face=f,bubble=b)
    elseif family===:hex
        return (vertex=0,edge=e,face=12*p*(p+1),bubble=3*p*p*(p+1))
    elseif family===:pri
        fq=6*p*(p+1)
        ft=p>=1 ? 6*(p-1)+2*(p-1)*(p-2) : 0
        b=3*max(p-1,0)*p+max(p-1,0)*max(p-2,0)*p+
          max(p-1,0)*p*(p+1)÷2
        return (vertex=0,edge=e,face=fq+ft,bubble=b)
    end
    return (vertex=0,edge=e,face=0,bubble=0)
end

# (offset,count,:tri|:quad) blocks in the concatenated Hcurl face table.
function _hcurl_face_blocks(family::Symbol,p::Int)
    blocks=Tuple{Int,Int,Symbol}[]
    offset=0
    if family===:qua
        push!(blocks,(0,2*p*(p+1),:quad))
        return blocks
    elseif family===:tri
        push!(blocks,(0,p>=1 ? 3*(p-1)+(p-1)*(p-2) : 0,:tri))
        return blocks
    elseif family===:tet
        per=p>=1 ? 3*(p-1)+(p-1)*(p-2) : 0
        for _ in 1:4
            push!(blocks,(offset,per,:tri))
            offset+=per
        end
        return blocks
    elseif family===:hex
        per=2*p*(p+1)
        for _ in 1:6
            push!(blocks,(offset,per,:quad))
            offset+=per
        end
        return blocks
    elseif family===:pri
        per=2*p*(p+1)
        for _ in 1:3
            push!(blocks,(offset,per,:quad))
            offset+=per
        end
        pert=p>=1 ? 3*(p-1)+(p-1)*(p-2) : 0
        for _ in 1:2
            push!(blocks,(offset,pert,:tri))
            offset+=pert
        end
        return blocks
    end
    return blocks
end

# --- Hcurl Line --------------------------------------------------------------
function _hcurl_tables(::Val{:lin},p::Int,u,v,w)
    psie0=1.0
    psie1=u
    sub=u
    e=_V3[(1.0,0.0,0.0)]
    if p>=1
        push!(e,(psie1,0.0,0.0))
        leg=[_gmsh_legendre(k,sub) for k in 0:max(p-1,0)]
        for iedge in 2:p
            val=Float64(Float32(2*iedge-1)/Float32(iedge))*leg[iedge]*psie1-
                Float64(Float32(iedge-1)/Float32(iedge))*leg[iedge-1]*psie0
            push!(e,(val,0.0,0.0))
        end
    end
    return e,_V3[],_V3[]
end

function _curl_tables(::Val{:lin},p::Int,u,v,w)
    return fill((0.0,0.0,0.0),p+1),_V3[],_V3[]
end

# --- Hcurl Quadrangle --------------------------------------------------------
function _hcurl_tables(::Val{:qua},p::Int,u,v,w)
    l1=0.5*(1.0+u)
    l2=0.5*(1.0-u)
    l3=0.5*(1.0+v)
    l4=0.5*(1.0-v)
    legu=[_gmsh_legendre(k,u) for k in 0:p]
    legv=[_gmsh_legendre(k,v) for k in 0:p]
    e=_V3[]
    sizehint!(e,4*(p+1))
    for k in 0:p
        push!(e,(l4*legu[k+1],0.0,0.0))
    end
    for k in 0:p
        push!(e,(0.0,l1*legv[k+1],0.0))
    end
    for k in 0:p
        push!(e,(l3*legu[k+1],0.0,0.0))
    end
    for k in 0:p
        push!(e,(0.0,l2*legv[k+1],0.0))
    end
    f=_V3[]
    sizehint!(f,2*p*(p+1))
    for n1 in 0:p
        for n2 in 2:p+1
            push!(f,(legu[n1+1]*_gmsh_lobatto(n2,v),0.0,0.0))
        end
    end
    for n1 in 2:p+1
        for n2 in 0:p
            push!(f,(0.0,legv[n2+1]*_gmsh_lobatto(n1,u),0.0))
        end
    end
    return e,f,_V3[]
end

function _curl_tables(::Val{:qua},p::Int,u,v,w)
    e=_V3[]
    sizehint!(e,4*(p+1))
    legu=[_gmsh_legendre(k,u) for k in 0:p]
    legv=[_gmsh_legendre(k,v) for k in 0:p]
    for k in 0:p
        push!(e,(0.0,0.0,0.5*legu[k+1]))
    end
    for k in 0:p
        push!(e,(0.0,0.0,0.5*legv[k+1]))
    end
    for k in 0:p
        push!(e,(0.0,0.0,-0.5*legu[k+1]))
    end
    for k in 0:p
        push!(e,(0.0,0.0,-0.5*legv[k+1]))
    end
    f=_V3[]
    sizehint!(f,2*p*(p+1))
    for n1 in 0:p
        for n2 in 2:p+1
            push!(f,(0.0,0.0,-legu[n1+1]*_gmsh_dlobatto(n2,v)))
        end
    end
    for n1 in 2:p+1
        for n2 in 0:p
            push!(f,(0.0,0.0,legv[n2+1]*_gmsh_dlobatto(n1,u)))
        end
    end
    return e,f,_V3[]
end

# --- Hcurl Triangle ----------------------------------------------------------
# Whitney edge vectors on the mapped triangle (jacob=2 rescales to the Gmsh
# reference domain). n/t use the 2-D dot product; component 3 stays zero.
function _hcurl_tri_psi(u::Float64,v::Float64)
    l1=v
    l2=1.0-u-v
    l3=u
    s05=sqrt(0.5)
    n1=(0.0,1.0,0.0)
    n2=(-s05,-s05,0.0)
    n3=(1.0,0.0,0.0)
    # dot(n2,t1)=-s05; dot(n3,t1)=1; dot(n3,t2)=-1; dot(n1,t2)=1;
    # dot(n1,t3)=-1; dot(n2,t3)=s05
    psie0=(ntuple(j->-l3*n2[j]/s05+l2*n3[j],3),
           ntuple(j->-l1*n3[j]+l3*n1[j],3),
           ntuple(j->-l2*n1[j]+l1*n2[j]/s05,3))
    psie1=(ntuple(j->-l3*n2[j]/s05-l2*n3[j],3),
           ntuple(j->-l1*n3[j]-l3*n1[j],3),
           ntuple(j->-l2*n1[j]-l1*n2[j]/s05,3))
    return psie0,psie1,(l3-l2,l1-l3,l2-l1)
end

function _hcurl_tables(::Val{:tri},p::Int,u,v,w)
    psie0,psie1,subs=_hcurl_tri_psi(u,v)
    jacob=2.0
    e=_V3[]
    sizehint!(e,3*(p+1))
    legs=[[_gmsh_legendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:3]
    for i in 1:3
        push!(e,jacob .* psie0[i])
        p>=1 || continue
        push!(e,jacob .* psie1[i])
        for iedge in 2:p
            val=jacob .* (Float64(Float32(2*iedge-1)/Float32(iedge))*legs[i][iedge] .* psie1[i] .-
                          Float64(Float32(iedge-1)/Float32(iedge))*legs[i][iedge-1] .* psie0[i])
            push!(e,val)
        end
    end
    l1=v
    l2=1.0-u-v
    l3=u
    f=_V3[]
    sizehint!(f,p>=1 ? 3*(p-1)+(p-1)*(p-2) : 0)
    # edge-based face functions (one group per edge)
    products=(l3*l2,l1*l3,l1*l2)
    nds=((0.0,0.5,0.0),(-0.5,-0.5,0.0),(0.5,0.0,0.0))
    for i in 1:3
        for i1 in 2:p
            val=jacob*products[i]*legs[i][i1-1] .* nds[i]
            push!(f,val)
        end
    end
    # genuine face functions: x group then mirrored y group
    product123=l1*l2*l3
    genuine=Float64[]
    for n1 in 0:p-3
        for n2 in 0:p-3-n1
            push!(genuine,
                  0.5*jacob*product123*legs[1][n1+1]*legs[3][n2+1])
            push!(f,(genuine[end],0.0,0.0))
        end
    end
    for value in genuine
        push!(f,(0.0,value,0.0))
    end
    return e,f,_V3[]
end

function _curl_tables(::Val{:tri},p::Int,u,v,w)
    l1=v
    l2=1.0-u-v
    l3=u
    det=4.0
    psie0,psie1,subs=_hcurl_tri_psi(u,v)
    curlpsie0=(1.0,1.0,1.0)
    curlpsie1=(0.0,0.0,0.0)
    # dlambda values on the mapped domain, matching the source's 2-D table
    dsub=((1.0,0.5),(-0.5,0.5),(-0.5,-1.0))
    legs=[[_gmsh_legendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:3]
    dlegs=[[_gmsh_dlegendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:3]
    e=_V3[]
    sizehint!(e,3*(p+1))
    for i in 1:3
        push!(e,(0.0,0.0,det*curlpsie0[i]))
        p>=1 || continue
        push!(e,(0.0,0.0,det*curlpsie1[i]))
        for iedge in 2:p
            val=det*(Float64(Float32(2*iedge-1)/Float32(iedge))*
                     (dsub[i][1]*dlegs[i][iedge]*psie1[i][2]-
                      dsub[i][2]*dlegs[i][iedge]*psie1[i][1])-
                     Float64(Float32(iedge-1)/Float32(iedge))*
                     (curlpsie0[i]*legs[i][iedge-1]+
                      dsub[i][1]*dlegs[i][iedge-1]*psie0[i][2]-
                      dsub[i][2]*dlegs[i][iedge-1]*psie0[i][1]))
            push!(e,(0.0,0.0,val))
        end
    end
    f=_V3[]
    sizehint!(f,p>=1 ? 3*(p-1)+(p-1)*(p-2) : 0)
    dlambda23=0.5*(l2-l3)
    prod32=l3*l2
    for n1 in 2:p
        push!(f,(0.0,0.0,
                 0.5*det*(dlambda23*legs[1][n1-1]+
                          prod32*dsub[1][1]*dlegs[1][n1-1])))
    end
    dlambda13U=0.5*l1
    dlambda13V=0.5*l3
    prod13=l3*l1
    for n1 in 2:p
        push!(f,(0.0,0.0,
                 -det*0.5*(dlambda13U*legs[2][n1-1]+
                           prod13*dsub[2][1]*dlegs[2][n1-1]-
                           (dlambda13V*legs[2][n1-1]+
                            prod13*dsub[2][2]*dlegs[2][n1-1]))))
    end
    dlambda12=0.5*(l2-l1)
    prod12=l2*l1
    for n1 in 2:p
        push!(f,(0.0,0.0,
                 -0.5*det*(dlambda12*legs[3][n1-1]+
                           prod12*dsub[3][2]*dlegs[3][n1-1])))
    end
    prod123=l1*l2*l3
    dlambda123U=0.5*l1*(l2-l3)
    dlambda123V=0.5*l3*(l2-l1)
    for n1 in 0:p-3
        for n2 in 0:p-3-n1
            push!(f,(0.0,0.0,
                     -0.5*det*(dlambda123V*legs[1][n1+1]*legs[3][n2+1]+
                               prod123*dsub[1][2]*dlegs[1][n1+1]*legs[3][n2+1]+
                               prod123*dsub[3][2]*legs[1][n1+1]*dlegs[3][n2+1])))
        end
    end
    for n1 in 0:p-3
        for n2 in 0:p-3-n1
            push!(f,(0.0,0.0,
                     0.5*det*(dlambda123U*legs[1][n1+1]*legs[3][n2+1]+
                              prod123*dsub[1][1]*dlegs[1][n1+1]*legs[3][n2+1]+
                              prod123*dsub[3][1]*legs[1][n1+1]*dlegs[3][n2+1])))
        end
    end
    return e,f,_V3[]
end

# --- Hcurl Tetrahedron -------------------------------------------------------
const _TET_N=((0.0,1.0,0.0),
              (-1.0/sqrt(3.0),-1.0/sqrt(3.0),-1.0/sqrt(3.0)),
              (1.0,0.0,0.0),(0.0,0.0,1.0))
const _TET_T=((1.0,0.0,0.0),(-1.0,1.0,0.0),(0.0,-1.0,0.0),
              (0.0,0.0,1.0),(-1.0,0.0,1.0),(0.0,-1.0,1.0))
# Whitney (lambda_a, n_a, lambda_b, n_b, t) per edge: psie0 = la*na/(na.t) +
# lb*nb/(nb.t); psie1 swaps + to -.
const _TET_EDGE_WHITNEY=
    ((3,2,2,3,1),(1,3,3,1,2),(2,1,1,2,3),
     (4,2,2,4,4),(4,1,1,4,6),(4,3,3,4,5))
const _TET_CURL_PSIE0=
    ((0.0,-1.0,1.0),(0.0,0.0,1.0),(-1.0,0.0,1.0),
     (-1.0,1.0,0.0),(1.0,0.0,0.0),(0.0,-1.0,0.0))

function _hcurl_tet_terms(u::Float64,v::Float64,w::Float64)
    l1=v
    l2=1.0-u-v-w
    l3=u
    l4=w
    lam=(l1,l2,l3,l4)
    dlams=((0.0,0.5,0.0),(-0.5,-0.5,-0.5),(0.5,0.0,0.0),(0.0,0.0,0.5))
    psie0=Vector{_V3}(undef,6)
    psie1=Vector{_V3}(undef,6)
    for (e,(la,na,lb,nb,t)) in enumerate(_TET_EDGE_WHITNEY)
        na_v=_TET_N[na]
        nb_v=_TET_N[nb]
        t_v=_TET_T[t]
        dna=na_v[1]*t_v[1]+na_v[2]*t_v[2]+na_v[3]*t_v[3]
        dnb=nb_v[1]*t_v[1]+nb_v[2]*t_v[2]+nb_v[3]*t_v[3]
        a=lam[la]
        b=lam[lb]
        psie0[e]=ntuple(j->a*na_v[j]/dna+b*nb_v[j]/dnb,3)
        psie1[e]=ntuple(j->a*na_v[j]/dna-b*nb_v[j]/dnb,3)
    end
    subs=(l3-l2,l1-l3,l2-l1,l4-l2,l4-l1,l4-l3,l2-l4,l1-l2,l3-l4)
    # dsubtraction per the source's 9x3 table
    dsub=(ntuple(j->dlams[3][j]-dlams[2][j],3),
          ntuple(j->dlams[1][j]-dlams[3][j],3),
          ntuple(j->dlams[2][j]-dlams[1][j],3),
          ntuple(j->dlams[4][j]-dlams[2][j],3),
          ntuple(j->dlams[4][j]-dlams[1][j],3),
          ntuple(j->dlams[4][j]-dlams[3][j],3),
          ntuple(j->dlams[2][j]-dlams[4][j],3),
          ntuple(j->dlams[1][j]-dlams[2][j],3),
          ntuple(j->dlams[3][j]-dlams[4][j],3))
    return lam,dlams,psie0,psie1,subs,dsub
end

function _hcurl_tables(::Val{:tet},p::Int,u,v,w)
    lam,dlams,psie0,psie1,subs,dsub=_hcurl_tet_terms(u,v,w)
    l1,l2,l3,l4=lam
    jacob=2.0
    legs=[[_gmsh_legendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:9]
    e=_V3[]
    sizehint!(e,6*(p+1))
    for i in 1:6
        push!(e,jacob .* psie0[i])
        p>=1 || continue
        push!(e,jacob .* psie1[i])
        for iedge in 2:p
            val=jacob .* (Float64(Float32(2*iedge-1)/Float32(iedge))*legs[i][iedge] .* psie1[i] .-
                          Float64(Float32(iedge-1)/Float32(iedge))*legs[i][iedge-1] .* psie0[i])
            push!(e,val)
        end
    end
    # per-face: 3 edge-based groups (product, nD, legendre-index), then two
    # genuine groups, then the face-based interior block.
    face_data=(
        (((l3*l2,(0.0,0.5,0.0),1),(l1*l3,(-0.5,-0.5,-0.5),2),
          (l1*l2,(0.5,0.0,0.0),3)),
         1,3,(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5),l1*l2*l3),
        (((l3*l2,(0.0,0.0,0.5),1),(l4*l3,(-0.5,-0.5,-0.5),6),
          (l4*l2,(0.5,0.0,0.0),7)),
         1,7,(0.5,0.0,0.0),(0.0,0.0,0.5),(0.0,0.5,0.0),l4*l2*l3),
        (((l1*l2,(0.0,0.0,0.5),8),(l4*l1,(-0.5,-0.5,-0.5),5),
          (l4*l2,(0.0,0.5,0.0),7)),
         8,7,(0.0,0.5,0.0),(0.0,0.0,0.5),(0.5,0.0,0.0),l1*l4*l2),
        (((l1*l3,(0.0,0.0,0.5),2),(l4*l1,(0.5,0.0,0.0),5),
          (l4*l3,(0.0,0.5,0.0),9)),
         2,9,(0.0,0.5,0.0),(0.0,0.0,0.5),(-0.5,-0.5,-0.5),l1*l4*l3))
    f=_V3[]
    sizehint!(f,p>=1 ? 12*(p-1)+4*(p-2)*(p-1) : 0)
    b=_V3[]
    sizehint!(b,p>=1 ? (p-1)*(p-2)*(p-3)÷2+2*(p-2)*(p-1) : 0)
    for face in 1:4
        edgefuncs,iv1,iv2,tan1,tan2,nit,faceprod=face_data[face]
        for (prod,nd,idx) in edgefuncs
            for i1 in 2:p
                push!(f,(jacob*prod*legs[idx][i1-1]) .* nd)
            end
        end
        copies=Float64[]
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                value=jacob*faceprod*legs[iv1][n1+1]*legs[iv2][n2+1]
                push!(copies,value)
                push!(f,value .* tan1)
            end
        end
        for value in copies
            push!(f,value .* tan2)
        end
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                push!(b,(jacob*faceprod*legs[iv1][n1+1]*
                         legs[iv2][n2+1]) .* nit)
            end
        end
    end
    if p>3
        product=l1*l2*l3*l4
        genuine=Float64[]
        for n1 in 0:p-4
            phi1=_gmsh_kernel(n1,subs[8])
            for n2 in 0:p-4-n1
                phi2=_gmsh_kernel(n2,subs[1])
                for n3 in 0:p-4-n2-n1
                    push!(genuine,
                          jacob*product*phi1*phi2*
                          _gmsh_kernel(n3,subs[4]))
                    push!(b,(genuine[end],0.0,0.0))
                end
            end
        end
        for value in genuine
            push!(b,(0.0,value,0.0))
        end
        for value in genuine
            push!(b,(0.0,0.0,value))
        end
    end
    return e,f,b
end

# Gmsh's curlFunction: a * (grad x nD).
@inline function _curl_func(a::Float64,nD::_V3,grad::_V3)
    return (a*(nD[3]*grad[2]-nD[2]*grad[3]),
            a*(nD[1]*grad[3]-nD[3]*grad[1]),
            a*(nD[2]*grad[1]-nD[1]*grad[2]))
end

function _curl_tables(::Val{:tet},p::Int,u,v,w)
    lam,dlams,psie0,psie1,subs,dsub=_hcurl_tet_terms(u,v,w)
    l1,l2,l3,l4=lam
    detjacob=4.0
    legs=[[_gmsh_legendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:9]
    dlegs=[[_gmsh_dlegendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:9]
    e=_V3[]
    sizehint!(e,6*(p+1))
    for i in 1:6
        push!(e,detjacob .* _TET_CURL_PSIE0[i])
        if p>=1
            push!(e,(0.0,0.0,0.0))
        end
        for iedge in 2:p
            c1=Float64(Float32(2*iedge-1)/Float32(iedge))
            c2=Float64(Float32(iedge-1)/Float32(iedge))
            cp=_TET_CURL_PSIE0[i]
            val=(detjacob*(c1*(dsub[i][2]*dlegs[i][iedge]*psie1[i][3]-
                              dsub[i][3]*dlegs[i][iedge]*psie1[i][2])-
                          c2*(cp[1]*legs[i][iedge-1]+
                              dsub[i][2]*dlegs[i][iedge-1]*psie0[i][3]-
                              dsub[i][3]*dlegs[i][iedge-1]*psie0[i][2])),
                 detjacob*(c1*(dsub[i][3]*dlegs[i][iedge]*psie1[i][1]-
                              dsub[i][1]*dlegs[i][iedge]*psie1[i][3])-
                          c2*(cp[2]*legs[i][iedge-1]+
                              dsub[i][3]*dlegs[i][iedge-1]*psie0[i][1]-
                              dsub[i][1]*dlegs[i][iedge-1]*psie0[i][3])),
                 detjacob*(c1*(dsub[i][1]*dlegs[i][iedge]*psie1[i][2]-
                              dsub[i][2]*dlegs[i][iedge]*psie1[i][1])-
                          c2*(cp[3]*legs[i][iedge-1]+
                              dsub[i][1]*dlegs[i][iedge-1]*psie0[i][2]-
                              dsub[i][2]*dlegs[i][iedge-1]*psie0[i][1])))
            push!(e,val)
        end
    end
    # dfaceProduct per face: d(la*lb*lc) of the face product
    dfaceprod=(
        ntuple(i->l3*l2*dlams[1][i]+l3*dlams[2][i]*l1+dlams[3][i]*l2*l1,3),
        ntuple(i->l3*l2*dlams[4][i]+l3*dlams[2][i]*l4+dlams[3][i]*l2*l4,3),
        ntuple(i->l1*l2*dlams[4][i]+l1*dlams[2][i]*l4+dlams[1][i]*l2*l4,3),
        ntuple(i->l1*l3*dlams[4][i]+l1*dlams[3][i]*l4+dlams[1][i]*l3*l4,3))
    # Per face: (product, dproduct, nD, legendre index) for each of the three
    # edge-based groups. dproduct = grad(la*lb) = dla*lb + dlb*la.
    grad2(a,da,b,db)=ntuple(i->da[i]*b+db[i]*a,3)
    face_data=(
        (((l3*l2,grad2(l3,dlams[3],l2,dlams[2]),(0.0,0.5,0.0),1),
          (l1*l3,grad2(l3,dlams[3],l1,dlams[1]),(-0.5,-0.5,-0.5),2),
          (l1*l2,grad2(l1,dlams[1],l2,dlams[2]),(0.5,0.0,0.0),3)),
         1,3,(0.5,0.0,0.0),(0.0,0.5,0.0),(0.0,0.0,0.5),l1*l2*l3),
        (((l3*l2,grad2(l3,dlams[3],l2,dlams[2]),(0.0,0.0,0.5),1),
          (l4*l3,grad2(l4,dlams[4],l3,dlams[3]),(-0.5,-0.5,-0.5),6),
          (l4*l2,grad2(l4,dlams[4],l2,dlams[2]),(0.5,0.0,0.0),7)),
         1,7,(0.5,0.0,0.0),(0.0,0.0,0.5),(0.0,0.5,0.0),l4*l2*l3),
        (((l1*l2,grad2(l1,dlams[1],l2,dlams[2]),(0.0,0.0,0.5),8),
          (l4*l1,grad2(l4,dlams[4],l1,dlams[1]),(-0.5,-0.5,-0.5),5),
          (l4*l2,grad2(l4,dlams[4],l2,dlams[2]),(0.0,0.5,0.0),7)),
         8,7,(0.0,0.5,0.0),(0.0,0.0,0.5),(0.5,0.0,0.0),l1*l4*l2),
        (((l1*l3,grad2(l3,dlams[3],l1,dlams[1]),(0.0,0.0,0.5),2),
          (l4*l1,grad2(l4,dlams[4],l1,dlams[1]),(0.5,0.0,0.0),5),
          (l4*l3,grad2(l4,dlams[4],l3,dlams[3]),(0.0,0.5,0.0),9)),
         2,9,(0.0,0.5,0.0),(0.0,0.0,0.5),(-0.5,-0.5,-0.5),l1*l4*l3))
    f=_V3[]
    sizehint!(f,p>=1 ? 12*(p-1)+4*(p-2)*(p-1) : 0)
    b=_V3[]
    sizehint!(b,p>=1 ? (p-1)*(p-2)*(p-3)÷2+2*(p-2)*(p-1) : 0)
    for face in 1:4
        edgefuncs,iv1,iv2,tan1,tan2,nit,faceprod=face_data[face]
        for (prod,dprod,nd,idx) in edgefuncs
            for i1 in 2:p
                gleg=ntuple(i->dsub[idx][i]*dlegs[idx][i1-1],3)
                grad=ntuple(i->dprod[i]*legs[idx][i1-1]+prod*gleg[i],3)
                push!(f,_curl_func(detjacob,nd,grad))
            end
        end
        copies=_V3[]
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                gface=ntuple(3) do i
                    dfaceprod[face][i]*legs[iv1][n1+1]*legs[iv2][n2+1]+
                    faceprod*dlegs[iv1][n1+1]*dsub[iv1][i]*legs[iv2][n2+1]+
                    faceprod*legs[iv1][n1+1]*dsub[iv2][i]*dlegs[iv2][n2+1]
                end
                push!(copies,gface)
                push!(f,_curl_func(detjacob,tan1,gface))
            end
        end
        for gface in copies
            push!(f,_curl_func(detjacob,tan2,gface))
        end
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                gface=ntuple(3) do i
                    dfaceprod[face][i]*legs[iv1][n1+1]*legs[iv2][n2+1]+
                    faceprod*dlegs[iv1][n1+1]*dsub[iv1][i]*legs[iv2][n2+1]+
                    faceprod*legs[iv1][n1+1]*dsub[iv2][i]*dlegs[iv2][n2+1]
                end
                push!(b,_curl_func(detjacob,nit,gface))
            end
        end
    end
    if p>3
        l1234=l1*l2*l3*l4
        dl1234=ntuple(3) do i
            l1*l2*l3*dlams[4][i]+l1*l2*dlams[3][i]*l4+
            l1*dlams[2][i]*l3*l4+dlams[1][i]*l2*l3*l4
        end
        genuine=_V3[]
        for n1 in 0:p-4
            phi1=_gmsh_kernel(n1,subs[8])
            dphi1=_gmsh_dkernel(n1,subs[8])
            for n2 in 0:p-4-n1
                phi2=_gmsh_kernel(n2,subs[1])
                dphi2=_gmsh_dkernel(n2,subs[1])
                for n3 in 0:p-4-n2-n1
                    phi3=_gmsh_kernel(n3,subs[4])
                    dphi3=_gmsh_dkernel(n3,subs[4])
                    g=ntuple(3) do i
                        dl1234[i]*phi1*phi2*phi3+
                        l1234*phi2*dsub[8][i]*dphi1*phi3+
                        l1234*dphi2*dsub[1][i]*phi1*phi3+
                        l1234*phi2*dsub[4][i]*phi1*dphi3
                    end
                    push!(genuine,g)
                    # x-directed genuine bubble: curl = (0, g3, -g2)
                    push!(b,(0.0,detjacob*g[3],-detjacob*g[2]))
                end
            end
        end
        for g in genuine
            push!(b,(-detjacob*g[3],0.0,detjacob*g[1]))
        end
        for g in genuine
            push!(b,(detjacob*g[2],-detjacob*g[1],0.0))
        end
    end
    return e,f,b
end

# --- Hcurl Hexahedron --------------------------------------------------------
function _hcurl_tables(::Val{:hex},p::Int,u,v,w)
    uvw=(u,v,w)
    lam=(0.5*(1.0+u),0.5*(1.0-u),0.5*(1.0+v),0.5*(1.0-v),
         0.5*(1.0+w),0.5*(1.0-w))
    legs=[[_gmsh_legendre(k,uvw[i]) for k in 0:p] for i in 1:3]
    lobs=[[_gmsh_lobatto(it,uvw[i]) for it in 2:p+1] for i in 1:3]
    prods=(lam[4]*lam[6],lam[2]*lam[6],lam[2]*lam[4],lam[1]*lam[6],
           lam[4]*lam[1],lam[3]*lam[6],lam[3]*lam[1],lam[3]*lam[2],
           lam[4]*lam[5],lam[5]*lam[2],lam[5]*lam[1],lam[5]*lam[3])
    dirs=((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),(0.0,1.0,0.0),
          (0.0,0.0,1.0),(1.0,0.0,0.0),(0.0,0.0,1.0),(0.0,0.0,1.0),
          (1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,1.0,0.0),(1.0,0.0,0.0))
    uvws=(1,2,3,2,3,1,3,3,1,2,2,1)
    e=_V3[]
    sizehint!(e,12*(p+1))
    for iedge in 1:12
        for k in 0:p
            push!(e,(legs[uvws[iedge]][k+1]*prods[iedge]) .* dirs[iedge])
        end
    end
    f=_V3[]
    sizehint!(f,12*p*(p+1))
    face_uvw=((6,1,2),(4,1,3),(2,2,3),(1,2,3),(3,1,3),(5,1,2))
    face_dir=((1.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),
              (0.0,1.0,0.0),(1.0,0.0,0.0),(1.0,0.0,0.0))
    face_dir2=((0.0,1.0,0.0),(0.0,0.0,1.0),(0.0,0.0,1.0),
               (0.0,0.0,1.0),(0.0,0.0,1.0),(0.0,1.0,0.0))
    for iface in 1:6
        idxl,v1,v2=face_uvw[iface]
        for i1 in 0:p
            for i2 in 0:p-1
                push!(f,(lam[idxl]*legs[v1][i1+1]*lobs[v2][i2+1]) .*
                         face_dir[iface])
            end
        end
        for i1 in 0:p-1
            for i2 in 0:p
                push!(f,(lam[idxl]*lobs[v1][i1+1]*legs[v2][i2+1]) .*
                         face_dir2[iface])
            end
        end
    end
    b=_V3[]
    sizehint!(b,3*p*p*(p+1))
    for i1 in 0:p
        for i2 in 0:p-1
            for i3 in 0:p-1
                push!(b,(legs[1][i1+1]*lobs[2][i2+1]*lobs[3][i3+1],
                         0.0,0.0))
            end
        end
    end
    for i1 in 0:p-1
        for i2 in 0:p
            for i3 in 0:p-1
                push!(b,(0.0,lobs[1][i1+1]*legs[2][i2+1]*lobs[3][i3+1],
                         0.0))
            end
        end
    end
    for i1 in 0:p-1
        for i2 in 0:p-1
            for i3 in 0:p
                push!(b,(0.0,0.0,
                         lobs[1][i1+1]*lobs[2][i2+1]*legs[3][i3+1]))
            end
        end
    end
    return e,f,b
end

function _curl_tables(::Val{:hex},p::Int,u,v,w)
    uvw=(u,v,w)
    lam=(0.5*(1.0+u),0.5*(1.0-u),0.5*(1.0+v),0.5*(1.0-v),
         0.5*(1.0+w),0.5*(1.0-w))
    dlam=(0.5,-0.5,0.5,-0.5,0.5,-0.5)
    legs=[[_gmsh_legendre(k,uvw[i]) for k in 0:p] for i in 1:3]
    lobs=[[_gmsh_lobatto(it,uvw[i]) for it in 2:p+1] for i in 1:3]
    dlobs=[[_gmsh_dlobatto(it,uvw[i]) for it in 2:p+1] for i in 1:3]
    cprod=((0.0,lam[4]*dlam[6],-dlam[4]*lam[6]),
           (-lam[2]*dlam[6],0.0,dlam[2]*lam[6]),
           (lam[2]*dlam[4],-dlam[2]*lam[4],0.0),
           (-lam[1]*dlam[6],0.0,dlam[1]*lam[6]),
           (lam[1]*dlam[4],-dlam[1]*lam[4],0.0),
           (0.0,lam[3]*dlam[6],-dlam[3]*lam[6]),
           (lam[1]*dlam[3],-dlam[1]*lam[3],0.0),
           (lam[2]*dlam[3],-dlam[2]*lam[3],0.0),
           (0.0,lam[4]*dlam[5],-dlam[4]*lam[5]),
           (-lam[2]*dlam[5],0.0,dlam[2]*lam[5]),
           (-lam[1]*dlam[5],0.0,dlam[1]*lam[5]),
           (0.0,lam[3]*dlam[5],-dlam[3]*lam[5]))
    uvws=(1,2,3,2,3,1,3,3,1,2,2,1)
    e=_V3[]
    sizehint!(e,12*(p+1))
    for iedge in 1:12
        for k in 0:p
            push!(e,(legs[uvws[iedge]][k+1]) .* cprod[iedge])
        end
    end
    # per face: (uvw1, uvw2, vec1, il1, jl1, vec2, il2, jl2)
    face_curl=(
        (1,2,(0.0,dlam[6],-lam[6]),3,2,(-dlam[6],0.0,lam[6]),3,1),
        (1,3,(0.0,lam[4],-dlam[4]),2,3,(dlam[4],-lam[4],0.0),2,1),
        (2,3,(-lam[2],0.0,dlam[2]),1,3,(lam[2],-dlam[2],0.0),1,2),
        (2,3,(-lam[1],0.0,dlam[1]),1,3,(lam[1],-dlam[1],0.0),1,2),
        (1,3,(0.0,lam[3],-dlam[3]),2,3,(dlam[3],-lam[3],0.0),2,1),
        (1,2,(0.0,dlam[5],-lam[5]),3,2,(-dlam[5],0.0,lam[5]),3,1))
    f=_V3[]
    sizehint!(f,12*p*(p+1))
    for iface in 1:6
        v1,v2,vec1,il1,jl1,vec2,il2,jl2=face_curl[iface]
        for i1 in 0:p
            for i2 in 0:p-1
                val=(0.0,0.0,0.0)
                lv=legs[v1][i1+1]
                val=(jl1==1 ? vec1[1]*lv*lobs[v2][i2+1] :
                     il1==1 ? vec1[1]*lv*dlobs[v2][i2+1] : 0.0,
                     jl1==2 ? vec1[2]*lv*lobs[v2][i2+1] :
                     il1==2 ? vec1[2]*lv*dlobs[v2][i2+1] : 0.0,
                     jl1==3 ? vec1[3]*lv*lobs[v2][i2+1] :
                     il1==3 ? vec1[3]*lv*dlobs[v2][i2+1] : 0.0)
                push!(f,val)
            end
        end
        for i1 in 0:p-1
            for i2 in 0:p
                lv=legs[v2][i2+1]
                val=(jl2==1 ? vec2[1]*lobs[v1][i1+1]*lv :
                     il2==1 ? vec2[1]*dlobs[v1][i1+1]*lv : 0.0,
                     jl2==2 ? vec2[2]*lobs[v1][i1+1]*lv :
                     il2==2 ? vec2[2]*dlobs[v1][i1+1]*lv : 0.0,
                     jl2==3 ? vec2[3]*lobs[v1][i1+1]*lv :
                     il2==3 ? vec2[3]*dlobs[v1][i1+1]*lv : 0.0)
                push!(f,val)
            end
        end
    end
    b=_V3[]
    sizehint!(b,3*p*p*(p+1))
    for i1 in 0:p
        for i2 in 0:p-1
            for i3 in 0:p-1
                push!(b,(0.0,
                         legs[1][i1+1]*lobs[2][i2+1]*dlobs[3][i3+1],
                         -legs[1][i1+1]*dlobs[2][i2+1]*lobs[3][i3+1]))
            end
        end
    end
    for i1 in 0:p-1
        for i2 in 0:p
            for i3 in 0:p-1
                push!(b,(-lobs[1][i1+1]*legs[2][i2+1]*dlobs[3][i3+1],
                         0.0,
                         dlobs[1][i1+1]*legs[2][i2+1]*lobs[3][i3+1]))
            end
        end
    end
    for i1 in 0:p-1
        for i2 in 0:p-1
            for i3 in 0:p
                push!(b,(lobs[1][i1+1]*dlobs[2][i2+1]*legs[3][i3+1],
                         -dlobs[1][i1+1]*lobs[2][i2+1]*legs[3][i3+1],
                         0.0))
            end
        end
    end
    return e,f,b
end

# --- Hcurl Prism -------------------------------------------------------------
# matrixVectorProductForMapping(a, v) = (2a*v1, 2a*v2, a*v3)
@inline function _pri_map(a::Float64,v::_V3)
    return (2.0*a*v[1],2.0*a*v[2],a*v[3])
end

function _hcurl_pri_terms(u::Float64,v::Float64,w::Float64)
    l1=v
    l2=1.0-u-v
    l3=u
    l4=0.5*(1.0+w)
    l5=0.5*(1.0-w)
    s05=sqrt(0.5)
    n1=(0.0,1.0,0.0)
    n2=(-s05,-s05,0.0)
    n3=(1.0,0.0,0.0)
    # 2-D dots: n2.t1=-s05, n3.t1=1, n3.t2=-1, n1.t2=1, n1.t3=1, n2.t3=-s05
    psie0=Vector{_V3}(undef,3)
    psie1=Vector{_V3}(undef,3)
    psie0[1]=ntuple(j->-l3*n2[j]/s05+l2*n3[j],3)
    psie0[2]=ntuple(j->-l1*n3[j]+l3*n1[j],3)
    psie0[3]=ntuple(j->l2*n1[j]-l1*n2[j]/s05,3)
    psie1[1]=ntuple(j->-l3*n2[j]/s05-l2*n3[j],3)
    psie1[2]=ntuple(j->-l1*n3[j]-l3*n1[j],3)
    psie1[3]=ntuple(j->-l2*n1[j]-l1*n2[j]/s05,3)
    subs=(l3-l2,l1-l3,l1-l2,l2-l1)
    return (l1,l2,l3,l4,l5),psie0,psie1,subs
end

function _hcurl_tables(::Val{:pri},p::Int,u,v,w)
    lam,psie0,psie1,subs=_hcurl_pri_terms(u,v,w)
    l1,l2,l3,l4,l5=lam
    legs=[[_gmsh_legendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:4]
    legw=[_gmsh_legendre(k,w) for k in 0:p]
    lobw=[_gmsh_lobatto(k,w) for k in 0:p+1]
    e=_V3[]
    sizehint!(e,9*(p+1))
    # horizontal edges: (edge index, psi index, lambda)
    for i in 1:9
        if i in (1,2,4,7,8,9)
            index=i in (1,7) ? 1 : (i in (2,8) ? 3 : 2)
            # C++ index: edges 0,6->0; 1,7->2; 3,8->1 (0-based)
            index=(i==1 || i==7) ? 1 : (i==2 || i==8) ? 3 : 2
            lambda=(i<=4) ? l5 : l4
            push!(e,_pri_map(lambda,psie0[index]))
            if p>=1
                push!(e,_pri_map(lambda,psie1[index]))
                for iedge in 2:p
                    val=Float64(Float32(2*iedge-1)/Float32(iedge))*legs[index][iedge] .*
                        psie1[index] .-
                        Float64(Float32(iedge-1)/Float32(iedge))*legs[index][iedge-1] .*
                        psie0[index]
                    push!(e,_pri_map(lambda,val))
                end
            end
        else
            lamb=i==3 ? l2 : (i==5 ? l3 : l1)
            for iedge in 0:p
                push!(e,_pri_map(1.0,(0.0,0.0,lamb*legw[iedge+1])))
            end
        end
    end
    f=_V3[]
    sizehint!(f,6*p*(p+1)+(p>=1 ? 6*(p-1)+2*(p-1)*(p-2) : 0))
    # quad faces
    for iFace in 1:3
        index=iFace==1 ? 1 : (iFace==2 ? 3 : 2)
        prod=iFace==1 ? l2*l3 : (iFace==2 ? l1*l2 : l1*l3)
        if p>0
            for n1 in 0:p
                facePsi=n1==0 ? psie0[index] :
                        n1==1 ? psie1[index] :
                        ntuple(j->Float64(Float32(2*n1-1)/Float32(n1))*legs[index][n1]*
                                  psie1[index][j]-
                                  Float64(Float32(n1-1)/Float32(n1))*legs[index][n1-1]*
                                  psie0[index][j],3)
                for n2 in 2:p+1
                    push!(f,_pri_map(lobw[n2+1],facePsi))
                end
            end
            for n1 in 2:p+1
                phie=prod*_gmsh_kernel(n1-2,subs[index])
                for n2 in 0:p
                    push!(f,_pri_map(1.0,(0.0,0.0,phie*legw[n2+1])))
                end
            end
        end
    end
    # tri faces (bottom iFace=0 uses lambda5, top uses lambda4)
    for iFace in 1:2
        lambda=iFace==1 ? l5 : l4
        tri_edge=((l2*l3,(0.0,0.5,0.0),1),
                  (l1*l3,(-0.5,-0.5,0.0),2),
                  (l1*l2,(0.5,0.0,0.0),4))
        for (prod,nd,idx) in tri_edge
            for n1 in 2:p
                push!(f,_pri_map(prod*legs[idx][n1-1]*lambda,nd))
            end
        end
        product=l1*l2*l3
        genuine=Float64[]
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                push!(genuine,0.5*lambda*product*legs[1][n1+1]*legs[4][n2+1])
                push!(f,_pri_map(1.0,(genuine[end],0.0,0.0)))
            end
        end
        for value in genuine
            push!(f,_pri_map(1.0,(0.0,value,0.0)))
        end
    end
    b=_V3[]
    sizehint!(b,3*max(p-1,0)*p+max(p-1,0)*max(p-2,0)*p+
              max(p-1,0)*p*(p+1)÷2)
    # quad-face-based bubbles
    qfb=((l2*l3,(0.0,0.5,0.0),1),(l1*l2,(0.5,0.0,0.0),4),
         (l1*l3,(-0.5,-0.5,0.0),2))
    for (prod,nd,index) in qfb
        for n1 in 2:p
            for n2 in 2:p+1
                push!(b,_pri_map(lobw[n2+1]*prod*legs[index][n1-1],nd))
            end
        end
    end
    # genuine bubbles
    product=l1*l2*l3
    genuine=_V3[]
    for n1 in 0:p-3
        for n2 in 0:p-3-n1
            term=(product*legs[1][n1+1]*legs[4][n2+1],0.0,0.0)
            for n3 in 2:p+1
                push!(genuine,_pri_map(lobw[n3+1],term))
                push!(b,genuine[end])
            end
        end
    end
    for value in genuine
        push!(b,(0.0,-value[1],0.0))
    end
    # tri-face-based bubbles
    for n1 in 0:p-2
        for n2 in 0:p-2-n1
            term=product*legs[1][n1+1]*legs[4][n2+1]
            for n3 in 0:p
                push!(b,(0.0,0.0,term*legw[n3+1]))
            end
        end
    end
    return e,f,b
end

function _curl_tables(::Val{:pri},p::Int,u,v,w)
    lam,psie0,psie1,subs=_hcurl_pri_terms(u,v,w)
    l1,l2,l3,l4,l5=lam
    dl4=0.5
    dl5=-0.5
    legs=[[_gmsh_legendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:4]
    dlegs=[[_gmsh_dlegendre(k,subs[i]) for k in 0:max(p-1,0)] for i in 1:4]
    legw=[_gmsh_legendre(k,w) for k in 0:p]
    lobw=[_gmsh_lobatto(k,w) for k in 0:p+1]
    dlobw=[_gmsh_dlobatto(k,w) for k in 0:p+1]
    curlpsie0=(1.0,1.0,-1.0)
    dsub=((1.0,0.5),(-0.5,0.5),(0.5,1.0),(-0.5,-1.0))
    @inline mvpm_curl(v)= (2.0*v[1],2.0*v[2],4.0*v[3])
    e=_V3[]
    sizehint!(e,9*(p+1))
    for i in 1:9
        if i in (1,2,4,7,8,9)
            index=(i==1 || i==7) ? 1 : (i==2 || i==8) ? 3 : 2
            lambda=(i<=4) ? l5 : l4
            dlambda=(i<=4) ? dl5 : dl4
            push!(e,mvpm_curl((-dlambda*psie0[index][2],
                               dlambda*psie0[index][1],
                               lambda*curlpsie0[index])))
            if p>=1
                push!(e,mvpm_curl((-dlambda*psie1[index][2],
                                   dlambda*psie1[index][1],0.0)))
                for iedge in 2:p
                    c1=Float64(Float32(2*iedge-1)/Float32(iedge))
                    c2=Float64(Float32(iedge-1)/Float32(iedge))
                    curlpsie=c1*(dsub[index][1]*dlegs[index][iedge]*
                                 psie1[index][2]-
                                 dsub[index][2]*dlegs[index][iedge]*
                                 psie1[index][1])-
                             c2*(curlpsie0[index]*legs[index][iedge-1]+
                                 dsub[index][1]*dlegs[index][iedge-1]*
                                 psie0[index][2]-
                                 dsub[index][2]*dlegs[index][iedge-1]*
                                 psie0[index][1])
                    psi1=c1*legs[index][iedge]*psie1[index][1]-
                         c2*legs[index][iedge-1]*psie0[index][1]
                    psi2=c1*legs[index][iedge]*psie1[index][2]-
                         c2*legs[index][iedge-1]*psie0[index][2]
                    push!(e,mvpm_curl((-dlambda*psi2,
                                       dlambda*psi1,
                                       lambda*curlpsie)))
                end
            end
        else
            dlambU=i==3 ? -0.5 : (i==5 ? 0.5 : 0.0)
            dlambV=i==3 ? -0.5 : (i==6 ? 0.5 : 0.0)
            for iedge in 0:p
                push!(e,mvpm_curl((dlambV*legw[iedge+1],
                                   -dlambU*legw[iedge+1],0.0)))
            end
        end
    end
    f=_V3[]
    sizehint!(f,6*p*(p+1)+(p>=1 ? 6*(p-1)+2*(p-1)*(p-2) : 0))
    # quad faces
    for iFace in 1:3
        index=iFace==1 ? 1 : (iFace==2 ? 3 : 2)
        prod=iFace==1 ? l2*l3 : (iFace==2 ? l1*l2 : l1*l3)
        dProd=iFace==1 ? (0.5*(l2-l3),-0.5*l3) :
              iFace==2 ? (-0.5*l1,0.5*(l2-l1)) : (0.5*l1,0.5*l3)
        if p>0
            for n1 in 0:p
                if n1==0
                    psie=(psie0[index][1],psie0[index][2])
                    curlpsie=curlpsie0[index]
                elseif n1==1
                    psie=(psie1[index][1],psie1[index][2])
                    curlpsie=0.0
                else
                    psie=ntuple(2) do j
                        Float64(Float32(2*n1-1)/Float32(n1))*legs[index][n1]*
                        psie1[index][j]-
                        Float64(Float32(n1-1)/Float32(n1))*legs[index][n1-1]*psie0[index][j]
                    end
                    curlpsie=Float64(Float32(2*n1-1)/Float32(n1))*
                             (dsub[index][1]*dlegs[index][n1]*
                              psie1[index][2]-
                              dsub[index][2]*dlegs[index][n1]*
                              psie1[index][1])-
                             Float64(Float32(n1-1)/Float32(n1))*
                             (curlpsie0[index]*legs[index][n1-1]+
                              dsub[index][1]*dlegs[index][n1-1]*
                              psie0[index][2]-
                              dsub[index][2]*dlegs[index][n1-1]*
                              psie0[index][1])
                end
                for n2 in 2:p+1
                    push!(f,mvpm_curl((-dlobw[n2+1]*psie[2],
                                       dlobw[n2+1]*psie[1],
                                       lobw[n2+1]*curlpsie)))
                end
            end
            for n1 in 2:p+1
                dphie=(dProd[1]*_gmsh_kernel(n1-2,subs[index])+
                       prod*dsub[index][1]*
                       _gmsh_dkernel(n1-2,subs[index]),
                       dProd[2]*_gmsh_kernel(n1-2,subs[index])+
                       prod*dsub[index][2]*
                       _gmsh_dkernel(n1-2,subs[index]))
                for n2 in 0:p
                    push!(f,mvpm_curl((dphie[2]*legw[n2+1],
                                       -dphie[1]*legw[n2+1],0.0)))
                end
            end
        end
    end
    # tri faces
    prod123=l1*l2*l3
    dl123U=0.5*l1*(l2-l3)
    dl123V=0.5*l3*(l2-l1)
    for iFace in 1:2
        lambda=iFace==1 ? l5 : l4
        dlambda=iFace==1 ? dl5 : dl4
        # edge-based groups
        prod=l2*l3
        dProdU=0.5*(l2-l3)
        for n1 in 2:p
            push!(f,mvpm_curl(
                      (-0.5*dlambda*prod*legs[1][n1-1],0.0,
                       0.5*lambda*(dProdU*legs[1][n1-1]+
                                   prod*dsub[1][1]*dlegs[1][n1-1]))))
        end
        prod=l1*l3
        dProdU=0.5*l1
        dProdV=0.5*l3
        for n1 in 2:p
            c0=0.5*dlambda*prod*legs[2][n1-1]
            push!(f,mvpm_curl(
                      (c0,-c0,
                       -0.5*lambda*(dProdU*legs[2][n1-1]+
                                    prod*dsub[2][1]*dlegs[2][n1-1]-
                                    (dProdV*legs[2][n1-1]+
                                     prod*dsub[2][2]*dlegs[2][n1-1])))))
        end
        prod=l1*l2
        dProdV=0.5*(l2-l1)
        for n1 in 2:p
            push!(f,mvpm_curl(
                      (0.0,0.5*dlambda*prod*legs[4][n1-1],
                       -0.5*(dProdV*legs[4][n1-1]+
                             prod*dsub[4][2]*dlegs[4][n1-1]))))
        end
        # genuine face functions
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                push!(f,mvpm_curl(
                          (0.0,
                           0.5*dlambda*prod123*legs[1][n1+1]*legs[4][n2+1],
                           -0.5*lambda*
                           (dl123V*legs[1][n1+1]*legs[4][n2+1]+
                            prod123*dsub[1][2]*dlegs[1][n1+1]*
                            legs[4][n2+1]+
                            prod123*dsub[4][2]*legs[1][n1+1]*
                            dlegs[4][n2+1]))))
            end
        end
        for n1 in 0:p-3
            for n2 in 0:p-3-n1
                push!(f,mvpm_curl(
                          (-0.5*dlambda*prod123*legs[1][n1+1]*legs[4][n2+1],
                           0.0,
                           0.5*lambda*
                           (dl123U*legs[1][n1+1]*legs[4][n2+1]+
                            prod123*dsub[1][1]*dlegs[1][n1+1]*
                            legs[4][n2+1]+
                            prod123*dsub[4][1]*legs[1][n1+1]*
                            dlegs[4][n2+1]))))
            end
        end
    end
    b=_V3[]
    sizehint!(b,3*max(p-1,0)*p+max(p-1,0)*max(p-2,0)*p+
              max(p-1,0)*p*(p+1)÷2)
    # quad-face-based bubbles
    for iFace in 1:3
        index=iFace==1 ? 1 : (iFace==2 ? 4 : 2)
        prod=iFace==1 ? l2*l3 : (iFace==2 ? l1*l2 : l1*l3)
        if iFace==1
            dProdU=0.5*(l2-l3)
            for n1 in 2:p
                phi=prod*legs[index][n1-1]
                dphiU=dProdU*legs[index][n1-1]+
                      prod*dsub[index][1]*dlegs[index][n1-1]
                for n2 in 2:p+1
                    push!(b,mvpm_curl((-0.5*dlobw[n2+1]*phi,0.0,
                                       0.5*lobw[n2+1]*dphiU)))
                end
            end
        elseif iFace==2
            dProdV=0.5*(l2-l1)
            for n1 in 2:p
                phi=prod*legs[index][n1-1]
                dphiV=dProdV*legs[index][n1-1]+
                      prod*dsub[index][2]*dlegs[index][n1-1]
                for n2 in 2:p+1
                    push!(b,mvpm_curl((0.0,0.5*dlobw[n2+1]*phi,
                                       -0.5*lobw[n2+1]*dphiV)))
                end
            end
        else
            dProdU=0.5*l1
            dProdV=0.5*l3
            for n1 in 2:p
                phi=prod*legs[index][n1-1]
                dphi=dProdU*legs[index][n1-1]+
                     prod*dsub[index][1]*dlegs[index][n1-1]-
                     (dProdV*legs[index][n1-1]+
                      prod*dsub[index][2]*dlegs[index][n1-1])
                for n2 in 2:p+1
                    push!(b,mvpm_curl((0.5*phi*dlobw[n2+1],
                                       -0.5*phi*dlobw[n2+1],
                                       -0.5*lobw[n2+1]*dphi)))
                end
            end
        end
    end
    # genuine bubbles
    for n1 in 0:p-3
        for n2 in 0:p-3-n1
            phi=prod123*legs[1][n1+1]*legs[4][n2+1]
            dphiV=-(dl123V*legs[1][n1+1]*legs[4][n2+1]+
                    prod123*dsub[1][2]*dlegs[1][n1+1]*legs[4][n2+1]+
                    prod123*dsub[4][2]*legs[1][n1+1]*dlegs[4][n2+1])
            for n3 in 2:p+1
                push!(b,mvpm_curl((0.0,dlobw[n3+1]*phi,
                                   lobw[n3+1]*dphiV)))
            end
        end
    end
    for n1 in 0:p-3
        for n2 in 0:p-3-n1
            phi=prod123*legs[1][n1+1]*legs[4][n2+1]
            dphiU=-(dl123U*legs[1][n1+1]*legs[4][n2+1]+
                    prod123*dsub[1][1]*dlegs[1][n1+1]*legs[4][n2+1]+
                    prod123*dsub[4][1]*legs[1][n1+1]*dlegs[4][n2+1])
            for n3 in 2:p+1
                push!(b,mvpm_curl((dlobw[n3+1]*phi,0.0,
                                   lobw[n3+1]*dphiU)))
            end
        end
    end
    # tri-face-based bubbles
    for n1 in 0:p-2
        for n2 in 0:p-2-n1
            dphiU=-(dl123U*legs[1][n1+1]*legs[4][n2+1]+
                    prod123*dsub[1][1]*dlegs[1][n1+1]*legs[4][n2+1]+
                    prod123*dsub[4][1]*legs[1][n1+1]*dlegs[4][n2+1])
            dphiV=dl123V*legs[1][n1+1]*legs[4][n2+1]+
                  prod123*dsub[1][2]*dlegs[1][n1+1]*legs[4][n2+1]+
                  prod123*dsub[4][2]*legs[1][n1+1]*dlegs[4][n2+1]
            for n3 in 0:p
                push!(b,mvpm_curl((dphiV*legw[n3+1],
                                   dphiU*legw[n3+1],0.0)))
            end
        end
    end
    return e,f,b
end

# ---------------------------------------------------------------------------
# Hcurl orientation machinery.
#
# Negative edge flag: functions at even local indices flip sign (all Hcurl
# families). Face orientation uses the same flag machinery as H1 but with the
# Hcurl-specific recompute formulas from the reference sources.
# ---------------------------------------------------------------------------
function _hcurl_negative_edge_table(e::Vector{_V3},edge_funcs::Int)
    out=copy(e)
    ne=length(e)÷edge_funcs
    for edge in 0:ne-1
        for k in 0:edge_funcs-1
            if k%2==0
                idx=edge*edge_funcs+k+1
                out[idx]=(-e[idx][1],-e[idx][2],-e[idx][3])
            end
        end
    end
    return out
end

# Shared flag-driven lambda/dlambda permutation for triangular faces.
function _permute_face_lams(lam,dlam,f1::Int,f2::Int)
    l=collect(lam)
    d=[collect(dlam[i]) for i in 1:3]
    if f1==1 && f2==-1
        l[1],l[2]=l[2],l[1]
        d[1],d[2]=d[2],d[1]
    elseif f1==0 && f2==-1
        l[2],l[3]=l[3],l[2]
        d[2],d[3]=d[3],d[2]
    elseif f1==2 && f2==-1
        l[1],l[3]=l[3],l[1]
        d[1],d[3]=d[3],d[1]
    elseif f1==1 && f2==1
        l[1],l[2],l[3]=l[2],l[3],l[1]
        d[1],d[2],d[3]=d[2],d[3],d[1]
    elseif f1==2 && f2==1
        l[1],l[2],l[3]=l[3],l[1],l[2]
        d[1],d[2],d[3]=d[3],d[1],d[2]
    end
    return l,d
end

# Oriented tri-face block via permuted lambdas, HcurlLegendre variant.
# jedge multiplies edge-based functions; jgenuine multiplies the genuine
# copy scalar; lam45 is the prism top/bottom lambda factor (1 otherwise).
function _hcurl_orient_tri_face(lam::NTuple{3,Float64},dlam::NTuple{3,_V3},
                                pf::Int,f1::Int,f2::Int,jedge::Float64,
                                jgenuine::Float64,lam45::Float64)
    l,d=_permute_face_lams(lam,dlam,f1,f2)
    sub=(l[2]-l[1],l[3]-l[2],l[1]-l[3])
    out=_V3[]
    sizehint!(out,3*(pf-1)+(pf-1)*(pf-2))
    for i in 1:3
        product2=i==1 ? l[2]*l[1] : (i==2 ? l[2]*l[3] : l[3]*l[1])
        normal=i==1 ? d[3] : (i==2 ? d[1] : d[2])
        for i1 in 2:pf
            value=jedge*product2*lam45*_gmsh_legendre(i1-2,sub[i])
            push!(out,(value*normal[1],value*normal[2],value*normal[3]))
        end
    end
    product=l[1]*l[2]*l[3]
    copies=Float64[]
    for n1_ in 0:pf-3
        ls1=_gmsh_legendre(n1_,sub[1])
        for n2_ in 0:pf-3-n1_
            push!(copies,
                  jgenuine*product*ls1*_gmsh_legendre(n2_,sub[3])*lam45)
            push!(out,copies[end] .* Tuple(d[2]))
        end
    end
    for value in copies
        push!(out,value .* Tuple(d[3]))
    end
    return out
end

# Oriented tri-face block, CurlHcurlLegendre variant. `full3d` selects the
# tetrahedron-style curlFunction form; otherwise the 2-D triangle/prism form
# is used with the prism's detjacob/det/lambda45/dlambda45 constants.
function _curl_orient_tri_face(lam::NTuple{3,Float64},dlam::NTuple{3,_V3},
                               pf::Int,f1::Int,f2::Int,full3d::Bool,
                               detjacob::Float64,det::Float64,
                               lam45::Float64,dlam45::Float64,
                               preperm::Bool=false)
    l,d=_permute_face_lams(lam,dlam,f1,f2)
    sub=(l[2]-l[1],l[3]-l[2],l[1]-l[3])
    dsub=(ntuple(j->d[2][j]-d[1][j],3),
          ntuple(j->d[3][j]-d[2][j],3),
          ntuple(j->d[1][j]-d[3][j],3))
    out=_V3[]
    sizehint!(out,3*(pf-1)+(pf-1)*(pf-2))
    for i in 1:3
        # C++ pairs: i=0 -> (lambda[0],lambda[1]); i=1 -> (lambda[2],lambda[1]);
        # i=2 -> (lambda[0],lambda[2])
        ia,ib=i==1 ? (1,2) : (i==2 ? (3,2) : (1,3))
        prod_=l[ia]*l[ib]
        dprod=ntuple(j->d[ia][j]*l[ib]+d[ib][j]*l[ia],3)
        normal=i==1 ? d[3] : (i==2 ? d[1] : d[2])
        for i1 in 2:pf
            leg=_gmsh_legendre(i1-2,sub[i])
            dleg=_gmsh_dlegendre(i1-2,sub[i])
            dphi=ntuple(j->dprod[j]*leg+prod_*dsub[i][j]*dleg,3)
            if full3d
                push!(out,_curl_func(detjacob,Tuple(normal),dphi))
            else
                c0=-detjacob*normal[2]*dlam45*prod_*leg
                c1=detjacob*normal[1]*dlam45*prod_*leg
                c2=det*lam45*(normal[2]*dphi[1]-normal[1]*dphi[2])
                push!(out,(c0,c1,c2))
            end
        end
    end
    subBA=sub[1]
    subAC=sub[3]
    dsubBA=dsub[1]
    dsubAC=dsub[3]
    lsub_ac=[_gmsh_legendre(k,subAC) for k in 0:max(pf-3,0)]
    dlsub_ac=[_gmsh_dlegendre(k,subAC) for k in 0:max(pf-3,0)]
    lsub_ba=[_gmsh_legendre(k,subBA) for k in 0:max(pf-3,0)]
    dlsub_ba=[_gmsh_dlegendre(k,subBA) for k in 0:max(pf-3,0)]
    dprod=preperm ?
        (0.5*lam[3]*(lam[1]-lam[2]),0.5*lam[2]*(lam[1]-lam[3]),0.0) :
        ntuple(3) do i
            d[1][i]*l[2]*l[3]+d[2][i]*l[1]*l[3]+d[3][i]*l[2]*l[1]
        end
    product=l[1]*l[2]*l[3]
    copies=_V3[]
    for n1_ in 0:pf-3
        for n2_ in 0:pf-3-n1_
            gface=ntuple(3) do i
                dprod[i]*lsub_ac[n2_+1]*lsub_ba[n1_+1]+
                product*dsubBA[i]*lsub_ac[n2_+1]*dlsub_ba[n1_+1]+
                product*dsubAC[i]*dlsub_ac[n2_+1]*lsub_ba[n1_+1]
            end
            push!(copies,gface)
            if full3d
                push!(out,_curl_func(detjacob,Tuple(d[2]),gface))
            else
                c0=-detjacob*d[2][2]*dlam45*product*lsub_ba[n1_+1]*
                    lsub_ac[n2_+1]
                c1=detjacob*d[2][1]*dlam45*product*lsub_ba[n1_+1]*
                    lsub_ac[n2_+1]
                c2=det*lam45*(gface[1]*d[2][2]-d[2][1]*gface[2])
                push!(out,(c0,c1,c2))
            end
        end
    end
    idx=1
    for n1_ in 0:pf-3
        for n2_ in 0:pf-3-n1_
            gface=copies[idx]
            idx+=1
            if full3d
                push!(out,_curl_func(detjacob,Tuple(d[3]),gface))
            else
                c0=-detjacob*d[3][2]*dlam45*product*lsub_ba[n1_+1]*
                    lsub_ac[n2_+1]
                c1=detjacob*d[3][1]*dlam45*product*lsub_ba[n1_+1]*
                    lsub_ac[n2_+1]
                c2=det*lam45*(gface[1]*d[3][2]-d[3][1]*gface[2])
                push!(out,(c0,c1,c2))
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Hcurl face descriptors and orientation drivers
# ---------------------------------------------------------------------------
# Tri-face descriptors: (lam tuple, dlam tuple, jedge, jgenuine, lam45,
# detjacob, det, dlam45, full3d)
function _hcurl_tri_face_desc(::Val{:tri},u::Float64,v::Float64)
    l1=v
    l2=1.0-u-v
    l3=u
    return ((l2,l3,l1),
            ((-0.5,-0.5,0.0),(0.5,0.0,0.0),(0.0,0.5,0.0)),
            2.0,2.0,1.0,0.0,4.0,0.0,false)
end

function _hcurl_tri_face_desc(::Val{:tet},face::Int,u::Float64,v::Float64,
                              w::Float64)
    l1=v
    l2=1.0-u-v-w
    l3=u
    l4=w
    lams=(l1,l2,l3,l4)
    dlams=((0.0,0.5,0.0),(-0.5,-0.5,-0.5),(0.5,0.0,0.0),(0.0,0.0,0.5))
    picks=((2,3,1),(2,3,4),(2,1,4),(3,1,4))
    pk=picks[face]
    return ((lams[pk[1]],lams[pk[2]],lams[pk[3]]),
            (dlams[pk[1]],dlams[pk[2]],dlams[pk[3]]),
            2.0,2.0,1.0,4.0,4.0,0.0,true)
end

function _hcurl_tri_face_desc(::Val{:pri},face::Int,u::Float64,v::Float64,
                              w::Float64)
    l1=v
    l2=1.0-u-v
    l3=u
    l4=0.5*(1.0+w)
    l5=0.5*(1.0-w)
    lam45=face==1 ? l5 : l4
    dlam45=face==1 ? -0.5 : 0.5
    return ((l2,l3,l1),
            ((-0.5,-0.5,0.0),(0.5,0.0,0.0),(0.0,0.5,0.0)),
            2.0,2.0,lam45,2.0,4.0,dlam45,false)
end

# Quad-face oriented block (HcurlQuad/HcurlBrick style).
# lam_factor: the face's normal lambda; v1,v2: coordinate components (1=u,2=v,
# 3=w); dir1/dir2: unit direction tuples; flag triple (f1,f2,f3).
# For the flag3==-1 branch the groups are recomputed in transposed order.
function _hcurl_orient_quad_face(canonical::Vector{_V3},p::Int,
                                 uvw::NTuple{3,Float64},lam_factor::Float64,
                                 v1::Int,v2::Int,dir1::_V3,dir2::_V3,
                                 f1::Int,f2::Int,f3::Int)
    if f3==1
        out=copy(canonical)
        it=1
        for it1 in 0:p
            for it2 in 0:p-1
                i1=(f1==-1 && it1%2==0) ? -1.0 : 1.0
                i2=(f2==-1 && it2%2!=0) ? -1.0 : 1.0
                out[it]=out[it] .* (i1*i2)
                it+=1
            end
        end
        for it1 in 0:p-1
            for it2 in 0:p
                i1=(f1==-1 && it1%2!=0) ? -1.0 : 1.0
                i2=(f2==-1 && it2%2==0) ? -1.0 : 1.0
                out[it]=out[it] .* (i1*i2)
                it+=1
            end
        end
        return out
    end
    lkv1=[_gmsh_lobatto(it,uvw[v2]) for it in 2:p+1]
    lkv2=[_gmsh_lobatto(it,uvw[v1]) for it in 2:p+1]
    leg1=[_gmsh_legendre(it,uvw[v1]) for it in 0:p]
    leg2=[_gmsh_legendre(it,uvw[v2]) for it in 0:p]
    out=_V3[]
    sizehint!(out,2*p*(p+1))
    for it1 in 0:p
        for it2 in 0:p-1
            i1=(f2==-1 && it1%2==0) ? -1.0 : 1.0
            i2=(f1==-1 && it2%2!=0) ? -1.0 : 1.0
            push!(out,(lam_factor*leg2[it1+1]*lkv2[it2+1]*i1*i2) .* dir2)
        end
    end
    for it1 in 0:p-1
        for it2 in 0:p
            i1=(f2==-1 && it1%2!=0) ? -1.0 : 1.0
            i2=(f1==-1 && it2%2==0) ? -1.0 : 1.0
            push!(out,(lam_factor*leg1[it2+1]*lkv1[it1+1]*i1*i2) .* dir1)
        end
    end
    return out
end

# Curl variant of the oriented quad-face block. Per face, vec1/vec2 carry the
# lambda/dlambda coefficients and (il,jl) select which component multiplies
# the differentiated versus plain Lobatto factor.
function _curl_orient_quad_face(p::Int,uvw::NTuple{3,Float64},
                                vec1::_V3,il1::Int,jl1::Int,
                                vec2::_V3,il2::Int,jl2::Int,
                                v1::Int,v2::Int,
                                f1::Int,f2::Int,f3::Int,
                                canonical::Vector{_V3})
    if f3==1
        out=copy(canonical)
        it=1
        for it1 in 0:p
            for it2 in 0:p-1
                i1=(f1==-1 && it1%2==0) ? -1.0 : 1.0
                i2=(f2==-1 && it2%2!=0) ? -1.0 : 1.0
                out[it]=out[it] .* (i1*i2)
                it+=1
            end
        end
        for it1 in 0:p-1
            for it2 in 0:p
                i1=(f1==-1 && it1%2!=0) ? -1.0 : 1.0
                i2=(f2==-1 && it2%2==0) ? -1.0 : 1.0
                out[it]=out[it] .* (i1*i2)
                it+=1
            end
        end
        return out
    end
    lkv1=[_gmsh_lobatto(it,uvw[v2]) for it in 2:p+1]
    lkv2=[_gmsh_lobatto(it,uvw[v1]) for it in 2:p+1]
    dlkv1=[_gmsh_dlobatto(it,uvw[v2]) for it in 2:p+1]
    dlkv2=[_gmsh_dlobatto(it,uvw[v1]) for it in 2:p+1]
    leg1=[_gmsh_legendre(it,uvw[v1]) for it in 0:p]
    leg2=[_gmsh_legendre(it,uvw[v2]) for it in 0:p]
    out=_V3[]
    sizehint!(out,2*p*(p+1))
    for it1 in 0:p
        for it2 in 0:p-1
            i1=(f2==-1 && it1%2==0) ? -1.0 : 1.0
            i2=(f1==-1 && it2%2!=0) ? -1.0 : 1.0
            val=ntuple(3) do j
                if j==jl2
                    vec2[j]*lkv2[it2+1]*leg2[it1+1]*i1*i2
                elseif j==il2
                    vec2[j]*dlkv2[it2+1]*leg2[it1+1]*i1*i2
                else
                    0.0
                end
            end
            push!(out,val)
        end
    end
    for it1 in 0:p-1
        for it2 in 0:p
            i1=(f2==-1 && it1%2!=0) ? -1.0 : 1.0
            i2=(f1==-1 && it2%2==0) ? -1.0 : 1.0
            val=ntuple(3) do j
                if j==jl1
                    vec1[j]*leg1[it2+1]*lkv1[it1+1]*i1*i2
                elseif j==il1
                    vec1[j]*leg1[it2+1]*dlkv1[it1+1]*i1*i2
                else
                    0.0
                end
            end
            push!(out,val)
        end
    end
    return out
end

# Quad-face descriptors for the oriented path:
# (lam_index into hex lambdas, v1, v2, dir1, dir2) for the value table and
# (vec1,il1,jl1,vec2,il2,jl2,v1,v2) for the curl table.
const _HCURL_HEX_FACE_DESC=
    ((6,1,2,(1.0,0.0,0.0),(0.0,1.0,0.0)),
     (4,1,3,(1.0,0.0,0.0),(0.0,0.0,1.0)),
     (2,2,3,(0.0,1.0,0.0),(0.0,0.0,1.0)),
     (1,2,3,(0.0,1.0,0.0),(0.0,0.0,1.0)),
     (3,1,3,(1.0,0.0,0.0),(0.0,0.0,1.0)),
     (5,1,2,(1.0,0.0,0.0),(0.0,1.0,0.0)))

function _curl_hex_face_desc(face::Int,lam,dlam)
    l6,l5,l4,l3=lam[6],lam[5],lam[4],lam[3]
    l2v,l1v=lam[2],lam[1]
    d6,d5,d4,d3=dlam[6],dlam[5],dlam[4],dlam[3]
    d2,d1=dlam[2],dlam[1]
    if face==1
        return ((0.0,d6,-l6),3,2,(-d6,0.0,l6),3,1,1,2)
    elseif face==2
        return ((0.0,l4,-d4),2,3,(d4,-l4,0.0),2,1,1,3)
    elseif face==3
        return ((-l2v,0.0,d2),1,3,(l2v,-d2,0.0),1,2,2,3)
    elseif face==4
        return ((-l1v,0.0,d1),1,3,(l1v,-d1,0.0),1,2,2,3)
    elseif face==5
        return ((0.0,l3,-d3),2,3,(d3,-l3,0.0),2,1,1,3)
    else
        return ((0.0,d5,-l5),3,2,(-d5,0.0,l5),3,1,1,2)
    end
end

# Prism quad-face oriented recompute (flag3==-1). Per-face 2-D Whitney
# vectors and prod/sub pairs from the pinned source; jacob=2 maps the
# transverse components.
function _hcurl_pri_orient_quad_face(p::Int,face::Int,u::Float64,
                                     v::Float64,w::Float64,f1::Int,f2::Int)
    l1=v
    l2=1.0-u-v
    l3=u
    if face==1
        prod=l2*l3
        sub=l3-l2
        psie0=(l3+l2,l3)
        psie1=(l3-l2,l3)
    elseif face==2
        prod=l1*l2
        sub=l1-l2
        psie0=(l1,l1+l2)
        psie1=(l1,l1-l2)
    else
        prod=l1*l3
        sub=l1-l3
        psie0=(-l1,l3)
        psie1=(-l1,-l3)
    end
    jacob=2.0
    lsub=[_gmsh_legendre(k,sub) for k in 0:max(p-1,0)]
    phi=[prod*_gmsh_kernel(k,sub) for k in 0:max(p-1,0)]
    out=_V3[]
    sizehint!(out,2*p*(p+1))
    for it1 in 0:p
        lw=_gmsh_legendre(it1,w)
        i1=(f2==-1 && it1%2==0) ? -1.0 : 1.0
        for it2 in 2:p+1
            i2=(f1==-1 && it2%2!=0) ? -1.0 : 1.0
            push!(out,(0.0,0.0,i1*i2*phi[it2-1]*lw))
        end
    end
    for it1 in 2:p+1
        lw=_gmsh_lobatto(it1,w)
        i1=(f2==-1 && it1%2!=0) ? -1.0 : 1.0
        if f1==-1
            push!(out,(-jacob*i1*lw*psie0[1],-jacob*i1*lw*psie0[2],0.0))
            push!(out,(jacob*i1*lw*psie1[1],jacob*i1*lw*psie1[2],0.0))
        else
            push!(out,(jacob*i1*lw*psie0[1],jacob*i1*lw*psie0[2],0.0))
            push!(out,(jacob*i1*lw*psie1[1],jacob*i1*lw*psie1[2],0.0))
        end
        for it2 in 2:p
            i2=(f1==-1 && it2%2==0) ? -1.0 : 1.0
            c1=Float64(Float32(2*it2-1)/Float32(it2))
            c2=Float64(Float32(it2-1)/Float32(it2))
            px=i1*i2*lw*jacob*(c1*lsub[it2]*psie1[1]-c2*lsub[it2-1]*psie0[1])
            py=i1*i2*lw*jacob*(c1*lsub[it2]*psie1[2]-c2*lsub[it2-1]*psie0[2])
            push!(out,(px,py,0.0))
        end
    end
    return out
end

# Curl variant of the prism quad-face oriented recompute.
function _curl_pri_orient_quad_face(p::Int,face::Int,u::Float64,
                                    v::Float64,w::Float64,f1::Int,f2::Int)
    l1=v
    l2=1.0-u-v
    l3=u
    detjacob=2.0
    det=4.0
    if face==1
        prod=l2*l3
        sub=l3-l2
        psie0=(l3+l2,l3)
        psie1=(l3-l2,l3)
        dprod=(0.5*(l2-l3),-0.5*l3)
        dsub=(1.0,0.5)
        cpsie0=1.0
    elseif face==2
        prod=l1*l2
        sub=l1-l2
        psie0=(l1,l1+l2)
        psie1=(l1,l1-l2)
        dprod=(-0.5*l1,0.5*(l2-l1))
        dsub=(0.5,1.0)
        cpsie0=-1.0
    else
        prod=l1*l3
        sub=l1-l3
        psie0=(-l1,l3)
        psie1=(-l1,-l3)
        dprod=(0.5*l1,0.5*l3)
        dsub=(-0.5,0.5)
        cpsie0=1.0
    end
    lsub=[_gmsh_legendre(k,sub) for k in 0:max(p-1,0)]
    dlsub=[_gmsh_dlegendre(k,sub) for k in 0:max(p-1,0)]
    phi=[_gmsh_kernel(k,sub) for k in 0:max(p-1,0)]
    dphi=[_gmsh_dkernel(k,sub) for k in 0:max(p-1,0)]
    out=_V3[]
    sizehint!(out,2*p*(p+1))
    for it1 in 0:p
        lw=_gmsh_legendre(it1,w)
        i1=(f2==-1 && it1%2==0) ? -1.0 : 1.0
        for it2 in 2:p+1
            i2=(f1==-1 && it2%2!=0) ? -1.0 : 1.0
            dphie=(dprod[1]*phi[it2-1]+prod*dsub[1]*dphi[it2-1],
                   dprod[2]*phi[it2-1]+prod*dsub[2]*dphi[it2-1])
            push!(out,(detjacob*i1*i2*lw*dphie[2],
                       -detjacob*i1*i2*lw*dphie[1],0.0))
        end
    end
    for it1 in 2:p+1
        lw=_gmsh_lobatto(it1,w)
        dlw=_gmsh_dlobatto(it1,w)
        i1=(f2==-1 && it1%2!=0) ? -1.0 : 1.0
        if f1==-1
            push!(out,(detjacob*i1*dlw*psie0[2],
                       -detjacob*i1*dlw*psie0[1],
                       -det*lw*cpsie0*i1))
            push!(out,(-detjacob*i1*dlw*psie1[2],
                       detjacob*i1*dlw*psie1[1],0.0))
        else
            push!(out,(-detjacob*i1*dlw*psie0[2],
                       detjacob*i1*dlw*psie0[1],
                       det*lw*cpsie0*i1))
            push!(out,(-detjacob*i1*dlw*psie1[2],
                       detjacob*i1*dlw*psie1[1],0.0))
        end
        for it2 in 2:p
            i2=(f1==-1 && it2%2==0) ? -1.0 : 1.0
            c1=Float64(Float32(2*it2-1)/Float32(it2))
            c2=Float64(Float32(it2-1)/Float32(it2))
            psie=(c1*lsub[it2]*psie1[1]-c2*lsub[it2-1]*psie0[1],
                  c1*lsub[it2]*psie1[2]-c2*lsub[it2-1]*psie0[2])
            cpsie=c1*(dsub[1]*dlsub[it2]*psie1[2]-
                      dsub[2]*dlsub[it2]*psie1[1])-
                  c2*(cpsie0*lsub[it2-1]+
                      dsub[1]*dlsub[it2-1]*psie0[2]-
                      dsub[2]*dlsub[it2-1]*psie0[1])
            push!(out,(-i1*i2*dlw*detjacob*psie[2],
                       i1*i2*dlw*detjacob*psie[1],
                       i1*i2*lw*det*cpsie))
        end
    end
    return out
end

# --- Hcurl drivers -----------------------------------------------------------
const _HCURL_FAMILIES=(:lin,:tri,:qua,:tet,:hex,:pri)

# Tri-face flag pairs in next_permutation order over sorted vertices.
function _hcurl_orient_one_face(family::Symbol,p::Int,u::Float64,v::Float64,
                                w::Float64,face::Int,isquad::Bool,
                                f1::Int,f2::Int,f3::Int,curl::Bool,
                                canonical::Vector{_V3})
    if isquad
        if family===:qua
            if curl
                if f3==1
                    out=copy(canonical)
                    it=1
                    for it1 in 0:p
                        for it2 in 0:p-1
                            i1=(f1==-1 && it1%2==0) ? -1.0 : 1.0
                            i2=(f2==-1 && it2%2!=0) ? -1.0 : 1.0
                            out[it]=out[it].*(i1*i2)
                            it+=1
                        end
                    end
                    for it1 in 0:p-1
                        for it2 in 0:p
                            i1=(f1==-1 && it1%2!=0) ? -1.0 : 1.0
                            i2=(f2==-1 && it2%2==0) ? -1.0 : 1.0
                            out[it]=out[it].*(i1*i2)
                            it+=1
                        end
                    end
                    return out
                end
                legv=[_gmsh_legendre(k,v) for k in 0:p]
                legu=[_gmsh_legendre(k,u) for k in 0:p]
                out=_V3[]
                sizehint!(out,2*p*(p+1))
                for it1 in 0:p
                    for it2 in 2:p+1
                        i1=(f2==-1 && it1%2==0) ? -1.0 : 1.0
                        i2=(f1==-1 && it2%2!=0) ? -1.0 : 1.0
                        push!(out,(0.0,0.0,
                                   legv[it1+1]*_gmsh_dlobatto(it2,u)*i1*i2))
                    end
                end
                for it1 in 2:p+1
                    for it2 in 0:p
                        i1=(f2==-1 && it1%2!=0) ? -1.0 : 1.0
                        i2=(f1==-1 && it2%2==0) ? -1.0 : 1.0
                        push!(out,(0.0,0.0,
                                   -legu[it2+1]*_gmsh_dlobatto(it1,v)*i1*i2))
                    end
                end
                return out
            end
            return _hcurl_orient_quad_face(canonical,p,(u,v,w),1.0,1,2,
                                           (1.0,0.0,0.0),(0.0,1.0,0.0),
                                           f1,f2,f3)
        elseif family===:hex
            lam=(0.5*(1.0+u),0.5*(1.0-u),0.5*(1.0+v),0.5*(1.0-v),
                 0.5*(1.0+w),0.5*(1.0-w))
            dlam=(0.5,-0.5,0.5,-0.5,0.5,-0.5)
            desc=_HCURL_HEX_FACE_DESC[face]
            if curl
                cdesc=_curl_hex_face_desc(face,lam,dlam)
                return _curl_orient_quad_face(p,(u,v,w),cdesc[1],cdesc[2],
                                              cdesc[3],cdesc[4],cdesc[5],
                                              cdesc[6],cdesc[7],cdesc[8],
                                              f1,f2,f3,canonical)
            end
            return _hcurl_orient_quad_face(canonical,p,(u,v,w),
                                           lam[desc[1]],desc[2],desc[3],
                                           desc[4],desc[5],f1,f2,f3)
        elseif family===:pri
            if f3==1
                # canonical block with sign flips
                out=copy(canonical)
                it=1
                for it1 in 0:p
                    for it2 in 2:p+1
                        i1=(f1==-1 && it1%2==0) ? -1.0 : 1.0
                        i2=(f2==-1 && it2%2!=0) ? -1.0 : 1.0
                        out[it]=out[it].*(i1*i2)
                        it+=1
                    end
                end
                for it1 in 2:p+1
                    for it2 in 0:p
                        i1=(f1==-1 && it1%2!=0) ? -1.0 : 1.0
                        i2=(f2==-1 && it2%2==0) ? -1.0 : 1.0
                        out[it]=out[it].*(i1*i2)
                        it+=1
                    end
                end
                return out
            end
            if curl
                return _curl_pri_orient_quad_face(p,face,u,v,w,f1,f2)
            end
            return _hcurl_pri_orient_quad_face(p,face,u,v,w,f1,f2)
        end
        error("no quad faces for Hcurl family $family")
    else
        # triangular face
        if family===:tri
            desc=_hcurl_tri_face_desc(Val(:tri),u,v)
        elseif family===:tet
            desc=_hcurl_tri_face_desc(Val(:tet),face,u,v,w)
        elseif family===:pri
            desc=_hcurl_tri_face_desc(Val(:pri),face-3,u,v,w)
        else
            error("no tri faces for Hcurl family $family")
        end
        lam,dlam,je,jg,lam45,dj,det,dl45,full3d=desc
        if curl
            return _curl_orient_tri_face(lam,dlam,p,f1,f2,full3d,dj,det,
                                         lam45,dl45,family===:tri)
        end
        return _hcurl_orient_tri_face(lam,dlam,p,f1,f2,je,jg,lam45)
    end
end

# ---------------------------------------------------------------------------
# Hcurl getKeysInfo / typeKeys / top-level evaluation
# ---------------------------------------------------------------------------
function _hcurl_keys_info(family::Symbol,p::Int)
    counts=_hcurl_counts(family,p)
    info=Vector{Tuple{Int32,Int32}}(undef,
        counts.vertex+counts.edge+counts.face+counts.bubble)
    it=1
    nedge=length(get(_SOLIN_EDGES,family,()))
    for _ in 1:nedge, i in 0:p
        info[it]=(Int32(1),Int32(i))
        it+=1
    end
    if family===:tri
        for _ in 1:3, i1 in 2:p
            info[it]=(Int32(2),Int32(i1))
            it+=1
        end
        for _ in 1:2, n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0)
            info[it]=(Int32(2),Int32(n1+n2+1))
            it+=1
        end
    elseif family===:qua
        for _ in 1:1
            for n1 in 0:p, n2 in 2:p+1
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
            for n1 in 2:p+1, n2 in 0:p
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
        end
    elseif family===:tet
        for _ in 1:4
            for _ in 1:3, i1 in 2:p
                info[it]=(Int32(2),Int32(i1))
                it+=1
            end
            for _ in 1:2, n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0)
                info[it]=(Int32(2),Int32(n1+n2+1))
                it+=1
            end
        end
        for _ in 1:4, n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0)
            info[it]=(Int32(3),Int32(n1+n2+1))
            it+=1
        end
        for _ in 1:3, n1 in 1:max(p-3,0), n2 in 1:max(p-2-n1,0),
            n3 in 1:max(p-1-n1-n2,0)
            info[it]=(Int32(3),Int32(n1+n2+n3+1))
            it+=1
        end
    elseif family===:hex
        for _ in 1:6
            for n1 in 0:p, n2 in 2:p+1
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
            for n1 in 2:p+1, n2 in 0:p
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
        end
        for n1 in 0:p, n2 in 2:p+1, n3 in 2:p+1
            info[it]=(Int32(3),Int32(max(n1,n2,n3)))
            it+=1
        end
        for n1 in 2:p+1, n2 in 0:p, n3 in 2:p+1
            info[it]=(Int32(3),Int32(max(n1,n2,n3)))
            it+=1
        end
        for n1 in 2:p+1, n2 in 2:p+1, n3 in 0:p
            info[it]=(Int32(3),Int32(max(n1,n2,n3)))
            it+=1
        end
    elseif family===:pri
        for _ in 1:3
            for n1 in 0:p, n2 in 2:p+1
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
            for n1 in 2:p+1, n2 in 0:p
                info[it]=(Int32(2),Int32(max(n1,n2)))
                it+=1
            end
        end
        for _ in 1:2
            for _ in 1:3, n1 in 2:p
                info[it]=(Int32(2),Int32(n1))
                it+=1
            end
            for _ in 1:2, n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0)
                info[it]=(Int32(2),Int32(n1+n2+1))
                it+=1
            end
        end
        for _ in 1:3, n1 in 2:p, n2 in 2:p+1
            info[it]=(Int32(3),Int32(max(n1,n2)))
            it+=1
        end
        for _ in 1:2, n1 in 1:max(p-2,0), n2 in 1:max(p-1-n1,0),
            n3 in 2:p+1
            info[it]=(Int32(3),Int32(max(n1+n2+1,n3)))
            it+=1
        end
        for n1 in 1:max(p-1,0), n2 in 1:max(p-n1,0), n3 in 0:p
            info[it]=(Int32(3),Int32(max(n1+n2+1,n3)))
            it+=1
        end
    end
    it==length(info)+1 || error(
        "HierarchicalBases: internal Hcurl keysInfo length mismatch")
    return info
end

function _hcurl_type_keys(family::Symbol,p::Int)
    counts=_hcurl_counts(family,p)
    per_edge=p+1
    blocks=_hcurl_face_blocks(family,p)
    per_quad=0
    per_tri=0
    for block in blocks
        if block[3]===:quad
            per_quad=block[2]
        else
            per_tri=block[2]
        end
    end
    const1=per_edge+1
    const2=const1+per_quad
    const3=const1+per_tri
    const4=counts.bubble+max(const2,const3)
    keys=Vector{Int32}(undef,
        counts.vertex+counts.edge+counts.face+counts.bubble)
    it=1
    for _ in 1:counts.vertex
        keys[it]=Int32(0)
        it+=1
    end
    for _ in 1:length(get(_SOLIN_EDGES,family,())), k in 1:per_edge
        keys[it]=Int32(k)
        it+=1
    end
    for block in blocks
        limit=block[3]===:quad ? const2 : const3
        for k in const1:limit-1
            keys[it]=Int32(k)
            it+=1
        end
    end
    for k in max(const3,const2):const4-1
        keys[it]=Int32(k)
        it+=1
    end
    it==length(keys)+1 || error(
        "HierarchicalBases: internal Hcurl typeKeys length mismatch")
    return keys
end

@inline function _hcurl_tables(family::Symbol,p::Int,u,v,w)
    family===:lin && return _hcurl_tables(Val(:lin),p,u,v,w)
    family===:tri && return _hcurl_tables(Val(:tri),p,u,v,w)
    family===:qua && return _hcurl_tables(Val(:qua),p,u,v,w)
    family===:tet && return _hcurl_tables(Val(:tet),p,u,v,w)
    family===:hex && return _hcurl_tables(Val(:hex),p,u,v,w)
    return _hcurl_tables(Val(:pri),p,u,v,w)
end

@inline function _curl_tables(family::Symbol,p::Int,u,v,w)
    family===:lin && return _curl_tables(Val(:lin),p,u,v,w)
    family===:tri && return _curl_tables(Val(:tri),p,u,v,w)
    family===:qua && return _curl_tables(Val(:qua),p,u,v,w)
    family===:tet && return _curl_tables(Val(:tet),p,u,v,w)
    family===:hex && return _curl_tables(Val(:hex),p,u,v,w)
    return _curl_tables(Val(:pri),p,u,v,w)
end

function _hcurl_all_oriented_faces(family::Symbol,p::Int,
                                   u::Float64,v::Float64,w::Float64,
                                   f::Vector{_V3},curl::Bool)
    blocks=_hcurl_face_blocks(family,p)
    quad_total=0
    tri_total=0
    for block in blocks
        block[3]===:quad ? (quad_total+=block[2]) : (tri_total+=block[2])
    end
    quad_all=Vector{Vector{_V3}}(undef,8)
    for o in 1:8
        work=copy(f)
        f1,f2,f3=_QUAD_ORIENTATION_FLAGS[o]
        for (face,block) in enumerate(blocks)
            block[3]===:quad || continue
            offset,count,_=block
            work[offset+1:offset+count]=_hcurl_orient_one_face(
                family,p,u,v,w,face,true,f1,f2,f3,curl,
                f[offset+1:offset+count])
        end
        quad_all[o]=work[1:quad_total]
    end
    tri_all=Vector{Vector{_V3}}(undef,6)
    for o in 1:6
        work=copy(f)
        f1,f2=_TRI_ORIENTATION_FLAGS[o]
        for (face,block) in enumerate(blocks)
            block[3]===:tri || continue
            offset,count,_=block
            work[offset+1:offset+count]=_hcurl_orient_one_face(
                family,p,u,v,w,face,false,f1,f2,1,curl,
                f[offset+1:offset+count])
        end
        tri_all[o]=work[quad_total+1:end]
    end
    return quad_all,tri_all,quad_total
end

# Write one oriented Hcurl/CurlHcurl point into `out`: function-major,
# 3 components per function.
function _hcurl_write_point!(out::Vector{Float64},cursor::Int,
                             et::Vector{_V3},eneg::Vector{_V3},
                             ftab::Vector{_V3},btab::Vector{_V3},
                             quad_all,tri_all,quad_total::Int,
                             family::Symbol,p::Int,perm::Vector{Int})
    blocks=_hcurl_face_blocks(family,p)
    per_edge=p+1
    eout=et
    edges=get(_SOLIN_EDGES,family,())
    if !isempty(edges)
        eout=copy(et)
        for (edge,(a,b)) in enumerate(edges)
            if perm[b+1]<perm[a+1]
                lo=(edge-1)*per_edge
                for k in 1:per_edge
                    eout[lo+k]=eneg[lo+k]
                end
            end
        end
    end
    fout=ftab
    if !isempty(blocks)
        fout=copy(ftab)
        faces=_SOLIN_FACES[family]
        for (face,verts) in enumerate(faces)
            offset,count,kind=blocks[face]
            count==0 && continue
            if kind===:tri
                nums=ntuple(i->perm[verts[i]+1],3)
                f1,f2=_tri_face_flags(nums)
                io=_number_orientation_tri(f1,f2)
                src=tri_all[io+1]
                base=offset-quad_total
                for k in 1:count
                    fout[offset+k]=src[base+k]
                end
            else
                nums=ntuple(i->perm[verts[i]+1],4)
                f1,f2,f3=_quad_face_flags(nums)
                io=_number_orientation_quad(f1,f2,f3)
                src=quad_all[io+1]
                for k in 1:count
                    fout[offset+k]=src[offset+k]
                end
            end
        end
    end
    for f in eout, c in 1:3
        out[cursor+=1]=f[c]
    end
    for f in fout, c in 1:3
        out[cursor+=1]=f[c]
    end
    for f in btab, c in 1:3
        out[cursor+=1]=f[c]
    end
    return cursor
end

function _hcurl_basis_eval(family::Symbol,p::Int,curl::Bool,
                           coords::Vector{Float64},point_count::Int,
                           orientations::Vector{Int},
                           result::Vector{Float64})
    counts=_hcurl_counts(family,p)
    nvertices=family===:lin ? 2 :
              (family===:tri || family===:qua) ? (family===:tri ? 3 : 4) :
              family===:tet ? 4 : family===:hex ? 8 : 6
    per_edge=p+1
    etabs=Vector{Vector{_V3}}(undef,point_count)
    enegs=Vector{Vector{_V3}}(undef,point_count)
    ftabs=Vector{Vector{_V3}}(undef,point_count)
    btabs=Vector{Vector{_V3}}(undef,point_count)
    quads=Vector{Any}(undef,point_count)
    tris=Vector{Any}(undef,point_count)
    quadtotals=Vector{Int}(undef,point_count)
    @inbounds for q in 1:point_count
        u=coords[3q-2]
        v=coords[3q-1]
        w=coords[3q]
        et,ft,bt=curl ? _curl_tables(family,p,u,v,w) :
                        _hcurl_tables(family,p,u,v,w)
        etabs[q]=et
        enegs[q]=isempty(et) ? et : _hcurl_negative_edge_table(et,per_edge)
        ftabs[q]=ft
        btabs[q]=bt
        if isempty(ft)
            quads[q]=nothing
            tris[q]=nothing
            quadtotals[q]=0
        else
            qa,ta,qt=_hcurl_all_oriented_faces(
                family,p,u,v,w,ft,curl)
            quads[q]=qa
            tris[q]=ta
            quadtotals[q]=qt
        end
    end
    cursor=0
    for orientation in orientations
        perm=_orientation_permutation(nvertices,orientation)
        for q in 1:point_count
            cursor=_hcurl_write_point!(
                result,cursor,etabs[q],enegs[q],ftabs[q],btabs[q],
                quads[q],tris[q],quadtotals[q],family,p,perm)
        end
    end
    cursor==length(result) || error(
        "HierarchicalBases: internal Hcurl result length mismatch")
    return result
end
