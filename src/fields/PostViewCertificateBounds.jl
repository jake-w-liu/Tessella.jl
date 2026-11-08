# Constructor-only enclosure of every Float64 query admitted by a scheme's
# intrinsic coordinates, artificial normal coordinates and residual certificate.
@inline _pv_up_add(a,b)=iszero(a) ? b : iszero(b) ? a : nextfloat(a+b)
@inline _pv_up_mul(a,b)=iszero(a) || iszero(b) ? 0.0 : nextfloat(a*b)
function _pv_up_pow(a,degree)
    result=1.0;factor=a;n=Int(degree)
    while n>0
        isodd(n) && (result=_pv_up_mul(result,factor))
        n>>=1;n>0 && (factor=_pv_up_mul(factor,factor))
    end
    return result
end
function _pv_down_pow(a,degree)
    result=1.0;factor=a;n=Int(degree)
    while n>0
        isodd(n) && (result=max(0.0,prevfloat(result*factor)))
        n>>=1;n>0 && (factor=max(0.0,prevfloat(factor*factor)))
    end
    return result
end

function _pv_certificate_canonical_basis(kind)
    if kind==1
        return [.5 -.5;.5 .5],[(0.,0.,0.),(1.,0.,0.)]
    elseif kind==2
        return [1. -1. -1.;0. 1. 0.;0. 0. 1.],[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.)]
    elseif kind==3
        return [1. -1. -1. -1.;0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1.],
            [(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)]
    elseif kind==7
        signs=((-1.,-1.),(1.,-1.),(1.,1.),(-1.,1.))
        c=[.25*(j==1 ? 1. : j==2 ? a : j==3 ? b : a*b) for (a,b) in signs,j in 1:4]
        return c,[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(1.,1.,0.)]
    elseif kind==4
        signs=((-1.,-1.,-1.),(1.,-1.,-1.),(1.,1.,-1.),(-1.,1.,-1.),
               (-1.,-1.,1.),(1.,-1.,1.),(1.,1.,1.),(-1.,1.,1.))
        exponents=[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),
                   (1.,1.,0.),(1.,0.,1.),(0.,1.,1.),(1.,1.,1.)]
        c=[.125*a^e[1]*b^e[2]*d^e[3] for (a,b,d) in signs,e in exponents]
        return c,exponents
    end
    return [.5 -.5 -.5 -.5 .5 .5;0. .5 0. 0. -.5 0.;0. 0. .5 0. 0. -.5;
            .5 -.5 -.5 .5 -.5 -.5;0. .5 0. 0. .5 0.;0. 0. .5 0. 0. .5],
        [(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(1.,0.,1.),(0.,1.,1.)]
end

# On a tolerance-extended collapsed pyramid, |u|,|v| <= 1-w+tol.
# For delta=1-w, (delta+tol)^s/max(epsilon,delta)^d has no interior
# maximum on [epsilon,1+tol]: its stationary point is a minimum. This
# preserves the collapse relation instead of using a whole tensor box.
function _pv_pyramid_factor(s,d,tol,epsilon=1e-14)
    first=_pv_up_pow(_pv_up_add(epsilon,tol),s)/_pv_down_pow(epsilon,d)
    delta=_pv_up_add(1.0,tol)
    last=_pv_up_pow(_pv_up_add(delta,tol),s)/_pv_down_pow(delta,d)
    return nextfloat(max(first,last))
end

function _pv_certificate_monomial(kind,e,tol)
    all(x->isfinite(x) && x>=0 && isinteger(x) && x<Float64(typemax(Int)),e) || return nothing
    _pv_geometry_zero_term(kind,e) && return ((0,0),0.0,(0.0,0.0,0.0))
    zero_range=_pv_scheme_monomial_range(kind,e)
    if kind==6
        m=min(e[1],e[2]);s=e[1]+e[2];radius=_pv_up_add(1.0,tol)
        magnitude=_pv_up_mul(_pv_pyramid_factor(s,m,tol),_pv_up_pow(radius,e[3]))
        signed=rem(e[1],2.)==1. || rem(e[2],2.)==1. ||
            (tol>0 && rem(e[3],2.)==1.)
        return (signed ? (-magnitude,magnitude) : (0.0,magnitude)),magnitude,(0.,0.,0.)
    end
    dimension=_pv_scheme_dim(kind)
    factors=ntuple(3) do axis
        axis>dimension && return 0.0
        simplex=kind in (2,3) || (kind==5 && axis<=2)
        return _pv_up_add(1.0,(simplex ? (kind==3 ? 3 : 2) : 1)*tol)
    end
    magnitude=1.0
    for axis in 1:3
        magnitude=_pv_up_mul(magnitude,_pv_up_pow(factors[axis],e[axis]))
    end
    derivatives=ntuple(3) do axis
        (axis>dimension || iszero(e[axis])) && return 0.0
        value=e[axis]
        for d in 1:3;value=_pv_up_mul(value,_pv_up_pow(factors[d],e[d]-(d==axis)));end
        return value
    end
    return zero_range,magnitude,derivatives
end

function _pv_scheme_certificate_bounds(cell,tol;normal_frame=nothing)
    extent=floatmax(Float64);full=((-extent,-extent,-extent),(extent,extent,extent))
    base=_pv_scheme_base(cell);scale=_pv_scheme_scale(cell)
    geometry=cell.geometry;kind=cell.kind;dim=_pv_scheme_dim(kind)
    intrinsic=_pv_up_add(tol,_pv_up_mul(8eps(Float64),max(1.0,1+tol)))
    lo=zeros(3);hi=zeros(3);mass=zeros(3);derivatives=zeros(3,3)
    operations=32.0
    if cell.flat && kind==6
        # Canonical flat pyramids use their separate unclamped w!=1 branch.
        # The nearest distinct Float64 below1 is1-2^-53, so this floating
        # program has a finite denominator bound even in the allowed strip.
        radius=_pv_up_add(1.0,2intrinsic)
        rational=_pv_up_mul(_pv_up_add(1.0,intrinsic),_pv_pyramid_factor(2.,1.,intrinsic,2.0^-53))
        shape=_pv_up_mul(.25,_pv_up_add(_pv_up_pow(_pv_up_add(1.0,radius),2),
            _pv_up_add(_pv_up_add(1.0,intrinsic),rational)))
        for axis in 1:3,node in axes(cell.coords,2)
            weight=node==5 ? _pv_up_add(1.0,intrinsic) : shape
            value=_pv_up_mul(abs(_pv_relative(cell,node,base,axis,scale)),weight)
            mass[axis]=_pv_up_add(mass[axis],value)
        end
        lo.=-mass;hi.=mass
        operations=256.0
    else
        coefficients,exponents=geometry!==nothing ?
            (geometry.coefficients,geometry.exponents) : cell.flat ?
            _pv_certificate_canonical_basis(kind) : (cell.coefgeo,cell.expgeo)
        nodes=geometry===nothing ? size(coefficients,1) : 1
        operations=nextfloat(32.0+8.0*nodes*length(exponents))
        for term in eachindex(exponents)
            bounded=_pv_certificate_monomial(kind,exponents[term],intrinsic)
            bounded===nothing && return full
            range,magnitude,gradient=bounded
            for degree in exponents[term];operations=nextfloat(operations+16degree);end
            all(isfinite,(magnitude,gradient...)) || return full
            for axis in 1:3
                alpha_lo=0.0;alpha_hi=0.0
                for node in 1:nodes
                    q=geometry===nothing ? _pv_relative(cell,node,base,axis,scale) : 1.0
                    c=geometry===nothing ? coefficients[node,term] : coefficients[axis,term]
                    (iszero(q) || iszero(c)) && continue
                    product=q*c
                    alpha_lo=prevfloat(alpha_lo+prevfloat(product))
                    alpha_hi=nextfloat(alpha_hi+nextfloat(product))
                    weight=_pv_up_mul(abs(q),abs(c))
                    mass[axis]=_pv_up_add(mass[axis],_pv_up_mul(weight,magnitude))
                    for d in 1:dim
                        derivatives[axis,d]=_pv_up_add(derivatives[axis,d],_pv_up_mul(weight,gradient[d]))
                    end
                end
                (iszero(alpha_lo) && iszero(alpha_hi)) && continue
                if kind==6
                    products=(alpha_lo*range[1],alpha_lo*range[2],alpha_hi*range[1],alpha_hi*range[2])
                    lower=prevfloat(minimum(products));upper=nextfloat(maximum(products))
                elseif range==(0,0)
                    lower=0.;upper=0.
                elseif range==(-1,1)
                    radius=max(abs(alpha_lo),abs(alpha_hi));lower=-radius;upper=radius
                elseif range==(1,1)
                    lower=alpha_lo;upper=alpha_hi
                else
                    lower=min(0.,alpha_lo);upper=max(0.,alpha_hi)
                end
                lo[axis]=prevfloat(lo[axis]+lower);hi[axis]=nextfloat(hi[axis]+upper)
            end
        end
    end
    mu=nextfloat(operations*(eps(Float64)/2))
    isfinite(mu) && mu<.5 || return full
    gamma=nextfloat(mu/prevfloat(1-mu))
    for axis in 1:3
        floor=_pv_up_mul(_pv_up_mul(nextfloat(0.0),operations),_pv_up_add(mass[axis],1.0))
        error=_pv_up_add(_pv_up_mul(gamma,mass[axis]),floor)
        lo[axis]=prevfloat(lo[axis]-error);hi[axis]=nextfloat(hi[axis]+error)
        for d in 1:dim
            value=derivatives[axis,d]
            derivatives[axis,d]=iszero(value) ? 0.0 : _pv_up_add(value,_pv_up_add(_pv_up_mul(gamma,value),floor))
        end
    end
    # Extend exact-domain map bounds by a mean-value bound. Pyramids above
    # already use their full tolerance-extended collapsed domain directly.
    if kind!=6
        delta=kind in (2,3) ? dim*intrinsic : intrinsic
        for axis in 1:3,d in 1:dim
            width=_pv_up_mul(derivatives[axis,d],kind==5 && d<=2 ? 2intrinsic : delta)
            lo[axis]=prevfloat(lo[axis]-width);hi[axis]=nextfloat(hi[axis]+width)
        end
    end
    ratio=geometry===nothing ? 1.0 : abs(geometry.normal_ratio)
    normal=zeros(3)
    if normal_frame!==nothing
        for axis in 1:3
            normal[axis]=_pv_up_mul(tol,_pv_up_add(normal_frame[1][axis],normal_frame[2][axis]))
        end
    elseif dim==1
        d=Tuple(derivatives[:,1]);a=(d[2],max(d[1],d[3]),d[2])
        for axis in 1:3
            i=mod1(axis+1,3);j=mod1(axis+2,3)
            cross=_pv_up_mul(ratio,_pv_up_add(_pv_up_mul(d[i],a[j]),_pv_up_mul(d[j],a[i])))
            normal[axis]=_pv_up_mul(tol,_pv_up_add(a[axis],cross))
        end
    elseif dim==2
        for axis in 1:3
            i=mod1(axis+1,3);j=mod1(axis+2,3)
            cross=_pv_up_add(_pv_up_mul(derivatives[i,1],derivatives[j,2]),
                             _pv_up_mul(derivatives[j,1],derivatives[i,2]))
            normal[axis]=_pv_up_mul(tol,_pv_up_mul(ratio,cross))
        end
    end
    for axis in 1:3
        lo[axis]=prevfloat(lo[axis]-normal[axis]);hi[axis]=nextfloat(hi[axis]+normal[axis])
    end
    all(isfinite,lo) && all(isfinite,hi) || return full
    # The inverse's final intrinsic residual is at most max(1e-10,256eps*Q).
    # Solve Q<=B+error(Q), also enclosing the few residual/query roundings.
    k=nextfloat(max(2048eps(Float64),gamma));k<.5 || return full
    bound=max(maximum(abs,lo),maximum(abs,hi))
    radius=nextfloat(_pv_up_add(bound,1e-10)/prevfloat(1-k))
    certificate=_pv_up_add(1e-10,_pv_up_mul(k,max(1.0,radius)))
    lo.=prevfloat.(lo.-certificate);hi.=nextfloat.(hi.+certificate)
    bounds=ntuple(3) do axis
        # Float queries between these adjacent physical coordinates do not
        # exist. This keeps huge constant translations from inflating every
        # axis of an otherwise tiny box through diagonal padding.
        lower_gap=abs((prevfloat(base[axis])-base[axis])/scale)
        upper_gap=abs((nextfloat(base[axis])-base[axis])/scale)
        if lo[axis]>-lower_gap && hi[axis]<upper_gap
            return (base[axis],base[axis])
        end
        error=_pv_up_mul(gamma,_pv_up_add(abs(base[axis]),_pv_up_mul(scale,mass[axis])))
        lower=prevfloat(base[axis]+prevfloat(scale*lo[axis]))
        upper=nextfloat(base[axis]+nextfloat(scale*hi[axis]))
        return (prevfloat(lower-error),nextfloat(upper+error))
    end
    all(b->all(isfinite,b),bounds) || return full
    return (ntuple(axis->bounds[axis][1],3),ntuple(axis->bounds[axis][2],3))
end

# Compact cells use the same intrinsic maps, with distinct artificial
# columns for Lines and Quadrangles. Bind those columns from their own
# evaluators; a canonical scheme's derivative-based frame is not equivalent.
function _pv_compact_certificate_bounds(coords,cells,column,kind,tol)
    C=coords[:,cells[:,column]]
    cell=_PVSchemeCell(kind,true,kind==6,zeros(0,0),NTuple{3,Float64}[],
        zeros(0,0),NTuple{3,Float64}[],C,zeros(0,0))
    scale=_pv_scheme_scale(cell);base=_pv_scheme_base(cell)
    frame=nothing
    if kind==1
        d=ntuple(axis->_pv_relative(cell,2,base,axis,scale),3)
        artificial=if (abs(d[1])>=abs(d[2]) && abs(d[1])>=abs(d[3])) ||
                      (abs(d[2])>=abs(d[1]) && abs(d[2])>=abs(d[3]))
            (abs(d[2]),abs(d[1]),0.0)
        else
            (0.0,abs(d[3]),abs(d[2]))
        end
        normal=ntuple(3) do axis
            i=mod1(axis+1,3);j=mod1(axis+2,3)
            _pv_up_mul(scale,_pv_up_add(_pv_up_mul(abs(d[i]),artificial[j]),
                                      _pv_up_mul(abs(d[j]),artificial[i])))
        end
        frame=(artificial,normal)
    elseif kind==2
        frame=((0.,0.,0.),(0.,0.,0.))
    elseif kind==7
        local_cells=reshape(Int32.(1:4),4,1)
        _,_,normal=_postview_reference_jacobian(C,local_cells,1,4,kind,0.,0.,0.,scale)
        frame=((0.,0.,0.),map(abs,normal))
    end
    lower,upper=_pv_scheme_certificate_bounds(cell,tol;normal_frame=frame)
    return lower,upper
end
