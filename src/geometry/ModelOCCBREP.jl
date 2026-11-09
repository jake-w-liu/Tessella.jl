# Native OCCT ASCII BREP interchange. Grammar reference:
# https://github.com/Open-Cascade-SAS/OCCT/wiki/brep_format
# Parse into an independent representation before allocating kernel entities.
mutable struct _OCCBrepCursor
    tokens::Vector{SubString{String}}
    position::Int
    caller::String
    version::Int
end

function _occ_brep_token!(c::_OCCBrepCursor)
    c.position<=length(c.tokens) || throw(ArgumentError(
        "$(c.caller): truncated BREP record"))
    t=c.tokens[c.position]
    c.position+=1
    return t
end

function _occ_brep_number!(c::_OCCBrepCursor,::Type{T}=Float64) where T
    token=_occ_brep_token!(c)
    value=tryparse(T,token)
    value!==nothing && isfinite(value) || throw(ArgumentError(
        "$(c.caller): invalid BREP $(T) value $(repr(token))"))
    return value::T
end

function _occ_brep_count!(c::_OCCBrepCursor;minimum=0)
    n=_occ_brep_number!(c,Int)
    minimum<=n<=length(c.tokens)-c.position+1 || throw(ArgumentError(
        "$(c.caller): BREP count $n exceeds the available records"))
    return n
end

function _occ_brep_flag!(c::_OCCBrepCursor)
    n=_occ_brep_number!(c,Int)
    n in (0,1) || throw(ArgumentError("$(c.caller): invalid BREP flag $n"))
    return n==1
end

_occ_brep_point!(c,n=3)=ntuple(_->_occ_brep_number!(c),n)

function _occ_brep_header!(c,name)
    _occ_brep_token!(c)==name || throw(ArgumentError(
        "$(c.caller): expected BREP $name section"))
    return _occ_brep_count!(c)
end

function _occ_brep_knots!(c,count;limit=length(c.tokens))
    knots=Float64[]
    last=-Inf
    for _ in 1:count
        value=_occ_brep_number!(c)
        multiplicity=_occ_brep_count!(c;minimum=1)
        value>last || throw(ArgumentError(
            "$(c.caller): BREP distinct knots must increase"))
        multiplicity<=limit-length(knots) || throw(ArgumentError(
            "$(c.caller): expanded BREP knot count exceeds its pole/degree bound"))
        append!(knots,Iterators.repeated(value,multiplicity))
        last=value
    end
    return knots
end

function _occ_brep_curve!(c,n,depth=0)
    depth<64 || throw(ArgumentError("$(c.caller): BREP curve nesting exceeds 64"))
    kind=_occ_brep_number!(c,Int)
    if kind==1
        return (kind=kind,origin=_occ_brep_point!(c,n),dir=_occ_brep_point!(c,n))
    elseif kind in 2:5
        center=_occ_brep_point!(c,n)
        normal=n==3 ? _occ_brep_point!(c,3) : nothing
        X=_occ_brep_point!(c,n); Y=_occ_brep_point!(c,n)
        a=_occ_brep_number!(c)
        b=kind in (3,5) ? _occ_brep_number!(c) : a
        return (kind=kind,center=center,normal=normal,X=X,Y=Y,a=a,b=b)
    elseif kind in (6,7)
        rational=_occ_brep_flag!(c)
        periodic=kind==7 ? _occ_brep_flag!(c) : false
        degree=_occ_brep_number!(c,Int)
        1<=degree<=25 || throw(ArgumentError(
            "$(c.caller): BREP curve degree must be in 1:25"))
        count=kind==6 ? degree+1 : _occ_brep_count!(c;minimum=2)
        nk=kind==7 ? _occ_brep_count!(c;minimum=2) : 2
        poles=Tuple[]; weights=Float64[]
        for _ in 1:count
            push!(poles,_occ_brep_point!(c,n))
            push!(weights,rational ? _occ_brep_number!(c) : 1.0)
        end
        knots=kind==6 ? vcat(zeros(degree+1),ones(degree+1)) :
            _occ_brep_knots!(c,nk;limit=count+degree+1)
        return (kind=kind,degree=degree,periodic=periodic,poles=poles,
                weights=weights,knots=knots)
    elseif kind==8
        lo=_occ_brep_number!(c); hi=_occ_brep_number!(c)
        lo<hi || throw(ArgumentError("$(c.caller): invalid BREP curve trim"))
        return (kind=kind,lo=lo,hi=hi,base=_occ_brep_curve!(c,n,depth+1))
    elseif kind==9
        distance=_occ_brep_number!(c)
        direction=n==3 ? _occ_brep_point!(c,3) : nothing
        return (kind=kind,distance=distance,direction=direction,
                base=_occ_brep_curve!(c,n,depth+1))
    end
    throw(ArgumentError("$(c.caller): unknown BREP $n-D curve kind $kind"))
end

function _occ_brep_surface!(c,depth=0)
    depth<64 || throw(ArgumentError("$(c.caller): BREP surface nesting exceeds 64"))
    kind=_occ_brep_number!(c,Int)
    if kind in 1:5
        center=_occ_brep_point!(c); axis=_occ_brep_point!(c)
        X=_occ_brep_point!(c); Y=_occ_brep_point!(c)
        a=kind==1 ? 0.0 : _occ_brep_number!(c)
        b=kind in (3,5) ? _occ_brep_number!(c) : 0.0
        return (kind=kind,center=center,axis=axis,X=X,Y=Y,a=a,b=b)
    elseif kind==6
        direction=_occ_brep_point!(c)
        return (kind=kind,direction=direction,curve=_occ_brep_curve!(c,3))
    elseif kind==7
        center=_occ_brep_point!(c); axis=_occ_brep_point!(c)
        return (kind=kind,center=center,axis=axis,curve=_occ_brep_curve!(c,3))
    elseif kind in (8,9)
        rational_u=_occ_brep_flag!(c); rational_v=_occ_brep_flag!(c)
        periodic_u=kind==9 ? _occ_brep_flag!(c) : false
        periodic_v=kind==9 ? _occ_brep_flag!(c) : false
        pu=_occ_brep_number!(c,Int); pv=_occ_brep_number!(c,Int)
        1<=pu<=25 && 1<=pv<=25 || throw(ArgumentError(
            "$(c.caller): BREP surface degrees must be in 1:25"))
        nu=kind==8 ? pu+1 : _occ_brep_count!(c;minimum=2)
        nv=kind==8 ? pv+1 : _occ_brep_count!(c;minimum=2)
        nku=kind==9 ? _occ_brep_count!(c;minimum=2) : 2
        nkv=kind==9 ? _occ_brep_count!(c;minimum=2) : 2
        nu<=div(length(c.tokens)-c.position+1,3nv) || throw(ArgumentError(
            "$(c.caller): BREP surface pole count exceeds available records"))
        poles=Matrix{NTuple{3,Float64}}(undef,nu,nv)
        weights=ones(nu,nv)
        for i in 1:nu,j in 1:nv
            poles[i,j]=_occ_brep_point!(c)
            (rational_u || rational_v) && (weights[i,j]=_occ_brep_number!(c))
        end
        U=kind==8 ? vcat(zeros(pu+1),ones(pu+1)) : _occ_brep_knots!(c,nku;limit=nu+pu+1)
        V=kind==8 ? vcat(zeros(pv+1),ones(pv+1)) : _occ_brep_knots!(c,nkv;limit=nv+pv+1)
        return (kind=kind,degree_u=pu,degree_v=pv,periodic_u=periodic_u,
            periodic_v=periodic_v,poles=poles,weights=weights,U=U,V=V)
    elseif kind==10
        bounds=_occ_brep_point!(c,4)
        bounds[1]<bounds[2] && bounds[3]<bounds[4] || throw(ArgumentError(
            "$(c.caller): invalid BREP surface trim"))
        return (kind=kind,bounds=bounds,base=_occ_brep_surface!(c,depth+1))
    elseif kind==11
        distance=_occ_brep_number!(c)
        return (kind=kind,distance=distance,base=_occ_brep_surface!(c,depth+1))
    end
    throw(ArgumentError("$(c.caller): unknown BREP surface kind $kind"))
end

function _occ_brep_reference!(c,nshapes,nlocations;current=0)
    token=_occ_brep_token!(c)
    occursin(r"^[+\-ie]\d+$",token) || throw(ArgumentError(
        "$(c.caller): invalid BREP topology reference $(repr(token))"))
    shape=tryparse(Int,SubString(token,2))
    shape!==nothing && current<shape<=nshapes || throw(ArgumentError(
        "$(c.caller): BREP topology reference must name an earlier record"))
    location=_occ_brep_number!(c,Int)
    0<=location<=nlocations || throw(ArgumentError(
        "$(c.caller): invalid BREP location reference $location"))
    return (shape=shape,orientation=first(token),location=location)
end

function _occ_brep_location_product(first,second,caller)
    isempty(second) && return first
    terms=collect(first)
    isempty(terms) && (terms=Tuple{Int,Int}[])
    for (datum,power) in second
        if !isempty(terms) && terms[end][1]==datum
            combined=try
                Base.checked_add(terms[end][2],power)
            catch
                throw(ArgumentError("$caller: BREP location exponent exceeds Int"))
            end
            pop!(terms)
            combined==0 || push!(terms,(datum,combined))
        elseif power!=0
            push!(terms,(datum,power))
        end
        length(terms)<=4096 || throw(ArgumentError(
            "$caller: expanded BREP location exceeds 4096 elementary terms"))
    end
    return Tuple(terms)
end

function _occ_brep_location_power(word,power,caller)
    power==1 && return word
    if power<0
        word=Tuple((datum,Base.checked_sub(0,exponent)) for
                   (datum,exponent) in reverse(word))
        power=Base.checked_sub(0,power)
    end
    result=()
    while power>0
        isodd(power) && (result=_occ_brep_location_product(result,word,caller))
        power>>=1
        power>0 && (word=_occ_brep_location_product(word,word,caller))
    end
    return result
end

function _occ_brep_location_matrix(word,locations,caller)
    result=Matrix{Float64}(LinearAlgebra.I,4,4)
    for (datum,power) in word
        result=result*(locations[datum]^power)
        all(isfinite,result) || throw(ArgumentError(
            "$caller: BREP location composition exceeds Float64"))
    end
    return result
end

function _occ_brep_parse(source::String,caller)
    version=match(r"CASCADE Topology V([123]),",source)
    version===nothing && throw(ArgumentError("$caller: invalid ASCII BREP header"))
    start=findnext("Locations",source,version.offset)
    start===nothing && throw(ArgumentError("$caller: missing BREP Locations section"))
    # `2CN` is legal in closed-surface edge representations. Splitting the
    # continuity suffix is the sole non-whitespace token boundary in the file.
    payload=replace(source[first(start):end],r"(?<=\d)([CG][0123N])"=>s" \1")
    c=_OCCBrepCursor(split(payload),1,String(caller),parse(Int,version.captures[1]))
    locations=Matrix{Float64}[]
    location_words=Tuple[]
    word_capacity=max(4096,16length(c.tokens))
    for location_index in 1:_occ_brep_header!(c,"Locations")
        kind=_occ_brep_number!(c,Int)
        matrix=Matrix{Float64}(LinearAlgebra.I,4,4)
        word=()
        if kind==1
            for i in 1:3,j in 1:4
                matrix[i,j]=_occ_brep_number!(c)
            end
            # Float products can underflow for a perfectly invertible
            # location. Nine exact Float64 rationals make this fixed-size
            # rank check independent of coordinate magnitude.
            a=Rational{BigInt}.(@view matrix[1:3,1:3])
            determinant=a[1,1]*(a[2,2]*a[3,3]-a[2,3]*a[3,2])-
                a[1,2]*(a[2,1]*a[3,3]-a[2,3]*a[3,1])+
                a[1,3]*(a[2,1]*a[3,2]-a[2,2]*a[3,1])
            !iszero(determinant) || throw(ArgumentError(
                "$caller: singular BREP location"))
            word=((location_index,1),)
        elseif kind==2
            while true
                ref=_occ_brep_number!(c,Int)
                ref==0 && break
                1<=ref<=length(locations) || throw(ArgumentError(
                    "$caller: invalid composed BREP location"))
                power=_occ_brep_number!(c,Int)
                typemin(Int32)<=power<=typemax(Int32) || throw(ArgumentError(
                    "$caller: BREP location power exceeds Int32"))
                word=_occ_brep_location_product(word,
                    _occ_brep_location_power(location_words[ref],power,caller),caller)
            end
            matrix=_occ_brep_location_matrix(word,locations,caller)
        else
            throw(ArgumentError("$caller: unknown BREP location kind $kind"))
        end
        length(word)<=word_capacity || throw(ArgumentError(
            "$caller: expanded BREP location tables exceed the source record budget"))
        word_capacity-=length(word)
        push!(locations,matrix)
        push!(location_words,word)
    end
    curves2=[_occ_brep_curve!(c,2) for _ in 1:_occ_brep_header!(c,"Curve2ds")]
    curves3=[_occ_brep_curve!(c,3) for _ in 1:_occ_brep_header!(c,"Curves")]
    # Stored tessellations are optional display caches. Consume and validate
    # their declared records; the CAD importer uses the analytic geometry.
    for _ in 1:_occ_brep_header!(c,"Polygon3D")
        n=_occ_brep_count!(c;minimum=2); parameters=_occ_brep_flag!(c)
        _occ_brep_number!(c)
        for _ in 1:n; _occ_brep_point!(c); end
        parameters && foreach(_->_occ_brep_number!(c),1:n)
    end
    for _ in 1:_occ_brep_header!(c,"PolygonOnTriangulations")
        n=_occ_brep_count!(c;minimum=2)
        for _ in 1:n; _occ_brep_number!(c,Int); end
        _occ_brep_token!(c)=="p" || throw(ArgumentError(
            "$caller: invalid BREP polygon parameters"))
        _occ_brep_number!(c); parameters=_occ_brep_flag!(c)
        parameters && foreach(_->_occ_brep_number!(c),1:n)
    end
    surfaces=[_occ_brep_surface!(c) for _ in 1:_occ_brep_header!(c,"Surfaces")]
    for _ in 1:_occ_brep_header!(c,"Triangulations")
        n=_occ_brep_count!(c;minimum=3); nt=_occ_brep_count!(c;minimum=1)
        uv=_occ_brep_flag!(c); normals=c.version==3 ? _occ_brep_flag!(c) : false
        _occ_brep_number!(c)
        for _ in 1:n; _occ_brep_point!(c); end
        uv && c.version>=2 && foreach(_->_occ_brep_point!(c,2),1:n)
        for _ in 1:nt,_ in 1:3
            node=_occ_brep_number!(c,Int)
            1<=node<=n || throw(ArgumentError("$caller: invalid BREP triangle node"))
        end
        normals && foreach(_->_occ_brep_point!(c),1:n)
    end
    count=_occ_brep_header!(c,"TShapes")
    shapes=Vector{NamedTuple}(undef,count)
    for index in count:-1:1
        kind=_occ_brep_token!(c)
        data=if kind=="Ve"
            tolerance=_occ_brep_number!(c); point=_occ_brep_point!(c)
            while true
                _occ_brep_number!(c)
                representation=_occ_brep_number!(c,Int)
                representation==0 && break
                representation in 1:3 || throw(ArgumentError(
                    "$caller: invalid BREP vertex representation"))
                representation==3 ? _occ_brep_number!(c) : _occ_brep_number!(c,Int)
                representation in (2,3) && _occ_brep_number!(c,Int)
                _occ_brep_number!(c,Int)
            end
            (tolerance=tolerance,point=point)
        elseif kind=="Ed"
            tolerance=_occ_brep_number!(c)
            same_parameter=_occ_brep_flag!(c); same_range=_occ_brep_flag!(c)
            degenerate=_occ_brep_flag!(c); representations=NamedTuple[]
            while true
                representation=_occ_brep_number!(c,Int)
                representation==0 && break
                if representation==1
                    curve=_occ_brep_number!(c,Int); location=_occ_brep_number!(c,Int)
                    lo=_occ_brep_number!(c); hi=_occ_brep_number!(c)
                    push!(representations,(kind=representation,curve=curve,
                        location=location,lo=lo,hi=hi))
                elseif representation in (2,3)
                    curve=_occ_brep_number!(c,Int)
                    curve2=representation==3 ? _occ_brep_number!(c,Int) : 0
                    representation==3 && _occ_brep_token!(c)
                    surface=_occ_brep_number!(c,Int); location=_occ_brep_number!(c,Int)
                    lo=_occ_brep_number!(c); hi=_occ_brep_number!(c)
                    c.version>=2 && _occ_brep_point!(c,4)
                    push!(representations,(kind=representation,curve=curve,
                        curve2=curve2,surface=surface,location=location,lo=lo,hi=hi))
                elseif representation==4
                    _occ_brep_token!(c)
                    for _ in 1:4; _occ_brep_number!(c,Int); end
                elseif representation in (5,6,7)
                    for _ in 1:(representation-3); _occ_brep_number!(c,Int); end
                else
                    throw(ArgumentError("$caller: invalid BREP edge representation"))
                end
            end
            (tolerance=tolerance,same_parameter=same_parameter,same_range=same_range,
             degenerate=degenerate,representations=representations)
        elseif kind=="Fa"
            natural=_occ_brep_flag!(c); tolerance=_occ_brep_number!(c)
            surface=_occ_brep_number!(c,Int); location=_occ_brep_number!(c,Int)
            if c.position<=length(c.tokens) && c.tokens[c.position]=="2"
                _occ_brep_token!(c); _occ_brep_number!(c,Int)
            end
            (natural=natural,tolerance=tolerance,surface=surface,location=location)
        elseif kind in ("Wi","Sh","So","CS","Co")
            NamedTuple()
        else
            throw(ArgumentError("$caller: unknown BREP topology kind $(repr(kind))"))
        end
        flags=_occ_brep_token!(c)
        occursin(r"^[01]{7}$",flags) || throw(ArgumentError(
            "$caller: invalid BREP shape flags"))
        children=NamedTuple[]
        while c.position<=length(c.tokens) && c.tokens[c.position]!="*"
            push!(children,_occ_brep_reference!(c,count,length(locations);current=index))
        end
        _occ_brep_token!(c)=="*" || throw(ArgumentError("$caller: missing BREP shape terminator"))
        shapes[index]=(kind=String(kind),data=data,flags=String(flags),children=children)
    end
    root=_occ_brep_reference!(c,count,length(locations))
    c.position==length(c.tokens)+1 || throw(ArgumentError(
        "$caller: unexpected data following the BREP root"))
    return (locations=locations,location_words=location_words,curves2=curves2,curves3=curves3,
            surfaces=surfaces,shapes=shapes,root=root)
end

function _occ_brep_lookup(items,index,caller,what)
    1<=index<=length(items) || throw(ArgumentError(
        "$caller: invalid BREP $what reference $index"))
    return items[index]
end

function _occ_brep_untrim(curve)
    while curve.kind==8; curve=curve.base; end
    if curve.kind==9
        base=_occ_brep_untrim(curve.base)
        d=curve.distance
        if base.kind==1
            direction=base.dir
            normal=length(direction)==2 ? (direction[2],-direction[1]) :
                _occ_cross(direction,curve.direction)
            magnitude=sqrt(sum(abs2,normal))
            magnitude>0 || throw(ArgumentError("BREP offset line has a singular normal"))
            return merge(base,(origin=base.origin.+(d/magnitude).*normal,))
        elseif base.kind==2
            handedness=length(base.center)==2 ?
                base.X[1]*base.Y[2]-base.X[2]*base.Y[1] :
                _dot3(base.normal,curve.direction)
            length(base.center)==3 && abs(abs(handedness)-1)>1e-10 && return curve
            radius=base.a+d*sign(handedness)
            radius>0 || throw(ArgumentError("BREP circle offset collapses its radius"))
            return merge(base,(a=radius,b=radius))
        end
    end
    return curve
end

function _occ_brep_nurbs_curve(curve,caller)
    curve=_occ_brep_untrim(curve)
    curve.kind in (6,7) || throw(ArgumentError(
        "$caller: BREP curve kind $(curve.kind) is not a spline"))
    poles=NTuple{3,Float64}[length(p)==2 ? (p[1],p[2],0.0) : p for p in curve.poles]
    U=copy(curve.knots); weights=copy(curve.weights)
    if curve.periodic
        # A periodic OCCT knot sequence lists one period, including its
        # repeated seam knot. Extend both the cyclic poles and knot pattern.
        period=U[end]-U[1]
        period>0 || throw(ArgumentError("$caller: invalid periodic BREP knots"))
        seam=count(==(U[1]),U)
        seam==count(==(U[end]),U) || throw(ArgumentError(
            "$caller: periodic BREP seam multiplicities disagree"))
        n=length(poles); p=curve.degree
        base=U[1:end-seam]
        length(base)==n || throw(ArgumentError(
            "$caller: periodic BREP knot and pole counts disagree"))
        shift=p+1-seam
        poles=vcat(poles,poles[1:p]); weights=vcat(weights,weights[1:p])
        # The cyclic extension above has n+p controls, hence n+2p+1 knots.
        U=Float64[base[mod(i,n)+1]+fld(i,n)*period for i in -shift:n+2p-shift]
    end
    return NURBS.NURBSCurve(curve.degree,U,poles,weights)
end

function _occ_brep_pcurve(curve,caller)
    curve=_occ_brep_untrim(curve)
    curve.kind==1 && return (lin2d=(curve.origin,curve.dir),)
    if curve.kind in (2,3)
        X,Y=curve.X,curve.Y
        base=curve.kind==2 ? (circ2d=(curve.center,X,curve.a),) :
            (ellips2d=(curve.center,X,curve.a,curve.b),)
        return X[1]*Y[2]-X[2]*Y[1]<0 ? (neg2d=(base,)) : base
    elseif curve.kind in (6,7)
        nc=_occ_brep_nurbs_curve(curve,caller)
        if all(==(nc.weights[1]),nc.weights)
            return (spl2d=(nc.knots,nc.degree,[(p[1],p[2]) for p in nc.controls]),)
        end
        return (brep_nurbs2d=_ModelOCCNURBSCurve(:bspline,nc),)
    elseif curve.kind in (4,5)
        return (brep_conic2d=curve,)
    end
    throw(ArgumentError("$caller: native BREP pcurve kind $(curve.kind) is not implemented"))
end

function _occ_pcurve_d2(pc::NamedTuple{(:brep_conic2d,)},t::Float64)
    curve=pc.brep_conic2d
    if curve.kind==4
        curve.a>0 || throw(ArgumentError("BREP parabola focal length must be positive"))
        x=t*t/(4curve.a); y=t
        dx=t/(2curve.a); dy=1.0
        ddx=1/(2curve.a); ddy=0.0
    else
        s=Base.sinh(t); c=Base.cosh(t)
        x=curve.a*c; y=curve.b*s
        dx=curve.a*s; dy=curve.b*c
        ddx=x; ddy=y
    end
    p=ntuple(i->curve.center[i]+x*curve.X[i]+y*curve.Y[i],2)
    d=ntuple(i->dx*curve.X[i]+dy*curve.Y[i],2)
    dd=ntuple(i->ddx*curve.X[i]+ddy*curve.Y[i],2)
    return p,d,dd
end
_occ_pcurve_d0(pc::NamedTuple{(:brep_conic2d,)},t::Float64)=first(_occ_pcurve_d2(pc,t))
function _occ_pcurve_d1(pc::NamedTuple{(:brep_conic2d,)},t::Float64)
    p,d,_=_occ_pcurve_d2(pc,t)
    return p,d
end

function _occ_pcurve_d0(pc::NamedTuple{(:brep_nurbs2d,)},t::Float64)
    p,_,_=_model_occ_curve_d2(pc.brep_nurbs2d,t)
    return p[1],p[2]
end
function _occ_pcurve_d1(pc::NamedTuple{(:brep_nurbs2d,)},t::Float64)
    p,d,_=_model_occ_curve_d2(pc.brep_nurbs2d,t)
    return (p[1],p[2]),(d[1],d[2])
end
function _occ_pcurve_d2(pc::NamedTuple{(:brep_nurbs2d,)},t::Float64)
    p,d,dd=_model_occ_curve_d2(pc.brep_nurbs2d,t)
    return (p[1],p[2]),(d[1],d[2]),(dd[1],dd[2])
end

# General affine maps in edge parameter and surface UV coordinates. This
# keeps rational traces exact and does not sample them into polylines.
function _occ_brep_pcurve_map(pc,alpha,beta,scale,shift)
    alpha==1.0 && beta==0.0 && scale==(1.0,1.0) && shift==(0.0,0.0) && return pc
    if hasproperty(pc,:lin2d)
        o,d=pc.lin2d
        return (lin2d=(ntuple(i->scale[i]*(o[i]+beta*d[i])+shift[i],2),
                       ntuple(i->scale[i]*alpha*d[i],2)),)
    end
    return (brep_map2d=(pc,alpha,beta,scale,shift),)
end

function _occ_pcurve_d0(pc::NamedTuple{(:brep_map2d,)},t::Float64)
    base,a,b,scale,shift=pc.brep_map2d
    p=_occ_pcurve_d0(base,a*t+b)
    return ntuple(i->scale[i]*p[i]+shift[i],2)
end
function _occ_pcurve_d1(pc::NamedTuple{(:brep_map2d,)},t::Float64)
    base,a,b,scale,shift=pc.brep_map2d
    p,d=_occ_pcurve_d1(base,a*t+b)
    return ntuple(i->scale[i]*p[i]+shift[i],2),ntuple(i->a*scale[i]*d[i],2)
end
function _occ_pcurve_d2(pc::NamedTuple{(:brep_map2d,)},t::Float64)
    base,a,b,scale,shift=pc.brep_map2d
    p,d,dd=_occ_pcurve_d2(base,a*t+b)
    return ntuple(i->scale[i]*p[i]+shift[i],2),
        ntuple(i->a*scale[i]*d[i],2),ntuple(i->a*a*scale[i]*dd[i],2)
end

function _occ_brep_affine(matrix)
    linear=ntuple(i->matrix[div(i-1,3)+1,mod(i-1,3)+1],9)
    offset=ntuple(i->matrix[i,4],3)
    step=ntuple(i->matrix[div(i-1,4)+1,mod(i-1,4)+1],16)
    return _AffineTransform(linear,offset,[step],(linear=linear,loc=offset))
end

function _occ_brep_chain_uses(kernel,uses,caller)
    isempty(uses) && throw(ArgumentError("$caller: empty BREP wire"))
    endpoints(use)=begin
        edge=kernel.shapes[use[1]]
        use[2] ? (edge.children[2],edge.children[1]) :
                 (edge.children[1],edge.children[2])
    end
    heads=first.(endpoints.(uses)); tails=last.(endpoints.(uses))
    adjacency=Dict{UInt64,Vector{Int}}()
    balance=Dict{UInt64,Int}()
    for i in eachindex(uses)
        push!(get!(adjacency,heads[i],Int[]),i)
        balance[heads[i]]=get(balance,heads[i],0)+1
        balance[tails[i]]=get(balance,tails[i],0)-1
    end
    starts=[v for (v,n) in balance if n==1]
    valid=all(n->abs(n)<=1,values(balance)) && length(starts)<=1 &&
        count(==(-1),values(balance))==length(starts)
    valid || throw(ArgumentError("$caller: BREP wire has inconsistent edge orientations"))
    start=isempty(starts) ? heads[1] : starts[1]
    # Hierholzer traversal also handles the pole/rim loops attached to a
    # seam. Greedily closing that seam first can strand a valid rim loop.
    cursors=Dict{UInt64,Int}()
    vertices=UInt64[start]; edges=Int[]; trail=Int[]
    while !isempty(vertices)
        vertex=last(vertices)
        indices=get(adjacency,vertex,Int[])
        cursor=get(cursors,vertex,1)
        if cursor<=length(indices)
            index=indices[cursor]; cursors[vertex]=cursor+1
            push!(edges,index); push!(vertices,tails[index])
        else
            pop!(vertices)
            isempty(edges) || push!(trail,pop!(edges))
        end
    end
    length(trail)==length(uses) || throw(ArgumentError(
        "$caller: BREP wire has disconnected edge components"))
    return uses[reverse(trail)]
end

function _occ_brep_curve_geometry(curve,lo,hi,caller)
    curve=_occ_brep_untrim(curve)
    lo<hi || throw(ArgumentError("$caller: BREP edge range must increase"))
    if curve.kind==1
        return (occ=:line,origin=curve.origin,dir=curve.dir,t0=lo,t1=hi)
    elseif curve.kind==2
        return (occ=:circle,center=curve.center,n=curve.normal,X=curve.X,Y=curve.Y,
                r=curve.a,t0=lo,t1=hi)
    elseif curve.kind==3
        return (occ=:ellipse,center=curve.center,n=curve.normal,X=curve.X,Y=curve.Y,
                rx=curve.a,ry=curve.b,t0=lo,t1=hi)
    elseif curve.kind in (6,7)
        nc=_occ_brep_nurbs_curve(curve,caller)
        cv=_ModelOCCNURBSCurve(curve.kind==6 ? :bezier : :bspline,nc)
        return (occ=:copied_curve,curve=cv,t0=lo,t1=hi)
    elseif curve.kind==4
        f=curve.a
        f>0 || throw(ArgumentError("$caller: BREP parabola focal length must be positive"))
        value=t->ntuple(i->curve.center[i]+(t*t/(4f))*curve.X[i]+t*curve.Y[i],3)
        a=value(lo); b=value(hi); middle=value((lo+hi)/2)
        pole=ntuple(i->2middle[i]-(a[i]+b[i])/2,3)
        nc=NURBS.NURBSCurve(2,[lo,lo,lo,hi,hi,hi],[a,pole,b])
        return (occ=:copied_curve,curve=_ModelOCCNURBSCurve(:bspline,nc),t0=lo,t1=hi)
    end
    throw(ArgumentError("$caller: native BREP 3-D curve kind $(curve.kind) is not implemented"))
end

function _occ_brep_generatrix(curve,lo,hi,caller)
    curve=_occ_brep_untrim(curve)
    if curve.kind in (2,3)
        values,scales,weights,U,reparam=_occ_nurbs_conic_axis(lo,hi)
        ry=curve.kind==2 ? curve.a : curve.b
        poles=NTuple{3,Float64}[ntuple(i->curve.center[i]+scales[j]*(
            curve.a*_gm_cos(t)*curve.X[i]+ry*_gm_sin(t)*curve.Y[i]),3)
            for (j,t) in enumerate(values)]
        return NURBS.NURBSCurve(2,U,poles,weights),reparam
    elseif curve.kind==1
        poles=NTuple{3,Float64}[curve.origin.+t.*curve.dir for t in (lo,hi)]
        return NURBS.NURBSCurve(1,[lo,lo,hi,hi],poles),nothing
    end
    g=_occ_brep_curve_geometry(curve,lo,hi,caller)
    g.occ===:copied_curve || throw(ArgumentError(
        "$caller: BREP generatrix has no native rational representation"))
    return g.curve.curve,nothing
end

function _occ_brep_surface_geometry(surface,uvb,caller)
    if surface.kind==10
        return _occ_brep_surface_geometry(surface.base,uvb,caller)
    elseif surface.kind==11
        base=surface.base
        while base.kind==10; base=base.base; end
        d=surface.distance
        if base.kind==1
            base=merge(base,(center=base.center.+d.*base.axis,))
        elseif base.kind in (2,4)
            base=merge(base,(a=base.a+d,))
        elseif base.kind==3
            s,co=_gm_sincos(base.b)
            base=merge(base,(center=base.center.-(d*s).*base.axis,a=base.a+d*co))
        elseif base.kind==5
            base=merge(base,(b=base.b+d,))
        else
            throw(ArgumentError("$caller: native BREP offset surface kind $(base.kind) is not implemented"))
        end
        return _occ_brep_surface_geometry(base,uvb,caller)
    elseif surface.kind in (6,7)
        lo,hi=uvb
        if surface.kind==6
            gen,reparam=_occ_brep_generatrix(surface.curve,lo[1],hi[1],caller)
            n=length(gen.controls)
            C=Matrix{NTuple{3,Float64}}(undef,n,2)
            W=Matrix{Float64}(undef,n,2)
            for i in 1:n,j in 1:2
                v=j==1 ? lo[2] : hi[2]
                C[i,j]=gen.controls[i].+v.*surface.direction
                W[i,j]=gen.weights[i]
            end
            ns=NURBS.NURBSSurface(gen.degree,1,gen.knots,
                [lo[2],lo[2],hi[2],hi[2]],C,W)
            reparams=(u=reparam,v=nothing)
        else
            gen,reparam=_occ_brep_generatrix(surface.curve,lo[2],hi[2],caller)
            ns=_occ_iges_revolve_surface(gen,surface.center,surface.axis,
                lo[1],hi[1],caller)
            _,_,_,_,urep=_occ_nurbs_conic_axis(lo[1],hi[1])
            reparams=(u=urep,v=reparam)
        end
        g=(occ=:tensor_bspline,patch=_ModelNURBSSurface(:bspline,ns),
            reparam=reparams,reversed=false)
        return (g,(1.0,1.0),(0.0,0.0))
    elseif surface.kind in (8,9)
        poles=surface.poles; weights=surface.weights
        U=surface.U; V=surface.V
        if surface.periodic_u
            columns=[_occ_brep_nurbs_curve((kind=7,degree=surface.degree_u,
                periodic=true,poles=collect(poles[:,j]),weights=collect(weights[:,j]),
                knots=U),caller) for j in axes(poles,2)]
            poles=hcat((column.controls for column in columns)...)
            weights=hcat((column.weights for column in columns)...)
            U=columns[1].knots
        end
        if surface.periodic_v
            rows=[_occ_brep_nurbs_curve((kind=7,degree=surface.degree_v,
                periodic=true,poles=collect(poles[i,:]),weights=collect(weights[i,:]),
                knots=V),caller) for i in axes(poles,1)]
            poles=permutedims(hcat((row.controls for row in rows)...))
            weights=permutedims(hcat((row.weights for row in rows)...))
            V=rows[1].knots
        end
        ns=NURBS.NURBSSurface(surface.degree_u,surface.degree_v,U,V,poles,weights)
        return ((occ=:tensor_bspline,patch=_ModelNURBSSurface(:bspline,ns),
                 reversed=false),(1.0,1.0),(0.0,0.0))
    end
    surface.kind in 1:5 || throw(ArgumentError(
        "$caller: native BREP surface kind $(surface.kind) is not implemented"))
    C,A,X,Y=surface.center,surface.axis,surface.X,surface.Y
    common=(center=C,axis=A,X=X,Y=Y,reversed=false)
    if surface.kind==1
        return (merge(common,(occ=:plane,uvb=uvb,)),(1.0,1.0),(0.0,0.0))
    elseif surface.kind==2
        v0,v1=uvb[1][2],uvb[2][2]
        g=merge(common,(occ=:cylinder,center=C.+v0.*A,
                        radius=surface.a,height=v1-v0))
        return (g,(1.0,1.0),(0.0,-v0))
    elseif surface.kind==3
        v0,v1=uvb[1][2],uvb[2][2]
        alpha=surface.b
        s,co=_gm_sincos(alpha)
        g=merge(common,(occ=:cone,center=C.+(v0*co).*A,
            r1=surface.a+v0*s,r2=surface.a+v1*s,height=(v1-v0)*co))
        return (g,(1.0,1.0),(0.0,-v0))
    elseif surface.kind==4
        return (merge(common,(occ=:sphere,radius=surface.a)),(1.0,1.0),(0.0,0.0))
    end
    return (merge(common,(occ=:torus,r1=surface.a,r2=surface.b,angle=_OCC_TWO_PI)),
            (1.0,1.0),(0.0,0.0))
end

_model_occ_import_brep(m::GeoModel,path;kwargs...)=
    _model_occ_import_brep_source(m,read(path,String);kwargs...)

function _model_occ_import_brep_source(m::GeoModel,source::String;
        highest_dim_only=true,caller="OCC.importShapes",topology_receipt=nothing)
    parsed=_occ_brep_parse(source,caller)
    capacity=100000+64min(ncodeunits(source),14062)
    # This exact dynamic-program count detects compact compound DAGs before
    # walking their expanded occurrence tree or copying any live state.
    work=zeros(Int,length(parsed.shapes))
    for index in length(parsed.shapes):-1:1
        preflight_shape=parsed.shapes[index]
        count=1
        if preflight_shape.kind in ("Co","CS","Sh","Wi")
            for child in preflight_shape.children
                count=min(capacity+1,count+work[child.shape])
                count>capacity && break
            end
        end
        work[index]=count
    end
    work[parsed.root.shape]<=capacity || throw(ArgumentError(
        "$caller: expanded BREP topology exceeds import capacity ($capacity occurrences)"))
    scratch=GeoModel()
    scratch.occ_kernel=deepcopy(m.occ_kernel)
    kernel=scratch.occ_kernel
    for d in 0:3
        kernel.maxima[d+3]=max(kernel.maxima[d+3],m.next_tag[d+1])
    end
    isempty(m.loops) || (kernel.maxima[2]=max(kernel.maxima[2],maximum(keys(m.loops))))
    isempty(m.surface_loops) ||
        (kernel.maxima[1]=max(kernel.maxima[1],maximum(keys(m.surface_loops))))
    identity=Matrix{Float64}(LinearAlgebra.I,4,4)
    cache=Dict{Tuple{Int,Tuple,Bool},UInt64}()
    sources=Dict{UInt64,Tuple{Int,Matrix{Float64},Tuple}}()
    parameter_maps=Dict{UInt64,Tuple{Float64,Float64}}()
    # Compound DAGs can describe exponentially many required occurrences.
    # Bound expanded work/output independently of recursion depth and avoid
    # merging distinct source instances to enforce the bound.
    remaining=Ref(capacity)
    consume()=begin
        remaining[]>0 || throw(ArgumentError(
            "$caller: expanded BREP topology exceeds import capacity ($capacity occurrences)"))
        remaining[]-=1
        nothing
    end
    located(key,index)=index==0 ? key : _occ_brep_location_product(key,
        _occ_brep_lookup(parsed.location_words,index,caller,"location"),caller)
    matrices=Dict{Tuple,Matrix{Float64}}(()=>identity)
    matrix_of(key)=get!(matrices,key) do
        _occ_brep_location_matrix(key,parsed.locations,caller)
    end
    placed(key,index)=matrix_of(located(key,index))
    function bind(d,geometry,children)
        tag=_model_occ_next_tag(kernel,d,-1,caller)
        id=only(_model_occ_reserve_shapes(kernel,1,caller))
        _model_occ_bind_shape!(kernel,d,tag,id,_ModelOCCShape(d,geometry,children))
        kernel.next_shape=id+UInt64(1)
        return id
    end
    function realize(index,matrix,location_key=(),reversed=false,depth=0)
        depth<256 || throw(ArgumentError("$caller: BREP topology nesting exceeds 256"))
        shape=_occ_brep_lookup(parsed.shapes,index,caller,"shape")
        # Edge/face orientation belongs to its parent's use, while wires
        # and shells store oriented use lists themselves.
        oriented=shape.kind in ("Wi","Sh") && reversed
        key=(index,location_key,oriented)
        haskey(cache,key) && return cache[key]
        consume()
        child(ref)=realize(ref.shape,placed(location_key,ref.location),
            located(location_key,ref.location),ref.orientation=='-',depth+1)
        data=shape.data
        id=if shape.kind=="Ve"
            isempty(shape.children) || throw(ArgumentError("$caller: BREP vertex has children"))
            p=_occ_trsf_apply(_occ_brep_affine(matrix),data.point)
            bind(0,(point=p,mesh_size=0.0),UInt64[])
        elseif shape.kind=="Ed"
            vertices=[(child(ref),ref.orientation) for ref in shape.children
                      if parsed.shapes[ref.shape].kind=="Ve"]
            isempty(vertices) && throw(ArgumentError("$caller: BREP edge has no vertices"))
            va=findfirst(v->v[2]=='+',vertices); vb=findfirst(v->v[2]=='-',vertices)
            a=vertices[something(va,1)][1]; b=vertices[something(vb,length(vertices))][1]
            records=data.representations
            ri=findfirst(r->r.kind==1,records)
            if data.degenerate
                r=findfirst(r->r.kind in (2,3),records)
                r===nothing && throw(ArgumentError("$caller: degenerate BREP edge lacks a pcurve"))
                g=(occ=:degenerate,t0=records[r].lo,t1=records[r].hi)
                result=bind(1,(kind=:degenerate,geometry=g),UInt64[a,b])
                parameter_maps[result]=(1.0,0.0)
                result
            else
                ri===nothing && throw(ArgumentError(
                    "$caller: native BREP edge without a 3-D curve is not implemented"))
                r=records[ri]
                curve=_occ_brep_lookup(parsed.curves3,r.curve,caller,"3-D curve")
                original=_occ_brep_curve_geometry(curve,r.lo,r.hi,caller)
                trsf=_occ_brep_affine(placed(location_key,r.location))
                g=_transform_occ_curve_record(scratch,0,original,trsf,caller,"BREP edge")
                result=bind(1,(kind=g.occ,geometry=g),UInt64[a,b])
                alpha=(r.hi-r.lo)/(g.t1-g.t0)
                reflected=original.occ in (:circle,:ellipse) &&
                    original.t1-original.t0>=_OCC_TWO_PI && _linear_det(trsf.linear)<0
                parameter_maps[result]=reflected ? (-1.0,0.0) :
                    (alpha,r.lo-alpha*g.t0)
                result
            end
        elseif shape.kind=="Wi"
            uses=Tuple{UInt64,Bool}[]
            for ref in shape.children
                parsed.shapes[ref.shape].kind=="Ed" || throw(ArgumentError(
                    "$caller: BREP wire child is not an edge"))
                push!(uses,(child(ref),ref.orientation=='-'))
            end
            uses=_occ_brep_chain_uses(kernel,uses,caller)
            oriented && (uses=[(eid,!rev) for (eid,rev) in reverse(uses)])
            bind(-1,(uses=uses,),UInt64[first(use) for use in uses])
        elseif shape.kind=="Fa"
            wire_ids=UInt64[child(ref) for ref in shape.children]
            all(id->kernel.shapes[id].dimension==-1,wire_ids) || throw(ArgumentError(
                "$caller: BREP face child is not a wire"))
            surface=_occ_brep_lookup(parsed.surfaces,data.surface,caller,"surface")
            pcurves=Dict{Int,NamedTuple}(); ranges=Dict{Int,Tuple{Float64,Float64}}()
            for wid in wire_ids,(eid,rev) in kernel.shapes[wid].geometry.uses
                ei,ematrix,edge_location=sources[eid]
                edge=parsed.shapes[ei]
                ri=findfirst(edge.data.representations) do r
                    r.kind in (2,3) && r.surface==data.surface &&
                        located(edge_location,r.location)==located(location_key,data.location)
                end
                if ri===nothing
                    plane,_,_=_occ_brep_surface_geometry(surface,([0.,0.],[1.,1.]),caller)
                    plane.occ===:plane || throw(ArgumentError(
                        "$caller: BREP face edge has no matching pcurve"))
                    plane=merge(plane,(pcurves=Dict{Int,NamedTuple}(),))
                    plane=_transform_occ_surface_record(plane,
                        _occ_brep_affine(placed(location_key,data.location)),caller,"BREP plane")
                    eg=kernel.shapes[eid].geometry.geometry
                    tag=_model_occ_shape_tag(kernel,eid,caller)
                    pcurves[tag]=(fwd=_occ_brep_plane_pcurve(eg,plane),rev=nothing)
                    ranges[tag]=(eg.t0,eg.t1)
                    continue
                end
                r=edge.data.representations[ri]
                first_pc=_occ_brep_pcurve(_occ_brep_lookup(parsed.curves2,r.curve,caller,"2-D curve"),caller)
                second_pc=r.kind==3 ? _occ_brep_pcurve(
                    _occ_brep_lookup(parsed.curves2,r.curve2,caller,"2-D curve"),caller) : nothing
                alpha,beta=parameter_maps[eid]
                eg=kernel.shapes[eid].geometry.geometry
                if !edge.data.same_range
                    curve_index=findfirst(x->x.kind==1,edge.data.representations)
                    if curve_index===nothing
                        alpha=(r.hi-r.lo)/(eg.t1-eg.t0)
                        beta=r.lo-alpha*eg.t0
                    else
                        curve_rep=edge.data.representations[curve_index]
                        factor=(r.hi-r.lo)/(curve_rep.hi-curve_rep.lo)
                        alpha,beta=factor*alpha,r.lo+factor*(beta-curve_rep.lo)
                    end
                end
                first_pc=_occ_brep_pcurve_map(first_pc,alpha,beta,(1.0,1.0),(0.0,0.0))
                second_pc===nothing || (second_pc=_occ_brep_pcurve_map(
                    second_pc,alpha,beta,(1.0,1.0),(0.0,0.0)))
                tag=_model_occ_shape_tag(kernel,eid,caller)
                pcurves[tag]=(fwd=first_pc,rev=second_pc)
                ranges[tag]=(eg.t0,eg.t1)
            end
            natural=isempty(wire_ids)
            isempty(pcurves) && !natural && throw(ArgumentError(
                "$caller: native BREP face has no boundary curves"))
            uvb=natural ? _occ_brep_natural_bounds(surface,caller) :
                _model_occ_prim_uvb(pcurves,ranges)
            base,scale,shift=_occ_brep_surface_geometry(surface,uvb,caller)
            base=merge(base,(pcurves=pcurves,))
            geometry=_transform_occ_surface_record(base,
                _occ_brep_affine(placed(location_key,data.location)),caller,"BREP face")
            vscale=_model_occ_pcurve_vscale(base,geometry)
            actual_pc=Dict{Int,NamedTuple}()
            for (tag,entry) in pcurves
                map_pc=pc->pc===nothing ? nothing : _occ_brep_pcurve_map(
                    pc,1.0,0.0,(scale[1],vscale*scale[2]),(shift[1],vscale*shift[2]))
                actual_pc[tag]=(fwd=map_pc(entry.fwd),rev=map_pc(entry.rev))
            end
            remapped=Set{UInt64}()
            for wid in wire_ids,(eid,_) in kernel.shapes[wid].geometry.uses
                eid in remapped && continue
                push!(remapped,eid)
                ei,_,_=sources[eid]
                edge=parsed.shapes[ei]
                edge.data.same_parameter && continue
                eg=kernel.shapes[eid].geometry.geometry
                eg.occ===:degenerate && continue
                tag=_model_occ_shape_tag(kernel,eid,caller)
                entry=actual_pc[tag]
                if geometry.occ===:plane
                    actual_pc[tag]=(fwd=_occ_brep_plane_pcurve(eg,geometry),rev=nothing)
                else
                    tolerance=8max(edge.data.tolerance,data.tolerance,_OCC_CONFUSION)
                    # The old trace supplies a branch hint, while the native
                    # 3-D curve and actual face supply the exact point map.
                    wrap=pc->pc===nothing ? nothing : (brep_inverse2d=(eg,geometry,pc,
                        eg.t0,eg.t1,tolerance),)
                    actual_pc[tag]=(fwd=wrap(entry.fwd),rev=wrap(entry.rev))
                end
            end
            geometry=merge(geometry,(pcurves=actual_pc,reversed=reversed,))
            if natural
                transform_uv=p->(p[1]*scale[1]+shift[1],
                    vscale*(p[2]*scale[2]+shift[2]))
                geometry=merge(geometry,(natural_bounds=(
                    transform_uv(uvb[1]),transform_uv(uvb[2])),))
                geometry=merge(geometry,(natural_collapsed=
                    _occ_natural_collapsed_sides(geometry),))
            end
            bounds=natural ? geometry.natural_bounds :
                _model_occ_prim_uvb(actual_pc,ranges)
            canonical=_occ_surface_bounds(geometry)
            uperiod=_brep_uperiod(geometry)
            vperiod=_occ_surface_isvperiodic(geometry) ? _OCC_TWO_PI : 0.0
            full_u=uperiod>0 ? bounds[2][1]-bounds[1][1]>=uperiod-1e-9 :
                bounds[1][1]<=canonical[1][1]+1e-9 && bounds[2][1]>=canonical[2][1]-1e-9
            full_v=vperiod>0 ? bounds[2][2]-bounds[1][2]>=vperiod-1e-9 :
                bounds[1][2]<=canonical[1][2]+1e-9 && bounds[2][2]>=canonical[2][2]-1e-9
            trimmed=length(wire_ids)>1 || !(full_u && full_v)
            geometry=merge(geometry,(trimmed=trimmed,))
            geometry.occ==:plane && (geometry=merge(geometry,
                (uvb=bounds,)))
            vuv=Dict{Int,NTuple{2,Float64}}()
            for wid in wire_ids,(eid,rev) in kernel.shapes[wid].geometry.uses
                eg=kernel.shapes[eid].geometry.geometry
                tag=_model_occ_shape_tag(kernel,eid,caller)
                entry=actual_pc[tag]
                pc=rev && entry.rev!==nothing ? entry.rev : entry.fwd
                for (vid,t) in zip(kernel.shapes[eid].children,(eg.t0,eg.t1))
                    get!(vuv,_model_occ_shape_tag(kernel,vid,caller),_occ_pcurve_d0(pc,t))
                end
                if eg.occ===:degenerate
                    kernel.shapes[eid]=_ModelOCCShape(1,
                        (kind=:degenerate,geometry=merge(eg,(face=geometry,))),
                        kernel.shapes[eid].children)
                end
            end
            _model_nurbs_surface(geometry) && (geometry=merge(geometry,(vertex_uv=vuv,)))
            bind(2,(kind=geometry.occ,geometry=geometry),wire_ids)
        elseif shape.kind=="Sh"
            uses=Tuple{UInt64,Bool}[]
            for ref in shape.children
                use_reversed=xor(ref.orientation=='-',oriented)
                fid=realize(ref.shape,placed(location_key,ref.location),
                    located(location_key,ref.location),use_reversed,depth+1)
                kernel.shapes[fid].dimension==2 || throw(ArgumentError(
                    "$caller: BREP shell child is not a face"))
                push!(uses,(fid,use_reversed))
            end
            bind(-2,(uses=uses,),UInt64[first(use) for use in uses])
        elseif shape.kind=="So"
            children=UInt64[realize(ref.shape,placed(location_key,ref.location),
                located(location_key,ref.location),xor(ref.orientation=='-',reversed),
                depth+1) for ref in shape.children]
            all(id->kernel.shapes[id].dimension==-2,children) || throw(ArgumentError(
                "$caller: BREP solid child is not a shell"))
            bind(3,(kind=:solid,),children)
        else
            throw(ArgumentError("$caller: internal BREP compound was realized as an entity"))
        end
        sources[id]=(index,matrix,location_key)
        cache[key]=id
        return id
    end
    roots=Tuple{Int,Int}[]
    function visit(ref,matrix,location_key=(),parent_reversed=false,depth=0)
        depth<256 || throw(ArgumentError("$caller: BREP compound nesting exceeds 256"))
        consume()
        location_key=located(location_key,ref.location)
        matrix=matrix_of(location_key)
        reversed=xor(parent_reversed,ref.orientation=='-')
        shape=parsed.shapes[ref.shape]
        if shape.kind in ("Co","CS","Sh","Wi")
            for sub in shape.children; visit(sub,matrix,location_key,reversed,depth+1); end
        else
            id=realize(ref.shape,matrix,location_key,reversed)
            push!(roots,kernel.shape_tags[id])
        end
    end
    visit(parsed.root,identity)
    isempty(roots) && throw(ArgumentError("$caller: BREP has no transferable native geometry"))
    if !highest_dim_only
        counts=Dict{UInt64,Int}()
        function count_result(id)
            get!(counts,id) do
                shape=kernel.shapes[id]
                subtotal=Int(shape.dimension>=0)
                for child in shape.children
                    child_count=count_result(child)
                    child_count<=capacity-subtotal || throw(ArgumentError(
                        "$caller: expanded BREP result exceeds import capacity ($capacity occurrences)"))
                    subtotal+=child_count
                end
                return subtotal
            end
        end
        returned_count=0
        for root in roots
            root_count=count_result(kernel.bindings[root])
            root_count<=capacity-returned_count || throw(ArgumentError(
                "$caller: expanded BREP result exceeds import capacity ($capacity occurrences)"))
            returned_count+=root_count
        end
    end
    result=_occ_import_results(scratch,roots,highest_dim_only,caller)
    references=Dict{Int,Tuple{Int,Int}}()
    if topology_receipt!==nothing
        # XAO uses TopExp_MapShapes indices: first-visit preorder over all
        # topological types, rather than BREP's reverse storage numbers.
        seen=Set{Tuple{Int,Tuple}}()
        counter=Ref(0)
        function record(ref,matrix,location_key=())
            location_key=located(location_key,ref.location)
            matrix=matrix_of(location_key)
            key=(ref.shape,location_key)
            key in seen && return nothing
            push!(seen,key); counter[]+=1
            shape=parsed.shapes[ref.shape]
            if shape.kind in ("Ve","Ed","Fa","So")
                id=cache[(ref.shape,location_key,false)]
                references[counter[]]=kernel.shape_tags[id]
            end
            for sub in shape.children; record(sub,matrix,location_key); end
            return nothing
        end
        record(parsed.root,identity)
    end
    m.occ_kernel=kernel
    topology_receipt===nothing || (topology_receipt[]=references)
    return result
end

function _occ_brep_retag_entities!(kernel,targets,imported;caller="XAO import")
    old=Dict{UInt64,Tuple{Int,Int}}()
    occupied=Set{Tuple{Int,Int}}()
    for (id,target) in targets
        id in imported || throw(ArgumentError("$caller: topology references retained CAD geometry"))
        haskey(kernel.shape_tags,id) || throw(ArgumentError("$caller: topology references unbound geometry"))
        dimension,tag=target
        dimension==kernel.shapes[id].dimension && 1<=tag<=typemax(Int32) ||
            throw(ArgumentError("$caller: invalid topology entity tag"))
        target in occupied && throw(ArgumentError("$caller: duplicate topology entity tag"))
        push!(occupied,target)
        previous=get(kernel.bindings,target,nothing)
        previous===nothing || haskey(targets,previous) || throw(ArgumentError(
            "$caller: topology entity tag collides with retained geometry"))
        old[id]=kernel.shape_tags[id]
    end
    for (id,(dimension,tag)) in old
        delete!(kernel.bindings,(dimension,tag))
        _model_occ_entity_map_unbind!(kernel.entity_maps[dimension+3],tag)
    end
    curves=Dict{Int,Int}(); points=Dict{Int,Int}()
    for (id,(dimension,tag)) in targets
        previous=old[id][2]
        kernel.bindings[(dimension,tag)]=id
        kernel.shape_tags[id]=(dimension,tag)
        _model_occ_entity_map_bind!(kernel.entity_maps[dimension+3],tag)
        kernel.maxima[dimension+3]=max(kernel.maxima[dimension+3],tag)
        dimension==0 && (points[previous]=tag)
        dimension==1 && (curves[previous]=tag)
    end
    shapes=Dict(id=>kernel.shapes[id] for id in imported)
    _model_occ_rekey_pcurves!(kernel,shapes,curves;tagmap0=points,caller=caller)
    kernel.changed=true
    return nothing
end

include("ModelOCCBREPTrace.jl")
include("ModelOCCBREPNatural.jl")
