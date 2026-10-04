# Included inside QuadTriNoNewRectGridCertificates. Literal twice-reference
# coordinates are pinned Gmsh4.15.2 MSH11--14 metadata, separate from production
# basis/elevation/certification code. All physical computations use exact
# rational representations of the actual Float64 payload.
const _P2_REFERENCE_TWICE=Dict(
    11=>[(0,0,0),(2,0,0),(0,2,0),(0,0,2),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(0,1,1),(1,0,1)],
    12=>[(-2,-2,-2),(2,-2,-2),(2,2,-2),(-2,2,-2),(-2,-2,2),(2,-2,2),(2,2,2),(-2,2,2),
        (0,-2,-2),(-2,0,-2),(-2,-2,0),(2,0,-2),(2,-2,0),(0,2,-2),(2,2,0),(-2,2,0),
        (0,-2,2),(-2,0,2),(2,0,2),(0,2,2),(0,0,-2),(0,-2,0),(-2,0,0),(2,0,0),(0,2,0),(0,0,2),(0,0,0)],
    13=>[(0,0,-2),(2,0,-2),(0,2,-2),(0,0,2),(2,0,2),(0,2,2),(1,0,-2),(0,1,-2),(0,0,0),
        (1,1,-2),(2,0,0),(0,2,0),(1,0,2),(0,1,2),(1,1,2),(1,0,0),(0,1,0),(1,1,0)],
    14=>[(-2,-2,0),(2,-2,0),(2,2,0),(-2,2,0),(0,0,2),(0,-2,0),(-2,0,0),(-1,-1,1),
        (2,0,0),(1,-1,1),(0,2,0),(1,1,1),(-1,1,1),(0,0,0)])
const _P2_POLYNOMIAL=Dict{NTuple{3,Int},Q}
const _P2_ZERO_POWER=(0,0,0)

function _p2_add(polynomials...)
    result=_P2_POLYNOMIAL()
    for polynomial in polynomials,(power,value) in polynomial
        result[power]=get(result,power,zero(Q))+value
    end
    filter!(pair->!iszero(last(pair)),result)
    return result
end
_p2_scale(polynomial,value)=_P2_POLYNOMIAL(power=>coefficient*value for (power,coefficient) in polynomial if !iszero(coefficient*value))
function _p2_multiply(a,b)
    result=_P2_POLYNOMIAL()
    for (first_power,x) in a,(second_power,y) in b
        power=ntuple(k->first_power[k]+second_power[k],3)
        result[power]=get(result,power,zero(Q))+x*y
    end
    filter!(pair->!iszero(last(pair)),result)
    return result
end
function _p2_power(polynomial,n)
    result=_P2_POLYNOMIAL(_P2_ZERO_POWER=>one(Q))
    for _ in 1:n;result=_p2_multiply(result,polynomial);end
    return result
end
function _p2_derivative(polynomial,axis)
    return _P2_POLYNOMIAL(ntuple(k->power[k]-(k==axis),3)=>value*power[axis]
        for (power,value) in polynomial if power[axis]>0)
end
function _p2_compose(polynomial,replacements)
    result=_P2_POLYNOMIAL()
    for (power,value) in polynomial
        term=_P2_POLYNOMIAL(_P2_ZERO_POWER=>value)
        for axis in 1:3;term=_p2_multiply(term,_p2_power(replacements[axis],power[axis]));end
        result=_p2_add(result,term)
    end
    return result
end
_p2_determinant(a,b,c)=_p2_add(
    _p2_multiply(a[1],_p2_add(_p2_multiply(b[2],c[3]),_p2_scale(_p2_multiply(b[3],c[2]),-1))),
    _p2_scale(_p2_multiply(a[2],_p2_add(_p2_multiply(b[1],c[3]),_p2_scale(_p2_multiply(b[3],c[1]),-1))),-1),
    _p2_multiply(a[3],_p2_add(_p2_multiply(b[1],c[2]),_p2_scale(_p2_multiply(b[2],c[1]),-1))))

function _p2_inverse(matrix)
    n=size(matrix,1);work=hcat(matrix,[Q(i==j) for i in 1:n,j in 1:n])
    for column in 1:n
        pivot=findfirst(row->!iszero(work[row,column]),column:n)
        require(pivot!==nothing,"singular literal nodal Vandermonde")
        row=column+Int(pivot)-1
        if row!=column;work[[column,row],:]=work[[row,column],:];end
        work[column,:]./=work[column,column]
        for other in 1:n
            other==column && continue
            scale=work[other,column]
            iszero(scale) || (work[other,:].-=scale.*work[column,:])
        end
    end
    return work[:,n+1:2n]
end

function _p2_reference(msh)
    points=[Tuple(Q(value)/2 for value in point) for point in _P2_REFERENCE_TWICE[msh]]
    powers=NTuple{3,Int}[Tuple(power) for power in Iterators.product(0:2,0:2,0:2)
        if msh==12 || msh==11 && sum(power)<=2 || msh==13 && power[1]+power[2]<=2 ||
            msh==14 && power[3]<=2-max(power[1],power[2])]
    require(length(points)==length(powers),"literal nodal space/reference width differs")
    function monomial(p,power)
        msh==14 && p[3]==1 && max(power[1],power[2])>0 && return zero(Q)
        value=prod(p[axis]^power[axis] for axis in 1:3)
        return msh==14 ? value/(1-p[3])^min(power[1],power[2]) : value
    end
    vandermonde=[monomial(p,power) for p in points,power in powers]
    return (;points,powers,inverse=_p2_inverse(vandermonde))
end
const _P2_SPACES=Dict(msh=>_p2_reference(msh) for msh in (11,12,13,14))

function _p2_pyramid_term(i,j,k,h;scale=1)
    require(min(i,j,k,h)>=0,"negative collapsed polynomial exponent")
    return _P2_POLYNOMIAL((i,j,k+degree)=>Q(scale)*binomial(h,degree)*(-1)^degree for degree in 0:h)
end
function _p2_map_derivatives(msh,vertices)
    space=_P2_SPACES[msh]
    anchored=[Q(vertices[node][dimension])-Q(vertices[1][dimension]) for node in eachindex(vertices),dimension in 1:3]
    coefficients=space.inverse*anchored
    if msh!=14
        maps=[_P2_POLYNOMIAL(power=>coefficients[row,dimension] for (row,power) in enumerate(space.powers)
            if !iszero(coefficients[row,dimension])) for dimension in 1:3]
        return [[_p2_derivative(map,axis) for map in maps] for axis in 1:3]
    end
    derivatives=[[_P2_POLYNOMIAL() for _ in 1:3] for _ in 1:3]
    for (row,(i,j,k)) in enumerate(space.powers)
        height=max(i,j)
        du=i>0 ? _p2_pyramid_term(i-1,j,k,height-1;scale=i) : _P2_POLYNOMIAL()
        dv=j>0 ? _p2_pyramid_term(i,j-1,k,height-1;scale=j) : _P2_POLYNOMIAL()
        dw=_p2_add(k>0 ? _p2_pyramid_term(i,j,k-1,height;scale=k) : _P2_POLYNOMIAL(),
            min(i,j)>0 ? _p2_pyramid_term(i,j,k,height-1;scale=min(i,j)) : _P2_POLYNOMIAL())
        for (axis,term) in enumerate((du,dv,dw)),dimension in 1:3
            derivatives[axis][dimension]=_p2_add(derivatives[axis][dimension],_p2_scale(term,coefficients[row,dimension]))
        end
    end
    return derivatives
end

function _p2_bernstein(polynomial)
    degrees=ntuple(axis->isempty(polynomial) ? 0 : maximum(power[axis] for power in keys(polynomial)),3)
    coefficients=Q[]
    for index in Iterators.product((0:degree for degree in degrees)...)
        value=zero(Q)
        for (power,coefficient) in polynomial
            if all(power[axis]<=index[axis] for axis in 1:3)
                value+=coefficient*prod(Q(binomial(index[axis],power[axis]))/binomial(degrees[axis],power[axis]) for axis in 1:3)
            end
        end
        push!(coefficients,value)
    end
    return degrees,coefficients
end

function quadratic_map_certificate(msh,vertices)
    require(haskey(_P2_SPACES,msh) && length(vertices)==length(_P2_REFERENCE_TWICE[msh]),"unexpected quadratic family/width")
    derivatives=_p2_map_derivatives(msh,vertices);polynomial=_p2_determinant(derivatives...)
    unit=_P2_POLYNOMIAL(_P2_ZERO_POWER=>one(Q));variables=[_P2_POLYNOMIAL(ntuple(k->Int(k==axis),3)=>one(Q)) for axis in 1:3]
    if msh in (11,13)
        r=variables[1];s=_p2_multiply(_p2_add(unit,_p2_scale(r,-1)),variables[2])
        t=msh==11 ? _p2_multiply(_p2_multiply(_p2_add(unit,_p2_scale(r,-1)),_p2_add(unit,_p2_scale(variables[2],-1))),variables[3]) :
            _p2_add(_p2_scale(variables[3],2),_p2_scale(unit,-1))
        replacements=(r,s,t)
    else
        replacements=Tuple(_p2_add(_p2_scale(variable,2),_p2_scale(unit,-1)) for variable in variables)
        msh==14 && (replacements=(replacements[1],replacements[2],variables[3]))
    end
    degrees,coefficients=_p2_bernstein(_p2_compose(polynomial,replacements))
    require(!isempty(coefficients) && minimum(coefficients)>0,"whole quadratic reference determinant is not certified positive")
    volume=zero(Q)
    for ((i,j,k),coefficient) in polynomial
        weight=if msh==11
            Q(factorial(big(i))*factorial(big(j))*factorial(big(k)))/factorial(big(i+j+k+3))
        elseif msh==13
            Q(factorial(big(i))*factorial(big(j)))/factorial(big(i+j+2))*(iseven(k) ? Q(2)/(k+1) : zero(Q))
        elseif msh==12
            prod(iseven(power) ? Q(2)/(power+1) : zero(Q) for power in (i,j,k))
        else
            (iseven(i) ? Q(2)/(i+1) : zero(Q))*(iseven(j) ? Q(2)/(j+1) : zero(Q))*
                Q(2factorial(big(k)))/factorial(big(k+3))
        end
        volume+=coefficient*weight
    end
    require(volume>0,"nonpositive exact quadratic map volume")
    return (;degrees,coefficients=length(coefficients),minimum=minimum(coefficients),maximum=maximum(coefficients),volume)
end

function certify_quadratic(mesh,f,source)
    primary_counts=Dict(11=>4,12=>8,13=>6,14=>5)
    primary=sort!(unique(Int(node) for block in blocks(mesh) for cell in eachcol(block.nodes)
        for node in cell[1:primary_counts[Int(block.msh)]]))
    positions=Dict(node=>position for (position,node) in enumerate(primary))
    linear_blocks=[ElementBlock(Int(block.msh)-7,Int32[positions[Int(block.nodes[row,column])]
        for row in 1:primary_counts[Int(block.msh)],column in axes(block.nodes,2)]) for block in blocks(mesh)]
    linear=MixedMesh(mesh.coords[:,primary],linear_blocks)
    underlying=certify(linear,f,source)
    support=Quad.certify_quadratic(mesh,linear)
    maps=NamedTuple[]
    for block in blocks(mesh),cell in eachcol(block.nodes)
        push!(maps,quadratic_map_certificate(Int(block.msh),Tuple(point(mesh,node) for node in cell)))
    end
    require(length(maps)==underlying.ncell,"quadratic cell set changed")
    return merge(underlying,(;total=sum(map.volume for map in maps),maps,nsupport=support.nsupport,linear))
end
