# Certify a tensor seam over its complete knot domain. Coordinate samples
# cannot establish a polynomial identity, and a world-position-dependent
# tolerance would change topology when a model is translated.
const _NURBSSeamRational=Rational{BigInt}
const _NURBSSeamInterpolation=Dict{Int,Matrix{_NURBSSeamRational}}()
const _NURBSSeamInterpolationLock=ReentrantLock()

function _model_nurbs_seam_basis(knots,degree,count,parameter,
        span=clamp(searchsortedlast(knots,parameter)-1,degree,count-1))
    basis=zeros(_NURBSSeamRational,degree+1)
    left=zeros(_NURBSSeamRational,degree)
    right=similar(left);basis[1]=1
    for level in 1:degree
        left[level]=parameter-knots[span+2-level]
        right[level]=knots[span+level+1]-parameter
        saved=zero(_NURBSSeamRational)
        for index in 0:level-1
            denominator=right[index+1]+left[level-index]
            term=iszero(denominator) ? zero(saved) : basis[index+1]/denominator
            basis[index+1]=saved+right[index+1]*term
            saved=left[level-index]*term
        end
        basis[level+1]=saved
    end
    return (span-degree+1:span+1),basis
end

function _model_nurbs_seam_iso(surface,axis,parameter)
    degree=axis==1 ? surface.degree_u : surface.degree_v
    knots=_NURBSSeamRational.(axis==1 ? surface.knots_u : surface.knots_v)
    indices,basis=_model_nurbs_seam_basis(knots,degree,
        size(surface.controls,axis),_NURBSSeamRational(parameter))
    count=size(surface.controls,3-axis)
    return [ntuple(component->sum(begin
        row,column=axis==1 ? (index,other) : (other,index)
        weight=_NURBSSeamRational(surface.weights[row,column])*basis[ordinal]
        component==4 ? weight : weight*_NURBSSeamRational(
            surface.controls[row,column][component])
    end for (ordinal,index) in enumerate(indices)),4) for other in 1:count]
end

function _model_nurbs_seam_evaluate(controls,degree,knots,parameter,span)
    indices,basis=_model_nurbs_seam_basis(knots,degree,length(controls),parameter,span)
    return ntuple(component->sum(basis[ordinal]*controls[index][component]
        for (ordinal,index) in enumerate(indices)),4)
end

function _model_nurbs_seam_interpolation(degree)
    lock(_NURBSSeamInterpolationLock) do
        get!(_NURBSSeamInterpolation,degree) do
            matrix=Matrix{_NURBSSeamRational}(undef,degree+1,degree+1)
            for row in 0:degree,column in 0:degree
                t=degree==0 ? zero(_NURBSSeamRational) :
                    _NURBSSeamRational(row,degree)
                matrix[row+1,column+1]=binomial(BigInt(degree),column)*
                    t^column*(1-t)^(degree-column)
            end
            inv(matrix)
        end
    end
end

function _model_nurbs_seam_bezier(controls,degree,knots,lower,upper)
    samples=Matrix{_NURBSSeamRational}(undef,degree+1,4)
    # Use this span's polynomial at both endpoints, including the left limit
    # of a repeated upper knot. A neighboring span can have another trace.
    span=clamp(searchsortedlast(knots,(lower+upper)/2)-1,degree,length(controls)-1)
    for index in 0:degree
        fraction=degree==0 ? zero(_NURBSSeamRational) :
            _NURBSSeamRational(index,degree)
        point=_model_nurbs_seam_evaluate(controls,degree,knots,
            lower+(upper-lower)*fraction,span)
        for component in 1:4;samples[index+1,component]=point[component];end
    end
    coefficients=_model_nurbs_seam_interpolation(degree)*samples
    return [ntuple(component->coefficients[index,component],4)
        for index in 1:degree+1]
end

function _model_nurbs_seam_split(controls)
    work=copy(controls);left=[first(work)];right=[last(work)]
    for level in length(work)-1:-1:1
        for index in 1:level
            work[index]=ntuple(component->
                (work[index][component]+work[index+1][component])/2,4)
        end
        push!(left,work[1]);push!(right,work[level])
    end
    return left,reverse!(right)
end

function _model_nurbs_seam_near(first,last,tolerance_squared)
    first[4]>0 && last[4]>0 || return false
    return sum((first[axis]/first[4]-last[axis]/last[4])^2
        for axis in 1:3)<=tolerance_squared
end

# A fixed parameter side is a constant rational curve exactly when its
# homogeneous control coefficients all dehomogenize to the same point.
# Derive that iso-curve from the active basis, including nonclamped ends.
function _model_nurbs_collapsed_iso_point(g,axis,value)
    reparam=get(g,:reparam,nothing)
    mapping=reparam===nothing ? nothing : axis==1 ? reparam.u : reparam.v
    parameter=mapping===nothing ? value : first(_occ_conic_reparam(mapping,value))
    coefficients=_model_nurbs_seam_iso(g.patch.surface,axis,parameter)
    point=first(coefficients)
    all(coefficient->_model_nurbs_seam_near(point,coefficient,0),coefficients) ||
        return nothing
    return ntuple(component->Float64(point[component]/point[4]),3)
end

function _model_nurbs_seam_panel(first,last,tolerance_squared,budget,depth=0)
    budget[]>0 || return false
    budget[]-=1
    _model_nurbs_seam_near(first[1],last[1],tolerance_squared) &&
        _model_nurbs_seam_near(first[end],last[end],tolerance_squared) || return false
    degree=length(first)-1
    numerators=zeros(_NURBSSeamRational,2degree+1,3)
    for i in 0:degree,j in 0:degree
        weight=_NURBSSeamRational(binomial(BigInt(degree),i)*
            binomial(BigInt(degree),j),binomial(BigInt(2degree),i+j))
        for axis in 1:3
            numerators[i+j+1,axis]+=weight*(first[i+1][axis]*last[j+1][4]-
                last[j+1][axis]*first[i+1][4])
        end
    end
    denominator=minimum(point[4] for point in first)*minimum(point[4] for point in last)
    bound=sum(maximum(abs,@view(numerators[:,axis]))^2 for axis in 1:3)
    denominator>0 && bound<=tolerance_squared*denominator^2 && return true
    depth<16 || return false
    fa,fb=_model_nurbs_seam_split(first);la,lb=_model_nurbs_seam_split(last)
    _model_nurbs_seam_near(fa[end],la[end],tolerance_squared) || return false
    return _model_nurbs_seam_panel(fa,la,tolerance_squared,budget,depth+1) &&
           _model_nurbs_seam_panel(fb,lb,tolerance_squared,budget,depth+1)
end

function _model_nurbs_certified_seam_period(g,axis)
    surface=g.patch.surface
    lo,hi=_model_nurbs_surface_bounds(g)
    width=hi[axis]-lo[axis]
    isfinite(width) && width>0 || return 0.0
    tolerance=Float64(get(g,:seam_tolerance,1e-7))
    isfinite(tolerance) && tolerance>=0 || return 0.0
    tolerance_squared=_NURBSSeamRational(tolerance)^2
    first=_model_nurbs_seam_iso(surface,axis,lo[axis])
    last=_model_nurbs_seam_iso(surface,axis,hi[axis])
    all(point->point[4]>0,first) && all(point->point[4]>0,last) || return 0.0
    # Proportional weights give a shared rational basis: its positive
    # normalized coefficients form a convex combination of pole differences.
    ratio=first[1][4]/last[1][4]
    if all(first[index][4]==ratio*last[index][4] for index in eachindex(first)) &&
            all(_model_nurbs_seam_near(first[index],last[index],tolerance_squared)
                for index in eachindex(first))
        return width
    end
    degree=axis==1 ? surface.degree_v : surface.degree_u
    # The native CAD surface constructors and BREP reader cap degree at25.
    degree<=25 || return 0.0
    knots=_NURBSSeamRational.(axis==1 ? surface.knots_v : surface.knots_u)
    cuts=sort!(unique!([_NURBSSeamRational(lo[3-axis]);
        knots[(_NURBSSeamRational(lo[3-axis]).<knots).&
              (knots.<_NURBSSeamRational(hi[3-axis]))];
        _NURBSSeamRational(hi[3-axis])]))
    budget=Ref(4096)
    for index in 1:length(cuts)-1
        lower,upper=cuts[index],cuts[index+1]
        fa=_model_nurbs_seam_bezier(first,degree,knots,lower,upper)
        la=_model_nurbs_seam_bezier(last,degree,knots,lower,upper)
        _model_nurbs_seam_panel(fa,la,tolerance_squared,budget) || return 0.0
    end
    return width
end
