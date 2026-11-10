# Positive line-measure weights permit a local weighted mean without forming
# tiny or overflowing world first moments. The rule's `first` field is this
# local mean; its absolute field is the matching mean absolute displacement.
function _model_occ_curve_integrals_rule(visit,panels,nodes,weights,origin,center)
    area=0.0;first=(0.,0.,0.);absolute_first=(0.,0.,0.)
    second=zeros(3,3);absolute_second=zeros(3,3)
    for (lower,upper) in panels
        half=upper/2-lower/2;middle=lower/2+upper/2
        for index in eachindex(nodes),side in (-1,1)
            parameter=middle+side*half*nodes[index]
            point,derivative=visit(parameter)
            weighted_length=(weights[index]*half)*hypot(derivative...)
            weighted_length==0 && continue
            shifted=_occ3_sub(point,origin)
            total=area+weighted_length
            fraction=weighted_length/total
            first=_occ3_add(first,_occ3_scale(_occ3_sub(shifted,first),fraction))
            absolute_first=_occ3_add(absolute_first,
                _occ3_scale(_occ3_sub(abs.(shifted),absolute_first),fraction))
            area=total
            if center!==nothing
                delta=_occ3_sub(shifted,center)
                for a in 1:3,b in 1:3
                    term=(weighted_length*delta[a])*delta[b]
                    second[a,b]+=term;absolute_second[a,b]+=abs(term)
                end
            end
        end
    end
    return (;area,first,second,
        absolute=(area=area,first=absolute_first,second=absolute_second))
end

function _model_occ_curve_integrals(visit,panels,origin,order;center=nothing,
                                    fields=(:area,:first))
    length(panels)<=65536 || throw(ArgumentError("OCC curve: quadrature capacity exceeded"))
    work=length(panels)*order
    work<=2097152 || throw(ArgumentError("OCC curve: quadrature capacity exceeded"))
    nodes,weights=_occ_natural_gauss_rule(order)
    previous=_model_occ_curve_integrals_rule(visit,panels,nodes,weights,origin,center)
    order*=2
    nodes,weights=_occ_natural_gauss_rule(order)
    for refinement in 0:8
        length(panels)<=div(2097152-work,order) || throw(ArgumentError(
            "OCC curve: adaptive quadrature capacity exceeded"))
        work+=length(panels)*order
        current=_model_occ_curve_integrals_rule(visit,panels,nodes,weights,origin,center)
        _occ_natural_integrals_converged(previous,current,fields) && return current
        refinement<8 || break
        previous=current
        refined=Tuple{Float64,Float64}[]
        for (lower,upper) in panels
            middle=_convex_coordinate(lower,upper,.5)
            lower<middle<upper || throw(ArgumentError(
                "OCC curve: no representable quadrature subdivision"))
            push!(refined,(lower,middle),(middle,upper))
        end
        length(refined)<=65536 || throw(ArgumentError("OCC curve: quadrature capacity exceeded"))
        panels=refined
    end
    throw(ArgumentError("OCC curve: quadrature did not converge"))
end

function _model_occ_curve_stats(rec;want_second=false)
    visit=rec.occ===:ellipse ?
        t->(_occ_ellipse_point(rec,t),_occ_ellipse_derivative(rec,t)) :
        rec.occ===:copied_curve ?
        t->begin values=_model_occ_curve_d2(rec.curve,t)
            (values[1],values[2]) end :
        rec.occ===:surface_iso ?
        t->(_model_surface_iso_point(rec,t),_model_surface_iso_d1(rec,t)) :
        throw(ArgumentError("OCC curve: unsupported measure kind $(rec.occ)"))
    bounds=Float64[rec.t0]
    rec.occ===:copied_curve && append!(bounds,
        _model_occ_curve_panels(rec.curve,rec.t0,rec.t1))
    push!(bounds,rec.t1)
    panels=Tuple{Float64,Float64}[(bounds[index],bounds[index+1])
        for index in 1:length(bounds)-1]
    origin=visit(rec.t0/2+rec.t1/2)[1]
    degree=rec.occ===:copied_curve ? rec.curve.curve.degree : 2
    order=16
    while order<cld(5degree,2);order*=2;end
    order<=64 || throw(ArgumentError("OCC curve: quadrature degree exceeds capacity"))
    integrals=_model_occ_curve_integrals(visit,panels,origin,order)
    centroid=integrals.first
    second=want_second ? _model_occ_curve_integrals(visit,panels,origin,order;
        center=centroid,fields=(:second,)).second : zeros(3,3)
    return (area=integrals.area,origin=origin,centroid=centroid,second=second)
end
