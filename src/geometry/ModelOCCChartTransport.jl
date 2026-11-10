# A converted conic's B-spline parameter is affine on each rational
# quadratic span; it is not the source angle. These maps are used by
# boundary traces and analytic measure formulas, never by patch evaluation.
function _model_occ_conic_chart(lower,upper)
    count=max(1,ceil(Int,(upper-lower)/(2π/3)))
    width=(upper-lower)/count
    return _ModelOCCConicReparam(
        Float64[lower+i*width for i in 0:count],fill(_gm_tan(width/4),count))
end

function _occ_parameter_inverse_d2(mapping::_ModelOCCConicReparam,s)
    lower,width,middle,q=_occ_conic_reparam_lookup(mapping,Float64(s))
    z=(2*(s-lower)/width-1)*q
    first=4q/(width*(1+z*z))
    return middle+2*_gm_atan(z),first,-first*first*z
end

struct _ModelOCCScaledReparam{R}
    base::R
    scale::Float64
end
function _occ_conic_reparam(mapping::_ModelOCCScaledReparam,t)
    value,first,second=mapping.base===nothing ? (t/mapping.scale,1.,0.) :
        _occ_conic_reparam(mapping.base,t/mapping.scale)
    return value,first/mapping.scale,second/(mapping.scale*mapping.scale)
end
function _occ_parameter_inverse_d2(mapping::_ModelOCCScaledReparam,s)
    value,first,second=mapping.base===nothing ? (s,1.,0.) :
        _occ_parameter_inverse_d2(mapping.base,s)
    return value*mapping.scale,first*mapping.scale,second*mapping.scale
end
_occ_conic_reparam_inv(mapping::_ModelOCCScaledReparam,s)=
    first(_occ_parameter_inverse_d2(mapping,s))

# Restore the analytic coordinates for closed-form integration. The
# geometry is the same exact conic; only the integration chart changes.
function _occ_pcurve_d2(pc::NamedTuple{(:brep_analytic2d,)},t::Float64)
    base,maps=pc.brep_analytic2d
    point,tangent,acceleration=_occ_pcurve_d2(base,t)
    states=ntuple(2) do axis
        map=axis==1 ? maps.u : maps.v
        map===nothing && return (point[axis],tangent[axis],acceleration[axis])
        value,first,second=_occ_parameter_inverse_d2(map,point[axis])
        return value,first*tangent[axis],second*tangent[axis]^2+first*acceleration[axis]
    end
    return ntuple(i->states[i][1],2),ntuple(i->states[i][2],2),ntuple(i->states[i][3],2)
end
_occ_pcurve_d0(pc::NamedTuple{(:brep_analytic2d,)},t::Float64)=first(_occ_pcurve_d2(pc,t))
function _occ_pcurve_d1(pc::NamedTuple{(:brep_analytic2d,)},t::Float64)
    point,first,_=_occ_pcurve_d2(pc,t)
    return point,first
end

function _model_occ_analytic_surface_record(record)
    maps=get(record,:analytic_chart,nothing)
    maps===nothing && return record
    pcurves=Dict{Int,NamedTuple}()
    for (tag,entry) in get(record,:pcurves,Dict{Int,NamedTuple}())
        pcurves[tag]=(fwd=entry.fwd===nothing ? nothing : (brep_analytic2d=(entry.fwd,maps),),
            rev=entry.rev===nothing ? nothing : (brep_analytic2d=(entry.rev,maps),))
    end
    return merge(Base.structdiff(record,(patch=nothing,analytic_chart=nothing)),
        (pcurves=pcurves,))
end

function _model_occ_full_patch_integrals(record;kwargs...)
    lower,upper=_model_nurbs_surface_bounds(record)
    return _occ_natural_surface_integrals(merge(record,
        (natural_bounds=(Tuple(lower),Tuple(upper)),));kwargs...)
end
