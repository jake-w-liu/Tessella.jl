# A face can use its surface's finite natural domain without any wires.
# The chart perimeter belongs to that Surface; it does not introduce CAD
# Point or Curve entities, including along a periodic seam or a pole row.
function _occ_brep_natural_curve_bounds(curve,caller)
    curve.kind==8 && return (curve.lo,curve.hi)
    curve.kind==9 && return _occ_brep_natural_curve_bounds(curve.base,caller)
    curve.kind in (2,3) && return (0.0,_OCC_TWO_PI)
    if curve.kind in (6,7)
        native=_occ_brep_nurbs_curve(curve,caller)
        return (native.knots[native.degree+1],native.knots[end-native.degree])
    end
    throw(ArgumentError("$caller: BREP surface has an unbounded generatrix"))
end

function _occ_brep_natural_bounds(surface,caller)
    if surface.kind==10
        a,b,c,d=surface.bounds
        return ((a,c),(b,d))
    elseif surface.kind==11
        return _occ_brep_natural_bounds(surface.base,caller)
    elseif surface.kind==4
        return ((0.0,-π/2),(_OCC_TWO_PI,π/2))
    elseif surface.kind==5
        return ((0.0,0.0),(_OCC_TWO_PI,_OCC_TWO_PI))
    elseif surface.kind==7
        lower,upper=_occ_brep_natural_curve_bounds(surface.curve,caller)
        return ((0.0,lower),(_OCC_TWO_PI,upper))
    elseif surface.kind in (8,9)
        g,_,_=_occ_brep_surface_geometry(surface,((0.0,0.0),(1.0,1.0)),caller)
        lo,hi=_model_nurbs_surface_bounds(g)
        return (Tuple(lo),Tuple(hi))
    end
    throw(ArgumentError("$caller: BREP face has an unbounded natural surface domain"))
end

function _occ_natural_surface_integrals(g;centered=false,center=nothing)
    lo,hi=g.natural_bounds
    function panels(axis)
        lower,upper=lo[axis],hi[axis]
        if _model_nurbs_surface(g)
            knots=axis==1 ? g.patch.surface.knots_u : g.patch.surface.knots_v
            return _model_occ_span_cuts(knots,lower,upper)
        end
        span=upper-lower
        isfinite(span) && span>0 || throw(ArgumentError(
            "OCC natural face: domain width exceeds Float64"))
        count=span<=π/2 ? 1 : ceil(Int,min(span/(π/2),4097.0))
        count<=4096 || throw(ArgumentError("OCC natural face: quadrature capacity exceeded"))
        cuts=Float64[(1-i/count)*lower+(i/count)*upper for i in 0:count]
        return collect(zip(cuts[1:end-1],cuts[2:end]))
    end
    upanels=panels(1);vpanels=panels(2)
    length(upanels)<=div(65536,max(1,length(vpanels))) || throw(ArgumentError(
        "OCC natural face: quadrature capacity exceeded"))
    area=0.0;first=(0.0,0.0,0.0);second=zeros(3,3)
    volume=0.0;div_first=(0.0,0.0,0.0);div_second=zeros(3,3)
    for (ua,ub) in upanels,(va,vb) in vpanels
        hu=ub/2-ua/2;mu=ua/2+ub/2;hv=vb/2-va/2;mv=va/2+vb/2
        for i in 1:4,si in (-1,1),j in 1:4,sj in (-1,1)
            u=mu+si*hu*_OCC_GL_T[i];v=mv+sj*hv*_OCC_GL_T[j]
            p,du,dv=_occ_surface_d1(g,u,v)
            normal=_occ_cross(du,dv)
            weight=_OCC_GL_W[i]*_OCC_GL_W[j]*hu*hv
            wa=weight*_occ_modulus(normal);rn=weight*_dot3(p,normal)
            area+=wa;volume+=rn/3
            first=_occ3_add(first,_occ3_scale(p,wa))
            div_first=_occ3_add(div_first,_occ3_scale(p,rn))
            for a in 1:3,b in 1:3
                second[a,b]+=wa*p[a]*p[b]
                div_second[a,b]+=rn*p[a]*p[b]/5
            end
        end
    end
    # Centered second moment and divergence tensor ∬(r-c)⊗(r-c),
    # ∬(r-c)⊗(r-c)((r-c)·n)dA, evaluated directly in a second pass so faces
    # far from the origin do not lose accuracy to subtraction.
    centered_second=zeros(3,3)
    centered_div_second=zeros(3,3)
    if centered && area>0
        c=center===nothing ? _occ3_scale(first,1/area) : center
        for (ua,ub) in upanels,(va,vb) in vpanels
            hu=ub/2-ua/2;mu=ua/2+ub/2;hv=vb/2-va/2;mv=va/2+vb/2
            for i in 1:4,si in (-1,1),j in 1:4,sj in (-1,1)
                u=mu+si*hu*_OCC_GL_T[i];v=mv+sj*hv*_OCC_GL_T[j]
                p,du,dv=_occ_surface_d1(g,u,v)
                weight=_OCC_GL_W[i]*_OCC_GL_W[j]*hu*hv
                normal=_occ_cross(du,dv)
                d=_occ3_sub(p,c);wa=weight*_occ_modulus(normal)
                rn=weight*_dot3(d,normal)
                for a in 1:3,b in 1:3
                    centered_second[a,b]+=wa*d[a]*d[b]
                    centered_div_second[a,b]+=rn*d[a]*d[b]/5
                end
            end
        end
    end
    return (;area,first,second,centered_second,centered_div_second,
        volume,div_first,div_second)
end

function _occ_natural_surface_bbox(g)
    if _model_nurbs_surface(g)
        return _model_occ_controls_box(g.patch.surface.controls)
    elseif g.occ===:sphere
        return (g.center[1]-g.radius,g.center[2]-g.radius,g.center[3]-g.radius,
                g.center[1]+g.radius,g.center[2]+g.radius,g.center[3]+g.radius)
    elseif g.occ===:torus
        return _occ_torus_bounding_box(g)
    end
    lo,hi=g.natural_bounds
    angles=Float64[lo[1],hi[1]]
    if g.occ in (:cylinder,:cone)
        Y=_occ_frame_y(g)
        for axis in 1:3,opposite in (0.0,π)
            angle=lo[1]+mod(_gm_atan2(Y[axis],g.X[axis])+opposite-lo[1],_OCC_TWO_PI)
            angle<=hi[1] && push!(angles,angle)
        end
    end
    return _model_occ_box_points(_occ_surface_point(g,u,v)
        for u in angles for v in (lo[2],hi[2]))
end

function _occ_natural_collapsed_sides(g)
    lo,hi=g.natural_bounds
    sides=Tuple{Int,Float64,NTuple{3,Float64}}[]
    if g.occ===:sphere
        for v in (lo[2],hi[2])
            abs(v)==π/2 || continue
            point=ntuple(i->g.center[i]+copysign(g.radius,v)*g.axis[i],3)
            push!(sides,(2,v,point))
        end
    elseif _model_nurbs_surface(g)
        surface=g.patch.surface
        for axis in 1:2,upper in (false,true)
            knots=axis==1 ? surface.knots_u : surface.knots_v
            degree=axis==1 ? surface.degree_u : surface.degree_v
            value=upper ? hi[axis] : lo[axis]
            terminal=upper ? knots[end-degree:end] : knots[1:degree+1]
            all(==(value),terminal) || continue
            index=upper ? size(surface.controls,axis) : 1
            poles=axis==1 ? @view(surface.controls[index,:]) :
                           @view(surface.controls[:,index])
            point=first(poles)
            all(==(point),poles) && push!(sides,(axis,value,point))
        end
    end
    return sides
end

function _occ_natural_surface_point(g,u,v,up=_brep_uperiod(g),vp=_brep_vperiod(g))
    lo,hi=g.natural_bounds
    for (axis,value,point) in g.natural_collapsed
        (axis==1 ? u : v)==value && return point
    end
    up>0 && hi[1]-lo[1]>=up-64eps(up) && u==hi[1] && (u=lo[1])
    vp>0 && hi[2]-lo[2]>=vp-64eps(vp) && v==hi[2] && (v=lo[2])
    if g.occ===:sphere && abs(v)==π/2
        return ntuple(i->g.center[i]+copysign(g.radius,v)*g.axis[i],3)
    end
    return _occ_surface_point(g,Float64(u),Float64(v))
end

function _occ_natural_uv_lerp(a,b,fraction)
    fraction==0 && return a
    fraction==1 && return b
    return ntuple(k->a[k]==b[k] ? a[k] :
        (1-fraction)*a[k]+fraction*b[k],2)
end

function _occ_natural_surface_scaffold!(m,t,g,emit,segs,options,caller)
    lo,hi=g.natural_bounds
    all(isfinite,(lo...,hi...)) && lo[1]<hi[1] && lo[2]<hi[2] ||
        throw(ArgumentError("$caller: Surface[$t] has no finite natural rectangle"))
    options===nothing && (options=_model_default_mesh1d_options(m,caller))
    background=_ModelAlgorithmBackgroundSize(m,2,t,options,nothing,caller)
    field=options.inner_field
    up=_brep_uperiod(g);vp=_brep_vperiod(g)
    lift=(u,v)->_occ_natural_surface_point(g,u,v,up,vp)
    function score(a,b)
        mid=_occ_natural_uv_lerp(a,b,.5)
        pa=lift(a...);pm=lift(mid...);pb=lift(b...)
        if field isa AbstractAnisoField
            return max(metric_edge_length(field,pa,pb;entity=(2,t)),
                metric_edge_length(field,pa,pm;entity=(2,t))+
                metric_edge_length(field,pm,pb;entity=(2,t)))
        end
        size=min(background(pa...),background(pm...),background(pb...))
        return (_b3dist(pa,pm)+_b3dist(pm,pb))/size
    end
    # Opposite periodic sides evaluate the same canonical world points, so
    # they receive the same subdivision and close with literal world rows.
    corners=((lo[1],lo[2]),(hi[1],lo[2]),(hi[1],hi[2]),(lo[1],hi[2]))
    count=Ref(0)
    function append_edge(a,b,depth)
        if score(a,b)>.6
            depth<20 || throw(ArgumentError(
                "$caller: natural Surface[$t] boundary exceeds subdivision capacity"))
            mid=_occ_natural_uv_lerp(a,b,.5)
            (mid==a || mid==b) && throw(ArgumentError(
                "$caller: natural Surface[$t] boundary has no representable subdivision"))
            append_edge(a,mid,depth+1);append_edge(mid,b,depth+1)
            return nothing
        end
        count[]<65536 || throw(ArgumentError(
            "$caller: natural Surface[$t] boundary exceeds 65536 segments"))
        count[]+=1
        i=emit(a,lift(a...),(2,Int32(t)),[a...])
        j=emit(b,lift(b...),(2,Int32(t)),[b...])
        push!(segs,(i,j))
        return nothing
    end
    for side in 1:4
        a=corners[side];b=corners[mod1(side+1,4)]
        for station in 0:3
            p=_occ_natural_uv_lerp(a,b,station/4)
            q=_occ_natural_uv_lerp(a,b,(station+1)/4)
            append_edge(p,q,0)
        end
    end
    return nothing
end

function _occ_natural_collapsed_groups(g,parameters,coordinates)
    groups=Dict{Int32,Int}()
    points=Dict{NTuple{3,Float64},Int}()
    for (axis,value,point) in g.natural_collapsed
        group=get!(points,point,length(points)+1)
        for row in eachindex(parameters)
            parameters[row][axis]==value && isequal(coordinates[row],point) || continue
            groups[Int32(row)]=group
        end
    end
    return groups
end

function _occ_natural_surface_closed(g)
    hasproperty(g,:natural_bounds) || return false
    lo,hi=g.natural_bounds
    periods=(_brep_uperiod(g),_brep_vperiod(g))
    for axis in 1:2
        period=periods[axis]
        period>0 && isapprox(hi[axis]-lo[axis],period;rtol=1e-12,atol=0) && continue
        all(any(side->side[1]==axis && side[2]==value,g.natural_collapsed)
            for value in (lo[axis],hi[axis])) || return false
    end
    return true
end

# Derive the expected chart boundary from the declared natural rectangle,
# independently of the emitted triangle incidence. Its periodic/degenerate
# images are checked by the existing final quotient certificate.
function _occ_natural_surface_boundary(g,parameters,numbers,caller)
    lo,hi=g.natural_bounds
    wanted=Set{NTuple{2,UInt64}}()
    for axis in 1:2,side in (lo[axis],hi[axis])
        rows=findall(p->p[axis]==side,parameters)
        sort!(rows;by=row->parameters[row][3-axis])
        length(rows)>=2 || throw(ArgumentError(
            "$caller: natural Surface boundary lacks its domain corners"))
        parameters[first(rows)][3-axis]==lo[3-axis] &&
            parameters[last(rows)][3-axis]==hi[3-axis] || throw(ArgumentError(
                "$caller: natural Surface boundary misses its declared domain"))
        for index in 1:length(rows)-1
            push!(wanted,minmax(numbers[rows[index]],numbers[rows[index+1]]))
        end
    end
    return wanted
end
