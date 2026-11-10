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

const _OCC_NATURAL_GAUSS_RULES=Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const _OCC_NATURAL_GAUSS_LOCK=ReentrantLock()

function _occ_natural_gauss_rule(order)
    order in (8,16,32,64,128) || throw(ArgumentError(
        "OCC natural face: unsupported quadrature order"))
    order==8 && return (_OCC_GL_T,_OCC_GL_W)
    lock(_OCC_NATURAL_GAUSS_LOCK) do
        get!(_OCC_NATURAL_GAUSS_RULES,order) do
            # Legendre's orthonormal three-term recurrence gives this
            # symmetric Jacobi matrix (Golub-Welsch quadrature).
            # https://dlmf.nist.gov/3.5.vi
            offdiagonal=Float64[index/sqrt(4index^2-1) for index in 1:order-1]
            decomposition=LinearAlgebra.eigen(LinearAlgebra.SymTridiagonal(
                zeros(order),offdiagonal))
            nodes=Float64[(decomposition.values[order+1-index]-
                decomposition.values[index])/2 for index in 1:order÷2]
            weights=Float64[decomposition.vectors[1,index]^2+
                decomposition.vectors[1,order+1-index]^2 for index in 1:order÷2]
            weights./=sum(weights)
            return nodes,weights
        end
    end
end

function _occ_natural_measure_panels(g,axis)
    lo,hi=g.natural_bounds
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
    cuts=Float64[_convex_coordinate(lower,upper,i/count) for i in 0:count]
    return collect(zip(cuts[1:end-1],cuts[2:end]))
end

function _occ_natural_surface_integrals_rule(g,upanels,vpanels,nodes,weights;
        centered=false,center=nothing)
    length(upanels)<=div(65536,max(1,length(vpanels))) || throw(ArgumentError(
        "OCC natural face: quadrature capacity exceeded"))
    area=0.0;first=(0.0,0.0,0.0);second=zeros(3,3)
    volume=0.0;div_first=(0.0,0.0,0.0);div_second=zeros(3,3)
    absolute_first=(0.0,0.0,0.0);absolute_second=zeros(3,3)
    absolute_volume=0.0;absolute_div_first=(0.0,0.0,0.0)
    absolute_div_second=zeros(3,3)
    for (ua,ub) in upanels,(va,vb) in vpanels
        hu=ub/2-ua/2;mu=ua/2+ub/2;hv=vb/2-va/2;mv=va/2+vb/2
        for i in eachindex(nodes),si in (-1,1),j in eachindex(nodes),sj in (-1,1)
            u=mu+si*hu*nodes[i];v=mv+sj*hv*nodes[j]
            p,du,dv=_occ_surface_d1(g,u,v)
            normal=_occ_cross(du,dv)
            weight=weights[i]*weights[j]*hu*hv
            wa=weight*_occ_modulus(normal);rn=weight*_dot3(p,normal)
            area+=wa;volume+=rn/3
            first=_occ3_add(first,_occ3_scale(p,wa))
            div_first=_occ3_add(div_first,_occ3_scale(p,rn))
            absolute_volume+=abs(rn/3)
            absolute_first=_occ3_add(absolute_first,abs.(_occ3_scale(p,wa)))
            absolute_div_first=_occ3_add(absolute_div_first,abs.(_occ3_scale(p,rn)))
            for a in 1:3,b in 1:3
                surface_term=wa*p[a]*p[b];volume_term=rn*p[a]*p[b]/5
                second[a,b]+=surface_term;div_second[a,b]+=volume_term
                absolute_second[a,b]+=abs(surface_term)
                absolute_div_second[a,b]+=abs(volume_term)
            end
        end
    end
    # Centered second moment and divergence tensor ∬(r-c)⊗(r-c),
    # ∬(r-c)⊗(r-c)((r-c)·n)dA, evaluated directly in a second pass so faces
    # far from the origin do not lose accuracy to subtraction.
    centered_second=zeros(3,3)
    centered_div_second=zeros(3,3)
    absolute_centered_second=zeros(3,3)
    absolute_centered_div_second=zeros(3,3)
    if centered && area>0
        c=center===nothing ? _occ3_scale(first,1/area) : center
        for (ua,ub) in upanels,(va,vb) in vpanels
            hu=ub/2-ua/2;mu=ua/2+ub/2;hv=vb/2-va/2;mv=va/2+vb/2
            for i in eachindex(nodes),si in (-1,1),j in eachindex(nodes),sj in (-1,1)
                u=mu+si*hu*nodes[i];v=mv+sj*hv*nodes[j]
                p,du,dv=_occ_surface_d1(g,u,v)
                weight=weights[i]*weights[j]*hu*hv
                normal=_occ_cross(du,dv)
                d=_occ3_sub(p,c);wa=weight*_occ_modulus(normal)
                rn=weight*_dot3(d,normal)
                for a in 1:3,b in 1:3
                    surface_term=wa*d[a]*d[b];volume_term=rn*d[a]*d[b]/5
                    centered_second[a,b]+=surface_term
                    centered_div_second[a,b]+=volume_term
                    absolute_centered_second[a,b]+=abs(surface_term)
                    absolute_centered_div_second[a,b]+=abs(volume_term)
                end
            end
        end
    end
    absolute=(area=area,first=absolute_first,second=absolute_second,
        centered_second=absolute_centered_second,
        centered_div_second=absolute_centered_div_second,
        volume=absolute_volume,div_first=absolute_div_first,
        div_second=absolute_div_second)
    return (;area,first,second,centered_second,centered_div_second,
        volume,div_first,div_second,absolute)
end

function _occ_natural_integrals_converged(previous,current,fields)
    for field in (:area,:first,:second,:volume,:div_first,:div_second,
                  :centered_second,:centered_div_second)
        fields===nothing || field in fields || continue
        a=getproperty(previous,field);b=getproperty(current,field)
        bound_a=getproperty(previous.absolute,field)
        bound_b=getproperty(current.absolute,field)
        for (x,y,ba,bb) in zip(a isa Number ? (a,) : a,b isa Number ? (b,) : b,
                bound_a isa Number ? (bound_a,) : bound_a,
                bound_b isa Number ? (bound_b,) : bound_b)
            # Equal signed infinities represent an out-of-range result;
            # NaNs and inconsistent rules never establish convergence.
            isinf(x) && x==y && continue
            isfinite(x) && isfinite(y) || return false
            isfinite(ba) && isfinite(bb) || return false
            # Absolute quadrature contributions supply a component's own
            # cancellation scale. A diameter power can overflow or swamp the
            # finite volume of a long, thin solid.
            tolerance=1e-12*max(ba,bb)+1e-11*max(abs(x),abs(y))
            abs(x-y)<=tolerance || return false
        end
    end
    return true
end

function _occ_natural_surface_integrals(g;centered=false,center=nothing,fields=nothing)
    upanels=_occ_natural_measure_panels(g,1);vpanels=_occ_natural_measure_panels(g,2)
    if !_model_nurbs_surface(g)
        return _occ_natural_surface_integrals_rule(g,upanels,vpanels,
            _OCC_GL_T,_OCC_GL_W;centered,center)
    end
    surface=g.patch.surface
    # A polynomial tensor's divergence second moments have degree at most
    # 5p-1 on each knot span. Rational patches additionally compare increasing
    # orders and refine their spans instead of trusting a fixed GL8 rule.
    required=cld(5max(surface.degree_u,surface.degree_v),2)
    rational=any(!=(first(surface.weights)),surface.weights)
    order=8
    while order<max(required,rational ? 16 : 8);order*=2;end
    pairs=length(upanels)*length(vpanels)
    pairs<=div(33554432,order^2) || throw(ArgumentError(
        "OCC natural face: adaptive quadrature capacity exceeded"))
    work=pairs*order^2
    nodes,weights=_occ_natural_gauss_rule(order)
    previous=_occ_natural_surface_integrals_rule(g,upanels,vpanels,nodes,weights;centered,center)
    order*=2
    nodes,weights=_occ_natural_gauss_rule(order)
    for refinement in 0:8
        pairs=length(upanels)*length(vpanels)
        pairs<=div(33554432-work,order^2) || throw(ArgumentError(
            "OCC natural face: adaptive quadrature capacity exceeded"))
        work+=pairs*order^2
        current=_occ_natural_surface_integrals_rule(g,upanels,vpanels,nodes,weights;centered,center)
        _occ_natural_integrals_converged(previous,current,fields) && return current
        refinement<8 || break
        previous=current
        for panels in (upanels,vpanels)
            refined=Tuple{Float64,Float64}[]
            for (lower,upper) in panels
                middle=_convex_coordinate(lower,upper,.5)
                lower<middle<upper || throw(ArgumentError(
                    "OCC natural face: no representable quadrature subdivision"))
                push!(refined,(lower,middle),(middle,upper))
            end
            empty!(panels);append!(panels,refined)
        end
    end
    throw(ArgumentError("OCC natural face: quadrature did not converge"))
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
        for axis in 1:2,upper in (false,true)
            value=upper ? hi[axis] : lo[axis]
            point=_model_nurbs_collapsed_iso_point(g,axis,value)
            point===nothing || push!(sides,(axis,value,point))
        end
    end
    return sides
end

function _occ_natural_surface_point(g,u,v,up=_brep_uperiod(g),vp=_brep_vperiod(g))
    lo,hi=g.natural_bounds
    for (axis,value,point) in g.natural_collapsed
        (axis==1 ? u : v)==value && return point
    end
    canonical_u=up>0 && hi[1]-lo[1]>=up-64eps(up) && u==hi[1] ? lo[1] : u
    canonical_v=vp>0 && hi[2]-lo[2]>=vp-64eps(vp) && v==hi[2] ? lo[2] : v
    if g.occ===:sphere && abs(canonical_v)==π/2
        return ntuple(i->g.center[i]+copysign(g.radius,canonical_v)*g.axis[i],3)
    end
    return _occ_surface_point(g,Float64(canonical_u),Float64(canonical_v))
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
    function append_edge(visit,a,b,depth)
        if score(a,b)>.6
            depth<20 || throw(ArgumentError(
                "$caller: natural Surface[$t] boundary exceeds subdivision capacity"))
            mid=_occ_natural_uv_lerp(a,b,.5)
            (mid==a || mid==b) && throw(ArgumentError(
                "$caller: natural Surface[$t] boundary has no representable subdivision"))
            visit(visit,a,mid,depth+1);visit(visit,mid,b,depth+1)
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
            append_edge(append_edge,p,q,0)
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

# A natural rectangle has only Surface-owned chart cuts. Certify its
# physical quotient before adaptation: repeated pole preimages disappear,
# while every surviving triangle must be nondegenerate and embedded.
function _occ_natural_cover_certificate(cells,world,caller)
    incidence=_metric_surface_incidence(cells,caller)
    boundary=Set(edge for (edge,hits) in incidence if length(hits)==1)
    _metric_surface_quotient_ok(cells,world,boundary) || return false
    lookup=Dict{NTuple{3,Float64},Int32}()
    rows=Int32[get!(lookup,p,Int32(length(lookup)+1)) for p in world]
    triangles=NTuple{3,Int32}[]
    for cell in cells
        a,b,c=(rows[row] for row in cell)
        (a==b || b==c || c==a) && continue
        triangle_area((world[row] for row in cell)...)>0 || return false
        push!(triangles,(a,b,c))
    end
    points=Matrix{Float64}(undef,3,length(lookup))
    for (point,row) in lookup;points[:,row].=point;end
    connectivity=Matrix{Int32}(undef,3,length(triangles))
    for (column,cell) in enumerate(triangles);connectivity[:,column].=cell;end
    return _plc_self_intersection3(Mesh(points;tris=connectivity))===nothing
end

function _occ_natural_grid_stations(g,axis,values)
    lo,hi=g.natural_bounds
    stations=Float64[values;[_convex_coordinate(lo[axis],hi[axis],i/4) for i in 0:4]]
    if _model_nurbs_surface(g)
        surface=g.patch.surface
        knots=axis==1 ? surface.knots_u : surface.knots_v
        reparam=get(g,:reparam,nothing)
        mapping=reparam===nothing ? nothing : axis==1 ? reparam.u : reparam.v
        bounds=_model_nurbs_surface_bounds(g)
        for knot in knots
            bounds[1][axis]<=knot<=bounds[2][axis] || continue
            value=mapping===nothing ? knot : _occ_conic_reparam_inv(mapping,knot)
            lo[axis]<value<hi[axis] && push!(stations,value)
        end
    end
    return sort!(unique!(stations))
end

# If a sparse CDT shortcuts a periodic or collapsed chart, rebuild genuine
# parameter rings. Aligned rings make polar fans and both seam partitions
# conform before the size-field producer performs its local adaptation.
# Existing embedded covers are retained verbatim.
function _occ_natural_surface_cover!(g,u,v,world,cells,owners,parameters,tag,caller)
    _occ_natural_cover_certificate(cells,world,caller) && return cells
    lo,hi=g.natural_bounds
    us=_occ_natural_grid_stations(g,1,
        Float64[u[row] for row in eachindex(u) if v[row]==lo[2] || v[row]==hi[2]])
    vs=_occ_natural_grid_stations(g,2,
        Float64[v[row] for row in eachindex(v) if u[row]==lo[1] || u[row]==hi[1]])
    up=_brep_uperiod(g);vp=_brep_vperiod(g)
    capacity=max(65536,length(world))
    for attempt in 0:16
        length(us)<=div(capacity,length(vs)) || break
        empty!(u);empty!(v);empty!(world);empty!(owners);empty!(parameters)
        for second in vs,first in us
            point=_occ_natural_surface_point(g,first,second,up,vp)
            all(isfinite,point) || throw(ArgumentError(
                "$caller: natural Surface chart exceeds Float64"))
            push!(u,first);push!(v,second);push!(world,point)
            push!(owners,(2,Int32(tag)));push!(parameters,[first,second])
        end
        cells=NTuple{3,Int32}[]
        count=length(us)
        for j in 1:length(vs)-1,i in 1:count-1
            a=Int32((j-1)*count+i);b=a+Int32(1);c=a+Int32(count);d=c+Int32(1)
            push!(cells,(a,b,d),(a,d,c))
        end
        _occ_natural_cover_certificate(cells,world,caller) && return cells
        attempt<16 || break
        for stations in (us,vs)
            midpoints=Float64[_convex_coordinate(stations[index],stations[index+1],.5)
                for index in 1:length(stations)-1]
            append!(stations,midpoints)
            sort!(unique!(stations))
        end
    end
    throw(ArgumentError("$caller: natural Surface cover exceeds quotient subdivision capacity"))
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
