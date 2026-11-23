# An OCC-style tensor-product patch owns its control geometry and the separate
# BRep boundary entities. Control Point identity is not boundary-vertex identity.
struct _ModelNURBSSurface
    family::Symbol
    surface::NURBS.NURBSSurface
    nets::NTuple{6,Union{Nothing,Matrix{NTuple{4,Float64}}}}
end

_model_nurbs_surface(g)=g!==nothing && hasproperty(g,:patch) &&
    g.patch isa _ModelNURBSSurface

function _model_nurbs_derivative_net(net,degree,knots,axis)
    degree==0 && return nothing
    nu,nv=size(net)
    out=Matrix{NTuple{4,Float64}}(undef,nu-(axis==1),nv-(axis==2))
    for j in axes(out,2),i in axes(out,1)
        k=axis==1 ? i : j
        denominator=knots[k+degree+1]-knots[k+1]
        factor=denominator==0.0 ? 0.0 : degree/denominator
        a=net[i,j];b=axis==1 ? net[i+1,j] : net[i,j+1]
        out[i,j]=ntuple(component->factor*(b[component]-a[component]),4)
    end
    return out
end

function _ModelNURBSSurface(family,surface::NURBS.NURBSSurface)
    net=map(surface.controls,surface.weights) do point,weight
        (weight*point[1],weight*point[2],weight*point[3],weight)
    end
    all(point->all(isfinite,point),net) || throw(ArgumentError(
        "NURBS Surface: homogeneous control points must be finite"))
    pu,pv=surface.degree_u,surface.degree_v
    U,V=surface.knots_u,surface.knots_v
    du=_model_nurbs_derivative_net(net,pu,U,1)
    dv=_model_nurbs_derivative_net(net,pv,V,2)
    duu=pu<2 ? nothing : _model_nurbs_derivative_net(du,pu-1,U[2:end-1],1)
    dvv=pv<2 ? nothing : _model_nurbs_derivative_net(dv,pv-1,V[2:end-1],2)
    duv=_model_nurbs_derivative_net(du,pv,V,2)
    return _ModelNURBSSurface(family,surface,(net,du,dv,duu,dvv,duv))
end

# De Boor's algorithm on homogeneous points, in a buffer the caller owns
# (`degree+1` entries are used). `knots` may be a view: the derivative nets
# use the knot vector with its ends trimmed.
function _model_nurbs_deboor_h!(buffer::Vector{NTuple{4,Float64}},points,degree::Int,
                                knots::AbstractVector{Float64},u::Float64)
    count=length(points)
    span=clamp(searchsortedlast(knots,u)-1,degree,count-1)
    # A periodic spline's active upper endpoint can precede repeated knots
    # in its cyclic extension. Evaluate its last nonzero span, including a
    # writer-rounded edge endpoint just beyond the active knot domain.
    while span>degree && knots[span+1]==knots[span+2]
        span-=1
    end
    @inbounds for j in 0:degree
        buffer[j+1]=points[span-degree+j+1]
    end
    @inbounds for level in 1:degree,j in degree:-1:level
        lower=knots[span-degree+j+1]
        upper=knots[span+j-level+2]
        factor=(u-lower)/(upper-lower)
        left=buffer[j];right=buffer[j+1]
        buffer[j+1]=((1.0-factor)*left[1]+factor*right[1],
                     (1.0-factor)*left[2]+factor*right[2],
                     (1.0-factor)*left[3]+factor*right[3],
                     (1.0-factor)*left[4]+factor*right[4])
    end
    return @inbounds buffer[degree+1]
end

function _model_nurbs_deboor_h(points,degree,knots,u)
    return _model_nurbs_deboor_h!(Vector{NTuple{4,Float64}}(undef,degree+1),
        points,Int(degree),knots,Float64(u))
end

# Scratch for one evaluation of a patch: the column of row results and one
# de Boor buffer.
function _model_nurbs_scratch(patch::_ModelNURBSSurface)
    surface=patch.surface
    net=patch.nets[1]::Matrix{NTuple{4,Float64}}
    return (Vector{NTuple{4,Float64}}(undef,size(net,1)),
            Vector{NTuple{4,Float64}}(undef,max(surface.degree_u,surface.degree_v)+1))
end

function _model_nurbs_homogeneous(patch::_ModelNURBSSurface,slot::Int,u::Float64,
        v::Float64,columns::Vector{NTuple{4,Float64}},
        buffer::Vector{NTuple{4,Float64}})
    net=patch.nets[slot]
    net===nothing && return (0.0,0.0,0.0,0.0)
    ou,ov=slot==1 ? (0,0) : slot==2 ? (1,0) : slot==3 ? (0,1) :
          slot==4 ? (2,0) : slot==5 ? (0,2) : (1,1)
    surface=patch.surface
    U=@view surface.knots_u[1+ou:end-ou]
    V=@view surface.knots_v[1+ov:end-ov]
    rows=size(net,1)
    for i in 1:rows
        columns[i]=_model_nurbs_deboor_h!(buffer,@view(net[i,:]),
            surface.degree_v-ov,V,v)
    end
    return _model_nurbs_deboor_h!(buffer,@view(columns[1:rows]),
        surface.degree_u-ou,U,u)
end

function _model_nurbs_homogeneous(patch,slot,u,v)
    columns,buffer=_model_nurbs_scratch(patch)
    return _model_nurbs_homogeneous(patch,Int(slot),Float64(u),Float64(v),
        columns,buffer)
end

@inline function _model_nurbs_weight_check(h,u,v)
    isfinite(h[4]) && !iszero(h[4]) || throw(ArgumentError(
        "NURBS Surface: zero or non-finite homogeneous weight at ($u,$v)"))
    return nothing
end

# A patch that closes on itself in one direction (clamped there, first and
# last control rows equal, weights included) is periodic in it: a parameter
# one period outside the knot domain names the surface point one period back.
# Evaluating it as given instead extends the last span's polynomial, which
# for a converted torus stays on the torus but at another point and with
# other derivatives. A face meshed on a sheet shifted by a period (the sheet
# is anchored where its seam curve's pcurve lies) then folded over itself.
# Excursions within rounding of the domain keep their one-sided evaluation.
@inline function _model_nurbs_closed_parameters(patch::_ModelNURBSSurface,u::Float64,v::Float64)
    s=patch.surface
    lo=s.knots_u[s.degree_u+1];hi=s.knots_u[end-s.degree_u]
    tolerance=1e-9*(hi-lo)
    (u<lo-tolerance || u>hi+tolerance) && (u=_model_nurbs_wrapped(s,1,u,lo,hi))
    lo=s.knots_v[s.degree_v+1];hi=s.knots_v[end-s.degree_v]
    tolerance=1e-9*(hi-lo)
    (v<lo-tolerance || v>hi+tolerance) && (v=_model_nurbs_wrapped(s,2,v,lo,hi))
    return u,v
end

@noinline function _model_nurbs_wrapped(s,axis::Int,value::Float64,lo::Float64,hi::Float64)
    width=hi-lo
    isfinite(value) && width>0 && _model_nurbs_closed(s,axis) || return value
    return lo+mod(value-lo,width)
end

function _model_nurbs_closed(s,axis::Int)
    knots=axis==1 ? s.knots_u : s.knots_v
    degree=axis==1 ? s.degree_u : s.degree_v
    for k in 1:degree
        knots[k+1]==knots[1] && knots[end-k]==knots[end] || return false
    end
    controls=s.controls;weights=s.weights
    last=size(controls,axis)
    last>1 || return false
    scale=0.0
    for point in controls
        scale=max(scale,abs(point[1]),abs(point[2]),abs(point[3]))
    end
    tolerance=1e-12*max(scale,floatmin(Float64))
    for other in axes(controls,3-axis)
        a=axis==1 ? controls[1,other] : controls[other,1]
        b=axis==1 ? controls[last,other] : controls[other,last]
        wa=axis==1 ? weights[1,other] : weights[other,1]
        wb=axis==1 ? weights[last,other] : weights[other,last]
        abs(a[1]-b[1])<=tolerance && abs(a[2]-b[2])<=tolerance &&
            abs(a[3]-b[3])<=tolerance && abs(wa-wb)<=1e-12*max(abs(wa),abs(wb)) ||
            return false
    end
    return true
end

# Point and derivatives of the bare patch. The evaluation used to be written
# in `_model_nurbs_surface_d2` itself, where the results were assigned in two
# branches and captured by closures: every intermediate was boxed, and each
# homogeneous evaluation sliced the knot vectors and allocated a buffer per
# control row. That was several hundred heap objects per evaluated point,
# and a caller that wanted the point alone still paid for all six nets.
function _model_nurbs_patch_d2(patch::_ModelNURBSSurface,u::Float64,v::Float64)
    u,v=_model_nurbs_closed_parameters(patch,u,v)
    columns,buffer=_model_nurbs_scratch(patch)
    h=_model_nurbs_homogeneous(patch,1,u,v,columns,buffer)
    hu=_model_nurbs_homogeneous(patch,2,u,v,columns,buffer)
    hv=_model_nurbs_homogeneous(patch,3,u,v,columns,buffer)
    huu=_model_nurbs_homogeneous(patch,4,u,v,columns,buffer)
    hvv=_model_nurbs_homogeneous(patch,5,u,v,columns,buffer)
    huv=_model_nurbs_homogeneous(patch,6,u,v,columns,buffer)
    _model_nurbs_weight_check(h,u,v)
    point=ntuple(axis->h[axis]/h[4],3)
    du=ntuple(axis->(hu[axis]-hu[4]*point[axis])/h[4],3)
    dv=ntuple(axis->(hv[axis]-hv[4]*point[axis])/h[4],3)
    duu=ntuple(axis->(huu[axis]-huu[4]*point[axis]-2.0*hu[4]*du[axis])/h[4],3)
    dvv=ntuple(axis->(hvv[axis]-hvv[4]*point[axis]-2.0*hv[4]*dv[axis])/h[4],3)
    duv=ntuple(axis->(huv[axis]-huv[4]*point[axis]-hu[4]*dv[axis]-hv[4]*du[axis])/h[4],3)
    result=(point,du,dv,duu,dvv,duv)
    all(value->all(isfinite,value),result) || throw(ArgumentError(
        "NURBS Surface: evaluated derivatives must be finite"))
    return result
end
_model_nurbs_patch_d2(patch::_ModelNURBSSurface,u,v)=
    _model_nurbs_patch_d2(patch,Float64(u),Float64(v))

# The first three results of `_model_nurbs_patch_d2`, from three nets.
function _model_nurbs_patch_d1(patch::_ModelNURBSSurface,u::Float64,v::Float64)
    u,v=_model_nurbs_closed_parameters(patch,u,v)
    columns,buffer=_model_nurbs_scratch(patch)
    h=_model_nurbs_homogeneous(patch,1,u,v,columns,buffer)
    hu=_model_nurbs_homogeneous(patch,2,u,v,columns,buffer)
    hv=_model_nurbs_homogeneous(patch,3,u,v,columns,buffer)
    _model_nurbs_weight_check(h,u,v)
    point=ntuple(axis->h[axis]/h[4],3)
    du=ntuple(axis->(hu[axis]-hu[4]*point[axis])/h[4],3)
    dv=ntuple(axis->(hv[axis]-hv[4]*point[axis])/h[4],3)
    all(isfinite,point) && all(isfinite,du) && all(isfinite,dv) ||
        throw(ArgumentError("NURBS Surface: evaluated derivatives must be finite"))
    return point,du,dv
end

# The point alone, from the control net.
function _model_nurbs_patch_point(patch::_ModelNURBSSurface,u::Float64,v::Float64)
    u,v=_model_nurbs_closed_parameters(patch,u,v)
    columns,buffer=_model_nurbs_scratch(patch)
    h=_model_nurbs_homogeneous(patch,1,u,v,columns,buffer)
    _model_nurbs_weight_check(h,u,v)
    point=ntuple(axis->h[axis]/h[4],3)
    all(isfinite,point) || throw(ArgumentError(
        "NURBS Surface: evaluated derivatives must be finite"))
    return point
end

# The conic reparametrization of a record that carries one: the patch
# parameters with their first and second derivatives.
function _model_nurbs_reparam(g,u,v)
    ru,rv=g.reparam.u,g.reparam.v
    su,dsu,ddsu=ru===nothing ? (u,1.0,0.0) : _occ_conic_reparam(ru,u)
    sv,dsv,ddsv=rv===nothing ? (v,1.0,0.0) : _occ_conic_reparam(rv,v)
    return su,dsu,ddsu,sv,dsv,ddsv
end
@inline _model_nurbs_reparametrized(g)=hasproperty(g,:reparam) && g.reparam!==nothing

function _model_nurbs_surface_d2(g,u,v)
    _model_nurbs_reparametrized(g) || return _model_nurbs_patch_d2(g.patch,u,v)
    su,dsu,ddsu,sv,dsv,ddsv=_model_nurbs_reparam(g,u,v)
    point,du,dv,duu,dvv,duv=_model_nurbs_patch_d2(g.patch,su,sv)
    return (point,
        (du[1]*dsu,du[2]*dsu,du[3]*dsu),
        (dv[1]*dsv,dv[2]*dsv,dv[3]*dsv),
        (duu[1]*dsu*dsu+du[1]*ddsu,duu[2]*dsu*dsu+du[2]*ddsu,
         duu[3]*dsu*dsu+du[3]*ddsu),
        (dvv[1]*dsv*dsv+dv[1]*ddsv,dvv[2]*dsv*dsv+dv[2]*ddsv,
         dvv[3]*dsv*dsv+dv[3]*ddsv),
        (duv[1]*dsu*dsv,duv[2]*dsu*dsv,duv[3]*dsu*dsv))
end

# `(point, du, dv)` of a NURBS record: `_model_nurbs_surface_d2` without the
# three second-derivative nets.
function _model_nurbs_surface_d1(g,u,v)
    _model_nurbs_reparametrized(g) ||
        return _model_nurbs_patch_d1(g.patch,Float64(u),Float64(v))
    su,dsu,_,sv,dsv,_=_model_nurbs_reparam(g,u,v)
    point,du,dv=_model_nurbs_patch_d1(g.patch,Float64(su),Float64(sv))
    return point,(du[1]*dsu,du[2]*dsu,du[3]*dsu),(dv[1]*dsv,dv[2]*dsv,dv[3]*dsv)
end

# The point of a NURBS record: one net instead of six.
function _model_nurbs_surface_point(g,u,v)
    _model_nurbs_reparametrized(g) ||
        return _model_nurbs_patch_point(g.patch,Float64(u),Float64(v))
    su,_,_,sv,_,_=_model_nurbs_reparam(g,u,v)
    return _model_nurbs_patch_point(g.patch,Float64(su),Float64(sv))
end

function _model_nurbs_surface_bounds(g)
    s=g.patch.surface
    return [s.knots_u[s.degree_u+1],s.knots_v[s.degree_v+1]],
           [s.knots_u[end-s.degree_u],s.knots_v[end-s.degree_v]]
end

# Complete-domain seam certificate for trim/Boolean unwrap support. The
# positive-weight homogeneous curves are compared through exact Bernstein
# bounds, including nonclamped/rational supports and intermediate knot spans.
include("ModelNURBSSeams.jl")
# The certificate runs in exact rational arithmetic and depends only on the
# surface, the axis and the seam tolerance, so it is computed once per
# surface: ray containment asks it for every face hit (t20's STEP slice
# spent a third of its time re-certifying the same seams). Keyed weakly by
# the control net, like the projection grids; a hit is checked against the
# rest of the surface.
const _MODEL_NURBS_SEAM_PERIODS=WeakKeyDict{Matrix{NTuple{3,Float64}},
    Vector{Tuple{Any,Float64}}}()
const _MODEL_NURBS_SEAM_LOCK=ReentrantLock()

function _model_nurbs_patch_seam_period(g,dim::Int)
    surface=g.patch.surface
    controls=surface.controls
    controls isa Matrix{NTuple{3,Float64}} ||
        return _model_nurbs_certified_seam_period(g,dim)
    key=(dim,Float64(get(g,:seam_tolerance,1e-7)),surface.degree_u,
         surface.degree_v,surface.knots_u,surface.knots_v,surface.weights)
    same(k)=k[1]==key[1] && k[2]===key[2] && k[3]==key[3] && k[4]==key[4] &&
        (k[5]===key[5] || k[5]==key[5]) && (k[6]===key[6] || k[6]==key[6]) &&
        (k[7]===key[7] || k[7]==key[7])
    hit=lock(_MODEL_NURBS_SEAM_LOCK) do
        entries=get(_MODEL_NURBS_SEAM_PERIODS,controls,nothing)
        entries===nothing && return nothing
        for (k,p) in entries
            same(k) && return p
        end
        return nothing
    end
    hit===nothing || return hit
    p=_model_nurbs_certified_seam_period(g,dim)
    lock(_MODEL_NURBS_SEAM_LOCK) do
        entries=get!(()->Tuple{Any,Float64}[],_MODEL_NURBS_SEAM_PERIODS,controls)
        any(e->same(e[1]),entries) || push!(entries,(key,p))
    end
    return p
end

function _model_nurbs_knots(count,degree,knots,multiplicities,caller,axis)
    mult=Int[]
    for value in multiplicities
        value isa Integer && !(value isa Bool) || throw(ArgumentError(
            "$caller: $axis multiplicities must be integers"))
        push!(mult,Int(value))
    end
    if isempty(knots)
        size=count-degree+1
        unique_knots=Float64.(0:size-1)
        # Native vector::resize retains existing interior multiplicities,
        # initializes new slots to one, then replaces the two endpoints.
        resize!(mult,min(length(mult),size))
        append!(mult,ones(Int,size-length(mult)))
        mult[1]=degree+1;mult[end]=degree+1
    else
        unique_knots=Float64[_finite_scalar(k,caller,"$axis knot") for k in knots]
        length(unique_knots)==length(mult) || throw(ArgumentError(
            "$caller: $axis knots and multiplicities must have equal length"))
    end
    all(>(0),mult) || throw(ArgumentError(
        "$caller: $axis multiplicities must be positive integers"))
    length(unique_knots)>=2 && issorted(unique_knots) && allunique(unique_knots) ||
        throw(ArgumentError("$caller: $axis knots must be strictly increasing"))
    for (position,value) in enumerate(mult)
        limit=position in (1,length(mult)) ? degree+1 : degree
        value<=limit || throw(ArgumentError(
            "$caller: $axis knot multiplicity exceeds the native degree limit"))
    end
    sum(mult)==count+degree+1 || throw(ArgumentError(
        "$caller: $axis knot multiplicity sum must be $(count+degree+1)"))
    expanded=Float64[]
    for (k,m) in zip(unique_knots,mult)
        append!(expanded,fill(k,m))
    end
    return expanded
end

function _model_tensor_surface!(m,points,num_points_u,family,tag,degree_u,degree_v,
                                weights,knots_u,knots_v,multiplicities_u,
                                multiplicities_v,wire_tags,wire_3d,caller)
    num_points_u isa Integer && !(num_points_u isa Bool) && num_points_u>=2 ||
        throw(ArgumentError("$caller: need at least two control points along U"))
    nu=Int(num_points_u)
    tags=Int[_tag(point,caller,0) for point in points]
    length(tags)%nu==0 || throw(ArgumentError("$caller: wrong control grid shape"))
    nv=length(tags)÷nu
    nv>=2 || throw(ArgumentError("$caller: need at least two control points along V"))
    for point in tags
        haskey(m.points,point) || throw(ArgumentError("$caller: unknown control Point[$point]"))
    end
    isempty(wire_tags) || throw(ArgumentError(
        "$caller: trimmed tensor-product Surface wires are not implemented"))
    wire_3d isa Bool || throw(ArgumentError("$caller: wire_3d must be Bool"))
    grid=reshape(tags,nu,nv)
    family===:bspline && (all(grid[1,j]==grid[end,j] for j in 1:nv) ||
     all(grid[i,1]==grid[i,end] for i in 1:nu)) && throw(ArgumentError(
        "$caller: periodic same-Point-tag tensor-product Surfaces are not implemented"))
    for degree in (degree_u,degree_v)
        degree isa Integer && !(degree isa Bool) || throw(ArgumentError(
            "$caller: Surface degrees must be integers"))
    end
    pu=family===:bezier ? nu-1 : min(degree_u<=0 ? 3 : Int(degree_u),nu-1)
    pv=family===:bezier ? nv-1 : min(degree_v<=0 ? 3 : Int(degree_v),nv-1)
    pu<=25 && pv<=25 || throw(ArgumentError(
        "$caller: Surface degree exceeds the native maximum of 25"))
    U=family===:bezier ? vcat(zeros(pu+1),ones(pu+1)) :
        _model_nurbs_knots(nu,pu,knots_u,multiplicities_u,caller,"U")
    V=family===:bezier ? vcat(zeros(pv+1),ones(pv+1)) :
        _model_nurbs_knots(nv,pv,knots_v,multiplicities_v,caller,"V")
    W=isempty(weights) ? ones(nu,nv) : begin
        length(weights)==length(tags) || throw(ArgumentError("$caller: wrong weight count"))
        reshape(Float64[_finite_scalar(w,caller,"weight") for w in weights],nu,nv)
    end
    controls=map(point->m.points[point],grid)
    surface=NURBS.NURBSSurface(pu,pv,U,V,controls,W)
    patch=_ModelNURBSSurface(family,surface)
    requested=_tag(tag,caller,2)
    (requested>0 && (haskey(m.surfaces,requested) || haskey(m.discrete,(2,requested)))) &&
        throw(ArgumentError("$caller: Surface[$requested] already exists"))
    staged=_model_mesh_transaction_copy(m)
    lo=(U[pu+1],V[pv+1]);hi=(U[end-pu],V[end-pv])
    point_tags=Int[]
    for uv in (lo,(lo[1],hi[2]),(hi[1],lo[2]),hi)
        # The BRep vertices lie on the active-domain corners. Non-clamped
        # knot vectors do not interpolate the corner control poles.
        point=first(_model_nurbs_surface_d2((patch=patch,),uv...))
        boundary_point=add_point!(staged,point...;mesh_size=1.0)
        staged.point_size[boundary_point]=0.0
        push!(point_tags,boundary_point)
    end
    vertex_uv=Dict(point_tags[1]=>lo,point_tags[2]=>(lo[1],hi[2]),
                   point_tags[3]=>(hi[1],lo[2]),point_tags[4]=>hi)
    curve_tags=Int[]
    # Each boundary edge is an iso-parameter line of the patch, so its pcurve
    # is that line of the (u,v) square. With the pcurves on the record the
    # face is meshed in its own parameter domain like every other curved OCC
    # face; without them it reached the plane mesher and was rejected
    # ("plane geometry requires a Plane surface").
    pcurves=Dict{Int,NamedTuple}()
    geometry=(occ=family===:bezier ? :tensor_bezier : :tensor_bspline,
              patch=patch,vertex_uv=vertex_uv,reversed=false)
    endpoints=((1,2),(1,3),(3,4),(2,4))
    for side in 1:4
        curve=_alloc_tag!(staged,1,0,caller)
        a,b=endpoints[side]
        staged.curves[curve]=(point_tags[a],point_tags[b])
        staged.curve_types[curve]=family===:bezier ? :tensor_bezier : :tensor_bspline
        axis=side in (1,3) ? 2 : 1
        fixed=side in (1,2) ? lo[3-axis] : hi[3-axis]
        staged.curve_geometry[curve]=(occ=:surface_iso,face=geometry,axis=axis,
                                      fixed=fixed,t0=lo[axis],t1=hi[axis])
        pcurves[curve]=(fwd=(lin2d=axis==1 ? ((0.0,Float64(fixed)),(1.0,0.0)) :
                                             ((Float64(fixed),0.0),(0.0,1.0)),),
                        rev=nothing)
        push!(curve_tags,curve)
    end
    loop=add_curve_loop!(staged,[-curve_tags[1],curve_tags[2],curve_tags[3],-curve_tags[4]])
    surface_tag=_alloc_tag!(staged,2,requested,caller)
    staged.surfaces[surface_tag]=[loop]
    staged.surface_types[surface_tag]=family===:bezier ? :tensor_bezier : :tensor_bspline
    staged.surface_geometry[surface_tag]=merge(geometry,
        (boundary_curves=curve_tags,pcurves=pcurves))
    for field in fieldnames(GeoModel)
        setfield!(m,field,getfield(staged,field))
    end
    return surface_tag
end

"""Add an untrimmed nonperiodic Bezier Surface from a U-fast public Point grid."""
function add_bezier_surface!(m::GeoModel,points,num_points_u;tag::Integer=0,
                             wire_tags=Int[],wire_3d=false)
    return _model_tensor_surface!(m,points,num_points_u,:bezier,tag,3,3,
        Float64[],Float64[],Float64[],Int[],Int[],wire_tags,wire_3d,"add_bezier_surface!")
end

"""Add an untrimmed nonperiodic rational BSpline Surface from a U-fast Point grid."""
function add_bspline_surface!(m::GeoModel,points,num_points_u;tag::Integer=0,
                              degree_u=3,degree_v=3,weights=Float64[],
                              knots_u=Float64[],knots_v=Float64[],
                              multiplicities_u=Int[],multiplicities_v=Int[],
                              wire_tags=Int[],wire_3d=false)
    return _model_tensor_surface!(m,points,num_points_u,:bspline,tag,degree_u,degree_v,
        weights,knots_u,knots_v,multiplicities_u,multiplicities_v,wire_tags,wire_3d,
        "add_bspline_surface!")
end

_model_surface_iso_uv(g,u)=g.axis==1 ? (u,g.fixed) : (g.fixed,u)
_model_surface_iso_point(g,u)=first(_model_nurbs_surface_d2(g.face,_model_surface_iso_uv(g,u)...))
function _model_surface_iso_d1(g,u)
    values=_model_nurbs_surface_d2(g.face,_model_surface_iso_uv(g,u)...)
    return values[g.axis+1]
end
function _model_surface_iso_d2(g,u)
    values=_model_nurbs_surface_d2(g.face,_model_surface_iso_uv(g,u)...)
    return values[g.axis+3]
end

function _model_nurbs_surface_project(g,point)
    result=_model_nurbs_projection_cached(g.patch,point)
    u,v=result.uv[1],result.uv[2]
    if hasproperty(g,:reparam) && g.reparam!==nothing
        ru,rv=g.reparam.u,g.reparam.v
        u=ru===nothing ? u : _occ_conic_reparam_inv(ru,u)
        v=rv===nothing ? v : _occ_conic_reparam_inv(rv,v)
    end
    return u,v,result.xyz
end

function _model_nurbs_surface_support(m,surface,nodes,key,caller;
                                      quadrangle_boundary=nothing)
    g=m.surface_geometry[surface]
    uv=(0.0,0.0)
    for (position,weight) in key
        node=nodes[position]
        owner=node.owner
        parameter=if owner==(2,Int32(surface)) && node.parameter!==nothing &&
                     length(node.parameter)==2
            node.parameter
        elseif owner[1]==0
            model_reparametrize_on_surface(m,0,Int(owner[2]),Float64[],surface)
        elseif owner[1]==1 && node.parameter!==nothing && !isempty(node.parameter)
            model_reparametrize_on_surface(m,1,Int(owner[2]),[node.parameter[1]],surface)
        else
            project=_model_nurbs_surface_project(g,node.point)
            (project[1],project[2])
        end
        uv=(uv[1]+weight*parameter[1],uv[2]+weight*parameter[2])
    end
    target=ntuple(axis->sum(weight*nodes[position].point[axis]
                           for (position,weight) in key),3)
    if quadrangle_boundary!==nothing
        length(quadrangle_boundary)==8 || throw(ArgumentError(
            "$caller: quadratic quadrangle interior needs its eight boundary births"))
        target=ntuple(3) do axis
            result=0.0
            for slot in 1:8
                coefficient=slot<=4 ? -0.25 : 0.5
                result+=coefficient*nodes[quadrangle_boundary[slot]].point[axis]
            end
            result
        end
    end
    u,v,point=_model_nurbs_surface_project(g,target)
    return point,Float64[u,v]
end

# Native TF spacing is measured in the actual boundary vertex world frame,
# independently of the UV coordinates interpolated by TRAN_QUA.
function _model_tensor_transfinite_metric(first,opposite,caller)
    length(first)==length(opposite) || throw(ArgumentError(
        "$caller: tensor Surface opposite metric sides have unequal births"))
    result=zeros(Float64,length(first));total=0.0
    for interval in 1:length(first)-1
        a,b=first[interval],first[interval+1]
        dx=a[1]-b[1];dy=a[2]-b[2];dz=a[3]-b[3]
        d1=sqrt(dx*dx+dy*dy+dz*dz)
        a,b=opposite[interval],opposite[interval+1]
        dx=a[1]-b[1];dy=a[2]-b[2];dz=a[3]-b[3]
        d2=sqrt(dx*dx+dy*dy+dz*dz)
        total+=0.5*(d1+d2)
        result[interval+1]=total
    end
    isfinite(total) && total>0.0 || throw(ArgumentError(
        "$caller: tensor Surface boundary metric has no finite positive length"))
    for interval in 2:length(result)-1
        result[interval]/=total
    end
    result[end]=1.0
    return result
end

@inline function _model_tensor_transfinite_coons(c1,c2,c3,c4,s1,s2,s3,s4,u,v)
    return (1.0-u)*c4+u*c2+(1.0-v)*c1+v*c3-
        ((1.0-u)*(1.0-v)*s1+u*(1.0-v)*s2+u*v*s3+(1.0-u)*v*s4)
end

function _model_nurbs_transfinite_surface_mesh(m,t,param_sizes,caller;
        size_field=nothing,_mesh_input=nothing,_mesh_node_high_water=nothing,
        _return_surface_production=false,_geometry_points=nothing)
    spec,signed_curves,world_sides,nside,_,junctions,_,_,_=
        _transfinite_surface_sides(m,t,caller;_mesh_input)
    nside==4 || throw(ArgumentError("$caller: tensor Surface[$t] needs four TF sides"))
    if isempty(spec.corners)
        junctions,world_sides=_transfinite_chained_tri(m,signed_curves,world_sides,t,caller)
    end
    side_keys=_model_transfinite_side_keys(m,signed_curves,junctions,caller)
    uv_sides=Vector{Vector{NTuple{3,Float64}}}(undef,4)
    g=m.surface_geometry[t]
    for side in 1:4
        signed=side_keys[side];curve=abs(signed)
        geometry=m.curve_geometry[curve]
        source=_model_mesh_curve_births(_mesh_input,m,curve,caller)
        uv=NTuple{3,Float64}[]
        if source===nothing
            for point in world_sides[side]
                projected=_model_nurbs_surface_project(g,point)
                push!(uv,(projected[1],projected[2],0.0))
            end
        else
            rows=signed>0 ? source.rows : reverse(source.rows)
            length(rows)==length(world_sides[side]) || throw(ArgumentError(
                "$caller: tensor Surface side disagrees with actual Curve[$curve] births"))
            data=source.mesh.entity_data
            for (ordinal,row) in enumerate(rows)
                parameter=if ordinal==1
                    signed>0 ? geometry.t0 : geometry.t1
                elseif ordinal==length(rows)
                    signed>0 ? geometry.t1 : geometry.t0
                else
                    original=data.node_parametric[row]
                    original!==nothing && length(original)==1 || throw(ArgumentError(
                        "$caller: tensor Surface Curve[$curve] birth has no actual parameter"))
                    original[1]
                end
                u,v=_model_surface_iso_uv(geometry,parameter)
                push!(uv,(u,v,0.0))
            end
        end
        uv_sides[side]=uv
    end
    bottom,right,top,left=uv_sides
    side_tags=map(abs,side_keys)
    kernel=_model_surface_recombined(m,t) ?
        mesh_transfinite_quad_patch(bottom,right,top,left;
            side_tags,arrangement=spec.arrangement) :
        mesh_transfinite_patch(bottom,right,top,left;
            side_tags,arrangement=spec.arrangement)
    parameters=copy(kernel.coords)
    width=length(bottom);horizontal=width-1;vertical=length(right)-1
    metric_u=_model_tensor_transfinite_metric(world_sides[1],reverse(world_sides[3]),caller)
    metric_v=_model_tensor_transfinite_metric(world_sides[2],reverse(world_sides[4]),caller)
    c1=bottom[1];c2=bottom[end];c3=right[end];c4=top[end]
    for i in 1:horizontal-1,j in 1:vertical-1
        row=i+1+j*width
        for axis in 1:2
            parameters[axis,row]=_model_tensor_transfinite_coons(
                bottom[i+1][axis],right[j+1][axis],top[horizontal-i+1][axis],
                left[vertical-j+1][axis],c1[axis],c2[axis],c3[axis],c4[axis],
                metric_u[i+1],metric_v[j+1])
        end
    end
    for row in axes(kernel.coords,2)
        kernel.coords[:,row].=first(_model_nurbs_surface_d2(g,parameters[1,row],parameters[2,row]))
    end
    # The actual Curve rows are existing mesh vertices. Restore them by tagged
    # kernel events rather than by coordinate coincidence or private rank.
    for side in 1:4
        rows=Int32[]
        for (nodes,tags) in _model_transfinite_boundary_events(kernel)
            for column in axes(nodes,2)
                tags[column]==side_tags[side] || continue
                a,b=nodes[:,column]
                isempty(rows) ? push!(rows,a) : rows[end]==a || throw(ArgumentError(
                    "$caller: tensor Surface side events are not consecutive"))
                push!(rows,b)
            end
        end
        length(rows)==length(world_sides[side]) || throw(ArgumentError(
            "$caller: tensor Surface kernel changed the Curve birth count"))
        for (row,point) in zip(rows,world_sides[side])
            kernel.coords[:,row].=point
        end
    end
    mesh=kernel isa MixedMesh ? _surface_part_strip(kernel,caller,t) :
         Mesh(kernel.coords;tris=kernel.tris)
    _transfinite_quad_frame_reversed(m,signed_curves,junctions,caller,t) &&
        (mesh=_reversed_surface_mesh(mesh,caller,t))
    mesh=_consume_surface_attributes(m,t,mesh,caller)
    if _return_surface_production && _mesh_input!==nothing
        product=_model_transfinite_surface_production(m,t,kernel,mesh,_mesh_input,
            _mesh_node_high_water,param_sizes,caller;size_field,side_keys,
            _geometry_points=nothing)
        # meshGFaceTransfinite allocates interiors with i outer, j inner.
        # Map those producer events explicitly onto the kernel's row storage.
        fresh_count=count(!,product.retained_nodes)
        next_tag=product.node_high_water-UInt64(fresh_count)
        birth_rows=Int32.(findall(product.retained_nodes))
        for i in 1:horizontal-1,j in 1:vertical-1
            row=i+1+j*width
            product.retained_nodes[row] && throw(ArgumentError(
                "$caller: tensor Surface interior birth already has a boundary identity"))
            next_tag+=1;product.node_tags[row]=next_tag
            push!(birth_rows,Int32(row))
        end
        next_tag==product.node_high_water || throw(ArgumentError(
            "$caller: tensor Surface interior birth count disagrees with its node high water"))
        for row in eachindex(product.node_parameters)
            product.node_owners[row]==(2,Int32(t)) || continue
            product.node_parameters[row]=Float64[parameters[1,row],parameters[2,row]]
        end
        return _ModelSurfaceProduction(product.mesh,product.node_tags,
            product.node_owners,product.node_parameters,product.retained_nodes,
            product.node_high_water,product.element_tags,product.element_high_water,
            birth_rows)
    end
    return mesh
end
