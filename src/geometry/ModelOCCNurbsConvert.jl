# `OCC_Internals::convertToNURBS` — `BRepBuilderAPI_NurbsConvert` applied per
# bound entity: geometry records are replaced by exact NURBS forms under the
# same tags and topology. Lines become degree-1 B-splines on the identical
# (arc-length) parametrization; circle and ellipse edges become quadratic-
# rational B-splines exposing their rational parameter. Conic and hyperbolic
# pcurves transport the edge parameter and surface chart exactly. Spline
# edges are already NURBS; vertices and degenerate/iso records pass through.
# Face, shell and solid closures stage their shared geometry together.

include("ModelOCCChartTransport.jl")

# Rational-quadratic image of a circle/ellipse arc: equal spans of at most
# 2π/3, interior poles on the tangent intersections `C + sec(Δ/2)·(p(m)−C)`,
# weights cos(Δ/2). `frame` resolves the eval-consistent record so
# native-backed circles convert through the frame evaluators actually use.
function _model_occ_nurbs_conic(geometry,caller)
    frame=_occ_2df_frame(geometry)
    C=frame.center
    eval=geometry.occ===:circle ? _occ_circle_point : _occ_ellipse_point
    t0,t1=geometry.t0,geometry.t1
    span=t1-t0
    span>0.0 || throw(ArgumentError("$caller: degenerate conic range"))
    n=max(1,ceil(Int,span/(2π/3)))
    δ=span/n
    weight=_gm_cos(δ/2)
    secant=1/weight
    poles=NTuple{3,Float64}[eval(geometry,t0)]
    weights=Float64[1.0]
    for i in 1:n
        a=t0+(i-1)*δ
        b=i==n ? t1 : t0+i*δ
        m=(a+b)/2
        pm=eval(geometry,m)
        push!(poles,ntuple(k->fma(pm[k]-C[k],secant,C[k]),3))
        push!(weights,weight)
        push!(poles,eval(geometry,b))
        push!(weights,1.0)
    end
    knots=Float64[t0,t0,t0]
    for i in 1:n-1
        knot=t0+i*δ
        push!(knots,knot,knot)
    end
    append!(knots,(t1,t1,t1))
    curve=NURBS.NURBSCurve(2,knots,poles,weights)
    return _ModelOCCNURBSCurve(:bspline,curve)
end

# Degree-1 B-spline on the identical arc-length parametrization — no
# reparametrization wrapper required.
function _model_occ_nurbs_line(geometry,caller)
    t0,t1=geometry.t0,geometry.t1
    t1>t0 || throw(ArgumentError("$caller: degenerate line range"))
    return _model_occ_curve_geometry(:bspline,
        NTuple{3,Float64}[_occ_line_point(geometry,t0),
                          _occ_line_point(geometry,t1)];
        degree=1,knots=[t0,t1],multiplicities=[2,2],caller=caller)
end

function _model_occ_nurbs_curve_record(geometry,caller)
    occ=get(geometry,:occ,nothing)
    if occ===:line
        # Analytic fields stay as inert metadata: eval dispatches on
        # `occ===:copied_curve` → `.curve`, while materialized-primitive
        # consistency checks can still read center/radius parameters.
        return merge(geometry,(occ=:copied_curve,converted_from=:line,
                curve=_model_occ_nurbs_line(geometry,caller)))
    elseif occ in (:circle,:ellipse)
        return merge(geometry,(occ=:copied_curve,converted_from=occ,
                curve=_model_occ_nurbs_conic(geometry,caller)))
    elseif occ===:copied_curve
        curve=geometry.curve
        if curve isa _ModelOCCHyperbolaCurve
            return merge(geometry,(curve=_model_occ_hyperbola_nurbs(curve,caller),))
        end
        curve isa _ModelOCCReparamCurve && return geometry
        # `GeomConvert` turns a Bezier into a single-span B-spline.
        curve.family===:bezier ||
            return geometry
        return merge(geometry,(curve=_ModelOCCNURBSCurve(:bspline,curve.curve),))
    elseif occ in (:degenerate,:surface_iso)
        # No standalone 3-D curve — the degenerate edge evaluates through its
        # owner face's pcurve and the iso edge through the surface record.
        return geometry
    end
    throw(ArgumentError(
        "$caller: OCC curve kind $(occ) cannot be converted to NURBS"))
end

function _model_occ_convert_to_nurbs!(m::GeoModel,dim_tags,caller;
        mesh_clear=false,mesh_topology=false)
    kernel=m.occ_kernel
    # Validate and stage everything before the first kernel mutation —
    # failures leave the registry byte-identical (a superset of upstream,
    # which converts greedily per entity).
    roots=UInt64[]
    for dt in dim_tags
        length(dt)==2 || throw(ArgumentError(
            "$caller: dimTag pairs must have two entries"))
        dim=_model_occ_dimension(dt[1],caller)
        tag=_model_occ_requested_tag(dt[2],caller)
        id=get(kernel.bindings,(dim,tag),nothing)
        id===nothing && throw(ArgumentError(
            "$caller: Unknown OpenCASCADE entity of dimension $dim with tag $tag"))
        push!(roots,id)
    end
    # Convert each edge/face record in the transitive shape closure once.
    staged=Dict{UInt64,NamedTuple}()
    curve_maps=Dict{UInt64,Union{_ModelOCCConicReparam,_ModelOCCHyperbolicReparam}}()
    surface_maps=Dict{UInt64,NamedTuple}()
    stack=copy(roots)
    seen=Set{UInt64}()
    order=UInt64[]
    while !isempty(stack)
        id=pop!(stack)
        id in seen && continue
        push!(seen,id)
        push!(order,id)
        shape=kernel.shapes[id]
        append!(stack,shape.children)
        if shape.dimension==1
            record=shape.geometry
            converted=_model_occ_nurbs_curve_record(record.geometry,caller)
            converted===record.geometry && continue
            staged[id]=merge(record,(kind=:tensor_bspline,geometry=converted))
            source=get(record.geometry,:curve,nothing)
            source isa _ModelOCCHyperbolaCurve &&
                (curve_maps[id]=_occ_hyperbolic_full_map(source.reparam))
            if record.geometry.occ in (:circle,:ellipse)
                curve_maps[id]=_model_occ_conic_chart(record.geometry.t0,record.geometry.t1)
            end
        elseif shape.dimension==2
            record=shape.geometry
            converted=_model_occ_nurbs_surface_record(record.geometry,caller)
            converted===record.geometry && continue
            staged[id]=merge(record,(geometry=converted,))
            mappings=get(converted,:analytic_chart,_occ_hyperbolic_surface_maps(record.geometry))
            mappings==(u=nothing,v=nothing) || (surface_maps[id]=mappings)
        end
    end
    # An explicitly converted edge exposes a new rational parameter. Every
    # referencing face keeps its physical trace, including a retained face
    # outside the requested closure. Converted surfaces also expose their
    # new rational UV coordinates.
    for (fid,face) in kernel.shapes
        face.dimension==2 || continue
        current=get(staged,fid,face.geometry)
        rec=current.geometry
        hasproperty(rec,:pcurves) || continue
        output=get(surface_maps,fid,(u=nothing,v=nothing))
        changed=output!=(u=nothing,v=nothing)
        pcurves=Dict{Int,NamedTuple}()
        for (tag,entry) in rec.pcurves
            eid=get(kernel.bindings,(1,abs(tag)),UInt64(0))
            input=get(curve_maps,eid,nothing)
            if input===nothing && eid!=0
                edge=kernel.shapes[eid].geometry.geometry
                if edge.occ===:surface_iso && hasproperty(edge,:face) &&
                        edge.face===face.geometry.geometry
                    input=edge.axis==1 ? output.u : output.v
                end
            end
            if input===nothing && output==(u=nothing,v=nothing)
                pcurves[tag]=entry
            else
                changed=true
                fwd=entry.fwd===nothing ? nothing : (brep_conversion2d=(entry.fwd,input,output),)
                rev=entry.rev===nothing ? nothing : (brep_conversion2d=(entry.rev,input,output),)
                pcurves[tag]=(fwd=fwd,rev=rev)
            end
        end
        if changed
            converted=merge(rec,(pcurves=pcurves,))
            if output!=(u=nothing,v=nothing) && hasproperty(rec,:vertex_uv)
                vertices=Dict{Int,NTuple{2,Float64}}()
                for (tag,uv) in rec.vertex_uv
                    vertices[tag]=ntuple(2) do axis
                        map=axis==1 ? output.u : output.v
                        map===nothing ? uv[axis] : first(_occ_conic_reparam(map,uv[axis]))
                    end
                end
                converted=merge(converted,(vertex_uv=vertices,))
            end
            staged[fid]=merge(current,(geometry=converted,))
        end
    end
    # Member degenerate/surface-iso edges embed their owner face record —
    # rebind each to the converted record so the pcurve lookup and the
    # surface evaluation stay on the same object.
    for id in order
        shape=kernel.shapes[id]
        shape.dimension==2 && haskey(staged,id) || continue
        newface=staged[id].geometry
        for wire_id in shape.children
            wire=kernel.shapes[wire_id]
            haskey(wire.geometry,:uses) || continue
            for (eid,_) in wire.geometry.uses
                eshape=kernel.shapes[eid]
                rec=eshape.geometry.geometry
                hasproperty(rec,:face) || continue
                entry=get(staged,eid,eshape.geometry)
                staged[eid]=merge(entry,(geometry=merge(
                    entry.geometry,(face=newface,)),))
            end
        end
    end
    for (id,record) in staged
        shape=kernel.shapes[id]
        kernel.shapes[id]=_ModelOCCShape(shape.dimension,record,shape.children)
        kernel.changed=true
        tag=get(kernel.shape_tags,id,nothing)
        tag===nothing && continue
        t=tag[2]
        if shape.dimension==1
            # The rebound edge is re-imported on the next synchronize — drop
            # the previously materialized model entries like the corner ops.
            delete!(m.curves,t);delete!(m.curve_types,t)
            delete!(m.curve_geometry,t)
            delete!(m.curve_params,t)
        elseif shape.dimension==2
            delete!(m.surfaces,t);delete!(m.surface_types,t)
            delete!(m.surface_geometry,t)
        end
    end
    return mesh_clear ? _model_occ_nurbs_mesh_selection(kernel,staged,order;topology=mesh_topology) :
        (entities=Tuple{Int,Int}[],upward=Dict{Tuple{Int,Int32},Vector{Tuple{Int,Int32}}}())
end

# Clear the converted closure and meshed owners of its changed subshapes.
# Other boundary products can survive when an unconverted owner still uses
# them; the native mesh-clear operation enforces that ownership rule.
function _model_occ_nurbs_mesh_selection(kernel,staged,order;topology=false)
    empty_upward=Dict{Tuple{Int,Int32},Vector{Tuple{Int,Int32}}}()
    isempty(staged) && return (entities=Tuple{Int,Int}[],upward=empty_upward)
    parents=Dict{UInt64,Vector{UInt64}}()
    for (id,shape) in kernel.shapes,child in shape.children
        push!(get!(()->UInt64[],parents,child),id)
    end
    selected=Set{UInt64}(order)
    pending=collect(keys(staged));visited=Set{UInt64}(pending)
    position=1
    while position<=length(pending)
        id=pending[position];position+=1
        push!(selected,id)
        for parent in get(parents,id,UInt64[])
            parent in visited && continue
            push!(visited,parent);push!(pending,parent)
        end
    end
    entities=Tuple{Int,Int}[]
    for id in selected
        entry=get(kernel.shape_tags,id,nothing)
        entry===nothing || !(0<=entry[1]<=3) || push!(entities,entry)
    end
    sort!(unique!(entities);by=entry->(-entry[1],entry[2]))
    topology || return (entities=entities,upward=empty_upward)
    upward=empty_upward
    for (dimension,tag) in entities
        id=kernel.bindings[(dimension,tag)]
        direct=Tuple{Int,Int32}[]
        queue=copy(get(parents,id,UInt64[]));seen=Set{UInt64}(queue)
        cursor=1
        while cursor<=length(queue)
            parent=queue[cursor];cursor+=1
            entry=get(kernel.shape_tags,parent,nothing)
            if entry!==nothing && entry[1]==dimension+1
                push!(direct,(entry[1],Int32(entry[2])))
                continue
            end
            for next in get(parents,parent,UInt64[])
                next in seen && continue
                push!(seen,next);push!(queue,next)
            end
        end
        upward[(dimension,Int32(tag))]=unique!(direct)
    end
    return (entities=entities,upward=upward)
end

# A CAD replacement has fresh meshing status. Mesh.clear itself deliberately
# retains status, so reset only replaced products that are now absent; shared
# retained products and unrelated entities keep their completed state.
function _model_occ_nurbs_mesh_status!(m,entities)
    authority=m.meshing.mesh_authority
    authority===nothing && return nothing
    for (dimension,tag) in entities
        owner=(dimension,Int32(tag))
        haskey(authority.products,owner) && continue
        delete!(authority.done,owner)
        delete!(authority.native_inputs,owner)
        delete!(authority.surface_constructor_inputs,owner)
        delete!(authority.producer_rows,owner)
    end
    return nothing
end

# ── surface conversion ─────────────────────────────────────────────────────
# Analytic quadric conversions retain their existing external chart through
# `patch`+`reparam`. Hyperbolic tensor axes expose the rational chart instead;
# their pcurves transport both the edge parameter and surface UV coordinates.

# One conic axis' pole layout over [t0,t1]: per-pole parameter value, radial
# factor (1 or sec(Δ/2)) and weight — shared by every quadratic-rational
# surface builder.
function _occ_nurbs_conic_axis(t0,t1)
    span=t1-t0
    n=max(1,ceil(Int,span/(2π/3)))
    δ=span/n
    secant=1/_gm_cos(δ/2); weight=_gm_cos(δ/2)
    vals=Float64[];qs=Float64[];ws=Float64[]
    for i in 0:2n
        if iseven(i)
            push!(vals,t0+(i>>1)*δ);push!(qs,1.0);push!(ws,1.0)
        else
            push!(vals,t0+(i>>1)*δ+δ/2);push!(qs,secant);push!(ws,weight)
        end
    end
    knots=Float64[t0,t0,t0]
    for i in 1:n-1
        knot=t0+i*δ
        push!(knots,knot,knot)
    end
    append!(knots,(t1,t1,t1))
    reparam=_ModelOCCConicReparam(Float64[t0+i*δ for i in 0:n],
        fill(_gm_tan(δ/4),n))
    return vals,qs,ws,knots,reparam
end

function _model_occ_nurbs_plane(g)
    lo,hi=g.uvb
    U=Float64[lo[1],lo[1],hi[1],hi[1]]
    V=Float64[lo[2],lo[2],hi[2],hi[2]]
    C,X,Y=g.center,g.X,_occ_frame_y(g)
    us=(lo[1],hi[1]);vs=(lo[2],hi[2])
    controls=Matrix{NTuple{3,Float64}}(undef,2,2)
    for i in 1:2,j in 1:2
        controls[i,j]=ntuple(k->fma(us[i],X[k],fma(vs[j],Y[k],C[k])),3)
    end
    patch=_ModelNURBSSurface(:bspline,
        NURBS.NURBSSurface(1,1,U,V,controls,ones(2,2)))
    return patch,nothing
end

function _model_occ_nurbs_cylinder(g)
    angle=get(g,:angle,2π)
    vals,qs,ws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    V=Float64[0.0,0.0,g.height,g.height]
    C,X,Y,Z=g.center,g.X,_occ_frame_y(g),g.axis
    r=g.radius
    nu=length(vals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,2)
    W=Matrix{Float64}(undef,nu,2)
    for i in 1:nu,j in 1:2
        v=j==1 ? 0.0 : g.height
        s,c=_gm_sincos(vals[i])
        controls[i,j]=ntuple(k->fma(qs[i]*r*c,X[k],
            fma(qs[i]*r*s,Y[k],fma(v,Z[k],C[k]))),3)
        W[i,j]=ws[i]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,1,U,V,controls,W))
    return patch,(u=urep,v=nothing)
end

function _model_occ_nurbs_cone(g)
    angle=get(g,:angle,2π)
    vals,qs,ws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    slant=sqrt(g.height*g.height+(g.r1-g.r2)*(g.r1-g.r2))
    V=Float64[0.0,0.0,slant,slant]
    sa=_gm_atan((g.r2-g.r1)/g.height);ssa,csa=_gm_sincos(sa)
    C,X,Y,Z=g.center,g.X,_occ_frame_y(g),g.axis
    nu=length(vals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,2)
    W=Matrix{Float64}(undef,nu,2)
    for i in 1:nu,j in 1:2
        v=j==1 ? 0.0 : slant
        R=fma(v,ssa,g.r1);A3=v*csa
        s,c=_gm_sincos(vals[i])
        controls[i,j]=ntuple(k->fma(qs[i]*R*c,X[k],
            fma(qs[i]*R*s,Y[k],fma(A3,Z[k],C[k]))),3)
        W[i,j]=ws[i]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,1,U,V,controls,W))
    return patch,(u=urep,v=nothing)
end

function _model_occ_nurbs_sphere(g)
    angle=get(g,:angle,2π)
    lat1,lat2=get(g,:lat1,-π/2),get(g,:lat2,π/2)
    uvals,uqs,uws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    vvals,vqs,vws,V,vrep=_occ_nurbs_conic_axis(lat1,lat2)
    C,X,Y,Z=g.center,g.X,g.Y,g.axis
    r=g.radius
    nu,nv=length(uvals),length(vvals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    W=Matrix{Float64}(undef,nu,nv)
    for i in 1:nu,j in 1:nv
        su,cu=_gm_sincos(uvals[i]);sv,cv=_gm_sincos(vvals[j])
        # Row j is a circle of radius r·qv_j·cos(v̂_j) centered at
        # C + r·qv_j·sin(v̂_j)·Z — the mid-v factor applies to both
        # components, the mid-u factor only to the radial direction.
        radial=r*vqs[j]*cv*uqs[i]
        controls[i,j]=ntuple(k->fma(radial*cu,X[k],
            fma(radial*su,Y[k],fma(r*vqs[j]*sv,Z[k],C[k]))),3)
        W[i,j]=uws[i]*vws[j]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,2,U,V,controls,W))
    return patch,(u=urep,v=vrep)
end

function _model_occ_nurbs_torus(g)
    angle=get(g,:angle,2π)
    uvals,uqs,uws,U,urep=_occ_nurbs_conic_axis(0.0,angle)
    vvals,vqs,vws,V,vrep=_occ_nurbs_conic_axis(0.0,2π)
    C,X,Y,Z=g.center,g.X,g.Y,g.axis
    r1,r2=g.r1,g.r2
    nu,nv=length(uvals),length(vvals)
    controls=Matrix{NTuple{3,Float64}}(undef,nu,nv)
    W=Matrix{Float64}(undef,nu,nv)
    for i in 1:nu,j in 1:nv
        su,cu=_gm_sincos(uvals[i]);sv,cv=_gm_sincos(vvals[j])
        rho=uqs[i]*(r1+r2*vqs[j]*cv)
        A3=r2*vqs[j]*sv
        controls[i,j]=ntuple(k->fma(rho*cu,X[k],
            fma(rho*su,Y[k],fma(A3,Z[k],C[k]))),3)
        W[i,j]=uws[i]*vws[j]
    end
    patch=_ModelNURBSSurface(:bspline,NURBS.NURBSSurface(2,2,U,V,controls,W))
    return patch,(u=urep,v=vrep)
end

function _model_occ_nurbs_surface_record(geometry,caller)
    occ=get(geometry,:occ,nothing)
    if hasproperty(geometry,:patch)
        mappings=_occ_hyperbolic_surface_maps(geometry)
        mappings==(u=nothing,v=nothing) ||
            return _model_occ_hyperbolic_surface_nurbs(geometry,caller)
        # Already-NURBS tensor record: a bezier patch renormalizes to a
        # B-spline like `GeomConvert`; anything else passes through.
        geometry.patch.family===:bezier || return geometry
        return merge(geometry,(occ=:tensor_bspline,
            patch=_ModelNURBSSurface(:bspline,geometry.patch.surface),))
    end
    patch,reparam=if occ===:plane
        _model_occ_nurbs_plane(geometry)
    elseif occ===:cylinder
        _model_occ_nurbs_cylinder(geometry)
    elseif occ===:cone
        _model_occ_nurbs_cone(geometry)
    elseif occ===:sphere
        _model_occ_nurbs_sphere(geometry)
    elseif occ===:torus
        _model_occ_nurbs_torus(geometry)
    else
        throw(ArgumentError(
            "$caller: OCC surface kind $(occ) cannot be converted to NURBS"))
    end
    return reparam===nothing ? merge(geometry,(patch=patch,)) :
        merge(geometry,(patch=patch,analytic_chart=reparam))
end
