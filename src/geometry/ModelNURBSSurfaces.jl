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

function _model_nurbs_deboor_h(points,degree,knots,u)
    count=length(points)
    span=clamp(searchsortedlast(knots,u)-1,degree,count-1)
    buffer=NTuple{4,Float64}[points[span-degree+j+1] for j in 0:degree]
    for level in 1:degree,j in degree:-1:level
        lower=knots[span-degree+j+1]
        upper=knots[span+j-level+2]
        factor=(u-lower)/(upper-lower)
        left=buffer[j];right=buffer[j+1]
        buffer[j+1]=ntuple(axis->(1.0-factor)*left[axis]+factor*right[axis],4)
    end
    return buffer[end]
end

function _model_nurbs_homogeneous(patch,slot,u,v)
    net=patch.nets[slot]
    net===nothing && return (0.0,0.0,0.0,0.0)
    ou,ov=slot==1 ? (0,0) : slot==2 ? (1,0) : slot==3 ? (0,1) :
          slot==4 ? (2,0) : slot==5 ? (0,2) : (1,1)
    surface=patch.surface
    U=surface.knots_u[1+ou:end-ou];V=surface.knots_v[1+ov:end-ov]
    columns=NTuple{4,Float64}[
        _model_nurbs_deboor_h(@view(net[i,:]),surface.degree_v-ov,V,v)
        for i in axes(net,1)]
    return _model_nurbs_deboor_h(columns,surface.degree_u-ou,U,u)
end

function _model_nurbs_surface_d2(g,u,v)
    patch=g.patch
    h=_model_nurbs_homogeneous(patch,1,u,v)
    hu=_model_nurbs_homogeneous(patch,2,u,v)
    hv=_model_nurbs_homogeneous(patch,3,u,v)
    huu=_model_nurbs_homogeneous(patch,4,u,v)
    hvv=_model_nurbs_homogeneous(patch,5,u,v)
    huv=_model_nurbs_homogeneous(patch,6,u,v)
    isfinite(h[4]) && !iszero(h[4]) || throw(ArgumentError(
        "NURBS Surface: zero or non-finite homogeneous weight at ($u,$v)"))
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

function _model_nurbs_surface_bounds(g)
    s=g.patch.surface
    return [s.knots_u[s.degree_u+1],s.knots_v[s.degree_v+1]],
           [s.knots_u[end-s.degree_u],s.knots_v[end-s.degree_v]]
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
        push!(curve_tags,curve)
    end
    loop=add_curve_loop!(staged,[-curve_tags[1],curve_tags[2],curve_tags[3],-curve_tags[4]])
    surface_tag=_alloc_tag!(staged,2,requested,caller)
    staged.surfaces[surface_tag]=[loop]
    staged.surface_types[surface_tag]=family===:bezier ? :tensor_bezier : :tensor_bspline
    staged.surface_geometry[surface_tag]=merge(geometry,(boundary_curves=curve_tags,))
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
    return result.uv[1],result.uv[2],result.xyz
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
