# `discreteEdge::_discretization`/`_pars` are geometry, separate from MVertex
# parameters. Capturing geometry never changes a mesh vertex's native data.
function _model_discrete_curve_geometry(m::GeoModel,tag::Integer)
    geometry=get(m.curve_geometry,Int(tag),nothing)
    return geometry!==nothing && hasproperty(geometry,:discrete_curve) &&
           geometry.discrete_curve ? geometry : nothing
end

# GEdge::length(..., 10) evaluates the unsplit piecewise derivative at the
# literal ten-point rule. Splitting at geometric knots changes HighOrder.cpp's
# equidistant-parameter Newton iteration.
const _DISCRETE_CURVE_GL10_POINTS = (
    -9.739065285171716e-01,-8.650633666889845e-01,-6.794095682990244e-01,
    -4.333953941292472e-01,-1.488743389816312e-01,1.488743389816312e-01,
    4.333953941292472e-01,6.794095682990244e-01,8.650633666889845e-01,
    9.739065285171716e-01)
const _DISCRETE_CURVE_GL10_WEIGHTS = (
    6.667134430868774e-02,1.494513491505805e-01,2.190863625159822e-01,
    2.692667193099962e-01,2.955242247147529e-01,2.955242247147529e-01,
    2.692667193099962e-01,2.190863625159822e-01,1.494513491505805e-01,
    6.667134430868774e-02)

function _model_discrete_curve_length10(geometry,u0,u1,caller)
    length=0.0
    jacobian=(u1-u0)*0.5
    for index in 1:10
        abscissa=_DISCRETE_CURVE_GL10_POINTS[index]
        parameter=u0*0.5*(1.0-abscissa)+u1*0.5*(1.0+abscissa)
        derivative=_discrete_curve_derivative(geometry.frames,parameter,caller)
        norm=sqrt(derivative[1]*derivative[1]+derivative[2]*derivative[2]+
                  derivative[3]*derivative[3])
        length+=norm*_DISCRETE_CURVE_GL10_WEIGHTS[index]*jacobian
    end
    return length
end

function _model_discrete_curve_equidistant(geometry,u0,u1,relaxation,caller)
    increment=(u1-u0)/2.0
    parameter=u0+increment
    upper=parameter+increment
    epsilon=(u1-u0)*1e-5
    for iteration in 2:50
        residual=_model_discrete_curve_length10(geometry,parameter,upper,caller)-
                 _model_discrete_curve_length10(geometry,u0,parameter,caller)
        parameter+=epsilon
        perturbed=_model_discrete_curve_length10(geometry,parameter,upper,caller)-
                  _model_discrete_curve_length10(geometry,u0,parameter,caller)
        jacobian=(perturbed-residual)/epsilon
        parameter-=epsilon
        correction=residual/jacobian
        parameter-=relaxation*correction
        (parameter<u0 || parameter>u1) && break
        abs(correction)<1e-6 && return true,parameter
    end
    return false,parameter
end

function _model_discrete_curve_p2_parameter(geometry,first,last,caller)
    lower,upper=minmax(first,last)
    relaxation=1.0
    while relaxation>1e-2
        converged,parameter=_model_discrete_curve_equidistant(
            geometry,lower,upper,relaxation,caller)
        converged && return parameter
        relaxation*=0.5
    end
    return lower+(upper-lower)/2.0
end

function _model_discrete_curve_segment_middle(points,index)
    first,last=points[index],points[index+1]
    return ((first[1]+last[1])*0.5,(first[2]+last[2])*0.5,(first[3]+last[3])*0.5)
end

function _model_discrete_curve_distance(first,last)
    x=first[1]-last[1];y=first[2]-last[2];z=first[3]-last[3]
    return sqrt(x*x+y*y+z*z)
end

function _model_discrete_curve_curvature(geometry,parameter,closed,caller)
    points=geometry.points
    length(points)<=3 && return 0.0
    edge=_discrete_curve_frame_index(geometry.frames,parameter)
    edge==0 && return 0.0
    previous=edge-1;following=edge+1
    if edge==1
        edge+=1;previous+=1;following+=1
        closed && (previous=length(points)-1)
    elseif edge==length(points)-1
        edge-=1;previous-=1;following-=1
        closed && (following=1)
    end
    a=_model_discrete_curve_segment_middle(points,previous)
    b=_model_discrete_curve_segment_middle(points,edge)
    c=_model_discrete_curve_segment_middle(points,following)
    A=_model_discrete_curve_distance(b,c)
    B=_model_discrete_curve_distance(c,a)
    C=_model_discrete_curve_distance(a,b)
    radicand=(A+B+C)*(-A+B+C)*(A-B+C)*(A+B-C)
    radius=A*B*C/(radicand<0.0 ? NaN : sqrt(radicand))
    return 1.0/radius
end

function _model_discrete_curve_second_derivative(geometry,parameter,caller)
    if parameter-1e-3<=geometry.t0
        first=_discrete_curve_derivative(geometry.frames,parameter,caller)
        last=_discrete_curve_derivative(geometry.frames,parameter+1e-3,caller)
        scale=1000.0
    elseif parameter+1e-3>=geometry.t1
        first=_discrete_curve_derivative(geometry.frames,parameter-1e-3,caller)
        last=_discrete_curve_derivative(geometry.frames,parameter,caller)
        scale=1000.0
    else
        first=_discrete_curve_derivative(geometry.frames,parameter-1e-3,caller)
        last=_discrete_curve_derivative(geometry.frames,parameter+1e-3,caller)
        scale=500.0
    end
    return (scale*(last[1]-first[1]),scale*(last[2]-first[2]),scale*(last[3]-first[3]))
end

function _model_discrete_curve_refine_projection(geometry,target,parameter,relaxation,lc,caller)
    point=_discrete_curve_point(geometry.frames,parameter,caller)
    difference=(point[1]-target[1],point[2]-target[2],point[3]-target[3])
    error=_model_discrete_curve_distance(point,target)
    threshold=1e-8*lc
    for iteration in 1:25
        error>threshold || break
        derivative=_discrete_curve_derivative(geometry.frames,parameter,caller)
        correction=(difference[1]*derivative[1]+difference[2]*derivative[2]+
                    difference[3]*derivative[3])/
                   (derivative[1]*derivative[1]+derivative[2]*derivative[2]+
                    derivative[3]*derivative[3])
        correction<1e-8 && error>threshold && (correction=1.0)
        parameter=min(geometry.t1,max(geometry.t0,parameter-relaxation*correction))
        point=_discrete_curve_point(geometry.frames,parameter,caller)
        difference=(point[1]-target[1],point[2]-target[2],point[3]-target[3])
        error=_model_discrete_curve_distance(point,target)
    end
    return error<=threshold,parameter,error
end

function _model_discrete_curve_xyz_to_u(geometry,target,lc,caller;relaxation=1.0)
    errors=Dict{Float64,Float64}()
    for index in 0:20
        parameter=geometry.t0+(geometry.t1-geometry.t0)/20*index
        converged,result,error=_model_discrete_curve_refine_projection(
            geometry,target,parameter,relaxation,lc,caller)
        converged && return true,result
        errors[error]=result
    end
    if relaxation>0.1
        converged,result=_model_discrete_curve_xyz_to_u(
            geometry,target,lc,caller;relaxation=0.75*relaxation)
        converged && return true,result
        point=_discrete_curve_point(geometry.frames,result,caller)
        errors[_model_discrete_curve_distance(point,target)]=result
    end
    return false,errors[minimum(keys(errors))]
end

function _model_discrete_curve_golden_section(geometry,target,first,middle,last,caller)
    golden2=2.0-(1.0+sqrt(5.0))/2.0
    trial=middle+golden2*(last-middle)
    abs(last-first)<1e-9*(abs(middle)+abs(trial)) && return (last+first)/2
    trial_distance=_model_discrete_curve_distance(
        target,_discrete_curve_point(geometry.frames,trial,caller))
    middle_distance=_model_discrete_curve_distance(
        target,_discrete_curve_point(geometry.frames,middle,caller))
    return trial_distance<middle_distance ?
        _model_discrete_curve_golden_section(geometry,target,middle,trial,last,caller) :
        _model_discrete_curve_golden_section(geometry,target,trial,middle,first,caller)
end

function _model_discrete_curve_closest_point(geometry,target,caller)
    lower,upper=minmax(geometry.t0,geometry.t1)
    increment=(upper-lower)/99.0
    distance=1e22;best=lower
    for index in 0:99
        parameter=lower+index*increment
        point=_discrete_curve_point(geometry.frames,parameter,caller)
        current=_model_discrete_curve_distance(target,point)
        if current<distance
            best=parameter;distance=current
        end
    end
    parameter=best==lower ? _model_discrete_curve_golden_section(
        geometry,target,best,best+increment/2,best+increment,caller) :
        best==upper ? _model_discrete_curve_golden_section(
        geometry,target,best-increment,best-increment/2,best,caller) :
        _model_discrete_curve_golden_section(
        geometry,target,best-increment,best,best+increment,caller)
    return _discrete_curve_point(geometry.frames,parameter,caller),parameter
end

function _model_discrete_curve_sources(m::GeoModel,caller)
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    owners=Dict{UInt64,Tuple{Int,Int}}()
    authority=m.meshing.mesh_authority
    if authority!==nothing && authority.mesh!==nothing
        actual=authority.mesh;data=actual.entity_data
        for (row,external) in enumerate(data.external_node_tags)
            tag=UInt64(external)
            coordinates[tag]=(actual.coords[1,row],actual.coords[2,row],actual.coords[3,row])
            owner=data.node_entities[row];owners[tag]=(owner[1],Int(owner[2]))
        end
    end
    # `authority.mesh` is a published snapshot — `classify_surfaces!`
    # redistributes nodes across new discrete entities without republishing
    # it, so its owner rows can lag the records. The record claims are the
    # current ownership and must win (the curve-parametrization plan reads
    # them for endpoint Point entities); coordinates still must agree
    # wherever the tag appears, and two records may not claim the same tag.
    claimed=Dict{UInt64,Tuple{Int,Int}}()
    for (owner,record) in _discrete_mesh_records_model(m)
        for (column,tag) in enumerate(record.node_tags)
            point=(record.node_coords[1,column],record.node_coords[2,column],record.node_coords[3,column])
            if haskey(claimed,tag)
                (coordinates[tag]==point && claimed[tag]==owner) ||
                    throw(ArgumentError(
                        "$caller: node $tag has conflicting actual geometry sources"))
            else
                haskey(coordinates,tag) && coordinates[tag]!=point &&
                    throw(ArgumentError(
                        "$caller: node $tag has conflicting actual geometry sources"))
                coordinates[tag]=point
                claimed[tag]=owner
                owners[tag]=owner
            end
        end
    end
    return coordinates,owners
end

# Literal `SortEdgeConsecutive`: the fallback comparator is the current public
# MVertex number. A consecutively oriented Line vector bypasses this sorting.
function _model_discrete_curve_chain(edges,caller,tag)
    chain=UInt64[first(edges)[1]]
    sorted=true
    for (a,b) in edges
        if a==last(chain)
            push!(chain,b)
        else
            sorted=false
            break
        end
    end
    sorted && return chain
    neighbors=Dict{UInt64,Vector{UInt64}}()
    for (a,b) in edges
        for (node,other) in ((a,b),(b,a))
            adjacent=get!(()->UInt64[],neighbors,node)
            length(adjacent)<2 || throw(ArgumentError(
                "$caller: discrete Curve[$tag] has a branched Line topology"))
            push!(adjacent,other)
        end
    end
    chains=Vector{UInt64}[]
    while !isempty(neighbors)
        ordered=sort!(collect(keys(neighbors)))
        endpoint=findfirst(node->length(neighbors[node])==1,ordered)
        start=endpoint===nothing ? first(ordered) : ordered[endpoint]
        adjacent=neighbors[start]
        previous=length(adjacent)==1 ? nothing : adjacent[2]==start ? adjacent[1] : adjacent[2]
        current=start;path=UInt64[]
        while true
            haskey(neighbors,current) || throw(ArgumentError(
                "$caller: discrete Curve[$tag] has an invalid Line chain"))
            push!(path,current)
            adjacent=pop!(neighbors,current)
            first_neighbor=adjacent[1]
            second_neighbor=length(adjacent)==1 ? nothing : adjacent[2]
            next=first_neighbor==previous ? second_neighbor :
                 second_neighbor==previous ? first_neighbor : nothing
            previous=current
            next===nothing && break
            current=next
            if current==start
                push!(path,current)
                break
            end
        end
        length(path)>2 && path[end-1]==path[end] && pop!(path)
        push!(chains,path)
    end
    length(chains)==1 || throw(ArgumentError(
        "$caller: discrete Curve[$tag] has disconnected Line chains"))
    return only(chains)
end

function _model_discrete_curve_geometry_plan(m,tag,record,coordinates,owners,caller)
    existing=_model_discrete_curve_geometry(m,tag)
    existing!==nothing && return existing,copy(record.boundary)
    edges=NTuple{2,UInt64}[]
    for (index,msh) in enumerate(record.element_types)
        msh_family(msh)===:lin || continue
        nodes=record.element_nodes[index]
        push!(edges,(nodes[1],nodes[2]))
    end
    isempty(edges) && return nothing,copy(record.boundary)
    chain=_model_discrete_curve_chain(edges,caller,tag)
    points=NTuple{3,Float64}[]
    for node in chain
        haskey(coordinates,node) || throw(ArgumentError(
            "$caller: discrete Curve[$tag] references unknown actual node $node"))
        push!(points,coordinates[node])
    end
    declared=Int[entity for (dimension,entity) in record.boundary if dimension==0]
    first_owner=get(owners,first(chain),(-1,0))
    last_owner=get(owners,last(chain),(-1,0))
    begin_point=first_owner[1]==0 ? first_owner[2] : isempty(declared) ? 0 : first(declared)
    end_point=last_owner[1]==0 ? last_owner[2] : isempty(declared) ? 0 : last(declared)
    if begin_point==0 && last_owner[1]==0 && _model_point_distance(first(points),last(points))<1e-14
        begin_point=last_owner[2]
    end
    if end_point==0 && first_owner[1]==0 && _model_point_distance(first(points),last(points))<1e-14
        end_point=first_owner[2]
    end
    begin_point!=0 && end_point!=0 || throw(ArgumentError(
        "$caller: discrete Curve[$tag] has no begin or end Point"))
    parameters=Float64.(0:length(points)-1)
    frames=Tuple{Vector{Vector{Float64}},Vector{NTuple{3,Float64}}}[
        ([Float64[parameters[index]],Float64[parameters[index+1]]],
         [points[index],points[index+1]]) for index in 1:length(points)-1]
    geometry=(discrete_curve=true,points=points,parameters=parameters,frames=frames,
              t0=first(parameters),t1=last(parameters))
    return geometry,NTuple{2,Int}[(0,begin_point),(0,end_point)]
end
