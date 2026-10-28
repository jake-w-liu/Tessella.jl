# Operation-owned HighOrder.cpp edge/face factories. Identity is the public
# vertex tag; reference coordinates only describe a slot within that topology.
struct _NativeOrderFace
    corners::Vector{Int32}
    nodes::Dict{Tuple,Int32}
end

function _native_order_edges(family)
    family===:lin && return ((0,1),)
    return getfield(Elements,Symbol("_",uppercase(string(family)),"_EDGES"))
end

function _native_order_faces(family)
    family===:tri && return ((0,1,2,-1),)
    family===:qua && return ((0,1,2,3),)
    return getfield(Elements,Symbol("_",uppercase(string(family)),"_FACES"))
end

@inline function _native_order_edge_key(assembly,first,last)
    return minmax(assembly.nodes[first].tag,assembly.nodes[last].tag)
end

function _native_order_face_key(assembly,corners)
    return Tuple(sort!(UInt64[assembly.nodes[node].tag for node in corners]))
end

function _native_order_face_slot_key(assembly,corners,coordinate,order)
    u,v=coordinate
    weights=length(corners)==3 ? (order-u-v,u,v) :
        ((order-u)*(order-v),u*(order-v),u*v,(order-u)*v)
    return Tuple(sort!([(assembly.nodes[node].tag,weight)
                        for (node,weight) in zip(corners,weights)];by=first))
end

function _native_order_face_coordinates(family,order)
    coordinates=getfield(Elements,Symbol("_",family,"_monomials"))(order,false)
    return coordinates[_native_inner_boundary_count(family,order,length(coordinates))+1:end]
end

function _native_order_seed_edges!(edges,assembly,cell)
    spec=msh_spec(cell.msh)
    primary=Elements.msh_num_nodes(msh_type(spec.family,1))
    for (index,(a,b)) in enumerate(_native_order_edges(spec.family))
        first,last=cell.nodes[a+1],cell.nodes[b+1]
        key=_native_order_edge_key(assembly,first,last)
        haskey(edges,key) && continue
        start=primary+(index-1)*(spec.order-1)+1
        supports=cell.nodes[start:start+spec.order-2]
        assembly.nodes[first].tag<=assembly.nodes[last].tag || reverse!(supports)
        edges[key]=supports
    end
    return nothing
end

function _native_order_seed_surface!(faces,assembly,cell)
    spec=msh_spec(cell.msh)
    primary=Elements.msh_num_nodes(msh_type(spec.family,1))
    corners=cell.nodes[1:primary]
    key=_native_order_face_key(assembly,corners)
    haskey(faces,key) && return nothing
    slots=Dict{Tuple,Int32}()
    if !spec.serendipity
        first=primary+length(_native_order_edges(spec.family))*(spec.order-1)+1
        for (coordinate,node) in zip(_native_order_face_coordinates(spec.family,spec.order),
                                      @view cell.nodes[first:end])
            slots[_native_order_face_slot_key(assembly,corners,coordinate,spec.order)]=node
        end
    end
    faces[key]=_NativeOrderFace(corners,slots)
    return nothing
end

function _native_order_surface_uv(m,assembly,node,surface)
    source=assembly.nodes[node]
    dim,entity=source.owner
    parameter=source.parameter
    dim==2 && entity==surface && parameter!==nothing && length(parameter)>=2 &&
        return parameter[1:2]
    if dim==0 && haskey(m.points,Int(entity))
        return model_reparametrize_on_surface(m,0,Int(entity),Float64[],Int(surface))
    elseif dim==1 && haskey(m.curves,Int(entity))
        parameter!==nothing && !isempty(parameter) || return nothing
        return model_reparametrize_on_surface(m,1,Int(entity),parameter[1:1],Int(surface))
    end
    return model_parametrization(m,2,Int(surface),collect(source.point))
end

# HighOrder.cpp calls reparamMeshVertexOnFace for FACE interiors, while
# EDGE placement calls reparamMeshEdgeOnFace/getAllParameters. For a native
# GmshModel Point on a Plane only the FACE path uses edited mesh XYZ.
function _native_order_surface_face_uv(m,assembly,node,surface)
    source=assembly.nodes[node]
    dimension,entity=source.owner
    if dimension==0 && haskey(m.points,Int(entity)) &&
            Model.model_entity_type(m,2,Int(surface))=="Plane"
        kernel=m.meshing.native_geo_kernel
        origin=get(kernel.point_aliases,Int(entity),Int(entity))
        if haskey(kernel.points,origin)
            return model_parametrization(m,2,Int(surface),collect(source.point))
        end
    end
    return _native_order_surface_uv(m,assembly,node,surface)
end

function _native_order_edge_supports!(edges,coefficients,assembly,m,cell,order,linear)
    spec=msh_spec(cell.msh)
    primary=Elements.msh_num_nodes(msh_type(spec.family,1))
    supports=Int32[]
    for (index,(a,b)) in enumerate(_native_order_edges(spec.family))
        first,last=cell.nodes[a+1],cell.nodes[b+1]
        key=_native_order_edge_key(assembly,first,last)
        increasing=assembly.nodes[first].tag<=assembly.nodes[last].tag
        existing=get(edges,key,nothing)
        if existing!==nothing
            length(existing)==order-1 || throw(ArgumentError(
                "$(assembly.caller): shared edge has a different polynomial order"))
            append!(supports,increasing ? existing : reverse(existing))
            continue
        end
        start=primary+(index-1)*(spec.order-1)+1
        source=vcat(Int32[first,last],cell.nodes[start:start+spec.order-2])
        coefficient=get!(coefficients,spec.order) do
            _dim01_native_line_coefficients(spec.order)
        end
        geometric=!linear && cell.owner[1]==2 &&
            !_dim01_fully_discrete(m,(2,Int(cell.owner[2])))
        first_uv=geometric ? _native_order_surface_uv(m,assembly,first,cell.owner[2]) : nothing
        last_uv=geometric ? _native_order_surface_uv(m,assembly,last,cell.owner[2]) : nothing
        local_supports=Int32[]
        for slot in 1:order-1
            parameter=nothing
            if geometric && first_uv!==nothing && last_uv!==nothing
                # GFace edge subdivision uses equal parameter increments.
                ratio=Float64(slot)*(1.0/order)
                parameter=Float64[first_uv[axis]+(last_uv[axis]-first_uv[axis])*ratio for axis in 1:2]
                if Model.model_entity_type(m,2,Int(cell.owner[2]))!="Plane"
                    a=assembly.nodes[first].point;b=assembly.nodes[last].point
                    chord=Float64[ratio*b[axis]+(1.0-ratio)*a[axis] for axis in 1:3]
                    _,parameter=model_closest_point(m,2,Int(cell.owner[2]),chord)
                end
                values=model_value(m,2,Int(cell.owner[2]),parameter)
                point=(values[1],values[2],values[3])
            else
                coordinate=Float64(slot)*(2.0/order)-1.0
                weights=_dim01_native_line_weights(spec.order,coordinate,coefficient)
                point=(0.0,0.0,0.0)
                for (node,weight) in zip(source,weights)
                    value=assembly.nodes[node].point
                    point=(point[1]+weight*value[1],point[2]+weight*value[2],point[3]+weight*value[3])
                end
            end
            push!(local_supports,_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),
                point,cell.owner,parameter))
        end
        edges[key]=increasing ? copy(local_supports) : reverse(local_supports)
        append!(supports,local_supports)
    end
    return supports
end

function _native_order_inner_nodes!(assembly,m,owner,boundary,rows,geometric)
    supports=Int32[]
    uv=geometric ? [_native_order_surface_face_uv(m,assembly,node,owner[2]) for node in boundary] : nothing
    use_parameters=geometric && all(value->value!==nothing,uv)
    for row in rows
        point=(0.0,0.0,0.0)
        parameter=use_parameters ? zeros(Float64,2) : nothing
        for (column,weight) in row
            value=assembly.nodes[boundary[column]].point
            point=(point[1]+weight*value[1],point[2]+weight*value[2],point[3]+weight*value[3])
            if use_parameters
                parameter[1]+=weight*uv[column][1]
                parameter[2]+=weight*uv[column][2]
            end
        end
        if geometric
            if use_parameters && Model.model_entity_type(m,2,Int(owner[2]))=="Plane"
                values=model_value(m,2,Int(owner[2]),parameter)
            else
                values,projected=model_closest_point(m,2,Int(owner[2]),collect(point))
                use_parameters && (parameter=projected)
            end
            point=(values[1],values[2],values[3])
        end
        push!(supports,_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),point,owner,parameter))
    end
    return supports
end

function _native_order_face_boundary(assembly,corners,edges,order)
    family=length(corners)==3 ? :tri : :qua
    boundary=copy(corners)
    for (a,b) in _native_order_edges(family)
        first,last=corners[a+1],corners[b+1]
        values=edges[_native_order_edge_key(assembly,first,last)]
        append!(boundary,assembly.nodes[first].tag<=assembly.nodes[last].tag ? values : reverse(values))
    end
    return family,boundary
end

function _native_order_face_supports!(faces,edges,placement,assembly,m,cell,corners,order,linear)
    family=length(corners)==3 ? :tri : :qua
    coordinates=_native_order_face_coordinates(family,order)
    key=_native_order_face_key(assembly,corners)
    existing=get(faces,key,nothing)
    if existing!==nothing
        return Int32[existing.nodes[_native_order_face_slot_key(assembly,corners,q,order)] for q in coordinates]
    end
    _,boundary=_native_order_face_boundary(assembly,corners,edges,order)
    rows=get!(placement,(family,order)) do
        _native_inner_rows(family,order)
    end
    geometric=!linear && cell.owner[1]==2 && !_dim01_fully_discrete(m,(2,Int(cell.owner[2])))
    supports=_native_order_inner_nodes!(assembly,m,cell.owner,boundary,rows,geometric)
    slots=Dict(_native_order_face_slot_key(assembly,corners,q,order)=>node
        for (q,node) in zip(coordinates,supports))
    faces[key]=_NativeOrderFace(copy(corners),slots)
    return supports
end

function _dim01_higher_order!(assembly,m,entities,order;only_visible=false,linear=false,incomplete=false)
    # The verified Curve provider owns its geometry/interpolation pass.
    _dim01_curve_order!(assembly,m,entities,order;only_visible,linear)
    first_orders=Dict{Tuple{Int,Int32},Int}()
    for cell in assembly.cells
        cell.owner[1]>=2 && get!(first_orders,cell.owner,msh_spec(cell.msh).order)
    end
    edges=Dict{Tuple{UInt64,UInt64},Vector{Int32}}()
    faces=Dict{Tuple,_NativeOrderFace}()
    coefficients=Dict{Int,Matrix{Float64}}()
    placement=Dict{Tuple{Symbol,Int},Vector{Vector{Pair{Int,Float64}}}}()
    for cell in assembly.cells
        cell.owner[1]==1 || continue
        only_visible && model_entity_visibility(m,1,Int(cell.owner[2]))==0 && continue
        _native_order_seed_edges!(edges,assembly,cell)
    end
    for cell in assembly.cells
        cell.owner[1]>=2 || continue
        only_visible && model_entity_visibility(m,cell.owner[1],Int(cell.owner[2]))==0 && continue
        if first_orders[cell.owner]==order
            if cell.owner[1]==2
                _native_order_seed_edges!(edges,assembly,cell)
                _native_order_seed_surface!(faces,assembly,cell)
            end
            continue
        end
        spec=msh_spec(cell.msh)
        target=msh_type(spec.family,order;serendipity=incomplete)
        primary=Elements.msh_num_nodes(msh_type(spec.family,1))
        corners=cell.nodes[1:primary]
        boundary=vcat(corners,_native_order_edge_supports!(edges,coefficients,assembly,m,cell,order,linear))
        if !msh_spec(target).serendipity
            for face in _native_order_faces(spec.family)
                local_corners=Int32[corners[index+1] for index in face if index>=0]
                append!(boundary,_native_order_face_supports!(faces,edges,placement,assembly,m,
                    cell,local_corners,order,linear))
            end
            if cell.owner[1]==3
                rows=get!(placement,(spec.family,order)) do
                    _native_inner_rows(spec.family,order)
                end
                append!(boundary,_native_order_inner_nodes!(assembly,m,cell.owner,boundary,rows,false))
            end
        end
        length(boundary)==Elements.msh_num_nodes(target) || throw(ArgumentError(
            "$(assembly.caller): native MSH $target factory produced invalid connectivity"))
        cell.nodes=boundary
        cell.msh=target
        cell.tag=_dim01_next_tag!(assembly,false)
        # The lifted cell only carries a curved map once geometric placement
        # moved its high-order nodes; the whole-domain Jacobian bound certifies
        # it exactly like the P2 construction gate. Families whose upstream
        # determinant space does not exist stay uncertified, as upstream.
        _ho_constructed_cell_certify(assembly.nodes,cell,assembly.caller)
    end
    return nothing
end
