"""
    MeshReferenceGeometry

Robust forward reference maps and Gmsh-shaped Jacobian arrays for finalized
linear-simplex [`Mesh`](@ref) values. Segment, triangle, and tetrahedron maps use
the Gmsh 4.15.2 reference conventions and retain signed tetrahedron determinants.
Fast floating filters fall back to exact-rational and `BigFloat` evaluation when
ordinary arithmetic cannot certify a finite, nondegenerate result.
"""
module MeshReferenceGeometry

using ..MeshTypes: Mesh, node
using ..MeshPointLocation: mesh_element_block, mesh_element_offsets
using ..Elements: msh_spec

export mesh_jacobian, mesh_jacobians

const _Exact=Rational{BigInt}
const _FILTER_FACTOR=32eps(Float64)
const _AFFINE_FILTER_FACTOR=64eps(Float64)

@inline _cross3(a,b)=(a[2]*b[3]-a[3]*b[2],
                      a[3]*b[1]-a[1]*b[3],
                      a[1]*b[2]-a[2]*b[1])
@inline _dot3(a,b)=muladd(a[1],b[1],muladd(a[2],b[2],a[3]*b[3]))
@inline _det3(a,b,c)=_dot3(a,_cross3(b,c))
@inline _finite3(a)=isfinite(a[1])&&isfinite(a[2])&&isfinite(a[3])

function _checked_integer(value,caller::AbstractString,name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    return try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError ||
         err isa MethodError) || rethrow()
        throw(ArgumentError(
            "$caller: $name exceeds the platform Int range"))
    end
end

function _checked_element_type(value,caller::AbstractString)
    element_type=_checked_integer(value,caller,"element_type")
    msh_spec(element_type)
    return element_type
end

function _checked_element_tag(mesh::Mesh,value,caller::AbstractString)
    tag=_checked_integer(value,caller,"element_tag")
    _,_,total=mesh_element_offsets(mesh)
    1<=tag<=total || throw(ArgumentError(
        "$caller: unknown element tag $value; expected a dense tag in 1:$total"))
    return tag
end

function _checked_local_coordinate(value,index::Int,caller::AbstractString)
    value isa Real || throw(ArgumentError(
        "$caller: local_coord[$index] must be a real number"))
    value isa Bool && throw(ArgumentError(
        "$caller: local_coord[$index] must not be Bool"))
    converted=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError ||
         err isa MethodError) || rethrow()
        throw(ArgumentError(
            "$caller: local_coord[$index] must be Float64-representable"))
    end
    isfinite(converted) || throw(ArgumentError(
        "$caller: local_coord[$index] must be finite"))
    return converted
end

function _checked_local_coordinates(values,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: local_coord must be a vector or tuple of real numbers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    coordinate_count=length(values)
    coordinate_count%3==0 || throw(ArgumentError(
        "$caller: local_coord length must be divisible by 3; got " *
        "$coordinate_count"))
    result=Vector{Float64}(undef,coordinate_count)
    for (index,value) in enumerate(values)
        result[index]=_checked_local_coordinate(value,index,caller)
    end
    return result,coordinate_count÷3
end

function _result_lengths(element_count::Int,point_count::Int,
                         caller::AbstractString)
    try
        evaluation_count=Base.checked_mul(element_count,point_count)
        return Base.checked_mul(9,evaluation_count),evaluation_count,
               Base.checked_mul(3,evaluation_count)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "$caller: requested Jacobian result exceeds the platform array range"))
    end
end

function _allocate_results(element_count::Int,point_count::Int,
                           caller::AbstractString)
    jacobian_count,determinant_count,coordinate_count=
        _result_lengths(element_count,point_count,caller)
    return Vector{Float64}(undef,jacobian_count),
           Vector{Float64}(undef,determinant_count),
           Vector{Float64}(undef,coordinate_count)
end

function _exact_to_float(value::_Exact,caller::AbstractString,
                         description::AbstractString)
    result=setprecision(BigFloat,256) do
        Float64(BigFloat(value))
    end
    isfinite(result) || throw(ArgumentError(
        "$caller: $description is not Float64-representable"))
    return result
end

function _derivative_description(kind::Symbol,tag::Int,component::Int)
    kind===:segment && return "element $tag segment derivative component $component"
    kind===:triangle_u && return "element $tag triangle u-derivative component $component"
    kind===:triangle_v && return "element $tag triangle v-derivative component $component"
    return "element $tag tetrahedron derivative component $component"
end

function _difference(a::Float64,b::Float64,caller::AbstractString,
                     tag::Int,kind::Symbol,component::Int)
    value=b-a
    (isfinite(value) && (value!=0.0 || a==b)) && return value
    exact=_Exact(b)-_Exact(a)
    description=_derivative_description(kind,tag,component)
    result=_exact_to_float(exact,caller,description)
    result==0.0 && exact!=0 && throw(ArgumentError(
        "$caller: $description is nonzero but below Float64 resolution"))
    return result
end

function _half_difference(a::Float64,b::Float64,caller::AbstractString,
                          tag::Int,component::Int)
    difference=b-a
    if isfinite(difference)
        value=difference/2
        (value!=0.0 || difference==0.0) && return value
    end
    exact=(_Exact(b)-_Exact(a))/2
    description=_derivative_description(:segment,tag,component)
    result=_exact_to_float(exact,caller,description)
    result==0.0 && exact!=0 && throw(ArgumentError(
        "$caller: $description is nonzero but below Float64 resolution"))
    return result
end

@inline function _stable_midpoint(a::Float64,b::Float64)
    return signbit(a)==signbit(b) ? a+(b-a)/2 : a/2+b/2
end

function _segment_frame(a,b,caller::AbstractString,tag::Int)
    derivative=ntuple(3) do axis
        _half_difference(
            a[axis],b[axis],caller,tag,axis)
    end
    determinant=hypot(derivative...)
    isfinite(determinant) || throw(ArgumentError(
        "$caller: element $tag segment determinant is not " *
        "Float64-representable"))
    determinant>0 || throw(ArgumentError(
        "$caller: element $tag is a degenerate segment"))

    x,y,z=derivative
    transverse=if ((abs(x)>=abs(y) && abs(x)>=abs(z)) ||
                   (abs(y)>=abs(x) && abs(y)>=abs(z)))
        (y,-x,0.0)
    else
        (0.0,z,-y)
    end
    transverse_norm=hypot(transverse...)
    transverse_norm>0 || error(
        "$caller: internal segment regularization failure for element $tag")
    transverse=ntuple(i->transverse[i]/transverse_norm,3)
    tangent=ntuple(i->derivative[i]/determinant,3)
    normal=_cross3(tangent,transverse)
    normal_norm=hypot(normal...)
    normal_norm>0 || error(
        "$caller: internal segment normal failure for element $tag")
    normal=ntuple(i->normal[i]/normal_norm,3)
    jacobian=(derivative...,transverse...,normal...)
    all(isfinite,jacobian) || error(
        "$caller: internal non-finite segment frame for element $tag")
    return jacobian,determinant,derivative
end

function _exact_triangle_normal(a,b,c,caller::AbstractString,tag::Int)
    ae=ntuple(i->_Exact(a[i]),3)
    edge1=ntuple(i->_Exact(b[i])-ae[i],3)
    edge2=ntuple(i->_Exact(c[i])-ae[i],3)
    cross=_cross3(edge1,edge2)
    squared=_dot3(cross,cross)
    squared>0 || throw(ArgumentError(
        "$caller: element $tag is a degenerate triangle"))
    return setprecision(BigFloat,256) do
        magnitude=sqrt(BigFloat(squared))
        determinant=Float64(magnitude)
        isfinite(determinant) || throw(ArgumentError(
            "$caller: element $tag triangle determinant is not " *
            "Float64-representable"))
        determinant==0.0 && throw(ArgumentError(
            "$caller: element $tag triangle determinant is nonzero but " *
            "below Float64 resolution"))
        normal=ntuple(3) do axis
            Float64(BigFloat(cross[axis])/magnitude)
        end
        _finite3(normal) || error(
            "$caller: internal non-finite triangle normal for element $tag")
        normal,determinant
    end
end

function _triangle_normal(a,b,c,edge1,edge2,
                          caller::AbstractString,tag::Int)
    cross=_cross3(edge1,edge2)
    safe=true
    @inbounds for (component,i,j,k,l) in
        ((1,2,3,3,2),(2,3,1,1,3),(3,1,2,2,1))
        permanent=abs(edge1[i]*edge2[j])+abs(edge1[k]*edge2[l])
        safe &= isfinite(permanent) &&
                (permanent==0.0 ||
                 abs(cross[component])>_FILTER_FACTOR*permanent)
    end
    magnitude=hypot(cross...)
    if safe && isfinite(magnitude) && magnitude>0
        return ntuple(i->cross[i]/magnitude,3),magnitude
    end
    return _exact_triangle_normal(a,b,c,caller,tag)
end

function _triangle_frame(a,b,c,caller::AbstractString,tag::Int)
    edge1=ntuple(3) do axis
        _difference(
            a[axis],b[axis],caller,tag,:triangle_u,axis)
    end
    edge2=ntuple(3) do axis
        _difference(
            a[axis],c[axis],caller,tag,:triangle_v,axis)
    end
    normal,determinant=
        _triangle_normal(a,b,c,edge1,edge2,caller,tag)
    return (edge1...,edge2...,normal...),determinant,(edge1,edge2)
end

function _exact_tetrahedron_determinant(a,b,c,d,
                                        caller::AbstractString,tag::Int)
    ae=ntuple(i->_Exact(a[i]),3)
    edges=ntuple(3) do local_node
        point=(b,c,d)[local_node]
        ntuple(i->_Exact(point[i])-ae[i],3)
    end
    determinant=_det3(edges...)
    determinant!=0 || throw(ArgumentError(
        "$caller: element $tag is a degenerate tetrahedron"))
    result=_exact_to_float(
        determinant,caller,"element $tag tetrahedron determinant")
    result==0.0 && throw(ArgumentError(
        "$caller: element $tag tetrahedron determinant is nonzero but " *
        "below Float64 resolution"))
    return result
end

function _tetrahedron_determinant(a,b,c,d,edge1,edge2,edge3,
                                  caller::AbstractString,tag::Int)
    determinant=_det3(edge1,edge2,edge3)
    permanent=
        abs(edge1[1])*(abs(edge2[2]*edge3[3])+abs(edge2[3]*edge3[2]))+
        abs(edge1[2])*(abs(edge2[1]*edge3[3])+abs(edge2[3]*edge3[1]))+
        abs(edge1[3])*(abs(edge2[1]*edge3[2])+abs(edge2[2]*edge3[1]))
    if isfinite(determinant) && isfinite(permanent) &&
       abs(determinant)>_FILTER_FACTOR*permanent
        return determinant
    end
    return _exact_tetrahedron_determinant(a,b,c,d,caller,tag)
end

function _tetrahedron_frame(a,b,c,d,caller::AbstractString,tag::Int)
    points=(b,c,d)
    edges=ntuple(3) do local_node
        point=points[local_node]
        ntuple(3) do axis
            _difference(
                a[axis],point[axis],caller,tag,:tetrahedron,
                3(local_node-1)+axis)
        end
    end
    determinant=_tetrahedron_determinant(
        a,b,c,d,edges...,caller,tag)
    return (edges[1]...,edges[2]...,edges[3]...),determinant,edges
end

function _affine_axis_fast(origin::Float64,edges,parameters)
    value=origin
    bound=abs(origin)
    @inbounds for index in eachindex(edges,parameters)
        product=parameters[index]*edges[index]
        isfinite(product) || return 0.0,false
        (product!=0.0 || parameters[index]==0.0 || edges[index]==0.0) ||
            return 0.0,false
        value=muladd(parameters[index],edges[index],value)
        bound+=abs(product)
    end
    isfinite(value) && isfinite(bound) &&
        (bound==0.0 || abs(value)>_AFFINE_FILTER_FACTOR*bound) ||
        return 0.0,false
    return value,true
end

function _exact_physical_axis(vertices,parameters,::Val{dimension},axis::Int) where dimension
    if dimension==1
        midpoint=(_Exact(vertices[1][axis])+_Exact(vertices[2][axis]))/2
        derivative=(_Exact(vertices[2][axis])-
                    _Exact(vertices[1][axis]))/2
        return midpoint+_Exact(parameters[1])*derivative
    end
    anchor=_Exact(vertices[1][axis])
    value=anchor
    @inbounds for local_node in 1:dimension
        edge=_Exact(vertices[local_node+1][axis])-anchor
        value+=_Exact(parameters[local_node])*edge
    end
    return value
end

function _physical_point(vertices,origin,edges,parameters,dimension::Val{D},
                         caller::AbstractString,tag::Int) where D
    if D==1
        parameters[1]==-1.0 && return vertices[1]
        parameters[1]==1.0 && return vertices[2]
    end
    static_dimension=Val(D)
    return ntuple(3) do axis
        axis_edges=ntuple(index->edges[index][axis],static_dimension)
        value,certified=_affine_axis_fast(
            origin[axis],axis_edges,parameters)
        certified && return value
        exact=_exact_physical_axis(vertices,parameters,dimension,axis)
        result=_exact_to_float(
            exact,caller,
            "element $tag physical coordinate component $axis")
        result==0.0 && exact!=0 && throw(ArgumentError(
            "$caller: element $tag physical coordinate component $axis " *
            "is nonzero but below Float64 resolution"))
        result
    end
end

function _write_evaluations!(jacobians,determinants,coordinates,
                             element_index::Int,point_count::Int,
                             local_coordinates,jacobian,determinant,
                             vertices,origin,edges,dimension::Val{D},
                             caller::AbstractString,tag::Int) where D
    @inbounds for point in 1:point_count
        evaluation=(element_index-1)*point_count+point
        jacobian_offset=9(evaluation-1)
        for component in 1:9
            jacobians[jacobian_offset+component]=jacobian[component]
        end
        determinants[evaluation]=determinant
        local_offset=3(point-1)
        parameters=ntuple(
            index->local_coordinates[local_offset+index],dimension)
        physical=_physical_point(
            vertices,origin,edges,parameters,dimension,caller,tag)
        coordinate_offset=3(evaluation-1)
        coordinates[coordinate_offset+1]=physical[1]
        coordinates[coordinate_offset+2]=physical[2]
        coordinates[coordinate_offset+3]=physical[3]
    end
    return nothing
end

function _write_cell!(jacobians,determinants,coordinates,
                      mesh::Mesh,cells::Matrix{Int32},cell::Int,
                      element_index::Int,tag::Int,point_count::Int,
                      local_coordinates,caller::AbstractString,
                      dimension::Val{1})
    @inbounds a=node(mesh,cells[1,cell]);b=node(mesh,cells[2,cell])
    jacobian,determinant,derivative=_segment_frame(a,b,caller,tag)
    origin=ntuple(i->_stable_midpoint(a[i],b[i]),3)
    return _write_evaluations!(
        jacobians,determinants,coordinates,element_index,point_count,
        local_coordinates,jacobian,determinant,(a,b),origin,
        (derivative,),dimension,caller,tag)
end

function _write_cell!(jacobians,determinants,coordinates,
                      mesh::Mesh,cells::Matrix{Int32},cell::Int,
                      element_index::Int,tag::Int,point_count::Int,
                      local_coordinates,caller::AbstractString,
                      dimension::Val{2})
    @inbounds a=node(mesh,cells[1,cell]);b=node(mesh,cells[2,cell])
    @inbounds c=node(mesh,cells[3,cell])
    jacobian,determinant,edges=_triangle_frame(a,b,c,caller,tag)
    return _write_evaluations!(
        jacobians,determinants,coordinates,element_index,point_count,
        local_coordinates,jacobian,determinant,(a,b,c),a,edges,
        dimension,caller,tag)
end

function _write_cell!(jacobians,determinants,coordinates,
                      mesh::Mesh,cells::Matrix{Int32},cell::Int,
                      element_index::Int,tag::Int,point_count::Int,
                      local_coordinates,caller::AbstractString,
                      dimension::Val{3})
    @inbounds a=node(mesh,cells[1,cell]);b=node(mesh,cells[2,cell])
    @inbounds c=node(mesh,cells[3,cell])
    @inbounds d=node(mesh,cells[4,cell])
    jacobian,determinant,edges=
        _tetrahedron_frame(a,b,c,d,caller,tag)
    return _write_evaluations!(
        jacobians,determinants,coordinates,element_index,point_count,
        local_coordinates,jacobian,determinant,(a,b,c,d),a,edges,
        dimension,caller,tag)
end

"""
    mesh_jacobians(mesh, element_type, local_coord)

Return detached `(jacobians, determinants, coordinates)` for every cached element
of one Gmsh type at concatenated `(u,v,w)` evaluation points. Results are ordered
by element and then point. Jacobian matrices use Gmsh's column-flattened layout.
The finalized cache contains only types 1, 2, and 4; another known type returns
three empty vectors. Malformed coordinates, degenerate maps, and non-finite or
unrepresentable results fail explicitly.
"""
function mesh_jacobians(mesh::Mesh,element_type,local_coord)
    caller="mesh_jacobians"
    msh=_checked_element_type(element_type,caller)
    local_coordinates,point_count=
        _checked_local_coordinates(local_coord,caller)
    block=mesh_element_block(mesh,msh)
    (block===nothing || point_count==0) &&
        return Float64[],Float64[],Float64[]
    offset,cells=block
    element_count=size(cells,2)
    element_count==0 && return Float64[],Float64[],Float64[]
    jacobians,determinants,coordinates=
        _allocate_results(element_count,point_count,caller)
    @inbounds for cell in axes(cells,2)
        _write_cell!(
            jacobians,determinants,coordinates,mesh,cells,cell,cell,
            offset+cell,point_count,local_coordinates,caller,
            Val(msh==1 ? 1 : msh==2 ? 2 : 3))
    end
    return jacobians,determinants,coordinates
end

"""
    mesh_jacobian(mesh, element_tag, local_coord)

Return detached Gmsh-shaped Jacobians, determinants, and physical coordinates for
one dense cached linear-simplex element tag. Evaluation-point and error contracts
match [`mesh_jacobians`](@ref).
"""
function mesh_jacobian(mesh::Mesh,element_tag,local_coord)
    caller="mesh_jacobian"
    tag=_checked_element_tag(mesh,element_tag,caller)
    local_coordinates,point_count=
        _checked_local_coordinates(local_coord,caller)
    point_count==0 && return Float64[],Float64[],Float64[]
    triangle_offset,tetrahedron_offset,_=mesh_element_offsets(mesh)
    if tag<=triangle_offset
        cells=mesh.segs;cell=tag
    elseif tag<=tetrahedron_offset
        cells=mesh.tris;cell=tag-triangle_offset
    else
        cells=mesh.tets;cell=tag-tetrahedron_offset
    end
    jacobians,determinants,coordinates=
        _allocate_results(1,point_count,caller)
    _write_cell!(
        jacobians,determinants,coordinates,mesh,cells,cell,1,tag,
        point_count,local_coordinates,caller,
        Val(tag<=triangle_offset ? 1 : tag<=tetrahedron_offset ? 2 : 3))
    return jacobians,determinants,coordinates
end

end # module MeshReferenceGeometry
