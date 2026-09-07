"""
    MeshFunctionSpaces

First-order reference finite-element function spaces and global degree-of-freedom
keys for finalized linear-simplex [`Mesh`](@ref) values. The implementation follows
Gmsh 4.15.2's reference coordinates, output layout, lexicographic orientation
indices, nodal keys, and lowest-order edge keys. Unsupported higher-order and mixed
families fail explicitly.
"""
module MeshFunctionSpaces

using ..MeshTypes: Mesh
using ..MeshEntityTopology: MeshEdgeTopology, _mesh_edges,
                            _simplex_edge_patterns
using ..MeshPointLocation: mesh_element_block, mesh_element_record
using ..MeshReferenceGeometry: _Exact, _checked_element_tag,
                               _checked_element_type,
                               _checked_local_coordinates, _exact_to_float,
                               _stable_midpoint

export mesh_basis_functions, mesh_basis_orientation,
       mesh_basis_orientations, mesh_key_dimension, mesh_keys,
       mesh_keys_for_element, mesh_keys_information,
       mesh_number_of_keys, mesh_number_of_orientations

struct _FunctionSpace
    kind::Symbol
    components::Int
    hierarchical::Bool
    key_dimension::Int
    key_order::Int
end

const _LAGRANGE=_FunctionSpace(:lagrange,1,false,0,-1)
const _LAGRANGE_1=_FunctionSpace(:lagrange,1,false,0,1)
const _GRAD_LAGRANGE=_FunctionSpace(:grad_lagrange,3,false,0,-1)
const _GRAD_LAGRANGE_1=_FunctionSpace(:grad_lagrange,3,false,0,1)
const _H1_1=_FunctionSpace(:lagrange,1,true,0,1)
const _GRAD_H1_1=_FunctionSpace(:grad_lagrange,3,true,0,1)
const _HCURL_0=_FunctionSpace(:hcurl,3,true,1,0)
const _CURL_HCURL_0=_FunctionSpace(:curl_hcurl,3,true,1,0)

function _function_space(value,caller::AbstractString)
    value isa AbstractString || throw(ArgumentError(
        "$caller: function_space_type must be a string"))
    name=String(value)
    occursin('\0',name) && throw(ArgumentError(
        "$caller: function_space_type must not contain NUL"))
    (name=="Lagrange" || name=="IsoParametric") && return _LAGRANGE
    name=="Lagrange1" && return _LAGRANGE_1
    (name=="GradLagrange" || name=="GradIsoParametric") &&
        return _GRAD_LAGRANGE
    name=="GradLagrange1" && return _GRAD_LAGRANGE_1
    name=="H1Legendre1" && return _H1_1
    name=="GradH1Legendre1" && return _GRAD_H1_1
    name=="HcurlLegendre0" && return _HCURL_0
    name=="CurlHcurlLegendre0" && return _CURL_HCURL_0
    throw(ArgumentError(
        "$caller: unsupported function_space_type $(repr(name)); supported " *
        "first-order spaces are Lagrange, IsoParametric, Lagrange1, " *
        "GradLagrange, GradIsoParametric, GradLagrange1, H1Legendre1, " *
        "GradH1Legendre1, HcurlLegendre0, and CurlHcurlLegendre0"))
end

function _linear_simplex_type(value,caller::AbstractString)
    element_type=_checked_element_type(value,caller)
    element_type in (1,2,4) || throw(ArgumentError(
        "$caller: element type $element_type is not a supported linear segment, " *
        "triangle, or tetrahedron"))
    return element_type
end

@inline function _vertex_count(element_type::Int)
    element_type==1 && return 2
    element_type==2 && return 3
    return 4
end

@inline function _edge_count(element_type::Int)
    element_type==1 && return 1
    element_type==2 && return 3
    return 6
end

@inline _orientation_count(element_type::Int)=
    element_type==1 ? 2 : element_type==2 ? 6 : 24

function _checked_orientation_sequence(values,norientations::Int,
                                       hierarchical::Bool,
                                       caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: wanted_orientations must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    if isempty(values)
        return hierarchical ? collect(0:norientations-1) : Int[0]
    end
    result=Vector{Int}(undef,length(values))
    seen=Set{Int}()
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: wanted_orientations[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: wanted_orientations[$index] must not be Bool"))
        orientation=try
            Int(value)
        catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: wanted_orientations[$index] exceeds Int bounds"))
        end
        0<=orientation<norientations || throw(ArgumentError(
            "$caller: orientation $orientation is outside " *
            "0:$(norientations-1)"))
        orientation in seen && throw(ArgumentError(
            "$caller: duplicate wanted orientation $orientation"))
        push!(seen,orientation)
        result[index]=orientation
    end
    return result
end

function _checked_result_length(factors::Int...;caller::AbstractString)
    result=1
    try
        for factor in factors
            result=Base.checked_mul(result,factor)
        end
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError(
            "$caller: requested basis result exceeds the platform array range"))
    end
    return result
end

function _origin_barycentric(values,caller::AbstractString,point::Int)
    fast=1.0
    scale=1.0
    for value in values
        fast-=value
        scale=max(scale,abs(value))
    end
    if isfinite(fast) && abs(fast)>64eps(Float64)*scale
        return fast
    end
    exact=_Exact(1)
    for value in values
        exact-=_Exact(value)
    end
    result=_exact_to_float(
        exact,caller,"evaluation point $point origin barycentric coordinate")
    result==0.0 && exact!=0 && throw(ArgumentError(
        "$caller: evaluation point $point origin barycentric coordinate is " *
        "nonzero but below Float64 resolution"))
    return result
end

function _barycentric_coordinates(element_type::Int,u::Float64,v::Float64,
                                  w::Float64,caller::AbstractString,point::Int)
    if element_type==1
        return ((1.0-u)/2,(1.0+u)/2)
    elseif element_type==2
        return (_origin_barycentric((u,v),caller,point),u,v)
    end
    return (_origin_barycentric((u,v,w),caller,point),u,v,w)
end

@inline function _reference_gradients(element_type::Int)
    element_type==1 && return ((-0.5,0.0,0.0),(0.5,0.0,0.0))
    element_type==2 && return (
        (-1.0,-1.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0))
    return ((-1.0,-1.0,-1.0),(1.0,0.0,0.0),
            (0.0,1.0,0.0),(0.0,0.0,1.0))
end

@inline function _canonical_edge_pairs(element_type::Int)
    element_type==1 && return ((1,2),)
    element_type==2 && return ((1,2),(2,3),(1,3))
    return ((1,2),(2,3),(1,3),(1,4),(3,4),(2,4))
end

function _orientation_permutation(vertex_count::Int,orientation::Int)
    remaining=collect(1:vertex_count)
    permutation=Vector{Int}(undef,vertex_count)
    rank=orientation
    @inbounds for position in 1:vertex_count
        block=factorial(vertex_count-position)
        selection=rank÷block+1
        rank%=block
        permutation[position]=remaining[selection]
        deleteat!(remaining,selection)
    end
    return permutation
end

function _exact_barycentric_coordinate(element_type::Int,vertex::Int,
                                       u::Float64,v::Float64,w::Float64)
    if element_type==1
        return vertex==1 ? (_Exact(1)-_Exact(u))/2 :
                           (_Exact(1)+_Exact(u))/2
    elseif element_type==2
        vertex==1 && return _Exact(1)-_Exact(u)-_Exact(v)
        return vertex==2 ? _Exact(u) : _Exact(v)
    end
    vertex==1 && return _Exact(1)-_Exact(u)-_Exact(v)-_Exact(w)
    vertex==2 && return _Exact(u)
    vertex==3 && return _Exact(v)
    return _Exact(w)
end

function _scaled_whitney_component(element_type::Int,
                                   first_vertex::Int,lambda_i::Float64,
                                   gradient_i::Float64,
                                   second_vertex::Int,lambda_j::Float64,
                                   gradient_j::Float64,
                                   u::Float64,v::Float64,w::Float64,
                                   caller::AbstractString,point::Int,
                                   edge::Int,component::Int)
    first_term=lambda_i*gradient_j
    second_term=lambda_j*gradient_i
    difference=first_term-second_term
    value=2.0*difference
    scale=max(abs(first_term),abs(second_term))
    if isfinite(value) && (scale==0.0 ||
                           abs(difference)>64eps(Float64)*scale)
        return value
    end
    exact_i=_exact_barycentric_coordinate(
        element_type,first_vertex,u,v,w)
    exact_j=_exact_barycentric_coordinate(
        element_type,second_vertex,u,v,w)
    exact=2*(exact_i*_Exact(gradient_j)-exact_j*_Exact(gradient_i))
    result=_exact_to_float(
        exact,caller,"evaluation point $point edge $edge component $component")
    result==0.0 && exact!=0 && throw(ArgumentError(
        "$caller: evaluation point $point edge $edge component $component is " *
        "nonzero but below Float64 resolution"))
    return result
end

@inline function _curl_component(gradient_i,gradient_j,component::Int)
    component==1 && return 4.0*(gradient_i[2]*gradient_j[3]-
                               gradient_i[3]*gradient_j[2])
    component==2 && return 4.0*(gradient_i[3]*gradient_j[1]-
                               gradient_i[1]*gradient_j[3])
    return 4.0*(gradient_i[1]*gradient_j[2]-
                gradient_i[2]*gradient_j[1])
end

"""Return first-order Gmsh-shaped reference basis values and orientation count."""
function mesh_basis_functions(element_type_value,local_coord,
                              function_space_type,wanted_orientations=Int32[];
                              caller::AbstractString="mesh_basis_functions")
    element_type=_linear_simplex_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    coordinates,point_count=_checked_local_coordinates(local_coord,caller)
    total_orientations=space.hierarchical ?
        _orientation_count(element_type) : 1
    orientations=_checked_orientation_sequence(
        wanted_orientations,total_orientations,space.hierarchical,caller)
    function_count=space.key_dimension==0 ?
        _vertex_count(element_type) : _edge_count(element_type)
    result_length=_checked_result_length(
        length(orientations),point_count,function_count,space.components;
        caller=caller)
    result=Vector{Float64}(undef,result_length)
    gradients=_reference_gradients(element_type)
    pairs=_canonical_edge_pairs(element_type)
    cursor=0
    @inbounds for orientation in orientations
        permutation=space.hierarchical ?
            _orientation_permutation(_vertex_count(element_type),orientation) :
            Int[]
        for point in 1:point_count
            base=3point-2
            lambda=_barycentric_coordinates(
                element_type,coordinates[base],coordinates[base+1],
                coordinates[base+2],caller,point)
            if space.kind==:lagrange
                for value in lambda
                    cursor+=1
                    result[cursor]=value
                end
            elseif space.kind==:grad_lagrange
                for gradient in gradients,component in 1:3
                    cursor+=1
                    result[cursor]=gradient[component]
                end
            else
                for (edge,(first_vertex,second_vertex)) in enumerate(pairs)
                    sign=permutation[first_vertex]<permutation[second_vertex] ?
                        1.0 : -1.0
                    for component in 1:3
                        cursor+=1
                        result[cursor]=sign*(space.kind==:hcurl ?
                            _scaled_whitney_component(
                                element_type,first_vertex,
                                lambda[first_vertex],
                                gradients[first_vertex][component],
                                second_vertex,
                                lambda[second_vertex],
                                gradients[second_vertex][component],
                                coordinates[base],coordinates[base+1],
                                coordinates[base+2],
                                caller,point,edge,component) :
                            _curl_component(
                                gradients[first_vertex],
                                gradients[second_vertex],component))
                    end
                end
            end
        end
    end
    return Int32(space.components),result,Int32(total_orientations)
end

"""Return the number of reference orientations for a supported function space."""
function mesh_number_of_orientations(element_type_value,function_space_type;
                                     caller::AbstractString=
                                         "mesh_number_of_orientations")
    element_type=_linear_simplex_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    return Int32(space.hierarchical ? _orientation_count(element_type) : 1)
end

function _orientation_rank(node_tags)
    count=length(node_tags)
    rank=0
    @inbounds for first in 1:count-1
        smaller=0
        for second in first+1:count
            smaller+=node_tags[second]<node_tags[first]
        end
        rank+=smaller*factorial(count-first)
    end
    return Int32(rank)
end

function _cell_orientation(cells::AbstractMatrix{Int32},cell::Int,
                           space::_FunctionSpace)
    space.hierarchical || return Int32(0)
    return _orientation_rank(@view cells[:,cell])
end

"""Return one orientation index per cached element of a linear-simplex type."""
function mesh_basis_orientations(mesh::Mesh,element_type_value,
                                 function_space_type;
                                 caller::AbstractString=
                                     "mesh_basis_orientations")
    element_type=_linear_simplex_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    block=mesh_element_block(mesh,element_type)
    block===nothing && return Int32[]
    _,cells=block
    result=Vector{Int32}(undef,size(cells,2))
    @inbounds for cell in axes(cells,2)
        result[cell]=_cell_orientation(cells,cell,space)
    end
    return result
end

"""Return the orientation index for one dense cached element tag."""
function mesh_basis_orientation(mesh::Mesh,element_tag_value,
                                function_space_type;
                                caller::AbstractString=
                                    "mesh_basis_orientation")
    tag=_checked_element_tag(mesh,element_tag_value,caller)
    space=_function_space(function_space_type,caller)
    record=mesh_element_record(mesh,tag)
    _linear_simplex_type(record.element_type,caller)
    return space.hierarchical ? _orientation_rank(record.node_tags) : Int32(0)
end

"""Return the node (0) or edge (1) dimension that owns keys for a space."""
function mesh_key_dimension(element_type_value,function_space_type;
                            caller::AbstractString="mesh_key_dimension")
    _linear_simplex_type(element_type_value,caller)
    return _function_space(function_space_type,caller).key_dimension
end

"""Return the number of keys per supported linear-simplex element."""
function mesh_number_of_keys(element_type_value,function_space_type;
                             caller::AbstractString="mesh_number_of_keys")
    element_type=_linear_simplex_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    return Int32(space.key_dimension==0 ?
                 _vertex_count(element_type) : _edge_count(element_type))
end

function _checked_bool(value,caller::AbstractString,name::AbstractString)
    value isa Bool || throw(ArgumentError("$caller: $name must be Bool"))
    return value
end

function _coordinate_result_length(key_count::Int,return_coord::Bool,
                                   caller::AbstractString)
    return_coord || return 0
    return _checked_result_length(3,key_count;caller=caller)
end

function _node_keys(mesh::Mesh,cells::AbstractMatrix{Int32},
                    return_coord::Bool,caller::AbstractString)
    key_count=length(cells)
    type_keys=fill(Int32(0),key_count)
    entity_keys=Vector{UInt64}(undef,key_count)
    coordinates=Vector{Float64}(
        undef,_coordinate_result_length(key_count,return_coord,caller))
    cursor=0
    @inbounds for cell in axes(cells,2),local_node in axes(cells,1)
        cursor+=1
        node=Int(cells[local_node,cell])
        entity_keys[cursor]=UInt64(node)
        if return_coord
            offset=3cursor-3
            coordinates[offset+1]=mesh.coords[1,node]
            coordinates[offset+2]=mesh.coords[2,node]
            coordinates[offset+3]=mesh.coords[3,node]
        end
    end
    return type_keys,entity_keys,coordinates
end

function _edge_keys(mesh::Mesh,cells::AbstractMatrix{Int32},
                    element_type::Int,topology::Union{Nothing,MeshEdgeTopology},
                    return_coord::Bool,caller::AbstractString)
    topology===nothing && throw(ArgumentError(
        "$caller: edge keys require a populated mesh edge catalog"))
    patterns=_simplex_edge_patterns(element_type)
    key_count=_checked_result_length(
        size(cells,2),length(patterns);caller=caller)
    node_pairs=Vector{UInt64}(undef,2key_count)
    cursor=0
    @inbounds for cell in axes(cells,2),pattern in patterns
        cursor+=1
        node_pairs[2cursor-1]=UInt64(cells[pattern[1],cell])
        node_pairs[2cursor]=UInt64(cells[pattern[2],cell])
    end
    entity_keys,_=_mesh_edges(topology,mesh,node_pairs,caller)
    type_keys=fill(Int32(1),key_count)
    coordinates=Vector{Float64}(
        undef,_coordinate_result_length(key_count,return_coord,caller))
    if return_coord
        @inbounds for key in 1:key_count
            first_node=Int(node_pairs[2key-1])
            second_node=Int(node_pairs[2key])
            offset=3key-3
            for component in 1:3
                coordinates[offset+component]=_stable_midpoint(
                    mesh.coords[component,first_node],
                    mesh.coords[component,second_node])
            end
        end
    end
    return type_keys,entity_keys,coordinates
end

function _keys_for_cells(mesh::Mesh,cells::AbstractMatrix{Int32},
                         element_type::Int,space::_FunctionSpace,
                         topology::Union{Nothing,MeshEdgeTopology},
                         return_coord::Bool,caller::AbstractString)
    return space.key_dimension==0 ?
        _node_keys(mesh,cells,return_coord,caller) :
        _edge_keys(
            mesh,cells,element_type,topology,return_coord,caller)
end

"""Return detached keys for all cached elements of one supported type."""
function mesh_keys(mesh::Mesh,element_type_value,function_space_type,
                   topology::Union{Nothing,MeshEdgeTopology}=nothing;
                   return_coord=true,caller::AbstractString="mesh_keys")
    element_type=_linear_simplex_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    coordinates_requested=_checked_bool(
        return_coord,caller,"return_coord")
    block=mesh_element_block(mesh,element_type)
    block===nothing && return Int32[],UInt64[],Float64[]
    _,cells=block
    return _keys_for_cells(
        mesh,cells,element_type,space,topology,coordinates_requested,caller)
end

"""Return detached keys for one dense cached element tag."""
function mesh_keys_for_element(
    mesh::Mesh,element_tag_value,function_space_type,
    topology::Union{Nothing,MeshEdgeTopology}=nothing;
    return_coord=true,caller::AbstractString="mesh_keys_for_element")
    tag=_checked_element_tag(mesh,element_tag_value,caller)
    space=_function_space(function_space_type,caller)
    coordinates_requested=_checked_bool(
        return_coord,caller,"return_coord")
    record=mesh_element_record(mesh,tag)
    element_type=_linear_simplex_type(record.element_type,caller)
    cells=reshape(Int32.(record.node_tags),length(record.node_tags),1)
    return _keys_for_cells(
        mesh,cells,element_type,space,topology,coordinates_requested,caller)
end

function _checked_type_keys(values,expected::Int,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: type_keys must be a vector or tuple of integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    result=Vector{Int32}(undef,length(values))
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: type_keys[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: type_keys[$index] must not be Bool"))
        converted=try Int(value) catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: type_keys[$index] exceeds Int bounds"))
        end
        converted==expected || throw(ArgumentError(
            "$caller: type_keys[$index] must be $expected for this function space"))
        result[index]=Int32(converted)
    end
    return result
end

function _checked_entity_keys(values,caller::AbstractString)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$caller: entity_keys must be a vector or tuple of positive integers"))
    values isa AbstractArray && Base.require_one_based_indexing(values)
    result=Vector{UInt64}(undef,length(values))
    for (index,value) in enumerate(values)
        value isa Integer || throw(ArgumentError(
            "$caller: entity_keys[$index] must be an integer"))
        value isa Bool && throw(ArgumentError(
            "$caller: entity_keys[$index] must not be Bool"))
        converted=try UInt64(value) catch err
            err isa InterruptException && rethrow()
            (err isa InexactError || err isa OverflowError ||
             err isa MethodError) || rethrow()
            throw(ArgumentError(
                "$caller: entity_keys[$index] must fit UInt64"))
        end
        converted>0 || throw(ArgumentError(
            "$caller: entity_keys[$index] must be positive"))
        result[index]=converted
    end
    return result
end

"""Return `(entity dimension, polynomial order)` metadata for complete key groups."""
function mesh_keys_information(type_keys,entity_keys,element_type_value,
                               function_space_type;
                               caller::AbstractString="mesh_keys_information")
    element_type=_linear_simplex_type(element_type_value,caller)
    space=_function_space(function_space_type,caller)
    expected_type=space.key_dimension
    types=_checked_type_keys(type_keys,expected_type,caller)
    entities=_checked_entity_keys(entity_keys,caller)
    length(types)==length(entities) || throw(ArgumentError(
        "$caller: type_keys and entity_keys must have equal lengths"))
    keys_per_element=space.key_dimension==0 ?
        _vertex_count(element_type) : _edge_count(element_type)
    length(types)%keys_per_element==0 || throw(ArgumentError(
        "$caller: key count $(length(types)) must be divisible by " *
        "$keys_per_element for element type $element_type"))
    return fill(
        (Int32(space.key_dimension),Int32(space.key_order)),length(types))
end

end # module MeshFunctionSpaces
