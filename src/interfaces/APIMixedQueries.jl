# Native reference maps for the linear and quadratic mixed-element cache. In
# particular, a hexahedron remains a hexahedron during inversion and evaluation.
import ..MeshPointLocation: mesh_element_record, _local_coordinates,
                            _locate_elements
import ..MeshReferenceGeometry: mesh_jacobian, mesh_jacobians
import ..MeshFunctionSpaces: mesh_basis_orientation, mesh_basis_orientations,
                             mesh_keys, mesh_keys_for_element
import ..MeshElementQuality: mesh_element_qualities, mesh_element_quality
using ..MeshPointLocation: MeshPointLocation
using ..MeshReferenceGeometry: MeshReferenceGeometry
using ..MeshEntityTopology: MeshEntityTopology
using ..MeshElementQuality: MeshElementQuality

function _mixed_reference_node_index(nodes,coordinate,caller)
    for index in axes(nodes,2)
        (nodes[1,index],nodes[2,index],nodes[3,index])==coordinate && return index
    end
    throw(ArgumentError("$caller: the element has no interpolation node at $coordinate"))
end

function _mesh_query_edge_patterns(msh::Int,primary::Bool,caller::AbstractString)
    spec=msh_spec(msh)
    first_order=spec.family===:pnt ? 15 : msh_type(spec.family,1)
    patterns=MeshEntityTopology._simplex_edge_patterns(first_order)
    (primary || spec.order<=1 || isempty(patterns)) && return patterns
    spec.order==2 || throw(ArgumentError(
        "$caller: full edge-node queries currently support linear and quadratic elements; " *
        "use primary=true for type-$msh order-$(spec.order) vertices"))
    reference=Elements.lagrange_nodes(msh)
    return Tuple((edge...,_mixed_reference_node_index(reference,
        ntuple(axis->(reference[axis,edge[1]]+reference[axis,edge[2]])/2,3),caller))
        for edge in patterns)
end

function _mesh_query_face_patterns(msh::Int,face_type::Int,primary::Bool,caller::AbstractString)
    spec=msh_spec(msh)
    first_order=spec.family===:pnt ? 15 : msh_type(spec.family,1)
    patterns=MeshEntityTopology._simplex_face_patterns(first_order,face_type)
    (primary || spec.order<=1 || isempty(patterns)) && return patterns
    spec.order==2 || throw(ArgumentError(
        "$caller: full face-node queries currently support linear and quadratic elements; " *
        "use primary=true for type-$msh order-$(spec.order) vertices"))
    reference=Elements.lagrange_nodes(msh)
    return Tuple(begin
        midpoints=ntuple(face_type) do edge
            first_node=face[edge];second_node=face[mod1(edge+1,face_type)]
            coordinate=ntuple(axis->(reference[axis,first_node]+reference[axis,second_node])/2,3)
            _mixed_reference_node_index(reference,coordinate,caller)
        end
        if face_type==4 && !spec.serendipity
            center=ntuple(axis->sum(reference[axis,vertex] for vertex in face)/4,3)
            (face...,midpoints...,_mixed_reference_node_index(reference,center,caller))
        else
            (face...,midpoints...)
        end
    end for face in patterns)
end

function _mixed_query_reference(mesh::MixedMesh,tag::Int,caller::AbstractString)
    offset=0
    for block in mesh.blocks
        block isa ElementBlock || throw(ArgumentError(
            "$caller: mixed reference queries require fixed-connectivity elements"))
        msh=block.msh;cells=block.nodes
        if offset<tag<=offset+size(cells,2)
            spec=msh_spec(msh)
            (spec.order in (1,2) || spec.family===:pnt) || throw(ArgumentError(
                "$caller: mixed reference queries require linear or quadratic elements; " *
                "type $msh has order $(spec.order)"))
            spec.family===:trih && throw(ArgumentError(
                "$caller: Trihedron has no reference mapping in Gmsh 4.15.2"))
            return Int(msh),spec.dim,cells,tag-offset,spec.family
        end
        offset+=size(cells,2)
    end
    throw(ArgumentError("$caller: unknown element tag $tag"))
end

function mesh_element_record(mesh::MixedMesh,element_tag)
    tag=_mesh_query_element_tag(mesh,element_tag,"mesh_element_record")
    msh,dim,cells,cell,_=_mixed_query_reference(mesh,tag,"mesh_element_record")
    return (element_type=Int32(msh),dimension=dim,
            node_tags=UInt64.(cells[:,cell]))
end

@inline _mixed_node(mesh::MixedMesh,cells,cell::Int,local_node::Int)=
    (mesh.coords[1,cells[local_node,cell]],
     mesh.coords[2,cells[local_node,cell]],
     mesh.coords[3,cells[local_node,cell]])

@inline _mixed_axpy3(value,delta,origin)=
    (muladd(value,delta[1],origin[1]),muladd(value,delta[2],origin[2]),
     muladd(value,delta[3],origin[3]))

function _mixed_map(mesh::MixedMesh,cells,cell::Int,family::Val,
                    uvw::NTuple{3,Float64},caller::AbstractString,tag::Int,msh::Int=0;
                    shape=nothing)
    u,v,w=uvw
    values,gradients=if shape!==nothing
        shape
    elseif msh==0 || msh_spec(msh).order<2
        (MeshFunctionSpaces._first_order_values(family,u,v,w,caller,tag),
         MeshFunctionSpaces._first_order_gradients(family,u,v,w,caller,tag))
    else
        local_coord=[u,v,w]
        _,nodal_values,_=mesh_basis_functions(msh,local_coord,"Lagrange";caller=caller)
        _,nodal_gradients,_=mesh_basis_functions(msh,local_coord,"GradLagrange";caller=caller)
        (nodal_values,[(nodal_gradients[3i-2],nodal_gradients[3i-1],
                        nodal_gradients[3i]) for i in eachindex(nodal_values)])
    end
    origin=_mixed_node(mesh,cells,cell,1)
    # Subtract one vertex before summation: translation must not contaminate
    # derivatives, and meshes far from the origin retain their local precision.
    mapped=origin
    du=(0.0,0.0,0.0);dv=du;dw=du
    @inbounds for local_node in 2:length(values)
        vertex=_mixed_node(mesh,cells,cell,local_node)
        delta=(vertex[1]-origin[1],vertex[2]-origin[2],vertex[3]-origin[3])
        value=values[local_node];gradient=gradients[local_node]
        mapped=_mixed_axpy3(value,delta,mapped)
        du=_mixed_axpy3(gradient[1],delta,du)
        dv=_mixed_axpy3(gradient[2],delta,dv)
        dw=_mixed_axpy3(gradient[3],delta,dw)
    end
    all(isfinite,(mapped...,du...,dv...,dw...)) || throw(ArgumentError(
        "$caller: reference mapping for element $tag is not Float64-representable"))
    return mapped,du,dv,dw
end

@inline function _mixed_initial_coordinates(family::Symbol)
    family===:tri && return (1/3,1/3,0.0)
    family===:tet && return (0.25,0.25,0.25)
    family===:pri && return (1/3,1/3,0.0)
    family===:pyr && return (0.0,0.0,0.25)
    return (0.0,0.0,0.0)
end

function _mixed_scale(mesh::MixedMesh,cells,cell::Int)
    origin=_mixed_node(mesh,cells,cell,1)
    scale=0.0
    @inbounds for local_node in 2:size(cells,1)
        vertex=_mixed_node(mesh,cells,cell,local_node)
        scale=max(scale,hypot(vertex[1]-origin[1],vertex[2]-origin[2],
                              vertex[3]-origin[3]))
    end
    return scale
end

function _mixed_newton_step(rhs,du,dv,dw,dimension::Int,
                            caller::AbstractString,tag::Int)
    dot=MeshPointLocation._dot3
    if dimension==1
        denominator=dot(du,du)
        isfinite(denominator) && denominator>0 || throw(ArgumentError(
            "$caller: element $tag has a degenerate reference mapping"))
        return (dot(du,rhs)/denominator,0.0,0.0)
    elseif dimension==2
        aa=dot(du,du);ab=dot(du,dv);bb=dot(dv,dv)
        denominator=aa*bb-ab*ab
        if !isfinite(denominator) || denominator<=sqrt(eps(Float64))*aa*bb
            coordinates,_=MeshPointLocation._triangle_coordinates_exact(
                (0.0,0.0,0.0),du,dv,rhs,caller,tag)
            return coordinates
        end
        ar=dot(du,rhs);br=dot(dv,rhs)
        return ((bb*ar-ab*br)/denominator,(aa*br-ab*ar)/denominator,0.0)
    end
    determinant=MeshPointLocation._det3(du,dv,dw)
    permanent=abs(du[1])*(abs(dv[2]*dw[3])+abs(dv[3]*dw[2]))+
              abs(du[2])*(abs(dv[1]*dw[3])+abs(dv[3]*dw[1]))+
              abs(du[3])*(abs(dv[1]*dw[2])+abs(dv[2]*dw[1]))
    if !isfinite(determinant) || abs(determinant)<=sqrt(eps(Float64))*permanent
        coordinates,_=MeshPointLocation._tetrahedron_coordinates_exact(
            (0.0,0.0,0.0),du,dv,dw,rhs,caller,tag)
        return coordinates
    end
    return (MeshPointLocation._det3(rhs,dv,dw)/determinant,
            MeshPointLocation._det3(du,rhs,dw)/determinant,
            MeshPointLocation._det3(du,dv,rhs)/determinant)
end

function _local_coordinates(mesh::MixedMesh,tag::Int,p::NTuple{3,Float64},
                             caller::AbstractString;location::Bool=false)
    msh,dimension,cells,cell,family=_mixed_query_reference(mesh,tag,caller)
    linear=msh_spec(msh).order<2
    a=_mixed_node(mesh,cells,cell,1)
    if family===:pnt
        return (0.0,0.0,0.0),hypot((p .- a)...),0
    elseif linear && family===:lin
        coordinates,residual=MeshPointLocation._segment_coordinates(
            a,_mixed_node(mesh,cells,cell,2),p,caller,tag)
        return coordinates,residual,1
    elseif linear && family===:tri
        coordinates,residual=MeshPointLocation._triangle_coordinates(
            a,_mixed_node(mesh,cells,cell,2),_mixed_node(mesh,cells,cell,3),
            p,caller,tag)
        return coordinates,residual,2
    elseif linear && family===:tet
        coordinates,residual=MeshPointLocation._tetrahedron_coordinates(
            a,_mixed_node(mesh,cells,cell,2),_mixed_node(mesh,cells,cell,3),
            _mixed_node(mesh,cells,cell,4),p,caller,tag)
        return coordinates,residual,3
    end
    scale=_mixed_scale(mesh,cells,cell)
    isfinite(scale) && scale>0 || throw(ArgumentError(
        "$caller: element $tag is degenerate or not Float64-representable"))
    # The pyramid tip has one physical point and one canonical reference point.
    if family===:pyr && p==_mixed_node(mesh,cells,cell,5)
        return (0.0,0.0,1.0),0.0,3
    end
    coordinates=_mixed_initial_coordinates(family)
    relative_scale=inv(scale)
    family_value=Val(family)
    for iteration in 1:64
        mapped,du,dv,dw=_mixed_map(
            mesh,cells,cell,family_value,coordinates,caller,tag,msh)
        rhs=((p[1]-mapped[1])*relative_scale,(p[2]-mapped[2])*relative_scale,
             (p[3]-mapped[3])*relative_scale)
        scaled_du=(du[1]*relative_scale,du[2]*relative_scale,du[3]*relative_scale)
        scaled_dv=(dv[1]*relative_scale,dv[2]*relative_scale,dv[3]*relative_scale)
        scaled_dw=(dw[1]*relative_scale,dw[2]*relative_scale,dw[3]*relative_scale)
        step=_mixed_newton_step(rhs,scaled_du,scaled_dv,scaled_dw,
                                 dimension,caller,tag)
        all(isfinite,step) || throw(ArgumentError(
            "$caller: local coordinates for element $tag are not Float64-representable"))
        updated=(coordinates[1]+step[1],coordinates[2]+step[2],
                 coordinates[3]+step[3])
        if maximum(abs,step)<=16eps(Float64)*max(1.0,maximum(abs,updated))
            final,_,_,_=_mixed_map(mesh,cells,cell,family_value,updated,caller,tag,msh)
            residual=hypot((p .- final)...)*relative_scale
            return updated,residual,dimension
        end
        coordinates=updated
    end
    location && return nothing,Inf,dimension
    throw(ArgumentError("$caller: reference inversion did not converge for element $tag"))
end

struct _MixedMeshLocator
    mesh::MixedMesh
    order::Vector{Int}
    lower::Vector{NTuple{3,Float64}}
    upper::Vector{NTuple{3,Float64}}
    left::Vector{Int}
    right::Vector{Int}
    first::Vector{Int}
    count::Vector{Int}
    dimensions::Vector{Int}
    families::Vector{Symbol}
end

function _MixedMeshLocator(mesh::MixedMesh)
    _,_,total=_mesh_element_offsets(mesh)
    order=collect(1:total)
    primitive_lower=Vector{NTuple{3,Float64}}(undef,total)
    primitive_upper=similar(primitive_lower);centroids=similar(primitive_lower)
    dimensions=Vector{Int}(undef,total);families=Vector{Symbol}(undef,total)
    for (msh,dim,offset,cells,_) in _cache_catalog(mesh),cell in axes(cells,2)
        tag=offset+cell
        _,_,_,_,family=_mixed_query_reference(mesh,tag,"mixed mesh locator")
        dimensions[tag]=dim;families[tag]=family
        lower=(Inf,Inf,Inf);upper=(-Inf,-Inf,-Inf)
        for vertex in axes(cells,1)
            point=_mixed_node(mesh,cells,cell,vertex)
            lower=(min(lower[1],point[1]),min(lower[2],point[2]),min(lower[3],point[3]))
            upper=(max(upper[1],point[1]),max(upper[2],point[2]),max(upper[3],point[3]))
        end
        scale=_mixed_scale(mesh,cells,cell)
        primitive_lower[tag],primitive_upper[tag],centroids[tag]=
            MeshPointLocation._finish_element_bounds(lower,upper,scale,dim)
        if msh_spec(msh).order>1
            # Curved Lagrange maps can leave the nodal convex hull. Without a
            # certified Bernstein enclosure, retain every high-order cell as
            # a candidate rather than silently excluding part of its geometry.
            primitive_lower[tag]=(-Inf,-Inf,-Inf)
            primitive_upper[tag]=(Inf,Inf,Inf)
        end
        # Check the center map before constructing the reusable locator.
        _mixed_frame(mesh,cells,cell,family,dim,
                      _mixed_initial_coordinates(family),"mixed mesh locator",tag,Int(msh))
    end
    lower=NTuple{3,Float64}[];upper=NTuple{3,Float64}[]
    left=Int[];right=Int[];firsts=Int[];counts=Int[]
    if total>0
        builder=MeshPointLocation._LocatorBuilder(order,lower,upper,left,right,
            firsts,counts,primitive_lower,primitive_upper,centroids)
        MeshPointLocation._locator_build!(builder,1,total)
    end
    return _MixedMeshLocator(mesh,order,lower,upper,left,right,firsts,counts,
                             dimensions,families)
end

const LAST_MIXED_MESH_LOCATOR=Ref{Union{Nothing,_MixedMeshLocator}}(nothing)

function _cached_mesh_locator_locked(mesh::MixedMesh)
    locator=LAST_MIXED_MESH_LOCATOR[]
    if locator===nothing || locator.mesh!==mesh
        locator=_MixedMeshLocator(mesh)
        LAST_MIXED_MESH_LOCATOR[]=locator
    end
    return locator
end

function _mixed_candidate_tags!(tags,locator::_MixedMeshLocator,tree_node::Int,p)
    tree_node==0 && return tags
    MeshPointLocation._bounds_contain(
        locator.lower[tree_node],locator.upper[tree_node],p) || return tags
    count=locator.count[tree_node]
    if count>0
        first=locator.first[tree_node]
        append!(tags,@view locator.order[first:first+count-1])
    else
        _mixed_candidate_tags!(tags,locator,locator.left[tree_node],p)
        _mixed_candidate_tags!(tags,locator,locator.right[tree_node],p)
    end
    return tags
end

@inline function _mixed_inside_reference(uvw,residual,family::Symbol,tolerance)
    u,v,w=uvw
    family===:pnt && return residual<=tolerance
    family===:lin && return abs(u)<=1+tolerance && residual<=tolerance
    family===:tri && return u>=-tolerance && v>=-tolerance &&
        u+v<=1+tolerance && residual<=tolerance
    family===:qua && return max(abs(u),abs(v))<=1+tolerance && residual<=tolerance
    family===:tet && return min(u,v,w)>=-tolerance && u+v+w<=1+tolerance
    family===:hex && return max(abs(u),abs(v),abs(w))<=1+tolerance && residual<=tolerance
    family===:pri && return u>=-tolerance && v>=-tolerance &&
        u+v<=1+tolerance && abs(w)<=1+tolerance && residual<=tolerance
    return w>=-tolerance && w<=1+tolerance &&
        max(abs(u),abs(v))<=1-w+tolerance && residual<=tolerance
end

function _mixed_matching_tags!(matches,candidates,locator,p,dimension,tolerance,caller)
    empty!(matches)
    for tag in candidates
        (dimension<0 || locator.dimensions[tag]==dimension) || continue
        uvw,residual,_=_local_coordinates(locator.mesh,tag,p,caller;location=true)
        # A curved candidate need not have an inverse outside its image. That
        # must not abort a query that belongs to another valid element.
        uvw===nothing && continue
        _mixed_inside_reference(uvw,residual,locator.families[tag],tolerance) &&
            push!(matches,tag)
    end
    return matches
end

function _locate_elements(locator::_MixedMeshLocator,p::NTuple{3,Float64},
                          dimension::Int,strict::Bool,caller::AbstractString)
    isempty(locator.order) && return UInt64[]
    candidates=Int[]
    _mixed_candidate_tags!(candidates,locator,1,p)
    matches=Int[]
    _mixed_matching_tags!(matches,candidates,locator,p,dimension,
        MeshPointLocation.STRICT_REFERENCE_TOLERANCE,caller)
    if isempty(matches) && !strict
        for tolerance in MeshPointLocation.RELAXED_REFERENCE_TOLERANCES
            _mixed_matching_tags!(matches,eachindex(locator.order),locator,p,
                                  dimension,tolerance,caller)
            isempty(matches) || break
        end
    end
    sort!(matches;by=tag->(-locator.dimensions[tag],tag))
    return UInt64.(matches)
end

function _mixed_frame(mesh::MixedMesh,cells,cell,family::Symbol,dimension::Int,
                      uvw,caller::AbstractString,tag::Int,msh::Int=0;shape=nothing)
    mapped,du,dv,dw=_mixed_map(mesh,cells,cell,Val(family),uvw,caller,tag,msh;shape=shape)
    if dimension==0
        return (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),1.0,mapped
    elseif dimension==1
        frame,determinant,_=MeshReferenceGeometry._segment_frame(
            (-du[1],-du[2],-du[3]),du,caller,tag)
        return frame,determinant,mapped
    elseif dimension==2
        normal=MeshReferenceGeometry._cross3(du,dv)
        determinant=hypot(normal...)
        isfinite(determinant) && determinant>0 || throw(ArgumentError(
            "$caller: element $tag has a degenerate surface mapping"))
        return (du...,dv...,(normal ./ determinant)...),determinant,mapped
    end
    determinant=MeshReferenceGeometry._tetrahedron_determinant(
        (0.0,0.0,0.0),du,dv,dw,du,dv,dw,caller,tag)
    return (du...,dv...,dw...),determinant,mapped
end

function _mixed_write_jacobians!(jacobians,determinants,coordinates,mesh,cells,
                                 selected,offset,msh,local_coordinates,point_count,
                                 caller)
    spec=msh_spec(msh)
    (spec.order in (1,2) || spec.family===:pnt) || throw(ArgumentError(
        "$caller: mixed reference queries require linear or quadratic elements"))
    spec.family===:trih && throw(ArgumentError(
        "$caller: Trihedron has no reference mapping in Gmsh 4.15.2"))
    # Reference shape data depends on type and local point, not the physical
    # cell. Evaluate the nodal engine once for a bulk query and reuse it.
    shapes=if spec.order==2
        _,values,_=mesh_basis_functions(msh,local_coordinates,"Lagrange";caller=caller)
        _,gradients,_=mesh_basis_functions(msh,local_coordinates,"GradLagrange";caller=caller)
        count=spec.nnodes
        [(@view(values[count*(point-1)+1:count*point]),
          [(gradients[3(count*(point-1)+local_node)-2],
            gradients[3(count*(point-1)+local_node)-1],
            gradients[3(count*(point-1)+local_node)]) for local_node in 1:count])
         for point in 1:point_count]
    else
        nothing
    end
    for (slot,cell) in enumerate(selected),point in 1:point_count
        base=3(point-1)
        uvw=(local_coordinates[base+1],local_coordinates[base+2],local_coordinates[base+3])
        frame,determinant,mapped=_mixed_frame(mesh,cells,cell,spec.family,
            spec.dim,uvw,caller,offset+cell,msh;
            shape=shapes===nothing ? nothing : shapes[point])
        evaluation=(slot-1)*point_count+point
        @inbounds for component in 1:9
            jacobians[9(evaluation-1)+component]=frame[component]
        end
        determinants[evaluation]=determinant
        @inbounds for component in 1:3
            coordinates[3(evaluation-1)+component]=mapped[component]
        end
    end
    return nothing
end

function mesh_jacobians(mesh::MixedMesh,element_type,local_coord,
                        positions::Union{Nothing,UnitRange{Int},AbstractVector{<:Integer}}=nothing)
    caller="mesh_jacobians"
    msh=MeshReferenceGeometry._checked_element_type(element_type,caller)
    local_coordinates,point_count=MeshReferenceGeometry._checked_local_coordinates(local_coord,caller)
    block=_mesh_element_block(mesh,msh)
    (block===nothing || point_count==0) && return Float64[],Float64[],Float64[]
    offset,cells=block
    selected=if positions===nothing
        1:size(cells,2)
    elseif positions isa UnitRange{Int}
        MeshReferenceGeometry._checked_element_range(positions,size(cells,2),caller)
        positions
    else
        MeshReferenceGeometry._checked_element_positions(positions,size(cells,2),caller)
    end
    isempty(selected) && return Float64[],Float64[],Float64[]
    results=MeshReferenceGeometry._allocate_results(length(selected),point_count,caller)
    _mixed_write_jacobians!(results...,mesh,cells,selected,offset,msh,
                            local_coordinates,point_count,caller)
    return results
end

function mesh_jacobian(mesh::MixedMesh,element_tag,local_coord)
    caller="mesh_jacobian"
    tag=_mesh_query_element_tag(mesh,element_tag,caller)
    local_coordinates,point_count=MeshReferenceGeometry._checked_local_coordinates(local_coord,caller)
    msh,_,cells,cell,_=_mixed_query_reference(mesh,tag,caller)
    point_count==0 && return Float64[],Float64[],Float64[]
    results=MeshReferenceGeometry._allocate_results(1,point_count,caller)
    _mixed_write_jacobians!(results...,mesh,cells,cell:cell,tag-cell,msh,
                            local_coordinates,point_count,caller)
    return results
end

function mesh_basis_orientations(mesh::MixedMesh,element_type,function_space_type,
        positions::Union{Nothing,UnitRange{Int},AbstractVector{<:Integer}}=nothing;
        caller::AbstractString="mesh_basis_orientations")
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    msh,family,_,_=MeshFunctionSpaces._basis_element_contract(element_type,space,caller)
    block=_mesh_element_block(mesh,msh)
    block===nothing && return Int32[]
    _,cells=block
    selected=if positions===nothing
        1:size(cells,2)
    elseif positions isa UnitRange{Int}
        MeshFunctionSpaces._checked_orientation_range(positions,size(cells,2),caller)
        positions
    else
        MeshReferenceGeometry._checked_element_positions(positions,size(cells,2),caller)
    end
    return Int32[_mixed_cell_orientation(mesh,cells,cell,space,family) for cell in selected]
end

@inline function _mixed_cell_orientation(mesh,cells,cell,space,family)
    space.hierarchical || return Int32(0)
    vertex_count=family===:pnt ? 1 : family===:lin ? 2 : family===:tri ? 3 :
        (family===:tet || family===:qua) ? 4 : family===:hex ? 8 :
        family===:pri ? 6 : 5
    # Hierarchical orientations are determined by public vertex tags. Sparse
    # tags can reverse a corner comparison even when dense indices increase.
    rank=0
    @inbounds for first in 1:vertex_count-1
        first_tag=_cache_node_tag(mesh,cells[first,cell])
        smaller=0
        for second in first+1:vertex_count
            smaller+=_cache_node_tag(mesh,cells[second,cell])<first_tag
        end
        rank+=smaller*factorial(vertex_count-first)
    end
    return Int32(rank)
end

function mesh_basis_orientation(mesh::MixedMesh,element_tag,function_space_type;
                                caller::AbstractString="mesh_basis_orientation")
    tag=_mesh_query_element_tag(mesh,element_tag,caller)
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    msh,_,cells,cell,family=_mixed_query_reference(mesh,tag,caller)
    MeshFunctionSpaces._basis_element_contract(msh,space,caller)
    return _mixed_cell_orientation(mesh,cells,cell,space,family)
end

function _mixed_keys_for_cells(mesh::MixedMesh,cells,msh::Int,space,
                               nodal_count,topology,face_topology,element_tags,
                               return_coord::Bool,caller::AbstractString)
    if !space.hierarchical && size(cells,1)!=nodal_count
        throw(ArgumentError("$caller: the cached element stores $(size(cells,1)) " *
            "nodes, but the requested basis requires $nodal_count"))
    end
    family=msh_spec(msh).family
    if !space.hierarchical
        counts=(vertex=size(cells,1),edge=0,face=0,bubble=0)
        type_pattern=zeros(Int32,size(cells,1))
        per_edge=0;face_blocks=Tuple{Int,Int,Symbol}[]
    else
        order=space.key_order
        hcurl=space.key_dimension==1
        counts=hcurl ? MeshFunctionSpaces._hcurl_counts(family,order) :
                       MeshFunctionSpaces._h1_counts(family,order)
        type_pattern=hcurl ? MeshFunctionSpaces._hcurl_type_keys(family,order) :
                             MeshFunctionSpaces._h1_type_keys(family,order)
        per_edge=hcurl ? order+1 : max(order-1,0)
        face_blocks=hcurl ? MeshFunctionSpaces._hcurl_face_blocks(family,order) :
                            MeshFunctionSpaces._h1_face_blocks(family,order)
    end
    function_count=length(type_pattern)
    key_count=MeshFunctionSpaces._checked_result_length(
        size(cells,2),function_count;caller=caller)
    type_keys=Vector{Int32}(undef,key_count)
    entity_keys=Vector{UInt64}(undef,key_count)
    coordinates=Vector{Float64}(undef,return_coord ? 3key_count : 0)
    primary=msh_spec(msh).family===:pnt ? 15 : msh_type(family,1)
    edge_patterns=MeshEntityTopology._simplex_edge_patterns(primary)
    # Gmsh's element keys follow its MElement face order. Prism face functions
    # pack the quadrangular faces first, followed by the triangular faces.
    face_patterns=family===:pri ?
        (MeshEntityTopology._simplex_face_patterns(primary,4)...,
         MeshEntityTopology._simplex_face_patterns(primary,3)...) :
        (MeshEntityTopology._simplex_face_patterns(primary,3)...,
         MeshEntityTopology._simplex_face_patterns(primary,4)...)
    cursor=0
    for cell in axes(cells,2)
        for local_node in 1:counts.vertex
            cursor+=1
            type_keys[cursor]=type_pattern[cursor-(cell-1)*function_count]
            node=Int(cells[local_node,cell]);entity_keys[cursor]=_cache_node_tag(mesh,node)
            if return_coord
                @inbounds for component in 1:3
                    coordinates[3(cursor-1)+component]=mesh.coords[component,node]
                end
            end
        end
        if counts.edge>0
            topology===nothing && throw(ArgumentError(
                "$caller: hierarchical keys require a populated mesh edge catalog"))
            for pattern in edge_patterns
                first_node=cells[pattern[1],cell];second_node=cells[pattern[2],cell]
                edge_tag=get(topology.tags,MeshEntityTopology._edge_key(first_node,second_node),nothing)
                edge_tag===nothing && throw(ArgumentError("$caller: unknown mesh edge"))
                for _ in 1:per_edge
                    cursor+=1
                    type_keys[cursor]=type_pattern[cursor-(cell-1)*function_count]
                    entity_keys[cursor]=edge_tag
                    if return_coord
                        @inbounds for component in 1:3
                            coordinates[3(cursor-1)+component]=
                                MeshReferenceGeometry._stable_midpoint(
                                    mesh.coords[component,first_node],mesh.coords[component,second_node])
                        end
                    end
                end
            end
        end
        if counts.face>0
            face_topology===nothing && throw(ArgumentError(
                "$caller: hierarchical keys require a populated mesh face catalog"))
            for (face,pattern) in enumerate(face_patterns)
                face_count=face_blocks[face][2]
                face_count==0 && continue
                face_tag=if length(pattern)==3
                    key=MeshEntityTopology._triangle_key(
                        cells[pattern[1],cell],cells[pattern[2],cell],cells[pattern[3],cell])
                    get(face_topology.triangle_tags,key,nothing)
                else
                    key=MeshEntityTopology._quadrangle_key(cells[pattern[1],cell],
                        cells[pattern[2],cell],cells[pattern[3],cell],cells[pattern[4],cell])
                    get(face_topology.quadrangle_tags,key,nothing)
                end
                face_tag===nothing && throw(ArgumentError("$caller: unknown mesh face"))
                for _ in 1:face_count
                    cursor+=1
                    type_keys[cursor]=type_pattern[cursor-(cell-1)*function_count]
                    entity_keys[cursor]=face_tag
                    if return_coord
                        @inbounds for component in 1:3
                            total=0.0
                            for local_node in pattern
                                total+=mesh.coords[component,cells[local_node,cell]]
                            end
                            coordinates[3(cursor-1)+component]=total/length(pattern)
                        end
                    end
                end
            end
        end
        for _ in 1:counts.bubble
            cursor+=1
            type_keys[cursor]=type_pattern[cursor-(cell-1)*function_count]
            entity_keys[cursor]=element_tags[cell]
            if return_coord
                @inbounds for component in 1:3
                    total=0.0
                    for local_node in axes(cells,1)
                        total+=mesh.coords[component,cells[local_node,cell]]
                    end
                    coordinates[3(cursor-1)+component]=total/size(cells,1)
                end
            end
        end
    end
    cursor==key_count || error("$caller: internal mixed key count mismatch")
    return type_keys,entity_keys,coordinates
end

_mesh_face_topology_for_key_cells(mesh::Mesh,topology,cells,msh,space)=
    _mesh_face_topology_for_cells(mesh,topology,cells,msh)

function _mesh_face_topology_for_key_cells(mesh::MixedMesh,topology,cells,msh,function_space_type)
    space=MeshFunctionSpaces._function_space(function_space_type,"API.mesh key query")
    family=msh_spec(msh).family
    primary=msh_type(family,1)
    blocks=space.key_dimension==1 ?
        MeshFunctionSpaces._hcurl_face_blocks(family,space.key_order) :
        MeshFunctionSpaces._h1_face_blocks(family,space.key_order)
    patterns=family===:pri ?
        (MeshEntityTopology._simplex_face_patterns(primary,4)...,
         MeshEntityTopology._simplex_face_patterns(primary,3)...) :
        (MeshEntityTopology._simplex_face_patterns(primary,3)...,
         MeshEntityTopology._simplex_face_patterns(primary,4)...)
    selected=Tuple(patterns[index] for index in eachindex(patterns) if blocks[index][2]>0)
    replacement=MeshEntityTopology._face_topology_copy(topology,size(mesh.coords,2))
    MeshEntityTopology._append_generated_faces!(replacement,cells,selected,primary)
    return replacement
end

function mesh_keys(mesh::MixedMesh,element_type,function_space_type,
                   topology::Union{Nothing,MeshEdgeTopology}=nothing,
                   face_topology::Union{Nothing,MeshFaceTopology}=nothing;
                   return_coord=true,caller::AbstractString="mesh_keys")
    return _mixed_mesh_keys(mesh,element_type,function_space_type,nothing,
                            topology,face_topology,return_coord,caller)
end

function mesh_keys(mesh::MixedMesh,element_type,function_space_type,
                   positions::AbstractVector{<:Integer},
                   topology::Union{Nothing,MeshEdgeTopology}=nothing,
                   face_topology::Union{Nothing,MeshFaceTopology}=nothing;
                   return_coord=true,caller::AbstractString="mesh_keys")
    return _mixed_mesh_keys(mesh,element_type,function_space_type,positions,
                            topology,face_topology,return_coord,caller)
end

function _mixed_mesh_keys(mesh::MixedMesh,element_type,function_space_type,
                          positions,topology,face_topology,return_coord,caller)
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    msh,_,nodal_count,_=MeshFunctionSpaces._basis_element_contract(element_type,space,caller)
    coordinate_requested=MeshFunctionSpaces._checked_bool(return_coord,caller,"return_coord")
    block=_mesh_element_block(mesh,msh)
    block===nothing && return Int32[],UInt64[],Float64[]
    offset,cells=block
    selected=positions===nothing ? collect(axes(cells,2)) :
        MeshReferenceGeometry._checked_element_positions(positions,size(cells,2),caller)
    isempty(selected) && return Int32[],UInt64[],Float64[]
    selected_cells=positions===nothing ? cells : cells[:,selected]
    return _mixed_keys_for_cells(mesh,selected_cells,msh,space,nodal_count,
        topology,face_topology,_cache_element_tags(mesh,offset .+ selected),coordinate_requested,caller)
end

function mesh_keys_for_element(mesh::MixedMesh,element_tag,function_space_type,
                   topology::Union{Nothing,MeshEdgeTopology}=nothing,
                   face_topology::Union{Nothing,MeshFaceTopology}=nothing;
                   return_coord=true,caller::AbstractString="mesh_keys_for_element")
    tag=_mesh_query_element_tag(mesh,element_tag,caller)
    msh,_,cells,cell,_=_mixed_query_reference(mesh,tag,caller)
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    _,_,nodal_count,_=MeshFunctionSpaces._basis_element_contract(msh,space,caller)
    coordinate_requested=MeshFunctionSpaces._checked_bool(return_coord,caller,"return_coord")
    return _mixed_keys_for_cells(mesh,cells[:,cell:cell],msh,space,nodal_count,
        topology,face_topology,UInt64[_cache_element_tag(mesh,tag)],coordinate_requested,caller)
end

function mesh_element_qualities(mesh::MixedMesh,element_tags,quality_name="minSICN")
    caller="mesh_element_qualities"
    name=MeshElementQuality._quality_name(quality_name,caller)
    (element_tags isa AbstractVector || element_tags isa Tuple) || throw(ArgumentError(
        "$caller: element_tags must be a vector or tuple of integers"))
    element_tags isa AbstractArray && Base.require_one_based_indexing(element_tags)
    result=Vector{Float64}(undef,length(element_tags))
    for (index,value) in enumerate(element_tags)
        tag=_mesh_query_element_tag(mesh,value,caller)
        msh,_,cells,cell,family=_mixed_query_reference(mesh,tag,caller)
        a=_mixed_node(mesh,cells,cell,1)
        if name=="minEdge" || name=="maxEdge"
            primary=msh_type(family,1)
            edges=MeshEntityTopology._simplex_edge_patterns(primary)
            isempty(edges) && throw(ArgumentError("$caller: $name is unsupported for $family elements"))
            quality=name=="minEdge" ? Inf : 0.0
            for (first_node,second_node) in edges
                length,_=MeshElementQuality._segment_geometry(
                    _mixed_node(mesh,cells,cell,first_node),_mixed_node(mesh,cells,cell,second_node))
                quality=name=="minEdge" ? min(quality,length) : max(quality,length)
            end
            result[index]=quality
        elseif msh_spec(msh).order==1 && family===:lin
            name in MeshElementQuality._UNDEFINED_SEGMENT_QUALITIES && throw(ArgumentError(
                "$caller: Gmsh 4.15.2 does not define $name reliably for segment elements"))
            length,valid=MeshElementQuality._segment_geometry(a,_mixed_node(mesh,cells,cell,2))
            result[index]=MeshElementQuality._segment_quality(name,length,valid)
        elseif msh_spec(msh).order==1 && family===:tri
            result[index]=MeshElementQuality._triangle_quality(name,a,
                _mixed_node(mesh,cells,cell,2),_mixed_node(mesh,cells,cell,3))
        elseif msh_spec(msh).order==1 && family===:tet
            result[index]=MeshElementQuality._tetrahedron_quality(name,a,
                _mixed_node(mesh,cells,cell,2),_mixed_node(mesh,cells,cell,3),_mixed_node(mesh,cells,cell,4))
        else
            throw(ArgumentError("$caller: quality $(repr(name)) is unsupported for native " *
                "type-$msh $family elements; supported native measures are minEdge and maxEdge"))
        end
    end
    return result
end

mesh_element_quality(mesh::MixedMesh,element_tag,quality_name="minSICN")=
    only(mesh_element_qualities(mesh,(element_tag,),quality_name))
