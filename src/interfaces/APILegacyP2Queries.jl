# Queries on the retained simplex overlay use its actual interpolation nodes.
# The linear cache continues to supply entity membership and public cell tags.
@inline _legacy_p2_cells(mesh::P2TriMesh) = mesh.tri6
@inline _legacy_p2_cells(mesh::P2Mesh) = mesh.tet10
@inline _legacy_p2_type(::P2TriMesh) = 9
@inline _legacy_p2_type(::P2Mesh) = 11
@inline function _cache_node_tag(mesh::Union{P2TriMesh,P2Mesh},index::Integer)
    cached=LAST_MESH[]
    return _high_order_overlay(cached)===mesh ?
        _cache_node_tag(_mesh_public_tags(cached),index) : UInt64(index)
end

function _legacy_p2_selected(mesh,selected,caller)
    count=size(_legacy_p2_cells(mesh),2)
    if selected isa UnitRange{Int}
        MeshReferenceGeometry._checked_element_range(selected,count,caller)
        return selected
    end
    return MeshReferenceGeometry._checked_element_positions(selected,count,caller)
end

function _legacy_p2_basis_orientations(mesh,selected,function_space_type;
        caller::AbstractString="API.mesh.get_basis_functions_orientation")
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    _,family,_,_=MeshFunctionSpaces._basis_element_contract(
        _legacy_p2_type(mesh),space,caller)
    positions=_legacy_p2_selected(mesh,selected,caller)
    cells=_legacy_p2_cells(mesh)
    return Int32[_mixed_cell_orientation(mesh,cells,cell,space,family)
                 for cell in positions]
end

function _legacy_p2_keys(mesh,selected,offset,function_space_type,
        topology=nothing,face_topology=nothing;
        return_coord=true,caller::AbstractString="API.mesh.get_keys")
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    msh,_,nodal_count,_=MeshFunctionSpaces._basis_element_contract(
        _legacy_p2_type(mesh),space,caller)
    coords=MeshFunctionSpaces._checked_bool(return_coord,caller,"return_coord")
    positions=_legacy_p2_selected(mesh,selected,caller)
    cells=@view _legacy_p2_cells(mesh)[:,positions]
    cached=LAST_MESH[]
    public=_high_order_overlay(cached)===mesh ? _mesh_public_tags(cached) : nothing
    tags=_cache_element_tags(public,offset .+ positions)
    return _mixed_keys_for_cells(mesh,cells,msh,space,nodal_count,
        topology,face_topology,tags,coords,caller)
end

function _legacy_p2_barycenters(mesh,selected,fast::Bool,primary::Bool,
        caller::AbstractString)
    positions=_legacy_p2_selected(mesh,selected,caller)
    cells=_legacy_p2_cells(mesh)
    count=primary ? msh_properties(_legacy_p2_type(mesh)).num_primary_nodes : size(cells,1)
    return _mesh_barycenters(mesh,@view(cells[1:count,positions]),fast,caller)
end

function _legacy_p2_pattern_nodes(mesh,selected,patterns,caller::AbstractString)
    positions=_legacy_p2_selected(mesh,selected,caller)
    return _mesh_pattern_nodes(@view(_legacy_p2_cells(mesh)[:,positions]),patterns)
end

# Simplex P2 derivatives are affine. Keep their values at the reference corners
# in tuples so repeated quality samples need neither nodal-basis arrays nor a
# temporary mesh. Translate and scale before forming products.
function _legacy_p2_quality_vertices(mesh,cell::Int,caller::AbstractString)
    cells=_legacy_p2_cells(mesh)
    1<=cell<=size(cells,2) || throw(ArgumentError("$caller: unknown quadratic cell $cell"))
    origin=_mixed_node(mesh,cells,cell,1)
    scale=0.0
    for local_node in axes(cells,1),axis in 1:3
        delta=mesh.coords[axis,cells[local_node,cell]]-origin[axis]
        isfinite(delta) || throw(ArgumentError("$caller: quadratic cell $cell has unrepresentable coordinate differences"))
        scale=max(scale,abs(delta))
    end
    scale>0 || return nothing,0.0
    element_scale=scale
    count=mesh isa P2TriMesh ? 3 : 4
    vertices=ntuple(count) do vertex
        q=vertex==1 ? (0.0,0.0,0.0) :
            vertex==2 ? (1.0,0.0,0.0) :
            vertex==3 ? (0.0,1.0,0.0) : (0.0,0.0,1.0)
        gradients=if mesh isa P2TriMesh
            dr,ds=HighOrder._p2tri_grads(q[1],q[2])
            ntuple(i->(dr[i],ds[i],0.0),6)
        else
            HighOrder._p2_grads(q...)
        end
        ntuple(3) do column
            ntuple(3) do axis
                total=0.0
                for local_node in 2:length(gradients)
                    delta=(mesh.coords[axis,cells[local_node,cell]]-origin[axis])/element_scale
                    total=muladd(gradients[local_node][column],delta,total)
                end
                total
            end
        end
    end
    return vertices,scale
end

@inline function _legacy_p2_quality_frame(vertices,lambda)
    return ntuple(column->ntuple(axis->sum(
        lambda[vertex]*vertices[vertex][column][axis] for vertex in eachindex(vertices)),3),3)
end

function _legacy_p2_primary_normal(mesh::P2TriMesh,cell,scale)
    cells=mesh.tri6
    origin=_mixed_node(mesh,cells,cell,1)
    edges=ntuple(vertex->ntuple(axis->
        (mesh.coords[axis,cells[vertex+1,cell]]-origin[axis])/scale,3),2)
    normal=MeshElementQuality._cross3(edges...)
    length=hypot(normal...)
    return length>0 ? normal ./ length : (0.0,0.0,0.0),length
end

const _LEGACY_P2_TRI_DET_GROUP=ntuple(i->ntuple(j->
    findfirst(==((count(==(1),(i,j)),count(==(2),(i,j)),count(==(3),(i,j)))),
              HighOrder._B2TRI_MULTI),3),3)

function _legacy_p2_det_coefficients(vertices::NTuple{3},normal)
    coefficients=zeros(Float64,6)
    for first in 1:3,second in 1:3
        group=_LEGACY_P2_TRI_DET_GROUP[first][second]
        coefficients[group]+=MeshElementQuality._dot3(normal,
            MeshElementQuality._cross3(vertices[first][1],vertices[second][2]))
    end
    for index in eachindex(coefficients)
        coefficients[index]/=HighOrder._B2TRI_WEIGHT[index]
    end
    return coefficients
end

function _legacy_p2_det_coefficients(vertices::NTuple{4},normal=nothing)
    coefficients=zeros(Float64,20)
    for first in 1:4,second in 1:4,third in 1:4
        group=HighOrder._B3_GROUP[first,second,third]
        coefficients[group]+=MeshElementQuality._det3(
            vertices[first][1],vertices[second][2],vertices[third][3])
    end
    for index in eachindex(coefficients)
        coefficients[index]/=HighOrder._B3_NPERM[index]
    end
    return coefficients
end

function _legacy_p2_scaled_quantity(value,scale,power)
    factor=scale^power
    result=value*factor
    if isfinite(result) && (result!=0 || value==0)
        return result
    end
    return setprecision(BigFloat,128) do
        MeshElementQuality._float_preserve(BigFloat(value)*BigFloat(scale)^power)
    end
end

function _legacy_p2_sampled_shape(mesh,vertices,normal,name)
    points=mesh isa P2TriMesh ?
        ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),
         (0.5,0.5,0.0),(0.0,0.5,0.5),(0.5,0.0,0.5)) :
        ((1.0,0.0,0.0,0.0),(0.0,1.0,0.0,0.0),(0.0,0.0,1.0,0.0),(0.0,0.0,0.0,1.0),
         (0.5,0.5,0.0,0.0),(0.0,0.5,0.5,0.0),(0.5,0.0,0.5,0.0),
         (0.5,0.0,0.0,0.5),(0.0,0.0,0.5,0.5),(0.0,0.5,0.0,0.5))
    minimum_value=Inf
    for lambda in points
        du,dv,dw=_legacy_p2_quality_frame(vertices,lambda)
        value=if mesh isa P2TriMesh
            shape=MeshElementQuality._triangle_quality(name,(0.0,0.0,0.0),du,dv)
            cross=MeshElementQuality._cross3(du,dv)
            signed=MeshElementQuality._dot3(normal,cross)
            # SICN uses singular values and only its sign uses the primary
            # normal. SIGE uses the projected signed determinant itself.
            if name=="minSIGE"
                magnitude=hypot(cross...)
                magnitude>0 ? shape*(signed/magnitude) : 0.0
            else
                copysign(shape,signed)
            end
        else
            MeshElementQuality._tetrahedron_quality(name,(0.0,0.0,0.0),du,dv,dw)
        end
        minimum_value=min(minimum_value,value)
    end
    return minimum_value
end

const _LEGACY_P2_TRI_AREA_RULE=let
    points,weights=mesh_integration_points(2,"Gauss5")
    Tuple(((1-points[3point-2]-points[3point-1],points[3point-2],points[3point-1]),
           weights[point]) for point in eachindex(weights))
end

function _legacy_p2_triangle_area(vertices,scale)
    # MElement::getVolume uses polynomial-order + 3 integration points.
    area=0.0
    for (lambda,weight) in _LEGACY_P2_TRI_AREA_RULE
        du,dv,_=_legacy_p2_quality_frame(vertices,lambda)
        area+=weight*hypot(MeshElementQuality._cross3(du,dv)...)
    end
    return _legacy_p2_scaled_quantity(area,scale,2)
end

function _legacy_p2_quality(mesh,cell::Int,quality_name,caller::AbstractString)
    name=MeshElementQuality._quality_name(quality_name,caller)
    cells=_legacy_p2_cells(mesh)
    1<=cell<=size(cells,2) || throw(ArgumentError("$caller: unknown quadratic cell $cell"))
    primary=ntuple(vertex->_mixed_node(mesh,cells,cell,vertex),mesh isa P2TriMesh ? 3 : 4)
    if !(name in ("minSICN","minSIGE","minSJ","minDetJac","maxDetJac","minIsotropy") ||
         (name=="volume" && mesh isa P2TriMesh))
        return mesh isa P2TriMesh ? MeshElementQuality._triangle_quality(name,primary...) :
            MeshElementQuality._tetrahedron_quality(name,primary...)
    end
    vertices,scale=_legacy_p2_quality_vertices(mesh,cell,caller)
    if vertices===nothing
        return mesh isa P2TriMesh ? MeshElementQuality._triangle_quality(name,primary...) :
            MeshElementQuality._tetrahedron_quality(name,primary...)
    end
    if mesh isa P2TriMesh
        normal,primary_determinant=_legacy_p2_primary_normal(mesh,cell,scale)
    else
        normal=nothing
        origin=primary[1]
        edges=ntuple(vertex->ntuple(axis->(primary[vertex+1][axis]-origin[axis])/scale,3),3)
        primary_determinant=abs(MeshElementQuality._det3(edges...))
    end
    name in ("minSICN","minSIGE") &&
        return _legacy_p2_sampled_shape(mesh,vertices,normal,name)
    name=="volume" && mesh isa P2TriMesh &&
        return _legacy_p2_triangle_area(vertices,scale)
    if name=="minSJ"
        coefficients=_legacy_p2_det_coefficients(vertices,normal)
        return primary_determinant>0 ? minimum(coefficients)/primary_determinant : 0.0
    end
    if name in ("minDetJac","maxDetJac","minIsotropy")
        return _legacy_p2_bounded_quality(mesh,vertices,normal,scale,name,caller)
    end
    return mesh isa P2TriMesh ? MeshElementQuality._triangle_quality(name,primary...) :
        MeshElementQuality._tetrahedron_quality(name,primary...)
end
