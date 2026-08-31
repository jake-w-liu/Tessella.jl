"""
    MeshPointLocation

Deterministic point location and reference-coordinate inversion for finalized
linear-simplex [`Mesh`](@ref) values. A reusable [`SimplexLocator`](@ref) owns a
compact AABB hierarchy over dense segment, triangle, and tetrahedron tags. The
geometric predicates use scaled `Float64` filters with exact-rational fallbacks for
degenerate or ill-conditioned affine maps.
"""
module MeshPointLocation

using ..MeshTypes: Mesh, node, nsegs, ntris, ntets, triangle_area, tet_volume

export SimplexLocator, mesh_element_offsets, mesh_element_block
export mesh_element_record, mesh_local_coordinates, locate_elements

const STRICT_REFERENCE_TOLERANCE=1.0e-6
const RELAXED_REFERENCE_TOLERANCES=(1.0e-5,1.0e-4,1.0e-3,1.0e-2,
                                    1.0e-1,1.0)
const _LOCATOR_LEAF_SIZE=8
# Below sqrt(eps), solving through a floating determinant can lose enough digits to
# change a 1e-6 reference-tolerance decision even when the determinant's sign is
# reliable. Route those affine maps through the exact-rational path.
const _FILTER_FACTOR=sqrt(eps(Float64))
const _Exact=Rational{BigInt}

"""Return the triangle offset, tetrahedron offset, and total dense element count."""
function mesh_element_offsets(mesh::Mesh)
    triangle_offset=nsegs(mesh)
    tetrahedron_offset=Base.checked_add(triangle_offset,ntris(mesh))
    total=Base.checked_add(tetrahedron_offset,ntets(mesh))
    return triangle_offset,tetrahedron_offset,total
end

"""Return `(dense_tag_offset, connectivity)` for integer MSH type 1, 2, or 4."""
function mesh_element_block(mesh::Mesh,element_type::Integer)
    element_type isa Bool && throw(ArgumentError(
        "mesh_element_block: element_type must not be Bool"))
    triangle_offset,tetrahedron_offset,_=mesh_element_offsets(mesh)
    element_type==1 && return 0,mesh.segs
    element_type==2 && return triangle_offset,mesh.tris
    element_type==4 && return tetrahedron_offset,mesh.tets
    return nothing
end

@inline function _element_dimension(mesh::Mesh,tag::Int)
    triangle_offset,tetrahedron_offset,_=mesh_element_offsets(mesh)
    tag<=triangle_offset && return 1
    tag<=tetrahedron_offset && return 2
    return 3
end

@inline function _element_reference(mesh::Mesh,tag::Int)
    triangle_offset,tetrahedron_offset,_=mesh_element_offsets(mesh)
    if tag<=triangle_offset
        return Int32(1),1,mesh.segs,tag
    elseif tag<=tetrahedron_offset
        return Int32(2),2,mesh.tris,tag-triangle_offset
    end
    return Int32(4),3,mesh.tets,tag-tetrahedron_offset
end

function _checked_element_tag(mesh::Mesh,value,caller::AbstractString)
    value isa Integer || throw(ArgumentError(
        "$caller: element_tag must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: element_tag must not be Bool"))
    tag=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError || err isa MethodError) ||
            rethrow()
        throw(ArgumentError(
            "$caller: element_tag exceeds the platform Int range"))
    end
    _,_,total=mesh_element_offsets(mesh)
    1<=tag<=total || throw(ArgumentError(
        "$caller: unknown element tag $value; expected a dense tag in 1:$total"))
    return tag
end

function _checked_coordinate(value,caller::AbstractString,name::AbstractString)
    value isa Real || throw(ArgumentError(
        "$caller: $name must be a real number"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    converted=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError || err isa MethodError) ||
            rethrow()
        throw(ArgumentError(
            "$caller: $name must be Float64-representable"))
    end
    isfinite(converted) || throw(ArgumentError(
        "$caller: $name must be finite"))
    return converted
end

function _checked_point(x,y,z,caller::AbstractString)
    return (_checked_coordinate(x,caller,"x"),
            _checked_coordinate(y,caller,"y"),
            _checked_coordinate(z,caller,"z"))
end

function _checked_dimension(value,caller::AbstractString)
    value isa Integer || throw(ArgumentError("$caller: dim must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: dim must not be Bool"))
    dimension=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        (err isa InexactError || err isa OverflowError || err isa MethodError) ||
            rethrow()
        throw(ArgumentError("$caller: dim exceeds the platform Int range"))
    end
    dimension in -1:3 || throw(ArgumentError(
        "$caller: dim must be -1 or in 0:3"))
    return dimension
end

"""
    mesh_element_record(mesh, element_tag)

Return detached type, dimension, and node-tag data for one dense simplex element.
The finalized `Mesh` does not own a model-entity classification tag.
"""
function mesh_element_record(mesh::Mesh,element_tag)
    caller="mesh_element_record"
    tag=_checked_element_tag(mesh,element_tag,caller)
    element_type,dimension,cells,cell=_element_reference(mesh,tag)
    node_tags=Vector{UInt64}(undef,size(cells,1))
    @inbounds for local_node in axes(cells,1)
        node_tags[local_node]=UInt64(cells[local_node,cell])
    end
    return (element_type=element_type,dimension=dimension,node_tags=node_tags)
end

@inline _sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _dot3(a,b)=muladd(a[1],b[1],muladd(a[2],b[2],a[3]*b[3]))
@inline _norm3(a)=hypot(a[1],a[2],a[3])
@inline _finite3(a)=isfinite(a[1])&&isfinite(a[2])&&isfinite(a[3])
@inline _maxabs3(a)=max(abs(a[1]),abs(a[2]),abs(a[3]))
@inline _scale3(a,s)=(a[1]*s,a[2]*s,a[3]*s)

@inline function _det3(a,b,c)
    return a[1]*(b[2]*c[3]-b[3]*c[2])-
           a[2]*(b[1]*c[3]-b[3]*c[1])+
           a[3]*(b[1]*c[2]-b[2]*c[1])
end

function _scaled_vectors(a,p,edges...)
    raw_edges=map(edge->_sub3(edge,a),edges)
    all(_finite3,raw_edges) || return nothing
    scale=maximum(_maxabs3,raw_edges)
    scale>0 || return nothing
    reciprocal=inv(scale)
    scaled_edges=map(edge->_scale3(edge,reciprocal),raw_edges)
    all(_finite3,scaled_edges) || return nothing
    rhs=_sub3(p,a)
    _finite3(rhs) || return nothing
    scaled_rhs=_scale3(rhs,reciprocal)
    _finite3(scaled_rhs) || return nothing
    return scaled_rhs,scaled_edges
end

@inline _exact3(a)=(_Exact(a[1]),_Exact(a[2]),_Exact(a[3]))
@inline _exact_sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _exact_dot3(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

@inline _exact_float(value)=Float64(value)

function _exact_residual(residual_squared,scale_squared,
                         caller::AbstractString,tag::Int)
    residual_squared==0 && return 0.0
    scale_squared>0 || throw(ArgumentError(
        "$caller: element $tag is degenerate"))
    result=setprecision(BigFloat,256) do
        sqrt(BigFloat(residual_squared)/BigFloat(scale_squared))
    end
    # An infinite residual is a valid internal no-match result. It must not abort a
    # relaxed scan just because another extremely small element is far from `p`.
    return Float64(result)
end

function _segment_coordinates_exact(a,b,p,caller::AbstractString,tag::Int)
    ae=_exact3(a);be=_exact3(b);pe=_exact3(p)
    edge=_exact_sub3(be,ae);rhs=_exact_sub3(pe,ae)
    denominator=_exact_dot3(edge,edge)
    denominator>0 || throw(ArgumentError(
        "$caller: element $tag is a degenerate segment"))
    parameter=_exact_dot3(rhs,edge)/denominator
    residual=ntuple(i->rhs[i]-parameter*edge[i],3)
    residual_squared=_exact_dot3(residual,residual)
    u=_exact_float(2parameter-1)
    return (u,0.0,0.0),
           _exact_residual(residual_squared,denominator,caller,tag)
end

function _triangle_coordinates_exact(a,b,c,p,caller::AbstractString,tag::Int)
    ae=_exact3(a);be=_exact3(b);ce=_exact3(c);pe=_exact3(p)
    edge1=_exact_sub3(be,ae);edge2=_exact_sub3(ce,ae)
    rhs=_exact_sub3(pe,ae)
    gram11=_exact_dot3(edge1,edge1)
    gram12=_exact_dot3(edge1,edge2)
    gram22=_exact_dot3(edge2,edge2)
    determinant=gram11*gram22-gram12*gram12
    determinant>0 || throw(ArgumentError(
        "$caller: element $tag is a degenerate triangle"))
    right1=_exact_dot3(rhs,edge1);right2=_exact_dot3(rhs,edge2)
    u=(right1*gram22-right2*gram12)/determinant
    v=(right2*gram11-right1*gram12)/determinant
    residual=ntuple(i->rhs[i]-u*edge1[i]-v*edge2[i],3)
    edge3=_exact_sub3(edge2,edge1)
    scale_squared=max(gram11,gram22,_exact_dot3(edge3,edge3))
    residual_squared=_exact_dot3(residual,residual)
    return (_exact_float(u),_exact_float(v),0.0),
           _exact_residual(residual_squared,scale_squared,caller,tag)
end

function _tetrahedron_coordinates_exact(a,b,c,d,p,
                                        caller::AbstractString,tag::Int)
    ae=_exact3(a);be=_exact3(b);ce=_exact3(c);de=_exact3(d);pe=_exact3(p)
    edge1=_exact_sub3(be,ae);edge2=_exact_sub3(ce,ae)
    edge3=_exact_sub3(de,ae);rhs=_exact_sub3(pe,ae)
    determinant=_det3(edge1,edge2,edge3)
    determinant!=0 || throw(ArgumentError(
        "$caller: element $tag is a degenerate tetrahedron"))
    u=_det3(rhs,edge2,edge3)/determinant
    v=_det3(edge1,rhs,edge3)/determinant
    w=_det3(edge1,edge2,rhs)/determinant
    return (_exact_float(u),_exact_float(v),_exact_float(w)),0.0
end

function _segment_coordinates(a,b,p,caller::AbstractString,tag::Int)
    scaled=_scaled_vectors(a,p,b)
    scaled===nothing && return _segment_coordinates_exact(a,b,p,caller,tag)
    rhs,edges=scaled;edge=only(edges)
    denominator=_dot3(edge,edge)
    denominator>0 || return _segment_coordinates_exact(a,b,p,caller,tag)
    parameter=_dot3(rhs,edge)/denominator
    residual=ntuple(i->rhs[i]-parameter*edge[i],3)
    u=muladd(2.0,parameter,-1.0)
    normalized_residual=_norm3(residual)/sqrt(denominator)
    all(isfinite,(u,normalized_residual)) ||
        return _segment_coordinates_exact(a,b,p,caller,tag)
    return (u,0.0,0.0),normalized_residual
end

function _triangle_coordinates(a,b,c,p,caller::AbstractString,tag::Int)
    scaled=_scaled_vectors(a,p,b,c)
    scaled===nothing &&
        return _triangle_coordinates_exact(a,b,c,p,caller,tag)
    rhs,edges=scaled;edge1,edge2=edges
    gram11=_dot3(edge1,edge1);gram12=_dot3(edge1,edge2)
    gram22=_dot3(edge2,edge2)
    product=gram11*gram22;square=gram12*gram12
    determinant=product-square
    permanent=abs(product)+abs(square)
    (isfinite(determinant)&&determinant>_FILTER_FACTOR*permanent) ||
        return _triangle_coordinates_exact(a,b,c,p,caller,tag)
    right1=_dot3(rhs,edge1);right2=_dot3(rhs,edge2)
    u=(right1*gram22-right2*gram12)/determinant
    v=(right2*gram11-right1*gram12)/determinant
    residual=ntuple(i->rhs[i]-u*edge1[i]-v*edge2[i],3)
    edge3=_sub3(edge2,edge1)
    scale=max(sqrt(gram11),sqrt(gram22),_norm3(edge3))
    normalized_residual=_norm3(residual)/scale
    all(isfinite,(u,v,normalized_residual)) ||
        return _triangle_coordinates_exact(a,b,c,p,caller,tag)
    return (u,v,0.0),normalized_residual
end

function _tetrahedron_coordinates(a,b,c,d,p,
                                  caller::AbstractString,tag::Int)
    scaled=_scaled_vectors(a,p,b,c,d)
    scaled===nothing &&
        return _tetrahedron_coordinates_exact(a,b,c,d,p,caller,tag)
    rhs,edges=scaled;edge1,edge2,edge3=edges
    determinant=_det3(edge1,edge2,edge3)
    permanent=
        abs(edge1[1])*(abs(edge2[2]*edge3[3])+abs(edge2[3]*edge3[2]))+
        abs(edge1[2])*(abs(edge2[1]*edge3[3])+abs(edge2[3]*edge3[1]))+
        abs(edge1[3])*(abs(edge2[1]*edge3[2])+abs(edge2[2]*edge3[1]))
    (isfinite(determinant)&&abs(determinant)>_FILTER_FACTOR*permanent) ||
        return _tetrahedron_coordinates_exact(a,b,c,d,p,caller,tag)
    u=_det3(rhs,edge2,edge3)/determinant
    v=_det3(edge1,rhs,edge3)/determinant
    w=_det3(edge1,edge2,rhs)/determinant
    all(isfinite,(u,v,w)) ||
        return _tetrahedron_coordinates_exact(a,b,c,d,p,caller,tag)
    return (u,v,w),0.0
end

function _local_coordinates(mesh::Mesh,tag::Int,p::NTuple{3,Float64},
                            caller::AbstractString)
    _,dimension,cells,cell=_element_reference(mesh,tag)
    @inbounds a=node(mesh,cells[1,cell])
    if dimension==1
        @inbounds b=node(mesh,cells[2,cell])
        coordinates,residual=_segment_coordinates(a,b,p,caller,tag)
    elseif dimension==2
        @inbounds b=node(mesh,cells[2,cell]);c=node(mesh,cells[3,cell])
        coordinates,residual=_triangle_coordinates(a,b,c,p,caller,tag)
    else
        @inbounds b=node(mesh,cells[2,cell]);c=node(mesh,cells[3,cell])
        @inbounds d=node(mesh,cells[4,cell])
        coordinates,residual=
            _tetrahedron_coordinates(a,b,c,d,p,caller,tag)
    end
    return coordinates,residual,dimension
end

function _require_local_coordinates(coordinates,caller::AbstractString,tag::Int)
    all(isfinite,coordinates) || throw(ArgumentError(
        "$caller: local coordinates for element $tag are not " *
        "Float64-representable"))
    return coordinates
end

"""
    mesh_local_coordinates(mesh, element_tag, x, y, z) -> (u, v, w)

Return affine reference coordinates for a dense linear-simplex element. Segment and
triangle inputs are orthogonally projected onto their affine span; their unused
reference coordinates are exactly zero. Degenerate and nonrepresentable inversions
fail explicitly.
"""
function mesh_local_coordinates(mesh::Mesh,element_tag,x,y,z)
    caller="mesh_local_coordinates"
    tag=_checked_element_tag(mesh,element_tag,caller)
    p=_checked_point(x,y,z,caller)
    coordinates,_,_=_local_coordinates(mesh,tag,p,caller)
    return _require_local_coordinates(coordinates,caller,tag)
end

@inline function _safe_edge_length(a,b)
    difference=_sub3(a,b)
    _finite3(difference) || return Inf
    return _norm3(difference)
end

@inline function _stable_midpoint(lower::Float64,upper::Float64)
    return signbit(lower)==signbit(upper) ?
           lower+(upper-lower)/2 : lower/2+upper/2
end

function _finish_element_bounds(lo::NTuple{3,Float64},
                                hi::NTuple{3,Float64},
                                maximum_edge::Float64,dimension::Int)
    coordinate_scale=max(abs(lo[1]),abs(lo[2]),abs(lo[3]),
                         abs(hi[1]),abs(hi[2]),abs(hi[3]))
    roundoff_padding=8eps(coordinate_scale)
    reference_padding=(dimension+2)*STRICT_REFERENCE_TOLERANCE*maximum_edge
    padding=max(roundoff_padding,reference_padding)
    padded_lo=(lo[1]-padding,lo[2]-padding,lo[3]-padding)
    padded_hi=(hi[1]+padding,hi[2]+padding,hi[3]+padding)
    centroid=(_stable_midpoint(lo[1],hi[1]),
              _stable_midpoint(lo[2],hi[2]),
              _stable_midpoint(lo[3],hi[3]))
    return padded_lo,padded_hi,centroid
end

function _segment_bounds(mesh::Mesh,cell::Int,tag::Int)
    @inbounds a=node(mesh,mesh.segs[1,cell]);b=node(mesh,mesh.segs[2,cell])
    a!=b || throw(ArgumentError(
        "SimplexLocator: element $tag is a degenerate segment"))
    lo=(min(a[1],b[1]),min(a[2],b[2]),min(a[3],b[3]))
    hi=(max(a[1],b[1]),max(a[2],b[2]),max(a[3],b[3]))
    return _finish_element_bounds(lo,hi,_safe_edge_length(a,b),1)
end

function _triangle_bounds(mesh::Mesh,cell::Int,tag::Int)
    @inbounds a=node(mesh,mesh.tris[1,cell]);b=node(mesh,mesh.tris[2,cell])
    @inbounds c=node(mesh,mesh.tris[3,cell])
    triangle_area(a,b,c)>0 || throw(ArgumentError(
        "SimplexLocator: element $tag is a degenerate triangle"))
    lo=(min(a[1],b[1],c[1]),min(a[2],b[2],c[2]),min(a[3],b[3],c[3]))
    hi=(max(a[1],b[1],c[1]),max(a[2],b[2],c[2]),max(a[3],b[3],c[3]))
    maximum_edge=max(_safe_edge_length(a,b),_safe_edge_length(a,c),
                     _safe_edge_length(b,c))
    return _finish_element_bounds(lo,hi,maximum_edge,2)
end

function _tetrahedron_bounds(mesh::Mesh,cell::Int,tag::Int)
    @inbounds a=node(mesh,mesh.tets[1,cell]);b=node(mesh,mesh.tets[2,cell])
    @inbounds c=node(mesh,mesh.tets[3,cell]);d=node(mesh,mesh.tets[4,cell])
    tet_volume(a,b,c,d)>0 || throw(ArgumentError(
        "SimplexLocator: element $tag is a degenerate tetrahedron"))
    lo=(min(a[1],b[1],c[1],d[1]),min(a[2],b[2],c[2],d[2]),
        min(a[3],b[3],c[3],d[3]))
    hi=(max(a[1],b[1],c[1],d[1]),max(a[2],b[2],c[2],d[2]),
        max(a[3],b[3],c[3],d[3]))
    maximum_edge=max(
        _safe_edge_length(a,b),_safe_edge_length(a,c),
        _safe_edge_length(a,d),_safe_edge_length(b,c),
        _safe_edge_length(b,d),_safe_edge_length(c,d))
    return _finish_element_bounds(lo,hi,maximum_edge,3)
end

"""
    SimplexLocator(mesh)

Build a deterministic AABB hierarchy for all dense linear-simplex elements in
`mesh`. Construction rejects any degenerate cell, so subsequent searches cannot
silently skip invalid geometry. The finalized mesh geometry must remain unchanged
for the lifetime of the locator.
"""
struct SimplexLocator
    mesh::Mesh
    order::Vector{Int}
    lower::Vector{NTuple{3,Float64}}
    upper::Vector{NTuple{3,Float64}}
    left::Vector{Int}
    right::Vector{Int}
    first::Vector{Int}
    count::Vector{Int}

    function SimplexLocator(mesh::Mesh)
        _,_,total=mesh_element_offsets(mesh)
        order=collect(1:total)
        primitive_lower=Vector{NTuple{3,Float64}}(undef,total)
        primitive_upper=similar(primitive_lower)
        centroids=similar(primitive_lower)
        tag=0
        @inbounds for cell in axes(mesh.segs,2)
            tag+=1
            primitive_lower[tag],primitive_upper[tag],centroids[tag]=
                _segment_bounds(mesh,cell,tag)
        end
        @inbounds for cell in axes(mesh.tris,2)
            tag+=1
            primitive_lower[tag],primitive_upper[tag],centroids[tag]=
                _triangle_bounds(mesh,cell,tag)
        end
        @inbounds for cell in axes(mesh.tets,2)
            tag+=1
            primitive_lower[tag],primitive_upper[tag],centroids[tag]=
                _tetrahedron_bounds(mesh,cell,tag)
        end
        lower=NTuple{3,Float64}[];upper=NTuple{3,Float64}[]
        left=Int[];right=Int[];firsts=Int[];counts=Int[]
        total==0 && return new(
            mesh,order,lower,upper,left,right,firsts,counts)

        function range_bounds(first::Int,last::Int)
            lo1=Inf;lo2=Inf;lo3=Inf
            hi1=-Inf;hi2=-Inf;hi3=-Inf
            @inbounds for position in first:last
                tag=order[position]
                primitive_lo=primitive_lower[tag]
                primitive_hi=primitive_upper[tag]
                lo1=min(lo1,primitive_lo[1]);lo2=min(lo2,primitive_lo[2])
                lo3=min(lo3,primitive_lo[3]);hi1=max(hi1,primitive_hi[1])
                hi2=max(hi2,primitive_hi[2]);hi3=max(hi3,primitive_hi[3])
            end
            return (lo1,lo2,lo3),(hi1,hi2,hi3)
        end

        function build!(first::Int,last::Int)
            tree_node=length(left)+1
            push!(lower,(0.0,0.0,0.0));push!(upper,(0.0,0.0,0.0))
            push!(left,0);push!(right,0);push!(firsts,0);push!(counts,0)
            range_count=last-first+1
            if range_count<=_LOCATOR_LEAF_SIZE
                lower[tree_node],upper[tree_node]=range_bounds(first,last)
                firsts[tree_node]=first;counts[tree_node]=range_count
                return tree_node
            end
            centroid_lo1=Inf;centroid_lo2=Inf;centroid_lo3=Inf
            centroid_hi1=-Inf;centroid_hi2=-Inf;centroid_hi3=-Inf
            @inbounds for position in first:last
                center=centroids[order[position]]
                centroid_lo1=min(centroid_lo1,center[1])
                centroid_lo2=min(centroid_lo2,center[2])
                centroid_lo3=min(centroid_lo3,center[3])
                centroid_hi1=max(centroid_hi1,center[1])
                centroid_hi2=max(centroid_hi2,center[2])
                centroid_hi3=max(centroid_hi3,center[3])
            end
            extent1=centroid_hi1-centroid_lo1
            extent2=centroid_hi2-centroid_lo2
            extent3=centroid_hi3-centroid_lo3
            axis=extent1>=extent2 ?
                 (extent1>=extent3 ? 1 : 3) :
                 (extent2>=extent3 ? 2 : 3)
            sort!(view(order,first:last);alg=QuickSort,
                  lt=(left_tag,right_tag)->begin
                      left_value=centroids[left_tag][axis]
                      right_value=centroids[right_tag][axis]
                      left_value<right_value ||
                          (left_value==right_value && left_tag<right_tag)
                  end)
            middle=(first+last)>>>1
            left_node=build!(first,middle)
            right_node=build!(middle+1,last)
            left[tree_node]=left_node;right[tree_node]=right_node
            left_lower=lower[left_node];right_lower=lower[right_node]
            left_upper=upper[left_node];right_upper=upper[right_node]
            lower[tree_node]=(min(left_lower[1],right_lower[1]),
                              min(left_lower[2],right_lower[2]),
                              min(left_lower[3],right_lower[3]))
            upper[tree_node]=(max(left_upper[1],right_upper[1]),
                              max(left_upper[2],right_upper[2]),
                              max(left_upper[3],right_upper[3]))
            return tree_node
        end

        build!(1,total)==1 || error("SimplexLocator: invalid hierarchy root")
        return new(mesh,order,lower,upper,left,right,firsts,counts)
    end
end

@inline function _bounds_contain(lower,upper,p)
    return lower[1]<=p[1]<=upper[1] &&
           lower[2]<=p[2]<=upper[2] &&
           lower[3]<=p[3]<=upper[3]
end

function _candidate_tags!(tags::Vector{Int},locator::SimplexLocator,
                          tree_node::Int,p::NTuple{3,Float64})
    tree_node==0 && return tags
    _bounds_contain(locator.lower[tree_node],locator.upper[tree_node],p) ||
        return tags
    count=locator.count[tree_node]
    if count>0
        first=locator.first[tree_node]
        append!(tags,@view locator.order[first:first+count-1])
        return tags
    end
    _candidate_tags!(tags,locator,locator.left[tree_node],p)
    _candidate_tags!(tags,locator,locator.right[tree_node],p)
    return tags
end

@inline function _inside_reference(coordinates,residual::Float64,
                                   dimension::Int,tolerance::Float64)
    u,v,w=coordinates
    if dimension==1
        return -(1+tolerance)<=u<=1+tolerance && residual<=tolerance
    elseif dimension==2
        return u>=-tolerance && v>=-tolerance &&
               u+v<=1+tolerance && residual<=tolerance
    end
    return u>=-tolerance && v>=-tolerance && w>=-tolerance &&
           u+v+w<=1+tolerance
end

@inline function _dimension_matches(mesh::Mesh,tag::Int,dimension::Int)
    return dimension<0 || _element_dimension(mesh,tag)==dimension
end

function _matching_tags!(matches::Vector{Int},candidates,
                         locator::SimplexLocator,p::NTuple{3,Float64},
                         dimension::Int,tolerance::Float64,
                         caller::AbstractString)
    empty!(matches)
    for tag in candidates
        _dimension_matches(locator.mesh,tag,dimension) || continue
        coordinates,residual,element_dimension=
            _local_coordinates(locator.mesh,tag,p,caller)
        _inside_reference(
            coordinates,residual,element_dimension,tolerance) && push!(matches,tag)
    end
    return matches
end

function _sort_matches!(matches::Vector{Int},mesh::Mesh)
    sort!(matches;by=tag->(-_element_dimension(mesh,tag),tag))
    return matches
end

function _locate_elements(locator::SimplexLocator,p::NTuple{3,Float64},
                          dimension::Int,strict::Bool,
                          caller::AbstractString)
    isempty(locator.order) && return UInt64[]
    candidates=Int[]
    _candidate_tags!(candidates,locator,1,p)
    matches=Int[]
    _matching_tags!(matches,candidates,locator,p,dimension,
                    STRICT_REFERENCE_TOLERANCE,caller)
    if isempty(matches) && !strict
        for tolerance in RELAXED_REFERENCE_TOLERANCES
            _matching_tags!(matches,eachindex(locator.order),locator,p,
                            dimension,tolerance,caller)
            isempty(matches) || break
        end
    end
    _sort_matches!(matches,locator.mesh)
    return UInt64.(matches)
end

"""
    locate_elements(locator, x, y, z; dim=-1, strict=false)

Return detached dense tags of all cached simplex elements containing or near a
point. Results are ordered by decreasing dimension and then increasing tag. Strict
search uses Gmsh 4.15.2's default `1e-6` reference tolerance. Relaxed search widens
that tolerance by decades through `1.0`, stopping at the first nonempty level.
"""
function locate_elements(locator::SimplexLocator,x,y,z;dim=-1,strict=false)
    caller="locate_elements"
    p=_checked_point(x,y,z,caller)
    dimension=_checked_dimension(dim,caller)
    strict isa Bool || throw(ArgumentError("$caller: strict must be Bool"))
    return _locate_elements(locator,p,dimension,strict,caller)
end

end # module
