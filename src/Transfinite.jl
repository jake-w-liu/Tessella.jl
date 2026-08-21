"""
    Transfinite

Validated four-sided planar transfinite patches. Boundary chains are supplied
already discretized; opposite chains must contain the same number of nodes.
The interior interpolation and triangle arrangements reproduce the documented
four-corner Gmsh transfinite-surface construction for an affine planar surface.
"""
module Transfinite

using ..MeshTypes: Mesh, boundary_edges, nnodes, nsegs, ntris, validate
using ..Predicates: orient2

export mesh_transfinite_patch

const _DEFAULT_MAX_NODES = 10_000_000
const _DEFAULT_MAX_TRIANGLES = 20_000_000
const _INT32_MAX = Int(typemax(Int32))
const _BOUNDARY_AUDIT_MULTIPLIER = 64
const _BOUNDARY_AUDIT_FLOOR = 4096

struct _PlaneFrame
    u::NTuple{3,Float64}
    v::NTuple{3,Float64}
    n::NTuple{3,Float64}
end

struct _SegmentBox
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    segment::Int
end

struct _BoundaryNode
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    segment::Int
    left::Int
    right::Int
    count::Int
end

@inline _sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _dot3(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _cross3(a,b)=(a[2]*b[3]-a[3]*b[2],
                      a[3]*b[1]-a[1]*b[3],
                      a[1]*b[2]-a[2]*b[1])
@inline _norm3(a)=hypot(a[1],a[2],a[3])
@inline _edge_key(a::Int32,b::Int32)=a<b ? (a,b) : (b,a)

function _checked_add(a::Int,b::Int,what::AbstractString)
    try
        return Base.checked_add(a,b)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError("mesh_transfinite_patch: $what count overflows Int"))
    end
end

function _checked_mul(a::Int,b::Int,what::AbstractString)
    try
        return Base.checked_mul(a,b)
    catch err
        err isa InterruptException && rethrow()
        err isa OverflowError || rethrow()
        throw(ArgumentError("mesh_transfinite_patch: $what count overflows Int"))
    end
end

function _limit(value::Integer,name::AbstractString)
    value isa Bool && throw(ArgumentError(
        "mesh_transfinite_patch: $name must not be Bool"))
    value>=0 || throw(ArgumentError(
        "mesh_transfinite_patch: $name must be non-negative"))
    value<=typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_patch: $name exceeds the Int32 topology limit"))
    return Int(value)
end

function _tag(value,name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "mesh_transfinite_patch: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "mesh_transfinite_patch: $name must not be Bool"))
    0<=value<=typemax(Int32) || throw(ArgumentError(
        "mesh_transfinite_patch: $name must lie in 0:$(typemax(Int32))"))
    return Int32(value)
end

function _arrangement(value)
    value isa Symbol || throw(ArgumentError(
        "mesh_transfinite_patch: arrangement must be a Symbol"))
    value in (:left,:right,:alternate_left,:alternate_right) ||
        throw(ArgumentError(
            "mesh_transfinite_patch: arrangement must be :left, :right, " *
            ":alternate_left, or :alternate_right"))
    return value
end

function _point3(raw,side::Int,index::Int)
    local count
    try
        count=length(raw)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_patch: side $side point $index is not indexable"))
    end
    count==3 || throw(ArgumentError(
        "mesh_transfinite_patch: side $side point $index must have exactly three coordinates"))
    point=try
        (Float64(raw[1]),Float64(raw[2]),Float64(raw[3]))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_patch: side $side point $index coordinates must " *
            "be Float64-representable: $(sprint(showerror,err))"))
    end
    all(isfinite,point) || throw(ArgumentError(
        "mesh_transfinite_patch: side $side point $index has a non-finite coordinate"))
    return point
end

function _convert_side(side::AbstractVector,number::Int)
    result=Vector{NTuple{3,Float64}}(undef,length(side))
    destination=1
    for raw in side
        destination<=length(result) || throw(ArgumentError(
            "mesh_transfinite_patch: side $number iteration produced more than " *
            "its declared length"))
        result[destination]=_point3(raw,number,destination)
        destination+=1
    end
    destination==length(result)+1 || throw(ArgumentError(
        "mesh_transfinite_patch: side $number iteration length changed during conversion"))
    return result
end

function _distance(a,b,description::AbstractString)
    dx=a[1]-b[1];dy=a[2]-b[2];dz=a[3]-b[3]
    (isfinite(dx)&&isfinite(dy)&&isfinite(dz)) || throw(ArgumentError(
        "mesh_transfinite_patch: $description coordinate span overflows Float64"))
    distance=hypot(dx,dy,dz)
    (isfinite(distance)&&distance>0) || throw(ArgumentError(
        "mesh_transfinite_patch: $description has zero or non-finite length"))
    return distance
end

function _validate_side_edges(side,number::Int)
    @inbounds for i in 1:length(side)-1
        _distance(side[i+1],side[i],"side $number segment $i")
    end
    return nothing
end

function _boundary_ring(sides)
    count=sum(length(side)-1 for side in sides;init=0)
    ring=Vector{NTuple{3,Float64}}(undef,count)
    cursor=0
    @inbounds for side in sides
        for i in 1:length(side)-1
            cursor+=1;ring[cursor]=side[i]
        end
    end
    cursor==count || throw(ErrorException(
        "mesh_transfinite_patch: internal boundary count invariant failed"))
    return ring
end

function _normalization(ring)
    origin=ring[1];scale=0.0
    @inbounds for (i,point) in pairs(ring)
        dx=point[1]-origin[1];dy=point[2]-origin[2];dz=point[3]-origin[3]
        (isfinite(dx)&&isfinite(dy)&&isfinite(dz)) || throw(ArgumentError(
            "mesh_transfinite_patch: boundary coordinate span overflows Float64 at node $i"))
        scale=max(scale,abs(dx),abs(dy),abs(dz))
    end
    (isfinite(scale)&&scale>0) || throw(ArgumentError(
        "mesh_transfinite_patch: boundary is geometrically degenerate"))
    return origin,scale
end

@inline function _normalize(point,origin,scale)
    ((point[1]-origin[1])/scale,
     (point[2]-origin[2])/scale,
     (point[3]-origin[3])/scale)
end

function _normalized_side(side,origin,scale)
    result=Vector{NTuple{3,Float64}}(undef,length(side))
    @inbounds for i in eachindex(side)
        point=_normalize(side[i],origin,scale)
        all(isfinite,point) || throw(ArgumentError(
            "mesh_transfinite_patch: normalized boundary coordinate is not finite"))
        result[i]=point
    end
    return result
end

function _plane_frame(ring,origin,scale)
    nx=0.0;ny=0.0;nz=0.0
    @inbounds for i in eachindex(ring)
        p=_normalize(ring[i],origin,scale)
        q=_normalize(ring[mod1(i+1,length(ring))],origin,scale)
        nx+=(p[2]-q[2])*(p[3]+q[3])
        ny+=(p[3]-q[3])*(p[1]+q[1])
        nz+=(p[1]-q[1])*(p[2]+q[2])
    end
    normal_length=hypot(nx,ny,nz)
    (isfinite(normal_length)&&normal_length>0) || throw(ArgumentError(
        "mesh_transfinite_patch: boundary has no representable plane normal"))
    normal=(nx/normal_length,ny/normal_length,nz/normal_length)
    reference=abs(normal[1])<=abs(normal[2]) ?
        (abs(normal[1])<=abs(normal[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
        (abs(normal[2])<=abs(normal[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    axis_u=_cross3(reference,normal);axis_length=_norm3(axis_u)
    (isfinite(axis_length)&&axis_length>0) || throw(ArgumentError(
        "mesh_transfinite_patch: could not construct an in-plane frame"))
    axis_u=(axis_u[1]/axis_length,axis_u[2]/axis_length,axis_u[3]/axis_length)
    axis_v=_cross3(normal,axis_u)
    frame=_PlaneFrame(axis_u,axis_v,normal)
    tolerance=256eps(Float64)
    @inbounds for (i,point) in pairs(ring)
        normalized=_normalize(point,origin,scale)
        distance=abs(_dot3(normalized,normal))
        (isfinite(distance)&&distance<=tolerance) || throw(ArgumentError(
            "mesh_transfinite_patch: boundary node $i is not coplanar " *
            "(normalized distance $distance exceeds $tolerance)"))
    end
    return frame
end

@inline _project(frame::_PlaneFrame,point)=
    (_dot3(point,frame.u),_dot3(point,frame.v))

@inline function _on_segment(a,b,p)
    orient2(a,b,p)==0 || return false
    return min(a[1],b[1])<=p[1]<=max(a[1],b[1]) &&
           min(a[2],b[2])<=p[2]<=max(a[2],b[2])
end

@inline function _segments_intersect(a,b,c,d)
    o1=orient2(a,b,c);o2=orient2(a,b,d)
    o3=orient2(c,d,a);o4=orient2(c,d,b)
    return (o1==0&&_on_segment(a,b,c)) || (o2==0&&_on_segment(a,b,d)) ||
           (o3==0&&_on_segment(c,d,a)) || (o4==0&&_on_segment(c,d,b)) ||
           (o1!=0&&o2!=0&&o3!=0&&o4!=0&&o1!=o2&&o3!=o4)
end

@inline _boundary_adjacent(i::Int,j::Int,n::Int)=
    j==i+1 || (i==1&&j==n)

function _adjacent_overlap(a,b,c,d)
    shared = a==c || a==d ? a : b==c || b==d ? b : nothing
    shared===nothing && return true
    for point in (a,b)
        point==shared || !_on_segment(c,d,point) || return true
    end
    for point in (c,d)
        point==shared || !_on_segment(a,b,point) || return true
    end
    return false
end

@inline _box_overlap(a,b)=a.xmin<=b.xmax&&b.xmin<=a.xmax&&
                                 a.ymin<=b.ymax&&b.ymin<=a.ymax

function _build_boundary_tree!(nodes,order,boxes,lo::Int,hi::Int)
    xmin=Inf;xmax=-Inf;ymin=Inf;ymax=-Inf
    @inbounds for position in lo:hi
        box=boxes[order[position]]
        xmin=min(xmin,box.xmin);xmax=max(xmax,box.xmax)
        ymin=min(ymin,box.ymin);ymax=max(ymax,box.ymax)
    end
    count=hi-lo+1;node_index=length(nodes)+1
    push!(nodes,_BoundaryNode(xmin,xmax,ymin,ymax,0,0,0,count))
    if lo==hi
        @inbounds segment=boxes[order[lo]].segment
        nodes[node_index]=_BoundaryNode(xmin,xmax,ymin,ymax,segment,0,0,1)
        return node_index
    end
    xspan=xmax-xmin;yspan=ymax-ymin
    if xspan>=yspan
        sort!(@view(order[lo:hi]);
              by=index->begin box=boxes[index];box.xmin/2+box.xmax/2 end,
              alg=QuickSort)
    else
        sort!(@view(order[lo:hi]);
              by=index->begin box=boxes[index];box.ymin/2+box.ymax/2 end,
              alg=QuickSort)
    end
    middle=lo+(hi-lo)÷2
    left=_build_boundary_tree!(nodes,order,boxes,lo,middle)
    right=_build_boundary_tree!(nodes,order,boxes,middle+1,hi)
    nodes[node_index]=_BoundaryNode(xmin,xmax,ymin,ymax,0,left,right,count)
    return node_index
end

function _audit_boundary_pair(points,i::Int,j::Int,count::Int)
    a=points[i];b=points[mod1(i+1,count)]
    c=points[j];d=points[mod1(j+1,count)]
    _segments_intersect(a,b,c,d) || return nothing
    adjacent=_boundary_adjacent(min(i,j),max(i,j),count)
    if !adjacent || _adjacent_overlap(a,b,c,d)
        throw(ArgumentError(
            "mesh_transfinite_patch: boundary segments $i and $j intersect"))
    end
    return nothing
end

function _validate_simple_boundary(points)
    count=length(points)
    boxes=Vector{_SegmentBox}(undef,count)
    @inbounds for i in 1:count
        a=points[i];b=points[mod1(i+1,count)]
        (isfinite(a[1])&&isfinite(a[2])) || throw(ArgumentError(
            "mesh_transfinite_patch: projected boundary node $i is not finite"))
        boxes[i]=_SegmentBox(min(a[1],b[1]),max(a[1],b[1]),
                             min(a[2],b[2]),max(a[2],b[2]),i)
    end
    order=collect(1:count);nodes=_BoundaryNode[]
    node_capacity=_checked_add(_checked_mul(2,count,"boundary audit node"),-1,
                               "boundary audit node")
    sizehint!(nodes,node_capacity)
    root=_build_boundary_tree!(nodes,order,boxes,1,count)
    limit=max(_BOUNDARY_AUDIT_FLOOR,
              _checked_mul(_BOUNDARY_AUDIT_MULTIPLIER,count,
                           "boundary-intersection audit"))
    stack=Tuple{Int,Int}[(root,root)]
    visits=0;candidates=0
    while !isempty(stack)
        first,second=pop!(stack);visits+=1
        visits<=limit || throw(ArgumentError(
            "mesh_transfinite_patch: boundary intersection audit exceeded " *
            "its bounded traversal limit $limit"))
        node1=nodes[first];node2=nodes[second]
        _box_overlap(node1,node2) || continue
        if first==second
            node1.segment!=0 && continue
            push!(stack,(node1.left,node1.left),(node1.left,node1.right),
                        (node1.right,node1.right))
        elseif node1.segment!=0&&node2.segment!=0
            candidates+=1
            candidates<=limit || throw(ArgumentError(
                "mesh_transfinite_patch: boundary intersection audit exceeded " *
                "its bounded candidate limit $limit"))
            _audit_boundary_pair(points,node1.segment,node2.segment,count)
        elseif node2.segment!=0 || (node1.segment==0&&node1.count>=node2.count)
            push!(stack,(node1.left,second),(node1.right,second))
        else
            push!(stack,(first,node2.left),(first,node2.right))
        end
    end
    return nothing
end

function _averaged_parameters(first,opposite,direction::AbstractString)
    length(first)==length(opposite) || throw(ErrorException(
        "mesh_transfinite_patch: internal opposite-side count invariant failed"))
    result=zeros(Float64,length(first));total=0.0
    @inbounds for i in 1:length(first)-1
        first_length=_distance(first[i+1],first[i],"$direction side interval $i")
        opposite_length=_distance(opposite[i+1],opposite[i],
                                  "$direction opposite-side interval $i")
        increment=0.5first_length+0.5opposite_length
        (isfinite(increment)&&increment>0) || throw(ArgumentError(
            "mesh_transfinite_patch: $direction averaged chord $i is not finite and positive"))
        total+=increment
        isfinite(total) || throw(ArgumentError(
            "mesh_transfinite_patch: $direction averaged chord sum overflows Float64"))
        result[i+1]=total
    end
    total>0 || throw(ArgumentError(
        "mesh_transfinite_patch: $direction averaged chord sum is zero"))
    @inbounds for i in 2:length(result)-1
        result[i]/=total
        (isfinite(result[i])&&result[i-1]<result[i]<1) || throw(ArgumentError(
            "mesh_transfinite_patch: $direction coordinates are not strictly increasing"))
    end
    result[end]=1.0
    return result
end

# Gmsh 4.15.2 `meshGFaceTransfinite.cpp` uses the average chord spacing of
# opposing sides for u/v, then the standard four-sided transfinite (Coons)
# interpolation. With c1 translated to zero, its bilinear corner correction has
# the compact form below.
@inline function _coons(left,right,bottom,top,c2,c3,c4,u,v)
    one_u=1-u;one_v=1-v
    ntuple(3) do coordinate
        one_u*left[coordinate]+u*right[coordinate]+
        one_v*bottom[coordinate]+v*top[coordinate]-
        (u*one_v*c2[coordinate]+u*v*c3[coordinate]+one_u*v*c4[coordinate])
    end
end

@inline function _physical_point(normalized,origin,scale)
    point=(origin[1]+scale*normalized[1],
           origin[2]+scale*normalized[2],
           origin[3]+scale*normalized[3])
    all(isfinite,point) || throw(ArgumentError(
        "mesh_transfinite_patch: generated coordinate is not finite"))
    return point
end

@inline _node(i::Int,j::Int,width::Int)=Int32(i+1+j*width)

@inline function _right_diagonal(arrangement::Symbol,i::Int,j::Int)
    # Exact zero-based parity from Gmsh 4.15.2: AlternateRight selects the
    # v1-v3 diagonal on odd i+j; AlternateLeft selects it on even i+j.
    arrangement===:right && return true
    arrangement===:alternate_right && return isodd(i+j)
    arrangement===:alternate_left && return iseven(i+j)
    return false
end

function _fill_segments!(segments,tags,width::Int,L::Int,H::Int,side_tags)
    cursor=0
    @inbounds for i in 0:L-1
        cursor+=1;segments[1,cursor]=_node(i,0,width)
        segments[2,cursor]=_node(i+1,0,width);tags[cursor]=side_tags[1]
    end
    @inbounds for j in 0:H-1
        cursor+=1;segments[1,cursor]=_node(L,j,width)
        segments[2,cursor]=_node(L,j+1,width);tags[cursor]=side_tags[2]
    end
    @inbounds for i in L:-1:1
        cursor+=1;segments[1,cursor]=_node(i,H,width)
        segments[2,cursor]=_node(i-1,H,width);tags[cursor]=side_tags[3]
    end
    @inbounds for j in H:-1:1
        cursor+=1;segments[1,cursor]=_node(0,j,width)
        segments[2,cursor]=_node(0,j-1,width);tags[cursor]=side_tags[4]
    end
    cursor==size(segments,2) || throw(ErrorException(
        "mesh_transfinite_patch: internal segment count invariant failed"))
    return nothing
end

function _fill_triangles!(triangles,width::Int,L::Int,H::Int,arrangement::Symbol)
    cursor=0
    @inbounds for i in 0:L-1,j in 0:H-1
        v1=_node(i,j,width);v2=_node(i+1,j,width)
        v3=_node(i+1,j+1,width);v4=_node(i,j+1,width)
        if _right_diagonal(arrangement,i,j)
            cursor+=1;triangles[1,cursor]=v1;triangles[2,cursor]=v2
            triangles[3,cursor]=v3
            cursor+=1;triangles[1,cursor]=v3;triangles[2,cursor]=v4
            triangles[3,cursor]=v1
        else
            cursor+=1;triangles[1,cursor]=v1;triangles[2,cursor]=v2
            triangles[3,cursor]=v4
            cursor+=1;triangles[1,cursor]=v4;triangles[2,cursor]=v2
            triangles[3,cursor]=v3
        end
    end
    cursor==size(triangles,2) || throw(ErrorException(
        "mesh_transfinite_patch: internal triangle count invariant failed"))
    return nothing
end

@inline function _project_output(coords,node_id::Int32,origin,scale,frame)
    node=Int(node_id)
    normalized=((coords[1,node]-origin[1])/scale,
                (coords[2,node]-origin[2])/scale,
                (coords[3,node]-origin[3])/scale)
    all(isfinite,normalized) || throw(ArgumentError(
        "mesh_transfinite_patch: generated coordinate cannot be projected"))
    return _project(frame,normalized)
end

function _validate_triangle_orientation(coords,triangles,origin,scale,frame)
    reference=0
    @inbounds for triangle in axes(triangles,2)
        a=_project_output(coords,triangles[1,triangle],origin,scale,frame)
        b=_project_output(coords,triangles[2,triangle],origin,scale,frame)
        c=_project_output(coords,triangles[3,triangle],origin,scale,frame)
        orientation=orient2(a,b,c)
        orientation!=0 || throw(ArgumentError(
            "mesh_transfinite_patch: cell triangle $triangle is folded or degenerate"))
        if reference==0
            reference=orientation
        elseif orientation!=reference
            throw(ArgumentError(
                "mesh_transfinite_patch: cell triangle $triangle reverses patch orientation"))
        end
    end
    return reference
end

function _validate_boundary_postcondition(mesh::Mesh)
    actual,max_incidence=boundary_edges(mesh.tris)
    max_incidence==2 || throw(ErrorException(
        "mesh_transfinite_patch: internal triangle incidence postcondition failed"))
    expected=Vector{NTuple{2,Int32}}(undef,nsegs(mesh))
    @inbounds for segment in 1:nsegs(mesh)
        expected[segment]=_edge_key(mesh.segs[1,segment],mesh.segs[2,segment])
    end
    sort!(actual);sort!(expected)
    actual==expected || throw(ErrorException(
        "mesh_transfinite_patch: triangle boundary does not equal emitted segments"))
    return nothing
end

"""
    mesh_transfinite_patch(side1, side2, side3, side4;
        arrangement=:left, face_tag=0, side_tags=(0,0,0,0),
        max_nodes=10_000_000, max_triangles=20_000_000) -> Mesh

Construct a four-sided planar transfinite triangle patch. Each side is an
already-discretized vector of finite 3-D points, oriented cyclically as
`c1→c2`, `c2→c3`, `c3→c4`, and `c4→c1`. Adjacent endpoints must match exactly
after conversion to `Float64`, and opposite sides must have equal node counts.

The supported Gmsh triangle arrangements are `:left`, `:right`,
`:alternate_left`, and `:alternate_right`. The returned mesh contains the four
boundary segment chains; `face_tag` is copied to every triangle and each entry
of `side_tags` to the corresponding chain. Resource counts and caller limits
are checked before output allocation. The function returns a validated simple
patch or throws a precise blocker; it never falls back to unstructured meshing.

This bounded operation does not discretize curves, apply size or quality fields,
smooth the grid, handle holes or three-sided/quasi-transfinite patches, map a
general CAD parameterization, generate quadrangles, or construct transfinite
volumes. A boundary whose spatial intersection audit exceeds its linear bounded
candidate budget is rejected instead of risking unbounded work.
"""
function mesh_transfinite_patch(side1::AbstractVector,side2::AbstractVector,
                                side3::AbstractVector,side4::AbstractVector;
                                arrangement=:left,face_tag::Integer=0,
                                side_tags=(0,0,0,0),
                                max_nodes::Integer=_DEFAULT_MAX_NODES,
                                max_triangles::Integer=_DEFAULT_MAX_TRIANGLES)::Mesh
    mode=_arrangement(arrangement)
    node_limit=_limit(max_nodes,"max_nodes")
    triangle_limit=_limit(max_triangles,"max_triangles")
    side_tags isa Tuple && length(side_tags)==4 || throw(ArgumentError(
        "mesh_transfinite_patch: side_tags must be a four-integer tuple"))
    physical_side_tags=ntuple(i->_tag(side_tags[i],"side_tags[$i]"),4)
    physical_face_tag=_tag(face_tag,"face_tag")

    lengths=(length(side1),length(side2),length(side3),length(side4))
    @inbounds for i in 1:4
        lengths[i]>=2 || throw(ArgumentError(
            "mesh_transfinite_patch: side $i needs at least two points"))
    end
    lengths[1]==lengths[3] || throw(ArgumentError(
        "mesh_transfinite_patch: opposite sides 1 and 3 have non-matching node counts " *
        "$(lengths[1]) and $(lengths[3])"))
    lengths[2]==lengths[4] || throw(ArgumentError(
        "mesh_transfinite_patch: opposite sides 2 and 4 have non-matching node counts " *
        "$(lengths[2]) and $(lengths[4])"))
    L=lengths[1]-1;H=lengths[2]-1
    nodes=_checked_mul(lengths[1],lengths[2],"node")
    logical_cells=_checked_mul(L,H,"logical-cell")
    triangles=_checked_mul(2,logical_cells,"triangle")
    segments=_checked_mul(2,_checked_add(L,H,"segment"),"segment")
    nodes<=_INT32_MAX || throw(ArgumentError(
        "mesh_transfinite_patch: $nodes nodes exceed the Int32 indexing limit"))
    triangles<=_INT32_MAX || throw(ArgumentError(
        "mesh_transfinite_patch: $triangles triangles exceed the Int32 topology limit"))
    segments<=_INT32_MAX || throw(ArgumentError(
        "mesh_transfinite_patch: $segments segments exceed the Int32 topology limit"))
    nodes<=node_limit || throw(ArgumentError(
        "mesh_transfinite_patch: $nodes nodes exceed max_nodes=$node_limit"))
    triangles<=triangle_limit || throw(ArgumentError(
        "mesh_transfinite_patch: $triangles triangles exceed max_triangles=$triangle_limit"))

    sides=(_convert_side(side1,1),_convert_side(side2,2),
           _convert_side(side3,3),_convert_side(side4,4))
    @inbounds for side in 1:4
        _validate_side_edges(sides[side],side)
        next=mod1(side+1,4)
        sides[side][end]==sides[next][1] || throw(ArgumentError(
            "mesh_transfinite_patch: side $side endpoint does not exactly match " *
            "side $next start point"))
    end
    corners=(sides[1][1],sides[2][1],sides[3][1],sides[4][1])
    length(Set(corners))==4 || throw(ArgumentError(
        "mesh_transfinite_patch: the four corners must be distinct"))

    ring=_boundary_ring(sides)
    origin,scale=_normalization(ring)
    frame=_plane_frame(ring,origin,scale)
    projected_ring=NTuple{2,Float64}[
        _project(frame,_normalize(point,origin,scale)) for point in ring]
    _validate_simple_boundary(projected_ring)

    bottom=_normalized_side(sides[1],origin,scale)
    right=_normalized_side(sides[2],origin,scale)
    top=reverse(_normalized_side(sides[3],origin,scale))
    left=reverse(_normalized_side(sides[4],origin,scale))
    u=_averaged_parameters(bottom,top,"u")
    v=_averaged_parameters(right,left,"v")
    c2=bottom[end];c3=top[end];c4=top[1]

    width=L+1
    coordinates=Matrix{Float64}(undef,3,nodes)
    @inbounds for j in 0:H,i in 0:L
        normalized = if j==0
            bottom[i+1]
        elseif i==L
            right[j+1]
        elseif j==H
            top[i+1]
        elseif i==0
            left[j+1]
        else
            _coons(left[j+1],right[j+1],bottom[i+1],top[i+1],
                   c2,c3,c4,u[i+1],v[j+1])
        end
        all(isfinite,normalized) || throw(ArgumentError(
            "mesh_transfinite_patch: transfinite interpolation generated a non-finite coordinate"))
        point = if j==0
            sides[1][i+1]
        elseif i==L
            sides[2][j+1]
        elseif j==H
            sides[3][L-i+1]
        elseif i==0
            sides[4][H-j+1]
        else
            _physical_point(normalized,origin,scale)
        end
        node=Int(_node(i,j,width))
        coordinates[1,node]=point[1];coordinates[2,node]=point[2]
        coordinates[3,node]=point[3]
    end

    segment_topology=Matrix{Int32}(undef,2,segments)
    segment_tags=Vector{Int32}(undef,segments)
    _fill_segments!(segment_topology,segment_tags,width,L,H,physical_side_tags)
    triangle_topology=Matrix{Int32}(undef,3,triangles)
    _fill_triangles!(triangle_topology,width,L,H,mode)
    _validate_triangle_orientation(coordinates,triangle_topology,origin,scale,frame)
    triangle_tags=fill(physical_face_tag,triangles)

    mesh=Mesh(coordinates;segs=segment_topology,tris=triangle_topology,
              seg_tag=segment_tags,tri_tag=triangle_tags)
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ErrorException(
        "mesh_transfinite_patch: internal output validation failed — " *
        join(diagnostic.messages,"; ")))
    (nnodes(mesh)==nodes&&nsegs(mesh)==segments&&ntris(mesh)==triangles) ||
        throw(ErrorException(
            "mesh_transfinite_patch: internal output count postcondition failed"))
    _validate_boundary_postcondition(mesh)
    return mesh
end

end # module Transfinite
