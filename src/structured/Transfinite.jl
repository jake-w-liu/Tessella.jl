"""
    Transfinite

Validated four-sided planar transfinite patches. Boundary chains are supplied
already discretized; opposite chains must contain the same number of nodes.
The interior interpolation and triangle arrangements reproduce the documented
four-corner Gmsh transfinite-surface construction for an affine planar surface.
"""
module Transfinite

using ..MeshTypes: Mesh, boundary_edges, nnodes, nsegs, ntris, validate
using ..Predicates: orient2, orient3

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

struct _SegmentBox3
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    zmin::Float64
    zmax::Float64
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

struct _BoundaryNode3
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    zmin::Float64
    zmax::Float64
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

function _limit(value,name::AbstractString)
    value isa Integer || throw(ArgumentError(
        "mesh_transfinite_patch: $name must be an integer"))
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
    values=try
        (raw[1],raw[2],raw[3])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "mesh_transfinite_patch: side $side point $index coordinates must " *
            "be Float64-representable: $(sprint(showerror,err))"))
    end
    any(value -> value isa Bool,values) && throw(ArgumentError(
        "mesh_transfinite_patch: side $side point $index coordinates must not be Bool"))
    point=try
        (Float64(values[1]),Float64(values[2]),Float64(values[3]))
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

function _frame_from_normal(normal)
    reference=abs(normal[1])<=abs(normal[2]) ?
        (abs(normal[1])<=abs(normal[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
        (abs(normal[2])<=abs(normal[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    axis_u=_cross3(reference,normal);axis_length=_norm3(axis_u)
    (isfinite(axis_length)&&axis_length>0) || throw(ArgumentError(
        "mesh_transfinite_patch: could not construct an in-plane frame"))
    axis_u=(axis_u[1]/axis_length,axis_u[2]/axis_length,axis_u[3]/axis_length)
    axis_v=_cross3(normal,axis_u)
    return _PlaneFrame(axis_u,axis_v,normal)
end

function _exact_plane_frame(ring,origin)
    exact_origin=ntuple(d->Rational{BigInt}(origin[d]),3)
    points=Vector{NTuple{3,Rational{BigInt}}}(undef,length(ring))
    @inbounds for i in eachindex(ring)
        points[i]=ntuple(d->Rational{BigInt}(ring[i][d])-exact_origin[d],3)
    end
    first=1;first_norm=_dot3(points[1],points[1])
    @inbounds for i in 2:length(points)
        point_norm=_dot3(points[i],points[i])
        if point_norm>first_norm
            first=i;first_norm=point_norm
        end
    end
    axis=points[first]
    zero_exact=Rational{BigInt}(0)
    best=(zero_exact,zero_exact,zero_exact);best_norm=zero_exact
    @inbounds for point in points
        candidate=_cross3(axis,point)
        candidate_norm=_dot3(candidate,candidate)
        if candidate_norm>best_norm
            best=candidate;best_norm=candidate_norm
        end
    end
    best_norm>0 || return nothing
    @inbounds for point in points
        _dot3(best,point)==0 || return nothing
    end
    normal=setprecision(BigFloat,256) do
        length_big=sqrt(best_norm)
        (Float64(BigFloat(best[1])/length_big),
         Float64(BigFloat(best[2])/length_big),
         Float64(BigFloat(best[3])/length_big))
    end
    all(isfinite,normal)&&_norm3(normal)>0 || return nothing
    return _frame_from_normal(normal)
end



# `allow_warped` admits a genuinely non-coplanar boundary: the patch then
# follows the 3-D Coons interpolation of its side chains instead of living in
# a single plane. The Newell normal of the closed ring still defines the
# orientation frame; it remains valid for any non-self-intersecting closed
# polygon. Exactly-coplanar rings take the planar frame unchanged so planar
# output is bit-identical either way.
#
# `_ring_coplanar_exact` distinguishes a truly warped ring from a coplanar
# ring with Float64 accumulation noise: three non-collinear anchors span a
# candidate plane and every remaining vertex must satisfy orient3 == 0 — an
# exact-sign predicate that allocates nothing on the ordinary float path.
function _ring_coplanar_exact(ring)
    n=length(ring)
    n<3 && return false
    a=ring[1];b_index=0;best_norm=-1.0
    @inbounds for i in 2:n
        d=_sub3(ring[i],a);norm2=_dot3(d,d)
        norm2>best_norm && (best_norm=norm2;b_index=i)
    end
    b_index==0 && return false
    b=ring[b_index];ab=_sub3(b,a)
    c_index=0;best_area=0.0
    @inbounds for i in 2:n
        i==b_index && continue
        cross=_cross3(ab,_sub3(ring[i],a));area2=_dot3(cross,cross)
        area2>best_area && (best_area=area2;c_index=i)
    end
    if c_index==0
        # Every Float64 cross product canceled to zero (exactly representable
        # near-collinear coordinates): resolve the third anchor exactly so a
        # valid affine patch still receives its coplanarity certificate.
        ea=ntuple(d->Rational{BigInt}(a[d]),3)
        eb=ntuple(d->Rational{BigInt}(b[d]),3)
        eab=_sub3(eb,ea)
        @inbounds for i in 2:n
            i==b_index && continue
            ep=ntuple(d->Rational{BigInt}(ring[i][d]),3)
            exact_cross=_cross3(eab,_sub3(ep,ea))
            _dot3(exact_cross,exact_cross)>0 || continue
            c_index=i;break
        end
        c_index==0 && return false
    end
    c=ring[c_index]
    @inbounds for i in eachindex(ring)
        orient3(a,b,c,ring[i])!=0 && return false
    end
    return true
end

function _patch_frame(ring,origin,scale,allow_warped::Bool)
    nx=0.0;ny=0.0;nz=0.0
    @inbounds for i in eachindex(ring)
        p=_normalize(ring[i],origin,scale)
        q=_normalize(ring[mod1(i+1,length(ring))],origin,scale)
        nx+=(p[2]-q[2])*(p[3]+q[3])
        ny+=(p[3]-q[3])*(p[1]+q[1])
        nz+=(p[1]-q[1])*(p[2]+q[2])
    end
    normal_length=hypot(nx,ny,nz)
    if !(isfinite(normal_length)&&normal_length>0)
        exact=_exact_plane_frame(ring,origin)
        exact===nothing && throw(ArgumentError(
            "mesh_transfinite_patch: boundary has no representable plane normal"))
        return exact,true
    end
    normal=(nx/normal_length,ny/normal_length,nz/normal_length)
    tolerance=256eps(Float64)
    warped=false
    @inbounds for point in ring
        normalized=_normalize(point,origin,scale)
        distance=abs(_dot3(normalized,normal))
        if !(isfinite(distance)&&distance<=tolerance)
            warped=true;break
        end
    end
    if warped
        _ring_coplanar_exact(ring) &&
            return _frame_from_normal(normal),true
        allow_warped || throw(ArgumentError(
            "mesh_transfinite_patch: boundary is not coplanar and the patch " *
            "is not allowed to be warped"))
    end
    return _frame_from_normal(normal),!warped
end

# In exact arithmetic, two non-coplanar open segments cannot intersect, so a
# 3-D intersection audit only needs the coplanar case plus the shared-endpoint
# handling of adjacent ring edges. Coplanar candidates are tested in their
# best-conditioned coordinate-plane projection; fully collinear candidates
# fall back to a scalar interval-overlap test along the dominant axis.
function _point_on_segment3(a,b,p)
    cross=_cross3(_sub3(b,a),_sub3(p,a))
    if !(cross[1]==0&&cross[2]==0&&cross[3]==0)
        ea=ntuple(d->Rational{BigInt}(a[d]),3)
        eb=ntuple(d->Rational{BigInt}(b[d]),3)
        ep=ntuple(d->Rational{BigInt}(p[d]),3)
        cross_exact=_cross3(_sub3(eb,ea),_sub3(ep,ea))
        (cross_exact[1]==0&&cross_exact[2]==0&&cross_exact[3]==0) ||
            return false
    end
    return min(a[1],b[1])<=p[1]<=max(a[1],b[1]) &&
           min(a[2],b[2])<=p[2]<=max(a[2],b[2]) &&
           min(a[3],b[3])<=p[3]<=max(a[3],b[3])
end

@inline _plane_axes(pair::Int)=pair==1 ? (2,3) : pair==2 ? (1,3) : (1,2)
@inline _proj_axes(p,axes)=(p[axes[1]],p[axes[2]])

function _segments_intersect_3d(a,b,c,d)
    orient3(a,b,c,d)!=0 && return false
    axes=0
    @inbounds for tri in ((a,b,c),(a,b,d),(c,d,a),(c,d,b))
        for pair in ((1,2),(1,3),(2,3))
            orient2(_proj_axes(tri[1],pair),_proj_axes(tri[2],pair),
                    _proj_axes(tri[3],pair))!=0 || continue
            axes=pair;break
        end
        axes!=0 && break
    end
    if axes==0
        # All four endpoints are collinear: compare interval ranges along the
        # dominant coordinate of the common direction.
        ab=_sub3(b,a)
        axis=abs(ab[1])>=abs(ab[2])&&abs(ab[1])>=abs(ab[3]) ? 1 :
              abs(ab[2])>=abs(ab[3]) ? 2 : 3
        lo1,hi1=minmax(a[axis],b[axis]);lo2,hi2=minmax(c[axis],d[axis])
        return max(lo1,lo2)<=min(hi1,hi2)
    end
    return _segments_intersect(_proj_axes(a,axes),_proj_axes(b,axes),
                               _proj_axes(c,axes),_proj_axes(d,axes))
end

function _adjacent_overlap_3d(a,b,c,d)
    shared = a==c || a==d ? a : b==c || b==d ? b : nothing
    shared===nothing && return true
    for point in (a,b)
        point==shared || !_point_on_segment3(c,d,point) || return true
    end
    for point in (c,d)
        point==shared || !_point_on_segment3(a,b,point) || return true
    end
    return false
end

@inline _box_overlap_3d(a,b)=a.xmin<=b.xmax&&b.xmin<=a.xmax&&
                                   a.ymin<=b.ymax&&b.ymin<=a.ymax&&
                                   a.zmin<=b.zmax&&b.zmin<=a.zmax

function _build_boundary_tree_3d!(nodes,order,boxes,lo::Int,hi::Int)
    xmin=Inf;xmax=-Inf;ymin=Inf;ymax=-Inf;zmin=Inf;zmax=-Inf
    @inbounds for position in lo:hi
        box=boxes[order[position]]
        xmin=min(xmin,box.xmin);xmax=max(xmax,box.xmax)
        ymin=min(ymin,box.ymin);ymax=max(ymax,box.ymax)
        zmin=min(zmin,box.zmin);zmax=max(zmax,box.zmax)
    end
    count=hi-lo+1;node_index=length(nodes)+1
    push!(nodes,_BoundaryNode3(xmin,xmax,ymin,ymax,zmin,zmax,0,0,0,count))
    if lo==hi
        @inbounds segment=boxes[order[lo]].segment
        nodes[node_index]=_BoundaryNode3(xmin,xmax,ymin,ymax,zmin,zmax,
                                         segment,0,0,1)
        return node_index
    end
    xspan=xmax-xmin;yspan=ymax-ymin;zspan=zmax-zmin
    axis=xspan>=yspan&&xspan>=zspan ? 1 : yspan>=zspan ? 2 : 3
    sort!(@view(order[lo:hi]);
          by=index->begin box=boxes[index]
              axis==1 ? box.xmin/2+box.xmax/2 :
              axis==2 ? box.ymin/2+box.ymax/2 : box.zmin/2+box.zmax/2 end,
          alg=QuickSort)
    middle=lo+(hi-lo)÷2
    left=_build_boundary_tree_3d!(nodes,order,boxes,lo,middle)
    right=_build_boundary_tree_3d!(nodes,order,boxes,middle+1,hi)
    nodes[node_index]=_BoundaryNode3(xmin,xmax,ymin,ymax,zmin,zmax,
                                     0,left,right,count)
    return node_index
end

function _audit_boundary_pair_3d(points,i::Int,j::Int,count::Int)
    a=points[i];b=points[mod1(i+1,count)]
    c=points[j];d=points[mod1(j+1,count)]
    _segments_intersect_3d(a,b,c,d) || return nothing
    adjacent=_boundary_adjacent(min(i,j),max(i,j),count)
    if !adjacent || _adjacent_overlap_3d(a,b,c,d)
        throw(ArgumentError(
            "mesh_transfinite_patch: boundary segments $i and $j intersect"))
    end
    return nothing
end

function _validate_simple_boundary_3d(points)
    count=length(points)
    boxes=Vector{_SegmentBox3}(undef,count)
    @inbounds for i in 1:count
        a=points[i];b=points[mod1(i+1,count)]
        (isfinite(a[1])&&isfinite(a[2])&&isfinite(a[3])) ||
            throw(ArgumentError(
                "mesh_transfinite_patch: boundary node $i is not finite"))
        boxes[i]=_SegmentBox3(min(a[1],b[1]),max(a[1],b[1]),
                              min(a[2],b[2]),max(a[2],b[2]),
                              min(a[3],b[3]),max(a[3],b[3]),i)
    end
    order=collect(1:count);nodes=_BoundaryNode3[]
    node_capacity=_checked_add(_checked_mul(2,count,"boundary audit node"),-1,
                               "boundary audit node")
    sizehint!(nodes,node_capacity)
    root=_build_boundary_tree_3d!(nodes,order,boxes,1,count)
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
        _box_overlap_3d(node1,node2) || continue
        if first==second
            node1.segment!=0 && continue
            push!(stack,(node1.left,node1.left),(node1.left,node1.right),
                        (node1.right,node1.right))
        elseif node1.segment!=0&&node2.segment!=0
            candidates+=1
            candidates<=limit || throw(ArgumentError(
                "mesh_transfinite_patch: boundary intersection audit exceeded " *
                "its bounded candidate limit $limit"))
            _audit_boundary_pair_3d(points,node1.segment,node2.segment,count)
        elseif node2.segment!=0 || (node1.segment==0&&node1.count>=node2.count)
            push!(stack,(node1.left,second),(node1.right,second))
        else
            push!(stack,(first,node2.left),(first,node2.right))
        end
    end
    return nothing
end

# Warped patches cannot be folded back to a reference plane, so orientation is
# audited directly in 3-D: every output triangle must carry a nonzero area
# (with an exact-arithmetic fallback when Float64 squares underflow), and the
# area-weighted normal field must agree with the ring's Newell normal — which
# rejects an inverted patch. Local folds between adjacent triangles are
# rejected by requiring every interior shared edge to keep both incident
# triangle normals on the same side of the dihedral (dot >= 0).
function _triangle_normal_3d(coords,i1,i2,i3)
    ax,ay,az=coords[1,i1],coords[2,i1],coords[3,i1]
    bx,by,bz=coords[1,i2],coords[2,i2],coords[3,i2]
    cx,cy,cz=coords[1,i3],coords[2,i3],coords[3,i3]
    return _cross3((bx-ax,by-ay,bz-az),(cx-ax,cy-ay,cz-az))
end

function _validate_warped_patch(coords,triangles,ring,origin,scale,L::Int,H::Int)
    nx=0.0;ny=0.0;nz=0.0
    @inbounds for i in eachindex(ring)
        p=_normalize(ring[i],origin,scale)
        q=_normalize(ring[mod1(i+1,length(ring))],origin,scale)
        nx+=(p[2]-q[2])*(p[3]+q[3])
        ny+=(p[3]-q[3])*(p[1]+q[1])
        nz+=(p[1]-q[1])*(p[2]+q[2])
    end
    normals=Vector{NTuple{3,Float64}}(undef,size(triangles,2))
    gx=0.0;gy=0.0;gz=0.0
    @inbounds for triangle in axes(triangles,2)
        n=_triangle_normal_3d(coords,Int(triangles[1,triangle]),
                              Int(triangles[2,triangle]),
                              Int(triangles[3,triangle]))
        n2=_dot3(n,n)
        if !(isfinite(n2)&&n2>0)
            ea=ntuple(d->Rational{BigInt}(coords[d,Int(triangles[1,triangle])]),3)
            eb=ntuple(d->Rational{BigInt}(coords[d,Int(triangles[2,triangle])]),3)
            ec=ntuple(d->Rational{BigInt}(coords[d,Int(triangles[3,triangle])]),3)
            en=_cross3(_sub3(eb,ea),_sub3(ec,ea))
            _dot3(en,en)>0 || throw(ArgumentError(
                "mesh_transfinite_patch: cell triangle $triangle is degenerate"))
            n=ntuple(d->Float64(en[d]),3)
            nl=sqrt(_dot3(n,n));n=(n[1]/nl,n[2]/nl,n[3]/nl)
        end
        normals[triangle]=n
        gx+=n[1];gy+=n[2];gz+=n[3]
    end
    gx*nx+gy*ny+gz*nz>0 || throw(ArgumentError(
        "mesh_transfinite_patch: warped patch reverses boundary orientation"))
    # Edge-sharing normal consistency enumerated analytically over the
    # regular cell grid — no auxiliary incidence storage beyond `normals`.
    # Cell (i,j) owns triangles 2(i·H+j)+1 and 2(i·H+j)+2, which share the
    # cell diagonal; cells (i,j) and (i+1,j) share the grid edge between
    # nodes (i+1,j) and (i+1,j+1); cells (i,j) and (i,j+1) share the edge
    # between nodes (i,j+1) and (i+1,j+1).
    @inbounds for i in 0:L-1,j in 0:H-1
        base=2(i*H+j)
        _dot3(normals[base+1],normals[base+2])>=0 || throw(ArgumentError(
            "mesh_transfinite_patch: warped patch folds across a shared edge"))
        if i<L-1
            a=_node(i+1,j,L+1);b=_node(i+1,j+1,L+1)
            _dot3(normals[_edge_owner(triangles,base+1,base+2,a,b)],
                  normals[_edge_owner(triangles,base+2H+1,base+2H+2,a,b)])>=0 ||
                throw(ArgumentError(
                    "mesh_transfinite_patch: warped patch folds across " *
                    "a shared edge"))
        end
        if j<H-1
            a=_node(i,j+1,L+1);b=_node(i+1,j+1,L+1)
            _dot3(normals[_edge_owner(triangles,base+1,base+2,a,b)],
                  normals[_edge_owner(triangles,base+3,base+4,a,b)])>=0 ||
                throw(ArgumentError(
                    "mesh_transfinite_patch: warped patch folds across " *
                    "a shared edge"))
        end
    end
    return 1
end

# The triangle (of the two owned by a cell) containing both endpoints of a
# shared grid edge.
@inline function _edge_owner(triangles,t1::Int,t2::Int,a,b)
    _edge_contains(triangles,t1,a,b) && return t1
    _edge_contains(triangles,t2,a,b) && return t2
    throw(ErrorException(
        "mesh_transfinite_patch: internal warped-audit adjacency failed"))
end

@inline _edge_contains(triangles,t,a,b)=
    _tri_contains(triangles,t,a)&&_tri_contains(triangles,t,b)
@inline _tri_contains(triangles,t,v)=
    triangles[1,t]==v||triangles[2,t]==v||triangles[3,t]==v

@inline _project(frame::_PlaneFrame,point)=
    (_dot3(point,frame.u),_dot3(point,frame.v))

# Orthogonal projection of `point` onto the plane through `anchor` with unit
# `normal` — the `planeSurface` reparametrization of a boundary vertex whose
# true position sits off the surface plane.
@inline function _project_onto_plane(point,anchor,normal)
    d=_dot3(_sub3(point,anchor),normal)
    return _sub3(point,(d*normal[1],d*normal[2],d*normal[3]))
end

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
        allow_warped=false,
        max_nodes=10_000_000, max_triangles=20_000_000) -> Mesh

Construct a four-sided transfinite triangle patch. Each side is an
already-discretized vector of finite 3-D points, oriented cyclically as
`c1→c2`, `c2→c3`, `c3→c4`, and `c4→c1`. Adjacent endpoints must match exactly
after conversion to `Float64`, and opposite sides must have equal node counts.
By default the boundary must be coplanar; `allow_warped=true` instead admits a
genuinely non-coplanar ring, whose interior follows the three-dimensional
Coons interpolation of the four side chains (the ruled-surface analogue).

The supported Gmsh triangle arrangements are `:left`, `:right`,
`:alternate_left`, and `:alternate_right`. The returned mesh contains the four
boundary segment chains; `face_tag` is copied to every triangle and each entry
of `side_tags` to the corresponding chain. Resource counts and caller limits
are checked before output allocation. The function returns a validated simple
patch or throws a precise blocker; it never falls back to unstructured meshing.
Warped patches are audited in 3-D: the boundary must be a non-intersecting
closed ring, every output triangle must have nonzero area, and the
area-weighted normal field must agree with the ring's orientation while no
pair of adjacent triangles may fold across their shared edge.

This bounded operation does not discretize curves, apply size or quality fields,
smooth the grid, handle holes or three-sided/quasi-transfinite patches, map a
general CAD parameterization, generate quadrangles, or construct transfinite
volumes. A boundary whose spatial intersection audit exceeds its linear bounded
candidate budget is rejected instead of risking unbounded work.
"""
function mesh_transfinite_patch(side1::AbstractVector,side2::AbstractVector,
                                side3::AbstractVector,side4::AbstractVector;
                                arrangement=:left,face_tag=0,
                                side_tags=(0,0,0,0),
                                allow_warped::Bool=false,
                                project_plane=nothing,
                                max_nodes=_DEFAULT_MAX_NODES,
                                max_triangles=_DEFAULT_MAX_TRIANGLES)::Mesh
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

    # `project_plane` is Gmsh's `planeSurface` transfinite semantics for a
    # non-coplanar boundary: the side chains are projected onto the declared
    # plane for the (u,v) bookkeeping and interior interpolation (the
    # interior stays exactly planar), while emitted boundary nodes keep the
    # true positions. Mutually exclusive with `allow_warped`.
    allow_warped && project_plane!==nothing && throw(ArgumentError(
        "mesh_transfinite_patch: project_plane and allow_warped are " *
        "mutually exclusive"))
    interp_sides=sides
    if project_plane!==nothing
        (project_plane isa Tuple && length(project_plane)==2) ||
            throw(ArgumentError("mesh_transfinite_patch: project_plane must " *
                                "be an (anchor, normal) tuple"))
        anchor,plane_normal=project_plane
        normal_norm=_norm3(plane_normal)
        (isfinite(normal_norm)&&normal_norm>0) || throw(ArgumentError(
            "mesh_transfinite_patch: project_plane normal is degenerate"))
        unit_normal=(plane_normal[1]/normal_norm,plane_normal[2]/normal_norm,
                     plane_normal[3]/normal_norm)
        interp_sides=ntuple(4) do i
            [_project_onto_plane(p,anchor,unit_normal) for p in sides[i]]
        end
    end
    ring=_boundary_ring(interp_sides)
    origin,scale=_normalization(ring)
    frame,planar=_patch_frame(ring,origin,scale,allow_warped)
    if planar
        projected_ring=NTuple{2,Float64}[
            _project(frame,_normalize(point,origin,scale)) for point in ring]
        _validate_simple_boundary(projected_ring)
    else
        normalized_ring=NTuple{3,Float64}[
            _normalize(point,origin,scale) for point in ring]
        _validate_simple_boundary_3d(normalized_ring)
    end
    emit_ring=ring
    if project_plane!==nothing
        emit_ring=_boundary_ring(sides)
        _validate_simple_boundary_3d(NTuple{3,Float64}[
            _normalize(point,origin,scale) for point in emit_ring])
    end

    bottom=_normalized_side(interp_sides[1],origin,scale)
    right=_normalized_side(interp_sides[2],origin,scale)
    top=reverse(_normalized_side(interp_sides[3],origin,scale))
    left=reverse(_normalized_side(interp_sides[4],origin,scale))
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
    if planar && project_plane===nothing
        _validate_triangle_orientation(coordinates,triangle_topology,
                                       origin,scale,frame)
    else
        _validate_warped_patch(coordinates,triangle_topology,emit_ring,
                               origin,scale,L,H)
    end
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
