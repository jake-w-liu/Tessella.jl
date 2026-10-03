# A conservative global certificate for the first NoNewVerts source-quad
# slice. Retained cells use their interval's eight column nodes. Recorded
# problems additionally use a rounded centroid; the subsequent local face-fan
# certificate must prove that actual point lies inside the boundary and hence
# its convex hull. Certifying these hulls therefore certifies the actual P1 cell
# complex, without making an inference from the continuous CAD sweep.
#
# Preconditions are finite, distinct columns, four distinct source indices,
# strictly increasing layer levels, and the subsequent local cell certificate.
# Adjacent hulls must meet only in their
# exactly planar shared cap. Nonadjacent hulls require a strict separating
# plane; unrelated touching boundaries remain unsupported. Supporting triangle
# planes are a sufficient certificate, not a complete intersection test, so an
# undecidable pair receives a precise pending-planner diagnostic.
#
# A balanced logical-interval BVH has 2N-1 fixed-size nodes. Certified supporting
# outer-cap planes prune fine angular ranges whose coordinate AABBs overlap.
# Explicit linear traversal and N*log(N) point-test budgets bound configurations
# whose broad phase would otherwise approach quadratic work. Predicate work per
# final candidate is constant: two sets of at most 56 triangle planes.

const _EXTRUDE_NONEW_GLOBAL_VISITS_PER_INTERVAL = 128

struct _ExtrudeNoNewHullNode
    bounds::NTuple{6,Float64}
    lower::Int32
    upper::Int32
    left::Int32
    right::Int32
    lower_side::Int8
    upper_side::Int8
end

@noinline function _extrude_nonew_global_pending(caller, reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts global P1 cell-hull " *
        "disjointness is not certified ($reason); the cyclic/global " *
        "self-intersection planner is pending"))
end

@inline function _extrude_nonew_hull_bounds(v)
    lo_x=v[1][1];lo_y=v[1][2];lo_z=v[1][3]
    hi_x=lo_x;hi_y=lo_y;hi_z=lo_z
    for i in 2:8
        p=v[i]
        lo_x=min(lo_x,p[1]);lo_y=min(lo_y,p[2]);lo_z=min(lo_z,p[3])
        hi_x=max(hi_x,p[1]);hi_y=max(hi_y,p[2]);hi_z=max(hi_z,p[3])
    end
    return (lo_x,lo_y,lo_z,hi_x,hi_y,hi_z)
end

@inline function _extrude_nonew_hull_union(a,b)
    return (min(a[1],b[1]),min(a[2],b[2]),min(a[3],b[3]),
            max(a[4],b[4]),max(a[5],b[5]),max(a[6],b[6]))
end

@inline function _extrude_nonew_hull_boxes_separate(a,b)
    return a[4]<b[1] || b[4]<a[1] || a[5]<b[2] || b[5]<a[2] ||
           a[6]<b[3] || b[6]<a[3]
end

function _extrude_nonew_hull_range_side(cols,cell,plane,lower,upper,
        point_budget,caller;required=0)
    p=cols[plane,cell[1]];q=cols[plane,cell[2]];r=cols[plane,cell[3]]
    side=required
    for layer in lower:upper,k in 1:4
        point_budget[]-=1
        point_budget[]>=0 || _extrude_nonew_global_pending(caller,
            "the bounded supporting-plane point-test budget is exceeded")
        sign=orient3(p,q,r,cols[layer,cell[k]])
        if required!=0
            sign==required || return Int8(0)
        elseif sign!=0
            if side==0
                side=sign
            elseif side!=sign
                return Int8(0)
            end
        end
    end
    return Int8(side)
end

function _extrude_nonew_hull_tree!(nodes,cols,cell,lower::Int,upper::Int,
        point_budget,caller)
    index=length(nodes)+1
    lower_side=_extrude_nonew_hull_range_side(cols,cell,lower,lower,upper+1,
                                             point_budget,caller)
    upper_side=_extrude_nonew_hull_range_side(cols,cell,upper+1,lower,upper+1,
                                             point_budget,caller)
    if lower==upper
        bounds=_extrude_nonew_hull_bounds(_extrude_nonew_corners(cols,cell,lower))
        push!(nodes,_ExtrudeNoNewHullNode(bounds,Int32(lower),Int32(upper),0,0,
                                        lower_side,upper_side))
        return Int32(index)
    end
    push!(nodes,_ExtrudeNoNewHullNode((0.0,0.0,0.0,0.0,0.0,0.0),
                                    Int32(lower),Int32(upper),0,0,
                                    lower_side,upper_side))
    middle=lower+(upper-lower)÷2
    left=_extrude_nonew_hull_tree!(nodes,cols,cell,lower,middle,point_budget,caller)
    right=_extrude_nonew_hull_tree!(nodes,cols,cell,middle+1,upper,point_budget,caller)
    bounds=_extrude_nonew_hull_union(nodes[left].bounds,nodes[right].bounds)
    nodes[index]=_ExtrudeNoNewHullNode(bounds,Int32(lower),Int32(upper),left,right,
                                    lower_side,upper_side)
    return Int32(index)
end

function _extrude_nonew_hull_caps_separate(a,b,cols,cell,point_budget,caller)
    for (plane,side) in ((Int(a.lower),a.lower_side),(Int(a.upper)+1,a.upper_side))
        side==0 && continue
        _extrude_nonew_hull_range_side(cols,cell,plane,b.lower,b.upper+1,
            point_budget,caller;required=-side)!=0 && return true
    end
    return false
end

# Use the shared cap as an exact supporting plane on both sides. Checking all
# eight nonshared corners also excludes adjacent interval interiors overlapping
# behind their nominal shared cap. Coplanarity uses exact orient3 signs.
function _extrude_nonew_hull_adjacent(cols,cell,layer,caller)
    shared=layer+1
    p=cols[shared,cell[1]];q=cols[shared,cell[2]]
    r=cols[shared,cell[3]];s=cols[shared,cell[4]]
    orient3(p,q,r,s)==0 || _extrude_nonew_global_pending(caller,
        "shared cap of intervals $layer and $(layer+1) is not exactly planar")
    sign=orient3(p,q,r,cols[layer,cell[1]])
    sign!=0 || _extrude_nonew_global_pending(caller,
        "interval $layer has no strict shared-cap halfspace")
    for k in 1:4
        orient3(p,q,r,cols[layer,cell[k]])==sign &&
        orient3(p,q,r,cols[layer+2,cell[k]])==-sign ||
            _extrude_nonew_global_pending(caller,
                "adjacent interval hulls $layer and $(layer+1) lack opposite strict cap halfspaces")
    end
    return nothing
end

# A triangle of one hull is a supporting plane only if all its vertices are
# on one side (or on the plane). Strict separation requires every vertex of
# the other hull to lie on the opposite side. Degenerate triangles provide no
# certificate and are skipped without floating normal/dot calculations.
function _extrude_nonew_hull_support_separates(v,w)
    for a in 1:6,b in a+1:7,c in b+1:8
        side=0
        supporting=true
        for k in 1:8
            sign=orient3(v[a],v[b],v[c],v[k])
            if sign!=0
                if side==0
                    side=sign
                elseif side!=sign
                    supporting=false
                    break
                end
            end
        end
        supporting && side!=0 || continue
        separates=true
        for k in 1:8
            if orient3(v[a],v[b],v[c],w[k])!=-side
                separates=false
                break
            end
        end
        separates && return true
    end
    return false
end

function _extrude_nonew_hull_pair!(nodes,cols,cell,first::Int32,second::Int32,
        remaining::Base.RefValue{Int},point_budget::Base.RefValue{Int},caller)
    remaining[]-=1
    remaining[]>=0 || _extrude_nonew_global_pending(caller,
        "the bounded hull-pair traversal budget is exceeded")
    a=nodes[first];b=nodes[second]
    if first==second
        a.left==0 && return nothing
        _extrude_nonew_hull_pair!(nodes,cols,cell,a.left,a.left,remaining,point_budget,caller)
        _extrude_nonew_hull_pair!(nodes,cols,cell,a.right,a.right,remaining,point_budget,caller)
        _extrude_nonew_hull_pair!(nodes,cols,cell,a.left,a.right,remaining,point_budget,caller)
        return nothing
    end
    _extrude_nonew_hull_boxes_separate(a.bounds,b.bounds) && return nothing
    (_extrude_nonew_hull_caps_separate(a,b,cols,cell,point_budget,caller) ||
     _extrude_nonew_hull_caps_separate(b,a,cols,cell,point_budget,caller)) && return nothing
    if a.left==0 && b.left==0
        abs(a.lower-b.lower)==1 && return nothing
        v=_extrude_nonew_corners(cols,cell,Int(a.lower))
        w=_extrude_nonew_corners(cols,cell,Int(b.lower))
        (_extrude_nonew_hull_support_separates(v,w) ||
         _extrude_nonew_hull_support_separates(w,v)) ||
            _extrude_nonew_global_pending(caller,
                "nonadjacent interval hulls $(a.lower) and $(b.lower) have no verified strict separating plane")
        return nothing
    end
    # Splitting the larger logical range keeps the traversal balanced without
    # sorting or copying columns. Each unordered leaf pair is visited once.
    if b.left==0 || (a.left!=0 && a.upper-a.lower>=b.upper-b.lower)
        _extrude_nonew_hull_pair!(nodes,cols,cell,a.left,second,remaining,point_budget,caller)
        _extrude_nonew_hull_pair!(nodes,cols,cell,a.right,second,remaining,point_budget,caller)
    else
        _extrude_nonew_hull_pair!(nodes,cols,cell,first,b.left,remaining,point_budget,caller)
        _extrude_nonew_hull_pair!(nodes,cols,cell,first,b.right,remaining,point_budget,caller)
    end
    return nothing
end

function _extrude_nonew_global_certify(cols,cell,spec,caller)
    intervals=size(cols,1)-1
    intervals>0 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts global certificate needs at least one interval"))
    intervals<=(_EXTRUDE_NONEW_MAX_NODES-5)÷4 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts global certificate exceeds the node limit"))
    for layer in 1:intervals-1
        _extrude_nonew_hull_adjacent(cols,cell,layer,caller)
    end
    intervals==1 && return nothing
    nodes=_ExtrudeNoNewHullNode[]
    sizehint!(nodes,2intervals-1)
    depth=ndigits(intervals;base=2)
    point_budget=Ref(128intervals*(depth+1))
    root=_extrude_nonew_hull_tree!(nodes,cols,cell,1,intervals,point_budget,caller)
    remaining=Ref(_EXTRUDE_NONEW_GLOBAL_VISITS_PER_INTERVAL*intervals)
    _extrude_nonew_hull_pair!(nodes,cols,cell,root,root,remaining,point_budget,caller)
    return nothing
end
