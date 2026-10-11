module QualityFactors

# Derived from a fixed 3-by-3 Gaussian factorization and barycentric gradients.
# The finite arithmetic policy is specified in the review contract: initial
# row scaling, first equal pivot, ordered products, and exact zero rejection.
# This is a separate tuple factorization, not a transcription of a general LU
# implementation. Reference implementations stay outside production inputs.
const Row3 = NTuple{3,Float64}
struct Factors3
    rows::NTuple{3,Row3}
    origins::NTuple{3,Int}
end

@inline rowmaximum(r::Row3) = max(abs(r[1]),abs(r[2]),abs(r[3]))
@inline function eliminate_first(r::Row3, pivot::Row3)
    multiplier = r[1] / pivot[1]
    return multiplier == 0.0 ? (multiplier,r[2],r[3]) :
        (multiplier,r[2]-multiplier*pivot[2],r[3]-multiplier*pivot[3])
end

@inline function factor3(rows::NTuple{3,Row3})
    maxima = (rowmaximum(rows[1]),rowmaximum(rows[2]),rowmaximum(rows[3]))
    any(iszero,maxima) && return nothing
    scales = (1.0/maxima[1],1.0/maxima[2],1.0/maxima[3])
    scores = (abs(rows[1][1])*scales[1],abs(rows[2][1])*scales[2],
              abs(rows[3][1])*scales[3])
    best = scores[1] > 0.0 ? scores[1] : 0.0
    pivot = 1
    if scores[2] > best
        pivot = 2; best = scores[2]
    end
    scores[3] > best && (pivot = 3; best = scores[3])
    best > 0.0 || return nothing
    # A swap of the selected row with row one retains the remaining row order.
    origins = pivot == 1 ? (1,2,3) : pivot == 2 ? (2,1,3) : (3,2,1)
    first = rows[origins[1]]
    second = eliminate_first(rows[origins[2]],first)
    third = eliminate_first(rows[origins[3]],first)
    raw_score2 = abs(second[2])*scales[origins[2]]
    score2 = raw_score2 > 0.0 ? raw_score2 : 0.0
    score3 = abs(third[2])*scales[origins[3]]
    if score3 > score2
        second,third = third,second
        origins = (origins[1],origins[3],origins[2])
        score2 = score3
    end
    score2 > 0.0 || return nothing
    multiplier = third[2] / second[2]
    third = (third[1],multiplier,multiplier == 0.0 ? third[3] :
                               third[3]-multiplier*second[3])
    third[3] == 0.0 && return nothing
    return Factors3((first,second,third),origins)
end

# Forward and reverse triangular equations, expanded for dimension three.
@inline function solve3(factors::Factors3, rhs::Row3)
    r1,r2,r3 = factors.rows
    p1,p2,p3 = factors.origins
    y1 = rhs[p1] - 0.0
    y2 = rhs[p2] - (0.0+r2[1]*y1)
    y3 = rhs[p3] - ((0.0+r3[1]*y1)+r3[2]*y2)
    x3 = (y3-0.0)/r3[3]
    x2 = (y2-(0.0+r2[3]*x3))/r2[2]
    x1 = (y1-((0.0+r1[2]*x2)+r1[3]*x3))/r1[1]
    return (x1,x2,x3)
end

@inline subtract3(a,b) = (a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline dot3(a,b) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline negative_cross(a,b) = (-(a[2]*b[3]-a[3]*b[2]),
    -(a[3]*b[1]-a[1]*b[3]),-(a[1]*b[2]-a[2]*b[1]))

function face_normal(a,b,c)
    edges = (subtract3(b,a),subtract3(a,c),subtract3(c,b))
    lengths = (dot3(edges[1],edges[1]),dot3(edges[2],edges[2]),
               dot3(edges[3],edges[3]))
    # Omit the first longest edge. The other two, in cyclic order, define the
    # same oriented face while avoiding an unnecessarily long-edge cross.
    omitted = lengths[2] > lengths[1] ? 2 : 1
    lengths[3] > lengths[omitted] && (omitted = 3)
    first,second = omitted == 1 ? (2,3) : omitted == 2 ? (3,1) : (1,2)
    return negative_cross(edges[first],edges[second])
end

@inline function normalize3(normal)
    length = sqrt(dot3(normal,normal))
    length == 0.0 && return (false,normal)
    return (true,(normal[1]/length,normal[2]/length,normal[3]/length))
end

@inline function normalize_four(n1,n2,n3,n4)
    v1,u1=normalize3(n1);v2,u2=normalize3(n2)
    v3,u3=normalize3(n3);v4,u4=normalize3(n4)
    return (v1 && v2 && v3 && v4,(u1,u2,u3,u4))
end

function tetrahedron_cosines(points)
    a,b,c,d = points
    factors = factor3((subtract3(a,d),subtract3(b,d),subtract3(c,d)))
    zero = (0.0,0.0,0.0)
    normals = (zero,zero,zero,zero)
    valid = false
    if factors !== nothing
        r1,r2,r3 = factors.rows
        volume = abs(r1[1]*r2[2]*r3[3])/6.0
        if volume > 0.0
            n1 = solve3(factors,(1.0,0.0,0.0))
            n2 = solve3(factors,(0.0,1.0,0.0))
            n3 = solve3(factors,(0.0,0.0,1.0))
            n4 = (-n1[1]-n2[1]-n3[1],-n1[2]-n2[2]-n3[2],
                  -n1[3]-n2[3]-n3[3])
            valid,normals = normalize_four(n1,n2,n3,n4)
        end
    end
    if !valid
        valid,normals = normalize_four(face_normal(c,b,d),face_normal(a,c,d),
                                      face_normal(b,a,d),face_normal(a,b,c))
    end
    valid || return (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
    n1,n2,n3,n4 = normals
    return (clamp(-dot3(n1,n2),-1.0,1.0),clamp(-dot3(n2,n3),-1.0,1.0),
            clamp(-dot3(n3,n4),-1.0,1.0),clamp(-dot3(n1,n4),-1.0,1.0),
            clamp(-dot3(n3,n1),-1.0,1.0),clamp(-dot3(n2,n4),-1.0,1.0))
end
end
