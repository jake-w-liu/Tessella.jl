"""
Independent exact-arithmetic oracles for CRC verification (DEVELOPMENT.md §CRC).

These deliberately do **not** reuse any production code path. Predicates are
checked by building the *full homogeneous* determinant matrix (3×3 / 4×4 / 5×5)
in exact `Rational{BigInt}` and evaluating it with a **generic recursive Laplace
expansion** — a structurally different computation from the hand-expanded minor
forms in `src/Predicates.jl`. Agreement across the two independent routes is the
oracle; a transcription/sign/index bug in either shows up as a mismatch.
"""
module Oracles

export exact_det, oracle_orient2, oracle_orient3, oracle_incircle, oracle_incircle3, oracle_insphere

const QBig = Rational{BigInt}

# Generic exact determinant by recursive Laplace (cofactor) expansion along the
# first row. Independent of any production expansion order. O(n!) — fine for the
# n ≤ 5 homogeneous matrices used here.
function exact_det(M::Matrix{QBig})::QBig
    n = size(M, 1)
    @assert size(M, 2) == n "matrix must be square"
    n == 1 && return M[1, 1]
    n == 2 && return M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]
    acc = zero(QBig)
    sgn = one(QBig)
    cols = collect(2:n)
    for j in 1:n
        keep = [c for c in 1:n if c != j]
        minor = M[2:n, keep]
        acc += sgn * M[1, j] * exact_det(minor)
        sgn = -sgn
    end
    return acc
end

_q(x)::QBig = convert(QBig, Float64(x))
_isign(x)::Int = x > 0 ? 1 : (x < 0 ? -1 : 0)

"""Sign of the 3×3 homogeneous orientation determinant of (a,b,c)."""
function oracle_orient2(a, b, c)::Int
    M = QBig[_q(a[1]) _q(a[2]) 1;
             _q(b[1]) _q(b[2]) 1;
             _q(c[1]) _q(c[2]) 1]
    return _isign(exact_det(M))
end

"""Sign of the 4×4 homogeneous orientation determinant of (a,b,c,d)."""
function oracle_orient3(a, b, c, d)::Int
    M = QBig[_q(a[1]) _q(a[2]) _q(a[3]) 1;
             _q(b[1]) _q(b[2]) _q(b[3]) 1;
             _q(c[1]) _q(c[2]) _q(c[3]) 1;
             _q(d[1]) _q(d[2]) _q(d[3]) 1]
    return _isign(exact_det(M))
end

"""Sign of the 4×4 homogeneous in-circle determinant of (a,b,c,d)."""
function oracle_incircle(a, b, c, d)::Int
    la = _q(a[1])^2 + _q(a[2])^2
    lb = _q(b[1])^2 + _q(b[2])^2
    lc = _q(c[1])^2 + _q(c[2])^2
    ld = _q(d[1])^2 + _q(d[2])^2
    M = QBig[_q(a[1]) _q(a[2]) la 1;
             _q(b[1]) _q(b[2]) lb 1;
             _q(c[1]) _q(c[2]) lc 1;
             _q(d[1]) _q(d[2]) ld 1]
    return _isign(exact_det(M))
end

"""Exact in-circle determinant for coplanar 3-D points using a non-degenerate projection and the full 3-D squared norm."""
function oracle_incircle3(a,b,c,d)::Int
    axes=oracle_orient2((a[1],a[2]),(b[1],b[2]),(c[1],c[2]))!=0 ? (1,2) :
         oracle_orient2((a[2],a[3]),(b[2],b[3]),(c[2],c[3]))!=0 ? (2,3) : (3,1)
    lift(p)=_q(p[1])^2+_q(p[2])^2+_q(p[3])^2
    M=QBig[_q(a[axes[1]]) _q(a[axes[2]]) lift(a) 1;
           _q(b[axes[1]]) _q(b[axes[2]]) lift(b) 1;
           _q(c[axes[1]]) _q(c[axes[2]]) lift(c) 1;
           _q(d[axes[1]]) _q(d[axes[2]]) lift(d) 1]
    return _isign(exact_det(M))
end

"""Sign of the 5×5 homogeneous in-sphere determinant of (a,b,c,d,e)."""
function oracle_insphere(a, b, c, d, e)::Int
    lift(p) = _q(p[1])^2 + _q(p[2])^2 + _q(p[3])^2
    M = QBig[_q(a[1]) _q(a[2]) _q(a[3]) lift(a) 1;
             _q(b[1]) _q(b[2]) _q(b[3]) lift(b) 1;
             _q(c[1]) _q(c[2]) _q(c[3]) lift(c) 1;
             _q(d[1]) _q(d[2]) _q(d[3]) lift(d) 1;
             _q(e[1]) _q(e[2]) _q(e[3]) lift(e) 1]
    return _isign(exact_det(M))
end

end # module Oracles
