"""
    MeshTypes

Compact, cache-friendly mesh container and its topology / quality / checksum
machinery (PLAN.md §3 "MeshTypes"; DEVELOPMENT.md "Mesh-CRC checksum").

Storage is column-major and index-compact:

* `coords :: Matrix{Float64}` is `3 × nnodes`; node `i` occupies the contiguous
  column `i`, so `node(m,i)` returns an `NTuple{3,Float64}` with **zero
  allocation**. Planar meshes carry `z = 0`.
* Cell connectivity is `Matrix{Int32}` (`2×`, `3×`, `4×` for segments, triangles,
  tets). `Int32` halves memory versus `Int` and still addresses 2.1 G nodes.

This is the *finalized* container used for I/O, validation and regression
checksums. Meshing kernels keep their own growable working structures and emit a
`Mesh` at the end.

Everything here is verified against independent oracles (analytic areas/volumes,
Euler characteristic of known complexes, hand-computed regular-tetrahedron
quality) in `test/meshtypes_test.jl`.
"""
module MeshTypes

import SHA

export Mesh, nnodes, nsegs, ntris, ntets, node, node2
export triangle_area, tet_volume, tet_signed_volume
export tet_dihedral_extrema, tet_circumradius, tet_radius_edge
export boundary_faces, boundary_edges, unique_edges, unique_faces
export euler_characteristic, boundary_euler, is_closed_manifold
export bounding_box, validate, MeshDiagnostic
export mesh_crc, crc_string

# ════════════════════════════════════════════════════════════════════════════════
# Container
# ════════════════════════════════════════════════════════════════════════════════

"""
    Mesh(coords; segs, tris, tets, seg_tag, tri_tag, tet_tag)

Finalized simplex mesh. `coords` is `3 × nnodes` (`Matrix{Float64}`); `segs`,
`tris`, `tets` are `2×`, `3×`, `4×` `Int32` connectivity (1-based). Tag vectors,
when supplied, are per-cell region/physical ids; otherwise zero-filled.
"""
struct Mesh
    coords::Matrix{Float64}   # 3 × nnodes
    segs::Matrix{Int32}       # 2 × nsegs
    tris::Matrix{Int32}       # 3 × ntris
    tets::Matrix{Int32}       # 4 × ntets
    seg_tag::Vector{Int32}
    tri_tag::Vector{Int32}
    tet_tag::Vector{Int32}

    function Mesh(coords::AbstractMatrix{<:Real};
                  segs::AbstractMatrix{<:Integer}=Matrix{Int32}(undef, 2, 0),
                  tris::AbstractMatrix{<:Integer}=Matrix{Int32}(undef, 3, 0),
                  tets::AbstractMatrix{<:Integer}=Matrix{Int32}(undef, 4, 0),
                  seg_tag::AbstractVector{<:Integer}=zeros(Int32, size(segs, 2)),
                  tri_tag::AbstractVector{<:Integer}=zeros(Int32, size(tris, 2)),
                  tet_tag::AbstractVector{<:Integer}=zeros(Int32, size(tets, 2)))
        size(coords, 1) == 3 || throw(ArgumentError("coords must be 3 × nnodes"))
        size(segs, 1) == 2 || throw(ArgumentError("segs must be 2 × nsegs"))
        size(tris, 1) == 3 || throw(ArgumentError("tris must be 3 × ntris"))
        size(tets, 1) == 4 || throw(ArgumentError("tets must be 4 × ntets"))
        length(seg_tag) == size(segs, 2) || throw(ArgumentError("seg_tag length mismatch"))
        length(tri_tag) == size(tris, 2) || throw(ArgumentError("tri_tag length mismatch"))
        length(tet_tag) == size(tets, 2) || throw(ArgumentError("tet_tag length mismatch"))
        nn = size(coords, 2)
        _check_ids(segs, nn, "segment"); _check_ids(tris, nn, "triangle"); _check_ids(tets, nn, "tet")
        new(Matrix{Float64}(coords), Matrix{Int32}(segs), Matrix{Int32}(tris),
            Matrix{Int32}(tets), Vector{Int32}(seg_tag), Vector{Int32}(tri_tag),
            Vector{Int32}(tet_tag))
    end
end

function _check_ids(cells::AbstractMatrix, nn::Integer, what::AbstractString)
    @inbounds for j in axes(cells, 2), i in axes(cells, 1)
        v = cells[i, j]
        (1 <= v <= nn) || throw(ArgumentError("$what $j references node $v outside 1:$nn"))
    end
end

@inline nnodes(m::Mesh) = size(m.coords, 2)
@inline nsegs(m::Mesh)  = size(m.segs, 2)
@inline ntris(m::Mesh)  = size(m.tris, 2)
@inline ntets(m::Mesh)  = size(m.tets, 2)

"""`node(m, i)` → the `(x,y,z)` of node `i` as an allocation-free `NTuple{3}`."""
@inline function node(m::Mesh, i::Integer)
    @inbounds (m.coords[1, i], m.coords[2, i], m.coords[3, i])
end

"""`node2(m, i)` → the `(x,y)` of node `i` (planar), allocation-free."""
@inline function node2(m::Mesh, i::Integer)
    @inbounds (m.coords[1, i], m.coords[2, i])
end

# ════════════════════════════════════════════════════════════════════════════════
# Geometry / quality (allocation-free, tuple in — scalar out)
# ════════════════════════════════════════════════════════════════════════════════

@inline _sub(a, b) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline _dot(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline _cross(a, b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline _norm(a) = sqrt(_dot(a, a))
@inline _scale(a, s) = (a[1]*s, a[2]*s, a[3]*s)
@inline _add(a, b) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])

"""Unsigned area of triangle `(a,b,c)` (3-D points)."""
@inline function triangle_area(a, b, c)
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))
end

"""Signed volume of tet `(a,b,c,d)` = det[b-a, c-a, d-a]/6 (sign = orientation)."""
@inline function tet_signed_volume(a, b, c, d)
    return _dot(_sub(b, a), _cross(_sub(c, a), _sub(d, a))) / 6.0
end

@inline tet_volume(a, b, c, d) = abs(tet_signed_volume(a, b, c, d))

"""
    tet_dihedral_extrema(a,b,c,d) -> (min_angle, max_angle) in radians

Min and max of the six dihedral angles of the tet. Slivers show a min near 0 or
a max near π. Computed from edge-orthogonal in-face vectors (numerically stable).
"""
function tet_dihedral_extrema(a, b, c, d)
    P = (a, b, c, d)
    # edges as vertex-index pairs and the two "apex" vertices of each edge
    edges = ((1,2,3,4), (1,3,2,4), (1,4,2,3), (2,3,1,4), (2,4,1,3), (3,4,1,2))
    mn = Inf; mx = -Inf
    @inbounds for e in edges
        i, j, k, l = e
        eij = _sub(P[j], P[i])
        len = _norm(eij)
        len == 0 && return (0.0, Float64(pi))
        u = _scale(eij, 1/len)                 # unit edge direction
        vk = _sub(P[k], P[i]); vl = _sub(P[l], P[i])
        nk = _sub(vk, _scale(u, _dot(vk, u)))  # component of vk ⊥ edge
        nl = _sub(vl, _scale(u, _dot(vl, u)))
        dk = _norm(nk); dl = _norm(nl)
        (dk == 0 || dl == 0) && return (0.0, Float64(pi))
        cosang = clamp(_dot(nk, nl) / (dk*dl), -1.0, 1.0)
        ang = acos(cosang)
        ang < mn && (mn = ang)
        ang > mx && (mx = ang)
    end
    return (mn, mx)
end

"""Circumradius of tet `(a,b,c,d)`. Returns `Inf` for a degenerate (flat) tet."""
function tet_circumradius(a, b, c, d)
    A = _sub(b, a); B = _sub(c, a); C = _sub(d, a)
    denom = 2.0 * _dot(A, _cross(B, C))       # 12·signed volume
    denom == 0 && return Inf
    la = _dot(A, A); lb = _dot(B, B); lc = _dot(C, C)
    O = _add(_add(_scale(_cross(B, C), la), _scale(_cross(C, A), lb)),
             _scale(_cross(A, B), lc))
    O = _scale(O, 1/denom)                     # circumcenter relative to a
    return _norm(O)
end

"""Radius-edge ratio = circumradius / shortest edge (Delaunay-refinement quality)."""
function tet_radius_edge(a, b, c, d)
    R = tet_circumradius(a, b, c, d)
    P = (a, b, c, d)
    mn = Inf
    @inbounds for i in 1:4, j in i+1:4
        l = _norm(_sub(P[i], P[j])); l < mn && (mn = l)
    end
    mn == 0 && return Inf
    return R / mn
end

# ════════════════════════════════════════════════════════════════════════════════
# Topology: facet incidence, boundary extraction, Euler, manifoldness
# ════════════════════════════════════════════════════════════════════════════════

@inline _sort3(a, b, c) = a <= b ? (a <= c ? (b <= c ? (a,b,c) : (a,c,b)) : (c,a,b)) :
                                    (b <= c ? (a <= c ? (b,a,c) : (b,c,a)) : (c,b,a))
@inline _sort2(a, b) = a <= b ? (a, b) : (b, a)

"""
    boundary_faces(tets) -> Vector{NTuple{3,Int32}}, max_incidence

Triangular faces incident to exactly one tet (the volume boundary), each as a
sorted vertex triple. Also returns the maximum face incidence: `>2` ⇒ the tet
complex is non-manifold.
"""
function boundary_faces(tets::AbstractMatrix{<:Integer})
    inc = Dict{NTuple{3,Int32}, Int}()
    sizehint!(inc, 2 * size(tets, 2))
    @inbounds for t in axes(tets, 2)
        v1 = Int32(tets[1,t]); v2 = Int32(tets[2,t]); v3 = Int32(tets[3,t]); v4 = Int32(tets[4,t])
        for f in (_sort3(v2,v3,v4), _sort3(v1,v3,v4), _sort3(v1,v2,v4), _sort3(v1,v2,v3))
            inc[f] = get(inc, f, 0) + 1
        end
    end
    maxinc = 0
    bnd = NTuple{3,Int32}[]
    for (f, c) in inc
        c > maxinc && (maxinc = c)
        c == 1 && push!(bnd, f)
    end
    return bnd, maxinc
end

"""
    boundary_edges(tris) -> Vector{NTuple{2,Int32}}, max_incidence

Edges incident to exactly one triangle (the surface boundary), each a sorted
vertex pair; plus the max edge incidence (`>2` ⇒ non-manifold surface).
"""
function boundary_edges(tris::AbstractMatrix{<:Integer})
    inc = Dict{NTuple{2,Int32}, Int}()
    sizehint!(inc, 2 * size(tris, 2))
    @inbounds for t in axes(tris, 2)
        v1 = Int32(tris[1,t]); v2 = Int32(tris[2,t]); v3 = Int32(tris[3,t])
        for e in (_sort2(v1,v2), _sort2(v2,v3), _sort2(v1,v3))
            inc[e] = get(inc, e, 0) + 1
        end
    end
    maxinc = 0
    bnd = NTuple{2,Int32}[]
    for (e, c) in inc
        c > maxinc && (maxinc = c)
        c == 1 && push!(bnd, e)
    end
    return bnd, maxinc
end

"""Set of unique undirected edges over triangles+tets (for Euler characteristic)."""
function unique_edges(tris::AbstractMatrix{<:Integer}, tets::AbstractMatrix{<:Integer})
    edges = Set{NTuple{2,Int32}}()
    @inbounds for t in axes(tris, 2)
        v1=Int32(tris[1,t]); v2=Int32(tris[2,t]); v3=Int32(tris[3,t])
        push!(edges, _sort2(v1,v2)); push!(edges, _sort2(v2,v3)); push!(edges, _sort2(v1,v3))
    end
    @inbounds for t in axes(tets, 2)
        v1=Int32(tets[1,t]); v2=Int32(tets[2,t]); v3=Int32(tets[3,t]); v4=Int32(tets[4,t])
        push!(edges, _sort2(v1,v2)); push!(edges, _sort2(v1,v3)); push!(edges, _sort2(v1,v4))
        push!(edges, _sort2(v2,v3)); push!(edges, _sort2(v2,v4)); push!(edges, _sort2(v3,v4))
    end
    return edges
end

"""Set of unique triangular faces over triangles+tets (for Euler characteristic)."""
function unique_faces(tris::AbstractMatrix{<:Integer}, tets::AbstractMatrix{<:Integer})
    faces = Set{NTuple{3,Int32}}()
    @inbounds for t in axes(tris, 2)
        push!(faces, _sort3(Int32(tris[1,t]), Int32(tris[2,t]), Int32(tris[3,t])))
    end
    @inbounds for t in axes(tets, 2)
        v1=Int32(tets[1,t]); v2=Int32(tets[2,t]); v3=Int32(tets[3,t]); v4=Int32(tets[4,t])
        push!(faces, _sort3(v2,v3,v4)); push!(faces, _sort3(v1,v3,v4))
        push!(faces, _sort3(v1,v2,v4)); push!(faces, _sort3(v1,v2,v3))
    end
    return faces
end

"""
    euler_characteristic(m) -> Int

`V − E + F − T` over the mesh's triangle+tet complex (unique edges/faces). For a
triangulated 3-ball this is `1`; for a triangulated 2-sphere (tris only, closed)
it is `2`.
"""
function euler_characteristic(m::Mesh)
    # V counts nodes referenced by the triangle+tet complex (isolated coords
    # excluded). Segment nodes are deliberately NOT counted: E/F/T span only
    # tris+tets, so adding segment vertices to V without their 1-cells to E would
    # inflate V − E + F − T (a lone segment would read χ=2 instead of 1).
    used = falses(nnodes(m))
    @inbounds for t in axes(m.tris, 2), i in 1:3; used[m.tris[i,t]] = true; end
    @inbounds for t in axes(m.tets, 2), i in 1:4; used[m.tets[i,t]] = true; end
    V = count(used)
    E = length(unique_edges(m.tris, m.tets))
    F = length(unique_faces(m.tris, m.tets))
    T = ntets(m)
    return V - E + F - T
end

"""Euler characteristic of the extracted boundary surface (`V − E + F`)."""
function boundary_euler(m::Mesh)
    bf, _ = boundary_faces(m.tets)
    isempty(bf) && return 0
    verts = Set{Int32}()
    edges = Set{NTuple{2,Int32}}()
    for f in bf
        push!(verts, f[1]); push!(verts, f[2]); push!(verts, f[3])
        push!(edges, _sort2(f[1],f[2])); push!(edges, _sort2(f[2],f[3])); push!(edges, _sort2(f[1],f[3]))
    end
    return length(verts) - length(edges) + length(bf)
end

"""
    is_closed_manifold(m) -> Bool

True iff the nonempty tet complex is a combinatorial 3-manifold with boundary.
Besides face incidence, every vertex link is certified as a connected triangulated
sphere (interior vertex) or disk (boundary vertex); this rejects cells that touch
only at an edge/vertex and other pinched complexes that a face-count test misses.
"""
function is_closed_manifold(m::Mesh)
    ok, _ = _tet_topology(m.tets)
    return ok
end

# A tetrahedral pseudocomplex can have face incidence <= 2 yet fail to be a
# 3-manifold (two solids touching only at one edge or vertex are the common case).
# The definitive local criterion is the link of every vertex: it must be a connected
# 2-manifold homeomorphic to S² (interior) or D² (boundary).  For each link we check
# edge incidence, triangle connectivity, boundary degree, and Euler characteristic.
# A link edge is the original face opposite one of its vertices, so these incidence
# checks also certify every tet face without a redundant global face table.  Global
# records are flat, compact arrays; per-link scratch is reused, avoiding one
# Set/Vector allocation per mesh vertex.
@inline function _uf_find!(p::Vector{Int32}, i::Int32)
    while @inbounds(p[i] != i)
        @inbounds p[i] = p[p[i]]
        @inbounds i = p[i]
    end
    return i
end
@inline function _uf_union!(p::Vector{Int32}, a::Int32, b::Int32)
    ra = _uf_find!(p, a); rb = _uf_find!(p, b)
    ra != rb && (@inbounds p[ra] = rb)
    return nothing
end
@inline function _nunique_sorted(v)
    isempty(v) && return 0
    n = 1
    @inbounds for i in 2:length(v); v[i] != v[i-1] && (n += 1); end
    return n
end

function _tet_topology(tets::AbstractMatrix{<:Integer})
    nt = size(tets, 2)
    nt == 0 && return (false, "tet complex is empty")
    nt <= typemax(Int32) ||
        return (false, "tet count $nt exceeds the Int32 topology-audit limit")
    n4 = try
        Base.checked_mul(4, nt)
    catch
        return (false, "tet topology record count overflow")
    end

    vertex_tets = Vector{NTuple{2,Int32}}(undef, n4)
    tetkeys = Vector{NTuple{4,Int32}}(undef, nt)
    vi = 0
    @inbounds for t in 1:nt
        a=Int32(tets[1,t]); b=Int32(tets[2,t]); c=Int32(tets[3,t]); d=Int32(tets[4,t])
        (a==b || a==c || a==d || b==c || b==d || c==d) &&
            return (false, "tet $t repeats a vertex")
        tetkeys[t] = _sort4(a,b,c,d)
        ti = Int32(t)
        for v in (a,b,c,d); vi += 1; vertex_tets[vi] = (v,ti); end
    end

    sort!(tetkeys; alg=QuickSort)
    @inbounds for i in 2:nt
        tetkeys[i] == tetkeys[i-1] && return (false, "duplicate tetrahedron")
    end

    sort!(vertex_tets; alg=QuickSort)
    linkverts = Int32[]
    linkedges = NTuple{3,Int32}[]       # (link edge endpoints, local triangle id)
    boundaryverts = Int32[]             # endpoints of incidence-1 link edges
    parent = Int32[]
    hasboundary = false
    i = 1
    while i <= n4
        j = i + 1
        @inbounds while j <= n4 && vertex_tets[j][1] == vertex_tets[i][1]; j += 1; end
        v = @inbounds vertex_tets[i][1]
        degree = j - i
        empty!(linkverts); empty!(linkedges); empty!(boundaryverts)
        resize!(parent, degree)
        @inbounds for k in 1:degree; parent[k] = Int32(k); end

        @inbounds for k in i:j-1
            t = Int(vertex_tets[k][2]); localid = Int32(k-i+1)
            a=Int32(tets[1,t]); b=Int32(tets[2,t]); c=Int32(tets[3,t]); d=Int32(tets[4,t])
            x,y,z = v==a ? (b,c,d) : (v==b ? (a,c,d) : (v==c ? (a,b,d) : (a,b,c)))
            push!(linkverts, x, y, z)
            e1=_sort2(x,y); e2=_sort2(y,z); e3=_sort2(x,z)
            push!(linkedges, (e1[1],e1[2],localid),
                             (e2[1],e2[2],localid),
                             (e3[1],e3[2],localid))
        end

        sort!(linkedges; alg=QuickSort)
        nedges = 0; k = 1
        @inbounds while k <= length(linkedges)
            l = k + 1
            while l <= length(linkedges) &&
                  linkedges[l][1] == linkedges[k][1] && linkedges[l][2] == linkedges[k][2]
                l += 1
            end
            incidence = l-k; nedges += 1
            incidence > 2 && return (false, "vertex $v has a non-manifold link edge")
            if incidence == 2
                _uf_union!(parent, linkedges[k][3], linkedges[k+1][3])
            else
                hasboundary = true
                push!(boundaryverts, linkedges[k][1], linkedges[k][2])
            end
            k = l
        end

        root = _uf_find!(parent, Int32(1))
        @inbounds for k in 2:degree
            _uf_find!(parent, Int32(k)) == root ||
                return (false, "vertex $v has a disconnected link")
        end

        sort!(linkverts; alg=QuickSort)
        nlinkverts = _nunique_sorted(linkverts)
        if !isempty(boundaryverts)
            sort!(boundaryverts; alg=QuickSort)
            k = 1
            @inbounds while k <= length(boundaryverts)
                l = k + 1
                while l <= length(boundaryverts) && boundaryverts[l] == boundaryverts[k]; l += 1; end
                l-k == 2 || return (false, "vertex $v has a non-circular link boundary")
                k = l
            end
        end
        chi = nlinkverts - nedges + degree
        expected = isempty(boundaryverts) ? 2 : 1
        chi == expected ||
            return (false, "vertex $v link has Euler characteristic $chi (expected $expected)")
        i = j
    end
    hasboundary || return (false, "tet complex has an empty boundary")
    return (true, "ok")
end

# ════════════════════════════════════════════════════════════════════════════════
# Bounding box, validation, checksum
# ════════════════════════════════════════════════════════════════════════════════

"""Axis-aligned bounding box `(lo, hi)` over referenced nodes (isolated coords
excluded, matching [`euler_characteristic`](@ref)). Since `mesh_crc` embeds this
box, an unreferenced far-away node must not be allowed to enlarge it."""
function bounding_box(m::Mesh)
    nnodes(m) == 0 && return ((0.0,0.0,0.0), (0.0,0.0,0.0))
    used = falses(nnodes(m))
    @inbounds for t in axes(m.segs, 2), i in 1:2; used[m.segs[i,t]] = true; end
    @inbounds for t in axes(m.tris, 2), i in 1:3; used[m.tris[i,t]] = true; end
    @inbounds for t in axes(m.tets, 2), i in 1:4; used[m.tets[i,t]] = true; end
    lo = (Inf, Inf, Inf); hi = (-Inf, -Inf, -Inf)
    any_used = false
    @inbounds for i in 1:nnodes(m)
        used[i] || continue
        any_used = true
        p = node(m, i)
        lo = (min(lo[1],p[1]), min(lo[2],p[2]), min(lo[3],p[3]))
        hi = (max(hi[1],p[1]), max(hi[2],p[2]), max(hi[3],p[3]))
    end
    any_used || return ((0.0,0.0,0.0), (0.0,0.0,0.0))
    return lo, hi
end

struct MeshDiagnostic
    ok::Bool
    messages::Vector{String}
end
Base.show(io::IO, d::MeshDiagnostic) =
    print(io, "MeshDiagnostic(ok=$(d.ok)" * (isempty(d.messages) ? "" : ", " * join(d.messages, "; ")) * ")")

"""
    validate(m; require_positive_tets=true) -> MeshDiagnostic

Completeness contract (DEVELOPMENT.md): a mesh is either validated or carries a
precise diagnostic — never a silent defect. Checks index bounds (already enforced
at construction), non-degenerate cells, positive tet volumes, finite coordinates,
and non-manifold faces. Empty-region detection is the caller's job (it must count
tets *per region*); `validate` reports the global tet count so "no empty volumes"
can never be vacuously true.
"""
function validate(m::Mesh; require_positive_tets::Bool=true)
    msgs = String[]
    # finite coordinates
    @inbounds for i in 1:nnodes(m)
        p = node(m, i)
        (isfinite(p[1]) && isfinite(p[2]) && isfinite(p[3])) ||
            (push!(msgs, "node $i has non-finite coordinate"); break)
    end
    # degenerate / inverted tets
    nneg = 0; nflat = 0; nbadvol = 0
    @inbounds for t in 1:ntets(m)
        a = node(m, m.tets[1,t]); b = node(m, m.tets[2,t])
        c = node(m, m.tets[3,t]); d = node(m, m.tets[4,t])
        v = tet_signed_volume(a, b, c, d)
        if !isfinite(v); nbadvol += 1
        elseif v == 0; nflat += 1
        elseif v < 0; nneg += 1
        end
    end
    nbadvol > 0 && push!(msgs, "$nbadvol tets have non-finite computed volume")
    nflat > 0 && push!(msgs, "$nflat degenerate (zero-volume) tets")
    if require_positive_tets && nneg > 0
        push!(msgs, "$nneg negatively-oriented tets")
    end
    # degenerate triangles
    ndegtri = 0; nbadtri = 0
    @inbounds for t in 1:ntris(m)
        a = node(m, m.tris[1,t]); b = node(m, m.tris[2,t]); c = node(m, m.tris[3,t])
        ar = triangle_area(a, b, c)
        !isfinite(ar) ? (nbadtri += 1) : (ar == 0 && (ndegtri += 1))
    end
    nbadtri > 0 && push!(msgs, "$nbadtri triangles have non-finite computed area")
    ndegtri > 0 && push!(msgs, "$ndegtri degenerate (zero-area) triangles")
    # Full combinatorial-manifold audit (faces alone miss edge/vertex pinches).
    if ntets(m) > 0
        topok, reason = _tet_topology(m.tets)
        topok || push!(msgs, "non-manifold tet complex: $reason")
    end
    return MeshDiagnostic(isempty(msgs), msgs)
end

"""
    mesh_crc(m) -> NamedTuple

Deterministic regression checksum (DEVELOPMENT.md "Mesh-CRC checksum"):
`(n_nodes, n_segs, n_tris, n_tets, bbox, dihedral (min,mean), radius-edge
(min,mean), n_boundary_faces, sha)` where `sha` is the SHA-256 of the mesh's
*canonically sorted* connectivity (invariant to per-cell vertex order and cell
order — so it tracks topology, not incidental labelling).
"""
function mesh_crc(m::Mesh)
    lo, hi = bounding_box(m)
    # tet quality aggregates
    dmin = Inf; dsum = 0.0; remin = Inf; resum = 0.0; nt = ntets(m)
    @inbounds for t in 1:nt
        a = node(m, m.tets[1,t]); b = node(m, m.tets[2,t])
        c = node(m, m.tets[3,t]); d = node(m, m.tets[4,t])
        dmn, _ = tet_dihedral_extrema(a, b, c, d)
        dmn < dmin && (dmin = dmn); dsum += dmn
        re = tet_radius_edge(a, b, c, d)
        re < remin && (remin = re)
        resum += re
    end
    dihedral = nt == 0 ? (0.0, 0.0) : (dmin, dsum/nt)
    radedge  = nt == 0 ? (0.0, 0.0) : (remin, resum/nt)
    nbf = nt == 0 ? 0 : length(first(boundary_faces(m.tets)))

    return (n_nodes = nnodes(m), n_segs = nsegs(m), n_tris = ntris(m), n_tets = nt,
            bbox = (lo, hi), dihedral = dihedral, radius_edge = radedge,
            n_boundary_faces = nbf, sha = _connectivity_sha(m))
end

"""SHA-256 (hex) of canonicalized connectivity: each cell's vertices sorted, then
cells sorted lexicographically, per cell type, hashed as little-endian Int32."""
function _connectivity_sha(m::Mesh)
    ctx = SHA.SHA2_256_CTX()
    _hash_cells!(ctx, m.segs, 2)
    _hash_cells!(ctx, m.tris, 3)
    _hash_cells!(ctx, m.tets, 4)
    return bytes2hex(SHA.digest!(ctx))
end

function _hash_cells!(ctx, cells::AbstractMatrix{<:Integer}, k::Int)
    n = size(cells, 2)
    rows = Vector{NTuple{4,Int32}}(undef, n)   # pad to 4 for a uniform key type
    @inbounds for j in 1:n
        a = Int32(cells[1,j]); b = Int32(cells[2,j])
        if k == 2
            (a, b) = _sort2(a, b); rows[j] = (a, b, Int32(0), Int32(0))
        elseif k == 3
            c = Int32(cells[3,j]); (a, b, c) = _sort3(a, b, c); rows[j] = (a, b, c, Int32(0))
        else
            c = Int32(cells[3,j]); d = Int32(cells[4,j])
            s = _sort4(a, b, c, d); rows[j] = s
        end
    end
    sort!(rows)
    buf = Vector{UInt8}(undef, 16)
    @inbounds for r in rows
        _put_i32!(buf, 1, r[1]); _put_i32!(buf, 5, r[2])
        _put_i32!(buf, 9, r[3]); _put_i32!(buf, 13, r[4])
        SHA.update!(ctx, buf)
    end
    # length tag so [] vs [] of different type can't collide
    tag = UInt8[k % UInt8, (n % 256) % UInt8, ((n >> 8) % 256) % UInt8, ((n >> 16) % 256) % UInt8]
    SHA.update!(ctx, tag)
end

@inline function _sort4(a, b, c, d)
    # sort 4 Int32 ascending (branchy but tiny/allocation-free)
    if a > b; a, b = b, a; end
    if c > d; c, d = d, c; end
    if a > c; a, c = c, a; end
    if b > d; b, d = d, b; end
    if b > c; b, c = c, b; end
    return (a, b, c, d)
end

@inline function _put_i32!(buf, off, v::Int32)
    u = reinterpret(UInt32, v)
    @inbounds buf[off]   = u % UInt8
    @inbounds buf[off+1] = (u >> 8)  % UInt8
    @inbounds buf[off+2] = (u >> 16) % UInt8
    @inbounds buf[off+3] = (u >> 24) % UInt8
    return nothing
end

"""One-line human-readable form of a `mesh_crc` result for STATUS.md logging."""
function crc_string(crc)
    lo, hi = crc.bbox
    return string(
        "nodes=", crc.n_nodes, " segs=", crc.n_segs, " tris=", crc.n_tris,
        " tets=", crc.n_tets,
        " bbox=[", round.(lo; sigdigits=6), "→", round.(hi; sigdigits=6), "]",
        " dihedral(min,mean)=", round.(crc.dihedral; sigdigits=5),
        " radedge(min,mean)=", round.(crc.radius_edge; sigdigits=5),
        " bfaces=", crc.n_boundary_faces,
        " sha=", crc.sha[1:16], "…")
end

end # module MeshTypes
