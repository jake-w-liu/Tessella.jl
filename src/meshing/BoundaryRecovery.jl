# Boundary-conforming tetrahedralization of a closed triangulated surface.
#
# The surface vertices are tetrahedralized inside an enclosing cube. Every
# surface edge and then every surface triangle missing from that
# tetrahedralization is recovered locally: by flips (2-3 flips and edge
# removals) where they suffice, otherwise by removing the cells the constraint
# crosses and refilling that cavity through an exhaustive search. Finally the
# cells outside the surface are discarded. A refill may add a vertex strictly
# inside its cavity; no vertex is ever placed on the surface, so the result
# always has exactly the input triangles as its boundary.
#
# Robustness rests on two invariants that every mutation re-establishes with
# exact orientation tests on the input coordinates:
#   * every live cell is strictly positively oriented, and
#   * a local replacement has the same oriented boundary as the cells it removes.
# Together they imply that the cells always tile the enclosing cube: a
# positively oriented chain with the boundary of a ball covers it exactly
# once. Nothing else is needed for validity, so an approximate in-sphere test
# may guide the insertion and no step tests cells for intersection. A recovery
# that cannot proceed returns a diagnostic instead of a mesh.
module BoundaryRecovery

using ..Predicates: orient3
using ..MeshTypes: Mesh

export recover_boundary_flips

const _O3D_ERRBOUND = (7.0 + 56.0 * 2.0^-53) * 2.0^-53
# Slots of the face opposite each vertex slot, ordered so that the opposite
# vertex lies on the positive side.
const _FV = ((2, 4, 3), (1, 3, 4), (1, 4, 2), (1, 2, 3))
const _EDGES = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
const _MAX_RING = 16

mutable struct Triangulation
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    tv::Vector{Int32}          # four vertices per cell
    tn::Vector{Int32}          # neighbour across the face opposite each slot (0: cube boundary)
    alive::Vector{Bool}
    free::Vector{Int32}
    v2t::Vector{Int32}         # one live cell incident to each vertex
    mark::Vector{UInt32}
    stamp::UInt32
    vmark::Vector{UInt32}
    vstamp::UInt32
    rng::UInt32
    cav::Vector{Int32}
    queue::Vector{Int32}
    bnd::Vector{NTuple{5,Int32}}
    ekeys::Vector{UInt64}
    evals::Vector{Int32}
    egens::Vector{UInt32}
    egen::UInt32
end

@inline _point(T::Triangulation, v::Integer) = @inbounds (T.x[v], T.y[v], T.z[v])

# Sign of the signed volume of (a,b,c,d): positive for a valid cell.
@inline function _vol(T::Triangulation, a::Integer, b::Integer, c::Integer, d::Integer)::Int
    x = T.x; y = T.y; z = T.z
    @inbounds begin
        dx = x[d]; dy = y[d]; dz = z[d]
        adx = x[a] - dx; ady = y[a] - dy; adz = z[a] - dz
        bdx = x[b] - dx; bdy = y[b] - dy; bdz = z[b] - dz
        cdx = x[c] - dx; cdy = y[c] - dy; cdz = z[c] - dz
    end
    bdxcdy = bdx * cdy; cdxbdy = cdx * bdy
    cdxady = cdx * ady; adxcdy = adx * cdy
    adxbdy = adx * bdy; bdxady = bdx * ady
    det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady)
    permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz) +
                (abs(cdxady) + abs(adxcdy)) * abs(bdz) +
                (abs(adxbdy) + abs(bdxady)) * abs(cdz)
    # The static filter is only trusted away from underflow and overflow.
    if 1.0e-250 < permanent < 1.0e250
        bound = _O3D_ERRBOUND * permanent
        det > bound && return -1
        det < -bound && return 1
    end
    return -orient3(_point(T, a), _point(T, b), _point(T, c), _point(T, d))
end

# Approximate in-sphere test for a positively oriented cell. It only steers
# which cells an insertion replaces; validity never depends on its answer.
@inline function _insphere_approx(T::Triangulation, t::Integer, p::Integer)::Bool
    x = T.x; y = T.y; z = T.z; base = 4 * (t - 1)
    @inbounds begin
        a = T.tv[base+1]; b = T.tv[base+2]; c = T.tv[base+3]; d = T.tv[base+4]
        ex = x[p]; ey = y[p]; ez = z[p]
        aex = x[a] - ex; aey = y[a] - ey; aez = z[a] - ez
        bex = x[b] - ex; bey = y[b] - ey; bez = z[b] - ez
        cex = x[c] - ex; cey = y[c] - ey; cez = z[c] - ez
        dex = x[d] - ex; dey = y[d] - ey; dez = z[d] - ez
    end
    ab = aex * bey - bex * aey; bc = bex * cey - cex * bey
    cd = cex * dey - dex * cey; da = dex * aey - aex * dey
    ac = aex * cey - cex * aey; bd = bex * dey - dex * bey
    abc = aez * bc - bez * ac + cez * ab
    bcd = bez * cd - cez * bd + dez * bc
    cda = cez * da + dez * ac + aez * cd
    dab = dez * ab + aez * bd + bez * da
    alift = aex * aex + aey * aey + aez * aez
    blift = bex * bex + bey * bey + bez * bez
    clift = cex * cex + cey * cey + cez * cez
    dlift = dex * dex + dey * dey + dez * dez
    return (dlift * abc - clift * dab) + (blift * cda - alift * bcd) < 0.0
end

# Whether cell (a,b,c,d) is flat to within rounding: its volume is negligible
# against the cube of its longest edge. Vertices that are coplanar only up to
# rounding (on a ruled or tilted planar face, say) form such cells; they are
# valid, but too thin to place a point in, so a refill never adds one and a
# cavity that contains one is widened around it.
@inline function _flat(T::Triangulation, a::Integer, b::Integer, c::Integer, d::Integer,
                       threshold::Float64=1.0e-10)::Bool
    x = T.x; y = T.y; z = T.z
    @inbounds begin
        ux = x[b] - x[a]; uy = y[b] - y[a]; uz = z[b] - z[a]
        vx = x[c] - x[a]; vy = y[c] - y[a]; vz = z[c] - z[a]
        wx = x[d] - x[a]; wy = y[d] - y[a]; wz = z[d] - z[a]
    end
    volume = ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)
    longest = max(ux * ux + uy * uy + uz * uz, vx * vx + vy * vy + vz * vz,
                  wx * wx + wy * wy + wz * wz,
                  (vx - ux)^2 + (vy - uy)^2 + (vz - uz)^2,
                  (wx - ux)^2 + (wy - uy)^2 + (wz - uz)^2,
                  (wx - vx)^2 + (wy - vy)^2 + (wz - vz)^2)
    return abs(volume) <= threshold * longest * sqrt(longest)
end

@inline _tv(T::Triangulation, t::Integer, i::Integer) = @inbounds T.tv[4*(t-1)+i]
@inline _tn(T::Triangulation, t::Integer, i::Integer) = @inbounds T.tn[4*(t-1)+i]

@inline function _face(T::Triangulation, t::Integer, s::Integer)
    base = 4 * (t - 1); f = _FV[s]
    return @inbounds (T.tv[base+f[1]], T.tv[base+f[2]], T.tv[base+f[3]])
end

@inline function _sort3(a::Int32, b::Int32, c::Int32)
    a > b && ((a, b) = (b, a))
    b > c && ((b, c) = (c, b))
    a > b && ((a, b) = (b, a))
    return (a, b, c)
end
@inline _sort3(f::NTuple{3,Int32}) = _sort3(f[1], f[2], f[3])
@inline _ekey(a::Integer, b::Integer) =
    a < b ? (UInt64(a) << 32) | UInt64(b) : (UInt64(b) << 32) | UInt64(a)

@inline function _slot(T::Triangulation, t::Integer, v::Integer)::Int
    base = 4 * (t - 1)
    @inbounds for i in 1:4
        T.tv[base+i] == v && return i
    end
    return 0
end

@inline function _neighbor_slot(T::Triangulation, t::Integer, other::Integer)::Int
    base = 4 * (t - 1)
    @inbounds for i in 1:4
        T.tn[base+i] == other && return i
    end
    return 0
end

@inline function _rand!(T::Triangulation)
    T.rng = T.rng * UInt32(1664525) + UInt32(1013904223)
    return T.rng >> 16
end

function _newcell!(T::Triangulation, a::Int32, b::Int32, c::Int32, d::Int32)::Int32
    if isempty(T.free)
        push!(T.tv, a, b, c, d); push!(T.tn, 0, 0, 0, 0)
        push!(T.alive, true); push!(T.mark, UInt32(0))
        t = Int32(length(T.alive))
    else
        t = pop!(T.free); base = 4 * (t - 1)
        @inbounds begin
            T.tv[base+1] = a; T.tv[base+2] = b; T.tv[base+3] = c; T.tv[base+4] = d
            T.tn[base+1] = 0; T.tn[base+2] = 0; T.tn[base+3] = 0; T.tn[base+4] = 0
            T.alive[t] = true; T.mark[t] = 0
        end
    end
    @inbounds begin
        T.v2t[a] = t; T.v2t[b] = t; T.v2t[c] = t; T.v2t[d] = t
    end
    return t
end

@inline function _kill!(T::Triangulation, t::Integer)
    @inbounds T.alive[t] = false
    push!(T.free, Int32(t))
    return nothing
end

# ── small open-addressed edge table used to stitch a refilled cavity ─────────

function _edges_reset!(T::Triangulation, capacity::Int)
    size = 64
    while size < 2 * capacity
        size <<= 1
    end
    if length(T.ekeys) < size
        resize!(T.ekeys, size); resize!(T.evals, size)
        resize!(T.egens, size); fill!(T.egens, UInt32(0)); T.egen = UInt32(0)
    end
    T.egen += UInt32(1)
    if T.egen == typemax(UInt32)
        fill!(T.egens, UInt32(0)); T.egen = UInt32(1)
    end
    return nothing
end

# Index of `key`, and whether it was already present.
@inline function _edges_find!(T::Triangulation, key::UInt64)
    mask = length(T.ekeys) - 1
    index = Int((key * 0x9e3779b97f4a7c15) >> 40) & mask
    @inbounds while true
        if T.egens[index+1] != T.egen
            T.egens[index+1] = T.egen; T.ekeys[index+1] = key
            return index + 1, false
        end
        T.ekeys[index+1] == key && return index + 1, true
        index = (index + 1) & mask
    end
end

# ── construction ────────────────────────────────────────────────────────────

function Triangulation(coords::AbstractMatrix{Float64}, count::Int)
    total = count + 8
    x = Vector{Float64}(undef, total); y = similar(x); z = similar(x)
    lo = (Inf, Inf, Inf); hi = (-Inf, -Inf, -Inf)
    @inbounds for i in 1:count
        px = coords[1, i]; py = coords[2, i]; pz = coords[3, i]
        x[i] = px; y[i] = py; z[i] = pz
        lo = (min(lo[1], px), min(lo[2], py), min(lo[3], pz))
        hi = (max(hi[1], px), max(hi[2], py), max(hi[3], pz))
    end
    extent = max(hi[1] - lo[1], hi[2] - lo[2], hi[3] - lo[3])
    (isfinite(extent) && extent > 0) || throw(ArgumentError(
        "boundary recovery: surface vertices have no finite positive extent"))
    pad = 0.3 * extent
    lo = (lo[1] - pad, lo[2] - pad, lo[3] - pad)
    hi = (hi[1] + pad, hi[2] + pad, hi[3] + pad)
    corners = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    for (k, corner) in enumerate(corners)
        x[count+k] = corner[1] == 1 ? hi[1] : lo[1]
        y[count+k] = corner[2] == 1 ? hi[2] : lo[2]
        z[count+k] = corner[3] == 1 ? hi[3] : lo[3]
    end
    T = Triangulation(x, y, z, Int32[], Int32[], Bool[], Int32[], zeros(Int32, total),
                      UInt32[], UInt32(0), zeros(UInt32, total), UInt32(0), UInt32(12345),
                      Int32[], Int32[], NTuple{5,Int32}[], UInt64[], Int32[], UInt32[], UInt32(0))
    sizehint!(T.tv, 32 * total); sizehint!(T.tn, 32 * total)
    sizehint!(T.alive, 8 * total); sizehint!(T.mark, 8 * total)
    for row in ((8, 3, 4, 2), (8, 1, 2, 4), (2, 7, 8, 3), (2, 1, 8, 5), (5, 2, 6, 8), (8, 2, 6, 7))
        a = Int32(count + row[1]); b = Int32(count + row[2])
        c = Int32(count + row[3]); d = Int32(count + row[4])
        sign = _vol(T, a, b, c, d)
        sign == 0 && throw(ArgumentError("boundary recovery: degenerate enclosing cube"))
        sign < 0 && ((a, b) = (b, a))
        _newcell!(T, a, b, c, d)
    end
    faces = Dict{NTuple{3,Int32},Tuple{Int32,Int32}}()
    for t in Int32(1):Int32(6), s in 1:4
        key = _sort3(_face(T, t, s))
        other = get(faces, key, nothing)
        if other === nothing
            faces[key] = (t, Int32(s))
        else
            T.tn[4*(t-1)+s] = other[1]
            T.tn[4*(other[1]-1)+other[2]] = t
            delete!(faces, key)
        end
    end
    return T
end

@inline function _part3(v::UInt64)
    v &= 0x00000000001fffff
    v = (v | v << 32) & 0x001f00000000ffff
    v = (v | v << 16) & 0x001f0000ff0000ff
    v = (v | v << 8) & 0x100f00f00f00f00f
    v = (v | v << 4) & 0x10c30c30c30c30c3
    v = (v | v << 2) & 0x1249249249249249
    return v
end

# Biased randomized insertion order: geometrically growing random rounds, each
# sorted along a Morton curve, so walks stay short and cavities stay small.
function _insertion_order(T::Triangulation, count::Int)
    order = collect(Int32, 1:count)
    state = UInt64(0x2545f4914f6cdd1d)
    @inbounds for i in count:-1:2
        state ⊻= state << 13; state ⊻= state >> 7; state ⊻= state << 17
        j = Int(state % UInt64(i)) + 1
        order[i], order[j] = order[j], order[i]
    end
    lo = (Inf, Inf, Inf); hi = (-Inf, -Inf, -Inf)
    @inbounds for i in 1:count
        lo = (min(lo[1], T.x[i]), min(lo[2], T.y[i]), min(lo[3], T.z[i]))
        hi = (max(hi[1], T.x[i]), max(hi[2], T.y[i]), max(hi[3], T.z[i]))
    end
    scale = ntuple(k -> hi[k] > lo[k] ? 2097151.0 / (hi[k] - lo[k]) : 0.0, 3)
    codes = Vector{UInt64}(undef, count)
    @inbounds for i in 1:count
        cx = unsafe_trunc(UInt64, (T.x[i] - lo[1]) * scale[1])
        cy = unsafe_trunc(UInt64, (T.y[i] - lo[2]) * scale[2])
        cz = unsafe_trunc(UInt64, (T.z[i] - lo[3]) * scale[3])
        codes[i] = _part3(cx) | (_part3(cy) << 1) | (_part3(cz) << 2)
    end
    stop = count
    while stop > 0
        start = stop > 256 ? stop - (stop * 7) ÷ 8 + 1 : 1
        sort!(view(order, start:stop); by=i -> @inbounds(codes[i]))
        stop = start - 1
    end
    return order
end

# Visibility walk to a cell whose closure contains vertex `p`.
function _locate(T::Triangulation, p::Int32, start::Int32)::Int32
    t = start
    limit = 4 * length(T.alive) + 64
    for _ in 1:limit
        offset = Int(_rand!(T) & 0x3)
        moved = false
        for k in 0:3
            s = ((k + offset) & 3) + 1
            a, b, c = _face(T, t, s)
            if _vol(T, a, b, c, p) < 0
                next = _tn(T, t, s)
                next == 0 && throw(ArgumentError(
                    "boundary recovery: a vertex lies outside the enclosing cube"))
                t = next; moved = true
                break
            end
        end
        moved || return t
    end
    for candidate in Int32(1):Int32(length(T.alive))
        T.alive[candidate] || continue
        inside = true
        for s in 1:4
            a, b, c = _face(T, candidate, s)
            _vol(T, a, b, c, p) < 0 && (inside = false; break)
        end
        inside && return candidate
    end
    throw(ArgumentError("boundary recovery: point location failed"))
end

# Insert vertex `p`. Returns a new cell containing it, or the negated index of
# an existing vertex with identical coordinates.
function _insert!(T::Triangulation, p::Int32, start::Int32)::Int32
    t0 = _locate(T, p, start)
    stamp = (T.stamp += UInt32(1))
    seeds = _seed_cavity!(T, p, t0, stamp)
    seeds < 0 && return Int32(seeds)
    cav = T.cav
    index = 1
    while index <= length(cav)
        t = cav[index]; index += 1
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || @inbounds(T.mark[next]) == stamp) && continue
            if _insphere_approx(T, next, p)
                @inbounds T.mark[next] = stamp
                push!(cav, next)
            end
        end
    end
    _repair_cavity!(T, p, seeds, stamp)
    return _fill_cavity!(T, p, seeds, stamp)
end

# Start the cavity of `p` with every cell whose closure contains it: that
# union is always a valid cavity, whatever the in-sphere guidance says
# afterwards. Returns the number of such cells, or the negated index of an
# existing vertex with the coordinates of `p`.
function _seed_cavity!(T::Triangulation, p::Int32, t0::Int32, stamp::UInt32)::Int
    cav = T.cav; empty!(cav)
    push!(cav, t0); @inbounds T.mark[t0] = stamp
    index = 1
    while index <= length(cav)
        t = cav[index]; index += 1
        zeros = 0
        for s in 1:4
            a, b, c = _face(T, t, s)
            if _vol(T, a, b, c, p) == 0
                zeros += 1
                next = _tn(T, t, s)
                next == 0 && throw(ArgumentError(
                    "boundary recovery: a vertex lies on the enclosing cube"))
                if @inbounds(T.mark[next]) != stamp
                    @inbounds T.mark[next] = stamp
                    push!(cav, next)
                end
            end
        end
        if zeros >= 3
            for s in 1:4
                v = _tv(T, t, s)
                (T.x[v] == T.x[p] && T.y[v] == T.y[p] && T.z[v] == T.z[p]) && return -Int(v)
            end
        end
    end
    return length(cav)
end

# Replace the marked cavity by the cells joining its boundary faces to `p`.
function _fill_cavity!(T::Triangulation, p::Int32, seeds::Int, stamp::UInt32)::Int32
    cav = T.cav
    bnd = T.bnd
    if !_cavity_boundary!(T, stamp)
        # A pinched boundary: fall back to the seed star, which is a ball.
        for k in seeds+1:length(cav)
            @inbounds T.mark[cav[k]] = 0
        end
        resize!(cav, seeds)
        _cavity_boundary!(T, stamp) || throw(ArgumentError(
            "boundary recovery: seed cavity is not a ball"))
    end
    for t in cav
        @inbounds(T.mark[t]) == stamp && _kill!(T, t)
    end
    _edges_reset!(T, 3 * length(bnd))
    last = Int32(0)
    for (a, b, c, outside, outside_slot) in bnd
        cell = _newcell!(T, a, b, c, p)
        base = 4 * (cell - 1)
        @inbounds T.tn[base+4] = outside
        outside != 0 && @inbounds(T.tn[4*(outside-1)+outside_slot] = cell)
        for (slot, key) in ((1, _ekey(b, c)), (2, _ekey(a, c)), (3, _ekey(a, b)))
            position, found = _edges_find!(T, key)
            if found
                code = @inbounds T.evals[position]
                other = code >> 2; other_slot = (code & 3) + 1
                @inbounds T.tn[base+slot] = other
                @inbounds T.tn[4*(other-1)+other_slot] = cell
            else
                @inbounds T.evals[position] = (cell << 2) | Int32(slot - 1)
            end
        end
        last = cell
    end
    return last
end

# Shrink the cavity until every boundary face is strictly visible from `p`,
# the cavity is face-connected to its seeds, and no vertex would be orphaned.
function _repair_cavity!(T::Triangulation, p::Int32, seeds::Int, stamp::UInt32;
                        removed::Bool=false)
    cav = T.cav
    while true
        changed = false
        for k in seeds+1:length(cav)
            t = cav[k]
            @inbounds(T.mark[t]) == stamp || continue
            for s in 1:4
                next = _tn(T, t, s)
                if next == 0 || @inbounds(T.mark[next]) != stamp
                    a, b, c = _face(T, t, s)
                    if _vol(T, a, b, c, p) <= 0
                        @inbounds T.mark[t] = 0
                        changed = true; removed = true
                        break
                    end
                end
            end
        end
        changed && continue
        if removed
            # Drop pieces no longer connected to the seeds: visited cells take
            # a fresh mark, then the survivors get the cavity mark back.
            removed = false
            seen = (T.stamp += UInt32(1))
            queue = T.queue; empty!(queue)
            for k in 1:seeds
                @inbounds T.mark[cav[k]] = seen
                push!(queue, cav[k])
            end
            head = 1
            while head <= length(queue)
                t = queue[head]; head += 1
                for s in 1:4
                    next = _tn(T, t, s)
                    next == 0 && continue
                    if @inbounds(T.mark[next]) == stamp
                        @inbounds T.mark[next] = seen
                        push!(queue, next)
                    end
                end
            end
            for t in cav
                m = @inbounds T.mark[t]
                if m == seen
                    @inbounds T.mark[t] = stamp
                elseif m == stamp
                    @inbounds T.mark[t] = 0
                end
            end
        end
        # No vertex may lose all of its cells: each must lie on the boundary.
        vstamp = (T.vstamp += UInt32(1))
        for t in cav
            @inbounds(T.mark[t]) == stamp || continue
            for s in 1:4
                next = _tn(T, t, s)
                if next == 0 || @inbounds(T.mark[next]) != stamp
                    a, b, c = _face(T, t, s)
                    @inbounds begin
                        T.vmark[a] = vstamp; T.vmark[b] = vstamp; T.vmark[c] = vstamp
                    end
                end
            end
        end
        for k in seeds+1:length(cav)
            t = cav[k]
            @inbounds(T.mark[t]) == stamp || continue
            for s in 1:4
                if @inbounds(T.vmark[_tv(T, t, s)]) != vstamp
                    @inbounds T.mark[t] = 0
                    changed = true; removed = true
                    break
                end
            end
        end
        changed || break
    end
    return nothing
end

# Collect the oriented boundary faces of the marked cavity into `T.bnd`.
# Returns false when an edge of the boundary is not shared by exactly two faces.
function _cavity_boundary!(T::Triangulation, stamp::UInt32)::Bool
    bnd = T.bnd; empty!(bnd)
    for t in T.cav
        @inbounds(T.mark[t]) == stamp || continue
        for s in 1:4
            next = _tn(T, t, s)
            if next == 0 || @inbounds(T.mark[next]) != stamp
                a, b, c = _face(T, t, s)
                slot = next == 0 ? 0 : _neighbor_slot(T, next, t)
                push!(bnd, (a, b, c, next, Int32(slot)))
            end
        end
    end
    _edges_reset!(T, 3 * length(bnd))
    open = 0
    for (a, b, c, _, _) in bnd
        for key in (_ekey(b, c), _ekey(a, c), _ekey(a, b))
            position, found = _edges_find!(T, key)
            if found
                count = @inbounds T.evals[position]
                count >= 2 && return false
                @inbounds T.evals[position] = 2
                open -= 1
            else
                @inbounds T.evals[position] = 1
                open += 1
            end
        end
    end
    return open == 0
end

# ── recovery state ──────────────────────────────────────────────────────────

# Scratch state of the exhaustive cavity refill (see `_refill!`).
mutable struct Solver
    vertices::Vector{Int32}                     # admissible cell vertices
    front::Vector{NTuple{3,Int32}}              # oriented faces still needing a cell on their positive side
    closed::Vector{NTuple{3,Int32}}             # sorted faces with cells already settled on both sides
    cells::Vector{NTuple{4,Int32}}
    nodes::Int
    limit::Int
    volume::Float64
    budget::Float64
    segments::Vector{Tuple{Int32,Int32}}        # segments that must be edges: no cell may cross them
    forbidden::Vector{Tuple{Int32,Int32}}       # vertex pairs no cell may join
    membranes::Vector{NTuple{3,Int32}}          # triangles that must be faces: no cell may cross them
    own::Int32                                  # vertices above this index were added for the refill
    memo::Dict{NTuple{4,Int32},Bool}            # cells already tested for orientation and crossing
end
Solver() = Solver(Int32[], NTuple{3,Int32}[], NTuple{3,Int32}[], NTuple{4,Int32}[], 0, 0, 0.0, 0.0,
                  Tuple{Int32,Int32}[], Tuple{Int32,Int32}[], NTuple{3,Int32}[], typemax(Int32),
                  Dict{NTuple{4,Int32},Bool}())

mutable struct Recovery
    T::Triangulation
    faces::Set{NTuple{3,Int32}}                 # sorted surface facets
    edges::Set{UInt64}                          # surface edges
    star::Vector{Int32}
    ring_vertices::Vector{Int32}
    ring_cells::Vector{Int32}
    olds::Vector{Int32}
    news::Vector{NTuple{4,Int32}}
    facemap::Dict{NTuple{3,Int32},Int32}
    outer::Vector{NTuple{5,Int32}}
    outer_matched::Vector{Bool}
    links::Vector{NTuple{4,Int32}}
    outer_links::Vector{NTuple{4,Int32}}
    kept_faces::Vector{NTuple{3,Int32}}
    kept_edges::Vector{UInt64}
    new_edges::Vector{UInt64}
    ids::Vector{Int32}
    ok::Matrix{Bool}
    split::Matrix{Int8}
    triangle::Array{Int8,3}
    sides::Vector{Int8}
    flips::Int
    steiner::Vector{Int32}                      # vertices added by recovery, after the cube corners
    cavity::Vector{Int32}
    crossings::Vector{NTuple{7,Float64}}        # nearest obstruction met by `_edge_cavity!`: point, normal, size
    crossed::Vector{Tuple{Int32,Int32}}         # edges `_edge_cavity!` found crossed within rounding
    wings::Dict{UInt64,NTuple{2,Int32}}         # surface edge -> vertices opposite it in its two facets
    solver::Solver
end

function Recovery(T::Triangulation)
    return Recovery(T, Set{NTuple{3,Int32}}(), Set{UInt64}(), Int32[], Int32[], Int32[],
                    Int32[], NTuple{4,Int32}[], Dict{NTuple{3,Int32},Int32}(),
                    NTuple{5,Int32}[], Bool[], NTuple{4,Int32}[], NTuple{4,Int32}[],
                    NTuple{3,Int32}[], UInt64[], UInt64[], Int32[],
                    Matrix{Bool}(undef, _MAX_RING, _MAX_RING),
                    Matrix{Int8}(undef, _MAX_RING, _MAX_RING),
                    Array{Int8,3}(undef, _MAX_RING, _MAX_RING, _MAX_RING), Int8[], 0,
                    Int32[], Int32[], NTuple{7,Float64}[], Tuple{Int32,Int32}[],
                    Dict{UInt64,NTuple{2,Int32}}(), Solver())
end

function _star!(R::Recovery, v::Int32)
    T = R.T; out = R.star; empty!(out)
    t0 = @inbounds T.v2t[v]
    (t0 != 0 && T.alive[t0] && _slot(T, t0, v) != 0) || throw(ArgumentError(
        "boundary recovery: stale vertex incidence"))
    stamp = (T.stamp += UInt32(1))
    @inbounds T.mark[t0] = stamp
    push!(out, t0); index = 1
    while index <= length(out)
        t = out[index]; index += 1; base = 4 * (t - 1)
        @inbounds for s in 1:4
            T.tv[base+s] == v && continue
            next = T.tn[base+s]
            (next == 0 || T.mark[next] == stamp) && continue
            T.mark[next] = stamp
            push!(out, next)
        end
    end
    return out
end

function _edge_cell(R::Recovery, a::Int32, b::Int32)::Int32
    for t in _star!(R, a)
        _slot(R.T, t, b) != 0 && return t
    end
    return Int32(0)
end

function _face_cell(R::Recovery, a::Int32, b::Int32, c::Int32)::Int32
    for t in _star!(R, a)
        (_slot(R.T, t, b) != 0 && _slot(R.T, t, c) != 0) && return t
    end
    return Int32(0)
end

# Vertices around edge (u,v), ordered so each cell (u,v,r[i],r[i+1]) is positive.
function _ring!(R::Recovery, u::Int32, v::Int32, t0::Int32)::Bool
    T = R.T; rv = R.ring_vertices; rt = R.ring_cells
    empty!(rv); empty!(rt)
    w1 = Int32(0); w2 = Int32(0)
    for s in 1:4
        w = _tv(T, t0, s)
        (w == u || w == v) && continue
        w1 == 0 ? (w1 = w) : (w2 = w)
    end
    x, y = _vol(T, u, v, w1, w2) > 0 ? (w1, w2) : (w2, w1)
    t = t0
    while true
        push!(rt, t); push!(rv, x)
        next = _tn(T, t, _slot(T, t, x))
        next == 0 && return false
        w = Int32(0)
        for s in 1:4
            candidate = _tv(T, next, s)
            (candidate == u || candidate == v || candidate == y) && continue
            w = candidate
        end
        x = y; y = w; t = next
        t == t0 && break
        length(rt) > 4096 && throw(ArgumentError("boundary recovery: edge ring does not close"))
    end
    return true
end

@inline function _even_permutation(f::NTuple{3,Int32}, g::NTuple{3,Int32})
    return g == f || g == (f[2], f[3], f[1]) || g == (f[3], f[1], f[2])
end

# Replace cells `olds` by `news`. The replacement is applied only when every new
# cell is strictly positive, the oriented boundary is unchanged, interior faces
# pair up, no vertex is lost, and no surface facet or edge disappears.
function _replace!(R::Recovery, olds::Vector{Int32}, news::Vector{NTuple{4,Int32}};
                   drop::Int32=Int32(0))::Bool
    T = R.T
    for (a, b, c, d) in news
        _vol(T, a, b, c, d) > 0 || return false
    end
    stamp = (T.stamp += UInt32(1))
    for t in olds
        @inbounds T.mark[t] = stamp
    end
    facemap = R.facemap; empty!(facemap)
    outer = R.outer; empty!(outer)
    kept_faces = R.kept_faces; empty!(kept_faces)
    kept_edges = R.kept_edges; empty!(kept_edges)
    for t in olds
        for s in 1:4
            f = _face(T, t, s); key = _sort3(f); next = _tn(T, t, s)
            if next != 0 && @inbounds(T.mark[next]) == stamp
                key in R.faces && push!(kept_faces, key)
            else
                slot = next == 0 ? 0 : _neighbor_slot(T, next, t)
                push!(outer, (f[1], f[2], f[3], next, Int32(slot)))
                facemap[key] = -Int32(length(outer))
            end
        end
        for (i, j) in _EDGES
            key = _ekey(_tv(T, t, i), _tv(T, t, j))
            key in R.edges && push!(kept_edges, key)
        end
    end
    matched = R.outer_matched; resize!(matched, length(outer)); fill!(matched, false)
    links = R.links; empty!(links)
    outer_links = R.outer_links; empty!(outer_links)
    new_edges = R.new_edges; empty!(new_edges)
    vstamp = (T.vstamp += UInt32(1))
    pending = 0
    for (index, cell) in enumerate(news)
        for s in 1:4
            f = (cell[_FV[s][1]], cell[_FV[s][2]], cell[_FV[s][3]]); key = _sort3(f)
            code = get(facemap, key, Int32(0))
            if code < 0 && code != typemin(Int32)
                record = outer[-code]
                (matched[-code] || !_even_permutation((record[1], record[2], record[3]), f)) && return false
                matched[-code] = true
                push!(outer_links, (Int32(index), Int32(s), record[4], record[5]))
            elseif code == 0
                facemap[key] = Int32(4 * (index - 1) + s)
                pending += 1
            elseif code == typemin(Int32)
                return false
            else
                other = (code - 1) ÷ 4 + 1; other_slot = (code - 1) % 4 + 1
                g = (news[other][_FV[other_slot][1]], news[other][_FV[other_slot][2]],
                     news[other][_FV[other_slot][3]])
                _even_permutation(g, f) && return false
                push!(links, (Int32(index), Int32(s), Int32(other), Int32(other_slot)))
                facemap[key] = typemin(Int32)
                pending -= 1
            end
        end
        for (i, j) in _EDGES
            push!(new_edges, _ekey(cell[i], cell[j]))
        end
        @inbounds for i in 1:4
            T.vmark[cell[i]] = vstamp
        end
    end
    (pending == 0 && all(matched)) || return false
    for key in kept_faces
        get(facemap, key, Int32(0)) == typemin(Int32) || return false
    end
    if !isempty(kept_edges)
        sort!(new_edges)
        for key in kept_edges
            isempty(searchsorted(new_edges, key)) && return false
        end
    end
    for t in olds, s in 1:4
        vertex = _tv(T, t, s)
        (vertex == drop || @inbounds(T.vmark[vertex]) == vstamp) || return false
    end
    for t in olds
        _kill!(T, t)
    end
    ids = R.ids; empty!(ids)
    for (a, b, c, d) in news
        push!(ids, _newcell!(T, a, b, c, d))
    end
    for (index, s, next, slot) in outer_links
        cell = ids[index]
        @inbounds T.tn[4*(cell-1)+s] = next
        next != 0 && @inbounds(T.tn[4*(next-1)+slot] = cell)
    end
    for (index, s, other, other_slot) in links
        cell = ids[index]; mate = ids[other]
        @inbounds T.tn[4*(cell-1)+s] = mate
        @inbounds T.tn[4*(mate-1)+other_slot] = cell
    end
    R.flips += 1
    return true
end

# Whether chord (i,k) of the ring polygon properly crosses chord (lo,hi).
@inline function _crosses(i::Int, k::Int, lo::Int, hi::Int)
    lo == 0 && return false
    (i == lo || i == hi || k == lo || k == hi) && return false
    return (lo < i < hi) != (lo < k < hi)
end

# Remove edge (u,v) by retriangulating its ring. `sides`, when given, holds a
# sign per ring vertex and forbids new edges joining strictly opposite signs;
# `require` names two ring vertices that must become joined by an edge.
function _remove_edge!(R::Recovery, u::Int32, v::Int32, t0::Int32;
                       plane::Union{Nothing,NTuple{3,Int32}}=nothing,
                       require::Tuple{Int32,Int32}=(Int32(0), Int32(0)))::Bool
    T = R.T
    _ekey(u, v) in R.edges && return false
    _ring!(R, u, v, t0) || return false
    rv = R.ring_vertices
    n = length(rv)
    3 <= n <= _MAX_RING || return false
    sides = R.sides; resize!(sides, n)
    if plane === nothing
        fill!(sides, Int8(0))
    else
        for i in 1:n
            w = rv[i]
            sides[i] = (w == plane[1] || w == plane[2] || w == plane[3]) ? Int8(0) :
                       Int8(_vol(T, plane[1], plane[2], plane[3], w))
        end
    end
    lo = 0; hi = 0
    if require[1] != 0
        for i in 1:n
            rv[i] == require[1] && (lo = i)
            rv[i] == require[2] && (hi = i)
        end
        (lo == 0 || hi == 0) && (lo = 0; hi = 0)
        lo > hi && ((lo, hi) = (hi, lo))
    end
    ok = R.ok; split = R.split; triangle = R.triangle
    @inbounds for i in 1:n, k in 1:n, j in 1:n
        triangle[i, k, j] = 0
    end
    @inbounds for i in 1:n-1
        ok[i, i+1] = true
    end
    @inbounds for len in 2:n-1, i in 1:n-len
        j = i + len
        ok[i, j] = false
        for k in i+1:j-1
            (ok[i, k] && ok[k, j]) || continue
            if k > i + 1
                (sides[i] * sides[k] < 0 || _crosses(i, k, lo, hi)) && continue
            end
            if j > k + 1
                (sides[k] * sides[j] < 0 || _crosses(k, j, lo, hi)) && continue
            end
            state = triangle[i, k, j]
            if state == 0
                good = _vol(T, rv[i], rv[k], rv[j], v) > 0 && _vol(T, rv[k], rv[i], rv[j], u) > 0
                state = good ? Int8(1) : Int8(2)
                triangle[i, k, j] = state
            end
            state == 1 || continue
            ok[i, j] = true; split[i, j] = k
            break
        end
    end
    @inbounds ok[1, n] || return false
    news = R.news; empty!(news)
    stack = Tuple{Int,Int}[(1, n)]
    while !isempty(stack)
        i, j = pop!(stack)
        j - i >= 2 || continue
        k = Int(@inbounds split[i, j])
        push!(news, (rv[i], rv[k], rv[j], v))
        push!(news, (rv[k], rv[i], rv[j], u))
        push!(stack, (i, k)); push!(stack, (k, j))
    end
    olds = R.olds; empty!(olds); append!(olds, R.ring_cells)
    return _replace!(R, olds, news)
end

# 2-3 flip of the face shared by cells `t` and `next`, creating the edge that
# joins their apices.
function _flip23!(R::Recovery, t::Int32, s::Int, next::Int32)::Bool
    T = R.T
    x, y, z = _face(T, t, s)
    a = _tv(T, t, s)
    apex = Int32(0)
    for k in 1:4
        w = _tv(T, next, k)
        (w == x || w == y || w == z) && continue
        apex = w
    end
    news = R.news; empty!(news)
    push!(news, (apex, a, x, y)); push!(news, (apex, a, y, z)); push!(news, (apex, a, z, x))
    olds = R.olds; empty!(olds); push!(olds, t); push!(olds, next)
    return _replace!(R, olds, news)
end

# Where the segment from vertex `a` towards vertex `b` leaves the star of `a`:
# kind 0, the edge is present; 1, through the interior of the face opposite
# `a` in `cell` (slot given); 2, through the edge (e1,e2) of `cell`; 3,
# through vertex e1; -1, undetermined.
function _exit(R::Recovery, a::Int32, b::Int32)
    T = R.T
    star = _star!(R, a)
    for t in star
        _slot(T, t, b) != 0 && return (0, t, 0, Int32(0), Int32(0))
    end
    for t in star
        s = _slot(T, t, a)
        x, y, z = _face(T, t, s)
        s1 = -_vol(T, a, x, y, b); s1 < 0 && continue
        s2 = -_vol(T, a, y, z, b); s2 < 0 && continue
        s3 = -_vol(T, a, z, x, b); s3 < 0 && continue
        zeros = (s1 == 0) + (s2 == 0) + (s3 == 0)
        zeros == 0 && return (1, t, s, Int32(0), Int32(0))
        if zeros == 1
            e1, e2 = s1 == 0 ? (x, y) : (s2 == 0 ? (y, z) : (z, x))
            return (2, t, 0, e1, e2)
        end
        w = (s1 == 0 && s3 == 0) ? x : ((s1 == 0 && s2 == 0) ? y : z)
        return (3, t, 0, w, Int32(0))
    end
    return (-1, Int32(0), 0, Int32(0), Int32(0))
end

@inline function _apex(T::Triangulation, t::Integer, x::Int32, y::Int32, z::Int32)::Int32
    for k in 1:4
        w = _tv(T, t, k)
        (w == x || w == y || w == z) || return w
    end
    return Int32(0)
end

# Recover edge (a,b) with flips. Returns :ok, :fail, or a diagnostic string
# when the configuration proves the input surface is not a valid boundary.
function _recover_edge!(R::Recovery, a::Int32, b::Int32)
    T = R.T
    for _ in 1:40
        kind, cell, slot, e1, e2 = _exit(R, a, b)
        kind == 0 && return :ok
        if kind == 1
            x, y, z = _face(T, cell, slot)
            _sort3(x, y, z) in R.faces &&
                return "surface edge ($a,$b) crosses surface facet ($x,$y,$z)"
            next = _tn(T, cell, slot)
            next == 0 && return :fail
            _flip23!(R, cell, slot, next) && continue
            apex = _apex(T, next, x, y, z)
            progressed = false
            for (p, q) in ((x, y), (y, z), (z, x))
                _vol(T, apex, a, p, q) <= 0 || continue
                if _remove_edge!(R, p, q, cell)
                    progressed = true
                    break
                end
            end
            progressed || return :fail
        elseif kind == 2
            _ekey(e1, e2) in R.edges &&
                return "surface edge ($a,$b) crosses surface edge ($e1,$e2)"
            _remove_edge!(R, e1, e2, cell; require=(a, b)) || _remove_edge!(R, e1, e2, cell) ||
                return :fail
        elseif kind == 3
            return "surface edge ($a,$b) passes through vertex $e1"
        else
            return :fail
        end
    end
    return :fail
end

# Recover triangle (a,b,c), whose three edges are present, with edge removals
# that never create an edge crossing the triangle's plane.
function _recover_face!(R::Recovery, a::Int32, b::Int32, c::Int32)
    T = R.T
    for _ in 1:40
        _face_cell(R, a, b, c) != 0 && return :ok
        progressed = false
        for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b))
            t0 = _edge_cell(R, p, q)
            t0 == 0 && return :fail
            _ring!(R, p, q, t0) || continue
            rv = R.ring_vertices; n = length(rv)
            u = Int32(0); v = Int32(0); wedge = Int32(0)
            for i in 1:n
                x = rv[i]; y = rv[i == n ? 1 : i + 1]
                (x == r || y == r) && continue
                if _vol(T, p, q, x, r) > 0 && _vol(T, p, q, r, y) > 0
                    u = x; v = y; wedge = R.ring_cells[i]
                    break
                end
            end
            u == 0 && continue
            _ekey(u, v) in R.edges &&
                return "surface edge ($u,$v) pierces surface facet ($a,$b,$c)"
            if _remove_edge!(R, u, v, wedge; plane=(a, b, c))
                progressed = true
                break
            end
        end
        progressed || return :fail
    end
    return :fail
end

# ── cavity refill ───────────────────────────────────────────────────────────
# A constraint that flips cannot recover is recovered by removing the cells it
# crosses and refilling that cavity by exhaustive search: cells are added one
# at a time on the positive side of a face that still lacks one, until no such
# face is left. Every cell is strictly positive and the faces pair up, so the
# result has the cavity's boundary and therefore tiles it exactly; no
# intersection test is involved. The search is bounded, and a cavity it cannot
# refill is reported as a failure.

# Approximate six-fold signed volume, positive for a valid cell.
@inline function _volf(T::Triangulation, a::Integer, b::Integer, c::Integer, d::Integer)
    @inbounds begin
        ax = T.x[a]; ay = T.y[a]; az = T.z[a]
        ux = T.x[b] - ax; uy = T.y[b] - ay; uz = T.z[b] - az
        vx = T.x[c] - ax; vy = T.y[c] - ay; vz = T.z[c] - az
        wx = T.x[d] - ax; wy = T.y[d] - ay; wz = T.z[d] - az
    end
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)
end

@inline function _dist2(T::Triangulation, a::Integer, b::Integer)
    @inbounds return (T.x[a] - T.x[b])^2 + (T.y[a] - T.y[b])^2 + (T.z[a] - T.z[b])^2
end

# Oriented face with its smallest vertex first.
@inline function _canonical(a::Int32, b::Int32, c::Int32)
    (a < b && a < c) && return (a, b, c)
    (b < a && b < c) && return (b, c, a)
    return (c, a, b)
end
@inline _reversed(f::NTuple{3,Int32}) = (f[1], f[3], f[2])

@inline _cell(T::Triangulation, t::Integer) = (_tv(T, t, 1), _tv(T, t, 2), _tv(T, t, 3), _tv(T, t, 4))

# Whether `target` lies in the closed cone that the positively oriented `cell`
# spans at its vertex in slot `i`.
@inline function _in_cone(T::Triangulation, cell::NTuple{4,Int32}, i::Int, target::Int32)::Bool
    apex = cell[i]; f = _FV[i]
    x = cell[f[1]]; y = cell[f[2]]; z = cell[f[3]]
    return _vol(T, apex, x, y, target) <= 0 && _vol(T, apex, y, z, target) <= 0 &&
           _vol(T, apex, z, x, target) <= 0
end

@inline _consistent(x::Int, y::Int, z::Int) = (x >= 0 && y >= 0 && z >= 0) || (x <= 0 && y <= 0 && z <= 0)

# Whether the open segment (a,b) meets the closed, positively oriented cell
# (p1,p2,p3,p4), none of whose vertices is `a` or `b`.
function _segment_meets_cell(T::Triangulation, a::Int32, b::Int32, p1::Int32, p2::Int32,
                             p3::Int32, p4::Int32)::Bool
    # A face plane with both endpoints strictly outside separates them.
    (_vol(T, p2, p4, p3, a) < 0 && _vol(T, p2, p4, p3, b) < 0) && return false
    (_vol(T, p1, p3, p4, a) < 0 && _vol(T, p1, p3, p4, b) < 0) && return false
    (_vol(T, p1, p4, p2, a) < 0 && _vol(T, p1, p4, p2, b) < 0) && return false
    (_vol(T, p1, p2, p3, a) < 0 && _vol(T, p1, p2, p3, b) < 0) && return false
    # Otherwise they meet exactly when the line does: seen along the line, its
    # trace must lie in the closed projection of one of the four faces.
    o12 = _vol(T, a, b, p1, p2); o13 = _vol(T, a, b, p1, p3); o14 = _vol(T, a, b, p1, p4)
    o23 = _vol(T, a, b, p2, p3); o24 = _vol(T, a, b, p2, p4); o34 = _vol(T, a, b, p3, p4)
    return _consistent(o12, o23, -o13) || _consistent(o12, o24, -o14) ||
           _consistent(o13, o34, -o14) || _consistent(o23, o34, -o24)
end

# Whether a positively oriented cell rules out (a,b) as an edge.
function _cell_blocks_edge(T::Triangulation, a::Int32, b::Int32, cell::NTuple{4,Int32})::Bool
    ia = 0; ib = 0
    for i in 1:4
        cell[i] == a && (ia = i)
        cell[i] == b && (ib = i)
    end
    (ia != 0 && ib != 0) && return false
    ia != 0 && return _in_cone(T, cell, ia, b)
    ib != 0 && return _in_cone(T, cell, ib, a)
    return _segment_meets_cell(T, a, b, cell[1], cell[2], cell[3], cell[4])
end

# Whether some edge of cell `t` strictly pierces the open triangle (a,b,c).
_cell_crosses_triangle(T::Triangulation, t::Integer, a::Int32, b::Int32, c::Int32) =
    _cell_crosses_triangle(T, _cell(T, t), a, b, c)

function _cell_crosses_triangle(T::Triangulation, cell::NTuple{4,Int32}, a::Int32, b::Int32,
                                c::Int32)::Bool
    sides = ntuple(i -> (v = cell[i]; (v == a || v == b || v == c) ? 0 : _vol(T, a, b, c, v)), 4)
    for (i, j) in _EDGES
        sides[i] * sides[j] < 0 || continue
        p = cell[i]; q = cell[j]
        o1 = _vol(T, p, q, a, b); o2 = _vol(T, p, q, b, c); o3 = _vol(T, p, q, c, a)
        ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0)) && return true
    end
    return false
end

# Whether a cell with apex `w` may be placed on the positive side of front face `f`.
function _admissible(R::Recovery, S::Solver, f::NTuple{3,Int32}, w::Int32)::Bool
    T = R.T; x, y, z = f
    (w == x || w == y || w == z) && return false
    # Orientation and crossing depend only on the cell, not on the search state.
    cell = (x, y, z, w)
    sound = get(S.memo, cell, nothing)
    if sound === nothing
        sound = _vol(T, x, y, z, w) > 0 && !_flat(T, x, y, z, w)
        if sound
            for (a, b) in S.forbidden
                if (a == x || a == y || a == z || a == w) && (b == x || b == y || b == z || b == w)
                    sound = false
                    break
                end
            end
        end
        if sound
            for (a, b) in S.segments
                if _cell_blocks_edge(T, a, b, cell)
                    sound = false
                    break
                end
            end
        end
        if sound
            for (a, b, c) in S.membranes
                if _cell_crosses_triangle(T, cell, a, b, c)
                    sound = false
                    break
                end
            end
        end
        S.memo[cell] = sound
    end
    sound || return false
    for g in ((x, y, w), (y, z, w), (z, x, w))
        face = _canonical(g[1], g[2], g[3])
        _reversed(face) in S.front && continue       # bonds to a cell already waiting there
        face in S.front && return false              # the waiting cell is on the other side
        _sort3(face) in S.closed && return false     # both sides are already settled
    end
    return true
end

@inline function _shape(T::Triangulation, x::Int32, y::Int32, z::Int32, w::Int32)
    total = _dist2(T, x, y) + _dist2(T, y, z) + _dist2(T, z, x) + _dist2(T, x, w) +
            _dist2(T, y, w) + _dist2(T, z, w)
    return total > 0 ? _volf(T, x, y, z, w) / (total * sqrt(total)) : 0.0
end

function _search!(R::Recovery, S::Solver, required::Vector{Int32})::Bool
    T = R.T
    if isempty(S.front)
        for v in required
            any(cell -> v in cell, S.cells) || return false
        end
        return true
    end
    (S.nodes += 1) > S.limit && return false
    # Branch on the face with the fewest admissible apices.
    best = 0; fewest = typemax(Int)
    for (k, f) in enumerate(S.front)
        count = 0
        for w in S.vertices
            if _admissible(R, S, f, w)
                count += 1
                count >= fewest && break
            end
        end
        count == 0 && return false
        if count < fewest
            fewest = count; best = k
            fewest == 1 && break
        end
    end
    f = S.front[best]
    x, y, z = f
    candidates = Tuple{Float64,Int32}[]
    for w in S.vertices
        _admissible(R, S, f, w) || continue
        # Added vertices come first, in the order offered: they were only
        # offered because the cavity's own vertices could not fill it. Those
        # follow, best shaped cell first.
        push!(candidates, (w > S.own ? -Inf : -_shape(T, x, y, z, w), w))
    end
    sort!(candidates)
    saved = copy(S.front); settled = length(S.closed)
    for (_, w) in candidates
        volume = _volf(T, x, y, z, w)
        S.volume + volume <= S.budget * (1 + 1.0e-9) + 1.0e-300 || continue
        deleteat!(S.front, best); push!(S.closed, _sort3(f))
        for g in ((x, y, w), (y, z, w), (z, x, w))
            face = _canonical(g[1], g[2], g[3])
            mate = findfirst(==(_reversed(face)), S.front)
            if mate === nothing
                push!(S.front, face)
            else
                deleteat!(S.front, mate); push!(S.closed, _sort3(face))
            end
        end
        push!(S.cells, (x, y, z, w)); S.volume += volume
        _search!(R, S, required) && return true
        pop!(S.cells); S.volume -= volume
        resize!(S.front, length(saved)); copyto!(S.front, saved)
        resize!(S.closed, settled)
        S.nodes > S.limit && return false
    end
    return false
end

# Replace `cavity` by cells found by search. `edge` names a segment that must
# become an edge; `membrane` an oriented triangle that must become a face;
# `extra` vertices that may be used in addition to the cavity's own.
function _refill!(R::Recovery, cavity::Vector{Int32};
                  edge::Tuple{Int32,Int32}=(Int32(0), Int32(0)),
                  required::Vector{NTuple{3,Int32}}=NTuple{3,Int32}[],
                  forbidden::Vector{Tuple{Int32,Int32}}=Tuple{Int32,Int32}[],
                  extra::Vector{Int32}=Int32[], limit::Int=3000)::Bool
    T = R.T; S = R.solver
    empty!(S.vertices); empty!(S.front); empty!(S.closed); empty!(S.cells)
    stamp = (T.stamp += UInt32(1))
    for t in cavity
        @inbounds T.mark[t] = stamp
    end
    vstamp = (T.vstamp += UInt32(1)); budget = 0.0
    membranes = NTuple{3,Int32}[]
    for t in cavity
        budget += _volf(T, _tv(T, t, 1), _tv(T, t, 2), _tv(T, t, 3), _tv(T, t, 4))
        for s in 1:4
            v = _tv(T, t, s)
            if @inbounds(T.vmark[v]) != vstamp
                @inbounds T.vmark[v] = vstamp
                push!(S.vertices, v)
            end
            next = _tn(T, t, s)
            f = _face(T, t, s)
            if next != 0 && @inbounds(T.mark[next]) == stamp
                # A surface facet inside the cavity stays, as a two-sided face.
                (next > t && _sort3(f) in R.faces) && push!(membranes, _canonical(f[1], f[2], f[3]))
            else
                push!(S.front, _canonical(f[1], f[2], f[3]))
            end
        end
    end
    # A vertex off the cavity boundary must be kept by the refill, and so must
    # a surface edge that runs through the cavity's interior.
    interior = Int32[]
    for v in S.vertices
        any(f -> v in f, S.front) || push!(interior, v)
    end
    empty!(S.forbidden); append!(S.forbidden, forbidden)
    empty!(S.segments)
    edge[1] != 0 && push!(S.segments, edge)
    for t in cavity, (i, j) in _EDGES
        p = _tv(T, t, i); q = _tv(T, t, j)
        _ekey(p, q) in R.edges || continue
        any(f -> p in f && q in f, S.front) && continue
        segment = p < q ? (p, q) : (q, p)
        segment in S.segments || push!(S.segments, segment)
    end
    S.own = isempty(extra) ? typemax(Int32) : minimum(extra) - Int32(1)
    append!(S.vertices, extra)
    sort!(S.vertices)
    for face in required
        push!(membranes, _canonical(face[1], face[2], face[3]))
    end
    empty!(S.membranes)
    for face in membranes
        push!(S.front, face); push!(S.front, _reversed(face)); push!(S.membranes, face)
    end
    S.nodes = 0; S.limit = limit; S.volume = 0.0; S.budget = budget
    empty!(S.memo)
    _search!(R, S, interior) || return false
    news = R.news; empty!(news); append!(news, S.cells)
    return _replace!(R, cavity, news)
end

function _add_vertex!(R::Recovery, px::Float64, py::Float64, pz::Float64)::Int32
    T = R.T
    push!(T.x, px); push!(T.y, py); push!(T.z, pz)
    push!(T.v2t, Int32(0)); push!(T.vmark, UInt32(0))
    push!(R.steiner, Int32(length(T.x)))
    return Int32(length(T.x))
end

# Collect in `R.cavity` the cells the open segment (a,b) meets, widened by the
# rings of edges it crosses within rounding of an exact intersection, so that
# a nearly planar quadrilateral can simply exchange its diagonals.
function _edge_cavity!(R::Recovery, a::Int32, b::Int32)::Bool
    T = R.T
    cavity = R.cavity; empty!(cavity); empty!(R.crossings); empty!(R.crossed)
    star = copy(_star!(R, a))
    stamp = (T.stamp += UInt32(1))
    for t in star
        if _in_cone(T, _cell(T, t), _slot(T, t, a), b)
            @inbounds T.mark[t] = stamp
            push!(cavity, t)
        end
    end
    isempty(cavity) && return false
    index = 1
    while index <= length(cavity)
        t = cavity[index]; index += 1
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || @inbounds(T.mark[next]) == stamp) && continue
            cell = _cell(T, next)
            if _cell_blocks_edge(T, a, b, cell)
                @inbounds T.mark[next] = stamp
                push!(cavity, next)
            end
        end
        length(cavity) > 400 && return false
    end
    ax = T.x[a]; ay = T.y[a]; az = T.z[a]
    dx = T.x[b] - ax; dy = T.y[b] - ay; dz = T.z[b] - az
    dd = dx * dx + dy * dy + dz * dz
    exact = length(cavity)
    nearest = Inf
    for k in 1:exact
        t = cavity[k]
        for (i, j) in _EDGES
            p = _tv(T, t, i); q = _tv(T, t, j)
            (p == a || p == b || q == a || q == b) && continue
            # Closest approach of the two segments.
            ex = T.x[q] - T.x[p]; ey = T.y[q] - T.y[p]; ez = T.z[q] - T.z[p]
            rx = T.x[p] - ax; ry = T.y[p] - ay; rz = T.z[p] - az
            ee = ex * ex + ey * ey + ez * ez
            de = dx * ex + dy * ey + dz * ez
            dr = dx * rx + dy * ry + dz * rz
            er = ex * rx + ey * ry + ez * rz
            denominator = dd * ee - de * de
            denominator > 1.0e-12 * dd * ee || continue
            along = (dr * ee - er * de) / denominator
            across = (dr * de - er * dd) / denominator
            (0.02 < along < 0.98 && 0.02 < across < 0.98) || continue
            gx = rx + across * ex - along * dx; gy = ry + across * ey - along * dy
            gz = rz + across * ez - along * dz
            gap = gx * gx + gy * gy + gz * gz
            if gap < nearest
                # The edge passing closest to the segment locates the vertices
                # a refill may add: beside the segment, off the surface.
                nx = dy * ez - dz * ey; ny = dz * ex - dx * ez; nz = dx * ey - dy * ex
                len = sqrt(nx * nx + ny * ny + nz * nz)
                if len > 0
                    nearest = gap
                    empty!(R.crossings)
                    push!(R.crossings, (ax + along * dx, ay + along * dy, az + along * dz,
                                        nx / len, ny / len, nz / len, sqrt(min(dd, ee))))
                end
            end
            gap <= 1.0e-18 * max(dd, ee) || continue
            pair = p < q ? (p, q) : (q, p)
            pair in R.crossed || push!(R.crossed, pair)
            _ring!(R, p, q, t) || continue
            for other in R.ring_cells
                if @inbounds(T.mark[other]) != stamp
                    @inbounds T.mark[other] = stamp
                    push!(cavity, other)
                end
            end
        end
    end
    return true
end

# Collect in `R.cavity` the cells that the open triangle (a,b,c) crosses.
function _face_cavity!(R::Recovery, a::Int32, b::Int32, c::Int32)::Bool
    T = R.T
    seeds = Int32[]
    for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b))
        t0 = _edge_cell(R, p, q)
        t0 == 0 && return false
        _ring!(R, p, q, t0) || return false
        rv = R.ring_vertices; n = length(rv)
        for i in 1:n
            x = rv[i]; y = rv[i == n ? 1 : i + 1]
            (x == r || y == r) && continue
            if _vol(T, p, q, x, r) > 0 && _vol(T, p, q, r, y) > 0
                push!(seeds, R.ring_cells[i])
                break
            end
        end
    end
    isempty(seeds) && return false
    cavity = R.cavity; empty!(cavity)
    stamp = (T.stamp += UInt32(1))
    for t in seeds
        if @inbounds(T.mark[t]) != stamp
            @inbounds T.mark[t] = stamp
            push!(cavity, t)
        end
    end
    index = 1
    while index <= length(cavity)
        t = cavity[index]; index += 1
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || @inbounds(T.mark[next]) == stamp) && continue
            if _cell_crosses_triangle(T, next, a, b, c)
                @inbounds T.mark[next] = stamp
                push!(cavity, next)
            end
        end
        length(cavity) > 400 && return false
    end
    return true
end

# Insert the off-surface vertex `p` by a cavity that respects the surface: it
# does not cross a surface facet that is already a face and does not swallow
# a surface edge that is already an edge. Returns a new cell at `p`, or 0 when
# `p` cannot be inserted that way.
function _insert_guarded!(R::Recovery, p::Int32, start::Int32)::Int32
    T = R.T
    t0 = _locate(T, p, start)
    stamp = (T.stamp += UInt32(1))
    seeds = _seed_cavity!(T, p, t0, stamp)
    cav = T.cav
    if seeds != 1
        # On a face or edge of a cell, or on a vertex: not a place for it.
        for t in cav
            @inbounds T.mark[t] = 0
        end
        return Int32(0)
    end
    index = 1
    while index <= length(cav)
        t = cav[index]; index += 1
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || @inbounds(T.mark[next]) == stamp) && continue
            _insphere_approx(T, next, p) || continue
            _sort3(_face(T, t, s)) in R.faces && continue
            @inbounds T.mark[next] = stamp
            push!(cav, next)
        end
    end
    kept = Set{UInt64}()
    removed = false
    while true
        _repair_cavity!(T, p, seeds, stamp; removed)
        empty!(kept)
        for t in cav
            @inbounds(T.mark[t]) == stamp || continue
            for s in 1:4
                next = _tn(T, t, s)
                if next == 0 || @inbounds(T.mark[next]) != stamp
                    a, b, c = _face(T, t, s)
                    push!(kept, _ekey(a, b)); push!(kept, _ekey(b, c)); push!(kept, _ekey(a, c))
                end
            end
        end
        removed = false
        for (k, t) in enumerate(cav)
            @inbounds(T.mark[t]) == stamp || continue
            for (i, j) in _EDGES
                key = _ekey(_tv(T, t, i), _tv(T, t, j))
                (key in R.edges && !(key in kept)) || continue
                if k <= seeds
                    for other in cav
                        @inbounds T.mark[other] = 0
                    end
                    return Int32(0)
                end
                @inbounds T.mark[t] = 0
                removed = true
                break
            end
        end
        removed || break
    end
    return _fill_cavity!(T, p, seeds, stamp)
end

# Insert a vertex on either side of the surface next to the edge (a,b) that
# could not be recovered, so that the cells around it stop spanning the whole
# region. The side that ends up outside is carved away with the exterior.
function _seed_edge!(R::Recovery, a::Int32, b::Int32)::Bool
    T = R.T
    _edge_cavity!(R, a, b) || return false
    isempty(R.crossings) && return false
    cx, cy, cz, nx, ny, nz, extent = R.crossings[1]
    inserted = false
    for side in (1.0, -1.0), factor in (0.3, 0.1, 0.03)
        offset = side * factor * extent
        vertex = _add_vertex!(R, cx + offset * nx, cy + offset * ny, cz + offset * nz)
        if _insert_guarded!(R, vertex, T.v2t[a]) != 0
            inserted = true
            break
        end
    end
    return inserted
end

# Add to `cavity` the cells that the open triangle (a,b,c) crosses and that
# connect to it, so that the triangle can be required as a face of the refill.
function _extend_cavity!(R::Recovery, cavity::Vector{Int32}, a::Int32, b::Int32, c::Int32)::Bool
    T = R.T
    stamp = (T.stamp += UInt32(1))
    for t in cavity
        @inbounds T.mark[t] = stamp
    end
    index = 1
    while index <= length(cavity)
        t = cavity[index]; index += 1
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || @inbounds(T.mark[next]) == stamp) && continue
            if _cell_crosses_triangle(T, next, a, b, c)
                @inbounds T.mark[next] = stamp
                push!(cavity, next)
            end
        end
        length(cavity) > 400 && return false
    end
    return true
end

# Add to `cavity` the neighbours of its thin cells. A thin cell whose faces
# bound the cavity is all but the only cell that fits between them, and it
# leaves no room for an added vertex to see both, so a refill has no freedom
# there until those neighbours are released too. Returns whether the cavity grew.
function _widen_cavity!(R::Recovery, cavity::Vector{Int32})::Bool
    T = R.T
    stamp = (T.stamp += UInt32(1))
    for t in cavity
        @inbounds T.mark[t] = stamp
    end
    before = length(cavity)
    for k in 1:before
        t = cavity[k]
        a, b, c, d = _cell(T, t)
        _flat(T, a, b, c, d, 1.0e-3) || continue
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || @inbounds(T.mark[next]) == stamp) && continue
            @inbounds T.mark[next] = stamp
            push!(cavity, next)
        end
    end
    return length(cavity) > before
end

# Try the refill with the cavity's own vertices, then with vertices added on
# either side of the point `(cx,cy,cz)` along the unit direction `(nx,ny,nz)`,
# and once more on the cavity widened around its flat cells.
function _refill_staged!(R::Recovery, cavity::Vector{Int32}, edge::Tuple{Int32,Int32},
                         required::Vector{NTuple{3,Int32}}, forbidden::Vector{Tuple{Int32,Int32}},
                         cx, cy, cz, nx, ny, nz, extent; widen::Bool=true)::Bool
    if _refill_offsets!(R, cavity, edge, required, forbidden, cx, cy, cz, nx, ny, nz, extent)
        return true
    end
    widen || return false
    wider = copy(cavity)
    for _ in 1:2
        _widen_cavity!(R, wider) || return false
        length(wider) <= 60 || return false
        if _refill_offsets!(R, wider, edge, required, forbidden, cx, cy, cz, nx, ny, nz, extent)
            return true
        end
    end
    return false
end

# The offset direction `(nx,ny,nz)` with its components across the cavity's
# thin cells removed. A vertex added beside the point `(cx,cy,cz)` sees both
# outer faces of a thin cell only from within the slab between them, so it has
# to move along that slab. Returns `nothing` when nothing changes or nothing
# is left of the direction.
function _slab_direction(R::Recovery, cavity::Vector{Int32}, nx, ny, nz)
    T = R.T
    dx = nx; dy = ny; dz = nz
    changed = false
    for _ in 1:3, t in cavity
        a, b, c, d = _cell(T, t)
        _flat(T, a, b, c, d, 1.0e-3) || continue
        # Normal of the cell's largest face stands in for the normal of its slab.
        mx = 0.0; my = 0.0; mz = 0.0; largest = 0.0
        for (p, q, r) in ((a, b, c), (a, b, d), (a, c, d), (b, c, d))
            ux = T.x[q] - T.x[p]; uy = T.y[q] - T.y[p]; uz = T.z[q] - T.z[p]
            vx = T.x[r] - T.x[p]; vy = T.y[r] - T.y[p]; vz = T.z[r] - T.z[p]
            wx = uy * vz - uz * vy; wy = uz * vx - ux * vz; wz = ux * vy - uy * vx
            area = wx * wx + wy * wy + wz * wz
            if area > largest
                largest = area; mx = wx; my = wy; mz = wz
            end
        end
        largest > 0 || continue
        across = (dx * mx + dy * my + dz * mz) / largest
        abs(across) * sqrt(largest) > 1.0e-6 || continue
        dx -= across * mx; dy -= across * my; dz -= across * mz
        changed = true
    end
    changed || return nothing
    len = sqrt(dx * dx + dy * dy + dz * dz)
    len > 1.0e-3 || return nothing
    return dx / len, dy / len, dz / len
end

function _refill_offsets!(R::Recovery, cavity::Vector{Int32}, edge::Tuple{Int32,Int32},
                          required::Vector{NTuple{3,Int32}}, forbidden::Vector{Tuple{Int32,Int32}},
                          cx, cy, cz, nx, ny, nz, extent)::Bool
    _refill!(R, cavity; edge, required, forbidden, limit=600) && return true
    slab = _slab_direction(R, cavity, nx, ny, nz)
    if slab !== nothing
        for factor in (0.3, 0.1, 0.03, 0.01, 0.001)
            offset = factor * extent
            extra = Int32[_add_vertex!(R, cx + offset * slab[1], cy + offset * slab[2], cz + offset * slab[3]),
                          _add_vertex!(R, cx - offset * slab[1], cy - offset * slab[2], cz - offset * slab[3])]
            _refill!(R, cavity; edge, required, forbidden, extra, limit=300) && return true
        end
    end
    for factor in (0.3, 0.1, 0.03, 0.01, 0.001)
        offset = factor * extent
        extra = Int32[_add_vertex!(R, cx + offset * nx, cy + offset * ny, cz + offset * nz),
                      _add_vertex!(R, cx - offset * nx, cy - offset * ny, cz - offset * nz)]
        _refill!(R, cavity; edge, required, forbidden, extra, limit=300) && return true
    end
    return false
end

# Recover edge (a,b) by refilling the cells it crosses. Its two facets are
# recovered with it when their other edges already exist: requiring them
# separates the two sides of the surface, which keeps the search small.
function _refill_edge!(R::Recovery, a::Int32, b::Int32)::Bool
    T = R.T
    _edge_cavity!(R, a, b) || return false
    isempty(R.crossings) && return false
    cx, cy, cz, nx, ny, nz, extent = R.crossings[1]
    forbidden = copy(R.crossed)
    pipe = copy(R.cavity)
    cavity = copy(pipe)
    required = NTuple{3,Int32}[]
    for c in get(R.wings, _ekey(a, b), (Int32(0), Int32(0)))
        c == 0 && continue
        (_edge_cell(R, a, c) != 0 && _edge_cell(R, b, c) != 0) || continue
        _face_cell(R, a, b, c) == 0 || continue
        any(t -> _slot(T, t, c) != 0, cavity) || continue
        _extend_cavity!(R, cavity, a, b, c) || continue
        push!(required, (a, b, c))
    end
    if !isempty(required)
        _refill_staged!(R, cavity, (a, b), required, forbidden, cx, cy, cz, nx, ny, nz, extent) &&
            return true
    end
    return _refill_staged!(R, pipe, (a, b), NTuple{3,Int32}[], forbidden, cx, cy, cz, nx, ny, nz,
                           extent)
end

# Recover facet (a,b,c), whose edges exist, by refilling the cells it crosses.
function _refill_face!(R::Recovery, a::Int32, b::Int32, c::Int32)::Bool
    _face_cavity!(R, a, b, c) || return false
    cavity = copy(R.cavity)
    T = R.T
    ux = T.x[b] - T.x[a]; uy = T.y[b] - T.y[a]; uz = T.z[b] - T.z[a]
    vx = T.x[c] - T.x[a]; vy = T.y[c] - T.y[a]; vz = T.z[c] - T.z[a]
    nx = uy * vz - uz * vy; ny = uz * vx - ux * vz; nz = ux * vy - uy * vx
    len = sqrt(nx * nx + ny * ny + nz * nz)
    len > 0 || return false
    extent = sqrt(min(_dist2(T, a, b), _dist2(T, b, c), _dist2(T, a, c)))
    return _refill_staged!(R, cavity, (Int32(0), Int32(0)), NTuple{3,Int32}[(a, b, c)],
                           Tuple{Int32,Int32}[],
                           (T.x[a] + T.x[b] + T.x[c]) / 3, (T.y[a] + T.y[b] + T.y[c]) / 3,
                           (T.z[a] + T.z[b] + T.z[c]) / 3, nx / len, ny / len, nz / len, extent)
end

# Remove added vertex `p` by merging it into a neighbour that sees its whole
# link. Returns false when no neighbour does.
function _collapse!(R::Recovery, p::Int32)::Bool
    T = R.T
    star = copy(_star!(R, p))
    candidates = Int32[]
    for t in star, s in 1:4
        w = _tv(T, t, s)
        if w == p
            continue
        elseif _tn(T, t, s) == 0
            return false
        end
        w in candidates || push!(candidates, w)
    end
    sort!(candidates; by=w -> (_dist2(T, w, p), w))
    news = NTuple{4,Int32}[]
    for w in candidates
        empty!(news)
        for t in star
            _slot(T, t, w) != 0 && continue
            push!(news, ntuple(i -> (vertex = _tv(T, t, i); vertex == p ? w : vertex), 4))
        end
        _replace!(R, star, news; drop=p) && return true
    end
    return false
end

# Cells enclosed by the surface: label the regions the surface facets
# separate, nest them from the outside in, and keep those at odd depth.
function _classify(R::Recovery, nfaces::Int)
    T = R.T
    ncells = length(T.alive)
    region = zeros(Int32, ncells)
    queue = Int32[]
    nregions = 0
    outside = Int32(0)
    for seed in Int32(1):Int32(ncells)
        (T.alive[seed] && region[seed] == 0) || continue
        nregions += 1
        region[seed] = nregions
        empty!(queue); push!(queue, seed); head = 1
        while head <= length(queue)
            t = queue[head]; head += 1
            for s in 1:4
                next = _tn(T, t, s)
                if next == 0
                    outside = Int32(nregions)
                    continue
                end
                region[next] == 0 || continue
                _sort3(_face(T, t, s)) in R.faces && continue
                region[next] = nregions
                push!(queue, next)
            end
        end
    end
    outside == 0 && return nothing, "no cell touches the enclosing cube"
    adjacency = [Int32[] for _ in 1:nregions]
    walls = 0
    for t in Int32(1):Int32(ncells)
        T.alive[t] || continue
        for s in 1:4
            next = _tn(T, t, s)
            (next == 0 || next < t) && continue
            _sort3(_face(T, t, s)) in R.faces || continue
            walls += 1
            first = region[t]; second = region[next]
            first == second &&
                return nothing, "the surface is not closed: a facet does not separate two regions"
            push!(adjacency[first], second); push!(adjacency[second], first)
        end
    end
    walls == nfaces || return nothing, "only $walls of $nfaces surface facets are cell faces"
    depth = fill(-1, nregions)
    depth[outside] = 0
    order = Int32[outside]; head = 1
    while head <= length(order)
        current = order[head]; head += 1
        for other in adjacency[current]
            if depth[other] < 0
                depth[other] = depth[current] + 1
                push!(order, other)
            elseif iseven(depth[other] - depth[current])
                return nothing, "the surface does not bound a consistent interior"
            end
        end
    end
    inside = falses(ncells)
    kept = 0
    for t in 1:ncells
        if T.alive[t] && isodd(depth[region[t]])
            inside[t] = true; kept += 1
        end
    end
    kept > 0 || return nothing, "the surface encloses no cells"
    return inside, ""
end

# ── driver ──────────────────────────────────────────────────────────────────

"""
    recover_boundary_flips(coords, tris; steiner=true) -> (mesh, reason, stats)

Tetrahedralize the region enclosed by the closed triangulated surface `tris`
(a `3 × F` matrix of rows into the `3 × N` matrix `coords`). Columns of
`coords` that no triangle references become interior vertices.

The vertices are tetrahedralized first; surface edges and facets missing from
that tetrahedralization are recovered with flips and, failing that, by
refilling the cells they cross. With `steiner=true` a refill may use vertices
added strictly inside the cells it replaces. Such vertices never lie on the
surface, so the boundary of the result is always exactly the input triangles;
those that end up inside the region and cannot be merged away again are
appended after the input columns.

On success `mesh` is a [`Mesh`](@ref) and `reason` is empty; otherwise `mesh`
is `nothing` and `reason` says which constraint could not be recovered or why
the input is not a valid boundary. `stats` reports the work done.
"""
function recover_boundary_flips(coords::AbstractMatrix{Float64}, tris::AbstractMatrix{<:Integer};
                                steiner::Bool=true)
    count = size(coords, 2); nfaces = size(tris, 2)
    missing_edges = 0; missing_faces = 0; refilled = 0; kept = 0; flips = 0
    report(cells=0) = (missing_edges=missing_edges, missing_faces=missing_faces, flips=flips,
                       refilled=refilled, steiner=kept, cells=cells)
    (size(coords, 1) == 3 && size(tris, 1) == 3) ||
        return nothing, "coordinates and triangles must have three rows", report()
    (count >= 4 && nfaces >= 4) ||
        return nothing, "a closed surface needs at least four vertices and facets", report()
    count <= typemax(Int32) ÷ 8 || return nothing, "too many vertices", report()
    all(isfinite, coords) || return nothing, "nonfinite coordinates", report()
    T = Triangulation(coords, count)
    R = Recovery(T)
    facets = Vector{NTuple{3,Int32}}(undef, nfaces)
    sizehint!(R.faces, nfaces); sizehint!(R.edges, 2 * nfaces)
    @inbounds for f in 1:nfaces
        a = Int32(tris[1, f]); b = Int32(tris[2, f]); c = Int32(tris[3, f])
        (1 <= a <= count && 1 <= b <= count && 1 <= c <= count && a != b && b != c && a != c) ||
            return nothing, "facet $f does not name three distinct vertices", report()
        key = _sort3(a, b, c)
        key in R.faces && return nothing, "duplicate facet ($a,$b,$c)", report()
        push!(R.faces, key); facets[f] = (a, b, c)
        for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b))
            edge = _ekey(p, q)
            push!(R.edges, edge)
            wing = get(R.wings, edge, (Int32(0), Int32(0)))
            R.wings[edge] = wing[1] == 0 ? (r, Int32(0)) : (wing[1], r)
        end
    end
    hint = Int32(1)
    for p in _insertion_order(T, count)
        result = _insert!(T, p, hint)
        result < 0 && return nothing, "vertices $(-result) and $p share coordinates", report()
        hint = result
    end
    # Edges first: a facet can only be recovered once its three edges exist.
    pending = Tuple{Int32,Int32}[]
    for key in R.edges
        a = Int32(key >> 32); b = Int32(key & 0xffffffff)
        _edge_cell(R, a, b) == 0 && push!(pending, (a, b))
    end
    sort!(pending)
    missing_edges = length(pending)
    # Stage 0 uses flips only, stage 1 adds refills, stage 2 first inserts
    # vertices beside each edge still missing.
    stage = 0
    while !isempty(pending)
        before = length(pending)
        remaining = Tuple{Int32,Int32}[]
        for (a, b) in pending
            stage == 2 && _edge_cell(R, a, b) == 0 && _seed_edge!(R, a, b)
            result = _recover_edge!(R, a, b)
            result === :fail && (result = _recover_edge!(R, b, a))
            result === :ok && continue
            flips = R.flips
            result isa String && return nothing, result, report()
            if stage > 0 && _refill_edge!(R, a, b)
                refilled += 1
                continue
            end
            push!(remaining, (a, b))
        end
        pending = remaining
        if length(pending) == before
            (stage == 2 || !steiner) && break
            stage += 1
        end
    end
    flips = R.flips
    isempty(pending) || return nothing, "$(length(pending)) of $missing_edges missing surface " *
        "edges could not be recovered (first: $(pending[1]))", report()
    waiting = NTuple{3,Int32}[]
    for (a, b, c) in facets
        _face_cell(R, a, b, c) == 0 && push!(waiting, (a, b, c))
    end
    missing_faces = length(waiting)
    refill = false
    while !isempty(waiting)
        before = length(waiting)
        remaining = NTuple{3,Int32}[]
        for (a, b, c) in waiting
            result = _recover_face!(R, a, b, c)
            result === :ok && continue
            flips = R.flips
            result isa String && return nothing, result, report()
            if refill && _refill_face!(R, a, b, c)
                refilled += 1
                continue
            end
            push!(remaining, (a, b, c))
        end
        waiting = remaining
        if length(waiting) == before
            (refill || !steiner) && break
            refill = true
        end
    end
    flips = R.flips
    isempty(waiting) || return nothing, "$(length(waiting)) of $missing_faces missing surface " *
        "facets could not be recovered (first: $(waiting[1]))", report()
    # Merge away the vertices a refill added, where a neighbour can take over.
    for p in R.steiner
        (T.v2t[p] != 0 && T.alive[T.v2t[p]] && _slot(T, T.v2t[p], p) != 0) || continue
        _collapse!(R, p) && (T.v2t[p] = 0)
    end
    flips = R.flips
    inside, reason = _classify(R, nfaces)
    inside === nothing && return nothing, reason, report()
    ncells = 0
    rows = Dict{Int32,Int32}()
    extra = Int32[]
    for t in Int32(1):Int32(length(T.alive))
        (T.alive[t] && inside[t]) || continue
        ncells += 1
        for s in 1:4
            vertex = _tv(T, t, s)
            vertex > count || continue
            vertex > count + 8 ||
                return nothing, "an interior cell touches the enclosing cube", report()
            if !haskey(rows, vertex)
                push!(extra, vertex); rows[vertex] = 0
            end
        end
    end
    sort!(extra)
    for (k, vertex) in enumerate(extra)
        rows[vertex] = Int32(count + k)
    end
    kept = length(extra)
    out = Matrix{Float64}(undef, 3, count + kept)
    @inbounds for i in 1:count
        out[1, i] = coords[1, i]; out[2, i] = coords[2, i]; out[3, i] = coords[3, i]
    end
    for (k, vertex) in enumerate(extra)
        out[1, count+k] = T.x[vertex]; out[2, count+k] = T.y[vertex]; out[3, count+k] = T.z[vertex]
    end
    tets = Matrix{Int32}(undef, 4, ncells)
    column = 0
    for t in Int32(1):Int32(length(T.alive))
        (T.alive[t] && inside[t]) || continue
        column += 1
        for s in 1:4
            vertex = _tv(T, t, s)
            tets[s, column] = vertex > count ? rows[vertex] : vertex
        end
    end
    return Mesh(out; tets=tets), "", report(ncells)
end

recover_boundary_flips(surface::Mesh; steiner::Bool=true) =
    recover_boundary_flips(surface.coords, surface.tris; steiner)

end # module BoundaryRecovery
