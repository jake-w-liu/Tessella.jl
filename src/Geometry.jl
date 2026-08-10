"""
    Geometry

Stage-5 native constructive primitives (PLAN.md §3 "Geometry", §5 native CSG
path). Each builder returns a **closed, manifold, outward-oriented** triangle
[`Mesh`](@ref) — a boundary surface ready for [`Mesh3D.tetrahedralize`](@ref) /
[`mesh_volume`](@ref). All are verified `Heal.is_meshable` and, when filled, to
reproduce the exact analytic volume.

These are the primitives the ASCENT `solid_model` emits (boxes, cylinders, and
holes/bores); a full Boolean CSG kernel that combines them into the literal
enclosure geometry is the remaining Stage-5 work.
"""
module Geometry

using ..MeshTypes: Mesh

export box_surface, cylinder_surface, box_tunnel_surface

# ── box ─────────────────────────────────────────────────────────────────────────
"""
    box_surface(x0,x1, y0,y1, z0,z1) -> Mesh

Closed surface of the axis-aligned box `[x0,x1]×[y0,y1]×[z0,z1]` (12 outward
triangles).
"""
function box_surface(x0::Real,x1::Real, y0::Real,y1::Real, z0::Real,z1::Real)
    (x0<x1 && y0<y1 && z0<z1) || throw(ArgumentError("box_surface: need x0<x1, y0<y1, z0<z1"))
    C = Float64[x0 x1 x1 x0 x0 x1 x1 x0; y0 y0 y1 y1 y0 y0 y1 y1; z0 z0 z0 z0 z1 z1 z1 z1]
    F = [(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),
         (2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
    t = Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    return Mesh(C; tris=t)
end

# ── cylinder (solid, watertight lateral + caps) ─────────────────────────────────
"""
    cylinder_surface(center, axis, radius, height; nθ=24, nz=2) -> Mesh

Closed surface of a solid cylinder: base `center`, unit-ish `axis`, `radius`,
`height`. `nθ` circumferential sectors, `nz` axial levels. Lateral wall + two cap
fans share the rim rings (watertight by construction).
"""
function cylinder_surface(center, axis, radius::Real, height::Real; nθ::Integer=24, nz::Integer=2)
    R=Float64(radius); H=Float64(height)
    (R>0 && H>0 && nθ>=3 && nz>=2) || throw(ArgumentError("cylinder_surface: R,H>0, nθ≥3, nz≥2"))
    c=(Float64(center[1]),Float64(center[2]),Float64(center[3]))
    ez=_unit((Float64(axis[1]),Float64(axis[2]),Float64(axis[3])))
    ax = abs(ez[1])<=abs(ez[2]) ? (abs(ez[1])<=abs(ez[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(ez[2])<=abs(ez[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    ex=_unit(_cross(ax,ez)); ey=_cross(ez,ex)
    on(θ,z) = (c[1]+R*cos(θ)*ex[1]+R*sin(θ)*ey[1]+z*ez[1],
               c[2]+R*cos(θ)*ex[2]+R*sin(θ)*ey[2]+z*ez[2],
               c[3]+R*cos(θ)*ex[3]+R*sin(θ)*ey[3]+z*ez[3])
    V=Tuple{Float64,Float64,Float64}[]
    for j in 0:nz-1, i in 0:nθ-1; push!(V, on(2π*i/nθ, H*j/(nz-1))); end
    ci=length(V)+1; push!(V, (c[1], c[2], c[3]))          # bottom centre
    cti=length(V)+1; push!(V, (c[1]+H*ez[1], c[2]+H*ez[2], c[3]+H*ez[3]))   # top centre
    idx(j,i)=(j-1)*nθ + mod(i,nθ) + 1
    Tr=NTuple{3,Int32}[]
    for j in 1:nz-1, i in 0:nθ-1                            # wall (outward)
        a=idx(j,i);b=idx(j,i+1);cc=idx(j+1,i+1);d=idx(j+1,i)
        push!(Tr,(Int32(a),Int32(b),Int32(cc))); push!(Tr,(Int32(a),Int32(cc),Int32(d)))
    end
    for i in 0:nθ-1; push!(Tr,(Int32(ci),Int32(idx(1,i+1)),Int32(idx(1,i)))); end   # bottom cap
    for i in 0:nθ-1; push!(Tr,(Int32(cti),Int32(idx(nz,i)),Int32(idx(nz,i+1)))); end # top cap
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    tm=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); tm[:,k]=Int32[f...]; end
    return Mesh(C; tris=tm)
end

# ── box with a rectangular through-tunnel (genus-1 bore) ────────────────────────
"""
    box_tunnel_surface(ox0,ox1, oy0,oy1, z0,z1, ix0,ix1, iy0,iy1) -> Mesh

Closed surface of the box `[ox0,ox1]×[oy0,oy1]×[z0,z1]` with a rectangular
through-tunnel `[ix0,ix1]×[iy0,iy1]` along `z` (a genus-1 solid — the coax-bore
analogue). The inner `[ix,iy]` rectangle must lie strictly inside the outer one.
"""
function box_tunnel_surface(ox0,ox1, oy0,oy1, z0,z1, ix0,ix1, iy0,iy1)
    (ox0<ix0<ix1<ox1 && oy0<iy0<iy1<oy1 && z0<z1) ||
        throw(ArgumentError("box_tunnel_surface: inner rectangle must be strictly inside the outer"))
    V=Tuple{Float64,Float64,Float64}[]
    for (x,y) in [(ox0,oy0),(ox1,oy0),(ox1,oy1),(ox0,oy1)]; push!(V,(Float64(x),Float64(y),Float64(z0))); end
    for (x,y) in [(ox0,oy0),(ox1,oy0),(ox1,oy1),(ox0,oy1)]; push!(V,(Float64(x),Float64(y),Float64(z1))); end
    for (x,y) in [(ix0,iy0),(ix1,iy0),(ix1,iy1),(ix0,iy1)]; push!(V,(Float64(x),Float64(y),Float64(z0))); end
    for (x,y) in [(ix0,iy0),(ix1,iy0),(ix1,iy1),(ix0,iy1)]; push!(V,(Float64(x),Float64(y),Float64(z1))); end
    Tr=NTuple{3,Int32}[]; q(a,b,c,d)=(push!(Tr,(Int32(a),Int32(b),Int32(c)));push!(Tr,(Int32(a),Int32(c),Int32(d))))
    q(1,9,10,2);q(2,10,11,3);q(3,11,12,4);q(4,12,9,1)      # bottom frame (-z)
    q(5,6,14,13);q(6,7,15,14);q(7,8,16,15);q(8,5,13,16)    # top frame (+z)
    q(1,2,6,5);q(2,3,7,6);q(3,4,8,7);q(4,1,5,8)            # outer walls
    q(9,13,14,10);q(10,14,15,11);q(11,15,16,12);q(12,16,13,9)  # inner tunnel walls
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    tm=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); tm[:,k]=Int32[f...]; end
    return Mesh(C; tris=tm)
end

@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline function _unit(a)
    l=sqrt(a[1]^2+a[2]^2+a[3]^2); l==0 && throw(ArgumentError("Geometry: zero-length axis"))
    (a[1]/l,a[2]/l,a[3]/l)
end

end # module Geometry
