"""
    Geometry

Stage-5 native constructive primitives (PLAN.md §3 "Geometry", §5 native CSG
path). Each builder returns a **closed, manifold, outward-oriented** triangle
[`Mesh`](@ref) — a boundary surface ready for [`Mesh3D.tetrahedralize`](@ref) /
[`mesh_volume`](@ref). All are verified `Heal.is_meshable` and, when filled, to
reproduce the exact analytic volume.

These are the primitives the ASCENT `solid_model` emits (boxes, cylinders, holes/
bores, and axis-aligned cavities via [`box_shell_surface`](@ref)); a full Boolean
CSG kernel that combines them into the literal enclosure geometry is the remaining
Stage-5 work.
"""
module Geometry

using ..MeshTypes: Mesh

export box_surface, cylinder_surface, box_tunnel_surface, box_shell_surface

@inline function _finite_real(x::Real, caller::AbstractString, name::AbstractString)
    y = try
        Float64(x)
    catch err
        throw(ArgumentError("$caller: $name must be Float64-representable: $(sprint(showerror, err))"))
    end
    isfinite(y) || throw(ArgumentError("$caller: $name must be finite (got $x)"))
    return y
end

@inline function _point3(p, caller::AbstractString, name::AbstractString)
    length(p) >= 3 || throw(ArgumentError("$caller: $name needs three coordinates"))
    return (_finite_real(p[1],caller,"$name[1]"),
            _finite_real(p[2],caller,"$name[2]"),
            _finite_real(p[3],caller,"$name[3]"))
end

# ── box ─────────────────────────────────────────────────────────────────────────
"""
    box_surface(x0,x1, y0,y1, z0,z1) -> Mesh

Closed surface of the axis-aligned box `[x0,x1]×[y0,y1]×[z0,z1]` (12 outward
triangles).
"""
function box_surface(x0::Real,x1::Real, y0::Real,y1::Real, z0::Real,z1::Real)
    x0=_finite_real(x0,"box_surface","x0"); x1=_finite_real(x1,"box_surface","x1")
    y0=_finite_real(y0,"box_surface","y0"); y1=_finite_real(y1,"box_surface","y1")
    z0=_finite_real(z0,"box_surface","z0"); z1=_finite_real(z1,"box_surface","z1")
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
    R=_finite_real(radius,"cylinder_surface","radius")
    H=_finite_real(height,"cylinder_surface","height")
    (R>0 && H>0 && nθ>=3 && nz>=2) || throw(ArgumentError("cylinder_surface: R,H>0, nθ≥3, nz≥2"))
    nθ <= typemax(Int) && nz <= typemax(Int) ||
        throw(ArgumentError("cylinder_surface: nθ and nz exceed the platform Int limit"))
    ntheta=Int(nθ); nlevels=Int(nz)
    nv = try
        Base.checked_add(Base.checked_mul(ntheta,nlevels), 2)
    catch
        throw(ArgumentError("cylinder_surface: requested node count overflows the platform Int limit"))
    end
    nv <= typemax(Int32) ||
        throw(ArgumentError("cylinder_surface: $nv nodes exceed the Int32 indexing limit"))
    try
        Base.checked_mul(Base.checked_mul(2,ntheta),nlevels)
    catch
        throw(ArgumentError("cylinder_surface: requested triangle count overflows the platform Int limit"))
    end
    c=_point3(center,"cylinder_surface","center")
    ez=_unit(_point3(axis,"cylinder_surface","axis"))
    ax = abs(ez[1])<=abs(ez[2]) ? (abs(ez[1])<=abs(ez[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
         (abs(ez[2])<=abs(ez[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    ex=_unit(_cross(ax,ez)); ey=_cross(ez,ex)
    on(θ,z) = (c[1]+R*cos(θ)*ex[1]+R*sin(θ)*ey[1]+z*ez[1],
               c[2]+R*cos(θ)*ex[2]+R*sin(θ)*ey[2]+z*ez[2],
               c[3]+R*cos(θ)*ex[3]+R*sin(θ)*ey[3]+z*ez[3])
    V=Tuple{Float64,Float64,Float64}[]
    sizehint!(V, nv)
    for j in 0:nlevels-1, i in 0:ntheta-1; push!(V, on(2π*i/ntheta, H*j/(nlevels-1))); end
    ci=length(V)+1; push!(V, (c[1], c[2], c[3]))          # bottom centre
    cti=length(V)+1; push!(V, (c[1]+H*ez[1], c[2]+H*ez[2], c[3]+H*ez[3]))   # top centre
    idx(j,i)=(j-1)*ntheta + mod(i,ntheta) + 1
    Tr=NTuple{3,Int32}[]
    for j in 1:nlevels-1, i in 0:ntheta-1                            # wall (outward)
        a=idx(j,i);b=idx(j,i+1);cc=idx(j+1,i+1);d=idx(j+1,i)
        push!(Tr,(Int32(a),Int32(b),Int32(cc))); push!(Tr,(Int32(a),Int32(cc),Int32(d)))
    end
    for i in 0:ntheta-1; push!(Tr,(Int32(ci),Int32(idx(1,i+1)),Int32(idx(1,i)))); end   # bottom cap
    for i in 0:ntheta-1; push!(Tr,(Int32(cti),Int32(idx(nlevels,i)),Int32(idx(nlevels,i+1)))); end # top cap
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
    raw = (ox0,ox1,oy0,oy1,z0,z1,ix0,ix1,iy0,iy1)
    names = ("ox0","ox1","oy0","oy1","z0","z1","ix0","ix1","iy0","iy1")
    vals = ntuple(i -> _finite_real(raw[i],"box_tunnel_surface",names[i]), 10)
    ox0,ox1,oy0,oy1,z0,z1,ix0,ix1,iy0,iy1 = vals
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

# ── hollow box (axis-aligned Boolean difference: outer box MINUS inner cavity) ──
"""
    box_shell_surface(ox0,ox1, oy0,oy1, oz0,oz1, ix0,ix1, iy0,iy1, iz0,iz1) -> Mesh

Closed boundary surface of the axis-aligned outer box `[ox0,ox1]×[oy0,oy1]×[oz0,oz1]`
MINUS a strictly-interior inner box `[ix0,ix1]×[iy0,iy1]×[iz0,iz1]` — a hollow box
(cavity / "case shell"), the native axis-aligned CSG difference. The inner box must
lie strictly inside the outer one (validated).

The boundary has **two** components: the outer box surface (12 outward triangles,
the [`box_surface`](@ref) pattern) and the inner cavity surface (the same pattern
with **reversed** triangle winding, so its normals face the cavity interior — i.e.
point out of the solid shell). Together they form a closed, manifold, 2-component
surface with zero open edges, ready for [`Mesh3D.tetrahedralize`](@ref) — which
fills the shell and reports the exact volume `outer − inner`.
"""
function box_shell_surface(ox0::Real,ox1::Real, oy0::Real,oy1::Real, oz0::Real,oz1::Real,
                           ix0::Real,ix1::Real, iy0::Real,iy1::Real, iz0::Real,iz1::Real)
    raw = (ox0,ox1,oy0,oy1,oz0,oz1,ix0,ix1,iy0,iy1,iz0,iz1)
    names = ("ox0","ox1","oy0","oy1","oz0","oz1","ix0","ix1","iy0","iy1","iz0","iz1")
    vals = ntuple(i -> _finite_real(raw[i],"box_shell_surface",names[i]), 12)
    ox0,ox1,oy0,oy1,oz0,oz1,ix0,ix1,iy0,iy1,iz0,iz1 = vals
    (ox0<ix0<ix1<ox1 && oy0<iy0<iy1<oy1 && oz0<iz0<iz1<oz1) ||
        throw(ArgumentError("box_shell_surface: inner box must be strictly inside the outer"))
    Co = Float64[ox0 ox1 ox1 ox0 ox0 ox1 ox1 ox0; oy0 oy0 oy1 oy1 oy0 oy0 oy1 oy1; oz0 oz0 oz0 oz0 oz1 oz1 oz1 oz1]
    Ci = Float64[ix0 ix1 ix1 ix0 ix0 ix1 ix1 ix0; iy0 iy0 iy1 iy1 iy0 iy0 iy1 iy1; iz0 iz0 iz0 iz0 iz1 iz1 iz1 iz1]
    C = hcat(Co, Ci)                                     # inner vertices are 9..16
    Fo = [(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),
          (2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]   # outer: outward normals
    Fi = [(a+8, c+8, b+8) for (a,b,c) in Fo]             # inner: +8 offset, reversed winding
    F = vcat(Fo, Fi)
    t = Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    return Mesh(C; tris=t)
end

@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline function _unit(a)
    l=sqrt(a[1]^2+a[2]^2+a[3]^2)
    (isfinite(l) && l > 0) || throw(ArgumentError("Geometry: axis must have finite positive length"))
    (a[1]/l,a[2]/l,a[3]/l)
end

end # module Geometry
