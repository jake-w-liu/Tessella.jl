"""
    Tessella

A Julia-native, robust, memory-efficient mesh generator. Goal: *design → mesh
always works*. See `PLAN.md` for the staged architecture, `DEVELOPMENT.md` for the
mandatory CRC (Correctness–Robustness–Completeness) discipline, and `STATUS.md` for
the live stage board.

Pipeline (0-D → 3-D): exact predicates ([`Predicates`](@ref)) underpin a robust
2-D Delaunay/CDT/refinement core ([`Mesh2D`](@ref)), 1-D edge + surface meshing
([`Mesh1D`](@ref), [`MeshSurface`](@ref)) under size fields ([`SizeField`](@ref)),
and a 3-D Delaunay kernel + volume filling ([`Mesh3D`](@ref)), with quality/
optimization ([`Optimize`](@ref)), surface healing ([`Heal`](@ref)), and gmsh
`.msh` / STL I/O ([`IO`](@ref)). The top-level [`mesh_volume`](@ref) ties them
together with the *validated-or-explicit-blocker* contract.
"""
module Tessella

# ── Version / capability banner ────────────────────────────────────────────────
# Highest development stage whose CRC gate (STATUS.md) is green. Stage 6 includes
# high-order elements and solver-consumable I/O; stages 0–5 are prerequisites.
const TESSELLA_STAGE = 6  # see STATUS.md stage board

# ── Submodules (PLAN.md §3) ────────────────────────────────────────────────────
include("Predicates.jl")     # Stage 0: adaptive exact orient/incircle/insphere + SoS
include("MeshTypes.jl")      # Stage 0: compact SoA mesh, topology, quality, CRC checksum
include("ExactMesh3D.jl")    # Stage 3: exact-coordinate (Rational{BigInt}) 3-D Delaunay
include("IO.jl")             # Stage 0: .msh v2/v4 read/write, STL, .geo scan
include("Mesh2D.jl")         # Stage 1: 2-D Delaunay + CDT + Ruppert refinement
include("SizeField.jl")      # Stage 2/4: size fields
include("Mesh1D.jl")         # Stage 2: graded 1-D edge meshing
include("MeshSurface.jl")    # Stage 2: planar / cylinder / parametric surface meshing
include("Mesh3D.jl")         # Stage 3: 3-D Delaunay + volume filling (+ multi-region)
include("RecoverCDT.jl")     # Stage 3: general conforming-Delaunay boundary recovery (exact)
include("Optimize.jl")       # Stage 4: tet quality report + Laplacian smoothing
include("Heal.jl")           # Stage 5: surface-defect detection ("heal, don't fail")
include("Geometry.jl")       # Stage 5: native constructive primitive surfaces
include("CAD.jl")            # Stage 5: native analytical geometry (surfaces + exact imprints), no OCC
include("HighOrder.jl")      # Stage 6: quadratic (P2) tet generation + type-11 I/O

using .MeshTypes: Mesh, validate, mesh_crc, tet_signed_volume, boundary_faces
using .Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using .Mesh3D: tetrahedralize, tetrahedralize_multi, tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion, recover_boundary, mesh_boolean, mesh_sized_conforming, mesh_cylinder, refine_to_size
using .RecoverCDT: recover_boundary_cdt, recover_partition_cdt, mesh_sized_cdt
using .Optimize: smooth_laplacian, smooth_odt, smooth_optimize, remove_slivers, mesh_quality
using .Heal: is_meshable

# Install Mesh3D's exact-recovery extension only after both modules and their
# public types are available, avoiding a Mesh3D ↔ RecoverCDT include cycle.
Mesh3D._recover_boundary_exact(surface::Mesh) = recover_boundary_cdt(surface)
Mesh3D._recover_partition_exact(surfaces::AbstractVector{Mesh}) = recover_partition_cdt(surfaces)

export mesh_volume, mesh_planar, mesh_sized_extrude, mesh_sized, refine_to_size, stage
# curated re-exports of the public API
export Mesh, validate, mesh_crc, mesh_quality, is_meshable
export tetrahedralize, tetrahedralize_multi, tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion, recover_boundary, recover_boundary_cdt, recover_partition_cdt, mesh_sized_cdt, mesh_boolean, mesh_sized_conforming, mesh_cylinder, smooth_laplacian, smooth_odt, smooth_optimize, remove_slivers

"""
    stage() -> Int

Current implemented development stage (see the `STATUS.md` stage board).
"""
stage() = TESSELLA_STAGE

"""
    mesh_planar(xs, ys, segments; min_angle_deg=25.0, max_area=Inf, rng_seed=1) -> Mesh

Quality 2-D triangle mesh of the planar straight-line graph (points `(xs,ys)` +
constraint `segments`, `(i,j)` index pairs): constrained-Delaunay triangulate,
Ruppert-refine to the angle/area bound, keep the interior of the constrained
domain, and return a validated 2-D [`Mesh`](@ref) (nodes carry `z = 0`). The
domain boundary must be closed constrained loops. This is the 2-D counterpart of
[`mesh_volume`](@ref).
"""
function mesh_planar(xs::Vector{Float64}, ys::Vector{Float64},
                     segments::AbstractVector{<:Tuple{Integer,Integer}};
                     min_angle_deg::Real=25.0, max_area::Real=Inf, rng_seed::Integer=1)
    T = constrained_delaunay(xs, ys, segments; rng_seed=rng_seed)
    interior = refine!(T; min_angle_deg=min_angle_deg, max_area=max_area)
    m = to_mesh(T; interior=interior)
    size(m.tris, 2) > 0 ||
        throw(ArgumentError("mesh_planar: constrained boundary encloses no triangles"))
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_planar: produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

"""
    mesh_sized_extrude(xs, ys, segments, z0, z1; hmax, min_angle_deg=25.0) -> Mesh

**Uniform size-controlled** tet mesh of the prismatic (extruded) domain whose
cross-section is the polygon `(xs, ys)` with boundary/hole `segments`, extruded from
`z=z0` to `z=z1`, with a **guaranteed maximum edge length `≤ hmax`**. The cross-section
is meshed by the size-controlled 2-D Ruppert mesher (edges `≤ hmax/√2`), then extruded
into `⌈(z1−z0)/(hmax/√2)⌉` layers; each triangular prism is split into three tets by a
**column-global-index diagonal rule**, so every shared quad face picks the same diagonal
in both incident prisms — the result is **conforming** (no face shared by >2 tets),
watertight, provably valid (each tet oriented positive), and of **exact volume**
(area(polygon)·(z1−z0), boundary preserved — no jitter). The `√2` factor bounds the
prism face-diagonal `√(edge²+layer²) ≤ hmax`.

Handles **non-convex** cross-sections and polygons with holes (via the `segments`
constrained-Delaunay + interior classification). This is the extruded/prismatic case of
uniform sizing on an arbitrary domain — complementing [`mesh_box`](@ref) (box) and
[`mesh_cylinder`](@ref) (cylinder). For a general closed faceted surface, use
[`mesh_sized`](@ref).
"""
function mesh_sized_extrude(xs::AbstractVector, ys::AbstractVector,
                            segments::AbstractVector{<:Tuple{Integer,Integer}},
                            z0::Real, z1::Real; hmax::Real, min_angle_deg::Real=25.0)
    hm=try Float64(hmax) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: hmax must be Float64-representable: $(sprint(showerror,err))"))
    end
    za=try Float64(z0) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: z0 must be Float64-representable: $(sprint(showerror,err))"))
    end
    zb=try Float64(z1) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: z1 must be Float64-representable: $(sprint(showerror,err))"))
    end
    (isfinite(hm)&&hm>0) || throw(ArgumentError("mesh_sized_extrude: hmax must be finite and positive (got $hmax)"))
    (isfinite(za)&&isfinite(zb)&&zb>za) ||
        throw(ArgumentError("mesh_sized_extrude: need finite z1 > z0 (got z0=$z0, z1=$z1)"))
    h2 = hm/sqrt(2.0)                                      # prism face-diagonal ≤ hmax
    (isfinite(h2)&&h2>0) || throw(ArgumentError("mesh_sized_extrude: hmax is below Float64 spacing resolution"))
    T2 = constrained_delaunay(Float64.(xs), Float64.(ys), segments)
    interior = refine!(T2; min_angle_deg=min_angle_deg, max_area=0.5*h2*h2, size=(cx,cy)->h2)
    m2 = to_mesh(T2; interior=interior)                    # 2-D mesh: coords (z=0), tris
    nv2 = size(m2.coords, 2)
    nv2 >= 3 && size(m2.tris, 2) >= 1 ||
        throw(ErrorException("mesh_sized_extrude: cross-section produced no interior triangles (check segments orientation/closure)"))
    ratio=(zb-za)/h2
    (isfinite(ratio)&&ratio<=typemax(Int)) ||
        throw(ArgumentError("mesh_sized_extrude: axial layer count exceeds the platform Int limit"))
    nz=max(1,ceil(Int,ratio));dz=(zb-za)/nz
    levels=try Base.checked_add(nz,1) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: axial layer count overflows Int"))
    end
    nout=try Base.checked_mul(levels,nv2) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: node count overflows Int"))
    end
    nout<=typemax(Int32) || throw(ArgumentError("mesh_sized_extrude: $nout nodes exceed Int32"))
    coords = Matrix{Float64}(undef, 3, nout)
    @inbounds for k in 0:nz, i in 1:nv2
        id = k*nv2 + i
        coords[1,id] = m2.coords[1,i]; coords[2,id] = m2.coords[2,i]; coords[3,id] = za + k*dz
    end
    _n(id) = @inbounds (coords[1,id], coords[2,id], coords[3,id])
    ntout=try Base.checked_mul(3,Base.checked_mul(size(m2.tris,2),nz)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized_extrude: tetrahedron count overflows Int"))
    end
    tets = Matrix{Int32}(undef, 4, ntout); col = 0
    @inbounds for ti in 1:size(m2.tris,2)
        v = sort!(Int[m2.tris[1,ti], m2.tris[2,ti], m2.tris[3,ti]])   # v[1]<v[2]<v[3] (column-index rule)
        for k in 0:nz-1
            b1=k*nv2+v[1]; b2=k*nv2+v[2]; b3=k*nv2+v[3]
            t1=(k+1)*nv2+v[1]; t2=(k+1)*nv2+v[2]; t3=(k+1)*nv2+v[3]
            for te in ((b1,b2,b3,t3), (b1,b2,t3,t2), (b1,t2,t3,t1))
                a,b,c,d = te
                tet_signed_volume(_n(a),_n(b),_n(c),_n(d)) < 0 && ((c,d) = (d,c))   # orient positive
                col += 1; tets[1,col]=a; tets[2,col]=b; tets[3,col]=c; tets[4,col]=d
            end
        end
    end
    m = Mesh(coords; tets=tets)
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_sized_extrude: produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

"""
    mesh_sized(surface::Mesh; hmax) -> Mesh

**General uniform size-controlled** tet mesh of the domain enclosed by an arbitrary
closed triangulated `surface` (curved or polyhedral, convex or non-convex), with a
**guaranteed maximum edge length `≤ hmax`**. Two stages:

1. **Fill** — [`tetrahedralize`](@ref) (Delaunay + ray-cast interior classification)
   gives a valid conforming tet mesh of the domain; if that is not watertight (a
   non-convex boundary the Delaunay restriction misses), [`recover_boundary`](@ref) is
   used instead.
2. **Size** — [`refine_to_size`](@ref) longest-edge-bisects the fill until every edge is
   `≤ hmax`. Bisection is conforming, validity/volume-preserving, and terminating, so the
   result is a valid, watertight, conforming mesh with `maxedge ≤ hmax` at the domain's
   own (faceted-surface) volume.

Verified on convex boxes, spheres, conical frusta, and a non-convex L-prism (guaranteed
`maxedge ≤ hmax`, valid, watertight, conforming). Returns the mesh, or throws an explicit
blocker if the surface cannot be filled watertight (never a silently bad mesh). For
axis-aligned boxes prefer [`mesh_box`](@ref) and for extrusions [`mesh_sized_extrude`](@ref)
(exact boundary, no faceting); `mesh_sized` is the general fallback for arbitrary surfaces.
"""
function mesh_sized(surface::Mesh; hmax::Real)
    hm=try Float64(hmax) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("mesh_sized: hmax must be Float64-representable: $(sprint(showerror,err))"))
    end
    (isfinite(hm)&&hm>0) || throw(ArgumentError("mesh_sized: hmax must be finite and positive (got $hmax)"))
    m0 = _conforming_fill(surface; caller="mesh_sized")
    m = refine_to_size(m0, hm)
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_sized: refinement produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

# Fill a closed PLC without ever accepting a convex-hull cap as its boundary.
# `tetrahedralize` owns the restriction, star-shaped fan, exact-CDT recovery, and
# bounded Float64 fallback, and returns only after Mesh3D's exact geometric gate.
function _conforming_fill(surface::Mesh; caller::AbstractString,
                          rng_seed::Integer=1, optimize::Bool=false)
    try
        return tetrahedralize(surface; rng_seed=rng_seed, optimize=optimize)
    catch err
        err isa InterruptException && rethrow()
        (err isa ArgumentError || err isa ErrorException) || rethrow()
        throw(ErrorException("$caller: could not construct a conforming volume mesh — " *
                             sprint(showerror, err)))
    end
end

"""
    mesh_volume(surface; smooth=true, smooth_iters=5, rng_seed=1, check=true) -> Mesh

Robust volume mesh of the closed triangulated `surface`, with the "always-valid or
explicit blocker" contract: the surface is first screened by [`is_meshable`](@ref)
(unless `check=false`); a defect raises an `ArgumentError` carrying the precise
[`Heal.SurfaceReport`](@ref) — never a silent bad mesh. Otherwise an initial
[`tetrahedralize`](@ref) fill must pass the exact PLC-boundary conformity gate; a
non-conforming convex-hull cap falls through [`recover_boundary`](@ref) and then
[`recover_boundary_cdt`](@ref). If `smooth`, the conforming fill is Laplacian-smoothed
([`smooth_laplacian`](@ref)). The result is `validate`-checked before return.
"""
function mesh_volume(surface::Mesh; smooth::Bool=true, smooth_iters::Integer=5,
                     optimize::Bool=false, rng_seed::Integer=1, check::Bool=true)
    if check
        ok, report = is_meshable(surface)
        ok || throw(ArgumentError("mesh_volume: input surface is not meshable — " *
                                  join(report.messages, "; ")))
    end
    m = check ? _conforming_fill(surface; caller="mesh_volume", rng_seed=rng_seed,
                                 optimize=optimize) :
                tetrahedralize(surface; rng_seed=rng_seed, optimize=optimize, check=false)
    smooth && (m = smooth_laplacian(m; iters=smooth_iters))
    # optimization-based (min-dihedral) smoothing targets slivers the mean-smoother and
    # topological flips leave behind; only in the optimize path (it is a local search).
    optimize && (m = smooth_optimize(m; iters=smooth_iters))
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_volume: produced an invalid mesh — " *
                                    join(diag.messages, "; ")))
    return m
end

end # module Tessella
