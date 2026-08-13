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
# Highest development stage whose CRC gate (STATUS.md) is green. 3 ⇒ Stage 2
# (1-D + surface meshing) complete; Stage 3 (3-D kernel + volume filling) + Stage 4
# (quality/smoothing) + Stage 5 (heal detection) landed incrementally on top.
const TESSELLA_STAGE = 3  # see STATUS.md stage board

# ── Submodules (PLAN.md §3) ────────────────────────────────────────────────────
include("Predicates.jl")     # Stage 0: adaptive exact orient/incircle/insphere + SoS
include("ExactMesh3D.jl")    # Stage 3: exact-coordinate (Rational{BigInt}) 3-D Delaunay
include("MeshTypes.jl")      # Stage 0: compact SoA mesh, topology, quality, CRC checksum
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

using .MeshTypes: Mesh, validate, mesh_crc
using .Mesh2D: constrained_delaunay, refine!, classify_interior, to_mesh
using .Mesh3D: tetrahedralize, tetrahedralize_multi, tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion, recover_boundary, mesh_boolean, mesh_sized_conforming, mesh_cylinder
using .RecoverCDT: recover_boundary_cdt, mesh_sized_cdt
using .Optimize: smooth_laplacian, smooth_odt, smooth_optimize, mesh_quality
using .Heal: is_meshable

export mesh_volume, mesh_planar, stage
# curated re-exports of the public API
export Mesh, validate, mesh_crc, mesh_quality, is_meshable
export tetrahedralize, tetrahedralize_multi, tetrahedralize_conforming, tetrahedralize_conforming_exact, tets_per_region, mesh_box, mesh_box_regions, BoxRegion, recover_boundary, recover_boundary_cdt, mesh_sized_cdt, mesh_boolean, mesh_sized_conforming, mesh_cylinder, smooth_laplacian, smooth_odt, smooth_optimize

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
    diag = validate(m)
    diag.ok || throw(ErrorException("mesh_planar: produced an invalid mesh — " * join(diag.messages, "; ")))
    return m
end

"""
    mesh_volume(surface; smooth=true, smooth_iters=5, rng_seed=1, check=true) -> Mesh

Robust volume mesh of the closed triangulated `surface`, with the "always-valid or
explicit blocker" contract: the surface is first screened by [`is_meshable`](@ref)
(unless `check=false`); a defect raises an `ArgumentError` carrying the precise
[`Heal.SurfaceReport`](@ref) — never a silent bad mesh. Otherwise it is filled
with tetrahedra ([`tetrahedralize`](@ref)) and, if `smooth`, Laplacian-smoothed
([`smooth_laplacian`](@ref)). The result is `validate`-checked before return.
"""
function mesh_volume(surface::Mesh; smooth::Bool=true, smooth_iters::Integer=5,
                     optimize::Bool=false, rng_seed::Integer=1, check::Bool=true)
    if check
        ok, report = is_meshable(surface)
        ok || throw(ArgumentError("mesh_volume: input surface is not meshable — " *
                                  join(report.messages, "; ")))
    end
    m = tetrahedralize(surface; rng_seed=rng_seed, optimize=optimize)
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
