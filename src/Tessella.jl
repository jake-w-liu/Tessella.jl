"""
    Tessella

A Julia-native, robust, memory-efficient mesh generator.

Goal: *design → mesh always works*. See `PLAN.md` for the staged architecture and
`DEVELOPMENT.md` for the mandatory CRC (Correctness–Robustness–Completeness)
discipline. Development is tracked in `STATUS.md`.

Nothing here is implemented yet — this is the Stage 0 skeleton created alongside
the project plan. Begin from `startup.md`.
"""
module Tessella

# ── Version / capability banner ────────────────────────────────────────────────
# Highest development stage whose CRC gate (STATUS.md) is green. 3 ⇒ Stage 2
# (1-D edge meshing + surface meshing) complete; 2 ⇒ Stage 1 (2-D); 1 ⇒ Stage 0.
const TESSELLA_STAGE = 3  # see STATUS.md stage board

# ── Submodules (filled per PLAN.md §3 as stages land) ──────────────────────────
include("Predicates.jl")
include("MeshTypes.jl")
include("IO.jl")
include("Mesh2D.jl")
include("SizeField.jl")
include("Mesh1D.jl")
include("MeshSurface.jl")
include("Mesh3D.jl")
include("Optimize.jl")

"""
    stage() -> Int

Current implemented development stage (0 = foundations only). Bumped as the
`STATUS.md` stage board advances and its CRC gate passes.
"""
stage() = TESSELLA_STAGE

# ── Submodule layout (to be filled per PLAN.md §3) ─────────────────────────────
# include("Predicates.jl")   # Stage 0: adaptive exact orient/incircle/insphere
# include("MeshTypes.jl")    # Stage 0: compact half-facet mesh + topology
# include("IO.jl")           # Stage 0: .msh v2/v4 read/write, gmsh .geo reader
# include("SizeField.jl")    # Stage 4
# include("Mesh1D.jl")       # Stage 2
# include("Mesh2D.jl")       # Stage 1/2
# include("Mesh3D.jl")       # Stage 3  (robust boundary recovery)
# include("Optimize.jl")     # Stage 4
# include("Heal.jl")         # Stage 5
# include("Geometry.jl")     # Stage 5

end # module Tessella
