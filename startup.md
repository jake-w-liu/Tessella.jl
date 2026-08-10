# Tessella.jl — Startup (open a fresh session from this file)

You are starting development of **Tessella.jl**, a Julia-native robust mesh
generator that replaces gmsh in the ASCENT / ascent-studio toolchain. This file
is your entry point. Read it, then read `PLAN.md`, `DEVELOPMENT.md`, `STATUS.md`
in that order before writing code.

## What this project is (and is not)

- **Is:** a from-scratch, idiomatic, optimized Julia mesh generator — geometry
  model + 1-D/2-D/3-D meshing + **robust boundary recovery** + healing + `.msh`
  I/O — informed by the gmsh source but a rewrite, not a transliteration.
- **Is not:** a literal gmsh port, and **not** a reimplementation of OpenCASCADE.
  Scope is stated honestly in `PLAN.md` §1. Do not silently expand it.

## Why it exists (the standing proof case)

gmsh 4.13.1 **and** 4.15.2-git both fail to mesh a shielded-enclosure coax
feed-through — `Invalid boundary mesh (overlapping facets) on surface 86`, all
volumes empty. It is **not** memory (peak 3.3 GB) and it resists every standard
remedy. The failing geometry is captured at
`test/fixtures/enclosure_coax_junction.geo`. **Tessella's Stage-3 acceptance test
is to mesh that file with all three volumes (pin/case/air) filled and validated.**

## Non-negotiable discipline (from `DEVELOPMENT.md`)

- Every function ships with an **independent oracle**, mutation-sensitive tests,
  and a **mesh-CRC checksum**. `spec → oracle → implement → CRC → reverify`.
- **Count actual elements.** "No empty volumes" is a false positive if there are
  no volumes (this exact trap sank `occ.healShapes()` — see `STATUS.md`). A valid
  volume mesh has tets in *every* region.
- Exact adaptive geometric predicates (Shewchuk + Simulation of Simplicity) are
  the foundation of 3-D robustness. Build and verify them first.
- No weakened tolerances, no `@test true`, no silent empty regions.

## Reference material

- **gmsh 4.13.1 source** (study, don't copy): download
  `https://gmsh.info/src/gmsh-4.13.1-source.tgz`. Key files to read:
  - `src/mesh/Generator.cpp` — the 0D→1D→2D→3D pipeline.
  - `src/mesh/meshGRegionBoundaryRecovery.cpp` — **the robustness core** (this is
    where the enclosure fails); study its recovery scheme.
  - `src/mesh/delaunay3d.cpp`, `meshGRegionDelaunayInsertion.cpp` — 3-D kernel.
  - `src/mesh/meshGFaceDelaunayInsertion.cpp`, `meshGFaceOptimize.cpp` — 2-D.
  - `src/numeric/` — predicates and robust numerics.
  - `contrib/hxt`, `contrib/Netgen` — reference algorithms (HXT = fast Delaunay;
    Netgen = robust advancing front). MMG = remesh/repair.
  - If a local copy exists under `reference/` (gitignored), use it directly.
- **gmsh 4.13.1 SDK** (Python API, handy for A/B cross-checks): downloadable from
  `https://gmsh.info/bin/macOS/gmsh-4.13.1-MacOSARM-sdk.tgz`.

## First task — Stage 0 (Foundations)

Work strictly to the `STATUS.md` stage board. Stage 0 deliverables, each behind
its CRC gate, in order:

1. **Exact predicates** (`src/Predicates.jl`): adaptive `orient2`, `orient3`,
   `incircle`, `insphere` + SoS. Oracle: exact `Rational{BigInt}`/`BigFloat`
   recomputation on an exhaustive set of degenerate configurations. This is the
   single most important correctness foundation — do it first, verify hardest.
2. **Mesh data structures** (`src/MeshTypes.jl`): compact SoA nodes + simplices,
   half-facet adjacency, topology queries (Euler check, boundary extraction).
3. **`.msh` I/O** (`src/IO.jl`): read/write MSH v2 and v4; round-trip a gmsh
   `.msh` byte-for-byte on connectivity; parse the fixture `.geo` enough to drive
   Stage 2+ (or ingest a boundary surface mesh to start).
4. Wire them into `Tessella`, add `test/` CRC suites, bump `stage()` to 1 only
   when the Stage-0 gate in `STATUS.md` is green.

Update `STATUS.md` (stage board + log + any CRC checksums) at the end of every
working session. Keep `main` green. Commit in small, oracle-justified steps.

## Housekeeping

- Local git is initialized with an initial commit.
- **GitHub remote is pending**: run `! gh auth login` once, then
  `gh repo create Tessella.jl --private --source=. --remote=origin --push`
  (or ask the assistant to do it after you authenticate).
- Julia ≥ 1.11. `julia --project=. -e 'using Pkg; Pkg.test()'` must stay green.
