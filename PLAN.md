# Tessella.jl — Architecture & Rigorous Port Plan

> A Julia-native, robust, memory-efficient mesh generator. Goal: **design → mesh
> always works** (the HFSS experience), so ASCENT/ascent-studio never fail at the
> geometry→mesh step. Informed by a study of the gmsh 4.13.1 source, but a
> *rewrite* — not a transliteration — in idiomatic, optimized Julia.

---

## 1. Honest scope (read first)

Porting gmsh literally is **not** the plan, and pretending otherwise would be
dishonest:

- gmsh 4.13.1 **core** (`src/`) is **~340k LOC of C++** across 759 files.
- With its bundled `contrib/` (OpenCASCADE geometry kernel, TetGen, Netgen, HXT,
  blossom, voro++, …) the full toolkit is **>1.5M LOC**.
- The geometry kernel gmsh *uses* is **OpenCASCADE (OCC)** — millions of lines of
  BREP/NURBS/Boolean C++. **A pure-Julia OCC replacement is explicitly out of
  near-term scope.** It would be its own multi-year project.

What gmsh actually *is* = **a geometry model (`GModel`) + meshing algorithms + an
OCC interface + I/O**. The pain we hit in the ASCENT campaign (the enclosure
"overlapping facets / no elements in volume") lives in the **meshing algorithms
(3-D boundary recovery)**, not in OCC. That layer we **can and will own** in
Julia.

**Therefore Tessella stages the work:** own the meshing pipeline + robustness
first (ASCENT's real need), consume geometry from OCC via interop initially, and
grow toward gmsh feature parity and (optionally, much later) a native geometry
kernel. Every stage ships something ASCENT can use.

## 2. gmsh architecture as studied (4.13.1 source)

`src/` modules and what each is (LOC approximate, from the extracted source):

| module | LOC | role | port priority |
|---|---:|---|---|
| `numeric` | 40k | robust predicates, Gauss rules, bases, `Numeric.cpp` | **core** (predicates only; bases belong to ASCENT) |
| `geo` | 83k | `GModel`, `GEntity/GVertex/GEdge/GFace/GRegion`, OCC interface (`GModelIO_OCC`), many `GModelIO_*` writers | **core** (model + MSH/geo I/O; skip exotic formats) |
| `mesh` | 60k | the algorithms — see below | **core** (the heart) |
| `common` | — | context, options, logging, OS glue | partial (infra only) |
| `parser` | — | `.geo` scripting language (lex/yacc) | later (compat) |
| `fltk`,`graphics`,`post`,`plugin`,`solver` | — | GUI, post-processing, FEM | **out of scope** (ASCENT is the solver) |

The meshing pipeline (`src/mesh/Generator.cpp`) is a strict dimensional cascade:

```
0D vertices → 1D edges (meshGEdge) → 2D faces (meshGFace: Delaunay insertion / BDS /
frontal / transfinite / pack) → 3D regions (meshGRegion: Delaunay insertion +
boundary recovery + HXT + local mesh modification) → high order → optimize
```

**Robustness-critical files** (where gmsh fails on hard geometry):
- `meshGRegionBoundaryRecovery.cpp` — recovering the 2-D boundary triangulation
  as constrained facets of the 3-D tetrahedralization. **This is exactly the
  "Invalid boundary mesh (overlapping facets) on surface 86" failure** on the
  enclosure coax feed-through.
- `delaunay3d.cpp`, `meshGRegionDelaunayInsertion.cpp` — 3-D Delaunay kernel.
- `meshGRegionLocalMeshMod.cpp`, `meshGRegionMMG.cpp` — sliver removal / repair.
- `numeric` robust predicates — `orient3d`/`insphere` must be *exact* (adaptive)
  or the kernel produces the inconsistencies that surface as overlapping facets.

Bundled fallback meshers we studied as references: **HXT** (fast parallel 3-D
Delaunay), **Netgen** (robust OCC advancing-front), **MMG** (remesh/repair),
**blossom** (quad recombination), **voro++**, **Revoropt** (CVT).

## 3. Tessella architecture (Julia)

Idiomatic, type-stable, allocation-conscious Julia. No transliterated C++.

```
Tessella
├── Predicates    exact adaptive orient2/orient3/incircle/insphere (Shewchuk), SoS
├── Geometry      GModel, GVertex/GEdge/GFace/GRegion; kernel interop (Stage 5)
├── MeshTypes     half-facet / compact SoA mesh; Node, Tri, Tet; topology queries
├── SizeField     background mesh, curvature-, distance-, and boundary-driven sizing
├── Mesh1D        edge meshing under a size field
├── Mesh2D        Delaunay + constrained Delaunay (CDT) + Ruppert/Chew refinement + frontal
├── Mesh3D        Delaunay + constrained + ROBUST BOUNDARY RECOVERY + sliver handling
├── Optimize      Laplacian/ODT smoothing, edge/face swaps, sliver removal
├── Heal          coincident-face / sliver / thin-feature detection + repair
├── IO            .msh (v2/v4) read/write, gmsh .geo compat reader, STL
└── Fallback      staged strategies + optional external backend shim (transition only)
```

**Design principles (the reason we do this at all):**
1. **Exactness where it counts.** Adaptive exact geometric predicates + Simulation
   of Simplicity. Floating-point predicate inconsistency is *the* root of
   overlapping-facet failures. This is non-negotiable and gets an exact-arithmetic
   oracle.
2. **Robust boundary recovery.** Constrained 3-D Delaunay with a recovery scheme
   that provably restores every boundary facet (segment/facet recovery via
   flips + Steiner points), so a valid input surface **always** yields a valid
   volume mesh. The enclosure coax-junction is the standing acceptance test.
3. **Heal, don't fail.** Detect near-coincident faces, zero-thickness slivers,
   and multi-way junctions up front; imprint/merge within tolerance; never emit
   an invalid BREP to the mesher.
4. **Always-valid or explicit blocker.** Every run either returns a validated
   mesh (all regions filled, positive volumes, manifold boundary, Euler check) or
   a precise diagnostic — never a silent empty volume.
5. **Memory efficiency.** Compact SoA, in-place kernels, streaming I/O; avoid the
   pointer-chasing `MVertex*` graph gmsh uses.

## 4. Staged roadmap

Each stage follows the repo's core loop: **spec → independent oracle → implement
→ CRC test → benchmark → gmsh cross-check**. No stage is "done" without a
mutation-sensitive test suite and a CRC-stamped regression artifact.

- **Stage 0 — Foundations.** Repo, CI, CRC discipline, mesh data structures,
  `.msh` v2/v4 read+write (cross-checked against gmsh output), **exact predicates
  with an exact-rational oracle**. *Exit:* round-trip a gmsh `.msh`; predicates
  pass exhaustive degenerate-configuration tests.
- **Stage 1 — 2-D core.** Delaunay (Bowyer–Watson / incremental), constrained
  Delaunay (CDT) for PSLGs, quality refinement (Ruppert + Chew). *Exit:* mesh
  arbitrary planar straight-line graphs; angle/area guarantees; CRC vs analytic.
- **Stage 2 — Surfaces & edges.** 1-D edge meshing under a size field; 2-D
  meshing of parametric/BREP faces (surface Delaunay in parameter space with the
  metric). *Exit:* mesh the flat/patch/coax faces from the HFSS cases.
- **Stage 3 — 3-D core + boundary recovery (the robustness milestone).** 3-D
  Delaunay, constrained tetrahedralization, **robust boundary recovery**, sliver
  handling. *Exit:* **mesh the enclosure (9.2) coax feed-through that gmsh
  4.13.1/4.15.2 cannot**, all volumes filled, validated.
- **Stage 4 — Sizing & optimization.** Background/curvature/distance size fields;
  smoothing, swaps, sliver removal; quality histograms. *Exit:* match/beat gmsh
  quality on the 22 HFSS geometries.
- **Stage 5 — Geometry kernel.** OCC interop (Booleans → BREP → faces) with the
  `Heal` layer, or a native CSG path for the primitives ASCENT emits. *Exit:*
  ingest ASCENT `solid_model` geometry directly.
- **Stage 6 — Integration & parity.** High-order (curved) elements; drop-in
  replacement for ASCENT's `gmsh -3` call; **regression across all 22 HFSS
  UserGuide cases** with CRC-stamped meshes. *Exit:* ASCENT solves every case
  through Tessella.

## 5. CRC & verification discipline (mandatory)

Applies to every function (see `DEVELOPMENT.md`):
- **Independent oracle** for each algorithm: exact-rational predicate, analytic
  mesh, Euler/Delaunay invariants, conservation of boundary, or cross-check vs a
  reference gmsh mesh.
- **Mutation-sensitive tests** (a sign flip / swapped index / dropped Steiner
  point must fail a test).
- **CRC-stamped regression artifacts**: every accepted mesh writes a checksum
  (node/element counts, bbox, quality stats, boundary hash) recorded in
  `STATUS.md`; changes require a justified diff.
- **Bug-first**: suspicious success (e.g. "no empty volumes" with zero tets — the
  exact false positive we hit with `occ.healShapes`) is a bug until disproven by
  checking the *actual* element counts.
- No weakened tolerances, no skipped degenerate cases, no `@test true`.

## 6. Non-goals (near-term, stated so scope stays honest)

- Pure-Julia OpenCASCADE / NURBS CAD kernel.
- GUI, post-processing, FEM assembly/solve (that is ASCENT), and the long tail of
  `GModelIO_*` formats (CGNS, MED, NIfTI, …).
- Beating HXT on raw parallel throughput before correctness/robustness is proven.

## 7. First concrete target (bridges the campaign)

The **enclosure coax feed-through** (surface 86: shield/case/air junction) is the
standing Stage-3 acceptance case: gmsh 4.13.1 and 4.15.2 both fail it; Tessella's
robust boundary recovery must mesh it with all three volumes (pin/case/air)
filled and validated. That single case is the proof that this project earns its
existence.
