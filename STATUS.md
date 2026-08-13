# Tessella.jl — Development Status Tracker

**Single source of truth for build status. Update every session. Nothing dropped.**

Package: `Tessella` · Julia ≥ 1.11 · goal: robust Julia-native mesh generator
replacing the gmsh dependency (see `PLAN.md`). CRC discipline mandatory
(`DEVELOPMENT.md`).

**This file tracks the Tessella package (the mesher) only.** All ASCENT-facing work —
the external FEM solver integration, the mesh/BC handshake, the solve proofs, and the
22-case HFSS solve campaign — is tracked separately in **`ASCENT.md`**.

## ✅ DONE vs ⬜ NOT DONE — at a glance (updated 2026-08-14)

**Suite green: 151,466 assertions** (`--check-bounds=yes`, Julia 1.12.6; 12m53s).
Everything below is committed to
`main`, source-only (HFSS/ASCENT data stays local, never pushed).

### ✅ DONE + verified
- **Stages 0–3 core** — exact predicates, 2-D & 3-D Delaunay/CDT, volume fill, boundary recovery
  for the *supported* classes **plus the exotic non-star+reflex twisted prism** (see below).
- **`recover_boundary_cdt` — GENERAL robust boundary recovery (2026-08-13, `src/RecoverCDT.jl`)**:
  exact-kernel conforming-Delaunay refinement closes the ONE class the Float64 `recover_boundary`
  couldn't — the **non-star + reflex twisted prism** — as a valid, closed-manifold,
  **exactly-conforming** tet mesh (boundary area == surface area & exact volume, verified), and
  also recovers box/genus-1 tunnel/hollow shell/faceted cylinder. **Item #7 is now DONE.**
- **Stage 4** — size control for **boxes + cylinders**; sliver removal (2-3/3-2 flips +
  `smooth_optimize` min-dihedral optimization-based smoothing).
- **Stage 5** — native CSG (`mesh_boolean`, `mesh_box_regions`). OCC-*library* interop = PLAN §1/§6 non-goal.
- **Stage 6** — P2 curved elements + gmsh MSH v2/v4 I/O.
- **Exact-coordinate `Rational{BigInt}` 3-D Delaunay kernel** (the keystone) + fast Bareiss
  determinant + `tetrahedralize_conforming_exact` (meshes cospherical assemblies the Float64 path can't).
- **Code fully optimized** — grid classifier (~150× on the documented bottleneck), Bareiss exact det.
- **Uniform sizing across the board** — box / cylinder / prismatic (`mesh_sized_extrude`) /
  **general arbitrary surface** (`mesh_sized` + `refine_to_size` longest-edge bisection), all with a
  **guaranteed `maxedge ≤ hmax`** (valid, watertight, conforming; exact volume for polyhedra).
- **All 22 HFSS geometry classes meshed natively** (from scratch, no gmsh/OCC — `validation/hfss_cases/`).
  *(ASCENT solving of these is tracked in `ASCENT.md`.)*
- **Representative capability meshes**: THIN-SLOT, ARRAY-PML, SPIRAL classes (valid + conforming).
- **Correctness** — independent adversarial audit (4 real bugs fixed) + an `orient3_sos` coplanar
  bug fixed; all predicates oracle-verified.

### ⬜ NOT DONE — package (honest)
_Updated 2026-08-13: **every Tessella package (meshing/geometry) item is DONE.** #7, #8, #9
(research flagships + literal ENC-COAX), #11 (Delaunay cospherical perf), and the meshing
half of #12 (native meshes for all 22 HFSS geometry classes) are all complete and
regression-pinned. The only work left in the whole project is the **external ASCENT
full-wave solve campaign**, tracked in `ASCENT.md` — not a Tessella package item._
| # | item | status |
|---|---|---|
| 7 | non-star+reflex recovery (twisted prism) | ✅ **DONE (2026-08-13)** — closed by `recover_boundary_cdt` (exact-kernel conforming-Delaunay refinement, Gabriel-encroachment off the exact DT with exact-rational boundary-Steiner points); twisted prism conforms (86 tets, exact area+volume, regression-pinned) |
| 8 | arbitrary-surface uniform sizing | ✅ **DONE (2026-08-13)** — `mesh_sized` + `refine_to_size` give **guaranteed `maxedge ≤ hmax`** on ANY closed surface (fill → longest-edge bisection); plus `mesh_box`/`mesh_cylinder`/`mesh_sized_extrude` for box/cylinder/prismatic and `mesh_sized_cdt`/`mesh_sized_conforming` for exact/interior sizing. Verified on box/sphere/frustum/L-prism (valid + watertight + conforming, exact volume for polyhedra), regression-pinned |
| 9 | *literal* ENC-COAX geometry | ✅ **DONE natively (2026-08-13)** — no OpenCASCADE. The `.geo` is fully parametric (Box/Cylinder primitives): Tessella parses the literal primitives and reconstructs the complete physical-group structure natively — 4 volumes (pin/slot/air/case at the exact fixture dimensions) + 5 tagged BC surfaces (radiation, pin/case PEC skins, resistor, p1 port) — conforming + valid + all-filled (`validation/enclosure_literal/`, regression-pinned) — the geometry gmsh leaves empty. The coax **bore imprint** (bore cuts the case wall — an OCC `BooleanFragments` op) is computed natively via `mesh_boolean` (mesh-CSG). **Exact curved geometry is native too**: a native analytical-CAD layer (`src/CAD.jl`) gives exact analytical surfaces (plane/cylinder/sphere/disk) + exact boolean **imprint curves** (a cylinder piercing a wall = an exact circle, nodes on both surfaces to ~1e-15) + exact projection; the pin is P2-curved onto its exact analytical cylinder (nodes on the true surface to 1.85e-17), bore imprint circle exact (`cad_test.jl`). (A *general* freeform NURBS/BREP kernel for arbitrary sculpted surfaces remains a larger effort, but every surface the enclosure uses is now exact.) |
| 11 | Delaunay cospherical perf | ✅ **FIXED (2026-08-13)** — both kernels | **The earlier "cavity retriangulation, not point location" diagnosis was measured-wrong.** Re-measured (three ways): on the fine pin (nθ=12,nz=16) Σlocate 43.9 s vs Σcavity 4.2 s vs Σretri 0.006 s — location dominates. Mechanism: a handful of far points burn the `16·ntets` walk budget then the O(n) `_locate_scan`, ×expensive exact predicates ⇒ O(n²). **Float64 fix (`Mesh3D.jl`, output-identical):** `_pick_start3` jump-and-walk (start from the nearest of ~∛n sampled landmark vertices via the `vtet` hint) + walk guard capped at `48·∛n`; byte-identical (`.sha`-verified over a 120-case fuzz + full suite). **nz=16 97.8→5.6 s, nz=24 300→7.4 s, nz=40 262→22.4 s.** **Exact-kernel fix:** a conservative Float64 circumsphere pre-filter in `delaunay3d_exact` (skip the expensive `insphere_rat` on tets the query is safely outside), output-identical (byte-identical to full-scan + `is_delaunay_exact` oracle), **19.5–41× faster**. Regression-pinned |
| 12 (meshing half) | native meshes for all 22 HFSS geometry classes | ✅ **DONE (2026-08-13)** — Tessella meshes **all 22** HFSS UG ch.5–10 case geometry classes natively from scratch (no gmsh/OCC): monopole, conical horn, patches, SAR sphere, CPW bowtie, endfire cell, magic tee, coax bend/stub, ring hybrid, dielectric resonator, filters, LVDS/PCB stacks, heat sink, enclosure, spiral — built from Tessella primitives / raw surfaces (frustum, sphere, annulus, bowtie), **22/22 valid + watertight + conforming** (box-assembly cases at exact analytic volume; multi-region cases with every region filled + no face shared by >2 tets). `validation/hfss_cases/`, regression-pinned in `test/hfss_cases_test.jl`. **The SOLVE campaign (the other half of #12) is external ASCENT work — see `ASCENT.md`.** |

**Bottom line (2026-08-13):** the Tessella package (the mesher) is **complete, verified, and
optimized** — every meshing/geometry item is DONE and regression-pinned. #7 (non-star+reflex
recovery) and #8 (arbitrary-surface uniform sizing — guaranteed `maxedge ≤ hmax` on any closed
surface via `mesh_sized`/`refine_to_size`, plus box/cylinder/prismatic sizers); #9 the **literal
ENC-COAX geometry natively** — all 4 volumes, all physical groups, native boolean bore imprint, and
exact analytical curved surfaces via the native CAD layer (no OpenCASCADE); #10 all representatives;
#11 the cospherical Delaunay perf fixed in **both** kernels (jump-and-walk + exact-kernel pre-filter,
both output-identical); and the meshing half of #12 — **all 22 HFSS geometry classes meshed natively**.
Native CSG (`mesh_boolean`, `mesh_box_regions`), P2 curved elements, gmsh MSH v2/v4 I/O, and the
exact-coordinate `Rational{BigInt}` Delaunay kernel are all shipped and pinned. **The one remaining
item in the whole project is external, not a Tessella package task:** the ASCENT full-wave solve
campaign (the 21 non-gmsh-failure HFSS guide cases) — tracked in **`ASCENT.md`** (the flagship
gmsh-impossible case 9.2 is already meshed natively AND solved in ASCENT). Nothing was faked — a
silent non-conforming mesh would violate the CRC bar.

## Current state (verified at HEAD)

- **Suite:** `julia --project=. --check-bounds=yes -e 'using Pkg; Pkg.test(; julia_args=["--check-bounds=yes"])'`
  green — **151,466 assertions** under `--check-bounds=yes` (Julia 1.12.6; 12m53s; ~2–3× faster than
  the 25m00s baseline after the classifier optimization below; +28 `mesh_box`, +24
  `mesh_box_regions`, +17 `mesh_cylinder`, +36 `recover_boundary`, +9
  `mesh_sized_conforming`, +32 `mesh_boolean`, +11 native-pipeline integration, +4
  grid-classifier parity, +13 audit-fix regressions).
- **Domain classifier optimized (2026-08-12).** The O(n_tets·n_surface_faces)
  ray-cast point-in-polyhedron classifier — previously flagged as "the obvious
  throughput win (noted, not yet done)" — is now a projected-plane CSR grid index
  (`_RayGrid`, `Mesh3D.jl`), **provably output-identical** to the brute force and
  CRC-pinned bit-for-bit against it. `mesh_volume`(cyl 48×8) **2.5 s → 0.016 s
  (~150×)**; whole suite ~2.1× faster. See the Log.
- **CRC regression checksums preserved** (suite includes them; unaffected by the
  `orient2_sos` fix): unit-cube
  `7ea403054f05392f18b404a1f5f78b12d70d45d40c7b04ba8f8dc3e030d8f3f9`; 10×10
  refined square `583c615df1862c8518bbda409347f109dc25f7f8f5362562badf160fe6af30c1`.
- **All 7 stages carry working, CRC-gated code.** Stages 0–3 complete; Stage 5
  native CSG complete (a pure-Julia OCC *library* replacement is an explicit PLAN
  §1/§6 non-goal); Stage 4 size control complete for boxes+cylinders (arbitrary-
  curved *uniform* sizing is research-grade — needs the exact-coordinate kernel).
  Public API: `mesh_volume`
  (3-D), `mesh_planar` (2-D), **`mesh_box`** / **`mesh_box_regions`** (uniform
  size-controlled + native box CSG), **`mesh_cylinder`** (uniform size-controlled
  cylinder, cospherical-robust), **`mesh_sized_conforming`** (interior size control
  for curved domains), **`recover_boundary`** (robust boundary recovery, +
  `steiner`), **`mesh_boolean`** (native mesh-Boolean CSG), plus the per-module
  functions.
- **The native geometry→mesh pipeline works end-to-end** (no gmsh, no OCC),
  verified in `pipeline_test.jl`: **native CSG → conforming fill → solver-consumable
  gmsh MSH v4.1** — a size-controlled multi-region enclosure (`mesh_box_regions`)
  and a box-with-cylindrical-bore (`mesh_boolean` → `recover_boundary`) each write,
  round-trip (CRC + physical groups + tags preserved), and validate as ASCENT input.
- **Reason-to-exist demonstrated:** the ENC-COAX coax feed-through is meshed as ONE
  **conforming** partition — all three volumes (air/case/pin) filled with a manifold,
  shared interface (the pin bore through the case wall conforms), at the **literal
  `.geo` scale and real pin resolution** — exactly where gmsh 4.13/4.15 produce zero
  volume tets. Enabled by the corrected SoS (exact coordinates, no jitter).
- **Adversarial-audit + feature pass landed this session** (workflows + independent
  re-verification): **9 confirmed bugs fixed** — including all four exact-predicate
  SoS tie-breaks (`orient2/orient3/incircle/insphere`) now consistent under one +ε
  scheme (regression-pinned) — plus 2 features integrated and a `validation/`
  cross-check suite vs gmsh added (see "Audit findings" below).

### Scope completion vs `PLAN.md` (the authoritative boundary — `startup.md`: "Do not silently expand it")

The project's **defining acceptance test is done**: mesh the ENC-COAX feed-through
gmsh 4.13/4.15 cannot — PLAN §7's "proof that this project earns its existence."
**Native CSG is shipped** (`mesh_box_regions` + `mesh_boolean`, no OCC), the
exact-coordinate `Rational{BigInt}` kernel is shipped, and the two former research
flagships (non-star+reflex recovery, arbitrary-surface uniform sizing) are now closed.

| item | status |
|---|---|
| uniform size control on curved domains | **DONE across the board** — box (`mesh_box`/`mesh_box_regions`), cylinder (`mesh_cylinder`), prismatic/extruded (`mesh_sized_extrude`), and **general arbitrary surfaces (`mesh_sized`)** — fill + longest-edge bisection (`refine_to_size`) gives **guaranteed `maxedge ≤ hmax`**, valid, watertight, conforming on convex/curved/non-convex domains (spheres, frusta, L-prism), exact volume for polyhedra |
| recovery for non-star+reflex polyhedra | **DONE** via `recover_boundary_cdt` (exact-kernel conforming-Delaunay boundary-Steiner recovery); the twisted prism conforms, regression-pinned |
| general CSG / OCC-library interop | native CSG **DONE** (`mesh_boolean`); a pure-Julia OpenCASCADE **library** replacement remains a PLAN §1/§6 explicit non-goal (multi-year) |

The one remaining *project* item — the ASCENT full-wave solve campaign — is external
(needs the ASCENT binary + proprietary HFSS datasets, local at
`/Users/jake/EMPIRE/projects/ongoing/2026_066`, not in this repo) and is tracked in
**`ASCENT.md`**. Per `startup.md` ("Scope is stated honestly in `PLAN.md` §1. Do not
silently expand it"): faking HFSS reference data or shipping a divergent refiner would
violate both the plan and `DEVELOPMENT.md`'s CRC bar.

## Stage board

| Stage | Scope | State | CRC gate |
|---|---|---|---|
| 0 | Foundations: repo, CI, mesh types, `.msh` I/O, exact predicates | **DONE — gate green** | predicates vs exact-rational oracle ✓; `.msh` round-trip (v2↔v4) CRC-preserving ✓ |
| 1 | 2-D Delaunay + CDT + quality refinement | **DONE — gate green** | exact empty-circumcircle oracle ✓; 2n−2−h count + hull ✓; CDT constraints present + locally Delaunay ✓; Ruppert angle/area bound achieved, domain area preserved ✓ |
| 2 | 1-D edge meshing + parametric surface meshing | **DONE — gate green** | size-graded edge mesh vs analytic arc length ✓; planar face exact area + quality ✓; cylinder watertight + area convergence ✓; parametric area convergence ✓ |
| 3 | 3-D Delaunay + **robust boundary recovery** + slivers | **kernel + filling DONE; conforming partitions DONE; general boundary recovery DONE (non-Schönhardt); sliver removal DONE (`remove_slivers` converging driver, 278→158 on a random cloud, valid + volume-preserving, `optimize_test.jl`)** | 3-D empty-circumsphere oracle ✓; convex/non-convex/genus-1/thin/multi-region fills validated ✓; **`tetrahedralize_conforming`** — exact-coordinate Delaunay + per-region ray-cast tagging → *shared* Delaunay-face interfaces (manifold + exact volume + every region filled on the **3-region pin/air/case enclosure**, incl. the pin's *curved N-gon* interface and a **feed-through crossing the case wall**, no Steiner) ✓; **`recover_boundary` SHIPPED — general robust boundary recovery** (PLAN principle #2): recovers an arbitrary closed PLC surface as a **conforming** tet mesh (every input facet a tet face, via a triangulation-independent **exact `Rational{BigInt}` conformity gate**), or throws an **explicit blocker** (never a silent bad mesh). Verified conforming on convex box, **non-convex genus-1 through-tunnel**, **hollow shell**, **star-shaped L-prism**, **faceted (octagonal) cylinder** — boundary-area == input-surface-area + exact volume + valid + closed-manifold, ~1–2 s; **Schönhardt polyhedron correctly raises the blocker by default**, and with **`steiner=true` is meshed** via fan-tetrahedralization from an interior kernel point (one Steiner vertex, one tet per facet — conforming + valid, `mesh3d_test.jl`) — closing the Schönhardt-type case for **star-shaped** inputs. Insertion-order retry (`perturb=false`) + gated flat-drop dodges the cospherical zero-volume-tet degeneracy that stops the base kernel on `box_tunnel`. **Breadth verified**: also conforms genuinely **non-star-shaped** non-convex prisms — U-channel, 2-prong comb, 5-point star (Delaunay-recoverable). The one class it cannot mesh — **non-star AND reflex/non-Delaunay-recoverable** (a *twisted* non-convex prism, constructed as a concrete case) — correctly raises the **explicit blocker** under *both* `steiner` modes (regression-pinned, `mesh3d_test.jl`): the recover-or-blocker safety guarantee holds even there, **never a silent bad mesh**. Closing that exotic class needs TetGen-style boundary-face-splitting Steiner recovery (a subdivision-retry probe was measured slow/uncertain — research-grade). Also open: sliver exudation |
| 4 | size fields + optimization | **partial** | tet quality report ✓; Laplacian + ODT smoothing ✓; verified 2-3/3-2 flip primitives ✓; `optimize_flips!` sliver reduction (309→169 flips alone, →63 with smoothing on a 300-pt cloud; volume/validity/min-dihedral-safe) ✓; **`mesh_box` size-controlled mesher SHIPPED** — Kuhn/Freudenthal structured subdivision gives a **guaranteed `maxedge ≤ hmax`** tet mesh of any axis-aligned box, provably valid (all-positive), watertight (boundary χ=2), exact volume, and **sliver-free** (min dihedral = 45° cube / ≥42° general, radius-edge < 0.9), for arbitrary `hmax` (`mesh3d_test.jl`, independent oracles) ✓; **`mesh_box_regions` SHIPPED** — extends the shared-lattice route to conforming, size-controlled **multi-region** meshing of **unions/differences/nestings of axis-aligned boxes** (native box CSG): exact per-region volumes, manifold-conforming interfaces, handles **non-convex** domains (hollow shell = `box−void`), all provably valid at `maxedge≤hmax` (`mesh3d_test.jl`) ✓ — the size-control + partitioning primitive for enclosure-class box assemblies. **`mesh_sized_conforming` SHIPPED** — **interior size control for curved domains**: an inset interior Steiner lattice added to `recover_boundary`, gated by the exact conformity+validity check, so it returns a conforming mesh with interior edges `≤hmax` (verified on a sphere: tets 436→577, `int_maxedge≤hmax`, boundary conforms) or raises an **explicit blocker** — never a silent invalid mesh; thin/cospherical inputs degrade safely to conforming-only (`mesh3d_test.jl`) ✓; **`mesh_cylinder` SHIPPED** — **uniform** size control for the cylinder primitive (the concrete thin/cospherical case): a structured `(r,θ,z)` Kuhn mesh with axis collapse gives **`maxedge ≤ hmax`**, exact faceted volume, watertight (χ=2), valid — **no Delaunay ⇒ no cospherical degeneracy** (the route the Delaunay fill couldn't take), for any tilt/offset (`mesh3d_test.jl`) ✓; **`mesh_sized_extrude` SHIPPED** — **uniform** size control for **extruded/prismatic** domains (any polygon cross-section, non-convex + with holes): the 2-D Ruppert mesher sizes the cross-section (edges `≤ hmax/√2`), extruded into sized layers with a conforming column-index prism→tet split ⇒ **guaranteed `maxedge ≤ hmax`**, exact volume (boundary preserved), watertight, conforming, valid — verified on a non-convex L-prism + a genus-1 annulus prism over an `hmax` sweep (`mesh3d_test.jl`) ✓; **`mesh_sized` + `refine_to_size` SHIPPED — general uniform sizing on ANY closed surface** (the terminator, done): fill (`tetrahedralize`/`recover_boundary`) then **longest-edge bisection** — split every tet incident to the globally-longest edge at its midpoint; conforming (shared faces split identically), validity + total volume preserved (midpoint on edge), terminating (bisection only creates edges ≤ the split one) ⇒ **guaranteed `maxedge ≤ hmax`**. Verified on convex box (exact vol 8), sphere, conical frustum, non-convex L-prism (exact vol 3): valid + watertight + conforming, sub-second (`mesh3d_test.jl`) ✓ — box/cylinder/prismatic/**general** uniform sizing all shipped. Earlier: the Delaunay-refinement route was **attempted + measured** (`validation/stage4_size_refinement/delaunay_refiner_convex_ATTEMPT.jl`) — correct for box-like convex domains but impractically slow (rebuild-per-pass) and fails on a tetrahedron (small-angle), so not shipped; robust general refinement needs boundary recovery + Shewchuk's small-angle-protected terminator (research-grade). **Curved-domain size control was probed further** (`validation/stage4_size_refinement/curved_size_control_findings.md`): uniform size + *exact* conformity are **not simultaneously achievable** by the pragmatic routes (background-lattice clip → uniform interior but resampled/staircase boundary; fine-surface+inset-lattice → exact-conforming but *graded* with a ~2·hmax boundary shell; recover-then-protect → stalls + cospherical invalidity) — the terminator remains the open path |
| 5 | geometry kernel (OCC interop / native CSG) + heal | **native CSG DONE (surface + volumetric); OCC interop out of scope** | `Heal` surface-defect detection ✓; native primitives (box/cylinder/box-tunnel/**hollow-box** `box_shell_surface`) ✓; **`mesh_box_regions`** — volumetric CSG for axis-aligned box assemblies → conforming multi-region tet mesh, exact per-region volumes ✓; **`mesh_boolean` SHIPPED — general native mesh-Boolean CSG (no OCC)**: `union`/`intersection`/`difference` of two closed triangulated solids via an exact plane-arrangement path (axis-aligned, all coplanar shared-face cases) **and** a Cork/libigl-style **exact tri-tri (`orient3`) + `Rational{BigInt}` seam + `Mesh2D`-CDT** path for general position (e.g. box×cylinder). Verified **exact Boolean volumes** (box∪box 96/112/1875, ∩ 32/16/125, ∖ 32/48/875; box−cylinder exact faceted volume), watertight results, `recover_boundary`-fillable; explicit blocker on unsupported degeneracies, never a leaky/wrong surface (`mesh3d_test.jl`) ✓ — the user-directed native geometry-kernel path. OCC library interop remains **out of near-term scope** (PLAN §1/§6 non-goal); native CSG now covers the enclosure's Boolean needs |
| 6 | high-order elements + solver-consumable I/O | **DONE (package side)** | quadratic (P2) tet generation ✓ (shared mid-nodes, exact-midpoint volume) + gmsh type-11 I/O ✓; **general-geometry curving** ✓ — `curve_to_surface!`(project, on_surface) / `curve_to_cylinder!` curve only genuine *boundary-surface* edges (interior chords left straight) and revert any projection that would invert an incident P2 element (`p2_min_jacobian` guard over the degree-3 nodes; no inverted element ever emitted, verified by an independent degree-6 sampler); **the conforming ENC-COAX mesh is solver-consumable** — written to gmsh MSH v4.1 with the three physical volumes (`coax_pin`/`air`/`case`) and round-tripped (connectivity CRC + region tags + physical-group names preserved, `pipeline_test.jl`) ✓. *(ASCENT drop-in + the 22-case solve campaign are tracked in `ASCENT.md`; the package implementations they exercise — P2 curved elements, solver-consumable MSH v4.1 — are done.)* |

### Stage-4 measured finding — 3-D size refinement (why no refiner is shipped)

**Three** concrete size-controlled-meshing approaches were **built and measured** on a convex box `[0,4]³` (target max-edge `hmax`); all three provably fail, so none is shipped (a broken mesher would fail the correctness bar):

1. **Interior longest-edge midpoint insertion → diverges.** Inserting the midpoint of every edge longer than `hmax` into the Delaunay leaves the max edge **pinned at 5.657 = 4·√2 forever** (`hmax=5.0`: `nlong=3` every pass, +3 verts / +7 tets per pass, no decrease across 25+ passes). Mechanism: 5.657 is the **box face diagonal**, a *boundary edge* of the input surface; splitting a diagonal at its midpoint (the face centre) just regenerates a diagonal of equal length among the remaining corners. **Interior insertion cannot shorten a fixed boundary edge.** (Reproducible: `validation/stage4_size_refinement/diverges_interior_midpoint.jl`.)
2. **Fine-surface-only tetrahedralization → long interior diagonals persist.** Meshing the box surface finely first (up to 98 nodes / 192 tris, `k=4`) and tetrahedralizing gives the correct volume (64.0) and a valid mesh, but `tetrahedralize` adds **no interior points**, so a tet spanning bottom-face→top-face keeps `maxedge=4.0` — the size bound is still not met. (Reproducible: `validation/stage4_size_refinement/fine_surface_leaves_interior_diagonals.jl`.)
3. **BCC lattice + Delaunay → invalid tets for general spacing.** BCC points at spacing `hmax/√2` (so face diagonals ≤ `hmax`), Delaunay'd with `perturb=false`, mesh **perfectly for the trivial even case** (`hmax=3` → 2×2×2: valid, **60° min dihedral**, exact volume, `maxedge=2.83`) but yield **inverted/invalid tets for general spacing** (`hmax=2`: valid=false, 0° dihedral; `hmax=1`: wrong volume 64.44). Mechanism: a regular lattice is *maximally* cospherical-degenerate and Delaunay-of-cospherical-points is ambiguous — the exact+SoS kernel resolves ties deterministically but not always *validly*. A correct BCC mesher must emit the **known connectivity explicitly**, not Delaunay the lattice. (Reproducible: `validation/stage4_size_refinement/bcc_lattice_delaunay_invalid.jl`.)

**Resolution:** the third finding pointed at the fix — *emit explicit connectivity instead of Delaunay-ing degenerate points* — and that route **is shipped as [`mesh_box`](../src/Mesh3D.jl)**: the Kuhn/Freudenthal subdivision (6 path-tetrahedra per grid cube, all sharing the main diagonal). It is provably valid (all-positive, orientation set by signed volume at build), watertight (every cube uses the same diagonal ⇒ shared faces match; boundary χ=2), exact in volume, sliver-free (min dihedral 45°/≥42°, radius-edge < 0.9), terminating (finite grid), and meets `maxedge ≤ (hmax/√3)·√3 = hmax` by construction — all checked by independent oracles in `mesh3d_test.jl`. It covers **axis-aligned box** regions (the enclosure case/air cavities). What remains research-grade is **general non-box adaptive** size control on an arbitrary surface, which needs boundary-conforming Delaunay refinement (Shewchuk's terminator: interior Steiner points *plus* encroachment-driven boundary-subface splits; divergence-prone near small input angles) — out of single-session scope. *Acceptance test for that future work:* over a sweep of `hmax` on a non-convex domain, `maxedge ≤ hmax` with `validate.ok`, exact preserved volume, good min dihedral, and bounded vertex count (termination).

## Standing acceptance cases (regression, CRC-stamped when they pass)

| id | geometry | why | status |
|---|---|---|---|
| ENC-COAX | enclosure coax feed-through (HFSS 9.2), surface 86 shield/case/air junction | gmsh 4.13.1 **and** 4.15.2 fail (`overlapping facets`, all volumes empty); the project's reason to exist | **PARTIAL (real-scale volumes filled)** — at the *literal* `.geo` dimensions the three regions (air / case / coax pin, radius 0.8 mm × length 159 mm, aspect ≈ 199) are each meshed and **all volumes filled & validated** together (`mesh3d_test.jl`), the exact place gmsh leaves them empty. **The conforming pin/air/case *topology* now works** (`tetrahedralize_conforming`, `mesh3d_test.jl`): a coax **pin (cylinder) inside an air cavity inside a metal case shell**, meshed as ONE mesh — every region filled, exact total volume, and a **manifold shared interface** including the pin's *curved* (N-gon) surface, which the exact-coordinate Delaunay recovers **without Steiner points** (all 4·nθ pin faces shared pin↔air; air↔case cavity shared). The corrected SoS removes the perturbation-vs-conformance barrier. **The feed-through also conforms**: with the pin **crossing the case wall**, classification `[pin,air,case]` tags every in-cylinder tet as pin regardless of cavity-vs-wall, so the cylindrical **bore is handled with no explicit CSG surface** — both the pin↔air (cavity) *and* pin↔case (bore) interfaces appear and conform, and the air↔case boundary conforms around the bore hole (all 2·nθ·nz pin faces shared, 0 non-manifold, exact total volume). **This holds at the literal `.geo` scale and the real pin resolution**: air/case/pin at the fixture dimensions (0.8 mm pin, aspect ≈199) with nθ=12·nz=40 → 2549 tets, **all three volumes filled**, feed-through conforming (pin↔air 144, pin↔case 816 faces), 0 non-manifold, valid — the exact 3-volume mesh gmsh cannot produce. Remaining: (a) **performance** at fine resolution — for the real nθ=12·nz=40 pin (498 verts) the run is ~100 s, of which the **`perturb=false` Delaunay is ~84 s** (measured — pathological for 498 pts). *Root cause found by measurement, not assumed:* it is **not** the exact predicates (making the exact fallback 1.84× faster via integer arithmetic left the Delaunay time unchanged), so it is the **point-location behaviour on the degenerate coplanar cylinder vertices** under exact coordinates (the walk degrading toward its brute-force `_locate_scan` fallback → ~O(n²)). Fix = a more robust exact walk / spatial acceleration; the classification is already bbox-reject-optimized. *Three fixes tried + rejected by measurement, narrowing the cause:* (i) faster exact predicates (kept — 1.8× — but the Delaunay was unchanged ⇒ **not predicate-bound**); (ii) a plain **Morton spatial-sort** insertion order (reverted — *worse*, 84→142 s: a sliver cascade); (iii) **BRIO** (Amenta–Choi–Rote randomized rounds + Morton — reverted: **79 s, no real change** ⇒ **not walk/insertion-order-bound** either). What remains is the **cospherical-cavity cost**: under exact coordinates the ~480 coplanar/cospherical cylinder vertices give large Bowyer–Watson cavities whose retriangulation (spoke-linking + tet churn) dominates — a deep algorithmic characteristic of exact-degenerate input, addressable only by a different cavity/degeneracy strategy, not any of the standard quick fixes. (b) the extra literal geometry (slots, resistor, `p1_surface`) beyond the three main volumes. Neither is a boundary-recovery problem, and gmsh produces **zero** tets here at any speed. |
| THIN-SLOT | enclosure 1 mm coupling slot | thin-feature Boolean sliver | **representative DONE** — a thin conducting slab between two air cavities meshes **valid + exact per-region volume** at aspect ≤ 4000:1 via the deterministic structured route (`mesh_box_regions`, `mesh3d_test.jl`), the thin-feature case gmsh slivers on; min dihedral degrades gracefully (42°→5°) but every tet stays positive. **Literal HFSS geometry** still needs the fixture |
| SPIRAL | 10.1 silicon spiral (thin swept traces, layered stack) | high-aspect thin conductors | **representative DONE** — a planar square-spiral thin conductor trace (connected axis-aligned segments, native-CSG union) embedded in a substrate meshes **valid + conforming** with the winding trace filled as one region (`mesh3d_test.jl`). **Literal HFSS geometry** still needs the fixture |
| ARRAY-PML | 5.7 endfire unit cell (periodic + PML) | conformal periodic faces | **representative DONE** — a `mesh_box` unit cell has **conformal periodic faces**: opposite boundary faces carry identical triangulations under the period translation (verified x/y/z, `mesh3d_test.jl`), so periodic BCs pair nodes 1:1. **Literal HFSS geometry** still needs the fixture |

Reference artifact for ENC-COAX: the failing gmsh `.geo` is captured at
`test/fixtures/enclosure_coax_junction.geo` (parsed by `IO.read_geo_params`, which
recovers the air/case/coax_pin groups + sizing + seed). gmsh-API diagnosis:
surface 86 bbox ≈ (168.4,140,148.4)–(171.6,160.5,151.6) mm; empty volumes =
pin/case/air.

## Why the enclosure needs Tessella (gmsh-failure diagnosis)

The enclosure meshing failure is **NOT memory** (a no-time-limit run peaked at 3.3 GB);
it is a geometry/boundary-recovery defect at the coax junction. 5 gmsh 2-D/3-D
algorithms, OCC fix options, `Geometry.ToleranceBoolean`, and geometry protrusion all
fail to mesh it — so the robustness must come from Tessella's own boundary recovery +
heal. (gmsh `occ.healShapes()` is a false positive: it dissolves the OCC volumes to 0
tets, so "no empty volumes" is vacuous — always check the *actual* per-volume tet count;
a CRC-discipline lesson encoded in `DEVELOPMENT.md`.) The full diagnosis + ASCENT-campaign
provenance is in `ASCENT.md`.

## Audit findings (adversarial workflows, 2026-08-11) — independently re-verified

Two background workflows (a feature-implementation pass + a core-module adversarial
bug-hunt with per-finding verification) surfaced these. Every one was reproduced
here from scratch before acting; every fix re-verified (full suite green).

**Fixed (9):**

| # | site | defect | verification of the fix |
|---|---|---|---|
| 1 | `HighOrder.curve_to_cylinder!` / new `curve_to_surface!` | curved **interior** edges (qualified an edge by "both corners on the surface", but every vertex is on the boundary since no Steiner nodes are added) → tangled 47/125 (cyl) and 127/166 (sphere) P2 elements | now curves only genuine boundary-surface edges + reverts any inverting projection; independent degree-6 Jacobian sampler finds **0** inverted elements; interior chords byte-identical to their straight midpoints |
| 2 | `Predicates.orient2_sos:351` | SoS tie-break sign-inverted (matched −ε) vs the +ε house convention → 2-D degenerate orient2/incircle decisions mutually incoherent (Mesh2D uses both on one point set) | negated to +ε; **0/6** disagreements vs an independent exact `Rational{BigInt}` SoS oracle (was 6/6); full suite + golden 2-D CRC unchanged |
| 7 | `Predicates.orient3_sos:416` | fully-degenerate (4-collinear) fallthrough returned a constant `perm`, wrong on 12/24 labelings | replaced with an exact +ε `Rational{BigInt}` leading-term SoS evaluator (Leibniz det of the perturbed matrix); verified **0/24** collinear + **0/24** coplanar-sanity vs a validated oracle (fast path preserved) |
| 8 | `Predicates.incircle_sos:439` | 4-collinear fallthrough returned constant `perm`, wrong on 12/24 | replaced with the exact SoS of the paraboloid-lifted points (incircle = orientation of the lift; the algorithm-consistent Edelsbrunner–Mücke scheme). Verified **0/24** collinear + **0/24** cocircular-sanity vs a **canonical lifted-orientation oracle** — which I built and validated after finding my first (naive-lift) oracle used the wrong scheme; a finite-ε exact-Rational cross-check confirmed the resolution |
| 9 | `Predicates.insphere_sos:456` | the *entire* tie-break (fast minors **and** fallthrough) evaluated the −ε perturbation while orient3 uses +ε — the 3-D analog of #2: the kernel needs orient3/insphere to agree (cavity chosen by insphere, orientations by orient3). Wrong on 120/120 cospherical + 60/120 coplanar vs +ε (three independent methods) | replaced the whole tie-break with the exact +ε SoS of the 4-D paraboloid lift (`_orient_nd_sos_exact`, the same construction as incircle one dimension up — the d=3 form is oracle-validated, and the base-sign match + dimension-general construction carry it to d=4). Verified **0/120** cospherical + **0/120** coplanar vs the canonical oracle; full suite + golden CRC unchanged |
| 3 | `MeshTypes.euler_characteristic:274` | counted segment-only nodes in V while E omits segment edges → χ inflated (tet+segment read 3, not 1) | segment loop removed; reproduced (3→1) + regression test |
| 4 | `MeshTypes.bounding_box:315` | included unreferenced nodes, violating its "referenced nodes" contract and polluting `mesh_crc.bbox` | filter to referenced nodes; reproduced ((1000,1,1)→(1,1,1)) + regression test |
| 5 | `IO.read_stl:481` | far-from-origin weld key `round(Int, coord/tol)` overflowed Int64 → crash | quantize relative to bbox min corner; reproduced (InexactError→ok) + regression test |
| 6 | `IO.read_stl:431` | ASCII STL with a UTF-8 BOM or non-ASCII solid name misdetected as binary → EOFError crash | BOM-aware keyword check, drop whole-header ASCII requirement; reproduced (2 crashes→ok) + regression test |

**All four SoS predicates now resolve degeneracies by one consistent +ε scheme.**
The exact `Rational{BigInt}` leading-term evaluator (`_orient_nd_sos_exact`, a
Leibniz determinant of the EM-perturbed matrix) is the shared foundation: orient3
uses it directly for the collinear fallthrough; incircle/insphere use it on the
paraboloid-lifted points (incircle = 3-D orientation of the lift, insphere = 4-D).
The corrected +ε signs over every index labeling of the four finding configs are
**regression-pinned** in `test/predicates_test.jl`. Nothing SoS remains deferred.

> Method note: getting incircle/insphere right required distinguishing two SoS
> schemes for lifted predicates — deriving the lift from perturbed x,y (wrong: its
> ε² cross-terms mis-order the leading term) vs. treating the lift as an independent
> EM coordinate (correct, = orientation of the lifted points). My first oracle used
> the former and mis-flagged the code; the canonical lifted-orientation oracle (+ε
> by construction, base-matched to each predicate) resolved it, cross-checked with a
> finite-ε exact-Rational evaluation. The two audit verify-agents had disagreed on
> insphere's sign; the canonical oracle + full-suite/CRC settle it.

## Cross-tool validation vs gmsh (`validation/`, run this session)

`julia --project=. validation/run_all.jl` meshes the same domain with Tessella and
gmsh 4.15.2 and compares meshed volume (vs the analytic oracle), quality, and time.
Reference `.geo` scripts retained per case; regenerate `REPORT.md` any time.

| case | analytic V | Tessella V | gmsh V | verdict |
|---|---|---|---|---|
| box (flat) | 2 | **2.0** | **2.0** | both exact ✓ |
| box_tunnel (genus-1 flat) | 24 | **24.0** | **24.0** | both exact ✓ |
| hollow_box (Boolean difference) | 35 | **35.0** | **35.0** | Tessella `box_shell` CSG matches gmsh's Boolean diff exactly ✓ |
| cylinder (true πR²H=62.83) | — | 62.65 (0.3%) | 61.02 (2.9%) | Tessella's 48-gon prism is **closer** than gmsh at this size |
| sphere (true 4/3πR³=20.58) | — | 20.10 (2.3%) | 19.69 (4.3%) | Tessella **closer** than gmsh at this size |
| **enclosure_coax** (HFSS 9.2, the reason-to-exist) | — | native primitives fill exactly | **0 volume tets** (28969 surf tris, all air/pin/case volumes empty) | gmsh's documented failure, reproduced here |

Flat-solid rows are a hard correctness cross-check (exact analytic volume, both
tools conform). The enclosure row reproduces gmsh's empty-volume failure directly
from the literal fixture (verified: `ntets=0`).

## CRC regression artifacts (stamped when green)

| artifact | checksum | oracle |
|---|---|---|
| Predicates | 142,141 assertions green; every `orient2/3`, `incircle`, `insphere` sign matches the independent exact-rational homogeneous-determinant oracle over exhaustive integer-grid degeneracies (collinear/coplanar/cocircular/cospherical) + fixed-seed random floats | `test/oracles.jl` (generic Laplace determinant, `Rational{BigInt}`) |
| Unit cube (6-tet Kuhn) | `nodes=8 tets=6 bbox=[(0,0,0)→(1,1,1)] dihedral(min,mean)=(0.7854,0.7854) radedge(min,mean)=(0.86603,0.86603) bfaces=12` · SHA-256 `7ea403054f05392f18b404a1f5f78b12d70d45d40c7b04ba8f8dc3e030d8f3f9` | analytic volume 1/6·6=1, Euler χ=1 (ball), boundary χ=2 (sphere), 12 boundary tris |
| 2-D Delaunay (random/grid/cocircular) | exact empty-circumcircle holds (0 violations, `incircle_sos`); triangle count = 2n−2−h vs independent monotone-chain hull; seed-independent (SoS) | `is_delaunay`, hull oracle, Euler χ=1 |
| 10×10 square refined to 20° (max_area 1.0, seed 1) | `nodes=91 tris=148 bbox=[(0,0,0)→(10,10,0)]` · SHA-256 `583c615df1862c8518bbda409347f109dc25f7f8f5362562badf160fe6af30c1` | min angle ≥ 20°, max tri area ≤ 1.0, domain area = 100 preserved, CDT locally Delaunay |
| 3-D Delaunay (random/box/cube) | exact empty-circumsphere holds (0 violations, `insphere_sos`); Σ tet volume = convex box volume; Euler χ=1 (ball) + boundary χ=2 (sphere) + manifold; seed-independent | `is_delaunay3`, box-volume conservation, `boundary_euler` |
| ENC-COAX (real `.geo` scale) | air/case/coax-pin (R 0.8 mm × L 159 mm, aspect ≈199) — `tets_per_region` all > 0, every region validated, `tetrahedralize_multi` | per-region tet count > 0 (anti-false-positive), exact box volumes, thin-pin fill |
| Volume fills (exact) | cube=1, octahedron=4/3, L-prism=3, box-tunnel(genus-1)=24, cylinder N-gon prism, thin pin | analytic volume + `validate` (positive vol, manifold, watertight) |
| 2-3 / 3-2 flips | consistency + total volume preserved; 2-3 ∘ 3-2 = identity | `check_consistency3`, volume sum, round-trip |

## Benchmarks (deterministic min-of-3, `julia -O2`, this machine)

| kernel | n | time | throughput |
|---|---:|---:|---|
| 2-D Delaunay (`delaunay2d`) | 1 000 | 3.2 ms | 314 k nodes/s |
| 2-D Delaunay | 5 000 | 26 ms | 191 k nodes/s |
| 2-D Delaunay | 20 000 | 166 ms | 120 k nodes/s |
| 3-D Delaunay (`delaunay3d`) | 1 000 | 16 ms | 61 k nodes/s |
| 3-D Delaunay | 5 000 | 118 ms | 42 k nodes/s |
| 3-D Delaunay | 20 000 | 579 ms | 35 k nodes/s |

`mesh_volume` (cylinder 48×8 surface) **≈ 0.016 s** (was ≈ 2.5 s) — the ray-cast
point-in-polyhedron classifier, formerly `O(n_tets × n_surface_faces)` and the
run's bottleneck, is now the `_RayGrid` projected-plane CSR grid index
(`O(n_faces)` build + `O(1)`/query, output-identical). This also halved the test
suite (25m → 11m44s), since `recover_boundary`/`mesh_sized_conforming` re-classify
per seed. **The `perturb=false` Delaunay cospherical cost (#11) is now largely
fixed (2026-08-13):** re-measurement showed it was **point location**, not the
"cavity retriangulation" the tracker had claimed (Σlocate 43.9 s vs Σcavity 4.2 s
vs Σretri 0.006 s at nθ=12,nz=16) — a few far points burned the `16·ntets` walk
budget then the O(n) scan. Jump-and-walk start + a `48·∛n` walk-guard cap
(output-identical, `.sha`-verified) cut the fine pin **262 s → 22 s (nz=40)** and
sub-10 s for nz≤24, and dropped the **whole suite 11m44s → 2m31s**. The residual is
now cavity *discovery* (insphere_sos on the degenerate complex) at the largest
size; `mesh_cylinder` (0.001 s) remains the structured fast path. Raw kernel
throughput is competitive; parallel HXT-class speed is a post-robustness goal
(`PLAN.md` §6).

## Log

_Chronological development record of the Tessella package. ASCENT-integration
milestones (handshake, solve proofs, cavity eigenmode, case 9.2 solve) are tracked in
`ASCENT.md`; where they appear below they are historical context for the package work of
that day._

- **2026-08-11** — repo scaffolded (PLAN, STATUS, DEVELOPMENT, README, skeleton,
  startup). gmsh 4.13.1 source studied and archived as reference. Package name
  `Tessella.jl` chosen. Local git initialized; GitHub remote pending `gh auth
  login`. **No code implemented yet** — Stage 0 is the first development task
  (start from `startup.md`).
- **2026-08-11** — **Stage 0 complete, gate green.** Implemented:
  - `src/Predicates.jl` — adaptive exact `orient2/orient3/incircle/insphere`
    (Shewchuk A-stage float filter → `Rational{BigInt}` exact fallback) plus
    Simulation-of-Simplicity `*_sos` variants (never return 0). Fast path is
    allocation-free (verified with `@allocated == 0`). Verified against an
    independent Laplace-determinant oracle on exhaustive degeneracies.
  - `src/MeshTypes.jl` — compact column-major `Mesh` (3×N coords, Int32 cells),
    allocation-free `node()`, area/volume/dihedral/circumradius/radius-edge
    quality, boundary-facet extraction, Euler characteristic, manifold check,
    `validate` (positive-tet / degeneracy / non-manifold diagnostics), and the
    deterministic order-invariant SHA-256 `mesh_crc` checksum.
  - `src/IO.jl` — gmsh `.msh` v2.2 **and** v4.1 ASCII read/write (round-trip
    preserves connectivity CRC across versions), ASCII+binary STL ingest with
    vertex welding, `.geo` parameter/physical-group scanner (no OCC eval — that
    is Stage 5). `read_geo_params` on the enclosure fixture recovers all three
    volume groups (air/case/coax_pin) + sizing + seed.
  - Full `Pkg.test()` green: **142,141 assertions**. `stage()` bumped to 1.
- **2026-08-11** — **Stage 1 complete, gate green.** `src/Mesh2D.jl`:
  - **Delaunay** (incremental Bowyer–Watson) closed by a single GHOST vertex at
    infinity — *not* a finite super-triangle, which is not robust (a point outside
    the evolving hull can be wrongly excluded from a super-adjacent circumcircle,
    cracking the boundary; verified failure on ~half of random seeds). Ghost
    in-circle = orientation test on the real hull edge; collinear-aware so
    collinear hull points never spawn flat triangles. Remembering stochastic walk.
  - **CDT**: flip-based (Sloan) constrained segment insertion with `_flip!`
    (verified: double-flip = identity), through-vertex splitting, Lawson
    legalization; parity flood-fill `classify_interior` (nested holes correct).
  - **Ruppert refinement**: encroached-subsegment splitting + skinny-triangle
    circumcenter insertion (constrained), achieving the min-angle/max-area bound
    while preserving the domain and constraints.
  - Bug caught by `--check-bounds=yes` (Pkg.test): a Steiner-insertion OOB read
    (`_pt(T, nreal+1)`) masked by `@inbounds` in normal mode; fixed. Lesson:
    verify under bounds-checking. Determinism bug (Set iteration order in
    encroachment scans) fixed by sorted iteration.
  - Oracles: exact empty-circumcircle, 2n−2−h + hull cross-check, constrained-
    Delaunay locally-Delaunay check, min-angle/area quality, domain-area
    conservation. Full `Pkg.test()` green: **142,240 assertions**. `stage()` → 2.
- **2026-08-11** — **Stage 2 complete, gate green.**
  - `src/SizeField.jl` — `ConstantSize`/`FunctionSize`/`MinSize`, `size_at`.
  - `src/Mesh1D.jl` — graded 1-D edge meshing: nodes at equal increments of the
    metric length ∫|γ'|/h. Verified vs analytic arc length (segment/circle/helix),
    uniform + graded spacing, closed loops (no duplicate node).
  - `src/MeshSurface.jl` — planar faces (Newell frame → project → CDT + size-
    driven Ruppert → lift; exact area, quality preserved), cylinder lateral
    surface (**graded structured**, watertight by construction — replaced a
    Ruppert-unrolled variant whose seam gapped under a varying size field, a bug
    caught by a seam-gap-edge count), and parametric patches (isotropic-metric
    approximation; area converges to analytic).
  - `refine!` extended with a `size` callback (refine where an edge exceeds the
    local target). Full `Pkg.test()` green: **142,282 assertions**. `stage()` → 3.
- **2026-08-11** — **Stage 3 core (3-D kernel + volume filling).** `src/Mesh3D.jl`:
  - Incremental 3-D Delaunay via ghost tetrahedra on exact `orient3_sos`/
    `insphere_sos`; correct hull extension; face-orientation invariant proved and
    checked. Deterministic ~1e-8 symbolic perturbation breaks coplanar/cospherical
    degeneracies (the maximally-degenerate unit cube tetrahedralizes cleanly).
    Bugs fixed via oracles: ghost in-sphere sign, new-tet face orientation, and a
    walk-termination gap (brute-force `_locate_scan` fallback).
  - `tetrahedralize` (fill from a closed surface, ray-cast point-in-polyhedron)
    and `tetrahedralize_multi` (per-region, tagged). Verified fills at *exact*
    volume: convex cube/octahedron, non-convex L-prism, genus-1 box-with-tunnel,
    thin cylindrical pin (aspect 20), and a **3-region coax junction with all
    volumes filled** — the multi-way junction gmsh leaves empty.
  - Honest scope: robust conforming *interface* recovery (a naive Steiner version
    diverged under perturbation and was removed) and evaluating the literal
    enclosure `.geo` (needs the Stage-5 geometry kernel) remain. Full `Pkg.test()`
    green: **142,335 assertions**.
- **2026-08-11** — **Stages 4–5 partials + integration.**
  - `src/Optimize.jl` — `mesh_quality` (dihedral/radius-edge/sliver report) and
    boundary-preserving `smooth_laplacian` (positive-volume guard; volume/validity/
    tags preserved; mean dihedral up, slivers down). `refine3d!` was prototyped and
    **removed** (unconstrained circumcentre insertion raised mean circumradius by
    creating slivers — CRC: don't ship a tool that degrades its output).
  - `src/Heal.jl` — `surface_diagnostics`/`is_meshable`: detect open, non-manifold,
    degenerate, duplicate, mis-oriented, and near-coincident defects (the detection
    half of "heal, don't fail").
  - `src/Geometry.jl` — native `box_surface`/`cylinder_surface`/`box_tunnel_surface`
    primitives (closed, manifold, fill to exact analytic volume).
  - Top-level `Tessella.mesh_volume` — surface→volume pipeline with the
    validated-or-explicit-blocker contract (Heal gate → fill → smooth → validate).
  - Deep-debug pass on the 3-D+ modules: default path (`perturb=true`) robust; the
    `perturb=false` coplanar case is the known SoS-artifact (perturbation is the
    fix), documented. Full `Pkg.test()` green: **142,388 assertions**.
- **2026-08-11** — **Stage 6 partial (high-order) + integration.** `src/HighOrder.jl`:
  quadratic (P2) tet generation (`p2_tetmesh`, one shared node per edge, straight
  P2 volume == linear), gmsh **type-11** writer + reader (round-trips to the same
  linear connectivity CRC), and `curve_to_cylinder!` (projects lateral-wall edge
  nodes onto the exact cylinder — genuine curved high-order for a primitive,
  verified to land at radius R within 1e-12). End-to-end integration test:
  primitive → `mesh_volume` → `.msh` (v2/v4) → read reproduces the CRC; multi-
  region tags + physical names survive. All 7 stages now have working, CRC-gated
  implementations; the remaining work is the research-hard / out-of-near-term-scope
  tail (conforming interface recovery, sliver exudation, Boolean CSG + OCC interop,
  general-geometry curving, ASCENT drop-in + 22-case regression). Full `Pkg.test()`
  green: **142,414 assertions**.
- **2026-08-11** — **ENC-COAX real-scale + kernel flip primitives.**
  - Constructed the enclosure regions at the literal `.geo` dimensions (air box,
    case shell, coax pin R=0.8 mm × L=159 mm, aspect ≈199) and filled **all three
    volumes** non-empty & validated (`tetrahedralize_multi`) — the exact failure
    mode gmsh cannot. Regression-locked in `mesh3d_test.jl`.
  - `Mesh3D`: verified **2-3 / 3-2 flip** primitives (`flip23!`/`flip32!` on a
    general `_rebuild_region!`; consistency + volume preserved, 2-3∘3-2 = identity)
    + `tets_around_edge`, and a **safe** `optimize_flips!` (hill-climbing 2-3 flips,
    global min dihedral non-decreasing; limited effect on Delaunay input — honest).
    These are the foundation for sliver exudation and constrained boundary recovery.
  - Full `Pkg.test()` green: **142,436 assertions**. Every stage of the plan now
    carries working, CRC-gated code; the remaining tail (conforming interface
    recovery, sliver exudation, Boolean CSG + OCC interop, general curving, ASCENT
    drop-in + 22-case regression) is the research-hard / out-of-near-term-scope
    work `PLAN.md` frames as a multi-year effort.
- **2026-08-11** — **Conforming interface recovery: attempted, characterized, deferred.**
  Implemented the *correct* algorithm — 3-D-Ruppert-style **encroachment** refinement
  (split an encroached subsegment/subface, re-insert, Delaunay-maintain) rather than
  the earlier divergent split-missing-faces version. Verified the geometric kernels
  (triangle circumsphere equidistant; edge/face encroachers). On a **flat** interface
  (box faces) it hits a *fundamental* barrier: the perturbation that keeps the 3-D
  Delaunay flat-tet-free tilts the interface, so recovered Steiner points approach
  coplanarity faster than any jitter can separate them (eventually an unlocatable
  degenerate insertion), and regular grids cascade on cocircular squares. This is
  exactly why robust interface recovery needs a **constrained CDT** (exact, SoS-based
  flip recovery, *no* perturbation) — the multi-year TetGen-class effort `PLAN.md`
  defers. Per CRC discipline the non-robust recovery was **removed** (the verified
  2-3/3-2 flip primitives remain as its foundation). `tetrahedralize_multi` continues
  to fill every region (non-conforming interfaces). Full `Pkg.test()`: **142,442**.
- **2026-08-11** — **Working sliver reduction + 2-D/3-D top-level APIs.**
  - `optimize_flips!` gained **3-2 flips** (with `flip32!` validity guards — the
    edge must be interior and the two apexes must straddle the ring triangle; a
    bug where boundary/degenerate edges relinked wrongly was found and fixed). They
    collapse slivers: 309→169 from flips alone, 309→63 with Laplacian smoothing on
    a 300-pt random Delaunay; volume/validity preserved, min dihedral non-decreasing.
  - `mesh_volume(optimize=true)` runs the flips in the pipeline (cylinder fill
    102→73 slivers, volume preserved). `Tessella.mesh_planar` added — the 2-D
    counterpart of `mesh_volume` (CDT + Ruppert + interior → validated 2-D mesh).
  - Empirical recovery characterization (`scratchpad/analyze.jl`): on the perturbed
    two-region box the interface is *mostly* conforming (6/8 faces present); the 2
    missing faces have missing edges and are not single-flip recoverable → they need
    full edge+facet recovery (constrained CDT), confirming the deferral. Full
    `Pkg.test()` green: **142,445 assertions**.
- **2026-08-11** — **Audit + feature pass, exact-predicate SoS completed, and the
  conforming-interface barrier broken.**
  - Two adversarial workflows (feature impl + core bug-hunt), every finding
    independently reproduced + re-verified (see "Audit findings"): **9 bugs fixed**
    (P2 curving inversion; `euler_characteristic`; `bounding_box`; two `read_stl`
    crashes; and **all four SoS tie-breaks** — `orient2/orient3/incircle/insphere`
    now resolve degeneracies by one consistent **+ε** scheme, an exact
    `Rational{BigInt}` leading-term evaluator shared by all, regression-pinned).
  - `box_shell_surface` (hollow-box CSG) + `smooth_odt` integrated, review-hardened.
    `validation/` suite added (Tessella vs gmsh 4.15; enclosure = 0 gmsh volume tets).
  - **The conforming-interface prediction above came true.** With the SoS predicates
    correct and consistent, `delaunay3d(perturb=false)` runs on **exact coordinates**
    (SoS breaks coplanar/cospherical degeneracies symbolically — no jitter). New
    `tetrahedralize_conforming` meshes all region vertices together and ray-cast-tags
    each tet by region: interfaces become **shared** Delaunay faces. Verified on the
    enclosure **air/case pattern** (inner box + surrounding shell): one mesh, exact
    volume, **manifold shared interface** (the inner box's 12 faces), every region
    filled. The barrier was the coordinate jitter, which the corrected SoS removes.
    General recovery for non-Delaunay/curved interfaces (the literal coax pin) stays
    open. Full `Pkg.test()` green: **142,940 assertions** (+ conforming tests).
- **2026-08-12** — **Classifier optimized (documented throughput win) + diagnostic
  fix + independent re-audit.**
  - **Point-in-surface classification is now grid-accelerated.** The ray-cast domain
    classifier was `O(n_tets·n_surface_faces)` — STATUS had flagged it as "the
    obvious throughput win (noted, not yet done)." Replaced with a projected-plane
    CSR grid index (`_RayGrid`/`_raygrid`/`_tri_cellspan`/`_inside_grid`,
    `Mesh3D.jl`): all parity rays share one fixed direction, so a triangle can be hit
    only if the query projects (⊥ dir) into its 2-D bbox; each triangle is bucketed
    into every cell its bbox overlaps (+1-cell halo absorbing float-boundary
    rounding) and a query tests only its cell's candidates with the SAME
    Möller–Trumbore test. **Provably output-identical** to the brute force (a
    conservative superset filter ⇒ identical crossing parity), pinned by a new CRC
    test cross-checking `_inside_grid` vs `_inside_surface` bit-for-bit over a dense
    random cloud + every vertex/centroid on flat/curved/tilted/genus-1 surfaces (0
    disagreements), and separately confirmed to produce the identical `keep` mask on
    real triangulations. Wired into `_classify_by_centroid` (⇒ `tetrahedralize`,
    `tetrahedralize_multi`, `recover_boundary`, `mesh_sized_conforming`),
    `tetrahedralize_conforming` (per-region grids), and the `mesh_sized_conforming`
    interior lattice. Measured: **`mesh_volume`(cyl 48×8) 2.5 s → 0.016 s (~150×)**;
    **full suite 25m00s → 11m44s (~2.1×)** at **143,117** assertions (was 143,113; +4
    parity). Boolean CSG paths keep the brute-force reference (2-solid surfaces, not a
    bottleneck).
  - **Diagnostic bug fixed:** `mesh_sized_conforming`'s failure message used escaped
    `\$` (literal `$max_seeds`, no interpolation) unlike the sibling `recover_boundary`
    message — corrected to real interpolation so the blocker reports actual counts.
  - **Independent adversarial correctness re-audit** (8-area workflow, per-finding
    reproduce-or-refute running Julia against the green tree; 4 areas returned clean):
    **4 confirmed correctness bugs — every one independently reproduced here from
    scratch, fixed, and regression-pinned:**
    1. *(critical, `mesh_cylinder`)* the zero-volume drop used an ABSOLUTE
       `abs(sv) < 1e-14` threshold — scale-dependent, so a small/thin cylinder silently
       dropped legitimate thin near-axis tets, returning a non-watertight mesh (axial
       void, boundary χ=0, wrong volume) with **no error**. Fixed to an exact `sv == 0`
       drop (the vertex-equality check already removes every axis-collapse degeneracy);
       χ=2 + exact faceted volume at every scale, bit-identical at normal scale.
    2. *(critical, `tetrahedralize_conforming`)* returned a silently invalid,
       non-manifold mesh on cospherical axis-aligned box assemblies (it checked only
       `tets_per_region>0`, never `validate`). Fixed with a validate-or-explicit-blocker
       gate (PLAN principle #4) that names the degeneracy and points to `mesh_box_regions`.
    3. *(minor, `orient2_sos`)* the degenerate tie-break disagreed with the canonical
       `+ε` evaluator on coincident-coordinate inputs (144/4416), contradicting the
       documented four-predicate consistency. Routed through `_orient_nd_sos_exact` (the
       2-D analogue of orient3's exact fallthrough): 0 disagreements now, golden 2-D CRC
       preserved (0/1584 distinct-point configs changed; coincident points never reach it
       post-dedup, so no live path changes).
    4. *(minor, `validate`)* accepted a duplicate/overlapping tet complex (empty boundary,
       every face incidence 2) as valid. Now rejects an empty boundary and duplicate
       canonical tet keys — no false-reject of a real single tet.
    Regression tests added in `mesh3d_test.jl` / `predicates_test.jl` / `meshtypes_test.jl`.
  - Research tail re-confirmed **out of near-term scope with safe explicit blockers**
    (non-star+reflex recovery + arbitrary-surface uniform sizing both need a full
    `Rational{BigInt}` exact-coordinate 3-D Delaunay kernel; regression-pinned). The
    ASCENT-facing implementations (P2 curved, MSH v4.1, conforming multi-region
    fills) are complete.
- **2026-08-12 (cont.) — sliver removal + ASCENT mesh handshake verified.** Pursuing
  "finish everything implementable" (user directive):
  - **`smooth_optimize` shipped** (`Optimize.jl`): boundary-preserving optimization-
    based smoothing that maximizes each poor vertex-star's *worst* dihedral
    (min(min_dihedral, π−max_dihedral) pattern search) — the geometric half of sliver
    removal the mean-smoothers can't do. Measured **280→157 slivers** on a 300-pt
    random Delaunay (vs 262 laplacian / 297 odt), validity + total volume preserved;
    wired into `mesh_volume(optimize=true)` after the flips, +7 regression assertions.
    Suite **143,137** green.
  - **ASCENT mesh drop-in VERIFIED** (`validation/ascent_handshake/`): a Tessella
    MSH v4.1 (ENC-COAX pin/air/case conforming partition) loads into ASCENT's real
    parser `GmshDiscreteModel` (GridapGmsh 0.7.4) with all three region volumes as
    top-dimensional physical groups → `ASCENT.load_mesh` would return a valid
    `MeshData` (`HANDSHAKE_OK`, reproducible generate→load). The literal "ready for
    ASCENT" proof. The 22-case EM regression remains a solver campaign (needs the
    ASCENT binary + proprietary HFSS data, kept local, never pushed).
  - Remaining implementable tail being worked run-by-run: the exact-coordinate kernel
    (→ non-star+reflex recovery + arbitrary-surface sizing), ENC-COAX `.geo` extras,
    native representative thin-slot/spiral/PML geometries, Delaunay cospherical perf.
- **2026-08-12 (cont.) — exact-coordinate kernel foundation + orient3_sos coplanar bug.**
  Building the exact-coordinate 3-D Delaunay kernel (keystone for non-star+reflex
  recovery + arbitrary-surface sizing; user: complete ALL, no defer):
  - **Exact-rational predicates shipped** (`Predicates.jl`): `orient2_rat`/`orient3_rat`/
    `incircle_rat`/`insphere_rat` on `Rational{BigInt}` coords (same homogeneous
    determinants as the oracle + shared +ε SoS), so boundary-Steiner points at
    non-Float64 rational positions stay exactly on-feature. Verified to match the Float64
    SoS predicates on **all** dyadic inputs (random + integer-grid degeneracies).
  - **Found + fixed a real `orient3_sos` bug** while building them: its degenerate branch
    used a hand-coded 8-minor sequence that MIScomputed the +ε SoS on some
    **coplanar-DISTINCT** configs (verified against an independent finite-ε exact-rational
    oracle: points on x=2 gave −1 where the true +ε limit is +1) — reachable in the
    perturb=false kernel, masked by the perturb default + the recover gate + seed retry.
    Routed its degenerate branch through the canonical `_orient3_sos_exact` (same fix
    pattern as the earlier orient2_sos audit fix); now all four predicates break ties by
    ONE canonical +ε evaluator. Regression-pinned; suite re-verified.
  - **Exact-coordinate 3-D Delaunay kernel SHIPPED + verified** (`src/ExactMesh3D.jl`):
    a valid Delaunay tetrahedralization on exact `Rational{BigInt}` coords via a bounded
    super-tetrahedron + exact Bowyer–Watson (only the verified `orient3_rat`/`insphere_rat`
    + shared +ε SoS; no ghost tets, no coplanar-in-circle sqrt, no jitter). The exact
    predicates evaluate their determinant sign by a **fraction-free Bareiss** integer
    algorithm (clear row denominators → exact O(n³) integer determinant, no gcd/fraction
    growth) — verified sign-identical to the rational Laplace form, much faster (n=80
    exact Delaunay ~1 s). Verified: exact
    empty-circumsphere (0 violations), valid closed-manifold positive-volume mesh (χ=2),
    exact box fill (volume 8), and it breaks the maximally-cospherical 3×3×3 grid cleanly
    where the Float64 perturb=false kernel emits degenerate tets. Regression-pinned
    (`mesh3d_test.jl`), suite green. **This is the foundation the tracker always named
    as the recovery blocker.** First payoff: **`tetrahedralize_conforming_exact`** — the
    exact-kernel conforming multi-region mesher **conformingly meshes the cospherical
    2×2×2 unit-box assembly** (48 tets, valid, all 8 regions filled, exact volume 8,
    max-face-incidence 2) that the Float64 `tetrahedralize_conforming` must raise its
    blocker on (regression-pinned).
  - Boundary-Steiner recovery ON that kernel is the remaining piece and is genuinely
    **TetGen-class** (measured, not assumed): the exact Delaunay of the twisted prism's 16
    vertices recovers 18/28 facets (10 missing lateral faces); **naive facet-centroid
    Steiner insertion diverges** (each split spawns 3 coplanar sub-facets that are also
    missing ⇒ growth), so it needs the encroachment-ordered CDT / column decomposition
    with robust per-column boundary-Steiner (the majority of the 6 columns are non-star,
    needing on-face Steiner points — exactly the `ExactMesh` output the design requires).
    **Ten distinct recovery approaches were built + measured** — naive facet-centroid
    Steiner (diverges), per-column fan/star (all 6 columns non-star, all throw),
    **Rivara longest-edge bisection with the correct triangulation-independent conformity
    check** (boundary-area vs surface-area *fluctuates* 29.9→38.0→36.7 around 37.2 but
    never converges over 6 refinements), and a **Gabriel-encroachment loop** (circumcenter/
    diametral split — the right *criterion*, but the simplified batch form diverges
    20→36→81→161→255→382 without subsegment protection + one-at-a-time insertion ordering +
    the anti-encroachment rule) — none conform. The divergence of the simplified encroachment
    loop pinpoints *why* only the full engine works. An **11th attempt implemented the
    CORRECT algorithm** — segment-first Gabriel-encroachment refinement, one feature at a
    time, reject-circumcenter-that-encroaches-a-subsegment (Shewchuk's provably-terminating
    conforming-Delaunay), with encroachment checked against the vertex set (no re-Delaunay
    in the loop). It runs without diverging, but the naive exact-arithmetic
    O(features·verts)/step form is too slow (>3.6 s/step) to confirm termination
    in-session. **So the class is closeable-in-principle with the now-identified correct
    algorithm; what remains is the shippable engine = spatial acceleration of the
    encroachment checks + integration (recover_boundary on the exact kernel → `ExactMesh`
    output + exact certificate + interior classification + tests) — the genuine
    multi-session build.** Four more approaches (12–15) were then built + measured: K-layer
    twist subdivision + Delaunay (boundary area diverges, non-conforming); vertical-slab
    fan-tetrahedralization (valid+closed but volume grows with K — the finer surface, not
    the original); per-column centroid fan; and a **decisive per-column star test via
    faceted-vs-fan volume — ALL 6 columns are non-star** (fan volume ≠ column faceted
    volume, e.g. col 1: 1.046 vs 1.000). So every column provably needs a *boundary*-Steiner
    point (split a twisted face), which is exactly the exact-coordinate boundary-Steiner
    engine. **Fifteen distinct recovery approaches, each implemented + verified with running
    code; the conclusion is firm — this class needs the full boundary-Steiner CDT engine
    (spatial-index-accelerated encroachment refinement + ExactMesh output), a bounded but
    genuine multi-session build on the shipped exact kernel.** Confirmed: this class needs the
    full TetGen constrained-Delaunay engine (cavity retriangulation + subsegment-
    encroachment protection), a genuine multi-session build. Until it lands, the class
    keeps its **safe explicit blocker** (never a silent bad mesh).
  - (earlier) Exact-rational Delaunay + boundary-Steiner recovery design vetted (multi-agent
    feasibility workflow). Verdict on the pinned U/40 twisted prism: closing it to the
    exact-conformity CRC bar is **NOT possible with a Float64-coordinate output** — the
    locked columns need *boundary*-Steiner splits, and there is no Float64-representable
    point exactly on a cos40°/sin40° twisted edge/face (measured: 4/8 twisted lateral
    edges have no Float64-exact on-edge midpoint), so a rounded boundary Steiner point
    leaves the feature and the exact gate correctly rejects it. It **IS closeable** via
    an **`ExactMesh`** output: Steiner nodes carry exact `Rational{BigInt}` coords and the
    conformity/validity certificate runs on those stored rationals (no rounding between
    certifying and returning ⇒ bit-for-bit certified). Concrete pieces vetted: the exact
    ghost/coplanar in-circle is done sqrt-free by dropping the dominant-normal axis then
    `incircle_rat`; interior Steiner points are safe (margin ≫ ulp, verified per point
    with `orient3_rat`). Remaining build ≈ exact kernel + Gabriel-encroachment Steiner
    recovery + `ExactMesh` type/certificate (~500–800 LOC exact-arithmetic), in verified
    increments. Until it lands, the class keeps its **safe explicit blocker** (never a
    silent bad mesh, regression-pinned).
- **2026-08-13 — general curved-surface uniform sizing SHIPPED (`mesh_sized` +
  `refine_to_size`).** The last open meshing item. After confirming (step-by-step) that
  the existing sizers error/blocker on a fine sphere (cospherical) and that a direct
  Delaunay+lattice fill gives valid/watertight/conforming meshes but `maxedge` ~1.5–3×
  `hmax`, the fix is a **longest-edge bisection terminator** (`refine_to_size`,
  `Mesh3D.jl`): repeatedly split every tet incident to the globally-longest edge at its
  midpoint. Splitting *all* tets around one edge keeps the mesh **conforming** (shared
  faces split identically); the midpoint lies on the edge ⇒ **validity + total volume
  preserved**; and a bisection only creates edges no longer than the one it splits
  (median inequality) ⇒ the global max is non-increasing and it **terminates** at
  `maxedge ≤ hmax`. `mesh_sized(surface; hmax)` = fill (`tetrahedralize`, or
  `recover_boundary` for non-convex) then `refine_to_size`. Verified over the battery —
  convex box (exact vol 8), sphere, conical frustum, non-convex L-prism (exact vol 3) —
  all **`maxedge ≤ hmax`, valid, watertight, conforming**, sub-second; regression-pinned
  (`mesh3d_test.jl`). Box/cylinder/prismatic/general uniform sizing are all now shipped.
- **2026-08-13 — "complete ALL" pass: #11 fixed, all 22 HFSS geometries meshed
  natively, sliver-exudation driver, sizing assessed.**
  - **#11 (cospherical Delaunay perf) fixed** — re-measured (the tracker's "cavity not
    location" root-cause was wrong: Σlocate 43.9 s vs Σcavity 4.2 s vs Σretri 0.006 s at
    nθ=12,nz=16), fixed by jump-and-walk start + a `48·∛n` walk-guard cap in
    `Mesh3D.locate3`/`_pick_start3` (output-identical, `.sha`-verified over a 120-case
    fuzz). Fine pin 262 s → 22 s (nz=40); whole suite 11m44s → ~3m. Regression-pinned.
  - **#12 meshing complete — all 22 HFSS case geometries meshed natively from scratch**
    (no gmsh/OCC): `validation/hfss_cases/hfss_case_meshes.jl` builds each with Tessella
    primitives / raw surfaces (frustum, sphere, annulus, bowtie added), **22/22 valid +
    watertight + conforming**, regression-pinned in `test/hfss_cases_test.jl`. Verified
    two ways: a 21-agent build workflow (20/21 ok, the SAR sphere honestly reported as a
    recover_boundary cospherical-hang) and an independent load-back-and-validate of every
    `.msh`, then consolidated into one reviewed module (SAR sphere fixed by filling via
    `tetrahedralize`, whose perturb=true Delaunay + ray-cast avoids the hang).
  - **`remove_slivers` shipped** (`Optimize.jl`): converging exudation driver
    (smooth_optimize to convergence, validity+volume-gated) — 278→158 slivers on a random
    cloud, valid, volume-preserving; pinned in `optimize_test.jl`.
  - **Exact-coordinate kernel sped up ~20–40× (`ExactMesh3D.delaunay3d_exact`).** It was
    O(n²): a brute-force cavity scan tested the expensive `insphere_rat` against *every*
    live tet per insertion, with no point location. Added a **conservative Float64
    circumsphere pre-filter** (`_fcircum`/`_maybe_in_sphere`): skip the exact test on a tet
    only when the query is *safely* outside its float circumsphere (borderline + degenerate
    tets keep the exact test), so the result is **provably unchanged** — verified
    byte-identical vs the full-scan reference AND `is_delaunay_exact` = 0 violations on
    random + degenerate-grid inputs; measured **19.5× (n=120), 41× (n=200)**, speedup
    growing with n. This accelerates every exact-kernel path (`tetrahedralize_conforming_
    exact`, `mesh_sized_cdt`, `recover_boundary_cdt`). Regression-pinned (`mesh3d_test.jl`,
    n=120 exercises the pruning; oracle gates correctness).
  - **Uniform sizing on extruded/prismatic domains SHIPPED (`mesh_sized_extrude`).**
    The general-surface route (exact-midpoint flat-facet refinement + `mesh_sized_cdt`)
    was measured **correct but impractically slow** — a refined polyhedral surface is
    heavily coplanar, the worst case for the exact kernel (flat tets have huge
    circumspheres, so even the new pre-filter can't prune) — so it is NOT shipped.
    Instead, for the large **extruded/prismatic** class the degeneracy is sidestepped
    entirely: the 2-D size-controlled Ruppert mesher meshes the cross-section (edges
    `≤ hmax/√2`), then it is extruded into sized layers and each prism split into 3 tets
    by a **column-global-index diagonal rule** (every shared quad face picks the same
    diagonal in both neighbours ⇒ conforming). Meets the acceptance test on the
    **non-convex L-prism** and a **holed (genus-1) square-annulus prism** over an `hmax`
    sweep: `maxedge ≤ hmax`, `validate.ok`, **exact volume** (boundary preserved, no
    jitter), watertight, conforming, sub-0.03 s — regression-pinned (`mesh3d_test.jl`).
    So box (`mesh_box`), cylinder (`mesh_cylinder`), and prismatic (`mesh_sized_extrude`)
    uniform sizing are all shipped; the general curved-surface terminator was then closed
    the same day (see the `mesh_sized`/`refine_to_size` entry above).
- **2026-08-13 — split ASCENT tracking into `ASCENT.md`.** Per user, all ASCENT-facing
  work (external FEM solver integration, mesh/BC handshake, solve proofs, cavity
  eigenmode, case 9.2 solve, the 22-case HFSS solve campaign, and the HFSS reference
  audit) moved out of `STATUS.md` into a dedicated `ASCENT.md`, so `STATUS.md` now tracks
  the Tessella package (the mesher) only. The current-state sections here are package-only
  with pointers to `ASCENT.md`; the Log retains historical context.
- **2026-08-13 — deep-debug audit: 8 real bugs found + fixed (every one reproduced,
  fixed, and CRC-regression-pinned).** A fresh 12-area adversarial audit (parallel
  agents + own verification; every finding reproduced against running Julia before any
  fix). Suite green after: **151,438 assertions**.
  1. *(critical, `ExactMesh3D`)* the Float64 circumsphere pre-filter was **non-conservative
     for far-from-origin coordinates** — `Float64(b)−Float64(a)` catastrophically cancels
     at magnitude ≳1e13 with fractional coords, so it wrongly skipped cavity tets and
     `delaunay3d_exact` emitted a non-Delaunay/non-conforming mesh (the exact "overlapping
     facets" failure the package exists to prevent). Fixed: run the pre-filter on **per-axis
     origin-shifted** coords (shift exact in rational ⇒ small, well-conditioned Float64),
     and disable it for absurd (>1e11) extents (full exact test). Verified 0 invalid to
     offset 1e16, byte-identical near-origin.
  2. *(major, `RecoverCDT.certify_exact`)* the per-region area certificate filtered boundary
     faces by **plane, not region membership**, so two disconnected coplanar regions (a
     U/plus/slot/multi-component surface) double-counted and `recover_boundary_cdt` /
     `mesh_sized_cdt` threw a **spurious area-mismatch blocker on valid conforming meshes**.
     Fixed with a centroid `in_region_closed` scope; U-channel now recovers.
  3. *(major, `Tessella.mesh_sized`)* silently returned an **over-filled** mesh for a
     non-convex (Schönhardt) domain — `tetrahedralize` caps a concavity, still valid +
     watertight but the wrong domain. Fixed: require the fill's **boundary area == input
     surface area**, else fall through `recover_boundary → recover_boundary_cdt`.
  4. *(minor, `mesh_sized`)* false-rejected a single-tet domain (max face incidence 1, not
     2); the manifold test is now `≤ 2`.
  5. *(minor, `refine_to_size`)* dropped `tet_tag` — multi-region meshes lost their
     partition. Now each child inherits its parent tet's tag.
  6. *(minor, `mesh_cylinder`/`cylinder_surface`)* a zero-length axis produced a silent
     NaN mesh; now throws `ArgumentError`.
  7. *(minor, `IO._read_physical_names!`)* collapsed runs of interior whitespace in a
     `$PhysicalNames` name; now parses the quoted substring verbatim.
  8. *(minor, `CAD.project_to(::CylinderS)`)* on-axis projection crashed for an x-parallel
     axis; now references the least-aligned coordinate axis.
  Not bugs: `mesh_boolean` on a tilted-axis cylinder throws its safe watertight blocker
  (never a leaky surface — correct fail-safe); `_ray_hits_tri`'s absolute `1e-12` cutoff
  degrades only at sub-picometer scale (out of the mm/m domain; documented, both classifiers
  degrade identically so their pinned equivalence holds). Predicates + MeshTypes oracles:
  audited clean.
- **2026-08-13 — deep-debug PASS 2 (audit the fixes) + 1 pre-existing bug, resolved in PASS 3.** A second
  fresh pass over the 8 fixes + a public-API sweep.
  - **Fixed: `CAD.project_to(::CylinderS)` — the pass-1 on-axis fix was INCOMPLETE.** For an
    axis with irrational unit components (e.g. `(1,1,1)/√3`), the foot-point rounding leaves a
    tiny axis-PARALLEL residual, so `rl==0` misses the on-axis case and the normal branch
    normalises that parallel residual → the projected point moves ALONG the axis (radial
    distance 0, residual −r). Fixed with a **scale-relative** on-axis tolerance
    (`rl ≤ 1e-10·(|axl|+r+1)`); verified residual <1e-11 for all axes at multiple on-axis
    points + off-axis still exact. `cad_test` green (8078). Pinned. Pass-2 otherwise clean
    (certify_exact fix, refine_to_size tags, ExactMesh3D shift edge cases, mesh_sized curved,
    mesh_box_regions/mesh_boolean/tetrahedralize_multi/p2 all re-verified correct).
  - **RESOLVED in PASS 3: `ExactMesh3D.delaunay3d_exact` returned 0 tets
    for a near-coplanar / high-aspect point set.** The finite super-tetrahedron (`K=10`,
    lines ~84-95) is too small: a flat/sliver tet has a circumsphere of radius ~base²/height
    that strictly contains the super-vertices, so the real tet is dropped by the final
    super-vertex filter ⇒ a **silent empty (vacuously "Delaunay") result**. Triggers at aspect
    ratio ≳ 20. Repro: `delaunay3d_exact([(0,0,0),(1000,0,0),(0,1000,0),(500,500,1)])` → 0
    (expect 1); aspect sweep 2/10 → 1 tet, 20/50/100/500/1000 → 0. **Reachability is limited**
    — a thin *box* (20×20×0.2, aspect 100) via `recover_boundary_cdt` gives the correct 6 tets
    (vol 80), so the common thin-feature paths route through configurations that work; it
    triggers on sparse-off-plane near-coplanar clouds. **Fix direction (bounded + correct):**
    a coplanarity pre-check (all-coplanar ⇒ return empty, which is correct) then **retry the
    build with a geometrically growing super-tet** (K×=4 each attempt, exact rational so no
    overflow) until **every input point appears in some output tet** (a correct 3-D Delaunay
    uses all points); cap the growths and raise an explicit blocker rather than ever return a
    silent partial/empty mesh for a non-coplanar input. PASS 3 implemented and verified this
    bounded construction and completeness gate.
- **2026-08-14 — deep-debug PASS 3, exact-kernel completeness fix.** Reproduced the
  PASS-2 bug before editing: aspect 2/10 produced
  one tet while 20/50/100/500/1000 produced **zero**, and the empty output passed
  `is_delaunay_exact` vacuously. The direct cause is now pinned in code: growing only
  the positive super-vertices cannot release the fixed negative super-vertex from a
  sliver's huge circumsphere. Fixed `delaunay3d_exact` with a first-pass-compatible
  K=1 build, exact affine-conditioning fallback, all-four-vertex geometric growth,
  bounded retries, physical-orientation/duplicate/face-incidence/vertex-completeness
  certification, and explicit rejection of duplicate inputs. The 4-point independent
  oracle passes at aspect 20, 1000, and 1,000,000 (9/9 assertions); coplanar input still
  correctly returns no volume cells and duplicates raise `ArgumentError`; the full
  `mesh3d_test.jl` passed **410/410** under `--check-bounds=yes` (9m23s). Allocation
  measurement on a 40-point exact build: **17,920,376 B before vs 18,161,840 B after**
  (+241,464 B / +1.35%), with cavity scratch dictionaries now reused per insertion.
  Full post-fix package gate: **151,464/151,464** assertions passed under
  `--check-bounds=yes` on Julia 1.12.6 (20m16s). The cospherical performance gate now
  measures process CPU time so unrelated scheduler contention cannot create a false
  regression; the previously failing run measured 76.13s wall but only 9.26s CPU.
- **2026-08-14 — deep-debug PASS 4, top-level PLC conformity.** Reproduced a silent
  wrong-domain result from `mesh_volume`: a reflex twisted prism returned a valid,
  watertight convex-hull cap of volume **1.1160253971** instead of the PLC volume
  **0.8660254038**. `mesh_volume` and `mesh_sized` now share a conforming-fill path:
  the fast restricted-Delaunay result must pass the exact region/crease/boundary gate,
  otherwise boundary recovery and exact conforming-Delaunay recovery are tried, with
  all blocker reasons retained. `check=false` intentionally preserves the documented
  caller-responsibility bypass. Regression result: the twisted prism validates at
  `sqrt(3)/2`; pipeline **42/42**, Mesh3D **410/410**, and the full package suite
  **151,466/151,466** pass under bounds checking on Julia 1.12.6 (12m53s).
