# Tessella.jl status

ASCENT integration and the external HFSS full-wave campaign are tracked separately in
[`ASCENT.md`](ASCENT.md).

## Status

The original Stage 0–6 simplex-mesher roadmap is complete and remains the regression
baseline. The active roadmap now targets independent Gmsh 4.15.2 parity and is **not
complete**. Work is ordered by ASCENT meshing value before UI and post-processing.

| Stage | Capability | State |
|---|---|---|
| 0 | exact predicates, mesh types, validation, MSH/STL I/O | DONE |
| 1 | 2-D Delaunay, constrained Delaunay, quality refinement | DONE |
| 2 | size fields, graded curves, planar/cylinder/parametric surfaces | DONE |
| 3 | 3-D Delaunay, exact-coordinate kernel, volume fill, conforming recovery | DONE |
| 4 | uniform sizing, quality metrics, flips, smoothing, sliver reduction | DONE |
| 5 | healing diagnostics, native primitives, analytical CAD, imprints, mesh CSG | DONE |
| 6 | globally certified P2 tetrahedra and solver-consumable I/O | DONE |

### Active parity increment

| Track | State | Verified implementation increment |
|---|---|---|
| P1 | **IN PROGRESS** | Native scalar/anisotropic catalog, strict `.geo` field graph with injected model/view context, Gmsh-style 1-D policy, and field/entity-aware 2-D, surface, and 3-D refinement |
| P2 | **IN PROGRESS** | 125 fixed-node Gmsh types plus ten serializable cut/border/child/sub-element records, mixed blocks/entities/classification metadata, structural validation/CRC, and ASCII/binary MSH v2.2/v4.1 read/write |
| P3 | **IN PROGRESS** | Native analytical surfaces/imprints, closed box/cylinder/cone/geodesic-sphere primitives, cavities, mesh Boolean CSG, finalized-mesh affine transforms, and bounded `.geo` constant expressions/ranges |
| P4 | **IN PROGRESS** | Recombination, uniform refinement, curve laws, planar triangle/quad transfinite patches, affine five-/six-face transfinite volumes, and recombined hexahedra |
| P5–P6 | **PENDING** | No state change |

P1 does not claim boundary-layer element topology, the full Gmsh automatic-sizing
pipeline, high-order/custom-interpolation, multiple-time-step, or mixed-component
`PostView`, materially warped quadrangles, `PostView` tensor-to-metric evaluation, direct tensor or
metric-meshing parity, full `.geo`/CAD-model execution, or exact CAD distance. P2 does
not claim general mixed-element generation or recombination beyond P4's first-order
surface pairing, MINI basis-selector tags 138/139 as mesh records, integration of
mixed blocks into the simplex meshing kernels, curved high-order
Jacobian certification, ancillary/unknown-section preservation (binary readers reject
unsupported sections explicitly), repeated `$Nodes`, non-8-byte binary data, internal
indices beyond `Int32`, or lossless multi-physical-group MSH v2.2 projection. Some
registered fixed tags and polygon-border type 69 require explicit Tessella-only output
because Gmsh 4.15.2 cannot consume them safely. MSH2 ASCII preserves variable records
and parent/domain links; binary MSH2 and MSH4 have explicitly narrower special-record
contracts. Pinned Gmsh 4.15.2 corrupts distinct parent links in its own binary MSH2
rewrite, and nonzero-physical special MSH4 requires compatible node/entity
classification metadata for a safe rewrite. P3 does not yet claim a general entity kernel,
OpenCASCADE/BREP/NURBS, CAD import/export, transformations of analytical/CAD
entities, or full `.geo` execution. Its scanner handles finite arithmetic constants,
pure numeric functions, prior scalar bindings, explicit field/physical tags, and
finite constant ranges in recognized numeric field lists/selectors. Entirely numeric
Physical memberships are range-checked but remain geometry data. It rejects loops,
macros, dynamic tags, option reads, stateful functions, dynamic/general ranges,
logical/ternary evaluation, CSG statements, and mixed geometry-derived
physical-group right-hand-side evaluation.

P4's uniform-refinement slice applies the exact Gmsh 4.15.2 linear segment, triangle,
and tetrahedron child templates while sharing edge midpoints, compacting unused nodes,
and preserving parent tags. Its straight-curve slice covers normalized affine-line
Progression/Power, Bump, and Beta parameters. Its surface transfinite slice covers
already-discretized, count-matched, three- and four-sided planar chains using Gmsh's
specific triangular and average-chord Coons interpolation. Four-sided grids can also
be emitted as first-order Gmsh type-3 quadrangles with exact projected
corner-Jacobian certification. Affine eight-corner blocks use Gmsh's unrecombined
six-tetrahedron transfinite volume subdivision; canonical affine triangular prisms
use its legacy collapsed-grid five-face tetrahedral path. P4 does not yet claim Gmsh's Blossom/full-quad
algorithms, non-affine CAD curve integration, FlexibleTransfinite or HWall/size-map
curve laws, quasi-transfinite or holed patches,
general CAD parameterizations, recombined three-sided patches, curved/warped or
compact-TransfiniteTri volumes, volume/hybrid
recombination, selective or high-order refinement, coarsening, boundary-layer element
topology, or periodic/embedded model constraints.

OpenCASCADE/NURBS, remaining algorithms/fields, broad formats and API, GUI, and
post-processing are unfinished parity tracks, not project non-goals.

## Current worktree verification

Re-measured on 2026-08-21 with Julia 1.12.7 after the transfinite-hex, geo-range,
entity-kernel, NURBS, `.geo` execution, boundary-layer, periodic, API/CLI/GUI/post,
and P6 box-API increment. Both bounds-checked package runs matched:

- `julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 163,096/163,096 assertions passed twice (11m08.3s, then 10m40.8s).
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  — exited 0 against Gmsh 4.15.2-git. Exact flat-model volumes box=2, tunnel=24,
  hollow box=35; cylinder prism 62.652572; enclosure gmsh empty solids reproduced.
- Size-field child: `SIZE_FIELD_DIFFERENTIAL_OK gmsh=4.15.2 plugin_calls=23
  direct_cases=23 direct_samples=63 mesh_cases=5 context_skips=5`.
- Geo-range child: `GEO_RANGE_DIFFERENTIAL_OK gmsh=4.15.2-git float_cases=13
  integer_cases=4 wrapped_cases=2 samples=58 bit_exact=1`.
- Transfinite hex child: `TRANSFINITE_HEX_DIFFERENTIAL_OK gmsh=4.15.2-git cases=8
  nodes=288 hexahedra=96 boundary_quadrangles=256 max_node_error=9.50e-12`.
- P6 box API child: `GMSH_PARITY_BOX_OK gmsh=4.15.2 tessella_volume=1
  gmsh_tets=1158 tessella_tets=12`.
- Fresh-process `using Tessella` plus `mesh_volume`/`size_at` succeeded twice on
  `box_surface(0,1,0,1,0,1)`: 9 nodes / 12 tets, `size_at==0.5`, both `validate` ok.
- Focused hex CRC 142/142; IO 305/305; NURBS 20/20; entity/`.geo` 15/15;
  boundary-layer 7/7; periodic 4/4; API/CLI/GUI/post 15/15.
- `git diff --check` passed. `.grok/` is gitignored and absent from the index.

Previous aggregate on the same day with Julia 1.12.7, kept as historical:

- `julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 161,183/161,183 assertions passed in 9m05.3s.
- `julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl`
  — exited 0 against Gmsh 4.15.2-git. Tessella preserved the exact flat-model
  volumes (box 2, tunnel 24, hollow box 35), completed the curved-model comparisons,
  and reproduced the literal enclosure's non-zero Gmsh exit with zero volume
  tetrahedra.
- The required size-field child reported
  `SIZE_FIELD_DIFFERENTIAL_OK gmsh=4.15.2 plugin_calls=23 direct_cases=23
  direct_samples=63 mesh_cases=5 context_skips=5`; the context skips are explicit
  non-claims listed in `validation/size_fields/STATUS.md`.
- Focused bounds-checked gates passed 6,915/6,915 size-field assertions,
  2,398/2,398 fixed-node/mixed-element assertions, and 74/74 Mesh1D assertions.
  The Mesh1D gate also passed under Julia 1.11.
- `detect_ambiguities(Tessella; recursive=true)` and
  `Base.Docs.undocumented_names(Tessella; private=false)` both returned zero.
- `git diff --check` passed.

After that stable aggregate gate, the isolated native-primitive increment passed
69/69 bounds-checked geometry assertions under both Julia 1.12.7 and Julia 1.11,
including deterministic CRCs, analytical/polyhedral volume checks, direct volume
meshing, resource/error paths, and linear allocation-growth ratchets. It will be
included in the next aggregate package gate. The subsequent finalized-mesh transform
increment passed 56/56 bounds-checked assertions under Julia 1.12.7 and Julia 1.11.9.
Its independent Gmsh 4.15.2 `model.mesh.affineTransform` oracle matched all four
fixture nodes exactly; 10,000- and 20,000-node translations allocated 504,984 and
996,504 bytes respectively, and the method-ambiguity scan remained empty. The binary
mixed-MSH increment then passed 2,579/2,579 bounds-checked assertions under both Julia
1.12.7 and Julia 1.11.9, including native and opposite-endian MSH 2.2/4.1, full-width
v4 tags, parametric nodes, atomic/resource failures, and Gmsh 4.15.2 acceptance.
Reading 2,000 and 4,000 nodes allocated 1,670,704 and 4,239,472 bytes respectively.
The following recombination increment passed 44/44 bounds-checked assertions under
Julia 1.12.7 and Julia 1.11.9. Its 12×12 grid produced 144 quadrangles with CRC
`dbb1bf17965d4e011e7f51a452c6a03e4018a628ffc7e7d9b33d9fc6b922439f`; 20×20 and
40×40 grids allocated 4,073,984 and 17,348,496 bytes. Gmsh 4.15.2 accepted both the
ASCII and binary recombined MSH 4.1 fixtures with `-check -parse_and_exit`.

The subsequent `.geo` constant-expression increment passed 231/231 bounds-checked IO
assertions under Julia 1.12.7 and Julia 1.11.9. An independent Gmsh 4.15.2 expression
oracle matched 34 accepted expressions exactly and matched seven error cases; an API
oracle also matched explicit expression-derived physical and field tags. Unsupported
control-flow contexts are rejected when they can affect relevant statements, and
their scalar bindings are invalidated before later use.

The next PostView increment passed 6,953/6,953 bounds-checked size-field assertions
under Julia 1.12.7 and Julia 1.11.9. It covers first-order scalar/vector point, line,
triangle, quadrangle, tetrahedron, hexahedron, prism, and pyramid list data; tensor
views preserve Gmsh's scalar `MAX_LC` result. The Gmsh 4.15.2 pointwise oracle covered
160 samples with maximum absolute error `6.66e-16`, plus eight exact closest-node
fallbacks. Warm 10,000-query loops for quadrangle, hexahedron, prism, pyramid, and
vector-quadrangle fields allocated at most 64 bytes in total.

The subsequent uniform-refinement increment passed 78/78 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 API differential matched
the ordered 2/4/8 segment/triangle/tetrahedron child templates and physical tags; the
combined fixture produced 10 nodes and CRC
`db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031`.
Refining 2,000 and 4,000 connected segments allocated 423,504 and 845,648 bytes,
respectively (1.99679×). The focused module ambiguity and public-doc scans both
returned zero.

The following special-element increment passed 2,822/2,822 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. It covers types 34/35/67/68/69/70/133–136,
compact variable connectivity, parent/domain references, validation and CRC, MSH2
ASCII and fixed-width binary records, and fixed unlinked MSH4 records. Pinned source
hashes and installed Gmsh 4.15.2 probes cover the accepted formats and their explicit
blockers. Rejected 30,000- and 300,000-entry variable records allocated 5,328 and
5,120 bytes after warm-up with `max_connectivity=0`; accepted 2,000- and 4,000-record
fixtures allocated 13,000,112 and 26,462,288 bytes. The tests pin Gmsh's type-69
normal-lifecycle crash and binary distinct-parent rewrite corruption as external
limitations instead of claiming unsafe compatibility.

The subsequent transfinite increment passed 91/91 bounds-checked assertions under
Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 API differential matched all
four triangle arrangements and 80 boundary/interior coordinate samples with maximum
absolute node error `4.44e-16`; it also pins Gmsh's mismatched-opposite-side fallback
at 18 nodes/23 triangles and its holed-surface error. The 64×64 and 128×64 allocation
fixtures used 1,771,392 and 3,342,496 bytes under Julia 1.12.7. Focused ambiguity and
public-documentation scans returned zero.

The subsequent straight-curve increment passed 310/310 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential covered
Progression/Power, Bump, and Beta in 21 cases and 219 coordinates with maximum
absolute error `5.49e-8`; its deterministic parameter SHA-256 is
`b525628945d55f49bc5151ad313d697f9c32f7b3325015e5d0080a69260ffd0e`.
The 20,000- and 40,000-node fixtures allocated 163,904 and 327,744 bytes for each
law; focused ambiguity and public-documentation scans returned zero.

The subsequent three-sided transfinite increment passed 103/103 bounds-checked
assertions under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential
matched exact boundary/triangle topology and 150 nodes across four arrangement names,
four resolutions, and planar/tilted geometries, with maximum coordinate error
`2.81e-15`. Its deterministic mesh SHA-256 is
`5f231f0c22f6812c247514103e21653e15281194911bf456b4c00a9d190a40df`.
The 64- and 128-division fixtures allocated 2,383,144 and 6,154,008 bytes under Julia
1.12.7; focused ambiguity and public-documentation scans returned zero.

The subsequent affine-volume increment passed 74/74 bounds-checked assertions under
Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential matched all 144
tetrahedra, 128 boundary triangles, and 72 mapped nodes across two affine blocks, with
maximum coordinate residual `9.51e-12`. The 8×8×4 and 16×8×4 fixtures allocated
385,984 and 774,128 bytes under Julia 1.12.7. Focused ambiguity and public-documentation
scans returned zero.

The subsequent recombined-quadrangle increment passed 128/128 bounds-checked
assertions under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential
matched ordered type-1/type-3 connectivity and 151 nodes across four arrangement
names, four resolutions, and planar/tilted geometries, with maximum coordinate error
`5.662137425588298e-15`. Its deterministic mixed-mesh SHA-256 is
`d05e9cbc57de975f1f88d8fb1da62e0636ea4bbe63a038d3ca060b16201b0ea2`.
The 64×64 and 128×64 fixtures allocated 3,430,528 and 6,274,688 bytes under Julia
1.12.7. Exact-rational fallback regressions cover valid thin affine patches whose
Float64 plane-normal or normalized-projection calculations cancel; focused ambiguity
and public-documentation scans returned zero.

The subsequent five-face-prism increment passed 130/130 bounds-checked assertions
under Julia 1.12.7 and Julia 1.11.9. Its required Gmsh 4.15.2 differential matched
the exact canonical tetrahedron vertex order and boundary topology for 117
tetrahedra, 106 boundary triangles, and 63 nodes across three affine prisms; maximum
coordinate error was `6.840577468974138e-12`. The deterministic mesh SHA-256 is
`a16de779890f62f8a09d928cbef67a6f13b09c6765a7d91ce8e86de78c14db6e`.
The 24×12×6 and 48×12×6 fixtures allocated 1,997,056 and 4,008,256 bytes under
Julia 1.12.7. Exact-predicate and exponent-scaled volume gates reject represented
folds or material loss while accepting verified subnormal affine prisms; focused
ambiguity and public-documentation scans returned zero.

Earlier measurements on 2026-08-14 with Julia 1.12.6, kept because they were
not re-run in this gate:

- `julia --project=. --check-bounds=yes validation/run_all.jl`
  — passed in 182.30 s against Gmsh 4.15.2-git, including exact box/tunnel/hollow-box
  volumes, curved-reference comparisons, and the literal enclosure's reproduced
  non-zero Gmsh exit with zero volume tetrahedra.
- Focused gates passed: 983/983 field/refinement assertions, 66/66 existing 1-D and
  surface assertions, and 439/439 Mesh3D assertions.
- Final-code enclosure field acceptance with the literal graph clamped to
  `size_min=0.002`, `size_max=0.012` returned a valid mesh with 706,962 nodes,
  3,838,729 tetrahedra, and 113,350 tagged triangles in 680.61 s, with a
  1,172,013,056-byte maximum resident set. The fixture's literal 0.2667 mm minimum
  remains a separate full-resolution scalability gate and is not claimed here.

## Previous complete baseline gate

Verified on 2026-08-14 with Julia 1.12.6:

- `julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'`
  — 151,787/151,787 assertions passed in 6m54.1s.
- `julia --project --check-bounds=yes validation/run_all.jl`
  — completed against gmsh 4.15.2-git. Tessella matched the exact model volumes for
  the box (2), box tunnel (24), hollow box (35), and 48-gon cylinder prism
  (62.652572). The external enclosure fixture reproduced gmsh's non-zero exit and
  zero tagged volume tetrahedra.
- A mixed segment/triangle/tetrahedron MSH v4.1 written by Tessella passed
  `gmsh <file> -check -v 3` with exit 0.
- `using Test, Tessella; detect_ambiguities(Tessella; recursive=true)` reported 0
  method ambiguities. `Base.Docs.undocumented_names(Tessella; private=false)`
  reported 0 undocumented public names.
- `git diff --check` passed, and the source/test scan found no `TODO`, `FIXME`,
  `WIP`, placeholder, or unimplemented markers.

The last runtime-bearing baseline head for this gate is
`27010e312455fe3dbb67ce0efd2f51800487c3ad`, and `origin/main` was checked to the
same SHA after its push. Later commits change documentation/docstrings only; their
package-load/focused checks and remote synchronization were also verified.

## Deep-debug closure

The cumulative audit repaired and regression-pinned the following classes:

- exact Delaunay completeness and fallback certification;
- top-level PLC/interface conformity and full triangle/tetrahedron vertex-link
  manifold checks, including pinches and duplicate cells;
- finite/range/resource contracts for public geometry, meshing, refinement,
  optimization, P2, and I/O APIs;
- exact coplanar 3-D circle decisions, consistent bounded SoS evaluation, robust
  extreme-magnitude area/volume signs, and Float64-resolution blockers;
- 2-D CDT/refinement midpoint, encroachment, duplicate-constraint, callback, and
  termination behavior;
- lower-dimensional topology/tag preservation during volume refinement;
- strict MSH section/entity/count/tag validation, atomic writers, and structural
  STL parsing/welding;
- exact global cubic-Bernstein P2 Jacobian certification and exact coefficient
  volume integration;
- analytical projection/imprint edge cases and generated-surface postconditions;
- interrupt propagation through conversion/overflow error paths.

The package's recover-or-block contract intentionally reports an explicit diagnostic
when a requested mesh exceeds representational or configured resource limits; that is
a safety contract, not a silent fallback. The closure statement applies only to the
historical Stage 0–6 scope, not to the active Gmsh parity objective.

## Memory evidence

The reproducible stress measurement uses
`refine_to_size(mesh_box(0,3,0,3,0,3; hmax=3), 0.25)` after one warm-up and `GC.gc()`.
Each measured snapshot produced 22,065 nodes, 98,304 tetrahedra, and a valid mesh:

| Version | Allocated bytes | Timed run |
|---|---:|---:|
| pre-fix `5d3e466` source snapshot | 75,038,688 | 0.3517 s |
| previous hardened implementation | 70,337,184 | 0.2972 s |
| previous field/refinement snapshot | 100,434,688 | 1.1318 s (concurrent enclosure load) |
| current P1/P2 worktree | 72,870,624 | 0.3027–0.3169 s |

The historical hardened version reduced allocation by 6.27% from `5d3e466`. The
previous field/refinement snapshot was 30,097,504 bytes (42.79%) above that version on
the same mesh. The current worktree is 27,564,064 bytes (27.44%) below that snapshot
and 2,533,440 bytes (3.60%) above the historical hardened implementation. All three
current runs produced 22,065 nodes, 98,304 tetrahedra, a valid mesh, and CRC
`7378bf6e460c596aa04f9f3b4a8cc9ce176a69f2386f253f7bd25057b551ceb9`.
The refinement queue keeps one live heap record per long edge via a companion set,
preventing duplicate queued records; the benchmark verifies the resulting mesh rather
than relying on allocation alone.

## Acceptance cases

- Exact predicates remain cross-checked against independent exact-rational oracles,
  including degenerate and randomized configurations.
- General boundary recovery covers the non-star/reflex twisted prism through the
  exact conforming-Delaunay path, while every returned mesh is independently checked
  for input-facet conformity and manifold topology.
- The literal enclosure/coax fixture is reconstructed natively with four filled
  material volumes and tagged boundary groups. The gmsh-impossible case is also
  solver-loadable and solved; those proofs live in [`ASCENT.md`](ASCENT.md).
- All 22 HFSS guide geometry classes have native, valid, watertight, conforming mesh
  regressions in `test/hfss_cases_test.jl`. Re-running the remaining literal cases as
  full-wave ASCENT studies is external solver work, not a missing Tessella feature.

## Provenance

The final hardening commits on `main` were pushed individually, including:

- `5ec9b4e` exact Delaunay completeness;
- `abbfbe1` top-level PLC conformity;
- `9776a51` safe surface diagnostics;
- `bc9f294` manifold topology and PLC fills;
- `a630ab4` 2-D kernel/refinement contracts;
- `9b3b0ac` curve/surface bounds;
- `0b1ec2e` geometry/P2/optimization/I/O contracts;
- `5d3e466` volume meshing, partition recovery, bounded SoS;
- `d4b0dc0` industrial closure hardening;
- `27010e3` bounded flip-optimizer passes and final source marker cleanup.

Historical closure claims above are tied to their dated commands and measurements.
Current P1/P2 claims are limited to the checked source, tests, differential, aggregate
validation, and memory evidence above; the explicit P1/P2 non-claims remain open.
