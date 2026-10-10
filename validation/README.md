# Validation — Tessella vs external mesh tools

Independent cross-validation uses pinned Gmsh 4.15.2 as a development oracle;
production meshing remains Julia-native. The geometric case report compares
represented domains, analytic volume, element quality and wall-clock time.
Focused differentials separately check public identities, ordered typed cells,
coordinates, ownership, raw and computed parameters, lifecycle state, diagnostics
and allocation behavior. Volume agreement alone does not establish those contracts.

Status reviewed **2026-10-10**: the authoritative [handoff](../HANDOFF.md#continuation-handoff-2026-10-10)
supersedes earlier completion statements. Published code and dated captures are
separate from **Root V24** and **Release V6**, which remain unqualified combined
candidates. Existing scoped Source proofs do not qualify the pending public,
resource, aggregate or package release. No validation was rerun for this documentation update.

## Layout

```
validation/
  run_all.jl             # aggregate: selected mandatory children and case report
  REPORT.md              # generated case report; git-ignored, not a release certificate
  support/
    common.jl            # helpers: gmsh runner, tet metrics, comparison
  size_fields/
    differential.jl      # required Gmsh 4.15.2 field differential
    STATUS.md            # exact coverage and explicit non-claims
  geo_ranges/
    differential.jl      # required bit-exact Gmsh 4.15.2 constant-range differential
  geo_control_flow/
    differential.jl      # required bit-exact Gmsh 4.15.2 .geo control-flow differential
  geo_transforms/
    differential.jl      # required Gmsh 4.15.2 .geo transform/Duplicata differential
  geo_extrude/
    differential.jl      # required Gmsh 4.15.2 .geo Extrude differential
  geo_curved/
    differential.jl      # required bit-exact Gmsh 4.15.2 arc/surface-filling differential
  geo_splines/
    differential.jl      # required bit-exact Gmsh 4.15.2 spline-family differential
  geo_primitives/
    differential.jl      # required Gmsh 4.15.2 OCC primitive-layout differential
  geo_constraints/
    differential.jl      # required Gmsh 4.15.2 .geo meshing-constraint/Delete differential
  uniform_refine/
    differential.jl      # required Gmsh 4.15.2 linear-simplex template differential
  mesh_affine_transform/
    differential.jl      # required whole-cache affine-transform differential
  mesh_data_queries/
    differential.jl      # required bulk and connectivity-derived query differential
  api_generate01/
    differential.jl      # required classified 0D/1D generation and sparse lifecycle differential
  api_point_duplicates/
    point_duplicate_refine_resources.jl # actual Point identity, refinement allocation and implementation AST
  mesh_entity_topology/
    differential.jl      # required automatic/manual global edge/face differential
  mesh_point_location/
    differential.jl      # required simplex location/reference-coordinate differential
  mesh_jacobians/
    differential.jl      # required linear-simplex Jacobian/reference-map differential
  mesh_quadrature/
    differential.jl      # required fixed-family reference-quadrature differential
  mesh_function_spaces/
    differential.jl      # required arbitrary-order nodal plus simplex hierarchy differential
  mesh_element_qualities/
    differential.jl      # required linear-simplex quality-measure differential
  element_catalog_queries/
    differential.jl      # required fixed type/property catalog differential
  high_order/
    differential.jl      # required Gmsh 4.15.2 type-11 tetrahedron differential
  transfinite/
    differential.jl      # required Gmsh 4.15.2 four-sided patch differential
  transfinite_curve/
    differential.jl      # required Gmsh 4.15.2 straight-curve-law differential
  transfinite_triangle/
    differential.jl      # required Gmsh 4.15.2 three-sided patch differential
  transfinite_quad/
    differential.jl      # required Gmsh 4.15.2 recombined-quad differential
  transfinite_volume/
    differential.jl      # required Gmsh 4.15.2 affine-volume differential
  transfinite_prism/
    differential.jl      # required Gmsh 4.15.2 five-face-prism differential
  transfinite_hex/
    differential.jl      # required Gmsh 4.15.2 recombined-hexahedron differential
  quadtri_nonew/
    differential.jl      # required isolated NoNew source-quad certificates and oracle
  gmsh_parity/
    box_api.jl           # Tessella API box volume vs analytic 1 and Gmsh 4.15.2
    geo_geometry_expressions.jl # bounded geometry-expression execution
    geo_list_variables.jl # bounded numeric-list and entity-reuse execution
    geo_mesh_sizes.jl # Point sizing and topology-derived Physical groups
    geo_dynamic_tags.jl # geometry/Physical allocation and lifecycle, SetMaxTag, OCC
    model_topology_queries.jl # entity/boundary/adjacency API differential
    model_entity_identity.jl # entity names and live-reference retagging
    model_entity_removal.jl # ordered dependency-safe recursive removal
    model_spatial_queries.jl # analytical bounds and containment queries
    model_entity_metadata.jl # native types, plane properties, partition ownership
    model_entity_evaluation.jl # evaluation, projection, and surface reparametrization
    model_entity_state.jl # visibility, color, Point coordinates, attributes
    boolean_boxes.jl      # Boolean snapshot ownership and Delete lifecycle
    nurbs_surface.jl      # OCC patch plus two-way IGES 126/128/144 interoperability
    periodic_translation.jl # native/projected periodic pairs vs Gmsh 4.15.2
    periodic_embedded_curve.jl # embedded periodic-curve MSH2/MSH4 lifecycle
    periodic_curve_graph.jl # dependency-graph/expression periodic lifecycle
    periodic_surface_volume.jl # planar periodic explicit-volume boundaries
    periodic_curve_branch.geo # one master reused by two embedded curves
    periodic_curve_chain.geo # acyclic master/slave dependency chain
    periodic_curve_expressions.geo # scalar/expression/range periodic chain
    periodic_curve_affine_expressions.geo # 16-entry expression affine map
    periodic_curve_rotate_expressions.geo # expression rotation map
    periodic_native.geo  # bounded native translation-periodic fixture
    periodic_two_direction.geo # shared-corner x/y-periodic fixture
    embed_point.jl        # classified Point-In-Surface MSH4 projection lifecycle
    embed_line.jl         # classified Line-In-Surface MSH4 projection lifecycle
    embed_sheet.jl        # Surface-In-Volume plus nested point/curve MSH lifecycle
    embed_sheet_hole.jl   # holed Surface-In-Volume MSH2/MSH4 lifecycle
    explicit_shell.jl     # planar Surface Loop/Volume MSH2/MSH4 lifecycle
    ...                   # focused API, CAD, and boundary-layer cases
  cases/
    01_box/box.geo               # reference gmsh script (retained)
    02_cylinder/cylinder.geo
    03_box_tunnel/box_tunnel.geo
    04_hollow_box/hollow_box.geo
    05_sphere/sphere.geo
    06_enclosure_coax/enclosure_coax_junction.geo   # ASCENT acceptance fixture
```

Each case folder keeps its reference `.geo` script. Generated `gmsh_out.msh` files
are git-ignored.

## Run

Tessella supports **Julia 1.12.x and 1.13.x** (`Project.toml`); the recorded
qualification runtimes are 1.12.7 and 1.13.1. From the repository root, use the
Juliaup selectors below, or substitute the corresponding installed executable:

```sh
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/run_all.jl
julia +1.13 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/run_all.jl
```

Follow the current handoff queue before launching: adopt and drain the existing
owned reader, keep one facade/package/aggregate reader, and require fresh physical
and available commit memory above 30 GiB. Bind the actual package path, exact
inputs, executable and driver; capture the actual exit, drain the process and
verify post-run hashes. Use normal `-O2` and default inlining. Debug runs with
`--compile=min` or disabled compiled modules are not release qualification.
Keep each resource test's original thread settings, including explicit four-thread
children; the aggregate already launches its Windows power probe with four threads.

Examples in the subdirectory READMEs use `+1.12`; repeat with `+1.13` unless both
executables are listed. The `sh` environment assignment and backslash continuation
syntax are for POSIX shells. In PowerShell set `$env:GMSH_JULIA_API` separately
and place the Julia invocation on one line. Put the pinned Gmsh CLI on `PATH`
and select its matching binding/DLL; individual drivers also accept their
documented `GMSH_EXECUTABLE` setting.

The aggregate gate requires the Gmsh 4.15.2 CLI and matching Julia API. It launches
the size-field, constant-range, uniform-refinement, whole-cache affine,
bulk/derived mesh-query, automatic/manual global edge/face topology, simplex
point-location, linear-simplex Jacobian/reference-map, fixed-family reference-quadrature,
actual- and explicit-order nodal and simplex hierarchical basis/orientation/key,
fixed-element catalog,
quadratic-tetrahedron, four-sided transfinite, straight transfinite
curve-law, three-sided transfinite, recombined-quadrangle, affine
transfinite-volume, five-face-prism, recombined-hexahedron, and `.geo`
meshing-constraint/Delete-lifecycle, mixed cache/query/refinement, 0D/1D generation,
transfinite QuadTri, AddVerts and NoNew differentials as
mandatory bounds-checked children. It also runs focused box, square, cone,
cylinder, Boolean snapshot/Delete lifecycle, NURBS/IGES, classified Point/Line-In-Surface and
Surface-In-Volume projection with nested sheet constraints and a hole, native `.geo`,
bounded expression- and numeric-list-backed geometry and entity lists,
point-local and explicit-topology `.geo` mesh-size constraints, API updates, and
spatial surface grading, bounded geometry and global Physical tag allocators,
the Physical-group API lifecycle, explicit model-topology, entity-identity,
entity-removal, spatial-query, native entity metadata including plane properties,
Point/Line/Plane evaluation and surface-reparametrization API differentials, and
entity visibility/color, Point-coordinate, and model-attribute state differentials,
and factory-aware `SetMaxTag` for tracked explicit and primitive topology,
projected MSH2/MSH4 single-/two-direction periodic surfaces, compact periodic node
pairs, embedded periodic curves, reusable-master, chained, and expression-backed
curve graphs, planar periodic boundaries of an explicit volume, and 2-D
boundary-layer parity cases.
The NURBS child both imports Gmsh-generated IGES and has Gmsh import and mesh
Tessella-generated type 126/128/144 records. Missing or wrong-version Gmsh,
failed probes, and parity mismatches make the aggregate command fail. Mesh-case
results print to the terminal and to `validation/REPORT.md`. The aggregate
contains selected probes; it does not replace every planned public body, family
resource gate or the complete package tests listed in the handoff.

## What each case checks

| Case | Geometry | Oracle | Point |
|------|----------|--------|-------|
| 01_box | axis-aligned box | V = 2 exact | both meshers must conform to a flat solid |
| 02_cylinder | solid cylinder | V = πR²H | curved-surface fidelity trade-off (Tessella inscribed N-gon vs gmsh true circle) |
| 03_box_tunnel | box with a through-tunnel | V = 24 exact | genus-1 flat solid |
| 04_hollow_box | box minus interior cavity | V = 35 exact | Boolean-difference (CSG) solid |
| 05_sphere | ball | V = 4/3·πR³ | curved-surface fidelity |
| 06_enclosure_coax | ASCENT coax feed-through | volumes non-empty | historical acceptance capture: tested gmsh 4.13/4.15 runs left air/case/pin volumes **empty**; literal CAD parity remains separate from the native reconstruction |

Flat solids give an *exact* analytic volume, so a passing row supplies a volume
cross-check. It must be accompanied by conformity, cell-map and identity audits.
Curved solids compare each represented surface model with the analytic curved
volume; neither a small volume error nor a successful report proves exact CAD parity.

## Adding a case

Drop a `.geo` in a new `cases/NN_name/` folder and add a matching Tessella surface
builder + analytic volume to the `cases` list in `run_all.jl`.

## Extending to other tools

The convention is tool-agnostic: add a `run_<tool>` helper in `support/common.jl` (mirroring
`run_gmsh`) that meshes to a Tessella-readable format, and a column in the report.
Retain each tool's input script alongside the case.
