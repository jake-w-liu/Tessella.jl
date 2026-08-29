# Validation — Tessella vs external mesh tools

Independent cross-validation of Tessella against established open-source meshers
(currently **gmsh**). The convention: for every geometry, mesh the *same domain*
with both tools and compare against a shared oracle — the **analytic volume** — plus
element quality and wall-clock time. This checks that Tessella is not merely
self-consistent but agrees with an independent implementation, and quantifies where
it is more accurate or faster.

## Layout

```
validation/
  run_all.jl             # driver: runs every case, writes REPORT.md
  REPORT.md              # generated results table (git-ignored until you run it)
  support/
    common.jl            # helpers: gmsh runner, tet metrics, comparison
  size_fields/
    differential.jl      # required Gmsh 4.15.2 field differential
    STATUS.md            # exact coverage and explicit non-claims
  geo_ranges/
    differential.jl      # required bit-exact Gmsh 4.15.2 constant-range differential
  uniform_refine/
    differential.jl      # required Gmsh 4.15.2 linear-simplex template differential
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
    model_entity_metadata.jl # native entity types and partition ownership
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

Use the supported Julia 1.12.x runtime:

```sh
julia --project=. --check-bounds=yes validation/run_all.jl
```

The aggregate gate requires the Gmsh 4.15.2 CLI and matching Julia API. It launches
the size-field, constant-range, uniform-refinement, quadratic-tetrahedron,
four-sided transfinite, straight transfinite
curve-law, three-sided transfinite, recombined-quadrangle, affine
transfinite-volume, five-face-prism, and recombined-hexahedron differentials as
mandatory bounds-checked children. It also runs focused box, square, cone,
cylinder, Boolean snapshot/Delete lifecycle, NURBS/IGES, classified Point/Line-In-Surface and
Surface-In-Volume projection with nested sheet constraints and a hole, native `.geo`,
bounded expression- and numeric-list-backed geometry and entity lists,
point-local and explicit-topology `.geo` mesh-size constraints, API updates, and
spatial surface grading, bounded geometry and global Physical tag allocators,
the Physical-group API lifecycle, explicit model-topology, entity-identity,
entity-removal, spatial-query, and native entity-metadata API differentials, and
factory-aware `SetMaxTag` for tracked explicit and primitive topology,
projected MSH2/MSH4 single-/two-direction periodic surfaces, compact periodic node
pairs, embedded periodic curves, reusable-master, chained, and expression-backed
curve graphs, planar periodic boundaries of an explicit volume, and 2-D
boundary-layer parity cases.
The NURBS child both imports Gmsh-generated IGES and has Gmsh import and mesh
Tessella-generated type 126/128/144 records. Missing or wrong-version Gmsh,
failed probes, and parity mismatches make the aggregate command fail. Mesh-case
results print to the terminal and to `validation/REPORT.md`.

## What each case checks

| Case | Geometry | Oracle | Point |
|------|----------|--------|-------|
| 01_box | axis-aligned box | V = 2 exact | both meshers must conform to a flat solid |
| 02_cylinder | solid cylinder | V = πR²H | curved-surface fidelity trade-off (Tessella inscribed N-gon vs gmsh true circle) |
| 03_box_tunnel | box with a through-tunnel | V = 24 exact | genus-1 flat solid |
| 04_hollow_box | box minus interior cavity | V = 35 exact | Boolean-difference (CSG) solid |
| 05_sphere | ball | V = 4/3·πR³ | curved-surface fidelity |
| 06_enclosure_coax | ASCENT coax feed-through | volumes non-empty | the acceptance case: gmsh 4.13/4.15 leave the air/case/pin volumes **empty**; Tessella's native pipeline targets filling them |

Flat solids give an *exact* analytic volume, so a passing row is a hard correctness
cross-check. Curved solids expose the geometric-fidelity trade-off honestly (each
tool is exact for its own surface model).

## Adding a case

Drop a `.geo` in a new `cases/NN_name/` folder and add a matching Tessella surface
builder + analytic volume to the `cases` list in `run_all.jl`.

## Extending to other tools

The convention is tool-agnostic: add a `run_<tool>` helper in `support/common.jl` (mirroring
`run_gmsh`) that meshes to a Tessella-readable format, and a column in the report.
Retain each tool's input script alongside the case.
