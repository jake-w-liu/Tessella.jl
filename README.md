# Tessella.jl

Tessella is a **Julia-native tetrahedral mesh generator** built around exact
geometric predicates, certified topology, conforming boundary recovery, and
explicit failure diagnostics. It targets public Gmsh 4.15.2 behavior parity —
Gmsh is used only as a reference oracle; all meshing runs natively in Julia.

It was created for the ASCENT electromagnetics workflow after gmsh 4.13.1 and
4.15.2-git left the enclosure/coax acceptance geometry with zero volume
elements. The simplex-mesher roadmap is complete; full Gmsh parity is in
progress across the P1–P6 tracks — see
[DEVIN_HANDOFF.md](DEVIN_HANDOFF.md) for the open-work ledger.

## Quickstart

Requires **Julia 1.12.x or 1.13.x**.

```julia
using Tessella
```

### Command line

```sh
julia -e 'using Tessella; Tessella.CLI.main(ARGS)' file.geo -3 -o mesh.msh
```

The CLI mirrors gmsh's batch interface: `file.geo -1|-2|-3`, `-o out.msh`,
`-format msh2|msh4`, `-bin`, `-clscale`/`-clmin`/`-clmax`, `-order`,
`-scale`, `-part N`, `-refine`, `-save_all`, `-log`, `-cpu`, `-statreport`
and more. `-help` and `-version` work standalone; `file.msh -check` validates
a mesh without generating. See [docs/cli_controls.md](docs/cli_controls.md).

### `.geo` scripts

Standard gmsh `.geo` files execute natively — geometry kernels
(`Point`/`Line`/`Circle`/`Ellipse`/`Spline`/`BSpline`/`Nurbs`,
`Extrude`/`Revolve`/`Fillet`/`Chamfer`, `Boolean*` CSG in both the built-in
and OpenCASCADE factories), `Physical` groups (including
`Physical GeoEntity{d}` forms), `Transfinite`/`Recombine`, `Field` size
specs, periodicity, `PartitionMesh`, `Plugin(...)` dispatch, and
`Save`/`Print` output.

### Public API

```julia
const A = Tessella.API
A.initialize()
occ = A.model.occ
occ.add_box(0, 0, 0, 1, 1, 1)
occ.synchronize()
A.mesh.generate(3)
A.write("box.msh")
A.finalize()
```

The API mirrors `gmsh/model`, `gmsh/model/occ`, `gmsh/model/geo`,
`gmsh/model/mesh`, `gmsh/mesh`, `gmsh/mesh/field`, `gmsh/option`, `gmsh/view`,
and `gmsh/plugin` surfaces, plus `API.open`/`API.merge`/`API.write` for
MSH/GEO/POS files. The oracle-registered plugin set is complete (CurvedBnd
distance fields, Bubbles, Crack, Homology, NearToFarField, SimplePartition,
and the full cut/transform/compute families).

## What Tessella does

- **Meshing**: Delaunay/frontal 2-D and 3-D with conforming boundary
  recovery, transfinite and structured grids, recombination to
  triangles/quads, uniform and barycentric refinement, high-order elements
  with whole-domain Jacobian certification, size-field evaluation
  (MathEval/PostView/octree `AutomaticMeshSizeField`), periodic and embedded
  constraints, partitioning with entity/ghost/topology layers.
- **Geometry**: built-in and OCC-kernel modeling, Booleans on solids and
  surfaces with generic curved-section tracing, STEP/IGES import (solids,
  surface sheets, mapped occurrences, assemblies, SI units), affine
  transforms, extrusion/lofts/fillets.
- **I/O**: MSH 2.2 + 4.1 (ASCII/binary, UInt64 tags, preserved sections),
  STL, UNV, Medit `.mesh`/`.meshb`, SU2, BDF, INP, NEU, P3D, MAIL, IR3,
  DIFF, VTK, OFF, OBJ, PLY2, TOCHNOG `.dat`, MATLAB `.m`, VRML, LS-DYNA
  `.key`, RADIOSS `.rad`, X3D, POS/CSV view data. Clear errors for
  dependency-gated formats (CGNS/MED/CDF).

Feature details and per-topic scope live in [`docs/`](docs/).

## Testing

```sh
julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

`Pkg.test()` selects minimal compilation automatically; for optimized-mode
runs set `TESSELLA_TEST_COMPILER=optimized` (see
[test/README.md](test/README.md)). The default suite is small real fixtures
and analytic invariants — no external tools. Extended sweeps, allocation
audits, benchmarks, and solver campaigns are opt-in:

```sh
julia --project=. --check-bounds=yes test/extended/runtests.jl
```

## Repository map

- [`src/`](src/) — `core`, `fields`, `geometry`, `meshing`, `structured`, `interfaces`
- [`test/`](test/) — matching suites plus `integration` end-to-end cases
- [`validation/`](validation/) — external-tool and solver-facing fixtures
- [`docs/`](docs/) — per-topic implementation scope
- [`PLAN.md`](PLAN.md) — architecture and scope boundary
- [`STATUS.md`](STATUS.md) / [`HANDOFF.md`](HANDOFF.md) — dated verification
  evidence and remaining work
- [`DEVIN_HANDOFF.md`](DEVIN_HANDOFF.md) — complete open-work checklist

## License

No license has been declared by the repository owner.
