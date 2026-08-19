# Tessella.jl

Tessella is a Julia-native tetrahedral mesh generator built around exact geometric
predicates, certified topology, conforming boundary recovery, and explicit failure
diagnostics. It was created for the ASCENT electromagnetics workflow after gmsh
4.13.1 and 4.15.2-git left the enclosure/coax acceptance geometry with zero volume
elements.

The original simplex-mesher roadmap is complete through Stage 6. Development has now
expanded toward independent Gmsh 4.15.2 feature and behavioral parity, with
ASCENT-relevant meshing capabilities implemented first. That parity target is **not
complete**. The current bounds-checked suite passes 152,789/152,789 assertions; the
live implementation and verification record is [`STATUS.md`](STATUS.md).
The separate external ASCENT solve campaign is recorded in [`ASCENT.md`](ASCENT.md).

```julia
using Tessella
using Tessella.IO: write_msh

m = mesh_volume(surface)                 # closed triangle surface → validated tets
ms = mesh_sized(surface; hmax=0.5)       # certified maximum-edge size bound
near = DistanceField(surface)
field = ThresholdField(near; dist_min=0.0, dist_max=0.2,
                       size_min=0.02, size_max=0.5)
mf = mesh_sized(surface; field=field)    # spatially graded volume mesh
write_msh("mesh.msh", ms; version=4.1)   # solver-consumable gmsh MSH
```

## Implemented scope

- adaptive exact 2-D/3-D predicates with one consistent Simulation-of-Simplicity
  convention and exact-rational test oracles;
- compact simplex meshes with finite-input, cell, tag, manifold-link, and quality
  validation;
- 2-D Delaunay/CDT and quality refinement, graded curves, and planar, cylindrical,
  and parametric surface meshing;
- Float64 and exact-coordinate 3-D Delaunay kernels, constrained boundary recovery,
  multi-region partitions, and explicit recovery blockers;
- uniform and field-driven sizing for boxes, cylinders, extrusions, and arbitrary
  closed faceted surfaces, including Gmsh-compatible `Distance`, `Threshold`, `Box`,
  `Ball`, finite `Cylinder`, `Frustum`, `Min`, `Max`, and final size-bound semantics;
  discrete distance queries use a deterministic AABB hierarchy, and
  lower-dimensional cells and tags are preserved through volume refinement;
- quality reporting, flips, Laplacian/ODT/targeted sliver smoothing, healing
  diagnostics, native primitives, analytical surfaces, imprints, and mesh Boolean CSG;
- globally certified quadratic tetrahedra, plus strict and atomic MSH v2.2/v4.1 and
  STL I/O.

Mixed element families, the remaining Gmsh field types and algorithms, complete CAD/
BREP and scripting support, broad file/API compatibility, GUI, and post-processing are
still pending parity work. See [`PLAN.md`](PLAN.md) rather than treating the completed
Stage 0–6 baseline as Gmsh completeness. ASCENT remains Tessella's primary solver
consumer.

## Verification

```sh
julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'
julia --project --check-bounds=yes validation/run_all.jl
```

The first command runs the full CRC suite. The second compares analytic volumes and
quality with the installed gmsh and reproduces the enclosure/coax failure. See
[`DEVELOPMENT.md`](DEVELOPMENT.md) for the correctness, robustness, and completeness
rules.

## Repository map

- [`src/`](src/) — package implementation
- [`test/`](test/) — package regression and independent-oracle tests
- [`validation/`](validation/) — external-tool and solver-facing validation fixtures
- [`PLAN.md`](PLAN.md) — architecture and scope boundary
- [`STATUS.md`](STATUS.md) — current package verification record
- [`ASCENT.md`](ASCENT.md) — external solver integration and HFSS campaign

## License

No license has been declared by the repository owner.
