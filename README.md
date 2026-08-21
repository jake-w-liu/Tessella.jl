# Tessella.jl

Tessella is a Julia-native tetrahedral mesh generator built around exact geometric
predicates, certified topology, conforming boundary recovery, and explicit failure
diagnostics. It was created for the ASCENT electromagnetics workflow after gmsh
4.13.1 and 4.15.2-git left the enclosure/coax acceptance geometry with zero volume
elements.

The original simplex-mesher roadmap is complete through Stage 6. Development has now
expanded toward independent Gmsh 4.15.2 feature and behavioral parity, with
ASCENT-relevant meshing capabilities implemented first. That parity target is **not
complete**. The last stable bounds-checked package gate passed 161,183/161,183
assertions; newer isolated increments and their focused gates are recorded in
[`STATUS.md`](STATUS.md).
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
- uniform and field-driven 1-D, 2-D, surface, and 3-D sizing, including geometric and
  composite (`Distance`, `Threshold`, `Box`, `Ball`, finite `Cylinder`, `Frustum`,
  `Min`, `Max`, bounds), analytic/derivative (`MathEval`, `Gradient`, `Laplacian`,
  `Mean`, `Curvature`, `MaxEigenHessian`), coordinate/storage (`LonLat`, `Parametric`,
  `Structured`), entity-aware (`Restrict`, `Constant`), sampled/context (`Extend`,
  `Octree`, first-order scalar/vector `PostView` on every standard point/line/surface/
  volume list element), scalar `BoundaryLayer`,
  discrete `AutomaticMeshSizeField` analogue, and `ExternalProcess` fields;
- anisotropic metric primitives and `MathEvalAniso`, `MinAniso`, `IntersectAniso`, and
  `AttractorAnisoCurve`; discrete distance queries use a deterministic AABB hierarchy,
  and lower-dimensional cells and tags are preserved through volume refinement;
- quality reporting, flips, Laplacian/ODT/targeted sliver smoothing, healing
  diagnostics, native box/cylinder/cone/geodesic-sphere primitives, analytical
  surfaces, imprints, mesh Boolean CSG, and validated affine transformations of
  finalized simplex meshes with orientation and physical-tag preservation;
- globally certified quadratic tetrahedra, plus strict and atomic simplex MSH
  v2.2/v4.1 and STL I/O;
- a resource-bounded `.geo` scanner for finite arithmetic constants, pure numeric
  functions, prior scalar bindings, sizing options, and explicit field/physical tags;
- a 125-type fixed-node Gmsh catalog plus ten serializable cut/border/child/
  sub-element records, compact variable connectivity and parent/domain references,
  mixed entity/classification metadata, structural validation/CRC, and ASCII/binary
  MSH v2.2/v4.1 mixed-element I/O with opposite-endian decoding;
- deterministic, physical-tag-preserving surface triangle-to-quadrangle
  recombination with convexity, resource-growth, CRC, and Gmsh-load gates;
- one-level uniform linear-simplex refinement using Gmsh 4.15.2's ordered 2/4/8
  child templates, shared deterministic edge midpoints, tag preservation, resource
  bounds, and output validation;
- validated four-sided planar transfinite triangle patches using average-chord Coons
  interpolation, all four Gmsh diagonal arrangements, physical tags, bounded
  intersection auditing, and exact orientation postconditions;
- normalized straight-curve transfinite parameters for Gmsh's Progression/Power,
  Bump, and Beta laws, with signed orientation and Float64 representability gates;
- validated three-sided planar structured triangle patches using Gmsh's specific
  triangular interpolation, compact topology, and exact boundary postconditions;
- affine six-face transfinite volumes using Gmsh's unrecombined six-tetrahedron
  subdivision, with boundary tags and exact topology/orientation postconditions.

P1 through P4 remain **in progress**. Current non-claims include boundary-layer element
construction, the full Gmsh automatic-sizing pipeline, broader `PostView` data,
including high-order/custom interpolation, multiple time steps, materially warped
quadrangles, mixed component counts, and tensor-to-metric evaluation,
general entity/OCC/BREP/NURBS and full `.geo` execution (including loops, macros,
dynamic tags, option reads, ranges, CSG statements, and physical-group RHSs),
mixed-element generation or recombination beyond first-order surface triangle pairing,
non-affine CAD curve integration and FlexibleTransfinite/HWall laws,
quasi-transfinite or holed patches, transfinite
quadrangles, five-face or curved/warped transfinite volumes,
selective/high-order refinement, simplex-kernel integration,
MINI basis-selector tags 138/139 as mesh records, curved-cell Jacobian certification,
non-8-byte binary data, and ancillary-section preservation. MSH2 ASCII is the lossless
format for variable connectivity and parent/domain links; binary MSH2 and MSH4 have
explicitly narrower special-record contracts. Type 69 and some registered fixed tags
require Tessella-only output because Gmsh 4.15.2 cannot consume them safely, and
nonzero-physical special MSH4 requires compatible classification metadata for a
Gmsh-safe rewrite. Complete CAD/BREP, broad file/API compatibility, GUI, and
post-processing remain pending. See
[`PLAN.md`](PLAN.md) rather than treating the completed Stage 0–6 baseline as Gmsh
completeness. ASCENT remains Tessella's primary solver consumer.

## Verification

```sh
julia --project --check-bounds=yes -e 'using Pkg; Pkg.test()'
julia --project --check-bounds=yes validation/run_all.jl
```

The first command runs the full package/CRC suite. The second compares analytic
volumes and quality with Gmsh, reproduces the enclosure/coax failure, and runs the
Gmsh 4.15.2 size-field, uniform-refinement, transfinite-patch, and straight
transfinite curve-law, three-sided transfinite, and affine transfinite-volume
differentials as mandatory bounds-checked children. Missing or wrong-version Gmsh and
differential failures make the aggregate command fail. See
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
