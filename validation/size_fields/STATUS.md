# Size-field differential status

Target: Gmsh 4.15.2. The runner requires both the Gmsh CLI and Julia API; a
missing runtime, wrong version, failed probe, or parity mismatch exits nonzero.
Intentional context gaps are printed as `CONTEXT_SKIP` and counted separately.

## Coverage classification

| Classification | Covered behavior | Meaning |
|---|---|---|
| Direct pointwise | `MathEval` (including an `F<n>` reference), `Gradient`, `Laplacian`, `Mean`, `Curvature`, `MaxEigenHessian`, `LonLat`, `Param`, `Structured`, entity-aware `Restrict`/`Constant`, scalar-point `PostView`, `ExternalProcess` | Gmsh values come from the official `MeshSizeFieldView` plugin on model-backed nodes and are compared with Tessella at the same coordinates/context. |
| Direct sampled/context | `Octree`, `Extend`, scalar `BoundaryLayer`, scalar `AttractorAnisoCurve` | The same sampled representation or entity context is supplied to both implementations. This does not claim broader CAD/topology behavior. |
| Scalar operator only | `MathEvalAniso`, `MinAniso`, `IntersectAniso` | Gmsh's public plugin exposes the scalar operator, not the underlying `SMetric3`. |
| Mesh-observed approximate | `Distance→Threshold`, `Box`, `Ball`, finite `Cylinder`, `Frustum` | Node counts and regional spacings are bounded; exact nodes are not expected from different curve discretizers. |
| Implemented, not directly covered | first-order scalar line/triangle/tetrahedron `PostView` interpolation | Unit-oracle tested, but this runner directly probes only scalar-point views; no Gmsh differential claim is made here. |
| Explicit context skip/non-claim | anisotropic metric tensors/metric-driven meshing, boundary-layer element topology, `AutomaticMeshSizeField`, and non-simplex/higher-order/vector/tensor `PostView` | No equivalent public oracle/shared model state exists for the first three; the listed broader `PostView` data is unsupported. No parity claim is made. |

## Known Gmsh probe constraints

- VERIFIED probe constraint: an `Octree` view probe before any Gmsh mesh pass
  crashed in `OctreeField::Cell::evaluate`. The runner performs a real mesh
  pass before evaluating that fixture.
- VERIFIED probe constraint: a `MathEval` field referring to another
  `MathEval` field did not return when driven by the view plugin. The runner
  exercises `F<n>` through a `Box` dependency instead.
- Gmsh 4.15.2 scalarizes the tested 1-D `MathEvalAniso` background and does not
  provide an honest directional 1-D metric comparison. It is reported as a
  context skip, not accepted with a loose tolerance.

## Last measured result

On 2026-08-21, with `/opt/homebrew/bin/gmsh` and the matching Homebrew Julia
API/library, this bounds-checked command exited 0:

```sh
GMSH_EXECUTABLE=/opt/homebrew/bin/gmsh \
GMSH_JULIA_API=/opt/homebrew/lib/gmsh.jl \
julia --project=. --startup-file=no --check-bounds=yes validation/size_fields/differential.jl
```

The terminal summary was:

```text
SIZE_FIELD_DIFFERENTIAL_OK gmsh=4.15.2 plugin_calls=23 direct_cases=23 direct_samples=63 mesh_cases=5 context_skips=5
```

The largest accepted direct discrepancy was `1.78e-15` for
`MaxEigenHessian`; all five mesh-observed cases satisfied their explicit count
and spacing bounds. The five `CONTEXT_SKIP` records above remain non-claims.

Aggregate failure propagation was also checked with an invalid explicit
runtime:

```sh
GMSH_EXECUTABLE=/definitely/missing/gmsh \
julia --project=. --startup-file=no --check-bounds=yes validation/run_all.jl
```

The child reported the invalid path and `validation/run_all.jl` exited 1 with a
`ProcessFailedException`, before entering its mesh-case/report-writing loop.
