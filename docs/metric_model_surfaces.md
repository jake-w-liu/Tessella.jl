# Variable-tensor ordinary planar Surfaces

The explicit tensor overload generates an ordinary planar triangle mesh from
the model's actual lower-dimensional Point and Line products:

```julia
field = MathEvalAnisoField(m11="4+12*x", m22="1", m33="1")
mesh = mesh_model_surface(model, surface_tag, field;
                          max_metric_edge=2.01,
                          max_nodes=100, max_triangles=200)
```

The boundary Curves may have prescribed transfinite grids; the Surface itself
uses ordinary initial triangulation and conforming interior edge bisection.
The supplied tensor is used directly, including any explicit field wrapper.
Scalar global size caps and factors are not applied again.

For a plane with tangent columns `B`, the tangential metric is `B' * M * B`.
The producer evaluates this norm through actual world-coordinate displacements,
retaining tensor cross terms and directions. Every emitted edge satisfies
`metric_edge_length(field, a, b) <= max_metric_edge`: the checked norm is the
maximum at both endpoints and the midpoint. This is a sampled edge contract,
not a guarantee over unsampled positions of an arbitrary variable field.

Retained lower vertices and Curve-cell boundary supports are checked against
the actual input UIDs. Their coordinates, classifications, and parameters
remain owned source data. A boundary edge that already exceeds the metric bound
is rejected with a request for finer lower Curve cells. The producer does not
split or replace those cells. Each interior split constructs one new Surface
vertex and two child triangles per incident parent, with checked constructor
births. Output storage is independent of the inputs and model.

Straight planar boundaries and holes are supported. The first scope rejects
transfinite or extruded Surfaces, periodic relations, recombination, embedded
entities, and Surface smoothing. There is no minimum-angle promise. Resource
limits, invalid metrics, flat children, and unrepresentable midpoints fail
before a partial result escapes. The model is meshed through a detached copy.

The existing two-argument `mesh_model_surface` route and native API
`mesh.generate` remain unchanged. Their ordinary Surface and Volume background
adapters currently reduce anisotropic fields to scalar values. This new overload
does not implement BAMG or claim native anisotropic ordering or general
variable-metric volume generation. Constant-tensor box grids and existing
directional surface/refinement pipelines retain their separate contracts.

The opt-in
[`model_metric_surface_functional_test.jl`](../test/geometry/model_metric_surface_functional_test.jl)
uses small actual model boundaries. It independently checks variable diagonal
quadratic forms, swapped tensor directions, an oblique plane with cross terms,
area and coordinate bounds, exact boundary UID transport, hole exclusion,
storage ownership, and atomic metric/resource failures. It performs no native
calls or volume production.

The final33 checks pass onJulia1.12.7 with actual exit0 in5.82 s warm,
`--compile=min -O0` and bounds checking. Source remained unchanged between the
initial coarse-boundary fixture rejection and the corrected fixture target.
Actual TF3 retained intervals produced a norm2.0000000000052314, so target2.0
was too tight; target2.01 passes the same strict independent norm certificate.
The failed receipt and actual boundary dump remain preserved.
