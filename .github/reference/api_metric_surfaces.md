# Explicit tensor Surface endpoint

`API.mesh.generate_metric_surface(tag; max_metric_edge=1.0,
max_nodes=10_000_000, max_triangles=20_000_000)` consumes the tensor background
selected with `API.mesh.field.set_as_background_mesh`. The existing session
size intersection, bounds and factor are applied once.

This endpoint requires the model's single planar Surface, no Regions, and
linear input and output elements. Its lower Curve cells must already satisfy
the requested sampled metric bound. Generate sufficiently fine lower products
first, for example with transfinite Curve counts and `API.mesh.generate(1)`;
the endpoint can also construct them within its detached transaction.

The shared [tensor producer](metric_model_surfaces.md) supports planar holes,
geometric curved boundaries, embedded Points and Curves, and recombination
algorithms 0–3. Lower UIDs, coordinates, classifications and literal Curve
incidence remain owned source data. The refined triangular product is checked
against the full tensor at each edge's endpoints and midpoint. Recombination
does not imply the same bound for every final quad edge.

Node and element UIDs, classifications, Curve cells and Surface constructors
are published together after production and ownership checks. `Mesh.Renumber=0`
preserves retained lower labels; enabling it requests dense renumbering.
`Mesh.SaveAll` retains the ordinary publication semantics. Returned meshes are
detached. A failed transaction leaves the held mesh and UID high-water marks
unchanged.

The endpoint's admission still rejects transfinite/extruded Surface
construction, periodic slaves, higher element orders, `MeshOnlyEmpty`,
`MeshOnlyVisible`, queued homology, and higher-order held input cells. Resource
limits count actual retained and newly born Surface rows.

Ordinary `API.mesh.generate(2)` supports algorithm 7 through the native Julia
metric producer and has its own broader dispatch. It does not require this
single-Surface endpoint. Neither path depends on an external mesher, and neither
claims identical BAMG ordering.

Public regressions:
[`api_metric_surface_functional_test.jl`](../../test/interfaces/api_metric_surface_functional_test.jl)
and [`api_metric_surface_embedded_test.jl`](../../test/interfaces/api_metric_surface_embedded_test.jl).
See [development evidence](../development/README.md) for qualified runtime scopes.
