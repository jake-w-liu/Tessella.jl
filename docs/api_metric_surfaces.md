# Explicit tensor Surface generation

The explicit `API.mesh.generate_metric_surface(tag; max_metric_edge=1.0,
max_nodes=10_000_000, max_triangles=20_000_000)` endpoint consumes the tensor
background selected with `API.mesh.field.set_as_background_mesh`. A
`MathEvalAniso` field can vary spatially; the selected background's existing
model size intersection, size bounds and factor are applied once.

This first scope requires one linear Plane Surface and no Regions. Boundary
Curves must be straight, and their actual retained Line2 cells must already
satisfy the requested metric bound. Generate a sufficiently fine lower mesh
first, for example using transfinite boundary Curve counts and
`API.mesh.generate(1)`. The endpoint can also construct the lower products as
part of its detached transaction. Interior edges are bisected conformingly.
Every final edge is checked using its actual displacement and the full tensor
at both endpoints and the midpoint. This certifies those three samples; it does
not certify arbitrary unsampled points or a minimum triangle angle.

Node and element UIDs, classifications, Curve cells and Surface constructor
metadata are published together only after production and ownership checks.
`Mesh.Renumber=0` preserves retained lower labels; enabling it requests the
usual dense renumbering. `Mesh.SaveAll` retains the ordinary publication
semantics. Returned meshes are detached. Failed production leaves the held
mesh and UID high-water marks unchanged.

Unsupported contexts include transfinite or recombined Surfaces, extrusion,
periodicity, embeddings, explicit Surface smoothing, higher element orders,
`MeshOnlyEmpty`, `MeshOnlyVisible`, queued homology, and higher-order held input
cells. Resource limits apply
to the actual Surface product, including its retained boundary rows.

This is a Tessella tensor meshing endpoint. Ordinary `API.mesh.generate` keeps
its existing algorithm dispatch. Native BAMG algorithm 7 remains guarded; this
endpoint does not claim that algorithm or its ordering. The
[official Gmsh tutorial t17](https://gmsh.info/doc/texinfo/#t17) uses BAMG for
anisotropic background meshing. Tiny opaque 4.15.2 observations show that
algorithm 6 can retain identical meshes when tensor directions change, so its
scalar behavior is not treated as a confirmed parity defect.

The opt-in [public functional fixture](../test/interfaces/api_metric_surface_functional_test.jl)
checks independent metric bounds, direction changes, lower UID ownership,
background factor/cap behavior, and failure atomicity. No default unstructured
generation workload is added.

The final47 checks passed onJulia1.13.1 with actual exit0 in49.80 s including
loading. Earlier failures caught retained Line cell UID renewal and regeneration
constructor rows; both were fixed in the tensor-only transport. Their original
receipts remain preserved.
