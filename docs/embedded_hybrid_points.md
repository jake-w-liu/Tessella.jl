# Embedded Point in a hybrid Region

One ordinary explicit Region with a certified convex planar Tri3/Quad4 shell
can consume one strictly interior embedded Point. The shell must contain at
least one quadrangle so generation uses the Tet4/Pyramid5 hybrid constructor.

After declaring the Region geometry, create and embed the Point through the
public API:

```julia
API.model.add_point(0.25, 0.3, 0.4; tag=9, meshSize=1.0)
API.mesh.embed(0, [9], 3, 1)
API.mesh.generate(3)
node_tags, coordinates, parameters = API.mesh.get_nodes(0, 9)
```

The constructor consumes the Point's actual generated Point15 cell and node
UID. During Volume generation, it retains that UID, Point ownership and exact
coordinates through the internal refinement, smoothing and optimization stages. The Point must be referenced by volume
cells with a closed connected interior Tet star. MSH4 export retains both the
Point-owned node and its incidence in the volume mesh.

Outside or boundary Points and missing Point15 provenance raise an error
without publishing a replacement mesh. Multiple Points, embedded Curves or
Surfaces, and embedded-Point multi-Region batches remain guarded. Existing
exterior identity and signed-volume certificates also apply to this context.

The opt-in regression is
[`api_hybrid_embedded_point_functional_test.jl`](../test/interfaces/api_hybrid_embedded_point_functional_test.jl).
