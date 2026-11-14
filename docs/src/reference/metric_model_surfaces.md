# Tensor Surface generation

The explicit tensor overload builds a planar mesh from the model's actual
Point and Curve products:

```julia
field = MathEvalAnisoField(m11="4+12*x", m22="1", m33="1")
mesh = mesh_model_surface(model, surface_tag, field;
    max_metric_edge=2.01, max_nodes=100, max_triangles=200)
```

The supplied tensor, including an explicit wrapper, is consumed directly.
Global scalar size caps and factors are not applied again. Ordinary
`mesh_model_surface(model, tag; size_field=field)` and `API.mesh.generate(2)`
also consume tensors when algorithm 7 is selected; those paths apply the
session's directional size bounds, factor and anisotropy cap.

The native producer uses metric-Delaunay Steiner insertion, conforming edge
bisection, flips and vertex removal. It retains full world-coordinate tensor
cross terms, equivalent to the tangential pullback `B' * M * B`. Every edge of
the refined triangular product satisfies
`metric_edge_length(field, a, b) <= max_metric_edge`, evaluated at its endpoints
and midpoint. This is a sampled contract; it does not bound arbitrary unsampled
points of a variable field or promise a minimum angle.

Actual lower vertices retain their UIDs, coordinates, classifications and
parameters. Planar faces support holes, geometric curved boundaries, embedded
Points and Curves, and recombination algorithms 0–3. An embedded Curve retains
literal edges with two incident Surface cells. A supplied Point on a Curve
splits that actual chain before full-quad coarsening. Retained lower cells must
satisfy the requested metric bound; the producer rejects under-resolution
instead of silently replacing owned cells. The triangular edge bound above is
checked before recombination and is not a claim about every final quad edge.

The ordinary algorithm-7 dispatcher also has ruled and OCC parametric chart
producers.
The explicit tensor overload remains planar. Transfinite and extruded Surface
construction uses its own producer. Algorithm-7 ordering is native Julia and
does not claim identical BAMG chronology or matching tie order.

Periodic Surface masters publish the final refined/recombined constructor to
the ordinary owned copy producer. Per-Surface elliptic smoothing is retained as
a transfinite attribute, as in the [Gmsh manual](https://www.gmsh.info/doc/texinfo/gmsh.html#Structured-grids).
It does not reject ordinary algorithm-7 generation.

New vertices and cells carry constructor births. Resource limits, invalid
metrics and representability failures return no partial mesh or model mutation;
the explicit overload works on a detached model copy. See the
[session endpoint](api_metric_surfaces.md) for publication and additional
session restrictions.
