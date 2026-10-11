# Constant-tensor box generation

The positional `Metric3` overload constructs a structured tetrahedral box using
three directional interval counts:

```julia
metric = Metric3(9.0, 1.0, 1.0)
mesh = mesh_box(0, 1, 0, 1, 0, 1, metric; max_metric_edge=1.0)
```

The tensor measures an edge displacement `d` by `sqrt(d' * M * d)`. Every
actual stored edge meets `max_metric_edge`. Diagonal metrics produce separate
axis counts; offdiagonal entries contribute to a conservative bound for all
seven Kuhn edge displacement classes. The grid remains aligned with the box
axes; it does not rotate to the metric eigenvectors.

Literal endpoint planes remain exact. Rounded intervals are checked before the
existing six-Tet Kuhn subdivision is mapped to physical coordinates. The result
is validated for positive, nondegenerate tetrahedra and manifold topology.
`max_nodes` and `max_tets` bound allocation; unrepresentable coordinates and
resource/count errors reject explicitly.

The positional `mesh_box_regions(regions, metric; max_metric_edge=...)`
overload uses one shared grid for all `BoxRegion` entries. Every literal Region
face becomes a grid plane. Adjacent materials share node and face rows;
`tet_tag` retains their tags. The last containing Region takes priority,
and `tag == 0` carves a void, matching the physical `hmax` overload.

```julia
regions = [BoxRegion(0, 1, 0, 1, 0, 1, 3),
           BoxRegion(1, 2, 0, 1, 0, 1, 7)]
mesh = mesh_box_regions(regions, metric)
```

Resource limits apply to the entire shared bounding grid before compaction,
including outside and void cells. Widely separated boxes can therefore exceed
a limit even when the final mesh would be small. Input Regions remain unchanged;
signed zero is treated as one geometric plane.

The six-argument `mesh_box(...; hmax=...)` retains its physical Euclidean bound.
This overload takes an explicit constant tensor, not a sampled field assumed to
be constant. Existing Surface meshing and Volume refinement already consume
directional field metrics. General native metric-mesher parity remains open.
