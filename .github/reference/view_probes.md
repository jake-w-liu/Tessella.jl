# View probes

The canonical query accepts list views and model-based NodeData, ElementData
and ElementNodeData:

```julia
values, distance = Tessella.API.view.probe(view_tag, x, y, z;
    step=-1, numComp=-1, gradient=false, distanceMax=0.0, dim=-1)
```

`step=-1` returns all allocated time steps, concatenated in step order. An
explicit step samples that frame's own geometry. All-step model queries use the
first frame's material geometry and bind later values by actual node/element
UIDs. Missing rows of a populated later frame retain preceding material values;
an allocated empty step emits zeros without erasing that history. Aliases retain
their source identities when another model is active.

`numComp` selects scalar (1), vector (3), or tensor (9) records. Automatic
selection prefers an exact scalar match, then vector, then tensor. `dim` filters
by element dimension. Full geometry and interpolation bases are retained for
supported high-order cells. Gradients are physical, grouped by component with
three coordinates per component; a curve or surface returns the minimum-norm
tangential gradient.

Exact element matches take precedence. A positive `distanceMax` enables
nearest-node sampling within that inclusive radius; a negative value enables it
without a radius limit. An ordinary automatic miss still reports the nearest
distance. A class-filtered miss with no radius request returns empty values and
distance -1. Supplying `xElemCoord`, `yElemCoord`, and `zElemCoord` restricts the
search to the ordered element coordinates supplied by the caller.

The legacy `API.view.probe(model_name, view_tag, x, y, z)` overload samples
step 0 of the named model and returns its component vector. It retains arbitrary
legacy component widths and reports an explicit error when no element contains
the point. The canonical overload returns `(values, distance)`.

The public Probe plugin uses the same all-step material query semantics. Mesh
sizing has a separate step-0 contract described in
[model view fields](model_view_fields.md).

Regression:
[`api_view_probe_basis_test.jl`](../../test/interfaces/api_view_probe_basis_test.jl).
Qualified supported-runtime and optimized checks are recorded in the
[development evidence](../development/README.md).
