# Model views as mesh-size fields

A `PostView` background field accepts model-based `NodeData`, `ElementData`,
and `ElementNodeData`, as well as list views. Keep the source mesh in its own
model while creating the target geometry in another model.

After adding the source mesh and its values:

```julia
view = Tessella.API.view.add("background sizes")
Tessella.API.view.add_model_data(view, 0, "source model", "NodeData",
    node_uids, [[0.25] for _ in node_uids])

# Select/create the target model before adding its background field.
field = Tessella.API.mesh.field.add("PostView")
Tessella.API.mesh.field.set_number(field, "ViewTag", view)
Tessella.API.mesh.field.set_as_background_mesh(field)
```

Sizing samples step 0, independently of the view's displayed time step. An
empty step 0 contributes no size constraint. Scalar data supplies a size;
vectors use their norm; nine-component data supplies a metric tensor through
the existing tensor field path. The mesher still validates positive sizes and
positive-definite metrics at consumption.

Fields retain source model identities and an owned geometry/value snapshot.
Aliases share their source data while field construction detaches the records.
Sparse cells without complete nodal values contribute no record; unknown UIDs
and malformed element-node rows are explicit errors.

Full element geometry and interpolation are retained for supported high-order
families, including complete pyramids and warped quadrangles. Curved geometry
search boxes enclose the parameterized maps rather than just their nodes.
Overlapping cells use the same dimensional/family precedence for linear and
high-order data. Unsupported element bases keep their explicit diagnostics.
