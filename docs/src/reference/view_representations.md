# List and model view representation changes

A session view has one active data representation. A successful
`API.view.add_model_data` call replaces existing list records with the admitted
model data. A successful `API.view.add_list_data` call on a model view replaces
its model rows and bindings with the admitted list records. The view keeps its
tag, name, numeric display options, and displayed `TimeStep`/`Time` until the
caller explicitly selects a step again.

The model-to-list transition admits the 24 plain linear list families
`SP` through `TY`, including an empty recognized list. It validates complete
coordinate/value layouts before removing model data. Unsupported high-order
or interpolation layouts remain outside this transition contract. Existing
list-to-list admission, physical-time history and export validation keep their
existing behavior.

`NbTimeStep` follows the new representation: model data retain their allocated
sparse prefix, ordinary list records determine their actual step count, and an
empty recognized replacement list resets the count to zero. Returning from
model data to a list starts a new implicit time history; previously imported
physical list times do not leak through the model representation. Updating a
list while it remains a list retains its physical-time history.

The inactive query reports a checked `ArgumentError`: `get_list_data` rejects
an active model view and `get_model_data` rejects an active list view. The
public query tuple shapes remain unchanged. Model/list export succeeds after
a valid transition; the writer still rejects internally corrupted mixed
storage. Removing one view's bindings does not mutate any held model slot or
another view's data.

Transition rejection preserves the previous representation, metadata and
model state. This is stronger than the observed Gmsh 4.15.2 unknown-list-type
failure, which removed model data before reporting an error. It does not claim
native partial-failure or undefined C++ exception parity.

These representation changes preserve PostView field selection, interpolation-list
import, generated meshes or broader visualization controls.
