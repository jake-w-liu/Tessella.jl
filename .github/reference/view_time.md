# Session view time metadata

`API.view.get_number(tag,"NbTimeStep")` reports actual data steps. Setter calls
for NbTimeStep are ignored. Plain list families use their minimum nonempty
step count; clearing all list families retains the last computed count, while
a fresh empty view reports zero. Model data use their allocated step prefix,
including actual empty gaps with default time zero. A data-kind change starts
the new actual prefix and no longer admits steps left over from the old kind.

TimeStep setters truncate finite fractions to a signed 32-bit integer. Negative
indices select the last legal step; indices at or above the count select zero.
For a fresh empty view, a negative selection yields -1 and a nonnegative one
yields zero. Time updates to the selected physical time. Nonnegative Time
setters select the nearest physical time, choosing the first index on ties;
negative Time requests leave display state unchanged. Other numeric options
retain their existing storage behavior.

Data updates preserve the current displayed TimeStep and Time until explicit
selection, even when replacing data makes the selection temporarily invalid.
Imported list TIME history is retained when the data count shrinks; growing
beyond existing history appends implicit index times. The lossless POS writer
keeps its strict TIME/value count certificate instead of emitting mismatched
payloads. Sparse model nearest-time lookup uses populated contexts and the
earliest zero-time gap without allocating arrays proportional to the maximum
step identifier. Geometry, UIDs, held model slots and data values are unchanged
by numeric selection.

The behavior is grounded in retained public Gmsh 4.15.2 controls and the
[official view option documentation](https://gmsh.info/doc/texinfo/gmsh.html#Post_002dprocessing-view-options).
The public SDK rejects a PostView field option called TimeStep and exposes no
direct field evaluator; this increment makes no field selection/parity change.
The existing indexed PostView field path remains separate. A retained public
one-dimensional control applied a two-step scalar Line view with sizes 0.5 and
0.25 as a background. Display selections zero and one both produced two Line
elements and the same three node coordinates. This supports keeping the
current field selection for that bounded control; it is not a general field
time-parity claim. Its opaque result is
`C:/tmp/tessella_postview_selection_20261011/native_postview_time_selection.json`.

Qualification is bounded and uses only ordinary list files and one caller-owned
Point cell. It covers truncation/wrapping, physical-time and empty metadata,
list updates, sparse gaps, overwritten model times, kind-reset bounds, finite
and integer-range failure atomicity, and retained UID/data/geometry. It includes
one million-step sparse range as a storage boundary check without constructing
a million rows. Two cheap default assertions query one two-step point list.
Actual runtime results will be recorded in the final packet.

The final reader passed 37/37 (35 functional checks plus two smoke assertions),
Julia 1.12.7 with `--compile=min -O0 --check-bounds=yes` and one thread, actual
exit zero, warm wall 1.0141333 seconds. Failed v1 is retained: the test used the
SDK's five-item model-data return shape instead of Tessella's existing
`(types,tags,flat_values)` query contract. Corrected expectations preserve all
metadata/UID/geometry checks; implementation Source is unchanged.

[List/model representation replacement](view_representations.md) now clears
inactive data and starts the new representation's step/time history while
retaining display state. Same-list updates preserve the qualified behavior above.
