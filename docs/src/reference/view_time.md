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

See the [view option documentation](https://gmsh.info/doc/texinfo/gmsh.html#Post_002dprocessing-view-options)
for Gmsh option meanings. Display time selection is separate from PostView
mesh-size field selection; changing a display step does not select a new
background-field step.

[List/model representation replacement](view_representations.md) now clears
inactive data and starts the new representation's step/time history while
retaining display state. Same-list updates preserve the time behavior above.
