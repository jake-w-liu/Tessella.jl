# List-based view exports

`API.view.write(tag, path)` accepts `.pos`, `.txt`, `.csv` and `.vtk` for
list-based views. List data store all x coordinates, then all y and all z
coordinates, followed by component/node values for every time step.

TXT writes one numeric row per element/time pair. Each corner contributes
`step time 0 element x y z value...`. Steps and element ordinals start at zero;
list times are the step indices. CSV separates the same fields with commas.
Coordinates and values use 17 significant digits. Scalar, vector and tensor
records retain 1, 3 and 9 components. TXT/CSV permit mixed component counts,
so row widths follow each element's shape and component count.

VTK writes an ASCII legacy unstructured grid. Every element retains its own
points: coincident coordinates can carry different values. Linear point,
line, triangle, quadrangle, tetrahedron, hexahedron, prism and pyramid cells
are supported. Global `TIME` field data records all step times; point arrays
`step_0`, `step_1`, ... retain every step. Scalars use `SCALARS`, vectors use
`VECTORS`, and tensors use nine-component `FIELD` arrays to preserve
nonsymmetric input. VTK requires one common component count and rejects mixed
scalar/vector/tensor storage. VTK append is rejected; TXT/CSV append adds rows.
`PostProcessing.Binary` produces an ASCII-only warning for these formats.
Payload and format validation occur before opening the destination.

Public Gmsh 4.15.2 observations write identical whitespace tables for `.txt`,
`.csv` and `.vtk`. Tessella follows the TXT record layout and deliberately
emits comma-separated CSV and valid VTK. Those two extension aliases therefore
do not have byte parity with the pinned SDK. Mixed list/model storage and
high-order list records remain explicit unsupported export contexts.

Model-based views use their own [CSV schema and VTK arrays](model_view_exports.md).

Format references: [Gmsh list-data tutorial](https://gmsh.info/doc/texinfo/gmsh.html#x3),
[VTK legacy format](https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html),
[VTK cell types](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkCellType.h).

[Session time metadata](view_time.md) derives ordinary list-step counts and
selects displayed physical times while retaining the actual data arrays.
Imported TIME history survives updates; the lossless writer rejects mismatched
TIME/value counts.
