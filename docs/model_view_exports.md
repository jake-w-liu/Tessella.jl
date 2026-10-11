# Model-based CSV and VTK exports

Model CSV writes a rectangular table with these columns:

`record,dataset,model,step,step_last,time,kind,components,uid,value_node,node_uid,entity_dim,entity_tag,x,y,z,msh,physical_tag,connectivity,value_1,...`

Records describe ranges, nodes, cells or data. Empty fields are blank. Dataset
ordinals start at zero in the first-occurrence order of the view's held model
step ranges. Every record includes its dataset ordinal and escaped model name,
so duplicate names and repeated per-model UIDs remain distinguishable.
Range records retain allocated first/last steps. Node records include
unreferenced nodes; cell records retain complete connectivity in Gmsh order,
owner, physical tag and element type. Connectivity uses semicolon-separated
node UIDs. Data records retain every stored component row, time and UID;
`value_node` is its zero-based row index within that UID. Float64 values use
17 significant digits; UIDs remain exact decimal integers.

Model VTK writes ASCII legacy unstructured grids. Linear Point, Line,
Triangle, Quadrangle, Tetrahedron, Hexahedron, Prism and Pyramid cells, plus
quadratic Line3, Triangle6, Tetrahedron10, Quad9, Hex27 and Prism18, are supported.
Complete Curve orders3–10 use VTK Lagrange Curve cells (type68), retaining every
endpoint and interior node. Other higher-order layouts fail before opening the destination. Complete quadratic
layouts use the VTK reference-node ordering, with connectivity and data bound to
the original UIDs. Pyramid14 still rejects: the standard quadratic VTK pyramid
has 13 nodes, so dropping its additional node would lose stored data.

NodeData uses point arrays, ElementData uses cell arrays, and ElementNodeData
uses separate per-cell points to preserve discontinuous values. Unreferenced
nodes remain in the grid. Absent values are NaN with a `present_STEP` mask of
zero. Stored values retain all components. VTK requires complete ElementNodeData
for each provided cell; CSV preserves any stored row width. `STEP`, `TIME`
and `ALLOCATED_RANGES` retain populated step/time bindings and sparse allocation
history. Value arrays are named `step_STEP`.

Multiple held meshes write separate `_0.vtk`, `_1.vtk`, ... files; CSV combines
them using dataset ordinals. All held payloads are validated before opening any
destination. Unknown UIDs and unsupported topology raise errors. Model CSV/VTK
append is rejected, and binary requests produce an ASCII-only warning.

These formats provide documented serialization rather than native extension
byte parity. Existing model POS/TXT behavior is retained. VTK follows the
[legacy format](https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html).
The [VTK reader](https://github.com/Kitware/VTK/blob/master/IO/Legacy/vtkDataReader.cxx)
supports `vtktypeuint64` for exact identities and NaN for missing values.
Its [quadratic tetrahedron topology](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkQuadraticTetra.cxx)
defines the target mid-edge order. Reader compatibility was reviewed in the
primary source; an external VTK reader was not run locally. Lagrange Curve
node ordering follows the two endpoints then increasing interior parameters,
as defined by [vtkHigherOrderCurve](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkHigherOrderCurve.cxx);
its degree is inferred from the complete point count. The [cell catalog](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkCellType.h)
assigns type68. Focused cubic/tenth-order NodeData and ElementNodeData checks
use independently specified curved geometry and vector values.

Model-data row tags accept checked UInt64 identities, including values beyond
signed Int64. Boolean, negative, nonintegral and overflowing tags are rejected
before mutation. This follows Tessella's owned UID storage; huge model-view tags
are not a pinned native-parity claim. Their public native observation faulted,
while the bounded Tessella regression verifies data and exact export identities.
