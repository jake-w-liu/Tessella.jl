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
quadratic Line3, Triangle6, Tetrahedron10, Quad8/9, Hex20/27, Prism15/18 and
Pyramid13, are supported.
Complete Curve orders3–10 use VTK Lagrange Curve cells (type68), retaining every
endpoint and interior node. Complete Triangle/Quadrangle/Tetrahedron orders3–10
and Hexahedron/Prism orders3–9 use VTK Lagrange types69–73, preserving every
node. Those are the existing catalog limits; no new element family is invented.
Incomplete higher-order layouts fail before opening the destination. Quadratic
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
The quadratic layouts use a bijection between canonical corner/edge/face
identities, preserving every stored node and value. Edge-only cells use
VTK types 23/25/26/27; full Quad9/Hex27/Prism18 use 28/29/32. The
[quadrangle](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkQuadraticQuad.cxx),
[hexahedron](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkQuadraticHexahedron.cxx),
[wedge](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkQuadraticWedge.cxx)
and [pyramid](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkQuadraticPyramid.cxx)
classes define the target edge order. Pyramid edge identities use the Euclidean
corner/edge lattice rather than VTK's collapsed parametric apex coordinates.
The new opt-in test checks all four edge-only families against independently
specified geometry and vector values; only two cheap Pyramid13 assertions join
the default suite. Its [quadratic tetrahedron topology](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkQuadraticTetra.cxx)
defines the target mid-edge order. For the earlier quadratic/Curve increment, reader compatibility was reviewed
in the primary source; no external reader was run for those fixtures. Lagrange Curve
node ordering follows the two endpoints then increasing interior parameters,
as defined by [vtkHigherOrderCurve](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkHigherOrderCurve.cxx);
its degree is inferred from the complete point count. The [cell catalog](https://github.com/Kitware/VTK/blob/master/Common/DataModel/vtkCellType.h)
assigns type68. Focused cubic/tenth-order NodeData and ElementNodeData checks
use independently specified curved geometry and vector values.

Complete higher-order Surface/Volume output uses VTK5.1 `OFFSETS` and
`CONNECTIVITY` cell arrays. Earlier VTK headers trigger a reader conversion
that swaps two Hex edge blocks, even when the stored order is correct. The
new writer enumerates the canonical integer reference lattice and certifies a
bijection to every owned reference node. Prism interiors use flat triangular
rows rather than recursive Triangle interior ordering. Existing linear,
quadratic and Curve-only files keep their earlier2.0 encoding.

Focused qualification checks40 public VTK reference lattices and NodeData/
ElementNodeData on10 supported catalog models. Actual VTK9.7.1 reads all20
exports with matching geometry, cell type, UInt64 node/cell identities and
every component in two vector time steps. This bounded reader result does
not establish all-family native byte parity. VTK is not a package or default
test dependency. Only two tiny cubic-Triangle checks enter the default suite;
the166 serialization/error assertions remain opt-in.

Model-data row tags accept checked UInt64 identities, including values beyond
signed Int64. Boolean, negative, nonintegral and overflowing tags are rejected
before mutation. This follows Tessella's owned UID storage; huge model-view tags
are not a pinned native-parity claim. Their public native observation faulted,
while the bounded Tessella regression verifies data and exact export identities.
