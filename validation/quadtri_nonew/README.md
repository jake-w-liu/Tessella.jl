# QuadTriNoNewVerts differential

Run against the pinned Gmsh 4.15.2 Julia binding and runtime:

```sh
GMSH_JULIA_API=/path/to/gmsh.jl julia --project=. --check-bounds=yes validation/quadtri_nonew/differential.jl
```

`QUADTRI_NONEW_CASE` filters fixture names. `QUADTRI_NONEW_ORACLE_ONLY=1`
certifies the Gmsh side while native integration is being completed.

The P1 matrix contains 20 cases with one nondegenerate source quadrangle,
positive uniform and grouped layer counts ending at normalized height 1.0,
free and recombined laterals, translation in both normal directions, tilted
translation, exterior-axis rotation by Pi/6, and twist by Pi/6. Four uniform
affine cases also exercise full P2 elevation and return to P1. Two independent
region cases preserve separate node identity, including coincident geometry
created with point setters that do not run coherence. The four non-unit final
height cases are expected upstream errors, counted separately from valid meshes.
Four native multi-turn helical cases use 52 layers, angle 13Pi/6 and separated
axial translations of +/-3. Their public products receive independent geometry
certificates; full Gmsh parity is not claimed. Six large-angle Gmsh probes are
recorded separately: free laterals at pitches +/-3 fail with
`FindDiagonalEdgeIndices`, while recombined laterals at pitches +/-3 and +/-0.1
take a different Tet/Prism/Pyramid route. Observed products contain 212 or 213
referenced volume nodes; one allocator-dependent outcome adds a volume-owned
centroid in the first interval. Their actual node and family counts and local
certificate outcomes are reported, without accepting them as valid parity
fixtures.
The triangle extension adds 36 P1 cases with one plain source triangle: both
source windings, positive and negative normal translation, tilted translation,
uniform and graded groups, and mild rotation/twist in both angle directions.
Free laterals produce three tetrahedra per interval; RecombLaterals produces
one six-node prism per interval. Both use exactly three column nodes per level,
with no centroid or collapsed fourth column. Source and cap remain triangles.
Four affine triangle cases certify full Tet10/Prism18 support and the actual
public P2/P1 roundtrip. Two more cases check separate and coincident independent
source identities. This extends the original quad matrix and leaves its 24
native CRC records unchanged.

The bounded two-Tri extension replays twelve already captured Gmsh 4.15.2
products from `test/artifacts/quadtri_nonew_two_tri_oracle.toml`; it does not
regenerate those oracle meshes. The exact inputs use one four-sided TF2 source
with two actual Tri3 cells, normal translation in both directions, free and
recombined laterals, and one, three uniform, or three graded intervals. The
catalog retains input hashes, all raw P1/P2 node and lower/volume cell data,
classified owners, primary identity remaps, interpolation supports, and stored
and computed parameter payloads. `QUADTRI_NONEW_CASE=nonew_two_tri` selects these
twelve cases. Separate terminal counters preserve the original 86-case evidence.
The companion `test/artifacts/quadtri_nonew_two_tri_crc.txt` pins 24 native API
P1/P2 products for the same matrix. It uses the existing strict public
node/cell/owner/parameter serializer, extended with sorted actual support
catalogs and complete computed surface own/closure queries. The 24 rows are
byte-identical on normal Julia 1.12.7 and 1.13.1 with bounds checks and are
asserted by `api_nonew_two_tri_boundary_test.jl`. The earlier 40 triangle and
24 quadrangle records remain unchanged.

Each two-Tri P1 product has `4(N+1)` column nodes and either `6N` tetrahedra or
`2N` prisms. Its actual source winding and diagonal must match the saved source;
both native and saved volume products receive independent typed-face,
opposite-incidence, convex-source/product and whole Tet4/Pri6 map certificates.
Free volume template choices may differ from the saved pointer-selected
connectivity. Every native mesh must match its own finalized external surfaces.
The full P2 graph has `18N+9` nodes, owned by Points/Curves/Surfaces/Volume in
counts `8`, `8N+4`, `8N-2`, `2N-1`. Each of the four laterals has `2N-1` owned
nodes and `6N+3` nodes including its boundary. Both cap diagonal midpoints are
Surface-owned; internal shared supports are Volume-owned without an invented
CAD surface. Native P2 support identity is checked against its own actual
projected P1 Line/Tri/Quad carriers and interpolation, followed by actual public
Jacobian quadrature and the return to P1. Saved raw tags are not native tag
allocation pins. Coordinate comparisons retain the absolute `2e-11` bound.

The six recombined saved fixtures contain 56 lateral face centers with empty
stored UVs. They have valid independently computed geometric inverses. The
native query contract returns complete computed own/closure UV pairs and must
evaluate them back to the actual node coordinates. Separate provenance counters
record this difference; the saved raw partial arrays remain in the artifact.

Run the independent allocation and payload gate with normal compilation:
`julia --project=. --check-bounds=yes validation/quadtri_nonew/two_tri_resources.jl`.
It needs no Gmsh installation. Both lateral policies pass through standalone,
GEO, classified projection and public API generation at 1,000, 2,000 and 4,000
intervals. Each measured output receives a complete linear audit of actual
column identities, cell maps, partition volumes and opposite internal face
cycles; the four paths must produce the same geometry digest. Allocation
doubling and lowered helper code are checked separately, with every production
and certificate input hashed before and after. Direct indexed construction
uses certified columns for the volume and its lateral surfaces, avoiding
temporary coordinate-cell arrays and coordinate welding dictionaries. API
classification and merge tables reserve their known payload bounds. Normal
mesh constructors and boundary/projection copies remain part of these measured
allocations.

The two recombined triangle P2 cases also record a parameter-provenance gap.
Both policies require five native owned nodes and 21 nodes with boundaries on
each lateral, with complete computed parameter arrays of lengths 10 and 42.
The free Tet10 Gmsh route has the same parameter lengths. Separate public units
check ordinary Tri6 boundary curves, independent coincident Tet10 carriers,
identity-preserving mutations and the refined child-face supports.
For the recombined Gmsh route, `getNodes` supplies only four stored parameter
values and `getNode` supplies no stored UV for its three Prism18 face centers
on each planar ruled lateral. Nine face centers per case therefore have empty stored
parameters. The native 2D/3D cache retains its existing computed-UV query
contract: every returned UV must evaluate back to the actual node. With boundary
nodes included, the fixture returns 21 nodes with 36 stored parameter values in
Gmsh and 42 computed values natively. Direct geometric parametrization is
checked independently. These two cases count the provenance difference
explicitly; they do not claim exact stored-parameter parity. Complete stored
2D/3D parameter provenance remains part of the broader public lifecycle work.

Other source grids, mixed roots, collapsed columns, shared neighbors,
copied-source chains and closed revolutions require their separate planner
phases and are outside these bounded slices. The two-Tri category additionally
requires a strictly convex source in an axis-aligned plane with straight CAD
curves and normal translation. Rotations, twists and tilted translation of this
source category remain separate planner work.

Independent certificates check the actual typed face complex, opposite
orientations of shared faces, a closed manifold boundary, Euler characteristic,
positive exact face fans, positive represented cell volumes, and exact agreement
with the emitted surface cells. Coordinates use an absolute 2e-11 bound.
Translation additionally has its analytic swept volume. Warped quadrangles use
their bilinear surface flux; no hidden diagonal or continuous curved-CAD volume
is substituted for the represented P1 mesh.
For triangles, the full Tet4 Jacobian is checked with exact rational arithmetic.
The Pri6 determinant is differentiated directly from its six reference shape
functions. It is affine over the reference triangle and quadratic axially;
exact endpoint and stationary-minimum checks at all three triangle vertices
therefore certify the entire reference prism, including warped cells.

Free volume cells contain no added vertices. In the quadrangle slice,
RecombLaterals keeps earlier hexes and
adds exactly one centroid in the final interval, producing two tetrahedra and
five pyramids there. The upstream NoNew top split uses vertex pointer ordering;
the native planner uses deterministic source-node precedence. Both opposite
top diagonals are admissible, and each mesh must match its own finalized cap
and lateral surface plan.

The recorded Gmsh oracle demonstrates real allocator dependence. The same
one-quad, three-layer, free-lateral input can contain 16 Tet + 1 Pyr,
12 Tet + 3 Pyr, or 8 Tet + 5 Pyr. Counts alone are insufficient; every outcome must pass
the independent certificates. Exact equality to a random oracle connectivity
is not required. RecombLaterals counts and centroid are stable, while its two
observed typed signatures correspond to the two admissible cap diagonals.
Native repeated execution must preserve its complete mixed-mesh CRC. Global
GEO meshes retain the existing CAD/control point parts; these are separate from
the volume's no-added-vertex rule.
Triangle free-lateral signatures can also differ between Gmsh allocations.
Their family counts, actual typed complex, finalized boundary and full P1 maps
remain mandatory certificates; a random oracle's diagonal choices are not pins.

Geometry units additionally cover standalone surface/volume and API/GEO
consistency, equivalent regroupings, affine full-P2 elevation and return to P1,
and explicit unsupported-input failures with preserved API state. The curved
cases test standalone source/cap/lateral queries, API generation 2, volume
projection with node renumbering, physical names, and atomic rejection of a
changed coordinate or winding. API generation 3 retains the established
top-dimensional actual-cell cache contract: volume cells and complete node
ownership are queried directly, while missing lower-dimensional cells are not
synthesized. Classified `model_to_mixed` and GEO products contain the actual
planned boundary cells.

P2 certificates evaluate independent P1 reference interpolation at each
quadratic support, verify shared node identity and complete primary cell
preservation, and integrate positive Gauss4 Jacobians to the analytic affine
volume. Gmsh and native P2 meshes are certified against their own P1 cells,
because upstream pointer choices can select different admissible complexes.
High-order `.geo SetOrder` and curved CAD P2 placement remain separate blockers.

Full or multiple pure revolutions require global intersection planning. Positive
local volume alone does not certify nonadjacent interiors are disjoint. The
native guard fails before
publishing discretizations, cache state or allocation counters. Terminal heights
0.5 and 2.0 instead fail in Gmsh because the CAD cap remains at u=1; native
generation rejects these precisely before allocation.
The former native thin-pitch meshes have a common strict-interior witness in
their actual Hex8 cells from intervals 1 and 49, verified by independent
trilinear inversion. This proof is separate from the different Gmsh cell route.
Both thin-pitch signs are blocked; both separated-pitch signs remain supported.
The conservative global hull certificate also rejects undecidable separation
without asserting that every such input necessarily overlaps.

A captured former native Hex8 also demonstrates why positive total volume and
valid face incidence do not establish a regular cell map: its exact rational
trilinear Jacobian is negative at a strictly interior reference point. Public
standalone, GEO and API tests reject this input atomically. The production
full-domain Jacobian certificates cover each emitted Hex8, Prism6 and Pyramid5;
undecidable positivity fails with a precise pending-planner diagnostic.

Native artifact builders live in
`test/geometry/quadtri_nonew_certificates.jl` as `artifact_records()`. Their
digests include exact coordinate bits, oriented typed connectivity, entity
ownership, physical names and metadata. Golden native records complement the
geometry certificates; random Gmsh connectivity is never used as a golden hash.
The 24 records in `test/artifacts/quadtri_nonew_crc.txt` cover 20 global GEO P1
products and four API volume P2 products. CAD point parts use sorted entity tags
for NoNew products, keeping their exact output stable across dictionary insertion
and rehashing. These pins agree on Julia 1.12.7 and 1.13.1 with bounds checks.

The separate 40 records in `test/artifacts/quadtri_nonew_triangle_crc.txt` cover
36 triangle global GEO P1 products and four actual public API P2 products. Build
them with `QuadTriNoNewTriangleCertificates.artifact_records()` from
`test/geometry/quadtri_nonew_triangle_certificates.jl`. Each tab-separated row
contains `name`, `kind`, `n_nodes`, `n_cells`, `families`, and `sha`. The global
digest combines the exact global mesh, every sorted entity mesh part and the
classified volume product. The public P2 digest reads actual node and element
queries, including external labels, coordinate bits, owner dimensions/tags,
computed parameters, oriented connectivity and physical names. It therefore
includes the simplex high-order overlay as well as native mixed cells.

Regeneration uses normal Julia compilation and `--check-bounds=yes` on both
Julia 1.12.7 and 1.13.1. Compare all 40 complete rows before replacing the pin
file; retain the independent geometry certificates and the original 24 quad
records. A changed digest is evidence to investigate, rather than a reason to
skip an ownership or parameter query.

The separate quad-patch replay reads
`test/artifacts/quadtri_nonew_quad_patch_oracle.toml`. It promotes twelve saved
Gmsh 4.15.2 rectangular products and four independently captured convex
trapezoids without remeshing the oracle. Inputs cover one or three uniform
intervals, three graded intervals, both normal signs and both lateral policies;
the trapezoids cover the graded three-interval profile. Full original P1/P2
coordinates and bits, actual lower and volume cells, owner tags, stored node
parameters, packed parameter queries, primary identity remaps and interpolation
support carriers remain in the artifact. The sixteen input and JSON hashes
identify the exact captures.

The native source must retain nine actual vertices and four Quad4 cells sharing
one interior pivot. The replay independently certifies the actual sampled source
partition, each emitted Tet4/Pyramid5/translated Hex8 map, macro volumes and
opposite internal typed faces, pivot cap diagonals and finalized exterior cells.
P2 construction checks actual primary identities and shared interpolation
supports, positive public quadrature Jacobians and integrated swept volume.
Every new P2 node is classified against its native actual P1 lower-cell carrier,
which also verifies the native three-node source and copied top Curve chains.
The recorded rounded trapezoid has two source Curve middle samples exactly
noncollinear with their endpoint chords; actual sampled topology is authoritative.
Raw pointer-selected oracle volume connectivity remains evidence rather than a
native tag or template pin.

For N intervals the bounded patch has 9(N+1) P1 nodes and 50N+25 P2 nodes. P2
owner counts are 8 Point, 8N+20 Curve, 24N+6 Surface and 18N-9 Volume nodes.
Each cap has nine owned nodes and 25 closure nodes; each lateral has 6N-3 owned
nodes and 10N+5 closure nodes. Native computed UVs must reevaluate all actual
nodes. Eight recombined saved fixtures have 160 surface nodes with empty upstream
stored UVs; these are counted provenance differences. The driver appends separate
quad-patch counters while preserving the original 86 oracle samples, twelve
two-triangle saved cases and earlier native artifact pins. Set
`QUADTRI_NONEW_CASE=next_quad` for the sixteen-case targeted replay.

Run the separate quad-patch resource gate from the repository root with normal
compilation and bounds checks on both runtimes. Replace the executable paths
below with the installed Julia 1.12.7 and 1.13.1 binaries; no Gmsh installation
or binding is needed. An optional final argument naming a file in an existing
directory writes the measured rows and input hashes as TOML.

```sh
/path/to/julia-1.12.7/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/quad_patch_resources.jl
/path/to/julia-1.13.1/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/quad_patch_resources.jl
```

Each run produces 24 rows: standalone volume, GEO, classified projection and
public API generation, for free and recombined laterals at 1,000, 2,000 and
4,000 intervals. After warming each operation, the helper takes the minimum of
three construction allocation measurements and three construction timings.
Output extraction and the independent geometry audit are outside those
measurements; normal constructors, returned-mesh copies and projection/merge
work inside the operation remain counted. Each path must preserve the same
actual geometry digest, complete column identities, typed families, macro
partition volumes and opposite internal face incidences.

For each path and policy, every doubling retains the unchanged gate
`allocated_next <= 2.15 * allocated_previous + 65536`. The lowered-code scan
recognizes `GlobalRef(Core, :Box)` as well as direct boxes, and must first detect
a deliberately boxed captured local in an actual lowered positive control.
Production and helper inputs are hashed before and after the run. The terminal
`QUAD_PATCH_RESOURCE_OK` marker follows the assertions; documenting these
commands does not establish final resource or full-validation success.
