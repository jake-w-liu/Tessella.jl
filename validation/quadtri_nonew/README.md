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

Run the separate quad-strip resource gate with normal compilation and bounds
checks on both runtimes; it needs no Gmsh installation or binding:

```sh
/path/to/julia-1.12.7/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/quad_strip_resources.jl
/path/to/julia-1.13.1/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/quad_strip_resources.jl
```

The strip gate produces 24 path/policy/scale rows: standalone volume, GEO,
classified projection and public API generation, with free and recombined
laterals at 1,000, 2,000 and 4,000 intervals. Warmed allocation measurements
cover construction only; extraction and the independent audit run afterward.
Constructor copies and projection or merge work performed during construction
remain counted. The audit uses the six actual source columns, independently
locates the two terminal centroids for recombined laterals, and certifies the
whole cell maps, bilinear quadrangle boundary flux, per-macro accounting and
opposite typed internal faces. A terminal centroid cell need not span both
slab planes; its actual extent must remain inside its assigned terminal macro.

Every doubling retains `allocated_next <= 2.15 * allocated_previous + 65536`.
The lowered-code check recognizes `GlobalRef(Core, :Box)` and includes an
actual boxed-local positive control. It also scans generated keyword bodies,
including optional positional wrappers, and verifies a boxed keyword control.
Input hashes must remain unchanged.
These reproducible commands describe the gate. Final strip resource and full
release results are recorded in the repository's `STATUS.md`.

The quad-strip replay reads sixteen captured Gmsh 4.15.2 P1/P2 products from
`test/artifacts/quadtri_nonew_quad_strip_oracle.toml`, without remeshing them.
Eight unit XY cases cover one or three uniform intervals, both normal signs
and both lateral policies. Eight independent variants add normalized graded
layers, rounded trapezoids, source Progression 4, YZ/XZ normal directions,
reversed Curve definitions, sparse tags and opposite transfinite pins. Exact
input hashes, raw classified lower/volume cells, coordinate and parameter bits,
primary remaps, interpolation supports, computed UV queries and independent
whole-reference-map certificates remain in the artifact. The saved actual
Progression samples and the unchanged absolute `2e-11` geometry tolerance are
authoritative. `QUADTRI_NONEW_CASE=quad_strip` selects this sixteen-case replay.

The source has six actual nodes and two Quad4 cells with one shared edge.
Free products use only the `6(N+1)` source columns. Recombined products add
exactly two Volume-owned terminal centroids and contain `2(N-1)` Hex8,
four Tet4 and ten Pyramid5 cells. Saved free family/connectivity choices are
pointer-dependent evidence, rather than native template or tag pins. Both
products must preserve actual source and cap states, all typed opposite shared
faces and exact per-macro slab volumes. A general Prism6 is integrated from
its actual reference map, including warped bilinear faces; a fixed flat
tetrahedral partition is not substituted for that map.

For N intervals the full P2 graph has `30N+15` nodes with free laterals and
`30N+31` with recombined laterals. Point, Curve and Surface ownership counts
are `8`, `8N+12` and `16N-2`; Volume counts are `6N-3` and `6N+13` respectively.
Each cap has three owned nodes and fifteen closure nodes. A lateral over a
three-node source Curve chain has `6N-3` owned nodes and `10N+5` closure nodes;
a two-node chain has `2N-1` and `6N+3`. Native interpolation nodes are checked
against the actual projected P1 carrier identities, with shared support
conformity, positive public Jacobian quadrature and the P2/P1 roundtrip.

The eight recombined captures contain 126 surface nodes with empty stored UVs.
Their independent computed inverses remain valid. Native complete computed UV
queries must reevaluate the actual coordinates; separate strip counters record
this stored-parameter provenance difference. The original 114 cases, their
serializers and earlier native CRC records remain unchanged. The new replay and
resource commands do not claim completed final release gates.

The three-Quad strip replay appends sixteen saved Gmsh 4.15.2 products from
`test/artifacts/quadtri_nonew_three_quad_strip_oracle.toml`. Eight original unit
XY fixtures cover one or three uniform intervals, both normal signs and both
lateral policies. Eight independent variants add graded normalized layers,
rounded trapezoids, skew parallelograms, YZ/XZ normals, native Progression 4,
reversed Curves, sparse tags and opposite transfinite pins. The artifact retains
all exact input hashes, raw coordinate and stored-parameter bits, classified
lower and volume cells, explicit primary remaps, interpolation support owners,
computed UV queries and 1,584 independent whole-map certificates at each order.
Saved source samples remain authoritative at the unchanged absolute `2e-11`
coordinate tolerance. Set `QUADTRI_NONEW_CASE=three_quad_strip` for this replay.

This source has eight actual boundary nodes and three strictly convex Quad4
cells in a path, with two shared edges and opposite two/four-node Curve chains.
The joined planner preserves original source cell indices and corner identities.
Free products use only `8(N+1)` retained column nodes. Recombined products add
exactly three Volume-owned terminal centroids and contain `3(N-1)` Hex8, six
Tet4 and fifteen Pyramid5 cells. Free family and connectivity choices in the
saved pointer-dependent products are evidence, not native pins. Actual cell
reference maps, opposite typed faces and each macro slab partition are checked;
the saved warped Prism6 witness uses its actual bilinear map and exact reference
integral, rather than a fixed tetrahedral proxy.

For N intervals and C centroids (zero free, three recombined), P2 has
`42N+21+8C` nodes. Point, Curve, Surface and Volume ownership counts are
`8`, `8N+20`, `24N-2` and `10N-5+8C`. Each cap has five owned nodes and 21
closure nodes. Four-node Curve laterals have `10N-5` owned and `14N+7` closure
nodes; two-node Curve laterals have `2N-1` and `6N+3`. Native supports are
compared against actual classified P1 carriers, and computed surface parameters
must reevaluate every queried node. The eight recombined captures contain 168
surface nodes with empty upstream stored UVs; separate three-strip counters
record those provenance gaps without replacing the native computed-UV contract.

The original 130 live/saved cases, serializers, counters and 136 earlier native
CRC rows remain unchanged. New three-strip counters are reported separately.
These fixture and command descriptions do not establish final release success.

Run the separate three-Quad strip resource gate with normal compilation and
bounds checks on both runtimes; it does not require Gmsh:

```sh
/path/to/julia-1.12.7/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/three_quad_strip_resources.jl
/path/to/julia-1.13.1/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/three_quad_strip_resources.jl
```

The gate measures 24 path/policy/scale rows at 1,000, 2,000 and 4,000 intervals
for standalone volume, GEO, classified projection and public API generation.
Warmed construction measurements retain constructor and downstream copies;
output extraction and independent audits follow the measurement. The audit
checks the eight actual source columns, exactly three recombined terminal
centroids, full cell reference maps, mapped bilinear face flux, every macro
partition and opposite typed faces. Each doubling keeps the unchanged bound
`allocated_next <= 2.15 * allocated_previous + 65536`. The lowered-code scan
recognizes `GlobalRef(Core, :Box)` and includes deliberately boxed positional
and keyword positive controls. Input hashes must remain stable. Final resource
and release results are recorded only after the corresponding gates complete.

The four-Quad strip extension replays 24 saved primary Gmsh 4.15.2 products
from `test/artifacts/quadtri_nonew_four_quad_strip_oracle.toml`. Eight unit
fixtures cover one/three intervals, both normal signs and both lateral policies;
sixteen additional fixtures cover rounded/skew sources, all coordinate planes,
represented graded/nonbinary levels, Progression 4, reversed curves, sparse
tags and opposite pins. The actual source has ten boundary nodes, four Quad4
cells, thirteen edges, three shared edges forming a path, and CAD chain widths
two/five. The earlier 146 scoped cases, serializers and artifact records remain
unchanged. `QUADTRI_NONEW_CASE=next_four_quad_strip` selects the 24 saved cases.

The artifact retains all exact inputs/hashes, raw coordinate/parameter bits,
actual classified lower and volume cells, primary P1/P2 identity remaps,
interpolation supports and stored/computed UV records. Promotion performs no
meshing and verifies a lossless raw JSON roundtrip. It contains 3,153 exact
whole-map records per order and five actual Pri6/Pri18 cells. Their bilinear
interfaces are certified using actual reference maps and exact mapped-flux
integration, rather than a flat tetrahedral split. Pointer-sensitive free
families/connectivity/public labels remain evidence of admissible products,
not a universal native pin. Native products must satisfy their own actual
source identities, top/lateral policy, whole maps and opposite typed faces.

For N intervals and C added centers, P1 nodes are `10(N+1)+C`; P2 nodes are
`54N+27+8C`. Free laterals add no center. RecombLaterals records exactly four
terminal problems, giving eight Tet4, twenty Pyr5 and `4(N-1)` retained Hex8.
P1 Point/Curve/Surface/Volume owner counts are `8`, `4N+8`, `6N-6`, `C`;
P2 owner counts are `8`, `8N+28`, `32N-2`, `14N-7+8C`. Actual interpolation
supports number `44N+17+7C`. Caps have seven owned/27 closure nodes; a long
lateral has `14N-7` owned/`18N+9` closure nodes, and a short lateral has
`2N-1` owned/`6N+3` closure nodes. Owner assertions follow actual classified
Line/Tri/Quad support identities, not endpoint or coordinate guesses.

The twelve recombined saved products contain 340 empty stored Surface UV
records (80 from the original eight products, 260 from the sixteen variants).
Computed inverses are retained separately. Native complete computed UV queries
must evaluate the actual nodes at the unchanged `2e-11` coordinate bound;
the dedicated counter records the stored-parameter provenance difference.
Actual Progression 4 samples are preserved, including measured deviations
from ideal analytic fractions; the coordinate bound is not relaxed.

The new terminal marker is `QUADTRI_NONEW_FOUR_QUAD_STRIP_DIFFERENTIAL_OK`,
with separate saved/P2/variant/provenance/empty-UV counters. Full execution
reports 24/24/16/12/340 for this extension. Native geometry, API, resource and
full release results require their completed guarded runs; fixture promotion
and primary certificates alone do not establish native parity.

The four-strip geometry suite independently checks 192 finished products
(all 24 original cell orders, four cyclic frame patterns and both policies).
It separately checks all 256 independent cyclic source-catalog frames for
one retained original order. The latter checks source identity and shared
edge reversal; they do not claim 24-by-256 finished-volume coverage.

Run the separate four-Quad strip resource gate with normal compilation and
bounds checks on both runtimes; it does not require Gmsh:

```sh
/path/to/julia-1.12.7/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/four_quad_strip_resources.jl
/path/to/julia-1.13.1/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/four_quad_strip_resources.jl
```

The 24 path/policy/scale rows cover 1,000, 2,000 and 4,000 intervals for
standalone volume, GEO, classified projection and public API generation.
Warmed allocations measure construction, including constructor and downstream
copies; output extraction and independent audits occur afterward. The audits
check ten actual source columns, exactly four recombined terminal centroids,
full cell reference maps, mapped bilinear face flux, every macro partition,
opposite typed faces and actual P2 carrier identities. The unchanged growth
bound is `allocated_next <= 2.15 * allocated_previous + 65536`. The lowered
code scan recognizes `GlobalRef(Core, :Box)` and uses deliberately boxed
positional and keyword positive controls. Input hashes must remain stable.
Final resource and release results require their completed guarded runs.

The rectangular-grid extension adds 24 saved primary Gmsh 4.15.2 P1/P2
products from `test/artifacts/quadtri_nonew_rect_grid_oracle.toml` and
`test/artifacts/quadtri_nonew_rect_grid_variants_oracle.toml`. The 16 original
controls cover 2-by-3 and 3-by-3 grids, N1/N3, both normal signs and both
lateral policies. Eight variants add the other coordinate planes, skew and
rounded trapezoid geometry, Progression 4/Bump 2 and grouped layer levels.
Raw Float64 bits, classified cells, actual source identities, interpolation
supports, Curve parameters and stored/computed Surface UV records are retained.
The strict extension preserves the earlier 170 cases, serializers and CRC rows.
`QUADTRI_NONEW_CASE=rect_grid` selects these 24 saved cases. The separate
`QUADTRI_NONEW_RECT_GRID_DIFFERENTIAL_OK` marker reports saved/P2/variant and
raw stored-parameter provenance counts.

The source certificate uses actual regular grid incidence, strict coherent
convex quadrangles and a Jordan boundary. Four actual native Line chains and
strictly ordered stored axis-normal columns establish the full source product.
Each physical source edge has one lateral choice per interval. B3 corner,
adjacent-B2 edge and B0 interior cells use fixed corner factories and linear
cap carry; no centroid is added. Native original source-column ordinal ranking
sets eligible terminal cap and interior lateral choices. Upstream uses primary
vertex pointer ordering, so the saved products require existence of one
consistent global rank. Their free cell families and connectivity are not
universal native pins.

For F=a*b, B=2(a+b) and N intervals, recombined products contain F(N-1) Hex8,
4F-2B Tet4 and F+B Pyr5 cells. P1 has (a+1)(b+1)(N+1) nodes and P2 has
(2a+1)(2b+1)(2N+1). Typed supports and carrier ownership come from actual
primary identities. Independent exact literal-reference Vandermonde maps
certify actual Tet10, Hex27, Prism18 and rational Pyramid14 nodal geometry
with full-reference Bernstein bounds and exact volume integrals. Deliberately
folded corner and interior support controls reject invalid quadratic maps.
Empty raw upstream stored Surface UVs are counted separately; native computed
parameters must reevaluate the queried coordinates at unchanged 2e-11.

Run the rectangular-grid resource gate with normal compilation and bounds
checks on both supported runtimes; it does not require Gmsh:

```sh
/path/to/julia-1.12.7/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/rect_grid_resources.jl
/path/to/julia-1.13.1/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/rect_grid_resources.jl
```

The gate measures 24 P1 path/policy/scale rows at 500, 1,000 and 2,000
intervals for standalone volume, GEO, classified projection and public API.
Separate actual source-size series measure 1,000, 2,000 and 4,000 quadrangles.
Construction allocations include constructor/downstream copies; extraction
and independent audits follow measurement. The unchanged doubling bound is
`allocated_next <= 2.15 * allocated_previous + 65536`. Six P2 audits cover
2-by-3 N1/N3 and 3-by-3 N3 with both policies, including actual B0 Hex27/Prism18
supports. Full-cell P1 maps, parent partitions, opposite typed faces, P2
ownership/closure/UVs and true `Core.Box` controls are checked. Input hashes
must remain stable. Release success requires completed guarded resource,
strict and package runs recorded in STATUS.md.

## Arbitrary-length recombined B4 strips

The strict driver adds twelve independently captured M5/M7/M9 all-boundary
strips, with one or three intervals, both normal signs, unit/skew/rounded
geometry and literal Progression/Bump controls. Select them with
`QUADTRI_NONEW_CASE=b4_rec_strip`. The earlier194 cases retain their contracts.
The saved artifact keeps exact recipes, source parameters, Float64 bits,
typed cells, real centers, full-reference map bounds, interpolation supports,
carrier proofs and stored/computed UV records. The pinned oracle's body-center
warnings remain recorded; the driver checks their exact expected form and
rejects Error diagnostics. The twelve captures contain84 real terminal
centers and384 empty stored Surface UV records. Native computed UVs must
reevaluate the actual queried coordinates within the unchanged2e-11 tolerance.

The source admits exactly one Quad across and at least five along, using the
same actual regular-disk incidence, strict convexity and Jordan boundary
proof as the rectangular path. Original source-column ordinals determine one
consistent terminal rank; public oracle node tags do not substitute for its
pointer order. Earlier macros are whole Hex8 cells. Each terminal macro emits
its actual strictly interior mean and a fully certified seven-cell fan.
For M source Quads and N intervals, family counts are2M Tet4, M(N-1) Hex8
and5M Pyr5; P1/P2 node counts are2(M+1)(N+1)+M and(12M+6)N+14M+3.
Actual primary/support identities determine carriers, including the real
Volume centers and their radial supports. Independent exact reference maps
certify every actual linear/quadratic cell and macro integral.

Two retained M7 Progression4 recipes exposed a preexisting source-sampling
error on a length-two offset Line. The native GEO Line density route now uses
the pinned bounded1e-5 first derivative and numerically integrated length.
It preserves public exact evaluation and the separate recombination-count
protocol. The strict native products pass the original coordinate tolerance;
the artifact preserves the earlier failed source-only evidence as provenance.
Other unfinished source categories remain work under the active full-parity
goal. Resource and final package release results
are recorded separately in STATUS.md after the candidate is frozen.

Run the recombined B4 resource gate with normal compilation and bounds checks
on both supported runtimes; it does not require Gmsh:

```sh
/path/to/julia-1.12.7/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/b4_rec_strip_resources.jl
/path/to/julia-1.13.1/bin/julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/b4_rec_strip_resources.jl
```

The gate measures twelve P1 construction rows: five source Quads with 500,
1,000 and 2,000 intervals through standalone volume, GEO, classified projection
and public API generation. A separate one-interval source series uses 1,000,
2,000 and 4,000 actual Quads and their real terminal centers. Independent
audits follow the measured construction and copies. The allocation bound stays
`allocated_next <= 2.15 * allocated_previous + 65536`. Three P2 products cover
M5/N1, M5/N3 and M7/N3, with actual support identities, carriers, closure and
UV reevaluation. Full-cell interval bounds, actual mean enclosures, macro
integrals, opposite typed internal faces and the complete shell are checked.
Actual named and generated keyword/closure bodies are scanned for true
`Core.Box`, with positive controls. Release success requires completed guarded
resource, strict and package runs with stable inputs and index.

The aggregate source-area certificate cancels exact dyadic signed triangle
areas; a fixed ulp bound on summing thousands of rounded cell areas is not
used. Positive per-cell areas and their independent map tolerances remain
checked. Confirmed container-growth failures are corrected through checked
capacity hints in tag validation, face topology/projection and API ownership
and merge construction, without changing mesh identities or the allocation
growth bound.

## Arbitrary-length free B4 strips

The strict driver replays twelve original free M5/M7/M9 primary captures with
one or three intervals and both normal signs. Select them with
`QUADTRI_NONEW_CASE=b4_free_strip`. Original recipes, coordinate/parameter bits,
typed cells, primary remaps, interpolation supports and stored/computed UV
payloads remain literal authority. These captures contain no retained center;
their pointer-dependent families and diagonals are independently certified
rather than used as native connectivity pins. The earlier 206 cases retain
their contracts and every prior CRC artifact remains unchanged.

The native planner acts on the actual certified source incidence. It finishes
all physical face decisions before choosing an existing-corner factory.
Retained problems append their actual eight-corner mean in original
cell/interval order and emit a complete fan, including nonterminal intervals.
The permanent suite separately checks all 729 final masks, the 253 supported
literal C++ factory complexes, and every retained-center fan with exact whole
maps, opposite internal faces and complete shells. Actual centerful Source
products exercise operation Scope, classified Projection and public seeded
mesh lifecycle; ordinary recipes separately exercise real generate(3).

Run the allocation gate with normal compilation and bounds checks:

```sh
julia --startup-file=no --project=. --check-bounds=yes validation/quadtri_nonew/b4_free_strip_resources.jl
```

It measures actual plan, emission, standalone, GEO, projection, classification,
merge and API paths at 500/1000/2000 intervals, with a separate 1000/2000/4000
source series at three intervals. Actual P2 supports, carriers, UVs and two
retained-center Source witnesses receive independent audits. Allocation
doubling retains `allocated_next <= 2.15 * allocated_previous + 65536`; actual
lowered method bodies include keyword and closure bodies with true Core.Box
positive controls, including generation reconciliation, independent-record
label validation, order changes, CAD boundary queries and the shared native
support-node parameter paths. Public lifecycle checks independently cover
raw native records, discrete records, aggregate closure repetitions, repeated
generation and atomic rejection of unsupported label collisions.
The joint inventory also scans public renaming, Point identity binding,
metadata/refinement/order adapters and API recombination bodies. Analytical
disjoint unit-right triangles measure actual label-table allocation and retained
size at 1000/2000/4000 cells, with unchanged `2.15*previous+65536` growth bounds.
These measurements are additional to the existing free geometry growth rows.
The focused joint AST check covers 234 actual bodies with zero Core.Box on
both supported runtimes; final resource runs must bind the final joint tree.
Repeated 2D Curve retention and global quadratic discrete-Curve elevation
remain separate confirmed generation-phase implementation work. The allocator
fixture uses independently proved stable discrete Point records and preserves
the original unresolved-reference rejection controls.
Release success requires completed guarded package,
strict and all family resource checks, recorded in STATUS.md.
