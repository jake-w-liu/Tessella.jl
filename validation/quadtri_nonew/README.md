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
Source grids, triangular or mixed roots, collapsed columns, shared neighbors,
copied-source chains and closed revolutions require their separate planner
phases and are outside this first slice.

Independent certificates check the actual typed face complex, opposite
orientations of shared faces, a closed manifold boundary, Euler characteristic,
positive exact face fans, positive represented cell volumes, and exact agreement
with the emitted surface cells. Coordinates use an absolute 2e-11 bound.
Translation additionally has its analytic swept volume. Warped quadrangles use
their bilinear surface flux; no hidden diagonal or continuous curved-CAD volume
is substituted for the represented P1 mesh.

Free volume cells contain no added vertices. RecombLaterals keeps earlier hexes and
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
