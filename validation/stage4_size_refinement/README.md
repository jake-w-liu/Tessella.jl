# Stage-4 size refinement — measured findings

**Resolution (shipped):** `works_mesh_box_kuhn.jl` demonstrates
[`mesh_box`](../../src/Mesh3D.jl) — a **correct, sliver-free, size-controlled**
tet mesher for axis-aligned boxes (Kuhn/Freudenthal explicit subdivision):
`maxedge ≤ hmax`, exact volume, `validate.ok`, watertight (boundary χ=2), min
dihedral 45°/≥42° — for arbitrary `hmax`. It covers the enclosure's box regions.
The three failures below are *why the naive routes don't work* and what pointed
at the explicit-connectivity fix `mesh_box` uses; general **non-box adaptive**
refinement (arbitrary surface) remains the research-grade open item.

Three naive approaches were built and **measured** to fail on a convex box `[0,4]³`:

| script | approach | result |
|---|---|---|
| `diverges_interior_midpoint.jl` | insert the midpoint of every edge > `hmax` into the Delaunay, repeat | **diverges** — max edge pinned at `4√2 ≈ 5.657` (the box face diagonal) forever; +3 verts/+7 tets per pass, no decrease |
| `fine_surface_leaves_interior_diagonals.jl` | mesh the box surface finely (k×k per face) first, then tetrahedralize | correct volume + valid, but **interior diagonals persist** (`maxedge ≈ 4.0` at k=4) — `tetrahedralize` adds no interior points |
| `bcc_lattice_delaunay_invalid.jl` | BCC lattice points (spacing `hmax/√2`) + Delaunay | works for the trivial even case (`hmax=3` → 2×2×2, valid, **60° dihedral**, exact volume) but **invalid/inverted tets for general spacing** (`hmax=2`: valid=false, 0° dihedral; `hmax=1`: volume wrong at 64.44) |

**Mechanisms.** (1)/(2) The face diagonal is a *boundary edge* of the input
surface; splitting it at its midpoint (the face centre) regenerates a diagonal of
equal length among the remaining corners — interior insertion can never shorten a
fixed boundary edge, and a fine boundary alone leaves long interior diagonals.
(3) A regular lattice is *maximally* cospherical-degenerate; Delaunay of
cospherical points is ambiguous and the exact+SoS kernel resolves ties
deterministically but not always into a *valid* tetrahedralization — a correct BCC
mesh must emit the **known BCC connectivity explicitly**, not Delaunay the lattice.

**Conclusion (measured, not assumed).** A genuine 3-D size bound requires
**boundary-conforming Delaunay refinement** — Steiner points inserted *and*
boundary sub-faces split under an encroachment rule (Shewchuk's terminator),
which is divergence-prone near small input angles and is a research-grade
component, out of single-session scope.

**Acceptance test for a future correct implementation.** On the box, for a sweep
of `hmax`: `refine → maxedge ≤ hmax`, with `validate.ok`, `is_delaunay3`, exact
preserved volume, and a bounded vertex count (termination).

Run either script with:

```
julia --project=. validation/stage4_size_refinement/<script>.jl
```
