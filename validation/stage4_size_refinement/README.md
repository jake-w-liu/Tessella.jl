# Stage-4 size refinement — measured findings

Why Tessella ships **no** 3-D interior size-refiner (a non-converging one would
fail the correctness bar). Two concrete approaches were built and **measured** to
fail the size bound on a convex box `[0,4]³`:

| script | approach | result |
|---|---|---|
| `diverges_interior_midpoint.jl` | insert the midpoint of every edge > `hmax` into the Delaunay, repeat | **diverges** — max edge pinned at `4√2 ≈ 5.657` (the box face diagonal) forever; +3 verts/+7 tets per pass, no decrease |
| `fine_surface_leaves_interior_diagonals.jl` | mesh the box surface finely (k×k per face) first, then tetrahedralize | correct volume + valid, but **interior diagonals persist** (`maxedge ≈ 4.0` at k=4) — `tetrahedralize` adds no interior points |

**Mechanism.** The face diagonal is a *boundary edge* of the input surface;
splitting it at its midpoint (the face centre) regenerates a diagonal of equal
length among the remaining corners. Interior point insertion can never shorten a
fixed boundary edge, and a fine boundary alone leaves long interior diagonals.

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
