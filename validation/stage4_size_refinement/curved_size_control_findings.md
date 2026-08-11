# Curved-domain size control — measured findings

**Goal:** a *uniform* `maxedge ≤ hmax` **and** *exact-conforming* tetrahedral mesh
of a **curved / non-box** domain (cylinder, sphere) given as a closed surface.

**Measured conclusion (design→probe workflow, grounded by running Julia):**
uniform size **and** exact input-facet conformity are **not simultaneously
achievable** for curved domains with the pragmatic routes — a fundamental tension
(it is exactly what Shewchuk's Delaunay-refinement *terminator* exists to manage,
and the terminator itself struggles on the massively-cospherical cylinder rings).
Three approaches were probed; each lands at a different, honestly-characterized
trade-off point:

| approach | size control | boundary conformity | measured / blocker |
|---|---|---|---|
| **(a) background-lattice clip** (mesh_box + centroid-keep) | **uniform** interior `≤hmax`, sliver-free (45°) | **resampled** — lattice verts aren't on the surface: no-snap ⇒ *staircase* (non-conforming), snap ⇒ overshoot ~1.3× (1.269 at hmax=1.0) and ~1% inverted tets | good interior, non-input-conforming boundary |
| **(b) fine-surface + inset interior lattice + `recover_boundary`** | **graded**: boundary edges `≤hmax` (fine input surface) and deep-interior `≤hmax` (lattice diagonal = hmax, measured 1.0 at hmax=1.0), **but a ~hmax-thick boundary shell at ~2·hmax** | **exact** (enforced by the `recover_boundary` `_rb_gate`) | convex-only; a non-convex concavity's inset point can pierce a region ⇒ the gate rejects |
| **(c) recover-then-protect** (insert interior circumcenters that don't encroach boundary faces) | **stalls**: boundary diametral balls blanket the interior — 48/48 circumcenters rejected on the cylinder; maxedge stays ~1.4–2.2 ≫ hmax=0.6 | exact at the seed | + **cospherical invalidity**: `delaunay3d(perturb=false)` on the circular rings emits positively-oriented **overlapping** tets (vol 6.0→9.5), *not* caught by `validate` — the BCC/lattice degeneracy the parent README documents |

**What is achievable and shipped instead:**
- **Axis-aligned** uniform size control + conforming multi-region + CSG — fully
  solved (`mesh_box` / `mesh_box_regions`), exact, sliver-free.
- **Curved boundary conformity at the input-surface resolution** — solved
  (`recover_boundary`): give it a finely-triangulated curved surface (`MeshSurface`)
  and the boundary is size-controlled and exactly conforming; only the *interior*
  lacks independent size control.

**Open (research-grade):** simultaneous *uniform interior* size control **and**
*exact* boundary conformity on curved domains — the boundary-protected Shewchuk
terminator with a cospherical-degeneracy-robust cavity strategy. This session
measured that all three pragmatic shortcuts fall short (above); a correct version
is a focused kernel effort, not a quick addition.
