# Curved-domain size control — measured findings

**Goal:** a *uniform* `maxedge ≤ hmax` **and** *exact-conforming* tetrahedral mesh
of a **curved / non-box** domain (cylinder, sphere) given as a closed surface.

**Measured conclusion (design→probe workflow, grounded by running Julia):**
uniform size **and** exact input-facet conformity are **not simultaneously
achievable** for curved domains with the pragmatic routes — a fundamental tension
(it is exactly what Shewchuk's Delaunay-refinement *terminator* exists to manage,
and the terminator itself struggles on the massively-cospherical cylinder rings).
Three approaches were probed; each lands at a different, honestly-characterized
trade-off point. Approach (b) was implemented **with the exact
conformity+validity gate + insertion-order retry** and shipped as
**`mesh_sized_conforming`**: it adds an inset interior Steiner lattice and accepts
only a *gated* (conforming + valid) result, else raises an **explicit blocker** —
so it can **never** emit the invalid mesh a raw (b) does on cospherical inputs.
Measured: on a **sphere** (thick curved) it conforms and genuinely reduces the
interior edge length (`int_maxedge ≤ hmax`, tets 436→577); on a **thin** cylinder
the inset removes the interior lattice so it degrades safely to conforming-only.
What remains genuinely open is **uniform** (not graded) interior size control on
**thin / maximally-cospherical** curved domains — the cospherical-robust Shewchuk
terminator.

| approach | size control | boundary conformity | measured / blocker |
|---|---|---|---|
| **(a) background-lattice clip** (mesh_box + centroid-keep) | **uniform** interior `≤hmax`, sliver-free (45°) | **resampled** — lattice verts aren't on the surface: no-snap ⇒ *staircase* (non-conforming), snap ⇒ overshoot ~1.3× (1.269 at hmax=1.0) and ~1% inverted tets | good interior, non-input-conforming boundary |
| **(b) fine-surface + inset interior lattice** | claimed graded (deep-interior `≤hmax`) | claimed exact | **directly re-measured: FAILS on a cylinder** — the inset lattice + the cylinder's cospherical rings make `delaunay3d(perturb=false)` emit **invalid, non-conforming** output (coarse 9-gon cylinder: `valid=false`, boundary-area ≠ surface-area; fine cylinder: also cospherical-slow, >5 min). Not shippable — it emits invalid meshes on the actual enclosure geometry |
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
