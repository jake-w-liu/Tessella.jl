# Curved-domain size control — measured findings

## Current continuation status (2026-10-10)

This is a historical record of three exploratory strategies and the gated
implementation they motivated; the capture date is not recorded here. All
captured table values are retained, and no probe was rerun for this review.
The failures are evidence about those strategies, not mathematical impossibility
or qualification of today's Root V24/release V6 candidates.

See the [current handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
for live reader ownership, remaining full-parity work, and the normal verified
push to `main`. Current verification must cover Julia 1.12.x and 1.13.x on the
same bound source. Keep uniform size, boundary conformity, validity, termination,
and thin/cospherical cases as separate acceptance obligations.

**Goal:** a *uniform* `maxedge ≤ hmax` **and** *exact-conforming* tetrahedral mesh
of a **curved / non-box** domain (cylinder, sphere) given as a closed surface.

**Measured conclusion (design→probe workflow, grounded by running Julia):**
the three tested pragmatic routes did **not simultaneously achieve** uniform
size and exact input-facet conformity on the recorded curved cases. That scoped
failure motivates further boundary-protected refinement and degeneracy handling;
it does not establish that the goal is mathematically impossible.
Three approaches were probed; each lands at a different, honestly-characterized
trade-off point. Approach (b) was implemented **with the exact
conformity+validity gate + insertion-order retry** and shipped as
**`mesh_sized_conforming`**: it adds an inset interior Steiner lattice and accepts
only a *gated* (conforming + valid) result, else raises an **explicit blocker** —
so it can **never** emit the invalid mesh a raw (b) does on cospherical inputs.
Measured: on a **sphere** (thick curved) it conforms and genuinely reduces the
interior edge length (`int_maxedge ≤ hmax`, tets 436→577); on a **thin** cylinder
the inset removes the interior lattice so it degrades safely to conforming-only.
This historical probe left **uniform** (not graded) interior size control on
**thin / maximally-cospherical** curved domains unfinished. Keep that requirement
active until a current implementation and independent qualification satisfy it.

| approach | size control | boundary conformity | measured / blocker |
|---|---|---|---|
| **(a) background-lattice clip** (mesh_box + centroid-keep) | **uniform** interior `≤hmax`, sliver-free (45°) | **resampled** — lattice verts aren't on the surface: no-snap ⇒ *staircase* (non-conforming), snap ⇒ overshoot ~1.3× (1.269 at hmax=1.0) and ~1% inverted tets | good interior, non-input-conforming boundary |
| **(b) fine-surface + inset interior lattice** | claimed graded (deep-interior `≤hmax`) | claimed exact | **directly re-measured: FAILS on a cylinder** — the inset lattice + the cylinder's cospherical rings make `delaunay3d(perturb=false)` emit **invalid, non-conforming** output (coarse 9-gon cylinder: `valid=false`, boundary-area ≠ surface-area; fine cylinder: also cospherical-slow, >5 min). Not shippable — it emits invalid meshes on the actual enclosure geometry |
| **(c) recover-then-protect** (insert interior circumcenters that don't encroach boundary faces) | **stalls**: boundary diametral balls blanket the interior — 48/48 circumcenters rejected on the cylinder; maxedge stays ~1.4–2.2 ≫ hmax=0.6 | exact at the seed | + **cospherical invalidity**: `delaunay3d(perturb=false)` on the circular rings emits positively-oriented **overlapping** tets (vol 6.0→9.5), *not* caught by `validate` — the BCC/lattice degeneracy the parent README documents |

The table describes the **raw experimental routes**. In particular, raw (b)'s
recorded failure is separate from the gated `mesh_sized_conforming` wrapper
described above; accepting only a valid conforming result does not itself prove
uniform size control.

**Recorded implemented alternatives:**
- **Axis-aligned** uniform size control + conforming multi-region + CSG — fully
  solved (`mesh_box` / `mesh_box_regions`), exact, sliver-free.
- **Curved boundary conformity at the input-surface resolution** — solved
  (`recover_boundary`): give it a finely-triangulated curved surface (`MeshSurface`)
  and the boundary is size-controlled and exactly conforming; only the *interior*
  lacks independent size control.

**Unfinished in the historical probes:** simultaneous *uniform interior* size
control and *exact* boundary conformity on curved domains. A boundary-protected
refinement strategy with robust cospherical handling is one candidate, not the
only possible algorithm. The three recorded shortcuts fell short; the active
roadmap still requires an implementation plus current independent correctness,
robustness, resource, and termination evidence.
