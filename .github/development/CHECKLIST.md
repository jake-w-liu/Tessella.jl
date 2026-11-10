# Implementation and release checklist

This is the current requirement index, not a declaration that every row is
unfinished or complete. Preserve the whole target while qualifying each
subdomain against current source and fresh evidence. Detailed subitems and
historical receipts remain in the
[complete previous ledger](https://github.com/jake-w-liu/Tessella.jl/blob/fbf34516da30408248e21a592ea599fdb9736b5c/.github/development/history/DEVIN_HANDOFF.md#complete-remaining-checklist)
and [P1–P6 plan](https://github.com/jake-w-liu/Tessella.jl/blob/fbf34516da30408248e21a592ea599fdb9736b5c/.github/development/history/PLAN.md#gmsh-parity-tracks). Dates and prior rejection
classifications do not establish the current implementation state.

## Required domains

| ID | Completion requirement |
| --- | --- |
| CAD01 | Nested product-usage DAGs, repeated paths, dimensional transforms, reachability/cycle/orphan checks, and separate occurrence ownership. |
| CAD02 | Mixed mapped/product workflows, multiple roots and solids per leaf, geometry plus subassembly usages, and general occurrence transforms with actual context/item provenance. |
| CAD03 | Ordinary and mapped-root unit conversion, conversion-based length/angular units, compatible uncertainty aliases, and target-unit contexts. |
| CAD04 | Curved trims, periodic/general NURBS topology, nonlinear UV charts, and complete imported analytic Face/Curve parameterization and registry. Mesh import alone is not a CAD kernel. |
| CAD05 | Contact/touching/overlapping assemblies with genuine interfaces; arbitrary CAD transformations, Booleans, imprints, and intersections beyond supported materialized primitives. |
| GEO01 | Remaining `.geo` language and genuine `model.geo`/`model.occ` frontend coverage, allocators, synchronization, lifecycle, and full CAD execution. Reuse existing Model operations. |
| PER01 | Proper-rigid rotated hybrid copies with certified derived cross-side Curve transport. |
| PER02 | Multiple periodic directions/pairs, general affine and cross-Region/shared Face contexts, consistent corner/edge relations, and embedding combinations. |
| EMB01 | Larger, subdivided, interior-vertex, holed, multiple, and curved sheets; supported mixed sheet families; nonlinear Curves and broader multi-Region embedded composition. |
| ALG01 | All successful ordinary surface/volume algorithms and mixed-block unstructured integration/recombination. Enumerate dispatch requirements instead of removing unsupported selectors. |
| MET01 | Native anisotropic algorithm-7/BAMG functionality. The explicit conforming tensor-bisection API is a distinct implemented algorithm and must not be renamed BAMG. |
| MET02 | Variable-metric Volume generation; tensor boundary adaptation and curved, periodic, recombined, and embedded Surface contexts with actual lower ownership transport. |
| MET03 | Full native automatic sizing/octree behavior, including curvature, gaps, bulk sizing, and gradation; demonstrate any HXT/p4est equivalence claims. |
| LAY01 | General NoNew source grids/categories and neighboring/shared-Region propagation; periodic, compound, reversed, constrained, and nonlinear retained sweep sources. |
| LAY02 | Boundary-layer Extrude and pipe frontends beyond ordinary Layers; remaining curved-layer and recombination contexts. |
| HO01 | Curved/mixed/hybrid high-order generation, projection/refinement placement, and global Jacobian certification across every required family/order. |
| OPT01 | Mixed Surface/nonlinear/primitive relocation, general high-order/hybrid optimization, and independent Netgen-equivalent behavior. |
| PART01 | Actual Point/Curve classification in partition children/interfaces/ghosts; high-order cells, periodic/parametric topology, triple junctions, and distinct coincident supports. |
| PAR01 | Actual parallel producer and partition workflows with deterministic ownership and output. Parallel query slicing alone is insufficient. |
| IO01 | True model filename import/merge/export and long-tail formats; bridge standalone STL/CAD readers into real session lifecycles. |
| IO02 | Lossless higher/incomplete-pyramid and remaining high-order target mappings where the formats permit them; pyramid model-view interpolation beyond P2. |
| VIEW01 | POS text/interpolation/high-order session storage, query, and writer transport. |
| VIEW02 | Homogeneous model data, list strings, interpolation matrices, aliases, combine/probe, per-window visibility, real owned data, and time semantics. |
| CLI01 | General conversion, partition, adaptation, and post-processing action frontends. |
| UI01 | Plugin option/run catalog, parser/ONELAB/graphics/FLTK services, and actual interactive GUI/post controls beyond headless command state. |
| MATCH01 | Remaining native matching/optimizer ordering and ties through independent algorithms and opaque public evidence. Preserve optimum certificates and seeded odd-cycle fixes; never inspect restricted matching implementations. |
| P601 | Requirement-by-requirement tutorial/API/option/format conformance across P1–P6, including quality, basis, queries, and lifecycle contexts. Record every missing endpoint or guard. |

CAD01–CAD03 and PER01 have historical implementation receipts. Other rows also
contain substantial completed domains. Preserve these implementations and
verify current behavior before classifying an old blocker as a missing feature.
Conversely, a classified rejection does not satisfy a required successful
upstream context.

## Residuals requiring current verification

Current source audits retain the following areas for investigation;
the full per-row ledger remains authoritative for additional subitems:

- Broader closed-Face native mean-plane frame contexts. The direct-model
  single-closed-spline frame gap is fixed in the current increment; see the
  fresh evidence in [development notes](README.md). Recheck historical
  minimal-compiler GmshLibm diagnoses against current source.
- Natural finite sphere/torus/B-spline/plane faces, closed natural solids,
  crossed-cylinder seams/fuse images, centroidal CAD measures and the wedge
  orientation are now implemented and default-wired. Native pcurve-only edges
  also retain their exact composed geometry, derivatives, affine images and
  trimmed planar witnesses. Natural-solid moments retain a nearby frame and
  rational faces use bounded adaptive quadrature. Remaining BREP families
  include hyperbolas/general offsets and internal/external orientations;
  general quadric-bounded shell inertia still rejects rather than approximating.
- Full-quad and anisotropic integration with embedded, transfinite, smoothing,
  periodic/compound, and curved/scalar-field contexts.
- NURBS trim/topology and general CAD Boolean/fillet contexts beyond the
  recently verified tutorial and coincident/tangent fixes.
- Parallel eligibility beyond the established detached producer classes,
  including RTREE recombined boundary write-back; stateful callbacks retain
  their documented serial contract.
- Full GUI/session services and native matching tie identity. A certified
  optimum proves cost correctness, not upstream chronological tie identity.
- Systematic option/format/API conformance, optimization/resource measurements,
  and a fresh independent audit after the final fix.

The two original designs remain in Git history:
[nested STEP assemblies](https://github.com/jake-w-liu/Tessella.jl/blob/fbf34516da30408248e21a592ea599fdb9736b5c/.github/development/history/designs/nested_step_assemblies.md) and
[meshing/periodic transport](https://github.com/jake-w-liu/Tessella.jl/blob/fbf34516da30408248e21a592ea599fdb9736b5c/.github/development/history/designs/meshing_remaining.md). Their initial
increments were subsequently implemented; follow-on requirements remain in
the domain rows above. Do not restart completed designs from their old prose.

## Version 0.1.0 release

- [ ] Complete all genuine required implementation work and fresh bug audits.
- [x] Support correctness/resource claims with evidence of the same scope.
- [x] Keep only README Markdown among tracked root files, retain meaningful API reference
  material, and preserve third-party notices and historical evidence.
- [x] Set version 0.1.0, declare dependency compatibility, retain the chosen MIT
  license and third-party notices, and enable supported-runtime CI.
- [ ] Run final ordinary `Pkg.test()` on the final source on Julia 1.12 and 1.13; separately record
  any optimized compiler checks. Push verified source normally to `main` and
  verify the exact remote SHA.
- [ ] Disclose substantial AI assistance and verify the maintainer's human
  review/understanding required by General's registration policy. AI assistance
  is disclosed in the root README; human review is not yet verified.
- [ ] Trigger Registrator for the final verified commit, resolve its actual
  checks, and verify the General registration merge and package installation.
  A version bump, comment, pull request, or tag alone is not registration.

Registration policy: [General](https://github.com/JuliaRegistries/General#registering-a-package-in-general),
[LLM contributions](https://github.com/JuliaRegistries/General#is-there-any-policy-regarding-the-use-of-llms-in-registered-packages),
[Registrator](https://github.com/JuliaRegistries/Registrator.jl#via-the-github-app),
and [AutoMerge checks](https://juliaregistries.github.io/RegistryCI.jl/stable/guidelines/).
