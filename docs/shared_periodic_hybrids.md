# Shared-Tri periodic hybrid components

Connected ordinary hybrid Regions can retain translated periodic boundary
copies while sharing literal Tri3 interfaces. Each Region has one declared
periodic Surface pair, and each master and slave Face is an exterior boundary
owned by that Region alone. Shared interfaces retain their actual lower Tri
UIDs, ordered nodes and opposite Region incidence.

Admission binds the complete actual copied Surface constructor, classified
source/target rows and explicit UID-pair events before any RNG or factory-clock
use. Repeated slave UIDs on neighboring periodic Faces must agree on their
current master UID and affine map. No relation is selected from coordinates;
coordinates verify the existing events. Conflicting shared copy relations
reject before production.

The connected producer retains its owned global recovery, shared-panel
multiplicity, Region partition, node/cell identity and final exterior/positive
volume certificates. Periodic constructor binding is checked again after each
final regional optimizer. The initial scope retains default
`Mesh.OptimizePyramids=0`, translations and one pair per Region; shared periodic
Faces, cross-Region pairs, incompatible copy relations, arbitrary affine maps
and general periodic component parity remain outside it.

The opt-in [public fixture](../test/interfaces/api_hybrid_shared_periodic_functional_test.jl)
uses two adjacent unit Regions with two actual shared Tri3 cells and two
translated periodic Quad4 copies. It checks independent volumes, opposite
interface sides, retained lower node UIDs/coordinates and cyclic incidence,
four primary and nine P2 periodic pairs, and complete owned MSH4 round-trip
records. Ordinary SetOrder1 cell factories may renew lower cell labels;
`Mesh.Renumber=0` preserves the fixture's literal node labels. No Volume
generation is added to default tests.

The negative fixtures cover the default Pyramid-mode boundary and an actual
mixed-direction lower-copy setup that rejects in the existing Surface planner
before a lower mesh is published. The joint copy-consistency rule remains an
admission certificate; no fabricated copy-event maps are used to exercise it.

The final bounded reader passed63/63 onJulia1.12.7, actual exit0 in36.98 s
with minimal compilation and bounds checks. It contains59 positive/P2/MSH4
checks and four actual atomic guard checks. Its first failed fixture receipt
remains separate; Source was unchanged when the lower-guard expectation was
corrected. No claimed joint-gate runtime coverage is inferred from that control.
