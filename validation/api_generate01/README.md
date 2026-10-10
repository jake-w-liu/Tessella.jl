# API generation 0/1 differential

Status reviewed **2026-10-10**: this page describes the bounded 0D/1D differential and its recorded artifact evidence.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

Run from the repository root with the installed Gmsh 4.15.2 Julia binding and
runtime DLL:

```powershell
$env:GMSH_JULIA_API = 'C:/path/to/gmsh.jl'
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=. validation/api_generate01/differential.jl
```

The driver checks the binding and loaded DLL versions. A successful run prints
`API_GENERATE01_DIFFERENTIAL_OK` with fixture, stage, and explicit blocker counts.
`--compiled-modules=no --compile=min` can be added for bounded debugging;
such runs do not qualify the normal-compilation release.

Coverage includes fresh dimension0, dimension1 Point15/Line cells, actual raw
Point cell selection, independent coincident entities and distinct coincident
raw node tags, globally shared record connectivity, sparse nonmonotonic public
labels, ownership and query ordering. Dimension0 preserves the actual existing
cell set; it does not synthesize Point15 cells. Dimension1 retains discrete
higher-dimensional blocks and removes native higher-dimensional cells.

Fixtures exercise ElementOrder versus immediate setOrder, native/discrete P1/P2
behavior, raw-only immediate order changes, OnlyEmpty/Visible, physical-group
renumber priority, grading laws, field clipping and size factor, callbacks,
periodic primary/full-P2 nodes and seven-output periodic keys. Lifecycle checks
include allocator high-water across clear and model switching, automatic tags,
refinement, duplicate survivors, explicit/default renumbering, authoritative
record suppression, visibility persistence, homology and endpoint lookup.
Periodic order reduction and partial clear verify that published node relations
contain only live pairs; invalid clear scopes must leave model and cache intact.
An unrelated hidden unmeshed Circle does not block refinement of a straight Line
whose actual mesh and boundary entities satisfy the linear CAD requirement.

Tags, typed oriented connectivity, node/cell ownership and query array order are
compared exactly. Stored raw fixture coordinates are exact. Ordinary coordinate
comparisons use `atol=5e-12, rtol=5e-12`; endpoint inverse coordinates use1e-10.
The three Progression/Bump/Beta grading fixtures and the periodic master
grading fixture retain their historical `atol=3e-7` comparison with Gmsh's
sampled adaptive density primitive. Native nonuniform Lines now use that same
density-integration route. Independent saved-coordinate unit tests enforce
`atol=2e-11, rtol=0` for the three laws and retain `atol=3e-12` for the periodic
master. Standalone analytic parameter helpers keep their separate exact-law
tests. Periodic pair dictionaries are exact; raw
pair array order depends on Gmsh's maps keyed by allocated MVertex pointers.

The published-baseline driver records a precise atomic source-state blocker: legacy
native dimension2/3 caches with nonempty attached mesh records do not retain
the allocation history needed for a lossless transition to dimension1. The
attached Point111/Point15cell777 cube fixture proves upstream node tags111:130
and allocator maxima130/846, and verifies that native rejection leaves all
source state intact. This fixture is counted separately in the final marker.

The driver does not claim coverage for native raw Line8 meshes lacking actual
Point boundary meshes: four such pinned Gmsh setter fixtures crashed the DLL in
the independent review. Fully discrete P2 and native raw P2 with actual Point
records are safe verified fixtures. This gap is separate from mixed Line1/Line8
families on one entity and from the legacy source-history blocker. This driver's
historical curved-CAD guard is not a current blanket support statement: newer
curved P2 placement/refinement candidates have scoped Source proofs, with their
whole public qualification still pending in the handoff.

`test/interfaces/api_generate01_artifacts.jl` builds nine compact public API
fixtures pinned in `test/artifacts/api_generate01_crc.txt`. Their `mixed_crc`
digests include coordinate bits, oriented typed connectivity, physical names,
entity ownership, public labels, parameters and periodic metadata. Allocator
maxima are pinned alongside each digest. The recorded nine records matched bitwise under
normal Julia 1.12.7 and 1.13.1 with bounds checks. Graded artifact hashes pin the
native density-integration samples, independently checked against Gmsh 4.15.2
coordinates. Exactly the three graded hashes changed for the native Line fix;
the other six records remain byte-identical. The differential retains the
historical bound stated above.
