# Session POS list import

`API.open(path)` and `API.merge(path)` now admit `.pos` files decoded by the
existing `Post.read_pos` parser. The session contract supports all 24 plain
linear scalar/vector/tensor families on Point, Line, Triangle, Quadrangle,
Tetrahedron, Hexahedron, Prism and Pyramid. It retains literal coordinates,
components, source order within a family, named/empty views and physical times.
Imported type groups use the order observed through public Gmsh 4.15.2 queries.

Open reuses an empty current model or appends a real empty model after a
nonempty current model; its name is the filename stem and its associated
filename is the POS path. Prior models and held views remain alive. Merge adds
only list views and preserves the current model, mesh and associated filename.
Comment-only POS files create no views. Parse, admission and tag allocation
failures publish nothing; this is intentionally stronger than native POS
parsing that can publish a valid prefix from a malformed file.

T2/T3 text records import as list string entries keyed `"T2"`/`"T3"`: one
list entry per record with the record's per-step strings `'\0'`-joined
(upstream's `T2C`/`T3C` char buffers) and four coords per record — x, y, z
(0 for T2), style (upstream `T2D`/`T3D` minus the char-buffer index). Their
numeric arguments are `FExpr` expressions evaluated at parse time —
`General.GraphicsWidth/4` and friends resolve through the session `OPTIONS`
map first and the generated `DefaultOptions` table afterwards — so the t8
`view4.pos` annotations import with evaluated positions. `view.write` on a
merged view re-emits the records verbatim. Interpolation schemes and
higher-order list suffixes still reject atomically. The richer standalone
POS decoder supports those structures; this increment does not silently
discard them or claim the full POS language.
Names use the actual session view for both public Name options and writer
headers. Imported `NbTimeStep` counts data steps, so empty views report zero
even if they declare TIME, matching a tiny public SDK observation.

Session list storage now retains a time vector. Existing caller-owned list
APIs keep implicit default times. POS, CSV and VTK writing preserve imported
times, including a single nondefault time and explicit times on empty views.
The POS serializer retains its previous output for an ordinary implicit
single-step time of zero. Existing PostView fields still select integer step
indices. [Session time metadata](view_time.md) now derives actual step counts
and resolves displayed TimeStep/Time selection. A tiny native Line2 field
control keeps the same first-step output despite changing display selection;
this observation does not establish broader field-context parity.

Qualification is a bounded file-only reader with no mesh generation. It checks
real session open/merge behavior, list values/detachment, names, time exports,
held binary model-view continuity, the 24-family public witness and atomic
rejection. Two default smoke assertions import one point record and merge an
empty file. The final qualification packet preserves exact Source and receipts.

The final bounded reader passed 30/30 (28 functional checks plus two default
smoke assertions), Julia 1.13.1 with `--compile=min -O0 --check-bounds=yes` and
one thread, actual exit zero, warm wall 3.9672326 seconds. The failed first
receipt remains preserved: its sole failed assertion incorrectly expected
opening the named empty current model to remove the prior initial empty slot.
The corrected fixture retains prior slots; implementation Source is unchanged.

Valid [list/model representation changes](view_representations.md) replace
inactive storage atomically. Switching through model data starts fresh implicit
list times; updating a list within its existing representation keeps imported TIME.
