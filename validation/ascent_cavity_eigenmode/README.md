# ASCENT cavity-eigenmode physics validation on a Tessella mesh

Status reviewed **2026-10-10**: this page describes the dated ASCENT cavity capture below.
The authoritative [handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
separates published code and historical results from **Root V24** and
**Release V6**, which remain unqualified combined candidates. No rerun
is claimed here. Tessella supports Julia **1.12.x and 1.13.x**; follow
the [validation launch guidance](../README.md#run) for both runtimes.

**Historical capture verified 2026-08-13.** ASCENT computed cavity resonant
frequencies on a Tessella-generated mesh and compared them with the closed-form
analytic reference. This captured case exercised geometry → mesh → solve → compare.
It has not been rerun on the current integration candidates.

A rectangular PEC cavity `a×b×d` has analytic resonances `f_mnp = (c/2)·√((m/a)²+(n/b)²+(p/d)²)`.
For `a=1, b=0.5, d=0.75 m` the dominant TE101 mode is `f = (c/2)·√(1/a²+1/d²) = 249.827 MHz`.

`generate_cavity.jl` (Tessella env) meshes the cavity (`mesh_box`, structured Kuhn, 3564 tets).
`solve_eigenmode.jl` (ASCENT env) loads it, builds the FEM cache (first-order Nedelec H(curl), PEC
walls), runs `solve_eigenmodes`, and matches each of the first 5 FEM modes to its nearest analytic
resonance (robust to degeneracies). Result:

```
FEM mode 1 = 249.711 MHz → analytic 249.827 MHz  (err 0.046 %)
FEM mode 2 = 334.435 MHz → analytic 335.178 MHz  (err 0.222 %)
FEM mode 3 = 359.241 MHz → analytic 360.306 MHz  (err 0.296 %)
FEM mode 4 = 359.814 MHz → analytic 360.306 MHz  (err 0.137 %)   # the degenerate pair, resolved
FEM mode 5 = 390.305 MHz → analytic 390.242 MHz  (err 0.016 %)
CAVITY_EIGENMODE_OK  (max err 0.296 %)
```

**The whole first-5-mode spectrum matches to <0.3 %** — within first-order-Nedelec discretization on
this mesh, including the correctly-resolved near-degenerate pair at ~360 MHz. So ASCENT solves the
Maxwell eigenvalue problem on a Tessella mesh and recovers the correct physical **spectrum** against
a first-principles reference.

This is a *complete physics regression case* on a Tessella mesh with an independent (analytic)
oracle — the same shape as an HFSS UserGuide cavity example. The literal 22 HFSS cases run this
pipeline on the guide's specific antenna/microwave geometries (OCC-built, proprietary reference
data) — the remaining external campaign; the dated capture demonstrates the mesh→solve→validate loop for this cavity.

## Reproduce

```
julia +1.12 -O2 --startup-file=no --history-file=no --check-bounds=yes --threads=1 --gcthreads=1 --heap-size-hint=2G --project=/path/to/Tessella.jl validation/ascent_cavity_eigenmode/generate_cavity.jl
julia --project=/path/to/ASCENT    validation/ascent_cavity_eigenmode/solve_eigenmode.jl
```

Use the Julia version and dependencies required by the external ASCENT project
for the solve command; Tessella's 1.12/1.13 support does not qualify either ASCENT
environment. The generated `cavity.msh` is a git-ignored build artifact.
