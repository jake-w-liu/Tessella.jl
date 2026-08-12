# ASCENT cavity-eigenmode physics validation on a Tessella mesh

**Verified 2026-08-13.** The strongest "ready for ASCENT to use" proof: ASCENT computes a real
**physical result** — the resonant frequency of a cavity — on a Tessella-generated mesh, and it
matches the **closed-form analytic reference**. This is exactly the structure of an HFSS regression
case (geometry → mesh → solve → compare to a known reference), run end-to-end here.

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
data) — the remaining external campaign; this proves the mesh→solve→validate loop it is built on.

## Reproduce

```
julia --project=<Tessella.jl>        validation/ascent_cavity_eigenmode/generate_cavity.jl
julia --project=<2026_066/ASCENT>    validation/ascent_cavity_eigenmode/solve_eigenmode.jl
```

The generated `cavity.msh` is a git-ignored build artifact.
