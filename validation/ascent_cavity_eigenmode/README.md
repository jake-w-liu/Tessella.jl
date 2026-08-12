# ASCENT cavity-eigenmode physics validation on a Tessella mesh

**Verified 2026-08-13.** The strongest "ready for ASCENT to use" proof: ASCENT computes a real
**physical result** — the resonant frequency of a cavity — on a Tessella-generated mesh, and it
matches the **closed-form analytic reference**. This is exactly the structure of an HFSS regression
case (geometry → mesh → solve → compare to a known reference), run end-to-end here.

A rectangular PEC cavity `a×b×d` has analytic resonances `f_mnp = (c/2)·√((m/a)²+(n/b)²+(p/d)²)`.
For `a=1, b=0.5, d=0.75 m` the dominant TE101 mode is `f = (c/2)·√(1/a²+1/d²) = 249.827 MHz`.

`generate_cavity.jl` (Tessella env) meshes the cavity (`mesh_box`, structured Kuhn, 3564 tets).
`solve_eigenmode.jl` (ASCENT env) loads it, builds the FEM cache (first-order Nedelec H(curl), PEC
walls), and runs `solve_eigenmodes`. Result:

```
lowest resonant f (FEM) = 249.711 MHz
analytic TE101          = 249.827 MHz
relative error          = 0.046 %
CAVITY_EIGENMODE_OK
```

**0.046 % error** — within first-order-Nedelec discretization on this mesh. So ASCENT solves the
Maxwell eigenvalue problem on a Tessella mesh and recovers the correct physics against a first-
principles reference.

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
