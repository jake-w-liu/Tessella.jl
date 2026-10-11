# Winslow Pyramid optimization

For the certified ordinary single-Region hybrid producer,
`Mesh.OptimizePyramids=1.0` runs regularized Winslow optimization on each
actual Pyramid-apex star. The objective combines shape and volume distortion
over the incident Tets and four Pyramid decompositions. Its analytic gradient
drives a constrained L-BFGS solve with a strong-Wolfe line search.

Each patch must have a connected, closed spherical boundary and exactly one
interior apex. The optimizer visits actual allocated apex UIDs in ascending
order and locks all other node positions. Its output must improve the patch
objective and leave every decomposition positive. Existing typed-cell,
ownership, exterior-incidence, volume, embedded-Point and periodic transport
certificates still apply. Nonzero optimization for multi-Region batches
remains unsupported.

The opt-in regression
`test/interfaces/api_hybrid_winslow_functional_test.jl` checks independent
finite-difference gradients, a simplex with a known optimum, and tiny cubes
with one and six Pyramid apices. The one-apex position is checked against an
independent Gmsh 4.15.2 observation. Both public fixtures retain their lower
node UIDs and positions, owned ordered boundary cells, positive typed volume
of one, conforming exterior and MSH4 UID/owner transport.

The 44 assertions passed in Julia 1.13 with bounds checks and
`--compile=min -O0`; the warm process took 18.41 seconds. These volume fixtures
remain outside the routine default suite. Platform-specific native pointer
visitation and exact multi-apex native chronology require separate validation.
