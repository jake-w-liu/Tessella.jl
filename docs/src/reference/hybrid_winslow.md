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
certificates still apply. [Region batches](hybrid_batch_pyramids.md)
also support this optimizer with their continuing global UID factory.

Exact Gmsh multi-apex visitation order and platform-specific pointer ordering are not guaranteed.
