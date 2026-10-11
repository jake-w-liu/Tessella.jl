# Negative pyramid factors

`API.option("Mesh.OptimizePyramids", value)` retains a finite floating-point
value. For one ordinary Delaunay Region with a certified convex shell and
planar Quad boundaries, a negative value places each Pyramid apex at the Quad
centroid plus `abs(value) * diagonal_length * unit_normal`. The diagonal joins
the first and third actual Quad vertices; the normal follows the actual Face
winding. Every apex must be strictly inside the shell. The producer preserves
its allocated node UID and position through insertion and Tet optimization,
and skips Pyramid relocation for this negative setting.

The volume still must pass the existing typed-cell, exterior-incidence,
ownership and analytic-volume certificates. A selected embedded Point keeps
its actual lower UID, position, owner and conforming Tet star. Invalid winding,
exterior apices and unsupported modes fail before publishing a replacement
mesh. Nonnegative settings other than exactly `1.0` use the existing Pyramid
relocation path. Exactly `1.0` requests the unsupported Winslow untangler;
nonzero settings for multi-Region batches remain unsupported.

The bounded opt-in regression is
`test/interfaces/api_hybrid_normal_apex_functional_test.jl`. It uses one cube
with an interior Point, analytic apex and cell-volume checks, an MSH4 round
trip, and small atomic-rejection checks. This scope does not establish curved
Quad support, Winslow optimization or public refinement behavior.
