# Pyramid options in GEO and CLI calls

`Mesh.OptimizePyramids` is captured for each generation call. API generation
uses its option snapshot; `execute_geo` and CLI generation use the numeric
option context of the current `.geo` program. A local negative factor stays
floating-point and reaches the hybrid constructor and its certificates. The
operation neither inherits nor changes an unrelated API option or cached mesh.
The constructor, optimizer and [batch producer](hybrid_batch_pyramids.md) use the same scoped
value.

The opt-in regression `test/interfaces/geo_hybrid_pyramid_options_functional_test.jl`
meshes one small explicit cube through `execute_geo` and CLI `-3`, with local
`-0.1` while the API retains `2.0`. It checks analytic apex positions, positive
typed cells, volume, exterior and MSH4 UID/ownership transport. An oversized
local negative factor rejects because its apices lie outside the volume,
while retaining the API mesh and the prior CLI output file. A local value of
`1.0` selects the [Winslow untangler](hybrid_winslow.md).
