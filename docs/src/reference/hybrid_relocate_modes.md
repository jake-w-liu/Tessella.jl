# Pyramid relocation settings

`Mesh.OptimizePyramids` is a finite floating-point option. In the ordinary
single-Region hybrid producer, nonnegative values other than exactly `1.0`
select the existing Pyramid relocation kernel, including its normal and
golden-search passes. This matches the native branch used by `0`, `2` and
fractional values such as `0.7`; the setting does not disable relocation.

Exactly `1.0` selects the separate [Winslow untangler](hybrid_winslow.md).
Negative values retain the certified normal-offset constructor. Region
batches with separate components or shared Triangles use the same kernels through the
[hybrid batch producer](hybrid_batch_pyramids.md). All source-geometry, ownership, exterior,
typed-cell and volume certificates still apply.
