# Pyramid relocation settings

`Mesh.OptimizePyramids` is a finite floating-point option. In the ordinary
single-Region hybrid producer, nonnegative values other than exactly `1.0`
select the existing Pyramid relocation kernel, including its normal and
golden-search passes. This matches the native branch used by `0`, `2` and
fractional values such as `0.7`; the setting does not disable relocation.

Exactly `1.0` selects the separate [Winslow untangler](hybrid_winslow.md).
Negative values retain the certified normal-offset constructor. Nonzero
multi-Region settings remain unsupported. All Source, ownership, exterior,
typed-cell and volume certificates still apply.

The opt-in regression `test/interfaces/api_hybrid_relocate_modes_functional_test.jl`
uses the same tiny cube with `0`, `2` and `0.7`, proves that centroid births
move into the volume, checks positive cells and analytic volume, and compares
the actual node and element identities, coordinates and connectivity. The
negative-factor regression separately preserves the embedded-Point coverage.
