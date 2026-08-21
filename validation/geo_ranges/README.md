# Finite `.geo` range differential

`differential.jl` compares Tessella's bounded constant-list expansion with the
installed Gmsh 4.15.2 parser through its Julia API. It checks bit-exact Float64
values for two- and three-term `start:end[:increment]` lists, repeated-addition
rounding, descending and empty ranges, Gmsh's signed-32-bit `%` operator semantics,
integer-list truncation, and outer negation of brace lists.

The implementation was traced against Gmsh tag `gmsh_4_15_2`, commit
`657c8e915f60405e6cad0c8ec7faf812bfff1a60`, especially
`src/parser/Gmsh.y` (`ListOfDouble`, `ListOfDoubleWithBraces`, `FExpr_Multi`, and
the `FExpr '%' FExpr` action). The gate requires both the Gmsh 4.15.2 CLI and
matching Julia API and exits nonzero on a missing runtime, version mismatch, or any
value mismatch.

Run it with either supported Julia runtime:

```sh
julia +1.12 --project=. --startup-file=no --check-bounds=yes \
  validation/geo_ranges/differential.jl
julia +1.11 --project=. --startup-file=no --check-bounds=yes \
  validation/geo_ranges/differential.jl
```

This bounded scanner is not a general `.geo` interpreter. It expands finite constant
ranges only in recognized numeric field options and field selectors. Entirely numeric
Physical memberships are checked for bounded range expansion but remain opaque
geometry data. Dynamic lists, loops, macros, geometry-derived or mixed Physical
right-hand sides, option reads, logical/ternary evaluation, and CSG statements remain
explicit non-claims.
