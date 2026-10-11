# Certified fill for retained rigid prisms

`import_iges(path)` and the existing surface volume-fill front door can now fill
a bounded rigid prism after ordinary and fan fills fail, before the existing
exact fallback. The retained arbitrarily rotated Plane190 through-hole fixture
uses 16 source vertices, 32 source triangles, and 24 positive tetrahedra. Its
material volume is 0.75 and its through-hole remains empty.

Recognition requires at most 128 source vertices, two separated planar caps,
and a unique one-to-one source vertex pairing by one extrusion vector. The two
caps must have identical triangle supports under that literal row pairing.
Every cap border edge must have exactly the two corresponding sidewall
triangles. Their shared diagonals impose a single acyclic precedence on the
base vertices. A deterministic topological order supplies three staircase
tetrahedra per cap triangle, with consistent diagonals across adjacent prisms.

Coordinate tolerance is used only to recognize the pairing. The fill retains
the original coordinates and source row references: it does not flatten,
snap, or replace vertices. The mesh constructor owns the output storage.
Acceptance requires positive mesh validation and the unchanged full surface
fill certificate, including the actual input boundary and PLC conformity.

Ambiguous pairing, different cap triangulations, cyclic sidewall precedence,
and unsupported shapes leave the existing fallback in place. In particular,
the retained quarter-turn through-hole fixture has different cap triangle
supports and correctly uses that fallback. This path does not add general
CAD, curved-prism, or arbitrary polyhedron support.

The opt-in regression is
[`prismatic_fill_functional_test.jl`](../test/meshing/prismatic_fill_functional_test.jl).
It checks the actual rotated IGES public fill, literal boundary supports,
positive analytic volume, hole exclusion, output storage independence,
translated/outer-box/triangular-prism controls, the unequal-cap public fallback,
and cheap unsupported/cyclic recognition cases. It performs no native calls.
