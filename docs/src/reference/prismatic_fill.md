# Certified fill for retained rigid prisms

`import_iges(path)` and surface volume filling support bounded rigid prisms.
This producer runs after ordinary and fan fills fail, before the exact fallback.

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
and unsupported shapes leave the existing fallback in place. This path does not add general
CAD, curved-prism, or arbitrary polyhedron support.
