# Licenses and source provenance

Tessella contains Julia translations and adaptations of upstream source code,
including Gmsh and OCCT routines. These are source derivations, even though
Tessella does not load or execute either project's meshing or CAD library.
Compatibility tests and manual-based reimplementations do not remove the
license obligations of the translated portions.

The combined package is distributed under **GPL-3.0-or-later**; see
[the package license](../LICENSE). Independently authored Tessella material
continues to carry the [MIT grant](NOTICES.md#tessella-independent-material) when used separately. The original
MIT notice did not grant rights to third-party material under MIT, and these
corrections do not revoke an existing valid MIT grant. Component licenses
apply to their respective portions; an `AND` expression in a mixed source
file records those component obligations, rather than offering a choice of
licenses for the whole file.

Gmsh-derived portions require GPL-2.0-or-later. The libstdc++ source adaptation
in `Mesh2DGmshInitialOrder.jl` is GPL-3.0-or-later with the GCC Runtime Library
Exception. Its copied prime table and hashtable policy are retained. That
exception permits specified combinations of compiled target code; it is not
a general permission to redistribute translated library source under MIT or
LGPL. GPL version 3 is therefore the common license for this distribution.
Simply choosing LGPL would not satisfy the Gmsh and libstdc++ source terms.
No additional restriction is imposed on the independently authored material.

OCCT 7.9.3 source is LGPL-2.1-only with the OCCT exception. Its notices and
exception are retained. LGPL 2.1 section 3 permits applying ordinary GPL terms
to a copy in a combined GPL work. The exception concerns object code using
OCCT headers; it does not place translated classifier source under MIT.
The Gmsh linking exception is retained with its original scope in
[NOTICES.md#gmsh](NOTICES.md#gmsh); it does not relicense other components.

## Source inventory

The following source relationships are recorded in the Julia files themselves
and were checked against the named releases. Other comments that mention an
upstream routine describe compatibility context and are not, by that fact
alone, evidence of a source translation. This inventory identifies confirmed
source adaptations; it is not a claim that every other file was developed
without consulting upstream source.

| Tessella source | Upstream source and retained terms |
| --- | --- |
| `src/core/Elements.jl` | Gmsh 4.15.2 `src/numeric/pointsGenerators.cpp`: local-node monomial lattices; GPL-2.0-or-later. |
| `src/core/HierarchicalBases.jl` | Gmsh 4.15.2 `src/numeric/OrthogonalPoly.cpp` and hierarchical basis routines: polynomial evaluation and orientation; GPL-2.0-or-later. |
| `src/meshing/Mesh1D.jl` | Gmsh 4.15.2 `src/mesh/meshGEdge.cpp`: `filterPoints`; GPL-2.0-or-later. |
| `src/geometry/ModelGmshCurveIntegration.jl` | Gmsh 4.15.2 `src/mesh/meshGEdge.cpp`: curve integration and smoothing; GPL-2.0-or-later. |
| `src/interfaces/APIPlugins.jl` | Gmsh 4.15.2 `src/geo/MEdge.cpp`: `SortEdgeConsecutive`, alongside plugin compatibility routines; GPL-2.0-or-later. |
| `src/geometry/ModelDiscreteCurveGeometry.jl` | Gmsh 4.15.2 `src/geo/MEdge.cpp`: `SortEdgeConsecutive`; GPL-2.0-or-later. |
| `src/geometry/ModelCurved.jl` | Gmsh 4.15.2 `src/geo/Geo.cpp` `EndCurve` and parameter-search routines; GPL-2.0-or-later. OCCT 7.9.3 geometric operations, extrema/projection, `BRepClass_FaceClassifier` chain and associated 2-D transition/intersection routines; LGPL-2.1-only WITH OCCT-exception-1.0. |
| `src/meshing/Mesh3DGmshOptimize.jl` | Gmsh 4.15.2 `src/mesh/meshGRegionLocalMeshMod.cpp`: ring-swap pattern ordering; GPL-2.0-or-later. The upstream SHA-256 is `abfb4c8cbda1356aeaf8cf017eea2ebf103acbbd2a22664a49664c0157bec658`. |
| `src/geometry/ModelNativeTetgenArena.jl` | Gmsh 4.15.2 `src/mesh/tetgenBR.cxx` and `tetgenBR.h`: topology tables and bounded flip routines. The Gmsh-distributed Tetgen/BR agreement notice is retained in [NOTICES.md#tetgenbr](NOTICES.md#tetgenbr); its referenced Gmsh GPL terms are in [NOTICES.md#gmsh](NOTICES.md#gmsh). |
| `src/geometry/ModelSplines.jl` | Gmsh 4.15.2 `src/geo/GeoInterpolation.cpp` and `GEdge` projection routines: spline-family evaluation, derivatives and projection; GPL-2.0-or-later. |
| `src/interfaces/GeoOptionTables.jl` | Gmsh 4.15.2 `src/common/DefaultOptions.h` and `Options.cpp`: mechanically generated option tables; GPL-2.0-or-later. |
| `src/interfaces/GeoColorNames.jl` | Gmsh 4.15.2 `src/common/Colors.h`: generated `ColorString` table; GPL-2.0-or-later. |
| `src/interfaces/IO.jl` | Gmsh 4.15.2 `src/parser/Gmsh.l`: verbatim lexer keyword table, alongside `DefaultOptions.h` option data; GPL-2.0-or-later. |
| `src/meshing/BoundaryLayer.jl` | Gmsh 4.15.2 `src/mesh/meshGFace.cpp` and `boundaryLayersData.cpp`: `fanType` cell emission loops; GPL-2.0-or-later. |
| `src/structured/TransfiniteTriangle.jl` | Gmsh 4.15.2 `src/mesh/meshGFaceTransfinite.cpp`: triangle and collapsed-grid node/cell emission; GPL-2.0-or-later. |
| `src/structured/TransfiniteVolume.jl` | Gmsh 4.15.2 `src/mesh/meshGRegionTransfinite.cpp`: `transfiniteHex` expression ordering; GPL-2.0-or-later. |
| `src/structured/TransfinitePrism.jl` | Gmsh 4.15.2 `src/mesh/meshGRegionTransfinite.cpp` and transfinite3 surface emission adapted onto prism slots; GPL-2.0-or-later. |
| `src/interfaces/APICellPlacement.jl` | Gmsh 4.15.2 `src/numeric/InnerVertexPlacement.cpp`: boundary-to-interior placement routines; GPL-2.0-or-later. |
| `src/meshing/Recombine3D.jl` | HXT in Gmsh 4.15.2 `contrib/hxt/Combine/src/hxt_quality.cpp`: Float64 ports of `computeJacobian`, `computeBezierCoefficients`, `hexExactValidity` and `hexSubdivisionBasedQuality`. Upstream authors: Amaury Johnen, Kilian Verhetsel, Jeanne Pellerin. GPL-2.0-or-later; see [HXT's grant and exception](NOTICES.md#hxt). The surrounding independently implemented recombination retains MIT. |
| `src/interfaces/APINativeLookup.jl` | Gmsh 4.15.2 `GModel` tag-query/cache paths: `getMeshElementByTag`, node lookup, cache rebuilding and error state; GPL-2.0-or-later. |
| `src/meshing/Mesh2DGmshInitialOrder.jl` | GCC 11.3.0 libstdc++ `include/bits/hashtable_policy.h`, `include/bits/hashtable.h`, `src/c++11/hashtable_c++0x.cc` and `src/shared/hashtable-aux.cc`: prime rehash policy, prime table and unique-key linked-list order; GPL-3.0-or-later WITH GCC-exception-3.1. |
| `src/meshing/WeightedMatching.jl` | NetworkX 3.6.1 weighted matching; BSD-3-Clause. See [the matching provenance](#weighted-matching-provenance) and [full notice](NOTICES.md#networkx). |
| `src/meshing/GmshPairSort.jl` | Original HP/SGI STL sorting and heap routines under their permissive notices; see [NOTICES.md#hpsgi-stl](NOTICES.md#hpsgi-stl). |
| `src/core/GmshLibm.jl` | MinGW-w64 native Windows math adaptations; ZPL-2.1 and public-domain portions, detailed in [NOTICES.md#mingw-w64](NOTICES.md#mingw-w64). |

Gmsh sources carry Copyright (C) 1997–2025 C. Geuzaine, J.-F. Remacle.
The OCCT sources used here carry Copyright (c) 1991–1999 Matra Datavision
(the exact initial year varies by upstream file) and Copyright (c) 1999–2014
OPEN CASCADE SAS. The libstdc++ sources carry Copyright (C) 2007–2021
Free Software Foundation, Inc. (the initial year varies by upstream file).
Gmsh's Tetgen/BR distribution carries Copyright 2015–2016 Weierstrass
Institute for Applied Analysis and Stochastics. Source-file notices preserve
these attributions alongside Tessella's 2026 Julia adaptation and changes.

## Upstream releases and full terms

- [Gmsh 4.15.2 sources](https://github.com/live-clones/gmsh/tree/gmsh_4_15_2),
  [official license](https://gmsh.info/LICENSE.txt),
  [retained license and exception](NOTICES.md#gmsh), and
  [GPL version 2](GPL-2.0-or-later.txt).
- [HXT as distributed with Gmsh 4.15.2](https://github.com/live-clones/gmsh/tree/gmsh_4_15_2/contrib/hxt),
  [grant and Gmsh combination exception](NOTICES.md#hxt).
- [OCCT 7.9.3 sources](https://github.com/Open-Cascade-SAS/OCCT/tree/V7_9_3),
  [LGPL 2.1](LGPL-2.1-only.txt), and
  [OCCT exception 1.0](NOTICES.md#occt-exception).
- [GCC 11.3.0 sources](https://github.com/gcc-mirror/gcc/tree/releases/gcc-11.3.0/libstdc%2B%2B-v3),
  [GPL version 3](../LICENSE), and
  [GCC Runtime Library Exception 3.1](NOTICES.md#gcc-runtime-library-exception).

The complete third-party notices and license texts must accompany redistributed
copies as their terms require. Generated mesh data is not relicensed merely
because Tessella generated it; the GPL governs output only when its contents
constitute a covered work.

## Weighted matching provenance

`WeightedMatching.jl` adapts the Galil/Edmonds primal-dual algorithm and
bookkeeping in NetworkX 3.6.1 `max_weight_matching`, using indexed state,
bounded reusable blossom slots and workspaces.

- [Upstream source](https://github.com/networkx/networkx/blob/networkx-3.6.1/networkx/algorithms/matching.py)
- Source SHA-256: `5e08becf2f42da936f1b04e2c53155eeeb7b5be97cb3bcd480985af6ff03d2cf`
- [BSD-3-Clause notice](NOTICES.md#networkx)

The minimum-cost graph interpretation follows Gmsh 4.15.2
`meshGFaceOptimize.cpp`. Production does not execute Gmsh or NetworkX.
`OrderedMatching.jl` and `OrderedFractionalMatching.jl` independently implement
singleton-dual and alternating-tree searches, with primal/dual certificates.
`SeededWeightedMatching.jl` continues fractional singleton state through the
BSD weighted solver. Minimum-cost correctness does not establish identical
chronological tie choices to Gmsh.

Production code does not incorporate Concorde or Blossom IV. Gmsh's
`contrib/blossom/README.txt` grants its special license only within the Gmsh system.
