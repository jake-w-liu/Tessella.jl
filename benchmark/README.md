# Meshing benchmarks

## Raw point-cloud benchmark

This benchmark compares Tessella's Julia mesh kernels with Gmsh's public C API
for planar triangulation and spatial tetrahedralization. Results apply to raw
point-cloud meshing.

On the measured shared Apple M3 host, Tessella was about 1.8–2.3 times faster
for hull-constrained 2D clouds. For raw 3D clouds, Tessella was faster at 16
points; Gmsh HXT was about 2–2.5 times faster at 64–1024 points.
The [retained results](results/apple_m3_julia_1_13_4_gmsh_4_15_2.toml) use Julia
1.13.4 and Gmsh 4.15.2 and include every configuration, timing sample, source
hash, and four small geometry fixtures.

Times are microseconds: the median of three input-seed medians, each from nine
interleaved repetitions. Both implementations passed the exact hull-completion
checks in these table groups. Empty circumspheres were checked exhaustively at
16 points.

| Task | Points | Gmsh (µs) | Tessella (µs) |
| --- | ---: | ---: | ---: |
| 2D, hull constraints | 16 | 55.83 | 30.88 |
| 2D, hull constraints | 64 | 265.42 | 128.79 |
| 2D, hull constraints | 256 | 1692.88 | 744.54 |
| 2D, hull constraints | 1024 | 8501.75 | 3992.96 |
| 3D, HXT | 16 | 66.04 | 30.71 |
| 3D, HXT | 64 | 195.62 | 389.67 |
| 3D, HXT | 256 | 837.88 | 2121.38 |
| 3D, HXT | 1024 | 4174.54 | 9779.67 |

The small constrained 2D Gmsh outputs filled the hull but failed the exact
empty-circle criterion with both default and zero random factors. Its raw
triangulation documentation does not explicitly promise that criterion.
Unconstrained default-algorithm outputs that failed this benchmark's convex-hull
completion checks remain in the data; HXT is the primary 3D comparison.

Install an external [Gmsh 4.15.2 SDK](https://gmsh.info/#Download). Gmsh is optional
and is not a Tessella dependency. From the repository root:

```sh
export GMSH_LIBRARY=/path/to/sdk/lib/libgmsh.so
julia --project=. --startup-file=no -O2 benchmark/raw_pointcloud.jl --check-driver
julia --project=. --startup-file=no --compiled-modules=existing --compile=yes --threads=1 -O2 benchmark/raw_pointcloud.jl results.toml
```

Use the SDK's `.dylib` on macOS or `.dll` on Windows. The driver also searches
common library locations when `GMSH_LIBRARY` is unset. `--check-driver` checks
case generation, the C ABI, and exact-oracle examples before loading Tessella.

Both adapters receive identical `D×N` Float64 coordinates and return owned
`3×N` Float64 coordinates with Int32 simplex connectivity. Timings include
coordinate extraction, output conversion, and release of Gmsh-owned buffers.
Tessella's exported mesh construction and validation costs remain included.
The [documented raw Gmsh APIs](https://gmsh.info/doc/texinfo/gmsh.html#Namespace-gmsh_002falgorithm)
accept interleaved coordinates and one-based constraint and element indices.

Julia, BLAS, and Gmsh use one thread. Each case primes both adapters, then runs
nine paired repetitions in a fixed randomized order. TOML results retain every
timing, quartiles, execution order, versions, build details, source hashes, and
host notes. Gmsh's default random factors are retained; additional small rows
disable both factors explicitly. The seed option is equal, but the two programs
do not share a random-number generator. Tessella 3D uses `perturb=false` to
retain supplied coordinates.

The ordinary families are uniform box clouds and random interior points inside
a fixed enclosing simplex, at 16, 64, 256, and 1024 points with three input seeds.
Additional 16-point clouds are thin, scaled, or translated. Simplex-enclosed
cases also run with the same explicit hull edges or faces supplied to both
adapters. Gmsh's constrained 3D endpoint requires HXT (`Mesh.Algorithm3D=10`);
that choice is recorded separately. Unconstrained 3D cases run both the default
Delaunay algorithm and HXT, retaining every outcome. Tessella's 3D adapter checks
face coverage before requesting recovery of an uncovered face.
The official [raw tetrahedralization example](https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_15_2/examples/api/raw_tetrahedralization.py)
selects HXT; use those rows as the primary 3D comparison.

Validation runs outside the timed region on every distinct output, including
the priming call. Garbage from preceding validation is collected before each
case's meshing measurements. All sizes use a common exact binary scale to check integer
determinants, opposite sides of shared faces, connected cells, a closed convex
boundary, area or volume conservation, and preservation and use of input points.
The 16-point cases also check every empty circumsphere using Rational{BigInt}
arithmetic. Larger rows are marked
`complete_mesh_pass`; they do not establish the empty-circumsphere criterion.
Outputs with added Steiner points are `oracle_unsupported`, because this driver
does not certify containment of added points in the original input hull.
Completion and empty-circumsphere results are separate: a complete triangulation
that fails the latter is `complete_mesh_non_delaunay`. Errors and
incomplete or invalid geometry remain separate outcomes. A faster incomplete
mesh is not a speed advantage for producing a complete mesh.

Measure fresh-process latency separately for the simplex-enclosed task with
explicit hull constraints:

```sh
julia --project=. --startup-file=no -O2 benchmark/raw_pointcloud.jl --latency Gmsh 3 64 gmsh-cold.toml
julia --project=. --startup-file=no -O2 benchmark/raw_pointcloud.jl --latency Tessella 3 64 tessella-cold.toml
```

The launcher measures a fresh Julia process through a flushed mesh-returned
marker, before independent validation. The TOML also separates package loading,
first meshing call, and driver time. Record whether the package cache was already
prepared: a first installation's precompilation differs from startup with an
existing cache. Other work
on a shared host can distort measurements; repeat on an idle host before
interpreting small timing differences or publishing speed ratios.

Fresh-process captures on the same host used 64 points in 3D, identical hull
constraints, and HXT for Gmsh. Times below are seconds; each row is one capture.
They include driver setup and compilation through the returned mesh.

| Implementation | Julia module cache | Through mesh | Package load | First mesh call |
| --- | --- | ---: | ---: | ---: |
| [Gmsh](results/apple_m3_julia_1_13_4_gmsh_cache_enabled.toml) | Enabled | 6.69 | <0.001 | 0.006 |
| [Tessella, first capture](results/apple_m3_julia_1_13_4_tessella_cache_enabled_first.toml) | Enabled | 104.15 | 89.39 | 0.463 |
| [Tessella, repeated import](results/apple_m3_julia_1_13_4_tessella_cache_enabled_repeat.toml) | Enabled | 13.57 | 0.246 | 0.488 |
| [Gmsh](results/apple_m3_julia_1_13_4_gmsh_cache_disabled.toml) | Disabled | 8.06 | <0.001 | 0.003 |
| [Tessella](results/apple_m3_julia_1_13_4_tessella_cache_disabled.toml) | Disabled | 111.92 | 92.46 | 0.524 |

Cache validity was unspecified for the first cache-enabled capture. The repeated
import followed a successful import of identical source. The disabled captures
used `--compiled-modules=no`; all captures used `--compile=yes -O2` and one
thread. Tessella's package loading and first-call compilation add startup cost
beyond the warm meshing measurements above.

## GEO startup

Run each case in a fresh Julia process from the repository root. Prepare the
dependencies and measure package precompilation separately:

```sh
JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. --startup-file=no -O2 --check-bounds=yes --threads=1 \
  -e 'using Pkg; Pkg.instantiate(); @time Pkg.precompile()'
```

Then measure a transfinite or unstructured unit square:

```sh
julia --project=. --startup-file=no --compiled-modules=existing --pkgimages=yes \
  --compile=yes -O2 --check-bounds=yes --threads=1 \
  benchmark/geo_startup.jl transfinite transfinite.toml
julia --project=. --startup-file=no --compiled-modules=existing --pkgimages=yes \
  --compile=yes -O2 --check-bounds=yes --threads=1 \
  benchmark/geo_startup.jl unstructured unstructured.toml
```

The TOML separates package loading, the first meshing call, and the median of
seven warm calls. It records cache availability, runtime options, source hash,
mesh size, and checksum. Each case checks positive triangle areas, unit-square
coverage, and identical coordinates and connectivity checksums across repeats.
Omit the output filename to print TOML to standard output.

The [example captures](results/apple_m3_julia_1_13_4_geo_startup.toml) report one
run of each case. Julia loaded a compatible cache with optimization level 2,
bounds checks, and one thread. The cache build was timed separately using
`Base.compilecache`. Re-run on an idle host for the source and Julia version
being used.

## Solid meshing

This benchmark meshes CAD solids end to end through the session API of both
programs: `generate(2)` for curves and surfaces, then `generate(3)` for the
tetrahedra, with one uniform mesh size and default options. It measures the
whole pipeline a user runs, including entity bookkeeping, not the kernels alone.

On the measured Windows host (Intel Core Ultra 7 265K, Julia 1.13.1, Gmsh
4.15.2, one thread each, Tessella 0.1.3), Tessella's volume pass took between
0.65 and 0.89 times Gmsh's on the six meshes above 100,000 tetrahedra, between
0.91 and 1.30 times on the six between 20,000 and 100,000, and between 0.95
and 1.13 times on the five below 20,000. Its surface pass took between 0.37
and 0.65 times Gmsh's in every case. Tessella is the faster program on large
meshes and on surfaces; on small and medium volumes the two are within a
third of each other, and Gmsh is ahead on seven of those eleven.

| Solid | Size | Tetrahedra (T / G) | Min. gamma (T / G) | Surface s (T / G) | Volume s (T / G) | Volume T/G |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| box | 0.1 | 4,977 / 4,979 | 0.30 / 0.30 | 0.007 / 0.014 | 0.030 / 0.032 | 0.95 |
| box | 0.05 | 37,210 / 37,262 | 0.27 / 0.30 | 0.028 / 0.053 | 0.281 / 0.216 | 1.30 |
| box | 0.03 | 178,858 / 178,869 | 0.26 / 0.28 | 0.099 / 0.159 | 1.190 / 1.337 | 0.89 |
| box | 0.02 | 560,539 / 560,591 | 0.26 / 0.26 | 0.222 / 0.342 | 3.908 / 6.047 | 0.65 |
| sphere | 0.2 | 2,702 / 2,583 | 0.30 / 0.27 | 0.006 / 0.014 | 0.016 / 0.016 | 1.04 |
| sphere | 0.1 | 20,663 / 20,564 | 0.30 / 0.29 | 0.021 / 0.053 | 0.105 / 0.112 | 0.93 |
| sphere | 0.06 | 88,575 / 91,035 | 0.25 / 0.30 | 0.059 / 0.144 | 0.516 / 0.569 | 0.91 |
| sphere | 0.04 | 299,836 / 294,952 | 0.21 / 0.19 | 0.119 / 0.324 | 1.804 / 2.607 | 0.69 |
| cylinder | 0.1 | 7,771 / 7,933 | 0.30 / 0.30 | 0.012 / 0.025 | 0.058 / 0.052 | 1.12 |
| cylinder | 0.05 | 62,101 / 59,581 | 0.30 / 0.21 | 0.042 / 0.100 | 0.446 / 0.414 | 1.08 |
| cylinder | 0.03 | 265,908 / 265,942 | 0.26 / 0.28 | 0.123 / 0.279 | 1.815 / 2.655 | 0.68 |
| torus | 0.15 | 14,484 / 14,468 | 0.30 / 0.26 | 0.027 / 0.060 | 0.130 / 0.134 | 0.97 |
| torus | 0.1 | 46,956 / 47,127 | 0.22 / 0.15 | 0.066 / 0.138 | 0.475 / 0.456 | 1.04 |
| torus | 0.07 | 133,270 / 133,188 | 0.21 / 0.12 | 0.144 / 0.285 | 1.246 / 1.482 | 0.84 |
| drilled block | 0.15 | 10,366 / 10,345 | 0.30 / 0.30 | 0.015 / 0.034 | 0.074 / 0.065 | 1.13 |
| drilled block | 0.08 | 59,612 / 59,861 | 0.21 / 0.30 | 0.051 / 0.109 | 0.467 / 0.400 | 1.17 |
| drilled block | 0.05 | 236,664 / 236,283 | 0.24 / 0.24 | 0.144 / 0.287 | 1.634 / 2.225 | 0.73 |

`T` is Tessella and `G` is Gmsh. Times are seconds; each driver meshes a case
twice in one process and reports the second run, which excludes compilation,
and the table gives the faster of two such driver runs. The Tessella driver
collects the garbage of earlier cases before it starts the clock. The
[retained results](results/) hold both runs of each program. Minimum gamma is
Gmsh's `gamma` quality of the worst tetrahedron: Tessella's is lower on five
rows, by up to 0.09 (the drilled block at 0.08), and higher on seven. The two
programs do not produce the same meshes: tetrahedron counts differ by under
5%, and the errors of the enclosed volumes against the exact solids agree to
within 1%.

The drilled block is a 2×2×2 box with a cylinder of radius 0.5 cut through it.

### Models with many faces

An n×n grid of unit plane faces at size 0.2 measures the cost that grows with
the number of CAD faces. Both programs produce identical node and triangle
counts. Times are the median of four driver runs:

| Faces | Triangles | Curves s (T / G) | Surfaces s (T / G) |
| ---: | ---: | ---: | ---: |
| 25 | 1,670 | 0.001 / 0.002 | 0.008 / 0.014 |
| 100 | 6,658 | 0.004 / 0.005 | 0.031 / 0.054 |
| 400 | 26,544 | 0.016 / 0.019 | 0.156 / 0.217 |

Tessella 0.1.2 took 9.8 s for the surfaces of the 400-face grid: it rebuilt
its mesh bookkeeping after every face. That is gone. What remains of the
bookkeeping is a record per mesh cell at several stages of a generation: in
the examples below, `tutorials/t10.geo` meshes its 46,000 triangles in 0.17 s
where Gmsh takes 0.40 s, and the `generate(3)` that follows, which has no
volume to mesh, takes 0.24 s and allocates 174 MiB where Gmsh takes 0.015 s.

### Gmsh's tutorials and examples

The third comparison runs Gmsh's own `.geo` files: every file of `tutorials`
and `examples/simple_geo` in the Gmsh 4.15.2 source tree, opened and meshed in
one, two and three dimensions by both programs. Four files that need a
post-processing view, a surface mesh or a background mesh file are left out.

Of the 35 files, Gmsh completes 33 (the two homology files stop for want of a
mesh, in both programs). Tessella completes 27 and stops on six more: a
compound surface (`t12`), two OpenCASCADE Boolean models (`t18`, `t19`), a
sliced STEP solid in two dimensions (`t20`), the volume of `cone.geo` and the
discrete `sphere-discrete.geo`. Its node counts equal Gmsh's on 28 files in
one dimension; node and cell counts are equal on 18 files in two dimensions
and on 15 in three. On the files that complete, a count that differs does so
by under 3%, except on `t5` (15% more triangles on its ruled sphere patches
and 28% more tetrahedra), `t15` (8% fewer tetrahedra) and `t21`, a
partitioned mesh (3,756 cells against 1,880).

End to end (open and three `generate` calls), on the 27 files that both
programs mesh, Tessella took a median of 1.5 times Gmsh's time (4.4 s against
2.9 s in total) and was the faster program on three (`splines.geo`, `t2`,
`t10`). The largest ratios are on the smallest structured files, where a few
milliseconds of bookkeeping are the whole run: `transfinite.geo` takes 0.052 s
against 0.005 s. These are small models, of at most 70,000 cells.

| File | 1-D nodes (T / G) | 2-D cells (T / G) | 3-D cells (T / G) | Seconds (T / G) | Tessella stops at |
| --- | ---: | ---: | ---: | ---: | --- |
| t1.geo | 80 / 80 | 728 / 728 | 0 / 0 | 0.01 / 0.01 |  |
| t10.geo | 575 / 575 | 45,762 / 45,762 | 0 / 0 | 0.43 / 0.46 |  |
| t11.geo | 66 / 66 | 3,484 / 3,485 | 0 / 0 | 0.16 / 0.16 |  |
| t12.geo | 67 / 111 | – / 308 | – / 0 | 0.00 / 0.02 | generate(2) |
| t14.geo | – / – | – / – | – / – | 0.01 / 0.00 | as Gmsh: generate(1) |
| t15.geo | 65 / 65 | 302 / 302 | 431 / 466 | 0.01 / 0.01 |  |
| t16.geo | 244 / 244 | 3,048 / 3,064 | 15,307 / 15,152 | 0.37 / 0.15 |  |
| t18.geo | – / 433 | – / 5,082 | – / 15,086 | 0.00 / 0.45 | open |
| t19.geo | 983 / 719 | – / 15,838 | – / 57,656 | 3.97 / 1.11 | generate(2) |
| t2.geo | 278 / 278 | 2,724 / 2,720 | 2,348 / 2,369 | 0.04 / 0.04 |  |
| t20.geo | 474 / 474 | – / 3,042 | – / 5,648 | 0.83 / 0.32 | generate(2) |
| t21.geo | 1,193 / 332 | 3,756 / 1,880 | 0 / 0 | 0.07 / 0.06 |  |
| t3.geo | 434 / 434 | 5,872 / 5,872 | 34,216 / 34,216 | 0.30 / 0.09 |  |
| t4.geo | 164 / 164 | 1,449 / 1,449 | 0 / 0 | 0.02 / 0.02 |  |
| t5.geo | 364 / 364 | 3,604 / 3,128 | 18,022 / 14,082 | 0.23 / 0.13 |  |
| t6.geo | 132 / 132 | 713 / 713 | 0 / 0 | 0.01 / 0.00 |  |
| t7.geo | 259 / 259 | 3,315 / 3,315 | 0 / 0 | 0.07 / 0.05 |  |
| antenna.geo | 5,110 / 5,110 | 21,518 / 21,518 | 0 / 0 | 0.81 / 0.36 |  |
| cone.geo | 124 / 124 | 1,160 / 1,248 | – / 2,460 | 0.02 / 0.03 | generate(3) |
| cube.geo | 44 / 44 | 260 / 260 | 375 / 373 | 0.01 / 0.01 |  |
| filter.geo | 570 / 570 | 4,588 / 4,588 | 11,592 / 11,369 | 0.19 / 0.13 |  |
| hex.geo | 132 / 132 | 288 / 288 | 192 / 192 | 0.01 / 0.00 |  |
| homology.geo | – / – | – / – | – / – | 0.01 / 0.00 | as Gmsh: generate(1) |
| indheat.geo | 1,807 / 1,807 | 7,362 / 7,412 | 37,893 / 37,741 | 0.60 / 0.42 |  |
| machine.geo | 1,296 / 1,296 | 14,823 / 14,863 | 0 / 0 | 0.20 / 0.15 |  |
| piece-extr-rec.geo | 1,215 / 1,215 | 6,146 / 6,146 | 12,230 / 12,230 | 0.13 / 0.10 |  |
| piece-extr.geo | 1,179 / 1,179 | 11,560 / 11,560 | 67,800 / 67,800 | 0.28 / 0.19 |  |
| piece.geo | 817 / 817 | 7,160 / 7,166 | 14,142 / 14,131 | 0.19 / 0.19 |  |
| pripyrtet.geo | 35 / 35 | 150 / 150 | 162 / 162 | 0.01 / 0.00 |  |
| sphere-discrete.geo | – / 326 | – / 648 | – / 1,639 | 0.02 / 0.01 | generate(1) |
| sphere.geo | 31 / 31 | 110 / 112 | 110 / 112 | 0.01 / 0.00 |  |
| splines.geo | 782 / 782 | 12,266 / 12,266 | 0 / 0 | 0.08 / 0.13 |  |
| square_regular.geo | 72 / 72 | 162 / 162 | 0 / 0 | 0.01 / 0.00 |  |
| tower.geo | 2,826 / 2,826 | 7,200 / 7,200 | 0 / 0 | 0.14 / 0.04 |  |
| transfinite.geo | 273 / 273 | 1,944 / 1,944 | 4,374 / 4,374 | 0.05 / 0.00 |  |

Times are seconds for the whole file, the faster of two driver runs; a dash
marks a step that was not reached.

### Running it

Install the Gmsh Python package (`pip install gmsh`); it is not a Tessella
dependency. From the repository root:

```sh
julia --project=. --startup-file=no -O2 --threads=1 benchmark/solid_meshing.jl tessella.toml
python benchmark/solid_meshing_gmsh.py gmsh.toml
julia --project=. --startup-file=no -O2 --threads=1 benchmark/solid_meshing.jl --faces tessella_faces.toml
python benchmark/solid_meshing_gmsh.py --faces gmsh_faces.toml
julia --project=. --startup-file=no -O2 --threads=1 benchmark/gmsh_examples.jl /path/to/gmsh-4.15.2-source tessella_examples.toml
python benchmark/gmsh_examples_gmsh.py /path/to/gmsh-4.15.2-source gmsh_examples.toml
```

Run the programs one after the other on an idle host. Each TOML lists, per
case, the node, triangle and tetrahedron counts, the enclosed volume, the
minimum and mean gamma, and the two times; the examples files list, per
`.geo` file, the time to open it, the time and the node and cell counts of
each `generate` call, and the step and message of an error. Omit the file
name to print to standard output.
