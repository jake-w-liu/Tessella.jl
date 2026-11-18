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
4.15.2, one thread each), Tessella's volume pass took between 0.92 and 1.30
times Gmsh's on the six meshes above 100,000 tetrahedra, and between 1.0 and
1.1 times on the largest mesh of each curved solid. It took 1.3 to 1.7 times
as long between 20,000 and 100,000 tetrahedra and 1.15 to 2.2 times below
20,000. The surfaces of the sphere, cylinder and torus meshed in 0.5 to 1.0
times Gmsh's time; the planar faces of the box took 1.4 to 2.1 times as long
and the drilled block 1.1 to 2.0 times. Tessella is faster than Gmsh in one
volume case, the 562,000-tetrahedron box.

| Solid | Size | Tetrahedra (T / G) | Min. gamma (T / G) | Surface s (T / G) | Volume s (T / G) | Volume T/G |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| box | 0.1 | 4,966 / 4,979 | 0.30 / 0.30 | 0.024 / 0.014 | 0.064 / 0.029 | 2.20 |
| box | 0.05 | 37,357 / 37,262 | 0.30 / 0.30 | 0.097 / 0.054 | 0.356 / 0.214 | 1.67 |
| box | 0.03 | 178,145 / 178,869 | 0.28 / 0.28 | 0.328 / 0.153 | 1.944 / 1.494 | 1.30 |
| box | 0.02 | 562,185 / 560,591 | 0.20 / 0.26 | 0.637 / 0.465 | 6.658 / 7.245 | 0.92 |
| sphere | 0.2 | 2,736 / 2,583 | 0.30 / 0.27 | 0.009 / 0.019 | 0.031 / 0.021 | 1.52 |
| sphere | 0.1 | 20,225 / 20,564 | 0.26 / 0.29 | 0.041 / 0.072 | 0.186 / 0.146 | 1.28 |
| sphere | 0.06 | 90,510 / 91,035 | 0.30 / 0.30 | 0.113 / 0.195 | 1.220 / 0.707 | 1.73 |
| sphere | 0.04 | 297,451 / 294,952 | 0.27 / 0.19 | 0.269 / 0.434 | 3.322 / 3.192 | 1.04 |
| cylinder | 0.1 | 8,136 / 7,933 | 0.30 / 0.30 | 0.020 / 0.035 | 0.075 / 0.066 | 1.15 |
| cylinder | 0.05 | 60,828 / 59,581 | 0.30 / 0.21 | 0.104 / 0.135 | 0.723 / 0.504 | 1.43 |
| cylinder | 0.03 | 265,991 / 265,942 | 0.29 / 0.28 | 0.341 / 0.357 | 3.262 / 3.226 | 1.01 |
| torus | 0.15 | 14,484 / 14,468 | 0.30 / 0.26 | 0.043 / 0.070 | 0.204 / 0.158 | 1.29 |
| torus | 0.1 | 46,956 / 47,127 | 0.22 / 0.15 | 0.143 / 0.157 | 0.738 / 0.574 | 1.28 |
| torus | 0.07 | 133,270 / 133,188 | 0.21 / 0.12 | 0.312 / 0.316 | 1.921 / 1.781 | 1.08 |
| drilled block | 0.15 | 10,275 / 10,345 | 0.30 / 0.30 | 0.043 / 0.039 | 0.122 / 0.081 | 1.51 |
| drilled block | 0.08 | 59,725 / 59,861 | 0.29 / 0.30 | 0.282 / 0.139 | 0.662 / 0.485 | 1.37 |
| drilled block | 0.05 | 236,509 / 236,283 | 0.26 / 0.24 | 0.421 / 0.366 | 2.727 / 2.668 | 1.02 |

`T` is Tessella and `G` is Gmsh. Times are seconds; each driver meshes a case
twice in one process and reports the second run, which excludes compilation,
and the table gives the faster of two such driver runs. The
[retained results](results/) hold both runs of each program. Minimum gamma is
Gmsh's `gamma` quality of the worst tetrahedron. The two programs do not produce
the same meshes: tetrahedron counts differ by under 3%, 6% on the coarsest
sphere, and the errors of the enclosed volumes against the exact solids agree
to within 1%.

The drilled block is a 2×2×2 box with a cylinder of radius 0.5 cut through it.

### Models with many faces

The cost that the table above does not show grows with the number of CAD
faces. Tessella rebuilds its mesh bookkeeping after every face it meshes, so
the surface pass grows with the square of the face count. For an n×n grid of
unit plane faces at size 0.2, where both programs produce identical node and
triangle counts:

| Faces | Triangles | Curves s (T / G) | Surfaces s (T / G) |
| ---: | ---: | ---: | ---: |
| 25 | 1,670 | 0.013 / 0.001 | 0.091 / 0.014 |
| 100 | 6,658 | 0.048 / 0.005 | 0.659 / 0.052 |
| 400 | 26,544 | 0.206 / 0.018 | 9.760 / 0.206 |

The triangulation itself takes about a quarter of a second of the 400-face
run; the rest is the per-face rebuild. Until that is made incremental, expect
models with hundreds of faces to mesh far more slowly than in Gmsh.

### Running it

Install the Gmsh Python package (`pip install gmsh`); it is not a Tessella
dependency. From the repository root:

```sh
julia --project=. --startup-file=no -O2 --threads=1 benchmark/solid_meshing.jl tessella.toml
python benchmark/solid_meshing_gmsh.py gmsh.toml
julia --project=. --startup-file=no -O2 --threads=1 benchmark/solid_meshing.jl --faces tessella_faces.toml
python benchmark/solid_meshing_gmsh.py --faces gmsh_faces.toml
```

Run the programs one after the other on an idle host. Each TOML lists, per
case, the node, triangle and tetrahedron counts, the enclosed volume, the
minimum and mean gamma, and the two times. Omit the file name to print to
standard output.
