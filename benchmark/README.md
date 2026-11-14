# Raw point-cloud benchmark

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
