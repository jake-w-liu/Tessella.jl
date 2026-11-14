# Tessella.jl

Tessella generates meshes in Julia, with exact geometric predicates, CAD
modeling, size fields, and boundary recovery. Production meshing requires no
external mesher or CAD executable. Julia 1.12 and later are accepted; CI tests
Julia 1.12 and 1.13 on Linux, Windows, and macOS.

## Install

While the first registration is under review, install from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/jake-w-liu/Tessella.jl")
```

## Generate a volume mesh

```julia
using Tessella

const A = Tessella.API
A.initialize()
try
    A.model.occ.add_box(0, 0, 0, 1, 1, 1)
    A.model.occ.synchronize()
    A.mesh.generate(3)
    A.write("box.msh")
finally
    A.finalize()
end
```

The session API follows Gmsh conventions with Julia-style names. The
[feature reference](reference/index.md) describes supported geometry, fields,
constraints, formats, and post-processing. The API pages publish docstrings
for the top-level functions and all submodules, including internal helpers.
Support varies by geometry and operation; unsupported contexts raise explicit
errors. Complete Gmsh parity is not guaranteed.

## Command line

```sh
julia -e 'using Tessella; Tessella.CLI.main(ARGS)' file.geo -3 -o mesh.msh
```

See the [CLI controls](reference/cli_controls.md), or pass `-help`.

Developed with substantial AI assistance. Maintainer: Jake Liu.

## License

See the [package license](https://github.com/jake-w-liu/Tessella.jl/blob/main/LICENSE)
and [third-party source notices](https://github.com/jake-w-liu/Tessella.jl/blob/main/LICENSES/README.md)
for redistribution terms.
