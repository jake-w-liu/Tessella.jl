# Tessella.jl

Tessella is a Julia-native mesh generator with exact geometric predicates,
conforming boundary recovery, CAD modeling, and explicit failure diagnostics.
Production meshing runs in Julia without an external mesher or CAD executable.
Julia 1.12.x and 1.13.x are supported.

Public Gmsh 4.15.2 compatibility remains an active implementation goal. The
[development checklist](.github/development/CHECKLIST.md) records the full scope;
a passing functional suite does not establish complete parity.

## Install

Until registration is complete, install from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/jake-w-liu/Tessella.jl")
```

## Generate a mesh

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

The model, mesh, option, field, view, and plugin APIs follow Gmsh conventions
with Julia-style names. Native capabilities include 2-D/3-D meshing,
structured grids, recombination, high-order elements, size fields, periodic
and embedded constraints, partitions, and STEP/IGES geometry. Unsupported
contexts return explicit errors. See the [feature reference](.github/reference/README.md)
for API examples and supported domains.

## Command line

```sh
julia -e 'using Tessella; Tessella.CLI.main(ARGS)' file.geo -3 -o mesh.msh
```

Native `.geo` execution supports modeling, sizing, physical groups, meshing,
and post-processing. The CLI accepts `-1`, `-2`, `-3`, `-format msh2|msh4`,
`-bin`, `-order`, and other [documented controls](.github/reference/cli_controls.md).
Use `-help` for the current option list or `file.msh -check` to validate a mesh.

## Test and contribute

```sh
julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

The default tests use small real fixtures and analytic/topological checks.
They select minimal compilation automatically and need no external mesher.
Optimized compiler checks and extended campaigns are described in
[test/README.md](test/README.md). Read the
[development instructions](.github/development/README.md) before changing code.

Developed with substantial AI assistance in implementation, tests, and
documentation. Maintainer: Jake Liu.

## License

Tessella is available under the [MIT license](LICENSE). Third-party code retains
its respective [license and provenance notices](LICENSES/).
