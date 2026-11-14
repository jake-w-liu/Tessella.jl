# Contributing

Keep production meshing independent of external meshers and CAD executables.
Preserve geometry, topology, entity ownership and file-format behavior when
changing the implementation. Add a small reproducer for a bug fix and test
the affected API, including failure paths that must leave the model unchanged.
Retain third-party copyright and license notices when adapting source code.

## Run tests

Run package tests with bounds checking from the repository root:

```sh
julia --project=. --startup-file=no --check-bounds=yes -e 'using Pkg; Pkg.test()'
```

The [test guide](https://github.com/jake-w-liu/Tessella.jl/blob/main/test/README.md)
describes the default suite, optimized compiler checks and extended tests.
Use the test relevant to a change before running the full suite. Changes to
numerical code should also run with normal optimization, since compiler mode
can affect floating-point behavior and allocation measurements.

## Build the documentation

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html`. The build includes source docstrings and checks
missing docstrings and unresolved cross-references. CI publishes documentation
from `main` to the site's `dev` version. Keep usage examples and supported-input
limits in the feature reference alongside the affected API docstrings.
