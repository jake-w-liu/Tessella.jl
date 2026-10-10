# Gmsh size-field differential

## Current continuation status (2026-10-10)

This is an oracle-only differential runner; production meshing stays independent
of Gmsh. Use supported Julia 1.12.x and 1.13.x for current verification, with the
same source and Gmsh 4.15.2 inputs. The historical 2026-08-21 measurements are in
[STATUS.md](STATUS.md), not a fresh qualification of Root V24/release V6.
Read the [current handoff](../../HANDOFF.md#continuation-handoff-2026-10-10)
before launching this package/aggregate reader: adopt the existing owner and
follow the fresh-memory/hash guards. A context skip remains unfinished coverage
under the full parity goal, not a completion waiver. Verified implementations
still require a normal push to `main`.

`differential.jl` requires the Gmsh 4.15.2 CLI and Julia API. It checks field
values directly with Gmsh's official `MeshSizeFieldView` plugin, then retains
the original line-mesh grading checks for `Distance→Threshold`, `Box`, `Ball`,
finite `Cylinder`, and `Frustum`.

Run the bounds-checked gate with:

```sh
julia --project=. --startup-file=no --history-file=no --check-bounds=yes -O2 --threads=1 --gcthreads=1 validation/size_fields/differential.jl
```

The runner searches `PATH`, executable-relative library paths, and historical
Homebrew locations. Bind both paths explicitly when those defaults do not match
the current machine; successful startup validates the actual CLI/API versions:

```sh
GMSH_EXECUTABLE=/path/to/gmsh \
GMSH_JULIA_API=/path/to/gmsh.jl \
julia --project=. --startup-file=no --history-file=no --check-bounds=yes -O2 --threads=1 --gcthreads=1 validation/size_fields/differential.jl
```

On Windows, the same explicit bindings can be supplied in PowerShell (replace
the placeholders with the pinned paths; run once with each supported Julia):

```powershell
$env:GMSH_EXECUTABLE = 'C:/path/to/gmsh.exe'
$env:GMSH_JULIA_API = 'C:/path/to/gmsh.jl'
& 'C:/path/to/julia.exe' --project=. --startup-file=no --history-file=no --check-bounds=yes -O2 --threads=1 --gcthreads=1 validation/size_fields/differential.jl
```

These are invocation examples, not commands executed by this documentation review.
Successful output contains all of the following:

- `GMSH_RUNTIME_OK` with CLI/API versions and resolved paths;
- one `DIRECT` line per pointwise comparison;
- five `MESH_APPROX` lines for mesh-observed grading;
- a counted `CONTEXT_DEPENDENT_SKIPS` block for behavior without an equivalent
  public oracle or shared model context;
- a final `SIZE_FIELD_DIFFERENTIAL_OK` summary.

Missing/wrong-version Gmsh, plugin failures, empty probes, and parity mismatches
exit nonzero. `CONTEXT_SKIP` means only that the named behavior is outside the
available differential oracle; it is never printed as a pass. See
[`STATUS.md`](STATUS.md) for the exact classification and limitations.

`validation/run_all.jl` launches this script as a required bounds-checked child
process, so its nonzero status propagates to the aggregate validation command.
