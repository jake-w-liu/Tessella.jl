# Gmsh size-field differential

`differential.jl` requires the Gmsh 4.15.2 CLI and Julia API. It checks field
values directly with Gmsh's official `MeshSizeFieldView` plugin, then retains
the original line-mesh grading checks for `Distance→Threshold`, `Box`, `Ball`,
finite `Cylinder`, and `Frustum`.

Run the bounds-checked gate with:

```sh
julia --project=. --startup-file=no --check-bounds=yes validation/size_fields/differential.jl
```

The Homebrew 4.15.2 installation is discovered automatically. Override either
path when needed:

```sh
GMSH_EXECUTABLE=/path/to/gmsh \
GMSH_JULIA_API=/path/to/gmsh.jl \
julia --project=. --startup-file=no --check-bounds=yes validation/size_fields/differential.jl
```

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
