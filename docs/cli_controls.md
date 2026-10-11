# CLI controls

`tessella file.geo -1|-2|-3` generates a mesh. Use `-o output.msh`,
`-format msh2|msh4` and `-bin` to choose its destination and MSH encoding.
Input/output aliases and duplicate/conflicting flags are rejected before
execution. Without a dimension flag, the CLI executes the script and returns
its input path; output-format flags require a dimension.

| Flag | Meaning | Accepted values |
| --- | --- | --- |
| `-clscale` | Mesh size factor after bounds | Finite, positive |
| `-clmin` | Minimum mesh size | Finite, nonnegative |
| `-clmax` | Maximum mesh size | Finite, positive, at least the minimum |
| `-order` | Generated element order | Integer 1–10; generated family support applies |

Command-line values retain priority over script assignments, including sizing
aliases and in-script mesh statements. This follows the pinned public Gmsh
4.15.2 initialization/open observation. Curve generation supports orders 1–10;
Surface/Volume orders retain the [implemented scope](../README.md).
Unsupported generated families/orders fail explicitly.

`-help`/`--help` prints and returns usage; `-version`/`--version` prints and
returns the package version. Use an informational flag alone, without a file.
Numeric values, duplicates and incompatible bounds are checked before parsing
the input or replacing output.

Flag meanings follow the [Gmsh command-line manual](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-command_002dline-interface).
The focused regression uses a single tiny public Line, analytic spacing/order,
option precedence and failure atomicity; it runs without an external mesher.

`file.msh -check` reads MSH2.2/4.1 ASCII or binary and returns the input path
without generating or writing a mesh. It checks finite coordinates, connectivity,
positive unique stored UIDs, aligned ownership/periodic/partition/data metadata,
and duplicate cells. Missing references and repeated cells raise errors.
This mode rejects generation, sizing and output flags. It checks structure;
curved Jacobian validity and native diagnostic chronology remain separate work.

The opt-in [mesh-check fixture](../test/interfaces/cli_check_functional_test.jl)
passed all20 fixed assertions onJulia1.13.1, actual exit0/35.63 s including
loading. It covers all four encodings, unchanged input bytes/directory contents,
duplicate and missing-node failures, and conflicting flags. Two no-generation
assertions enter the concise default suite. Tiny public native4.15.2 CLI controls
confirm duplicate-cell failure; a zero-volume warning is outside this structural
certificate. The observation is retained under
`C:/tmp/tessella_cli_check_20261011/native_cli_check_observation.json`.
