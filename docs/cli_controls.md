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
