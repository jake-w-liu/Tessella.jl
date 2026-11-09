# Tessella continuation record — 2026-11-10

This record supersedes the 2026-11-09 stop handoff. Work resumed from it and
completed the concrete pending correctness tracks except registration and the
remaining long-tail checklist families listed at the end.

## Published state

Code HEAD and verified remote `main`: `461f84dec75b77c2f942864db2849885702814a2`.
Later commits in this working tree are documentation/checklist only under
`.github/development/`. The pre-resume safety snapshot remains at
`C:/Users/User/Tessella-handoff-20261109-1716` with its hash manifest and
patch; it was not overwritten.

## Completed in this continuation

- Native ONELAB process server: real named clients, nested CONNECT, complete
  vector/metadata stores, detached session staging, owned process-tree cleanup,
  `DefineString` Gmsh value/`Name` semantics (including empty server values),
  finite run/connect deadlines enforced by an independent watchdog, and
  kill-before-socket-close cleanup for non-reading clients.
- Field inverse normalization: the FMA equality check now evaluates the true
  `(point-base)-query*scale` residual with exact tails and conservatively uses
  the existing precise solver for overflow/underflow. Thin-map values are
  recovered.
- Boolean seams and fuse images: crossed-cylinder seam anchors, periodic N-path
  snapping, original-source image tags, retained operands and internal-wall
  tracking. The 44-case public fuse fixture and crossed-cylinder seam/mass/
  centroid/inertia controls are default-wired.
- Natural BREP: finite natural sphere/torus/B-spline/plane faces and closed
  natural solids consume the actual Face mesh without invented Point/Curve
  entities. Natural solid centroid and centroidal inertia are supported.
- CAD measures: `get_matrix_of_inertia` now returns Gmsh-compatible centroidal
  tensors for curves, plane/conic faces, primitives, plane/cylinder/natural
  solids; cylinder-bounded solid centers of mass are implemented; conic plane
  boundary moments support disk inertia. `add_wedge` now matches the public
  OpenCASCADE/Gmsh xy-trapezoid profile extruded along z.
- ONELAB Unix portability: `Cmd(Vector{String}; dir, detach)` is wrapped in the
  keyword-bearing constructor, spawn failures become `ArgumentError`, GEO
  command embedding is POSIX-safe, and short nested clients get a bounded
  reader grace before a non-connect diagnosis.

## Qualification

Final ordinary `Pkg.test()` on clean `461f84d`:

- Julia 1.12.7: **4976/4976 passed**, 91m40.5s, actual exit 0.
  `C:/Users/User/Tessella-final-461f84d-112-20261110.log` (+ `.exit`).
- Julia 1.13.3: **4976/4976 passed**, 96m60.0s reported, actual exit 0.
  `C:/Users/User/Tessella-final-461f84d-1133-20261110.log` (+ `.exit`).

Linux WSL Ubuntu (Julia 1.12.7, minimal) reproduced the original Unix ONELAB
failure and passed after the fix: ONELAB 135/135, inertia 34/34, wedge 720/720,
natural BREP 182/182, crossed-cylinder seam 140/140, fuse images 134/134.
Those focused files cover every new default test on Linux. GitHub Actions run
`42685928378` passed all four Ubuntu/Windows × Julia 1.12/1.13 jobs on the
final `d75d8ef` tree (docs/checklist tree; production code is `461f84d`).

Other focused receipts are listed in
[development notes](README.md#local-verification-receipts). The exact residual
checks were: field FMA 44/44 on both runtimes, ONELAB 218/218 on both runtimes,
and the confirmed red witness now matches the exact rational thin-map value.

## Remaining work

- Registration is **not** complete. The release checklist still needs the
  maintainer's verified human review/understanding and a Registrator trigger;
  no GitHub token, `gh`, SSH key or Registrator action was available in this
  environment.
- The 27-domain checklist still has broad long-tail work. Confirmed remaining
  BREP families include nondegenerate pcurve-only edges, hyperbolas/general
  offsets and internal/external orientations; general quadric-bounded shell
  inertia still rejects rather than approximating. The wider residuals
  (full-quad/anisotropic contexts, broader NURBS Boolean/fillet, parallel
  eligibility, full GUI/matching identity and systematic conformance) remain
  as recorded in `CHECKLIST.md`.
- No final independent release-readiness audit or complete-optimization claim is
  made. The full-package results above qualify the tested source and scopes.

Registration policy: [General](https://github.com/JuliaRegistries/General#registering-a-package-in-general),
[LLM contributions](https://github.com/JuliaRegistries/General#is-there-any-policy-regarding-the-use-of-llms-in-registered-packages),
[Registrator](https://github.com/JuliaRegistries/Registrator.jl#via-the-github-app).
