# IGES plane parameter charts

`import_iges` maps UV-only CurveOnSurface142 wires on an explicit
PlaneSurface190 Form1 through its declared chart: `C + u*x + v*y`.
Point116 gives the origin; Direction123 records give the normal and reference
direction. The importer normalizes them, certifies their approximate
perpendicularity and removes the small normal projection from the reference
axis. The resulting right-handed frame also certifies Face coplanarity.

Supported UV wires are bounded Line110 Form0 and ordered CompositeCurve102
Form0 records. Literal endpoints must join and close with exactly zero UV Z.
The wrapper must reference the Face's actual plane. Existing model-space
CurveOnSurface records retain priority and their vertices stay in model space.

A finite chain of proper rigid TransformationMatrix124 Form0 records may be
referenced by the plane's directory entry. Each matrix's directory field7
references its later parent transform. Composition uses `Rparent*Rchild` and
`Tparent + Rparent*Tchild`; the result maps the origin as `R*C + T` and the
axes as `R*x` and `R*y`.
Every matrix and composed rotation must be finite and orthonormal with
positive determinant. Cycles, missing/invalid parents, accumulated rigidity
loss, reflections, scaling and transforms on the Point/Direction/UV wire
records reject. Implicit Plane108 charts and nonlinear UV/NURBS
charts remain unsupported. This mesh importer does not create an analytic CAD
registry or general trimmed/periodic NURBS topology.

These field meanings follow sections 4.50 and 4.21 of the
[IGES 5.3 standard](https://paulbourke.net/dataformats/iges/IGES.pdf)
(printed pages 203 and 113–115), and the public
[PlaneSurface schema](https://occt3d.com/dev/doc/refman/html/class_i_g_e_s_solid___plane_surface.html).
Mapped rings retain the existing hole-domain, orientation, closed-shell and
positive-volume certificates.

The untransformed chart checks passed 49 assertions on Julia 1.12 in 22.43 s;
the transformed chart checks passed 38 in 23.11 s, actual exit zero, using
minimal compilation and bounds checks. They independently verify mapped
corners, signed shell volume, Face areas, void witnesses, positive Tet volume,
boundary incidence, assembly transport, model-space priority and failure
atomicity. Only five inexpensive Surface assertions enter the default suite;
full tests are in `test/geometry/brep_plane190_{uv,transform}_functional_test.jl`.
Opaque Gmsh 4.15.2 CAD imports confirmed Face geometry and areas; they generated
no meshes and establish no native meshing or element-order parity.

The positive transformed-hole fill uses an exact quarter-turn. Arbitrary
rotation is verified for Surface import. A diagnostic fill of an arbitrary
rotated nonstar-shaped hole was stopped after 415.44 s in unchanged generic
conforming-Delaunay recovery; its failed receipt and verified process identity
remain in the isolated transform packet. Arbitrary-rotation hole fill and its
performance remain unqualified.

The chain increment passed 26 assertions on Julia 1.12, actual exit zero in
32.48 s including compilation; body 2.7 s. Two cheap Surface assertions join
the default suite. Two/three-matrix fixtures match independent sequential
point mapping and a collapsed single-matrix control. A reproduced failure
where individually tolerated matrices accumulated nonrigid scaling is fixed
by certifying each composite as well.

The opaque Gmsh 4.15.2 CAD import applies only the child matrix to these chain
fixtures, while its collapsed single-matrix control agrees with the standard
composition. Tessella follows IGES section4.21's declared parent/child rule;
this chain support establishes standard transport rather than native chain
import parity. No volume meshing occurs in the chain tests.
