# Pyramid options for Region batches

The ordinary hybrid batch producer forwards the generation's finite
`Mesh.OptimizePyramids` value to each Region component. Negative
values construct the actual normal-offset apices and certify their strict
containment before Delaunay recovery. Exactly `1.0` runs Winslow optimization
over actual Region-owned apex stars. Other nonnegative values run the existing
Pyramid relocation kernel.

All components continue the same global node and element UID factories.
The batch builds and binds every component's Pyramid products before deferred
Tet optimization. Each completed Region must preserve its literal lower
constructor, owned typed exterior, positive cells and certified volume.
Negative-factor apices retain their actual allocated UID and exact position
through that final optimization. Regions connected by shared Triangles also
support these settings. Their actual common Surface incidence partitions the
Tet cells into declared Regions, and each optimizer sees only its Region's
stars. For negative factors, every recorded global birth binds to the actual
Region constructor by Quad UID, ordered corner UIDs and literal Source XYZ;
the producer certifies each Region's containment before recovery. Shared Quad
recovery still requires a two-sided Pyramid constructor.

The opt-in regression
`test/interfaces/api_hybrid_batch_pyramids_functional_test.jl` uses two separated
unit cubes with `-0.1` and `1.0`, plus bounded atomic failures. It checks
analytic normal positions, independent Winslow one-apex positions, both unit
volumes, conforming exteriors, global UID uniqueness, Region ownership and
MSH4 transport. The preceding separated-component Source passed 74 assertions
in a warm Julia 1.13 process with
`--compile=min -O0` and bounds checks in 59.37 seconds. The default suite is
unchanged. Its later shared-Quad atomic test correction is not a new full
74-check qualification. The final shared-component Source is qualified by a
second opt-in regression,
`test/interfaces/api_hybrid_shared_tri_pyramids_functional_test.jl`, checks
adjacent unit cubes with a shared triangulated Surface. Its 61 assertions
passed with the same Julia flags in 85.56 seconds including 35.73 seconds of
cold loading. Both negative and Winslow cases produce positive unit volumes,
retain the literal shared incidence and preserve NodeUID/Region-owner
transport through MSH4. The existing 35-assertion default0 hybrid regression
also passed once against the final batch Source in 36.94 seconds warm.
