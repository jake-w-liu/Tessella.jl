using Test
using Tessella

# Upstream `meshGFaceTransfinite` applies `Mesh.Smoothing - 1` (or the
# surface's `Smoother`/`setSmoothing` count) sweeps of its elliptic smoother
# to four-corner transfinite faces, in the face's (u,v) chart: rigid for a
# Plane Surface, the ruled chart for the extruded lateral `Surface`s.
# Expected interior nodes are Gmsh 4.15.2's for this script.
const _SMOOTHED_INTERIOR = [
    (0.23013991143936594, 0.7714440806691288), (0.23136096203492068, 0.23136402865901526),
    (0.24253300984315274, 0.5023326755602727), (0.5010830588092773, 0.7604748101941715),
    (0.5022491434581428, 0.5022501297246098), (0.5023773006665808, 0.24257869983491806),
    (0.7605013799020273, 0.5011093232225087), (0.7710917161950613, 0.7710892017838424),
    (0.7715452298638017, 0.23024195944912806)]

function _smoothed_face_interior(A, dim, tag, plane)
    tags, xyz, _ = A.mesh.get_nodes(dim, tag, false, false)
    pts = [Tuple(xyz[3i-2:3i]) for i in eachindex(tags)]
    return sort([plane(p) for p in pts])
end

@testset "transfinite elliptic smoothing matches Gmsh" begin
    A = Tessella.API
    path = joinpath(mktempdir(), "smooth.geo")
    write(path, """
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
        Extrude {0,0,1} { Surface{1}; }
        Transfinite Curve{:} = 5 Using Bump 0.3;
        Transfinite Surface{:};
        Transfinite Volume{1};
        Mesh.Smoothing=5;
        """)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(3)
        expected = sort(_SMOOTHED_INTERIOR)
        plane = _smoothed_face_interior(A, 2, 1, p -> (p[1], p[2]))
        ruled = _smoothed_face_interior(A, 2, 13, p -> (p[1], p[3]))
        @test length(plane) == 9 && length(ruled) == 9
        @test maximum(maximum(abs.(a .- b)) for (a, b) in zip(plane, expected)) < 1e-12
        @test maximum(maximum(abs.(a .- b)) for (a, b) in zip(ruled, expected)) < 1e-10
        # The default `Mesh.Smoothing = 1` leaves the Bump grid unsmoothed.
        A.option("Mesh.Smoothing", 1)
        A.mesh.clear()
        A.mesh.generate(2)
        raw = _smoothed_face_interior(A, 2, 1, p -> (p[1], p[2]))
        @test maximum(maximum(abs.(a .- b)) for (a, b) in zip(raw, expected)) > 1e-3
    finally
        A.finalize()
    end
end
