module GeoNegativePlaneLoopTest

using Test, Tessella

# `Plane Surface(1) = {-1}`: upstream `SetSurfaceGeneratrices` walks a
# negative exterior wire backwards with every curve negated, so the face runs
# the other way. Gmsh 4.15.2 on the counter-clockwise unit square below:
# normal (0,0,-1), oriented boundary {-1,-2,-3,-4}, and with a transfinite
# surface all 8 triangles facing -z (tutorial t6 uses the same form).
function face(extra::AbstractString, loop::AbstractString)
    path = joinpath(mktempdir(), "negative.geo")
    write(path, """
        Point(1)={0,0,0,0.5}; Point(2)={1,0,0,0.5}; Point(3)={1,1,0,0.5}; Point(4)={0,1,0,0.5};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};
        Plane Surface(1)={$loop};
        $extra
        """)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(2)
        tags, xyz, _ = A.mesh.get_nodes()
        pos = Dict(tags[i] => (xyz[3i-2], xyz[3i-1]) for i in eachindex(tags))
        _, _, nodes = A.mesh.get_elements(2, 1)
        up = down = 0
        for c in 1:3:length(nodes[1])
            a, b, d = pos[nodes[1][c]], pos[nodes[1][c+1]], pos[nodes[1][c+2]]
            nz = (b[1]-a[1])*(d[2]-a[2]) - (b[2]-a[2])*(d[1]-a[1])
            nz > 0 ? (up += 1) : (down += 1)
        end
        return (normal = A.model.get_normal(1, [0.5, 0.5]),
                boundary = Set(t for (_, t) in A.model.get_boundary([(2, 1)], false, true)),
                up = up, down = down)
    finally
        A.finalize()
    end
end

@testset "a negative exterior loop reverses the plane surface" begin
    positive = face("", "1")
    @test positive.normal[3] ≈ 1.0
    @test positive.boundary == Set([1, 2, 3, 4])
    negative = face("", "-1")
    @test negative.normal[3] ≈ -1.0                  # Gmsh 4.15.2: (0,0,-1)
    @test negative.boundary == Set([-1, -2, -3, -4]) # Gmsh 4.15.2
    transfinite = face("Transfinite Curve{:}=3; Transfinite Surface{1};", "-1")
    @test (transfinite.up, transfinite.down) == (0, 8) # Gmsh 4.15.2: 8 facing -z
end

end # module
