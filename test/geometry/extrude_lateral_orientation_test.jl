module ExtrudeLateralOrientationTest

using Test, Tessella

# Upstream `orientMeshGFace` orients every swept lateral surface mesh along the
# lateral surface's normal, (signed generatrix) × (sweep direction). For a
# counter-clockwise source loop that is outward, including generatrices used
# reversed (curve 2 below) and the triangulated laterals of an unrecombined
# sweep (whose constrained-diagonal triangles are created inward and flipped).
function lateral_normals_outward(recombine::Bool)
    path = joinpath(mktempdir(), "sweep.geo")
    write(path, """
        Point(1)={0,0,0,0.5}; Point(2)={1,0,0,0.5}; Point(3)={1,1,0,0.5}; Point(4)={0,1,0,0.5};
        Line(1)={1,2}; Line(2)={3,2}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1,-2,3,4}; Plane Surface(1)={1};
        Transfinite Curve{1,2,3,4}=3; Transfinite Surface{1};
        $(recombine ? "Recombine Surface{1};" : "")
        Extrude {0,0,1} { Surface{1}; Layers{2}; $(recombine ? "Recombine;" : "") }
        """)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(3)
        tags, xyz, _ = A.mesh.get_nodes()
        at = Dict(t => (xyz[3i-2], xyz[3i-1], xyz[3i]) for (i, t) in enumerate(tags))
        outward = true
        for surface in (13, 17, 21, 25)
            types, _, nodes = A.mesh.get_elements(2, surface)
            for (type, conn) in zip(types, nodes)
                width = type == 3 ? 4 : 3
                for q in 1:width:length(conn)
                    a, b, c = at[conn[q]], at[conn[q+1]], at[conn[q+2]]
                    u = b .- a; v = c .- a
                    n = (u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])
                    centre = (a .+ b .+ c) ./ 3 .- (0.5, 0.5, 0.5)
                    sum(n .* centre) > 0 || (outward = false)
                end
            end
        end
        return outward
    finally
        A.finalize()
    end
end

@testset "swept lateral surfaces follow the lateral normal" begin
    @test lateral_normals_outward(true)
    @test lateral_normals_outward(false)
end

end # module
