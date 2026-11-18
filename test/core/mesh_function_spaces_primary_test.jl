# The suite includes every test file into Main; this file keeps its aliases to itself.
module FunctionSpacesPrimaryTests
using Test
using Tessella
import Tessella.MeshFunctionSpaces as FS

# The first-order H1 and lowest-order H(curl) functions belong to the element
# family: a 3-node Line or a 6-node Triangle has those of its straight primary
# element, as in Gmsh 4.15.2. A higher-order type used to select the
# Tetrahedron's tables and write past the result sized for its own family.

@testset "First-order hierarchical bases on higher-order simplices" begin
    points = Dict(1 => [-0.3, 0.0, 0.0, 0.6, 0.0, 0.0],
                  2 => [0.2, 0.3, 0.0, 0.6, 0.1, 0.0],
                  4 => [0.1, 0.2, 0.3, 0.3, 0.1, 0.2])
    # getBasisFunctions(type, points, "H1Legendre1", [0]) in Gmsh 4.15.2.
    oracle = Dict(1 => [0.65, 0.35, 0.2, 0.8],
                  2 => [0.5, 0.2, 0.3, 0.3, 0.6, 0.1],
                  4 => [0.4, 0.1, 0.2, 0.3, 0.4, 0.3, 0.1, 0.2])
    orientations = Dict(1 => 2, 2 => 6, 4 => 24)
    for (higher, primary) in ((8, 1), (26, 1), (9, 2), (21, 2), (11, 4), (29, 4))
        components, values, count = FS.mesh_basis_functions(higher, points[primary], "H1Legendre1", [0])
        @test components == 1 && count == orientations[primary]
        @test values ≈ oracle[primary] atol = 1e-14
        for space in ("H1Legendre1", "GradH1Legendre1", "HcurlLegendre0", "CurlHcurlLegendre0")
            @test FS.mesh_basis_functions(higher, points[primary], space, [0]) ==
                  FS.mesh_basis_functions(primary, points[primary], space, [0])
            every = collect(0:orientations[primary]-1)
            @test FS.mesh_basis_functions(higher, points[primary], space, every) ==
                  FS.mesh_basis_functions(primary, points[primary], space, every)
        end
    end
end
end
