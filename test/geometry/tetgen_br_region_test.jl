module TetgenBRRegionTest

using Test, Tessella

const API = Tessella.API

# Gmsh 4.15.2 runs part of TetGen on every Delaunay Region before refining it
# (meshGRegionBoundaryRecovery.cpp: reconstruct the tetrahedralization,
# recover the boundary, carve the exterior, recover Delaunayness, optimize).
# Mesh.Algorithm3D = 3 stops after that pass. The references are Gmsh's: node
# and tetrahedron counts, and a checksum of the tetrahedra in Gmsh's element
# order with their node tags in Gmsh's vertex order
# (sum over i of i*(a+7b+49c+343d) mod 1e9+7), so order and vertex order of
# the initial mesh, of the refinement (insertVerticesInRegion) and of the
# optimizer are all held to Gmsh's.
function mesh_checksum(build; initial=false)
    API.initialize()
    try
        API.option("General.Terminal", 0)
        build()
        if initial
            API.option("Mesh.Algorithm3D", 3)
            API.option("Mesh.Optimize", 0)
        end
        API.mesh.generate(3)
        tags = API.mesh.get_nodes()[1]
        types, _, nodes = API.mesh.get_elements(3)
        tets = Int[n for (t, ns) in zip(types, nodes) if t == 4 for n in ns]
        checksum = 0
        for i in 1:length(tets) ÷ 4
            a, b, c, d = tets[4i-3:4i]
            checksum = mod(checksum + i * (a + 7b + 49c + 343d), 1_000_000_007)
        end
        return (nodes=length(tags), tets=length(tets) ÷ 4, checksum=checksum)
    finally
        API.finalize()
    end
end

occ_box() = () -> begin
    API.model.occ.add_box(0, 0, 0, 1, 1, 1)
    API.model.occ.synchronize()
end

# The references are Gmsh's on x86-64. The surface mesh, the inserted points
# and the optimizer's quality comparisons come out of long floating-point
# chains that differ in the last bits on other architectures (AArch64 fuses
# multiply-adds), so there the counts are held to a band only.
function check_mesh(result, nodes, tets, checksum)
    if Sys.ARCH === :x86_64
        @test result == (nodes=nodes, tets=tets, checksum=checksum)
    else
        @test abs(result.nodes - nodes) <= max(3, nodes ÷ 50)
        @test abs(result.tets - tets) <= max(10, tets ÷ 25)
    end
end

@testset "TetGen/BR region pass and refinement in Gmsh's order" begin
    # Initial mesh: flips of dihedral angles above 179 degrees after the
    # exterior of the enclosing box is carved.
    check_mesh(mesh_checksum(occ_box(); initial=true), 272, 769, 993886645)
    # Refinement (breadth-first cavities, Gmsh's shell and cell order) and the
    # optimizer on that initial mesh.
    check_mesh(mesh_checksum(occ_box()), 339, 1125, 287225513)
end

end
