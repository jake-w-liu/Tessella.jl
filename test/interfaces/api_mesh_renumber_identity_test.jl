module ApiMeshRenumberIdentityTest

using Test, Tessella

const API = Tessella.API

# renumber_nodes()/renumber_elements() right after generate: the mesh is
# already numbered in entity order, so Gmsh 4.15.2 changes nothing and the
# tag counters keep their values (98/499 here, 99 with a physical group,
# whose partial-save parking generate already applied). Tessella takes the
# no-copy path for this case; the observable state must be the same.
function square(physical)
    API.initialize()
    API.option("General.Terminal", 0)
    g = API.model.geo
    p = [g.add_point(x, y, 0, 0) for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1))]
    l = [g.add_line(p[i], p[mod1(i + 1, 4)]) for i in 1:4]
    s = g.add_plane_surface([g.add_curve_loop(l)])
    g.synchronize()
    physical && API.model.add_physical_group(2, [s]; tag=7)
    API.option("Mesh.MeshSizeMax", 0.5)
    API.mesh.generate(2)
end

snapshot() = (API.mesh.get_nodes()[1:2], API.mesh.get_elements()[2:3],
              API.mesh.get_max_node_tag(), API.mesh.get_max_element_tag())

@testset "renumber after generate is the identity, as in Gmsh" begin
    for physical in (false, true)
        try
            square(physical)
            before = snapshot()
            API.mesh.renumber_nodes()
            @test snapshot() == before
            API.mesh.renumber_elements()
            @test snapshot() == before
            # Gmsh 4.15.2's counters for this square.
            if Sys.ARCH === :x86_64
                @test length(before[1][1]) == 98
                @test before[3] == (physical ? 99 : 98)
                @test before[4] == 499
            end
            # A real permutation still goes through the full commit.
            tags = before[1][1]
            first_xyz = API.mesh.get_node(tags[1])[1]
            API.mesh.renumber_nodes(tags, reverse(tags))
            @test API.mesh.get_node(tags[end])[1] == first_xyz
            @test sort(API.mesh.get_nodes()[1]) == sort(tags)
        finally
            API.finalize()
        end
    end
end

end
