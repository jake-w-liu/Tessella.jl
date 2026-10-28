@testset "public API reference queries, session state and mesh lifecycle" begin
    api = Tessella.API
    api.finalize()
    @test_throws ArgumentError api.mesh.get()
    @test_throws ArgumentError api.model.add_point(0, 0, 0)

    coordinates, weights = api.mesh.get_integration_points(2, "Gauss2")
    points = reshape(coordinates, 3, :)
    @test sum(weights) ≈ 1 / 2
    @test sum(weights .* points[1, :]) ≈ 1 / 6
    @test sum(weights .* points[2, :] .^ 2) ≈ 1 / 12
    components, values, orientations = api.mesh.get_basis_functions(
        2, [0.2, 0.3, 0.0], "Lagrange")
    @test components == 1 && orientations == 1
    @test values ≈ [0.5, 0.2, 0.3] && sum(values) ≈ 1.0
    @test msh_num_nodes(11) == 10 && msh_dimension(11) == 3 && msh_order(11) == 2

    api.initialize()
    try
        @test api.option("Mesh.MeshSizeFactor", 2.0) == 2.0
        @test_throws ArgumentError api.option("Mesh.MeshSizeFactor", 0.0)
        @test api.option("Mesh.MeshSizeFactor") == 2.0
        api.model.add_point(0, 0, 0; tag=1)
        api.model.add_point(1, 0, 0; tag=2)
        api.model.add_line(1, 2; tag=1)
        api.mesh.set_transfinite_curve(1, 3)
        api.model.add_physical_group(1, [1]; tag=7, name="wire")
        @test api.model.get_entities(1) == [(1, 1)]
        @test api.mesh.set_reverse(3, 999, true) === nothing
        @test api.mesh.set_smoothing(3, 999, 5) === nothing
        @test Set(api.model.get_boundary([(1, 1)], false, false, false)) ==
              Set([(0, 1), (0, 2)])

        api.mesh.generate(1)
        node_tags, flat_coordinates, _ = api.mesh.get_nodes()
        line_tags, line_nodes = api.mesh.get_elements_by_type(1)
        @test length(node_tags) == 3 && length(line_tags) == 2
        @test length(line_nodes) == 4 && Set(line_nodes) ⊆ Set(node_tags)
        # Allow 1e-10 absolute placement error from native curve inversion.
        @test sort(reshape(flat_coordinates, 3, :)[1, :]) ≈ [0.0, 0.5, 1.0] atol=1e-10 rtol=0
        @test length(api.mesh.get_nodes_for_physical_group(1, 7)[1]) == 3
        baseline = api.mesh.get_nodes()
        flat_coordinates[1] += 100
        @test api.mesh.get_nodes() == baseline
        detached = api.mesh.get()
        detached.coords[1, 1] += 100
        @test api.mesh.get_nodes() == baseline
        @test_throws ArgumentError api.mesh.generate(4)
        @test api.mesh.get_nodes() == baseline

        api.model.set_visibility([(1, 1)], 0, true)
        @test api.model.get_visibility(0, 1) == 0
        api.model.set_attribute("owner", ["test"])
        detached_attribute = api.model.get_attribute("owner")
        detached_attribute[1] = "changed"
        @test api.model.get_attribute("owner") == ["test"]
        @test api.mesh.get_nodes() == baseline

        api.mesh.clear()
        @test_throws ArgumentError api.mesh.get()
        api.mesh.generate(1)
        @test length(api.mesh.get_nodes()[1]) == 3
        point_uid = only(api.mesh.get_nodes(0, 2)[1])
        retained_node_tags, retained_coordinates, _ = api.mesh.get_nodes()
        expected_coordinates = reshape(copy(retained_coordinates), 3, :)
        expected_coordinates[:, findfirst(==(point_uid), retained_node_tags)] = [2.0, 0.0, 0.0]
        retained_elements = api.mesh.get_elements()
        api.model.set_coordinates(2, 2, 0, 0)
        @test validate(api.mesh.get()).ok
        @test api.mesh.get_node(point_uid)[1] == [2.0, 0.0, 0.0]
        @test api.mesh.get_nodes()[1:2] == (retained_node_tags, vec(expected_coordinates))
        @test api.mesh.get_elements() == retained_elements
        @test api.model.get_value(0, 2, []) == [2.0, 0.0, 0.0]

        api.model.add("second")
        @test isempty(api.model.get_entities())
        @test isempty(api.model.get_attribute_names())
        api.model.set_current("")
        @test api.model.get_entities(1) == [(1, 1)]

        # MSH4 node tags span the full `size_t` range upstream — verified
        # against live Gmsh 4.15.2 (node tag 5000000003 stores verbatim).
        api.model.add("wide")
        api.model.add_discrete_entity(2, 7)
        wide_nodes = [2_147_483_748, 3_000_000_005, 5_000_000_003]
        api.mesh.add_nodes(2, 7, wide_nodes,
                           [0.0, 0, 0, 1, 0, 0, 0, 1, 0])
        @test api.mesh.get_max_node_tag() == 5_000_000_003
        @test api.mesh.get_nodes(2, 7)[1] == UInt64.(wide_nodes)
        # Element tags route through the upstream signed 32-bit constructor
        # slot (`MElement(int)` sign-extended into size_t): live Gmsh 4.15.2
        # returns these exact widened values for the same inputs.
        api.mesh.add_elements_by_type(7, 2, [777], [2_147_483_748, 3_000_000_005, 5_000_000_003])
        api.mesh.add_elements_by_type(7, 2, [3_000_000_007], [2_147_483_748, 3_000_000_005, 5_000_000_003])
        api.mesh.add_elements_by_type(7, 2, [5_000_000_000], [2_147_483_748, 3_000_000_005, 5_000_000_003])
        _, element_tags, _ = api.mesh.get_elements(2, 7)
        @test reduce(vcat, element_tags) ==
              UInt64[777, 18_446_744_072_414_584_327, 705_032_704]
        @test api.mesh.get_max_element_tag() == 18_446_744_072_414_584_327
    finally
        api.finalize()
    end
    @test_throws ArgumentError api.mesh.get_nodes()
end
