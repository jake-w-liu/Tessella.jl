module MshOpenGmshParityTest

using Test, Tessella

# Gmsh 4.15.2 opens every file below; each expectation was taken from Gmsh itself
# (gmsh.open + the matching gmsh.model query on the same text).
const API = Tessella.API

const NEGATIVE_TET = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 0 1
1 0 0 0 1 1 1 1 1 0
\$EndEntities
\$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
\$EndNodes
\$Elements
1 1 1 1
3 1 4 1
1 1 3 2 4
\$EndElements
"""

# Two tetrahedra, the second with its own nodes lying on the first one's face
# (a hanging node) and one inverted element.
const HANGING = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 0 1
1 0 0 -1 1 1 1 1 1 0
\$EndEntities
\$Nodes
1 7 1 7
3 1 0 7
1
2
3
4
5
6
7
0 0 0
1 0 0
0 1 0
0 0 1
1 1 1
0.25 0.25 0
0.3 0.2 -1
\$EndNodes
\$Elements
1 2 1 2
3 1 4 2
1 1 2 3 4
2 1 5 6 7
\$EndElements
"""

# A surface whose boundary lists its two curves twice (seam style), with a
# Curve that carries no element.
const SEAM = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
3 2 1 0
1 0 0 0 0
2 1 0 0 0
3 0 1 0 0
1 0 0 0 1 0 0 0 2 -1 2
2 0 0 0 0 1 0 0 2 -3 1
1 0 0 0 1 1 0 0 4 1 2 -1 -2
\$EndEntities
\$Nodes
4 3 1 3
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
0 1 0
2 1 0 0
\$EndNodes
\$Elements
1 1 1 1
2 1 2 1
1 1 2 3
\$EndElements
"""

const PERIODIC = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
4 2 0 0
1 0 0 0 0
2 0 1 0 0
3 1 0 0 0
4 1 1 0 0
1 0 0 0 0 1 0 0 2 -1 2
2 1 0 0 1 1 0 0 2 -3 4
\$EndEntities
\$Nodes
2 4 1 4
1 1 0 2
1
2
0 0 0
0 1 0
1 2 0 2
3
4
1 0 0
1 1 0
\$EndNodes
\$Elements
2 2 1 2
1 1 1 1
1 1 2
1 2 1 1
2 3 4
\$EndElements
\$Periodic
3
0 3 1
0
1
3 1
0 4 2
0
1
4 2
1 2 1
0
2
3 1
4 2
\$EndPeriodic
"""

# Old-style partitioned MSH 2.2: tags are physical, elementary, #partitions, partition.
partitioned(version, tags) = """
\$MeshFormat
$version 0 8
\$EndMeshFormat
\$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
\$EndNodes
\$Elements
2
1 2 $(length(tags[1])) $(join(tags[1], ' ')) 1 2 3
2 2 $(length(tags[2])) $(join(tags[2], ' ')) 1 3 4
\$EndElements
"""

function with_file(f, text)
    path = tempname() * ".msh"
    write(path, text)
    try
        return f(path)
    finally
        rm(path; force = true)
    end
end

function opened(f, text)
    with_file(text) do path
        API.initialize()
        try
            API.open(path)
            return f()
        finally
            API.finalize()
        end
    end
end

@testset "A.open accepts what Gmsh opens" begin
    @testset "inverted and hanging tetrahedra" begin
        opened(NEGATIVE_TET) do
            types, tags, nodes = API.mesh.get_elements(3)
            @test types == [4] && nodes == [[1, 3, 2, 4]]
        end
        opened(HANGING) do
            _, tags, nodes = API.mesh.get_elements(3)
            @test length(tags[1]) == 2 && length(API.mesh.get_nodes()[1]) == 7
            out = tempname() * ".msh"
            API.option("Mesh.MshFileVersion", 4.1)
            API.write(out)
            API.initialize()
            API.open(out)
            @test API.mesh.get_elements(3)[3] == nodes
            rm(out)
        end
    end

    @testset "format versions 2.0, 2.1 and 4.2" begin
        for version in ("2.0", "2.1", "2.2")
            # Before 2.2 the third tag is the partition number itself.
            tags = version == "2.2" ? ([1, 1, 1, 1], [1, 1, 1, 2]) : ([1, 1, 1], [1, 1, 2])
            opened(partitioned(version, tags)) do
                @test API.model.get_number_of_partitions() == 2
            end
        end
        opened(replace(NEGATIVE_TET, "4.1 0 8" => "4.2 0 8")) do
            @test length(API.mesh.get_nodes()[1]) == 4
        end
    end

    @testset "old-style MSH 2 partitions become partition entities" begin
        opened(partitioned("2.2", ([1, 1, 1, 1], [1, 1, 1, 2]))) do
            @test API.model.get_number_of_partitions() == 2
            @test API.model.get_partitions(2, 2) == [1]
            @test API.model.get_partitions(2, 3) == [2]
            @test API.model.get_parent(2, 2) == (2, 1)
            @test API.mesh.get_elements(2, 2)[2] == [UInt64[1]]
            @test API.mesh.get_elements(2, 3)[2] == [UInt64[2]]
        end
    end

    @testset "boundary entities are listed once" begin
        opened(SEAM) do
            @test [b[2] for b in API.model.get_boundary([(2, 1)], false, false, false)] == [1, 2]
        end
    end

    @testset "bounding boxes of loaded entities" begin
        opened(SEAM) do
            @test API.model.get_bounding_box(0, 3) == (0.0, 1.0, 0.0, 0.0, 1.0, 0.0)
            @test API.model.get_bounding_box(2, 1) == (0.0, 0.0, 0.0, 1.0, 1.0, 0.0)
            # Gmsh: a Curve without elements has an empty box.
            err = try API.model.get_bounding_box(1, 1); nothing catch e; e end
            @test err isa ArgumentError && occursin("Empty bounding box", err.msg)
            @test API.model.get_bounding_box(-1, -1) == (0.0, 0.0, 0.0, 1.0, 1.0, 0.0)
            @test API.model.get_entities_in_bounding_box(-1, -1, -1, 2, 2, 2, 1) == [(1, 1), (1, 2)]
            @test API.model.get_entities_in_bounding_box(-1, -1, -1, 2, 2, 2, 0) ==
                  [(0, 1), (0, 2), (0, 3)]
        end
    end

    @testset "periodic Points and Curves of a loaded file" begin
        opened(PERIODIC) do
            master, slaves, masters, affine = API.mesh.get_periodic_nodes(0, 3, false)
            @test master == 1 && slaves == [3] && masters == [1] && isempty(affine)
            master, slaves, masters, _ = API.mesh.get_periodic_nodes(1, 2, false)
            @test master == 1
            @test sort(collect(zip(slaves, masters))) == [(3, 1), (4, 2)]
            @test API.mesh.get_periodic_nodes(1, 1, false)[1] == 1
        end
    end

    @testset "Gmsh-shaped model data" begin
        API.initialize()
        try
            API.model.add("m")
            API.model.geo.add_point(0, 0, 0; meshSize = 1.0)
            API.model.geo.add_point(1, 0, 0; meshSize = 1.0)
            API.model.geo.add_line(1, 2)
            API.model.geo.synchronize()
            API.mesh.generate(1)
            nodes = API.mesh.get_nodes()[1]
            view = API.view.add("v")
            API.view.add_model_data(view, 0, "m", "NodeData", nodes,
                [[Float64(i), 2.0i] for i in eachindex(nodes)], 0.5, 2)
            kind, tags, data, time, components = API.view.get_model_data_full(view, 0)
            @test kind == "NodeData" && tags == nodes && time == 0.5 && components == 2
            @test data == [[Float64(i), 2.0i] for i in eachindex(nodes)]
            @test API.view.get_model_data(view, 0)[3] == reduce(vcat, data)
        finally
            API.finalize()
        end
    end
end

end # module
