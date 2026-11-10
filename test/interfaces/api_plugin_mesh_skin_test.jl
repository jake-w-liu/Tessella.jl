module PluginMeshSkinChecks
using Test,Tessella
const A=Tessella.API
function setup(dim,typ,nodes;coords=[0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.,0.,0.,1.,1.,0.,1.,1.,1.,1.,0.,1.,1.])
    A.initialize()
    A.option("General.Terminal",0)
    A.model.add_discrete_entity(dim,17)
    A.mesh.add_nodes(dim,17,collect(1:length(coords)÷3),coords)
    width=Tessella.Elements.msh_num_nodes(typ)
    A.mesh.add_elements_by_type(17,typ,collect(1:length(nodes)÷width),nodes)
    A.plugin.set_number("Skin","FromMesh",1)
end
@testset "Skin publishes actual discrete mesh boundaries" begin
    try
        setup(2,2,[1,2,3])
        source=A.mesh.get_elements(2,17)
        before=A.mesh.get_nodes()
        @test A.plugin.run("Skin")===nothing
        @test A.view.get_tags()==Int32[]
        @test A.model.get_entities(1)==[(1,1)]
        types,uids,cells=A.mesh.get_elements(1,1)
        @test types==[1]
        @test uids==[[2,3,4]]
        @test cells==[[1,2,3,1,2,3]]
        @test A.mesh.get_elements(2,17)==source
        @test A.mesh.get_nodes()==before
        @test isempty(A.mesh.get_nodes(1,1)[1])
        @test isempty(A.model.get_boundary([(1,1)]))
        A.plugin.run("Skin")
        @test A.model.get_entities(1)==[(1,1),(1,2)]
        @test A.mesh.get_elements(1,2)[2]==[[5,6,7]]
    finally
        A.finalize()
    end
end

# Expected oriented boundary records captured through the public Gmsh 4.15.2
# API. No reference mesher is needed by these regression tests.
@testset "Skin family factories retain oriented primary facets" begin
    cases=((3,[1,2,3,4],[1],[[1,2,4,1,2,3,3,4]]),
           (4,[1,2,3,4],[2],[[1,3,2,1,2,4,1,4,3,4,2,3]]),
           (5,collect(1:8),[3],[[1,4,3,2,1,2,6,5,1,5,8,4,2,3,7,6,3,4,8,7,5,6,7,8]]),
           (6,collect(1:6),[2,3],[[1,3,2,4,5,6],[1,2,5,4,1,4,6,3,2,3,6,5]]),
           (7,collect(1:5),[2,3],[[1,2,5,4,1,5,2,3,5,3,4,5],[1,4,3,2]]),
           (9,[1,2,3,5,6,4],[1],[[1,2,3,1,2,3]]),
           (140,collect(1:4),[2,3],[[1,4,2,2,4,3],[1,2,3,4]]))
    for (typ,nodes,expected_types,expected_cells) in cases
        dimension=Int(Tessella.Elements.msh_spec(typ).dim)
        try
            setup(dimension,typ,nodes)
            source=A.mesh.get_elements(dimension,17)
            A.plugin.run("Skin")
            types,uids,cells=A.mesh.get_elements(dimension-1,1)
            @test types==expected_types
            @test cells==expected_cells
            @test vcat(uids...)==collect(2:1+sum(length(c)÷Tessella.Elements.msh_num_nodes(t) for (t,c) in zip(types,cells)))
            @test A.mesh.get_elements(dimension,17)==source
            @test isempty(A.mesh.get_nodes(dimension-1,1)[1])
        finally
            A.finalize()
        end
    end
end

@testset "Skin cancellation, visibility, and identity" begin
    for separate in (false,true)
        try
            setup(2,2,[1,2,3])
            owner=separate ? A.model.add_discrete_entity(2,18) : 17
            A.mesh.add_elements_by_type(owner,2,[2],[2,4,3])
            # Element visibility does not filter the entity's boundary scan.
            A.mesh.set_visibility([1],0)
            A.plugin.set_number("Skin","View",123456.)
            A.plugin.run("Skin")
            @test A.mesh.get_elements(1,1)==([1],[[3,4,5,6]],[[1,2,3,1,2,4,4,3]])
        finally
            A.finalize()
        end
    end
    for visible in (0.,1.)
        try
            setup(2,2,[1,2,3])
            A.model.set_visibility([(2,17)],0)
            A.plugin.set_number("Skin","Visible",visible)
            A.plugin.run("Skin")
            @test A.model.get_entities(1)==[(1,1)]
            types,_,cells=A.mesh.get_elements(1,1)
            @test visible==0 ? types==[1] && cells==[[1,2,3,1,2,3]] : isempty(types)
        finally
            A.finalize()
        end
    end
    try
        setup(2,2,[1,2,3])
        A.mesh.add_nodes(2,17,[9,10,11],[0.,0.,0.,1.,0.,0.,0.,1.,0.])
        A.mesh.add_elements_by_type(17,2,[2],[9,10,11])
        A.plugin.run("Skin")
        @test A.mesh.get_elements(1,1)[3]==[[1,2,3,1,2,3,9,10,11,9,10,11]]
    finally
        A.finalize()
    end
    try
        setup(2,2,[1,2,3,1,3,2])
        A.plugin.run("Skin")
        @test A.model.get_entities(1)==[(1,1)]
        @test isempty(A.mesh.get_elements(1,1)[1])
    finally
        A.finalize()
    end
    try
        setup(2,2,[1,2,3])
        A.model.add_discrete_entity(3,23)
        A.plugin.run("Skin")
        @test A.model.get_entities(2)==[(2,17),(2,18)]
        @test isempty(A.mesh.get_elements(2,18)[1])
        @test isempty(A.model.get_entities(1))
    finally
        A.finalize()
    end
    for populated in (false,true)
        try
            A.initialize();A.option("General.Terminal",0)
            if populated
                A.model.add_discrete_entity(1,17)
                A.mesh.add_nodes(1,17,[1,2],[0.,0.,0.,1.,0.,0.])
                A.mesh.add_elements_by_type(17,1,[1],[1,2])
            end
            before=A.model.get_entities()
            A.plugin.set_number("Skin","FromMesh",1)
            A.plugin.run("Skin")
            @test A.model.get_entities()==before
            @test isempty(A.view.get_tags())
        finally
            A.finalize()
        end
    end
end

@testset "Skin attaches boundaries beside native imported cache" begin
    mktempdir() do dir
        path=joinpath(dir,"skin-input.msh")
        write(path,"""\$MeshFormat
2.2 0 8
\$EndMeshFormat
\$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 1 1 0
\$EndNodes
\$Elements
2
41 2 2 0 17 1 2 3
42 2 2 0 17 2 4 3
\$EndElements
""")
        try
            A.initialize();A.option("General.Terminal",0);A.open(path)
            @test A.LAST_MESH[]!==nothing
            source=A.mesh.get_elements(2,17)
            nodes=A.mesh.get_nodes()
            cache=A.LAST_MESH[]
            A.plugin.set_number("Skin","FromMesh",1)
            A.plugin.run("Skin")
            @test A.LAST_MESH[]===cache
            @test A.mesh.get_elements(2,17)==source
            @test A.mesh.get_nodes()==nodes
            @test A.mesh.get_elements(1,1)==([1],[[43,44,45,46]],[[1,2,3,1,2,4,4,3]])
            @test isempty(A.mesh.get_nodes(1,1)[1])
            result=joinpath(dir,"skin-output.msh")
            A.write(result)
            @test isfile(result)
            A.open(result)
            @test A.mesh.get_elements(1,1)[3]==[[1,2,3,1,2,4,4,3]]
        finally
            A.finalize()
        end
    end
end
end
