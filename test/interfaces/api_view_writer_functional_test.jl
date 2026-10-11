module ViewWriterFunctionalChecks
using Test, Tessella
const A = Tessella.API

function triangle_model(name, first_node, element)
    A.model.add(name)
    A.model.add_discrete_entity(2, 7)
    tags = collect(first_node:first_node + 2)
    A.mesh.add_nodes(2, 7, tags, [0., 0., 0., 1., 0., 0., 0., 1., 0.])
    A.mesh.add_elements_by_type(7, 2, [element], tags)
    return tags
end

@testset "view model rows, numbering, export and stream cleanup" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.option("General.Verbosity", 99)
        tags = triangle_model("first", 101, 301)
        view = A.view.add("scalar")
        A.view.add_model_data(view, 0, "first", "NodeData", [tags[1]], [[1.25]])
        A.view.add_model_data(view, 0, "first", "NodeData", [tags[2]], [[2.5]])
        A.view.add_model_data(view, 2, "first", "NodeData", [tags[3]], [[3.75]])
        @test A.view.get_model_data(view, 0) == (["NodeData"], UInt64[101, 102], [1.25, 2.5])
        @test A.view.get_model_data(view, 1) == (["NodeData"], UInt64[], Float64[])
        A.mesh.renumber_nodes(tags, [1, 2, 3])
        @test A.view.get_model_data(view, 0) == (["NodeData"], UInt64[1, 2], [1.25, 2.5])
        @test A.mesh.get_max_node_tag() == 103

        mktempdir() do directory
            path = joinpath(directory, "single.msh")
            A.logger.start()
            A.view.write(view, path)
            logs = A.logger.get()
            A.logger.stop()
            @test !any(message -> occursin("dense node numbering", message), logs)
            written = Tessella.read_mixed_msh(path)
            @test size(written.coords, 2) == 3
            @test [section.values for section in written.data_sections] == [[1.25, 2.5], [3.75]]
            bytes = read(path)
            A.view.write(view, path, true)
            @test read(path) == vcat(bytes, bytes)

            triangle_model("second", 201, 401)
            A.view.add_model_data(view, 3, "second", "NodeData", [201], [[4.5]])
            multi = joinpath(directory, "multi.msh")
            A.logger.start()
            A.view.write(view, multi)
            logs = A.logger.get()
            A.logger.stop()
            exporting = findfirst(==("Info: Exporting multi-mesh view in separate files"), logs)
            rebuilding = findfirst(message -> occursin("Rebuilding mesh node cache", message), logs)
            @test exporting !== nothing && rebuilding !== nothing && exporting < rebuilding
            @test !isfile(multi)
            @test [section.values for section in Tessella.read_mixed_msh(joinpath(directory, "multi_0.msh")).data_sections] ==
                [[1.25, 2.5], [3.75]]
            @test only(Tessella.read_mixed_msh(joinpath(directory, "multi_1.msh")).data_sections).values == [4.5]

            bad = A.view.add("orphan row")
            A.view.add_model_data(bad, 0, "second", "NodeData", [999], [[9.]])
            failure = joinpath(directory, "failure.msh")
            before = A.mesh.get_nodes()
            @test_throws ArgumentError A.view.write(bad, failure)
            @test A.mesh.get_nodes() == before
            @test occursin("\$NodeData", read(failure, String))
            rm(failure)
            @test !ispath(failure)
        end
    finally
        A.finalize()
    end
end

@testset "quadratic pyramid interpolation exports with its model data" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.model.add("pyramid")
        A.model.add_discrete_entity(3, 7)
        _, dimension, _, count, coordinates, _ = A.mesh.get_element_properties(14)
        @test dimension == 3 && count == 14
        tags = collect(101:100 + count)
        A.mesh.add_nodes(3, 7, tags, coordinates)
        A.mesh.add_elements_by_type(7, 14, [301], tags)
        view = A.view.add("pyramid values")
        values = collect(1.:count)
        A.view.add_model_data(view, 0, "pyramid", "ElementNodeData", [301], [values])
        mktempdir() do directory
            path = joinpath(directory, "pyramid.msh")
            A.view.write(view, path)
            @test occursin("\$InterpolationScheme", read(path, String))
            section = only(Tessella.read_mixed_msh(path).data_sections)
            @test section.row_nodes == Int32[14] && section.values == values
            A.option("PostProcessing.Binary", 1)
            A.view.write(view, path)
            @test only(Tessella.read_mixed_msh(path).data_sections).values == values
        end
    finally
        A.finalize()
    end
end

@testset "caller-owned records survive a covered cache's mutations" begin
    A.initialize()
    try
        A.option("General.Terminal", 0)
        mktempdir() do directory
            geometry=joinpath(directory,"triangle.geo")
            write(geometry,"Point(1)={0,0,0,1}; Point(2)={1,0,0,1}; Point(3)={0,1,0,1};\n" *
                "Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1}; Curve Loop(1)={1,2,3};\n" *
                "Plane Surface(1)={1}; Transfinite Curve{1,2,3}=2; Transfinite Surface{1};\n")
            A.open_geo!(geometry;mesh_dim=0)
            A.mesh.generate(2)
            corner=first(first(A.mesh.get_elements(2,1)[3]))
            extra_coordinates=[2.,0.,0.,2.,1.,0.]
            parameters=A.model.get_parametrization(2,1,extra_coordinates)
            A.mesh.add_nodes(2,1,[900,901],extra_coordinates,parameters)
            A.mesh.add_elements_by_type(1,2,[901],[corner,900,901])
            @test 900 in first(A.mesh.get_nodes())
            @test 901 in reduce(vcat,A.mesh.get_elements(2,1)[2];init=UInt64[])
            A.mesh.set_node(901,[2.,1.,2.])
            A.mesh.relocate_nodes(2,1)
            @test first(A.mesh.get_node(901))≈[2.,1.,0.] atol=1e-12
            A.mesh.set_node(900,[3.,0.,0.])
            @test first(A.mesh.get_node(900))==[3.,0.,0.]
            view=A.view.add("added cell")
            A.view.add_model_data(view,0,"","ElementData",[901],[[2.]])
            path=joinpath(directory,"added.msh")
            A.view.write(view,path)
            @test only(Tessella.read_mixed_msh(path).data_sections).values==[2.]
            cell_tags=first(A.mesh.get_elements_by_type(2,1))
            @test_throws ArgumentError A.mesh.reorder_elements(2,1,fill(0,length(cell_tags)))
            @test first(A.mesh.get_elements_by_type(2,1))==cell_tags
            A.mesh.reorder_elements(2,1,collect(length(cell_tags)-1:-1:0))
            @test first(A.mesh.get_elements_by_type(2,1))==reverse(cell_tags)
            A.mesh.renumber_nodes([900],[1000])
            @test first(A.mesh.get_node(1000))==[3.,0.,0.]
            A.mesh.affine_transform([1.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.,0.])
            @test first(A.mesh.get_node(1000))==[4.,0.,0.]
            A.mesh.remove_elements(2,1,[901])
            @test !(901 in reduce(vcat,A.mesh.get_elements(2,1)[2];init=UInt64[]))
            @test first(A.mesh.get_node(1000))==[4.,0.,0.]
        end
    finally
        A.finalize()
    end
end

@testset "list POS layout and PostView field interpolation" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        point=A.view.add("point")
        A.view.add_list_data(point,"SP",1,[2.,3.,4.,7.])
        triangle=A.view.add("linear triangle")
        blocked=[0.,2.,0.,0.,0.,3.,0.,0.,0.,1.,3.,4.]
        A.view.add_list_data(triangle,"ST",1,blocked)
        field_tag=A.mesh.field.add("PostView")
        A.mesh.field.set_number(field_tag,"ViewTag",triangle)
        A.mesh.field.set_as_background_mesh(field_tag)
        field=A._session_size_field_locked(A.CURRENT[];raw=true)
        @test Tessella.field_value(field,.5,.75,0.)≈2.25 atol=1e-12
        for invalid in (2.5,NaN,1e100)
            A.mesh.field.set_string(field_tag,"ViewTag",string(invalid))
            @test_throws ArgumentError A._session_size_field_locked(A.CURRENT[];raw=true)
        end
        A.mesh.field.set_number(field_tag,"ViewTag",-1)
        A.mesh.field.set_number(field_tag,"ViewIndex",0)
        @test Tessella.field_value(A._session_size_field_locked(A.CURRENT[];raw=true),2.,3.,4.)==7.
        for invalid in (1.5,NaN,1e100)
            A.mesh.field.set_string(field_tag,"ViewIndex",string(invalid))
            @test_throws ArgumentError A._session_size_field_locked(A.CURRENT[];raw=true)
        end
        A.mesh.field.set_number(field_tag,"ViewIndex",1)
        @test Tessella.field_value(A._session_size_field_locked(A.CURRENT[];raw=true),.5,.75,0.)≈2.25 atol=1e-12
        @test A.view.get_list_data(triangle)==(["ST"],[blocked])
        tetrahedron=A.view.add("tetrahedron")
        A.view.add_list_data(tetrahedron,"SS",1,
            [0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,2.,3.,4.])
        tetrahedral_field=A._post_view_field(A.SESSION_VIEWS[Int(tetrahedron)],2)
        @test Tessella.field_value(tetrahedral_field,.1,.2,.3)≈2.4 atol=1e-12
        mktempdir() do directory
            path=joinpath(directory,"list.pos")
            A.view.write(point,path)
            element=only(only(Tessella.Post.read_pos(path)).elements)
            @test element.kind==:point && vec(element.coords)==[2.,3.,4.]
            @test vec(element.values)==[7.]
            A.view.write(triangle,path)
            element=only(only(Tessella.Post.read_pos(path)).elements)
            @test element.coords==[0. 2. 0.;0. 0. 3.;0. 0. 0.]
            @test vec(element.values)==[1.,3.,4.]
            A.view.write(point,path,true)
            @test length(Tessella.Post.read_pos(path))==2
            unsupported=joinpath(directory,"unsupported.msh")
            @test_throws ArgumentError A.view.write(point,unsupported)
            @test !ispath(unsupported)
            A.view.add_list_data(point,"ST",1,[1.,2.])
            before=read(path)
            @test_throws ArgumentError A.view.write(point,path)
            @test read(path)==before
        end

        A.model.add("cross layout")
        A.model.add_discrete_entity(2,7)
        A.mesh.add_nodes(2,7,[11,12,13],[0.,0.,0.,2.,0.,0.,0.,3.,0.])
        A.mesh.add_elements_by_type(7,2,[31],[11,12,13])
        views=A.mesh.compute_cross_field()
        for view in views
            _,data=A.view.get_list_data(view)
            @test first(data)[1:9]==blocked[1:9]
        end
    finally
        A.finalize()
    end
end

@testset "outward orientation writes covered caller cells" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            geometry=joinpath(directory,"tetrahedron.geo")
            write(geometry,"""
                Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
                Point(3)={0,1,0,1}; Point(4)={0,0,1,1};
                Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};
                Line(4)={1,4}; Line(5)={2,4}; Line(6)={3,4};
                Curve Loop(1)={1,2,3}; Plane Surface(1)={1};
                Curve Loop(2)={1,5,-4}; Plane Surface(2)={2};
                Curve Loop(3)={2,6,-5}; Plane Surface(3)={3};
                Curve Loop(4)={3,4,-6}; Plane Surface(4)={4};
                Surface Loop(1)={-1,2,3,4}; Volume(1)={1};
                Transfinite Curve{1:6}=2; Transfinite Surface{1:4};
                """)
            A.open_geo!(geometry;mesh_dim=0)
            A.mesh.generate(2)
            cached_tag=first(first(A.mesh.get_elements(2,1)[2]))
            nodes=A.mesh.get_element(cached_tag)[2]
            A.mesh.add_elements_by_type(1,2,[9001],nodes)
            before=A.mesh.get_element(9001)[2]
            A.mesh.set_outward_orientation(1)
            @test A.mesh.get_element(9001)[2]==before[[1,3,2]]
            @test A.mesh.get_element(cached_tag)[2]==before[[1,3,2]]
        end
    finally
        A.finalize()
    end
end
end
