module X3dOrderChecks
using Test,Tessella,LinearAlgebra
const F=Tessella.IO
const E=Tessella.Elements
const A=Tessella.API

function specimen(msh)
    coords=E.lagrange_nodes(msh)
    return E.MixedMesh(coords,[E.ElementBlock(msh,reshape(Int32.(1:size(coords,2)),:,1))])
end

function scene_points(path)
    body=match(r"<Coordinate point=\"(.*?)\"></Coordinate>"s,read(path,String))
    return reshape(parse.(Float64,split(body.captures[1])),3,:)
end

function make_source(dir)
    source=joinpath(dir,"source.msh")
    coords=[-1. 1. 1. -1. 2.;-1. -1. 1. 1. 2.;0. 0. 0. 0. 0.]
    mesh=E.MixedMesh(coords,[E.ElementBlock(3,reshape(Int32[1,2,3,4],4,1)),
        E.ElementBlock(1,reshape(Int32[1,2],2,1)),E.ElementBlock(15,reshape(Int32[5],1,1))];
        elementary_entities=[Int32[1],Int32[2],Int32[3]])
    Tessella.write_mixed_msh(source,mesh;version=2.2)
    return source
end

@testset "X3D retains all triangle and quadrangle corners at every order" begin
    mktempdir() do dir
        path=joinpath(dir,"surface.x3d")
        for family in (:tri,:qua)
            types=unique(E.msh_type(family,p;serendipity=incomplete)
                         for p in 1:10 for incomplete in (false,true))
            for msh in types
                mesh=specimen(msh)
                expected=family===:tri ? [1,2,3] : [1,2,3,1,3,4]
                @test F.write_x3d(path,mesh)==path
                @test scene_points(path)==mesh.coords[:,expected]
                @test isfile(joinpath(dir,"surface_x3d.html"))
            end
        end
        # The quadratic unit reference square has area four, not the two
        # produced by dropping the first-third-fourth corner triangle.
        F.write_x3d(path,specimen(10))
        points=scene_points(path)
        @test sum(norm(cross(points[:,i+1]-points[:,i],points[:,i+2]-points[:,i]))/2
                  for i in 1:3:size(points,2))==4.
        warped=specimen(10)
        warped.coords[3,3]=1.
        F.write_x3d(path,warped)
        @test scene_points(path)==warped.coords[:,[1,2,3,1,3,4]]
    end
end

@testset "X3D source and destination checks preserve held files" begin
    mktempdir() do dir
        path=joinpath(dir,"held.x3d")
        html=joinpath(dir,"held_x3d.html")
        for corrupt in (m->(m.blocks[1].nodes[1,1]=99),
                        m->(m.coords[1,1]=NaN),m->pop!(m.blocks[1].tags))
            mesh=specimen(10)
            corrupt(mesh)
            write(path,"held geometry");write(html,"held launcher")
            before=readdir(dir)
            @test_throws ArgumentError F.write_x3d(path,mesh)
            @test read(path,String)=="held geometry"
            @test read(html,String)=="held launcher"
            @test readdir(dir)==before
        end
        for family in (:tri,:qua)
            write(path,"held geometry");write(html,"held launcher")
            @test_throws ArgumentError F.write_x3d(path,specimen(E.msh_type(family,0)))
            @test read(path,String)=="held geometry"
            @test read(html,String)=="held launcher"
        end
        collision=joinpath(dir,"collision.x3d")
        write(collision,"held geometry")
        html_dir=joinpath(dir,"collision_x3d.html")
        mkdir(html_dir)
        before=readdir(dir)
        @test_throws ArgumentError F.write_x3d(collision,specimen(10))
        @test read(collision,String)=="held geometry"
        @test isdir(html_dir) && isempty(readdir(html_dir))
        @test readdir(dir)==before
        # An unselected order-zero surface does not need polygon corners.
        F.write_x3d(path,specimen(E.msh_type(:tri,0));surfaces=false,edges=true)
        @test !occursin("<TriangleSet",read(path,String))
        @test occursin("<Group DEF=\"edges\">",read(path,String))
    end
end

@testset "public X3D export retains higher-order quadrangle area and source state" begin
    mktempdir() do dir
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.SaveAll",1)
            A.option("Print.X3dSurfaces",1)
            A.model.add("x3d order")
            A.model.add_discrete_entity(2,1)
            msh=E.msh_type(:qua,6)
            coords=E.lagrange_nodes(msh)
            tags=collect(1:size(coords,2))
            A.mesh.add_nodes(2,1,tags,vec(coords))
            A.mesh.add_elements_by_type(1,msh,[1],tags)
            held=(nodes=A.mesh.get_nodes(),elements=A.mesh.get_elements())
            path=joinpath(dir,"public.x3d")
            @test A.write(path)===nothing
            @test scene_points(path)==coords[:,[1,2,3,1,3,4]]
            @test A.mesh.get_nodes()==held.nodes && A.mesh.get_elements()==held.elements
            collision=joinpath(dir,"collision.x3d")
            write(collision,"held geometry")
            mkdir(joinpath(dir,"collision_x3d.html"))
            @test_throws ArgumentError A.write(collision)
            @test read(collision,String)=="held geometry"
            @test A.mesh.get_nodes()==held.nodes && A.mesh.get_elements()==held.elements
        finally
            A.finalize()
        end
    end
end

@testset "canonical X3D API options select geometry and retain lifecycle semantics" begin
    mktempdir() do dir
        A.initialize()
        try
            A.option("General.Terminal",0)
            for (key,default) in (("Print.X3dSurfaces",1.),("Print.X3dEdges",0.),
                                  ("Print.X3dVertices",0.))
                @test A.option(key)==default
                @test A.option(key,2.9)==2.
                for invalid in (true,Inf,NaN)
                    @test_throws ArgumentError A.option(key,invalid)
                    @test A.option(key)==2.
                end
                A.option(key,default)
            end
            A.open(make_source(dir))
            A.option("Mesh.SaveAll",1)
            held=(nodes=A.mesh.get_nodes(),elements=A.mesh.get_elements())
            path=joinpath(dir,"selection.x3d")
            for (surfaces,edges,vertices) in ((1,0,0),(0,1,0),(0,0,1),(1,1,1))
                A.option("Print.X3dSurfaces",surfaces)
                A.option("Print.X3dEdges",edges)
                A.option("Print.X3dVertices",vertices)
                A.write(path)
                body=read(path,String)
                @test (occursin("<TriangleSet",body),occursin("<LineSet",body),
                       occursin("<PointSet",body))==(surfaces!=0,edges!=0,vertices!=0)
            end
            for key in ("Print.X3dSurfaces","Print.X3dEdges","Print.X3dVertices")
                A.option(key,0)
            end
            geometry=read(path);launcher=read(joinpath(dir,"selection_x3d.html"))
            @test_throws ArgumentError A.write(path)
            @test read(path)==geometry && read(joinpath(dir,"selection_x3d.html"))==launcher
            @test A.mesh.get_nodes()==held.nodes && A.mesh.get_elements()==held.elements
            A.clear()
            @test (A.option("Print.X3dSurfaces"),A.option("Print.X3dEdges"),
                   A.option("Print.X3dVertices"))==(0.,0.,0.)
        finally
            A.finalize()
        end
        A.initialize()
        try
            @test (A.option("Print.X3dSurfaces"),A.option("Print.X3dEdges"),
                   A.option("Print.X3dVertices"))==(1.,0.,0.)
        finally
            A.finalize()
        end
    end
end

@testset "GEO Save and CLI setnumber transport canonical X3D selection" begin
    mktempdir() do dir
        source=make_source(dir)
        held=read(source)
        output=joinpath(dir,"geo.x3d")
        geo=joinpath(dir,"source.geo")
        write(geo,"Print.X3dSurfaces=0; Print.X3dEdges=1; Print.X3dVertices=0; Merge \"$(replace(source,'\\'=>'/'))\"; Save \"$(replace(output,'\\'=>'/'))\";")
        executed=Tessella.execute_geo(geo;mesh_dim=0)
        @test isempty(executed.errors)
        @test isfile(output)
        body=read(output,String)
        @test !occursin("<TriangleSet",body) && occursin("<LineSet",body)
        output=joinpath(dir,"cli.x3d")
        @test Tessella.CLI.main([source,"-save","-format","x3d","-save_all",
            "-v","0","-setnumber","General.Terminal","0",
            "-setnumber","Print.X3dSurfaces","0",
            "-setnumber","Print.X3dEdges","0",
            "-setnumber","Print.X3dVertices","1","-o",output])==output
        body=read(output,String)
        @test !occursin("<TriangleSet",body) && occursin("<PointSet",body)
        @test read(source)==held
    end
end
end
