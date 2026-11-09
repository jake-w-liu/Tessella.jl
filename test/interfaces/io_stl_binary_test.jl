module StlBinaryChecks
using Test,Tessella
const F=Tessella.IO
const A=Tessella.API

function binary_facets(path)
    bytes=read(path)
    stream=IOBuffer(bytes)
    header=read(stream,80)
    count=Int(ltoh(read(stream,UInt32)))
    @test length(bytes)==84+50count
    records=map(1:count) do _
        values=ntuple(_->Float64(reinterpret(Float32,ltoh(read(stream,UInt32)))),12)
        attribute=ltoh(read(stream,UInt16))
        (normal=values[1:3],points=ntuple(i->values[3i+1:3i+3],3),attribute=attribute)
    end
    @test eof(stream)
    return header,records
end

function triangle(;scale=1.)
    Mesh([0. scale 0.;0. 0. scale;0. 0. 0.];tris=reshape(Int32[1,2,3],3,1))
end

function normal_checksum(points,count)
    result=0.
    for _ in 1:count
        normal=F._stl_unit_normal(points...)
        result+=sum(normal)
    end
    return result
end

@testset "native binary STL preserves its literal facet layout" begin
    mktempdir() do dir
        coordinates=[0. 1. 1. 0. 2.;0. 0. 1. 1. 0.;0. 0. 0. 0. 0.]
        mesh=MixedMesh(coordinates,[
            ElementBlock(2,reshape(Int32[2,5,3],3,1)),
            ElementBlock(3,reshape(Int32[1,2,3,4],4,1)),
            ElementBlock(1,reshape(Int32[1,2],2,1))])
        path=joinpath(dir,"facets.stl")
        @test F.write_stl(path,mesh;binary=true)==path
        header,records=binary_facets(path)
        @test startswith(String(header),"Tessella binary STL")
        @test length(records)==3
        @test all(record->record.attribute==0 && record.normal==(0.,0.,1.),records)
        expected=[ntuple(i->Tuple(coordinates[:,cell[i]]),3) for cell in
            ((2,5,3),(1,2,3),(1,3,4))]
        @test getproperty.(records,:points)==expected
        decoded=F.read_stl(path;merge_tol=0.)
        @test size(decoded.tris,2)==3
        @test validate(decoded).ok
        @test Set(Tuple(Tuple(decoded.coords[:,v]) for v in cell) for cell in eachcol(decoded.tris))==Set(expected)
        @test mesh.coords==coordinates
        @test F._stl_binary_count(UInt64(typemax(UInt32)))==typemax(UInt32)
        for count in (0,-1,true,UInt64(typemax(UInt32))+1,typemax(UInt64))
            @test_throws ArgumentError F._stl_binary_count(count)
        end
        @test_throws ArgumentError F.write_stl(path,mesh;binary=1)
        @test length(read(path))==234
    end
end

@testset "STL normals retain scale and winding" begin
    mktempdir() do dir
        for scale in (1e-300,1e-100,1.,1e100,1e300)
            mesh=triangle(;scale)
            path=joinpath(dir,"normal.stl")
            F.write_stl(path,mesh)
            line=only(filter(line->occursin("facet normal",line),readlines(path)))
            @test parse.(Float64,split(line)[3:5])==[0.,0.,1.]
            points=ntuple(i->Tuple(mesh.coords[:,i]),3)
            normal_checksum(points,2)
            if Base.JLOptions().compile_enabled!=3
                @test (@allocated normal_checksum(points,100))==0
            end
        end
        @test F._stl_unit_normal((1e308,0.,0.),(-1e308,0.,0.),(1e308,1e308,0.))==(0.,0.,-1.)
        @test F._stl_unit_normal((0.,0.,0.),(1.,1.,1.),(2.,2.,2.))==(0.,0.,0.)
        reversed=Mesh([0. 1. 0.;0. 0. 1.;0. 0. 0.];tris=reshape(Int32[1,3,2],3,1))
        path=joinpath(dir,"reversed.stl")
        F.write_stl(path,reversed;binary=true)
        _,records=binary_facets(path)
        @test only(records).normal==(0.,0.,-1.)
    end
end

@testset "binary STL failures preserve the destination" begin
    mktempdir() do dir
        path=joinpath(dir,"existing.stl")
        for mesh in (triangle(;scale=1e100),triangle(;scale=1e-100),
                Mesh(zeros(3,0)),
                Mesh([0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1.];
                    tets=reshape(Int32[1,2,3,4],4,1)))
            write(path,"held destination")
            @test_throws ArgumentError F.write_stl(path,mesh;binary=true)
            @test read(path,String)=="held destination"
            @test readdir(dir)==["existing.stl"]
        end
        corrupted=MixedMesh([0. 1. 0.;0. 0. 1.;0. 0. 0.],[ElementBlock(2,reshape(Int32[1,2,3],3,1))])
        corrupted.blocks[1].nodes[1,1]=99
        @test_throws ArgumentError F.write_stl(path,corrupted;binary=true)
        @test read(path,String)=="held destination"
    end
end

@testset "API, GEO Save and CLI transport binary STL selection" begin
    mktempdir() do dir
        ascii=joinpath(dir,"source.stl")
        F.write_stl(ascii,triangle())
        source=read(ascii)
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.open(ascii)
            nodes=A.mesh.get_nodes();elements=A.mesh.get_elements()
            A.option("Mesh.Binary",1)
            out=joinpath(dir,"api.stl")
            A.write(out)
            _,records=binary_facets(out)
            @test length(records)==1
            @test only(records).points==((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))
            @test A.mesh.get_nodes()==nodes && A.mesh.get_elements()==elements
            A.option("Mesh.Binary",0)
            A.write(out)
            @test startswith(read(out,String),"solid tessella")
        finally
            A.finalize()
        end
        geo=joinpath(dir,"source.geo");geo_out=joinpath(dir,"geo.stl")
        write(geo,"Mesh.Binary=1; Merge \"$(replace(ascii,'\\'=>'/'))\"; Save \"$(replace(geo_out,'\\'=>'/'))\";")
        executed=execute_geo(geo;mesh_dim=0)
        @test isempty(executed.errors)
        _,records=binary_facets(geo_out)
        @test length(records)==1
        cli_out=joinpath(dir,"cli.stl")
        @test Tessella.CLI.main([geo,"-2","-bin","-o",cli_out])==cli_out
        _,records=binary_facets(cli_out)
        @test length(records)==1
        @test read(ascii)==source
    end
end
end
