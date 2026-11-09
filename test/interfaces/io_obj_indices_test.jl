module ObjIndexChecks
using Test,Tessella
const F=Tessella.IO
const A=Tessella.API
const PREFIX="v 100 100 100\nv 0 0 0\nv 1 0 0\nv 0 1 0\n"
const FACE=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))

function face_points(path)
    _,flat,records=F._read_obj(path,"OBJ index regression")
    coordinates=reshape(flat,3,:)
    return [Tuple(Tuple(coordinates[:,row]) for row in record.nodes) for record in records]
end

@testset "OBJ relative indices use one-based current vertex positions" begin
    mktempdir() do dir
        path=joinpath(dir,"relative.obj")
        for face in ("f -3 -2 -1", "f 2 -2 -1", "f 2 3 4")
            write(path,PREFIX*face*"\n")
            @test face_points(path)==[FACE]
        end
        write(path,PREFIX*"vt 0 0\nvt 1 0\nvt 0 1\nvn 0 0 1\nf -3/1/1 -2/2/1 -1/3/1\n")
        @test face_points(path)==[FACE]
        write(path,"v 0 0 0\nv 1 0 0\nv 0 1 0\nf -3 -2 -1\nv 100 100 100\n")
        @test face_points(path)==[FACE]
        write(path,"v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf -4 -3 -2 -1\n")
        @test face_points(path)==[((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))]
        for face in ("f 0 2 3", "f -5 -2 -1", "f 2 3 5", "f -999999999999999999999 2 3")
            write(path,PREFIX*face*"\n")
            @test_throws ArgumentError face_points(path)
        end
    end
end

@testset "OBJ API import and merge preserve geometry and failed imports" begin
    mktempdir() do dir
        path=joinpath(dir,"relative.obj")
        write(path,PREFIX*"f -3 -2 -1\n")
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.open(path)
            _,xyz,_=A.mesh.get_nodes()
            @test Set(Tuple(p) for p in eachcol(reshape(xyz,3,:)))==Set(FACE)
            @test sum(length,A.mesh.get_elements()[2])==1
            A.merge(path)
            @test length(A.mesh.get_nodes()[1])==6
            @test sum(length,A.mesh.get_elements()[2])==2
            held=(name=A.model.get_current(),names=A.model.list(),
                nodes=A.mesh.get_nodes(),elements=A.mesh.get_elements())
            bad=joinpath(dir,"bad.obj")
            for face in ("f 0 2 3", "f -5 -2 -1", "f 2 3 5")
                write(bad,PREFIX*face*"\n")
                @test_throws ArgumentError A.open(bad)
                @test A.model.get_current()==held.name && A.model.list()==held.names
                @test A.mesh.get_nodes()==held.nodes && A.mesh.get_elements()==held.elements
            end
        finally
            A.finalize()
        end
    end
end
end
