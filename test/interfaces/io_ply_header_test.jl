module PlyHeaderChecks
using Test,Tessella
const F=Tessella.IO
const A=Tessella.API
const VERTICES="0 0 0\n1 0 0\n0 1 0\n"

function header(encoding,nv,nf)
    faces=nf===nothing ? "" : "element face $nf\nproperty list uchar int vertex_indices\n"
    return "ply\nformat $encoding 1.0\nelement vertex $nv\nproperty float x\nproperty float y\nproperty float z\n"*faces*"end_header\n"
end

function write_count_file(path,encoding,nv,nf)
    open(path,"w") do io
        write(io,header(encoding,nv,nf))
        if encoding=="ascii"
            write(io,VERTICES)
        else
            pack=encoding=="binary_little_endian" ? htol : hton
            for value in (0.,0.,0.,1.,0.,0.,0.,1.,0.)
                write(io,pack(reinterpret(UInt32,Float32(value))))
            end
        end
    end
end

@testset "PLY declared counts distinguish missing, zero and negative values" begin
    mktempdir() do dir
        path=joinpath(dir,"counts.ply")
        for encoding in ("ascii","binary_little_endian","binary_big_endian")
            for count in (-1,-2,typemin(Int))
                write_count_file(path,encoding,3,count)
                @test_throws ArgumentError F._read_foreign_records(path,"PLY count regression")
                write_count_file(path,encoding,count,0)
                @test_throws ArgumentError F._read_foreign_records(path,"PLY count regression")
            end
        end
        for count in (nothing,0)
            write(path,header("ascii",3,count)*VERTICES)
            @test isempty(F._read_foreign_records(path,"PLY count regression")[3])
        end
        write(path,header("ascii",3,1)*VERTICES*"3 0 1 2\n")
        @test only(F._read_foreign_records(path,"PLY count regression")[3]).nodes==[1,2,3]
        # A duplicate declaration remains invalid even when its count is zero.
        write(path,replace(header("ascii",3,0),"end_header"=>"element face 0\nend_header"))
        @test_throws ArgumentError F._read_foreign_records(path,"PLY count regression")
    end
end

@testset "PLY invalid counts preserve the current model on open and merge" begin
    mktempdir() do dir
        path=joinpath(dir,"valid.ply")
        write(path,header("ascii",3,1)*VERTICES*"3 0 1 2\n")
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.open(path)
            held=(name=A.model.get_current(),names=A.model.list(),
                  nodes=A.mesh.get_nodes(),elements=A.mesh.get_elements())
            bad=joinpath(dir,"negative.ply")
            for count in (-1,-2),operation in (A.open,A.merge)
                write(bad,header("ascii",3,count)*VERTICES)
                @test_throws ArgumentError operation(bad)
                @test A.model.get_current()==held.name && A.model.list()==held.names
                @test A.mesh.get_nodes()==held.nodes && A.mesh.get_elements()==held.elements
            end
        finally
            A.finalize()
        end
    end
end
end
