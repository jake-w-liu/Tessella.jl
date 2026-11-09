module SurfaceRecordChecks
using Test,Tessella
const F=Tessella.IO
const A=Tessella.API
const VERTICES="0 0 0\n1 0 0\n0 1 0\n"
const PLY_HEADER="ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n"

function rejected(path)
    try
        F._read_foreign_records(path,"surface record regression")
        return false
    catch err
        err isa InterruptException && rethrow()
        return err isa ArgumentError
    end
end

@testset "OFF separated headers, comments, face counts and indices" begin
    mktempdir() do dir
        path=joinpath(dir,"records.off")
        for header in ("OFF 3 1 0\n","OFF\n3 1 0\n",
                       "# preface\nOFF # magic\n\n# counts\n3 1 0\n")
            write(path,header*VERTICES*"3 0 1 2 # triangle\n")
            labels,flat,records=F._read_off(path,"OFF record regression")
            @test labels==[1,2,3]
            @test reshape(flat,3,:)==[0. 1. 0.;0. 0. 1.;0. 0. 0.]
            @test length(records)==1 && records[1].nodes==[1,2,3]
        end
        # Optional colors after a supported face remain accepted.
        write(path,"OFF 3 1 0\n"*VERTICES*"3 0 1 2 0.5 0.5 0.5\n")
        @test length(F._read_off(path,"OFF record regression")[3])==1
        for body in ("OFF 3 1 0\n"*VERTICES,
                     "OFF 3 0 0\n"*VERTICES*"3 0 1 2\n",
                     "OFF 3 1 0\n"*VERTICES*"3 0 1\n",
                     "OFF 3 1 0\n"*VERTICES*"3 -1 1 2\n",
                     "OFF 3 1 0\n"*VERTICES*"3 0 q 2\n",
                     "OFF 3 1 0\n"*VERTICES*"3 0 1 3\n",
                     "OFF 3 1 0\n"*VERTICES*"3 0 1 9223372036854775807\n",
                     "OFF 3 1 0\n"*VERTICES*"9223372036854775807 0 1 2\n",
                     "OFF\n","OFF\nq 1 0\n","OFF -1 0 0\n")
            write(path,body)
            @test rejected(path)
        end
        ply=joinpath(dir,"records.ply2")
        for face in ("3 -1 1 2","3 0 1 3","3 0 q 2","3 0 1",
                     "3 0 1 9223372036854775807")
            write(ply,"3\n1\n"*VERTICES*face*"\n")
            @test rejected(ply)
        end
        write(ply,"# vertices\n3\n# faces\n1\n"*VERTICES*"3 0 1 2 # triangle\n")
        @test length(F._read_ply2(ply,"PLY2 record regression")[3])==1
    end
end

@testset "truncated declared vertex counts have bounded allocations" begin
    mktempdir() do dir
        for (ext,makebody) in (
            ("off",n->"OFF $n 0 0\n"),
            ("ply2",n->"$n\n0\n"),
            ("ply",n->replace(PLY_HEADER,"vertex 3"=>"vertex $n")),
            ("vtk",n->"# vtk DataFile Version 2.0\nresource\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS $n float\n"))
            path=joinpath(dir,"bounded."*ext)
            for count in (100_000,1_000_000_000,typemax(Int))
                write(path,makebody(count))
                @test rejected(path)
                rejected(path) # warm the exact malformed path
                @test (@allocated rejected(path))<100_000
            end
        end
        # Binary PLY must use the actual scalar widths for its feasibility
        # bound, including properties beyond the three coordinate fields.
        path=joinpath(dir,"binary.ply")
        for encoding in ("binary_little_endian","binary_big_endian")
            write(path,replace(PLY_HEADER,"ascii"=>encoding,
                               "vertex 3"=>"vertex 100000"))
            @test rejected(path)
            rejected(path)
            @test (@allocated rejected(path))<100_000
        end
        # VTK's token iterator may buffer coordinates on the POINTS line.
        vtk=joinpath(dir,"inline.vtk")
        write(vtk,"# vtk DataFile Version 2.0\ninline\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 1 float 0 0 0\nCELLS 1 2\n1 0\nCELL_TYPES 1\n1\n")
        @test length(F._read_foreign_records(vtk,"VTK inline regression")[1])==1
    end
end

@testset "OFF public import handles separated headers and rolls back failures" begin
    mktempdir() do dir
        valid=joinpath(dir,"valid.off")
        write(valid,"OFF\n# counts\n3 1 0\n"*VERTICES*"3 0 1 2\n")
        A.initialize()
        try
            A.option("General.Terminal",0)
            @test A.open(valid)===nothing
            @test length(A.mesh.get_nodes()[1])==3
            @test sum(length,A.mesh.get_elements()[2])==1
            held=(name=A.model.get_current(),names=A.model.list(),
                  nodes=A.mesh.get_nodes(),elements=A.mesh.get_elements())
            bad=joinpath(dir,"bad.off")
            for body in ("OFF 3 1 0\n"*VERTICES,
                         "OFF 3 1 0\n"*VERTICES*"3 -1 1 2\n",
                         "OFF 1000000000 0 0\n")
                write(bad,body)
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
