module ForeignResourceChecks
using Test,Tessella
const F=Tessella.IO

function rejected(path)
    try
        F._read_foreign_records(path,"native reader resource regression")
        return false
    catch err
        err isa InterruptException && rethrow()
        return err isa ArgumentError
    end
end

function diffhead(nv,ne)
    return "Number of space dim. = 3\nNumber of elements   = $ne\nNumber of nodes      = $nv\nMax number of nodes in an element: 3\nOnly one material/subdomain\n0 Boundary indicators:\n#\n"
end

@testset "native reader counts cannot reserve impossible payloads" begin
    mktempdir() do dir
        for (ext,body) in (
            ("p3d",n->"1\n$n 1 1\n"),
            ("mail",n->"$n 0\n"),("mail",n->"0 $n\n"),
            ("ir3",n->"33\n$n 0\n"),("ir3",n->"33\n0 $n\n"),
            ("unv",n->"    -1\n  2412\n1 91 0 0 7 $n\n"),
            ("diff",n->diffhead(n,0)),("diff",n->diffhead(0,n)*"#\n"))
            path=joinpath(dir,"bounded."*ext)
            for count in (100_000,1_000_000_000,typemax(Int))
                write(path,body(count))
                @test rejected(path)
                rejected(path)
                @test (@allocated rejected(path))<100_000
            end
        end
        # Unsupported UNV descriptors still consume their label payloads;
        # their untrusted widths need the same bound before sizehint!.
        path=joinpath(dir,"unknown.unv")
        write(path,"    -1\n  2412\n1 999 0 0 7 100000\n")
        @test rejected(path)
        rejected(path)
        @test (@allocated rejected(path))<100_000

        # Header arithmetic is bounded before multiplying dimensions or
        # summing element sections, including values at the integer limit.
        path=joinpath(dir,"overflow.p3d")
        write(path,"1\n$(typemax(Int)) $(typemax(Int)) 2\n")
        @test_throws ArgumentError F.read_p3d(path)
        path=joinpath(dir,"overflow.ir3")
        write(path,"33\n0 $(typemax(Int)) $(typemax(Int))\n")
        @test_throws ArgumentError F.read_ir3(path)
        path=joinpath(dir,"indicators.diff")
        write(path,replace(diffhead(0,0),"0 Boundary"=>"$(typemax(Int)) Boundary"))
        @test_throws ArgumentError F.read_diff(path)
        @test (@allocated rejected(path))<100_000
    end
end

@testset "bounded native readers retain compact positive payloads" begin
    mktempdir() do dir
        path=joinpath(dir,"inline.p3d")
        write(path,"1 2 2 1 0 1 0 1 0 0 1 1 0 0 0 0\n")
        mesh=F.read_p3d(path).mesh
        @test size(mesh.coords,2)==4 && mesh.blocks[1].msh==3
        path=joinpath(dir,"inline.mail")
        write(path,"3 1 0 0 0 1 0 0 0 1 0 1 2 3\n")
        mesh=F.read_mail(path).mesh
        @test mesh.blocks[1].nodes==reshape(Int32[1,2,3],3,1)
        path=joinpath(dir,"inline.ir3")
        write(path,"33\n3 1\n1 0 0 0 2 1 0 0 3 0 1 0 1 7 3 1 2 3\n")
        mesh=F.read_ir3(path).mesh
        @test mesh.blocks[1].nodes==reshape(Int32[1,2,3],3,1)
        @test mesh.elementary_entities==[Int32[7]]
        for (ext,reader,body) in (("p3d",F.read_p3d,"1\n100000 1 1\n"),
                                  ("mail",F.read_mail,"100000 0\n"),
                                  ("ir3",F.read_ir3,"33\n100000 0\n"),
                                  ("unv",F.read_unv,"    -1\n  2412\n1 91 0 0 7 100000\n"),
                                  ("diff",F.read_diff,diffhead(100000,0)))
            path=joinpath(dir,"public."*ext)
            write(path,body)
            @test_throws ArgumentError reader(path)
        end
    end
end
end
