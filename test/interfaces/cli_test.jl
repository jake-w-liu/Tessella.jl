using Test
using Tessella
using Tessella.CLI: main
using Tessella.IO: read_msh
using Tessella.MeshTypes: mesh_crc, ntris, validate

const _SQUARE_GEO="""
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
"""

@testset "bounded non-destructive CLI" begin
    mktempdir() do directory
        input=joinpath(directory,"square")
        write(input,_SQUARE_GEO)

        @test main([input])==input
        @test read(input,String)==_SQUARE_GEO

        output=main([input,"-2"])
        @test output==input*".msh"
        @test isfile(output)
        @test read(input,String)==_SQUARE_GEO
        mesh=read_msh(output).mesh
        @test validate(mesh).ok && ntris(mesh)>0
        @test mesh_crc(mesh).sha==
              "92e578bac6d8feb3f0f845f100665dcc145edf924965ad72be716da933f34461"

        explicit=joinpath(directory,"explicit.msh")
        wrapped="x"*input*"x"
        input_view=SubString(wrapped,nextind(wrapped,firstindex(wrapped)),
                            prevind(wrapped,lastindex(wrapped)))
        @test main([input_view])==input
        @test main([input,"-o",explicit,"-2"])==explicit
        @test mesh_crc(read_msh(explicit).mesh)==mesh_crc(mesh)

        @test_throws ArgumentError main([input,"-2","-o",input])
        @test read(input,String)==_SQUARE_GEO
        link=joinpath(directory,"input-link")
        hardlink(input,link)
        @test_throws ArgumentError main([input,"-2","-o",link])
        @test read(input,String)==_SQUARE_GEO

        @test_throws ArgumentError main(String[])
        @test_throws ArgumentError main([input,input])
        @test_throws ArgumentError main([input,"--bad"])
        @test_throws ArgumentError main([input,"-2","-2"])
        @test_throws ArgumentError main([input,"-2","-3"])
        @test_throws ArgumentError main([input,"-2","-o","a","-o","b"])
        @test_throws ArgumentError main([input,"-2","-o"])
        @test_throws ArgumentError main([input,"-2","-o","-3"])
        @test_throws ArgumentError main([input,"-o","unused.msh"])
        @test_throws ArgumentError main([""])
        @test_throws ArgumentError main([input*"\0bad"])
        @test_throws ArgumentError main(fill("x",10_001))
        @test_throws ArgumentError main([repeat("x",1_000_001)])
    end

    @test isempty(Docs.undocumented_names(Tessella.CLI;private=false))
end
