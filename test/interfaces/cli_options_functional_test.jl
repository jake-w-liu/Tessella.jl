module CLIOptionsFunctionalChecks
using Test, Tessella
const E=Tessella.Elements
const LINE="Point(1)={0,0,0,0.25};Point(2)={1,0,0,0.25};Line(1)={1,2};"
curve_cells(m)=sum(size(b.nodes,2) for b in m.blocks if E.msh_dimension(Int(b.msh))==1;init=0)

@testset "CLI numeric options and information" begin
    mktempdir() do directory
        source=joinpath(directory,"line.geo");output=joinpath(directory,"line.msh")
        write(source,LINE)
        Tessella.CLI.main([source,"-1","-o",output])
        baseline=E.read_mixed_msh(output)
        Tessella.CLI.main([source,"-1","-clscale","0.5","-o",output])
        smaller=E.read_mixed_msh(output)
        @test Tessella.validate(smaller).ok
        @test curve_cells(smaller)>curve_cells(baseline)
        @test isapprox(maximum(diff(sort(smaller.coords[1,:]))),0.125;rtol=1e-10)
        Tessella.CLI.main([source,"-1","-clmin","0.5","-clmax","0.5","-o",output])
        bounded=E.read_mixed_msh(output)
        @test Tessella.validate(bounded).ok
        @test curve_cells(bounded)==2
        Tessella.CLI.main([source,"-1","-order","3","-o",output])
        cubic=E.read_mixed_msh(output)
        @test Tessella.validate(cubic).ok
        @test only(filter(b->E.msh_dimension(Int(b.msh))==1,cubic.blocks)).msh==26
        @test curve_cells(cubic)==curve_cells(baseline)
        # Alias writes and command execution in the script still observe forced values.
        write(source,LINE*"Mesh.CharacteristicLengthFactor=4;Mesh.ElementOrder=1;Mesh 1;")
        Tessella.CLI.main([source,"-1","-clscale","0.5","-order","3","-o",output])
        forced=E.read_mixed_msh(output)
        @test curve_cells(forced)==curve_cells(smaller)
        @test only(filter(b->E.msh_dimension(Int(b.msh))==1,forced.blocks)).msh==26
        saved=read(output);original=read(source)
        for flags in (["-clscale"],["-clscale","0"],["-clscale","Inf"],
                      ["-clscale","NaN"],["-clmin","-1"],["-clmax","0"],
                      ["-clmax","bad"],["-order","1.5"],["-order","11"],
                      ["-order","0"],["-clmin","2","-clmax","1"],
                      ["-clscale","1","-clscale","2"],["-help"],["-version"])
            @test_throws ArgumentError Tessella.CLI.main(vcat([source,"-1","-o",output],flags))
        end
        @test read(output)==saved && read(source)==original
        help=redirect_stdout(devnull) do
            Tessella.CLI.main(["--help"])
        end
        version=redirect_stdout(devnull) do
            Tessella.CLI.main(["-version"])
        end
        @test occursin("-clscale",help) && occursin("-order",help)
        @test version==string(Base.pkgversion(Tessella))
    end
end
end
