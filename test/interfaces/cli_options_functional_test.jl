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

const SQUARE="""
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
# Scalar triangle coverage of [-1,2]^2 at h=0.1 — the PostView background
# field must densify the unit square well past the default h=1 grading.
const BGM_POS="""
View "h" {
ST(-1,-1,0, 2,-1,0, -1,2,0){0.1,0.1,0.1};
ST(2,-1,0, 2,2,0, -1,2,0){0.1,0.1,0.1};
};
"""

@testset "CLI option flags wired to stores" begin
    mktempdir() do directory
        source=joinpath(directory,"read.geo")
        output=joinpath(directory,"read.msh")
        # Each flag's value must land in the .geo option store; `Printf`
        # reads observe the forced value with priority over assignments.
        function readback(option,flags)
            target=replace(joinpath(directory,"value.txt"),'\\'=>'/')
            write(source,LINE*"Printf(\"%.17g\",$option) > \"$target\";")
            Tessella.CLI.main(vcat([source,"-1","-o",output],flags))
            return parse(Float64,read(joinpath(directory,"value.txt"),String))
        end
        @test readback("Mesh.AnisoMax",["-aniso_max","7"])==7.0
        @test readback("Mesh.AnisoMax",["-anisoMax","9"])==9.0
        @test readback("Mesh.RandomFactor",["-rand","0.5"])==0.5
        @test readback("Mesh.SmoothRatio",["-smooth_ratio","2.2"])==2.2
        @test readback("Mesh.SmoothRatio",["-smoothRatio","2.5"])==2.5
        @test readback("Mesh.IgnorePeriodicity",["-ignore_periodicity"])==1.0
        @test readback("Mesh.IgnorePeriodicity",["-ignorePeriodicity"])==1.0
        @test readback("Mesh.MeshSizeFromCurvatureIsotropic",
                       ["-clcurviso"])==1.0
        # Script assignments lose to the forced flag value.
        write(source,"Mesh.AnisoMax=3;"*LINE)
        Tessella.CLI.main([source,"-1","-aniso_max","8","-o",output])
        @test isfile(output)
        target=replace(joinpath(directory,"value.txt"),'\\'=>'/')
        write(source,"Mesh.AnisoMax=3;"*LINE*"Printf(\"%.17g\","*
                     "Mesh.AnisoMax) > \"$target\";")
        Tessella.CLI.main([source,"-1","-aniso_max","8","-o",output])
        @test parse(Float64,read(joinpath(directory,"value.txt"),String))==8.0
        # Both spellings on one option reject like other option aliases.
        @test_throws ArgumentError Tessella.CLI.main(
            [source,"-1","-aniso_max","5","-anisoMax","6"])
        # Flags that never existed upstream reject as obsolete spellings.
        for flag in ("-ignore_periocity","-claniso")
            @test occursin("not an upstream",
                try Tessella.CLI.main([source,"-1",flag]);"" catch e e.msg end)
        end
    end
end

@testset "CLI -bgm background mesh" begin
    mktempdir() do directory
        source=joinpath(directory,"sq.geo");write(source,SQUARE)
        pos=joinpath(directory,"bg.pos");write(pos,BGM_POS)
        output=joinpath(directory,"sq.msh")
        Tessella.CLI.main([source,"-2","-o",output])
        baseline=Tessella.Elements.read_mixed_msh(output)
        @test size(baseline.coords,2)<20
        Tessella.CLI.main([source,"-2","-bgm",pos,"-o",output])
        refined=Tessella.Elements.read_mixed_msh(output)
        @test Tessella.validate(refined).ok
        @test size(refined.coords,2)>10*size(baseline.coords,2)
        # Session path (-part forces the API session): identical effect.
        partitioned=joinpath(directory,"sqp.msh")
        Tessella.CLI.main([source,"-2","-part","2","-bgm",pos,"-o",partitioned])
        through_session=Tessella.Elements.read_mixed_msh(partitioned)
        @test size(through_session.coords,2)>10*size(baseline.coords,2)
        # A missing file warns and meshes on — upstream's non-fatal arm.
        Tessella.CLI.main([source,"-2","-bgm",
            joinpath(directory,"absent.pos"),"-o",output])
        absent=Tessella.Elements.read_mixed_msh(output)
        @test size(absent.coords,2)==size(baseline.coords,2)
        # Non-.pos files and nonexistent parses reject loudly.
        nonpos=joinpath(directory,"bg.msh")
        Tessella.CLI.main([source,"-2","-o",nonpos])
        @test_throws ArgumentError Tessella.CLI.main(
            [source,"-2","-bgm",nonpos,"-o",output])
        @test_throws ArgumentError Tessella.CLI.main([source,"-2","-bgm"])
        @test_throws ArgumentError Tessella.CLI.main(
            [output,"-check","-bgm",pos])
    end
end
end
