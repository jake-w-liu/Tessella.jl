# CLI `-0` — verified against the Gmsh 4.15.2 oracle: the model executes
# without meshing and serializes as an unrolled `.geo` script to
# `<input>_unrolled` (or the `-o` path). Distinct point characteristic
# lengths emit as `cl__N` variables in first-use order; entity statements
# interleave Transfinite/Recombine/embeds after their owner; periodic
# constraints and non-built-in entity kinds are dropped upstream.
module CliUnrolledGeoChecks
using Test, Tessella
using Tessella.CLI: main
using Tessella.GeoUnrolled: write_geo_unrolled
using Tessella.GeoExec: execute_geo
using Tessella.Elements: read_mixed_msh, write_mixed_msh

const _UNROLL_GEO="""
Point(1) = {0, 0, 0, 0.3};
Point(2) = {1, 0, 0, 0.3};
Point(3) = {0, 1, 0, 0.5};
Point(4) = {0.5, 0.5, 0, 0.2};
Point(5) = {0.3, 0.3, 0, 0.4};
Line(1) = {1, 2};
Line(2) = {2, 3};
Circle(3) = {3, 4, 1};
Transfinite Curve {1} = 5 Using Progression 1.2;
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};
Recombine Surface {1};
Point {5} In Surface {1};
Physical Surface("dom", 7) = {1};
"""

@testset "CLI -0: unrolled .geo output" begin
    mktempdir() do dir
        geo=joinpath(dir,"m.geo");write(geo,_UNROLL_GEO)
        dest=main([geo,"-0"])
        @test dest==geo*"_unrolled"
        @test isfile(dest)
        text=read(dest,String)
        # cl__N variables in first-use order — 0.3 appears before 0.5.
        @test occursin("cl__1 = 0.3;",text)
        @test occursin("cl__2 = 0.5;",text)
        @test occursin("Point(1) = {0, 0, 0, cl__1};",text)
        @test occursin("Point(3) = {0, 1, 0, cl__2};",text)
        @test occursin("Line(1) = {1, 2};",text)
        @test occursin("Circle(3) = {3, 4, 1};",text)
        # Transfinite interleaves immediately after its curve; the surface's
        # loop emits lazily before the surface statement.
        @test occursin("Transfinite Curve {1} = 5 Using Progression 1.2;",text)
        @test occursin("Curve Loop(1) = {1, 2, 3};",text)
        @test occursin("Plane Surface(1) = {1};",text)
        @test occursin("Recombine Surface {1};",text)
        @test occursin("Point {5} In Surface {1};",text)
        @test occursin("Physical Surface(\"dom\") = {1};",text)
        # Oracle drops periodic constraints entirely.
        @test !occursin("Periodic",text)
    end
end

@testset "CLI -0: -o destination and alias rejection" begin
    mktempdir() do dir
        geo=joinpath(dir,"m.geo");write(geo,_UNROLL_GEO)
        out=joinpath(dir,"named.geo")
        dest=main([geo,"-0","-o",out])
        @test dest==out
        @test isfile(out)
        err=try
            main([geo,"-0","-o",geo]);""
        catch e;sprint(showerror,e);end
        @test occursin("alias",err)
        err=try
            main([geo,"-0","-2"]);""
        catch e;sprint(showerror,e);end
        @test occursin("conflicting batch action",err)
    end
end

@testset "CLI -0: .msh input emits discrete points" begin
    mktempdir() do dir
        geo=joinpath(dir,"line.geo")
        write(geo,"""
Point(1) = {0,0,0,0.4}; Point(2) = {1,0,0,0.4};
Line(1) = {1,2};
""")
        msh=joinpath(dir,"line.msh")
        main([geo,"-1","-o",msh])
        dest=main([msh,"-0"])
        @test dest==splitext(msh)[1]*".geo_unrolled"
        text=read(dest,String)
        # Only point records serialize from a mesh input — one cl__1 size
        # variable at the upstream unset value.
        @test occursin("cl__1 = 1e+22;",text)
        @test occursin("Point(1) = {0, 0, 0, cl__1};",text)
        @test occursin("Point(2) = {1, 0, 0, cl__1};",text)
        @test !occursin("Line(",text)
        @test !occursin("Surface",text)
    end
end

@testset "CLI -convert: .geo input unrolls" begin
    mktempdir() do dir
        geo=joinpath(dir,"m.geo");write(geo,_UNROLL_GEO)
        dest=main(["-convert",geo])
        @test dest==geo*"_new"
        text=read(dest,String)
        @test occursin("cl__1 = 0.3;",text)
        @test occursin("Plane Surface(1) = {1};",text)
    end
end
end
