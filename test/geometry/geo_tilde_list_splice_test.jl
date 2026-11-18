# `l~{i}()` / `l~{i}[]` splice the list the `~{}`-suffixed name resolves to,
# as in examples/boolean/twist.geo (`Curve{l~{i-1}()}`). Gmsh 4.15.2 prints
# `n=4 3`, `m=5 7` for the equivalent Printf script.
using Test
using Tessella

@testset "geo tilde list splices" begin
    dir=mktempdir()
    path=joinpath(dir,"tilde.geo")
    write(path,"""
        l~{1}() = {1:4};
        i = 1;
        l~{i+1}() = {l~{i}(), 7};
        m[] = {l~{2}[]};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test ex.lists["l_1"]==[1.0,2.0,3.0,4.0]
    @test ex.lists["l_2"]==[1.0,2.0,3.0,4.0,7.0]
    @test ex.lists["m"]==[1.0,2.0,3.0,4.0,7.0]
end
