module GeoOCCThruSectionsShellChecks
using Test,Tessella
const A=Tessella.API

# The bare `ThruSections{wires}` term lofts a shell — Gmsh.y calls
# `addThruSections(-1, wires, false, …)` — and only the tagged definition
# `ThruSections(tag) = {wires};` makes a solid. Tessella asked for a solid in
# both, so open section wires failed with "Making solid requires closed
# wires" (thrusections.geo lines 49 and 51, hyperboloid.geo). Gmsh 4.15.2
# has 4 points, 4 curves and 1 surface for the two arcs below.
@testset "geo bare ThruSections lofts a shell" begin
    dir=mktempdir()
    path=joinpath(dir,"loft.geo")
    write(path,"""
        SetFactory("OpenCASCADE");
        Circle(1) = {0,0,0, 0.2, Pi/3};
        Circle(2) = {0.1,0.04,1, 0.3, Pi/2};
        Wire(1) = 1; Wire(2) = 2;
        out[] = ThruSections{1,2};
        """)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    @test length(ex.lists["out"])==1
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(4,4,1,0)
    finally
        A.finalize()
    end
end
end
