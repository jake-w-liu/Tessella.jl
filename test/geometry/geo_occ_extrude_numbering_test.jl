module GeoOccExtrudeNumberingTest

using Test, Tessella

# Under `SetFactory("OpenCASCADE")` the `.geo` `Extrude` statement runs the
# OCC sweep (upstream `OCC_Internals::extrude`/`revolve`): OCC numbering for
# the created entities and no GEO axis-centre points. Gmsh 4.15.2 numbers the
# surfaces 1..6 and keeps 8 points for both scripts below.
function entities_after(sweep::AbstractString)
    path = joinpath(mktempdir(), "occ.geo")
    write(path, """
        SetFactory("OpenCASCADE");
        Rectangle(1)={1,0,0,1,0.5};
        Rotate {{1,0,0},{0,0,0},Pi/2} { Surface{1}; }
        Transfinite Curve{:}=3; Transfinite Surface{1}; Recombine Surface{1};
        out[] = $sweep;
        """)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(3)
        types, tags, _ = A.mesh.get_elements(3)
        return ([t for (_, t) in A.model.get_entities(0)],
                [t for (_, t) in A.model.get_entities(2)],
                Dict(t => length(e) for (t, e) in zip(types, tags)))
    finally
        A.finalize()
    end
end

@testset ".geo OCC-factory extrusions use OCC numbering" begin
    points, surfaces, cells = entities_after(
        "Extrude {0,0.7,0} { Surface{1}; Layers{2}; Recombine; }")
    @test points == collect(1:8)
    @test surfaces == collect(1:6)
    @test cells == Dict(5 => 8)
    points, surfaces, cells = entities_after(
        "Extrude {{0,0,1},{0,0,0},Pi/2} { Surface{1}; Layers{3}; Recombine; }")
    @test points == collect(1:8)
    @test surfaces == collect(1:6)
    @test cells == Dict(5 => 12)
end

end # module
