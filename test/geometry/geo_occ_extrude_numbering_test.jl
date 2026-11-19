module GeoOccExtrudeNumberingTest

using Test, Tessella
const A = Tessella.API

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

# The sweep's new entities take each dimension's next free tag in the
# result's exploration order, so every tag has to sit where Gmsh puts it: the
# built-in sweep numbered a rectangle's prism 13, 17, 21, 25 and 26, and a
# script that addresses the top as `Surface{6}` hit nothing. Tags, the
# returned list and each entity's position are Gmsh 4.15.2's.
function sweep(text)
    dir=mktempdir()
    path=joinpath(dir,"sweep.geo")
    write(path,"SetFactory(\"OpenCASCADE\");\n"*text)
    ex=execute_geo(path)
    @test isempty(ex.errors)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        centre(d,t)=(b=A.model.get_bounding_box(d,t);
            round.(((b[1]+b[4])/2,(b[2]+b[5])/2,(b[3]+b[6])/2);digits=6).+0.0)
        return Int.(ex.lists["out"]),
            Dict((d,t)=>centre(d,t) for d in 1:2 for (_,t) in A.model.get_entities(d)),
            ex.values
    finally
        A.finalize()
    end
end

@testset "geo OCC Extrude numbers its entities like OCC" begin
    out,at,values=sweep("""
        Rectangle(1) = {0,0,0,1,1};
        out[] = Extrude {0,0,1} { Surface{1}; Layers{2}; Recombine; };
        """)
    @test out==[6,1,2,3,4,5]
    @test at[(2,6)]==(0.5,0.5,1.0)
    @test [at[(2,t)] for t in 2:5]==
        [(0.5,0.0,0.5),(1.0,0.5,0.5),(0.5,1.0,0.5),(0.0,0.5,0.5)]
    # per wall: the two swept-vertex edges, then the translated copy
    @test [at[(1,t)] for t in 5:12]==[(0.0,0.0,0.5),(1.0,0.0,0.5),
        (0.5,0.0,1.0),(1.0,1.0,0.5),(1.0,0.5,1.0),(0.0,1.0,0.5),
        (0.5,1.0,1.0),(0.0,0.5,1.0)]

    # a quarter turn: arcs before the rotated copy, the walls before the cap
    out,at,values=sweep("""
        Rectangle(1) = {1,0,0,1,1};
        out[] = Extrude {{0,1,0},{0,0,0},Pi/2} { Surface{1}; };
        """)
    @test out==[6,1,2,3,4,5]
    @test at[(2,6)]==(0.0,0.5,-1.5)
    @test [at[(1,t)] for t in 5:9]==[(0.5,0.0,-0.5),(1.0,0.0,-1.0),
        (0.0,0.0,-1.5),(1.0,1.0,-1.0),(0.0,0.5,-2.0)]

    # a prism along a direction that is not the profile's normal used to be
    # rejected by the OCC sweep ("lateral face winding is inconsistent")
    out,at,values=sweep("""
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={0,1,0};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};
        Curve Loop(1)={1,2,3}; Plane Surface(1)={1};
        out[] = Extrude {0.2,0.1,1} { Surface{1}; };
        v = Mass Volume{1};
        """)
    @test out==[5,1,2,3,4]
    @test at[(2,5)]==(0.7,0.6,1.0)
    @test values["v"]≈0.5 rtol=1e-12
end

end # module
