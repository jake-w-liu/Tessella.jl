# `plugin("Bubbles")` — Gmsh 4.15.2 `Bubbles.cpp`: writes one `.geo` file
# containing the model's points/lines/surfaces followed by the "bubbles" of
# each face's triangulation vertex-Voronoi (cell-boundary points =
# incident-triangle barycenters plus the vertex's own position for nodes on
# dim<2 entities, angle-sorted about their barycenter and shrunk by
# `ShrinkFactor`). Each input face contributes one Physical Surface of
# bubble plane-surfaces plus a leftover `Plane Surface` holding the face
# loop and every bubble loop. On a unit square the file must be valid .geo
# and round-trip through `open_geo`.
module PluginBubblesChecks
using Test,Tessella
const A=Tessella.API

function _square()
    A.model.add("sq")
    for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
        A.model.occ.add_point(x,y,0.)
    end
    for i in 1:4; A.model.occ.add_line(i,mod(i,4)+1); end
    A.model.occ.add_curve_loop([1,2,3,4])
    A.model.occ.add_plane_surface([1])
    A.model.occ.synchronize()
    A.mesh.generate(2)
end

@testset "Plugin(Bubbles) writes a round-trippable .geo" begin
    A.initialize()
    try
        _square()
        path=joinpath(mktempdir(),"bubbles.geo")
        A.plugin.set_string("Bubbles","OutputFile",path)
        A.plugin.run("Bubbles")
        @test isfile(path)
        text=read(path,String)
        @test occursin(r"^lc = [\d.]+;"m,text)
        # The original square's entities are emitted verbatim.
        @test occursin("Point(1) = {0, 0, 0, lc};",text)
        @test occursin("Line(4) = {4, 1};",text)
        @test occursin("Curve Loop(1) = {1, 2, 3, 4};",text)
        @test occursin("Plane Surface(1) = {1};",text)
        @test occursin("Delete { Surface {1}; }",text)
        # Every BSpline closes on its first point (a bubble ring).
        for m in eachmatch(r"BSpline\(\d+\) = \{([\d, ]+)\};",text)
            ids=parse.(Int,split(strip(m.captures[1]),r",\s*"))
            @test first(ids)==last(ids)
            @test length(ids)>=4
        end
        @test !isempty(collect(eachmatch(r"BSpline\(\d+\)",text)))
        # The bubbles physical group plus the leftover face group.
        phys=[m.captures[1] for m in
            eachmatch(r"Physical Surface\(\d+\) = \{([^}]*)\};",text)]
        @test length(phys)==2
        @test occursin(":",phys[1])  # contiguous bubble surface range
        # The leftover surface carries the face loop + every bubble loop.
        leftover=match(r"Plane Surface\(\d+\) = \{1, (\d+):(\d+)\};",text)
        @test leftover!==nothing
        @test parse(Int,leftover.captures[1])<=parse(Int,leftover.captures[2])
        # The emitted file is well-formed .geo: it re-imports cleanly.
        A.open_geo!(path)
        @test !isempty(A.model.get_entities(2))
    finally
        A.finalize()
    end
end

@testset "Plugin(Bubbles) ShrinkFactor shrinks cells" begin
    A.initialize()
    try
        _square()
        dir=mktempdir()
        p0=joinpath(dir,"b0.geo");p1=joinpath(dir,"b1.geo")
        A.plugin.set_string("Bubbles","OutputFile",p0)
        A.plugin.run("Bubbles")
        A.plugin.set_number("Bubbles","ShrinkFactor",0.5)
        A.plugin.set_string("Bubbles","OutputFile",p1)
        A.plugin.run("Bubbles")
        @test isfile(p0) && isfile(p1)
        @test read(p0,String)!=read(p1,String)
    finally
        A.finalize()
    end
end
end
