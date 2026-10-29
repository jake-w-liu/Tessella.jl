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
        out=joinpath(dir,"named.geo_unrolled")
        dest=main([geo,"-0","-o",out])
        @test dest==out
        @test isfile(out)
        # Upstream `-0 -o x.geo` answers `Unknown output file format` — the
        # name's extension selects the payload format, validated at
        # `CreateOutputFile` after the script runs.
        err=try
            main([geo,"-0","-o",joinpath(dir,"named.geo")]);""
        catch e;sprint(showerror,e);end
        @test occursin("unknown output file format",err)
        # `-2 -o` to a `.geo` name rejects the same way (after generation).
        err=try
            main([geo,"-2","-o",joinpath(dir,"named.geo")]);""
        catch e;sprint(showerror,e);end
        @test occursin("unknown output file format",err)
        # Output aliasing the input path still rejects.
        err=try
            main([geo,"-2","-o",geo]);""
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

@testset "CLI -convert: upstream _new payloads" begin
    mktempdir() do dir
        # Geometry-only `.geo`: upstream `OpenProject` accumulates no views
        # and no mesh, so `-convert` creates no `_new` at all.
        geo=joinpath(dir,"m.geo");write(geo,_UNROLL_GEO)
        dest=main(["-convert",geo])
        @test dest==geo*"_new"
        @test !isfile(dest)
        # A `.geo` whose script generates a mesh emits binary MSH 4.1.
        meshed=joinpath(dir,"meshed.geo")
        write(meshed,"""
Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
Mesh 2;
""")
        dest=main(["-convert",meshed])
        @test dest==meshed*"_new"
        @test isfile(dest)
        bytes=read(dest)
        # Binary MSH 4.1 header: `$MeshFormat`, `4.1 1 8`, then the int32
        # endianness sentinel.
        @test startswith(String(bytes[1:30]),"\$MeshFormat\n4.1 1 8\n")
        @test reinterpret(Int32,bytes[21:24])==Int32[1]
        mesh=read_mixed_msh(dest)
        @test mesh.coords isa Matrix{Float64}
        # A `.pos` input converts to the binary `.pos` 1.3 layout upstream
        # writes for `-convert` (`PView::write(name, 1, append)`).
        pos=joinpath(dir,"v.pos")
        write(pos,"""View "sample" {
SP(0,0,0){1.5};
SL(0,0,0,1,0,0){2,3};
ST(0,0,0,1,0,0,0,1,0){4,5,6};
T3(0,0,0,17){"note"};
};
""")
        dest=main(["-convert",pos])
        @test dest==pos*"_new"
        raw=read(dest)
        # `String(::Vector{UInt8})` takes ownership of the buffer and leaves
        # the vector empty — copy before each conversion.
        text=String(copy(raw))
        @test startswith(text,"\$PostFormat /* Gmsh 1.3, binary */")
        @test occursin("\$View /* sample */",text)
        @test occursin("\$EndView",text)
        # Header: name + 1 step + 24 family counts (SP=SL=ST=1 in the
        # scalar slot of their triples) + T2/T3 record/byte counts
        # (`0 0 1 5` — one T3 text of 5 chars), verified byte-for-byte
        # against the 4.15.2 oracle's `v.pos_new`.
        counts=vcat([1,0,0, 1,0,0, 1,0,0],zeros(Int,15),[0,0,1,5])
        hdr="sample "*join(vcat([1],counts)," ")*"\n"
        @test occursin(hdr,text)
        # A mesh-bearing foreign input rewrites as binary MSH 4.1.
        stl=joinpath(dir,"t.stl")
        write(stl,"""solid t
facet normal 0 0 1
outer loop
vertex 0 0 0
vertex 1 0 0
vertex 0 1 0
endloop
endfacet
endsolid t
""")
        dest=main(["-convert",stl])
        @test dest==stl*"_new"
        @test read_mixed_msh(dest).coords isa Matrix{Float64}
        # Unsupported extensions keep the named rejection.
        err=try
            main(["-convert",joinpath(dir,"x.step")]);""
        catch e;sprint(showerror,e);end
        @test occursin("does not exist",err)
    end
end

@testset "CLI -format/-o: foreign output payloads" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo")
        write(geo,"""
Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
""")
        # `-o <name>.stl` on a `.geo` input writes real STL bytes — not an
        # MSH payload under a foreign name.
        stl=joinpath(dir,"sq.stl")
        dest=main([geo,"-2","-o",stl])
        @test dest==stl
        @test startswith(read(stl,String),"solid")
        # `-format` picks the payload format upstream-style: `-o` naming
        # does not override it.
        odd=joinpath(dir,"odd.msh")
        dest=main([geo,"-2","-format","vtk","-o",odd])
        @test dest==odd
        @test startswith(read(odd,String),"# vtk DataFile")
        # No `-o`: `-format` defaults the destination to the canonical
        # extension (`GetDefaultFileName(Mesh.FileFormat)`).
        main([geo,"-2","-format","off"])
        @test isfile(joinpath(dir,"sq.off"))
        @test occursin("OFF",read(joinpath(dir,"sq.off"),String))
        # `-format geo_unrolled` and `-o <name>.geo_unrolled` emit the
        # unrolled model script (`FORMAT_GEO` upstream).
        unrolled=joinpath(dir,"sq.geo_unrolled")
        dest=main([geo,"-2","-format","geo_unrolled","-o",unrolled])
        @test dest==unrolled
        @test occursin("Plane Surface(1) = {1};",read(unrolled,String))
        # Unknown format names and unserializable `-o` extensions reject
        # like upstream's `Unknown output file format`.
        err=try
            main([geo,"-2","-format","bogus"]); ""
        catch e;sprint(showerror,e);end
        @test occursin("unknown output format",err)
        err=try
            main([geo,"-2","-o",joinpath(dir,"x.pos")]); ""
        catch e;sprint(showerror,e);end
        @test occursin("POS output is not implemented",err)
        err=try
            main([geo,"-2","-o",joinpath(dir,"x.bogus")]); ""
        catch e;sprint(showerror,e);end
        @test occursin("unknown output file format",err)
        err=try
            main([geo,"-2","-o",joinpath(dir,"noext")]); ""
        catch e;sprint(showerror,e);end
        @test occursin("unknown output file format",err)
    end
end

@testset "CLI flags: oracle-verified tail" begin
    mktempdir() do dir
        geo=joinpath(dir,"sq.geo")
        write(geo,"Point(1) = {0,0,0,0.5};\n")
        # `-noenv`/`-nolocale` are startup environment controls upstream —
        # batch no-ops here, accepted silently.
        @test main([geo,"-noenv","-nolocale","-0"])==
            joinpath(dir,"sq.geo_unrolled")
        # `-hybrid` sets Recombine3DAll + ALGO_3D_RTREE + ALGO_2D_PACK_PRLGRMS
        # upstream; none has a native producer, so the named rejection stays.
        for flag in ("-hybrid","-dual","-voronoi","-part_weight")
            err=try
                main([geo,flag]); ""
            catch e;sprint(showerror,e);end
            @test occursin("not supported",err)
        end
        err=try
            main([geo,"-partWeight","tet","2"]); ""
        catch e;sprint(showerror,e);end
        @test occursin("weighted graph partition",err)
        # Every upstream `-algo` spelling parses and stores the real option
        # value; unknown names keep the named rejection.
        for name in ("bamg","front3d","delquad","quadqs","auto","hxt")
            err=try
                main([geo,"-algo",name]); "ok"
            catch e;sprint(showerror,e);end
            @test !occursin("unknown mesh algorithm",err)
        end
        err=try
            main([geo,"-algo","nosuchalg"]); ""
        catch e;sprint(showerror,e);end
        @test occursin("unknown mesh algorithm",err)
        # `-optimize` keeps working but warns like upstream's obsolete flag.
        # `redirect_stderr` needs an OS-level stream — a `Pipe`, not an
        # `IOBuffer`.
        pipe=Pipe()
        redirect_stderr(pipe) do
            main([geo,"-optimize"])
        end
        close(pipe.in)
        @test occursin("obsolete",String(read(pipe.out)))
        # `-convert` with no following file exits silently like upstream.
        @test main(["-convert"])==""
    end
end
end
