module GmshIOCompatibilityChecks
using Test
using Tessella.MeshTypes
using Tessella.IO

function _cube()
    coords = Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    raw = [ (1,2,3,7), (1,2,6,7), (1,4,3,7), (1,4,8,7), (1,5,6,7), (1,5,8,7) ]
    tets = Matrix{Int32}(undef, 4, length(raw))
    tags = Int32[]
    for (k, t) in enumerate(raw)
        f(i) = (coords[1,i],coords[2,i],coords[3,i])
        v = tet_signed_volume(f(t[1]),f(t[2]),f(t[3]),f(t[4]))
        tt = v < 0 ? (t[1],t[2],t[4],t[3]) : t
        tets[:,k] = Int32[tt...]
        push!(tags, Int32(4))     # physical tag 4 ("case")
    end
    return Mesh(coords; tets=tets, tet_tag=tags)
end

@testset "Live Gmsh I/O compatibility" begin
    mktempdir() do dir
        @testset "Gmsh v4 rewrite preserves nodes, tags, and CRC" begin
            gmsh=Sys.which("gmsh")
            @test gmsh!==nothing
            if gmsh!==nothing
                version=strip(read(`$gmsh --version`,String))
                @test version=="4.15.2" || startswith(version,"4.15.2-")
                m=_cube()
                m.tet_tag[1]=Int32(5)
                names=Dict((3,4)=>"case\tregion",(3,5)=>"insert")
                source=joinpath(dir,"gmsh_rewrite_source.msh")
                rewritten=joinpath(dir,"gmsh_rewrite_result.msh")
                write_msh(source,m;version=4.1,physical_names=names)
                output=IOBuffer()
                process=run(pipeline(ignorestatus(
                    `$gmsh $source -0 -save -format msh41 -o $rewritten -v 5`);
                    stdout=output,stderr=output))
                log=String(take!(output))
                @test success(process)
                @test !occursin("Error",log)
                back=read_msh(rewritten)
                @test mesh_crc(back.mesh).sha==mesh_crc(m).sha
                @test back.physical_names==names
                @test_throws ArgumentError read_msh(rewritten;max_blocks=1)
            end
        end

        @testset "Gmsh numeric list ranges" begin
            gmsh=Sys.which("gmsh")
            @test gmsh!==nothing
            if gmsh!==nothing
                gmsh_version=strip(read(`$gmsh --version`,String))
                @test gmsh_version=="4.15.2" ||
                      startswith(gmsh_version,"4.15.2-")
                function gmsh_range_output(code)
                    output=IOBuffer()
                    devnull=Sys.iswindows() ? "NUL" : "/dev/null"
                    command=ignorestatus(`$gmsh $devnull -parse_and_exit -string $code`)
                    run(pipeline(command;stdout=output,stderr=output))
                    return String(take!(output))
                end
                values=gmsh_range_output(raw"""
                    a[] = {1:5}; b[] = {5:1:-2};
                    c[] = {5:-2:1}; d[] = {5:-2};
                    Printf("R %g %g %g %g %g %g %g %g %g",
                           #a[], a[0], a[4], #b[], b[0], b[2], #c[], #d[], d[7]);
                    """)
                @test occursin("R 5 1 5 3 5 1 0 8 -2",values)
                accepted=gmsh_range_output(raw"""
                    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
                    Point(3)={2,0,0,1}; Point(4)={3,0,0,1};
                    Point(5)={4,0,0,1};
                    Physical Point("ranged", 11)={1:5};
                    Field[1]=Box; Field[2]=Box; Field[3]=Box;
                    Field[4]=Box; Field[5]=Box; Field[6]=Min;
                    Field[6].FieldsList={1:5};
                    Field[7]=BoundaryLayer; Field[8]=BoundaryLayer;
                    BoundaryLayer Field={7:8};
                    """)
                @test !occursin("Error",accepted)
                wrapped_accepted=gmsh_range_output(raw"""
                    Field[1]=BoundaryLayer;
                    Field[1].PointsList=-{-1:-3};
                    Field[1].SizesList=-{-0.1:-0.3:-0.1};
                    Field[3]=Box; Background Field=(1+2)*{1:1};
                    """)
                @test !occursin("Error",wrapped_accepted)
                modulo_values=gmsh_range_output(raw"""
                    a[]={(5.5%2.2):(5.5%2.2)};
                    Printf("MOD %g %g", #a[], a[0]);
                    """)
                @test occursin("MOD 1 1",modulo_values)
                precedence_rejected=gmsh_range_output(
                    "Field[1]=Box; Field[3]=Box; " *
                    "Background Field=1+2*{1:1};")
                @test occursin("syntax error",precedence_rejected)
                cardinality_rejected=gmsh_range_output(
                    "Field[1]=Box; Background Field=0.4*{3:4};")
                @test occursin("Only 1 field",cardinality_rejected)
                comment_rejected=gmsh_range_output(
                    "Field[1]=BoundaryLayer; " *
                    "Field[1].SizesList={1/* comment */2:12};")
                @test occursin("syntax error",comment_rejected)
                for source in (
                    "Field[1]=BoundaryLayer; BoundaryLayer Field=-{};",
                    "Field[1]=BoundaryLayer; BoundaryLayer Field=2*{};",
                    "Field[1]=Distance; Field[1].PointsList=-{};",
                    "Field[1]=BoundaryLayer; BoundaryLayer Field=;")
                    @test occursin("syntax error",gmsh_range_output(source))
                end
                rejected=gmsh_range_output("a[] = {1:5:0};")
                @test occursin("Wrong increment in '1:5:0'",rejected)
            end

        end
    end
end
end
