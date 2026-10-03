using Test
using Tessella
using Tessella.IO: read_geo_params, read_stl

@testset "failed file readers release their streams immediately" begin
    facet="""
        solid bounded
        facet normal 0 0 1
        outer loop
        vertex 0 0 0
        vertex 1 0 0
        vertex 0 1 0
        endloop
        endfacet
        endsolid bounded
        """
    cases=(
        (execute_geo,"}\n"),
        (execute_geo,"If (1\nEndIf\n"),
        (read_geo_params,"If (1\nEndIf\n"),
        (read_stl,"solid malformed\nvertex 0 0 0\nendsolid malformed\n"),
        (read_stl,"solid malformed\nfacet normal invalid 0 1\n"),
        (path->read_stl(path;max_facets=0),facet),
        (path->Tessella.IO._scan_geo_statements(
            _->throw(ArgumentError("aborted statement consumer")),path),
         "a = 1;\nb = 2;\n"))
    mktempdir() do directory
        for (index,(reader,source)) in enumerate(cases)
            path=joinpath(directory,"rejected_$index.input")
            write(path,source)
            @test_throws ArgumentError reader(path)
            # Windows refuses this unlink while the reader still owns a file
            # handle. Do not run GC between rejection and removal.
            @test begin
                rm(path)
                !ispath(path)
            end
        end
    end
end
