using Test
using Tessella

@testset "Public documentation and dispatch" begin
    for mod in (
        Tessella,
        Tessella.API,
        Tessella.API.mesh,
        Tessella.BoundaryLayer,
        Tessella.CLI,
        Tessella.Elements,
        Tessella.GUI,
        Tessella.GeoExec,
        Tessella.HighOrder,
        Tessella.Mesh1D,
        Tessella.Mesh2D,
        Tessella.MeshElementQuality,
        Tessella.MeshEntityTopology,
        Tessella.MeshFunctionSpaces,
        Tessella.MeshPointLocation,
        Tessella.MeshReferenceGeometry,
        Tessella.MeshSurface,
        Tessella.MeshTypes,
        Tessella.Model,
        Tessella.Optimize,
        Tessella.Periodic,
        Tessella.Post,
        Tessella.Recombine,
        Tessella.Refine,
        Tessella.TransfinitePrism,
        Tessella.TransfiniteVolume,
        Tessella.Transform,
    )
        @testset "$mod" begin
            @test isempty(Base.Docs.undocumented_names(mod; private=false))
        end
    end
    @test isempty(Test.detect_ambiguities(
        Tessella.API,
        Tessella.Elements,
        Tessella.GeoExec,
        Tessella.HighOrder,
        Tessella.Mesh1D,
        Tessella.MeshElementQuality,
        Tessella.MeshEntityTopology,
        Tessella.MeshFunctionSpaces,
        Tessella.MeshPointLocation,
        Tessella.MeshQuadrature,
        Tessella.MeshReferenceGeometry,
        Tessella.MeshSurface,
        Tessella.Model,
        Tessella.Optimize,
        Tessella.Periodic,
        Tessella.Recombine,
        Tessella.Refine,
        Tessella.TransfiniteCurve,
        Tessella.TransfiniteHex,
        Tessella.TransfinitePrism,
        Tessella.TransfiniteQuad,
        Tessella.TransfiniteTriangle,
        Tessella.TransfiniteVolume,
        Tessella.Transform;
        recursive=true))
end
