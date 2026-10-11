module ExplicitPLCIncidenceChecks
using Test,Tessella
const M=Tessella.Model
function fixtures()
    coords=Float64[0 1 1 0 0 1 1 0;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1]
    tets=Int32[1 2 5 2 4 6;2 4 6 4 8 8;4 5 4 6 6 7;5 6 8 3 3 3]
    triangles=Int32[1 2 5 6 1 2 4 8 1 4 2 3;
                    4 4 6 7 2 6 8 7 5 5 3 7;
                    2 3 8 8 5 5 3 3 4 8 6 6]
    opposite=copy(triangles)
    opposite[:,1]=Int32[1,3,2];opposite[:,2]=Int32[1,4,3]
    return Mesh(coords;tris=triangles),Mesh(coords;tris=opposite),Mesh(coords;tets=tets)
end
@testset "actual PLC facet receipt distinguishes equal geometric sheets" begin
    source,opposite,volume=fixtures()
    @test validate(source).ok && validate(opposite).ok && validate(volume).ok
    # These are identical unit-cube geometric sheets. Their cell diagonals
    # differ, so only the recorded actual constructor facets can be retained.
    @test first(M._certify_surface_fill(source,volume))
    @test first(M._certify_surface_fill(opposite,volume))
    @test M._model_native_boundary_facets_match(source,volume)
    @test !M._model_native_boundary_facets_match(opposite,volume)
    snapshot=(copy(source.coords),copy(source.tris),copy(volume.coords),copy(volume.tets))
    changed=Mesh(source.coords;tris=source.tris[:,2:end])
    @test !M._model_native_boundary_facets_match(changed,volume)
    repeated=Mesh(source.coords;tris=hcat(source.tris,source.tris[:,1]))
    @test !M._model_native_boundary_facets_match(repeated,volume)
    invalid=deepcopy(volume);invalid.tets[1,1]=0
    @test !M._model_native_boundary_facets_match(source,invalid)
    invalid=deepcopy(volume);invalid.tets[1,1]=invalid.tets[2,1]
    @test !M._model_native_boundary_facets_match(source,invalid)
    invalid=Mesh(volume.coords;tets=hcat(volume.tets,volume.tets[:,1]))
    @test !M._model_native_boundary_facets_match(source,invalid)
    permuted=volume.tets[[2,1,3,4],end:-1:1]
    @test M._model_native_boundary_facets_match(source,Mesh(volume.coords;tets=permuted))
    @test snapshot==(source.coords,source.tris,volume.coords,volume.tets)
end
end
