module NativeHybridMshStatusChecks
using Test,Tessella
const M=Tessella.Model

function box_model()
    m=M.GeoModel()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points);M.add_point!(m,p...;tag,mesh_size=.4);end
    curves=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(curves);M.add_line!(m,ends...;tag);end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,curves) in enumerate(loops)
        M.add_curve_loop!(m,collect(curves);tag);M.add_plane_surface!(m,[tag];tag)
    end
    M.add_surface_loop!(m,collect(1:6);tag=1);M.add_volume!(m,[1];tag=1)
    return m
end
@testset "Clearing mesh authority preserves completed stages" begin
    @test M._model_mesh_cleared_authority(nothing)===nothing
    for history in ("cold","generate1","generate2_clear")
        m=box_model()
        done=Set{Tuple{Int,Int32}}()
        if history!="cold"
            union!(done,((1,Int32(tag)) for tag in 1:12))
            history=="generate2_clear" && union!(done,((2,Int32(tag)) for tag in 1:6))
        end
        previous=M._ModelMeshAuthority(nothing,Dict{Tuple{Int,Int32},M._ModelMeshProduct}(),done,UInt64[])
        cleared=M._model_mesh_cleared_authority(previous)
        @test cleared!==previous && cleared.done==done && cleared.done!==done
        @test cleared.mesh===nothing && isempty(cleared.products) && isempty(cleared.native_inputs)
        @test isempty(cleared.surface_constructor_inputs) && isempty(cleared.producer_rows)
        @test cleared.emitted_geometry===nothing && M._model_mesh_status(m,cleared)==-1
    end
end
end
