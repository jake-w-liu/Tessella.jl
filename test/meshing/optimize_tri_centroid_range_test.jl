module SurfaceCentroidRangeTests
using Test, Tessella
const MT=Tessella.MeshTypes
const A=Tessella.API

function fixture(offset)
    width=iszero(offset) ? 1.0 : 16eps(abs(offset))
    coords=Float64[offset offset+width offset+width offset offset+width/4;
                   0 0 1 1 .25; 0 0 0 0 0]
    tris=Int32[1 2 3 4; 2 3 4 1; 5 5 5 5]
    return MT.Mesh(coords;tris,tri_tag=Int32[7,8,9,10])
end

@testset "Surface Laplacian centroids retain finite coordinate range" begin
    cancellation=MT.Mesh(Float64[1e308 1e308 -1e308 -1e308 .5e308;
                                  0 .25 .25 0 .125; 0 0 0 0 0];
        tris=Int32[1 2 3 4;2 3 4 1;5 5 5 5],tri_tag=Int32[7,8,9,10])
    for source in (fixture(0.0),fixture(1e308),fixture(-1e308),cancellation)
        snapshot=copy(source.coords)
        @test MT.validate(source).ok
        exact=ntuple(3) do axis
            Float64(sum(Rational{BigInt}(source.coords[axis,node]) for node in 1:4)/4)
        end
        for iterations in (1,3)
            out=Tessella.Optimize._laplacian_smooth_tri_cache(source,iterations,nothing)
            @test MT.validate(out).ok
            @test out.coords[:,1:4]==snapshot[:,1:4]
            @test out.coords[1,5]==exact[1]
            @test Tuple(out.coords[2:3,5])==exact[2:3]
            @test out.tris==source.tris && out.tri_tag==source.tri_tag
            @test source.coords==snapshot
        end
        held=Tessella.Optimize._laplacian_smooth_tri_cache(source,1,falses(5))
        @test held.coords==snapshot
        for method in ("Laplace2D","Relocate2D")
            A.initialize()
            try
                A._replace_mesh_cache_locked!(source)
                @test A.mesh.optimize(method)===nothing
                tags,coords,_=A.mesh.get_nodes()
                @test length(tags)==5 && all(isfinite,coords)
                @test reshape(coords,3,:)[:,1:4]==snapshot[:,1:4]
                @test MT.validate(A.LAST_MESH[]).ok
            finally
                A.finalize()
            end
        end
    end
    ordinary=fixture(0.0)
    ordinary.coords[1,:].=[.1,1.7,1.9,-.3,.4]
    sums=zeros(3)
    for node in 1:4,axis in 1:3
        sums[axis]+=ordinary.coords[axis,node]
    end
    out=Tessella.Optimize._laplacian_smooth_tri_cache(ordinary,1,nothing)
    @test reinterpret(UInt64,out.coords[:,5])==reinterpret(UInt64,sums./4)
end
end
