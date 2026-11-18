using Test
using Tessella
import Tessella.API as A
import Tessella.Mesh3D as M3

# The vertices of a sphere mesh are cospherical, so their Delaunay
# tetrahedralization is decided by rounding and misses surface edges at some
# sizes. A missing edge that runs under a hull edge sits in a thin twisted
# cavity no refill can fill; it is recovered once a vertex is seeded on the
# inner side. Recovery used to seed one side only, could pick the outer one,
# and then left the whole fill to the slow float recovery.

function sphere_surface(h)
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.model.occ.add_sphere(0, 0, 0, 1)
        A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin", h)
        A.option("Mesh.MeshSizeMax", h)
        A.mesh.generate(2)
        tags, xyz, _ = A.mesh.get_nodes()
        row = Dict(tag => Int32(i) for (i, tag) in enumerate(tags))
        types, _, nodes = A.mesh.get_elements(2)
        triangles = Int32[]
        for i in eachindex(types)
            types[i] == 2 && append!(triangles, (row[node] for node in nodes[i]))
        end
        return Matrix{Float64}(reshape(xyz, 3, :)), Matrix{Int32}(reshape(triangles, 3, :))
    finally
        A.finalize()
    end
end

@testset "Sphere surfaces are recovered without the float fallback" begin
    for h in (0.5, 0.4, 0.35, 0.25, 0.22, 0.2)
        coords, tris = sphere_surface(h)
        fill, reason, stats = Tessella.BoundaryRecovery.recover_boundary_flips(coords, tris)
        @test fill !== nothing
        fill === nothing && (@info "sphere recovery" h reason stats; continue)
        @test validate(fill).ok
        @test M3._fill_boundary_is_surface(Mesh(coords; tris=tris), fill)
        # Added vertices are interior: the surface rows are kept as given.
        @test fill.coords[:, 1:size(coords, 2)] == coords
        @test size(fill.coords, 2) - size(coords, 2) <= 4
    end
end
