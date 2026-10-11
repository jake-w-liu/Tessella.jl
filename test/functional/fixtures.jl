function unit_tetrahedron()
    return Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                tets=reshape(Int32[1, 2, 3, 4], 4, 1), tet_tag=Int32[7])
end

function unit_square()
    return Mesh(Float64[0 1 1 0; 0 0 1 1; 0 0 0 0];
                segs=Int32[1 2 3 4; 2 3 4 1],
                tris=Int32[1 1; 2 3; 3 4],
                seg_tag=Int32[1, 2, 3, 4], tri_tag=Int32[7, 7])
end

surface_area(mesh) = sum(MT.triangle_area(MT.node(mesh, a), MT.node(mesh, b),
                                        MT.node(mesh, c))
                         for (a, b, c) in eachcol(mesh.tris); init=0.0)

signed_volume(mesh) = sum(MT.tet_signed_volume(MT.node(mesh, a), MT.node(mesh, b),
                                            MT.node(mesh, c), MT.node(mesh, d))
                         for (a, b, c, d) in eachcol(mesh.tets); init=0.0)
