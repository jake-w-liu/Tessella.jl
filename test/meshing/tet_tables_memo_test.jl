# The suite includes every test file into Main; this file keeps its aliases to itself.
module TetTablesMemoTests
using Test
using Tessella
import Tessella.MeshTypes as MT
import Tessella.Mesh3D as M3

# While a Volume is meshed, the sorted face table and the topology certificate
# of a tetrahedron complex are kept by a fingerprint of its connectivity, so
# the many callers that ask for them share one computation. The memo must
# never hand out the table of another complex.

# A Kuhn grid of n^3 cubes, six positively oriented tetrahedra each.
function kuhn(n)
    id(i, j, k) = Int32(i + (n + 1) * (j + (n + 1) * k) + 1)
    coords = Matrix{Float64}(undef, 3, (n + 1)^3)
    for k in 0:n, j in 0:n, i in 0:n
        coords[:, id(i, j, k)] .= (i, j, k)
    end
    tets = Matrix{Int32}(undef, 4, 6n^3)
    column = 0
    for k in 0:n-1, j in 0:n-1, i in 0:n-1
        corner(bits) = id(i + (bits & 1), j + ((bits >> 1) & 1), k + ((bits >> 2) & 1))
        for (a, b) in ((1, 3), (1, 5), (2, 3), (2, 6), (4, 5), (4, 6))
            cell = (corner(0), corner(a), corner(b), corner(7))
            p = ntuple(v -> (coords[1, cell[v]], coords[2, cell[v]], coords[3, cell[v]]), 4)
            MT.tet_signed_volume(p...) < 0 && (cell = (cell[2], cell[1], cell[3], cell[4]))
            tets[:, column+=1] .= cell
        end
    end
    return coords, tets
end

@testset "Face tables and certificates are kept per complex while a Volume is meshed" begin
    coords, tets = kuhn(3)
    plain = M3._sorted_tet_faces(tets)
    @test issorted(plain) && length(plain) == 4 * size(tets, 2)
    # Outside a Volume nothing is kept.
    @test MT._tet_tables() === nothing
    @test M3._sorted_tet_faces(tets) !== plain

    began = MT._tet_tables_begin!()
    try
        @test began
        @test !MT._tet_tables_begin!()                    # nested: the outer call owns the memo
        first = M3._sorted_tet_faces(tets)
        @test first == plain
        @test M3._sorted_tet_faces(tets) === first        # same complex, same table
        @test M3._sorted_tet_faces(copy(tets)) === first  # by content, not by array
        @test M3._single_tet_faces(tets) == M3._single_tet_faces(plain)

        # Another complex gets its own table, and both stay available.
        fewer = tets[:, 1:end-6]
        other = M3._sorted_tet_faces(fewer)
        @test other == M3._build_sorted_tet_faces(fewer) && other !== first
        @test M3._sorted_tet_faces(tets) === first
        @test M3._sorted_tet_faces(fewer) === other

        # Editing a complex in place changes what is returned for it.
        edited = copy(tets)
        @test M3._sorted_tet_faces(edited) === first
        edited[1, 1], edited[2, 1] = edited[2, 1], edited[1, 1]      # same faces, other orientation
        swapped = M3._sorted_tet_faces(edited)
        @test swapped == first && swapped !== first
        edited[3, 5] = Int32(size(coords, 2))                        # other faces: the far corner of the grid
        @test M3._sorted_tet_faces(edited) == M3._build_sorted_tet_faces(edited) != first

        # The certificate is remembered for the complex that passed it and for no other.
        tables = MT._tet_tables()
        @test MT._tet_topology(tets) == (true, "ok")
        @test (size(tets, 2), MT._tet_content_key(tets)) in tables.certified
        @test MT._tet_topology(tets) == (true, "ok")
        broken = copy(tets); broken[:, 2] .= broken[:, 1]            # one cell listed twice
        @test !MT._tet_topology(broken)[1]
        @test !((size(broken, 2), MT._tet_content_key(broken)) in tables.certified)
        @test MT.validate(MT.Mesh(coords; tets=tets)).ok
        @test !MT.validate(MT.Mesh(coords; tets=broken)).ok
        # A pinch table changes the question, so the memo is not consulted.
        @test MT._tet_topology(tets; pinch_link_chi=Dict{Int32,Int}(Int32(1) => 0))[1] isa Bool
    finally
        MT._tet_tables_end!(began)
    end
    @test MT._tet_tables() === nothing

    # Every entry enters the fingerprint, in order.
    key = MT._tet_content_key(tets)
    for index in (1, 2, length(tets) ÷ 2, length(tets))
        changed = copy(tets); changed[index] += Int32(1)
        @test MT._tet_content_key(changed) != key
    end
    reordered = tets[:, reverse(axes(tets, 2))]
    @test MT._tet_content_key(reordered) != key
end

@testset "A Volume leaves no memo behind" begin
    api = Tessella.API
    api.initialize()
    try
        api.option("General.Terminal", 0)
        api.model.occ.add_box(0, 0, 0, 1, 1, 1); api.model.occ.synchronize()
        api.option("Mesh.MeshSizeMin", 0.4); api.option("Mesh.MeshSizeMax", 0.4)
        api.mesh.generate(3)
        @test MT._tet_tables() === nothing
    finally
        api.finalize()
    end
end
end
