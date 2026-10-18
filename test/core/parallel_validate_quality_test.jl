# PAR01 slice — parallel O(cells) validation/quality scans under
# `General.NumThreads`. Every parallel scan writes disjoint per-element slots
# and reduces serially in traversal order, so diagnostics, qualities and
# checksums must be identical for any requested worker count. Run standalone
# (e.g. `julia -t4 test/core/parallel_validate_quality_test.jl`); the default
# suite does not register this file.

using Test
using Tessella
using Tessella.MeshTypes: Mesh, validate, mesh_crc, nsegs, ntris,
    ntets, unique_edges
using Tessella.Mesh3D: tetrahedralize
using Tessella.Geometry: box_surface
using Tessella.MeshElementQuality: mesh_element_qualities, ELEMENT_QUALITY_NAMES

A=Tessella.API
const _SEG_UNDEFINED=("minDetJac","maxDetJac","minSIGE","minIsotropy")

# A real produced volume mesh: the unit-cube boundary tetrahedralized by the
# native 3-D producer, plus a segment block rebuilt from its boundary edges so
# every cell kind is exercised.
function volume_fixture!(workers)
    A.finalize();A.initialize();A.model.add("parallel validate")
    A.option("General.Terminal",0)
    A.option("General.NumThreads",workers)
    vol=tetrahedralize(box_surface(0,1,0,1,0,1))
    @test validate(vol).ok && ntets(vol)>0
    edges=sort!(collect(unique_edges(vol.tris,vol.tets)))
    segs=Matrix{Int32}(undef,2,min(16,length(edges)))
    for (i,edge) in enumerate(edges)
        i>size(segs,2) && break
        segs[1,i],segs[2,i]=edge
    end
    return Mesh(vol.coords;segs=segs,tris=vol.tris,tets=vol.tets)
end

function all_qualities(mesh)
    tri_off=nsegs(mesh);tet_off=tri_off+ntris(mesh)
    total=tet_off+ntets(mesh)
    out=Dict{String,Vector{Float64}}()
    for name in ELEMENT_QUALITY_NAMES
        seg=name in _SEG_UNDEFINED ? Float64[] :
            mesh_element_qualities(mesh,collect(1:tri_off),name)
        rest=mesh_element_qualities(mesh,collect(tri_off+1:total),name)
        out[name]=vcat(seg,rest)
    end
    return out
end

@testset "parallel validation/quality stage stays deterministic" begin
    try
        serial=volume_fixture!(1)
        @test serial isa Mesh
        @test Tessella.MeshTypes._parallel_workers()==1
        serial_diag=validate(serial)
        @test serial_diag.ok
        serial_crc=mesh_crc(serial)
        serial_quality=all_qualities(serial)
        threaded=volume_fixture!(4)
        @test Tessella.MeshTypes._parallel_workers()==4
        @test threaded isa Mesh
        # Byte-identical product: the parallel scans cannot perturb generation.
        @test mesh_crc(threaded)==serial_crc
        threaded_diag=validate(threaded)
        @test threaded_diag.ok==serial_diag.ok
        @test threaded_diag.messages==serial_diag.messages
        threaded_quality=all_qualities(threaded)
        for name in ELEMENT_QUALITY_NAMES
            @test threaded_quality[name]==serial_quality[name]
        end
        # The same mesh evaluated under both worker counts reports identical
        # per-element quality values, bit for bit.
        tri_off=nsegs(serial);total=tri_off+ntris(serial)+ntets(serial)
        for name in ("minSICN","minDetJac","volume","gamma","minEdge")
            A.option("General.NumThreads",1)
            one=mesh_element_qualities(serial,collect(tri_off+1:total),name)
            A.option("General.NumThreads",4)
            two=mesh_element_qualities(serial,collect(tri_off+1:total),name)
            @test two==one
        end
    finally
        A.option("General.NumThreads",1)
        A.finalize()
    end
end

# A corrupted fixture exercises every flag class in the parallel scans:
# non-finite coordinate, negative tag, out-of-range node index, non-finite and
# zero measure cells, a negative tet, and duplicate segment/triangle keys. The
# Mesh constructor validates, so corruption is applied post-construction.
function corrupted_fixture()
    coords=Float64[0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 -1 0]
    segs=Int32[1 2;2 3]
    tris=Int32[1 2;2 3;3 4]
    tets=Int32[1 2;2 3;3 4;4 5]
    m=Mesh(coords;segs=segs,tris=tris,tets=tets)
    # Non-finite coordinate on the dangling node.
    m.coords[3,6]=NaN
    # Duplicate undirected segment: (1,2) already exists.
    m.segs[1,2]=Int32(1);m.segs[2,2]=Int32(2)
    # Zero-area duplicate triangle on the second slot.
    m.tris[1,2]=Int32(1);m.tris[2,2]=Int32(2);m.tris[3,2]=Int32(3)
    # Out-of-range node index in the first triangle.
    m.tris[1,1]=Int32(40)
    # Flip the second tet's orientation and corrupt its tag.
    m.tets[3,2],m.tets[4,2]=m.tets[4,2],m.tets[3,2]
    m.tet_tag[2]=Int32(-3)
    return m
end

# Structure-valid mesh that still trips the per-element scans: duplicate
# undirected segments, a zero-length segment, a zero-area duplicate triangle,
# a zero-volume tet, and a negatively-oriented tet. All ids stay in range and
# tags stay nonnegative so `_mesh_structure_messages!` does not return early.
function flagged_fixture()
    # Node 5 coincides with node 2, so cells referencing it measure zero.
    coords=Float64[0 1 0 0 1;0 0 1 0 0;0 0 0 1 0]
    # cols: (1,2) twice → one duplicate; (3,3) → zero length.
    segs=Int32[1 1 3;2 2 3]
    # cols: (1,2,3) twice → one duplicate; (1,2,5) → zero area.
    tris=Int32[1 1 1;2 2 2;3 3 5]
    # cols: (1,2,3,4) valid; (1,2,3,5) zero volume; (1,2,4,3) negative.
    tets=Int32[1 1 1;2 2 2;3 3 4;4 5 3]
    return Mesh(coords;segs=segs,tris=tris,tets=tets)
end

@testset "parallel diagnostics keep serial ordering and content" begin
    try
        A.finalize();A.initialize()
        m=corrupted_fixture()
        A.option("General.NumThreads",1)
        serial=validate(m;require_manifold_tris=false)
        @test !serial.ok && length(serial.messages)>=3
        A.option("General.NumThreads",4)
        threaded=validate(m;require_manifold_tris=false)
        @test threaded.ok==serial.ok
        @test threaded.messages==serial.messages
        # Degenerate-segment serialization flag follows the same ordering.
        serial2=validate(m;require_manifold_tris=false,allow_degenerate_segs=true)
        threaded2=validate(m;require_manifold_tris=false,allow_degenerate_segs=true)
        @test threaded2.messages==serial2.messages
        # Element-level flags: duplicate keys, zero measures, negative tets.
        f=flagged_fixture()
        A.option("General.NumThreads",1)
        serial3=validate(f;require_manifold_tris=false)
        @test !serial3.ok && length(serial3.messages)>=4
        A.option("General.NumThreads",4)
        threaded3=validate(f;require_manifold_tris=false)
        @test threaded3.ok==serial3.ok
        @test threaded3.messages==serial3.messages
        serial4=validate(f;require_manifold_tris=false,allow_degenerate_segs=true)
        threaded4=validate(f;require_manifold_tris=false,allow_degenerate_segs=true)
        @test threaded4.messages==serial4.messages
        @test length(serial4.messages)<length(serial3.messages)
    finally
        A.option("General.NumThreads",1)
        A.finalize()
    end
end
