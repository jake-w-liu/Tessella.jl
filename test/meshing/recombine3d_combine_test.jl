# Recombine3D star fixtures — `Recombine3D.combine_volume_to_hexes` probed
# directly on hand-built tetrahedral stars (the combine arm itself, decoupled
# from the Delaunay fill):
#
#   * 2-tet pair and 3-tet edge star — too few cells for any hexahedron;
#   * 5-tet and 6-tet cube splits — the canonical star the vertex walk must
#     recover whole, with boundary-quadrangle receipts on a corner star;
#   * two adjacent hex stars — the compatible-selection facet sharing;
#   * mixed-quality prism — a sheared star rejected by the quality gate while
#     its unit-cube neighbour still combines (driving the odd-incidence
#     trihedron conformalizer);
#   * a folded (inverted) star — the trilinear-Bezier certificate rejects;
#   * forbidden (embedded-sheet) half-triangles — the colour guard kills the
#     candidate;
#   * retained stray tets — unconsumed cells carry over untouched.
#
# Oracle reference points (Gmsh 4.15.2, `Mesh.Algorithm3D=9` + `Recombine3D`):
# only type-5/type-4/type-140 cells are ever emitted — upstream's
# `selectCellsGreedy` enables hexes and retained tets only — and a tetrahedron
# domain at lc=0.4 reproduces 1 hex + 47 tets + 5 trihedra + 1 boundary quad
# bit-for-bit in counts on the Tessella side.
module Recombine3DCombineChecks

using Test, Tessella
const R3=Tessella.Recombine3D
const E=Tessella.Elements

const _CUBE=Float64[0 0 0; 1 0 0; 1 1 0; 0 1 0;
                    0 0 1; 1 0 1; 1 1 1; 0 1 1]'

# Canonical 5-tet and 6-tet cube splits (see STATUS): both cover the cube
# conformally; the 6-tet split fans the 1-7 body diagonal.
const _TETS5=Int32[1 2 3 6; 1 3 4 8; 1 5 6 8; 3 6 7 8; 1 3 6 8]'
const _TETS6=Int32[1 2 3 7; 1 3 4 7; 1 4 8 7; 1 8 5 7; 1 5 6 7; 1 6 2 7]'

_orient6(a,b,c,d)=((b[1]-a[1])*((c[2]-a[2])*(d[3]-a[3])-(c[3]-a[3])*(d[2]-a[2]))-
                  (b[2]-a[2])*((c[1]-a[1])*(d[3]-a[3])-(c[3]-a[3])*(d[1]-a[1]))+
                  (b[3]-a[3])*((c[1]-a[1])*(d[2]-a[2])-(c[2]-a[2])*(d[1]-a[1])))

function _positive_tets(coords,tets)
    out=copy(tets)
    for j in axes(out,2)
        a,b,c,d=(coords[:,out[i,j]] for i in 1:4)
        _orient6(a,b,c,d)<0 && ((out[3,j],out[4,j])=(out[4,j],out[3,j]))
    end
    return out
end

function _fixture(tets)
    t=_positive_tets(_CUBE,tets)
    return Tessella.Mesh(_CUBE;tets=t)
end

# Boundary colouring: map every mesh-boundary triangle (sorted triple) to a
# surface tag — one tag per constant-coordinate face, or a single tag for all.
function _owner_by_face(mesh,per_face::Bool)
    m=R3._rtet_mesh(mesh,nothing,nothing)
    out=Dict{NTuple{3,Int32},Int32}()
    face_id=0;seen=Dict{Tuple,Int}()
    for tri in m.boundary_tris
        if per_face
            p=mesh.coords[:,tri[1]]
            c=mesh.coords[:,tri[2]]
            q=mesh.coords[:,tri[3]]
            # Plane the triangle lies on: constant (axis, value) pair —
            # x=0 and x=1 faces must not share a key.
            key=Tuple((i,p[i]) for i in 1:3 if p[i]==c[i]==q[i])
            fid=get!(seen,key,length(seen)+1)
            out[tri]=Int32(fid)
        else
            out[tri]=Int32(1)
        end
    end
    return out
end

function _inventory(m::Tessella.MixedMesh)
    counts=Dict{Int,Int}()
    for b in m.blocks
        counts[b.msh]=get(counts,b.msh,0)+size(b.nodes,2)
    end
    return counts
end

@testset "pairs and small stars cannot form a hexahedron" begin
    # 2-tet pair sharing a face.
    m=Tessella.Mesh(_CUBE;tets=_positive_tets(_CUBE,Int32[1 2 3 4; 2 3 4 6]'))
    @test R3.combine_volume_to_hexes(m)===nothing
    # 3-tet star around edge (1,2).
    m=Tessella.Mesh(_CUBE;tets=_positive_tets(_CUBE,
        Int32[1 2 3 4; 1 2 4 6; 1 2 6 3]'))
    @test R3.combine_volume_to_hexes(m)===nothing
    # A prism's 3-tet split — pyramids/prisms are never emitted either.
    prism=Float64[0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1]'
    t=Int32[1 2 3 5; 1 3 4 5; 3 4 5 6]'
    m=Tessella.Mesh(prism;tets=_positive_tets(prism,t))
    @test R3.combine_volume_to_hexes(m)===nothing
end

@testset "5-tet and 6-tet cube stars recover the hexahedron" begin
    for tets in (_TETS5,_TETS6)
        mesh=_fixture(tets)
        owners=_owner_by_face(mesh,true)
        res=R3.combine_volume_to_hexes(mesh,owners)
        @test res!==nothing
        inv=_inventory(res.mesh)
        @test inv==Dict(5=>1)
        # Every corner survives; the hex is the certified trilinear cell.
        @test length(res.kept_rows)==8
        hex=[res.mesh.blocks[1].nodes[:,i] for i in 1:1][1]
        @test sort(hex)==collect(Int32(1):Int32(8))
        # All six boundary facets pair to quadrangles — six distinct faces.
        @test length(res.quads)==6
        @test sort(res.quad_surfaces)==collect(Int32(1):Int32(6))
        @test length(res.halves)==6
        @test res.unpaired_boundary_quads==0
        @test Tessella.validate(res.mesh).ok
    end
    # Same star, uncoloured boundary: candidates combine but no quadrangle
    # receipt can be emitted — the halves count as unpaired.
    mesh=_fixture(_TETS6)
    res=R3.combine_volume_to_hexes(mesh)
    @test res!==nothing
    @test _inventory(res.mesh)==Dict(5=>1)
    @test isempty(res.quads)
    @test res.unpaired_boundary_quads==6
    # One mismatched owner pair: the `checkQuadFacetColor` legality audit
    # rejects the whole candidate — a hex quad facet may never span two
    # differently-coloured surface halves (upstream semantics, unlike two
    # uncoloured halves which merely stay triangulated).
    owners=_owner_by_face(mesh,false)
    one=first(keys(owners))
    owners[one]=Int32(9)
    res=R3.combine_volume_to_hexes(mesh,owners)
    @test res===nothing
end

@testset "two adjacent hex stars select compatibly" begin
    # 2x1x1 prism, twelve nodes: left hex (1..8) + right hex (5..8,9..12)
    # sharing the (5,6,7,8) quad facet with matching diagonal 5-7.
    coords=Float64[0 0 0; 1 0 0; 1 1 0; 0 1 0;
                   0 0 1; 1 0 1; 1 1 1; 0 1 1;
                   0 0 2; 1 0 2; 1 1 2; 0 1 2]'
    tets=Int32[1 2 3 7; 1 3 4 7; 1 4 8 7; 1 8 5 7; 1 5 6 7; 1 6 2 7;
               5 6 7 11; 5 7 8 11; 5 8 12 11; 5 12 9 11; 5 9 10 11; 5 10 6 11]'
    mesh=Tessella.Mesh(coords;tets=_positive_tets(coords,tets))
    res=R3.combine_volume_to_hexes(mesh,_owner_by_face(mesh,false))
    @test res!==nothing
    inv=_inventory(res.mesh)
    @test inv==Dict(5=>2)
    # Ten boundary facets pair (12 facets minus the two coplanar shared faces).
    @test length(res.quads)==10
    @test res.unpaired_boundary_quads==0
    @test Tessella.validate(res.mesh).ok
end

@testset "mixed-quality prism: gate selects only the certified star" begin
    # Left unit-cube star (corner quality 1.0) + right hex star whose top quad
    # is sheared by 1.5 in x — its corner scaled Jacobian falls under the
    # gate at min_quality=0.9, so only the cube combines and the shared facet
    # is conformalized by a type-140 trihedron.
    coords=Float64[0 0 0; 1 0 0; 1 1 0; 0 1 0;
                   0 0 1; 1 0 1; 1 1 1; 0 1 1;
                   1.5 0 2; 2.5 0 2; 2.5 1 2; 1.5 1 2]'
    tets=Int32[1 2 3 7; 1 3 4 7; 1 4 8 7; 1 8 5 7; 1 5 6 7; 1 6 2 7;
               5 6 7 11; 5 7 8 11; 5 8 12 11; 5 12 9 11; 5 9 10 11; 5 10 6 11]'
    mesh=Tessella.Mesh(coords;tets=_positive_tets(coords,tets))
    res=R3.combine_volume_to_hexes(mesh,_owner_by_face(mesh,false);
                                   min_quality=0.9)
    @test res!==nothing
    inv=_inventory(res.mesh)
    @test get(inv,5,0)==1
    @test get(inv,4,0)>0
    @test get(inv,140,0)==1
    @test Tessella.validate(res.mesh).ok
    # At the default threshold the sheared star is certified too — both
    # combine (the certificate is the strong gate, the threshold only prunes
    # near-degenerate corners).
    res2=R3.combine_volume_to_hexes(mesh,_owner_by_face(mesh,false))
    @test res2!==nothing
    @test get(_inventory(res2.mesh),5,0)==2
end

@testset "folded star fails the Bezier certificate" begin
    # Vertex 7 folded below the base plane — every star tet inverts.
    coords=copy(_CUBE);coords[3,7]=-0.2
    mesh=Tessella.Mesh(coords;tets=_positive_tets(coords,_TETS6))
    @test R3.combine_volume_to_hexes(mesh)===nothing
end

@testset "forbidden halves kill the covering candidate" begin
    mesh=_fixture(_TETS6)
    m=R3._rtet_mesh(mesh,nothing,nothing)
    forbidden=Set{NTuple{3,Int32}}([first(m.boundary_tris)])
    res=R3.combine_volume_to_hexes(mesh,nothing;forbidden_tris=forbidden)
    @test res===nothing
end

@testset "unconsumed tets carry over untouched" begin
    # Hex star plus a disjoint stray tet the census cannot fold.
    coords=hcat(_CUBE,Float64[2 0 0; 3 0 0; 3 1 0; 2 0 1]')
    tets=hcat(_TETS6,Int32[9 10 11 12]')
    mesh=Tessella.Mesh(coords;tets=_positive_tets(coords,tets))
    res=R3.combine_volume_to_hexes(mesh,_owner_by_face(mesh,false))
    @test res!==nothing
    inv=_inventory(res.mesh)
    @test get(inv,5,0)==1
    @test get(inv,4,0)==1
    @test Tessella.validate(res.mesh).ok
end

end # module
