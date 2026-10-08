# ── Stage-4 CRC suite: tet-mesh quality reporting + Laplacian & ODT smoothing ───
#
# Correctness  : quality metrics on a known-good mesh (regular-ish tets); smoothing
#                preserves total volume (boundary fixed) and validity; improves the
#                mean dihedral and does not increase the sliver count. ODT is checked
#                against an INDEPENDENT volume-weighted-circumcenter oracle that a
#                centroid/Lloyd update would fail (it discriminates the algorithm).
# Robustness   : slivery Delaunay-of-random-cloud input, boundary nodes pinned.
# Completeness  : no tet inverted (positive-volume guard), tags preserved.

using Test
import Tessella
using Tessella.Mesh3D
using Tessella.MeshTypes
using Tessella.Optimize

mvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)

mutable struct _RO; s::UInt64; end
_nfo(r::_RO) = (r.s ⊻= r.s<<13; r.s ⊻= r.s>>7; r.s ⊻= r.s<<17; (r.s>>11)/Float64(2^53))

_dist(p, q) = sqrt((p[1]-q[1])^2 + (p[2]-q[2])^2 + (p[3]-q[3])^2)

# INDEPENDENT circumcenter (Cramer linear solve of |x−vᵢ|² equal), distinct from the
# production scale-normalized `_odt_candidate` solve — the oracle for the ODT update.
function cramer_circumcenter(a,b,c,d)
    v1=collect(Float64,a); v2=collect(Float64,b); v3=collect(Float64,c); v4=collect(Float64,d)
    A=[ (v2.-v1)'; (v3.-v1)'; (v4.-v1)' ].*2
    rhs=[sum(v2.^2)-sum(v1.^2), sum(v3.^2)-sum(v1.^2), sum(v4.^2)-sum(v1.^2)]
    d3(M)=M[1,1]*(M[2,2]*M[3,3]-M[2,3]*M[3,2])-M[1,2]*(M[2,1]*M[3,3]-M[2,3]*M[3,1])+M[1,3]*(M[2,1]*M[3,2]-M[2,2]*M[3,1])
    D=d3(A); (d3([rhs A[:,2] A[:,3]])/D, d3([A[:,1] rhs A[:,3]])/D, d3([A[:,1] A[:,2] rhs])/D)
end

# Twelve cones over the outward-oriented unit-cube faces, with an off-center
# interior vertex. Reversing each face supplies positive tetrahedron order
# without using a production orientation predicate to build the oracle fixture.
function _optimize_tagged_cube_star()
    coordinates=Float64[0 1 1 0 0 1 1 0 .12;
                        0 0 1 1 0 0 1 1 .34;
                        0 0 0 0 1 1 1 1 .71]
    faces=Int32[1 1 5 5 1 1 2 2 3 3 4 4;
                3 4 6 7 2 6 3 7 4 8 1 5;
                2 3 7 8 6 5 7 6 8 7 5 8]
    edges=Int32[1 2 3 4 5 6 7 8 1 2 3 4;
                2 3 4 1 6 7 8 5 5 6 7 8]
    cells=vcat(faces[[1,3,2],:],fill(Int32(9),1,12))
    return Mesh(coordinates;segs=edges,tris=faces,tets=cells,
                seg_tag=Int32.(1:12),tri_tag=Int32.(21:32),
                tet_tag=Int32.(41:52))
end

# Independent rational determinant sum, rather than the production volume
# function: fixed cube boundary and positive cone cells must total exactly one.
function _optimize_exact_volume(mesh)
    R=Rational{BigInt}
    total=zero(R)
    for cell in eachcol(mesh.tets)
        a,b,c,d=cell
        u=R.(mesh.coords[:,b])-R.(mesh.coords[:,a])
        v=R.(mesh.coords[:,c])-R.(mesh.coords[:,a])
        w=R.(mesh.coords[:,d])-R.(mesh.coords[:,a])
        total+=(u[1]*(v[2]*w[3]-v[3]*w[2])-
                u[2]*(v[1]*w[3]-v[3]*w[1])+
                u[3]*(v[1]*w[2]-v[2]*w[1]))/6
    end
    return total
end

@testset "optimization returns independent buffers and preserves the unit cube" begin
    input=_optimize_tagged_cube_star()
    buffers=fieldnames(Mesh)
    snapshot=map(field->copy(getfield(input,field)),buffers)
    @test validate(input).ok
    @test _optimize_exact_volume(input)==1
    variants=(
        (mesh->smooth_laplacian(mesh;iters=0)),
        (mesh->smooth_laplacian(mesh;iters=1,relax=.5)),
        (mesh->smooth_laplacian(mesh;iters=3,relax=.9)),
        (mesh->smooth_odt(mesh;iters=0)),
        (mesh->smooth_odt(mesh;iters=1)),
        (mesh->smooth_odt(mesh;iters=3)),
        (mesh->smooth_optimize(mesh;iters=0)),
        (mesh->smooth_optimize(mesh;iters=1,sliver_deg=35.)),
        (mesh->smooth_optimize(mesh;iters=2,sliver_deg=35.,
                               movable=trues(nnodes(mesh)))),
        (mesh->smooth_optimize(mesh;iters=2,sliver_deg=35.,
                               movable=falses(nnodes(mesh)))),
        (mesh->first(remove_slivers(mesh;max_rounds=0))),
        (mesh->first(remove_slivers(mesh;max_rounds=1,sliver_deg=35.))),
        (mesh->first(remove_slivers(mesh;max_rounds=3,sliver_deg=35.))),
    )
    for optimize in variants
        result=optimize(input)
        @test validate(result).ok
        @test _optimize_exact_volume(result)==1
        @test result.coords[:,1:8]==input.coords[:,1:8]
        # Independent SHA-256 of the literal sorted, padded little-endian
        # connectivity records and family/count delimiters of this fixture.
        @test mesh_crc(result).sha==
            "b52840703e69c865b608b1671cf1c290616f62b9c49e7d424f68567fd7e2524d"
        for field in buffers
            @test getfield(result,field)!==getfield(input,field)
            field===:coords || (@test getfield(result,field)==getfield(input,field))
        end
        @test all(getfield(input,field)==saved for (field,saved) in zip(buffers,snapshot))
        # Exercise every nonempty returned buffer, including metadata. Object
        # identity alone would miss distinct array headers sharing storage.
        result.coords[1,1]+=.125
        for field in (:segs,:tris,:tets)
            matrix=getfield(result,field)
            matrix[1,1]=mod(matrix[1,1],nnodes(result))+1
        end
        for field in (:seg_tag,:tri_tag,:tet_tag)
            getfield(result,field)[1]+=1
        end
        @test all(getfield(input,field)==saved for (field,saved) in zip(buffers,snapshot))
    end
end

# A tagged segment and triangle in each disjoint right tetrahedron keep every
# serialized buffer nonempty. All vertices are on the boundary, so these cases
# isolate public validation/setup/output costs while still using real kernels.
function _optimize_allocation_mesh(count::Int)
    coordinates=zeros(3,4count)
    for index in 0:count-1
        first=4index+1;x=3.0index
        coordinates[:,first].=(x,0.,0.)
        coordinates[:,first+1].=(x+1,0.,0.)
        coordinates[:,first+2].=(x,1.,0.)
        coordinates[:,first+3].=(x,0.,1.)
    end
    cells=reshape(Int32.(1:4count),4,:)
    return Mesh(coordinates;tets=cells,segs=cells[1:2,:],tris=cells[1:3,:],
                seg_tag=fill(Int32(7),count),tri_tag=fill(Int32(8),count),
                tet_tag=fill(Int32(9),count))
end

function _optimize_warmed_allocation(operation,mesh)
    operation(mesh)
    minimum_bytes=typemax(Int)
    for _ in 1:5
        minimum_bytes=min(minimum_bytes,@allocated operation(mesh))
    end
    return minimum_bytes
end

@testset "optimization allocation budgets for tagged boundary meshes" begin
    for count in (1000,2000,4000)
        mesh=_optimize_allocation_mesh(count)
        @test validate(mesh).ok
        # One detached payload plus 10% and 4 KiB for array headers/setup.
        # The former redundant copy used approximately twice the payload.
        payload=sum(sizeof(getfield(mesh,field)) for field in fieldnames(Mesh))
        @test _optimize_warmed_allocation(Tessella.Optimize._copy_mesh,mesh)<=
              ceil(Int,1.1payload)+4096
        # Warmed Julia 1.12/1.13 measurements at 1k/2k/4k cells establish these
        # linear budgets with fixed setup slack. Reintroducing the old full
        # connectivity/tag pre-copy exceeds each unchanged budget at every size.
        @test _optimize_warmed_allocation(
            m->smooth_laplacian(m;iters=1),mesh)<=2660count+16384
        @test _optimize_warmed_allocation(
            m->smooth_odt(m;iters=1),mesh)<=1260count+16384
        @test _optimize_warmed_allocation(
            m->smooth_optimize(m;iters=1),mesh)<=1260count+16384
        @test _optimize_warmed_allocation(
            m->remove_slivers(m;max_rounds=1),mesh)<=300count+16384
    end
end

@testset "Optimize (Stage 4)" begin

    @testset "public parameter and metadata contracts" begin
        C = Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
        m = Mesh(C; segs=reshape(Int32[1,2],2,1),
                 tris=reshape(Int32[1,2,3],3,1),
                 tets=reshape(Int32[1,2,3,4],4,1),
                 seg_tag=Int32[7], tri_tag=Int32[8], tet_tag=Int32[9])
        @test_throws ArgumentError mesh_quality(m; sliver_deg=NaN)
        @test_throws ArgumentError mesh_quality(m; sliver_deg=-1)
        @test_throws ArgumentError mesh_quality(m; sliver_deg=true)
        @test_throws ArgumentError smooth_laplacian(m; iters=-1)
        @test_throws ArgumentError smooth_laplacian(m; iters=true)
        @test_throws ArgumentError smooth_laplacian(m; relax=NaN)
        @test_throws ArgumentError smooth_laplacian(m; relax=1.1)
        @test_throws ArgumentError smooth_laplacian(m; relax=true)
        @test_throws ArgumentError smooth_odt(m; iters=-1)
        @test_throws ArgumentError smooth_odt(m; iters=true)
        @test_throws ArgumentError smooth_optimize(m; iters=-1)
        @test_throws ArgumentError smooth_optimize(m; iters=true)
        @test_throws ArgumentError smooth_optimize(m; sliver_deg=Inf)
        @test_throws ArgumentError smooth_optimize(m; sliver_deg=true)
        @test_throws ArgumentError remove_slivers(m; max_rounds=-1)
        @test_throws ArgumentError remove_slivers(m; max_rounds=true)
        @test_throws ArgumentError remove_slivers(m; sliver_deg=true)
        @test_throws ArgumentError smooth_laplacian(m; iters=big(typemax(Int))+1)
        @test_throws ArgumentError smooth_odt(m; iters=big(typemax(Int))+1)
        @test_throws ArgumentError smooth_optimize(m; iters=big(typemax(Int))+1)
        @test_throws ArgumentError remove_slivers(m; max_rounds=big(typemax(Int))+1)
        @test Tessella.Optimize._convex_combine(-floatmax(Float64),floatmax(Float64),0.5) == 0.0
        near_cancel=Tessella.Optimize._convex_combine(
            nextfloat(-1e308),1e308,0.5)
        cancel_reference=setprecision(BigFloat,256) do
            Float64((BigFloat(nextfloat(-1e308))+BigFloat(1e308))/2)
        end
        @test near_cancel==cancel_reference
        @test TetQuality(1,10,20,30,1,2,0.1,0).n_tets==1
        @test_throws ArgumentError TetQuality(true,0,0,0,0,0,0,0)
        @test_throws ArgumentError TetQuality(1,true,20,30,1,2,0.1,0)
        @test_throws ArgumentError TetQuality(1,10,20,30,1,2,NaN,0)
        @test_throws ArgumentError TetQuality(big(typemax(Int))+1,10,20,30,1,2,0.1,0)
        @test_throws ArgumentError TetQuality(0,0,0,0,0,0,1,0)
        @test_throws ArgumentError TetQuality(1,30,20,40,1,2,0.1,0)
        @test_throws ArgumentError TetQuality(1,10,20,30,2,1,0.1,0)
        @test_throws ArgumentError TetQuality(1,10,20,30,1,2,0.1,2)
        for smoothed in (smooth_laplacian(m; iters=0), smooth_odt(m; iters=0),
                         smooth_optimize(m; iters=0))
            @test smoothed.segs == m.segs && smoothed.seg_tag == m.seg_tag
            @test smoothed.tris == m.tris && smoothed.tri_tag == m.tri_tag
            @test smoothed.tets == m.tets && smoothed.tet_tag == m.tet_tag
        end
        invalid = Mesh(C; tets=reshape(Int32[1,2,4,3],4,1))
        @test_throws ArgumentError mesh_quality(invalid)
        @test_throws ArgumentError smooth_laplacian(invalid)
        @test_throws ArgumentError smooth_odt(invalid)
        @test_throws ArgumentError smooth_optimize(invalid)
        nonfinite = Mesh(C; tets=reshape(Int32[1,2,3,4],4,1))
        nonfinite.coords[1,1] = NaN
        @test_throws ArgumentError smooth_laplacian(nonfinite)
        invalid_notets=Mesh(Float64[0 0;0 0;0 0];segs=reshape(Int32[1,2],2,1))
        @test_throws ArgumentError mesh_quality(invalid_notets)
        @test_throws ArgumentError smooth_odt(invalid_notets)
        detached,report=remove_slivers(m;max_rounds=0)
        @test detached!==m
        @test all(getfield(detached,field)!==getfield(m,field) for field in
                  (:coords,:segs,:tris,:tets,:seg_tag,:tri_tag,:tet_tag))
        @test report.slivers_before==report.slivers_after
    end

    @testset "mesh_quality on a single unit tet" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,4],4,1))
        q = mesh_quality(m)
        @test q.n_tets == 1
        @test q.min_volume ≈ 1/6 rtol=1e-12
        @test 0 < q.min_dihedral_deg <= q.mean_dihedral_deg <= q.max_dihedral_deg < 180
        @test q.n_slivers == 0                      # a unit right tet is not a sliver
    end

    @testset "smoothing preserves volume + validity, improves mean dihedral" begin
        r = _RO(0xABCDEF)
        n = 300
        xs=Float64[_nfo(r) for _ in 1:n]; ys=Float64[_nfo(r) for _ in 1:n]; zs=Float64[_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        q0 = mesh_quality(m)
        ms = smooth_laplacian(m; iters=10, relax=0.8)
        q1 = mesh_quality(ms)
        @test mvol(ms) ≈ mvol(m) rtol=1e-9          # boundary fixed ⇒ volume preserved
        @test validate(ms).ok                       # no inverted tets
        @test ntets(ms) == ntets(m)                 # topology unchanged
        @test q1.mean_dihedral_deg > q0.mean_dihedral_deg      # mean quality up
        @test q1.n_slivers <= q0.n_slivers          # slivers not increased
    end

    @testset "boundary nodes are pinned (cube unchanged)" begin
        # cube surface fill has only boundary nodes → smoothing is a no-op geometrically
        C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
        F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
        ct=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); ct[:,k]=Int32[f...]; end
        m = tetrahedralize(Mesh(C; tris=ct))
        ms = smooth_laplacian(m; iters=5)
        @test ms.coords ≈ m.coords                  # all nodes on the boundary ⇒ fixed
        @test mvol(ms) ≈ 1.0 rtol=1e-6
    end

    @testset "tags preserved through smoothing" begin
        r = _RO(0x1234)
        n = 120
        xs=Float64[_nfo(r) for _ in 1:n]; ys=Float64[_nfo(r) for _ in 1:n]; zs=Float64[_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        # slap arbitrary tags on
        tags = Int32[(t % 3) + 1 for t in 1:ntets(m)]
        mt = Mesh(m.coords; tets=m.tets, tet_tag=tags)
        ms = smooth_laplacian(mt; iters=3)
        @test ms.tet_tag == tags
    end

    @testset "ODT update = volume-weighted circumcenter average (discriminates from Lloyd)" begin
        # A stretched box (8 corners) with ONE off-centre interior node. One ODT step
        # must move that node to the volume-weighted average of its incident tets'
        # circumcenters — a target the centroid/Lloyd update does NOT hit.
        C = Float64[0 3 3 0 0 3 3 0 1.1; 0 0 1 1 0 0 1 1 0.35; 0 0 0 0 1 1 1 1 0.65]
        m = to_mesh3(delaunay3d(C[1,:], C[2,:], C[3,:]; rng_seed=1))
        bf,_ = boundary_faces(m.tets); onb = falses(nnodes(m))
        for f in bf; onb[f[1]]=true; onb[f[2]]=true; onb[f[3]]=true; end
        interior = findall(!, onb)
        @test length(interior) == 1                       # exactly one movable node
        v = interior[1]
        inc = [t for t in 1:ntets(m) if v in (m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])]
        # independent ODT target: Σ volₜ·circumcenterₜ / Σ volₜ
        sx=0.0; sy=0.0; sz=0.0; w=0.0
        for t in inc
            a=node(m,m.tets[1,t]); b=node(m,m.tets[2,t]); c=node(m,m.tets[3,t]); d=node(m,m.tets[4,t])
            vol=tet_volume(a,b,c,d); cc=cramer_circumcenter(a,b,c,d)
            sx+=vol*cc[1]; sy+=vol*cc[2]; sz+=vol*cc[3]; w+=vol
        end
        odt_target = (sx/w, sy/w, sz/w)
        # centroid (Lloyd) target, for contrast
        neigh=Set{Int}(); for t in inc, i in 1:4; push!(neigh, Int(m.tets[i,t])); end; delete!(neigh, Int(v))
        cen = (sum(node(m,n)[1] for n in neigh)/length(neigh),
               sum(node(m,n)[2] for n in neigh)/length(neigh),
               sum(node(m,n)[3] for n in neigh)/length(neigh))
        ms = smooth_odt(m; iters=1)
        got = (ms.coords[1,v], ms.coords[2,v], ms.coords[3,v])
        @test _dist(got, node(m,v)) > 0.1                 # the node genuinely moved
        @test collect(got) ≈ collect(odt_target) rtol=1e-11   # ODT update matched to machine precision
        @test _dist(odt_target, cen) > 1e-10              # ODT target ≠ centroid target
        @test _dist(got, cen) > 1e3 * _dist(got, odt_target)  # got is orders closer to ODT than to centroid
        @test validate(ms).ok
        @test mesh_crc(ms).sha==
            "31a280b6a063b428a11772b752ca1d9a64f70a8d4728029c23064a78e977ee00"
        for scale in (1e-100,1e100)
            scaled=Mesh(m.coords.*scale;tets=m.tets)
            scaled_result=smooth_odt(scaled;iters=1)
            scaled_target=ntuple(d->scaled_result.coords[d,v]/scale,3)
            @test collect(scaled_target)≈collect(got) atol=1e-12 rtol=1e-12
        end
    end

    @testset "ODT preserves volume + validity, lifts mean dihedral (slivery cloud)" begin
        r = _RO(0xC0FFEE); n = 300
        xs=Float64[_nfo(r) for _ in 1:n]; ys=Float64[_nfo(r) for _ in 1:n]; zs=Float64[_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        q0 = mesh_quality(m)
        ms = smooth_odt(m; iters=5)
        q1 = mesh_quality(ms)
        @test mvol(ms) ≈ mvol(m) rtol=1e-9                # boundary fixed ⇒ volume preserved
        @test validate(ms).ok                             # no inverted tets
        @test ntets(ms) == ntets(m)                       # topology unchanged
        @test ms.coords != m.coords                       # interior nodes moved
        @test q1.mean_dihedral_deg > q0.mean_dihedral_deg # ODT lifts the mean dihedral
    end

    @testset "ODT pins boundary nodes (cube unchanged)" begin
        C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
        F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
        ct=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); ct[:,k]=Int32[f...]; end
        m = tetrahedralize(Mesh(C; tris=ct))
        ms = smooth_odt(m; iters=5)
        @test ms.coords == m.coords                       # all nodes on the boundary ⇒ fixed
        @test mvol(ms) ≈ 1.0 rtol=1e-6
    end

    @testset "smooth_optimize targets slivers (min-dihedral), preserves volume/validity" begin
        r = _RO(UInt64(12345)); n = 300
        xs = [_nfo(r) for _ in 1:n]; ys = [_nfo(r) for _ in 1:n]; zs = [_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        q0 = mesh_quality(m)
        @test q0.n_slivers > 0                            # the random-cloud input is genuinely slivery
        mo = smooth_optimize(m; iters=10)
        q1 = mesh_quality(mo)
        @test validate(mo).ok                             # no inverted/degenerate tet
        @test mvol(mo) ≈ mvol(m) rtol=1e-9                # boundary fixed ⇒ total volume preserved
        @test ntets(mo) == ntets(m)                       # geometry-only move: topology unchanged
        @test q1.n_slivers < q0.n_slivers                 # strictly fewer slivers (it works)
        # and it beats the mean-optimizing Laplacian on the sliver count (targets the min angle)
        @test q1.n_slivers <= mesh_quality(smooth_laplacian(m; iters=10)).n_slivers
    end

    @testset "remove_slivers: converging exudation driver — fewer slivers, valid, volume-preserving" begin
        r = _RO(UInt64(12345)); n = 300
        xs = [_nfo(r) for _ in 1:n]; ys = [_nfo(r) for _ in 1:n]; zs = [_nfo(r) for _ in 1:n]
        m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=1))
        q0 = mesh_quality(m)
        @test q0.n_slivers > 0
        mo, rep = remove_slivers(m)
        @test validate(mo).ok                                 # never emits an invalid mesh
        @test mvol(mo) ≈ mvol(m) rtol=1e-9                    # boundary-preserving ⇒ volume conserved
        @test ntets(mo) == ntets(m)                           # geometric route: topology unchanged
        @test rep.slivers_after < rep.slivers_before          # it strictly reduces slivers
        @test rep.slivers_after == mesh_quality(mo).n_slivers # report is the measured truth
        @test rep.min_dihedral_after >= rep.min_dihedral_before - 1e-9   # worst angle non-worsening
    end

    @testset "smooth_optimize pins boundary nodes (all-boundary cube unchanged)" begin
        C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
        F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
        ct=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); ct[:,k]=Int32[f...]; end
        m = tetrahedralize(Mesh(C; tris=ct))
        mo = smooth_optimize(m; iters=5)
        @test mo.coords == m.coords                       # every node on the boundary ⇒ nothing moves
    end

    @testset "public documentation and deterministic CRC" begin
        C=Float64[0 1 0 0;0 0 1 0;0 0 0 1]
        m=Mesh(C;tets=reshape(Int32[1,2,3,4],4,1),tet_tag=Int32[7])
        @test mesh_crc(smooth_odt(m;iters=0)).sha==mesh_crc(m).sha
        @test isempty(Base.Docs.undocumented_names(Tessella.Optimize;private=false))
        @test isempty(Test.detect_ambiguities(Tessella.Optimize;recursive=true))
    end
end
