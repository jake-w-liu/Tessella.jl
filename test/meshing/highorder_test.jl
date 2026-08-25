# ── Stage-6 CRC suite: quadratic (P2) tet generation + curving + type-11 I/O ────
#
# Correctness  : exactly one shared mid-node per edge (node count = corners + edges,
#                cross-checked against MeshTypes.unique_edges); corners preserved;
#                edge nodes are exact midpoints; straight P2 volume == linear volume.
# Curving      : only genuine BOUNDARY-surface edges are curved (interior chords are
#                left untouched — regression against the "both corners on surface"
#                bug that tangled the volume); curved nodes land exactly on the true
#                surface; NO element is inverted (checked by an INDEPENDENT P2
#                Jacobian sampler in addition to the exact global guard); the guard
#                reverts a projection that would invert an incident element.
# Robustness   : single tet, filled volume mesh, empty mesh.
# Completeness : gmsh type-11 file writes and reads back to the same linear
#                connectivity (CRC) — solver-consumable.

using Test
import SHA
import Tessella
using Tessella.MeshTypes
using Tessella.Mesh3D
using Tessella.Geometry
using Tessella.Heal
using Tessella.HighOrder
using Tessella.IO

lvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)

# ── INDEPENDENT P2 Jacobian oracle (loop-form gradients; distinct from the module's
#    hand-expanded _p2_grads and its 20-node sample set) ───────────────────────────
function _oracle_grads(r, s, t)
    L = (1-r-s-t, r, s, t)
    dL = ((-1.0,-1.0,-1.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
    g = Vector{NTuple{3,Float64}}(undef, 10)
    for i in 1:4; c = 4L[i]-1; g[i] = (c*dL[i][1], c*dL[i][2], c*dL[i][3]); end
    for (k,(a,b)) in enumerate(((1,2),(2,3),(3,1),(1,4),(3,4),(2,4)))
        g[4+k] = (4*(L[a]*dL[b][1]+L[b]*dL[a][1]), 4*(L[a]*dL[b][2]+L[b]*dL[a][2]), 4*(L[a]*dL[b][3]+L[b]*dL[a][3]))
    end
    return g
end
function _oracle_detJ(p::P2Mesh, t, r, s, u)
    g = _oracle_grads(r, s, u); J = zeros(3,3)
    for k in 1:10
        v = p.tet10[k,t]
        for d in 1:3
            J[d,1]+=g[k][1]*p.coords[d,v]; J[d,2]+=g[k][2]*p.coords[d,v]; J[d,3]+=g[k][3]*p.coords[d,v]
        end
    end
    J[1,1]*(J[2,2]*J[3,3]-J[2,3]*J[3,2]) - J[1,2]*(J[2,1]*J[3,3]-J[2,3]*J[3,1]) + J[1,3]*(J[2,1]*J[3,2]-J[2,2]*J[3,1])
end
# Independent DENSER sample set: the degree-6 barycentric lattice (i,j,k)/6 → 84
# points, distinct from and finer than the module guard's degree-3 (20-node) set,
# so it catches any between-node Jacobian dip the guard could miss. Deterministic.
const _ORACLE_PTS = let S = NTuple{3,Float64}[]
    for j in 0:6, k in 0:6-j, l in 0:6-j-k; push!(S, (j/6, k/6, l/6)); end
    S
end
function _oracle_min_detJ(p::P2Mesh)
    mn = Inf
    for t in 1:ntets(p), (r,s,u) in _ORACLE_PTS; mn = min(mn, _oracle_detJ(p, t, r, s, u)); end
    mn
end
_straight_mid(p, a, b) = ((p.coords[1,a]+p.coords[1,b])/2, (p.coords[2,a]+p.coords[2,b])/2, (p.coords[3,a]+p.coords[3,b])/2)

# mid-node → its two corner endpoints (independent bookkeeping)
function _endpoints(p::P2Mesh)
    slots=((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,3,4),(10,2,4)); e=Dict{Int32,Tuple{Int32,Int32}}()
    for t in 1:ntets(p), (sl,i,j) in slots; e[p.tet10[sl,t]]=(p.tet10[i,t],p.tet10[j,t]); end
    e
end
function _boundary_edge_set(m::Mesh)
    bf,_ = boundary_faces(m.tets); s=Set{Tuple{Int32,Int32}}()
    for f in bf; push!(s,minmax(f[1],f[2]));push!(s,minmax(f[2],f[3]));push!(s,minmax(f[1],f[3])); end
    s
end

# closed sphere of radius R: octahedron subdivided twice, verts snapped to |p|=R.
function _sphere_surface(R)
    octaV=[(R,0.,0.),(-R,0.,0.),(0.,R,0.),(0.,-R,0.),(0.,0.,R),(0.,0.,-R)]
    octaF=[(1,3,5),(1,6,3),(1,5,4),(1,4,6),(2,5,3),(2,3,6),(2,4,5),(2,6,4)]
    snap(q)=(sc=R/sqrt(q[1]^2+q[2]^2+q[3]^2);(q[1]*sc,q[2]*sc,q[3]*sc))
    function subdivide(V,F)
        verts=collect(V); mid=Dict{Tuple{Int,Int},Int}()
        gm(a,b)=get!(mid,(min(a,b),max(a,b))) do
            pa=verts[a];pb=verts[b]; push!(verts,snap(((pa[1]+pb[1])/2,(pa[2]+pb[2])/2,(pa[3]+pb[3])/2))); length(verts)
        end
        nF=Tuple{Int,Int,Int}[]
        for (a,b,c) in F; ab=gm(a,b);bc=gm(b,c);ca=gm(c,a)
            push!(nF,(a,ab,ca));push!(nF,(ab,b,bc));push!(nF,(ca,bc,c));push!(nF,(ab,bc,ca)); end
        verts,nF
    end
    V1,F1=subdivide(octaV,octaF); V2,F2=subdivide(V1,F1)
    C=Matrix{Float64}(undef,3,length(V2)); for (k,v) in enumerate(V2); C[:,k]=[v...]; end
    Tm=Matrix{Int32}(undef,3,length(F2)); for (k,f) in enumerate(F2); Tm[:,k]=Int32[f...]; end
    Mesh(C; tris=Tm)
end

@testset "HighOrder P2 (Stage 6)" begin

    @testset "P2 container and curving contracts" begin
        @test_throws ArgumentError P2Mesh(zeros(3), zeros(Int32,10,0))
        @test_throws ArgumentError P2Mesh(zeros(3,0), zeros(Int32,10))
        @test_throws ArgumentError P2Mesh(zeros(2,10), zeros(Int32,10,0))
        @test_throws ArgumentError P2Mesh(zeros(3,10), zeros(Int32,9,0))
        @test_throws ArgumentError P2Mesh(trues(3,0), zeros(Int32,10,0))
        @test_throws ArgumentError P2Mesh(zeros(3,0), falses(10,0))
        @test_throws ArgumentError P2Mesh(
            zeros(3,0), zeros(Int32,10,0); tet_tag=Bool[])
        @test_throws ArgumentError P2Mesh(
            zeros(3,10), reshape(Int32.(1:10),10,1); tet_tag=Int32[])
        @test_throws ArgumentError P2Mesh(fill(NaN,3,10), reshape(Int32.(1:10),10,1))
        @test_throws ArgumentError P2Mesh(zeros(3,10), reshape(Int32[1,2,3,4,5,6,7,8,9,11],10,1))
        @test_throws ArgumentError P2Mesh(zeros(3,10), reshape(Int32[1,2,3,4,5,6,7,8,9,9],10,1))

        empty = P2Mesh(zeros(3,0), zeros(Int32,10,0))
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.), (0.,0.,1.), -1)
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.), (0.,0.,0.), 1)
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.), (0.,0.,1.), 1; rtol=NaN)
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.), (0.,0.,1.), true)
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.), (0.,0.,1.), "1")
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.,0.), (0.,0.,1.), 1)
        @test_throws ArgumentError curve_to_cylinder!(empty, (0.,0.,0.), (0.,0.,1.,0.), 1)
        @test_throws ArgumentError curve_to_cylinder!(empty, (false,0.,0.), (0.,0.,1.), 1)
        @test_throws ArgumentError curve_to_surface!(empty, identity, (x,y,z)->true; rtol=-1)
        @test_throws ArgumentError curve_to_surface!(empty, identity, (x,y,z)->true; rtol=true)
        @test_throws ArgumentError curve_to_surface!(empty, identity, (x,y,z)->true; rtol="1")

        @test Tessella.P2Mesh === P2Mesh
        @test Tessella.p2_tetmesh === p2_tetmesh
    end

    @testset "single tet → 10 nodes, midpoints, volume preserved" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                 tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[23])
        p = p2_tetmesh(m)
        @test size(p.tet10) == (10, 1)
        @test size(p.coords, 2) == 4 + 6            # 4 corners + 6 edges
        @test p.coords[:,1:4] == m.coords           # corners unchanged
        @test p.coords[:, p.tet10[5,1]] ≈ [0.5,0.0,0.0]   # mid of edge (1,2)
        @test p.tet_tag == Int32[23]
        @test p2_volume(p) ≈ 1/6 rtol=1e-12
        @test p2_min_jacobian(p) ≈ 6*(1/6) rtol=1e-12     # straight ⇒ detJ = 6·vol > 0
        @test validate(p).ok

        direct = P2Mesh(p.coords, p.tet10; tet_tag=p.tet_tag)
        @test !Base.mightalias(direct.coords, p.coords)
        @test !Base.mightalias(direct.tet10, p.tet10)
        @test !Base.mightalias(direct.tet_tag, p.tet_tag)
        @test !Base.mightalias(p.coords, m.coords)
        @test !Base.mightalias(p.tet10, m.tets)
        @test !Base.mightalias(p.tet_tag, m.tet_tag)

        @test_throws ArgumentError p2_tetmesh(m; max_nodes=true)
        @test_throws ArgumentError p2_tetmesh(m; max_tets=false)
        @test_throws ArgumentError p2_tetmesh(m; max_nodes=10.0)
        @test_throws ArgumentError p2_tetmesh(m; max_tets="1")
        @test_throws ArgumentError p2_tetmesh(m; max_nodes=-1)
        @test_throws ArgumentError p2_tetmesh(m; max_tets=-1)
        @test_throws ArgumentError p2_tetmesh(
            m; max_nodes=big(typemax(Int32)) + 1)
        @test_throws ArgumentError p2_tetmesh(
            m; max_tets=big(typemax(Int32)) + 1)
        @test_throws ArgumentError p2_tetmesh(m; max_nodes=9)
        @test_throws ArgumentError p2_tetmesh(m; max_tets=0)
        bounded = p2_tetmesh(m; max_nodes=UInt(10), max_tets=BigInt(1))
        @test (nnodes(bounded), ntets(bounded)) == (10, 1)
    end

    @testset "shared mid-nodes: node count = corners + unique edges" begin
        s = box_surface(0,1,0,1,0,1)
        m = tetrahedralize(s)
        p = p2_tetmesh(m)
        nedges = length(unique_edges(Matrix{Int32}(undef,3,0), m.tets))
        @test size(p.coords, 2) == nnodes(m) + nedges     # one shared node per edge
        @test p2_volume(p) ≈ lvol(m) rtol=1e-9            # straight P2 ⇒ same volume
    end

    @testset "gmsh type-11 round-trip (reads back to linear connectivity)" begin
        m = tetrahedralize(cylinder_surface((0.,0,0),(0.,0,1),1.0,2.0; nθ=16))
        p = p2_tetmesh(m)
        dir = mktempdir(); path = joinpath(dir, "p2.msh")
        write_msh_p2(path, p; tet_tag=fill(Int32(7), ntets(m)))
        f = read_msh(path)                                # type-11 read as 4-corner tets
        @test ntets(f.mesh) == ntets(m)
        @test mesh_crc(f.mesh).sha == mesh_crc(m).sha     # linear connectivity preserved
        @test all(f.mesh.tet_tag .== 7)
        @test validate(f.mesh).ok

        tagged = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                      tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[23])
        tagged_p2 = p2_tetmesh(tagged)
        tagged_path = joinpath(dir, "tagged-p2.msh")
        write_msh_p2(tagged_path, tagged_p2)
        tagged_read = read_msh(tagged_path)
        @test tagged_read.mesh.tet_tag == Int32[23]
        @test mesh_crc(tagged_read.mesh).sha == mesh_crc(tagged).sha
        @test bytes2hex(SHA.sha256(read(tagged_path))) ==
              "5a83ebe0386bda71c6761148ed3fe2f964f16c2da2f0b66b6951ef558f4927ab"
    end

    @testset "curve_to_cylinder!: curved nodes on the true cylinder, NO inversion" begin
        R = 2.0
        m = tetrahedralize(cylinder_surface((0.,0,0),(0.,0,1),R,5.0; nθ=16, nz=3))
        p = p2_tetmesh(m)
        @test p2_min_jacobian(p) > 0                       # straight mesh is valid
        nc = curve_to_cylinder!(p, (0.,0,0), (0.,0,1), R)
        @test nc > 0                                       # boundary-wall edges were curved
        # CORE: no element inverted — verified by BOTH the module guard's node set
        # and an INDEPENDENT, finer sampler (distinct gradients + degree-6 lattice).
        @test p2_min_jacobian(p) > 0
        @test _oracle_min_detJ(p) > 0
        # Each boundary-wall mid is EITHER curved exactly onto the cylinder OR left on
        # its straight chord (the validity guard reverts a curve that would invert an
        # incident element — a legitimate outcome). Count the curved ones; it is `nc`.
        e = _endpoints(p); be = _boundary_edge_set(m)
        rad(v) = sqrt(p.coords[1,v]^2 + p.coords[2,v]^2)
        onwall(v) = abs(rad(v) - R) <= 1e-6*R              # tol > SoS perturbation (~5e-8)
        curved = 0
        for (mid,(a,b)) in e
            (onwall(a) && onwall(b) && minmax(a,b) in be) || continue
            on_cyl = abs(rad(mid) - R) < 1e-9
            straight = (p.coords[1,mid],p.coords[2,mid],p.coords[3,mid]) == _straight_mid(p,a,b)
            @test on_cyl || straight                       # curved onto the wall, or reverted
            on_cyl && (curved += 1)
        end
        # every displaced mid lands on the cylinder, so the on-wall count ≥ nc (it
        # also includes vertical wall edges whose straight midpoint is already at R).
        @test curved >= nc && nc > 0
    end

    @testset "curve_to_surface! (sphere): interior chords untouched, boundary curved" begin
        R = 1.7
        sphere = _sphere_surface(R)
        @test is_meshable(sphere)[1]                       # watertight, oriented, manifold
        # Exercise the interior-chord regression on the raw restricted Delaunay.
        # The certified public fill may repair this faceted sphere with a one-centre
        # boundary fan, which correctly has no surface-to-surface interior chords.
        m = tetrahedralize(sphere; check=false)
        @test validate(m).ok && is_closed_manifold(m)
        p = p2_tetmesh(m)
        @test p2_min_jacobian(p) > 0                       # straight mesh valid

        proj(x,y,z)   = (s = R/sqrt(x^2+y^2+z^2); (x*s, y*s, z*s))
        onsurf(x,y,z) = abs(sqrt(x^2+y^2+z^2) - R) <= 1e-6*R
        nc = curve_to_surface!(p, proj, onsurf)
        @test nc > 0

        # CORE fix: NO element inverted (was 127/166 with the "both corners" bug),
        # by the independent finer oracle AND the module's own metric.
        @test _oracle_min_detJ(p) > 0
        @test p2_min_jacobian(p) > 0

        # Partition mid-nodes whose BOTH corner endpoints are on the sphere into
        # boundary-surface edges vs interior chords, using an INDEPENDENT recompute.
        e = _endpoints(p); be = _boundary_edge_set(m)
        rad(v) = sqrt(p.coords[1,v]^2 + p.coords[2,v]^2 + p.coords[3,v]^2)
        onS(v) = abs(rad(v) - R) <= 1e-6*R
        curved_bnd = 0; interior_seen = 0
        for (mid,(a,b)) in e
            (onS(a) && onS(b)) || continue
            here = (p.coords[1,mid],p.coords[2,mid],p.coords[3,mid])
            if minmax(a,b) in be
                # boundary edge: curved exactly onto the sphere, or reverted (left straight)
                @test abs(rad(mid) - R) < 1e-9 || here == _straight_mid(p,a,b)
                abs(rad(mid) - R) < 1e-9 && (curved_bnd += 1)
            else
                interior_seen += 1
                # REGRESSION: an interior chord's mid is byte-identical to its straight
                # midpoint — NEVER moved (the old "both corners on surface" code snapped
                # every such mid out onto the sphere, tangling 127/166 elements).
                @test here == _straight_mid(p,a,b)
            end
        end
        @test interior_seen > 0                            # interior chords do exist here
        @test curved_bnd == nc && curved_bnd > 0           # all curves are boundary-surface edges
    end

    @testset "validity guard reverts an inverting projection" begin
        # A single reference tet: all 4 faces are boundary ⇒ all 6 edges qualify.
        # `onedge12` marks only the two corners of edge (1,2) (they lie on the x-axis);
        # `far` would drag that edge's mid out to (5,0,0), folding the element.
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]; tets=reshape(Int32[1,2,3,4],4,1))
        onedge12(x,y,z) = abs(y) < 1e-12 && abs(z) < 1e-12
        far(x,y,z) = (5.0, 0.0, 0.0)
        # self-validating: placing the mid at (5,0,0) really does invert the element.
        bad = p2_tetmesh(m); bad.coords[:, bad.tet10[5,1]] = [5.0, 0.0, 0.0]
        @test p2_min_jacobian(bad) < 0                         # the projection WOULD invert
        # the guard must therefore refuse it, leaving the mesh straight and valid.
        p = p2_tetmesh(m)
        nc = curve_to_surface!(p, far, onedge12)
        @test nc == 0                                          # only qualifying move is reverted
        @test p.coords[:, p.tet10[5,1]] == [0.5, 0.0, 0.0]     # edge-(1,2) mid unchanged
        @test p2_min_jacobian(p) > 0                           # still valid (≈ straight)
    end

    @testset "curving is transactional across callback failures" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                 tets=reshape(Int32[1,2,3,4],4,1))
        p = p2_tetmesh(m)
        pristine = copy(p.coords)
        calls = Ref(0)
        function failing_project(x, y, z)
            calls[] += 1
            calls[] == 2 && error("deliberate projection failure")
            return (x, y, z + 0.01)
        end
        @test_throws ErrorException curve_to_surface!(
            p, failing_project, (x,y,z)->true)
        @test calls[] == 2
        @test p.coords == pristine
        @test validate(p).ok

        @test_throws ArgumentError curve_to_surface!(
            p, (x,y,z)->(x,y), (x,y,z)->true)
        @test p.coords == pristine
        @test_throws ArgumentError curve_to_surface!(
            p, (x,y,z)->(true,y,z), (x,y,z)->true)
        @test p.coords == pristine
        @test_throws ArgumentError curve_to_surface!(
            p, (x,y,z)->(x,y,z), (x,y,z)->1)
        @test p.coords == pristine

        @test curve_to_surface!(
            p, (x,y,z)->(Inf,y,z), (x,y,z)->true) == 0
        @test p.coords == pristine
    end

    @testset "p2_volume integrates curved geometry, not only corner tets" begin
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                 tets=reshape(Int32[1,2,3,4],4,1))
        p = p2_tetmesh(m)
        linear = p2_volume(p)
        # Move one mid-edge node without folding the element. The corner tet is
        # unchanged, but the isoparametric P2 volume must respond to the curvature.
        p.coords[3,p.tet10[5,1]] += 0.05
        @test p2_min_jacobian(p) > 0
        @test p2_volume(p) != linear
        # Independent degree-3 quadrature oracle for the exact cubic determinant.
        vals = [_oracle_detJ(p,1,0.25,0.25,0.25),
                _oracle_detJ(p,1,1/6,1/6,1/6),
                _oracle_detJ(p,1,0.5,1/6,1/6),
                _oracle_detJ(p,1,1/6,0.5,1/6),
                _oracle_detJ(p,1,1/6,1/6,0.5)]
        oracle = abs(((-4/5)*vals[1] + (9/20)*sum(vals[2:5])) / 6)
        @test p2_volume(p) ≈ oracle rtol=1e-13
    end

    @testset "global Jacobian certificate catches a between-sample fold" begin
        # This deterministic P2 tet is positive at all 20 nodes used by the old
        # degree-3 lattice guard (minimum 0.0843), but det(J)=-1.197 at the interior
        # point (1/4,1/2,1/4). A sampled guard silently accepted the fold.
        m = Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                 tets=reshape(Int32[1,2,3,4],4,1))
        p = p2_tetmesh(m)
        p.coords[:,p.tet10[5:10,1]] = Float64[
             0.598244   1.71802  -0.674877  -0.150606  -1.12732    0.694353;
             0.152545   0.407208  1.05877   -0.0461052  0.525971   0.280946;
             0.0303981 -1.00104   0.0559645  0.276108   0.689428   0.981785]
        old_samples = [(j/3,k/3,l/3) for j in 0:3 for k in 0:3-j for l in 0:3-j-k]
        @test minimum(_oracle_detJ(p,1,r,s,t) for (r,s,t) in old_samples) > 0.08
        @test _oracle_detJ(p,1,0.25,0.5,0.25) < -1.19
        @test p2_min_jacobian(p) < 0
        @test_throws ArgumentError p2_volume(p)
    end

    @testset "midpoint and P2 writer contracts" begin
        huge = Mesh(Float64[1e308 1e308; 0 1; 0 0];
                    segs=reshape(Int32[1,2],2,1))
        @test p2_tetmesh(huge).coords == huge.coords # no tets: no midpoint needed
        collapsed = Mesh(Float64[1e308 nextfloat(1e308) 0 0;
                                  0 0 1 0; 0 0 0 1];
                         tets=reshape(Int32[1,2,3,4],4,1))
        @test_throws ArgumentError p2_tetmesh(collapsed)

        tiny = nextfloat(0.0)
        subnormal = Mesh(Float64[tiny 5tiny tiny tiny;
                                  0 0 1 0; 0 0 0 1];
                         tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[9])
        subnormal_p2 = p2_tetmesh(subnormal)
        @test subnormal_p2.coords[:, subnormal_p2.tet10[5,1]] ==
              Float64[3tiny, 0, 0]
        @test p2_min_jacobian(subnormal_p2) == 4tiny
        @test p2_volume(subnormal_p2) == tiny
        @test subnormal_p2.tet_tag == Int32[9]
        @test validate(subnormal_p2).ok

        p = p2_tetmesh(Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                            tets=reshape(Int32[1,2,3,4],4,1)))
        @test_throws ArgumentError P2Mesh(p.coords, p.tet10; tet_tag=Int32[-1])
        @test_throws ArgumentError P2Mesh(p.coords, p.tet10; tet_tag=Bool[true])
        @test_throws ArgumentError write_msh_p2("ignored.msh",p;tet_tag=Int64[typemax(Int64)])
        @test_throws ArgumentError write_msh_p2("ignored.msh",p;tet_tag=Int32[-1])
        @test_throws ArgumentError write_msh_p2("ignored.msh",p;tet_tag=Bool[true])
        @test_throws ArgumentError write_msh_p2("ignored.msh",p;tet_tag=[1.0])
        @test_throws ArgumentError write_msh_p2("",p)
        @test_throws ArgumentError write_msh_p2(nothing,p)

        atomic_dir = mktempdir()
        atomic_path = joinpath(atomic_dir, "atomic.msh")
        write(atomic_path, "preserve-existing-output")
        folded = P2Mesh(p.coords, p.tet10; tet_tag=p.tet_tag)
        folded.coords[:, folded.tet10[5,1]] = [5.0, 0.0, 0.0]
        @test !validate(folded).ok
        @test_throws ArgumentError write_msh_p2(atomic_path, folded)
        @test read(atomic_path, String) == "preserve-existing-output"
        @test_throws ArgumentError write_msh_p2(atomic_dir, p)

        two=Mesh(Float64[0 1 0 0 0;0 0 1 0 0;0 0 0 1 -1];
                 tets=Int32[1 1;2 3;3 2;4 5])
        q=p2_tetmesh(two);badT=copy(q.tet10);badC=hcat(q.coords,q.coords[:,badT[5,1]])
        edge_slots=((5,1,2),(6,2,3),(7,3,1),(8,1,4),(9,3,4),(10,2,4))
        slot=only(s for (s,i,j) in edge_slots
                    if minmax(badT[i,2],badT[j,2])==(Int32(1),Int32(2)))
        badT[slot,2]=Int32(size(badC,2))
        @test_throws ArgumentError P2Mesh(badC,badT)
    end

    @testset "scale, translation, and mutable-storage validation" begin
        function curved_reference(scale)
            m = Mesh(Float64[scale 0 0 0;
                             0 scale 0 0;
                             0 0 0 scale];
                     tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[17])
            p = p2_tetmesh(m)
            @test curve_to_cylinder!(
                p, (0.,0.,0.), (0.,0.,1.), scale; rtol=1e-12) == 1
            return p
        end
        reference = curved_reference(1.0)
        for scale in (1e-100, 1.0, 1e100)
            p = curved_reference(scale)
            @test p.tet10 == reference.tet10
            @test p.tet_tag == reference.tet_tag
            @test p.coords ./ scale ≈ reference.coords atol=0 rtol=2eps()
            @test p2_min_jacobian(p) / scale^3 ≈
                  p2_min_jacobian(reference) atol=0 rtol=8eps()
            @test validate(p).ok
        end

        offset = 1e100
        width = 16eps(offset)
        translated = Mesh(Float64[offset + width offset offset offset;
                                   offset offset + width offset offset;
                                   offset offset offset offset + width];
                          tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[17])
        translated_p2 = p2_tetmesh(translated)
        @test curve_to_cylinder!(translated_p2, (offset,offset,offset),
                                 (0.,0.,1.), width; rtol=1e-12) == 1
        @test translated_p2.tet10 == reference.tet10
        @test validate(translated_p2).ok

        huge = floatmax(Float64)
        narrow = 1e-154
        extreme = Mesh(Float64[-huge huge -huge -huge;
                                0 0 narrow 0;
                                0 0 0 narrow];
                       tets=reshape(Int32[1,2,3,4],4,1))
        extreme_p2 = p2_tetmesh(extreme)
        @test curve_to_surface!(extreme_p2, (x,y,z)->(x,y,z),
                                (x,y,z)->false; rtol=1e-308) == 0
        @test validate(extreme_p2).ok

        function fresh()
            p2_tetmesh(Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                            tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[4]))
        end
        corrupt_coordinate = fresh()
        corrupt_coordinate.coords[1,1] = NaN
        @test !validate(corrupt_coordinate).ok
        @test_throws ArgumentError p2_volume(corrupt_coordinate)
        @test_throws ArgumentError p2_min_jacobian(corrupt_coordinate)

        corrupt_connectivity = fresh()
        corrupt_connectivity.tet10[1,1] = 0
        @test !validate(corrupt_connectivity).ok
        callback_calls = Ref(0)
        @test_throws ArgumentError curve_to_surface!(
            corrupt_connectivity,
            (x,y,z)->(callback_calls[] += 1; (x,y,z)),
            (x,y,z)->true)
        @test callback_calls[] == 0

        corrupt_tag = fresh()
        corrupt_tag.tet_tag[1] = -1
        @test !validate(corrupt_tag).ok
        corrupt_path = joinpath(mktempdir(), "corrupt.msh")
        @test_throws ArgumentError write_msh_p2(corrupt_path, corrupt_tag)
        @test !ispath(corrupt_path)
    end

    @testset "empty mesh" begin
        e = Mesh(Matrix{Float64}(undef,3,0))
        p = p2_tetmesh(e)
        @test size(p.tet10, 2) == 0
        @test p2_volume(p) == 0.0
        @test p2_min_jacobian(p) == 0.0
        @test curve_to_cylinder!(p, (0.,0,0), (0.,0,1), 1.0) == 0
        @test curve_to_surface!(p, (x,y,z)->(x,y,z), (x,y,z)->true) == 0
    end

    @testset "public documentation" begin
        @test isempty(Base.Docs.undocumented_names(Tessella.HighOrder; private=false))
        @test isempty(Test.detect_ambiguities(Tessella.HighOrder; recursive=true))
    end
end
