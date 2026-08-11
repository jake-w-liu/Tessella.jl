# ── Stage-3 CRC suite: 3-D Delaunay kernel ──────────────────────────────────────
#
# Correctness  : exact empty-circumsphere (is_delaunay3 via insphere_sos) is the
#                defining oracle; convex-region volume conservation (Σ tet vol =
#                bounding-box volume for a box-filling point set); Euler χ=1 (ball);
#                boundary faces form a closed manifold (2-sphere), χ=2.
# Robustness   : random clouds (multiple seeds), the maximally-degenerate cube
#                (8 cospherical corners) resolved by symbolic perturbation, and
#                a structured box grid.
# Completeness : exported tet meshes validate() (positive volumes, manifold).

using Test
using Tessella
using Tessella.Mesh3D
using Tessella.MeshTypes
using Tessella.Geometry
using Tessella.Heal: is_meshable

mesh_vol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                  for t in 1:ntets(m); init=0.0)

mutable struct _R3; s::UInt64; end
_nf(r::_R3) = (r.s ⊻= r.s<<13; r.s ⊻= r.s>>7; r.s ⊻= r.s<<17; (r.s>>11)/Float64(2^53))

@testset "Mesh3D Delaunay (Stage 3)" begin

    @testset "random clouds: exact empty-circumsphere + valid + manifold" begin
        for seed in (1, 7, 42)
            r = _R3(UInt64(seed)*0x9E3779B97F4A7C15 + 1)
            n = 90
            xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]; zs=Float64[_nf(r) for _ in 1:n]
            T = delaunay3d(xs, ys, zs; rng_seed=seed)
            @test check_consistency3(T)[1]
            dok, nv = is_delaunay3(T)
            @test dok
            @test nv == 0
            m = to_mesh3(T)
            @test validate(m).ok                       # positive volumes, manifold
            @test euler_characteristic(m) == 1         # tetrahedralized ball
            _, maxinc = boundary_faces(m.tets)
            @test maxinc == 2                          # closed manifold boundary
            @test boundary_euler(m) == 2               # boundary is a 2-sphere
        end
    end

    @testset "convex box: Σ tet volume = box volume" begin
        # a grid of points filling [0,2]×[0,3]×[0,1]; convex hull = the box.
        xs=Float64[]; ys=Float64[]; zs=Float64[]
        for i in 0:3, j in 0:4, k in 0:2
            push!(xs, 2*i/3); push!(ys, 3*j/4); push!(zs, 1*k/2)
        end
        T = delaunay3d(xs, ys, zs; rng_seed=3)
        @test check_consistency3(T)[1]
        m = to_mesh3(T)
        @test validate(m).ok
        @test mesh_vol(m) ≈ 6.0 rtol=1e-6              # 2·3·1, perturbation-tolerant
    end

    @testset "degenerate unit cube (8 cospherical corners) resolved" begin
        cx=Float64[0,1,1,0,0,1,1,0]; cy=Float64[0,0,1,1,0,0,1,1]; cz=Float64[0,0,0,0,1,1,1,1]
        T = delaunay3d(cx, cy, cz)
        @test check_consistency3(T)[1]
        @test is_delaunay3(T)[1]
        m = to_mesh3(T)
        @test validate(m).ok                           # NO flat tets
        @test mesh_vol(m) ≈ 1.0 rtol=1e-6
        @test euler_characteristic(m) == 1
    end

    @testset "single tet (4 points)" begin
        T = delaunay3d([0.0,1,0,0],[0.0,0,1,0],[0.0,0,0,1])
        m = to_mesh3(T)
        @test ntets(m) == 1
        @test mesh_vol(m) ≈ 1/6 rtol=1e-6
        @test validate(m).ok
    end

    @testset "seed-independence (deterministic perturbation ⇒ same mesh)" begin
        r = _R3(0xDEAD_BEEF_1234_5678)
        n = 80
        xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]; zs=Float64[_nf(r) for _ in 1:n]
        shas = String[]
        for seed in (1, 5, 99)
            m = to_mesh3(delaunay3d(xs, ys, zs; rng_seed=seed))
            push!(shas, mesh_crc(m).sha)
        end
        @test all(==(shas[1]), shas)
    end

    @testset "tetrahedralize: fill a domain from its boundary surface" begin
        # closed cube surface (12 outward triangles) → filled volume 1
        C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
        F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
        ct=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); ct[:,k]=Int32[f...]; end
        cs=Mesh(C; tris=ct)
        @test boundary_edges(cs.tris)[1] |> isempty       # closed surface
        m=tetrahedralize(cs)
        @test mesh_vol(m) ≈ 1.0 rtol=1e-6
        @test validate(m).ok
        _, mi = boundary_faces(m.tets); @test mi == 2      # watertight fill

        # non-convex L-prism (2×2×1 minus a 1×1×1 corner) → volume 3
        base=[(0.0,0.0),(2.0,0.0),(2.0,1.0),(1.0,1.0),(1.0,2.0),(0.0,2.0)]; nb=length(base)
        LC=Matrix{Float64}(undef,3,2nb)
        for (i,(x,y)) in enumerate(base); LC[:,i]=[x,y,0.0]; LC[:,i+nb]=[x,y,1.0]; end
        bt=[(1,2,3),(1,3,4),(1,4,5),(1,5,6)]
        lt=NTuple{3,Int32}[]
        for (a,b,c) in bt; push!(lt,(Int32(a),Int32(c),Int32(b))); end
        for (a,b,c) in bt; push!(lt,(Int32(a+nb),Int32(b+nb),Int32(c+nb))); end
        for i in 1:nb; j=i%nb+1; push!(lt,(Int32(i),Int32(j),Int32(j+nb))); push!(lt,(Int32(i),Int32(j+nb),Int32(i+nb))); end
        ltm=Matrix{Int32}(undef,3,length(lt)); for (k,f) in enumerate(lt); ltm[:,k]=Int32[f...]; end
        ls=Mesh(LC; tris=ltm)
        ml=tetrahedralize(ls)
        @test mesh_vol(ml) ≈ 3.0 rtol=1e-6                 # concave corner excluded
        @test validate(ml).ok

        # convex octahedron (vertices ±1 on axes) → volume 4/3
        OC=Float64[1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1]
        OF=[(1,3,5),(3,2,5),(2,4,5),(4,1,5),(3,1,6),(2,3,6),(4,2,6),(1,4,6)]
        ot=Matrix{Int32}(undef,3,length(OF)); for (k,f) in enumerate(OF); ot[:,k]=Int32[f...]; end
        os=Mesh(OC; tris=ot)
        mo=tetrahedralize(os)
        @test mesh_vol(mo) ≈ 4/3 rtol=1e-6
        @test validate(mo).ok
    end

    @testset "genus-1, thin features, multi-region (coax junction)" begin
        # axis-z cylinder volume surface, N-gon rims, K axial levels
        function cyl_vol(cx,cy,cz,R,H,N,K)
            V=Tuple{Float64,Float64,Float64}[]
            for j in 0:K-1, i in 0:N-1; push!(V,(cx+R*cospi(2i/N),cy+R*sinpi(2i/N),cz+H*j/(K-1))); end
            ci=length(V)+1; push!(V,(cx,cy,cz)); cti=length(V)+1; push!(V,(cx,cy,cz+H))
            idx(j,i)=(j-1)*N+mod(i,N)+1; Tr=NTuple{3,Int32}[]
            for j in 1:K-1, i in 0:N-1
                a=idx(j,i);b=idx(j,i+1);c=idx(j+1,i+1);d=idx(j+1,i)
                push!(Tr,(Int32(a),Int32(b),Int32(c))); push!(Tr,(Int32(a),Int32(c),Int32(d)))
            end
            for i in 0:N-1; push!(Tr,(Int32(ci),Int32(idx(1,i+1)),Int32(idx(1,i)))); end
            for i in 0:N-1; push!(Tr,(Int32(cti),Int32(idx(K,i)),Int32(idx(K,i+1)))); end
            C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
            tm=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); tm[:,k]=Int32[f...]; end
            Mesh(C; tris=tm)
        end
        function box_surf(x0,x1,y0,y1,z0,z1)
            C=Float64[x0 x1 x1 x0 x0 x1 x1 x0; y0 y0 y1 y1 y0 y0 y1 y1; z0 z0 z0 z0 z1 z1 z1 z1]
            F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
            t=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end; Mesh(C; tris=t)
        end
        function box_tun(ox0,ox1,oy0,oy1,z0,z1,ix0,ix1,iy0,iy1)
            V=Tuple{Float64,Float64,Float64}[]
            for (x,y) in [(ox0,oy0),(ox1,oy0),(ox1,oy1),(ox0,oy1)]; push!(V,(x,y,z0)); end
            for (x,y) in [(ox0,oy0),(ox1,oy0),(ox1,oy1),(ox0,oy1)]; push!(V,(x,y,z1)); end
            for (x,y) in [(ix0,iy0),(ix1,iy0),(ix1,iy1),(ix0,iy1)]; push!(V,(x,y,z0)); end
            for (x,y) in [(ix0,iy0),(ix1,iy0),(ix1,iy1),(ix0,iy1)]; push!(V,(x,y,z1)); end
            Tr=NTuple{3,Int32}[]; q(a,b,c,d)=(push!(Tr,(Int32(a),Int32(b),Int32(c)));push!(Tr,(Int32(a),Int32(c),Int32(d))))
            q(1,9,10,2);q(2,10,11,3);q(3,11,12,4);q(4,12,9,1); q(5,6,14,13);q(6,7,15,14);q(7,8,16,15);q(8,5,13,16)
            q(1,2,6,5);q(2,3,7,6);q(3,4,8,7);q(4,1,5,8); q(9,13,14,10);q(10,14,15,11);q(11,15,16,12);q(12,16,13,9)
            C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
            tm=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); tm[:,k]=Int32[f...]; end; Mesh(C; tris=tm)
        end

        # thin cylindrical pin (aspect 20) fills to the exact 12-gon prism volume
        pin = cyl_vol(0.0,0.0,0.0, 0.05, 1.0, 12, 4)
        mp = tetrahedralize(pin)
        @test validate(mp).ok
        @test mesh_vol(mp) ≈ 0.5*12*0.05^2*sinpi(2/12)*1.0 rtol=1e-6

        # box with a rectangular through-tunnel (genus-1 bore) → exact volume 24
        case = box_tun(1,5,1,5,1,3, 2,4,2,4)
        mc = tetrahedralize(case)
        @test validate(mc).ok
        @test mesh_vol(mc) ≈ 24.0 rtol=1e-6
        @test boundary_faces(mc.tets)[2] == 2                # watertight fill

        # THE STANDING CASE (representative): 3-region coax junction where the pin,
        # metal case and air gap meet at the bore walls — gmsh 4.13/4.15 leave every
        # volume empty here. Tessella fills all three, each with a positive tet count.
        pinb = box_surf(2.4,3.6, 2.4,3.6, 0,4)
        gap  = box_tun(2,4, 2,4, 1,3, 2.4,3.6,2.4,3.6)
        mm = tetrahedralize_multi([pinb, case, gap])
        tpr = tets_per_region(mm)
        @test length(tpr) == 3
        @test all(v > 0 for v in values(tpr))                # NO empty volume
        @test mesh_vol(mm) ≈ (1.2^2*4 + 24.0 + (2.0^2*2 - 1.2^2*2)) rtol=1e-6  # 5.76+24+5.12
        @test validate(mm).ok
    end

    @testset "ENC-COAX at real scale: all three volumes filled (gmsh leaves empty)" begin
        # dimensions from test/fixtures/enclosure_coax_junction.geo (metres)
        air  = box_surface(-0.0405, 0.2810, -0.0405, 0.2010, -0.0405, 0.3810)   # radiation/air box
        case = box_surface(-0.0005, 0.2205, -0.0005, 0.1405, -0.0005, 0.3005)   # case outer shell
        # coax pin: centre (0.17,0.1605,0.15), axis −y, length 0.1589, radius 8e-4 (aspect ≈ 199)
        pin  = cylinder_surface((0.17,0.1605,0.15), (0.0,-1.0,0.0), 0.0008, 0.1589; nθ=12, nz=40)
        for s in (air, case, pin); @test is_meshable(s)[1]; end
        @test mesh_vol(tetrahedralize(air))  ≈ 0.3215*0.2415*0.4215 rtol=1e-6     # exact box
        @test mesh_vol(tetrahedralize(case)) ≈ 0.221*0.141*0.301   rtol=1e-6      # exact box
        mpin = tetrahedralize(pin)
        @test validate(mpin).ok                                                   # thin pin fills
        @test mesh_vol(mpin) > 0
        # the standing acceptance: all three regions non-empty & validated together
        mm = tetrahedralize_multi([air, case, pin])
        tpr = tets_per_region(mm)
        @test length(tpr) == 3 && all(v > 0 for v in values(tpr))                 # NO empty volume
        @test validate(mm).ok
    end

    @testset "tetrahedralize_conforming: one partition, interfaces SHARED across regions" begin
        # Independent conformance oracle: map each triangular face to the tet-tags
        # incident to it. A conforming partition has (a) no face in >2 tets (manifold),
        # and (b) interface faces incident to two DIFFERENT region tags.
        facetags(m) = begin
            ft = Dict{NTuple{3,Int32},Vector{Int32}}()
            for t in 1:ntets(m), k in 1:4
                vs = (m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
                f = Tuple(sort(Int32[vs[j] for j in 1:4 if j != k]))
                push!(get!(ft, f, Int32[]), m.tet_tag[t])
            end
            ft
        end

        # two unit boxes sharing the plane x=1 → ONE mesh; the interface faces are shared
        # between region 1 and region 2 (with perturb-jitter they would be duplicated/off-plane).
        mm = tetrahedralize_conforming([box_surface(0,1,0,1,0,1), box_surface(1,2,0,1,0,1)])
        @test validate(mm).ok
        @test mesh_vol(mm) ≈ 2.0 rtol=1e-9                       # EXACT — no perturbation drift
        tprc = tets_per_region(mm)
        @test length(tprc) == 2 && all(v > 0 for v in values(tprc))
        ft = facetags(mm)
        @test all(length(v) <= 2 for v in values(ft))           # manifold: no face in >2 tets
        @test count(v -> length(v) == 2 && v[1] != v[2], values(ft)) > 0   # ≥1 shared interface face

        # ENCLOSURE air/case pattern: inner solid box (region 1) inside a surrounding shell
        # (region 2). Meshing all vertices together makes the inner box's 12 faces SHARED
        # between the two regions — a conforming interface, the exact thing coordinate jitter
        # (and hence gmsh's failure mode) destroys.
        inner = box_surface(1,2,1,2,1,2)
        shell = box_shell_surface(0,3,0,3,0,3, 1,2,1,2,1,2)
        me = tetrahedralize_conforming([inner, shell])
        @test validate(me).ok
        @test mesh_vol(me) ≈ 27.0 rtol=1e-9                      # inner(1) + shell(26) = outer(27)
        tpe = tets_per_region(me)
        @test length(tpe) == 2 && all(v > 0 for v in values(tpe))
        fe = facetags(me)
        @test all(length(v) <= 2 for v in values(fe))           # manifold
        @test count(v -> length(v) == 2 && v[1] != v[2], values(fe)) == 12   # the inner box's 12 faces, shared

        # THE ENC-COAX TOPOLOGY, conforming: a coax PIN (cylinder) inside an AIR cavity
        # inside a metal CASE shell. Classified innermost→outermost. The pin's CURVED
        # (N-gon) lateral+cap faces are recovered as shared interface faces too — the
        # exact-coordinate Delaunay conforms to the cylinder without any Steiner points.
        outer  = box_surface(0,6,0,6,0,6)                        # case outer (vol 216)
        cavity = box_surface(1,5,1,5,1,5)                        # air+pin live here (vol 64)
        pin    = cylinder_surface((3.,3.,1.5), (0.,0.,1.), 0.7, 3.0; nθ=8, nz=2)  # fully inside the cavity
        mc = tetrahedralize_conforming([pin, cavity, outer])     # order = pin, air, case
        @test validate(mc).ok
        @test mesh_vol(mc) ≈ 216.0 rtol=1e-9                     # pin + air + case = outer box, EXACT
        tpc = tets_per_region(mc)
        @test length(tpc) == 3 && all(v > 0 for v in values(tpc))   # every region filled (pin/air/case)
        fc = facetags(mc)
        @test all(length(v) <= 2 for v in values(fc))           # manifold everywhere
        # the pin/air interface — the cylinder's 4·nθ surface faces (nθ wall quads×2 + 2 caps×nθ),
        # since the pin sits fully inside the air — are EACH shared between region 1 (pin) and
        # region 2 (air): the CURVED interface is recovered conformingly, no Steiner points.
        @test count(v -> sort(v) == Int32[1,2], values(fc)) == 4*8       # all 32 pin faces conform to air
        @test count(v -> sort(v) == Int32[1,3], values(fc)) == 0         # pin does not touch the case
        @test count(v -> sort(v) == Int32[2,3], values(fc)) > 0          # air/case cavity interface conforms too
    end

    @testset "2-3 / 3-2 flips: consistent, volume-preserving, round-trip identity" begin
        M3 = Tessella.Mesh3D
        totvol(T) = begin
            s=0.0
            for t in eachindex(T.alive)
                (T.alive[t] && !M3._is_ghost_tet(T,t)) || continue
                a=M3._pt(T,M3._vert(T,t,1)); b=M3._pt(T,M3._vert(T,t,2))
                c=M3._pt(T,M3._vert(T,t,3)); d=M3._pt(T,M3._vert(T,t,4))
                s += abs(M3._dot(M3._subn(b,a), M3._cross(M3._subn(c,a), M3._subn(d,a)))/6)
            end; s
        end
        r = _R3(0x999); n = 40
        xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]; zs=Float64[_nf(r) for _ in 1:n]
        T = delaunay3d(xs, ys, zs; rng_seed=1)
        v0 = totvol(T); nt0 = ntets_live(T)
        # find a flippable interior face
        flipped = false
        for t1 in 1:length(T.alive)
            (T.alive[t1] && !M3._is_ghost_tet(T,t1)) || continue
            for k1 in 1:4
                t2 = M3._nbr(T,t1,k1)
                (t2 != 0 && !M3._is_ghost_tet(T,t2)) || continue
                d = M3._vert(T,t1,k1); e = M3._vert(T,t2, M3._nslot(T,t2,t1))
                if flip23!(T, t1, k1)
                    @test check_consistency3(T)[1]
                    @test ntets_live(T) == nt0 + 1                     # 2 → 3
                    @test totvol(T) ≈ v0 rtol=1e-9                     # volume preserved
                    ring = tets_around_edge(T, d, e)
                    @test length(ring) == 3
                    @test flip32!(T, d, e, ring)                       # reverse it
                    @test check_consistency3(T)[1]
                    @test ntets_live(T) == nt0                         # back to 2 → identity
                    @test totvol(T) ≈ v0 rtol=1e-9
                    flipped = true
                    break
                end
            end
            flipped && break
        end
        @test flipped                                                  # a flippable face existed
    end

    @testset "optimize_flips! is safe and reduces slivers (2-3 + 3-2 flips)" begin
        r = _R3(0xABCDEF); n = 300
        xs=Float64[_nf(r) for _ in 1:n]; ys=Float64[_nf(r) for _ in 1:n]; zs=Float64[_nf(r) for _ in 1:n]
        T = delaunay3d(xs, ys, zs; rng_seed=1)
        m0 = to_mesh3(T)
        v0 = mesh_vol(m0); q0 = mesh_quality(m0)
        nf = optimize_flips!(T; passes=6)
        @test nf > 0
        @test check_consistency3(T)[1]                                 # still a valid complex
        m1 = to_mesh3(T)
        @test mesh_vol(m1) ≈ v0 rtol=1e-9                              # volume preserved (flips)
        @test validate(m1).ok
        q1 = mesh_quality(m1)
        @test q1.min_dihedral_deg >= q0.min_dihedral_deg - 1e-9        # min never worsens (safe)
        @test q1.n_slivers < q0.n_slivers                              # 3-2 flips collapse slivers
        @test q1.mean_dihedral_deg > q0.mean_dihedral_deg              # mean quality up
    end

    @testset "degenerate input handled gracefully" begin
        # exact-degenerate (perturb=false): no non-coplanar 4-tuple ⇒ empty mesh
        @test ntets(to_mesh3(delaunay3d([0.0,1,2],[0.0,1,2],[0.0,1,2]; perturb=false))) == 0     # collinear
        @test ntets(to_mesh3(delaunay3d([0.0,1,0,1],[0.0,0,1,1],[0.0,0,0,0]; perturb=false))) == 0 # coplanar
        @test ntets(to_mesh3(delaunay3d([0.0],[0.0],[0.0]))) == 0                                 # 1 point
    end
end
