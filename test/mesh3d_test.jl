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
using Tessella.ExactMesh3D
using Tessella.MeshTypes
using Tessella.Geometry
using Tessella.Heal: is_meshable
using Tessella.Mesh2D: constrained_delaunay, to_mesh, classify_interior

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

    @testset "ENC-COAX LITERAL: parsed from the .geo, native reconstruction, conforming" begin
        # The enclosure .geo is fully parametric — every solid is a gmsh Box/Cylinder.
        # Parse those LITERAL primitives directly (no OCC eval) and reconstruct + mesh the
        # three main volumes at the exact fixture dimensions, as ONE conforming partition
        # (tetrahedralize_conforming_exact — robust to the thin-pin cosphericity). This is
        # the geometry gmsh 4.13/4.15 leave empty.
        geo = joinpath(@__DIR__, "fixtures", "enclosure_coax_junction.geo")
        boxes = Dict{String,NTuple{6,Float64}}(); cyls = Dict{String,NTuple{7,Float64}}()
        for ln in eachline(geo)
            m = match(r"(\w+)\s*=\s*newv;\s*Box\(\w+\)\s*=\s*\{([^}]+)\}", ln)
            m !== nothing && (boxes[m.captures[1]] = Tuple(parse.(Float64, split(m.captures[2],","))))
            m = match(r"(\w+)\s*=\s*newv;\s*Cylinder\(\w+\)\s*=\s*\{([^}]+)\}", ln)
            m !== nothing && (cyls[m.captures[1]] = Tuple(parse.(Float64, split(m.captures[2],","))))
        end
        for k in ("sm_air_inside","sm_case_outer","sm_air","sm_slot"); @test haskey(boxes,k); end
        @test haskey(cyls,"sm_coax_pin")
        ab = boxes["sm_air_inside"]; co = boxes["sm_case_outer"]; pn = cyls["sm_coax_pin"]
        ao = boxes["sm_air"]; sl = boxes["sm_slot"]
        bsurf(b) = box_surface(b[1],b[1]+b[4], b[2],b[2]+b[5], b[3],b[3]+b[6])
        airs  = bsurf(ab); slots = bsurf(sl)
        shell = box_shell_surface(co[1],co[1]+co[4], co[2],co[2]+co[5], co[3],co[3]+co[6],
                                  ab[1],ab[1]+ab[4], ab[2],ab[2]+ab[5], ab[3],ab[3]+ab[6])
        pinl  = sqrt(pn[4]^2+pn[5]^2+pn[6]^2)
        pins  = cylinder_surface((pn[1],pn[2],pn[3]), (pn[4],pn[5],pn[6]), pn[7], pinl; nθ=8, nz=2)
        # pin(1)/slot(2)/air(3)/case(4) — four literal regions, all native Box/Cylinder primitives
        ml = tetrahedralize_conforming_exact([pins, slots, airs, shell])
        @test validate(ml).ok                                                     # the mesh gmsh cannot produce
        tprl = tets_per_region(ml)
        @test length(tprl) == 4 && all(v -> v > 0, values(tprl))                  # all four literal volumes filled (NO empty volume)
        @test is_closed_manifold(ml)
        # total filled volume is positive and bounded by the air outer box (sm_air) that
        # contains the whole assembly — the pin feed-through extends past the case wall.
        @test 0 < mesh_vol(ml) <= ao[4]*ao[5]*ao[6] + 1e-9
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

        # THE FEED-THROUGH: the coax pin now CROSSES the case wall (z=5) — from the air
        # cavity into/through the metal case. Classification [pin,air,case] tags any
        # in-cylinder tet as `pin` whether it is in the cavity or the wall, so the
        # cylindrical BORE is handled with NO explicit CSG bore surface. Both the
        # pin↔air (cavity) and pin↔case (bore) interfaces appear and conform, and the
        # air↔case cavity boundary conforms around the bore hole.
        nθ, nz = 8, 3
        pinf = cylinder_surface((3.,3.,2.), (0.,0.,1.), 0.7, 3.5; nθ=nθ, nz=nz)   # z 2..5.5 crosses z=5
        mf = tetrahedralize_conforming([pinf, cavity, outer])
        @test validate(mf).ok
        @test mesh_vol(mf) ≈ 216.0 rtol=1e-9                     # total = outer box, EXACT
        tpf = tets_per_region(mf)
        @test length(tpf) == 3 && all(v > 0 for v in values(tpf))
        ff = facetags(mf)
        @test all(length(v) <= 2 for v in values(ff))           # manifold everywhere
        pin_air  = count(v -> sort(v) == Int32[1,2], values(ff))
        pin_case = count(v -> sort(v) == Int32[1,3], values(ff))
        @test pin_air > 0 && pin_case > 0                        # pin passes THROUGH the wall — both interfaces
        @test pin_air + pin_case == 2*nθ*nz                      # every pin face conforms (none lost at the wall)
        @test count(v -> sort(v) == Int32[2,3], values(ff)) > 0  # air/case boundary conforms around the bore

        # THE ACCEPTANCE CASE at the LITERAL `.geo` scale (metres): air / metal case /
        # coax pin (0.8 mm radius, aspect ≈199), the coax feed-through where gmsh 4.13/
        # 4.15 leave every volume empty. Meshed as ONE conforming partition — all three
        # volumes filled, the pin bore through the case wall conforming, exact total
        # volume (= the air box, which contains everything). (Coarse pin here to keep the
        # test quick; the real nθ=12·nz=40 resolution is verified out-of-band.)
        airL  = box_surface(-0.0405, 0.2810, -0.0405, 0.2010, -0.0405, 0.3810)
        caseL = box_surface(-0.0005, 0.2205, -0.0005, 0.1405, -0.0005, 0.3005)
        pinL  = cylinder_surface((0.17,0.1605,0.15), (0.0,-1.0,0.0), 0.0008, 0.1589; nθ=6, nz=3)
        mL = tetrahedralize_conforming([pinL, caseL, airL])      # classify pin → case → air
        @test validate(mL).ok
        @test mesh_vol(mL) ≈ 0.3215*0.2415*0.4215 rtol=1e-9      # = air-box volume, EXACT
        tpL = tets_per_region(mL)
        @test length(tpL) == 3 && all(v > 0 for v in values(tpL))   # air/case/pin all filled — NOT empty
        fL = facetags(mL)
        @test all(length(v) <= 2 for v in values(fL))            # manifold
        pinL_air  = count(v -> sort(v) == Int32[1,3], values(fL))
        pinL_case = count(v -> sort(v) == Int32[1,2], values(fL))
        @test pinL_air > 0 && pinL_case > 0                      # feed-through: pin in air AND through the case
        @test pinL_air + pinL_case == 2*6*3                      # all pin faces conform at literal scale
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

    # ── Stage-4: mesh_box — size-controlled structured (Kuhn/Freudenthal) box mesh ──
    # Independent oracles: the max-edge GUARANTEE (recomputed from scratch), exact
    # volume + boundary-area (analytic box formulas), a watertight boundary via
    # boundary Euler χ=2 (sphere), sliver-free quality via tet_dihedral_extrema, and
    # CRC determinism. This is a *provably*-correct capability (explicit connectivity,
    # not Delaunay-of-degenerate-lattice — see validation/stage4_size_refinement/).
    @testset "mesh_box: guaranteed max-edge, exact fill, watertight, sliver-free" begin
        boxvol(bx)  = (bx[2]-bx[1])*(bx[4]-bx[3])*(bx[6]-bx[5])
        boxarea(bx) = 2*((bx[2]-bx[1])*(bx[4]-bx[3]) + (bx[4]-bx[3])*(bx[6]-bx[5]) + (bx[2]-bx[1])*(bx[6]-bx[5]))
        maxedge(m) = maximum(begin
            vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]); e=0.0
            for i in 1:4, j in i+1:4
                p=node(m,vs[i]); q=node(m,vs[j]); e=max(e, hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3]))
            end; e end for t in 1:ntets(m))

        @testset "cube: max-edge bound + exact volume + valid + watertight, sweep hmax" begin
            bx = (0.0,4.0,0.0,4.0,0.0,4.0)
            for hmax in (3.0, 2.0, 1.0, 0.6)
                m = mesh_box(bx...; hmax=hmax)
                @test maxedge(m) <= hmax + 1e-9            # THE guarantee, recomputed
                @test mesh_vol(m) ≈ boxvol(bx) rtol=1e-9   # exact fill (independent sum)
                @test validate(m).ok                       # positive tets, non-degenerate, manifold
                @test boundary_euler(m) == 2               # boundary is a 2-sphere ⇒ watertight
            end
        end

        @testset "general (non-cube, shifted) box: volume + boundary-area + bbox oracle" begin
            bx = (-1.0, 2.0, 0.0, 5.0, 1.0, 3.0)           # 3×5×2 = 30, offset from origin
            m = mesh_box(bx...; hmax=0.8)
            @test mesh_vol(m) ≈ boxvol(bx) rtol=1e-9
            @test validate(m).ok
            @test boundary_euler(m) == 2
            bf, _ = boundary_faces(m.tets)
            area = sum(triangle_area(node(m,f[1]),node(m,f[2]),node(m,f[3])) for f in bf)
            @test area ≈ boxarea(bx) rtol=1e-9             # boundary == the 6 box faces exactly
            lo, hi = bounding_box(m)
            @test all(abs.(lo .- (bx[1],bx[3],bx[5])) .< 1e-12)
            @test all(abs.(hi .- (bx[2],bx[4],bx[6])) .< 1e-12)
        end

        @testset "quality: no slivers (bounded min dihedral, radius-edge)" begin
            m = mesh_box(0,4,0,4,0,4; hmax=1.0)
            dmin = minimum(tet_dihedral_extrema(node(m,m.tets[1,t]),node(m,m.tets[2,t]),
                            node(m,m.tets[3,t]),node(m,m.tets[4,t]))[1] for t in 1:ntets(m))
            remax = maximum(tet_radius_edge(node(m,m.tets[1,t]),node(m,m.tets[2,t]),
                            node(m,m.tets[3,t]),node(m,m.tets[4,t])) for t in 1:ntets(m))
            @test dmin > 0.6            # radians (~34°); Kuhn cube min dihedral = π/4 = 45°
            @test remax < 1.0           # radius-edge ratio well below the sliver regime
        end

        @testset "determinism (CRC) + error paths" begin
            @test mesh_crc(mesh_box(0,2,0,2,0,2; hmax=0.7)).sha ==
                  mesh_crc(mesh_box(0,2,0,2,0,2; hmax=0.7)).sha
            @test_throws ArgumentError mesh_box(0,0,0,1,0,1; hmax=0.5)   # zero x-extent
            @test_throws ArgumentError mesh_box(0,1,0,1,0,1; hmax=0.0)   # non-positive hmax
            @test_throws ArgumentError mesh_box(0,1,0,1,0,1; hmax=-1.0)
        end

        @testset "conformal periodic faces (ARRAY-PML class): opposite faces match under period" begin
            # A periodic unit cell needs opposite boundary faces to carry IDENTICAL
            # triangulations so periodic BCs pair nodes 1:1. The Kuhn subdivision is
            # translation-symmetric per cell, so each pair of opposite box faces matches
            # exactly under the period translation — the conformal-periodic-face capability.
            m = mesh_box(0.,2., 0.,2., 0.,2.; hmax=0.7)
            @test validate(m).ok
            bf, _ = boundary_faces(m.tets)
            cd(i) = (m.coords[1,i], m.coords[2,i], m.coords[3,i])
            onplane(f, ax, val) = all(abs(cd(v)[ax]-val) < 1e-12 for v in f)
            key(f, ax, sh) = sort([ntuple(k -> round(cd(v)[k] + (k==ax ? sh : 0.0); digits=9), 3) for v in f])
            for (ax, lo, hi) in ((1,0.,2.),(2,0.,2.),(3,0.,2.))
                loset = Set(key(f, ax, hi-lo) for f in bf if onplane(f, ax, lo))
                hiset = Set(key(f, ax, 0.0)  for f in bf if onplane(f, ax, hi))
                @test !isempty(loset)
                @test loset == hiset            # opposite faces conformal under the period
            end
        end
    end

    # ── Stage-4/5: mesh_box_regions — conforming multi-region + native box CSG ──────
    # Shared-grid Kuhn subdivision with per-cell region classification: union,
    # difference (void), nesting, multi-material — all at guaranteed edge size.
    # Independent oracles: exact per-region and total volumes (analytic), manifold
    # conformity (max face incidence == 2 everywhere ⇒ no gaps/T-junctions across
    # region interfaces), boundary Euler χ (counts boundary components), CRC.
    @testset "mesh_box_regions: conforming multi-region + native box CSG" begin
        maxedge(m) = maximum(begin
            vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]); e=0.0
            for i in 1:4, j in i+1:4
                p=node(m,vs[i]); q=node(m,vs[j]); e=max(e, hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3]))
            end; e end for t in 1:ntets(m))
        voltag(m,tag) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                            for t in 1:ntets(m) if m.tet_tag[t]==tag; init=0.0)

        @testset "nested: air inside case — conforming, exact per-region volumes" begin
            R = [BoxRegion(0,6,0,6,0,6,1), BoxRegion(1,5,1,5,1,5,2)]  # air (tag2) higher priority
            m = mesh_box_regions(R; hmax=1.5)
            @test maxedge(m) <= 1.5 + 1e-9
            @test validate(m).ok                            # positive + manifold ⇒ conforming
            @test sort(unique(m.tet_tag)) == Int32[1,2]
            @test voltag(m,2) ≈ 64.0 rtol=1e-9              # air = 4^3
            @test voltag(m,1) ≈ 216.0 - 64.0 rtol=1e-9      # case shell around the air
            @test mesh_vol(m) ≈ 216.0 rtol=1e-9             # exact union fill
            @test boundary_euler(m) == 2                    # air/case interface is INTERNAL ⇒ one boundary
        end

        @testset "interface conformity: no non-manifold face across a region split" begin
            R = [BoxRegion(0,4,0,4,0,4,1), BoxRegion(0,2,0,4,0,4,2)]  # split at x=2 (a grid plane)
            m = mesh_box_regions(R; hmax=1.0)
            _, maxinc = boundary_faces(m.tets)
            @test maxinc == 2                               # every internal face shared by exactly 2 tets
            @test 1 in m.tet_tag && 2 in m.tet_tag
            @test mesh_vol(m) ≈ 64.0 rtol=1e-9
            @test validate(m).ok
        end

        @testset "difference: hollow shell (box − void) is watertight with a cavity" begin
            R = [BoxRegion(0,6,0,6,0,6,1), BoxRegion(2,4,2,4,2,4,0)]  # inner box is void (tag 0)
            m = mesh_box_regions(R; hmax=1.5)
            @test mesh_vol(m) ≈ 216.0 - 8.0 rtol=1e-9       # 216 − 2^3
            @test validate(m).ok
            @test boundary_euler(m) == 4                    # outer sphere (2) + inner cavity sphere (2)
            @test maxedge(m) <= 1.5 + 1e-9
        end

        @testset "union of two adjacent boxes fuses conformingly" begin
            R = [BoxRegion(0,2,0,2,0,2,1), BoxRegion(2,4,0,2,0,2,1)]  # share the face x=2
            m = mesh_box_regions(R; hmax=1.0)
            @test mesh_vol(m) ≈ 16.0 rtol=1e-9              # 8 + 8, shared face internal
            @test validate(m).ok
            @test boundary_euler(m) == 2                    # one fused solid
        end

        @testset "single region matches mesh_box; determinism; error paths" begin
            @test mesh_vol(mesh_box_regions([BoxRegion(0,4,0,4,0,4,1)]; hmax=1.0)) ≈ 64.0 rtol=1e-9
            @test mesh_crc(mesh_box_regions([BoxRegion(0,2,0,2,0,2,7)]; hmax=0.7)).sha ==
                  mesh_crc(mesh_box_regions([BoxRegion(0,2,0,2,0,2,7)]; hmax=0.7)).sha
            @test_throws ArgumentError mesh_box_regions(BoxRegion[]; hmax=1.0)                  # no regions
            @test_throws ArgumentError mesh_box_regions([BoxRegion(0,0,0,1,0,1,1)]; hmax=0.5)   # degenerate box
            @test_throws ArgumentError mesh_box_regions([BoxRegion(0,1,0,1,0,1,1)]; hmax=0.0)   # bad hmax
            @test_throws ArgumentError mesh_box_regions([BoxRegion(0,2,0,2,0,2,0)]; hmax=1.0)   # all void → empty
        end

        @testset "thin coupling-slot (THIN-SLOT class): valid + exact at extreme aspect" begin
            # Representative of the enclosure 1 mm coupling slot — a THIN conducting slab
            # (the slot/wall) between two air cavities. This is the thin-feature Boolean
            # case gmsh slivers on; the deterministic structured route meshes it VALID with
            # EXACT per-region volumes at aspect ratios up to ~1000:1 (min dihedral degrades
            # gracefully with thinness but every tet stays positive).
            for thick in (1.0, 0.2)
                xc = 10.0
                R = [BoxRegion(0.0,20.0, 0.0,10.0, 0.0,10.0, 1),                 # air enclosure
                     BoxRegion(xc-thick/2, xc+thick/2, 0.0,10.0, 0.0,10.0, 2)]   # thin slot (priority)
                m = mesh_box_regions(R; hmax=1.0)
                @test validate(m).ok                              # valid despite the thin feature
                @test maxedge(m) <= 1.0 + 1e-9                    # size bound still met
                rv = Dict{Int32,Float64}()
                for t in 1:ntets(m)
                    rv[m.tet_tag[t]] = get(rv, m.tet_tag[t], 0.0) +
                        tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                end
                @test rv[Int32(2)] ≈ thick*10*10 rtol=1e-9                  # exact slot volume
                @test rv[Int32(1)] ≈ 20*10*10 - thick*10*10 rtol=1e-9       # exact air volume
                @test tets_per_region(m)[Int32(2)] > 0                      # thin region genuinely filled
            end
        end

        @testset "spiral conductor trace (SPIRAL class): thin winding conductor in a substrate" begin
            # Representative of the silicon-spiral thin swept conductor — a planar square-
            # spiral trace (thin high-aspect winding, connected axis-aligned segments) as a
            # native-CSG union embedded in a substrate. The structured route meshes the
            # winding conductor VALID and CONFORMING with the trace filled as one region.
            w=0.5; h=0.5; L=4.0
            segs=[(0.0,L,0.0,w,0.0,h),(L-w,L,0.0,L,0.0,h),(w,L,L-w,L,0.0,h),(w,2w,w,L,0.0,h)]
            R=BoxRegion[BoxRegion(-0.5,L+0.5,-0.5,L+0.5,-0.5,h+0.5,2)]   # substrate
            for s in segs; push!(R, BoxRegion(s..., 1)); end            # conductor segments
            m=mesh_box_regions(R; hmax=1.0)
            @test validate(m).ok                                        # valid winding
            @test maxedge(m) <= 1.0 + 1e-9
            tpr=tets_per_region(m)
            @test Set(keys(tpr)) == Set(Int32[1,2]) && all(v->v>0, values(tpr))  # conductor + substrate filled
            ftc=Dict{NTuple{3,Int32},Int}()
            for t in 1:ntets(m), k in 1:4
                vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]); f=Tuple(sort(Int32[vs[j] for j in 1:4 if j!=k])); ftc[f]=get(ftc,f,0)+1
            end
            @test maximum(values(ftc)) <= 2                             # conforming (conductor/substrate interface shared)
        end
    end

    # ── Stage-4: mesh_cylinder — uniform size-controlled cylinder (no Delaunay) ──────
    # Structured (r,θ,z) Kuhn mesh with axis collapse: uniform maxedge ≤ hmax, exact
    # faceted volume, watertight, valid — the cospherical-robust route the Delaunay
    # fill can't take. (Near-axis tets are thinner; size/validity guarantees hold.)
    @testset "mesh_cylinder: uniform size-controlled cylinder (cospherical-robust)" begin
        maxedge(m) = maximum(begin
            vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t]); e=0.0
            for i in 1:4, j in i+1:4
                p=node(m,vs[i]); q=node(m,vs[j]); e=max(e, hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3]))
            end; e end for t in 1:ntets(m))
        facetvol(R,H,nt) = 0.5*nt*R^2*sin(2π/nt)*H

        @testset "uniform maxedge ≤ hmax + exact volume + watertight, sweep hmax" begin
            for hmax in (2.0, 1.0, 0.6)
                m = mesh_cylinder((0.,0.,0.),(0.,0.,1.),2.0,4.0; hmax=hmax)
                @test maxedge(m) <= hmax + 1e-9               # THE uniform guarantee
                @test validate(m).ok                          # all-positive, manifold
                @test boundary_euler(m) == 2                  # watertight (cylinder ≅ sphere)
                nt = max(3, ceil(Int, 2π*2.0/(hmax/sqrt(3.0))))
                @test mesh_vol(m) ≈ facetvol(2.0,4.0,nt) rtol=1e-9   # exact faceted-prism volume
            end
        end
        @testset "tilted/offset cylinder valid + watertight; error paths" begin
            m = mesh_cylinder((1.,2.,3.),(1.,1.,1.),1.5,3.0; hmax=1.0)
            @test validate(m).ok
            @test boundary_euler(m) == 2
            @test maxedge(m) <= 1.0 + 1e-9
            @test_throws ArgumentError mesh_cylinder((0.,0.,0.),(0.,0.,1.),0.0,1.0; hmax=0.5)   # R=0
            @test_throws ArgumentError mesh_cylinder((0.,0.,0.),(0.,0.,1.),1.0,1.0; hmax=0.0)   # hmax=0
        end
    end

    # ── Stage-3: recover_boundary — robust boundary recovery / conforming tetra ──────
    # Recover an arbitrary closed PLC surface (possibly NON-convex) as a conforming
    # tet mesh, or throw an explicit blocker (never a silent bad mesh). Independent
    # conformity oracle: the tet-mesh boundary area == the input surface area (the
    # boundary IS the input surface, triangulation-independently), plus exact domain
    # volume, validity, and closed-manifold. Genuinely non-tetrahedralizable inputs
    # (Schönhardt) MUST raise the explicit blocker.
    @testset "recover_boundary: conforming tetrahedralization of arbitrary PLCs" begin
        surfarea(s) = sum(triangle_area(node(s,s.tris[1,t]),node(s,s.tris[2,t]),node(s,s.tris[3,t]))
                          for t in 1:ntris(s); init=0.0)
        bndarea(m)  = (bf = first(boundary_faces(m.tets));
                       sum(triangle_area(node(m,f[1]),node(m,f[2]),node(m,f[3])) for f in bf; init=0.0))
        conforms(s, exactvol) = begin
            m = recover_boundary(s)
            @test validate(m).ok                          # positive, non-degenerate, manifold
            @test is_closed_manifold(m)
            @test bndarea(m) ≈ surfarea(s) rtol=1e-9       # tet-mesh boundary == input surface
            @test mesh_vol(m) ≈ exactvol rtol=1e-9         # exact domain volume ⇒ boundary conformed
            m
        end
        # star-shaped L-prism (square minus a corner, extruded): area 12 × height 2 = 24
        function _Lprism()
            P=[(0.0,0.0),(4.0,0.0),(4.0,2.0),(2.0,2.0),(2.0,4.0),(0.0,4.0)]; n=length(P)
            xs=Float64[];ys=Float64[];zs=Float64[]
            for z in (0.0,2.0), (x,y) in P; push!(xs,x);push!(ys,y);push!(zs,z); end
            tr=NTuple{3,Int}[]; bot(i)=i; top(i)=n+i
            for i in 1:n; j=i%n+1; push!(tr,(bot(i),bot(j),top(j))); push!(tr,(bot(i),top(j),top(i))); end
            for i in 2:n-1; push!(tr,(bot(1),bot(i+1),bot(i))); push!(tr,(top(1),top(i),top(i+1))); end
            C=Matrix{Float64}(undef,3,length(xs)); for k in 1:length(xs); C[:,k]=[xs[k],ys[k],zs[k]]; end
            T=Matrix{Int32}(undef,3,length(tr)); for (t,f) in enumerate(tr); T[:,t]=Int32[f...]; end
            Mesh(C; tris=T)
        end
        # regular octagonal prism (a faceted cylinder), exact octagon area (1/2)·ns·r²·sin(2π/ns)·h
        function _octprism(r,h,ns)
            xs=Float64[];ys=Float64[];zs=Float64[]
            for z in (0.0,h), k in 0:ns-1; push!(xs,r*cos(2pi*k/ns)); push!(ys,r*sin(2pi*k/ns)); push!(zs,z); end
            tr=NTuple{3,Int}[]; bot(i)=i; top(i)=ns+i
            for k in 0:ns-1; a=k+1; b=(k+1)%ns+1; push!(tr,(bot(a),bot(b),top(b))); push!(tr,(bot(a),top(b),top(a))); end
            for k in 1:ns-2; push!(tr,(bot(1),bot(k+2),bot(k+1))); push!(tr,(top(1),top(k+1),top(k+2))); end
            C=Matrix{Float64}(undef,3,length(xs)); for i in 1:length(xs); C[:,i]=[xs[i],ys[i],zs[i]]; end
            T=Matrix{Int32}(undef,3,length(tr)); for (t,f) in enumerate(tr); T[:,t]=Int32[f...]; end
            Mesh(C; tris=T), 0.5*ns*r^2*sin(2pi/ns)*h
        end
        # Schönhardt polyhedron (twist θ, reflex diagonal) — NOT tetrahedralizable w/o Steiner
        function _schonhardt(θ)
            C=Matrix{Float64}(undef,3,6)
            for k in 0:2; C[:,k+1]=[cos(2pi*k/3),sin(2pi*k/3),0.0]; end
            for k in 0:2; C[:,k+4]=[cos(2pi*k/3+θ),sin(2pi*k/3+θ),1.0]; end
            tr=NTuple{3,Int}[(1,3,2),(4,5,6)]
            for i in 1:3
                ai=i; an=i%3+1; bi=i+3; bn=(i%3)+1+3
                push!(tr,(ai,an,bn)); push!(tr,(ai,bn,bi))     # reflex diagonal ai-bn
            end
            T=Matrix{Int32}(undef,3,length(tr)); for (t,f) in enumerate(tr); T[:,t]=Int32[f...]; end
            Mesh(C; tris=T)
        end

        @testset "convex box"                              begin conforms(box_surface(0,4,0,4,0,4), 64.0) end
        @testset "non-convex genus-1 (through-tunnel)"     begin conforms(box_tunnel_surface(0,6,0,6,0,6, 2,4,2,4), 216.0-24.0) end
        @testset "non-convex hollow shell"                 begin conforms(box_shell_surface(0,6,0,6,0,6, 1,5,1,5,1,5), 216.0-64.0) end
        @testset "non-convex star-shaped L-prism"          begin conforms(_Lprism(), 24.0) end
        @testset "faceted cylinder (octagonal prism)"      begin s,v=_octprism(2.0,3.0,8); conforms(s, v) end

        # NON-STAR-SHAPED non-convex prisms (extruded U-channel / comb / star): these have
        # no single kernel point, yet recover_boundary conforms them (Delaunay-recoverable).
        # Demonstrates the breadth — the only unhandled class is non-star AND non-Delaunay-
        # recoverable (exotic). Caps triangulated by Mesh2D CDT (handles non-convex loops).
        _prism(poly, h) = begin
            n=length(poly); xs=Float64[p[1] for p in poly]; ys=Float64[p[2] for p in poly]
            cap = to_mesh(constrained_delaunay(xs, ys, [(i, i%n+1) for i in 1:n]);
                          interior=classify_interior(constrained_delaunay(xs, ys, [(i, i%n+1) for i in 1:n])))
            capnode=[(cap.coords[1,i],cap.coords[2,i]) for i in 1:size(cap.coords,2)]
            pid=Dict((poly[i][1],poly[i][2])=>i for i in 1:n)
            coords=vcat([(p[1],p[2],0.0) for p in poly], [(p[1],p[2],h) for p in poly])
            tris=NTuple{3,Int}[]
            for i in 1:n; j=i%n+1; push!(tris,(i,j,n+j)); push!(tris,(i,n+j,n+i)); end
            for t in 1:ntris(cap)
                a=pid[capnode[cap.tris[1,t]]]; b=pid[capnode[cap.tris[2,t]]]; c=pid[capnode[cap.tris[3,t]]]
                push!(tris,(a,c,b)); push!(tris,(n+a,n+b,n+c))
            end
            C=Matrix{Float64}(undef,3,length(coords)); for (i,p) in enumerate(coords); C[:,i]=[p...]; end
            Tm=Matrix{Int32}(undef,3,length(tris)); for (t,f) in enumerate(tris); Tm[:,t]=Int32[f...]; end
            Mesh(C; tris=Tm)
        end
        _parea(poly)=abs(sum(poly[i][1]*poly[i%length(poly)+1][2]-poly[i%length(poly)+1][1]*poly[i][2]
                             for i in 1:length(poly)))/2
        @testset "non-star-shaped non-convex prisms (U-channel, comb, star)" begin
            U   = [(0.,0.),(3.,0.),(3.,3.),(2.,3.),(2.,1.),(1.,1.),(1.,3.),(0.,3.)]
            comb= [(0.,0.),(5.,0.),(5.,2.),(4.,2.),(4.,1.),(3.,1.),(3.,2.),(2.,2.),(2.,1.),(1.,1.),(1.,2.),(0.,2.)]
            star= [(2.,0.),(2.6,1.4),(4.,1.5),(3.,2.6),(3.4,4.),(2.,3.2),(0.6,4.),(1.,2.6),(0.,1.5),(1.4,1.4)]
            for poly in (U, comb, star)
                conforms(_prism(poly, 1.0), _parea(poly)*1.0)
            end
        end

        @testset "Schönhardt: blocks by default, meshes with steiner=true" begin
            s = _schonhardt(pi/6)
            @test_throws ErrorException recover_boundary(s)          # default: explicit blocker
            m = recover_boundary(s; steiner=true)                    # Steiner fan-tetrahedralization
            @test validate(m).ok
            @test is_closed_manifold(m)
            @test bndarea(m) ≈ surfarea(s) rtol=1e-9                 # boundary == input surface ⇒ conforming
            @test ntets(m) == ntris(s)                               # one tet per facet (fan)
        end

        @testset "exotic non-star + reflex (twisted prism): recover_boundary_cdt CONFORMS it" begin
            # A non-convex polygon extruded WITH A TWIST is both non-star-shaped AND
            # reflex/non-Delaunay-recoverable — the one class the Float64 `recover_boundary`
            # cannot mesh (it still correctly raises its blocker, even with steiner=true).
            # `recover_boundary_cdt` (exact-kernel conforming-Delaunay refinement) CLOSES it:
            # a valid, closed-manifold, exactly-conforming tet mesh via boundary-Steiner
            # points at exact rational positions (impossible with a Float64 kernel).
            U = [(0.,0.),(3.,0.),(3.,3.),(2.,3.),(2.,1.),(1.,1.),(1.,3.),(0.,3.)]
            twisted = begin
                n=length(U); cx=sum(p[1] for p in U)/n; cy=sum(p[2] for p in U)/n; θ=deg2rad(40)
                rot(p)=((p[1]-cx)*cos(θ)-(p[2]-cy)*sin(θ)+cx, (p[1]-cx)*sin(θ)+(p[2]-cy)*cos(θ)+cy)
                xs=Float64[p[1] for p in U]; ys=Float64[p[2] for p in U]
                cdt=constrained_delaunay(xs,ys,[(i,i%n+1) for i in 1:n])
                cap=to_mesh(cdt; interior=classify_interior(cdt))
                capn=[(cap.coords[1,i],cap.coords[2,i]) for i in 1:size(cap.coords,2)]
                pid=Dict((U[i][1],U[i][2])=>i for i in 1:n)
                coords=vcat([(p[1],p[2],0.0) for p in U], [(rot(p)[1],rot(p)[2],1.0) for p in U])
                tr=NTuple{3,Int}[]
                for i in 1:n; j=i%n+1; push!(tr,(i,j,n+j)); push!(tr,(i,n+j,n+i)); end
                for t in 1:ntris(cap)
                    a=pid[capn[cap.tris[1,t]]]; b=pid[capn[cap.tris[2,t]]]; c=pid[capn[cap.tris[3,t]]]
                    push!(tr,(a,c,b)); push!(tr,(n+a,n+b,n+c))
                end
                C=Matrix{Float64}(undef,3,length(coords)); for (i,p) in enumerate(coords); C[:,i]=[p...]; end
                Tm=Matrix{Int32}(undef,3,length(tr)); for (t,f) in enumerate(tr); Tm[:,t]=Int32[f...]; end
                Mesh(C; tris=Tm)
            end
            # the Float64 recover_boundary still declines (its safe blocker) ...
            @test_throws ErrorException recover_boundary(twisted; max_seeds=12)
            @test_throws ErrorException recover_boundary(twisted; max_seeds=12, steiner=true)
            # ... but recover_boundary_cdt CONFORMS it (exact-kernel Steiner recovery):
            mc = recover_boundary_cdt(twisted)
            @test validate(mc).ok
            @test is_closed_manifold(mc)
            @test bndarea(mc) ≈ surfarea(twisted) rtol=1e-9      # boundary == input surface ⇒ conforming
            @test ntets(mc) > ntris(twisted)                      # genuine interior fill (Steiner points added)
            # exact twist-preserved volume: cross-section area 7 stays, faceted solid ≈ 4.4186
            @test mesh_vol(mc) ≈ 4.418609603270241 rtol=1e-9
        end
        @testset "recover_boundary_cdt is general (supported classes too)" begin
            # the same exact CDT recovery conforms the supported classes with the exact volume.
            for (s, v) in ((box_surface(0,4,0,4,0,4), 64.0),
                           (box_tunnel_surface(0,6,0,6,0,6,2,4,2,4), 192.0),
                           (box_shell_surface(0,6,0,6,0,6,1,5,1,5,1,5), 152.0))
                m = recover_boundary_cdt(s)
                @test validate(m).ok && is_closed_manifold(m)
                @test bndarea(m) ≈ surfarea(s) rtol=1e-9
                @test mesh_vol(m) ≈ v rtol=1e-9
            end
        end
        @testset "error paths" begin
            @test_throws ArgumentError recover_boundary(Mesh(Float64[0 1 0; 0 0 1; 0 0 0]; tris=reshape(Int32[1,2,3],3,1)))
        end
    end

    # ── Stage-4: mesh_sized_conforming — interior size control for curved domains ────
    # Adds an inset interior Steiner lattice to recover_boundary, gated by the exact
    # conformity+validity check. For well-conditioned (thick) curved domains it gives
    # a conforming mesh with interior edges ≤ hmax; thin/cospherical inputs degrade to
    # conforming-only or raise an explicit blocker — never a silent invalid mesh.
    @testset "mesh_sized_conforming: interior size control for curved domains" begin
        _sfa(s)=sum(triangle_area(node(s,s.tris[1,t]),node(s,s.tris[2,t]),node(s,s.tris[3,t])) for t in 1:ntris(s); init=0.0)
        _bfa(m)=(bf=first(boundary_faces(m.tets)); sum(triangle_area(node(m,f[1]),node(m,f[2]),node(m,f[3])) for f in bf; init=0.0))
        _intmax(m)=begin
            isb=falses(nnodes(m)); for f in first(boundary_faces(m.tets)); isb[f[1]]=isb[f[2]]=isb[f[3]]=true; end; mx=0.0
            for t in 1:ntets(m); vs=(m.tets[1,t],m.tets[2,t],m.tets[3,t],m.tets[4,t])
                for i in 1:4,j in i+1:4; (isb[vs[i]]||isb[vs[j]])&&continue; p=node(m,vs[i]);q=node(m,vs[j]); mx=max(mx,hypot(p[1]-q[1],p[2]-q[2],p[3]-q[3])); end
            end; mx
        end
        _icosphere(R,k)=begin
            t=(1+sqrt(5))/2
            V=[(-1.,t,0.),(1.,t,0.),(-1.,-t,0.),(1.,-t,0.),(0.,-1.,t),(0.,1.,t),(0.,-1.,-t),(0.,1.,-t),(t,0.,-1.),(t,0.,1.),(-t,0.,-1.),(-t,0.,1.)]
            F=[(1,12,6),(1,6,2),(1,2,8),(1,8,11),(1,11,12),(2,6,10),(6,12,5),(12,11,3),(11,8,7),(8,2,9),(4,10,5),(4,5,3),(4,3,7),(4,7,9),(4,9,10),(5,10,6),(3,5,12),(7,3,11),(9,7,8),(10,9,2)]
            coords=[v for v in V]; idx=Dict{NTuple{3,Float64},Int}(); for (i,v) in enumerate(coords); idx[v]=i; end; faces=F
            for _ in 1:k
                vid(p)=get!(idx,p) do; push!(coords,p); length(coords); end; nf=NTuple{3,Int}[]
                for (a,b,c) in faces; pa=coords[a];pb=coords[b];pc=coords[c]; md(x,y)=((x[1]+y[1])/2,(x[2]+y[2])/2,(x[3]+y[3])/2)
                    ab=vid(md(pa,pb));bc=vid(md(pb,pc));ca=vid(md(pc,pa)); push!(nf,(a,ab,ca));push!(nf,(ab,b,bc));push!(nf,(ca,bc,c));push!(nf,(ab,bc,ca)); end
                faces=nf; end
            pr=[(v[1]/hypot(v...)*R,v[2]/hypot(v...)*R,v[3]/hypot(v...)*R) for v in coords]
            C=Matrix{Float64}(undef,3,length(pr)); for (i,p) in enumerate(pr); C[:,i]=[p...]; end
            Tm=Matrix{Int32}(undef,3,length(faces)); for (t,f) in enumerate(faces); Tm[:,t]=Int32[f...]; end
            Mesh(C; tris=Tm)
        end
        @testset "sphere (thick curved): conforming + genuine interior size reduction" begin
            s=_icosphere(3.0,2); base=recover_boundary(s); m=mesh_sized_conforming(s; hmax=1.5)
            @test validate(m).ok
            @test is_closed_manifold(m)
            @test _bfa(m) ≈ _sfa(s) rtol=1e-9            # boundary conforms to the input surface
            @test _intmax(m) <= 1.5 + 1e-9               # interior edges size-controlled
            @test ntets(m) > ntets(base)                 # genuinely refined vs the no-interior baseline
        end
        @testset "thin/cospherical input degrades safely (conforming, never silent-invalid)" begin
            cyl=cylinder_surface((0.,0,0),(0.,0,1),2.0,4.0; nθ=16, nz=4)
            m=mesh_sized_conforming(cyl; hmax=1.5)       # thin ⇒ few/no interior pts, but still conforming+valid
            @test validate(m).ok
            @test is_closed_manifold(m)
            @test _bfa(m) ≈ _sfa(cyl) rtol=1e-9
        end
        @testset "error path" begin
            @test_throws ArgumentError mesh_sized_conforming(box_surface(0,1,0,1,0,1); hmax=0.0)
        end

        @testset "mesh_sized_cdt: arbitrary-surface interior sizing on the exact CDT engine" begin
            # #8 — uniform interior size control on an arbitrary curved domain, built on the
            # exact conforming-Delaunay recovery. Conforming (exact), valid, closed, with
            # interior edges driven toward ≤ hmax; the certify gate guarantees it never
            # returns a non-conforming mesh (falls back to the conforming baseline if the
            # interior lattice would break conformity).
            s = _icosphere(3.0, 1)
            base = recover_boundary_cdt(s)
            m = mesh_sized_cdt(s; hmax=2.5)
            @test validate(m).ok
            @test is_closed_manifold(m)
            @test _bfa(m) ≈ _sfa(s) rtol=1e-9                 # boundary conforms exactly
            @test _intmax(m) <= 2.5 + 1e-9                    # interior edges size-controlled
            @test ntets(m) >= ntets(base)                     # refined (or safe baseline), never worse
            @test_throws ArgumentError mesh_sized_cdt(s; hmax=0.0)
        end
    end

    # ── Stage-5: mesh_boolean — native mesh-CSG (union/intersection/difference) ──────
    # No OCC. Independent oracle: divergence-theorem volume of the RESULT surface ==
    # the analytic Boolean volume; the result is watertight; and it fills into a valid
    # tet mesh via recover_boundary. Axis-aligned inputs (exact plane-arrangement path,
    # handles coplanar shared faces); box×cylinder exercises the tri-tri + CDT path.
    @testset "mesh_boolean: native mesh-CSG (no OCC)" begin
        divvol(s) = abs(sum(begin
            a=node(s,s.tris[1,t]); b=node(s,s.tris[2,t]); c=node(s,s.tris[3,t])
            a[1]*(b[2]*c[3]-b[3]*c[2]) - a[2]*(b[1]*c[3]-b[3]*c[1]) + a[3]*(b[1]*c[2]-b[2]*c[1])
        end for t in 1:ntris(s); init=0.0) / 6)
        watertight(s) = isempty(first(boundary_edges(s.tris)))
        csg(A,B,op,exactvol) = begin
            R = mesh_boolean(A,B,op)
            @test watertight(R)                          # closed result surface
            @test divvol(R) ≈ exactvol rtol=1e-9         # exact Boolean volume
            @test validate(recover_boundary(R)).ok       # result fills into a valid tet mesh
            R
        end

        @testset "box ∪/∩/∖ box, x-offset (coplanar shared faces)" begin
            A=box_surface(0,4,0,4,0,4); B=box_surface(2,6,0,4,0,4)   # ∩ = [2,4]×[0,4]×[0,4] = 32
            csg(A,B,:union,96.0); csg(A,B,:intersection,32.0); csg(A,B,:difference,32.0)
        end
        @testset "box ∪/∩/∖ box, xy-offset" begin
            A=box_surface(0,4,0,4,0,4); B=box_surface(2,6,2,6,0,4)   # ∩ = [2,4]×[2,4]×[0,4] = 16
            csg(A,B,:union,112.0); csg(A,B,:intersection,16.0); csg(A,B,:difference,48.0)
        end
        @testset "box ∪/∩/∖ box, general position" begin
            A=box_surface(0,10,0,10,0,10); B=box_surface(5,15,5,15,5,15)  # ∩ = [5,10]³ = 125
            csg(A,B,:union,1875.0); csg(A,B,:intersection,125.0); csg(A,B,:difference,875.0)
        end
        @testset "box − cylinder (tri-tri + CDT path): watertight + exact faceted volume" begin
            A=box_surface(0,4,0,4,0,4)
            cyl=cylinder_surface((2.,2.,-1.),(0.,0.,1.),1.0,6.0; nθ=16)   # r=1 through the box
            removed = 0.5*16*1.0^2*sin(2pi/16)*4.0                        # inscribed 16-gon × height 4
            R=mesh_boolean(A,cyl,:difference)
            @test watertight(R)
            @test divvol(R) ≈ 64.0 - removed rtol=1e-9
            @test validate(recover_boundary(R)).ok
        end
        @testset "error paths" begin
            A=box_surface(0,4,0,4,0,4); B=box_surface(2,6,2,6,0,4)
            @test_throws ArgumentError mesh_boolean(A,B,:foo)                       # bad op
            @test_throws ArgumentError mesh_boolean(Mesh(A.coords; tris=A.tris[:,1:end-1]), B, :union)  # open input
        end
    end

    # CRC gate for the grid-accelerated point-in-surface classifier: the grid filter
    # (_inside_grid) MUST return bit-identical parity to the brute-force reference
    # (_inside_surface) on every query, or the domain classification (and thus which
    # tets are kept) would silently change. Cross-checked on flat/curved/genus-1
    # surfaces over a dense random cloud that straddles inside/outside/near-face plus
    # every surface vertex and face centroid (the boundary/degenerate cases).
    @testset "raygrid classifier == brute-force parity (bit-identical)" begin
        dir = Mesh3D._CLASSIFY_DIR
        for surf in (box_surface(0.0,3.0, 0.0,2.0, 0.0,1.0),
                     cylinder_surface((0.0,0.0,0.0),(0.0,0.0,1.0), 1.0, 2.0; nθ=24),
                     cylinder_surface((0.3,-0.2,0.1),(0.2,0.3,1.0), 0.8, 5.0; nθ=40, nz=6),
                     box_tunnel_surface(0.0,4.0, 0.0,3.0, 0.0,2.0, 1.0,3.0, 1.0,2.0))
            g = Mesh3D._raygrid(surf, dir)
            lo = (minimum(surf.coords[1,:]), minimum(surf.coords[2,:]), minimum(surf.coords[3,:]))
            hi = (maximum(surf.coords[1,:]), maximum(surf.coords[2,:]), maximum(surf.coords[3,:]))
            pad = 0.15 .* (hi .- lo) .+ 1e-6
            rr = _R3(UInt64(20260812))
            disagree = 0
            for _ in 1:20000
                p = (lo[1]-pad[1] + _nf(rr)*(hi[1]-lo[1]+2pad[1]),
                     lo[2]-pad[2] + _nf(rr)*(hi[2]-lo[2]+2pad[2]),
                     lo[3]-pad[3] + _nf(rr)*(hi[3]-lo[3]+2pad[3]))
                Mesh3D._inside_surface(p, dir, surf) == Mesh3D._inside_grid(p, g) || (disagree += 1)
            end
            for f in 1:size(surf.tris,2)
                a=surf.tris[1,f]; b=surf.tris[2,f]; c=surf.tris[3,f]
                for pid in (a,b,c)
                    p=(surf.coords[1,pid],surf.coords[2,pid],surf.coords[3,pid])
                    Mesh3D._inside_surface(p,dir,surf)==Mesh3D._inside_grid(p,g) || (disagree+=1)
                end
                cx=(surf.coords[1,a]+surf.coords[1,b]+surf.coords[1,c])/3
                cy=(surf.coords[2,a]+surf.coords[2,b]+surf.coords[2,c])/3
                cz=(surf.coords[3,a]+surf.coords[3,b]+surf.coords[3,c])/3
                Mesh3D._inside_surface((cx,cy,cz),dir,surf)==Mesh3D._inside_grid((cx,cy,cz),g) || (disagree+=1)
            end
            @test disagree == 0
        end
    end

    @testset "audit-fix regressions (2026-08-12): silent-invalid-mesh guards" begin
        # F1: mesh_cylinder's zero-volume drop must be scale-independent. A small/thin
        # cylinder must stay watertight (boundary χ=2) with the exact faceted-prism
        # volume — an ABSOLUTE volume threshold silently dropped legitimate near-axis
        # tets at small scale, leaving an axial void (χ=0, wrong volume, no error).
        facet_vol(R,H,nθ) = 0.5*nθ*R^2*sin(2π/nθ)*H
        for (R,H,hmax) in ((1.0, 2.0, 0.2), (5e-4, 9e-5, 5.2e-5))
            m = mesh_cylinder((0.,0.,0.),(0.,0.,1.), R, H; hmax=hmax)
            a = hmax/sqrt(3.0); nθ = max(3, ceil(Int, 2π*R/a))
            @test validate(m).ok
            @test boundary_euler(m) == 2                          # watertight at every scale
            @test mesh_vol(m) ≈ facet_vol(R,H,nθ) rtol=1e-9       # exact faceted volume
        end

        # F2: tetrahedralize_conforming must raise an explicit blocker — never return a
        # silently invalid (zero-volume / non-manifold) mesh — on a cospherical
        # axis-aligned box assembly the exact kernel cannot break into positive tets.
        boxes2 = [box_surface(Float64(ix),Float64(ix+1),Float64(iy),Float64(iy+1),Float64(iz),Float64(iz+1))
                  for ix in 0:1 for iy in 0:1 for iz in 0:1]
        @test_throws ErrorException tetrahedralize_conforming(boxes2)
        # ...and the EXACT-coordinate path MESHES the same cospherical assembly (the exact
        # kernel breaks the degeneracy validly where the Float64 perturb=false path must
        # raise its blocker): valid, all 8 regions filled, conforming, exact volume.
        me = tetrahedralize_conforming_exact(boxes2)
        @test validate(me).ok
        tpe = tets_per_region(me)
        @test length(tpe) == 8 && all(v -> v > 0, values(tpe))
        @test mesh_vol(me) ≈ 8.0 rtol=1e-9
        ftc = Dict{NTuple{3,Int32},Int}()
        for t in 1:ntets(me), k in 1:4
            vs=(me.tets[1,t],me.tets[2,t],me.tets[3,t],me.tets[4,t]); f=Tuple(sort(Int32[vs[j] for j in 1:4 if j != k]))
            ftc[f] = get(ftc, f, 0) + 1
        end
        @test maximum(values(ftc)) <= 2                       # conforming: no face shared by >2 tets
    end

    @testset "exact-coordinate Delaunay kernel (Rational{BigInt})" begin
        # A valid Delaunay tetrahedralization on exact rational coords (empty circumsphere
        # exact; valid closed-manifold positive-volume mesh), breaking cospherical/coplanar
        # degeneracies deterministically with NO jitter — the foundation for boundary
        # recovery with Steiner points at non-Float64 rational positions.
        Q(x) = Rational{BigInt}(Float64(x))
        r = _R3(UInt64(12345))
        for n in (30, 60)
            pts = [(Q(_nf(r)), Q(_nf(r)), Q(_nf(r))) for _ in 1:n]
            tets = delaunay3d_exact(pts)
            ok, nv = is_delaunay_exact(pts, tets)
            @test ok && nv == 0                                   # exact empty circumsphere
            C = Matrix{Float64}(undef, 3, n)
            for i in 1:n; C[1,i]=Float64(pts[i][1]); C[2,i]=Float64(pts[i][2]); C[3,i]=Float64(pts[i][3]); end
            Tm = Matrix{Int32}(undef, 4, length(tets))
            for (j,t) in enumerate(tets); Tm[1,j]=t[1]; Tm[2,j]=t[2]; Tm[3,j]=t[3]; Tm[4,j]=t[4]; end
            m = Mesh(C; tets=Tm)
            @test validate(m).ok                                  # positive volumes + manifold
            @test boundary_euler(m) == 2                          # convex-hull boundary is a 2-sphere
        end
        # maximally cospherical: a 3×3×3 integer box grid — exact fill to the exact volume,
        # zero empty-circumsphere violations (the case the Float64 perturb=false kernel can
        # muddle with degenerate tets).
        boxpts = NTuple{3,Rational{BigInt}}[]
        for x in 0:2, y in 0:2, z in 0:2; push!(boxpts, (Q(x), Q(y), Q(z))); end
        bt = delaunay3d_exact(boxpts)
        okb, _ = is_delaunay_exact(boxpts, bt)
        @test okb
        vol6 = sum(begin
            a=boxpts[t[1]]; b=boxpts[t[2]]; c=boxpts[t[3]]; d=boxpts[t[4]]
            abs((b[1]-a[1])*((c[2]-a[2])*(d[3]-a[3])-(c[3]-a[3])*(d[2]-a[2])) -
                (b[2]-a[2])*((c[1]-a[1])*(d[3]-a[3])-(c[3]-a[3])*(d[1]-a[1])) +
                (b[3]-a[3])*((c[1]-a[1])*(d[2]-a[2])-(c[2]-a[2])*(d[1]-a[1])))
        end for t in bt)
        @test vol6 == 48                                          # exact box volume 8 (×6)
    end
end
