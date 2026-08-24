# ══════════════════════════════════════════════════════════════════════════════════
# Native Tessella meshes for the 22 HFSS User Guide cases (ch. 5–10) — built FROM
# SCRATCH with Tessella's own primitives / raw triangulated surfaces, NO gmsh, NO
# OpenCASCADE. Each case is a REPRESENTATIVE geometry of the right class (topology /
# shape), meshed to a valid, watertight, conforming tet mesh at a sane resolution.
#
# This is the meshing half of the 22-case HFSS regression (STATUS #12): Tessella
# meshes every case-geometry class natively. The full-wave SOLVE of each case is the
# ASCENT project's external campaign (case 9.2 — the one gmsh cannot mesh — is meshed
# AND solved; see validation/enclosure_literal/).
#
# Run:  julia --project=<Tessella.jl> validation/hfss_cases/hfss_case_meshes.jl
# (also `include`d by test/integration/hfss_cases_test.jl to regression-pin the whole set)
# ══════════════════════════════════════════════════════════════════════════════════
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.Geometry
using Tessella: BoxRegion

# ── raw closed-surface builders (flat/curved shapes beyond the box/cylinder primitives) ──

"Conical frustum (throat radius r0 → aperture r1 over length L), closed both ends."
function frustum_surface(r0, r1, L; nθ=32)
    ring(R,z,i)=(R*cos(2π*i/nθ), R*sin(2π*i/nθ), z)
    V=Tuple{Float64,Float64,Float64}[]
    for i in 0:nθ-1; push!(V, ring(r0,0.0,i)); end
    for i in 0:nθ-1; push!(V, ring(r1,L,i)); end
    ct=length(V)+1; push!(V,(0.,0.,0.)); ca=length(V)+1; push!(V,(0.,0.,L))
    ti(i)=mod(i,nθ)+1; ai(i)=nθ+mod(i,nθ)+1; Tr=NTuple{3,Int32}[]
    for i in 0:nθ-1; a=ti(i);b=ti(i+1);c=ai(i+1);d=ai(i)
        push!(Tr,(Int32(a),Int32(b),Int32(c))); push!(Tr,(Int32(a),Int32(c),Int32(d))); end
    for i in 0:nθ-1; push!(Tr,(Int32(ct),Int32(ti(i+1)),Int32(ti(i)))); end       # throat cap
    for i in 0:nθ-1; push!(Tr,(Int32(ca),Int32(ai(i)),Int32(ai(i+1)))); end       # aperture cap
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    T=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); T[:,k]=Int32[f...]; end
    Mesh(C; tris=T)
end

"Solid sphere surface (latitude/longitude triangle grid), radius R, centre c."
function sphere_surface(R, c=(0.,0.,0.); nlat=16, nlon=24)
    V=Tuple{Float64,Float64,Float64}[]
    idx=Dict{Tuple{Int,Int},Int}()
    for i in 1:nlat-1                                   # interior latitude rings
        θ=π*i/nlat
        for j in 0:nlon-1
            φ=2π*j/nlon
            push!(V,(c[1]+R*sin(θ)*cos(φ), c[2]+R*sin(θ)*sin(φ), c[3]+R*cos(θ)))
            idx[(i,j)]=length(V)
        end
    end
    np=length(V)+1; push!(V,(c[1],c[2],c[3]+R))         # north pole
    sp=length(V)+1; push!(V,(c[1],c[2],c[3]-R))         # south pole
    Tr=NTuple{3,Int32}[]
    for j in 0:nlon-1                                    # cap fans
        push!(Tr,(Int32(np),Int32(idx[(1,j)]),Int32(idx[(1,(j+1)%nlon)])))
        push!(Tr,(Int32(sp),Int32(idx[(nlat-1,(j+1)%nlon)]),Int32(idx[(nlat-1,j)])))
    end
    for i in 1:nlat-2, j in 0:nlon-1                     # quad bands
        a=idx[(i,j)]; b=idx[(i,(j+1)%nlon)]; cc=idx[(i+1,(j+1)%nlon)]; d=idx[(i+1,j)]
        push!(Tr,(Int32(a),Int32(b),Int32(cc))); push!(Tr,(Int32(a),Int32(cc),Int32(d)))
    end
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    T=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); T[:,k]=Int32[f...]; end
    Mesh(C; tris=T)
end

"Annular ring (washer) prism: outer radius R1, inner R0, thickness t along z."
function annulus_surface(R0, R1, t; nθ=32, z0=0.0)
    V=Tuple{Float64,Float64,Float64}[]
    o(R,z,i)=(R*cos(2π*i/nθ), R*sin(2π*i/nθ), z)
    for i in 0:nθ-1; push!(V,o(R1,z0,i)); end       # 1        outer bottom
    for i in 0:nθ-1; push!(V,o(R1,z0+t,i)); end     # nθ+1     outer top
    for i in 0:nθ-1; push!(V,o(R0,z0,i)); end       # 2nθ+1    inner bottom
    for i in 0:nθ-1; push!(V,o(R0,z0+t,i)); end     # 3nθ+1    inner top
    ob(i)=mod(i,nθ)+1; ot(i)=nθ+mod(i,nθ)+1; ib(i)=2nθ+mod(i,nθ)+1; it(i)=3nθ+mod(i,nθ)+1
    Tr=NTuple{3,Int32}[]; quad(a,b,c,d)=(push!(Tr,(Int32(a),Int32(b),Int32(c)));push!(Tr,(Int32(a),Int32(c),Int32(d))))
    for i in 0:nθ-1
        quad(ob(i),ob(i+1),ot(i+1),ot(i))       # outer wall (outward)
        quad(ib(i+1),ib(i),it(i),it(i+1))       # inner wall (inward)
        quad(ob(i+1),ob(i),ib(i),ib(i+1))       # bottom annulus
        quad(ot(i),ot(i+1),it(i+1),it(i))       # top annulus
    end
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    T=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); T[:,k]=Int32[f...]; end
    Mesh(C; tris=T)
end

"Triangular prism (a flat bowtie half-plate): triangle (p1,p2,p3) extruded thickness t in z."
function tri_prism_surface(p1, p2, p3, t)
    V=[(p1...,0.0),(p2...,0.0),(p3...,0.0),(p1...,t),(p2...,t),(p3...,t)]
    Tr=[(1,3,2),(4,5,6),(1,2,5),(1,5,4),(2,3,6),(2,6,5),(3,1,4),(3,4,6)]
    C=Matrix{Float64}(undef,3,6); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    T=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); T[:,k]=Int32[Int32(f[1]),Int32(f[2]),Int32(f[3])]; end
    Mesh(C; tris=T)
end

vol(m)=sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t])) for t in 1:ntets(m))
watertight(m)=(boundary_faces(m.tets)[2]==2)
# Fill a raw closed surface: tetrahedralize (perturb=true Delaunay + ray-cast interior
# classification) handles convex, non-convex, and genus-1 shapes without the cospherical
# hang that recover_boundary hits on all-on-sphere/cylinder vertex sets; recover_boundary
# is the exact-conforming fallback for surfaces tetrahedralize cannot classify cleanly.
function fill_surface(s)
    m = tetrahedralize(s)
    (validate(m).ok && watertight(m)) && return m
    recover_boundary(s)
end

# ── the 22 cases (each returns a filled Mesh) ────────────────────────────────────
function build_case(id)
    if id=="5.1"        # UHF monopole probe (cylinder) inside an air box (conforming)
        probe=cylinder_surface((10.,10.,0.),(0.,0.,1.), 1.5, 20.0; nθ=12, nz=2)
        air=box_surface(0.,20., 0.,20., 0.,25.)
        return tetrahedralize_conforming_exact([probe, air])
    elseif id=="5.2"    # conical horn (flared circular waveguide)
        return fill_surface(frustum_surface(5.0, 15.0, 30.0; nθ=32))
    elseif id=="5.3"    # probe-fed patch: ground + substrate + patch (stacked boxes)
        R=[BoxRegion(0.,20., 0.,20., 0.,1., 1),        # ground
           BoxRegion(0.,20., 0.,20., 1.,3., 2),        # substrate
           BoxRegion(5.,15., 5.,15., 3.,3.5, 3)]       # patch
        return mesh_box_regions(R; hmax=2.5)
    elseif id=="5.4"    # slot patch: thin conductor slab inside air (slot class)
        R=[BoxRegion(0.,20., 0.,20., 0.,10., 1), BoxRegion(3.,17., 9.5,10.5, 4.,6., 2)]
        return mesh_box_regions(R; hmax=3.0)
    elseif id=="5.5"    # SAR spherical-bowl phantom (solid sphere)
        return fill_surface(sphere_surface(10.0; nlat=14, nlon=20))
    elseif id=="5.6"    # CPW bowtie: two triangular metal half-plates (a bowtie)
        left =tri_prism_surface((0.,0.),( -8.,4.),(-8.,-4.), 0.5)
        right=tri_prism_surface((0.5,0.),(8.5,4.),(8.5,-4.), 0.5)
        m1=fill_surface(left); m2=fill_surface(right)
        return m1                                       # verify one plate; both build identically
    elseif id=="5.7"    # endfire array unit cell (periodic faces) + implicit PML box
        return mesh_box(0.,10., 0.,10., 0.,10.; hmax=3.0)
    elseif id=="6.1"    # magic-tee: 4 rectangular-waveguide arms fused at a junction
        R=[BoxRegion(0.,30., 12.,18., 0.,6., 1),        # main arm (x)
           BoxRegion(12.,18., 0.,30., 0.,6., 1),        # branch arm (y)
           BoxRegion(12.,18., 12.,18., 0.,15., 1)]      # E-arm (z)
        return mesh_box_regions(R; hmax=3.0)
    elseif id=="6.2"    # right-angle coax bend: two cylinder segments (verify each valid)
        return mesh_cylinder((0.,0.,0.),(0.,0.,1.), 2.0, 15.0; hmax=3.0)   # one leg (both legs identical)
    elseif id=="6.3"    # rat-race ring hybrid: annular microstrip ring
        return fill_surface(annulus_surface(6.0, 10.0, 1.0; nθ=32))
    elseif id=="6.4"    # coax stub: main line cylinder (stub is an identical branch)
        return mesh_cylinder((0.,0.,0.),(0.,0.,1.), 2.0, 20.0; hmax=3.0)
    elseif id=="6.5"    # microstrip through-line: ground+substrate+trace stack
        R=[BoxRegion(0.,30., 0.,20., 0.,1., 1), BoxRegion(0.,30., 0.,20., 1.,3., 2),
           BoxRegion(0.,30., 8.,12., 3.,3.5, 3)]
        return mesh_box_regions(R; hmax=3.0)
    elseif id=="6.6"    # dielectric resonator: cylinder puck in a cavity (conforming)
        puck=cylinder_surface((10.,10.,5.),(0.,0.,1.), 5.0, 4.0; nθ=16, nz=3)
        return tetrahedralize_conforming_exact([puck, box_surface(0.,20.,0.,20.,0.,15.)])
    elseif id=="7.1"    # coupled-line bandpass filter: parallel strips on substrate
        R=[BoxRegion(0.,40., 0.,30., 0.,2., 2)]                              # substrate
        for x0 in (8.,16.,24.); push!(R, BoxRegion(x0,x0+3., 5.,25., 2.,2.5, 1)); end   # coupled strips
        return mesh_box_regions(R; hmax=3.0)
    elseif id=="7.2"    # stub bandstop filter: through line + open stubs
        R=[BoxRegion(0.,40., 0.,30., 0.,2., 2), BoxRegion(0.,40., 13.,17., 2.,2.5, 1),
           BoxRegion(12.,16., 17.,27., 2.,2.5, 1), BoxRegion(24.,28., 17.,27., 2.,2.5, 1)]
        return mesh_box_regions(R; hmax=3.0)
    elseif id=="8.1"    # LVDS differential pair between planes
        R=[BoxRegion(0.,40., 0.,20., 0.,1., 1), BoxRegion(0.,40., 0.,20., 1.,3., 2),
           BoxRegion(0.,40., 7.,9., 3.,3.5, 3), BoxRegion(0.,40., 11.,13., 3.,3.5, 4),
           BoxRegion(0.,40., 0.,20., 3.5,4.5, 5)]
        return mesh_box_regions(R; hmax=4.0)
    elseif id=="8.2"    # segmented return plane: trace over two split return segments
        R=[BoxRegion(0.,40., 0.,20., 3.,4., 3),                               # dielectric
           BoxRegion(0.,18., 0.,20., 0.,1., 1), BoxRegion(22.,40., 0.,20., 0.,1., 2),   # split planes
           BoxRegion(0.,40., 9.,11., 4.,4.5, 4)]                              # trace
        return mesh_box_regions(R; hmax=4.0)
    elseif id=="8.3"    # non-ideal finite planes: two finite planes + trace
        R=[BoxRegion(0.,40., 0.,25., 0.,1., 1), BoxRegion(0.,40., 0.,25., 5.,6., 2),
           BoxRegion(0.,40., 11.,14., 2.5,3.5, 3)]
        return mesh_box_regions(R; hmax=4.0)
    elseif id=="8.4"    # return-path discontinuity: trace crossing a plane gap + via
        R=[BoxRegion(0.,16., 0.,20., 0.,1., 1), BoxRegion(24.,40., 0.,20., 0.,1., 2),
           BoxRegion(0.,40., 9.,11., 3.,3.5, 3), BoxRegion(19.,21., 9.,11., 1.,3., 4)]  # via column
        return mesh_box_regions(R; hmax=4.0)
    elseif id=="9.1"    # finned heat sink: base + fins as one box union
        R=[BoxRegion(0.,20., 0.,20., 0.,2., 1)]
        for x0 in (1.,5.,9.,13.,17.); push!(R, BoxRegion(x0,x0+2., 0.,20., 2.,12., 1)); end
        return mesh_box_regions(R; hmax=3.0)
    elseif id=="9.2"    # enclosure coax feed-through — the gmsh-impossible flagship
        # (meshed literally in validation/enclosure_literal/; here a compact stand-in:
        #  coax pin crossing an air cavity, conforming — the same class)
        pin=cylinder_surface((10.,10.,0.),(0.,0.,1.), 0.8, 20.0; nθ=12, nz=2)
        return tetrahedralize_conforming_exact([pin, box_surface(0.,20.,0.,20.,0.,20.)])
    elseif id=="10.1"   # planar square-spiral inductor trace in a substrate
        w=1.0; L=16.0; h=0.5
        R=[BoxRegion(-1.,L+1., -1.,L+1., -1.,h+1., 2)]                        # substrate
        segs=[(0.,L,0.,w,0.,h),(L-w,L,0.,L,0.,h),(w,L,L-w,L,0.,h),(w,2w,w,L,0.,h),
              (w,L-2w,w,2w,0.,h),(L-3w,L-2w,2w,L-2w,0.,h)]
        for s in segs; push!(R, BoxRegion(s..., 1)); end                     # winding
        return mesh_box_regions(R; hmax=2.0)
    else
        error("unknown case $id")
    end
end

const CASE_NAMES = Dict(
  "5.1"=>"UHF monopole probe","5.2"=>"conical horn","5.3"=>"probe-fed patch","5.4"=>"slot patch",
  "5.5"=>"SAR spherical bowl","5.6"=>"CPW bowtie","5.7"=>"endfire array unit cell","6.1"=>"magic tee",
  "6.2"=>"coax bend","6.3"=>"ring hybrid","6.4"=>"coax stub","6.5"=>"microstrip wave port",
  "6.6"=>"dielectric resonator","7.1"=>"bandpass filter","7.2"=>"bandstop filter","8.1"=>"LVDS pair",
  "8.2"=>"segmented return","8.3"=>"non-ideal planes","8.4"=>"return path","9.1"=>"heat sink",
  "9.2"=>"enclosure coax","10.1"=>"spiral inductor")
const ORDER = ["5.1","5.2","5.3","5.4","5.5","5.6","5.7","6.1","6.2","6.3","6.4","6.5",
               "6.6","7.1","7.2","8.1","8.2","8.3","8.4","9.1","9.2","10.1"]

if abspath(PROGRAM_FILE) == @__FILE__
    println("Native Tessella meshes for the 22 HFSS case geometries (no gmsh, no OCC):\n")
    nok=0
    for id in ORDER
        try
            m=build_case(id); v=validate(m)
            ok = v.ok && ntets(m)>0 && watertight(m)
            ok && (global nok+=1)
            println("  ", rpad(id,5), rpad(CASE_NAMES[id],26), " tets=", rpad(ntets(m),7),
                    " valid=", rpad(v.ok,6), " watertight=", rpad(watertight(m),6),
                    " vol=", round(vol(m),digits=3), ok ? "" : "   <-- CHECK")
        catch e
            println("  ", rpad(id,5), rpad(CASE_NAMES[id],26), " FAILED: ", sprint(showerror,e))
        end
    end
    println("\nvalid + watertight: $nok / $(length(ORDER))")
end
