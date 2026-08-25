# ── Integration CRC suite: the top-level mesh_volume pipeline ───────────────────
#
# Correctness  : end-to-end surface → volume mesh at the exact expected volume.
# Robustness   : the "validated or explicit blocker" contract — a defective input
#                raises a precise error, never a silent bad mesh.
# Completeness : output validates; smoothing preserves volume.

using Test
using Tessella
using Tessella.MeshTypes
using Tessella.IO
using Tessella.Geometry
using Tessella.Mesh3D: tetrahedralize_multi, tetrahedralize_conforming

mvpvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                for t in 1:ntets(m); init=0.0)

function _cube_surface()
    C=Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    F=[(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
    t=Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    Mesh(C; tris=t)
end

struct _UnreadPipelinePoints <: AbstractVector{Float64}
    count::Int
end
Base.size(values::_UnreadPipelinePoints) = (values.count,)
Base.IndexStyle(::Type{_UnreadPipelinePoints}) = IndexLinear()
Base.getindex(::_UnreadPipelinePoints, index::Int) =
    throw(ErrorException("pipeline resource preflight read point $index"))

@testset "mesh_volume pipeline (integration)" begin
    @testset "clean surface → validated volume mesh" begin
        m = mesh_volume(_cube_surface())
        @test validate(m).ok
        @test mvpvol(m) ≈ 1.0 rtol=1e-6
        @test mesh_quality(m).n_tets == ntets(m)
    end

    @testset "non-convex surface cannot be silently convex-hull capped" begin
        # This twisted prism's unrestricted Delaunay fill is valid and watertight but
        # has volume 1.116...: it caps the reflex boundary.  The PLC itself has volume
        # sqrt(3)/2 and requires exact boundary recovery.
        ang = deg2rad.((90, 210, 330))
        C = Matrix{Float64}(undef, 3, 6)
        for i in 1:3
            C[:,i] = [cos(ang[i]), sin(ang[i]), 0.0]
            C[:,i+3] = [cos(ang[i] + deg2rad(30)), sin(ang[i] + deg2rad(30)), 1.0]
        end
        faces = NTuple{3,Int32}[(1,3,2), (4,5,6)]
        for i in 1:3
            j = i % 3 + 1
            push!(faces, (i,j,j+3), (i,j+3,i+3))
        end
        tris = reduce(hcat, (Int32[f...] for f in faces))
        m = mesh_volume(Mesh(C; tris=tris); smooth=false)
        @test validate(m).ok
        @test mvpvol(m) ≈ sqrt(3)/2 rtol=1e-9
    end

    @testset "defective surface → explicit blocker (not silent)" begin
        cube = _cube_surface()
        openm = Mesh(cube.coords; tris=cube.tris[:, 1:end-1])   # missing a face
        @test_throws ArgumentError mesh_volume(openm)
        # check=false bypasses the gate (caller takes responsibility)
        @test mesh_volume(openm; check=false) isa Mesh
    end

    @testset "smoothing preserves volume through the pipeline" begin
        m1 = mesh_volume(_cube_surface(); smooth=false)
        m2 = mesh_volume(_cube_surface(); smooth=true, smooth_iters=8)
        @test mvpvol(m1) ≈ mvpvol(m2) rtol=1e-6
    end

    @testset "mesh_planar: 2-D PSLG → quality mesh (validated)" begin
        # unit square domain with an interior point; boundary is a closed loop
        xs = Float64[0,10,10,0]; ys = Float64[0,0,10,10]
        segs = [(1,2),(2,3),(3,4),(4,1)]
        m = mesh_planar(xs, ys, segs; min_angle_deg=25.0, max_area=2.0)
        @test validate(m).ok
        area = sum(triangle_area(node(m,m.tris[1,t]),node(m,m.tris[2,t]),node(m,m.tris[3,t]))
                   for t in 1:ntris(m); init=0.0)
        @test area ≈ 100.0 atol=1e-8            # domain area
        @test all(m.coords[3,:] .== 0.0)        # planar
    end

    @testset "top-level pipeline input, resource, and CRC contracts" begin
        square_segments = [(1,2),(2,3),(3,4),(4,1)]

        planar = mesh_planar([0,1,1,0], [0,0,1,1], square_segments;
                             min_angle_deg=0)
        @test validate(planar).ok
        @test mesh_crc(planar).sha ==
              "850fe31fb8b9c7946d716633cfabdfaf13850456a1b53474d21edfcfa9f194f4"
        @test mesh_crc(mesh_planar(Float64[0,1,1,0], Float64[0,0,1,1],
                                   square_segments; min_angle_deg=0)).sha ==
              mesh_crc(planar).sha
        @test_throws ArgumentError mesh_planar(
            Any[0,1,true,0], [0,0,1,1], square_segments; min_angle_deg=0)
        @test_throws ArgumentError mesh_planar(
            Any[0,1,"1",0], [0,0,1,1], square_segments; min_angle_deg=0)
        @test_throws ArgumentError mesh_planar(
            [0,1,1], [0,0,1,1], square_segments; min_angle_deg=0)
        @test_throws ArgumentError mesh_planar(
            [0,1,1,0], [0,0,1,1],
            Any[(1,2),(2,3),[3,4],(4,1)]; min_angle_deg=0)
        @test_throws ArgumentError mesh_planar(
            [0,1,1,0], [0,0,1,1],
            Any[(1,2),(2,true),(3,4),(4,1)]; min_angle_deg=0)
        @test_throws ArgumentError mesh_planar(
            [0,1,1,0], [0,0,1,1], square_segments; min_angle_deg=true)
        @test_throws ArgumentError mesh_planar(
            [0,1,1,0], [0,0,1,1], square_segments; max_area=-Inf)
        @test_throws ArgumentError mesh_planar(
            [0,1,1,0], [0,0,1,1], square_segments; rng_seed=true)
        @test_throws ArgumentError mesh_planar(
            [0,1,1,0], [0,0,1,1], square_segments; field=1)

        extruded = mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, 0, 1;
            hmax=2, min_angle_deg=0)
        @test validate(extruded).ok
        @test size(extruded.coords,2) == 10
        @test size(extruded.tets,2) == 12
        @test mesh_crc(extruded).sha ==
              "c7783021725d2dfd0b60b83536b5489f556b564af35fe66ef487e0bce15d9e3e"
        @test_throws ArgumentError mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, 0, 1;
            hmax=2, min_angle_deg=0, max_nodes=9)
        @test_throws ArgumentError mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, 0, 1;
            hmax=2, min_angle_deg=0, max_tets=11)
        unread = _UnreadPipelinePoints(4)
        preflight_error = try
            mesh_sized_extrude(unread, unread, square_segments, 0, 1;
                               hmax=2, min_angle_deg=0, max_nodes=5)
            nothing
        catch err
            err
        end
        @test preflight_error isa ArgumentError
        @test preflight_error isa ArgumentError &&
              occursin("max_nodes=5", sprint(showerror, preflight_error))
        @test_throws ArgumentError mesh_sized_extrude(
            unread, unread, square_segments, 0, 1;
            hmax=2, min_angle_deg=0, max_tets=2)
        @test_throws ArgumentError mesh_sized_extrude(
            Any[0,1,true,0], [0,0,1,1], square_segments, 0, 1; hmax=2)
        @test_throws ArgumentError mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, true, 1; hmax=2)
        @test_throws ArgumentError mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, 0, 1; hmax=true)
        @test_throws ArgumentError mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, 0, 1;
            hmax=2, min_angle_deg=true)
        @test_throws ArgumentError mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, 0, 1;
            hmax=2, max_nodes=true)
        @test_throws ArgumentError mesh_sized_extrude(
            [0.0,1e160,1e160,0.0], [0.0,0.0,1e160,1e160],
            square_segments, 0.0, 1.0; hmax=2e160, min_angle_deg=0)

        remote_z0 = 4.097470032826895e-162
        remote_z1 = 2.5149445970871698e185
        remote_hmax = 4.710819546778298e182
        remote = mesh_sized_extrude(
            [0,1,1,0], [0,0,1,1], square_segments, remote_z0, remote_z1;
            hmax=remote_hmax, min_angle_deg=0)
        @test all(==(remote_z0), view(remote.coords,3,1:4))
        @test all(==(remote_z1), view(remote.coords,3,size(remote.coords,2)-3:size(remote.coords,2)))
        @test length(unique(view(remote.coords,3,:))) == 756

        surface = _cube_surface()
        @test_throws ArgumentError mesh_sized(surface; hmax=true)
        @test_throws ArgumentError mesh_sized(surface; field=1)
        @test_throws ArgumentError mesh_sized(surface, 1)
        @test_throws ArgumentError mesh_sized(Mesh(zeros(3,3)); hmax=1)
        @test_throws ArgumentError mesh_volume(surface; smooth=1)
        @test_throws ArgumentError mesh_volume(surface; optimize=1)
        @test_throws ArgumentError mesh_volume(surface; check=1)
        @test_throws ArgumentError mesh_volume(surface; smooth_iters=true)
        @test_throws ArgumentError mesh_volume(surface; smooth_iters=-1)
        @test_throws ArgumentError mesh_volume(surface; rng_seed=true)
    end

    @testset "optimize=true is quality-monotone and preserves volume + validity" begin
        s = cylinder_surface((0.,0,0),(0.,0,1),2.0,5.0; nθ=20, nz=4)
        m0 = mesh_volume(s; optimize=false, smooth=false)
        m1 = mesh_volume(s; optimize=true,  smooth=true)
        @test mvpvol(m1) ≈ mvpvol(m0) rtol=1e-6      # flips+smoothing preserve volume
        @test validate(m1).ok
        # A boundary-recovered mesh can have no movable interior vertices; in that
        # constrained case optimization is correctly a no-op, never a regression.
        @test mesh_quality(m1).n_slivers <= mesh_quality(m0).n_slivers
    end

    @testset "end-to-end: primitive → volume → .msh → read (solver-consumable)" begin
        dir = mktempdir()
        for surf in (box_surface(0,2,0,1,0,1), cylinder_surface((0.,0,0),(0.,0,1),1.0,2.0; nθ=16))
            m = mesh_volume(surf; smooth=false)
            for ver in (2.2, 4.1)
                p = joinpath(dir, "vol_$(ver).msh")
                write_msh(p, m; version=ver)
                back = read_msh(p).mesh
                @test mesh_crc(back).sha == mesh_crc(m).sha    # connectivity preserved
                @test validate(back).ok
            end
        end
        # multi-region tags + physical names survive the round-trip
        mm = tetrahedralize_multi([box_surface(0,1,0,1,0,1), box_surface(1,2,0,1,0,1)])
        p = joinpath(dir, "multi.msh")
        write_msh(p, mm; version=4.1, physical_names=Dict((3,1)=>"regA",(3,2)=>"regB"))
        f = read_msh(p)
        @test mesh_crc(f.mesh).sha == mesh_crc(mm).sha
        @test sort(unique(f.mesh.tet_tag)) == Int32[1,2]
        @test f.physical_names[(3,1)] == "regA" && f.physical_names[(3,2)] == "regB"
    end

    @testset "conforming ENC-COAX mesh is solver-consumable (.msh with physical volumes)" begin
        # The whole point of the acceptance case: a CONFORMING air/case/pin mesh a
        # solver (ASCENT) can read. Build it, write gmsh MSH v4.1 with the three
        # physical volumes, read it back, and confirm connectivity (CRC), region tags,
        # and physical-group names all survive — i.e. the acceptance mesh is usable.
        dir = mktempdir()
        pin  = cylinder_surface((3.,3.,1.5), (0.,0.,1.), 0.7, 3.0; nθ=8, nz=2)   # inside the cavity
        air  = box_surface(1,5,1,5,1,5)
        case = box_surface(0,6,0,6,0,6)
        enc = tetrahedralize_conforming([pin, air, case])            # tags 1=pin, 2=air, 3=case
        @test validate(enc).ok
        @test sort(unique(enc.tet_tag)) == Int32[1,2,3]              # all three volumes present
        p = joinpath(dir, "enclosure.msh")
        names = Dict((3,1)=>"coax_pin", (3,2)=>"air", (3,3)=>"case")
        write_msh(p, enc; version=4.1, physical_names=names)
        f = read_msh(p)
        @test mesh_crc(f.mesh).sha == mesh_crc(enc).sha             # conforming connectivity preserved
        @test sort(unique(f.mesh.tet_tag)) == Int32[1,2,3]
        @test f.physical_names[(3,1)] == "coax_pin" &&
              f.physical_names[(3,2)] == "air" && f.physical_names[(3,3)] == "case"
        @test validate(f.mesh).ok                                   # round-tripped mesh still valid
    end

    @testset "native pipeline end-to-end: CSG → conforming fill → solver-consumable .msh" begin
        # The project's "design → mesh always works" goal, exercised entirely through
        # the native geometry+meshing stack (no OpenCASCADE): build geometry with the
        # native CSG operators, mesh it conformingly at a controlled size, and confirm
        # the result is a validated, ASCENT-consumable gmsh MSH v4.1 with physical groups.
        dir = mktempdir()

        # (A) size-controlled conforming MULTI-REGION enclosure via mesh_box_regions
        #     (native volumetric CSG): case shell (tag 1) around an air cavity (tag 2).
        enc = mesh_box_regions([BoxRegion(0,6,0,6,0,6,1), BoxRegion(1,5,1,5,1,5,2)]; hmax=1.5)
        @test validate(enc).ok
        @test sort(unique(enc.tet_tag)) == Int32[1,2]
        @test mvpvol(enc) ≈ 216.0 rtol=1e-9                          # exact filled volume
        pA = joinpath(dir, "enclosure_csg.msh")
        write_msh(pA, enc; version=4.1, physical_names=Dict((3,1)=>"case", (3,2)=>"air"))
        fA = read_msh(pA)
        @test mesh_crc(fA.mesh).sha == mesh_crc(enc).sha            # connectivity preserved
        @test sort(unique(fA.mesh.tet_tag)) == Int32[1,2]
        @test fA.physical_names[(3,1)] == "case" && fA.physical_names[(3,2)] == "air"
        @test validate(fA.mesh).ok

        # (B) general native CSG → boundary recovery → volume mesh → .msh: a box with a
        #     cylindrical bore (mesh_boolean difference), filled conformingly and written.
        bored = mesh_boolean(box_surface(0,4,0,4,0,4),
                             cylinder_surface((2.,2.,-1.),(0.,0.,1.),1.0,6.0; nθ=12), :difference)
        vol = recover_boundary(bored)
        @test validate(vol).ok
        removed = 0.5*12*1.0^2*sin(2pi/12)*4.0                       # inscribed 12-gon × height 4
        @test mvpvol(vol) ≈ 64.0 - removed rtol=1e-6                 # exact faceted-bore volume, conformed
        pB = joinpath(dir, "bored_block.msh")
        write_msh(pB, vol; version=4.1)
        fB = read_msh(pB)
        @test mesh_crc(fB.mesh).sha == mesh_crc(vol).sha
        @test validate(fB.mesh).ok
    end
end
