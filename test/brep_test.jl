using Test
using Tessella
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: nnodes, ntris, ntets, tet_volume, node, validate

const FIXTURES = joinpath(@__DIR__, "fixtures")

gvol(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),
                         node(m,m.tets[3,t]),node(m,m.tets[4,t]))
              for t in 1:ntets(m); init=0.0)

function signed_surface_volume(s)
    v = 0.0
    for k in 1:ntris(s)
        a = node(s, s.tris[1,k]); b = node(s, s.tris[2,k]); c = node(s, s.tris[3,k])
        cx = a[2]*b[3]-a[3]*b[2]; cy = a[3]*b[1]-a[1]*b[3]; cz = a[1]*b[2]-a[2]*b[1]
        v += (cx*c[1] + cy*c[2] + cz*c[3])
    end
    v/6
end

@testset "native STEP/IGES classified-solid import" begin
    @testset "ISO-10303-21 axis-aligned block" begin
        path=joinpath(FIXTURES,"unit_box.step")
        ents=parse_step_entities(read(path,String))
        @test count(e->e.kind=="CARTESIAN_POINT", values(ents))==8
        surface=import_step(path; fill=false)
        @test validate(surface).ok
        @test nnodes(surface)==8
        @test ntris(surface)==12
        @test signed_surface_volume(surface)≈2.0 atol=1e-12
        mesh=import_step(path)
        @test validate(mesh).ok
        @test ntets(mesh)>0
        @test gvol(mesh)≈2.0 atol=1e-12
    end

    @testset "ISO-10303-21 SPHERE and cylinder" begin
        sph=import_step(joinpath(FIXTURES,"sphere.step"))
        @test validate(sph).ok
        @test ntets(sph)>0
        oracle=Tessella.Geometry.sphere_surface((0.5,0.25,0.0),1.5)
        @test gvol(sph)≈signed_surface_volume(oracle) atol=1e-12

        cyl=import_step(joinpath(FIXTURES,"cylinder.step"))
        @test validate(cyl).ok
        @test ntets(cyl)>0
        nθ=24
        @test gvol(cyl)≈0.5*nθ*1.0^2*sin(2π/nθ)*2.0 atol=1e-12
        @test gvol(cyl)≈signed_surface_volume(
            Tessella.Geometry.cylinder_surface((0.,0.,0.),(0.,0.,1.),1.0,2.0)) atol=1e-12
    end

    @testset "IGES type 150 block and 158 sphere" begin
        box=import_iges(joinpath(FIXTURES,"unit_box.iges"))
        @test validate(box).ok
        @test ntets(box)>0
        @test gvol(box)≈2.0 atol=1e-12

        sph=import_iges(joinpath(FIXTURES,"sphere.iges"))
        @test validate(sph).ok
        @test ntets(sph)>0
        @test gvol(sph)≈signed_surface_volume(
            Tessella.Geometry.sphere_surface((0.5,0.25,0.0),1.5)) atol=1e-12
    end

    @testset "explicit blockers" begin
        @test_throws ArgumentError import_step("/no/such/part.step")
        @test_throws ArgumentError import_iges("/no/such/part.iges")
        err=try
            import_step(joinpath(FIXTURES,"unsupported.step")); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("MANIFOLD_SOLID_BREP", sprint(showerror, err))
        err=try
            import_iges(joinpath(FIXTURES,"unsupported.iges")); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("100", sprint(showerror, err))
        junk=mktemp() do path,io
            write(io,"not a cad file\n"); close(io); path
        end
        @test_throws ArgumentError import_step(junk)
    end
end
