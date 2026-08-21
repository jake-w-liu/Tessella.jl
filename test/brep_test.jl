using Test
using Tessella
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: nnodes, ntris, ntets, tet_volume, node, validate
using Tessella.NURBS: NURBSCurve, NURBSSurface, nurbs_eval

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

    @testset "classified cone solids" begin
        nθ=24; h=2.0; r1=1.0; r2=0.4
        oracle=Tessella.Geometry.cone_surface((0.,0.,0.),(0.,0.,1.),r1,r2,h)
        expected=signed_surface_volume(oracle)
        polygon_factor=0.5*nθ*sinpi(2/nθ)
        @test expected≈polygon_factor*h*(r1^2+r1*r2+r2^2)/3 rtol=2e-14

        step=import_step(joinpath(FIXTURES,"cone.step"))
        @test validate(step).ok
        @test ntets(step)>0
        @test gvol(step)≈expected atol=1e-12

        iges=import_iges(joinpath(FIXTURES,"cone.iges"))
        @test validate(iges).ok
        @test ntets(iges)>0
        @test gvol(iges)≈expected atol=1e-12

        cyl=import_iges(joinpath(FIXTURES,"cylinder.iges"))
        @test validate(cyl).ok
        @test gvol(cyl)≈0.5*nθ*1.0^2*sin(2π/nθ)*2.0 atol=1e-12
    end

    @testset "STEP/IGES NURBS import and IGES round-trip" begin
        curves=import_nurbs_step(joinpath(FIXTURES,"bezier.step"))
        @test length(curves)==1
        bezier=only(curves)
        @test bezier isa NURBSCurve
        P=[(0.0,0.0,0.0),(1.0,1.0,0.0),(2.0,0.0,0.0)]
        bernstein(t)=((1-t)^2 .* P[1] .+ 2(1-t)*t .* P[2] .+ t^2 .* P[3])
        for t in (0.0,0.25,0.5,0.75,1.0)
            q=nurbs_eval(bezier,t); b=bernstein(t)
            @test hypot(q[1]-b[1],q[2]-b[2],q[3]-b[3])<=1e-14
        end

        patches=import_nurbs_step(joinpath(FIXTURES,"bilinear.step"))
        @test length(patches)==1
        patch=only(patches)
        @test patch isa NURBSSurface
        q=nurbs_eval(patch,0.5,0.5)
        @test q[1]≈0.5 && q[2]≈0.5 && q[3]≈0.5

        w=1/sqrt(2)
        ctrls=[(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0),(-1.0,1.0,0.0),
               (-1.0,0.0,0.0),(-1.0,-1.0,0.0),(0.0,-1.0,0.0),(1.0,-1.0,0.0),(1.0,0.0,0.0)]
        knots=[0,0,0,1,1,2,2,3,3,4,4,4]./4
        weights=[1,w,1,w,1,w,1,w,1]
        circle=NURBSCurve(2,knots,ctrls,weights)
        imported=mktemp() do p,io
            close(io)
            export_iges_nurbs(p, [circle, NURBSSurface(1,1,[0,0,1,1],[0,0,1,1],
                [(0.0,0.0,0.0) (0.0,1.0,0.0); (1.0,0.0,0.0) (1.0,1.0,2.0)])])
            import_nurbs_iges(p)
        end
        @test length(imported)==2
        ic=imported[1]; isurf=imported[2]
        @test ic isa NURBSCurve && isurf isa NURBSSurface
        for (u,target) in ((0.0,(1.0,0.0,0.0)),(0.25,(0.0,1.0,0.0)),
                           (0.5,(-1.0,0.0,0.0)),(0.75,(0.0,-1.0,0.0)),(1.0,(1.0,0.0,0.0)))
            q=nurbs_eval(ic,u)
            @test hypot(q[1]-target[1],q[2]-target[2],q[3]-target[3])<=1e-12
        end
        q=nurbs_eval(isurf,0.5,0.5)
        @test q[1]≈0.5 && q[2]≈0.5 && q[3]≈0.5

        complex_src="""
ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('complex NURBS'),'2;1');
FILE_NAME('complex.step','2026-08-21T00:00:00',('Tessella'),('Tessella'),
          'Tessella.jl','Tessella.jl','');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(1.,1.,0.));
#3 = CARTESIAN_POINT('',(2.,0.,0.));
#4 = ( BOUNDED_CURVE() B_SPLINE_CURVE(2,(#1,#2,#3),.UNSPECIFIED.,.F.,.F.)
       B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)
       CURVE() GEOMETRIC_REPRESENTATION_ITEM()
       RATIONAL_B_SPLINE_CURVE((1.,1.,1.))
       REPRESENTATION_ITEM('') );
ENDSEC;
END-ISO-10303-21;
"""
        cc=mktemp() do p,io
            write(io,complex_src); close(io)
            only(import_nurbs_step(p))
        end
        q=nurbs_eval(cc,0.5); b=bernstein(0.5)
        @test hypot(q[1]-b[1],q[2]-b[2],q[3]-b[3])<=1e-14
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
