using Test
using SHA
using Tessella
using Tessella.BRep: parse_step_entities
using Tessella.MeshTypes: nnodes, ntris, ntets, tet_volume, node, validate
using Tessella.NURBS: NURBSCurve, NURBSSurface, nurbs_eval

const FIXTURES = joinpath(@__DIR__, "..", "fixtures")

_brep_volume(m) = sum(tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),
                                 node(m,m.tets[3,t]),node(m,m.tets[4,t]))
                      for t in 1:ntets(m); init=0.0)

function _brep_signed_surface_volume(s)
    v = 0.0
    for k in 1:ntris(s)
        a = node(s, s.tris[1,k]); b = node(s, s.tris[2,k]); c = node(s, s.tris[3,k])
        cx = a[2]*b[3]-a[3]*b[2]; cy = a[3]*b[1]-a[1]*b[3]; cz = a[1]*b[2]-a[2]*b[1]
        v += (cx*c[1] + cy*c[2] + cz*c[3])
    end
    v/6
end

function _brep_iges_parameter_line(payload; pointer=1, sequence=1)
    ncodeunits(payload)<=64 || throw(ArgumentError("test IGES payload exceeds 64 columns"))
    return rpad(payload,64)*lpad(pointer,8)*"P"*lpad(sequence,7)*"\n"
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
        @test _brep_signed_surface_volume(surface)≈2.0 atol=1e-12
        mesh=import_step(path)
        @test validate(mesh).ok
        @test ntets(mesh)>0
        @test _brep_volume(mesh)≈2.0 atol=1e-12
    end

    @testset "ISO-10303-21 SPHERE and cylinder" begin
        sph=import_step(joinpath(FIXTURES,"sphere.step"))
        @test validate(sph).ok
        @test ntets(sph)>0
        oracle=Tessella.Geometry.sphere_surface((0.5,0.25,0.0),1.5)
        @test _brep_volume(sph)≈_brep_signed_surface_volume(oracle) atol=1e-12

        cyl=import_step(joinpath(FIXTURES,"cylinder.step"))
        @test validate(cyl).ok
        @test ntets(cyl)>0
        nθ=24
        @test _brep_volume(cyl)≈0.5*nθ*1.0^2*sin(2π/nθ)*2.0 atol=1e-12
        @test _brep_volume(cyl)≈_brep_signed_surface_volume(
            Tessella.Geometry.cylinder_surface((0.,0.,0.),(0.,0.,1.),1.0,2.0)) atol=1e-12
    end

    @testset "IGES type 150 block and 158 sphere" begin
        box=import_iges(joinpath(FIXTURES,"unit_box.iges"))
        @test validate(box).ok
        @test ntets(box)>0
        @test _brep_volume(box)≈2.0 atol=1e-12

        sph=import_iges(joinpath(FIXTURES,"sphere.iges"))
        @test validate(sph).ok
        @test ntets(sph)>0
        @test _brep_volume(sph)≈_brep_signed_surface_volume(
            Tessella.Geometry.sphere_surface((0.5,0.25,0.0),1.5)) atol=1e-12
    end

    @testset "classified cone solids" begin
        nθ=24; h=2.0; r1=1.0; r2=0.4
        oracle=Tessella.Geometry.cone_surface((0.,0.,0.),(0.,0.,1.),r1,r2,h)
        expected=_brep_signed_surface_volume(oracle)
        polygon_factor=0.5*nθ*sinpi(2/nθ)
        @test expected≈polygon_factor*h*(r1^2+r1*r2+r2^2)/3 rtol=2e-14

        step=import_step(joinpath(FIXTURES,"cone.step"))
        @test validate(step).ok
        @test ntets(step)>0
        @test _brep_volume(step)≈expected atol=1e-12

        iges=import_iges(joinpath(FIXTURES,"cone.iges"))
        @test validate(iges).ok
        @test ntets(iges)>0
        @test _brep_volume(iges)≈expected atol=1e-12

        cyl=import_iges(joinpath(FIXTURES,"cylinder.iges"))
        @test validate(cyl).ok
        @test _brep_volume(cyl)≈0.5*nθ*1.0^2*sin(2π/nθ)*2.0 atol=1e-12
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
            lines=readlines(p)
            @test all(ncodeunits(line)==80 for line in lines)
            @test count(line->line[73]=='D',lines)==6
            @test count(line->line[73]=='P',lines)>=3
            @test bytes2hex(sha256(read(p)))==
                  "ae515df934189f3d0b3cf5614bd39427cd6d015ceb07a9858c44fc32ea7986f5"
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

        complex_surface_src="""
ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('complex rational surface'),'2;1');
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=CARTESIAN_POINT('',(0.,1.,0.));
#3=CARTESIAN_POINT('',(1.,0.,0.));
#4=CARTESIAN_POINT('',(1.,1.,2.));
#5=( BOUNDED_SURFACE()
     B_SPLINE_SURFACE(1,1,((#1,#2),(#3,#4)),.UNSPECIFIED.,.F.,.F.,.F.)
     B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)
     GEOMETRIC_REPRESENTATION_ITEM()
     RATIONAL_B_SPLINE_SURFACE(((1.,2.),(3.,4.)))
     REPRESENTATION_ITEM('') SURFACE() );
ENDSEC;
END-ISO-10303-21;
"""
        complex_surface=mktemp() do p,io
            write(io,complex_surface_src); close(io)
            only(import_nurbs_step(p))
        end
        @test complex_surface isa NURBSSurface
        @test complex_surface.weights==[1.0 2.0; 3.0 4.0]
        rational_mid=nurbs_eval(complex_surface,0.5,0.5)
        @test all(isapprox(rational_mid[i],(0.7,0.6,0.8)[i];atol=1e-14)
                  for i in 1:3)
    end

    @testset "strict parser, resource, and atomic-write contracts" begin
        duplicate="""
ISO-10303-21;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#1=CARTESIAN_POINT('',(1.,0.,0.));
ENDSEC;
END-ISO-10303-21;
"""
        @test_throws ArgumentError parse_step_entities(duplicate)
        nonfinite=replace(duplicate,
                          "#1=CARTESIAN_POINT('',(1.,0.,0.));"=>
                          "#2=CARTESIAN_POINT('',(1e999,0.,0.));")
        @test_throws ArgumentError parse_step_entities(nonfinite)
        oversized_id=replace(duplicate,
                             "#1=CARTESIAN_POINT('',(1.,0.,0.));"=>
                             "#999999999999999999999999=CARTESIAN_POINT('',(1.,0.,0.));")
        @test_throws ArgumentError parse_step_entities(oversized_id)

        corners=join(("#$(i)=CARTESIAN_POINT('',($(p[1]),$(p[2]),$(p[3])));"
                      for (i,p) in enumerate(((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                                               (0,0,1),(1,0,1),(1,1,1),(0,1,1)))),"\n")
        conflicting="ISO-10303-21;\nDATA;\n"*corners*
                    "\n#9=MANIFOLD_SOLID_BREP('',#1);\nENDSEC;\nEND-ISO-10303-21;\n"
        err=mktemp() do path,io
            write(io,conflicting); close(io)
            try import_step(path); nothing catch e e end
        end
        @test err isa ArgumentError
        @test occursin("MANIFOLD_SOLID_BREP",sprint(showerror,err))

        multiple_step="""
ISO-10303-21;
DATA;
#1=SPHERE('',\$,1.);
#2=SPHERE('',\$,2.);
ENDSEC;
END-ISO-10303-21;
"""
        multiple_step_error=mktemp() do path,io
            write(io,multiple_step); close(io)
            try import_step(path); nothing catch e e end
        end
        @test multiple_step_error isa ArgumentError
        @test occursin("multiple recognized solids",sprint(showerror,multiple_step_error))

        standard=_brep_iges_parameter_line("150,2D0,1D0,1D0,0D0,0D0,0D0;";
                                            pointer=12345678)
        standard_surface=mktemp() do path,io
            write(io,standard); close(io)
            import_iges(path; fill=false)
        end
        @test validate(standard_surface).ok
        @test _brep_signed_surface_volume(standard_surface)≈2.0 atol=1e-12
        malformed=_brep_iges_parameter_line("150,2.,bad,1.,0.,0.,0.;")
        mktemp() do path,io
            write(io,malformed); close(io)
            @test_throws ArgumentError import_iges(path)
        end
        mktemp() do path,io
            write(io,_brep_iges_parameter_line("150,2.,1.,1.,0.,0.,0.")); close(io)
            @test_throws ArgumentError import_iges(path)
        end
        multiple_iges=_brep_iges_parameter_line("158,1.,0.,0.,0.;";sequence=1)*
                      _brep_iges_parameter_line("158,2.,0.,0.,0.;";sequence=2)
        multiple_iges_error=mktemp() do path,io
            write(io,multiple_iges); close(io)
            try import_iges(path); nothing catch e e end
        end
        @test multiple_iges_error isa ArgumentError
        @test occursin("multiple recognized solids",sprint(showerror,multiple_iges_error))

        bad_count=_brep_iges_parameter_line("126,1000001,1,0,0,0,0,0;")
        mktemp() do path,io
            write(io,bad_count); close(io)
            @test_throws ArgumentError import_nurbs_iges(path)
        end
        @test_throws ArgumentError Tessella.BRep._expand_knots([1_000_001],[0.0],
                                                               "test")
        @test_throws ArgumentError Tessella.BRep._as_int(true,"test","count")

        curve=NURBSCurve(1,[0,0,1,1],[(0.0,0.0,0.0),(1.0,0.0,0.0)])
        mktemp() do path,io
            write(io,"preserve-me"); close(io)
            @test_throws ArgumentError export_iges_nurbs(path,[curve,42])
            @test read(path,String)=="preserve-me"
            curve.knots[1]=NaN
            @test_throws ArgumentError export_iges_nurbs(path,[curve])
            @test read(path,String)=="preserve-me"
        end
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
        mktemp() do path,io
            write(io,"not a cad file\n"); close(io)
            @test_throws ArgumentError import_step(path)
        end
    end
end
