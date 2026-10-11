module StepTorusSegmentChecks
using Test,Tessella,LinearAlgebra
const E=Tessella.Elements
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function _signed_surface_volume(mesh)
    volume=0.0
    for k in 1:size(mesh.tris,2)
        a,b,c=(mesh.coords[:,mesh.tris[i,k]] for i in 1:3)
        volume+=dot(a,cross(b,c))/6
    end
    return volume
end

function _with_source(body,replacements...)
    source=read(joinpath(FIXTURES,"brep_occ_torus_segment.step"),String)
    for (old,new) in replacements
        occursin(old,source) || throw(ArgumentError("fixture rewrite anchor missing: $old"))
        source=replace(source,old=>new)
    end
    return mktemp() do path,io
        write(io,source);close(io);body(path)
    end
end

@testset "rim-certified TOROIDAL_SURFACE shell imports as a torus segment" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_torus_segment.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0 && size(surface.tets,2)==0
    lo=vec(minimum(surface.coords,dims=2)); hi=vec(maximum(surface.coords,dims=2))
    @test lo≈[0.,0.,-1.] && hi≈[4.,4.,1.]
    # Lateral nodes satisfy the tube distance law; the two cap fan centres
    # lie on the tube centre circle instead.
    @test count(i->abs(sqrt((hypot(surface.coords[1,i],surface.coords[2,i])-3)^2+
              surface.coords[3,i]^2)-1)<=1e-9,axes(surface.coords,2))==
              size(surface.coords,2)-2
    # Quarter segment volume: (sweep/2pi)*2*pi^2*R*r^2 = 1.5*pi^2.
    exact=1.5π^2
    @test _signed_surface_volume(surface)>0.9exact
    @test _signed_surface_volume(surface)<=exact
end

@testset "torus segment certificate honours unit context" begin
    _with_source("SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(\$,.METRE.)") do path
        surface=import_step(path;fill=false)
        @test validate(surface).ok
        @test vec(maximum(surface.coords,dims=2))≈[4000.,4000.,1000.]
    end
end

@testset "trimmed or malformed toroidal topology rejects explicitly" begin
    # Lateral rim circles split into arcs are genuine trims.
    _with_source("#19 = EDGE_LOOP('',(#20,#24,#26,#28));"=>
        "#19 = EDGE_LOOP('',(#20,#90,#91,#26,#28));"*
        "#90 = ORIENTED_EDGE('',*,*,#92,.T.);"*
        "#91 = ORIENTED_EDGE('',*,*,#93,.T.);"*
        "#92 = EDGE_CURVE('',#23,#94,#34,.T.);"*
        "#93 = EDGE_CURVE('',#94,#23,#34,.T.);"*
        "#94 = VERTEX_POINT('',#95);"*
        "#95 = CARTESIAN_POINT('',(0.,4.,1.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("surface parallel",sprint(showerror,err))
    end

    # A rim circle whose radius disagrees with the minor radius is malformed.
    _with_source("#36 = CIRCLE('',#44,1.);"=>"#36 = CIRCLE('',#44,1.5);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("does not match the toroidal surface",sprint(showerror,err))
    end

    # A rim circle centred off the tube circle is malformed.
    _with_source("#46 = CARTESIAN_POINT('',(3.,0.,0.));"=>"#46 = CARTESIAN_POINT('',(2.5,0.,0.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("off the tube circle",sprint(showerror,err))
    end

    # A rim circle in a non-meridian plane is a genuine trim.
    _with_source("#47 = DIRECTION('',(0.,1.,0.));"=>"#47 = DIRECTION('',(0.,0.70710678118654752,0.70710678118654752));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("meridian plane",sprint(showerror,err))
    end

    # A meridian-normal rim circle whose plane misses the axis is malformed.
    _with_source("#46 = CARTESIAN_POINT('',(3.,0.,0.));"=>"#46 = CARTESIAN_POINT('',(0.,3.,0.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("does not contain the torus axis",sprint(showerror,err))
    end

    # A seam circle lifted off the surface parallels is a genuine trim.
    _with_source("#39 = AXIS2_PLACEMENT_3D('',#12,#35,#42);"=>"#39 = AXIS2_PLACEMENT_3D('',#90,#35,#42);".*
        "#90 = CARTESIAN_POINT('',(0.,0.,1.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("off the toroidal surface",sprint(showerror,err))
    end

    # A cap plane that is not its rim's meridian plane is malformed.
    _with_source("#57 = DIRECTION('',(0.,-1.,0.));"=>"#57 = DIRECTION('',(0.,0.,-1.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("not a meridian plane",sprint(showerror,err))
    end

    # Flipping one cap's orientation cannot certify either segment side.
    _with_source("#64 = DIRECTION('',(-1.,0.,0.));"=>"#64 = DIRECTION('',(1.,0.,0.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("orientations are inconsistent",sprint(showerror,err))
    end

    # An extra shell face is a different solid, not a segment certificate.
    _with_source("#16 = CLOSED_SHELL('',(#17,#40,#50));"=>"#16 = CLOSED_SHELL('',(#17,#40,#50,#40));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("degenerate",sprint(showerror,err))
    end
end

@testset "single-face full torus path is unchanged" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_torus.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
end

end # module
