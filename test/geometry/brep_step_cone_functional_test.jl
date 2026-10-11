module StepConicalChecks
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

function _tet_volumes(mesh)
    return [det(hcat((mesh.coords[:,cell[i]]-mesh.coords[:,cell[1]] for i in 2:4)...))/6
        for cell in eachcol(mesh.tets)]
end

function _with_source(body,replacements...;fixture="brep_occ_frustum.step")
    source=read(joinpath(FIXTURES,fixture),String)
    for (old,new) in replacements
        occursin(old,source) || throw(ArgumentError("fixture rewrite anchor missing: $old"))
        source=replace(source,old=>new)
    end
    return mktemp() do path,io
        write(io,source);close(io);body(path)
    end
end

@testset "rim-certified CONICAL_SURFACE shell imports as a full frustum" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_frustum.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0 && size(surface.tets,2)==0
    @test vec(minimum(surface.coords,dims=2))≈[-2.,-2.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[2.,2.,3.]
    # Frustum volume pi*h*(r1^2+r1*r2+r2^2)/3 with h=3, r1=2, r2=1.
    exact=7π
    @test _signed_surface_volume(surface)>0.9exact
    @test _signed_surface_volume(surface)<=exact

    filled=import_step(joinpath(FIXTURES,"brep_occ_frustum.step"))
    @test filled isa E.Mesh && validate(filled).ok
    @test size(filled.tets,2)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test isapprox(sum(volumes),exact;rtol=0.1)
end

@testset "apex-certified CONICAL_SURFACE shell imports as a full cone" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_cone.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0 && size(surface.tets,2)==0
    @test vec(minimum(surface.coords,dims=2))≈[-2.,-2.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[2.,2.,4.]
    @test any(all(isapprox.(surface.coords[:,i],[0.,0.,0.];atol=1e-12))
              for i in axes(surface.coords,2))
    # Full cone volume pi*h*r^2/3 with h=4, r=2.
    exact=16π/3
    @test _signed_surface_volume(surface)>0.9exact
    @test _signed_surface_volume(surface)<=exact

    filled=import_step(joinpath(FIXTURES,"brep_occ_cone.step"))
    @test filled isa E.Mesh && validate(filled).ok
    @test size(filled.tets,2)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test isapprox(sum(volumes),exact;rtol=0.1)
end

@testset "conical certificate honours unit context" begin
    _with_source("SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(\$,.METRE.)") do path
        surface=import_step(path;fill=false)
        @test validate(surface).ok
        @test vec(minimum(surface.coords,dims=2))≈[-2000.,-2000.,0.] &&
              vec(maximum(surface.coords,dims=2))≈[2000.,2000.,3000.]
    end
    _with_source("SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(\$,.METRE.)";
                 fixture="brep_occ_cone.step") do path
        surface=import_step(path;fill=false)
        @test validate(surface).ok
        @test vec(maximum(surface.coords,dims=2))≈[2000.,2000.,4000.]
    end
end

@testset "trimmed or malformed conical topology rejects explicitly" begin
    # Lateral rim circles split into arcs are genuine trims.
    _with_source("#19 = EDGE_LOOP('',(#20,#24,#26,#28));"=>
        "#19 = EDGE_LOOP('',(#20,#90,#91,#26,#28));"*
        "#90 = ORIENTED_EDGE('',*,*,#92,.T.);"*
        "#91 = ORIENTED_EDGE('',*,*,#93,.T.);"*
        "#92 = EDGE_CURVE('',#23,#94,#34,.T.);"*
        "#93 = EDGE_CURVE('',#94,#23,#34,.T.);"*
        "#94 = VERTEX_POINT('',#95);"*
        "#95 = CARTESIAN_POINT('',(-1.,0.,3.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("partial arc",sprint(showerror,err))
    end

    # A rim circle whose radius disagrees with the semi-angle is malformed.
    _with_source("#36 = CIRCLE('',#44,2.);"=>"#36 = CIRCLE('',#44,1.5);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("does not match the conical surface",sprint(showerror,err))
    end

    # A rim vertex pulled along the generator leaves the circle.
    _with_source("#32 = CARTESIAN_POINT('',(2.,0.,0.));"=>"#32 = CARTESIAN_POINT('',(1.5,0.,1.5));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("off the circle",sprint(showerror,err))
    end

    # A rim circle centred off the cone axis is malformed.
    _with_source("#45 = CARTESIAN_POINT('',(0.,0.,0.));"=>"#45 = CARTESIAN_POINT('',(0.5,0.,0.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("off the cone axis",sprint(showerror,err))
    end

    # A seam whose carrier line misses the apex is not a cone generator.
    _with_source("#33 = CARTESIAN_POINT('',(1.,0.,3.));"=>"#33 = CARTESIAN_POINT('',(1.,0.,4.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("not a cone generator",sprint(showerror,err))
    end

    # A cap plane that does not contain its rim circle is malformed.
    _with_source("#61 = AXIS2_PLACEMENT_3D('',#53,#47,#42);"=>"#61 = AXIS2_PLACEMENT_3D('',#46,#47,#42);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("does not contain its rim circle",sprint(showerror,err))
    end

    # Dropping a cap leaves a two-rim shell one cap short.
    _with_source("#16 = CLOSED_SHELL('',(#17,#40,#50));"=>"#16 = CLOSED_SHELL('',(#17,#40));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("rim(s) but",sprint(showerror,err))
    end

    # A conical face bounded by one rim without its apex is a genuine trim.
    _with_source("#19 = EDGE_LOOP('',(#20,#24,#26,#28));"=>"#19 = EDGE_LOOP('',(#20,#26,#28));";
                 fixture="brep_occ_cone.step") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("full cone or frustum",sprint(showerror,err))
    end

    # A degenerate semi-angle cannot certify a cone.
    _with_source(",-0.32175055439664219);"=>",0.);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("semi-angle is degenerate",sprint(showerror,err))
    end

    # An extra shell face is a different solid, not a cone certificate.
    _with_source("#16 = CLOSED_SHELL('',(#17,#40,#50));"=>"#16 = CLOSED_SHELL('',(#17,#40,#50,#40));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("degenerate",sprint(showerror,err))
    end

    # A lone conical face stays an unsupported-topology blocker.
    _with_source("#16 = CLOSED_SHELL('',(#17,#40,#50));"=>"#16 = CLOSED_SHELL('',(#17));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("unsupported BREP face or edge topology",sprint(showerror,err))
    end
end

@testset "primitive RIGHT_CIRCULAR_CONE path is unchanged" begin
    primitive=import_step(joinpath(FIXTURES,"cone.step");fill=false)
    @test primitive isa E.Mesh && validate(primitive).ok
end

end # module
