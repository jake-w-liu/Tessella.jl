module StepCylinderChecks
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

function _with_source(body,replacements...)
    source=read(joinpath(FIXTURES,"brep_occ_cylinder.step"),String)
    for (old,new) in replacements
        occursin(old,source) || throw(ArgumentError("fixture rewrite anchor missing: $old"))
        source=replace(source,old=>new)
    end
    return mktemp() do path,io
        write(io,source);close(io);body(path)
    end
end

@testset "rim-certified CYLINDRICAL_SURFACE shell imports as a full cylinder" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_cylinder.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0 && size(surface.tets,2)==0
    @test vec(minimum(surface.coords,dims=2))≈[-1.,-1.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[1.,1.,4.]
    # Lateral band nodes sit on the R=1 wall; cap fans keep every node inside.
    @test all(norm(surface.coords[1:2,i])<=1.0+1e-12 for i in axes(surface.coords,2))
    @test any(isapprox(norm(surface.coords[1:2,i]),1.0;atol=1e-12) for i in axes(surface.coords,2))
    exact=π*4.0
    @test _signed_surface_volume(surface)>0.9exact
    @test _signed_surface_volume(surface)<=exact

    filled=import_step(joinpath(FIXTURES,"brep_occ_cylinder.step"))
    @test filled isa E.Mesh && validate(filled).ok
    @test size(filled.tets,2)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test isapprox(sum(volumes),exact;rtol=0.1)
end

@testset "cylinder certificate honours unit context" begin
    _with_source("SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(\$,.METRE.)") do path
        surface=import_step(path;fill=false)
        @test validate(surface).ok
        @test vec(minimum(surface.coords,dims=2))≈[-1000.,-1000.,0.] &&
              vec(maximum(surface.coords,dims=2))≈[1000.,1000.,4000.]
    end
end

@testset "trimmed or malformed cylindrical topology rejects explicitly" begin
    # A cap bound containing a partial arc is a genuine trim.
    _with_source("#52 = EDGE_LOOP('',(#55));"=>
        "#52 = EDGE_LOOP('',(#55,#90,#91));"*
        "#90 = ORIENTED_EDGE('',*,*,#92,.T.);"*
        "#91 = ORIENTED_EDGE('',*,*,#93,.T.);"*
        "#92 = EDGE_CURVE('',#23,#94,#34,.T.);"*
        "#93 = EDGE_CURVE('',#94,#23,#34,.T.);"*
        "#94 = VERTEX_POINT('',#95);"*
        "#95 = CARTESIAN_POINT('',(-1.,0.,4.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("partial arc",sprint(showerror,err))
    end

    # A rim circle whose radius disagrees with the surface is malformed.
    _with_source("#36 = CIRCLE('',#44,1.);"=>"#36 = CIRCLE('',#44,1.5);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("radius does not match",sprint(showerror,err))
    end

    # A rim vertex pulled off its circle is a geometric inconsistency.
    _with_source("#33 = CARTESIAN_POINT('',(1.,0.,4.));"=>"#33 = CARTESIAN_POINT('',(1.,0.,3.5));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("off the circle",sprint(showerror,err))
    end

    # A rim circle centred off the cylinder axis is malformed.
    _with_source("#37 = CARTESIAN_POINT('',(0.,0.,4.));"=>"#37 = CARTESIAN_POINT('',(1.,0.,4.));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("off the cylinder axis",sprint(showerror,err))
    end

    # Both caps on the bottom rim leave the top open.
    _with_source("#55 = ORIENTED_EDGE('',*,*,#25,.T.);"=>"#55 = ORIENTED_EDGE('',*,*,#29,.T.);",
                 "#61 = AXIS2_PLACEMENT_3D('',#53,#47,#42);"=>"#61 = AXIS2_PLACEMENT_3D('',#45,#47,#42);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("do not close both rims",sprint(showerror,err))
    end

    # A cap plane that does not contain its rim circle is malformed.
    _with_source("#61 = AXIS2_PLACEMENT_3D('',#53,#47,#42);"=>"#61 = AXIS2_PLACEMENT_3D('',#32,#47,#42);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("does not contain its rim circle",sprint(showerror,err))
    end

    # An extra shell face is a different solid, not a cylinder certificate.
    _with_source("#16 = CLOSED_SHELL('',(#17,#40,#50));"=>"#16 = CLOSED_SHELL('',(#17,#40,#50,#40));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("degenerate",sprint(showerror,err))
    end

    # A lone cylindrical face stays an unsupported-topology blocker.
    _with_source("#16 = CLOSED_SHELL('',(#17,#40,#50));"=>"#16 = CLOSED_SHELL('',(#17));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("unsupported BREP face or edge topology",sprint(showerror,err))
    end
end

@testset "primitive RIGHT_CIRCULAR_CYLINDER path is unchanged" begin
    primitive=import_step(joinpath(FIXTURES,"cylinder.step");fill=false)
    @test primitive isa E.Mesh && validate(primitive).ok
end

end # module
