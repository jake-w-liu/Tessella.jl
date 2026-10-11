module StepSphereChecks
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

function _with_source(body,fixture,replacements...)
    source=read(joinpath(FIXTURES,fixture),String)
    for (old,new) in replacements
        occursin(old,source) || throw(ArgumentError("fixture rewrite anchor missing: $old"))
        source=replace(source,old=>new)
    end
    return mktemp() do path,io
        write(io,source);close(io);body(path)
    end
end

@testset "closed SPHERICAL_SURFACE face imports as a full sphere" begin
    surface=import_step(joinpath(FIXTURES,"brep_occ_sphere.step");fill=false)
    @test surface isa E.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0 && size(surface.tets,2)==0
    @test all(isapprox(norm(surface.coords[:,i]),2.5;atol=1e-12) for i in axes(surface.coords,2))
    exact=4π*2.5^3/3
    @test _signed_surface_volume(surface)>0.9exact
    @test _signed_surface_volume(surface)<=exact

    filled=import_step(joinpath(FIXTURES,"brep_occ_sphere.step"))
    @test filled isa E.Mesh && validate(filled).ok
    @test size(filled.tets,2)>0
    volumes=_tet_volumes(filled)
    @test all(>(0),volumes)
    @test isapprox(sum(volumes),exact;rtol=0.1)
end

@testset "sphere cover certificate honours placement axis and unit context" begin
    # Flipping the placement axis swaps the computed poles; the vertex set is
    # unchanged so the same boundary must still certify.
    _with_source("brep_occ_sphere.step",
                 "#34 = DIRECTION('',(0.,0.,1.));"=>"#34 = DIRECTION('',(0.,0.,-1.));") do path
        surface=import_step(path;fill=false)
        @test validate(surface).ok
        @test all(isapprox(norm(surface.coords[:,i]),2.5;atol=1e-12) for i in axes(surface.coords,2))
    end

    # A metre-context sphere surface scales to millimetres like every other root.
    _with_source("brep_occ_sphere.step",
                 "SI_UNIT(.MILLI.,.METRE.)"=>"SI_UNIT(\$,.METRE.)") do path
        surface=import_step(path;fill=false)
        @test validate(surface).ok
        @test all(isapprox(norm(surface.coords[:,i]),2500.0;atol=1e-9) for i in axes(surface.coords,2))
    end
end

@testset "trimmed or malformed spherical topology rejects explicitly" begin
    # A boundary vertex off both poles is a genuine trim (cap/equator).
    _with_source("brep_occ_sphere.step",
                 "#26 = CARTESIAN_POINT('',(0.,0.,2.5));"=>"#26 = CARTESIAN_POINT('',(0.,0.,2.4));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("trimmed SPHERICAL_SURFACE",sprint(showerror,err))
    end

    # Poles must both appear — a ring that only touches the south pole does not
    # cover the sphere.
    _with_source("brep_occ_sphere.step",
                 "#26 = CARTESIAN_POINT('',(0.,0.,2.5));"=>"#26 = CARTESIAN_POINT('',(0.,0.,-2.5));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("trimmed SPHERICAL_SURFACE",sprint(showerror,err))
    end

    # A non-positive radius is malformed, not unsupported.
    _with_source("brep_occ_sphere.step",
                 "#31 = SPHERICAL_SURFACE('',#32,2.5);"=>"#31 = SPHERICAL_SURFACE('',#32,-2.5);") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("radius must be positive",sprint(showerror,err))
    end

    # Two spherical faces in one shell exceed the closed-face special case.
    _with_source("brep_occ_sphere.step",
                 "#16 = CLOSED_SHELL('',(#17));"=>"#16 = CLOSED_SHELL('',(#17,#17));") do path
        err=try import_step(path;fill=false); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("unsupported BREP face or edge topology",sprint(showerror,err))
    end
end

@testset "primitive SPHERE path is unchanged" begin
    # The classified SPHERE primitive (no BRep topology) keeps its own path.
    primitive=import_step(joinpath(FIXTURES,"sphere.step");fill=false)
    @test primitive isa E.Mesh && validate(primitive).ok
end

end # module
