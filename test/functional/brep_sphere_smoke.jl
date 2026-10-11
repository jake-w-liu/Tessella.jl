module BRepSphereSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_sphere.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0 &&
        all(isapprox(sum(abs2,surface.coords[:,i]),6.25;atol=1e-12) for i in axes(surface.coords,2))
end
let filled=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_sphere.step"))
    @test filled isa Tessella.Mesh && validate(filled).ok && size(filled.tets,2)>0
end
end
