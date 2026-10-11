module BRepConicalSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_frustum.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test vec(minimum(surface.coords,dims=2))≈[-2.,-2.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[2.,2.,3.]
end
let filled=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_cone.step"))
    @test filled isa Tessella.Mesh && validate(filled).ok && size(filled.tets,2)>0
end
end
