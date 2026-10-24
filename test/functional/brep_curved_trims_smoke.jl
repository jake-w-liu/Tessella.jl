module BRepCurvedTrimsSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_drilled_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==56 && Tessella.MeshTypes.ntris(surface)==112
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[4.,4.,1.]
end
let filled=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_drilled_plate.step"))
    @test filled isa Tessella.Mesh && validate(filled).ok && size(filled.tets,2)>0
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_hollow_cylinder.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==96 && Tessella.MeshTypes.ntris(surface)==192
end
end
