module BRepTorusSegmentSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_torus_segment.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test size(surface.tris,2)>0
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,-1.] &&
          vec(maximum(surface.coords,dims=2))≈[4.,4.,1.]
end
end
