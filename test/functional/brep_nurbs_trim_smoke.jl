module BRepNurbsTrimSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_bump_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==85 && Tessella.MeshTypes.ntris(surface)==166
    @test vec(minimum(surface.coords,dims=2))≈[0.,0.,0.] atol=1e-12
    @test vec(maximum(surface.coords,dims=2))≈[2.,1.,1.125] atol=1e-12
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_pipe.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==112 && Tessella.MeshTypes.ntris(surface)==220
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_pillow.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==130 && Tessella.MeshTypes.ntris(surface)==256
end
# Full pillow recovery stays in geometry/brep_nurbs_trim_test.jl (opt-in).
end
