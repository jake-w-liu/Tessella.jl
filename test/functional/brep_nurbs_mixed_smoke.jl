module BRepNurbsMixedSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_cyl_dome.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==97 && Tessella.MeshTypes.ntris(surface)==190
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_sph_cap.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==440 && Tessella.MeshTypes.ntris(surface)==876
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_sph_seam.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==760 && Tessella.MeshTypes.ntris(surface)==1516
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_nurbs_arb.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==158 && Tessella.MeshTypes.ntris(surface)==312
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_pcurve_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==85 && Tessella.MeshTypes.ntris(surface)==166
end
# Full curved-cap recovery stays in geometry/brep_nurbs_mixed_test.jl (opt-in).
end
