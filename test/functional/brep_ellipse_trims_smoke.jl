module BRepEllipseTrimsSmokeChecks
using Test, Tessella
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_elliptical_pipe.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==96 && Tessella.MeshTypes.ntris(surface)==188
    @test vec(minimum(surface.coords,dims=2))≈[-2.,-1.,0.] &&
          vec(maximum(surface.coords,dims=2))≈[2.,1.,1.]
end
let filled=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_elliptical_pipe.step"))
    @test filled isa Tessella.Mesh && validate(filled).ok && size(filled.tets,2)>0
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_ellipse_arc_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==54 && Tessella.MeshTypes.ntris(surface)==104
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_ellipse_hole_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==104 && Tessella.MeshTypes.ntris(surface)==208
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_spline_partial_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok && size(surface.tets,2)==0
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_prism_multi_rim.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==38 && Tessella.MeshTypes.ntris(surface)==72
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_rational_bulge_plate.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok
    @test Tessella.MeshTypes.nnodes(surface)==38 && Tessella.MeshTypes.ntris(surface)==72
end
let surface=import_step(joinpath(@__DIR__,"..","fixtures","brep_occ_nurbs_rational_bump.step");fill=false)
    @test surface isa Tessella.Mesh && validate(surface).ok && size(surface.tets,2)==0
end
end
