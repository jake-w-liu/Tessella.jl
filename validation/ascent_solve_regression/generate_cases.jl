# Generate a suite of distinct Tessella meshes for the ASCENT solve regression — the
# transferable core of a multi-case HFSS regression: prove ASCENT assembles+solves the
# Maxwell FEM system on DIVERSE Tessella geometries (single/multi-region, box/cylinder/
# nested), not just one. Run in the Tessella env:
#   julia --project=<Tessella.jl> validation/ascent_solve_regression/generate_cases.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.Geometry, Tessella.IO

function emit(name, surfaces, names)
    m = tetrahedralize_conforming_exact(surfaces)
    validate(m).ok || error("$name: invalid mesh")
    pn = Dict((3,Int32(i)) => names[i] for i in 1:length(names))
    p = joinpath(@__DIR__, "$name.msh")
    write_msh(p, m; version=4.1, physical_names=pn)
    println("  $name: tets=$(ntets(m)) regions=$(length(tets_per_region(m)))  →  $(basename(p))")
end

emit("box_cavity", [box_surface(0.,1.,0.,1.,0.,1.)], ["dielectric"])
emit("nested_box", [box_surface(0.25,0.75,0.25,0.75,0.25,0.75),
                    box_shell_surface(0.,1.,0.,1.,0.,1., 0.25,0.75,0.25,0.75,0.25,0.75)], ["core","shell"])
emit("coax3",      [cylinder_surface((0.5,0.5,0.3),(0.,0.,1.),0.12,0.6; nθ=8,nz=2),
                    box_surface(0.2,0.8,0.2,0.8,0.2,0.8),
                    box_shell_surface(0.,1.,0.,1.,0.,1., 0.2,0.8,0.2,0.8,0.2,0.8)], ["pin","air","case"])
emit("cyl_cavity", [cylinder_surface((0.,0.,0.),(0.,0.,1.),0.5,1.0; nθ=12)], ["dielectric"])
println("generated 4 Tessella cases for the ASCENT solve regression")
