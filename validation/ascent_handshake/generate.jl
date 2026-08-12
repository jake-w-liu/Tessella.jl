# Tessella side of the ASCENT mesh handshake: build the ENC-COAX conforming
# partition (pin/air/case) and write it as gmsh MSH v4.1 with the three physical
# volumes — the exact format ASCENT's load_mesh ingests. Run in the Tessella env.
#   julia --project=<Tessella.jl> validation/ascent_handshake/generate.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.Geometry, Tessella.IO

outer  = box_surface(0.,6., 0.,6., 0.,6.)                       # case outer (vol 216)
cavity = box_surface(1.,5., 1.,5., 1.,5.)                       # air+pin cavity (vol 64)
pin    = cylinder_surface((3.,3.,1.5), (0.,0.,1.), 0.7, 3.0; nθ=8, nz=2)  # coax pin
mc = tetrahedralize_conforming([pin, cavity, outer])            # one conforming mesh
diag = validate(mc)
diag.ok || error("generate.jl: conforming mesh invalid — $(diag.messages)")

out = joinpath(@__DIR__, "ascent_coax.msh")
write_msh(out, mc; version=4.1,
          physical_names=Dict((3,Int32(1))=>"coax_pin", (3,Int32(2))=>"air", (3,Int32(3))=>"case"))
println("wrote $out  (tets=$(ntets(mc)), regions=$(tets_per_region(mc)), $(filesize(out)) bytes)")
