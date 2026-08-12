# Generate a rectangular PEC cavity mesh in Tessella for an ASCENT eigenmode (resonance)
# solve validated against the closed-form analytic frequency — a real physics regression
# case (the kind HFSS UserGuide cavity examples validate against), end-to-end on a
# Tessella mesh. Run in the Tessella env:
#   julia --project=<Tessella.jl> validation/ascent_cavity_eigenmode/generate_cavity.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella, Tessella.Mesh3D, Tessella.MeshTypes, Tessella.IO

const A, B, D = 1.0, 0.5, 0.75            # cavity dimensions (m)
m0 = mesh_box(0.0,A, 0.0,B, 0.0,D; hmax=0.16)     # structured Kuhn mesh, fine enough for 1st-order Nedelec
m  = Mesh(m0.coords; tets=m0.tets, tet_tag=fill(Int32(1), ntets(m0)))   # one physical volume group
validate(m).ok || error("cavity mesh invalid")

out = joinpath(@__DIR__, "cavity.msh")
write_msh(out, m; version=4.1, physical_names=Dict((3,Int32(1))=>"vacuum"))

c0 = 299792458.0
f101 = (c0/2) * sqrt(1/A^2 + 1/D^2)       # TE101 analytic resonance
println("cavity a×b×d = $A×$B×$D  tets=$(ntets(m))  →  $(basename(out))")
println("analytic TE101 = ", round(f101/1e6, digits=3), " MHz")
