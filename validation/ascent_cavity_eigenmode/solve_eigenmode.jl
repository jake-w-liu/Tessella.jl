# ASCENT eigenmode (resonance) solve on the Tessella cavity mesh, validated against the
# closed-form analytic frequency — a real end-to-end physics regression case.
# Run in the ASCENT env (after generate_cavity.jl in the Tessella env):
#   julia --project=<2026_066/ASCENT> validation/ascent_cavity_eigenmode/solve_eigenmode.jl
using ASCENT
import ASCENT: load_mesh, Material, build_assembly_cache, solve_eigenmodes

path = joinpath(@__DIR__, "cavity.msh")
isfile(path) || error("run generate_cavity.jl (Tessella env) first to produce $path")
md = load_mesh(path)
ntags = length(md.volume_tags)
mat = Material(zeros(ntags), zeros(ntags); tag_names=md.volume_tags,
               εr=fill(ComplexF64(1.0),ntags), μr=fill(ComplexF64(1.0),ntags))  # vacuum cavity, PEC walls
cache = build_assembly_cache(md.model, mat, 1)
modes = solve_eigenmodes(cache, mat; n_modes=4)

fs = filter(f -> f > 1e6, sort(real.(modes.freqs)))     # drop spurious/near-zero modes
fnum = fs[1]
c0 = 299792458.0; A, B, D = 1.0, 0.5, 0.75
f101 = (c0/2) * sqrt(1/A^2 + 1/D^2)
err = abs(fnum - f101) / f101
println("ASCENT eigenmode solve on the Tessella cavity mesh:")
println("  lowest resonant f (FEM) = ", round(fnum/1e6, digits=3), " MHz")
println("  analytic TE101          = ", round(f101/1e6, digits=3), " MHz")
println("  relative error          = ", round(err*100, digits=3), " %")
err < 0.05 || error("cavity resonance error too large: $err")
println("CAVITY_EIGENMODE_OK — ASCENT computes the correct cavity resonance (physics) on a Tessella mesh")
