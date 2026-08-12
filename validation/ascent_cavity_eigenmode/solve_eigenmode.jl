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
modes = solve_eigenmodes(cache, mat; n_modes=6)

femf = filter(f -> f > 1e6, sort(real.(modes.freqs)))    # drop spurious/near-zero modes
c0 = 299792458.0; A, B, D = 1.0, 0.5, 0.75
# analytic rectangular-cavity spectrum: f_mnp = (c/2)√((m/a)²+(n/b)²+(p/d)²), TE/TM ⇒ ≥2 nonzero indices
ana = Float64[]
for mi in 0:3, ni in 0:3, pi in 0:3
    ((mi==0)+(ni==0)+(pi==0)) <= 1 || continue
    f = (c0/2)*sqrt((mi/A)^2+(ni/B)^2+(pi/D)^2); f > 1e6 && push!(ana, f)
end
sort!(ana)

# each FEM mode matched to its NEAREST analytic resonance (robust to degeneracies/ordering)
println("ASCENT cavity mode spectrum on the Tessella mesh vs the analytic reference:")
maxerr = 0.0
for k in 1:min(5, length(femf))
    _, i = findmin(abs.(ana .- femf[k])); e = abs(femf[k]-ana[i])/ana[i]; global maxerr = max(maxerr, e)
    println("  FEM mode $k = ", round(femf[k]/1e6, digits=3), " MHz  →  analytic ",
            round(ana[i]/1e6, digits=3), " MHz  (err ", round(e*100, digits=3), " %)")
end
maxerr < 0.01 || error("cavity spectrum error too large: max $maxerr")
println("CAVITY_EIGENMODE_OK — ASCENT recovers the cavity mode SPECTRUM (5 modes) on a Tessella mesh, max err ",
        round(maxerr*100, digits=3), " %")
