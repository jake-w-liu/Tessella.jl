# ASCENT solve regression across the Tessella case suite: for each mesh, assemble the
# Maxwell FEM operator and solve it (manufactured solution → machine-precision recovery),
# checking complex-symmetry too. Proves ASCENT solve-usability is ROBUST across diverse
# Tessella geometries — the transferable core of the 22-case HFSS regression.
#
# Run in the ASCENT env (after generate_cases.jl in the Tessella env):
#   julia --project=<2026_066/ASCENT> validation/ascent_solve_regression/solve_all.jl
using ASCENT
using ASCENT: load_mesh, cell_sigma_tensor, fe_spaces, assemble_diffusive_matrix
using Gridap.FESpaces: num_free_dofs
using LinearAlgebra, SparseArrays, Random

cases = ["box_cavity", "nested_box", "coax3", "cyl_cavity"]
allok = true
for name in cases
    path = joinpath(@__DIR__, "$name.msh")
    isfile(path) || error("missing $path — run generate_cases.jl (Tessella env) first")
    md = load_mesh(path)
    nreg = length(md.volume_tags)
    σh = fill(1.0e3, nreg); σh[1] = 1.0e7             # first region conductive
    σ = cell_sigma_tensor(md.model, md.volume_tags, σh, σh)
    U, V = fe_spaces(md.model, 1); n = num_free_dofs(V)
    A = assemble_diffusive_matrix(md.model, U, V, σ, 2π*1.0e9)
    relsym = norm(A - transpose(A), Inf) / norm(A, Inf)
    Random.seed!(7); xt = randn(ComplexF64, n); b = A*xt; x = A \ b
    err = norm(x - xt) / norm(xt)
    ok = (n > 0 && relsym < 1e-10 && err < 1e-6)
    global allok &= ok
    println("  $(rpad(name,12)) regions=$nreg ndof=$(lpad(n,4)) relsym=$(round(relsym,sigdigits=2)) solve_err=$(round(err,sigdigits=2)) $(ok ? "OK" : "FAIL")")
end
allok || error("SOLVE REGRESSION FAILED")
println("SOLVE_REGRESSION_OK — ASCENT assembles + solves Maxwell on every Tessella case")
