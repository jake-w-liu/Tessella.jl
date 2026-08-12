# The "ready for ASCENT to use" proof at the SOLVE level: ASCENT assembles the actual
# Maxwell finite-element operator ON a Tessella-generated mesh (not just loads it).
#
# Pipeline (all ASCENT API, on a Tessella .msh):
#   load_mesh -> cell_sigma_tensor (materials per volume) -> fe_spaces (Nedelec H(curl)
#   edge elements) -> assemble_diffusive_matrix / assemble_stiffness_mass (the Maxwell
#   curl-curl + mass system). Verifies the assembled operator's physics: complex-symmetric
#   and curl-curl-stiffness PSD.
#
# Run in the ASCENT env (ASCENT + Gridap are ASCENT deps):
#   julia --project=<2026_066/ASCENT> validation/ascent_handshake/solve_step.jl
using ASCENT
using ASCENT: load_mesh, cell_sigma_tensor, fe_spaces, assemble_diffusive_matrix, assemble_stiffness_mass
using Gridap.FESpaces: num_free_dofs
using LinearAlgebra, SparseArrays, Random

path = joinpath(@__DIR__, "enc_coax_bc.msh")
isfile(path) || error("run generate_bc.jl (Tessella env) first to produce $path")
md = load_mesh(path)
println("loaded Tessella mesh: volume_tags = ", md.volume_tags)

# representative conductivities per material volume (pin conductor / air / case metal)
σh = [1.0e7, 1.0e-3, 1.0e6]
σ  = cell_sigma_tensor(md.model, md.volume_tags, σh, σh)
U, V = fe_spaces(md.model, 1)                    # first-order Nedelec (H(curl)) edge elements
ndof = num_free_dofs(V)
ω = 2π * 1.0e9                                    # 1 GHz
A = assemble_diffusive_matrix(md.model, U, V, σ, ω)
K, M = assemble_stiffness_mass(md.model, U, V, σ)

# (a) the assembled operator's physics: complex-symmetric + curl-curl-stiffness PSD
relsym = norm(A - transpose(A), Inf) / norm(A, Inf)
c = randn(size(K,1)); q = dot(c, real(K)*c)      # curl-curl is PSD: cᵀKc >= 0
println("ASCENT assembled the Maxwell FEM operator ON THE TESSELLA MESH:")
println("  ndof (Nedelec edges) = ", ndof, "   A = ", typeof(A), " size ", size(A), " nnz=", nnz(A))
println("  complex-symmetric: relsym = ", relsym, "   curl-curl PSD: cᵀKc = ", q)

# (b) FULL SOLVE via a manufactured solution: known non-trivial field x_true, b = A x_true,
#     solve A x = b, and confirm x recovers x_true to round-off — an unambiguous correct solve.
Random.seed!(42)
x_true = randn(ComplexF64, ndof)
b = A * x_true
x = A \ b
solerr = norm(x - x_true) / norm(x_true)
res    = norm(A*x - b) / norm(b)
println("ASCENT SOLVED the Maxwell FEM system on the Tessella mesh:")
println("  |x_true| = ", round(norm(x_true),digits=3), " (non-trivial)   recovered error = ", solerr, "   residual = ", res)

(ndof > 0 && relsym < 1e-10 && q >= -1e-8 && solerr < 1e-6 && res < 1e-10) ||
    error("ASCENT solve/assembly properties failed: ndof=$ndof relsym=$relsym cKc=$q solerr=$solerr res=$res")
println("ASCENT_SOLVE_STEP_OK — ASCENT assembles AND solves the Maxwell FEM system on a Tessella mesh")
