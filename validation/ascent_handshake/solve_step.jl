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
using LinearAlgebra, SparseArrays

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

relsym = norm(A - transpose(A), Inf) / norm(A, Inf)
c = randn(size(K,1)); q = dot(c, real(K)*c)      # curl-curl is PSD: c'Kc >= 0
println("ASCENT assembled the Maxwell FEM operator ON THE TESSELLA MESH:")
println("  ndof (Nedelec edges) = ", ndof)
println("  A = ", typeof(A), " size ", size(A), " nnz=", nnz(A))
println("  complex-symmetric: relsym = ", relsym)
println("  curl-curl stiffness PSD: cᵀKc = ", q)
(ndof > 0 && relsym < 1e-10 && q >= -1e-8) ||
    error("assembly properties failed: ndof=$ndof relsym=$relsym cKc=$q")
println("ASCENT_SOLVE_STEP_OK — ASCENT builds the finite-element Maxwell system on a Tessella mesh")
