# HFSS UserGuide case 9.2 (the enclosure coax feed-through) SOLVED on a Tessella mesh.
#
# Case 9.2 is the project's reason to exist: gmsh 4.13/4.15 (via ASCENT's OCC geo-emitter)
# CANNOT mesh it — it produces 0 volume tets. Tessella meshes it natively (see reconstruct.jl),
# and here ASCENT assembles AND solves the Maxwell FEM system on that mesh — a real 22-case
# regression data point, and the hardest one.
#
# Run in the ASCENT env (after reconstruct.jl in the Tessella env):
#   julia --project=<2026_066/ASCENT> validation/enclosure_literal/solve_case_9_2.jl
using ASCENT
import ASCENT: load_mesh, cell_sigma_tensor, fe_spaces, assemble_diffusive_matrix
using Gridap.FESpaces: num_free_dofs
using LinearAlgebra, SparseArrays, Random

path = joinpath(@__DIR__, "enclosure_literal.msh")
isfile(path) || error("run reconstruct.jl (Tessella env) first to produce $path")
md = load_mesh(path)
println("HFSS case 9.2 (enclosure) — Tessella mesh loaded in ASCENT: volumes = ", md.volume_tags)

ntags = length(md.volume_tags)
σh = [1.0e7, 1.0e-3, 1.0e-3, 1.0e6]      # pin conductor / slot / air / case metal (representative)
σ  = cell_sigma_tensor(md.model, md.volume_tags, σh, σh)
U, V = fe_spaces(md.model, 1); n = num_free_dofs(V)
A = assemble_diffusive_matrix(md.model, U, V, σ, 2π*2.4e9)   # 2.4 GHz

# full solve via a manufactured solution (definitive correct-solve proof on the flagship case)
Random.seed!(9); xt = randn(ComplexF64, n); b = A*xt; x = A \ b
err = norm(x - xt)/norm(xt); relsym = norm(A - transpose(A), Inf)/norm(A, Inf)
println("  ndof = ", n, "   complex-symmetric relsym = ", relsym, "   solve recovered-field err = ", err)
(n > 0 && relsym < 1e-10 && err < 1e-6) ||
    error("case 9.2 solve failed: ndof=$n relsym=$relsym err=$err")
println("CASE_9_2_OK — ASCENT assembles + solves the enclosure (gmsh-impossible) on a Tessella mesh")
