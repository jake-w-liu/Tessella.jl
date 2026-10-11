module APINativeMsgResetCanonicalPublicOptimizerTests
using Tessella
# Native Msg errors survive initialize/finalize. GenerateMesh resets them;
# perform that public operation before fixtures captured with clean Msg state.
# Exact native evidence is sealed in native_optimizer_dirty_msg_order_witness_v1.
Tessella.API.initialize()
try
    Tessella.API.mesh.generate(0)
finally
    Tessella.API.finalize()
end
include(joinpath(@__DIR__,"api_native_public_optimizer_test.jl"))
include(joinpath(@__DIR__,"api_native_public_optimizer_quad_neighbors_test.jl"))
end
