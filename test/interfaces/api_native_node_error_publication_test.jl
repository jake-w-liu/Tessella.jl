module APINativeNodeErrorPublicationTests
using Tessella,Test
const VALIDATION=joinpath(@__DIR__,"..","..","validation","native_geometry")
include(joinpath(VALIDATION,"native_node_error_publication_v2.jl"))
const PREFIX=get(ENV,"TESSELLA_NATIVE_NODE_ERROR_OUTPUT_PREFIX",joinpath(mktempdir(),"native_node_error"))
@testset "Native node rejection return, diagnostic and Msg lifetime" begin
    NativeNodeErrorPublicationV2.run(PREFIX*".errors.toml")
end
end
