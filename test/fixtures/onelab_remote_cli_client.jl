# The actual remote CLI is the child under test; its solver child uses only
# Julia stdlibs in onelab_native_client.jl.
using Tessella
Tessella.CLI.main(ARGS)
