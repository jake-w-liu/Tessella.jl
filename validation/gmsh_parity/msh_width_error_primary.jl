# The pinned oracle retains the FILE handle for a rejected binary-width
# header. Process exit owns cleanup; no Tessella code participates here.
length(ARGS)==2 || error("expected pinned Gmsh API and narrow MSH fixture")
include(ARGS[1])
isfile(ARGS[2]) || error("narrow MSH fixture does not exist")
gmsh.GMSH_API_VERSION=="4.15.2" || error("expected Gmsh API 4.15.2")
gmsh.initialize(["gmsh","-v","0"])
try
    version=gmsh.option.getString("General.Version")
    (version=="4.15.2" || startswith(version,"4.15.2-")) ||
        error("expected Gmsh runtime 4.15.2, got $version")
    rejected=false
    try
        gmsh.clear()
        gmsh.open(ARGS[2])
    catch err
        err isa ErrorException || rethrow()
        expected="Binary file has sizeof(size_t) = 4, not matching machine sizeof(size_t) = 8"
        err.msg==expected || error("unexpected narrow-file oracle error: $(err.msg)")
        rejected=true
    end
    rejected || error("Gmsh 4.15.2 accepted a 4-byte size_t file")
    println("MSH_WIDTH_ERROR_PRIMARY_OK size_t_bytes=4 expected_machine_bytes=8")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
