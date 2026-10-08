# Invalid GEO input can leave parser FILE handles open in the pinned oracle,
# even after gmsh.clear/parser.clear/finalize. A child process owns those
# handles so its exit releases them before the caller removes fixture files.
length(ARGS)>=2 || error("expected Gmsh API path and invalid GEO fixtures")
include(ARGS[1])
gmsh.GMSH_API_VERSION=="4.15.2" || error("expected Gmsh API 4.15.2")
gmsh.initialize(String[],false,false)
try
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.setNumber("General.NumThreads",1)
    gmsh.option.setNumber("Mesh.ElementOrder",1)
    gmsh.option.setNumber("Mesh.FlexibleTransfinite",0)
    version=gmsh.option.getString("General.Version")
    (version=="4.15.2" || startswith(version,"4.15.2-")) || error(
        "expected Gmsh runtime 4.15.2, got $version")
    for path in ARGS[2:end]
        gmsh.clear()
        gmsh.parser.clear()
        rejected=false
        try
            gmsh.merge(path)
        catch err
            err isa InterruptException && rethrow()
            rejected=true
        end
        rejected || error("Gmsh accepted invalid source $(basename(path))")
    end
    println("GEO_CONSTRAINTS_ERROR_PRIMARY_OK cases=$(length(ARGS)-1)")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
