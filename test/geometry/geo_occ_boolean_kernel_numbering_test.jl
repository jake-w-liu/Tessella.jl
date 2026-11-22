module GeoOCCBooleanKernelNumberingChecks
using Test,Tessella
const A=Tessella.API

# `.geo` volume Booleans run through the kernel publication like the API's:
# `.geo` volume Booleans run through the kernel publication like the API's:
# a fragment of two abutting boxes and a box across their common face
# (compsolid2.geo with a box for the sphere) numbers its entities from 1 as
# Gmsh 4.15.2 does: 24/40/22/4, masses 7.5, 0.5, 0.5, 7.5.
@testset ".geo fragment numbering through the kernel" begin
    path=tempname()*".geo"
    write(path,"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 2,2,2};
        Box(2) = {1.5,0.5,0.5, 1,1,1};
        Box(3) = {2,0,0, 2,2,2};
        f() = BooleanFragments { Volume{1:3}; Delete; }{};
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test [[t for (_,t) in A.model.get_entities(d)] for d in 0:3]==
            [collect(1:24),collect(1:40),collect(1:22),collect(1:4)]
        # `getMass` of `.geo` volumes binds them into the kernel lazily
        Tessella.API._with_model(m->Tessella.GeoExec._geo_occ_bind_entities!(m,
            [(3,t) for t in 1:4],Tessella.GeoExec._geo_occ_prepare_kernel!(m),"test"))
        @test [A.model.occ.get_mass(3,t) for t in 1:4]≈[7.5,0.5,0.5,7.5] rtol=1e-12
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
