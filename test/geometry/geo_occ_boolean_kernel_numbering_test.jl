module GeoOCCBooleanKernelNumberingChecks
using Test,Tessella
const A=Tessella.API

# `.geo` volume Booleans run through the kernel publication like the API's:
# compsolid2.geo's fragment of two boxes and a sphere numbers its entities
# from 1 as Gmsh 4.15.2 does, and a vertex is stored exactly where its edge
# ends, so the closed circle bounding the sphere cap meshes conformingly.
@testset ".geo fragment numbering and exact vertices" begin
    path=tempname()*".geo"
    write(path,"""
        SetFactory("OpenCASCADE");
        Box(1) = {0,0,0, 2,2,2};
        Sphere(2) = {2.2, 1, 1, 0.5};
        Box(3) = {2,0,0, 2,2,2};
        f() = BooleanFragments { Volume{1:3}; Delete; }{};
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test [[t for (_,t) in A.model.get_entities(d)] for d in 0:3]==
            [collect(1:15),collect(1:24),collect(1:14),collect(1:4)]
        for (_,c) in A.model.get_entities(1)
            a,b=A.model.get_adjacencies(1,c)[2]
            a==b || continue
            lo,_=A.model.get_parametrization_bounds(1,c)
            @test A.model.get_value(1,c,[lo[1]])==A.model.get_value(0,a,Float64[])
        end
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
