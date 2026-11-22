# examples/boolean/periodic.geo: `Periodic Surface{9} = {8} Rotate{...}` pairs
# the two end disks of a torus segment, each bounded by one closed circle
# with a single point. The declaration stopped with "Surface[9] needs at
# least three distinct point tags"; Gmsh 4.15.2 pairs them ("Setting curve
# master 2 - 3" here) and so does Tessella now. Meshing the pair is a
# separate item (the 1D copy of a closed periodic curve).
module GeoOccPeriodicDiskTest
using Test, Tessella

@testset "Periodic Surface between one-point disks (Gmsh 4.15.2)" begin
    A=Tessella.API
    path=tempname()*".geo"
    write(path,"""
        SetFactory("OpenCASCADE");
        Torus(1) = {3,0,0, 1.5, 1, Pi/3};
        Periodic Surface{3} = {2} Rotate{{0,0,1}, {3,0,0}, Pi/3 };
        """)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(path)
        @test A.mesh.get_periodic(2,[3])==[2]
        @test A.mesh.get_periodic(1,[2])==[3]
    finally
        A.finalize()
        rm(path;force=true)
    end
end
end
