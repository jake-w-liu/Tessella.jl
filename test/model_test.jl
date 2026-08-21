using Test
using Tessella
using Tessella.MeshTypes: ntris, ntets, validate, tet_volume, node

@testset "entity model and .geo execution" begin
    m=GeoModel()
    p1=add_point!(m,0,0,0; tag=1, mesh_size=0.5)
    p2=add_point!(m,1,0,0; tag=2, mesh_size=0.5)
    p3=add_point!(m,1,1,0; tag=3, mesh_size=0.5)
    p4=add_point!(m,0,1,0; tag=4, mesh_size=0.5)
    @test (p1,p2,p3,p4)==(1,2,3,4)
    add_line!(m,1,2; tag=1); add_line!(m,2,3; tag=2)
    add_line!(m,3,4; tag=3); add_line!(m,4,1; tag=4)
    add_curve_loop!(m,[1,2,3,4]; tag=1)
    add_plane_surface!(m,[1]; tag=1)
    add_physical_group!(m,2,[1]; tag=10, name="front")
    surf=mesh_model_surface(m,1)
    @test validate(surf).ok
    @test ntris(surf)>0
    area=0.0
    for t in 1:ntris(surf)
        a=node(surf,surf.tris[1,t]); b=node(surf,surf.tris[2,t]); c=node(surf,surf.tris[3,t])
        area+=abs((b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2]))/2
    end
    @test area≈1.0 atol=1e-12
    @test_throws ArgumentError add_line!(m,1,1)
    @test_throws ArgumentError mesh_model_surface(m,99)

    box=GeoModel()
    add_box!(box,0,0,0,2,1,1; tag=1)
    vol=mesh_model_volume(box,1)
    @test validate(vol).ok
    @test ntets(vol)>0
    V=sum(tet_volume(node(vol,vol.tets[1,t]),node(vol,vol.tets[2,t]),
                     node(vol,vol.tets[3,t]),node(vol,vol.tets[4,t])) for t in 1:ntets(vol))
    @test V≈2.0 atol=1e-12

    geo=mktemp() do path,io
        write(io, """
            Point(1) = {0, 0, 0, 0.5};
            Point(2) = {1, 0, 0, 0.5};
            Point(3) = {1, 1, 0, 0.5};
            Point(4) = {0, 1, 0, 0.5};
            Line(1) = {1, 2};
            Line(2) = {2, 3};
            Line(3) = {3, 4};
            Line(4) = {4, 1};
            Line Loop(1) = {1, 2, 3, 4};
            Plane Surface(1) = {1};
            Physical Surface("front", 10) = {1};
            """)
        close(io)
        execute_geo(path; mesh_dim=2)
    end
    @test geo.mesh!==nothing
    @test validate(geo.mesh).ok
    @test ntris(geo.mesh)>0
    @test_throws ArgumentError execute_geo("/no/such/file.geo")
    @test_throws ArgumentError import_step("part.step")
    @test_throws ArgumentError import_iges("part.iges")

    cyl=GeoModel()
    add_cylinder!(cyl,0,0,0,0,0,2,1; tag=1)
    cvol=mesh_model_volume(cyl,1)
    @test validate(cvol).ok
    @test ntets(cvol)>0
    CV=sum(tet_volume(node(cvol,cvol.tets[1,t]),node(cvol,cvol.tets[2,t]),
                      node(cvol,cvol.tets[3,t]),node(cvol,cvol.tets[4,t]))
           for t in 1:ntets(cvol))
    @test CV≈0.5*24*sin(2π/24)*2 atol=1e-12

    sph=GeoModel()
    add_sphere!(sph,0,0,0,1; tag=1)
    svol=mesh_model_volume(sph,1)
    @test validate(svol).ok
    @test ntets(svol)>0

    bool=GeoModel()
    add_box!(bool,0,0,0,2,1,1; tag=1)
    add_box!(bool,0,0,0,1,1,1; tag=2)
    boolean_volumes!(bool,:difference,1,2; tag=3)
    bvol=mesh_model_volume(bool,3)
    @test validate(bvol).ok
    @test ntets(bvol)>0
    BV=sum(tet_volume(node(bvol,bvol.tets[1,t]),node(bvol,bvol.tets[2,t]),
                      node(bvol,bvol.tets[3,t]),node(bvol,bvol.tets[4,t]))
           for t in 1:ntets(bvol))
    @test BV≈1.0 atol=1e-12

    cylgeo=mktemp() do path,io
        write(io, """
            SetFactory("OpenCASCADE");
            Cylinder(1) = {0, 0, 0, 0, 0, 2, 1};
            """)
        close(io)
        execute_geo(path; mesh_dim=3)
    end
    @test cylgeo.mesh!==nothing
    @test validate(cylgeo.mesh).ok
    @test ntets(cylgeo.mesh)>0
    CG=sum(tet_volume(node(cylgeo.mesh,cylgeo.mesh.tets[1,t]),
                      node(cylgeo.mesh,cylgeo.mesh.tets[2,t]),
                      node(cylgeo.mesh,cylgeo.mesh.tets[3,t]),
                      node(cylgeo.mesh,cylgeo.mesh.tets[4,t]))
           for t in 1:ntets(cylgeo.mesh))
    @test CG≈0.5*24*sin(2π/24)*2 atol=1e-12

    boolgeo=mktemp() do path,io
        write(io, """
            SetFactory("OpenCASCADE");
            Box(1) = {
              0, 0, 0, 2, 1, 1
            };
            Box(2) = {0, 0, 0, 1, 1, 1};
            BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };
            """)
        close(io)
        execute_geo(path; mesh_dim=3)
    end
    @test length(boolgeo.model.volumes)==1
    @test validate(boolgeo.mesh).ok
    @test ntets(boolgeo.mesh)>0
    BG=sum(tet_volume(node(boolgeo.mesh,boolgeo.mesh.tets[1,t]),
                      node(boolgeo.mesh,boolgeo.mesh.tets[2,t]),
                      node(boolgeo.mesh,boolgeo.mesh.tets[3,t]),
                      node(boolgeo.mesh,boolgeo.mesh.tets[4,t]))
           for t in 1:ntets(boolgeo.mesh))
    @test BG≈1.0 atol=1e-12

    shifted=mktemp() do path,io
        write(io, """
            Box(1) = {0, 0, 0, 1, 1, 1};
            Translate {2, 0, 0} { Volume{1}; };
            """)
        close(io)
        execute_geo(path; mesh_dim=3)
    end
    @test validate(shifted.mesh).ok
    xs=(shifted.mesh.coords[1,i] for i in 1:size(shifted.mesh.coords,2))
    @test minimum(xs)≈2.0 atol=1e-12
    @test maximum(xs)≈3.0 atol=1e-12
    SV=sum(tet_volume(node(shifted.mesh,shifted.mesh.tets[1,t]),
                      node(shifted.mesh,shifted.mesh.tets[2,t]),
                      node(shifted.mesh,shifted.mesh.tets[3,t]),
                      node(shifted.mesh,shifted.mesh.tets[4,t]))
           for t in 1:ntets(shifted.mesh))
    @test SV≈1.0 atol=1e-12

    two=mktemp() do path,io
        write(io, "Box(1) = {0, 0, 0, 1, 1, 1};\nBox(2) = {2, 0, 0, 1, 1, 1};\n")
        close(io)
        path
    end
    @test_throws ArgumentError execute_geo(two; mesh_dim=3)
    extrude=mktemp() do path,io
        write(io, "Extrude {0, 0, 1} { Surface{1}; };\n")
        close(io)
        path
    end
    @test_throws ArgumentError execute_geo(extrude)
end
