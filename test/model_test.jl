using Test
using Tessella
using Tessella.MeshTypes: ntris, ntets, nnodes, validate, tet_volume, node

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

    hole=GeoModel()
    add_point!(hole,0,0,0; tag=1, mesh_size=0.5)
    add_point!(hole,1,0,0; tag=2, mesh_size=0.5)
    add_point!(hole,1,1,0; tag=3, mesh_size=0.5)
    add_point!(hole,0,1,0; tag=4, mesh_size=0.5)
    add_point!(hole,0.25,0.25,0; tag=5, mesh_size=0.5)
    add_point!(hole,0.75,0.25,0; tag=6, mesh_size=0.5)
    add_point!(hole,0.75,0.75,0; tag=7, mesh_size=0.5)
    add_point!(hole,0.25,0.75,0; tag=8, mesh_size=0.5)
    add_line!(hole,1,2; tag=1); add_line!(hole,2,3; tag=2)
    add_line!(hole,3,4; tag=3); add_line!(hole,4,1; tag=4)
    add_line!(hole,5,6; tag=5); add_line!(hole,6,7; tag=6)
    add_line!(hole,7,8; tag=7); add_line!(hole,8,5; tag=8)
    add_curve_loop!(hole,[1,2,3,4]; tag=1)
    add_curve_loop!(hole,[5,6,7,8]; tag=2)
    add_plane_surface!(hole,[1,2]; tag=1)
    hmesh=mesh_model_surface(hole,1)
    @test validate(hmesh).ok
    harea=sum(abs((node(hmesh,hmesh.tris[2,t])[1]-node(hmesh,hmesh.tris[1,t])[1])*
                  (node(hmesh,hmesh.tris[3,t])[2]-node(hmesh,hmesh.tris[1,t])[2])-
                  (node(hmesh,hmesh.tris[3,t])[1]-node(hmesh,hmesh.tris[1,t])[1])*
                  (node(hmesh,hmesh.tris[2,t])[2]-node(hmesh,hmesh.tris[1,t])[2]))/2
              for t in 1:ntris(hmesh))
    @test harea≈0.75 atol=1e-12

    emb=GeoModel()
    add_point!(emb,0,0,0; tag=1, mesh_size=0.5)
    add_point!(emb,1,0,0; tag=2, mesh_size=0.5)
    add_point!(emb,1,1,0; tag=3, mesh_size=0.5)
    add_point!(emb,0,1,0; tag=4, mesh_size=0.5)
    add_point!(emb,0.5,0.5,0; tag=5, mesh_size=0.5)
    add_line!(emb,1,2; tag=1); add_line!(emb,2,3; tag=2)
    add_line!(emb,3,4; tag=3); add_line!(emb,4,1; tag=4)
    add_curve_loop!(emb,[1,2,3,4]; tag=1)
    add_plane_surface!(emb,[1]; tag=1)
    embed!(emb,0,[5],2,1)
    emesh=mesh_model_surface(emb,1)
    @test validate(emesh).ok
    @test any(i->hypot(emesh.coords[1,i]-0.5,emesh.coords[2,i]-0.5,emesh.coords[3,i])<=1e-12,
              1:nnodes(emesh))
    earea=sum(abs((node(emesh,emesh.tris[2,t])[1]-node(emesh,emesh.tris[1,t])[1])*
                  (node(emesh,emesh.tris[3,t])[2]-node(emesh,emesh.tris[1,t])[2])-
                  (node(emesh,emesh.tris[3,t])[1]-node(emesh,emesh.tris[1,t])[1])*
                  (node(emesh,emesh.tris[2,t])[2]-node(emesh,emesh.tris[1,t])[2]))/2
              for t in 1:ntris(emesh))
    @test earea≈1.0 atol=1e-12
    outside=GeoModel()
    add_point!(outside,0,0,0; tag=1, mesh_size=0.5)
    add_point!(outside,1,0,0; tag=2, mesh_size=0.5)
    add_point!(outside,1,1,0; tag=3, mesh_size=0.5)
    add_point!(outside,0,1,0; tag=4, mesh_size=0.5)
    add_point!(outside,2,2,0; tag=5, mesh_size=0.5)
    add_line!(outside,1,2; tag=1); add_line!(outside,2,3; tag=2)
    add_line!(outside,3,4; tag=3); add_line!(outside,4,1; tag=4)
    add_curve_loop!(outside,[1,2,3,4]; tag=1)
    add_plane_surface!(outside,[1]; tag=1)
    embed!(outside,0,[5],2,1)
    @test_throws Exception mesh_model_surface(outside,1)

    cone=GeoModel()
    add_cone!(cone,0,0,0,0,0,2,1,0.5; tag=1)
    conevol=mesh_model_volume(cone,1)
    @test validate(conevol).ok
    @test ntets(conevol)>0
    nθ=24
    expected=0.5*nθ*sin(2π/nθ)*2*(1+0.5+0.25)/3
    CEV=sum(tet_volume(node(conevol,conevol.tets[1,t]),node(conevol,conevol.tets[2,t]),
                       node(conevol,conevol.tets[3,t]),node(conevol,conevol.tets[4,t]))
            for t in 1:ntets(conevol))
    @test CEV≈expected rtol=1e-12

    dilated=GeoModel()
    add_box!(dilated,0,0,0,1,1,1; tag=1)
    dilate_volume!(dilated,1,(0,0,0),2)
    dvol=mesh_model_volume(dilated,1)
    DV=sum(tet_volume(node(dvol,dvol.tets[1,t]),node(dvol,dvol.tets[2,t]),
                      node(dvol,dvol.tets[3,t]),node(dvol,dvol.tets[4,t]))
           for t in 1:ntets(dvol))
    @test DV≈8.0 atol=1e-12

    rotated=GeoModel()
    add_box!(rotated,0,0,0,1,1,1; tag=1)
    rotate_volume!(rotated,1,(0,0,1),(0,0,0),π/2)
    rvol=mesh_model_volume(rotated,1)
    RV=sum(tet_volume(node(rvol,rvol.tets[1,t]),node(rvol,rvol.tets[2,t]),
                      node(rvol,rvol.tets[3,t]),node(rvol,rvol.tets[4,t]))
           for t in 1:ntets(rvol))
    @test RV≈1.0 atol=1e-12
    rx=(rvol.coords[1,i] for i in 1:nnodes(rvol))
    ry=(rvol.coords[2,i] for i in 1:nnodes(rvol))
    @test minimum(rx)≈-1.0 atol=1e-12
    @test maximum(rx)≈0.0 atol=1e-12
    @test minimum(ry)≈0.0 atol=1e-12
    @test maximum(ry)≈1.0 atol=1e-12

    xform=mktemp() do path,io
        write(io, """
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 1, 1, 1};
            Dilate {{0, 0, 0}, 2} { Volume{1}; };
            Rotate {{0, 0, 1}, {0, 0, 0}, $(π/2)} { Volume{1}; };
            """)
        close(io)
        execute_geo(path; mesh_dim=3)
    end
    @test validate(xform.mesh).ok
    XV=sum(tet_volume(node(xform.mesh,xform.mesh.tets[1,t]),
                      node(xform.mesh,xform.mesh.tets[2,t]),
                      node(xform.mesh,xform.mesh.tets[3,t]),
                      node(xform.mesh,xform.mesh.tets[4,t]))
           for t in 1:ntets(xform.mesh))
    @test XV≈8.0 atol=1e-12

    embedgeo=mktemp() do path,io
        write(io, """
            Point(1) = {0, 0, 0, 0.5};
            Point(2) = {1, 0, 0, 0.5};
            Point(3) = {1, 1, 0, 0.5};
            Point(4) = {0, 1, 0, 0.5};
            Point(5) = {0.5, 0.5, 0, 0.5};
            Line(1) = {1, 2};
            Line(2) = {2, 3};
            Line(3) = {3, 4};
            Line(4) = {4, 1};
            Line Loop(1) = {1, 2, 3, 4};
            Plane Surface(1) = {1};
            Point{5} In Surface{1};
            """)
        close(io)
        execute_geo(path; mesh_dim=2)
    end
    @test any(i->hypot(embedgeo.mesh.coords[1,i]-0.5,embedgeo.mesh.coords[2,i]-0.5)<=1e-12,
              1:nnodes(embedgeo.mesh))

    lineemb=GeoModel()
    add_point!(lineemb,0,0,0; tag=1, mesh_size=0.5)
    add_point!(lineemb,1,0,0; tag=2, mesh_size=0.5)
    add_point!(lineemb,1,1,0; tag=3, mesh_size=0.5)
    add_point!(lineemb,0,1,0; tag=4, mesh_size=0.5)
    add_point!(lineemb,0.25,0.5,0; tag=5, mesh_size=0.5)
    add_point!(lineemb,0.75,0.5,0; tag=6, mesh_size=0.5)
    add_line!(lineemb,1,2; tag=1); add_line!(lineemb,2,3; tag=2)
    add_line!(lineemb,3,4; tag=3); add_line!(lineemb,4,1; tag=4)
    add_line!(lineemb,5,6; tag=5)
    add_curve_loop!(lineemb,[1,2,3,4]; tag=1)
    add_plane_surface!(lineemb,[1]; tag=1)
    embed!(lineemb,1,[5],2,1)
    lmesh=mesh_model_surface(lineemb,1)
    @test validate(lmesh).ok
    larea=sum(abs((node(lmesh,lmesh.tris[2,t])[1]-node(lmesh,lmesh.tris[1,t])[1])*
                  (node(lmesh,lmesh.tris[3,t])[2]-node(lmesh,lmesh.tris[1,t])[2])-
                  (node(lmesh,lmesh.tris[3,t])[1]-node(lmesh,lmesh.tris[1,t])[1])*
                  (node(lmesh,lmesh.tris[2,t])[2]-node(lmesh,lmesh.tris[1,t])[2]))/2
              for t in 1:ntris(lmesh))
    @test larea≈1.0 atol=1e-12
    @test Tessella.Model._mesh_covers_segment(lmesh,(0.25,0.5,0.0),(0.75,0.5,0.0))

    volpt=GeoModel()
    add_box!(volpt,0,0,0,1,1,1; tag=1)
    add_point!(volpt,0.2,0.3,0.4; tag=10)
    embed!(volpt,0,[10],3,1)
    vmesh=mesh_model_volume(volpt,1)
    @test validate(vmesh).ok
    @test ntets(vmesh)>0
    VV=sum(tet_volume(node(vmesh,vmesh.tets[1,t]),node(vmesh,vmesh.tets[2,t]),
                      node(vmesh,vmesh.tets[3,t]),node(vmesh,vmesh.tets[4,t]))
           for t in 1:ntets(vmesh))
    @test VV≈1.0 atol=1e-12
    @test any(i->hypot(vmesh.coords[1,i]-0.2,vmesh.coords[2,i]-0.3,vmesh.coords[3,i]-0.4)<=1e-12,
              1:nnodes(vmesh))

    linegeo=mktemp() do path,io
        write(io, """
            Point(1) = {0, 0, 0, 0.5};
            Point(2) = {1, 0, 0, 0.5};
            Point(3) = {1, 1, 0, 0.5};
            Point(4) = {0, 1, 0, 0.5};
            Point(5) = {0.25, 0.5, 0, 0.5};
            Point(6) = {0.75, 0.5, 0, 0.5};
            Line(1) = {1, 2};
            Line(2) = {2, 3};
            Line(3) = {3, 4};
            Line(4) = {4, 1};
            Line(5) = {5, 6};
            Line Loop(1) = {1, 2, 3, 4};
            Plane Surface(1) = {1};
            Line{5} In Surface{1};
            """)
        close(io)
        execute_geo(path; mesh_dim=2)
    end
    @test Tessella.Model._mesh_covers_segment(linegeo.mesh,(0.25,0.5,0.0),(0.75,0.5,0.0))
    lgarea=sum(abs((node(linegeo.mesh,linegeo.mesh.tris[2,t])[1]-node(linegeo.mesh,linegeo.mesh.tris[1,t])[1])*
                   (node(linegeo.mesh,linegeo.mesh.tris[3,t])[2]-node(linegeo.mesh,linegeo.mesh.tris[1,t])[2])-
                   (node(linegeo.mesh,linegeo.mesh.tris[3,t])[1]-node(linegeo.mesh,linegeo.mesh.tris[1,t])[1])*
                   (node(linegeo.mesh,linegeo.mesh.tris[2,t])[2]-node(linegeo.mesh,linegeo.mesh.tris[1,t])[2]))/2
               for t in 1:ntris(linegeo.mesh))
    @test lgarea≈1.0 atol=1e-12

    volgeo=mktemp() do path,io
        write(io, """
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 1, 1, 1};
            Point(10) = {0.2, 0.3, 0.4, 0.5};
            Point{10} In Volume{1};
            """)
        close(io)
        execute_geo(path; mesh_dim=3)
    end
    @test any(i->hypot(volgeo.mesh.coords[1,i]-0.2,volgeo.mesh.coords[2,i]-0.3,
                       volgeo.mesh.coords[3,i]-0.4)<=1e-12, 1:nnodes(volgeo.mesh))
    VG=sum(tet_volume(node(volgeo.mesh,volgeo.mesh.tets[1,t]),
                      node(volgeo.mesh,volgeo.mesh.tets[2,t]),
                      node(volgeo.mesh,volgeo.mesh.tets[3,t]),
                      node(volgeo.mesh,volgeo.mesh.tets[4,t]))
           for t in 1:ntets(volgeo.mesh))
    @test VG≈1.0 atol=1e-12
end
