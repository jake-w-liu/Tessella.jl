using Test
using Tessella
using Tessella.MeshTypes: Mesh, ntris, ntets, nnodes, validate, tet_volume, node,
                          mesh_crc

@testset "entity model validation and tag semantics" begin
    oriented=GeoModel()
    for (tag,(x,y)) in enumerate(((0,0),(1,0),(1,1),(0,1)))
        add_point!(oriented,x,y,0; tag=tag, mesh_size=0.5)
    end
    add_line!(oriented,1,2; tag=1)
    add_line!(oriented,2,3; tag=2)
    add_line!(oriented,3,4; tag=3)
    add_line!(oriented,4,1; tag=4)
    @test add_curve_loop!(oriented,[-4,-3,-2,-1]; tag=7)==7
    @test Tessella.Model._loop_points(oriented,7)==[1,4,3,2]
    @test add_plane_surface!(oriented,[7]; tag=7)==7
    reversed_mesh=mesh_model_surface(oriented,7)
    reversed_area=sum(abs((node(reversed_mesh,reversed_mesh.tris[2,t])[1]-
                           node(reversed_mesh,reversed_mesh.tris[1,t])[1])*
                          (node(reversed_mesh,reversed_mesh.tris[3,t])[2]-
                           node(reversed_mesh,reversed_mesh.tris[1,t])[2])-
                          (node(reversed_mesh,reversed_mesh.tris[3,t])[1]-
                           node(reversed_mesh,reversed_mesh.tris[1,t])[1])*
                          (node(reversed_mesh,reversed_mesh.tris[2,t])[2]-
                           node(reversed_mesh,reversed_mesh.tris[1,t])[2]))/2
                      for t in 1:ntris(reversed_mesh))
    @test validate(reversed_mesh).ok
    @test reversed_area≈1.0 atol=1e-12
    @test mesh_crc(reversed_mesh).sha==
          "d88d6244f73026b3450f9a4be9a0402160cf86ef41ddc047586c0f889b0955b0"

    loops_before=copy(oriented.loops)
    @test_throws ArgumentError add_curve_loop!(oriented,[1,3,2,4]; tag=8)
    @test_throws ArgumentError add_curve_loop!(oriented,[1,2,3]; tag=8)
    @test_throws ArgumentError add_curve_loop!(oriented,[0,2,3,4]; tag=8)
    @test_throws ArgumentError add_curve_loop!(oriented,(true,2,3,4); tag=8)
    @test_throws ArgumentError add_curve_loop!(oriented,[-(big(2)^100),2,3,4]; tag=8)
    @test oriented.loops==loops_before

    exhausted=GeoModel()
    @test add_point!(exhausted,0,0,0; tag=typemax(Int32))==typemax(Int32)
    points_before=copy(exhausted.points)
    @test_throws ArgumentError add_point!(exhausted,1,0,0)
    @test_throws ArgumentError add_point!(GeoModel(),0,0,0; tag=big(2)^100)
    @test exhausted.points==points_before

    invalid_primitives=GeoModel()
    @test_throws ArgumentError add_box!(invalid_primitives,0,0,0,Inf,1,1)
    @test_throws ArgumentError add_cylinder!(invalid_primitives,0,0,0,0,0,1,Inf)
    @test_throws ArgumentError add_cylinder!(invalid_primitives,0,0,0,
                                             floatmax(Float64),floatmax(Float64),0,1)
    @test_throws ArgumentError add_sphere!(invalid_primitives,0,0,0,NaN)
    @test_throws ArgumentError add_cone!(invalid_primitives,0,0,0,0,0,1,1,Inf)
    @test isempty(invalid_primitives.volumes)

    transformed=GeoModel()
    add_box!(transformed,0,0,0,1,2,3; tag=1)
    box_before=transformed.box_extents[1]
    @test_throws ArgumentError translate_volume!(transformed,1,(0,0))
    @test_throws ArgumentError translate_volume!(transformed,1,(Inf,0,0))
    @test_throws ArgumentError dilate_volume!(transformed,1,(0,0),2)
    @test_throws ArgumentError dilate_volume!(transformed,1,(0,0,0),Inf)
    @test_throws ArgumentError rotate_volume!(transformed,1,(0,0),(0,0,0),π/2)
    @test_throws ArgumentError rotate_volume!(transformed,1,(0,0,1),(0,0),π/2)
    @test_throws ArgumentError rotate_volume!(transformed,1,(0,0,1),(0,0,0),1e300)
    @test transformed.box_extents[1]==box_before
    @test translate_volume!(transformed,1,(2,-1,0.5))==1
    @test transformed.box_extents[1]==(2.0,-1.0,0.5,1.0,2.0,3.0)
    oversized=GeoModel()
    add_sphere!(oversized,floatmax(Float64),0,0,floatmax(Float64); tag=1)
    sphere_before=oversized.spheres[1]
    @test_throws ArgumentError translate_volume!(oversized,1,(floatmax(Float64),0,0))
    @test oversized.spheres[1]==sphere_before
    @test_throws ArgumentError dilate_volume!(oversized,1,(0,0,0),2)
    @test oversized.spheres[1]==sphere_before

    translated_primitives=GeoModel()
    add_cylinder!(translated_primitives,0,0,0,0,0,2,1;tag=1)
    add_sphere!(translated_primitives,0,0,0,1;tag=2)
    add_cone!(translated_primitives,0,0,0,0,0,2,1,0.5;tag=3)
    for tag in 1:3
        @test translate_volume!(translated_primitives,tag,(1,2,3))==tag
    end
    @test translated_primitives.cylinders[1].center==(1.0,2.0,3.0)
    @test translated_primitives.spheres[2].center==(1.0,2.0,3.0)
    @test translated_primitives.cones[3].center==(1.0,2.0,3.0)
    boolean_volumes!(translated_primitives,:union,1,2;tag=4)
    @test_throws ArgumentError translate_volume!(translated_primitives,4,(1,0,0))
    @test translated_primitives.booleans[4]==(op=:union,a=1,b=2)

    embedded=GeoModel()
    add_box!(embedded,0,0,0,1,1,1; tag=1)
    add_point!(embedded,0.2,0.2,0.2; tag=1)
    add_point!(embedded,0.3,0.3,0.3; tag=2)
    @test_throws ArgumentError embed!(embedded,true,[1],3,1)
    @test_throws ArgumentError embed!(embedded,0,[1],true,1)
    @test_throws ArgumentError embed!(embedded,0,[1,99],3,1)
    @test !haskey(embedded.embeds,(3,1))
    @test embed!(embedded,0,[1],3,1)==1
    embedded_before=copy(embedded.embeds[(3,1)])
    @test_throws ArgumentError embed!(embedded,0,[2,2],3,1)
    @test_throws ArgumentError embed!(embedded,0,[1],3,1)
    @test embedded.embeds[(3,1)]==embedded_before

    physical=GeoModel()
    add_box!(physical,0,0,0,1,1,1; tag=123)
    @test add_physical_group!(physical,3,[123])==1
    @test Tessella.Model.model_physical_tags(physical,3,1)==[123]
    returned=Tessella.Model.model_physical_tags(physical,3,1)
    push!(returned,999)
    @test Tessella.Model.model_physical_tags(physical,3,1)==[123]
    groups_before=copy(physical.physical)
    @test_throws ArgumentError add_physical_group!(physical,3,[999])
    @test_throws ArgumentError add_physical_group!(physical,3,[123,123])
    @test_throws ArgumentError add_physical_group!(physical,true,[123])
    @test physical.physical==groups_before
    @test_throws ArgumentError Tessella.Model.model_entity(physical,true,123)
    @test_throws ArgumentError Tessella.Model.model_physical_tags(physical,3,big(2)^100)
    @test_throws ArgumentError Tessella.Model.set_physical_name!(physical,true,1,"bad")
end

@testset "persistent native straight-curve periodic constraints" begin
    function periodic_square()
        model=GeoModel()
        for (tag,(x,y)) in enumerate(((0.0,0.0),(1.0,0.0),
                                      (1.0,1.0),(0.0,1.0)))
            add_point!(model,x,y,0;tag=tag,mesh_size=0.25)
        end
        for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            add_line!(model,first,last;tag=tag)
        end
        add_curve_loop!(model,[1,2,3,4];tag=1)
        add_plane_surface!(model,[1];tag=1)
        return model
    end

    function periodic_curve_graph()
        model=GeoModel()
        coordinates=(
            (1,0.0,0.0),(2,1.0,0.0),(3,1.0,1.0),(4,0.0,1.0),
            (101,0.2,0.2),(102,0.8,0.2),
            (103,0.2,0.5),(104,0.8,0.5),
            (105,0.2,0.8),(106,0.8,0.8),(107,0.425,0.8))
        for (tag,x,y) in coordinates
            add_point!(model,x,y,0;tag=tag,mesh_size=0.4)
        end
        endpoints=((1,1,2),(2,2,3),(3,3,4),(4,4,1),
                   (30,101,102),(20,103,104),(10,106,105))
        for (tag,first,last) in endpoints
            add_line!(model,first,last;tag=tag)
        end
        add_curve_loop!(model,[1,2,3,4];tag=1)
        add_plane_surface!(model,[1];tag=1)
        embed!(model,0,[107],2,1)
        embed!(model,1,[30,20,10],2,1)
        return model
    end

    translation=Float64[
        1,0,0,1,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1,
    ]
    model=periodic_square()
    @test set_periodic!(model,1,[2],[4],translation)===nothing
    translation[4]=99
    constraints=model_periodic_constraints(model)
    @test length(constraints)==1
    constraint=only(constraints)
    @test constraint.slave_entity==2
    @test constraint.master_entity==4
    @test constraint.reversed
    @test constraint.affine[4]==1
    empty!(constraints)
    @test length(model_periodic_constraints(model))==1
    invalid_mesh=Mesh([0.0 1.0 2.0;0.0 0.0 0.0;0.0 0.0 0.0];
                      tris=reshape(Int32[1,2,3],3,1))
    @test_throws ArgumentError model_periodic_nodes(model,invalid_mesh,1,2)

    mesh=mesh_model_surface(model,1)
    @test validate(mesh).ok
    @test mesh_crc(mesh).sha==
          "6ea713b4493eeb5b31e7c70ea5312290ef698424fd787bb04f1df4ef55f894cf"
    mapping=model_periodic_nodes(model,mesh,1,2)
    @test mapping.master_entity==4
    @test length(mapping.slave_nodes)==length(mapping.master_nodes)==7
    for (slave,master) in zip(mapping.slave_nodes,mapping.master_nodes)
        @test Tuple(mesh.coords[:,slave])==
              (mesh.coords[1,master]+1,mesh.coords[2,master],mesh.coords[3,master])
    end
    unsynchronized=deepcopy(mesh)
    unsynchronized.coords[3,mapping.slave_nodes[2]]+=1e-6
    @test validate(unsynchronized).ok
    @test_throws ArgumentError model_periodic_nodes(
        model,unsynchronized,1,2)
    mapping.slave_nodes[1]=1
    @test first(model_periodic_nodes(model,mesh,1,2).slave_nodes)!=1
    @test_throws ArgumentError model_periodic_nodes(model,mesh,1,4)
    @test_throws ArgumentError model_periodic_nodes(model,mesh,2,2)
    @test_throws ArgumentError mesh_model_surface(
        model,1;max_periodic_passes=false)
    @test_throws ArgumentError mesh_model_surface(
        model,1;max_periodic_passes=0)
    @test_throws ArgumentError mesh_model_surface(
        model,1;max_periodic_passes=65)
    @test_throws ErrorException mesh_model_surface(
        model,1;max_periodic_passes=1)

    # A surface containing only one entity from a relation is an explicit
    # blocker: silently meshing it would lose the required correspondence.
    add_point!(model,2,0,0;tag=5,mesh_size=0.25)
    add_point!(model,2,1,0;tag=6,mesh_size=0.25)
    add_line!(model,2,5;tag=5)
    add_line!(model,5,6;tag=6)
    add_line!(model,6,3;tag=7)
    add_curve_loop!(model,[5,6,7,-2];tag=2)
    add_plane_surface!(model,[2];tag=2)
    @test_throws ArgumentError mesh_model_surface(model,2)

    rotation_model=GeoModel()
    for (tag,(x,y)) in enumerate(((1.0,0.0),(2.0,0.0),
                                  (0.0,1.0),(0.0,2.0)))
        add_point!(rotation_model,x,y,0;tag=tag,mesh_size=0.5)
    end
    add_line!(rotation_model,1,2;tag=1)
    add_line!(rotation_model,2,4;tag=2)
    add_line!(rotation_model,3,4;tag=3)
    add_line!(rotation_model,3,1;tag=4)
    add_curve_loop!(rotation_model,[1,2,-3,4];tag=1)
    add_plane_surface!(rotation_model,[1];tag=1)
    rotation=(0.0,-1.0,0.0,0.0,
              1.0, 0.0,0.0,0.0,
              0.0, 0.0,1.0,0.0,
              0.0, 0.0,0.0,1.0)
    set_periodic!(rotation_model,1,[3],[1],rotation)
    @test !only(model_periodic_constraints(rotation_model)).reversed
    rotation_mesh=mesh_model_surface(rotation_model,1)
    @test mesh_crc(rotation_mesh).sha==
          "f6ad616e56d52d7e10a598a4079db2de9b3d5f2a777f492f5a2366946d8ea990"
    rotation_mapping=model_periodic_nodes(rotation_model,rotation_mesh,1,3)
    @test length(rotation_mapping.master_nodes)==3
    for (slave,master) in zip(rotation_mapping.slave_nodes,
                              rotation_mapping.master_nodes)
        @test Tuple(rotation_mesh.coords[:,slave])==
              (-rotation_mesh.coords[2,master],
                rotation_mesh.coords[1,master],
                rotation_mesh.coords[3,master])
    end

    double_periodic=periodic_square()
    set_periodic!(double_periodic,1,[2],[4],(
        1.0,0.0,0.0,1.0,
        0.0,1.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0))
    set_periodic!(double_periodic,1,[3],[1],(
        1.0,0.0,0.0,0.0,
        0.0,1.0,0.0,1.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0))
    double_mesh=mesh_model_surface(double_periodic,1)
    @test validate(double_mesh).ok
    @test mesh_crc(double_mesh).sha==
          "b82c9f0f4e235e90a754f2ec50b3a373ef0a2d514a79194d9b922873e35f8dd1"
    @test length(model_periodic_constraints(double_periodic))==2
    for slave_entity in (2,3)
        double_mapping=model_periodic_nodes(
            double_periodic,double_mesh,1,slave_entity)
        @test length(double_mapping.slave_nodes)==9
    end

    translate_y03=(1.0,0.0,0.0,0.0,
                   0.0,1.0,0.0,0.3,
                   0.0,0.0,1.0,0.0,
                   0.0,0.0,0.0,1.0)
    translate_y06=(1.0,0.0,0.0,0.0,
                   0.0,1.0,0.0,0.6,
                   0.0,0.0,1.0,0.0,
                   0.0,0.0,0.0,1.0)

    branch=periodic_curve_graph()
    set_periodic!(branch,1,[20],[30],translate_y03)
    set_periodic!(branch,1,[10],[30],translate_y06)
    branch_constraints=model_periodic_constraints(branch)
    @test Int.(getproperty.(branch_constraints,:slave_entity))==[10,20]
    @test Int.(getproperty.(branch_constraints,:master_entity))==[30,30]
    @test branch_constraints[1].reversed
    branch_mesh=mesh_model_surface(branch,1)
    @test validate(branch_mesh).ok
    @test mesh_crc(branch_mesh).sha==
          "dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338"
    for (slave_entity,offset) in ((10,0.6),(20,0.3))
        branch_mapping=model_periodic_nodes(
            branch,branch_mesh,1,slave_entity)
        @test branch_mapping.master_entity==30
        @test length(branch_mapping.slave_nodes)==9
        for (slave,master) in zip(branch_mapping.slave_nodes,
                                  branch_mapping.master_nodes)
            @test Tuple(branch_mesh.coords[:,slave])==
                  (branch_mesh.coords[1,master],
                   branch_mesh.coords[2,master]+offset,
                   branch_mesh.coords[3,master])
        end
    end

    chain=periodic_curve_graph()
    @test set_periodic!(
        chain,1,[20,10],[30,20],translate_y03)===nothing
    chain_constraints=model_periodic_constraints(chain)
    @test Int.(getproperty.(chain_constraints,:slave_entity))==[10,20]
    @test Int.(getproperty.(chain_constraints,:master_entity))==[20,30]
    cycle_error=try
        set_periodic!(chain,1,[30],[10],(
            1.0,0.0,0.0,0.0,
            0.0,1.0,0.0,-0.6,
            0.0,0.0,1.0,0.0,
            0.0,0.0,0.0,1.0))
        nothing
    catch err
        err
    end
    @test cycle_error isa ArgumentError
    @test occursin(
        "cyclic periodic dependency Curve[10] -> Curve[20] -> Curve[30] -> Curve[10]",
        sprint(showerror,cycle_error))
    @test model_periodic_constraints(chain)==chain_constraints
    chain_mesh=mesh_model_surface(chain,1)
    @test validate(chain_mesh).ok
    @test mesh_crc(chain_mesh).sha==
          "dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338"
    for (slave_entity,master_entity) in ((10,20),(20,30))
        chain_mapping=model_periodic_nodes(
            chain,chain_mesh,1,slave_entity)
        @test chain_mapping.master_entity==master_entity
        @test length(chain_mapping.slave_nodes)==9
        for (slave,master) in zip(chain_mapping.slave_nodes,
                                  chain_mapping.master_nodes)
            @test Tuple(chain_mesh.coords[:,slave])==
                  (chain_mesh.coords[1,master],
                   chain_mesh.coords[2,master]+0.3,
                   chain_mesh.coords[3,master])
        end
    end

    inconsistent=periodic_square()
    set_periodic!(inconsistent,1,[2],[4],(
        1.0,1e-4,0.0,1.0,
        0.0,1.0, 0.0,0.0,
        0.0,0.0, 1.0,0.0,
        0.0,0.0, 0.0,1.0);atol=1e-3)
    set_periodic!(inconsistent,1,[3],[1],(
        1.0,0.0,0.0,0.0,
        0.0,1.0,0.0,1.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0))
    @test_throws ErrorException mesh_model_surface(inconsistent,1)

    identity=(1.0,0.0,0.0,0.0,
              0.0,1.0,0.0,0.0,
              0.0,0.0,1.0,0.0,
              0.0,0.0,0.0,1.0)
    invalid=periodic_square()
    @test_throws ArgumentError set_periodic!(invalid,2,[2],[4],identity)
    @test_throws ArgumentError set_periodic!(invalid,false,[2],[4],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,Int[],Int[],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,[2],[4,1],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,[2,2],[4,1],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,[2],[9],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,[9],[4],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,[1],[2],identity)
    @test_throws ArgumentError set_periodic!(invalid,1,[2],[4],identity)
    @test_throws ArgumentError set_periodic!(
        invalid,1,[2],[4],zeros(4,4))
    @test_throws ArgumentError set_periodic!(
        invalid,1,[2],[4],identity;atol=true)
    @test_throws ArgumentError set_periodic!(
        invalid,1,[2],[4],identity;atol=-1)
    @test isempty(model_periodic_constraints(invalid))
    @test set_periodic!(invalid,1,[2],[4],(
        1.0,0.0,0.0,1.0,
        0.0,1.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0))===nothing
    @test_throws ArgumentError set_periodic!(invalid,1,[1],[3],(
        1.0,0.0,0.0,0.0,
        0.0,1.0,0.0,1.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0))
    @test length(model_periodic_constraints(invalid))==1

    degenerate=GeoModel()
    for tag in 1:4
        add_point!(degenerate,0,0,0;tag=tag)
    end
    add_line!(degenerate,1,2;tag=1)
    add_line!(degenerate,3,4;tag=2)
    @test_throws ArgumentError set_periodic!(
        degenerate,1,[2],[1],identity)
    @test isempty(model_periodic_constraints(degenerate))
    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end

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

    selective=mktemp() do path,io
        write(io, """
            /* A block comment is lexical whitespace. */
            Box(1) = {0, 0, 0, 1, 1, 1};
            Box(2) = {0.5, 0, 0, 1, 1, 1};
            BooleanUnion(3) = { Volume{1}; Delete; }{ Volume{2}; };
            """)
        close(io)
        execute_geo(path)
    end
    @test sort!(collect(keys(selective.model.volumes)))==[2,3]
    mktemp() do path,io
        write(io,"""
            Box(1) = {0, 0, 0, 1, 1, 1};
            Box(2) = {0.5, 0, 0, 1, 1, 1};
            BooleanUnion(3) = { Volume{1}; Remove; }{ Volume{2}; };
            """)
        close(io)
        @test_throws ArgumentError execute_geo(path)
    end

    mktemp() do path,io
        write(io, "Box(1) = {0, 0, 0, 1, 1, 1};\nBox(2) = {2, 0, 0, 1, 1, 1};\n")
        close(io)
        @test_throws ArgumentError execute_geo(path; mesh_dim=3)
        @test_throws ArgumentError execute_geo(path; mesh_dim=false)
        @test_throws ArgumentError execute_geo(path; mesh_dim=big(2)^100)
    end
    mktemp() do path,io
        write(io, "Extrude {0, 0, 1} { Surface{1}; };\n")
        close(io)
        @test_throws ArgumentError execute_geo(path)
    end
    mktemp() do path,io
        write(io,"/* unterminated\n"); close(io)
        @test_throws ArgumentError execute_geo(path)
    end
    mktemp() do path,io
        write(io,"Point(1) = {0, 0, 0, 1}};\n"); close(io)
        @test_throws ArgumentError execute_geo(path)
    end

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
    add_point!(emb,0.5,0.5,0; tag=5, mesh_size=0.25)
    add_line!(emb,1,2; tag=1); add_line!(emb,2,3; tag=2)
    add_line!(emb,3,4; tag=3); add_line!(emb,4,1; tag=4)
    add_curve_loop!(emb,[1,2,3,4]; tag=1)
    add_plane_surface!(emb,[1]; tag=1)
    embed!(emb,0,[5],2,1)
    emesh=mesh_model_surface(emb,1)
    @test validate(emesh).ok
    @test mesh_crc(emesh).sha==
          "13917dad18b19e8376640a379a2f1cd338aacca5b01f66e100b5bc372fc91371"
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

    linevol=GeoModel()
    add_box!(linevol,0,0,0,1,1,1; tag=1)
    add_point!(linevol,0.2,0.3,0.4; tag=1)
    add_point!(linevol,0.7,0.6,0.5; tag=2)
    add_line!(linevol,1,2; tag=1)
    embed!(linevol,1,[1],3,1)
    lv=mesh_model_volume(linevol,1)
    @test validate(lv).ok
    @test ntets(lv)>0
    LV=sum(tet_volume(node(lv,lv.tets[1,t]),node(lv,lv.tets[2,t]),
                      node(lv,lv.tets[3,t]),node(lv,lv.tets[4,t]))
           for t in 1:ntets(lv))
    @test LV≈1.0 atol=1e-12
    @test Tessella.Mesh3D.mesh_covers_segment3(lv,(0.2,0.3,0.4),(0.7,0.6,0.5))

    sheet=GeoModel()
    add_box!(sheet,0,0,0,1,1,1; tag=1)
    add_point!(sheet,0.2,0.2,0.5; tag=1)
    add_point!(sheet,0.8,0.2,0.5; tag=2)
    add_point!(sheet,0.5,0.8,0.5; tag=3)
    add_line!(sheet,1,2; tag=1); add_line!(sheet,2,3; tag=2); add_line!(sheet,3,1; tag=3)
    add_curve_loop!(sheet,[1,2,3]; tag=1)
    add_plane_surface!(sheet,[1]; tag=1)
    embed!(sheet,2,[1],3,1)
    sv=mesh_model_volume(sheet,1)
    @test validate(sv).ok
    @test ntets(sv)>0
    SV=sum(tet_volume(node(sv,sv.tets[1,t]),node(sv,sv.tets[2,t]),
                      node(sv,sv.tets[3,t]),node(sv,sv.tets[4,t]))
           for t in 1:ntets(sv))
    @test SV≈1.0 atol=1e-12
    @test Tessella.Mesh3D.mesh_covers_triangle3(sv,(0.2,0.2,0.5),(0.8,0.2,0.5),(0.5,0.8,0.5))

    sheetgeo=mktemp() do path,io
        write(io, """
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 1, 1, 1};
            Point(1) = {0.2, 0.2, 0.5, 0.5};
            Point(2) = {0.8, 0.2, 0.5, 0.5};
            Point(3) = {0.5, 0.8, 0.5, 0.5};
            Line(1) = {1, 2};
            Line(2) = {2, 3};
            Line(3) = {3, 1};
            Line Loop(1) = {1, 2, 3};
            Plane Surface(1) = {1};
            Surface{1} In Volume{1};
            """)
        close(io)
        execute_geo(path; mesh_dim=3)
    end
    @test validate(sheetgeo.mesh).ok
    @test Tessella.Mesh3D.mesh_covers_triangle3(sheetgeo.mesh,(0.2,0.2,0.5),(0.8,0.2,0.5),(0.5,0.8,0.5))
    SG=sum(tet_volume(node(sheetgeo.mesh,sheetgeo.mesh.tets[1,t]),
                      node(sheetgeo.mesh,sheetgeo.mesh.tets[2,t]),
                      node(sheetgeo.mesh,sheetgeo.mesh.tets[3,t]),
                      node(sheetgeo.mesh,sheetgeo.mesh.tets[4,t]))
           for t in 1:ntets(sheetgeo.mesh))
    @test SG≈1.0 atol=1e-12
end
