using Test
using Tessella
using Tessella.MeshTypes: bounding_box, mesh_crc, node, ntets, tet_volume, validate

function _snapshot_volume(mesh)
    validate(mesh).ok || error("Boolean snapshot test produced an invalid volume mesh")
    return sum(tet_volume(node(mesh,mesh.tets[1,t]),
                          node(mesh,mesh.tets[2,t]),
                          node(mesh,mesh.tets[3,t]),
                          node(mesh,mesh.tets[4,t]))
               for t in 1:ntets(mesh))
end

@testset "owned Boolean operands and volume deletion" begin
    model=GeoModel()
    add_box!(model,0,0,0,2,1,1;tag=1)
    add_box!(model,0,0,0,1,1,1;tag=2)
    @test boolean_volumes!(model,:difference,1,2;tag=3)==3
    @test model.booleans[3]==(op=:difference,a=1,b=2)
    @test haskey(model.boolean_operands,3)

    original=mesh_model_volume(model,3)
    original_crc=mesh_crc(original)
    @test _snapshot_volume(original)≈1.0 atol=1e-12
    @test original_crc.bbox==((1.0,0.0,0.0),(2.0,1.0,1.0))

    translate_volume!(model,1,(10,0,0))
    retained=mesh_model_volume(model,3)
    @test _snapshot_volume(retained)≈1.0 atol=1e-12
    @test mesh_crc(retained)==original_crc

    add_box!(model,1,0,0,0.5,1,1;tag=4)
    @test boolean_volumes!(model,:difference,3,4;tag=5)==5
    nested=mesh_model_volume(model,5)
    nested_crc=mesh_crc(nested)
    @test _snapshot_volume(nested)≈0.5 atol=1e-12
    @test nested_crc.bbox==((1.5,0.0,0.0),(2.0,1.0,1.0))

    add_box!(model,20,0,0,1,1,1;tag=6)
    add_point!(model,10.5,0.5,0.5;tag=90)
    add_point!(model,1.75,0.5,0.5;tag=91)
    embed!(model,0,[90],3,1)
    embed!(model,0,[91],3,3)
    @test add_physical_group!(model,3,[1,6];tag=10,name="mixed")==10
    @test add_physical_group!(model,3,[2];tag=11,name="tool")==11
    @test add_physical_group!(model,3,[3,6];tag=12,name="nested source")==12

    @test Tessella.Model._remove_volume_entity!(model,1)
    @test !haskey(model.volumes,1)
    @test !haskey(model.box_extents,1)
    @test !haskey(model.embeds,(3,1))
    @test model.physical[(3,10)]==[6]
    @test model.physical_names[(3,10)]=="mixed"

    @test Tessella.Model._remove_volume_entity!(model,2)
    @test !haskey(model.physical,(3,11))
    @test !haskey(model.physical_names,(3,11))
    @test !Tessella.Model._remove_volume_entity!(model,2)

    @test Tessella.Model._remove_volume_entity!(model,3)
    @test !haskey(model.booleans,3)
    @test !haskey(model.boolean_operands,3)
    @test !haskey(model.embeds,(3,3))
    @test model.physical[(3,12)]==[6]
    @test model.physical_names[(3,12)]=="nested source"

    translate_volume!(model,4,(10,0,0))
    add_sphere!(model,30,0,0,1;tag=1)
    add_box!(model,40,0,0,4,1,1;tag=2)
    add_cone!(model,50,0,0,0,0,1,1,0;tag=3)
    @test !haskey(model.box_extents,1)
    @test haskey(model.spheres,1)
    @test bounding_box(Tessella.Model._volume_surface(model,1))==
          ((29.0,-1.0,-1.0),(31.0,1.0,1.0))

    reused_source=mesh_model_volume(model,5)
    @test _snapshot_volume(reused_source)≈0.5 atol=1e-12
    @test mesh_crc(reused_source)==nested_crc

    stable_volumes=deepcopy(model.volumes)
    stable_booleans=copy(model.booleans)
    stable_operands=copy(model.boolean_operands)
    stable_next_tag=copy(model.next_tag)
    @test_throws ArgumentError boolean_volumes!(model,:difference,1,999;tag=7)
    @test_throws ArgumentError boolean_volumes!(model,:difference,1,2;tag=6)
    @test_throws ArgumentError boolean_volumes!(model,:xor,1,2;tag=7)
    @test model.volumes==stable_volumes
    @test model.booleans==stable_booleans
    @test model.boolean_operands==stable_operands
    @test model.next_tag==stable_next_tag

    damaged=deepcopy(model)
    delete!(damaged.boolean_operands,5)
    damaged_next_tag=copy(damaged.next_tag)
    err=try
        mesh_model_volume(damaged,5)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("Boolean Volume[5] has no owned operand snapshot",
                   sprint(showerror,err))
    err=try
        boolean_volumes!(damaged,:union,5,6;tag=7)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("recreate the Boolean in a fresh GeoModel",sprint(showerror,err))
    @test !haskey(damaged.volumes,7)
    @test !haskey(damaged.booleans,7)
    @test !haskey(damaged.boolean_operands,7)
    @test damaged.next_tag==damaged_next_tag
end

@testset "bounded .geo Boolean snapshots and selective Delete" begin
    execution=mktemp() do path,io
        write(io,"""
            SetFactory("OpenCASCADE");
            Box(1) = {0, 0, 0, 2, 1, 1};
            Box(2) = {0, 0, 0, 1, 1, 1};
            Point(90) = {0.5, 0.5, 0.5, 1};
            Point{90} In Volume{1};
            Physical Volume("retained", 10) = {1, 2};
            BooleanDifference(3) = {Volume{1}; Delete;}{Volume{2};};
            Translate {10, 0, 0} {Volume{2};};
            Sphere(1) = {30, 0, 0, 1};
            """)
        close(io)
        execute_geo(path)
    end
    @test sort!(collect(keys(execution.model.volumes)))==[1,2,3]
    @test !haskey(execution.model.box_extents,1)
    @test haskey(execution.model.spheres,1)
    @test !haskey(execution.model.embeds,(3,1))
    @test execution.model.physical[(3,10)]==[2]
    @test execution.model.physical_names[(3,10)]=="retained"
    result=mesh_model_volume(execution.model,3)
    @test _snapshot_volume(result)≈1.0 atol=1e-12
    @test mesh_crc(result).bbox==((1.0,0.0,0.0),(2.0,1.0,1.0))

    nested=mktemp() do path,io
        write(io,"""
            Box(1) = {0, 0, 0, 2, 1, 1};
            Box(2) = {0, 0, 0, 1, 1, 1};
            Box(4) = {1, 0, 0, 0.5, 1, 1};
            BooleanDifference(3) = {Volume{1}; Delete;}{Volume{2}; Delete;};
            BooleanDifference(5) = {Volume{3}; Delete;}{Volume{4}; Delete;};
            """)
        close(io)
        execute_geo(path;mesh_dim=3)
    end
    @test sort!(collect(keys(nested.model.volumes)))==[5]
    @test sort!(collect(keys(nested.model.booleans)))==[5]
    @test sort!(collect(keys(nested.model.boolean_operands)))==[5]
    @test _snapshot_volume(nested.mesh)≈0.5 atol=1e-12
    @test mesh_crc(nested.mesh).bbox==((1.5,0.0,0.0),(2.0,1.0,1.0))
end
