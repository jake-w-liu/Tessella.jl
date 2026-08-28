using Test
using Tessella
using Tessella.MeshTypes: mesh_crc
using Tessella.Model: model_boundary, model_entities,
                      model_entities_for_physical_group, model_entity_name,
                      model_periodic_constraints, model_set_tag!,
                      remove_entity_name!, set_entity_name!

function _identity_tetrahedron()
    model=GeoModel()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        add_point!(model,x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        add_line!(model,first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        add_curve_loop!(model,curves;tag=tag)
        add_plane_surface!(model,[tag];tag=tag)
    end
    add_surface_loop!(model,[21,-22,23,-24];tag=30)
    add_volume!(model,[30];tag=40)
    return model
end

function _identity_state(model)
    names=fieldnames(typeof(model))
    return NamedTuple{names}(map(name->deepcopy(getfield(model,name)),names))
end

function _identity_curve_chain()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0),
                      (5,0.0,2.0),(6,1.0,2.0))
        add_point!(model,x,y,0;tag=tag)
    end
    add_line!(model,1,2;tag=1)
    add_line!(model,3,4;tag=2)
    add_line!(model,5,6;tag=3)
    translate_down=(1.0,0.0,0.0,0.0,
                    0.0,1.0,0.0,-1.0,
                    0.0,0.0,1.0,0.0,
                    0.0,0.0,0.0,1.0)
    set_periodic!(model,1,[1,2],[2,3],translate_down)
    return model
end

function _identity_periodic_surfaces()
    model=GeoModel()
    coordinates=((0.0,0.0,0.0),(0.0,1.0,0.0),
                 (0.0,1.0,1.0),(0.0,0.0,1.0),
                 (1.0,0.0,0.0),(1.0,1.0,0.0),
                 (1.0,1.0,1.0),(1.0,0.0,1.0))
    for (tag,coordinate) in pairs(coordinates)
        add_point!(model,coordinate...;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),
             (5,5,6),(6,6,7),(7,7,8),(8,8,5))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_curve_loop!(model,[5,6,7,8];tag=2)
    add_plane_surface!(model,[1];tag=1)
    add_plane_surface!(model,[2];tag=2)
    translate_x=(1.0,0.0,0.0,1.0,
                 0.0,1.0,0.0,0.0,
                 0.0,0.0,1.0,0.0,
                 0.0,0.0,0.0,1.0)
    set_periodic!(model,2,[2],[1],translate_x)
    return model
end

@testset "owned model entity names" begin
    model=GeoModel()
    add_point!(model,0,0,0;tag=1)
    add_point!(model,1,0,0;tag=2)
    add_line!(model,1,2;tag=1)

    @test model_entity_name(model,0,1)==""
    @test model_entity_name(model,1,99)==""
    @test set_entity_name!(model,1,99,"missing")==""
    @test !haskey(model.entity_names,(1,99))
    @test set_entity_name!(model,0,1,"shared")=="shared"
    @test set_entity_name!(model,1,1,"shared")=="shared"
    @test model_entity_name(model,0,1)=="shared"
    @test model_entity_name(model,1,1)=="shared"
    @test set_entity_name!(model,1,1,"edge")=="edge"
    @test model_entity_name(model,1,1)=="edge"
    @test set_entity_name!(model,1,1,"")==""
    @test !haskey(model.entity_names,(1,1))
    set_entity_name!(model,0,1,"shared")
    set_entity_name!(model,0,2,"shared")
    set_entity_name!(model,1,1,"shared")
    @test remove_entity_name!(model,"shared")==3
    @test isempty(model.entity_names)
    @test remove_entity_name!(model,"shared")==0
    @test remove_entity_name!(model,"")==0

    stable=_identity_state(model)
    @test_throws ArgumentError model_entity_name(model,4,1)
    @test_throws ArgumentError model_entity_name(model,0,0)
    @test_throws ArgumentError model_entity_name(model,0,true)
    @test_throws ArgumentError set_entity_name!(model,0,1,7)
    @test_throws ArgumentError remove_entity_name!(model,7)
    @test _identity_state(model)==stable
end

@testset "atomic topology and metadata retagging" begin
    model=_identity_tetrahedron()
    add_physical_group!(model,0,[10,7];tag=101,name="vertices")
    add_physical_group!(model,1,[8,6];tag=102,name="edges")
    add_physical_group!(model,2,[21,22];tag=103,name="faces")
    add_physical_group!(model,3,[40];tag=104,name="domain")
    set_entity_name!(model,0,10,"origin")
    set_entity_name!(model,1,8,"base edge")
    set_entity_name!(model,1,6,"vertical edge")
    set_entity_name!(model,2,21,"base face")
    set_entity_name!(model,2,22,"side face")
    set_entity_name!(model,3,40,"tetrahedron")
    before=mesh_model_volume(model,40)
    before_crc=mesh_crc(before)
    before_tags=copy(before.tet_tag)

    @test model_set_tag!(model,0,10,110)==110
    @test model_set_tag!(model,1,8,108)==108
    @test model_set_tag!(model,1,6,106)==106
    @test model_set_tag!(model,2,21,121)==121
    @test model_set_tag!(model,2,22,122)==122
    @test model_set_tag!(model,3,40,140)==140

    @test model.curves[108]==(110,2)
    @test model.curves[11]==(7,110)
    @test model.curves[106]==(110,5)
    @test model.loops[21]==[108,3,11]
    @test model.loops[22]==[108,12,-106]
    @test model.loops[24]==[11,106,-4]
    @test model.surface_loops[30]==[121,-122,23,-24]
    @test model.volumes[140]==[30]
    @test model.point_size[110]==1.0
    @test model_entities_for_physical_group(model,0,101)==[7,110]
    @test model_entities_for_physical_group(model,1,102)==[106,108]
    @test model_entities_for_physical_group(model,2,103)==[121,122]
    @test model_entities_for_physical_group(model,3,104)==[140]
    @test model_entity_name(model,0,110)=="origin"
    @test model_entity_name(model,1,108)=="base edge"
    @test model_entity_name(model,1,106)=="vertical edge"
    @test model_entity_name(model,2,121)=="base face"
    @test model_entity_name(model,2,122)=="side face"
    @test model_entity_name(model,3,140)=="tetrahedron"
    @test model_entity_name(model,0,10)==""
    @test model_boundary(model,[(1,108)],false,true,false)==[(0,110),(0,2)]
    @test model_boundary(model,[(2,122)],false,true,false)==
          [(1,108),(1,12),(1,-106)]
    @test model_boundary(model,[(3,140)],false,true,false)==
          [(2,121),(2,-122),(2,23),(2,-24)]
    @test model.next_tag==[110,108,122,140]
    @test add_point!(model,2,2,2;tag=0)==111

    after=mesh_model_volume(model,140)
    @test mesh_crc(after)==before_crc
    @test after.tet_tag==before_tags
    @test !any(entity->entity in ((0,10),(1,8),(1,6),(2,21),(2,22),(3,40)),
               model_entities(model))
end

@testset "embedding source and target retagging" begin
    model=_identity_tetrahedron()
    add_point!(model,0.2,0.2,0.0;tag=30)
    add_point!(model,0.4,0.2,0.0;tag=31)
    add_line!(model,30,31;tag=30)
    embed!(model,0,[30],2,21)
    embed!(model,1,[30],2,21)
    embed!(model,2,[22],3,40)
    embed!(model,0,[31],3,40)

    model_set_tag!(model,0,30,130)
    model_set_tag!(model,1,30,130)
    model_set_tag!(model,2,21,121)
    model_set_tag!(model,2,22,122)
    model_set_tag!(model,3,40,140)
    @test model.embeds[(2,121)]==[(0,130),(1,130)]
    @test model.embeds[(3,140)]==[(2,122),(0,31)]
    @test !haskey(model.embeds,(2,21))
    @test !haskey(model.embeds,(3,40))
end

@testset "periodic identity follows curve and surface tags" begin
    curves=_identity_curve_chain()
    add_physical_group!(curves,1,[1,2,3];tag=10,name="periodic curves")
    for tag in 1:3
        set_entity_name!(curves,1,tag,"curve $tag")
    end
    model_set_tag!(curves,1,2,102)
    constraints=model_periodic_constraints(curves)
    @test [(Int(c.slave_entity),Int(c.master_entity)) for c in constraints]==
          [(1,102),(102,3)]
    model_set_tag!(curves,1,1,101)
    model_set_tag!(curves,1,3,103)
    constraints=model_periodic_constraints(curves)
    @test [(Int(c.slave_entity),Int(c.master_entity)) for c in constraints]==
          [(101,102),(102,103)]
    @test sort!(collect(keys(curves.periodic)))==[(1,101),(1,102)]
    @test model_entities_for_physical_group(curves,1,10)==[101,102,103]
    @test [model_entity_name(curves,1,tag) for tag in 101:103]==
          ["curve 1","curve 2","curve 3"]

    surfaces=_identity_periodic_surfaces()
    add_physical_group!(surfaces,2,[1,2];tag=20,name="periodic faces")
    set_entity_name!(surfaces,2,1,"master")
    set_entity_name!(surfaces,2,2,"slave")
    model_set_tag!(surfaces,2,1,101)
    model_set_tag!(surfaces,2,2,102)
    constraint=only(model_periodic_constraints(surfaces))
    @test (Int(constraint.slave_entity),Int(constraint.master_entity))==(102,101)
    @test only(keys(surfaces.periodic))==(2,102)
    @test model_entities_for_physical_group(surfaces,2,20)==[101,102]
    @test model_entity_name(surfaces,2,101)=="master"
    @test model_entity_name(surfaces,2,102)=="slave"
end

@testset "primitive and Boolean volume identity" begin
    model=GeoModel()
    add_box!(model,0,0,0,2,1,1;tag=1)
    add_box!(model,0,0,0,1,1,1;tag=2)
    boolean_volumes!(model,:difference,1,2;tag=3)
    add_point!(model,1.5,0.5,0.5;tag=10)
    embed!(model,0,[10],3,1)
    add_physical_group!(model,3,[1,3];tag=10,name="solids")
    set_entity_name!(model,3,1,"object")
    set_entity_name!(model,3,2,"tool")
    set_entity_name!(model,3,3,"cut")
    primitive_before=mesh_model_volume(model,1)
    result_before=mesh_model_volume(model,3)

    model_set_tag!(model,3,1,101)
    @test haskey(model.box_extents,101) && !haskey(model.box_extents,1)
    @test model.embeds[(3,101)]==[(0,10)]
    @test model.booleans[3]==(op=:difference,a=1,b=2)
    @test model_entity_name(model,3,101)=="object"
    @test mesh_crc(mesh_model_volume(model,101))==mesh_crc(primitive_before)
    @test mesh_crc(mesh_model_volume(model,3))==mesh_crc(result_before)

    model_set_tag!(model,3,3,103)
    @test haskey(model.booleans,103) && !haskey(model.booleans,3)
    @test haskey(model.boolean_operands,103) && !haskey(model.boolean_operands,3)
    @test model.booleans[103]==(op=:difference,a=1,b=2)
    @test model_entity_name(model,3,103)=="cut"
    @test mesh_crc(mesh_model_volume(model,103))==mesh_crc(result_before)

    model_set_tag!(model,3,2,102)
    @test haskey(model.box_extents,102) && !haskey(model.box_extents,2)
    @test model.booleans[103]==(op=:difference,a=1,b=2)
    @test model_entities_for_physical_group(model,3,10)==[101,103]
    @test model.next_tag[4]==103
    @test add_sphere!(model,20,0,0,1;tag=0)==104
end

@testset "retag validation and failure atomicity" begin
    model=GeoModel()
    add_point!(model,0,0,0;tag=1)
    add_point!(model,1,0,0;tag=2)
    add_line!(model,1,2;tag=1)
    calls=(
        ()->model_set_tag!(model,true,1,3),
        ()->model_set_tag!(model,0,true,3),
        ()->model_set_tag!(model,0,1,true),
        ()->model_set_tag!(model,0,0,3),
        ()->model_set_tag!(model,0,1,0),
        ()->model_set_tag!(model,0,99,3),
        ()->model_set_tag!(model,0,1,2),
        ()->model_set_tag!(model,0,1,1),
        ()->model_set_tag!(model,0,1,big(2)^100),
    )
    for call in calls
        before=_identity_state(model)
        @test_throws ArgumentError call()
        @test _identity_state(model)==before
    end

    corrupt_name=deepcopy(model)
    corrupt_name.entity_names[(0,9)]="stale"
    before=_identity_state(corrupt_name)
    @test_throws ArgumentError model_set_tag!(corrupt_name,0,1,9)
    @test _identity_state(corrupt_name)==before

    corrupt_boolean=GeoModel()
    add_box!(corrupt_boolean,0,0,0,2,1,1;tag=1)
    add_box!(corrupt_boolean,0,0,0,1,1,1;tag=2)
    boolean_volumes!(corrupt_boolean,:difference,1,2;tag=3)
    delete!(corrupt_boolean.boolean_operands,3)
    before=_identity_state(corrupt_boolean)
    err=try
        model_set_tag!(corrupt_boolean,3,3,4)
        nothing
    catch caught
        caught
    end
    @test err isa ErrorException
    @test occursin("inconsistent snapshot ownership",sprint(showerror,err))
    @test _identity_state(corrupt_boolean)==before

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
