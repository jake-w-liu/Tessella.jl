using Test
using Tessella
using Tessella.MeshTypes: mesh_crc
using Tessella.Model: model_entities, model_entities_for_physical_group,
                      model_entity_name, model_periodic_constraints,
                      model_physical_groups, model_physical_name,
                      remove_entities!, set_entity_name!

function _removal_two_triangles()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0))
        add_point!(model,x,y,0.0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,1),(4,2,4),(5,4,3))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3];tag=1)
    add_plane_surface!(model,[1];tag=1)
    add_curve_loop!(model,[4,5,-2];tag=2)
    add_plane_surface!(model,[2];tag=2)
    return model
end

function _removal_tetrahedron(;embedded=false)
    model=GeoModel()
    for (tag,x,y,z) in ((1,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (3,0.0,1.0,0.0),(4,0.0,0.0,1.0))
        add_point!(model,x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,1),(4,1,4),(5,2,4),(6,3,4))
        add_line!(model,first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((1,[1,2,3]),(2,[1,5,-4]),
                         (3,[2,6,-5]),(4,[3,4,-6]))
        add_curve_loop!(model,curves;tag=tag)
        add_plane_surface!(model,[tag];tag=tag)
    end
    add_surface_loop!(model,[1,-2,3,-4];tag=1)
    add_volume!(model,[1];tag=1)
    if embedded
        add_point!(model,0.2,0.2,0.2;tag=9)
        add_point!(model,0.3,0.2,0.2;tag=10)
        add_line!(model,9,10;tag=9)
        embed!(model,0,[9],3,1)
        embed!(model,1,[9],3,1)
    end
    return model
end

function _removal_embedded_surface()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,1.0,1.0),(4,0.0,1.0),
                      (5,0.25,0.25),(6,0.75,0.25))
        add_point!(model,x,y,0.0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),(5,5,6))
        add_line!(model,first_point,last_point;tag=tag)
    end
    add_curve_loop!(model,[1,2,3,4];tag=1)
    add_plane_surface!(model,[1];tag=1)
    embed!(model,0,[5],2,1)
    embed!(model,1,[5],2,1)
    return model
end

function _removal_periodic_curves()
    model=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,1.0,0.0),
                      (3,0.0,1.0),(4,1.0,1.0),
                      (5,0.0,2.0),(6,1.0,2.0))
        add_point!(model,x,y,0.0;tag=tag)
    end
    for tag in 1:3
        add_line!(model,2tag-1,2tag;tag=tag)
    end
    translate_down=(1.0,0.0,0.0,0.0,
                    0.0,1.0,0.0,-1.0,
                    0.0,0.0,1.0,0.0,
                    0.0,0.0,0.0,1.0)
    set_periodic!(model,1,[1,2],[2,3],translate_down)
    return model
end

function _removal_state(model)
    names=fieldnames(typeof(model))
    return NamedTuple{names}(map(name->deepcopy(getfield(model,name)),names))
end

@testset "ordered dependency-safe entity removal" begin
    model=_removal_two_triangles()
    retained_crc=mesh_crc(mesh_model_surface(model,2))
    initial=_removal_state(model)
    @test remove_entities!(model,[(0,1)])==0
    @test remove_entities!(model,[(1,1)])==0
    @test _removal_state(model)==initial

    @test remove_entities!(model,[(2,1)])==1
    @test model_entities(model,2)==[(2,2)]
    @test haskey(model.loops,1)
    @test remove_entities!(model,[(1,1)])==1
    @test !haskey(model.loops,1)
    @test remove_entities!(model,[(1,2)])==0
    @test mesh_crc(mesh_model_surface(model,2))==retained_crc

    low_first=_removal_two_triangles()
    @test remove_entities!(low_first,[(1,1),(2,1)])==1
    @test (1,1) in model_entities(low_first)
    @test !((2,1) in model_entities(low_first))

    high_first=_removal_two_triangles()
    @test remove_entities!(high_first,[(2,1),(1,1)])==2
    @test !((1,1) in model_entities(high_first))
    @test !((2,1) in model_entities(high_first))
    @test remove_entities!(high_first,[(2,99),(2,99)])==0
end

@testset "recursive explicit boundaries and embeddings" begin
    model=_removal_two_triangles()
    retained_crc=mesh_crc(mesh_model_surface(model,2))
    @test remove_entities!(model,[(2,1)],true)==4
    @test model_entities(model)==
          [(0,2),(0,3),(0,4),(1,2),(1,4),(1,5),(2,2)]
    @test !haskey(model.loops,1)
    @test haskey(model.loops,2)
    @test mesh_crc(mesh_model_surface(model,2))==retained_crc
    @test model.next_tag==[4,5,2,0]
    @test add_point!(model,2,2,0;tag=0)==5

    embedded=_removal_embedded_surface()
    before=_removal_state(embedded)
    @test remove_entities!(embedded,[(0,5),(1,5)])==0
    @test _removal_state(embedded)==before
    @test remove_entities!(embedded,[(2,1)],true)==9
    @test model_entities(embedded)==[(0,5),(0,6),(1,5)]
    @test isempty(embedded.embeds)

    source_first=_removal_embedded_surface()
    @test remove_entities!(source_first,[(0,5),(2,1)])==1
    @test (0,5) in model_entities(source_first)
    target_first=_removal_embedded_surface()
    @test remove_entities!(target_first,[(2,1),(0,5)])==1
    @test (0,5) in model_entities(target_first)
end

@testset "explicit volume recursive closure" begin
    model=_removal_tetrahedron()
    original_crc=mesh_crc(mesh_model_volume(model,1))
    for dimension in 0:3
        add_physical_group!(
            model,dimension,[first(last.(model_entities(model,dimension)))];
            tag=dimension+1,name="dimension $dimension")
    end
    for (dimension,tag) in model_entities(model)
        set_entity_name!(model,dimension,tag,"entity $dimension $tag")
    end
    nonrecursive=deepcopy(model)
    @test remove_entities!(nonrecursive,[(3,1)])==1
    @test isempty(model_entities(nonrecursive,3))
    @test length(model_entities(nonrecursive))==14
    @test isempty(model_physical_groups(nonrecursive,3))
    @test model_physical_name(nonrecursive,3,4)==""
    @test model_entity_name(nonrecursive,3,1)==""

    @test remove_entities!(model,[(3,1)],true)==15
    @test isempty(model_entities(model))
    @test isempty(model.loops)
    @test isempty(model.surface_loops)
    @test isempty(model.entity_names)
    @test isempty(model.physical)
    @test isempty(model.physical_names)
    @test model.next_tag==[4,6,4,1]

    embedded=_removal_tetrahedron(embedded=true)
    @test remove_entities!(embedded,[(3,1)],true)==15
    @test model_entities(embedded)==[(0,9),(0,10),(1,9)]
    @test isempty(embedded.embeds)
end

@testset "entity metadata and periodic cleanup" begin
    model=GeoModel()
    add_point!(model,0,0,0;tag=1)
    add_point!(model,1,0,0;tag=2)
    set_entity_name!(model,0,1,"first")
    set_entity_name!(model,0,2,"second")
    add_physical_group!(model,0,[1,2];tag=10,name="nodes")
    add_physical_group!(model,0,[1];tag=11,name="first only")
    @test remove_entities!(model,[(0,1)])==1
    @test model_entities_for_physical_group(model,0,10)==[2]
    @test model_physical_groups(model)==[(0,10)]
    @test model_physical_name(model,0,11)==""
    @test model_entity_name(model,0,1)==""
    @test remove_entities!(model,[(0,2)])==1
    @test isempty(model_physical_groups(model))
    @test model_physical_name(model,0,10)==""
    @test isempty(model.entity_names)

    periodic=_removal_periodic_curves()
    @test length(model_periodic_constraints(periodic))==2
    @test remove_entities!(periodic,[(1,2)])==1
    @test isempty(model_periodic_constraints(periodic))
    @test model_entities(periodic,1)==[(1,1),(1,3)]
end

@testset "primitive and Boolean removal ownership" begin
    model=GeoModel()
    add_box!(model,0,0,0,2,1,1;tag=1)
    add_box!(model,0,0,0,1,1,1;tag=2)
    boolean_volumes!(model,:difference,1,2;tag=3)
    result_crc=mesh_crc(mesh_model_volume(model,3))
    add_point!(model,1.5,0.5,0.5;tag=20)
    embed!(model,0,[20],3,1)
    add_physical_group!(model,3,[1,3];tag=10,name="solids")
    set_entity_name!(model,3,1,"object")
    set_entity_name!(model,3,3,"cut")

    @test remove_entities!(model,[(3,1)],true)==1
    @test !haskey(model.box_extents,1)
    @test !haskey(model.embeds,(3,1))
    @test model_entities_for_physical_group(model,3,10)==[3]
    @test model_entity_name(model,3,1)==""
    @test mesh_crc(mesh_model_volume(model,3))==result_crc

    @test remove_entities!(model,[(3,3)])==1
    @test !haskey(model.booleans,3)
    @test !haskey(model.boolean_operands,3)
    @test isempty(model_physical_groups(model))
    @test model_physical_name(model,3,10)==""

    add_cylinder!(model,10,0,0,0,0,1,1;tag=4)
    add_sphere!(model,20,0,0,1;tag=5)
    add_cone!(model,30,0,0,0,0,1,1,0;tag=6)
    @test remove_entities!(model,[(3,4),(3,5),(3,6)])==3
    @test isempty(model.cylinders)
    @test isempty(model.spheres)
    @test isempty(model.cones)
    @test model.next_tag[4]==6
    @test add_box!(model,40,0,0,1,1,1;tag=0)==7
end

@testset "removal validation and atomicity" begin
    calls=(
        model->remove_entities!(model,[(4,1)]),
        model->remove_entities!(model,[(0,0)]),
        model->remove_entities!(model,[(0,-1)]),
        model->remove_entities!(model,[(true,1)]),
        model->remove_entities!(model,[(0,true)]),
        model->remove_entities!(model,[(0,big(2)^100)]),
        model->remove_entities!(model,[(0,1,2)]),
        model->remove_entities!(model,[1]),
        model->remove_entities!(model,"bad"),
        model->remove_entities!(model,[(0,1)],1),
    )
    for call in calls
        model=_removal_two_triangles()
        before=_removal_state(model)
        @test_throws ArgumentError call(model)
        @test _removal_state(model)==before
    end

    corrupt=_removal_two_triangles()
    delete!(corrupt.loops,1)
    before=_removal_state(corrupt)
    err=try
        remove_entities!(corrupt,[(1,1)])
        nothing
    catch caught
        caught
    end
    @test err isa ErrorException
    @test occursin("references missing Loop",sprint(showerror,err))
    @test _removal_state(corrupt)==before

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
