using Test
using Tessella

function _topology_query_tetrahedron()
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

function _add_topology_query_tetra_shell!(model,offset,scale)
    coordinates=((0.0,0.0,0.0),(scale,0.0,0.0),
                 (0.0,scale,0.0),(0.0,0.0,scale))
    for (index,point) in pairs(coordinates)
        add_point!(model,point...;tag=offset+index)
    end
    for (index,(first_point,last_point)) in
            pairs(((1,2),(2,3),(3,1),(1,4),(2,4),(3,4)))
        add_line!(model,offset+first_point,offset+last_point;tag=offset+index)
    end
    for (index,curves) in pairs(((1,2,3),(1,5,-4),(2,6,-5),(3,4,-6)))
        signed_curves=Int[sign(curve)*(offset+abs(curve)) for curve in curves]
        add_curve_loop!(model,signed_curves;tag=offset+index)
        add_plane_surface!(model,[offset+index];tag=offset+index)
    end
    add_surface_loop!(
        model,[offset+1,-(offset+2),offset+3,-(offset+4)];tag=offset+1)
    return offset+1
end

@testset "deterministic explicit model topology queries" begin
    empty_model=GeoModel()
    @test Tessella.Model.model_entities(empty_model)==Tuple{Int,Int}[]
    @test Tessella.Model.model_entities(empty_model,2)==Tuple{Int,Int}[]
    @test Tessella.Model.model_dimension(empty_model)==-1
    @test Tessella.Model.model_boundary(empty_model,[])==Tuple{Int,Int}[]

    model=_topology_query_tetrahedron()
    expected_entities=[
        (0,2),(0,5),(0,7),(0,10),
        (1,3),(1,4),(1,6),(1,8),(1,11),(1,12),
        (2,21),(2,22),(2,23),(2,24),(3,40)]
    @test Tessella.Model.model_entities(model)==expected_entities
    @test Tessella.Model.model_entities(model,0)==
          [(0,2),(0,5),(0,7),(0,10)]
    @test Tessella.Model.model_entities(model,1)==
          [(1,3),(1,4),(1,6),(1,8),(1,11),(1,12)]
    @test Tessella.Model.model_entities(model,2)==
          [(2,21),(2,22),(2,23),(2,24)]
    @test Tessella.Model.model_entities(model,3)==[(3,40)]
    @test Tessella.Model.model_dimension(model)==3
    detached_entities=Tessella.Model.model_entities(model)
    push!(detached_entities,(3,99))
    @test Tessella.Model.model_entities(model)==expected_entities

    @test Tessella.Model.model_boundary(model,[(0,10)],false,false,false)==
          Tuple{Int,Int}[]
    @test Tessella.Model.model_boundary(model,[(0,10)],false,false,true)==
          [(0,10)]
    @test Tessella.Model.model_boundary(model,[(0,-10)],true,true,true)==
          [(0,10)]
    @test Tessella.Model.model_boundary(model,[(1,8)],false,false,false)==
          [(0,10),(0,2)]
    @test Tessella.Model.model_boundary(model,[1=>8],false,true,true)==
          [(0,10),(0,2)]
    @test Tessella.Model.model_boundary(model,[(1,-8)],false,false,false)==
          [(0,2),(0,10)]
    @test Tessella.Model.model_boundary(model,[(1,8)],true,false,false)==
          [(0,2),(0,10)]

    @test Tessella.Model.model_boundary(model,[(2,22)],false,false,false)==
          [(1,8),(1,12),(1,6)]
    @test Tessella.Model.model_boundary(model,[(2,22)],false,true,false)==
          [(1,8),(1,12),(1,-6)]
    @test Tessella.Model.model_boundary(model,[(2,-22)],false,true,false)==
          [(1,8),(1,12),(1,-6)]
    @test Tessella.Model.model_boundary(model,[(2,22)],true,true,false)==
          [(1,-6),(1,8),(1,12)]
    @test Tessella.Model.model_boundary(model,[(2,22)],false,false,true)==
          [(0,2),(0,5),(0,10)]

    two_surfaces=[(2,21),(2,22)]
    @test Tessella.Model.model_boundary(model,two_surfaces,false,false,false)==
          [(1,8),(1,3),(1,11),(1,8),(1,12),(1,6)]
    @test Tessella.Model.model_boundary(model,two_surfaces,false,true,false)==
          [(1,8),(1,3),(1,11),(1,8),(1,12),(1,-6)]
    @test Tessella.Model.model_boundary(model,two_surfaces,true,false,false)==
          [(1,3),(1,6),(1,11),(1,12)]
    @test Tessella.Model.model_boundary(model,two_surfaces,true,true,false)==
          [(1,3),(1,-6),(1,11),(1,12)]
    @test Tessella.Model.model_boundary(model,two_surfaces,false,false,true)==
          [(0,2),(0,7),(0,10),(0,2),(0,5),(0,10)]
    @test Tessella.Model.model_boundary(model,two_surfaces,true,false,true)==
          [(0,5),(0,7)]
    @test Tessella.Model.model_boundary(
        model,[(1,8),(2,21)],true,false,true)==[(0,7)]
    @test isempty(Tessella.Model.model_boundary(
        model,[(2,21),(2,21)],true,true,false))
    @test Tessella.Model.model_boundary(
        model,[(2,21),(2,21),(2,21)],true,true,false)==
        [(1,3),(1,8),(1,11)]

    @test Tessella.Model.model_boundary(model,[(3,40)],false,false,false)==
          [(2,21),(2,22),(2,23),(2,24)]
    @test Tessella.Model.model_boundary(model,[(3,40)],false,true,false)==
          [(2,21),(2,-22),(2,23),(2,-24)]
    @test Tessella.Model.model_boundary(model,[(3,-40)],false,true,false)==
          [(2,21),(2,-22),(2,23),(2,-24)]
    @test Tessella.Model.model_boundary(model,[(3,40)],false,false,true)==
          [(0,2),(0,5),(0,7),(0,10)]
    detached_boundary=Tessella.Model.model_boundary(
        model,[(2,21),(2,22)],true,true,false)
    push!(detached_boundary,(1,99))
    @test Tessella.Model.model_boundary(
        model,[(2,21),(2,22)],true,true,false)==
        [(1,3),(1,-6),(1,11),(1,12)]

    @test Tessella.Model.model_adjacencies(model,0,10)==([6,8,11],Int[])
    @test Tessella.Model.model_adjacencies(model,0,2)==([3,8,12],Int[])
    @test Tessella.Model.model_adjacencies(model,1,8)==([21,22],[10,2])
    @test Tessella.Model.model_adjacencies(model,1,6)==([22,24],[10,5])
    @test Tessella.Model.model_adjacencies(model,2,21)==([40],[8,3,11])
    @test Tessella.Model.model_adjacencies(model,2,22)==([40],[8,12,6])
    @test Tessella.Model.model_adjacencies(model,3,40)==(Int[],[21,22,23,24])
    detached_upward,detached_downward=
        Tessella.Model.model_adjacencies(model,1,8)
    push!(detached_upward,99);push!(detached_downward,99)
    @test Tessella.Model.model_adjacencies(model,1,8)==([21,22],[10,2])

    embedded=deepcopy(model)
    for (tag,x,y,z) in ((30,0.2,0.2,0.0),(31,0.4,0.2,0.0),
                        (32,0.2,0.4,0.0))
        add_point!(embedded,x,y,z;tag=tag)
    end
    add_line!(embedded,30,31;tag=30)
    embed!(embedded,0,[32],2,21)
    embed!(embedded,1,[30],2,21)
    @test Tessella.Model.model_adjacencies(embedded,0,30)==([30],Int[])
    @test Tessella.Model.model_adjacencies(embedded,0,32)==(Int[],Int[])
    @test Tessella.Model.model_adjacencies(embedded,1,30)==(Int[],[30,31])
    @test Tessella.Model.model_adjacencies(embedded,2,21)==([40],[8,3,11])

    holed=GeoModel()
    for (tag,x,y) in ((1,0.0,0.0),(2,2.0,0.0),(3,2.0,2.0),(4,0.0,2.0),
                      (5,0.5,0.5),(6,1.5,0.5),(7,1.5,1.5),(8,0.5,1.5))
        add_point!(holed,x,y,0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((1,1,2),(2,2,3),(3,3,4),(4,4,1),
             (5,5,6),(6,6,7),(7,7,8),(8,8,5))
        add_line!(holed,first_point,last_point;tag=tag)
    end
    add_curve_loop!(holed,[1,2,3,4];tag=1)
    add_curve_loop!(holed,[5,6,7,8];tag=2)
    add_plane_surface!(holed,[1,2];tag=1)
    @test Tessella.Model.model_boundary(holed,[(2,1)],false,true,false)==
          [(1,1),(1,2),(1,3),(1,4),(1,-8),(1,-7),(1,-6),(1,-5)]
    @test Tessella.Model.model_adjacencies(holed,2,1)==
          (Int[],[1,2,3,4,8,7,6,5])

    cavity=GeoModel()
    outer=_add_topology_query_tetra_shell!(cavity,0,2.0)
    inner=_add_topology_query_tetra_shell!(cavity,100,1.0)
    add_volume!(cavity,[outer,inner];tag=1)
    @test Tessella.Model.model_boundary(cavity,[(3,1)],false,true,false)==
          [(2,1),(2,-2),(2,3),(2,-4),
           (2,-101),(2,102),(2,-103),(2,104)]
    @test Tessella.Model.model_boundary(cavity,[(3,1)],false,false,true)==
          [(0,1),(0,2),(0,3),(0,4),(0,101),(0,102),(0,103),(0,104)]
    @test Tessella.Model.model_adjacencies(cavity,3,1)==
          (Int[],[1,2,3,4,101,102,103,104])

    primitive=GeoModel()
    add_box!(primitive,0,0,0,1,1,1;tag=1)
    add_box!(primitive,2,0,0,1,1,1;tag=2)
    boolean_volumes!(primitive,:union,1,2;tag=3)
    @test Tessella.Model.model_entities(primitive)==[(3,1),(3,2),(3,3)]
    @test Tessella.Model.model_dimension(primitive)==3
    @test_throws ArgumentError Tessella.Model.model_boundary(primitive,[(3,1)])
    @test_throws ArgumentError Tessella.Model.model_boundary(primitive,[(3,3)])
    @test_throws ArgumentError Tessella.Model.model_adjacencies(primitive,3,1)

    @test_throws ArgumentError Tessella.Model.model_entities(model,true)
    @test_throws ArgumentError Tessella.Model.model_entities(model,4)
    @test_throws ArgumentError Tessella.Model.model_boundary(model,nothing)
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[1])
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(4,1)])
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(1,0)])
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(1,true)])
    @test_throws ArgumentError Tessella.Model.model_boundary(
        model,[(1,-(big(2)^100))])
    @test_throws ArgumentError Tessella.Model.model_boundary(
        model,[(1,8)],1,false,false)
    @test_throws ArgumentError Tessella.Model.model_boundary(
        model,[(1,8)],true,1,false)
    @test_throws ArgumentError Tessella.Model.model_boundary(
        model,[(1,8)],true,false,1)
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(0,99)])
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(1,99)])
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(2,99)])
    @test_throws ArgumentError Tessella.Model.model_boundary(model,[(3,99)])
    @test_throws ArgumentError Tessella.Model.model_adjacencies(model,true,1)
    @test_throws ArgumentError Tessella.Model.model_adjacencies(model,1,0)
    @test_throws ArgumentError Tessella.Model.model_adjacencies(model,1,99)
    @test_throws ArgumentError Tessella.Model.model_adjacencies(
        model,1,big(2)^100)

    @test isempty(Docs.undocumented_names(Tessella.Model;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.Model;recursive=true))
end
