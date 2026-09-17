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
    # Box face loops allocate from the dedicated curve-loop namespace
    # (Gmsh's _maxLineLoopNum), so the second box's curves continue at 13.
    # Boolean results materialize explicit OCC-style topology: this disjoint
    # union yields one Volume per solid — 3 (operand A) and 4 (operand B) —
    # each a six-face box whose boundary entities are fresh result entities.
    @test Tessella.Model.model_entities(primitive)==
          [Tuple{Int,Int}[(0,point) for point in 1:32];
           Tuple{Int,Int}[(1,curve) for curve in 1:48];
           Tuple{Int,Int}[(2,surface) for surface in 1:24];
           [(3,1),(3,2),(3,3),(3,4)]]
    @test Tessella.Model.model_dimension(primitive)==3
    @test Tessella.Model.model_boundary(primitive,[(3,1)],false,true,false)==
          [(2,-1),(2,2),(2,-3),(2,4),(2,-5),(2,6)]
    @test Tessella.Model.model_boundary(primitive,[(3,2)],false,true,false)==
          [(2,-7),(2,8),(2,-9),(2,10),(2,-11),(2,12)]
    @test Tessella.Model.model_boundary(primitive,[(3,1)],false,false,true)==
          Tuple{Int,Int}[(0,point) for point in 1:8]
    @test Tessella.Model.model_adjacencies(primitive,3,1)==
          (Int[],[1,2,3,4,5,6])
    @test Tessella.Model.model_boundary(primitive,[(3,3)],false,true,false)==
          [(2,-13),(2,14),(2,-15),(2,16),(2,-17),(2,18)]
    @test Tessella.Model.model_boundary(primitive,[(3,4)],false,true,false)==
          [(2,-19),(2,20),(2,-21),(2,22),(2,-23),(2,24)]
    @test Tessella.Model.model_boundary(primitive,[(3,3)],false,false,true)==
          Tuple{Int,Int}[(0,point) for point in 17:24]
    @test Tessella.Model.model_adjacencies(primitive,3,3)==
          (Int[],[13,14,15,16,17,18])
    @test Tessella.Model.model_adjacencies(primitive,3,4)==
          (Int[],[19,20,21,22,23,24])

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

# Every materialized Boolean face must store pcurves whose parameter-range
# endpoints evaluate to the uv of the curve's endpoint vertices on that face
# — a reparametrization slip here breaks nested Booleans and CurveOnSurface
# queries while leaving 3-D counts intact.
function _check_materialized_pcurves(model,first_surface)
    for (s,loops) in model.surfaces
        s<first_surface && continue
        sg=model.surface_geometry[s]
        uper=sg.occ in (:cylinder,:cone,:sphere) ? 2π : 0.0
        function uvdist(p,q)
            du=p[1]-q[1]
            uper>0 && (du=mod(du+π,2π)-π)
            hypot(du,p[2]-q[2])
        end
        for lt in loops, sc in model.loops[lt]
            c=abs(sc)
            haskey(sg.pcurves,c) || continue   # degenerate (pole) edges bind
                                               # their 2-D trace separately
            ent=sg.pcurves[c]
            pc=sc>0 ? ent.fwd : ent.rev
            @test pc!==nothing
            pc===nothing && continue
            cg=model.curve_geometry[c]
            a,b=model.curves[c]
            pa=Tessella.Model._brep_surface_uv(sg,model.points[a])
            pb=Tessella.Model._brep_surface_uv(sg,model.points[b])
            # stored pcurves are edge-parametrized: t0/t1 must land on the
            # uv of the edge's endpoint vertices — the nested-Boolean
            # reparametrization regression check
            @test uvdist(Tessella.Model._brep_pc_eval(pc,cg.t0),pa)<1e-9
            @test uvdist(Tessella.Model._brep_pc_eval(pc,cg.t1),pb)<1e-9
            # the mid-parameter must agree with the curve's 3-D midpoint —
            # endpoints alone cannot catch an interior warp
            crv,=Tessella.Model._brep_edge_curve(model,c)
            crv.kind===:degenerate && continue
            tm=(cg.t0+cg.t1)/2
            qm=Tessella.Model._brep_surface_uv(sg,
                Tessella.Model._brep_curve_eval(crv,tm))
            @test uvdist(Tessella.Model._brep_pc_eval(pc,tm),qm)<1e-9
        end
    end
end

function _result_curve_incidence(model,first_surface)
    incidence=Dict{Int,Int}()
    for (s,loops) in model.surfaces
        s<first_surface && continue
        for lt in loops, c in model.loops[lt]
            incidence[abs(c)]=get(incidence,abs(c),0)+1
        end
    end
    return incidence
end

@testset "materialized Boolean result topology" begin
    # overlapping fuse: OCC same-domain gluing yields a clean box boundary
    fuse=GeoModel()
    add_box!(fuse,0,0,0,1,1,1;tag=1)
    add_box!(fuse,0.5,0,0,1,1,1;tag=2)
    @test boolean_volumes!(fuse,:union,1,2;tag=3)==3
    @test sort(collect(keys(fuse.volumes)))==[1,2,3]
    # totals include both operand boxes (16/24/12); the materialized
    # result adds 8/12/6 — the same boundary OCC reports for the fuse
    @test length(fuse.points)==24
    @test length(fuse.curves)==36
    @test length(fuse.surfaces)==18
    @test length(fuse.volumes[3])==1
    @test all(==(2),values(_result_curve_incidence(fuse,13)))
    _check_materialized_pcurves(fuse,13)

    # corner fuse: an L-shaped region, two T-split merged faces, and the
    # internal coincident wall — matching Gmsh 4.15.2's 14/21/9 result
    corner=GeoModel()
    add_box!(corner,0,0,0,1,1,1;tag=1)
    add_box!(corner,1,0.5,0.5,0.5,0.5,0.5;tag=2)
    @test boolean_volumes!(corner,:union,1,2;tag=3)==3
    @test length(corner.points)==30
    @test length(corner.curves)==45
    @test length(corner.surfaces)==21
    @test all(==(2),values(_result_curve_incidence(corner,13)))
    _check_materialized_pcurves(corner,13)

    # corner cut keeps the imprint splits: 11/16/7 like Gmsh
    cut=GeoModel()
    add_box!(cut,0,0,0,1,1,1;tag=1)
    add_box!(cut,1,0.5,0.5,0.5,0.5,0.5;tag=2)
    @test boolean_volumes!(cut,:difference,1,2;tag=3)==3
    @test length(cut.points)==27
    @test length(cut.curves)==40
    @test length(cut.surfaces)==19
    @test all(==(2),values(_result_curve_incidence(cut,13)))
    _check_materialized_pcurves(cut,13)

    # zero-volume face contact intersects to nothing: OCC binds no output,
    # so the preallocated tag is rolled back and the call returns 0
    empty=GeoModel()
    add_box!(empty,0,0,0,1,1,1;tag=1)
    add_box!(empty,1,0.5,0.5,0.5,0.5,0.5;tag=2)
    @test boolean_volumes!(empty,:intersection,1,2;tag=3)==0
    @test sort(collect(keys(empty.volumes)))==[1,2]
    @test !haskey(empty.booleans,3)
    @test !haskey(empty.boolean_operands,3)
    @test !haskey(empty.boolean_components,3)

    # a slab cut splits the box into two disjoint solids: OCC returns one
    # volume per solid, linked through the component table
    multi=GeoModel()
    add_box!(multi,0,0,0,1,1,1;tag=1)
    add_box!(multi,0.4,-0.5,-0.5,0.2,2,2;tag=2)
    @test boolean_volumes!(multi,:difference,1,2;tag=3)==3
    @test sort(collect(keys(multi.volumes)))==[1,2,3,4]
    @test multi.boolean_components[3]==3
    @test multi.boolean_components[4]==3
    @test multi.booleans[4]==(op=:difference,a=1,b=2)
    @test all(==(2),values(_result_curve_incidence(multi,13)))

    # nested Booleans read the materialized boundary of the earlier result
    nested=GeoModel()
    add_box!(nested,0,0,0,1,1,1;tag=1)
    add_box!(nested,0.5,0,0,1,1,1;tag=2)
    @test boolean_volumes!(nested,:union,1,2;tag=3)==3
    add_box!(nested,0.25,0.25,-0.5,0.5,0.5,2;tag=4)
    @test boolean_volumes!(nested,:difference,3,4;tag=5)==5
    @test nested.booleans[5]==(op=:difference,a=3,b=4)
    @test all(==(2),values(_result_curve_incidence(nested,25)))
    _check_materialized_pcurves(nested,25)

    # curved sections: box minus a through-hole cylinder — Gmsh 4.15.2
    # produces 10 points, 15 curves, and 7 surfaces for the result
    curved=GeoModel()
    add_box!(curved,0,0,0,1,1,1;tag=1)
    add_cylinder!(curved,0.5,0.5,-0.5,0,0,2,0.3;tag=2)
    np=length(curved.points);nc=length(curved.curves);ns=length(curved.surfaces)
    @test boolean_volumes!(curved,:difference,1,2;tag=3)==3
    @test (length(curved.points)-np,length(curved.curves)-nc,
           length(curved.surfaces)-ns)==(10,15,7)
    @test all(==(2),values(_result_curve_incidence(curved,ns+1)))
    _check_materialized_pcurves(curved,ns+1)

    # unsupported intersections fail explicitly rather than approximating
    unsupported=GeoModel()
    add_cylinder!(unsupported,0,0,0,0,0,1,0.3;tag=1)
    add_cylinder!(unsupported,0.2,0,0,0,1,0,0.3;tag=2)
    @test_throws ArgumentError boolean_volumes!(unsupported,:union,1,2;tag=3)
    @test !haskey(unsupported.volumes,3)
    @test !haskey(unsupported.booleans,3)
end
