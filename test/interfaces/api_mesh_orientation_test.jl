using Test,Tessella
const _MBO=Tessella.API

function _mbo_cube(;occ=false,extra=true)
    _MBO.initialize();_MBO.option("Mesh.Renumber",0)
    if occ
        _MBO.model.add_box(0,0,0,1,1,1;tag=1)
    else
        for (tag,point) in enumerate(((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                                      (0,0,1),(1,0,1),(1,1,1),(0,1,1)))
            _MBO.model.add_point(point...;tag)
        end
        for (tag,edge) in enumerate(((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
                                     (7,8),(8,5),(1,5),(2,6),(3,7),(4,8)))
            _MBO.model.add_line(edge...;tag)
        end
        loops=([1,2,3,4],[5,6,7,8],[1,10,-5,-9],
               [2,11,-6,-10],[3,12,-7,-11],[4,9,-8,-12])
        for (tag,curves) in enumerate(loops)
            _MBO.model.add_curve_loop(curves;tag)
            _MBO.model.add_plane_surface([tag];tag)
        end
        _MBO.model.add_surface_loop([-1,2,3,4,5,6];tag=1)
        _MBO.model.add_volume([1];tag=1)
    end
    if extra
        first=occ ? 9 : 21;curve=occ ? 13 : 21
        for (offset,point) in enumerate(((3,0,0),(4,0,0),(4,1,0),(3,1,0)))
            _MBO.model.add_point(point...;tag=first+offset-1)
        end
        for offset in 0:3
            _MBO.model.add_line(first+offset,first+mod(offset+1,4);tag=curve+offset)
        end
        _MBO.model.add_curve_loop(collect(curve:curve+3);tag=7)
        _MBO.model.add_plane_surface([7];tag=7)
    end
end


@testset "Native Box face sense and parameter frame" begin
    try
        _MBO.initialize();_MBO.model.add_box(2,-3,5,2,3,4;tag=1)
        normals=([-1.,0,0],[1.,0,0],[0.,-1,0],[0.,1,0],[0.,0,-1],[0.,0,1])
        points=([2.,-2.25,8.],[4.,-2.25,8.],[3.5,-3.,6.],[3.5,0.,6.],[2.5,-.75,5.],[2.5,-.75,9.])
        parameters=([3.,-.75],[3.,-.75],[1.,1.5],[1.,1.5],[.5,2.25],[.5,2.25])
        for surface in 1:6
            @test _MBO.model.get_normal(surface,parameters[surface])==normals[surface]
            @test _MBO.model.get_parametrization(2,surface,points[surface])≈parameters[surface] atol=2e-11
            @test _MBO.model.get_value(2,surface,parameters[surface])≈points[surface] atol=2e-11
        end
        _MBO.mesh.set_outward_orientation(1)
        @test all(!get(_MBO.CURRENT[].meshing.reverse,(2,s),false) for s in 1:6)
        @test length(_MBO.CURRENT[].meshing.orientation_stl)==6
        @test isempty(_MBO.CURRENT[].meshing.outward_orientation)
    finally
        _MBO.finalize()
    end
end

@testset "Edited native P2 support row preservation" begin
    try
        _mbo_cube(;extra=false)
        for curve in 1:12;_MBO.mesh.set_transfinite_curve(curve,2);end
        for surface in 1:6;_MBO.mesh.set_transfinite_surface(surface);end
        _MBO.mesh.generate(2);_MBO.mesh.set_order(2);_MBO.mesh.renumber_nodes()
        @test _MBO.LAST_MESH[] isa Tessella.MixedMesh
        @test _MBO.LAST_MESH_HIGH_ORDER[]===nothing
        labels=_MBO.mesh.get_elements(2,1)[2][1]
        element=_MBO.mesh.get_element(first(labels))
        node=element[2][5]
        old=_MBO.mesh.get_node(node)
        coordinates=copy(old[1]);coordinates[3]=.07
        _MBO.mesh.set_node(node,coordinates,old[2])
        before=_MBO.LAST_MESH[]
        cells=Dict(label=>_MBO.mesh.get_element(label) for labels in
            _MBO.mesh.get_elements(2,-1)[2] for label in labels)
        points=copy(before.coords)
        node_labels=copy(_MBO.LAST_MESH_CLASS[].public_tags.node_tags)
        parameters=deepcopy(_MBO.LAST_MESH_CLASS[].node_parametric)
        _MBO.mesh.set_outward_orientation(1)
        after=_MBO.LAST_MESH[]
        @test after.coords==points
        @test _MBO.LAST_MESH_CLASS[].public_tags.node_tags==node_labels
        @test _MBO.LAST_MESH_CLASS[].node_parametric==parameters
        @test _MBO.mesh.get_node(node)[1]==coordinates
        for (label,cell) in cells
            actual=_MBO.mesh.get_element(label)
            @test actual[1]==cell[1]==9
            @test actual[3:4]==cell[3:4]
            @test actual[2]==(cell[4]==1 ? cell[2][[1,3,2,6,5,4]] : cell[2])
        end
    finally
        _MBO.finalize()
    end
end


@testset "Net-zero orientation retains explicit carrier identity" begin
    carrier=Tessella.Mesh([0. 1 0;0 0 1;0 0 0];tris=reshape(Int32[1,2,3],3,1))
    actual=Tessella.simplex_to_mixed(carrier)
    changes=Dict(1=>Int32[1,2,3])
    replacement=Tessella.Model._model_orientation_mixed_replacement(actual,changes)
    @test replacement===actual
    @test Tessella.Model._model_orientation_carrier(actual,replacement,carrier,Int32[1,2,3],changes)===carrier
    reversed=Tessella.Model._model_orientation_mixed_replacement(actual,Dict(1=>Int32[1,3,2]))
    rebound=Tessella.Model._model_orientation_carrier(actual,reversed,carrier,Int32[1,2,3],Dict(1=>Int32[1,3,2]))
    @test rebound.tris==reshape(Int32[1,3,2],3,1)
    @test rebound.coords==carrier.coords
end


@testset "Pinned Box endpoint collapse rejects before entity allocation" begin
    try
        _MBO.initialize();_MBO.model.add_point(0,0,0;tag=21)
        model=_MBO.CURRENT[]
        before=Dict(name=>deepcopy(getfield(model,name)) for name in fieldnames(typeof(model)))
        @test_throws ArgumentError _MBO.model.add_box(1e16,0,0,1,1,1;tag=1)
        @test _MBO.CURRENT[]===model
        @test all(name->getfield(model,name)==before[name],fieldnames(typeof(model)))
        @test _MBO.model.get_entities()==[(0,21)]
    finally
        _MBO.finalize()
    end
    # Actual primary admits the nearby representable widths. The rejected
    # OCC1e-9 absolute-tolerance case remains outside this exact collapse fix.
    for (origin,width) in ((1e16,4.0),(1e15,1.0))
        try
            _MBO.initialize();_MBO.model.add_box(origin,0,0,width,1,1;tag=1)
            @test length(_MBO.model.get_entities())==27
            for curve in 9:12
                lower,upper=_MBO.model.get_parametrization_bounds(1,curve)
                @test [only(lower),only(upper)]==[0.0,width]
                @test all(isfinite,_MBO.model.get_value(1,curve,[width/2]))
            end
        finally
            _MBO.finalize()
        end
    end
end
