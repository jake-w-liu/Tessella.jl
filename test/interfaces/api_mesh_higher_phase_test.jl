using Test, Tessella

const _SCP=Tessella.API

@testset "Completed native Curve phase supplies actual generation-two input" begin
    # Independent 4.15.2 rows: square_generate1_to2_primary_v1.json,
    # SHA256 36485c10b3da5d8e2535c85e6e43a09e1b2fcd8896e7323ed897cc2f81f5cd1f.
    for only_empty in (0,1),order in (1,2)
        try
            _SCP.initialize();_SCP.option("Mesh.Renumber",0)
            _SCP.option("Mesh.MeshOnlyEmpty",only_empty)
            for (x,y) in ((0,0),(1,0),(1,1),(0,1))
                _SCP.model.add_point(x,y,0)
            end
            for (a,b) in ((1,2),(2,3),(3,4),(4,1));_SCP.model.add_line(a,b);end
            _SCP.model.add_curve_loop([1,2,3,4]);_SCP.model.add_plane_surface([1])
            for curve in 1:4;_SCP.mesh.set_transfinite_curve(curve,3);end
            _SCP.mesh.set_transfinite_surface(1)
            _SCP.mesh.generate(1)
            points=Dict(point=>_SCP.mesh.get_nodes(0,point) for point in 1:4)
            primaries=Dict(curve=>only(_SCP.mesh.get_nodes(1,curve)[1]) for curve in 1:4)
            @test primaries==Dict(curve=>UInt64(curve+4) for curve in 1:4)
            authority=_SCP.CURRENT[].meshing.mesh_authority
            @test Tessella.Model._model_mesh_status(_SCP.CURRENT[],authority)==1
            @test _SCP.mesh.get_elements_by_type(1)[1]==UInt64.(5:12)
            _SCP.option("Mesh.ElementOrder",order);_SCP.mesh.generate(2)
            @test Dict(point=>_SCP.mesh.get_nodes(0,point) for point in 1:4)==points
            @test all(curve->primaries[curve] in _SCP.mesh.get_nodes(1,curve)[1],1:4)
            line_type=order==1 ? 1 : 8;line_start=order==1 ? 13 : 29
            @test _SCP.mesh.get_elements_by_type(line_type)[1]==UInt64.(line_start:line_start+7)
            for curve in 1:4
                endpoints=(curve,mod1(curve+1,4));primary=curve+4
                expected=order==1 ? [endpoints[1],primary,primary,endpoints[2]] :
                    [endpoints[1],primary,8+2curve,primary,endpoints[2],9+2curve]
                @test _SCP.mesh.get_elements_by_type(line_type,curve)[2]==UInt64.(expected)
            end
            triangle_type=order==1 ? 2 : 9;triangle_start=order==1 ? 21 : 37
            @test _SCP.mesh.get_elements_by_type(triangle_type,1)[1]==UInt64.(triangle_start:triangle_start+7)
            @test length(_SCP.mesh.get_nodes()[1])==(order==1 ? 9 : 25)
            @test sum(_SCP.mesh.get_element_qualities(UInt64.(triangle_start:triangle_start+7),"volume"))≈1 atol=2e-11
        finally
            _SCP.finalize()
        end
    end
end

# Queued higher-phase/global-order parity, independently observed in Gmsh 4.15.2.
# public_labels_square_regeneration_primary_v2.json SHA256
# 655558baa12280d45202480e2c0df42485647a505b8401e466001b593fdc5442.
# This fixture intentionally remains separate from the accepted label allocator
# controls: native DONE Curve products survive a repeated Surface phase, while
# the global order pass elevates an unrelated fully discrete Curve.
# Element comparisons use each implementation's historical maxima because the
# legacy native Surface cache does not publish primary's lower-dimensional cells.

@testset "Repeated native Square phase retains Curve primaries and orders discrete Curve" begin
    try
        _SCP.initialize();_SCP.option("Mesh.Renumber",0)
        for (x,y) in ((0,0),(1,0),(1,1),(0,1))
            _SCP.model.add_point(x,y,0)
        end
        for (a,b) in ((1,2),(2,3),(3,4),(4,1))
            _SCP.model.add_line(a,b)
        end
        _SCP.model.add_curve_loop([1,2,3,4]);_SCP.model.add_plane_surface([1])
        _SCP.model.add_physical_group(2,[1];tag=99)
        for curve in 1:4;_SCP.mesh.set_transfinite_curve(curve,3);end
        _SCP.mesh.set_transfinite_surface(1)
        _SCP.mesh.generate(2)
        # The standing virgin cache uses its original grid row labels. Establish
        # the independent primary owner labels through an explicit public edit;
        # this leaves the virgin mesh and its CRC expectations unchanged.
        primaries=_SCP.mesh.get_nodes()[1]
        sort!(primaries;by=tag->begin
            value=_SCP.mesh.get_node(tag)
            (value[3],value[4])
        end,alg=Base.Sort.MergeSort)
        _SCP.mesh.renumber_nodes(primaries,UInt64.(1:length(primaries)))
        _SCP.mesh.set_order(2)
        old_elements=vcat(_SCP.mesh.get_elements()[2]...)
        _SCP.mesh.renumber_elements(old_elements,old_elements.+1000)
        _SCP.option("Mesh.ElementOrder",2)
        _SCP.model.add_discrete_entity(1,90)
        _SCP.mesh.add_nodes(1,90,[6001,6002],[10.,0,0,11,0,0])
        _SCP.mesh.add_elements_by_type(90,1,[901],[6001,6002])
        before_points=Dict(point=>sort(_SCP.mesh.get_nodes(0,point)[1]) for point in 1:4)
        @test before_points==Dict(point=>UInt64[point] for point in 1:4)
        @test Dict(curve=>minimum(_SCP.mesh.get_nodes(1,curve)[1]) for curve in 1:4)==
            Dict(curve=>UInt64(curve+4) for curve in 1:4)
        @test _SCP.mesh.get_element(901)==(1,UInt64[6001,6002],1,90)
        node_highwater=_SCP.mesh.get_max_node_tag()
        element_highwater=_SCP.mesh.get_max_element_tag()
        @test node_highwater==6002
        _SCP.mesh.generate(2)

        @test Dict(point=>sort(_SCP.mesh.get_nodes(0,point)[1]) for point in 1:4)==before_points
        expected_curves=Dict(1=>UInt64[5,6004,6005],2=>UInt64[6,6006,6007],
            3=>UInt64[7,6008,6009],4=>UInt64[8,6010,6011])
        @test Dict(curve=>sort(_SCP.mesh.get_nodes(1,curve)[1]) for curve in 1:4)==expected_curves
        @test sort(_SCP.mesh.get_nodes(2,1)[1])==UInt64[6003,6013,6014,6015,6016,6017,6018,6019,6020]
        native_nodes=vcat([_SCP.mesh.get_nodes(dim,entity)[1]
            for (dim,entity) in ((0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(1,3),(1,4),(2,1))]...)
        @test minimum(filter(tag->tag>node_highwater,native_nodes))==node_highwater+1
        types,surface_elements,surface_connections=_SCP.mesh.get_elements(2,1)
        @test types==[9]
        @test sum(length,surface_elements)==8
        @test all(tag->tag>element_highwater,Iterators.flatten(surface_elements))
        foreign_tags,foreign_nodes=_SCP.mesh.get_elements_by_type(8,90)
        @test length(foreign_tags)==1
        @test only(foreign_tags)>element_highwater
        @test foreign_nodes==UInt64[6001,6002,6012]
        @test _SCP.mesh.get_element(only(foreign_tags))==(8,UInt64[6001,6002,6012],1,90)
        @test sort(_SCP.mesh.get_nodes(1,90)[1])==UInt64[6001,6002,6012]
        @test _SCP.mesh.get_node(6001)[1]==[10.,0,0]
        @test _SCP.mesh.get_node(6002)[1]==[11.,0,0]
        @test _SCP.mesh.get_node(6012)[1]==[10.5,0,0]
        all_nodes=_SCP.mesh.get_nodes()[1]
        _,all_elements,all_connections=_SCP.mesh.get_elements()
        @test length(unique(all_nodes))==length(all_nodes)
        @test length(unique(vcat(all_elements...)))==sum(length,all_elements)
        known_nodes=Set(all_nodes)
        @test all(node->node in known_nodes,Iterators.flatten(all_connections))
        @test sum(_SCP.mesh.get_element_qualities(vcat(surface_elements...),"volume"))≈1 atol=2e-11
    finally
        _SCP.finalize()
    end
end

@testset "Actual mixed P2 simplex quality retains mapped area" begin
    # X(u,v)=(u,v-4d*u*(1-u-v),0), detJ=1+4d*u > 0.
    # Its exact area is 1/2+2d/3; this differs from the corner area for d>0.
    for displacement in (0.0,3/16)
        try
            _SCP.initialize();_SCP.option("Mesh.Renumber",0)
            _SCP.option("Mesh.MeshOnlyEmpty",1)
            _SCP.model.add_discrete_entity(2,91)
            coordinates=[0.0 1 0 .5 .5 0;0.0 0 1 -displacement .5 .5;zeros(1,6)]
            _SCP.mesh.add_nodes(2,91,1:6,vec(coordinates))
            _SCP.mesh.add_elements_by_type(91,9,[101],collect(1:6))
            _SCP.mesh.generate(2)
            @test _SCP.mesh.get_element(101)==(9,UInt64.(1:6),2,91)
            expected=.5+2displacement/3
            @test only(_SCP.mesh.get_element_qualities([101],"volume"))≈expected atol=2e-14
            # The retained overlay evaluator remains a separate control on
            # exactly these interpolation nodes, without publishing a copy.
            overlay=Tessella.HighOrder.P2TriMesh(coordinates,reshape(Int32.(1:6),6,1))
            @test _SCP._legacy_p2_quality(overlay,1,"volume","P2 area control")≈expected atol=2e-14
        finally
            _SCP.finalize()
        end
    end
end
