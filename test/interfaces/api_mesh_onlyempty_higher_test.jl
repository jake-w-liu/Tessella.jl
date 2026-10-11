using Test, Tessella
using Tessella.Elements: MixedMesh, ElementBlock
const _OE=Tessella.API

function _oe_square(;nx=7,ny=1)
    _OE.initialize();_OE.option("Mesh.Renumber",0)
    for (x,y) in ((0,0),(1,0),(1,1),(0,1));_OE.model.add_point(x,y,0);end
    for (a,b) in ((1,2),(2,3),(3,4),(4,1));_OE.model.add_line(a,b);end
    _OE.model.add_curve_loop([1,2,3,4]);_OE.model.add_plane_surface([1])
    for curve in (1,3);_OE.mesh.set_transfinite_curve(curve,nx+1);end
    for curve in (2,4);_OE.mesh.set_transfinite_curve(curve,ny+1);end
    _OE.mesh.set_transfinite_surface(1);_OE.mesh.set_recombine(2,1)
    return nothing
end

function _oe_area(nodes,coordinates,connections,width)
    points=Dict(tag=>coordinates[3i-2:3i] for (i,tag) in enumerate(nodes))
    total=0.0
    for cell in eachcol(reshape(connections,width,:))
        area=0.0
        for position in 1:width
            a=points[cell[position]];b=points[cell[mod1(position+1,width)]]
            area+=a[1]*b[2]-a[2]*b[1]
        end
        @test area>0
        total+=area/2
    end
    return total
end

@testset "Native 2D OnlyEmpty retains actual source vectors" begin
    for (nx,ny) in ((7,1),(2,2))
        try
            _oe_square(;nx,ny)
            _OE.mesh.generate(2)
            nodes=_OE.mesh.get_nodes()[1]
            labels=UInt64.(10003 .+ 7 .* reverse(0:length(nodes)-1))
            _OE.mesh.renumber_nodes(nodes,labels)
            old_cells,connections=_OE.mesh.get_elements_by_type(3,1)
            _OE.mesh.renumber_elements(old_cells,UInt64.(20003 .+ 11 .* (0:length(old_cells)-1)))
            _OE.mesh.reorder_elements(3,1,reverse(0:length(old_cells)-1))
            old_cells,connections=_OE.mesh.get_elements_by_type(3,1)
            before=_OE.mesh.get_nodes()
            source_coords=copy(_OE.LAST_MESH[].coords)
            records=Dict(tag=>_OE.mesh.get_node(tag) for tag in labels)
            highwater=_OE.mesh.get_max_element_tag()
            _OE.option("Mesh.MeshOnlyEmpty",1)
            actual=_OE.mesh.generate(2)
            @test actual isa MixedMesh
            cells,after=_OE.mesh.get_elements_by_type(3,1)
            @test after==connections
            @test length(cells)==nx*ny
            @test minimum(cells)>highwater
            @test isempty(intersect(cells,old_cells))
            @test Set(_OE.mesh.get_nodes()[1])==Set(labels)
            @test Dict(tag=>_OE.mesh.get_node(tag) for tag in labels)==records
            @test _oe_area(before[1],before[2],after,4)≈1 atol=8eps(Float64)
            authority=_OE.CURRENT[].meshing.mesh_authority
            source=Tessella.Model._model_mesh_product(authority,2,1)
            @test source.coords==source_coords
            @test source.entity_data.external_element_tags==[cells]
            @test source.entity_data.external_node_tags[vec(only(source.blocks).nodes)]==after
            # Mesh2D resets faces to PENDING before the OnlyEmpty return.
            @test Tessella.Model._model_mesh_status(_OE.CURRENT[],authority)==1
            # A second request consumes the published actual source, and the
            # native SetOrder1 factory allocates a fresh element table again.
            _OE.mesh.generate(2)
            second,again=_OE.mesh.get_elements_by_type(3,1)
            @test again==after
            @test minimum(second)>maximum(cells)
            @test Set(_OE.mesh.get_nodes()[1])==Set(labels)
        finally
            _OE.finalize()
        end
    end
end

@testset "OnlyEmpty phase and native status distinguish raw/discrete cells" begin
    @test Tessella.Model._model_mesh_generation_phases(2,2)==(one=false,two=true,three=false)
    @test Tessella.Model._model_mesh_generation_phases(3,2)==(one=false,two=false,three=true)
    @test Tessella.Model._model_mesh_generation_phases(3,0)==(one=true,two=true,three=true)
    try
        _oe_square(;nx=2,ny=2)
        _OE.option("Mesh.MeshOnlyEmpty",1)
        _OE.mesh.generate(2)
        @test length(_OE.mesh.get_elements_by_type(3,1)[1])==4
        @test length(_OE.mesh.get_elements_by_type(1)[1])==8
        @test _oe_area(_OE.mesh.get_nodes()[1],_OE.mesh.get_nodes()[2],
            _OE.mesh.get_elements_by_type(3,1)[2],4)≈1 atol=8eps(Float64)
        _OE.option("Mesh.ElementOrder",2)
        _OE.mesh.generate(2)
        @test length(_OE.mesh.get_elements_by_type(10,1)[1])==4
        @test all(isfinite,_OE.mesh.get_nodes()[2])
        @test Set(_OE.mesh.get_elements()[1])==Set(Int32[15,8,10])
    finally
        _OE.finalize()
    end
end
