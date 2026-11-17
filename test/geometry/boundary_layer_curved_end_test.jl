module BoundaryLayerCurvedEndChecks
using Test,Tessella
const A=Tessella.API
const B=Tessella.BoundaryLayer
const M=Tessella.Model
const fixture=joinpath(@__DIR__,"..","fixtures","geo_boundary_layer_curved_end.geo")

@testset "Boundary-layer columns retain curved adjacent-edge nodes" begin
    wall=Tessella.Mesh([0. 1.;0. 0.;0. 0.];segs=reshape(Int32[1,2],2,1))
    positions=Dict((1,1)=>(-.05,.1,0.),(1,2)=>(-.2,.2,0.),
                   (2,1)=>(1.05,.1,0.),(2,2)=>(1.2,.2,0.))
    snapshot=deepcopy(positions)
    strip=B._bl_2d_offsets(wall,[.1,.2];vertex_positions=positions)
    @test Tessella.Elements.validate(strip.mesh).ok
    for ((v,k),point) in positions
        @test Tuple(strip.mesh.coords[:,strip.id_of[(v,k,0)]])==point
    end
    @test positions==snapshot
    @test strip.tops[(0.,0.,0.)]==positions[(1,2)]
    @test strip.tops[(1.,0.,0.)]==positions[(2,2)]
    @test_throws ArgumentError B._bl_2d_offsets(wall,[.1];
        vertex_positions=Dict((3,1)=>(0.,.1,0.)))
    @test_throws ArgumentError B._bl_2d_offsets(wall,[.1];
        vertex_positions=Dict((1,2)=>(0.,.2,0.)))
    @test_throws ArgumentError B._bl_2d_offsets(wall,[.1];
        vertex_positions=Dict((1,1)=>(NaN,.1,0.)))

    # Public Gmsh 4.15.2 successfully meshes this Circle adjacent to a BL
    # wall (44 triangles,15 quads). Protect the actual curve-node support
    # and ownership without imposing another mesher's interior triangulation.
    A.initialize()
    try
        A.option("General.Terminal",0)
        @test A.option("Mesh.BoundaryLayerFanElements")==5.
        @test A.option("Mesh.BoundaryLayerFanElements",-2.5)==0.
        A.option("Mesh.BoundaryLayerFanElements",5)
        A.open_geo!(fixture)
        A.mesh.generate(1)
        tags,xyz,_=A.mesh.get_nodes(1,2,true)
        @test length(tags)>3
        before=Dict(tag=>Tuple(xyz[3i-2:3i]) for (i,tag) in enumerate(tags))
        A.mesh.generate(2)
        after_tags,after_xyz,_=A.mesh.get_nodes(1,2,true)
        after=Dict(tag=>Tuple(after_xyz[3i-2:3i]) for (i,tag) in enumerate(after_tags))
        @test after==before
        kinds,_,connectivity=A.mesh.get_elements(2,1)
        @test 3 in kinds
        @test length(only(connectivity[findall(==(3),kinds)]))>0
        used=Set(Iterators.flatten(connectivity))
        @test all(in(used),tags)
        surface_tags,surface_xyz,_=A.mesh.get_nodes(2,1)
        @test isempty(intersect(tags,surface_tags))
        # No duplicated Surface-owned point may sit on the curved end column.
        @test all(point->all(i->Tuple(surface_xyz[3i-2:3i])!=point,
            eachindex(surface_tags)),values(before))
        @test all(point->isapprox((point[1]+1)^2+point[2]^2,4.;atol=2e-14),
            values(after))
        saved_nodes=deepcopy(A.mesh.get_nodes())
        saved_elements=deepcopy(A.mesh.get_elements())
        A.mesh.field.set_number(1,"Thickness",-1)
        @test_throws ArgumentError A.mesh.generate(2)
        @test A.mesh.get_nodes()==saved_nodes
        @test A.mesh.get_elements()==saved_elements
        A.mesh.field.remove(1)
        A.mesh.clear()
        A.mesh.generate(2)
        @test !(3 in A.mesh.get_elements(2,1)[1])
    finally
        A.finalize()
    end
    for transform in ("", "Rotate {{1,0,0},{0,0,0},Pi/3} { Surface{1}; }")
        mktemp() do path,io
            write(io,read(fixture,String)*transform);close(io)
            result=execute_geo(path;mesh_dim=2)
            @test isempty(result.errors)
            @test Tessella.Elements.validate(result.mesh).ok
        end
    end
end

@testset "Empty endpoint layers retain native curve births" begin
    model=Tessella.GeoModel()
    Tessella.add_point!(model,0.,0.,0.;tag=1,mesh_size=.25)
    Tessella.add_point!(model,1.,0.,0.;tag=2,mesh_size=.25)
    Tessella.add_line!(model,1,2;tag=1)
    caller="empty endpoint layer regression"
    options=M._model_default_mesh1d_options(model,caller)
    ordinary=M._model_mesh_curve!(model,1,options,caller;_native_producer=true)
    spec=Tessella.IO.GeoFieldSpec(1,"BoundaryLayer",Dict(
        "PointsList"=>"{1,2}","Size"=>"0.5","Thickness"=>"1"))
    layer=M._ModelBoundaryLayerField(1,
        Tessella.SizeField._geo_boundary_layer_config(spec,5),caller)
    for (first,last) in ((nothing,nothing),(layer,nothing),(nothing,layer),(layer,layer))
        graded=M._model_curve_bl_grade_params(model,1,first,last,0.,1.,1.,
            options,caller;_native_producer=true)
        @test graded isa M._ModelNativeCurveGrade
        @test (graded.parameters,graded.rows,graded.count)==
            (ordinary.parameters,ordinary.rows,ordinary.count)
    end

    square="""
    Point(1)={0,0,0,0.25}; Point(2)={1,0,0,0.25};
    Point(3)={1,1,0,0.25}; Point(4)={0,1,0,0.25};
    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
    """
    plain,layered=map(("","""
    Field[1]=BoundaryLayer; Field[1].PointsList={1};
    Field[1].Size=0.5; Field[1].Thickness=1; BoundaryLayer Field=1;
    """)) do fields
        mktemp() do path,io
            write(io,square*fields);close(io)
            execute_geo(path;mesh_dim=2)
        end
    end
    @test isempty(plain.errors) && isempty(layered.errors)
    mesh=layered.mesh
    @test Tessella.Elements.validate(mesh).ok
    @test isequal(mesh.coords,plain.mesh.coords)
    @test mesh.segs==plain.mesh.segs && mesh.tris==plain.mesh.tris
    areas=[((mesh.coords[1,b]-mesh.coords[1,a])*(mesh.coords[2,c]-mesh.coords[2,a])-
            (mesh.coords[2,b]-mesh.coords[2,a])*(mesh.coords[1,c]-mesh.coords[1,a]))/2
           for (a,b,c) in eachcol(mesh.tris)]
    @test all(>(0),areas)
    @test isapprox(sum(areas),1.;rtol=64eps(Float64))
    actual=layered.model.meshing.mesh_authority.mesh
    expected=plain.model.meshing.mesh_authority.mesh
    @test Tessella.Elements.validate(actual).ok
    @test isequal(actual.coords,expected.coords)
    @test [(b.msh,b.nodes,b.tags) for b in actual.blocks]==
        [(b.msh,b.nodes,b.tags) for b in expected.blocks]
    data=actual.entity_data;expected_data=expected.entity_data
    @test isequal((data.external_node_tags,data.node_entities,data.node_parametric,
                   data.external_element_tags,data.block_entities),
                  (expected_data.external_node_tags,expected_data.node_entities,
                   expected_data.node_parametric,expected_data.external_element_tags,
                   expected_data.block_entities))
end
end
