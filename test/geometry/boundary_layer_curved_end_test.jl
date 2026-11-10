module BoundaryLayerCurvedEndChecks
using Test,Tessella
const A=Tessella.API
const B=Tessella.BoundaryLayer
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
end
