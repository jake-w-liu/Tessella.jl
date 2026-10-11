module GeoHybridPyramidOptionsFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements
const FIXTURE=joinpath(@__DIR__,"..","fixtures","hybrid_normal_apex.geo")

snapshot()=(A.model.get_current(),A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),
            A.option("Mesh.OptimizePyramids"))

function check_volume(mesh)
    data=mesh.entity_data
    @test data!==nothing
    cells=Dict(Int(block.msh)=>block.nodes for block in mesh.blocks
               if E.msh_dimension(Int(block.msh))==3)
    @test Set(keys(cells))==Set([4,7]) && size(cells[7],2)==6
    function tet_volume(nodes)
        a,b,c,d=(mesh.coords[:,row] for row in nodes)
        det(hcat(b-a,c-a,d-a))/6
    end
    tet_volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
    pyramid_volumes=[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                     for cell in eachcol(cells[7])]
    @test all(>(0),tet_volumes) && all(>(0),pyramid_volumes)
    @test sum(tet_volumes)+sum(pyramid_volumes)≈1.0 atol=2e-10
    offset=0.1sqrt(2.0)
    expected=([0.5,0.5,offset],[0.5,0.5,1-offset],
              [0.5,offset,0.5],[1-offset,0.5,0.5],
              [0.5,1-offset,0.5],[offset,0.5,0.5])
    apices=vec(cells[7][5,:])
    @test all(count(row->isapprox(mesh.coords[:,row],point;atol=8eps(),rtol=0),apices)==1
              for point in expected)
    @test all(isapprox(volume,offset/3;atol=8eps(),rtol=0) for volume in pyramid_volumes)
    @test all(data.node_entities[row]==(3,Int32(1)) for row in apices)
    @test all(all(==(Int32(1)),data.block_entities[index])
              for (index,block) in enumerate(mesh.blocks) if E.msh_dimension(Int(block.msh))==3)
    incidence=Dict{Tuple,Int}()
    for (kind,nodes) in cells
        faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
            ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
        for cell in eachcol(nodes),face in faces
            key=Tuple(sort(cell[collect(face)]))
            incidence[key]=get(incidence,key,0)+1
        end
    end
    quads=Set(Tuple(sort(block.nodes[:,column])) for block in mesh.blocks
              if block.msh==3 for column in axes(block.nodes,2))
    @test Set(face for (face,count) in incidence if count==1)==quads && length(quads)==6
    @test all(count->count in (1,2),values(incidence))
    @test Tessella.validate(mesh).ok
    return nothing
end

function cell_records(mesh)
    data=mesh.entity_data
    Dict(uid=>(Int(block.msh),data.block_entities[index][column],
               Tuple(data.external_node_tags[block.nodes[:,column]]))
         for (index,block) in enumerate(mesh.blocks)
         for (column,uid) in enumerate(data.external_element_tags[index]))
end

@testset "GEO and CLI transport their own floating-point Pyramid setting" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add_point(0,0,0;tag=101,meshSize=1.0)
        A.model.add_point(1,0,0;tag=102,meshSize=1.0)
        A.model.add_line(101,102;tag=101)
        A.mesh.generate(1)
        A.option("Mesh.OptimizePyramids",2.0)
        before=snapshot()
        execution=Tessella.execute_geo(FIXTURE;mesh_dim=3)
        @test execution.values["stored_factor"]==-0.1
        volume=Tessella.geo_entity_mesh(execution,3,1)
        projected=Tessella.Model.model_to_mixed(execution.model,volume,3,1)
        check_volume(projected)
        @test snapshot()==before
        @test projected.physical_names==Dict((2,7)=>"boundary",(3,3)=>"solid")
        mktempdir() do directory
            source=joinpath(directory,"cube.geo")
            output=joinpath(directory,"cube.msh")
            original=read(FIXTURE,String)
            write(source,original)
            @test Tessella.CLI.main([source,"-3","-o",output])==output
            recovered=E.read_mixed_msh(output)
            check_volume(recovered)
            @test Set(recovered.entity_data.external_node_tags)==Set(projected.entity_data.external_node_tags)
            @test cell_records(recovered)==cell_records(projected)
            @test recovered.physical_names==projected.physical_names
            @test recovered.entity_data.entities[(3,1)].boundaries==Int32[1,2,3,4,5,6]
            @test snapshot()==before && read(source,String)==original
            saved=read(output)
            write(source,replace(original,"Mesh.OptimizePyramids=-0.1;"=>"Mesh.OptimizePyramids=-1;"))
            @test_throws ArgumentError Tessella.execute_geo(source;mesh_dim=3)
            @test snapshot()==before
            @test_throws ArgumentError Tessella.CLI.main([source,"-3","-o",output])
            @test read(output)==saved && snapshot()==before
        end
    finally
        A.finalize()
    end
end
end
