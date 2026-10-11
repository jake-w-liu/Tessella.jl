module ModelImportGenerationChecks
using Test,Tessella
const A=Tessella.API

@testset "atomic file import with a tiny actual GEO mesh" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"generated-line.geo")
            Base.write(path,"Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(3)={1,2};Transfinite Curve{3}=3;Mesh 1;\n")
            @test A.open(path)===nothing
            node_uids,coordinates,_=A.mesh.get_nodes()
            xyz=reshape(coordinates,3,:)
            @test Set(node_uids)==Set(UInt64[1,2,3]) && sort(vec(xyz[1,:]))==[0.,.5,1.] && all(iszero,xyz[2:3,:])
            @test Set(A.model.get_entities())==Set([(0,1),(0,2),(1,3)])
            types,element_uids,connectivity=A.mesh.get_elements(1,3)
            @test types==[1] && element_uids==[UInt64[3,4]] && connectivity==[UInt64[1,3,3,2]]
            types,element_uids,_=A.mesh.get_elements()
            @test Set(types)==Set([1,15]) && Set(reduce(vcat,element_uids))==Set(UInt64[1,2,3,4])
            mid_uids,mid_coordinates,parameters=A.mesh.get_nodes(1,3)
            @test mid_uids==[3] && mid_coordinates==[.5,0,0] && parameters==[.5]

            A.model.add("legacy generated");count=length(A.model.list())
            result=A.open_geo!(path;mesh_dim=1)
            @test length(A.model.list())==count && A.model.get_current()=="legacy generated" && A.model.get_file_name()==path
            node_uids,coordinates,_=A.mesh.get_nodes();xyz=reshape(coordinates,3,:)
            println("Legacy generated nodes: ",(uids=node_uids,xyz=xyz,result_xyz=result.mesh.coords))
            @test length(node_uids)==3 && isapprox(sort(vec(xyz[1,:])),[0.,.5,1.];rtol=0,atol=1e-11) && all(iszero,xyz[2:3,:])
            types,tags,_=A.mesh.get_elements(1,3)
            @test types==[1] && length(only(tags))==2 && length(A.mesh.get_nodes(1,3)[1])==1
            @test result.mesh!==nothing && size(result.mesh.coords,2)==3
            before=A.mesh.get_nodes();result.mesh.coords[:,1].=9.;result.model.points[1]=(8.,8.,8.)
            @test A.mesh.get_nodes()==before && A.model.get_value(0,1,Float64[])==[0.,0,0]
        end
    finally;A.finalize();end
end
end
