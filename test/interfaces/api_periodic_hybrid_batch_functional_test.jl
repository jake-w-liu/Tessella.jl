module PeriodicHybridBatchFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements
const M=Tessella.Model
const TRANSLATION=[1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
function cube(region,offset,shift;copy_slave=true,affine=TRANSLATION,cross_region=false)
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points)
        A.model.add_point(p[1]+shift,p[2],p[3];tag=tag+offset,meshSize=1.)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        A.model.add_line((point+offset for point in ends)...;tag=tag+offset)
        A.mesh.set_transfinite_curve(tag+offset,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop([sign(line)*(abs(line)+offset) for line in lines];tag=tag+offset)
        A.model.add_plane_surface([tag+offset];tag=tag+offset)
        (tag!=(cross_region ? 1 : 2) || !copy_slave) && A.mesh.set_transfinite_surface(tag+offset)
        A.mesh.set_recombine(2,tag+offset)
    end
    A.model.add_surface_loop(collect(1:6).+offset;tag=region)
    A.model.add_volume([region];tag=region)
    A.mesh.set_periodic(2,[cross_region ? 1+offset : 2+offset],[cross_region ? 1 : 1+offset],affine)
end

function two_cubes(;copy_second=true,second_affine=TRANSLATION,cross_region=false)
    cube(1,0,0.)
    cube(2,100,3.;copy_slave=copy_second,affine=second_affine,cross_region)
end
snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
function cells(dim,tag)
    kinds,labels,vertices=A.mesh.get_elements(dim,tag)
    return Dict(uid=>(Int(kind),Tuple(matrix[:,column]))
        for (kind,uids,nodes) in zip(kinds,labels,vertices)
        for matrix in (reshape(nodes,E.msh_num_nodes(kind),length(uids)),)
        for (column,uid) in enumerate(uids))
end
function correspondence(slave,master,high,expected)
    link=A.mesh.get_periodic_nodes(2,slave,high)
    @test link.master_entity==master && link.affine==TRANSLATION
    @test length(link.slave_nodes)==length(link.master_nodes)==expected
    @test allunique(link.slave_nodes) && allunique(link.master_nodes)
    @test all(first(A.mesh.get_node(s))≈first(A.mesh.get_node(t))+[0.,0.,1.]
        for (s,t) in zip(link.slave_nodes,link.master_nodes))
    return Dict(zip(link.master_nodes,link.slave_nodes))
end

@testset "two separated periodic hybrid Regions consume actual retained copy events" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.OptimizePyramids",0.)
        two_cubes()
        A.mesh.generate(2)
        lower=Dict(owner=>cells(owner...) for owner in A.model.get_entities() if first(owner)<3)
        nodes_before=A.mesh.get_nodes()[1:2]
        lower_tags,lower_xyz=nodes_before
        lower_coordinates=Dict(uid=>lower_xyz[3row-2:3row] for (row,uid) in enumerate(lower_tags))
        original=Dict(slave=>correspondence(slave,master,false,4) for (slave,master) in ((2,1),(102,101)))
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        @test all(cells(owner...)==records for (owner,records) in lower)
        tags,xyz=nodes_before
        @test all(first(A.mesh.get_node(uid))==xyz[3row-2:3row] for (row,uid) in enumerate(tags))
        @test correspondence(2,1,false,4)==original[2]
        @test correspondence(102,101,false,4)==original[102]
        ids,xyz,_=A.mesh.get_nodes()
        coordinates=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        tet_volume(nodes)=begin
            a,b,c,d=(coordinates[uid] for uid in nodes)
            det(hcat(b-a,c-a,d-a))/6
        end
        regions=Dict{Int,Dict}()
        volumes=Float64[]
        for (region,offset) in ((1,0),(2,100))
            kinds,labels,vertices=A.mesh.get_elements(3,region)
            @test Set(kinds)==Set(Int32[4,7])
            blocks=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(uids))
                for (kind,uids,nodes) in zip(kinds,labels,vertices))
            @test size(blocks[4],2)>0 && size(blocks[7],2)==6
            contributions=[tet_volume(cell) for cell in eachcol(blocks[4])]
            append!(contributions,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                for cell in eachcol(blocks[7])])
            @test all(>(0),contributions)
            @test sum(contributions)≈1.0 atol=2e-10
            push!(volumes,sum(contributions))
            incidence=Dict{Tuple,Int}()
            for (kind,block) in blocks
                roles=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
                    ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
                for cell in eachcol(block),face in roles
                    key=Tuple(sort(cell[collect(face)]));incidence[key]=get(incidence,key,0)+1
                end
            end
            exterior=Set(Tuple(sort(only(last(A.mesh.get_elements(2,surface+offset))))) for surface in 1:6)
            @test all(count->count in (1,2),values(incidence))
            @test Set(face for (face,count) in incidence if count==1)==exterior
            @test Set(Tuple(sort(cell[1:4])) for cell in eachcol(blocks[7]))==exterior
            regions[region]=cells(3,region)
        end
        @test isempty(intersect(Set(keys(regions[1])),Set(keys(regions[2]))))
        @test Tessella.validate(A.mesh.get()).ok
        println("PERIODIC_HYBRID_BATCH Regions=2 volumes=",volumes);flush(stdout)
        A.mesh.set_order(2)
        @test all(Set(first(A.mesh.get_elements(3,region)))==Set(Int32[11,14]) for region in (1,2))
        # Default native renumbering gives old volume UIDs new roles at P2;
        # the retained lower corner identities and coordinates stay fixed.
        @test all(first(A.mesh.get_node(uid))==point for (uid,point) in lower_coordinates)
        expected=Dict(slave=>correspondence(slave,master,true,9) for (slave,master) in ((2,1),(102,101)))
        @test correspondence(2,1,false,4)==original[2]
        @test correspondence(102,101,false,4)==original[102]
        current=Dict((3,region)=>cells(3,region) for region in (1,2))
        mktempdir() do directory
            path=joinpath(directory,"periodic-regions.msh")
            A.write(path)
            restored=E.read_mixed_msh(path);data=restored.entity_data
            @test Tessella.validate(restored).ok
            @test Set(data.external_node_tags)==Set(first(A.mesh.get_nodes()))
            @test length(restored.periodic_links)>=2
            for (slave,master) in ((2,1),(102,101))
                link=only(link for link in restored.periodic_links if link.dim==2 && link.slave_entity==slave)
                @test link.master_entity==master && collect(link.affine)==TRANSLATION
                @test Dict(zip(data.external_node_tags[link.master_nodes],data.external_node_tags[link.slave_nodes]))==original[slave]
            end
            for region in (1,2)
                recovered=Dict(uid=>(Int(block.msh),Tuple(data.external_node_tags[block.nodes[:,column]]))
                    for (index,block) in enumerate(restored.blocks) if E.msh_dimension(block.msh)==3
                    for (column,uid) in enumerate(data.external_element_tags[index])
                    if data.block_entities[index][column]==region)
                @test recovered==current[(3,region)]
            end
        end
    finally
        A.finalize()
    end
end

@testset "unsupported periodic batch copy contexts reject atomically before construction" begin
    shear=[1.,0.,.25,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
    across=[1.,0.,0.,3.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
    for settings in ((copy_second=false,second_affine=TRANSLATION,cross_region=false),
                     (copy_second=true,second_affine=shear,cross_region=false),
                     (copy_second=true,second_affine=across,cross_region=true))
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.OptimizePyramids",0.)
            two_cubes(;settings...)
            # The cross-Region relation's derived Surface copies are valid
            # dim-2 products like any other periodic slave; the Region
            # certificate rejects the foreign pair at volume time.
            A.mesh.generate(2)
            A.option("Mesh.MeshOnlyEmpty",1)
            before=snapshot()
            @test_throws ArgumentError A.mesh.generate(3)
            @test snapshot()==before
        finally
            A.finalize()
        end
    end
end

@testset "incomplete actual copy event inventory cannot certify an otherwise valid lower mesh" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        two_cubes()
        A.mesh.generate(2)
        model=A.CURRENT[];authority=model.meshing.mesh_authority
        groups,receipts=M._model_native_hybrid_groups(model,[1,2],authority,"periodic inventory regression")
        events=Dict(slave=>(Int(history.master),[(pair.slave,pair.master) for pair in history.primary])
            for ((dim,slave),history) in model.meshing.periodic_uid_history if dim==2)
        @test Set(keys(events))==Set((2,102)) && all(length(last(record))==4 for record in values(events))
        incomplete=copy(events);master,pairs=incomplete[102];incomplete[102]=(master,pairs[1:3])
        before=snapshot()
        @test_throws ArgumentError M._model_native_hybrid_batch_periodic_certificate(model,receipts[2],incomplete,0.,"periodic inventory regression")
        @test snapshot()==before
    finally
        A.finalize()
    end
end
end
