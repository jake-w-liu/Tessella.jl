module HybridSharedPeriodicFunctionalChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const E=Tessella.Elements
const Z_TRANSLATION=[1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]
const Y_TRANSLATION=[1.,0.,0.,0.,0.,1.,0.,1.,0.,0.,1.,0.,0.,0.,0.,1.]

function adjacent_cubes(;conflicting=false)
    corners=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(corners);A.model.add_point(point...;tag,meshSize=1.);end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,pair) in enumerate(edges);A.model.add_line(pair...;tag);A.mesh.set_transfinite_curve(tag,2);end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,loop) in enumerate(loops)
        A.model.add_curve_loop(tag==1 ? collect(loop) : [-line for line in reverse(loop)];tag)
        A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1);A.model.add_volume([1];tag=1)
    for (tag,point) in zip(101:104,((2,0,0),(2,1,0),(2,0,1),(2,1,1)))
        A.model.add_point(point...;tag,meshSize=1.)
    end
    more=((2,101),(101,102),(102,3),(6,103),(103,104),(104,7),(101,103),(102,104))
    for (tag,pair) in zip(101:108,more);A.model.add_line(pair...;tag);A.mesh.set_transfinite_curve(tag,2);end
    right=((101,102,103,-2),(104,105,106,-6),(101,107,-104,-10),(102,108,-105,-107),(103,11,-106,-108))
    for (index,loop) in enumerate(right)
        tag=100+index
        A.model.add_curve_loop(index==1 ? collect(loop) : [-line for line in reverse(loop)];tag)
        A.model.add_plane_surface([tag];tag)
    end
    A.model.add_surface_loop([101,102,103,104,105,-4];tag=101);A.model.add_volume([101];tag=2)
    slaves=conflicting ? (2,105) : (2,102)
    for surface in vcat(collect(1:6),collect(101:105))
        surface in slaves || A.mesh.set_transfinite_surface(surface)
        surface==4 || A.mesh.set_recombine(2,surface)
    end
    A.mesh.set_periodic(2,[2],[1],Z_TRANSLATION)
    A.mesh.set_periodic(2,[conflicting ? 105 : 102],[conflicting ? 103 : 101],
        conflicting ? Y_TRANSLATION : Z_TRANSLATION)
    A.model.add_physical_group(3,[1];tag=201,name="left")
    A.model.add_physical_group(3,[2];tag=202,name="right")
    A.model.add_physical_group(2,[4];tag=301,name="interface")
end
function records(dim,tag)
    kinds,tags,flat=A.mesh.get_elements(dim,tag)
    Dict(uid=>(Int(kind),Tuple(matrix[:,column])) for (kind,labels,nodes) in zip(kinds,tags,flat)
        for matrix in (reshape(nodes,E.msh_num_nodes(kind),length(labels)),)
        for (column,uid) in enumerate(labels))
end
snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),
    A.mesh.get_periodic_nodes(2,2),A.mesh.get_periodic_nodes(2,102))
function session(body)
    A.initialize()
    try
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0);A.option("Mesh.SaveAll",1)
        A.option("Mesh.OptimizePyramids",0.)
        body()
    finally;A.finalize();end
end
function correspondence(slave,master,high,count)
    link=A.mesh.get_periodic_nodes(2,slave,high)
    @test link.master_entity==master && link.affine==Z_TRANSLATION
    @test length(link.slave_nodes)==length(link.master_nodes)==count &&
        allunique(link.slave_nodes) && allunique(link.master_nodes)
    @test all(first(A.mesh.get_node(s))≈first(A.mesh.get_node(m))+[0.,0.,1.]
        for (s,m) in zip(link.slave_nodes,link.master_nodes))
    Dict(zip(link.master_nodes,link.slave_nodes))
end

@testset "connected periodic hybrid Regions retain their shared Tri interface" begin
    session() do
        adjacent_cubes();A.mesh.generate(2)
        lower=Dict(owner=>records(owner...) for owner in A.model.get_entities() if owner[1]<3)
        shared=Tuple[last(cell) for cell in values(lower[(2,4)])]
        @test length(shared)==2 && all(first(cell)==2 for cell in values(lower[(2,4)]))
        lower_ids,lower_xyz,_=A.mesh.get_nodes()
        original=Dict(slave=>correspondence(slave,master,false,4) for (slave,master) in ((2,1),(102,101)))
        A.option("Mesh.MeshOnlyEmpty",1);A.mesh.generate(3)
        @test all(Set(values(records(owner...)))==Set(values(cells)) for (owner,cells) in lower)
        @test all(first(A.mesh.get_node(uid))==lower_xyz[3row-2:3row] for (row,uid) in enumerate(lower_ids))
        @test correspondence(2,1,false,4)==original[2]
        @test correspondence(102,101,false,4)==original[102]
        ids,xyz,_=A.mesh.get_nodes()
        positions=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        tet_volume(cell)=begin;a,b,c,d=(positions[uid] for uid in cell);det(hcat(b-a,c-a,d-a))/6;end
        incidences=Dict{Int,Dict}();regions=Dict{Int,Dict}();tets=Dict{Int,Matrix{UInt64}}()
        for region in (1,2)
            kinds,tags,flat=A.mesh.get_elements(3,region)
            @test Set(kinds)==Set(Int32[4,7])
            blocks=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(labels)) for (kind,labels,nodes) in zip(kinds,tags,flat))
            @test size(blocks[4],2)>0 && size(blocks[7],2)==5
            volumes=[tet_volume(cell) for cell in eachcol(blocks[4])]
            append!(volumes,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]]) for cell in eachcol(blocks[7])])
            @test all(>(0),volumes) && isapprox(sum(volumes),1.;atol=2e-10)
            incidence=Dict{Tuple,Int}()
            for (kind,block) in blocks
                roles=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) : ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
                for cell in eachcol(block),slots in roles
                    face=Tuple(sort(cell[collect(slots)]));incidence[face]=get(incidence,face,0)+1
                end
            end
            surfaces=region==1 ? collect(1:6) : [101,102,103,104,105,4]
            boundary=Set(Tuple(sort(collect(nodes))) for surface in surfaces for (_,nodes) in values(lower[(2,surface)]))
            @test all(count->count in (1,2),values(incidence)) && Set(face for (face,count) in incidence if count==1)==boundary
            @test all(incidence[Tuple(sort(collect(face)))]==1 for face in shared)
            @test A.model.get_physical_groups_for_entity(3,region)==[200+region]
            incidences[region]=incidence;regions[region]=records(3,region);tets[region]=blocks[4]
            println("SHARED_PERIODIC_HYBRID Region=",region," Tet4=",size(blocks[4],2)," Pyramid5=5 volume=",sum(volumes));flush(stdout)
        end
        @test isempty(intersect(Set(keys(regions[1])),Set(keys(regions[2]))))
        for face in shared
            opposite=UInt64[only(node for node in cell if !(node in face)) for region in (1,2)
                for cell in eachcol(tets[region]) if all(node->node in cell,face)]
            @test length(opposite)==2
            a,b,c=(positions[node] for node in face)
            sides=[dot(cross(b-a,c-a),positions[node]-a) for node in opposite]
            @test sides[1]*sides[2]<0
        end
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.set_order(2)
        @test all(Set(first(A.mesh.get_elements(3,region)))==Set(Int32[11,14]) for region in (1,2))
        @test all(first(A.mesh.get_node(uid))==lower_xyz[3row-2:3row] for (row,uid) in enumerate(lower_ids))
        for (slave,master) in ((2,1),(102,101))
            correspondence(slave,master,true,9)
            @test correspondence(slave,master,false,4)==original[slave]
        end
        current=Dict(owner=>records(owner...) for owner in A.model.get_entities())
        mktempdir() do directory
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",0)
            path=joinpath(directory,"shared-periodic.msh");A.write(path)
            restored=E.read_mixed_msh(path);data=restored.entity_data
            @test Tessella.validate(restored).ok
            @test Set(data.external_node_tags)==Set(first(A.mesh.get_nodes()))
            for (slave,master) in ((2,1),(102,101))
                link=only(link for link in restored.periodic_links if link.dim==2 && link.slave_entity==slave)
                @test link.master_entity==master && collect(link.affine)==Z_TRANSLATION
                @test Dict(zip(data.external_node_tags[link.master_nodes],data.external_node_tags[link.slave_nodes]))==original[slave]
            end
            recovered=Dict{Tuple{Int,Int},Dict}()
            for (index,block) in enumerate(restored.blocks),column in axes(block.nodes,2)
                owner=(E.msh_dimension(block.msh),Int(data.block_entities[index][column]))
                get!(Dict,recovered,owner)[data.external_element_tags[index][column]]=
                    (Int(block.msh),Tuple(data.external_node_tags[block.nodes[:,column]]))
            end
            @test recovered==current
            @test A.model.get_physical_groups_for_entity(2,4)==[301]
        end
    end
end
@testset "unsupported periodic Pyramid mode rejects before publication" begin
    session() do
        adjacent_cubes();A.mesh.generate(2);A.option("Mesh.MeshOnlyEmpty",1);A.option("Mesh.OptimizePyramids",-.1)
        before=snapshot()
        error=try;A.mesh.generate(3);nothing;catch failure;failure isa InterruptException && rethrow();failure;end
        @test error isa ArgumentError && occursin("require default Pyramid optimization",sprint(showerror,error))
        @test snapshot()==before
    end
end
@testset "unsupported mixed-direction joint copies reject atomically" begin
    session() do
        adjacent_cubes(;conflicting=true)
        # The mixed-direction slave copies are themselves valid dim-2
        # products; the shared interface slave UID is bound by both Regions
        # under different masters/affines, so the joint periodic certificate
        # rejects the conflicting copy at volume time.
        A.mesh.generate(2)
        before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
        error=try;A.mesh.generate(3);nothing;catch failure;failure isa InterruptException && rethrow();failure;end
        @test error isa ArgumentError && occursin("conflicting current master UID or affine",sprint(showerror,error))
        @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==before
    end
end
end
