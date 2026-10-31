module PeriodicHybridFunctionalChecks
using Test,LinearAlgebra,Tessella
const A=Tessella.API
const E=Tessella.Elements
const TRANSLATION=[1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.]

function cube(;copy_slave=true,affine=TRANSLATION)
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points)
        A.model.add_point(p...;tag,meshSize=1.)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        A.model.add_line(ends...;tag)
        A.mesh.set_transfinite_curve(tag,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        A.model.add_curve_loop(collect(lines);tag)
        A.model.add_plane_surface([tag];tag)
        (tag!=2 || !copy_slave) && A.mesh.set_transfinite_surface(tag)
        A.mesh.set_recombine(2,tag)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    A.mesh.set_periodic(2,[2],[1],affine)
end

snapshot()=(A.mesh.get_nodes(),A.mesh.get_elements(),
            A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())

function correspondence(high_order,expected)
    link=A.mesh.get_periodic_nodes(2,2,high_order)
    @test link.master_entity==1 && link.affine==TRANSLATION
    @test length(link.slave_nodes)==length(link.master_nodes)==expected
    @test allunique(link.slave_nodes) && allunique(link.master_nodes)
    @test all(first(A.mesh.get_node(slave))≈first(A.mesh.get_node(master))+[0.,0.,1.]
              for (slave,master) in zip(link.slave_nodes,link.master_nodes))
    return Dict(zip(link.master_nodes,link.slave_nodes))
end

@testset "translated copied Quad boundary admits positive hybrid cube" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube()
        A.mesh.generate(3)
        families,labels,vertices=A.mesh.get_elements(3,1)
        @test Set(families)==Set(Int32[4,7])
        cells=Dict(Int(kind)=>reshape(nodes,E.msh_num_nodes(kind),length(ids))
                   for (kind,ids,nodes) in zip(families,labels,vertices))
        @test size(cells[4],2)>0 && size(cells[7],2)==6
        ids,xyz,_=A.mesh.get_nodes()
        points=Dict(uid=>xyz[3row-2:3row] for (row,uid) in enumerate(ids))
        function tet_volume(nodes)
            a,b,c,d=(points[uid] for uid in nodes)
            det(hcat(b-a,c-a,d-a))/6
        end
        volumes=[tet_volume(cell) for cell in eachcol(cells[4])]
        append!(volumes,[tet_volume(cell[[1,2,3,5]])+tet_volume(cell[[1,3,4,5]])
                         for cell in eachcol(cells[7])])
        @test all(>(0),volumes)
        @test sum(volumes)≈1.0 atol=2e-10
        incidence=Dict{Tuple,Int}()
        for (kind,nodes) in cells
            faces=kind==4 ? ((1,2,3),(1,2,4),(1,3,4),(2,3,4)) :
                ((1,2,3,4),(1,2,5),(2,3,5),(3,4,5),(4,1,5))
            for cell in eachcol(nodes),face in faces
                key=Tuple(sort(cell[collect(face)]))
                incidence[key]=get(incidence,key,0)+1
            end
        end
        quads=Tuple[]
        for surface in 1:6
            kinds,tags,nodes=A.mesh.get_elements(2,surface)
            @test kinds==Int32[3] && length(only(tags))==1
            push!(quads,Tuple(sort(only(nodes))))
        end
        @test all(count->count in (1,2),values(incidence))
        @test Set(face for (face,count) in incidence if count==1)==Set(quads)
        @test Set(Tuple(sort(cell[1:4])) for cell in eachcol(cells[7]))==Set(quads)
        pairs=correspondence(false,4)
        _,master_tags,master_nodes=A.mesh.get_elements(2,1)
        _,slave_tags,slave_nodes=A.mesh.get_elements(2,2)
        @test map(uid->pairs[uid],only(master_nodes))==only(slave_nodes)
        @test isempty(intersect(only(master_tags),only(slave_tags)))
        lower=Dict(uid=>points[uid] for uid in vcat(only(master_nodes),only(slave_nodes)))
        A.mesh.set_order(2)
        @test Set(first(A.mesh.get_elements(3,1)))==Set(Int32[11,14])
        @test all(first(A.mesh.get_node(uid))==point for (uid,point) in lower)
        correspondence(false,4)
        correspondence(true,9)
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end

@testset "explicitly transfinite periodic slave leaves actual mesh unchanged" begin
    # An explicit transfinite slave bypasses the periodic copy path — the
    # hybrid planner cannot certify its boundary against actual copy events,
    # so volume generation rejects atomically without touching the mesh.
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube(;copy_slave=false,affine=TRANSLATION)
        A.mesh.generate(2)
        before=snapshot()
        @test_throws ArgumentError A.mesh.generate(3)
        @test snapshot()==before
    finally
        A.finalize()
    end
end

@testset "retained periodic Surface copy history feeds one hybrid Region" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        cube()
        A.mesh.generate(2)
        original=correspondence(false,4)
        lower=Dict((dim,tag)=>A.mesh.get_elements(dim,tag) for (dim,tag) in A.model.get_entities()
                   if dim<3)
        lower_nodes=A.mesh.get_nodes()
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        @test Set(first(A.mesh.get_elements(3,1)))==Set(Int32[4,7])
        @test correspondence(false,4)==original
        @test all(A.mesh.get_elements(owner...)==cells for (owner,cells) in lower)
        tags,xyz,_=lower_nodes
        @test all(first(A.mesh.get_node(uid))==xyz[3row-2:3row]
                  for (row,uid) in enumerate(tags))
        @test Tessella.validate(A.mesh.get()).ok
    finally
        A.finalize()
    end
end
end
